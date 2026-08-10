#!/usr/bin/env python3
"""Exact higher single-trace audit on the Hopf--Box module.

The finite word space and the privileged sixth-moment test were frozen in
commit aea13cb.  Commit 35566e8 corrected the decision boundary, before the
full sixth moment was evaluated, so that complete dynamic selection rather
than cubic purity is decisive.

All word values are reconstructed exactly from six modular computations.  A
proved operator-norm bound makes the CRT modulus more than twice the largest
possible absolute value, so signed reconstruction is unique rather than
probabilistic.
"""

from itertools import combinations_with_replacement, product
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.ntheory.modular import crt

from verify_hopf_fibration_invariants import (
    build_2I,
    build_adjacency,
    build_fiber_adjacency,
    find_all_hopf_fibrations,
)


OUTPUT = Path(__file__).with_name("hopf_higher_single_trace.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def trace_product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = result @ matrix
    return int(np.trace(result))


def projective_key(left, right):
    divisor = math.gcd(abs(left), abs(right))
    left //= divisor
    right //= divisor
    if left < 0 or (left == 0 and right < 0):
        left = -left
        right = -right
    return left, right


print("="*78)
print("EXACT HIGHER SINGLE-TRACE HOPF AUDIT")
print("="*78)

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
basis = [boxes[index]-boxes[5] for index in range(5)]
gram = sp.Matrix([
    [trace_product(basis[row], basis[col]) for col in range(5)]
    for row in range(5)
])
check("the fixed geometry reconstructs the five-dimensional Hopf--Box space",
      len(boxes) == 6 and np.array_equal(sum(boxes), np.zeros_like(adjacency))
      and gram.rank() == 5)

# Exact minimal polynomial of the adjacency.  Horner evaluation is entirely
# integral and remains safely inside int64 for this degree.
z = sp.symbols("z")
minimal_polynomial = sp.Poly(
    sp.expand(z*(z-12)*(z-3)*(z+2)*(z+3)
              *(z**2-6*z-36)*(z**2-4*z-16)),
    z,
    domain=sp.ZZ,
)
expected_coefficients = [
    1, -20, 39, 940, -1584, -18360, -3456, 103680, 124416, 0
]
identity = np.eye(120, dtype=np.int64)
residual = minimal_polynomial.all_coeffs()[0]*identity
for coefficient in minimal_polynomial.all_coeffs()[1:]:
    residual = residual@adjacency+int(coefficient)*identity

# Prove I,A,...,A^8 independent over Q.  Modular row reduction identifies a
# 9-by-9 integer minor; its nonzero exact determinant is the certificate.
powers = [identity]
for _ in range(8):
    powers.append(powers[-1]@adjacency)
flat_power_matrix = np.column_stack([matrix.reshape(-1) for matrix in powers])
rank_prime = 1000003
selected_rows = []
modular_basis = []
for row_index, row in enumerate(flat_power_matrix):
    candidate = [list(vector) for vector in modular_basis]
    vector = [int(value % rank_prime) for value in row]
    candidate.append(vector)
    if int(sp.polys.matrices.DomainMatrix.from_list(
        candidate, sp.GF(rank_prime)
    ).rank()) > len(modular_basis):
        modular_basis.append(vector)
        selected_rows.append(row_index)
    if len(selected_rows) == 9:
        break
independence_minor = sp.Matrix(flat_power_matrix[selected_rows, :].tolist())
check("m_A has degree nine and annihilates A coefficientwise",
      minimal_polynomial.all_coeffs() == expected_coefficients
      and np.count_nonzero(residual) == 0)
check("no polynomial of degree below nine annihilates A",
      len(selected_rows) == 9 and independence_minor.det() != 0,
      f"exact independence minor det={independence_minor.det()}")

# The exact A5 character calculation gives two invariant cubics.  Values by
# class are chi_W=(5,1,-1,0); the symmetric-cube character values are
# (35,3,2,0), with class sizes (1,15,20,24).
sym3_character_sum = 1*35+15*3+20*2+24*0
invariant_dimension = sym3_character_sum//60
check("the independently evaluated invariant-cubic dimension is two",
      sym3_character_sum == 120 and invariant_dimension == 2)

# Two exact evaluation points.  The first is Box_0.  The second is the integer
# direction from the fourth-moment counterexample, before normalization.
X_first = boxes[0]
X_second = sum((coefficient*direction for coefficient, direction in
                zip((-1, 1, 1, -1, 1), basis)), np.zeros_like(adjacency))

def G0_value(X):
    return trace_product(X, X, X)


def C_box_value(X):
    return sum(trace_product(X, box)**3 for box in boxes)


evaluation_basis = sp.Matrix([
    [G0_value(X_first), C_box_value(X_first)],
    [G0_value(X_second), C_box_value(X_second)],
])
evaluation_determinant = evaluation_basis.det()
check("the two-point evaluation map is exact and nondegenerate",
      evaluation_basis == sp.Matrix([
          [14400, 358318080000],
          [622080, 0],
      ])
      and evaluation_determinant == -222902511206400000,
      f"det={evaluation_determinant}")

triples = list(combinations_with_replacement(range(9), 3))
check("the complete reduced individual-word list has N=165",
      len(triples) == math.comb(11, 3) == 165)

# Exact word values by CRT.  For symmetric A and X,
# |Tr(A^a X A^b X A^c X)| <= 120*12^(a+b+c)*||X||^3.
# The maximum absolute row sum bounds ||X||, and a+b+c<=24.
primes = [1000003, 1000033, 1000037, 1000039, 1000081, 1000099]
modulus = math.prod(primes)
row_bounds = [
    int(np.max(np.sum(np.abs(X_first), axis=1))),
    int(np.max(np.sum(np.abs(X_second), axis=1))),
]
absolute_bounds = [120*12**24*row_bound**3 for row_bound in row_bounds]
check("the CRT modulus proves unique signed reconstruction for both points",
      row_bounds == [20, 72]
      and all(modulus > 2*bound for bound in absolute_bounds),
      f"modulus={modulus}; max bound={max(absolute_bounds)}")


def modular_word_values(X, prime):
    adjacency_mod = adjacency % prime
    powers_mod = [np.eye(120, dtype=np.int64)]
    for _ in range(8):
        powers_mod.append(powers_mod[-1]@adjacency_mod % prime)
    left = [matrix@(X % prime) % prime for matrix in powers_mod]
    pairs = [[left[first]@left[second] % prime for second in range(9)]
             for first in range(9)]
    return {
        triple: int(np.sum((pairs[triple[0]][triple[1]]
                            *left[triple[2]].T) % prime) % prime)
        for triple in triples
    }


residues_first = {triple: [] for triple in triples}
residues_second = {triple: [] for triple in triples}
for prime in primes:
    first_mod = modular_word_values(X_first, prime)
    second_mod = modular_word_values(X_second, prime)
    for triple in triples:
        residues_first[triple].append(first_mod[triple])
        residues_second[triple].append(second_mod[triple])


def signed_crt(residues):
    value = int(crt(primes, residues, check=True)[0])
    return value if value <= modulus//2 else value-modulus


word_values = {
    triple: (signed_crt(residues_first[triple]),
             signed_crt(residues_second[triple]))
    for triple in triples
}

# One unused prime provides an independent residue cross-check after signed
# reconstruction.
extra_prime = 1000117
extra_first = modular_word_values(X_first, extra_prime)
extra_second = modular_word_values(X_second, extra_prime)
check("all 330 reconstructed integers pass an unused-prime cross-check",
      all(word_values[triple][0] % extra_prime == extra_first[triple]
              and word_values[triple][1] % extra_prime == extra_second[triple]
              for triple in triples))

projective_lines = {}
selector_hits = []
old_line_hits = []
g_first, c_first = map(int, evaluation_basis.row(0))
g_second, c_second = map(int, evaluation_basis.row(1))
for triple, (first, second) in word_values.items():
    projective_lines.setdefault(projective_key(first, second), []).append(triple)
    if first*c_second-second*c_first == 0:
        selector_hits.append(triple)
    if first*g_second-second*g_first == 0:
        old_line_hits.append(triple)
check("the 165 words occupy exactly 24 distinct invariant cubic lines",
      len(projective_lines) == 24
      and sorted(len(items) for items in projective_lines.values())
      == [1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 7, 7, 8, 8,
          10, 10, 11, 11, 12, 12, 12, 12, 13])
check("zero of 165 individual canonical words equals the Hopf selector",
      selector_hits == [],
      "hit fraction=0/165")
check("only two words remain on the old cubic line",
      old_line_hits == [(0, 0, 0), (0, 0, 1)],
      f"old-line triples={old_line_hits}")
check("the complete adjacency-polynomial word span reaches both invariant lines",
      len(projective_lines) > 1 and invariant_dimension == 2,
      "span rank=2; arbitrary linear combinations would introduce fitting")

# For n=p-3, cyclically mark one of the three X insertions.  Every length-p
# word is thereby counted three times, giving the exact noncommutative formula
# K_p=(p/3) sum_{a+b+c=n} T_abc.
def K_evaluations(moment):
    total_first = 0
    total_second = 0
    degree = moment-3
    for first in range(degree+1):
        for second in range(degree-first+1):
            third = degree-first-second
            values = word_values[tuple(sorted((first, second, third)))]
            total_first += values[0]
            total_second += values[1]
    factor = sp.Rational(moment, 3)
    return sp.Matrix([factor*total_first, factor*total_second])


moment_coordinates = {}
for moment in range(3, 7):
    evaluations = K_evaluations(moment)
    coordinates = evaluation_basis.inv()*evaluations
    moment_coordinates[moment] = tuple(map(sp.factor, coordinates))
expected_coordinates = {
    3: (sp.Integer(1), sp.Integer(0)),
    4: (sp.Integer(-8), sp.Integer(0)),
    5: (sp.Integer(40), sp.Rational(1, 1244160)),
    6: (sp.Integer(-280), sp.Rational(1, 124416)),
}
check("K_3,...,K_6 have the exact preregistered basis coordinates",
      moment_coordinates == expected_coordinates,
      str(moment_coordinates))
check("the sixth-moment cubic is mixed rather than the pure selector",
      moment_coordinates[6] == (-280, sp.Rational(1, 124416)),
      "K6=-280 Tr(X^3)+C_box/124416")
posthoc_cancellation = tuple(
    sp.simplify(moment_coordinates[6][index]
                - 35*moment_coordinates[4][index])
    for index in range(2)
)
check("the post-comparison coefficient 35 cancels the unwanted cubic exactly",
      posthoc_cancellation == (0, sp.Rational(1, 124416)),
      "K6-35*K4=C_box/124416; recorded as PATTERN, not a derivation")

# Full sixth moment on the fixed q-sphere.  First test all desired signed
# vertices directly in integer arithmetic.
vertex_coordinates = []
for index in range(5):
    coordinate = [sp.Rational(-1, 6)]*5
    coordinate[index] = sp.Rational(5, 6)
    vertex_coordinates.append(sp.Matrix(coordinate))
vertex_coordinates.append(sp.Matrix([sp.Rational(-1, 6)]*5))

signed_vertex_data = {}
for sign, expected_multiplier in ((1, 38880), (-1, 69360)):
    values = []
    stationary = []
    for index, box in enumerate(boxes):
        D = adjacency+sign*box
        D_fifth = np.linalg.matrix_power(D, 5)
        values.append(trace_product(D_fifth, D))
        gradient_six = sp.Matrix([
            6*int(np.sum(D_fifth*direction.T)) for direction in basis
        ])
        gradient_q = 2*gram*(sign*vertex_coordinates[index])
        stationary.append(gradient_six == expected_multiplier*gradient_q)
    signed_vertex_data[sign] = (values, stationary, expected_multiplier)
check("all +/-Box_i are exact stationary points of the full sixth moment",
      all(signed_vertex_data[1][1]) and all(signed_vertex_data[-1][1]),
      "multipliers: +38880 and -69360")
check("the full sixth-moment values at signed vertices are exact",
      set(signed_vertex_data[1][0]) == {111974400}
      and set(signed_vertex_data[-1][0]) == {200678400})

# The cubic cancellation corresponds to the order-six truncated heat ratio
# t=3/35: up to an overall positive factor its angular functional is
# S4-S6/35.  Although its cubic part has the desired -C_box sign, the complete
# polynomial orders the signed vertices in the wrong direction.
truncated_heat_positive = sp.Integer(933120)-sp.Rational(111974400, 35)
truncated_heat_negative = sp.Integer(1163520)-sp.Rational(200678400, 35)
check("the seductive t=3/35 truncated-heat cancellation fails in full",
      truncated_heat_positive == -sp.Rational(15863040, 7)
      and truncated_heat_negative == -sp.Rational(31991040, 7)
      and truncated_heat_negative < truncated_heat_positive,
      "despite pure cubic sign, -Box_i is lower than +Box_i by 2304000")

# The same integer direction used as the second evaluation point has
# q=51840, so s=sqrt(5)/6 reaches q=7200.  Expand the noncommutative sixth
# power by all 64 words, then check both signs and their gradients exactly.
scale = sp.sqrt(5)/6


def trace_expansion(direction, degree):
    coefficients = [0]*(degree+1)
    for choices in product((0, 1), repeat=degree):
        matrices = [direction if choice else adjacency for choice in choices]
        coefficients[sum(choices)] += trace_product(*matrices)
    return coefficients


sixth_coefficients = trace_expansion(X_second, 6)


def sixth_value(signed_scale):
    return sp.expand(sum(coefficient*signed_scale**power
                         for power, coefficient
                         in enumerate(sixth_coefficients)))


lower_value = sixth_value(scale)
higher_value = sixth_value(-scale)
check("the exact competitor pair is normalized to q=7200",
      trace_product(X_second, X_second) == 51840
      and sp.simplify(scale**2*51840) == 7200)
check("a stationary competitor lies below the desired vertices for +S6",
      lower_value == 165542400-38592000*sp.sqrt(5)
      and sp.simplify(lower_value-111974400) < 0,
      f"S6={lower_value}")
check("the opposite competitor beats them for the fixed heat sign -S6",
      higher_value == 165542400+38592000*sp.sqrt(5)
      and sp.simplify(higher_value-111974400) > 0,
      f"S6={higher_value}; hence -S6 is lower there")


def sixth_gradient(signed_scale):
    result = []
    for test_direction in basis:
        coefficients = [0]*6
        for choices in product((0, 1), repeat=5):
            matrices = [X_second if choice else adjacency
                        for choice in choices]
            coefficients[sum(choices)] += 6*trace_product(
                *matrices, test_direction
            )
        result.append(sp.expand(sum(
            coefficient*signed_scale**power
            for power, coefficient in enumerate(coefficients)
        )))
    return sp.Matrix(result)


direction_coordinates = sp.Matrix([-1, 1, 1, -1, 1])
stationary_competitors = []
competitor_multipliers = []
for signed_scale, multiplier in (
    (scale, 57240-12840*sp.sqrt(5)),
    (-scale, 57240+12840*sp.sqrt(5)),
):
    gradient_six = sixth_gradient(signed_scale)
    gradient_q = 2*gram*(signed_scale*direction_coordinates)
    residuals = [sp.simplify(left-multiplier*right)
                 for left, right in zip(gradient_six, gradient_q)]
    stationary_competitors.append(all(residual == 0 for residual in residuals))
    competitor_multipliers.append(multiplier)
check("both exact sixth-moment competitors are constrained stationary points",
      all(stationary_competitors),
      f"multipliers={competitor_multipliers}")

payload = {
    "protocol_commits": ["aea13cb", "35566e8"],
    "minimal_polynomial": str(minimal_polynomial.as_expr()),
    "complete_word_space": {
        "exponent_range": [0, 8],
        "N_words": len(triples),
        "distinct_projective_cubic_lines": len(projective_lines),
        "selector_hits": len(selector_hits),
        "selector_hit_fraction": "0/165",
        "old_line_words": [list(triple) for triple in old_line_hits],
        "span_rank": 2,
        "CRT_primes": primes,
        "CRT_modulus": modulus,
        "unused_check_prime": extra_prime,
    },
    "moment_cubic_coordinates_in_G0_Cbox": {
        str(moment): [str(value) for value in coordinates]
        for moment, coordinates in moment_coordinates.items()
    },
    "sixth_moment": {
        "K6": "-280*Tr(X^3)+C_box/124416",
        "posthoc_cubic_cancellation": "K6-35*K4=C_box/124416",
        "posthoc_status": "PATTERN; coefficient found after target comparison",
        "truncated_heat_t": "3/35",
        "truncated_heat_failure": "S4-S6/35 is lower on -Box_i by 2304000",
        "S6_positive_Box": 111974400,
        "S6_negative_Box": 200678400,
        "desired_vertices_stationary": True,
        "competitor_lower_for_positive_S6": str(lower_value),
        "competitor_higher_for_heat_sign_minus_S6": str(higher_value),
        "competitors_stationary": True,
    },
    "verdict": (
        "DERIVED NEGATIVE: no individual reduced single-trace word equals "
        "the Hopf selector; K6 is a fixed mixture, and the complete sixth "
        "moment has exact stationary competitors defeating the desired six "
        "for both +S6 and the formal heat-trace sign -S6"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact higher-single-trace audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: 165 reduced words occupy 24 cubic lines and span both")
print("         A5-invariant directions.")
print("DERIVED NEGATIVE: selector hits = 0/165; K6 is a mixed cubic.")
print("DERIVED NEGATIVE: exact stationary competitors defeat both +S6 and -S6.")
print("KILL: the canonical sixth single-trace moment does not select six Hopf")
print("      fibrations; arbitrary word combinations would be fitting.")
raise SystemExit(0 if passed == tests else 1)
