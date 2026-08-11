#!/usr/bin/env python3
"""STEP 2 target comparison for the full Hopf-label Hessian spectrum.

The primitive list was frozen in commit 21a988e before this file compared it
with any distinguished Hopf direction.  Exact target, stationarity and
extremum tests are performed here without fitted combinations.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp

from verify_hopf_fibration_invariants import (
    build_2I,
    build_adjacency,
    build_fiber_adjacency,
    find_all_hopf_fibrations,
)


ENUMERATION = Path(__file__).with_name(
    "hopf_full_hessian_spectral_enumeration.json"
)
OUTPUT = Path(__file__).with_name(
    "hopf_full_hessian_spectral_target.json"
)
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


def span_coordinates(target, basis_polynomials, variables):
    polynomials = [sp.Poly(sp.expand(value), *variables)
                   for value in basis_polynomials]
    target_poly = sp.Poly(sp.expand(target), *variables)
    monomials = sorted(set(target_poly.monoms()).union(*(
        set(polynomial.monoms()) for polynomial in polynomials
    )))
    matrix = sp.Matrix([
        [polynomial.coeff_monomial(monomial)
         for polynomial in polynomials]
        for monomial in monomials
    ])
    vector = sp.Matrix([
        target_poly.coeff_monomial(monomial) for monomial in monomials
    ])
    if matrix.rank() != matrix.row_join(vector).rank():
        return None
    solution = sp.linsolve((matrix, vector))
    values = next(iter(solution))
    if any(value.free_symbols for value in values):
        return None
    return tuple(sp.simplify(value) for value in values)


def exact_sign(value):
    value = sp.simplify(value)
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    if value.is_zero:
        return 0
    raise RuntimeError(f"undecided exact sign: {value}")


def constrained_data(polynomial, point, q, variables):
    gradient = sp.Matrix([sp.diff(polynomial, variable)
                          for variable in variables]).subs(dict(zip(variables,
                                                                   point)))
    q_gradient = sp.Matrix([sp.diff(q, variable)
                            for variable in variables]).subs(dict(zip(variables,
                                                                     point)))
    ratios = [sp.simplify(gradient[index]/q_gradient[index])
              for index in range(len(variables)) if q_gradient[index] != 0]
    stationary = bool(ratios and len(set(ratios)) == 1
                      and gradient == ratios[0]*q_gradient)
    lagrange_multiplier = ratios[0] if stationary else None
    tangent_columns = sp.Matrix([list(q_gradient)]).nullspace()
    tangent_basis = sp.Matrix.hstack(*tangent_columns)
    if not stationary:
        return {
            "stationary": False,
            "lambda": None,
            "signature": None,
            "restricted_hessian": None,
        }
    lagrangian_hessian = (
        sp.hessian(polynomial, variables)
        - lagrange_multiplier*sp.hessian(q, variables)
    ).subs(dict(zip(variables, point)))
    restricted = sp.simplify(tangent_basis.T*lagrangian_hessian*tangent_basis)
    eigenvalues = restricted.eigenvals()
    signature = [0, 0, 0]
    for eigenvalue, multiplicity in eigenvalues.items():
        sign = exact_sign(eigenvalue)
        signature[{1: 0, -1: 1, 0: 2}[sign]] += multiplicity
    return {
        "stationary": True,
        "lambda": lagrange_multiplier,
        "signature": tuple(signature),
        "restricted_hessian": restricted,
    }


print("="*78)
print("FULL-HESSIAN SPECTRAL TARGET COMPARISON -- STEP 2")
print("="*78)

enumeration = json.loads(ENUMERATION.read_text())
check("the committed blind enumeration precedes target comparison",
      enumeration["protocol_commit"] == "3767638"
      and enumeration["target_comparison_performed"] is False
      and enumeration["primitive_spectral_count_N"] == 5)

u = sp.symbols("u0:5")
elementary = {
    degree: sp.sympify(value, locals={str(variable): variable for variable in u})
    for degree, value in enumeration["elementary_coefficients"].items()
}
e3, e4, e5 = elementary["3"], elementary["4"], elementary["5"]

vertices = build_2I()
adjacency = np.rint(build_adjacency(vertices)).astype(np.int64)
fibrations = find_all_hopf_fibrations(vertices)
fiber_adjacencies = [
    np.rint(build_fiber_adjacency(adjacency, fibration)).astype(np.int64)
    for fibration in fibrations
]
boxes = [6*fiber-adjacency for fiber in fiber_adjacencies]
field_basis = [boxes[index]-boxes[5] for index in range(5)]
field_gram = sp.Matrix([
    [trace_product(field_basis[row], field_basis[col])
     for col in range(5)]
    for row in range(5)
])
q = sp.expand((sp.Matrix(u).T*field_gram*sp.Matrix(u))[0])
check("the target audit independently reconstructs the committed norm",
      str(q) == enumeration["arena"]["field_norm_polynomial"])

old_cubic = sp.expand(sum(
    u[a]*u[b]*u[c]*trace_product(field_basis[a], field_basis[b],
                                field_basis[c])
    for a in range(5) for b in range(5) for c in range(5)
))
overlaps = [sp.expand(sum(
    u[a]*trace_product(field_basis[a], box) for a in range(5)
)) for box in boxes]
projector_cubic = sp.expand(sum(overlap**3 for overlap in overlaps))
projector_quartic = sp.expand(sum(overlap**4 for overlap in overlaps))
projector_quintic = sp.expand(sum(overlap**5 for overlap in overlaps))
check("the old and equal-weight cubics are independent",
      span_coordinates(old_cubic, [projector_cubic], u) is None
      and old_cubic != 0 and projector_cubic != 0)

e3_coordinates = span_coordinates(e3, [old_cubic, projector_cubic], u)
e4_coordinates = span_coordinates(e4, [q**2, projector_quartic], u)
e5_coordinates = span_coordinates(e5, [q*old_cubic,
                                        q*projector_cubic,
                                        projector_quintic], u)
print(f"  e3 in (old,C_box) = {e3_coordinates}")
print(f"  e4 in (q^2,sum overlap^4) = {e4_coordinates}")
print(f"  e5 in (q*old,q*C_box,sum overlap^5) = {e5_coordinates}")
check("the cubic primitive lies in the complete two-cubic invariant span",
      e3_coordinates == (0, sp.Integer(-7776)),
      "e3=-7776*C_box, with zero old-cubic component")
check("the quartic primitive lies in the symmetric norm/quartic span",
      e4_coordinates
      == (sp.Integer(15672832819200), sp.Integer(139968)))
check("the quintic is independent of the declared inherited degree-five span",
      e5_coordinates is None,
      "e5 is a new A5-invariant quintic, not q times either known cubic")

# The overlap map is an exact isomorphism from W to the zero-sum hyperplane.
# Its two constraints turn the global e3/e4 extremum questions into finite
# Lagrange-multiplier enumerations with no numerical search.
overlap_matrix = sp.Matrix(6, 5, lambda row, col:
                           sp.diff(overlaps[row], u[col]))
overlap_norm = sp.expand(sum(overlap**2 for overlap in overlaps))
check("overlaps identify W with the full zero-sum label hyperplane",
      overlap_matrix.rank() == 5
      and sp.expand(sum(overlaps)) == 0
      and sp.expand(overlap_norm-8640*q) == 0,
      "sum s_i=0 and sum s_i^2=8640*q")

coordinate_vertices = []
for index in range(5):
    coordinate = sp.Matrix([sp.Rational(-1, 6)]*5)
    coordinate[index] = sp.Rational(5, 6)
    coordinate_vertices.append(coordinate)
coordinate_vertices.append(sp.Matrix([sp.Rational(-1, 6)]*5))
signed_vertices = coordinate_vertices+[-point for point in coordinate_vertices]

primitive_polynomials = {3: e3, 4: e4, 5: e5}
vertex_data = {}
for degree, polynomial in primitive_polynomials.items():
    records = []
    for point in signed_vertices:
        data = constrained_data(polynomial, point, q, u)
        data["value"] = sp.simplify(polynomial.subs(dict(zip(u, point))))
        records.append(data)
    vertex_data[degree] = records
    print(f"  e{degree} +Box value={records[0]['value']}")
    print(f"  e{degree} +Box tangent signature={records[0]['signature']}")
    print(f"  e{degree} -Box tangent signature={records[6]['signature']}")

check("all three nonconstant primitives are stationary at all signed vertices",
      all(record["stationary"]
          for records in vertex_data.values() for record in records))
check("the signed-orbit values obey the exact parity of each degree",
      all(vertex_data[degree][index+6]["value"]
          == (-1)**degree*vertex_data[degree][index]["value"]
          for degree in primitive_polynomials for index in range(6)))

# Global cubic certificate.  At a constrained stationary point of sum s_i^3,
# every coordinate is a root of the same quadratic, so there are at most two
# values.  Enumerating their multiplicity k=1,...,5 is exhaustive.
R2 = sp.Integer(8640*7200)
cubic_stationary_values = {}
quartic_two_value_ratios = {}
for multiplicity in range(1, 6):
    positive_value = sp.sqrt(
        R2*sp.Rational(6-multiplicity, 6*multiplicity)
    )
    negative_value = -sp.Rational(multiplicity, 6-multiplicity)*positive_value
    cubic_stationary_values[multiplicity] = sp.simplify(
        multiplicity*positive_value**3
        +(6-multiplicity)*negative_value**3
    )
    # Scale-free fourth-moment ratio for a two-value stationary point.
    a = sp.Integer(6-multiplicity)
    b = sp.Integer(-multiplicity)
    quartic_two_value_ratios[multiplicity] = sp.factor(
        (multiplicity*a**4+(6-multiplicity)*b**4)
        /(multiplicity*a**2+(6-multiplicity)*b**2)**2
    )

target_cubic = sp.simplify(projector_cubic.subs(
    dict(zip(u, coordinate_vertices[0]))
))
check("the exhaustive cubic stationary census has unique k=1 maximum",
      cubic_stationary_values[1] == target_cubic
      and cubic_stationary_values[5] == -target_cubic
      and all(exact_sign(target_cubic-cubic_stationary_values[index]) == 1
              for index in range(2, 6)),
      f"five multiplicity values={cubic_stationary_values}")

target_overlap_vectors = {
    tuple(overlap.subs(dict(zip(u, point))) for overlap in overlaps)
    for point in signed_vertices
}
expected_target_overlap_vectors = set()
for sign in (1, -1):
    for exceptional in range(6):
        expected_target_overlap_vectors.add(tuple(
            sp.Integer(sign*(7200 if index == exceptional else -1440))
            for index in range(6)
        ))
check("the cubic equality cases pull back exactly to the signed Box simplex",
      target_overlap_vectors == expected_target_overlap_vectors
      and overlap_matrix.rank() == 5)
check("e3 and Tr(Hhat^3) have exactly the six +Box global minima",
      e3_coordinates == (0, -7776)
      and sp.expand(sp.sympify(
          enumeration["power_sums"]["3"],
          locals={str(variable): variable for variable in u},
      )-3*e3) == 0,
      "Tr(Hhat_X^3)=3e3=-23328*C_box")

# Global quartic certificate.  Three used roots of the stationary cubic must
# sum to zero.  Their multiplicities exhaust the positive ordered triples
# p+q+r=6.  The exceptional p=q=r=2 family has fixed ratio 1/4.
quartic_three_value_ratios = []
for p in range(1, 5):
    for r in range(1, 6-p):
        s = 6-p-r
        if s < 1:
            continue
        constraint = sp.Matrix([[1, 1, 1], [p, r, s]])
        nullspace = constraint.nullspace()
        if len(nullspace) == 2:
            # This occurs only for (2,2,2); a+b+c=0 gives
            # a^4+b^4+c^4=(a^2+b^2+c^2)^2/2.
            quartic_three_value_ratios.append(sp.Rational(1, 4))
            continue
        root_vector = nullspace[0]
        if len(set(root_vector)) < 3:
            # A degenerate root assignment is already in the two-value list.
            continue
        numerator = (p*root_vector[0]**4+r*root_vector[1]**4
                     +s*root_vector[2]**4)
        denominator = (p*root_vector[0]**2+r*root_vector[1]**2
                       +s*root_vector[2]**2)**2
        quartic_three_value_ratios.append(sp.factor(numerator/denominator))

quartic_target_ratio = sp.Rational(7, 10)
check("the exhaustive quartic stationary census has sharp ratio 7/10",
      quartic_two_value_ratios[1] == quartic_target_ratio
      and quartic_two_value_ratios[5] == quartic_target_ratio
      and all(value < quartic_target_ratio
              for index, value in quartic_two_value_ratios.items()
              if index not in (1, 5))
      and quartic_three_value_ratios
      and all(value < quartic_target_ratio
              for value in quartic_three_value_ratios),
      "equality only for one-versus-five overlap values")
check("e4 has exactly the twelve signed Box points as global maxima",
      e4_coordinates[1] > 0
      and target_overlap_vectors == expected_target_overlap_vectors,
      "e4=15672832819200*q^2+139968*sum_i s_i^4")

# e5 is locally sharp at the target, but the exact global problem is not
# inferred from reconnaissance.  Its independence above prevents reducing it
# to either already certified cubic inequality.
check("the exact local signatures classify all three primitive target orbits",
      all(record["signature"] == (4, 0, 0)
          for record in vertex_data[3][:6])
      and all(record["signature"] == (0, 4, 0)
              for record in vertex_data[3][6:])
      and all(record["signature"] == (0, 4, 0)
              for record in vertex_data[4])
      and all(record["signature"] == (0, 4, 0)
              for record in vertex_data[5][:6])
      and all(record["signature"] == (4, 0, 0)
              for record in vertex_data[5][6:]))

# Gaussian framing audit: reconstruct the actual physical Hessian at a target
# point.  Its (3+,2-) signature makes a real bosonic Euclidean Gaussian
# divergent; squaring it or declaring fermions would be a new input.
box_products = [[boxes[row]@boxes[col] for col in range(6)]
                for row in range(6)]


def label_hessian(X):
    return sp.Matrix(6, 6, lambda row, col:
                     3*(int(np.sum(X*box_products[row][col].T))
                        + int(np.sum(X*box_products[col][row].T))))


label_basis = sp.Matrix(6, 5, lambda row, col:
                        (1 if row == col else (-1 if row == 5 else 0)))
target_hessian = label_hessian(boxes[0])
target_restriction = (target_hessian*label_basis)[:5, :]
target_eigenvalues = target_restriction.eigenvals()
target_spectral_signature = [0, 0, 0]
for eigenvalue, multiplicity in target_eigenvalues.items():
    target_spectral_signature[{1: 0, -1: 1, 0: 2}[
        exact_sign(eigenvalue)
    ]] += multiplicity
check("the physical Hessian at every target is indefinite",
      target_spectral_signature == [3, 2, 0]
      and target_hessian*sp.ones(6, 1) == sp.zeros(6, 1),
      f"Hhat_Box signature={tuple(target_spectral_signature)}")
check("a convergent bosonic Gaussian cannot generate the target determinant",
      target_spectral_signature[1] > 0,
      "determinant/fermionic or H^2 interpretations require new input")

certified_geometric_hits = [3, 4]
positive_minimum_hits = [3]
check("the preregistered look-elsewhere counts are explicit",
      len(certified_geometric_hits) == 2
      and len(positive_minimum_hits) == 1,
      "certified extrema: 2/5 primitives (2/3 nonconstant); +Box minima: 1/5")

payload = {
    "protocol_commits": {
        "preregistration": "3767638",
        "blind_enumeration": "21a988e",
    },
    "primitive_count_N": 5,
    "nonconstant_primitive_count": 3,
    "exact_relations": {
        "e3": "-7776*C_box",
        "Tr_Hhat_cubed": "-23328*C_box",
        "e4": "15672832819200*q^2+139968*sum_i overlap_i^4",
        "e5_in_span_q_old_q_Cbox_overlap5": False,
    },
    "global_target_results": {
        "e3": "global minima exactly +Box_i; maxima exactly -Box_i",
        "e4": "global maxima exactly +/-Box_i",
        "e5": "OPEN globally; exact strict local maxima at +Box_i",
    },
    "hit_fractions": {
        "certified_geometric_extremal_loci_among_N": "2/5",
        "certified_geometric_extremal_loci_among_nonconstant": "2/3",
        "positive_action_global_minimum_plus_Box_among_N": "1/5",
    },
    "Gaussian_audit": {
        "Hhat_Box_signature": target_spectral_signature,
        "convergent_real_bosonic_Gaussian": False,
    },
    "verdict": (
        "STRUCTURAL ADVANCE: the full off-diagonal Hessian has "
        "Tr(Hhat_X^3)=-23328*C_box, whose global minima on q=7200 are "
        "exactly the six +Box_i.  No label superselection or fitted "
        "coefficient is used.  Physical action status remains OPEN because "
        "the extended Hessian field and its positive cubic spectral moment "
        "are not yet derived dynamics; the Hessian is indefinite, so this is "
        "not a convergent bosonic Gaussian one-loop action."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the exact STEP 2 target-comparison JSON was written", OUTPUT.exists())

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("STRUCTURAL ADVANCE: Tr(Hhat_X^3)=-23328*C_box.")
print("DERIVED: its global minima on q=7200 are exactly the six +Box_i.")
print("DERIVED: e4 has the twelve signed Box_i as exact global maxima.")
print("OPEN: e5 is a strict local target extremum; global uniqueness unproved.")
print("OPEN PHYSICS: Hhat is indefinite and no bosonic Gaussian is convergent.")
raise SystemExit(0 if passed == tests else 1)
