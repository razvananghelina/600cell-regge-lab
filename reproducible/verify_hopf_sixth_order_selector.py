#!/usr/bin/env python3
"""Exact exhaustive audit of the first canonical higher Hopf-axis moment.

The exhaustive protocol was frozen in commit 66d89d3.  Candidate values had
already been seen exploratorily and are declared in that preregistration; this
script supplies the missing Groebner and quotient-dimension proof.
"""

from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp


OUTPUT = Path(__file__).with_name("hopf_sixth_order_selector.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def exact_zero(expression):
    return sp.simplify(sp.expand(expression)) == 0


def dot(left, right):
    return sum(a*b for a, b in zip(left, right))


def projector(vector):
    norm_square = sp.simplify(dot(vector, vector))
    return tuple(tuple(sp.radsimp(sp.simplify(left*right/norm_square))
                       for right in vector) for left in vector)


def same_matrix(left, right):
    return all(exact_zero(left[row][col]-right[row][col])
               for row in range(len(left))
               for col in range(len(left[0])))


def deduplicate_lines(vectors):
    result = []
    for vector in vectors:
        candidate = projector(vector)
        if not any(same_matrix(candidate, existing[1])
                   for existing in result):
            result.append((vector, candidate))
    return result


def quadratic_form(matrix, vector):
    return sp.expand(sum(vector[row]*matrix[row][col]*vector[col]
                         for row in range(3) for col in range(3)))


def value_on_line(polynomial, variables, vector):
    norm_square = sp.simplify(dot(vector, vector))
    value = polynomial.subs(dict(zip(variables, vector)))
    return sp.simplify(value/norm_square**3)


def line_is_critical(polynomial, variables, vector, lagrange_value):
    norm_square = sp.simplify(dot(vector, vector))
    gradient = [sp.diff(polynomial, variable).subs(dict(zip(variables, vector)))
                for variable in variables]
    # For a homogeneous degree-six polynomial and n=v/sqrt(v.v), the unit
    # Lagrange equation is equivalent to grad S(v)=2 lambda (v.v)^2 v.
    return all(exact_zero(gradient[index]
                          - 2*lagrange_value*norm_square**2*vector[index])
               for index in range(3))


def quotient_dimension(groebner_basis):
    """Count standard monomials of a zero-dimensional lex basis exactly."""
    leading = [polynomial.LM(order=groebner_basis.order).exponents
               for polynomial in groebner_basis.polys]
    dimension = len(groebner_basis.gens)
    bounds = []
    for variable in range(dimension):
        pure_powers = [monomial[variable] for monomial in leading
                       if monomial[variable] > 0
                       and all(monomial[other] == 0
                               for other in range(dimension)
                               if other != variable)]
        if not pure_powers:
            raise RuntimeError("basis is not visibly zero-dimensional")
        bounds.append(min(pure_powers))
    standard = []
    for exponent in product(*(range(bound) for bound in bounds)):
        divisible = any(all(exponent[index] >= monomial[index]
                            for index in range(dimension))
                        for monomial in leading)
        if not divisible:
            standard.append(exponent)
    return len(standard), leading


print("="*78)
print("EXACT SIXTH-ORDER HOPF SELECTOR")
print("="*78)

sqrt5 = sp.sqrt(5)
phi = (1+sqrt5)/2
x, y, z, lagrange = sp.symbols("x y z lambda", real=True)
variables = (x, y, z)
radius_square = x*x+y*y+z*z

vertices = []
for first, second in product((1, -1), repeat=2):
    vertices.extend((
        (sp.Integer(0), sp.Integer(first), sp.Integer(second)*phi),
        (sp.Integer(first), sp.Integer(second)*phi, sp.Integer(0)),
        (sp.Integer(first)*phi, sp.Integer(0), sp.Integer(second)),
    ))
vertices = sorted(set(vertices), key=lambda item: tuple(map(str, item)))
norms = {sp.simplify(dot(vertex, vertex)) for vertex in vertices}
check("the exact coordinate carrier has twelve equal-radius vertices",
      len(vertices) == 12 and len(norms) == 1,
      f"radius squared={next(iter(norms))}")

edges = []
for left, right in combinations(range(len(vertices)), 2):
    if exact_zero(dot(vertices[left], vertices[right])-phi):
        edges.append((left, right))
edge_set = {tuple(sorted(edge)) for edge in edges}
faces = [face for face in combinations(range(len(vertices)), 3)
         if all(tuple(sorted(edge)) in edge_set
                for edge in combinations(face, 2))]
degrees = [sum(index in edge for edge in edges)
           for index in range(len(vertices))]
check("maximal-dot incidence is the icosahedron",
      len(edges) == 30 and len(faces) == 20 and set(degrees) == {5},
      f"f=(12,{len(edges)},{len(faces)}), degree={set(degrees)}")

vertex_lines = deduplicate_lines(vertices)
edge_vectors = [tuple(vertices[left][coordinate]
                      + vertices[right][coordinate] for coordinate in range(3))
                for left, right in edges]
edge_lines = deduplicate_lines(edge_vectors)
face_vectors = [tuple(sum(vertices[index][coordinate] for index in face)
                      for coordinate in range(3)) for face in faces]
face_lines = deduplicate_lines(face_vectors)
check("vertex/edge/face centres give the 6/15/10 symmetry-axis lines",
      (len(vertex_lines), len(edge_lines), len(face_lines)) == (6, 15, 10),
      f"line counts={(len(vertex_lines), len(edge_lines), len(face_lines))}")

axis_quadratics = [quadratic_form(matrix, variables)
                   for _, matrix in vertex_lines]
S2 = sp.expand(sum(axis_quadratics))
S4 = sp.expand(sum(value**2 for value in axis_quadratics))
S6 = sp.expand(sum(value**3 for value in axis_quadratics))
check("the quadratic equal-weight moment is exactly radial",
      exact_zero(S2-2*radius_square))
check("the quartic equal-weight moment is exactly radial",
      exact_zero(S4-sp.Rational(6, 5)*radius_square**2))
sixth_radial_residual = sp.expand(
    S6-sp.Rational(6, 7)*radius_square**3
)
check("the sixth moment has a nonzero anisotropic part",
      not exact_zero(sixth_radial_residual))

lagrange_equations = [sp.diff(S6, variable)-2*lagrange*variable
                      for variable in variables]
lagrange_equations.append(radius_square-1)
critical_basis = sp.groebner(
    lagrange_equations, x, y, z, lagrange,
    extension=sqrt5, order="grevlex"
)
check("the exact Lagrange ideal is zero-dimensional",
      critical_basis.is_zero_dimensional)

eliminants = []
for polynomial in critical_basis.polys:
    expression = polynomial.as_expr()
    if expression.free_symbols <= {lagrange} and expression != 0:
        eliminants.append(sp.Poly(expression, lagrange,
                                  extension=sqrt5).monic().as_expr())
elimination_polynomial = min(eliminants, key=lambda expression:
                             sp.Poly(expression, lagrange).degree())
expected_eliminant = sp.expand(
    (lagrange-sp.Rational(78, 25))
    *(lagrange-sp.Rational(12, 5))
    *(lagrange-sp.Rational(34, 15))
)
check("the Lagrange multiplier has exactly the preregistered cubic",
      exact_zero(elimination_polynomial-expected_eliminant),
      f"monic eliminant={sp.factor(elimination_polynomial)}")

euler_residual = sum(variable*sp.diff(S6, variable)
                     for variable in variables)-6*S6
check("Euler homogeneity gives lambda=3*S6 on the unit sphere",
      exact_zero(euler_residual))

orbit_data = (
    ("C10/Hopf fivefold", vertex_lines,
     sp.Rational(26, 25), sp.Rational(78, 25), 12),
    ("C4 twofold", edge_lines,
     sp.Rational(4, 5), sp.Rational(12, 5), 30),
    ("C6 threefold", face_lines,
     sp.Rational(34, 45), sp.Rational(34, 15), 20),
)
quotient_records = {}
for label, lines, expected_value, expected_lambda, expected_points in orbit_data:
    values = {value_on_line(S6, variables, vector) for vector, _ in lines}
    critical = all(line_is_critical(S6, variables, vector, expected_lambda)
                   for vector, _ in lines)
    check(f"every {label} line is critical with the exact predicted value",
          values == {expected_value} and critical,
          f"unoriented lines={len(lines)}, S6={next(iter(values))}")

    specialized_equations = [
        sp.diff(S6, variable)-2*expected_lambda*variable
        for variable in variables
    ]+[radius_square-1]
    specialized_basis = sp.groebner(
        specialized_equations, x, y, z,
        extension=sqrt5, order="lex"
    )
    quotient_dim, leading = quotient_dimension(specialized_basis)
    quotient_records[label] = {
        "S6": str(expected_value),
        "lambda": str(expected_lambda),
        "unoriented_lines": len(lines),
        "exhibited_real_points": 2*len(lines),
        "quotient_dimension": quotient_dim,
        "leading_monomials": [list(exponents) for exponents in leading],
    }
    check(f"the {label} orbit exhausts its complex critical fibre",
          specialized_basis.is_zero_dimensional
          and quotient_dim == expected_points == 2*len(lines),
          f"real points={2*len(lines)}, quotient dimension={quotient_dim}")

critical_values = [record[2] for record in orbit_data]
check("compactness and the exhaustive critical values fix global extrema",
      min(critical_values) == sp.Rational(34, 45)
      and max(critical_values) == sp.Rational(26, 25),
      "min=34/45 on ten C6 lines; max=26/25 on six C10/Hopf lines")

payload = {
    "protocol_commit": "66d89d3",
    "provenance": "candidate values declared before exhaustive verifier",
    "icosahedron": {"vertices": len(vertices), "edges": len(edges),
                    "faces": len(faces)},
    "axis_lines": {"C10_fivefold": len(vertex_lines),
                   "C4_twofold": len(edge_lines),
                   "C6_threefold": len(face_lines)},
    "moments": {
        "S2": "2*(x^2+y^2+z^2)",
        "S4": "6/5*(x^2+y^2+z^2)^2",
        "S6": str(sp.collect(S6, sqrt5)),
        "S6_radial_residual_nonzero": not exact_zero(sixth_radial_residual),
    },
    "lagrange": {
        "zero_dimensional": critical_basis.is_zero_dimensional,
        "monic_elimination_polynomial": str(sp.factor(elimination_polynomial)),
        "lambda_equals_3_S6": exact_zero(euler_residual),
    },
    "critical_orbits": quotient_records,
    "global": {
        "minimum": "34/45",
        "minimum_orbit": "10 unoriented C6/threefold lines",
        "maximum": "26/25",
        "maximum_orbit": "6 unoriented C10/Hopf/fivefold lines",
    },
    "conditional_selector": (
        "-g*S6 with g>0 has the six unoriented Hopf axes as its exact "
        "degenerate global minima; neither g nor its sign is derived"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("exact sixth-order audit JSON was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: equal quadratic and quartic moments are radial; the sixth is not.")
print("DERIVED: S6 has only three critical orbits and values on S^2:")
print("         C10/Hopf max 26/25, C4 intermediate 4/5, C6 min 34/45.")
print("DERIVED CONDITIONAL: -g*S6 with g>0 selects six degenerate Hopf axes.")
print("OPEN: the theory supplies neither this potential, its sign, nor g.")
raise SystemExit(0 if passed == tests else 1)
