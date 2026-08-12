#!/usr/bin/env python3
"""Certified asymmetric Lorentzian Regge tent witness.

Protocol commit: 4a63f66.  The numerical witness and every decision boundary
were disclosed before this verifier was written.  A passing run proves an
admissible stationary pole for the frozen asymmetric boundary data; it does
not turn the target-found boundary data into a selected physical clock.
"""

import itertools
import json
import math
from pathlib import Path

from flint import arb, ctx
import mpmath as mp
import networkx as nx
import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_lorentzian_asymmetric_tent.json"
PROTOCOL_COMMIT = "4a63f66"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


class Dual:
    """First-order dual number whose coefficients may be Arb balls."""

    def __init__(self, value, derivative=0):
        self.v = value if isinstance(value, arb) else arb(value)
        self.d = derivative if isinstance(derivative, arb) else arb(derivative)

    @staticmethod
    def lift(other):
        return other if isinstance(other, Dual) else Dual(other)

    def __add__(self, other):
        other = self.lift(other)
        return Dual(self.v + other.v, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        return Dual(self.v * other.v, self.d * other.v + self.v * other.d)

    __rmul__ = __mul__

    def reciprocal(self):
        return Dual(1 / self.v, -self.d / (self.v * self.v))

    def __truediv__(self, other):
        return self * self.lift(other).reciprocal()

    def __rtruediv__(self, other):
        return self.lift(other) / self

    def sqrt(self):
        root = self.v.sqrt()
        return Dual(root, self.d / (2 * root))

    def acos(self):
        return Dual(self.v.acos(), -self.d / (1 - self.v * self.v).sqrt())


def projected_components(rho, q1, q2, q3):
    """Projected normal-plane Gram entries for one ordered hinge."""
    c1 = (1 - rho - q1) / 2
    c2 = (1 - rho - q2) / 2
    c3 = (1 - rho - q3) / 2
    denominator = 4 * (c1 * c1 + rho)
    p11 = (4*c1*c1 - 4*c1*c2 + 4*c2*c2 + 3*rho) / denominator
    p22 = (4*c1*c1 - 4*c1*c3 + 4*c3*c3 + 3*rho) / denominator
    p12 = (
        2*c1*c1 - 2*c1*c2 - 2*c1*c3 + 4*c2*c3 + rho
    ) / denominator
    return p11, p22, p12


def theta_dual(rho, q1, q2, q3):
    p11, p22, p12 = projected_components(rho, q1, q2, q3)
    return (p12 / (p11*p22).sqrt()).acos()


def weight_dual(rho, q):
    radicand = 4*rho + (1-q-rho)*(1-q-rho)
    return (1+q+rho) / (4*radicand.sqrt())


def determinant_3(matrix):
    a, b, c = matrix[0]
    _, d, e = matrix[1]
    _, _, f = matrix[2]
    return a*d*f + 2*b*c*e - a*e*e - d*c*c - f*b*b


def final_tetra_gram(values):
    return [
        [
            values[i] if i == j else (values[i]+values[j]-1)/2
            for j in range(3)
        ]
        for i in range(3)
    ]


def principal_minors(matrix):
    result = []
    for size in (1, 2, 3):
        for subset in itertools.combinations(range(3), size):
            if size == 1:
                value = matrix[subset[0]][subset[0]]
            elif size == 2:
                i, j = subset
                value = matrix[i][i]*matrix[j][j] - matrix[i][j]*matrix[j][i]
            else:
                value = determinant_3(matrix)
            result.append(value)
    return result


def minkowski_control(rho, q1, q2, q3):
    """Reconstruct one simplex and its facet normals without P formulas."""
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    q_values = np.array([q1, q2, q3], dtype=float)
    c = (1.0-rho-q_values)/2.0
    root_rho = math.sqrt(rho)
    time_components = -c/root_rho
    old_link_gram = np.full((3, 3), 0.5)
    np.fill_diagonal(old_link_gram, 1.0)
    spatial_gram = old_link_gram + np.outer(time_components, time_components)
    spatial_rows = np.linalg.cholesky(spatial_gram)
    y = np.array([root_rho, 0.0, 0.0, 0.0])
    u_rows = np.column_stack((time_components, spatial_rows))
    edges = np.vstack((y, u_rows))

    expected = np.empty((4, 4), dtype=float)
    expected[0, 0] = -rho
    expected[0, 1:] = c
    expected[1:, 0] = c
    expected[1:, 1:] = old_link_gram
    measured = edges @ eta @ edges.T

    def outward_normal(included, omitted):
        constraints = edges[list(included)] @ eta
        _, _, right_vectors = np.linalg.svd(constraints)
        normal = right_vectors[-1]
        norm_squared = float(normal @ eta @ normal)
        if norm_squared <= 0:
            raise RuntimeError("facet normal is not spacelike")
        normal /= math.sqrt(norm_squared)
        if normal @ eta @ edges[omitted] > 0:
            normal *= -1
        return normal

    normal_u2 = outward_normal((0, 1, 3), 2)
    normal_u3 = outward_normal((0, 1, 2), 3)
    outward_cosine = float(normal_u2 @ eta @ normal_u3)
    measured_theta = math.pi-math.acos(np.clip(outward_cosine, -1.0, 1.0))

    p11, p22, p12 = projected_components(float(rho), q1, q2, q3)
    formula_theta = math.acos(p12/math.sqrt(p11*p22))
    return {
        "q": [q1, q2, q3],
        "gram_residual": float(np.max(np.abs(measured-expected))),
        "angle_residual": float(measured_theta-formula_theta),
        "normal_norm_residual": max(
            abs(float(normal_u2 @ eta @ normal_u2)-1.0),
            abs(float(normal_u3 @ eta @ normal_u3)-1.0),
        ),
    }


print("=" * 78)
print("ASYMMETRIC LORENTZIAN 600-CELL TENT WITNESS")
print("=" * 78)

# -------------------------------------------------------------------------
# Exact one-simplex algebra
# -------------------------------------------------------------------------
rho, c1, c2, c3 = sp.symbols("rho c1 c2 c3", real=True)
gram = sp.Matrix([
    [-rho, c1, c2, c3],
    [c1, 1, sp.Rational(1, 2), sp.Rational(1, 2)],
    [c2, sp.Rational(1, 2), 1, sp.Rational(1, 2)],
    [c3, sp.Rational(1, 2), sp.Rational(1, 2), 1],
])
Q = (
    c1**2+c2**2+c3**2
    +(c1-c2)**2+(c1-c3)**2+(c2-c3)**2+2*rho
)
check(
    "the asymmetric four-simplex determinant is the preregistered sum of squares",
    sp.simplify(gram.det()+Q/4) == 0,
    "det G=-Q/4, with Q>0 for rho>0",
)
link = gram[1:, 1:]
schur = sp.factor(gram.det()/link.det())
check(
    "the positive link block and negative Schur complement give inertia (3+,1-)",
    link.eigenvals() == {sp.Rational(1, 2): 2, sp.Integer(2): 1}
    and sp.simplify(schur+Q/2) == 0,
    f"Schur complement={schur}",
)

hinge = gram.extract([0, 1], [0, 1])
others = gram.extract([0, 1], [2, 3])
raw_others = gram.extract([2, 3], [2, 3])
projected = sp.simplify(raw_others-others.T*hinge.inv()*others)
denominator = 4*(c1**2+rho)
expected_projected = sp.Matrix([
    [(4*c1**2-4*c1*c2+4*c2**2+3*rho)/denominator,
     (2*c1**2-2*c1*c2-2*c1*c3+4*c2*c3+rho)/denominator],
    [(2*c1**2-2*c1*c2-2*c1*c3+4*c2*c3+rho)/denominator,
     (4*c1**2-4*c1*c3+4*c3**2+3*rho)/denominator],
])
check(
    "orthogonal projection derives the asymmetric dihedral-angle formula exactly",
    all(sp.simplify(projected[i, j]-expected_projected[i, j]) == 0
        for i in range(2) for j in range(2)),
)

q_symbol = sp.symbols("q", positive=True, real=True)
area_radicand = 4*rho+(1-q_symbol-rho)**2
area = sp.sqrt(area_radicand)/4
weight = sp.factor(sp.diff(area, rho))
expected_weight = (1+q_symbol+rho)/(4*sp.sqrt(area_radicand))
check(
    "the timelike hinge area and positive pole weight are exact",
    sp.simplify(weight-expected_weight) == 0,
    "A/a^2=sqrt(4rho+(1-q-rho)^2)/4 and dA/drho>0",
)

# -------------------------------------------------------------------------
# Independent combinatorial icosahedron and symmetry audit
# -------------------------------------------------------------------------
graph = nx.icosahedral_graph()
triangles = sorted(
    tuple(sorted(clique))
    for clique in nx.enumerate_all_cliques(graph)
    if len(clique) == 3
)
check(
    "the link is the combinatorial icosahedron",
    graph.number_of_nodes() == 12
    and graph.number_of_edges() == 30
    and len(triangles) == 20
    and all(graph.degree[node] == 5 for node in graph),
    "f-vector of the link: 12 vertices, 30 edges, 20 faces",
)
face_incidence = {node: 0 for node in graph}
for triangle in triangles:
    for node in triangle:
        face_incidence[node] += 1
check(
    "every internal tent hinge belongs to exactly five four-simplices",
    set(face_incidence.values()) == {5},
)

automorphisms = list(nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter())
u0 = 0
stabilizer = [mapping for mapping in automorphisms if mapping[u0] == u0]
distances = nx.single_source_shortest_path_length(graph, u0)

def group_orbits(group):
    remaining = set(graph.nodes)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = {mapping[seed] for mapping in group}
        orbits.append(orbit)
        remaining -= orbit
    return sorted(orbits, key=lambda orbit: (len(orbit), min(orbit)))


full_orbits = group_orbits(automorphisms)
stabilizer_orbits = group_orbits(stabilizer)
shells = [set(node for node in graph if distances[node] == distance) for distance in range(4)]
check(
    "the full icosahedral graph automorphism group is transitive of order 120",
    len(automorphisms) == 120 and len(full_orbits) == 1,
)
check(
    "a selected vertex has order-ten stabilizer and four shells 1,5,5,1",
    len(stabilizer) == 10
    and sorted(stabilizer_orbits, key=lambda orbit: min(distances[v] for v in orbit)) == shells
    and [len(shell) for shell in shells] == [1, 5, 5, 1],
)
triangle_shell_types = {}
for triangle in triangles:
    shell_type = tuple(sorted(distances[node] for node in triangle))
    triangle_shell_types[shell_type] = triangle_shell_types.get(shell_type, 0)+1
expected_shell_types = {
    (0, 1, 1): 5,
    (1, 1, 2): 5,
    (1, 2, 2): 5,
    (2, 2, 3): 5,
}
check(
    "the 20 link faces split into the four preregistered shell types",
    triangle_shell_types == expected_shell_types,
    str(triangle_shell_types),
)

# -------------------------------------------------------------------------
# Arb interval certificate for the disclosed witness
# -------------------------------------------------------------------------
ctx.prec = 200
rho_ball = arb(1)/4
q1_ball = arb(3)/2
q2_ball = arb(4)/5
q3_ball = arb(3)/2
x_lo = arb(11)/25
x_hi = arb(9)/20
x_interval = arb(arb(89)/200, arb(1)/200)


def q_by_node_dual(x):
    shell_values = [x, Dual(q1_ball), Dual(q2_ball), Dual(q3_ball)]
    return {node: shell_values[distances[node]] for node in graph}


def global_E_dual(x):
    rho_dual = Dual(rho_ball)
    q_values = q_by_node_dual(x)
    deficits = {node: Dual(2*arb.pi()) for node in graph}
    for triangle in triangles:
        for node in triangle:
            other = [vertex for vertex in triangle if vertex != node]
            deficits[node] = deficits[node] - theta_dual(
                rho_dual, q_values[node], q_values[other[0]], q_values[other[1]]
            )
    total = Dual(0)
    for node in graph:
        total += deficits[node]*weight_dual(rho_dual, q_values[node])
    return total, deficits


lo_eval, _ = global_E_dual(Dual(x_lo, 1))
hi_eval, _ = global_E_dual(Dual(x_hi, 1))
interval_eval, interval_deficits = global_E_dual(Dual(x_interval, 1))
check(
    "Arb certifies opposite endpoint signs for the pole equation",
    lo_eval.v < 0 and hi_eval.v > 0,
    f"E(11/25)={lo_eval.v}; E(9/20)={hi_eval.v}",
)
check(
    "Arb certifies strict monotonicity throughout the complete root bracket",
    interval_eval.d > 0,
    f"dE/dx={interval_eval.d}",
)

# Causal and real-angle domains are checked on the whole x interval, not just
# at the numerical root.  Use every actual face so this also audits the shell
# assignment and incidence loop.
q_shell_balls = [x_interval, q1_ball, q2_ball, q3_ball]
all_minors = []
angle_branch_determinants = []
interval_angles = []
for triangle in triangles:
    q_values = [q_shell_balls[distances[node]] for node in triangle]
    all_minors.extend(principal_minors(final_tetra_gram(q_values)))
    for index in range(3):
        others_indices = [j for j in range(3) if j != index]
        p11, p22, p12 = projected_components(
            rho_ball,
            q_values[index],
            q_values[others_indices[0]],
            q_values[others_indices[1]],
        )
        angle_branch_determinants.append(p11*p22-p12*p12)
        interval_angles.append((p12/(p11*p22).sqrt()).acos())

check(
    "all final tetrahedra stay strictly spacelike on the complete bracket",
    all(value > 0 for value in all_minors),
    f"smallest lower bound={min(float(value.lower()) for value in all_minors):.6g}",
)
check(
    "all 60 dihedral angles remain on a real nonsingular branch",
    all(value > 0 for value in angle_branch_determinants),
    "every projected 2x2 Gram determinant is strictly positive",
)

all_q_balls = [x_interval, q1_ball, q2_ball, q3_ball]
area_radicands = [4*rho_ball+(1-q-rho_ball)**2 for q in all_q_balls]
weights = [(1+q+rho_ball)/(4*radicand.sqrt())
           for q, radicand in zip(all_q_balls, area_radicands)]
check(
    "all boundary lengths, timelike hinge areas and pole weights stay admissible",
    all(q > 0 for q in all_q_balls)
    and all(radicand > 0 for radicand in area_radicands)
    and all(value > 0 for value in weights),
)
check(
    "the witness contains an admissible angle strictly above the symmetric 72-degree bound",
    any(angle > 2*arb.pi()/5 for angle in interval_angles),
    f"largest certified interval lower endpoint="
    f"{max(float(angle.lower()) for angle in interval_angles)*180/math.pi:.6f} deg",
)

# -------------------------------------------------------------------------
# High-precision root and independent Minkowski controls
# -------------------------------------------------------------------------
mp.mp.dps = 90


def theta_mp(rho_value, qa, qb, qc):
    c_a = (1-rho_value-qa)/2
    c_b = (1-rho_value-qb)/2
    c_c = (1-rho_value-qc)/2
    den = 4*(c_a*c_a+rho_value)
    p11 = (4*c_a*c_a-4*c_a*c_b+4*c_b*c_b+3*rho_value)/den
    p22 = (4*c_a*c_a-4*c_a*c_c+4*c_c*c_c+3*rho_value)/den
    p12 = (2*c_a*c_a-2*c_a*c_b-2*c_a*c_c+4*c_b*c_c+rho_value)/den
    return mp.acos(p12/mp.sqrt(p11*p22))


def global_E_mp(x):
    rho_value = mp.mpf(1)/4
    shell_values = [x, mp.mpf(3)/2, mp.mpf(4)/5, mp.mpf(3)/2]
    q_values = {node: shell_values[distances[node]] for node in graph}
    deficits = {node: 2*mp.pi for node in graph}
    for triangle in triangles:
        for node in triangle:
            other = [vertex for vertex in triangle if vertex != node]
            deficits[node] -= theta_mp(
                rho_value, q_values[node], q_values[other[0]], q_values[other[1]]
            )
    weight_values = {
        node: (1+q_values[node]+rho_value)
        /(4*mp.sqrt(4*rho_value+(1-q_values[node]-rho_value)**2))
        for node in graph
    }
    total = mp.fsum(deficits[node]*weight_values[node] for node in graph)
    return total, deficits, weight_values


root_lower = mp.mpf(11)/25
root_upper = mp.mpf(9)/20
for _ in range(280):
    midpoint = (root_lower+root_upper)/2
    if global_E_mp(midpoint)[0] < 0:
        root_lower = midpoint
    else:
        root_upper = midpoint
root = (root_lower+root_upper)/2
root_residual, root_deficits, root_weights = global_E_mp(root)
check(
    "high-precision bisection isolates the unique preregistered root",
    root_lower < root < root_upper
    and abs(root_residual) < mp.mpf("1e-80")
    and abs(root-mp.mpf("0.44333089835748125745")) < mp.mpf("1e-19"),
    f"x*={mp.nstr(root, 60)}, |E|={mp.nstr(abs(root_residual), 5)}",
)

def orbit_value(values, shell):
    entries = [values[node] for node in graph if distances[node] == shell]
    if max(entries)-min(entries) > mp.mpf("1e-70"):
        raise RuntimeError("stabilizer orbit values are not constant")
    return entries[0]


deficit_orbits = [orbit_value(root_deficits, shell) for shell in range(4)]
weight_orbits = [orbit_value(root_weights, shell) for shell in range(4)]
check(
    "the pole is a weighted cancellation, not twelve locally flat hinges",
    deficit_orbits[0] > 0
    and deficit_orbits[1] < 0
    and deficit_orbits[2] > 0
    and deficit_orbits[3] > 0
    and max(abs(value) for value in deficit_orbits) > mp.mpf("0.1"),
    "deficit orbits=" + str([mp.nstr(value, 16) for value in deficit_orbits]),
)

root_float = float(root)
shell_float = [root_float, 1.5, 0.8, 1.5]
ordered_types = sorted({
    (
        shell_float[distances[node]],
        *sorted(shell_float[distances[other]] for other in triangle if other != node),
    )
    for triangle in triangles for node in triangle
})
controls = [minkowski_control(0.25, *ordered_type) for ordered_type in ordered_types]
check(
    "independent Minkowski coordinates reconstruct every ordered shell simplex",
    max(control["gram_residual"] for control in controls) < 2e-14
    and max(control["normal_norm_residual"] for control in controls) < 2e-14,
    f"{len(controls)} ordered types; max Gram residual="
    f"{max(control['gram_residual'] for control in controls):.3e}",
)
check(
    "independent facet normals reproduce every projected-angle value",
    max(abs(control["angle_residual"]) for control in controls) < 3e-14,
    f"max angle residual={max(abs(control['angle_residual']) for control in controls):.3e}",
)

# Selection is a separate, deliberately negative gate.
check(
    "bare symmetry has one length orbit whereas the witness needs a chosen direction",
    len(full_orbits) == 1 and len(stabilizer_orbits) == 4,
    "full-invariant dimension=1; selected-vertex-stabilizer dimension=4",
)
check(
    "one pole equation cannot select the disclosed four shell lengths and pole",
    interval_eval.d > 0 and len(stabilizer_orbits) == 4,
    "the derivative proves local root regularity, but the boundary ratios were fitted",
)

result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "scope": {
        "carrier": "[v,v']*L_v with L_v the combinatorial icosahedron",
        "action": "ordinary real Lorentzian Regge action on timelike hinges",
        "volume_coefficient": 0,
        "variation": "pole only; boundary q values held fixed",
    },
    "combinatorics": {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "faces": len(triangles),
        "faces_per_internal_hinge": 5,
        "automorphism_group_order": len(automorphisms),
        "selected_vertex_stabilizer_order": len(stabilizer),
        "shell_sizes": [len(shell) for shell in shells],
        "triangle_shell_types": {str(key): value for key, value in triangle_shell_types.items()},
    },
    "frozen_witness": {
        "rho": "1/4",
        "q_shells": ["x", "3/2", "4/5", "3/2"],
        "root_bracket": ["11/25", "9/20"],
        "endpoint_E": [str(lo_eval.v), str(hi_eval.v)],
        "derivative_on_bracket": str(interval_eval.d),
        "root": mp.nstr(root, 80),
        "root_enclosure_width": mp.nstr(root_upper-root_lower, 8),
        "root_residual": mp.nstr(root_residual, 8),
        "deficit_orbits": [mp.nstr(value, 40) for value in deficit_orbits],
        "weight_orbits": [mp.nstr(value, 40) for value in weight_orbits],
    },
    "selection_audit": {
        "full_automorphism_invariant_length_dimension": 1,
        "selected_vertex_stabilizer_length_dimension": 4,
        "pole_equations": 1,
        "status": "STRUCTURAL/FITTED NON-SELECTION",
        "reason": "the chosen link direction and three frozen shell ratios are not selected",
    },
    "verdict": {
        "asymmetric_vacuum_existence": "DERIVED",
        "universal_symmetric_angle_bound": "REFUTED outside the symmetric ansatz",
        "physical_clock": "OPEN",
    },
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(result, indent=2)+"\n")

print("-" * 78)
print(f"Result: {passed}/{tests} checks passed")
print(f"Wrote {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
