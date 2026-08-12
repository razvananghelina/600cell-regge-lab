#!/usr/bin/env python3
"""Exact round-S3 audit of the de Rham heat response on Hopf spin two.

Protocol commit 31ecea7 froze the hypotheses, observable and scope.  This is
a post-recognition structural audit, not a blind physical discovery.
"""

from itertools import product
from math import comb
import json
from pathlib import Path

import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "hopf_kahler_induced_gravity.json"
PROTOCOL_COMMIT = "31ecea7"

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def zero(value):
    return sy.simplify(sy.expand(value)) == 0


def zero_matrix(matrix):
    return all(zero(value) for value in matrix)


print("Hopf spin-two response of the continuum de Rham heat coefficient")

# A generic left-invariant metric on SU(2), in a fixed Lie-algebra frame.
a, b, c, d, e, f = sy.symbols("a b c d e f", real=True)
variables = (a, b, c, d, e, f)
metric = sy.Matrix(((a, d, e), (d, b, f), (e, f, c)))
metric_inverse = metric.inv()

# [e_i,e_j] = C_ij^k e_k = 2 epsilon_ijk e_k.
structure = [
    [
        [2 * sy.LeviCivita(i, j, k) for k in range(3)]
        for j in range(3)
    ]
    for i in range(3)
]

# Koszul formula for constant left-invariant vector fields:
# 2 <nabla_i e_j,e_k> = <[i,j],k>-<[j,k],i>+<[k,i],j>.
connection = [[[] for _ in range(3)] for _ in range(3)]
for i in range(3):
    for j in range(3):
        lowered = sy.Matrix([
            sy.Rational(1, 2) * sum(
                structure[i][j][m] * metric[m, k]
                - structure[j][k][m] * metric[m, i]
                + structure[k][i][m] * metric[m, j]
                for m in range(3)
            )
            for k in range(3)
        ])
        connection[i][j] = [
            sy.factor(value) for value in metric_inverse * lowered
        ]

torsion_controls = []
metric_controls = []
for i in range(3):
    gamma_i = sy.Matrix(3, 3, lambda row, col: connection[i][col][row])
    metric_controls.append(zero_matrix(gamma_i.T * metric + metric * gamma_i))
    for j in range(3):
        torsion_controls.extend(
            zero(
                connection[i][j][k]
                - connection[j][i][k]
                - structure[i][j][k]
            )
            for k in range(3)
        )
check(
    "Koszul connection is metric compatible and torsion free",
    all(metric_controls) and all(torsion_controls),
)

# R(e_i,e_j)e_k = nabla_i nabla_j e_k - nabla_j nabla_i e_k
#                  - nabla_[e_i,e_j] e_k.
curvature = [[[[sy.S.Zero] * 3 for _ in range(3)]
              for _ in range(3)] for _ in range(3)]
for i in range(3):
    for j in range(3):
        for k in range(3):
            for ell in range(3):
                curvature[i][j][k][ell] = sy.factor(sum(
                    connection[j][k][m] * connection[i][m][ell]
                    - connection[i][k][m] * connection[j][m][ell]
                    - structure[i][j][m] * connection[m][k][ell]
                    for m in range(3)
                ))

ricci = sy.Matrix(3, 3, lambda j, k: sy.factor(sum(
    curvature[i][j][k][i] for i in range(3)
)))
scalar_curvature = sy.factor(sum(
    metric_inverse[j, k] * ricci[j, k]
    for j in range(3) for k in range(3)
))

closed_scalar = sy.factor(
    2 * (sy.trace(metric) ** 2 - 2 * sy.trace(metric * metric))
    / metric.det()
)
check(
    "direct Koszul curvature gives the exact SU(2) scalar formula",
    zero(scalar_curvature - closed_scalar),
    "R(G)=2((Tr G)^2-2 Tr(G^2))/det(G)",
)

round_point = {a: 1, b: 1, c: 1, d: 0, e: 0, f: 0}
round_ricci = ricci.subs(round_point)
round_scalar = sy.simplify(scalar_curvature.subs(round_point))
check(
    "round normalization is Ric=2g and R=6",
    round_ricci == 2 * sy.eye(3) and round_scalar == 6,
)

# In dimension three, the scale-invariant fixed-volume Einstein functional
# per fixed coordinate-volume factor is Y=R det(G)^(1/3).
yamabe_density = sy.simplify(scalar_curvature * metric.det() ** sy.Rational(1, 3))
gradient = sy.Matrix([
    sy.simplify(sy.diff(yamabe_density, variable).subs(round_point))
    for variable in variables
])
hessian = sy.Matrix([
    [
        sy.simplify(
            sy.diff(yamabe_density, left, right).subs(round_point)
        )
        for right in variables
    ]
    for left in variables
])
check(
    "the unit round metric is stationary for the normalized functional",
    gradient == sy.zeros(6, 1),
)

sqrt2 = sy.sqrt(2)
sqrt3 = sy.sqrt(3)
sqrt6 = sy.sqrt(6)
# Coordinates are (G11,G22,G33,G12,G13,G23).  These five matrices are an
# orthonormal Frobenius basis of Sym^2_0(R^3).
tracefree_coordinates = (
    sy.Matrix((1 / sqrt2, -1 / sqrt2, 0, 0, 0, 0)),
    sy.Matrix((1 / sqrt6, 1 / sqrt6, -2 / sqrt6, 0, 0, 0)),
    sy.Matrix((0, 0, 0, 1 / sqrt2, 0, 0)),
    sy.Matrix((0, 0, 0, 0, 1 / sqrt2, 0)),
    sy.Matrix((0, 0, 0, 0, 0, 1 / sqrt2)),
)
scale_coordinate = sy.Matrix((1 / sqrt3, 1 / sqrt3, 1 / sqrt3, 0, 0, 0))
restricted_hessian = sy.Matrix([
    [sy.simplify(left.dot(hessian * right))
     for right in tracefree_coordinates]
    for left in tracefree_coordinates
])
check(
    "normalized curvature Hessian is -4 times the STF norm",
    restricted_hessian == -4 * sy.eye(5),
    "delta^2 Y(H,H)=-4 Tr(H^2)",
)
check(
    "the removed overall scale direction is exactly null",
    zero_matrix(hessian * scale_coordinate),
)
check(
    "the curvature response has exact STF inertia (0 positive,5 negative,0 zero)",
    restricted_hessian.rank() == 5
    and restricted_hessian.eigenvals() == {-4: 5},
)

# Connect this full STF space back to the six fivefold Hopf axes.
sqrt5 = sy.sqrt(5)
phi = (1 + sqrt5) / 2
vertices = []
for first, second in product((1, -1), repeat=2):
    vertices.extend((
        sy.Matrix((0, first, second * phi)),
        sy.Matrix((first, second * phi, 0)),
        sy.Matrix((first * phi, 0, second)),
    ))


def projector(vector):
    return sy.simplify(vector * vector.T / vector.dot(vector))


projectors = {}
for vertex in vertices:
    candidate = projector(vertex)
    key = tuple(sy.radsimp(value) for value in candidate)
    projectors.setdefault(key, candidate)
centered = tuple(
    sy.simplify(candidate - sy.eye(3) / 3)
    for candidate in projectors.values()
)
hopf_columns = [sy.Matrix(tensor).reshape(9, 1) for tensor in centered]
hopf_responses = [
    sy.simplify(-4 * sy.trace(tensor * tensor)) for tensor in centered
]
check(
    "the six Hopf tensors span every tested STF direction",
    len(centered) == 6 and sy.Matrix.hstack(*hopf_columns).rank() == 5,
)
check(
    "all six Hopf rays have the same nonzero curvature Hessian",
    hopf_responses == [-sy.Rational(8, 3)] * 6,
    "delta^2Y(T_i,T_i)=-8/3",
)

# Heat coefficient for the Hodge Laplacian.  In the convention
# Delta_p=-(nabla^2+E_p), E_p=-q_p(R), the local a2 density is
# tr(E_p+R I/6).  The Weitzenbock trace is
# tr q_p = binomial(n-2,p-1) R.
dimension = 3
degree_ranks = [comb(dimension, degree) for degree in range(dimension + 1)]
weitzenbock_multipliers = [
    comb(dimension - 2, degree - 1)
    if 1 <= degree <= dimension - 1 else 0
    for degree in range(dimension + 1)
]
heat_a2_multipliers = [
    sy.Rational(rank, 6) - q_multiplier
    for rank, q_multiplier in zip(degree_ranks, weitzenbock_multipliers)
]
check(
    "degreewise de Rham a2 curvature multipliers are exact",
    degree_ranks == [1, 3, 3, 1]
    and weitzenbock_multipliers == [0, 1, 1, 0]
    and heat_a2_multipliers == [
        sy.Rational(1, 6), -sy.Rational(1, 2),
        -sy.Rational(1, 2), sy.Rational(1, 6),
    ],
)
ordinary_heat_multiplier = sy.simplify(sum(heat_a2_multipliers))
supertrace_multiplier = sy.simplify(sum(
    (-1) ** degree * heat_a2_multipliers[degree]
    for degree in range(4)
))
check(
    "ordinary full-exterior heat trace has curvature coefficient -2R/3",
    ordinary_heat_multiplier == -sy.Rational(2, 3),
)
check(
    "the graded supertrace curvature coefficient cancels",
    supertrace_multiplier == 0,
    "ordinary trace is load-bearing; the index supertrace sees no metric response",
)

heat_restricted_hessian = sy.simplify(
    ordinary_heat_multiplier * restricted_hessian
)
check(
    "the normalized heat a2 response is positive definite on Hopf spin two",
    heat_restricted_hessian == sy.Rational(8, 3) * sy.eye(5),
    "delta^2 a2(H,H)=(8/3) Tr(H^2), up to universal volume/(4pi)^(3/2)",
)
check(
    "the heat-coefficient response has exact inertia (5 positive,0 zero,0 negative)",
    heat_restricted_hessian.rank() == 5
    and heat_restricted_hessian.eigenvals() == {sy.Rational(8, 3): 5},
)

# Direct path control: G(s)=exp(sH) through second order has determinant one
# through second order for tracefree H.  This verifies that multiplying the
# normalized Hessian above is genuinely the fixed-volume integrated heat
# response, not an accidental normalization identity.
s = sy.symbols("s", real=True)
h0, h1, h2, h3, h4 = sy.symbols("h0:5", real=True)
generic_hopf_tensor = sy.Matrix((
    (h0, h2, h3),
    (h2, h1, h4),
    (h3, h4, -h0 - h1),
))
metric_path = sy.eye(3) + s * generic_hopf_tensor + sy.Rational(1, 2) * s**2 * generic_hopf_tensor**2
path_scalar = sy.factor(
    2 * (sy.trace(metric_path) ** 2 - 2 * sy.trace(metric_path**2))
    / metric_path.det()
)
path_volume = sy.sqrt(metric_path.det())
path_heat_a2 = ordinary_heat_multiplier * path_scalar * path_volume
path_volume_second = sy.simplify(sy.diff(path_volume, s, 2).subs(s, 0))
path_heat_second = sy.simplify(sy.diff(path_heat_a2, s, 2).subs(s, 0))
check(
    "an explicit determinant-one Hopf metric path gives the same heat Hessian",
    zero(path_volume_second)
    and zero(
        path_heat_second
        - sy.Rational(8, 3) * sy.trace(generic_hopf_tensor**2)
    ),
    "G(s)=exp(sH)+O(s^3), delta^2 A2=(8/3) Tr(H^2)",
)

# Source-scope audit.  These checks deliberately concern the authoritative
# constructors, not every historical exploratory file.
kahler_source = (HERE / "verify_kahler_dirac.py").read_text()
whitney_source = (HERE / "verify_whitney_kahler_induction.py").read_text()
spectral_source = (HERE / "verify_spectral_action.py").read_text()
check(
    "the certified finite incidence D has no continuous metric parameter",
    "D = sp.bmat" in kahler_source
    and "local_whitney_mass" not in kahler_source
    and "metric_parameter" not in kahler_source,
)
check(
    "the Whitney D is metric-derived but evaluated on fixed vertex data",
    "def local_whitney_mass(points, degree):" in whitney_source
    and "mass_coarse = [assemble_mass(reference_vertices" in whitney_source
    and "metric_parameter" not in whitney_source,
)
check(
    "the finite spectral certificate computes moments, not a variable-metric action",
    "c0 = N_total" in spectral_source
    and "c1 = np.sum(all_evals_D2)" in spectral_source
    and "c2 = 0.5 * np.sum(all_evals_D2**2)" in spectral_source
    and "metric_parameter" not in spectral_source,
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": "post-recognition hostile structural audit",
    "phenomenological_target_used": False,
    "hypotheses": {
        "carrier": "unit round S3=SU(2)",
        "metric_family": "positive left-invariant Gram matrices G",
        "perturbations": "homogeneous H in Sym^2_0(R^3)",
        "operator": "continuum de Rham Kahler-Dirac D_g=d+d*_g",
        "observable": "ordinary heat trace a2 curvature coefficient",
        "not_assumed": [
            "Einstein equation", "Regge action", "Lorentzian time",
            "diffeomorphism quotient", "source coupling", "Newton scale",
        ],
    },
    "exact_results": {
        "scalar_curvature": "2((Tr G)^2-2Tr(G^2))/det(G)",
        "round_ricci": "2g",
        "round_scalar_curvature": 6,
        "normalized_curvature_hessian_on_stf": "-4 I_5",
        "ordinary_derham_heat_a2_curvature_multiplier": "-2/3",
        "supertrace_a2_curvature_multiplier": "0",
        "normalized_heat_a2_hessian_on_stf": "8/3 I_5",
        "heat_hessian_inertia": [5, 0, 0],
        "hopf_span_rank": 5,
    },
    "verdicts": [
        {
            "label": "DERIVED CONTINUUM KINEMATICS",
            "claim": "the exact round-S3 curvature and de Rham heat identities",
        },
        {
            "label": "STRUCTURAL INDUCED-GRAVITY ADVANCE",
            "claim": (
                "the continuum operator has a nondegenerate curvature response "
                "on every homogeneous Hopf spin-two direction"
            ),
        },
        {
            "label": "OPEN PHYSICAL SELECTION",
            "claim": (
                "the finite/refined theory does not yet select variable metric, "
                "spectral functional, normalization, Lorentzian gauge dynamics "
                "or universal source coupling"
            ),
        },
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("machine-readable induced-gravity certificate was written", OUTPUT.exists())

print(f"\nRESULT: {passed}/{tests} checks passed")
print("DERIVED: the de Rham heat a2 coefficient sees all five Hopf TT directions.")
print("STRUCTURAL: this is an induced-curvature mechanism, not selected GR.")
print("OPEN: variable metric, spectral functional, time/gauge/source and scale.")
raise SystemExit(0 if passed == tests else 1)
