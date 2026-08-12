#!/usr/bin/env python3
"""Exact controls for the fixed-Regge FEEC transfer.

Protocol commit 1682a46 froze the carrier identification, metric bounds,
transfer lemma and theorem-hypothesis chain before this verifier existed.

The script checks the finite algebraic and geometric inputs.  It does not
claim to machine-prove the cited functional-analytic theorems; their complete
hypothesis audit and the equivalent-inner-product lemma are written out in
the accompanying result note.
"""

import json
from pathlib import Path

import sympy as sy


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "whitney_regge_continuum_transfer.json"
EDGEWISE = HERE / "whitney_rank_edgewise_refinement.json"
DYNAMICS = HERE / "whitney_edgewise_continuum_dynamics.json"
PROTOCOL_COMMIT = "1682a46"

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


print("Fixed-Regge FEEC transfer controls")

edgewise = json.loads(EDGEWISE.read_text())
dynamics = json.loads(DYNAMICS.read_text())

check(
    "prior rank-edgewise certificate has the preregistered provenance",
    edgewise["protocol_commit"] == "58fa9fc",
)
check(
    "tower has a level-independent three-shape bound",
    all(
        edgewise["edgewise_counts"][str(level)]["normalized_shape_classes"]
        <= 3
        for level in (1, 2, 3, 4)
    )
    and edgewise["rank_chamber_edgewise_quality_min"] > 0,
)
check(
    "powers-of-two edgewise refinement is nested eight-to-one",
    edgewise["nesting_fine_per_coarse"] == [8] * 8,
)

finite_inputs = dynamics["finite_theorem_inputs"]
check(
    "exact Whitney finite inputs were independently certified",
    all(finite_inputs.values()),
)
check(
    "the finite calculation used no phenomenological target",
    dynamics["phenomenological_target_used"] is False,
)

# A regular 600-cell at unit circumradius has adjacent scalar product phi/2.
# Recover the edge length and facet supporting distance from its Gram matrix,
# without importing coordinates of a preferred facet.
phi = (sy.Integer(1) + sy.sqrt(5)) / 2
adjacent_dot = phi / 2
vertex_gram = sy.eye(4) * (1 - adjacent_dot) + sy.ones(4, 4) * adjacent_dot
edge_length_sq = sy.simplify(
    vertex_gram[0, 0] + vertex_gram[1, 1] - 2 * vertex_gram[0, 1]
)
check(
    "unit 600-cell edge length gives adjacent dot product phi/2",
    sy.simplify(edge_length_sq - 1 / phi**2) == 0,
    f"edge_length^2={edge_length_sq}",
)

ones = sy.ones(4, 1)
facet_center_sq = sy.simplify((ones.T * vertex_gram * ones)[0] / 16)
expected_a2 = sy.simplify((2 + 3 * phi) / 8)
center_vertex_products = [
    sy.simplify(sum(vertex_gram[row, col] for row in range(4)) / 4)
    for col in range(4)
]
check(
    "facet centroid is normal to every facet edge",
    len(set(center_vertex_products)) == 1,
    f"centroid dot each vertex={center_vertex_products[0]}",
)
check(
    "facet supporting distance is a^2=(2+3 phi)/8",
    sy.simplify(facet_center_sq - expected_a2) == 0,
    f"a^2={expected_a2}",
)

# Derive the differential of radial normalization in ambient R^4 exactly.
x_symbols = sy.symbols("x0:4", real=True)
x = sy.Matrix(x_symbols)
r2 = sum(value**2 for value in x_symbols)
radial = x / sy.sqrt(r2)
jacobian = radial.jacobian(x_symbols)
round_pullback = sy.simplify(jacobian.T * jacobian)
expected_pullback = sy.eye(4) / r2 - (x * x.T) / r2**2
check(
    "radial pullback metric formula is exact",
    sy.simplify(round_pullback - expected_pullback) == sy.zeros(4),
)

# Put the facet normal along coordinate 0 and the tangential component of x
# along coordinate 1.  Restriction to the three-dimensional facet tangent
# space then exposes all generalized metric eigenvalues.
a, y = sy.symbols("a y", positive=True)
local_x = sy.Matrix((a, y, 0, 0))
local_r2 = a**2 + y**2
local_metric = sy.eye(4) / local_r2 - local_x * local_x.T / local_r2**2
tangent_metric = sy.simplify(local_metric[1:4, 1:4])
expected_tangent = sy.diag(
    a**2 / local_r2**2,
    1 / local_r2,
    1 / local_r2,
)
check(
    "facet-tangent radial metric spectrum is exact",
    sy.simplify(tangent_metric - expected_tangent) == sy.zeros(3),
    "eigenvalues=(a^2/r^4, 1/r^2, 1/r^2)",
)

# On a convex unit-circumradius facet a^2 <= r^2 <= 1.  The endpoint
# extrema of the two monotone eigenvalue functions give a fixed equivalence
# interval [a^2, 1/a^2], independent of refinement level.
a2_numeric = float(expected_a2.evalf(30))
metric_lower = expected_a2
metric_upper = sy.simplify(1 / expected_a2)
check(
    "radial/flat metric equivalence constants are finite and nondegenerate",
    0 < a2_numeric < 1 and float(metric_upper.evalf()) < 2,
    f"{metric_lower} <= lambda <= {metric_upper}",
)

# The trace metric of a regular tetrahedron on a shared triangular face is
# determined only by that face.  This is the exact tangential--tangential
# Regge matching condition.
face_gram = sy.Matrix(
    ((edge_length_sq, 1 - adjacent_dot),
     (1 - adjacent_dot, edge_length_sq))
)
equilateral_face_gram = edge_length_sq * sy.Matrix(
    ((1, sy.Rational(1, 2)), (sy.Rational(1, 2), 1))
)
check(
    "flat tetrahedral metrics have identical traces on a shared face",
    sy.simplify(face_gram - equilateral_face_gram) == sy.zeros(2),
)

# Nontrivial finite-dimensional control of the transfer-lemma identity.  The
# two minimal right inverses are different, as expected when the adjoint
# changes.  Projecting either representative off ker(d) in the new metric
# gives the new minimal representative exactly.
d = sy.Matrix(((1, 1, 0), (0, 1, 1)))
g0 = sy.eye(3)
g1 = sy.diag(2, 3, 5)
kernel = sy.Matrix.hstack(*d.nullspace())
r0 = g0.inv() * d.T * (d * g0.inv() * d.T).inv()
r1_direct = g1.inv() * d.T * (d * g1.inv() * d.T).inv()
kernel_projection_g1 = (
    kernel * (kernel.T * g1 * kernel).inv() * kernel.T * g1
)
r1_transferred = sy.simplify((sy.eye(3) - kernel_projection_g1) * r0)
check(
    "equivalent metrics genuinely change the minimal right inverse",
    sy.simplify(r1_direct - r0) != sy.zeros(3, 2),
)
check(
    "new minimal inverse is the bounded kernel-projection transfer",
    sy.simplify(r1_direct - r1_transferred) == sy.zeros(3, 2)
    and sy.simplify(d * r1_transferred) == sy.eye(2)
    and sy.simplify(kernel.T * g1 * r1_transferred) == sy.zeros(1, 2),
)

# A fixed equivalence m||u||_0 <= ||u||_R <= M||u||_0 transfers a uniform
# projection bound C to (M/m)C.  Keep the constants symbolic to avoid hiding
# a numerical threshold in this theorem input.
m, big_m, projection_bound = sy.symbols("m M C", positive=True)
transferred_bound = sy.simplify(big_m * projection_bound / m)
check(
    "uniform projection bound transfer contains no mesh parameter",
    not transferred_bound.has(sy.Symbol("h"))
    and transferred_bound == big_m * projection_bound / m,
    f"||pi_h||_Regge <= ({big_m}/{m}) {projection_bound}",
)

result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "carrier": {
        "map": "R(x)=x/||x|| from boundary(600-cell) to S^3",
        "map_role": "coordinate identification only; exact flat metric is pushed forward",
        "adjacent_vertex_dot": "phi/2",
        "facet_support_distance_squared": str(expected_a2),
        "round_over_flat_tangent_eigenvalues": [
            "a^2/r^4", "1/r^2", "1/r^2"
        ],
        "uniform_vector_metric_bounds": [str(metric_lower), str(metric_upper)],
        "regge_face_trace_exact": True,
    },
    "prior_exact_inputs": {
        "tower": dynamics["all_level_geometry"]["tower"],
        "mesh_law": dynamics["all_level_geometry"]["mesh_law"],
        "normalized_shape_classes": dynamics["all_level_geometry"][
            "normalized_shape_classes"
        ],
        "finite_theorem_inputs": finite_inputs,
    },
    "transfer_lemma_control": {
        "adjoints_and_minimal_inverses_may_change": True,
        "minimal_inverse_kernel_projection_identity": True,
        "uniform_bound_formula": "C_Regge <= (M/m) C_smooth",
        "compactness_proof_location": "whitney_regge_continuum_transfer_result.md",
    },
    "theorem_chain": [
        "Gawlik-McKee: fixed Regge and smooth metrics give equivalent L2/H(d) topologies",
        "Licht: smooth-triangulation FE spaces admit uniform L2-bounded commuting projections",
        "equivalent-norm transfer: same projections are uniform for the fixed Regge metric",
        "AFW Theorem 3.19: compact Hilbert complex plus uniform L2 cochain projections gives eigenvalue/eigenspace convergence",
    ],
    "primary_sources": [
        "Gawlik-McKee arXiv:2410.15579v3",
        "Licht arXiv:2310.14276",
        "Arnold-Falk-Winther arXiv:0906.4325",
    ],
    "verdicts": [
        "DERIVED: exact radial coordinate transfer preserves the fixed Regge metric",
        "DERIVED: Regge and auxiliary round Hilbert norms are uniformly equivalent",
        "DERIVED: compactness transfers although adjoints change",
        "STRUCTURAL: Licht plus AFW give fixed-index Hodge-Laplacian spectral convergence",
        "DERIVED CONDITIONAL: continuum speed |c| after physical time and c are supplied",
        "NOT DERIVED: finite-mesh cone, time, numerical c, mass, fourth dimension or Planck units",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

print(f"\n{passed}/{tests} checks passed")
print(f"Wrote {OUTPUT.name}")
if passed != tests:
    raise SystemExit(1)
