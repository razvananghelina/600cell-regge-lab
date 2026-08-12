#!/usr/bin/env python3
"""Lorentzian Regge gate for the symmetric 600-cell tent carrier.

Protocol commit: 14dd9aa.  The preliminary no-go was disclosed before this
verifier was written.  This script certifies a scoped result under the
ordinary real Lorentzian Regge convention for timelike hinges; it does not
claim that this action or a cosmological coefficient is selected.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_lorentzian_tent.json"
EUCLIDEAN_RESULT = HERE / "gravity_tent_move_regge.json"
PROTOCOL_COMMIT = "14dd9aa"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def minkowski_control(rho_value, q_value):
    """Build one simplex in R^(1,3), then obtain facet normals by nullspaces.

    This does not call the inverse-Gram angle formula being tested.
    Vertices other than the origin are returned as the rows (y,u1,u2,u3).
    """
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    root_rho = np.sqrt(rho_value)
    y = np.array([root_rho, 0.0, 0.0, 0.0])
    y_dot_u = (1.0 - q_value - rho_value) / 2.0
    u_time = -y_dot_u / root_rho

    spatial_gram = np.full((3, 3), 0.5 + u_time**2)
    np.fill_diagonal(spatial_gram, 1.0 + u_time**2)
    spatial_rows = np.linalg.cholesky(spatial_gram)
    u_rows = np.column_stack((np.full(3, u_time), spatial_rows))
    edges = np.vstack((y, u_rows))
    measured_gram = edges @ eta @ edges.T

    expected_gram = np.full((4, 4), 0.5)
    expected_gram[0, :] = y_dot_u
    expected_gram[:, 0] = y_dot_u
    expected_gram[0, 0] = -rho_value
    np.fill_diagonal(expected_gram[1:, 1:], 1.0)

    def outward_normal(included, omitted):
        constraints = edges[list(included)] @ eta
        _, _, right_vectors = np.linalg.svd(constraints)
        normal = right_vectors[-1]
        norm_squared = float(normal @ eta @ normal)
        if norm_squared <= 0:
            raise RuntimeError("facet normal is not spacelike")
        normal /= np.sqrt(norm_squared)
        if normal @ eta @ edges[omitted] > 0:
            normal *= -1
        return normal

    # The internal triangle is span(y,u1).  Its two incident facets in this
    # simplex omit u2 and u3 respectively.
    normal_u2 = outward_normal((0, 1, 3), 2)
    normal_u3 = outward_normal((0, 1, 2), 3)
    outward_cosine = float(normal_u2 @ eta @ normal_u3)
    theta = np.pi - np.arccos(np.clip(outward_cosine, -1.0, 1.0))

    expected_cosine = (
        q_value**2 + 2*q_value*rho_value - 2*q_value
        + rho_value**2 + 1
    ) / (
        2 * (
            q_value**2 + 2*q_value*rho_value - 2*q_value
            + rho_value**2 + rho_value + 1
        )
    )
    return {
        "rho": rho_value,
        "q": q_value,
        "gram_residual": float(np.max(np.abs(measured_gram-expected_gram))),
        "normal_norms": [
            float(normal_u2 @ eta @ normal_u2),
            float(normal_u3 @ eta @ normal_u3),
        ],
        "theta": float(theta),
        "expected_theta": float(np.arccos(expected_cosine)),
        "angle_residual": float(theta-np.arccos(expected_cosine)),
        "gram_eigenvalues": np.linalg.eigvalsh(measured_gram).tolist(),
    }


print("=" * 78)
print("LORENTZIAN 600-CELL REGGE TENT GATE")
print("=" * 78)

with EUCLIDEAN_RESULT.open() as stream:
    carrier_certificate = json.load(stream)
carrier = carrier_certificate["tent_carrier"]
check(
    "the inherited canonical carrier has 12 fivefold internal hinges",
    carrier["internal_triangles"] == 12
    and carrier["four_simplices_per_internal_triangle"] == 5
    and carrier["four_simplices"] == 20,
)

rho, q = sp.symbols("rho q", positive=True, real=True)
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
c = (1-q-rho) / 2
gram = sp.Matrix([
    [-rho, c, c, c],
    [c, 1, sp.Rational(1, 2), sp.Rational(1, 2)],
    [c, sp.Rational(1, 2), 1, sp.Rational(1, 2)],
    [c, sp.Rational(1, 2), sp.Rational(1, 2), 1],
])
link_block = gram[1:, 1:]
square = (q+rho-1)**2
positive_polynomial = 3*square + 8*rho
expected_determinant = -positive_polynomial/16
gram_determinant = sp.factor(gram.det())
check(
    "the signed Gram determinant is exact and strictly negative",
    sp.simplify(gram_determinant-expected_determinant) == 0,
    f"det G = {expected_determinant}",
)

link_eigenvalues = sorted(
    [sp.simplify(value) for value in link_block.eigenvals().keys()],
    key=float,
)
schur_complement = sp.factor(gram_determinant/link_block.det())
check(
    "the link block plus negative Schur complement prove inertia (3+,1-)",
    link_eigenvalues == [sp.Rational(1, 2), sp.Integer(2)]
    and link_block.eigenvals()[sp.Rational(1, 2)] == 2
    and sp.simplify(schur_complement + positive_polynomial/8) == 0,
    f"link eigenvalues with multiplicity={{1/2: 2, 2: 1}}, Schur={schur_complement}",
)
check(
    "Lorentzian realizability has no hidden upper bound on rho or q",
    sp.expand(positive_polynomial)
    == sp.expand(3*(q+rho-1)**2+8*rho),
    "the determinant is nonzero for every rho>0 and q>0",
)

# Pairwise spacelike new edges do not by themselves make the final cone
# tetrahedron a spacelike hypersurface.  Its Gram eigenvalues expose the
# extra q>1/3 condition.  The no-go below holds on the larger q>0 set and
# therefore also on this physically stricter two-spacelike-boundary sector.
final_boundary_gram = sp.Matrix([
    [q, q-sp.Rational(1, 2), q-sp.Rational(1, 2)],
    [q-sp.Rational(1, 2), q, q-sp.Rational(1, 2)],
    [q-sp.Rational(1, 2), q-sp.Rational(1, 2), q],
])
check(
    "a genuinely spacelike final tetrahedral star additionally requires q>1/3",
    final_boundary_gram.eigenvals()
    == {sp.Rational(1, 2): 2, 3*q-1: 1},
    "pairwise spacelike edges q>0 were not sufficient; the vacuum no-go still covers q>1/3",
)

# The internal triangle is spanned by y and u1.  Its signed two-Gram
# determinant is negative, so use the real magnitude of the timelike area.
triangle_gram = gram.extract([0, 1], [0, 1])
area_radicand = sp.factor(4*rho+(1-q-rho)**2)
area = sp.sqrt(area_radicand)/4
check(
    "each internal triangle is timelike with the exact real area magnitude",
    sp.simplify(triangle_gram.det()+area_radicand/4) == 0,
    f"A_t/a^2 = {area}",
)
area_derivative_rho = sp.factor(sp.diff(area, rho))
expected_area_derivative = (q+rho+1)/(4*sp.sqrt(area_radicand))
check(
    "the timelike hinge area is strictly increasing with proper pole time",
    sp.simplify(area_derivative_rho-expected_area_derivative) == 0,
    "dA/drho=(q+rho+1)/(4 sqrt(F))>0 and drho/dtau>0",
)

gram_inverse = sp.simplify(gram.inv())
D = sp.factor(square+3*rho)
normal_norm = sp.factor(gram_inverse[2, 2])
expected_normal_norm = 4*D/positive_polynomial
check(
    "the facets meeting at a timelike hinge have spacelike normals",
    sp.simplify(normal_norm-expected_normal_norm) == 0,
    f"normal squared = {expected_normal_norm} > 0",
)

inverse_gram_cosine = sp.factor(
    -gram_inverse[2, 3]
    / sp.sqrt(gram_inverse[2, 2]*gram_inverse[3, 3])
)
expected_cosine = sp.factor(
    (q**2+2*q*rho-2*q+rho**2+1)
    /
    (2*(q**2+2*q*rho-2*q+rho**2+rho+1))
)
check(
    "inverse-Gram normals give the preregistered Lorentzian dihedral cosine",
    # SymPy retains Abs(positive polynomial) after sqrt(x**2).  Equality of
    # squares plus N=(q+rho-1)^2+2rho>0 and D>0 fixes the positive branch.
    sp.simplify(inverse_gram_cosine**2-expected_cosine**2) == 0
    and sp.simplify(
        q**2+2*q*rho-2*q+rho**2+1-(square+2*rho)
    ) == 0,
    f"cos(theta) = {expected_cosine}",
)

range_form = sp.Rational(1, 2)-rho/(2*D)
check(
    "the angle formula reduces to 1/2-rho/(2[(q+rho-1)^2+3rho])",
    sp.simplify(expected_cosine-range_form) == 0,
)
one_third_gap = sp.factor(sp.Rational(1, 3)-(sqrt5-1)/4)
check(
    "the full angle range lies strictly below the fivefold target 2*pi/5",
    sp.simplify(one_third_gap-(7-3*sqrt5)/12) == 0
    and bool(sp.N(one_third_gap) > 0),
    "1/3<=cos(theta)<1/2 and 1/3>cos(2*pi/5)",
)

target_cosine = (sqrt5-1)/4
positive_target_gap = sp.factor(
    (3-sqrt5)*(square+rho/phi**2)/(4*D),
    extension=sqrt5,
)
check(
    "the target gap is an exact positive sum of squares",
    sp.simplify(expected_cosine-target_cosine-positive_target_gap) == 0,
    "cos(theta)-cos(72 deg)=(3-sqrt(5))[(q+rho-1)^2+rho/phi^2]/(4D)>0",
)

minimum_deficit = 2*sp.pi-5*sp.acos(sp.Rational(1, 3))
check(
    "the fivefold deficit is uniformly positive",
    float(sp.N(minimum_deficit, 30)) > 0,
    f"epsilon >= {float(sp.N(minimum_deficit, 16)):.15f} rad",
)
check(
    "the zero-volume Lorentzian pole equation has no stationary point",
    float(sp.N(minimum_deficit, 30)) > 0
    and bool(sp.N(expected_area_derivative.subs({rho: 1, q: 1})) > 0),
    "dS/dtau=12 epsilon dA/dtau has one nonzero sign on rho,q>0",
)

# Independent coordinate construction in Minkowski space.  These samples
# include static/nonstatic and small/large proper-time points.
control_points = [
    (0.2, 1.0),
    (0.7, 1.0),
    (4.0, 1.0),
    (0.4, 0.3),
    (3.0, 5.0),
]
controls = [minkowski_control(*point) for point in control_points]
check(
    "explicit Minkowski coordinates reconstruct every prescribed edge interval",
    max(item["gram_residual"] for item in controls) < 2e-14,
    f"max residual={max(item['gram_residual'] for item in controls):.3e}",
)
check(
    "explicit facet normals are spacelike at all deterministic controls",
    max(
        abs(value-1.0)
        for item in controls for value in item["normal_norms"]
    ) < 2e-14
    and all(
        sum(value < -1e-12 for value in item["gram_eigenvalues"]) == 1
        for item in controls
    ),
)
check(
    "Minkowski-normal angles independently match the symbolic formula",
    max(abs(item["angle_residual"]) for item in controls) < 2e-14,
    f"max angle residual={max(abs(item['angle_residual']) for item in controls):.3e}",
)

# Static q=1 specialization and generic volume coefficient.
static_determinant = sp.factor(gram_determinant.subs(q, 1))
static_cosine = sp.factor(expected_cosine.subs(q, 1))
static_area = sp.factor(area.subs(q, 1))
static_volume = sp.sqrt(rho*(8+3*rho))/96
check(
    "the static Lorentzian determinant and angle reduce exactly",
    sp.simplify(static_determinant+rho*(8+3*rho)/16) == 0
    and sp.simplify(static_cosine-(2+rho)/(2*(3+rho))) == 0,
)
check(
    "the static timelike area and four-volume reduce exactly",
    sp.simplify(static_area-sp.sqrt(rho*(4+rho))/4) == 0
    and sp.simplify(96*static_volume-sp.sqrt(rho*(8+3*rho))) == 0,
)

static_theta = sp.acos(static_cosine)
static_epsilon = 2*sp.pi-5*static_theta
# S = 12 A epsilon - lambda * 20 V.  Schlaefli removes derivatives of
# angles, and ell=lambda*a^2 is dimensionless.
required_ell = sp.factor(
    sp.Rational(12, 20)
    * static_epsilon * sp.diff(static_area, rho)
    / sp.diff(static_volume, rho)
)
expected_ell = sp.factor(
    sp.Rational(72, 5)*static_epsilon*(rho+2)/(3*rho+4)
    * sp.sqrt((3*rho+8)/(rho+4))
)
check(
    "a generic volume term gives the exact unselected root coefficient curve",
    sp.simplify(required_ell-expected_ell) == 0,
    "ell=lambda*a^2 is a function of the desired rho, not a selected constant",
)
ell_zero = sp.simplify(sp.limit(required_ell, rho, 0, dir="+"))
ell_infinity = sp.simplify(sp.limit(required_ell, rho, sp.oo))
expected_zero = sp.Rational(36, 5)*sp.sqrt(2)*minimum_deficit
expected_infinity = sp.Rational(8, 5)*sp.sqrt(3)*sp.pi
check(
    "the volume-coefficient curve has the two preregistered exact limits",
    sp.simplify(ell_zero-expected_zero) == 0
    and sp.simplify(ell_infinity-expected_infinity) == 0,
    f"ell(0+)={float(sp.N(ell_zero)):.12f}, ell(infinity)={float(sp.N(ell_infinity)):.12f}",
)
check(
    "the coefficient needed for a root is genuinely rho-dependent",
    abs(float(sp.N(required_ell.subs(rho, 1)))-float(sp.N(ell_zero))) > 1.0,
    f"ell(1)={float(sp.N(required_ell.subs(rho, 1))):.12f}",
)

golden_rho = sp.simplify(phi**-2)
golden_lorentzian_deficit = sp.N(
    static_epsilon.subs(rho, golden_rho), 30
)
check(
    "the Euclidean golden ratio is not a Lorentzian vacuum root",
    float(golden_lorentzian_deficit) > float(sp.N(minimum_deficit)),
    f"epsilon_L(phi^-2)={float(golden_lorentzian_deficit):.12f} rad",
)

verdict = (
    "DERIVED SCOPED LORENTZIAN VACUUM NO-GO"
    if passed == tests else "REFUTED OR INCOMPLETE"
)
result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "external_primary_sources": [
        "https://arxiv.org/abs/1110.5694",
        "https://arxiv.org/abs/1810.09042",
        "https://arxiv.org/abs/1108.1974",
    ],
    "hypotheses": {
        "carrier": "[v,v'] joined with the icosahedral link L_v",
        "old_and_link_squared_length": "a^2",
        "new_cone_squared_length": "q*a^2, q>0 fixed",
        "tent_pole_squared_interval": "-rho*a^2, rho>0 varied",
        "signature": "(-,+,+,+)",
        "action": "ordinary real Lorentzian Regge action on timelike hinges",
        "volume_coefficient": 0,
        "matter_or_higher_terms": False,
    },
    "carrier": {
        "internal_timelike_hinges": 12,
        "four_simplices_per_hinge": 5,
    },
    "exact_geometry": {
        "gram_determinant": "-[3*(q+rho-1)^2+8*rho]/16",
        "inertia": [3, 1],
        "inertia_order": ["positive", "negative"],
        "spacelike_final_boundary_condition": "q>1/3",
        "timelike_hinge_area_over_a2": "sqrt(4*rho+(1-q-rho)^2)/4",
        "dihedral_cosine": (
            "(q^2+2*q*rho-2*q+rho^2+1)/"
            "(2*(q^2+2*q*rho-2*q+rho^2+rho+1))"
        ),
        "dihedral_cosine_range": "1/3 <= cos(theta) < 1/2",
        "angle_range": "pi/3 < theta <= acos(1/3) < 2*pi/5",
        "minimum_deficit_radians": float(sp.N(minimum_deficit, 17)),
        "stationary_points_zero_volume": 0,
    },
    "minkowski_coordinate_controls": controls,
    "static_branch": {
        "gram_determinant": "-rho*(8+3*rho)/16",
        "dihedral_cosine": "(2+rho)/(2*(3+rho))",
        "hinge_area_over_a2": "sqrt(rho*(4+rho))/4",
        "simplex_volume_over_a4": "sqrt(rho*(8+3*rho))/96",
        "golden_rho_deficit_radians": float(golden_lorentzian_deficit),
    },
    "volume_control": {
        "coefficient_definition": "ell=lambda*a^2 for action term -lambda*V_total",
        "required_ell": (
            "(72/5)*epsilon*(rho+2)/(3*rho+4)"
            "*sqrt((3*rho+8)/(rho+4))"
        ),
        "ell_zero_limit": str(ell_zero),
        "ell_zero_limit_float": float(sp.N(ell_zero, 17)),
        "ell_infinity_limit": str(ell_infinity),
        "ell_infinity_limit_float": float(sp.N(ell_infinity, 17)),
        "ell_at_rho_one": float(sp.N(required_ell.subs(rho, 1), 17)),
        "root_is_selected_without_lambda": False,
    },
    "derived": [
        "all q>0,rho>0 simplices in the symmetric family are Lorentzian and nondegenerate",
        "every internal hinge is timelike and has an ordinary Euclidean normal plane",
        "the fivefold deficit is strictly positive throughout the real Lorentzian family",
        "the zero-volume-coefficient symmetric pole equation has no stationary point",
        "releasing the final symmetric spatial radius q does not restore a vacuum root",
    ],
    "structural": [
        "ordinary Lorentzian Regge calculus is adopted as the dynamics",
        "the local symmetric tent ansatz is the physical sector being tested",
    ],
    "open": [
        "a theory-selected volume/cosmological coefficient",
        "nonsymmetric Lorentzian boundary data",
        "matter or higher-curvature terms",
        "a non-overlapping global causal update carrier",
        "c, G, Planck time and Planck mass",
    ],
    "not_claimed": [
        "a no-go for all Lorentzian Regge geometries",
        "a no-go for dynamics outside the local symmetric tent carrier",
        "the generic volume coefficient equals a measured cosmological constant",
    ],
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
with OUTPUT.open("w") as stream:
    json.dump(result, stream, indent=2)
    stream.write("\n")

print("-" * 78)
print(f"VERDICT: {verdict}")
print(f"RESULT: {passed}/{tests} checks passed")
print(f"Machine-readable result: {OUTPUT}")

raise SystemExit(0 if passed == tests else 1)
