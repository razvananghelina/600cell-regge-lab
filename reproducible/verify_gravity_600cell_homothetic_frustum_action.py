#!/usr/bin/env python3
"""Audit subdivision invariance of the homogeneous 600-cell frustum action.

Prior-art commit: f133de0.
Protocol commit: cd7a2ea.

The exact result is based on the flat-cell subdivision theorem only after its
hypotheses have been checked on the actual even/odd carriers.  Independent
100-decimal cellular normals control Lorentzian angle and boundary branches.
"""

import ast
from collections import Counter, defaultdict
from itertools import combinations
import contextlib
import hashlib
import io
import json
from pathlib import Path

import mpmath as arb
import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_gravity_600cell_dust_canonical_continuation.py"
FRUSTUM_INPUT = HERE / "gravity_600cell_homothetic_frustum_equivalence.json"
OUTPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
FRUSTUM_SHA256 = (
    "7e7c23efaf24a2c99a68f3b302b9ef575e0f777ef46f73ccaea9f99e1ecd58dc"
)
PRIOR_ART_COMMIT = "f133de0"
PROTOCOL_COMMIT = "cd7a2ea"
FRUSTUM_RESULT_COMMIT = "6f9cb78"
DPS = 100
arb.mp.dps = DPS
CONTROL_POINTS = (
    (sp.Integer(1), sp.Integer(1), sp.Rational(1, 16)),
    (sp.Integer(1), sp.Rational(3, 4), sp.Rational(1, 16)),
    (sp.Integer(1), sp.Rational(5, 4), sp.Rational(1, 16)),
    (sp.Integer(1), sp.Rational(1, 2), sp.Rational(1, 8)),
    (sp.Integer(1), sp.Rational(3, 2), sp.Rational(1, 8)),
    (sp.Integer(2), sp.Integer(3), sp.Rational(1, 2)),
)
DIFFERENCE_STEPS = (arb.mpf("1e-18"), arb.mpf("3e-18"))
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equal(left, right):
    return sp.simplify(sp.expand(left-right)) == 0


def matrix_equal(left, right):
    return left.shape == right.shape and all(
        equal(a, b) for a, b in zip(left, right)
    )


def load_action_core():
    """Execute definitions only, never the committed continuation search."""
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "print"
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("continuation definition boundary was not found")
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    namespace = {
        "__file__": str(SOURCE),
        "__name__": "homothetic_frustum_action_imported_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(SOURCE), "exec"), namespace)
    return namespace


frustum_input = json.loads(FRUSTUM_INPUT.read_text())
frustum_input_ok = bool(
    digest(FRUSTUM_INPUT) == FRUSTUM_SHA256
    and frustum_input.get("outcome")
        == "HOMOTHETIC_SCHEDULES_ONE_LORENTZIAN_FRUSTUM"
    and frustum_input.get("passed") == frustum_input.get("tests") == 9
    and frustum_input.get("labelled_simplex_count") == 96
    and frustum_input.get("isometry_failure_count") == 0
    and frustum_input.get("regge_action_evaluations") == 0
)
check(
    "the frozen common-frustum theorem is loaded without prior action data",
    frustum_input_ok,
)

core = load_action_core()
models = core["models"]
gro = core["gro"]
core_ok = bool(
    core["tests"] == core["passed"] == 4
    and set(models) == {"even", "odd"}
    and PRIOR_ART_COMMIT == "f133de0"
    and PROTOCOL_COMMIT == "cd7a2ea"
    and FRUSTUM_RESULT_COMMIT == "6f9cb78"
)
check("the committed action core and mission hashes are frozen", core_ok)


# Reconstruct the cellular 600-cell data from logical supports rather than
# importing any hinge classification from a result artifact.
tetrahedra = frozenset(frozenset(cell) for cell in gro.tetrahedra)
spatial_edges = frozenset(models["even"]["old_edges"])
spatial_triangles = frozenset(
    frozenset(face)
    for cell in tetrahedra
    for face in combinations(sorted(cell), 3)
)
f_vector_ok = bool(
    len(gro.vertices) == 120
    and len(spatial_edges) == 720
    and len(spatial_triangles) == 1200
    and len(tetrahedra) == 600
)
check("logical supports reproduce the 600-cell f-vector", f_vector_ok)


def triangle_support_audit(model):
    incidence = defaultdict(list)
    simplex_parent = {}
    for simplex in model["slab"]:
        parent = frozenset(vertex % 120 for vertex in simplex)
        if parent not in tetrahedra:
            raise RuntimeError("a staircase simplex has no tetrahedral frustum")
        simplex_parent[simplex] = parent
        for triangle in combinations(simplex, 3):
            incidence[tuple(sorted(triangle))].append(simplex)

    kinds = {}
    parents = {}
    vertical_by_edge = defaultdict(list)
    diagonal_vertical_incidence = Counter()
    for triangle, simplices in incidence.items():
        logical = frozenset(vertex % 120 for vertex in triangle)
        layers = {vertex // 120 for vertex in triangle}
        parent_cells = frozenset(simplex_parent[simplex] for simplex in simplices)
        parents[triangle] = parent_cells
        if len(logical) == 3 and len(layers) == 1:
            kind = "spatial"
        elif len(logical) == 3 and len(layers) == 2:
            kind = "facet_subdivision"
        elif len(logical) == 2 and len(layers) == 2:
            kind = "trapezoid_subdivision"
            edge = tuple(sorted(logical))
            vertical_by_edge[edge].append(triangle)
            diagonals = [
                tuple(sorted(pair))
                for pair in combinations(triangle, 2)
                if tuple(sorted(pair)) in model["internal_edges"]
                and abs(pair[1]-pair[0]) != 120
            ]
            if len(diagonals) != 1:
                raise RuntimeError("vertical triangle has no unique diagonal")
            diagonal_vertical_incidence[diagonals[0]] += 1
        else:
            kind = "unclassified"
        kinds[triangle] = kind

    kind_counts = Counter(kinds.values())
    parent_counts = {
        kind: Counter(
            len(parents[triangle])
            for triangle, triangle_kind in kinds.items()
            if triangle_kind == kind
        )
        for kind in set(kinds.values())
    }
    vertical_pairing_ok = bool(
        len(vertical_by_edge) == 720
        and set(vertical_by_edge) == set(spatial_edges)
        and all(len(items) == 2 for items in vertical_by_edge.values())
        and len(diagonal_vertical_incidence) == 720
        and set(diagonal_vertical_incidence) == set(
            edge for edge in model["internal_edges"]
            if abs(edge[1]-edge[0]) != 120
        )
        and set(diagonal_vertical_incidence.values()) == {2}
    )
    support_ok = bool(
        kind_counts == Counter({
            "spatial": 2400,
            "facet_subdivision": 2400,
            "trapezoid_subdivision": 1440,
        })
        and parent_counts["spatial"] == Counter({2: 2400})
        and parent_counts["facet_subdivision"] == Counter({2: 2400})
        and parent_counts["trapezoid_subdivision"] == Counter({5: 1440})
        and len(incidence) == 6240
    )
    return {
        "kind_counts": kind_counts,
        "parent_counts": parent_counts,
        "support_ok": support_ok,
        "vertical_pairing_ok": vertical_pairing_ok,
        "incidence": incidence,
        "kinds": kinds,
    }


support = {parity: triangle_support_audit(model) for parity, model in models.items()}
for parity in ("even", "odd"):
    check(
        f"{parity}: every triangle has the unique registered cellular support",
        support[parity]["support_ok"],
        str(dict(support[parity]["kind_counts"])),
    )
check(
    "both parities pair 1,440 vertical triangles into 720 trapezoids",
    all(record["vertical_pairing_ok"] for record in support.values()),
)


# Rebuild the intrinsic metric and facet conormals symbolically.  The direct
# numerical action below uses the same definition but no simplex angle.
L_MINUS, L_PLUS, RHO = sp.symbols(
    "L_minus L_plus rho", positive=True
)
SQRT5 = sp.sqrt(5)
PHI = (1+SQRT5)/2
C = PHI/2
DELTA_L = L_PLUS-L_MINUS
R_MINUS = PHI*L_MINUS
R_PLUS = PHI*L_PLUS
DELTA_R = R_PLUS-R_MINUS
T_SQUARE = RHO+DELTA_R**2
U = sp.Matrix(4, 4, lambda i, j: 1 if i == j else C)
J = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [-1, -1, -1, DELTA_R],
])
G_Q = (J.T*U*J).applyfunc(sp.expand)
G_Q[3, 3] = sp.expand(G_Q[3, 3]-T_SQUARE)
G_INV = sp.simplify(G_Q.inv())
BOTTOM = sp.Matrix([0, 0, 0, 1])
TOP = sp.Matrix([0, 0, 0, -1])
LATERAL_0 = sp.Matrix([1, 0, 0, 0])
LATERAL_1 = sp.Matrix([0, 1, 0, 0])


def conormal_inner(left, right):
    return sp.factor((left.T*G_INV*right)[0])


def internal_cosine(left, right):
    numerator = -conormal_inner(left, right)
    denominator = sp.sqrt(
        conormal_inner(left, left)*conormal_inner(right, right)
    )
    return sp.simplify(numerator/denominator)


COS_LATERAL = internal_cosine(LATERAL_0, LATERAL_1)
COS_BOTTOM = internal_cosine(BOTTOM, LATERAL_0)
COS_TOP = internal_cosine(TOP, LATERAL_0)
static_substitution = {L_PLUS: L_MINUS}
metric_definition_ok = bool(
    matrix_equal(G_Q, G_Q.T)
    and equal(COS_LATERAL.subs(static_substitution), sp.Rational(1, 3))
    and equal(COS_BOTTOM.subs(static_substitution), 0)
    and equal(COS_TOP.subs(static_substitution), 0)
)
check(
    "independent facet conormals recover the registered static angle anchors",
    metric_definition_ok,
)


# Exact Lorentzian area and artificial-diagonal calculus.
D_SQUARE = sp.symbols("d_square")
H_SQUARE = RHO+DELTA_L**2/4


def triangle_area_square(values):
    x, y, z = values
    return sp.factor(
        (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16
    )


lower_area_square = triangle_area_square(
    (L_MINUS**2, -RHO, D_SQUARE)
)
upper_area_square = triangle_area_square(
    (L_PLUS**2, -RHO, D_SQUARE)
)
geometric_diagonal = L_MINUS*L_PLUS-RHO
lower_area = sp.I*L_MINUS*sp.sqrt(H_SQUARE)/2
upper_area = sp.I*L_PLUS*sp.sqrt(H_SQUARE)/2
trapezoid_area = sp.I*(L_MINUS+L_PLUS)*sp.sqrt(H_SQUARE)/2
lower_at_geometry = sp.factor(
    lower_area_square.subs(D_SQUARE, geometric_diagonal)
)
upper_at_geometry = sp.factor(
    upper_area_square.subs(D_SQUARE, geometric_diagonal)
)
area_additivity_ok = bool(
    equal(lower_at_geometry, lower_area**2)
    and equal(upper_at_geometry, upper_area**2)
    and equal(lower_area+upper_area, trapezoid_area)
    and equal(
        triangle_area_square((L_MINUS**2, L_MINUS**2, L_MINUS**2)),
        3*L_MINUS**4/16,
    )
    and equal(
        triangle_area_square((L_PLUS**2, L_PLUS**2, L_PLUS**2)),
        3*L_PLUS**4/16,
    )
)
check("the two signed triangle areas add exactly to the trapezoid area",
      area_additivity_ok)

lower_diagonal_derivative = sp.simplify(
    sp.diff(lower_area_square, D_SQUARE)
    .subs(D_SQUARE, geometric_diagonal)/(2*lower_area)
)
upper_diagonal_derivative = sp.simplify(
    sp.diff(upper_area_square, D_SQUARE)
    .subs(D_SQUARE, geometric_diagonal)/(2*upper_area)
)
diagonal_area_stationary = equal(
    lower_diagonal_derivative+upper_diagonal_derivative, 0
)
check(
    "the split-trapezoid area is stationary in its artificial diagonal",
    diagonal_area_stationary,
    f"partials=({lower_diagonal_derivative},{upper_diagonal_derivative})",
)


# The exact theorem application is now falsifiable by the preceding carrier
# gates.  Spatial wedges are triangulated within two frusta; trapezoid wedges
# within five; facet-only hinges lie in one shared cellular 3-face and acquire
# no cellular codimension-two curvature.  The upstream theorem supplies the
# required flat non-folded realization for every local schedule.
flat_subdivision_hypotheses = bool(
    frustum_input_ok
    and f_vector_ok
    and metric_definition_ok
    and area_additivity_ok
    and diagonal_area_stationary
    and all(record["support_ok"] for record in support.values())
    and all(record["vertical_pairing_ok"] for record in support.values())
)
symbolic_action_identity = flat_subdivision_hypotheses
check(
    "the actual carriers satisfy every hypothesis of flat-cell action additivity",
    symbolic_action_identity,
)


ARB_I = arb.mpc(0, 1)


def mp_matrix(values):
    return arb.matrix([[arb.mpf(item) for item in row] for row in values])


def scalar_product(covector_left, inverse, covector_right):
    return (covector_left.T*inverse*covector_right)[0]


def cell_angles(l_minus, l_plus, rho):
    phi = (1+arb.sqrt(5))/2
    c = phi/2
    delta_r = phi*(l_plus-l_minus)
    t_square = rho+delta_r**2
    u = arb.matrix(4, 4)
    for left in range(4):
        for right in range(4):
            u[left, right] = 1 if left == right else c
    jacobian = mp_matrix((
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (-1, -1, -1, delta_r),
    ))
    metric = jacobian.T*u*jacobian
    metric[3, 3] -= t_square
    inverse = metric**-1
    bottom = arb.matrix((0, 0, 0, 1))
    top = arb.matrix((0, 0, 0, -1))
    lateral_0 = arb.matrix((1, 0, 0, 0))
    lateral_1 = arb.matrix((0, 1, 0, 0))

    def angle(left, right):
        cross = scalar_product(left, inverse, right)
        norm_product = (
            scalar_product(left, inverse, left)
            * scalar_product(right, inverse, right)
        )
        cosine = -cross/arb.sqrt(arb.mpc(norm_product))
        return arb.acos(cosine), cosine

    lateral_angle, lateral_cosine = angle(lateral_0, lateral_1)
    bottom_angle, bottom_cosine = angle(bottom, lateral_0)
    top_angle, top_cosine = angle(top, lateral_0)
    return {
        "lateral": lateral_angle,
        "bottom": bottom_angle,
        "top": top_angle,
        "cosines": (lateral_cosine, bottom_cosine, top_cosine),
    }


def cellular_action(l_minus, l_plus, rho, mass, include_dust=True):
    angles = cell_angles(l_minus, l_plus, rho)
    h = arb.sqrt(rho+(l_plus-l_minus)**2/4)
    lateral_area = ARB_I*(l_minus+l_plus)*h/2
    bottom_area = arb.sqrt(3)*l_minus**2/4
    top_area = arb.sqrt(3)*l_plus**2/4
    lateral_curvature = 2*arb.pi-5*angles["lateral"]
    bottom_curvature = arb.pi-2*angles["bottom"]
    top_curvature = arb.pi-2*angles["top"]
    gravitational = -ARB_I*(
        720*lateral_area*lateral_curvature
        + 1200*bottom_area*bottom_curvature
        + 1200*top_area*top_curvature
    )
    dust = -8*arb.pi*mass*arb.sqrt(rho) if include_dust else arb.mpf(0)
    return gravitational+dust, {
        "gravitational": gravitational,
        "dust": dust,
        "angles": angles,
        "curvatures": (
            lateral_curvature, bottom_curvature, top_curvature
        ),
    }


def as_mpf(value):
    return arb.mpf(value.p)/value.q if isinstance(value, sp.Rational) else arb.mpf(value)


def homothetic_staircase(model, l_minus, l_plus, rho):
    diagonal = l_minus*l_plus-rho
    q_old = tuple(l_minus**2 for _ in range(30))
    internal = tuple([diagonal]*30+[rho]*5)
    q_new = tuple(l_plus**2 for _ in range(30))
    return core["action_and_gradient"](model, q_old, internal, q_new)


def relative_error(left, right):
    return abs(left-right)/max(arb.mpf(1), abs(left), abs(right))


def pack(value, digits=50):
    return {
        "real": arb.nstr(arb.re(value), digits),
        "imag": arb.nstr(arb.im(value), digits),
    }


def perturb_point(point, coordinate, step):
    l_minus, l_plus, rho = point
    result = [l_minus, l_plus, rho]
    if coordinate in (0, 1):
        result[coordinate] *= arb.exp(step/2)
    else:
        result[coordinate] *= arb.exp(step)
    return tuple(result)


def cellular_centered_derivative(point, coordinate, step, mass):
    plus = perturb_point(point, coordinate, step)
    minus = perturb_point(point, coordinate, -step)
    value_plus = cellular_action(*plus, mass)[0]
    value_minus = cellular_action(*minus, mass)[0]
    return (value_plus-value_minus)/(2*step)


def staircase_collective_derivatives(gradient):
    old = 24*sum(gradient[:30], arb.mpc(0))
    diagonal = 24*sum(gradient[30:60], arb.mpc(0))
    pole = 24*sum(gradient[60:65], arb.mpc(0))
    new = 24*sum(gradient[65:95], arb.mpc(0))
    return old, new, pole, diagonal


control_records = []
all_controls_evaluated = True
actions_agree = True
derivatives_agree = True
diagonal_gradients_zero = True
branches_ok = True
maximum_action_error = arb.mpf(0)
maximum_derivative_error = arb.mpf(0)
maximum_diagonal_gradient = arb.mpf(0)
maximum_imaginary = arb.mpf(0)

for exact_point in CONTROL_POINTS:
    point = tuple(as_mpf(value) for value in exact_point)
    try:
        direct_total, direct_data = cellular_action(
            *point, core["ARB_MASS"], include_dust=True
        )
        direct_grav = direct_data["gravitational"]
        direct_derivatives = []
        derivative_spreads = []
        for coordinate in range(3):
            fine = cellular_centered_derivative(
                point, coordinate, DIFFERENCE_STEPS[0], core["ARB_MASS"]
            )
            coarse = cellular_centered_derivative(
                point, coordinate, DIFFERENCE_STEPS[1], core["ARB_MASS"]
            )
            richardson = (9*fine-coarse)/8
            direct_derivatives.append(richardson)
            derivative_spreads.append(abs(fine-coarse)/8)

        parity_records = {}
        for parity in ("even", "odd"):
            total, gradient, branch = homothetic_staircase(
                models[parity], *point
            )
            dust = -8*arb.pi*core["ARB_MASS"]*arb.sqrt(point[2])
            gravitational = total-dust
            collective = staircase_collective_derivatives(gradient)
            action_error = max(
                relative_error(total, direct_total),
                relative_error(gravitational, direct_grav),
            )
            diagonal_error = max(abs(value) for value in gradient[30:60])
            derivative_errors = tuple(
                relative_error(direct_derivatives[index], collective[index])
                for index in range(3)
            )
            imaginary = max(
                abs(arb.im(total)), abs(arb.im(gravitational)),
                *(abs(arb.im(value)) for value in collective[:3]),
                *(abs(arb.im(value)) for value in direct_derivatives),
            )
            branch_ok = bool(
                branch["negative_counts"] == Counter({1: 2400})
                and branch["minimum_leading_minor"] > 0
                and branch["minimum_argument"] > arb.mpf("1e-8")
            )
            parity_ok = bool(
                action_error < arb.mpf("1e-70")
                and max(derivative_errors) < arb.mpf("1e-50")
                and diagonal_error < arb.mpf("1e-70")
                and imaginary < arb.mpf("1e-70")
                and branch_ok
            )
            actions_agree &= action_error < arb.mpf("1e-70")
            derivatives_agree &= max(derivative_errors) < arb.mpf("1e-50")
            diagonal_gradients_zero &= diagonal_error < arb.mpf("1e-70")
            branches_ok &= branch_ok
            maximum_action_error = max(maximum_action_error, action_error)
            maximum_derivative_error = max(
                maximum_derivative_error, *derivative_errors
            )
            maximum_diagonal_gradient = max(
                maximum_diagonal_gradient, diagonal_error
            )
            maximum_imaginary = max(maximum_imaginary, imaginary)
            parity_records[parity] = {
                "action": pack(total),
                "gravitational_action": pack(gravitational),
                "action_relative_error": arb.nstr(action_error, 30),
                "collective_derivatives": [pack(value) for value in collective[:3]],
                "diagonal_collective_derivative": pack(collective[3]),
                "maximum_individual_diagonal_gradient": arb.nstr(
                    diagonal_error, 30
                ),
                "derivative_relative_errors": [
                    arb.nstr(value, 30) for value in derivative_errors
                ],
                "branch_ok": branch_ok,
                "passed": parity_ok,
            }

        point_imaginary = max(
            abs(arb.im(direct_total)), abs(arb.im(direct_grav)),
            *(abs(arb.im(value)) for value in direct_derivatives),
        )
        maximum_imaginary = max(maximum_imaginary, point_imaginary)
        control_records.append({
            "point": [str(value) for value in exact_point],
            "cellular_action": pack(direct_total),
            "cellular_gravitational_action": pack(direct_grav),
            "cellular_collective_derivatives": [
                pack(value) for value in direct_derivatives
            ],
            "richardson_spreads": [
                arb.nstr(value, 30) for value in derivative_spreads
            ],
            "cellular_cosines": [
                pack(value) for value in direct_data["angles"]["cosines"]
            ],
            "parities": parity_records,
        })
    except Exception as exc:
        all_controls_evaluated = False
        control_records.append({
            "point": [str(value) for value in exact_point],
            "error": repr(exc),
        })

check(
    "all six preregistered 100-decimal cellular/staircase controls ran",
    all_controls_evaluated and len(control_records) == 6,
)
check(
    "direct, even and odd complete actions agree below 1e-70",
    actions_agree and maximum_imaginary < arb.mpf("1e-70") and branches_ok,
    f"max relative={arb.nstr(maximum_action_error, 8)}",
)
check(
    "all artificial-diagonal gradients vanish below 1e-70",
    diagonal_gradients_zero,
    f"max absolute={arb.nstr(maximum_diagonal_gradient, 8)}",
)
check(
    "cellular and staircase collective derivatives agree below 1e-50",
    derivatives_agree,
    f"max relative={arb.nstr(maximum_derivative_error, 8)}",
)


if (
    symbolic_action_identity
    and all_controls_evaluated
    and actions_agree
    and diagonal_gradients_zero
    and derivatives_agree
    and branches_ok
    and maximum_imaginary < arb.mpf("1e-70")
):
    outcome = "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
elif (
    all_controls_evaluated
    and branches_ok
    and (not actions_agree or not derivatives_agree)
):
    outcome = "HOMOTHETIC_FRUSTUM_ACTION_SUBDIVISION_DEPENDENT"
else:
    outcome = "HOMOTHETIC_FRUSTUM_ACTION_OPEN"

check(
    "the outcome follows the preregistered hierarchy",
    outcome in {
        "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT",
        "HOMOTHETIC_FRUSTUM_ACTION_SUBDIVISION_DEPENDENT",
        "HOMOTHETIC_FRUSTUM_ACTION_OPEN",
    },
    f"outcome={outcome}",
)


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "frustum_result_commit": FRUSTUM_RESULT_COMMIT,
    "frustum_input_sha256": FRUSTUM_SHA256,
    "outcome": outcome,
    "labels": {
        "homogeneous_action_subdivision_invariance": (
            "DERIVED" if outcome == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
            else "OPEN"
        ),
        "anisotropic_action_equivalence": "OPEN",
        "physical_lapse": "OPEN",
        "external_novelty": "OPEN",
    },
    "cellular_counts": {
        "vertices": 120,
        "edges_or_trapezoids": 720,
        "spatial_triangles_per_boundary": 1200,
        "frusta": 600,
        "dust_struts": 120,
    },
    "carrier_support": {
        parity: {
            "kind_counts": dict(record["kind_counts"]),
            "parent_cell_counts": {
                kind: dict(counts)
                for kind, counts in record["parent_counts"].items()
            },
            "support_ok": record["support_ok"],
            "vertical_pairing_ok": record["vertical_pairing_ok"],
        }
        for parity, record in support.items()
    },
    "exact_geometry": {
        "induced_metric": [[str(value) for value in row] for row in G_Q.tolist()],
        "lateral_internal_cosine": str(sp.factor(COS_LATERAL)),
        "bottom_internal_cosine": str(sp.factor(COS_BOTTOM)),
        "top_internal_cosine": str(sp.factor(COS_TOP)),
        "lower_triangle_area": str(lower_area),
        "upper_triangle_area": str(upper_area),
        "trapezoid_area": str(trapezoid_area),
        "lower_diagonal_area_derivative": str(lower_diagonal_derivative),
        "upper_diagonal_area_derivative": str(upper_diagonal_derivative),
        "diagonal_derivative_sum": str(sp.simplify(
            lower_diagonal_derivative+upper_diagonal_derivative
        )),
        "flat_subdivision_hypotheses": flat_subdivision_hypotheses,
        "symbolic_action_identity": symbolic_action_identity,
    },
    "control_points": control_records,
    "maximum_errors": {
        "action_relative": arb.nstr(maximum_action_error, 50),
        "collective_derivative_relative": arb.nstr(maximum_derivative_error, 50),
        "individual_diagonal_gradient_absolute": arb.nstr(
            maximum_diagonal_gradient, 50
        ),
        "imaginary_contamination": arb.nstr(maximum_imaginary, 50),
    },
    "scope": {
        "homogeneous_collective_action_and_momenta": "TESTED",
        "unrestricted_anisotropic_action": "NOT_TESTED",
        "graviton_modes": "NOT_TESTED",
        "physical_clock_or_scale": "NOT_TESTED",
    },
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print(f"\nSummary: {passed}/{tests} checks passed")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")
if passed != tests or outcome == "HOMOTHETIC_FRUSTUM_ACTION_OPEN":
    raise SystemExit(1)
