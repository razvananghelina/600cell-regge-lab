#!/usr/bin/env python3
"""Exact regular-lapse identity on the 600-cell dust slab.

Prior-art commit: 898d123.
Protocol commit: cf492b9.

The theorem candidate is decided by exact SymPy algebra plus independent
100-digit evaluations of the already certified geometric action core.  A
finite numerical sample cannot promote the outcome to PROVED.
"""

import ast
from collections import Counter, defaultdict
import contextlib
import io
from itertools import combinations, permutations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_gravity_600cell_dust_canonical_continuation.py"
OUTPUT = HERE / "gravity_600cell_dust_regular_lapse_identity.json"
PRIOR_ART_COMMIT = "898d123"
PROTOCOL_COMMIT = "cf492b9"
CONTINUATION_RESULT_COMMIT = "fcfe0c9"
CONTROL_RATIOS = (
    sp.Rational(1),
    sp.Rational(3, 4),
    sp.Rational(1, 2),
    sp.Rational(1, 4),
    sp.Rational(1, 16),
    sp.Rational(1, 256),
)
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


def load_continuation_core():
    """Execute only the definition/control prefix of the committed verifier."""
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
        raise RuntimeError("continuation main marker was not found")
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    namespace = {
        "__file__": str(SOURCE),
        "__name__": "regular_lapse_imported_continuation_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(SOURCE), "exec"), namespace)
    return namespace


core = load_continuation_core()
arb = core["arb"]
models = core["models"]

check(
    "the committed continuation definition prefix retains its four controls",
    core["tests"] == core["passed"] == 4,
)
check(
    "the mission commits and two parity carriers are fixed",
    PRIOR_ART_COMMIT == "898d123"
    and PROTOCOL_COMMIT == "cf492b9"
    and CONTINUATION_RESULT_COMMIT == "fcfe0c9"
    and set(models) == {"even", "odd"},
)


PAIRS5 = tuple(combinations(range(5), 2))
PERMUTATIONS5 = tuple(permutations(range(5)))


def edge_kind(model, edge):
    edge = tuple(sorted(edge))
    if edge in model["old_edges"] or edge in model["new_edges"]:
        return "S"
    if edge[1]-edge[0] == 120:
        return "P"
    if edge in model["internal_edges"]:
        return "D"
    raise ValueError(f"edge outside slab: {edge}")


def raw_angle_key(model, simplex, triangle):
    simplex = tuple(sorted(simplex))
    labels = tuple(
        edge_kind(model, (simplex[left], simplex[right]))
        for left, right in PAIRS5
    )
    marked = tuple(
        index for index, vertex in enumerate(simplex) if vertex in triangle
    )
    return labels, marked


def canonical_angle_string(raw_key):
    labels, marked = raw_key
    matrix = [[None]*5 for _ in range(5)]
    for label, (left, right) in zip(labels, PAIRS5):
        matrix[left][right] = matrix[right][left] = label
    marked = set(marked)
    best = None
    for permutation in PERMUTATIONS5:
        bits = tuple(int(permutation[index] in marked) for index in range(5))
        edges = tuple(
            matrix[permutation[left]][permutation[right]]
            for left, right in PAIRS5
        )
        candidate = bits, edges
        if best is None or candidate < best:
            best = candidate
    if best[0] != (0, 0, 1, 1, 1):
        raise RuntimeError("canonical hinge marking is not 00111")
    return "".join(best[1])


def triangle_and_angle_patterns(model):
    incidence = defaultdict(list)
    for simplex in model["slab"]:
        for triangle in combinations(simplex, 3):
            incidence[tuple(sorted(triangle))].append(simplex)

    raw_keys = {
        raw_angle_key(model, simplex, triangle)
        for triangle, simplices in incidence.items()
        for simplex in simplices
    }
    canonical = {key: canonical_angle_string(key) for key in raw_keys}
    triangle_kinds = Counter()
    patterns = Counter()
    vertical_triangles = []
    for triangle, simplices in incidence.items():
        labels = "".join(sorted(
            edge_kind(model, edge) for edge in combinations(triangle, 2)
        ))
        boundary = triangle in model["boundary_triangles"]
        kind = ("B" if boundary else "I", labels, len(simplices))
        triangle_kinds[kind] += 1
        pattern = tuple(sorted(
            canonical[raw_angle_key(model, simplex, triangle)]
            for simplex in simplices
        ))
        patterns[(("B" if boundary else "I"), labels, pattern)] += 1
        if labels == "DPS":
            vertical_triangles.append(triangle)

    edge_incidence = {"S": Counter(), "D": Counter(), "P": Counter()}
    for triangle in vertical_triangles:
        for edge in combinations(triangle, 2):
            edge = tuple(sorted(edge))
            edge_incidence[edge_kind(model, edge)][edge] += 1
    return incidence, triangle_kinds, patterns, edge_incidence


topology = {}
for parity, model in models.items():
    incidence, kinds, patterns, edge_incidence = triangle_and_angle_patterns(model)
    topology[parity] = {
        "incidence": incidence,
        "kinds": kinds,
        "patterns": patterns,
        "edge_incidence": edge_incidence,
    }

expected_triangle_kinds = Counter({
    ("B", "SSS", 2): 1584,
    ("B", "SSS", 3): 432,
    ("B", "SSS", 4): 384,
    ("I", "DDS", 4): 1584,
    ("I", "DDS", 5): 432,
    ("I", "DDS", 6): 384,
    ("I", "DPS", 5): 1440,
})

topology_ok = all(
    data["kinds"] == expected_triangle_kinds
    and Counter(data["edge_incidence"]["S"].values()) == Counter({1: 1440})
    and Counter(data["edge_incidence"]["D"].values()) == Counter({2: 720})
    and Counter(data["edge_incidence"]["P"].values()) == Counter({12: 120})
    for data in topology.values()
)
parity_patterns_equal = topology["even"]["patterns"] == topology["odd"]["patterns"]
check(
    "both slabs have the exact product-hinge and edge-incidence census",
    topology_ok,
    "2400 boundary SSS, 2400 flat internal DDS, 1440 vertical DPS; "
    "S:D:P incidences are 1:2:12 per edge",
)
check(
    "the complete angle-incidence pattern is identical in both parities",
    parity_patterns_equal and len(topology["even"]["patterns"]) == 15,
)


u = sp.symbols("u", positive=True)
EDGE_VALUE = {"S": sp.Integer(1), "D": 1-u, "P": -u}


def signed_volume_square_symbolic(squared, vertices):
    vertices = tuple(vertices)
    if len(vertices) == 1:
        return sp.Integer(1)
    base = vertices[0]
    others = vertices[1:]
    gram = sp.Matrix([
        [
            (squared[base, left]+squared[base, right]-squared[left, right])/2
            for right in others
        ]
        for left in others
    ])
    return sp.factor(gram.det()/sp.factorial(len(others))**2)


def symbolic_angle_data(edge_string):
    squared = sp.zeros(5)
    for label, (left, right) in zip(edge_string, PAIRS5):
        squared[left, right] = squared[right, left] = EDGE_VALUE[label]
    gram = sp.Matrix([
        [
            (squared[0, left]+squared[0, right]-squared[left, right])/2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram.inv()
    simplex_volume = signed_volume_square_symbolic(squared, range(5))
    facet_a = signed_volume_square_symbolic(squared, (1, 2, 3, 4))
    facet_b = signed_volume_square_symbolic(squared, (0, 2, 3, 4))
    hinge_volume = signed_volume_square_symbolic(squared, (2, 3, 4))
    derivative = sp.zeros(4)
    opposite = {0, 1}
    for left in range(1, 5):
        for right in range(1, 5):
            derivative[left-1, right-1] = (
                int({0, left} == opposite)
                + int({0, right} == opposite)
                - int(left != right and {left, right} == opposite)
            )/2
    volume_derivative = sp.factor(
        simplex_volume*sp.trace(inverse*derivative)
    )
    denominator = sp.sqrt(facet_a)*sp.sqrt(facet_b)
    cosine = sp.factor(sp.simplify(16*volume_derivative/denominator))
    sine = sp.factor(sp.simplify(
        -sp.Rational(4, 3)*sp.sqrt(hinge_volume)
        * sp.sqrt(simplex_volume)/denominator
    ))
    argument = sp.factor(sp.simplify(cosine+sp.I*sine))
    leading_minors = tuple(
        sp.factor(gram[:size, :size].det()) for size in range(1, 5)
    )
    return {
        "simplex_volume": simplex_volume,
        "facet_a": facet_a,
        "facet_b": facet_b,
        "hinge_volume": hinge_volume,
        "cosine": cosine,
        "sine": sine,
        "argument": argument,
        "leading_minors": leading_minors,
    }


angle_strings = sorted({
    string
    for _, _, pattern in topology["even"]["patterns"]
    for string in pattern
})
angle_data = {string: symbolic_angle_data(string) for string in angle_strings}

angle_metric_ok = bool(
    len(angle_strings) == 11
    and all(data["simplex_volume"] == -u/sp.Integer(1152)
            for data in angle_data.values())
    and all(sp.simplify(data["cosine"]**2+data["sine"]**2) == 1
            for data in angle_data.values())
)


def no_root_in_open_half(expression):
    polynomial = sp.Poly(sp.factor(expression), u)
    if polynomial.degree() > 1:
        return False
    for root in sp.solve(polynomial.as_expr(), u):
        if root.is_real is True and bool(0 < root < sp.Rational(1, 2)):
            return False
    return True


raw_simplex_strings = set()
for model in models.values():
    for simplex in model["slab"]:
        simplex = tuple(sorted(simplex))
        raw_simplex_strings.add("".join(
            edge_kind(model, (simplex[left], simplex[right]))
            for left, right in PAIRS5
        ))

raw_simplex_data = {
    string: symbolic_angle_data(string) for string in raw_simplex_strings
}
gram_branch_ok = True
minor_forms = set()
for data in raw_simplex_data.values():
    minors = data["leading_minors"]
    minor_forms.add(tuple(map(str, minors)))
    if not all(no_root_in_open_half(value) for value in minors):
        gram_branch_ok = False
        continue
    signs = [1]
    for value in minors:
        at_quarter = sp.sign(sp.simplify(value.subs(u, sp.Rational(1, 4))))
        signs.append(int(at_quarter))
    gram_branch_ok &= all(sign != 0 for sign in signs)
    gram_branch_ok &= sum(left != right for left, right in zip(signs, signs[1:])) == 1

volume_branch_ok = all(
    no_root_in_open_half(value)
    for data in angle_data.values()
    for value in (data["facet_a"], data["facet_b"], data["hinge_volume"])
)

check(
    "exact Gram and angle algebra yields eleven nondegenerate local types",
    angle_metric_ok and gram_branch_ok and volume_branch_ok
    and len(raw_simplex_strings) == 20 and len(minor_forms) == 10,
    "every 4-simplex has volume square -u/1152 and inertia (3,1) on 0<u<1/2",
)


vertical_strings = {
    "DDDSSSDSDP",
    "SDDSDDSSDP",
    "SDSSDSSDPS",
}
positive_real_strings = {"PDDDSSSSSS", "PDDSSSDSDD"}
negative_imaginary_strings = set(angle_strings)-vertical_strings-positive_real_strings
vertical_argument = sp.Rational(1, 3)-2*sp.sqrt(2)*sp.I/3

vertical_arguments_ok = all(
    sp.simplify(angle_data[string]["argument"]-vertical_argument) == 0
    for string in vertical_strings
)
positive_real_form_ok = all(
    sp.simplify(sp.im(angle_data[string]["argument"])) == 0
    for string in positive_real_strings
)
negative_imaginary_form_ok = all(
    sp.simplify(sp.re(angle_data[string]["argument"])) == 0
    for string in negative_imaginary_strings
)

# Exact sign certificates.  Every radical below is real and nonzero on
# 0<u<1/2.  The two nontrivial differences are positive because their
# squared differences are 9-18u and 18-27u respectively.
branch_sign_polynomials = (
    u,
    1-u,
    1-2*u,
    2-3*u,
    3-4*u,
    9-18*u,
    18-27*u,
)
branch_sign_ok = all(
    sp.simplify(poly.subs(u, sp.Rational(1, 4))) > 0
    and no_root_in_open_half(poly)
    for poly in branch_sign_polynomials
)

angle_branch_ok = bool(
    vertical_arguments_ok
    and positive_real_form_ok
    and negative_imaginary_form_ok
    and branch_sign_ok
)
check(
    "all exact angle arguments remain off the logarithm cut on 0<u<1/2",
    angle_branch_ok,
    "each is positive real, negative imaginary, or the fixed fourth-quadrant vertical angle",
)


pattern_products_ok = True
branch_anchors_ok = True
alpha = sp.acos(sp.Rational(1, 3))
pattern_records = []
for (boundary_code, labels, pattern), multiplicity in sorted(
    topology["even"]["patterns"].items(), key=str
):
    product = sp.factor(sp.radsimp(sp.simplify(sp.prod(
        angle_data[string]["argument"] for string in pattern
    ))))
    if boundary_code == "B" and labels == "SSS":
        expected_product = -1
        expected_anchor = -sp.pi
    elif boundary_code == "I" and labels == "DDS":
        expected_product = 1
        expected_anchor = -2*sp.pi
    elif boundary_code == "I" and labels == "DPS":
        expected_product = vertical_argument**5
        expected_anchor = -5*alpha
    else:
        expected_product = sp.nan
        expected_anchor = sp.nan
        pattern_products_ok = False

    pattern_products_ok &= sp.simplify(product-expected_product) == 0
    anchor = sp.Integer(0)
    for string in pattern:
        if string in vertical_strings:
            anchor -= alpha
        elif string in positive_real_strings:
            anchor += 0
        else:
            anchor -= sp.pi/2
    branch_anchors_ok &= sp.simplify(anchor-expected_anchor) == 0
    pattern_records.append({
        "boundary": boundary_code == "B",
        "labels": labels,
        "multiplicity": multiplicity,
        "angle_types": dict(Counter(pattern)),
        "product": str(product),
        "anchor_sum": str(anchor),
    })

curvature_identity_ok = bool(pattern_products_ok and branch_anchors_ok)
check(
    "exact products and anchored branches give only epsilon3 curvature",
    curvature_identity_ok,
    "boundary SSS: 0; internal DDS: 0; vertical DPS: 2*pi-5*acos(1/3)",
)


L, tau, epsilon = sp.symbols("L tau epsilon", positive=True)
mass = sp.Rational(90)*epsilon*L/sp.pi
spatial_square = L**2
diagonal_square = L**2-tau**2
pole_square = -tau**2


def triangle_area_square_symbolic(values):
    x, y, z = values
    return sp.factor((2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16)


def triangle_partials_symbolic(values):
    x, y, z = values
    return (
        sp.factor((y+z-x)/8),
        sp.factor((x+z-y)/8),
        sp.factor((x+y-z)/8),
    )


vertical_values = (spatial_square, diagonal_square, pole_square)
vertical_area_square = triangle_area_square_symbolic(vertical_values)
vertical_area = sp.I*L*tau/2
partials = triangle_partials_symbolic(vertical_values)
log_derivatives = (
    sp.simplify(partials[0]*spatial_square/(2*vertical_area)),
    sp.simplify(partials[1]*diagonal_square/(2*vertical_area)),
    sp.simplify(partials[2]*pole_square/(2*vertical_area)),
)
expected_log_derivatives = (sp.I*L*tau/4, 0, sp.I*L*tau/4)
area_derivative_ok = bool(
    vertical_area_square == -L**2*tau**2/4
    and all(sp.simplify(left-right) == 0
            for left, right in zip(log_derivatives, expected_log_derivatives))
)

gravitational_action = sp.simplify(
    -sp.I*sp.Integer(1440)*vertical_area*epsilon
)
dust_action = sp.simplify(-8*sp.pi*mass*tau)
total_action = sp.simplify(gravitational_action+dust_action)
boundary_gradient = sp.simplify(-sp.I*epsilon*log_derivatives[0])
diagonal_gradient = sp.simplify(
    -sp.I*sp.Integer(2)*epsilon*log_derivatives[1]
)
pole_gravitational_gradient = sp.simplify(
    -sp.I*sp.Integer(12)*epsilon*log_derivatives[2]
)
pole_dust_gradient = sp.simplify(
    -(4*sp.pi*mass/sp.Integer(5))*tau/sp.Integer(24)
)
pole_total_gradient = sp.simplify(
    pole_gravitational_gradient+pole_dust_gradient
)
pre_momentum = sp.simplify(-boundary_gradient)
post_momentum = sp.simplify(boundary_gradient)

formula_flags = {
    "vertical_area_and_log_derivatives": area_derivative_ok,
    "gravitational_action": sp.simplify(
        gravitational_action-720*epsilon*L*tau
    ) == 0,
    "dust_action": sp.simplify(dust_action+720*epsilon*L*tau) == 0,
    "total_action": total_action == 0,
    "diagonal_gradient": diagonal_gradient == 0,
    "pole_gradient": pole_total_gradient == 0,
    "pre_momentum": sp.simplify(
        pre_momentum+epsilon*L*tau/4
    ) == 0,
    "post_momentum": sp.simplify(
        post_momentum-epsilon*L*tau/4
    ) == 0,
}
exact_formula_ok = bool(all(formula_flags.values()))
check(
    "the exact hinge census gives the registered action and all gradients",
    exact_formula_ok,
    "S_grav=720 epsilon3 L tau, S_dust=-S_grav, p_pre=-epsilon3 L tau/4",
)


# The committed base lies inside the wider symbolic domain.  A handwritten
# exact bound used in the result note is rho0<1, zeta>1/2 and R0>4 (from
# 3<pi<10/3), hence L>2 and u0<1/4<1/2.  This numerical line is only a
# consistency control, not the proof of that elementary inequality.
base_u = core["ARB_RHO"]/core["ARB_L0_SQUARE"]
check(
    "the committed positive-lapse interval lies inside 0<u<1/2",
    0 < base_u < arb.mpf("0.5"),
    f"u0={arb.nstr(base_u, 30)}",
)


numerical_records = []
numerical_formula_ok = True
numerical_evaluation_complete = True
for parity in ("even", "odd"):
    model = models[parity]
    for ratio_exact in CONTROL_RATIOS:
        ratio = arb.mpf(ratio_exact.p)/ratio_exact.q
        rho = ratio*core["ARB_RHO"]
        slant = core["ARB_L0_SQUARE"]-rho
        x = tuple([slant]*30+[rho]*5)
        q_new = tuple([core["ARB_L0_SQUARE"]]*30)
        try:
            total, gradient, branch = core["action_and_gradient"](
                model, core["ARB_BASE_OLD"], x, q_new
            )
        except Exception as exc:  # recorded as numerical openness
            numerical_evaluation_complete = False
            numerical_records.append({
                "parity": parity,
                "rho_ratio": str(ratio_exact),
                "error": repr(exc),
            })
            continue

        tau_value = arb.sqrt(rho)
        epsilon_value = core["ARB_EPSILON_3"]
        L_value = core["ARB_L0"]
        expected_pre = -epsilon_value*L_value*tau_value/4
        expected_post = -expected_pre
        dust = -(8*arb.pi*core["ARB_MASS"])*tau_value
        gravitational = total-dust
        expected_gravitational = 720*epsilon_value*L_value*tau_value
        internal_error = max(abs(value) for value in gradient[30:65])
        pre_error = max(
            abs(-gradient[index]-expected_pre)/abs(expected_pre)
            for index in range(30)
        )
        post_error = max(
            abs(gradient[65+index]-expected_post)/abs(expected_post)
            for index in range(30)
        )
        action_error = max(
            abs(gravitational-expected_gravitational),
            abs(dust+expected_gravitational),
            abs(total),
        )/max(arb.mpf(1), abs(gravitational), abs(dust))
        maximum_imaginary = max(
            abs(arb.im(total)), *(abs(arb.im(value)) for value in gradient)
        )
        branch_ok = core["branch_pass"](branch, maximum_imaginary)
        record_ok = bool(
            branch_ok
            and maximum_imaginary < arb.mpf("1e-70")
            and internal_error < arb.mpf("1e-60")
            and pre_error < arb.mpf("1e-60")
            and post_error < arb.mpf("1e-60")
            and action_error < arb.mpf("1e-60")
        )
        numerical_formula_ok &= record_ok
        numerical_records.append({
            "parity": parity,
            "rho_ratio": str(ratio_exact),
            "branch_ok": bool(branch_ok),
            "internal_gradient_error": arb.nstr(internal_error, 40),
            "pre_relative_error": arb.nstr(pre_error, 40),
            "post_relative_error": arb.nstr(post_error, 40),
            "action_relative_error": arb.nstr(action_error, 40),
            "maximum_imaginary": arb.nstr(maximum_imaginary, 40),
            "passed": record_ok,
        })

check(
    "all twelve independent 100-digit geometric controls were evaluated",
    numerical_evaluation_complete and len(numerical_records) == 12,
)
check(
    "the original geometric action agrees with every exact formula control",
    numerical_formula_ok,
    "six preregistered rho/rho0 values in each parity",
)


stored_maps_ok = all(
    sorted(core["gluing_input"]["parities"][parity]["geometry"][
        "old_to_final_orbit_map"
    ]) == list(range(30))
    for parity in ("even", "odd")
)
time_reflection_ok = bool(
    exact_formula_ok
    and pre_momentum == -post_momentum
    and stored_maps_ok
)
check(
    "uniform exact pre/post momenta obey the derived orbit-map reflection",
    time_reflection_ok,
)


symbolic_complete = bool(
    topology_ok
    and parity_patterns_equal
    and angle_metric_ok
    and gram_branch_ok
    and volume_branch_ok
    and angle_branch_ok
    and curvature_identity_ok
)
if symbolic_complete and exact_formula_ok and numerical_formula_ok:
    outcome = "REGULAR_LAPSE_IDENTITY_PROVED"
elif symbolic_complete and not exact_formula_ok:
    outcome = "REGULAR_LAPSE_IDENTITY_REFUTED"
elif numerical_evaluation_complete and not numerical_formula_ok:
    outcome = "REGULAR_LAPSE_IDENTITY_REFUTED"
elif numerical_formula_ok and not symbolic_complete:
    outcome = "REGULAR_LAPSE_PATTERN_ONLY"
else:
    outcome = "REGULAR_LAPSE_IDENTITY_NUMERICALLY_OPEN"

zero_lapse_consequence = bool(
    outcome == "REGULAR_LAPSE_IDENTITY_PROVED" and time_reflection_ok
)
check(
    "the scientific outcome follows the preregistered hierarchy",
    outcome in {
        "REGULAR_LAPSE_IDENTITY_PROVED",
        "REGULAR_LAPSE_IDENTITY_REFUTED",
        "REGULAR_LAPSE_PATTERN_ONLY",
        "REGULAR_LAPSE_IDENTITY_NUMERICALLY_OPEN",
    },
    f"outcome={outcome}",
)


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "continuation_result_commit": CONTINUATION_RESULT_COMMIT,
    "labels": {
        "identity": "DERIVED" if outcome == "REGULAR_LAPSE_IDENTITY_PROVED" else "OPEN",
        "external_novelty": "OPEN",
        "physical_clock": "OPEN",
        "expansion": "NOT_TESTED",
    },
    "domain": {
        "rho": "0 < rho <= rho0",
        "normalized_symbolic_domain": "0 < u=rho/L^2 < 1/2",
        "base_u": arb.nstr(base_u, 50),
    },
    "topology": {
        "triangle_kinds": {
            f"{boundary}:{labels}:incidence={incidence}": count
            for (boundary, labels, incidence), count
            in sorted(expected_triangle_kinds.items())
        },
        "vertical_edge_incidence": {
            "boundary_spatial": "1440 edges x 1 triangle",
            "diagonal": "720 edges x 2 triangles",
            "pole": "120 edges x 12 triangles",
        },
        "angle_patterns": pattern_records,
        "parities_identical": parity_patterns_equal,
    },
    "exact_algebra": {
        "simplex_volume_square": "-u/1152",
        "local_angle_types": len(angle_strings),
        "raw_simplex_types": len(raw_simplex_strings),
        "leading_minor_forms": len(minor_forms),
        "vertical_angle_argument": str(vertical_argument),
        "vertical_deficit": "2*pi-5*acos(1/3)",
        "formula_flags": formula_flags,
        "gravitational_action": str(gravitational_action),
        "dust_action": str(dust_action),
        "total_action": str(total_action),
        "diagonal_gradient": str(diagonal_gradient),
        "pole_gravitational_gradient": str(pole_gravitational_gradient),
        "pole_dust_gradient": str(pole_dust_gradient),
        "pole_total_gradient": str(pole_total_gradient),
        "pre_momentum_per_edge": str(pre_momentum),
        "post_momentum_per_edge": str(post_momentum),
    },
    "numerical_controls": numerical_records,
    "homotopy_consequence": {
        "time_reflection_exact": time_reflection_ok,
        "target": "p(lambda)=(1-2*lambda)*p_pre(rho0)",
        "positive_lapse_solution": "rho=rho0*(1-2*lambda)^2 for 0<=lambda<1/2",
        "zero_lapse_boundary": "lambda=1/2" if zero_lapse_consequence else None,
        "lambda_is_physical_time": False,
    },
    "outcome": outcome,
    "zero_lapse_consequence": zero_lapse_consequence,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
if zero_lapse_consequence:
    print(
        "DERIVED: the same-orientation regular connected branch reaches the "
        "zero-lapse boundary at lambda=1/2 and supplies no nondegenerate "
        "forward spatial frame."
    )
else:
    print(
        "OPEN/NEGATIVE: the exact zero-lapse consequence is not established."
    )

if passed != tests:
    raise SystemExit(1)

