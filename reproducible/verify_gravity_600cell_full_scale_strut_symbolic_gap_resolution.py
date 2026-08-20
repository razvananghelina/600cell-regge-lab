#!/usr/bin/env python3
"""Coordinate-free wedge resolution of the symbolic scale--strut gap."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_symbolic_gap_resolution_protocol.md"
GAP_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_symbolic_gap_result.md"
FIRST_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_symbolic_adversarial_protocol.md"
FIRST_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_symbolic_adversarial.py"
GAP_JSON = HERE / "gravity_600cell_full_scale_strut_symbolic_adversarial.json"
PRECISION_JSON = HERE / "gravity_600cell_full_scale_strut_precision.json"

PROTOCOL_COMMIT = "d6c62a7"
EXPECTED_HASHES = {
    "protocol": "48b1eab544807e7b8724e312144783f5d4e161e84793dc22bb05b382487f9088",
    "gap_result": "6e81e90427996557e2d9eab62dd1b01d098121d971d66c14f03899d8ed5a017d",
    "first_protocol": "ab1b82702f4acb95998f35eda91b2472870a0b730e46f6c391e73320cc16fea1",
    "first_source": "abf83f50fb5cda2ce6a4820b2d105ce3895da7f91504e886c76b559ad6964f2e",
    "gap_json": "96a0d91528c7cf69e9286a50bed1f81184ef58adc3a6dd4b8e76b62c6c5c223a",
    "precision_json": "2a2a79271a92fc2ddde343a9d0651402df6eeb4a90efa2697e26f54cafcdf60f",
}
INPUTS = {
    "protocol": PROTOCOL,
    "gap_result": GAP_RESULT,
    "first_protocol": FIRST_PROTOCOL,
    "first_source": FIRST_SOURCE,
    "gap_json": GAP_JSON,
    "precision_json": PRECISION_JSON,
}

lam, tau = sp.symbols("lambda tau", real=True)
A, B, C, D = sp.symbols("A B C D")
ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
POINTS = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
SHARED = (0, 1, 2)
GENERIC_FIELD = sp.QQ.frac_field(lam, tau)

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def homogeneous(point):
    return point.col_join(sp.ones(1, 1))


def covector(moving, fixed):
    return sp.simplify(2 * (moving - fixed).T * ETA)


def endpoint_row(source, target):
    result = sp.zeros(1, 8)
    result[0, source] = A
    result[0, target] = B
    result[0, 4 + source] = C
    result[0, 4 + target] = D
    return result


def vertex_response(scale, lapse):
    top = tuple(scale * point + lapse * NORMAL for point in POINTS)
    response = sp.zeros(16, 8)
    determinants = []
    matrices = []
    for target in range(4):
        rows = []
        rhs = []
        for source in range(4):
            if source == target:
                continue
            rows.append(covector(top[target], POINTS[source]))
            rhs.append(endpoint_row(source, target))
        rows.append(covector(top[target], POINTS[target]))
        strut = sp.zeros(1, 8)
        strut[0, 4 + target] = 1
        rhs.append(strut)
        matrix = sp.Matrix.vstack(*rows)
        matrices.append(matrix)
        determinants.append(sp.factor(matrix.det()))
        response[4 * target:4 * target + 4, :] = matrix.inv(method="DM") * sp.Matrix.vstack(*rhs)
    return top, sp.simplify(response), determinants, matrices


def affine_reflection(point, hyperplane):
    origin = hyperplane[0]
    tangents = sp.Matrix.hstack(*(value - origin for value in hyperplane[1:]))
    normals = (tangents.T * ETA).nullspace()
    if len(normals) != 1:
        raise RuntimeError("lateral normal is not one-dimensional")
    normal = normals[0]
    norm = sp.factor((normal.T * ETA * normal)[0])
    coefficient = sp.cancel(2 * ((point - origin).T * ETA * normal)[0] / norm)
    return sp.simplify(point - coefficient * normal), normal, norm


def lateral_transition(scale, lapse):
    top = tuple(scale * point + lapse * NORMAL for point in POINTS)
    hyperplane = (POINTS[0], POINTS[1], POINTS[2], top[0])
    reflected_lower, normal, normal_norm = affine_reflection(POINTS[3], hyperplane)
    reflected_upper, _, _ = affine_reflection(top[3], hyperplane)
    domain = (POINTS[0], POINTS[1], POINTS[2], top[0], POINTS[3])
    codomain = (POINTS[0], POINTS[1], POINTS[2], top[0], reflected_lower)
    domain_h = sp.Matrix.hstack(*(homogeneous(point) for point in domain))
    codomain_h = sp.Matrix.hstack(*(homogeneous(point) for point in codomain))
    domain_det = sp.factor(domain_h.det())
    transition = sp.simplify(codomain_h * domain_h.inv(method="DM"))
    linear = transition[:4, :4]
    mapping = []
    for index in SHARED:
        mapping.extend(transition * homogeneous(POINTS[index]) - homogeneous(POINTS[index]))
        mapping.extend(transition * homogeneous(top[index]) - homogeneous(top[index]))
    mapping.extend(transition * homogeneous(top[3]) - homogeneous(reflected_upper))
    control = bool(
        all(sp.cancel(value) == 0 for value in mapping)
        and sp.simplify(linear.T * ETA * linear) == ETA
        and sp.factor(linear.det() ** 2 - 1) == 0
        and transition[4, :] == sp.Matrix([[0, 0, 0, 0, 1]])
    )
    return {
        "top": top,
        "transition": transition,
        "linear": linear,
        "domain_determinant": domain_det,
        "normal": normal,
        "normal_norm": normal_norm,
        "control": control,
    }


def lorentz_basis():
    result = []
    for left, right in combinations(range(4), 2):
        generator = sp.zeros(4)
        generator[left, right] = 1
        generator[right, left] = -ETA[left, left] / ETA[right, right]
        result.append(generator)
    return tuple(result)


def poincare_evaluation(points):
    columns = [
        sp.Matrix.vstack(*(generator * point for point in points))
        for generator in lorentz_basis()
    ]
    for axis in range(4):
        vector = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(vector for _ in points)))
    return sp.Matrix.hstack(*columns)


def selector():
    result = sp.zeros(12, 16)
    for position, vertex in enumerate(SHARED):
        result[4 * position:4 * position + 4, 4 * vertex:4 * vertex + 4] = sp.eye(4)
    return result


def embedding(vertices):
    result = sp.zeros(8, 10)
    for local, global_vertex in enumerate(vertices):
        result[local, global_vertex] = 1
        result[4 + local, 5 + global_vertex] = 1
    return result


def wedge_equations(connection, residual):
    equations = []
    for column in range(residual.cols):
        for left, right in combinations(range(connection.rows), 2):
            equations.append(sp.cancel(
                connection[left] * residual[right, column]
                - connection[right] * residual[left, column]
            ))
    return equations


def matrix_entries(matrix):
    return [matrix[row, column] for row in range(matrix.rows) for column in range(matrix.cols)]


def build_case(scale, lapse):
    top, response, determinants, vertex_matrices = vertex_response(scale, lapse)
    transition = lateral_transition(scale, lapse)
    shared_selector = selector()
    source = embedding((0, 1, 2, 3))
    target = embedding((0, 1, 2, 4))
    block = sp.diag(*([transition["linear"]] * 3))
    difference = sp.simplify(
        shared_selector * response * source
        - block * shared_selector * response * target
    )
    lower = poincare_evaluation(tuple(POINTS[index] for index in SHARED))
    lower_kernel = sp.Matrix.hstack(*lower.nullspace())
    upper = poincare_evaluation(tuple(top[index] for index in SHARED))
    connection = sp.simplify(upper * lower_kernel)
    substituted = sp.simplify(difference.subs({A: 8 - B, C: 1 - D}))
    wedge = wedge_equations(connection, substituted)
    fixed = matrix_entries(substituted)
    connection_norm = sp.factor(sum(value**2 for value in connection))
    return {
        "top": top,
        "response": response,
        "determinants": determinants,
        "vertex_matrices": vertex_matrices,
        "transition": transition,
        "difference": difference,
        "substituted": substituted,
        "lower_rank": lower.rank(),
        "lower_kernel_dimension": lower_kernel.cols,
        "connection": connection,
        "connection_rank": connection.rank(),
        "connection_norm": connection_norm,
        "wedge": wedge,
        "fixed": fixed,
    }


def numerators(expressions):
    result = []
    for expression in expressions:
        numerator, _ = sp.fraction(sp.cancel(expression))
        numerator = sp.factor(numerator)
        if numerator:
            result.append(numerator)
    return result


def groebner_basis(expressions, generic=True):
    values = numerators(expressions)
    domain = GENERIC_FIELD if generic else sp.QQ
    basis = sp.groebner(values, B, D, order="lex", domain=domain)
    monic = [sp.cancel(polynomial.as_expr()) for polynomial in basis.polys]
    return monic, values


def normalize_factor(expression):
    return sp.factor(sp.Poly(expression, lam, tau, domain=sp.QQ).monic().as_expr())


def factors_from_part(part):
    factors = set()
    if part == 0:
        return factors
    _, values = sp.factor_list(part, gens=(A, B, C, D, lam, tau))
    for factor, _ in values:
        if factor.free_symbols and factor.free_symbols <= {lam, tau}:
            factors.add(normalize_factor(factor))
    return factors


def denominator_factors(expressions):
    result = set()
    for expression in expressions:
        _, denominator = sp.fraction(sp.cancel(expression))
        result |= factors_from_part(denominator)
    return result


def numerator_factors(expressions):
    result = set()
    for expression in expressions:
        numerator, _ = sp.fraction(sp.cancel(expression))
        result |= factors_from_part(numerator)
    return result


def basis_strings(values):
    return [str(sp.factor(value)) for value in values]


def special_record(scale, lapse):
    case = build_case(sp.Integer(scale), sp.Integer(lapse))
    basis, residual_numerators = groebner_basis(case["wedge"], generic=False)
    expected = sp.groebner([
        B - (sp.Rational(2) + sp.Rational(2) * lapse**2 / (scale - 1) ** 2),
        D - sp.Rational(scale, scale - 1),
    ], B, D, order="lex", domain=sp.QQ)
    expected_basis = [sp.cancel(value.as_expr()) for value in expected.polys]
    ranks = [matrix.rank() for matrix in case["vertex_matrices"]]
    ok = bool(
        ranks == [4, 4, 4, 4]
        and case["transition"]["control"]
        and case["connection_rank"] == 1
        and basis == expected_basis
    )
    return ok, {
        "lambda": scale,
        "tau": lapse,
        "vertex_solve_ranks": ranks,
        "transition_control": case["transition"]["control"],
        "connection_rank": case["connection_rank"],
        "wedge_groebner_basis": basis_strings(basis),
        "expected_groebner_basis": basis_strings(expected_basis),
        "nonzero_wedge_numerator_count": len(residual_numerators),
    }


print("=" * 78)
print("COORDINATE-FREE SYMBOLIC GAP RESOLUTION")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all gap-resolution inputs retain frozen provenance", provenance_ok, str(hashes))

gap = json.loads(GAP_JSON.read_text())
precision = json.loads(PRECISION_JSON.read_text())
gap_frozen = bool(
    gap.get("outcome") == "FULL_SCALE_STRUT_SYMBOLIC_HYPOTHESIS_GAP"
    and gap.get("passed") == gap.get("tests") == 14
    and gap.get("formula_agrees") is True
    and set(gap.get("unexpected_exception_factors", [])) == {
        "lambda + tau - 1",
        "lambda - tau - 1",
        "lambda**2 - 2*lambda + 3*tau**2 + 1",
    }
    and precision.get("outcome") == "FULL_SCALE_STRUT_PRECISION_RESOLVED"
    and precision.get("passed") == precision.get("tests") == 9
)
check("the 14/14 hypothesis gap and resolved precision are preserved literally", gap_frozen)

print("building generic pivot-free wedge system", flush=True)
generic = build_case(lam, tau)
geometry_ok = bool(
    generic["transition"]["control"]
    and generic["lower_rank"] == 9
    and generic["lower_kernel_dimension"] == 1
    and generic["connection_rank"] == 1
    and all(value != 0 for value in generic["determinants"])
)
check("the generic construction has four invertible vertex solves and one connection line", geometry_ok)

wedge_basis, wedge_numerators = groebner_basis(generic["wedge"], generic=True)
expected_generators = [
    (lam - 1) ** 2 * B - 2 * (lam - 1) ** 2 - 2 * tau**2,
    (lam - 1) * D - lam,
]
expected_object = sp.groebner(
    expected_generators, B, D, order="lex", domain=GENERIC_FIELD
)
expected_basis = [sp.cancel(value.as_expr()) for value in expected_object.polys]
wedge_agrees = wedge_basis == expected_basis
check("the pivot-free wedge ideal is evaluated and classified", bool(wedge_basis), f"agrees={wedge_agrees}, basis={basis_strings(wedge_basis)}")

fixed_basis, fixed_numerators = groebner_basis(generic["fixed"], generic=True)
fixed_rejected = fixed_basis != expected_basis
check("pivot-free analysis still rejects fixed-frame gluing", fixed_rejected, str(basis_strings(fixed_basis)))

B_expected = 2 + 2 * tau**2 / (lam - 1) ** 2
D_expected = lam / (lam - 1)
corrupted = [sp.cancel(value.subs({B: B_expected, D: D_expected + 1})) for value in generic["wedge"]]
corrupted_nonzero = [value for value in corrupted if value != 0]
corruption_ok = bool(corrupted_nonzero)
check("D+1 leaves nonzero coordinate-free wedge residuals", corruption_ok, f"nonzero={len(corrupted_nonzero)}")

all_rational_expressions = (
    matrix_entries(generic["response"])
    + matrix_entries(generic["transition"]["transition"])
    + matrix_entries(generic["connection"])
    + generic["wedge"]
    + expected_basis
)
actual_denominators = denominator_factors(all_rational_expressions)
rank_expressions = (
    generic["determinants"]
    + [
        generic["transition"]["domain_determinant"],
        generic["transition"]["normal_norm"],
        generic["connection_norm"],
    ]
)
rank_factors = numerator_factors(rank_expressions) | denominator_factors(rank_expressions)
coefficient_numerator_factors = numerator_factors(all_rational_expressions)

tau_factor = normalize_factor(tau)
lambda_factor = normalize_factor(lam - 1)
lateral_factor = normalize_factor((lam - 1) ** 2 - 3 * tau**2)
positive_factor = normalize_factor((lam - 1) ** 2 + 3 * tau**2)
allowed_actual = {tau_factor, lambda_factor, lateral_factor, positive_factor}
genuine_extras = (actual_denominators | rank_factors) - allowed_actual
factors_classified = True
check(
    "actual denominators, rank factors and coefficient numerators are separated",
    factors_classified,
    f"denominators={sorted(map(str, actual_denominators))}, rank={sorted(map(str, rank_factors))}, extras={sorted(map(str, genuine_extras))}",
)

positive_identity = sp.expand(
    ((lam - 1) ** 2 + 3 * tau**2)
    - ((lam - 1) ** 2 + (sp.sqrt(3) * tau) ** 2)
) == 0
connection_norm_numerator, connection_norm_denominator = sp.fraction(
    sp.cancel(generic["connection_norm"])
)
positive_ok = bool(
    positive_identity
    and positive_factor in numerator_factors([connection_norm_numerator])
    and tau_factor in numerator_factors([connection_norm_numerator])
)
check("the extra quadratic is a positive real connection-norm certificate", positive_ok, f"norm={sp.factor(generic['connection_norm'])}")

special_records = []
special_ok = True
for scale, lapse in ((-1, 2), (3, 2)):
    ok, record = special_record(scale, lapse)
    special_ok &= ok
    special_records.append(record)
check("both lambda±tau-1 real strata retain full ranks and the disclosed ideal", special_ok, str(special_records))

firewall_ok = bool(
    set(hashes) == set(INPUTS)
    and not any("action" in str(path).lower() or "hessian" in str(path).lower() or "strong" in str(path).lower() for path in INPUTS.values())
)
check("the gap-resolution firewall loaded no action, Hessian or strong-equation source", firewall_ok)

controls_ok = bool(
    provenance_ok and gap_frozen and geometry_ok and fixed_rejected
    and corruption_ok and factors_classified and positive_ok and special_ok
    and firewall_ok
)
genuine_gap = bool(genuine_extras or not special_ok)
if not controls_ok:
    outcome = "FULL_SCALE_STRUT_GAP_CONTROL_FAILED"
elif genuine_gap:
    outcome = "FULL_SCALE_STRUT_GAP_GENUINE"
elif not wedge_agrees:
    outcome = "FULL_SCALE_STRUT_GAP_DISAGREEMENT"
elif wedge_agrees:
    outcome = "FULL_SCALE_STRUT_GAP_REAL_RESOLVED"
else:
    outcome = "FULL_SCALE_STRUT_GAP_OPEN"

allowed = {
    "FULL_SCALE_STRUT_GAP_CONTROL_FAILED",
    "FULL_SCALE_STRUT_GAP_GENUINE",
    "FULL_SCALE_STRUT_GAP_REAL_RESOLVED",
    "FULL_SCALE_STRUT_GAP_DISAGREEMENT",
    "FULL_SCALE_STRUT_GAP_OPEN",
}
check("the preregistered gap-resolution hierarchy assigns one outcome", outcome in allowed, outcome)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "real_domain": [
        "lambda != 1",
        "tau != 0",
        "(lambda-1)^2-3*tau^2 != 0",
    ],
    "wedge_nonzero_numerator_count": len(wedge_numerators),
    "wedge_groebner_basis": basis_strings(wedge_basis),
    "expected_groebner_basis": basis_strings(expected_basis),
    "wedge_agrees": wedge_agrees,
    "fixed_frame_groebner_basis": basis_strings(fixed_basis),
    "fixed_frame_nonzero_numerator_count": len(fixed_numerators),
    "corrupted_D_plus_one_nonzero_wedge_count": len(corrupted_nonzero),
    "actual_denominator_factors": sorted(map(str, actual_denominators)),
    "rank_factors": sorted(map(str, rank_factors)),
    "coefficient_numerator_factors": sorted(map(str, coefficient_numerator_factors)),
    "allowed_actual_or_positive_factors": sorted(map(str, allowed_actual)),
    "genuine_extra_factors": sorted(map(str, genuine_extras)),
    "connection_component_norm": str(sp.factor(generic["connection_norm"])),
    "positive_connection_norm_certificate": positive_ok,
    "special_linear_strata": special_records,
    "complex_specializations": "OPEN; uniform theorem claimed only on the stated real domain",
    "classification": {
        "generic_real_endpoint_formula": (
            "DERIVED EXACT; SYMBOLIC GAP RESOLVED"
            if outcome == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED" else "OPEN"
        ),
        "accepted_curved_240_carrier": (
            "DERIVED KINEMATIC"
            if outcome == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED" else "OPEN"
        ),
        "action_or_physical_dynamics": "NOT EVALUATED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
if passed != tests or outcome == "FULL_SCALE_STRUT_GAP_CONTROL_FAILED":
    raise SystemExit(1)
