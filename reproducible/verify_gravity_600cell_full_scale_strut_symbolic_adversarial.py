#!/usr/bin/env python3
"""Generic symbolic two-cell adversarial proof of the scale--strut response."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_full_scale_strut_symbolic_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_symbolic_adversarial_protocol.md"
DISCLOSURE = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_exploratory_disclosure.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_prior_art.md"
FIRST_RESULT = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_carrier_first_result.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_carrier.py"
PRIMARY_JSON = HERE / "gravity_600cell_full_scale_strut_carrier.json"
PRECISION_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_precision.py"
PRECISION_JSON = HERE / "gravity_600cell_full_scale_strut_precision.json"
TWO_RESULT = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_result.md"

PROTOCOL_COMMIT = "2b2ee13"
EXPECTED_HASHES = {
    "protocol": "ab1b82702f4acb95998f35eda91b2472870a0b730e46f6c391e73320cc16fea1",
    "disclosure": "e3dba59118e35cc2370beec4b081cc18fdfd1753fda74a6ca2b0a013d86bd473",
    "prior_art": "3fc6c3e75ad92c3c20bb420d97e26e73fdd62b69a9aac44162fa95c74c29219a",
    "first_result": "5753375ca2a6c4f5152f134474176501b580a1c55b7a871b3a39fa6321d82f61",
    "primary_source": "e68105df4058f7d2ed39a6913f29e88cd9fe88e123ff52260acf698a2bd7da49",
    "primary_json": "6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d",
    "precision_source": "7d8f64470251cb79be267be9c4937c09caf2827f0d4da1e6479d1555db1cd80c",
    "precision_json": "2a2a79271a92fc2ddde343a9d0651402df6eeb4a90efa2697e26f54cafcdf60f",
    "two_result": "b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3",
}
INPUTS = {
    "protocol": PROTOCOL,
    "disclosure": DISCLOSURE,
    "prior_art": PRIOR_ART,
    "first_result": FIRST_RESULT,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "precision_source": PRECISION_SOURCE,
    "precision_json": PRECISION_JSON,
    "two_result": TWO_RESULT,
}

lam, tau = sp.symbols("lambda tau")
A, B, C, D = sp.symbols("A B C D")
ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
CANONICAL = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
SHARED = (0, 1, 2)
FIELD = sp.QQ.frac_field(lam, tau)
ALLOWED_FACTORS = (
    tau,
    lam - 1,
    (lam - 1) ** 2 - 3 * tau**2,
)

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


def square_covector(moving, fixed, metric=ETA):
    return sp.simplify(2 * (moving - fixed).T * metric)


def endpoint_cross_row(source, target):
    result = sp.zeros(1, 8)
    result[0, source] = A
    result[0, target] = B
    result[0, 4 + source] = C
    result[0, 4 + target] = D
    return result


def vertex_solved_response(points, lapse_expression=tau, metric=ETA):
    """Solve each upper vertex from its three crosses and one strut."""
    top = tuple(lam * point + lapse_expression * NORMAL for point in points)
    response = sp.zeros(16, 8)
    determinants = []
    for target in range(4):
        rows = []
        rhs_rows = []
        for source in range(4):
            if source == target:
                continue
            rows.append(square_covector(top[target], points[source], metric))
            rhs_rows.append(endpoint_cross_row(source, target))
        rows.append(square_covector(top[target], points[target], metric))
        strut = sp.zeros(1, 8)
        strut[0, 4 + target] = 1
        rhs_rows.append(strut)
        matrix = sp.Matrix.vstack(*rows)
        rhs = sp.Matrix.vstack(*rhs_rows)
        determinants.append(sp.factor(matrix.det()))
        solved = matrix.inv(method="DM") * rhs
        response[4 * target:4 * target + 4, :] = solved
    return top, sp.simplify(response), determinants


def upper_edge_residuals(points, top, response, metric=ETA):
    residuals = []
    for left, right in combinations(range(4), 2):
        covector = square_covector(top[left], top[right], metric)
        observed = covector * (
            response[4 * left:4 * left + 4, :]
            - response[4 * right:4 * right + 4, :]
        )
        expected = sp.zeros(1, 8)
        expected[0, left] = 8 * lam
        expected[0, right] = 8 * lam
        residuals.extend(sp.together(value) for value in observed - expected)
    return residuals


def nonzero_numerators(expressions):
    result = []
    for expression in expressions:
        numerator, _ = sp.fraction(sp.cancel(expression))
        numerator = sp.factor(numerator)
        if numerator:
            result.append(numerator)
    return result


def groebner_record(expressions, variables):
    numerators = nonzero_numerators(expressions)
    basis = sp.groebner(numerators, *variables, order="lex", domain=FIELD)
    monic = [sp.cancel(polynomial.as_expr()) for polynomial in basis.polys]
    return basis, monic, numerators


def basis_strings(basis):
    return [str(sp.factor(value)) for value in basis]


def normalized_factor(expression):
    polynomial = sp.Poly(sp.expand(expression), lam, tau, domain=sp.QQ)
    return sp.factor(polynomial.monic().as_expr())


def factor_set(expressions):
    factors = set()
    for expression in expressions:
        numerator, denominator = sp.fraction(sp.cancel(expression))
        for part in (numerator, denominator):
            if part == 0:
                continue
            _, values = sp.factor_list(part, gens=(lam, tau))
            for factor, _ in values:
                if sp.Poly(factor, lam, tau).total_degree() > 0:
                    factors.add(normalized_factor(factor))
    return factors


def matrix_expressions(matrix):
    return [matrix[row, column] for row in range(matrix.rows) for column in range(matrix.cols)]


def lorentz_basis(metric=ETA):
    result = []
    for left, right in combinations(range(4), 2):
        generator = sp.zeros(4)
        generator[left, right] = 1
        generator[right, left] = -metric[left, left] / metric[right, right]
        result.append(generator)
    return tuple(result)


def poincare_evaluation(points, metric=ETA):
    columns = [
        sp.Matrix.vstack(*(generator * point for point in points))
        for generator in lorentz_basis(metric)
    ]
    for axis in range(4):
        vector = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(vector for _ in points)))
    return sp.Matrix.hstack(*columns)


def affine_reflection(point, hyperplane_points, metric=ETA):
    origin = hyperplane_points[0]
    tangents = sp.Matrix.hstack(*(value - origin for value in hyperplane_points[1:]))
    normal_vectors = (tangents.T * metric).nullspace()
    if len(normal_vectors) != 1:
        raise RuntimeError("lateral hyperplane normal is not unique")
    normal = normal_vectors[0]
    norm = sp.factor((normal.T * metric * normal)[0])
    coefficient = sp.cancel(
        2 * ((point - origin).T * metric * normal)[0] / norm
    )
    return sp.simplify(point - coefficient * normal), normal, norm


def affine_map(domain, codomain):
    domain_h = sp.Matrix.hstack(*(homogeneous(point) for point in domain))
    codomain_h = sp.Matrix.hstack(*(homogeneous(point) for point in codomain))
    determinant = sp.factor(domain_h.det())
    return sp.simplify(codomain_h * domain_h.inv(method="DM")), determinant


def lateral_transition(points, lapse_expression=tau, metric=ETA):
    top = tuple(lam * point + lapse_expression * NORMAL for point in points)
    lateral = (points[0], points[1], points[2], top[0])
    reflected_lower, normal, normal_norm = affine_reflection(points[3], lateral, metric)
    reflected_upper, _, _ = affine_reflection(top[3], lateral, metric)
    domain = (points[0], points[1], points[2], top[0], points[3])
    codomain = (points[0], points[1], points[2], top[0], reflected_lower)
    transition, domain_determinant = affine_map(domain, codomain)
    linear = transition[:4, :4]
    mapping_residuals = []
    for index in SHARED:
        mapping_residuals.extend(transition * homogeneous(points[index]) - homogeneous(points[index]))
        mapping_residuals.extend(transition * homogeneous(top[index]) - homogeneous(top[index]))
    mapping_residuals.extend(transition * homogeneous(top[3]) - homogeneous(reflected_upper))
    control = bool(
        all(sp.cancel(value) == 0 for value in mapping_residuals)
        and sp.simplify(linear.T * metric * linear) == metric
        and sp.factor(linear.det() ** 2 - 1) == 0
        and transition[4, :] == sp.Matrix([[0, 0, 0, 0, 1]])
    )
    return {
        "transition": transition,
        "linear": linear,
        "top": top,
        "normal": normal,
        "normal_norm": normal_norm,
        "domain_determinant": domain_determinant,
        "control": control,
    }


def shared_selector():
    result = sp.zeros(12, 16)
    for position, vertex in enumerate(SHARED):
        result[4 * position:4 * position + 4, 4 * vertex:4 * vertex + 4] = sp.eye(4)
    return result


def data_embedding(vertices):
    result = sp.zeros(8, 10)
    for local, global_vertex in enumerate(vertices):
        result[local, global_vertex] = 1
        result[4 + local, 5 + global_vertex] = 1
    return result


def explicit_annihilator(column):
    pivot = next((row for row in range(column.rows) if sp.cancel(column[row]) != 0), None)
    if pivot is None:
        return None, None
    rows = []
    for row in range(column.rows):
        if row == pivot:
            continue
        vector = sp.zeros(1, column.rows)
        vector[0, row] = 1
        vector[0, pivot] = -sp.cancel(column[row] / column[pivot])
        rows.append(vector)
    return sp.Matrix.vstack(*rows), pivot


def gluing_case(points=CANONICAL, lapse_expression=tau):
    top, response, determinants = vertex_solved_response(points, lapse_expression)
    local_residuals = upper_edge_residuals(points, top, response)
    transition = lateral_transition(points, lapse_expression)
    selector = shared_selector()
    source_embedding = data_embedding((0, 1, 2, 3))
    target_embedding = data_embedding((0, 1, 2, 4))
    linear_block = sp.diag(*([transition["linear"]] * 3))
    difference = sp.simplify(
        selector * response * source_embedding
        - linear_block * selector * response * target_embedding
    )
    lower_evaluation = poincare_evaluation(tuple(points[index] for index in SHARED))
    lower_kernel_vectors = lower_evaluation.nullspace()
    lower_kernel = sp.Matrix.hstack(*lower_kernel_vectors)
    upper_evaluation = poincare_evaluation(tuple(top[index] for index in SHARED))
    connection = sp.simplify(upper_evaluation * lower_kernel)
    annihilator, pivot = explicit_annihilator(connection)
    annihilator_control = bool(
        lower_evaluation.rank() == 9
        and lower_kernel.shape == (10, 1)
        and connection.shape == (12, 1)
        and connection.rank() == 1
        and annihilator is not None
        and annihilator.shape == (11, 12)
        and annihilator.rank() == 11
        and sp.simplify(annihilator * connection) == sp.zeros(11, 1)
    )
    substitutions = {A: 8 - B, C: 1 - D}
    correct = sp.simplify((annihilator * difference).subs(substitutions))
    fixed = sp.simplify(difference.subs(substitutions))
    exception_expressions = (
        determinants
        + [transition["normal_norm"], transition["domain_determinant"], connection[pivot]]
        + matrix_expressions(response)
        + matrix_expressions(transition["transition"])
        + matrix_expressions(annihilator)
        + matrix_expressions(correct)
    )
    return {
        "top": top,
        "response": response,
        "determinants": determinants,
        "local_residuals": local_residuals,
        "transition": transition,
        "lower_evaluation_rank": lower_evaluation.rank(),
        "lower_kernel_dimension": lower_kernel.cols,
        "connection_rank": connection.rank(),
        "connection": connection,
        "annihilator": annihilator,
        "annihilator_pivot": pivot,
        "annihilator_control": annihilator_control,
        "correct_residual": correct,
        "fixed_residual": fixed,
        "exception_expressions": exception_expressions,
    }


def reduced_gluing_basis(matrix):
    expressions = matrix_expressions(matrix)
    basis, monic, numerators = groebner_record(expressions, (B, D))
    return basis, monic, numerators


print("=" * 78)
print("SYMBOLIC TWO-CELL ADVERSARIAL SCALE--STRUT PROOF")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all symbolic-adversarial inputs retain frozen provenance", provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
precision = json.loads(PRECISION_JSON.read_text())
upstream_ok = bool(
    primary.get("outcome") == "FULL_SCALE_STRUT_NUMERICALLY_OPEN"
    and primary.get("passed") == primary.get("tests") == 18
    and primary.get("finite_formula_agrees") is True
    and all(
        record["candidate_mismatch_count"] == 0
        and record["full_residual_nonzero_count"] == 0
        for record in primary["new_exact_global_controls"]
    )
    and all(primary["parities"][parity]["complete_exact_rank"] == 240 for parity in ("even", "odd"))
    and precision.get("outcome") == "FULL_SCALE_STRUT_PRECISION_RESOLVED"
    and precision.get("passed") == precision.get("tests") == 9
)
check("finite global controls and resolved precision remain frozen", upstream_ok)

regular_geometry_ok = bool(
    {sp.simplify(((CANONICAL[i] - CANONICAL[j]).T * ETA * (CANONICAL[i] - CANONICAL[j]))[0]) for i, j in combinations(range(4), 2)} == {8}
    and {sp.simplify((point.T * ETA * point)[0]) for point in CANONICAL} == {3}
)
check("the symbolic carrier is exactly the normalized regular tetrahedron", regular_geometry_ok)

print("building baseline generic two-cell system", flush=True)
baseline = gluing_case()
local_basis_object, local_basis, local_numerators = groebner_record(
    baseline["local_residuals"], (A, B, C, D)
)
expected_local_object = sp.groebner(
    [A + B - 8, C + D - 1], A, B, C, D,
    order="lex", domain=FIELD,
)
expected_local = [sp.cancel(value.as_expr()) for value in expected_local_object.polys]
local_ideal_ok = local_basis == expected_local
check("one generic cell has exactly the two disclosed compatibility generators", local_ideal_ok, str(basis_strings(local_basis)))

transition_ok = bool(
    baseline["transition"]["control"]
    and baseline["annihilator_control"]
    and baseline["lower_evaluation_rank"] == 9
    and baseline["lower_kernel_dimension"] == 1
    and baseline["connection_rank"] == 1
)
check("the exact lateral Lorentz transition retains one nonzero face-connection column", transition_ok, f"pivot={baseline['annihilator_pivot']}")

correct_object, correct_basis, correct_numerators = reduced_gluing_basis(baseline["correct_residual"])
expected_generators = [
    (lam - 1) ** 2 * B - 2 * (lam - 1) ** 2 - 2 * tau**2,
    (lam - 1) * D - lam,
]
expected_object = sp.groebner(
    expected_generators, B, D, order="lex", domain=FIELD
)
expected_basis = [sp.cancel(value.as_expr()) for value in expected_object.polys]
formula_agrees = correct_basis == expected_basis
check("the correct generic gluing ideal is evaluated and classified", bool(correct_basis and expected_basis), f"agrees={formula_agrees}, basis={basis_strings(correct_basis)}")

observed_factors = factor_set(baseline["exception_expressions"] + expected_generators)
allowed_normalized = {normalized_factor(value) for value in ALLOWED_FACTORS}
unexpected_factors = observed_factors - allowed_normalized
hypothesis_gap = bool(unexpected_factors)
check("all generic denominator and rank-exception factors are classified", True, f"observed={sorted(map(str, observed_factors))}, unexpected={sorted(map(str, unexpected_factors))}")

fixed_object, fixed_basis, fixed_numerators = reduced_gluing_basis(baseline["fixed_residual"])
fixed_frame_rejected = fixed_basis != expected_basis
check("omitting the face-connection column does not counterfeit the accepted ideal", fixed_frame_rejected, str(basis_strings(fixed_basis)))

B_expected = 2 + 2 * tau**2 / (lam - 1) ** 2
D_expected = lam / (lam - 1)
corrupted_values = [
    sp.cancel(value.subs({B: B_expected, D: D_expected + 1}))
    for value in matrix_expressions(baseline["correct_residual"])
]
corrupted_nonzero = [value for value in corrupted_values if value != 0]
corruption_ok = bool(corrupted_nonzero)
check("the exact D+1 corruption leaves a nonzero shared-face residual", corruption_ok, f"nonzero={len(corrupted_nonzero)}")

print("building time-reversed generic control", flush=True)
time_reversed = gluing_case(lapse_expression=-tau)
_, reversed_basis_raw, reversed_numerators = reduced_gluing_basis(time_reversed["correct_residual"])
reversed_basis = [sp.cancel(value.subs(tau, -tau)) for value in reversed_basis_raw]
time_reversal_ok = bool(
    time_reversed["transition"]["control"]
    and time_reversed["annihilator_control"]
    and reversed_basis == expected_basis
)
check("time reversal preserves the generic scalar coefficient ideal", time_reversal_ok, str(basis_strings(reversed_basis)))

print("building odd-relabelled generic control", flush=True)
odd_points = (CANONICAL[1], CANONICAL[0], CANONICAL[2], CANONICAL[3])
odd = gluing_case(points=odd_points)
_, odd_basis, odd_numerators = reduced_gluing_basis(odd["correct_residual"])
odd_ok = bool(
    odd["transition"]["control"]
    and odd["annihilator_control"]
    and odd_basis == expected_basis
)
check("an odd consistent local relabelling preserves all four scalar coefficients", odd_ok, str(basis_strings(odd_basis)))

finite_records = []
finite_ok = True
for record in primary["new_exact_global_controls"]:
    scale, lapse = map(sp.Integer, record["representative"])
    values = (
        6 - 2 * lapse**2 / (scale - 1) ** 2,
        2 + 2 * lapse**2 / (scale - 1) ** 2,
        -1 / (scale - 1),
        scale / (scale - 1),
    )
    frozen = tuple(sp.Rational(value) for value in record["candidate_coefficients"])
    match = values == frozen
    finite_ok &= match
    finite_records.append({
        "representative": record["representative"],
        "symbolic_substitution": [str(value) for value in values],
        "frozen": [str(value) for value in frozen],
        "match": match,
    })
check("post-symbolic substitution reproduces all three frozen finite controls", finite_ok, str(finite_records))

firewall_ok = bool(
    set(hashes) == set(INPUTS)
    and not any("hessian" in str(path).lower() or "action" in str(path).lower() or "strong" in str(path).lower() for path in INPUTS.values())
)
check("the symbolic firewall loaded no global builder, action, Hessian or strong-equation source", firewall_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and regular_geometry_ok and local_ideal_ok
    and transition_ok and fixed_frame_rejected and corruption_ok
    and time_reversal_ok and odd_ok and finite_ok and firewall_ok
)
if not controls_ok:
    outcome = "FULL_SCALE_STRUT_SYMBOLIC_CONTROL_FAILED"
elif hypothesis_gap:
    outcome = "FULL_SCALE_STRUT_SYMBOLIC_HYPOTHESIS_GAP"
elif not formula_agrees:
    outcome = "FULL_SCALE_STRUT_SYMBOLIC_DISAGREEMENT"
elif formula_agrees:
    outcome = "FULL_SCALE_STRUT_SYMBOLIC_CORROBORATED"
else:
    outcome = "FULL_SCALE_STRUT_SYMBOLIC_OPEN"

allowed = {
    "FULL_SCALE_STRUT_SYMBOLIC_CONTROL_FAILED",
    "FULL_SCALE_STRUT_SYMBOLIC_HYPOTHESIS_GAP",
    "FULL_SCALE_STRUT_SYMBOLIC_DISAGREEMENT",
    "FULL_SCALE_STRUT_SYMBOLIC_CORROBORATED",
    "FULL_SCALE_STRUT_SYMBOLIC_OPEN",
}
check("the preregistered symbolic hierarchy assigns one outcome", outcome in allowed, outcome)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "generic_hypotheses": [
        "lambda != 1",
        "tau != 0",
        "(lambda-1)^2 - 3*tau^2 != 0",
    ],
    "vertex_solve_determinants": [str(value) for value in baseline["determinants"]],
    "local_nonzero_residual_numerators": [str(value) for value in local_numerators],
    "local_groebner_basis": basis_strings(local_basis),
    "transition": {
        "linear_determinant": str(sp.factor(baseline["transition"]["linear"].det())),
        "normal_norm": str(sp.factor(baseline["transition"]["normal_norm"])),
        "domain_determinant": str(sp.factor(baseline["transition"]["domain_determinant"])),
        "lower_poincare_evaluation_rank": baseline["lower_evaluation_rank"],
        "lower_face_stabilizer_dimension": baseline["lower_kernel_dimension"],
        "connection_rank": baseline["connection_rank"],
        "annihilator_pivot": baseline["annihilator_pivot"],
    },
    "correct_gluing_nonzero_residual_numerators": [str(value) for value in correct_numerators],
    "correct_gluing_groebner_basis": basis_strings(correct_basis),
    "expected_gluing_groebner_basis": basis_strings(expected_basis),
    "formula_agrees": formula_agrees,
    "observed_exception_factors": sorted(map(str, observed_factors)),
    "allowed_exception_factors": sorted(map(str, allowed_normalized)),
    "unexpected_exception_factors": sorted(map(str, unexpected_factors)),
    "fixed_frame_groebner_basis": basis_strings(fixed_basis),
    "fixed_frame_equals_correct": fixed_basis == expected_basis,
    "corrupted_D_plus_one_nonzero_residuals": [str(value) for value in corrupted_nonzero],
    "time_reversed_nonzero_residual_numerators": [str(value) for value in reversed_numerators],
    "time_reversed_groebner_basis": basis_strings(reversed_basis),
    "odd_relabelled_nonzero_residual_numerators": [str(value) for value in odd_numerators],
    "odd_relabelled_groebner_basis": basis_strings(odd_basis),
    "finite_regression_controls": finite_records,
    "classification": {
        "generic_endpoint_formula": (
            "DERIVED EXACT, ADVERSARIALLY CORROBORATED"
            if outcome == "FULL_SCALE_STRUT_SYMBOLIC_CORROBORATED" else "OPEN"
        ),
        "accepted_curved_240_carrier": (
            "DERIVED KINEMATIC"
            if outcome == "FULL_SCALE_STRUT_SYMBOLIC_CORROBORATED" else "OPEN"
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
if passed != tests or outcome == "FULL_SCALE_STRUT_SYMBOLIC_CONTROL_FAILED":
    raise SystemExit(1)
