#!/usr/bin/env python3
"""Independent irregular-frustum audit of relative-Poincare stratification."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_cellular_frustum_relative_poincare_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_adversarial_protocol.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_prior_art.md"
CORRECTION_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_covariance_correction_protocol.md"
CORRECTION_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_relative_poincare_covariance_correction.py"
CORRECTION_JSON = HERE / "gravity_600cell_cellular_frustum_relative_poincare_covariance_correction.json"
IRREGULAR_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_anisotropic_rigidity_adversarial.py"

PROTOCOL_COMMIT = "6250871"
EXPECTED_HASHES = {
    "protocol": "84444067b369471c52ffe36ce63e1aef863bfc52dd4408595547067c37cb0803",
    "prior_art": "a8811a441fecd137b37085e4018fea7abb3f365750dfc426d16f1b46e5282e7c",
    "correction_protocol": "ad569b4c4ecbfb4b6d3db7d0225dd26f0426e7e973ddb9c82c21daff5b404b3c",
    "correction_source": "e85e2df690234e19e0343183499c6f8465bc149bc66a87c5246dfe0bda4c1d61",
    "correction_json": "f571869be3341b74b2341c2bf776e99b21174f9f0fb0c5d02e42585c2f3ebaa2",
    "irregular_source": "ecc5e0cb5f8913325f00137245f33299c7607b395d219161bf7e0e806068c18a",
}

ETA = sp.diag(1, 1, 1, -1)
BOTTOM = tuple(sp.Matrix(point) for point in (
    (5, 0, 0, 0),
    (0, 5, 0, 0),
    (0, 0, 5, 0),
    (3, 4, 0, 0),
))
REPRESENTATIVES = ((1, 7), (2, 7), (3, 13))
NORMAL = sp.Matrix((0, 0, 0, 1))
PAIRS = tuple(combinations(range(4), 2))

A_SYMBOLS = sp.symbols("a0:16")
B_SYMBOLS = sp.symbols("b0:4")
PARAMETERS = A_SYMBOLS + B_SYMBOLS
A_MATRIX = sp.Matrix(4, 4, A_SYMBOLS)
B_VECTOR = sp.Matrix(B_SYMBOLS)
Y_SYMBOLS = sp.symbols("y0:16")
Y_POINTS = tuple(sp.Matrix(Y_SYMBOLS[4 * i:4 * i + 4]) for i in range(4))

tests = 0
passed = 0


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
    return sha256(path.read_bytes()).hexdigest()


def nullspace_matrix(matrix):
    vectors = matrix.nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.cols, 0)


def equal_span(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


def square(left, right, metric):
    delta = left - right
    return sp.expand((delta.T * metric * delta)[0])


def top(scale, lapse, bottom=BOTTOM):
    return tuple(sp.Matrix((
        scale * point[0], scale * point[1], scale * point[2], lapse
    )) for point in bottom)


def polynomial_jacobian(bottom, points, metric):
    polynomials = [
        square(Y_POINTS[left], Y_POINTS[right], metric)
        for left, right in PAIRS
    ]
    polynomials.extend(
        square(Y_POINTS[index], bottom[index], metric)
        for index in range(4)
    )
    symbolic = sp.Matrix(polynomials).jacobian(Y_SYMBOLS)
    values = {
        Y_SYMBOLS[4 * index + axis]: points[index][axis]
        for index in range(4) for axis in range(4)
    }
    return symbolic.subs(values)


def lorentz_constraint(metric):
    symmetric = A_MATRIX.T * metric + metric * A_MATRIX
    equations = [symmetric[row, column]
                 for row in range(4) for column in range(row, 4)]
    return sp.Matrix(equations).jacobian(A_SYMBOLS)


def affine_displacement(points):
    expressions = sp.Matrix.vstack(*(A_MATRIX * q + B_VECTOR for q in points))
    return expressions.jacobian(PARAMETERS)


def stabilizer(metric, normal):
    lorentz = lorentz_constraint(metric)
    fixes_normal = (A_MATRIX * normal).jacobian(A_SYMBOLS)
    return nullspace_matrix(lorentz.col_join(fixes_normal))


def analyse(bottom, points, metric):
    direct = polynomial_jacobian(bottom, points, metric)
    lorentz_a = lorentz_constraint(metric)
    lorentz_full = lorentz_a.row_join(sp.zeros(10, 4))
    displacement = affine_displacement(points)
    strut = direct[6:10, :]
    translation = strut * displacement[:, 16:20]
    system = lorentz_full.col_join(strut * displacement)
    parameter_kernel = nullspace_matrix(system)
    displacement_image = displacement * parameter_kernel
    direct_kernel = nullspace_matrix(direct)
    a_image = parameter_kernel[:16, :]
    all_lorentz = nullspace_matrix(lorentz_a)
    return {
        "direct": direct,
        "direct_rank": direct.rank(),
        "direct_kernel": direct_kernel,
        "lorentz_rank": lorentz_a.rank(),
        "all_lorentz": all_lorentz,
        "rigid_displacement_rank": (
            displacement * nullspace_matrix(lorentz_full)
        ).rank(),
        "system_rank": system.rank(),
        "parameter_kernel": parameter_kernel,
        "displacement_image": displacement_image,
        "kernel_equality": equal_span(displacement_image, direct_kernel),
        "a_image": a_image,
        "a_image_rank": a_image.rank(),
        "translation": translation,
        "translation_rank": translation.rank(),
        "pure_translation_dimension": 4 - translation.rank(),
    }


paths = {
    "protocol": PROTOCOL,
    "prior_art": PRIOR_ART,
    "correction_protocol": CORRECTION_PROTOCOL,
    "correction_source": CORRECTION_SOURCE,
    "correction_json": CORRECTION_JSON,
    "irregular_source": IRREGULAR_SOURCE,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all adversarial inputs have exact frozen provenance",
      provenance_ok, str(hashes))

correction = json.loads(CORRECTION_JSON.read_text())
correction_ok = bool(
    correction["outcome"]
    == "STATIC_STABILIZER_AND_EXPANDING_LORENTZ_CHART_CORROBORATED"
    and correction["passed"] == correction["tests"] == 13
    and correction["preserved_original_outcome"]
    == "RELATIVE_POINCARE_CONTROL_FAILED"
)
check("the corrected result and preserved original failure both persist",
      correction_ok)

affine = sp.Matrix.hstack(*(BOTTOM[index] - BOTTOM[0]
                             for index in range(1, 4)))
norms = tuple(square(point, sp.zeros(4, 1), ETA) for point in BOTTOM)
geometry_ok = bool(affine.rank() == 3 and set(norms) == {25})
check("the adversarial tetrahedron is irregular, equinorm and nondegenerate",
      geometry_ok)

lorentz_control = lorentz_constraint(ETA)
minus_lorentz_control = lorentz_constraint(-ETA)
algebra_ok = bool(
    lorentz_control.shape == (10, 16)
    and lorentz_control.rank() == 10
    and nullspace_matrix(lorentz_control).shape == (16, 6)
    and minus_lorentz_control == -lorentz_control
)
check("the redundant 16-entry matrix constraints leave exact so(3,1)",
      algebra_ok)

symbol_lambda, symbol_tau = sp.symbols("lambda tau")
symbol_top = top(symbol_lambda, symbol_tau)
symbol_direct = polynomial_jacobian(BOTTOM, symbol_top, ETA)
symbol_displacement = affine_displacement(symbol_top)
symbol_translation = symbol_direct[6:10, :] * symbol_displacement[:, 16:20]
symbol_determinant = sp.factor(symbol_translation.det())
symbol_quotient = sp.factor(
    symbol_determinant / (symbol_tau * (symbol_lambda - 1) ** 3)
)
symbolic_ok = bool(
    symbol_determinant != 0 and symbol_quotient.is_number
    and symbol_quotient != 0
    and sp.expand(
        symbol_determinant
        - symbol_quotient * symbol_tau * (symbol_lambda - 1) ** 3
    ) == 0
)
check("the irregular translation determinant has the same exact zero set",
      symbolic_ok, f"det(T)={symbol_determinant}")

boost = sp.eye(4)
boost[1, 1] = sp.Rational(13, 12)
boost[1, 3] = sp.Rational(5, 12)
boost[3, 1] = sp.Rational(5, 12)
boost[3, 3] = sp.Rational(13, 12)
boost_ok = bool(boost.T * ETA * boost == ETA and boost.det() == 1)
check("the adversarial y-time boost is exact and proper Lorentz", boost_ok)

static_stabilizer = stabilizer(ETA, NORMAL)
boosted_normal = boost * NORMAL
boosted_stabilizer = stabilizer(ETA, boosted_normal)
boosted_bottom = tuple(boost * point for point in BOTTOM)

records = []
representatives_ok = True
direct_reconstruction_ok = True
static_ok = True
expanding_ok = True
metric_sign_ok = True
boosted_static_ok = True

for scale, lapse in REPRESENTATIVES:
    points = top(scale, lapse)
    strut_values = tuple(square(points[i], BOTTOM[i], ETA)
                         for i in range(4))
    representatives_ok &= bool(
        len(set(strut_values)) == 1 and strut_values[0] < 0
    )

    result = analyse(BOTTOM, points, ETA)
    negated = analyse(BOTTOM, points, -ETA)
    local_direct = bool(
        result["direct_rank"] == 10
        and result["direct_kernel"].shape == (16, 6)
        and result["lorentz_rank"] == 10
        and result["rigid_displacement_rank"] == 10
        and result["system_rank"] == 14
        and result["parameter_kernel"].shape == (20, 6)
        and result["kernel_equality"]
    )
    direct_reconstruction_ok &= local_direct

    local_metric_sign = bool(
        negated["direct"] == -result["direct"]
        and negated["direct_rank"] == result["direct_rank"]
        and negated["system_rank"] == result["system_rank"]
        and negated["a_image_rank"] == result["a_image_rank"]
        and negated["translation_rank"] == result["translation_rank"]
        and negated["pure_translation_dimension"]
        == result["pure_translation_dimension"]
        and equal_span(negated["direct_kernel"], result["direct_kernel"])
    )
    metric_sign_ok &= local_metric_sign

    if scale == 1:
        local_static = bool(
            result["a_image_rank"] == 3
            and result["translation_rank"] == 1
            and result["pure_translation_dimension"] == 3
            and equal_span(result["a_image"], static_stabilizer)
        )
        transformed_points = tuple(boost * point for point in points)
        transformed = analyse(boosted_bottom, transformed_points, ETA)
        local_boosted_static = bool(
            transformed["kernel_equality"]
            and transformed["a_image_rank"] == 3
            and transformed["translation_rank"] == 1
            and transformed["pure_translation_dimension"] == 3
            and equal_span(transformed["a_image"], boosted_stabilizer)
        )
        static_ok &= local_static
        boosted_static_ok &= local_boosted_static
        local_expanding = False
    else:
        local_expanding = bool(
            result["a_image_rank"] == 6
            and result["translation_rank"] == 4
            and result["pure_translation_dimension"] == 0
            and equal_span(result["a_image"], result["all_lorentz"])
        )
        expanding_ok &= local_expanding
        local_static = False
        local_boosted_static = False

    records.append({
        "scale": scale,
        "lapse": lapse,
        "strut_squared_length": str(strut_values[0]),
        "direct_polynomial_rank": result["direct_rank"],
        "redundant_system_rank": result["system_rank"],
        "parameter_nullity": result["parameter_kernel"].cols,
        "direct_kernel_equality": result["kernel_equality"],
        "a_image_rank": result["a_image_rank"],
        "translation_block_rank": result["translation_rank"],
        "pure_translation_dimension": result["pure_translation_dimension"],
        "static_stabilizer": local_static,
        "boosted_static_stabilizer": local_boosted_static,
        "expanding_full_lorentz": local_expanding,
        "metric_sign_control": local_metric_sign,
    })

check("all irregular representatives have equal timelike struts",
      representatives_ok)
check("the redundant affine-Lorentz kernels equal the polynomial flex kernels",
      direct_reconstruction_ok)
check("the irregular static A-image is the normal stabilizer with 3 translations",
      static_ok)
check("both irregular expanding A-images are all of so(3,1)", expanding_ok)
check("negating the Lorentz metric changes no decision space or rank",
      metric_sign_ok)
check("the independently boosted static image fixes the boosted normal",
      boosted_static_ok)

controls_ok = bool(
    provenance_ok and correction_ok and geometry_ok and algebra_ok
    and symbolic_ok and boost_ok and representatives_ok and metric_sign_ok
)
corroborated = bool(
    controls_ok and direct_reconstruction_ok and static_ok
    and expanding_ok and boosted_static_ok
)
disagreement = bool(
    controls_ok and not corroborated
)

if not controls_ok:
    outcome = "ADVERSARIAL_POINCARE_STRATIFICATION_CONTROL_FAILED"
elif corroborated:
    outcome = "ADVERSARIAL_POINCARE_STRATIFICATION_CORROBORATED"
elif disagreement:
    outcome = "ADVERSARIAL_POINCARE_STRATIFICATION_DISAGREEMENT"
else:
    outcome = "ADVERSARIAL_POINCARE_STRATIFICATION_OPEN"

allowed = {
    "ADVERSARIAL_POINCARE_STRATIFICATION_CONTROL_FAILED",
    "ADVERSARIAL_POINCARE_STRATIFICATION_CORROBORATED",
    "ADVERSARIAL_POINCARE_STRATIFICATION_DISAGREEMENT",
    "ADVERSARIAL_POINCARE_STRATIFICATION_OPEN",
}
check("the adversarial hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "symbolic_translation_determinant": str(symbol_determinant),
    "symbolic_determinant_quotient": str(symbol_quotient),
    "records": records,
    "classification": {
        "local_relative_poincare_stratification": (
            "ADVERSARIALLY CORROBORATED DERIVED EXACT"
            if corroborated else "OPEN"
        ),
        "uniform_lorentz_frame_through_static_slice": "REFUTED",
        "expanding_lorentz_chart": (
            "ADVERSARIALLY CORROBORATED LOCAL"
            if corroborated else "OPEN"
        ),
        "connection_extrinsic_curvature_gluing_action_or_dynamics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("det(T) =", symbol_determinant)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"A-image={record['a_image_rank']}, "
        f"rank(T)={record['translation_block_rank']}, "
        f"pure translations={record['pure_translation_dimension']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
