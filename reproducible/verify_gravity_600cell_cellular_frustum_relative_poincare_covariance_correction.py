#!/usr/bin/env python3
"""Invariant correction of the static Poincare covariance classifier."""

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_cellular_frustum_relative_poincare_covariance_correction.json"
ORIGINAL_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_protocol.md"
ORIGINAL_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_relative_poincare.py"
FAILED_JSON = HERE / "gravity_600cell_cellular_frustum_relative_poincare.json"
FAILURE_NOTE = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_covariance_failure.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_covariance_correction_protocol.md"

PROTOCOL_COMMIT = "5bd3b5b"
EXPECTED_HASHES = {
    "original_protocol": "f88599f2b23d3459f95cba1be12401cd404d98d56a8615f3fa103d1112cc9c7a",
    "original_source": "308d97cc0b057d3ac79cbc8a4706a63fe1b3d76792d84d8576a84df7d7d63514",
    "failed_json": "3ac5cce9db2b2f828e0ced2114f301f761dd9371847b712ad47119709396cf7d",
    "failure_note": "de81078a7a44af6c49461b0d22d3a59a4edff0e2e4c29432529e00c9f4e39109",
    "protocol": "ad569b4c4ecbfb4b6d3db7d0225dd26f0426e7e973ddb9c82c21daff5b404b3c",
}

ETA = sp.diag(1, 1, 1, -1)
P = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
NORMAL = sp.Matrix((0, 0, 0, 1))

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


def kernel(matrix):
    vectors = matrix.nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.cols, 0)


def same_image(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


def make_generators():
    generators = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        matrix = sp.zeros(4)
        matrix[a, b] = 1
        matrix[b, a] = -ETA[a, a] / ETA[b, b]
        generators.append(matrix)
    return tuple(generators)


GENERATORS = make_generators()
GENERATOR_COORDINATES = sp.Matrix.hstack(
    *(matrix.reshape(16, 1) for matrix in GENERATORS)
)


def coordinates(matrix):
    solution, free = GENERATOR_COORDINATES.gauss_jordan_solve(
        matrix.reshape(16, 1)
    )
    if free.rows:
        raise RuntimeError("Lorentz basis coordinate ambiguity")
    return solution


def upper(scale, lapse):
    return tuple(sp.Matrix((
        scale * point[0], scale * point[1], scale * point[2], lapse
    )) for point in P)


def displacement_design(points):
    columns = [sp.Matrix.vstack(*(generator * q for q in points))
               for generator in GENERATORS]
    for axis in range(4):
        vector = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(vector for _ in points)))
    return sp.Matrix.hstack(*columns)


def strut_derivative(bottom, top):
    derivative = sp.zeros(4, 16)
    for index, (p, q) in enumerate(zip(bottom, top)):
        derivative[index, 4 * index:4 * index + 4] = (
            2 * ETA * (q - p)
        ).T
    return derivative


def data(bottom, top):
    design = displacement_design(top)
    strut = strut_derivative(bottom, top)
    constraint = strut * design
    constraint_kernel = kernel(constraint)
    lorentz_image = constraint_kernel[:6, :]
    translation = constraint[:, 6:10]
    return {
        "design": design,
        "constraint": constraint,
        "kernel": constraint_kernel,
        "lorentz_image": lorentz_image,
        "lorentz_rank": lorentz_image.rank(),
        "translation_rank": translation.rank(),
        "pure_translation_dimension": 4 - translation.rank(),
        "coordinate_rotation_rank": lorentz_image[:3, :].rank(),
        "coordinate_boost_rank": lorentz_image[3:6, :].rank(),
    }


def stabilizer(normal):
    action = sp.Matrix.hstack(*(generator * normal for generator in GENERATORS))
    return kernel(action)


def repeated_map(matrix):
    result = sp.zeros(16)
    for index in range(4):
        result[4 * index:4 * index + 4, 4 * index:4 * index + 4] = matrix
    return result


paths = {
    "original_protocol": ORIGINAL_PROTOCOL,
    "original_source": ORIGINAL_SOURCE,
    "failed_json": FAILED_JSON,
    "failure_note": FAILURE_NOTE,
    "protocol": PROTOCOL,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("the failed execution and correction protocol have frozen provenance",
      provenance_ok, str(hashes))

failed = json.loads(FAILED_JSON.read_text())
failed_records = failed["records"]
failed_preserved = bool(
    failed["outcome"] == "RELATIVE_POINCARE_CONTROL_FAILED"
    and failed["passed"] == 12 and failed["tests"] == 13
    and len(failed_records) == 3
    and [record["lorentz_covariance"] for record in failed_records]
    == [False, True, True]
)
check("the original 12/13 static-only covariance failure is preserved",
      failed_preserved)

basis_ok = bool(
    GENERATOR_COORDINATES.rank() == 6
    and all(A.T * ETA + ETA * A == sp.zeros(4) for A in GENERATORS)
)
check("the independent six-generator basis is exact so(3,1)", basis_ok)

L = sp.eye(4)
L[0, 0] = sp.Rational(5, 4)
L[0, 3] = sp.Rational(3, 4)
L[3, 0] = sp.Rational(3, 4)
L[3, 3] = sp.Rational(5, 4)
boost_ok = bool(L.T * ETA * L == ETA and L.det() == 1)
check("the correction uses the same exact proper Lorentz boost", boost_ok)

adjoint = sp.Matrix.hstack(
    *(coordinates(L * generator * L.inv()) for generator in GENERATORS)
)
parameter_transport = sp.zeros(10)
parameter_transport[:6, :6] = adjoint
parameter_transport[6:10, 6:10] = L
vertex_transport = repeated_map(L)
transport_ok = bool(
    adjoint.rank() == 6 and parameter_transport.det() != 0
)
check("adjoint and Poincare parameter transports are exact isomorphisms",
      transport_ok)

original_stabilizer = stabilizer(NORMAL)
boosted_normal = L * NORMAL
boosted_stabilizer = stabilizer(boosted_normal)
abstract_stabilizer_covariance = same_image(
    boosted_stabilizer, adjoint * original_stabilizer
)
check("observer-normal stabilizers transform by exact Lorentz conjugation",
      abstract_stabilizer_covariance)

records = []
intertwining_ok = True
kernel_transport_ok = True
invariant_ranks_ok = True
static_stabilizer_ok = True
noninvariant_split_control = True
expanding_graph_ok = True

boosted_bottom = tuple(L * point for point in P)
for scale, lapse in REPRESENTATIVES:
    top = upper(scale, lapse)
    transformed_top = tuple(L * point for point in top)
    original = data(P, top)
    transformed = data(boosted_bottom, transformed_top)

    local_intertwining = bool(
        transformed["design"] * parameter_transport
        == vertex_transport * original["design"]
        and transformed["constraint"] * parameter_transport
        == original["constraint"]
    )
    local_kernel_transport = same_image(
        transformed["kernel"], parameter_transport * original["kernel"]
    )
    local_invariant_ranks = bool(
        transformed["lorentz_rank"] == original["lorentz_rank"]
        and transformed["translation_rank"] == original["translation_rank"]
        and transformed["pure_translation_dimension"]
        == original["pure_translation_dimension"]
    )
    intertwining_ok &= local_intertwining
    kernel_transport_ok &= local_kernel_transport
    invariant_ranks_ok &= local_invariant_ranks

    original_split = (
        original["coordinate_rotation_rank"],
        original["coordinate_boost_rank"],
    )
    transformed_split = (
        transformed["coordinate_rotation_rank"],
        transformed["coordinate_boost_rank"],
    )

    if scale == 1:
        local_static = bool(
            original["lorentz_rank"] == transformed["lorentz_rank"] == 3
            and original["pure_translation_dimension"]
            == transformed["pure_translation_dimension"] == 3
            and same_image(original["lorentz_image"], original_stabilizer)
            and same_image(transformed["lorentz_image"], boosted_stabilizer)
            and same_image(
                transformed["lorentz_image"],
                adjoint * original["lorentz_image"],
            )
        )
        local_split_control = bool(
            original_split != transformed_split
            and transformed["coordinate_boost_rank"] > 0
        )
        static_stabilizer_ok &= local_static
        noninvariant_split_control &= local_split_control
        local_expanding = False
    else:
        local_expanding = bool(
            original["lorentz_rank"] == transformed["lorentz_rank"] == 6
            and original["translation_rank"]
            == transformed["translation_rank"] == 4
            and original["pure_translation_dimension"]
            == transformed["pure_translation_dimension"] == 0
        )
        expanding_graph_ok &= local_expanding
        local_static = False
        local_split_control = True

    records.append({
        "scale": scale,
        "lapse": lapse,
        "intertwining_exact": local_intertwining,
        "kernel_transport_exact": local_kernel_transport,
        "invariant_ranks_equal": local_invariant_ranks,
        "original_lorentz_rank": original["lorentz_rank"],
        "transformed_lorentz_rank": transformed["lorentz_rank"],
        "original_pure_translation_dimension": original["pure_translation_dimension"],
        "transformed_pure_translation_dimension": transformed["pure_translation_dimension"],
        "original_coordinate_rotation_boost_ranks": list(original_split),
        "transformed_coordinate_rotation_boost_ranks": list(transformed_split),
        "static_observer_stabilizer": local_static,
        "noninvariant_coordinate_split_control": local_split_control,
        "expanding_full_lorentz_graph": local_expanding,
    })

check("all displacement and constraint matrices intertwine exactly",
      intertwining_ok)
check("all Poincare kernels transport exactly under the boost",
      kernel_transport_ok)
check("all invariant projection and translation ranks are preserved",
      invariant_ranks_ok)
check("the static Lorentz images equal the corresponding observer stabilizers",
      static_stabilizer_ok)
check("the old coordinate rotation/boost split changes under the boost",
      noninvariant_split_control)
check("both expanding strata remain full graphs over so(3,1)",
      expanding_graph_ok)

controls_ok = bool(
    provenance_ok and failed_preserved and basis_ok and boost_ok
    and transport_ok and abstract_stabilizer_covariance
)
real_disagreement = bool(
    controls_ok and (not intertwining_ok or not kernel_transport_ok)
)
corroborated = bool(
    controls_ok and intertwining_ok and kernel_transport_ok
    and invariant_ranks_ok and static_stabilizer_ok
    and noninvariant_split_control and expanding_graph_ok
)

if not controls_ok:
    outcome = "POINCARE_COVARIANCE_CORRECTION_CONTROL_FAILED"
elif real_disagreement:
    outcome = "POINCARE_COVARIANCE_REAL_DISAGREEMENT"
elif corroborated:
    outcome = "STATIC_STABILIZER_AND_EXPANDING_LORENTZ_CHART_CORROBORATED"
else:
    outcome = "POINCARE_COVARIANCE_CORRECTION_OPEN"

allowed = {
    "POINCARE_COVARIANCE_CORRECTION_CONTROL_FAILED",
    "POINCARE_COVARIANCE_REAL_DISAGREEMENT",
    "STATIC_STABILIZER_AND_EXPANDING_LORENTZ_CHART_CORROBORATED",
    "POINCARE_COVARIANCE_CORRECTION_OPEN",
}
check("the correction hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "preserved_original_outcome": failed["outcome"],
    "records": records,
    "classification": {
        "original_static_covariance_failure": "NONINVARIANT CLASSIFIER",
        "relative_poincare_kernel": (
            "DERIVED EXACT" if corroborated else "OPEN"
        ),
        "static_lorentz_image": (
            "OBSERVER-NORMAL STABILIZER, DIMENSION 3"
            if corroborated else "OPEN"
        ),
        "static_pure_translation_sector": (
            "DIMENSION 3" if corroborated else "OPEN"
        ),
        "expanding_lorentz_graph": (
            "DERIVED EXACT LOCAL" if corroborated else "OPEN"
        ),
        "uniform_lorentz_frame_through_static_slice": "REFUTED",
        "connection_extrinsic_curvature_gluing_or_dynamics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"Lorentz rank={record['original_lorentz_rank']}, "
        f"pure translations={record['original_pure_translation_dimension']}, "
        f"coordinate split "
        f"{record['original_coordinate_rotation_boost_ranks']} -> "
        f"{record['transformed_coordinate_rotation_boost_ranks']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
