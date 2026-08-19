#!/usr/bin/env python3
"""Exact Poincare decomposition of cellular tetrahedral-frustum flexes."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_cellular_frustum_relative_poincare.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_relative_poincare_protocol.md"
RIGIDITY_RESULT = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_anisotropic_rigidity_result.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_anisotropic_rigidity.py"
PRIMARY_JSON = HERE / "gravity_600cell_cellular_frustum_anisotropic_rigidity.json"
ADVERSARIAL_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_anisotropic_rigidity_adversarial.py"
ADVERSARIAL_JSON = HERE / "gravity_600cell_cellular_frustum_anisotropic_rigidity_adversarial.json"

PROTOCOL_COMMIT = "f2dc6d3"
EXPECTED_HASHES = {
    "prior_art": "a8811a441fecd137b37085e4018fea7abb3f365750dfc426d16f1b46e5282e7c",
    "protocol": "f88599f2b23d3459f95cba1be12401cd404d98d56a8615f3fa103d1112cc9c7a",
    "rigidity_result": "7221fd948c6e21aa59ea2738d4ef6a224f13b674689a051837366fcfa76f203b",
    "primary_source": "2f766503296aa43f6192d2cce6ce44faac3b7fb57ba131ba0fbf393a2da80f60",
    "primary_json": "c55f98313121018ff5ca1fc834260e8f2f075248a21fd7b99a356d89b2d18255",
    "adversarial_source": "ecc5e0cb5f8913325f00137245f33299c7607b395d219161bf7e0e806068c18a",
    "adversarial_json": "7763287a12075a911134b24e5f23c3c682198923bda1ab8f75ac1d9541540fc1",
}

ETA = sp.diag(1, 1, 1, -1)
BOTTOM = tuple(sp.Matrix(point) for point in (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
))
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
TOP_PAIRS = tuple(combinations(range(4), 2))
ORIGIN_SHIFT = sp.Matrix((2, -1, 3, 4))
VERTEX_PERMUTATION = (0, 2, 3, 1)

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


def kernel_matrix(matrix):
    basis = matrix.nullspace()
    return sp.Matrix.hstack(*basis) if basis else sp.zeros(matrix.cols, 0)


def same_column_space(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


def lorentz_basis():
    result = []
    for left, right in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4, 4)
        generator[left, right] = 1
        generator[right, left] = -ETA[left, left] / ETA[right, right]
        if generator.T * ETA + ETA * generator != sp.zeros(4, 4):
            raise RuntimeError("invalid Lorentz generator")
        result.append(generator)
    return tuple(result)


LORENTZ_BASIS = lorentz_basis()
LORENTZ_FLAT = sp.Matrix.hstack(
    *(sp.Matrix(generator).reshape(16, 1) for generator in LORENTZ_BASIS)
)


def top_points(scale, lapse, bottom=BOTTOM):
    return tuple(sp.Matrix((
        scale * point[0],
        scale * point[1],
        scale * point[2],
        lapse,
    )) for point in bottom)


def top_edge_jacobian(top):
    matrix = sp.zeros(6, 16)
    for row, (left, right) in enumerate(TOP_PAIRS):
        gradient = 2 * ETA * (top[left] - top[right])
        matrix[row, 4 * left:4 * left + 4] = gradient.T
        matrix[row, 4 * right:4 * right + 4] = -gradient.T
    return matrix


def strut_jacobian(bottom, top):
    matrix = sp.zeros(4, 16)
    for index in range(4):
        gradient = 2 * ETA * (top[index] - bottom[index])
        matrix[index, 4 * index:4 * index + 4] = gradient.T
    return matrix


def poincare_design(top):
    columns = []
    for generator in LORENTZ_BASIS:
        columns.append(sp.Matrix.vstack(*(generator * point for point in top)))
    for axis in range(4):
        direction = sp.eye(4)[:, axis]
        columns.append(sp.Matrix.vstack(*(direction for _ in top)))
    return sp.Matrix.hstack(*columns)


def analyse(bottom, top):
    edge = top_edge_jacobian(top)
    strut = strut_jacobian(bottom, top)
    full = edge.col_join(strut)
    design = poincare_design(top)
    constraint = strut * design
    translation = constraint[:, 6:10]
    parameter_kernel = kernel_matrix(constraint)
    displacement_kernel = design * parameter_kernel
    direct_kernel = kernel_matrix(full)
    projection = parameter_kernel[:6, :]
    rotation_projection = projection[:3, :]
    boost_projection = projection[3:6, :]
    return {
        "edge": edge,
        "strut": strut,
        "full": full,
        "design": design,
        "constraint": constraint,
        "translation": translation,
        "parameter_kernel": parameter_kernel,
        "displacement_kernel": displacement_kernel,
        "direct_kernel": direct_kernel,
        "edge_rank": edge.rank(),
        "design_rank": design.rank(),
        "full_rank": full.rank(),
        "constraint_rank": constraint.rank(),
        "translation_rank": translation.rank(),
        "parameter_nullity": parameter_kernel.cols,
        "lorentz_projection_rank": projection.rank(),
        "rotation_projection_rank": rotation_projection.rank(),
        "boost_projection_rank": boost_projection.rank(),
        "pure_translation_dimension": 4 - translation.rank(),
        "edge_annihilation_rank": (edge * design).rank(),
        "design_equals_edge_kernel": same_column_space(
            design, kernel_matrix(edge)
        ),
        "displacement_equals_direct_kernel": same_column_space(
            displacement_kernel, direct_kernel
        ),
    }


def classification_tuple(record):
    return (
        record["edge_rank"],
        record["design_rank"],
        record["full_rank"],
        record["constraint_rank"],
        record["translation_rank"],
        record["parameter_nullity"],
        record["lorentz_projection_rank"],
        record["rotation_projection_rank"],
        record["boost_projection_rank"],
        record["pure_translation_dimension"],
    )


def block_vertex_map(matrix):
    result = sp.zeros(16, 16)
    for index in range(4):
        result[4 * index:4 * index + 4, 4 * index:4 * index + 4] = matrix
    return result


def permutation_map(permutation, block_size):
    result = sp.zeros(block_size * len(permutation))
    identity = sp.eye(block_size)
    for new_index, old_index in enumerate(permutation):
        result[
            block_size * new_index:block_size * new_index + block_size,
            block_size * old_index:block_size * old_index + block_size,
        ] = identity
    return result


def lorentz_coordinates(generator):
    solution, parameters = LORENTZ_FLAT.gauss_jordan_solve(
        sp.Matrix(generator).reshape(16, 1)
    )
    if parameters.rows:
        raise RuntimeError("nonunique Lorentz coordinates")
    return solution


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "rigidity_result": RIGIDITY_RESULT,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "adversarial_source": ADVERSARIAL_SOURCE,
    "adversarial_json": ADVERSARIAL_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all local Poincare inputs have frozen provenance", provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
adversarial = json.loads(ADVERSARIAL_JSON.read_text())
upstream_ok = bool(
    primary["outcome"] == "CELLULAR_FRUSTUM_SIX_SHAPES_UNDERDETERMINED"
    and primary["passed"] == primary["tests"] == 11
    and all(record["fixed_bottom_flex_dimension"] == 6
            for record in primary["records"])
    and adversarial["outcome"]
    == "ADVERSARIAL_CELLULAR_FRUSTUM_SIX_SHAPES_CORROBORATED"
    and adversarial["passed"] == adversarial["tests"] == 9
    and all(record["base_nullity"] == 6
            for record in adversarial["records"])
)
check("both independent upstream six-flex certificates persist", upstream_ok)

basis_ok = bool(
    len(LORENTZ_BASIS) == 6
    and LORENTZ_FLAT.rank() == 6
    and all(generator.T * ETA + ETA * generator == sp.zeros(4, 4)
            for generator in LORENTZ_BASIS)
)
check("the disclosed rotation/boost basis spans exact so(3,1)", basis_ok)

symbol_lambda, symbol_tau = sp.symbols("lambda tau")
symbol_top = top_points(symbol_lambda, symbol_tau)
symbol_translation = analyse(BOTTOM, symbol_top)["translation"]
symbol_determinant = sp.factor(symbol_translation.det())
determinant_quotient = sp.factor(
    symbol_determinant / (symbol_tau * (symbol_lambda - 1) ** 3)
)
determinant_ok = bool(
    symbol_determinant != 0
    and determinant_quotient.is_number
    and determinant_quotient != 0
    and sp.expand(
        symbol_determinant
        - determinant_quotient * symbol_tau * (symbol_lambda - 1) ** 3
    ) == 0
)
check(
    "the translation determinant exposes the static cubic degeneration",
    determinant_ok,
    f"det(T)={symbol_determinant}",
)

boost = sp.eye(4)
boost[0, 0] = sp.Rational(5, 4)
boost[0, 3] = sp.Rational(3, 4)
boost[3, 0] = sp.Rational(3, 4)
boost[3, 3] = sp.Rational(5, 4)
boost_ok = boost.T * ETA * boost == ETA and boost.det() == 1
check("the covariance control is an exact proper Lorentz boost", boost_ok)

records = []
rigid_motion_controls = True
constraint_controls = True
static_prediction = True
expanding_prediction = True
graph_controls = True
origin_covariance = True
lorentz_covariance = True
permutation_covariance = True

for scale, lapse in REPRESENTATIVES:
    top = top_points(scale, lapse)
    record = analyse(BOTTOM, top)

    rigid_motion_controls &= bool(
        record["edge_rank"] == 6
        and record["design_rank"] == 10
        and record["edge_annihilation_rank"] == 0
        and record["design_equals_edge_kernel"]
        and record["full_rank"] == 10
    )
    constraint_controls &= bool(
        record["constraint_rank"] == 4
        and record["parameter_nullity"] == 6
        and record["displacement_equals_direct_kernel"]
        and record["direct_kernel"].cols == 6
    )

    if scale == 1:
        local_static = bool(
            record["translation_rank"] == 1
            and record["pure_translation_dimension"] == 3
            and record["lorentz_projection_rank"] == 3
            and record["rotation_projection_rank"] == 3
            and record["boost_projection_rank"] == 0
        )
        static_prediction &= local_static
        local_graph = False
    else:
        graph = sp.eye(6).col_join(
            -record["translation"].inv() * record["constraint"][:, :6]
        )
        local_graph = bool(
            record["translation_rank"] == 4
            and record["pure_translation_dimension"] == 0
            and record["lorentz_projection_rank"] == 6
            and record["constraint"] * graph == sp.zeros(4, 6)
            and same_column_space(graph, record["parameter_kernel"])
        )
        expanding_prediction &= local_graph
        graph_controls &= local_graph

    shifted_bottom = tuple(point + ORIGIN_SHIFT for point in BOTTOM)
    shifted_top = tuple(point + ORIGIN_SHIFT for point in top)
    shifted = analyse(shifted_bottom, shifted_top)
    shift_parameter_map = sp.eye(10)
    for column, generator in enumerate(LORENTZ_BASIS):
        shift_parameter_map[6:10, column] = -generator * ORIGIN_SHIFT
    local_origin_covariance = bool(
        shifted["design"] * shift_parameter_map == record["design"]
        and shifted["constraint"] * shift_parameter_map
        == record["constraint"]
        and classification_tuple(shifted) == classification_tuple(record)
    )
    origin_covariance &= local_origin_covariance

    boosted_bottom = tuple(boost * point for point in BOTTOM)
    boosted_top = tuple(boost * point for point in top)
    boosted = analyse(boosted_bottom, boosted_top)
    lorentz_parameter_map = sp.zeros(10, 10)
    for column, generator in enumerate(LORENTZ_BASIS):
        transformed = boost * generator * boost.inv()
        lorentz_parameter_map[:6, column] = lorentz_coordinates(transformed)
    lorentz_parameter_map[6:10, 6:10] = boost
    vertex_boost = block_vertex_map(boost)
    local_lorentz_covariance = bool(
        boosted["design"] * lorentz_parameter_map
        == vertex_boost * record["design"]
        and boosted["constraint"] * lorentz_parameter_map
        == record["constraint"]
        and classification_tuple(boosted) == classification_tuple(record)
    )
    lorentz_covariance &= local_lorentz_covariance

    permuted_bottom = tuple(BOTTOM[index] for index in VERTEX_PERMUTATION)
    permuted_top = tuple(top[index] for index in VERTEX_PERMUTATION)
    permuted = analyse(permuted_bottom, permuted_top)
    vertex_permutation = permutation_map(VERTEX_PERMUTATION, 4)
    strut_permutation = permutation_map(VERTEX_PERMUTATION, 1)
    local_permutation_covariance = bool(
        permuted["design"] == vertex_permutation * record["design"]
        and permuted["constraint"] == strut_permutation * record["constraint"]
        and classification_tuple(permuted) == classification_tuple(record)
    )
    permutation_covariance &= local_permutation_covariance

    records.append({
        "scale": scale,
        "lapse": lapse,
        "edge_rank": record["edge_rank"],
        "design_rank": record["design_rank"],
        "full_constraint_rank": record["full_rank"],
        "poincare_constraint_rank": record["constraint_rank"],
        "poincare_constraint_nullity": record["parameter_nullity"],
        "translation_block_rank": record["translation_rank"],
        "pure_translation_flex_dimension": record["pure_translation_dimension"],
        "lorentz_projection_rank": record["lorentz_projection_rank"],
        "rotation_projection_rank": record["rotation_projection_rank"],
        "boost_projection_rank": record["boost_projection_rank"],
        "is_lorentz_graph": local_graph,
        "origin_covariance": local_origin_covariance,
        "lorentz_covariance": local_lorentz_covariance,
        "vertex_permutation_covariance": local_permutation_covariance,
    })

check("all top-edge flexes are exactly restricted Poincare motions",
      rigid_motion_controls)
check("the four struts cut the Poincare space to the direct six-flex kernel",
      constraint_controls)
check("the static flexes are exactly rotations plus spatial translations",
      static_prediction)
check("both expanding flex spaces are exact graphs over so(3,1)",
      expanding_prediction and graph_controls)
check("the classification is covariant under an exact origin shift",
      origin_covariance)
check("the classification is covariant under a rational Lorentz boost",
      lorentz_covariance)
check("the classification survives simultaneous paired-vertex relabelling",
      permutation_covariance)

controls_ok = bool(
    provenance_ok and upstream_ok and basis_ok and determinant_ok and boost_ok
    and rigid_motion_controls and constraint_controls
    and origin_covariance and lorentz_covariance and permutation_covariance
)
uniform = bool(
    controls_ok
    and all(record["lorentz_projection_rank"] == 6
            and record["pure_translation_flex_dimension"] == 0
            for record in records)
)
stratified = bool(
    controls_ok and static_prediction and expanding_prediction
)
refuted = bool(
    controls_ok and any(
        record["scale"] != 1
        and (record["lorentz_projection_rank"] != 6
             or record["pure_translation_flex_dimension"] != 0)
        for record in records
    )
)

if not controls_ok:
    outcome = "RELATIVE_POINCARE_CONTROL_FAILED"
elif uniform:
    outcome = "UNIFORM_RELATIVE_LORENTZ_CHART"
elif stratified:
    outcome = "STRATIFIED_RELATIVE_LORENTZ_CHART"
elif refuted:
    outcome = "RELATIVE_LORENTZ_INTERPRETATION_REFUTED"
else:
    outcome = "RELATIVE_LORENTZ_INTERPRETATION_OPEN"

allowed = {
    "RELATIVE_POINCARE_CONTROL_FAILED",
    "UNIFORM_RELATIVE_LORENTZ_CHART",
    "STRATIFIED_RELATIVE_LORENTZ_CHART",
    "RELATIVE_LORENTZ_INTERPRETATION_REFUTED",
    "RELATIVE_LORENTZ_INTERPRETATION_OPEN",
}
check("the preregistered hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "metric": [int(ETA[index, index]) for index in range(4)],
    "lorentz_basis": ["J01", "J02", "J12", "K03", "K13", "K23"],
    "symbolic_translation_determinant": str(symbol_determinant),
    "symbolic_determinant_quotient": str(determinant_quotient),
    "records": records,
    "classification": {
        "six_flexes_equal_uniform_lorentz_frame": (
            "DERIVED EXACT" if uniform else "REFUTED DERIVED EXACT"
        ),
        "expanding_local_lorentz_chart": (
            "DERIVED EXACT" if expanding_prediction and controls_ok else "OPEN"
        ),
        "static_flex_content": (
            "THREE ROTATIONS PLUS THREE RELATIVE SPATIAL TRANSLATIONS"
            if static_prediction and controls_ok else "OPEN"
        ),
        "physical_connection_extrinsic_curvature_or_action": "NOT TESTED",
        "global_gluing_closure_and_shape_matching": "OPEN",
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
        f"rank(T)={record['translation_block_rank']}, "
        f"pure translations={record['pure_translation_flex_dimension']}, "
        f"Lorentz projection={record['lorentz_projection_rank']}, "
        f"rotation/boost={record['rotation_projection_rank']}/"
        f"{record['boost_projection_rank']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
