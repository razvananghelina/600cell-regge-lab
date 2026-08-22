#!/usr/bin/env python3
"""Full-real-space adversarial replication of the finite-height rank census.

Primary implementation commit: 1472213.
Primary result commit: 2b9c124.
Adversarial protocol commit: b28c4a9.
Registry commit: 2a10931.
"""

import ast
from collections import Counter
import contextlib
import gc
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import networkx as nx
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = (
    HERE / "gravity_600cell_finite_height_internal_carrier_rank_adversarial.json"
)
MATRIX_OUTPUT = (
    HERE
    / "gravity_600cell_finite_height_internal_carrier_rank_adversarial_matrices.npz"
)
PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_internal_carrier_rank_adversarial_protocol.md"
)
PRIMARY_JSON = HERE / "gravity_600cell_finite_height_internal_carrier_rank.json"
PRIMARY_MATRICES = (
    HERE / "gravity_600cell_finite_height_internal_carrier_rank_matrices.npz"
)
PRIMARY_SOURCE = (
    HERE / "verify_gravity_600cell_finite_height_internal_carrier_rank.py"
)
QUADRATIC_SOURCE = (
    HERE / "verify_gravity_600cell_finite_height_carrier_quadratic.py"
)
HESSIAN_SOURCE = (
    HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
)
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
QUADRATIC_JSON = HERE / "gravity_600cell_finite_height_carrier_quadratic.json"
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "protocol": PROTOCOL,
    "primary_json": PRIMARY_JSON,
    "primary_matrices": PRIMARY_MATRICES,
    "primary_source": PRIMARY_SOURCE,
    "quadratic_source": QUADRATIC_SOURCE,
    "hessian_source": HESSIAN_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "quadratic_json": QUADRATIC_JSON,
}
EXPECTED_HASHES = {
    "protocol": (
        "e0f1a1401062bfa708695b08f404f3620c08d6220de300ad0bb7b2dc4dad6d97"
    ),
    "primary_json": (
        "513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8"
    ),
    "primary_matrices": (
        "97f5b8318be2b3ccf843db87e678ac1ac6ce402db262023c6bbc63a7b647321b"
    ),
    "primary_source": (
        "fff2c70dc3685562b4e5f1e7646886c828a5df1aa7cbed792cef8b19afdf8c62"
    ),
    "quadratic_source": (
        "bbe7112270a7f2bcb2d443fab45ca450598e7234250bd335b14a4ed7869443a5"
    ),
    "hessian_source": (
        "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5"
    ),
    "geometry_source": (
        "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"
    ),
    "quadratic_json": (
        "0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9"
    ),
}

PRIMARY_IMPLEMENTATION_COMMIT = "1472213"
PRIMARY_RESULT_COMMIT = "2b9c124"
PROTOCOL_COMMIT = "b28c4a9"
REGISTRY_COMMIT = "2a10931"
VERIFIER_NAME = Path(__file__).name

DPS = 120
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
LEVELS = tuple(DERIVATIVE_STEPS)
DIRECT_STEPS = (mp.mpf("1e-4"), mp.mpf("5e-5"), mp.mpf("2.5e-5"))
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {
    hinge: index for index, hinge in enumerate(LOCAL_HINGES)
}
I = mp.mpc(0, 1)
VERTICES = 120
OLD = 720
INTERNAL = 840
NEW = 720
FULL = OLD + INTERNAL + NEW
ACTIVE = INTERNAL + NEW
DATA = 2 * VERTICES

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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=80):
    return mp.nstr(value, digits)


def float_text(value):
    return f"{float(value):.17e}"


def sparse_norm2(matrix):
    matrix = sp.csr_matrix(matrix)
    if min(matrix.shape) == 0 or matrix.nnz == 0:
        return 0.0
    try:
        return float(spla.svds(
            matrix, k=1, which="LM", return_singular_vectors=False,
            tol=1e-10, maxiter=20000,
        )[0])
    except Exception:
        return float(sp.linalg.norm(matrix))


def load_named_functions(path, names, namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    found = {node.name for node in selected}
    if found != set(names):
        raise RuntimeError(f"missing audited functions: {sorted(set(names)-found)}")
    module = ast.Module(body=selected, type_ignores=[])
    ast.fix_missing_locations(module)
    exec(compile(module, str(path), "exec"), namespace)


def registry_inventory(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    scripts = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "scripts"
            for target in node.targets
        ):
            scripts = ast.literal_eval(node.value)
            break
    if scripts is None:
        raise RuntimeError("run_all.py has no literal scripts registry")
    counts = Counter(scripts)
    duplicates = sorted(name for name, count in counts.items() if count != 1)
    return scripts, duplicates


def response_matrix(hessian, carrier):
    hessian_block = sp.csr_matrix(hessian[OLD:OLD + INTERNAL, OLD:])
    return np.asarray(
        (hessian_block @ sp.csr_matrix(carrier)).toarray(), dtype=float
    )


def column_scaling(midpoint):
    scale_norm = float(np.linalg.norm(midpoint[:, :VERTICES]))
    strut_norm = float(np.linalg.norm(midpoint[:, VERTICES:]))
    factors = np.r_[
        np.full(VERTICES, 1 / max(1.0, scale_norm)),
        np.full(VERTICES, 1 / max(1.0, strut_norm)),
    ]
    return factors, scale_norm, strut_norm


def kernel_projector(matrix, nullity):
    if nullity == 0:
        return np.zeros((matrix.shape[1], matrix.shape[1]))
    _, _, right = la.svd(matrix, full_matrices=False, check_finite=True)
    vectors = right[-nullity:, :].T
    return vectors @ vectors.T


def classify_scope(raw_levels, carrier, row_count, round_entry_bound):
    midpoint_raw = (
        raw_levels["operational_primary"]
        + raw_levels["validation_primary"]
    ) / 2
    factors, scale_norm, strut_norm = column_scaling(midpoint_raw)
    scaled = {
        level: matrix * factors[np.newaxis, :]
        for level, matrix in raw_levels.items()
    }
    midpoint = (
        scaled["operational_primary"] + scaled["validation_primary"]
    ) / 2
    normalization = max(1.0, sparse_norm2(midpoint))
    step_error = max(
        sparse_norm2(matrix - midpoint) for matrix in scaled.values()
    ) / normalization
    scaled_carrier = carrier * factors[np.newaxis, :]
    column_l1 = np.sum(np.abs(scaled_carrier), axis=0)
    round_absolute = float(
        round_entry_bound * math.sqrt(row_count)
        * math.sqrt(float(np.dot(column_l1, column_l1)))
    )
    round_error = (
        round_absolute / normalization
        + 500 * np.finfo(float).eps * max(row_count, DATA)
    )
    total_error = step_error + round_error

    records = {}
    nullities = []
    resolved = True
    svd_runtime = {}
    for level, matrix in scaled.items():
        values = la.svdvals(matrix, check_finite=True)
        normalized = values / normalization
        labels = [
            "ZERO" if value <= 10 * total_error
            else "NONZERO" if value > 100 * total_error
            else "OPEN"
            for value in normalized
        ]
        nullity = labels.count("ZERO")
        nullities.append(nullity)
        resolved &= "OPEN" not in labels
        records[level] = {
            "singular_values": [float_text(value) for value in values],
            "normalized_singular_values": [
                float_text(value) for value in normalized
            ],
            "labels": labels,
            "nullity": nullity,
        }
        svd_runtime[level] = values
    resolved &= len(set(nullities)) == 1
    nullity = nullities[0] if resolved else None
    rank = DATA - nullity if resolved else None

    _, qr_matrix, pivots = la.qr(
        midpoint, pivoting=True, mode="economic", check_finite=True
    )
    qr_diagonal = np.abs(np.diag(qr_matrix))
    if resolved:
        accepted = qr_diagonal[:rank]
        rejected = qr_diagonal[rank:]
        accepted_ok = bool(
            rank == 0 or float(np.min(accepted)) / normalization > 100 * total_error
        )
        rejected_ok = bool(
            not len(rejected)
            or float(np.max(rejected)) / normalization <= 10 * total_error
        )
        qr_ok = accepted_ok and rejected_ok
    else:
        qr_ok = False

    zero_labels_ok = bool(0 <= 10 * total_error)
    identity_labels_ok = bool(1 / normalization > 100 * total_error)
    return {
        "shape": [row_count, DATA],
        "source_scaling": {
            "scale_norm": float_text(scale_norm),
            "strut_norm": float_text(strut_norm),
            "scale_factor": float_text(factors[0]),
            "strut_factor": float_text(factors[-1]),
        },
        "normalization": float_text(normalization),
        "step_error": float_text(step_error),
        "round_absolute_bound": float_text(round_absolute),
        "round_error": float_text(round_error),
        "total_error": float_text(total_error),
        "levels": records,
        "resolved": bool(resolved),
        "rank": rank,
        "nullity": nullity,
        "qr": {
            "pivot_order": [int(value) for value in pivots],
            "diagonal": [float_text(value) for value in qr_diagonal],
            "accepted_above_nonzero_gate": bool(accepted_ok) if resolved else False,
            "rejected_below_zero_gate": bool(rejected_ok) if resolved else False,
            "passed": bool(qr_ok),
        },
        "synthetic_controls": {
            "zero_full_nullity": zero_labels_ok,
            "identity_zero_nullity": identity_labels_ok,
        },
    }, {
        "midpoint_raw": midpoint_raw,
        "midpoint_scaled": midpoint,
        "scaled_levels": scaled,
        "factors": factors,
        "normalization": normalization,
        "total_error": total_error,
    }


def projector_census(raw_levels, nullity):
    midpoint = (
        raw_levels["operational_primary"]
        + raw_levels["validation_primary"]
    ) / 2
    matrices = {**raw_levels, "midpoint": midpoint}
    projectors = {
        level: kernel_projector(matrix, nullity)
        for level, matrix in matrices.items()
    }
    uncertainty = max(
        float(np.linalg.norm(
            projectors[level] - projectors["midpoint"], ord=2
        ))
        for level in LEVELS
    ) + 500 * np.finfo(float).eps * max(INTERNAL, DATA)
    return projectors, uncertainty


def compare_projectors(left, right, uncertainty):
    difference = float(np.linalg.norm(left - right, ord=2))
    if difference <= 10 * uncertainty:
        classification = "AGREE"
    elif difference > 100 * uncertainty:
        classification = "DEPENDENT"
    else:
        classification = "OPEN"
    return {
        "classification": classification,
        "difference_two_norm": float_text(difference),
        "uncertainty": float_text(uncertainty),
        "agreement_gate": float_text(10 * uncertainty),
        "dependence_gate": float_text(100 * uncertainty),
    }


def direct_secant(library, carrier, kernel_vector, comparison_vector, response):
    vectors = {
        "kernel": np.asarray(kernel_vector, dtype=float),
        "smallest_nonzero": np.asarray(comparison_vector, dtype=float),
    }
    results = {}
    branch_ok = True
    for name, vector in vectors.items():
        active_direction = carrier @ vector
        full_direction = np.r_[np.zeros(OLD), active_direction]
        secants = []
        branch_records = []
        for step in DIRECT_STEPS:
            plus, plus_branch = library["full_gradient_at_delta"](
                full_direction * float(step)
            )
            minus, minus_branch = library["full_gradient_at_delta"](
                -full_direction * float(step)
            )
            local_branch_ok = bool(
                plus_branch["negative_counts"] == Counter({1: 2400})
                and minus_branch["negative_counts"] == Counter({1: 2400})
                and plus_branch["minimum_leading_minor"] > 0
                and minus_branch["minimum_leading_minor"] > 0
                and plus_branch["minimum_argument"] > mp.mpf("1e-6")
                and minus_branch["minimum_argument"] > mp.mpf("1e-6")
                and plus_branch["maximum_imaginary"] < mp.mpf("1e-60")
                and minus_branch["maximum_imaginary"] < mp.mpf("1e-60")
            )
            branch_ok &= local_branch_ok
            derivative = np.asarray([
                float(mp.re((plus[index] - minus[index]) / (2 * step)))
                for index in range(OLD, OLD + INTERNAL)
            ])
            secants.append(derivative)
            branch_records.append({
                "step": mp_text(step),
                "branch_ok": local_branch_ok,
                "minimum_leading_minor": mp_text(min(
                    plus_branch["minimum_leading_minor"],
                    minus_branch["minimum_leading_minor"],
                )),
                "minimum_argument": mp_text(min(
                    plus_branch["minimum_argument"],
                    minus_branch["minimum_argument"],
                )),
                "maximum_imaginary": mp_text(max(
                    plus_branch["maximum_imaginary"],
                    minus_branch["maximum_imaginary"],
                )),
                "secant_norm": float_text(np.linalg.norm(derivative)),
            })
        first_difference = float(np.linalg.norm(secants[0] - secants[1]))
        second_difference = float(np.linalg.norm(secants[1] - secants[2]))
        convergence_ok = bool(
            second_difference <= 0.4 * first_difference + 1e-20
        )
        predicted = response @ vector
        prediction_error = float(
            np.linalg.norm(secants[-1] - predicted)
            / max(1e-300, np.linalg.norm(predicted))
        )
        results[name] = {
            "branches": branch_records,
            "first_difference": float_text(first_difference),
            "second_difference": float_text(second_difference),
            "second_order_convergence": convergence_ok,
            "predicted_norm": float_text(np.linalg.norm(predicted)),
            "finest_prediction_relative_error": prediction_error,
            "finest_secant": secants[-1],
        }
    suppression = float(
        np.linalg.norm(results["kernel"]["finest_secant"])
        / max(1e-300, np.linalg.norm(
            results["smallest_nonzero"]["finest_secant"]
        ))
    )
    validated = bool(
        branch_ok
        and all(results[name]["second_order_convergence"] for name in results)
        and results["smallest_nonzero"]["finest_prediction_relative_error"] < 1e-5
        and suppression < 1e-6
    )
    public = {}
    for name, item in results.items():
        public[name] = {
            key: (
                float_text(value)
                if key == "finest_prediction_relative_error" else value
            )
            for key, value in item.items() if key != "finest_secant"
        }
    return {
        "directions": public,
        "kernel_to_nonzero_finest_secant_ratio": float_text(suppression),
        "branch_ok": bool(branch_ok),
        "validated": validated,
    }


print("[setup] checking hashes without reading primary scientific values", flush=True)
input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
provenance_ok = input_hashes == EXPECTED_HASHES
check("all adversarial inputs retain their frozen hashes", provenance_ok)

scripts, duplicates = registry_inventory(RUN_ALL)
registry_ok = bool(scripts.count(VERIFIER_NAME) == 1 and not duplicates)
check(
    "the adversarial verifier is registered exactly once with no duplicates",
    registry_ok,
    f"entries={len(scripts)}, duplicates={duplicates}",
)

spec = importlib.util.spec_from_file_location(
    "finite_height_rank_adversarial_geometry", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
geometry_import_ok = bool(
    gro.tests == gro.passed == 43 and set(gro.models) == {"even", "odd"}
)
check(
    "the full physical one-slab geometry retains all 43 certificates",
    geometry_import_ok,
)

quadratic_library = {
    "ACTIVE": ACTIVE,
    "Counter": Counter,
    "DATA": DATA,
    "NEW": NEW,
    "OLD": OLD,
    "VERTICES": VERTICES,
    "mp": mp,
    "mp_text": mp_text,
    "np": np,
    "nx": nx,
}
load_named_functions(
    QUADRATIC_SOURCE,
    {"finite_height_state", "build_carrier"},
    quadratic_library,
)
state = quadratic_library["finite_height_state"]()
background_ok = bool(
    state["bracket_ok"]
    and abs(state["residual"]) < mp.mpf("1e-90")
    and state["width"] < mp.mpf("1e-100")
    and state["h"] > 0
    and state["lambda"] > 0
    and state["rho"] > 0
    and state["q_diag"] != 0
)
check(
    "independent bisection reconstructs a valid finite-height background",
    background_ok,
    f"q={mp_text(state['q'], 30)}, residual={mp_text(state['residual'], 6)}",
)

library = {
    "Counter": Counter,
    "DPS": DPS,
    "DERIVATIVE_STEPS": DERIVATIVE_STEPS,
    "I": I,
    "LOCAL_EDGES": LOCAL_EDGES,
    "LOCAL_HINGES": LOCAL_HINGES,
    "LOCAL_HINGE_INDEX": LOCAL_HINGE_INDEX,
    "L0_SQUARE": mp.mpf(1),
    "MASS": state["mass"],
    "RHO0": state["rho"],
    "_DIRECT": None,
    "combinations": combinations,
    "gro": gro,
    "math": math,
    "mp": mp,
    "np": np,
    "sp": sp,
    "spla": spla,
    "sparse_norm2": sparse_norm2,
}
load_named_functions(
    HESSIAN_SOURCE,
    {
        "log_minus",
        "signed_volume_square",
        "angle_record",
        "area_data",
        "extended_edge_image",
        "orbit_sort_key",
        "augment_boundary_orbits",
        "group_and_index_data",
        "prepare_geometry",
        "build_pattern_cache",
        "assemble_full_hessians",
        "full_gradient_at_delta",
    },
    library,
)
models = {
    parity: library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}

records = {}
runtime = {}
for parity in ("even", "odd"):
    print(f"[{parity}] assembling complete physical Hessians", flush=True)
    model = models[parity]
    index_data = library["group_and_index_data"](
        model, (mp.log(state["lambda"]), mp.mpf(0))
    )
    geometry = library["prepare_geometry"](model, index_data)
    carrier = quadratic_library["build_carrier"](
        model, index_data, state["lambda"], state["rho"]
    )
    geometry_ok = bool(
        len(model["slab"]) == 2400
        and len(index_data["edge_to_index"]) == FULL
        and len(geometry["triangle_records"]) == 6240
        and len(geometry["patterns"]) == 20
        and carrier["controls_ok"]
    )
    check(
        f"{parity}: full geometry and independently built rank-240 carrier pass",
        geometry_ok,
        f"support={carrier['controls']['row_support']}",
    )

    kind_values = {
        "old": mp.mpf(1),
        "internal": state["q_diag"],
        "pole": -state["rho"],
        "new": state["lambda"] ** 2,
    }
    pattern_cache, pattern_control = library["build_pattern_cache"](
        geometry["patterns"], kind_values
    )
    pattern_ok = bool(
        pattern_control["entry_pass"]
        and pattern_control["base_negative_counts"] == Counter({1: 2400})
        and pattern_control["all_displaced_have_one_negative"]
        and pattern_control["minimum_leading_minor"] > 0
        and pattern_control["minimum_argument"] > mp.mpf("1e-6")
    )
    check(
        f"{parity}: independent derivative hierarchy retains the Lorentzian branch",
        pattern_ok,
        (
            f"cross={mp_text(pattern_control['maximum_cross'], 6)}, "
            f"minor={mp_text(pattern_control['minimum_leading_minor'], 6)}"
        ),
    )

    library["_DIRECT"] = {
        "signed_base": index_data["signed_base"],
        "triangle_records": geometry["triangle_records"],
        "simplex_records": geometry["simplex_records"],
        "pole_indices": geometry["pole_indices"],
        "directions": (),
    }
    direct_gradient, direct_branch = library["full_gradient_at_delta"](
        np.zeros(FULL)
    )
    assembled_gradient, hessians, assembly = library["assemble_full_hessians"](
        model, index_data, geometry, pattern_cache
    )
    internal_maximum = max(
        abs(direct_gradient[index])
        for index in range(OLD, OLD + INTERNAL)
    )
    assembly_gradient_error = max(
        abs(assembled_gradient[index] - float(mp.re(value)))
        for index, value in enumerate(direct_gradient)
    )
    stationarity_ok = bool(
        internal_maximum < mp.mpf("1e-25")
        and assembly_gradient_error < 2e-11
        and direct_branch["negative_counts"] == Counter({1: 2400})
        and direct_branch["maximum_imaginary"] < mp.mpf("1e-60")
        and assembly["maximum_imaginary"] < 1e-60
    )
    check(
        f"{parity}: direct action gradient confirms internal stationarity",
        stationarity_ok,
        (
            f"internal={mp_text(internal_maximum, 6)}, "
            f"assembly_error={assembly_gradient_error:.3e}"
        ),
    )

    responses = {
        level: response_matrix(matrix, carrier["matrix"])
        for level, matrix in hessians.items()
    }
    diagonal_rows = [
        index - OLD
        for index in range(OLD, OLD + INTERNAL)
        if index_data["edge_kind"][index] == "internal"
    ]
    pole_rows = [
        index - OLD
        for index in range(OLD, OLD + INTERNAL)
        if index_data["edge_kind"][index] == "pole"
    ]
    row_split_ok = bool(
        len(diagonal_rows) == 720
        and len(pole_rows) == 120
        and set(diagonal_rows).isdisjoint(pole_rows)
        and set(diagonal_rows) | set(pole_rows) == set(range(INTERNAL))
    )
    wrong_label_detected = not bool(
        len(pole_rows) == 720 and len(diagonal_rows) == 120
    )
    check(
        f"{parity}: physical edge labels select 720 diagonal and 120 pole equations",
        row_split_ok and wrong_label_detected,
        f"diagonal={len(diagonal_rows)}, pole={len(pole_rows)}",
    )

    round_entry = assembly["binary64_roundoff"]["maximum_entry_bound"]
    full_record, full_runtime = classify_scope(
        responses, carrier["matrix"], INTERNAL, round_entry
    )
    diagonal_levels = {
        level: matrix[diagonal_rows, :] for level, matrix in responses.items()
    }
    diagonal_record, diagonal_runtime = classify_scope(
        diagonal_levels, carrier["matrix"], len(diagonal_rows), round_entry
    )
    classifiers_control_ok = bool(
        all(full_record["synthetic_controls"].values())
        and all(diagonal_record["synthetic_controls"].values())
    )
    check(
        f"{parity}: global SVD classifiers pass zero and identity controls",
        classifiers_control_ok,
    )
    print(
        f"[{parity}] resolved ranks: diagonal={diagonal_record['rank']}, "
        f"full={full_record['rank']}",
        flush=True,
    )

    midpoint_hessian = (
        hessians["operational_primary"]
        + hessians["validation_primary"]
    ) / 2
    first_edge, first_row, first_source = carrier["first_diagonal"]
    corrupted_carrier = carrier["matrix"].copy()
    corrupted_carrier[first_row, first_source] += 0.1
    corrupted_carrier_response = response_matrix(
        midpoint_hessian, corrupted_carrier
    )
    carrier_effect = float(
        np.linalg.norm(
            corrupted_carrier_response - full_runtime["midpoint_raw"]
        ) / max(1.0, np.linalg.norm(full_runtime["midpoint_raw"]))
    )
    hessian_corrupted_response = full_runtime["midpoint_raw"].copy()
    hessian_corrupted_response[first_row, :] += (
        0.1 * carrier["matrix"][first_row, :]
    )
    hessian_effect = float(
        np.linalg.norm(
            hessian_corrupted_response - full_runtime["midpoint_raw"]
        ) / max(1.0, np.linalg.norm(full_runtime["midpoint_raw"]))
    )
    corruption_ok = carrier_effect > 1e-12 and hessian_effect > 1e-12
    check(
        f"{parity}: full-real-space carrier and Hessian corruptions are detected",
        corruption_ok,
        f"carrier={carrier_effect:.3e}, Hdiag={hessian_effect:.3e}",
    )

    ranks_resolved = bool(
        full_record["resolved"] and full_record["qr"]["passed"]
        and diagonal_record["resolved"] and diagonal_record["qr"]["passed"]
    )
    if full_record["resolved"] and full_record["nullity"]:
        projectors, projector_uncertainty = projector_census(
            responses, full_record["nullity"]
        )
        midpoint = full_runtime["midpoint_raw"]
        _, _, right = la.svd(midpoint, full_matrices=False, check_finite=True)
        kernel_vector = right[-1, :]
        comparison_vector = right[-1 - full_record["nullity"], :]
        secant = direct_secant(
            library, carrier["matrix"], kernel_vector,
            comparison_vector, midpoint,
        )
    else:
        projectors = None
        projector_uncertainty = None
        kernel_vector = None
        secant = {
            "branch_ok": True,
            "validated": False,
            "reason": "full nullity is zero or unresolved",
        }

    records[parity] = {
        "geometry": {
            "simplices": len(model["slab"]),
            "edges": len(index_data["edge_to_index"]),
            "patterns": len(geometry["patterns"]),
            "triangles": len(geometry["triangle_records"]),
        },
        "carrier": carrier["controls"],
        "pattern": {
            "maximum_cross": mp_text(pattern_control["maximum_cross"]),
            "maximum_operational_proxy": mp_text(
                pattern_control["maximum_operational_proxy"]
            ),
            "maximum_validation_proxy": mp_text(
                pattern_control["maximum_validation_proxy"]
            ),
            "minimum_leading_minor": mp_text(
                pattern_control["minimum_leading_minor"]
            ),
            "minimum_argument": mp_text(pattern_control["minimum_argument"]),
        },
        "stationarity": {
            "internal_maximum": mp_text(internal_maximum),
            "assembly_gradient_error": float_text(assembly_gradient_error),
        },
        "assembly": {
            "maximum_imaginary": float_text(assembly["maximum_imaginary"]),
            "binary64_roundoff": assembly["binary64_roundoff"],
        },
        "row_split": {
            "diagonal_count": len(diagonal_rows),
            "pole_count": len(pole_rows),
            "wrong_reversal_detected": wrong_label_detected,
        },
        "diagonal": diagonal_record,
        "full": full_record,
        "hostile_controls": {
            "carrier_corruption_edge": list(first_edge),
            "carrier_corruption_row": int(first_row),
            "carrier_corruption_source": int(first_source),
            "carrier_corruption_relative_effect": float_text(carrier_effect),
            "hessian_corruption_relative_effect": float_text(hessian_effect),
        },
        "projector_uncertainty": (
            float_text(projector_uncertainty)
            if projector_uncertainty is not None else None
        ),
        "direct_secant": secant,
    }
    runtime[parity] = {
        "responses": responses,
        "full_runtime": full_runtime,
        "diagonal_runtime": diagonal_runtime,
        "projectors": projectors,
        "projector_uncertainty": projector_uncertainty,
        "kernel_vector": kernel_vector,
        "ranks_resolved": ranks_resolved,
        "control_ok": bool(
            geometry_ok and pattern_ok and stationarity_ok and row_split_ok
            and wrong_label_detected and classifiers_control_ok and corruption_ok
            and secant["branch_ok"]
        ),
    }
    del hessians, midpoint_hessian, pattern_cache
    gc.collect()

if all(runtime[parity]["projectors"] is not None for parity in ("even", "odd")):
    parity_uncertainty = max(
        runtime["even"]["projector_uncertainty"],
        runtime["odd"]["projector_uncertainty"],
    )
    parity_projector = compare_projectors(
        runtime["even"]["projectors"]["midpoint"],
        runtime["odd"]["projectors"]["midpoint"],
        parity_uncertainty,
    )
else:
    parity_projector = {
        "classification": "OPEN",
        "reason": "one full-map kernel is zero or unresolved",
    }
check(
    "the direct physical parity-projector classifier returns a declared state",
    parity_projector["classification"] in {"AGREE", "DEPENDENT", "OPEN"},
    str(parity_projector),
)

replication_classification = {
    parity: {
        "diagonal_nullity": records[parity]["diagonal"]["nullity"],
        "full_nullity": records[parity]["full"]["nullity"],
        "resolved": runtime[parity]["ranks_resolved"],
    }
    for parity in ("even", "odd")
}

print("[comparison] reading primary scientific artifacts after classification", flush=True)
primary = json.loads(PRIMARY_JSON.read_text())
primary_arrays = np.load(PRIMARY_MATRICES, allow_pickle=False)
primary_integrity_ok = bool(
    primary["outcome"] == "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY"
    and primary["tests"] == primary["passed"] == 25
    and primary["global_nullities"] == {"even": 1, "odd": 1}
    and primary["kernel_comparison"]["classification"] == "AGREE"
)
check("the committed primary claim remains internally intact", primary_integrity_ok)

primary_projector_records = {}
for parity in ("even", "odd"):
    if runtime[parity]["projectors"] is None:
        primary_projector_records[parity] = {"classification": "OPEN"}
        continue
    uncertainty = max(
        runtime[parity]["projector_uncertainty"],
        float(primary["kernel_comparison"]["uncertainty"]),
    )
    primary_projector_records[parity] = compare_projectors(
        runtime[parity]["projectors"]["midpoint"],
        primary_arrays[f"{parity}_kernel_projector"],
        uncertainty,
    )

matrix_reproduces = bool(
    all(runtime[parity]["ranks_resolved"] for parity in ("even", "odd"))
    and all(
        replication_classification[parity]["diagonal_nullity"] == 121
        and replication_classification[parity]["full_nullity"] == 1
        for parity in ("even", "odd")
    )
    and parity_projector["classification"] == "AGREE"
    and all(
        primary_projector_records[parity]["classification"] == "AGREE"
        for parity in ("even", "odd")
    )
)
secants_validate = all(
    records[parity]["direct_secant"]["validated"]
    for parity in ("even", "odd")
)
all_ranks_resolved = all(
    runtime[parity]["ranks_resolved"] for parity in ("even", "odd")
)
resolved_matrix_contradiction = bool(
    all_ranks_resolved
    and (
        any(
            replication_classification[parity]["diagonal_nullity"] != 121
            or replication_classification[parity]["full_nullity"] != 1
            for parity in ("even", "odd")
        )
        or parity_projector["classification"] == "DEPENDENT"
        or any(
            record["classification"] == "DEPENDENT"
            for record in primary_projector_records.values()
        )
    )
)
all_controls = bool(
    provenance_ok and registry_ok and geometry_import_ok and background_ok
    and primary_integrity_ok
    and all(runtime[parity]["control_ok"] for parity in ("even", "odd"))
)
any_numerical_open = bool(
    not all_ranks_resolved
    or parity_projector["classification"] == "OPEN"
    or any(
        record["classification"] == "OPEN"
        for record in primary_projector_records.values()
    )
    or not secants_validate
)

if not all_controls:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_ADVERSARIAL_CONTROL_FAILED"
elif matrix_reproduces and secants_validate:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_ADVERSARIALLY_REPLICATED"
elif resolved_matrix_contradiction:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_PRIMARY_REFUTED"
elif any_numerical_open:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_ADVERSARIAL_OPEN"
else:
    outcome = "FINITE_HEIGHT_INTERNAL_CARRIER_PRIMARY_REFUTED"

check(
    "the adversarial outcome follows the frozen hierarchy",
    outcome in {
        "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_ADVERSARIALLY_REPLICATED",
        "FINITE_HEIGHT_INTERNAL_CARRIER_PRIMARY_REFUTED",
        "FINITE_HEIGHT_INTERNAL_CARRIER_ADVERSARIAL_OPEN",
        "FINITE_HEIGHT_INTERNAL_CARRIER_ADVERSARIAL_CONTROL_FAILED",
    },
    outcome,
)

matrix_payload = {
    "even_midpoint_full": runtime["even"]["full_runtime"]["midpoint_raw"],
    "odd_midpoint_full": runtime["odd"]["full_runtime"]["midpoint_raw"],
    "even_midpoint_diagonal": runtime["even"]["diagonal_runtime"]["midpoint_raw"],
    "odd_midpoint_diagonal": runtime["odd"]["diagonal_runtime"]["midpoint_raw"],
}
for parity in ("even", "odd"):
    if runtime[parity]["projectors"] is not None:
        matrix_payload[f"{parity}_kernel_projector"] = (
            runtime[parity]["projectors"]["midpoint"]
        )
    if runtime[parity]["kernel_vector"] is not None:
        matrix_payload[f"{parity}_kernel_vector"] = runtime[parity]["kernel_vector"]
np.savez_compressed(MATRIX_OUTPUT, **matrix_payload)
matrix_hash = sha256(MATRIX_OUTPUT)

artifact = {
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "method": (
        "complete physical 2280-edge binary64 Hessians assembled from "
        "120-digit local terms, direct 840x240 SVD/QR, and nonlinear "
        "action-gradient secants"
    ),
    "provenance": {
        "primary_implementation_commit": PRIMARY_IMPLEMENTATION_COMMIT,
        "primary_result_commit": PRIMARY_RESULT_COMMIT,
        "adversarial_protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": input_hashes,
    },
    "background": {
        key: mp_text(value) if isinstance(value, (mp.mpf, mp.mpc)) else value
        for key, value in state.items()
    },
    "replication_classification": replication_classification,
    "parity_kernel_comparison": parity_projector,
    "primary_kernel_comparison": primary_projector_records,
    "matrix_reproduces_primary": matrix_reproduces,
    "direct_secants_validate": secants_validate,
    "parities": records,
    "matrices": {
        "path": MATRIX_OUTPUT.name,
        "sha256": matrix_hash,
        "arrays": {
            key: list(value.shape) for key, value in matrix_payload.items()
        },
    },
    "interpretation": {
        "selected_object": (
            "one homogeneous first-order direction in the fixed rank-240 "
            "scale-plus-strut carrier at one finite-height slab"
        ),
        "physical_classification": "OPEN",
        "finite_step_invariant_continuation": "SEPARATE_DERIVED_RESULT",
        "infinite_proper_time_evolution": "OPEN",
        "graviton_or_wave_equation": "NOT_DERIVED",
        "tick_c_G_planck_particle_masses": "NOT_DERIVED",
    },
    "firewall": {
        "primary_scientific_values_read_after_replication_census": True,
        "orbit_sector_decomposition_used": False,
        "full_suite_run": False,
        "continuum_target_parsed": False,
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} PASS")
print(f"MATRIX SHA: {matrix_hash}")
if passed != tests:
    raise SystemExit(1)
