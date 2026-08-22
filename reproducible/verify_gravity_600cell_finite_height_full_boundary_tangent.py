#!/usr/bin/env python3
"""Target-free full canonical tangent at the first finite-height dust slab.

Prior-art gate commit: c3fc22d.
Protocol commit: 2f6a4a5.
Registry commit: aae40c7.

No eigenvalue, continuum harmonic, speed, Planck or particle target is used.
"""

import ast
from collections import Counter, defaultdict
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import zipfile

from flint import acb, acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_finite_height_full_boundary_tangent.json"
NUMERIC_OUTPUT = (
    HERE / "gravity_600cell_finite_height_full_boundary_tangent.npz"
)
PRIOR_ART = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_prior_art.md"
)
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_carrier_quadratic.json"
ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_adversarial.json"
)
INTERNAL_RANK_INPUT = (
    HERE / "gravity_600cell_finite_height_internal_carrier_rank.json"
)
RECONCILIATION_INPUT = (
    HERE / "gravity_600cell_finite_height_internal_kernel_canonical_reconciliation.json"
)
OLD_TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
OLD_TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
OLD_TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "prior_art": PRIOR_ART,
    "primary": PRIMARY_INPUT,
    "adversarial": ADVERSARIAL_INPUT,
    "internal_rank": INTERNAL_RANK_INPUT,
    "reconciliation": RECONCILIATION_INPUT,
    "old_tangent": OLD_TANGENT_INPUT,
    "old_tangent_numeric": OLD_TANGENT_NUMERIC,
    "old_tangent_source": OLD_TANGENT_SOURCE,
    "rank_source": RANK_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
}
EXPECTED_HASHES = {
    "prior_art": "6fe3e10daf97fd60849a837e56716ced594e19c77117ecc14f862822edc10074",
    "primary": "0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9",
    "adversarial": "54915cf364c36af6bbc8e1dbd36433079269d293453478bfdf589e547d462ad6",
    "internal_rank": "513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8",
    "reconciliation": "81ec0379247023451e82ab42f5beb026ee2d1b083aa5e2553e42b894554266f6",
    "old_tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "old_tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "old_tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

PRIOR_ART_COMMIT = "c3fc22d"
PROTOCOL_COMMIT = "2f6a4a5"
REGISTRY_COMMIT = "aae40c7"
VERIFIER_NAME = Path(__file__).name

DPS = 180
BALL_DPS = 140
mp.mp.dps = DPS
ctx.dps = BALL_DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-25"),
    "operational_shadow": mp.mpf("5e-26"),
    "validation_primary": mp.mpf("2.5e-26"),
    "validation_shadow": mp.mpf("1.25e-26"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
LEVELS = ("K01", "K12", "K23")
ARITHMETIC_FLOOR = mp.mpf("1e-150")
CLASSIFIER_FLOOR = mp.mpf("1e-135")
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
ORBIT_ORDER = (30, 35, 30)

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


def complex_record(value):
    return {
        "real": float_text(np.real(value)),
        "imaginary": float_text(np.imag(value)),
    }


def load_named_functions(path, names, namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    found = {node.name for node in selected}
    missing = set(names) - found
    if missing:
        raise RuntimeError(f"missing audited functions: {sorted(missing)}")
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


def finite_height_state(committed):
    v = mp.mpf(3) / 2

    def epsilon(value):
        square = value**2
        return 2 * mp.pi - 5 * mp.acos(
            (square + 2) / (2 * (square + 3))
        )

    def mu(value):
        return 180 * epsilon(value) / (mp.pi * mp.sqrt(value**2 + 4))

    def momentum(value):
        square = value**2
        return (
            180 * value * epsilon(value) / mp.sqrt(square + 4)
            - 600 * mp.sqrt(3)
            * mp.asinh(value / mp.sqrt(8 * (square + 3)))
        )

    mass = mu(v)
    incoming = momentum(v)

    def elimination(value):
        return (
            4 * mp.pi * (mu(value) - mass)
            + value * (momentum(value) - incoming)
        )

    left = mp.mpf(9)
    right = mp.mpf(10)
    left_value = elimination(left)
    right_value = elimination(right)
    bracket_ok = left_value * right_value < 0
    for _ in range(600):
        middle = (left + right) / 2
        middle_value = elimination(middle)
        if middle_value == 0:
            left = right = middle
            break
        if left_value * middle_value < 0:
            right = middle
            right_value = middle_value
        else:
            left = middle
            left_value = middle_value
        if right - left < mp.mpf("1e-160"):
            break
    q = (left + right) / 2
    h = (momentum(q) - incoming) / (2 * mp.pi * mu(q))
    lam = 1 + h * q
    rho = h**2
    committed_error = max(
        abs(q - mp.mpf(committed["q"])),
        abs(h - mp.mpf(committed["h"])),
        abs(lam - mp.mpf(committed["lambda"])),
        abs(rho - mp.mpf(committed["rho"])),
        abs(mass - mp.mpf(committed["mass"])),
    )
    return {
        "v": v,
        "q": q,
        "h": h,
        "lambda": lam,
        "rho": rho,
        "mass": mass,
        "elimination": elimination(q),
        "bracket_width": right - left,
        "bracket_ok": bracket_ok,
        "committed_error": committed_error,
    }


def mp_frobenius(matrix):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in matrix))


def mp_difference_frobenius(left, right):
    if left.rows != right.rows or left.cols != right.cols:
        raise ValueError("matrix shape mismatch")
    return mp.sqrt(mp.fsum(
        abs(left[row, column] - right[row, column]) ** 2
        for row in range(left.rows) for column in range(left.cols)
    ))


def mp_submatrix(matrix, rows, columns):
    return mp.matrix([
        [matrix[row, column] for column in columns]
        for row in rows
    ])


def mp_to_numpy(matrix):
    return np.asarray([
        [
            complex(
                float(mp.re(matrix[row, column])),
                float(mp.im(matrix[row, column])),
            )
            for column in range(matrix.cols)
        ]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def richardson_kernel(coarse, fine):
    keys = set(coarse) | set(fine)
    return {
        key: (4 * fine.get(key, 0) - coarse.get(key, 0)) / 3
        for key in keys
    }


def mp_to_acb(value):
    return acb(
        mp.nstr(mp.re(value), 170),
        mp.nstr(mp.im(value), 170),
    )


def mp_matrix_to_acb(matrix):
    return acb_mat(
        matrix.rows,
        matrix.cols,
        [mp_to_acb(value) for value in matrix],
    )


def acb_midpoint_and_radii(matrix):
    midpoint = np.empty((matrix.nrows(), matrix.ncols()), dtype=np.complex128)
    radii = np.empty((matrix.nrows(), matrix.ncols()), dtype=float)
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            value = matrix[row, column]
            midpoint[row, column] = complex(
                float(value.real.mid()),
                float(value.imag.mid()),
            )
            radii[row, column] = math.hypot(
                float(value.real.rad().upper()),
                float(value.imag.rad().upper()),
            )
    return midpoint, radii


def deterministic_npz(path, arrays):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with zipfile.ZipFile(
        temporary, "w", compression=zipfile.ZIP_DEFLATED
    ) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(
                buffer, np.asarray(arrays[name]), allow_pickle=False
            )
            info = zipfile.ZipInfo(
                f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0)
            )
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, buffer.getvalue())
    temporary.replace(path)


def expanded_types(start, stop, dimension):
    return [
        kind * dimension + component
        for kind in range(start, stop)
        for component in range(dimension)
    ]


def pre_legendre_matrix(block, dimension):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    xd = len(internal)
    nd = len(new)
    result = mp.matrix(xd + nd, xd + nd)
    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    for row in range(xd):
        for column in range(xd):
            result[row, column] = k_xx[row, column]
        for column in range(nd):
            result[row, xd + column] = k_xn[row, column]
    for row in range(nd):
        for column in range(xd):
            result[xd + row, column] = -k_ox[row, column]
        for column in range(nd):
            result[xd + row, xd + column] = -k_on[row, column]
    return result


def canonical_tangent_ball(
    block,
    old,
    internal,
    new,
    output_final_indices,
    common_indices=None,
    pre_momentum_sign=1,
    include_k_no=True,
):
    """Apply the frozen implicit Legendre formula to one Hessian block."""
    od = len(old)
    xd = len(internal)
    nd = len(new)
    if od != nd or len(output_final_indices) != od:
        raise ValueError("boundary dimensions or output map are inconsistent")

    j_matrix = mp.matrix(xd + nd, xd + nd)
    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    for row in range(xd):
        for column in range(xd):
            j_matrix[row, column] = k_xx[row, column]
        for column in range(nd):
            j_matrix[row, xd + column] = k_xn[row, column]
    for row in range(nd):
        for column in range(xd):
            j_matrix[xd + row, column] = -k_ox[row, column]
        for column in range(nd):
            j_matrix[xd + row, xd + column] = -k_on[row, column]

    rhs = mp.matrix(xd + nd, 2 * od)
    k_xo = mp_submatrix(block, internal, old)
    k_oo = mp_submatrix(block, old, old)
    for row in range(xd):
        for column in range(od):
            rhs[row, column] = -k_xo[row, column]
    for row in range(nd):
        for column in range(od):
            rhs[xd + row, column] = k_oo[row, column]
        rhs[xd + row, od + row] = pre_momentum_sign

    j_ball = mp_matrix_to_acb(j_matrix)
    rhs_ball = mp_matrix_to_acb(rhs)
    det_j = j_ball.det()
    if det_j.contains(0):
        return {
            "j_matrix": j_matrix,
            "det_j": det_j,
            "tangent": None,
            "defect": None,
        }
    solved = j_ball.solve(rhs_ball)
    y_x = acb_mat(xd, 2 * od, [
        solved[row, column]
        for row in range(xd)
        for column in range(2 * od)
    ])
    y_n = acb_mat(nd, 2 * od, [
        solved[xd + row, column]
        for row in range(nd)
        for column in range(2 * od)
    ])

    k_nx = mp_matrix_to_acb(mp_submatrix(block, new, internal))
    k_nn = mp_matrix_to_acb(mp_submatrix(block, new, new))
    direct = acb_mat(nd, 2 * od)
    if include_k_no:
        k_no = mp_matrix_to_acb(mp_submatrix(block, new, old))
        for row in range(nd):
            for column in range(od):
                direct[row, column] = k_no[row, column]
    p_post = direct + k_nx * y_x + k_nn * y_n
    raw = acb_mat(2 * nd, 2 * od)
    for row in range(nd):
        for column in range(2 * od):
            raw[row, column] = y_n[row, column]
            raw[nd + row, column] = p_post[row, column]

    local = acb_mat(2 * od, 2 * od)
    for old_position, final_position in enumerate(output_final_indices):
        for column in range(2 * od):
            local[old_position, column] = raw[final_position, column]
            local[od + old_position, column] = raw[nd + final_position, column]

    if common_indices is None:
        common_indices = tuple(range(od))
    phase_indices = tuple(common_indices) + tuple(
        od + index for index in common_indices
    )
    tangent = acb_mat(2 * od, 2 * od)
    for row, local_row in enumerate(phase_indices):
        for column, local_column in enumerate(phase_indices):
            tangent[row, column] = local[local_row, local_column]

    omega = acb_mat(2 * od, 2 * od)
    for index in range(od):
        omega[index, od + index] = 1
        omega[od + index, index] = -1
    defect = tangent.transpose().conjugate() * omega * tangent - omega
    return {
        "j_matrix": j_matrix,
        "det_j": det_j,
        "tangent": tangent,
        "defect": defect,
    }


def scalar_formula_control():
    block = mp.matrix([
        [7, 3, 11],
        [3, 2, 5],
        [11, 5, 13],
    ])
    expected = np.asarray([
        [-5 / 7, -2 / 7],
        [22 / 7, -1 / 7],
    ], dtype=float)
    good = canonical_tangent_ball(
        block, [0], [1], [2], [0]
    )
    bad_direct = canonical_tangent_ball(
        block, [0], [1], [2], [0], include_k_no=False
    )
    bad_sign = canonical_tangent_ball(
        block, [0], [1], [2], [0], pre_momentum_sign=-1
    )
    good_mid, good_rad = acb_midpoint_and_radii(good["tangent"])
    good_defect, good_defect_rad = acb_midpoint_and_radii(good["defect"])
    bad_direct_mid, _ = acb_midpoint_and_radii(bad_direct["tangent"])
    bad_direct_defect, _ = acb_midpoint_and_radii(bad_direct["defect"])
    bad_sign_mid, _ = acb_midpoint_and_radii(bad_sign["tangent"])
    expected_bad = np.asarray([
        [-5 / 7, -2 / 7],
        [-55 / 7, -1 / 7],
    ], dtype=float)
    return {
        "passed": bool(
            not good["det_j"].contains(0)
            and abs(float(good["det_j"].real.mid()) + 7) < 1e-14
            and np.max(np.abs(good_mid.real - expected)) < 1e-14
            and np.max(np.abs(good_mid.imag)) < 1e-14
            and np.max(good_rad) < 1e-100
            and np.linalg.norm(good_defect) < 1e-14
            and np.max(good_defect_rad) < 1e-100
            and np.max(np.abs(bad_direct_mid.real - expected_bad)) < 1e-14
            and abs(bad_direct_defect[0, 1] + 22 / 7) < 1e-14
            and np.linalg.norm(bad_sign_mid - good_mid) > 1
        ),
        "determinant": str(good["det_j"]),
        "maximum_good_matrix_error": float(
            np.max(np.abs(good_mid.real - expected))
        ),
        "good_defect_norm": float(np.linalg.norm(good_defect)),
        "bad_direct_defect_norm": float(np.linalg.norm(bad_direct_defect)),
        "bad_sign_distance": float(np.linalg.norm(bad_sign_mid - good_mid)),
    }


def boundary_identification(index_data):
    old_seeds = tuple(
        min(index_data["orbit_edges"][orbit_type])
        for orbit_type in range(30)
    )
    mapping = []
    group_preserved = True
    for old_type in range(30):
        candidates = set()
        for group, edge in enumerate(index_data["orbit_edges"][old_type]):
            shifted = tuple(int(vertex) + VERTICES for vertex in edge)
            global_index = index_data["edge_to_index"].get(shifted)
            if global_index is None:
                group_preserved = False
                continue
            final_type, final_group = divmod(global_index, 24)
            candidates.add(final_type - 65)
            group_preserved &= final_group == group
        if len(candidates) != 1:
            mapping.append(-1)
        else:
            mapping.append(candidates.pop())
    mapping = tuple(mapping)
    ok = bool(
        sorted(mapping) == list(range(30))
        and group_preserved
        and len(set(old_seeds)) == 30
    )
    return {
        "mapping": mapping,
        "old_seeds": old_seeds,
        "group_preserved": group_preserved,
        "ok": ok,
    }


def expanded_common_indices(local_old_seeds, common_seeds, dimension):
    by_seed = {seed: index for index, seed in enumerate(local_old_seeds)}
    if set(by_seed) != set(common_seeds):
        raise RuntimeError("old boundary orbit seed sets differ")
    return tuple(
        by_seed[seed] * dimension + component
        for seed in common_seeds
        for component in range(dimension)
    )


def expanded_output_mapping(mapping, dimension):
    return tuple(
        mapping[old_type] * dimension + component
        for old_type in range(30)
        for component in range(dimension)
    )


def scaled_frobenius(matrix, scales):
    return mp.sqrt(mp.fsum(
        abs(matrix[row, column] * scales[column]) ** 2
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    ))


def scaled_difference_frobenius(left, right, scales):
    return mp.sqrt(mp.fsum(
        abs((left[row, column] - right[row, column]) * scales[column]) ** 2
        for row in range(left.rows)
        for column in range(left.cols)
    ))


def pre_legendre_analysis(blocks, dimension):
    matrices = {
        level: pre_legendre_matrix(block, dimension)
        for level, block in blocks.items()
    }
    middle_numpy = mp_to_numpy(matrices["K12"])
    column_norms = np.linalg.norm(middle_numpy, axis=0)
    scales = np.asarray(
        [1 / max(1.0, float(value)) for value in column_norms],
        dtype=float,
    )
    scale_mp = tuple(mp.mpf(float(value)) for value in scales)
    normalization = max(
        mp.mpf(1), scaled_frobenius(matrices["K12"], scale_mp)
    )
    step_error = max(
        scaled_difference_frobenius(
            matrices["K01"], matrices["K12"], scale_mp
        ),
        scaled_difference_frobenius(
            matrices["K12"], matrices["K23"], scale_mp
        ),
    ) / normalization

    singular = {}
    conditions = {}
    for level, matrix in matrices.items():
        scaled = mp_to_numpy(matrix) * scales[np.newaxis, :]
        values = la.svd(
            scaled, compute_uv=False, lapack_driver="gesvd"
        )
        singular[level] = values
        conditions[level] = float(values[0] / values[-1])
    middle_gesdd = la.svd(
        middle_numpy * scales[np.newaxis, :],
        compute_uv=False,
        lapack_driver="gesdd",
    )
    driver_difference = float(np.max(np.abs(
        singular["K12"] - middle_gesdd
    )))
    middle_norm2 = float(singular["K12"][0])
    normalization_float = float(normalization)
    svd_error = max(
        driver_difference / normalization_float,
        (
            10 * np.finfo(float).eps * max(1.0, middle_norm2)
            / normalization_float
        ),
    )
    total_error = float(step_error) + svd_error + 1e-135
    normalized_minima = {
        level: float(values[-1] / normalization_float)
        for level, values in singular.items()
    }

    determinant_balls = {}
    determinant_flags = {}
    for level, matrix in matrices.items():
        determinant = mp_matrix_to_acb(matrix).det()
        determinant_balls[level] = str(determinant)
        determinant_flags[level] = not determinant.contains(0)
    gap_ok = all(
        value > 100 * total_error
        for value in normalized_minima.values()
    )
    regular = bool(gap_ok and all(determinant_flags.values()))
    return {
        "matrices": matrices,
        "regular": regular,
        "record": {
            "classification": "REGULAR" if regular else "NUMERICALLY_OPEN",
            "shape": list(middle_numpy.shape),
            "column_scale_minimum": float_text(np.min(scales)),
            "column_scale_maximum": float_text(np.max(scales)),
            "normalization": mp_text(normalization),
            "step_error": mp_text(step_error),
            "svd_driver_difference": float_text(driver_difference),
            "svd_error": float_text(svd_error),
            "total_error": float_text(total_error),
            "normalized_smallest_singular_values": {
                key: float_text(value)
                for key, value in normalized_minima.items()
            },
            "condition_estimates": {
                key: float_text(value) for key, value in conditions.items()
            },
            "singular_values": {
                key: [float_text(value) for value in values]
                for key, values in singular.items()
            },
            "determinants_exclude_zero": determinant_flags,
            "determinant_balls": determinant_balls,
            "gap_gate": float_text(100 * total_error),
        },
    }


def tangent_analysis(ball_records):
    midpoints = {
        level: item["midpoint"] for level, item in ball_records.items()
    }
    radii = {level: item["radii"] for level, item in ball_records.items()}
    defects = {
        level: item["defect_midpoint"]
        for level, item in ball_records.items()
    }
    defect_radii = {
        level: item["defect_radii"]
        for level, item in ball_records.items()
    }
    middle = midpoints["K12"]
    tangent_variation = (
        la.norm(midpoints["K01"] - middle, "fro")
        + la.norm(middle - midpoints["K23"], "fro")
    )
    maximum_tangent_radius = max(
        la.norm(value, "fro") for value in radii.values()
    )
    defect_variation = (
        la.norm(defects["K01"] - defects["K12"], "fro")
        + la.norm(defects["K12"] - defects["K23"], "fro")
    )
    maximum_defect_radius = max(
        la.norm(value, "fro") for value in defect_radii.values()
    )
    middle_norm2 = float(la.norm(middle, 2))
    tangent_effect = (
        2 * max(1.0, middle_norm2) * tangent_variation
        + tangent_variation**2
    )
    round_effect = (
        100 * np.finfo(float).eps * max(1.0, middle_norm2) ** 2
    )
    epsilon_symplectic = (
        tangent_effect
        + defect_variation
        + maximum_tangent_radius
        + maximum_defect_radius
        + round_effect
        + 1e-135
    )
    symplectic_norm = float(la.norm(defects["K12"], 2))
    symplectic_ok = bool(symplectic_norm <= 10 * epsilon_symplectic)

    log_moduli = {}
    for level, matrix in midpoints.items():
        _, log_modulus = np.linalg.slogdet(matrix)
        log_moduli[level] = float(log_modulus)
    condition = float(np.linalg.cond(middle))
    radius_log_effect = (
        middle.shape[0]
        * max(1.0, condition)
        * maximum_tangent_radius
        / max(1.0, middle_norm2)
    )
    epsilon_log_det = (
        abs(log_moduli["K01"] - log_moduli["K12"])
        + abs(log_moduli["K12"] - log_moduli["K23"])
        + radius_log_effect
        + 100 * np.finfo(float).eps * middle.shape[0] * max(1.0, condition)
        + 1e-135
    )
    determinant_ok = bool(
        abs(log_moduli["K12"]) <= 10 * epsilon_log_det
    )
    return {
        "midpoints": midpoints,
        "radii": radii,
        "defects": defects,
        "defect_radii": defect_radii,
        "canonicality_ok": bool(symplectic_ok and determinant_ok),
        "record": {
            "tangent_variation_frobenius": float_text(tangent_variation),
            "maximum_tangent_ball_radius_frobenius": float_text(
                maximum_tangent_radius
            ),
            "defect_variation_frobenius": float_text(defect_variation),
            "maximum_defect_ball_radius_frobenius": float_text(
                maximum_defect_radius
            ),
            "symplectic_defect_two_norm": float_text(symplectic_norm),
            "epsilon_symplectic": float_text(epsilon_symplectic),
            "symplectic_ok": symplectic_ok,
            "log_determinant_moduli": {
                key: float_text(value) for key, value in log_moduli.items()
            },
            "epsilon_log_determinant_modulus": float_text(epsilon_log_det),
            "determinant_modulus_ok": determinant_ok,
            "tangent_condition_estimate": float_text(condition),
            "canonicality_ok": bool(symplectic_ok and determinant_ok),
        },
    }


def schedule_compare(even, odd):
    middle_even = even["midpoints"]["K12"]
    middle_odd = odd["midpoints"]["K12"]
    normalization = max(
        1.0,
        float(la.norm(middle_even, "fro")),
        float(la.norm(middle_odd, "fro")),
    )

    def within(analysis):
        middle = analysis["midpoints"]["K12"]
        return (
            la.norm(analysis["midpoints"]["K01"] - middle, "fro")
            + la.norm(middle - analysis["midpoints"]["K23"], "fro")
            + max(la.norm(value, "fro") for value in analysis["radii"].values())
            + 100 * np.finfo(float).eps * max(1.0, la.norm(middle, 2))
        )

    uncertainty = (
        float(within(even) + within(odd)) / normalization + 1e-135
    )
    distance = float(
        la.norm(middle_even - middle_odd, "fro") / normalization
    )
    if distance <= 10 * uncertainty:
        label = "SCHEDULE_ROBUST"
    elif distance > 100 * uncertainty:
        label = "SCHEDULE_DEPENDENT"
    else:
        label = "SCHEDULE_OPEN"

    synthetic = middle_even.copy()
    synthetic[0, 0] += 1e-3
    synthetic_distance = float(
        la.norm(synthetic - middle_even, "fro") / normalization
    )
    synthetic_detected = bool(synthetic_distance > 100 * uncertainty)

    nd = middle_even.shape[0] // 2
    dimension = nd // 30
    cyclic_configuration = np.asarray([
        ((orbit + 1) % 30) * dimension + component
        for orbit in range(30)
        for component in range(dimension)
    ], dtype=int)
    cyclic_phase = np.concatenate((
        cyclic_configuration,
        nd + cyclic_configuration,
    ))
    cyclic = middle_even[cyclic_phase, :]
    cyclic_distance = float(
        la.norm(cyclic - middle_even, "fro") / normalization
    )
    cyclic_detected = bool(cyclic_distance > 100 * uncertainty)
    return {
        "label": label,
        "distance": distance,
        "uncertainty": uncertainty,
        "synthetic_distance": synthetic_distance,
        "synthetic_detected": synthetic_detected,
        "cyclic_output_map_distance": cyclic_distance,
        "cyclic_output_map_detected": cyclic_detected,
    }


print("=" * 78)
print("FINITE-HEIGHT FULL 1440-DIMENSIONAL CANONICAL TANGENT")
print("=" * 78)

input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
payloads = {
    "primary": json.loads(PRIMARY_INPUT.read_text()),
    "adversarial": json.loads(ADVERSARIAL_INPUT.read_text()),
    "internal_rank": json.loads(INTERNAL_RANK_INPUT.read_text()),
    "reconciliation": json.loads(RECONCILIATION_INPUT.read_text()),
    "old_tangent": json.loads(OLD_TANGENT_INPUT.read_text()),
}
provenance_ok = input_hashes == EXPECTED_HASHES
check(
    "all preregistered inputs and implementation sources retain frozen hashes",
    provenance_ok,
    str(input_hashes),
)

registered_scripts, registry_duplicates = registry_inventory(RUN_ALL)
registry_ok = bool(
    registered_scripts.count(VERIFIER_NAME) == 1
    and not registry_duplicates
)
check(
    "the verifier is registered exactly once and the registry has no duplicates",
    registry_ok,
    f"entries={len(registered_scripts)}, duplicates={registry_duplicates}",
)

accepted_inputs_ok = bool(
    payloads["primary"]["outcome"]
    == "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY"
    and payloads["primary"]["passed"] == payloads["primary"]["tests"] == 22
    and payloads["adversarial"]["outcome"]
    == "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_ADVERSARIALLY_REPLICATED"
    and payloads["adversarial"]["passed"]
    == payloads["adversarial"]["tests"] == 18
    and payloads["internal_rank"]["outcome"]
    == "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY"
    and payloads["internal_rank"]["passed"]
    == payloads["internal_rank"]["tests"] == 25
    and payloads["reconciliation"]["outcome"]
    == "INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_FIXED_INPUT_REMOVES_IT"
    and payloads["reconciliation"]["passed"]
    == payloads["reconciliation"]["tests"] == 14
    and payloads["old_tangent"]["outcome"]
    == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and payloads["old_tangent"]["passed"]
    == payloads["old_tangent"]["tests"] == 19
    and payloads["old_tangent"]["numeric_archive_sha256"]
    == EXPECTED_HASHES["old_tangent_numeric"]
)
check(
    "the accepted finite-height controls and old formula control remain intact",
    accepted_inputs_ok,
)

state = finite_height_state(payloads["primary"]["background"])
background_ok = bool(
    state["bracket_ok"]
    and abs(state["elimination"]) < mp.mpf("1e-140")
    and state["bracket_width"] < mp.mpf("1e-150")
    and state["committed_error"] < mp.mpf("1e-70")
    and state["h"] > 0
    and state["lambda"] > 0
    and state["rho"] > 0
    and state["lambda"] - state["rho"] > 0
)
check(
    "deterministic bisection independently reconstructs the finite-height background",
    background_ok,
    (
        f"q={mp_text(state['q'], 35)}, "
        f"E={mp_text(state['elimination'], 6)}, "
        f"width={mp_text(state['bracket_width'], 6)}"
    ),
)

scalar_control = scalar_formula_control()
check(
    "the exact scalar generating-function control passes and both formula corruptions fail",
    scalar_control["passed"],
    str(scalar_control),
)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_finite_height_full_tangent", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
geometry_import_ok = bool(gro.tests == gro.passed == 43)
check(
    "the direct complete one-slab geometry import retains all 43 certificates",
    geometry_import_ok,
)

library = {
    "ARITHMETIC_FLOOR": ARITHMETIC_FLOOR,
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
    "VARIANTS": VARIANTS,
    "cluster_sorted": None,
    "combinations": combinations,
    "defaultdict": defaultdict,
    "gro": gro,
    "math": math,
    "mp": mp,
    "mp_frobenius": mp_frobenius,
    "mp_submatrix": mp_submatrix,
    "np": np,
}
load_named_functions(
    RANK_SOURCE,
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
    },
    library,
)
load_named_functions(
    OLD_TANGENT_SOURCE,
    {
        "cluster_sorted",
        "high_precision_sector_bases",
        "high_precision_pattern_cache",
        "assemble_full_representative_kernels",
        "project_full_kernel",
    },
    library,
)
models = {
    parity: library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}

# Freeze both common representation carriers before constructing a Hessian.
prepared = {}
for parity in ("even", "odd"):
    model = models[parity]
    index_data = library["group_and_index_data"](
        model, (mp.log(state["lambda"]), mp.mpf(0))
    )
    geometry = library["prepare_geometry"](model, index_data)
    boundary = boundary_identification(index_data)
    sectors, sector_control = library["high_precision_sector_bases"](
        index_data
    )
    carrier_ok = bool(
        len(model["old_edges"]) == OLD
        and len(model["internal_edges"]) == INTERNAL
        and len(model["new_edges"]) == NEW
        and len(model["slab"]) == 2400
        and len(geometry["triangle_records"]) == 6240
        and len(geometry["simplex_records"]) == 2400
        and len(geometry["patterns"]) == 20
        and len(index_data["edge_to_index"]) == FULL
        and len(index_data["orbit_edges"]) == sum(ORBIT_ORDER)
        and boundary["ok"]
    )
    basis_ok = bool(
        sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sum(sector["dimension"] ** 2 for sector in sectors) == 24
        and sum(60 * sector["dimension"] ** 2 for sector in sectors) == 1440
        and all(
            sector_control[key] < mp.mpf("1e-140")
            for key in (
                "maximum_orthonormal",
                "maximum_center",
                "maximum_splitter",
                "maximum_right_leakage",
                "maximum_conjugate_pair",
            )
        )
    )
    check(
        f"{parity}: complete carrier and literal 720-edge boundary shift pass",
        carrier_ok,
        f"mapping={boundary['mapping']}",
    )
    check(
        f"{parity}: seven target-free minimal sectors exhaust the phase carrier",
        basis_ok,
        f"dimensions={[sector['dimension'] for sector in sectors]}",
    )
    prepared[parity] = {
        "model": model,
        "index_data": index_data,
        "geometry": geometry,
        "boundary": boundary,
        "sectors": sectors,
        "sector_control": sector_control,
        "carrier_ok": carrier_ok,
        "basis_ok": basis_ok,
    }

table_equal = bool(np.array_equal(
    prepared["even"]["index_data"]["table"],
    prepared["odd"]["index_data"]["table"],
))
actions_equal = bool(
    len(prepared["even"]["index_data"]["actions"])
    == len(prepared["odd"]["index_data"]["actions"]) == 24
    and all(
        np.array_equal(left, right)
        for left, right in zip(
            prepared["even"]["index_data"]["actions"],
            prepared["odd"]["index_data"]["actions"],
        )
    )
)
seed_sets_equal = bool(
    set(prepared["even"]["boundary"]["old_seeds"])
    == set(prepared["odd"]["boundary"]["old_seeds"])
)
basis_distance = mp.mpf(0)
basis_signature_equal = True
for even_sector, odd_sector in zip(
    prepared["even"]["sectors"], prepared["odd"]["sectors"]
):
    basis_signature_equal &= bool(
        even_sector["dimension"] == odd_sector["dimension"]
        and even_sector["splitter"] == odd_sector["splitter"]
        and abs(
            even_sector["old_central_eigenvalue"]
            - odd_sector["old_central_eigenvalue"]
        ) < mp.mpf("1e-140")
    )
    basis_distance = max(
        basis_distance,
        mp_difference_frobenius(
            even_sector["basis"], odd_sector["basis"]
        ),
    )
common_carrier_ok = bool(
    table_equal
    and actions_equal
    and seed_sets_equal
    and basis_signature_equal
    and basis_distance < mp.mpf("1e-140")
)
check(
    "both parities have literally common group, boundary and minimal-basis coordinates",
    common_carrier_ok,
    f"basis_distance={mp_text(basis_distance, 6)}",
)
common_seeds = tuple(sorted(prepared["even"]["boundary"]["old_seeds"]))

records = {}
runtime = {}
numeric_arrays = {}
for parity in ("even", "odd"):
    item = prepared[parity]
    index_data = item["index_data"]
    geometry = item["geometry"]
    kind_values = {
        "old": mp.mpf(1),
        "internal": state["lambda"] - state["rho"],
        "pole": -state["rho"],
        "new": state["lambda"] ** 2,
    }
    print(f"[{parity}] differentiating 20 local Lorentzian patterns", flush=True)
    pattern_cache, branch_control = library["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    branch_ok = bool(
        branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and branch_control["minimum_leading_minor"] > 0
        and branch_control["minimum_argument"] > mp.mpf("1e-6")
    )
    check(
        f"{parity}: all derivative levels retain the Lorentzian branch and hierarchy",
        branch_ok,
        (
            f"cross={mp_text(branch_control['maximum_cross'], 6)}, "
            f"proxy={mp_text(branch_control['maximum_proxy'], 6)}"
        ),
    )
    raw_kernels, kernel_control = library[
        "assemble_full_representative_kernels"
    ](index_data, geometry, pattern_cache)
    kernels = {
        "K01": richardson_kernel(
            raw_kernels["operational_primary"],
            raw_kernels["operational_shadow"],
        ),
        "K12": richardson_kernel(
            raw_kernels["operational_shadow"],
            raw_kernels["validation_primary"],
        ),
        "K23": richardson_kernel(
            raw_kernels["validation_primary"],
            raw_kernels["validation_shadow"],
        ),
    }
    richardson_imaginary = max(
        (
            abs(mp.im(value))
            for kernel in kernels.values()
            for value in kernel.values()
        ),
        default=mp.mpf(0),
    )
    kernel_ok = bool(
        set(raw_kernels) == set(VARIANTS)
        and all(len(kernel) > 0 for kernel in raw_kernels.values())
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-140")
        and richardson_imaginary < mp.mpf("1e-140")
    )
    check(
        f"{parity}: complete orbit-convolution Hessian kernels are physically real",
        kernel_ok,
        (
            f"entries={kernel_control['nonzero_entries']}, "
            f"imag={mp_text(richardson_imaginary, 6)}"
        ),
    )

    sector_records = []
    sector_runtime = []
    all_reciprocity = True
    all_regular = True
    all_canonical = True
    for sector_index, sector in enumerate(item["sectors"]):
        dimension = int(sector["dimension"])
        print(
            f"[{parity}] sector {sector_index + 1}/7 d={dimension}: "
            "projecting K and classifying J",
            flush=True,
        )
        blocks = {
            level: library["project_full_kernel"](kernel, sector)
            for level, kernel in kernels.items()
        }
        normalization = max(mp.mpf(1), mp_frobenius(blocks["K12"]))
        reciprocity_variation = max(
            mp_difference_frobenius(blocks["K01"], blocks["K12"]),
            mp_difference_frobenius(blocks["K12"], blocks["K23"]),
        ) / normalization
        reciprocity_residual = max(
            mp_difference_frobenius(block, block.H) / normalization
            for block in blocks.values()
        )
        reciprocity_epsilon = reciprocity_variation + CLASSIFIER_FLOOR
        reciprocity_ok = bool(
            reciprocity_residual <= 10 * reciprocity_epsilon
        )
        all_reciprocity &= reciprocity_ok

        rank_analysis = pre_legendre_analysis(blocks, dimension)
        all_regular &= rank_analysis["regular"]
        tangent_runtime = None
        tangent_record = None
        if rank_analysis["regular"]:
            old = expanded_types(0, 30, dimension)
            internal = expanded_types(30, 65, dimension)
            new = expanded_types(65, 95, dimension)
            output_mapping = expanded_output_mapping(
                item["boundary"]["mapping"], dimension
            )
            common_indices = expanded_common_indices(
                item["boundary"]["old_seeds"], common_seeds, dimension
            )
            ball_records = {}
            tangent_build_ok = True
            for level, block in blocks.items():
                built = canonical_tangent_ball(
                    block,
                    old,
                    internal,
                    new,
                    output_mapping,
                    common_indices,
                )
                tangent_build_ok &= bool(
                    built["tangent"] is not None
                    and not built["det_j"].contains(0)
                )
                if not tangent_build_ok:
                    continue
                midpoint, radii = acb_midpoint_and_radii(built["tangent"])
                defect_midpoint, defect_radii = acb_midpoint_and_radii(
                    built["defect"]
                )
                ball_records[level] = {
                    "midpoint": midpoint,
                    "radii": radii,
                    "defect_midpoint": defect_midpoint,
                    "defect_radii": defect_radii,
                }
                prefix = f"{parity}_sector{sector_index}_{level}"
                numeric_arrays[f"{prefix}_tangent_midpoint"] = midpoint
                numeric_arrays[f"{prefix}_tangent_radii"] = radii
                numeric_arrays[f"{prefix}_defect_midpoint"] = defect_midpoint
                numeric_arrays[f"{prefix}_defect_radii"] = defect_radii
            if tangent_build_ok and set(ball_records) == set(LEVELS):
                tangent_runtime = tangent_analysis(ball_records)
                tangent_record = tangent_runtime["record"]
                all_canonical &= tangent_runtime["canonicality_ok"]
            else:
                all_canonical = False
                tangent_record = {
                    "canonicality_ok": False,
                    "reason": "Flint solve failed after regular classification",
                }
        else:
            all_canonical = False

        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "constant_overlap": mp_text(sector["constant_overlap"]),
            "center_value": mp_text(sector["old_central_eigenvalue"]),
            "splitter": sector["splitter"],
            "reciprocity": {
                "residual": mp_text(reciprocity_residual),
                "variation": mp_text(reciprocity_variation),
                "epsilon": mp_text(reciprocity_epsilon),
                "passed": reciprocity_ok,
            },
            "pre_legendre": rank_analysis["record"],
            "tangent": tangent_record,
        })
        sector_runtime.append({
            "rank": rank_analysis,
            "tangent": tangent_runtime,
        })

    check(
        f"{parity}: all projected Hessian blocks retain calibrated reciprocity",
        all_reciprocity,
    )
    check(
        f"{parity}: all seven pre-Legendre sectors receive the frozen classifier",
        len(sector_records) == 7
        and all(
            record["pre_legendre"]["classification"]
            in {"REGULAR", "NUMERICALLY_OPEN"}
            for record in sector_records
        ),
        f"regular={sum(r['pre_legendre']['classification']=='REGULAR' for r in sector_records)}/7",
    )
    records[parity] = {
        "controls_ok": bool(
            item["carrier_ok"]
            and item["basis_ok"]
            and branch_ok
            and kernel_ok
            and all_reciprocity
        ),
        "all_regular": bool(all_regular),
        "all_canonical": bool(all_regular and all_canonical),
        "boundary_mapping": list(item["boundary"]["mapping"]),
        "old_orbit_seeds": [list(seed) for seed in item["boundary"]["old_seeds"]],
        "branch": {
            "minimum_leading_minor": mp_text(
                branch_control["minimum_leading_minor"]
            ),
            "minimum_argument": mp_text(branch_control["minimum_argument"]),
            "maximum_cross": mp_text(branch_control["maximum_cross"]),
            "maximum_proxy": mp_text(branch_control["maximum_proxy"]),
            "maximum_raw_imaginary": mp_text(
                branch_control["maximum_raw_angle_or_derivative_imaginary"]
            ),
        },
        "kernel": {
            "raw_nonzero_entries": kernel_control["nonzero_entries"],
            "raw_maximum_imaginary": mp_text(
                kernel_control["maximum_imaginary"]
            ),
            "richardson_maximum_imaginary": mp_text(richardson_imaginary),
        },
        "sector_control": {
            key: (
                value
                if isinstance(value, list)
                else mp_text(value)
            )
            for key, value in item["sector_control"].items()
        },
        "sectors": sector_records,
    }
    runtime[parity] = sector_runtime

all_regular = bool(all(records[p]["all_regular"] for p in ("even", "odd")))
all_canonical = bool(
    all(records[p]["all_canonical"] for p in ("even", "odd"))
)
schedule_records = []
actual_hostile_ok = True
if all_regular and all_canonical:
    for sector_index in range(7):
        comparison = schedule_compare(
            runtime["even"][sector_index]["tangent"],
            runtime["odd"][sector_index]["tangent"],
        )
        comparison["sector_index"] = sector_index
        comparison["dimension"] = records["even"]["sectors"][sector_index][
            "dimension"
        ]
        actual_hostile_ok &= bool(
            comparison["synthetic_detected"]
            and comparison["cyclic_output_map_detected"]
        )
        schedule_records.append(comparison)
    labels = Counter(item["label"] for item in schedule_records)
    if labels["SCHEDULE_DEPENDENT"]:
        schedule_outcome = "SCHEDULE_DEPENDENT"
    elif labels["SCHEDULE_OPEN"]:
        schedule_outcome = "SCHEDULE_OPEN"
    else:
        schedule_outcome = "SCHEDULE_ROBUST"
else:
    schedule_outcome = "NOT_EVALUATED"

check(
    "the actual schedule classifier detects its synthetic and boundary-map corruptions when applicable",
    actual_hostile_ok,
    f"schedule={schedule_outcome}",
)

base_controls_ok = bool(
    provenance_ok
    and registry_ok
    and accepted_inputs_ok
    and background_ok
    and scalar_control["passed"]
    and geometry_import_ok
    and common_carrier_ok
    and actual_hostile_ok
    and all(records[parity]["controls_ok"] for parity in ("even", "odd"))
)
if not base_controls_ok:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_CONTROL_FAILED"
elif not all_regular:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_RANK_OPEN"
elif not all_canonical:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED"
elif schedule_outcome == "SCHEDULE_DEPENDENT":
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_DEPENDENT_PRIMARY"
elif schedule_outcome == "SCHEDULE_OPEN":
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_OPEN"
else:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY"

check(
    "the result follows the preregistered outcome hierarchy",
    outcome in {
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_CONTROL_FAILED",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_RANK_OPEN",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_DEPENDENT_PRIMARY",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_OPEN",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY",
    },
    outcome,
)

deterministic_npz(NUMERIC_OUTPUT, numeric_arrays)
numeric_hash = sha256(NUMERIC_OUTPUT)
artifact = {
    "outcome": outcome,
    "status": "PRIMARY_ONLY_ADVERSARIAL_REPLICATION_REQUIRED_FOR_MATERIAL_RESULT",
    "tests": tests,
    "passed": passed,
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": input_hashes,
    },
    "background": {
        key: mp_text(value)
        for key, value in state.items()
        if key != "bracket_ok"
    },
    "full_boundary_configuration_dimension": 720,
    "full_boundary_phase_dimension": 1440,
    "minimal_sector_dimensions": [
        record["dimension"] for record in records["even"]["sectors"]
    ],
    "common_basis_distance": mp_text(basis_distance),
    "scalar_formula_control": scalar_control,
    "all_pre_legendre_sectors_regular": all_regular,
    "all_tangent_sectors_canonical": all_canonical,
    "schedule_outcome": schedule_outcome,
    "schedule_sectors": [
        {
            key: (
                float_text(value)
                if isinstance(value, float)
                else value
            )
            for key, value in record.items()
        }
        for record in schedule_records
    ],
    "parities": records,
    "numeric_archive": {
        "path": NUMERIC_OUTPUT.name,
        "sha256": numeric_hash,
        "arrays": {
            key: list(value.shape) for key, value in numeric_arrays.items()
        },
    },
    "classification": {
        "canonical_map": (
            "DERIVED_COMPUTATIONAL_PRIMARY"
            if outcome.endswith("_PRIMARY")
            else "OPEN_OR_FAILED"
        ),
        "schedule_comparison": schedule_outcome,
        "one_step_eigenvalue_spectrum": "NOT_COMPUTED",
        "physical_gauge_invariant_modes": "OPEN",
        "second_anisotropic_tick": "OPEN",
        "refinement_and_continuum": "OPEN",
        "wave_equation_limiting_speed_G_planck_particles": "NOT_DERIVED",
        "external_novelty": "OPEN",
    },
    "firewall": {
        "eigenvalues_computed": False,
        "continuum_target_parsed": False,
        "desired_rank_or_schedule_result_parsed": False,
        "full_suite_run": False,
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"OUTCOME: {outcome}")
print(f"SCHEDULE: {schedule_outcome}")
for parity in ("even", "odd"):
    regular_count = sum(
        record["pre_legendre"]["classification"] == "REGULAR"
        for record in records[parity]["sectors"]
    )
    canonical_count = sum(
        bool(record["tangent"] and record["tangent"].get("canonicality_ok"))
        for record in records[parity]["sectors"]
    )
    print(
        f"{parity}: regular={regular_count}/7, canonical={canonical_count}/7",
        flush=True,
    )
print(f"RESULT: {passed}/{tests} PASS")
print(f"NUMERIC SHA: {numeric_hash}")
if passed != tests:
    raise SystemExit(1)
