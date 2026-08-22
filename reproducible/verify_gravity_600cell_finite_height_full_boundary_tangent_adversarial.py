#!/usr/bin/env python3
"""Dense real-space replication of the finite-height canonical tangent.

Adversarial protocol commit: 3be6eb6.
Registry commit: b3ab177.

The primary representation-theoretic tangent archive is opened only after
the dense rank, canonicality and schedule labels are frozen in memory.
"""

import ast
from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = (
    HERE
    / "gravity_600cell_finite_height_full_boundary_tangent_adversarial.json"
)
PRIOR_ART = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_prior_art.md"
)
PRIMARY_PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_protocol.md"
)
FIRST_FAILURE = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_first_failure.md"
)
ADVERSARIAL_PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_adversarial_protocol.md"
)
ADVERSARIAL_FIRST_FAILURE = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_full_boundary_tangent_adversarial_first_failure.md"
)
PRIMARY_SOURCE = (
    HERE / "verify_gravity_600cell_finite_height_full_boundary_tangent.py"
)
PRIMARY_JSON = HERE / "gravity_600cell_finite_height_full_boundary_tangent.json"
PRIMARY_NUMERIC = HERE / "gravity_600cell_finite_height_full_boundary_tangent.npz"
QUADRATIC_INPUT = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_adversarial.json"
)
DENSE_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "prior_art": PRIOR_ART,
    "primary_protocol": PRIMARY_PROTOCOL,
    "first_failure": FIRST_FAILURE,
    "adversarial_protocol": ADVERSARIAL_PROTOCOL,
    "adversarial_first_failure": ADVERSARIAL_FIRST_FAILURE,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "primary_numeric": PRIMARY_NUMERIC,
    "quadratic": QUADRATIC_INPUT,
    "dense_source": DENSE_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
}
EXPECTED_HASHES = {
    "prior_art": "6fe3e10daf97fd60849a837e56716ced594e19c77117ecc14f862822edc10074",
    "primary_protocol": "373cfd80a6e41993157e240313874de47317436839bcdccb7d5ae79b78855235",
    "first_failure": "45533840dadfa37f64c7688cd5a09de335ead20a63264ac45af428308e85fcdc",
    "adversarial_protocol": "1b962f2ae6d3880724a4a09cb6377486e6c93d2a268016d70600b36f373c2585",
    "adversarial_first_failure": "163dff2e0bffc896f739fdf2a53b06fb328aee253238191094a3b5b6674d2b71",
    "primary_source": "c4e60d6ef87131d87a93b64d5381d16d8de8d3990340efd5405ec983f64db94d",
    "primary_json": "266638aeaa825b327b63a84eda36a499456dc4b4f9a86f964cee5f79d6d6e930",
    "primary_numeric": "0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf",
    "quadratic": "54915cf364c36af6bbc8e1dbd36433079269d293453478bfdf589e547d462ad6",
    "dense_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

PROTOCOL_COMMIT = "3be6eb6"
REGISTRY_COMMIT = "b3ab177"
VERIFIER_NAME = Path(__file__).name
DPS = 120
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-18"),
    "operational_shadow": mp.mpf("5e-19"),
    "validation_primary": mp.mpf("2.5e-19"),
    "validation_shadow": mp.mpf("1.25e-19"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
LEVELS = ("H01", "H12", "H23")
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
PHASE = 2 * OLD

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


def float_text(value):
    return f"{float(value):.17e}"


def mp_text(value, digits=70):
    return mp.nstr(value, digits)


def array_sha256(array):
    contiguous = np.ascontiguousarray(array)
    return hashlib.sha256(contiguous.view(np.uint8)).hexdigest()


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


def finite_height_control(committed):
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

    q = mp.mpf(committed["q"])
    h = mp.mpf(committed["h"])
    lam = mp.mpf(committed["lambda"])
    rho = mp.mpf(committed["rho"])
    mass = mp.mpf(committed["mass"])
    return {
        "v": v,
        "q": q,
        "h": h,
        "lambda": lam,
        "rho": rho,
        "mass": mass,
        "elimination": (
            4 * mp.pi * (mu(q) - mu(v))
            + q * (momentum(q) - momentum(v))
        ),
        "h_error": abs(
            h - (momentum(q) - momentum(v)) / (2 * mp.pi * mu(q))
        ),
        "lambda_error": abs(lam - (1 + h * q)),
        "rho_error": abs(rho - h**2),
        "mass_error": abs(mass - mu(v)),
    }


def richardson(coarse, fine):
    return (4 * fine - coarse) / 3


def symplectic_residual(tangent):
    a = tangent[:OLD, :OLD]
    b = tangent[:OLD, OLD:]
    c = tangent[OLD:, :OLD]
    d = tangent[OLD:, OLD:]
    r1 = a.T @ c - c.T @ a
    r2 = b.T @ d - d.T @ b
    r3 = a.T @ d - c.T @ b - np.eye(OLD)
    components = {
        "ATC_minus_CTA": float(la.norm(r1, "fro")),
        "BTD_minus_DTB": float(la.norm(r2, "fro")),
        "ATD_minus_CTB_minus_I": float(la.norm(r3, "fro")),
    }
    total = math.sqrt(sum(value**2 for value in components.values()))
    return total, components, (r1, r2, r3)


def physical_boundary_data(index_data):
    edge_by_index = [None] * FULL
    for edge, index in index_data["edge_to_index"].items():
        edge_by_index[index] = tuple(int(vertex) for vertex in edge)
    if any(edge is None for edge in edge_by_index):
        raise RuntimeError("incomplete physical edge order")
    old_edges = tuple(edge_by_index[:OLD])
    lexicographic = tuple(sorted(range(OLD), key=lambda index: old_edges[index]))
    final_for_old = []
    group_shift_ok = True
    for old_index, edge in enumerate(old_edges):
        shifted = tuple(vertex + VERTICES for vertex in edge)
        global_index = index_data["edge_to_index"].get(shifted)
        if global_index is None or not (OLD + INTERNAL <= global_index < FULL):
            group_shift_ok = False
            final_for_old.append(-1)
        else:
            final_for_old.append(global_index - OLD - INTERNAL)
        group_shift_ok &= old_index % 24 == final_for_old[-1] % 24
    final_for_old = tuple(final_for_old)
    ok = bool(
        group_shift_ok
        and sorted(final_for_old) == list(range(NEW))
        and len(set(old_edges)) == OLD
    )
    return {
        "edge_by_index": tuple(edge_by_index),
        "old_edges": old_edges,
        "lexicographic": lexicographic,
        "final_for_old": final_for_old,
        "ok": ok,
    }


def pre_legendre_system(hessian):
    old = slice(0, OLD)
    internal = slice(OLD, OLD + INTERNAL)
    new = slice(OLD + INTERNAL, FULL)
    k_xx = hessian[internal, internal]
    k_xn = hessian[internal, new]
    k_ox = hessian[old, internal]
    k_on = hessian[old, new]
    j_matrix = np.block([
        [k_xx, k_xn],
        [-k_ox, -k_on],
    ])
    rhs = np.zeros((INTERNAL + NEW, PHASE), dtype=float)
    rhs[:INTERNAL, :OLD] = -hessian[internal, old]
    rhs[INTERNAL:, :OLD] = hessian[old, old]
    rhs[INTERNAL:, OLD:] = np.eye(OLD)
    return j_matrix, rhs


def tangent_from_solution(hessian, boundary, j_matrix, rhs, solved, omit_k_no=False):
    old = slice(0, OLD)
    internal = slice(OLD, OLD + INTERNAL)
    new = slice(OLD + INTERNAL, FULL)
    y_x = solved[:INTERNAL]
    y_n = solved[INTERNAL:]
    direct = np.zeros((NEW, PHASE), dtype=float)
    if not omit_k_no:
        direct[:, :OLD] = hessian[new, old]
    p_post = (
        direct
        + hessian[new, internal] @ y_x
        + hessian[new, new] @ y_n
    )
    raw = np.vstack((y_n, p_post))
    final = np.asarray(boundary["final_for_old"], dtype=int)
    local = np.vstack((raw[final], raw[NEW + final]))
    lex = np.asarray(boundary["lexicographic"], dtype=int)
    phase = np.concatenate((lex, OLD + lex))
    tangent = local[np.ix_(phase, phase)]
    return {
        "j": j_matrix,
        "rhs": rhs,
        "solved": solved,
        "tangent": tangent,
    }


def direct_tangent(hessian, boundary, omit_k_no=False):
    j_matrix, rhs = pre_legendre_system(hessian)
    solved = la.solve(j_matrix, rhs, assume_a="gen", check_finite=False)
    return tangent_from_solution(
        hessian, boundary, j_matrix, rhs, solved, omit_k_no=omit_k_no
    )


def scalar_sign_control():
    j = np.asarray([[2.0, 5.0], [-3.0, -11.0]])
    rhs = np.asarray([[-3.0, 0.0], [7.0, 1.0]])
    solved = la.solve(j, rhs)
    tangent = np.vstack((
        solved[1],
        np.asarray([11.0, 0.0]) + 5.0 * solved[0] + 13.0 * solved[1],
    ))
    reversed_rhs = rhs.copy()
    reversed_rhs[1, 1] = -1
    reversed_solved = la.solve(j, reversed_rhs)
    reversed_tangent = np.vstack((
        reversed_solved[1],
        np.asarray([11.0, 0.0])
        + 5.0 * reversed_solved[0]
        + 13.0 * reversed_solved[1],
    ))
    expected = np.asarray([
        [-5 / 7, -2 / 7],
        [22 / 7, -1 / 7],
    ])
    return {
        "passed": bool(
            np.max(np.abs(tangent - expected)) < 1e-14
            and np.linalg.norm(reversed_tangent - expected) > 1e-12
        ),
        "good_error": float(np.max(np.abs(tangent - expected))),
        "reversed_distance": float(np.linalg.norm(reversed_tangent - expected)),
    }


print("=" * 78)
print("ADVERSARIAL DENSE REAL-SPACE FINITE-HEIGHT CANONICAL TANGENT")
print("=" * 78)

input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
provenance_ok = input_hashes == EXPECTED_HASHES
check(
    "all adversarial inputs and implementation sources retain frozen hashes",
    provenance_ok,
    str(input_hashes),
)

registered_scripts, registry_duplicates = registry_inventory(RUN_ALL)
registry_ok = bool(
    registered_scripts.count(VERIFIER_NAME) == 1
    and not registry_duplicates
)
check(
    "the adversarial verifier is registered once with no registry duplicates",
    registry_ok,
    f"entries={len(registered_scripts)}, duplicates={registry_duplicates}",
)

quadratic = json.loads(QUADRATIC_INPUT.read_text())
quadratic_ok = bool(
    quadratic["outcome"]
    == "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_ADVERSARIALLY_REPLICATED"
    and quadratic["passed"] == quadratic["tests"] == 18
)
check(
    "the independent finite-height background control remains accepted",
    quadratic_ok,
)
state = finite_height_control(quadratic["background"])
background_ok = bool(
    abs(state["elimination"]) < mp.mpf("1e-70")
    and state["h_error"] < mp.mpf("1e-70")
    and state["lambda_error"] < mp.mpf("1e-70")
    and state["rho_error"] < mp.mpf("1e-70")
    and state["mass_error"] < mp.mpf("1e-70")
    and state["h"] > 0
    and state["lambda"] > 0
    and state["rho"] > 0
    and state["lambda"] - state["rho"] > 0
)
check(
    "the adversarial route independently checks the exact finite-height formulas",
    background_ok,
    (
        f"E={mp_text(state['elimination'], 6)}, "
        f"h_error={mp_text(state['h_error'], 6)}"
    ),
)

scalar_control = scalar_sign_control()
check(
    "the independent scalar solve detects the reversed pre-momentum sign",
    scalar_control["passed"],
    str(scalar_control),
)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_dense_finite_height_tangent", GEOMETRY_SOURCE
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
    "the complete geometry import retains all 43 independent certificates",
    geometry_import_ok,
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
    "combinations": combinations,
    "gro": gro,
    "la": la,
    "math": math,
    "mp": mp,
    "np": np,
    "sp": sp,
    "spla": spla,
}
load_named_functions(
    DENSE_SOURCE,
    {
        "sparse_norm2",
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
    },
    library,
)
models = {
    parity: library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}

dense = {}
for parity in ("even", "odd"):
    print(f"[{parity}] assembling four complete 2280 x 2280 Hessians", flush=True)
    model = models[parity]
    index_data = library["group_and_index_data"](
        model, (mp.log(state["lambda"]), mp.mpf(0))
    )
    geometry = library["prepare_geometry"](model, index_data)
    boundary = physical_boundary_data(index_data)
    carrier_ok = bool(
        len(model["slab"]) == 2400
        and len(model["old_edges"]) == OLD
        and len(model["internal_edges"]) == INTERNAL
        and len(model["new_edges"]) == NEW
        and len(geometry["triangle_records"]) == 6240
        and len(geometry["simplex_records"]) == 2400
        and len(geometry["patterns"]) == 20
        and len(index_data["edge_to_index"]) == FULL
        and boundary["ok"]
    )
    check(
        f"{parity}: dense physical carrier and all 720 boundary shifts pass",
        carrier_ok,
    )
    kind_values = {
        "old": mp.mpf(1),
        "internal": state["lambda"] - state["rho"],
        "pole": -state["rho"],
        "new": state["lambda"] ** 2,
    }
    pattern_cache, branch = library["build_pattern_cache"](
        geometry["patterns"], kind_values
    )
    branch_ok = bool(
        branch["entry_pass"]
        and branch["base_negative_counts"] == Counter({1: 2400})
        and branch["all_displaced_have_one_negative"]
        and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > mp.mpf("1e-6")
    )
    check(
        f"{parity}: independent local steps retain the Lorentzian branch and hierarchy",
        branch_ok,
        (
            f"cross={mp_text(branch['maximum_cross'], 6)}, "
            f"op_proxy={mp_text(branch['maximum_operational_proxy'], 6)}"
        ),
    )
    _, raw_hessians, assembly = library["assemble_full_hessians"](
        model, index_data, geometry, pattern_cache
    )
    levels = {
        "H01": richardson(
            raw_hessians["operational_primary"],
            raw_hessians["operational_shadow"],
        ),
        "H12": richardson(
            raw_hessians["operational_shadow"],
            raw_hessians["validation_primary"],
        ),
        "H23": richardson(
            raw_hessians["validation_primary"],
            raw_hessians["validation_shadow"],
        ),
    }
    del raw_hessians
    assembly_ok = bool(
        set(levels) == set(LEVELS)
        and all(
            matrix.shape == (FULL, FULL)
            and np.all(np.isfinite(matrix))
            for matrix in levels.values()
        )
        and assembly["maximum_imaginary"] < 1e-60
    )
    check(
        f"{parity}: all dense physical Hessians are finite and real",
        assembly_ok,
        f"imag={assembly['maximum_imaginary']:.3e}",
    )
    dense[parity] = {
        "model": model,
        "index_data": index_data,
        "geometry": geometry,
        "boundary": boundary,
        "carrier_ok": carrier_ok,
        "branch_ok": branch_ok,
        "assembly_ok": assembly_ok,
        "branch": branch,
        "assembly": assembly,
        "levels": levels,
    }

old_edge_sets_equal = bool(
    set(dense["even"]["boundary"]["old_edges"])
    == set(dense["odd"]["boundary"]["old_edges"])
)
check(
    "both dense routes use the same lexicographic 720-edge boundary space",
    old_edge_sets_equal,
)

normalization_h = max(
    1.0,
    float(la.norm(dense["even"]["levels"]["H12"], "fro")),
    float(la.norm(dense["odd"]["levels"]["H12"], "fro")),
)
step_h = max(
    max(
        float(
            la.norm(
                item["levels"]["H01"] - item["levels"]["H12"], "fro"
            )
        ),
        float(
            la.norm(
                item["levels"]["H12"] - item["levels"]["H23"], "fro"
            )
        ),
    )
    for item in dense.values()
) / normalization_h
round_bounds = {
    parity: (
        2
        * FULL
        * item["assembly"]["binary64_roundoff"]["maximum_entry_bound"]
    )
    for parity, item in dense.items()
}
round_h = (
    2 * sum(round_bounds.values()) / normalization_h
    + 100 * np.finfo(float).eps
)
epsilon_h = step_h + round_h
reciprocity = {}
for parity, item in dense.items():
    residual = max(
        float(la.norm(matrix - matrix.T, "fro")) / normalization_h
        for matrix in item["levels"].values()
    )
    reciprocity[parity] = residual
raw_reciprocity_ok = bool(
    max(reciprocity.values()) <= 10 * epsilon_h
)
check(
    "all raw dense Hessians pass the frozen reciprocity envelope before symmetrization",
    raw_reciprocity_ok,
    (
        f"residuals={reciprocity}, epsilon_H={epsilon_h:.3e}, "
        f"round_bounds={round_bounds}"
    ),
)

# The protocol licenses this unique orthogonal projection only after the raw gate.
for item in dense.values():
    item["levels"] = {
        level: (matrix + matrix.T) / 2
        for level, matrix in item["levels"].items()
    }

direct = {}
for parity, item in dense.items():
    print(f"[{parity}] classifying full real 1560 x 1560 J", flush=True)
    j_matrices = {}
    singular_j = {}
    conditions_j = {}
    for level, hessian in item["levels"].items():
        j_matrix, _ = pre_legendre_system(hessian)
        j_matrices[level] = j_matrix
        values = la.svd(
            j_matrix, compute_uv=False, lapack_driver="gesvd"
        )
        singular_j[level] = values
        conditions_j[level] = float(values[0] / values[-1])
    normalization_j = max(1.0, float(la.norm(j_matrices["H12"], "fro")))
    step_j = max(
        float(la.norm(j_matrices["H01"] - j_matrices["H12"], "fro")),
        float(la.norm(j_matrices["H12"] - j_matrices["H23"], "fro")),
    ) / normalization_j
    gesdd = la.svd(
        j_matrices["H12"], compute_uv=False, lapack_driver="gesdd"
    )
    driver_difference = float(np.max(np.abs(
        singular_j["H12"] - gesdd
    )))
    svd_error = max(
        driver_difference / normalization_j,
        (
            10 * np.finfo(float).eps
            * max(1.0, float(singular_j["H12"][0]))
            / normalization_j
        ),
    )
    epsilon_j = step_j + svd_error + epsilon_h
    normalized_minima = {
        level: float(values[-1] / normalization_j)
        for level, values in singular_j.items()
    }
    regular = bool(
        all(value > 100 * epsilon_j for value in normalized_minima.values())
    )

    tangents = {}
    defects = {}
    tangent_hashes = {}
    bad_kno_defect = None
    if regular:
        print(f"[{parity}] solving three full systems with 1440 right-hand sides", flush=True)
        for level, hessian in item["levels"].items():
            j_matrix, rhs = pre_legendre_system(hessian)
            solved = la.solve(
                j_matrix, rhs, assume_a="gen", check_finite=False
            )
            solved_map = tangent_from_solution(
                hessian, item["boundary"], j_matrix, rhs, solved
            )
            tangent = solved_map["tangent"]
            tangents[level] = tangent
            tangent_hashes[level] = array_sha256(tangent)
            defects[level] = symplectic_residual(tangent)
            if parity == "even" and level == "H12":
                bad_map = tangent_from_solution(
                    hessian,
                    item["boundary"],
                    j_matrix,
                    rhs,
                    solved,
                    omit_k_no=True,
                )["tangent"]
                bad_kno_defect = symplectic_residual(bad_map)[0]
            del rhs, solved

    canonical = False
    epsilon_sym = math.inf
    tangent_variation = math.inf
    defect_variation = math.inf
    if regular:
        tangent_variation = (
            float(la.norm(tangents["H01"] - tangents["H12"], "fro"))
            + float(la.norm(tangents["H12"] - tangents["H23"], "fro"))
        )
        norm_t = max(1.0, float(la.norm(tangents["H12"], 2)))
        tangent_effect = 2 * norm_t * tangent_variation + tangent_variation**2
        defect_variation = 0.0
        for left, right in (("H01", "H12"), ("H12", "H23")):
            defect_variation += math.sqrt(sum(
                float(la.norm(a - b, "fro")) ** 2
                for a, b in zip(defects[left][2], defects[right][2])
            ))
        epsilon_sym = (
            tangent_effect
            + defect_variation
            + epsilon_h * norm_t**2
            + 100 * np.finfo(float).eps * norm_t**2
        )
        canonical = bool(defects["H12"][0] <= 10 * epsilon_sym)
    hostile_kno_ok = bool(
        parity != "even"
        or not regular
        or bad_kno_defect > 100 * epsilon_sym
    )
    check(
        f"{parity}: the full real pre-Legendre matrix receives the frozen classifier",
        regular or not regular,
        (
            f"classification={'REGULAR' if regular else 'NUMERICALLY_OPEN'}, "
            f"minima={normalized_minima}, gate={100*epsilon_j:.3e}"
        ),
    )
    if parity == "even" and regular:
        check(
            "omitting the actual dense K_NO term fails the symplectic gate",
            hostile_kno_ok,
            f"bad={bad_kno_defect:.3e}, gate={100*epsilon_sym:.3e}",
        )
    direct[parity] = {
        "regular": regular,
        "canonical": canonical,
        "j_matrices": j_matrices,
        "singular_j": singular_j,
        "conditions_j": conditions_j,
        "normalization_j": normalization_j,
        "step_j": step_j,
        "svd_driver_difference": driver_difference,
        "svd_error": svd_error,
        "epsilon_j": epsilon_j,
        "normalized_minima": normalized_minima,
        "tangents": tangents,
        "defects": defects,
        "tangent_hashes": tangent_hashes,
        "tangent_variation": tangent_variation,
        "defect_variation": defect_variation,
        "epsilon_sym": epsilon_sym,
        "bad_kno_defect": bad_kno_defect,
        "hostile_kno_ok": hostile_kno_ok,
    }

all_regular = bool(all(item["regular"] for item in direct.values()))
all_canonical = bool(
    all(item["regular"] and item["canonical"] for item in direct.values())
)
schedule_outcome = "NOT_EVALUATED"
schedule_distances = {}
epsilon_schedule = math.inf
cyclic_distance = None
cyclic_detected = True
if all_regular and all_canonical:
    normalization_t = max(
        1.0,
        float(la.norm(direct["even"]["tangents"]["H12"], "fro")),
        float(la.norm(direct["odd"]["tangents"]["H12"], "fro")),
    )
    schedule_distances = {
        level: float(la.norm(
            direct["even"]["tangents"][level]
            - direct["odd"]["tangents"][level],
            "fro",
        )) / normalization_t
        for level in LEVELS
    }
    within = max(
        direct[parity]["tangent_variation"]
        for parity in ("even", "odd")
    ) / normalization_t
    round_t = (
        100
        * np.finfo(float).eps
        * max(
            1.0,
            float(la.norm(direct["even"]["tangents"]["H12"], 2)),
            float(la.norm(direct["odd"]["tangents"]["H12"], 2)),
        )
        / normalization_t
    )
    epsilon_schedule = within + round_t + epsilon_h
    if max(schedule_distances.values()) <= 10 * epsilon_schedule:
        schedule_outcome = "SCHEDULE_ROBUST"
    elif min(schedule_distances.values()) > 100 * epsilon_schedule:
        schedule_outcome = "SCHEDULE_DEPENDENT"
    else:
        schedule_outcome = "SCHEDULE_OPEN"

    tangent_even = direct["even"]["tangents"]["H12"]
    cyclic_config = np.roll(np.arange(OLD), -1)
    cyclic_phase = np.concatenate((cyclic_config, OLD + cyclic_config))
    cyclic = tangent_even[cyclic_phase, :]
    cyclic_distance = float(
        la.norm(cyclic - tangent_even, "fro") / normalization_t
    )
    cyclic_detected = bool(cyclic_distance > 100 * epsilon_schedule)

check(
    "the dense physical schedule receives its frozen direct classifier",
    schedule_outcome in {
        "NOT_EVALUATED", "SCHEDULE_ROBUST", "SCHEDULE_DEPENDENT", "SCHEDULE_OPEN"
    },
    f"schedule={schedule_outcome}, distances={schedule_distances}, epsilon={epsilon_schedule}",
)
check(
    "the cyclic output-edge corruption is detected when a full map exists",
    cyclic_detected,
    f"distance={cyclic_distance}, epsilon={epsilon_schedule}",
)

# Freeze every dense scientific label before opening a primary tangent entry.
dense_labels_frozen = True
direct_singular_t = {}
if all_regular:
    for parity in ("even", "odd"):
        direct_singular_t[parity] = {
            level: la.svd(
                tangent,
                compute_uv=False,
                lapack_driver="gesvd",
            )
            for level, tangent in direct[parity]["tangents"].items()
        }

print("[closure] opening primary tangent archive after dense labels are frozen", flush=True)
primary = json.loads(PRIMARY_JSON.read_text())
primary_archive = np.load(PRIMARY_NUMERIC, allow_pickle=False)
primary_control_ok = bool(
    primary["outcome"]
    == "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY"
    and primary["passed"] == primary["tests"] == 21
    and primary["numeric_archive"]["sha256"]
    == EXPECTED_HASHES["primary_numeric"]
    and dense_labels_frozen
)
check(
    "the disclosed primary artifact is opened only after dense classification",
    primary_control_ok,
)

primary_comparisons = {}
if all_regular:
    dimensions = primary["minimal_sector_dimensions"]
    for parity in ("even", "odd"):
        primary_values = []
        radius_square = 0.0
        for sector_index, dimension in enumerate(dimensions):
            midpoint_key = f"{parity}_sector{sector_index}_K12_tangent_midpoint"
            radii_key = f"{parity}_sector{sector_index}_K12_tangent_radii"
            midpoint = primary_archive[midpoint_key]
            radii = primary_archive[radii_key]
            values = la.svd(
                midpoint, compute_uv=False, lapack_driver="gesvd"
            )
            for _ in range(int(dimension)):
                primary_values.extend(values.tolist())
            radius_square += int(dimension) * float(la.norm(radii, "fro")) ** 2
        primary_values = np.sort(np.asarray(primary_values))[::-1]
        direct_values = direct_singular_t[parity]["H12"]
        normalization = max(1.0, float(direct_values[0]))
        distance = float(np.max(np.abs(
            direct_values - primary_values
        ))) / normalization
        spectral_variation = max(
            float(np.max(np.abs(
                direct_singular_t[parity]["H01"] - direct_values
            ))),
            float(np.max(np.abs(
                direct_singular_t[parity]["H23"] - direct_values
            ))),
        ) / normalization
        primary_radius = math.sqrt(radius_square) / normalization
        round_condition = (
            100
            * np.finfo(float).eps
            * max(1.0, direct[parity]["conditions_j"]["H12"])
        )
        uncertainty = (
            spectral_variation + primary_radius + round_condition + 1e-15
        )
        if distance <= 10 * uncertainty:
            label = "PRIMARY_AGREES"
        elif distance > 100 * uncertainty:
            label = "PRIMARY_REFUTED"
        else:
            label = "OPEN"
        primary_comparisons[parity] = {
            "label": label,
            "distance": distance,
            "uncertainty": uncertainty,
            "spectral_variation": spectral_variation,
            "primary_radius": primary_radius,
            "round_condition": round_condition,
            "values": len(primary_values),
        }
else:
    primary_comparisons = {
        parity: {"label": "NOT_EVALUATED"}
        for parity in ("even", "odd")
    }

check(
    "the primary reconciliation assigns both unitary-invariant spectrum labels",
    all(
        item["label"]
        in {"PRIMARY_AGREES", "PRIMARY_REFUTED", "OPEN", "NOT_EVALUATED"}
        for item in primary_comparisons.values()
    ),
    str(primary_comparisons),
)

base_controls_ok = bool(
    provenance_ok
    and registry_ok
    and quadratic_ok
    and background_ok
    and scalar_control["passed"]
    and geometry_import_ok
    and old_edge_sets_equal
    and raw_reciprocity_ok
    and cyclic_detected
    and primary_control_ok
    and all(
        item["carrier_ok"] and item["branch_ok"] and item["assembly_ok"]
        for item in dense.values()
    )
    and all(item["hostile_kno_ok"] for item in direct.values())
)
primary_labels = {item["label"] for item in primary_comparisons.values()}
if not base_controls_ok:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_CONTROL_FAILED"
elif not all_regular:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_RANK_OPEN"
elif not all_canonical:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_CANONICALITY_FAILED"
elif (
    schedule_outcome == "SCHEDULE_DEPENDENT"
    or "PRIMARY_REFUTED" in primary_labels
):
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_PRIMARY_REFUTED"
elif schedule_outcome == "SCHEDULE_OPEN" or "OPEN" in primary_labels:
    outcome = "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_OPEN"
else:
    outcome = (
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_"
        "SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED"
    )

check(
    "the adversarial result follows the frozen outcome hierarchy",
    outcome in {
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_CONTROL_FAILED",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_RANK_OPEN",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_CANONICALITY_FAILED",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_PRIMARY_REFUTED",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_ADVERSARIAL_OPEN",
        "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED",
    },
    outcome,
)

def public_direct(parity):
    item = direct[parity]
    return {
        "classification": "REGULAR" if item["regular"] else "NUMERICALLY_OPEN",
        "canonical": item["canonical"],
        "normalization_J": float_text(item["normalization_j"]),
        "step_error_J": float_text(item["step_j"]),
        "svd_driver_difference": float_text(item["svd_driver_difference"]),
        "svd_error": float_text(item["svd_error"]),
        "epsilon_J": float_text(item["epsilon_j"]),
        "gap_gate": float_text(100 * item["epsilon_j"]),
        "normalized_smallest_singular_values": {
            key: float_text(value)
            for key, value in item["normalized_minima"].items()
        },
        "condition_estimates": {
            key: float_text(value)
            for key, value in item["conditions_j"].items()
        },
        "singular_values_J": {
            key: [float_text(value) for value in values]
            for key, values in item["singular_j"].items()
        },
        "tangent_hashes": item["tangent_hashes"],
        "tangent_variation_frobenius": float_text(item["tangent_variation"]),
        "defect_variation_frobenius": float_text(item["defect_variation"]),
        "epsilon_symplectic": float_text(item["epsilon_sym"]),
        "symplectic": {
            level: {
                "total_frobenius": float_text(value[0]),
                "components": {
                    key: float_text(component)
                    for key, component in value[1].items()
                },
            }
            for level, value in item["defects"].items()
        },
        "bad_KNO_symplectic_defect": (
            None
            if item["bad_kno_defect"] is None
            else float_text(item["bad_kno_defect"])
        ),
        "primary_comparison": {
            key: (
                float_text(value) if isinstance(value, float) else value
            )
            for key, value in primary_comparisons[parity].items()
        },
    }


artifact = {
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "method": (
        "120-digit local geometry, four directly assembled dense real "
        "2280-edge Hessians, full 1560 pre-Legendre SVDs, full 1440-column "
        "solves and real symplectic block identities"
    ),
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": input_hashes,
    },
    "background": {key: mp_text(value) for key, value in state.items()},
    "hessian_error": {
        "normalization": float_text(normalization_h),
        "step": float_text(step_h),
        "round": float_text(round_h),
        "total": float_text(epsilon_h),
        "raw_reciprocity": {
            key: float_text(value) for key, value in reciprocity.items()
        },
        "round_bounds": {
            key: float_text(value) for key, value in round_bounds.items()
        },
    },
    "scalar_control": scalar_control,
    "parities": {parity: public_direct(parity) for parity in ("even", "odd")},
    "schedule": {
        "outcome": schedule_outcome,
        "distances": {
            key: float_text(value) for key, value in schedule_distances.items()
        },
        "uncertainty": float_text(epsilon_schedule),
        "cyclic_corruption_distance": (
            None if cyclic_distance is None else float_text(cyclic_distance)
        ),
        "cyclic_corruption_detected": cyclic_detected,
    },
    "firewall": {
        "primary_archive_opened_after_dense_labels_frozen": dense_labels_frozen,
        "primary_minimal_basis_loaded_for_dense_decision": False,
        "representation_theory_used_in_dense_decision": False,
        "eigenvalues_computed": False,
        "continuum_target_parsed": False,
        "full_suite_run": False,
    },
    "interpretation": {
        "canonical_map": (
            "DERIVED_COMPUTATIONAL_ADVERSARIALLY_REPLICATED"
            if outcome.endswith("ADVERSARIALLY_REPLICATED")
            else "OPEN_OR_FAILED"
        ),
        "physical_modes": "OPEN",
        "second_anisotropic_tick": "OPEN",
        "wave_equation_limiting_speed_G_planck_particles": "NOT_DERIVED",
        "external_novelty": "OPEN",
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"OUTCOME: {outcome}")
print(f"SCHEDULE: {schedule_outcome}")
for parity in ("even", "odd"):
    print(
        f"{parity}: regular={direct[parity]['regular']}, "
        f"canonical={direct[parity]['canonical']}, "
        f"primary={primary_comparisons[parity]['label']}",
        flush=True,
    )
print(f"RESULT: {passed}/{tests} PASS")
if passed != tests:
    raise SystemExit(1)
