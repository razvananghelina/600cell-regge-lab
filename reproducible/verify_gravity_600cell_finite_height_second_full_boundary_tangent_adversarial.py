#!/usr/bin/env python3
"""Dense adversarial replication of the second finite-height tangent.

Protocol commit: eceda30.
Registry commit: 4ce8ee9.

The scientific decision is made in the complete lexicographically ordered
720-edge boundary basis. Minimal representation sectors and primary tangent
entries are opened only after every dense label has been frozen in memory.
No tangent spectrum is computed.
"""

import ast
from collections import Counter
from itertools import combinations
import contextlib
import gc
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_finite_height_second_full_boundary_tangent_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_second_full_boundary_tangent_adversarial_protocol.md"
FIRST_RUN_FAILURE = ROOT / "docs/gravity/gravity_600cell_second_full_boundary_tangent_adversarial_first_run_failure.md"
FIRST_RUN_CORRECTION = ROOT / "docs/gravity/gravity_600cell_second_full_boundary_tangent_adversarial_first_run_correction.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_finite_height_second_full_boundary_tangent.py"
PRIMARY_JSON = HERE / "gravity_600cell_finite_height_second_full_boundary_tangent.json"
PRIMARY_NUMERIC = HERE / "gravity_600cell_finite_height_second_full_boundary_tangent.npz"
FIRST_DENSE_SOURCE = HERE / "verify_gravity_600cell_finite_height_full_boundary_tangent_adversarial.py"
FIRST_DENSE_JSON = HERE / "gravity_600cell_finite_height_full_boundary_tangent_adversarial.json"
FIRST_PRIMARY_SOURCE = HERE / "verify_gravity_600cell_finite_height_full_boundary_tangent.py"
FIRST_PRIMARY_NUMERIC = HERE / "gravity_600cell_finite_height_full_boundary_tangent.npz"
COMPOSITION_ADVERSARIAL = HERE / "gravity_600cell_finite_height_composition_adversarial.json"
COMPOSITION_ADVERSARIAL_SOURCE = HERE / "verify_gravity_600cell_finite_height_composition_adversarial.py"
DENSE_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
OLD_TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "protocol": PROTOCOL,
    "first_run_failure": FIRST_RUN_FAILURE,
    "first_run_correction": FIRST_RUN_CORRECTION,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "primary_numeric": PRIMARY_NUMERIC,
    "first_dense_source": FIRST_DENSE_SOURCE,
    "first_dense_json": FIRST_DENSE_JSON,
    "first_primary_source": FIRST_PRIMARY_SOURCE,
    "first_primary_numeric": FIRST_PRIMARY_NUMERIC,
    "composition_adversarial": COMPOSITION_ADVERSARIAL,
    "composition_adversarial_source": COMPOSITION_ADVERSARIAL_SOURCE,
    "dense_source": DENSE_SOURCE,
    "old_tangent_source": OLD_TANGENT_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
}
EXPECTED_HASHES = {
    "protocol": "d6f3a27771dc258f3c1c0cd84d99780a45436bb7dcd544e43289a1fe188e97fc",
    "first_run_failure": "56bf30a4046c8b66162c2cd7aaf77acecc92edd6d88c655123b81b9a10ddeab5",
    "first_run_correction": "725acff64492b65b5d3dc7fe466755b0c1ac34ea6a80c7442106c8b90bcb32f3",
    "primary_source": "8c29c66eb4ec253229685cdbe56eb0371fb00860f0f2e15804fca0c7c64ec536",
    "primary_json": "f97db13031fd366b74e7d327abf61d8d23c24ee2889e500a2cfc747ed7dd6990",
    "primary_numeric": "c80a1af93deddc9526c362a8a76eec5dfc4b8360440a1a69f73ca61a419aac9a",
    "first_dense_source": "6e4fdbe6822ac024d2bad9aa22fe73f85aa9805a018146482ace47cfbefc43b6",
    "first_dense_json": "ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7",
    "first_primary_source": "c4e60d6ef87131d87a93b64d5381d16d8de8d3990340efd5405ec983f64db94d",
    "first_primary_numeric": "0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf",
    "composition_adversarial": "d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9",
    "composition_adversarial_source": "8395e921ab1c1f518abb567a114f1eb8bfdf2068be031bff55c8d2f0cff56c2b",
    "dense_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "old_tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

PROTOCOL_COMMIT = "eceda30"
REGISTRY_COMMIT = "4ce8ee9"
VERIFIER_NAME = Path(__file__).name
DPS = 120
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-18"),
    "operational_shadow": mp.mpf("5e-19"),
    "validation_primary": mp.mpf("2.5e-19"),
    "validation_shadow": mp.mpf("1.25e-19"),
}
RAW_ORDER = tuple(DERIVATIVE_STEPS)
LEVELS = ("H01", "H12", "H23")
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
VERTICES = 120
OLD = 720
INTERNAL = 840
NEW = 720
FULL = OLD + INTERNAL + NEW
PHASE = 2 * OLD
MACHINE_EPS = np.finfo(float).eps
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


def array_sha256(array):
    value = np.ascontiguousarray(array)
    return hashlib.sha256(value.view(np.uint8)).hexdigest()


def float_text(value):
    return f"{float(value):.17e}"


def mp_text(value, digits=75):
    return mp.nstr(value, digits)


def load_named_functions(path, names, namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
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
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "scripts" for target in node.targets
        ):
            scripts = ast.literal_eval(node.value)
            break
    if scripts is None:
        raise RuntimeError("run_all.py has no literal scripts registry")
    counts = Counter(scripts)
    return scripts, sorted(name for name, count in counts.items() if count != 1)


def epsilon(value):
    square = value**2
    return 2 * mp.pi - 5 * mp.acos((square + 2) / (2 * (square + 3)))


def mu(value):
    return 180 * epsilon(value) / (mp.pi * mp.sqrt(value**2 + 4))


def momentum(value):
    square = value**2
    return (
        180 * value * epsilon(value) / mp.sqrt(square + 4)
        - 600 * mp.sqrt(3) * mp.asinh(value / mp.sqrt(8 * (square + 3)))
    )


def bisect_root(function, left, right):
    left, right = mp.mpf(left), mp.mpf(right)
    left_value, right_value = function(left), function(right)
    bracket_ok = bool(left_value * right_value < 0)
    for _ in range(650):
        middle = (left + right) / 2
        middle_value = function(middle)
        if middle_value == 0:
            left = right = middle
            break
        if left_value * middle_value < 0:
            right, right_value = middle, middle_value
        else:
            left, left_value = middle, middle_value
        if right - left < mp.mpf("1e-105"):
            break
    return (left + right) / 2, right - left, bracket_ok


def reconstruct_history(adversarial):
    v = mp.mpf(3) / 2
    m0, pi0 = mu(v), momentum(v)

    def advance(mass, incoming, bracket):
        equation = lambda q: 4 * mp.pi * (mu(q) - mass) + q * (momentum(q) - incoming)
        q, width, bracket_ok = bisect_root(equation, *bracket)
        h = (momentum(q) - incoming) / (2 * mp.pi * mu(q))
        ratio = 1 + h * q
        return {
            "q": q, "h": h, "r": ratio, "m_in": mass, "pi_in": incoming,
            "m_out": mass / ratio,
            "pi_out": momentum(q) + 2 * mp.pi * h * mu(q) / ratio,
            "residual": equation(q), "width": width, "bracket_ok": bracket_ok,
        }

    first = advance(m0, pi0, (9, 10))
    second = advance(first["m_out"], first["pi_out"], (31, 32))
    run = next(item for item in adversarial["precision_runs"] if item["target_precision"] == 180)
    branch = run["branches"]["B"]
    committed_error = max(
        abs(first["q"] - mp.mpf(run["q1"])), abs(first["h"] - mp.mpf(run["h1"])),
        abs(first["r"] - mp.mpf(run["L1"])), abs(first["m_out"] - mp.mpf(run["m1"])),
        abs(first["pi_out"] - mp.mpf(run["pi1"])), abs(second["q"] - mp.mpf(branch["q2"])),
        abs(second["h"] - mp.mpf(branch["h2"])), abs(second["r"] - mp.mpf(branch["ratio"])),
    )
    junction_error = max(
        abs(first["m_out"] - second["m_in"]), abs(first["pi_out"] - second["pi_in"]),
        abs(m0 - first["r"] * first["m_out"]),
    )
    return {"v": v, "m0": m0, "pi0": pi0, "first": first, "second": second,
            "committed_error": committed_error, "junction_error": junction_error}


def richardson(coarse, fine):
    return (4 * fine - coarse) / 3


def physical_boundary_data(index_data):
    edge_by_index = [None] * FULL
    for edge, index in index_data["edge_to_index"].items():
        edge_by_index[index] = tuple(int(vertex) for vertex in edge)
    if any(edge is None for edge in edge_by_index):
        raise RuntimeError("incomplete physical edge order")
    old_edges = tuple(edge_by_index[:OLD])
    lexicographic = tuple(sorted(range(OLD), key=lambda index: old_edges[index]))
    final_for_old, group_shift_ok = [], True
    for old_index, edge in enumerate(old_edges):
        shifted = tuple(vertex + VERTICES for vertex in edge)
        global_index = index_data["edge_to_index"].get(shifted)
        if global_index is None or not (OLD + INTERNAL <= global_index < FULL):
            group_shift_ok, final = False, -1
        else:
            final = global_index - OLD - INTERNAL
        final_for_old.append(final)
        group_shift_ok &= old_index % 24 == final % 24
    return {
        "old_edges": old_edges, "lexicographic": lexicographic,
        "final_for_old": tuple(final_for_old),
        "ok": bool(group_shift_ok and sorted(final_for_old) == list(range(NEW)) and len(set(old_edges)) == OLD),
    }


def pre_legendre_system(hessian):
    old, internal, new = slice(0, OLD), slice(OLD, OLD + INTERNAL), slice(OLD + INTERNAL, FULL)
    j_matrix = np.block([
        [hessian[internal, internal], hessian[internal, new]],
        [-hessian[old, internal], -hessian[old, new]],
    ])
    rhs = np.zeros((INTERNAL + NEW, PHASE), dtype=float)
    rhs[:INTERNAL, :OLD] = -hessian[internal, old]
    rhs[INTERNAL:, :OLD] = hessian[old, old]
    rhs[INTERNAL:, OLD:] = np.eye(OLD)
    return j_matrix, rhs


def tangent_from_solution(hessian, boundary, solved, omit_k_no=False):
    old, internal, new = slice(0, OLD), slice(OLD, OLD + INTERNAL), slice(OLD + INTERNAL, FULL)
    y_x, y_n = solved[:INTERNAL], solved[INTERNAL:]
    direct = np.zeros((NEW, PHASE), dtype=float)
    if not omit_k_no:
        direct[:, :OLD] = hessian[new, old]
    p_post = direct + hessian[new, internal] @ y_x + hessian[new, new] @ y_n
    raw = np.vstack((y_n, p_post))
    final = np.asarray(boundary["final_for_old"], dtype=int)
    local = np.vstack((raw[final], raw[NEW + final]))
    lex = np.asarray(boundary["lexicographic"], dtype=int)
    phase = np.concatenate((lex, OLD + lex))
    return local[np.ix_(phase, phase)]


def symplectic_residual(tangent):
    a, b = tangent[:OLD, :OLD], tangent[:OLD, OLD:]
    c, d = tangent[OLD:, :OLD], tangent[OLD:, OLD:]
    pieces = (a.T @ c - c.T @ a, b.T @ d - d.T @ b, a.T @ d - c.T @ b - np.eye(OLD))
    components = {name: float(la.norm(value, "fro")) for name, value in zip(
        ("ATC_minus_CTA", "BTD_minus_DTB", "ATD_minus_CTB_minus_I"), pieces
    )}
    return math.sqrt(sum(value**2 for value in components.values())), components


def scale_tangent(tangent, factor):
    result = tangent.copy()
    result[:OLD, OLD:] /= float(factor)
    result[OLD:, :OLD] *= float(factor)
    return result


def classify(distance, uncertainty):
    if distance <= 10 * uncertainty:
        return "AGREES"
    if distance > 100 * uncertainty:
        return "REFUTED"
    return "OPEN"


def family_variation(family):
    middle = family["H12"]
    return float(la.norm(family["H01"] - middle, "fro") + la.norm(middle - family["H23"], "fro"))


def make_dense_library(gro, mass, rho0, length_square):
    namespace = {
        "Counter": Counter, "DPS": DPS, "DERIVATIVE_STEPS": DERIVATIVE_STEPS, "I": I,
        "LOCAL_EDGES": LOCAL_EDGES, "LOCAL_HINGES": LOCAL_HINGES,
        "LOCAL_HINGE_INDEX": LOCAL_HINGE_INDEX, "L0_SQUARE": mp.mpf(length_square),
        "MASS": mp.mpf(mass), "RHO0": mp.mpf(rho0), "combinations": combinations,
        "gro": gro, "la": la, "math": math, "mp": mp, "np": np, "sp": sp, "spla": spla,
    }
    load_named_functions(DENSE_SOURCE, {
        "sparse_norm2", "log_minus", "signed_volume_square", "angle_record", "area_data",
        "extended_edge_image", "orbit_sort_key", "augment_boundary_orbits", "group_and_index_data",
        "prepare_geometry", "build_pattern_cache", "assemble_full_hessians",
    }, namespace)
    return namespace


def assemble_stage(name, parity, library, model, ratio):
    print(f"[{parity}/{name}] preparing complete carrier and four Hessians", flush=True)
    started = time.time()
    index_data = library["group_and_index_data"](model, (mp.log(ratio), mp.mpf(0)))
    geometry = library["prepare_geometry"](model, index_data)
    boundary = physical_boundary_data(index_data)
    carrier_ok = bool(
        len(model["slab"]) == 2400 and len(model["old_edges"]) == OLD
        and len(model["internal_edges"]) == INTERNAL and len(model["new_edges"]) == NEW
        and len(geometry["triangle_records"]) == 6240 and len(geometry["simplex_records"]) == 2400
        and len(geometry["patterns"]) == 20 and len(index_data["edge_to_index"]) == FULL and boundary["ok"]
    )
    values = {
        "old": library["L0_SQUARE"],
        "internal": ratio * library["L0_SQUARE"] - library["RHO0"],
        "pole": -library["RHO0"], "new": ratio**2 * library["L0_SQUARE"],
    }
    pattern_cache, branch = library["build_pattern_cache"](geometry["patterns"], values)
    branch_ok = bool(
        branch["entry_pass"] and branch["base_negative_counts"] == Counter({1: 2400})
        and branch["all_displaced_have_one_negative"] and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > mp.mpf("1e-6")
    )
    _, raw, assembly = library["assemble_full_hessians"](model, index_data, geometry, pattern_cache)
    levels = {
        "H01": richardson(raw[RAW_ORDER[0]], raw[RAW_ORDER[1]]),
        "H12": richardson(raw[RAW_ORDER[1]], raw[RAW_ORDER[2]]),
        "H23": richardson(raw[RAW_ORDER[2]], raw[RAW_ORDER[3]]),
    }
    normalization = max(1.0, float(la.norm(levels["H12"], "fro")))
    step = max(float(la.norm(levels["H01"] - levels["H12"], "fro")),
               float(la.norm(levels["H12"] - levels["H23"], "fro"))) / normalization
    round_error = 2 * FULL * assembly["binary64_roundoff"]["maximum_entry_bound"] / normalization + 100 * MACHINE_EPS
    uncertainty = step + round_error
    raw_reciprocity = {key: float(la.norm(value - value.T, "fro")) / normalization for key, value in raw.items()}
    richardson_reciprocity = {key: float(la.norm(value - value.T, "fro")) / normalization for key, value in levels.items()}
    reciprocity_ok = max((*raw_reciprocity.values(), *richardson_reciprocity.values())) <= 10 * uncertainty
    return {
        "name": name, "parity": parity, "index_data": index_data, "boundary": boundary,
        "carrier_ok": carrier_ok, "branch_ok": branch_ok, "branch": branch, "assembly": assembly,
        "assembly_ok": bool(all(value.shape == (FULL, FULL) and np.all(np.isfinite(value))
                                for value in (*raw.values(), *levels.values())) and assembly["maximum_imaginary"] < 1e-60),
        "raw": raw, "levels": levels, "normalization_h": normalization, "step_h": step,
        "round_h": round_error, "epsilon_h": uncertainty, "raw_reciprocity": raw_reciprocity,
        "richardson_reciprocity": richardson_reciprocity, "reciprocity_ok": bool(reciprocity_ok),
        "seconds": time.time() - started,
    }


def scale_hessian_records(normalized, physical, factor):
    records = []
    def one(kind, key, left_family, right_family, order):
        left, right = left_family[key], right_family[key]
        normalization = max(1.0, float(la.norm(left, "fro")), abs(float(factor)) * float(la.norm(right, "fro")))
        distance = float(la.norm(left - float(factor) * right, "fro")) / normalization
        position = order.index(key)
        neighbor = order[position + 1] if position + 1 < len(order) else order[position - 1]
        variation = (float(la.norm(left_family[key] - left_family[neighbor], "fro"))
                     + abs(float(factor)) * float(la.norm(right_family[key] - right_family[neighbor], "fro"))) / normalization
        uncertainty = variation + physical["round_h"] + normalized["round_h"] + 100 * MACHINE_EPS
        return {"kind": kind, "level": key, "distance": distance, "uncertainty": uncertainty,
                "label": classify(distance, uncertainty)}
    records.extend(one("raw", key, physical["raw"], normalized["raw"], RAW_ORDER) for key in RAW_ORDER)
    records.extend(one("richardson", key, physical["levels"], normalized["levels"], LEVELS) for key in LEVELS)
    labels = {record["label"] for record in records}
    outcome = "SCALE_LIFT_REFUTED" if "REFUTED" in labels else "SCALE_LIFT_OPEN" if "OPEN" in labels else "SCALE_LIFT_CONFIRMED"
    return records, outcome


def analyze_stage(stage, hostile_kno=False):
    stage["levels"] = {key: (value + value.T) / 2 for key, value in stage["levels"].items()}
    j_matrices, singular, singular_gesdd, conditions = {}, {}, {}, {}
    for level, hessian in stage["levels"].items():
        print(f"[{stage['parity']}/{stage['name']}/{level}] two J SVD drivers", flush=True)
        j_matrix, _ = pre_legendre_system(hessian)
        j_matrices[level] = j_matrix
        singular[level] = la.svd(j_matrix, compute_uv=False, lapack_driver="gesvd")
        singular_gesdd[level] = la.svd(j_matrix, compute_uv=False, lapack_driver="gesdd")
        conditions[level] = float(singular[level][0] / singular[level][-1])
    normalization_j = max(1.0, float(la.norm(j_matrices["H12"], "fro")))
    step_j = max(float(la.norm(j_matrices["H01"] - j_matrices["H12"], "fro")),
                 float(la.norm(j_matrices["H12"] - j_matrices["H23"], "fro"))) / normalization_j
    driver_error = max(float(np.max(np.abs(singular[key] - singular_gesdd[key]))) / normalization_j for key in LEVELS)
    binary_error = 10 * MACHINE_EPS * max(float(singular[key][0]) for key in LEVELS) / normalization_j
    epsilon_j = step_j + driver_error + binary_error + stage["epsilon_h"]
    minima = {key: float(singular[key][-1] / normalization_j) for key in LEVELS}
    regular = all(value > 100 * epsilon_j for value in minima.values())
    tangents, defects, hashes, bad_kno_defect = {}, {}, {}, None
    if regular:
        for level, hessian in stage["levels"].items():
            print(f"[{stage['parity']}/{stage['name']}/{level}] solving 1440 RHS", flush=True)
            j_matrix, rhs = pre_legendre_system(hessian)
            solved = la.solve(j_matrix, rhs, assume_a="gen", check_finite=False)
            tangent = tangent_from_solution(hessian, stage["boundary"], solved)
            tangents[level], defects[level], hashes[level] = tangent, symplectic_residual(tangent), array_sha256(tangent)
            if hostile_kno and level == "H12":
                bad_kno_defect = symplectic_residual(tangent_from_solution(hessian, stage["boundary"], solved, True))[0]
            del rhs, solved
    tangent_variation, epsilon_symplectic, canonical = math.inf, math.inf, False
    if regular:
        tangent_variation = family_variation(tangents)
        norm_t = max(1.0, *(float(la.norm(value, 2)) for value in tangents.values()))
        epsilon_symplectic = (2 * norm_t * tangent_variation + tangent_variation**2
                              + stage["epsilon_h"] * norm_t**2
                              + 100 * MACHINE_EPS * max(conditions.values()) * norm_t**2)
        canonical = all(value[0] <= 10 * epsilon_symplectic for value in defects.values())
    hostile_kno_ok = not hostile_kno or not regular or bad_kno_defect > 100 * epsilon_symplectic
    stage["analysis"] = {
        "regular": bool(regular), "canonical": bool(canonical), "normalization_j": normalization_j,
        "step_j": step_j, "driver_error": driver_error, "binary_error": binary_error,
        "epsilon_j": epsilon_j, "minima": minima, "conditions": conditions, "tangents": tangents,
        "defects": defects, "tangent_hashes": hashes, "tangent_variation": tangent_variation,
        "epsilon_symplectic": epsilon_symplectic, "bad_kno_defect": bad_kno_defect,
        "hostile_kno_ok": bool(hostile_kno_ok),
    }
    del j_matrices, singular, singular_gesdd
    gc.collect()
    return stage["analysis"]


def compare_tangent_families(
    left,
    right,
    left_transform=None,
    conditioning=1.0,
    extra_uncertainty=0.0,
):
    transformed = {key: left_transform(value) if left_transform else value for key, value in left.items()}
    normalization = max(1.0, float(la.norm(transformed["H12"], "fro")), float(la.norm(right["H12"], "fro")))
    distances = {key: float(la.norm(transformed[key] - right[key], "fro")) / normalization for key in LEVELS}
    uncertainty = (
        (family_variation(transformed) + family_variation(right))
        / normalization
        + float(extra_uncertainty)
        + 100 * MACHINE_EPS * max(1.0, float(conditioning))
    )
    return {"distances": distances, "uncertainty": uncertainty,
            "label": classify(max(distances.values()), uncertainty)}


def product_family(second, first):
    return {key: second[key] @ first[key] for key in LEVELS}


def product_conditioning(first_condition, second_condition):
    return (
        float(first_condition)
        + float(second_condition)
        + MACHINE_EPS * float(first_condition) * float(second_condition)
    )


def canonical_family(family, conditioning):
    defects = {key: symplectic_residual(value) for key, value in family.items()}
    variation = family_variation(family)
    norm_t = max(1.0, *(float(la.norm(value, 2)) for value in family.values()))
    uncertainty = 2 * norm_t * variation + variation**2 + 100 * MACHINE_EPS * max(1.0, conditioning) * norm_t**2
    return {"passed": bool(all(value[0] <= 10 * uncertainty for value in defects.values())),
            "defects": defects, "variation": variation, "uncertainty": uncertainty}


def make_sector_basis(index_data, boundary, sector, common_seeds):
    dimension, selected = sector["dimension"], sector["basis"]
    old_seeds = tuple(min(index_data["orbit_edges"][orbit_type]) for orbit_type in range(30))
    type_by_seed = {seed: orbit_type for orbit_type, seed in enumerate(old_seeds)}
    inverse_lex = {old_index: position for position, old_index in enumerate(boundary["lexicographic"])}
    q = np.zeros((OLD, 30 * dimension), dtype=np.complex128)
    for common_type, seed in enumerate(common_seeds):
        orbit_type = type_by_seed[seed]
        for group in range(24):
            dense_row = inverse_lex[24 * orbit_type + group]
            for component in range(dimension):
                q[dense_row, common_type * dimension + component] = complex(selected[group, component])
    phase = np.zeros((PHASE, 60 * dimension), dtype=np.complex128)
    phase[:OLD, :30 * dimension], phase[OLD:, 30 * dimension:] = q, q
    return phase


def project_family(family, basis):
    return {key: basis.conj().T @ value @ basis for key, value in family.items()}


def compare_primary_entry(projected, midpoint, radii, conditioning):
    middle = projected["H12"]
    normalization = max(1.0, float(np.max(np.abs(middle))), float(np.max(np.abs(midpoint))))
    distance = float(np.max(np.abs(middle - midpoint))) / normalization
    variation = max(
        float(np.max(np.abs(projected["H01"] - middle))),
        float(np.max(np.abs(projected["H23"] - middle))),
    ) / normalization
    primary_radius = float(np.max(radii)) / normalization
    binary = 100 * MACHINE_EPS * max(1.0, conditioning)
    uncertainty = variation + primary_radius + binary
    return {"label": classify(distance, uncertainty), "distance": distance, "uncertainty": uncertainty,
            "dense_variation": variation, "primary_radius": primary_radius, "binary": binary}


def public_stage(stage):
    analysis = stage["analysis"]
    return {
        "carrier_ok": stage["carrier_ok"], "branch_ok": stage["branch_ok"],
        "assembly_ok": stage["assembly_ok"], "reciprocity_ok": stage["reciprocity_ok"],
        "seconds": float_text(stage["seconds"]),
        "hessian": {
            "normalization": float_text(stage["normalization_h"]), "step": float_text(stage["step_h"]),
            "round": float_text(stage["round_h"]), "uncertainty": float_text(stage["epsilon_h"]),
            "raw_reciprocity": {key: float_text(value) for key, value in stage["raw_reciprocity"].items()},
            "richardson_reciprocity": {key: float_text(value) for key, value in stage["richardson_reciprocity"].items()},
            "richardson_hashes": {key: array_sha256(value) for key, value in stage["levels"].items()},
        },
        "pre_legendre": {
            "classification": "REGULAR" if analysis["regular"] else "NUMERICALLY_OPEN",
            "normalized_smallest_singular_values": {key: float_text(value) for key, value in analysis["minima"].items()},
            "gap_gate": float_text(100 * analysis["epsilon_j"]),
            "condition_estimates": {key: float_text(value) for key, value in analysis["conditions"].items()},
            "gesvd_gesdd_error": float_text(analysis["driver_error"]),
        },
        "tangent": {
            "canonical": analysis["canonical"], "hashes": analysis["tangent_hashes"],
            "variation": float_text(analysis["tangent_variation"]),
            "symplectic_uncertainty": float_text(analysis["epsilon_symplectic"]),
            "symplectic_defects": {key: {"total": float_text(value[0]),
                "components": {name: float_text(component) for name, component in value[1].items()}}
                for key, value in analysis["defects"].items()},
            "bad_KNO_defect": None if analysis["bad_kno_defect"] is None else float_text(analysis["bad_kno_defect"]),
        },
    }


print("=" * 78)
print("DENSE ADVERSARIAL SECOND FINITE-HEIGHT TANGENT AND TWO-STEP MAP")
print("=" * 78)
input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
provenance_ok = input_hashes == EXPECTED_HASHES
check("all frozen protocol, artifacts and audited sources retain hashes", provenance_ok, str(input_hashes))
registered, duplicates = registry_inventory(RUN_ALL)
registry_ok = registered.count(VERIFIER_NAME) == 1 and not duplicates
check("the adversarial verifier is registered once with no duplicates", registry_ok,
      f"entries={len(registered)}, duplicates={duplicates}")

primary_header = json.loads(PRIMARY_JSON.read_text())
first_dense_header = json.loads(FIRST_DENSE_JSON.read_text())
composition = json.loads(COMPOSITION_ADVERSARIAL.read_text())
accepted_inputs_ok = bool(
    primary_header["outcome"] == "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY"
    and primary_header["tests"] == primary_header["passed"] == 31
    and primary_header["numeric_archive"]["array_count"] == 672
    and primary_header["numeric_archive"]["sha256"] == EXPECTED_HASHES["primary_numeric"]
    and first_dense_header["outcome"] == "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED"
    and first_dense_header["tests"] == first_dense_header["passed"] == 22
    and composition["outcome"] == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED"
    and composition["tests"] == composition["passed"] == 9
)
check("the primary header and accepted dense/history controls retain outcomes", accepted_inputs_ok)
history = reconstruct_history(composition)
first_step, second_step = history["first"], history["second"]
history_ok = bool(
    first_step["bracket_ok"] and second_step["bracket_ok"]
    and abs(first_step["residual"]) < mp.mpf("1e-100") and abs(second_step["residual"]) < mp.mpf("1e-100")
    and first_step["width"] < mp.mpf("1e-100") and second_step["width"] < mp.mpf("1e-100")
    and history["committed_error"] < mp.mpf("1e-65") and history["junction_error"] < mp.mpf("1e-100")
    and min(first_step["q"], first_step["h"], first_step["r"], first_step["m_out"],
            second_step["q"], second_step["h"], second_step["r"], second_step["m_out"],
            second_step["r"] - second_step["h"]**2) > 0
)
check("independent bisection reconstructs the frozen 180-digit branch-B history", history_ok,
      f"q1={mp_text(first_step['q'],25)}, q2={mp_text(second_step['q'],25)}, error={mp_text(history['committed_error'],6)}")

spec = importlib.util.spec_from_file_location("global_regge_orbits_dense_second_tangent", GEOMETRY_SOURCE)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
geometry_ok = gro.tests == gro.passed == 43
check("the complete geometric carrier retains all 43 certificates", geometry_ok)

r1, r2, h1, h2 = first_step["r"], second_step["r"], first_step["h"], second_step["h"]
m0, m1, c = history["m0"], first_step["m_out"], r1**2
first_library = make_dense_library(gro, m0, h1**2, mp.mpf(1))
second_normalized_library = make_dense_library(gro, m1, h2**2, mp.mpf(1))
second_physical_library = make_dense_library(gro, m0, c * h2**2, c)
models = {parity: first_library["augment_boundary_orbits"](model) for parity, model in gro.models.items()}

stages = {parity: {} for parity in ("even", "odd")}
scale_hessian = {}
for parity in ("even", "odd"):
    model = models[parity]
    stages[parity]["first_physical"] = assemble_stage("first_physical", parity, first_library, model, r1)
    del stages[parity]["first_physical"]["raw"]
    gc.collect()
    analyze_stage(stages[parity]["first_physical"])
    stages[parity]["second_normalized"] = assemble_stage("second_normalized", parity, second_normalized_library, model, r2)
    stages[parity]["second_physical"] = assemble_stage("second_physical", parity, second_physical_library, model, r2)
    records, label = scale_hessian_records(stages[parity]["second_normalized"], stages[parity]["second_physical"], c)
    scale_hessian[parity] = {"records": records, "label": label}
    check(f"{parity}: all three dense carriers, branches and raw reciprocity gates pass",
          all(item["carrier_ok"] and item["branch_ok"] and item["assembly_ok"] and item["reciprocity_ok"]
              for item in stages[parity].values()))
    check(f"{parity}: seven raw/Richardson Hessian scale comparisons receive labels",
          len(records) == 7 and all(item["label"] in {"AGREES", "REFUTED", "OPEN"} for item in records),
          f"outcome={label}, labels={dict(Counter(x['label'] for x in records))}")
    for stage_name in ("second_normalized", "second_physical"):
        del stages[parity][stage_name]["raw"]
    gc.collect()
    for stage_name in ("second_normalized", "second_physical"):
        item = stages[parity][stage_name]
        analyze_stage(item, hostile_kno=(parity == "even" and stage_name == "second_physical"))
    for stage_name, item in stages[parity].items():
        check(f"{parity}/{stage_name}: the full J rank and canonical map receive labels",
              isinstance(item["analysis"]["regular"], bool)
              and isinstance(item["analysis"]["canonical"], bool),
              f"minima={item['analysis']['minima']}, gate={100*item['analysis']['epsilon_j']:.3e}")

first_hash_control = {}
for parity in ("even", "odd"):
    first_hash_control[parity] = (
        stages[parity]["first_physical"]["analysis"]["tangent_hashes"]
        == first_dense_header["parities"][parity]["tangent_hashes"]
    )
check("the independently recomputed first slab reproduces all six accepted dense hashes",
      all(first_hash_control.values()), str(first_hash_control))
scale_labels = {item["label"] for item in scale_hessian.values()}
scale_lift_outcome = "SCALE_LIFT_REFUTED" if "SCALE_LIFT_REFUTED" in scale_labels else "SCALE_LIFT_OPEN" if "SCALE_LIFT_OPEN" in scale_labels else "SCALE_LIFT_CONFIRMED"
check("both complete real Hessian families freeze one scale-lift outcome",
      scale_lift_outcome in {"SCALE_LIFT_CONFIRMED", "SCALE_LIFT_REFUTED", "SCALE_LIFT_OPEN"}, scale_lift_outcome)

maps_available = all(
    item["analysis"]["regular"]
    for parity in stages.values() for item in parity.values()
)
tangent_scale, identity_hostile, linear_hostile = {}, {}, {}
if maps_available:
    for parity in ("even", "odd"):
        normalized = stages[parity]["second_normalized"]["analysis"]["tangents"]
        physical = stages[parity]["second_physical"]["analysis"]["tangents"]
        tangent_condition = max(
            *stages[parity]["second_normalized"]["analysis"]["conditions"].values(),
            *stages[parity]["second_physical"]["analysis"]["conditions"].values(),
        )
        tangent_extra = (
            stages[parity]["second_normalized"]["epsilon_h"]
            + stages[parity]["second_physical"]["epsilon_h"]
        )
        tangent_scale[parity] = compare_tangent_families(
            normalized, physical, lambda value: scale_tangent(value, c),
            tangent_condition, tangent_extra,
        )
        identity_hostile[parity] = compare_tangent_families(
            normalized, physical, conditioning=tangent_condition,
            extra_uncertainty=tangent_extra,
        )
        linear_hostile[parity] = compare_tangent_families(
            normalized, physical, lambda value: scale_tangent(value, r1),
            tangent_condition, tangent_extra,
        )
else:
    for parity in ("even", "odd"):
        empty = {"label": "NOT_EVALUATED", "distances": {}, "uncertainty": math.inf}
        tangent_scale[parity] = dict(empty)
        identity_hostile[parity] = dict(empty)
        linear_hostile[parity] = dict(empty)
check("the direct physical c=r1^2 tangent comparison receives both labels",
      all(item["label"] in {"AGREES", "REFUTED", "OPEN", "NOT_EVALUATED"}
          for item in tangent_scale.values()),
      str({key: value["label"] for key, value in tangent_scale.items()}))
hostile_tangents_ok = bool(
    not maps_available
    or (
        any(item["label"] == "REFUTED" for item in identity_hostile.values())
        and any(item["label"] == "REFUTED" for item in linear_hostile.values())
        and stages["even"]["second_physical"]["analysis"]["hostile_kno_ok"]
    )
)
check("identity, r1 and omitted-K_NO hostile tangent constructions are rejected",
      hostile_tangents_ok)

if maps_available:
    second_schedule_condition = max(
        *stages["even"]["second_physical"]["analysis"]["conditions"].values(),
        *stages["odd"]["second_physical"]["analysis"]["conditions"].values(),
    )
    second_schedule_extra = (
        stages["even"]["second_physical"]["epsilon_h"]
        + stages["odd"]["second_physical"]["epsilon_h"]
    )
    second_schedule = compare_tangent_families(
        stages["even"]["second_physical"]["analysis"]["tangents"],
        stages["odd"]["second_physical"]["analysis"]["tangents"],
        conditioning=second_schedule_condition,
        extra_uncertainty=second_schedule_extra,
    )
else:
    second_schedule = {
        "label": "NOT_EVALUATED", "distances": {}, "uncertainty": math.inf,
    }
check("the two direct physical second-slab schedules receive the dense classifier",
      second_schedule["label"] in {"AGREES", "REFUTED", "OPEN", "NOT_EVALUATED"},
      f"label={second_schedule['label']}, distances={second_schedule['distances']}")

products, product_canonical, product_conditions = {}, {}, {}
if maps_available:
    for p1 in ("even", "odd"):
        for p2 in ("even", "odd"):
            key = f"{p1}_{p2}"
            products[key] = product_family(
                stages[p2]["second_physical"]["analysis"]["tangents"],
                stages[p1]["first_physical"]["analysis"]["tangents"],
            )
            conditioning = product_conditioning(
                max(stages[p1]["first_physical"]["analysis"]["conditions"].values()),
                max(stages[p2]["second_physical"]["analysis"]["conditions"].values()),
            )
            product_conditions[key] = conditioning
            product_canonical[key] = canonical_family(products[key], conditioning)
check("all four dense products receive the real symplectic block classifier",
      (not maps_available)
      or (len(product_canonical) == 4
          and all(isinstance(item["passed"], bool) for item in product_canonical.values())))
product_schedule = []
for left, right in combinations(tuple(products), 2):
    product_schedule.append({"left": left, "right": right,
                             **compare_tangent_families(
                                 products[left], products[right],
                                 conditioning=max(product_conditions[left], product_conditions[right]),
                             )})
product_labels = {item["label"] for item in product_schedule}
product_schedule_outcome = (
    "NOT_EVALUATED" if not maps_available else
    "TWO_STEP_SCHEDULE_DEPENDENT" if "REFUTED" in product_labels else
    "TWO_STEP_SCHEDULE_OPEN" if "OPEN" in product_labels else
    "TWO_STEP_SCHEDULE_ROBUST"
)
check("all six dense pairs among the four two-step products receive labels",
      (not maps_available) or len(product_schedule) == 6,
      f"outcome={product_schedule_outcome}, labels={dict(Counter(x['label'] for x in product_schedule))}")

if maps_available:
    synthetic = {key: value.copy() for key, value in products["even_even"].items()}
    for value in synthetic.values():
        value[0, 0] += 1e-3
    synthetic_record = compare_tangent_families(
        products["even_even"], synthetic,
        conditioning=product_conditions["even_even"],
    )
    first_even = stages["even"]["first_physical"]["analysis"]["tangents"]
    normalized_even = stages["even"]["second_normalized"]["analysis"]["tangents"]
    wrong_identity_product = product_family(normalized_even, first_even)
    wrong_linear_product = product_family(
        {key: scale_tangent(value, r1) for key, value in normalized_even.items()},
        first_even,
    )
    identity_product = compare_tangent_families(
        products["even_even"], wrong_identity_product,
        conditioning=product_conditions["even_even"],
    )
    linear_product = compare_tangent_families(
        products["even_even"], wrong_linear_product,
        conditioning=product_conditions["even_even"],
    )
    hostile_products_ok = bool(
        synthetic_record["label"] == "REFUTED"
        and identity_product["label"] == "REFUTED"
        and linear_product["label"] == "REFUTED"
    )
else:
    synthetic_record = identity_product = linear_product = {
        "label": "NOT_EVALUATED", "distances": {}, "uncertainty": math.inf,
    }
    hostile_products_ok = True
check("the 1e-3 corruption and both wrong product lifts are detected", hostile_products_ok,
      f"synthetic={synthetic_record['label']}, identity={identity_product['label']}, linear={linear_product['label']}")

dense_labels = {
    "scale_lift": scale_lift_outcome,
    "all_regular": all(item["analysis"]["regular"] for parity in stages.values() for item in parity.values()),
    "all_canonical": all(item["analysis"]["canonical"] for parity in stages.values() for item in parity.values()),
    "correct_tangent_scale": all(item["label"] == "AGREES" for item in tangent_scale.values()),
    "second_schedule": second_schedule["label"], "product_schedule": product_schedule_outcome,
    "products_canonical": bool(
        maps_available and all(item["passed"] for item in product_canonical.values())
    ),
}
dense_labels_frozen = True
print("[closure] dense labels frozen; opening primary bases and tangent entries", flush=True)
basis_library = {
    "I": I, "math": math, "mp": mp, "np": np,
    "mp_frobenius": lambda matrix: mp.sqrt(mp.fsum(abs(value) ** 2 for value in matrix)),
    "mp_submatrix": lambda matrix, rows, columns: mp.matrix([[matrix[row, column] for column in columns] for row in rows]),
}
load_named_functions(OLD_TANGENT_SOURCE, {"cluster_sorted", "high_precision_sector_bases"}, basis_library)
primary_numeric = np.load(PRIMARY_NUMERIC, allow_pickle=False)
first_primary_numeric = np.load(FIRST_PRIMARY_NUMERIC, allow_pickle=False)
archive_schema_ok = len(primary_numeric.files) == 672 and len(first_primary_numeric.files) == 168 and dense_labels_frozen
check("both primary archives open only after dense classification with exact array counts", archive_schema_ok,
      f"second={len(primary_numeric.files)}, first={len(first_primary_numeric.files)}")
reference_index = stages["even"]["second_normalized"]["index_data"]
reference_boundary = stages["even"]["second_normalized"]["boundary"]
with mp.workdps(180):
    sectors, sector_control = basis_library["high_precision_sector_bases"](
        reference_index
    )
common_seeds = tuple(sorted(min(reference_index["orbit_edges"][orbit_type]) for orbit_type in range(30)))
basis_ok = (sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
            and sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
            and sum(60 * sector["dimension"]**2 for sector in sectors) == PHASE)
check("the delayed deterministic minimal bases exhaust the 1440-dimensional phase carrier", basis_ok,
      f"dimensions={[sector['dimension'] for sector in sectors]}")

primary_comparisons = {"first": [], "second": [], "products": []}
if maps_available:
    for sector_index, sector in enumerate(sectors):
        basis = make_sector_basis(reference_index, reference_boundary, sector, common_seeds)
        dimension = sector["dimension"]
        for parity in ("even", "odd"):
            projected_first = project_family(
                stages[parity]["first_physical"]["analysis"]["tangents"], basis
            )
            first_key = f"{parity}_sector{sector_index}_K12"
            record = compare_primary_entry(
                projected_first,
                first_primary_numeric[f"{first_key}_tangent_midpoint"],
                first_primary_numeric[f"{first_key}_tangent_radii"],
                stages[parity]["first_physical"]["analysis"]["conditions"]["H12"],
            )
            primary_comparisons["first"].append({
                "parity": parity, "sector": sector_index,
                "dimension": dimension, **record,
            })
            projected_second = project_family(
                stages[parity]["second_physical"]["analysis"]["tangents"], basis
            )
            second_key = f"{parity}_sector{sector_index}_physical_K12"
            record = compare_primary_entry(
                projected_second,
                primary_numeric[f"{second_key}_tangent_midpoint"],
                primary_numeric[f"{second_key}_tangent_radii"],
                stages[parity]["second_physical"]["analysis"]["conditions"]["H12"],
            )
            primary_comparisons["second"].append({
                "parity": parity, "sector": sector_index,
                "dimension": dimension, **record,
            })
        for pair, family in products.items():
            projected_product = project_family(family, basis)
            product_key = f"product_{pair}_sector{sector_index}_K12"
            p1, p2 = pair.split("_")
            conditioning = product_conditioning(
                stages[p1]["first_physical"]["analysis"]["conditions"]["H12"],
                stages[p2]["second_physical"]["analysis"]["conditions"]["H12"],
            )
            record = compare_primary_entry(
                projected_product,
                primary_numeric[f"{product_key}_tangent_midpoint"],
                primary_numeric[f"{product_key}_tangent_radii"],
                conditioning,
            )
            primary_comparisons["products"].append({
                "pair": pair, "sector": sector_index,
                "dimension": dimension, **record,
            })
        del basis
comparison_counts_ok = bool(
    not maps_available
    or (
        len(primary_comparisons["first"]) == 14
        and len(primary_comparisons["second"]) == 14
        and len(primary_comparisons["products"]) == 28
    )
)
all_primary_labels = [item["label"] for family in primary_comparisons.values() for item in family]
check("delayed entrywise reconciliation covers 14 first, 14 second and 28 product blocks", comparison_counts_ok,
      str({key: len(value) for key, value in primary_comparisons.items()}))
check("every delayed entrywise primary comparison receives the frozen classifier",
      (not maps_available)
      or all(label in {"AGREES", "REFUTED", "OPEN"} for label in all_primary_labels),
      str(dict(Counter(all_primary_labels))))

base_controls_ok = bool(
    provenance_ok and registry_ok and accepted_inputs_ok and history_ok and geometry_ok
    and all(first_hash_control.values()) and archive_schema_ok and basis_ok and comparison_counts_ok
    and hostile_tangents_ok and hostile_products_ok
    and all(item["carrier_ok"] and item["branch_ok"] and item["assembly_ok"] and item["reciprocity_ok"]
            and item["analysis"]["hostile_kno_ok"] for parity in stages.values() for item in parity.values())
)
if not base_controls_ok:
    outcome = "SECOND_FULL_TANGENT_DENSE_CONTROL_FAILED"
elif scale_lift_outcome == "SCALE_LIFT_REFUTED":
    outcome = "SECOND_FULL_TANGENT_DENSE_SCALE_LIFT_REFUTED"
elif scale_lift_outcome == "SCALE_LIFT_OPEN":
    outcome = "SECOND_FULL_TANGENT_DENSE_SCALE_LIFT_OPEN"
elif not dense_labels["all_regular"]:
    outcome = "SECOND_FULL_TANGENT_DENSE_RANK_OPEN"
elif not (dense_labels["all_canonical"] and dense_labels["correct_tangent_scale"] and dense_labels["products_canonical"]):
    outcome = "SECOND_FULL_TANGENT_DENSE_CANONICALITY_FAILED"
elif dense_labels["second_schedule"] == "REFUTED" or dense_labels["product_schedule"] == "TWO_STEP_SCHEDULE_DEPENDENT":
    outcome = "SECOND_FULL_TANGENT_DENSE_SCHEDULE_DEPENDENT"
elif dense_labels["second_schedule"] == "OPEN" or dense_labels["product_schedule"] == "TWO_STEP_SCHEDULE_OPEN":
    outcome = "SECOND_FULL_TANGENT_DENSE_OPEN"
elif "REFUTED" in all_primary_labels:
    outcome = "SECOND_FULL_TANGENT_PRIMARY_REFUTED"
elif "OPEN" in all_primary_labels:
    outcome = "SECOND_FULL_TANGENT_DENSE_OPEN"
else:
    outcome = "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED"
allowed_outcomes = {
    "SECOND_FULL_TANGENT_DENSE_CONTROL_FAILED", "SECOND_FULL_TANGENT_DENSE_SCALE_LIFT_REFUTED",
    "SECOND_FULL_TANGENT_DENSE_SCALE_LIFT_OPEN", "SECOND_FULL_TANGENT_DENSE_RANK_OPEN",
    "SECOND_FULL_TANGENT_DENSE_CANONICALITY_FAILED", "SECOND_FULL_TANGENT_DENSE_SCHEDULE_DEPENDENT",
    "SECOND_FULL_TANGENT_DENSE_OPEN", "SECOND_FULL_TANGENT_PRIMARY_REFUTED",
    "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED",
}
check("the dense adversarial result follows the preregistered nine-level hierarchy", outcome in allowed_outcomes, outcome)


def public_record(record):
    return {key: float_text(value) if isinstance(value, float) else value for key, value in record.items()}


artifact = {
    "outcome": outcome, "tests": tests, "passed": passed,
    "method": "six direct 2280-edge Hessian assemblies, eighteen full pre-Legendre SVD pairs and 1440-column solves, dense real symplectic identities and delayed entrywise closure",
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "first_run_correction_commit": "88e3ef3",
        "input_sha256": input_hashes,
    },
    "history": {
        "v": mp_text(history["v"]), "m0": mp_text(history["m0"]), "pi0": mp_text(history["pi0"]),
        "q1": mp_text(first_step["q"]), "h1": mp_text(first_step["h"]), "r1": mp_text(first_step["r"]),
        "m1": mp_text(first_step["m_out"]), "pi1": mp_text(first_step["pi_out"]),
        "q2": mp_text(second_step["q"]), "h2": mp_text(second_step["h"]), "r2": mp_text(second_step["r"]),
        "m2": mp_text(second_step["m_out"]), "pi2": mp_text(second_step["pi_out"]),
        "committed_error": mp_text(history["committed_error"]), "junction_error": mp_text(history["junction_error"]),
    },
    "scale_factor_r1_squared": mp_text(c), "dense_labels_frozen_before_primary_entries": dense_labels,
    "hessian_scale_lift": {parity: {"label": item["label"], "records": [public_record(record) for record in item["records"]]}
                           for parity, item in scale_hessian.items()},
    "tangent_scale_lift": {parity: public_record(record) for parity, record in tangent_scale.items()},
    "second_schedule": public_record(second_schedule), "product_schedule_outcome": product_schedule_outcome,
    "product_schedule": [public_record(record) for record in product_schedule],
    "product_canonicality": {key: {
        "passed": value["passed"], "variation": float_text(value["variation"]),
        "uncertainty": float_text(value["uncertainty"]),
        "defects": {level: float_text(record[0]) for level, record in value["defects"].items()},
        "hashes": {level: array_sha256(products[key][level]) for level in LEVELS},
    } for key, value in product_canonical.items()},
    "hostile_controls": {
        "identity_tangent": {key: public_record(value) for key, value in identity_hostile.items()},
        "linear_r1_tangent": {key: public_record(value) for key, value in linear_hostile.items()},
        "synthetic_product": public_record(synthetic_record), "identity_product": public_record(identity_product),
        "linear_r1_product": public_record(linear_product),
    },
    "first_dense_hash_control": first_hash_control,
    "stages": {parity: {name: public_stage(item) for name, item in family.items()} for parity, family in stages.items()},
    "delayed_primary_reconciliation": {
        "basis_decimal_digits": 180,
        "product_binary_conditioning": "kappa_1+kappa_2+u*kappa_1*kappa_2",
        "comparisons": {
            key: [public_record(record) for record in value]
            for key, value in primary_comparisons.items()
        },
    },
    "firewall": {
        "primary_entries_opened_after_dense_labels_frozen": dense_labels_frozen,
        "primary_minimal_basis_used_in_dense_decision": False, "orbit_convolution_used_in_dense_decision": False,
        "tangent_eigenvalues_computed": False, "tangent_singular_values_computed": False,
        "continuum_or_particle_target_parsed": False, "full_suite_run": False,
    },
    "classification": {
        "two_step_finite_height_linear_response": "DERIVED_COMPUTATIONAL_ADVERSARIALLY_REPLICATED" if outcome.endswith("ADVERSARIALLY_REPLICATED") else "OPEN_OR_REFUTED_BY_OUTCOME",
        "physical_constraint_quotient": "OPEN", "physical_mode_spectrum": "NOT_COMPUTED",
        "wave_equation_limiting_speed_G_planck_particles": "NOT_DERIVED", "external_novelty": "OPEN",
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
print("-" * 78)
print(f"OUTCOME: {outcome}")
print(f"HESSIAN SCALE: {scale_lift_outcome}")
print(f"SECOND SCHEDULE: {second_schedule['label']}")
print(f"PRODUCT SCHEDULE: {product_schedule_outcome}")
print(f"PRIMARY LABELS: {dict(Counter(all_primary_labels))}")
print(f"RESULT: {passed}/{tests} PASS")
if passed != tests:
    raise SystemExit(1)
