#!/usr/bin/env python3
"""Primary second finite-height full tangent and physical two-step map.

Prior-art commit: eae691a.
Protocol commit: f80e174.
Pre-evaluation repository-prior-art correction: ddb3f6e.
Pre-evaluation history-provenance correction: 8f096be.
Preserved first control failure: 38f59f5.
Preserved converter failure: b263afa.
Registry commit: f0b8177.

No eigenvalue spectrum, tangent singular-value census, continuum harmonic,
speed, Planck or particle target is evaluated.
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
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_finite_height_second_full_boundary_tangent.json"
NUMERIC_OUTPUT = (
    HERE / "gravity_600cell_finite_height_second_full_boundary_tangent.npz"
)
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_second_full_boundary_tangent_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_second_full_boundary_tangent_protocol.md"
)
PROTOCOL_CORRECTION = (
    ROOT
    / "docs/gravity/gravity_600cell_second_full_boundary_tangent_protocol_correction.md"
)
HISTORY_CORRECTION = (
    ROOT
    / "docs/gravity/gravity_600cell_second_full_boundary_tangent_history_provenance_correction.md"
)
FIRST_FAILURE = (
    ROOT
    / "docs/gravity/gravity_600cell_second_full_boundary_tangent_first_run_failure.md"
)
SECOND_FAILURE = (
    ROOT
    / "docs/gravity/gravity_600cell_second_full_boundary_tangent_second_run_failure.md"
)
HISTORY_INPUT = HERE / "gravity_600cell_finite_height_asymptotic_map.json"
HISTORY_SOURCE = HERE / "verify_gravity_600cell_finite_height_asymptotic_map.py"
FIRST_INPUT = HERE / "gravity_600cell_finite_height_full_boundary_tangent.json"
FIRST_NUMERIC = HERE / "gravity_600cell_finite_height_full_boundary_tangent.npz"
FIRST_SOURCE = HERE / "verify_gravity_600cell_finite_height_full_boundary_tangent.py"
FIRST_ADVERSARIAL = (
    HERE / "gravity_600cell_finite_height_full_boundary_tangent_adversarial.json"
)
COMPOSITION_INPUT = HERE / "gravity_600cell_finite_height_composition.json"
COMPOSITION_SOURCE = HERE / "verify_gravity_600cell_finite_height_composition.py"
COMPOSITION_ADVERSARIAL = (
    HERE / "gravity_600cell_finite_height_composition_adversarial.json"
)
COMPOSITION_ADVERSARIAL_SOURCE = (
    HERE / "verify_gravity_600cell_finite_height_composition_adversarial.py"
)
OLD_TWO_PRIOR = (
    ROOT / "docs/gravity/gravity_600cell_dust_two_step_full_tangent_prior_art.md"
)
OLD_TWO_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_dust_two_step_full_tangent_protocol.md"
)
OLD_TWO_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_dust_two_step_full_tangent_result.md"
)
OLD_TWO_SOURCE = HERE / "verify_gravity_600cell_dust_two_step_full_tangent.py"
OLD_TWO_INPUT = HERE / "gravity_600cell_dust_two_step_full_tangent.json"
OLD_TWO_NUMERIC = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"
OLD_TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "protocol_correction": PROTOCOL_CORRECTION,
    "history_correction": HISTORY_CORRECTION,
    "first_failure": FIRST_FAILURE,
    "second_failure": SECOND_FAILURE,
    "history": HISTORY_INPUT,
    "history_source": HISTORY_SOURCE,
    "first": FIRST_INPUT,
    "first_numeric": FIRST_NUMERIC,
    "first_source": FIRST_SOURCE,
    "first_adversarial": FIRST_ADVERSARIAL,
    "composition": COMPOSITION_INPUT,
    "composition_source": COMPOSITION_SOURCE,
    "composition_adversarial": COMPOSITION_ADVERSARIAL,
    "composition_adversarial_source": COMPOSITION_ADVERSARIAL_SOURCE,
    "old_two_prior": OLD_TWO_PRIOR,
    "old_two_protocol": OLD_TWO_PROTOCOL,
    "old_two_result": OLD_TWO_RESULT,
    "old_two_source": OLD_TWO_SOURCE,
    "old_two": OLD_TWO_INPUT,
    "old_two_numeric": OLD_TWO_NUMERIC,
    "old_tangent_source": OLD_TANGENT_SOURCE,
    "rank_source": RANK_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
}
EXPECTED_HASHES = {
    "prior_art": "d3740e0b08b2f3ec6adf2c69c762e5e5dc0cdd87a571d6d27bc62e78518e70be",
    "protocol": "c6615e5011f0a07e5ddaccd00c63d7ed9419e4058451afc3e2243423098c7024",
    "protocol_correction": "cd8b99a64e3c88c91781188db0e08e5f8bf2cfc1d450e305949f7350288cda80",
    "history_correction": "fc321949e07de6ec551743abf5274c337c9ed0adb8628142fa1e83f218f1164b",
    "first_failure": "d9cdc909aeddb8a01b0051420c2dead72007905d3e8420fe259b843d1aefd0ef",
    "second_failure": "58eda42f4babb4ea3c47d3ade1d5e5d327595ce1edb55127c4e488229bf804a6",
    "history": "a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6",
    "history_source": "3aafdb326eb9299d9e69ef79c0726eeb09f214b9dee1dc848848e34e0920b208",
    "first": "266638aeaa825b327b63a84eda36a499456dc4b4f9a86f964cee5f79d6d6e930",
    "first_numeric": "0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf",
    "first_source": "c4e60d6ef87131d87a93b64d5381d16d8de8d3990340efd5405ec983f64db94d",
    "first_adversarial": "ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7",
    "composition": "d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447",
    "composition_source": "cb4cf619dc54922d3a64d5e000a6cd0d3c19f71cda2a18b66d69f68271496422",
    "composition_adversarial": "d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9",
    "composition_adversarial_source": "8395e921ab1c1f518abb567a114f1eb8bfdf2068be031bff55c8d2f0cff56c2b",
    "old_two_prior": "e7d865b9e72a411eee61e0ce091cde0d912fd9d9f773708f00a9c7046a6785f9",
    "old_two_protocol": "d5dd44ece724b65351b35fd18e6d334dbf4b68e9f2757484b76fdaf6c42fe0cf",
    "old_two_result": "014a86460433e9e8ab72a2aae029bed774306c2095aef1a543522d780c783038",
    "old_two_source": "c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717",
    "old_two": "f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc",
    "old_two_numeric": "ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d",
    "old_tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

PRIOR_ART_COMMIT = "eae691a"
PROTOCOL_COMMIT = "f80e174"
CORRECTION_COMMIT = "ddb3f6e"
HISTORY_CORRECTION_COMMIT = "8f096be"
FIRST_FAILURE_COMMIT = "38f59f5"
SECOND_FAILURE_COMMIT = "b263afa"
REGISTRY_COMMIT = "f0b8177"
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


def load_named_functions(path, names, namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [
        node
        for node in tree.body
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


def bisect_root(function, left, right):
    left = mp.mpf(left)
    right = mp.mpf(right)
    left_value = function(left)
    right_value = function(right)
    bracket_ok = left_value * right_value < 0
    for _ in range(700):
        middle = (left + right) / 2
        middle_value = function(middle)
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
    root = (left + right) / 2
    return root, right - left, bracket_ok


def reconstruct_history(first_background, composition, adversarial):
    v = mp.mpf(3) / 2
    m0 = mu(v)
    pi0 = momentum(v)

    def advance(mass, incoming, bracket):
        function = lambda q: (
            4 * mp.pi * (mu(q) - mass)
            + q * (momentum(q) - incoming)
        )
        q, width, bracket_ok = bisect_root(function, *bracket)
        h = (momentum(q) - incoming) / (2 * mp.pi * mu(q))
        ratio = 1 + h * q
        next_mass = mass / ratio
        next_momentum = momentum(q) + 2 * mp.pi * h * mu(q) / ratio
        return {
            "q": q,
            "h": h,
            "r": ratio,
            "m_in": mass,
            "pi_in": incoming,
            "m_out": next_mass,
            "pi_out": next_momentum,
            "residual": function(q),
            "width": width,
            "bracket_ok": bracket_ok,
        }

    first = advance(m0, pi0, (9, 10))
    second = advance(first["m_out"], first["pi_out"], (31, 32))
    primary_state = next(
        item for item in composition["composition"]
        if mp.mpf(item["v"]) == v
    )
    primary_branch = next(
        root for root in primary_state["roots"]
        if root["physical"] and mp.mpf(31) < mp.mpf(root["q2"]) < mp.mpf(32)
    )
    adversarial_run = next(
        item for item in adversarial["precision_runs"]
        if item["target_precision"] == 180
    )
    adversarial_branch = adversarial_run["branches"]["B"]
    primary_serialized_matches = {
        "q2": mp.nstr(second["q"], 60) == primary_branch["q2"],
        "h2": mp.nstr(second["h"], 60) == primary_branch["h2"],
        "r2": mp.nstr(second["r"], 60) == primary_branch["scale_ratio"],
    }
    committed_error = max(
        abs(first["q"] - mp.mpf(first_background["q"])),
        abs(first["h"] - mp.mpf(first_background["h"])),
        abs(first["r"] - mp.mpf(first_background["lambda"])),
        abs(first["q"] - mp.mpf(adversarial_run["q1"])),
        abs(first["h"] - mp.mpf(adversarial_run["h1"])),
        abs(first["r"] - mp.mpf(adversarial_run["L1"])),
        abs(first["m_out"] - mp.mpf(adversarial_run["m1"])),
        abs(first["pi_out"] - mp.mpf(adversarial_run["pi1"])),
        abs(second["q"] - mp.mpf(adversarial_branch["q2"])),
        abs(second["h"] - mp.mpf(adversarial_branch["h2"])),
        abs(second["r"] - mp.mpf(adversarial_branch["ratio"])),
    )
    junction_error = max(
        abs(first["m_out"] - second["m_in"]),
        abs(first["pi_out"] - second["pi_in"]),
        abs(m0 - first["r"] * first["m_out"]),
    )
    return {
        "v": v,
        "m0": m0,
        "pi0": pi0,
        "first": first,
        "second": second,
        "committed_error": committed_error,
        "junction_error": junction_error,
        "primary_branch_q2": mp.mpf(primary_branch["q2"]),
        "adversarial_branch_q2": mp.mpf(adversarial_branch["q2"]),
        "primary_serialized_matches": primary_serialized_matches,
    }


def exact_scale_control():
    c, a, b, cc, d = sp.symbols("c a b cc d", nonzero=True)
    omega = sp.Matrix([[0, 1], [-1, 0]])
    dilation = sp.diag(1, c)
    tangent = sp.Matrix([[a, b], [cc, d]])
    conjugated = dilation * tangent * dilation.inv()
    conformal = sp.simplify(dilation.T * omega * dilation - c * omega)
    conjugate_defect = (conjugated.T * omega * conjugated - omega).applyfunc(
        lambda value: sp.factor(value.subs(d, (1 + b * cc) / a))
    )
    s, beta = sp.symbols("s beta", nonzero=True)
    canonical_scale = sp.diag(s, 1 / s)
    shear = sp.Matrix([[1, 0], [beta, 1]])
    family_ok = bool(
        canonical_scale.T * omega * canonical_scale == omega
        and shear.T * omega * shear == omega
    )
    return {
        "passed": bool(
            conformal == sp.zeros(2)
            and conjugate_defect == sp.zeros(2)
            and family_ok
        ),
        "conformal_identity": str(dilation.T * omega * dilation),
        "conjugated": str(conjugated),
        "canonical_families": family_ok,
    }


def richardson(coarse, fine):
    keys = set(coarse) | set(fine)
    return {
        key: (4 * fine.get(key, 0) - coarse.get(key, 0)) / 3
        for key in keys
    }


def kernel_frobenius(kernel):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in kernel.values()))


def kernel_difference(left, right, right_scale=mp.mpf(1)):
    keys = set(left) | set(right)
    return mp.sqrt(mp.fsum(
        abs(left.get(key, 0) - right_scale * right.get(key, 0)) ** 2
        for key in keys
    ))


def scale_tangent_arrays(midpoint, radii, factor):
    size = midpoint.shape[0] // 2
    row = np.concatenate((np.ones(size), np.full(size, float(factor))))
    column = np.concatenate((np.ones(size), np.full(size, 1 / float(factor))))
    multiplier = row[:, np.newaxis] * column[np.newaxis, :]
    return midpoint * multiplier, radii * np.abs(multiplier)


def omega_array(size):
    omega = np.zeros((2 * size, 2 * size), dtype=float)
    omega[:size, size:] = np.eye(size)
    omega[size:, :size] = -np.eye(size)
    return omega


def product_with_radius(left, left_radius, right, right_radius):
    midpoint = left @ right
    radius = (
        np.abs(left) @ right_radius
        + left_radius @ np.abs(right)
        + left_radius @ right_radius
    )
    return midpoint, radius


def defect_with_radius(midpoint, radius):
    size = midpoint.shape[0] // 2
    omega = omega_array(size)
    defect = midpoint.conj().T @ omega @ midpoint - omega
    bound = (
        radius.T @ np.abs(omega @ midpoint)
        + np.abs(midpoint.conj().T @ omega) @ radius
        + radius.T @ np.abs(omega) @ radius
    )
    return defect, bound


def ball_records_from_built(built_by_level, converter):
    records = {}
    for level, built in built_by_level.items():
        midpoint, radii = converter(built["tangent"])
        defect_midpoint, defect_radii = converter(built["defect"])
        records[level] = {
            "midpoint": midpoint,
            "radii": radii,
            "defect_midpoint": defect_midpoint,
            "defect_radii": defect_radii,
        }
    return records


def comparison_record(left, right):
    normalization = max(
        1.0,
        float(la.norm(left["midpoints"]["K12"], "fro")),
        float(la.norm(right["midpoints"]["K12"], "fro")),
    )

    def within(item):
        middle = item["midpoints"]["K12"]
        return (
            la.norm(item["midpoints"]["K01"] - middle, "fro")
            + la.norm(middle - item["midpoints"]["K23"], "fro")
            + max(la.norm(value, "fro") for value in item["radii"].values())
            + 100 * np.finfo(float).eps * max(1.0, la.norm(middle, 2))
        )

    uncertainty = float(within(left) + within(right)) / normalization + 1e-135
    distances = {
        level: float(la.norm(
            left["midpoints"][level] - right["midpoints"][level], "fro"
        )) / normalization
        for level in LEVELS
    }
    maximum = max(distances.values())
    if maximum <= 10 * uncertainty:
        label = "AGREES"
    elif min(distances.values()) > 100 * uncertainty:
        label = "REFUTED"
    else:
        label = "OPEN"
    return {
        "label": label,
        "distances": distances,
        "uncertainty": uncertainty,
    }


def make_library(gro, mass, rho0, length_square):
    namespace = {
        "ARITHMETIC_FLOOR": ARITHMETIC_FLOOR,
        "Counter": Counter,
        "DPS": DPS,
        "DERIVATIVE_STEPS": DERIVATIVE_STEPS,
        "I": I,
        "LOCAL_EDGES": LOCAL_EDGES,
        "LOCAL_HINGES": LOCAL_HINGES,
        "LOCAL_HINGE_INDEX": LOCAL_HINGE_INDEX,
        "L0_SQUARE": mp.mpf(length_square),
        "MASS": mp.mpf(mass),
        "RHO0": mp.mpf(rho0),
        "VARIANTS": VARIANTS,
        "cluster_sorted": None,
        "combinations": combinations,
        "defaultdict": defaultdict,
        "gro": gro,
        "math": math,
        "mp": mp,
        "mp_frobenius": lambda matrix: mp.sqrt(
            mp.fsum(abs(value) ** 2 for value in matrix)
        ),
        "mp_submatrix": lambda matrix, rows, columns: mp.matrix([
            [matrix[row, column] for column in columns]
            for row in rows
        ]),
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
        namespace,
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
        namespace,
    )
    return namespace


def kernel_family_variation(family, order):
    return mp.fsum(
        kernel_difference(family[left], family[right])
        for left, right in zip(order, order[1:])
    )


def kernel_scale_records(physical, normalized, factor, order):
    physical_variation = kernel_family_variation(physical, order)
    normalized_variation = factor * kernel_family_variation(normalized, order)
    key_count = max(
        len(set(physical[level]) | set(normalized[level]))
        for level in order
    )
    flint_conversion = (
        mp.sqrt(max(1, key_count))
        * (ARITHMETIC_FLOOR + CLASSIFIER_FLOOR)
        * max(mp.mpf(1), abs(factor))
    )
    family_magnitude = max(
        mp.mpf(1),
        *(kernel_frobenius(physical[level]) for level in order),
        *(abs(factor) * kernel_frobenius(normalized[level]) for level in order),
    )
    binary_conversion = (
        100 * mp.mpf(np.finfo(float).eps) * family_magnitude
    )
    arithmetic = flint_conversion + binary_conversion
    records = {}
    for level in order:
        normalization = max(
            mp.mpf(1),
            kernel_frobenius(physical[level]),
            abs(factor) * kernel_frobenius(normalized[level]),
        )
        distance = kernel_difference(
            physical[level], normalized[level], right_scale=factor
        ) / normalization
        uncertainty = (
            physical_variation + normalized_variation + arithmetic
        ) / normalization
        if distance <= 10 * uncertainty:
            label = "AGREES"
        elif distance > 100 * uncertainty:
            label = "REFUTED"
        else:
            label = "OPEN"
        records[level] = {
            "label": label,
            "distance": mp_text(distance),
            "uncertainty": mp_text(uncertainty),
            "normalization": mp_text(normalization),
            "flint_conversion_bound": mp_text(flint_conversion / normalization),
            "binary_conversion_bound": mp_text(binary_conversion / normalization),
        }
    return records


def family_label(records, agree, refuted, open_label):
    labels = [record["label"] for record in records]
    if any(label == "REFUTED" for label in labels):
        return refuted
    if any(label == "OPEN" for label in labels):
        return open_label
    return agree


def scaled_analysis(analysis, factor):
    midpoints = {}
    radii = {}
    for level in LEVELS:
        midpoint, radius = scale_tangent_arrays(
            analysis["midpoints"][level], analysis["radii"][level], factor
        )
        midpoints[level] = midpoint
        radii[level] = radius
    return {"midpoints": midpoints, "radii": radii}


def product_analysis_record(left, right, core):
    midpoints = {}
    radii = {}
    defects = {}
    defect_radii = {}
    for level in LEVELS:
        midpoint, radius = product_with_radius(
            left["midpoints"][level],
            left["radii"][level],
            right["midpoints"][level],
            right["radii"][level],
        )
        defect, bound = defect_with_radius(midpoint, radius)
        midpoints[level] = midpoint
        radii[level] = radius
        defects[level] = defect
        defect_radii[level] = bound
    analysis = core["tangent_analysis"]({
        level: {
            "midpoint": midpoints[level],
            "radii": radii[level],
            "defect_midpoint": defects[level],
            "defect_radii": defect_radii[level],
        }
        for level in LEVELS
    })
    return analysis


def public_comparison(record):
    return {
        "label": record["label"],
        "distances": {
            key: float_text(value) for key, value in record["distances"].items()
        },
        "uncertainty": float_text(record["uncertainty"]),
    }


print("=" * 78)
print("SECOND FINITE-HEIGHT FULL TANGENT AND PHYSICAL TWO-STEP MAP")
print("=" * 78)

input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
payloads = {
    "history": json.loads(HISTORY_INPUT.read_text()),
    "first": json.loads(FIRST_INPUT.read_text()),
    "first_adversarial": json.loads(FIRST_ADVERSARIAL.read_text()),
    "composition": json.loads(COMPOSITION_INPUT.read_text()),
    "composition_adversarial": json.loads(COMPOSITION_ADVERSARIAL.read_text()),
    "old_two": json.loads(OLD_TWO_INPUT.read_text()),
}
provenance_ok = input_hashes == EXPECTED_HASHES
check(
    "all frozen inputs, corrections and implementation sources retain exact hashes",
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
    payloads["history"]["outcome"]
    == "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB"
    and payloads["first"]["outcome"]
    == "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY"
    and payloads["first"]["passed"] == payloads["first"]["tests"] == 21
    and payloads["first"]["numeric_archive"]["sha256"]
    == EXPECTED_HASHES["first_numeric"]
    and len(payloads["first"]["numeric_archive"]["arrays"]) == 168
    and payloads["first_adversarial"]["outcome"]
    == "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED"
    and payloads["first_adversarial"]["passed"]
    == payloads["first_adversarial"]["tests"] == 22
    and payloads["composition"]["outcome"]
    == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE"
    and payloads["composition"]["passed"]
    == payloads["composition"]["tests"] == 10
    and payloads["composition_adversarial"]["outcome"]
    == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED"
    and payloads["composition_adversarial"]["passed"]
    == payloads["composition_adversarial"]["tests"] == 9
    and payloads["old_two"]["outcome"]
    == "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED"
    and payloads["old_two"]["passed"] == payloads["old_two"]["tests"] == 16
    and payloads["old_two"]["target_comparisons_performed"] is False
    and payloads["old_two"]["numeric_archive_sha256"]
    == EXPECTED_HASHES["old_two_numeric"]
)
check(
    "the accepted first tangent, branch-B history and old method controls remain intact",
    accepted_inputs_ok,
)

history = reconstruct_history(
    payloads["first"]["background"],
    payloads["composition"],
    payloads["composition_adversarial"],
)
first_step = history["first"]
second_step = history["second"]
history_ok = bool(
    first_step["bracket_ok"]
    and second_step["bracket_ok"]
    and abs(first_step["residual"]) < mp.mpf("1e-140")
    and abs(second_step["residual"]) < mp.mpf("1e-140")
    and first_step["width"] < mp.mpf("1e-150")
    and second_step["width"] < mp.mpf("1e-150")
    and history["committed_error"] < mp.mpf("1e-65")
    and all(history["primary_serialized_matches"].values())
    and history["junction_error"] < mp.mpf("1e-140")
    and min(
        first_step["q"], first_step["h"], first_step["r"],
        first_step["m_out"], second_step["q"], second_step["h"],
        second_step["r"], second_step["m_out"],
        second_step["r"] - second_step["h"] ** 2,
    ) > 0
)
check(
    "deterministic bisection reconstructs the committed branch-B two-slab history",
    history_ok,
    (
        f"q1={mp_text(first_step['q'], 28)}, "
        f"q2={mp_text(second_step['q'], 28)}, "
        f"committed_error={mp_text(history['committed_error'], 6)}"
    ),
)

scale_control = exact_scale_control()
check(
    "the exact conformal-symplectic scale lift and non-unique canonical families pass",
    scale_control["passed"],
    str(scale_control),
)

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_second_finite_height_tangent", GEOMETRY_SOURCE
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

core = {
    "VERTICES": VERTICES,
    "acb": acb,
    "acb_mat": acb_mat,
    "ast": ast,
    "Counter": Counter,
    "io": io,
    "la": la,
    "math": math,
    "mp": mp,
    "np": np,
    "Path": Path,
    "zipfile": zipfile,
    "float_text": float_text,
    "mp_text": mp_text,
}
load_named_functions(
    FIRST_SOURCE,
    {
        "mp_frobenius",
        "mp_difference_frobenius",
        "mp_submatrix",
        "mp_to_numpy",
        "mp_to_acb",
        "mp_matrix_to_acb",
        "acb_midpoint_and_radii",
        "deterministic_npz",
        "expanded_types",
        "pre_legendre_matrix",
        "canonical_tangent_ball",
        "boundary_identification",
        "expanded_common_indices",
        "expanded_output_mapping",
        "scaled_frobenius",
        "scaled_difference_frobenius",
        "pre_legendre_analysis",
        "tangent_analysis",
    },
    core,
)

r1 = first_step["r"]
r2 = second_step["r"]
h1 = first_step["h"]
h2 = second_step["h"]
m0 = history["m0"]
m1 = first_step["m_out"]
scale_factor = r1**2

first_library = make_library(gro, m0, h1**2, mp.mpf(1))
normalized_library = make_library(gro, m1, h2**2, mp.mpf(1))
physical_library = make_library(
    gro, m0, scale_factor * h2**2, scale_factor
)
models = {
    parity: normalized_library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}

# Freeze the literal representation carrier at r1 and r2 before a Hessian.
prepared = {}
for parity in ("even", "odd"):
    model = models[parity]
    first_index = first_library["group_and_index_data"](
        model, (mp.log(r1), mp.mpf(0))
    )
    normalized_index = normalized_library["group_and_index_data"](
        model, (mp.log(r2), mp.mpf(0))
    )
    physical_index = physical_library["group_and_index_data"](
        model, (mp.log(r2), mp.mpf(0))
    )
    geometry_normalized = normalized_library["prepare_geometry"](
        model, normalized_index
    )
    geometry_physical = physical_library["prepare_geometry"](
        model, physical_index
    )
    first_boundary = core["boundary_identification"](first_index)
    boundary = core["boundary_identification"](normalized_index)
    physical_boundary = core["boundary_identification"](physical_index)
    first_sectors, first_sector_control = first_library[
        "high_precision_sector_bases"
    ](first_index)
    sectors, sector_control = normalized_library[
        "high_precision_sector_bases"
    ](normalized_index)

    carrier_ok = bool(
        len(model["old_edges"]) == OLD
        and len(model["internal_edges"]) == INTERNAL
        and len(model["new_edges"]) == NEW
        and len(model["slab"]) == 2400
        and len(geometry_normalized["triangle_records"]) == 6240
        and len(geometry_normalized["simplex_records"]) == 2400
        and len(geometry_normalized["patterns"]) == 20
        and geometry_normalized["patterns"] == geometry_physical["patterns"]
        and len(normalized_index["edge_to_index"]) == FULL
        and len(normalized_index["orbit_edges"]) == sum(ORBIT_ORDER)
        and first_boundary["ok"] and boundary["ok"] and physical_boundary["ok"]
    )
    basis_quality = bool(
        first_sector_control["irrep_dimensions"]
        == sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and first_sector_control["isotypic_dimensions"]
        == sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sum(sector["dimension"] ** 2 for sector in sectors) == 24
        and sum(60 * sector["dimension"] ** 2 for sector in sectors) == 1440
        and all(
            control[key] < mp.mpf("1e-140")
            for control in (first_sector_control, sector_control)
            for key in (
                "maximum_orthonormal",
                "maximum_center",
                "maximum_splitter",
                "maximum_right_leakage",
                "maximum_conjugate_pair",
            )
        )
    )
    basis_distance = max(
        core["mp_difference_frobenius"](left["basis"], right["basis"])
        for left, right in zip(first_sectors, sectors)
    )
    action_equal = bool(
        np.array_equal(first_index["table"], normalized_index["table"])
        and np.array_equal(normalized_index["table"], physical_index["table"])
        and all(
            np.array_equal(left, right)
            for left, right in zip(first_index["actions"], normalized_index["actions"])
        )
        and all(
            np.array_equal(left, right)
            for left, right in zip(normalized_index["actions"], physical_index["actions"])
        )
    )
    signature_equal = all(
        left["dimension"] == right["dimension"]
        and left["splitter"] == right["splitter"]
        and abs(
            left["old_central_eigenvalue"]
            - right["old_central_eigenvalue"]
        ) < mp.mpf("1e-140")
        for left, right in zip(first_sectors, sectors)
    )
    physical_edge_scale_error = max(
        abs(physical - scale_factor * normalized)
        for physical, normalized in zip(
            physical_index["signed_base"], normalized_index["signed_base"]
        )
    )
    common_ok = bool(
        action_equal
        and signature_equal
        and basis_distance < mp.mpf("1e-140")
        and first_boundary["mapping"] == boundary["mapping"]
        == physical_boundary["mapping"]
        and set(first_boundary["old_seeds"]) == set(boundary["old_seeds"])
        == set(physical_boundary["old_seeds"])
        and physical_edge_scale_error < mp.mpf("1e-140")
    )
    check(
        f"{parity}: r1/r2 carriers and the 720-edge boundary shift are literal",
        carrier_ok and common_ok,
        (
            f"basis_distance={mp_text(basis_distance, 6)}, "
            f"edge_scale_error={mp_text(physical_edge_scale_error, 6)}"
        ),
    )
    check(
        f"{parity}: seven target-free minimal sectors exhaust 1440 phase dimensions",
        basis_quality,
        f"dimensions={[sector['dimension'] for sector in sectors]}",
    )
    prepared[parity] = {
        "model": model,
        "first_index": first_index,
        "index": normalized_index,
        "physical_index": physical_index,
        "geometry": geometry_normalized,
        "physical_geometry": geometry_physical,
        "boundary": boundary,
        "sectors": sectors,
        "carrier_ok": carrier_ok,
        "basis_ok": basis_quality,
        "common_ok": common_ok,
        "basis_distance": basis_distance,
        "physical_edge_scale_error": physical_edge_scale_error,
    }

parity_basis_distance = max(
    core["mp_difference_frobenius"](left["basis"], right["basis"])
    for left, right in zip(
        prepared["even"]["sectors"], prepared["odd"]["sectors"]
    )
)
parity_carrier_ok = bool(
    np.array_equal(
        prepared["even"]["index"]["table"],
        prepared["odd"]["index"]["table"],
    )
    and all(
        np.array_equal(left, right)
        for left, right in zip(
            prepared["even"]["index"]["actions"],
            prepared["odd"]["index"]["actions"],
        )
    )
    and set(prepared["even"]["boundary"]["old_seeds"])
    == set(prepared["odd"]["boundary"]["old_seeds"])
    and parity_basis_distance < mp.mpf("1e-140")
)
check(
    "both staircase parities have literally common group and minimal-basis coordinates",
    parity_carrier_ok,
    f"basis_distance={mp_text(parity_basis_distance, 6)}",
)
common_seeds = tuple(sorted(prepared["even"]["boundary"]["old_seeds"]))

records = {}
runtime = {}
numeric_arrays = {}
for parity in ("even", "odd"):
    item = prepared[parity]
    kind_normalized = {
        "old": mp.mpf(1),
        "internal": r2 - h2**2,
        "pole": -h2**2,
        "new": r2**2,
    }
    kind_physical = {
        key: scale_factor * value for key, value in kind_normalized.items()
    }
    print(
        f"[{parity}] independently differentiating normalized and physical patterns",
        flush=True,
    )
    normalized_cache, normalized_branch = normalized_library[
        "high_precision_pattern_cache"
    ](item["geometry"]["patterns"], kind_normalized)
    physical_cache, physical_branch = physical_library[
        "high_precision_pattern_cache"
    ](item["physical_geometry"]["patterns"], kind_physical)
    branch_ok = bool(
        normalized_branch["entry_pass"] and physical_branch["entry_pass"]
        and normalized_branch["base_negative_counts"] == Counter({1: 2400})
        and physical_branch["base_negative_counts"] == Counter({1: 2400})
        and normalized_branch["displaced_negative_counts"] == Counter({1: 1600})
        and physical_branch["displaced_negative_counts"] == Counter({1: 1600})
        and min(
            normalized_branch["minimum_leading_minor"],
            physical_branch["minimum_leading_minor"],
        ) > 0
        and min(
            normalized_branch["minimum_argument"],
            physical_branch["minimum_argument"],
        ) > mp.mpf("1e-6")
    )
    check(
        f"{parity}: both direct derivative assemblies retain the Lorentzian branch",
        branch_ok,
        (
            f"normalized_cross={mp_text(normalized_branch['maximum_cross'], 6)}, "
            f"physical_cross={mp_text(physical_branch['maximum_cross'], 6)}"
        ),
    )

    raw_normalized, normalized_kernel_control = normalized_library[
        "assemble_full_representative_kernels"
    ](item["index"], item["geometry"], normalized_cache)
    raw_physical, physical_kernel_control = physical_library[
        "assemble_full_representative_kernels"
    ](item["physical_index"], item["physical_geometry"], physical_cache)
    richardson_normalized = {
        "K01": richardson(
            raw_normalized["operational_primary"],
            raw_normalized["operational_shadow"],
        ),
        "K12": richardson(
            raw_normalized["operational_shadow"],
            raw_normalized["validation_primary"],
        ),
        "K23": richardson(
            raw_normalized["validation_primary"],
            raw_normalized["validation_shadow"],
        ),
    }
    richardson_physical = {
        "K01": richardson(
            raw_physical["operational_primary"],
            raw_physical["operational_shadow"],
        ),
        "K12": richardson(
            raw_physical["operational_shadow"],
            raw_physical["validation_primary"],
        ),
        "K23": richardson(
            raw_physical["validation_primary"],
            raw_physical["validation_shadow"],
        ),
    }
    maximum_imaginary = max(
        normalized_kernel_control["maximum_imaginary"],
        physical_kernel_control["maximum_imaginary"],
        max(
            abs(mp.im(value))
            for family in (richardson_normalized, richardson_physical)
            for kernel in family.values()
            for value in kernel.values()
        ),
    )
    kernel_ok = bool(
        set(raw_normalized) == set(raw_physical) == set(VARIANTS)
        and all(len(kernel) > 0 for kernel in raw_normalized.values())
        and all(len(kernel) > 0 for kernel in raw_physical.values())
        and maximum_imaginary < mp.mpf("1e-140")
    )
    raw_scale = kernel_scale_records(
        raw_physical, raw_normalized, scale_factor, VARIANTS
    )
    richardson_scale = kernel_scale_records(
        richardson_physical, richardson_normalized, scale_factor, LEVELS
    )
    parity_scale_label = family_label(
        list(raw_scale.values()) + list(richardson_scale.values()),
        "SCALE_LIFT_CONFIRMED",
        "SCALE_LIFT_REFUTED",
        "SCALE_LIFT_OPEN",
    )
    check(
        f"{parity}: independently assembled physical kernels satisfy the frozen scale classifier",
        kernel_ok and parity_scale_label in {
            "SCALE_LIFT_CONFIRMED", "SCALE_LIFT_REFUTED", "SCALE_LIFT_OPEN"
        },
        f"label={parity_scale_label}, imaginary={mp_text(maximum_imaginary, 6)}",
    )

    sector_records = []
    sector_runtime = []
    all_reciprocity = True
    all_regular = True
    all_canonical = True
    for sector_index, sector in enumerate(item["sectors"]):
        dimension = int(sector["dimension"])
        print(
            f"[{parity}] sector {sector_index + 1}/7 d={dimension}: J and T",
            flush=True,
        )
        normalized_blocks = {
            level: normalized_library["project_full_kernel"](kernel, sector)
            for level, kernel in richardson_normalized.items()
        }
        physical_blocks = {
            level: physical_library["project_full_kernel"](kernel, sector)
            for level, kernel in richardson_physical.items()
        }
        reciprocity = {}
        sector_reciprocity = True
        for label, blocks in (
            ("normalized", normalized_blocks), ("physical", physical_blocks)
        ):
            normalization = max(
                mp.mpf(1), core["mp_frobenius"](blocks["K12"])
            )
            variation = max(
                core["mp_difference_frobenius"](blocks["K01"], blocks["K12"]),
                core["mp_difference_frobenius"](blocks["K12"], blocks["K23"]),
            ) / normalization
            residual = max(
                core["mp_difference_frobenius"](block, block.H) / normalization
                for block in blocks.values()
            )
            uncertainty = variation + CLASSIFIER_FLOOR
            ok = bool(residual <= 10 * uncertainty)
            sector_reciprocity &= ok
            reciprocity[label] = {
                "residual": mp_text(residual),
                "variation": mp_text(variation),
                "uncertainty": mp_text(uncertainty),
                "passed": ok,
            }
        all_reciprocity &= sector_reciprocity

        normalized_rank = core["pre_legendre_analysis"](
            normalized_blocks, dimension
        )
        physical_rank = core["pre_legendre_analysis"](
            physical_blocks, dimension
        )
        sector_regular = bool(
            normalized_rank["regular"] and physical_rank["regular"]
        )
        all_regular &= sector_regular
        tangent_records = {}
        tangent_runtime = {}
        scale_comparison = {
            "label": "OPEN", "distances": {level: math.inf for level in LEVELS},
            "uncertainty": math.inf,
        }
        identity_comparison = scale_comparison
        linear_comparison = scale_comparison
        if sector_regular:
            old = core["expanded_types"](0, 30, dimension)
            internal = core["expanded_types"](30, 65, dimension)
            new = core["expanded_types"](65, 95, dimension)
            output_mapping = core["expanded_output_mapping"](
                item["boundary"]["mapping"], dimension
            )
            common_indices = core["expanded_common_indices"](
                item["boundary"]["old_seeds"], common_seeds, dimension
            )
            for label, blocks in (
                ("normalized", normalized_blocks), ("physical", physical_blocks)
            ):
                built_by_level = {
                    level: core["canonical_tangent_ball"](
                        block, old, internal, new, output_mapping, common_indices
                    )
                    for level, block in blocks.items()
                }
                build_ok = all(
                    built["tangent"] is not None
                    and not built["det_j"].contains(0)
                    for built in built_by_level.values()
                )
                if build_ok:
                    ball_records = ball_records_from_built(
                        built_by_level, core["acb_midpoint_and_radii"]
                    )
                    analysis = core["tangent_analysis"](ball_records)
                    tangent_runtime[label] = analysis
                    tangent_records[label] = analysis["record"]
                    for level in LEVELS:
                        prefix = (
                            f"{parity}_sector{sector_index}_{label}_{level}"
                        )
                        numeric_arrays[f"{prefix}_tangent_midpoint"] = (
                            ball_records[level]["midpoint"]
                        )
                        numeric_arrays[f"{prefix}_tangent_radii"] = (
                            ball_records[level]["radii"]
                        )
                        numeric_arrays[f"{prefix}_defect_midpoint"] = (
                            ball_records[level]["defect_midpoint"]
                        )
                        numeric_arrays[f"{prefix}_defect_radii"] = (
                            ball_records[level]["defect_radii"]
                        )
                else:
                    tangent_records[label] = {
                        "canonicality_ok": False,
                        "reason": "Flint solve failed after regular classification",
                    }
            if set(tangent_runtime) == {"normalized", "physical"}:
                scale_comparison = comparison_record(
                    tangent_runtime["physical"],
                    scaled_analysis(tangent_runtime["normalized"], scale_factor),
                )
                identity_comparison = comparison_record(
                    tangent_runtime["physical"],
                    scaled_analysis(tangent_runtime["normalized"], mp.mpf(1)),
                )
                linear_comparison = comparison_record(
                    tangent_runtime["physical"],
                    scaled_analysis(tangent_runtime["normalized"], r1),
                )
                sector_canonical = bool(
                    tangent_runtime["normalized"]["canonicality_ok"]
                    and tangent_runtime["physical"]["canonicality_ok"]
                    and scale_comparison["label"] == "AGREES"
                )
            else:
                sector_canonical = False
        else:
            sector_canonical = False
        all_canonical &= sector_canonical

        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "constant_overlap": mp_text(sector["constant_overlap"]),
            "center_value": mp_text(sector["old_central_eigenvalue"]),
            "splitter": sector["splitter"],
            "reciprocity": reciprocity,
            "normalized_pre_legendre": normalized_rank["record"],
            "physical_pre_legendre": physical_rank["record"],
            "tangents": tangent_records,
            "correct_scale_lift": public_comparison(scale_comparison),
            "identity_lift_hostile": public_comparison(identity_comparison),
            "linear_r1_lift_hostile": public_comparison(linear_comparison),
        })
        sector_runtime.append({
            "tangents": tangent_runtime,
            "scale": scale_comparison,
            "identity": identity_comparison,
            "linear": linear_comparison,
            "regular": sector_regular,
            "canonical": sector_canonical,
        })

    check(
        f"{parity}: all projected normalized and physical Hessians retain reciprocity",
        all_reciprocity,
    )
    check(
        f"{parity}: all fourteen pre-Legendre blocks receive the frozen classifier",
        len(sector_records) == 7
        and all(
            record[key]["classification"] in {"REGULAR", "NUMERICALLY_OPEN"}
            for record in sector_records
            for key in ("normalized_pre_legendre", "physical_pre_legendre")
        ),
        (
            f"regular={sum(r['regular'] for r in sector_runtime) * 2}/14"
        ),
    )
    check(
        f"{parity}: every regular direct tangent receives canonical and scale-lift labels",
        all(
            record["scale"]["label"] in {"AGREES", "REFUTED", "OPEN"}
            for record in sector_runtime
        ),
        f"canonical={sum(r['canonical'] for r in sector_runtime)}/7",
    )
    records[parity] = {
        "controls_ok": bool(
            item["carrier_ok"] and item["basis_ok"] and item["common_ok"]
            and branch_ok and kernel_ok and all_reciprocity
        ),
        "branch_ok": branch_ok,
        "kernel_ok": kernel_ok,
        "scale_lift_label": parity_scale_label,
        "raw_kernel_scale": raw_scale,
        "richardson_kernel_scale": richardson_scale,
        "all_regular": all_regular,
        "all_canonical": all_canonical,
        "sectors": sector_records,
    }
    runtime[parity] = sector_runtime

kernel_scale_labels = [records[parity]["scale_lift_label"] for parity in ("even", "odd")]
if "SCALE_LIFT_REFUTED" in kernel_scale_labels:
    scale_lift_outcome = "SCALE_LIFT_REFUTED"
elif "SCALE_LIFT_OPEN" in kernel_scale_labels:
    scale_lift_outcome = "SCALE_LIFT_OPEN"
else:
    scale_lift_outcome = "SCALE_LIFT_CONFIRMED"
check(
    "the seven raw/Richardson scale comparisons per parity freeze one global label",
    scale_lift_outcome in {
        "SCALE_LIFT_CONFIRMED", "SCALE_LIFT_REFUTED", "SCALE_LIFT_OPEN"
    },
    f"parities={kernel_scale_labels}, global={scale_lift_outcome}",
)

all_regular = bool(all(records[parity]["all_regular"] for parity in ("even", "odd")))
all_canonical = bool(
    all(records[parity]["all_canonical"] for parity in ("even", "odd"))
)
maps_available = bool(
    all(
        set(record["tangents"]) == {"normalized", "physical"}
        for parity in ("even", "odd")
        for record in runtime[parity]
    )
)
correct_tangent_scale_ok = bool(
    maps_available
    and all(
        record["scale"]["label"] == "AGREES"
        for parity in ("even", "odd")
        for record in runtime[parity]
    )
)
identity_tangent_refuted = bool(
    maps_available
    and any(
        record["identity"]["label"] == "REFUTED"
        for parity in ("even", "odd")
        for record in runtime[parity]
    )
)
linear_tangent_refuted = bool(
    maps_available
    and any(
        record["linear"]["label"] == "REFUTED"
        for parity in ("even", "odd")
        for record in runtime[parity]
    )
)
check(
    "the direct physical tangent classifier can reject both wrong unit lifts when applicable",
    (not maps_available) or (identity_tangent_refuted and linear_tangent_refuted),
    (
        f"available={maps_available}, identity={identity_tangent_refuted}, "
        f"linear={linear_tangent_refuted}"
    ),
)

second_schedule_records = []
if maps_available:
    for sector_index in range(7):
        comparison = comparison_record(
            runtime["even"][sector_index]["tangents"]["physical"],
            runtime["odd"][sector_index]["tangents"]["physical"],
        )
        second_schedule_records.append({
            "sector_index": sector_index,
            "dimension": records["even"]["sectors"][sector_index]["dimension"],
            **public_comparison(comparison),
        })
second_labels = [record["label"] for record in second_schedule_records]
if not maps_available:
    second_schedule_outcome = "NOT_EVALUATED"
elif "REFUTED" in second_labels:
    second_schedule_outcome = "SCHEDULE_DEPENDENT"
elif "OPEN" in second_labels:
    second_schedule_outcome = "SCHEDULE_OPEN"
else:
    second_schedule_outcome = "SCHEDULE_ROBUST"
check(
    "all seven second-slab schedule comparisons receive the frozen 10/100 label",
    (not maps_available) or len(second_schedule_records) == 7,
    f"outcome={second_schedule_outcome}, labels={dict(Counter(second_labels))}",
)

# This assignment is the protocol firewall: no first tangent entry is opened
# above this line.
frozen_second_labels = {
    "kernel_scale": scale_lift_outcome,
    "all_regular": all_regular,
    "all_canonical": all_canonical,
    "correct_tangent_scale": correct_tangent_scale_ok,
    "second_schedule": second_schedule_outcome,
}
first_numeric = np.load(FIRST_NUMERIC)
first_expected_keys = {
    f"{parity}_sector{sector}_{level}_{suffix}"
    for parity in ("even", "odd")
    for sector in range(7)
    for level in LEVELS
    for suffix in (
        "tangent_midpoint", "tangent_radii",
        "defect_midpoint", "defect_radii",
    )
}
first_archive_ok = bool(
    len(first_numeric.files) == 168
    and set(first_numeric.files) == first_expected_keys
)
check(
    "the first tangent archive is opened only after second labels and has its exact schema",
    bool(frozen_second_labels) and first_archive_ok,
    f"arrays={len(first_numeric.files)}",
)

first_runtime = {parity: [] for parity in ("even", "odd")}
if first_archive_ok:
    for parity in ("even", "odd"):
        for sector_index in range(7):
            first_runtime[parity].append({
                "midpoints": {
                    level: first_numeric[
                        f"{parity}_sector{sector_index}_{level}_tangent_midpoint"
                    ]
                    for level in LEVELS
                },
                "radii": {
                    level: first_numeric[
                        f"{parity}_sector{sector_index}_{level}_tangent_radii"
                    ]
                    for level in LEVELS
                },
            })

schedule_pairs = tuple(
    (first_parity, second_parity)
    for first_parity in ("even", "odd")
    for second_parity in ("even", "odd")
)
product_runtime = {
    f"{first_parity}_{second_parity}": []
    for first_parity, second_parity in schedule_pairs
}
products_available = bool(maps_available and first_archive_ok)
if products_available:
    for first_parity, second_parity in schedule_pairs:
        pair_key = f"{first_parity}_{second_parity}"
        for sector_index in range(7):
            analysis = product_analysis_record(
                runtime[second_parity][sector_index]["tangents"]["physical"],
                first_runtime[first_parity][sector_index],
                core,
            )
            product_runtime[pair_key].append(analysis)
            for level in LEVELS:
                prefix = f"product_{pair_key}_sector{sector_index}_{level}"
                numeric_arrays[f"{prefix}_tangent_midpoint"] = analysis[
                    "midpoints"
                ][level]
                numeric_arrays[f"{prefix}_tangent_radii"] = analysis[
                    "radii"
                ][level]
                numeric_arrays[f"{prefix}_defect_midpoint"] = analysis[
                    "defects"
                ][level]
                numeric_arrays[f"{prefix}_defect_radii"] = analysis[
                    "defect_radii"
                ][level]

all_products_canonical = bool(
    products_available
    and all(
        analysis["canonicality_ok"]
        for analyses in product_runtime.values()
        for analysis in analyses
    )
)
check(
    "all four physical two-step products satisfy propagated canonical bounds when applicable",
    (not products_available) or all_products_canonical,
    f"available={products_available}",
)

product_schedule_records = []
if products_available:
    pair_keys = tuple(product_runtime)
    for sector_index in range(7):
        for left_key, right_key in combinations(pair_keys, 2):
            comparison = comparison_record(
                product_runtime[left_key][sector_index],
                product_runtime[right_key][sector_index],
            )
            product_schedule_records.append({
                "sector_index": sector_index,
                "dimension": records["even"]["sectors"][sector_index]["dimension"],
                "left": left_key,
                "right": right_key,
                **public_comparison(comparison),
            })
product_labels = [record["label"] for record in product_schedule_records]
if not products_available:
    product_schedule_outcome = "NOT_EVALUATED"
elif "REFUTED" in product_labels:
    product_schedule_outcome = "TWO_STEP_SCHEDULE_DEPENDENT"
elif "OPEN" in product_labels:
    product_schedule_outcome = "TWO_STEP_SCHEDULE_OPEN"
else:
    product_schedule_outcome = "TWO_STEP_SCHEDULE_ROBUST"
check(
    "all 42 pairwise comparisons among four schedule products receive labels",
    (not products_available) or len(product_schedule_records) == 42,
    f"outcome={product_schedule_outcome}, labels={dict(Counter(product_labels))}",
)

synthetic_detections = []
identity_product_records = []
linear_product_records = []
if products_available:
    for sector_index in range(7):
        correct = product_runtime["even_even"][sector_index]
        synthetic = {
            "midpoints": {
                level: correct["midpoints"][level].copy() for level in LEVELS
            },
            "radii": {
                level: correct["radii"][level].copy() for level in LEVELS
            },
        }
        for level in LEVELS:
            synthetic["midpoints"][level][0, 0] += 1e-3
        synthetic_record = comparison_record(correct, synthetic)
        synthetic_detections.append(synthetic_record["label"] == "REFUTED")

        normalized_second = runtime["even"][sector_index]["tangents"][
            "normalized"
        ]
        wrong_identity = product_analysis_record(
            scaled_analysis(normalized_second, mp.mpf(1)),
            first_runtime["even"][sector_index],
            core,
        )
        wrong_linear = product_analysis_record(
            scaled_analysis(normalized_second, r1),
            first_runtime["even"][sector_index],
            core,
        )
        identity_record = comparison_record(correct, wrong_identity)
        linear_record = comparison_record(correct, wrong_linear)
        identity_product_records.append({
            "sector_index": sector_index,
            **public_comparison(identity_record),
        })
        linear_product_records.append({
            "sector_index": sector_index,
            **public_comparison(linear_record),
        })

synthetic_detected = bool(
    products_available and synthetic_detections and all(synthetic_detections)
)
identity_product_refuted = bool(
    products_available
    and any(record["label"] == "REFUTED" for record in identity_product_records)
)
linear_product_refuted = bool(
    products_available
    and any(record["label"] == "REFUTED" for record in linear_product_records)
)
check(
    "the product classifier detects the synthetic schedule corruption when applicable",
    (not products_available) or synthetic_detected,
    f"detected={sum(synthetic_detections)}/{len(synthetic_detections)}",
)
check(
    "the direct physical product rejects identity and r1 hostile lifts in an actual sector",
    (not products_available)
    or (identity_product_refuted and linear_product_refuted),
    f"identity={identity_product_refuted}, linear={linear_product_refuted}",
)

post_map_controls_ok = bool(
    (not maps_available)
    or (identity_tangent_refuted and linear_tangent_refuted)
)
post_product_controls_ok = bool(
    (not products_available)
    or (synthetic_detected and identity_product_refuted and linear_product_refuted)
)
base_controls_ok = bool(
    provenance_ok and registry_ok and accepted_inputs_ok and history_ok
    and scale_control["passed"] and geometry_import_ok and parity_carrier_ok
    and first_archive_ok
    and all(records[parity]["controls_ok"] for parity in ("even", "odd"))
    and post_map_controls_ok and post_product_controls_ok
)

if not base_controls_ok:
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_CONTROL_FAILED"
elif scale_lift_outcome == "SCALE_LIFT_REFUTED":
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_SCALE_LIFT_REFUTED"
elif scale_lift_outcome == "SCALE_LIFT_OPEN":
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_SCALE_LIFT_OPEN"
elif not all_regular:
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_RANK_OPEN"
elif not (
    all_canonical and correct_tangent_scale_ok and all_products_canonical
):
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED"
elif (
    second_schedule_outcome == "SCHEDULE_DEPENDENT"
    or product_schedule_outcome == "TWO_STEP_SCHEDULE_DEPENDENT"
):
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_SCHEDULE_DEPENDENT"
elif (
    second_schedule_outcome == "SCHEDULE_OPEN"
    or product_schedule_outcome == "TWO_STEP_SCHEDULE_OPEN"
):
    outcome = "SECOND_FULL_BOUNDARY_TANGENT_SCHEDULE_OPEN"
else:
    outcome = "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY"

allowed_outcomes = {
    "SECOND_FULL_BOUNDARY_TANGENT_CONTROL_FAILED",
    "SECOND_FULL_BOUNDARY_TANGENT_SCALE_LIFT_REFUTED",
    "SECOND_FULL_BOUNDARY_TANGENT_SCALE_LIFT_OPEN",
    "SECOND_FULL_BOUNDARY_TANGENT_RANK_OPEN",
    "SECOND_FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED",
    "SECOND_FULL_BOUNDARY_TANGENT_SCHEDULE_DEPENDENT",
    "SECOND_FULL_BOUNDARY_TANGENT_SCHEDULE_OPEN",
    "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY",
}
check(
    "the result follows the frozen eight-level outcome hierarchy",
    outcome in allowed_outcomes,
    outcome,
)

expected_numeric_arrays = 672
archive_schema_ok = bool(
    (not products_available) or len(numeric_arrays) == expected_numeric_arrays
)
core["deterministic_npz"](NUMERIC_OUTPUT, numeric_arrays)
numeric_hash = sha256(NUMERIC_OUTPUT)
check(
    "the deterministic archive has the complete schema when all maps exist",
    archive_schema_ok,
    f"arrays={len(numeric_arrays)}, expected_if_complete={expected_numeric_arrays}",
)

artifact = {
    "outcome": outcome,
    "status": "PRIMARY_ONLY_ADVERSARIAL_REPLICATION_REQUIRED_FOR_MATERIAL_RESULT",
    "tests": tests,
    "passed": passed,
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "repository_prior_art_correction_commit": CORRECTION_COMMIT,
        "history_provenance_correction_commit": HISTORY_CORRECTION_COMMIT,
        "first_control_failure_commit": FIRST_FAILURE_COMMIT,
        "converter_failure_commit": SECOND_FAILURE_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": input_hashes,
    },
    "history": {
        "v": mp_text(history["v"]),
        "m0": mp_text(history["m0"]),
        "pi0": mp_text(history["pi0"]),
        "q1": mp_text(first_step["q"]),
        "h1": mp_text(first_step["h"]),
        "r1": mp_text(first_step["r"]),
        "m1": mp_text(first_step["m_out"]),
        "pi1": mp_text(first_step["pi_out"]),
        "q2": mp_text(second_step["q"]),
        "h2": mp_text(second_step["h"]),
        "r2": mp_text(second_step["r"]),
        "m2": mp_text(second_step["m_out"]),
        "pi2": mp_text(second_step["pi_out"]),
        "committed_error": mp_text(history["committed_error"]),
        "junction_error": mp_text(history["junction_error"]),
        "primary_60_digit_serialized_matches": history[
            "primary_serialized_matches"
        ],
    },
    "scale_factor_r1_squared": mp_text(scale_factor),
    "exact_scale_control": scale_control,
    "frozen_second_labels_before_first_tangent_open": frozen_second_labels,
    "kernel_scale_lift_outcome": scale_lift_outcome,
    "all_pre_legendre_sectors_regular": all_regular,
    "all_direct_tangents_canonical_and_scale_consistent": all_canonical,
    "second_schedule_outcome": second_schedule_outcome,
    "two_step_schedule_outcome": product_schedule_outcome,
    "second_schedule_comparisons": second_schedule_records,
    "two_step_schedule_comparisons": product_schedule_records,
    "hostile_controls": {
        "identity_tangent_lift_refuted": identity_tangent_refuted,
        "linear_r1_tangent_lift_refuted": linear_tangent_refuted,
        "synthetic_schedule_corruption_detected_all_sectors": synthetic_detected,
        "identity_product_lift_refuted": identity_product_refuted,
        "linear_r1_product_lift_refuted": linear_product_refuted,
        "identity_product_records": identity_product_records,
        "linear_r1_product_records": linear_product_records,
    },
    "parities": records,
    "product_canonicality": {
        pair: [analysis["record"] for analysis in analyses]
        for pair, analyses in product_runtime.items()
    },
    "numeric_archive": {
        "path": NUMERIC_OUTPUT.name,
        "sha256": numeric_hash,
        "array_count": len(numeric_arrays),
        "arrays": {key: list(value.shape) for key, value in numeric_arrays.items()},
    },
    "classification": {
        "two_step_finite_height_response": (
            "DERIVED_COMPUTATIONAL_PRIMARY"
            if outcome == "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY"
            else "OPEN_OR_REFUTED_BY_OUTCOME"
        ),
        "method_novelty": "REFUTED_REPOSITORY_KNOWN",
        "new_background_coefficients_external_novelty": "OPEN",
        "physical_constraint_quotient": "OPEN",
        "physical_mode_spectrum": "NOT_COMPUTED",
        "wave_equation_limiting_speed_G_planck_particles": "NOT_DERIVED",
        "homogeneous_convergence_duration_completeness": "OPEN",
    },
    "firewall": {
        "first_tangent_entries_opened_before_second_labels": False,
        "tangent_eigenvalues_computed": False,
        "tangent_singular_values_computed": False,
        "continuum_or_particle_target_parsed": False,
        "old_two_step_spectrum_parsed": False,
        "full_suite_run": False,
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"OUTCOME: {outcome}")
print(f"KERNEL SCALE: {scale_lift_outcome}")
print(f"SECOND SCHEDULE: {second_schedule_outcome}")
print(f"TWO-STEP SCHEDULE: {product_schedule_outcome}")
for parity in ("even", "odd"):
    print(
        f"{parity}: regular={sum(r['regular'] for r in runtime[parity])}/7, "
        f"canonical={sum(r['canonical'] for r in runtime[parity])}/7",
        flush=True,
    )
print(f"RESULT: {passed}/{tests} PASS")
print(f"NUMERIC SHA: {numeric_hash}")
if passed != tests:
    raise SystemExit(1)
