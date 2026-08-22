#!/usr/bin/env python3
"""Target-free quadratic schedule test on the finite-height 240 carrier.

Prior-art commit: 5fab3f5.
Protocol commit: 01382fb.
Provenance clarifications: a020978, 2f9361b.
Registry commit: a6c93d4.

No quadratic spectrum, continuum mode, speed, or desired eigenvalue is read.
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
OUTPUT = HERE / "gravity_600cell_finite_height_carrier_quadratic.json"
MATRIX_OUTPUT = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_matrices.npy"
)
RUN_ALL = HERE / "run_all.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
HESSIAN_SOURCE = (
    HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
)
CARRIER_SOURCE = (
    HERE / "verify_gravity_600cell_full_scale_strut_canonical_intersection.py"
)

INPUTS = {
    "finite_height": HERE / "gravity_600cell_finite_height_fourth_slab.json",
    "homothetic_action": HERE / "gravity_600cell_homothetic_frustum_action.json",
    "carrier": HERE / "gravity_600cell_full_scale_strut_carrier.json",
    "carrier_precision": HERE / "gravity_600cell_full_scale_strut_precision.json",
    "canonical_precision": (
        HERE / "gravity_600cell_full_scale_strut_canonical_precision.json"
    ),
    "symbolic_gap": (
        HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
    ),
    "schedule_census": (
        HERE / "gravity_600cell_staircase_orientation_selector.json"
    ),
    "geometry_source": GEOMETRY_SOURCE,
    "hessian_source": HESSIAN_SOURCE,
    "carrier_source": CARRIER_SOURCE,
}

EXPECTED_HASHES = {
    "finite_height": (
        "cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468"
    ),
    "homothetic_action": (
        "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
    ),
    "carrier": (
        "6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d"
    ),
    "carrier_precision": (
        "2a2a79271a92fc2ddde343a9d0651402df6eeb4a90efa2697e26f54cafcdf60f"
    ),
    "canonical_precision": (
        "75351ae4dfde26dd75ed8faa927b0a49cd725d83c7629d4545268030b54e2706"
    ),
    "symbolic_gap": (
        "ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179"
    ),
    "schedule_census": (
        "6a0bf4112baf3868a728beb45f04e9bdf1420cfbb93b26e6b8680041cb5d37f2"
    ),
    "geometry_source": (
        "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"
    ),
    "hessian_source": (
        "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5"
    ),
    "carrier_source": (
        "a2d5390d39c725a5fb586fefce9da34cede3a1fb84bbe36791f8b0599b3eae42"
    ),
}

PRIOR_ART_COMMIT = "5fab3f5"
PROTOCOL_COMMIT = "01382fb"
SCHEDULE_PROVENANCE_COMMIT = "a020978"
COMPLETE_PROVENANCE_COMMIT = "2f9361b"
REGISTRY_COMMIT = "a6c93d4"
VERIFIER_NAME = Path(__file__).name

DPS = 120
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
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


def norm2(matrix):
    values = la.svdvals(np.asarray(matrix), check_finite=True)
    return float(values[0]) if len(values) else 0.0


def frobenius_difference(left, right):
    difference = np.asarray(left) - np.asarray(right)
    result = float(np.linalg.norm(difference))
    del difference
    return result


def load_named_functions(path, names, namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name in names
    ]
    found = {node.name for node in selected}
    missing = set(names) - found
    if missing:
        raise RuntimeError(f"missing AST-loaded functions: {sorted(missing)}")
    module = ast.Module(body=selected, type_ignores=[])
    ast.fix_missing_locations(module)
    exec(compile(module, str(path), "exec"), namespace)


def registry_inventory(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    scripts = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "scripts"
               for target in node.targets):
            scripts = ast.literal_eval(node.value)
            break
    if scripts is None:
        raise RuntimeError("run_all.py has no literal scripts registry")
    counts = Counter(scripts)
    duplicates = sorted(name for name, count in counts.items() if count != 1)
    return scripts, duplicates


def finite_height_state():
    def epsilon(value):
        square = value**2
        argument = (square + 2) / (2 * (square + 3))
        return 2 * mp.pi - 5 * mp.acos(argument)

    def mu(value):
        return 180 * epsilon(value) / (mp.pi * mp.sqrt(value**2 + 4))

    def momentum(value):
        square = value**2
        return (
            180 * value * epsilon(value) / mp.sqrt(square + 4)
            - 600 * mp.sqrt(3)
            * mp.asinh(value / mp.sqrt(8 * (square + 3)))
        )

    v = mp.mpf(3) / 2
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
    if not bracket_ok:
        return {
            "bracket_ok": False,
            "left_value": left_value,
            "right_value": right_value,
        }
    for _ in range(500):
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
        if right - left < mp.mpf("1e-110"):
            break
    q = (left + right) / 2
    residual = elimination(q)
    h = (momentum(q) - incoming) / (2 * mp.pi * mu(q))
    lam = 1 + h * q
    rho = h**2
    q_diag = lam - rho
    tau_square = 3 * (lam - 1) ** 2 + 8 * rho
    exceptional = (lam - 1) ** 2 - 3 * tau_square
    return {
        "bracket_ok": True,
        "v": v,
        "q": q,
        "h": h,
        "lambda": lam,
        "rho": rho,
        "mass": mass,
        "incoming_momentum": incoming,
        "residual": residual,
        "width": right - left,
        "q_diag": q_diag,
        "tau_square": tau_square,
        "exceptional_factor": exceptional,
    }


def build_carrier(model, index_data, lam, rho):
    q_diag = lam - rho
    coefficient_a = -16 * rho / (lam - 1) ** 2
    coefficient_b = 8 + 16 * rho / (lam - 1) ** 2
    a_value = coefficient_a / (8 * q_diag)
    b_value = coefficient_b / (8 * q_diag)
    kappa = rho / ((lam - 1) * q_diag)

    matrix = np.zeros((ACTIVE, DATA), dtype=float)
    edge_by_row = [None] * ACTIVE
    diagonal_rows = []
    pole_rows = []
    for edge, global_index in index_data["edge_to_index"].items():
        if global_index < OLD:
            continue
        row = global_index - OLD
        edge = tuple(int(value) for value in edge)
        edge_by_row[row] = edge
        kind = index_data["edge_kind"][global_index]
        if kind == "pole":
            lower, upper = edge
            if upper != lower + VERTICES:
                raise RuntimeError("pole endpoint convention failed")
            matrix[row, VERTICES + lower] = 1.0
            pole_rows.append((lower, row))
        elif kind == "internal":
            lower, upper = edge
            target = upper - VERTICES
            if not (0 <= lower < VERTICES and 0 <= target < VERTICES):
                raise RuntimeError("oriented diagonal endpoint convention failed")
            matrix[row, lower] = float(a_value)
            matrix[row, target] = float(b_value)
            matrix[row, VERTICES + lower] = float(kappa)
            matrix[row, VERTICES + target] = float(-lam * kappa)
            diagonal_rows.append((edge, row, lower))
        elif kind == "new":
            left, right = (edge[0] - VERTICES, edge[1] - VERTICES)
            matrix[row, left] = float(1 / lam)
            matrix[row, right] = float(1 / lam)
        else:
            raise RuntimeError(f"unexpected active edge kind: {kind}")

    if any(edge is None for edge in edge_by_row):
        raise RuntimeError("carrier row coverage is incomplete")

    row_support = Counter(np.count_nonzero(matrix, axis=1).tolist())
    scale_support = np.count_nonzero(matrix[:, :VERTICES], axis=0)
    strut_support = np.count_nonzero(matrix[:, VERTICES:], axis=0)

    pole_rows.sort()
    pole_block = matrix[
        np.asarray([row for _, row in pole_rows]), VERTICES:
    ]
    pole_identity = bool(np.array_equal(pole_block, np.eye(VERTICES)))

    graph = nx.Graph()
    graph.add_nodes_from(range(VERTICES))
    for edge in model["new_edges"]:
        graph.add_edge(edge[0] - VERTICES, edge[1] - VERTICES)
    graph_control = {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "connected": nx.is_connected(graph),
        "bipartite": nx.is_bipartite(graph),
    }
    exact_rank_proof = bool(
        pole_identity
        and graph_control["vertices"] == VERTICES
        and graph_control["edges"] == NEW
        and graph_control["connected"]
        and not graph_control["bipartite"]
    )

    collective_exact_error = max(
        abs(a_value + b_value - 1 / q_diag),
        abs(kappa - lam * kappa + rho / q_diag),
    )
    collective_numeric_error = mp.mpf(0)
    scale_sum = matrix[:, :VERTICES].sum(axis=1)
    strut_sum = matrix[:, VERTICES:].sum(axis=1)
    for row, edge in enumerate(edge_by_row):
        global_index = OLD + row
        kind = index_data["edge_kind"][global_index]
        if kind == "pole":
            expected_scale, expected_strut = mp.mpf(0), mp.mpf(1)
        elif kind == "internal":
            expected_scale = 1 / q_diag
            expected_strut = -rho / q_diag
        else:
            expected_scale, expected_strut = 2 / lam, mp.mpf(0)
        collective_numeric_error = max(
            collective_numeric_error,
            abs(mp.mpf(scale_sum[row]) - expected_scale),
            abs(mp.mpf(strut_sum[row]) - expected_strut),
        )

    controls = {
        "row_support": dict(sorted(row_support.items())),
        "scale_support_minimum": int(scale_support.min()),
        "scale_support_maximum": int(scale_support.max()),
        "strut_support_minimum": int(strut_support.min()),
        "strut_support_maximum": int(strut_support.max()),
        "pole_identity": pole_identity,
        "graph": graph_control,
        "exact_rank_240": exact_rank_proof,
        "collective_exact_error": mp_text(collective_exact_error),
        "collective_numeric_error": mp_text(collective_numeric_error),
        "coefficients": {
            "a": mp_text(a_value),
            "b": mp_text(b_value),
            "kappa": mp_text(kappa),
        },
    }
    controls_ok = bool(
        row_support == Counter({1: 120, 2: 720, 4: 720})
        and np.all(scale_support == 24)
        and np.all(strut_support == 13)
        and exact_rank_proof
        and collective_exact_error < mp.mpf("1e-100")
        and collective_numeric_error < mp.mpf("1e-13")
    )
    first_diagonal = min(diagonal_rows, key=lambda item: item[0])
    return {
        "matrix": matrix,
        "edge_by_row": tuple(edge_by_row),
        "first_diagonal": first_diagonal,
        "controls": controls,
        "controls_ok": controls_ok,
    }


def quadratic_form(hessian, carrier):
    h_sparse = sp.csr_matrix(hessian)
    g_sparse = sp.csr_matrix(carrier)
    result = (g_sparse.T @ (h_sparse @ g_sparse)).toarray()
    return np.asarray(result, dtype=float)


hashes = {name: sha256(path) for name, path in INPUTS.items()}
payloads = {
    name: json.loads(path.read_text())
    for name, path in INPUTS.items()
    if name not in {"geometry_source", "hessian_source", "carrier_source"}
}
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and payloads["finite_height"]["outcome"]
    == "SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB"
    and payloads["homothetic_action"]["outcome"]
    == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and payloads["carrier"]["tests"] == payloads["carrier"]["passed"] == 18
    and payloads["carrier_precision"]["outcome"]
    == "FULL_SCALE_STRUT_PRECISION_RESOLVED"
    and payloads["canonical_precision"]["tests"]
    == payloads["canonical_precision"]["passed"] == 17
    and payloads["symbolic_gap"]["outcome"]
    == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED"
    and payloads["schedule_census"]["outcome"]
    == "ORIENTATION_DOES_NOT_SELECT_PARITY"
)
check(
    "all preregistered artifacts and implementation libraries have frozen provenance",
    provenance_ok,
    str(hashes),
)

registered_scripts, registry_duplicates = registry_inventory(RUN_ALL)
registry_ok = bool(
    registered_scripts.count(VERIFIER_NAME) == 1
    and not registry_duplicates
)
check(
    "the verifier was registered exactly once before execution",
    registry_ok,
    f"registry={len(registered_scripts)}, duplicates={registry_duplicates}",
)


# Import only the already certified combinatorial model.  The large Hessian
# verifier and its old scientific state are loaded below by AST function
# extraction, so their bottom-level calculation is never executed.
spec = importlib.util.spec_from_file_location(
    "finite_height_quadratic_geometry", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
with contextlib.redirect_stdout(io.StringIO()):
    try:
        spec.loader.exec_module(gro)
    except SystemExit as upstream_exit:
        if upstream_exit.code not in (None, 0):
            raise
geometry_import_ok = bool(
    gro.tests == gro.passed == 43
    and set(gro.models) == {"even", "odd"}
)
check(
    "the complete one-slab combinatorial geometry retains all 43 certificates",
    geometry_import_ok,
)


schedule = payloads["schedule_census"]
orbits = schedule["h4_cover_action"]["order_orbits"]
schedule_ok = bool(
    schedule["schedule_census"]["orders"] == 120
    and schedule["schedule_census"]["unique_slabs"] == 120
    and [item["size"] for item in orbits] == [60, 60]
    and all(item["all_controls_pass"] for item in orbits)
    and schedule["h4_cover_action"]["distinct_induced_permutations"] == 60
    and schedule["h4_cover_action"]["kernel_size"] == 24
)
check(
    "the two representatives cover the two certified H4 schedule orbits of size 60",
    schedule_ok,
    f"orbit_sizes={[item['size'] for item in orbits]}",
)


state = finite_height_state()
committed_state = payloads["finite_height"]["history"]["slab_1"]
if state["bracket_ok"]:
    committed_error = max(
        abs(state["q"] - mp.mpf(committed_state["q"])),
        abs(state["h"] - mp.mpf(committed_state["h"])),
        abs(state["lambda"] - mp.mpf(committed_state["ratio"])),
    )
else:
    committed_error = mp.inf
root_ok = bool(
    state["bracket_ok"]
    and abs(state["residual"]) < mp.mpf("1e-90")
    and committed_error < mp.mpf("1e-70")
    and state["width"] < mp.mpf("1e-100")
)
check(
    "deterministic bisection reconstructs the frozen q in (9,10) independently",
    root_ok,
    (
        f"q={mp_text(state.get('q', mp.nan), 40)}, "
        f"residual={mp_text(state.get('residual', mp.inf), 6)}, "
        f"committed_error={mp_text(committed_error, 6)}"
    ),
)

domain_ok = bool(
    state["h"] > 0
    and state["lambda"] > 0
    and state["rho"] > 0
    and state["lambda"] != 1
    and state["tau_square"] > 0
    and state["exceptional_factor"] != 0
    and state["q_diag"] != 0
)
check(
    "the finite-height carrier lies on the complete generic real domain",
    domain_ok,
    (
        f"lambda={mp_text(state['lambda'], 30)}, "
        f"rho={mp_text(state['rho'], 30)}, "
        f"exceptional={mp_text(state['exceptional_factor'], 8)}"
    ),
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
}
load_named_functions(
    HESSIAN_SOURCE,
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
        "full_gradient_at_delta",
    },
    library,
)
models = {
    parity: library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}


runtime = {}
records = {}
for parity in ("even", "odd"):
    print(f"[{parity}] preparing the complete finite-height geometry", flush=True)
    model = models[parity]
    index_data = library["group_and_index_data"](
        model, (mp.log(state["lambda"]), mp.mpf(0))
    )
    geometry = library["prepare_geometry"](model, index_data)
    carrier = build_carrier(model, index_data, state["lambda"], state["rho"])

    geometry_ok = bool(
        len(model["slab"]) == 2400
        and len(model["old_edges"]) == OLD
        and len(model["internal_edges"]) == INTERNAL
        and len(model["new_edges"]) == NEW
        and len(index_data["edge_to_index"]) == FULL
        and len(geometry["triangle_records"]) == 6240
        and len(geometry["simplex_records"]) == 2400
        and len(geometry["patterns"]) == 20
        and all(
            index_data["edge_kind"][index] == "old"
            for index in range(OLD)
        )
        and all(
            index_data["edge_kind"][index] in {"internal", "pole"}
            for index in range(OLD, OLD + INTERNAL)
        )
        and all(
            index_data["edge_kind"][index] == "new"
            for index in range(OLD + INTERNAL, FULL)
        )
    )
    check(
        f"{parity}: the full 2280-edge geometry has the frozen active ordering",
        geometry_ok,
        f"patterns={len(geometry['patterns'])}, triangles={len(geometry['triangle_records'])}",
    )

    kind_values = {
        "old": mp.mpf(1),
        "internal": state["q_diag"],
        "pole": -state["rho"],
        "new": state["lambda"] ** 2,
    }
    print(f"[{parity}] differentiating 20 local simplex patterns", flush=True)
    pattern_cache, local_control = library["build_pattern_cache"](
        geometry["patterns"], kind_values
    )
    local_ok = bool(
        local_control["entry_pass"]
        and local_control["base_negative_counts"] == Counter({1: 2400})
        and local_control["all_displaced_have_one_negative"]
        and local_control["minimum_leading_minor"] > 0
        and local_control["minimum_argument"] > mp.mpf("1e-6")
    )
    check(
        f"{parity}: all base and displaced local simplices retain the Lorentzian branch",
        local_ok,
        (
            f"cross={mp_text(local_control['maximum_cross'], 6)}, "
            f"minor={mp_text(local_control['minimum_leading_minor'], 6)}, "
            f"argument={mp_text(local_control['minimum_argument'], 6)}"
        ),
    )

    direct_payload = {
        "signed_base": index_data["signed_base"],
        "triangle_records": geometry["triangle_records"],
        "simplex_records": geometry["simplex_records"],
        "pole_indices": geometry["pole_indices"],
        "directions": (),
    }
    library["_DIRECT"] = direct_payload
    direct_gradient, direct_branch = library["full_gradient_at_delta"](
        np.zeros(FULL)
    )

    print(f"[{parity}] assembling four complete 2280 x 2280 Hessians", flush=True)
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
    gradient_ok = bool(
        internal_maximum < mp.mpf("1e-25")
        and assembly_gradient_error < 2e-11
        and direct_branch["negative_counts"] == Counter({1: 2400})
        and direct_branch["minimum_leading_minor"] > 0
        and direct_branch["minimum_argument"] > mp.mpf("1e-6")
        and direct_branch["maximum_imaginary"] < mp.mpf("1e-60")
        and assembly["maximum_imaginary"] < 1e-60
    )
    check(
        f"{parity}: all artificial equations vanish and the direct gradient is real",
        gradient_ok,
        (
            f"internal={mp_text(internal_maximum, 6)}, "
            f"assembly={assembly_gradient_error:.3e}, "
            f"imaginary={mp_text(direct_branch['maximum_imaginary'], 6)}"
        ),
    )

    check(
        f"{parity}: the geometry-selected scale-plus-strut carrier passes exact controls",
        carrier["controls_ok"],
        (
            f"support={carrier['controls']['row_support']}, "
            f"rank240={carrier['controls']['exact_rank_240']}"
        ),
    )

    operational = hessians["operational_primary"]
    shadow = hessians["operational_shadow"]
    validation = hessians["validation_primary"]
    validation_shadow = hessians["validation_shadow"]
    derivative_frobenius = (
        frobenius_difference(operational, shadow)
        + frobenius_difference(validation, validation_shadow)
        + frobenius_difference(operational, validation)
    )
    round_entry = assembly["binary64_roundoff"]["maximum_entry_bound"]
    round_antisymmetry = 2 * FULL * round_entry
    reciprocity_allowance = 10 * (derivative_frobenius + round_antisymmetry)
    antisymmetry = {
        name: frobenius_difference(matrix, matrix.T)
        for name, matrix in hessians.items()
    }
    cross = np.abs(operational - validation)
    proxy = (
        np.abs(operational - shadow)
        + np.abs(validation - validation_shadow)
        + 1e-70
    )
    entrywise_ok = bool(np.all(cross <= 10 * proxy))
    del cross, proxy
    reciprocity_ok = bool(
        entrywise_ok
        and max(antisymmetry.values()) <= reciprocity_allowance
    )

    active_hessians = {
        name: matrix[OLD:, OLD:]
        for name, matrix in hessians.items()
    }
    forms = {
        name: quadratic_form(matrix, carrier["matrix"])
        for name, matrix in active_hessians.items()
    }
    symmetric_forms = {
        name: (matrix + matrix.T) / 2
        for name, matrix in forms.items()
    }
    form_antisymmetry = {
        name: norm2(matrix - matrix.T)
        for name, matrix in forms.items()
    }

    permutation = np.asarray(sorted(
        range(ACTIVE), key=lambda row: carrier["edge_by_row"][row]
    ), dtype=int)
    reordered_hessian = active_hessians["operational_primary"][
        np.ix_(permutation, permutation)
    ]
    reordered_carrier = carrier["matrix"][permutation, :]
    reordered_form = quadratic_form(reordered_hessian, reordered_carrier)
    reordered_form = (reordered_form + reordered_form.T) / 2
    reorder_absolute_error = norm2(
        reordered_form - symmetric_forms["operational_primary"]
    )
    del reordered_hessian, reordered_carrier, reordered_form

    reciprocity_and_reorder_precontrol = bool(
        reciprocity_ok
        and all(math.isfinite(value) for value in form_antisymmetry.values())
        and math.isfinite(reorder_absolute_error)
    )
    check(
        f"{parity}: Hessian reciprocity and physical-edge reordering controls are finite",
        reciprocity_and_reorder_precontrol,
        (
            f"antisym_max={max(antisymmetry.values()):.3e}, "
            f"allowance={reciprocity_allowance:.3e}, "
            f"reorder_abs={reorder_absolute_error:.3e}"
        ),
    )

    old_gradient = {}
    new_gradient = {}
    for edge, index in index_data["edge_to_index"].items():
        kind = index_data["edge_kind"][index]
        if kind == "old":
            old_gradient[tuple(edge)] = mp.re(direct_gradient[index])
        elif kind == "new":
            logical = tuple(value - VERTICES for value in edge)
            new_gradient[logical] = mp.re(direct_gradient[index])

    carrier_norm = norm2(carrier["matrix"])
    first_edge, first_row, first_source = carrier["first_diagonal"]
    corrupted_carrier = carrier["matrix"].copy()
    corrupted_carrier[first_row, first_source] += 0.1
    corrupted_form = quadratic_form(
        active_hessians["operational_primary"], corrupted_carrier
    )
    corrupted_form = (corrupted_form + corrupted_form.T) / 2
    carrier_corruption_absolute = norm2(
        corrupted_form - symmetric_forms["operational_primary"]
    )
    del corrupted_carrier, corrupted_form

    records[parity] = {
        "geometry": {
            "simplices": len(model["slab"]),
            "edges": len(index_data["edge_to_index"]),
            "triangles": len(geometry["triangle_records"]),
            "patterns": len(geometry["patterns"]),
        },
        "branch": {
            "minimum_leading_minor": mp_text(
                local_control["minimum_leading_minor"]
            ),
            "minimum_argument": mp_text(local_control["minimum_argument"]),
            "maximum_imaginary": mp_text(
                direct_branch["maximum_imaginary"]
            ),
        },
        "gradient": {
            "internal_maximum": mp_text(internal_maximum),
            "assembly_error": float_text(assembly_gradient_error),
        },
        "carrier": carrier["controls"],
        "hessian": {
            "derivative_frobenius": float_text(derivative_frobenius),
            "round_entry_bound": float_text(round_entry),
            "reciprocity_allowance": float_text(reciprocity_allowance),
            "antisymmetry_frobenius": {
                name: float_text(value)
                for name, value in antisymmetry.items()
            },
            "form_antisymmetry_two_norm": {
                name: float_text(value)
                for name, value in form_antisymmetry.items()
            },
            "reorder_absolute_error": float_text(reorder_absolute_error),
        },
        "carrier_norm": float_text(carrier_norm),
        "carrier_corruption": {
            "edge": list(first_edge),
            "row": int(first_row),
            "source_column": int(first_source),
            "absolute_form_change": float_text(carrier_corruption_absolute),
        },
    }
    runtime[parity] = {
        "forms": symmetric_forms,
        "gradient_old": old_gradient,
        "gradient_new": new_gradient,
        "carrier": carrier["matrix"],
        "carrier_norm": carrier_norm,
        "round_entry": round_entry,
        "reorder_absolute": reorder_absolute_error,
        "carrier_corruption_absolute": carrier_corruption_absolute,
        "control_ok": bool(
            geometry_ok
            and local_ok
            and gradient_ok
            and carrier["controls_ok"]
            and reciprocity_and_reorder_precontrol
        ),
    }

    del hessians, active_hessians, forms, pattern_cache
    gc.collect()


old_keys_match = set(runtime["even"]["gradient_old"]) == set(
    runtime["odd"]["gradient_old"]
)
new_keys_match = set(runtime["even"]["gradient_new"]) == set(
    runtime["odd"]["gradient_new"]
)
old_gradient_difference = (
    max(
        abs(runtime["even"]["gradient_old"][edge]
            - runtime["odd"]["gradient_old"][edge])
        for edge in runtime["even"]["gradient_old"]
    )
    if old_keys_match else mp.inf
)
new_gradient_difference = (
    max(
        abs(runtime["even"]["gradient_new"][edge]
            - runtime["odd"]["gradient_new"][edge])
        for edge in runtime["even"]["gradient_new"]
    )
    if new_keys_match else mp.inf
)
boundary_gradient_ok = bool(
    old_keys_match
    and new_keys_match
    and old_gradient_difference < mp.mpf("1e-25")
    and new_gradient_difference < mp.mpf("1e-25")
)
check(
    "the even and odd first gradients agree edgewise on both physical boundaries",
    boundary_gradient_ok,
    (
        f"old={mp_text(old_gradient_difference, 6)}, "
        f"new={mp_text(new_gradient_difference, 6)}"
    ),
)


scheme_names = tuple(DERIVATIVE_STEPS)
op_name = "operational_primary"
normalization = max(
    1.0,
    norm2(runtime["even"]["forms"][op_name]),
    norm2(runtime["odd"]["forms"][op_name]),
)
difference_matrices = {
    name: runtime["even"]["forms"][name]
    - runtime["odd"]["forms"][name]
    for name in scheme_names
}
normalized_differences = {
    name: norm2(matrix) / normalization
    for name, matrix in difference_matrices.items()
}
e_step = max(
    norm2(runtime[parity]["forms"][name]
          - runtime[parity]["forms"][op_name]) / normalization
    for parity in ("even", "odd")
    for name in scheme_names
)
e_round = (
    ACTIVE * runtime["even"]["round_entry"]
    * runtime["even"]["carrier_norm"] ** 2
    + ACTIVE * runtime["odd"]["round_entry"]
    * runtime["odd"]["carrier_norm"] ** 2
) / normalization
e_total = e_step + e_round + 100 * np.finfo(float).eps
pairwise_difference_agreement = max(
    norm2(difference_matrices[left] - difference_matrices[right])
    / normalization
    for index, left in enumerate(scheme_names)
    for right in scheme_names[index + 1:]
)
numerical_error_ok = bool(
    all(math.isfinite(value) for value in normalized_differences.values())
    and math.isfinite(e_step)
    and math.isfinite(e_round)
    and math.isfinite(e_total)
    and e_total > 0
)
check(
    "the preregistered target-independent error budget is finite",
    numerical_error_ok,
    (
        f"e_step={e_step:.3e}, e_round={e_round:.3e}, "
        f"e_total={e_total:.3e}"
    ),
)


reorder_relative = max(
    runtime[parity]["reorder_absolute"] / normalization
    for parity in ("even", "odd")
)
reorder_ok = reorder_relative <= 10 * e_total

uniform = np.zeros((DATA, 2), dtype=float)
uniform[:VERTICES, 0] = 1 / math.sqrt(VERTICES)
uniform[VERTICES:, 1] = 1 / math.sqrt(VERTICES)
homogeneous_relative = norm2(
    uniform.T @ difference_matrices[op_name] @ uniform
) / normalization
homogeneous_ok = homogeneous_relative <= 10 * e_total
check(
    "the known homogeneous subdivision and edge-reordering controls agree",
    homogeneous_ok and reorder_ok,
    (
        f"homogeneous={homogeneous_relative:.3e}, "
        f"reorder={reorder_relative:.3e}, gate={10*e_total:.3e}"
    ),
)


odd_carrier = runtime["odd"]["carrier"]
synthetic_u = odd_carrier[:, 0] / np.linalg.norm(odd_carrier[:, 0])
synthetic_image = odd_carrier.T @ synthetic_u
synthetic_relative = (
    float(np.dot(synthetic_image, synthetic_image)) / normalization
)
synthetic_ok = synthetic_relative > 100 * e_total
check(
    "the frozen symmetric rank-one Hessian corruption is resolved",
    synthetic_ok,
    f"effect={synthetic_relative:.3e}, gate={100*e_total:.3e}",
)


carrier_corruption_relative = (
    runtime["odd"]["carrier_corruption_absolute"] / normalization
)
carrier_corruption_ok = carrier_corruption_relative > 100 * e_total
check(
    "the frozen +1/10 incidence-coefficient corruption is resolved",
    carrier_corruption_ok,
    f"effect={carrier_corruption_relative:.3e}, gate={100*e_total:.3e}",
)


all_controls = bool(
    provenance_ok
    and registry_ok
    and geometry_import_ok
    and schedule_ok
    and root_ok
    and domain_ok
    and all(runtime[parity]["control_ok"] for parity in ("even", "odd"))
    and boundary_gradient_ok
    and numerical_error_ok
    and homogeneous_ok
    and reorder_ok
    and synthetic_ok
    and carrier_corruption_ok
)

maximum_difference = max(normalized_differences.values())
minimum_difference = min(normalized_differences.values())
if not all_controls:
    outcome = "FINITE_HEIGHT_QUADRATIC_CONTROL_FAILED"
elif maximum_difference <= 10 * e_total:
    outcome = "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY"
elif (
    minimum_difference > 100 * e_total
    and pairwise_difference_agreement <= 10 * e_total
):
    outcome = "FINITE_HEIGHT_QUADRATIC_PARITY_DEPENDENT_PRIMARY"
else:
    outcome = "FINITE_HEIGHT_QUADRATIC_PARITY_OPEN"

outcome_ok = outcome in {
    "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY",
    "FINITE_HEIGHT_QUADRATIC_PARITY_DEPENDENT_PRIMARY",
    "FINITE_HEIGHT_QUADRATIC_PARITY_OPEN",
    "FINITE_HEIGHT_QUADRATIC_CONTROL_FAILED",
}
check(
    "the primary outcome follows the preregistered hierarchy without a spectrum",
    outcome_ok,
    (
        f"outcome={outcome}, min_d={minimum_difference:.3e}, "
        f"max_d={maximum_difference:.3e}, agreement={pairwise_difference_agreement:.3e}"
    ),
)


matrix_stack = np.stack([
    runtime["even"]["forms"][op_name],
    runtime["odd"]["forms"][op_name],
    difference_matrices[op_name],
])
np.save(MATRIX_OUTPUT, matrix_stack, allow_pickle=False)
matrix_hash = sha256(MATRIX_OUTPUT)

artifact = {
    "outcome": outcome,
    "status": "PRIMARY_ONLY_PENDING_MECHANICALLY_DIFFERENT_ADVERSARIAL_REPLICATION",
    "tests": tests,
    "passed": passed,
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "schedule_provenance_commit": SCHEDULE_PROVENANCE_COMMIT,
        "complete_provenance_commit": COMPLETE_PROVENANCE_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": hashes,
        "source_sha256": sha256(Path(__file__)),
    },
    "background": {
        "v": "3/2",
        "q": mp_text(state["q"]),
        "h": mp_text(state["h"]),
        "lambda": mp_text(state["lambda"]),
        "rho": mp_text(state["rho"]),
        "mass": mp_text(state["mass"]),
        "elimination_residual": mp_text(state["residual"]),
        "root_bracket_width": mp_text(state["width"]),
        "q_diagonal": mp_text(state["q_diag"]),
        "tau_square": mp_text(state["tau_square"]),
        "exceptional_factor": mp_text(state["exceptional_factor"]),
        "committed_reconstruction_error": mp_text(committed_error),
    },
    "schedule_scope": {
        "orders": 120,
        "h4_orbit_sizes": [item["size"] for item in orbits],
        "representatives": ["even", "odd"],
    },
    "first_gradient_covariance": {
        "old_boundary_maximum_difference": mp_text(old_gradient_difference),
        "new_boundary_maximum_difference": mp_text(new_gradient_difference),
    },
    "parities": records,
    "comparison": {
        "normalization": float_text(normalization),
        "normalized_difference_two_norms": {
            name: float_text(value)
            for name, value in normalized_differences.items()
        },
        "pairwise_difference_agreement": float_text(
            pairwise_difference_agreement
        ),
        "e_step": float_text(e_step),
        "e_round": float_text(e_round),
        "e_total": float_text(e_total),
        "independent_gate": float_text(10 * e_total),
        "dependent_gate": float_text(100 * e_total),
        "homogeneous_relative_difference": float_text(homogeneous_relative),
        "reorder_relative_difference": float_text(reorder_relative),
        "synthetic_rank_one_relative_effect": float_text(synthetic_relative),
        "carrier_corruption_relative_effect": float_text(
            carrier_corruption_relative
        ),
    },
    "matrices": {
        "path": MATRIX_OUTPUT.name,
        "sha256": matrix_hash,
        "shape": list(matrix_stack.shape),
        "order": ["even_operational", "odd_operational", "even_minus_odd"],
    },
    "firewall": {
        "quadratic_spectrum_retained_or_inspected": False,
        "matrix_two_norms_used_for_error_classification_only": True,
        "continuum_mode_parsed": False,
        "speed_target_parsed": False,
        "desired_eigenvalue_parsed": False,
        "full_suite_run": False,
    },
    "interpretation": {
        "primary_only": True,
        "nonlinear_integrability": "OPEN",
        "boundary_evolution": "NOT_TESTED",
        "gravitons": "NOT_TESTED",
        "limiting_speed": "NOT_TESTED",
        "c_G_Planck_scales": "NOT_DERIVED",
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("\n" + "=" * 78)
print(f"OUTCOME: {outcome}")
print(f"MATRIX SHA-256: {matrix_hash}")
print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)
