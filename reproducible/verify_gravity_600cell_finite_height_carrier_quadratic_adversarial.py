#!/usr/bin/env python3
"""High-precision orbit-kernel replication of finite-height parity equality.

Primary result commit: d2796de.
Adversarial protocol commit: 0235162.
Registry commit: 6fde152.

The complete binary64 ambient Hessian assembler is deliberately not loaded.
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

import mpmath as mp
import networkx as nx
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_adversarial.json"
)
MATRIX_OUTPUT = (
    HERE
    / "gravity_600cell_finite_height_carrier_quadratic_adversarial_matrices.npy"
)
PRIMARY_JSON = HERE / "gravity_600cell_finite_height_carrier_quadratic.json"
PRIMARY_MATRICES = (
    HERE / "gravity_600cell_finite_height_carrier_quadratic_matrices.npy"
)
PRIMARY_SOURCE = (
    HERE / "verify_gravity_600cell_finite_height_carrier_quadratic.py"
)
HIGH_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = (
    HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
)
SYMBOLIC_INPUT = (
    HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
)
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "primary_json": PRIMARY_JSON,
    "primary_matrices": PRIMARY_MATRICES,
    "primary_source": PRIMARY_SOURCE,
    "high_source": HIGH_SOURCE,
    "geometry_source": GEOMETRY_SOURCE,
    "symbolic_input": SYMBOLIC_INPUT,
}
EXPECTED_HASHES = {
    "primary_json": (
        "0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9"
    ),
    "primary_matrices": (
        "e01bfb28d4c5313b466118315f8ca22c16c2cdc4e94ab05f30c730a136d81cb2"
    ),
    "primary_source": (
        "bbe7112270a7f2bcb2d443fab45ca450598e7234250bd335b14a4ed7869443a5"
    ),
    "high_source": (
        "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571"
    ),
    "geometry_source": (
        "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf"
    ),
    "symbolic_input": (
        "ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179"
    ),
}

PRIMARY_RESULT_COMMIT = "d2796de"
PROTOCOL_COMMIT = "0235162"
REGISTRY_COMMIT = "6fde152"
VERIFIER_NAME = Path(__file__).name

DPS = 180
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-25"),
    "operational_shadow": mp.mpf("5e-26"),
    "validation_primary": mp.mpf("2.5e-26"),
    "validation_shadow": mp.mpf("1.25e-26"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
ARITHMETIC_FLOOR = mp.mpf("1e-150")
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
        if any(isinstance(target, ast.Name) and target.id == "scripts"
               for target in node.targets):
            scripts = ast.literal_eval(node.value)
            break
    if scripts is None:
        raise RuntimeError("run_all.py has no literal scripts registry")
    counts = Counter(scripts)
    duplicates = sorted(name for name, count in counts.items() if count != 1)
    return scripts, duplicates


def mp_frobenius(matrix):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in matrix))


def mp_symmetric(matrix):
    rows, columns = matrix.rows, matrix.cols
    result = mp.matrix(rows, columns)
    for row in range(rows):
        for column in range(columns):
            result[row, column] = (
                matrix[row, column] + matrix[column, row]
            ) / 2
    return result


def mp_difference(left, right):
    return mp.matrix([
        [left[row, column] - right[row, column]
         for column in range(left.cols)]
        for row in range(left.rows)
    ])


def mp_linear_combination(left_factor, left, right_factor, right):
    return mp.matrix([
        [
            left_factor * left[row, column]
            + right_factor * right[row, column]
            for column in range(left.cols)
        ]
        for row in range(left.rows)
    ])


def richardson(coarse, fine):
    return mp_linear_combination(mp.mpf(-1) / 3, coarse, mp.mpf(4) / 3, fine)


def quadratic_value(matrix, vector):
    return mp.fsum(
        vector[row] * matrix[row, column] * vector[column]
        for row in range(DATA)
        for column in range(DATA)
        if vector[row] != 0 and vector[column] != 0
    )


def mp_to_numpy(matrix):
    return np.asarray([
        [float(mp.re(matrix[row, column])) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.float64)


def finite_height_formula_control(primary):
    q = mp.mpf(primary["background"]["q"])
    h = mp.mpf(primary["background"]["h"])
    lam = mp.mpf(primary["background"]["lambda"])
    rho = mp.mpf(primary["background"]["rho"])
    mass = mp.mpf(primary["background"]["mass"])
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

    elimination = (
        4 * mp.pi * (mu(q) - mu(v))
        + q * (momentum(q) - momentum(v))
    )
    reconstructed_h = (
        momentum(q) - momentum(v)
    ) / (2 * mp.pi * mu(q))
    return {
        "q": q,
        "h": h,
        "lambda": lam,
        "rho": rho,
        "mass": mass,
        "elimination": elimination,
        "h_error": abs(h - reconstructed_h),
        "lambda_error": abs(lam - (1 + h * q)),
        "rho_error": abs(rho - h**2),
        "mass_error": abs(mass - mu(v)),
    }


def build_sparse_carrier(model, index_data, lam, rho):
    q_diag = lam - rho
    A = -16 * rho / (lam - 1) ** 2
    B = 8 + 16 * rho / (lam - 1) ** 2
    a_value = A / (8 * q_diag)
    b_value = B / (8 * q_diag)
    kappa = rho / ((lam - 1) * q_diag)

    rows = [dict() for _ in range(ACTIVE)]
    edge_by_row = [None] * ACTIVE
    poles = {}
    first_diagonal = None
    for edge, global_index in index_data["edge_to_index"].items():
        if global_index < OLD:
            continue
        edge = tuple(int(value) for value in edge)
        row = global_index - OLD
        edge_by_row[row] = edge
        kind = index_data["edge_kind"][global_index]
        if kind == "pole":
            lower, upper = edge
            if upper != lower + VERTICES:
                raise RuntimeError("invalid pole endpoint convention")
            rows[row][VERTICES + lower] = mp.mpf(1)
            poles[lower] = row
        elif kind == "internal":
            lower, upper = edge
            target = upper - VERTICES
            rows[row] = {
                lower: a_value,
                target: b_value,
                VERTICES + lower: kappa,
                VERTICES + target: -lam * kappa,
            }
            candidate = (edge, row, lower)
            if first_diagonal is None or candidate[0] < first_diagonal[0]:
                first_diagonal = candidate
        elif kind == "new":
            left, right = edge[0] - VERTICES, edge[1] - VERTICES
            rows[row] = {left: 1 / lam, right: 1 / lam}
        else:
            raise RuntimeError(f"unexpected active edge kind {kind}")

    support = Counter(len(row) for row in rows)
    scale_support = Counter()
    strut_support = Counter()
    for row in rows:
        for column in row:
            if column < VERTICES:
                scale_support[column] += 1
            else:
                strut_support[column - VERTICES] += 1
    pole_identity = bool(
        len(poles) == VERTICES
        and all(rows[poles[vertex]] == {VERTICES + vertex: mp.mpf(1)}
                for vertex in range(VERTICES))
    )

    graph = nx.Graph()
    graph.add_nodes_from(range(VERTICES))
    for edge in model["new_edges"]:
        graph.add_edge(edge[0] - VERTICES, edge[1] - VERTICES)
    graph_ok = bool(
        graph.number_of_nodes() == VERTICES
        and graph.number_of_edges() == NEW
        and nx.is_connected(graph)
        and not nx.is_bipartite(graph)
    )

    algebra_error = max(
        abs(a_value + b_value - 1 / q_diag),
        abs(kappa - lam * kappa + rho / q_diag),
    )
    collective_error = mp.mpf(0)
    for row_index, row in enumerate(rows):
        kind = index_data["edge_kind"][OLD + row_index]
        scale_sum = mp.fsum(value for column, value in row.items()
                            if column < VERTICES)
        strut_sum = mp.fsum(value for column, value in row.items()
                            if column >= VERTICES)
        if kind == "pole":
            expected = (mp.mpf(0), mp.mpf(1))
        elif kind == "internal":
            expected = (1 / q_diag, -rho / q_diag)
        else:
            expected = (2 / lam, mp.mpf(0))
        collective_error = max(
            collective_error,
            abs(scale_sum - expected[0]),
            abs(strut_sum - expected[1]),
        )

    controls_ok = bool(
        support == Counter({1: 120, 2: 720, 4: 720})
        and len(scale_support) == VERTICES
        and len(strut_support) == VERTICES
        and all(scale_support[vertex] == 24 for vertex in range(VERTICES))
        and all(strut_support[vertex] == 13 for vertex in range(VERTICES))
        and pole_identity
        and graph_ok
        and algebra_error < mp.mpf("1e-160")
        and collective_error < mp.mpf("1e-160")
    )
    return {
        "rows": rows,
        "edge_by_row": tuple(edge_by_row),
        "first_diagonal": first_diagonal,
        "controls_ok": controls_ok,
        "controls": {
            "support": dict(sorted(support.items())),
            "scale_support": [
                min(scale_support.values()), max(scale_support.values())
            ],
            "strut_support": [
                min(strut_support.values()), max(strut_support.values())
            ],
            "pole_identity": pole_identity,
            "connected_nonbipartite_upper_graph": graph_ok,
            "exact_rank_240": pole_identity and graph_ok,
            "algebra_error": mp_text(algebra_error),
            "collective_error": mp_text(collective_error),
            "coefficients": {
                "a": mp_text(a_value),
                "b": mp_text(b_value),
                "kappa": mp_text(kappa),
            },
        },
    }


def pull_orbit_kernel(kernel, index_data, carrier_rows, reverse=False):
    table = index_data["table"]
    result = mp.matrix(DATA, DATA)
    maximum_imaginary = mp.mpf(0)
    active_keys = 0
    expanded_entries = 0
    for (row_type, column_type, relative), value in kernel.items():
        if row_type < 30 or column_type < 30:
            continue
        active_keys += 1
        maximum_imaginary = max(maximum_imaginary, abs(mp.im(value)))
        real_value = mp.re(value)
        for row_group in range(24):
            if reverse:
                column_group = int(table[relative, row_group])
            else:
                column_group = int(table[row_group, relative])
            row = 24 * (row_type - 30) + row_group
            column = 24 * (column_type - 30) + column_group
            left = carrier_rows[row]
            right = carrier_rows[column]
            if not left or not right or real_value == 0:
                continue
            expanded_entries += 1
            for left_column, left_value in left.items():
                for right_column, right_value in right.items():
                    result[left_column, right_column] += (
                        real_value * left_value * right_value
                    )
    return result, {
        "active_kernel_keys": active_keys,
        "expanded_nonzero_hessian_entries": expanded_entries,
        "maximum_kernel_imaginary": maximum_imaginary,
    }


def combine_kernels(coarse, fine):
    keys = set(coarse) | set(fine)
    return {
        key: (4 * fine.get(key, 0) - coarse.get(key, 0)) / 3
        for key in keys
    }


def action_at_active_delta(
    model, index_data, geometry, angle_record, area_data, mass, delta
):
    signed = list(index_data["signed_base"])
    for row, value in enumerate(delta):
        signed[OLD + row] *= mp.exp(value)

    curvature = [
        mp.pi if record["boundary"] else 2 * mp.pi
        for record in geometry["triangle_records"]
    ]
    negative_counts = Counter()
    minimum_minor = mp.inf
    minimum_argument = mp.inf
    for simplex in geometry["simplex_records"]:
        values = [signed[index] for index in simplex["indices"]]
        angles, branch = angle_record(values)
        negative_counts[branch["negative_directions"]] += 1
        minimum_minor = min(
            minimum_minor, branch["minimum_leading_minor"]
        )
        minimum_argument = min(minimum_argument, branch["minimum_argument"])
        for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
            curvature[triangle] += angles[hinge_index]

    action = mp.mpc(0)
    for triangle_index, record in enumerate(geometry["triangle_records"]):
        values = [signed[index] for index in record["indices"]]
        area, _, _ = area_data(values)
        action += -I * area * curvature[triangle_index]
    for pole_index in geometry["pole_indices"]:
        action += -(8 * mp.pi * mass / VERTICES) * mp.sqrt(-signed[pole_index])
    return action, {
        "negative_counts": negative_counts,
        "minimum_leading_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": abs(mp.im(action)),
    }


def carrier_direction(rows, vector):
    return [
        mp.fsum(value * vector[column] for column, value in row.items())
        for row in rows
    ]


hashes = {name: sha256(path) for name, path in INPUTS.items()}
primary = json.loads(PRIMARY_JSON.read_text())
symbolic = json.loads(SYMBOLIC_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and primary["outcome"]
    == "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY"
    and primary["tests"] == primary["passed"] == 22
    and symbolic["outcome"] == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED"
)
check(
    "the primary result and independent implementation sources have frozen provenance",
    provenance_ok,
    str(hashes),
)

registered_scripts, registry_duplicates = registry_inventory(RUN_ALL)
registry_ok = bool(
    registered_scripts.count(VERIFIER_NAME) == 1
    and not registry_duplicates
)
check(
    "the adversarial verifier was registered exactly once before execution",
    registry_ok,
    f"registry={len(registered_scripts)}, duplicates={registry_duplicates}",
)


spec = importlib.util.spec_from_file_location(
    "finite_height_quadratic_adversarial_geometry", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
with contextlib.redirect_stdout(io.StringIO()):
    try:
        spec.loader.exec_module(gro)
    except SystemExit as upstream_exit:
        if upstream_exit.code not in (None, 0):
            raise
geometry_ok = bool(gro.tests == gro.passed == 43 and set(gro.models) == {"even", "odd"})
check(
    "the independently imported complete geometry retains all 43 certificates",
    geometry_ok,
)


state = finite_height_formula_control(primary)
background_ok = bool(
    abs(state["elimination"]) < mp.mpf("1e-70")
    and state["h_error"] < mp.mpf("1e-70")
    and state["lambda_error"] < mp.mpf("1e-70")
    and state["rho_error"] < mp.mpf("1e-70")
    and state["mass_error"] < mp.mpf("1e-70")
    and state["h"] > 0
    and state["lambda"] > 0
    and state["rho"] > 0
)
check(
    "the primary background independently satisfies the exact finite-height equations",
    background_ok,
    (
        f"E={mp_text(state['elimination'], 6)}, "
        f"h_error={mp_text(state['h_error'], 6)}, "
        f"lambda_error={mp_text(state['lambda_error'], 6)}"
    ),
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
    "combinations": combinations,
    "defaultdict": defaultdict,
    "gro": gro,
    "mp": mp,
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
    HIGH_SOURCE,
    {"high_precision_pattern_cache", "assemble_full_representative_kernels"},
    library,
)
models = {
    parity: library["augment_boundary_orbits"](model)
    for parity, model in gro.models.items()
}


records = {}
runtime = {}
for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing the high-precision orbit kernel", flush=True)
    model = models[parity]
    index_data = library["group_and_index_data"](
        model, (mp.log(state["lambda"]), mp.mpf(0))
    )
    geometry = library["prepare_geometry"](model, index_data)
    carrier = build_sparse_carrier(
        model, index_data, state["lambda"], state["rho"]
    )
    carrier_and_geometry_ok = bool(
        len(model["slab"]) == 2400
        and len(index_data["edge_to_index"]) == FULL
        and len(geometry["patterns"]) == 20
        and carrier["controls_ok"]
    )
    check(
        f"{parity}: sparse carrier and orbit geometry pass independent exact controls",
        carrier_and_geometry_ok,
        f"support={carrier['controls']['support']}, rank240={carrier['controls']['exact_rank_240']}",
    )

    kind_values = {
        "old": mp.mpf(1),
        "internal": state["lambda"] - state["rho"],
        "pole": -state["rho"],
        "new": state["lambda"] ** 2,
    }
    pattern_cache, pattern_control = library["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    pattern_ok = bool(
        pattern_control["entry_pass"]
        and pattern_control["base_negative_counts"] == Counter({1: 2400})
        and set(pattern_control["displaced_negative_counts"]) == {1}
        and pattern_control["minimum_leading_minor"] > 0
        and pattern_control["minimum_argument"] > mp.mpf("1e-6")
        and pattern_control["maximum_raw_angle_or_derivative_imaginary"]
        < mp.mpf("1e-140")
    )
    check(
        f"{parity}: the four 180-digit derivative levels retain the Lorentzian branch",
        pattern_ok,
        (
            f"cross={mp_text(pattern_control['maximum_cross'], 6)}, "
            f"proxy={mp_text(pattern_control['maximum_proxy'], 6)}, "
            f"imag={mp_text(pattern_control['maximum_raw_angle_or_derivative_imaginary'], 6)}"
        ),
    )

    kernels, kernel_control = library["assemble_full_representative_kernels"](
        index_data, geometry, pattern_cache
    )
    kernel_ok = bool(
        set(kernels) == set(VARIANTS)
        and all(len(kernel) > 0 for kernel in kernels.values())
        and kernel_control["maximum_imaginary"] < mp.mpf("1e-140")
    )
    check(
        f"{parity}: identity-row 2T kernels assemble with negligible imaginary residue",
        kernel_ok,
        (
            f"nonzero={kernel_control['nonzero_entries']}, "
            f"imag={mp_text(kernel_control['maximum_imaginary'], 6)}"
        ),
    )

    raw_forms = {}
    pull_controls = {}
    for name in VARIANTS:
        raw, pull_control = pull_orbit_kernel(
            kernels[name], index_data, carrier["rows"]
        )
        raw_forms[name] = mp_symmetric(raw)
        pull_control["antisymmetry_frobenius"] = mp_text(
            mp_frobenius(mp_difference(raw, mp.matrix(raw).T))
        )
        pull_controls[name] = pull_control

    R01 = richardson(
        raw_forms["operational_primary"],
        raw_forms["operational_shadow"],
    )
    R12 = richardson(
        raw_forms["operational_shadow"],
        raw_forms["validation_primary"],
    )
    R23 = richardson(
        raw_forms["validation_primary"],
        raw_forms["validation_shadow"],
    )
    richardson_finite = all(
        mp.isfinite(value)
        for matrix in (R01, R12, R23)
        for value in matrix
    )
    check(
        f"{parity}: direct orbit convolution yields three finite Richardson forms",
        richardson_finite,
    )

    records[parity] = {
        "carrier": carrier["controls"],
        "pattern": {
            "minimum_leading_minor": mp_text(
                pattern_control["minimum_leading_minor"]
            ),
            "minimum_argument": mp_text(
                pattern_control["minimum_argument"]
            ),
            "maximum_cross": mp_text(pattern_control["maximum_cross"]),
            "maximum_proxy": mp_text(pattern_control["maximum_proxy"]),
            "maximum_imaginary": mp_text(
                pattern_control["maximum_raw_angle_or_derivative_imaginary"]
            ),
        },
        "kernel": {
            "nonzero_entries": kernel_control["nonzero_entries"],
            "maximum_imaginary": mp_text(
                kernel_control["maximum_imaginary"]
            ),
            "pull_controls": pull_controls,
        },
    }
    runtime[parity] = {
        "model": model,
        "index_data": index_data,
        "geometry": geometry,
        "carrier": carrier,
        "kernels": kernels,
        "forms": {"R01": R01, "R12": R12, "R23": R23},
        "control_ok": bool(
            carrier_and_geometry_ok and pattern_ok and kernel_ok
            and richardson_finite
        ),
    }


normalization = max(
    mp.mpf(1),
    mp_frobenius(runtime["even"]["forms"]["R12"]),
    mp_frobenius(runtime["odd"]["forms"]["R12"]),
)
differences = {
    level: mp_difference(
        runtime["even"]["forms"][level],
        runtime["odd"]["forms"][level],
    )
    for level in ("R01", "R12", "R23")
}
d_values = {
    level: mp_frobenius(matrix) / normalization
    for level, matrix in differences.items()
}
e_step = max(
    mp_frobenius(mp_difference(
        runtime[parity]["forms"]["R01"],
        runtime[parity]["forms"]["R12"],
    )) / normalization
    for parity in ("even", "odd")
)
e_step = max(
    e_step,
    max(
        mp_frobenius(mp_difference(
            runtime[parity]["forms"]["R12"],
            runtime[parity]["forms"]["R23"],
        )) / normalization
        for parity in ("even", "odd")
    ),
)
e_total = e_step + mp.mpf("1e-145")
difference_agreement = max(
    mp_frobenius(mp_difference(differences[left], differences[right]))
    / normalization
    for index, left in enumerate(("R01", "R12", "R23"))
    for right in ("R01", "R12", "R23")[index + 1:]
)
high_precision_independent = max(d_values.values()) <= 10 * e_total
high_precision_dependent = bool(
    min(d_values.values()) > 100 * e_total
    and difference_agreement <= 10 * e_total
)
classification_finite = bool(
    all(mp.isfinite(value) for value in d_values.values())
    and mp.isfinite(e_step)
    and e_total > 0
)
check(
    "the target-free high-precision parity classification is resolved or honestly open",
    classification_finite,
    (
        f"d={[mp_text(d_values[x], 8) for x in ('R01','R12','R23')]}, "
        f"e_step={mp_text(e_step, 8)}, e_total={mp_text(e_total, 8)}"
    ),
)


# Carrier corruption: rebuild only the finest odd Richardson kernel with the
# preregistered single coefficient change.
odd = runtime["odd"]
corrupted_rows = [dict(row) for row in odd["carrier"]["rows"]]
first_edge, first_row, first_source = odd["carrier"]["first_diagonal"]
corrupted_rows[first_row][first_source] += mp.mpf(1) / 10
kernel_R12_odd = combine_kernels(
    odd["kernels"]["operational_shadow"],
    odd["kernels"]["validation_primary"],
)
corrupted_R12, _ = pull_orbit_kernel(
    kernel_R12_odd, odd["index_data"], corrupted_rows
)
corrupted_R12 = mp_symmetric(corrupted_R12)
carrier_corruption_relative = mp_frobenius(mp_difference(
    corrupted_R12, odd["forms"]["R12"]
)) / normalization
carrier_corruption_ok = carrier_corruption_relative > mp.mpf("1e-12")
check(
    "the frozen +1/10 carrier corruption is resolved at high precision",
    carrier_corruption_ok,
    (
        f"edge={first_edge}, effect={mp_text(carrier_corruption_relative, 8)}"
    ),
)


# Rank-one hostile Hessian control, evaluated without modifying any kernel.
column_zero = [row.get(0, mp.mpf(0)) for row in odd["carrier"]["rows"]]
column_norm = mp.sqrt(mp.fsum(value**2 for value in column_zero))
u = [value / column_norm for value in column_zero]
image = [
    mp.fsum(
        odd["carrier"]["rows"][row].get(column, mp.mpf(0)) * u[row]
        for row in range(ACTIVE)
    )
    for column in range(DATA)
]
synthetic_relative = mp.fsum(value**2 for value in image) / normalization
synthetic_ok = synthetic_relative > mp.mpf("1e-6")
check(
    "the rank-one hostile Hessian perturbation remains visible",
    synthetic_ok,
    f"effect={mp_text(synthetic_relative, 8)}",
)


# Convention diagnostic only: reverse the noncommutative product in the
# finest odd kernel convolution.
wrong_R12, _ = pull_orbit_kernel(
    kernel_R12_odd, odd["index_data"], odd["carrier"]["rows"], reverse=True
)
wrong_R12 = mp_symmetric(wrong_R12)
wrong_convention_relative = mp_frobenius(mp_difference(
    wrong_R12, odd["forms"]["R12"]
)) / normalization


# The primary matrices are loaded only now, after the independent forms and
# parity classification exist in memory.
primary_matrices = np.load(PRIMARY_MATRICES, allow_pickle=False)
primary_shape_ok = primary_matrices.shape == (3, DATA, DATA)
primary_errors = {}
adversarial_numpy = {}
for index, parity in enumerate(("even", "odd")):
    adversarial_numpy[parity] = mp_to_numpy(runtime[parity]["forms"]["R12"])
    primary_errors[parity] = (
        np.linalg.norm(primary_matrices[index] - adversarial_numpy[parity])
        / max(
            1.0,
            np.linalg.norm(primary_matrices[index]),
            np.linalg.norm(adversarial_numpy[parity]),
        )
    )
primary_matrix_ok = bool(
    primary_shape_ok and max(primary_errors.values()) < 1e-10
)
check(
    "both frozen primary matrices agree with the independent 180-digit reconstruction",
    primary_matrix_ok,
    str(primary_errors),
)


# Two direct scalar-action second derivatives, fixed before this run.
directions = {}
uniform = [mp.mpf(0) for _ in range(DATA)]
for vertex in range(VERTICES):
    uniform[vertex] = 1 / mp.sqrt(VERTICES)
directions["uniform_scale"] = uniform
local = [mp.mpf(0) for _ in range(DATA)]
local[0] = mp.mpf(1) / 2
local[1] = -mp.mpf(1) / 2
local[VERTICES + 2] = mp.mpf(1) / 2
local[VERTICES + 3] = -mp.mpf(1) / 2
directions["local_scale_strut"] = local

action_records = {}
action_controls_ok = True
primary_action_errors = []
for parity in ("even", "odd"):
    print(f"[{parity}] running two direct scalar-action second derivatives", flush=True)
    item = runtime[parity]
    zero_delta = [mp.mpf(0) for _ in range(ACTIVE)]
    base_action, base_branch = action_at_active_delta(
        item["model"], item["index_data"], item["geometry"],
        library["angle_record"], library["area_data"], state["mass"], zero_delta,
    )
    action_records[parity] = {}
    for name, vector in directions.items():
        edge_direction = carrier_direction(item["carrier"]["rows"], vector)
        derivatives = []
        branch_records = [base_branch]
        for step in (mp.mpf("1e-20"), mp.mpf("5e-21")):
            plus, plus_branch = action_at_active_delta(
                item["model"], item["index_data"], item["geometry"],
                library["angle_record"], library["area_data"], state["mass"],
                [step * value for value in edge_direction],
            )
            minus, minus_branch = action_at_active_delta(
                item["model"], item["index_data"], item["geometry"],
                library["angle_record"], library["area_data"], state["mass"],
                [-step * value for value in edge_direction],
            )
            derivatives.append((plus - 2 * base_action + minus) / step**2)
            branch_records.extend((plus_branch, minus_branch))
        direct = (4 * derivatives[1] - derivatives[0]) / 3
        predicted = quadratic_value(item["forms"]["R12"], vector)
        relative = abs(direct - predicted) / max(
            mp.mpf(1), abs(direct), abs(predicted)
        )
        branch_ok = bool(
            all(record["negative_counts"] == Counter({1: 2400})
                for record in branch_records)
            and min(record["minimum_leading_minor"] for record in branch_records) > 0
            and min(record["minimum_argument"] for record in branch_records)
            > mp.mpf("1e-6")
            and max(record["maximum_imaginary"] for record in branch_records)
            < mp.mpf("1e-140")
        )
        direction_ok = bool(relative < mp.mpf("1e-55") and branch_ok)
        action_controls_ok &= direction_ok

        primary_prediction = mp.mpf(str(
            float(np.asarray(vector, dtype=float)
                  @ primary_matrices[0 if parity == "even" else 1]
                  @ np.asarray(vector, dtype=float))
        ))
        primary_relative = abs(direct - primary_prediction) / max(
            mp.mpf(1), abs(direct), abs(primary_prediction)
        )
        primary_action_errors.append(primary_relative)
        action_records[parity][name] = {
            "direct": mp_text(direct),
            "predicted": mp_text(predicted),
            "relative_error": mp_text(relative),
            "primary_relative_error": mp_text(primary_relative),
            "branch_ok": branch_ok,
        }
check(
    "all four direct scalar-action controls reproduce the orbit-kernel forms",
    action_controls_ok,
    str({
        parity: {
            name: record["relative_error"]
            for name, record in rows.items()
        }
        for parity, rows in action_records.items()
    }),
)


all_controls = bool(
    provenance_ok
    and registry_ok
    and geometry_ok
    and background_ok
    and all(runtime[parity]["control_ok"] for parity in ("even", "odd"))
    and classification_finite
    and carrier_corruption_ok
    and synthetic_ok
    and primary_shape_ok
    and action_controls_ok
)

if not all_controls:
    outcome = "FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_CONTROL_FAILED"
elif high_precision_independent and primary_matrix_ok:
    outcome = (
        "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_"
        "ADVERSARIALLY_REPLICATED"
    )
elif (
    high_precision_dependent
    and not primary_matrix_ok
    and max(primary_action_errors) > mp.mpf("1e-10")
):
    outcome = "FINITE_HEIGHT_QUADRATIC_PRIMARY_REFUTED"
else:
    outcome = "FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_OPEN"

check(
    "the adversarial outcome follows the frozen hierarchy",
    outcome in {
        "FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_ADVERSARIALLY_REPLICATED",
        "FINITE_HEIGHT_QUADRATIC_PRIMARY_REFUTED",
        "FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_OPEN",
        "FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_CONTROL_FAILED",
    },
    outcome,
)


matrix_stack = np.stack([
    adversarial_numpy["even"],
    adversarial_numpy["odd"],
    adversarial_numpy["even"] - adversarial_numpy["odd"],
])
np.save(MATRIX_OUTPUT, matrix_stack, allow_pickle=False)
matrix_hash = sha256(MATRIX_OUTPUT)

artifact = {
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "method": (
        "180-digit identity-row 2T orbit kernels, exact group convolution, "
        "sparse direct carrier pullback, and scalar-action controls"
    ),
    "provenance": {
        "primary_result_commit": PRIMARY_RESULT_COMMIT,
        "adversarial_protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": hashes,
        "source_sha256": sha256(Path(__file__)),
    },
    "background": {
        "q": mp_text(state["q"]),
        "h": mp_text(state["h"]),
        "lambda": mp_text(state["lambda"]),
        "rho": mp_text(state["rho"]),
        "mass": mp_text(state["mass"]),
        "exact_elimination_residual": mp_text(state["elimination"]),
    },
    "derivative_steps": {
        name: mp_text(value) for name, value in DERIVATIVE_STEPS.items()
    },
    "parities": records,
    "classification": {
        "normalization_frobenius": mp_text(normalization),
        "relative_parity_differences": {
            level: mp_text(value) for level, value in d_values.items()
        },
        "e_step": mp_text(e_step),
        "e_total": mp_text(e_total),
        "difference_agreement": mp_text(difference_agreement),
        "independent_gate": mp_text(10 * e_total),
        "dependent_gate": mp_text(100 * e_total),
        "high_precision_independent": high_precision_independent,
        "high_precision_dependent": high_precision_dependent,
        "primary_matrix_relative_errors": {
            parity: f"{value:.17e}" for parity, value in primary_errors.items()
        },
        "carrier_corruption_relative_effect": mp_text(
            carrier_corruption_relative
        ),
        "synthetic_rank_one_relative_effect": mp_text(synthetic_relative),
        "wrong_group_product_diagnostic": mp_text(
            wrong_convention_relative
        ),
    },
    "direct_action_controls": action_records,
    "matrices": {
        "path": MATRIX_OUTPUT.name,
        "sha256": matrix_hash,
        "shape": list(matrix_stack.shape),
        "order": ["even_R12", "odd_R12", "even_minus_odd"],
    },
    "firewall": {
        "complete_binary64_hessian_assembler_loaded": False,
        "quadratic_eigenmodes_inspected": False,
        "continuum_target_parsed": False,
        "speed_target_parsed": False,
        "full_suite_run": False,
    },
    "interpretation": {
        "scope": "one finite-height background, one-sided quadratic tangent form",
        "nonlinear_integrability": "OPEN",
        "boundary_evolution": "NOT_TESTED",
        "gravitons": "NOT_TESTED",
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
