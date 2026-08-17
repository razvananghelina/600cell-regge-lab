#!/usr/bin/env python3
"""Internal Regge-curvature response of the 119 strong tangent modes.

Prior-art commit: 564b3a9.
Protocol commit: 2629f5e.
No continuum, speed, Planck or desired response target is loaded.
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

from flint import acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
TANGENT_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_internal_curvature_response.json"

PRIOR_ART_COMMIT = "564b3a9"
PROTOCOL_COMMIT = "2629f5e"
EXPECTED_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "tangent": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

DPS = 100
BALL_DPS = 80
mp.mp.dps = DPS
ctx.dps = BALL_DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
ARITHMETIC_FLOOR = mp.mpf("1e-70")
tests = passed = 0


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


def norm2(matrix):
    singular = la.svdvals(np.asarray(matrix))
    return float(singular[0]) if len(singular) else 0.0


hashes = {
    "tick": sha256(TICK_INPUT),
    "tangent": sha256(TANGENT_INPUT),
    "tangent_numeric": sha256(TANGENT_NUMERIC),
    "tangent_source": sha256(TANGENT_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
}
tick = json.loads(TICK_INPUT.read_text())
tangent_input = json.loads(TANGENT_INPUT.read_text())
numeric = np.load(TANGENT_NUMERIC)
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tangent_input["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and tangent_input["passed"] == tangent_input["tests"] == 19
    and tangent_input["numeric_archive_arrays"] == len(numeric.files) == 224
    and tangent_input["numeric_archive_sha256"] == hashes["tangent_numeric"]
)
check("all preregistered inputs have exact frozen provenance", provenance_ok, str(hashes))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_internal_curvature", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check(
    "the direct one-slab geometry import retains all 43 certificates",
    gro.tests == gro.passed == 43,
)


M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2


def load_named_functions(source, wanted):
    tree = ast.parse(source.read_text(), filename=str(source))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch in {source.name}: {wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"),
        globals(),
    )


load_named_functions(RANK_SOURCE, {
    "orbit_sort_key",
    "augment_boundary_orbits",
    "log_minus",
    "signed_volume_square",
    "angle_record",
    "area_data",
    "extended_edge_image",
    "group_and_index_data",
    "prepare_geometry",
})
load_named_functions(TANGENT_SOURCE, {
    "mp_frobenius",
    "mp_submatrix",
    "cluster_sorted",
    "high_precision_sector_bases",
    "high_precision_pattern_cache",
    "assemble_full_representative_kernels",
    "project_full_kernel",
    "mp_to_numpy",
    "mp_to_acb",
    "mp_matrix_to_acb",
    "acb_midpoint_and_radii",
    "expanded_types",
})

models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}


def triangle_area_square(index_data, record):
    x = [index_data["signed_base"][index] for index in record["indices"]]
    return (
        2 * (x[0] * x[1] + x[0] * x[2] + x[1] * x[2])
        - x[0] ** 2 - x[1] ** 2 - x[2] ** 2
    ) / 16


def extended_triangle_image(action, triangle):
    return tuple(sorted(gro.extended_image(action, vertex) for vertex in triangle))


def group_inverses(index_data):
    table = index_data["table"]
    identity = index_data["identity"]
    return tuple(
        next(
            right for right in range(24)
            if table[left, right] == identity and table[right, left] == identity
        )
        for left in range(24)
    )


def triangle_response_data(model, index_data, geometry, pattern_cache):
    """Build the direct sparse D_kappa and its representative orbit kernel."""
    records_by_triangle = {
        record["triangle"]: record for record in geometry["triangle_records"]
    }
    actions = index_data["actions"]
    identity = index_data["identity"]
    table = index_data["table"]
    inverses = group_inverses(index_data)

    unseen = set(geometry["triangles"])
    triangle_orbits = []
    triangle_location = {}
    while unseen:
        seed = min(unseen)
        ordered = tuple(extended_triangle_image(action, seed) for action in actions)
        if len(set(ordered)) != 24:
            raise RuntimeError("triangle orbit is not free")
        orbit_type = len(triangle_orbits)
        triangle_orbits.append(ordered)
        for group, triangle in enumerate(ordered):
            triangle_location[triangle] = (orbit_type, group)
        unseen -= set(ordered)

    causal = {}
    row_phase = {}
    minimum_area_square = mp.inf
    maximum_area_imaginary = mp.mpf(0)
    for record in geometry["triangle_records"]:
        q = triangle_area_square(index_data, record)
        maximum_area_imaginary = max(maximum_area_imaginary, abs(mp.im(q)))
        real_q = mp.re(q)
        minimum_area_square = min(minimum_area_square, abs(real_q))
        if real_q > 0:
            causal[record["triangle"]] = "spacelike"
            row_phase[record["triangle"]] = -I
        elif real_q < 0:
            causal[record["triangle"]] = "timelike"
            row_phase[record["triangle"]] = mp.mpf(1)
        else:
            causal[record["triangle"]] = "null"
            row_phase[record["triangle"]] = mp.nan

    base_deficits = {
        record["triangle"]: mp.pi if record["boundary"] else 2 * mp.pi
        for record in geometry["triangle_records"]
    }
    incidence = Counter()
    for simplex in geometry["simplex_records"]:
        angles = pattern_cache[simplex["pattern"]]["base_angles"]
        for hinge_index, triangle_index in enumerate(simplex["hinge_triangles"]):
            triangle = geometry["triangle_records"][triangle_index]["triangle"]
            base_deficits[triangle] += angles[hinge_index]
            incidence[triangle] += 1

    physical_counts = Counter()
    orbit_counts = Counter()
    maximum_base_imaginary = mp.mpf(0)
    base_values = defaultdict(list)
    for triangle, deficit in base_deficits.items():
        record = records_by_triangle[triangle]
        value = row_phase[triangle] * deficit
        maximum_base_imaginary = max(maximum_base_imaginary, abs(mp.im(value)))
        zero = abs(value) < mp.mpf("1e-40")
        key = (
            "boundary" if record["boundary"] else "internal",
            causal[triangle],
            "zero" if zero else "nonzero",
        )
        physical_counts[key] += 1
        base_values[key].append(mp.re(value))

    orbit_labels_constant = True
    for orbit in triangle_orbits:
        keys = []
        for triangle in orbit:
            record = records_by_triangle[triangle]
            value = row_phase[triangle] * base_deficits[triangle]
            keys.append((
                "boundary" if record["boundary"] else "internal",
                causal[triangle],
                "zero" if abs(value) < mp.mpf("1e-40") else "nonzero",
            ))
        orbit_labels_constant &= len(set(keys)) == 1
        orbit_counts[keys[0]] += 1

    internal_global_types = [
        orbit_type for orbit_type, orbit in enumerate(triangle_orbits)
        if not records_by_triangle[orbit[identity]]["boundary"]
    ]
    internal_position = {
        orbit_type: position
        for position, orbit_type in enumerate(internal_global_types)
    }

    direct = {name: defaultdict(lambda: mp.mpc(0)) for name in VARIANTS}
    for simplex in geometry["simplex_records"]:
        derivatives = pattern_cache[simplex["pattern"]]["derivatives"]
        for hinge_index, triangle_index in enumerate(simplex["hinge_triangles"]):
            record = geometry["triangle_records"][triangle_index]
            if record["boundary"]:
                continue
            triangle = record["triangle"]
            phase = row_phase[triangle]
            for column, edge_index in enumerate(simplex["indices"]):
                for name in VARIANTS:
                    direct[name][(triangle, edge_index)] += (
                        phase * derivatives[name][hinge_index, column]
                    )

    kernels = {name: defaultdict(lambda: mp.mpc(0)) for name in VARIANTS}
    for name in VARIANTS:
        for (triangle, edge_index), value in direct[name].items():
            global_type, row_group = triangle_location[triangle]
            if row_group != identity:
                continue
            edge_type, edge_group = divmod(edge_index, 24)
            kernels[name][(
                internal_position[global_type], edge_type, edge_group
            )] += value

    maximum_equivariance_residual = {}
    maximum_derivative_imaginary = {}
    for name in VARIANTS:
        predicted = {}
        for (internal_type, edge_type, relative_group), value in kernels[name].items():
            global_type = internal_global_types[internal_type]
            for row_group in range(24):
                triangle = triangle_orbits[global_type][row_group]
                edge_group = int(table[row_group, relative_group])
                predicted[(triangle, 24 * edge_type + edge_group)] = value
        keys = set(predicted) | set(direct[name])
        maximum_equivariance_residual[name] = max(
            (abs(direct[name][key] - predicted.get(key, 0)) for key in keys),
            default=mp.mpf(0),
        )
        maximum_derivative_imaginary[name] = max(
            (abs(mp.im(value)) for value in direct[name].values()),
            default=mp.mpf(0),
        )

        # Independently verify the preregistered relative-group convention.
        convention_residual = mp.mpf(0)
        for (triangle, edge_index), value in direct[name].items():
            global_type, row_group = triangle_location[triangle]
            edge_type, edge_group = divmod(edge_index, 24)
            relative = int(table[inverses[row_group], edge_group])
            expected = kernels[name].get((
                internal_position[global_type], edge_type, relative
            ), 0)
            convention_residual = max(convention_residual, abs(value - expected))
        maximum_equivariance_residual[name] = max(
            maximum_equivariance_residual[name], convention_residual
        )

    return {
        "triangle_orbits": tuple(triangle_orbits),
        "triangle_location": triangle_location,
        "internal_global_types": tuple(internal_global_types),
        "kernels": kernels,
        "physical_counts": physical_counts,
        "orbit_counts": orbit_counts,
        "incidence_counts": Counter(
            (
                "boundary" if records_by_triangle[triangle]["boundary"] else "internal",
                causal[triangle],
                count,
            )
            for triangle, count in incidence.items()
        ),
        "base_values": base_values,
        "minimum_absolute_area_square": minimum_area_square,
        "maximum_area_square_imaginary": maximum_area_imaginary,
        "maximum_base_kappa_imaginary": maximum_base_imaginary,
        "maximum_derivative_imaginary": maximum_derivative_imaginary,
        "maximum_equivariance_residual": maximum_equivariance_residual,
        "orbit_labels_constant": bool(orbit_labels_constant),
        "direct_nonzero_entries": {
            name: sum(abs(value) > ARITHMETIC_FLOOR for value in values.values())
            for name, values in direct.items()
        },
        "kernel_nonzero_entries": {
            name: sum(abs(value) > ARITHMETIC_FLOOR for value in values.values())
            for name, values in kernels.items()
        },
    }


def project_curvature_kernel(kernel, sector, internal_types=160):
    dimension = sector["dimension"]
    block = mp.matrix(internal_types * dimension, 95 * dimension)
    right = sector["right_representations"]
    for (row_type, column_type, group), value in kernel.items():
        representation = right[group]
        for row in range(dimension):
            for column in range(dimension):
                block[
                    row_type * dimension + row,
                    column_type * dimension + column,
                ] += value * representation[row, column]
    return block


def canonical_response_ball(block, dimension):
    old = expanded_types(0, 30, dimension)
    internal = expanded_types(30, 65, dimension)
    new = expanded_types(65, 95, dimension)
    nd = 30 * dimension
    xd = 35 * dimension
    size = xd + nd

    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    j_matrix = mp.matrix(size, size)
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

    rhs = mp.matrix(size, 2 * nd)
    k_xo = mp_submatrix(block, internal, old)
    k_oo = mp_submatrix(block, old, old)
    for row in range(xd):
        for column in range(nd):
            rhs[row, column] = -k_xo[row, column]
    for row in range(nd):
        for column in range(nd):
            rhs[xd + row, column] = k_oo[row, column]
        rhs[xd + row, nd + row] = 1

    j_ball = mp_matrix_to_acb(j_matrix)
    response_ball = j_ball.solve(mp_matrix_to_acb(rhs))
    midpoint, radii = acb_midpoint_and_radii(response_ball)
    return {
        "determinant": j_ball.det(),
        "midpoint": midpoint,
        "radii": radii,
    }


def orthonormal_columns(matrix):
    q, _ = la.qr(matrix, mode="economic")
    return q


def subspace_distance(left, right):
    q_left = orthonormal_columns(left)
    q_right = orthonormal_columns(right)
    overlap = la.svdvals(q_left.conj().T @ q_right)
    minimum = min(1.0, max(0.0, float(np.min(overlap))))
    return math.sqrt(max(0.0, 1 - minimum**2))


def extreme_subspace(matrix, count, branch):
    eigenvalues, eigenvectors = la.eig(matrix)
    moduli = np.abs(eigenvalues)
    if branch == "plus":
        order = np.argsort(moduli)[::-1]
        selected = order[:count]
        boundary = math.sqrt(moduli[order[count - 1]] * moduli[order[count]])
        selector = lambda value: abs(value) > boundary
        gap = float(moduli[order[count - 1]] / moduli[order[count]])
    else:
        order = np.argsort(moduli)
        selected = order[:count]
        boundary = math.sqrt(moduli[order[count - 1]] * moduli[order[count]])
        selector = lambda value: abs(value) < boundary
        gap = float(moduli[order[count]] / moduli[order[count - 1]])
    _, schur_vectors, selected_count = la.schur(
        matrix, output="complex", sort=selector
    )
    schur_basis = schur_vectors[:, :selected_count]
    direct_basis = orthonormal_columns(eigenvectors[:, selected])
    unselected = np.setdiff1d(np.arange(len(eigenvalues)), selected)
    separation = float(np.min(np.abs(
        eigenvalues[selected, None] - eigenvalues[None, unselected]
    )))
    return {
        "basis": schur_basis,
        "direct_basis": direct_basis,
        "selected_count": int(selected_count),
        "gap": gap,
        "boundary": boundary,
        "direct_distance": subspace_distance(schur_basis, direct_basis),
        "eigenvector_condition": float(np.linalg.cond(eigenvectors)),
        "spectral_separation": separation,
        "selected_moduli": np.sort(moduli[selected]),
    }


def singular_record(variant_data, branch=None):
    singular = {}
    direct_singular = {}
    ball_bounds = []
    tangent_bounds = []
    matrices = {}
    for name, data in variant_data.items():
        if branch is None:
            matrix = data["F"]
            direct_matrix = None
        else:
            extreme = data[f"extreme_{branch}"]
            matrix = data["F"] @ extreme["basis"]
            direct_matrix = data["F"] @ extreme["direct_basis"]
        matrices[name] = matrix
        singular[name] = la.svdvals(matrix)
        if direct_matrix is not None:
            direct_singular[name] = la.svdvals(direct_matrix)

        response_ball = data["D_norm"] * data["Z_radius"]
        ball_bounds.append(response_ball)
        if branch is not None:
            extreme = data[f"extreme_{branch}"]
            angle_bound = (
                data["tangent_radius"] * extreme["eigenvector_condition"]
                / max(1e-300, extreme["spectral_separation"])
            )
            tangent_bounds.append(norm2(data["F"]) * angle_bound)

    op = singular["operational_primary"]
    epsilon_step = (
        float(np.max(np.abs(op - singular["operational_shadow"])))
        + float(np.max(np.abs(
            singular["validation_primary"] - singular["validation_shadow"]
        )))
        + float(np.max(np.abs(op - singular["validation_primary"])))
    )
    epsilon_direct = max(
        (
            float(np.max(np.abs(singular[name] - direct_singular[name])))
            for name in VARIANTS
        ),
        default=0.0,
    )
    epsilon_ball = max(ball_bounds, default=0.0)
    epsilon_tangent = max(tangent_bounds, default=0.0)
    epsilon_binary = (
        10 * np.finfo(float).eps * max(matrices["operational_primary"].shape)
        * max(1.0, float(op[0]))
    )
    epsilon = (
        epsilon_step + epsilon_direct + epsilon_ball + epsilon_tangent
        + epsilon_binary + 1e-70
    )
    zero = op <= 10 * epsilon
    nonzero = op > 100 * epsilon
    open_flags = ~(zero | nonzero)
    if np.all(zero):
        label = "ZERO"
    elif np.all(nonzero):
        label = "INJECTIVE"
    else:
        label = "PARTIAL_OR_OPEN"
    return {
        "label": label,
        "resolved_rank": int(np.sum(nonzero)),
        "zero_count": int(np.sum(zero)),
        "open_count": int(np.sum(open_flags)),
        "columns": int(matrices["operational_primary"].shape[1]),
        "minimum_singular": float(op[-1]),
        "maximum_singular": float(op[0]),
        "condition": float(op[0] / op[-1]) if op[-1] else math.inf,
        "epsilon_singular": epsilon,
        "epsilon_step": epsilon_step,
        "epsilon_direct": epsilon_direct,
        "epsilon_ball": epsilon_ball,
        "epsilon_tangent": epsilon_tangent,
        "epsilon_binary": epsilon_binary,
        "variant_minimum_singular": {
            name: float(values[-1]) for name, values in singular.items()
        },
        "variant_maximum_singular": {
            name: float(values[0]) for name, values in singular.items()
        },
    }


def sf(value):
    return f"{float(value):.17e}"


def serialize_counter(counter):
    return {
        "|".join(str(part) for part in key): int(value)
        for key, value in sorted(counter.items())
    }


def serialize_singular(record):
    return {
        "label": record["label"],
        "resolved_rank": record["resolved_rank"],
        "zero_count": record["zero_count"],
        "open_count": record["open_count"],
        "columns": record["columns"],
        "minimum_singular": sf(record["minimum_singular"]),
        "maximum_singular": sf(record["maximum_singular"]),
        "condition": sf(record["condition"]),
        "epsilon_singular": sf(record["epsilon_singular"]),
        "epsilon_step": sf(record["epsilon_step"]),
        "epsilon_direct": sf(record["epsilon_direct"]),
        "epsilon_ball": sf(record["epsilon_ball"]),
        "epsilon_tangent": sf(record["epsilon_tangent"]),
        "epsilon_binary": sf(record["epsilon_binary"]),
        "variant_minimum_singular": {
            name: sf(value)
            for name, value in record["variant_minimum_singular"].items()
        },
        "variant_maximum_singular": {
            name: sf(value)
            for name, value in record["variant_maximum_singular"].items()
        },
    }


print("=" * 78)
print("INTERNAL REGGE-CURVATURE RESPONSE OF THE 119 STRONG TANGENT MODES")
print("=" * 78)

records = {}
all_full_records = []
all_strong_records = []
global_controls = provenance_ok and gro.tests == gro.passed == 43

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing triangle carrier and direct curvature Jacobian", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)
    sectors, sector_control = high_precision_sector_bases(index_data)

    s = mp.mpf(state[0])
    kind_values = {
        "old": L0_SQUARE,
        "internal": mp.exp(s) * L0_SQUARE - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s) * L0_SQUARE,
    }
    pattern_cache, branch_control = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    curvature_data = triangle_response_data(
        model, index_data, geometry, pattern_cache
    )

    expected_physical = Counter({
        ("boundary", "spacelike", "nonzero"): 2400,
        ("internal", "spacelike", "zero"): 2400,
        ("internal", "timelike", "nonzero"): 1440,
    })
    expected_orbits = Counter({
        ("boundary", "spacelike", "nonzero"): 100,
        ("internal", "spacelike", "zero"): 100,
        ("internal", "timelike", "nonzero"): 60,
    })
    carrier_ok = bool(
        len(geometry["triangles"]) == 6240
        and len(curvature_data["triangle_orbits"]) == 260
        and all(len(orbit) == 24 for orbit in curvature_data["triangle_orbits"])
        and len(curvature_data["internal_global_types"]) == 160
        and curvature_data["physical_counts"] == expected_physical
        and curvature_data["orbit_counts"] == expected_orbits
        and curvature_data["orbit_labels_constant"]
        and curvature_data["minimum_absolute_area_square"] > mp.mpf("1e-4")
        and curvature_data["maximum_area_square_imaginary"] < ARITHMETIC_FLOOR
        and curvature_data["maximum_base_kappa_imaginary"] < ARITHMETIC_FLOOR
    )
    check(
        f"{parity}: all triangle, causal and free-orbit census gates pass",
        carrier_ok,
        f"physical={dict(curvature_data['physical_counts'])}",
    )

    local_ok = bool(
        branch_control["entry_pass"]
        and branch_control["base_negative_counts"] == Counter({1: 2400})
        and branch_control["displaced_negative_counts"] == Counter({1: 1600})
        and all(
            value < ARITHMETIC_FLOOR
            for value in curvature_data["maximum_derivative_imaginary"].values()
        )
    )
    check(
        f"{parity}: local angle branches and physicalized derivative reality pass",
        local_ok,
        "max imag=" + mp.nstr(
            max(curvature_data["maximum_derivative_imaginary"].values()), 8
        ),
    )

    equivariance_ok = all(
        value < ARITHMETIC_FLOOR
        for value in curvature_data["maximum_equivariance_residual"].values()
    )
    check(
        f"{parity}: direct 3840-row incidence equals the orbit kernel exactly",
        equivariance_ok,
        "max residual=" + mp.nstr(
            max(curvature_data["maximum_equivariance_residual"].values()), 8
        ),
    )

    hessian_kernels, hessian_control = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    basis_ok = bool(
        sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and all(
            value < ARITHMETIC_FLOOR
            for key, value in sector_control.items() if key.startswith("maximum_")
        )
        and hessian_control["maximum_imaginary"] < ARITHMETIC_FLOOR
    )

    determinant_ok = True
    selection_ok = True
    stored_ok = True
    sector_records = []
    strong_weight = 0
    homogeneous_open = None

    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        trivial = bool(sector["constant_overlap"] > mp.mpf("0.5"))
        count = 4 if trivial else 5 * dimension
        strong_weight += dimension * count
        print(
            f"[{parity}] sector {sector_index + 1}/7 d={dimension}, k={count}: "
            "Flint solve and curvature ranks",
            flush=True,
        )

        stored = tangent_input["parities"][parity]["sectors"][sector_index]
        stored_ok &= bool(
            stored["sector_index"] == sector_index
            and stored["dimension"] == dimension
            and abs(float(stored["constant_overlap"]) - float(sector["constant_overlap"]))
            < 1e-12
            and stored["resolved_off_unit_count_minimal"] == 2 * count
        )

        hessian_blocks = {
            name: project_full_kernel(kernel, sector)
            for name, kernel in hessian_kernels.items()
        }
        curvature_blocks = {
            name: project_curvature_kernel(kernel, sector)
            for name, kernel in curvature_data["kernels"].items()
        }
        variant_data = {}
        for name in VARIANTS:
            response = canonical_response_ball(hessian_blocks[name], dimension)
            determinant_ok &= not response["determinant"].contains(0)
            nd = 30 * dimension
            old_projection = np.zeros((nd, 2 * nd), dtype=np.complex128)
            old_projection[:, :nd] = np.eye(nd)
            z_matrix = np.vstack((old_projection, response["midpoint"]))
            z_radii = np.vstack((
                np.zeros_like(old_projection, dtype=float), response["radii"]
            ))
            d_matrix = mp_to_numpy(curvature_blocks[name])
            f_matrix = d_matrix @ z_matrix
            tangent = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_midpoint"
            ]
            tangent_radii = numeric[
                f"{parity}_sector{sector_index}_{name}_tangent_radii"
            ]
            plus = extreme_subspace(tangent, count, "plus")
            minus = extreme_subspace(tangent, count, "minus")
            selection_ok &= bool(
                plus["selected_count"] == count
                and minus["selected_count"] == count
                and plus["gap"] > 2
                and minus["gap"] > 2
                and plus["direct_distance"] < 1e-4
                and minus["direct_distance"] < 1e-4
            )
            variant_data[name] = {
                "F": f_matrix,
                "D_norm": norm2(d_matrix),
                "Z_radius": float(la.norm(z_radii, "fro")),
                "tangent_radius": float(la.norm(tangent_radii, "fro")),
                "extreme_plus": plus,
                "extreme_minus": minus,
            }

        full_record = singular_record(variant_data)
        plus_record = singular_record(variant_data, "plus")
        minus_record = singular_record(variant_data, "minus")
        all_full_records.append({
            "parity": parity, "sector_index": sector_index, **full_record
        })
        all_strong_records.extend((
            {"parity": parity, "sector_index": sector_index,
             "branch": "plus", **plus_record},
            {"parity": parity, "sector_index": sector_index,
             "branch": "minus", **minus_record},
        ))

        if trivial:
            op_tangent = numeric[
                f"{parity}_sector{sector_index}_operational_primary_tangent_midpoint"
            ]
            moduli = np.sort(np.abs(la.eigvals(op_tangent)))
            homogeneous_open = {
                "contracting_fifth_modulus": float(moduli[4]),
                "expanding_fifth_modulus": float(moduli[-5]),
                "contracting_strong_gap": variant_data[
                    "operational_primary"
                ]["extreme_minus"]["gap"],
                "expanding_strong_gap": variant_data[
                    "operational_primary"
                ]["extreme_plus"]["gap"],
            }

        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "trivial": trivial,
            "selected_count": count,
            "multiplicity_weighted_selected_count": dimension * count,
            "minimum_plus_gap": min(
                data["extreme_plus"]["gap"] for data in variant_data.values()
            ),
            "minimum_minus_gap": min(
                data["extreme_minus"]["gap"] for data in variant_data.values()
            ),
            "maximum_plus_direct_distance": max(
                data["extreme_plus"]["direct_distance"]
                for data in variant_data.values()
            ),
            "maximum_minus_direct_distance": max(
                data["extreme_minus"]["direct_distance"]
                for data in variant_data.values()
            ),
            "full": full_record,
            "plus": plus_record,
            "minus": minus_record,
        })

    reconstruction_ok = bool(
        basis_ok and determinant_ok and selection_ok and stored_ok
        and strong_weight == 119 and homogeneous_open is not None
    )
    check(
        f"{parity}: all sector bases, Flint solves and frozen 119-mode selections pass",
        reconstruction_ok,
        f"weight={strong_weight}, determinant={determinant_ok}, selection={selection_ok}",
    )
    classification_ok = bool(
        len(sector_records) == 7
        and all(
            item[key]["resolved_rank"] + item[key]["zero_count"]
            + item[key]["open_count"] == item[key]["columns"]
            for item in sector_records for key in ("full", "plus", "minus")
        )
    )
    check(
        f"{parity}: all 7 full and 14 restricted spectra receive complete rank ledgers",
        classification_ok,
    )

    controls_ok = bool(
        carrier_ok and local_ok and equivariance_ok and reconstruction_ok
        and classification_ok
    )
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "curvature_data": curvature_data,
        "strong_weight": strong_weight,
        "homogeneous_open_pair": homogeneous_open,
        "sectors": sector_records,
    }


full_labels = Counter(item["label"] for item in all_full_records)
strong_labels = Counter(item["label"] for item in all_strong_records)
all_full_injective = len(all_full_records) == 14 and full_labels == Counter({"INJECTIVE": 14})
all_strong_injective = (
    len(all_strong_records) == 28 and strong_labels == Counter({"INJECTIVE": 28})
)
all_strong_zero = len(all_strong_records) == 28 and strong_labels == Counter({"ZERO": 28})

if not global_controls:
    outcome = "INTERNAL_CURVATURE_RESPONSE_CONTROL_FAILED"
elif all_full_injective and all_strong_injective:
    outcome = "FULL_BOUNDARY_PHASE_CURVATURE_INJECTIVE"
elif all_strong_injective:
    outcome = "STRONG_TANGENT_CURVATURE_INJECTIVE"
elif all_strong_zero:
    outcome = "STRONG_TANGENT_CURVATURE_ZERO"
else:
    outcome = "STRONG_TANGENT_CURVATURE_PARTIAL_OR_OPEN"

check(
    "the preregistered census contains exactly 14 full maps and 28 strong restrictions",
    len(all_full_records) == 14 and len(all_strong_records) == 28,
    f"full={dict(full_labels)}, strong={dict(strong_labels)}",
)
check(
    "the frozen hierarchy assigns the internal-curvature outcome",
    outcome in {
        "INTERNAL_CURVATURE_RESPONSE_CONTROL_FAILED",
        "FULL_BOUNDARY_PHASE_CURVATURE_INJECTIVE",
        "STRONG_TANGENT_CURVATURE_INJECTIVE",
        "STRONG_TANGENT_CURVATURE_ZERO",
        "STRONG_TANGENT_CURVATURE_PARTIAL_OR_OPEN",
    },
    f"outcome={outcome}",
)


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "planck_target_parsed": False,
    "outcome": outcome,
    "attempt_counts": {
        "full_phase_maps": len(all_full_records),
        "strong_restrictions": len(all_strong_records),
    },
    "full_label_counts": dict(full_labels),
    "strong_label_counts": dict(strong_labels),
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "strong_count_with_multiplicity": record["strong_weight"],
            "triangle_census": {
                "triangles": len(record["curvature_data"]["triangle_location"]),
                "free_orbits": len(record["curvature_data"]["triangle_orbits"]),
                "internal_orbits": len(record["curvature_data"]["internal_global_types"]),
                "physical_counts": serialize_counter(
                    record["curvature_data"]["physical_counts"]
                ),
                "orbit_counts": serialize_counter(
                    record["curvature_data"]["orbit_counts"]
                ),
                "incidence_counts": serialize_counter(
                    record["curvature_data"]["incidence_counts"]
                ),
                "minimum_absolute_area_square": mp.nstr(
                    record["curvature_data"]["minimum_absolute_area_square"], 70
                ),
                "maximum_base_kappa_imaginary": mp.nstr(
                    record["curvature_data"]["maximum_base_kappa_imaginary"], 8
                ),
                "maximum_derivative_imaginary": {
                    name: mp.nstr(value, 8)
                    for name, value in record["curvature_data"][
                        "maximum_derivative_imaginary"
                    ].items()
                },
                "maximum_equivariance_residual": {
                    name: mp.nstr(value, 8)
                    for name, value in record["curvature_data"][
                        "maximum_equivariance_residual"
                    ].items()
                },
                "direct_nonzero_entries": record["curvature_data"][
                    "direct_nonzero_entries"
                ],
                "representative_nonzero_entries": record["curvature_data"][
                    "kernel_nonzero_entries"
                ],
            },
            "homogeneous_open_pair": {
                key: sf(value)
                for key, value in record["homogeneous_open_pair"].items()
            },
            "sectors": [
                {
                    "sector_index": sector["sector_index"],
                    "irrep_dimension": sector["dimension"],
                    "trivial": sector["trivial"],
                    "selected_count": sector["selected_count"],
                    "multiplicity_weighted_selected_count": sector[
                        "multiplicity_weighted_selected_count"
                    ],
                    "minimum_expanding_gap": sf(sector["minimum_plus_gap"]),
                    "minimum_contracting_gap": sf(sector["minimum_minus_gap"]),
                    "maximum_expanding_direct_schur_distance": sf(
                        sector["maximum_plus_direct_distance"]
                    ),
                    "maximum_contracting_direct_schur_distance": sf(
                        sector["maximum_minus_direct_distance"]
                    ),
                    "full_phase_response": serialize_singular(sector["full"]),
                    "expanding_response": serialize_singular(sector["plus"]),
                    "contracting_response": serialize_singular(sector["minus"]),
                }
                for sector in record["sectors"]
            ],
        }
        for parity, record in records.items()
    },
    "classification": {
        "zero_nonzero_and_rank": "DERIVED COMPUTATIONAL",
        "response_norm_magnitudes": "PATTERN ONLY",
        "exact_curvature_preserving_lapse_interpretation": (
            "REFUTED"
            if outcome in {
                "FULL_BOUNDARY_PHASE_CURVATURE_INJECTIVE",
                "STRONG_TANGENT_CURVATURE_INJECTIVE",
            }
            else "OPEN"
        ),
        "physical_graviton_interpretation": "OPEN",
        "constraint_satisfaction": "OPEN",
        "continuum_dispersion": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Full maps: {dict(full_labels)} / 14")
print(f"Strong restrictions: {dict(strong_labels)} / 28")
if all_full_records:
    print(
        "Full minimum singular range: "
        f"{min(item['minimum_singular'] for item in all_full_records):.6g} ... "
        f"{max(item['minimum_singular'] for item in all_full_records):.6g}"
    )
if all_strong_records:
    print(
        "Strong minimum singular range: "
        f"{min(item['minimum_singular'] for item in all_strong_records):.6g} ... "
        f"{max(item['minimum_singular'] for item in all_strong_records):.6g}"
    )
print(f"Results: {passed}/{tests} tests passed.")
sys.exit(0 if passed == tests else 1)
