#!/usr/bin/env python3
"""Full 2T-resolved canonical-rank census of the dynamic 600-cell dust slab.

Prior-art commit: f266cc2.
Protocol commits: 6b61e70, ce13a0e.
No continuum spectrum, speed target, or full tangent spectrum is evaluated.
"""

from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
import multiprocessing as mp_pool
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_dust_full_anisotropic_legendre_rank.json"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
TANGENT_INPUT = HERE / "gravity_600cell_dust_dynamic_tangent.json"
GLUING_INPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
PRIOR_ART_COMMIT = "f266cc2"
PROTOCOL_COMMIT = "6b61e70"
CLARIFICATION_COMMIT = "ce13a0e"
ROUNDING_CORRECTION_COMMIT = "5bd2cd1"
EXPECTED_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "tangent": "1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5",
}
DPS = 100
mp.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
DIRECT_STEPS = (mp.mpf("1e-8"), mp.mpf("5e-9"))
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
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
    values = la.svdvals(np.asarray(matrix))
    return float(values[0]) if len(values) else 0.0


def sparse_norm2(matrix):
    matrix = sp.csr_matrix(matrix)
    if matrix.nnz == 0:
        return 0.0
    try:
        return float(spla.svds(
            matrix, k=1, which="LM", return_singular_vectors=False,
            tol=1e-10, maxiter=20000,
        )[0])
    except Exception:
        return float(sp.linalg.norm(matrix))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_full_anisotropic", HERE / "verify_gravity_global_regge_orbits.py"
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
    "the imported complete one-slab geometry retains all 43 certificates",
    gro.tests == gro.passed == 43,
)


def orbit_sort_key(orbit, phase):
    representative = min(orbit)
    logical = tuple(vertex % 120 for vertex in representative)
    phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
    return phase_pair, tuple(sorted(orbit))


def augment_boundary_orbits(base):
    old_orbits = tuple(sorted(
        gro.orbit_partition(base["old_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    final_orbits = tuple(sorted(
        gro.orbit_partition(base["new_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    return {
        **base,
        "old_orbits": old_orbits,
        "final_orbits": final_orbits,
    }


models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}

tick = json.loads(TICK_INPUT.read_text())
tangent_input = json.loads(TANGENT_INPUT.read_text())
gluing = json.loads(GLUING_INPUT.read_text())
hashes = {"tick": sha256(TICK_INPUT), "tangent": sha256(TANGENT_INPUT)}
input_ok = bool(
    hashes == EXPECTED_HASHES
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tangent_input["outcome"] in {
        "DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT",
        "DYNAMIC_SHAPE_TANGENT_SCHEDULE_ROBUST",
        "DYNAMIC_SHAPE_TANGENT_SCHEDULE_OPEN",
    }
    and set(tick["solutions"]) == set(tangent_input["parities"]) == {"even", "odd"}
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
)
check("the two frozen dynamic inputs have exact committed provenance", input_ok, str(hashes))


M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2


def log_minus(value):
    scale = max(mp.mpf(1), abs(value))
    if abs(mp.im(value)) < mp.mpf("1e-80") * scale:
        real = mp.re(value)
        if real < 0:
            return mp.log(-real) - I * mp.pi
        return mp.log(real)
    return mp.log(value)


def signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices) - 1
    if dimension == 0:
        return mp.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = mp.matrix([
        [
            (squared[base][left] + squared[base][right] - squared[left][right]) / 2
            for right in others
        ]
        for left in others
    ])
    return mp.det(gram) / (mp.factorial(dimension) ** 2)


def angle_record(edge_values):
    squared = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
    for value, (left, right) in zip(edge_values, LOCAL_EDGES):
        squared[left][right] = squared[right][left] = value
    gram = mp.matrix([
        [
            (squared[0][left] + squared[0][right] - squared[left][right]) / 2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram**-1
    simplex_volume_square = signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    leading = []
    for size in range(1, 5):
        principal = mp.matrix([
            [gram[left, right] for right in range(size)]
            for left in range(size)
        ])
        leading.append(mp.det(principal))
    signs = [1] + [1 if value > 0 else -1 if value < 0 else 0 for value in leading]
    negative_directions = None if 0 in signs else sum(
        left != right for left, right in zip(signs, signs[1:])
    )

    angles = [mp.mpc(0) for _ in LOCAL_HINGES]
    minimum_argument = mp.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge = tuple(
            vertex for vertex in range(5) if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = signed_volume_square(squared, hinge)
        gram_derivative = mp.matrix(4, 4)
        opposite_edge = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                gram_derivative[left - 1, right - 1] = (
                    int({0, left} == opposite_edge)
                    + int({0, right} == opposite_edge)
                    - int(left != right and {left, right} == opposite_edge)
                ) / 2
        product = inverse * gram_derivative
        volume_derivative = simplex_volume_square * sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            mp.sqrt(mp.mpc(facet_volume_squares[omitted_a]))
            * mp.sqrt(mp.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16 * volume_derivative / denominator
        sine = -(mp.mpf(4) / 3) * (
            mp.sqrt(mp.mpc(hinge_volume_square))
            * mp.sqrt(mp.mpc(simplex_volume_square))
        ) / denominator
        argument = cosine + I * sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[LOCAL_HINGE_INDEX[hinge]] = -I * log_minus(argument)
    return angles, {
        "negative_directions": negative_directions,
        "minimum_leading_minor": min(abs(value) for value in leading),
        "minimum_argument": minimum_argument,
    }


def area_data(edge_values):
    x = list(edge_values)
    q = (
        2 * (x[0] * x[1] + x[0] * x[2] + x[1] * x[2])
        - x[0] ** 2 - x[1] ** 2 - x[2] ** 2
    ) / 16
    area = mp.sqrt(mp.mpc(q))
    q_plain = (
        (x[1] + x[2] - x[0]) / 8,
        (x[0] + x[2] - x[1]) / 8,
        (x[0] + x[1] - x[2]) / 8,
    )
    q_log = [q_plain[index] * x[index] for index in range(3)]
    q_hessian = [[mp.mpf(0) for _ in range(3)] for _ in range(3)]
    for row in range(3):
        for column in range(3):
            plain_second = mp.mpf(-1) / 8 if row == column else mp.mpf(1) / 8
            q_hessian[row][column] = plain_second * x[row] * x[column]
            if row == column:
                q_hessian[row][column] += q_plain[row] * x[row]
    gradient = [value / (2 * area) for value in q_log]
    hessian = [[
        q_hessian[row][column] / (2 * area)
        - q_log[row] * q_log[column] / (4 * area**3)
        for column in range(3)
    ] for row in range(3)]
    return area, gradient, hessian


def extended_edge_image(action, edge):
    return tuple(sorted(gro.extended_image(action, vertex) for vertex in edge))


def group_and_index_data(model, state):
    actions = tuple(sorted(tuple(int(value) for value in action) for action in model["stabilizer"]))
    action_arrays = tuple(np.asarray(action, dtype=np.int16) for action in actions)
    action_index = {action: index for index, action in enumerate(actions)}
    identity = action_index[tuple(range(120))]
    table = np.empty((24, 24), dtype=np.int16)
    for left, action in enumerate(action_arrays):
        for right, other in enumerate(action_arrays):
            table[left, right] = action_index[tuple(action[other])]
    inverses = tuple(
        next(
            right for right in range(24)
            if table[left, right] == identity and table[right, left] == identity
        )
        for left in range(24)
    )
    orders = []
    for element in range(24):
        product = identity
        order = 0
        while True:
            order += 1
            product = int(table[product, element])
            if product == identity:
                break
        orders.append(order)
    unseen = set(range(24))
    classes = []
    while unseen:
        seed = min(unseen)
        conjugacy_class = tuple(sorted({
            int(table[table[group, seed], inverses[group]]) for group in range(24)
        }))
        classes.append(conjugacy_class)
        unseen -= set(conjugacy_class)
    classes = tuple(sorted(classes, key=lambda item: (len(item), orders[item[0]], item)))

    global_orbits = model["old_orbits"] + model["edge_orbits"] + model["final_orbits"]
    orbit_edges = []
    edge_to_index = {}
    for orbit_type, orbit in enumerate(global_orbits):
        seed = min(orbit)
        ordered = tuple(extended_edge_image(action, seed) for action in action_arrays)
        if len(set(ordered)) != 24 or set(ordered) != set(orbit):
            raise RuntimeError("edge orbit is not free regular action")
        orbit_edges.append(ordered)
        for group, edge in enumerate(ordered):
            edge_to_index[edge] = 24 * orbit_type + group

    s = mp.mpf(state[0])
    r = mp.mpf(state[1])
    rho = RHO0 * mp.exp(r)
    values = {
        "old": L0_SQUARE,
        "internal": mp.exp(s) * L0_SQUARE - rho,
        "pole": -rho,
        "new": mp.exp(2 * s) * L0_SQUARE,
    }
    edge_kind = [None] * 2280
    signed_base = [None] * 2280
    for edge, index in edge_to_index.items():
        if edge in model["old_edges"]:
            kind = "old"
        elif edge in model["new_edges"]:
            kind = "new"
        elif edge[1] - edge[0] == 120:
            kind = "pole"
        else:
            kind = "internal"
        edge_kind[index] = kind
        signed_base[index] = values[kind]

    return {
        "actions": action_arrays,
        "table": table,
        "identity": identity,
        "orders": tuple(orders),
        "classes": classes,
        "global_orbits": global_orbits,
        "orbit_edges": tuple(orbit_edges),
        "edge_to_index": edge_to_index,
        "edge_kind": tuple(edge_kind),
        "signed_base": tuple(signed_base),
        "rho": rho,
    }


def prepare_geometry(model, index_data):
    edge_to_index = index_data["edge_to_index"]
    triangles = tuple(sorted(model["faces"][3]))
    triangle_index = {triangle: index for index, triangle in enumerate(triangles)}
    triangle_records = []
    for triangle in triangles:
        edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
        indices = tuple(edge_to_index[edge] for edge in edges)
        triangle_records.append({
            "triangle": triangle,
            "indices": indices,
            "boundary": triangle in model["boundary_triangles"],
        })
    simplex_records = []
    patterns = Counter()
    for simplex in sorted(model["slab"]):
        indices = tuple(
            edge_to_index[tuple(sorted((simplex[left], simplex[right])))]
            for left, right in LOCAL_EDGES
        )
        pattern = tuple(index_data["edge_kind"][index] for index in indices)
        hinge_triangles = tuple(
            triangle_index[tuple(sorted(simplex[position] for position in hinge))]
            for hinge in LOCAL_HINGES
        )
        simplex_records.append({
            "simplex": simplex,
            "indices": indices,
            "pattern": pattern,
            "hinge_triangles": hinge_triangles,
        })
        patterns[pattern] += 1
    return {
        "triangles": triangles,
        "triangle_records": tuple(triangle_records),
        "simplex_records": tuple(simplex_records),
        "patterns": patterns,
        "pole_indices": tuple(
            index for index, kind in enumerate(index_data["edge_kind"]) if kind == "pole"
        ),
    }


def build_pattern_cache(patterns, kind_values):
    cache = {}
    maximum_cross = mp.mpf(0)
    maximum_op_proxy = mp.mpf(0)
    maximum_val_proxy = mp.mpf(0)
    entry_pass = True
    minimum_minor = mp.inf
    minimum_argument = mp.inf
    negative_counts = Counter()
    for pattern in sorted(patterns):
        values = [kind_values[kind] for kind in pattern]
        base_angles, branch = angle_record(values)
        negative_counts[branch["negative_directions"]] += patterns[pattern]
        minimum_minor = min(minimum_minor, branch["minimum_leading_minor"])
        minimum_argument = min(minimum_argument, branch["minimum_argument"])
        derivatives_mp = {}
        for name, step in DERIVATIVE_STEPS.items():
            matrix = [[mp.mpc(0) for _ in LOCAL_EDGES] for _ in LOCAL_HINGES]
            for column in range(10):
                plus = list(values)
                minus = list(values)
                plus[column] *= mp.exp(step)
                minus[column] *= mp.exp(-step)
                plus_angles, plus_branch = angle_record(plus)
                minus_angles, minus_branch = angle_record(minus)
                for displaced_branch in (plus_branch, minus_branch):
                    negative_counts[(name, displaced_branch["negative_directions"])] += 1
                    minimum_minor = min(
                        minimum_minor, displaced_branch["minimum_leading_minor"]
                    )
                    minimum_argument = min(
                        minimum_argument, displaced_branch["minimum_argument"]
                    )
                for row in range(10):
                    matrix[row][column] = (
                        plus_angles[row] - minus_angles[row]
                    ) / (2 * step)
            derivatives_mp[name] = matrix

        for row in range(10):
            for column in range(10):
                op = derivatives_mp["operational_primary"][row][column]
                op_shadow = derivatives_mp["operational_shadow"][row][column]
                val = derivatives_mp["validation_primary"][row][column]
                val_shadow = derivatives_mp["validation_shadow"][row][column]
                cross = abs(op - val)
                op_proxy = abs(op - op_shadow)
                val_proxy = abs(val - val_shadow)
                maximum_cross = max(maximum_cross, cross)
                maximum_op_proxy = max(maximum_op_proxy, op_proxy)
                maximum_val_proxy = max(maximum_val_proxy, val_proxy)
                entry_pass &= bool(
                    cross <= 10 * (op_proxy + val_proxy + mp.mpf("1e-70"))
                )
        cache[pattern] = {
            "base_angles_mp": tuple(base_angles),
            "derivatives": {
                name: np.array([
                    [complex(value) for value in row] for row in matrix
                ], dtype=np.complex128)
                for name, matrix in derivatives_mp.items()
            },
        }
    return cache, {
        "entry_pass": bool(entry_pass),
        "maximum_cross": maximum_cross,
        "maximum_operational_proxy": maximum_op_proxy,
        "maximum_validation_proxy": maximum_val_proxy,
        "minimum_leading_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "base_negative_counts": Counter({
            key: value for key, value in negative_counts.items() if isinstance(key, int)
        }),
        "all_displaced_have_one_negative": all(
            not isinstance(key, tuple) or key[1] == 1 for key in negative_counts
        ),
    }


def assemble_full_hessians(model, index_data, geometry, pattern_cache):
    signed_base = index_data["signed_base"]
    curvature = [mp.pi if record["boundary"] else 2 * mp.pi
                 for record in geometry["triangle_records"]]
    for simplex in geometry["simplex_records"]:
        angles = pattern_cache[simplex["pattern"]]["base_angles_mp"]
        for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
            curvature[triangle] += angles[hinge_index]

    triangle_local = []
    gradient = np.zeros(2280, dtype=np.complex128)
    common = np.zeros((2280, 2280), dtype=np.complex128)
    common_absolute_terms = np.zeros((2280, 2280), dtype=np.float64)
    common_term_counts = np.zeros((2280, 2280), dtype=np.uint16)
    action = mp.mpc(0)
    for triangle_index, record in enumerate(geometry["triangle_records"]):
        indices = np.asarray(record["indices"], dtype=int)
        values = [signed_base[index] for index in indices]
        area, area_gradient, area_hessian = area_data(values)
        ag = np.asarray([complex(value) for value in area_gradient])
        ah = np.asarray([[complex(value) for value in row] for row in area_hessian])
        epsilon = curvature[triangle_index]
        coefficient = complex(-I * epsilon)
        area_term = coefficient * ah
        gradient[indices] += coefficient * ag
        common[np.ix_(indices, indices)] += area_term
        common_absolute_terms[np.ix_(indices, indices)] += np.abs(area_term)
        common_term_counts[np.ix_(indices, indices)] += 1
        action += -I * area * epsilon
        triangle_local.append({"indices": indices, "area_gradient": ag})

    dust_gradient = -(4 * mp.pi * MASS / 120) * mp.sqrt(index_data["rho"])
    dust_hessian = -(2 * mp.pi * MASS / 120) * mp.sqrt(index_data["rho"])
    dust_action = -(8 * mp.pi * MASS / 120) * 120 * mp.sqrt(index_data["rho"])
    action += dust_action
    for index in geometry["pole_indices"]:
        gradient[index] += complex(dust_gradient)
        common[index, index] += complex(dust_hessian)
        common_absolute_terms[index, index] += abs(complex(dust_hessian))
        common_term_counts[index, index] += 1

    matrices = {}
    assembly_roundoff = None
    maximum_imaginary = max(abs(mp.im(action)), float(np.max(np.abs(gradient.imag))))
    for name in DERIVATIVE_STEPS:
        matrix = common.copy()
        if name == "operational_primary":
            absolute_terms = common_absolute_terms.copy()
            term_counts = common_term_counts.copy()
        for simplex in geometry["simplex_records"]:
            simplex_indices = np.asarray(simplex["indices"], dtype=int)
            derivatives = pattern_cache[simplex["pattern"]]["derivatives"][name]
            for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
                area_record = triangle_local[triangle]
                rows = area_record["indices"]
                contribution = -1j * np.outer(
                    area_record["area_gradient"], derivatives[hinge_index]
                )
                matrix[np.ix_(rows, simplex_indices)] += contribution
                if name == "operational_primary":
                    absolute_terms[np.ix_(rows, simplex_indices)] += np.abs(contribution)
                    term_counts[np.ix_(rows, simplex_indices)] += 1
        maximum_imaginary = max(maximum_imaginary, float(np.max(np.abs(matrix.imag))))
        if name == "operational_primary":
            unit_roundoff = np.finfo(np.float64).eps / 2
            rounding_operations = term_counts.astype(np.float64) + 16
            gamma = (
                rounding_operations * unit_roundoff
                / (1 - rounding_operations * unit_roundoff)
            )
            envelope = gamma * absolute_terms
            assembly_roundoff = {
                "antisymmetry_bound": sparse_norm2(envelope + envelope.T),
                "maximum_entry_bound": float(np.max(envelope)),
                "maximum_term_count": int(np.max(term_counts)),
            }
        matrices[name] = matrix.real.copy()
        del matrix
    return gradient.real.copy(), matrices, {
        "action_real": mp.re(action),
        "action_imaginary": mp.im(action),
        "maximum_imaginary": maximum_imaginary,
        "curvature_maximum_imaginary": max(abs(mp.im(value)) for value in curvature),
        "binary64_roundoff": assembly_roundoff,
    }


def canonical_from_full(matrix):
    row_types = tuple(range(30, 65)) + tuple(range(0, 30))
    column_types = tuple(range(30, 65)) + tuple(range(65, 95))
    rows = np.concatenate([np.arange(24 * kind, 24 * (kind + 1)) for kind in row_types])
    columns = np.concatenate([
        np.arange(24 * kind, 24 * (kind + 1)) for kind in column_types
    ])
    result = matrix[np.ix_(rows, columns)].copy()
    result[35 * 24:, :] *= -1
    return result


def trivial_hessian(matrix):
    constant = np.ones((24, 1)) / math.sqrt(24)
    basis = np.kron(np.eye(95), constant)
    return basis.T @ matrix @ basis


def canonical_from_reduced(matrix):
    result = np.empty((65, 65), dtype=float)
    result[:35, :35] = matrix[30:65, 30:65]
    result[:35, 35:] = matrix[30:65, 65:95]
    result[35:, :35] = -matrix[0:30, 30:65]
    result[35:, 35:] = -matrix[0:30, 65:95]
    return result


def build_tangent_reduced(matrix, mapping):
    canonical = canonical_from_reduced(matrix)
    rhs = np.zeros((65, 60))
    rhs[:35, :30] = -matrix[30:65, 0:30]
    rhs[35:, :30] = matrix[0:30, 0:30]
    rhs[35:, 30:] = np.eye(30)
    solved = la.solve(canonical, rhs)
    raw = np.zeros((60, 60))
    raw[:30, :] = solved[35:, :]
    raw[30:, :] = (
        np.c_[matrix[65:95, 0:30], np.zeros((30, 30))]
        + matrix[65:95, 30:65] @ solved[:35, :]
        + matrix[65:95, 65:95] @ solved[35:, :]
    )
    permutation = np.asarray(mapping, dtype=int)
    return canonical, np.r_[raw[permutation, :], raw[30 + permutation, :]]


def sector_bases(index_data):
    table = index_data["table"]
    regular = []
    for element in range(24):
        matrix = np.zeros((24, 24), dtype=np.complex128)
        for column in range(24):
            matrix[int(table[element, column]), column] = 1
        regular.append(matrix)
    central = np.zeros((24, 24), dtype=np.complex128)
    for class_index, conjugacy_class in enumerate(index_data["classes"]):
        for element in conjugacy_class:
            central += (2**class_index) * regular[element]
    eigenvalues, eigenvectors = la.eig(central)
    centers = []
    for value in eigenvalues:
        if not any(abs(value - center) < 1e-7 for center in centers):
            centers.append(value)
    centers.sort(key=lambda value: (round(value.real, 8), round(value.imag, 8)))
    sectors = []
    constant = np.ones(24, dtype=np.complex128) / math.sqrt(24)
    for center in centers:
        indices = np.flatnonzero(np.abs(eigenvalues - center) < 1e-7)
        isotypic = la.qr(eigenvectors[:, indices], mode="economic")[0]
        dimension = int(round(math.sqrt(len(indices))))
        if dimension**2 != len(indices):
            raise RuntimeError("central eigenspace is not a square dimension")
        if dimension == 1:
            selected = isotypic
            splitter = None
        else:
            selected = None
            splitter = None
            for element in range(24):
                if element == index_data["identity"]:
                    continue
                hermitian = 1j * (regular[element] - regular[element].conj().T)
                restricted = isotypic.conj().T @ hermitian @ isotypic
                values, vectors = la.eigh((restricted + restricted.conj().T) / 2)
                clusters = []
                for position, value in enumerate(values):
                    if not clusters or abs(value - values[clusters[-1][0]]) > 1e-7:
                        clusters.append([position])
                    else:
                        clusters[-1].append(position)
                if len(clusters) == dimension and all(
                    len(cluster) == dimension for cluster in clusters
                ):
                    splitter = element
                    selected = isotypic @ vectors[:, clusters[0]]
                    break
            if selected is None:
                raise RuntimeError("deterministic representation splitter failed")
        selected = la.qr(selected, mode="economic")[0]
        constant_overlap = float(np.linalg.norm(selected.conj().T @ constant) ** 2)
        sectors.append({
            "central_eigenvalue": center,
            "isotypic_dimension": len(indices),
            "irrep_dimension": dimension,
            "splitter": splitter,
            "isotypic_basis": isotypic,
            "basis": selected,
            "constant_overlap": constant_overlap,
        })
    return sectors, {
        "central_normal_defect": float(np.linalg.norm(
            central @ central.conj().T - central.conj().T @ central
        )),
        "isotypic_dimensions": sorted(item["isotypic_dimension"] for item in sectors),
    }


def analyse_sector(sector, canonical_matrices):
    dimension = sector["irrep_dimension"]
    basis = np.kron(np.eye(65), sector["basis"])
    isotypic_basis = np.kron(np.eye(65), sector["isotypic_basis"])
    blocks = {
        name: basis.conj().T @ matrix @ basis
        for name, matrix in canonical_matrices.items()
    }
    operational = blocks["operational_primary"]
    block_image = canonical_matrices["operational_primary"] @ basis
    leakage_matrix = block_image - basis @ operational
    basis_leakage = norm2(leakage_matrix)
    isotypic_operational = (
        isotypic_basis.conj().T
        @ canonical_matrices["operational_primary"]
        @ isotypic_basis
    )
    isotypic_image = (
        canonical_matrices["operational_primary"] @ isotypic_basis
    )
    isotypic_leakage = norm2(
        isotypic_image - isotypic_basis @ isotypic_operational
    )
    orthonormal_error = float(np.linalg.norm(
        basis.conj().T @ basis - np.eye(65 * dimension), 2
    ))
    isotypic_orthonormal_error = float(np.linalg.norm(
        isotypic_basis.conj().T @ isotypic_basis
        - np.eye(65 * dimension**2), 2
    ))
    epsilon_basis = (
        max(basis_leakage, isotypic_leakage)
        + max(orthonormal_error, isotypic_orthonormal_error) * norm2(operational)
    )

    op_difference = blocks["operational_primary"] - blocks["operational_shadow"]
    val_difference = blocks["validation_primary"] - blocks["validation_shadow"]
    cross_difference = blocks["operational_primary"] - blocks["validation_primary"]
    epsilon_matrix = norm2(op_difference) + norm2(val_difference) + norm2(cross_difference)

    left, singular, right_h = la.svd(
        operational, full_matrices=False, lapack_driver="gesdd"
    )
    singular_gesvd = la.svd(
        operational, compute_uv=False, lapack_driver="gesvd"
    )
    isotypic_singular = la.svd(
        isotypic_operational, compute_uv=False, lapack_driver="gesvd"
    )
    expected_isotypic_singular = np.repeat(singular_gesvd, dimension)
    repetition_error = float(np.max(np.abs(
        isotypic_singular - expected_isotypic_singular
    )))
    driver_difference = float(np.max(np.abs(singular - singular_gesvd)))
    residuals = []
    for index, value in enumerate(singular):
        u = left[:, index]
        v = right_h[index, :].conj()
        residuals.append(float(
            np.linalg.norm(operational @ v - value * u)
            + np.linalg.norm(operational.conj().T @ u - value * v)
        ))
    epsilon_svd = max(residuals, default=0.0) + driver_difference
    epsilon_global = epsilon_matrix + epsilon_basis + epsilon_svd

    other_singular = {
        name: la.svd(matrix, compute_uv=False, lapack_driver="gesvd")
        for name, matrix in blocks.items() if name != "operational_primary"
    }
    resolved = singular > 100 * epsilon_global
    zero = np.ones(len(singular), dtype=bool)
    for index in range(len(singular)):
        estimates = [singular[index]] + [
            values[index] for values in other_singular.values()
        ]
        zero[index] = max(estimates) < 10 * epsilon_global
    open_flags = ~(resolved | zero)
    directional = []
    for index in range(max(0, len(singular) - 10), len(singular)):
        u = left[:, index]
        v = right_h[index, :].conj()
        directional.append({
            "index_descending": index,
            "singular_value": singular[index],
            "operational_proxy": abs(np.vdot(u, op_difference @ v)),
            "validation_proxy": abs(np.vdot(u, val_difference @ v)),
            "cross_proxy": abs(np.vdot(u, cross_difference @ v)),
            "svd_residual": residuals[index],
        })

    projection_singular = []
    if np.any(zero):
        null_rows = right_h[np.flatnonzero(zero), :]
        projection_singular = la.svd(
            null_rows[:, 35 * dimension:], compute_uv=False
        ).tolist()

    return {
        "blocks": blocks,
        "singular_values": singular,
        "resolved": resolved,
        "zero": zero,
        "open": open_flags,
        "epsilon_global": epsilon_global,
        "epsilon_matrix": epsilon_matrix,
        "epsilon_basis": epsilon_basis,
        "epsilon_svd": epsilon_svd,
        "basis_leakage": basis_leakage,
        "isotypic_leakage": isotypic_leakage,
        "orthonormal_error": orthonormal_error,
        "isotypic_orthonormal_error": isotypic_orthonormal_error,
        "repetition_error": repetition_error,
        "driver_difference": driver_difference,
        "directional": directional,
        "projection_singular": projection_singular,
        "relative_ranks": {
            f"{threshold:.0e}": int(np.sum(singular > threshold * singular[0]))
            for threshold in (1e-7, 1e-9, 1e-11, 1e-13, 1e-15)
        },
        "pseudoconstraints": int(np.sum(
            resolved & (singular < 1e-6 * singular[0])
        )),
    }


_DIRECT = None


def initialize_direct_worker(payload):
    global _DIRECT
    mp.mp.dps = DPS
    _DIRECT = payload


def full_gradient_at_delta(delta):
    payload = _DIRECT
    signed = [
        base * mp.exp(mp.mpf(str(delta[index])))
        for index, base in enumerate(payload["signed_base"])
    ]
    curvature = [
        mp.pi if record["boundary"] else 2 * mp.pi
        for record in payload["triangle_records"]
    ]
    minimum_minor = mp.inf
    minimum_argument = mp.inf
    negative_counts = Counter()
    for simplex in payload["simplex_records"]:
        values = [signed[index] for index in simplex["indices"]]
        angles, branch = angle_record(values)
        negative_counts[branch["negative_directions"]] += 1
        minimum_minor = min(minimum_minor, branch["minimum_leading_minor"])
        minimum_argument = min(minimum_argument, branch["minimum_argument"])
        for hinge_index, triangle in enumerate(simplex["hinge_triangles"]):
            curvature[triangle] += angles[hinge_index]
    gradient = [mp.mpc(0) for _ in range(2280)]
    action = mp.mpc(0)
    for triangle_index, record in enumerate(payload["triangle_records"]):
        values = [signed[index] for index in record["indices"]]
        area, area_gradient, _ = area_data(values)
        epsilon = curvature[triangle_index]
        action += -I * area * epsilon
        for index, derivative in zip(record["indices"], area_gradient):
            gradient[index] += -I * epsilon * derivative
    for index in payload["pole_indices"]:
        rho = -signed[index]
        action += -(8 * mp.pi * MASS / 120) * mp.sqrt(rho)
        gradient[index] += -(4 * mp.pi * MASS / 120) * mp.sqrt(rho)
    maximum_imaginary = max(
        abs(mp.im(action)), *(abs(mp.im(value)) for value in gradient)
    )
    return gradient, {
        "negative_counts": negative_counts,
        "minimum_leading_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
    }


def direct_worker(task):
    direction_index, step_text = task
    step = mp.mpf(step_text)
    direction = _DIRECT["directions"][direction_index]
    plus_delta = np.asarray(direction, dtype=float) * float(step)
    minus_delta = -plus_delta
    plus, plus_branch = full_gradient_at_delta(plus_delta)
    minus, minus_branch = full_gradient_at_delta(minus_delta)
    derivative = np.array([
        float(mp.re((left - right) / (2 * step)))
        for left, right in zip(plus, minus)
    ])
    branch_ok = bool(
        plus_branch["negative_counts"] == Counter({1: 2400})
        and minus_branch["negative_counts"] == Counter({1: 2400})
        and plus_branch["minimum_leading_minor"] > 0
        and minus_branch["minimum_leading_minor"] > 0
        and plus_branch["minimum_argument"] > mp.mpf("1e-6")
        and minus_branch["minimum_argument"] > mp.mpf("1e-6")
        and plus_branch["maximum_imaginary"] < mp.mpf("1e-60")
        and minus_branch["maximum_imaginary"] < mp.mpf("1e-60")
    )
    return {
        "direction": direction_index,
        "step": step_text,
        "derivative": derivative,
        "branch_ok": branch_ok,
        "minimum_leading_minor": str(min(
            plus_branch["minimum_leading_minor"], minus_branch["minimum_leading_minor"]
        )),
        "minimum_argument": str(min(
            plus_branch["minimum_argument"], minus_branch["minimum_argument"]
        )),
        "maximum_imaginary": str(max(
            plus_branch["maximum_imaginary"], minus_branch["maximum_imaginary"]
        )),
    }


def directional_controls(index_data, geometry, operational_hessian):
    directions = []
    names = []
    old_scale = np.zeros(2280)
    old_scale[:720] = 1 / math.sqrt(720)
    directions.append(old_scale)
    names.append("old_boundary_scale")

    internal_lapse = np.zeros(2280)
    rho = float(index_data["rho"])
    diagonal = float(mp.exp(mp.mpf(_DIRECT_STATE[0])) * L0_SQUARE - index_data["rho"])
    for index, kind in enumerate(index_data["edge_kind"]):
        if kind == "internal":
            internal_lapse[index] = -rho / diagonal
        elif kind == "pole":
            internal_lapse[index] = 1
    internal_lapse /= np.linalg.norm(internal_lapse)
    directions.append(internal_lapse)
    names.append("internal_collective_lapse")

    new_scale = np.zeros(2280)
    new_scale[1560:] = 1 / math.sqrt(720)
    directions.append(new_scale)
    names.append("new_boundary_scale")

    for name, orbit_type in (
        ("old_orbit_contrast", 0),
        ("internal_orbit_contrast", 30),
        ("new_orbit_contrast", 65),
    ):
        direction = np.zeros(2280)
        ordered = sorted(index_data["orbit_edges"][orbit_type])
        first = index_data["edge_to_index"][ordered[0]]
        second = index_data["edge_to_index"][ordered[1]]
        direction[first] = 1 / math.sqrt(2)
        direction[second] = -1 / math.sqrt(2)
        directions.append(direction)
        names.append(name)

    payload = {
        "signed_base": index_data["signed_base"],
        "triangle_records": geometry["triangle_records"],
        "simplex_records": geometry["simplex_records"],
        "pole_indices": geometry["pole_indices"],
        "directions": tuple(directions),
    }
    tasks = [
        (direction, mp.nstr(step, 20))
        for direction in range(6) for step in DIRECT_STEPS
    ]
    fork = mp_pool.get_context("fork")
    with fork.Pool(
        processes=min(8, len(tasks)), initializer=initialize_direct_worker,
        initargs=(payload,),
    ) as pool:
        raw = pool.map(direct_worker, tasks, chunksize=1)
    records = {(item["direction"], item["step"]): item for item in raw}
    controls = []
    all_pass = True
    for direction_index, (name, direction) in enumerate(zip(names, directions)):
        coarse = records[(direction_index, mp.nstr(DIRECT_STEPS[0], 20))]
        fine = records[(direction_index, mp.nstr(DIRECT_STEPS[1], 20))]
        richardson = (4 * fine["derivative"] - coarse["derivative"]) / 3
        assembled = operational_hessian @ direction
        relative = float(
            np.linalg.norm(richardson - assembled)
            / max(1.0, np.linalg.norm(richardson), np.linalg.norm(assembled))
        )
        passed_direction = bool(
            coarse["branch_ok"] and fine["branch_ok"] and relative < 2e-6
        )
        all_pass &= passed_direction
        controls.append({
            "name": name,
            "relative_error": relative,
            "branch_ok": bool(coarse["branch_ok"] and fine["branch_ok"]),
            "passed": passed_direction,
            "richardson_change": float(np.linalg.norm(richardson - fine["derivative"])),
            "minimum_leading_minor": min(
                coarse["minimum_leading_minor"], fine["minimum_leading_minor"],
                key=mp.mpf,
            ),
            "minimum_argument": min(
                coarse["minimum_argument"], fine["minimum_argument"], key=mp.mpf,
            ),
            "maximum_imaginary": max(
                coarse["maximum_imaginary"], fine["maximum_imaginary"], key=mp.mpf,
            ),
        })
    return bool(all_pass), controls


def serialize_float(value):
    return f"{float(value):.17e}"


print("=" * 78)
print("FULL ANISOTROPIC 600-CELL DUST CANONICAL-RANK CENSUS")
print("=" * 78)

records = {}
global_controls = True

for parity in ("even", "odd"):
    print(f"[{parity}] preparing the complete 2280-edge carrier", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    _DIRECT_STATE = state
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)

    carrier_ok = bool(
        len(model["slab"]) == 2400
        and len(model["old_edges"]) == 720
        and len(model["internal_edges"]) == 840
        and len(model["new_edges"]) == 720
        and len(model["faces"][3]) == 6240
        and len(index_data["edge_to_index"]) == 2280
        and Counter(index_data["orders"]) == Counter({6: 8, 3: 8, 4: 6, 2: 1, 1: 1})
        and sorted(map(len, index_data["classes"])) == [1, 1, 4, 4, 4, 4, 6]
        and len(geometry["patterns"]) == 20
    )
    check(
        f"{parity}: the full carrier and target-free 2T census match the protocol",
        carrier_ok,
        f"edges=2280, triangles={len(model['faces'][3])}, patterns={len(geometry['patterns'])}",
    )

    s = mp.mpf(state[0])
    kind_values = {
        "old": L0_SQUARE,
        "internal": mp.exp(s) * L0_SQUARE - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s) * L0_SQUARE,
    }
    print(f"[{parity}] differentiating {len(geometry['patterns'])} local simplex patterns", flush=True)
    pattern_cache, local_control = build_pattern_cache(geometry["patterns"], kind_values)
    local_ok = bool(
        local_control["entry_pass"]
        and local_control["base_negative_counts"] == Counter({1: 2400})
        and local_control["all_displaced_have_one_negative"]
        and local_control["minimum_leading_minor"] > 0
        and local_control["minimum_argument"] > mp.mpf("1e-6")
    )
    check(
        f"{parity}: all local angle derivatives pass the frozen branch and step gates",
        local_ok,
        f"cross={mp.nstr(local_control['maximum_cross'], 7)}, "
        f"op={mp.nstr(local_control['maximum_operational_proxy'], 7)}, "
        f"val={mp.nstr(local_control['maximum_validation_proxy'], 7)}",
    )

    print(f"[{parity}] assembling four complete 2280 x 2280 Hessians", flush=True)
    gradient, hessians, assembly = assemble_full_hessians(
        model, index_data, geometry, pattern_cache
    )
    op = hessians["operational_primary"]
    op_shadow = hessians["operational_shadow"]
    val = hessians["validation_primary"]
    val_shadow = hessians["validation_shadow"]
    cross = np.abs(op - val)
    proxy = np.abs(op - op_shadow) + np.abs(val - val_shadow) + 1e-70
    entry_ok = bool(np.all(cross <= 10 * proxy))

    base_payload = {
        "signed_base": index_data["signed_base"],
        "triangle_records": geometry["triangle_records"],
        "simplex_records": geometry["simplex_records"],
        "pole_indices": geometry["pole_indices"],
        "directions": (),
    }
    _DIRECT = base_payload
    base_gradient_mp, base_branch = full_gradient_at_delta(np.zeros(2280))
    internal_mp = base_gradient_mp[720:1560]
    old_mp = base_gradient_mp[:720]
    new_mp = base_gradient_mp[1560:]
    stored_pre_mp = tuple(
        mp.mpf(value) for value in tick["solutions"][parity]["pre_momentum"]
    )
    stored_post_mp = tuple(
        mp.mpf(value) for value in tick["solutions"][parity]["post_momentum"]
    )
    orbit_pre_mp = tuple(
        mp.fsum(-old_mp[24 * orbit + group] for group in range(24)) / 24
        for orbit in range(30)
    )
    orbit_post_mp = tuple(
        mp.fsum(new_mp[24 * orbit + group] for group in range(24)) / 24
        for orbit in range(30)
    )
    orbit_spread_mp = mp.mpf(0)
    for values, count in ((old_mp, 30), (internal_mp, 35), (new_mp, 30)):
        for orbit in range(count):
            block = values[24 * orbit:24 * (orbit + 1)]
            mean = mp.fsum(block) / 24
            orbit_spread_mp = max(
                orbit_spread_mp, *(abs(value - mean) for value in block)
            )
    internal_maximum_mp = max(abs(value) for value in internal_mp)
    momentum_error_mp = max(
        *(abs(value - target) for value, target in zip(orbit_pre_mp, stored_pre_mp)),
        *(abs(value - target) for value, target in zip(orbit_post_mp, stored_post_mp)),
    )
    assembly_gradient_error = max(
        abs(gradient[index] - float(mp.re(value)))
        for index, value in enumerate(base_gradient_mp)
    )
    gradient_ok = bool(
        internal_maximum_mp < mp.mpf("1e-25")
        and momentum_error_mp < mp.mpf("1e-40")
        and orbit_spread_mp < mp.mpf("1e-40")
        and assembly_gradient_error < 2e-11
        and base_branch["negative_counts"] == Counter({1: 2400})
        and base_branch["minimum_leading_minor"] > 0
        and base_branch["minimum_argument"] > mp.mpf("1e-6")
        and base_branch["maximum_imaginary"] < mp.mpf("1e-60")
        and assembly["maximum_imaginary"] < 1e-60
    )
    check(
        f"{parity}: all 840 internal equations and 1440 boundary momenta reproduce edgewise",
        gradient_ok,
        f"internal={mp.nstr(internal_maximum_mp, 6)}, "
        f"momentum={mp.nstr(momentum_error_mp, 6)}, "
        f"spread={mp.nstr(orbit_spread_mp, 6)}, assembly={assembly_gradient_error:.3e}",
    )

    antisymmetric = op - op.T
    derivative_error = (
        sparse_norm2(op - op_shadow)
        + sparse_norm2(val - val_shadow)
        + sparse_norm2(op - val)
    )
    antisymmetric_norm = sparse_norm2(antisymmetric)
    assembly_roundoff_bound = assembly["binary64_roundoff"]["antisymmetry_bound"]
    reciprocity_allowance = 10 * derivative_error + assembly_roundoff_bound
    reciprocity_ok = bool(antisymmetric_norm <= reciprocity_allowance)
    check(
        f"{parity}: the full Hessian is reciprocal inside its calibrated derivative error",
        reciprocity_ok and entry_ok,
        f"antisym={antisymmetric_norm:.3e}, derivative={derivative_error:.3e}, "
        f"roundoff={assembly_roundoff_bound:.3e}, entry={entry_ok}",
    )

    reduced_hessians = {name: trivial_hessian(matrix) for name, matrix in hessians.items()}
    mapping = tuple(gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"])
    reduced_canonical, reduced_tangent = build_tangent_reduced(
        reduced_hessians["operational_primary"], mapping
    )
    stored_canonical_singular = np.asarray([
        float(value) for value in tangent_input["parities"][parity]["canonical_singular_values"]
    ])
    current_canonical_singular = la.svd(
        reduced_canonical, compute_uv=False, lapack_driver="gesvd"
    )
    singular_error = float(np.max(
        np.abs(current_canonical_singular - stored_canonical_singular)
        / np.maximum(1.0, np.abs(stored_canonical_singular))
    ))
    stored_tangent = np.asarray(
        [[float(value) for value in row]
         for row in tangent_input["parities"][parity]["tangent_matrix"]]
    )
    tangent_error = float(
        np.linalg.norm(reduced_tangent - stored_tangent)
        / max(1.0, np.linalg.norm(stored_tangent), np.linalg.norm(reduced_tangent))
    )
    reduced_ok = bool(singular_error < 2e-8 and tangent_error < 2e-8)
    check(
        f"{parity}: the invariant restriction reproduces the committed reduced map",
        reduced_ok,
        f"singular={singular_error:.3e}, tangent={tangent_error:.3e}",
    )

    print(f"[{parity}] running six preregistered full-action directional controls", flush=True)
    directions_ok, direction_records = directional_controls(index_data, geometry, op)
    check(
        f"{parity}: all six independent full-gradient directions reproduce the Hessian",
        directions_ok,
        "max relative=" + f"{max(item['relative_error'] for item in direction_records):.3e}",
    )

    sectors, sector_geometry = sector_bases(index_data)
    sector_geometry_ok = bool(
        sector_geometry["central_normal_defect"] < 1e-10
        and sector_geometry["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sorted(item["irrep_dimension"] for item in sectors) == [1, 1, 1, 2, 2, 2, 3]
        and sum(65 * item["irrep_dimension"] ** 2 for item in sectors) == 1560
        and sum(item["constant_overlap"] > 1 - 1e-12 for item in sectors) == 1
    )
    canonical_matrices = {
        name: canonical_from_full(matrix) for name, matrix in hessians.items()
    }
    print(f"[{parity}] evaluating the seven frozen 2T canonical blocks", flush=True)
    sector_records = [analyse_sector(sector, canonical_matrices) for sector in sectors]
    sector_control_ok = bool(
        sector_geometry_ok
        and all(
            max(record["basis_leakage"], record["isotypic_leakage"])
            <= 10 * (record["epsilon_matrix"] + record["epsilon_svd"] + 1e-70)
            and record["repetition_error"] <= 10 * (
                record["epsilon_matrix"]
                + record["epsilon_basis"]
                + record["epsilon_svd"]
                + 1e-70
            )
            for record in sector_records
        )
    )
    check(
        f"{parity}: seven deterministic 2T blocks exhaust the complete canonical matrix",
        sector_control_ok,
        f"dims={[65*item['irrep_dimension'] for item in sectors]}, "
        f"max leakage={max(max(item['basis_leakage'], item['isotypic_leakage']) for item in sector_records):.3e}, "
        f"max repetition={max(item['repetition_error'] for item in sector_records):.3e}",
    )

    controls_ok = bool(
        carrier_ok and local_ok and gradient_ok and entry_ok and reciprocity_ok
        and reduced_ok and directions_ok and sector_control_ok
    )
    any_open = any(np.any(record["open"]) for record in sector_records)
    any_zero = any(np.any(record["zero"]) for record in sector_records)
    all_resolved = all(np.all(record["resolved"]) for record in sector_records)
    if not controls_ok:
        outcome = "FULL_CANONICAL_ASSEMBLY_CONTROL_FAILED"
    elif any_open:
        outcome = "FULL_CANONICAL_RANK_NUMERICALLY_OPEN"
    elif all_resolved:
        outcome = "FULL_CANONICAL_LEGENDRE_REGULAR"
    elif any_zero:
        outcome = "FULL_CANONICAL_DEGENERACY"
    else:
        outcome = "FULL_CANONICAL_RANK_NUMERICALLY_OPEN"
    full_rank = sum(
        int(np.sum(record["resolved"])) * sector["irrep_dimension"]
        for sector, record in zip(sectors, sector_records)
    )
    full_nullity = sum(
        int(np.sum(record["zero"])) * sector["irrep_dimension"]
        for sector, record in zip(sectors, sector_records)
    )
    full_open = sum(
        int(np.sum(record["open"])) * sector["irrep_dimension"]
        for sector, record in zip(sectors, sector_records)
    )
    check(
        f"{parity}: the frozen hierarchy assigns the complete rank outcome",
        outcome in {
            "FULL_CANONICAL_ASSEMBLY_CONTROL_FAILED",
            "FULL_CANONICAL_RANK_NUMERICALLY_OPEN",
            "FULL_CANONICAL_LEGENDRE_REGULAR",
            "FULL_CANONICAL_DEGENERACY",
        },
        f"outcome={outcome}, rank={full_rank}, zero={full_nullity}, open={full_open}",
    )
    global_controls &= controls_ok

    serialized_sectors = []
    for sector, record in zip(sectors, sector_records):
        serialized_sectors.append({
            "central_eigenvalue": {
                "real": serialize_float(sector["central_eigenvalue"].real),
                "imaginary": serialize_float(sector["central_eigenvalue"].imag),
            },
            "isotypic_dimension": sector["isotypic_dimension"],
            "irrep_dimension": sector["irrep_dimension"],
            "minimal_block_size": 65 * sector["irrep_dimension"],
            "splitter_group_index": sector["splitter"],
            "constant_overlap": serialize_float(sector["constant_overlap"]),
            "singular_values": [serialize_float(value) for value in record["singular_values"]],
            "resolved_count_minimal": int(np.sum(record["resolved"])),
            "zero_count_minimal": int(np.sum(record["zero"])),
            "open_count_minimal": int(np.sum(record["open"])),
            "epsilon_global": serialize_float(record["epsilon_global"]),
            "epsilon_matrix": serialize_float(record["epsilon_matrix"]),
            "epsilon_basis": serialize_float(record["epsilon_basis"]),
            "epsilon_svd": serialize_float(record["epsilon_svd"]),
            "basis_leakage": serialize_float(record["basis_leakage"]),
            "isotypic_leakage": serialize_float(record["isotypic_leakage"]),
            "singular_repetition_error": serialize_float(record["repetition_error"]),
            "condition": serialize_float(
                record["singular_values"][0] / record["singular_values"][-1]
                if record["singular_values"][-1] else math.inf
            ),
            "relative_ranks": record["relative_ranks"],
            "pseudoconstraint_candidates_minimal": record["pseudoconstraints"],
            "new_boundary_projection_singular_values": [
                serialize_float(value) for value in record["projection_singular"]
            ],
            "weak_direction_diagnostics": [
                {key: value if isinstance(value, int) else serialize_float(value)
                 for key, value in item.items()}
                for item in record["directional"]
            ],
        })
    records[parity] = {
        "outcome": outcome,
        "controls_pass": controls_ok,
        "full_resolved_rank": full_rank,
        "full_error_consistent_nullity": full_nullity,
        "full_numerically_open_count": full_open,
        "pseudoconstraint_candidates": sum(
            record["pseudoconstraints"] * sector["irrep_dimension"]
            for sector, record in zip(sectors, sector_records)
        ),
        "carrier": {
            "edge_variables": 2280,
            "old": 720,
            "internal": 840,
            "new": 720,
            "four_simplices": 2400,
            "triangles": 6240,
            "simplex_patterns": len(geometry["patterns"]),
            "group_orders": dict(sorted(Counter(index_data["orders"]).items())),
            "conjugacy_class_sizes": sorted(map(len, index_data["classes"])),
        },
        "assembly": {
            "maximum_internal_residual": mp.nstr(internal_maximum_mp, 40),
            "maximum_momentum_error": mp.nstr(momentum_error_mp, 40),
            "maximum_orbit_gradient_spread": mp.nstr(orbit_spread_mp, 40),
            "binary_assembly_gradient_error": serialize_float(assembly_gradient_error),
            "maximum_imaginary": serialize_float(assembly["maximum_imaginary"]),
            "entrywise_calibration_pass": entry_ok,
            "antisymmetric_spectral_norm": serialize_float(antisymmetric_norm),
            "derivative_error_spectral": serialize_float(derivative_error),
            "binary64_assembly_roundoff_antisymmetry_bound": serialize_float(
                assembly_roundoff_bound
            ),
            "binary64_assembly_roundoff_maximum_entry_bound": serialize_float(
                assembly["binary64_roundoff"]["maximum_entry_bound"]
            ),
            "binary64_assembly_maximum_term_count": assembly[
                "binary64_roundoff"
            ]["maximum_term_count"],
            "reduced_singular_error": serialize_float(singular_error),
            "reduced_tangent_error": serialize_float(tangent_error),
        },
        "directional_controls": direction_records,
        "sector_geometry": sector_geometry,
        "sectors": serialized_sectors,
    }

if not global_controls:
    combined_outcome = "FULL_CANONICAL_ASSEMBLY_CONTROL_FAILED"
elif any(item["outcome"] == "FULL_CANONICAL_RANK_NUMERICALLY_OPEN" for item in records.values()):
    combined_outcome = "FULL_CANONICAL_RANK_NUMERICALLY_OPEN"
elif all(item["outcome"] == "FULL_CANONICAL_LEGENDRE_REGULAR" for item in records.values()):
    combined_outcome = "FULL_CANONICAL_LEGENDRE_REGULAR"
elif any(item["outcome"] == "FULL_CANONICAL_DEGENERACY" for item in records.values()):
    combined_outcome = "FULL_CANONICAL_DEGENERACY"
else:
    combined_outcome = "FULL_CANONICAL_RANK_NUMERICALLY_OPEN"

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "clarification_commit": CLARIFICATION_COMMIT,
    "rounding_correction_commit": ROUNDING_CORRECTION_COMMIT,
    "input_sha256": hashes,
    "full_720_edge_boundary": True,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "tangent_map_constructed": False,
    "outcome": combined_outcome,
    "parities": records,
    "classification": {
        "rank": "DERIVED COMPUTATIONAL" if combined_outcome == "FULL_CANONICAL_LEGENDRE_REGULAR" else "OPEN_OR_DEGENERATE",
        "pseudoconstraints": "STRUCTURAL CANDIDATES ONLY",
        "gauge_identification": "OPEN",
        "graviton_interpretation": "OPEN",
        "limiting_speed": "OPEN",
        "refinement": "OPEN",
    },
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Summary: {passed}/{tests} checks passed")
print(f"Outcome: {combined_outcome}")
for parity, record in records.items():
    print(
        f"  {parity}: rank={record['full_resolved_rank']}, "
        f"zero={record['full_error_consistent_nullity']}, "
        f"open={record['full_numerically_open_count']}, "
        f"pseudo={record['pseudoconstraint_candidates']}"
    )
print(f"Artifact: {OUTPUT}")

if passed != tests:
    sys.exit(1)
