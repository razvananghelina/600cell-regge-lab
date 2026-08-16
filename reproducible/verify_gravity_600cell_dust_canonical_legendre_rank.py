#!/usr/bin/env python3
"""Canonical pre-Legendre rank census for the 600-cell dust slab.

Prior-art commit: 31b1690.
Protocol commit: 4837a16; tolerance clarification: be145e4.
No nonlinear root search is performed.
"""

from collections import Counter
import contextlib
import importlib.util
import io
from itertools import combinations
import json
import math
import multiprocessing as mp_pool
from pathlib import Path
import sys

import mpmath as arb
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_dust_canonical_legendre_rank.json"
PRIOR_ART_COMMIT = "31b1690"
PROTOCOL_COMMIT = "4837a16"
TOLERANCE_COMMIT = "be145e4"
CLARIFICATION_COMMIT = "92182e1"
PROJECTION_RULE_COMMIT = "9684918"
GLUING_RESULT_COMMIT = "a766740"
GAUGE_INPUT = HERE / "gravity_600cell_dust_gauge_quotient_precision.json"
GLUING_INPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
PUBLISHED_INPUT = HERE / "gravity_600cell_published_dust_control.json"
DPS = 100
arb.mp.dps = DPS
DERIVATIVE_STEPS = {
    "operational": (arb.mpf("1e-20"), arb.mpf("1e-15")),
    "validation": (arb.mpf("3e-20"), arb.mpf("3e-15")),
}
ARITHMETIC_FLOOR = arb.mpf("1e-70")
ENTRY_GATE_FACTOR = arb.mpf(10)
NONZERO_FACTOR = arb.mpf(100)
ZERO_FACTOR = arb.mpf(10)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def relative_error(left, right):
    return abs(left-right)/max(arb.mpf(1), abs(left), abs(right))


spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_canonical_rank",
    HERE / "verify_gravity_global_regge_orbits.py",
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
    "the imported one-slab geometric/action core retains all 43 certificates",
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
        "old_lookup": {
            edge: index for index, orbit in enumerate(old_orbits) for edge in orbit
        },
        "final_lookup": {
            edge: index
            for index, orbit in enumerate(final_orbits) for edge in orbit
        },
    }


models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}

for parity, model in models.items():
    check(
        f"{parity}: the canonical carrier has 30 old, 35 internal and 30 new orbits",
        len(model["old_orbits"]) == 30
        and len(model["edge_orbits"]) == 35
        and len(model["pole_orbits"]) == 5
        and len(model["final_orbits"]) == 30
        and Counter(map(len, model["old_orbits"])) == Counter({24: 30})
        and Counter(map(len, model["edge_orbits"])) == Counter({24: 35})
        and Counter(map(len, model["final_orbits"])) == Counter({24: 30}),
    )


gauge_input = json.loads(GAUGE_INPUT.read_text())
gluing_input = json.loads(GLUING_INPUT.read_text())
published_input = json.loads(PUBLISHED_INPUT.read_text())
check(
    "all committed control artifacts have the required provenance and parities",
    gauge_input["protocol_commit"] == "da34272"
    and gluing_input["protocol_commit"] == "29dcfa5"
    and gluing_input["branch_precision_correction_commit"] == "ab75d91"
    and gluing_input["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and published_input["protocol_commit"] == "cc0902b"
    and set(gauge_input["parities"])
        == set(gluing_input["parities"])
        == set(published_input["parities"])
        == {"even", "odd"},
)


ARB_I = arb.mpc(0, 1)
ARB_M_STAR = arb.mpf(10)
ARB_ZETA = (arb.pi**2*arb.sqrt(2)/50)**(arb.mpf(1)/3)
ARB_R0 = 4*ARB_M_STAR/(3*arb.pi)
ARB_L0 = ARB_ZETA*ARB_R0
ARB_L0_SQUARE = ARB_L0**2
ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
ARB_MASS = (90/arb.pi)*ARB_EPSILON_3*ARB_L0
ARB_TAU = arb.mpf("0.0102")
ARB_RHO = ARB_TAU**2
ARB_SLANT_SQUARE = ARB_L0_SQUARE-ARB_RHO
ARB_BASE_OLD = tuple(ARB_L0_SQUARE for _ in range(30))
ARB_BASE_X = tuple(
    [ARB_SLANT_SQUARE for _ in range(30)]
    + [ARB_RHO for _ in range(5)]
)
ARB_BASE_NEW = tuple(ARB_L0_SQUARE for _ in range(30))


def arb_log_minus(value):
    scale = max(arb.mpf(1), abs(value))
    if abs(arb.im(value)) < arb.mpf("1e-80")*scale:
        real = arb.re(value)
        if real < 0:
            return arb.log(-real)-ARB_I*arb.pi
        return arb.log(real)
    return arb.log(value)


def arb_signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices)-1
    if dimension == 0:
        return arb.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = arb.matrix([
        [
            (squared[base][left]+squared[base][right]-squared[left][right])/2
            for right in others
        ]
        for left in others
    ])
    return arb.det(gram)/(arb.factorial(dimension)**2)


def arb_angle_record(squared):
    gram = arb.matrix([
        [
            (squared[0][left]+squared[0][right]-squared[left][right])/2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram**-1
    simplex_volume_square = arb_signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: arb_signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    leading_minors = []
    for size in range(1, 5):
        principal = arb.matrix([
            [gram[left, right] for right in range(size)]
            for left in range(size)
        ])
        leading_minors.append(arb.det(principal))
    signs = [1]
    for value in leading_minors:
        signs.append(1 if value > 0 else -1 if value < 0 else 0)
    negative_directions = None if 0 in signs else sum(
        left != right for left, right in zip(signs, signs[1:])
    )

    angles = {}
    minimum_argument = arb.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge_vertices = tuple(
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = arb_signed_volume_square(squared, hinge_vertices)
        gram_derivative = arb.matrix(4, 4)
        opposite_edge = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                gram_derivative[left-1, right-1] = (
                    int({0, left} == opposite_edge)
                    + int({0, right} == opposite_edge)
                    - int(left != right and {left, right} == opposite_edge)
                )/2
        product = inverse*gram_derivative
        volume_derivative = simplex_volume_square*sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            arb.sqrt(arb.mpc(facet_volume_squares[omitted_a]))
            * arb.sqrt(arb.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16*volume_derivative/denominator
        sine = -arb.mpf(4)/3*(
            arb.sqrt(arb.mpc(hinge_volume_square))
            * arb.sqrt(arb.mpc(simplex_volume_square))
        )/denominator
        argument = cosine+ARB_I*sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[hinge_vertices] = -ARB_I*arb_log_minus(argument)
    return {
        "angles": angles,
        "negative_directions": negative_directions,
        "minimum_leading_minor": min(abs(value) for value in leading_minors),
        "minimum_argument": minimum_argument,
    }


def triangle_area_square(values):
    x, y, z = values
    return (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16


def triangle_area_square_partials(values):
    x, y, z = values
    return ((y+z-x)/8, (x+z-y)/8, (x+y-z)/8)


def edge_data(model, edge, q_old, x, q_new):
    edge = tuple(sorted(edge))
    if edge in model["old_lookup"]:
        index = model["old_lookup"][edge]
        return q_old[index], index
    if edge in model["edge_to_variable"]:
        internal = model["edge_to_variable"][edge]
        value = int(model["edge_jacobian"][edge])*x[internal]
        return value, 30+internal
    if edge in model["final_lookup"]:
        index = model["final_lookup"][edge]
        return q_new[index], 65+index
    raise ValueError(f"edge absent from 95-variable slab: {edge}")


def simplex_squared(model, simplex, q_old, x, q_new):
    squared = [[arb.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value, _ = edge_data(
            model, (simplex[left], simplex[right]), q_old, x, q_new
        )
        squared[left][right] = squared[right][left] = value
    return squared


def action_and_gradient(model, q_old, x, q_new):
    curvature = [
        arb.pi if min(orbit) in model["boundary_triangles"] else 2*arb.pi
        for orbit in model["triangle_orbits"]
    ]
    negative_counts = Counter()
    minimum_leading_minor = arb.inf
    minimum_argument = arb.inf
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        record = arb_angle_record(simplex_squared(model, simplex, q_old, x, q_new))
        negative_counts[record["negative_directions"]] += len(simplex_orbit)
        minimum_leading_minor = min(
            minimum_leading_minor, record["minimum_leading_minor"]
        )
        minimum_argument = min(minimum_argument, record["minimum_argument"])
        for local_hinge, angle in record["angles"].items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            curvature[model["triangle_to_orbit"][triangle]] += angle

    action_sum = arb.mpc(0)
    gradient = [arb.mpc(0) for _ in range(95)]
    for orbit_index, orbit in enumerate(model["triangle_orbits"]):
        triangle = min(orbit)
        edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
        data = tuple(edge_data(model, edge, q_old, x, q_new) for edge in edges)
        values = tuple(item[0] for item in data)
        area_square = triangle_area_square(values)
        area = arb.sqrt(arb.mpc(area_square))
        triangle_curvature = curvature[orbit_index]
        multiplicity = arb.mpf(len(orbit))
        action_sum += multiplicity*area*triangle_curvature
        for (squared_value, variable), partial in zip(
            data, triangle_area_square_partials(values)
        ):
            gradient[variable] += (
                -ARB_I*(multiplicity/24)*triangle_curvature
                * partial*squared_value/(2*area)
            )
    gravitational = -ARB_I*action_sum
    dust = -(8*arb.pi*ARB_MASS/5)*sum(arb.sqrt(value) for value in x[30:35])
    for internal in range(30, 35):
        gradient[30+internal] += (
            -(4*arb.pi*ARB_MASS/5)*arb.sqrt(x[internal])/24
        )
    return gravitational+dust, tuple(gradient), {
        "negative_counts": negative_counts,
        "minimum_leading_minor": minimum_leading_minor,
        "minimum_argument": minimum_argument,
        "representative_simplices": len(model["simplex_orbits"]),
    }


def perturb_base(coordinate, delta):
    q_old = list(ARB_BASE_OLD)
    x = list(ARB_BASE_X)
    q_new = list(ARB_BASE_NEW)
    if coordinate < 30:
        q_old[coordinate] *= arb.exp(delta)
    elif coordinate < 65:
        x[coordinate-30] *= arb.exp(delta)
    else:
        q_new[coordinate-65] *= arb.exp(delta)
    return tuple(q_old), tuple(x), tuple(q_new)


def pack_complex(value):
    return arb.nstr(arb.re(value), DPS+5), arb.nstr(arb.im(value), DPS+5)


def unpack_complex(value):
    return arb.mpc(arb.mpf(value[0]), arb.mpf(value[1]))


def pack_branch(record):
    return {
        "negative_counts": {
            "NONE" if key is None else str(key): value
            for key, value in record["negative_counts"].items()
        },
        "minimum_leading_minor": arb.nstr(record["minimum_leading_minor"], DPS+5),
        "minimum_argument": arb.nstr(record["minimum_argument"], DPS+5),
        "representative_simplices": record["representative_simplices"],
    }


def unpack_branch(record):
    return {
        "negative_counts": Counter({
            None if key == "NONE" else int(key): value
            for key, value in record["negative_counts"].items()
        }),
        "minimum_leading_minor": arb.mpf(record["minimum_leading_minor"]),
        "minimum_argument": arb.mpf(record["minimum_argument"]),
        "representative_simplices": record["representative_simplices"],
    }


_WORKER_MODEL = None


def initialize_worker(model):
    global _WORKER_MODEL
    arb.mp.dps = DPS
    _WORKER_MODEL = model


def gradient_worker(task):
    coordinate, delta_text = task
    q_old, x, q_new = perturb_base(coordinate, arb.mpf(delta_text))
    action, gradient, branch = action_and_gradient(_WORKER_MODEL, q_old, x, q_new)
    return {
        "action": pack_complex(action),
        "gradient": [pack_complex(value) for value in gradient],
        "branch": pack_branch(branch),
    }


def matrix_from_gradients(values, h):
    matrix = arb.matrix(95, 95)
    positive = arb.nstr(h, 20)
    negative = arb.nstr(-h, 20)
    for column in range(95):
        plus = values[(column, positive)]
        minus = values[(column, negative)]
        for row in range(95):
            matrix[row, column] = arb.re((plus[row]-minus[row])/(2*h))
    return matrix


def matrix_sub(left, right):
    return left-right


def max_abs_entry(matrix):
    return max(abs(matrix[row, column]) for row in range(matrix.rows)
               for column in range(matrix.cols))


def frobenius_norm(matrix):
    return arb.sqrt(sum(abs(matrix[row, column])**2
                        for row in range(matrix.rows)
                        for column in range(matrix.cols)))


def spectral_norm(matrix):
    singular = arb.svd_r(matrix, compute_uv=False)
    return singular[0] if singular.rows else arb.mpf(0)


def extract_canonical(matrix):
    result = arb.matrix(65, 65)
    # Rows 0:35 are internal equations; rows 35:65 are p_pre=-g_old.
    # Columns 0:35 are internal unknowns; columns 35:65 are q_new.
    for row in range(35):
        for column in range(35):
            result[row, column] = matrix[30+row, 30+column]
        for column in range(30):
            result[row, 35+column] = matrix[30+row, 65+column]
    for row in range(30):
        for column in range(35):
            result[35+row, column] = -matrix[row, 30+column]
        for column in range(30):
            result[35+row, 35+column] = -matrix[row, 65+column]
    return result


def matrix_to_numpy(matrix):
    return np.array([
        [float(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ])


def relative_frobenius(left, right):
    return np.linalg.norm(left-right)/max(1.0, np.linalg.norm(left), np.linalg.norm(right))


def householder_complement(vector):
    vector = np.asarray(vector, dtype=float)
    vector = vector/np.linalg.norm(vector)
    sign = 1.0 if vector[0] >= 0 else -1.0
    reflector_vector = vector.copy()
    reflector_vector[0] += sign
    reflector = (
        np.eye(len(vector))
        - 2*np.outer(reflector_vector, reflector_vector)
        / float(reflector_vector@reflector_vector)
    )
    return reflector[:, 1:]


raw_lapse = np.r_[
    np.full(30, -float(ARB_RHO/ARB_SLANT_SQUARE)),
    np.ones(5),
]
lapse = raw_lapse/np.linalg.norm(raw_lapse)
internal_complement = householder_complement(lapse)
scale = np.ones(30)/math.sqrt(30)
shape_complement = householder_complement(scale)
internal_basis = np.column_stack((lapse, internal_complement))
boundary_basis = np.column_stack((scale, shape_complement))
sector_basis = np.zeros((65, 65))
sector_basis[:35, :35] = internal_basis
sector_basis[35:, 35:] = boundary_basis
check(
    "the frozen lapse/scale Householder bases are orthonormal",
    np.linalg.norm(sector_basis.T@sector_basis-np.eye(65), 2) < 3e-14,
)


def sector_overlaps(vector):
    coefficients = sector_basis.T@np.asarray(vector, dtype=float)
    return {
        "internal_lapse": float(np.sum(coefficients[0:1]**2)),
        "internal_transverse": float(np.sum(coefficients[1:35]**2)),
        "boundary_scale": float(np.sum(coefficients[35:36]**2)),
        "boundary_shape": float(np.sum(coefficients[36:65]**2)),
    }


def block_norms(matrix):
    transformed = sector_basis.T@matrix@sector_basis
    cuts = ((0, 1), (1, 35), (35, 36), (36, 65))
    names = ("lapse", "transverse", "scale", "shape")
    return {
        row_name: {
            column_name: float(np.linalg.norm(
                transformed[row_start:row_end, column_start:column_end]
            ))
            for column_name, (column_start, column_end) in zip(names, cuts)
        }
        for row_name, (row_start, row_end) in zip(names, cuts)
    }


def mp_bilinear(u, matrix, v):
    return abs((u.T*matrix*v)[0])


def complex_record(value, digits=70):
    return {
        "real": arb.nstr(arb.re(value), digits),
        "imaginary": arb.nstr(arb.im(value), 20),
    }


print("="*78)
print("CANONICAL PRE-LEGENDRE RANK CENSUS")
print("="*78)

fork_context = mp_pool.get_context("fork")
records = {}

for parity in ("even", "odd"):
    model = models[parity]
    print(f"[{parity}] evaluating the 100-decimal base gradient", flush=True)
    base_action, base_gradient, base_branch = action_and_gradient(
        model, ARB_BASE_OLD, ARB_BASE_X, ARB_BASE_NEW
    )
    maximum_imaginary = max(
        abs(arb.im(base_action)), *(abs(arb.im(value)) for value in base_gradient)
    )

    published_action = arb.mpc(
        str(published_input["parities"][parity]["total_action_real"]),
        str(published_input["parities"][parity]["total_action_imaginary"]),
    )
    action_error = relative_error(base_action, published_action)
    internal_residuals = base_gradient[30:65]
    internal_maximum = max(abs(value) for value in internal_residuals)
    stored_pre = tuple(
        arb.mpf(item["real"])
        for item in gluing_input["parities"][parity]["momenta"]["pre"]
    )
    stored_post = tuple(
        arb.mpf(item["real"])
        for item in gluing_input["parities"][parity]["momenta"]["post"]
    )
    pre = tuple(-base_gradient[index] for index in range(30))
    post = tuple(base_gradient[65+index] for index in range(30))
    momentum_error = max(
        *(abs(value-target) for value, target in zip(pre, stored_pre)),
        *(abs(value-target) for value, target in zip(post, stored_post)),
    )
    base_control = bool(
        action_error < arb.mpf("5e-8")
        and internal_maximum < arb.mpf("1e-7")
        and momentum_error < arb.mpf("1e-20")
        and base_branch["negative_counts"] == Counter({1: 2400})
        and base_branch["minimum_leading_minor"] > 0
        and base_branch["minimum_argument"] > arb.mpf("1e-6")
        and maximum_imaginary < arb.mpf("1e-70")
    )
    check(
        f"{parity}: the independent base action, equations and momenta reproduce",
        base_control,
        f"action={float(action_error):.3e}, internal={float(internal_maximum):.3e}, "
        f"momentum={float(momentum_error):.3e}",
    )

    tasks = []
    for coordinate in range(95):
        for pair in DERIVATIVE_STEPS.values():
            for h in pair:
                for sign in (1, -1):
                    tasks.append((coordinate, arb.nstr(sign*h, 20)))
    print(f"[{parity}] evaluating {len(tasks)} calibrated gradient points", flush=True)
    with fork_context.Pool(
        processes=8,
        initializer=initialize_worker,
        initargs=(model,),
    ) as pool:
        raw_results = pool.map(gradient_worker, tasks, chunksize=1)

    gradient_values = {}
    branch_pass = True
    minimum_leading_minor = base_branch["minimum_leading_minor"]
    minimum_argument = base_branch["minimum_argument"]
    representative_simplices = base_branch["representative_simplices"]
    for task, raw in zip(tasks, raw_results):
        action = unpack_complex(raw["action"])
        gradient = tuple(unpack_complex(value) for value in raw["gradient"])
        branch = unpack_branch(raw["branch"])
        gradient_values[task] = gradient
        maximum_imaginary = max(
            maximum_imaginary,
            abs(arb.im(action)),
            *(abs(arb.im(value)) for value in gradient),
        )
        minimum_leading_minor = min(
            minimum_leading_minor, branch["minimum_leading_minor"]
        )
        minimum_argument = min(minimum_argument, branch["minimum_argument"])
        representative_simplices += branch["representative_simplices"]
        branch_pass &= (
            branch["negative_counts"] == Counter({1: 2400})
            and branch["minimum_leading_minor"] > 0
            and branch["minimum_argument"] > arb.mpf("1e-6")
        )
    branch_pass = bool(branch_pass and maximum_imaginary < arb.mpf("1e-70"))
    check(
        f"{parity}: all 761 gradients stay on the certified Lorentzian branch",
        branch_pass,
        f"representatives={representative_simplices}, "
        f"minor={arb.nstr(minimum_leading_minor, 8)}, "
        f"argument={arb.nstr(minimum_argument, 8)}, "
        f"imag={arb.nstr(maximum_imaginary, 8)}",
    )

    matrices = {
        name: {
            "primary": matrix_from_gradients(gradient_values, pair[0]),
            "shadow": matrix_from_gradients(gradient_values, pair[1]),
        }
        for name, pair in DERIVATIVE_STEPS.items()
    }
    k_op = matrices["operational"]["primary"]
    k_op_shadow = matrices["operational"]["shadow"]
    k_val = matrices["validation"]["primary"]
    k_val_shadow = matrices["validation"]["shadow"]
    d_k_op = k_op-k_op_shadow
    d_k_val = k_val-k_val_shadow
    d_k_cross = k_op-k_val
    entry_pass = bool(all(
        abs(d_k_cross[row, column])
        <= ENTRY_GATE_FACTOR*(
            abs(d_k_op[row, column])+abs(d_k_val[row, column])+ARITHMETIC_FLOOR
        )
        for row in range(95) for column in range(95)
    ))
    check(
        f"{parity}: all 9025 Hessian entries pass operational/validation calibration",
        entry_pass,
        f"max cross={arb.nstr(max_abs_entry(d_k_cross), 8)}, "
        f"max op proxy={arb.nstr(max_abs_entry(d_k_op), 8)}, "
        f"max val proxy={arb.nstr(max_abs_entry(d_k_val), 8)}",
    )

    # The 95-dimensional calibration norms are used only for reciprocity.
    k_error_spectral = (
        spectral_norm(d_k_op)+spectral_norm(d_k_val)+spectral_norm(d_k_cross)
        + ARITHMETIC_FLOOR
    )
    antisymmetric = k_op-k_op.T
    antisymmetric_spectral = spectral_norm(antisymmetric)
    antisymmetric_frobenius = frobenius_norm(antisymmetric)
    reciprocity_pass = bool(
        antisymmetric_spectral <= ENTRY_GATE_FACTOR*k_error_spectral
    )
    check(
        f"{parity}: the independently differentiated 95-variable Hessian is reciprocal",
        reciprocity_pass,
        f"antisym2={arb.nstr(antisymmetric_spectral, 8)}, "
        f"error2={arb.nstr(k_error_spectral, 8)}",
    )

    k_np = matrix_to_numpy(k_op)
    stored = gauge_input["parities"][parity]["correction"]
    stored_h = np.array(stored["corrected_hessian"], dtype=float)
    stored_b = np.array(stored["corrected_boundary_block"], dtype=float)
    new_h = k_np[30:65, 30:65]
    new_b = k_np[30:65, 65:95]
    h_error = relative_frobenius(new_h, stored_h)
    b_error = relative_frobenius(new_b, stored_b)
    upstream_block_pass = bool(h_error < 1e-6 and b_error < 1e-6)
    check(
        f"{parity}: new internal/final blocks reproduce the corrected upstream blocks",
        upstream_block_pass,
        f"H={h_error:.3e}, B={b_error:.3e}",
    )

    j_op = extract_canonical(k_op)
    j_op_shadow = extract_canonical(k_op_shadow)
    j_val = extract_canonical(k_val)
    j_val_shadow = extract_canonical(k_val_shadow)
    d_op = j_op-j_op_shadow
    d_val = j_val-j_val_shadow
    d_cross = j_op-j_val
    epsilon_global = (
        spectral_norm(d_op)+spectral_norm(d_val)+spectral_norm(d_cross)
        + ARITHMETIC_FLOOR
    )

    print(f"[{parity}] computing the 100-decimal canonical SVD", flush=True)
    u_op, singular_op, v_op = arb.svd_r(j_op)
    u_val, singular_val, v_val = arb.svd_r(j_val)
    singular_values = [singular_op[index] for index in range(65)]
    directional_errors = []
    for index in range(65):
        u = u_op[:, index]
        v = v_op[index, :].T
        directional_errors.append(
            mp_bilinear(u, d_op, v)
            + mp_bilinear(u, d_val, v)
            + mp_bilinear(u, d_cross, v)
            + ARITHMETIC_FLOOR
        )
    resolved = [value > NONZERO_FACTOR*epsilon_global for value in singular_values]
    zero_consistent = [value < ZERO_FACTOR*epsilon_global for value in singular_values]
    open_flags = [not a and not b for a, b in zip(resolved, zero_consistent)]
    resolved_rank = sum(resolved)
    nullity = sum(zero_consistent)

    j_np = matrix_to_numpy(j_op)
    np_u, np_s, np_vh = np.linalg.svd(j_np)
    binary_spectrum_error = max(
        abs(arb.mpf(str(np_s[index]))-singular_values[index])
        / max(arb.mpf(1), singular_values[index])
        for index in range(65)
    )
    relative_ranks = {
        f"{threshold:.0e}": int(np.sum(np_s > threshold*np_s[0]))
        for threshold in (1e-7, 1e-9, 1e-11, 1e-13, 1e-15)
    }

    # Principal angles between the weakest five left/right subspaces.
    op_right_weak = np.array([
        [float(v_op[index, column]) for column in range(65)]
        for index in range(60, 65)
    ])
    val_right_weak = np.array([
        [float(v_val[index, column]) for column in range(65)]
        for index in range(60, 65)
    ])
    op_left_weak = np.array([
        [float(u_op[row, index]) for index in range(60, 65)]
        for row in range(65)
    ])
    val_left_weak = np.array([
        [float(u_val[row, index]) for index in range(60, 65)]
        for row in range(65)
    ])
    right_cosines = np.linalg.svd(
        op_right_weak@val_right_weak.T, compute_uv=False
    )
    left_cosines = np.linalg.svd(
        op_left_weak.T@val_left_weak, compute_uv=False
    )
    right_angles = np.arccos(np.clip(right_cosines, -1, 1))
    left_angles = np.arccos(np.clip(left_cosines, -1, 1))

    transformed_block_norms = block_norms(j_np)
    vector_diagnostics = []
    diagnostic_indices = sorted(set(
        list(range(60, 65))
        + [index for index, flag in enumerate(open_flags) if flag]
        + [index for index, flag in enumerate(zero_consistent) if flag]
    ))
    for index in diagnostic_indices:
        right = np.array([float(v_op[index, column]) for column in range(65)])
        left = np.array([float(u_op[row, index]) for row in range(65)])
        vector_diagnostics.append({
            "index_descending": index,
            "singular_value": arb.nstr(singular_values[index], 70),
            "left_sector_overlaps": sector_overlaps(left),
            "right_sector_overlaps": sector_overlaps(right),
        })

    null_projection_singular = []
    null_projection_rank = 0
    null_projection_tolerance = None
    if nullity:
        null_indices = [
            index for index, flag in enumerate(zero_consistent) if flag
        ]
        projected = arb.matrix(nullity, 30)
        for row, index in enumerate(null_indices):
            for column in range(30):
                projected[row, column] = v_op[index, 35+column]
        projected_svd = arb.svd_r(projected, compute_uv=False)
        null_projection_singular = [
            projected_svd[index] for index in range(projected_svd.rows)
        ]
        resolved_values = [
            value for value, flag in zip(singular_values, resolved) if flag
        ]
        if resolved_values:
            null_projection_tolerance = min(
                arb.mpf(1),
                NONZERO_FACTOR*epsilon_global/min(resolved_values),
            )
            null_projection_rank = sum(
                value > null_projection_tolerance
                for value in null_projection_singular
            )
        else:
            null_projection_rank = None

    controls_pass = bool(
        base_control and branch_pass and entry_pass and reciprocity_pass
        and upstream_block_pass
    )
    if not controls_pass:
        outcome = "CANONICAL_GRADIENT_CONTROL_FAILED"
    elif any(open_flags):
        outcome = "CANONICAL_RANK_NUMERICALLY_OPEN"
    elif resolved_rank == 65:
        outcome = "CANONICAL_LEGENDRE_REGULAR"
    elif nullity == 1:
        null_index = zero_consistent.index(True)
        null_right = np.array([
            float(v_op[null_index, column]) for column in range(65)
        ])
        null_overlap = sector_overlaps(null_right)
        if (
            null_projection_rank == 0
            and null_overlap["internal_lapse"] > 0.999999
        ):
            outcome = "ONE_CANONICAL_LAPSE_NULL"
        else:
            outcome = "ADDITIONAL_CANONICAL_DEGENERACY"
    else:
        outcome = "ADDITIONAL_CANONICAL_DEGENERACY"

    resolved_pseudoconstraints = [
        index for index, flag in enumerate(resolved)
        if flag and singular_values[index] < arb.mpf("1e-6")*singular_values[0]
    ]
    check(
        f"{parity}: the canonical rank outcome is assigned by calibrated bands",
        outcome in {
            "CANONICAL_GRADIENT_CONTROL_FAILED",
            "CANONICAL_RANK_NUMERICALLY_OPEN",
            "CANONICAL_LEGENDRE_REGULAR",
            "ONE_CANONICAL_LAPSE_NULL",
            "ADDITIONAL_CANONICAL_DEGENERACY",
        },
        f"outcome={outcome}, resolved={resolved_rank}, null={nullity}, "
        f"open={sum(open_flags)}, smin/error={arb.nstr(singular_values[-1]/epsilon_global, 8)}",
    )

    records[parity] = {
        "outcome": outcome,
        "controls_pass": controls_pass,
        "base": {
            "action": complex_record(base_action),
            "relative_action_error": float(action_error),
            "internal_residuals": [complex_record(value) for value in internal_residuals],
            "maximum_internal_residual": arb.nstr(internal_maximum, 30),
            "pre_momenta": [complex_record(value) for value in pre],
            "post_momenta": [complex_record(value) for value in post],
            "maximum_momentum_control_error": arb.nstr(momentum_error, 30),
        },
        "branch": {
            "precision_digits": DPS,
            "gradient_evaluations": 761,
            "representative_simplices": representative_simplices,
            "minimum_absolute_leading_principal_minor": arb.nstr(
                minimum_leading_minor, 30
            ),
            "minimum_angle_argument_modulus": arb.nstr(minimum_argument, 30),
            "maximum_imaginary_contamination": arb.nstr(maximum_imaginary, 30),
            "pass": branch_pass,
        },
        "calibration": {
            "steps": {
                name: [arb.nstr(value, 10) for value in pair]
                for name, pair in DERIVATIVE_STEPS.items()
            },
            "entrywise_pass": entry_pass,
            "maximum_operational_proxy": arb.nstr(max_abs_entry(d_k_op), 30),
            "maximum_validation_proxy": arb.nstr(max_abs_entry(d_k_val), 30),
            "maximum_cross_difference": arb.nstr(max_abs_entry(d_k_cross), 30),
            "hessian_error_spectral": arb.nstr(k_error_spectral, 30),
            "hessian_antisymmetry_spectral": arb.nstr(antisymmetric_spectral, 30),
            "hessian_antisymmetry_frobenius": arb.nstr(antisymmetric_frobenius, 30),
            "reciprocity_pass": reciprocity_pass,
            "upstream_internal_relative_frobenius": float(h_error),
            "upstream_boundary_relative_frobenius": float(b_error),
            "upstream_blocks_pass": upstream_block_pass,
        },
        "canonical_matrix": [
            [arb.nstr(j_op[row, column], 70) for column in range(65)]
            for row in range(65)
        ],
        "spectrum": {
            "singular_values": [arb.nstr(value, 70) for value in singular_values],
            "epsilon_global": arb.nstr(epsilon_global, 70),
            "directional_errors": [
                arb.nstr(value, 70) for value in directional_errors
            ],
            "ratios_to_global_error": [
                arb.nstr(value/epsilon_global, 30) for value in singular_values
            ],
            "ratios_to_directional_error": [
                arb.nstr(value/error, 30)
                for value, error in zip(singular_values, directional_errors)
            ],
            "resolved_nonzero": resolved,
            "zero_consistent": zero_consistent,
            "numerically_open": open_flags,
            "resolved_rank": resolved_rank,
            "nullity": nullity,
            "relative_threshold_ranks": relative_ranks,
            "condition_2": (
                arb.nstr(singular_values[0]/singular_values[-1], 30)
                if singular_values[-1] else None
            ),
            "binary64_normalized_spectrum_error": float(binary_spectrum_error),
            "weakest_five_right_principal_angles": right_angles.tolist(),
            "weakest_five_left_principal_angles": left_angles.tolist(),
            "resolved_pseudoconstraint_indices_at_relative_1e-6": (
                resolved_pseudoconstraints
            ),
        },
        "sector_decomposition": {
            "row_column_order": ["lapse", "transverse", "scale", "shape"],
            "block_frobenius_norms": transformed_block_norms,
            "weak_vector_diagnostics": vector_diagnostics,
            "null_final_boundary_projection_singular_values": [
                arb.nstr(value, 70) for value in null_projection_singular
            ],
            "null_final_boundary_projection_tolerance": (
                arb.nstr(null_projection_tolerance, 70)
                if null_projection_tolerance is not None else None
            ),
            "null_final_boundary_projection_rank": null_projection_rank,
        },
    }


outcomes = {record["outcome"] for record in records.values()}
continuation_accepted = outcomes <= {
    "CANONICAL_LEGENDRE_REGULAR", "ONE_CANONICAL_LAPSE_NULL"
}
if "CANONICAL_GRADIENT_CONTROL_FAILED" in outcomes:
    verdict = (
        "CONTROL FAILURE: at least one canonical gradient certificate failed; "
        "no rank or nonlinear-evolution claim follows."
    )
elif outcomes == {"CANONICAL_LEGENDRE_REGULAR"}:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL: the pre-Legendre inversion is regular "
        "in both order-24 schedule carriers.  Fixed pre-momentum lifts the "
        "Dirichlet lapse degeneracy.  No nonlinear next frame is solved."
    )
elif outcomes == {"ONE_CANONICAL_LAPSE_NULL"}:
    verdict = (
        "DERIVED COMPUTATIONAL LOCAL: both pre-Legendre maps retain exactly "
        "the preregistered internal collective lapse null; geometry evolution "
        "is locally unique modulo that generator.  No nonlinear frame is solved."
    )
else:
    verdict = (
        "DERIVED/OPEN CANONICAL CENSUS: the two parity outcomes are "
        f"{sorted(outcomes)}; the nonlinear next-frame solve remains blocked."
    )

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "tolerance_commit": TOLERANCE_COMMIT,
    "clarification_commit": CLARIFICATION_COMMIT,
    "projection_rule_commit": PROJECTION_RULE_COMMIT,
    "gluing_result_commit": GLUING_RESULT_COMMIT,
    "prediction_frozen_before_evaluation": {
        "rank": 65,
        "weakest_sector": "collective-lapse / homogeneous-scale",
    },
    "coordinates": {
        "columns": "log(x[35]), log(q_new[30])",
        "rows": "per-edge log internal equations[35], pre-momenta[30]",
        "complete_hessian_order": "log(q_old[30]), log(x[35]), log(q_new[30])",
    },
    "parities": records,
    "outcomes": sorted(outcomes),
    "nonlinear_continuation_accepted": continuation_accepted,
    "labels": {
        "rank_census": (
            "CONTROL FAILURE"
            if "CANONICAL_GRADIENT_CONTROL_FAILED" in outcomes
            else "OPEN NUMERIC"
            if "CANONICAL_RANK_NUMERICALLY_OPEN" in outcomes
            else "DERIVED COMPUTATIONAL LOCAL"
        ),
        "nonlinear_next_frame": "NOT SOLVED",
        "continuum_degrees_of_freedom": "NOT ESTABLISHED",
        "inflation_or_physical_clock": "NOT ESTABLISHED",
    },
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOMES: {sorted(outcomes)}")
print(verdict)
print(f"NONLINEAR CONTINUATION ACCEPTED: {continuation_accepted}")
raise SystemExit(0 if passed == tests else 1)
