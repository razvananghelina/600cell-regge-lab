#!/usr/bin/env python3
"""Blind canonical tangent census about the first dynamic 600-cell dust tick.

Prior-art commit: 25722d9.
Protocol commit: 0bceb9b.
No continuum harmonic or speed target is loaded.
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
import multiprocessing as mp_pool
from pathlib import Path
import sys

import mpmath as arb
import numpy as np
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_dust_dynamic_tangent.json"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
GLUING_INPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
CANONICAL_SOURCE = HERE / "verify_gravity_600cell_dust_canonical_legendre_rank.py"
ACTION_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"

PRIOR_ART_COMMIT = "25722d9"
PROTOCOL_COMMIT = "0bceb9b"
INPUT_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "gluing": "a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77",
    "canonical_source": "396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e",
    "action_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}

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
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_audited_functions():
    wanted = {
        "orbit_sort_key",
        "augment_boundary_orbits",
        "arb_log_minus",
        "arb_signed_volume_square",
        "arb_angle_record",
        "triangle_area_square",
        "triangle_area_square_partials",
        "edge_data",
        "simplex_squared",
        "action_and_gradient",
        "perturb_base",
        "pack_complex",
        "unpack_complex",
        "pack_branch",
        "unpack_branch",
        "initialize_worker",
        "gradient_worker",
        "matrix_from_gradients",
        "max_abs_entry",
        "frobenius_norm",
        "spectral_norm",
        "extract_canonical",
        "matrix_to_numpy",
    }
    tree = ast.parse(CANONICAL_SOURCE.read_text(), filename=str(CANONICAL_SOURCE))
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch: missing={wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(CANONICAL_SOURCE), "exec"),
        globals(),
    )


hashes = {
    "tick": file_digest(TICK_INPUT),
    "gluing": file_digest(GLUING_INPUT),
    "canonical_source": file_digest(CANONICAL_SOURCE),
    "action_source": file_digest(ACTION_SOURCE),
}
load_audited_functions()

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_dynamic_tangent", ACTION_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise

print("Blind dynamic canonical-tangent census", flush=True)
check(
    "the imported one-slab action core retains all 43 certificates",
    gro.tests == gro.passed == 43,
)

models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}
carrier_ok = all(
    len(model["old_orbits"]) == 30
    and len(model["edge_orbits"]) == 35
    and len(model["pole_orbits"]) == 5
    and len(model["final_orbits"]) == 30
    and Counter(map(len, model["old_orbits"])) == Counter({24: 30})
    and Counter(map(len, model["edge_orbits"])) == Counter({24: 35})
    and Counter(map(len, model["final_orbits"])) == Counter({24: 30})
    for model in models.values()
)
check("both dynamic tangent carriers have the frozen 30+35+30 orbits", carrier_ok)

tick = json.loads(TICK_INPUT.read_text())
gluing = json.loads(GLUING_INPUT.read_text())
maps = {
    parity: tuple(
        gluing["parities"][parity]["geometry"]["old_to_final_orbit_map"]
    )
    for parity in ("even", "odd")
}
provenance_ok = bool(
    hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "25722d9"
    and PROTOCOL_COMMIT == "0bceb9b"
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tick["passed"] == tick["tests"] == 7
    and set(tick["solutions"]) == {"even", "odd"}
    and gluing["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and all(sorted(mapping) == list(range(30)) for mapping in maps.values())
)
check(
    "all source hashes, dynamic states and gluing maps are frozen",
    provenance_ok,
    str(hashes),
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


def mp_matrix_to_strings(matrix, digits=50):
    return [
        [arb.nstr(matrix[row, column], digits) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def make_omega(size):
    omega = arb.matrix(2*size, 2*size)
    for index in range(size):
        omega[index, size+index] = 1
        omega[size+index, index] = -1
    return omega


OMEGA = make_omega(30)


def build_tangent(matrix, mapping):
    j_matrix = extract_canonical(matrix)
    rhs = arb.matrix(65, 60)
    for row in range(35):
        for column in range(30):
            rhs[row, column] = -matrix[30+row, column]
    for row in range(30):
        for column in range(30):
            rhs[35+row, column] = matrix[row, column]
        rhs[35+row, 30+row] = 1
    solved = (j_matrix**-1)*rhs
    raw = arb.matrix(60, 60)
    for row in range(30):
        for column in range(60):
            raw[row, column] = solved[35+row, column]
            direct = matrix[65+row, column] if column < 30 else arb.mpf(0)
            induced = sum(
                matrix[65+row, 30+internal]*solved[internal, column]
                for internal in range(35)
            ) + sum(
                matrix[65+row, 65+new]*solved[35+new, column]
                for new in range(30)
            )
            raw[30+row, column] = direct+induced
    tangent = arb.matrix(60, 60)
    for next_index, final_index in enumerate(mapping):
        for column in range(60):
            tangent[next_index, column] = raw[final_index, column]
            tangent[30+next_index, column] = raw[30+final_index, column]
    return j_matrix, tangent


def symplectic_defect(matrix):
    return matrix.T*OMEGA*matrix-OMEGA


def singular_values_mp(matrix):
    values = arb.svd_r(matrix, compute_uv=False)
    return tuple(values[index] for index in range(values.rows))


def reciprocal_vector(values):
    size = len(values)
    return tuple(values[index]*values[size-1-index]-1 for index in range(size//2))


def infinity_vector(values):
    return max(abs(value) for value in values) if values else arb.mpf(0)


def vector_difference(left, right):
    return tuple(a-b for a, b in zip(left, right))


def householder_boundary_basis():
    scale = arb.matrix([1/arb.sqrt(30) for _ in range(30)])
    vector = scale.copy()
    vector[0] += 1
    denominator = (vector.T*vector)[0]
    reflector = arb.eye(30)-2*(vector*vector.T)/denominator
    basis = arb.matrix(30, 30)
    for row in range(30):
        basis[row, 0] = scale[row]
        for column in range(1, 30):
            basis[row, column] = reflector[row, column]
    return basis


BOUNDARY_BASIS = householder_boundary_basis()
PHASE_BASIS = arb.matrix(60, 60)
for block in (0, 30):
    for row in range(30):
        for column in range(30):
            PHASE_BASIS[block+row, block+column] = BOUNDARY_BASIS[row, column]


def sector_matrix(matrix):
    return PHASE_BASIS.T*matrix*PHASE_BASIS


SCALE_INDICES = (0, 30)
SHAPE_INDICES = tuple(list(range(1, 30))+list(range(31, 60)))


def select_matrix(matrix, rows, columns):
    return arb.matrix([[matrix[row, column] for column in columns] for row in rows])


def mixing_vector(matrix):
    transformed = sector_matrix(matrix)
    values = []
    for row in SCALE_INDICES:
        for column in SHAPE_INDICES:
            values.append(transformed[row, column])
    for row in SHAPE_INDICES:
        for column in SCALE_INDICES:
            values.append(transformed[row, column])
    return tuple(values)


def vector_frobenius(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def optimal_eigen_match(left, right):
    left = np.asarray(left, dtype=np.complex128)
    right = np.asarray(right, dtype=np.complex128)
    rows, columns = linear_sum_assignment(np.abs(left[:, None]-right[None, :]))
    distances = np.abs(left[rows]-right[columns])
    reordered = np.empty_like(left)
    reordered[rows] = right[columns]
    return float(np.max(distances)), reordered


def numerical_spectrum(matrix):
    array = matrix_to_numpy(matrix)
    singular = np.linalg.svd(array, compute_uv=False)
    eigenvalues, eigenvectors = np.linalg.eig(array)
    eigenvector_condition = float(np.linalg.cond(eigenvectors))
    sign, logabsdet = np.linalg.slogdet(array)
    determinant = float(sign*np.exp(logabsdet)) if logabsdet < 700 else math.copysign(math.inf, sign)
    return {
        "array": array,
        "singular": singular,
        "eigenvalues": eigenvalues,
        "eigenvector_condition": eigenvector_condition,
        "determinant": determinant,
    }


def calibrated_spectrum(matrices):
    spectra = {name: numerical_spectrum(matrix) for name, matrix in matrices.items()}
    operational = spectra["operational"]
    scale_eig = max(1.0, float(np.max(np.abs(operational["eigenvalues"]))))
    eigen_errors = []
    for name in ("operational_shadow", "validation", "validation_shadow"):
        distance, _ = optimal_eigen_match(
            operational["eigenvalues"], spectra[name]["eigenvalues"]
        )
        eigen_errors.append(distance)
    epsilon_eig = max(eigen_errors)+1e-15*scale_eig
    scale_svd = max(1.0, float(operational["singular"][0]))
    epsilon_svd = max(
        float(np.max(np.abs(operational["singular"]-spectra[name]["singular"])))
        for name in ("operational_shadow", "validation", "validation_shadow")
    )+1e-15*scale_svd
    moduli = np.abs(operational["eigenvalues"])
    distances_to_unit = np.abs(moduli-1)
    unit_consistent = int(np.sum(distances_to_unit <= 10*epsilon_eig))
    off_unit = int(np.sum(distances_to_unit > 100*epsilon_eig))
    open_count = len(moduli)-unit_consistent-off_unit
    return {
        "raw": spectra,
        "epsilon_eig": epsilon_eig,
        "epsilon_svd": epsilon_svd,
        "unit_consistent": unit_consistent,
        "off_unit": off_unit,
        "open": open_count,
    }


fork_context = mp_pool.get_context("fork")
records = {}
all_base_ok = True
all_branch_ok = True
all_entry_ok = True
all_reciprocity_ok = True

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing the accepted dynamic background", flush=True)
    state = tuple(arb.mpf(value) for value in tick["solutions"][parity]["state"])
    scale_log, rho_log = state
    rho = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(scale_log)*ARB_L0_SQUARE-rho
    ARB_BASE_OLD = tuple(ARB_L0_SQUARE for _ in range(30))
    ARB_BASE_X = tuple([diagonal for _ in range(30)]+[rho for _ in range(5)])
    ARB_BASE_NEW = tuple(arb.exp(2*scale_log)*ARB_L0_SQUARE for _ in range(30))

    model = models[parity]
    base_action, base_gradient, base_branch = action_and_gradient(
        model, ARB_BASE_OLD, ARB_BASE_X, ARB_BASE_NEW
    )
    stored_pre = tuple(
        arb.mpf(value) for value in tick["solutions"][parity]["pre_momentum"]
    )
    stored_post = tuple(
        arb.mpf(value) for value in tick["solutions"][parity]["post_momentum"]
    )
    pre = tuple(-arb.re(base_gradient[index]) for index in range(30))
    post = tuple(arb.re(base_gradient[65+index]) for index in range(30))
    internal_maximum = max(abs(value) for value in base_gradient[30:65])
    momentum_error = max(
        *(abs(value-target) for value, target in zip(pre, stored_pre)),
        *(abs(value-target) for value, target in zip(post, stored_post)),
    )
    maximum_imaginary = max(
        abs(arb.im(base_action)), *(abs(arb.im(value)) for value in base_gradient)
    )
    base_ok = bool(
        internal_maximum < arb.mpf("1e-25")
        and momentum_error < arb.mpf("1e-45")
        and base_branch["negative_counts"] == Counter({1: 2400})
        and base_branch["minimum_leading_minor"] > 0
        and base_branch["minimum_argument"] > arb.mpf("1e-6")
        and maximum_imaginary < arb.mpf("1e-70")
    )
    all_base_ok &= base_ok
    check(
        f"{parity}: the dynamic equations and committed momenta reproduce",
        base_ok,
        f"internal={arb.nstr(internal_maximum, 8)}, momentum={arb.nstr(momentum_error, 8)}",
    )

    tasks = []
    for coordinate in range(95):
        for pair in DERIVATIVE_STEPS.values():
            for step in pair:
                for sign in (1, -1):
                    tasks.append((coordinate, arb.nstr(sign*step, 20)))
    print(f"[{parity}] evaluating {len(tasks)} calibrated dynamic gradients", flush=True)
    with fork_context.Pool(
        processes=8, initializer=initialize_worker, initargs=(model,)
    ) as pool:
        raw_results = pool.map(gradient_worker, tasks, chunksize=1)

    gradient_values = {}
    branch_ok = True
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
        branch_ok &= (
            branch["negative_counts"] == Counter({1: 2400})
            and branch["minimum_leading_minor"] > 0
            and branch["minimum_argument"] > arb.mpf("1e-6")
        )
    branch_ok = bool(branch_ok and maximum_imaginary < arb.mpf("1e-70"))
    all_branch_ok &= branch_ok
    check(
        f"{parity}: all 761 dynamic gradients retain the Lorentzian branch",
        branch_ok,
        f"representatives={representative_simplices}, minor={arb.nstr(minimum_leading_minor, 8)}, "
        f"argument={arb.nstr(minimum_argument, 8)}, imag={arb.nstr(maximum_imaginary, 8)}",
    )

    matrices = {
        "operational": matrix_from_gradients(
            gradient_values, DERIVATIVE_STEPS["operational"][0]
        ),
        "operational_shadow": matrix_from_gradients(
            gradient_values, DERIVATIVE_STEPS["operational"][1]
        ),
        "validation": matrix_from_gradients(
            gradient_values, DERIVATIVE_STEPS["validation"][0]
        ),
        "validation_shadow": matrix_from_gradients(
            gradient_values, DERIVATIVE_STEPS["validation"][1]
        ),
    }
    d_op = matrices["operational"]-matrices["operational_shadow"]
    d_val = matrices["validation"]-matrices["validation_shadow"]
    d_cross = matrices["operational"]-matrices["validation"]
    entry_ok = bool(all(
        abs(d_cross[row, column])
        <= ENTRY_GATE_FACTOR*(
            abs(d_op[row, column])+abs(d_val[row, column])+ARITHMETIC_FLOOR
        )
        for row in range(95) for column in range(95)
    ))
    all_entry_ok &= entry_ok
    check(
        f"{parity}: all 9025 dynamic Hessian entries pass calibration",
        entry_ok,
        f"cross={arb.nstr(max_abs_entry(d_cross), 8)}, "
        f"op={arb.nstr(max_abs_entry(d_op), 8)}, val={arb.nstr(max_abs_entry(d_val), 8)}",
    )

    hessian_error = (
        spectral_norm(d_op)+spectral_norm(d_val)+spectral_norm(d_cross)
        + ARITHMETIC_FLOOR
    )
    antisymmetric = matrices["operational"]-matrices["operational"].T
    antisymmetric_norm = spectral_norm(antisymmetric)
    reciprocity_ok = bool(antisymmetric_norm <= 10*hessian_error)
    all_reciprocity_ok &= reciprocity_ok
    check(
        f"{parity}: the dynamic complete Hessian is reciprocal",
        reciprocity_ok,
        f"antisym={arb.nstr(antisymmetric_norm, 8)}, error={arb.nstr(hessian_error, 8)}",
    )

    canonical = {}
    tangents = {}
    for name, matrix in matrices.items():
        canonical[name], tangents[name] = build_tangent(matrix, maps[parity])
    j_op = canonical["operational"]
    epsilon_j = (
        spectral_norm(j_op-canonical["operational_shadow"])
        + spectral_norm(canonical["validation"]-canonical["validation_shadow"])
        + spectral_norm(j_op-canonical["validation"])
        + ARITHMETIC_FLOOR
    )
    j_singular = singular_values_mp(j_op)
    j_resolved = tuple(value > 100*epsilon_j for value in j_singular)
    j_zero = tuple(value < 10*epsilon_j for value in j_singular)
    j_open = tuple(not nonzero and not zero for nonzero, zero in zip(j_resolved, j_zero))

    epsilon_t = (
        spectral_norm(tangents["operational"]-tangents["operational_shadow"])
        + spectral_norm(tangents["validation"]-tangents["validation_shadow"])
        + spectral_norm(tangents["operational"]-tangents["validation"])
        + ARITHMETIC_FLOOR
    )
    defects = {name: symplectic_defect(value) for name, value in tangents.items()}
    epsilon_sym = (
        spectral_norm(defects["operational"]-defects["operational_shadow"])
        + spectral_norm(defects["validation"]-defects["validation_shadow"])
        + spectral_norm(defects["operational"]-defects["validation"])
        + ARITHMETIC_FLOOR
    )
    symplectic_norm = spectral_norm(defects["operational"])
    symplectic_ok = bool(symplectic_norm <= 10*epsilon_sym)

    tangent_singular = {
        name: singular_values_mp(value) for name, value in tangents.items()
    }
    reciprocal = {
        name: reciprocal_vector(value) for name, value in tangent_singular.items()
    }
    epsilon_reciprocal = (
        infinity_vector(vector_difference(
            reciprocal["operational"], reciprocal["operational_shadow"]
        ))
        + infinity_vector(vector_difference(
            reciprocal["validation"], reciprocal["validation_shadow"]
        ))
        + infinity_vector(vector_difference(
            reciprocal["operational"], reciprocal["validation"]
        ))
        + ARITHMETIC_FLOOR
    )
    reciprocal_norm = infinity_vector(reciprocal["operational"])
    reciprocal_ok = bool(reciprocal_norm <= 10*epsilon_reciprocal)

    mixing = {name: mixing_vector(value) for name, value in tangents.items()}
    epsilon_mix = (
        vector_frobenius(vector_difference(
            mixing["operational"], mixing["operational_shadow"]
        ))
        + vector_frobenius(vector_difference(
            mixing["validation"], mixing["validation_shadow"]
        ))
        + vector_frobenius(vector_difference(
            mixing["operational"], mixing["validation"]
        ))
        + ARITHMETIC_FLOOR
    )
    mixing_norm = vector_frobenius(mixing["operational"])
    if mixing_norm <= 10*epsilon_mix:
        mixing_label = "SCALE_SHAPE_INVARIANT"
    elif mixing_norm > 100*epsilon_mix:
        mixing_label = "SCALE_SHAPE_MIXED"
    else:
        mixing_label = "SCALE_SHAPE_OPEN"

    full_spectrum = calibrated_spectrum(tangents)
    shape_matrices = None
    shape_spectrum = None
    scale_spectrum = None
    if mixing_label == "SCALE_SHAPE_INVARIANT":
        sector_matrices = {name: sector_matrix(value) for name, value in tangents.items()}
        shape_matrices = {
            name: select_matrix(value, SHAPE_INDICES, SHAPE_INDICES)
            for name, value in sector_matrices.items()
        }
        scale_matrices = {
            name: select_matrix(value, SCALE_INDICES, SCALE_INDICES)
            for name, value in sector_matrices.items()
        }
        shape_spectrum = calibrated_spectrum(shape_matrices)
        scale_spectrum = calibrated_spectrum(scale_matrices)

    print(
        f"[{parity}] J smin={arb.nstr(j_singular[-1], 10)} "
        f"epsJ={arb.nstr(epsilon_j, 8)} sym={arb.nstr(symplectic_norm, 8)} "
        f"epsSym={arb.nstr(epsilon_sym, 8)} mix={arb.nstr(mixing_norm, 8)} "
        f"epsMix={arb.nstr(epsilon_mix, 8)} {mixing_label}",
        flush=True,
    )

    records[parity] = {
        "base": {
            "state": state,
            "internal_maximum": internal_maximum,
            "momentum_error": momentum_error,
            "maximum_imaginary": maximum_imaginary,
            "minimum_leading_minor": minimum_leading_minor,
            "minimum_argument": minimum_argument,
            "passed": base_ok,
        },
        "branch_ok": branch_ok,
        "entry_ok": entry_ok,
        "reciprocity_ok": reciprocity_ok,
        "hessian_error": hessian_error,
        "antisymmetric_norm": antisymmetric_norm,
        "canonical": canonical,
        "epsilon_j": epsilon_j,
        "j_singular": j_singular,
        "j_resolved": j_resolved,
        "j_zero": j_zero,
        "j_open": j_open,
        "tangents": tangents,
        "epsilon_t": epsilon_t,
        "symplectic_norm": symplectic_norm,
        "epsilon_sym": epsilon_sym,
        "symplectic_ok": symplectic_ok,
        "tangent_singular": tangent_singular,
        "reciprocal_norm": reciprocal_norm,
        "epsilon_reciprocal": epsilon_reciprocal,
        "reciprocal_ok": reciprocal_ok,
        "mixing_norm": mixing_norm,
        "epsilon_mix": epsilon_mix,
        "mixing_label": mixing_label,
        "full_spectrum": full_spectrum,
        "shape_matrices": shape_matrices,
        "shape_spectrum": shape_spectrum,
        "scale_spectrum": scale_spectrum,
    }


all_j_resolved = all(all(record["j_resolved"]) for record in records.values())
any_j_open = any(any(record["j_open"]) for record in records.values())
any_j_zero = any(any(record["j_zero"]) for record in records.values())
all_symplectic = all(
    record["symplectic_ok"] and record["reciprocal_ok"]
    for record in records.values()
)
mixing_labels = {record["mixing_label"] for record in records.values()}

schedule_label = "NOT_AVAILABLE"
schedule_record = None
if mixing_labels == {"SCALE_SHAPE_INVARIANT"}:
    even = records["even"]["shape_spectrum"]
    odd = records["odd"]["shape_spectrum"]
    eigen_distance, _ = optimal_eigen_match(
        even["raw"]["operational"]["eigenvalues"],
        odd["raw"]["operational"]["eigenvalues"],
    )
    singular_distance = float(np.max(np.abs(
        even["raw"]["operational"]["singular"]
        - odd["raw"]["operational"]["singular"]
    )))
    eigen_uncertainty = even["epsilon_eig"]+odd["epsilon_eig"]
    singular_uncertainty = even["epsilon_svd"]+odd["epsilon_svd"]
    if (
        eigen_distance <= 10*eigen_uncertainty
        and singular_distance <= 10*singular_uncertainty
    ):
        schedule_label = "SCHEDULE_ROBUST"
    elif (
        eigen_distance > 100*eigen_uncertainty
        or singular_distance > 100*singular_uncertainty
    ):
        schedule_label = "SCHEDULE_DEPENDENT"
    else:
        schedule_label = "SCHEDULE_OPEN"
    schedule_record = {
        "eigen_distance": eigen_distance,
        "eigen_uncertainty": eigen_uncertainty,
        "singular_distance": singular_distance,
        "singular_uncertainty": singular_uncertainty,
        "label": schedule_label,
    }

controls_ok = bool(
    provenance_ok and carrier_ok and all_base_ok and all_branch_ok
    and all_entry_ok and all_reciprocity_ok
)
if not controls_ok:
    outcome = "DYNAMIC_TANGENT_CONTROL_FAILED"
elif any_j_open:
    outcome = "DYNAMIC_CANONICAL_RANK_OPEN"
elif any_j_zero or not all_j_resolved:
    outcome = "DYNAMIC_CANONICAL_DEGENERATE"
elif not all_symplectic:
    outcome = "DYNAMIC_TANGENT_SYMPLECTICITY_FAILED"
elif "SCALE_SHAPE_MIXED" in mixing_labels:
    outcome = "DYNAMIC_SCALE_SHAPE_MIXED"
elif "SCALE_SHAPE_OPEN" in mixing_labels:
    outcome = "DYNAMIC_SCALE_SHAPE_OPEN"
elif schedule_label == "SCHEDULE_DEPENDENT":
    outcome = "DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT"
elif schedule_label == "SCHEDULE_OPEN":
    outcome = "DYNAMIC_SHAPE_TANGENT_SCHEDULE_OPEN"
else:
    outcome = "DYNAMIC_SHAPE_TANGENT_BLIND_CENSUS_CERTIFIED"

allowed_outcomes = {
    "DYNAMIC_TANGENT_CONTROL_FAILED",
    "DYNAMIC_CANONICAL_RANK_OPEN",
    "DYNAMIC_CANONICAL_DEGENERATE",
    "DYNAMIC_TANGENT_SYMPLECTICITY_FAILED",
    "DYNAMIC_SCALE_SHAPE_MIXED",
    "DYNAMIC_SCALE_SHAPE_OPEN",
    "DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT",
    "DYNAMIC_SHAPE_TANGENT_SCHEDULE_OPEN",
    "DYNAMIC_SHAPE_TANGENT_BLIND_CENSUS_CERTIFIED",
}
check(
    "the frozen hierarchy assigns one target-independent tangent outcome",
    outcome in allowed_outcomes,
    outcome,
)


def serialize_spectrum(record):
    operational = record["raw"]["operational"]
    return {
        "singular_values": [f"{value:.17e}" for value in operational["singular"]],
        "eigenvalues": [
            {"real": f"{value.real:.17e}", "imaginary": f"{value.imag:.17e}"}
            for value in operational["eigenvalues"]
        ],
        "eigenvalue_moduli": [
            f"{abs(value):.17e}" for value in operational["eigenvalues"]
        ],
        "spectral_radius": f"{max(abs(value) for value in operational['eigenvalues']):.17e}",
        "determinant": f"{operational['determinant']:.17e}",
        "eigenvector_condition": f"{operational['eigenvector_condition']:.17e}",
        "epsilon_eig": f"{record['epsilon_eig']:.17e}",
        "epsilon_svd": f"{record['epsilon_svd']:.17e}",
        "unit_consistent_count": record["unit_consistent"],
        "off_unit_count": record["off_unit"],
        "open_unit_count": record["open"],
    }


artifact_records = {}
for parity, record in records.items():
    artifact_records[parity] = {
        "base": {
            key: (
                value if isinstance(value, bool)
                else [arb.nstr(item, 50) for item in value]
                if isinstance(value, tuple)
                else arb.nstr(value, 50)
            )
            for key, value in record["base"].items()
        },
        "branch_ok": record["branch_ok"],
        "entry_ok": record["entry_ok"],
        "reciprocity_ok": record["reciprocity_ok"],
        "hessian_error": arb.nstr(record["hessian_error"], 40),
        "antisymmetric_norm": arb.nstr(record["antisymmetric_norm"], 40),
        "epsilon_j": arb.nstr(record["epsilon_j"], 40),
        "canonical_singular_values": [arb.nstr(value, 50) for value in record["j_singular"]],
        "canonical_resolved_rank": sum(record["j_resolved"]),
        "canonical_zero_count": sum(record["j_zero"]),
        "canonical_open_count": sum(record["j_open"]),
        "epsilon_t": arb.nstr(record["epsilon_t"], 40),
        "symplectic_norm": arb.nstr(record["symplectic_norm"], 40),
        "epsilon_sym": arb.nstr(record["epsilon_sym"], 40),
        "symplectic_ok": record["symplectic_ok"],
        "reciprocal_norm": arb.nstr(record["reciprocal_norm"], 40),
        "epsilon_reciprocal": arb.nstr(record["epsilon_reciprocal"], 40),
        "reciprocal_ok": record["reciprocal_ok"],
        "mixing_norm": arb.nstr(record["mixing_norm"], 40),
        "epsilon_mix": arb.nstr(record["epsilon_mix"], 40),
        "mixing_label": record["mixing_label"],
        "tangent_matrix": mp_matrix_to_strings(record["tangents"]["operational"], 50),
        "full_spectrum": serialize_spectrum(record["full_spectrum"]),
        "shape_spectrum": (
            serialize_spectrum(record["shape_spectrum"])
            if record["shape_spectrum"] is not None else None
        ),
        "scale_spectrum": (
            serialize_spectrum(record["scale_spectrum"])
            if record["scale_spectrum"] is not None else None
        ),
    }

n_distinct = (
    1 if schedule_label == "SCHEDULE_ROBUST"
    else 2 if schedule_label == "SCHEDULE_DEPENDENT"
    else "OPEN"
)
artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "full_720_edge_carrier": False,
    "number_of_maps": 2,
    "number_of_distinct_maps": n_distinct,
    "parities": artifact_records,
    "schedule": schedule_record,
    "classification": {
        "canonical_rank_resolved": all_j_resolved,
        "symplectic": all_symplectic,
        "scale_shape_labels": sorted(mixing_labels),
        "schedule_label": schedule_label,
        "graviton_interpretation": "NOT ESTABLISHED",
        "long_time_stability": "NOT TESTED",
        "refinement": "NOT TESTED",
        "emergent_time": "OPEN",
        "limiting_speed": "OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

for parity, record in records.items():
    full = record["full_spectrum"]
    print(
        f"[{parity}] full spectral radius="
        f"{max(abs(value) for value in full['raw']['operational']['eigenvalues']):.12g} "
        f"unit/off/open={full['unit_consistent']}/{full['off_unit']}/{full['open']}",
        flush=True,
    )
    if record["shape_spectrum"] is not None:
        shape = record["shape_spectrum"]
        print(
            f"[{parity}] shape spectral radius="
            f"{max(abs(value) for value in shape['raw']['operational']['eigenvalues']):.12g} "
            f"unit/off/open={shape['unit_consistent']}/{shape['off_unit']}/{shape['open']}",
            flush=True,
        )
if schedule_record is not None:
    print(
        "schedule {} eig={:.6e}/{:.6e} svd={:.6e}/{:.6e}".format(
            schedule_label,
            schedule_record["eigen_distance"],
            schedule_record["eigen_uncertainty"],
            schedule_record["singular_distance"],
            schedule_record["singular_uncertainty"],
        ),
        flush=True,
    )
print(f"OUTCOME: {outcome}", flush=True)
print(f"Tests passed: {passed}/{tests}", flush=True)
raise SystemExit(0 if passed == tests else 1)
