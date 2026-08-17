#!/usr/bin/env python3
"""Freeze calibrated anisotropic canonical input/seed pairs before nonlinear solves.

Prior-art commit: 526a202.
Protocol commit: 05f76c3.
"""

import ast
from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import multiprocessing as mp_pool
from pathlib import Path
import sys

import mpmath as arb
import numpy as np


HERE = Path(__file__).resolve().parent
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
TANGENT_INPUT = HERE / "gravity_600cell_dust_dynamic_tangent.json"
ENUMERATION_INPUT = HERE / "gravity_600cell_dust_dynamic_tangent_conjugacy_enumeration.json"
DIRECTION_INPUT = HERE / "gravity_600cell_dust_nonlinear_relative.json"
CANONICAL_SOURCE = HERE / "verify_gravity_600cell_dust_canonical_legendre_rank.py"
ACTION_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
OUTPUT = HERE / "gravity_600cell_dust_nonlinear_boundary_covariance_seeds.json"

PRIOR_ART_COMMIT = "526a202"
PROTOCOL_COMMIT = "05f76c3"
INPUT_HASHES = {
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "tangent": "1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5",
    "enumeration": "51b52457eba84ca1e41926b6e4fb1c51032f788b70bde916a3fb755d0323cb3e",
    "directions": "6e7d108ec7b1a2c80b412134a301084aea14f9457fedba1fd840820ad6f558dd",
    "canonical_source": "396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e",
    "action_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}
DPS = 100
DERIVATIVE_STEPS = {
    "operational": (arb.mpf("1e-20"), arb.mpf("1e-15")),
    "validation": (arb.mpf("3e-20"), arb.mpf("3e-15")),
}
ETA = arb.mpf("1e-4")
ARITHMETIC_FLOOR = arb.mpf("1e-70")

arb.mp.dps = DPS


def digest(path):
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
    }
    tree = ast.parse(CANONICAL_SOURCE.read_text(), filename=str(CANONICAL_SOURCE))
    body = [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in wanted]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch: missing={wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]),
                str(CANONICAL_SOURCE), "exec"),
        globals(),
    )


def vector_norm(vector):
    return arb.sqrt(sum(abs(value)**2 for value in vector))


def matrix_strings(matrix, digits=60):
    return [[arb.nstr(matrix[row, column], digits)
             for column in range(matrix.cols)]
            for row in range(matrix.rows)]


def vector_strings(vector, digits=60):
    return [arb.nstr(value, digits) for value in vector]


def response_record(matrix):
    canonical = extract_canonical(matrix)
    rhs = arb.matrix(65, 60)
    for row in range(35):
        for column in range(30):
            rhs[row, column] = -matrix[30+row, column]
    for row in range(30):
        for column in range(30):
            rhs[35+row, column] = matrix[row, column]
        rhs[35+row, 30+row] = 1
    response = (canonical**-1)*rhs
    return canonical, response


def permute_vector(permutation, vector):
    result = [arb.mpf(0) for _ in range(30)]
    for source, target in enumerate(permutation):
        result[target] = vector[source]
    return tuple(result)


paths = {
    "tick": TICK_INPUT,
    "tangent": TANGENT_INPUT,
    "enumeration": ENUMERATION_INPUT,
    "directions": DIRECTION_INPUT,
    "canonical_source": CANONICAL_SOURCE,
    "action_source": ACTION_SOURCE,
}
hashes = {name: digest(path) for name, path in paths.items()}
tick = json.loads(TICK_INPUT.read_text())
tangent_input = json.loads(TANGENT_INPUT.read_text())
enumeration = json.loads(ENUMERATION_INPUT.read_text())
direction_input = json.loads(DIRECTION_INPUT.read_text())

load_audited_functions()
spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_nonlinear_covariance_seeds", ACTION_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise

models = {
    parity: augment_boundary_orbits(model)
    for parity, model in gro.models.items()
}

physical_candidates = [
    record for record in enumeration["boundary_candidates"]
    if "IDENTICAL_PHYSICAL_EDGE_SETS" in record["sources"]
]
physical_map = tuple(physical_candidates[0]["permutation"])
raw_directions = direction_input["direction_enumeration"]["directions"]
directions = []
for record in raw_directions:
    vector = [arb.mpf(str(value)) for value in record["boundary_vector"]]
    mean = sum(vector)/30
    vector = [value-mean for value in vector]
    norm = vector_norm(vector)
    vector = tuple(value/norm for value in vector)
    directions.append({"index": int(record["index"]), "even": vector,
                       "odd": permute_vector(physical_map, vector)})

provenance_ok = bool(
    hashes == INPUT_HASHES
    and tick.get("outcome") == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and tick.get("passed") == tick.get("tests") == 7
    and tangent_input.get("passed") == tangent_input.get("tests") == 12
    and enumeration.get("outcome") == "GEOMETRIC_CONJUGACY_CANDIDATES_ENUMERATED"
    and enumeration.get("passed") == enumeration.get("tests") == 7
    and direction_input.get("passed") == direction_input.get("tests") == 10
    and direction_input.get("eta") == 0.0001
)
direction_ok = bool(
    len(physical_candidates) == 1
    and len(directions) == 4
    and direction_input["direction_enumeration"]["absolute_rank_above_1e-8"] == 4
    and all(
        abs(sum(record[parity])) < arb.mpf("1e-80")
        and abs(vector_norm(record[parity])-1) < arb.mpf("1e-80")
        for record in directions for parity in ("even", "odd")
    )
)
carrier_ok = bool(
    gro.tests == gro.passed == 43
    and all(
        len(model["old_orbits"]) == 30
        and len(model["edge_orbits"]) == 35
        and len(model["final_orbits"]) == 30
        for model in models.values()
    )
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

fork_context = mp_pool.get_context("fork")
records = {}
all_base_ok = True
all_branch_ok = True
all_entry_ok = True
all_reciprocity_ok = True
all_rank_ok = True

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing dynamic response Hessians", flush=True)
    scale_log, rho_log = (
        arb.mpf(value) for value in tick["solutions"][parity]["state"]
    )
    rho = ARB_RHO*arb.exp(rho_log)
    diagonal = arb.exp(scale_log)*ARB_L0_SQUARE-rho
    ARB_BASE_OLD = tuple(ARB_L0_SQUARE for _ in range(30))
    ARB_BASE_X = tuple([diagonal]*30+[rho]*5)
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

    tasks = []
    for coordinate in range(95):
        for pair in DERIVATIVE_STEPS.values():
            for step in pair:
                for sign in (1, -1):
                    tasks.append((coordinate, arb.nstr(sign*step, 20)))
    with fork_context.Pool(
        processes=8, initializer=initialize_worker, initargs=(model,)
    ) as pool:
        raw_results = pool.map(gradient_worker, tasks, chunksize=1)

    gradient_values = {}
    branch_ok = True
    minimum_minor = base_branch["minimum_leading_minor"]
    minimum_argument = base_branch["minimum_argument"]
    for task, raw in zip(tasks, raw_results):
        action = unpack_complex(raw["action"])
        gradient = tuple(unpack_complex(value) for value in raw["gradient"])
        branch = unpack_branch(raw["branch"])
        gradient_values[task] = gradient
        maximum_imaginary = max(
            maximum_imaginary, abs(arb.im(action)),
            *(abs(arb.im(value)) for value in gradient),
        )
        minimum_minor = min(minimum_minor, branch["minimum_leading_minor"])
        minimum_argument = min(minimum_argument, branch["minimum_argument"])
        branch_ok &= bool(
            branch["negative_counts"] == Counter({1: 2400})
            and branch["minimum_leading_minor"] > 0
            and branch["minimum_argument"] > arb.mpf("1e-6")
        )
    branch_ok &= maximum_imaginary < arb.mpf("1e-70")
    all_branch_ok &= branch_ok

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
    entry_ok = all(
        abs(d_cross[row, column]) <= 10*(
            abs(d_op[row, column])+abs(d_val[row, column])+ARITHMETIC_FLOOR
        )
        for row in range(95) for column in range(95)
    )
    all_entry_ok &= entry_ok
    hessian_error = (
        spectral_norm(d_op)+spectral_norm(d_val)+spectral_norm(d_cross)
        + ARITHMETIC_FLOOR
    )
    antisymmetric_norm = spectral_norm(
        matrices["operational"]-matrices["operational"].T
    )
    reciprocity_ok = antisymmetric_norm <= 10*hessian_error
    all_reciprocity_ok &= reciprocity_ok

    canonical = {}
    response = {}
    for name, matrix in matrices.items():
        canonical[name], response[name] = response_record(matrix)
    epsilon_j = (
        spectral_norm(canonical["operational"]-canonical["operational_shadow"])
        + spectral_norm(canonical["validation"]-canonical["validation_shadow"])
        + spectral_norm(canonical["operational"]-canonical["validation"])
        + ARITHMETIC_FLOOR
    )
    singular = tuple(arb.svd_r(canonical["operational"], compute_uv=False))
    rank_ok = all(value > 100*epsilon_j for value in singular)
    all_rank_ok &= rank_ok
    response_error = (
        spectral_norm(response["operational"]-response["operational_shadow"])
        + spectral_norm(response["validation"]-response["validation_shadow"])
        + spectral_norm(response["operational"]-response["validation"])
        + ARITHMETIC_FLOOR
    )

    records[parity] = {
        "base_old": ARB_BASE_OLD,
        "base_x": ARB_BASE_X,
        "base_new": ARB_BASE_NEW,
        "base_pre": pre,
        "base_post": post,
        "canonical": {
            "operational": canonical["operational"],
            "validation": canonical["validation"],
        },
        "response": {
            "operational": response["operational"],
            "validation": response["validation"],
        },
        "base_ok": base_ok,
        "branch_ok": branch_ok,
        "entry_ok": entry_ok,
        "reciprocity_ok": reciprocity_ok,
        "rank_ok": rank_ok,
        "epsilon_j": epsilon_j,
        "smallest_singular": singular[-1],
        "response_error": response_error,
        "minimum_minor": minimum_minor,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": maximum_imaginary,
        "hessian_error": hessian_error,
        "antisymmetric_norm": antisymmetric_norm,
    }
    print(
        f"[{parity}] smin={arb.nstr(singular[-1], 9)} "
        f"epsJ={arb.nstr(epsilon_j, 6)} response_err={arb.nstr(response_error, 6)}",
        flush=True,
    )

p_star = abs(sum(records["even"]["base_pre"])/30)
cases = []
for direction in directions:
    for sector in ("POSITION", "MOMENTUM"):
        rays = {}
        response_vectors = {}
        gains = []
        for parity in ("even", "odd"):
            d = direction[parity]
            ray = arb.matrix(60, 1)
            if sector == "POSITION":
                for index in range(30):
                    ray[index] = d[index]
            else:
                for index in range(30):
                    ray[30+index] = p_star*d[index]
            rays[parity] = ray
            response_vectors[parity] = {}
            for variant in ("operational", "validation"):
                value = records[parity]["response"][variant]*ray
                response_vectors[parity][variant] = value
                gains.append(vector_norm(value))
        gain = max(gains)
        amplitude = ETA/gain
        for sign in (-1, 1):
            for level_name, level in (("HALF", arb.mpf("0.5")),
                                      ("FULL", arb.mpf(1))):
                factor = sign*level*amplitude
                case = {
                    "id": f"d{direction['index']}_{sector.lower()}_s{sign:+d}_{level_name.lower()}",
                    "direction_index": direction["index"],
                    "sector": sector,
                    "sign": sign,
                    "level": arb.nstr(level, 10),
                    "gain": arb.nstr(gain, 60),
                    "amplitude": arb.nstr(amplitude, 60),
                    "factor": arb.nstr(factor, 60),
                    "parities": {},
                }
                for parity in ("even", "odd"):
                    ray = tuple(rays[parity][index] for index in range(60))
                    case["parities"][parity] = {
                        "input_ray": vector_strings(ray),
                        "unknown_seed_delta_operational": vector_strings(
                            tuple(factor*response_vectors[parity]["operational"][index]
                                  for index in range(65))
                        ),
                        "unknown_seed_delta_validation": vector_strings(
                            tuple(factor*response_vectors[parity]["validation"][index]
                                  for index in range(65))
                        ),
                    }
                cases.append(case)

case_ok = bool(
    len(cases) == 32
    and len({case["id"] for case in cases}) == 32
    and all(arb.mpf(case["amplitude"]) > 0 for case in cases)
    and all(abs(arb.mpf(case["gain"])*arb.mpf(case["amplitude"])-ETA)
            < arb.mpf("1e-55") for case in cases)
)

tests = [
    ("all frozen input hashes and artifact provenance", provenance_ok),
    ("imported carrier retains 43/43 and 30+35+30 dimensions", carrier_ok),
    ("unique physical-edge map and four normalized directions", direction_ok),
    ("both dynamic backgrounds reproduce", all_base_ok),
    ("all 1522 Hessian evaluations retain the Lorentzian branch", all_branch_ok),
    ("all Hessian entries pass operational/validation calibration", all_entry_ok),
    ("both complete Hessians pass reciprocity", all_reciprocity_ok),
    ("both 65 by 65 canonical matrices have resolved full rank", all_rank_ok),
    ("exactly 32 derived finite paired cases", case_ok),
    ("no nonlinear perturbed case was evaluated or compared", True),
]
passed = sum(bool(ok) for _, ok in tests)

payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "eta": arb.nstr(ETA, 20),
    "p_star": arb.nstr(p_star, 60),
    "physical_edge_permutation": list(physical_map),
    "number_of_directions": 4,
    "number_of_paired_cases": len(cases),
    "nonlinear_perturbed_action_evaluations": 0,
    "nonlinear_outputs_compared": False,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "full_720_edge_carrier": False,
    "directions": [{
        "index": record["index"],
        "even": vector_strings(record["even"]),
        "odd": vector_strings(record["odd"]),
    } for record in directions],
    "parities": {
        parity: {
            "base_old": vector_strings(records[parity]["base_old"]),
            "base_x": vector_strings(records[parity]["base_x"]),
            "base_new": vector_strings(records[parity]["base_new"]),
            "base_pre": vector_strings(records[parity]["base_pre"]),
            "base_post": vector_strings(records[parity]["base_post"]),
            "canonical_operational": matrix_strings(
                records[parity]["canonical"]["operational"]
            ),
            "canonical_validation": matrix_strings(
                records[parity]["canonical"]["validation"]
            ),
            "epsilon_j": arb.nstr(records[parity]["epsilon_j"], 60),
            "smallest_singular": arb.nstr(
                records[parity]["smallest_singular"], 60
            ),
            "response_error": arb.nstr(records[parity]["response_error"], 60),
            "minimum_leading_minor": arb.nstr(
                records[parity]["minimum_minor"], 60
            ),
            "minimum_angle_argument": arb.nstr(
                records[parity]["minimum_argument"], 60
            ),
            "maximum_imaginary": arb.nstr(
                records[parity]["maximum_imaginary"], 60
            ),
        } for parity in ("even", "odd")
    },
    "cases": cases,
    "passed": passed,
    "tests": len(tests),
    "outcome": (
        "NONLINEAR_BOUNDARY_COVARIANCE_CASES_FROZEN"
        if passed == len(tests)
        else "NONLINEAR_BOUNDARY_COVARIANCE_SEED_CONTROL_FAILED"
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(f"cases={len(cases)} p_star={arb.nstr(p_star, 10)}")
print(f"OUTCOME: {payload['outcome']}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) else 1)
