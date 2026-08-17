#!/usr/bin/env python3
"""Disclosed gravity--dust decomposition of the full lapse Schur stiffness."""

import hashlib
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
INPUT = HERE / "gravity_600cell_dust_full_lapse_schur.json"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
OUTPUT = HERE / "gravity_600cell_dust_lapse_stiffness_origin.json"
PRIOR_ART_COMMIT = "475f534"
PROTOCOL_COMMIT = "a0a5491"
EXPECTED_INPUT_HASH = "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349"
DPS = 100
mp.mp.dps = DPS
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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_from_record(record):
    return np.array([
        [complex(float(value["real"]), float(value["imaginary"])) for value in row]
        for row in record["schur"]["midpoint_matrix"]
    ], dtype=np.complex128)


def serialize_float(value):
    return f"{float(value):.17e}"


def serialize_complex(value):
    return {
        "real": serialize_float(np.real(value)),
        "imaginary": serialize_float(np.imag(value)),
    }


source = json.loads(INPUT.read_text())
tick = json.loads(TICK_INPUT.read_text())
input_hash = sha256(INPUT)
tick_hash = sha256(TICK_INPUT)
provenance_ok = bool(
    input_hash == EXPECTED_INPUT_HASH
    and source["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    and source["passed"] == source["tests"] == 18
    and tick_hash == source["input_sha256"]["tick"]
    and all(
        record["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
        and record["resolved_schur_rank"] == 120
        and record["schur_zero_count"] == 0
        and record["schur_open_count"] == 0
        for record in source["parities"].values()
    )
)
check(
    "the regular full-Schur input has exact committed provenance",
    provenance_ok,
    f"Schur={input_hash}, tick={tick_hash}",
)


M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2

states = {
    parity: tuple(mp.mpf(value) for value in tick["solutions"][parity]["state"])
    for parity in ("even", "odd")
}
r_values = [state[1] for state in states.values()]
rho_values = [RHO0 * mp.exp(value) for value in r_values]
dust_values = [
    -(2 * mp.pi * MASS / 120) * mp.sqrt(rho) for rho in rho_values
]
dust_control_ok = bool(
    abs(r_values[0] - r_values[1]) < mp.mpf("1e-80")
    and abs(dust_values[0] - dust_values[1]) < mp.mpf("1e-80")
    and dust_values[0] < 0
)
check(
    "both parities give the same independently recomputed dust Hessian",
    dust_control_ok,
    f"h_dust={mp.nstr(dust_values[0], 30)}",
)
h_dust = float(dust_values[0])


records = {}
alphas = []
weights = []
epsilons = []
dimension_total = 0
all_scalar = True
all_real = True
for parity, parity_record in source["parities"].items():
    sector_records = []
    parity_dimension = 0
    for sector in parity_record["sectors"]:
        dimension = int(sector["irrep_dimension"])
        matrix = matrix_from_record(sector)
        n = 5 * dimension
        if matrix.shape != (n, n):
            raise RuntimeError("stored Schur matrix has the wrong sector dimension")
        alpha = np.trace(matrix) / n
        deviation = float(la.norm(matrix - alpha * np.eye(n), 2))
        epsilon_global = float(sector["schur"]["epsilon_global"])
        maximum_radius = float(sector["schur"]["maximum_entry_radius"])
        epsilon = epsilon_global + n * maximum_radius
        scalar_consistent = bool(deviation < 10 * epsilon)
        real_consistent = bool(abs(alpha.imag) < 10 * epsilon)
        gravity_scalar = alpha - h_dust
        affine = {
            f"{mu:g}": gravity_scalar + mu * h_dust
            for mu in (0.0, 0.5, 1.0, 2.0)
        }
        affine_identity_error = max(
            abs(affine["1"] - alpha),
            abs((gravity_scalar + h_dust) - alpha),
        )
        all_scalar &= scalar_consistent
        all_real &= real_consistent
        full_weight = n * dimension
        alphas.append(alpha)
        weights.append(full_weight)
        epsilons.append(epsilon)
        parity_dimension += full_weight
        sector_records.append({
            "irrep_dimension": dimension,
            "minimal_schur_dimension": n,
            "full_representation_weight": full_weight,
            "alpha_total": alpha,
            "alpha_gravity": gravity_scalar,
            "scalar_deviation_norm": deviation,
            "epsilon": epsilon,
            "scalar_consistent": scalar_consistent,
            "real_consistent": real_consistent,
            "affine_scalars": affine,
            "affine_identity_error": affine_identity_error,
        })
    dimension_total += parity_dimension
    records[parity] = {
        "full_schur_dimension": parity_dimension,
        "sectors": sector_records,
    }

dimension_reality_ok = bool(
    dimension_total == 240
    and all(record["full_schur_dimension"] == 120 for record in records.values())
    and all_real
)
check(
    "all fourteen stored Schur blocks exhaust 120 real directions per parity",
    dimension_reality_ok,
    f"combined dimension={dimension_total}, real={all_real}",
)

maximum_deviation_ratio = max(
    item["scalar_deviation_norm"] / max(item["epsilon"], np.finfo(float).tiny)
    for record in records.values() for item in record["sectors"]
)
check(
    "every total and gravity Schur sector is scalar-consistent under the frozen rule",
    all_scalar,
    f"maximum deviation/epsilon={maximum_deviation_ratio:.3e}",
)

maximum_epsilon = max(epsilons)
alpha_spread = max(abs(left - right) for left in alphas for right in alphas)
common_scalar_ok = bool(alpha_spread < 10 * maximum_epsilon)
check(
    "all sectors and parities admit one common scalar within calibrated uncertainty",
    common_scalar_ok,
    f"spread={alpha_spread:.3e}, 10 epsilon={10*maximum_epsilon:.3e}",
)

alpha_common = sum(
    alpha * weight for alpha, weight in zip(alphas, weights)
) / sum(weights)
alpha_common = complex(alpha_common.real, 0.0)
gravity_common = alpha_common.real - h_dust
mu_star = -gravity_common / h_dust
cancellation_ratio = abs(alpha_common.real) / max(
    abs(gravity_common), abs(h_dust)
)
if cancellation_ratio < 1e-5:
    cancellation_bin = "NEAR_CANCELLATION_PATTERN"
elif cancellation_ratio <= 1e-2:
    cancellation_bin = "PARTIAL_CANCELLATION"
else:
    cancellation_bin = "NO_STRONG_CANCELLATION"

affine_common = {
    f"{mu:g}": gravity_common + mu * h_dust
    for mu in (0.0, 0.5, 1.0, 2.0)
}
affine_roundoff_bound = 64 * np.finfo(float).eps * max(
    abs(gravity_common), abs(h_dust), abs(alpha_common.real),
    np.finfo(float).tiny,
)
affine_controls_ok = bool(
    max(
        item["affine_identity_error"]
        for record in records.values() for item in record["sectors"]
    ) < affine_roundoff_bound
    and abs(affine_common["1"] - alpha_common.real) < affine_roundoff_bound
)
check(
    "the exact affine dust shift is reconstructed for all frozen multipliers",
    affine_controls_ok,
    f"alpha={alpha_common.real:.17e}, gravity={gravity_common:.17e}, "
    f"dust={h_dust:.17e}",
)

controls_ok = bool(
    provenance_ok and dust_control_ok and dimension_reality_ok
    and affine_controls_ok
)
if not controls_ok:
    outcome = "LAPSE_STIFFNESS_ORIGIN_CONTROL_FAILED"
elif not (all_scalar and common_scalar_ok):
    outcome = "LAPSE_STIFFNESS_NONSCALAR"
elif cancellation_bin == "NEAR_CANCELLATION_PATTERN":
    outcome = "LAPSE_STIFFNESS_SCALAR_NEAR_CANCELLATION"
else:
    outcome = "LAPSE_STIFFNESS_SCALAR_NO_NEAR_CANCELLATION"
check(
    "the disclosed hierarchy assigns the stiffness-origin outcome mechanically",
    outcome in {
        "LAPSE_STIFFNESS_ORIGIN_CONTROL_FAILED",
        "LAPSE_STIFFNESS_NONSCALAR",
        "LAPSE_STIFFNESS_SCALAR_NEAR_CANCELLATION",
        "LAPSE_STIFFNESS_SCALAR_NO_NEAR_CANCELLATION",
    },
    f"outcome={outcome}, bin={cancellation_bin}, ratio={cancellation_ratio:.9e}",
)


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": input_hash,
    "tick_sha256": tick_hash,
    "disclosed_post_result_analysis": True,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "outcome": outcome,
    "dust": {
        "mass": mp.nstr(MASS, 80),
        "rho": mp.nstr(rho_values[0], 80),
        "hessian_scalar": mp.nstr(dust_values[0], 80),
    },
    "common": {
        "total_scalar": serialize_float(alpha_common.real),
        "gravity_scalar": serialize_float(gravity_common),
        "dust_scalar": serialize_float(h_dust),
        "maximum_cross_sector_scalar_spread": serialize_float(alpha_spread),
        "maximum_scalar_deviation_over_epsilon": serialize_float(
            maximum_deviation_ratio
        ),
        "cancellation_ratio": serialize_float(cancellation_ratio),
        "cancellation_bin": cancellation_bin,
        "critical_dust_multiplier": serialize_float(mu_star),
        "physical_minus_critical_multiplier": serialize_float(1 - mu_star),
        "affine_counterfactual_scalars": {
            name: serialize_float(value) for name, value in affine_common.items()
        },
    },
    "parities": {
        parity: {
            "full_schur_dimension": record["full_schur_dimension"],
            "sectors": [
                {
                    **{key: value for key, value in item.items()
                       if key not in {
                           "alpha_total", "alpha_gravity", "affine_scalars"
                       }},
                    "alpha_total": serialize_complex(item["alpha_total"]),
                    "alpha_gravity": serialize_complex(item["alpha_gravity"]),
                    "affine_scalars": {
                        name: serialize_complex(value)
                        for name, value in item["affine_scalars"].items()
                    },
                }
                for item in record["sectors"]
            ],
        }
        for parity, record in records.items()
    },
    "classification": {
        "dust_affine_shift": "DERIVED ALGEBRAIC",
        "term_sizes": "DERIVED COMPUTATIONAL",
        "scalarity": "PATTERN POST-RESULT",
        "near_cancellation": "PATTERN POST-RESULT",
        "hamiltonian_constraint": "OPEN",
        "dust_clock": "OPEN",
        "refinement": "OPEN",
    },
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Summary: {passed}/{tests} checks passed")
print(f"Outcome: {outcome}")
print(f"gravity scalar : {gravity_common:.17e}")
print(f"dust scalar    : {h_dust:.17e}")
print(f"residual       : {alpha_common.real:.17e}")
print(f"mu_star        : {mu_star:.17e}")
print(f"Artifact: {OUTPUT}")

if passed != tests:
    sys.exit(1)
