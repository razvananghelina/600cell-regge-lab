#!/usr/bin/env python3
"""High-precision correction audit of the blind dynamic shape spectra.

Prior-art commit: 16f4310.
Protocol commit: 0e8f561.
"""

import hashlib
import json
import math
from pathlib import Path

import mpmath as arb
import numpy as np
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
INPUT = HERE / "gravity_600cell_dust_dynamic_tangent.json"
OUTPUT = HERE / "gravity_600cell_dust_dynamic_tangent_precision.json"

PRIOR_ART_COMMIT = "16f4310"
PROTOCOL_COMMIT = "0e8f561"
INPUT_SHA256 = "1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5"
LOW_DPS = 100
HIGH_DPS = 160
NORMALIZATION = 2**20
STORE_ERROR = arb.mpf("3e-49")
CONSISTENT_FACTOR = arb.mpf(10)
DEPENDENT_FACTOR = arb.mpf(100)
EIG_RESIDUAL_TOLERANCE = arb.mpf("1e-80")

arb.mp.dps = HIGH_DPS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frobenius(matrix):
    return arb.sqrt(sum(abs(matrix[row, column])**2
                        for row in range(matrix.rows)
                        for column in range(matrix.cols)))


def boundary_basis():
    scale = arb.matrix([1/arb.sqrt(30) for _ in range(30)])
    vector = scale.copy()
    vector[0] += 1
    reflector = arb.eye(30)-2*(vector*vector.T)/(vector.T*vector)[0]
    basis = arb.matrix(30, 30)
    for row in range(30):
        basis[row, 0] = scale[row]
        for column in range(1, 30):
            basis[row, column] = reflector[row, column]
    return basis


def phase_basis():
    boundary = boundary_basis()
    phase = arb.matrix(60, 60)
    for block in (0, 30):
        for row in range(30):
            for column in range(30):
                phase[block+row, block+column] = boundary[row, column]
    return phase


PHASE_BASIS = phase_basis()
SCALE = (0, 30)
SHAPE = tuple(range(1, 30))+tuple(range(31, 60))


def read_matrix(rows):
    return arb.matrix([[arb.mpf(value) for value in row] for row in rows])


def select(matrix, indices):
    return arb.matrix([[matrix[row, column] for column in indices]
                       for row in indices])


def reconstruct_shape(tangent):
    sectors = PHASE_BASIS.T*tangent*PHASE_BASIS
    shape = select(sectors, SHAPE)
    mixing = []
    for row in SCALE:
        for column in SHAPE:
            mixing.append(sectors[row, column])
    for row in SHAPE:
        for column in SCALE:
            mixing.append(sectors[row, column])
    return shape, arb.sqrt(sum(abs(value)**2 for value in mixing))


def power_traces(matrix, dps):
    with arb.workdps(dps):
        normalized = matrix/arb.mpf(NORMALIZATION)
        power = arb.eye(matrix.rows)
        values = []
        for _ in range(matrix.rows):
            power = power*normalized
            values.append(+sum(power[index, index]
                               for index in range(matrix.rows)))
        return values


def optimal_match(left, right):
    cost = np.asarray([
        [float(abs(a-b)) for b in right]
        for a in left
    ], dtype=np.float64)
    rows, columns = linear_sum_assignment(cost)
    pairs = sorted(zip(rows.tolist(), columns.tolist()))
    distances = [abs(left[row]-right[column]) for row, column in pairs]
    return max(distances), pairs


def eigensystem(matrix, dps, diagnostics=False):
    with arb.workdps(dps):
        eigenvalues, vectors = arb.eig(matrix, left=False, right=True)
        values = [+value for value in eigenvalues]
        vectors = arb.matrix(vectors)
        record = {"values": values, "vectors": vectors}
        if diagnostics:
            diagonal = arb.diag(values)
            residual = frobenius(matrix*vectors-vectors*diagonal)
            residual /= max(arb.mpf(1), frobenius(matrix)*frobenius(vectors))
            inverse = vectors**-1
            condition = frobenius(vectors)*frobenius(inverse)
            record.update({
                "residual": +residual,
                "condition": +condition,
                "invertible": True,
            })
        return record


def classify(distance, uncertainty, consistent_label, dependent_label,
             open_label):
    if distance <= CONSISTENT_FACTOR*uncertainty:
        return consistent_label
    if distance > DEPENDENT_FACTOR*uncertainty:
        return dependent_label
    return open_label


def decimal(value, digits=80):
    return arb.nstr(value, digits, strip_zeros=False)


def complex_record(value):
    return {
        "real": decimal(arb.re(value)),
        "imaginary": decimal(arb.im(value)),
    }


artifact = json.loads(INPUT.read_text())
provenance_ok = bool(
    digest(INPUT) == INPUT_SHA256
    and artifact.get("prior_art_commit") == "25722d9"
    and artifact.get("protocol_commit") == "0bceb9b"
    and artifact.get("outcome") == "DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT"
    and artifact.get("passed") == 12
    and artifact.get("tests") == 12
    and artifact.get("number_of_maps") == 2
    and artifact.get("continuum_target_parsed") is False
    and artifact.get("speed_target_parsed") is False
    and artifact.get("full_720_edge_carrier") is False
)

records = {}
dimensions_ok = True
mixing_ok = True
for parity in ("even", "odd"):
    source = artifact["parities"][parity]
    rows = source["tangent_matrix"]
    dimensions_ok &= bool(
        len(rows) == 60 and all(len(row) == 60 for row in rows)
    )
    tangent = read_matrix(rows)
    shape, mixing = reconstruct_shape(tangent)
    epsilon_mix = arb.mpf(source["epsilon_mix"])
    mixing_ok &= mixing <= 10*epsilon_mix
    delta = arb.mpf(source["epsilon_t"])+STORE_ERROR
    sigma = (
        arb.mpf(source["shape_spectrum"]["singular_values"][0])
        + arb.mpf(source["shape_spectrum"]["epsilon_svd"])
    )
    print(f"[{parity}] computing 58 normalized power traces", flush=True)
    traces_low = power_traces(shape, LOW_DPS)
    traces_high = power_traces(shape, HIGH_DPS)
    records[parity] = {
        "shape": shape,
        "mixing": mixing,
        "epsilon_mix": epsilon_mix,
        "delta": delta,
        "sigma": sigma,
        "traces_low": traces_low,
        "traces_high": traces_high,
    }

trace_rows = []
trace_precision_ok = True
for index in range(58):
    k = index+1
    even = records["even"]
    odd = records["odd"]
    difference = abs(even["traces_high"][index]-odd["traces_high"][index])
    uncertainty = arb.mpf("1e-140")
    for record in (even, odd):
        alpha = (record["sigma"]+record["delta"])/NORMALIZATION
        uncertainty += (
            58*k*alpha**(k-1)*(record["delta"]/NORMALIZATION)
            + abs(record["traces_low"][index]-record["traces_high"][index])
        )
    trace_precision_ok &= bool(arb.isfinite(difference)
                               and arb.isfinite(uncertainty)
                               and uncertainty > 0)
    trace_rows.append({
        "k": k,
        "difference": difference,
        "uncertainty": uncertainty,
        "ratio": difference/uncertainty,
    })

if all(row["difference"] <= 10*row["uncertainty"] for row in trace_rows):
    trace_label = "TRACE_SPECTRUM_CONSISTENT"
elif any(row["difference"] > 100*row["uncertainty"] for row in trace_rows):
    trace_label = "TRACE_SPECTRUM_DEPENDENT"
else:
    trace_label = "TRACE_SPECTRUM_OPEN"

for parity in ("even", "odd"):
    print(f"[{parity}] computing the 100-digit eigensystem", flush=True)
    low = eigensystem(records[parity]["shape"], LOW_DPS)
    print(f"[{parity}] computing the 160-digit eigensystem", flush=True)
    high = eigensystem(records[parity]["shape"], HIGH_DPS, diagnostics=True)
    convergence, convergence_pairs = optimal_match(low["values"], high["values"])
    eigen_uncertainty = (
        high["condition"]*records[parity]["delta"]
        + convergence+arb.mpf("1e-100")
    )
    records[parity].update({
        "eig_low": low,
        "eig_high": high,
        "eig_convergence": convergence,
        "eig_convergence_pairs": convergence_pairs,
        "eig_uncertainty": eigen_uncertainty,
    })
    print(
        f"[{parity}] residual={decimal(high['residual'], 8)} "
        f"kappa_F={decimal(high['condition'], 8)}",
        flush=True,
    )

eigen_distance, eigen_pairs = optimal_match(
    records["even"]["eig_high"]["values"],
    records["odd"]["eig_high"]["values"],
)
eigen_uncertainty = (
    records["even"]["eig_uncertainty"]
    + records["odd"]["eig_uncertainty"]
)
eigen_label = classify(
    eigen_distance,
    eigen_uncertainty,
    "EIG_SPECTRUM_CONSISTENT",
    "EIG_SPECTRUM_DEPENDENT",
    "EIG_SPECTRUM_OPEN",
)

if (trace_label == "TRACE_SPECTRUM_DEPENDENT"
        and eigen_label == "EIG_SPECTRUM_DEPENDENT"):
    outcome = "SCHEDULE_DEPENDENCE_CONFIRMED_HIGH_PRECISION"
elif (trace_label == "TRACE_SPECTRUM_CONSISTENT"
      and eigen_label == "EIG_SPECTRUM_CONSISTENT"):
    outcome = "SCHEDULE_SPECTRUM_NOT_RESOLVED"
else:
    outcome = "SCHEDULE_SPECTRUM_OPEN_NONNORMAL"

residuals_ok = all(
    records[parity]["eig_high"]["residual"] < EIG_RESIDUAL_TOLERANCE
    for parity in ("even", "odd")
)
conditions_ok = all(
    records[parity]["eig_high"]["invertible"]
    and arb.isfinite(records[parity]["eig_high"]["condition"])
    and records[parity]["eig_high"]["condition"] >= 1
    for parity in ("even", "odd")
)
classification_ok = bool(
    trace_label in {
        "TRACE_SPECTRUM_CONSISTENT",
        "TRACE_SPECTRUM_DEPENDENT",
        "TRACE_SPECTRUM_OPEN",
    }
    and eigen_label in {
        "EIG_SPECTRUM_CONSISTENT",
        "EIG_SPECTRUM_DEPENDENT",
        "EIG_SPECTRUM_OPEN",
    }
)

tests = [
    ("frozen input provenance", provenance_ok),
    ("two 60 by 60 tangent matrices", dimensions_ok),
    ("scale-shape reconstruction", mixing_ok),
    ("all 58 trace comparisons finite", trace_precision_ok),
    ("even 160-digit eigenpair residual", records["even"]["eig_high"]["residual"]
     < EIG_RESIDUAL_TOLERANCE),
    ("odd 160-digit eigenpair residual", records["odd"]["eig_high"]["residual"]
     < EIG_RESIDUAL_TOLERANCE),
    ("eigenvector matrices invertible", conditions_ok),
    ("frozen classifiers evaluated", classification_ok),
]
passed = sum(bool(ok) for _, ok in tests)
if passed != len(tests):
    outcome = "PRECISION_AUDIT_CONTROL_FAILED"

output = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": INPUT_SHA256,
    "low_dps": LOW_DPS,
    "high_dps": HIGH_DPS,
    "normalization": NORMALIZATION,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "full_720_edge_carrier": False,
    "passed": passed,
    "tests": len(tests),
    "outcome": outcome,
    "trace": {
        "label": trace_label,
        "maximum_ratio": decimal(max(row["ratio"] for row in trace_rows)),
        "rows": [{
            "k": row["k"],
            "difference": decimal(row["difference"]),
            "uncertainty": decimal(row["uncertainty"]),
            "ratio": decimal(row["ratio"]),
        } for row in trace_rows],
    },
    "eigen": {
        "label": eigen_label,
        "distance": decimal(eigen_distance),
        "uncertainty": decimal(eigen_uncertainty),
        "ratio": decimal(eigen_distance/eigen_uncertainty),
        "matching": [[int(row), int(column)] for row, column in eigen_pairs],
    },
    "parities": {},
}

for parity in ("even", "odd"):
    record = records[parity]
    high = record["eig_high"]
    output["parities"][parity] = {
        "mixing_norm": decimal(record["mixing"]),
        "epsilon_mix": decimal(record["epsilon_mix"]),
        "delta_matrix": decimal(record["delta"]),
        "largest_singular_bound": decimal(record["sigma"]),
        "eigen_residual": decimal(high["residual"]),
        "eigenvector_condition_frobenius": decimal(high["condition"]),
        "eigen_convergence": decimal(record["eig_convergence"]),
        "eigen_uncertainty": decimal(record["eig_uncertainty"]),
        "eigenvalues": [complex_record(value) for value in high["values"]],
    }

OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")

for label, ok in tests:
    print(f"{'PASS' if ok else 'FAIL'}: {label}")
print(
    "trace {} max_ratio={}".format(
        trace_label, decimal(max(row["ratio"] for row in trace_rows), 10)
    )
)
print(
    "eigen {} distance={} uncertainty={} ratio={}".format(
        eigen_label,
        decimal(eigen_distance, 10),
        decimal(eigen_uncertainty, 10),
        decimal(eigen_distance/eigen_uncertainty, 10),
    )
)
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) else 1)
