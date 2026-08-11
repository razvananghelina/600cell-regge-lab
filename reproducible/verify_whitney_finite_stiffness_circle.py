#!/usr/bin/env python3
"""Known-answer S1 topology gate for finite Whitney trace stiffness.

Protocol commit d499ca0 froze the exact Bloch formula, grids, stiffness laws,
numerical gates, and kill boundary before evaluation.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy import linalg
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_finite_stiffness_circle.json")
PROTOCOL_COMMIT = "d499ca0"
LEVELS = (8, 16, 32, 64, 128)
NUMERICAL_LEVELS = (8, 16)
NUMERICAL_KAPPAS = (1.0, 2.0, 4.0)
MATCH_GATE = 1e-11
ZERO_TOLERANCE = 1e-9
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def local_circle_pencil(vertices, kappa):
    h = 1.0 / vertices
    local_metric = np.asarray((
        (h / 3.0, h / 6.0, 0.0),
        (h / 6.0, h / 3.0, 0.0),
        (0.0, 0.0, 1.0 / h),
    ))
    local_weak = np.asarray((
        (0.0, 0.0, -1.0 / h),
        (0.0, 0.0, 1.0 / h),
        (-1.0 / h, 1.0 / h, 0.0),
    ))
    metric = linalg.block_diag(*([local_metric] * vertices))
    weak = linalg.block_diag(*([local_weak] * vertices))
    jump = np.zeros((vertices, 3 * vertices))
    for vertex in range(vertices):
        jump[vertex, 3 * vertex] = 1.0
        jump[vertex, 3 * ((vertex - 1) % vertices) + 1] = -1.0
    penalty = jump.T @ jump
    return weak + kappa * penalty, metric, weak, penalty, jump


def assembled_circle_pencil(vertices):
    h = 1.0 / vertices
    differential = np.zeros((vertices, vertices))
    for edge in range(vertices):
        differential[edge, edge] = -1.0
        differential[edge, (edge + 1) % vertices] = 1.0
    mass_zero = np.zeros((vertices, vertices))
    for edge in range(vertices):
        left, right = edge, (edge + 1) % vertices
        mass_zero[left, left] += h / 3.0
        mass_zero[right, right] += h / 3.0
        mass_zero[left, right] += h / 6.0
        mass_zero[right, left] += h / 6.0
    mass_one = np.eye(vertices) / h
    metric = linalg.block_diag(mass_zero, mass_one)
    forward = mass_one @ differential
    weak = np.block([
        [np.zeros_like(mass_zero), forward.T],
        [forward, np.zeros_like(mass_one)],
    ])
    return weak, metric


def q_zero_basis(vertices):
    basis = np.zeros((3 * vertices, 3))
    normalization = 1.0 / math.sqrt(vertices)
    for edge in range(vertices):
        basis[3 * edge, 0] = normalization
        basis[3 * edge + 1, 1] = normalization
        basis[3 * edge + 2, 2] = normalization
    return basis


def exact_q_zero(vertices, kappa):
    root = math.sqrt(36.0 * kappa * kappa + 12.0)
    return np.asarray((
        (6.0 * kappa - root) * vertices,
        0.0,
        (6.0 * kappa + root) * vertices,
    ))


def harmonic_mass(vertices, kappa):
    # Rationalized to avoid cancellation at large kappa.
    return (
        12.0 * vertices
        / (math.sqrt(36.0 * kappa * kappa + 12.0) + 6.0 * kappa)
    )


print("=" * 78)
print("FINITE WHITNEY TRACE STIFFNESS ON THE UNIT CIRCLE")
print("=" * 78)

kappa_symbol = sy.symbols("kappa", positive=True)
predicted_block = sy.Matrix((
    (12 * kappa_symbol, -sy.sqrt(12)),
    (-sy.sqrt(12), 0),
))
predicted_eigenvalues = predicted_block.eigenvals()
expected_eigenvalues = {
    6 * kappa_symbol - sy.sqrt(36 * kappa_symbol ** 2 + 12): 1,
    6 * kappa_symbol + sy.sqrt(36 * kappa_symbol ** 2 + 12): 1,
}
spectral_variable = sy.symbols("lambda")
expected_characteristic = sy.expand(
    (spectral_variable - next(iter(expected_eigenvalues)))
    * (spectral_variable - tuple(expected_eigenvalues)[1])
)
observed_characteristic = predicted_block.charpoly(
    spectral_variable
).as_expr()
check("the preregistered mass-orthonormal q=0 block has the exact roots",
      sy.simplify(observed_characteristic - expected_characteristic) == 0,
      str(predicted_eigenvalues))

numerical_records = []
all_hermitian_positive = True
all_zero_counts = True
all_matches = True
all_assembled_zeros = True
maximum_match_error = 0.0
for vertices in NUMERICAL_LEVELS:
    assembled_weak, assembled_metric = assembled_circle_pencil(vertices)
    assembled_values = linalg.eigvalsh(assembled_weak, assembled_metric)
    assembled_zero_count = int(np.count_nonzero(
        np.abs(assembled_values) < ZERO_TOLERANCE
    ))
    all_assembled_zeros &= assembled_zero_count == 2
    for kappa in NUMERICAL_KAPPAS:
        pencil, metric, weak, penalty, jump = local_circle_pencil(
            vertices, kappa
        )
        all_hermitian_positive &= (
            np.array_equal(pencil, pencil.T)
            and np.min(np.linalg.eigvalsh(metric)) > 0.0
            and np.min(np.linalg.eigvalsh(penalty)) > -1e-12
        )
        values = linalg.eigvalsh(pencil, metric)
        zero_count = int(np.count_nonzero(np.abs(values) < ZERO_TOLERANCE))
        all_zero_counts &= zero_count == 1

        basis = q_zero_basis(vertices)
        compressed_weak = basis.T @ pencil @ basis
        compressed_metric = basis.T @ metric @ basis
        q_values = linalg.eigvalsh(compressed_weak, compressed_metric)
        expected = exact_q_zero(vertices, kappa)
        relative_error = float(np.max(
            np.abs(q_values - expected) / np.maximum(1.0, np.abs(expected))
        ))
        maximum_match_error = max(maximum_match_error, relative_error)
        full_match_error = max(
            min(abs(values - expected_value))
            / max(1.0, abs(expected_value))
            for expected_value in expected
        )
        all_matches &= (
            relative_error < MATCH_GATE and full_match_error < MATCH_GATE
        )
        numerical_records.append({
            "vertices": vertices,
            "kappa": kappa,
            "local_zero_count": zero_count,
            "assembled_zero_count": assembled_zero_count,
            "q_zero_eigenvalues": q_values.tolist(),
            "exact_q_zero_eigenvalues": expected.tolist(),
            "q_zero_maximum_relative_error": relative_error,
            "full_spectrum_match_relative_error": full_match_error,
        })

check("all complete local pencils are Hermitian with positive metric/penalty",
      all_hermitian_positive)
check("every finite-kappa local pencil has exactly one zero mode",
      all_zero_counts)
check("all complete spectra contain the exact q=0 Bloch eigenvalues",
      all_matches,
      f"maximum compressed relative error={maximum_match_error:.3e}")
check("exact conforming assembly independently restores two Betti zero modes",
      all_assembled_zeros)

asymptotic_records = []
for vertices in LEVELS:
    laws = {
        "fixed_kappa_1": 1.0,
        "kappa_equals_N": float(vertices),
        "kappa_equals_N_squared": float(vertices ** 2),
    }
    asymptotic_records.append({
        "vertices": vertices,
        "h": 1.0 / vertices,
        **{
            label: {
                "kappa": kappa,
                "kappa_times_h": kappa / vertices,
                "lifted_harmonic_absolute_eigenvalue": harmonic_mass(
                    vertices, kappa
                ),
                "positive_stiff_eigenvalue": exact_q_zero(
                    vertices, kappa
                )[2],
            }
            for label, kappa in laws.items()
        },
    })

fixed_masses = [record["fixed_kappa_1"][
    "lifted_harmonic_absolute_eigenvalue"] for record in asymptotic_records]
linear_masses = [record["kappa_equals_N"][
    "lifted_harmonic_absolute_eigenvalue"] for record in asymptotic_records]
quadratic_masses = [record["kappa_equals_N_squared"][
    "lifted_harmonic_absolute_eigenvalue"] for record in asymptotic_records]

N = sy.symbols("N", positive=True)
def exact_mass(kappa):
    return sy.simplify(
        12 * N / (sy.sqrt(36 * kappa ** 2 + 12) + 6 * kappa)
    )

limits = {
    "fixed_kappa_1": sy.limit(exact_mass(1), N, sy.oo),
    "kappa_equals_N": sy.limit(exact_mass(N), N, sy.oo),
    "kappa_equals_N_squared": sy.limit(exact_mass(N ** 2), N, sy.oo),
}
check("fixed finite kappa makes the lifted harmonic mass diverge exactly",
      limits["fixed_kappa_1"] == sy.oo
      and all(fixed_masses[index + 1] > fixed_masses[index]
              for index in range(len(fixed_masses) - 1)),
      f"masses={fixed_masses}")
check("kappa proportional to 1/h leaves a nonzero unit harmonic mass",
      limits["kappa_equals_N"] == 1
      and all(abs(value - 1.0) < 0.01 for value in linear_masses[-2:]),
      f"masses={linear_masses}")
check("only kappa*h -> infinity recovers the harmonic zero in this family",
      limits["kappa_equals_N_squared"] == 0
      and all(quadratic_masses[index + 1] < quadratic_masses[index]
              for index in range(len(quadratic_masses) - 1)),
      f"masses={quadratic_masses}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "exact_q_zero": {
        "mass_orthonormal_block_times_h": (
            "[[12*kappa,-sqrt(12)],[-sqrt(12),0]]"
        ),
        "eigenvalues_times_h": (
            "6*kappa +/- sqrt(36*kappa^2+12), plus scalar zero"
        ),
        "lifted_harmonic_mass": (
            "12/(h*(sqrt(36*kappa^2+12)+6*kappa))"
        ),
    },
    "symbolic_limits": {key: str(value) for key, value in limits.items()},
    "numerical_records": numerical_records,
    "asymptotic_records": asymptotic_records,
    "verdicts": [
        "DERIVED NEGATIVE: fixed finite kappa loses the S1 harmonic one-form",
        "DERIVED NEGATIVE: kappa proportional to 1/h leaves a nonzero topological mass",
        "DERIVED: harmonic recovery requires kappa*h tending to infinity or new consistency fluxes",
        "OPEN: a geometry-selected flux completion and its 3D Kähler-Dirac limit",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured circle certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("SYMBOLIC_LIMITS=" + str(limits))
print("FIXED_MASSES=" + str(fixed_masses))
print("LINEAR_MASSES=" + str(linear_masses))
print("QUADRATIC_MASSES=" + str(quadratic_masses))
raise SystemExit(0 if passed == tests else 1)
