#!/usr/bin/env python3
"""Unique copy-symmetric central Whitney flux on the unit circle.

Protocol commit 2b83f31 froze the averaging rule, support, grids, kappa set,
known-answer comparison, and acceptance boundary before evaluation.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy import linalg
import sympy as sy


OUTPUT = Path(__file__).with_name("whitney_central_flux_circle.json")
PROTOCOL_COMMIT = "2b83f31"
LEVELS = (8, 16, 32, 64, 128, 256)
FULL_LEVELS = (8, 16)
KAPPAS = (0.5, 1.0, 2.0, 4.0)
ZERO_TOLERANCE = 1e-9
FULL_MATCH_GATE = 1e-10
FINAL_IMPROVEMENT_GATE = 16.0
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def circle_maps(vertices, symbolic=False):
    dtype = object if symbolic else np.float64
    half = sy.Rational(1, 2) if symbolic else 0.5
    one = sy.Integer(1) if symbolic else 1.0
    zero = sy.Integer(0) if symbolic else 0.0
    injection = sy.zeros(2 * vertices, vertices) if symbolic else np.zeros(
        (2 * vertices, vertices), dtype=dtype
    )
    average = sy.zeros(vertices, 2 * vertices) if symbolic else np.zeros(
        (vertices, 2 * vertices), dtype=dtype
    )
    jump = sy.zeros(vertices, 2 * vertices) if symbolic else np.zeros(
        (vertices, 2 * vertices), dtype=dtype
    )
    differential = sy.zeros(vertices, vertices) if symbolic else np.zeros(
        (vertices, vertices), dtype=dtype
    )
    for vertex in range(vertices):
        injection[2 * vertex, vertex] = one
        injection[2 * vertex + 1, (vertex + 1) % vertices] = one
        average[vertex, 2 * vertex] = half
        average[vertex, 2 * ((vertex - 1) % vertices) + 1] = half
        jump[vertex, 2 * vertex] = one
        jump[vertex, 2 * ((vertex - 1) % vertices) + 1] = -one
        differential[vertex, vertex] = -one
        differential[vertex, (vertex + 1) % vertices] = one
    return injection, average, jump, differential


def central_pencil(vertices, kappa):
    h = 1.0 / vertices
    injection, average, jump, differential = circle_maps(vertices)
    vertex_metric = np.zeros((2 * vertices, 2 * vertices))
    local_zero = h * np.asarray(((1.0 / 3.0, 1.0 / 6.0),
                                 (1.0 / 6.0, 1.0 / 3.0)))
    for edge in range(vertices):
        indices = (2 * edge, 2 * edge + 1)
        vertex_metric[np.ix_(indices, indices)] += local_zero
    edge_metric = np.eye(vertices) / h
    metric = linalg.block_diag(vertex_metric, edge_metric)
    assembled_forward = edge_metric @ differential
    central_forward = assembled_forward @ average
    weak = np.block([
        [np.zeros_like(vertex_metric), central_forward.T],
        [central_forward, np.zeros_like(edge_metric)],
    ])
    vertex_penalty = jump.T @ jump
    penalty = linalg.block_diag(vertex_penalty, np.zeros_like(edge_metric))
    pencil = weak + kappa * penalty
    return {
        "pencil": pencil,
        "metric": metric,
        "weak": weak,
        "penalty": penalty,
        "injection": injection,
        "average": average,
        "jump": jump,
        "differential": differential,
        "assembled_forward": assembled_forward,
        "central_forward": central_forward,
    }


def no_flux_pencil(vertices, kappa, data):
    h = 1.0 / vertices
    local_forward = np.zeros((vertices, 2 * vertices))
    for edge in range(vertices):
        local_forward[edge, 2 * edge] = -1.0 / h
        local_forward[edge, 2 * edge + 1] = 1.0 / h
    vertex_size = 2 * vertices
    edge_size = vertices
    weak = np.block([
        [np.zeros((vertex_size, vertex_size)), local_forward.T],
        [local_forward, np.zeros((edge_size, edge_size))],
    ])
    return weak + kappa * data["penalty"]


def bloch_basis(vertices, momentum):
    phase = np.exp(1j * momentum * np.arange(vertices)) / math.sqrt(vertices)
    basis = np.zeros((3 * vertices, 3), dtype=np.complex128)
    for edge in range(vertices):
        basis[2 * edge, 0] = phase[edge]
        basis[2 * edge + 1, 1] = phase[edge]
        basis[2 * vertices + edge, 2] = phase[edge]
    return basis


def bloch_spectrum(data, momentum):
    vertices = data["differential"].shape[0]
    basis = bloch_basis(vertices, momentum)
    weak = basis.conj().T @ data["pencil"] @ basis
    metric = basis.conj().T @ data["metric"] @ basis
    return linalg.eigvalsh(weak, metric)


print("=" * 78)
print("COPY-SYMMETRIC CENTRAL WHITNEY FLUX ON THE UNIT CIRCLE")
print("=" * 78)

left_weight, right_weight = sy.symbols("left_weight right_weight")
unique_weights = sy.solve((
    left_weight - right_weight,
    left_weight + right_weight - 1,
), (left_weight, right_weight), dict=True)
check("copy exchange symmetry and left inversion uniquely force 1/2,1/2",
      unique_weights == [{left_weight: sy.Rational(1, 2),
                          right_weight: sy.Rational(1, 2)}],
      str(unique_weights))

exact_gates = True
for vertices in FULL_LEVELS:
    injection, average, jump, differential = circle_maps(
        vertices, symbolic=True
    )
    forward = vertices * differential
    exact_gates &= (
        average * injection == sy.eye(vertices)
        and jump * injection == sy.zeros(vertices, vertices)
        and (forward * average) * injection == forward
    )
check("LJ=I, RJ=0 and central F times J=assembled F exactly",
      exact_gates)

support_ok = True
odd_ok = True
zero_vectors_ok = True
zero_counts_ok = True
negative_control_ok = True
q_zero_ok = True
full_match_error = 0.0
full_records = []
for vertices in FULL_LEVELS:
    for kappa in KAPPAS:
        data = central_pencil(vertices, kappa)
        vertex_size = 2 * vertices
        parity = np.diag(np.r_[np.ones(vertex_size), -np.ones(vertices)])
        odd_ok &= np.max(np.abs(
            parity @ data["weak"] + data["weak"] @ parity
        )) < 1e-13

        for edge in range(vertices):
            observed = set(np.flatnonzero(
                np.abs(data["central_forward"][edge]) > 0.0
            ))
            allowed = {
                2 * edge,
                2 * ((edge - 1) % vertices) + 1,
                2 * ((edge + 1) % vertices),
                2 * edge + 1,
            }
            support_ok &= observed <= allowed

        constant_scalar = np.r_[np.ones(vertex_size), np.zeros(vertices)]
        constant_one_form = np.r_[np.zeros(vertex_size), np.ones(vertices)]
        zero_vectors_ok &= (
            np.max(np.abs(data["pencil"] @ constant_scalar)) < 1e-12
            and np.max(np.abs(data["pencil"] @ constant_one_form)) < 1e-12
        )

        values = linalg.eigvalsh(data["pencil"], data["metric"])
        zero_count = int(np.count_nonzero(np.abs(values) < ZERO_TOLERANCE))
        zero_counts_ok &= zero_count == 2
        negative_values = linalg.eigvalsh(
            no_flux_pencil(vertices, kappa, data), data["metric"]
        )
        negative_zero_count = int(np.count_nonzero(
            np.abs(negative_values) < ZERO_TOLERANCE
        ))
        negative_control_ok &= negative_zero_count == 1

        q_zero = bloch_spectrum(data, 0.0)
        expected_q_zero = np.asarray((0.0, 0.0, 12.0 * kappa * vertices))
        q_zero_error = float(np.max(
            np.abs(q_zero - expected_q_zero)
            / np.maximum(1.0, np.abs(expected_q_zero))
        ))
        q_zero_ok &= q_zero_error < 1e-11

        q_first = bloch_spectrum(data, 2.0 * math.pi / vertices)
        positive_bloch = q_first[q_first > ZERO_TOLERANCE][0]
        positive_full = values[values > ZERO_TOLERANCE][0]
        match_error = abs(positive_full / positive_bloch - 1.0)
        full_match_error = max(full_match_error, match_error)
        full_records.append({
            "vertices": vertices,
            "kappa": kappa,
            "zero_count": zero_count,
            "negative_control_zero_count": negative_zero_count,
            "q_zero_eigenvalues": q_zero.tolist(),
            "first_positive_full": float(positive_full),
            "first_positive_bloch": float(positive_bloch),
            "full_bloch_relative_error": float(match_error),
        })

check("the central flux has only endpoint-star support", support_ok)
check("the flux weak operator is exactly odd before adding the even penalty",
      odd_ok)
check("constant scalar and constant one-form are explicit exact zero vectors",
      zero_vectors_ok)
check("every complete finite-kappa central pencil has exactly two zero modes",
      zero_counts_ok)
check("the original no-flux negative control still has exactly one zero",
      negative_control_ok)
check("the q=0 mismatch eigenvalue is exactly 12*kappa/h",
      q_zero_ok)
check("full first-positive spectra match the q=2*pi/N Bloch branch",
      full_match_error < FULL_MATCH_GATE,
      f"maximum relative error={full_match_error:.3e}")

convergence_records = []
all_monotone = True
all_improve = True
all_velocity = True
for kappa in KAPPAS:
    values = []
    errors = []
    velocities = []
    for vertices in LEVELS:
        data = central_pencil(vertices, kappa)
        spectrum = bloch_spectrum(data, 2.0 * math.pi / vertices)
        positive = float(spectrum[spectrum > ZERO_TOLERANCE][0])
        values.append(positive)
        errors.append(abs(positive - 2.0 * math.pi))
        velocities.append(positive / (2.0 * math.pi))
    error_ratios = [errors[index] / errors[index + 1]
                    for index in range(len(errors) - 1)]
    monotone = all(errors[index + 1] < errors[index]
                   for index in range(len(errors) - 1))
    improvement = errors[0] / errors[-1]
    velocity_toward_one = abs(velocities[-1] - 1.0) < abs(
        velocities[0] - 1.0
    )
    all_monotone &= monotone
    all_improve &= improvement > FINAL_IMPROVEMENT_GATE
    all_velocity &= velocity_toward_one
    convergence_records.append({
        "kappa": kappa,
        "levels": list(LEVELS),
        "first_positive_eigenvalues": values,
        "absolute_errors_to_2pi": errors,
        "consecutive_error_ratios": error_ratios,
        "initial_to_final_error_improvement": improvement,
        "dimensionless_velocities": velocities,
        "monotone": monotone,
    })

check("all four fixed-kappa low branches converge monotonically toward 2*pi",
      all_monotone,
      "final errors=" + str([
          record["absolute_errors_to_2pi"][-1]
          for record in convergence_records
      ]))
check("every frozen error improves by more than the preregistered factor 16",
      all_improve,
      "improvements=" + str([
          record["initial_to_final_error_improvement"]
          for record in convergence_records
      ]))
check("every low-mode velocity moves toward the unit continuum value",
      all_velocity,
      "final velocities=" + str([
          record["dimensionless_velocities"][-1]
          for record in convergence_records
      ]))

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "candidate_count": 1,
    "averaging_weights": [0.5, 0.5],
    "full_matrix_records": full_records,
    "convergence_records": convergence_records,
    "verdicts": [
        "DERIVED: unique copy-symmetric central flux preserves both S1 Betti modes at finite kappa",
        "DERIVED NUMERICAL: all frozen fixed-kappa low branches converge to 2*pi",
        "DERIVED: mismatch branch remains at 12*kappa/h",
        "OPEN: canonical 3D Whitney flux, finite-cutoff chirality and causal time",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured central-flux certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for record in convergence_records:
    print(
        f"kappa={record['kappa']}: final="
        f"{record['first_positive_eigenvalues'][-1]:.12g}, "
        f"error={record['absolute_errors_to_2pi'][-1]:.3e}, "
        f"improvement={record['initial_to_final_error_improvement']:.3f}"
    )
raise SystemExit(0 if passed == tests else 1)

