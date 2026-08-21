#!/usr/bin/env python3
"""Direct-2400-simplex adversarial audit of tick scale covariance.

Protocol commit: fd1f8a8.
Primary result commit: 0ac0aba.

This verifier does not call the primary 100-decimal orbit action.  It first
constructs and evaluates every direct full-simplex state with raw squared-
length derivatives, and only then reads the primary artifact.
"""

import ast
from collections import Counter
import contextlib
import io
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_gravity_global_boundary_legendre.py"
PRIMARY = HERE / "gravity_600cell_tick_scale_covariance.json"
OUTPUT = HERE / "gravity_600cell_tick_scale_covariance_adversarial.json"
PROTOCOL_COMMIT = "fd1f8a8"
PRIMARY_RESULT_COMMIT = "0ac0aba"
ALPHAS = (3.0 / 5.0, 7.0 / 4.0)
HOSTILE_MINIMUM = 1e-8
SUPPORT_MINIMUM = 1e-12
EPSILON = float(np.finfo(float).eps)
OPERATION_FACTOR = 128.0
SIMPLEX_COUNT = 2400
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


def relative_error(left, right):
    return abs(left - right) / max(1.0, abs(left), abs(right))


def load_direct_core():
    """Load definitions through full_evaluation, before upstream main work."""
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "CONTROLS"
                for target in node.targets
            )
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("direct full-simplex evaluator cutoff not found")
    namespace = {
        "__file__": str(SOURCE),
        "__name__": "tick_scale_covariance_direct_core",
    }
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(SOURCE), "exec"), namespace)
    return namespace


print("Adversarial direct-simplex tick scale-covariance audit", flush=True)
core = load_direct_core()
models = core["models"]
full_evaluation = core["full_evaluation"]
provenance_ok = bool(
    PROTOCOL_COMMIT == "fd1f8a8"
    and PRIMARY_RESULT_COMMIT == "0ac0aba"
    and core["tests"] == core["passed"] == 7
    and set(models) == {"even", "odd"}
)
check(
    "the independent direct-2400-simplex core and frozen provenance are intact",
    provenance_ok,
)


# Reconstruct the physical constants independently in binary64.
M_STAR = 10.0
ZETA = (math.pi**2 * math.sqrt(2.0) / 50.0) ** (1.0 / 3.0)
R0 = 4.0 * M_STAR / (3.0 * math.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON_3 = 2.0 * math.pi - 5.0 * math.acos(1.0 / 3.0)
MASS = (90.0 / math.pi) * EPSILON_3 * L0
TAU = 0.0102
RHO = TAU**2
SLANT_SQUARE = L0_SQUARE - RHO


def perturbed(values, modulus, center):
    values = np.asarray(values, dtype=float)
    factors = np.array([
        math.exp(1e-6 * ((index % modulus) - center))
        for index in range(len(values))
    ])
    return values * factors


base_old = perturbed(np.full(30, L0_SQUARE), 7, 3)
base_variables = perturbed(
    np.r_[np.full(30, SLANT_SQUARE), np.full(5, RHO), np.full(30, L0_SQUARE)],
    5,
    2,
)

# Construct all states before the primary artifact is read.
scaled_states = {
    alpha: {
        "r": alpha**2,
        "old": (alpha**2) * base_old,
        "variables": (alpha**2) * base_variables,
    }
    for alpha in ALPHAS
}


def direct_total(model, old_values, variables, mass):
    action, gradient, old_gradient, branch = full_evaluation(
        model, variables, old_values
    )
    dust_action = -(8.0 * math.pi * mass / 5.0) * float(
        np.sum(np.sqrt(variables[30:35]))
    )
    dust_gradient = np.zeros(65, dtype=complex)
    dust_gradient[30:35] = (
        -(4.0 * math.pi * mass / 5.0) / np.sqrt(variables[30:35])
    )
    total_action = action + dust_action
    total_gradient = gradient + dust_gradient
    imaginary = max(
        abs(total_action.imag),
        float(np.max(np.abs(total_gradient.imag))),
        float(np.max(np.abs(old_gradient.imag))),
    )
    branch_ok = bool(
        branch["negative_counts"] == Counter({1: SIMPLEX_COUNT})
        and branch["minimum_gram"] > 1e-8
        and branch["minimum_argument"] > 1e-6
        and imaginary < 2e-7
    )
    return {
        "action": total_action,
        "gradient": total_gradient,
        "old_gradient": old_gradient,
        "branch": branch,
        "imaginary": imaginary,
        "branch_ok": branch_ok,
    }


def propagated_bound(*evaluations):
    minimum_gram = min(item["branch"]["minimum_gram"] for item in evaluations)
    conditioning = max(1.0, 1.0 / minimum_gram)
    return OPERATION_FACTOR * EPSILON * SIMPLEX_COUNT * conditioning


records = {}
constructed = {}
covariance_ok = True
hostile_ok = True
support_ok = True
envelope_ok = True
for parity, model in models.items():
    base = direct_total(model, base_old, base_variables, MASS)
    internal_support = int(np.sum(np.abs(base["gradient"][:35]) > SUPPORT_MINIMUM))
    boundary_support = int(
        np.sum(np.abs(base["gradient"][35:]) > SUPPORT_MINIMUM)
        + np.sum(np.abs(base["old_gradient"]) > SUPPORT_MINIMUM)
    )
    support_ok &= bool(
        base["branch_ok"] and internal_support >= 20 and boundary_support >= 20
    )
    parity_records = {
        "base": {
            "action_real": float(base["action"].real),
            "action_imaginary": float(base["action"].imag),
            "internal_raw_gradient_support": internal_support,
            "boundary_raw_gradient_support": boundary_support,
            "minimum_gram": float(base["branch"]["minimum_gram"]),
            "minimum_argument": float(base["branch"]["minimum_argument"]),
            "maximum_imaginary": float(base["imaginary"]),
        },
        "scales": {},
    }
    constructed[parity] = {"base": base, "scales": {}}
    for alpha, state in scaled_states.items():
        r = state["r"]
        scaled = direct_total(
            model, state["old"], state["variables"], alpha * MASS
        )
        fixed_mass = direct_total(
            model, state["old"], state["variables"], MASS
        )
        bound = propagated_bound(base, scaled, fixed_mass)
        action_error = relative_error(scaled["action"], r * base["action"])
        raw_gradient_error = max(
            relative_error(left, right)
            for left, right in zip(scaled["gradient"], base["gradient"])
        )
        old_gradient_error = max(
            relative_error(left, right)
            for left, right in zip(scaled["old_gradient"], base["old_gradient"])
        )
        fixed_action_error = relative_error(
            fixed_mass["action"], r * base["action"]
        )
        fixed_pole_error = max(
            relative_error(fixed_mass["gradient"][index], base["gradient"][index])
            for index in range(30, 35)
        )
        one_envelope_ok = bool(bound < 1e-3)
        one_covariance_ok = bool(
            scaled["branch_ok"]
            and action_error <= bound
            and raw_gradient_error <= bound
            and old_gradient_error <= bound
        )
        one_hostile_ok = bool(
            fixed_mass["branch_ok"]
            and fixed_action_error > max(HOSTILE_MINIMUM, 100.0 * bound)
            and fixed_pole_error > max(HOSTILE_MINIMUM, 100.0 * bound)
        )
        envelope_ok &= one_envelope_ok
        covariance_ok &= one_covariance_ok
        hostile_ok &= one_hostile_ok
        constructed[parity]["scales"][alpha] = {
            "scaled": scaled,
            "fixed_mass": fixed_mass,
        }
        parity_records["scales"][format(alpha, ".17g")] = {
            "r": r,
            "propagated_binary64_bound": bound,
            "action_error": action_error,
            "maximum_raw_internal_and_final_gradient_error": raw_gradient_error,
            "maximum_raw_old_gradient_error": old_gradient_error,
            "fixed_mass_action_error": fixed_action_error,
            "fixed_mass_maximum_pole_error": fixed_pole_error,
            "branch_ok": bool(scaled["branch_ok"]),
            "covariance_pass": one_covariance_ok,
            "hostile_control_pass": one_hostile_ok,
        }
        check(
            f"{parity}, alpha={alpha:.8g}: direct action scales and all 95 raw derivatives are invariant",
            one_covariance_ok and one_envelope_ok,
            f"bound={bound:.3e}, action={action_error:.3e}, grad={raw_gradient_error:.3e}, old={old_gradient_error:.3e}",
        )
        check(
            f"{parity}, alpha={alpha:.8g}: the independent fixed-mass control rejects covariance",
            one_hostile_ok,
            f"action={fixed_action_error:.3e}, pole={fixed_pole_error:.3e}",
        )
    records[parity] = parity_records

check(
    "the direct off-shell state is branch-valid and the propagated envelopes are nonvacuous",
    support_ok and envelope_ok,
    "; ".join(
        f"{parity}: internal={item['base']['internal_raw_gradient_support']}, boundary={item['base']['boundary_raw_gradient_support']}"
        for parity, item in records.items()
    ),
)


# Only now may the primary result enter the process.
primary = json.loads(PRIMARY.read_text())
primary_ok = bool(
    primary.get("outcome") == "TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED"
    and primary.get("tests") == primary.get("passed") == 12
    and tuple(float(value) for value in primary.get("alphas", ())) == ALPHAS
)
primary_action_errors = {}
for parity in models:
    primary_text = primary["parities"][parity]["base"]["action"].replace(" ", "")
    primary_action = complex(primary_text)
    primary_action_errors[parity] = relative_error(
        constructed[parity]["base"]["action"], primary_action
    )
implementations_agree = bool(
    primary_ok
    and covariance_ok
    and max(primary_action_errors.values()) < 2e-8
)
check(
    "the independently completed direct construction agrees with the frozen primary artifact",
    implementations_agree,
    "; ".join(
        f"{parity} base-action error={error:.3e}"
        for parity, error in primary_action_errors.items()
    ),
)

if not (
    provenance_ok and support_ok and envelope_ok and hostile_ok and primary_ok
):
    outcome = "TICK_SCALE_COVARIANCE_ADVERSARIAL_CONTROL_FAILED"
elif not (covariance_ok and implementations_agree):
    outcome = "TICK_SCALE_COVARIANCE_IMPLEMENTATIONS_DISAGREE"
else:
    outcome = "ABSOLUTE_CLASSICAL_TICK_NO_GO_CORROBORATED"

check(
    "the adversarial mechanical outcome corroborates the conditional absolute-tick no-go",
    outcome == "ABSOLUTE_CLASSICAL_TICK_NO_GO_CORROBORATED",
    outcome,
)

artifact = {
    "title": "Adversarial direct-simplex audit of classical tick scale covariance",
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_result_commit": PRIMARY_RESULT_COMMIT,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "independence": {
        "geometry_action": "direct 2400-simplex binary64 evaluator",
        "derivative_coordinates": "raw squared lengths",
        "primary_geometry_action_not_called": True,
        "primary_artifact_read_after_all_direct_evaluations": True,
        "shared_upstream": (
            "the same Regge action definition, Lorentzian angle convention, "
            "carrier, and physical dust normalization"
        ),
    },
    "propagated_bound": (
        "128*binary64_epsilon*2400*max(1,1/minimum_Gram_modulus)"
    ),
    "alphas": list(ALPHAS),
    "parities": records,
    "primary_action_errors": primary_action_errors,
    "interpretation": {
        "derived_exact_conditional": (
            "under q->alpha^2 q and mass->alpha mass, the scale-free classical "
            "Regge-dust equations and canonical data occur in global scale families"
        ),
        "derived_negative": (
            "the stated classical action cannot select an absolute nonzero tick"
        ),
        "open": [
            "dimensionless tick ratio selection under refinement",
            "relational dust time",
            "a separately justified scale-breaking ingredient",
            "quantum dimensional transmutation",
        ],
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

print(f"RESULT: {passed}/{tests}", flush=True)
print(f"OUTCOME: {outcome}", flush=True)
if passed != tests:
    raise SystemExit(1)

