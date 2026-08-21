#!/usr/bin/env python3
"""Exact-action-first adversarial audit of the half-step obstruction."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_relational_step_composition.json"
OUTPUT = HERE / "gravity_600cell_relational_step_composition_adversarial.json"
PRIMARY_SHA256 = (
    "b7c4d46b96eb0fbc266237390462474c132a0852a1e39e42e7bbf009ca15a2e6"
)
PRIMARY_IMPLEMENTATION_COMMIT = "55cb3a2"
ADVERSARIAL_PROTOCOL_COMMIT = "177232b"
CORRECTION_PROTOCOL_COMMIT = "1b0a024"
SECOND_CORRECTION_PROTOCOL_COMMIT = "e548944"
mp.mp.dps = 100
E_VALUES = (mp.mpf(1) / 100, mp.mpf(1) / 200, mp.mpf(1) / 400)
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_zero(expression):
    return sp.cancel(sp.expand(expression), extension=sp.sqrt(2)) == 0


def text(value, digits=40):
    return mp.nstr(value, digits)


primary = json.loads(PRIMARY_INPUT.read_text())
provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and primary["outcome"] == "SAME_STATE_HALF_STEP_BRANCH_ABSENT"
    and primary["passed"] == primary["tests"] == 8
    and PRIMARY_IMPLEMENTATION_COMMIT == "55cb3a2"
    and ADVERSARIAL_PROTOCOL_COMMIT == "177232b"
    and CORRECTION_PROTOCOL_COMMIT == "1b0a024"
    and SECOND_CORRECTION_PROTOCOL_COMMIT == "e548944"
)
check("the primary result and adversarial protocol are frozen", provenance_ok)


# Mechanically independent order: differentiate the unexpanded action first.
L_MINUS, L_PLUS, RHO = sp.symbols("L_minus L_plus rho", positive=True)
ALPHA = sp.acos(sp.Rational(1, 3))
EPSILON = 2 * sp.pi - 5 * ALPHA
MASS = sp.Rational(90) / sp.pi * EPSILON
DELTA = L_PLUS - L_MINUS
HEIGHT = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
ACTION = (
    360 * (L_MINUS + L_PLUS) * HEIGHT * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600 * sp.sqrt(3) * (L_MINUS**2 - L_PLUS**2) * sp.asinh(BOOST)
    - 8 * sp.pi * MASS * sp.sqrt(RHO)
)
F_EXACT = sp.factor(RHO * sp.diff(ACTION, RHO))
P_EXACT = sp.factor(L_MINUS * sp.diff(ACTION, L_MINUS) / 2)

TAU = sp.symbols("tau", positive=True)
static_ok = bool(
    sp.simplify(ACTION.subs({L_MINUS: 1, L_PLUS: 1, RHO: TAU**2})) == 0
    and sp.simplify(
        P_EXACT.subs({L_MINUS: 1, L_PLUS: 1, RHO: TAU**2})
        - 180 * EPSILON * TAU
    ) == 0
)
check("the exact-action-first route reproduces the static control", static_ok)


e = sp.symbols("e", positive=True)
A = sp.symbols("A", real=True)

# exp(A e^2)=1+A e^2+O(e^4).  The rescaled limits below depend only on
# (L_plus-1)/e^2, so the tangent path is exactly equivalent for these limits.
exp_tangent = sp.limit((sp.exp(A * e**2) - 1) / e**2, e, 0, dir="+")
tangent_ok = exp_tangent == A
check("the linear proxy has the exact exponential half-step tangent", tangent_ok)

tangent_substitution = {
    L_MINUS: 1,
    L_PLUS: 1 + A * e**2,
    RHO: e**2 / 4,
}
raw_f_limit = sp.factor(
    sp.limit(F_EXACT.subs(tangent_substitution) / e**3, e, 0, dir="+")
)
raw_p_limit = sp.factor(
    sp.limit(
        (P_EXACT.subs(tangent_substitution) + 180 * EPSILON * e) / e,
        e,
        0,
        dir="+",
    )
)

# SymPy expresses the real asinh branch through this principal complex log.
# Its exact polar form follows from
# (2 sqrt(2)+i)/3 = exp(i(pi/2-acos(1/3))).
complex_log = sp.log(2 * sp.sqrt(2) + sp.I)
branch_log_value = sp.log(3) + sp.I * (sp.pi / 2 - ALPHA)
unit_complex_identity = sp.simplify(
    sp.expand_complex(
        sp.exp(sp.I * (sp.pi / 2 - ALPHA))
        - (2 * sp.sqrt(2) + sp.I) / 3
    )
)
branch_range_ok = bool(
    sp.N(sp.pi / 2 - ALPHA, 80) > 0
    and sp.N(sp.pi / 2 - ALPHA, 80) < sp.pi / 2
)
branch_identity_ok = bool(unit_complex_identity == 0 and branch_range_ok)
reduced_p_limit = sp.factor(raw_p_limit.subs(complex_log, branch_log_value))

D_PHYSICAL = 5 * sp.sqrt(2) / 3 - EPSILON
expected_f_limit = 90 * A * (D_PHYSICAL * A + EPSILON)
expected_p_limit = 180 * (D_PHYSICAL * A + sp.Rational(3, 2) * EPSILON)
limits_ok = bool(
    exact_zero(raw_f_limit - expected_f_limit)
    and exact_zero(reduced_p_limit - expected_p_limit)
)
check(
    "exact-action-first limits reproduce both independent leading equations",
    branch_identity_ok and limits_ok,
    f"F={sp.factor(raw_f_limit)}; P={sp.factor(reduced_p_limit)}",
)


lapse_root = -EPSILON / D_PHYSICAL
momentum_root = -sp.Rational(3, 2) * EPSILON / D_PHYSICAL
resultant = sp.factor(sp.resultant(expected_f_limit, expected_p_limit, A))
physical_signs_ok = bool(
    sp.N(EPSILON, 100) > 0
    and sp.N(D_PHYSICAL, 100) > 0
    and sp.N(resultant, 100) != 0
)
root_obstructions = (
    sp.simplify(expected_p_limit.subs(A, lapse_root)),
    sp.simplify(expected_f_limit.subs(A, momentum_root)),
)
no_common_root = bool(
    physical_signs_ok
    and all(sp.simplify(value) != 0 for value in root_obstructions)
)
check(
    "physical epsilon gives a nonzero resultant and two resolved obstructions",
    no_common_root,
    f"resultant={resultant}; obstructions={root_obstructions}",
)


changed_state_limit = sp.factor(
    expected_p_limit - 90 * EPSILON
)
changed_state_ok = bool(
    sp.simplify(expected_f_limit.subs(A, lapse_root)) == 0
    and sp.simplify(changed_state_limit.subs(A, lapse_root)) == 0
)
check("the changed-state hostile control restores the lapse branch", changed_state_ok)


# Direct 100-decimal evaluations use exp(A e^2), not the tangent proxy.
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO), F_EXACT, "mpmath")
p_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO), P_EXACT, "mpmath")
epsilon_mp = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
d_mp = 5 * mp.sqrt(2) / 3 - epsilon_mp
a_lapse_mp = -epsilon_mp / d_mp
a_momentum_mp = -mp.mpf(3) * epsilon_mp / (2 * d_mp)
expected_p_obstruction_mp = 90 * epsilon_mp
expected_f_obstruction_mp = mp.mpf(135) * epsilon_mp**2 / (2 * d_mp)

numeric_records = []
for e_value in E_VALUES:
    rho_value = e_value**2 / 4
    lp_lapse = mp.exp(a_lapse_mp * e_value**2)
    lp_momentum = mp.exp(a_momentum_mp * e_value**2)
    numeric_records.append({
        "e": e_value,
        "lapse_root_scaled_lapse": f_numeric(1, lp_lapse, rho_value) / e_value**3,
        "lapse_root_scaled_same_state_momentum": (
            p_numeric(1, lp_lapse, rho_value) + 180 * epsilon_mp * e_value
        ) / e_value,
        "lapse_root_scaled_changed_state_momentum": (
            p_numeric(1, lp_lapse, rho_value) + 90 * epsilon_mp * e_value
        ) / e_value,
        "momentum_root_scaled_lapse": (
            f_numeric(1, lp_momentum, rho_value) / e_value**3
        ),
        "momentum_root_scaled_momentum": (
            p_numeric(1, lp_momentum, rho_value) + 180 * epsilon_mp * e_value
        ) / e_value,
    })


def halving_orders(values):
    return tuple(
        mp.log(abs(values[index] / values[index + 1])) / mp.log(2)
        for index in range(2)
    )


zero_sequences = {
    "lapse_at_lapse_root": [
        record["lapse_root_scaled_lapse"] for record in numeric_records
    ],
    "changed_momentum_at_lapse_root": [
        record["lapse_root_scaled_changed_state_momentum"]
        for record in numeric_records
    ],
    "momentum_at_momentum_root": [
        record["momentum_root_scaled_momentum"] for record in numeric_records
    ],
}
zero_orders = {name: halving_orders(values) for name, values in zero_sequences.items()}
zero_convergence_ok = all(
    abs(values[2]) < abs(values[1]) < abs(values[0])
    and all(mp.mpf("1.8") <= order <= mp.mpf("2.2") for order in zero_orders[name])
    for name, values in zero_sequences.items()
)
check(
    "all three direct defining residuals vanish with quadratic order",
    zero_convergence_ok,
    "; ".join(
        f"{name} orders=({text(orders[0], 10)},{text(orders[1], 10)})"
        for name, orders in zero_orders.items()
    ),
)

p_obstructions = [
    record["lapse_root_scaled_same_state_momentum"] for record in numeric_records
]
f_obstructions = [record["momentum_root_scaled_lapse"] for record in numeric_records]
p_error = abs(p_obstructions[-1] - expected_p_obstruction_mp)
f_error = abs(f_obstructions[-1] - expected_f_obstruction_mp)
p_drift = abs(p_obstructions[-1] - p_obstructions[-2])
f_drift = abs(f_obstructions[-1] - f_obstructions[-2])
p_signal_ratio = abs(p_obstructions[-1]) / max(p_drift, mp.mpf("1e-99"))
f_signal_ratio = abs(f_obstructions[-1]) / max(f_drift, mp.mpf("1e-99"))
obstructions_ok = bool(
    p_error < p_drift
    and f_error < f_drift
    and p_signal_ratio > 100
    and f_signal_ratio > 100
)
check(
    "both direct nonzero obstructions dominate their precision drifts",
    obstructions_ok,
    f"momentum signal/drift={text(p_signal_ratio, 12)}; lapse signal/drift={text(f_signal_ratio, 12)}",
)


all_controls = bool(
    provenance_ok
    and static_ok
    and tangent_ok
    and branch_identity_ok
    and limits_ok
    and no_common_root
    and changed_state_ok
    and zero_convergence_ok
    and obstructions_ok
)
outcome = (
    "SAME_STATE_HALF_STEP_ABSENCE_ADVERSARIALLY_CORROBORATED"
    if all_controls
    else "HALF_STEP_ADVERSARIAL_DISAGREEMENT"
)
check(
    "the adversarial hierarchy assigns the corroborated verdict",
    outcome == "SAME_STATE_HALF_STEP_ABSENCE_ADVERSARIALLY_CORROBORATED",
    outcome,
)

artifact = {
    "primary_input_sha256": digest(PRIMARY_INPUT),
    "primary_implementation_commit": PRIMARY_IMPLEMENTATION_COMMIT,
    "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
    "correction_protocol_commit": CORRECTION_PROTOCOL_COMMIT,
    "second_correction_protocol_commit": SECOND_CORRECTION_PROTOCOL_COMMIT,
    "method": "differentiate_exact_action_then_take_rescaled_limits",
    "exact_control_booleans": {
        "branch_identity_ok": branch_identity_ok,
        "limits_ok": limits_ok,
    },
    "exact_limits": {
        "lapse": str(sp.factor(expected_f_limit)),
        "same_state_momentum": str(sp.factor(expected_p_limit)),
        "changed_state_momentum": str(sp.factor(changed_state_limit)),
        "lapse_root": str(sp.factor(lapse_root)),
        "same_state_momentum_root": str(sp.factor(momentum_root)),
        "resultant": str(resultant),
        "root_obstructions": [str(sp.factor(value)) for value in root_obstructions],
    },
    "numeric": {
        "e_values": [text(value, 30) for value in E_VALUES],
        "records": [
            {key: text(value, 50) for key, value in record.items()}
            for record in numeric_records
        ],
        "zero_orders": {
            name: [text(value, 30) for value in orders]
            for name, orders in zero_orders.items()
        },
        "expected_momentum_obstruction": text(expected_p_obstruction_mp, 50),
        "expected_lapse_obstruction": text(expected_f_obstruction_mp, 50),
        "momentum_signal_to_drift": text(p_signal_ratio, 30),
        "lapse_signal_to_drift": text(f_signal_ratio, 30),
    },
    "labels": {
        "same_state_half_step_absence": "DERIVED_NEGATIVE_ADVERSARIALLY_CORROBORATED_SCOPED",
        "fundamental_tick": "NOT_DERIVED",
        "generic_nonzero_velocity_composition": "OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nOutcome: {outcome}")
print(f"Checks: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
