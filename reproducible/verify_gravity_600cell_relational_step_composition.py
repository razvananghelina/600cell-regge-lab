#!/usr/bin/env python3
"""Exact same-state half-step branch test for the cellular dust map.

Prior-art commit: 4f4b5c2.
Protocol commit: a025951.
"""

import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ACTION_INPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
OUTPUT = HERE / "gravity_600cell_relational_step_composition.json"
ACTION_SHA256 = (
    "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
)
PRIOR_ART_COMMIT = "4f4b5c2"
PROTOCOL_COMMIT = "a025951"
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


action_artifact = json.loads(ACTION_INPUT.read_text())
provenance_ok = bool(
    digest(ACTION_INPUT) == ACTION_SHA256
    and action_artifact["outcome"] == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and action_artifact["passed"] == action_artifact["tests"] == 16
    and PRIOR_ART_COMMIT == "4f4b5c2"
    and PROTOCOL_COMMIT == "a025951"
)
check("the exact cellular action input and preregistration are frozen", provenance_ok)


# Reconstruct the complete action rather than importing a weak-lapse equation.
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
F_RHO = RHO * sp.diff(ACTION, RHO)
P_MINUS = L_MINUS * sp.diff(ACTION, L_MINUS) / 2
P_PLUS = L_PLUS * sp.diff(ACTION, L_PLUS) / 2

TAU = sp.symbols("tau", positive=True)
static_action = sp.simplify(
    ACTION.subs({L_MINUS: 1, L_PLUS: 1, RHO: TAU**2})
)
static_p_minus = sp.simplify(
    P_MINUS.subs({L_MINUS: 1, L_PLUS: 1, RHO: TAU**2})
)
static_p_plus = sp.simplify(
    P_PLUS.subs({L_MINUS: 1, L_PLUS: 1, RHO: TAU**2})
)
static_ok = bool(
    static_action == 0
    and exact_zero(static_p_minus - 180 * EPSILON * TAU)
    and exact_zero(static_p_plus - 180 * EPSILON * TAU)
)
check("the reconstructed action has the exact all-lapse static family", static_ok)


# Constant-force calibration: this detects state misalignment independently of
# every Regge formula.
q0, p0, force, step = sp.symbols("q0 p0 force step", real=True)


def constant_force_flow(q, p, duration):
    return (
        q + duration * p + force * duration**2 / 2,
        p + force * duration,
    )


coarse_control = constant_force_flow(q0, p0, step)
fine_mid = constant_force_flow(q0, p0, step / 2)
fine_control = constant_force_flow(*fine_mid, step / 2)
constant_force_ok = all(
    sp.expand(left - right) == 0
    for left, right in zip(coarse_control, fine_control)
)
check("the exact constant-force calibration composes from the same state", constant_force_ok)


# Build S/e as a series in x=e^2 for a general nominal duration factor s.
# Differentiating this series gives F/e and P/e without importing the earlier
# all-index coefficient equations.
x = sp.symbols("x", positive=True)
a_minus, a_plus, relative_rho, duration_factor = sp.symbols(
    "a_minus a_plus relative_rho duration_factor", real=True
)
lm = sp.exp(a_minus * x)
lp = sp.exp(a_plus * x)
rho_bar = duration_factor**2 * sp.exp(relative_rho * x)
delta = lp - lm
delta_over_x = sp.series(delta / x, x, 0, 3).removeO()
delta_square_over_x = sp.series(delta**2 / x, x, 0, 3).removeO()
height_bar = sp.sqrt(rho_bar + delta_square_over_x / 4)
cosine_bar = (
    delta_square_over_x + 2 * rho_bar
) / (2 * (delta_square_over_x + 3 * rho_bar))
boost_bar = delta_over_x / sp.sqrt(
    8 * (delta_square_over_x + 3 * rho_bar)
)
asinh_over_sqrt_x = boost_bar - x * boost_bar**3 / 6 + 3 * x**2 * boost_bar**5 / 40
s_bar = sp.series(
    360 * (lm + lp) * height_bar * (2 * sp.pi - 5 * sp.acos(cosine_bar))
    + 600 * sp.sqrt(3) * (lm**2 - lp**2) * asinh_over_sqrt_x
    - 8 * sp.pi * MASS * sp.sqrt(rho_bar),
    x,
    0,
    3,
).removeO().expand()
f_scaled = sp.series(sp.diff(s_bar, relative_rho) / x, x, 0, 2).removeO().expand()
p_minus_scaled = sp.series(sp.diff(s_bar, a_minus) / (2 * x), x, 0, 2).removeO().expand()

E = sp.symbols("epsilon", positive=True, real=True)
alpha_symbol = sp.symbols("alpha", positive=True, real=True)


def compact(expression):
    with_alpha = sp.expand(expression).xreplace({ALPHA: alpha_symbol})
    return sp.factor(
        with_alpha.subs(alpha_symbol, (2 * sp.pi - E) / 5),
        extension=sp.sqrt(2),
    )


A = sp.symbols("A", real=True)
D = 5 * sp.sqrt(2) / 3 - E
base_substitution = {a_minus: 0, a_plus: A, relative_rho: 0}

coarse_f = compact(
    f_scaled.coeff(x, 1).subs({**base_substitution, duration_factor: 1})
)
coarse_p = compact(
    p_minus_scaled.coeff(x, 0).subs(
        {**base_substitution, duration_factor: 1}
    )
    + 180 * E
)
coarse_expected_f = 45 * A * (D * A + 4 * E)
coarse_expected_p = 90 * (D * A + 4 * E)
coarse_factorization_ok = bool(
    exact_zero(coarse_f - coarse_expected_f)
    and exact_zero(coarse_p - coarse_expected_p)
)
coarse_root = -4 * E / D
known_coarse_root = -12 * E / (5 * sp.sqrt(2) - 3 * E)
coarse_branch_ok = bool(
    coarse_factorization_ok
    and exact_zero(coarse_expected_f.subs(A, coarse_root))
    and exact_zero(coarse_expected_p.subs(A, coarse_root))
    and exact_zero(coarse_root - known_coarse_root)
)
check(
    "the coarse equations recover the blind-derived contracting coefficient",
    coarse_branch_ok,
    f"A_coarse={sp.factor(coarse_root)}",
)


fine_f = compact(
    f_scaled.coeff(x, 1).subs(
        {**base_substitution, duration_factor: sp.Rational(1, 2)}
    )
)
fine_p_same = compact(
    p_minus_scaled.coeff(x, 0).subs(
        {**base_substitution, duration_factor: sp.Rational(1, 2)}
    )
    + 180 * E
)
fine_expected_f = 90 * A * (D * A + E)
fine_expected_p_same = 180 * (D * A + sp.Rational(3, 2) * E)
fine_factorization_ok = bool(
    exact_zero(fine_f - fine_expected_f)
    and exact_zero(fine_p_same - fine_expected_p_same)
)
check(
    "the half-lapse lapse and same-state momentum equations factor exactly",
    fine_factorization_ok,
    f"F={fine_expected_f}; P={fine_expected_p_same}",
)

lapse_roots = (sp.Integer(0), -E / D)
same_state_momentum_root = -sp.Rational(3, 2) * E / D
no_common_root = bool(
    D != 0
    and E != 0
    and all(
        not exact_zero(root - same_state_momentum_root)
        for root in lapse_roots
    )
)
resultant = sp.factor(sp.resultant(fine_expected_f, fine_expected_p_same, A))
resultant_nonzero = bool(resultant != 0)
check(
    "the same-state half-step leading equations have zero common roots",
    no_common_root and resultant_nonzero,
    f"lapse roots={lapse_roots}; momentum root={same_state_momentum_root}; resultant={resultant}",
)


# Hostile control: the old lambda family halves the initial momentum together
# with the lapse.  That altered state restores the nonzero lapse root exactly.
fine_p_changed_state = compact(
    p_minus_scaled.coeff(x, 0).subs(
        {**base_substitution, duration_factor: sp.Rational(1, 2)}
    )
    + 90 * E
)
fine_expected_p_changed = 180 * (D * A + E)
changed_state_root = -E / D
changed_state_restores_branch = bool(
    exact_zero(fine_p_changed_state - fine_expected_p_changed)
    and exact_zero(fine_expected_f.subs(A, changed_state_root))
    and exact_zero(fine_expected_p_changed.subs(A, changed_state_root))
)
changed_state_difference = sp.factor((180 * E - 90 * E) / (180 * E))
check(
    "halving the momentum restores a branch but fails the same-state gate",
    changed_state_restores_branch and changed_state_difference == sp.Rational(1, 2),
    f"restored A={changed_state_root}; relative state change={changed_state_difference}",
)


N_FINE = 0 if no_common_root and resultant_nonzero else None
if not (provenance_ok and static_ok and constant_force_ok and coarse_branch_ok):
    outcome = "RELATIONAL_STEP_COMPOSITION_CONTROL_FAILED"
elif N_FINE == 0:
    outcome = "SAME_STATE_HALF_STEP_BRANCH_ABSENT"
else:
    outcome = "TEMPORAL_REFINEMENT_NUMERICALLY_UNRESOLVED"

check(
    "the frozen outcome hierarchy assigns the exact branch verdict",
    outcome == "SAME_STATE_HALF_STEP_BRANCH_ABSENT",
    outcome,
)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "action_input_sha256": digest(ACTION_INPUT),
    "same_state": {
        "L0": "1",
        "coarse_and_fine_initial_momentum_over_e": "180*epsilon",
        "mass": "(90/pi)*epsilon",
    },
    "coarse": {
        "lapse_leading": str(coarse_expected_f),
        "momentum_leading": str(coarse_expected_p),
        "common_contracting_root": str(sp.factor(coarse_root)),
    },
    "fine_first_slab": {
        "nominal_duration_ratio": "1/2",
        "lapse_leading": str(fine_expected_f),
        "same_state_momentum_leading": str(fine_expected_p_same),
        "lapse_roots": [str(sp.factor(value)) for value in lapse_roots],
        "same_state_momentum_root": str(sp.factor(same_state_momentum_root)),
        "resultant": str(resultant),
        "N_fine": N_FINE,
    },
    "hostile_changed_state": {
        "incoming_momentum_over_e": "90*epsilon",
        "relative_initial_momentum_change": str(changed_state_difference),
        "restored_root": str(sp.factor(changed_state_root)),
        "branch_restored": changed_state_restores_branch,
    },
    "nonlinear_stage_executed": False,
    "nonlinear_stage_reason": (
        "The preregistered exact leading branch count is zero; a nonlinear "
        "seed search cannot restore a missing weak-lapse branch."
    ),
    "labels": {
        "same_state_half_step_absence": "DERIVED_NEGATIVE_SCOPED",
        "old_lambda_comparison": "INVALID_CHANGED_INITIAL_STATE",
        "fundamental_tick": "NOT_DERIVED",
        "absolute_tick": "DERIVED_NEGATIVE_UNDER_SCALE_COVARIANCE_HYPOTHESES",
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

