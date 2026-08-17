#!/usr/bin/env python3
"""Stage-B comparison of the frozen blind cellular jet with four tick artifacts.

The analytic prediction was committed first in 76a09ab.  This verifier does
not solve or tune a state; it only evaluates disclosed frozen artifacts.
"""

import hashlib
import json
from pathlib import Path

import mpmath as arb
import sympy as sp


HERE = Path(__file__).resolve().parent
BLIND_INPUT = HERE / "gravity_600cell_cellular_weak_lapse_blind.json"
WEAK_INPUT = HERE / "gravity_600cell_dust_weak_lapse_recurrence.json"
FOURTH_INPUT = HERE / "gravity_600cell_dust_fourth_tick.json"
SCALE_INPUT = HERE / "gravity_600cell_dust_homothetic_mass_conservation.json"
OUTPUT = HERE / "gravity_600cell_cellular_weak_lapse_comparison.json"

BLIND_SHA256 = "6d39e9a4594d9c9ead102f94cf9115d8474132ecce511fe7359826dcc73b9de0"
WEAK_SHA256 = "500be1c4e2d7ec4104b9773bc1cfc71065c9d930607119eb616367d18fa5d8f9"
FOURTH_SHA256 = "4d8d03957675a6f454c1ad05102ffd1711f48c2e5a19f09b2898a60c9f07020d"
SCALE_SHA256 = "72225b1ca17de18f6d77aac43972f4fdca18e24575c8640e8be5e5636316fad0"
BLIND_RESULT_COMMIT = "76a09ab"
PROTOCOL_COMMIT = "71d10b4"
DPS = 100
arb.mp.dps = DPS
ASYMPTOTIC_FLOOR = arb.mpf("1e-40")
RESIDUAL_TOLERANCE = arb.mpf("1e-20")
LAMBDA_TEXTS = ("0.5", "0.25", "0.125")
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mp(value):
    return arb.mpf(str(value))


def text(value, digits=50):
    return arb.nstr(value, digits)


blind = json.loads(BLIND_INPUT.read_text())
weak = json.loads(WEAK_INPUT.read_text())
fourth = json.loads(FOURTH_INPUT.read_text())
scale = json.loads(SCALE_INPUT.read_text())
hashes = {
    "blind": digest(BLIND_INPUT),
    "weak": digest(WEAK_INPUT),
    "fourth": digest(FOURTH_INPUT),
    "scale": digest(SCALE_INPUT),
}
provenance_ok = bool(
    hashes == {
        "blind": BLIND_SHA256,
        "weak": WEAK_SHA256,
        "fourth": FOURTH_SHA256,
        "scale": SCALE_SHA256,
    }
    and BLIND_RESULT_COMMIT == "76a09ab"
    and blind["protocol_commit"] == PROTOCOL_COMMIT
    and blind["outcome"] == "CELLULAR_WEAK_LAPSE_JET_DERIVED"
    and blind["passed"] == blind["tests"] == 8
    and blind["tick_artifacts_parsed"] is False
    and blind["labels"]["comparison_to_committed_ticks"] == "NOT_PERFORMED"
    and weak["outcome"] == "WEAK_LAPSE_QUADRATIC_INTEGER_LAW"
    and weak["passed"] == weak["tests"] == 5
    and fourth["outcome"] == "FOURTH_TICK_WEAK_LAPSE_PREDICTION_CONFIRMED"
    and fourth["passed"] == fourth["tests"] == 5
    and tuple(weak["lambdas"]) == tuple(fourth["lambdas"]) == LAMBDA_TEXTS
)
check(
    "the blind result precedes and freezes every disclosed comparison input",
    provenance_ok,
    f"blind commit={BLIND_RESULT_COMMIT}; hashes={hashes}",
)


# Map the frozen target-blind outputs to the names already stored by the
# numerical recurrence audits.  No numerical target is written here.
ratios = blind["blind_ratios"]
weak_predictions = {
    "u2_over_u1": mp(ratios["u_over_u1"]["2"]),
    "u3_over_u1": mp(ratios["u_over_u1"]["3"]),
    "a2_over_u1": mp(ratios["a_over_u1"]["2"]),
    "a3_over_u1": mp(ratios["a_over_u1"]["3"]),
    "v2_over_v1": mp(ratios["v_over_v1"]["2"]),
    "v3_over_v1": mp(ratios["v_over_v1"]["3"]),
    "r2_over_v1": mp(ratios["r_over_v1"]["2"]),
    "r3_over_v1": mp(ratios["r_over_v1"]["3"]),
    "post1_over_k": mp(ratios["p_post_over_k"]["1"]),
    "post2_over_k": mp(ratios["p_post_over_k"]["2"]),
    "post3_over_k": mp(ratios["p_post_over_k"]["3"]),
}
fourth_predictions = {
    "u4_over_u1": mp(ratios["u_over_u1"]["4"]),
    "a4_over_u1": mp(ratios["a_over_u1"]["4"]),
    "v4_over_v1": mp(ratios["v_over_v1"]["4"]),
    "r4_over_v1": mp(ratios["r_over_v1"]["4"]),
    "post4_over_k": mp(ratios["p_post_over_k"]["4"]),
}

normalized_records = []
normalized_ok = True
for name, prediction in weak_predictions.items():
    record = weak["asymptotic"][name]
    intercept = mp(record["richardson_fine"])
    band = 10*mp(record["richardson_epsilon"])
    error = abs(intercept-prediction)
    accepted = error <= band
    normalized_ok &= accepted
    normalized_records.append({
        "source": "weak", "parity": "even=odd", "name": name,
        "blind_prediction": text(prediction),
        "richardson_intercept": text(intercept),
        "frozen_band": text(band), "absolute_error": text(error),
        "accepted": accepted,
    })

training_band = mp(fourth["training_band_from_artifact"])
for parity in ("even", "odd"):
    for name, prediction in fourth_predictions.items():
        record = fourth["asymptotic"][parity][name]
        intercept = mp(record["richardson_fine"])
        internal_band = 10*mp(record["epsilon4"])
        error = abs(intercept-prediction)
        accepted = error <= internal_band and error <= training_band
        normalized_ok &= accepted
        normalized_records.append({
            "source": "fourth", "parity": parity, "name": name,
            "blind_prediction": text(prediction),
            "richardson_intercept": text(intercept),
            "internal_band": text(internal_band),
            "external_training_band": text(training_band),
            "absolute_error": text(error), "accepted": accepted,
        })
check(
    "all 21 disclosed normalized intercepts contain the blind predictions",
    normalized_ok and len(normalized_records) == 21,
    f"comparisons={len(normalized_records)}; max error="
    f"{text(max(mp(item['absolute_error']) for item in normalized_records), 12)}",
)


# Absolute leading coefficients.  Parse the exact blind expressions, then
# convert e=lambda*tau0/L0 using only the committed carrier constants.
epsilon_symbol = sp.symbols("epsilon_3", positive=True)
locals_exact = {"epsilon_3": epsilon_symbol, "sqrt": sp.sqrt}
epsilon_value = 2*sp.pi-5*sp.acos(sp.Rational(1, 3))


def evaluate_blind_exact(value):
    exact = sp.sympify(value, locals=locals_exact).subs(
        epsilon_symbol, epsilon_value
    )
    return mp(sp.N(exact, DPS))


A1 = evaluate_blind_exact(blind["coefficients"]["A"]["1"])
R1 = evaluate_blind_exact(blind["coefficients"]["R"]["1"])
L0 = mp(scale["frozen_parameters"]["L0"])
rho0 = mp(scale["frozen_parameters"]["rho0"])
tau0 = arb.sqrt(rho0)
e0_square = rho0/L0**2
absolute_predictions = {
    "u1_over_lambda2": A1*e0_square,
    "v1_over_lambda2": R1*e0_square,
}
absolute_records = []
absolute_ok = bool(abs(tau0**2-rho0) < arb.mpf("1e-90"))
for name, prediction in absolute_predictions.items():
    values = [mp(weak["leading_coefficients"][lam][name]) for lam in LAMBDA_TEXTS]
    coarse = (4*values[1]-values[0])/3
    fine = (4*values[2]-values[1])/3
    band = 10*(abs(fine-coarse)+ASYMPTOTIC_FLOOR)
    error = abs(fine-prediction)
    accepted = error <= band
    absolute_ok &= accepted
    absolute_records.append({
        "name": name, "blind_prediction": text(prediction, 70),
        "values": [text(value, 70) for value in values],
        "richardson_coarse": text(coarse, 70),
        "richardson_fine": text(fine, 70),
        "frozen_band": text(band, 70),
        "absolute_error": text(error, 70), "accepted": accepted,
    })
check(
    "the committed tau0/L0 converts both absolute blind coefficients correctly",
    absolute_ok,
    f"e0^2={text(e0_square, 18)}; A1 e0^2={text(absolute_predictions['u1_over_lambda2'], 18)}; "
    f"R1 e0^2={text(absolute_predictions['v1_over_lambda2'], 18)}",
)


# Reconstruct the cellular action and feed every stored state through it.
# This is evaluation only: no Newton step, root selection or staircase model.
L_MINUS, L_PLUS, RHO = sp.symbols("L_minus L_plus rho", positive=True)
EPSILON3 = 2*sp.pi-5*sp.acos(sp.Rational(1, 3))
MU = sp.Rational(90)/sp.pi*EPSILON3
DELTA = L_PLUS-L_MINUS
H = sp.sqrt(RHO+DELTA**2/4)
C = (DELTA**2+2*RHO)/(2*(DELTA**2+3*RHO))
BETA = DELTA/sp.sqrt(8*(DELTA**2+3*RHO))
S = (
    360*(L_MINUS+L_PLUS)*H*(2*sp.pi-5*sp.acos(C))
    + 600*sp.sqrt(3)*(L_MINUS**2-L_PLUS**2)*sp.asinh(BETA)
    - 8*sp.pi*MU*sp.sqrt(RHO)
)
F = sp.lambdify((L_MINUS, L_PLUS, RHO), RHO*sp.diff(S, RHO), "mpmath")
P_MINUS = sp.lambdify(
    (L_MINUS, L_PLUS, RHO), L_MINUS*sp.diff(S, L_MINUS)/2, "mpmath"
)
P_PLUS = sp.lambdify(
    (L_MINUS, L_PLUS, RHO), L_PLUS*sp.diff(S, L_PLUS)/2, "mpmath"
)

rho0_dimensionless = rho0/L0**2
residual_records = []
maximum_f = arb.mpf(0)
maximum_g = arb.mpf(0)
maximum_parity_state_difference = arb.mpf(0)
for lam_text in LAMBDA_TEXTS:
    lam = mp(lam_text)
    parity_states = {}
    for parity in ("even", "odd"):
        parity_states[parity] = (
            [(arb.mpf(0), arb.mpf(0))]
            + [tuple(mp(value) for value in item["state"])
               for item in weak["solves"][lam_text][parity]]
            + [tuple(mp(value) for value in fourth["solves"][lam_text][parity]["state"])]
        )
    for even_state, odd_state in zip(parity_states["even"], parity_states["odd"]):
        maximum_parity_state_difference = max(
            maximum_parity_state_difference,
            *(abs(left-right) for left, right in zip(even_state, odd_state)),
        )
    for parity in ("even", "odd"):
        states = parity_states[parity]
        slabs = []
        for n in range(1, 5):
            l_minus = arb.exp(states[n-1][0])
            l_plus = arb.exp(states[n][0])
            rho_n = lam**2*rho0_dimensionless*arb.exp(states[n][1])
            slabs.append((l_minus, l_plus, rho_n))
        static_slab = (arb.mpf(1), arb.mpf(1), lam**2*rho0_dimensionless)
        for n, slab in enumerate(slabs, start=1):
            previous = static_slab if n == 1 else slabs[n-2]
            f_value = F(*slab)
            g_value = P_PLUS(*previous)+P_MINUS(*slab)
            maximum_f = max(maximum_f, abs(f_value))
            maximum_g = max(maximum_g, abs(g_value))
            residual_records.append({
                "lambda": lam_text, "parity": parity, "n": n,
                "F_abs": text(abs(f_value), 50),
                "G_abs": text(abs(g_value), 50),
            })
residual_ok = bool(
    len(residual_records) == 24
    and maximum_f < RESIDUAL_TOLERANCE
    and maximum_g < RESIDUAL_TOLERANCE
)
check(
    "all 24 committed finite-lapse slabs satisfy the cellular lapse and seam equations",
    residual_ok,
    f"max |F|={text(maximum_f, 12)}; max |G|={text(maximum_g, 12)}; no re-solving",
)

cellular_only_ok = bool(
    maximum_parity_state_difference < arb.mpf("1e-50")
    and blind["action_input_sha256"]
        == "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
)
check(
    "the comparison is parity independent and uses no staircase carrier",
    cellular_only_ok,
    f"max even/odd state difference={text(maximum_parity_state_difference, 8)}",
)

if provenance_ok and normalized_ok and absolute_ok and residual_ok and cellular_only_ok:
    outcome = "CELLULAR_JET_EXPLAINS_FOUR_TICKS"
else:
    outcome = "CELLULAR_JET_DOES_NOT_EXPLAIN_FOUR_TICKS"
check(
    "the disclosed comparison follows the preregistered outcome rule",
    outcome in {
        "CELLULAR_JET_EXPLAINS_FOUR_TICKS",
        "CELLULAR_JET_DOES_NOT_EXPLAIN_FOUR_TICKS",
    },
    f"outcome={outcome}",
)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "blind_result_commit": BLIND_RESULT_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "labels": {
        "analytic_explanation_of_four_ticks": (
            "DERIVED" if outcome == "CELLULAR_JET_EXPLAINS_FOUR_TICKS" else "DERIVED NEGATIVE"
        ),
        "continuum_convergence": "OPEN",
        "physical_clock": "OPEN",
        "external_novelty": "OPEN",
    },
    "normalized_comparisons": normalized_records,
    "absolute_comparisons": absolute_records,
    "scale": {
        "L0": text(L0, 70), "rho0": text(rho0, 70),
        "tau0": text(tau0, 70), "e0_square": text(e0_square, 70),
    },
    "finite_state_residuals": residual_records,
    "maximum_cellular_lapse_residual": text(maximum_f, 60),
    "maximum_cellular_seam_residual": text(maximum_g, 60),
    "maximum_even_odd_state_difference": text(maximum_parity_state_difference, 30),
    "staircase_carrier_used": False,
    "states_resolved_again": False,
    "scope": {
        "ticks": 4,
        "lambdas": list(LAMBDA_TEXTS),
        "spatial_refinement": "NOT TESTED",
        "anisotropic_modes": "NOT TESTED",
        "absolute_clock": "NOT DERIVED",
    },
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print(f"\nSummary: {passed}/{tests} checks passed")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")
if passed != tests or outcome != "CELLULAR_JET_EXPLAINS_FOUR_TICKS":
    raise SystemExit(1)

