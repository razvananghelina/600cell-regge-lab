#!/usr/bin/env python3
"""Audit causal and local-Legendre selectors after two-slab branching."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
COMPOSITION_INPUT = HERE / "gravity_600cell_finite_height_composition.json"
OUTPUT = HERE / "gravity_600cell_finite_height_selector_audit.json"

COMPOSITION_SHA256 = (
    "d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447"
)
SELECTOR_PROTOCOL_COMMIT = "e3c8927"
COMPOSITION_RESULT_COMMIT = "5e93d05"

mp.mp.dps = 100
tests = 0
passed = 0


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


def text(value, digits=60):
    return mp.nstr(value, digits)


composition = json.loads(COMPOSITION_INPUT.read_text())
provenance_ok = bool(
    digest(COMPOSITION_INPUT) == COMPOSITION_SHA256
    and composition["outcome"] == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE"
    and composition["passed"] == composition["tests"] == 10
)
check(
    "the accepted branching artifact and post-result selector protocol are frozen",
    provenance_ok,
    f"protocol={SELECTOR_PROTOCOL_COMMIT}",
)


# Derive the determinant before reading either physical root.  Write u=mu(q),
# up=mu'(q), pp=p'(q).  The state-function identity is polynomial and does
# not require division by q, so q=0 is retained.
h, q, u, up, pp = sp.symbols("h q u up pp", real=True)
pi = sp.pi
c_h = 4 * pi * q * u
c_q = 8 * pi * up + 4 * pi * h * (u + q * up)
p_h = -2 * pi * u
p_q = pp - 2 * pi * h * up
jacobian = sp.expand(c_h * p_q - c_q * p_h)
jacobian_remainder = sp.factor(jacobian - 8 * pi**2 * h * u**2)
expected_remainder = 4 * pi * u * (4 * pi * up + q * pp)
determinant_reduction_ok = sp.expand(
    jacobian_remainder - expected_remainder
) == 0
check(
    "the exact two-equation determinant reduces to the state derivative identity",
    determinant_reduction_ok,
    f"remainder={jacobian_remainder}",
)


# Reproduce 4*pi*mu'(q)+q*p'(q)=0 from the already-proved closed derivative
# formulas.  K is left arbitrary: the cancellation is coefficient-level.
k, denominator = sp.symbols("K denominator", real=True, nonzero=True)
mu_prime = 180 * q * k / (pi * denominator)
p_prime = -720 * k / denominator
state_derivative_identity = sp.simplify(4 * pi * mu_prime + q * p_prime)
exact_jacobian = sp.simplify(
    jacobian.subs({up: mu_prime, pp: p_prime})
)
regularity_identity_ok = bool(
    state_derivative_identity == 0
    and sp.simplify(exact_jacobian - 8 * pi**2 * h * u**2) == 0
)
check(
    "every positive-height root has determinant 8*pi^2*h*mu(q)^2",
    regularity_identity_ok,
)


# Positivity of mu is analytic on the complete real q line.  The cosine in
# epsilon lies in [1/3,1/2), and 1/3 exceeds cos(2*pi/5).
golden_cos = (sp.sqrt(5) - 1) / 4
epsilon_lower_bound_ok = sp.simplify(sp.Rational(1, 3) - golden_cos) > 0
check(
    "mu(q) is strictly positive on the complete real branch",
    epsilon_lower_bound_ok,
    "1/3 > cos(2*pi/5), hence 2*pi-5*acos(c(q)) > 0",
)


# The tetrahedral carrier uses proper timelike struts.  Its central coordinate
# time changes with the radial displacement, so no finite q crosses light
# speed.  This is not the cubical central-height parametrisation.
phi = (1 + sp.sqrt(5)) / 2
beta_squared = sp.factor(phi**2 * q**2 / (1 + phi**2 * q**2))
causal_margin = sp.simplify(1 - beta_squared)
causal_identity_ok = sp.simplify(
    causal_margin - 1 / (1 + phi**2 * q**2)
) == 0
check(
    "the tetrahedral proper-strut parametrisation is subluminal for every finite q",
    causal_identity_ok,
    f"1-beta^2={causal_margin}",
)


def epsilon(value):
    return 2 * mp.pi - 5 * mp.acos(
        (value**2 + 2) / (2 * (value**2 + 3))
    )


def mu(value):
    return 180 * epsilon(value) / (mp.pi * mp.sqrt(value**2 + 4))


record = next(row for row in composition["composition"] if row["v"] == "1.5")
physical_roots = [root for root in record["roots"] if root["physical"]]
same_input_two_roots_ok = bool(
    record["physical_root_count"] == 2
    and len(physical_roots) == 2
    and physical_roots[0]["q2"] != physical_roots[1]["q2"]
)
check(
    "the frozen canonical input has two distinct physical second slabs",
    same_input_two_roots_ok,
    f"physical_roots={len(physical_roots)}",
)


root_audits = []
all_future = True
all_causal = True
all_regular = True
for root in physical_roots:
    q2 = mp.mpf(root["q2"])
    h2 = mp.mpf(root["h2"])
    ratio = mp.mpf(root["scale_ratio"])
    stored_jacobian = mp.mpf(root["jacobian"])
    phi_mp = (1 + mp.sqrt(5)) / 2
    beta2 = phi_mp**2 * q2**2 / (1 + phi_mp**2 * q2**2)
    determinant = 8 * mp.pi**2 * h2 * mu(q2) ** 2
    determinant_error = abs(determinant - stored_jacobian)
    future = h2 > 0
    causal = bool(
        ratio > 0
        and -h2**2 < 0
        and beta2 >= 0
        and beta2 < 1
    )
    regular = bool(
        determinant > 0
        and determinant_error < mp.mpf("1e-55")
    )
    all_future &= future
    all_causal &= causal
    all_regular &= regular
    root_audits.append(
        {
            "q": text(q2),
            "h": text(h2),
            "scale_ratio": text(ratio),
            "beta_squared": text(beta2),
            "causal_margin": text(1 - beta2),
            "strut_squared": text(-h2**2),
            "legendre_determinant": text(determinant),
            "stored_jacobian_error": text(determinant_error, 20),
            "future_oriented": future,
            "same_connected_action_domain": causal,
            "locally_regular": regular,
        }
    )

check(
    "both branches have the same future time orientation",
    all_future,
)
check(
    "both branches remain in the same connected causal cellular-action domain",
    all_causal,
    "rho=h^2>0, L_next/L_current>0 and beta^2<1",
)
check(
    "both branches are locally regular sheets of the discrete Legendre relation",
    all_regular,
    "nonzero local Jacobians do not imply global injectivity",
)


selector_count = sum(
    1
    for row in root_audits
    if row["future_oriented"]
    and row["same_connected_action_domain"]
    and row["locally_regular"]
)
outcome = (
    "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH"
    if provenance_ok
    and determinant_reduction_ok
    and regularity_identity_ok
    and epsilon_lower_bound_ok
    and causal_identity_ok
    and same_input_two_roots_ok
    and all_future
    and all_causal
    and all_regular
    and selector_count == 2
    else "SELECTOR_AUDIT_OPEN"
)
check(
    "the frozen hierarchy rejects causality or regularity as a branch selector",
    outcome == "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH",
    f"surviving_branches={selector_count}/2",
)


artifact = {
    "provenance": {
        "composition_sha256": COMPOSITION_SHA256,
        "composition_result_commit": COMPOSITION_RESULT_COMMIT,
        "selector_protocol_commit": SELECTOR_PROTOCOL_COMMIT,
    },
    "exact_certificates": {
        "state_derivative_identity": "4*pi*mu'(q)+q*p'(q)=0",
        "legendre_determinant": "8*pi^2*h*mu(q)^2",
        "proper_strut_square": "-h^2",
        "beta_squared": "phi^2*q^2/(1+phi^2*q^2)",
        "causal_margin": "1/(1+phi^2*q^2)",
    },
    "frozen_state": "v=3/2",
    "root_audits": root_audits,
    "surviving_branches": selector_count,
    "interpretation": {
        "label": "DERIVED NEGATIVE, SELECTOR-SCOPED",
        "local_legendre_regularity": True,
        "global_legendre_injectivity": False,
        "causal_selector": False,
        "future_orientation_selector": False,
        "additional_selector": "OPEN",
        "fundamental_tick": False,
    },
    "tests": tests,
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)
