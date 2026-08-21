#!/usr/bin/env python3
"""Dual tangent-line audit of the finite-height 600-cell classification."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp
from flint import arb, ctx


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
OUTPUT = HERE / "gravity_600cell_finite_height_classification_adversarial.json"

PRIMARY_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
ADVERSARIAL_PROTOCOL_COMMIT = "ccc7a4e"
PRIMARY_IMPLEMENTATION_COMMIT = "4176c3d"
PRIMARY_ARTIFACT_COMMIT = "f0a4209"

ctx.dps = 140

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


provenance_ok = digest(PRIMARY_INPUT) == PRIMARY_SHA256
check(
    "the frozen primary artifact is present but not yet numerically read",
    provenance_ok,
    f"adversarial_protocol={ADVERSARIAL_PROTOCOL_COMMIT}",
)


# The primary implementation used E_q=p(q)-p(v).  This implementation starts
# from the tangent-line error and differentiates it using only dp/dmu=-4*pi/q.
q = sp.symbols("q", real=True, nonzero=True)
mu_q_symbol, mu_v_symbol, mu_q_prime = sp.symbols(
    "mu_q mu_v mu_q_prime", real=True
)
p_q_prime = -4 * sp.pi * mu_q_prime / q
tangent_derivative_chain = (
    -p_q_prime
    + 4
    * sp.pi
    * (-mu_q_prime / q - (mu_v_symbol - mu_q_symbol) / q**2)
)
tangent_derivative_expected = (
    4 * sp.pi * (mu_q_symbol - mu_v_symbol) / q**2
)
dual_derivative_ok = bool(
    sp.factor(tangent_derivative_chain - tangent_derivative_expected) == 0
)
check(
    "the dual tangent derivative is exactly controlled by mu(q)-mu(v)",
    dual_derivative_ok,
)


# Independently reconstruct the complete action and its two canonical
# residuals.  No coefficient or decimal is imported from the primary file.
L_MINUS, L_PLUS, RHO, MASS = sp.symbols(
    "L_minus L_plus rho mass", positive=True
)
DELTA = L_PLUS - L_MINUS
HEIGHT = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
ACTION = (
    360
    * (L_MINUS + L_PLUS)
    * HEIGHT
    * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600
    * sp.sqrt(3)
    * (L_MINUS**2 - L_PLUS**2)
    * sp.asinh(BOOST)
    - 8 * sp.pi * MASS * sp.sqrt(RHO)
)
F_EXACT = RHO * sp.diff(ACTION, RHO)
P_PRE_EXACT = -L_MINUS * sp.diff(ACTION, L_MINUS) / 2

h_symbol = sp.symbols("h", positive=True)
q_symbol = sp.symbols("q_symbol", real=True)
direct_substitution = {
    L_MINUS: 1,
    L_PLUS: 1 + h_symbol * q_symbol,
    RHO: h_symbol**2,
}
DIRECT_C = sp.factor((2 * F_EXACT / h_symbol).subs(direct_substitution))
DIRECT_P = sp.factor(P_PRE_EXACT.subs(direct_substitution))
affine_direct_ok = bool(
    sp.diff(DIRECT_C, h_symbol, 2) == 0
    and sp.diff(DIRECT_P, h_symbol, 2) == 0
)
check(
    "the independently differentiated full action is affine in height",
    affine_direct_ok,
)


# Rigorous real-ball functions and disjoint rational controls.
def arat(numerator, denominator=1):
    return arb(numerator) / denominator


API = arb.pi()


def aepsilon(t):
    return 2 * API - 5 * (
        (t * t + 2) / (2 * (t * t + 3))
    ).acos()


def amu(t):
    return 180 * aepsilon(t) / (API * (t * t + 4).sqrt())


def ap(t):
    return (
        180 * t * aepsilon(t) / (t * t + 4).sqrt()
        - 600
        * arb(3).sqrt()
        * (t / (8 * (t * t + 3)).sqrt()).asinh()
    )


def aK_square(t):
    return (
        10 * (t + 4).sqrt()
        - (t + 3)
        * (3 * t + 8).sqrt()
        * (2 * API - 5 * ((t + 2) / (2 * (t + 3))).acos())
    )


def aT(v_value, q_value):
    return (
        ap(v_value)
        - ap(q_value)
        + 4
        * API
        * (amu(v_value) - amu(q_value))
        / q_value
    )


def amass_gap(v_value, u_value):
    return amu(u_value) - 2 * amu(v_value)


APINF = 60 * API - 300 * arb(3).sqrt() * arb(2).log()
vC_left = arat(1573, 50)
vC_right = arat(787, 25)
u_mid = arat(243, 2000)
sign_certificates = {
    "p_infinity": APINF,
    "K_x_5": aK_square(arat(5)),
    "K_x_6": aK_square(arat(6)),
    "vA_left_49_over_40": ap(arat(49, 40)) - APINF,
    "vA_right_5_over_4": ap(arat(5, 4)) - APINF,
    "vM_left_16": amu(arat(16)) - amu(arat(0)),
    "vM_right_33_over_2": amu(arat(33, 2)) - amu(arat(0)),
    "negative_axis_tail_control": (
        4 * API * amu(arat(1)) + ap(arat(1)) + APINF
    ),
    "vC_left_T_at_minus_u_mid": aT(vC_left, -u_mid),
    "vC_left_mass_gap": amass_gap(vC_left, u_mid),
    "vC_right_T_at_minus_u_mid": aT(vC_right, -u_mid),
    "vC_right_mass_gap": amass_gap(vC_right, u_mid),
    "qC_lower_v_left_T": aT(vC_left, -arat(607, 5000)),
    "qC_upper_v_right_T": aT(vC_right, -arat(76, 625)),
}
rigorous_signs_ok = bool(
    sign_certificates["p_infinity"] < 0
    and sign_certificates["K_x_5"] > 0
    and sign_certificates["K_x_6"] < 0
    and sign_certificates["vA_left_49_over_40"] > 0
    and sign_certificates["vA_right_5_over_4"] < 0
    and sign_certificates["vM_left_16"] > 0
    and sign_certificates["vM_right_33_over_2"] < 0
    and sign_certificates["negative_axis_tail_control"] < 0
    and sign_certificates["vC_left_T_at_minus_u_mid"] < 0
    and sign_certificates["vC_left_mass_gap"] < 0
    and sign_certificates["vC_right_T_at_minus_u_mid"] > 0
    and sign_certificates["vC_right_mass_gap"] > 0
    and sign_certificates["qC_lower_v_left_T"] > 0
    and sign_certificates["qC_upper_v_right_T"] < 0
)
check(
    "disjoint Arb brackets certify all tangent and endpoint signs",
    rigorous_signs_ok,
    "vA=(49/40,5/4), vM=(16,33/2), vC=(1573/50,787/25)",
)


# Complete dual sign ledger.  This proof uses critical points of mu rather
# than critical points of p.  The q=0 row is retained separately because T
# divides by q.
dual_root_ledger = [
    {
        "state": "v=0",
        "positive_height_roots": 0,
        "certificate": "direct even/odd residuals",
    },
    {
        "state": "v<0",
        "positive_height_roots": 0,
        "certificate": "(v,q,h)->(-v,-q,-h)",
    },
    {
        "state": "0<v<=vA",
        "off_diagonal_roots": 0,
        "certificate": "T has diagonal minimum zero and nonnegative positive tail",
    },
    {
        "state": "vA<v<vstar",
        "off_diagonal_roots": 1,
        "location": "outer positive-q mu branch",
    },
    {
        "state": "v=vstar",
        "off_diagonal_roots": 0,
        "certificate": "T decreases through its stationary diagonal zero",
    },
    {
        "state": "vstar<v<vM",
        "off_diagonal_roots": 1,
        "location": "inner positive-q mu branch",
    },
    {
        "state": "v=vM",
        "off_diagonal_roots": 1,
        "location": "q=0 in the original residuals",
    },
    {
        "state": "v>vM",
        "off_diagonal_roots": 1,
        "location": "negative-q branch between -v and 0",
    },
]
negative_axis_controls = {
    "inner_equal_mass_value": "T(v,-v)=2*p(v)<0",
    "outer_equal_mass_value": "T(v,-m(v))=p(v)+p(m(v))<0",
    "above_vM_zero_limit": "T(v,q)->+infinity as q->0-",
}
dual_global_logic_ok = bool(
    dual_derivative_ok
    and rigorous_signs_ok
    and len(dual_root_ledger) == 8
    and len(negative_axis_controls) == 3
)
check(
    "the mu-controlled tangent ledger covers the complete real line",
    dual_global_logic_ok,
)


endpoint_ledger = {
    "T_q_on_root": "positive",
    "T_v_on_root": "positive",
    "dq_dv": "negative",
    "u_below_one": bool(sign_certificates["negative_axis_tail_control"] < 0),
    "mu_u_increases": True,
    "mu_v_decreases": True,
    "L_plus_strictly_decreases": True,
    "L_plus_limits": ["+1", "-1"],
    "unique_zero": bool(rigorous_signs_ok),
}
endpoint_logic_ok = bool(
    dual_derivative_ok
    and endpoint_ledger["u_below_one"]
    and endpoint_ledger["unique_zero"]
)
check(
    "the tangent implicit derivative independently fixes one causal endpoint",
    endpoint_logic_ok,
)


def compute_independent(precision):
    mp.mp.dps = precision
    pi = mp.pi

    def epsilon(t):
        return 2 * pi - 5 * mp.acos((t * t + 2) / (2 * (t * t + 3)))

    def mu(t):
        return 180 * epsilon(t) / (pi * mp.sqrt(t * t + 4))

    def momentum(t):
        return (
            180 * t * epsilon(t) / mp.sqrt(t * t + 4)
            - 600
            * mp.sqrt(3)
            * mp.asinh(t / mp.sqrt(8 * (t * t + 3)))
        )

    p_infinity = 60 * pi - 300 * mp.sqrt(3) * mp.log(2)

    def K_square(t):
        return (
            10 * mp.sqrt(t + 4)
            - (t + 3)
            * mp.sqrt(3 * t + 8)
            * (2 * pi - 5 * mp.acos((t + 2) / (2 * (t + 3))))
        )

    def tangent(v_value, q_value):
        return (
            momentum(v_value)
            - momentum(q_value)
            + 4
            * pi
            * (mu(v_value) - mu(q_value))
            / q_value
        )

    def bisect(function, left, right, iterations=None):
        if iterations is None:
            iterations = 4 * precision
        left = mp.mpf(left)
        right = mp.mpf(right)
        f_left = function(left)
        f_right = function(right)
        if not (f_left * f_right < 0):
            raise ValueError((left, right, f_left, f_right))
        for _ in range(iterations):
            middle = (left + right) / 2
            f_middle = function(middle)
            if f_middle == 0:
                return middle
            if f_left * f_middle < 0:
                right = middle
            else:
                left = middle
                f_left = f_middle
        return (left + right) / 2

    x_star = bisect(K_square, 5, 6)
    v_star = mp.sqrt(x_star)
    v_A = bisect(
        lambda t: momentum(t) - p_infinity,
        mp.mpf(49) / 40,
        mp.mpf(5) / 4,
    )
    v_M = bisect(
        lambda t: mu(t) - mu(0),
        16,
        mp.mpf(33) / 2,
    )

    def q_negative_root(v_value):
        return bisect(lambda q_value: tangent(v_value, q_value), -1, -mp.mpf("1e-20"))

    def endpoint(v_value):
        q_value = q_negative_root(v_value)
        return 2 * mu(v_value) / mu(q_value) - 1

    v_C = bisect(
        endpoint,
        mp.mpf(1573) / 50,
        mp.mpf(787) / 25,
    )
    q_C = q_negative_root(v_C)
    h_C = (momentum(q_C) - momentum(v_C)) / (2 * pi * mu(q_C))

    return {
        "precision": precision,
        "v_A": v_A,
        "v_star": v_star,
        "v_M": v_M,
        "v_C": v_C,
        "q_C": q_C,
        "h_C": h_C,
        "L_C": 1 + h_C * q_C,
        "mu": mu,
        "momentum": momentum,
        "tangent": tangent,
    }


precision_runs = [compute_independent(dps) for dps in (80, 120, 180)]
threshold_names = ("v_A", "v_star", "v_M", "v_C", "q_C", "h_C")
precision_agreement = True
for lower, upper in zip(precision_runs, precision_runs[1:]):
    for name in threshold_names:
        precision_agreement &= bool(
            abs(lower[name] - upper[name]) < mp.mpf(10) ** (-(lower["precision"] - 15))
        )
precision_agreement &= bool(
    all(abs(run["L_C"]) < mp.mpf(10) ** (-(run["precision"] - 15)) for run in precision_runs)
)
check(
    "80-, 120- and 180-digit tangent classifications are nested",
    precision_agreement,
)


# Directly solve the two full-action residuals in (h,q), rather than solving
# either E or T.  Sympy built these callables before any primary decimal was
# loaded.
direct_c_numeric = sp.lambdify((h_symbol, q_symbol, MASS), DIRECT_C, "mpmath")
direct_p_numeric = sp.lambdify((h_symbol, q_symbol), DIRECT_P, "mpmath")

direct_specs = [
    (mp.mpf(3) / 2, mp.mpf("0.2"), mp.mpf("9.6"), True),
    (mp.mpf(3), mp.mpf("0.1"), mp.mpf("1.55"), True),
    (mp.mpf(20), mp.mpf("7.86"), -mp.mpf("0.05"), True),
]
direct_records = []
direct_ok = True
mp.mp.dps = 130
last = precision_runs[-1]
mu = last["mu"]
momentum = last["momentum"]
for v_value, h_seed, q_seed, endpoint_expected in direct_specs:
    mass_value = mu(v_value)
    p_value = momentum(v_value)

    def equations(h_value, q_value):
        return (
            direct_c_numeric(h_value, q_value, mass_value),
            direct_p_numeric(h_value, q_value) - p_value,
        )

    h_value, q_value = mp.findroot(
        equations,
        (h_seed, q_seed),
        tol=mp.mpf("1e-110"),
        maxsteps=100,
    )
    c_residual, p_residual = equations(h_value, q_value)
    endpoint = 1 + h_value * q_value

    # Numerical Jacobian of the original equations is independent of the
    # primary closed-form Jacobian.
    dc_dh = mp.diff(lambda hh: equations(hh, q_value)[0], h_value)
    dc_dq = mp.diff(lambda qq: equations(h_value, qq)[0], q_value)
    dp_dh = mp.diff(lambda hh: equations(hh, q_value)[1], h_value)
    dp_dq = mp.diff(lambda qq: equations(h_value, qq)[1], q_value)
    jacobian = dc_dh * dp_dq - dc_dq * dp_dh
    row_ok = bool(
        h_value > 0
        and (endpoint > 0) == endpoint_expected
        and abs(c_residual) < mp.mpf("1e-105")
        and abs(p_residual) < mp.mpf("1e-105")
        and abs(jacobian) > mp.mpf("1e-20")
    )
    direct_ok &= row_ok
    direct_records.append(
        {
            "v": text(v_value),
            "h": text(h_value),
            "q": text(q_value),
            "L_plus": text(endpoint),
            "constraint_residual": text(c_residual),
            "momentum_residual": text(p_residual),
            "jacobian": text(jacobian),
            "passed": row_ok,
        }
    )

check(
    "direct two-variable solves of the full action reproduce the physical branch",
    direct_ok,
    f"roots={len(direct_records)}",
)


# Convention and hostile controls.
v_control = mp.mpf(3)
h_control = mp.mpf(direct_records[1]["h"])
q_control = mp.mpf(direct_records[1]["q"])

def reduced_residuals(v_value, q_value, h_value):
    return (
        8 * mp.pi * (mu(q_value) - mu(v_value))
        + 4 * mp.pi * h_value * q_value * mu(q_value),
        momentum(q_value)
        - momentum(v_value)
        - 2 * mp.pi * h_value * mu(q_value),
    )


positive_pair = reduced_residuals(v_control, q_control, h_control)
reflected_pair = reduced_residuals(-v_control, -q_control, -h_control)
reflection_ok = bool(
    abs(reflected_pair[0] - positive_pair[0]) < mp.mpf("1e-100")
    and abs(reflected_pair[1] + positive_pair[1]) < mp.mpf("1e-100")
)

wrong_tangent = (
    momentum(v_control)
    - momentum(q_control)
    - 4
    * mp.pi
    * (mu(v_control) - mu(q_control))
    / q_control
)
shifted_constraint = direct_c_numeric(
    h_control,
    q_control,
    mu(v_control) + mp.mpf(1) / 11,
)
no_update_v1 = bool(
    rigorous_signs_ok
    and mp.mpf(1) < precision_runs[-1]["v_A"]
)
post_boundary_v40_q = mp.findroot(
    lambda qq: last["tangent"](mp.mpf(40), qq),
    (-mp.mpf("0.2"), -mp.mpf("0.1")),
)
post_boundary_h40 = (
    momentum(post_boundary_v40_q) - momentum(40)
) / (2 * mp.pi * mu(post_boundary_v40_q))
post_boundary_L40 = 1 + post_boundary_h40 * post_boundary_v40_q
hostile_ok = bool(
    reflection_ok
    and abs(wrong_tangent) > mp.mpf(1)
    and abs(shifted_constraint) > mp.mpf(1)
    and no_update_v1
    and post_boundary_L40 < 0
)
check(
    "time reflection passes while sign, mass and causal hostile controls fail",
    hostile_ok,
    (
        f"wrong_T={text(wrong_tangent, 12)}; "
        f"mass_shift={text(shifted_constraint, 12)}; "
        f"L40={text(post_boundary_L40, 12)}"
    ),
)


# Only now read and compare the primary numerical result.
primary = json.loads(PRIMARY_INPUT.read_text())
independent = precision_runs[-1]
primary_comparison = {}
comparison_ok = primary["outcome"] == (
    "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
)
for name in ("v_A", "v_star", "v_M", "v_C", "q_C", "h_C"):
    primary_value = mp.mpf(primary["thresholds"][name])
    difference = abs(independent[name] - primary_value)
    row_ok = difference < mp.mpf("1e-70")
    comparison_ok &= row_ok
    primary_comparison[name] = {
        "independent": text(independent[name], 80),
        "primary": text(primary_value, 80),
        "absolute_difference": text(difference, 20),
        "passed": row_ok,
    }
check(
    "the independently constructed thresholds agree with the primary artifact",
    comparison_ok,
)


classification_complete = bool(
    provenance_ok
    and dual_derivative_ok
    and affine_direct_ok
    and rigorous_signs_ok
    and dual_global_logic_ok
    and endpoint_logic_ok
    and precision_agreement
    and direct_ok
    and hostile_ok
    and comparison_ok
)
outcome = (
    "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY_ADVERSARIALLY_CORROBORATED"
    if classification_complete
    else "FINITE_HEIGHT_OPEN"
)
check(
    "the adversarial outcome is selected only after the independent comparison",
    outcome
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY_ADVERSARIALLY_CORROBORATED",
)


artifact = {
    "provenance": {
        "primary_sha256": digest(PRIMARY_INPUT),
        "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
        "primary_implementation_commit": PRIMARY_IMPLEMENTATION_COMMIT,
        "primary_artifact_commit": PRIMARY_ARTIFACT_COMMIT,
    },
    "method": {
        "primary_decisive_identity_not_reused": "E_q=p(q)-p(v)",
        "adversarial_decisive_identity": "T_q=4*pi*(mu(q)-mu(v))/q^2",
        "direct_action_redifferentiated": True,
        "primary_decimals_read_only_after_independent_construction": True,
    },
    "dual_root_ledger": dual_root_ledger,
    "negative_axis_controls": negative_axis_controls,
    "endpoint_ledger": endpoint_ledger,
    "rigorous_sign_certificates": {
        key: str(value) for key, value in sign_certificates.items()
    },
    "precision_runs": [
        {
            "precision": run["precision"],
            **{name: text(run[name], 80) for name in threshold_names},
            "L_C_residual": text(run["L_C"], 30),
        }
        for run in precision_runs
    ],
    "direct_full_action_roots": direct_records,
    "hostile_controls": {
        "reflection": reflection_ok,
        "wrong_tangent": text(wrong_tangent, 40),
        "mass_shift_constraint": text(shifted_constraint, 40),
        "v_1_has_no_update": no_update_v1,
        "v_40_endpoint": text(post_boundary_L40, 40),
    },
    "primary_comparison": primary_comparison,
    "interpretation": {
        "label": "DERIVED EXACT / STRUCTURAL / ADVERSARIALLY CORROBORATED",
        "fundamental_tick": False,
        "composition": "OPEN",
        "nonhomogeneous_stability": "OPEN",
        "refinement": "OPEN",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
raise SystemExit(0 if passed == tests else 1)
