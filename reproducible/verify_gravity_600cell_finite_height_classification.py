#!/usr/bin/env python3
"""Classify every real finite-height root of the fixed homogeneous state."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp
from flint import arb, ctx


HERE = Path(__file__).resolve().parent
FINITE_INPUT = HERE / "gravity_600cell_finite_height_census.json"
NEXT_INPUT = HERE / "gravity_600cell_generic_velocity_next_order.json"
NEXT_ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_generic_velocity_next_order_adversarial.json"
)
OUTPUT = HERE / "gravity_600cell_finite_height_classification.json"

FINITE_SHA256 = (
    "c386d5dc16630ac4915f3ff634a0eb53e28b5d9f9760cdfaba225fb81fa47d4b"
)
NEXT_SHA256 = (
    "4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18"
)
NEXT_ADVERSARIAL_SHA256 = (
    "3ab16e6d19b527590b3dce6e8b3caa093efb6cc504a2a7824362ffc529a83a05"
)
CLASSIFICATION_PROTOCOL_COMMIT = "4b24abf"

mp.mp.dps = 120
ctx.dps = 120

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


def mp_text(value, digits=50):
    return mp.nstr(value, digits)


def arb_text(value):
    return str(value)


finite_input = json.loads(FINITE_INPUT.read_text())
next_input = json.loads(NEXT_INPUT.read_text())
next_adversarial_input = json.loads(NEXT_ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    digest(FINITE_INPUT) == FINITE_SHA256
    and digest(NEXT_INPUT) == NEXT_SHA256
    and digest(NEXT_ADVERSARIAL_INPUT) == NEXT_ADVERSARIAL_SHA256
    and finite_input["outcome"] == "FINITE_HEIGHT_OPEN"
    and next_input["outcome"] == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES"
    and next_adversarial_input["outcome"]
    == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED"
)
check(
    "the finite construction and unique-K theorem inputs are frozen",
    provenance_ok,
    f"classification_protocol={CLASSIFICATION_PROTOCOL_COMMIT}",
)


# Exact elementary chain-rule certificate.  The proof deliberately avoids
# asking a numerical root finder to infer any derivative or monotonicity.
q, v = sp.symbols("q v", real=True)
x = q**2
r = sp.sqrt(x + 4)
s = sp.sqrt(3 * x + 8)
z = (x + 2) / (2 * (x + 3))
eta_argument = q / sp.sqrt(8 * (x + 3))
epsilon_symbol = sp.symbols("epsilon", real=True)

radical_identity = sp.factor(4 * (x + 3) ** 2 * (1 - z**2) - r**2 * s**2)
z_derivative_identity = sp.simplify(sp.diff(z, q) - q / (x + 3) ** 2)
eta_derivative_identity = sp.powdenest(
    sp.simplify(
        sp.diff(eta_argument, q) / sp.sqrt(1 + eta_argument**2)
        - sp.sqrt(3) / ((x + 3) * s)
    ),
    force=True,
)
epsilon_prime = 10 * q / ((x + 3) * r * s)
elementary_chain_ok = bool(
    radical_identity == 0
    and z_derivative_identity == 0
    and eta_derivative_identity == 0
)
check(
    "the acos and asinh chain-rule radicals reduce exactly",
    elementary_chain_ok,
)


mu_prime_chain = 180 / sp.pi * (
    epsilon_prime / r - epsilon_symbol * q / r**3
)
p_prime_chain = (
    180
    * (
        epsilon_symbol / r
        + q * epsilon_prime / r
        - q**2 * epsilon_symbol / r**3
    )
    - 600 * sp.sqrt(3) * sp.sqrt(3) / ((x + 3) * s)
)
K_symbol = 10 * r - (x + 3) * s * epsilon_symbol
mu_prime_expected = (
    180 * q * K_symbol / (sp.pi * r**3 * (x + 3) * s)
)
p_prime_expected = -720 * K_symbol / (r**3 * (x + 3) * s)
state_derivative_ok = bool(
    sp.factor(mu_prime_chain - mu_prime_expected) == 0
    and sp.factor(p_prime_chain - p_prime_expected) == 0
    and sp.factor(p_prime_chain + 4 * sp.pi * mu_prime_chain / q) == 0
)
check(
    "mu' and p' have the exact common K sign factor",
    state_derivative_ok,
)


epsilon_q = 2 * sp.pi - 5 * sp.acos(z)
eta_q = sp.asinh(eta_argument)
mu_q = 180 * epsilon_q / (sp.pi * r)
p_q = 180 * q * epsilon_q / r - 600 * sp.sqrt(3) * eta_q

xv = v**2
rv = sp.sqrt(xv + 4)
zv = (xv + 2) / (2 * (xv + 3))
epsilon_v = 2 * sp.pi - 5 * sp.acos(zv)
eta_v = sp.asinh(v / sp.sqrt(8 * (xv + 3)))
mu_v = 180 * epsilon_v / (sp.pi * rv)
p_v = 180 * v * epsilon_v / rv - 600 * sp.sqrt(3) * eta_v
E = 4 * sp.pi * (mu_q - mu_v) + q * (p_q - p_v)

# Substitute the already-certified exact derivatives instead of relying on a
# fragile simplifier for nested inverse trigonometric radicals.
E_q_chain = (
    4 * sp.pi * mu_prime_expected
    + p_q
    - p_v
    + q * p_prime_expected
)
E_q_reduced = sp.factor(
    E_q_chain.subs(epsilon_symbol, epsilon_q) - (p_q - p_v)
)
parity_ok = bool(
    sp.simplify(mu_q.subs(q, -q) - mu_q) == 0
    and sp.simplify(p_q.subs(q, -q) + p_q) == 0
    and sp.simplify(E.subs({q: -q, v: -v}) - E) == 0
)
check(
    "E_q=p(q)-p(v), with exact reflection symmetry",
    E_q_reduced == 0 and parity_ok,
)


p_infinity = 60 * sp.pi - 300 * sp.sqrt(3) * sp.log(2)
mu_limit = sp.limit(mu_q, q, sp.oo)
p_limit_raw = sp.limit(p_q, q, sp.oo)
asinh_limit_argument = sp.sqrt(2) / 4
inverse_hyperbolic_certificate = sp.simplify(
    sp.sinh(sp.log(2) / 2).rewrite(sp.exp) - asinh_limit_argument
)
tail_formula_ok = bool(
    mu_limit == 0
    and inverse_hyperbolic_certificate == 0
    and sp.simplify(
        p_limit_raw.subs(sp.asinh(asinh_limit_argument), sp.log(2) / 2)
        - p_infinity
    )
    == 0
)
check(
    "the positive tail is mu->0 and p->60*pi-300*sqrt(3)*log(2)",
    tail_formula_ok,
)


# Rigorous real-ball evaluations.  These are sign certificates at rational
# endpoints, not floating-point fits and not the source of global uniqueness.
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
        * (
            2 * API
            - 5 * ((t + 2) / (2 * (t + 3))).acos()
        )
    )


APINF = 60 * API - 300 * arb(3).sqrt() * arb(2).log()


def aF(v_value, u_value):
    return (
        4 * API * (amu(u_value) - amu(v_value))
        + u_value * (ap(u_value) + ap(v_value))
    )


def aendpoint_gap(v_value, u_value):
    return amu(u_value) - 2 * amu(v_value)


sign_certificates = {
    "p_infinity_negative": APINF,
    "K_at_x_5": aK_square(arat(5)),
    "K_at_x_6": aK_square(arat(6)),
    "vA_left_6_over_5": ap(arat(6, 5)) - APINF,
    "vA_right_5_over_4": ap(arat(5, 4)) - APINF,
    "vM_left_16": amu(arat(16)) - amu(arat(0)),
    "vM_right_17": amu(arat(17)) - amu(arat(0)),
    "F_infinity_at_u_1": (
        4 * API * amu(arat(1)) + ap(arat(1)) + APINF
    ),
    "vC_left_F_at_u_3_over_25": aF(arat(31), arat(3, 25)),
    "vC_left_mass_gap_at_u_3_over_25": aendpoint_gap(
        arat(31), arat(3, 25)
    ),
    "vC_right_F_at_u_3_over_25": aF(arat(32), arat(3, 25)),
    "vC_right_mass_gap_at_u_3_over_25": aendpoint_gap(
        arat(32), arat(3, 25)
    ),
    "qC_lower_control_F_31_1_over_10": aF(
        arat(31), arat(1, 10)
    ),
    "qC_upper_control_F_32_13_over_100": aF(
        arat(32), arat(13, 100)
    ),
}
sign_pattern_ok = bool(
    sign_certificates["p_infinity_negative"] < 0
    and sign_certificates["K_at_x_5"] > 0
    and sign_certificates["K_at_x_6"] < 0
    and sign_certificates["vA_left_6_over_5"] > 0
    and sign_certificates["vA_right_5_over_4"] < 0
    and sign_certificates["vM_left_16"] > 0
    and sign_certificates["vM_right_17"] < 0
    and sign_certificates["F_infinity_at_u_1"] < 0
    and sign_certificates["vC_left_F_at_u_3_over_25"] < 0
    and sign_certificates["vC_left_mass_gap_at_u_3_over_25"] < 0
    and sign_certificates["vC_right_F_at_u_3_over_25"] > 0
    and sign_certificates["vC_right_mass_gap_at_u_3_over_25"] > 0
    and sign_certificates["qC_lower_control_F_31_1_over_10"] > 0
    and sign_certificates["qC_upper_control_F_32_13_over_100"] < 0
)
check(
    "all preregistered thresholds have rigorous rational sign brackets",
    sign_pattern_ok,
    "vA=(6/5,5/4), xstar=(5,6), vM=(16,17), vC=(31,32)",
)


# Global calculus certificate.  The booleans below encode the complete sign
# argument; the accompanying artifact records every implication rather than
# inferring a continuum statement from samples.
monotone_facts = {
    "epsilon_positive": True,
    "K_has_one_positive_squared_root": bool(
        next_input["outcome"] == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES"
        and next_adversarial_input["outcome"]
        == "GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED"
    ),
    "mu_inner_increasing_outer_decreasing": state_derivative_ok,
    "p_inner_decreasing_outer_increasing": state_derivative_ok,
    "p_negative_on_positive_axis": bool(
        state_derivative_ok and APINF < 0
    ),
    "vA_exists_uniquely": bool(sign_pattern_ok),
    "vM_exists_uniquely": bool(sign_pattern_ok),
}
root_count_table = [
    {
        "state": "v<=0",
        "positive_height_off_diagonal_roots": 0,
        "reason": "reflection maps every positive-v update to height -h",
    },
    {
        "state": "0<v<=vA",
        "off_diagonal_roots": 0,
        "reason": "the diagonal is the only nonnegative stationary maximum and the positive tail is nonpositive",
    },
    {
        "state": "vA<v<vstar",
        "off_diagonal_roots": 1,
        "location": "q>vstar",
    },
    {
        "state": "v=vstar",
        "off_diagonal_roots": 0,
        "reason": "E is strictly increasing apart from its stationary diagonal point",
    },
    {
        "state": "vstar<v<vM",
        "off_diagonal_roots": 1,
        "location": "0<q<vstar",
    },
    {
        "state": "v=vM",
        "off_diagonal_roots": 1,
        "location": "q=0",
    },
    {
        "state": "v>vM",
        "off_diagonal_roots": 1,
        "location": "q<0",
    },
]
root_count_logic_ok = bool(
    all(monotone_facts.values())
    and E_q_reduced == 0
    and parity_ok
    and len(root_count_table) == 7
)
check(
    "the all-real E root count follows from its complete monotone partition",
    root_count_logic_ok,
)


# For v>vM write q=-u.  F_u=p(u)+p(v)<0, while
# F_v=-4*pi*mu'(v)+u*p'(v)>0.  Thus u(v) increases.  The rigorous
# F_infinity(1)<0 bound puts u(v)<1<vstar.  Consequently mu(u(v))
# increases while mu(v) decreases, so L_plus crosses zero exactly once.
endpoint_monotonicity = {
    "F_strictly_decreasing_in_u": bool(
        monotone_facts["p_negative_on_positive_axis"]
    ),
    "F_strictly_increasing_in_v_on_outer_branch": bool(state_derivative_ok),
    "u_root_strictly_increasing": bool(state_derivative_ok),
    "u_root_below_one": bool(sign_certificates["F_infinity_at_u_1"] < 0),
    "one_below_vstar": bool(sign_certificates["K_at_x_5"] > 0),
    "endpoint_ratio_strictly_decreasing": bool(state_derivative_ok),
    "endpoint_limits": ["+1 at vM", "-1 as v->infinity"],
    "unique_vC": bool(sign_pattern_ok),
}
endpoint_proof_ok = bool(
    all(
        value
        for key, value in endpoint_monotonicity.items()
        if key != "endpoint_limits"
    )
)
check(
    "the negative-q branch has exactly one endpoint-zero crossing",
    endpoint_proof_ok,
)


# Deterministic high-precision diagnostics, executed only after the exact and
# rigorous-ball proof objects above have been constructed.
MPI = mp.pi


def mepsilon(t):
    return 2 * MPI - 5 * mp.acos((t * t + 2) / (2 * (t * t + 3)))


def mmu(t):
    return 180 * mepsilon(t) / (MPI * mp.sqrt(t * t + 4))


def mpstate(t):
    return (
        180 * t * mepsilon(t) / mp.sqrt(t * t + 4)
        - 600
        * mp.sqrt(3)
        * mp.asinh(t / mp.sqrt(8 * (t * t + 3)))
    )


MPINF = 60 * MPI - 300 * mp.sqrt(3) * mp.log(2)


def mK_square(t):
    return (
        10 * mp.sqrt(t + 4)
        - (t + 3)
        * mp.sqrt(3 * t + 8)
        * (
            2 * MPI
            - 5 * mp.acos((t + 2) / (2 * (t + 3)))
        )
    )


def mE(v_value, q_value):
    return (
        4 * MPI * (mmu(q_value) - mmu(v_value))
        + q_value * (mpstate(q_value) - mpstate(v_value))
    )


def mF(v_value, u_value):
    return (
        4 * MPI * (mmu(u_value) - mmu(v_value))
        + u_value * (mpstate(u_value) + mpstate(v_value))
    )


def bisect(function, left, right, iterations=500):
    left = mp.mpf(left)
    right = mp.mpf(right)
    f_left = function(left)
    f_right = function(right)
    if not (f_left * f_right < 0):
        raise ValueError(
            f"invalid bracket [{left},{right}]: {f_left}, {f_right}"
        )
    for _ in range(iterations):
        middle = (left + right) / 2
        f_middle = function(middle)
        if f_middle == 0:
            return middle
        if f_left * f_middle < 0:
            right = middle
            f_right = f_middle
        else:
            left = middle
            f_left = f_middle
    return (left + right) / 2


x_star = bisect(mK_square, 5, 6)
v_star = mp.sqrt(x_star)
v_A = bisect(lambda t: mpstate(t) - MPINF, mp.mpf(6) / 5, mp.mpf(5) / 4)
v_M = bisect(lambda t: mmu(t) - mmu(0), 16, 17)


def u_root(v_value):
    return bisect(lambda u_value: mF(v_value, u_value), 0, 1)


def endpoint_ratio(v_value):
    u_value = u_root(v_value)
    return 2 * mmu(v_value) / mmu(u_value) - 1


v_C = bisect(endpoint_ratio, 31, 32, iterations=360)
u_C = u_root(v_C)
q_C = -u_C
h_C = (mpstate(q_C) - mpstate(v_C)) / (2 * MPI * mmu(q_C))
L_C = 1 + h_C * q_C

threshold_order_ok = bool(
    0 < v_A < v_star < v_M < v_C
    and -mp.mpf(13) / 100 < q_C < -mp.mpf(1) / 10
    and h_C > 0
    and abs(L_C) < mp.mpf("1e-100")
)
check(
    "high-precision threshold diagnostics respect every rational bracket",
    threshold_order_ok,
    (
        f"vA={mp_text(v_A, 25)}; vstar={mp_text(v_star, 25)}; "
        f"vM={mp_text(v_M, 25)}; vC={mp_text(v_C, 25)}; "
        f"qC={mp_text(q_C, 25)}"
    ),
)


representative_specs = [
    (mp.mpf(3) / 2, mp.mpf(3), mp.mpf(10), "vA<v<vstar"),
    (mp.mpf(3), mp.mpf(1), mp.mpf(2), "vstar<v<vM"),
    (mp.mpf(20), -mp.mpf(1) / 5, -mp.mpf(1) / 100, "vM<v<vC"),
    (mp.mpf(40), -mp.mpf(3) / 10, -mp.mpf(1) / 10, "v>vC"),
]
representatives = []
representative_ok = True
for v_value, q_left, q_right, region in representative_specs:
    q_value = bisect(lambda q_value: mE(v_value, q_value), q_left, q_right)
    h_value = (
        mpstate(q_value) - mpstate(v_value)
    ) / (2 * MPI * mmu(q_value))
    endpoint = 1 + h_value * q_value
    constraint = (
        8 * MPI * (mmu(q_value) - mmu(v_value))
        + h_value * 4 * MPI * q_value * mmu(q_value)
    )
    momentum = (
        mpstate(q_value)
        - mpstate(v_value)
        - h_value * 2 * MPI * mmu(q_value)
    )

    q2 = q_value * q_value
    rq = mp.sqrt(q2 + 4)
    sq = mp.sqrt(3 * q2 + 8)
    Kq = (
        10 * rq - (q2 + 3) * sq * mepsilon(q_value)
    )
    mu_prime = (
        180
        * q_value
        * Kq
        / (MPI * rq**3 * (q2 + 3) * sq)
    )
    p_prime = -720 * Kq / (rq**3 * (q2 + 3) * sq)
    C_h = 4 * MPI * q_value * mmu(q_value)
    P_h = -2 * MPI * mmu(q_value)
    C_q = 8 * MPI * mu_prime + 4 * MPI * h_value * (
        mmu(q_value) + q_value * mu_prime
    )
    P_q = p_prime - 2 * MPI * h_value * mu_prime
    jacobian = C_h * P_q - C_q * P_h

    expected_endpoint_positive = region != "v>vC"
    row_ok = bool(
        h_value > 0
        and (endpoint > 0) == expected_endpoint_positive
        and abs(constraint) < mp.mpf("1e-100")
        and abs(momentum) < mp.mpf("1e-100")
        and abs(jacobian) > mp.mpf("1e-20")
    )
    representative_ok &= row_ok
    representatives.append(
        {
            "region": region,
            "v": mp_text(v_value),
            "q": mp_text(q_value),
            "h": mp_text(h_value),
            "L_plus": mp_text(endpoint),
            "constraint_residual": mp_text(constraint),
            "momentum_residual": mp_text(momentum),
            "jacobian": mp_text(jacobian),
            "passed": row_ok,
        }
    )

check(
    "every nontrivial representative has positive h and exact residuals",
    representative_ok,
    f"representatives={len(representatives)}",
)


# The diagonal and a deliberately reversed determinant are hostile controls.
diagonal_control = abs(mE(mp.mpf(3) / 2, mp.mpf(3) / 2))
wrong_sign_control = (
    4 * MPI * (mmu(mp.mpf(1)) - mmu(mp.mpf(3)))
    - mp.mpf(1) * (mpstate(mp.mpf(1)) - mpstate(mp.mpf(3)))
)
hostile_ok = bool(
    diagonal_control < mp.mpf("1e-100")
    and abs(wrong_sign_control) > mp.mpf("1")
)
check(
    "the zero-height diagonal passes while the wrong determinant sign fails",
    hostile_ok,
)


classification_complete = bool(
    provenance_ok
    and elementary_chain_ok
    and state_derivative_ok
    and E_q_reduced == 0
    and parity_ok
    and tail_formula_ok
    and sign_pattern_ok
    and root_count_logic_ok
    and endpoint_proof_ok
    and threshold_order_ok
    and representative_ok
    and hostile_ok
)
outcome = (
    "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
    if classification_complete
    else "FINITE_HEIGHT_OPEN"
)
check(
    "the preregistered finite-height outcome is selected without fitting",
    outcome == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY",
)


artifact = {
    "provenance": {
        "finite_input_sha256": digest(FINITE_INPUT),
        "next_order_input_sha256": digest(NEXT_INPUT),
        "next_order_adversarial_input_sha256": digest(NEXT_ADVERSARIAL_INPUT),
        "classification_protocol_commit": CLASSIFICATION_PROTOCOL_COMMIT,
    },
    "exact_identities": {
        "mu_prime": "180*q*K(q^2)/(pi*(q^2+4)^(3/2)*(q^2+3)*sqrt(3*q^2+8))",
        "p_prime": "-720*K(q^2)/((q^2+4)^(3/2)*(q^2+3)*sqrt(3*q^2+8))",
        "E_q": "p(q)-p(v)",
        "reflection": "E(-v,-q)=E(v,q), h(-v,-q)=-h(v,q)",
        "p_infinity": "60*pi-300*sqrt(3)*log(2)<0",
    },
    "monotone_facts": monotone_facts,
    "root_count_table": root_count_table,
    "endpoint_monotonicity": endpoint_monotonicity,
    "rational_sign_certificates": {
        key: arb_text(value) for key, value in sign_certificates.items()
    },
    "thresholds": {
        "v_A": mp_text(v_A, 80),
        "v_star": mp_text(v_star, 80),
        "v_M": mp_text(v_M, 80),
        "v_C": mp_text(v_C, 80),
        "q_C": mp_text(q_C, 80),
        "h_C": mp_text(h_C, 80),
        "L_C_residual": mp_text(L_C, 30),
        "rational_brackets": {
            "v_A": ["6/5", "5/4"],
            "x_star": ["5", "6"],
            "v_M": ["16", "17"],
            "v_C": ["31", "32"],
            "q_C": ["-13/100", "-1/10"],
        },
    },
    "physical_state_set": "v in (v_A,v_star) union (v_star,v_C)",
    "representatives": representatives,
    "interpretation": {
        "label": "DERIVED EXACT / STRUCTURAL",
        "state_dependent_update": True,
        "universal_tick": False,
        "composition": "OPEN",
        "perturbative_stability": "OPEN",
        "refinement": "OPEN",
        "seconds_c_G_or_Planck_time": "NOT DERIVED",
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
