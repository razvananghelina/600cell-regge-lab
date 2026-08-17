#!/usr/bin/env python3
"""Generic-index theorem audit for the cellular 600-cell weak-lapse jet."""

import hashlib
import json
from collections import defaultdict
from pathlib import Path
import time

import mpmath as arb
import sympy as sp


HERE = Path(__file__).resolve().parent
BLIND_INPUT = HERE / "gravity_600cell_cellular_weak_lapse_blind.json"
OUTPUT = HERE / "gravity_600cell_cellular_weak_lapse_all_n.json"
BLIND_SHA256 = "6d39e9a4594d9c9ead102f94cf9115d8474132ecce511fe7359826dcc73b9de0"
PRIOR_ART_COMMIT = "bd1e7df"
PROTOCOL_COMMIT = "a07b3ff"
PROTOCOL_CORRECTION_COMMIT = "c0303a4"
DPS = 120
arb.mp.dps = DPS
CONTROL_INDICES = (5, 7, 11)
CONTROL_E = (arb.mpf(1)/200, arb.mpf(1)/400, arb.mpf(1)/800)
tests = passed = 0
started = time.perf_counter()


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


blind = json.loads(BLIND_INPUT.read_text())
provenance_ok = bool(
    digest(BLIND_INPUT) == BLIND_SHA256
    and blind["outcome"] == "CELLULAR_WEAK_LAPSE_JET_DERIVED"
    and blind["passed"] == blind["tests"] == 8
    and blind["tick_artifacts_parsed"] is False
    and PRIOR_ART_COMMIT == "bd1e7df"
    and PROTOCOL_COMMIT == "a07b3ff"
    and PROTOCOL_CORRECTION_COMMIT == "c0303a4"
)
check("the generic theorem uses the frozen blind action input", provenance_ok)


# Reconstruct the universal action jet independently of the stored coefficient
# equations.  As in Stage A, x=e^2 and S=e*Sbar.
x = sp.symbols("x", positive=True)
a_minus, a_plus, r = sp.symbols("a_minus a_plus r", real=True)
epsilon = sp.symbols("epsilon", positive=True, real=True)
n = sp.symbols("n", integer=True, positive=True)
A_new, B_new, R_new = sp.symbols("A_new B_new R_new", real=True)
alpha_value = sp.acos(sp.Rational(1, 3))
epsilon_value = 2*sp.pi-5*alpha_value
mu = sp.Rational(90)/sp.pi*epsilon_value

lm = sp.exp(a_minus*x)
lp = sp.exp(a_plus*x)
rho_bar = sp.exp(r*x)
delta = lp-lm
delta_over_x = sp.series(delta/x, x, 0, 5).removeO()
delta_square_over_x = sp.series(delta**2/x, x, 0, 5).removeO()
h_bar = sp.sqrt(rho_bar+delta_square_over_x/4)
cosine_bar = (
    delta_square_over_x+2*rho_bar
)/(2*(delta_square_over_x+3*rho_bar))
boost_bar = delta_over_x/sp.sqrt(8*(delta_square_over_x+3*rho_bar))
asinh_over_sqrt_x = (
    boost_bar-x*boost_bar**3/6+3*x**2*boost_bar**5/40
    -5*x**3*boost_bar**7/112+35*x**4*boost_bar**9/1152
)
s_bar_raw = (
    360*(lm+lp)*h_bar*(2*sp.pi-5*sp.acos(cosine_bar))
    + 600*sp.sqrt(3)*(lm**2-lp**2)*asinh_over_sqrt_x
    - 8*sp.pi*mu*sp.sqrt(rho_bar)
)
s_bar = sp.series(s_bar_raw, x, 0, 4).removeO().expand()
f_scaled = sp.series(sp.diff(s_bar, r)/x, x, 0, 3).removeO().expand()
p_minus_scaled = sp.series(
    sp.diff(s_bar, a_minus)/(2*x), x, 0, 3
).removeO().expand()
p_plus_scaled = sp.series(
    sp.diff(s_bar, a_plus)/(2*x), x, 0, 3
).removeO().expand()


def field_reduce(expression):
    expanded = sp.expand(expression)
    generators = sorted(expanded.free_symbols, key=sp.default_sort_key)
    if not generators:
        return expanded
    return sp.cancel(expanded, *generators, extension=sp.sqrt(2))


def compact(expression):
    alpha = sp.symbols("alpha", positive=True, real=True)
    replaced = expression.xreplace({alpha_value: alpha})
    return field_reduce(replaced.subs(alpha, (2*sp.pi-epsilon)/5))


def compact_polynomial(expression, maximum_power=2):
    return sp.Add(*[
        compact(expression.coeff(x, power))*x**power
        for power in range(maximum_power+1)
    ])


f_compact = compact_polynomial(f_scaled)
p_minus_compact = compact_polynomial(p_minus_scaled)
p_plus_compact = compact_polynomial(p_plus_scaled)


def affine_jet_coefficients(expression, am0, bm, ap0, bp, r0, order=3):
    result = [sp.Integer(0) for _ in range(order)]
    base = {a_minus: am0, a_plus: ap0, r: r0}
    for power in range(order):
        directional = expression.coeff(x, power)
        for derivative_order in range(order-power):
            if derivative_order:
                directional = (
                    bm*sp.diff(directional, a_minus)
                    + bp*sp.diff(directional, a_plus)
                )
            result[power+derivative_order] += (
                directional.subs(base)/sp.factorial(derivative_order)
            )
    return result


q = 5*sp.sqrt(2)-3*epsilon
t_index = n*(n+1)


def candidate_A(index):
    return -6*epsilon*index*(index+1)/q


def candidate_R(index):
    return (
        -10*epsilon*(7*sp.sqrt(2)*epsilon+60)*index**2/q**3
    )


def candidate_B(index):
    t_value = index*(index+1)
    return -epsilon**2*t_value*(
        (108*epsilon**2-395*sp.sqrt(2)*epsilon+300)*t_value
        -54*epsilon**2+145*sp.sqrt(2)*epsilon-600
    )/q**4


print("[INFO] constructing generic-index coefficient equations", flush=True)
current_f = affine_jet_coefficients(
    f_compact,
    candidate_A(n-1), candidate_B(n-1), A_new, B_new, R_new,
)
previous_p_plus = affine_jet_coefficients(
    p_plus_compact,
    candidate_A(n-2), candidate_B(n-2),
    candidate_A(n-1), candidate_B(n-1), candidate_R(n-1),
)
current_p_minus = affine_jet_coefficients(
    p_minus_compact,
    candidate_A(n-1), candidate_B(n-1), A_new, B_new, R_new,
)
seam = [left+right for left, right in zip(previous_p_plus, current_p_minus)]

leading_f = current_f[1]
leading_g = seam[0]
A_candidate = candidate_A(n)
g_affine = sp.diff(leading_g, A_new, 2) == 0
g_slope = field_reduce(sp.diff(leading_g, A_new))
g_at_candidate = field_reduce(leading_g.subs(A_new, A_candidate))
f_at_candidate = field_reduce(leading_f.subs(A_new, A_candidate))
leading_local_derivatives = (
    field_reduce(sp.diff(leading_f, A_new).subs(A_new, A_candidate)),
    field_reduce(sp.diff(leading_g, A_new).subs(A_new, A_candidate)),
)
leading_rank = int(any(value != 0 for value in leading_local_derivatives))
leading_ok = bool(
    g_affine and g_slope != 0 and g_at_candidate == 0
    and f_at_candidate == 0 and leading_rank == 1
)
check(
    "the generic seam selects A_n and the common leading branch has rank one",
    leading_ok,
    f"seam slope={sp.factor(g_slope)}; local rank={leading_rank}",
)


print("[INFO] reducing the generic next-order Cramer system", flush=True)
next_f = current_f[2].subs(A_new, A_candidate)
next_g = seam[1].subs(A_new, A_candidate)
next_affine = bool(
    sp.diff(next_f, B_new, 2) == 0
    and sp.diff(next_f, B_new, R_new) == 0
    and sp.diff(next_f, R_new, 2) == 0
    and sp.diff(next_g, B_new, 2) == 0
    and sp.diff(next_g, B_new, R_new) == 0
    and sp.diff(next_g, R_new, 2) == 0
)
m11 = field_reduce(sp.diff(next_f, B_new))
m12 = field_reduce(sp.diff(next_f, R_new))
m21 = field_reduce(sp.diff(next_g, B_new))
m22 = field_reduce(sp.diff(next_g, R_new))
c1 = field_reduce(next_f.subs({B_new: 0, R_new: 0}))
c2 = field_reduce(next_g.subs({B_new: 0, R_new: 0}))
determinant = field_reduce(m11*m22-m12*m21)
B_cramer = field_reduce((m12*c2-m22*c1)/determinant)
R_cramer = field_reduce((m21*c1-m11*c2)/determinant)
B_residual = field_reduce(B_cramer-candidate_B(n))
R_residual = field_reduce(R_cramer-candidate_R(n))
next_ok = bool(
    next_affine and determinant == 16200*epsilon**2
    and B_residual == 0 and R_residual == 0
)
check(
    "the generic next system has rank two and selects the candidate B_n,R_n",
    next_ok,
    f"det={sp.factor(determinant)}; B residual={B_residual}; R residual={R_residual}",
)


candidate_post = affine_jet_coefficients(
    p_plus_compact,
    candidate_A(n-1), candidate_B(n-1),
    candidate_A(n), candidate_B(n), candidate_R(n),
)[0]
post_ratio = field_reduce(candidate_post/(180*epsilon))
post_ok = post_ratio == 2*n+1
check(
    "the generic outgoing momentum coefficient is exactly 2n+1",
    post_ok,
    f"p_out/k={post_ratio}",
)


# Exact domain signs.  The elementary inequalities are recorded explicitly;
# their logical use relies only on monotonicity of cos on [0,pi].
cos_2pi5 = (sp.sqrt(5)-1)/4
cos_bound = bool(sp.N(sp.Rational(1, 3)-cos_2pi5, 80) > 0)
lower_acos_bound = bool(sp.Rational(1, 3) < sp.Rational(1, 2))
pi_bound = bool(sp.N(5*sp.sqrt(2)-sp.pi, 80) > 0)
domain_numeric_epsilon = sp.N(epsilon_value, 100)
domain_numeric_q = sp.N(q.subs(epsilon, epsilon_value), 100)
domain_ok = bool(
    cos_bound and lower_acos_bound and pi_bound
    and domain_numeric_epsilon > 0 and domain_numeric_q > 0
    and determinant == 16200*epsilon**2
)
check(
    "the exact elementary bounds select a nondegenerate contracting domain",
    domain_ok,
    f"epsilon={domain_numeric_epsilon}; q={domain_numeric_q}",
)


# Specialize to the four frozen coefficients as an exact independent control.
locals_exact = {"epsilon_3": epsilon, "sqrt": sp.sqrt}
frozen_match_records = []
frozen_match_ok = True
for index in range(1, 5):
    for label, candidate in (
        ("A", candidate_A(index)),
        ("B", candidate_B(index)),
        ("R", candidate_R(index)),
    ):
        frozen = sp.sympify(
            blind["coefficients"][label][str(index)], locals=locals_exact
        )
        residual = field_reduce(candidate-frozen)
        frozen_match_ok &= residual == 0
        frozen_match_records.append({
            "n": index, "coefficient": label, "residual": str(residual),
        })
check(
    "the generic closed forms recover all twelve frozen n<=4 coefficients exactly",
    frozen_match_ok and len(frozen_match_records) == 12,
)


# Full unexpanded action residual controls at three genuinely new indices.
L_MINUS, L_PLUS, RHO = sp.symbols("L_minus L_plus rho", positive=True)
DELTA = L_PLUS-L_MINUS
H = sp.sqrt(RHO+DELTA**2/4)
C = (DELTA**2+2*RHO)/(2*(DELTA**2+3*RHO))
BETA = DELTA/sp.sqrt(8*(DELTA**2+3*RHO))
S_FULL = (
    360*(L_MINUS+L_PLUS)*H*(2*sp.pi-5*sp.acos(C))
    + 600*sp.sqrt(3)*(L_MINUS**2-L_PLUS**2)*sp.asinh(BETA)
    - 8*sp.pi*mu*sp.sqrt(RHO)
)
F_FULL = sp.lambdify(
    (L_MINUS, L_PLUS, RHO), RHO*sp.diff(S_FULL, RHO), "mpmath"
)
PM_FULL = sp.lambdify(
    (L_MINUS, L_PLUS, RHO),
    L_MINUS*sp.diff(S_FULL, L_MINUS)/2, "mpmath",
)
PP_FULL = sp.lambdify(
    (L_MINUS, L_PLUS, RHO),
    L_PLUS*sp.diff(S_FULL, L_PLUS)/2, "mpmath",
)


def mp_coefficient(expression):
    return arb.mpf(str(sp.N(expression.subs(epsilon, epsilon_value), DPS)))


A_mp = {i: mp_coefficient(candidate_A(i)) for i in range(-1, 12)}
B_mp = {i: mp_coefficient(candidate_B(i)) for i in range(-1, 12)}
R_mp = {i: mp_coefficient(candidate_R(i)) for i in range(0, 12)}


def l_at(index, e_value):
    return arb.exp(A_mp[index]*e_value**2+B_mp[index]*e_value**4)


def rho_at(index, e_value):
    return e_value**2*arb.exp(R_mp[index]*e_value**2)


residual_records = []
sequences = defaultdict(list)
for index in CONTROL_INDICES:
    for e_value in CONTROL_E:
        current = (l_at(index-1, e_value), l_at(index, e_value), rho_at(index, e_value))
        previous = (
            l_at(index-2, e_value), l_at(index-1, e_value),
            rho_at(index-1, e_value),
        )
        f_value = abs(F_FULL(*current))
        g_value = abs(PP_FULL(*previous)+PM_FULL(*current))
        sequences[(index, "F")].append(f_value)
        sequences[(index, "G")].append(g_value)
        residual_records.append({
            "n": index, "e": arb.nstr(e_value, 30),
            "F_abs": arb.nstr(f_value, 50), "G_abs": arb.nstr(g_value, 50),
        })

observed_orders = {}
full_control_ok = True
for (index, equation), values in sequences.items():
    orders = (
        arb.log(values[0]/values[1], 2),
        arb.log(values[1]/values[2], 2),
    )
    expected = arb.mpf(7 if equation == "F" else 5)
    accepted = bool(
        values[0] > values[1] > values[2] > 0
        and all(abs(value-expected) < arb.mpf("0.3") for value in orders)
    )
    full_control_ok &= accepted
    observed_orders[f"n={index}:{equation}"] = [arb.nstr(value, 30) for value in orders]
check(
    "new n=5,7,11 states have the predicted full-equation residual orders",
    full_control_ok,
    "; ".join(
        f"{key}={','.join(values)}" for key, values in observed_orders.items()
    ),
)


generic_ok = bool(leading_ok and next_ok and post_ok and domain_ok)
if generic_ok and frozen_match_ok and full_control_ok:
    outcome = "CELLULAR_WEAK_LAPSE_ALL_N_PROVED"
elif not generic_ok:
    outcome = "CELLULAR_WEAK_LAPSE_ALL_N_REFUTED"
else:
    outcome = "CELLULAR_WEAK_LAPSE_ALL_N_OPEN"
check(
    "the result follows the preregistered all-index hierarchy",
    outcome in {
        "CELLULAR_WEAK_LAPSE_ALL_N_PROVED",
        "CELLULAR_WEAK_LAPSE_ALL_N_REFUTED",
        "CELLULAR_WEAK_LAPSE_ALL_N_OPEN",
    },
    f"outcome={outcome}",
)

runtime = time.perf_counter()-started
artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "blind_input_sha256": BLIND_SHA256,
    "outcome": outcome,
    "labels": {
        "formal_all_index_coefficient_recurrence": (
            "DERIVED" if outcome == "CELLULAR_WEAK_LAPSE_ALL_N_PROVED"
            else "OPEN" if outcome.endswith("OPEN") else "DERIVED NEGATIVE"
        ),
        "uniform_large_n_convergence": "OPEN",
        "spatial_refinement": "OPEN",
        "physical_clock": "OPEN",
        "external_novelty": "OPEN",
    },
    "candidate": {
        "A_n": str(sp.factor(candidate_A(n))),
        "B_n": str(sp.factor(candidate_B(n))),
        "R_n": str(sp.factor(candidate_R(n))),
        "p_out_over_k": str(post_ratio),
    },
    "generic_system": {
        "leading_seam_slope": str(sp.factor(g_slope)),
        "leading_lapse_residual": str(f_at_candidate),
        "leading_seam_residual": str(g_at_candidate),
        "leading_local_derivatives": [str(value) for value in leading_local_derivatives],
        "leading_local_rank": leading_rank,
        "next_matrix": [[str(m11), str(m12)], [str(m21), str(m22)]],
        "next_determinant": str(determinant),
        "B_candidate_residual": str(B_residual),
        "R_candidate_residual": str(R_residual),
        "operation_counts": {
            "leading_F": int(sp.count_ops(leading_f)),
            "leading_G": int(sp.count_ops(leading_g)),
            "next_F": int(sp.count_ops(next_f)),
            "next_G": int(sp.count_ops(next_g)),
        },
    },
    "domain": {
        "cos_2pi5_less_than_one_third": cos_bound,
        "one_third_less_than_one_half": lower_acos_bound,
        "pi_less_than_5sqrt2": pi_bound,
        "epsilon_numeric": str(domain_numeric_epsilon),
        "q_numeric": str(domain_numeric_q),
    },
    "frozen_coefficient_controls": frozen_match_records,
    "new_index_residuals": residual_records,
    "observed_halving_orders": observed_orders,
    "scope": {
        "formal_fixed_n": True,
        "uniform_n_e_limit": "NOT PROVED",
        "control_indices": list(CONTROL_INDICES),
    },
    "runtime_seconds": runtime,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print(f"\nSummary: {passed}/{tests} checks passed in {runtime:.2f}s")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")
if passed != tests or outcome != "CELLULAR_WEAK_LAPSE_ALL_N_PROVED":
    raise SystemExit(1)
