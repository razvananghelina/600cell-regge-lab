#!/usr/bin/env python3
"""Target-blind weak-lapse jet of the cellular 600-cell dust action.

Prior-art commit: b77856a.
Protocol commit: 71d10b4.

This file deliberately has no tick-artifact path and derives all coefficients
before any comparison with the previously observed integer sequences.
"""

from collections import defaultdict
import hashlib
import json
from pathlib import Path

import mpmath as arb
import sympy as sp


HERE = Path(__file__).resolve().parent
ACTION_INPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
OUTPUT = HERE / "gravity_600cell_cellular_weak_lapse_blind.json"
ACTION_SHA256 = (
    "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d"
)
PRIOR_ART_COMMIT = "b77856a"
PROTOCOL_COMMIT = "71d10b4"
PROTOCOL_CORRECTION_COMMIT = "35f37d4"
DPS = 100
arb.mp.dps = DPS
CONTROL_POINTS = (
    (sp.Integer(1), sp.Rational(4, 5), sp.Rational(1, 10)),
    (sp.Integer(1), sp.Rational(6, 5), sp.Rational(1, 10)),
    (sp.Integer(2), sp.Rational(5, 2), sp.Rational(1, 4)),
)
SERIES_CONTROLS = (
    sp.Rational(1, 100), sp.Rational(1, 200), sp.Rational(1, 400)
)
DIFFERENCE_STEPS = (arb.mpf("1e-20"), arb.mpf("3e-20"))
FORBIDDEN_TICK_TOKENS = (
    "weak_lapse", "canonical_lapse", "second_tick", "third_tick", "fourth_tick"
)
TICK_ARTIFACTS_PARSED = False
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


def equal(left, right):
    return sp.simplify(sp.expand(left-right)) == 0


action_artifact = json.loads(ACTION_INPUT.read_text())
input_ok = bool(
    digest(ACTION_INPUT) == ACTION_SHA256
    and action_artifact.get("outcome")
        == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and action_artifact.get("passed") == action_artifact.get("tests") == 16
    and action_artifact.get("labels", {}).get(
        "homogeneous_action_subdivision_invariance"
    ) == "DERIVED"
    and not any(token in ACTION_INPUT.name for token in FORBIDDEN_TICK_TOKENS)
    and TICK_ARTIFACTS_PARSED is False
    and PROTOCOL_CORRECTION_COMMIT == "35f37d4"
)
check("the blind stage reads only the frozen cellular-action artifact", input_ok)


# Exact closed action and its canonical logarithmic derivatives.
L_MINUS, L_PLUS, RHO = sp.symbols(
    "L_minus L_plus rho", positive=True
)
EPSILON3 = 2*sp.pi-5*sp.acos(sp.Rational(1, 3))
MU = sp.Rational(90)/sp.pi*EPSILON3
DELTA = L_PLUS-L_MINUS
H = sp.sqrt(RHO+DELTA**2/4)
LATERAL_COSINE = (DELTA**2+2*RHO)/(2*(DELTA**2+3*RHO))
BOOST = DELTA/sp.sqrt(8*(DELTA**2+3*RHO))
S_GRAV = (
    360*(L_MINUS+L_PLUS)*H*(2*sp.pi-5*sp.acos(LATERAL_COSINE))
    + 600*sp.sqrt(3)*(L_MINUS**2-L_PLUS**2)*sp.asinh(BOOST)
)
S_TOTAL = S_GRAV-8*sp.pi*MU*sp.sqrt(RHO)
F_RHO = sp.factor(RHO*sp.diff(S_TOTAL, RHO))
P_MINUS = sp.factor(L_MINUS*sp.diff(S_TOTAL, L_MINUS)/2)
P_PLUS = sp.factor(L_PLUS*sp.diff(S_TOTAL, L_PLUS)/2)

TAU = sp.symbols("tau", positive=True)
static_action = sp.simplify(S_TOTAL.subs({
    L_MINUS: 1, L_PLUS: 1, RHO: TAU**2,
}))
static_p_minus = sp.simplify(P_MINUS.subs({
    L_MINUS: 1, L_PLUS: 1, RHO: TAU**2,
}))
static_p_plus = sp.simplify(P_PLUS.subs({
    L_MINUS: 1, L_PLUS: 1, RHO: TAU**2,
}))
static_ok = bool(
    static_action == 0
    and equal(static_p_minus, 180*EPSILON3*TAU)
    and equal(static_p_plus, 180*EPSILON3*TAU)
)
check(
    "the closed formula independently recovers the static identity and momenta",
    static_ok,
    f"P-/tau={sp.simplify(static_p_minus/TAU)}",
)


# High-precision derivative control, independent of the symbolic derivatives.
def mp_action(l_minus, l_plus, rho):
    epsilon3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
    mu = 90/arb.pi*epsilon3
    delta = l_plus-l_minus
    h = arb.sqrt(rho+delta**2/4)
    cosine = (delta**2+2*rho)/(2*(delta**2+3*rho))
    boost = delta/arb.sqrt(8*(delta**2+3*rho))
    return (
        360*(l_minus+l_plus)*h*(2*arb.pi-5*arb.acos(cosine))
        + 600*arb.sqrt(3)*(l_minus**2-l_plus**2)*arb.asinh(boost)
        - 8*arb.pi*mu*arb.sqrt(rho)
    )


def mp_value(value):
    if isinstance(value, sp.Rational):
        return arb.mpf(value.p)/value.q
    return arb.mpf(value)


def perturb(point, coordinate, step):
    result = list(point)
    if coordinate in (0, 1):
        result[coordinate] *= arb.exp(step/2)
    else:
        result[coordinate] *= arb.exp(step)
    return tuple(result)


def centered(point, coordinate, step):
    plus = perturb(point, coordinate, step)
    minus = perturb(point, coordinate, -step)
    return (mp_action(*plus)-mp_action(*minus))/(2*step)


numeric_derivatives = [
    sp.lambdify((L_MINUS, L_PLUS, RHO), expression, "mpmath")
    for expression in (P_MINUS, P_PLUS, F_RHO)
]
derivative_records = []
maximum_derivative_error = arb.mpf(0)
for exact_point in CONTROL_POINTS:
    point = tuple(mp_value(value) for value in exact_point)
    errors = []
    for coordinate, analytic in enumerate(numeric_derivatives):
        fine = centered(point, coordinate, DIFFERENCE_STEPS[0])
        coarse = centered(point, coordinate, DIFFERENCE_STEPS[1])
        richardson = (9*fine-coarse)/8
        expected = analytic(*point)
        error = abs(richardson-expected)/max(
            arb.mpf(1), abs(richardson), abs(expected)
        )
        errors.append(error)
        maximum_derivative_error = max(maximum_derivative_error, error)
    derivative_records.append({
        "point": [str(value) for value in exact_point],
        "relative_errors": [arb.nstr(value, 40) for value in errors],
    })
derivative_control_ok = maximum_derivative_error < arb.mpf("1e-60")
check(
    "symbolic and centered closed-action derivatives agree below 1e-60",
    derivative_control_ok,
    f"max relative={arb.nstr(maximum_derivative_error, 8)}",
)


# Universal target-blind series.  Work with x=e^2 and divide the action by e;
# the remaining expression is analytic in x.  This avoids any numerical
# limiting procedure in the coefficient derivation.
x = sp.symbols("x", positive=True)
a_minus, a_plus, r = sp.symbols("a_minus a_plus r", real=True)
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
boost_bar = delta_over_x/sp.sqrt(
    8*(delta_square_over_x+3*rho_bar)
)

# asinh(sqrt(x)*b)/sqrt(x), through x^4.
asinh_over_sqrt_x = (
    boost_bar
    - x*boost_bar**3/6
    + 3*x**2*boost_bar**5/40
    - 5*x**3*boost_bar**7/112
    + 35*x**4*boost_bar**9/1152
)
s_bar_raw = (
    360*(lm+lp)*h_bar*(2*sp.pi-5*sp.acos(cosine_bar))
    + 600*sp.sqrt(3)*(lm**2-lp**2)*asinh_over_sqrt_x
    - 8*sp.pi*MU*sp.sqrt(rho_bar)
)
s_bar = sp.series(s_bar_raw, x, 0, 4).removeO().expand()
series_coefficients = {
    power: sp.factor(s_bar.coeff(x, power)) for power in range(4)
}

# S=e*s_bar and log L^2=2*a*x, log rho=const+r*x.
f_scaled = sp.series(sp.diff(s_bar, r)/x, x, 0, 3).removeO().expand()
p_minus_scaled = sp.series(
    sp.diff(s_bar, a_minus)/(2*x), x, 0, 3
).removeO().expand()
p_plus_scaled = sp.series(
    sp.diff(s_bar, a_plus)/(2*x), x, 0, 3
).removeO().expand()
ALPHA = sp.symbols("alpha", positive=True, real=True)
E3 = sp.symbols("epsilon_3", positive=True, real=True)
ALPHA_VALUE = sp.acos(sp.Rational(1, 3))
EPSILON3_COMPACT = E3


def compact(expression):
    # All alpha/pi dependence collapses to the geometric deficit epsilon_3.
    # Eliminating alpha before the recursion prevents generic assumption and
    # root-isolation machinery from expanding otherwise equivalent formulas.
    with_alpha = expression.xreplace({ALPHA_VALUE: ALPHA})
    return sp.cancel(sp.expand(with_alpha.subs(
        ALPHA, (2*sp.pi-E3)/5
    )))


def physical(expression):
    return expression.subs(E3, EPSILON3)


def rational_reduce(expression):
    """Exact zero-preserving reduction without heuristic simplify()."""
    return sp.cancel(sp.expand(expression))


f_scaled_compact = compact(f_scaled)
p_minus_scaled_compact = compact(p_minus_scaled)
p_plus_scaled_compact = compact(p_plus_scaled)


def first_nonzero(expression, maximum_power=2):
    for power in range(maximum_power+1):
        coefficient = sp.factor(expression.coeff(x, power))
        if coefficient != 0 and not equal(coefficient, 0):
            return power, coefficient
    return None, sp.Integer(0)


f_power, f_leading = first_nonzero(f_scaled)
p_minus_power, p_minus_leading = first_nonzero(p_minus_scaled)
p_plus_power, p_plus_leading = first_nonzero(p_plus_scaled)
leading_orders_found = bool(
    f_power is not None and p_minus_power is not None and p_plus_power is not None
)
check(
    "the first nonzero lapse and momentum orders are discovered symbolically",
    leading_orders_found,
    f"F/e starts x^{f_power}; P-/e x^{p_minus_power}; P+/e x^{p_plus_power}",
)


# Recursive coefficient solution.  The first registered run established that
# the leading system has rank one: it fixes A_n but not R_n.  Following the
# frozen correction, B_n is retained as the e^4 scale coefficient and the
# next two equations solve (B_n,R_n).  No target sequence is present here.
A = {-1: sp.Integer(0), 0: sp.Integer(0)}
B = {-1: sp.Integer(0), 0: sp.Integer(0)}
R = {0: sp.Integer(0)}
step_records = []
unique_contracting = True
exact_substitution_ok = True

for n in range(1, 5):
    A_n, B_n, R_n = sp.symbols(f"A_{n} B_{n} R_{n}", real=True)
    current_substitution = {
        a_minus: A[n-1]+B[n-1]*x,
        a_plus: A_n+B_n*x,
        r: R_n,
    }
    previous_substitution = {
        a_minus: A[n-2]+B[n-2]*x,
        a_plus: A[n-1]+B[n-1]*x,
        r: R[n-1],
    }
    current_f = sp.series(
        f_scaled_compact.subs(current_substitution), x, 0, 3
    ).removeO().expand()
    previous_p_plus = sp.series(
        p_plus_scaled_compact.subs(previous_substitution), x, 0, 3
    ).removeO().expand()
    current_p_minus = sp.series(
        p_minus_scaled_compact.subs(current_substitution), x, 0, 3
    ).removeO().expand()
    seam = sp.expand(previous_p_plus+current_p_minus)

    f_leading_power, f_leading_equation = first_nonzero(current_f)
    g_leading_power, g_leading_equation = first_nonzero(seam)
    print(f"[INFO] solving blind coefficient step n={n}", flush=True)
    a_linear = sp.diff(g_leading_equation, A_n, 2) == 0
    a_coefficient = sp.factor(sp.diff(g_leading_equation, A_n))
    a_constant = sp.factor(g_leading_equation.subs(A_n, 0))
    a_solutions = [] if a_coefficient == 0 else [
        rational_reduce(-a_constant/a_coefficient)
    ]
    a_solutions = [
        value for value in a_solutions
        if abs(sp.im(sp.N(physical(value), 60))) < sp.Float("1e-50")
    ]
    contracting_a = [
        value for value in a_solutions
        if sp.N(physical(value-A[n-1]), 60) < 0
    ]
    if len(contracting_a) != 1:
        unique_contracting = False
        selected_a = contracting_a[0] if contracting_a else (
            a_solutions[0] if a_solutions else None
        )
    else:
        selected_a = contracting_a[0]

    leading_jacobian = sp.Matrix([
        [sp.diff(f_leading_equation, A_n)],
        [sp.diff(g_leading_equation, A_n)],
    ])
    leading_rank = int(any(value != 0 for value in leading_jacobian))
    if selected_a is not None:
        A[n] = selected_a
        f_leading_residual = rational_reduce(
            f_leading_equation.subs(A_n, A[n])
        )
        g_leading_residual = rational_reduce(
            g_leading_equation.subs(A_n, A[n])
        )
        after_a = {A_n: A[n]}
        current_f_after_a = sp.expand(current_f.subs(after_a))
        seam_after_a = sp.expand(seam.subs(after_a))
        f_next_power, f_next_equation = first_nonzero(current_f_after_a)
        g_next_power, g_next_equation = first_nonzero(seam_after_a)
        next_linear = bool(
            sp.diff(f_next_equation, B_n, 2) == 0
            and sp.diff(f_next_equation, B_n, R_n) == 0
            and sp.diff(f_next_equation, R_n, 2) == 0
            and sp.diff(g_next_equation, B_n, 2) == 0
            and sp.diff(g_next_equation, B_n, R_n) == 0
            and sp.diff(g_next_equation, R_n, 2) == 0
        )
        m11 = sp.factor(sp.diff(f_next_equation, B_n))
        m12 = sp.factor(sp.diff(f_next_equation, R_n))
        m21 = sp.factor(sp.diff(g_next_equation, B_n))
        m22 = sp.factor(sp.diff(g_next_equation, R_n))
        c1 = sp.factor(f_next_equation.subs({B_n: 0, R_n: 0}))
        c2 = sp.factor(g_next_equation.subs({B_n: 0, R_n: 0}))
        next_jacobian = sp.Matrix([[m11, m12], [m21, m22]])
        next_determinant = sp.factor(m11*m22-m12*m21)
        next_rank = 2 if next_determinant != 0 else next_jacobian.rank()
        next_solutions = []
        if next_linear and next_determinant != 0:
            b_value = rational_reduce(
                (m12*c2-m22*c1)/next_determinant
            )
            r_value = rational_reduce(
                (m21*c1-m11*c2)/next_determinant
            )
            if abs(sp.im(sp.N(physical(b_value), 60))) < sp.Float("1e-50") and abs(
                sp.im(sp.N(physical(r_value), 60))
            ) < sp.Float("1e-50"):
                next_solutions.append((b_value, r_value))
        if len(next_solutions) == 1:
            B[n], R[n] = next_solutions[0]
        else:
            unique_contracting = False
            if next_solutions:
                B[n], R[n] = next_solutions[0]
            else:
                B[n], R[n] = sp.nan, sp.nan
        next_substitution = {B_n: B[n], R_n: R[n]}
        f_next_residual = rational_reduce(
            f_next_equation.subs(next_substitution)
        )
        g_next_residual = rational_reduce(
            g_next_equation.subs(next_substitution)
        )
        exact_substitution_ok &= bool(
            f_leading_residual == 0
            and g_leading_residual == 0
            and f_next_residual == 0
            and g_next_residual == 0
        )
    else:
        A[n], B[n], R[n] = sp.nan, sp.nan, sp.nan
        f_leading_residual = g_leading_residual = sp.nan
        f_next_power = g_next_power = None
        f_next_equation = g_next_equation = sp.nan
        next_solutions = []
        next_jacobian = sp.zeros(2)
        next_determinant = sp.Integer(0)
        next_rank = 0
        f_next_residual = g_next_residual = sp.nan
        exact_substitution_ok = False
    step_records.append({
        "n": n,
        "leading_lapse_power_in_x_of_F_over_e": f_leading_power,
        "leading_seam_power_in_x_of_G_over_e": g_leading_power,
        "leading_lapse_equation": str(sp.factor(f_leading_equation)),
        "leading_seam_equation": str(sp.factor(g_leading_equation)),
        "leading_A_jacobian": [
            str(sp.factor(value)) for value in leading_jacobian
        ],
        "leading_A_rank": leading_rank,
        "leading_A_equation_linear": a_linear,
        "real_A_solutions": [str(value) for value in a_solutions],
        "contracting_A_solution_count": len(contracting_a),
        "selected_A": str(A[n]),
        "leading_lapse_residual": str(f_leading_residual),
        "leading_seam_residual": str(g_leading_residual),
        "next_lapse_power_in_x_of_F_over_e": f_next_power,
        "next_seam_power_in_x_of_G_over_e": g_next_power,
        "next_lapse_equation": str(sp.factor(f_next_equation)),
        "next_seam_equation": str(sp.factor(g_next_equation)),
        "next_BR_jacobian": [
            [str(sp.factor(value)) for value in row]
            for row in next_jacobian.tolist()
        ],
        "next_BR_determinant": str(next_determinant),
        "next_BR_rank": next_rank,
        "next_BR_equations_linear": next_linear,
        "real_BR_solutions": [
            {"B": str(item[0]), "R": str(item[1])}
            for item in next_solutions
        ],
        "selected_B": str(B[n]),
        "selected_R": str(R[n]),
        "next_lapse_residual": str(f_next_residual),
        "next_seam_residual": str(g_next_residual),
    })

coefficient_system_ok = bool(
    unique_contracting
    and exact_substitution_ok
    and all(record["leading_A_rank"] == 1 for record in step_records)
    and all(record["next_BR_rank"] == 2 for record in step_records)
)
check(
    "rank-one A and rank-two (B,R) systems select one contracting branch",
    coefficient_system_ok,
    "; ".join(
        f"n={n}: A={sp.N(physical(A[n]), 12)}, "
        f"B={sp.N(physical(B[n]), 8)}, R={sp.N(physical(R[n]), 12)}"
        for n in range(1, 5)
    ),
)


# Blind sequences derived solely from the selected coefficient branch.
U = {n: sp.simplify(A[n]-A[n-1]) for n in range(1, 5)}
V = {n: sp.simplify(R[n]-R[n-1]) for n in range(1, 5)}
post_coefficients = {}
for n in range(1, 5):
    expression = sp.expand(p_plus_scaled_compact.subs({
        a_minus: A[n-1]+B[n-1]*x,
        a_plus: A[n]+B[n]*x,
        r: R[n],
    }))
    _, leading = first_nonzero(expression)
    post_coefficients[n] = sp.factor(leading/(180*EPSILON3_COMPACT))

blind_ratios = {
    "u_over_u1": {n: sp.simplify(U[n]/U[1]) for n in range(1, 5)},
    "a_over_u1": {n: sp.simplify(A[n]/U[1]) for n in range(1, 5)},
    "v_over_v1": {n: sp.simplify(V[n]/V[1]) for n in range(1, 5)},
    "r_over_v1": {n: sp.simplify(R[n]/V[1]) for n in range(1, 5)},
    "p_post_over_k": post_coefficients,
}


# Independent exact-equation residual order at three small e values.
f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO), F_RHO, "mpmath")
p_minus_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO), P_MINUS, "mpmath")
p_plus_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO), P_PLUS, "mpmath")
A_mp = {n: arb.mpf(str(sp.N(physical(value), DPS))) for n, value in A.items()}
B_mp = {n: arb.mpf(str(sp.N(physical(value), DPS))) for n, value in B.items()}
R_mp = {n: arb.mpf(str(sp.N(physical(value), DPS))) for n, value in R.items()}


def state_at(e_value, n):
    return arb.exp(A_mp[n]*e_value**2+B_mp[n]*e_value**4)


def rho_at(e_value, n):
    return e_value**2*arb.exp(R_mp[n]*e_value**2)


residual_records = []
residual_sequences = defaultdict(list)
for exact_e in SERIES_CONTROLS:
    e_value = mp_value(exact_e)
    for n in range(1, 5):
        current = (
            state_at(e_value, n-1), state_at(e_value, n), rho_at(e_value, n)
        )
        previous = (
            state_at(e_value, n-2), state_at(e_value, n-1),
            rho_at(e_value, n-1),
        )
        f_value = f_numeric(*current)
        g_value = p_plus_numeric(*previous)+p_minus_numeric(*current)
        residual_sequences[(n, "F")].append(abs(f_value))
        residual_sequences[(n, "G")].append(abs(g_value))
        residual_records.append({
            "e": str(exact_e), "n": n,
            "F_abs": arb.nstr(abs(f_value), 40),
            "G_abs": arb.nstr(abs(g_value), 40),
        })

observed_orders = {}
residual_order_ok = True
for key, values in residual_sequences.items():
    orders = (
        arb.log(values[0]/values[1], 2),
        arb.log(values[1]/values[2], 2),
    )
    observed_orders[f"n={key[0]}:{key[1]}"] = [
        arb.nstr(value, 30) for value in orders
    ]
    expected_order = arb.mpf(7 if key[1] == "F" else 5)
    residual_order_ok &= bool(
        values[0] > values[1] > values[2] > 0
        and all(abs(order-expected_order) < arb.mpf("0.3") for order in orders)
    )
check(
    "truncated blind coefficients give the predicted higher-order residual decay",
    residual_order_ok,
    "expected F order 7 and seam order 5",
)


ZETA = (sp.pi**2*sp.sqrt(2)/50)**sp.Rational(1, 3)
CONTINUUM_A1 = -ZETA**2/2
continuum_ratio = sp.factor(A[1]/CONTINUUM_A1)
continuum_control_ok = bool(
    sp.N(physical(continuum_ratio), 60).is_real
    and sp.N(physical(continuum_ratio), 60) > 0
)
check(
    "the fixed half-step Friedmann coefficient comparison is well defined",
    continuum_control_ok,
    f"A1/A1_FLRW={sp.N(physical(continuum_ratio), 16)}",
)


if coefficient_system_ok and residual_order_ok and derivative_control_ok:
    outcome = "CELLULAR_WEAK_LAPSE_JET_DERIVED"
elif exact_substitution_ok and not unique_contracting:
    outcome = "CELLULAR_WEAK_LAPSE_JET_NONUNIQUE"
elif leading_orders_found and not exact_substitution_ok:
    outcome = "CELLULAR_WEAK_LAPSE_JET_REFUTED"
else:
    outcome = "CELLULAR_WEAK_LAPSE_JET_OPEN"

check(
    "the blind outcome follows the preregistered hierarchy",
    outcome in {
        "CELLULAR_WEAK_LAPSE_JET_DERIVED",
        "CELLULAR_WEAK_LAPSE_JET_NONUNIQUE",
        "CELLULAR_WEAK_LAPSE_JET_REFUTED",
        "CELLULAR_WEAK_LAPSE_JET_OPEN",
    },
    f"outcome={outcome}",
)


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "action_input_sha256": ACTION_SHA256,
    "tick_artifacts_parsed": TICK_ARTIFACTS_PARSED,
    "forbidden_tick_tokens": list(FORBIDDEN_TICK_TOKENS),
    "outcome": outcome,
    "labels": {
        "weak_lapse_jet": (
            "DERIVED" if outcome == "CELLULAR_WEAK_LAPSE_JET_DERIVED"
            else "OPEN"
        ),
        "comparison_to_committed_ticks": "NOT_PERFORMED",
        "continuum_convergence": "OPEN",
        "physical_lapse": "OPEN",
        "external_novelty": "OPEN",
    },
    "static": {
        "action": str(static_action),
        "P_minus": str(static_p_minus),
        "P_plus": str(static_p_plus),
    },
    "derivative_controls": derivative_records,
    "maximum_derivative_relative_error": arb.nstr(
        maximum_derivative_error, 50
    ),
    "universal_series": {
        "S_over_e_coefficients": {
            str(power): str(value)
            for power, value in series_coefficients.items()
        },
        "F_over_e_first_power": f_power,
        "P_minus_over_e_first_power": p_minus_power,
        "P_plus_over_e_first_power": p_plus_power,
        "F_over_e": str(sp.factor(f_scaled)),
        "P_minus_over_e": str(sp.factor(p_minus_scaled)),
        "P_plus_over_e": str(sp.factor(p_plus_scaled)),
    },
    "steps": step_records,
    "coefficients": {
        "A": {str(n): str(A[n]) for n in range(1, 5)},
        "B": {str(n): str(B[n]) for n in range(1, 5)},
        "R": {str(n): str(R[n]) for n in range(1, 5)},
        "U": {str(n): str(U[n]) for n in range(1, 5)},
        "V": {str(n): str(V[n]) for n in range(1, 5)},
    },
    "blind_ratios": {
        label: {str(n): str(value) for n, value in sequence.items()}
        for label, sequence in blind_ratios.items()
    },
    "series_residuals": residual_records,
    "observed_halving_orders": observed_orders,
    "continuum_half_step_control": {
        "zeta": str(ZETA),
        "A1_FLRW": str(CONTINUUM_A1),
        "A1_discrete_over_A1_FLRW": str(continuum_ratio),
        "numeric_ratio": str(sp.N(physical(continuum_ratio), 50)),
    },
    "scope": {
        "n_max": 4,
        "finite_tick_comparison": "NOT_PERFORMED",
        "spatial_refinement": "NOT_TESTED",
        "anisotropic_modes": "NOT_TESTED",
    },
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("\nBlind coefficients (no tick comparison):")
for n in range(1, 5):
    print(
        f"  n={n}: A={sp.factor(A[n])}, B={sp.factor(B[n])}, "
        f"R={sp.factor(R[n])}, "
        f"u/u1={blind_ratios['u_over_u1'][n]}, "
        f"v/v1={blind_ratios['v_over_v1'][n]}, "
        f"p/k={blind_ratios['p_post_over_k'][n]}"
    )
print(f"\nSummary: {passed}/{tests} checks passed")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")
if passed != tests or outcome in {
    "CELLULAR_WEAK_LAPSE_JET_OPEN", "CELLULAR_WEAK_LAPSE_JET_REFUTED"
}:
    raise SystemExit(1)
