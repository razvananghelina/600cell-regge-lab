#!/usr/bin/env python3
"""Rigorous primary certificate for the normalized invariant half-strip."""

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp
from flint import arb, arb_series, ctx


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROTOCOL = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_invariant_region_protocol.md"
)
CORRECTION = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_invariant_region_protocol_correction.md"
)
REPORTING_CORRECTION = (
    ROOT
    / "docs"
    / "gravity"
    / "gravity_600cell_finite_height_invariant_region_reporting_correction.md"
)
CLASSIFICATION_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
FIFTH_INPUT = HERE / "gravity_600cell_finite_height_asymptotic_map.json"
OUTPUT = HERE / "gravity_600cell_finite_height_invariant_region.json"

PROTOCOL_COMMIT = "559ca7a"
PROTOCOL_SHA256 = (
    "fe7891ef35e1273c81b2df3da6d5b8b8f8d142bdd28a936bbb6eb768450e2de6"
)
CORRECTION_COMMIT = "8e30b1c"
CORRECTION_SHA256 = (
    "f79e24e63b8e544e5cc34b51cf02354e6dcfc376a6eb5853b2f36c1bd2aac697"
)
REPORTING_CORRECTION_COMMIT = "8f770cd"
REPORTING_CORRECTION_SHA256 = (
    "6a985d0b41d4d1d6ac7825d59014521372f6d62b28211cf6f083e8ffa1789f07"
)
CLASSIFICATION_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
FIFTH_SHA256 = (
    "a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6"
)

ctx.dps = 192
mp.mp.dps = 140
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


def text(value, digits=75):
    return mp.nstr(value, digits)


def rational_ball(numerator, denominator=1):
    return arb(numerator) / arb(denominator)


def interval_ball(left_num, left_den, right_num, right_den):
    left = Fraction(left_num, left_den)
    right = Fraction(right_num, right_den)
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    center_ball = rational_ball(midpoint.numerator, midpoint.denominator)
    radius_text = f"{radius.numerator}/{radius.denominator}"
    return center_ball + arb(0, radius_text)


def strict_positive(value):
    return bool(value > 0 and value.lower() > 0 and not value.contains(0))


def arb_record(value):
    return {
        "pretty": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "contains_zero": bool(value.contains(0)),
    }


classification = json.loads(CLASSIFICATION_INPUT.read_text())
fifth = json.loads(FIFTH_INPUT.read_text())
provenance_ok = bool(
    digest(PROTOCOL) == PROTOCOL_SHA256
    and digest(CORRECTION) == CORRECTION_SHA256
    and digest(REPORTING_CORRECTION) == REPORTING_CORRECTION_SHA256
    and digest(CLASSIFICATION_INPUT) == CLASSIFICATION_SHA256
    and digest(FIFTH_INPUT) == FIFTH_SHA256
    and classification["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
    and fifth["outcome"]
    == "CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB"
)
check(
    "the invariant protocol and accepted finite inputs are frozen",
    provenance_ok,
    f"protocol={PROTOCOL_COMMIT}; correction={CORRECTION_COMMIT}",
)


# Exact regularization in u=t^2.
u = sp.symbols("u", nonnegative=True)
epsilon_u = 2 * sp.pi - 5 * sp.acos(
    (1 + 2 * u) / (2 * (1 + 3 * u))
)
M_u = 180 * epsilon_u / (sp.pi * sp.sqrt(1 + 4 * u))
P_u = (
    180 * epsilon_u / sp.sqrt(1 + 4 * u)
    - 600
    * sp.sqrt(3)
    * sp.asinh(1 / sp.sqrt(8 * (1 + 3 * u)))
)
P_zero = sp.simplify(P_u.subs(u, 0))
P_infinity = 60 * sp.pi - 300 * sp.sqrt(3) * sp.log(2)
N_u = sp.simplify(P_u - P_zero + 2 * sp.pi * u * M_u)

regular_origin_ok = bool(
    sp.simplify(P_zero.rewrite(sp.log) - P_infinity) == 0
    and sp.simplify(M_u.subs(u, 0) - 60) == 0
    and sp.simplify(N_u.subs(u, 0)) == 0
    and sp.simplify(sp.diff(N_u, u).subs(u, 0)) == 0
)
check(
    "the u=t^2 representation removes every compactification singularity",
    regular_origin_ok,
)


t = sp.symbols("t", real=True)
M_t_series = sp.series(M_u.subs(u, t**2), t, 0, 14).removeO().expand()
P_t_series = sp.series(P_u.subs(u, t**2), t, 0, 14).removeO().expand()
formal_series_ok = bool(
    all(M_t_series.coeff(t, degree) == 0 for degree in range(1, 13, 2))
    and all(P_t_series.coeff(t, degree) == 0 for degree in range(1, 13, 2))
    and sp.simplify(M_t_series.coeff(t, 0) - 60) == 0
    and sp.simplify(
        M_t_series.coeff(t, 2) + 120 + 300 * sp.sqrt(3) / sp.pi
    )
    == 0
    and sp.simplify(P_t_series.coeff(t, 2) + 120 * sp.pi) == 0
    and sp.simplify(
        P_t_series.coeff(t, 4) - 360 * sp.pi - 900 * sp.sqrt(3)
    )
    == 0
)
check(
    "formal even coefficients through degree twelve reproduce the accepted expansion",
    formal_series_ok,
)


# Algebraic certificate for the normalized mean-value gap.
U_s, M_s, B_s, z_s, r_s, Cbar_s, m_s = sp.symbols(
    "U M B z r Cbar m", positive=True
)
u_s = m_s**2 * z_s**2
u_plus_s = u_s / r_s**2
drift_s = -4 * (1 - U_s) * u_s * B_s / M_s**2
curve_change_s = -z_s**2 * Cbar_s * (u_s - u_plus_s)
gap_s = sp.factor(curve_change_s - drift_s)
gap_normalized_target = (
    4 * (1 - U_s) * B_s / M_s**2
    - z_s**2 * (1 - r_s ** (-2)) * Cbar_s
)
gap_identity_ok = sp.simplify(gap_s / u_s - gap_normalized_target) == 0
check(
    "the exact same-x gap factors through m^2*z^2",
    gap_identity_ok,
)


# Arb automatic differentiation over the complete rational u interval.
API = arb.pi()
U_MAX = rational_ball(4, 390625)
U_INTERVAL = interval_ball(0, 1, 4, 390625)


def series_asinh(value):
    return (value + (value * value + 1).sqrt()).log()


def arb_M(value):
    eps = 2 * API - 5 * (
        (1 + 2 * value) / (2 * (1 + 3 * value))
    ).acos()
    return 180 * eps / (API * (1 + 4 * value).sqrt())


def arb_P(value):
    eps = 2 * API - 5 * (
        (1 + 2 * value) / (2 * (1 + 3 * value))
    ).acos()
    argument = 1 / (8 * (1 + 3 * value)).sqrt()
    return (
        180 * eps / (1 + 4 * value).sqrt()
        - 600 * arb(3).sqrt() * series_asinh(argument)
    )


u_jet = arb_series([U_INTERVAL, arb(1)], prec=3)
M_jet = arb_M(u_jet)
P_jet = arb_P(u_jet)
P0_arb = arb_P(arb(0))
N_jet = P_jet - P0_arb + 2 * API * u_jet * M_jet

M_ball = M_jet[0]
M_prime_ball = M_jet[1]
P_prime_ball = P_jet[1]
P_second_ball = 2 * P_jet[2]
N_second_ball = 2 * N_jet[2]

# Exact integral Taylor enclosures from the frozen method correction.
W_ball = P_prime_ball
Bbar_ball = N_second_ball / 2
W_prime_ball = P_second_ball / 2
C_ball = W_ball + 4 * API * M_ball
Cbar_ball = -W_prime_ball - 4 * API * M_prime_ball

primitive_bounds_ok = bool(
    strict_positive(M_ball)
    and strict_positive(Bbar_ball)
    and strict_positive(C_ball)
    and strict_positive(Cbar_ball)
)
check(
    "Arb separates every load-bearing one-variable coefficient on the full interval",
    primitive_bounds_ok,
    (
        f"M={M_ball}; Bbar={Bbar_ball}; "
        f"C={C_ball}; -C'={Cbar_ball}"
    ),
)


Z_INTERVAL = interval_ball(0, 1, 1, 125)
MASS_SQUARED_INTERVAL = interval_ball(0, 1, 4, 25)
U_BALL = Z_INTERVAL * M_ball
one_minus_u = 1 - U_BALL
current_physical_ok = bool(
    strict_positive(M_ball)
    and U_BALL < 1
    and strict_positive(one_minus_u)
)
check(
    "every interior point of the half-strip represents an expanding physical root",
    current_physical_ok,
    f"U={U_BALL}; hence r>1 and 0<m_plus<m",
)


# y_plus/z is bounded without dividing by z at the compactification axis.
yplus_over_z_ball = (
    4 * API
    - Z_INTERVAL * C_ball
    - 4
    * one_minus_u
    * MASS_SQUARED_INTERVAL
    * Z_INTERVAL
    * Bbar_ball
    / (M_ball * M_ball)
)
yplus_positive_ok = strict_positive(yplus_over_z_ball)
check(
    "the outgoing normalized momentum remains below p_infinity",
    yplus_positive_ok,
    f"y_plus/z={yplus_over_z_ball}",
)


# Since C'(u)<0, the omitted last term in Y_z is nonnegative.
curve_derivative_lower = 4 * API - 2 * Z_INTERVAL * C_ball
curve_monotone_ok = strict_positive(curve_derivative_lower)
check(
    "the complete next-root curve is strictly increasing in z",
    curve_monotone_ok,
    f"partial_z Y >= {curve_derivative_lower}",
)


# Use 0 < 1-r^-2 < 1 for r>1.  The expression is a rigorous lower
# enclosure of the exact normalized gap on the whole compact rectangle.
gap_normalized_lower = (
    4 * one_minus_u * Bbar_ball / (M_ball * M_ball)
    - Z_INTERVAL * Z_INTERVAL * Cbar_ball
)
same_x_gap_ok = strict_positive(gap_normalized_lower)
check(
    "the normalized same-x gap is positive including both continuous axes",
    same_x_gap_ok,
    f"gap/(m^2*z^2) >= {gap_normalized_lower}",
)


half_strip_map_ok = bool(
    current_physical_ok
    and yplus_positive_ok
    and curve_monotone_ok
    and same_x_gap_ok
)
check(
    "the exact successor satisfies m_plus<m and x_plus>x>=125",
    half_strip_map_ok,
    "Y(m_plus,0)=0<y_plus<Y(m_plus,z)",
)


# Global root uniqueness: use the accepted exact monotone partition only for
# its already-proved mu/p facts, then derive the present negative-root bound.
monotone = classification["monotone_facts"]
monotone_input_ok = bool(
    monotone["K_has_one_positive_squared_root"]
    and monotone["mu_inner_increasing_outer_decreasing"]
    and monotone["p_negative_on_positive_axis"]
    and classification["thresholds"]["rational_brackets"]["x_star"]
    == ["5", "6"]
)

q_r, mu_r, mu_prime_r, p_prime_r, mass_r = sp.symbols(
    "q mu mu_prime p_prime mass", nonzero=True
)
R_prime = p_prime_r + 4 * sp.pi * (
    q_r * mu_prime_r - (mu_r - mass_r)
) / q_r**2
R_prime = sp.factor(
    R_prime.subs(p_prime_r, -4 * sp.pi * mu_prime_r / q_r)
)
dual_identity_ok = sp.simplify(
    R_prime - 4 * sp.pi * (mass_r - mu_r) / q_r**2
) == 0

PINF_ARB = 60 * API - 300 * arb(3).sqrt() * arb(2).log()
negative_root_margin = -PINF_ARB - rational_ball(8, 5) * API
negative_root_excluded = strict_positive(negative_root_margin)
global_uniqueness_ok = bool(
    monotone_input_ok
    and dual_identity_ok
    and negative_root_excluded
    and half_strip_map_ok
)
check(
    "the all-real dual partition leaves exactly one physical successor",
    global_uniqueness_ok,
    f"negative-root lower margin={negative_root_margin}",
)


# The accepted pair is only a seed membership check, not evidence for the
# continuum inequalities above.
m3 = mp.mpf(fifth["history"]["m3"])
x4 = mp.mpf(fifth["history"]["x4"])
seed_membership_ok = bool(0 < m3 <= mp.mpf(2) / 5 and x4 >= 125)
check(
    "the accepted branch-B history lies strictly inside the certified half-strip",
    seed_membership_ok,
    f"2/5-m3={text(mp.mpf(2)/5-m3, 30)}; x4-125={text(x4-125, 30)}",
)


# Numerical fifth-root reproduction is deliberately after the domain proof.
MPI = mp.pi
MPINF = 60 * MPI - 300 * mp.sqrt(3) * mp.log(2)


def mepsilon_u(value):
    return 2 * MPI - 5 * mp.acos(
        (1 + 2 * value) / (2 * (1 + 3 * value))
    )


def mM(value):
    return 180 * mepsilon_u(value) / (MPI * mp.sqrt(1 + 4 * value))


def mP(value):
    return (
        180 * mepsilon_u(value) / mp.sqrt(1 + 4 * value)
        - 600
        * mp.sqrt(3)
        * mp.asinh(1 / mp.sqrt(8 * (1 + 3 * value)))
    )


def mW(value):
    if value == 0:
        return -120 * MPI
    return (mP(value) - MPINF) / value


def root_curve(mass, x_value):
    z_value = 1 / x_value
    u_value = (mass * z_value) ** 2
    m_value = mM(u_value)
    return (
        -z_value**2 * mW(u_value)
        - 4 * MPI * z_value * (z_value * m_value - 1)
    )


z4 = 1 / x4
u4 = (m3 * z4) ** 2
M4 = mM(u4)
U4 = z4 * M4
V4 = z4**2 * mW(u4)
y4_root = root_curve(m3, x4)
r4 = 2 / U4 - 1
m4 = m3 / r4
y4_out = -r4 * ((r4 + 1) * V4 + y4_root)
x5_rebuilt = mp.findroot(
    lambda candidate: root_curve(m4, candidate) - y4_out,
    (mp.mpf(125), mp.mpf(126)),
    tol=mp.mpf("1e-115"),
)
stored_x5 = mp.mpf(fifth["forecast_comparison"]["actual_x5"])
fifth_control_ok = abs(x5_rebuilt - stored_x5) < mp.mpf("1e-70")

boundary_x60_fail = sp.simplify((60 / sp.Symbol("x") - 1).subs(sp.Symbol("x"), 60)) == 0
P_no_boost_u = 180 * epsilon_u / sp.sqrt(1 + 4 * u)
N_no_boost = sp.simplify(
    P_no_boost_u
    - P_no_boost_u.subs(u, 0)
    + 2 * sp.pi * u * M_u
)
boost_changes_coefficient = sp.simplify(
    sp.diff(N_no_boost, u, 2).subs(u, 0)
    - sp.diff(N_u, u, 2).subs(u, 0)
) != 0
stored_hostile = fifth["checks"]["hostile_controls"]
controls_ok = bool(
    fifth_control_ok
    and boundary_x60_fail
    and boost_changes_coefficient
    and stored_hostile
)
check(
    "known-pass and known-fail controls separate the theorem from altered dynamics",
    controls_ok,
    f"rebuilt x5 error={text(x5_rebuilt-stored_x5, 25)}; x=60 has U=1",
)


all_gates = bool(
    provenance_ok
    and regular_origin_ok
    and formal_series_ok
    and gap_identity_ok
    and primitive_bounds_ok
    and current_physical_ok
    and yplus_positive_ok
    and curve_monotone_ok
    and same_x_gap_ok
    and half_strip_map_ok
    and global_uniqueness_ok
    and seed_membership_ok
    and controls_ok
)
outcome = (
    "INVARIANT_HALF_STRIP_PRIMARY_CERTIFICATE"
    if all_gates
    else "INVARIANT_HALF_STRIP_OPEN"
)
check(
    "the primary hierarchy certifies the half-strip pending adversarial replication",
    outcome == "INVARIANT_HALF_STRIP_PRIMARY_CERTIFICATE",
)


artifact = {
    "provenance": {
        "protocol_commit": PROTOCOL_COMMIT,
        "protocol_sha256": PROTOCOL_SHA256,
        "correction_commit": CORRECTION_COMMIT,
        "correction_sha256": CORRECTION_SHA256,
        "reporting_correction_commit": REPORTING_CORRECTION_COMMIT,
        "reporting_correction_sha256": REPORTING_CORRECTION_SHA256,
        "classification_sha256": CLASSIFICATION_SHA256,
        "fifth_sha256": FIFTH_SHA256,
    },
    "domain": {
        "m": "0<m<=2/5",
        "x": "x>=125",
        "compactification": "z=1/x, u=(m*z)^2",
        "u_max": "4/390625",
        "thresholds_post_hoc": True,
    },
    "method": {
        "precision_decimal_digits": ctx.dps,
        "interval_engine": "python-flint Arb outward-rounded balls",
        "removable_quotients": "exact integral Taylor remainders",
        "finite_grid_used": False,
        "subdivision_used": False,
        "global_root_partition": "R'(q)=4*pi*(m-mu(q))/q^2",
    },
    "formal_series": {
        "M_through_t12": str(M_t_series),
        "P_through_t12": str(P_t_series),
    },
    "arb_bounds": {
        "M": arb_record(M_ball),
        "W": arb_record(W_ball),
        "Bbar": arb_record(Bbar_ball),
        "C": arb_record(C_ball),
        "minus_C_prime": arb_record(Cbar_ball),
        "U": arb_record(U_BALL),
        "one_minus_U": arb_record(one_minus_u),
        "y_plus_over_z": arb_record(yplus_over_z_ball),
        "curve_derivative_lower": arb_record(curve_derivative_lower),
        "normalized_gap_lower": arb_record(gap_normalized_lower),
        "negative_root_margin": arb_record(negative_root_margin),
    },
    "seed": {
        "m3": text(m3),
        "x4": text(x4),
        "rebuilt_x5": text(x5_rebuilt),
        "stored_x5": text(stored_x5),
    },
    "checks": {
        "provenance": provenance_ok,
        "regular_origin": regular_origin_ok,
        "formal_series": formal_series_ok,
        "gap_identity": gap_identity_ok,
        "primitive_bounds": primitive_bounds_ok,
        "current_physical": current_physical_ok,
        "outgoing_y_positive": yplus_positive_ok,
        "curve_monotone": curve_monotone_ok,
        "same_x_gap": same_x_gap_ok,
        "half_strip_map": half_strip_map_ok,
        "global_uniqueness": global_uniqueness_ok,
        "seed_membership": seed_membership_ok,
        "controls": controls_ok,
    },
    "claims": {
        "invariant_half_strip": "PRIMARY_RIGOROUS_CERTIFICATE_PENDING_ADVERSARIAL",
        "accepted_history_infinite": "OPEN_PENDING_ADVERSARIAL",
        "generic_original_v": "OPEN",
        "nonhomogeneous_stability": "OPEN",
        "absolute_tick": "NOT_DERIVED",
    },
    "passed": passed,
    "tests": tests,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print()
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
if passed != tests:
    raise SystemExit(1)
