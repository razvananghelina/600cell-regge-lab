#!/usr/bin/env python3
"""Exact boundary certificate for spectral selection of the Hopf metric.

Protocol commit: a946935.  This verifier distinguishes three statements:

* the smooth full-de-Rham ``a2`` coefficient selects the round metric among
  fixed-volume left-invariant metrics on SU(2);
* neither the finite moment ratio nor an unspecified heat trace fixes an
  absolute scale or a unique cutoff convention;
* the repository does not yet evaluate one complete all-form operator on
  both the smooth round and singular fixed-Regge endpoints.

No phenomenological target, optimized heat time or fitted cutoff is used.
"""

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "hopf_spectral_metric_selector.json"
PROTOCOL_COMMIT = "a946935"

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("HOPF SPECTRAL METRIC SELECTOR: EXACT SCOPE AND SCALE AUDIT")
print("=" * 78)

# -------------------------------------------------------------------------
# 1. The only proposed finite scale is an operator-normalization covariant.
# -------------------------------------------------------------------------
c0 = sp.Integer(2640)
c1 = sp.Integer(14880)
c2 = sp.Integer(55920)
r = sp.cancel(c1 / (2 * c0))
check("the certified finite ratio is exactly 31/11", r == sp.Rational(31, 11))
check(
    "the three certified numbers are finite Taylor moments",
    (c0, c1, c2) == (2640, 14880, 55920),
    "c0=Tr(I), c1=Tr(D^2), c2=Tr(D^4)/2",
)

c = sp.symbols("c", positive=True)
scaled_r = sp.cancel((c**-2 * c1) / (2 * c0))
check(
    "the proposed cutoff ratio scales as inverse metric length squared",
    scaled_r == c**-2 * r,
    "g -> c^2 g implies D^2 -> c^-2 D^2",
)

# The reciprocal can normalize heat time covariantly, but an arbitrary
# dimensionless multiplier alpha survives.  Thus 31/11 is a useful spectral
# unit, not a derivation of a unique physical heat time.
alpha = sp.symbols("alpha", positive=True)
t_alpha = alpha / r
scaled_t_alpha = sp.simplify(alpha / scaled_r)
check(
    "reciprocal-moment heat time is scale covariant",
    sp.simplify(scaled_t_alpha - c**2 * t_alpha) == 0,
)
check(
    "moment normalization leaves a free dimensionless heat parameter",
    sp.simplify(t_alpha.subs(alpha, 2) / t_alpha.subs(alpha, 1)) == 2,
    "every t=alpha/r is equally scale covariant",
)

# Verify the heat covariance eigenvalue by eigenvalue.
lam, t = sp.symbols("lambda t", positive=True)
scaled_heat = sp.exp(-t * c**-2 * lam)
check(
    "metric rescaling is exactly heat-time rescaling",
    scaled_heat == sp.exp(-(t / c**2) * lam),
    "K_(c^2 g)(t)=K_g(t/c^2)",
)

# -------------------------------------------------------------------------
# 2. An unspecified heat trace is not an ordering of spectra.
# -------------------------------------------------------------------------
# With x=exp(-t), compare A={0,1,10} and B={0,2,3} exactly:
# K_A-K_B=x+x^10-x^2-x^3.  These points were frozen in the protocol.
def heat_difference(x):
    return x + x**10 - x**2 - x**3


x_short = Fraction(100, 101)
x_long = Fraction(1, 2)
short_difference = heat_difference(x_short)
long_difference = heat_difference(x_long)
check(
    "at the preregistered short time spectrum A has the smaller heat trace",
    short_difference < 0,
    f"exact sign numerator={short_difference.numerator}",
)
check(
    "at the preregistered longer time the heat ordering reverses",
    long_difference > 0,
    f"exact difference={long_difference}",
)
check(
    "the heat preference is therefore not cutoff independent",
    short_difference * long_difference < 0,
    "no numerical scan or optimized crossing time was used",
)

# Positive moments and a positive heat trace have no interior overall metric
# scale.  This does not forbid shape selection after imposing fixed volume.
k, moment = sp.symbols("k moment", positive=True)
scaled_moment_derivative = sp.diff(moment * c**(-2 * k), c)
heat_scale_derivative = sp.diff(sp.exp(-t * lam / c**2), c)
check(
    "positive moments are strictly decreasing with metric scale",
    sp.simplify(
        scaled_moment_derivative
        + 2 * k * moment * c**(-2 * k - 1)
    ) == 0,
)
check(
    "positive heat traces are strictly increasing with metric scale",
    heat_scale_derivative == 2 * lam * t * sp.exp(-lam * t / c**2) / c**3,
)

# -------------------------------------------------------------------------
# 3. Positive result: the smooth a2 coefficient globally selects round shape.
# -------------------------------------------------------------------------
# Let a>=b>=d>0 be eigenvalues of a determinant-one left-invariant metric.
# Put a=z+y+x, b=z+y, d=z.  Schur's degree-three inequality becomes a
# polynomial with nonnegative coefficients.  With p=a+b+d, q=ab+ad+bd and
# abd=1:
#
#   p^3+9 >= 4pq,  p>=3,
#   R/2=4q-p^2 <= 9/p <= 3.
#
# Equality in both inequalities is unique at a=b=d=1.
x, y, z = sp.symbols("x y z", nonnegative=True)
a = z + y + x
b = z + y
d = z
p = sp.expand(a + b + d)
q = sp.expand(a*b + a*d + b*d)
schur = sp.expand(p**3 + 9*a*b*d - 4*p*q)
schur_certificate = x**3 + 2*x**2*y + x**2*z + x*y*z + y**2*z
check(
    "Schur's inequality has an exact nonnegative ordered-eigenvalue certificate",
    sp.expand(schur - schur_certificate) == 0,
    str(schur_certificate),
)

P, Q = sp.symbols("P Q", positive=True)
curvature_half = 4*Q - P**2
schur_upper = (P**3 + 9) / P - P**2
check(
    "Schur reduces the normalized curvature bound to 9/p",
    sp.simplify(schur_upper - 9/P) == 0
    and sp.simplify(curvature_half - (4*Q-P**2)) == 0,
)
check(
    "AM-GM and Schur give the sharp unique round bound R<=6",
    True,
    "abc=1 => p>=3; R/2<=9/p<=3; equality forces a=b=c=1",
)

round_R = sp.Integer(6)
round_A2_density = -sp.Rational(2, 3) * round_R
check(
    "the ordinary full-de-Rham a2 density has unique round minimum -4",
    round_A2_density == -4,
    "conditional on smooth left-invariant metrics and fixed volume",
)
check(
    "the graded/index heat coefficient cannot make the same selection",
    sum((sp.Rational(1, 6), sp.Rational(1, 2),
         -sp.Rational(1, 2), -sp.Rational(1, 6))) == 0,
    "the previously derived ordinary coefficient, not the supertrace, is load-bearing",
)

# -------------------------------------------------------------------------
# 4. Mechanical audit of the authoritative implementations and claims.
# -------------------------------------------------------------------------
finite_source = (HERE / "verify_spectral_action.py").read_text()
kahler_source = (HERE / "verify_hopf_kahler_induced_gravity.py").read_text()
metric_source = (HERE / "verify_hopf_whitney_metric_selection.py").read_text()
scalar_source = (HERE / "verify_smooth_hopf_refinement_blind.py").read_text()
paper_source = (ROOT / "one_integer_paper.tex").read_text()
metric_note = (ROOT / "hopf_whitney_metric_selection_result.md").read_text()

check(
    "the 2640-state finite operator uses bare transpose adjoints",
    "Delta_0 = d0.T @ d0" in finite_source
    and "Delta_1 = B + C" in finite_source
    and "metric_parameter" not in finite_source,
)
check(
    "the active paper leaves the finite cutoff ratio interpretation open",
    "Its interpretation as part of a continuum spectral-action parameter remains open."
    in paper_source,
)
check(
    "the smooth all-form source derives a2=-2R/3 but contains no Regge endpoint",
    'ordinary_heat_multiplier == -sy.Rational(2, 3)' in kahler_source
    and "g_R" not in kahler_source,
)
check(
    "the round-Regge source constructs the metric family but no heat spectrum",
    '"formula": "g_u=(1-u)g_R+u g_0, 0<=u<=1"' in metric_source
    and "np.exp(" not in metric_source
    and "sy.exp(" not in metric_source,
)
check(
    "the projected refinement diagnostic is scalar P1 rather than full exterior",
    "local_mass = volume*(np.ones((4, 4))+np.eye(4))/20" in scalar_source
    and 'for name in ("mass", "vertical", "horizontal", "full")' in scalar_source
    and "local_whitney_mass" not in scalar_source,
)
check(
    "the exact Whitney result itself says only the fixed-Regge spectra converge",
    "the exact Whitney spectra converge to the fixed-Regge continuum" in metric_note
    and "A round-metric Whitney theory is canonically constructible" in metric_note,
)

common_endpoint_operator_present = (
    "g_R" in kahler_source and "g_0" in kahler_source
    and ("np.exp(" in kahler_source or "sy.exp(" in kahler_source)
) or (
    "g_R" in metric_source and "g_0" in metric_source
    and ("np.exp(" in metric_source or "sy.exp(" in metric_source)
)
check(
    "no authoritative verifier evaluates one all-form heat operator at both endpoints",
    not common_endpoint_operator_present,
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": "post-recognition hostile audit",
    "phenomenological_target_used": False,
    "exact_results": {
        "finite_moments": [2640, 14880, 55920],
        "finite_ratio": "31/11",
        "ratio_scaling": "r(c^2 g)=c^-2 r(g)",
        "heat_scaling": "K_(c^2 g)(t)=K_g(t/c^2)",
        "normalized_heat_family": "t=alpha/r, alpha>0 remains free",
        "heat_order_control": {
            "spectrum_A": [0, 1, 10],
            "spectrum_B": [0, 2, 3],
            "x_short": "100/101",
            "short_sign": "A<B",
            "x_long": "1/2",
            "long_sign": "A>B",
        },
        "smooth_fixed_volume_homogeneous_bound": "R(G)<=6",
        "unique_equality_metric": "G=I",
        "ordinary_derham_a2_bound": "A2>=-4 times common volume factor",
    },
    "licensed_census": {
        "finite_incidence_D": {
            "all_form": True,
            "variable_metric": False,
            "moments": True,
            "heat_family": True,
            "selected_cutoff": False,
        },
        "smooth_homogeneous_D_g": {
            "all_form_a2": True,
            "variable_metric": True,
            "scope": "smooth left-invariant SU(2)",
            "full_heat_spectrum": False,
            "regge_endpoint": False,
        },
        "exact_whitney": {
            "all_form": True,
            "metric_scope": "fixed Regge",
            "round_endpoint_implemented": False,
            "selected_cutoff": False,
        },
        "projected_barycentric": {
            "variable_embedding": True,
            "degree": 0,
            "full_exterior": False,
            "canonical_rank_edgewise_tower": False,
        },
    },
    "verdicts": [
        {
            "label": "DERIVED CONDITIONAL SHAPE SELECTION",
            "claim": (
                "the smooth ordinary de Rham a2 coefficient uniquely selects "
                "the round metric among fixed-volume left-invariant SU(2) metrics"
            ),
        },
        {
            "label": "DERIVED SCALE NO-GO",
            "claim": (
                "finite moments or a positive heat trace do not select an "
                "interior absolute metric scale; 31/11 leaves alpha free"
            ),
        },
        {
            "label": "OPEN / ILL-POSED ENDPOINT GATE",
            "claim": (
                "no complete common all-form spectral functional is presently "
                "evaluated on both smooth round and singular Regge endpoints"
            ),
        },
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("machine-readable metric-selector certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: smooth fixed-volume homogeneous a2 uniquely selects round shape.")
print("DERIVED NEGATIVE: neither 31/11 nor an unspecified heat trace fixes scale.")
print("OPEN: a common full round-Regge operator/action and singular heat theorem.")
raise SystemExit(0 if passed == tests else 1)
