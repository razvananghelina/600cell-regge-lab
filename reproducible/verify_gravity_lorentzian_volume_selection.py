#!/usr/bin/env python3
"""Audit whether existing theory selects the Lorentzian tent volume term.

Protocol commit: dc7a6ab.  The preliminary non-selection result and the
optimistic spatial spectral formula were disclosed before implementation.
This is a current-repository selection audit, not a no-go for future actions.
"""

import ast
from collections import Counter
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_lorentzian_volume_selection.json"
PROTOCOL_COMMIT = "dc7a6ab"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def load(name):
    return json.loads((HERE / name).read_text())


def registered_scripts():
    tree = ast.parse((HERE / "run_all.py").read_text())
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "scripts"
               for target in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("could not read scripts registry")


print("=" * 78)
print("LORENTZIAN TENT VOLUME-COEFFICIENT SELECTION AUDIT")
print("=" * 78)

tent = load("gravity_lorentzian_tent.json")
metric_selector = load("hopf_spectral_metric_selector.json")
induced_gravity = load("hopf_kahler_induced_gravity.json")
spectral_sign = load("round_regge_spectral_action_sign.json")
slab = load("gravity_time_slab_canonicity.json")
box_hessian = load("gravity_box4_full_hessian.json")

check(
    "the input is the committed 23/23 Lorentzian tent no-go certificate",
    tent["verdict"] == "DERIVED SCOPED LORENTZIAN VACUUM NO-GO"
    and tent["tests"] == tent["passed"] == 23,
)
volume_control = tent["volume_control"]
target_ells = [
    volume_control["ell_zero_limit_float"],
    volume_control["ell_at_rho_one"],
    volume_control["ell_infinity_limit_float"],
]
check(
    "the required static volume coefficients are positive and unselected",
    all(value > 0 for value in target_ells)
    and volume_control["root_is_selected_without_lambda"] is False,
    "controls ell={:.12f}, {:.12f}, {:.12f}".format(*target_ells),
)

# -------------------------------------------------------------------------
# Fixed finite spectral data: exact but not continuum heat coefficients.
# -------------------------------------------------------------------------
finite = metric_selector["exact_results"]
check(
    "the fixed finite moments and ratio are the certified exact values",
    finite["finite_moments"] == [2640, 14880, 55920]
    and finite["finite_ratio"] == "31/11"
    and sp.Rational(14880, 2*2640) == sp.Rational(31, 11),
)
finite_census = metric_selector["licensed_census"]["finite_incidence_D"]
check(
    "the finite incidence operator has neither variable metric nor selected cutoff",
    finite_census["all_form"] is True
    and finite_census["variable_metric"] is False
    and finite_census["selected_cutoff"] is False,
)
dimension_note = (ROOT / "dimension_reconciliation.md").read_text()
discrete_note = (ROOT / "spectral_action_discrete_theorem.md").read_text()
check(
    "the authoritative dimension audit forbids a Seeley--DeWitt volume reading",
    "matrix moments, not Seeley--DeWitt coefficients" in dimension_note
    and "assigning them" in dimension_note
    and "cosmological, Einstein--Hilbert" in dimension_note
    and "These numbers are not Seeley--DeWitt coefficients" in discrete_note,
)

# r is a useful inverse-length unit but does not select the dimensionless heat
# time beta.  Under g->c^2 g, r->r/c^2 and t=beta/r; beta remains arbitrary.
c_scale, beta = sp.symbols("c_scale beta", positive=True)
r0 = sp.Rational(31, 11)
r_scaled = r0 / c_scale**2
t_scaled = beta / r_scaled
check(
    "31/11 supplies only a scale-covariant unit with a free heat parameter",
    sp.simplify(t_scaled/c_scale**2-beta/r0) == 0
    and finite["normalized_heat_family"]
    == "t=alpha/r, alpha>0 remains free",
    "t/c^2=beta/r0 for every beta>0; no certificate fixes beta",
)

# -------------------------------------------------------------------------
# Strongest optimistic source: the smooth three-dimensional de Rham action.
# -------------------------------------------------------------------------
exact_induced = induced_gravity["exact_results"]
check(
    "the continuum ordinary de Rham coefficients are A0=8Vol and A2=-2R/3",
    exact_induced["ordinary_derham_heat_a2_curvature_multiplier"] == "-2/3"
    and induced_gravity["hypotheses"]["carrier"] == "unit round S3=SU(2)",
    "A0=rank Lambda*(T*S3)=8 times Vol follows from ranks 1+3+3+1",
)
check(
    "the registered positive-cutoff action leaves its cutoff and scale open",
    spectral_sign["status"]
    == "DERIVED CONDITIONAL SIGN; COMPLETE FINITE ACTION OPEN"
    and "selection of chi and Lambda" in spectral_sign["open"]
    and spectral_sign["gaussian_family"]["C2_over_C0"] == "a",
)

Lambda, C0, C2, V3, R3 = sp.symbols(
    "Lambda C0 C2 V3 R3", positive=True
)
A0 = 8*V3
A2 = -sp.Rational(2, 3)*R3
spatial_action = Lambda**3*C0*A0 + Lambda*C2*A2
curvature_normalization = -sp.Rational(3, 2)/(Lambda*C2)
normalized_spatial_action = sp.expand(curvature_normalization*spatial_action)
lambda3 = sp.simplify(12*Lambda**2*C0/C2)
check(
    "normalizing the spatial curvature term derives lambda_3=12 Lambda^2 C0/C2",
    sp.simplify(normalized_spatial_action-(R3-lambda3*V3)) == 0,
    "the common (4*pi)^(-3/2) factor cancels",
)

v, s = sp.symbols("v s", positive=True)
chi = sp.exp(-s*v**2)
I0 = sp.integrate(chi, (v, 0, sp.oo))
I2 = sp.integrate(chi*v**2, (v, 0, sp.oo))
gaussian_C0 = sp.simplify(4*I2/sp.sqrt(sp.pi))
gaussian_C2 = sp.simplify(2*I0/sp.sqrt(sp.pi))
check(
    "the admissible Gaussian family has the exact registered Mellin moments",
    sp.simplify(gaussian_C0-s**(-sp.Rational(3, 2))) == 0
    and sp.simplify(gaussian_C2-s**(-sp.Rational(1, 2))) == 0,
)

edge_scale, ell = sp.symbols("edge_scale ell", positive=True)
ell_spatial = sp.simplify(
    12*(Lambda*edge_scale)**2*gaussian_C0/gaussian_C2
)
effective_cutoff = Lambda*edge_scale/sp.sqrt(s)
check(
    "the optimistic dimensionless coefficient is 12 times effective cutoff squared",
    sp.simplify(ell_spatial-12*effective_cutoff**2) == 0,
    "ell_3=12*(Lambda*a)^2/s",
)
inverse_cutoff = sp.sqrt(ell/12)
check(
    "the positive-cutoff family covers every positive coefficient exactly",
    sp.simplify(12*inverse_cutoff**2-ell) == 0
    and inverse_cutoff.is_positive,
    "for any ell>0 choose (Lambda*a)/sqrt(s)=sqrt(ell/12)",
)

# The finite internal unit does not remove the same freedom.  If the Gaussian
# heat time is beta/r, then mu^2=r/beta and ell=12r/beta.
ell_from_finite_unit = 12*r0/beta
check(
    "using r=31/11 still leaves a continuous ell=372/(11 beta) family",
    sp.simplify(ell_from_finite_unit-sp.Rational(372, 11)/beta) == 0
    and sp.limit(ell_from_finite_unit, beta, 0, dir="+") == sp.oo
    and sp.limit(ell_from_finite_unit, beta, sp.oo) == 0,
)

# Nor does the positive Gaussian finite heat trace extremize its own cutoff.
# Every nonzero eigenvalue term is strictly increasing with Lambda.
spectral_eigenvalue = sp.symbols("spectral_eigenvalue", positive=True)
gaussian_mode = sp.exp(-spectral_eigenvalue/Lambda**2)
check(
    "the Gaussian finite spectral action cannot select an interior cutoff by extremization",
    sp.diff(gaussian_mode, Lambda).is_positive,
    "each positive mode has d exp(-lambda/Lambda^2)/dLambda>0",
)

# -------------------------------------------------------------------------
# Carrier/bridge census.
# -------------------------------------------------------------------------
check(
    "the certified spectral curvature response is explicitly three-dimensional",
    induced_gravity["hypotheses"]["carrier"] == "unit round S3=SU(2)"
    and "Lorentzian time" in induced_gravity["hypotheses"]["not_assumed"]
    and any(
        item["label"] == "OPEN PHYSICAL SELECTION"
        for item in induced_gravity["verdicts"]
    ),
)
check(
    "the canonical four-dimensional slab still has no Lorentzian Regge action",
    slab["tests"] == slab["passed"] == 23
    and "DERIVED CURRENT SLAB-DYNAMICS GAP" in slab["verdicts"]
    and "Lorentzian edge data and a Regge action on the cylinder"
    in slab["open"],
)
check(
    "the finite Box Hessian is certified only as static stiffness",
    any(
        verdict["label"] == "OPEN"
        and "no graviton follows" in verdict["claim"]
        for verdict in box_hessian["verdicts"]
    ),
)

alpha_source = (HERE / "verify_alpha_spectral.py").read_text().lower()
legacy_gravity_source = (HERE / "verify_gravity.py").read_text().lower()
bridge_terms = (
    "tent", "slab", "four-volume", "cosmological volume",
    "cutoff-to-edge", "lambda*a^2",
)
check(
    "the registered alpha verifier defines no tent/slab volume bridge",
    all(term not in alpha_source for term in bridge_terms),
    "its alpha equation is a Box/gauge normalization, not a cutoff-to-edge map",
)
check(
    "the legacy registered gravity verifier defines no tent/slab volume bridge",
    all(term not in legacy_gravity_source for term in bridge_terms),
    "its outputs are graph curvature and static Box stiffness",
)

# Historical experiments are reported rather than silently ignored.  They
# are not registered verifiers, and the CC experiment itself states that the
# heat cutoff is an input.
scripts = registered_scripts()
experiment_cc = (ROOT / "exp513_cc_spectral_action.py").read_text()
experiment_old_action = (ROOT / "exp261_spectral_action.py").read_text()
check(
    "the historical cosmological experiment explicitly leaves its cutoff free",
    "Lambda is an INPUT, not an output." in experiment_cc
    and "This is NOT the cosmological constant." in experiment_cc
    and "exp513_cc_spectral_action.py" not in scripts,
)
check(
    "the old full-bosonic-action assertion is exploratory, unregistered, and superseded",
    "This IS the full bosonic Lagrangian" in experiment_old_action
    and "exp261_spectral_action.py" not in scripts
    and "matrix moments, not Seeley--DeWitt coefficients"
    in dimension_note,
)

# Six frozen acceptance fields.  A True means that the candidate actually
# supplies that field, not that the field can be guessed externally.
acceptance_fields = (
    "same_carrier_4d_action_or_proved_transfer",
    "curvature_and_four_volume_in_one_action",
    "relative_sign_and_normalization_selected",
    "unique_dimensionless_ell",
    "cutoff_function_and_scale_selected_or_not_applicable",
    "native_4d_or_derived_3plus1_transfer",
)
candidates = {
    "generic_lorentzian_tent_extension": (
        True, True, False, False, True, True,
    ),
    "spatial_derham_spectral_action": (
        False, False, False, False, False, False,
    ),
    "fixed_finite_incidence_spectral_moments": (
        False, False, False, False, False, False,
    ),
    "canonical_cw_slab": (
        False, False, False, False, True, True,
    ),
    "box_edge_stiffness": (
        False, False, False, False, True, False,
    ),
    "box_alpha_equation": (
        False, False, False, False, True, False,
    ),
}
passing_candidates = [
    name for name, fields in candidates.items() if all(fields)
]
check(
    "no frozen candidate supplies all six coefficient-selection fields",
    not passing_candidates
    and all(len(fields) == len(acceptance_fields)
            for fields in candidates.values()),
    "passes=0/{} candidates".format(len(candidates)),
)

# Only now compare with the already committed tent requirement.  Every
# positive control can be matched by the free spectral scale, which is
# anti-selection rather than prediction.
fitted_effective_cutoffs = [sp.sqrt(sp.Float(value, 17)/12)
                            for value in target_ells]
fit_residuals = [
    abs(float(12*cutoff**2)-value)
    for cutoff, value in zip(fitted_effective_cutoffs, target_ells)
]
check(
    "every displayed tent coefficient can be reproduced only by choosing the free cutoff",
    max(fit_residuals) < 2e-15,
    "effective cutoffs={}".format(
        [float(value) for value in fitted_effective_cutoffs]
    ),
)
check(
    "the spectral candidate fixes a favorable sign but not a value",
    C0.is_positive and C2.is_positive
    and not passing_candidates,
    "STRUCTURAL SIGN ONLY; magnitude and 3D-to-4D transfer remain absent",
)

verdict = (
    "DERIVED CURRENT VOLUME-SELECTION ABSENCE; STRUCTURAL SIGN ONLY"
    if passed == tests else "REFUTED OR INCOMPLETE"
)
result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "question": (
        "does the committed theory select a unique ell=lambda*a^2 for the "
        "same Lorentzian tent without consulting a desired rho?"
    ),
    "acceptance_fields": list(acceptance_fields),
    "candidate_census": {
        name: dict(zip(acceptance_fields, fields))
        for name, fields in candidates.items()
    },
    "candidate_count": len(candidates),
    "fully_passing_candidates": passing_candidates,
    "finite_spectral_data": {
        "moments": finite["finite_moments"],
        "r": finite["finite_ratio"],
        "heat_time": "t=beta/r with beta>0 unselected",
        "dimension_status": "finite Taylor moments, not Seeley--DeWitt coefficients",
    },
    "optimistic_spatial_spectral_transfer": {
        "carrier": "three-dimensional S3, not the Lorentzian tent",
        "A0": "8*Vol_3",
        "A2": "-(2/3)*integral R_3",
        "normalized_lambda3": "12*Lambda^2*C0/C2",
        "gaussian_family": "ell_3=12*(Lambda*a)^2/s",
        "range": "all positive real ell_3",
        "relative_sign": "favorable after curvature normalization",
        "status": "STRUCTURAL SIGN ONLY; NO 3+1 TRANSFER",
    },
    "tent_comparison_after_census": {
        "ell_controls": target_ells,
        "fitted_effective_cutoffs": [
            float(value) for value in fitted_effective_cutoffs
        ],
        "interpretation": "free-cutoff fit, not prediction",
    },
    "derived": [
        "the fixed finite moments do not define a four-dimensional volume coefficient",
        "the strongest spatial de Rham spectral action fixes only a favorable relative sign",
        "its admissible cutoff family spans every positive dimensionless coefficient",
        "the canonical four-dimensional slab has no selected Lorentzian operator or action",
        "none of the six frozen current candidate classes passes all acceptance fields",
    ],
    "structural": [
        "if the spatial A0/A2 ratio were transferred to the tent, its sign could support a positive lambda",
        "31/11 can serve as an internal inverse-length unit after a heat parameter is supplied",
    ],
    "open": [
        "a selected four-dimensional Lorentzian spectral/Regge action",
        "a cutoff function and cutoff-to-edge normalization",
        "a derived 3+1 transfer including kinetic/extrinsic-curvature terms",
        "nonsymmetric tent data, matter and higher-curvature alternatives",
    ],
    "not_claimed": [
        "no future theory can generate a cosmological term",
        "the three-dimensional spectral coefficient is a four-dimensional cosmological constant",
        "the required tent coefficient is experimentally excluded",
    ],
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

print("-" * 78)
print(f"VERDICT: {verdict}")
print(f"RESULT: {passed}/{tests} checks passed")
print(f"Machine-readable result: {OUTPUT}")
raise SystemExit(0 if passed == tests else 1)
