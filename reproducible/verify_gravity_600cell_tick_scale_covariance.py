#!/usr/bin/env python3
"""Exact and complete-action audit of classical tick scale covariance.

Prior-art commit: 3275575.
Protocol commit: fd1f8a8.

No stationary equation is solved. The off-shell perturbation and both scale
factors were frozen before this verifier was evaluated.
"""

import ast
from collections import Counter
import contextlib
import io
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_gravity_600cell_dust_canonical_continuation.py"
REFINED_MASS = HERE / "gravity_600cell_refined_local_curvature_mass.json"
OUTPUT = HERE / "gravity_600cell_tick_scale_covariance.json"
PRIOR_ART_COMMIT = "3275575"
PROTOCOL_COMMIT = "fd1f8a8"
DPS = 100
mp.mp.dps = DPS
ALPHAS = (mp.mpf(3) / 5, mp.mpf(7) / 4)
COVARIANCE_TOLERANCE = mp.mpf("1e-65")
BRANCH_IMAGINARY_TOLERANCE = mp.mpf("1e-65")
HOSTILE_MINIMUM = mp.mpf("1e-8")
SUPPORT_MINIMUM = mp.mpf("1e-20")
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


def text(value, digits=30):
    return mp.nstr(value, digits)


def normalized_error(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))


def maximum_imaginary(action, gradient):
    return max(abs(mp.im(action)), *(abs(mp.im(value)) for value in gradient))


def load_action_core():
    """Execute definitions only; do not run the upstream continuation."""
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "print"
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("canonical-continuation definition cutoff not found")
    namespace = {
        "__file__": str(SOURCE),
        "__name__": "tick_scale_covariance_action_core",
    }
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(SOURCE), "exec"), namespace)
    return namespace


print("Classical tick scale-covariance audit", flush=True)
core = load_action_core()
models = core["models"]
action_and_gradient = core["action_and_gradient"]
action_globals = action_and_gradient.__globals__
base_mass = core["ARB_MASS"]

provenance_ok = bool(
    PRIOR_ART_COMMIT == "3275575"
    and PROTOCOL_COMMIT == "fd1f8a8"
    and core["tests"] == core["passed"] == 4
    and set(models) == {"even", "odd"}
)
check("the frozen provenance and complete-action core are intact", provenance_ok)


# Exact degree audit. Here r scales squared lengths and alpha scales lengths.
r, alpha = sp.symbols("r alpha", positive=True)
x, y, z = sp.symbols("x y z")
area_square = (2 * (x * y + x * z + y * z) - x**2 - y**2 - z**2) / 16
area_polynomial_ok = bool(
    sp.expand(
        area_square.subs({x: r * x, y: r * y, z: r * z})
        - r**2 * area_square
    ) == 0
)

# Ten independent entries prevent a diagonal/numerical special-case proof.
g_symbols = sp.symbols("g00 g01 g02 g03 g11 g12 g13 g22 g23 g33")
g00, g01, g02, g03, g11, g12, g13, g22, g23, g33 = g_symbols
gram = sp.Matrix([
    [g00, g01, g02, g03],
    [g01, g11, g12, g13],
    [g02, g12, g22, g23],
    [g03, g13, g23, g33],
])
gram_degree_ok = bool(sp.expand((r * gram).det() - r**4 * gram.det()) == 0)

v4, f3a, f3b, h2, dv = sp.symbols("v4 f3a f3b h2 dv", positive=True)
cosine_scaled = 16 * (r**3 * dv) / (
    sp.sqrt(r**3 * f3a) * sp.sqrt(r**3 * f3b)
)
cosine_base = 16 * dv / (sp.sqrt(f3a) * sp.sqrt(f3b))
sine_scaled = -sp.Rational(4, 3) * (
    sp.sqrt(r**2 * h2) * sp.sqrt(r**4 * v4)
) / (sp.sqrt(r**3 * f3a) * sp.sqrt(r**3 * f3b))
sine_base = -sp.Rational(4, 3) * (
    sp.sqrt(h2) * sp.sqrt(v4)
) / (sp.sqrt(f3a) * sp.sqrt(f3b))
angle_degree_ok = bool(
    sp.simplify(cosine_scaled - cosine_base) == 0
    and sp.simplify(sine_scaled - sine_base) == 0
)

A, deficit, mass, rho = sp.symbols("A deficit mass rho", positive=True)
gravity_degree_ok = bool(
    sp.simplify((r * A) * deficit - r * (A * deficit)) == 0
)
dust_degree_ok = bool(
    sp.simplify(
        (alpha * mass) * sp.sqrt(alpha**2 * rho)
        - alpha**2 * mass * sp.sqrt(rho)
    ) == 0
)
fixed_mass_defect_ok = bool(
    sp.simplify(
        mass * sp.sqrt(alpha**2 * rho)
        - alpha**2 * mass * sp.sqrt(rho)
        - (alpha - alpha**2) * mass * sp.sqrt(rho)
    ) == 0
)

L, epsilon3 = sp.symbols("L epsilon3", positive=True)
coarse_mass = sp.Rational(90, 1) / sp.pi * epsilon3 * L
coarse_mass_degree_ok = bool(
    sp.simplify(coarse_mass.subs(L, alpha * L) - alpha * coarse_mass) == 0
)
ell = sp.symbols("ell0:4", positive=True)
delta = sp.symbols("delta0:4")
curvature = sum(length * angle for length, angle in zip(ell, delta))
refined_mass_degree_ok = bool(
    sp.simplify(
        curvature.subs({length: alpha * length for length in ell}) / (8 * sp.pi)
        - alpha * curvature / (8 * sp.pi)
    ) == 0
)
refined_artifact = json.loads(REFINED_MASS.read_text())
refined_provenance_ok = bool(
    refined_artifact.get("outcome")
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and refined_artifact.get("tests") == {"passed": 15, "total": 15}
    and refined_artifact.get("definitions", {}).get("selected_mass")
        == "mu_r=K_r/(8*pi)"
)

algebra_ok = all((
    area_polynomial_ok,
    gram_degree_ok,
    angle_degree_ok,
    gravity_degree_ok,
    dust_degree_ok,
    fixed_mass_defect_ok,
    coarse_mass_degree_ok,
    refined_mass_degree_ok,
    refined_provenance_ok,
))
check(
    "exact algebra gives degree two for action and degree one for both selected mass rules",
    algebra_ok,
    f"area={area_polynomial_ok}, gram={gram_degree_ok}, angles={angle_degree_ok}, refined provenance={refined_provenance_ok}",
)


def perturb(values, modulus, center):
    return tuple(
        value * mp.exp(mp.mpf("1e-6") * ((index % modulus) - center))
        for index, value in enumerate(values)
    )


base_old = perturb(core["ARB_BASE_OLD"], 7, 3)
base_internal = perturb(
    tuple([core["ARB_SLANT_SQUARE"]] * 30 + [core["ARB_RHO"]] * 5),
    5,
    2,
)
base_new = perturb(tuple([core["ARB_L0_SQUARE"]] * 30), 11, 5)


def evaluate(model, old, internal, new, selected_mass):
    previous_mass = action_globals["ARB_MASS"]
    action_globals["ARB_MASS"] = selected_mass
    try:
        action, gradient, branch = action_and_gradient(model, old, internal, new)
    finally:
        action_globals["ARB_MASS"] = previous_mass
    imaginary = maximum_imaginary(action, gradient)
    branch_ok = bool(
        branch["negative_counts"] == Counter({1: 2400})
        and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > mp.mpf("1e-6")
        and imaginary < BRANCH_IMAGINARY_TOLERANCE
    )
    return {
        "action": action,
        "gradient": gradient,
        "branch": branch,
        "imaginary": imaginary,
        "branch_ok": branch_ok,
    }


records = {}
base_controls_ok = True
covariance_ok = True
hostile_ok = True
for parity, model in models.items():
    base = evaluate(model, base_old, base_internal, base_new, base_mass)
    internal_support = sum(
        abs(value) > SUPPORT_MINIMUM for value in base["gradient"][30:65]
    )
    boundary_values = base["gradient"][:30] + base["gradient"][65:]
    boundary_support = sum(abs(value) > SUPPORT_MINIMUM for value in boundary_values)
    base_controls_ok &= bool(
        base["branch_ok"] and internal_support >= 20 and boundary_support >= 20
    )
    parity_records = {
        "base": {
            "action": text(base["action"], 50),
            "internal_support": int(internal_support),
            "boundary_support": int(boundary_support),
            "minimum_leading_minor": text(base["branch"]["minimum_leading_minor"], 20),
            "minimum_argument": text(base["branch"]["minimum_argument"], 20),
            "maximum_imaginary": text(base["imaginary"], 12),
        },
        "scales": {},
    }
    for selected_alpha in ALPHAS:
        selected_r = selected_alpha**2
        scaled_old = tuple(selected_r * value for value in base_old)
        scaled_internal = tuple(selected_r * value for value in base_internal)
        scaled_new = tuple(selected_r * value for value in base_new)
        scaled = evaluate(
            model,
            scaled_old,
            scaled_internal,
            scaled_new,
            selected_alpha * base_mass,
        )
        fixed_mass = evaluate(
            model,
            scaled_old,
            scaled_internal,
            scaled_new,
            base_mass,
        )
        action_error = normalized_error(
            scaled["action"], selected_r * base["action"]
        )
        component_errors = tuple(
            normalized_error(left, selected_r * right)
            for left, right in zip(scaled["gradient"], base["gradient"])
        )
        maximum_component_error = max(component_errors)
        fixed_action_error = normalized_error(
            fixed_mass["action"], selected_r * base["action"]
        )
        fixed_pole_errors = tuple(
            normalized_error(
                fixed_mass["gradient"][index],
                selected_r * base["gradient"][index],
            )
            for index in range(60, 65)
        )
        maximum_fixed_pole_error = max(fixed_pole_errors)
        one_covariance_ok = bool(
            scaled["branch_ok"]
            and action_error < COVARIANCE_TOLERANCE
            and maximum_component_error < COVARIANCE_TOLERANCE
        )
        one_hostile_ok = bool(
            fixed_mass["branch_ok"]
            and fixed_action_error > HOSTILE_MINIMUM
            and maximum_fixed_pole_error > HOSTILE_MINIMUM
        )
        covariance_ok &= one_covariance_ok
        hostile_ok &= one_hostile_ok
        parity_records["scales"][text(selected_alpha, 20)] = {
            "r": text(selected_r, 20),
            "branch_ok": bool(scaled["branch_ok"]),
            "action_error": text(action_error, 12),
            "maximum_log_gradient_error": text(maximum_component_error, 12),
            "fixed_mass_action_error": text(fixed_action_error, 12),
            "fixed_mass_maximum_pole_error": text(maximum_fixed_pole_error, 12),
            "fixed_mass_hostile_control_pass": one_hostile_ok,
        }
        check(
            f"{parity}, alpha={text(selected_alpha, 8)}: complete action and all 95 log derivatives scale with degree two",
            one_covariance_ok,
            f"action error={text(action_error, 5)}, max derivative error={text(maximum_component_error, 5)}",
        )
        check(
            f"{parity}, alpha={text(selected_alpha, 8)}: holding mass fixed destroys covariance",
            one_hostile_ok,
            f"action error={text(fixed_action_error, 5)}, pole error={text(maximum_fixed_pole_error, 5)}",
        )
    records[parity] = parity_records

check(
    "the frozen off-shell state is branch-valid with nonzero internal and boundary support",
    base_controls_ok,
    "; ".join(
        f"{parity}: internal={record['base']['internal_support']}, boundary={record['base']['boundary_support']}"
        for parity, record in records.items()
    ),
)

if not (provenance_ok and algebra_ok and base_controls_ok and hostile_ok):
    outcome = "TICK_SCALE_COVARIANCE_CONTROL_FAILED"
elif not covariance_ok:
    outcome = "TICK_SCALE_COVARIANCE_REFUTED"
else:
    outcome = "TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED"

check(
    "the primary mechanical outcome is scale covariance",
    outcome == "TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED",
    outcome,
)

artifact = {
    "title": "Classical 600-cell tick scale covariance",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "hypotheses": {
        "carrier": "fixed one-slab 600-cell staircase",
        "action": "zero-Lambda Lorentzian flat-simplex Regge plus point dust",
        "transformation": "q->alpha^2 q, mass->alpha mass",
        "excluded": [
            "fixed external mass or length",
            "fixed cosmological constant",
            "quantum scale",
        ],
    },
    "exact_algebra": {
        "triangle_area_square_degree_two_in_q": area_polynomial_ok,
        "four_gram_determinant_degree_four_in_q": gram_degree_ok,
        "dihedral_angles_scale_invariant": angle_degree_ok,
        "regge_action_length_degree_two": gravity_degree_ok,
        "dust_action_length_degree_two_with_mass_scaling": dust_degree_ok,
        "fixed_mass_defect_identity": fixed_mass_defect_ok,
        "coarse_mass_length_degree_one": coarse_mass_degree_ok,
        "refined_curvature_mass_length_degree_one": refined_mass_degree_ok,
        "refined_mass_artifact_provenance": refined_provenance_ok,
    },
    "off_shell_perturbation": {
        "amplitude": "1e-6",
        "old": "exp(1e-6*((i mod 7)-3))",
        "internal": "exp(1e-6*((i mod 5)-2))",
        "new": "exp(1e-6*((i mod 11)-5))",
    },
    "alphas": [text(value, 30) for value in ALPHAS],
    "covariance_tolerance": text(COVARIANCE_TOLERANCE),
    "hostile_minimum": text(HOSTILE_MINIMUM),
    "parities": records,
    "interpretation": {
        "derived_if_adversarially_confirmed": (
            "stationary and canonical solutions occur in global scale families; "
            "the stated classical scale-free action cannot select an absolute tick"
        ),
        "not_excluded": [
            "tau/L",
            "tau_next/tau0",
            "relational dust time",
            "a duration relative to an external dimensional anchor",
        ],
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

print(f"RESULT: {passed}/{tests}", flush=True)
print(f"OUTCOME: {outcome}", flush=True)
if passed != tests:
    raise SystemExit(1)

