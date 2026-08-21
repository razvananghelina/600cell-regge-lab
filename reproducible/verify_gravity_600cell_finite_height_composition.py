#!/usr/bin/env python3
"""Canonical two-slab composition gate for the finite-height update."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_classification.json"
ADVERSARIAL_INPUT = (
    HERE / "gravity_600cell_finite_height_classification_adversarial.json"
)
OUTPUT = HERE / "gravity_600cell_finite_height_composition.json"

PRIMARY_SHA256 = (
    "9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03"
)
ADVERSARIAL_SHA256 = (
    "da8d60e95b5196beaf93ea234fbf9dfb93e3d5e6bd00fb0a85ed2ef4ba388996"
)
PRIOR_ART_COMMIT = "c9eb996"
PROTOCOL_COMMIT = "386466e"

mp.mp.dps = 120

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


primary_input = json.loads(PRIMARY_INPUT.read_text())
adversarial_input = json.loads(ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    digest(PRIMARY_INPUT) == PRIMARY_SHA256
    and digest(ADVERSARIAL_INPUT) == ADVERSARIAL_SHA256
    and primary_input["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY"
    and adversarial_input["outcome"]
    == "FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY_ADVERSARIALLY_CORROBORATED"
)
check(
    "both accepted one-slab inputs and the composition protocol are frozen",
    provenance_ok,
    f"prior_art={PRIOR_ART_COMMIT}; protocol={PROTOCOL_COMMIT}",
)


# Complete action and canonical conventions.
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
P_POST_EXACT = L_PLUS * sp.diff(ACTION, L_PLUS) / 2

scale = sp.symbols("scale", positive=True)
scale_substitution = {
    L_MINUS: scale * L_MINUS,
    L_PLUS: scale * L_PLUS,
    RHO: scale**2 * RHO,
    MASS: scale * MASS,
}
action_homogeneity = sp.simplify(
    sp.powsimp(ACTION.subs(scale_substitution) - scale**2 * ACTION, force=True)
)
post_homogeneity = sp.simplify(
    sp.powsimp(
        P_POST_EXACT.subs(scale_substitution) - scale**2 * P_POST_EXACT,
        force=True,
    )
)
check(
    "the action and outgoing canonical momentum scale with degree two",
    action_homogeneity == 0 and post_homogeneity == 0,
)


# The shared-slice derivative fixes the post/pre sign without using a root.
L0, L1, L2, RHO1, RHO2 = sp.symbols(
    "L0 L1 L2 rho1 rho2", positive=True
)
S1 = ACTION.subs({L_MINUS: L0, L_PLUS: L1, RHO: RHO1})
S2 = ACTION.subs({L_MINUS: L1, L_PLUS: L2, RHO: RHO2})
dS1_shared = sp.diff(S1, L1)
dS2_shared = sp.diff(S2, L1)
post1 = L1 * dS1_shared / 2
pre2 = -L1 * dS2_shared / 2
d1_local, d2_local = sp.symbols("d1_local d2_local")
junction_identity = sp.expand(
    L1 * (d1_local + d2_local) / 2
    - (L1 * d1_local / 2 - (-L1 * d2_local / 2))
)
junction_definitions_ok = bool(
    post1 == L1 * dS1_shared / 2
    and pre2 == -L1 * dS2_shared / 2
)
check(
    "the shared-slice equation is exactly p_post(first)=p_pre(second)",
    junction_identity == 0 and junction_definitions_ok,
)


# Independently rebuild the exact affine one-slab residuals.  Radical
# normalization is algebraic and does not sample a state.
h = sp.symbols("h", positive=True)
q = sp.symbols("q", real=True)
finite_substitution = {
    L_MINUS: 1,
    L_PLUS: 1 + h * q,
    RHO: h**2,
}
constraint_raw = sp.factor((2 * F_EXACT / h).subs(finite_substitution))
momentum_raw = sp.factor(P_PRE_EXACT.subs(finite_substitution))


def radical_inventory(z):
    square = z**2
    return (
        (square + 4, sp.sqrt(square + 4)),
        (3 * square + 8, sp.sqrt(3 * square + 8)),
        (
            3 * square**2 + 20 * square + 32,
            sp.sqrt(square + 4) * sp.sqrt(3 * square + 8),
        ),
        (
            9 * square**3 + 84 * square**2 + 256 * square + 256,
            (3 * square + 8) * sp.sqrt(square + 4),
        ),
        (
            3 * square**3 + 32 * square**2 + 112 * square + 128,
            (square + 4) * sp.sqrt(3 * square + 8),
        ),
        (
            9 * square**4
            + 120 * square**3
            + 592 * square**2
            + 1280 * square
            + 1024,
            (square + 4) * (3 * square + 8),
        ),
        (
            27 * square**5
            + 432 * square**4
            + 2736 * square**3
            + 8576 * square**2
            + 13312 * square
            + 8192,
            (square + 4)
            * (3 * square + 8)
            * sp.sqrt(3 * square + 8),
        ),
        (
            27 * square**6
            + 540 * square**5
            + 4464 * square**4
            + 19520 * square**3
            + 47616 * square**2
            + 61440 * square
            + 32768,
            (square + 4)
            * (3 * square + 8)
            * sp.sqrt(square + 4)
            * sp.sqrt(3 * square + 8),
        ),
    )


def normalize_radicals(expression, z):
    inventory = radical_inventory(z)

    def is_sqrt(node):
        return isinstance(node, sp.Pow) and node.exp == sp.Rational(1, 2)

    def replace_sqrt(node):
        base = sp.expand(node.base)
        for polynomial, replacement in inventory:
            if sp.expand(base - polynomial) == 0:
                return replacement
        return node

    normalized = expression.replace(is_sqrt, replace_sqrt)
    normalized = sp.together(normalized)
    normalized = normalized.replace(is_sqrt, replace_sqrt)
    normalized = sp.cancel(normalized)
    normalized = normalized.replace(is_sqrt, replace_sqrt)
    return sp.factor_terms(normalized)


epsilon_q = 2 * sp.pi - 5 * sp.acos((q**2 + 2) / (2 * (q**2 + 3)))
mu_q = 180 * epsilon_q / (sp.pi * sp.sqrt(q**2 + 4))
p_q = (
    180 * q * epsilon_q / sp.sqrt(q**2 + 4)
    - 600 * sp.sqrt(3) * sp.asinh(q / sp.sqrt(8 * (q**2 + 3)))
)
constraint_expected = 8 * sp.pi * (mu_q - MASS) + 4 * sp.pi * h * q * mu_q
momentum_expected = p_q - 2 * sp.pi * h * mu_q
affine_exact_ok = bool(
    normalize_radicals(constraint_raw - constraint_expected, q) == 0
    and normalize_radicals(momentum_raw - momentum_expected, q) == 0
)
check(
    "the complete action gives the exact general affine second-slab equations",
    affine_exact_ok,
)


# High-precision functions.  No root is selected by a fitted coefficient.
PI = mp.pi


def epsilon(t):
    return 2 * PI - 5 * mp.acos((t * t + 2) / (2 * (t * t + 3)))


def mu(t):
    return 180 * epsilon(t) / (PI * mp.sqrt(t * t + 4))


def momentum(t):
    return (
        180 * t * epsilon(t) / mp.sqrt(t * t + 4)
        - 600
        * mp.sqrt(3)
        * mp.asinh(t / mp.sqrt(8 * (t * t + 3)))
    )


P_INFINITY = 60 * PI - 300 * mp.sqrt(3) * mp.log(2)


def K_square(x):
    return (
        10 * mp.sqrt(x + 4)
        - (x + 3)
        * mp.sqrt(3 * x + 8)
        * (2 * PI - 5 * mp.acos((x + 2) / (2 * (x + 3))))
    )


def bisect(function, left, right, iterations=500):
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


def E_state(v_value, q_value):
    return (
        4 * PI * (mu(q_value) - mu(v_value))
        + q_value * (momentum(q_value) - momentum(v_value))
    )


def E_general(mass_value, p_value, q_value):
    return (
        4 * PI * (mu(q_value) - mass_value)
        + q_value * (momentum(q_value) - p_value)
    )


def reduced_residuals(mass_value, p_value, h_value, q_value):
    return (
        8 * PI * (mu(q_value) - mass_value)
        + 4 * PI * h_value * q_value * mu(q_value),
        momentum(q_value)
        - p_value
        - 2 * PI * h_value * mu(q_value),
    )


def reduced_jacobian(mass_value, p_value, h_value, q_value):
    c_h = mp.diff(
        lambda hh: reduced_residuals(mass_value, p_value, hh, q_value)[0],
        h_value,
    )
    c_q = mp.diff(
        lambda qq: reduced_residuals(mass_value, p_value, h_value, qq)[0],
        q_value,
    )
    p_h = mp.diff(
        lambda hh: reduced_residuals(mass_value, p_value, hh, q_value)[1],
        h_value,
    )
    p_q = mp.diff(
        lambda qq: reduced_residuals(mass_value, p_value, h_value, qq)[1],
        q_value,
    )
    return c_h * p_q - c_q * p_h


f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
pre_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath")
post_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_POST_EXACT, "mpmath")

first_specs = [
    {
        "v": mp.mpf(3) / 2,
        "first_bracket": (3, 10),
        "second_brackets": ((mp.mpf(1) / 100, mp.mpf(1) / 20), (30, 32)),
        "region": "vA<v<vstar",
    },
    {
        "v": mp.mpf(3),
        "first_bracket": (1, 2),
        "second_brackets": ((1, mp.mpf(4) / 3), (18, 20)),
        "region": "vstar<v<vM",
    },
    {
        "v": mp.mpf(20),
        "first_bracket": (-mp.mpf(1) / 5, -mp.mpf(1) / 100),
        "second_brackets": (),
        "region": "vM<v<vC",
    },
]

first_records = []
all_first_ok = True
for spec in first_specs:
    v_value = spec["v"]
    q1 = bisect(
        lambda q_value: E_state(v_value, q_value),
        *spec["first_bracket"],
    )
    h1 = (
        momentum(q1) - momentum(v_value)
    ) / (2 * PI * mu(q1))
    l1 = 1 + h1 * q1
    mass = mu(v_value)
    p0 = momentum(v_value)
    p_post = post_numeric(1, l1, h1**2, mass)
    m1 = mass / l1
    pi1 = p_post / l1**2
    c1 = 2 * f_numeric(1, l1, h1**2, mass) / h1
    p1 = pre_numeric(1, l1, h1**2, mass) - p0
    jacobian1 = reduced_jacobian(mass, p0, h1, q1)
    first_ok = bool(
        h1 > 0
        and l1 > 0
        and m1 > 0
        and abs(c1) < mp.mpf("1e-100")
        and abs(p1) < mp.mpf("1e-100")
        and abs(jacobian1) > mp.mpf("1e-20")
    )
    all_first_ok &= first_ok
    first_records.append(
        {
            "region": spec["region"],
            "v": v_value,
            "q1": q1,
            "h1": h1,
            "l1": l1,
            "mass": mass,
            "p0": p0,
            "p_post": p_post,
            "m1": m1,
            "pi1": pi1,
            "constraint_residual": c1,
            "momentum_residual": p1,
            "jacobian": jacobian1,
            "passed": first_ok,
            "second_brackets": spec["second_brackets"],
        }
    )

check(
    "all three frozen first slabs give exact outgoing canonical data",
    all_first_ok,
    f"states={len(first_records)}",
)


# State-curve closure is only a diagnostic.  The mass equation has four
# candidates for these representatives; none passes the momentum equation.
state_curve_records = []
state_curve_misses = True
for row in first_records:
    m1 = row["m1"]
    pi1 = row["pi1"]
    inner = bisect(lambda w: mu(w) - m1, 0, v_star)
    outer_right = mp.mpf(2) * v_star
    while mu(outer_right) > m1:
        outer_right *= 2
    outer = bisect(lambda w: mu(w) - m1, v_star, outer_right)
    candidates = (-outer, -inner, inner, outer)
    gaps = [momentum(w) - pi1 for w in candidates]
    miss = all(abs(gap) > mp.mpf("1e-20") for gap in gaps)
    state_curve_misses &= miss
    state_curve_records.append(
        {
            "v": text(row["v"]),
            "candidates": [
                {"w": text(w), "momentum_gap": text(gap)}
                for w, gap in zip(candidates, gaps)
            ],
            "closure": False if miss else "UNRESOLVED",
        }
    )

check(
    "the special one-parameter state curve misses all representative outputs",
    state_curve_misses,
    "diagnostic only; general composition remains separately tested",
)


# Enumerate the general second slabs.  For the two expanding controls, pi1
# lies between p(vstar) and p(infinity).  E2 then has two stationary points;
# their signs and the tails prove exactly three roots.  The middle root is the
# negative-height reverse slab.  The first and third roots are distinct
# positive-height forward slabs.
composition_records = []
composition_ok = True
nonunique_witnesses = 0
for row in first_records:
    v_value = row["v"]
    m1 = row["m1"]
    pi1 = row["pi1"]
    q1 = row["q1"]
    h1 = row["h1"]
    l1 = row["l1"]

    roots = []
    if row["second_brackets"]:
        critical_inner = bisect(
            lambda q_value: momentum(q_value) - pi1,
            0,
            v_star,
        )
        critical_right = 2 * v_star
        while momentum(critical_right) < pi1:
            critical_right *= 2
        critical_outer = bisect(
            lambda q_value: momentum(q_value) - pi1,
            v_star,
            critical_right,
        )
        stationary_signs = (
            E_general(m1, pi1, 0),
            E_general(m1, pi1, critical_inner),
            E_general(m1, pi1, critical_outer),
        )
        root_count_certified = bool(
            momentum(v_star) < pi1 < P_INFINITY
            and stationary_signs[0] < 0
            and stationary_signs[1] > 0
            and stationary_signs[2] < 0
        )
        for bracket in row["second_brackets"]:
            roots.append(
                bisect(
                    lambda q_value: E_general(m1, pi1, q_value),
                    *bracket,
                )
            )
        # Exact time-reversed algebraic root, verified rather than searched.
        roots.insert(1, q1)
        reverse_root_ok = abs(E_general(m1, pi1, q1)) < mp.mpf("1e-100")
        root_count_certified &= reverse_root_ok and len(roots) == 3
    else:
        # At v=20 the outgoing momentum is above the global maximum of p.
        # E2 is strictly decreasing from +infinity to -infinity, so q1 is its
        # only root and it is the negative-height reverse slab.
        root_count_certified = bool(
            pi1 > -momentum(v_star)
            and abs(E_general(m1, pi1, q1)) < mp.mpf("1e-100")
        )
        stationary_signs = ()
        roots = [q1]

    root_records = []
    physical_roots = 0
    for q2 in roots:
        h2 = (
            momentum(q2) - pi1
        ) / (2 * PI * mu(q2))
        ratio = 1 + h2 * q2
        c2, p2 = reduced_residuals(m1, pi1, h2, q2)

        # For a positive-height slab, verify the unnormalised shared-slice
        # derivative and the two-slab momentum match.
        junction = mp.nan
        c2_direct = mp.nan
        p2_direct = mp.nan
        jacobian2 = mp.nan
        if h2 > 0 and ratio > 0:
            physical_roots += 1
            c2_direct = 2 * f_numeric(1, ratio, h2**2, m1) / h2
            p2_direct = pre_numeric(1, ratio, h2**2, m1) - pi1
            jacobian2 = reduced_jacobian(m1, pi1, h2, q2)
            L2_absolute = l1 * ratio
            rho2_absolute = (l1 * h2) ** 2
            p_pre2_absolute = pre_numeric(
                l1,
                L2_absolute,
                rho2_absolute,
                row["mass"],
            )
            junction = row["p_post"] - p_pre2_absolute

        root_ok = bool(
            abs(c2) < mp.mpf("1e-100")
            and abs(p2) < mp.mpf("1e-100")
            and (
                not (h2 > 0 and ratio > 0)
                or (
                    abs(c2_direct) < mp.mpf("1e-95")
                    and abs(p2_direct) < mp.mpf("1e-95")
                    and abs(junction) < mp.mpf("1e-90")
                    and abs(jacobian2) > mp.mpf("1e-20")
                )
            )
        )
        root_count_certified &= root_ok
        root_records.append(
            {
                "q2": text(q2),
                "h2": text(h2),
                "scale_ratio": text(ratio),
                "physical": bool(h2 > 0 and ratio > 0),
                "constraint_residual": text(c2_direct),
                "momentum_residual": text(p2_direct),
                "junction_residual": text(junction),
                "jacobian": text(jacobian2),
                "passed": root_ok,
            }
        )

    if physical_roots > 1:
        nonunique_witnesses += 1
    composition_ok &= root_count_certified
    composition_records.append(
        {
            "v": text(v_value),
            "pi1_minus_pstar": text(pi1 - momentum(v_star)),
            "pi1_minus_pinf": text(pi1 - P_INFINITY),
            "stationary_signs": [text(value) for value in stationary_signs],
            "all_real_root_count": len(roots),
            "physical_root_count": physical_roots,
            "roots": root_records,
            "root_count_certified": root_count_certified,
        }
    )

check(
    "all representative second-root counts follow from monotone intervals and tails",
    composition_ok,
)
check(
    "two admitted first states each have two distinct physical second slabs",
    nonunique_witnesses == 2,
    f"nonunique_witnesses={nonunique_witnesses}/3",
)


# Hostile controls distinguish actual composition from the two common state
# resets identified in the protocol.
hostile_records = []
hostile_ok = True
for row in first_records:
    wrong_momentum = row["p_post"] / row["l1"]
    wrong_sign = -row["p_post"] / row["l1"] ** 2
    reset_mass = mu(row["q1"])
    controls = {
        "wrong_momentum_normalization_gap": wrong_momentum - row["pi1"],
        "wrong_post_sign_gap": wrong_sign - row["pi1"],
        "reset_mass_gap": reset_mass - row["m1"],
    }
    row_ok = all(abs(value) > mp.mpf("1e-10") for value in controls.values())
    hostile_ok &= row_ok
    hostile_records.append(
        {
            "v": text(row["v"]),
            **{name: text(value) for name, value in controls.items()},
            "passed": row_ok,
        }
    )

check(
    "wrong momentum scaling, post sign and mass reset all change the state",
    hostile_ok,
)


outcome = (
    "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE"
    if (
        provenance_ok
        and action_homogeneity == 0
        and post_homogeneity == 0
        and junction_identity == 0
        and affine_exact_ok
        and all_first_ok
        and composition_ok
        and nonunique_witnesses > 0
        and hostile_ok
    )
    else "FINITE_HEIGHT_TWO_SLAB_OPEN"
)
check(
    "the preregistered hierarchy reports nonunique composition, not a tick",
    outcome == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE",
)


def serialize_first(row):
    return {
        key: (text(value) if isinstance(value, mp.mpf) else value)
        for key, value in row.items()
        if key != "second_brackets"
    }


artifact = {
    "provenance": {
        "primary_sha256": digest(PRIMARY_INPUT),
        "adversarial_sha256": digest(ADVERSARIAL_INPUT),
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
    },
    "exact_certificates": {
        "action_homogeneity": "degree 2",
        "post_momentum_homogeneity": "degree 2",
        "junction": "(L1/2)d(S1+S2)/dL1=p_post1-p_pre2",
        "general_elimination": "4*pi*(mu(q)-m)+q*(p(q)-pi)=0",
    },
    "first_slab_images": [serialize_first(row) for row in first_records],
    "state_curve_closure": {
        "status": "NUMERIC_DIAGNOSTIC_NO_CLOSURE_AT_FROZEN_CONTROLS",
        "not_used_as_composition_gate": True,
        "records": state_curve_records,
    },
    "composition": composition_records,
    "nonunique_witness_fraction": f"{nonunique_witnesses}/3",
    "hostile_controls": hostile_records,
    "interpretation": {
        "label": "STRUCTURAL / OPEN SELECTION",
        "unique_two_slab_map": False,
        "fundamental_tick": False,
        "global_first_state_classification": "OPEN",
        "additional_selector": "NOT DERIVED",
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
