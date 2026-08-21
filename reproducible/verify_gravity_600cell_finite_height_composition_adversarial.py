#!/usr/bin/env python3
"""Direct two-equation replication of finite-height composition nonuniqueness."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
PRIMARY_INPUT = HERE / "gravity_600cell_finite_height_composition.json"
OUTPUT = HERE / "gravity_600cell_finite_height_composition_adversarial.json"

PRIMARY_SHA256 = (
    "d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447"
)
ADVERSARIAL_PROTOCOL_COMMIT = "e8628bc"
PRIMARY_ARTIFACT_COMMIT = "0f832c6"
PRIMARY_REGISTERED_COMMIT = "c7bffc3"
COMPARISON_CORRECTION_COMMIT = "25c4fed"

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


def text(value, digits=70):
    return mp.nstr(value, digits)


provenance_ok = digest(PRIMARY_INPUT) == PRIMARY_SHA256
check(
    "the frozen primary artifact is present but its roots are not yet read",
    provenance_ok,
    f"adversarial_protocol={ADVERSARIAL_PROTOCOL_COMMIT}",
)


# Reconstruct and redifferentiate the full action without importing any
# scalar elimination coefficient from the primary implementation.
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
scaling = {
    L_MINUS: scale * L_MINUS,
    L_PLUS: scale * L_PLUS,
    RHO: scale**2 * RHO,
    MASS: scale * MASS,
}
redifferentiated_ok = bool(
    sp.simplify(
        sp.powsimp(ACTION.subs(scaling) - scale**2 * ACTION, force=True)
    )
    == 0
    and sp.simplify(
        sp.powsimp(
            P_POST_EXACT.subs(scaling) - scale**2 * P_POST_EXACT,
            force=True,
        )
    )
    == 0
)
check(
    "the independently redifferentiated action has the required scaling",
    redifferentiated_ok,
)


f_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), F_EXACT, "mpmath")
pre_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_PRE_EXACT, "mpmath")
post_numeric = sp.lambdify((L_MINUS, L_PLUS, RHO, MASS), P_POST_EXACT, "mpmath")


def state_functions():
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

    return mu, momentum


def direct_residuals(h_value, q_value, mass_value, p_value):
    endpoint = 1 + h_value * q_value
    return (
        2 * f_numeric(1, endpoint, h_value**2, mass_value) / h_value,
        pre_numeric(1, endpoint, h_value**2, mass_value) - p_value,
    )


def direct_jacobian(h_value, q_value, mass_value, p_value):
    c_h = mp.diff(
        lambda hh: direct_residuals(hh, q_value, mass_value, p_value)[0],
        h_value,
    )
    c_q = mp.diff(
        lambda qq: direct_residuals(h_value, qq, mass_value, p_value)[0],
        q_value,
    )
    p_h = mp.diff(
        lambda hh: direct_residuals(hh, q_value, mass_value, p_value)[1],
        h_value,
    )
    p_q = mp.diff(
        lambda qq: direct_residuals(h_value, qq, mass_value, p_value)[1],
        q_value,
    )
    return c_h * p_q - c_q * p_h


def solve_direct(seed, mass_value, p_value):
    return mp.findroot(
        lambda hh, qq: direct_residuals(hh, qq, mass_value, p_value),
        seed,
        tol=mp.mpf("1e-115"),
        maxsteps=100,
    )


TARGET_PRECISIONS = (80, 120, 180)
FIRST_SEED = (mp.mpf(1) / 5, mp.mpf(10))
BRANCH_SEEDS = {
    "A": (mp.mpf(7), mp.mpf(1) / 50),
    "B": (mp.mpf(1) / 14, mp.mpf(31)),
}
precision_records = []
precision_ok = True
v_value_text = "1.5"

for target_precision in TARGET_PRECISIONS:
    # Guard digits allow every target-precision reconstruction, including the
    # 80-digit one, to be checked against the frozen 1e-90 residual threshold.
    mp.mp.dps = target_precision + 60
    mu, momentum = state_functions()
    v_value = mp.mpf(3) / 2
    mass = mu(v_value)
    p0 = momentum(v_value)

    h1, q1 = solve_direct(FIRST_SEED, mass, p0)
    L1 = 1 + h1 * q1
    p_post = post_numeric(1, L1, h1**2, mass)
    m1 = mass / L1
    pi1 = p_post / L1**2
    first_residuals = direct_residuals(h1, q1, mass, p0)
    first_jacobian = direct_jacobian(h1, q1, mass, p0)

    branches = {}
    row_ok = bool(
        h1 > 0
        and L1 > 0
        and max(abs(value) for value in first_residuals) < mp.mpf("1e-90")
        and abs(first_jacobian) > mp.mpf("1e-20")
    )
    for name, seed in BRANCH_SEEDS.items():
        h2, q2 = solve_direct(seed, m1, pi1)
        ratio = 1 + h2 * q2
        residuals = direct_residuals(h2, q2, m1, pi1)
        jacobian = direct_jacobian(h2, q2, m1, pi1)

        L2 = L1 * ratio
        rho2 = (L1 * h2) ** 2
        pre2_absolute = pre_numeric(L1, L2, rho2, mass)
        junction = p_post - pre2_absolute
        branch_ok = bool(
            h2 > 0
            and ratio > 0
            and max(abs(value) for value in residuals) < mp.mpf("1e-90")
            and abs(jacobian) > mp.mpf("1e-20")
            and abs(junction) < mp.mpf("1e-90")
        )
        row_ok &= branch_ok
        branches[name] = {
            "h2": h2,
            "q2": q2,
            "ratio": ratio,
            "constraint_residual": residuals[0],
            "momentum_residual": residuals[1],
            "jacobian": jacobian,
            "junction_residual": junction,
            "passed": branch_ok,
        }

    distinct = bool(
        abs(branches["A"]["q2"] - branches["B"]["q2"]) > 1
        and abs(branches["A"]["ratio"] - branches["B"]["ratio"]) > 1
    )
    row_ok &= distinct
    precision_ok &= row_ok
    precision_records.append(
        {
            "target_precision": target_precision,
            "h1": h1,
            "q1": q1,
            "L1": L1,
            "p_post": p_post,
            "m1": m1,
            "pi1": pi1,
            "first_jacobian": first_jacobian,
            "branches": branches,
            "distinct": distinct,
            "passed": row_ok,
        }
    )

check(
    "both direct full-action second slabs pass at 80, 120 and 180 digits",
    precision_ok,
)


precision_nesting_ok = True
for lower, upper in zip(precision_records, precision_records[1:]):
    for name in ("h1", "q1", "L1", "p_post", "m1", "pi1"):
        precision_nesting_ok &= bool(
            abs(lower[name] - upper[name]) < mp.mpf("1e-70")
        )
    for branch in ("A", "B"):
        for name in ("h2", "q2", "ratio"):
            precision_nesting_ok &= bool(
                abs(
                    lower["branches"][branch][name]
                    - upper["branches"][branch][name]
                )
                < mp.mpf("1e-70")
            )
check(
    "all independently solved roots are nested beyond 60 decimal digits",
    precision_nesting_ok,
)


# Perturb each rational seed in both directions.  This attacks accidental
# convergence to a seed-specific numerical point without using E2.
mp.mp.dps = 180
mu, momentum = state_functions()
v_value = mp.mpf(3) / 2
mass = mu(v_value)
p0 = momentum(v_value)
h1, q1 = solve_direct(FIRST_SEED, mass, p0)
L1 = 1 + h1 * q1
p_post = post_numeric(1, L1, h1**2, mass)
m1 = mass / L1
pi1 = p_post / L1**2
reference = precision_records[-1]["branches"]

seed_records = []
seed_stability_ok = True
for name, seed in BRANCH_SEEDS.items():
    for h_factor, q_factor in (
        (mp.mpf("0.95"), mp.mpf("0.95")),
        (mp.mpf("0.95"), mp.mpf("1.05")),
        (mp.mpf("1.05"), mp.mpf("0.95")),
        (mp.mpf("1.05"), mp.mpf("1.05")),
    ):
        perturbed = (seed[0] * h_factor, seed[1] * q_factor)
        h2, q2 = solve_direct(perturbed, m1, pi1)
        errors = (
            abs(h2 - reference[name]["h2"]),
            abs(q2 - reference[name]["q2"]),
        )
        row_ok = max(errors) < mp.mpf("1e-70")
        seed_stability_ok &= row_ok
        seed_records.append(
            {
                "branch": name,
                "h_factor": text(h_factor),
                "q_factor": text(q_factor),
                "h_error": text(errors[0]),
                "q_error": text(errors[1]),
                "passed": row_ok,
            }
        )

check(
    "both roots survive every preregistered +-5 percent seed perturbation",
    seed_stability_ok,
    f"perturbations={len(seed_records)}",
)


# The signed affine reverse control is deliberately not evaluated through
# sqrt(h^2) in the full action; it is not a positive-height geometry.
def reduced_signed(h_value, q_value):
    pi = mp.pi
    return (
        8 * pi * (mu(q_value) - m1)
        + 4 * pi * h_value * q_value * mu(q_value),
        momentum(q_value) - pi1 - 2 * pi * h_value * mu(q_value),
    )


h_reverse, q_reverse = mp.findroot(
    reduced_signed,
    (-mp.mpf(1) / 14, mp.mpf(10)),
    tol=mp.mpf("1e-130"),
    maxsteps=100,
)
reverse_ok = bool(
    h_reverse < 0
    and max(abs(value) for value in reduced_signed(h_reverse, q_reverse))
    < mp.mpf("1e-110")
)
check(
    "the reverse algebraic solution has negative height and is rejected",
    reverse_ok,
    f"h_reverse={text(h_reverse, 20)}",
)


# Convention attacks.  They must change the incoming second state or the
# shared derivative; no alternative is silently substituted into the solve.
wrong_scaled_pi = p_post / L1
wrong_sign_pi = -p_post / L1**2
reset_mass = mu(q1)
wrong_sign_junction = p_post - (-p_post)
hostile = {
    "wrong_scale_gap": wrong_scaled_pi - pi1,
    "wrong_sign_gap": wrong_sign_pi - pi1,
    "reset_mass_gap": reset_mass - m1,
    "wrong_sign_junction": wrong_sign_junction,
}
hostile_ok = all(abs(value) > mp.mpf("1e-10") for value in hostile.values())
check(
    "wrong scaling, sign and mass-reset conventions all fail distinctly",
    hostile_ok,
)


# Only now load and compare the primary roots.
primary = json.loads(PRIMARY_INPUT.read_text())
primary_v = next(row for row in primary["composition"] if row["v"] == v_value_text)
primary_physical = [root for root in primary_v["roots"] if root["physical"]]
primary_physical.sort(key=lambda root: mp.mpf(root["q2"]))
comparison = {}
comparison_ok = bool(
    primary["outcome"] == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE"
    and len(primary_physical) == 2
)
for branch, primary_root in zip(("A", "B"), primary_physical):
    differences = {
        "h2": abs(reference[branch]["h2"] - mp.mpf(primary_root["h2"])),
        "q2": abs(reference[branch]["q2"] - mp.mpf(primary_root["q2"])),
        "ratio": abs(
            reference[branch]["ratio"] - mp.mpf(primary_root["scale_ratio"])
        ),
    }
    serialized_matches = {
        "h2": text(reference[branch]["h2"], 60) == primary_root["h2"],
        "q2": text(reference[branch]["q2"], 60) == primary_root["q2"],
        "ratio": (
            text(reference[branch]["ratio"], 60)
            == primary_root["scale_ratio"]
        ),
    }
    row_ok = all(serialized_matches.values())
    comparison_ok &= row_ok
    comparison[branch] = {
        **{name: text(value, 20) for name, value in differences.items()},
        "serialized_60_digit_matches": serialized_matches,
        "passed": row_ok,
    }
check(
    "the direct roots agree with the primary elimination only after construction",
    comparison_ok,
)


complete = bool(
    provenance_ok
    and redifferentiated_ok
    and precision_ok
    and precision_nesting_ok
    and seed_stability_ok
    and reverse_ok
    and hostile_ok
    and comparison_ok
)
outcome = (
    "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED"
    if complete
    else "FINITE_HEIGHT_TWO_SLAB_OPEN"
)
check(
    "the adversarial hierarchy rejects a unique tick map",
    outcome == "FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED",
)


def serialize_precision(row):
    return {
        "target_precision": row["target_precision"],
        "h1": text(row["h1"]),
        "q1": text(row["q1"]),
        "L1": text(row["L1"]),
        "p_post": text(row["p_post"]),
        "m1": text(row["m1"]),
        "pi1": text(row["pi1"]),
        "first_jacobian": text(row["first_jacobian"]),
        "distinct": row["distinct"],
        "passed": row["passed"],
        "branches": {
            name: {
                key: (text(value) if isinstance(value, mp.mpf) else value)
                for key, value in branch.items()
            }
            for name, branch in row["branches"].items()
        },
    }


artifact = {
    "provenance": {
        "primary_sha256": digest(PRIMARY_INPUT),
        "adversarial_protocol_commit": ADVERSARIAL_PROTOCOL_COMMIT,
        "primary_artifact_commit": PRIMARY_ARTIFACT_COMMIT,
        "primary_registered_commit": PRIMARY_REGISTERED_COMMIT,
        "comparison_correction_commit": COMPARISON_CORRECTION_COMMIT,
    },
    "method": {
        "full_action_redifferentiated": True,
        "decisive_solver": "direct two-equation solve in (h,q)",
        "primary_scalar_elimination_used": False,
        "primary_roots_read_only_after_direct_construction": True,
        "state_curve_closure_used": False,
    },
    "precision_runs": [serialize_precision(row) for row in precision_records],
    "seed_perturbations": seed_records,
    "reverse_control": {
        "h": text(h_reverse),
        "q": text(q_reverse),
        "rejected": True,
    },
    "hostile_controls": {name: text(value) for name, value in hostile.items()},
    "primary_comparison": comparison,
    "interpretation": {
        "label": "DERIVED NEGATIVE FOR UNIQUENESS / STRUCTURAL OPEN SELECTION",
        "unique_deterministic_composition": False,
        "additional_selector": "NOT DERIVED",
        "fundamental_tick": False,
        "global_first_state_classification": "OPEN",
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
