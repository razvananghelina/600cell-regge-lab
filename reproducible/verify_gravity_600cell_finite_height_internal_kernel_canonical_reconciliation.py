#!/usr/bin/env python3
"""Exact canonical reconciliation of the finite-height internal kernel."""

import ast
from collections import Counter
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


mp.mp.dps = 180

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = (
    HERE
    / "gravity_600cell_finite_height_internal_kernel_canonical_reconciliation.json"
)
PRIOR_ART = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_internal_kernel_canonical_reconciliation_prior_art.md"
)
PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_finite_height_internal_kernel_canonical_reconciliation_protocol.md"
)
PRIMARY_JSON = HERE / "gravity_600cell_finite_height_internal_carrier_rank.json"
PRIMARY_MATRICES = (
    HERE / "gravity_600cell_finite_height_internal_carrier_rank_matrices.npz"
)
ADVERSARIAL_JSON = (
    HERE / "gravity_600cell_finite_height_internal_carrier_rank_adversarial.json"
)
ADVERSARIAL_MATRICES = (
    HERE
    / "gravity_600cell_finite_height_internal_carrier_rank_adversarial_matrices.npz"
)
SELECTOR_JSON = HERE / "gravity_600cell_finite_height_selector_audit.json"
SELECTOR_SOURCE = HERE / "verify_gravity_600cell_finite_height_selector_audit.py"
CARRIER_PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_finite_height_carrier_quadratic_protocol.md"
)
BACKGROUND_JSON = HERE / "gravity_600cell_finite_height_carrier_quadratic.json"
RUN_ALL = HERE / "run_all.py"

INPUTS = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "primary_json": PRIMARY_JSON,
    "primary_matrices": PRIMARY_MATRICES,
    "adversarial_json": ADVERSARIAL_JSON,
    "adversarial_matrices": ADVERSARIAL_MATRICES,
    "selector_json": SELECTOR_JSON,
    "selector_source": SELECTOR_SOURCE,
    "carrier_protocol": CARRIER_PROTOCOL,
    "background_json": BACKGROUND_JSON,
}
EXPECTED_HASHES = {
    "prior_art": (
        "9e474113a426c61e44fcce000ff6e4a0262bc00c951a09dd72eb3cca347d0c4b"
    ),
    "protocol": (
        "3c0437313cc52ab786318c2ce2c2ded2987aa2bc71a49e9a57a509dc1372c37b"
    ),
    "primary_json": (
        "513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8"
    ),
    "primary_matrices": (
        "97f5b8318be2b3ccf843db87e678ac1ac6ce402db262023c6bbc63a7b647321b"
    ),
    "adversarial_json": (
        "ddd4704b7d1deb6360e752b2ebfe5cc0b66d03819c9f8df7b74e24373aa98fb5"
    ),
    "adversarial_matrices": (
        "45ee09642485dc0e18c6378b9454882414bb07d7e5ad9dbd0c3a8896fd8a7f74"
    ),
    "selector_json": (
        "956cd655b8b3a5106029fb852df74b85bb59f922a4984542bc2e089f54799676"
    ),
    "selector_source": (
        "aca44fa0cf0f6a464ca1a8eaa61356941e8408950ae3801085b85ce314741503"
    ),
    "carrier_protocol": (
        "f73ee892258e33d43991fc8c74bc6f44e6c7f2ae57be56f057050c86ff646fad"
    ),
    "background_json": (
        "0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9"
    ),
}

PRIOR_ART_COMMIT = "5c30454"
PROTOCOL_COMMIT = "07f2b96"
REGISTRY_COMMIT = "e009970"
VERIFIER_NAME = Path(__file__).name
VERTICES = 120
DATA = 240

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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=80):
    return mp.nstr(value, digits)


def float_text(value):
    return f"{float(value):.17e}"


def registry_inventory(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    scripts = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "scripts"
            for target in node.targets
        ):
            scripts = ast.literal_eval(node.value)
            break
    if scripts is None:
        raise RuntimeError("run_all.py has no literal scripts registry")
    counts = Counter(scripts)
    duplicates = sorted(name for name, count in counts.items() if count != 1)
    return scripts, duplicates


def analytic_projector(ratio):
    vector = np.r_[
        np.ones(VERTICES), np.full(VERTICES, float(ratio))
    ]
    vector /= np.linalg.norm(vector)
    return np.outer(vector, vector), vector


def projector_classification(left, right, uncertainty):
    difference = float(np.linalg.norm(left - right, ord=2))
    if difference <= 10 * uncertainty:
        classification = "AGREE"
    elif difference > 100 * uncertainty:
        classification = "DEPENDENT"
    else:
        classification = "OPEN"
    return {
        "classification": classification,
        "difference_two_norm": float_text(difference),
        "uncertainty": float_text(uncertainty),
        "agreement_gate": float_text(10 * uncertainty),
        "dependence_gate": float_text(100 * uncertainty),
    }


def numerical_bridge(q_text, h_text, digits):
    with mp.workdps(digits):
        q_value = mp.mpf(q_text)
        h_value = mp.mpf(h_text)

        def epsilon(value):
            return 2 * mp.pi - 5 * mp.acos(
                (value**2 + 2) / (2 * (value**2 + 3))
            )

        def mu(value):
            return 180 * epsilon(value) / (
                mp.pi * mp.sqrt(value**2 + 4)
            )

        def momentum(value):
            return (
                180 * value * epsilon(value) / mp.sqrt(value**2 + 4)
                - 600 * mp.sqrt(3) * mp.asinh(
                    value / mp.sqrt(8 * (value**2 + 3))
                )
            )

        u_value = mu(q_value)
        up_value = mp.diff(mu, q_value)
        pp_value = mp.diff(momentum, q_value)
        c_h = 4 * mp.pi * q_value * u_value
        c_q = (
            8 * mp.pi * up_value
            + 4 * mp.pi * h_value * (u_value + q_value * up_value)
        )
        p_h = -2 * mp.pi * u_value
        p_q = pp_value - 2 * mp.pi * h_value * up_value
        c_sigma = c_q / h_value
        c_c = (h_value * c_h - q_value * c_q) / 2
        p_sigma = p_q / h_value
        p_c = (h_value * p_h - q_value * p_q) / 2
        ratio = -c_sigma / c_c
        d_c = c_sigma + c_c * ratio
        d_p = p_sigma + p_c * ratio
        determinant_hq = c_h * p_q - c_q * p_h
        expected_hq = 8 * mp.pi**2 * h_value * u_value**2
        determinant_sc = c_sigma * p_c - c_c * p_sigma
        expected_sc = -4 * mp.pi**2 * h_value * u_value**2
        return {
            "digits": digits,
            "q": q_value,
            "h": h_value,
            "mu": u_value,
            "mu_prime": up_value,
            "p_prime": pp_value,
            "C_h": c_h,
            "C_q": c_q,
            "P_h": p_h,
            "P_q": p_q,
            "C_sigma": c_sigma,
            "C_c": c_c,
            "P_sigma": p_sigma,
            "P_c": p_c,
            "ratio": ratio,
            "dC": d_c,
            "dP": d_p,
            "determinant_hq": determinant_hq,
            "expected_hq": expected_hq,
            "determinant_sc": determinant_sc,
            "expected_sc": expected_sc,
        }


input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
provenance_ok = input_hashes == EXPECTED_HASHES
check("all reconciliation inputs retain their frozen hashes", provenance_ok)

scripts, duplicates = registry_inventory(RUN_ALL)
registry_ok = bool(scripts.count(VERIFIER_NAME) == 1 and not duplicates)
check(
    "the reconciliation verifier is registered exactly once with no duplicates",
    registry_ok,
    f"entries={len(scripts)}, duplicates={duplicates}",
)

selector = json.loads(SELECTOR_JSON.read_text())
selector_ok = bool(
    selector["outcome"] == "STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH"
    and selector["tests"] == selector["passed"] == 10
    and selector["exact_certificates"]["state_derivative_identity"]
    == "4*pi*mu'(q)+q*p'(q)=0"
    and selector["exact_certificates"]["legendre_determinant"]
    == "8*pi^2*h*mu(q)^2"
)
check("the exact selector determinant certificate remains intact", selector_ok)

background = json.loads(BACKGROUND_JSON.read_text())["background"]
q_value = mp.mpf(background["q"])
h_value = mp.mpf(background["h"])
lambda_value = mp.mpf(background["lambda"])
rho_value = mp.mpf(background["rho"])
background_ok = bool(
    q_value > 0 and h_value > 0 and lambda_value > 0 and rho_value > 0
    and abs(lambda_value - (1 + h_value * q_value)) < mp.mpf("1e-70")
    and abs(rho_value - h_value**2) < mp.mpf("1e-70")
)
check(
    "the frozen positive-height background obeys lambda=1+h*q and rho=h^2",
    background_ok,
)

# Exact symbolic calculation before any stored kernel vector is loaded.
h, q, u, up, pp = sp.symbols("h q u up pp", positive=True, real=True)
pi = sp.pi
c_h = 4 * pi * q * u
c_q = 8 * pi * up + 4 * pi * h * (u + q * up)
p_h = -2 * pi * u
p_q = pp - 2 * pi * h * up
jacobian_hq = sp.Matrix([[c_h, c_q], [p_h, p_q]])
remainder = sp.factor(jacobian_hq.det() - 8 * pi**2 * h * u**2)
expected_remainder = 4 * pi * u * (4 * pi * up + q * pp)
determinant_reduction_ok = sp.expand(remainder - expected_remainder) == 0
check(
    "the exact (h,q) determinant reduces to the state derivative identity",
    determinant_reduction_ok,
    f"remainder={remainder}",
)

k, denominator = sp.symbols("K denominator", real=True, nonzero=True)
mu_prime = 180 * q * k / (pi * denominator)
p_prime = -720 * k / denominator
state_identity = sp.simplify(4 * pi * mu_prime + q * p_prime)
exact_hq = sp.simplify(
    jacobian_hq.det().subs({up: mu_prime, pp: p_prime})
)
state_identity_ok = bool(
    state_identity == 0
    and sp.simplify(exact_hq - 8 * pi**2 * h * u**2) == 0
)
check(
    "the independently reconstructed derivative identity gives a positive determinant",
    state_identity_ok,
)

coordinate = sp.Matrix([[0, h / 2], [1 / h, -q / 2]])
coordinate_determinant_ok = sp.simplify(coordinate.det() + sp.Rational(1, 2)) == 0
jacobian_sc = sp.simplify(jacobian_hq * coordinate)
exact_sc = sp.simplify(
    jacobian_sc.det().subs({up: mu_prime, pp: p_prime})
)
coordinate_identity_ok = bool(
    coordinate_determinant_ok
    and sp.simplify(exact_sc + 4 * pi**2 * h * u**2) == 0
)
check(
    "the exact carrier conversion has determinant -1/2 and preserves regularity",
    coordinate_identity_ok,
    f"det={coordinate.det()}",
)

c_sigma, c_c = jacobian_sc[0, 0], jacobian_sc[0, 1]
p_sigma, p_c = jacobian_sc[1, 0], jacobian_sc[1, 1]
tangent = sp.Matrix([c_c, -c_sigma])
symbolic_dc = sp.simplify((jacobian_sc[0, :] * tangent)[0])
symbolic_dp = sp.simplify(
    (jacobian_sc[1, :] * tangent)[0].subs({up: mu_prime, pp: p_prime})
)
tangent_identity_ok = bool(
    symbolic_dc == 0
    and sp.simplify(symbolic_dp - 4 * pi**2 * h * u**2) == 0
)
check(
    "the exact C tangent has nonzero canonical-momentum response",
    tangent_identity_ok,
    f"dP={symbolic_dp}",
)

wrong_half = sp.Matrix([[0, h], [1 / h, -q]])
wrong_half_detected = bool(
    sp.simplify(wrong_half.det() + sp.Rational(1, 2)) != 0
    and sp.simplify(wrong_half.det() + 1) == 0
)
check(
    "omitting the one-half in delta h=h*c/2 is rejected exactly",
    wrong_half_detected,
)

bridges = {
    digits: numerical_bridge(background["q"], background["h"], digits)
    for digits in (120, 180)
}
bridge_120 = bridges[120]
bridge_180 = bridges[180]
ratio_precision_error = abs(
    bridge_120["ratio"] - bridge_180["ratio"]
) / abs(bridge_180["ratio"])
determinant_error = max(
    abs(bridge_180["determinant_hq"] - bridge_180["expected_hq"])
    / abs(bridge_180["expected_hq"]),
    abs(bridge_180["determinant_sc"] - bridge_180["expected_sc"])
    / abs(bridge_180["expected_sc"]),
)
bridge_ok = bool(
    ratio_precision_error < mp.mpf("1e-100")
    and abs(bridge_180["dC"]) < mp.mpf("1e-100")
    and determinant_error < mp.mpf("1e-100")
    and abs(bridge_180["dP"]) > mp.mpf("1e-3")
    and bridge_180["mu"] > 0
)
check(
    "120/180-digit definitions reproduce a stable C tangent and nonzero dP",
    bridge_ok,
    (
        f"c/sigma={mp_text(bridge_180['ratio'], 40)}, "
        f"dP={mp_text(bridge_180['dP'], 12)}, "
        f"precision={mp_text(ratio_precision_error, 6)}"
    ),
)

# Stored rank artifacts are read only after the exact/numerical tangent is fixed.
primary = json.loads(PRIMARY_JSON.read_text())
adversarial = json.loads(ADVERSARIAL_JSON.read_text())
rank_inputs_ok = bool(
    primary["outcome"] == "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY"
    and primary["tests"] == primary["passed"] == 25
    and primary["global_nullities"] == {"even": 1, "odd": 1}
    and primary["kernel_comparison"]["classification"] == "AGREE"
    and adversarial["outcome"]
    == "FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_ADVERSARIALLY_REPLICATED"
    and adversarial["tests"] == adversarial["passed"] == 19
    and adversarial["matrix_reproduces_primary"]
    and adversarial["direct_secants_validate"]
    and all(
        adversarial["replication_classification"][parity]
        == {
            "diagonal_nullity": 121,
            "full_nullity": 1,
            "resolved": True,
        }
        for parity in ("even", "odd")
    )
)
check("both independently built rank-one kernel artifacts remain intact", rank_inputs_ok)

primary_arrays = np.load(PRIMARY_MATRICES, allow_pickle=False)
adversarial_arrays = np.load(ADVERSARIAL_MATRICES, allow_pickle=False)
analytic, analytic_vector = analytic_projector(bridge_180["ratio"])
uncertainty = max(
    float(primary["kernel_comparison"]["uncertainty"]),
    float(adversarial["parity_kernel_comparison"]["uncertainty"]),
)
projector_records = {}
for method, arrays in (
    ("primary", primary_arrays), ("adversarial", adversarial_arrays)
):
    projector_records[method] = {
        parity: projector_classification(
            analytic, arrays[f"{parity}_kernel_projector"], uncertainty
        )
        for parity in ("even", "odd")
    }
projectors_resolved = bool(
    all(
        record["classification"] in {"AGREE", "DEPENDENT"}
        for method in projector_records.values() for record in method.values()
    )
)
projectors_agree = bool(
    all(
        record["classification"] == "AGREE"
        for method in projector_records.values() for record in method.values()
    )
)
check(
    "the analytic line is resolved against four independently stored projectors",
    projectors_resolved,
    str({
        method: {parity: record["difference_two_norm"]
                 for parity, record in records.items()}
        for method, records in projector_records.items()
    }),
)

wrong_log_lambda, _ = analytic_projector(
    lambda_value * bridge_180["ratio"]
)
wrong_linear_rho, _ = analytic_projector(
    rho_value * bridge_180["ratio"]
)
hostile_projectors = {
    "sigma_is_delta_log_lambda": projector_classification(
        wrong_log_lambda, adversarial_arrays["even_kernel_projector"], uncertainty
    ),
    "c_is_delta_rho": projector_classification(
        wrong_linear_rho, adversarial_arrays["even_kernel_projector"], uncertainty
    ),
}
hostile_ok = all(
    record["classification"] == "DEPENDENT"
    for record in hostile_projectors.values()
)
check(
    "wrong logarithmic/linear carrier coordinate conventions are strongly rejected",
    hostile_ok,
    str(hostile_projectors),
)

controls_ok = bool(
    provenance_ok and registry_ok and selector_ok and background_ok
    and determinant_reduction_ok and state_identity_ok
    and coordinate_identity_ok and tangent_identity_ok
    and wrong_half_detected and bridge_ok and rank_inputs_ok and hostile_ok
)
if not controls_ok:
    outcome = "CANONICAL_RECONCILIATION_CONTROL_FAILED"
elif not projectors_resolved:
    outcome = "CANONICAL_RECONCILIATION_OPEN"
elif not projectors_agree:
    outcome = "INTERNAL_KERNEL_NOT_THE_HOMOGENEOUS_CONSTRAINT_TANGENT"
else:
    outcome = (
        "INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_"
        "FIXED_INPUT_REMOVES_IT"
    )

check(
    "the canonical reconciliation outcome follows the frozen hierarchy",
    outcome in {
        "CANONICAL_RECONCILIATION_CONTROL_FAILED",
        "CANONICAL_RECONCILIATION_OPEN",
        "INTERNAL_KERNEL_NOT_THE_HOMOGENEOUS_CONSTRAINT_TANGENT",
        "INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_FIXED_INPUT_REMOVES_IT",
    },
    outcome,
)

artifact = {
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "registry_commit": REGISTRY_COMMIT,
        "input_sha256": input_hashes,
    },
    "exact_certificates": {
        "state_derivative_identity": "4*pi*mu'(q)+q*p'(q)=0",
        "hq_determinant": "8*pi^2*h*mu(q)^2",
        "carrier_coordinate_matrix": "[[0,h/2],[1/h,-q/2]]",
        "carrier_coordinate_determinant": "-1/2",
        "sigma_c_determinant": "-4*pi^2*h*mu(q)^2",
        "fixed_input_intersection": "ker(R_p) intersect ker(dP) = {0}",
    },
    "background": {
        "q": background["q"],
        "h": background["h"],
        "lambda": background["lambda"],
        "rho": background["rho"],
    },
    "numerical_bridge": {
        "ratio_c_over_sigma": mp_text(bridge_180["ratio"]),
        "dC_on_normalized_sigma_tangent": mp_text(bridge_180["dC"]),
        "dP_on_normalized_sigma_tangent": mp_text(bridge_180["dP"]),
        "mu": mp_text(bridge_180["mu"]),
        "hq_determinant": mp_text(bridge_180["determinant_hq"]),
        "sigma_c_determinant": mp_text(bridge_180["determinant_sc"]),
        "ratio_120_180_relative_error": mp_text(ratio_precision_error),
        "determinant_relative_error": mp_text(determinant_error),
        "analytic_vector_scale_mean": float_text(np.mean(analytic_vector[:VERTICES])),
        "analytic_vector_strut_mean": float_text(np.mean(analytic_vector[VERTICES:])),
    },
    "projector_uncertainty": float_text(uncertainty),
    "projector_comparisons": projector_records,
    "hostile_coordinate_controls": hostile_projectors,
    "interpretation": {
        "label": "DERIVED EXACT/COMPUTATIONAL, BOUNDED NEGATIVE",
        "internal_kernel": "homogeneous lapse-constraint tangent",
        "fixed_incoming_canonical_tangent": "zero",
        "forced_response_to_varying_incoming_data": "NOT_TESTED",
        "physical_tick": "NOT_DERIVED",
        "gravitons_wave_equation_c_G_planck_particle_masses": "NOT_DERIVED",
        "finite_step_invariant_continuation": "SEPARATE_DERIVED_RESULT",
        "infinite_proper_time_evolution": "OPEN",
    },
    "firewall": {
        "exploratory_ratio_disclosed": True,
        "stored_kernel_vectors_loaded_after_exact_ratio": True,
        "full_suite_run": False,
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} PASS")
if passed != tests:
    raise SystemExit(1)
