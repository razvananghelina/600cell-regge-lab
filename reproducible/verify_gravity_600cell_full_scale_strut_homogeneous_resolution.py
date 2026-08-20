#!/usr/bin/env python3
"""Resolve the exact homogeneous weak-pole canonical intersection line."""

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIOR = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_prior_art.md"
DISCLOSURE = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_disclosure.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_protocol.md"
ACTION_SOURCE = HERE / "verify_gravity_600cell_homothetic_frustum_action.py"
ACTION_INPUT = HERE / "gravity_600cell_homothetic_frustum_action.json"
LAPSE_SOURCE = HERE / "verify_gravity_600cell_dust_homothetic_canonical_lapse.py"
LAPSE_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
CARRIER_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_carrier.py"
CARRIER_INPUT = HERE / "gravity_600cell_full_scale_strut_carrier.json"
SYMBOLIC_INPUT = HERE / "gravity_600cell_full_scale_strut_symbolic_gap_resolution.json"
BINARY_INPUT = HERE / "gravity_600cell_full_scale_strut_canonical_intersection.json"
PRECISION_SOURCE = HERE / "verify_gravity_600cell_full_scale_strut_canonical_precision.py"
PRECISION_INPUT = HERE / "gravity_600cell_full_scale_strut_canonical_precision.json"
ADVERSARIAL_INPUT = HERE / "gravity_600cell_full_scale_strut_canonical_precision_adversarial.json"
OUTPUT = HERE / "gravity_600cell_full_scale_strut_homogeneous_resolution.json"

PRIOR_COMMIT = "16e28ec"
PROTOCOL_COMMIT = "26697d7"
EXPECTED_HASHES = {
    "prior": "661455e95e33c5c7f5d18627e8b31892c3ca85944d3904323e2b53de7a94cf81",
    "disclosure": "8f936f6fcd87dff18aff8c409821df9f6a740507ec9dc1fa2e10aaa8ad1720e8",
    "protocol": "94260ae7c1a11d4d58fb5f829cb4f67ff39c8d893df86a991c0f0ef52359c4ad",
    "action_source": "61de49dd88f614044d2c24fcbbe02e6fb3bc39c05bd963c63c5bc419c4bbf0cd",
    "action": "c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d",
    "lapse_source": "8ae83004dcdeadfde27b91947a1c517915fa59af60807c9b128406a20c63508c",
    "lapse": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "carrier_source": "e68105df4058f7d2ed39a6913f29e88cd9fe88e123ff52260acf698a2bd7da49",
    "carrier": "6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d",
    "symbolic": "ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179",
    "binary": "b29cc33a9effeb2087fb6133359ee747d100d203778586372a7ceeebc2e4f070",
    "precision_source": "b836cad394fe8a54644d514b6f31cb899ff5c3697b6c2a5f4edfc2b5f0ac5d62",
    "precision": "75351ae4dfde26dd75ed8faa927b0a49cd725d83c7629d4545268030b54e2706",
    "adversarial": "ecf02fd76b0c1d4d95cd206c639a027400c2053bdb1850018d57ff2721861db3",
}
INPUTS = {
    "prior": PRIOR,
    "disclosure": DISCLOSURE,
    "protocol": PROTOCOL,
    "action_source": ACTION_SOURCE,
    "action": ACTION_INPUT,
    "lapse_source": LAPSE_SOURCE,
    "lapse": LAPSE_INPUT,
    "carrier_source": CARRIER_SOURCE,
    "carrier": CARRIER_INPUT,
    "symbolic": SYMBOLIC_INPUT,
    "binary": BINARY_INPUT,
    "precision_source": PRECISION_SOURCE,
    "precision": PRECISION_INPUT,
    "adversarial": ADVERSARIAL_INPUT,
}

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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(value, digits=70):
    return mp.nstr(value, digits)


def relative_change(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))


def vector_norm(vector):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in vector))


def normalized_projector(vector):
    norm = vector_norm(vector)
    unit = [value / norm for value in vector]
    return [[left * mp.conj(right) for right in unit] for left in unit]


def projector_distance(left, right):
    return mp.sqrt(mp.fsum(
        abs(left[row][column] - right[row][column]) ** 2
        for row in range(len(left)) for column in range(len(left))
    ))


print("=" * 78)
print("EXACT HOMOGENEOUS WEAK-POLE CANONICAL LINE")
print("=" * 78)

hashes = {name: digest(path) for name, path in INPUTS.items()}
action = json.loads(ACTION_INPUT.read_text())
lapse = json.loads(LAPSE_INPUT.read_text())
carrier = json.loads(CARRIER_INPUT.read_text())
symbolic = json.loads(SYMBOLIC_INPUT.read_text())
binary = json.loads(BINARY_INPUT.read_text())
precision = json.loads(PRECISION_INPUT.read_text())
adversarial = json.loads(ADVERSARIAL_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and action["outcome"] == "HOMOTHETIC_FRUSTUM_ACTION_INVARIANT"
    and action["passed"] == action["tests"] == 16
    and lapse["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and lapse["passed"] == lapse["tests"] == 7
    and carrier["passed"] == carrier["tests"] == 18
    and symbolic["outcome"] == "FULL_SCALE_STRUT_GAP_REAL_RESOLVED"
    and symbolic["passed"] == symbolic["tests"] == 11
    and binary["outcome"] == "FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN"
    and binary["passed"] == binary["tests"] == 13
    and precision["outcome"] == "FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN"
    and precision["passed"] == precision["tests"] == 17
    and adversarial["outcome"] == "NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED"
    and adversarial["passed"] == adversarial["tests"] == 7
)
check("all homogeneous-resolution inputs retain frozen provenance", provenance_ok)

# Independently reconstruct the closed action with new symbols.
L_MINUS, L_PLUS, RHO = sp.symbols("L_minus L_plus rho", positive=True)
DELTA = L_PLUS - L_MINUS
H = sp.sqrt(RHO + DELTA**2 / 4)
COSINE = (DELTA**2 + 2 * RHO) / (2 * (DELTA**2 + 3 * RHO))
BOOST = DELTA / sp.sqrt(8 * (DELTA**2 + 3 * RHO))
S = sp.factor(
    360 * (L_MINUS + L_PLUS) * H * (2 * sp.pi - 5 * sp.acos(COSINE))
    + 600 * sp.sqrt(3) * (L_MINUS**2 - L_PLUS**2) * sp.asinh(BOOST)
)
formula_ok = bool(
    str(S) == action["exact_geometry"]["closed_gravitational_action"]
    and action["exact_geometry"]["symbolic_action_identity"]
    and action["exact_geometry"]["flat_subdivision_hypotheses"]
    and action["exact_geometry"]["diagonal_derivative_sum"] == "0"
    and action["scope"]["homogeneous_collective_action_and_momenta"] == "TESTED"
)
check("the independent closed action matches the frozen subdivision identity", formula_ok)

P_MINUS = sp.factor((L_MINUS / 2) * sp.diff(S, L_MINUS))
P_S = sp.diff(P_MINUS, L_PLUS) * L_PLUS
P_Z = sp.diff(P_MINUS, RHO) * RHO
LAMBDA = L_PLUS / L_MINUS
SIGMA = -LAMBDA * P_Z
C = P_S
correct_momentum = sp.factor(P_S * (SIGMA / LAMBDA) + P_Z * C)

Q_DIAG = L_MINUS * L_PLUS - RHO
delta_log_diagonal = sp.factor(
    (sp.diff(Q_DIAG, L_PLUS) * L_PLUS * (SIGMA / LAMBDA)
     + sp.diff(Q_DIAG, RHO) * RHO * C) / Q_DIAG
)
expected_diagonal = sp.factor(
    (L_MINUS**2 * SIGMA - RHO * C) / Q_DIAG
)
delta_log_upper = sp.factor(
    sp.diff(sp.log(L_PLUS**2), L_PLUS) * L_PLUS * (SIGMA / LAMBDA)
)
delta_log_pole = sp.factor(
    sp.diff(sp.log(RHO), RHO) * RHO * C
)
DUST_MASS = sp.symbols("M", positive=True)
DUST = -8 * sp.pi * DUST_MASS * sp.sqrt(RHO)
symbolic_ok = bool(
    correct_momentum == 0
    and sp.factor(delta_log_diagonal - expected_diagonal) == 0
    and sp.factor(delta_log_upper - 2 * SIGMA / LAMBDA) == 0
    and sp.factor(delta_log_pole - C) == 0
    and sp.diff(DUST, L_MINUS) == 0
)
check("the exact momentum and all three carrier differentials vanish identically", symbolic_ok)

wrong_conversion = sp.factor(P_S * SIGMA + P_Z * C)
wrong_expected = sp.factor(P_S * P_Z * (1 - LAMBDA))
wrong_sign = sp.factor(P_S * P_Z + P_Z * P_S)
corruption_symbolic_ok = bool(
    sp.factor(wrong_conversion - wrong_expected) == 0
    and sp.factor(wrong_sign - 2 * P_S * P_Z) == 0
)
check("the missing-lambda and wrong-sign corruptions retain nonzero formulas", corruption_symbolic_ok)

PS_FUNCTION = sp.lambdify((L_MINUS, L_PLUS, RHO), P_S, modules="mpmath")
PZ_FUNCTION = sp.lambdify((L_MINUS, L_PLUS, RHO), P_Z, modules="mpmath")


def evaluate_background(parity, dps):
    mp.mp.dps = dps
    state = lapse["solutions"][parity]["state"]
    s = mp.mpf(state[0])
    z = mp.mpf(state[1])
    m_star = mp.mpf(10)
    zeta = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
    r0 = 4 * m_star / (3 * mp.pi)
    l0 = zeta * r0
    lam = mp.exp(s)
    rho = mp.mpf("0.0102") ** 2 * mp.exp(z)
    ps = PS_FUNCTION(l0, lam * l0, rho)
    pz = PZ_FUNCTION(l0, lam * l0, rho)
    ratio = -lam * pz / ps
    wrong = ps * pz * (1 - lam)
    sign_corrupt = 2 * ps * pz
    return {
        "s": s,
        "z": z,
        "l0": l0,
        "lambda": lam,
        "rho": rho,
        "p_s": ps,
        "p_z": pz,
        "ratio": ratio,
        "wrong_conversion": wrong,
        "wrong_sign": sign_corrupt,
    }


evaluations = {
    parity: {
        level: evaluate_background(parity, dps)
        for level, dps in (("P160", 160), ("P220", 220))
    }
    for parity in ("even", "odd")
}
numeric_convergence_ok = True
for parity in ("even", "odd"):
    low = evaluations[parity]["P160"]
    high = evaluations[parity]["P220"]
    numeric_convergence_ok &= bool(
        abs(low["p_s"]) > mp.mpf("1e5")
        and abs(low["p_z"]) > mp.mpf("1e-1")
        and abs(high["p_s"]) > mp.mpf("1e5")
        and abs(high["p_z"]) > mp.mpf("1e-1")
        and relative_change(low["p_s"], high["p_s"]) < mp.mpf("1e-140")
        and relative_change(low["p_z"], high["p_z"]) < mp.mpf("1e-140")
        and relative_change(low["ratio"], high["ratio"]) < mp.mpf("1e-140")
        and abs(high["wrong_conversion"]) > mp.mpf("1e-3")
        and abs(high["wrong_sign"]) > mp.mpf("1e5")
    )
check("P160/P220 automatic derivatives converge and both corruptions fail", numeric_convergence_ok)

parity_input_ok = bool(
    lapse["parity_gate"]["passed"]
    and mp.mpf(lapse["parity_gate"]["s_difference"]) < mp.mpf("1e-80")
    and mp.mpf(lapse["parity_gate"]["z_difference"]) < mp.mpf("1e-80")
)
generator_projectors = {}
bridge_records = {}
endpoint_bridge_ok = True
candidate_bridge_ok = True

for parity in ("even", "odd"):
    high = evaluations[parity]["P220"]
    endpoint = lapse["solutions"][parity]["endpoint_jacobian"]
    jacobian = endpoint["matrices"]["operational_primary"]
    j_s = mp.mpf(jacobian[1][0])
    j_z = mp.mpf(jacobian[1][1])
    epsilon = mp.mpf(endpoint["epsilon"])
    jacobian_ratio = -high["lambda"] * j_z / j_s
    denominator_margin = abs(j_s) - epsilon
    ratio_error_bound = high["lambda"] * (
        epsilon / denominator_margin
        + abs(j_z) * epsilon / (abs(j_s) * denominator_margin)
    )
    endpoint_error = abs(high["ratio"] - jacobian_ratio)
    endpoint_bridge_ok &= bool(
        endpoint["resolved"]
        and denominator_margin > 0
        and endpoint_error <= ratio_error_bound
    )

    binary_sector = binary["parities"][parity]["sectors"][6]
    scale_norm = mp.mpf(binary_sector["scaling"]["scale"])
    strut_norm = mp.mpf(binary_sector["scaling"]["strut"])
    frozen = precision["no_refit_validation"][parity][
        "candidate_P100_rounded_70_digits"
    ]
    scale_values = [mp.mpf(value["real"]) for value in frozen[:5]]
    strut_values = [mp.mpf(value["real"]) for value in frozen[5:]]
    scale_mean = mp.fsum(scale_values) / 5
    strut_mean = mp.fsum(strut_values) / 5
    physical_candidate_ratio = scale_norm * scale_mean / (strut_norm * strut_mean)
    scale_spread = max(abs(value - scale_mean) for value in scale_values)
    strut_spread = max(abs(value - strut_mean) for value in strut_values)
    candidate_error = abs(physical_candidate_ratio - high["ratio"])
    candidate_bridge_ok &= bool(
        candidate_error < mp.mpf("1e-30")
        and mp.mpf(0) < scale_spread < mp.mpf("1e-30")
        and mp.mpf(0) < strut_spread < mp.mpf("1e-30")
    )

    sigma = -high["lambda"] * high["p_z"]
    c_value = high["p_s"]
    generator = [sigma / scale_norm] * 5 + [c_value / strut_norm] * 5
    generator_projectors[parity] = normalized_projector(generator)
    bridge_records[parity] = {
        "p_s": text(high["p_s"]),
        "p_z": text(high["p_z"]),
        "lambda": text(high["lambda"]),
        "rho": text(high["rho"]),
        "analytic_sigma_over_c": text(high["ratio"]),
        "endpoint_sigma_over_c": text(jacobian_ratio),
        "endpoint_absolute_error": text(endpoint_error),
        "endpoint_error_bound": text(ratio_error_bound),
        "P100_candidate_sigma_over_c": text(physical_candidate_ratio),
        "P100_candidate_absolute_error": text(candidate_error),
        "P100_scale_spread": text(scale_spread),
        "P100_strut_spread": text(strut_spread),
        "wrong_conversion_residual": text(high["wrong_conversion"]),
        "wrong_sign_residual": text(high["wrong_sign"]),
    }

projector_difference = projector_distance(
    generator_projectors["even"], generator_projectors["odd"]
)
parity_ok = bool(parity_input_ok and projector_difference < mp.mpf("1e-70"))
check("even and odd analytic generator projectors agree", parity_ok,
      f"difference={text(projector_difference, 12)}")
check("the closed-action ratio agrees with the independent endpoint Jacobian", endpoint_bridge_ok)
check("the disclosed P100 vector identifies the analytic line without fitting", candidate_bridge_ok)

rank_records = {}
rank_closure_ok = True
for parity in ("even", "odd"):
    homogeneous = [
        record for record in precision["levels"]["P160"]["parities"][parity]
        if record["homogeneous"]
    ]
    if len(homogeneous) != 1:
        rank_closure_ok = False
        continue
    record = homogeneous[0]
    d_minors = record["D_single_column_deleted_minors"]
    k_minors = record["K_single_column_deleted_minors"]
    d_rank_lower = len(d_minors) == 10 and all(
        not item["contains_zero"] for item in d_minors
    )
    k_rank_lower = len(k_minors) == 15 and all(
        not item["contains_zero"] for item in k_minors
    )
    rank_closure_ok &= bool(record["controls_pass"] and d_rank_lower and k_rank_lower)
    rank_records[parity] = {
        "sector_index": record["sector_index"],
        "D_certified_rank_at_least": 9 if d_rank_lower else None,
        "K_certified_rank_at_least": 14 if k_rank_lower else None,
        "D_exact_rank": 9 if d_rank_lower and symbolic_ok else None,
        "K_exact_rank": 14 if k_rank_lower and symbolic_ok else None,
        "D_exact_nullity": 1 if d_rank_lower and symbolic_ok else None,
        "K_exact_nullity": 1 if k_rank_lower and symbolic_ok else None,
    }
check("the exact line and frozen minors close D/K nullity to one", rank_closure_ok)

controls_ok = bool(
    provenance_ok and formula_ok and corruption_symbolic_ok
    and numeric_convergence_ok and parity_ok
)
if not controls_ok:
    outcome = "HOMOGENEOUS_LINE_CONTROL_FAILED"
elif not symbolic_ok:
    outcome = "HOMOGENEOUS_LINE_SYMBOLIC_DISAGREEMENT"
elif not endpoint_bridge_ok or not candidate_bridge_ok:
    outcome = "HOMOGENEOUS_LINE_NUMERICAL_BRIDGE_DISAGREEMENT"
elif not rank_closure_ok:
    outcome = "HOMOGENEOUS_LINE_NULLITY_OPEN"
else:
    outcome = "HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE"

allowed = {
    "HOMOGENEOUS_LINE_CONTROL_FAILED",
    "HOMOGENEOUS_LINE_SYMBOLIC_DISAGREEMENT",
    "HOMOGENEOUS_LINE_NUMERICAL_BRIDGE_DISAGREEMENT",
    "HOMOGENEOUS_LINE_NULLITY_OPEN",
    "HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE",
}
check("the preregistered homogeneous hierarchy assigns one verdict", outcome in allowed, outcome)

payload = {
    "prior_commit": PRIOR_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "source_sha256": digest(Path(__file__)),
    "symbolic": {
        "closed_action_text_matches": formula_ok,
        "old_momentum_identity": str(correct_momentum),
        "diagonal_response_identity": str(sp.factor(delta_log_diagonal - expected_diagonal)),
        "upper_response_identity": str(sp.factor(delta_log_upper - 2 * SIGMA / LAMBDA)),
        "pole_response_identity": str(sp.factor(delta_log_pole - C)),
        "wrong_conversion_formula": str(wrong_conversion),
        "wrong_sign_formula": str(wrong_sign),
        "dust_old_momentum_derivative": str(sp.diff(DUST, L_MINUS)),
    },
    "bridges": bridge_records,
    "parity_generator_projector_difference": text(projector_difference),
    "rank_closure": rank_records,
    "classification": {
        "nonhomogeneous_canonical_intersection": "ZERO; ADVERSARIALLY REPLICATED",
        "homogeneous_weak_pole_canonical_intersection": (
            "ONE DIMENSION; PRIMARY SYMBOLIC/COMPUTATIONAL"
            if outcome == "HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE" else "OPEN"
        ),
        "omitted_pole_equation": "NOT EVALUATED",
        "gauge_or_physical": "OPEN",
        "tick_c_G_planck_mass": "NOT EVALUATED",
    },
    "outcome": outcome,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"TOTAL: {passed}/{tests} tests PASSED")
print(f"Artifact: {OUTPUT.name}")
if passed != tests:
    raise SystemExit(1)
