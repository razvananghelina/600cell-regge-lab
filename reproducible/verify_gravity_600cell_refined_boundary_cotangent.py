#!/usr/bin/env python3
"""All-schedule boundary covector of the curvature-matched refined seed."""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
CURVATURE_ADVERSARIAL = HERE / "gravity_600cell_refined_local_curvature_mass_adversarial.json"
COTANGENT = HERE / "gravity_600cell_refined_homogeneous_cotangent_lift.json"
COTANGENT_SOURCE = HERE / "verify_gravity_600cell_refined_homogeneous_cotangent_lift.py"
COARSE_IDENTITY = HERE / "gravity_600cell_dust_regular_lapse_identity.json"
ACCELERATION = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_blind.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_refined_boundary_cotangent_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_boundary_cotangent_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_boundary_cotangent.json"

PRIOR_ART_COMMIT = "e7a1545"
PROTOCOL_COMMIT = "5ee4f3e"
EXPECTED_HASHES = {
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "curvature": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "curvature_adversarial": "c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e",
    "cotangent": "93dd857bff3b406e86d41a8a4b05d6441cb0e3e1c11e4f53d098555b1218924b",
    "cotangent_source": "154081b12f74ed8597a4b72b37a99219d64c0905829da1566e840fc562b1c20c",
    "coarse_identity": "5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51",
    "acceleration": "2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2",
    "prior_art": "12a32f97129e7fa388261ffd5b16fa4463778ec67af82f676d0365d96300a672",
    "protocol": "38ccd8ab639c974da5bf58235582620a261f5e5b55ff4989ba14c9c7147cd6a6",
}
PAIR4 = tuple(combinations(range(4), 2))
LOCAL_TRIANGLES = np.asarray(tuple(combinations(range(5), 3)), dtype=np.int8)
VARIABLES = (
    tuple(("old",) + pair for pair in PAIR4)
    + tuple(("new",) + pair for pair in PAIR4)
    + tuple(("cross",) + pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
INTERNAL_VARIABLES = (
    tuple(("cross",) + pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
TAU_TEXT = "0.0102"

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
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=75):
    return mp.nstr(value, digits, strip_zeros=False)


def load_action_definitions():
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "mp": mp,
        "np": np,
        "json": json,
        "Path": Path,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "HERE": HERE,
        "ROOT": ROOT,
        "PAIR4": PAIR4,
        "LOCAL_TRIANGLES": LOCAL_TRIANGLES,
        "TAU_TEXT": TAU_TEXT,
        "VARIABLES": VARIABLES,
        "INTERNAL_VARIABLES": INTERNAL_VARIABLES,
        "FD_STEP_TEXTS": ("1e-15", "5e-16"),
        "FD_GATE_TEXT": "1e-24",
        "EXPECTED_F": (2640, 17040, 28800, 14400),
        "tests": 0,
        "passed": 0,
    }
    exec(
        compile(ast.Module(body=definitions, type_ignores=[]), str(ACTION_SOURCE), "exec"),
        namespace,
    )
    return namespace


paths = {
    "action_source": ACTION_SOURCE,
    "curvature": CURVATURE,
    "curvature_adversarial": CURVATURE_ADVERSARIAL,
    "cotangent": COTANGENT,
    "cotangent_source": COTANGENT_SOURCE,
    "coarse_identity": COARSE_IDENTITY,
    "acceleration": ACCELERATION,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all action, curvature, cotangent and coarse inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIOR_ART_COMMIT == "e7a1545"
    and PROTOCOL_COMMIT == "5ee4f3e",
)

curvature = json.loads(CURVATURE.read_text())
curvature_adversarial = json.loads(CURVATURE_ADVERSARIAL.read_text())
cotangent = json.loads(COTANGENT.read_text())
coarse_identity = json.loads(COARSE_IDENTITY.read_text())
acceleration = json.loads(ACCELERATION.read_text())
upstream_ok = check(
    "all frozen upstream artifacts carry their accepted scoped outcomes",
    curvature["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and curvature["tests"] == {"passed": 15, "total": 15}
    and curvature_adversarial["outcome"]
        == "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"
    and curvature_adversarial["tests"] == {"passed": 16, "total": 16}
    and cotangent["outcome"]
        == "REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNDERDETERMINED"
    and cotangent["orbit_total"]["rank"] == 1
    and cotangent["orbit_total"]["nullity"] == 5
    and coarse_identity["outcome"] == "REGULAR_LAPSE_IDENTITY_PROVED"
    and coarse_identity["passed"] == coarse_identity["tests"] == 13
    and acceleration["outcome"]
        == "CANONICAL_CARRIER_ACCELERATION_COEFFICIENTS_DERIVED"
    and acceleration["passed"] == acceleration["tests"] == 10,
)

pullback = tuple(
    mp.mpf(value) for value in cotangent["orbit_total"]["pullback"][0]
)
pullback_ok = check(
    "the frozen homogeneous orbit-total pullback is exactly (2,2,2,2,2,2)",
    pullback == (mp.mpf(2),) * 6,
)

defs = load_action_definitions()
definition_ok = check(
    "only the frozen stationary-fill action definitions were loaded",
    all(name in defs for name in (
        "tetrahedra_from_adjacency",
        "barycentric_chambers",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "evaluate_schedule",
    ))
    and "OUTPUT" not in defs,
)

_, adjacency, _ = build_600cell()
coarse_top = defs["tetrahedra_from_adjacency"](adjacency)
_, top, colours = defs["barycentric_chambers"](coarse_top)
spatial_cells = defs["all_simplices"](tuple(map(tuple, top)))
topology_ok = check(
    "the direct carrier rebuild has the exact projected barycentric f-vector",
    tuple(len(layer) for layer in spatial_cells)
        == (2640, 17040, 28800, 14400),
)

orders = tuple(permutations(range(4)))
selected_masses_text = curvature["selected_rank_matter"]["total_masses"]
schedule_records = []
retained_combinatorics = {}
maximum_cross_residual = mp.mpf(0)
maximum_vertical_residual = mp.mpf(0)
maximum_branch_identity = mp.mpf(0)
maximum_branch_imaginary = mp.mpf(0)
minimum_branch_argument = mp.inf
maximum_boundary_imaginary = mp.mpf(0)

print("[INFO] evaluating all 24 direct boundary covectors", flush=True)
with mp.workdps(100):
    tau = mp.mpf(TAU_TEXT)
    selected_masses = tuple(mp.mpf(value) for value in selected_masses_text)
    geometry = defs["exact_geometry"](100)
    geometry["mass"] = mp.mpf(0)
    coordinates = defs["base_coordinates"](geometry)
    for index, order in enumerate(orders):
        combinatorics = defs["schedule_combinatorics"](top, colours, order)
        evaluation = defs["evaluate_schedule"](combinatorics, geometry, coordinates)
        if order in ((0, 1, 2, 3), (3, 2, 1, 0)):
            retained_combinatorics[order] = combinatorics
        cross = tuple(evaluation["gradient"][("cross",) + pair] for pair in PAIR4)
        vertical = tuple(
            evaluation["gradient"]["rho", rank]
            - 4 * mp.pi * selected_masses[rank] * tau
            for rank in range(4)
        )
        old_gradient = tuple(
            evaluation["gradient"][("old",) + pair] for pair in PAIR4
        )
        new_gradient = tuple(
            evaluation["gradient"][("new",) + pair] for pair in PAIR4
        )
        pre = tuple(-value for value in old_gradient)
        post = new_gradient
        schedule_records.append({
            "order": order,
            "pre": pre,
            "post": post,
            "old_gradient": old_gradient,
            "new_gradient": new_gradient,
        })
        maximum_cross_residual = max(
            maximum_cross_residual, *(abs(value) for value in cross)
        )
        maximum_vertical_residual = max(
            maximum_vertical_residual, *(abs(value) for value in vertical)
        )
        maximum_branch_identity = max(
            maximum_branch_identity,
            evaluation["maximum_angle_identity_residual"],
        )
        maximum_branch_imaginary = max(
            maximum_branch_imaginary,
            evaluation["maximum_imaginary_curvature"],
        )
        minimum_branch_argument = min(
            minimum_branch_argument, evaluation["minimum_angle_argument"]
        )
        maximum_boundary_imaginary = max(
            maximum_boundary_imaginary,
            *(abs(mp.im(value)) for value in old_gradient + new_gradient),
        )
        if index in (5, 11, 17, 23):
            print(f"[INFO] schedules completed: {index + 1}/24", flush=True)

on_shell_ok = check(
    "curvature masses make all direct internal orbit equations stationary",
    maximum_cross_residual < mp.mpf("1e-60")
    and maximum_vertical_residual < mp.mpf("1e-60"),
    f"cross={mp_text(maximum_cross_residual, 8)}, "
    f"vertical={mp_text(maximum_vertical_residual, 8)}",
)
branch_ok = check(
    "all boundary evaluations stay finite on the accepted Lorentzian branch",
    maximum_branch_identity < mp.mpf("1e-80")
    and maximum_branch_imaginary < mp.mpf("1e-80")
    and maximum_boundary_imaginary < mp.mpf("1e-80")
    and minimum_branch_argument > mp.mpf("1e-8"),
    f"identity={mp_text(maximum_branch_identity, 8)}, "
    f"imag={mp_text(maximum_branch_imaginary, 8)}, "
    f"boundary imag={mp_text(maximum_boundary_imaginary, 8)}",
)

with mp.workdps(100):
    reference_pre = schedule_records[0]["pre"]
    reference_post = schedule_records[0]["post"]
    pre_spread = max(
        abs(record["pre"][index] - reference_pre[index])
        for record in schedule_records for index in range(6)
    )
    post_spread = max(
        abs(record["post"][index] - reference_post[index])
        for record in schedule_records for index in range(6)
    )
    reflection_error = max(
        abs(reference_post[index] + reference_pre[index]) for index in range(6)
    )
schedule_ok = check(
    "all 24 schedules select one complete six-component pre/post covector",
    pre_spread < mp.mpf("1e-60") and post_spread < mp.mpf("1e-60"),
    f"pre spread={mp_text(pre_spread, 8)}, post spread={mp_text(post_spread, 8)}",
)
reflection_ok = check(
    "the selected post covector is the negative of the pre covector",
    reflection_error < mp.mpf("1e-60"),
    f"max error={mp_text(reflection_error, 8)}",
)

finite_difference_errors = []
with mp.workdps(100):
    for order in ((0, 1, 2, 3), (3, 2, 1, 0)):
        combinatorics = retained_combinatorics[order]
        evaluation = defs["evaluate_schedule"](combinatorics, geometry, coordinates)
        for layer in ("old", "new"):
            for pair in PAIR4:
                key = (layer,) + pair
                derivatives = []
                for step_text in ("1e-12", "5e-13"):
                    step = mp.mpf(step_text)
                    plus = dict(coordinates)
                    minus = dict(coordinates)
                    plus[key] *= mp.exp(step)
                    minus[key] *= mp.exp(-step)
                    plus_action = defs["evaluate_schedule"](
                        combinatorics, geometry, plus
                    )["action"]
                    minus_action = defs["evaluate_schedule"](
                        combinatorics, geometry, minus
                    )["action"]
                    derivatives.append((plus_action - minus_action) / (2 * step))
                richardson = (4 * derivatives[1] - derivatives[0]) / 3
                analytic = evaluation["gradient"][key]
                finite_difference_errors.append(
                    abs(richardson - analytic) / max(mp.mpf(1), abs(analytic))
                )
    finite_difference_error = max(finite_difference_errors)
finite_difference_ok = check(
    "independent Richardson checks reproduce all 24 selected boundary gradients",
    finite_difference_error < mp.mpf("1e-36"),
    f"max relative error={mp_text(finite_difference_error, 8)}",
)

with mp.workdps(100):
    p_pre_fine = mp.fsum(
        pullback[index] * reference_pre[index] for index in range(6)
    )
    p_post_fine = mp.fsum(
        pullback[index] * reference_post[index] for index in range(6)
    )
    k_fine = mp.mpf(curvature["geometry"]["total_curvature"])
    m_fine = mp.mpf(curvature["selected_rank_matter"]["total_mass"])
    fine_pre_error = abs(p_pre_fine + tau * k_fine / 2)
    fine_post_error = abs(p_post_fine - tau * k_fine / 2)
    fine_mass_error = max(
        abs(p_pre_fine + 4 * mp.pi * tau * m_fine),
        abs(p_post_fine - 4 * mp.pi * tau * m_fine),
    )
fine_identity_ok = check(
    "the refined pullbacks equal plus/minus tau*K_fine/2 and 4*pi*tau*M_fine",
    max(fine_pre_error, fine_post_error, fine_mass_error) < mp.mpf("1e-60"),
    f"max error={mp_text(max(fine_pre_error, fine_post_error, fine_mass_error), 8)}",
)

with mp.workdps(100):
    zeta = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
    epsilon3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
    k_coarse = 720 * zeta * epsilon3
    m_coarse = k_coarse / (8 * mp.pi)
    p_pre_coarse = -tau * k_coarse / 2
    p_post_coarse = +tau * k_coarse / 2
    raw_ratio = p_pre_fine / p_pre_coarse
    curvature_ratio = k_fine / k_coarse
    mass_ratio = m_fine / m_coarse
    normalized_ratio = (p_pre_fine / m_fine) / (p_pre_coarse / m_coarse)
    ratio_error = max(
        abs(raw_ratio - curvature_ratio),
        abs(raw_ratio - mass_ratio),
        abs(normalized_ratio - 1),
    )
    raw_mismatch = abs(raw_ratio - 1)
ratio_ok = check(
    "raw momentum tracks finite curvature/mass renormalization and normalized momentum agrees",
    ratio_error < mp.mpf("1e-60") and raw_mismatch > mp.mpf("1e-4"),
    f"raw={mp_text(raw_ratio, 22)}, normalized={mp_text(normalized_ratio, 22)}",
)

with mp.workdps(100):
    corrupted_pre = list(reference_pre)
    corrupted_pre[0] += mp.mpf("1e-10")
    corrupted_spread = max(
        abs(corrupted_pre[index] - schedule_records[1]["pre"][index])
        for index in range(6)
    )
    corrupted_pullback = mp.fsum(
        pullback[index] * corrupted_pre[index] for index in range(6)
    )
corruption_ok = check(
    "a deliberate boundary-component corruption fails vector and pullback gates",
    corrupted_spread > mp.mpf("1e-10") / 2
    and abs(corrupted_pullback - p_pre_fine) > mp.mpf("1e-10"),
)

with mp.workdps(100):
    kernel_perturbation = (
        mp.mpf("1e-6"), -mp.mpf("1e-6"), mp.mpf(0),
        mp.mpf(0), mp.mpf(0), mp.mpf(0),
    )
    kernel_vector = tuple(
        reference_pre[index] + kernel_perturbation[index] for index in range(6)
    )
    kernel_pullback = mp.fsum(
        pullback[index] * kernel_vector[index] for index in range(6)
    )
    kernel_distance = max(
        abs(kernel_vector[index] - reference_pre[index]) for index in range(6)
    )
    kernel_control_condition = (
        abs(kernel_pullback - p_pre_fine) < mp.mpf("1e-90")
        and abs(kernel_distance - mp.mpf("1e-6")) < mp.mpf("1e-90")
    )
kernel_control_ok = check(
    "a distinct kernel-shifted covector has exactly the same scalar pullback",
    kernel_control_condition,
)

with mp.workdps(100):
    uniform_vector = (p_pre_fine / 12,) * 6
    uniform_distance = max(
        abs(uniform_vector[index] - reference_pre[index]) for index in range(6)
    )

scope = {
    "stationary_root_or_nested_census_executed": False,
    "hessian_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
}
scope_ok = check(
    "the calculation stays inside the boundary-cotangent scope",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    pullback_ok,
    definition_ok,
    topology_ok,
    on_shell_ok,
    branch_ok,
    reflection_ok,
    finite_difference_ok,
    ratio_ok,
    corruption_ok,
    kernel_control_ok,
    scope_ok,
))
if not controls_ok:
    outcome = "REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED"
elif not schedule_ok:
    outcome = "REFINED_BOUNDARY_COTANGENT_SCHEDULE_DEPENDENT"
elif not fine_identity_ok:
    outcome = "REFINED_BOUNDARY_COTANGENT_IDENTITY_REFUTED"
elif raw_mismatch < mp.mpf("1e-60"):
    outcome = "REFINED_BOUNDARY_COTANGENT_SELECTED_EXACT_COARSE"
else:
    outcome = "REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one boundary-cotangent outcome",
    outcome in {
        "REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED",
        "REFINED_BOUNDARY_COTANGENT_SCHEDULE_DEPENDENT",
        "REFINED_BOUNDARY_COTANGENT_IDENTITY_REFUTED",
        "REFINED_BOUNDARY_COTANGENT_SELECTED_EXACT_COARSE",
        "REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED",
    },
    outcome,
)

artifact_schedules = []
for record in schedule_records:
    artifact_schedules.append({
        "order": list(record["order"]),
        "pre": [mp_text(mp.re(value)) for value in record["pre"]],
        "post": [mp_text(mp.re(value)) for value in record["post"]],
    })

artifact = {
    "title": "All-schedule boundary cotangent of the refined on-shell seed",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "coordinate": "log squared edge length",
        "orbit_order": [f"{left}{right}" for left, right in PAIR4],
        "pullback": [mp_text(value) for value in pullback],
        "pre_sign": "P_pre=-dS/dlog(q_old)",
        "post_sign": "P_post=+dS/dlog(q_new)",
        "tau0": TAU_TEXT,
    },
    "on_shell_and_branch": {
        "maximum_cross_residual": mp_text(maximum_cross_residual),
        "maximum_vertical_residual": mp_text(maximum_vertical_residual),
        "maximum_angle_identity_residual": mp_text(maximum_branch_identity),
        "maximum_imaginary_curvature": mp_text(maximum_branch_imaginary),
        "minimum_angle_argument": mp_text(minimum_branch_argument),
        "maximum_boundary_imaginary": mp_text(maximum_boundary_imaginary),
    },
    "boundary_covectors": {
        "schedules": artifact_schedules,
        "selected_pre": [mp_text(mp.re(value)) for value in reference_pre],
        "selected_post": [mp_text(mp.re(value)) for value in reference_post],
        "maximum_pre_schedule_spread": mp_text(pre_spread),
        "maximum_post_schedule_spread": mp_text(post_spread),
        "time_reflection_error": mp_text(reflection_error),
        "maximum_richardson_relative_error": mp_text(finite_difference_error),
        "uniform_same_pullback_distance": mp_text(uniform_distance),
    },
    "pullback_comparison": {
        "p_pre_fine": mp_text(mp.re(p_pre_fine)),
        "p_post_fine": mp_text(mp.re(p_post_fine)),
        "K_fine": mp_text(k_fine),
        "M_fine": mp_text(m_fine),
        "p_pre_coarse": mp_text(p_pre_coarse),
        "p_post_coarse": mp_text(p_post_coarse),
        "K_coarse": mp_text(k_coarse),
        "M_coarse": mp_text(m_coarse),
        "raw_ratio": mp_text(mp.re(raw_ratio)),
        "curvature_ratio": mp_text(curvature_ratio),
        "mass_ratio": mp_text(mass_ratio),
        "normalized_ratio": mp_text(mp.re(normalized_ratio)),
        "raw_fractional_mismatch": mp_text(raw_mismatch),
        "maximum_ratio_identity_error": mp_text(ratio_error),
    },
    "controls": {
        "corrupted_vector_spread": mp_text(corrupted_spread),
        "corrupted_pullback_change": mp_text(abs(corrupted_pullback - p_pre_fine)),
        "kernel_shift_distance": mp_text(kernel_distance),
        "kernel_shift_pullback_change": mp_text(abs(kernel_pullback - p_pre_fine)),
    },
    "scope": scope,
    "status_labels": {
        "complete_refined_covector": "COMPUTED_BEFORE_SCALAR_COMPARISON",
        "fixed_radius_bare_transport": "TESTED_SEPARATELY",
        "mass_normalized_transport": "TESTED_SEPARATELY",
        "perfect_action": "NOT_ESTABLISHED",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print("P_pre=" + ", ".join(mp_text(mp.re(value), 18) for value in reference_pre))
print(f"raw ratio={mp_text(mp.re(raw_ratio), 18)}")
print(f"normalized ratio={mp_text(mp.re(normalized_ratio), 18)}")

raise SystemExit(0 if passed == tests else 1)
