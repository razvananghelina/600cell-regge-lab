#!/usr/bin/env python3
"""Corrected-state adversarial audit of classical tick scale covariance.

Original protocol commit: fd1f8a8.
State-correction protocol commit: 8e0f8a6.
Exact state/incidence diagnosis commit: 0fa8947.

The previous failed artifacts are inputs, not overwritten.
"""

import ast
from collections import Counter
import contextlib
from itertools import combinations
import io
import json
import math
from pathlib import Path
import re

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
PRIMARY_ACTION_SOURCE = HERE / "verify_gravity_600cell_dust_canonical_continuation.py"
DIRECT_ACTION_SOURCE = HERE / "verify_gravity_global_boundary_legendre.py"
PRIMARY_SCALE_SOURCE = HERE / "verify_gravity_600cell_tick_scale_covariance.py"
FAILED_ADVERSARIAL_SOURCE = HERE / "verify_gravity_600cell_tick_scale_covariance_adversarial.py"
PRECISION_SOURCE = HERE / "verify_gravity_600cell_tick_scale_covariance_precision.py"
PRIMARY_ARTIFACT = HERE / "gravity_600cell_tick_scale_covariance.json"
INCIDENCE_ARTIFACT = HERE / "gravity_600cell_orbit_action_incidence.json"
OUTPUT = HERE / "gravity_600cell_tick_scale_covariance_state_correction.json"
ORIGINAL_PROTOCOL_COMMIT = "fd1f8a8"
CORRECTION_PROTOCOL_COMMIT = "8e0f8a6"
DIAGNOSIS_COMMIT = "0fa8947"
DPS = 80
mp.mp.dps = DPS
ALPHAS_MP = (mp.mpf(3) / 5, mp.mpf(7) / 4)
ALPHAS_FLOAT = (3.0 / 5.0, 7.0 / 4.0)
SCALE_TOLERANCE = mp.mpf("1e-55")
PRIMARY_TOLERANCE = mp.mpf("1e-45")
HOSTILE_MINIMUM = 1e-8
SUPPORT_MINIMUM = 1e-12
EPSILON = float(np.finfo(float).eps)
OPERATION_FACTOR = 128.0
SIMPLEX_COUNT = 2400
I = mp.mpc(0, 1)
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


def text(value, digits=35):
    return mp.nstr(value, digits)


def relative_error(left, right):
    return abs(left - right) / max(1, abs(left), abs(right))


def definition_prefix(path, cutoff_target):
    tree = ast.parse(path.read_text(), filename=str(path))
    cut = None
    for index, node in enumerate(tree.body):
        if cutoff_target == "first_print":
            if (
                isinstance(node, ast.Expr)
                and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name)
                and node.value.func.id == "print"
            ):
                cut = index
                break
        elif (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == cutoff_target
                for target in node.targets
            )
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError(f"definition cutoff {cutoff_target} absent in {path}")
    namespace = {
        "__file__": str(path),
        "__name__": f"tick_state_correction_{path.stem}",
    }
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(path), "exec"), namespace)
    return namespace


def perturb_signature(path, names):
    tree = ast.parse(path.read_text(), filename=str(path))
    result = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name) or target.id not in names:
            continue
        if not isinstance(node.value, ast.Call) or len(node.value.args) < 3:
            continue
        result[target.id] = (
            ast.literal_eval(node.value.args[-2]),
            ast.literal_eval(node.value.args[-1]),
        )
    return result


def precision_functions():
    wanted = {
        "edge_square",
        "simplex_squared",
        "signed_volume_square",
        "log_minus",
        "angle_record",
        "triangle_area_square",
        "direct_action",
    }
    tree = ast.parse(PRECISION_SOURCE.read_text(), filename=str(PRECISION_SOURCE))
    selected = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    if {node.name for node in selected} != wanted:
        raise RuntimeError("arbitrary-precision direct action definitions incomplete")
    namespace = {
        "mp": mp,
        "I": I,
        "Counter": Counter,
        "combinations": combinations,
    }
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(PRECISION_SOURCE), "exec"), namespace)
    return namespace


def action_key(action):
    return tuple(int(value) for value in action)


def image_vertex(action, vertex):
    return int(action[vertex] if vertex < 120 else action[vertex - 120] + 120)


def image_edge(action, edge):
    return tuple(sorted(image_vertex(action, vertex) for vertex in edge))


print("Corrected-state adversarial tick scale-covariance audit", flush=True)
primary_core = definition_prefix(PRIMARY_ACTION_SOURCE, "first_print")
direct_core = definition_prefix(DIRECT_ACTION_SOURCE, "CONTROLS")
primary_models = primary_core["models"]
direct_models = direct_core["models"]
full_evaluation = direct_core["full_evaluation"]
direct_high = precision_functions()["direct_action"]

primary_signature = perturb_signature(
    PRIMARY_SCALE_SOURCE, {"base_old", "base_internal", "base_new"}
)
failed_signature = perturb_signature(
    FAILED_ADVERSARIAL_SOURCE, {"base_old", "base_variables"}
)
source_diagnosis_ok = bool(
    primary_signature
        == {"base_old": (7, 3), "base_internal": (5, 2), "base_new": (11, 5)}
    and failed_signature == {"base_old": (7, 3), "base_variables": (5, 2)}
)
provenance_ok = bool(
    ORIGINAL_PROTOCOL_COMMIT == "fd1f8a8"
    and CORRECTION_PROTOCOL_COMMIT == "8e0f8a6"
    and DIAGNOSIS_COMMIT == "0fa8947"
    and primary_core["tests"] == primary_core["passed"] == 4
    and direct_core["tests"] == direct_core["passed"] == 7
    and set(primary_models) == set(direct_models) == {"even", "odd"}
)
check(
    "the frozen source structurally exposes the new-boundary perturbation error",
    source_diagnosis_ok,
    f"primary={primary_signature}, failed={failed_signature}",
)
check("the corrected protocol and both action carriers retain their provenance", provenance_ok)


# Independent constants and correctly separated state pieces.
M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON_3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON_3 * L0
TAU = mp.mpf("0.0102")
RHO = TAU**2
SLANT_SQUARE = L0_SQUARE - RHO


def perturbed_mp(values, modulus, center):
    return tuple(
        value * mp.exp(mp.mpf("1e-6") * ((index % modulus) - center))
        for index, value in enumerate(values)
    )


base_old = perturbed_mp(tuple([L0_SQUARE] * 30), 7, 3)
base_internal = perturbed_mp(tuple([SLANT_SQUARE] * 30 + [RHO] * 5), 5, 2)
base_new = perturbed_mp(tuple([L0_SQUARE] * 30), 11, 5)
base_variables = base_internal + base_new


def primary_edge_value(model, edge):
    edge = tuple(sorted(edge))
    if edge in model["old_lookup"]:
        return base_old[model["old_lookup"][edge]]
    if edge in model["edge_to_variable"]:
        index = model["edge_to_variable"][edge]
        return mp.mpf(model["edge_jacobian"][edge]) * base_internal[index]
    if edge in model["final_lookup"]:
        return base_new[model["final_lookup"][edge]]
    raise ValueError(f"edge absent from primary carrier: {edge}")


def direct_edge_value(model, edge):
    edge = tuple(sorted(edge))
    if edge in model["old_to_orbit"]:
        return base_old[model["old_to_orbit"][edge]]
    if edge in model["all_edge_to_variable"]:
        index = model["all_edge_to_variable"][edge]
        return mp.mpf(model["all_edge_jacobian"][edge]) * base_variables[index]
    raise ValueError(f"edge absent from direct carrier: {edge}")


carrier_state_ok = True
for parity in primary_models:
    primary_model = primary_models[parity]
    direct_model = direct_models[parity]
    all_edges = (
        direct_model["old_edges"]
        | direct_model["internal_edges"]
        | direct_model["new_edges"]
    )
    final_direct = {
        edge: variable - 35
        for edge, variable in direct_model["final_to_variable"].items()
    }
    carrier_state_ok &= bool(
        primary_model["slab"] == direct_model["slab"]
        and tuple(sorted(map(action_key, primary_model["stabilizer"])))
            == tuple(sorted(map(action_key, direct_model["stabilizer"])))
        and primary_model["old_edges"] == direct_model["old_edges"]
        and primary_model["internal_edges"] == direct_model["internal_edges"]
        and primary_model["new_edges"] == direct_model["new_edges"]
        and primary_model["old_lookup"] == direct_model["old_to_orbit"]
        and primary_model["edge_to_variable"] == direct_model["edge_to_variable"]
        and primary_model["edge_jacobian"] == direct_model["edge_jacobian"]
        and primary_model["final_lookup"] == final_direct
    )
    for edge in all_edges:
        carrier_state_ok &= bool(
            primary_edge_value(primary_model, edge) == direct_edge_value(direct_model, edge)
        )
    for action in direct_model["stabilizer"]:
        for edge in all_edges:
            carrier_state_ok &= bool(
                direct_edge_value(direct_model, edge)
                == direct_edge_value(direct_model, image_edge(action, edge))
            )
check(
    "the corrected primary and direct states agree on every labelled edge and are stabilizer-invariant",
    carrier_state_ok,
)


incidence = json.loads(INCIDENCE_ARTIFACT.read_text())
incidence_ok = bool(
    incidence.get("outcome") == "ORBIT_ACTION_DIRECT_CONSTRUCTION_SUSPECT"
    and incidence.get("tests") == 12
    and incidence.get("passed") == 9
)
for parity in ("even", "odd"):
    data = incidence.get("enumerations", {}).get(parity, {})
    incidence_ok &= bool(
        len(data.get("triangle_orbit_sizes", ())) == 260
        and len(data.get("simplex_orbit_sizes", ())) == 100
        and len(data.get("flag_orbit_sizes", ())) == 1000
        and sum(data.get("flag_orbit_sizes", ())) == 24000
        and data.get("mismatches") == []
    )
check("the exact flag-incidence artifact has zero shortcut multiplicity mismatches", incidence_ok)


def direct_total_binary64(model, old_values, variables, mass):
    action, gradient, old_gradient, branch = full_evaluation(model, variables, old_values)
    dust_action = -(8.0 * math.pi * mass / 5.0) * float(
        np.sum(np.sqrt(variables[30:35]))
    )
    dust_gradient = np.zeros(65, dtype=complex)
    dust_gradient[30:35] = (
        -(4.0 * math.pi * mass / 5.0) / np.sqrt(variables[30:35])
    )
    total_action = action + dust_action
    total_gradient = gradient + dust_gradient
    imaginary = max(
        abs(total_action.imag),
        float(np.max(np.abs(total_gradient.imag))),
        float(np.max(np.abs(old_gradient.imag))),
    )
    branch_ok = bool(
        branch["negative_counts"] == Counter({1: SIMPLEX_COUNT})
        and branch["minimum_gram"] > 1e-8
        and branch["minimum_argument"] > 1e-6
        and imaginary < 2e-7
    )
    return {
        "action": total_action,
        "gradient": total_gradient,
        "old_gradient": old_gradient,
        "branch": branch,
        "imaginary": imaginary,
        "branch_ok": branch_ok,
    }


def propagated_bound(*evaluations):
    minimum_gram = min(item["branch"]["minimum_gram"] for item in evaluations)
    return OPERATION_FACTOR * EPSILON * SIMPLEX_COUNT * max(1.0, 1.0 / minimum_gram)


old_float = np.array([float(value) for value in base_old])
variables_float = np.array([float(value) for value in base_variables])
mass_float = float(MASS)
binary_records = {}
binary_ok = True
hostile_ok = True
support_ok = True
envelope_ok = True
for parity, model in direct_models.items():
    base = direct_total_binary64(model, old_float, variables_float, mass_float)
    internal_support = int(np.sum(np.abs(base["gradient"][:35]) > SUPPORT_MINIMUM))
    boundary_support = int(
        np.sum(np.abs(base["gradient"][35:]) > SUPPORT_MINIMUM)
        + np.sum(np.abs(base["old_gradient"]) > SUPPORT_MINIMUM)
    )
    support_ok &= bool(
        base["branch_ok"] and internal_support >= 20 and boundary_support >= 20
    )
    parity_record = {
        "base_action_real": float(base["action"].real),
        "base_action_imaginary": float(base["action"].imag),
        "internal_support": internal_support,
        "boundary_support": boundary_support,
        "scales": {},
    }
    for alpha in ALPHAS_FLOAT:
        r = alpha**2
        scaled = direct_total_binary64(
            model, r * old_float, r * variables_float, alpha * mass_float
        )
        fixed = direct_total_binary64(
            model, r * old_float, r * variables_float, mass_float
        )
        bound = propagated_bound(base, scaled, fixed)
        action_error = relative_error(scaled["action"], r * base["action"])
        gradient_error = max(
            relative_error(left, right)
            for left, right in zip(scaled["gradient"], base["gradient"])
        )
        old_error = max(
            relative_error(left, right)
            for left, right in zip(scaled["old_gradient"], base["old_gradient"])
        )
        fixed_action_error = relative_error(fixed["action"], r * base["action"])
        fixed_pole_error = max(
            relative_error(fixed["gradient"][index], base["gradient"][index])
            for index in range(30, 35)
        )
        one_binary_ok = bool(
            scaled["branch_ok"]
            and action_error <= bound
            and gradient_error <= bound
            and old_error <= bound
        )
        one_hostile_ok = bool(
            fixed["branch_ok"]
            and fixed_action_error > max(HOSTILE_MINIMUM, 100 * bound)
            and fixed_pole_error > max(HOSTILE_MINIMUM, 100 * bound)
        )
        binary_ok &= one_binary_ok
        hostile_ok &= one_hostile_ok
        envelope_ok &= bound < 1e-3
        parity_record["scales"][format(alpha, ".17g")] = {
            "bound": bound,
            "action_error": action_error,
            "gradient_error": gradient_error,
            "old_gradient_error": old_error,
            "fixed_mass_action_error": fixed_action_error,
            "fixed_mass_pole_error": fixed_pole_error,
            "covariance_pass": one_binary_ok,
            "hostile_pass": one_hostile_ok,
        }
        check(
            f"{parity}, alpha={alpha:.8g}: corrected direct binary64 action and all raw derivatives are covariant",
            one_binary_ok and bound < 1e-3,
            f"bound={bound:.3e}, action={action_error:.3e}, grad={gradient_error:.3e}, old={old_error:.3e}",
        )
        check(
            f"{parity}, alpha={alpha:.8g}: corrected fixed-mass hostile control fails covariance",
            one_hostile_ok,
            f"action={fixed_action_error:.3e}, pole={fixed_pole_error:.3e}",
        )
    binary_records[parity] = {"evaluation": base, "artifact": parity_record}

check(
    "the corrected off-shell state has nonzero internal and boundary raw-gradient support",
    support_ok,
)


# All direct arbitrary-precision actions are completed before reading primary data.
high_records = {parity: {} for parity in direct_models}
high_branch_ok = True
high_scale_ok = True
for parity, model in direct_models.items():
    for alpha in (mp.mpf(1),) + ALPHAS_MP:
        print(f"  evaluating direct 80d {parity}, alpha={text(alpha, 8)}", flush=True)
        r = alpha**2
        record = direct_high(
            model,
            tuple(r * value for value in base_old),
            tuple(r * value for value in base_variables),
            alpha * MASS,
        )
        state_branch_ok = bool(
            record["negative_counts"] == Counter({1: SIMPLEX_COUNT})
            and record["minimum_leading_minor"] > mp.mpf("1e-20")
            and record["minimum_argument"] > mp.mpf("1e-6")
            and record["maximum_imaginary"] < mp.mpf("1e-60")
        )
        high_branch_ok &= state_branch_ok
        high_records[parity][alpha] = record
    base_high = high_records[parity][mp.mpf(1)]["action"]
    for alpha in ALPHAS_MP:
        scale_error = relative_error(
            high_records[parity][alpha]["action"], alpha**2 * base_high
        )
        high_records[parity][alpha]["scale_error"] = scale_error
        high_scale_ok &= scale_error < SCALE_TOLERANCE
        check(
            f"{parity}, alpha={text(alpha, 8)}: corrected direct 80-decimal action has degree two",
            scale_error < SCALE_TOLERANCE,
            f"relative error={text(scale_error, 8)}",
        )

check("all six corrected direct high-precision states retain the Lorentzian branch", high_branch_ok)


def parse_mpc(value):
    match = re.fullmatch(r"\((\S+)\s+([+-])\s+(\S+)j\)", value.strip())
    if match is None:
        raise ValueError(f"cannot parse mpc text: {value}")
    imaginary = mp.mpf(match.group(3))
    if match.group(2) == "-":
        imaginary = -imaginary
    return mp.mpc(mp.mpf(match.group(1)), imaginary)


primary_artifact = json.loads(PRIMARY_ARTIFACT.read_text())
primary_artifact_ok = bool(
    primary_artifact.get("outcome") == "TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED"
    and primary_artifact.get("tests") == primary_artifact.get("passed") == 12
)
agreement_ok = True
agreement_records = {}
for parity in direct_models:
    stored = parse_mpc(primary_artifact["parities"][parity]["base"]["action"])
    high = high_records[parity][mp.mpf(1)]["action"]
    direct_primary_error = relative_error(high, stored)
    fresh_action, _, _ = primary_core["action_and_gradient"](
        primary_models[parity], base_old, base_internal, base_new
    )
    fresh_primary_error = relative_error(fresh_action, stored)
    binary_action = mp.mpc(
        binary_records[parity]["evaluation"]["action"].real,
        binary_records[parity]["evaluation"]["action"].imag,
    )
    binary_high_error = relative_error(binary_action, high)
    maximum_bound = max(
        item["bound"] for item in binary_records[parity]["artifact"]["scales"].values()
    )
    one_ok = bool(
        direct_primary_error < PRIMARY_TOLERANCE
        and fresh_primary_error < PRIMARY_TOLERANCE
        and binary_high_error <= maximum_bound
    )
    agreement_ok &= one_ok
    agreement_records[parity] = {
        "direct_high_vs_stored_primary_error": text(direct_primary_error, 15),
        "fresh_primary_vs_stored_primary_error": text(fresh_primary_error, 15),
        "binary64_direct_vs_high_precision_error": text(binary_high_error, 15),
        "binary64_maximum_propagated_bound": maximum_bound,
        "pass": one_ok,
    }
    check(
        f"{parity}: corrected direct, fresh primary and stored primary actions agree",
        one_ok,
        f"direct/primary={text(direct_primary_error, 6)}, fresh/stored={text(fresh_primary_error, 6)}, binary/high={text(binary_high_error, 6)}",
    )

controls_ok = bool(
    source_diagnosis_ok
    and provenance_ok
    and carrier_state_ok
    and incidence_ok
    and support_ok
    and envelope_ok
    and hostile_ok
    and high_branch_ok
    and primary_artifact_ok
)
if not controls_ok:
    outcome = "TICK_SCALE_STATE_CORRECTION_CONTROL_FAILED"
elif not (binary_ok and high_scale_ok and agreement_ok):
    outcome = "TICK_SCALE_CORRECTED_IMPLEMENTATIONS_DISAGREE"
else:
    outcome = "ABSOLUTE_CLASSICAL_TICK_NO_GO_ADVERSARIALLY_CORROBORATED"

check(
    "the corrected adversarial audit corroborates the conditional absolute classical tick no-go",
    outcome == "ABSOLUTE_CLASSICAL_TICK_NO_GO_ADVERSARIALLY_CORROBORATED",
    outcome,
)

artifact = {
    "title": "Corrected-state adversarial audit of classical tick scale covariance",
    "original_protocol_commit": ORIGINAL_PROTOCOL_COMMIT,
    "correction_protocol_commit": CORRECTION_PROTOCOL_COMMIT,
    "diagnosis_commit": DIAGNOSIS_COMMIT,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "source_diagnosis": {
        "primary": primary_signature,
        "failed_adversarial": failed_signature,
        "pass": source_diagnosis_ok,
    },
    "state": {
        "old": "exp(1e-6*((i mod 7)-3))",
        "internal": "exp(1e-6*((i mod 5)-2))",
        "new": "exp(1e-6*((i mod 11)-5))",
        "labelled_carrier_equality_and_invariance": carrier_state_ok,
    },
    "incidence_control": incidence_ok,
    "binary64": {
        parity: data["artifact"] for parity, data in binary_records.items()
    },
    "high_precision": {
        parity: {
            "states": {
                text(alpha, 20): {
                    "action": text(record["action"], 65),
                    "minimum_leading_minor": text(record["minimum_leading_minor"], 20),
                    "minimum_argument": text(record["minimum_argument"], 20),
                    "maximum_imaginary": text(record["maximum_imaginary"], 12),
                    "scale_error": text(record.get("scale_error", mp.mpf(0)), 15),
                }
                for alpha, record in records.items()
            }
        }
        for parity, records in high_records.items()
    },
    "agreement": agreement_records,
    "interpretation": {
        "derived_exact_adversarially_corroborated": (
            "under q->alpha^2 q and mass->alpha mass, the stated scale-free "
            "classical Regge-dust action and equations occur in global scale families"
        ),
        "derived_negative": (
            "the current classical theory cannot select an absolute nonzero tick"
        ),
        "not_excluded": [
            "tau/L",
            "tau_next/tau0",
            "relational dust time",
            "independently motivated scale-breaking physics",
        ],
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

print(f"RESULT: {passed}/{tests}", flush=True)
print(f"OUTCOME: {outcome}", flush=True)
if passed != tests:
    raise SystemExit(1)
