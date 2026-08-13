#!/usr/bin/env python3
"""Frozen 35-variable stationary-root search on the global Regge slab.

Protocol commit: 5b687a3.  The six starts, causal gates, logarithmic variables,
Levenberg--Marquardt damping list, stopping rules and validation tolerances were
fixed before any root iteration was run.
"""

import contextlib
from collections import Counter
import importlib.util
import io
import json
import math
import multiprocessing as mp
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_global_regge_roots.json"
PROTOCOL_COMMIT = "5b687a3"
CORRECTION_COMMIT = "82ab7b8"
CLASSIFICATION_CORRECTION_COMMIT = "6594398"
EVALUATOR_COMMIT = "d9fe159"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


# Import the certified evaluator as a library.  Its module-level audit still
# runs, but suppress its console transcript; its own counter is checked below.
evaluator_path = HERE / "verify_gravity_global_regge_orbits.py"
spec = importlib.util.spec_from_file_location("global_regge_orbits", evaluator_path)
gro = importlib.util.module_from_spec(spec)
upstream_transcript = io.StringIO()
with contextlib.redirect_stdout(upstream_transcript):
    spec.loader.exec_module(gro)


check(
    "the imported full-versus-orbit evaluator retains all 43 certificates",
    gro.tests == gro.passed == 43,
)
check(
    "both parity models expose exactly 35 honest internal-edge orbits",
    set(gro.models) == {"even", "odd"}
    and all(len(model["edge_orbits"]) == 35 for model in gro.models.values()),
)


def frozen_starts():
    j = np.arange(1, 31, dtype=float)
    k = np.arange(1, 6, dtype=float)
    return {
        "S0": np.r_[np.ones(30), np.full(5, 0.25)],
        "S1": np.r_[1+j/1000, 0.25+k/1000],
        "S2": np.r_[1-j/2000, 0.25+(6-k)/1500],
        "S3": np.r_[
            1+0.02*np.cos(2*np.pi*j/31),
            0.25+0.01*np.cos(2*np.pi*k/5),
        ],
        "S4": np.r_[
            1+0.02*np.sin(2*np.pi*j/31),
            0.25+0.01*np.sin(2*np.pi*k/5),
        ],
        "S5": np.r_[
            1+0.01*(-1.0)**(j-1),
            0.25+0.005*(-1.0)**(k-1),
        ],
    }


STARTS = frozen_starts()
LOWER_Y = -6.0
UPPER_Y = 6.0
JACOBIAN_STEP = 2e-5
DAMPING = (0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1.0, 100.0)
MAX_ACCEPTED_ITERATIONS = 80


def state(model, y, need_full=False):
    """Evaluate one logarithmic point and all frozen causal gates."""
    if not np.all(np.isfinite(y)) or np.any(y <= LOWER_Y) or np.any(y >= UPPER_Y):
        return {"valid": False, "reason": "logarithmic box"}
    x = np.exp(y)
    try:
        action, orbit_gradient, data = gro.reduced_evaluation(model, x)
    except Exception as error:
        return {"valid": False, "reason": f"evaluator: {type(error).__name__}: {error}"}
    residual = orbit_gradient/24
    action_imaginary = abs(action.imag)/max(1.0, abs(action.real))
    gradient_imaginary = float(np.max(np.abs(residual.imag)))/max(
        1.0, float(np.max(np.abs(residual.real)))
    )
    causal = (
        data["negative_counts"] == gro.Counter({1: 100})
        and data["minimum_absolute_gram_eigenvalue"] > 1e-8
        and data["minimum_argument"] > 1e-6
        and action_imaginary < 2e-7
        and gradient_imaginary < 2e-7
    )
    result = {
        "valid": bool(causal),
        "reason": "ok" if causal else "causal/branch gate",
        "x": x,
        "y": y,
        "action": action,
        "orbit_gradient": orbit_gradient,
        "residual": residual.real,
        "residual_inf": float(np.linalg.norm(residual.real, ord=np.inf)),
        "residual_2": float(np.linalg.norm(residual.real)),
        "action_imaginary_relative": float(action_imaginary),
        "gradient_imaginary_relative": float(gradient_imaginary),
        "minimum_absolute_gram_eigenvalue": data["minimum_absolute_gram_eigenvalue"],
        "minimum_angle_argument": data["minimum_argument"],
    }
    if need_full and causal:
        full_action, full_gradient, full_data = gro.full_evaluation(model, x)
        result["full_action"] = full_action
        result["full_restricted_gradient"] = full_gradient
        result["full_data"] = full_data
    return result


_WORKER_MODEL = None


def initialize_worker(model):
    global _WORKER_MODEL
    _WORKER_MODEL = model


def worker_state(y):
    return state(_WORKER_MODEL, y)


def jacobian_y(model, y, pool=None):
    """Frozen centered Jacobian of the per-edge residual in log variables."""
    jacobian = np.empty((35, 35))
    points = []
    for variable in range(35):
        plus = y.copy()
        minus = y.copy()
        plus[variable] += JACOBIAN_STEP
        minus[variable] -= JACOBIAN_STEP
        points.extend((plus, minus))
    evaluated = (
        pool.map(worker_state, points, chunksize=2)
        if pool is not None else [state(model, point) for point in points]
    )
    for variable in range(35):
        plus_state = evaluated[2*variable]
        minus_state = evaluated[2*variable+1]
        if not plus_state["valid"] or not minus_state["valid"]:
            reasons = {plus_state["reason"], minus_state["reason"]}
            if "logarithmic box" in reasons:
                return None, "artificial box boundary"
            return None, "noncausal centered Jacobian point"
        jacobian[:, variable] = (
            plus_state["residual"]-minus_state["residual"]
        )/(2*JACOBIAN_STEP)
    return jacobian, "ok"


def solve_start(model, start_name, start_x, pool=None):
    y = np.log(start_x)
    history = []
    ill_conditioned_streak = 0
    accepted_iterations = 0
    total_candidates = 0
    while True:
        current = state(model, y)
        if not current["valid"]:
            return {
                "status": "invalid",
                "reason": current["reason"],
                "accepted_iterations": accepted_iterations,
                "total_candidates": total_candidates,
                "history": history,
            }
        history.append({
            "iteration": accepted_iterations,
            "residual_inf": current["residual_inf"],
            "residual_2": current["residual_2"],
            "x_min": float(np.min(current["x"])),
            "x_max": float(np.max(current["x"])),
            "minimum_gram_eigenvalue": current["minimum_absolute_gram_eigenvalue"],
            "minimum_angle_argument": current["minimum_angle_argument"],
        })
        if current["residual_inf"] < 1e-10 and current["residual_2"] < 3e-10:
            return {
                "status": "converged",
                "reason": "frozen residual tolerance",
                "accepted_iterations": accepted_iterations,
                "total_candidates": total_candidates,
                "y": y,
                "x": current["x"],
                "state": current,
                "history": history,
            }
        if accepted_iterations >= MAX_ACCEPTED_ITERATIONS:
            return {
                "status": "iteration_limit",
                "reason": "80 accepted iterations",
                "accepted_iterations": accepted_iterations,
                "total_candidates": total_candidates,
                "y": y,
                "x": current["x"],
                "state": current,
                "history": history,
            }
        jacobian, jacobian_reason = jacobian_y(model, y, pool=pool)
        if jacobian is None:
            status = (
                "boundary_contact"
                if jacobian_reason == "artificial box boundary"
                else "jacobian_failure"
            )
            return {
                "status": status,
                "reason": jacobian_reason,
                "accepted_iterations": accepted_iterations,
                "total_candidates": total_candidates,
                "y": y,
                "x": current["x"],
                "state": current,
                "history": history,
            }
        singular_values = np.linalg.svd(jacobian, compute_uv=False)
        condition = (
            math.inf if singular_values[-1] == 0
            else float(singular_values[0]/singular_values[-1])
        )
        history[-1]["jacobian_condition"] = condition
        if condition > 1e14:
            ill_conditioned_streak += 1
        else:
            ill_conditioned_streak = 0
        if ill_conditioned_streak >= 3:
            return {
                "status": "ill_conditioned",
                "reason": "condition >1e14 for three iterations",
                "accepted_iterations": accepted_iterations,
                "total_candidates": total_candidates,
                "y": y,
                "x": current["x"],
                "state": current,
                "history": history,
            }
        normal = jacobian.T @ jacobian
        right_hand_side = -jacobian.T @ current["residual"]
        sigma_scale = max(1.0, singular_values[0]**2)
        accepted = None
        for relative_mu in DAMPING:
            try:
                step = np.linalg.solve(
                    normal + relative_mu*sigma_scale*np.eye(35),
                    right_hand_side,
                )
            except np.linalg.LinAlgError:
                continue
            for backtrack in range(21):
                alpha = 2.0**(-backtrack)
                candidate_y = y + alpha*step
                total_candidates += 1
                candidate = state(model, candidate_y)
                if not candidate["valid"]:
                    continue
                if candidate["residual_2"] <= (
                    1-1e-4*alpha
                )*current["residual_2"]:
                    accepted = (
                        candidate_y, relative_mu, backtrack,
                        candidate["residual_2"],
                    )
                    break
            if accepted is not None:
                break
        if accepted is None:
            return {
                "status": "no_step",
                "reason": "no causal Armijo step in frozen list",
                "accepted_iterations": accepted_iterations,
                "total_candidates": total_candidates,
                "y": y,
                "x": current["x"],
                "state": current,
                "history": history,
            }
        y, relative_mu, backtrack, new_norm = accepted
        history[-1]["accepted_relative_mu"] = relative_mu
        history[-1]["accepted_backtrack"] = backtrack
        history[-1]["accepted_new_residual_2"] = new_norm
        accepted_iterations += 1


def serializable_search(result):
    payload = {
        key: value for key, value in result.items()
        if key not in {"y", "x", "state"}
    }
    if "x" in result:
        payload["x"] = result["x"].tolist()
    if "state" in result:
        state_result = result["state"]
        payload["terminal"] = {
            "residual_inf": state_result["residual_inf"],
            "residual_2": state_result["residual_2"],
            "action_real": float(state_result["action"].real),
            "action_imaginary": float(state_result["action"].imag),
            "minimum_absolute_gram_eigenvalue": state_result["minimum_absolute_gram_eigenvalue"],
            "minimum_angle_argument": state_result["minimum_angle_argument"],
        }
    return payload


print("=" * 78)
print("GLOBAL 35-VARIABLE LORENTZIAN REGGE ROOT SEARCH")
print("=" * 78)

searches = {parity: {} for parity in ("even", "odd")}
fork_context = mp.get_context("fork")
for parity in ("even", "odd"):
    model = gro.models[parity]
    with fork_context.Pool(
        processes=8, initializer=initialize_worker, initargs=(model,)
    ) as pool:
        for start_name, start_x in STARTS.items():
            initial = state(model, np.log(start_x))
            check(
                f"{parity}/{start_name}: frozen start satisfies every causal gate",
                initial["valid"],
                initial.get("reason", ""),
            )
            if not initial["valid"]:
                searches[parity][start_name] = {
                    "status": "invalid_start", "reason": initial["reason"]
                }
                continue
            print(
                f"[RUN ] {parity}/{start_name}: "
                f"||r||2={initial['residual_2']:.9g}, "
                f"||r||inf={initial['residual_inf']:.9g}",
                flush=True,
            )
            result = solve_start(model, start_name, start_x, pool=pool)
            searches[parity][start_name] = result
            terminal = result.get("state", {})
            detail = (
                f"status={result['status']}, accepted={result['accepted_iterations']}, "
                f"candidates={result['total_candidates']}"
            )
            if terminal:
                detail += (
                    f", ||r||2={terminal['residual_2']:.3e}, "
                    f"||r||inf={terminal['residual_inf']:.3e}"
                )
            print(f"[DONE] {parity}/{start_name}: {detail}", flush=True)

# Cluster all converged candidates before the expensive post-root controls.
converged = [
    (parity, start_name, result)
    for parity, parity_results in searches.items()
    for start_name, result in parity_results.items()
    if result["status"] == "converged"
]


clusters = []
for parity, start_name, result in converged:
    assigned = False
    for cluster in clusters:
        if cluster["parity"] != parity:
            continue
        representative = cluster["representative"]["x"]
        distance = float(np.linalg.norm(result["x"]-representative))/max(
            1.0, float(np.linalg.norm(result["x"])), float(np.linalg.norm(representative))
        )
        if distance < 1e-8:
            cluster["members"].append((start_name, result))
            assigned = True
            break
    if not assigned:
        clusters.append({
            "parity": parity,
            "representative": result,
            "members": [(start_name, result)],
        })


def hessian_x_at_relative_step(model, x, relative_step):
    """Differentiate orbit gradients and return the 35x35 x-Hessian."""
    hessian = np.empty((35, 35))
    for variable in range(35):
        delta = relative_step*x[variable]
        plus_x = x.copy()
        minus_x = x.copy()
        plus_x[variable] += delta
        minus_x[variable] -= delta
        plus_state = state(model, np.log(plus_x))
        minus_state = state(model, np.log(minus_x))
        if not plus_state["valid"] or not minus_state["valid"]:
            return None, "noncausal Hessian control"
        hessian[:, variable] = (
            plus_state["orbit_gradient"].real-minus_state["orbit_gradient"].real
        )/(2*delta)
    return hessian, "ok"


validation_records = []
validated_clusters = []
for cluster_index, cluster in enumerate(clusters):
    parity = cluster["parity"]
    result = cluster["representative"]
    start_names = [name for name, _ in cluster["members"]]
    model = gro.models[parity]
    y = result["y"]
    x = result["x"]
    reduced = state(model, y, need_full=True)
    full_data = reduced["full_data"]
    individual_gradient = full_data["individual_gradient"]
    max_individual = float(np.max(np.abs(individual_gradient)))
    full_action_error = gro.relative_error(reduced["action"], reduced["full_action"])
    full_restricted_error = max(
        gro.relative_error(left, right)
        for left, right in zip(
            reduced["orbit_gradient"], reduced["full_restricted_gradient"]
        )
    )
    check(
        f"{parity}/root{cluster_index}: candidate passes unreduced 840-gradient gate",
        max_individual < 2e-10
        and full_action_error < 2e-10
        and full_restricted_error < 2e-10,
        f"max individual={max_individual:.3e}, action={full_action_error:.3e}, "
        f"restricted={full_restricted_error:.3e}",
    )
    validation_ok = (
        max_individual < 2e-10
        and full_action_error < 2e-10
        and full_restricted_error < 2e-10
    )

    hessian_steps = (1e-4, 3e-5, 1e-5)
    hessian_records = {}
    hessian_control_ok = True
    for relative_step in hessian_steps:
        hessian, reason = hessian_x_at_relative_step(model, x, relative_step)
        if hessian is None:
            hessian_control_ok = False
            hessian_records[str(relative_step)] = {"valid": False, "reason": reason}
            continue
        symmetry = float(
            np.linalg.norm(hessian-hessian.T)/max(1.0, np.linalg.norm(hessian))
        )
        singular_values = np.linalg.svd((hessian+hessian.T)/2, compute_uv=False)
        ranks = {
            str(threshold): int(np.sum(
                singular_values > threshold*singular_values[0]
            ))
            for threshold in (1e-7, 1e-9, 1e-11)
        }
        hessian_records[str(relative_step)] = {
            "valid": True,
            "symmetry_residual": symmetry,
            "singular_values": singular_values.tolist(),
            "ranks": ranks,
        }
        hessian_control_ok &= symmetry < 3e-5
    rank_values = [
        tuple(record["ranks"].values())
        for record in hessian_records.values() if record.get("valid")
    ]
    stable_ranks = (
        len(rank_values) == 3
        and len({value for row in rank_values for value in row}) == 1
    )
    hessian_control_ok &= stable_ranks
    check(
        f"{parity}/root{cluster_index}: Hessian rank is stable at three steps and thresholds",
        hessian_control_ok,
        f"ranks={rank_values}",
    )
    validation_ok &= hessian_control_ok

    action_difference_records = {}
    action_difference_ok = True
    for relative_step in (1e-4, 3e-5):
        direct_gradient = np.empty(35, dtype=complex)
        causal_directions = True
        for variable in range(35):
            delta = relative_step*x[variable]
            plus_x = x.copy()
            minus_x = x.copy()
            plus_x[variable] += delta
            minus_x[variable] -= delta
            plus_state = state(model, np.log(plus_x))
            minus_state = state(model, np.log(minus_x))
            if not plus_state["valid"] or not minus_state["valid"]:
                causal_directions = False
                break
            plus_action = gro.full_evaluation(model, plus_x)[0]
            minus_action = gro.full_evaluation(model, minus_x)[0]
            direct_gradient[variable] = (plus_action-minus_action)/(2*delta)
        if causal_directions:
            errors = [
                gro.relative_error(left, right)
                for left, right in zip(direct_gradient, reduced["orbit_gradient"])
            ]
            maximum_error = max(errors)
            imaginary_residual = float(np.max(np.abs(direct_gradient.imag)))
        else:
            maximum_error = math.inf
            imaginary_residual = math.inf
        action_difference_records[str(relative_step)] = {
            "causal": causal_directions,
            "maximum_relative_error": maximum_error,
            "imaginary_residual": imaginary_residual,
        }
        action_difference_ok &= (
            causal_directions and maximum_error < 3e-5
            and imaginary_residual < 3e-5
        )
    check(
        f"{parity}/root{cluster_index}: complete-action differences validate all 35 derivatives",
        action_difference_ok,
        f"errors={[record['maximum_relative_error'] for record in action_difference_records.values()]}",
    )
    validation_ok &= action_difference_ok

    rounded_x = np.round(x, 14)
    rounded_restart = solve_start(model, "rounded", rounded_x)
    if rounded_restart["status"] == "converged":
        restart_distance = float(np.linalg.norm(rounded_restart["x"]-x))/max(
            1.0, float(np.linalg.norm(x))
        )
    else:
        restart_distance = math.inf
    restart_ok = (
        rounded_restart["status"] == "converged" and restart_distance < 2e-10
    )
    check(
        f"{parity}/root{cluster_index}: 14-digit rounded restart returns to the same root",
        restart_ok,
        f"status={rounded_restart['status']}, distance={restart_distance:.3e}",
    )
    validation_ok &= restart_ok

    record = {
        "parity": parity,
        "cluster": cluster_index,
        "starts": start_names,
        "x": x.tolist(),
        "max_individual_gradient": max_individual,
        "full_action_relative_error": full_action_error,
        "full_restricted_gradient_relative_error": full_restricted_error,
        "root_action_real": float(reduced["action"].real),
        "hessian_controls": hessian_records,
        "complete_action_difference_controls": action_difference_records,
        "rounded_restart_status": rounded_restart["status"],
        "rounded_restart_relative_distance": restart_distance,
        "validated": bool(validation_ok),
    }
    validation_records.append(record)
    if validation_ok:
        validated_clusters.append(record)


status_counts = {
    parity: dict(Counter(result["status"] for result in results.values()))
    for parity, results in searches.items()
}
listed_statuses = {
    "converged", "boundary_contact", "jacobian_failure", "iteration_limit",
    "no_step", "ill_conditioned", "invalid", "invalid_start",
}
all_searches = [result for results in searches.values() for result in results.values()]
check(
    "all twelve frozen searches terminate in a protocol-listed state",
    len(all_searches) == 12
    and all(result["status"] in listed_statuses for result in all_searches),
    f"status counts={status_counts}",
)

first_run_terminal_r2 = {
    "even": {
        "S0": 2.5059651478046785, "S1": 2.7212104112995137,
        "S2": 2.5057927432693554, "S3": 2.7160164213061537,
        "S4": 2.573114204835357, "S5": 2.603435539860999,
    },
    "odd": {
        "S0": 2.6832154785213906, "S1": 2.8502777897902614,
        "S2": 2.5895949945436465, "S3": 2.8578022096555213,
        "S4": 2.6056393737342565, "S5": 2.5924775428825457,
    },
}
reproduction_errors = []
for parity, parity_results in searches.items():
    for start_name, result in parity_results.items():
        terminal_r2 = result["history"][-1]["residual_2"]
        expected = first_run_terminal_r2[parity][start_name]
        reproduction_errors.append(abs(terminal_r2-expected)/max(1.0, expected))
check(
    "parallel rerun reproduces every disclosed terminal residual",
    max(reproduction_errors) < 2e-8
    and status_counts == {
        "even": {
            "boundary_contact": 3,
            "iteration_limit": 2,
            "jacobian_failure": 1,
        },
        "odd": {"boundary_contact": 4, "iteration_limit": 2},
    },
    f"max relative error={max(reproduction_errors):.3e}",
)
verdict = (
    "DERIVED NUMERICAL WITNESS: at least one stationary root in a frozen "
    "35-variable invariant subspace passed the complete float64 validation; "
    "an arbitrary-precision certificate is still required."
    if validated_clusters else
    "DERIVED SCOPED NEGATIVE: no stationary point was found from the six frozen "
    "starts in the stated causal box and 35-variable invariant subspaces, or "
    "no converged candidate passed validation.  This is not a global "
    "nonexistence result."
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "correction_commit": CORRECTION_COMMIT,
    "classification_correction_commit": CLASSIFICATION_CORRECTION_COMMIT,
    "evaluator_commit": EVALUATOR_COMMIT,
    "disclosed_first_run": {
        "passed": 14,
        "tests": 14,
        "validated_roots": 0,
        "generic_jacobian_failures_reclassified_as_boundary_contacts": 7,
        "causal_branch_jacobian_failures": 1,
        "iteration_limits": 4,
        "terminal_residual_2_min": min(
            value for row in first_run_terminal_r2.values() for value in row.values()
        ),
        "terminal_residual_2_max": max(
            value for row in first_run_terminal_r2.values() for value in row.values()
        ),
    },
    "starts_per_parity": len(STARTS),
    "searches": {
        parity: {
            start_name: serializable_search(result)
            for start_name, result in parity_results.items()
        }
        for parity, parity_results in searches.items()
    },
    "status_counts": status_counts,
    "converged_candidates": len(converged),
    "distinct_converged_clusters": len(clusters),
    "validated_clusters": len(validated_clusters),
    "validations": validation_records,
    "verdict": verdict,
    "labels": {
        "roots_found": "DERIVED NUMERICAL WITNESS" if validated_clusters else "NONE VALIDATED FROM FROZEN STARTS",
        "no_root_in_full_35_variable_box": "NOT PROVED",
        "no_root_in_full_840_variable_space": "NOT PROVED",
        "phase_parity_selection": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)
