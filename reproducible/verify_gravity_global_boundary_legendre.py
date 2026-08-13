#!/usr/bin/env python3
"""Global 65-variable boundary-Legendre and rectangular-rank audit.

Protocol commit: 8c2482b.  Final-boundary lengths are variables; their
derivatives are post-momenta, while only the 35 internal derivatives are
constraints.  This verifier performs no stationary-root search.
"""

from collections import Counter, defaultdict
import contextlib
import importlib.util
import io
from itertools import combinations
import json
import math
import multiprocessing as mp
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_global_boundary_legendre.json"
PROTOCOL_COMMIT = "8c2482b"
UPSTREAM_COMMIT = "d9fe159"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def relative_error(left, right):
    return abs(left-right)/max(1.0, abs(left), abs(right))


# Load the already certified action branch and slab reconstruction.
spec = importlib.util.spec_from_file_location(
    "global_regge_orbits", HERE / "verify_gravity_global_regge_orbits.py"
)
gro = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(gro)


check(
    "the imported action branch retains all 43 upstream certificates",
    gro.tests == gro.passed == 43,
)


def orbit_sort_key(orbit, phase):
    representative = min(orbit)
    logical = tuple(vertex % 120 for vertex in representative)
    phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
    return phase_pair, tuple(sorted(orbit))


def augment_model(base):
    old_orbits = tuple(sorted(
        gro.orbit_partition(base["old_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    final_orbits = tuple(sorted(
        gro.orbit_partition(base["new_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    old_to_orbit = {
        edge: index for index, orbit in enumerate(old_orbits) for edge in orbit
    }
    final_to_variable = {
        edge: 35+index
        for index, orbit in enumerate(final_orbits) for edge in orbit
    }
    all_edge_to_variable = dict(base["edge_to_variable"])
    all_edge_to_variable.update(final_to_variable)
    all_edge_jacobian = dict(base["edge_jacobian"])
    all_edge_jacobian.update({edge: 1.0 for edge in base["new_edges"]})
    all_edges = tuple(sorted(
        base["internal_edges"] | base["old_edges"] | base["new_edges"]
    ))
    all_edge_position = {edge: index for index, edge in enumerate(all_edges)}
    return {
        **base,
        "old_orbits": old_orbits,
        "final_orbits": final_orbits,
        "old_to_orbit": old_to_orbit,
        "final_to_variable": final_to_variable,
        "all_edge_to_variable": all_edge_to_variable,
        "all_edge_jacobian": all_edge_jacobian,
        "all_edges": all_edges,
        "all_edge_position": all_edge_position,
    }


models = {parity: augment_model(base) for parity, base in gro.models.items()}

for parity, model in models.items():
    final_phase_pairs = Counter()
    for orbit in model["final_orbits"]:
        representative = min(orbit)
        logical = tuple(vertex-120 for vertex in representative)
        final_phase_pairs[tuple(sorted(
            model["phase"][vertex] for vertex in logical
        ))] += 1
    boundary_triangle_orbits = tuple(
        orbit for orbit in model["triangle_orbits"]
        if min(orbit) in model["boundary_triangles"]
    )
    check(
        f"{parity}: old and final boundaries each have 30 edge orbits of 24",
        len(model["old_orbits"]) == len(model["final_orbits"]) == 30
        and Counter(map(len, model["old_orbits"])) == Counter({24: 30})
        and Counter(map(len, model["final_orbits"])) == Counter({24: 30}),
    )
    check(
        f"{parity}: every phase pair has three final-boundary edge orbits",
        final_phase_pairs == Counter({pair: 3 for pair in combinations(range(5), 2)}),
    )
    check(
        f"{parity}: 2400 boundary triangles form 100 orbits of 24",
        len(boundary_triangle_orbits) == 100
        and Counter(map(len, boundary_triangle_orbits)) == Counter({24: 100}),
    )


def edge_square(model, edge, variables, old_values=None):
    edge = tuple(sorted(edge))
    if edge in model["all_edge_to_variable"]:
        variable = model["all_edge_to_variable"][edge]
        return model["all_edge_jacobian"][edge]*variables[variable]
    if edge in model["old_to_orbit"]:
        if old_values is None:
            return 1.0
        return old_values[model["old_to_orbit"][edge]]
    raise ValueError(f"edge absent from augmented slab: {edge}")


def simplex_squared(model, simplex, variables, old_values=None):
    squared = np.zeros((5, 5))
    for left, right in combinations(range(5), 2):
        value = edge_square(
            model, (simplex[left], simplex[right]), variables, old_values
        )
        squared[left, right] = squared[right, left] = value
    return squared


def triangle_area_data(model, triangle, variables, old_values=None):
    edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
    values = tuple(edge_square(model, edge, variables, old_values) for edge in edges)
    area_square = gro.triangle_area_square(values)
    variable_derivatives = {}
    old_derivatives = {}
    for edge, partial in zip(edges, gro.triangle_area_partials(values)):
        if edge in model["all_edge_to_variable"]:
            variable = model["all_edge_to_variable"][edge]
            variable_derivatives[variable] = (
                partial*model["all_edge_jacobian"][edge]
            )
        elif edge in model["old_to_orbit"]:
            old_derivatives[model["old_to_orbit"][edge]] = partial
    return area_square, variable_derivatives, old_derivatives


def reduced_evaluation(model, variables, old_values=None):
    triangle_orbits = model["triangle_orbits"]
    curvature = np.array([
        gro.hinge_constant(model, min(orbit)) for orbit in triangle_orbits
    ], dtype=complex)
    minimum_argument = math.inf
    minimum_gram = math.inf
    negative_counts = Counter()
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        squared = simplex_squared(model, simplex, variables, old_values)
        angles, arguments, eigenvalues = gro.angle_data(squared)
        negative_counts[int(np.sum(eigenvalues < -1e-10))] += 1
        minimum_gram = min(minimum_gram, float(np.min(np.abs(eigenvalues))))
        for local_hinge, angle in angles.items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            orbit_index = model["triangle_to_orbit"][triangle]
            curvature[orbit_index] += angle
            minimum_argument = min(minimum_argument, abs(arguments[local_hinge]))
    action_sum = 0j
    gradient = np.zeros(65, dtype=complex)
    old_gradient = np.zeros(30, dtype=complex)
    for orbit_index, orbit in enumerate(triangle_orbits):
        triangle = min(orbit)
        area_square, variable_derivatives, old_derivatives = triangle_area_data(
            model, triangle, variables, old_values
        )
        root_area = np.lib.scimath.sqrt(area_square)
        triangle_curvature = curvature[orbit_index]
        action_sum += 24*root_area*triangle_curvature
        for variable, derivative in variable_derivatives.items():
            gradient[variable] += -1j*24*triangle_curvature*derivative/(2*root_area)
        for variable, derivative in old_derivatives.items():
            old_gradient[variable] += -1j*24*triangle_curvature*derivative/(2*root_area)
    return -1j*action_sum, gradient, old_gradient, {
        "curvature": curvature,
        "minimum_argument": minimum_argument,
        "minimum_gram": minimum_gram,
        "negative_counts": negative_counts,
    }


def full_evaluation(model, variables, old_values=None):
    angle_incidence = defaultdict(list)
    minimum_argument = math.inf
    minimum_gram = math.inf
    negative_counts = Counter()
    for simplex in model["slab"]:
        squared = simplex_squared(model, simplex, variables, old_values)
        angles, arguments, eigenvalues = gro.angle_data(squared)
        negative_counts[int(np.sum(eigenvalues < -1e-10))] += 1
        minimum_gram = min(minimum_gram, float(np.min(np.abs(eigenvalues))))
        for local_hinge, angle in angles.items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            angle_incidence[triangle].append(angle)
            minimum_argument = min(minimum_argument, abs(arguments[local_hinge]))
    curvature = {
        triangle: gro.hinge_constant(model, triangle)+sum(angles)
        for triangle, angles in angle_incidence.items()
    }
    action_sum = 0j
    edge_gradient = np.zeros(len(model["all_edges"]), dtype=complex)
    for triangle, triangle_curvature in curvature.items():
        edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
        values = tuple(edge_square(model, edge, variables, old_values) for edge in edges)
        area_square = gro.triangle_area_square(values)
        root_area = np.lib.scimath.sqrt(area_square)
        action_sum += root_area*triangle_curvature
        for edge, partial in zip(edges, gro.triangle_area_partials(values)):
            jacobian = -1.0 if edge in model["internal_edges"] and edge[1]-edge[0] == 120 else 1.0
            edge_gradient[model["all_edge_position"][edge]] += (
                -1j*triangle_curvature*partial*jacobian/(2*root_area)
            )
    gradient = np.zeros(65, dtype=complex)
    for variable, orbit in enumerate(model["edge_orbits"]):
        gradient[variable] = sum(
            edge_gradient[model["all_edge_position"][edge]] for edge in orbit
        )
    for final_index, orbit in enumerate(model["final_orbits"]):
        gradient[35+final_index] = sum(
            edge_gradient[model["all_edge_position"][edge]] for edge in orbit
        )
    old_gradient = np.array([
        sum(edge_gradient[model["all_edge_position"][edge]] for edge in orbit)
        for orbit in model["old_orbits"]
    ])
    orbit_curvature = np.empty(len(model["triangle_orbits"]), dtype=complex)
    maximum_curvature_spread = 0.0
    for orbit_index, orbit in enumerate(model["triangle_orbits"]):
        values = [curvature[triangle] for triangle in orbit]
        orbit_curvature[orbit_index] = sum(values)/len(values)
        maximum_curvature_spread = max(
            maximum_curvature_spread,
            max(abs(value-orbit_curvature[orbit_index]) for value in values),
        )
    return -1j*action_sum, gradient, old_gradient, {
        "curvature": orbit_curvature,
        "minimum_argument": minimum_argument,
        "minimum_gram": minimum_gram,
        "negative_counts": negative_counts,
        "maximum_curvature_spread": maximum_curvature_spread,
        "edge_gradient": edge_gradient,
    }


def controls():
    j = np.arange(1, 31, dtype=float)
    k = np.arange(1, 6, dtype=float)
    l = np.arange(1, 31, dtype=float)
    return {
        "B0": np.r_[np.ones(30), np.full(5, 0.25), np.ones(30)],
        "B1": np.r_[1+j/1000, 0.25+k/1000, 1+l/1500],
        "B2": np.r_[1-j/2000, 0.25+(6-k)/1500, 1-l/2500],
    }


CONTROLS = controls()
control_records = {parity: {} for parity in models}

print("=" * 78)
print("GLOBAL REGGE BOUNDARY-LEGENDRE RANK")
print("=" * 78)

for parity, model in models.items():
    for name, variables in CONTROLS.items():
        reduced_action, reduced_gradient, reduced_old, reduced_data = reduced_evaluation(
            model, variables
        )
        full_action, full_gradient, full_old, full_data = full_evaluation(
            model, variables
        )
        action_error = relative_error(reduced_action, full_action)
        gradient_error = max(
            relative_error(left, right)
            for left, right in zip(reduced_gradient, full_gradient)
        )
        old_gradient_error = max(
            relative_error(left, right)
            for left, right in zip(reduced_old, full_old)
        )
        curvature_error = max(
            relative_error(left, right)
            for left, right in zip(reduced_data["curvature"], full_data["curvature"])
        )
        real_residual = max(
            abs(reduced_action.imag),
            float(np.max(np.abs(reduced_gradient.imag))),
            float(np.max(np.abs(reduced_old.imag))),
        )
        check(
            f"{parity}/{name}: 100-orbit action and all 95 derivatives match the full slab",
            action_error < 2e-8 and gradient_error < 2e-8
            and old_gradient_error < 2e-8 and curvature_error < 2e-9
            and full_data["maximum_curvature_spread"] < 2e-9,
            f"action={action_error:.2e}, grad65={gradient_error:.2e}, "
            f"old30={old_gradient_error:.2e}, curvature={curvature_error:.2e}",
        )
        check(
            f"{parity}/{name}: all controls remain on the real Lorentzian branch",
            reduced_data["negative_counts"] == Counter({1: 100})
            and full_data["negative_counts"] == Counter({1: 2400})
            and reduced_data["minimum_gram"] > 1e-8
            and reduced_data["minimum_argument"] > 1e-6
            and real_residual < 2e-7,
            f"real={real_residual:.2e}, gram={reduced_data['minimum_gram']:.3e}, "
            f"argument={reduced_data['minimum_argument']:.3e}",
        )
        control_records[parity][name] = {
            "action": reduced_action,
            "gradient": reduced_gradient,
            "old_gradient": reduced_old,
            "minimum_gram": reduced_data["minimum_gram"],
            "minimum_argument": reduced_data["minimum_argument"],
            "action_error": action_error,
            "gradient_error": gradient_error,
            "old_gradient_error": old_gradient_error,
            "curvature_error": curvature_error,
            "real_residual": real_residual,
        }


# B0 internal components must be identical to the upstream fixed-boundary result.
for parity, model in models.items():
    upstream_gradient = gro.control_results[parity]["R0"]["gradient"]
    current_gradient = control_records[parity]["B0"]["gradient"][:35]
    error = max(relative_error(a, b) for a, b in zip(upstream_gradient, current_gradient))
    check(
        f"{parity}: B0 internal gradient reproduces the certified 35-variable result",
        error < 2e-10,
        f"relative error={error:.3e}",
    )


# Direct full-action derivatives.  Fork workers share the read-only models.
_WORKER_MODEL = None


def initialize_worker(model):
    global _WORKER_MODEL
    _WORKER_MODEL = model


def full_action_worker(variables):
    return full_evaluation(_WORKER_MODEL, variables)[0]


direct_gradient_records = {parity: {} for parity in models}
fork_context = mp.get_context("fork")
direct_step = 3e-5
for parity, model in models.items():
    with fork_context.Pool(
        processes=8, initializer=initialize_worker, initargs=(model,)
    ) as pool:
        for name in ("B0", "B1"):
            variables = CONTROLS[name]
            points = []
            deltas = []
            for variable in range(65):
                delta = direct_step*variables[variable]
                plus = variables.copy()
                minus = variables.copy()
                plus[variable] += delta
                minus[variable] -= delta
                points.extend((plus, minus))
                deltas.append(delta)
            actions = pool.map(full_action_worker, points, chunksize=2)
            direct_gradient = np.array([
                (actions[2*variable]-actions[2*variable+1])/(2*deltas[variable])
                for variable in range(65)
            ])
            analytic = control_records[parity][name]["gradient"]
            maximum_error = max(
                relative_error(left, right)
                for left, right in zip(direct_gradient, analytic)
            )
            imaginary_residual = float(np.max(np.abs(direct_gradient.imag)))
            direct_gradient_records[parity][name] = {
                "gradient": direct_gradient,
                "maximum_error": maximum_error,
                "imaginary_residual": imaginary_residual,
            }
            check(
                f"{parity}/{name}: complete-action differences reproduce all 65 derivatives",
                maximum_error < 3e-5 and imaginary_residual < 3e-5,
                f"relative error={maximum_error:.3e}, imag={imaginary_residual:.3e}",
            )


# Rectangular internal-residual Jacobian at B0.
rank_records = {}
rank_step = 5e-4
for parity, model in models.items():
    variables = CONTROLS["B0"]
    jacobian = np.empty((35, 65))
    for variable in range(65):
        delta_y = rank_step
        plus = variables.copy()
        minus = variables.copy()
        plus[variable] *= math.exp(delta_y)
        minus[variable] *= math.exp(-delta_y)
        plus_gradient = reduced_evaluation(model, plus)[1][:35].real/24
        minus_gradient = reduced_evaluation(model, minus)[1][:35].real/24
        jacobian[:, variable] = (plus_gradient-minus_gradient)/(2*delta_y)
    internal_block = jacobian[:, :35]
    final_block = jacobian[:, 35:]
    spectra = {
        "combined": np.linalg.svd(jacobian, compute_uv=False),
        "internal": np.linalg.svd(internal_block, compute_uv=False),
        "final": np.linalg.svd(final_block, compute_uv=False),
    }
    ranks = {
        block: {
            str(threshold): int(np.sum(values > threshold*values[0]))
            for threshold in (1e-7, 1e-9, 1e-11)
        }
        for block, values in spectra.items()
    }
    residual = control_records[parity]["B0"]["gradient"][:35].real/24
    u_final, s_final, _ = np.linalg.svd(final_block, full_matrices=True)
    final_rank = ranks["final"][str(1e-9)]
    projection = u_final[:, :final_rank] @ (u_final[:, :final_rank].T @ residual)
    orthogonal = residual-projection
    combined_step = np.linalg.lstsq(jacobian, -residual, rcond=None)[0]
    final_step, _, _, _ = np.linalg.lstsq(final_block, -residual, rcond=None)
    final_linear_residual = residual+final_block@final_step
    rank_records[parity] = {
        "jacobian": jacobian,
        "spectra": spectra,
        "ranks": ranks,
        "internal_condition": float(spectra["internal"][0]/spectra["internal"][-1]),
        "residual_norm": float(np.linalg.norm(residual)),
        "projection_norm": float(np.linalg.norm(projection)),
        "orthogonal_norm": float(np.linalg.norm(orthogonal)),
        "combined_step_norm": float(np.linalg.norm(combined_step)),
        "combined_linear_residual": float(np.linalg.norm(residual+jacobian@combined_step)),
        "final_step_norm": float(np.linalg.norm(final_step)),
        "final_linear_residual": float(np.linalg.norm(final_linear_residual)),
    }
    check(
        f"{parity}: rectangular/internal/final ranks are stable at all frozen thresholds",
        all(len(set(block.values())) == 1 for block in ranks.values()),
        f"ranks={ranks}",
    )
    check(
        f"{parity}: the combined linearized system solves every internal residual direction",
        ranks["combined"] == {"1e-07": 35, "1e-09": 35, "1e-11": 35}
        and rank_records[parity]["combined_linear_residual"] < 1e-10,
        f"linear residual={rank_records[parity]['combined_linear_residual']:.3e}",
    )


# At the regular control, compare boundary derivative multisets and direct values.
momentum_records = {}
for parity in models:
    old_derivative = control_records[parity]["B0"]["old_gradient"].real
    final_derivative = control_records[parity]["B0"]["gradient"][35:].real
    old_direct = np.empty(30)
    # Direct old-boundary differences are cheap in the orbit evaluator and
    # independent of the analytic area derivative.
    model = models[parity]
    base_variables = CONTROLS["B0"]
    old_base = np.ones(30)
    for variable in range(30):
        delta = direct_step
        plus_old = old_base.copy()
        minus_old = old_base.copy()
        plus_old[variable] += delta
        minus_old[variable] -= delta
        plus_action = reduced_evaluation(model, base_variables, plus_old)[0]
        minus_action = reduced_evaluation(model, base_variables, minus_old)[0]
        old_direct[variable] = ((plus_action-minus_action)/(2*delta)).real
    old_direct_error = max(
        relative_error(left, right) for left, right in zip(old_direct, old_derivative)
    )
    final_direct = direct_gradient_records[parity]["B0"]["gradient"][35:].real
    final_direct_error = max(
        relative_error(left, right) for left, right in zip(final_direct, final_derivative)
    )
    multiset_error = max(
        abs(left-right) for left, right in zip(
            np.sort(old_derivative), np.sort(final_derivative)
        )
    )/max(1.0, float(np.max(np.abs(old_derivative))), float(np.max(np.abs(final_derivative))))
    # Canonical signs: p_pre=-dS/d(old), p_post=+dS/d(final).  Equal derivative
    # multisets therefore mean p_pre and p_post are opposite under reversal.
    check(
        f"{parity}: old and final derivative multisets obey regular time reversal",
        multiset_error < 2e-10,
        f"relative multiset error={multiset_error:.3e}",
    )
    check(
        f"{parity}: direct differences independently reproduce all boundary momenta",
        old_direct_error < 3e-5 and final_direct_error < 3e-5
        and np.max(np.abs(old_derivative.imag if np.iscomplexobj(old_derivative) else 0)) < 2e-10,
        f"old={old_direct_error:.3e}, final={final_direct_error:.3e}",
    )
    momentum_records[parity] = {
        "old_action_derivatives": old_derivative,
        "final_action_derivatives": final_derivative,
        "pre_momenta": -old_derivative,
        "post_momenta": final_derivative,
        "multiset_error": multiset_error,
        "old_direct_error": old_direct_error,
        "final_direct_error": final_direct_error,
    }


# Parity comparison is reported, not assumed equal.
parity_rank_spectrum_errors = {}
for block in ("combined", "internal", "final"):
    even = rank_records["even"]["spectra"][block]
    odd = rank_records["odd"]["spectra"][block]
    parity_rank_spectrum_errors[block] = float(
        np.max(np.abs(even-odd))/max(1.0, even[0], odd[0])
    )

verdict = (
    "DERIVED COMPUTATIONAL: the global slab defines a real 65-variable "
    "boundary-Legendre action with 35 internal constraints and 30 final "
    "post-momenta.  Its rectangular rank and boundary-response subspace are "
    "reported without a root search."
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "upstream_commit": UPSTREAM_COMMIT,
    "variables": {
        "internal": 35,
        "final_boundary": 30,
        "total": 65,
        "old_boundary_momentum_orbits": 30,
    },
    "controls": {
        parity: {
            name: {
                "action_real": float(record["action"].real),
                "action_imaginary": float(record["action"].imag),
                "internal_gradient_norm": float(np.linalg.norm(record["gradient"][:35].real)),
                "post_momentum_norm": float(np.linalg.norm(record["gradient"][35:].real)),
                "pre_derivative_norm": float(np.linalg.norm(record["old_gradient"].real)),
                "minimum_gram": record["minimum_gram"],
                "minimum_angle_argument": record["minimum_argument"],
                "full_action_error": record["action_error"],
                "full_gradient_error": record["gradient_error"],
                "full_old_gradient_error": record["old_gradient_error"],
            }
            for name, record in parity_records.items()
        }
        for parity, parity_records in control_records.items()
    },
    "rank": {
        parity: {
            "ranks": record["ranks"],
            "singular_values": {
                block: values.tolist() for block, values in record["spectra"].items()
            },
            "internal_condition": record["internal_condition"],
            "regular_residual_norm": record["residual_norm"],
            "final_image_projection_norm": record["projection_norm"],
            "final_left_null_projection_norm": record["orthogonal_norm"],
            "combined_minimum_norm_step": record["combined_step_norm"],
            "combined_linear_residual": record["combined_linear_residual"],
            "final_only_minimum_norm_step": record["final_step_norm"],
            "final_only_linear_residual": record["final_linear_residual"],
        }
        for parity, record in rank_records.items()
    },
    "boundary_momenta": {
        parity: {
            key: value.tolist() if isinstance(value, np.ndarray) else value
            for key, value in record.items()
        }
        for parity, record in momentum_records.items()
    },
    "phase_parity_rank_spectrum_errors": parity_rank_spectrum_errors,
    "labels": {
        "65_variable_boundary_action": "DERIVED COMPUTATIONAL",
        "boundary_momenta": "DERIVED COMPUTATIONAL OFF SHELL",
        "stationary_internal_solution": "NOT SEARCHED",
        "canonical_evolution": "OPEN UNTIL ON SHELL",
    },
    "verdict": verdict,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)
