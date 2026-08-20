#!/usr/bin/env python3
"""Bounded positive-root search for the refined H4 internal equations.

Prior-art commit: 29162db.
Protocol commit: 722fb3c.
"""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
FILL = HERE / "gravity_600cell_refined_h4_stationary_fill.json"
JACOBIAN_SOURCE = HERE / "verify_gravity_600cell_refined_h4_internal_jacobian.py"
JACOBIAN = HERE / "gravity_600cell_refined_h4_internal_jacobian.json"
JACOBIAN_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_internal_jacobian_result.md"
)
PRIOR_ART = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_stationary_root_prior_art.md"
)
PROTOCOL = (
    ROOT / "docs/gravity/gravity_600cell_refined_h4_stationary_root_protocol.md"
)
OUTPUT = HERE / "gravity_600cell_refined_h4_stationary_root.json"
PRIOR_ART_COMMIT = "29162db"
PROTOCOL_COMMIT = "722fb3c"
INPUT_HASHES = {
    "reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py":
        "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "reproducible/gravity_600cell_refined_h4_stationary_fill.json":
        "283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018",
    "reproducible/verify_gravity_600cell_refined_h4_internal_jacobian.py":
        "6f74f0a73d15b1e50e61e0afe56d74b162d4a98ac87368979f4bd52fe86b6b4e",
    "reproducible/gravity_600cell_refined_h4_internal_jacobian.json":
        "b900021c21df67c1de1ae18929be302b0d47d2f267c4a919388711a0a0bf5eaa",
    "docs/gravity/gravity_600cell_refined_h4_internal_jacobian_result.md":
        "8a6603c810a5615956c8cc2ec8a9b8c3a6e015a0f6495f2f8a58811585418ee5",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_prior_art.md":
        "450bb7ad0bee2107c7b80de652aa2dc04fd0e51e10f3798a7b717ec64e251478",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_protocol.md":
        "fccb4b6fd28b4e24dfdb99e657cc31f2a98ca0e5b0d500c099db5b3b3ec619d5",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_control_correction.md":
        "fc31a55898e2c7fb357e386ba9020e26e8f26096e80ddd982ba2d2d2c0ae3caf",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_lorentzian_correction.md":
        "230f7389c97f7356ca732240cf5a6ded82685c3dac44173883e6942be2143775",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_performance_interruption.md":
        "4a945ba4609a282bf0f65979a3f557fcca2214333868f913a41785cb6f258117",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_fast_evaluator_protocol.md":
        "f86d9bb13a88566930b78ef9e64b2bf026bfe7bd5d6a3cf62944ee50fdf1b9c3",
    "docs/gravity/gravity_600cell_refined_h4_stationary_root_determinism_correction.md":
        "f1a5422d4e036afe6c900ae7be8f1499367e1b2c68e85b23bd7d1c32774ef88d",
}
PAIR4 = tuple(combinations(range(4), 2))
VARIABLES = (
    tuple(("old",)+pair for pair in PAIR4)
    + tuple(("new",)+pair for pair in PAIR4)
    + tuple(("cross",)+pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
INTERNAL_VARIABLES = (
    tuple(("cross",)+pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
TAU_TEXT = "0.0102"
MAIN_LOWER = np.asarray([-0.35]*6+[-8.0]*4)
MAIN_UPPER = np.asarray([0.35]*6+[2.0]*4)
LADDER_LOWERS = (-4.0, -8.0, -12.0, -16.0)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def sf(value, digits=17):
    value = float(value)
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return format(value, f".{digits}g")


def mp_text(value, digits=60):
    return mp.nstr(value, digits)


def load_action_definitions():
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
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
        "LOCAL_TRIANGLES": np.asarray(
            tuple(combinations(range(5), 3)), dtype=np.int8
        ),
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
        compile(ast.Module(body=definitions, type_ignores=[]), str(SOURCE), "exec"),
        namespace,
    )
    return namespace


def parse_states(label):
    return tuple(
        (int(item.split("t")[0][1:]), int(item.split("t")[1]))
        for item in label.split("|")
    )


def parse_combinatorics(record):
    return {
        "order": tuple(record["order"]),
        "pentachora": record["pentachora"],
        "distinct_pentachora": record["distinct_pentachora"],
        "triangles": record["triangles"],
        "boundary_triangles": record["boundary_triangles"],
        "mixed_triangle_types": record["mixed_triangle_types"],
        "simplex_types": [
            {"states": parse_states(item["states"]), "count": item["count"]}
            for item in record["simplex_types"]
        ],
        "triangle_types": [
            {
                "states": parse_states(item["states"]),
                "count": item["count"],
                "incidence_values": item["incidence_values"],
                "contributions": [
                    {
                        "simplex": parse_states(value["simplex"]),
                        "triangle": parse_states(value["triangle"]),
                        "multiplicity": value["multiplicity"],
                    }
                    for value in item["signature"]
                ],
            }
            for item in record["triangle_types"]
        ],
    }


def coordinates_from_y(base, y):
    result = dict(base)
    for index, key in enumerate(INTERNAL_VARIABLES):
        result[key] *= mp.exp(mp.mpf(y[index]))
    return result


def fast_coordinates_from_y(base, y):
    result = dict(base)
    for index, key in enumerate(INTERNAL_VARIABLES):
        result[key] *= math.exp(float(y[index]))
    return result


def fast_edge_coordinate(left, right):
    rank_left, layer_left = left
    rank_right, layer_right = right
    if layer_left == layer_right:
        if rank_left == rank_right:
            raise ValueError("zero boundary edge")
        return (("old" if layer_left == 0 else "new",)
                + tuple(sorted((rank_left, rank_right)))), 1.0
    if rank_left == rank_right:
        return ("rho", rank_left), -1.0
    return ("cross",)+tuple(sorted((rank_left, rank_right))), 1.0


def fast_simplex_squared(states, coordinates):
    squared = np.zeros((5, 5), dtype=float)
    for left, right in combinations(range(5), 2):
        key, sign = fast_edge_coordinate(states[left], states[right])
        squared[left, right] = squared[right, left] = sign*coordinates[key]
    return squared


def fast_signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices)-1
    if dimension == 0:
        return 1.0
    base = vertices[0]
    others = vertices[1:]
    gram = np.asarray([
        [(squared[base, left]+squared[base, right]-squared[left, right])/2
         for right in others] for left in others
    ])
    return float(np.linalg.det(gram))/(math.factorial(dimension)**2)


def fast_log_minus(value):
    scale = max(1.0, abs(value))
    if abs(value.imag) < 1e-13*scale:
        real = value.real
        if real < 0:
            return math.log(-real)-1j*math.pi
        return complex(math.log(real), 0.0)
    return np.log(value)


def fast_angle_record(squared):
    gram = np.asarray([
        [(squared[0, left]+squared[0, right]-squared[left, right])/2
         for right in range(1, 5)] for left in range(1, 5)
    ])
    inverse = np.linalg.inv(gram)
    simplex_volume_square = fast_signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: fast_signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        ) for omitted in range(5)
    }
    angles = {}
    maximum_identity = 0.0
    minimum_argument = math.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge = tuple(vertex for vertex in range(5)
                      if vertex not in (omitted_a, omitted_b))
        hinge_volume_square = fast_signed_volume_square(squared, hinge)
        derivative = np.zeros((4, 4), dtype=float)
        opposite = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                derivative[left-1, right-1] = (
                    int({0, left} == opposite)
                    + int({0, right} == opposite)
                    - int(left != right and {left, right} == opposite)
                )/2
        volume_derivative = simplex_volume_square*np.trace(inverse@derivative)
        denominator = (
            np.sqrt(complex(facet_volume_squares[omitted_a]))
            * np.sqrt(complex(facet_volume_squares[omitted_b]))
        )
        cosine = 16*volume_derivative/denominator
        sine = (-4/3)*(
            np.sqrt(complex(hinge_volume_square))
            * np.sqrt(complex(simplex_volume_square))
        )/denominator
        maximum_identity = max(maximum_identity, abs(cosine*cosine+sine*sine-1))
        argument = cosine+1j*sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[hinge] = -1j*fast_log_minus(argument)
    return angles, maximum_identity, minimum_argument


def fast_triangle_area_and_derivatives(states, coordinates):
    edge_records = []
    for left, right in combinations(range(3), 2):
        key, sign = fast_edge_coordinate(states[left], states[right])
        edge_records.append((key, sign, sign*coordinates[key]))
    x, y, z = (record[2] for record in edge_records)
    area_square = (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16
    partials = ((y+z-x)/8, (x+z-y)/8, (x+y-z)/8)
    area = np.sqrt(complex(area_square))
    derivatives = [
        (key, partial*sign/(2*area))
        for (key, sign, _), partial in zip(edge_records, partials)
    ]
    return area, derivatives


def fast_evaluate_schedule(record, geometry, coordinates):
    angle_lookup = {}
    maximum_identity = 0.0
    minimum_argument = math.inf
    for simplex_record in record["simplex_types"]:
        states = simplex_record["states"]
        angles, identity, argument = fast_angle_record(
            fast_simplex_squared(states, coordinates)
        )
        maximum_identity = max(maximum_identity, identity)
        minimum_argument = min(minimum_argument, argument)
        for local_triangle, angle in angles.items():
            triangle_states = tuple(sorted(
                (states[index] for index in local_triangle),
                key=lambda state: state[0]+4*state[1],
            ))
            angle_lookup[(states, triangle_states)] = angle
    gravitational_sum = 0j
    gradient = {key: 0j for key in VARIABLES}
    curvatures = []
    for triangle_record in record["triangle_types"]:
        states = triangle_record["states"]
        boundary = len({layer for _, layer in states}) == 1
        curvature = math.pi if boundary else 2*math.pi
        for contribution in triangle_record["contributions"]:
            curvature += contribution["multiplicity"]*angle_lookup[
                (contribution["simplex"], contribution["triangle"])
            ]
        area, derivatives = fast_triangle_area_and_derivatives(states, coordinates)
        multiplicity = triangle_record["count"]
        gravitational_sum += multiplicity*area*curvature
        for key, area_derivative in derivatives:
            gradient[key] += (
                -1j*multiplicity*curvature*area_derivative*coordinates[key]
            )
        curvatures.append(curvature)
    gravitational = -1j*gravitational_sum
    dust = 0.0
    rank_mass = geometry["mass"]/4
    for rank in range(4):
        rho = coordinates["rho", rank]
        dust -= 8*math.pi*rank_mass*math.sqrt(rho)
        gradient["rho", rank] -= 4*math.pi*rank_mass*math.sqrt(rho)
    return {
        "action": gravitational+dust,
        "gravitational": gravitational,
        "dust": dust,
        "gradient": np.asarray([gradient[key] for key in INTERNAL_VARIABLES]),
        "maximum_angle_identity_residual": maximum_identity,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary_curvature": max(abs(value.imag) for value in curvatures),
    }


def fast_equation_state(record, geometry, base, y):
    try:
        result = fast_evaluate_schedule(
            record, geometry, fast_coordinates_from_y(base, y)
        )
        gradient = result["gradient"]
        maximum_imaginary = max(
            abs(result["action"].imag), abs(result["gravitational"].imag),
            float(np.max(np.abs(gradient.imag))),
        )
        finite = bool(
            np.all(np.isfinite(gradient))
            and np.isfinite(result["action"])
            and math.isfinite(result["minimum_angle_argument"])
            and math.isfinite(result["maximum_angle_identity_residual"])
            and math.isfinite(result["maximum_imaginary_curvature"])
        )
        branch = bool(
            finite and result["minimum_angle_argument"] > 1e-10
            and result["maximum_angle_identity_residual"] < 5e-11
            and maximum_imaginary < 5e-9
        )
        return {
            "finite": finite,
            "branch": branch,
            "result": result,
            "gradient": gradient,
            "maximum_imaginary": maximum_imaginary,
            "maximum_imaginary_curvature": result["maximum_imaginary_curvature"],
        }
    except Exception as error:
        return {
            "finite": False,
            "branch": False,
            "error": repr(error),
            "gradient": None,
            "maximum_imaginary": math.inf,
            "maximum_imaginary_curvature": math.inf,
        }


def evaluate_equations(record, geometry, base, y):
    try:
        coordinates = coordinates_from_y(base, y)
        result = actions["evaluate_schedule"](record, geometry, coordinates)
        gradient = mp.matrix([result["gradient"][key] for key in INTERNAL_VARIABLES])
        maximum_imaginary = max(
            [abs(mp.im(result["action"])), abs(mp.im(result["gravitational"]))]
            + [abs(mp.im(value)) for value in gradient]
        )
        maximum_imaginary_curvature = abs(result["maximum_imaginary_curvature"])
        finite = all(
            mp.isfinite(value)
            for value in list(gradient) + [
                result["action"], result["minimum_angle_argument"],
                result["maximum_angle_identity_residual"], maximum_imaginary,
                maximum_imaginary_curvature,
            ]
        )
        branch = bool(
            finite
            and result["minimum_angle_argument"] > mp.mpf("1e-10")
            and result["maximum_angle_identity_residual"] < mp.mpf("1e-20")
            and maximum_imaginary < mp.mpf("1e-18")
        )
        return {
            "finite": finite,
            "branch": branch,
            "result": result,
            "gradient": gradient,
            "maximum_imaginary": maximum_imaginary,
            "maximum_imaginary_curvature": maximum_imaginary_curvature,
        }
    except Exception as error:
        return {
            "finite": False,
            "branch": False,
            "error": repr(error),
            "gradient": None,
            "maximum_imaginary": mp.inf,
            "maximum_imaginary_curvature": mp.inf,
        }


def real_numpy_gradient(state):
    return np.asarray([complex(value).real for value in state["gradient"]])


def solver_attempt(record, geometry, base, hessian, seed, lower, upper,
                   fast=False):
    started = time.perf_counter()
    diagnostics = {
        "evaluations": 0,
        "invalid_evaluations": 0,
        "minimum_angle_argument": math.inf,
        "maximum_angle_identity_residual": 0.0,
        "maximum_imaginary": 0.0,
        "maximum_imaginary_curvature": 0.0,
    }

    def residual(y):
        diagnostics["evaluations"] += 1
        state = (fast_equation_state if fast else evaluate_equations)(
            record, geometry, base, y
        )
        if not state["finite"] or not state["branch"]:
            diagnostics["invalid_evaluations"] += 1
            return np.full(10, 1e6)
        result = state["result"]
        diagnostics["minimum_angle_argument"] = min(
            diagnostics["minimum_angle_argument"],
            float(result["minimum_angle_argument"]),
        )
        diagnostics["maximum_angle_identity_residual"] = max(
            diagnostics["maximum_angle_identity_residual"],
            float(result["maximum_angle_identity_residual"]),
        )
        diagnostics["maximum_imaginary"] = max(
            diagnostics["maximum_imaginary"], float(state["maximum_imaginary"])
        )
        diagnostics["maximum_imaginary_curvature"] = max(
            diagnostics["maximum_imaginary_curvature"],
            float(state["maximum_imaginary_curvature"]),
        )
        return np.linalg.solve(hessian, real_numpy_gradient(state))

    clipped = np.minimum(np.maximum(np.asarray(seed, dtype=float), lower+1e-6), upper-1e-6)
    optimization = least_squares(
        residual,
        clipped,
        bounds=(lower, upper),
        method="trf",
        jac="3-point",
        diff_step=1e-5,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=1200,
        x_scale=1.0,
    )
    endpoint = np.asarray(optimization.x)
    final = (fast_equation_state if fast else evaluate_equations)(
        record, geometry, base, endpoint
    )
    if final["finite"] and final["branch"]:
        gradient = real_numpy_gradient(final)
        preconditioned = np.linalg.solve(hessian, gradient)
        gradient_norm = float(np.linalg.norm(gradient))
        residual_norm = float(np.linalg.norm(preconditioned))
        minimum_argument = float(final["result"]["minimum_angle_argument"])
        maximum_identity = float(final["result"]["maximum_angle_identity_residual"])
        maximum_imaginary = float(final["maximum_imaginary"])
        maximum_imaginary_curvature = float(final["maximum_imaginary_curvature"])
    else:
        gradient_norm = residual_norm = math.inf
        minimum_argument = 0.0
        maximum_identity = maximum_imaginary = math.inf
        maximum_imaginary_curvature = math.inf
    lower_distances = endpoint-lower
    upper_distances = upper-endpoint
    distance = float(min(np.min(lower_distances), np.min(upper_distances)))
    lower_active = [int(i) for i, value in enumerate(lower_distances) if value <= 1e-5]
    upper_active = [int(i) for i, value in enumerate(upper_distances) if value <= 1e-5]
    return {
        "seed": [sf(value) for value in seed],
        "clipped_seed": [sf(value) for value in clipped],
        "seed_was_clipped": bool(np.max(np.abs(clipped-np.asarray(seed))) > 0),
        "endpoint": [sf(value) for value in endpoint],
        "optimizer_success": bool(optimization.success),
        "optimizer_status": int(optimization.status),
        "optimizer_message": str(optimization.message),
        "nfev": int(optimization.nfev),
        "njev": None if optimization.njev is None else int(optimization.njev),
        "raw_gradient_norm": sf(gradient_norm),
        "preconditioned_residual_norm": sf(residual_norm),
        "branch_valid": bool(final["finite"] and final["branch"]),
        "distance_from_box": sf(distance),
        "active_lower_indices": lower_active,
        "active_upper_indices": upper_active,
        "minimum_angle_argument": sf(minimum_argument),
        "maximum_angle_identity_residual": sf(maximum_identity),
        "maximum_imaginary": sf(maximum_imaginary),
        "maximum_imaginary_curvature": sf(maximum_imaginary_curvature),
        "refinement_eligible": bool(
            final["finite"] and final["branch"]
            and distance > 1e-5 and residual_norm < 1e-7
        ),
        "_elapsed_seconds": sf(time.perf_counter()-started),
        "diagnostics": {
            key: sf(value) if isinstance(value, float) else value
            for key, value in diagnostics.items()
        },
        "_endpoint_array": endpoint,
        "_residual_norm_float": residual_norm,
    }


def synthetic_controls():
    diagonal = np.diag(np.arange(1.0, 11.0))
    target = np.asarray([.1, -.1, .05, -.05, .02, -.02, .2, -.2, .1, -.1])
    positive = least_squares(
        lambda x: diagonal@(x-target), np.zeros(10),
        bounds=(MAIN_LOWER, MAIN_UPPER), method="trf", jac="3-point",
        diff_step=1e-5, xtol=1e-12, ftol=1e-12, gtol=1e-12,
        max_nfev=1200, x_scale=1.0,
    )
    negative = least_squares(
        lambda x: np.exp(x)+1, np.zeros(10),
        bounds=(MAIN_LOWER, MAIN_UPPER), method="trf", jac="3-point",
        diff_step=1e-5, xtol=1e-12, ftol=1e-12, gtol=1e-12,
        max_nfev=1200, x_scale=1.0,
    )
    positive_norm = float(np.linalg.norm(diagonal@(positive.x-target)))
    negative_norm = float(np.linalg.norm(np.exp(negative.x)+1))
    return {
        "positive_root_error": sf(np.linalg.norm(positive.x-target)),
        "positive_residual_norm": sf(positive_norm),
        "negative_residual_norm": sf(negative_norm),
        "negative_classified_as_root": bool(negative_norm < 1e-10),
        "positive_ok": positive_norm < 1e-10,
        "negative_ok": negative_norm >= 1e-10,
    }


def mp_preconditioned_norm(gradient, hessian):
    return mp.norm(mp.lu_solve(hessian, mp.matrix([mp.re(value) for value in gradient])))


def mp_jacobian(record, geometry, base, y, step):
    matrix = mp.matrix(10, 10)
    for column in range(10):
        plus = list(y)
        minus = list(y)
        plus[column] += step
        minus[column] -= step
        upper = evaluate_equations(record, geometry, base, plus)
        lower = evaluate_equations(record, geometry, base, minus)
        if not (upper["finite"] and upper["branch"] and lower["finite"] and lower["branch"]):
            raise ValueError("Jacobian displacement left the branch")
        derivative = (upper["gradient"]-lower["gradient"])/(2*step)
        for row in range(10):
            matrix[row, column] = mp.re(derivative[row])
    return matrix


def in_main_interior(y, margin=mp.mpf("1e-5")):
    return all(
        mp.mpf(MAIN_LOWER[index])+margin < y[index]
        < mp.mpf(MAIN_UPPER[index])-margin
        for index in range(10)
    )


def refine_candidate(record, reverse_record, hessian_strings, endpoint):
    with mp.workdps(100):
        geometry = actions["exact_geometry"](100)
        base = actions["base_coordinates"](geometry)
        hessian = mp.matrix([[mp.mpf(value) for value in row] for row in hessian_strings])
        y = [mp.mpf(str(value)) for value in endpoint]
        history = []
        converged = False
        for iteration in range(40):
            state = evaluate_equations(record, geometry, base, y)
            if not (state["finite"] and state["branch"]):
                break
            norm_before = mp_preconditioned_norm(state["gradient"], hessian)
            history.append(mp_text(norm_before, 30))
            if norm_before < mp.mpf("1e-35"):
                converged = True
                break
            jacobian = mp_jacobian(record, geometry, base, y, mp.mpf("1e-18"))
            try:
                proposal = mp.lu_solve(
                    jacobian, -mp.matrix([mp.re(value) for value in state["gradient"]])
                )
            except Exception:
                break
            accepted = False
            for power in range(11):
                damping = mp.mpf(1)/(2**power)
                trial = [y[index]+damping*proposal[index] for index in range(10)]
                if not in_main_interior(trial):
                    continue
                trial_state = evaluate_equations(record, geometry, base, trial)
                if not (trial_state["finite"] and trial_state["branch"]):
                    continue
                trial_norm = mp_preconditioned_norm(trial_state["gradient"], hessian)
                if trial_norm < norm_before:
                    y = trial
                    accepted = True
                    break
            if not accepted:
                break

    with mp.workdps(140):
        geometry = actions["exact_geometry"](140)
        base = actions["base_coordinates"](geometry)
        hessian = mp.matrix([[mp.mpf(value) for value in row] for row in hessian_strings])
        state = evaluate_equations(record, geometry, base, y)
        reverse = evaluate_equations(reverse_record, geometry, base, y)
        analytic = state["gradient"] if state["finite"] else mp.matrix([mp.inf]*10)
        action_derivatives = []
        for index in range(10):
            values = []
            for step in (mp.mpf("1e-10"), mp.mpf("5e-11")):
                plus = list(y)
                minus = list(y)
                plus[index] += step
                minus[index] -= step
                upper = evaluate_equations(record, geometry, base, plus)
                lower = evaluate_equations(record, geometry, base, minus)
                values.append((upper["result"]["action"]-lower["result"]["action"])/(2*step))
            action_derivatives.append((4*values[1]-values[0])/3)
        action_derivatives = mp.matrix(action_derivatives)
        gradient_scale = max(mp.mpf(1), max(abs(value) for value in analytic))
        action_difference = max(
            abs(analytic[index]-action_derivatives[index]) for index in range(10)
        )/gradient_scale
        jacobian_a = mp_jacobian(record, geometry, base, y, mp.mpf("1e-15"))
        jacobian_b = mp_jacobian(record, geometry, base, y, mp.mpf("5e-16"))
        jacobian = mp.matrix([
            [mp.re(jacobian_b[row, column]+jacobian_b[column, row])/2
             for column in range(10)] for row in range(10)
        ])
        step_difference = max(
            abs(jacobian_a[row, column]-jacobian_b[row, column])
            for row in range(10) for column in range(10)
        )
        matrix_scale = max(
            mp.mpf(1), max(abs(jacobian[row, column]) for row in range(10) for column in range(10))
        )
        spectral_error = 10*(100*step_difference+mp.mpf("1e-50")*matrix_scale)
        eigenvalues, _ = mp.eigsy(jacobian)
        rank = sum(abs(eigenvalues[index]) > spectral_error for index in range(10))
        reverse_difference = max(
            abs(analytic[index]-reverse["gradient"][index]) for index in range(10)
        ) if reverse["finite"] else mp.inf
        preconditioned_norm = mp_preconditioned_norm(analytic, hessian)
        result = state["result"]
        accepted = bool(
            converged
            and state["finite"] and state["branch"] and reverse["finite"] and reverse["branch"]
            and preconditioned_norm < mp.mpf("1e-30")
            and action_difference < mp.mpf("1e-25")
            and state["maximum_imaginary"] < mp.mpf("1e-50")
            and result["minimum_angle_argument"] > mp.mpf("1e-8")
            and result["maximum_angle_identity_residual"] < mp.mpf("1e-50")
            and in_main_interior(y)
            and reverse_difference < mp.mpf("1e-40")
        )
        return {
            "endpoint": [mp_text(value, 65) for value in y],
            "newton_history": history,
            "newton_converged": converged,
            "preconditioned_residual_norm": mp_text(preconditioned_norm, 45),
            "analytic_action_relative_difference": mp_text(action_difference, 45),
            "maximum_imaginary": mp_text(state["maximum_imaginary"], 45),
            "maximum_imaginary_curvature": mp_text(
                state["maximum_imaginary_curvature"], 45
            ),
            "minimum_angle_argument": mp_text(result["minimum_angle_argument"], 45),
            "maximum_angle_identity_residual": mp_text(
                result["maximum_angle_identity_residual"], 45
            ),
            "reverse_equation_difference": mp_text(reverse_difference, 45),
            "root_jacobian_rank": rank,
            "root_jacobian_spectral_error": mp_text(spectral_error, 35),
            "root_jacobian_eigenvalues": [mp_text(eigenvalues[i], 45) for i in range(10)],
            "accepted_finite_positive_root": accepted,
        }


print("="*78)
print("REFINED H4 STATIONARY POSITIVE-ROOT SEARCH")
print("="*78)

previous_artifact = json.loads(OUTPUT.read_text()) if OUTPUT.exists() else None

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "all frozen action, Jacobian, prior-art and protocol inputs have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "29162db" and PROTOCOL_COMMIT == "722fb3c",
)

fill = json.loads(FILL.read_text())
jacobian_artifact = json.loads(JACOBIAN.read_text())
upstream_ok = check(
    "the two upstream outcomes and every schedule rank/inertia are exactly as required",
    fill["outcome"] == "REFINED_H4_INDUCED_FILL_OFF_SHELL"
    and fill["tests"] == {"passed": 12, "total": 12}
    and jacobian_artifact["outcome"]
        == "REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_MULTIPLE_CLASSES"
    and jacobian_artifact["tests"] == {"passed": 9, "total": 9}
    and jacobian_artifact["census"]["matrix_class_count"] == 12
    and all(item["rank"] == 10
            and item["inertia_positive_zero_negative"] == [8, 0, 2]
            for item in jacobian_artifact["census"]["schedules"]),
)
corrupted = np.asarray(
    [[float(value) for value in row]
     for row in jacobian_artifact["census"]["schedules"][0]["matrix"]]
)
corrupted[-1] = corrupted[0]
corruption_ok = check(
    "the required full-rank gate rejects a deliberately rank-nine matrix",
    np.linalg.matrix_rank(corrupted) < 10,
)

actions = load_action_definitions()
definitions_ok = check(
    "only frozen AST definitions are loaded and the required action functions exist",
    {"exact_geometry", "base_coordinates", "evaluate_schedule"} <= set(actions)
    and "OUTPUT" not in actions,
)
combinatorics = [parse_combinatorics(record) for record in fill["combinatorics"]]
orders = tuple(record["order"] for record in combinatorics)
order_to_index = {order: index for index, order in enumerate(orders)}
combinatorics_ok = check(
    "all 24 committed reduced carriers reconstruct without a new schedule choice",
    orders == tuple(permutations(range(4)))
    and all(record["pentachora"] == 57600 and record["triangles"] == 149280
            and record["mixed_triangle_types"] == 0
            for record in combinatorics),
)

schedule_by_order = {
    tuple(item["order"]): item for item in jacobian_artifact["census"]["schedules"]
}
class_records = []
for item in jacobian_artifact["census"]["matrix_classes"]:
    representative = tuple(item["orders"][0])
    reverse = tuple(reversed(representative))
    class_records.append({
        "class": item["class"],
        "order": representative,
        "reverse": reverse,
        "record": combinatorics[order_to_index[representative]],
        "reverse_record": combinatorics[order_to_index[reverse]],
        "matrix_record": schedule_by_order[representative],
    })

anchors = (
    [0.0]*10,
    [1e-6, -1e-6, 2e-6, -2e-6, 3e-6, -3e-6,
     -1e-4, -1e-4, -1e-4, -1e-4],
    [-3e-6, 2e-6, -1e-6, 1e-6, -2e-6, 3e-6,
     -5e-5, -1e-4, -1.5e-4, -2e-4],
)
with mp.workdps(80):
    geometry80 = actions["exact_geometry"](80)
    base80 = actions["base_coordinates"](geometry80)
    reversal_differences = []
    reversal_branches = []
    reversal_anchor_diagnostics = []
    for class_record in class_records:
        for anchor_index, anchor in enumerate(anchors):
            left = evaluate_equations(class_record["record"], geometry80, base80, anchor)
            right = evaluate_equations(class_record["reverse_record"], geometry80, base80, anchor)
            reversal_branches.append(
                left["finite"] and left["branch"] and right["finite"] and right["branch"]
            )
            if reversal_branches[-1]:
                difference = max(
                    [abs(left["result"]["action"]-right["result"]["action"])]
                    + [abs(left["gradient"][i]-right["gradient"][i]) for i in range(10)]
                )
            else:
                difference = mp.inf
            reversal_differences.append(difference)
            reversal_anchor_diagnostics.append({
                "class": class_record["class"],
                "anchor": anchor_index,
                "left_branch": bool(left["finite"] and left["branch"]),
                "right_branch": bool(right["finite"] and right["branch"]),
                "left_minimum_argument": mp_text(
                    left["result"]["minimum_angle_argument"], 25
                ) if left["finite"] else "nonfinite",
                "right_minimum_argument": mp_text(
                    right["result"]["minimum_angle_argument"], 25
                ) if right["finite"] else "nonfinite",
                "left_maximum_imaginary": mp_text(left["maximum_imaginary"], 25),
                "right_maximum_imaginary": mp_text(right["maximum_imaginary"], 25),
                "left_maximum_imaginary_curvature": mp_text(
                    left["maximum_imaginary_curvature"], 25
                ),
                "right_maximum_imaginary_curvature": mp_text(
                    right["maximum_imaginary_curvature"], 25
                ),
                "difference": mp_text(difference, 25),
            })
time_reversal_ok = check(
    "all 12 time-reversal pairs agree at all three preregistered nonlinear anchors",
    all(reversal_branches) and max(reversal_differences) < mp.mpf("1e-60"),
    f"max difference={mp_text(max(reversal_differences), 8)}",
)

fast_geometry = {"mass": float(geometry80["mass"])}
fast_base = {key: float(value) for key, value in base80.items()}
fast_comparisons = []
fast_control_values = []
with mp.workdps(80):
    for class_record in class_records:
        proposal = [
            mp.mpf(value) for value in
            class_record["matrix_record"]["unapplied_newton_proposal"]
        ]
        comparison_anchors = (
            anchors[0], anchors[1], anchors[2],
            [mp.mpf("0.25")*value for value in proposal],
            [mp.mpf("0.5")*value for value in proposal],
            proposal,
        )
        for anchor_index, anchor in enumerate(comparison_anchors):
            precise = evaluate_equations(
                class_record["record"], geometry80, base80, anchor
            )
            fast = fast_equation_state(
                class_record["record"], fast_geometry, fast_base, anchor
            )
            fast_reverse = fast_equation_state(
                class_record["reverse_record"], fast_geometry, fast_base, anchor
            )
            if precise["finite"] and fast["finite"] and fast_reverse["finite"]:
                action_error = abs(
                    fast["result"]["action"]-complex(precise["result"]["action"])
                )/max(1, abs(complex(precise["result"]["action"])))
                precise_gradient = np.asarray([
                    complex(value) for value in precise["gradient"]
                ])
                gradient_error = float(np.max(np.abs(
                    fast["gradient"]-precise_gradient
                )))/max(1, float(np.max(np.abs(precise_gradient))))
                argument_error = abs(
                    fast["result"]["minimum_angle_argument"]
                    - float(precise["result"]["minimum_angle_argument"])
                )
                reverse_error = max(
                    abs(fast["result"]["action"]-fast_reverse["result"]["action"]),
                    float(np.max(np.abs(
                        fast["gradient"]-fast_reverse["gradient"]
                    ))),
                )
            else:
                action_error = gradient_error = argument_error = reverse_error = math.inf
            passed_fast = bool(
                precise["branch"] and fast["branch"] and fast_reverse["branch"]
                and action_error < 5e-9 and gradient_error < 5e-9
                and argument_error < 5e-10
                and fast["result"]["maximum_angle_identity_residual"] < 5e-11
                and fast["maximum_imaginary"] < 5e-9
                and reverse_error < 5e-9
            )
            fast_control_values.append(passed_fast)
            fast_comparisons.append({
                "class": class_record["class"],
                "anchor": anchor_index,
                "precise_branch": bool(precise["finite"] and precise["branch"]),
                "fast_branch": bool(fast["finite"] and fast["branch"]),
                "action_relative_error": sf(action_error),
                "gradient_relative_error": sf(gradient_error),
                "minimum_argument_error": sf(argument_error),
                "fast_angle_identity_residual": sf(
                    fast["result"]["maximum_angle_identity_residual"]
                    if fast["finite"] else math.inf
                ),
                "fast_physical_imaginary": sf(fast["maximum_imaginary"]),
                "fast_time_reversal_difference": sf(reverse_error),
                "passed": passed_fast,
            })
fast_evaluator_ok = check(
    "the binary64 search evaluator matches the 80-decimal action at all 72 anchors",
    len(fast_comparisons) == 72 and all(fast_control_values),
    "max errors action={}, gradient={}, argument={}".format(
        sf(max(float(item["action_relative_error"]) for item in fast_comparisons)),
        sf(max(float(item["gradient_relative_error"]) for item in fast_comparisons)),
        sf(max(float(item["minimum_argument_error"]) for item in fast_comparisons)),
    ),
)

synthetic = synthetic_controls()
synthetic_positive_ok = check(
    "the frozen solver recovers its synthetic known interior root",
    synthetic["positive_ok"], synthetic["positive_residual_norm"],
)
synthetic_negative_ok = check(
    "the frozen solver does not call its synthetic no-root system a root",
    synthetic["negative_ok"], synthetic["negative_residual_norm"],
)

controls_before_search = all((
    provenance_ok, upstream_ok, corruption_ok, definitions_ok,
    combinatorics_ok, time_reversal_ok, synthetic_positive_ok,
    synthetic_negative_ok, fast_evaluator_ok,
))

main_attempts = []
ladder_attempts = []
refinements = []


def write_checkpoint(stage):
    def public(item):
        return {key: value for key, value in item.items() if not key.startswith("_")}
    checkpoint = {
        "title": "INCOMPLETE checkpoint: refined H4 stationary root search",
        "date": "2026-08-20",
        "stage": stage,
        "protocol_commit": PROTOCOL_COMMIT,
        "main_attempts": [public(item) for item in main_attempts],
        "boundary_ladder_attempts": [public(item) for item in ladder_attempts],
        "outcome": "INCOMPLETE",
    }
    OUTPUT.write_text(json.dumps(checkpoint, indent=2, sort_keys=True)+"\n")


if controls_before_search:
    print("[INFO] running 72 preregistered main-box attempts", flush=True)
    for class_position, class_record in enumerate(class_records):
        matrix = np.asarray(
            [[float(value) for value in row]
             for row in class_record["matrix_record"]["matrix"]]
        )
        proposal = np.asarray([
            float(value) for value in
            class_record["matrix_record"]["unapplied_newton_proposal"]
        ])
        s4 = np.asarray([.05, -.05, .05, -.05, .05, -.05, -1.5, -2, -2.5, -3])
        s5 = np.asarray([-.05, .05, -.05, .05, -.05, .05, -3, -2.5, -2, -1.5])
        seeds = (
            np.zeros(10), .5*proposal, proposal,
            np.asarray([0.0]*6+[-4.0]*4), s4, s5,
        )
        for seed_index, seed in enumerate(seeds):
            attempt = solver_attempt(
                class_record["record"], fast_geometry, fast_base, matrix,
                seed, MAIN_LOWER, MAIN_UPPER, fast=True,
            )
            attempt.update({"class": class_record["class"], "seed_index": seed_index})
            main_attempts.append(attempt)
        write_checkpoint(f"main_class_{class_position}")
        print(f"[INFO] main classes completed: {class_position+1}/12", flush=True)

    print("[INFO] running 48 zero-lapse boundary-ladder attempts", flush=True)
    for class_position, class_record in enumerate(class_records):
        matrix = np.asarray(
            [[float(value) for value in row]
             for row in class_record["matrix_record"]["matrix"]]
        )
        proposal = .5*np.asarray([
            float(value) for value in
            class_record["matrix_record"]["unapplied_newton_proposal"]
        ])
        for lower_lapse in LADDER_LOWERS:
            lower = np.asarray([-0.35]*6+[lower_lapse]*4)
            attempt = solver_attempt(
                class_record["record"], fast_geometry, fast_base, matrix,
                proposal, lower, MAIN_UPPER, fast=True,
            )
            attempt.update({
                "class": class_record["class"],
                "lower_lapse_bound": lower_lapse,
            })
            ladder_attempts.append(attempt)
        write_checkpoint(f"ladder_class_{class_position}")
        print(f"[INFO] ladder classes completed: {class_position+1}/12", flush=True)

    print("[INFO] refining every eligible interior endpoint", flush=True)
    for attempt in main_attempts:
        if attempt["refinement_eligible"]:
            class_record = class_records[attempt["class"]]
            refinement = refine_candidate(
                class_record["record"], class_record["reverse_record"],
                class_record["matrix_record"]["matrix"],
                attempt["_endpoint_array"],
            )
            refinement.update({"class": attempt["class"], "seed_index": attempt["seed_index"]})
            refinements.append(refinement)

search_count_ok = check(
    "the complete frozen look-elsewhere census contains exactly 72+48 attempts",
    (not controls_before_search)
    or (len(main_attempts) == 72 and len(ladder_attempts) == 48),
    f"main={len(main_attempts)}, ladder={len(ladder_attempts)}",
)

accepted_by_class = {index: [] for index in range(12)}
for refinement in refinements:
    if refinement["accepted_finite_positive_root"]:
        accepted_by_class[refinement["class"]].append(refinement)

boundary_patterns = {}
for class_index in range(12):
    records = sorted(
        (item for item in ladder_attempts if item["class"] == class_index),
        key=lambda item: item["lower_lapse_bound"], reverse=True,
    )
    norms = [item["_residual_norm_float"] for item in records]
    active = [any(index >= 6 for index in item["active_lower_indices"]) for item in records]
    slopes = []
    for left, right in zip(records, records[1:]):
        left_norm = left["_residual_norm_float"]
        right_norm = right["_residual_norm_float"]
        if left_norm > 0 and right_norm > 0 and math.isfinite(left_norm+right_norm):
            slopes.append(math.log(right_norm/left_norm)/(
                right["lower_lapse_bound"]-left["lower_lapse_bound"]
            ))
        else:
            slopes.append(math.nan)
    boundary_patterns[class_index] = {
        "passes": bool(
            not accepted_by_class[class_index]
            and len(records) == 4 and all(active)
            and all(norms[index+1] < norms[index] for index in range(3))
        ),
        "norms_by_lower_bound": {
            str(item["lower_lapse_bound"]): sf(item["_residual_norm_float"])
            for item in records
        },
        "pairwise_log_slopes": [sf(value) for value in slopes],
    }

accepted_classes = [index for index, values in accepted_by_class.items() if values]
classification_ok = check(
    "every accepted finite root passes all preregistered high-precision gates",
    all(item["accepted_finite_positive_root"]
        for values in accepted_by_class.values() for item in values),
    f"accepted classes={accepted_classes}, refined endpoints={len(refinements)}",
)

all_controls = controls_before_search and search_count_ok and classification_ok
if not all_controls:
    outcome = "REFINED_H4_STATIONARY_ROOT_CONTROL_FAILED"
elif len(accepted_classes) == 12:
    outcome = "REFINED_H4_FINITE_ROOTS_ALL_CLASSES"
elif accepted_classes:
    outcome = "REFINED_H4_FINITE_ROOTS_SOME_CLASSES"
elif all(item["passes"] for item in boundary_patterns.values()):
    outcome = "REFINED_H4_NO_FINITE_ROOT_FOUND_ZERO_LAPSE_PATTERN"
else:
    outcome = "REFINED_H4_NO_FINITE_ROOT_FOUND_OTHER"
outcome_ok = check(
    "the frozen outcome hierarchy classifies the bounded root search",
    outcome in {
        "REFINED_H4_STATIONARY_ROOT_CONTROL_FAILED",
        "REFINED_H4_FINITE_ROOTS_ALL_CLASSES",
        "REFINED_H4_FINITE_ROOTS_SOME_CLASSES",
        "REFINED_H4_NO_FINITE_ROOT_FOUND_ZERO_LAPSE_PATTERN",
        "REFINED_H4_NO_FINITE_ROOT_FOUND_OTHER",
    },
    outcome,
)


def deterministic_science(record):
    def clean(value):
        if isinstance(value, dict):
            return {
                key: clean(item) for key, item in value.items()
                if key != "elapsed_seconds" and not key.startswith("_elapsed")
            }
        if isinstance(value, list):
            return [clean(item) for item in value]
        return value
    return clean({
        key: record[key] for key in
        ("definitions", "controls", "search", "scope", "outcome")
    })

def strip_private(attempt):
    return {key: value for key, value in attempt.items() if not key.startswith("_")}


artifact = {
    "title": "Bounded positive-root search for the refined H4 slab",
    "date": "2026-08-20",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "classes": 12,
        "main_attempts": 72,
        "boundary_ladder_attempts": 48,
        "look_elsewhere_denominator": 120,
        "main_cross_bounds": [-.35, .35],
        "main_lapse_bounds": [-8, 2],
        "boundary_ladder_lower_bounds": list(LADDER_LOWERS),
        "absolute_tick_target_loaded": False,
        "effective_boundary_hessian_computed": False,
    },
    "controls": {
        "synthetic_solver": synthetic,
        "fast_evaluator_comparisons": fast_comparisons,
        "maximum_time_reversal_anchor_difference": mp_text(
            max(reversal_differences), 45
        ),
        "time_reversal_anchor_count": len(reversal_differences),
        "time_reversal_anchor_diagnostics": reversal_anchor_diagnostics,
    },
    "search": {
        "main_attempts": [strip_private(item) for item in main_attempts],
        "boundary_ladder_attempts": [strip_private(item) for item in ladder_attempts],
        "refinements": refinements,
        "accepted_class_count": len(accepted_classes),
        "accepted_classes": accepted_classes,
        "accepted_attempt_hit_fraction": {
            "numerator": sum(
                item["accepted_finite_positive_root"] for item in refinements
            ),
            "denominator": 120,
        },
        "zero_lapse_boundary_patterns": {
            str(key): value for key, value in boundary_patterns.items()
        },
    },
    "scope": {
        "bounded_search_is_nonexistence_proof": False,
        "zero_lapse_pattern_is_physical_tick": False,
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
previous_science_ok = check(
    "all scientific fields reproduce the first complete artifact exactly",
    previous_artifact is not None
    and previous_artifact.get("outcome") == outcome
    and deterministic_science(previous_artifact) == deterministic_science(artifact),
)
artifact["tests"] = {"passed": passed, "total": tests}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"ACCEPTED ROOT CLASSES: {len(accepted_classes)}/12")
print(f"REFINED ENDPOINTS: {len(refinements)}")
print(f"ZERO-LAPSE PATTERN CLASSES: {sum(v['passes'] for v in boundary_patterns.values())}/12")
print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} checks passed")
sys.exit(0 if outcome_ok and passed == tests else 1)
