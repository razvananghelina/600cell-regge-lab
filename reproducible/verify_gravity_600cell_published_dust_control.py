#!/usr/bin/env python3
"""External control against the published 600-cell dust sandwich.

Protocol commit: cc0902b.  The published constants are reconstructed from
their formulas.  No parameter is optimized and no stationary-root search is
performed.
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
import sys

import mpmath as arb
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_published_dust_control.json"
PROTOCOL_COMMIT = "cc0902b"
PRIOR_ART_COMMIT = "e7d8bd5"
UPSTREAM_COMMIT = "b11185e"
PRECISION_CORRECTION_COMMIT = "3056c7d"
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


# Recheck the complete boundary-Legendre construction before extending it.
spec = importlib.util.spec_from_file_location(
    "global_boundary_legendre",
    HERE / "verify_gravity_global_boundary_legendre.py",
)
bl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = bl
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(bl)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise

check(
    "the imported 65-variable action retains all 33 upstream certificates",
    bl.tests == bl.passed == 33,
)


print("=" * 78)
print("PUBLISHED 600-CELL DUST-SANDWICH CONTROL")
print("=" * 78)

# De Felice--Fabri data reconstructed from the displayed formulas.
M_STAR = 10.0
ZETA = (math.pi**2*math.sqrt(2.0)/50.0)**(1.0/3.0)
R0 = 4.0*M_STAR/(3.0*math.pi)
L0 = ZETA*R0
L0_SQUARE = L0**2
EPSILON_3 = 2.0*math.pi-5.0*math.acos(1.0/3.0)
MASS = (90.0/math.pi)*EPSILON_3*L0
TAU = 0.0102
TAU_SQUARE = TAU**2
SLANT_SQUARE = L0_SQUARE-TAU_SQUARE

published_controls = {
    "l0_square": 7.69379990138304,
    "l0_square_minus_tau_square": 7.69369586138304,
    "solved_slant_square": 7.69369586138301,
    "solved_final_square": 7.69379990138297,
}
published_errors = {
    "l0_square": abs(L0_SQUARE-published_controls["l0_square"]),
    "l0_square_minus_tau_square": abs(
        SLANT_SQUARE-published_controls["l0_square_minus_tau_square"]
    ),
    "solved_slant_square": abs(
        SLANT_SQUARE-published_controls["solved_slant_square"]
    ),
    "solved_final_square": abs(
        L0_SQUARE-published_controls["solved_final_square"]
    ),
}
check(
    "the unrounded source formulas reproduce the paper's printed sandwich",
    max(published_errors.values()) < 2e-13 and 10.202 < MASS < 10.203,
    f"max squared-length error={max(published_errors.values()):.3e}, M={MASS:.12f}",
)
check(
    "the published slant relation is imposed without fitting",
    abs(SLANT_SQUARE-(L0_SQUARE-TAU_SQUARE)) < 1e-15,
    f"l0^2={L0_SQUARE:.14f}, d^2={SLANT_SQUARE:.14f}, tau^2={TAU_SQUARE:.8f}",
)


def dust_action_and_gradient(variables):
    poles = variables[30:35]
    action = -(8.0*math.pi*MASS/5.0)*float(np.sum(np.sqrt(poles)))
    gradient = np.zeros(65)
    gradient[30:35] = -(4.0*math.pi*MASS)/(5.0*np.sqrt(poles))
    return action, gradient


def total_reduced_evaluation(model, variables, old_values):
    action, gradient, old_gradient, data = bl.reduced_evaluation(
        model, variables, old_values
    )
    dust_action, dust_gradient = dust_action_and_gradient(variables)
    return action+dust_action, gradient+dust_gradient, old_gradient, data


def phase_pair_variables(model):
    grouped = defaultdict(list)
    for orbit in model["diagonal_orbits"]:
        representative = min(orbit)
        variable = model["edge_to_variable"][representative]
        logical = tuple(vertex % 120 for vertex in representative)
        pair = tuple(sorted(model["phase"][vertex] for vertex in logical))
        grouped[pair].append(variable)
    return {pair: tuple(sorted(indices)) for pair, indices in grouped.items()}


# Independent 60-decimal action-only evaluator frozen in correction 3056c7d.
arb.mp.dps = 60
ARB_I = arb.mpc(0, 1)
ARB_TAU = arb.mpf("0.0102")
ARB_M_STAR = arb.mpf(10)
ARB_ZETA = (arb.pi**2*arb.sqrt(2)/50)**(arb.mpf(1)/3)
ARB_R0 = 4*ARB_M_STAR/(3*arb.pi)
ARB_L0 = ARB_ZETA*ARB_R0
ARB_L0_SQUARE = ARB_L0**2
ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
ARB_MASS = (90/arb.pi)*ARB_EPSILON_3*ARB_L0
ARB_TAU_SQUARE = ARB_TAU**2
ARB_SLANT_SQUARE = ARB_L0_SQUARE-ARB_TAU_SQUARE
ARB_OLD_VALUES = [ARB_L0_SQUARE]*30
ARB_BASE_VARIABLES = (
    [ARB_SLANT_SQUARE]*30
    + [ARB_TAU_SQUARE]*5
    + [ARB_L0_SQUARE]*30
)


def arb_edge_square(model, edge, variables, old_values=ARB_OLD_VALUES):
    edge = tuple(sorted(edge))
    if edge in model["all_edge_to_variable"]:
        variable = model["all_edge_to_variable"][edge]
        jacobian = arb.mpf(int(model["all_edge_jacobian"][edge]))
        return jacobian*variables[variable]
    if edge in model["old_to_orbit"]:
        return old_values[model["old_to_orbit"][edge]]
    raise ValueError(f"edge absent from arbitrary-precision slab: {edge}")


def arb_simplex_squared(model, simplex, variables):
    squared = [[arb.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value = arb_edge_square(
            model, (simplex[left], simplex[right]), variables
        )
        squared[left][right] = squared[right][left] = value
    return squared


def arb_signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices)-1
    if dimension == 0:
        return arb.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = arb.matrix([
        [
            (squared[base][left]+squared[base][right]-squared[left][right])/2
            for right in others
        ]
        for left in others
    ])
    return arb.det(gram)/(arb.factorial(dimension)**2)


def arb_log_minus(value):
    scale = max(arb.mpf(1), abs(value))
    if abs(arb.im(value)) < arb.mpf("1e-50")*scale:
        real = arb.re(value)
        if real < 0:
            return arb.log(-real)-ARB_I*arb.pi
        return arb.log(real)
    return arb.log(value)


def arb_angle_data(squared):
    gram = arb.matrix([
        [
            (squared[0][left]+squared[0][right]-squared[left][right])/2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram**-1
    simplex_volume_square = arb_signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: arb_signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    angles = {}
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge_vertices = tuple(
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = arb_signed_volume_square(
            squared, hinge_vertices
        )
        gram_derivative = arb.matrix(4, 4)
        opposite_edge = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                hit_0_left = int({0, left} == opposite_edge)
                hit_0_right = int({0, right} == opposite_edge)
                hit_pair = int(left != right and {left, right} == opposite_edge)
                gram_derivative[left-1, right-1] = (
                    hit_0_left+hit_0_right-hit_pair
                )/2
        product = inverse*gram_derivative
        volume_derivative = simplex_volume_square*sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            arb.sqrt(arb.mpc(facet_volume_squares[omitted_a]))
            * arb.sqrt(arb.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16*volume_derivative/denominator
        sine = -arb.mpf(4)/3*(
            arb.sqrt(arb.mpc(hinge_volume_square))
            * arb.sqrt(arb.mpc(simplex_volume_square))
        )/denominator
        argument = cosine+ARB_I*sine
        angles[hinge_vertices] = -ARB_I*arb_log_minus(argument)
    return angles


def arb_triangle_area_square(values):
    x, y, z = values
    return (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16


def arb_gravitational_action(model, variables):
    triangle_orbits = model["triangle_orbits"]
    curvature = [
        arb.pi if min(orbit) in model["boundary_triangles"] else 2*arb.pi
        for orbit in triangle_orbits
    ]
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        squared = arb_simplex_squared(model, simplex, variables)
        for local_hinge, angle in arb_angle_data(squared).items():
            triangle = tuple(
                sorted(simplex[position] for position in local_hinge)
            )
            curvature[model["triangle_to_orbit"][triangle]] += angle
    action_sum = arb.mpc(0)
    for orbit_index, orbit in enumerate(triangle_orbits):
        triangle = min(orbit)
        edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
        values = tuple(
            arb_edge_square(model, edge, variables) for edge in edges
        )
        area = arb.sqrt(arb.mpc(arb_triangle_area_square(values)))
        action_sum += 24*area*curvature[orbit_index]
    return -ARB_I*action_sum


def arb_action_components(model, variables):
    gravitational = arb_gravitational_action(model, variables)
    dust = -(8*arb.pi*ARB_MASS/5)*sum(
        arb.sqrt(value) for value in variables[30:35]
    )
    return gravitational, dust, gravitational+dust


old_values = np.full(30, L0_SQUARE)
base_variables = np.concatenate((
    np.full(30, SLANT_SQUARE),
    np.full(5, TAU_SQUARE),
    np.full(30, L0_SQUARE),
))

records = {}
for parity, model in bl.models.items():
    reduced_action, reduced_gradient, reduced_old, reduced_data = (
        bl.reduced_evaluation(model, base_variables, old_values)
    )
    full_action, full_gradient, full_old, full_data = bl.full_evaluation(
        model, base_variables, old_values
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
    imaginary_residual = max(
        abs(reduced_action.imag),
        float(np.max(np.abs(reduced_gradient.imag))),
        float(np.max(np.abs(reduced_old.imag))),
        abs(full_action.imag),
        float(np.max(np.abs(full_gradient.imag))),
        float(np.max(np.abs(full_old.imag))),
    )
    check(
        f"{parity}: full and 100-orbit gravitational actions agree at the published point",
        action_error < 3e-8
        and gradient_error < 3e-8
        and old_gradient_error < 3e-8,
        f"action={action_error:.3e}, grad65={gradient_error:.3e}, old30={old_gradient_error:.3e}",
    )
    check(
        f"{parity}: the published sandwich is nondegenerate Lorentzian on the real branch",
        reduced_data["negative_counts"] == Counter({1: 100})
        and full_data["negative_counts"] == Counter({1: 2400})
        and reduced_data["minimum_gram"] > 1e-8
        and reduced_data["minimum_argument"] > 1e-6
        and imaginary_residual < 3e-7,
        f"gram={reduced_data['minimum_gram']:.3e}, argument={reduced_data['minimum_argument']:.3e}, imag={imaginary_residual:.3e}",
    )

    dust_action, dust_gradient = dust_action_and_gradient(base_variables)
    total_action = reduced_action+dust_action
    total_gradient = reduced_gradient+dust_gradient
    per_edge = total_gradient[:35].real/24.0
    diagonal = per_edge[:30]
    poles = per_edge[30:35]
    gravitational_pole_length_lhs = (
        2.0*TAU*reduced_gradient[30:35].real/24.0
    )
    pole_length_rhs = math.pi*MASS/15.0
    pole_balance_error = max(
        relative_error(value, pole_length_rhs)
        for value in gravitational_pole_length_lhs
    )

    grouped = phase_pair_variables(model)
    phase_pairs = {}
    for pair, variables in sorted(grouped.items()):
        values = np.array([diagonal[variable] for variable in variables])
        phase_pairs[f"{pair[0]}-{pair[1]}"] = {
            "variables": list(variables),
            "residuals": values,
            "sum": float(np.sum(values)),
            "spread": float(np.max(values)-np.min(values)),
        }

    records[parity] = {
        "model": model,
        "reduced_action": reduced_action,
        "total_action": total_action,
        "reduced_gradient": reduced_gradient,
        "total_gradient": total_gradient,
        "per_edge": per_edge,
        "diagonal": diagonal,
        "poles": poles,
        "phase_pairs": phase_pairs,
        "action_error": action_error,
        "gradient_error": gradient_error,
        "old_gradient_error": old_gradient_error,
        "imaginary_residual": imaginary_residual,
        "minimum_gram": reduced_data["minimum_gram"],
        "minimum_argument": reduced_data["minimum_argument"],
        "pole_length_lhs": gravitational_pole_length_lhs,
        "pole_length_rhs": pole_length_rhs,
        "pole_balance_error": pole_balance_error,
    }

check(
    "each parity has exactly three diagonal orbits for every phase pair",
    all(
        len(record["phase_pairs"]) == 10
        and {len(item["variables"]) for item in record["phase_pairs"].values()} == {3}
        for record in records.values()
    ),
)


# Centered differences of the complete total restricted action.  The direct
# 2400-simplex comparison above proves that this orbit sum is the complete
# action restriction; the differences are independent of analytic area
# derivatives.
_WORKER_MODEL = None


def initialize_worker(model):
    global _WORKER_MODEL
    _WORKER_MODEL = model


def total_action_worker(variables):
    return total_reduced_evaluation(_WORKER_MODEL, variables, old_values)[0]


direct_step = 3e-6
fork_context = mp.get_context("fork")
for parity, record in records.items():
    model = record["model"]
    points = []
    deltas = []
    for variable in range(35):
        delta = direct_step*base_variables[variable]
        plus = base_variables.copy()
        minus = base_variables.copy()
        plus[variable] += delta
        minus[variable] -= delta
        points.extend((plus, minus))
        deltas.append(delta)
    with fork_context.Pool(
        processes=8, initializer=initialize_worker, initargs=(model,)
    ) as pool:
        actions = pool.map(total_action_worker, points, chunksize=2)
    direct_gradient = np.array([
        (actions[2*variable]-actions[2*variable+1])/(2.0*deltas[variable])
        for variable in range(35)
    ])
    analytic = record["total_gradient"][:35]
    direct_error = max(
        relative_error(left, right)
        for left, right in zip(direct_gradient, analytic)
    )
    direct_imaginary = float(np.max(np.abs(direct_gradient.imag)))
    record["direct_gradient"] = direct_gradient
    record["direct_error"] = direct_error
    record["direct_imaginary"] = direct_imaginary
    binary64_pass = bool(
        direct_error < 5e-5 and direct_imaginary < 3e-7
    )
    record["binary64_direct_pass"] = binary64_pass
    print(
        f"[{'PASS' if binary64_pass else 'RECORDED PRECISION LIMIT'}] "
        f"{parity}: binary64 total-action differences"
    )
    print(f"       relative error={direct_error:.3e}, imag={direct_imaginary:.3e}")


# Correction 3056c7d: repeat the same action differences, at the same points
# and step, with an independent action-only implementation at 60 decimals.
_ARB_WORKER_MODEL = None


def initialize_arb_worker(model):
    global _ARB_WORKER_MODEL
    _ARB_WORKER_MODEL = model
    arb.mp.dps = 60


def arb_action_worker(variables):
    return arb_action_components(_ARB_WORKER_MODEL, variables)


def arb_relative_error(left, right):
    return abs(left-right)/max(arb.mpf(1), abs(left), abs(right))


arb_step = arb.mpf("3e-6")
for parity, record in records.items():
    model = record["model"]
    base_components = arb_action_components(model, ARB_BASE_VARIABLES)
    double_total = arb.mpc(
        str(record["total_action"].real), str(record["total_action"].imag)
    )
    base_action_error = float(arb_relative_error(base_components[2], double_total))
    check(
        f"{parity}: the 60-decimal action reproduces the certified base action",
        base_action_error < 3e-8,
        f"relative error={base_action_error:.3e}",
    )

    points = []
    deltas = []
    for variable in range(35):
        delta = arb_step*ARB_BASE_VARIABLES[variable]
        plus = list(ARB_BASE_VARIABLES)
        minus = list(ARB_BASE_VARIABLES)
        plus[variable] += delta
        minus[variable] -= delta
        points.extend((plus, minus))
        deltas.append(delta)
    with fork_context.Pool(
        processes=8, initializer=initialize_arb_worker, initargs=(model,)
    ) as pool:
        components = pool.map(arb_action_worker, points, chunksize=1)

    arb_gravitational_gradient = []
    arb_total_gradient = []
    for variable, delta in enumerate(deltas):
        plus = components[2*variable]
        minus = components[2*variable+1]
        arb_gravitational_gradient.append((plus[0]-minus[0])/(2*delta))
        arb_total_gradient.append((plus[2]-minus[2])/(2*delta))

    analytic = [
        arb.mpc(str(value.real), str(value.imag))
        for value in record["total_gradient"][:35]
    ]
    analytic_error = float(max(
        arb_relative_error(left, right)
        for left, right in zip(arb_total_gradient, analytic)
    ))
    imaginary_derivative = float(max(
        abs(arb.im(value)) for value in arb_total_gradient
    ))
    direct_per_edge_maximum = float(max(
        abs(value)/24 for value in arb_total_gradient
    ))
    check(
        f"{parity}: 60-decimal differences reproduce all 35 sourced derivatives",
        analytic_error < 5e-8 and imaginary_derivative < 1e-30,
        f"analytic error={analytic_error:.3e}, imag={imaginary_derivative:.3e}",
    )
    check(
        f"{parity}: direct 60-decimal residuals pass the frozen stationarity gate",
        direct_per_edge_maximum <= 1e-7,
        f"max per-edge residual={direct_per_edge_maximum:.3e}",
    )

    dust_pole_derivative = [
        -(4*arb.pi*ARB_MASS)/(5*arb.sqrt(ARB_BASE_VARIABLES[variable]))
        for variable in range(30, 35)
    ]
    cancellation_ratios = [
        abs(arb_total_gradient[variable])
        / max(
            abs(arb_gravitational_gradient[variable]),
            abs(dust_pole_derivative[variable-30]),
        )
        for variable in range(30, 35)
    ]
    record.update({
        "arb_base_action": base_components[2],
        "arb_base_action_error": base_action_error,
        "arb_gravitational_gradient": arb_gravitational_gradient,
        "arb_total_gradient": arb_total_gradient,
        "arb_analytic_error": analytic_error,
        "arb_imaginary_derivative": imaginary_derivative,
        "arb_direct_per_edge_maximum": direct_per_edge_maximum,
        "arb_dust_pole_derivative": dust_pole_derivative,
        "arb_cancellation_ratios": cancellation_ratios,
    })


stationarity_threshold = 1e-7
for parity, record in records.items():
    arb_per_edge = [value/24 for value in record["arb_total_gradient"]]
    arb_diagonal = arb_per_edge[:30]
    arb_poles = arb_per_edge[30:35]
    pole_pass = bool(max(abs(value) for value in arb_poles) <= stationarity_threshold)
    diagonal_pass = bool(
        max(abs(value) for value in arb_diagonal) <= stationarity_threshold
    )
    arb_phase_sums = {}
    for pair, item in record["phase_pairs"].items():
        arb_phase_sums[pair] = sum(
            arb_diagonal[variable] for variable in item["variables"]
        )
    phase_sum_pass = all(
        abs(value) <= stationarity_threshold for value in arb_phase_sums.values()
    )
    if pole_pass and diagonal_pass:
        outcome = "FULL_REPRODUCTION"
    elif pole_pass:
        outcome = "POLE_ONLY_REPRODUCTION"
    else:
        outcome = "SOURCE_NORMALIZATION_MISMATCH"
    restricted_label = (
        "STRUCTURAL RESTRICTED-ACTION CANCELLATION"
        if phase_sum_pass and not diagonal_pass
        else "ABSENT"
    )
    record.update({
        "pole_pass": pole_pass,
        "diagonal_pass": diagonal_pass,
        "phase_sum_pass": phase_sum_pass,
        "outcome": outcome,
        "restricted_label": restricted_label,
        "arb_per_edge": arb_per_edge,
        "arb_diagonal": arb_diagonal,
        "arb_poles": arb_poles,
        "arb_phase_sums": arb_phase_sums,
    })
    print(f"[{parity}] outcome={outcome}")
    print(
        "       "
        f"max|pole|={float(max(abs(v) for v in arb_poles)):.12e}, "
        f"max|diagonal|={float(max(abs(v) for v in arb_diagonal)):.12e}, "
        f"||R35||={float(arb.sqrt(sum(abs(v)**2 for v in arb_per_edge))):.12e}, "
        f"phase-pair sums pass={phase_sum_pass}"
    )


parity_residual_multiset_error = max(
    abs(left-right)
    for left, right in zip(
        np.sort(records["even"]["per_edge"]),
        np.sort(records["odd"]["per_edge"]),
    )
)

outcomes = {record["outcome"] for record in records.values()}
if outcomes == {"FULL_REPRODUCTION"}:
    verdict = (
        "DERIVED EXTERNAL CONTROL: the published time-symmetric dust sandwich "
        "solves all 35 complete one-slab orbit equations in both parities."
    )
elif "SOURCE_NORMALIZATION_MISMATCH" not in outcomes:
    verdict = (
        "DERIVED PARTIAL EXTERNAL CONTROL: the published dust source balances "
        "all pole equations, but at least one complete one-slab diagonal "
        "equation remains nonstationary."
    )
else:
    verdict = (
        "DERIVED CONTROL MISMATCH: the published dust normalization does not "
        "balance every one-slab pole equation under the frozen conventions."
    )


def serializable_record(record):
    return {
        "gravitational_action_real": float(record["reduced_action"].real),
        "gravitational_action_imaginary": float(record["reduced_action"].imag),
        "total_action_real": float(record["total_action"].real),
        "total_action_imaginary": float(record["total_action"].imag),
        "per_edge_residuals": record["per_edge"].tolist(),
        "diagonal_residuals": record["diagonal"].tolist(),
        "pole_residuals": record["poles"].tolist(),
        "diagonal_maximum_absolute": float(np.max(np.abs(record["diagonal"]))),
        "diagonal_norm": float(np.linalg.norm(record["diagonal"])),
        "pole_maximum_absolute": float(np.max(np.abs(record["poles"]))),
        "pole_norm": float(np.linalg.norm(record["poles"])),
        "full_residual_norm": float(np.linalg.norm(record["per_edge"])),
        "phase_pairs": {
            pair: {
                "variables": item["variables"],
                "residuals": item["residuals"].tolist(),
                "sum": item["sum"],
                "spread": item["spread"],
            }
            for pair, item in record["phase_pairs"].items()
        },
        "pole_length_equation_lhs": record["pole_length_lhs"].tolist(),
        "pole_length_equation_rhs": record["pole_length_rhs"],
        "pole_balance_relative_error": record["pole_balance_error"],
        "full_reduced_action_error": record["action_error"],
        "full_reduced_gradient_error": record["gradient_error"],
        "full_reduced_old_gradient_error": record["old_gradient_error"],
        "imaginary_residual": record["imaginary_residual"],
        "minimum_absolute_gram_eigenvalue": record["minimum_gram"],
        "minimum_angle_argument_modulus": record["minimum_argument"],
        "direct_total_gradient_error": record["direct_error"],
        "direct_total_gradient_imaginary": record["direct_imaginary"],
        "binary64_direct_pass": record["binary64_direct_pass"],
        "arbitrary_precision_digits": 60,
        "arbitrary_precision_base_action_real": arb.nstr(
            arb.re(record["arb_base_action"]), 50
        ),
        "arbitrary_precision_base_action_imaginary": arb.nstr(
            arb.im(record["arb_base_action"]), 50
        ),
        "arbitrary_precision_base_action_error": record["arb_base_action_error"],
        "arbitrary_precision_total_gradient": [
            {
                "real": arb.nstr(arb.re(value), 50),
                "imaginary": arb.nstr(arb.im(value), 50),
            }
            for value in record["arb_total_gradient"]
        ],
        "arbitrary_precision_gravitational_gradient": [
            {
                "real": arb.nstr(arb.re(value), 50),
                "imaginary": arb.nstr(arb.im(value), 50),
            }
            for value in record["arb_gravitational_gradient"]
        ],
        "arbitrary_precision_dust_pole_derivative": [
            arb.nstr(value, 50) for value in record["arb_dust_pole_derivative"]
        ],
        "arbitrary_precision_per_edge_residuals": [
            {
                "real": arb.nstr(arb.re(value), 50),
                "imaginary": arb.nstr(arb.im(value), 50),
            }
            for value in record["arb_per_edge"]
        ],
        "arbitrary_precision_analytic_gradient_error": record["arb_analytic_error"],
        "arbitrary_precision_imaginary_derivative": record["arb_imaginary_derivative"],
        "arbitrary_precision_maximum_per_edge_residual": record[
            "arb_direct_per_edge_maximum"
        ],
        "arbitrary_precision_pole_cancellation_ratios": [
            arb.nstr(value, 50) for value in record["arb_cancellation_ratios"]
        ],
        "arbitrary_precision_phase_pair_sums": {
            pair: {
                "real": arb.nstr(arb.re(value), 50),
                "imaginary": arb.nstr(arb.im(value), 50),
            }
            for pair, value in record["arb_phase_sums"].items()
        },
        "pole_pass": record["pole_pass"],
        "diagonal_pass": record["diagonal_pass"],
        "phase_pair_sum_pass": record["phase_sum_pass"],
        "outcome": record["outcome"],
        "restricted_label": record["restricted_label"],
    }


payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "upstream_commit": UPSTREAM_COMMIT,
    "precision_correction_commit": PRECISION_CORRECTION_COMMIT,
    "source": "De Felice and Fabri, arXiv:gr-qc/0009093",
    "published_data": {
        "M_star": M_STAR,
        "zeta": ZETA,
        "R0": R0,
        "l0": L0,
        "l0_square": L0_SQUARE,
        "epsilon_3": EPSILON_3,
        "M": MASS,
        "tau": TAU,
        "tau_square": TAU_SQUARE,
        "slant_square": SLANT_SQUARE,
        "printed_control_errors": published_errors,
    },
    "stationarity_threshold": stationarity_threshold,
    "parities": {
        parity: serializable_record(record)
        for parity, record in records.items()
    },
    "parity_residual_multiset_error": parity_residual_multiset_error,
    "verdict": verdict,
    "labels": {
        "published_geometry_and_source": "EXTERNAL CONTROL",
        "binary64_total_action_differences": "RECORDED NUMERICAL FAILURE",
        "arbitrary_precision_total_action_differences": "CORRECTED CONTROL",
        "full_one_slab_stationarity": (
            "DERIVED" if outcomes == {"FULL_REPRODUCTION"} else "NOT ESTABLISHED"
        ),
        "multi_tick_reproduction": "NOT TESTED",
        "parameter_selection": "NOT CLAIMED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} implementation checks passed")
print(verdict)
raise SystemExit(0 if passed == tests else 1)
