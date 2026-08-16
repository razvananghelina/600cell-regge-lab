#!/usr/bin/env python3
"""Preregistered two-slab 600-cell dust gluing and momentum-sign control.

Protocol commit: 29dcfa5.  This verifier constructs the three-layer complex
directly, derives its orbit maps from vertices, and compares the direct
two-slab action with two independently evaluated one-slab factors.  It is a
structural control only and performs no root search.
"""

from collections import Counter
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
OUTPUT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
PROTOCOL_COMMIT = "29dcfa5"
PRIOR_ART_COMMIT = "8c45290"
FRAMING_CORRECTION_COMMIT = "620461d"
SCHEDULE_CORRECTION_COMMIT = "6c4a377"
DPS = 100
ACTION_RELATIVE_TOLERANCE = arb.mpf("5e-8")
IMAGINARY_TOLERANCE = arb.mpf("1e-70")
DERIVATIVE_FLOOR = arb.mpf("1e-60")
DERIVATIVE_GATE_FACTOR = arb.mpf(10)
NONZERO_GATE_FACTOR = arb.mpf(100)
AUDIT_LOG_STEP = arb.mpf("1e-6")
DERIVATIVE_STEPS = {
    "operational": (arb.mpf("1e-20"), arb.mpf("1e-15")),
    "validation": (arb.mpf("3e-20"), arb.mpf("3e-15")),
}
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def relative_error(left, right):
    return abs(left-right)/max(arb.mpf(1), abs(left), abs(right))


# Import only the certified geometric/action core.  The new action and all
# boundary differences below are reconstructed here, not read from its JSON.
spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_two_slab", HERE / "verify_gravity_global_regge_orbits.py"
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise

check(
    "the imported one-slab geometric/action core retains all 43 certificates",
    gro.tests == gro.passed == 43,
)


def orbit_sort_key(orbit, phase):
    representative = min(orbit)
    logical = tuple(vertex % 120 for vertex in representative)
    phase_pair = tuple(sorted(phase[vertex] for vertex in logical))
    return phase_pair, tuple(sorted(orbit))


def augment_boundary_orbits(base):
    old_orbits = tuple(sorted(
        gro.orbit_partition(base["old_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    final_orbits = tuple(sorted(
        gro.orbit_partition(base["new_edges"], base["stabilizer"]),
        key=lambda orbit: orbit_sort_key(orbit, base["phase"]),
    ))
    old_lookup = {
        edge: index for index, orbit in enumerate(old_orbits) for edge in orbit
    }
    final_lookup = {
        edge: index for index, orbit in enumerate(final_orbits) for edge in orbit
    }
    return {
        **base,
        "old_orbits": old_orbits,
        "final_orbits": final_orbits,
        "old_lookup": old_lookup,
        "final_lookup": final_lookup,
    }


models = {
    parity: augment_boundary_orbits(model) for parity, model in gro.models.items()
}


def layered_image(action, vertex):
    return int(action[vertex % 120]) + 120*(vertex//120)


def layered_orbit_partition(items, stabilizer):
    item_set = set(items)
    unseen = set(items)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = frozenset(
            tuple(sorted(layered_image(action, vertex) for vertex in seed))
            for action in stabilizer
        )
        if not orbit <= item_set:
            raise RuntimeError("three-layer stabilizer left the complex")
        orbits.append(orbit)
        unseen -= orbit
    return tuple(orbits)


def shifted(item, amount):
    return tuple(vertex+amount for vertex in item)


def derive_maps(model):
    final_orbit_by_set = {
        frozenset(orbit): index for index, orbit in enumerate(model["final_orbits"])
    }
    old_to_final = []
    for orbit in model["old_orbits"]:
        image = frozenset(shifted(edge, 120) for edge in orbit)
        old_to_final.append(final_orbit_by_set[image])
    if sorted(old_to_final) != list(range(30)):
        raise RuntimeError("consecutive old-to-final orbit map is not bijective")
    final_to_old = [None]*30
    for old_index, final_index in enumerate(old_to_final):
        final_to_old[final_index] = old_index

    reversal_actions = []
    induced_maps = []
    ordering = model["ordering"]
    for action in gro.h4_actions:
        if not all(
            frozenset(int(action[vertex]) for vertex in ordering[index])
            == ordering[4-index]
            for index in range(5)
        ):
            continue
        reversal_actions.append(action)

        def reverse_vertex(vertex):
            return int(action[vertex % 120]) + 120*(1-vertex//120)

        reversed_slab = frozenset(
            tuple(sorted(reverse_vertex(vertex) for vertex in simplex))
            for simplex in model["slab"]
        )
        if reversed_slab != model["slab"]:
            raise RuntimeError("phase-reversing action does not reverse the slab")
        final_to_old_reversal = []
        for orbit in model["final_orbits"]:
            representative = min(orbit)
            image = tuple(sorted(reverse_vertex(vertex) for vertex in representative))
            final_to_old_reversal.append(model["old_lookup"][image])
        induced_maps.append(tuple(final_to_old_reversal))
    unique_reversal_maps = sorted(set(induced_maps))
    if len(unique_reversal_maps) != 1:
        raise RuntimeError(
            f"time reversal induced {len(unique_reversal_maps)} quotient maps"
        )
    return {
        "old_to_final": tuple(old_to_final),
        "final_to_old": tuple(final_to_old),
        "reversal_action_count": len(reversal_actions),
        "reversal_unique_map_count": len(unique_reversal_maps),
        "reversal_final_to_old": unique_reversal_maps[0],
    }


def construct_combined(model):
    first_slab = frozenset(model["slab"])
    second_slab = frozenset(shifted(simplex, 120) for simplex in model["slab"])
    slab = first_slab | second_slab
    if first_slab & second_slab:
        raise RuntimeError("the two slabs share a four-simplex")
    edges = frozenset(
        tuple(sorted(edge))
        for simplex in slab for edge in combinations(simplex, 2)
    )
    triangles = frozenset(
        tuple(sorted(triangle))
        for simplex in slab for triangle in combinations(simplex, 3)
    )
    facets = Counter(
        tuple(sorted(facet))
        for simplex in slab for facet in combinations(simplex, 4)
    )
    triangle_orbits = tuple(sorted(
        layered_orbit_partition(triangles, model["stabilizer"]),
        key=gro.canonical_key,
    ))
    simplex_orbits = tuple(sorted(
        layered_orbit_partition(slab, model["stabilizer"]),
        key=gro.canonical_key,
    ))
    triangle_lookup = {
        triangle: index
        for index, orbit in enumerate(triangle_orbits) for triangle in orbit
    }

    edge_data = {}

    def insert(edge, value):
        if edge in edge_data and edge_data[edge] != value:
            raise RuntimeError(f"ambiguous combined edge assignment: {edge}")
        edge_data[edge] = value

    for edge, index in model["old_lookup"].items():
        insert(edge, ("q0", index, 1))
    for edge, index in model["final_lookup"].items():
        insert(edge, ("q1", index, 1))
    for edge, index in model["final_lookup"].items():
        insert(shifted(edge, 120), ("q2", index, 1))
    for edge, index in model["edge_to_variable"].items():
        insert(edge, ("x1", index, int(model["edge_jacobian"][edge])))
        insert(
            shifted(edge, 120),
            ("x2", index, int(model["edge_jacobian"][edge])),
        )
    if set(edge_data) != set(edges):
        missing = sorted(set(edges)-set(edge_data))
        extra = sorted(set(edge_data)-set(edges))
        raise RuntimeError(f"combined edge cover mismatch: missing={missing[:3]}, extra={extra[:3]}")

    outer_triangles = frozenset(
        triangle for triangle in triangles
        if len({vertex//120 for vertex in triangle}) == 1
        and next(iter({vertex//120 for vertex in triangle})) in (0, 2)
    )
    shared_facets = frozenset(shifted(tetrahedron, 120) for tetrahedron in gro.tetrahedra)
    return {
        "slab": slab,
        "first_slab": first_slab,
        "second_slab": second_slab,
        "edges": edges,
        "triangles": triangles,
        "facets": facets,
        "shared_facets": shared_facets,
        "triangle_orbits": triangle_orbits,
        "triangle_lookup": triangle_lookup,
        "simplex_orbits": simplex_orbits,
        "edge_data": edge_data,
        "outer_triangles": outer_triangles,
    }


maps = {parity: derive_maps(model) for parity, model in models.items()}
combined_models = {
    parity: construct_combined(model) for parity, model in models.items()
}


for parity in ("even", "odd"):
    model = models[parity]
    mapping = maps[parity]
    combined = combined_models[parity]
    check(
        f"{parity}: incidence derives a bijective consecutive orbit map",
        sorted(mapping["old_to_final"]) == list(range(30))
        and sorted(mapping["final_to_old"]) == list(range(30)),
    )
    check(
        f"{parity}: all 24 phase-reversing automorphisms induce one quotient map",
        mapping["reversal_action_count"] == 24
        and mapping["reversal_unique_map_count"] == 1
        and sorted(mapping["reversal_final_to_old"]) == list(range(30)),
    )
    shared_multiplicities = Counter(
        combined["facets"][facet] for facet in combined["shared_facets"]
    )
    check(
        f"{parity}: the direct carrier is a 360-vertex, 4800-simplex gluing",
        len({vertex for simplex in combined["slab"] for vertex in simplex}) == 360
        and len(combined["first_slab"]) == len(combined["second_slab"]) == 2400
        and len(combined["slab"]) == 4800
        and len(combined["shared_facets"]) == 600
        and shared_multiplicities == Counter({2: 600})
        and set(combined["edge_data"]) == set(combined["edges"]),
        f"triangles={len(combined['triangles'])}, triangle orbits={len(combined['triangle_orbits'])}, "
        f"simplex orbits={len(combined['simplex_orbits'])}",
    )
    check(
        f"{parity}: every direct simplex and triangle orbit has size 24",
        Counter(map(len, combined["simplex_orbits"])) == Counter({24: 200})
        and sum(map(len, combined["triangle_orbits"])) == len(combined["triangles"])
        and set(map(len, combined["triangle_orbits"])) == {24},
    )


# Published external-control values reconstructed from their displayed formulas.
arb.mp.dps = DPS
ARB_M_STAR = arb.mpf(10)
ARB_ZETA = (arb.pi**2*arb.sqrt(2)/50)**(arb.mpf(1)/3)
ARB_R0 = 4*ARB_M_STAR/(3*arb.pi)
ARB_L0 = ARB_ZETA*ARB_R0
ARB_L0_SQUARE = ARB_L0**2
ARB_EPSILON_3 = 2*arb.pi-5*arb.acos(arb.mpf(1)/3)
ARB_MASS = (90/arb.pi)*ARB_EPSILON_3*ARB_L0
ARB_TAU = arb.mpf("0.0102")
ARB_TAU_SQUARE = ARB_TAU**2
ARB_SLANT_SQUARE = ARB_L0_SQUARE-ARB_TAU_SQUARE
ARB_BASE_Q = tuple(ARB_L0_SQUARE for _ in range(30))
ARB_BASE_X = tuple(
    [ARB_SLANT_SQUARE for _ in range(30)]
    + [ARB_TAU_SQUARE for _ in range(5)]
)
ARB_I = arb.mpc(0, 1)


def arb_log_minus(value):
    scale = max(arb.mpf(1), abs(value))
    if abs(arb.im(value)) < arb.mpf("1e-80")*scale:
        real = arb.re(value)
        if real < 0:
            return arb.log(-real)-ARB_I*arb.pi
        return arb.log(real)
    return arb.log(value)


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
        hinge_volume_square = arb_signed_volume_square(squared, hinge_vertices)
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
        angles[hinge_vertices] = -ARB_I*arb_log_minus(cosine+ARB_I*sine)
    return angles


def arb_triangle_area_square(values):
    x, y, z = values
    return (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16


def one_edge_square(model, edge, x, q_old, q_final):
    edge = tuple(sorted(edge))
    if edge in model["edge_to_variable"]:
        index = model["edge_to_variable"][edge]
        return int(model["edge_jacobian"][edge])*x[index]
    if edge in model["old_lookup"]:
        return q_old[model["old_lookup"][edge]]
    if edge in model["final_lookup"]:
        return q_final[model["final_lookup"][edge]]
    raise ValueError(f"edge absent from one-slab evaluator: {edge}")


def one_simplex_squared(model, simplex, x, q_old, q_final):
    squared = [[arb.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value = one_edge_square(
            model, (simplex[left], simplex[right]), x, q_old, q_final
        )
        squared[left][right] = squared[right][left] = value
    return squared


def arb_one_action(model, x, q_old, q_final):
    curvature = [
        arb.pi if min(orbit) in model["boundary_triangles"] else 2*arb.pi
        for orbit in model["triangle_orbits"]
    ]
    for simplex_orbit in model["simplex_orbits"]:
        simplex = min(simplex_orbit)
        squared = one_simplex_squared(model, simplex, x, q_old, q_final)
        for local_hinge, angle in arb_angle_data(squared).items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            curvature[model["triangle_to_orbit"][triangle]] += angle
    action_sum = arb.mpc(0)
    for index, orbit in enumerate(model["triangle_orbits"]):
        triangle = min(orbit)
        values = tuple(
            one_edge_square(model, edge, x, q_old, q_final)
            for edge in combinations(triangle, 2)
        )
        area = arb.sqrt(arb.mpc(arb_triangle_area_square(values)))
        action_sum += len(orbit)*area*curvature[index]
    gravitational = -ARB_I*action_sum
    dust = -(8*arb.pi*ARB_MASS/5)*sum(arb.sqrt(value) for value in x[30:35])
    return gravitational+dust


def direct_edge_square(combined, edge, x1, x2, q0, q1, q2):
    kind, index, sign = combined["edge_data"][tuple(sorted(edge))]
    values = {"x1": x1, "x2": x2, "q0": q0, "q1": q1, "q2": q2}[kind]
    return sign*values[index]


def direct_simplex_squared(combined, simplex, x1, x2, q0, q1, q2):
    squared = [[arb.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value = direct_edge_square(
            combined, (simplex[left], simplex[right]), x1, x2, q0, q1, q2
        )
        squared[left][right] = squared[right][left] = value
    return squared


def arb_direct_action(combined, x1, x2, q0, q1, q2):
    curvature = [
        arb.pi if min(orbit) in combined["outer_triangles"] else 2*arb.pi
        for orbit in combined["triangle_orbits"]
    ]
    for simplex_orbit in combined["simplex_orbits"]:
        simplex = min(simplex_orbit)
        squared = direct_simplex_squared(combined, simplex, x1, x2, q0, q1, q2)
        for local_hinge, angle in arb_angle_data(squared).items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            curvature[combined["triangle_lookup"][triangle]] += angle
    action_sum = arb.mpc(0)
    for index, orbit in enumerate(combined["triangle_orbits"]):
        triangle = min(orbit)
        values = tuple(
            direct_edge_square(combined, edge, x1, x2, q0, q1, q2)
            for edge in combinations(triangle, 2)
        )
        area = arb.sqrt(arb.mpc(arb_triangle_area_square(values)))
        action_sum += len(orbit)*area*curvature[index]
    gravitational = -ARB_I*action_sum
    dust = -(8*arb.pi*ARB_MASS/5)*(
        sum(arb.sqrt(value) for value in x1[30:35])
        + sum(arb.sqrt(value) for value in x2[30:35])
    )
    return gravitational+dust


def perturb(values, index, delta):
    result = list(values)
    result[index] *= arb.exp(delta)
    return tuple(result)


def old_q_from_shared(q_shared, mapping):
    return tuple(q_shared[mapping["old_to_final"][index]] for index in range(30))


def pack_complex(value):
    return arb.nstr(arb.re(value), DPS+5), arb.nstr(arb.im(value), DPS+5)


def unpack_complex(value):
    return arb.mpc(arb.mpf(value[0]), arb.mpf(value[1]))


_WORKER_MODEL = None
_WORKER_COMBINED = None
_WORKER_MAP = None


def initialize_worker(model, combined, mapping):
    global _WORKER_MODEL, _WORKER_COMBINED, _WORKER_MAP
    arb.mp.dps = DPS
    _WORKER_MODEL = model
    _WORKER_COMBINED = combined
    _WORKER_MAP = mapping


def action_worker(task):
    kind = task[0]
    if kind == "one":
        _, boundary, index, delta_text = task
        delta = arb.mpf(delta_text)
        q_old = ARB_BASE_Q
        q_final = ARB_BASE_Q
        if boundary == "old":
            q_old = perturb(q_old, index, delta)
        else:
            q_final = perturb(q_final, index, delta)
        return pack_complex(arb_one_action(
            _WORKER_MODEL, ARB_BASE_X, q_old, q_final
        ))
    if kind == "direct":
        _, index, delta_text = task
        q1 = perturb(ARB_BASE_Q, index, arb.mpf(delta_text))
        return pack_complex(arb_direct_action(
            _WORKER_COMBINED,
            ARB_BASE_X, ARB_BASE_X,
            ARB_BASE_Q, q1, ARB_BASE_Q,
        ))
    if kind == "glue":
        _, index, delta_text = task
        q1 = ARB_BASE_Q if index < 0 else perturb(
            ARB_BASE_Q, index, arb.mpf(delta_text)
        )
        direct = arb_direct_action(
            _WORKER_COMBINED,
            ARB_BASE_X, ARB_BASE_X,
            ARB_BASE_Q, q1, ARB_BASE_Q,
        )
        factor_sum = (
            arb_one_action(_WORKER_MODEL, ARB_BASE_X, ARB_BASE_Q, q1)
            + arb_one_action(
                _WORKER_MODEL,
                ARB_BASE_X,
                old_q_from_shared(q1, _WORKER_MAP),
                ARB_BASE_Q,
            )
        )
        return pack_complex(direct), pack_complex(factor_sum)
    raise ValueError(f"unknown action task {kind}")


def float_one_edge_square(model, edge, q_old, q_final):
    edge = tuple(sorted(edge))
    if edge in model["edge_to_variable"]:
        index = model["edge_to_variable"][edge]
        return int(model["edge_jacobian"][edge])*float(ARB_BASE_X[index])
    if edge in model["old_lookup"]:
        return q_old[model["old_lookup"][edge]]
    return q_final[model["final_lookup"][edge]]


def float_direct_edge_square(combined, edge, q1):
    kind, index, sign = combined["edge_data"][tuple(sorted(edge))]
    if kind in ("x1", "x2"):
        return sign*float(ARB_BASE_X[index])
    if kind == "q1":
        return q1[index]
    return float(ARB_BASE_Q[index])


def branch_signature_one(model, boundary=None, index=None, delta=0.0):
    q_old = np.full(30, float(ARB_L0_SQUARE))
    q_final = np.full(30, float(ARB_L0_SQUARE))
    if boundary == "old":
        q_old[index] *= math.exp(delta)
    elif boundary == "final":
        q_final[index] *= math.exp(delta)
    negative = Counter()
    minimum_gram = math.inf
    minimum_argument = math.inf
    for orbit in model["simplex_orbits"]:
        simplex = min(orbit)
        squared = np.zeros((5, 5))
        for left, right in combinations(range(5), 2):
            value = float_one_edge_square(
                model, (simplex[left], simplex[right]), q_old, q_final
            )
            squared[left, right] = squared[right, left] = value
        _, arguments, eigenvalues = gro.angle_data(squared)
        negative[int(np.sum(eigenvalues < -1e-10))] += len(orbit)
        minimum_gram = min(minimum_gram, float(np.min(np.abs(eigenvalues))))
        minimum_argument = min(minimum_argument, *(abs(v) for v in arguments.values()))
    return negative, minimum_gram, minimum_argument


def branch_signature_direct(combined, index=None, delta=0.0):
    q1 = np.full(30, float(ARB_L0_SQUARE))
    if index is not None:
        q1[index] *= math.exp(delta)
    negative = Counter()
    minimum_gram = math.inf
    minimum_argument = math.inf
    for orbit in combined["simplex_orbits"]:
        simplex = min(orbit)
        squared = np.zeros((5, 5))
        for left, right in combinations(range(5), 2):
            value = float_direct_edge_square(
                combined, (simplex[left], simplex[right]), q1
            )
            squared[left, right] = squared[right, left] = value
        _, arguments, eigenvalues = gro.angle_data(squared)
        negative[int(np.sum(eigenvalues < -1e-10))] += len(orbit)
        minimum_gram = min(minimum_gram, float(np.min(np.abs(eigenvalues))))
        minimum_argument = min(minimum_argument, *(abs(v) for v in arguments.values()))
    return negative, minimum_gram, minimum_argument


def vector_norm(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def derivative_vectors(action_values, kind, size):
    result = {}
    for name, (primary_h, shadow_h) in DERIVATIVE_STEPS.items():
        rows = []
        for h in (primary_h, shadow_h):
            row = []
            for index in range(size):
                plus = action_values[(kind, index, arb.nstr(h, 20))]
                minus = action_values[(kind, index, arb.nstr(-h, 20))]
                row.append((plus-minus)/(48*h))
            rows.append(tuple(row))
        result[name] = {
            "primary": rows[0],
            "shadow": rows[1],
            "proxy": tuple(a-b for a, b in zip(rows[0], rows[1])),
        }
    result["agreement"] = tuple(
        result["operational"]["primary"][index]
        - result["validation"]["primary"][index]
        for index in range(size)
    )
    result["uncertainty"] = tuple(
        abs(result["operational"]["proxy"][index])
        + abs(result["validation"]["proxy"][index])
        + abs(result["agreement"][index])
        + DERIVATIVE_FLOOR
        for index in range(size)
    )
    result["calibration_pass"] = all(
        abs(result["agreement"][index])
        <= DERIVATIVE_GATE_FACTOR*(
            abs(result["operational"]["proxy"][index])
            + abs(result["validation"]["proxy"][index])
            + DERIVATIVE_FLOOR
        )
        for index in range(size)
    )
    return result


def complex_record(value):
    return {
        "real": arb.nstr(arb.re(value), 60),
        "imaginary": arb.nstr(arb.im(value), 20),
    }


def vector_record(values):
    return [complex_record(value) for value in values]


print("="*78)
print("TWO-SLAB 600-CELL DUST GLUING CONTROL")
print("="*78)

records = {}
overall_outcome = "TWO_SLAB_GLUING_CONTROL_PASSED"
fork_context = mp.get_context("fork")

for parity in ("even", "odd"):
    print(f"[{parity}] auditing every frozen branch point")
    model = models[parity]
    combined = combined_models[parity]
    mapping = maps[parity]

    branch_minimum_gram = math.inf
    branch_minimum_argument = math.inf
    branch_pass = True
    branch_evaluations = 0

    # Base and every +/-1e-6 action-gluing audit point.  The direct signature
    # must equal the two factor signatures after the incidence-derived map.
    gluing_branch_match = True
    base_one = branch_signature_one(model)
    base_direct = branch_signature_direct(combined)
    branch_evaluations += 2
    branch_minimum_gram = min(branch_minimum_gram, base_one[1], base_direct[1])
    branch_minimum_argument = min(
        branch_minimum_argument, base_one[2], base_direct[2]
    )
    branch_pass &= base_one[0] == Counter({1: 2400})
    branch_pass &= base_direct[0] == Counter({1: 4800})
    gluing_branch_match &= base_direct[0] == base_one[0]+base_one[0]
    for index in range(30):
        old_index = mapping["final_to_old"][index]
        for delta in (float(AUDIT_LOG_STEP), -float(AUDIT_LOG_STEP)):
            first = branch_signature_one(model, "final", index, delta)
            second = branch_signature_one(model, "old", old_index, delta)
            direct = branch_signature_direct(combined, index, delta)
            branch_evaluations += 3
            branch_minimum_gram = min(
                branch_minimum_gram, first[1], second[1], direct[1]
            )
            branch_minimum_argument = min(
                branch_minimum_argument, first[2], second[2], direct[2]
            )
            branch_pass &= (
                first[0] == second[0] == Counter({1: 2400})
                and direct[0] == Counter({1: 4800})
            )
            gluing_branch_match &= direct[0] == first[0]+second[0]

    # Audit every tiny finite-difference point as frozen, rather than assuming
    # that it inherits the branch of the larger +/-1e-6 envelope.
    derivative_deltas = sorted({
        float(sign*h)
        for pair in DERIVATIVE_STEPS.values() for h in pair for sign in (-1, 1)
    })
    for boundary in ("old", "final"):
        for index in range(30):
            for delta in derivative_deltas:
                signature = branch_signature_one(model, boundary, index, delta)
                branch_evaluations += 1
                branch_minimum_gram = min(branch_minimum_gram, signature[1])
                branch_minimum_argument = min(branch_minimum_argument, signature[2])
                branch_pass &= signature[0] == Counter({1: 2400})
    for index in range(30):
        for delta in derivative_deltas:
            signature = branch_signature_direct(combined, index, delta)
            branch_evaluations += 1
            branch_minimum_gram = min(branch_minimum_gram, signature[1])
            branch_minimum_argument = min(branch_minimum_argument, signature[2])
            branch_pass &= signature[0] == Counter({1: 4800})

    branch_pass &= branch_minimum_gram > 1e-8 and branch_minimum_argument > 1e-6
    check(
        f"{parity}: all {branch_evaluations} frozen evaluations retain one timelike direction",
        branch_pass,
        f"min Gram={branch_minimum_gram:.3e}, min argument={branch_minimum_argument:.3e}",
    )
    check(
        f"{parity}: direct and factor branch counts agree at all 61 gluing points",
        gluing_branch_match,
    )
    if not (branch_pass and gluing_branch_match):
        overall_outcome = "TWO_SLAB_DERIVATIVE_CONTROL_FAILED"

    tasks = [("glue", -1, "0")]
    tasks.extend(
        ("glue", index, arb.nstr(sign*AUDIT_LOG_STEP, 20))
        for index in range(30) for sign in (1, -1)
    )
    for boundary in ("old", "final"):
        for index in range(30):
            for pair in DERIVATIVE_STEPS.values():
                for h in pair:
                    for sign in (1, -1):
                        tasks.append((
                            "one", boundary, index, arb.nstr(sign*h, 20)
                        ))
    for index in range(30):
        for pair in DERIVATIVE_STEPS.values():
            for h in pair:
                for sign in (1, -1):
                    tasks.append(("direct", index, arb.nstr(sign*h, 20)))

    print(f"[{parity}] evaluating {len(tasks)} arbitrary-precision action tasks")
    with fork_context.Pool(
        processes=8,
        initializer=initialize_worker,
        initargs=(model, combined, mapping),
    ) as pool:
        raw_results = pool.map(action_worker, tasks, chunksize=1)

    gluing_errors = []
    maximum_imaginary = arb.mpf(0)
    action_values = {}
    for task, raw in zip(tasks, raw_results):
        if task[0] == "glue":
            direct = unpack_complex(raw[0])
            factor_sum = unpack_complex(raw[1])
            gluing_errors.append(relative_error(direct, factor_sum))
            maximum_imaginary = max(
                maximum_imaginary,
                abs(arb.im(direct)),
                abs(arb.im(factor_sum)),
            )
        else:
            value = unpack_complex(raw)
            boundary_kind = task[1] if task[0] == "one" else "shared"
            index = task[2] if task[0] == "one" else task[1]
            delta_text = task[3] if task[0] == "one" else task[2]
            action_values[(boundary_kind, index, delta_text)] = value
            maximum_imaginary = max(maximum_imaginary, abs(arb.im(value)))

    action_gluing_pass = (
        len(gluing_errors) == 61
        and max(gluing_errors) < ACTION_RELATIVE_TOLERANCE
    )
    check(
        f"{parity}: direct and summed actions agree at the base and all 60 orbit audits",
        action_gluing_pass,
        f"max relative error={float(max(gluing_errors)):.3e}",
    )
    if not action_gluing_pass:
        overall_outcome = "TWO_SLAB_ACTION_GLUING_FAILED"

    old_derivatives = derivative_vectors(action_values, "old", 30)
    final_derivatives = derivative_vectors(action_values, "final", 30)
    shared_derivatives = derivative_vectors(action_values, "shared", 30)
    derivative_calibration_pass = (
        old_derivatives["calibration_pass"]
        and final_derivatives["calibration_pass"]
        and shared_derivatives["calibration_pass"]
        and maximum_imaginary < IMAGINARY_TOLERANCE
    )
    check(
        f"{parity}: all 90 derivatives pass the frozen operational/validation calibration",
        derivative_calibration_pass,
        f"max imaginary={arb.nstr(maximum_imaginary, 8)}",
    )
    if not derivative_calibration_pass:
        overall_outcome = "TWO_SLAB_DERIVATIVE_CONTROL_FAILED"

    p_pre = tuple(-value for value in old_derivatives["operational"]["primary"])
    p_post = final_derivatives["operational"]["primary"]
    direct_shared = shared_derivatives["operational"]["primary"]

    reversal_post = [arb.mpc(0) for _ in range(30)]
    reversal_uncertainty = [arb.mpf(0) for _ in range(30)]
    for final_index, old_index in enumerate(mapping["reversal_final_to_old"]):
        reversal_post[old_index] = p_post[final_index]
        reversal_uncertainty[old_index] = final_derivatives["uncertainty"][final_index]
    time_reversal_residual = tuple(
        p_pre[index]+reversal_post[index] for index in range(30)
    )
    time_reversal_pass = all(
        abs(time_reversal_residual[index])
        <= DERIVATIVE_GATE_FACTOR*(
            old_derivatives["uncertainty"][index]
            + reversal_uncertainty[index]
        )
        for index in range(30)
    )
    check(
        f"{parity}: the vertex-derived reversal map gives p_pre + R p_post = 0 orbitwise",
        time_reversal_pass,
        f"max residual={float(max(map(abs, time_reversal_residual))):.3e}",
    )

    mapped_pre = [arb.mpc(0) for _ in range(30)]
    mapped_pre_uncertainty = [arb.mpf(0) for _ in range(30)]
    for old_index, final_index in enumerate(mapping["old_to_final"]):
        mapped_pre[final_index] = p_pre[old_index]
        mapped_pre_uncertainty[final_index] = old_derivatives["uncertainty"][old_index]
    cusp = tuple(
        p_post[index]-mapped_pre[index] for index in range(30)
    )
    shared_identity_residual = tuple(
        direct_shared[index]-cusp[index] for index in range(30)
    )
    shared_error = tuple(
        shared_derivatives["uncertainty"][index]
        + final_derivatives["uncertainty"][index]
        + mapped_pre_uncertainty[index]
        for index in range(30)
    )
    shared_identity_pass = all(
        abs(shared_identity_residual[index])
        <= DERIVATIVE_GATE_FACTOR*shared_error[index]
        for index in range(30)
    )
    cusp_uncertainty = tuple(
        final_derivatives["uncertainty"][index]
        + mapped_pre_uncertainty[index]
        for index in range(30)
    )
    cusp_norm = vector_norm(cusp)
    direct_shared_norm = vector_norm(direct_shared)
    cusp_error_norm = vector_norm(cusp_uncertainty)
    shared_error_norm = vector_norm(shared_derivatives["uncertainty"])
    nonzero_pass = (
        cusp_norm > NONZERO_GATE_FACTOR*max(cusp_error_norm, DERIVATIVE_FLOOR)
        and direct_shared_norm
        > NONZERO_GATE_FACTOR*max(shared_error_norm, DERIVATIVE_FLOOR)
    )
    check(
        f"{parity}: the direct shared derivative equals p_post - P p_pre orbitwise",
        shared_identity_pass,
        f"max residual={float(max(map(abs, shared_identity_residual))):.3e}",
    )
    check(
        f"{parity}: the repeated sandwich has the preregistered nonzero momentum cusp",
        nonzero_pass,
        f"||cusp||={float(cusp_norm):.12e}, ||direct||={float(direct_shared_norm):.12e}, "
        f"error={float(max(cusp_error_norm, shared_error_norm)):.3e}",
    )
    momentum_pass = time_reversal_pass and shared_identity_pass and nonzero_pass
    if not momentum_pass:
        overall_outcome = "TWO_SLAB_MOMENTUM_SIGN_CONTROL_FAILED"

    records[parity] = {
        "geometry": {
            "vertices": 360,
            "four_simplices": len(combined["slab"]),
            "shared_tetrahedra": len(combined["shared_facets"]),
            "triangles": len(combined["triangles"]),
            "triangle_orbits": len(combined["triangle_orbits"]),
            "simplex_orbits": len(combined["simplex_orbits"]),
            "old_to_final_orbit_map": list(mapping["old_to_final"]),
            "reversal_final_to_old_orbit_map": list(mapping["reversal_final_to_old"]),
            "reversal_action_count": mapping["reversal_action_count"],
            "unique_reversal_quotient_maps": mapping["reversal_unique_map_count"],
        },
        "branch": {
            "evaluations": branch_evaluations,
            "minimum_absolute_gram_eigenvalue": branch_minimum_gram,
            "minimum_angle_argument_modulus": branch_minimum_argument,
            "pass": branch_pass and gluing_branch_match,
        },
        "action_gluing": {
            "points": len(gluing_errors),
            "maximum_relative_error": float(max(gluing_errors)),
            "pass": action_gluing_pass,
        },
        "derivative_calibration": {
            "digits": DPS,
            "steps": {
                name: [arb.nstr(value, 10) for value in pair]
                for name, pair in DERIVATIVE_STEPS.items()
            },
            "maximum_imaginary_contamination": arb.nstr(maximum_imaginary, 20),
            "old_pass": old_derivatives["calibration_pass"],
            "final_pass": final_derivatives["calibration_pass"],
            "shared_pass": shared_derivatives["calibration_pass"],
            "pass": derivative_calibration_pass,
        },
        "momenta": {
            "pre": vector_record(p_pre),
            "post": vector_record(p_post),
            "direct_shared": vector_record(direct_shared),
            "predicted_cusp": vector_record(cusp),
            "time_reversal_residual": vector_record(time_reversal_residual),
            "shared_identity_residual": vector_record(shared_identity_residual),
            "cusp_norm": arb.nstr(cusp_norm, 60),
            "direct_shared_norm": arb.nstr(direct_shared_norm, 60),
            "cusp_uncertainty_norm": arb.nstr(cusp_error_norm, 20),
            "direct_shared_uncertainty_norm": arb.nstr(shared_error_norm, 20),
            "time_reversal_pass": time_reversal_pass,
            "shared_identity_pass": shared_identity_pass,
            "nonzero_cusp_pass": nonzero_pass,
        },
        "pass": (
            branch_pass and gluing_branch_match and action_gluing_pass
            and derivative_calibration_pass and momentum_pass
        ),
    }


if not all(record["pass"] for record in records.values()):
    if any(not record["branch"]["pass"] for record in records.values()):
        overall_outcome = "TWO_SLAB_DERIVATIVE_CONTROL_FAILED"
    elif any(not record["action_gluing"]["pass"] for record in records.values()):
        overall_outcome = "TWO_SLAB_ACTION_GLUING_FAILED"
    elif any(not record["derivative_calibration"]["pass"] for record in records.values()):
        overall_outcome = "TWO_SLAB_DERIVATIVE_CONTROL_FAILED"
    else:
        overall_outcome = "TWO_SLAB_MOMENTUM_SIGN_CONTROL_FAILED"

verdict = (
    "DERIVED CONTROL: direct two-slab action gluing, canonical orbit maps, "
    "boundary-momentum signs and the nonzero repeated-sandwich cusp pass in "
    "both schedule parities.  No evolving next frame has been solved."
    if overall_outcome == "TWO_SLAB_GLUING_CONTROL_PASSED"
    else "DERIVED CONTROL FAILURE: the preregistered two-slab gluing chain "
         f"stopped at {overall_outcome}."
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "prior_art_commit": PRIOR_ART_COMMIT,
    "framing_correction_commit": FRAMING_CORRECTION_COMMIT,
    "schedule_correction_commit": SCHEDULE_CORRECTION_COMMIT,
    "attempts": ["even repeated schedule", "odd repeated schedule"],
    "attempt_count": 2,
    "outcome": overall_outcome,
    "labels": {
        "two_slab_gluing": "DERIVED CONTROL" if all(
            record["pass"] for record in records.values()
        ) else "DERIVED CONTROL FAILURE",
        "next_frame": "NOT SOLVED",
        "canonical_evolution": "OPEN",
        "physical_time_scale": "OPEN",
    },
    "parities": records,
    "tests": tests,
    "passed": passed,
    "verdict": verdict,
}
OUTPUT.write_text(json.dumps(payload, indent=2)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {overall_outcome}")
print(verdict)
raise SystemExit(0 if passed == tests and all(r["pass"] for r in records.values()) else 1)
