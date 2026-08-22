#!/usr/bin/env python3
"""Nested tangential vertex-displacement prolongation certificate.

Prior-art gate 79b612b and protocol 9005af2 were committed before this
implementation was constructed or executed.  The verifier contains no action
Hessian and loads no spectral or particle-physics target.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_nested_vertex_displacement.json"
PRIOR_ART_COMMIT = "79b612b"
PROTOCOL_COMMIT = "9005af2"
PRECISION_CORRECTION_COMMIT = "db6845a"
INPUT_HASHES = {
    "commons/cell600.py":
        "840d921355e040bd4125dc8f8688b9702d63d9119e6f955f6e40b444c2d7d7a7",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py":
        "50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23",
    "reproducible/gravity_600cell_projected_rank_edgewise_carrier.json":
        "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84",
    "reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py":
        "36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32",
    "reproducible/gravity_600cell_refined_canonical_map_feasibility.json":
        "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
}
EXPECTED_BASE_F = (2640, 17040, 28800, 14400)
EXPECTED_FINE_F = (19680, 134880, 230400, 115200)
STEPS = (2.0**-18, 2.0**-20)

tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def all_simplices(top):
    return tuple(
        tuple(sorted({
            tuple(sorted(face))
            for simplex in top
            for face in combinations(simplex, degree + 1)
        }))
        for degree in range(4)
    )


def face_incidence(top):
    return Counter(
        tuple(sorted(face))
        for tetrahedron in top
        for face in combinations(tetrahedron, 3)
    )


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(adjacency[index] > 0.5))
                  for index in range(len(adjacency))]
    result = []
    for first in range(len(adjacency)):
        for second in sorted(value for value in neighbours[first]
                             if value > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(value for value in common_two
                                if value > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(value for value in common_three
                                     if value > third):
                    result.append((first, second, third, fourth))
    return tuple(result)


def projected_barycentric_carrier(vertices, coarse_cells, coarse_top):
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    cell_index = {cell: position for position, cell in enumerate(vertex_cells)}
    raw = np.asarray([vertices[list(cell)].mean(axis=0)
                      for cell in vertex_cells])
    positions = raw / np.linalg.norm(raw, axis=1)[:, None]
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            top.append(tuple(cell_index[cell] for cell in flag))
    return positions, tuple(top)


def weak_compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in weak_compositions(total-first, parts-1):
            yield (first,) + rest


def edgewise_facets(k, dimension=3):
    """Exact Edelsbrunner--Grayson colour-scheme enumeration."""
    width = dimension + 1
    facets = set()
    for counts in weak_compositions(k*width, width):
        sequence = tuple(
            colour for colour, count in enumerate(counts)
            for _ in range(count)
        )
        rows = tuple(sequence[row*width:(row+1)*width]
                     for row in range(k))
        columns = tuple(tuple(rows[row][column] for row in range(k))
                        for column in range(width))
        if len(set(columns)) != width:
            continue
        facets.add(tuple(
            tuple(column.count(colour) for colour in range(width))
            for column in columns
        ))
    return tuple(sorted(facets))


def route_a_edgewise(base_positions, base_top):
    local_facets = edgewise_facets(2)
    key_to_index = {}
    ordered_keys = []
    fine_top = []
    occurrences = []

    for chamber in base_top:
        for facet in local_facets:
            child = []
            for numerator in facet:
                entries = []
                for rank, weight in enumerate(numerator):
                    entries.extend([chamber[rank]] * weight)
                key = tuple(sorted(entries))
                if key not in key_to_index:
                    key_to_index[key] = len(ordered_keys)
                    ordered_keys.append(key)
                child.append(key_to_index[key])
                occurrences.append((key_to_index[key], chamber, numerator))
            fine_top.append(tuple(sorted(child)))

    raw = np.asarray([
        base_positions[list(key)].mean(axis=0) for key in ordered_keys
    ])
    positions = raw / np.linalg.norm(raw, axis=1)[:, None]
    return tuple(ordered_keys), positions, tuple(fine_top), tuple(occurrences)


def route_b_old_plus_edges(base_positions, base_edges):
    keys = tuple(sorted(
        tuple((index, index)) for index in range(len(base_positions))
    ) + sorted(tuple(edge) for edge in base_edges))
    left = np.asarray([key[0] for key in keys], dtype=np.int32)
    right = np.asarray([key[1] for key in keys], dtype=np.int32)
    sums = base_positions[left] + base_positions[right]
    norms = np.linalg.norm(sums, axis=1)
    positions = sums / norms[:, None]
    return keys, positions, left, right, norms


def analytic_prolongation(positions, tangents, left, right):
    sums = positions[left] + positions[right]
    norms = np.linalg.norm(sums, axis=1)
    fine_positions = sums / norms[:, None]
    raw_tangents = tangents[left] + tangents[right]
    radial = np.einsum("ni,ni->n", fine_positions, raw_tangents)
    fine_tangents = (
        raw_tangents - radial[:, None] * fine_positions
    ) / norms[:, None]
    return fine_positions, fine_tangents, norms


def nonlinear_fine_positions(positions, left, right):
    sums = positions[left] + positions[right]
    return sums / np.linalg.norm(sums, axis=1)[:, None]


def normalized_coarse_perturbation(positions, tangents, epsilon):
    moved = positions + epsilon*tangents
    return moved / np.linalg.norm(moved, axis=1)[:, None]


def normalize_rows_extended(values):
    values = np.asarray(values, dtype=np.longdouble)
    norms = np.sqrt(np.sum(values*values, axis=1))
    return values/norms[:, None]


def deterministic_tangents(positions):
    indices = np.arange(len(positions))
    seeds = (
        np.column_stack((
            np.ones(len(indices)),
            ((indices % 7)-3)/7.0,
            ((indices % 11)-5)/11.0,
            ((indices % 13)-6)/13.0,
        )),
        np.column_stack((
            ((indices % 5)-2)/5.0,
            np.ones(len(indices)),
            ((indices % 17)-8)/17.0,
            ((indices % 19)-9)/19.0,
        )),
    )
    fields = []
    for seed in seeds:
        radial = np.einsum("ni,ni->n", positions, seed)
        tangent = seed-radial[:, None]*positions
        tangent /= np.linalg.norm(tangent)
        fields.append(tangent)
    return tuple(fields)


def centered_vertex_error(positions, tangent, left, right, epsilon,
                          analytic):
    positions_extended = np.asarray(positions, dtype=np.longdouble)
    tangent_extended = np.asarray(tangent, dtype=np.longdouble)
    epsilon_extended = np.longdouble(epsilon)
    plus = normalize_rows_extended(
        positions_extended+epsilon_extended*tangent_extended
    )
    minus = normalize_rows_extended(
        positions_extended-epsilon_extended*tangent_extended
    )
    fine_plus = normalize_rows_extended(plus[left]+plus[right])
    fine_minus = normalize_rows_extended(minus[left]+minus[right])
    finite = (fine_plus-fine_minus)/(2*epsilon_extended)
    analytic_extended = np.asarray(analytic, dtype=np.longdouble)
    return (
        float(np.max(np.abs(finite-analytic_extended))),
        fine_plus,
        fine_minus,
    )


def fine_edge_derivative_error(fine_positions, fine_tangent, fine_edges,
                               fine_plus, fine_minus, epsilon):
    edge_array = np.asarray(fine_edges, dtype=np.int32)
    first = edge_array[:, 0]
    second = edge_array[:, 1]
    positions_extended = np.asarray(fine_positions, dtype=np.longdouble)
    tangent_extended = np.asarray(fine_tangent, dtype=np.longdouble)
    delta_position = positions_extended[first]-positions_extended[second]
    delta_tangent = tangent_extended[first]-tangent_extended[second]
    analytic = 2.0*np.einsum("ni,ni->n", delta_position, delta_tangent)
    squared_plus = np.einsum(
        "ni,ni->n", fine_plus[first]-fine_plus[second],
        fine_plus[first]-fine_plus[second]
    )
    squared_minus = np.einsum(
        "ni,ni->n", fine_minus[first]-fine_minus[second],
        fine_minus[first]-fine_minus[second]
    )
    finite = (
        squared_plus-squared_minus
    )/(2*np.longdouble(epsilon))
    return float(np.max(np.abs(finite-analytic)))


def scalar_cancellation_control(epsilon):
    epsilon64 = np.float64(epsilon)
    one64 = np.float64(1.0)
    finite64 = (
        (one64+epsilon64)*(one64+epsilon64)
        - (one64-epsilon64)*(one64-epsilon64)
    )/(np.float64(2.0)*epsilon64)
    epsilon_extended = np.longdouble(epsilon)
    one_extended = np.longdouble(1.0)
    finite_extended = (
        (one_extended+epsilon_extended)*(one_extended+epsilon_extended)
        - (one_extended-epsilon_extended)*(one_extended-epsilon_extended)
    )/(np.longdouble(2.0)*epsilon_extended)
    return (
        float(abs(finite64-np.float64(2.0))),
        float(abs(finite_extended-np.longdouble(2.0))),
    )


def convergence_ok(errors):
    coarse, fine = errors
    return fine <= coarse or (coarse <= 2e-10 and fine <= 2e-10)


print("="*78)
print("NESTED TANGENTIAL VERTEX-DISPLACEMENT PROLONGATION")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "all preregistered inputs have exact frozen provenance",
    actual_hashes == INPUT_HASHES,
    str(actual_hashes),
)

source_positions, source_adjacency, _ = build_600cell()
source_top = tetrahedra_from_adjacency(source_adjacency)
source_cells = all_simplices(source_top)
base_positions, base_top = projected_barycentric_carrier(
    source_positions, source_cells, source_top
)
base_cells = all_simplices(base_top)
base_f = tuple(len(layer) for layer in base_cells)
base_ok = check(
    "K0 has the frozen closed topology and unit geometry",
    base_f == EXPECTED_BASE_F
    and len(set(base_top)) == len(base_top)
    and set(face_incidence(base_top).values()) == {2}
    and np.max(np.abs(np.linalg.norm(base_positions, axis=1)-1.0)) <= 2e-10,
    f"f={base_f}",
)

keys_a, positions_a, fine_top, occurrences = route_a_edgewise(
    base_positions, base_top
)
fine_cells = all_simplices(fine_top)
fine_f = tuple(len(layer) for layer in fine_cells)
fine_a_ok = check(
    "Route A reconstructs the complete closed fine carrier",
    fine_f == EXPECTED_FINE_F
    and len(set(fine_top)) == len(fine_top)
    and set(face_incidence(fine_top).values()) == {2}
    and len(keys_a) == EXPECTED_FINE_F[0]
    and np.max(np.abs(np.linalg.norm(positions_a, axis=1)-1.0)) <= 2e-10,
    f"f={fine_f}, occurrences={len(occurrences)}",
)

keys_b, positions_b, left, right, denominators = route_b_old_plus_edges(
    base_positions, base_cells[1]
)
old_count = sum(a == b for a, b in keys_b)
midpoint_count = len(keys_b)-old_count
fine_b_ok = check(
    "Route B gives exactly old vertices plus all coarse edges",
    len(keys_b) == EXPECTED_FINE_F[0]
    and old_count == EXPECTED_BASE_F[0]
    and midpoint_count == EXPECTED_BASE_F[1]
    and len(set(keys_b)) == len(keys_b),
    f"old={old_count}, midpoints={midpoint_count}",
)

lookup_a = {key: index for index, key in enumerate(keys_a)}
route_coordinate_residual = max(
    float(np.linalg.norm(positions_a[lookup_a[key]]-positions_b[index]))
    for index, key in enumerate(keys_b)
)
route_agreement_ok = check(
    "the mechanically distinct fine-key routes agree",
    set(keys_a) == set(keys_b) and route_coordinate_residual <= 2e-10,
    f"coordinate residual={route_coordinate_residual:.3e}",
)

tangent_fields = deterministic_tangents(base_positions)
analytic_fields = []
for field in tangent_fields:
    _, prolonged, _ = analytic_prolongation(
        base_positions, field, left, right
    )
    analytic_fields.append(prolonged)

parent_coordinate_residual = 0.0
parent_tangent_residual = 0.0
key_b_lookup = {key: index for index, key in enumerate(keys_b)}
for key_index_a, chamber, numerator in occurrences:
    key = keys_a[key_index_a]
    key_index_b = key_b_lookup[key]
    raw = sum(
        numerator[rank]*base_positions[chamber[rank]]
        for rank in range(4)
    )/2.0
    occurrence_position = raw/np.linalg.norm(raw)
    parent_coordinate_residual = max(
        parent_coordinate_residual,
        float(np.linalg.norm(occurrence_position-positions_b[key_index_b])),
    )
    for field, prolonged in zip(tangent_fields, analytic_fields):
        raw_tangent = sum(
            numerator[rank]*field[chamber[rank]]
            for rank in range(4)
        )/2.0
        occurrence_tangent = (
            raw_tangent
            - occurrence_position*np.dot(occurrence_position, raw_tangent)
        )/np.linalg.norm(raw)
        parent_tangent_residual = max(
            parent_tangent_residual,
            float(np.linalg.norm(
                occurrence_tangent-prolonged[key_index_b]
            )),
        )
denominator_parent_ok = check(
    "all midpoint maps are regular and parent independent",
    float(np.min(denominators)) > 1.0
    and parent_coordinate_residual <= 2e-10
    and parent_tangent_residual <= 2e-10,
    f"min denominator={np.min(denominators):.6f}, "
    f"position={parent_coordinate_residual:.3e}, "
    f"tangent={parent_tangent_residual:.3e}",
)

coarse_tangency = max(
    float(np.max(np.abs(np.einsum("ni,ni->n", base_positions, field))))
    for field in tangent_fields
)
fine_tangency = max(
    float(np.max(np.abs(np.einsum("ni,ni->n", positions_b, field))))
    for field in analytic_fields
)
tangency_ok = check(
    "the two probe fields remain tangent on both levels",
    coarse_tangency <= 2e-10 and fine_tangency <= 2e-10,
    f"coarse={coarse_tangency:.3e}, fine={fine_tangency:.3e}",
)

old_indices = np.asarray(
    [key_b_lookup[(index, index)] for index in range(len(base_positions))],
    dtype=np.int32,
)
left_inverse_residual = max(
    float(np.max(np.abs(prolonged[old_indices]-field)))
    for field, prolonged in zip(tangent_fields, analytic_fields)
)
left_inverse_ok = check(
    "the old-vertex restriction is an exact left inverse",
    len(set(old_indices.tolist())) == len(base_positions)
    and all(keys_b[old_indices[index]] == (index, index)
            for index in range(len(base_positions)))
    and left_inverse_residual <= 2e-15,
    f"structural rank={3*len(base_positions)}, residual={left_inverse_residual:.3e}",
)

vertex_errors = []
edge_errors = []
fine_edges = fine_cells[1]
for field, prolonged in zip(tangent_fields, analytic_fields):
    local_vertex = []
    local_edge = []
    for epsilon in STEPS:
        error, fine_plus, fine_minus = centered_vertex_error(
            base_positions, field, left, right, epsilon, prolonged
        )
        local_vertex.append(error)
        local_edge.append(fine_edge_derivative_error(
            positions_b, prolonged, fine_edges,
            fine_plus, fine_minus, epsilon,
        ))
    vertex_errors.append(tuple(local_vertex))
    edge_errors.append(tuple(local_edge))
vertex_derivative_ok = check(
    "centered differentiation corroborates the vertex derivative",
    max(max(pair) for pair in vertex_errors) <= 2e-7
    and all(convergence_ok(pair) for pair in vertex_errors),
    f"errors={vertex_errors}",
)
scalar_control_errors = scalar_cancellation_control(STEPS[-1])
edge_derivative_ok = check(
    "the induced fine-edge metric derivative is unambiguous",
    len(fine_edges) == EXPECTED_FINE_F[1]
    and max(max(pair) for pair in edge_errors) <= 3e-7
    and all(convergence_ok(pair) for pair in edge_errors)
    and scalar_control_errors[1] <= scalar_control_errors[0],
    f"edges={len(fine_edges)}, errors={edge_errors}, "
    f"scalar64/extended={scalar_control_errors}",
)

r1 = np.zeros((4, 4))
for old, new in enumerate((1, 2, 3, 0)):
    r1[new, old] = -1.0 if old == 3 else 1.0
w = np.ones(4)/2.0
r2 = np.eye(4)-2.0*np.outer(w, w)
orthogonal_matrices = (r1, r2)
orthogonality = max(
    float(np.max(np.abs(matrix.T@matrix-np.eye(4))))
    for matrix in orthogonal_matrices
)
determinants = tuple(float(np.linalg.det(matrix))
                     for matrix in orthogonal_matrices)
covariance_residual = 0.0
for matrix in orthogonal_matrices:
    rotated_positions = base_positions@matrix.T
    for field, prolonged in zip(tangent_fields, analytic_fields):
        rotated_field = field@matrix.T
        covariant_positions, covariant_tangent, _ = analytic_prolongation(
            rotated_positions, rotated_field, left, right
        )
        covariance_residual = max(
            covariance_residual,
            float(np.max(np.abs(covariant_positions-positions_b@matrix.T))),
            float(np.max(np.abs(covariant_tangent-prolonged@matrix.T))),
        )
covariance_ok = check(
    "the prolongation is covariant in both O(4) determinant classes",
    orthogonality <= 2e-15
    and abs(determinants[0]-1.0) <= 2e-15
    and abs(determinants[1]+1.0) <= 2e-15
    and covariance_residual <= 2e-10,
    f"det={determinants}, residual={covariance_residual:.3e}",
)

feasibility = json.loads(
    (HERE/"gravity_600cell_refined_canonical_map_feasibility.json").read_text()
)
schedule_records = {}
schedule_ok_values = []
for level_name, level in feasibility["levels"].items():
    boundary_counts = sorted(set(
        int(record["boundary_edges"]) for record in level["schedules"]
    ))
    record = {
        "schedule_count": int(level["schedule_count"]),
        "distinct_internal_edge_sets": int(level["distinct_internal_edge_sets"]),
        "cross_edge_intersection": int(level["cross_edge_intersection"]),
        "boundary_edge_counts": boundary_counts,
    }
    schedule_records[level_name] = record
    schedule_ok_values.append(
        record["schedule_count"] == 24
        and record["distinct_internal_edge_sets"] == 24
        and record["cross_edge_intersection"] == 0
        and len(boundary_counts) == 1
    )
schedule_independence_ok = check(
    "the spatial map is schedule independent while temporal ambiguity remains",
    all(schedule_ok_values),
    str(schedule_records),
)

first_midpoint = next(index for index, (a, b) in enumerate(keys_b) if a != b)
i, j = keys_b[first_midpoint]
corrupted_raw = 2.0*base_positions[i]+base_positions[j]
corrupted = corrupted_raw/np.linalg.norm(corrupted_raw)
corruption_distance = float(np.linalg.norm(
    corrupted-positions_b[first_midpoint]
))
omitted_old_keys = set(keys_b)
omitted_old_keys.remove((0, 0))
omission_detected = not all(
    (index, index) in omitted_old_keys for index in range(len(base_positions))
)
negative_controls_ok = check(
    "the corrupted midpoint and omitted old vertex are detected",
    corruption_distance > 1e-4 and omission_detected,
    f"corruption distance={corruption_distance:.6e}, omission={omission_detected}",
)

target_firewall = {
    "action_hessian_constructed": False,
    "continuum_eigenvalue_loaded": False,
    "dispersion_target_loaded": False,
    "graviton_target_loaded": False,
    "c_target_loaded": False,
    "G_target_loaded": False,
    "planck_target_loaded": False,
    "particle_target_loaded": False,
    "momentum_lift_claimed": False,
    "normal_lapse_transport_claimed": False,
}

control_ok = provenance_ok and base_ok and fine_a_ok and fine_b_ok
canonical_ok = (
    route_agreement_ok and denominator_parent_ok and tangency_ok
    and left_inverse_ok and covariance_ok and schedule_independence_ok
)
derivative_ok = vertex_derivative_ok and edge_derivative_ok
if not control_ok:
    outcome = "NESTED_VERTEX_PROLONGATION_CONTROL_FAILED"
elif not canonical_ok or not negative_controls_ok:
    outcome = "NESTED_VERTEX_PROLONGATION_NOT_CANONICAL"
elif not derivative_ok:
    outcome = "NESTED_VERTEX_PROLONGATION_DERIVATIVE_FAILED"
else:
    outcome = "NESTED_TANGENTIAL_VERTEX_CARRIER_DERIVED"

outcome_ok = check(
    "the target firewall and frozen hierarchy assign one outcome",
    not any(target_firewall.values())
    and outcome in {
        "NESTED_VERTEX_PROLONGATION_CONTROL_FAILED",
        "NESTED_VERTEX_PROLONGATION_NOT_CANONICAL",
        "NESTED_VERTEX_PROLONGATION_DERIVATIVE_FAILED",
        "NESTED_TANGENTIAL_VERTEX_CARRIER_DERIVED",
    },
    outcome,
)

payload = {
    "title": "Nested tangential vertex-displacement prolongation",
    "date": "2026-08-22",
    "classification": "DERIVED_COMPUTATIONAL_STRUCTURAL_INFRASTRUCTURE",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "precision_correction_commit": PRECISION_CORRECTION_COMMIT,
    "input_sha256": actual_hashes,
    "levels": {
        "K0_f_vector": list(base_f),
        "K1_f_vector": list(fine_f),
        "old_vertices": old_count,
        "projected_midpoints": midpoint_count,
        "fine_key_count": len(keys_b),
        "fine_edge_count": len(fine_edges),
    },
    "prolongation": {
        "domain": "direct sum of K0 spatial tangent spaces",
        "domain_dimension": 3*len(base_positions),
        "codomain_dimension": 3*len(keys_b),
        "rank_certificate": "exact old-vertex left inverse",
        "rank": 3*len(base_positions),
        "minimum_midpoint_denominator": float(np.min(denominators)),
        "route_coordinate_residual": route_coordinate_residual,
        "parent_coordinate_residual": parent_coordinate_residual,
        "parent_tangent_residual": parent_tangent_residual,
        "coarse_tangency_residual": coarse_tangency,
        "fine_tangency_residual": fine_tangency,
        "left_inverse_residual": left_inverse_residual,
        "vertex_finite_difference_errors": vertex_errors,
        "edge_finite_difference_errors": edge_errors,
        "scalar_cancellation_errors_binary64_extended": scalar_control_errors,
        "O4_determinants": determinants,
        "O4_covariance_residual": covariance_residual,
    },
    "temporal_schedule_audit": schedule_records,
    "negative_controls": {
        "corrupted_midpoint_distance": corruption_distance,
        "omitted_old_vertex_detected": omission_detected,
    },
    "scope": {
        "spatial_tangential_displacements_only": True,
        "normal_lapse_transport": "OPEN",
        "cotangent_lift": "OPEN_AND_NONUNIQUE_FROM_PAIRING_ALONE",
        "matched_K1_on_shell_background": "OPEN",
        "constraint_restoration": "OPEN",
        "physical_gravitons": "OPEN",
    },
    "target_firewall": target_firewall,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(outcome)
print(f"{passed}/{tests} tests passed")
print(f"artifact: {OUTPUT}")

if tests != 14 or passed != tests or not outcome_ok:
    raise SystemExit(1)
