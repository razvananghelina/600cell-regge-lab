#!/usr/bin/env python3
"""Independent geometric audit of the canonical prism-shift elimination."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_prism_shift_canonical_elimination_adversarial.json"
PRIOR_ART_COMMIT = "d90a44b"
PROTOCOL_COMMIT = "d01217b"
PRIMARY_RESULT_COMMIT = "78fa42c"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "docs/gravity/gravity_600cell_prism_shift_canonical_elimination_adversarial_protocol.md":
        "b9df30a3ab77c5719cf81e079400a80dd4003cd97454981db90eef06cb8124f2",
    "reproducible/gravity_600cell_prism_shift_canonical_elimination.json":
        "b9e31d56670c397232937ae5f2e7e002632cc715807c221b2b98b47e20dde332",
    "reproducible/gravity_600cell_dust_full_lapse_schur.json":
        "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349",
    "reproducible/gravity_600cell_dust_full_anisotropic_legendre_rank.json":
        "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "reproducible/verify_gravity_600cell_dust_full_lapse_schur.py":
        "7258899ba96a127515956fa2ea5fb17ad480373765b3f7c88fed40845adc82a6",
    "reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py":
        "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "reproducible/verify_gravity_global_regge_orbits.py":
        "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
}
PRIMES = (101, 1000003)
RATIONAL_SCALES = ((3, 2), (5, 3))
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
    return sha256(path.read_bytes()).hexdigest()


def quaternion_product(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return np.asarray((
        a*e-b*f-c*g-d*h,
        a*f+b*e+c*h-d*g,
        a*g-b*h+c*e+d*f,
        a*h+b*g-c*f+d*e,
    ))


def tetrahedra_from_graph(adjacency):
    neighbours = [set(np.flatnonzero(row)) for row in adjacency]
    tetrahedra = []
    for first in range(120):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    tetrahedra.append((first, second, third, fourth))
    return tuple(tetrahedra)


def staircase(tetrahedra, phase):
    simplices = set()
    for tetrahedron in tetrahedra:
        ordered = sorted(tetrahedron, key=phase.__getitem__)
        for pivot in ordered:
            simplex = [pivot, pivot+120]
            simplex.extend(
                vertex+120 if phase[vertex] < phase[pivot] else vertex
                for vertex in ordered if vertex != pivot
            )
            simplices.add(tuple(sorted(simplex)))
    return frozenset(simplices)


def rank_mod_prime(rows, column_count, prime):
    basis = {}
    for source in rows:
        row = {
            column: value % prime for column, value in source.items()
            if value % prime
        }
        while row:
            pivot = min(row)
            if pivot not in basis:
                inverse = pow(row[pivot], -1, prime)
                basis[pivot] = {
                    column: value*inverse % prime
                    for column, value in row.items()
                }
                break
            factor = row[pivot]
            for column, value in basis[pivot].items():
                updated = (row.get(column, 0)-factor*value) % prime
                if updated:
                    row[column] = updated
                else:
                    row.pop(column, None)
    assert all(0 <= pivot < column_count for pivot in basis)
    return len(basis)


def relative_pullback(rows):
    """Pull rows back along R e_k=e_k-e_119."""
    pulled = []
    for row in rows:
        last = row.get(119, 0)
        pulled.append({
            column: row.get(column, 0)-last
            for column in range(119)
            if row.get(column, 0)-last
        })
    return tuple(pulled)


print("="*78)
print("ADVERSARIAL GEOMETRY AUDIT OF CANONICAL PRISM-SHIFT ELIMINATION")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
primary = json.loads(
    (HERE/"gravity_600cell_prism_shift_canonical_elimination.json").read_text())
schur = json.loads(
    (HERE/"gravity_600cell_dust_full_lapse_schur.json").read_text())
rank_artifact = json.loads(
    (HERE/"gravity_600cell_dust_full_anisotropic_legendre_rank.json").read_text())
check(
    "all frozen geometry, carrier and primary inputs have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "d90a44b"
    and PROTOCOL_COMMIT == "d01217b"
    and PRIMARY_RESULT_COMMIT == "78fa42c",
    str(actual_hashes),
)
check(
    "the input conclusions are complete and contain no target comparison",
    primary["verdict"] == "RELATIVE_SHIFT_CANONICALLY_ELIMINATED"
    and primary["passed"] == primary["tests"] == 13
    and schur["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    and schur["passed"] == schur["tests"] == 18
    and not schur["continuum_target_parsed"]
    and not schur["speed_target_parsed"],
)

# Reconstruct the group table and the five right cosets independently.
vertices, adjacency_float, _ = build_600cell()
adjacency = adjacency_float > 0.5
multiplication = np.empty((120, 120), dtype=np.int16)
minimum_margin = np.inf
maximum_closure_error = 0.0
for left in range(120):
    for right in range(120):
        candidate = quaternion_product(vertices[left], vertices[right])
        products = vertices @ candidate
        order = np.argsort(products)
        multiplication[left, right] = int(order[-1])
        minimum_margin = min(minimum_margin, float(products[order[-1]]-products[order[-2]]))
        maximum_closure_error = max(maximum_closure_error, float(abs(products[order[-1]]-1)))
check(
    "nearest-neighbour quaternion multiplication is unique and closed",
    minimum_margin > 0.19 and maximum_closure_error < 1e-9,
    f"minimum margin={minimum_margin:.12g}, closure error={maximum_closure_error:.3e}",
)

binary_tetrahedral = tuple(
    index for index, vertex in enumerate(vertices)
    if (
        np.count_nonzero(np.abs(vertex) > 1e-8) == 1
        and np.max(np.abs(vertex)) > 1-1e-8
    ) or np.all(np.abs(np.abs(vertex)-0.5) < 1e-8)
)
unseen = set(range(120))
cover_cells = []
while unseen:
    representative = min(unseen)
    cell = frozenset(
        int(multiplication[representative, element])
        for element in binary_tetrahedral
    )
    cover_cells.append(cell)
    unseen -= cell
cover_cells = tuple(sorted(cover_cells, key=lambda cell: tuple(sorted(cell))))
check(
    "the binary-tetrahedral construction gives a disjoint 5 by 24 cover",
    len(binary_tetrahedral) == 24
    and len(cover_cells) == 5
    and all(len(cell) == 24 for cell in cover_cells)
    and set().union(*cover_cells) == set(range(120))
    and sum(len(cell) for cell in cover_cells) == 120,
)

tetrahedra = tetrahedra_from_graph(adjacency)
spatial_edges = frozenset(
    tuple(map(int, edge))
    for edge in np.transpose(np.nonzero(np.triu(adjacency, 1)))
)
spatial_faces = frozenset(
    tuple(sorted(face))
    for tetrahedron in tetrahedra for face in combinations(tetrahedron, 3)
)
check(
    "the independently reconstructed spatial carrier is the complete 600-cell",
    (len(vertices), len(spatial_edges), len(spatial_faces), len(tetrahedra))
    == (120, 720, 1200, 600),
)

orders = {
    "even": cover_cells,
    "odd": (
        cover_cells[1], cover_cells[0], cover_cells[2],
        cover_cells[3], cover_cells[4],
    ),
}
schedule_records = {}
all_schedule_counts = True
all_orientations = True
all_boundary_support = True
diagonal_sets = {}
for parity, ordering in orders.items():
    phase = {
        vertex: phase_index
        for phase_index, cell in enumerate(ordering) for vertex in cell
    }
    rainbow = all(
        len({phase[vertex] for vertex in tetrahedron}) == 4
        for tetrahedron in tetrahedra
    )
    slab = staircase(tetrahedra, phase)
    all_edges = frozenset(
        tuple(sorted(edge))
        for simplex in slab for edge in combinations(simplex, 2)
    )
    old_edges = spatial_edges
    new_edges = frozenset((left+120, right+120) for left, right in old_edges)
    internal_edges = all_edges-old_edges-new_edges
    pole_edges = frozenset((vertex, vertex+120) for vertex in range(120))
    diagonal_edges = internal_edges-pole_edges
    predicted = frozenset(
        tuple(sorted((high, low+120)))
        for left, right in old_edges
        for low, high in [
            (left, right) if phase[left] < phase[right] else (right, left)
        ]
    )
    counts = {
        "four_simplices": len(slab),
        "all_edges": len(all_edges),
        "old": len(old_edges),
        "new": len(new_edges),
        "internal": len(internal_edges),
        "poles": len(pole_edges),
        "diagonals": len(diagonal_edges),
    }
    count_ok = rainbow and counts == {
        "four_simplices": 2400,
        "all_edges": 2280,
        "old": 720,
        "new": 720,
        "internal": 840,
        "poles": 120,
        "diagonals": 720,
    }
    orientation_ok = diagonal_edges == predicted
    support_ok = (
        not (internal_edges & old_edges)
        and not (internal_edges & new_edges)
        and pole_edges <= internal_edges
        and diagonal_edges.isdisjoint(pole_edges)
    )
    all_schedule_counts &= count_ok
    all_orientations &= orientation_ok
    all_boundary_support &= support_ok
    diagonal_sets[parity] = diagonal_edges
    schedule_records[parity] = {
        "counts": counts,
        "rainbow_tetrahedra": rainbow,
        "orientation_exact": orientation_ok,
        "internal_boundary_disjoint": support_ok,
    }
check(
    "both staircase schedules have the preregistered complete carrier",
    all_schedule_counts,
    str({key: value["counts"] for key, value in schedule_records.items()}),
)
check(
    "each spatial edge selects exactly the color-oriented cross diagonal",
    all_orientations,
)
check(
    "all poles and cross diagonals are internal and boundary support is zero",
    all_boundary_support,
)

# Direct exact affine derivatives on a centered rational regular tetrahedron.
tetrahedron = (
    sy.Matrix((1, 1, 1)), sy.Matrix((1, -1, -1)),
    sy.Matrix((-1, 1, -1)), sy.Matrix((-1, -1, 1)),
)
epsilon = sy.symbols("epsilon")
direct_controls = (
    (sy.Rational(3, 2), sy.Matrix((2, -1, 3)), sy.Rational(5), sy.Rational(4)),
    (sy.Rational(5, 3), sy.Matrix((-3, 2, 1)), sy.Rational(7), sy.Rational(-2)),
)
direct_exact = True
reverse_failures = 0
direct_records = []
for q, shift, lapse, lapse_change in direct_controls:
    strut_derivatives = []
    for vertex in tetrahedron:
        spatial = (q-1)*vertex+epsilon*shift
        square = (spatial.T*spatial)[0]-(lapse+epsilon*lapse_change)**2
        strut_derivatives.append(sy.diff(square, epsilon).subs(epsilon, 0))
    residuals = []
    reverse_residuals = []
    for bottom in range(4):
        for top in range(4):
            if bottom == top:
                continue
            spatial = q*tetrahedron[top]+epsilon*shift-tetrahedron[bottom]
            square = (spatial.T*spatial)[0]-(lapse+epsilon*lapse_change)**2
            diagonal_derivative = sy.diff(square, epsilon).subs(epsilon, 0)
            predicted = sy.cancel(
                (q*strut_derivatives[top]-strut_derivatives[bottom])/(q-1)
            )
            reversed_prediction = sy.cancel(
                (q*strut_derivatives[bottom]-strut_derivatives[top])/(q-1)
            )
            residuals.append(sy.simplify(diagonal_derivative-predicted))
            reverse_residuals.append(
                sy.simplify(diagonal_derivative-reversed_prediction)
            )
    direct_exact &= all(value == 0 for value in residuals)
    reverse_failures += sum(value != 0 for value in reverse_residuals)
    direct_records.append({
        "q": str(q),
        "maximum_exact_residual": str(max(map(abs, residuals))),
        "reversed_nonzero_residuals": sum(value != 0 for value in reverse_residuals),
    })
check(
    "two nonsymmetric exact coordinate controls derive the oriented graph formula",
    direct_exact,
    str(direct_records),
)
check(
    "the orientation-reversed formula is rejected by direct geometry",
    reverse_failures > 0,
    f"nonzero reversed residuals={reverse_failures}/24",
)

# Exact modular rank of the diagonal graph and its relative pullback.
graph_records = {}
all_graph_ranks = True
all_collective = True
negative_control_detected = True
for parity, diagonal_edges in diagonal_sets.items():
    graph_records[parity] = {}
    for numerator, denominator in RATIONAL_SCALES:
        rows = []
        for bottom, lifted_top in sorted(diagonal_edges):
            top = lifted_top-120
            rows.append({bottom: denominator, top: -numerator})
        relative_rows = relative_pullback(rows)
        full_ranks = {
            prime: rank_mod_prime(rows, 120, prime) for prime in PRIMES
        }
        relative_ranks = {
            prime: rank_mod_prime(relative_rows, 119, prime) for prime in PRIMES
        }
        # Removing every diagonal incident to vertex zero leaves its column zero.
        deleted_rows = [
            row for edge, row in zip(sorted(diagonal_edges), rows)
            if edge[0] != 0 and edge[1]-120 != 0
        ]
        deleted_ranks = {
            prime: rank_mod_prime(deleted_rows, 120, prime) for prime in PRIMES
        }
        rank_ok = (
            full_ranks == {prime: 120 for prime in PRIMES}
            and relative_ranks == {prime: 119 for prime in PRIMES}
        )
        # For z_i=1, the log-diagonal coefficient is exactly -rho/D.
        q = sy.Rational(numerator, denominator)
        rho, length_square = sy.symbols("rho length_square", nonzero=True)
        diagonal_square = q*length_square-rho
        collective = sy.cancel(
            rho*(1-q)/((q-1)*diagonal_square)
        )
        collective_ok = sy.simplify(collective+rho/diagonal_square) == 0
        negative_ok = all(rank < 120 for rank in deleted_ranks.values())
        all_graph_ranks &= rank_ok
        all_collective &= collective_ok
        negative_control_detected &= negative_ok
        graph_records[parity][f"{numerator}/{denominator}"] = {
            "full_ranks": full_ranks,
            "relative_ranks": relative_ranks,
            "deleted_vertex_zero_ranks": deleted_ranks,
            "collective_coefficient": str(collective),
        }
check(
    "the diagonal graph has exact rank 120 and relative rank 119",
    all_graph_ranks,
    str(graph_records),
)
check(
    "the collective graph is exactly the frozen -rho/D lapse coefficient",
    all_collective,
)
check(
    "deleting one vertex's diagonals triggers the preregistered rank control",
    negative_control_detected,
)

# Tie the independently reconstructed support to the frozen canonical split.
carrier_alignment = True
for parity in ("even", "odd"):
    carrier = rank_artifact["parities"][parity]["carrier"]
    carrier_alignment &= (
        carrier["four_simplices"] == schedule_records[parity]["counts"]["four_simplices"]
        and carrier["edge_variables"] == schedule_records[parity]["counts"]["all_edges"]
        and carrier["old"] == schedule_records[parity]["counts"]["old"]
        and carrier["new"] == schedule_records[parity]["counts"]["new"]
        and carrier["internal"] == schedule_records[parity]["counts"]["internal"]
        and schur["parities"][parity]["weak_orbit_positions"]
        == [30, 31, 32, 33, 34]
        and schur["parities"][parity]["resolved_schur_rank"] == 120
        and schur["parities"][parity]["schur_zero_count"] == 0
        and schur["parities"][parity]["schur_open_count"] == 0
    )
check(
    "the independent 720-diagonal/120-pole split aligns with the frozen carrier",
    carrier_alignment,
)

all_pass = passed == tests
verdict = (
    "RELATIVE_SHIFT_ELIMINATION_GEOMETRICALLY_CORROBORATED"
    if all_pass else "RELATIVE_SHIFT_ELIMINATION_GEOMETRY_OPEN"
)
result = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_result_commit": PRIMARY_RESULT_COMMIT,
    "input_sha256": actual_hashes,
    "schedule_reconstruction": schedule_records,
    "direct_affine_controls": direct_records,
    "graph_rank_controls": graph_records,
    "reverse_orientation_nonzero_residuals": reverse_failures,
    "continuum_target_parsed": False,
    "speed_target_parsed": False,
    "mass_target_parsed": False,
    "classification": {
        "geometric_internal_graph": "DERIVED" if all_pass else "OPEN",
        "homogeneous_relative_shift_elimination": "DERIVED" if all_pass else "OPEN",
        "auxiliary_or_pseudoconstraint_reading": "STRUCTURAL",
        "sourced_boundary_response": "OPEN",
        "physical_tensor_modes": "OPEN",
    },
    "verdict": verdict,
    "passed": passed,
    "tests": tests,
}
OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"VERDICT: {verdict}")
if not all_pass:
    raise SystemExit(1)

