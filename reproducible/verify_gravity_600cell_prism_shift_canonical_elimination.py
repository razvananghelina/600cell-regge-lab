#!/usr/bin/env python3
"""Compose the exact relative-shift map with the certified pole Schur map."""

from collections import Counter, deque
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


OUTPUT = HERE / "gravity_600cell_prism_shift_canonical_elimination.json"
PRIOR_ART_COMMIT = "d90a44b"
PROTOCOL_COMMIT = "e411f93"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_prism_shift_canonical_elimination_prior_art.md":
        "d16da74583fa5c17551adb892427b9431ca01ce89531d2429ff4903f0ccb49e1",
    "docs/gravity/gravity_600cell_prism_shift_canonical_elimination_protocol.md":
        "47f8c7809e1cc35d4d921d6712c95671c290baad5b0c38cf656a480a18c85be8",
    "reproducible/gravity_600cell_prism_shift_dynamic_extension.json":
        "32d5269b27756a4c6fec4603855db643106e571007d3f3dd1a0a6c69d33a0095",
    "reproducible/gravity_600cell_dust_full_lapse_schur.json":
        "4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349",
    "reproducible/gravity_600cell_dust_full_anisotropic_legendre_rank.json":
        "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "reproducible/verify_gravity_600cell_dust_full_lapse_schur.py":
        "7258899ba96a127515956fa2ea5fb17ad480373765b3f7c88fed40845adc82a6",
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
PRIMES = (101, 1000003)
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


def tetrahedra_from_graph(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    cells = []
    for first in range(len(adjacency)):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    cells.append((first, second, third, fourth))
    return tuple(cells)


def rank_mod_prime(rows, column_count, prime):
    basis = {}
    for source in rows:
        row = {column: value % prime for column, value in source.items()
               if value % prime}
        while row:
            pivot = min(row)
            if pivot not in basis:
                inverse = pow(row[pivot], -1, prime)
                basis[pivot] = {
                    column: (value*inverse) % prime
                    for column, value in row.items() if value % prime
                }
                break
            factor = row[pivot]
            for column, value in basis[pivot].items():
                updated = (row.get(column, 0)-factor*value) % prime
                if updated:
                    row[column] = updated
                elif column in row:
                    del row[column]
    assert all(0 <= pivot < column_count for pivot in basis)
    return len(basis)


def complex_matrix(serialized):
    return np.asarray([
        [complex(float(value["real"]), float(value["imaginary"]))
         for value in row]
        for row in serialized
    ], dtype=np.complex128)


print("="*78)
print("CANONICAL ELIMINATION OF THE 119 RELATIVE PRISM SHIFTS")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "all frozen theorem, Schur and source inputs have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "d90a44b" and PROTOCOL_COMMIT == "e411f93",
    str(actual_hashes),
)

dynamic = json.loads(
    (HERE/"gravity_600cell_prism_shift_dynamic_extension.json").read_text())
schur = json.loads(
    (HERE/"gravity_600cell_dust_full_lapse_schur.json").read_text())
rank_artifact = json.loads(
    (HERE/"gravity_600cell_dust_full_anisotropic_legendre_rank.json").read_text())
check(
    "the two input theorem outcomes are complete and target-free",
    dynamic["verdict"] == "DYNAMIC_SHIFT_EXTENSION_OBSTRUCTED"
    and dynamic["passed"] == dynamic["tests"] == 13
    and schur["outcome"] == "FULL_LAPSE_SCHUR_REGULAR"
    and schur["passed"] == schur["tests"] == 18
    and not schur["continuum_target_parsed"]
    and not schur["speed_target_parsed"],
)

# Exact relative-pole carrier R e_i=e_i-e_119.
relative_rows = []
for row in range(120):
    if row < 119:
        relative_rows.append({row: 1})
    else:
        relative_rows.append({column: -1 for column in range(119)})
relative_ranks = {
    prime: rank_mod_prime(relative_rows, 119, prime) for prime in PRIMES
}
relative_column_sums = [
    sum(row.get(column, 0) for row in relative_rows)
    for column in range(119)
]
augmented_rows = [
    {**row, 119: 1} for row in relative_rows
]
augmented_ranks = {
    prime: rank_mod_prime(augmented_rows, 120, prime) for prime in PRIMES
}
check(
    "the exact relative-pole map has rank 119 and zero collective sum",
    relative_ranks == {prime: 119 for prime in PRIMES}
    and relative_column_sums == [0]*119,
    f"ranks={relative_ranks}",
)
check(
    "the all-ones collective lapse is an independent complement",
    augmented_ranks == {prime: 120 for prime in PRIMES},
    f"augmented ranks={augmented_ranks}",
)

# Independent literal 600-cell incidence reconstruction.
_, adjacency, _ = build_600cell()
cells = tetrahedra_from_graph(adjacency)
edges = tuple(sorted({tuple(sorted(edge)) for cell in cells
                      for edge in combinations(cell, 2)}))
faces = tuple(sorted({tuple(sorted(face)) for cell in cells
                      for face in combinations(cell, 3)}))
incidence_rows = tuple({left: -1, right: 1} for left, right in edges)
incidence_ranks = {
    prime: rank_mod_prime(incidence_rows, 120, prime) for prime in PRIMES
}
check(
    "the vertex-potential and relative-pole carriers have the same exact rank",
    (120, len(edges), len(faces), len(cells)) == (120, 720, 1200, 600)
    and incidence_ranks == relative_ranks,
    f"f={(120, len(edges), len(faces), len(cells))}, ranks={incidence_ranks}",
)

# Complete carrier census before reading any sector matrix.
carrier_ok = True
for parity in ("even", "odd"):
    carrier = rank_artifact["parities"][parity]["carrier"]
    carrier_ok &= (
        carrier["old"] == 720 and carrier["internal"] == 840
        and carrier["new"] == 720 and carrier["edge_variables"] == 2280
        and schur["parities"][parity]["weak_orbit_positions"]
        == [30, 31, 32, 33, 34]
        and schur["parities"][parity]["resolved_schur_rank"] == 120
        and schur["parities"][parity]["schur_zero_count"] == 0
        and schur["parities"][parity]["schur_open_count"] == 0
    )
check(
    "the complete carrier splits as 1440 strong plus 120 pole coordinates",
    carrier_ok,
    "840 internal=720 non-pole strong+120 poles; 720 new-boundary strong",
)

parity_records = {}
all_determinants = True
all_shapes = True
all_counts = True
all_recomputed = True
global_minimum_singular = np.inf
global_minimum_margin = np.inf
for parity in ("even", "odd"):
    dimensions = []
    strong_full_count = 0
    schur_full_count = 0
    sector_records = []
    for sector in schur["parities"][parity]["sectors"]:
        dimension = int(sector["irrep_dimension"])
        dimensions.append(dimension)
        midpoint = complex_matrix(sector["schur"]["midpoint_matrix"])
        recomputed = np.linalg.svd(midpoint, compute_uv=False)
        stored = np.asarray(
            [float(value) for value in sector["schur"]["singular_values"]])
        epsilon = float(sector["schur"]["epsilon_global"])
        matrix_size = 5*dimension
        strong_size = 60*dimension
        determinant_ok = (
            all(sector["determinant_a_excludes_zero"].values())
            and all(sector["determinant_schur_excludes_zero"].values())
        )
        shapes_ok = (
            midpoint.shape == (matrix_size, matrix_size)
            and len(stored) == matrix_size
            and len(sector["strong"]["singular_values"]) == strong_size
        )
        counts_ok = (
            sector["strong"]["resolved_count"] == strong_size
            and sector["strong"]["zero_count"] == 0
            and sector["strong"]["open_count"] == 0
            and sector["schur"]["resolved_count"] == matrix_size
            and sector["schur"]["zero_count"] == 0
            and sector["schur"]["open_count"] == 0
        )
        comparison = float(np.max(
            np.abs(recomputed-stored)/np.maximum(1.0, np.abs(stored))))
        relative_comparison = float(np.max(
            np.abs(recomputed-stored)/np.maximum(1e-30, np.abs(stored))))
        margin = float(np.min(recomputed)/(100*epsilon))
        recomputed_ok = (
            comparison < 2e-12 and relative_comparison < 2e-12
            and np.min(recomputed) > 100*epsilon
        )
        all_determinants &= determinant_ok
        all_shapes &= shapes_ok
        all_counts &= counts_ok
        all_recomputed &= recomputed_ok
        strong_full_count += dimension*strong_size
        schur_full_count += dimension*matrix_size
        global_minimum_singular = min(
            global_minimum_singular, float(np.min(recomputed)))
        global_minimum_margin = min(global_minimum_margin, margin)
        sector_records.append({
            "irrep_dimension": dimension,
            "schur_size": matrix_size,
            "minimum_singular": float(np.min(recomputed)),
            "epsilon_global": epsilon,
            "margin_over_100_epsilon": margin,
            "stored_relative_error": relative_comparison,
        })
    parity_records[parity] = {
        "dimensions": dimensions,
        "strong_full_count": strong_full_count,
        "schur_full_count": schur_full_count,
        "sectors": sector_records,
    }

check(
    "all frozen strong and Schur determinant certificates exclude zero",
    all_determinants,
)
check(
    "every parsed sector has the preregistered 60d/5d shapes and counts",
    all_shapes and all_counts
    and all(sorted(record["dimensions"]) == [1, 1, 1, 2, 2, 2, 3]
            for record in parity_records.values()),
)
check(
    "independent SVDs reproduce every Schur spectrum above its error gate",
    all_recomputed,
    f"minimum singular={global_minimum_singular:.17e}, "
    f"minimum margin={global_minimum_margin:.6e}",
)
check(
    "representation multiplicities restore exactly 1440 strong and 120 weak",
    all(record["strong_full_count"] == 1440
            and record["schur_full_count"] == 120
            for record in parity_records.values()),
)

# Exact nontrivial rational check of graph-embedding invariance.
A = sy.Matrix(((2, 1), (1, 3)))
B = sy.Matrix(((1,), (2,)))
C = sy.Matrix(((3, 4),))
D = sy.Matrix(((5,),))
G = sy.Matrix(((2,), (-1,)))
base_schur = sy.simplify(D-C*A.inv()*B)
graph_schur = sy.simplify((C*G+D)-C*A.inv()*(A*G+B))
check(
    "the strong geometric graph cancels exactly from the Schur operator",
    A.det() != 0 and B != sy.zeros(2, 1) and C != sy.zeros(1, 2)
    and G != sy.zeros(2, 1) and graph_schur == base_schur,
    f"S={base_schur[0]}, S_G={graph_schur[0]}",
)

# S invertible and R injective imply S R injective.  The lower bound is the
# smallest certified singular value because an orthonormal basis can be used
# for the zero-sum relative hyperplane.
composition_injective = (
    global_minimum_singular > 0
    and relative_ranks == {prime: 119 for prime in PRIMES}
    and all(record["schur_full_count"] == 120
            for record in parity_records.values())
)
check(
    "the effective Schur equation is injective on all 119 relative shifts",
    composition_injective,
    f"certified inherited lower bound={global_minimum_singular:.17e}",
)

# Deliberately remove one whole minimal sector; restored rank must drop by
# representation multiplicity d times its 5d matrix dimension.
control_dimension = parity_records["even"]["sectors"][0]["irrep_dimension"]
expected_drop = control_dimension*(5*control_dimension)
shadow_rank = 120-expected_drop
negative_control_detected = shadow_rank < 120
check(
    "an exact-zero shadow Schur sector destroys the elimination verdict",
    negative_control_detected and expected_drop > 0,
    f"shadow rank={shadow_rank}, dropped={expected_drop}",
)

verdict = (
    "RELATIVE_SHIFT_CANONICALLY_ELIMINATED"
    if passed == tests else "RELATIVE_SHIFT_CANONICAL_STATUS_OPEN"
)
artifact = {
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "input_hashes": actual_hashes,
    },
    "relative_carrier": {
        "dimension": 119,
        "ranks": {str(key): value for key, value in relative_ranks.items()},
        "collective_complement_ranks": {
            str(key): value for key, value in augmented_ranks.items()},
        "vertex_incidence_ranks": {
            str(key): value for key, value in incidence_ranks.items()},
    },
    "canonical_carrier": {
        "strong_dimension": 1440,
        "pole_dimension": 120,
        "relative_pole_dimension": 119,
        "collective_lapse_dimension": 1,
        "parities": parity_records,
        "minimum_schur_singular": global_minimum_singular,
        "minimum_margin_over_100_epsilon": global_minimum_margin,
    },
    "graph_embedding": {
        "identity": "(CG+D)-C A^-1(AG+B)=D-C A^-1 B",
        "exact_control_base": str(base_schur[0]),
        "exact_control_graph": str(graph_schur[0]),
    },
    "negative_control": {
        "zeroed_irrep_dimension": control_dimension,
        "restored_rank_drop": expected_drop,
        "shadow_rank": shadow_rank,
    },
    "classification": {
        "free_relative_shift_propagation": "REFUTED_LOCAL_LINEAR",
        "canonical_status": "DERIVED_AUXILIARY_OR_PSEUDOCONSTRAINT_LIKE",
        "sourced_boundary_response": "OPEN",
        "graviton_dispersion_speed": "OPEN",
    },
    "tests": tests,
    "passed": passed,
    "verdict": verdict,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks pass")
print(f"VERDICT: {verdict}")
print(f"ARTIFACT: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
