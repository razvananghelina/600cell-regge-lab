#!/usr/bin/env python3
"""Exact homogeneous coarse-to-refined cotangent-lift rank certificate."""

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FEASIBILITY = HERE / "gravity_600cell_refined_canonical_map_feasibility.json"
FEASIBILITY_RESULT = (
    ROOT / "docs/gravity/gravity_600cell_refined_canonical_map_feasibility_result.md"
)
PRIOR_ART = (
    ROOT
    / "docs/gravity/gravity_600cell_refined_homogeneous_cotangent_lift_prior_art.md"
)
PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_refined_homogeneous_cotangent_lift_protocol.md"
)
OUTPUT = HERE / "gravity_600cell_refined_homogeneous_cotangent_lift.json"

PRIOR_ART_COMMIT = "3e188e7"
PROTOCOL_COMMIT = "3aefd29"
EXPECTED_HASHES = {
    "feasibility": "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
    "feasibility_result": "222739a82680c35c337b127c455fc0a8a7c24c2bd4a2f6c8f9953eeb3251e681",
    "prior_art": "fa391b36501f08dfee9a0bd588f651e772cebc38643acd4c852fc95bb8cd6f21",
    "protocol": "0ed274170526e1e84901e5425f2773f372bbc04cfe1072a50ec1452ea0688fc7",
}
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EXPECTED_POPULATIONS = (1440, 3600, 2400, 3600, 3600, 2400)

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def matrix_rows(matrix):
    return [[str(matrix[row, column]) for column in range(matrix.cols)]
            for row in range(matrix.rows)]


def exact_rank_geometry():
    phi = (1 + sp.sqrt(5)) / 2
    c_value = phi / 2

    def norm_square(size):
        return size * (1 + (size - 1) * c_value)

    squares = []
    positivity_gaps = []
    for left, right in PAIRS:
        a, b = left + 1, right + 1
        dot = sp.simplify(
            a * (1 + (b - 1) * c_value)
            / sp.sqrt(norm_square(a) * norm_square(b))
        )
        squares.append(sp.simplify(2 - 2 * dot))
        positivity_gaps.append(sp.simplify(1 - dot**2))
    return tuple(squares), tuple(positivity_gaps)


paths = {
    "feasibility": FEASIBILITY,
    "feasibility_result": FEASIBILITY_RESULT,
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all frozen inputs match their preregistered hashes",
    actual_hashes == EXPECTED_HASHES,
)

feasibility = json.loads(FEASIBILITY.read_text())
level = feasibility["levels"]["projected_barycentric"]
upstream_ok = check(
    "the accepted refined feasibility census supplies K0",
    feasibility["outcome"] == "REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED"
    and level["spatial_f_vector"] == [2640, 17040, 28800, 14400]
    and level["distinct_internal_edge_sets"] == 24,
)

population_map = level["colour_pair_edge_populations"]
populations = tuple(population_map[f"{left}-{right}"] for left, right in PAIRS)
population_ok = check(
    "the six positive orbit populations and total are reproduced",
    populations == EXPECTED_POPULATIONS
    and sum(populations) == 17040
    and all(value > 0 for value in populations),
    f"N={populations}, sum={sum(populations)}",
)

unit_squares, positivity_gaps = exact_rank_geometry()
geometry_ok = check(
    "the exact golden-ratio rank geometry has six positive edge squares",
    all(value.is_positive is True for value in unit_squares)
    and all(value.is_positive is True for value in positivity_gaps),
    "squares=" + ", ".join(str(value) for value in unit_squares),
)

h_squared = sp.Matrix([[2] * 6])
d_population = sp.diag(*populations)
pullback_total = h_squared
pullback_edge = h_squared * d_population

rank_total = pullback_total.rank()
rank_edge = pullback_edge.rank()
null_total = 6 - rank_total
null_edge = 6 - rank_edge
total_ok = check(
    "the orbit-total cotangent pullback has exact rank one and nullity five",
    rank_total == 1 and null_total == 5,
    f"rank={rank_total}, nullity={null_total}",
)
edge_ok = check(
    "the per-edge cotangent pullback has exact rank one and nullity five",
    rank_edge == 1 and null_edge == 5,
    f"rank={rank_edge}, nullity={null_edge}",
)

convention_ok = check(
    "positive population scaling exactly intertwines the two conventions",
    d_population.det() != 0
    and pullback_edge == pullback_total * d_population,
    f"det(D)={d_population.det()}",
)

kernel_total = pullback_total.nullspace()
kernel_edge = pullback_edge.nullspace()
kernel_ok = check(
    "both exact kernel bases have five independent null vectors",
    len(kernel_total) == 5
    and len(kernel_edge) == 5
    and sp.Matrix.hstack(*kernel_total).rank() == 5
    and sp.Matrix.hstack(*kernel_edge).rank() == 5
    and all(pullback_total * vector == sp.zeros(1, 1) for vector in kernel_total)
    and all(pullback_edge * vector == sp.zeros(1, 1) for vector in kernel_edge),
)

synthetic = sp.Matrix([[2]])
synthetic_ok = check(
    "the one-orbit control has a unique lift",
    synthetic.rank() == 1
    and len(synthetic.nullspace()) == 0
    and synthetic.inv() * sp.Matrix([1]) == sp.Matrix([sp.Rational(1, 2)]),
)

length_coordinate = sp.Matrix([[1] * 6])
reversal = sp.zeros(6, 6)
for index in range(6):
    reversal[index, 5 - index] = 1
coordinate_ok = check(
    "coordinate rescaling and orbit reversal preserve rank and nullity",
    length_coordinate.rank() == rank_total
    and len(length_coordinate.nullspace()) == null_total
    and (pullback_total * reversal).rank() == rank_total
    and len((pullback_total * reversal).nullspace()) == null_total
    and (pullback_edge * reversal).rank() == rank_edge
    and len((pullback_edge * reversal).nullspace()) == null_edge,
)

singular_control = sp.diag(0, *populations[1:])
negative_control_ok = check(
    "a zero-population corruption destroys convention invertibility",
    singular_control.det() == 0,
)

controls_ok = all((
    provenance_ok, upstream_ok, population_ok, geometry_ok, convention_ok,
    kernel_ok, synthetic_ok, coordinate_ok, negative_control_ok,
))
if not controls_ok:
    outcome = "REFINED_HOMOGENEOUS_COTANGENT_CONTROL_FAILED"
elif null_total == 0 and null_edge == 0:
    outcome = "REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNIQUE"
elif rank_total == rank_edge == 1 and null_total == null_edge == 5:
    outcome = "REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNDERDETERMINED"
else:
    outcome = "REFINED_HOMOGENEOUS_COTANGENT_LIFT_OPEN"

outcome_ok = check(
    "the frozen hierarchy assigns the cotangent-lift outcome",
    outcome in {
        "REFINED_HOMOGENEOUS_COTANGENT_CONTROL_FAILED",
        "REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNIQUE",
        "REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNDERDETERMINED",
        "REFINED_HOMOGENEOUS_COTANGENT_LIFT_OPEN",
    },
    outcome,
)

artifact = {
    "title": "Homogeneous coarse-to-refined cotangent-lift rank certificate",
    "date": "2026-08-21",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "definitions": {
        "coarse_homogeneous_configuration_dimension": 1,
        "refined_h4_edge_orbits": [f"{left}{right}" for left, right in PAIRS],
        "refined_h4_orbit_count": 6,
        "coordinate": "log squared edge length",
        "homothetic_tangent": [2] * 6,
        "preferred_particular_lift_computed": False,
        "action_hessian_used": False,
        "refined_slab_solved": False,
    },
    "controls": {
        "populations": list(populations),
        "population_sum": sum(populations),
        "unit_edge_squares_exact": [str(value) for value in unit_squares],
        "unit_edge_squares_decimal": [sp.N(value, 30).__str__()
                                          for value in unit_squares],
        "one_orbit_rank": synthetic.rank(),
        "one_orbit_nullity": len(synthetic.nullspace()),
        "one_orbit_unit_momentum_lift": "1/2",
        "zero_population_control_singular": singular_control.det() == 0,
    },
    "orbit_total": {
        "pullback": matrix_rows(pullback_total),
        "rank": rank_total,
        "nullity": null_total,
        "kernel_basis": [matrix_rows(vector) for vector in kernel_total],
        "affine_free_parameters_for_fixed_coarse_momentum": null_total,
    },
    "per_edge": {
        "pullback": matrix_rows(pullback_edge),
        "rank": rank_edge,
        "nullity": null_edge,
        "kernel_basis": [matrix_rows(vector) for vector in kernel_edge],
        "affine_free_parameters_for_fixed_coarse_momentum": null_edge,
    },
    "convention_intertwiner": {
        "diagonal": list(populations),
        "determinant": str(d_population.det()),
        "invertible": d_population.det() != 0,
    },
    "interpretation": {
        "status": "STRUCTURAL" if "UNDERDETERMINED" in outcome else "OPEN",
        "excluded_claim": (
            "geometry, H4 invariance and canonical pairing alone select a "
            "unique refined homogeneous momentum"
        ),
        "not_excluded": (
            "an action-selected, supermetric-selected or perfect-action "
            "coarse-to-fine transport"
        ),
        "tick_c_G_planck": "NOT COMPUTED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
print(f"ARTIFACT: {OUTPUT}")
raise SystemExit(0 if passed == tests and outcome_ok else 1)
