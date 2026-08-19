#!/usr/bin/env python3
"""Exact local rigidity test for an anisotropic cellular tetrahedral frustum."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_cellular_frustum_anisotropic_rigidity.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_anisotropic_rigidity_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_anisotropic_rigidity_protocol.md"
BALANCED_RESULT = ROOT / "docs/gravity/gravity_600cell_projected_rank_edgewise_balanced_slab_result.md"
BALANCED_SOURCE = HERE / "verify_gravity_600cell_projected_rank_edgewise_balanced_slab.py"
BALANCED_JSON = HERE / "gravity_600cell_projected_rank_edgewise_balanced_slab.json"
TOURNAMENT_SOURCE = HERE / "verify_gravity_600cell_projected_rank_edgewise_balanced_slab_adversarial.py"
TOURNAMENT_JSON = HERE / "gravity_600cell_projected_rank_edgewise_balanced_slab_adversarial.json"
CELLULAR_SOURCE = HERE / "verify_gravity_600cell_projected_refinement_acceleration_blind.py"
CELLULAR_JSON = HERE / "gravity_600cell_projected_refinement_acceleration_blind.json"

PROTOCOL_COMMIT = "adc7243"
EXPECTED_HASHES = {
    "prior_art": "92c88042e8233a542b9f21e96a99bc0d09cf13cff89a8e243354f97984baaaab",
    "protocol": "ac20410bced8408c9cc8ec609653c3036a029b8e1d439a84c3acc3d5960eb1e8",
    "balanced_result": "6cd15954c73129fad4ac5905bdbe4440e9ef2a748b8a41c712662a34da3599bc",
    "balanced_source": "f59b8fc89106b42077eca281ff3d956a5a5d6fb4be70b73465133035b1ce0f57",
    "balanced_json": "0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5",
    "tournament_source": "794797d15f37887eb09c9aa168db73286705c47e1421c38e58a71d607099f2a8",
    "tournament_json": "dd1043a8cb712adb4f0717f95024b9ce62132501198938bb997e7ab3dad8bf65",
    "cellular_source": "e88111adaeb333abf80b68e06e23d7840ef14399238ada9d0f3cd722d7934e50",
    "cellular_json": "640bc0dd3d6f1ae727f8113bf29514878874effffd14f539f5a43e3c3b18d069",
}
METRICS = {
    "lorentz": sp.diag(1, 1, 1, -1),
    "euclidean": sp.eye(4),
}
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
BOTTOM = (
    (1, 1, 1, 0),
    (1, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, -1, 1, 0),
)
BOTTOM_EDGES = tuple(combinations(range(4), 2))
TOP_EDGES = tuple((left + 4, right + 4) for left, right in BOTTOM_EDGES)
STRUTS = tuple((index, index + 4) for index in range(4))
CELLULAR_EDGES = BOTTOM_EDGES + TOP_EDGES + STRUTS
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


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def points_for(scale, lapse):
    bottom = tuple(sp.Matrix(point) for point in BOTTOM)
    top = tuple(sp.Matrix((
        scale * point[0], scale * point[1], scale * point[2], lapse
    )) for point in BOTTOM)
    return bottom + top


def squared_length(points, edge, metric):
    delta = points[edge[0]] - points[edge[1]]
    return sp.expand((delta.T * metric * delta)[0])


def rigidity(points, edges, metric):
    matrix = sp.zeros(len(edges), 4 * len(points))
    for row, (left, right) in enumerate(edges):
        gradient = 2 * metric * (points[left] - points[right])
        for axis in range(4):
            matrix[row, 4 * left + axis] = gradient[axis]
            matrix[row, 4 * right + axis] = -gradient[axis]
    return matrix


def lorentz_generators(points, metric):
    columns = []
    for axis in range(4):
        column = sp.zeros(4 * len(points), 1)
        for vertex in range(len(points)):
            column[4 * vertex + axis] = 1
        columns.append(column)
    for left, right in combinations(range(4), 2):
        generator = sp.zeros(4, 4)
        generator[left, right] = 1
        generator[right, left] = -metric[left, left] / metric[right, right]
        if generator.T * metric + metric * generator != sp.zeros(4, 4):
            raise RuntimeError("invalid metric-isometry generator")
        columns.append(sp.Matrix.vstack(*(generator * point for point in points)))
    return sp.Matrix.hstack(*columns)


def diagonals(order):
    position = {colour: index for index, colour in enumerate(order)}
    result = []
    for left, right in BOTTOM_EDGES:
        if position[left] < position[right]:
            result.append((left, right + 4))
        else:
            result.append((right, left + 4))
    return tuple(sorted(tuple(sorted(edge)) for edge in result))


def reverse_edge(edge):
    def reverse_vertex(vertex):
        return vertex + 4 if vertex < 4 else vertex - 4
    return tuple(sorted((reverse_vertex(edge[0]), reverse_vertex(edge[1]))))


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "balanced_result": BALANCED_RESULT,
    "balanced_source": BALANCED_SOURCE,
    "balanced_json": BALANCED_JSON,
    "tournament_source": TOURNAMENT_SOURCE,
    "tournament_json": TOURNAMENT_JSON,
    "cellular_source": CELLULAR_SOURCE,
    "cellular_json": CELLULAR_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all exact frustum-rigidity inputs have frozen provenance",
      provenance_ok, str(hashes))

balanced = json.loads(BALANCED_JSON.read_text())
tournament = json.loads(TOURNAMENT_JSON.read_text())
cellular = json.loads(CELLULAR_JSON.read_text())
upstream_ok = bool(
    balanced["tests"] == {"passed": 15, "total": 15}
    and balanced["selection"] == {
        "canonical_selection_passes": False,
        "classification": "STRUCTURAL",
        "existence_passes": True,
        "h4_invariant_ordered_slab_alternatives": 24,
        "ordered_slab_alternatives": 24,
        "time_reversal_fixed_orders": 0,
        "time_reversal_orbit_sizes": {"2": 12},
    }
    and tournament["tests"] == {"passed": 13, "total": 13}
    and cellular["outcome"] == "PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_DERIVED"
    and cellular["passed"] == cellular["tests"] == 9
)
check("the 24-fold slab ambiguity and homogeneous cellular control persist",
      upstream_ok)

orders = tuple(permutations(range(4)))
diagonal_sets = {order: diagonals(order) for order in orders}
distinct_diagonals = set(diagonal_sets.values())
reversal_ok = True
reversal_pairs = set()
for order, edge_set in diagonal_sets.items():
    reversed_edges = tuple(sorted(reverse_edge(edge) for edge in edge_set))
    reversed_order = order[::-1]
    reversal_ok &= reversed_edges == diagonal_sets[reversed_order]
    reversal_pairs.add(frozenset((order, reversed_order)))
schedule_control = bool(
    len(orders) == len(distinct_diagonals) == 24
    and all(len(edge_set) == 6 and len(set(edge_set)) == 6
            for edge_set in diagonal_sets.values())
    and reversal_ok and len(reversal_pairs) == 12
    and all(len(pair) == 2 for pair in reversal_pairs)
)
check("the six diagonals encode 24 schedules and twelve reversal pairs",
      schedule_control)

records = []
representative_control = True
cellular_rank_ok = True
isometry_ok = True
fixed_flex_ok = True
signature_ok = True
completion_ok = True
restricted_completion_ok = True

for scale, lapse in REPRESENTATIVES:
    points = points_for(scale, lapse)
    affine = sp.Matrix.hstack(*(point - points[0] for point in points[1:]))
    bottom_gram = sp.Matrix.hstack(*(points[i] - points[0]
                                     for i in range(1, 4)))
    top_gram = sp.Matrix.hstack(*(points[i + 4] - points[4]
                                  for i in range(1, 4)))
    strut_squares = tuple(
        squared_length(points, edge, METRICS["lorentz"]) for edge in STRUTS
    )
    local_control = bool(
        affine.rank() == 4
        and bottom_gram.rank() == top_gram.rank() == 3
        and len(set(strut_squares)) == 1
        and strut_squares[0] < 0
    )
    representative_control &= local_control

    lorentz = rigidity(points, CELLULAR_EDGES, METRICS["lorentz"])
    euclidean = rigidity(points, CELLULAR_EDGES, METRICS["euclidean"])
    rank_lorentz = lorentz.rank()
    rank_euclidean = euclidean.rank()
    kernel_dimension = 32 - rank_lorentz
    generators = lorentz_generators(points, METRICS["lorentz"])
    generator_rank = generators.rank()
    annihilation_rank = (lorentz * generators).rank()
    quotient_flex = kernel_dimension - generator_rank
    cellular_rank_ok &= bool(
        rank_lorentz == 16 and kernel_dimension == 16 and quotient_flex == 6
    )
    signature_ok &= rank_euclidean == rank_lorentz
    isometry_ok &= generator_rank == 10 and annihilation_rank == 0

    fixed_edges = TOP_EDGES + STRUTS
    fixed_full = rigidity(points, fixed_edges, METRICS["lorentz"])
    fixed = fixed_full[:, 16:]
    fixed_rank = fixed.rank()
    fixed_nullspace = fixed.nullspace()
    fixed_kernel = (
        sp.Matrix.hstack(*fixed_nullspace)
        if fixed_nullspace else sp.zeros(16, 0)
    )
    fixed_flex_ok &= bool(
        fixed.shape == (10, 16) and fixed_rank == 10
        and fixed_kernel.shape == (16, 6)
    )

    completion_ranks = Counter()
    fixed_completion_ranks = Counter()
    flex_action_ranks = Counter()
    euclidean_completion_ranks = Counter()
    for order in orders:
        added = diagonal_sets[order]
        full_edges = CELLULAR_EDGES + added
        full_lorentz = rigidity(points, full_edges, METRICS["lorentz"])
        full_euclidean = rigidity(points, full_edges, METRICS["euclidean"])
        diagonal_full = rigidity(points, added, METRICS["lorentz"])
        diagonal_top = diagonal_full[:, 16:]
        full_fixed = fixed.col_join(diagonal_top)
        completion_ranks[full_lorentz.rank()] += 1
        euclidean_completion_ranks[full_euclidean.rank()] += 1
        fixed_completion_ranks[full_fixed.rank()] += 1
        flex_action_ranks[(diagonal_top * fixed_kernel).rank()] += 1

    completion_ok &= bool(
        completion_ranks == {22: 24}
        and fixed_completion_ranks == {16: 24}
        and euclidean_completion_ranks == {22: 24}
    )
    restricted_completion_ok &= flex_action_ranks == {6: 24}
    records.append({
        "scale": scale,
        "lapse": lapse,
        "strut_squared_length": str(strut_squares[0]),
        "affine_dimension": affine.rank(),
        "cellular_edges": len(CELLULAR_EDGES),
        "cellular_lorentz_rank": rank_lorentz,
        "cellular_euclidean_rank": rank_euclidean,
        "full_coordinate_kernel_dimension": kernel_dimension,
        "isometry_generator_rank": generator_rank,
        "isometry_annihilation_rank": annihilation_rank,
        "nonisometric_flex_dimension": quotient_flex,
        "fixed_bottom_constraint_rank": fixed_rank,
        "fixed_bottom_flex_dimension": fixed_kernel.shape[1],
        "completion_rank_counts": {str(k): v for k, v in completion_ranks.items()},
        "euclidean_completion_rank_counts": {
            str(k): v for k, v in euclidean_completion_ranks.items()
        },
        "fixed_completion_rank_counts": {
            str(k): v for k, v in fixed_completion_ranks.items()
        },
        "diagonal_action_on_flex_rank_counts": {
            str(k): v for k, v in flex_action_ranks.items()
        },
    })

check("all three rational frusta are nondegenerate with timelike struts",
      representative_control)
check("the cellular graph has exact rank 16 and six quotient flexes",
      cellular_rank_ok)
check("all ten Lorentz isometries are independent exact null directions",
      isometry_ok)
check("fixing the bottom leaves exactly six top-frustum flexes",
      fixed_flex_ok)
check("Euclidean and Lorentz signatures give the same exact ranks",
      signature_ok)
check("every six-diagonal staircase completion has full exact rank",
      completion_ok)
check("every diagonal set acts with rank six on the flex kernel",
      restricted_completion_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and schedule_control
    and representative_control and isometry_ok and signature_ok
)
exact_six = bool(
    controls_ok and cellular_rank_ok and fixed_flex_ok
    and completion_ok and restricted_completion_ok
)
worse = bool(
    controls_ok and (
        any(record["nonisometric_flex_dimension"] > 6 for record in records)
        or not completion_ok or not restricted_completion_ok
    )
)
if not controls_ok:
    outcome = "CELLULAR_FRUSTUM_RIGIDITY_CONTROL_FAILED"
elif exact_six:
    outcome = "CELLULAR_FRUSTUM_SIX_SHAPES_UNDERDETERMINED"
elif worse:
    outcome = "CELLULAR_FRUSTUM_UNDERDETERMINATION_WORSE"
else:
    outcome = "CELLULAR_FRUSTUM_RIGIDITY_OPEN"
allowed = {
    "CELLULAR_FRUSTUM_RIGIDITY_CONTROL_FAILED",
    "CELLULAR_FRUSTUM_SIX_SHAPES_UNDERDETERMINED",
    "CELLULAR_FRUSTUM_UNDERDETERMINATION_WORSE",
    "CELLULAR_FRUSTUM_RIGIDITY_OPEN",
}
check("the preregistered exact hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "metric_signatures": {
        name: [int(metric[i, i]) for i in range(4)]
        for name, metric in METRICS.items()
    },
    "cellular_graph": {
        "vertices": 8,
        "edges": len(CELLULAR_EDGES),
        "coordinate_dimension": 32,
        "isometry_dimension": 10,
        "geometric_dimension_mod_isometry": 22,
    },
    "schedule_count": len(distinct_diagonals),
    "time_reversal_orbits": len(reversal_pairs),
    "records": records,
    "classification": {
        "homogeneous_cellular_action": "PRESERVED DERIVED COMPUTATIONAL",
        "anisotropic_metric_from_16_lengths": (
            "REFUTED DERIVED EXACT / STRUCTURAL" if exact_six else "OPEN"
        ),
        "six_staircase_diagonals": (
            "COMPLETE THE LOCAL METRIC, BUT 24 CHOICES REMAIN"
            if exact_six else "OPEN"
        ),
        "refined_anisotropic_hessian": (
            "NOT AUTHORIZED WITHOUT NEW SHAPE DATA"
            if exact_six else "OPEN"
        ),
        "ensemble_area_angle_or_physics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"rank={record['cellular_lorentz_rank']}, "
        f"nonisometric flexes={record['nonisometric_flex_dimension']}, "
        f"completions={record['completion_rank_counts']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
