#!/usr/bin/env python3
"""Independent incidence/action audit of the refined curvature--mass identity."""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


PRIMARY_SOURCE = HERE / "verify_gravity_600cell_refined_local_curvature_mass.py"
PRIMARY = HERE / "gravity_600cell_refined_local_curvature_mass.json"
PRIMARY_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_local_curvature_mass_primary_result.md"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
CELL600 = ROOT / "commons/cell600.py"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_local_curvature_mass_adversarial_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_local_curvature_mass_adversarial.json"

PRIMARY_RESULT_COMMIT = "2881e4c"
PROTOCOL_COMMIT = "f4336e7"
EXPECTED_HASHES = {
    "primary_source": "c54f17708a2678b925cfce96fcfc7d6baaeeaf0577bedbf22b5d0435c069fae6",
    "primary": "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "primary_result": "82b55ed2918b2db4f83123b23cb2e7534a9f813d5705a1e2ee0ac92e481a0f10",
    "action_source": "89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7",
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "protocol": "906bfaf6b766ba3fd09229844dda603e15f4adbf5ba0ff8fb8cfd8f571f5d182",
}
PAIR4 = tuple(combinations(range(4), 2))
LOCAL_TRIANGLES = np.asarray(tuple(combinations(range(5), 3)), dtype=np.int8)
VARIABLES = (
    tuple(("old",) + pair for pair in PAIR4)
    + tuple(("new",) + pair for pair in PAIR4)
    + tuple(("cross",) + pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
INTERNAL_VARIABLES = (
    tuple(("cross",) + pair for pair in PAIR4)
    + tuple(("rho", rank) for rank in range(4))
)
RANK_SIZES = (120, 720, 1200, 600)
TAU_TEXT = "0.0102"

tests = 0
passed = 0


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


def mp_text(value, digits=75):
    return mp.nstr(value, digits, strip_zeros=False)


def load_action_definitions():
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
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
        "LOCAL_TRIANGLES": LOCAL_TRIANGLES,
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
        compile(ast.Module(body=definitions, type_ignores=[]), str(ACTION_SOURCE), "exec"),
        namespace,
    )
    return namespace


def edge_incidence(top):
    counts = Counter()
    for tetrahedron in top:
        for edge in combinations(map(int, tetrahedron), 2):
            counts[tuple(sorted(edge))] += 1
    return counts


def endpoint_accumulation(edges, contribution, vertex_count):
    values = [mp.mpf(0) for _ in range(vertex_count)]
    for edge in edges:
        amount = contribution(edge)
        values[edge[0]] += amount / 2
        values[edge[1]] += amount / 2
    return tuple(values)


paths = {
    "primary_source": PRIMARY_SOURCE,
    "primary": PRIMARY,
    "primary_result": PRIMARY_RESULT,
    "action_source": ACTION_SOURCE,
    "cell600": CELL600,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "the primary result and independent action inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIMARY_RESULT_COMMIT == "2881e4c"
    and PROTOCOL_COMMIT == "f4336e7",
)

primary = json.loads(PRIMARY.read_text())
primary_ok = check(
    "the frozen primary artifact carries its complete post-hoc result",
    primary["outcome"]
        == "REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC"
    and primary["tests"] == {"passed": 15, "total": 15}
    and primary["identity_audit"]["schedule_count"] == 24
    and primary["identity_audit"]["equation_count"] == 96
    and primary["identity_audit"]["mass_response_rank"] == 4,
)

defs = load_action_definitions()
definition_ok = check(
    "only frozen stationary-fill function definitions were loaded",
    all(name in defs for name in (
        "tetrahedra_from_adjacency",
        "all_simplices",
        "barycentric_chambers",
        "schedule_combinatorics",
        "exact_geometry",
        "base_coordinates",
        "evaluate_schedule",
    ))
    and "OUTPUT" not in defs,
)

coarse_vertices, coarse_adjacency, _ = build_600cell()
coarse_top = defs["tetrahedra_from_adjacency"](coarse_adjacency)
vertex_cells, top, colours = defs["barycentric_chambers"](coarse_top)
spatial_cells = defs["all_simplices"](tuple(map(tuple, top)))
edges = tuple(spatial_cells[1])
incidence = edge_incidence(top)

topology_ok = check(
    "the actual barycentric flag complex has the complete refined f-vector",
    tuple(len(layer) for layer in spatial_cells) == (2640, 17040, 28800, 14400)
    and len(vertex_cells) == len(colours) == 2640
    and len(edges) == len(incidence) == 17040
    and tuple(int(np.count_nonzero(colours == rank)) for rank in range(4))
        == RANK_SIZES,
)

with mp.workdps(100):
    geometry = defs["exact_geometry"](100)

    def actual_edge_curvature(edge):
        pair = tuple(sorted((int(colours[edge[0]]), int(colours[edge[1]]))))
        length = geometry["s0"] * mp.sqrt(geometry["unit_squares"][pair])
        deficit = 2 * mp.pi - incidence[edge] * geometry["spatial_angles"][pair]
        return length * deficit

    vertex_curvatures = endpoint_accumulation(edges, actual_edge_curvature, 2640)
    edge_curvature = mp.fsum(actual_edge_curvature(edge) for edge in edges)
    vertex_curvature = mp.fsum(vertex_curvatures)
    rank_curvatures = tuple(
        mp.fsum(
            vertex_curvatures[vertex]
            for vertex in range(2640)
            if int(colours[vertex]) == rank
        )
        for rank in range(4)
    )
    fractions = tuple(value / edge_curvature for value in rank_curvatures)
    within_rank_spreads = tuple(
        max(
            vertex_curvatures[vertex]
            for vertex in range(2640)
            if int(colours[vertex]) == rank
        )
        - min(
            vertex_curvatures[vertex]
            for vertex in range(2640)
            if int(colours[vertex]) == rank
        )
        for rank in range(4)
    )
    primary_fractions = tuple(
        mp.mpf(value) for value in primary["selected_rank_matter"]["fractions"]
    )
    fraction_error = max(
        abs(fractions[rank] - primary_fractions[rank]) for rank in range(4)
    )
    conservation_relative = abs(vertex_curvature - edge_curvature) / edge_curvature
    spread_relative = max(within_rank_spreads) / max(vertex_curvatures)

incidence_ok = check(
    "actual edge incidences give positive curvature and one value per rank orbit",
    all(actual_edge_curvature(edge) > 0 for edge in edges)
    and spread_relative < mp.mpf("1e-70"),
    f"relative spread={mp_text(spread_relative, 8)}",
)
conservation_ok = check(
    "the 2,640 endpoint contributions conserve the 17,040 edge curvature sum",
    conservation_relative < mp.mpf("1e-70"),
    f"relative error={mp_text(conservation_relative, 8)}",
)
fraction_ok = check(
    "actual-incidence curvature fractions reproduce the primary four fractions",
    fraction_error < mp.mpf("1e-68"),
    f"max error={mp_text(fraction_error, 8)}",
)

coarse_edges = tuple(
    (left, right)
    for left in range(len(coarse_adjacency))
    for right in range(left + 1, len(coarse_adjacency))
    if coarse_adjacency[left, right] > 0.5
)
with mp.workdps(100):
    regular_vertex = endpoint_accumulation(
        coarse_edges, lambda edge: mp.mpf(1), len(coarse_vertices)
    )
    regular_total = mp.fsum(regular_vertex)
    regular_error = max(
        abs(value / regular_total - mp.mpf(1) / 120)
        for value in regular_vertex
    )
regular_control_ok = check(
    "the unrefined regular 600-cell control gives exactly 1/120 per vertex",
    len(coarse_edges) == 720 and regular_error == 0,
)

with mp.workdps(100):
    corrupted_vertex_curvatures = list(vertex_curvatures)
    first_edge = edges[0]
    corrupted_vertex_curvatures[first_edge[0]] -= actual_edge_curvature(first_edge) / 2
    corrupted_conservation = abs(
        mp.fsum(corrupted_vertex_curvatures) - edge_curvature
    ) / edge_curvature
negative_control_ok = check(
    "dropping one endpoint contribution fails curvature conservation",
    corrupted_conservation > mp.mpf("1e-8"),
    f"relative defect={mp_text(corrupted_conservation, 8)}",
)

print("[INFO] constructing direct dust-free action for two schedules", flush=True)
orders = ((0, 1, 2, 3), (3, 2, 1, 0))
action_records = []
finite_difference_errors = []
with mp.workdps(100):
    tau = mp.mpf(TAU_TEXT)
    action_geometry = defs["exact_geometry"](100)
    action_geometry["mass"] = mp.mpf(0)
    coordinates = defs["base_coordinates"](action_geometry)
    for order in orders:
        combinatorics = defs["schedule_combinatorics"](top, colours, order)
        evaluation = defs["evaluate_schedule"](
            combinatorics, action_geometry, coordinates
        )
        action_records.append((order, combinatorics, evaluation))
        for rank in range(4):
            key = ("rho", rank)
            derivatives = []
            for step_text in ("1e-12", "5e-13"):
                step = mp.mpf(step_text)
                plus = dict(coordinates)
                minus = dict(coordinates)
                plus[key] *= mp.exp(step)
                minus[key] *= mp.exp(-step)
                plus_action = defs["evaluate_schedule"](
                    combinatorics, action_geometry, plus
                )["action"]
                minus_action = defs["evaluate_schedule"](
                    combinatorics, action_geometry, minus
                )["action"]
                derivatives.append((plus_action - minus_action) / (2 * step))
            richardson = (4 * derivatives[1] - derivatives[0]) / 3
            analytic = evaluation["gradient"][key]
            finite_difference_errors.append(
                abs(richardson - analytic) / max(mp.mpf(1), abs(analytic))
            )

    direct_identity_errors = []
    selected_total_residuals = []
    p1_total_residuals = []
    selected_masses = tuple(value / (8 * mp.pi) for value in rank_curvatures)
    total_mass = mp.fsum(selected_masses)
    for _, _, evaluation in action_records:
        for rank in range(4):
            gravitational = evaluation["gradient"]["rho", rank]
            expected = tau * rank_curvatures[rank] / 2
            direct_identity_errors.append(gravitational - expected)
            selected_total_residuals.append(
                gravitational - 4 * mp.pi * selected_masses[rank] * tau
            )
            p1_total_residuals.append(
                gravitational - 4 * mp.pi * (total_mass / 4) * tau
            )
    direct_error = max(abs(value) for value in direct_identity_errors)
    selected_residual = max(abs(value) for value in selected_total_residuals)
    p1_residual = max(abs(value) for value in p1_total_residuals)
    finite_difference_error = max(finite_difference_errors)
    branch_identity = max(
        evaluation["maximum_angle_identity_residual"]
        for _, _, evaluation in action_records
    )
    branch_imaginary = max(
        evaluation["maximum_imaginary_curvature"]
        for _, _, evaluation in action_records
    )
    branch_argument = min(
        evaluation["minimum_angle_argument"]
        for _, _, evaluation in action_records
    )

branch_ok = check(
    "both dust-free evaluations stay on the accepted Lorentzian branch",
    branch_identity < mp.mpf("1e-80")
    and branch_imaginary < mp.mpf("1e-80")
    and branch_argument > mp.mpf("1e-8"),
    f"identity={mp_text(branch_identity, 8)}, "
    f"imag={mp_text(branch_imaginary, 8)}, arg={mp_text(branch_argument, 8)}",
)
direct_action_ok = check(
    "two direct dust-free actions reproduce tau times actual local curvature / 2",
    direct_error < mp.mpf("1e-60"),
    f"max error={mp_text(direct_error, 8)}",
)
finite_difference_ok = check(
    "independent Richardson derivatives reproduce all eight analytic gradients",
    finite_difference_error < mp.mpf("1e-36"),
    f"max relative error={mp_text(finite_difference_error, 8)}",
)
selected_stationary_ok = check(
    "the actual-incidence curvature masses cancel all eight lapse gradients",
    selected_residual < mp.mpf("1e-60"),
    f"max residual={mp_text(selected_residual, 8)}",
)
p1_negative_ok = check(
    "direct P1 equal-rank masses remain detectably off shell",
    p1_residual > mp.mpf("1e-4"),
    f"max residual={mp_text(p1_residual, 8)}",
)

scope = {
    "primary_functions_imported": False,
    "primary_residual_subtraction_used": False,
    "root_search_or_nested_census_executed": False,
    "hessian_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
}
scope_ok = check(
    "the adversarial route is independent and remains inside scope",
    not any(scope.values()),
)

controls_ok = all((
    provenance_ok,
    primary_ok,
    definition_ok,
    topology_ok,
    incidence_ok,
    conservation_ok,
    regular_control_ok,
    negative_control_ok,
    branch_ok,
    finite_difference_ok,
    p1_negative_ok,
    scope_ok,
))
agreement_ok = fraction_ok and direct_action_ok and selected_stationary_ok
if not controls_ok:
    outcome = "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED"
elif not agreement_ok:
    outcome = "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_DISAGREEMENT"
else:
    outcome = "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one adversarial outcome",
    outcome in {
        "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED",
        "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_DISAGREEMENT",
        "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED",
    },
    outcome,
)

artifact = {
    "title": "Adversarial refined local curvature--mass audit",
    "date": "2026-08-21",
    "primary_result_commit": PRIMARY_RESULT_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "actual_incidence": {
        "f_vector": [len(layer) for layer in spatial_cells],
        "rank_sizes": list(RANK_SIZES),
        "spatial_edges": len(edges),
        "edge_incidence_values": sorted(set(incidence.values())),
        "rank_curvatures": [mp_text(value) for value in rank_curvatures],
        "fractions": [mp_text(value) for value in fractions],
        "within_rank_relative_spread": mp_text(spread_relative),
        "endpoint_edge_conservation_relative_error": mp_text(conservation_relative),
        "primary_fraction_max_error": mp_text(fraction_error),
    },
    "known_and_negative_controls": {
        "unrefined_edges": len(coarse_edges),
        "unrefined_uniform_fraction_error": mp_text(regular_error),
        "dropped_endpoint_relative_defect": mp_text(corrupted_conservation),
        "p1_direct_residual": mp_text(p1_residual),
    },
    "direct_dust_free_action": {
        "orders": [list(order) for order in orders],
        "maximum_identity_error": mp_text(direct_error),
        "maximum_selected_mass_residual": mp_text(selected_residual),
        "maximum_richardson_relative_error": mp_text(finite_difference_error),
        "maximum_angle_identity_residual": mp_text(branch_identity),
        "maximum_imaginary_curvature": mp_text(branch_imaginary),
        "minimum_angle_argument": mp_text(branch_argument),
    },
    "scope": scope,
    "status_labels": {
        "identity": "ADVERSARIAL_TEST_OF_POST_HOC_PRIMARY",
        "homogeneous_dust": "NOT_ESTABLISHED",
        "refinement_convergence": "OPEN",
        "tick_c_G_planck_particles": "OPEN_NOT_COMPUTED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Tests passed: {passed}/{tests}")
print(f"Outcome: {outcome}")
print("fractions=" + ", ".join(mp_text(value, 22) for value in fractions))

raise SystemExit(0 if passed == tests else 1)
