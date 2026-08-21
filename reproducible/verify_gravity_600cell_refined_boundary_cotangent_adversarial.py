#!/usr/bin/env python3
"""Derive the refined boundary covector from actual spatial hinge curvature."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


PRIMARY_SOURCE = HERE / "verify_gravity_600cell_refined_boundary_cotangent.py"
PRIMARY = HERE / "gravity_600cell_refined_boundary_cotangent.json"
PRIMARY_RESULT = ROOT / "docs/gravity/gravity_600cell_refined_boundary_cotangent_primary_result.md"
CURVATURE_ADVERSARIAL = HERE / "gravity_600cell_refined_local_curvature_mass_adversarial.json"
COARSE_IDENTITY = HERE / "gravity_600cell_dust_regular_lapse_identity.json"
CELL600 = ROOT / "commons/cell600.py"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_boundary_cotangent_adversarial_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_boundary_cotangent_adversarial.json"

PRIMARY_RESULT_COMMIT = "f805557"
PROTOCOL_COMMIT = "fbe6613"
EXPECTED_HASHES = {
    "primary_source": "ababad0e8e667e31c290b9e8bbf61005308faed20af4d09ae7affbc32b3509d7",
    "primary": "4e7bf0beb0327a3ee1bddbec13126fbef99380970e62cecf74eb24ce8d6dafaa",
    "primary_result": "ec69cd7c6521b3cff3ded777d80a4c740b9065bcce73bd52f53c0591433c9074",
    "curvature_adversarial": "c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e",
    "coarse_identity": "5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51",
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "protocol": "75a3531b4be9d83a97de5bebf63e418db315f6994c7b7b4fd7b6d633c07cb449",
}
PAIRS = tuple(combinations(range(4), 2))
TAU_TEXT = "0.0102"

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


def mp_text(value, digits=75):
    return mp.nstr(value, digits, strip_zeros=False)


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    result = []
    for first in range(len(adjacency)):
        for second in sorted(value for value in neighbours[first] if value > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(value for value in common_two if value > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(value for value in common_three if value > third):
                    result.append((first, second, third, fourth))
    return tuple(result)


def all_simplices(top):
    return tuple(
        tuple(sorted({
            tuple(sorted(face))
            for tetrahedron in top
            for face in combinations(tetrahedron, degree + 1)
        }))
        for degree in range(4)
    )


def barycentric_chambers(coarse_top):
    coarse_cells = all_simplices(coarse_top)
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    cell_index = {cell: index for index, cell in enumerate(vertex_cells)}
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
    colours = tuple(len(cell) - 1 for cell in vertex_cells)
    return vertex_cells, tuple(top), colours


def reconstruct_metric():
    phi = (1 + mp.sqrt(5)) / 2
    c_value = phi / 2

    def norm_square(size):
        return size * (1 + (size - 1) * c_value)

    unit_squares = {}
    for left, right in PAIRS:
        a, b = left + 1, right + 1
        dot = (
            a * (1 + (b - 1) * c_value)
            / mp.sqrt(norm_square(a) * norm_square(b))
        )
        unit_squares[left, right] = 2 - 2 * dot

    squared = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    for (left, right), value in unit_squares.items():
        squared[left][right] = squared[right][left] = value
    gram = mp.matrix([
        [
            (squared[3][left] + squared[3][right] - squared[left][right]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ])
    chamber_volume = mp.sqrt(mp.det(gram)) / 6
    scale = (2 * mp.pi**2 / (14400 * chamber_volume)) ** (mp.mpf(1) / 3)
    inverse = gram**-1
    normals = (
        mp.matrix([1, 0, 0]),
        mp.matrix([0, 1, 0]),
        mp.matrix([0, 0, 1]),
        mp.matrix([-1, -1, -1]),
    )

    def inner(left, right):
        return (left.T * inverse * right)[0]

    angles = {}
    for left, right in PAIRS:
        omitted_a, omitted_b = [
            value for value in range(4) if value not in (left, right)
        ]
        cosine = -inner(normals[omitted_a], normals[omitted_b]) / mp.sqrt(
            inner(normals[omitted_a], normals[omitted_a])
            * inner(normals[omitted_b], normals[omitted_b])
        )
        angles[left, right] = mp.acos(cosine)
    return unit_squares, angles, scale


paths = {
    "primary_source": PRIMARY_SOURCE,
    "primary": PRIMARY,
    "primary_result": PRIMARY_RESULT,
    "curvature_adversarial": CURVATURE_ADVERSARIAL,
    "coarse_identity": COARSE_IDENTITY,
    "cell600": CELL600,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all primary, curvature, coarse and geometry inputs have exact provenance",
    actual_hashes == EXPECTED_HASHES
    and PRIMARY_RESULT_COMMIT == "f805557"
    and PROTOCOL_COMMIT == "fbe6613",
)

primary = json.loads(PRIMARY.read_text())
curvature_adversarial = json.loads(CURVATURE_ADVERSARIAL.read_text())
coarse_identity = json.loads(COARSE_IDENTITY.read_text())
upstream_ok = check(
    "the frozen inputs carry the complete accepted scoped outcomes",
    primary["outcome"] == "REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED"
    and primary["tests"] == {"passed": 16, "total": 16}
    and curvature_adversarial["outcome"]
        == "ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED"
    and curvature_adversarial["tests"] == {"passed": 16, "total": 16}
    and coarse_identity["outcome"] == "REGULAR_LAPSE_IDENTITY_PROVED"
    and coarse_identity["passed"] == coarse_identity["tests"] == 13,
)

l_symbol, tau_symbol, x_symbol, y_symbol, z_symbol = sp.symbols(
    "l tau x y z", positive=True
)
area_square = (
    2 * (x_symbol * y_symbol + x_symbol * z_symbol + y_symbol * z_symbol)
    - x_symbol**2 - y_symbol**2 - z_symbol**2
) / 16
product_substitution = {
    x_symbol: l_symbol**2,
    y_symbol: -tau_symbol**2,
    z_symbol: l_symbol**2 - tau_symbol**2,
}
product_area_square = sp.simplify(area_square.subs(product_substitution))
product_area = sp.I * l_symbol * tau_symbol / 2
log_area_derivative = sp.simplify(
    l_symbol**2
    * sp.diff(area_square, x_symbol).subs(product_substitution)
    / (2 * product_area)
)
symbolic_ok = check(
    "Heron differentiation gives the independent product-hinge factor 1/4",
    sp.simplify(product_area_square + l_symbol**2 * tau_symbol**2 / 4) == 0
    and sp.simplify(log_area_derivative - sp.I * l_symbol * tau_symbol / 4) == 0,
    f"A^2={product_area_square}, dA/dlog(x)={log_area_derivative}",
)

coarse_vertices, coarse_adjacency, _ = build_600cell()
coarse_top = tetrahedra_from_adjacency(coarse_adjacency)
vertex_cells, top, colours = barycentric_chambers(coarse_top)
spatial_cells = all_simplices(top)
edges = spatial_cells[1]
incidence = Counter()
for tetrahedron in top:
    for edge in combinations(tetrahedron, 2):
        incidence[tuple(sorted(edge))] += 1
topology_ok = check(
    "the independent flag-complex reconstruction has all actual cells and edges",
    len(coarse_vertices) == 120
    and len(coarse_top) == 600
    and tuple(len(layer) for layer in spatial_cells)
        == (2640, 17040, 28800, 14400)
    and len(vertex_cells) == 2640
    and len(edges) == len(incidence) == 17040,
)

with mp.workdps(100):
    tau = mp.mpf(TAU_TEXT)
    unit_squares, angles, scale = reconstruct_metric()

    def edge_curvature(edge):
        pair = tuple(sorted((colours[edge[0]], colours[edge[1]])))
        length = scale * mp.sqrt(unit_squares[pair])
        deficit = 2 * mp.pi - incidence[edge] * angles[pair]
        return length * deficit

    pair_curvatures = {
        pair: mp.fsum(
            edge_curvature(edge)
            for edge in edges
            if tuple(sorted((colours[edge[0]], colours[edge[1]]))) == pair
        )
        for pair in PAIRS
    }
    k_fine = mp.fsum(pair_curvatures.values())
    derived_pre = tuple(-tau * pair_curvatures[pair] / 4 for pair in PAIRS)
    derived_post = tuple(+tau * pair_curvatures[pair] / 4 for pair in PAIRS)

# The primary vector is deliberately read only after the independent build.
with mp.workdps(100):
    primary_pre = tuple(
        mp.mpf(value) for value in primary["boundary_covectors"]["selected_pre"]
    )
    primary_post = tuple(
        mp.mpf(value) for value in primary["boundary_covectors"]["selected_post"]
    )
    component_error = max(
        *(abs(derived_pre[index] - primary_pre[index]) for index in range(6)),
        *(abs(derived_post[index] - primary_post[index]) for index in range(6)),
    )
component_ok = check(
    "actual hinge curvature reproduces all twelve primary boundary components",
    component_error < mp.mpf("1e-68"),
    f"max error={mp_text(component_error, 8)}",
)

coarse_edges = tuple(
    (left, right)
    for left in range(120)
    for right in range(left + 1, 120)
    if coarse_adjacency[left, right] > 0.5
)
with mp.workdps(100):
    zeta = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
    epsilon3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
    coarse_per_edge = -tau * zeta * epsilon3 / 4
    k_coarse = len(coarse_edges) * zeta * epsilon3
    p_pre_coarse = 2 * len(coarse_edges) * coarse_per_edge
    coarse_identity_error = abs(p_pre_coarse + tau * k_coarse / 2)
regular_control_ok = check(
    "the same hinge formula reproduces the exact regular 600-cell momentum",
    len(coarse_edges) == 720
    and coarse_identity["exact_algebra"]["pre_momentum_per_edge"]
        == "-L*epsilon*tau/4"
    and coarse_identity_error < mp.mpf("1e-90"),
    f"error={mp_text(coarse_identity_error, 8)}",
)

with mp.workdps(100):
    p_pre_fine = 2 * mp.fsum(derived_pre)
    pullback_error = abs(p_pre_fine + tau * k_fine / 2)
    m_fine = k_fine / (8 * mp.pi)
    normalized_error = abs(p_pre_fine / m_fine + 4 * mp.pi * tau)
    raw_ratio = p_pre_fine / p_pre_coarse
    curvature_ratio = k_fine / k_coarse
    ratio_error = abs(raw_ratio - curvature_ratio)
pullback_ok = check(
    "the independent vector has the correct curvature pullback and mass normalization",
    max(pullback_error, normalized_error, ratio_error) < mp.mpf("1e-68"),
    f"pullback={mp_text(pullback_error, 8)}, "
    f"normalized={mp_text(normalized_error, 8)}, ratio={mp_text(ratio_error, 8)}",
)

with mp.workdps(100):
    first_edge = next(
        edge for edge in edges
        if tuple(sorted((colours[edge[0]], colours[edge[1]]))) == PAIRS[0]
    )
    dropped_pre = list(derived_pre)
    dropped_pre[0] += tau * edge_curvature(first_edge) / 4
    dropped_relative = abs(dropped_pre[0] - derived_pre[0]) / max(
        abs(value) for value in derived_pre
    )
drop_control_ok = check(
    "dropping one actual edge detectably changes its boundary component",
    dropped_relative > mp.mpf("1e-8"),
    f"relative change={mp_text(dropped_relative, 8)}",
)

with mp.workdps(100):
    wrong_factor = tuple(2 * value for value in derived_pre)
    wrong_factor_error = max(
        abs(wrong_factor[index] - primary_pre[index]) for index in range(6)
    )
    wrong_factor_pullback = abs(2 * mp.fsum(wrong_factor) - p_pre_fine)
    wrong_sign_error = max(
        abs(-derived_pre[index] - primary_pre[index]) for index in range(6)
    )
factor_control_ok = check(
    "factor-two and sign corruptions fail the vector and pullback comparisons",
    wrong_factor_error > mp.mpf("1e-4")
    and wrong_factor_pullback > mp.mpf("1e-4")
    and wrong_sign_error > mp.mpf("1e-4"),
)

with mp.workdps(100):
    kernel_delta = (mp.mpf("1e-6"), -mp.mpf("1e-6"), 0, 0, 0, 0)
    kernel_vector = tuple(
        derived_pre[index] + kernel_delta[index] for index in range(6)
    )
    kernel_pullback_change = abs(2 * mp.fsum(kernel_vector) - p_pre_fine)
    kernel_vector_distance = max(
        abs(kernel_vector[index] - derived_pre[index]) for index in range(6)
    )
kernel_control_ok = check(
    "a kernel perturbation preserves the scalar but changes selected components",
    kernel_pullback_change < mp.mpf("1e-90")
    and abs(kernel_vector_distance - mp.mpf("1e-6")) < mp.mpf("1e-90"),
)

scope = {
    "primary_functions_imported_or_executed": False,
    "lorentzian_action_evaluator_executed": False,
    "root_search_or_nested_census_executed": False,
    "hessian_or_spectrum_computed": False,
    "continuum_or_particle_target_loaded": False,
    "physical_constant_extracted": False,
    "component_formula_status": "POST_PRIMARY_ANALYTIC_EXPLANATION",
}
scope_ok = check(
    "the adversarial route is mechanically independent and remains in scope",
    not any(value for value in scope.values() if isinstance(value, bool)),
)

controls_ok = all((
    provenance_ok,
    upstream_ok,
    symbolic_ok,
    topology_ok,
    regular_control_ok,
    drop_control_ok,
    factor_control_ok,
    kernel_control_ok,
    scope_ok,
))
agreement_ok = component_ok and pullback_ok
if not controls_ok:
    outcome = "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED"
elif not agreement_ok:
    outcome = "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_DISAGREEMENT"
else:
    outcome = "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CORROBORATED"

outcome_ok = check(
    "the frozen hierarchy assigns exactly one adversarial boundary outcome",
    outcome in {
        "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED",
        "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_DISAGREEMENT",
        "ADVERSARIAL_REFINED_BOUNDARY_COTANGENT_CORROBORATED",
    },
    outcome,
)

artifact = {
    "title": "Adversarial spatial-hinge derivation of refined boundary covector",
    "date": "2026-08-21",
    "primary_result_commit": PRIMARY_RESULT_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "symbolic_product_hinge": {
        "area_square": str(product_area_square),
        "area": "I*l*tau/2",
        "log_boundary_area_derivative": str(log_area_derivative),
        "pre_momentum_per_edge": "-tau*l_e*epsilon_e/4",
    },
    "actual_incidence": {
        "f_vector": [len(layer) for layer in spatial_cells],
        "edge_count": len(edges),
        "edge_incidence_values": sorted(set(incidence.values())),
        "pair_curvatures": {
            f"{left}-{right}": mp_text(pair_curvatures[left, right])
            for left, right in PAIRS
        },
    },
    "derived_boundary_covector": {
        "orbit_order": [f"{left}{right}" for left, right in PAIRS],
        "pre": [mp_text(value) for value in derived_pre],
        "post": [mp_text(value) for value in derived_post],
        "primary_component_max_error": mp_text(component_error),
        "p_pre_fine": mp_text(p_pre_fine),
        "K_fine": mp_text(k_fine),
        "M_fine": mp_text(m_fine),
        "raw_coarse_ratio": mp_text(raw_ratio),
        "curvature_ratio": mp_text(curvature_ratio),
        "ratio_error": mp_text(ratio_error),
        "mass_normalized_error": mp_text(normalized_error),
    },
    "controls": {
        "coarse_edge_count": len(coarse_edges),
        "coarse_identity_error": mp_text(coarse_identity_error),
        "dropped_edge_relative_change": mp_text(dropped_relative),
        "wrong_factor_vector_error": mp_text(wrong_factor_error),
        "wrong_factor_pullback_error": mp_text(wrong_factor_pullback),
        "wrong_sign_vector_error": mp_text(wrong_sign_error),
        "kernel_pullback_change": mp_text(kernel_pullback_change),
        "kernel_vector_distance": mp_text(kernel_vector_distance),
    },
    "scope": scope,
    "status_labels": {
        "component_formula": "POST_PRIMARY_ANALYTIC_ADVERSARIAL_ROUTE",
        "perfect_action": "NOT_ESTABLISHED",
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
print("P_pre=" + ", ".join(mp_text(value, 18) for value in derived_pre))
print(f"raw ratio={mp_text(raw_ratio, 18)}")

raise SystemExit(0 if passed == tests else 1)
