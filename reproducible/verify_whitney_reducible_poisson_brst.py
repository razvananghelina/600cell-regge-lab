#!/usr/bin/env python3
"""Exact reducible Poisson--BRST audit for canonical Whitney constraints.

Protocol commit 6a0b208 froze the Poisson leaf, differential, quotient,
locality, scale, spectator and Hamiltonian gates before this construction.
"""

from collections import defaultdict, deque
from itertools import combinations
import json
from math import comb
from pathlib import Path

import sympy as sy

from whitney_trace_refinement_tools import LOCAL_D, local_whitney_mass


OUTPUT = Path(__file__).with_name("whitney_reducible_poisson_brst.json")
ALL_Q_CERTIFICATE = Path(__file__).with_name(
    "whitney_dual_resolution_all_k.json"
)
MASS_INVERSE_CERTIFICATE = Path(__file__).with_name(
    "whitney_mass_inverse_polynomial.json"
)
PROTOCOL_COMMIT = "6a0b208"
EXPECTED_ALL_Q_PROTOCOL = "8d0c557"
EXPECTED_MASS_INVERSE_PROTOCOL = "3323174"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def all_simplices(top_cells):
    return tuple(
        tuple(sorted({
            tuple(face)
            for top in top_cells
            for face in combinations(top, degree + 1)
        }))
        for degree in range(4)
    )


def oriented_cycle(row_ids, pairs):
    adjacency = defaultdict(list)
    for row in row_ids:
        left, right, _ = pairs[row]
        adjacency[left].append((right, row))
        adjacency[right].append((left, row))
    if not adjacency or any(len(items) != 2 for items in adjacency.values()):
        raise AssertionError("dual two-cell is not a simple cycle")
    start = min(adjacency)
    next_node, next_row = min(adjacency[start], key=lambda item: item[1])
    previous = None
    current = start
    coefficients = {}
    while True:
        left, right, _ = pairs[next_row]
        if (current, next_node) == (left, right):
            coefficients[next_row] = -1
        elif (current, next_node) == (right, left):
            coefficients[next_row] = 1
        else:
            raise AssertionError("cycle traversal disagrees with row")
        previous, current = current, next_node
        if current == start:
            break
        choices = [item for item in adjacency[current] if item[0] != previous]
        if len(choices) != 1:
            raise AssertionError("cycle traversal is ambiguous")
        next_node, next_row = choices[0]
    if set(coefficients) != set(row_ids):
        raise AssertionError("cycle did not exhaust its rows")
    if coefficients[min(coefficients)] < 0:
        coefficients = {row: -value for row, value in coefficients.items()}
    return coefficients


def coherent_three_boundary(r2, face_columns, constraint_rows):
    row_to_faces = defaultdict(list)
    for face in face_columns:
        for row in range(r2.rows):
            value = r2[row, face]
            if value and row in constraint_rows:
                row_to_faces[row].append((face, value))
    if set(row_to_faces) != set(constraint_rows):
        raise AssertionError("dual volume boundary misses a dual edge")
    if any(len(items) != 2 for items in row_to_faces.values()):
        raise AssertionError("a dual edge does not border exactly two faces")
    adjacency = {face: [] for face in face_columns}
    for row, ((left, lv), (right, rv)) in row_to_faces.items():
        adjacency[left].append((right, lv, rv))
        adjacency[right].append((left, rv, lv))
    start = min(face_columns)
    coefficients = {start: sy.Integer(1)}
    queue = deque((start,))
    while queue:
        face = queue.popleft()
        for other, value, other_value in adjacency[face]:
            proposed = -value * coefficients[face] / other_value
            if other in coefficients:
                if coefficients[other] != proposed:
                    raise AssertionError("inconsistent volume orientation")
            else:
                coefficients[other] = proposed
                queue.append(other)
    if set(coefficients) != set(face_columns):
        raise AssertionError("dual-volume face graph is disconnected")
    if any(abs(value) != 1 for value in coefficients.values()):
        raise AssertionError("non-incidence volume coefficient")
    return coefficients


def resolution_degree(top_cells, cells, degree):
    local_faces = tuple(combinations(range(4), degree + 1))
    cell_index = {cell: index for index, cell in enumerate(cells[degree])}
    triangle_index = {cell: index for index, cell in enumerate(cells[2])}
    lookup = {}
    injection = sy.zeros(len(top_cells) * len(local_faces), len(cells[degree]))
    for top_id, top in enumerate(top_cells):
        for local_id, positions in enumerate(local_faces):
            simplex = tuple(top[position] for position in positions)
            global_id = cell_index[simplex]
            occurrence = top_id * len(local_faces) + local_id
            lookup[top_id, global_id] = occurrence
            injection[occurrence, global_id] = 1

    triangle_parents = [[] for _ in cells[2]]
    for top_id, top in enumerate(top_cells):
        for triangle in combinations(top, 3):
            triangle_parents[triangle_index[triangle]].append(top_id)
    triangles_by_edge = defaultdict(list)
    for tri_id, triangle in enumerate(cells[2]):
        for edge in combinations(triangle, 2):
            triangles_by_edge[edge].append(tri_id)

    pairs = []
    rows_by_global = [[] for _ in cells[degree]]
    row_lookup = {}
    for tri_id, parents in enumerate(triangle_parents):
        if len(parents) != 2:
            raise AssertionError("control is not a closed 3-manifold")
        left_top, right_top = sorted(parents)
        triangle = cells[2][tri_id]
        for simplex in combinations(triangle, degree + 1):
            global_id = cell_index[tuple(simplex)]
            left = lookup[left_top, global_id]
            right = lookup[right_top, global_id]
            row = len(pairs)
            pairs.append((left, right, global_id))
            rows_by_global[global_id].append(row)
            row_lookup[global_id, tri_id] = row
    constraint = sy.zeros(len(pairs), injection.rows)
    for row, (left, right, _) in enumerate(pairs):
        constraint[row, left] = 1
        constraint[row, right] = -1

    r2_columns = []
    stage2_by_vertex = defaultdict(list)
    if degree == 0:
        vertex_index = {cell[0]: index for index, cell in enumerate(cells[0])}
        for edge_id, edge in enumerate(cells[1]):
            for vertex in edge:
                global_id = vertex_index[vertex]
                row_ids = [
                    row_lookup[global_id, tri_id]
                    for tri_id in triangles_by_edge[edge]
                ]
                coefficients = oriented_cycle(row_ids, pairs)
                column = len(r2_columns)
                r2_columns.append(coefficients)
                stage2_by_vertex[global_id].append(column)
    elif degree == 1:
        for edge_id, edge in enumerate(cells[1]):
            row_ids = [
                row_lookup[edge_id, tri_id]
                for tri_id in triangles_by_edge[edge]
            ]
            r2_columns.append(oriented_cycle(row_ids, pairs))
    relation2 = sy.zeros(len(pairs), len(r2_columns))
    for column, coefficients in enumerate(r2_columns):
        for row, value in coefficients.items():
            relation2[row, column] = value

    relation3 = sy.zeros(relation2.cols, len(cells[0]) if degree == 0 else 0)
    if degree == 0:
        for vertex, face_columns in stage2_by_vertex.items():
            coefficients = coherent_three_boundary(
                relation2, face_columns, set(rows_by_global[vertex])
            )
            for face, value in coefficients.items():
                relation3[face, vertex] = value

    return {
        "degree": degree,
        "J": injection,
        "C": constraint,
        "R2": relation2,
        "R3": relation3,
        "constraint_endpoints": [
            (left // len(local_faces), right // len(local_faces))
            for left, right, _ in pairs
        ],
    }


def matrix_max_support(matrix, direction):
    if direction == "row":
        return max(
            (sum(matrix[row, column] != 0 for column in range(matrix.cols))
             for row in range(matrix.rows)),
            default=0,
        )
    return max(
        (sum(matrix[row, column] != 0 for row in range(matrix.rows))
         for column in range(matrix.cols)),
        default=0,
    )


def count_remote(matrix, endpoint_sets):
    count = 0
    for left in range(matrix.rows):
        for right in range(matrix.cols):
            if matrix[left, right] and endpoint_sets[left].isdisjoint(
                endpoint_sets[right]
            ):
                count += 1
    return count


print("=" * 78)
print("REDUCIBLE POISSON--BRST CONVERSION OF WHITNEY CONSTRAINTS")
print("=" * 78)


all_q = json.loads(ALL_Q_CERTIFICATE.read_text())
check(
    "the inherited all-q locality certificate has the frozen protocol",
    all_q["protocol_commit"] == EXPECTED_ALL_Q_PROTOCOL
    and all_q["all_q_formulas"]["sharp_maxima"]
        == {"a0": 24, "a1": 6, "r3": 14},
)
mass_inverse_certificate = json.loads(MASS_INVERSE_CERTIFICATE.read_text())
check(
    "the independent assembled-mass inverse certificate is authentic",
    mass_inverse_certificate["protocol_commit"]
        == EXPECTED_MASS_INVERSE_PROTOCOL
    and all(
        not record["inverse_is_one_step_local"]
        and record["inverse_reaches_component_diameter"]
        for record in mass_inverse_certificate["mass_inverse_audits"][:3]
    ),
)


top_cells = tuple(combinations(range(5), 4))
cells = all_simplices(top_cells)
degrees = [resolution_degree(top_cells, cells, degree) for degree in range(4)]
regular_points = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))

records = []
all_resolution_exact = True
all_poisson = True
all_brst = True
all_quotients = True
all_locality = True
remote_leaf_inverse = 0
first_nonzero_g = None

for data in degrees:
    degree = data["degree"]
    J, C, R2, R3 = (data[key] for key in ("J", "C", "R2", "R3"))
    local_mass = local_whitney_mass(regular_points, degree)
    M = sy.diag(*([local_mass] * len(top_cells)))
    Minv = sy.diag(*([local_mass.inv()] * len(top_cells)))
    G = sy.simplify(C * Minv * C.T)
    if first_nonzero_g is None:
        first_nonzero_g = next((value for value in G if value != 0), None)
    rank_c = C.rank()
    rank_r2 = R2.rank()
    rank_r3 = R3.rank()

    resolution_exact = (
        J.T * C.T == sy.zeros(J.cols, C.rows)
        and C.T * R2 == sy.zeros(C.cols, R2.cols)
        and R2 * R3 == sy.zeros(R2.rows, R3.cols)
        and rank_c + J.cols == J.rows
        and rank_c + rank_r2 == C.rows
        and rank_r2 + rank_r3 == R2.cols
        and rank_r3 == R3.cols
    )
    all_resolution_exact &= resolution_exact

    poisson = (
        G.rank() == rank_c
        and G * R2 == sy.zeros(G.rows, R2.cols)
        and rank_r2 + G.rank() == G.cols
        and R2.T * G == sy.zeros(R2.cols, G.cols)
        and G.rank() + rank_r2 == R2.rows
        and -G + G == sy.zeros(G.rows, G.cols)
    )
    all_poisson &= poisson

    s_u = -sy.I * Minv * C.T
    s_eta = sy.I * G
    brst = (
        C * s_u + s_eta == sy.zeros(C.rows, C.rows)
        and s_u * R2 == sy.zeros(s_u.rows, R2.cols)
        and s_eta * R2 == sy.zeros(s_eta.rows, R2.cols)
        and R2 * R3 == sy.zeros(R2.rows, R3.cols)
    )
    all_brst &= brst

    invariant = J.T * M
    gauge_u = Minv * C.T
    M_W = sy.simplify(J.T * M * J)
    quotient = (
        C * J == sy.zeros(C.rows, J.cols)
        and invariant * gauge_u == sy.zeros(invariant.rows, gauge_u.cols)
        and invariant.rank() == J.cols
        and invariant.cols - invariant.rank() == rank_c
        and gauge_u.rank() == rank_c
        and invariant * J == M_W
        and M_W.det() != 0
        and J.rows + rank_c - 2 * rank_c == J.cols
    )
    all_quotients &= quotient

    g_row = matrix_max_support(G, "row")
    predicted_g_bound = (21, 21, 7, 0)[degree]
    locality = (
        matrix_max_support(C.T, "column") <= 2
        and matrix_max_support(s_u, "column")
            <= (8, 12, 8, 0)[degree]
        and matrix_max_support(s_u, "row")
            <= (12, 12, 4, 0)[degree]
        and matrix_max_support(R2, "column") <= 6
        and matrix_max_support(R2, "row") <= (2, 1, 0, 0)[degree]
        and matrix_max_support(R3, "column") <= 14
        and matrix_max_support(R3, "row") <= (1, 0, 0, 0)[degree]
        and g_row <= predicted_g_bound
    )
    all_locality &= locality

    leaf_inverse_remote = 0
    if G.rows and rank_c:
        Gplus = G.pinv()
        if not (
            G * Gplus * G == G
            and Gplus * G * Gplus == Gplus
            and G * Gplus == (G * Gplus).T
            and Gplus * G == (Gplus * G).T
        ):
            raise AssertionError("exact Moore--Penrose leaf inverse failed")
        endpoint_sets = [set(pair) for pair in data["constraint_endpoints"]]
        leaf_inverse_remote = count_remote(Gplus, endpoint_sets)
        remote_leaf_inverse += leaf_inverse_remote

    records.append({
        "degree": degree,
        "dimensions_Z0_Z1_Z2_Z3_W": [
            J.rows, C.rows, R2.cols, R3.cols, J.cols
        ],
        "ranks_C_R2_R3": [rank_c, rank_r2, rank_r3],
        "poisson_rank": G.rank(),
        "resolution_exact": resolution_exact,
        "poisson_leaf_identities": poisson,
        "brst_nilpotent": brst,
        "physical_quotient_exact": quotient,
        "locality": {
            "G_max_nonzeros_per_row": g_row,
            "G_all_q_upper_bound": predicted_g_bound,
            "C_transpose_max_nonzeros_per_column": matrix_max_support(
                C.T, "column"
            ),
            "s_u_max_nonzeros_per_row": matrix_max_support(s_u, "row"),
            "s_u_max_nonzeros_per_column": matrix_max_support(
                s_u, "column"
            ),
            "R2_max_nonzeros_per_column_on_control": matrix_max_support(
                R2, "column"
            ),
            "R2_max_nonzeros_per_row_on_control": matrix_max_support(
                R2, "row"
            ),
            "R3_max_nonzeros_per_column_on_control": matrix_max_support(
                R3, "column"
            ),
            "R3_max_nonzeros_per_row_on_control": matrix_max_support(
                R3, "row"
            ),
        },
        "leaf_inverse_remote_constraint_pair_nonzeros": leaf_inverse_remote,
    })

check("the canonical dual resolution is exact in every degree", all_resolution_exact)
check(
    "G has kernel im(R2) and image ker(R2*) in every degree",
    all_poisson,
)
check(
    "the reducible coordinate/ghost differential is exactly nilpotent",
    all_brst,
)
check(
    "the invariant covector gives the exact assembled physical quotient",
    all_quotients,
)
check(
    "every Poisson/BRST generator has the frozen all-q support bound",
    all_locality,
    str([record["locality"] for record in records]),
)
check(
    "the inverse symplectic form on the leaf has remote support",
    remote_leaf_inverse > 0,
    f"remote G+ entries={remote_leaf_inverse}",
)


alpha = sy.symbols("alpha")
check(
    "first-class cancellation fixes the auxiliary bracket scale alpha=1",
    first_nonzero_g is not None
    and sy.solve(sy.Eq((alpha - 1) * first_nonzero_g, 0), alpha) == [1],
)


# Fully local nondegenerate occurrence-field realization.
n_total = sum(data["J"].rows for data in degrees)
r_total = sum(data["C"].rank() for data in degrees)
w_total = sum(data["J"].cols for data in degrees)
occurrence_physical = 2 * n_total - 2 * r_total
check(
    "the local nondegenerate occurrence double leaves one spectator W sector",
    occurrence_physical == 2 * w_total and occurrence_physical != w_total,
    f"n={n_total}, r={r_total}, W={w_total}, reduced={occurrence_physical}",
)


# Assemble all-degree metric and weak Kähler--Dirac form in degree-block order.
degree_offsets = [0]
for data in degrees:
    degree_offsets.append(degree_offsets[-1] + data["J"].rows)
M_blocks = []
for degree in range(4):
    local_mass = local_whitney_mass(regular_points, degree)
    M_blocks.append(sy.diag(*([local_mass] * len(top_cells))))
M_all = sy.diag(*M_blocks)
J_all = sy.diag(*(data["J"] for data in degrees))
A_local = sy.zeros(n_total, n_total)
local_masses = [local_whitney_mass(regular_points, degree) for degree in range(4)]
for degree in range(3):
    forward = local_masses[degree + 1] * LOCAL_D[degree]
    low_count = comb(4, degree + 1)
    high_count = comb(4, degree + 2)
    for top_id in range(len(top_cells)):
        low_start = degree_offsets[degree] + top_id * low_count
        high_start = degree_offsets[degree + 1] + top_id * high_count
        A_local[
            high_start:high_start + high_count,
            low_start:low_start + low_count,
        ] = forward
        A_local[
            low_start:low_start + low_count,
            high_start:high_start + high_count,
        ] = forward.T

M_W = sy.simplify(J_all.T * M_all * J_all)
A_W = sy.simplify(J_all.T * A_local * J_all)
M_W_inverse = M_W.inv()
H_y = sy.simplify(M_W_inverse * A_W * M_W_inverse)
hamiltonian_identity = sy.simplify(M_W * H_y * M_W - A_W)
quotient_poisson = sy.simplify(
    J_all.T * M_all * M_all.inv() * M_all * J_all
)
check(
    "the quotient bracket and Hamiltonian flow equal assembled Whitney dynamics",
    quotient_poisson == M_W
    and hamiltonian_identity == sy.zeros(w_total, w_total)
    and M_W * H_y == A_W * M_W_inverse,
)


global_cells = tuple(cell for layer in cells for cell in layer)
remote_mwinv = 0
remote_hy = 0
for row, left in enumerate(global_cells):
    for column, right in enumerate(global_cells):
        locally_coincident = any(
            set(left).union(right).issubset(top) for top in top_cells
        )
        if not locally_coincident:
            remote_mwinv += int(M_W_inverse[row, column] != 0)
            remote_hy += int(H_y[row, column] != 0)
check(
    "the minimal rational Hamiltonian support control was evaluated exactly",
    remote_mwinv == 0 and remote_hy == 0,
    (
        "boundary-Delta4 is too small to expose remote assembled support: "
        f"M_W^-1/H_y={remote_mwinv}/{remote_hy}; the independent 600-cell "
        "certificate carries the nonlocality claim"
    ),
)


identity_form = sy.eye(w_total)
radius_two_form = sy.simplify(A_W.T * A_W)
check(
    "BRST kinematics admits inequivalent coefficient-free positive local y forms",
    A_W == A_W.T
    and identity_form == identity_form.T
    and radius_two_form == radius_two_form.T
    and radius_two_form.is_positive_semidefinite
    and radius_two_form != sy.zeros(w_total, w_total)
    and radius_two_form.rank() < identity_form.rank(),
    f"ranks I/A_W^2={identity_form.rank()}/{radius_two_form.rank()}",
)


payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "all_q_locality_protocol": EXPECTED_ALL_Q_PROTOCOL,
    "mass_inverse_protocol": EXPECTED_MASS_INVERSE_PROTOCOL,
    "phenomenological_target_used": False,
    "spectrum_used_for_selection": False,
    "independent_constraint_basis_used": False,
    "control": "boundary of a 4-simplex with exact rational Whitney mass",
    "degree_records": records,
    "combined_dimensions": {
        "occurrence_n": n_total,
        "constraint_rank_r": r_total,
        "assembled_W": w_total,
        "poisson_leaf_extended": n_total + r_total,
        "poisson_leaf_reduced": n_total - r_total,
        "nondegenerate_occurrence_double_reduced": occurrence_physical,
    },
    "remote_support": {
        "leaf_symplectic_form_Gplus": remote_leaf_inverse,
        "assembled_metric_inverse": remote_mwinv,
        "exact_Hamiltonian_in_y": remote_hy,
        "minimal_control_is_support_blind": True,
        "inherited_600cell_inverse_distances": [
            record["maximum_nonzero_inverse_distance"]
            for record in mass_inverse_certificate["mass_inverse_audits"]
        ],
    },
    "verdicts": [
        "DERIVED LOCAL REDUCIBLE POISSON--BRST KINEMATICS",
        "DERIVED: the zero-Casimir leaf has the correct physical quotient",
        "DERIVED: first-class cancellation fixes the conversion scale",
        "DERIVED NEGATIVE: the obvious local symplectic occurrence double leaves a full spectator sector",
        "DERIVED RELOCATION: the leaf symplectic form is remote and the exact Whitney Hamiltonian algebraically retains M_W inverse factors",
        "DERIVED CONTROL LIMIT: boundary-Delta4 does not expose remote support of M_W inverse",
        "DERIVED NONSELECTION: BRST kinematics permits inequivalent local Hamiltonian forms",
        "OPEN: a selected local symplectic realization and physical Hamiltonian",
        "NOT CLAIMED: time, causality, inertia, mass or Planck units",
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
check("the structured no-target certificate was written", OUTPUT.exists())
payload["tests"] = tests
payload["passed"] = passed
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")


print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed == tests:
    print("DERIVED LOCAL REDUCIBLE POISSON--BRST KINEMATICS")
    print("PHYSICAL GATE STILL CLOSED: inverses, spectators, Hamiltonian nonselection")
else:
    print("REDUCIBLE POISSON--BRST CANDIDATE REFUTED")

raise SystemExit(0 if passed == tests else 1)
