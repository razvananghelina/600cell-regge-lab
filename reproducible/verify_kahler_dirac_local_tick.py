#!/usr/bin/env python3
"""Canonical signed Grover--Szegedy tick on the 600-cell cochains.

The protocol was frozen in commit ff6b2ce before this spectrum-free gate was
evaluated.  The verifier uses no particle, Planck-scale or A1 target.
"""

from collections import Counter, defaultdict, deque
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sparse


OUTPUT = Path(__file__).with_name("kahler_dirac_local_tick.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def permutation_parity(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left+1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def mod2_rank(matrix):
    """Exact rank over F2 using Python-integer row bitsets."""
    matrix = matrix.tocsr()
    pivots = {}
    rank = 0
    for row in range(matrix.shape[0]):
        bits = 0
        for column in matrix.indices[matrix.indptr[row]:matrix.indptr[row+1]]:
            bits ^= 1 << int(column)
        while bits:
            pivot = bits.bit_length()-1
            if pivot in pivots:
                bits ^= pivots[pivot]
            else:
                pivots[pivot] = bits
                rank += 1
                break
    return rank


def sparse_max_abs(matrix):
    return 0.0 if matrix.nnz == 0 else float(np.max(np.abs(matrix.data)))


def build_walk(simplex_count, incidences, simplex_degrees):
    """Build arcs, signed local reflections, flip shift and U=S C.

    Each incidence is (lower_simplex, higher_simplex, incidence_sign).
    Arc phase is the incidence sign in the upward direction and +1 downward.
    """
    arc_tail = []
    arc_head = []
    arc_phase = []
    arc_lookup = {}
    for lower, higher, sign in incidences:
        for tail, head, phase in (
            (lower, higher, sign),
            (higher, lower, 1),
        ):
            arc_lookup[(tail, head)] = len(arc_tail)
            arc_tail.append(tail)
            arc_head.append(head)
            arc_phase.append(phase)
    arc_tail = np.asarray(arc_tail, dtype=np.int32)
    arc_head = np.asarray(arc_head, dtype=np.int32)
    arc_phase = np.asarray(arc_phase, dtype=np.int8)
    arc_count = len(arc_tail)

    outgoing = [[] for _ in range(simplex_count)]
    for arc, tail in enumerate(arc_tail):
        outgoing[int(tail)].append(arc)
    degree_array = np.asarray([len(arcs) for arcs in outgoing], dtype=np.int32)
    if not np.array_equal(degree_array, np.asarray(simplex_degrees)):
        raise RuntimeError("arc degrees do not match the supplied Hasse degrees")

    reverse = np.asarray([
        arc_lookup[(int(head), int(tail))]
        for tail, head in zip(arc_tail, arc_head)
    ], dtype=np.int32)
    shift = sparse.csr_matrix(
        (np.ones(arc_count), (reverse, np.arange(arc_count))),
        shape=(arc_count, arc_count),
    )

    coin_rows = []
    coin_columns = []
    coin_data = []
    exact_coin_blocks = True
    for tail, arcs in enumerate(outgoing):
        q = len(arcs)
        phases = arc_phase[arcs].astype(np.int64)
        numerator = 2*np.outer(phases, phases)-q*np.eye(q, dtype=np.int64)
        exact_coin_blocks &= np.array_equal(
            numerator @ numerator, q*q*np.eye(q, dtype=np.int64)
        )
        for local_row, row in enumerate(arcs):
            for local_column, column in enumerate(arcs):
                value = numerator[local_row, local_column]/q
                if value:
                    coin_rows.append(row)
                    coin_columns.append(column)
                    coin_data.append(value)
    coin = sparse.csr_matrix(
        (coin_data, (coin_rows, coin_columns)),
        shape=(arc_count, arc_count),
    )
    walk = shift @ coin

    embedding = sparse.csr_matrix(
        (arc_phase/np.sqrt(degree_array[arc_tail]),
         (np.arange(arc_count), arc_tail)),
        shape=(arc_count, simplex_count),
    )
    return {
        "arc_tail": arc_tail,
        "arc_head": arc_head,
        "arc_phase": arc_phase,
        "arc_lookup": arc_lookup,
        "outgoing": outgoing,
        "reverse": reverse,
        "shift": shift,
        "coin": coin,
        "walk": walk,
        "embedding": embedding,
        "exact_coin_blocks": exact_coin_blocks,
    }


print("="*78)
print("SIGNED KAEHLER--DIRAC LOCAL UNITARY TICK")
print("="*78)

# -------------------------------------------------------------------------
# Reconstruct the same complete oriented 600-cell cochain complex.
# -------------------------------------------------------------------------
sqrt5 = np.sqrt(5.0)
phi = (1+sqrt5)/2
vertex_set = set()
for coordinate in range(4):
    for sign in (-1.0, 1.0):
        vertex = [0.0]*4
        vertex[coordinate] = sign
        vertex_set.add(tuple(vertex))
for signs in product((-0.5, 0.5), repeat=4):
    vertex_set.add(signs)
base = [phi/2, 0.5, 1/(2*phi), 0.0]
even_permutations = [
    permutation for permutation in permutations(range(4))
    if permutation_parity(permutation) == 1
]
for permutation in even_permutations:
    vertex = [base[permutation[index]] for index in range(4)]
    nonzero = [index for index, value in enumerate(vertex)
               if abs(value) > 1e-12]
    for signs in product((-1, 1), repeat=3):
        signed = vertex[:]
        for index, sign in zip(nonzero, signs):
            signed[index] *= sign
        vertex_set.add(tuple(round(value, 10) for value in signed))
vertices = np.asarray(sorted(vertex_set))
dot_products = vertices @ vertices.T
edges = [
    (left, right)
    for left in range(120)
    for right in range(left+1, 120)
    if abs(dot_products[left, right]-phi/2) < 1e-3
]
adjacency = defaultdict(set)
for left, right in edges:
    adjacency[left].add(right)
    adjacency[right].add(left)
triangles = []
for left, right in edges:
    for third in adjacency[left] & adjacency[right]:
        if right < third:
            triangles.append((left, right, third))
tetrahedra = []
for first, second, third in triangles:
    for fourth in adjacency[first] & adjacency[second] & adjacency[third]:
        if third < fourth:
            tetrahedra.append((first, second, third, fourth))
cells = [
    [(index,) for index in range(120)],
    edges,
    triangles,
    tetrahedra,
]
dims = tuple(map(len, cells))
offsets = np.cumsum((0,)+dims)
check("the oriented cochain carrier has f-vector (120,720,1200,600)",
      dims == (120, 720, 1200, 600), str(dims))

indices = [{simplex: index for index, simplex in enumerate(layer)}
           for layer in cells]
coboundaries = []
incidences = []
for degree in range(3):
    rows = []
    columns = []
    values = []
    for high_index, simplex in enumerate(cells[degree+1]):
        for omitted in range(degree+2):
            face = simplex[:omitted]+simplex[omitted+1:]
            low_index = indices[degree][face]
            sign = (-1)**omitted
            rows.append(high_index)
            columns.append(low_index)
            values.append(sign)
            incidences.append((
                int(offsets[degree]+low_index),
                int(offsets[degree+1]+high_index),
                sign,
            ))
    coboundaries.append(sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(dims[degree+1], dims[degree]),
        dtype=np.int8,
    ))
check("the two consecutive coboundary products vanish exactly over Z",
      (coboundaries[1] @ coboundaries[0]).nnz == 0
      and (coboundaries[2] @ coboundaries[1]).nnz == 0)

blocks = [[None]*4 for _ in range(4)]
for degree in range(3):
    blocks[degree+1][degree] = coboundaries[degree]
    blocks[degree][degree+1] = coboundaries[degree].T
D = sparse.bmat(blocks, format="csr", dtype=np.int8)
form_degree = np.concatenate([
    np.full(size, degree, dtype=np.int8) for degree, size in enumerate(dims)
])
grading = np.where(form_degree % 2 == 0, 1, -1).astype(np.int8)
check("D is the exact nonzero signed Hasse adjacency and is grading-odd",
      D.shape == (2640, 2640)
      and D.nnz == 2*len(incidences)
      and (D-D.T).nnz == 0
      and sparse_max_abs(
          sparse.diags(grading, dtype=np.int8) @ D
          + D @ sparse.diags(grading, dtype=np.int8)
      ) == 0)

hasse_degree = np.asarray(np.abs(D).sum(axis=1)).ravel().astype(np.int32)
degree_multiset = Counter(map(int, hasse_degree))
check("the Hasse degrees are derived and nowhere zero",
      degree_multiset == Counter({12: 120, 7: 720, 5: 1200, 4: 600}),
      str(dict(sorted(degree_multiset.items()))))

# Exact mod-2 minors attain the maximal ranks allowed by d^2=0 and the
# constant 0-cocycle.  Therefore the same ranks hold over Q.
mod2_ranks = tuple(mod2_rank(matrix) for matrix in coboundaries)
betti = (
    dims[0]-mod2_ranks[0],
    dims[1]-mod2_ranks[0]-mod2_ranks[1],
    dims[2]-mod2_ranks[1]-mod2_ranks[2],
    dims[3]-mod2_ranks[2],
)
check("exact finite-field ranks imply Betti numbers (1,0,0,1)",
      mod2_ranks == (119, 601, 599)
      and betti == (1, 0, 0, 1),
      f"ranks={mod2_ranks}; betti={betti}")

# -------------------------------------------------------------------------
# Frozen signed Grover--Szegedy construction.
# -------------------------------------------------------------------------
construction = build_walk(len(form_degree), incidences, hasse_degree)
arc_tail = construction["arc_tail"]
arc_head = construction["arc_head"]
arc_phase = construction["arc_phase"]
shift = construction["shift"]
coin = construction["coin"]
walk = construction["walk"]
embedding = construction["embedding"]
reverse = construction["reverse"]
arc_count = len(arc_tail)

check("the walk carrier contains every directed incidence exactly once",
      arc_count == D.nnz == 14880
      and len(set(zip(map(int, arc_tail), map(int, arc_head)))) == arc_count)
check("every local coin is an exact reflection, with no adjustable angle",
      construction["exact_coin_blocks"],
      "for B_x=2 eta eta^T-q_x I, B_x^2=q_x^2 I over Z")
check("the flip shift is an exact involutive permutation",
      np.array_equal(reverse[reverse], np.arange(arc_count))
      and np.all(shift.data == 1)
      and np.array_equal(np.asarray(shift.sum(axis=0)).ravel(),
                         np.ones(arc_count)))

# Exact local identities above prove unitarity.  The full sparse product is
# retained only as a floating audit of the assembled 14880-dimensional U.
unitarity_residual = sparse_max_abs(
    walk.T @ walk-sparse.eye(arc_count, format="csr")
)
check("the assembled 14880-dimensional tick is unitary [numerical audit]",
      unitarity_residual < 1e-14,
      f"max |U^T U-I|={unitarity_residual:.3e}")

walk_rows, walk_columns = walk.nonzero()
locality = all(
    D[int(arc_tail[row]), int(arc_tail[column])] != 0
    for row, column in zip(walk_rows, walk_columns)
)
check("every one-tick transition crosses exactly one Hasse incidence",
      locality,
      f"nonzero walk transitions={walk.nnz}")
check("one micro-tick exchanges even and odd form degree exactly",
      all(grading[arc_tail[row]] == -grading[arc_tail[column]]
          for row, column in zip(walk_rows, walk_columns)))

# The phase convention was frozen so the two phases on an undirected
# incidence multiply to its signed coboundary coefficient.
phase_products_are_incidence = all(
    arc_phase[construction["arc_lookup"][(lower, higher)]]
    * arc_phase[construction["arc_lookup"][(higher, lower)]] == sign
    for lower, higher, sign in incidences
)
check("arc phases factor every signed incidence exactly",
      phase_products_are_incidence)

discriminant = (embedding.T @ shift @ embedding).tocsr()
normalizer = sparse.diags(1/np.sqrt(hasse_degree))
normalized_dirac = (normalizer @ D @ normalizer).tocsr()
discriminant_residual = sparse_max_abs(discriminant-normalized_dirac)
check("the walk discriminant is the normalized signed Kahler--Dirac D",
      discriminant_residual < 1e-14,
      f"max |A^T S A-Q^(-1/2) D Q^(-1/2)|="
      f"{discriminant_residual:.3e}")

# These two full operator identities imply the standard two-dimensional
# spectral map mu^2-2 lambda mu+1=0, hence mu=e^{+-i arccos(lambda)}.
spectral_identity_one = sparse_max_abs(walk @ embedding-shift @ embedding)
spectral_identity_two = sparse_max_abs(
    walk @ shift @ embedding
    - 2*shift @ embedding @ discriminant
    + embedding
)
check("the full spectral-mapping identities hold [numerical audit]",
      max(spectral_identity_one, spectral_identity_two) < 1e-14,
      f"residuals=({spectral_identity_one:.3e},"
      f"{spectral_identity_two:.3e})")
check("the two S3 harmonic modes remain zero modes of the discriminant",
      sum(betti) == 2 and np.all(hasse_degree > 0),
      "ker(Q^-1/2 D Q^-1/2)=Q^1/2 ker(D), so nullity is exactly 2")

# Hostile carrier audit.  A connected signed normalized graph has eigenvalue
# +1 exactly when its signs are balanced.  The Hasse graph is bipartite, so
# -1 has the same obstruction.  A boundary-of-boundary square supplies an
# explicit negative signed cycle, excluding both endpoints.  Consequently
# [A,S A] has full column rank 2*2640, leaving an exact 9600-dimensional
# walk-only complement in the arc dilation.
triangle = cells[2][0]
vertex = triangle[0]
edge_left = (triangle[0], triangle[1])
edge_right = (triangle[0], triangle[2])
cycle_simplices = (
    int(offsets[0]+indices[0][(vertex,)]),
    int(offsets[1]+indices[1][edge_left]),
    int(offsets[2]+indices[2][triangle]),
    int(offsets[1]+indices[1][edge_right]),
)
negative_cycle_product = 1
for index in range(4):
    negative_cycle_product *= int(D[
        cycle_simplices[index], cycle_simplices[(index+1) % 4]
    ])
spectral_sector_dimension = 2*sum(dims)
walk_only_dimension = arc_count-spectral_sector_dimension
check("an exact negative Hasse 4-cycle excludes discriminant eigenvalues +-1",
      negative_cycle_product == -1,
      f"signed cycle product={negative_cycle_product}")
check("the local arc dilation contains an exact 9600-dimensional extra sector",
      spectral_sector_dimension == 5280 and walk_only_dimension == 9600,
      "the invariant own-operator sector span(AH,S AH) has dimension 5280")

# Finite-cone witness.  The one-step support theorem above proves the general
# statement by induction; this sparse evolution audits four successive ticks.
initial_arc = construction["arc_lookup"][(incidences[0][0], incidences[0][1])]
state = np.zeros(arc_count)
state[initial_arc] = 1
initial_tail = int(arc_tail[initial_arc])
hasse_graph = [[] for _ in range(len(form_degree))]
for lower, higher, _ in incidences:
    hasse_graph[lower].append(higher)
    hasse_graph[higher].append(lower)
distances = {initial_tail: 0}
queue = deque((initial_tail,))
while queue:
    current = queue.popleft()
    for neighbor in hasse_graph[current]:
        if neighbor not in distances:
            distances[neighbor] = distances[current]+1
            queue.append(neighbor)
check("the full Hasse incidence graph is connected",
      len(distances) == len(form_degree))
cone_ok = True
norms = [float(np.dot(state, state))]
support_sizes = [1]
for tick in range(1, 5):
    state = walk @ state
    support = np.flatnonzero(np.abs(state) > 1e-14)
    cone_ok &= all(distances[int(arc_tail[arc])] <= tick for arc in support)
    norms.append(float(np.dot(state, state)))
    support_sizes.append(len(support))
check("localized evolution obeys the exact finite cone through four ticks",
      cone_ok and max(abs(value-1) for value in norms) < 1e-13,
      f"support sizes={support_sizes}; norms={norms}")
check("the local unitary dynamics is nontrivial",
      support_sizes[1] > 1 and not np.array_equal(
          np.flatnonzero(np.abs(state) > 1e-14), np.asarray((initial_arc,))
      ))

# -------------------------------------------------------------------------
# Known-answer calibration: the Hasse graph of an oriented cycle is itself
# a cycle.  Degree-two Grover reflections are signed swaps, so each chiral
# arc translates deterministically by one Hasse edge per micro-tick.
# -------------------------------------------------------------------------
cycle_vertices = 16
cycle_incidences = []
for edge in range(cycle_vertices):
    edge_simplex = cycle_vertices+edge
    cycle_incidences.append((edge, edge_simplex, -1))
    cycle_incidences.append(((edge+1) % cycle_vertices, edge_simplex, 1))
cycle = build_walk(
    2*cycle_vertices,
    cycle_incidences,
    np.full(2*cycle_vertices, 2, dtype=np.int32),
)
cycle_walk = cycle["walk"]


def deterministic_tail_after(start_arc, ticks):
    state = np.zeros(len(cycle["arc_tail"]))
    state[start_arc] = 1
    for _ in range(ticks):
        state = cycle_walk @ state
    support = np.flatnonzero(np.abs(state) > 1e-14)
    if len(support) != 1 or abs(np.dot(state, state)-1) > 1e-14:
        return None
    return int(cycle["arc_tail"][support[0]])


left_start = cycle["arc_lookup"][(0, cycle_vertices)]
right_start = cycle["arc_lookup"][(0, 2*cycle_vertices-1)]
left_positions = [deterministic_tail_after(left_start, 2*step)
                  for step in range(5)]
right_positions = [deterministic_tail_after(right_start, 2*step)
                   for step in range(5)]
expected_left = [(-step) % cycle_vertices for step in range(5)]
expected_right = [step % cycle_vertices for step in range(5)]
check("cycle calibration gives exact ballistic left/right translations",
      left_positions == expected_left and right_positions == expected_right,
      f"left={left_positions}; right={right_positions}")
check("calibrated speed is one Hasse edge per micro-tick without fitting",
      all(deterministic_tail_after(left_start, tick) is not None
          for tick in range(9)),
      "equivalently one original cycle edge per two incidence micro-ticks")

payload = {
    "protocol_commit": "ff6b2ce",
    "target_comparison_performed": False,
    "carrier": {
        "f_vector": list(dims),
        "cochain_dimension": int(sum(dims)),
        "directed_incidence_arc_dimension": arc_count,
        "hasse_degree_multiset": {
            str(key): value for key, value in sorted(degree_multiset.items())
        },
        "betti_numbers": list(betti),
        "harmonic_modes": int(sum(betti)),
    },
    "walk": {
        "definition": "U=S(2 sum_x |s_x><s_x|-I)",
        "adjustable_coin_angles": 0,
        "exact_local_coin_reflections": construction["exact_coin_blocks"],
        "exact_shift_involution": bool(np.array_equal(
            reverse[reverse], np.arange(arc_count)
        )),
        "unitarity_numerical_residual": unitarity_residual,
        "one_hasse_edge_per_micro_tick": locality,
        "exchanges_form_parity": True,
        "discriminant": "Q^(-1/2) D Q^(-1/2)",
        "discriminant_numerical_residual": discriminant_residual,
        "spectral_phase_map": "exp(+- i arccos(lambda))",
        "own_operator_invariant_sector_dimension": spectral_sector_dimension,
        "walk_only_complement_dimension": walk_only_dimension,
        "finite_cone_support_sizes_ticks_0_to_4": support_sizes,
    },
    "cycle_calibration": {
        "vertices": cycle_vertices,
        "left_vertex_positions_after_even_micro_ticks": left_positions,
        "right_vertex_positions_after_even_micro_ticks": right_positions,
        "speed_hasse_edges_per_micro_tick": 1,
    },
    "status": {
        "mathematical": (
            "DERIVED LOCAL UNITARY LIFT: a dimensionless, strictly local "
            "unitary micro-tick of the normalized signed Kahler--Dirac "
            "incidence is constructed."
        ),
        "physical": (
            "OPEN: the construction introduces a 9600-dimensional walk-only "
            "complement and does not select its own quantization functor, a "
            "refinement limit, Lorentzian units, c, G, hbar, Planck time or "
            "Planck mass."
        ),
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the structured target-blind certificate was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED_LOCAL_UNITARY_LIFT=True")
print("STRICT_CONE=1_HASSE_EDGE_PER_MICRO_TICK")
print("PLANCK_SCALE_DERIVED=False")
print("OPEN: refinement stability and fundamental physical selection.")
raise SystemExit(0 if passed == tests else 1)
