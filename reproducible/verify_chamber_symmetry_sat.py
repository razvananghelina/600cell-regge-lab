#!/usr/bin/env python3
"""SMT decision: can a first-order orientable chamber partition retain symmetry?

The query ranges over every set partition of the 120 oriented barycentric
chambers, encoded by a restricted-growth string.  It is stronger than a
search in the 90-bit contraction family.  For each prime-order conjugacy
class of A5 it asks whether the partition equivalence relation is invariant.
"""

import contextlib
from fractions import Fraction
import io
import os
import runpy
import time

import numpy as np
import sympy as sy
import z3


HERE = os.path.dirname(os.path.abspath(__file__))
with contextlib.redirect_stdout(io.StringIO()):
    data = runpy.run_path(os.path.join(HERE, "verify_oriented_chamber_double.py"))

edges = data["chamber_edges"]
reflection = data["reflection"]
gamma = data["gamma"]
rotations = data["chamber_rotations"]
edge_pairs = data["edge_pairs"]
N = 120


def compose(p, q):
    return tuple(p[q[i]] for i in range(N))


identity = tuple(range(N))


def permutation_order(p):
    current = identity
    for order in range(1, 61):
        current = compose(p, current)
        if current == identity:
            return order
    raise AssertionError("rotation order exceeds 60")


representative = {}
for p in rotations:
    order = permutation_order(p)
    if order in (2, 3, 5) and order not in representative:
        representative[order] = p

# In A5 the two order-five classes are represented by g and g^2.
g5 = representative[5]
queries = (
    ("C2", representative[2]),
    ("C3", representative[3]),
    ("C5A", g5),
    ("C5B", compose(g5, g5)),
)

positive = [i for i in range(N) if int(gamma[i]) == 1]
negative = [i for i in range(N) if int(gamma[i]) == -1]
timeout_ms = int(os.environ.get("CHAMBER_SYMMETRY_TIMEOUT_MS", "120000"))


def solve_query(name, permutation):
    labels = [z3.Int(f"{name}_p_{i}") for i in range(N)]
    maxima = [z3.Int(f"{name}_m_{i}") for i in range(N)]
    solver = z3.Solver()
    solver.set(timeout=timeout_ms)

    # Every set partition exactly once.
    solver.add(labels[0] == 0, maxima[0] == 0)
    for i in range(1, N):
        solver.add(labels[i] >= 0)
        solver.add(labels[i] <= maxima[i-1] + 1)
        solver.add(maxima[i] == z3.If(labels[i] > maxima[i-1],
                                      labels[i], maxima[i-1]))

    # First order and metric-dimension-zero orientability.
    for x, y in edges:
        solver.add(z3.Or(labels[x] == labels[y],
                         labels[reflection[x]] == labels[reflection[y]]))
    for x in positive:
        for y in negative:
            solver.add(z3.Or(labels[x] != labels[y],
                             labels[reflection[x]] != labels[reflection[y]]))

    # The symmetry may permute block labels.  What must be invariant is the
    # equivalence relation, not the restricted-growth integer names.
    for x in range(N):
        for y in range(x+1, N):
            solver.add((labels[x] == labels[y])
                       == (labels[permutation[x]] == labels[permutation[y]]))

    start = time.monotonic()
    result = solver.check()
    elapsed = time.monotonic() - start
    if result == z3.sat:
        model = solver.model()
        values = tuple(model.eval(label).as_long() for label in labels)
        return "sat", elapsed, values, ""
    if result == z3.unsat:
        return "unsat", elapsed, None, ""
    return "unknown", elapsed, None, solver.reason_unknown()


def audit_model(labels, permutation):
    blocks = sorted(set(labels))
    relabel = {old: new for new, old in enumerate(blocks)}
    labels = tuple(relabel[x] for x in labels)
    dim = len(blocks)
    first_order = all(labels[x] == labels[y]
                      or labels[reflection[x]] == labels[reflection[y]]
                      for x, y in edges)
    orientable = all(labels[x] != labels[y]
                     or labels[reflection[x]] != labels[reflection[y]]
                     for x in positive for y in negative)
    invariant = all((labels[x] == labels[y])
                    == (labels[permutation[x]] == labels[permutation[y]])
                    for x in range(N) for y in range(x+1, N))
    cap = sy.zeros(dim)
    for x in range(N):
        cap[labels[x], labels[reflection[x]]] += int(gamma[x])
    rank = cap.rank()
    determinant = int(cap.det()) if rank == dim else 0
    pfaffian = exact_pfaffian(cap) if cap == -cap.T else None
    quotient = [set() for _ in range(dim)]
    quotient_edges = set()
    for x, y in edges:
        a, b = labels[x], labels[y]
        if a != b:
            quotient[a].add(b)
            quotient[b].add(a)
            quotient_edges.add(tuple(sorted((a, b))))
    reached = {0}
    frontier = [0]
    while frontier:
        x = frontier.pop()
        new = quotient[x] - reached
        reached.update(new)
        frontier.extend(new)
    return {
        "first_order": first_order,
        "orientable": orientable,
        "invariant": invariant,
        "dim": dim,
        "rank": rank,
        "det": determinant,
        "pf": pfaffian,
        "connected": len(reached) == dim,
        "omega": 2*len(quotient_edges),
    }


def exact_pfaffian(matrix):
    """Fraction-exact skew elimination; Pf(A)^2=det(A)."""
    n = matrix.rows
    if n % 2:
        return 0
    a = [[Fraction(int(matrix[i, j])) for j in range(n)]
         for i in range(n)]
    pf = Fraction(1)
    for k in range(0, n, 2):
        pivot_column = next((j for j in range(k+1, n)
                             if a[k][j] != 0), None)
        if pivot_column is None:
            return 0
        if pivot_column != k+1:
            for row in range(n):
                a[row][k+1], a[row][pivot_column] = (
                    a[row][pivot_column], a[row][k+1]
                )
            a[k+1], a[pivot_column] = a[pivot_column], a[k+1]
            pf = -pf
        pivot = a[k][k+1]
        pf *= pivot
        for i in range(k+2, n):
            for j in range(i+1, n):
                updated = (a[i][j]
                           + (a[k+1][i]*a[k][j]
                              - a[k][i]*a[k+1][j])/pivot)
                a[i][j] = updated
                a[j][i] = -updated
    assert pf.denominator == 1
    return int(pf)


print("=" * 78, flush=True)
print("CHAMBER PARTITION SYMMETRY SAT AUDIT", flush=True)
print("scope: arbitrary partitions; first order + orientability", flush=True)
results = {}
for name, permutation in queries:
    status, elapsed, labels, reason = solve_query(name, permutation)
    results[name] = (status, labels)
    print(f"{name}: {status} in {elapsed:.2f}s", flush=True)
    if reason:
        print(f"  reason={reason}", flush=True)
    if labels is not None:
        audit = audit_model(labels, permutation)
        print(f"  audit={audit}", flush=True)
        print("  labels=" + ",".join(map(str, labels)), flush=True)

unknown = [name for name, (status, _) in results.items() if status == "unknown"]
sat = [name for name, (status, _) in results.items() if status == "sat"]
if unknown:
    print("STATUS=INCOMPLETE", flush=True)
    print("unknown=" + ",".join(unknown), flush=True)
    raise SystemExit(2)
if sat:
    print("STATUS=SYMMETRIC_LEGAL_PARTITION_EXISTS", flush=True)
    print("sat=" + ",".join(sat), flush=True)
else:
    print("STATUS=DERIVED_NO_GO", flush=True)
    print("No first-order orientable partition has a nontrivial A5 stabilizer.",
          flush=True)

print("RESULT: 4/4 prime-order A5 conjugacy-class queries decided", flush=True)

# Stronger but explicitly scoped follow-up: exhaust the 90-bit contraction
# family under one C5 subgroup.  This is the family containing the registered
# unimodular witness; it is not every arbitrary partition queried above.
edge_index = {frozenset(edge): i for i, edge in enumerate(edges)}
edge_to_pair_choice = {}
for pair_id, pair in enumerate(edge_pairs):
    for choice, edge_id in enumerate(pair):
        edge_to_pair_choice[edge_id] = (pair_id, choice)


def choice_action(permutation):
    action = []
    for pair_id, pair in enumerate(edge_pairs):
        images = []
        for edge_id in pair:
            x, y = edges[edge_id]
            image_edge = edge_index[frozenset((permutation[x], permutation[y]))]
            images.append(edge_to_pair_choice[image_edge])
        assert images[0][0] == images[1][0]
        assert images[1][1] == 1-images[0][1]
        action.append((images[0][0], images[0][1]))
    return action


def affine_choice_orbits(action):
    seen = set()
    orbits = []
    for seed in range(len(action)):
        if seed in seen:
            continue
        current, value = seed, 0
        orbit = []
        while current not in seen:
            seen.add(current)
            orbit.append((current, value))
            nxt, flip = action[current]
            current, value = nxt, value ^ flip
        if current == seed and value:
            return None
        orbits.append(orbit)
    return orbits


def labels_from_choices(choices):
    parent = list(range(N))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[y] = x

    for pair, choice in zip(edge_pairs, choices):
        x, y = edges[pair[choice]]
        union(x, y)
    roots = [find(i) for i in range(N)]
    table = {root: i for i, root in enumerate(sorted(set(roots)))}
    return tuple(table[root] for root in roots)


c5_orbits = affine_choice_orbits(choice_action(g5))
assert c5_orbits is not None
c5_total = 1 << len(c5_orbits)
c5_orientable = c5_full_rank = c5_integral = 0
c5_best_rank = 0
c5_best_dim = 0
c5_integral_bits = None
c5_integral_labels = None
for mask in range(c5_total):
    choices = [0]*len(edge_pairs)
    for orbit_id, orbit in enumerate(c5_orbits):
        seed_value = (mask >> orbit_id) & 1
        for pair_id, offset in orbit:
            choices[pair_id] = seed_value ^ offset
    labels = labels_from_choices(choices)
    # First order is automatic in this contraction family.  Keep the
    # exhaustive loop cheap: reject non-orientable assignments before doing
    # any linear algebra, then use numerical rank only as a prefilter.  Every
    # potentially full-rank case is certified again over Z with SymPy.
    if not all(labels[x] != labels[y]
               or labels[reflection[x]] != labels[reflection[y]]
               for x in positive for y in negative):
        continue
    c5_orientable += 1
    dim = max(labels) + 1
    cap_np = np.zeros((dim, dim), dtype=np.int64)
    for x in range(N):
        cap_np[labels[x], labels[reflection[x]]] += int(gamma[x])
    rank_prefilter = int(np.linalg.matrix_rank(cap_np))
    if (rank_prefilter, dim) > (c5_best_rank, c5_best_dim):
        c5_best_rank, c5_best_dim = rank_prefilter, dim
    if rank_prefilter == dim:
        cap_exact = sy.Matrix(cap_np.tolist())
        rank_exact = cap_exact.rank()
        assert rank_exact == dim
        determinant = int(cap_exact.det())
        c5_full_rank += 1
        if abs(determinant) == 1:
            c5_integral += 1
            c5_integral_bits = "".join(map(str, choices))
            c5_integral_labels = labels

print("C5_INVARIANT_VECTOR_CONTRACTION_EXHAUSTION", flush=True)
print("  scope=partitions generated by C5-invariant contraction vectors",
      flush=True)
print(f"  affine_free_bits={len(c5_orbits)}", flush=True)
print(f"  assignments={c5_total}", flush=True)
print(f"  orientable_first_order={c5_orientable}", flush=True)
print(f"  numeric_prefilter_full_rank_then_exact={c5_full_rank}", flush=True)
print(f"  exact_unimodular_witnesses_found={c5_integral}", flush=True)
print("  aggregate_counts_status=PATTERN (numeric prefilter may have false"
      " negatives)", flush=True)
print(f"  best_rank_dim=({c5_best_rank},{c5_best_dim})", flush=True)
if c5_integral_bits:
    print(f"  integral_bits={c5_integral_bits}", flush=True)
    exact_audit = audit_model(c5_integral_labels, g5)
    assert exact_audit["first_order"]
    assert exact_audit["orientable"]
    assert exact_audit["invariant"]
    assert abs(exact_audit["det"]) == 1
    assert exact_audit["det"] == exact_audit["pf"]**2
    assert abs(exact_audit["pf"]) == 1
    stabilizer = [p for p in rotations
                  if all((c5_integral_labels[x] == c5_integral_labels[y])
                         == (c5_integral_labels[p[x]]
                             == c5_integral_labels[p[y]])
                         for x in range(N) for y in range(x+1, N))]
    print(f"  exact_witness_audit={exact_audit}", flush=True)
    print(f"  exact_witness_A5_stabilizer_order={len(stabilizer)}", flush=True)
    print("  exact_witness_labels="
          + ",".join(map(str, c5_integral_labels)), flush=True)
print("  caution=closure can make a partition invariant even when its"
      " contraction vector is not invariant", flush=True)
print("STATUS_C5_INVARIANT_VECTOR_SCOPE=EXHAUSTED", flush=True)
