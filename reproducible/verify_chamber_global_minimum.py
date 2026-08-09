#!/usr/bin/env python3
"""
Exact branch-and-bound attempt for the global minimum of dim Omega^1.

This loads the certified chamber geometry, then enumerates contraction
choices on the 90 J-pairs of chamber edges.  Orientability means exactly
one edge in each J-pair may be internal to a partition.  Union-find closure
propagates choices and rejects a branch as soon as both paired edges become
internal.

A completed branch is tested for integral Poincare duality.  The search
either exhausts the tree (global certificate) or exits honestly with an
INCOMPLETE status and a rigorous explored-node count.
"""

import os
import runpy
import time

import numpy as np
import sympy as sy


HERE = os.path.dirname(os.path.abspath(__file__))
data = runpy.run_path(os.path.join(HERE, "verify_oriented_chamber_double.py"))

edge_pairs = data["edge_pairs"]
chamber_edges = data["chamber_edges"]
reflection = data["reflection"]
gamma = data["gamma"]
rows = data["rows"]
cols = data["cols"]

N = 120
M = len(edge_pairs)
TIME_LIMIT = float(os.environ.get("CHAMBER_GLOBAL_SECONDS", "180"))
NODE_LIMIT = int(os.environ.get("CHAMBER_GLOBAL_NODES", "10000000"))
TARGET = int(os.environ.get("CHAMBER_GLOBAL_TARGET", "125"))


def find(parent, x):
    while parent[x] != x:
        x = parent[x]
    return x


def union(parent, size, a, b):
    a, b = find(parent, a), find(parent, b)
    if a == b:
        return
    if size[a] < size[b]:
        a, b = b, a
    parent[b] = a
    size[a] += size[b]


def labels_from(parent):
    roots = [find(parent, x) for x in range(N)]
    table = {root: i for i, root in enumerate(sorted(set(roots)))}
    return tuple(table[root] for root in roots)


def omega_dimension(labels):
    return len({
        (labels[x], labels[y])
        for x, y in zip(rows, cols)
        if labels[x] != labels[y]
    })


def cap_matrix(labels):
    dim = len(set(labels))
    cap = np.zeros((dim, dim), dtype=np.int64)
    for x in range(N):
        cap[labels[x], labels[reflection[x]]] += int(gamma[x])
    return cap


start = time.monotonic()
nodes = leaves = conflicts = integral_leaves = 0
best_omega = 126
best_bits = None
timed_out = False
seen_partitions = set()


def propagate(parent, size, state):
    """Return False on conflict; otherwise force choices visible in closure."""
    changed = True
    while changed:
        changed = False
        for i, pair in enumerate(edge_pairs):
            a0, b0 = chamber_edges[pair[0]]
            a1, b1 = chamber_edges[pair[1]]
            c0 = find(parent, a0) == find(parent, b0)
            c1 = find(parent, a1) == find(parent, b1)
            if c0 and c1:
                return False
            if state[i] == 0:
                if c1:
                    return False
                if not c0:
                    union(parent, size, a0, b0)
                    changed = True
            elif state[i] == 1:
                if c0:
                    return False
                if not c1:
                    union(parent, size, a1, b1)
                    changed = True
            elif c0:
                state[i] = 0
                changed = True
            elif c1:
                state[i] = 1
                changed = True
    return True


def dfs(parent, size, state):
    global nodes, leaves, conflicts, integral_leaves
    global best_omega, best_bits, timed_out

    nodes += 1
    if nodes >= NODE_LIMIT or time.monotonic() - start >= TIME_LIMIT:
        timed_out = True
        return
    if not propagate(parent, size, state):
        conflicts += 1
        return

    undecided = [i for i, value in enumerate(state) if value < 0]
    if not undecided:
        leaves += 1
        labels = labels_from(parent)
        signature = tuple(labels)
        if signature in seen_partitions:
            return
        seen_partitions.add(signature)
        dim = len(set(labels))
        if dim <= 1 or dim % 2:
            return
        cap = cap_matrix(labels)
        if int(np.linalg.matrix_rank(cap)) != dim:
            return
        if abs(int(sy.Matrix(cap).det())) != 1:
            return
        integral_leaves += 1
        omega = omega_dimension(labels)
        if omega < best_omega:
            best_omega = omega
            best_bits = "".join(str(x) for x in state)
            print(f"NEW_BEST omega={omega} dim={dim} bits={best_bits}", flush=True)
        return

    # Branch on the pair whose two possible unions join the largest current
    # components, a deterministic fail-first heuristic.
    def branch_score(i):
        scores = []
        for edge_id in edge_pairs[i]:
            a, b = chamber_edges[edge_id]
            ra, rb = find(parent, a), find(parent, b)
            scores.append(0 if ra == rb else size[ra] * size[rb])
        return max(scores)

    choice_index = max(undecided, key=branch_score)
    for choice in (0, 1):
        child_parent = parent.copy()
        child_size = size.copy()
        child_state = state.copy()
        child_state[choice_index] = choice
        a, b = chamber_edges[edge_pairs[choice_index][choice]]
        union(child_parent, child_size, a, b)
        dfs(child_parent, child_size, child_state)
        if timed_out:
            return


parent = list(range(N))
size = [1] * N
state = [-1] * M

# J exchanges the two choices globally, so fix the first choice without
# losing partition orbits.
state[0] = 0
a, b = chamber_edges[edge_pairs[0][0]]
union(parent, size, a, b)
dfs(parent, size, state)

elapsed = time.monotonic() - start
print("=" * 78)
print("GLOBAL CHAMBER MINIMIZATION")
print(f"nodes={nodes} conflicts={conflicts} leaves={leaves}")
print(f"unique_partitions={len(seen_partitions)} integral_leaves={integral_leaves}")
print(f"best_omega={best_omega} elapsed={elapsed:.2f}s")
if best_bits is not None:
    print(f"BEST_BITS={best_bits}")
if timed_out:
    print("STATUS=INCOMPLETE -- no global minimum claim is licensed")
    raise SystemExit(2)
print("STATUS=EXHAUSTED")
print(f"GLOBAL_MINIMUM={best_omega}")
