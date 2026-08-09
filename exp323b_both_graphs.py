"""
EXP-323b: McKay Mass Weights - Compare CORRECT vs WRONG graph
==============================================================

The weight_search_script uses edges that give legs (1,3,4) from branch at dim 5.
The CORRECT McKay graph (computed from rho_2 x rho_i) gives legs (1,2,5) from
branch at dim 6 = the standard affine E8 Dynkin diagram.

This experiment runs the search on BOTH graphs.

Author: Computational exploration for R.-C. Anghelina's framework
"""

import numpy as np
from itertools import combinations, permutations
from collections import deque

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N_gen = 3
N_eig = 9
rank_E8 = 8

# Irrep dimensions (same for both)
DIMS = np.array([1, 2, 2, 3, 3, 4, 4, 5, 6])
# Node index: 0=rho1(1), 1=rho2(2), 2=rho3(2'), 3=rho4(3), 4=rho5(3'),
#             5=rho6(4), 6=rho7(4'), 7=rho8(5), 8=rho9(6)

TARGET_SORTED = np.array([0, 3, 5, 11, 11, 16, 17, 19, 26])
N_TO_FERMION = {0: 'e', 3: 'u', 5: 'd', 11: 'mu/s',
                16: 'c', 17: 'tau', 19: 'b', 26: 't'}

# ==================================================
# TWO GRAPH DEFINITIONS
# ==================================================

# WRONG graph from weight_search_script (legs 1,3,4 - NOT E8!)
EDGES_WRONG = [(0, 1), (1, 3), (3, 5), (5, 7), (7, 8), (2, 4), (4, 6), (6, 7)]
# Tree: rho1-rho2-rho4-rho6-rho8-rho9, branch at rho8: rho8-rho7-rho5-rho3

# CORRECT McKay graph (legs 1,2,5 - standard affine E8)
# Computed from rho_2 x rho_i tensor product decomposition
EDGES_CORRECT = [(0, 1), (1, 3), (3, 6), (6, 7), (7, 8), (8, 4), (8, 5), (5, 2)]
# Tree: rho1-rho2-rho4-rho7-rho8-rho9, branch at rho9: rho9-rho5, rho9-rho6-rho3


def search_solutions(edges, graph_name):
    """Search for all non-negative integer weight solutions on given graph."""
    print(f"\n{'='*70}")
    print(f"SEARCHING: {graph_name}")
    print(f"{'='*70}")

    # Print graph structure
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    # Find branch node
    for node in range(9):
        nbrs = adj.get(node, [])
        deg = len(nbrs)
        role = " [BRANCH]" if deg == 3 else (" [ENDPOINT]" if deg == 1 else "")
        print(f"  Node {node} (rho_{node+1}, dim {DIMS[node]}): degree {deg}{role}")

    branch_nodes = [n for n in range(9) if len(adj.get(n, [])) == 3]
    if branch_nodes:
        br = branch_nodes[0]
        # Find leg lengths
        def leg_length(start, avoid):
            length = 0
            cur = start
            prev = avoid
            while True:
                nbrs = [n for n in adj[cur] if n != prev]
                if len(nbrs) == 0:
                    return length + 1
                if len(nbrs) == 1:
                    prev = cur
                    cur = nbrs[0]
                    length += 1
                else:
                    return length + 1  # hit another branch
            return length + 1

        legs = []
        for nbr in adj[br]:
            legs.append(leg_length(nbr, br))
        legs.sort()
        print(f"  Branch at node {br} (rho_{br+1}, dim {DIMS[br]}), legs: {tuple(legs)}")

    # BFS from node 0
    parent = {0: None}
    parent_edge = {0: None}
    queue = deque([0])
    visited = {0}
    n_edges = len(edges)

    while queue:
        v = queue.popleft()
        for e_idx, (a, b) in enumerate(edges):
            for u in [a, b]:
                other = b if u == a else a
                if other == v and u not in visited:
                    visited.add(u)
                    parent[u] = v
                    parent_edge[u] = e_idx
                    queue.append(u)

    # Path matrix
    remaining_nodes = list(range(1, 9))
    P = np.zeros((n_edges, n_edges))
    for idx, node in enumerate(remaining_nodes):
        v = node
        while parent[v] is not None:
            e = parent_edge[v]
            P[idx, e] = 1
            v = parent[v]

    rank = np.linalg.matrix_rank(P)
    if rank < n_edges:
        print(f"  WARNING: Path matrix rank {rank} < {n_edges}, not invertible!")
        return []

    P_inv = np.linalg.inv(P)

    # Search
    remaining_targets = [3, 5, 11, 11, 16, 17, 19, 26]
    non_degen = [3, 5, 16, 17, 19, 26]
    total_tested = 0
    solutions = []

    for degen_pair in combinations(range(8), 2):
        other_nodes_idx = [i for i in range(8) if i not in degen_pair]

        for perm in permutations(range(6)):
            total_tested += 1

            t = np.zeros(8)
            t[degen_pair[0]] = 11
            t[degen_pair[1]] = 11
            for k in range(6):
                t[other_nodes_idx[k]] = non_degen[perm[k]]

            w = P_inv @ t

            if np.min(w) < -1e-10:
                continue

            w_rounded = np.round(w)
            if not np.allclose(w, w_rounded, atol=1e-10):
                continue

            w_int = w_rounded.astype(int)

            assignment = {0: 0}
            for idx in range(8):
                assignment[remaining_nodes[idx]] = int(t[idx])

            solutions.append({
                'assignment': assignment,
                'weights': w_int.tolist(),
                'degen_pair': (remaining_nodes[degen_pair[0]],
                               remaining_nodes[degen_pair[1]]),
            })

    print(f"\n  Total tested: {total_tested}")
    print(f"  Solutions found: {len(solutions)}")

    # Deduplicate (mu/s swap)
    unique = []
    seen = set()
    for sol in solutions:
        key = tuple(sorted(sol['assignment'].items()))
        # Create swapped version (swap mu and s at their nodes)
        swapped = dict(sol['assignment'])
        dp = sol['degen_pair']
        swapped[dp[0]], swapped[dp[1]] = swapped[dp[1]], swapped[dp[0]]
        swap_key = tuple(sorted(swapped.items()))
        if key not in seen and swap_key not in seen:
            seen.add(key)
            unique.append(sol)
    print(f"  Unique (mod mu/s swap): {len(unique)}")

    # Display solutions
    weight_names = {
        0: "0", 1: "1", 2: "F(3)", 3: "N_gen", 5: "a1",
        6: "b1", 7: "a1+2", 8: "rank(E8)", 9: "N_eig",
        10: "degree", 15: "N_gen*a1"
    }

    for sol_idx, sol in enumerate(unique):
        assign = sol['assignment']
        w = sol['weights']

        print(f"\n  --- Solution {sol_idx + 1} ---")

        # Assignment
        for node in range(9):
            n = assign[node]
            fermion = N_TO_FERMION.get(n, '?')
            print(f"    rho_{node+1}(d={DIMS[node]}): n={n:2d} ({fermion})")

        # Edge weights
        print(f"    Edge weights:")
        for e_idx, (a, b) in enumerate(edges):
            wi = w[e_idx]
            name = weight_names.get(wi, str(wi))
            print(f"      rho_{a+1}--rho_{b+1}: w={wi} ({name})")

        # Verify distances
        ok = True
        for node in range(9):
            v = node
            d = 0
            while parent.get(v) is not None:
                e = parent_edge[v]
                d += w[e]
                v = parent[v]
            if d != assign[node]:
                print(f"    MISMATCH at rho_{node+1}: dist={d}, n={assign[node]}")
                ok = False
        if ok:
            print(f"    All distances match: OK")

        # Weight set analysis
        w_sorted = sorted(w)
        w_set = sorted(set(w))
        framework = sorted([0, 1, 2, N_gen, a1, b1, rank_E8, N_eig])
        match = w_sorted == framework
        print(f"    Weight multiset: {w_sorted}")
        print(f"    Framework set:   {framework}")
        print(f"    EXACT MATCH: {match}")

        # Galois pairs
        galois_pairs = [(1, 2), (3, 4), (5, 6)]  # (rho2,rho3), (rho4,rho5), (rho6,rho7)
        print(f"    Galois pairs:")
        for a_node, b_node in galois_pairs:
            na, nb = assign[a_node], assign[b_node]
            fa = N_TO_FERMION.get(na, '?')
            fb = N_TO_FERMION.get(nb, '?')
            print(f"      rho_{a_node+1}<->rho_{b_node+1}: {fa}(n={na}) <-> {fb}(n={nb}), |diff|={abs(na-nb)}")

        # Sum
        total = sum(w)
        print(f"    Total weight sum: {total}")

    return unique


# ==================================================
# RUN BOTH SEARCHES
# ==================================================
print("=" * 72)
print("EXP-323b: McKAY MASS WEIGHTS - CORRECT vs WRONG GRAPH")
print("=" * 72)
print(f"\nTarget exponents: {list(TARGET_SORTED)}")
print(f"Framework weight set: [0, 1, 2, 3, 5, 6, 8, 9]")

print("\n\n" + "#" * 72)
print("# GRAPH A: WRONG (from weight_search_script, legs 1,3,4)")
print("#" * 72)
sols_wrong = search_solutions(EDGES_WRONG, "WRONG GRAPH (legs 1,3,4)")

print("\n\n" + "#" * 72)
print("# GRAPH B: CORRECT McKAY (computed from tensor products, legs 1,2,5)")
print("#" * 72)
sols_correct = search_solutions(EDGES_CORRECT, "CORRECT MCKAY GRAPH (legs 1,2,5)")


# ==================================================
# COMPARISON
# ==================================================
print("\n\n" + "=" * 72)
print("COMPARISON")
print("=" * 72)

print(f"\n  WRONG graph (legs 1,3,4): {len(sols_wrong)} unique solutions")
print(f"  CORRECT graph (legs 1,2,5): {len(sols_correct)} unique solutions")

if sols_wrong:
    print(f"\n  WRONG graph weight sets:")
    for i, sol in enumerate(sols_wrong):
        w = sorted(sol['weights'])
        print(f"    Solution {i+1}: {w}")

if sols_correct:
    print(f"\n  CORRECT graph weight sets:")
    for i, sol in enumerate(sols_correct):
        w = sorted(sol['weights'])
        print(f"    Solution {i+1}: {w}")
else:
    print(f"\n  CORRECT graph: NO solutions with non-negative integer weights!")
    print(f"  This means the exact McKay distance approach needs modification")
    print(f"  (e.g., signed weights, or different metric on the graph)")

# ==================================================
# PART 2: CORRECT GRAPH - RELAXED SEARCH
# ==================================================
print("\n\n" + "=" * 72)
print("CORRECT GRAPH: RELAXED SEARCH (allow ANY non-negative integer weights)")
print("=" * 72)

edges = EDGES_CORRECT
adj = {}
for a, b in edges:
    adj.setdefault(a, []).append(b)
    adj.setdefault(b, []).append(a)

parent = {0: None}
parent_edge = {0: None}
queue = deque([0])
visited = {0}

while queue:
    v = queue.popleft()
    for e_idx, (a, b) in enumerate(edges):
        for u in [a, b]:
            other = b if u == a else a
            if other == v and u not in visited:
                visited.add(u)
                parent[u] = v
                parent_edge[u] = e_idx
                queue.append(u)

remaining_nodes = list(range(1, 9))
P = np.zeros((8, 8))
for idx, node in enumerate(remaining_nodes):
    v = node
    while parent[v] is not None:
        e = parent_edge[v]
        P[idx, e] = 1
        v = parent[v]
P_inv = np.linalg.inv(P)

non_degen = [3, 5, 16, 17, 19, 26]
all_solutions = []

for degen_pair in combinations(range(8), 2):
    other_idx = [i for i in range(8) if i not in degen_pair]
    for perm in permutations(range(6)):
        t = np.zeros(8)
        t[degen_pair[0]] = 11
        t[degen_pair[1]] = 11
        for k in range(6):
            t[other_idx[k]] = non_degen[perm[k]]

        w = P_inv @ t
        if np.min(w) < -1e-10:
            continue
        w_r = np.round(w)
        if not np.allclose(w, w_r, atol=1e-10):
            continue

        w_int = w_r.astype(int).tolist()
        assignment = {0: 0}
        for idx in range(8):
            assignment[remaining_nodes[idx]] = int(t[idx])

        all_solutions.append({
            'assignment': assignment,
            'weights': w_int,
            'weights_sorted': sorted(w_int),
            'weight_set': sorted(set(w_int))
        })

# Deduplicate
unique_relaxed = []
seen = set()
for sol in all_solutions:
    key = tuple(sol['weights_sorted'])
    if key not in seen:
        seen.add(key)
        unique_relaxed.append(sol)

print(f"\n  Total solutions with non-negative integer weights: {len(all_solutions)}")
print(f"  Unique weight multisets: {len(unique_relaxed)}")

for i, sol in enumerate(unique_relaxed):
    assign = sol['assignment']
    w = sol['weights']
    print(f"\n  Weight multiset {i+1}: {sol['weights_sorted']}")
    print(f"    Weight set (unique): {sol['weight_set']}")
    print(f"    Sum: {sum(w)}")

    # How many framework constants in the weight set?
    fw_consts = {0, 1, 2, 3, 5, 6, 8, 9}
    n_fw = sum(1 for x in w if x in fw_consts)
    print(f"    Framework constants in weights: {n_fw}/{len(w)}")

    # Assignment
    for node in range(9):
        n = assign[node]
        fermion = N_TO_FERMION.get(n, '?')
        print(f"      rho_{node+1}(d={DIMS[node]}): n={n:2d} ({fermion})")

    # Edge weights
    for e_idx, (a, b) in enumerate(edges):
        wi = w[e_idx]
        print(f"      rho_{a+1}--rho_{b+1}: w={wi}")


# ==================================================
# PART 3: CORRECT GRAPH - SIGNED WEIGHTS
# ==================================================
print("\n\n" + "=" * 72)
print("CORRECT GRAPH: SIGNED WEIGHTS (allow negative)")
print("=" * 72)

signed_solutions = []
for degen_pair in combinations(range(8), 2):
    other_idx = [i for i in range(8) if i not in degen_pair]
    for perm in permutations(range(6)):
        t = np.zeros(8)
        t[degen_pair[0]] = 11
        t[degen_pair[1]] = 11
        for k in range(6):
            t[other_idx[k]] = non_degen[perm[k]]

        w = P_inv @ t
        w_r = np.round(w)
        if not np.allclose(w, w_r, atol=1e-10):
            continue

        w_int = w_r.astype(int).tolist()
        w_abs = sorted([abs(x) for x in w_int])

        assignment = {0: 0}
        for idx in range(8):
            assignment[remaining_nodes[idx]] = int(t[idx])

        signed_solutions.append({
            'assignment': assignment,
            'weights': w_int,
            'abs_sorted': w_abs,
        })

# Find solutions where |weights| = {0,1,2,3,5,6,8,9}
target_abs = [0, 1, 2, 3, 5, 6, 8, 9]
matching_signed = [s for s in signed_solutions if s['abs_sorted'] == target_abs]

print(f"\n  Total solutions with integer weights (any sign): {len(signed_solutions)}")
print(f"  Solutions where |w| = {{0,1,2,3,5,6,8,9}}: {len(matching_signed)}")

if matching_signed:
    for i, sol in enumerate(matching_signed[:10]):
        assign = sol['assignment']
        w = sol['weights']
        print(f"\n  Signed solution {i+1}: weights = {w}")
        for node in range(9):
            n = assign[node]
            fermion = N_TO_FERMION.get(n, '?')
            print(f"    rho_{node+1}(d={DIMS[node]}): n={n:2d} ({fermion})")
        for e_idx, (a, b) in enumerate(edges):
            print(f"    rho_{a+1}--rho_{b+1}: w={w[e_idx]}")

print(f"\n{'='*72}")
print("EXP-323b COMPLETE")
print(f"{'='*72}")
