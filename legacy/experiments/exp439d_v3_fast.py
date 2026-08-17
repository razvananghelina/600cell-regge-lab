"""
exp439d_v3_fast.py - Fast phi-power edge search with aggressive pruning.

We know C2 has EXACTLY 1 solution (proven in paper).
Question: how many C1-valid assignments have ALL edges phi-power?
If this number is small, phi-power is a STRONG independent constraint.
"""
import numpy as np
import time

phi = (1 + np.sqrt(5)) / 2
phi_p = -1 / phi
LOG_PHI = np.log(phi)

names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
n_bare = [0, 3, 5, 11, 11, 16, 17, 19, 26]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

def norm_zphi(a, b):
    return a**2 + a*b - b**2

def is_irreducible_norm(N):
    absN = abs(N)
    if absN <= 1: return True
    if absN < 2: return False
    for p in range(2, int(absN**0.5) + 1):
        if absN % p == 0: return False
    if absN == 5: return True
    return absN % 5 in (1, 4)

def is_phi_power_edge(ai, bi, aj, bj):
    da, db = aj - ai, bj - bi
    dz = da + db * phi
    dzg = da + db * phi_p
    if abs(dz) < 1e-12 and abs(dzg) < 1e-12: return True
    if abs(dzg) < 1e-12 or abs(dz) < 1e-12: return True
    log_r = np.log(abs(dz / dzg)) / LOG_PHI
    return abs(log_r - round(log_r)) < 0.01

def find_valid_ab(n, t_range=15):
    valid = []
    for t in range(-t_range, t_range + 1):
        a = -n - 6*t
        b = n + 5*t
        if is_irreducible_norm(norm_zphi(a, b)):
            valid.append((a, b))
    return valid

valid = [find_valid_ab(n) for n in n_bare]

print("C1-valid assignments per node:")
total_c1 = 1
for i in range(9):
    print(f"  {names[i]:3s} (n={n_bare[i]:2d}): {len(valid[i])}")
    total_c1 *= len(valid[i])
print(f"Total C1 space: {total_c1:,}")

# Tree: parent[node] = parent_node
parent = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:5}
dfs_order = [0, 1, 2, 3, 4, 5, 6, 7, 8]

# STEP 1: For each (parent, child) edge, precompute which child
# assignments are phi-power-compatible with which parent assignments.
# This turns the DFS into a chain of lookups.

print("\nPrecomputing phi-power compatibility per edge...")
# For each edge (p, c), for each parent assignment, list compatible child assignments
compat = {}
for node in dfs_order[1:]:
    p = parent[node]
    compat[node] = {}
    for pa, pb in valid[p]:
        key = (pa, pb)
        compat[node][key] = []
        for ca, cb in valid[node]:
            if is_phi_power_edge(pa, pb, ca, cb):
                compat[node][key].append((ca, cb))

# Show compatibility stats
print("\nCompatibility statistics:")
for node in dfs_order[1:]:
    p = parent[node]
    total_compat = sum(len(v) for v in compat[node].values())
    avg = total_compat / len(compat[node]) if compat[node] else 0
    print(f"  edge ({names[p]:>3s},{names[node]:>3s}): "
          f"avg {avg:.1f} compatible children per parent assignment")

# STEP 2: DFS using precomputed compatibility
print("\nRunning DFS with precomputed phi-power pruning...")
t0 = time.time()

count_phi = 0
solutions = []

def dfs(node_idx, assign):
    global count_phi
    if node_idx == 9:
        count_phi += 1
        if len(solutions) < 50:
            solutions.append(tuple(assign))
        return

    node = dfs_order[node_idx]
    if node == 0:
        # Root: try all
        for a, b in valid[0]:
            assign[0] = (a, b)
            dfs(1, assign)
    else:
        p = parent[node]
        key = assign[p]
        for ca, cb in compat[node].get(key, []):
            assign[node] = (ca, cb)
            dfs(node_idx + 1, assign)

assign = [None] * 9
dfs(0, assign)

elapsed = time.time() - t0
print(f"Done in {elapsed:.2f}s")

print(f"\n{'='*70}")
print(f"RESULTS")
print(f"{'='*70}")
print(f"  C1-valid total:              {total_c1:>12,}")
print(f"  C1 + ALL edges phi-power:    {count_phi:>12,}")
print(f"  Reduction factor:            {total_c1/max(count_phi,1):>12,.0f}x")
print(f"  Phi-power fraction:          {100*count_phi/total_c1:>12.4f}%")

# STEP 3: Check C2 for phi-power solutions
print(f"\nChecking C2 (edge sum = -6) for phi-power solutions:")
TARGET = -6
actual = ((0,0),(3,-2),(1,0),(1,1),(1,1),(2,1),(1,2),(4,1),(-1,4))

count_c2_in_phi = 0
for sol in solutions:
    edge_sum = 0
    for i, j in edges:
        ai, bi = sol[i]
        aj, bj = sol[j]
        edge_sum += norm_zphi(aj-ai, bj-bi)
    is_c2 = (edge_sum == TARGET)
    is_actual = (sol == actual)
    if is_c2:
        count_c2_in_phi += 1
    if is_c2 or is_actual or len(solutions) <= 30:
        c2_mark = " C2=YES" if is_c2 else ""
        act_mark = " *** ACTUAL ***" if is_actual else ""
        norms = [norm_zphi(a,b) for a,b in sol]
        print(f"  norms={norms}, edge_sum={edge_sum:+d}{c2_mark}{act_mark}")

if count_phi > 50:
    # We only stored first 50, need to check all
    print(f"\n  (Checked first {len(solutions)} of {count_phi} solutions)")
    print("  Re-running to check all for C2...")

    count_c2_in_phi = 0
    def dfs_check(node_idx, assign):
        global count_c2_in_phi
        if node_idx == 9:
            edge_sum = 0
            for i, j in edges:
                ai, bi = assign[i]
                aj, bj = assign[j]
                edge_sum += norm_zphi(aj-ai, bj-bi)
            if edge_sum == TARGET:
                count_c2_in_phi += 1
                norms = [norm_zphi(a,b) for a,b in assign]
                is_act = tuple(assign) == actual
                act_mark = " *** ACTUAL ***" if is_act else ""
                print(f"  FOUND: norms={norms}, edge_sum={edge_sum}{act_mark}")
            return
        node = dfs_order[node_idx]
        if node == 0:
            for a, b in valid[0]:
                assign[0] = (a, b)
                dfs_check(1, assign)
        else:
            p = parent[node]
            key = assign[p]
            for ca, cb in compat[node].get(key, []):
                assign[node] = (ca, cb)
                dfs_check(node_idx + 1, assign)

    assign2 = [None] * 9
    dfs_check(0, assign2)
    print(f"  Total phi-power solutions with C2: {count_c2_in_phi}")

# STEP 4: Show some phi-power solutions with their edge structure
print(f"\n{'='*70}")
print("SAMPLE PHI-POWER SOLUTIONS (first 10)")
print(f"{'='*70}")

for idx, sol in enumerate(solutions[:10]):
    edge_sum = 0
    edge_info = []
    for i, j in edges:
        ai, bi = sol[i]
        aj, bj = sol[j]
        da, db = aj-ai, bj-bi
        dz = da + db*phi
        dzg = da + db*phi_p
        en = norm_zphi(da, db)
        edge_sum += en
        if abs(dzg) > 1e-12 and abs(dz) > 1e-12:
            lp = round(np.log(abs(dz/dzg)) / LOG_PHI)
            edge_info.append(f"{lp:+d}")
        else:
            edge_info.append("dg")

    norms = [norm_zphi(a,b) for a,b in sol]
    is_c2 = (edge_sum == TARGET)
    is_act = (sol == actual)

    marks = ""
    if is_c2: marks += " C2"
    if is_act: marks += " ***ACTUAL***"

    print(f"\n  [{idx+1}] Norms: {norms}")
    print(f"       (a,b): {list(sol)}")
    print(f"       Edge phi-powers: [{', '.join(edge_info)}]")
    print(f"       Edge N-sum: {edge_sum:+d} (target: {TARGET}){marks}")

print(f"\n{'='*70}")
print("FINAL VERDICT")
print(f"{'='*70}")
print(f"""
  Total C1 space:               {total_c1:>12,}
  Phi-power solutions:          {count_phi:>12,} ({100*count_phi/total_c1:.4f}%)
  C2 solutions (known):                    1
  Phi-power + C2:               {count_c2_in_phi:>12,}
""")

if count_phi < total_c1 * 0.01:
    print("  => Phi-power is a STRONG constraint (< 1% of C1 space)")
    if count_c2_in_phi == 1:
        print("  => Phi-power + C2 still gives the UNIQUE assignment")
        if count_phi > 1:
            print(f"  => But phi-power alone is NOT sufficient ({count_phi} solutions)")
            print(f"  => Phi-power is INDEPENDENT of C2 (not equivalent)")
            print(f"  => It's a NEW structural property!")
else:
    print("  => Phi-power is a WEAK constraint (consequence of structure)")
