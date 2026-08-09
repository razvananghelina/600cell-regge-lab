"""
exp439d_v2_optimized.py
========================
Optimized search: prune at each tree level.
McKay graph is a tree, so we process level by level.

Question: among all C1-valid assignments, how many have the
"all edges phi-power" property? Is it rare or common?
"""
import numpy as np

phi = (1 + np.sqrt(5)) / 2
phi_p = -1 / phi
a1 = 5

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
    if absN % 5 in (1, 4): return True
    return False

def is_phi_power_edge(a_i, b_i, a_j, b_j, tol=0.01):
    da, db = a_j - a_i, b_j - b_i
    dz = da + db * phi
    dzg = da + db * phi_p
    if abs(dz) < 1e-12 and abs(dzg) < 1e-12: return True
    if abs(dzg) < 1e-12 or abs(dz) < 1e-12: return True
    log_phi = np.log(abs(dz / dzg)) / np.log(phi)
    return abs(log_phi - round(log_phi)) < tol

def find_valid_ab(n, t_range=15):
    valid = []
    for t in range(-t_range, t_range + 1):
        a = -n - 6*t
        b = n + 5*t
        N = norm_zphi(a, b)
        if is_irreducible_norm(N):
            valid.append((a, b))
    return valid

# Get valid assignments per node
valid = [find_valid_ab(n) for n in n_bare]

print("Valid assignments per node (C1):")
total_c1 = 1
for i in range(9):
    print(f"  {names[i]:3s} (n={n_bare[i]:2d}): {len(valid[i])}")
    total_c1 *= len(valid[i])
print(f"Total C1 combinations: {total_c1:,}")

# Tree processing order (BFS from node 0):
# Level 0: node 0 (e)
# Level 1: node 1 (u)  [parent: 0]
# Level 2: node 2 (d)  [parent: 1]
# Level 3: node 3 (s)  [parent: 2]
# Level 4: node 4 (mu) [parent: 3]
# Level 5: node 5 (c)  [parent: 4]
# Level 6: node 6 (tau) [parent: 5], node 8 (b) [parent: 5]
# Level 7: node 7 (t)  [parent: 6]
#
# Process as chain: 0,1,2,3,4,5 then branch 6,7 and 8

# Strategy 1: DFS with pruning
# At each node, check if edge to parent is phi-power

# Tree adjacency (parent of each node)
parent = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:5}
# DFS order
dfs_order = [0, 1, 2, 3, 4, 5, 6, 7, 8]

TARGET_EDGE_SUM = -6

# Counts
count_phi_power = 0
count_c2 = 0
count_both = 0
solutions_both = []

def dfs(node_idx, assignment, edge_sum_so_far):
    """DFS through tree, pruning on phi-power at each edge."""
    global count_phi_power, count_c2, count_both

    if node_idx == 9:
        # All nodes assigned - check totals
        count_phi_power += 1
        if edge_sum_so_far == TARGET_EDGE_SUM:
            count_both += 1
            if len(solutions_both) < 30:
                solutions_both.append([x for x in assignment])
        return

    node = dfs_order[node_idx]

    for a, b in valid[node]:
        # Check edge to parent (if not root)
        if node in parent:
            p = parent[node]
            ap, bp = assignment[p]
            if not is_phi_power_edge(ap, bp, a, b):
                continue
            # Compute edge norm for C2
            da, db = a - ap, b - bp
            e_norm = norm_zphi(da, db)
            new_edge_sum = edge_sum_so_far + e_norm
        else:
            new_edge_sum = edge_sum_so_far

        assignment[node] = (a, b)
        dfs(node_idx + 1, assignment, new_edge_sum)

print("\nRunning DFS with phi-power pruning...")
assignment = [None] * 9
dfs(0, assignment, 0)

print(f"\nRESULTS:")
print(f"  C1-valid with ALL edges phi-power: {count_phi_power:,}")
print(f"  Of those, also satisfying C2:      {count_both}")
print(f"  Selectivity: {count_phi_power}/{total_c1} = "
      f"{100*count_phi_power/total_c1:.4f}% of C1-space")

if count_phi_power > 0:
    print(f"  P(C2 | phi-power) = {count_both}/{count_phi_power} = "
          f"{count_both/count_phi_power:.6f}")

# Also count C2 without phi-power constraint (separate DFS)
print("\nCounting C2 solutions without phi-power constraint...")
count_c2_only = 0
solutions_c2 = []

def dfs_c2(node_idx, assignment, edge_sum_so_far):
    global count_c2_only
    if node_idx == 9:
        if edge_sum_so_far == TARGET_EDGE_SUM:
            count_c2_only += 1
            if len(solutions_c2) < 30:
                solutions_c2.append([x for x in assignment])
        return

    node = dfs_order[node_idx]
    for a, b in valid[node]:
        if node in parent:
            p = parent[node]
            ap, bp = assignment[p]
            da, db = a - ap, b - bp
            e_norm = norm_zphi(da, db)
            new_edge_sum = edge_sum_so_far + e_norm
        else:
            new_edge_sum = edge_sum_so_far
        assignment[node] = (a, b)
        dfs_c2(node_idx + 1, assignment, new_edge_sum)

assignment2 = [None] * 9
dfs_c2(0, assignment2, 0)
print(f"  C2-valid (no phi-power): {count_c2_only}")
print(f"  P(phi-power | C2) = {count_both}/{count_c2_only}")

# ============================================================
print("\n" + "=" * 80)
print("SOLUTIONS WITH BOTH PROPERTIES")
print("=" * 80)

for idx, sol in enumerate(solutions_both):
    print(f"\nSolution {idx+1}:")
    norms = []
    for i in range(9):
        a, b = sol[i]
        N = norm_zphi(a, b)
        norms.append(N)
        z_val = a + b * phi
        print(f"  {names[i]:3s}: ({a:+3d},{b:+3d}), N={N:+4d}, z={z_val:+8.4f}")

    # Edge phi-powers
    edge_sum = 0
    for ii, jj in edges:
        ai, bi = sol[ii]
        aj, bj = sol[jj]
        da, db = aj-ai, bj-bi
        dz = da + db*phi
        dzg = da + db*phi_p
        e_norm = norm_zphi(da, db)
        edge_sum += e_norm
        if abs(dzg) > 1e-12 and abs(dz) > 1e-12:
            lp = np.log(abs(dz/dzg)) / np.log(phi)
            label = f"phi^{round(lp):+d}"
        elif abs(dz) < 1e-12 and abs(dzg) < 1e-12:
            label = "0/0"
        else:
            label = "0 or inf"
        print(f"    edge ({names[ii]:>3s},{names[jj]:>3s}): N(dz)={e_norm:+3d}, dz/dz'={label}")
    print(f"  Edge sum: {edge_sum}")

    # Check if this is the ACTUAL assignment
    actual = [(0,0),(3,-2),(1,0),(1,1),(1,1),(2,1),(1,2),(4,1),(-1,4)]
    is_actual = all(sol[i] == actual[i] for i in range(9))
    if is_actual:
        print("  *** THIS IS THE ACTUAL UNIQUE ASSIGNMENT ***")

# ============================================================
print("\n" + "=" * 80)
print("C2 SOLUTIONS (without phi-power constraint) - first few")
print("=" * 80)

for idx, sol in enumerate(solutions_c2[:5]):
    print(f"\nC2 solution {idx+1}:")
    all_pp = True
    for i in range(9):
        a, b = sol[i]
        N = norm_zphi(a, b)
        print(f"  {names[i]:3s}: ({a:+3d},{b:+3d}), N={N:+4d}")

    for ii, jj in edges:
        ai, bi = sol[ii]
        aj, bj = sol[jj]
        pp, _, label = True, 0, "?"
        da, db = aj-ai, bj-bi
        dz = da+db*phi
        dzg = da+db*phi_p
        if abs(dzg) > 1e-12 and abs(dz) > 1e-12:
            lp = np.log(abs(dz/dzg)) / np.log(phi)
            pp = abs(lp - round(lp)) < 0.01
            label = f"phi^{round(lp):+d}" if pp else f"~phi^{lp:.2f}"
        else:
            label = "degen"
        if not pp: all_pp = False
        print(f"    ({names[ii]:>3s},{names[jj]:>3s}): dz/dz'={label:>12s}  phi-power: {pp}")
    print(f"  All phi-power: {all_pp}")

# ============================================================
print("\n" + "=" * 80)
print("FINAL VERDICT")
print("=" * 80)
print(f"""
  C1-valid total:          {total_c1:>12,}
  C1 + phi-power:          {count_phi_power:>12,}  ({100*count_phi_power/total_c1:.4f}%)
  C1 + C2:                 {count_c2_only:>12,}
  C1 + phi-power + C2:     {count_both:>12,}

  Phi-power selectivity: reduces C1 space by factor {total_c1/max(count_phi_power,1):.0f}x
  C2 selectivity:        reduces C1 space by factor {total_c1/max(count_c2_only,1):.0f}x
""")

if count_both == 1:
    print("  ==> PHI-POWER + C2 selects UNIQUE assignment!")
    print("  ==> But is phi-power INDEPENDENT of C2?")
    if count_phi_power > count_c2_only:
        print("  ==> phi-power is WEAKER than C2 (more solutions)")
    elif count_phi_power < count_c2_only:
        print("  ==> phi-power is STRONGER than C2 (fewer solutions)")
    else:
        print("  ==> phi-power and C2 have same selectivity")
elif count_both == count_c2_only:
    print("  ==> ALL C2 solutions have phi-power => it's a CONSEQUENCE of C1+C2")
elif count_both == 0:
    print("  ==> NO overlap! phi-power and C2 are incompatible?!")
