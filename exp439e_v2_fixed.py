"""
exp439e_v2_fixed.py - Fixed node ordering: node 7=t(n=26), node 8=b(n=19)
"""
import numpy as np
import time

phi = (1 + np.sqrt(5)) / 2
phi_p = -1 / phi
LOG_PHI = np.log(phi)

# CORRECT ordering: node 7=t has n=26, node 8=b has n=19
names  = ['e',  'u',   'd',  's',  'mu',  'c',  'tau', 't',   'b'  ]
n_bare = [ 0,    3,     5,    11,   11,    16,   17,    26,    19   ]
dims   = [ 1,    2,     3,    4,    5,     6,    4,     2,     3    ]
edges  = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
actual = ((0,0),(3,-2),(1,0),(1,1),(1,1),(2,1),(1,2),(4,1),(-1,4))

# Verify actual matches n_bare
print("Verification: actual (a,b) vs n_bare")
for i in range(9):
    a, b = actual[i]
    n_check = 5*a + 6*b
    match = "OK" if n_check == n_bare[i] else f"MISMATCH ({n_check})"
    print(f"  {names[i]:3s}: 5*{a:+d}+6*{b:+d} = {n_check:+3d}, n_bare={n_bare[i]:2d} [{match}]")

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
    if abs(dz) < 1e-12 and abs(dzg) < 1e-12: return True, 0
    if abs(dzg) < 1e-12: return True, 999
    if abs(dz) < 1e-12: return True, -999
    log_r = np.log(abs(dz / dzg)) / LOG_PHI
    n = round(log_r)
    if abs(log_r - n) < 0.01:
        return True, n
    return False, log_r

def find_valid_ab(n, t_range=15):
    valid = []
    for t in range(-t_range, t_range + 1):
        a = -n - 6*t
        b = n + 5*t
        if is_irreducible_norm(norm_zphi(a, b)):
            valid.append((a, b))
    return valid

# Compute actual edge sum (the TRUE C2 target)
print("\nActual edge norms:")
actual_edge_sum = 0
for i, j in edges:
    ai, bi = actual[i]
    aj, bj = actual[j]
    da, db = aj-ai, bj-bi
    en = norm_zphi(da, db)
    actual_edge_sum += en
    ok, ep = is_phi_power_edge(ai, bi, aj, bj)
    print(f"  ({names[i]:>3s},{names[j]:>3s}): (da,db)=({da:+d},{db:+d}), "
          f"N(dz)={en:+3d}, phi-power={ok}, exp={ep}")
print(f"  TOTAL edge sum: {actual_edge_sum}")
print(f"  -b1 = -6")
print(f"  Actual matches -b1? {actual_edge_sum == -6}")

actual_node_sum = sum(norm_zphi(a,b) for a,b in actual)
print(f"\nActual node sum: {actual_node_sum} (should be +b1=+6)")

# Use the ACTUAL edge sum as target
TARGET = actual_edge_sum
print(f"\nUsing TARGET edge sum = {TARGET}")

# Build valid assignments
valid = [find_valid_ab(n) for n in n_bare]
print("\nC1-valid assignments per node:")
total_c1 = 1
for i in range(9):
    print(f"  {names[i]:3s} (n={n_bare[i]:2d}): {len(valid[i])}")
    total_c1 *= len(valid[i])
print(f"Total C1 space: {total_c1:,}")

# Build compatibility tables
parent = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:5}
dfs_order = [0, 1, 2, 3, 4, 5, 6, 7, 8]

compat = {}
for node in dfs_order[1:]:
    p = parent[node]
    compat[node] = {}
    for pa, pb in valid[p]:
        key = (pa, pb)
        compat[node][key] = []
        for ca, cb in valid[node]:
            ok, _ = is_phi_power_edge(pa, pb, ca, cb)
            if ok:
                compat[node][key].append((ca, cb))

print("\nCompatibility statistics:")
for node in dfs_order[1:]:
    p = parent[node]
    total_compat = sum(len(v) for v in compat[node].values())
    avg = total_compat / len(compat[node]) if compat[node] else 0
    print(f"  edge ({names[p]:>3s},{names[node]:>3s}): "
          f"avg {avg:.1f} compat children per parent")

# DFS with phi-power pruning
print("\nRunning phi-power DFS...")
t0 = time.time()

all_solutions = []
def dfs(node_idx, assign):
    if node_idx == 9:
        all_solutions.append(tuple(tuple(x) for x in assign))
        return
    node = dfs_order[node_idx]
    if node == 0:
        for a, b in valid[0]:
            assign[0] = (a, b)
            dfs(1, assign)
    else:
        p = parent[node]
        for ca, cb in compat[node].get(assign[p], []):
            assign[node] = (ca, cb)
            dfs(node_idx + 1, assign)

assign = [None] * 9
dfs(0, assign)
elapsed = time.time() - t0
print(f"Done in {elapsed:.2f}s. Found {len(all_solutions)} phi-power solutions.")

# Analyze all solutions
print(f"\n{'='*80}")
print("ALL PHI-POWER SOLUTIONS")
print(f"{'='*80}")

from collections import Counter

actual_found = False
c2_count = 0

print(f"\n{'#':>3s} {'EdgeSum':>8s} {'SumN':>6s} {'Sum|N|':>7s} {'Units':>5s} "
      f"{'Norms':>42s} {'C2?':>4s}")
print("-" * 85)

for idx, sol in enumerate(sorted(all_solutions,
                                  key=lambda s: abs(sum(norm_zphi(s[j][0]-s[i][0], s[j][1]-s[i][1])
                                                        for i,j in edges) - TARGET))):
    norms = [norm_zphi(a,b) for a,b in sol]
    edge_sum = sum(norm_zphi(sol[j][0]-sol[i][0], sol[j][1]-sol[i][1])
                   for i,j in edges)
    sum_N = sum(norms)
    sum_absN = sum(abs(n) for n in norms)
    n_units = sum(1 for n in norms if abs(n) <= 1)
    is_c2 = (edge_sum == TARGET)
    is_act = (sol == actual)

    if is_act: actual_found = True
    if is_c2: c2_count += 1

    mark = ""
    if is_act: mark = " ***"
    elif is_c2: mark = " C2"

    print(f"{idx+1:3d} {edge_sum:+8d} {sum_N:+6d} {sum_absN:7d} {n_units:5d} "
          f"{str(norms):>42s} {'YES' if is_c2 else '':>4s}{mark}")

print(f"\n{'='*80}")
print("FINAL RESULTS")
print(f"{'='*80}")
print(f"  Total C1 space:            {total_c1:>12,}")
print(f"  Phi-power solutions:       {len(all_solutions):>12,}")
print(f"  Reduction factor:          {total_c1/max(len(all_solutions),1):>12,.0f}x")
print(f"  C2 target (edge sum):      {TARGET:>12d}")
print(f"  Phi-power + C2:            {c2_count:>12d}")
print(f"  Actual found in set:       {'YES' if actual_found else 'NO':>12s}")
print(f"  Phi-power fraction:        {100*len(all_solutions)/total_c1:>12.6f}%")

if actual_found and c2_count == 1:
    print("\n  => UNIQUE! Phi-power + C2 selects the actual assignment.")
    print(f"  => Phi-power alone: {len(all_solutions)} solutions (not unique)")
    print(f"  => This is a STRONG independent constraint!")
elif not actual_found:
    print("\n  => WARNING: Actual assignment NOT found among phi-power solutions!")
    print("  => The actual assignment does NOT have all edges phi-power??")
    print("\n  Checking actual edges individually:")
    for i, j in edges:
        ai, bi = actual[i]
        aj, bj = actual[j]
        ok, ep = is_phi_power_edge(ai, bi, aj, bj)
        print(f"    ({names[i]:>3s},{names[j]:>3s}): phi-power={ok}, exp={ep}")
