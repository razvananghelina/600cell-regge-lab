"""
exp439e_68_solutions.py
========================
Analyze ALL 68 phi-power solutions.
What do they share? What makes the actual one unique?
"""
import numpy as np
from collections import Counter

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

valid = [find_valid_ab(n) for n in n_bare]

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

# Collect ALL 68 solutions
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
print(f"Found {len(all_solutions)} phi-power solutions.\n")

# ============================================================
print("=" * 80)
print("PART 1: ALL 68 SOLUTIONS - NORM PATTERNS")
print("=" * 80)

actual = ((0,0),(3,-2),(1,0),(1,1),(1,1),(2,1),(1,2),(4,1),(-1,4))

# For each solution, compute key properties
sol_data = []
for sol in all_solutions:
    norms = tuple(norm_zphi(a, b) for a, b in sol)
    abs_norms = tuple(abs(n) for n in norms)

    edge_sum = 0
    edge_powers = []
    edge_norms = []
    for i, j in edges:
        ai, bi = sol[i]
        aj, bj = sol[j]
        da, db = aj - ai, bj - bi
        en = norm_zphi(da, db)
        edge_sum += en
        edge_norms.append(en)
        _, ep = is_phi_power_edge(ai, bi, aj, bj)
        edge_powers.append(ep)

    # Count units vs primes
    n_units = sum(1 for n in abs_norms if n <= 1)
    n_primes = sum(1 for n in abs_norms if n > 1)

    # Sum of |N|
    sum_absN = sum(abs_norms)

    is_actual = (sol == actual)

    sol_data.append({
        'sol': sol,
        'norms': norms,
        'abs_norms': abs_norms,
        'edge_sum': edge_sum,
        'edge_powers': tuple(edge_powers),
        'edge_norms': tuple(edge_norms),
        'n_units': n_units,
        'n_primes': n_primes,
        'sum_absN': sum_absN,
        'is_actual': is_actual,
    })

# Sort by edge_sum distance from -6
sol_data.sort(key=lambda d: abs(d['edge_sum'] - (-6)))

print(f"\n{'#':>3s} {'Edge sum':>9s} {'Sum|N|':>7s} {'Units':>5s} "
      f"{'Norms':>40s} {'Edge powers':>25s}")
print("-" * 100)
for idx, d in enumerate(sol_data):
    mark = " ***" if d['is_actual'] else ""
    ep_str = str(list(d['edge_powers']))
    print(f"{idx+1:3d} {d['edge_sum']:+9d} {d['sum_absN']:7d} {d['n_units']:5d} "
          f"{str(list(d['norms'])):>40s} {ep_str:>25s}{mark}")

# ============================================================
print("\n" + "=" * 80)
print("PART 2: WHAT PROPERTIES ARE SHARED?")
print("=" * 80)

# Edge sum distribution
edge_sums = [d['edge_sum'] for d in sol_data]
print(f"\nEdge sum (C2) distribution:")
for val, cnt in sorted(Counter(edge_sums).items()):
    mark = " <-- TARGET" if val == -6 else ""
    print(f"  edge_sum = {val:+5d}: {cnt:3d} solutions{mark}")

# Unit count distribution
unit_counts = [d['n_units'] for d in sol_data]
print(f"\nUnit count distribution:")
for val, cnt in sorted(Counter(unit_counts).items()):
    mark = " <-- ACTUAL" if val == 7 else ""
    print(f"  n_units = {val}: {cnt:3d} solutions{mark}")

# Sum |N| distribution
sum_absN_vals = [d['sum_absN'] for d in sol_data]
print(f"\nSum |N(z)| distribution:")
for val, cnt in sorted(Counter(sum_absN_vals).items()):
    mark = " <-- ACTUAL (=Tr(A^4)=48)" if val == 48 else ""
    print(f"  sum|N| = {val:5d}: {cnt:3d} solutions{mark}")

# Edge power patterns
ep_patterns = [d['edge_powers'] for d in sol_data]
print(f"\nEdge power pattern distribution:")
for pat, cnt in Counter(ep_patterns).most_common(15):
    mark = ""
    if any(d['is_actual'] and d['edge_powers'] == pat for d in sol_data):
        mark = " <-- ACTUAL"
    print(f"  {str(list(pat)):>30s}: {cnt:3d}{mark}")

# ============================================================
print("\n" + "=" * 80)
print("PART 3: NODE-BY-NODE ANALYSIS")
print("=" * 80)

# Which (a,b) values appear at each node?
for node in range(9):
    vals = Counter(d['sol'][node] for d in sol_data)
    actual_val = actual[node]
    print(f"\n  {names[node]:3s} (n={n_bare[node]:2d}):")
    for (a, b), cnt in vals.most_common():
        N = norm_zphi(a, b)
        mark = " <-- ACTUAL" if (a, b) == actual_val else ""
        print(f"    ({a:+3d},{b:+3d}) N={N:+5d}: appears {cnt:3d}/68{mark}")

# ============================================================
print("\n" + "=" * 80)
print("PART 4: THE ACTUAL SOLUTION vs CLOSEST COMPETITORS")
print("=" * 80)

# Show solutions closest to edge_sum = -6
print("\nSolutions ranked by |edge_sum - (-6)|:")
for idx, d in enumerate(sol_data[:15]):
    dist = abs(d['edge_sum'] + 6)
    mark = " ***ACTUAL***" if d['is_actual'] else ""
    print(f"  [{idx+1:2d}] edge_sum={d['edge_sum']:+5d} (dist={dist:4d}), "
          f"sum|N|={d['sum_absN']:4d}, units={d['n_units']}{mark}")
    if d['is_actual'] or idx < 5:
        for i in range(9):
            a, b = d['sol'][i]
            N = norm_zphi(a, b)
            z = a + b*phi
            print(f"       {names[i]:3s}: ({a:+3d},{b:+3d}) N={N:+4d} z={z:+8.4f}")

# ============================================================
print("\n" + "=" * 80)
print("PART 5: WHAT MAKES THE ACTUAL SOLUTION SPECIAL?")
print("=" * 80)

actual_data = [d for d in sol_data if d['is_actual']][0]

# Properties of the actual solution
print(f"\nActual solution properties:")
print(f"  Edge sum:     {actual_data['edge_sum']:+d} (=C2 target -b_1=-6)")
print(f"  Sum |N|:      {actual_data['sum_absN']} (=Tr(A^4)=48)")
print(f"  Units:        {actual_data['n_units']}/9")
print(f"  Edge powers:  {list(actual_data['edge_powers'])}")
print(f"  Edge norms:   {list(actual_data['edge_norms'])}")

# How many others share each property?
same_edge_sum = sum(1 for d in sol_data if d['edge_sum'] == actual_data['edge_sum'])
same_sum_absN = sum(1 for d in sol_data if d['sum_absN'] == actual_data['sum_absN'])
same_units = sum(1 for d in sol_data if d['n_units'] == actual_data['n_units'])
same_ep = sum(1 for d in sol_data if d['edge_powers'] == actual_data['edge_powers'])

print(f"\nSelectivity of each property (among 68 phi-power solutions):")
print(f"  Edge sum = -6:        {same_edge_sum}/68")
print(f"  Sum |N| = 48:         {same_sum_absN}/68")
print(f"  n_units = {actual_data['n_units']}:          {same_units}/68")
print(f"  Same edge powers:     {same_ep}/68")

# Combinations
both_es_saN = sum(1 for d in sol_data
                  if d['edge_sum'] == -6 and d['sum_absN'] == 48)
print(f"  Edge sum=-6 AND sum|N|=48: {both_es_saN}/68")

# ============================================================
print("\n" + "=" * 80)
print("PART 6: STRUCTURAL HIERARCHY OF CONSTRAINTS")
print("=" * 80)

# How much does each constraint reduce the space?
# C1: 1.6B -> 1.6B
# C1 + phi-power: 1.6B -> 68
# C1 + phi-power + max_units: 68 -> ?
# C1 + phi-power + C2: 68 -> 1

# Alternative: what single property (beyond C1+phi-power) selects unique?
print("\nSingle-property selectors (from 68 phi-power solutions):")

# Test: edge_sum alone
for target in sorted(set(d['edge_sum'] for d in sol_data)):
    cnt = sum(1 for d in sol_data if d['edge_sum'] == target)
    if cnt == 1:
        print(f"  edge_sum = {target:+d}: UNIQUE (selects 1/68)")
        if target == -6:
            print(f"    ^ This is the actual! C2 selects uniquely from phi-power set.")
        break

# Test: sum_absN alone
for target in sorted(set(d['sum_absN'] for d in sol_data)):
    cnt = sum(1 for d in sol_data if d['sum_absN'] == target)
    if cnt == 1:
        print(f"  sum|N| = {target}: UNIQUE (selects 1/68)")

# Test: norms tuple
norm_counter = Counter(d['norms'] for d in sol_data)
unique_norms = sum(1 for cnt in norm_counter.values() if cnt == 1)
print(f"  Unique norm tuples: {unique_norms}/68")
actual_norm_cnt = norm_counter[actual_data['norms']]
print(f"  Actual norm tuple count: {actual_norm_cnt}")

# ============================================================
print("\n" + "=" * 80)
print("PART 7: GALOIS COMPATIBILITY AS GAUGE FLATNESS")
print("=" * 80)

print("""
INTERPRETATION:

The phi-power constraint says: for every edge (i,j) of the McKay graph,
the ratio dz/dz' = phi^{n_ij} for some integer n_ij.

This means the "Galois phase" along each edge is QUANTIZED in units of
log(phi). In gauge theory language: the Wilson line around any edge
takes values in the DISCRETE group generated by phi.

Since the McKay graph is a tree (no loops), there are no non-trivial
holonomies to constrain. But the phi-power condition still reduces
1.6 billion -> 68. The remaining selection is by C2 (flatness), which
among these 68 Galois-quantized connections picks the unique FLAT one.

HIERARCHY:
  C1 (irreducibility):     1,629,388,800 solutions
  + Galois quantization:              68 solutions (24M x reduction)
  + C2 (flatness):                     1 solution  (68 x reduction)

The Galois quantization is doing 99.999996% of the work!
C2 just picks one from the remaining 68.
""")

# Actual edge powers and their meaning
print("Edge phi-powers of the actual solution:")
for (i, j), ep, en in zip(edges, actual_data['edge_powers'],
                           actual_data['edge_norms']):
    da, db = actual[j][0]-actual[i][0], actual[j][1]-actual[i][1]
    dz = da + db*phi
    dzg = da + db*phi_p
    sign = '+' if dz/dzg > 0 else '-'
    print(f"  ({names[i]:>3s},{names[j]:>3s}): dz/dz' = {sign}phi^{ep:+d}, "
          f"N(dz) = {en:+3d}, (da,db)=({da:+d},{db:+d})")

# Sum of edge powers
ep_sum = sum(actual_data['edge_powers'])
print(f"\n  Sum of edge phi-powers: {ep_sum}")
print(f"  (= sum of log_phi|dz/dz'| for all edges)")

# Check if edge power sum is special among 68
ep_sums = [sum(d['edge_powers']) for d in sol_data]
print(f"\n  Edge power sum distribution:")
for val, cnt in sorted(Counter(ep_sums).items()):
    mark = " <-- ACTUAL" if val == ep_sum else ""
    print(f"    sum = {val:+4d}: {cnt:3d}/68{mark}")
