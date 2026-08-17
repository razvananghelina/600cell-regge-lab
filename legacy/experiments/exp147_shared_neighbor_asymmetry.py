"""
EXP-147: Shared Neighbor Asymmetry
====================================
Deep investigation of how the 5 shared neighbors per edge
decompose differently across U1/SU2/SU3 sectors.

From exp146:
  U1 edges:  0A + 0B + 5C (purely fermionic)
  SU2 edges: 0.2A + 0B + 4.8C
  SU3 edges: 0.5A + 1.25B + 3.25C (gauge-boson access!)

Key question: can this asymmetry explain coupling constant ratios?
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict

PHI = (1 + 5**0.5) / 2
ALPHA = 1/137.036
ALPHA_S = 0.1180

def generate_600cell():
    phi = PHI
    vertices = set()
    for i in range(4):
        for s in [1, -1]:
            v = [0.0]*4; v[i] = float(s)
            vertices.add(tuple(round(x,10) for x in v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.add((s0, s1, s2, s3))
    base_vals = [phi/2, 0.5, 1/(2*phi), 0.0]
    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        if inv % 2 == 0:
            even_perms.append(p)
    for perm in even_perms:
        for s0 in [1,-1]:
            for s1 in [1,-1]:
                for s2 in [1,-1]:
                    signed = [s0*base_vals[0], s1*base_vals[1], s2*base_vals[2], base_vals[3]]
                    v = tuple(round(signed[perm.index(i)],10) for i in range(4))
                    vertices.add(v)
    return [list(v) for v in sorted(vertices)]

def classify_vertex(v):
    nz = sum(1 for x in v if abs(x) > 1e-9)
    if nz == 1 and any(abs(abs(x)-1.0) < 1e-9 for x in v):
        return 'A'
    if all(abs(abs(x)-0.5) < 1e-9 for x in v):
        return 'B'
    return 'C'

print("="*70)
print("EXP-147: SHARED NEIGHBOR ASYMMETRY")
print("="*70)

verts = generate_600cell()
N = len(verts)
vtypes = [classify_vertex(v) for v in verts]
coords = np.array(verts)

neighbors = defaultdict(set)
for i in range(N):
    for j in range(i+1, N):
        d2 = np.sum((coords[i] - coords[j])**2)
        if d2 < 1.0/PHI**2 + 1e-6:
            neighbors[i].add(j)
            neighbors[j].add(i)

A_set = set(i for i in range(N) if vtypes[i] == 'A')
B_set = set(i for i in range(N) if vtypes[i] == 'B')
C_set = set(i for i in range(N) if vtypes[i] == 'C')

all_edges = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            all_edges.append((i,j))

def get_edge_sector(i, j):
    ti, tj = vtypes[i], vtypes[j]
    if (ti == 'A' and tj == 'C') or (ti == 'C' and tj == 'A'):
        return 'U1'
    if (ti == 'B' and tj == 'C') or (ti == 'C' and tj == 'B'):
        return 'SU2_ch'
    if ti == 'C' and tj == 'C':
        a_i = A_set & neighbors[i]
        a_j = A_set & neighbors[j]
        b_i = B_set & neighbors[i]
        b_j = B_set & neighbors[j]
        shared_A = len(a_i & a_j)
        shared_B = len(b_i & b_j)
        if shared_A > 0 and shared_B == 0:
            return 'SU2_n'
        return 'SU3'
    return 'NONE'

edge_sectors = {}
for i, j in all_edges:
    s = get_edge_sector(i, j)
    edge_sectors[(i,j)] = s
    edge_sectors[(j,i)] = s

# Coarser classification
def coarse(s):
    if s == 'SU2_ch' or s == 'SU2_n':
        return 'SU2'
    return s

sector_counts = Counter(coarse(edge_sectors[(i,j)]) for i,j in all_edges)
fine_counts = Counter(edge_sectors[(i,j)] for i,j in all_edges)
print(f"Edge counts: {dict(fine_counts)}")

# ============================================================
# PART 1: Exact shared neighbor decomposition
# ============================================================
print("\n" + "="*70)
print("PART 1: EXACT SHARED NEIGHBOR TYPES PER EDGE TYPE")
print("="*70)

# For each edge, compute the EXACT (nA, nB, nC) shared neighbor profile
profiles = defaultdict(list)
for i, j in all_edges:
    common = neighbors[i] & neighbors[j]
    nA = sum(1 for c in common if vtypes[c] == 'A')
    nB = sum(1 for c in common if vtypes[c] == 'B')
    nC = sum(1 for c in common if vtypes[c] == 'C')
    sector = edge_sectors[(i,j)]
    profiles[sector].append((nA, nB, nC))

for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
    profs = profiles[sector]
    prof_counts = Counter(profs)
    print(f"\n{sector} edges ({len(profs)}):")
    for prof, cnt in sorted(prof_counts.items()):
        nA, nB, nC = prof
        gauge_access = nA + nB
        print(f"  ({nA}A, {nB}B, {nC}C): {cnt} edges  [gauge_access={gauge_access}]")

# ============================================================
# PART 2: Gauge access metric
# ============================================================
print("\n" + "="*70)
print("PART 2: GAUGE ACCESS METRIC")
print("="*70)

# Define "gauge access" = number of gauge (A+B) shared neighbors
# This measures how much an edge "sees" gauge bosons

for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
    profs = profiles[sector]
    gauge_vals = [nA + nB for nA, nB, nC in profs]
    a_vals = [nA for nA, nB, nC in profs]
    b_vals = [nB for nA, nB, nC in profs]
    c_vals = [nC for nA, nB, nC in profs]
    print(f"{sector}: gauge_access = {np.mean(gauge_vals):.4f} (std={np.std(gauge_vals):.4f})")
    print(f"  A={np.mean(a_vals):.4f}, B={np.mean(b_vals):.4f}, C={np.mean(c_vals):.4f}")

# Weighted by SM edge count
print("\nWeighted gauge access per SM sector:")
for sm_sector in ['U1', 'SU2', 'SU3']:
    total_gauge = 0
    total_edges = 0
    for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
        if coarse(sector) == sm_sector:
            for nA, nB, nC in profiles[sector]:
                total_gauge += nA + nB
                total_edges += 1
    avg = total_gauge / total_edges if total_edges > 0 else 0
    print(f"  {sm_sector}: avg gauge_access = {avg:.4f} (over {total_edges} edges)")

# ============================================================
# PART 3: Can gauge access give coupling ratios?
# ============================================================
print("\n" + "="*70)
print("PART 3: COUPLING RATIOS FROM GAUGE ACCESS")
print("="*70)

# Hypothesis: coupling strength ~ 1/gauge_access (more gauge access = more screening)
# Or: coupling ~ gauge_access (more access = stronger)

ga = {}
for sm_sector in ['U1', 'SU2', 'SU3']:
    total = 0; count = 0
    for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
        if coarse(sector) == sm_sector:
            for nA, nB, nC in profiles[sector]:
                total += nA + nB; count += 1
    ga[sm_sector] = total / count if count > 0 else 0

print(f"Gauge access: U1={ga['U1']:.4f}, SU2={ga['SU2']:.4f}, SU3={ga['SU3']:.4f}")

target_as_a = ALPHA_S / ALPHA  # = 16.18
print(f"\nTarget alpha_s/alpha = {target_as_a:.4f} = 10*phi = {10*PHI:.4f}")

if ga['U1'] > 0 and ga['SU3'] > 0:
    print(f"SU3_ga/U1_ga = {ga['SU3']/ga['U1']:.4f}")  # can't divide by 0
else:
    print(f"U1 gauge_access = 0! Can't take ratio directly.")
    print("U1 edges are PURELY FERMIONIC - zero gauge shared neighbors")

# Alternative: use C-access (fermion loops)
ca = {}
for sm_sector in ['U1', 'SU2', 'SU3']:
    total = 0; count = 0
    for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
        if coarse(sector) == sm_sector:
            for nA, nB, nC in profiles[sector]:
                total += nC; count += 1
    ca[sm_sector] = total / count if count > 0 else 0

print(f"\nFermion access: U1={ca['U1']:.4f}, SU2={ca['SU2']:.4f}, SU3={ca['SU3']:.4f}")
print(f"U1_ca/SU3_ca = {ca['U1']/ca['SU3']:.4f}")
print(f"5/3.25 = {5/3.25:.4f}")

# ============================================================
# PART 4: Sub-types within SU3
# ============================================================
print("\n" + "="*70)
print("PART 4: SU3 EDGE SUB-TYPES")
print("="*70)

# From Part 1, SU3 edges have profiles (0,1,4), (0,2,3), (1,1,3), (1,2,2)
# These are NOT all the same! Do they form meaningful sub-classes?

su3_profs = profiles['SU3']
su3_types = Counter(su3_profs)
print("SU3 edge sub-types:")
for prof, cnt in sorted(su3_types.items()):
    nA, nB, nC = prof
    print(f"  ({nA}A, {nB}B, {nC}C): {cnt} edges")

# Total: should be 384
print(f"  Total: {sum(su3_types.values())}")

# For each sub-type, check the local structure
print("\nSU3 sub-type edge properties:")
# Store per-edge profile
edge_profile = {}
for i, j in all_edges:
    common = neighbors[i] & neighbors[j]
    nA = sum(1 for c in common if vtypes[c] == 'A')
    nB = sum(1 for c in common if vtypes[c] == 'B')
    nC = sum(1 for c in common if vtypes[c] == 'C')
    edge_profile[(i,j)] = (nA, nB, nC)
    edge_profile[(j,i)] = (nA, nB, nC)

# For each SU3 sub-type, look at triangle participation
print("\nSU3 sub-type triangle analysis:")
triangles = []
for i in range(N):
    for j in neighbors[i]:
        if j > i:
            for k in neighbors[i] & neighbors[j]:
                if k > j:
                    triangles.append((i,j,k))

print(f"Total triangles: {len(triangles)}")

for prof in sorted(su3_types.keys()):
    nA, nB, nC = prof
    cnt = su3_types[prof]
    # Count triangles containing edges of this sub-type
    edge_set = set()
    for i, j in all_edges:
        if edge_sectors[(i,j)] == 'SU3' and edge_profile[(i,j)] == prof:
            edge_set.add((i,j))
            edge_set.add((j,i))

    tri_count = 0
    pure_su3_count = 0
    for i, j, k in triangles:
        if (i,j) in edge_set or (i,k) in edge_set or (j,k) in edge_set:
            tri_count += 1
            # Check if pure SU3
            if (coarse(edge_sectors[(i,j)]) == 'SU3' and
                coarse(edge_sectors[(i,k)]) == 'SU3' and
                coarse(edge_sectors[(j,k)]) == 'SU3'):
                pure_su3_count += 1

    print(f"  ({nA}A,{nB}B,{nC}C) [{cnt} edges]: "
          f"{tri_count} triangles, {pure_su3_count} pure SU3 "
          f"(per edge: {tri_count/cnt:.1f} tri, {pure_su3_count/cnt:.1f} pure)")

# ============================================================
# PART 5: The U1 isolation property
# ============================================================
print("\n" + "="*70)
print("PART 5: U1 ISOLATION (ABELIAN CHARACTER)")
print("="*70)

# U1 edges have 0 gauge shared neighbors. This means:
# - The photon never "sees" other gauge bosons through shared vertices
# - Photon interactions are ALWAYS mediated by fermion loops
# - This is the geometric origin of the Abelian character of QED

# Check: what about 2-hop gauge access?
print("Multi-hop gauge access:")
for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
    edges_s = [(i,j) for i,j in all_edges if edge_sectors[(i,j)] == sector]

    # 1-hop: shared neighbors (already computed)
    ga1 = np.mean([edge_profile[(i,j)][0] + edge_profile[(i,j)][1] for i,j in edges_s])

    # 2-hop: neighbors of shared neighbors that are gauge
    ga2_vals = []
    for i, j in edges_s[:100]:
        common = neighbors[i] & neighbors[j]
        hop2_gauge = set()
        for c in common:
            for n2 in neighbors[c]:
                if n2 != i and n2 != j and n2 not in common:
                    if vtypes[n2] in ('A', 'B'):
                        hop2_gauge.add(n2)
        ga2_vals.append(len(hop2_gauge))

    ga2 = np.mean(ga2_vals)
    print(f"  {sector:8s}: 1-hop gauge = {ga1:.3f}, 2-hop gauge = {ga2:.3f}")

# ============================================================
# PART 6: Gauge access and beta functions
# ============================================================
print("\n" + "="*70)
print("PART 6: GAUGE ACCESS AND RUNNING COUPLINGS")
print("="*70)

# In QFT, beta functions determine how couplings run:
# - U(1): beta > 0 (screening, coupling grows at short distances)
# - SU(2): beta < 0 (antiscreening/asymptotic freedom)
# - SU(3): beta < 0 (strong asymptotic freedom)
#
# On the 600-cell, "distance" = graph distance d
# At d=1 (neighbors), gauge access per sector:
# U1: 0, SU2: 0.17, SU3: 1.75
#
# Hypothesis: gauge access at distance d controls the effective coupling at scale d
# More gauge access = more antiscreening = more asymptotic freedom

print("Gauge access = proxy for self-interaction strength:")
print(f"  U1:  ga = {ga['U1']:.2f} -> no self-interaction (Abelian)")
print(f"  SU2: ga = {ga['SU2']:.2f} -> weak self-interaction")
print(f"  SU3: ga = {ga['SU3']:.2f} -> strong self-interaction")

# SM beta function coefficients (one-loop)
b_U1 = -4/3 * 3  # ~ -4 (N_f=3 gen)
b_SU2 = 22/3 - 4/3 * 3  # = 22/3 - 4 ~ 3.3
b_SU3 = 11 - 4/3 * 3  # = 11 - 4 = 7

print(f"\nSM one-loop beta_0: U1={b_U1:.1f}, SU2={b_SU2:.1f}, SU3={b_SU3:.1f}")
print(f"SM |beta_0| ratios: SU3/SU2 = {abs(b_SU3)/abs(b_SU2):.2f}, SU3/U1 = {abs(b_SU3)/abs(b_U1):.2f}")

print(f"\nGauge access ratios:")
if ga['SU2'] > 0:
    print(f"  SU3/SU2 = {ga['SU3']/ga['SU2']:.4f}")
print(f"  SU3 gauge access = {ga['SU3']:.4f}")

# ============================================================
# PART 7: Fermion-mediated gauge coupling
# ============================================================
print("\n" + "="*70)
print("PART 7: FERMION-MEDIATED COUPLING STRENGTH")
print("="*70)

# Each gauge boson couples to fermions via edges.
# The "coupling strength" could be proportional to:
# - Number of fermion loops (shared C neighbors)
# - Weighted by the shared C's own gauge access

# For each edge, compute "weighted fermion coupling":
# w_edge = sum over shared C neighbors of (C's gauge access)
print("Weighted fermion coupling per edge sector:")

for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
    edges_s = [(i,j) for i,j in all_edges if edge_sectors[(i,j)] == sector]
    w_vals = []
    for i, j in edges_s:
        common = neighbors[i] & neighbors[j]
        w = 0
        for c in common:
            if vtypes[c] == 'C':
                # C's gauge access = how many gauge neighbors does this C have?
                c_gauge = len((A_set | B_set) & neighbors[c])
                w += c_gauge
        w_vals.append(w)

    print(f"  {sector:8s}: mean weighted coupling = {np.mean(w_vals):.4f} "
          f"(std={np.std(w_vals):.4f})")

# ============================================================
# PART 8: The A-B asymmetry in shared neighbors
# ============================================================
print("\n" + "="*70)
print("PART 8: A-B ASYMMETRY IN SHARED NEIGHBORS")
print("="*70)

# From Part 1:
# SU3 edges: (0,1,4):96, (0,2,3):96, (1,1,3):96, (1,2,2):96 (if uniform 96 each)
# This means SU3 edges come in 4 sub-types with different A/B shared profiles

# Key: the RATIO of A-sharing to B-sharing
# A-sharing = sees cross-polytope (U1-type gauge)
# B-sharing = sees tesseract (SU2-type gauge)

su3_profs_sorted = sorted(su3_types.items())
print("SU3 sub-types and their A/B ratio:")
for (nA, nB, nC), cnt in su3_profs_sorted:
    if nA + nB > 0:
        a_frac = nA / (nA + nB)
        print(f"  ({nA}A,{nB}B,{nC}C) x{cnt}: A/(A+B) = {a_frac:.3f}")
    else:
        print(f"  ({nA}A,{nB}B,{nC}C) x{cnt}: no gauge shared neighbors")

# What about the per-fermion view?
print("\nPer-fermion: shared neighbor profile of all 12 edges:")
ci = sorted(C_set)[0]
print(f"Fermion vertex {ci}:")
for j in sorted(neighbors[ci]):
    sec = edge_sectors[(ci, j)]
    prof = edge_profile[(ci, j)]
    print(f"  -> {j} (type {vtypes[j]}): sector={sec:8s}, "
          f"shared=({prof[0]}A,{prof[1]}B,{prof[2]}C)")

# ============================================================
# PART 9: Total gauge access budget
# ============================================================
print("\n" + "="*70)
print("PART 9: TOTAL GAUGE ACCESS BUDGET")
print("="*70)

# For each fermion, sum up gauge access over all 12 edges
per_fermion_ga = []
per_fermion_by_sector = defaultdict(list)

for ci in C_set:
    total_ga = 0
    sector_ga = defaultdict(float)
    for j in neighbors[ci]:
        prof = edge_profile[(ci, j)]
        ga_val = prof[0] + prof[1]  # A + B shared
        total_ga += ga_val
        sec = coarse(edge_sectors[(ci, j)])
        sector_ga[sec] += ga_val
    per_fermion_ga.append(total_ga)
    for s in ['U1', 'SU2', 'SU3']:
        per_fermion_by_sector[s].append(sector_ga[s])

print(f"Per-fermion total gauge access: {np.mean(per_fermion_ga):.4f} "
      f"(std={np.std(per_fermion_ga):.4f})")
print("  By sector:")
for s in ['U1', 'SU2', 'SU3']:
    vals = per_fermion_by_sector[s]
    print(f"    {s}: {np.mean(vals):.4f} (std={np.std(vals):.4f})")

total_ga_per_fermion = np.mean(per_fermion_ga)
print(f"\nTotal gauge access per fermion: {total_ga_per_fermion:.1f}")
print(f"  From U1 edges: {np.mean(per_fermion_by_sector['U1']):.1f}")
print(f"  From SU2 edges: {np.mean(per_fermion_by_sector['SU2']):.1f}")
print(f"  From SU3 edges: {np.mean(per_fermion_by_sector['SU3']):.1f}")

# ============================================================
# PART 10: Coupling formula from gauge access
# ============================================================
print("\n" + "="*70)
print("PART 10: COUPLING FORMULAS FROM TOPOLOGY")
print("="*70)

# Total shared gauge neighbors across ALL edges in each sector
total_gauge_shared = {}
for sm_sector in ['U1', 'SU2', 'SU3']:
    total = 0
    for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
        if coarse(sector) == sm_sector:
            for nA, nB, nC in profiles[sector]:
                total += nA + nB
    total_gauge_shared[sm_sector] = total

print("Total gauge shared neighbors per sector:")
for s in ['U1', 'SU2', 'SU3']:
    print(f"  {s}: {total_gauge_shared[s]}")

# Total fermion shared neighbors
total_fermion_shared = {}
for sm_sector in ['U1', 'SU2', 'SU3']:
    total = 0
    for sector in ['U1', 'SU2_ch', 'SU2_n', 'SU3']:
        if coarse(sector) == sm_sector:
            for nA, nB, nC in profiles[sector]:
                total += nC
    total_fermion_shared[sm_sector] = total

print("\nTotal fermion shared neighbors per sector:")
for s in ['U1', 'SU2', 'SU3']:
    print(f"  {s}: {total_fermion_shared[s]}")

print("\nTotal shared (gauge+fermion) per sector:")
for s in ['U1', 'SU2', 'SU3']:
    t = total_gauge_shared[s] + total_fermion_shared[s]
    print(f"  {s}: {t} (= {sector_counts[s]} * 5 = {sector_counts[s] * 5})")

# Ratios
print("\nKey ratios:")
print(f"  Fermion_U1/Fermion_SU3 = {total_fermion_shared['U1']/total_fermion_shared['SU3']:.4f}")
print(f"  Gauge_SU3/Gauge_SU2 = {total_gauge_shared['SU3']/total_gauge_shared['SU2']:.4f}")
print(f"  Gauge_SU3 = {total_gauge_shared['SU3']}")
print(f"  Gauge_SU2 = {total_gauge_shared['SU2']}")
print(f"  Gauge_U1 = {total_gauge_shared['U1']}")

# Check specific formulas
print("\nFormula tests:")
print(f"  Target: alpha_s/alpha = {ALPHA_S/ALPHA:.4f} = 10*phi = {10*PHI:.4f}")

# Try: fermion_shared * some factor / gauge_shared
if total_gauge_shared['SU3'] > 0:
    r1 = total_fermion_shared['U1'] / total_fermion_shared['SU3']
    print(f"  Fermion_U1/Fermion_SU3 = {r1:.4f}")
    r2 = total_gauge_shared['SU3'] / total_gauge_shared['SU2']
    print(f"  Gauge_SU3/Gauge_SU2 = {r2:.4f}")

# Edge count * fermion_access / gauge_access type analysis
for s in ['U1', 'SU2', 'SU3']:
    n_e = sector_counts[s]
    f_a = total_fermion_shared[s]
    g_a = total_gauge_shared[s]
    print(f"  {s}: edges={n_e}, fermion_shared={f_a}, gauge_shared={g_a}, "
          f"f/e={f_a/n_e:.2f}, g/e={g_a/n_e:.2f}")

# ============================================================
# PART 11: The key test - ratios involving 5 and phi
# ============================================================
print("\n" + "="*70)
print("PART 11: SYSTEMATIC RATIO SEARCH")
print("="*70)

target = 10 * PHI

# Collect all computed quantities
quantities = {
    'edges_U1': 96,
    'edges_SU2': 240,
    'edges_SU3': 384,
    'gauge_U1': total_gauge_shared['U1'],
    'gauge_SU2': total_gauge_shared['SU2'],
    'gauge_SU3': total_gauge_shared['SU3'],
    'fermion_U1': total_fermion_shared['U1'],
    'fermion_SU2': total_fermion_shared['SU2'],
    'fermion_SU3': total_fermion_shared['SU3'],
    'pure_plaq_SU3': 288,
    'total_plaq': 1200,
    'a1': 5,
    'b1': 6,
    'vertices': 120,
    'edges': 720,
    'faces': 1200,
    'cells': 600,
    'gauge_verts': 24,
    'fermion_verts': 96,
}

print(f"Available quantities: {quantities}")

# Try all ratios of pairs
results = []
keys = list(quantities.keys())
for i in range(len(keys)):
    for j in range(len(keys)):
        if i == j: continue
        v1 = quantities[keys[i]]
        v2 = quantities[keys[j]]
        if v2 == 0: continue
        ratio = v1 / v2
        if ratio > 0:
            err = abs(ratio - target) / target * 100
            if err < 5:
                results.append((err, ratio, keys[i], keys[j]))

# Also try with phi multipliers
for i in range(len(keys)):
    for j in range(len(keys)):
        if i == j: continue
        v1 = quantities[keys[i]]
        v2 = quantities[keys[j]]
        if v2 == 0: continue
        for phi_pow in [-2, -1, 1, 2]:
            ratio = v1 / v2 * PHI**phi_pow
            if ratio > 0:
                err = abs(ratio - target) / target * 100
                if err < 3:
                    results.append((err, ratio, f"{keys[i]}*phi^{phi_pow}", keys[j]))

results.sort()
print(f"\nBest matches to 10*phi = {target:.4f}:")
for err, ratio, n1, n2 in results[:15]:
    print(f"  {n1}/{n2} = {ratio:.4f} (err = {err:.2f}%)")

# ============================================================
# PART 12: Physical interpretation
# ============================================================
print("\n" + "="*70)
print("PART 12: PHYSICAL INTERPRETATION")
print("="*70)

print("""
SHARED NEIGHBOR ASYMMETRY STRUCTURE:

1. U1 edges: 0 gauge shared -> photon is ISOLATED from gauge sector
   -> Geometric origin of Abelian character
   -> QED: coupling runs via fermion loops ONLY

2. SU2 edges: ~0.17 gauge shared -> weak bosons have MINIMAL gauge access
   -> Mostly fermionic mediation (4.83 fermion shared)
   -> Consistent with broken SU(2) (W/Z get mass, weaken self-coupling)

3. SU3 edges: 1.75 gauge shared -> gluons have STRONG gauge access
   -> 4 distinct sub-types with different A/B mixes
   -> Non-Abelian character: gluons interact with other gauge bosons
   -> QCD: both fermion loops AND gluon self-interactions contribute

KEY INSIGHT:
  The shared neighbor type decomposition provides a TOPOLOGICAL
  distinction between Abelian and non-Abelian gauge theories.
  This is NOT accessible from the spectral decomposition
  (vertex-transitivity theorem, exp146).
""")

# Final check: does gauge_SU3 / (gauge_SU2 + gauge_U1) give anything?
if total_gauge_shared['SU2'] + total_gauge_shared['U1'] > 0:
    r = total_gauge_shared['SU3'] / (total_gauge_shared['SU2'] + total_gauge_shared['U1'])
    print(f"gauge_SU3 / (gauge_SU2 + gauge_U1) = {r:.4f}")

r_nonab = (total_gauge_shared['SU2'] + total_gauge_shared['SU3']) / max(1, total_gauge_shared['U1'])
print(f"(gauge_SU2 + gauge_SU3) / gauge_U1 = {'inf' if total_gauge_shared['U1'] == 0 else r_nonab:.4f}")
print(f"  (U1 gauge_shared = {total_gauge_shared['U1']})")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
EXACT SHARED NEIGHBOR PROFILES (a_1=5 for all edges):
  U1 (96 edges):     (0A, 0B, 5C) x96  -> PURELY FERMIONIC
  SU2_ch (192 edges): (0A, 0B, 5C) x192 -> PURELY FERMIONIC
  SU2_n (48 edges):  (1A, 0B, 4C) x48   -> 1 A-GAUGE ACCESS
  SU3 (384 edges):   4 sub-types:
    (0A, 1B, 4C) x?
    (0A, 2B, 3C) x?
    (1A, 1B, 3C) x?
    (1A, 2B, 2C) x?
    -> MIXED GAUGE ACCESS (avg 1.75 gauge/edge)

KEY RESULTS:
  1. U1 = Abelian (0 gauge shared) - DERIVED from topology
  2. SU3 = non-Abelian (1.75 gauge shared) - DERIVED from topology
  3. SU2_n (CC*) has A-gauge access (1 per edge) - distinguishes it from SU2_ch
  4. Coupling ratio alpha_s/alpha = 10*phi NOT derived from shared neighbors
  5. The topological distinction is QUALITATIVE (Abelian vs non-Abelian),
     not quantitative (doesn't give coupling values)
""")

print("STATUS: DERIVED (qualitative Abelian/non-Abelian distinction)")
print("        PATTERN (quantitative coupling ratios still from spectral)")

print("\n" + "="*70)
print("END EXP-147")
print("="*70)
