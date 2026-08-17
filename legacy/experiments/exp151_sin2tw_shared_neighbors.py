"""
EXP-151: sin^2(theta_W) from Shared Neighbor Structure
======================================================
Can we derive sin^2(tW) = 6/26 from the shared-neighbor asymmetry?
Connect the edge-based structure (exp147/148) to the Weinberg angle.

Parts:
1. Recap edge counts and shared-neighbor profiles
2. Various ratio combinations
3. Connection to phi-sector dimension 26
4. Per-fermion gauge budget analysis
5. Gauge boson "strength" from shared neighbors
6. Projective weight analysis
7. B-vertex role (16 = 2*8, relate to weak isospin)
8. Cross-ratios between sectors
9. Intersection number b_1=6 and its role
10. Systematic search for 6/26 from combinatorics
"""

import numpy as np
from itertools import permutations, product
from collections import defaultdict, Counter

print("="*70)
print("EXP-151: sin^2(theta_W) from Shared Neighbor Structure")
print("="*70)

PHI = (1 + np.sqrt(5)) / 2

# ============================================================
# Build 600-cell (same as exp149)
# ============================================================
def build_600cell():
    verts = set()
    for i in range(4):
        for s in [+1, -1]:
            v = [0,0,0,0]
            v[i] = 2*s
            verts.add(tuple(v))
    for signs in product([1,-1], repeat=4):
        verts.add(signs)
    base = [PHI, 1, 1/PHI, 0]
    def even_perms(lst):
        perms = []
        for p in permutations(range(4)):
            inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
            if inv % 2 == 0:
                perms.append(tuple(lst[p[i]] for i in range(4)))
        return perms
    for ep in even_perms(base):
        for s0 in ([1,-1] if ep[0]!=0 else [1]):
            for s1 in ([1,-1] if ep[1]!=0 else [1]):
                for s2 in ([1,-1] if ep[2]!=0 else [1]):
                    for s3 in ([1,-1] if ep[3]!=0 else [1]):
                        verts.add((s0*ep[0], s1*ep[1], s2*ep[2], s3*ep[3]))
    verts = sorted(verts)
    assert len(verts) == 120
    coords = np.array(verts)
    dist2 = np.sum((coords[:,None]-coords[None,:])**2, axis=2)
    edge_len2 = 4/PHI**2
    adj = (np.abs(dist2-edge_len2)<1e-10) & (dist2>0.01)
    edges = set()
    for i in range(120):
        for j in range(i+1,120):
            if adj[i,j]: edges.add((i,j))
    assert len(edges)==720
    typeA,typeB,typeC=[],[],[]
    for i in range(120):
        c=coords[i]
        nz=np.sum(np.abs(c)>0.01)
        if nz==1: typeA.append(i)
        elif nz==4 and all(abs(abs(c[k])-1)<0.01 for k in range(4)): typeB.append(i)
        else: typeC.append(i)
    return coords,edges,adj,typeA,typeB,typeC,set(typeA),set(typeB),set(typeC)

coords,edges,adj_matrix,typeA,typeB,typeC,setA,setB,setC = build_600cell()

neighbors = defaultdict(set)
for i,j in edges:
    neighbors[i].add(j)
    neighbors[j].add(i)

def get_profile(i,j):
    shared = neighbors[i]&neighbors[j]
    return (len(shared&setA), len(shared&setB), len(shared&setC))

# Classify edges
edge_types = {}
edge_profiles = {}
for i,j in edges:
    ti = 'A' if i in setA else ('B' if i in setB else 'C')
    tj = 'A' if j in setA else ('B' if j in setB else 'C')
    edge_types[(i,j)] = ''.join(sorted([ti,tj]))
    edge_profiles[(i,j)] = get_profile(i,j)

# Count by type
ac_edges = [(i,j) for i,j in edges if edge_types[(i,j)]=='AC']
bc_edges = [(i,j) for i,j in edges if edge_types[(i,j)]=='BC']
cc_edges = [(i,j) for i,j in edges if edge_types[(i,j)]=='CC']
cc_star  = [(i,j) for i,j in cc_edges if edge_profiles[(i,j)]==(1,0,4)]

# SU3 sub-types
T1 = [(i,j) for i,j in cc_edges if edge_profiles[(i,j)]==(0,1,4)]
T2 = [(i,j) for i,j in cc_edges if edge_profiles[(i,j)]==(0,2,3)]
T3 = [(i,j) for i,j in cc_edges if edge_profiles[(i,j)]==(1,1,3)]

n_U1 = len(ac_edges)  # 96
n_SU2_ch = len(bc_edges)  # 192
n_SU2_n = len(cc_star)  # 48
n_SU3 = len(T1)+len(T2)+len(T3)  # 384

print(f"\nEdge counts:")
print(f"  U(1): AC = {n_U1}")
print(f"  SU(2)_ch: BC = {n_SU2_ch}")
print(f"  SU(2)_n: CC* = {n_SU2_n}")
print(f"  SU(3): T1={len(T1)}, T2={len(T2)}, T3={len(T3)}, total={n_SU3}")
print(f"  Total: {n_U1+n_SU2_ch+n_SU2_n+n_SU3} = 720")
print(f"\n  Per fermion: {n_U1//96} AC + {n_SU2_ch//96} BC + {n_SU2_n//96} CC* + {n_SU3//96} CC = {(n_U1+n_SU2_ch+n_SU2_n+n_SU3)//96}")

# ============================================================
# Part 1: Basic ratios
# ============================================================
print("\n" + "="*70)
print("PART 1: Basic edge count ratios")
print("="*70)

target = 6/26
print(f"Target: sin^2(tW) = 6/26 = {target:.6f}")
print(f"Experimental: 0.23122")

# Try all simple ratios
print(f"\nSimple ratios:")
quantities = {
    'n_U1': n_U1, 'n_SU2': n_SU2_ch+n_SU2_n, 'n_SU2_ch': n_SU2_ch,
    'n_SU2_n': n_SU2_n, 'n_SU3': n_SU3, 'n_total': 720,
    '|A|': 8, '|B|': 16, '|C|': 96, 'N': 120,
    'edges_per_C': 12, 'degree': 12, 'dim': 4, 'diameter': 5,
    'T1': len(T1), 'T2': len(T2), 'T3': len(T3),
    'a1': 5, 'b1': 6, 'phi_dim': 26,
}

close_ratios = []
for name1, v1 in quantities.items():
    for name2, v2 in quantities.items():
        if v2 != 0 and v1 != v2:
            r = v1/v2
            err = abs(r - target)/target * 100
            if err < 5:
                close_ratios.append((name1, name2, r, err))

close_ratios.sort(key=lambda x: x[3])
for n1, n2, r, err in close_ratios[:15]:
    print(f"  {n1}/{n2} = {quantities[n1]}/{quantities[n2]} = {r:.6f} (err {err:.3f}%)")

# ============================================================
# Part 2: Known derivation path: b_1/phi_dim = 6/26
# ============================================================
print("\n" + "="*70)
print("PART 2: Known spectral derivation")
print("="*70)

print(f"b_1 = 6 (second intersection number of 600-cell)")
print(f"phi-sector dim = 4+9+9+4 = 26 (multiplicities of phi-eigenvalues)")
print(f"sin^2(tW) = b_1/phi_dim = 6/26 = {6/26:.6f}")
print(f"\nThis is SPECTRAL. Can we find a TOPOLOGICAL path?")

# ============================================================
# Part 3: Shared neighbor based derivation attempts
# ============================================================
print("\n" + "="*70)
print("PART 3: Shared neighbor approaches")
print("="*70)

# Per-fermion gauge budget: 0(U1) + 1(SU2) + 14(SU3) = 15
# Per-edge shared C neighbors: U1=5, SU2_ch=5, SU2_n=4, SU3_avg=10/3
print("Per-fermion analysis:")
print(f"  Gauge budget: 0 U1 + 1 SU2 + 14 SU3 = 15")
print(f"  Edge budget: 1 AC + 2 BC + 0.5 CC* + 4 CC = 7.5 per fermion")

# SU2 gauge shared / total gauge shared
su2_gauge = 1  # per fermion
total_gauge = 15
print(f"\n  SU2_gauge/total_gauge = {su2_gauge}/{total_gauge} = {su2_gauge/total_gauge:.6f}")
print(f"  (not 6/26)")

# Try: weak sector edges / total non-trivial
weak = n_SU2_ch + n_SU2_n  # 240
strong = n_SU3  # 384
em = n_U1  # 96
print(f"\n  weak/(weak+strong) = {weak}/({weak}+{strong}) = {weak/(weak+strong):.6f}")
print(f"  weak/(weak+strong+em) = {weak}/{720} = {weak/720:.6f}")
print(f"  em/(em+weak) = {em}/({em}+{weak}) = {em/(em+weak):.6f}")

# Connection through shared B vertices
# SU2_n shares 1A, 0B. SU2_ch shares 0A, 0B.
# B-vertices are weak isospin doublets?
print(f"\n  |B|/(|A|+|B|) = 16/(8+16) = {16/24:.6f}")
print(f"  |A|/(|A|+|B|+|C|) = 8/120 = {8/120:.6f}")
print(f"  |B|/|C| = 16/96 = {16/96:.6f}")

# ============================================================
# Part 4: Combinatorial formulas for 6/26
# ============================================================
print("\n" + "="*70)
print("PART 4: Combinatorial search for 6/26")
print("="*70)

# 6/26 = 3/13
# Can we get 3 and 13 from the structure?
print("6/26 = 3/13")
print(f"3 = shared gauge per fermion (1A+2B) per SU2_n profile (1,0,4)")
print(f"  Actually per fermion: 3 gauge neighbors (1A + 2B)")
print(f"13 = ?")

# From exp147: per-fermion breakdown
# Each C vertex has: 1 A-neighbor, 2 B-neighbors, 9 C-neighbors = 12 total
# 3 gauge neighbors out of 12 total = 3/12 = 1/4
# Not 6/26

# Shared C neighbors across all edge types per fermion:
# U1: 5C shared, SU2_ch: 5C shared, SU2_n: 4C shared
# SU3: T1=4C, T2=3C, T3=3C (weighted avg)
# Total shared C per fermion via all 12 edges?

# Per fermion: 1*5(U1) + 2*5(SU2_ch) + 0.5*4(SU2_n) + 1*4(T1) + 1*3(T2) + 2*3(T3) = ?
shared_C_per_fermion = 1*5 + 2*5 + 0.5*4 + 1*4 + 1*3 + 2*3
print(f"\nShared C per fermion (all edges): {shared_C_per_fermion}")

# Total shared neighbors per fermion across all 12 edges:
# U1=(0,0,5): 5; SU2_ch=(0,0,5): 5*2=10; SU2_n=(1,0,4): 5*0.5=2.5
# T1=(0,1,4): 5*1=5; T2=(0,2,3): 5*1=5; T3=(1,1,3): 5*2=10
total_shared = 1*5 + 2*5 + 0.5*5 + 1*5 + 1*5 + 2*5
print(f"Total shared per fermion: {total_shared}")

# Gauge shared per fermion = 0+0+1+1+2+4 = 8? No...
# From exp147: per-fermion gauge budget = 15 (total gauge shared counting multiplicity)
# Let me compute more carefully
print(f"\nPer-fermion gauge budget breakdown:")
print(f"  1 AC edge x (0A+0B shared gauge) = 0")
print(f"  2 BC edges x (0A+0B shared gauge) = 0")
print(f"  0.5 CC* edges x (1A+0B=1 gauge) = 0.5")
# Wait, CC* has profile (1,0,4) -> shared 1A, 0B
# But per fermion we have exactly 1 CC* edge (from exp143: 48 total / 96 C verts...
# actually each CC* has 2 endpoints so 48*2/96 = 1 per fermion)
print(f"  Correction: each C vertex has exactly 1 CC* edge (48 edges, 96 C-verts, 2 endpoints)")
# Wait actually 48 CC* edges with 2 endpoints each = 96 endpoint slots / 96 C verts = 1.0
# So per fermion: 1 CC* edge sharing 1A = 1 gauge shared from CC*

# T3: 192 edges, (1A,1B,3C) -> gauge = 2 per edge
# Per fermion: T3 count = 192*2/96 = 4 T3 edges per fermion
# T1: 96*2/96 = 2 T1 edges per fermion, gauge = 1B per edge -> 2 gauge
# T2: 96*2/96 = 2 T2 edges per fermion, gauge = 2B per edge -> 4 gauge

gauge_from_CCstar = 1 * 1  # 1 edge * 1 gauge shared (1A)
gauge_from_T1 = 2 * 1  # 2 edges * 1B shared
gauge_from_T2 = 2 * 2  # 2 edges * 2B shared
gauge_from_T3 = 4 * 2  # 4 edges * (1A+1B) shared
gauge_total = gauge_from_CCstar + gauge_from_T1 + gauge_from_T2 + gauge_from_T3
print(f"\n  CC* (1 edge x 1 gauge): {gauge_from_CCstar}")
print(f"  T1 (2 edges x 1 gauge): {gauge_from_T1}")
print(f"  T2 (2 edges x 2 gauge): {gauge_from_T2}")
print(f"  T3 (4 edges x 2 gauge): {gauge_from_T3}")
print(f"  Total: {gauge_total}")

# Per-fermion gauge shared split by A/B
gauge_A = 1*1 + 4*1  # CC* shares 1A, T3 shares 1A -> 1+4=5A
gauge_B = 2*1 + 2*2 + 4*1  # T1 shares 1B, T2 shares 2B, T3 shares 1B -> 2+4+4=10B
print(f"\n  Gauge-A shared per fermion: {gauge_A}")
print(f"  Gauge-B shared per fermion: {gauge_B}")
print(f"  Total: {gauge_A + gauge_B}")

# Can we get 6/26?
print(f"\n  gauge_A / gauge_total = {gauge_A}/{gauge_A+gauge_B} = {gauge_A/(gauge_A+gauge_B):.6f}")
print(f"  gauge_B / gauge_total = {gauge_B}/{gauge_A+gauge_B} = {gauge_B/(gauge_A+gauge_B):.6f}")
print(f"  6/26 = {6/26:.6f}")

# Hmm, 5/15 = 1/3 and 10/15 = 2/3

# ============================================================
# Part 5: b_1 from topology
# ============================================================
print("\n" + "="*70)
print("PART 5: b_1 = 6 from topological perspective")
print("="*70)

# b_1 is the number of common neighbors at distance 2
# For 600-cell: if d(u,v)=2, then |N(u) cap N(v)| = b_1 = 6
# Is b_1 related to shared neighbor structure?

# Compute b_1 directly
sample_pairs = []
for i in range(120):
    for j in range(i+1, 120):
        if adj_matrix[i,j] == 0:  # not neighbors
            shared = neighbors[i] & neighbors[j]
            if len(shared) > 0:  # distance 2
                sample_pairs.append((i,j, len(shared)))
                if len(sample_pairs) > 100:
                    break
    if len(sample_pairs) > 100:
        break

b1_vals = [x[2] for x in sample_pairs]
print(f"b_1 (distance-2 common neighbors): {sorted(set(b1_vals))}")
# Should be uniformly 6 (distance-regular)

# 2-hop gauge access (from exp147): total ~ 6 = b_1
print(f"\n2-hop gauge access per edge ~ 6 = b_1 (from exp147)")
print(f"This is because at distance 2, every pair shares b_1=6 neighbors")
print(f"b_1 is a global topological invariant of the 600-cell")

# ============================================================
# Part 6: phi-sector from topology
# ============================================================
print("\n" + "="*70)
print("PART 6: Can phi-sector dimension 26 be derived topologically?")
print("="*70)

# phi-sector = eigenspaces with phi-irrational eigenvalues
# Multiplicities: 4+9+9+4 = 26
# Rational sector: 1+16+25+36+16 = 94 (but not 94!)
# Wait: 9 eigenvalues with mults 1,4,9,16,25,36,9,16,4
# phi-eigenvalues: lambda=6phi(m=4), 4phi(m=9), 4-4phi(m=9), 6-6phi(m=4) => dim=26
# Rational: lambda=12(m=1), 3(m=16), 0(m=25), -2(m=36), -3(m=16) => dim=94
# 26+94=120 check

print(f"Phi-sector: eigenspaces with irrational (phi-dependent) eigenvalues")
print(f"  6*phi (m=4), 4*phi (m=9), 4-4*phi (m=9), 6-6*phi (m=4)")
print(f"  Total: 4+9+9+4 = 26")
print(f"\nRational sector:")
print(f"  12 (m=1), 3 (m=16), 0 (m=25), -2 (m=36), -3 (m=16)")
print(f"  Total: 1+16+25+36+16 = 94")
print(f"\n26 + 94 = 120 (check)")

# 26 is a topological invariant (doesn't depend on labeling)
# But it comes from the SPECTRUM, not directly from edge counting

# Connection: 26 = |A|+|B|+2? No, 8+16+2=26!
print(f"\n26 = |A| + |B| + 2 = 8 + 16 + 2 = 26")
print(f"What is the +2?")
print(f"  26 - |A| - |B| = 26 - 24 = 2")
print(f"  |A| + |B| = 24 = |gauge vertices|")
print(f"  phi_dim = |gauge| + 2")

# Or: 26 = 120 - 94 = N - |rational sector|
# Or: 26 eigenvalue slots out of 120

# Try: 6/26 = b_1/(N-94) = b_1/(|A|+|B|+2)
print(f"\n6/26 = b_1 / phi_dim = b_1 / (|gauge| + 2)")
print(f"     = 6 / (24 + 2)")
print(f"     = 6 / 26")

# Can we understand the "+2"?
# 26 = dim of phi-sector. Galois pairs: (6phi, 6-6phi) and (4phi, 4-4phi)
# 4+4 = 8, 9+9 = 18. 8+18=26.
# 8 = |A|, 18 = |B| + 2? No, 18 != 18.
# Actually: mult 4 eigenspaces: 2 of them (one phi, one conjugate)
# mult 9 eigenspaces: 2 of them. Total = 2*4 + 2*9 = 26

print(f"\nGalois structure:")
print(f"  Pair 1: lambda=6*phi (m=4) <-> lambda=6-6*phi (m=4)")
print(f"  Pair 2: lambda=4*phi (m=9) <-> lambda=4-4*phi (m=9)")
print(f"  phi-sector = 2*(4+9) = 2*13 = 26")
print(f"\n6/26 = 6/(2*13) = 3/13")
print(f"  3 = gauge neighbors per fermion (1A + 2B)")
print(f"  13 = sum of Galois-paired multiplicities (4+9)")

# ============================================================
# Part 7: The 3/13 decomposition
# ============================================================
print("\n" + "="*70)
print("PART 7: sin^2(tW) = 3/13 deep structure")
print("="*70)

print(f"sin^2(tW) = 6/26 = 3/13")
print(f"\nNumerator 3:")
print(f"  - Each fermion has 3 gauge neighbors (1A + 2B)")
print(f"  - M_vacancy = 3 (from exp109)")
print(f"  - F/E = 1200/720 = 5/3 -> denominator")
print(f"  - 3 = number of gauge nbrs per fermion EXACTLY")

print(f"\nDenominator 13:")
print(f"  - 13 = sum of paired multiplicities (4+9)")
print(f"  - 13 is a Fibonacci number (F_7)")
print(f"  - 13 = |A| + |B|/2 - 1 = 8 + 8 - 3 = 13... no")
print(f"  - 13*2 = 26 = phi-sector dim")

# Connection to gauge budget
# Per fermion: 3 gauge neighbors, 9 C neighbors = 12 total
# Weak: 3 gauge / (3+9+1?)
print(f"\nPer-fermion: 3 gauge + 9 fermionic = 12 total neighbors")
print(f"  gauge/total = 3/12 = 1/4 = 0.250 (close to 6/26=0.231 but off)")

# Actually can we get 13 from geometry?
# 13 = number of distinct triangulation types? No (3 types)
# 13 = some vertex count in a substructure?

# Check: 4+9 = 13. These are m=4^2 and m=3^2. Dimensions of H4 irreps.
# dim_1 + dim_2 of the phi-sector H4 irreps
print(f"\n4 = 2^2 (H4 irrep dim for Galois pair 1)")
print(f"9 = 3^2 (H4 irrep dim for Galois pair 2)")
print(f"4+9 = 13 = phi-sector 'one-sided' dimension")

# ============================================================
# Part 8: Systematic formula search
# ============================================================
print("\n" + "="*70)
print("PART 8: Systematic search for 6/26")
print("="*70)

# Use all known integers from the 600-cell
vals = {
    '|A|': 8, '|B|': 16, '|C|': 96, 'N': 120,
    'E': 720, 'F': 1200, 'T': 600,
    'deg': 12, 'diam': 5,
    'a1': 5, 'b1': 6, 'c1': 1,  # intersection numbers
    'a2': -18, 'b2': -6, 'c2': 25,
    'gauge_per_f': 3,
    'CC_per_f': 9,
    'n_U1': 96, 'n_SU2': 240, 'n_SU3': 384,
    'n_CC*': 48, 'n_T1': 96, 'n_T2': 96, 'n_T3': 192,
    '|H4|': 14400,
    '24cell': 24,
    'phi_dim': 26,
}

# Search: (a*x + b*y) / (c*w + d*z) = 6/26 where a,b,c,d in {-2,-1,0,1,2}
print("Looking for exact 6/26 = 3/13 from simple combinations...")
found = []
keys = list(vals.keys())
for i, k1 in enumerate(keys):
    for j, k2 in enumerate(keys):
        if i == j: continue
        num = vals[k1]
        den = vals[k2]
        if den != 0 and num != 0:
            r = num/den
            if abs(r - 6/26) < 0.001:
                found.append(f"{k1}/{k2} = {num}/{den} = {r:.6f}")

# Also try a/b + c/d and a*b/(c*d)
for i, k1 in enumerate(keys):
    for j, k2 in enumerate(keys):
        if i >= j: continue
        prod = vals[k1] * vals[k2]
        for m, k3 in enumerate(keys):
            for n, k4 in enumerate(keys):
                if m >= n: continue
                den = vals[k3] * vals[k4]
                if den != 0:
                    r = prod / den
                    if abs(r - 6/26) < 0.0001 and r != 1:
                        found.append(f"{k1}*{k2}/({k3}*{k4}) = {vals[k1]}*{vals[k2]}/({vals[k3]}*{vals[k4]}) = {r:.6f}")

for f in sorted(set(found)):
    print(f"  {f}")

# ============================================================
# Part 9: The b_1 connection
# ============================================================
print("\n" + "="*70)
print("PART 9: b_1 = 6 as the numerator")
print("="*70)

print(f"b_1 = 6 is the SECOND intersection number of the 600-cell")
print(f"Physical meaning: at distance 2, any two vertices share 6 common neighbors")
print(f"\nIn the formula n = a_1*a + b_1*b = 5a + 6b:")
print(f"  a_1 = 5 controls the 'diameter' contribution")
print(f"  b_1 = 6 controls the 'decagon' contribution")
print(f"\nsin^2(tW) = b_1/phi_dim = 6/26")
print(f"  b_1 measures LATERAL connectivity (distance-2 sharing)")
print(f"  phi_dim measures SPECTRAL extent (phi-eigenvalue space)")

# Another path: b_1 = 6 = 2*3 = 2 * (gauge nbrs per fermion)
print(f"\nb_1 = 2 * (gauge_per_fermion) = 2 * 3 = 6")
print(f"phi_dim = 2 * 13")
print(f"So: sin^2(tW) = (2*gauge_per_f) / (2*13) = gauge_per_f / 13")

# ============================================================
# Part 10: Cross-check with vertex counting
# ============================================================
print("\n" + "="*70)
print("PART 10: New topological path?")
print("="*70)

# SU(2) has 3 generators. In the edge model: 1 AC + 2 BC + 1 CC* = ...
# Actually: SU(2) = 2 BC + 1 CC* = 3 edges per fermion
# SU(3) = 8 CC edges per fermion
# U(1) = 1 AC edge per fermion

# Total non-U(1) per fermion = 3 + 8 = 11
# Total = 12
# U(1) fraction = 1/12

# Weinberg angle: sin^2(tW) measures the U(1) component of the photon
# In SM: photon = cos(tW)*B + sin(tW)*W_3
# sin^2(tW) = g'^2/(g^2+g'^2)

# In edge model:
# Electroweak = U(1) + SU(2) = 1 + 3 = 4 edges per fermion
# sin^2(tW) should measure U(1)/electroweak ratio somehow

print(f"Electroweak edges per fermion: 1(U1) + 3(SU2) = 4")
print(f"Strong edges per fermion: 8(SU3)")
print(f"Total: 12")
print(f"\nRatios:")
print(f"  U1/EW = 1/4 = 0.250 (close but wrong)")
print(f"  U1/total = 1/12 = 0.083 (too small)")
print(f"  EW/total = 4/12 = 1/3 (too big)")

# GUT normalization factor: 5/3
# sin^2(tW)_GUT = 3/8 (at unification)
# sin^2(tW)_phys = (3/5)*g1^2/(g1^2 + g2^2)
# With GUT norm: sin^2(tW) = (3/5)/(1+g2^2/g1^2)

# From edge model with GUT normalization:
# g1^2 ~ n_U1 = 96, g2^2 ~ n_SU2 = 240
# sin^2 = (3/5) * 96 / (96 + 240) = (3/5) * 96/336 = (3/5)*(2/7) = 6/35 = 0.171
# Nope

# Or: hypercharge edges = AC edges carry hypercharge
# Weak edges = BC + CC* carry weak isospin
# GUT normalization: multiply hypercharge by sqrt(5/3)

# Try: n_U1 * (5/3) / (n_U1*(5/3) + n_SU2)
gut_sin2 = n_U1*(5/3) / (n_U1*(5/3) + n_SU2_ch + n_SU2_n)
print(f"\nGUT-normalized: n_U1*(5/3) / (n_U1*(5/3) + n_SU2)")
print(f"  = 96*(5/3) / (160 + 240) = {96*5/3:.1f} / {96*5/3 + 240:.1f} = {gut_sin2:.6f}")
print(f"  Target: {target:.6f}, err: {abs(gut_sin2-target)/target*100:.2f}%")

# With just charged SU2:
gut_sin2b = n_U1*(5/3) / (n_U1*(5/3) + n_SU2_ch)
print(f"\n  Only SU2_ch: {96*5/3:.1f} / ({96*5/3:.1f} + {n_SU2_ch}) = {gut_sin2b:.6f}")

# Per fermion with GUT norm:
pf_gut = 1*(5/3) / (1*(5/3) + 3)
print(f"\n  Per fermion: 1*(5/3) / (5/3 + 3) = {5/3:.4f}/{5/3+3:.4f} = {pf_gut:.6f}")
print(f"  Target: {target:.6f}, err: {abs(pf_gut-target)/target*100:.2f}%")

# 5/3 * 1 / (5/3 + 3) = 5/3 / 14/3 = 5/14 = 0.357 nope

# Try with 600-cell specific factor
# sin^2(tW) = U1_weight / (U1_weight + SU2_weight)
# where weights come from shared-neighbor structure

# From exp147:
# 10*phi = 480/48 * phi = (fermion_shared_U1/gauge_shared_SU2) * phi
# alpha_s/alpha = 10*phi
# sin^2(tW) = alpha/alpha_2 where alpha_2 = alpha/sin^2(tW)
# In terms of coupling ratios: sin^2(tW) = alpha_em/alpha_2

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
sin^2(tW) = 6/26 derivation paths:

PATH 1 (SPECTRAL - known from exp111):
  b_1/phi_dim = 6/26
  b_1 = second intersection number
  phi_dim = dim of phi-irrational eigenvalue sector

PATH 2 (SEMI-TOPOLOGICAL - new):
  sin^2(tW) = (2 * gauge_per_fermion) / (2 * (m_Galois_1 + m_Galois_2))
           = (2 * 3) / (2 * (4 + 9))
           = 6 / 26
  gauge_per_fermion = 3 (topological: 1A + 2B per C vertex)
  4, 9 = multiplicities of phi-Galois paired eigenspaces (spectral)

PATH 3 (NOT FOUND):
  Pure topological (edge-counting only) derivation NOT achieved.
  The denominator 26 requires spectral information.

KEY INSIGHT:
  The numerator b_1 = 6 = 2 * (gauge_nbrs_per_fermion)
  This connects b_1 to the A/B/C vertex structure topologically.
  But 26 = phi-sector dim remains spectral.

STATUS: PATTERN (second path semi-topological, no pure topological derivation)
""")
