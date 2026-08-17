# EXP-225b: QUANTITATIVE TEST OF GALOIS COMPLEMENTARITY
# ======================================================================
# Build the 600-cell, decompose A = A_rat + A_phi + A_phi' via spectral
# projectors, classify edges by coset structure, and test whether
# different gauge channels project onto different Galois sectors.
# ======================================================================

import numpy as np
from itertools import product

PHI = (1 + np.sqrt(5)) / 2
a_1 = 5; b_1 = 6; deg = 12; N_vert = 120
N_bag = 13; h_E8 = 30

print("="*70)
print("EXP-225b: QUANTITATIVE GALOIS SECTOR TEST")
print("="*70)

# ===================================================================
# PART A: BUILD 600-CELL
# ===================================================================
print("\n--- PART A: Build 600-cell ---")

vertices = []

# Type A: 8 vertices - permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        vertices.append(tuple(v))

# Type B: 16 vertices - (+-1/2, +-1/2, +-1/2, +-1/2)
for signs in product([0.5, -0.5], repeat=4):
    vertices.append(signs)

# Type C: 96 vertices - even permutations of (0, +-1/2, +-phi/2, +-1/(2*phi))
base_vals = [0.0, 0.5, PHI/2, 1/(2*PHI)]
# Even permutations of (0,1,2,3): (0123),(0231),(0312),(1032),(1203),(1320),
# (2013),(2130),(2301),(3021),(3102),(3210)
# Actually there are 12 even permutations of 4 elements
even_perms = [
    (0,1,2,3), (0,2,3,1), (0,3,1,2),
    (1,0,3,2), (1,2,0,3), (1,3,2,0),
    (2,0,1,3), (2,1,3,0), (2,3,0,1),
    (3,0,2,1), (3,1,0,2), (3,2,1,0)
]

for perm in even_perms:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                # Apply signs to non-zero positions
                v = [0.0]*4
                signs = [1, s1, s2, s3]
                for i in range(4):
                    v[perm[i]] = signs[i] * base_vals[i]
                vt = tuple(v)
                # Check unit norm
                norm = sum(x*x for x in vt)
                if abs(norm - 1.0) < 0.01:
                    # Check if already present
                    is_dup = False
                    for existing in vertices:
                        if all(abs(vt[k] - existing[k]) < 1e-6 for k in range(4)):
                            is_dup = True
                            break
                    if not is_dup:
                        vertices.append(vt)

print(f"  Total vertices: {len(vertices)}")

# Classify vertex types
type_A = []
type_B = []
type_C = []
for i, v in enumerate(vertices):
    n_zero = sum(1 for x in v if abs(x) < 0.01)
    if n_zero == 3:
        type_A.append(i)
    elif n_zero == 0 and all(abs(abs(x) - 0.5) < 0.01 for x in v):
        type_B.append(i)
    else:
        type_C.append(i)

print(f"  Type A: {len(type_A)}, Type B: {len(type_B)}, Type C: {len(type_C)}")

# Build adjacency matrix
verts = np.array(vertices)
inner = verts @ verts.T
adj = (inner > PHI/2 - 0.01).astype(int)
np.fill_diagonal(adj, 0)

degrees = adj.sum(axis=1)
print(f"  Degree check: min={degrees.min()}, max={degrees.max()}, mean={degrees.mean():.1f}")
n_edges = adj.sum() // 2
print(f"  Edges: {n_edges}")

if len(vertices) != 120 or n_edges != 720:
    print("  WARNING: vertex/edge count mismatch! Adjusting threshold...")
    # Try different thresholds
    for thresh in [0.805, 0.808, 0.809, 0.810, 0.8085, 0.8090]:
        adj_test = (inner > thresh).astype(int)
        np.fill_diagonal(adj_test, 0)
        ne = adj_test.sum() // 2
        degs = adj_test.sum(axis=1)
        if ne == 720 and degs.min() == 12 and degs.max() == 12:
            adj = adj_test
            n_edges = ne
            print(f"  Fixed with threshold {thresh}: {ne} edges, degree {degs.min()}")
            break

# ===================================================================
# PART B: COSET STRUCTURE (2T = A+B, 4 C-cosets)
# ===================================================================
print(f"\n--- PART B: Coset structure ---")

# The 2T coset = {A vertices} + {B vertices} = 24 elements
coset_2T = set(type_A + type_B)
print(f"  2T coset (A+B): {len(coset_2T)} vertices")

# Find 4 C-cosets using quaternion multiplication
# q1^{-1} * q2 in 2T <=> same coset
# For unit quaternion: q^{-1} = conjugate = (a, -b, -c, -d)

def quat_mult(q1, q2):
    """Multiply two quaternions q = (a,b,c,d) = a + bi + cj + dk"""
    a1,b1,c1,d1 = q1
    a2,b2,c2,d2 = q2
    return (a1*a2 - b1*b2 - c1*c2 - d1*d2,
            a1*b2 + b1*a2 + c1*d2 - d1*c2,
            a1*c2 - b1*d2 + c1*a2 + d1*b2,
            a1*d2 + b1*c2 - c1*b2 + d1*a2)

def quat_conj(q):
    return (q[0], -q[1], -q[2], -q[3])

def is_2T_element(q, tol=0.01):
    """Check if quaternion q is in 2T = {+-1,+-i,+-j,+-k, (+-1+-i+-j+-k)/2}"""
    a,b,c,d = q
    # Check type A: one component +-1, rest 0
    for idx in range(4):
        coords = [a,b,c,d]
        if abs(abs(coords[idx]) - 1.0) < tol:
            others = [coords[j] for j in range(4) if j != idx]
            if all(abs(x) < tol for x in others):
                return True
    # Check type B: all +-1/2
    if all(abs(abs(x) - 0.5) < tol for x in [a,b,c,d]):
        return True
    return False

# Classify C vertices into cosets
c_cosets = {}  # vertex_idx -> coset_id
coset_id = 0
unassigned = set(type_C)

while unassigned:
    # Pick a vertex
    v0_idx = min(unassigned)
    v0 = vertices[v0_idx]
    current_coset = set()

    for c_idx in list(unassigned):
        vc = vertices[c_idx]
        # Check if v0^{-1} * vc in 2T
        q_prod = quat_mult(quat_conj(v0), vc)
        if is_2T_element(q_prod):
            current_coset.add(c_idx)
            c_cosets[c_idx] = coset_id

    unassigned -= current_coset
    print(f"  C-coset {coset_id}: {len(current_coset)} vertices")
    coset_id += 1

# Full coset assignment
coset_of = {}
for idx in coset_2T:
    coset_of[idx] = 0  # 2T coset = coset 0
for idx, cid in c_cosets.items():
    coset_of[idx] = cid + 1  # C cosets = 1,2,3,4

print(f"  Total cosets: {max(coset_of.values())+1}")
for cid in range(5):
    count = sum(1 for v in coset_of.values() if v == cid)
    print(f"    Coset {cid}: {count} vertices")

# ===================================================================
# PART C: EDGE CLASSIFICATION BY COSET STRUCTURE
# ===================================================================
print(f"\n--- PART C: Edge classification ---")

# Count edges by (coset_i, coset_j) pair
from collections import Counter

edge_type_count = Counter()
edge_list_by_type = {}

for i in range(120):
    for j in range(i+1, 120):
        if adj[i,j]:
            ci, cj = coset_of[i], coset_of[j]
            key = (min(ci,cj), max(ci,cj))
            edge_type_count[key] += 1
            if key not in edge_list_by_type:
                edge_list_by_type[key] = []
            edge_list_by_type[key].append((i,j))

print(f"  Edge types by coset pair:")
intra_total = 0
inter_total = 0
for key in sorted(edge_type_count.keys()):
    count = edge_type_count[key]
    label = "INTRA" if key[0] == key[1] else "INTER"
    if key[0] == key[1]:
        intra_total += count
    else:
        inter_total += count
    print(f"    ({key[0]},{key[1]}): {count:>4d} edges  [{label}]")

print(f"  Total intra-coset: {intra_total}")
print(f"  Total inter-coset: {inter_total}")
print(f"  Total: {intra_total + inter_total}")

# Classify by vertex type pairs
vtype = {}
for i in type_A: vtype[i] = 'A'
for i in type_B: vtype[i] = 'B'
for i in type_C: vtype[i] = 'C'

edge_vtype_count = Counter()
for i in range(120):
    for j in range(i+1, 120):
        if adj[i,j]:
            t1, t2 = sorted([vtype[i], vtype[j]])
            edge_vtype_count[(t1,t2)] += 1

print(f"\n  Edge types by vertex type:")
for key in sorted(edge_vtype_count.keys()):
    print(f"    {key[0]}-{key[1]}: {edge_vtype_count[key]:>4d} edges")

# ===================================================================
# PART D: EIGENSYSTEM AND GALOIS DECOMPOSITION
# ===================================================================
print(f"\n--- PART D: Spectral decomposition ---")

eigenvalues_full, eigenvectors = np.linalg.eigh(adj.astype(float))

# Sort by eigenvalue (descending)
idx_sort = np.argsort(-eigenvalues_full)
eigenvalues_full = eigenvalues_full[idx_sort]
eigenvectors = eigenvectors[:, idx_sort]

# Identify the 9 distinct eigenvalues
from collections import defaultdict
eig_groups = defaultdict(list)
for i, ev in enumerate(eigenvalues_full):
    found = False
    for key in eig_groups:
        if abs(ev - key) < 0.01:
            eig_groups[key].append(i)
            found = True
            break
    if not found:
        eig_groups[ev] = [i]

print(f"  Distinct eigenvalues: {len(eig_groups)}")

# Known eigenvalues and their Galois sectors
# PHI eigenvalues: 6*phi (~9.708), 4*phi (~6.472)
# PHI' eigenvalues: 4-4*phi (~-2.472), 6-6*phi (~-3.708)
# Rational: 12, 3, 0, -2, -3

def galois_sector(lam):
    if abs(lam - round(lam)) < 0.05:
        return "RAT"
    elif lam > 0:
        return "PHI"
    else:
        return "PHI'"

print(f"\n  Eigenvalue classification:")
sector_indices = {"RAT": [], "PHI": [], "PHI'": []}
for ev in sorted(eig_groups.keys(), reverse=True):
    indices = eig_groups[ev]
    sector = galois_sector(ev)
    sector_indices[sector].extend(indices)
    print(f"    lambda = {ev:>8.4f}, mult = {len(indices):>2d}, sector = {sector}")

print(f"\n  Sector mode counts:")
for sector in ["RAT", "PHI", "PHI'"]:
    print(f"    {sector:>4s}: {len(sector_indices[sector]):>3d} modes")

# ===================================================================
# PART E: SPECTRAL PROJECTORS
# ===================================================================
print(f"\n--- PART E: Spectral projectors ---")

# Build A_rat, A_phi, A_phi' by summing lambda_k * |v_k><v_k|
# for eigenvalues in each sector

# Actually, we want the PROJECTOR P_sector = sum |v_k><v_k| (no lambda)
# And also the WEIGHTED projector A_sector = sum lambda_k |v_k><v_k|

P_rat = np.zeros((120, 120))
P_phi = np.zeros((120, 120))
P_phip = np.zeros((120, 120))

A_rat = np.zeros((120, 120))
A_phi = np.zeros((120, 120))
A_phip = np.zeros((120, 120))

for sector, P, A_sec in [("RAT", P_rat, A_rat), ("PHI", P_phi, A_phi), ("PHI'", P_phip, A_phip)]:
    for idx in sector_indices[sector]:
        v = eigenvectors[:, idx]
        P += np.outer(v, v)
        A_sec += eigenvalues_full[idx] * np.outer(v, v)

# Verify: P_rat + P_phi + P_phip = I
P_total = P_rat + P_phi + P_phip
print(f"  P_rat + P_phi + P_phip = I? max|P_total - I| = {np.max(np.abs(P_total - np.eye(120))):.2e}")

# Verify: A_rat + A_phi + A_phip = A
A_total = A_rat + A_phi + A_phip
print(f"  A_rat + A_phi + A_phip = A? max|A_total - A| = {np.max(np.abs(A_total - adj)):.2e}")

# ===================================================================
# PART F: GALOIS SECTOR WEIGHTS PER EDGE
# ===================================================================
print(f"\n--- PART F: Galois sector weights per edge ---")

print("""
  For each edge (i,j) with A_{ij}=1:
  A_rat(i,j) + A_phi(i,j) + A_phip(i,j) = 1
  These give the "Galois sector decomposition" of each edge.
""")

# Compute sector weights for each edge
edge_data = []
for i in range(120):
    for j in range(i+1, 120):
        if adj[i,j]:
            w_rat = A_rat[i,j]
            w_phi = A_phi[i,j]
            w_phip = A_phip[i,j]
            ci, cj = coset_of[i], coset_of[j]
            same_coset = (ci == cj)
            t1, t2 = sorted([vtype[i], vtype[j]])
            edge_data.append({
                'i': i, 'j': j,
                'w_rat': w_rat, 'w_phi': w_phi, 'w_phip': w_phip,
                'coset_pair': (min(ci,cj), max(ci,cj)),
                'same_coset': same_coset,
                'vtype': (t1, t2)
            })

# Verify sum = 1 for each edge
sums = [e['w_rat'] + e['w_phi'] + e['w_phip'] for e in edge_data]
print(f"  Edge weight sum: mean = {np.mean(sums):.6f}, std = {np.std(sums):.2e}")

# ===================================================================
# PART G: SECTOR WEIGHTS BY EDGE TYPE
# ===================================================================
print(f"\n--- PART G: Average sector weights by edge type ---")

# Group by same/different coset
intra_edges = [e for e in edge_data if e['same_coset']]
inter_edges = [e for e in edge_data if not e['same_coset']]

def print_sector_stats(edges, label):
    if not edges:
        print(f"  {label}: no edges")
        return
    w_rat = np.mean([e['w_rat'] for e in edges])
    w_phi = np.mean([e['w_phi'] for e in edges])
    w_phip = np.mean([e['w_phip'] for e in edges])
    print(f"  {label:>25s} ({len(edges):>3d} edges): "
          f"RAT={w_rat:>7.4f}  PHI={w_phi:>7.4f}  PHI'={w_phip:>7.4f}")

print_sector_stats(intra_edges, "INTRA-COSET")
print_sector_stats(inter_edges, "INTER-COSET")

print()

# Group by vertex type pair
for vt in sorted(set(e['vtype'] for e in edge_data)):
    subset = [e for e in edge_data if e['vtype'] == vt]
    label = f"{vt[0]}-{vt[1]}"
    print_sector_stats(subset, label)

print()

# Group by coset pair
for cp in sorted(set(e['coset_pair'] for e in edge_data)):
    subset = [e for e in edge_data if e['coset_pair'] == cp]
    same = "intra" if cp[0]==cp[1] else "inter"
    label = f"coset({cp[0]},{cp[1]}) [{same}]"
    print_sector_stats(subset, label)

# ===================================================================
# PART H: INTRA-COSET SUBSTRUCTURE
# ===================================================================
print(f"\n--- PART H: Intra-coset sub-classification ---")

# Within the 2T coset (coset 0):
coset0_edges = [e for e in edge_data if e['coset_pair'] == (0,0)]
print_sector_stats(coset0_edges, "2T coset (A+B)")

# AB edges within coset 0
ab_edges_0 = [e for e in coset0_edges if e['vtype'] == ('A','B')]
bb_edges_0 = [e for e in coset0_edges if e['vtype'] == ('B','B')]
print_sector_stats(ab_edges_0, "A-B within 2T")
print_sector_stats(bb_edges_0, "B-B within 2T")

# C cosets
for cid in range(1, 5):
    subset = [e for e in edge_data if e['coset_pair'] == (cid,cid)]
    print_sector_stats(subset, f"C-coset {cid} (intra)")

# ===================================================================
# PART I: THE KEY TEST - PROJECTOR WEIGHTS ON VERTICES
# ===================================================================
print(f"\n--- PART I: Projector diagonal by vertex type ---")

print("  P_sector(v,v) = probability that vertex v is in sector")
print()

for sector_name, P in [("RAT", P_rat), ("PHI", P_phi), ("PHI'", P_phip)]:
    diag = np.diag(P)
    # By vertex type
    a_avg = np.mean(diag[type_A])
    b_avg = np.mean(diag[type_B])
    c_avg = np.mean(diag[type_C])
    print(f"  P_{sector_name:>4s}(v,v): A={a_avg:.4f}  B={b_avg:.4f}  C={c_avg:.4f}")

print(f"\n  (Should all be: RAT={94/120:.4f}, PHI={13/120:.4f}, PHI'={13/120:.4f})")
print(f"  (Because vertex-transitive => uniform diagonal)")

# ===================================================================
# PART J: SECTOR WEIGHT DISTRIBUTION FOR EDGES
# ===================================================================
print(f"\n--- PART J: Distribution of phi/phi' weights ---")

# For each edge, compute the asymmetry: (w_phi - w_phip)
# If Galois complementarity is real, different edge types should
# have different asymmetries.

for label, subset in [
    ("AB (intra-2T)", [e for e in edge_data if e['vtype']==('A','B') and e['same_coset']]),
    ("BB (intra-2T)", [e for e in edge_data if e['vtype']==('B','B') and e['same_coset']]),
    ("AC (inter)", [e for e in edge_data if e['vtype']==('A','C')]),
    ("BC (inter)", [e for e in edge_data if e['vtype']==('B','C')]),
    ("CC (intra-C)", [e for e in edge_data if e['vtype']==('C','C') and e['same_coset']]),
    ("CC (inter-C)", [e for e in edge_data if e['vtype']==('C','C') and not e['same_coset']]),
]:
    if not subset:
        print(f"  {label:>20s}: no edges")
        continue
    asym = [e['w_phi'] - e['w_phip'] for e in subset]
    phi_vals = [e['w_phi'] for e in subset]
    phip_vals = [e['w_phip'] for e in subset]
    rat_vals = [e['w_rat'] for e in subset]

    n_unique_asym = len(set(round(a, 6) for a in asym))
    print(f"  {label:>20s} ({len(subset):>3d}): "
          f"phi={np.mean(phi_vals):>7.4f}  phi'={np.mean(phip_vals):>7.4f}  "
          f"asym={np.mean(asym):>8.5f}  rat={np.mean(rat_vals):>7.4f}  "
          f"({n_unique_asym} distinct asym values)")

# ===================================================================
# PART K: GAUGE CHANNEL IDENTIFICATION
# ===================================================================
print(f"\n--- PART K: Gauge channel attempt ---")

# In the framework:
# U(1): 48 edges forming a perfect matching on CC edges
# SU(2): 288 edges
# SU(3): 384 edges

# The shared-neighbor count might help classify:
print("  Shared neighbor profile for each edge type:")

for label, subset in [
    ("AB", [e for e in edge_data if e['vtype']==('A','B')]),
    ("AC", [e for e in edge_data if e['vtype']==('A','C')]),
    ("BC", [e for e in edge_data if e['vtype']==('B','C')]),
    ("CC intra", [e for e in edge_data if e['vtype']==('C','C') and e['same_coset']]),
    ("CC inter", [e for e in edge_data if e['vtype']==('C','C') and not e['same_coset']]),
    ("BB", [e for e in edge_data if e['vtype']==('B','B')]),
]:
    if not subset:
        continue
    shared_counts = []
    for e in subset[:20]:  # sample
        i, j = e['i'], e['j']
        shared = sum(1 for k in range(120) if adj[i,k] and adj[j,k])
        shared_counts.append(shared)
    counts = Counter(shared_counts)
    print(f"    {label:>12s}: shared neighbors = {dict(counts)}")

# ===================================================================
# PART L: GALOIS ASYMMETRY BY SHARED NEIGHBORS
# ===================================================================
print(f"\n--- PART L: Galois asymmetry vs shared neighbors ---")

# Classify ALL edges by shared neighbor count
shared_neighbor_data = {}
for e in edge_data:
    i, j = e['i'], e['j']
    shared = sum(1 for k in range(120) if adj[i,k] and adj[j,k])
    if shared not in shared_neighbor_data:
        shared_neighbor_data[shared] = []
    shared_neighbor_data[shared].append(e)

print(f"  {'Shared':>7s} {'Count':>6s} {'phi':>8s} {'phi_p':>8s} {'asym':>9s} {'rat':>8s}")
print(f"  {'-'*50}")
for shared in sorted(shared_neighbor_data.keys()):
    edges = shared_neighbor_data[shared]
    w_phi = np.mean([e['w_phi'] for e in edges])
    w_phip = np.mean([e['w_phip'] for e in edges])
    w_rat = np.mean([e['w_rat'] for e in edges])
    asym = w_phi - w_phip
    print(f"  {shared:>7d} {len(edges):>6d} {w_phi:>8.5f} {w_phip:>8.5f} {asym:>9.5f} {w_rat:>8.5f}")

# ===================================================================
# PART M: THE ADJACENCY MATRIX IN GALOIS SECTORS
# ===================================================================
print(f"\n--- PART M: Sector adjacency matrices ---")

# Trace of A_sector over edges
print("  Sector contributions to total adjacency:")
for name, A_sec in [("A_rat", A_rat), ("A_phi", A_phi), ("A_phip", A_phip)]:
    # Sum over all edges
    total = 0
    for i in range(120):
        for j in range(i+1, 120):
            if adj[i,j]:
                total += A_sec[i,j]
    print(f"    sum_edges {name}(i,j) = {total:.4f}")

# Frobenius norm of off-diagonal parts
print(f"\n  Off-diagonal norm of sector matrices:")
for name, A_sec in [("A_rat", A_rat), ("A_phi", A_phi), ("A_phip", A_phip)]:
    off_diag = A_sec - np.diag(np.diag(A_sec))
    frob = np.sqrt(np.sum(off_diag**2))
    print(f"    ||{name}||_F = {frob:.4f}")

# ===================================================================
# PART N: PROJECTOR ON EDGE SUBSETS
# ===================================================================
print(f"\n--- PART N: P_sector projectors on edge subsets ---")

# For each edge type, compute sum of P_sector(i,j)
# This tells us how much each edge type "lives" in each sector

print(f"  {'Edge type':>20s} {'N':>5s} {'P_rat':>8s} {'P_phi':>8s} {'P_phip':>8s} "
      f"{'phi-phip':>9s} {'ratio':>7s}")
print(f"  {'-'*65}")

all_results = []
for label, subset in [
    ("AB intra-2T", [e for e in edge_data if e['vtype']==('A','B') and e['same_coset']]),
    ("BB intra-2T", [e for e in edge_data if e['vtype']==('B','B') and e['same_coset']]),
    ("AC inter", [e for e in edge_data if e['vtype']==('A','C')]),
    ("BC inter", [e for e in edge_data if e['vtype']==('B','C')]),
    ("CC intra-C", [e for e in edge_data if e['vtype']==('C','C') and e['same_coset']]),
    ("CC inter-C", [e for e in edge_data if e['vtype']==('C','C') and not e['same_coset']]),
    ("ALL intra", intra_edges),
    ("ALL inter", inter_edges),
    ("ALL edges", edge_data),
]:
    if not subset:
        continue
    p_rat = np.mean([P_rat[e['i'],e['j']] for e in subset])
    p_phi = np.mean([P_phi[e['i'],e['j']] for e in subset])
    p_phip = np.mean([P_phip[e['i'],e['j']] for e in subset])
    diff = p_phi - p_phip
    ratio = p_phi / p_phip if abs(p_phip) > 1e-10 else float('inf')
    print(f"  {label:>20s} {len(subset):>5d} {p_rat:>8.5f} {p_phi:>8.5f} {p_phip:>8.5f} "
          f"{diff:>9.5f} {ratio:>7.3f}")
    all_results.append((label, len(subset), p_rat, p_phi, p_phip))

# ===================================================================
# PART O: SUMMARY AND INTERPRETATION
# ===================================================================
print(f"\n{'='*70}")
print("SUMMARY AND INTERPRETATION")
print("="*70)

print(f"""
  KEY RESULT: The spectral decomposition A = A_rat + A_phi + A_phip
  gives the Galois sector weight of each edge.

  For the 600-cell, which is VERTEX-TRANSITIVE:
  - The diagonal P_sector(v,v) is uniform (each vertex sees all sectors equally)
  - The EDGE weights carry the sector information

  If the weights are uniform across all edge types:
  - This means vertex-transitivity implies EDGE-uniformity too
  - Then Galois complementarity cannot be explained at the spectral level alone
  - It would need to come from the DYNAMICS (soliton/Potts model), not the spectrum

  If the weights differ by edge type:
  - Different gauge channels genuinely project onto different Galois sectors
  - This would be a SPECTRAL DERIVATION of Galois complementarity
""")

# Check uniformity
all_phis = [e['w_phi'] for e in edge_data]
all_phips = [e['w_phip'] for e in edge_data]

phi_std = np.std(all_phis)
phip_std = np.std(all_phips)
print(f"  Edge weight variation:")
print(f"  w_phi:  mean={np.mean(all_phis):.5f}, std={phi_std:.5f}")
print(f"  w_phip: mean={np.mean(all_phips):.5f}, std={phip_std:.5f}")

if phi_std < 0.001:
    print(f"\n  RESULT: Weights are UNIFORM across edges (std < 0.001)")
    print(f"  => Vertex-transitivity forces edge-uniformity")
    print(f"  => Galois complementarity is NOT spectral, it's DYNAMIC")
    print(f"  => The complementarity comes from soliton/gauge DYNAMICS, not from")
    print(f"     the static spectrum of the adjacency matrix.")
    print(f"\n  This is actually EXPECTED: vertex-transitivity means the graph")
    print(f"  'looks the same from every vertex', so edges cannot distinguish")
    print(f"  between Galois sectors at the spectral level.")
    print(f"\n  The complementarity must arise from:")
    print(f"  - The BREAKING of vertex-transitivity by soliton placement")
    print(f"  - The specific GAUGE CHANNEL (edge subset) each particle couples to")
    print(f"  - The coset structure (which breaks the symmetry)")
else:
    print(f"\n  RESULT: Weights VARY across edges (std = {phi_std:.5f})")
    print(f"  => Edge types DO project onto different Galois sectors!")
