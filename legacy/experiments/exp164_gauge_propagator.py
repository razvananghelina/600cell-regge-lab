"""
EXP-164: Gauge-Decomposed Propagator & Interaction Amplitudes
=============================================================
How do fermions interact through the gauge sector?

Strategy:
1. Build full 600-cell graph with vertex types (A, B, C) and edge types (U1, SU2, SU3)
2. Compute Green's function restricted to different channels
3. Compute effective interaction: G_gauge(C_i, C_j) = sum over paths through gauge vertices
4. Compare U(1), SU(2), SU(3) mediated amplitudes
5. Check if coupling ratios emerge from graph structure

Edge types (from exp143):
- U(1): 96 edges (1 per fermion, CC-type)
- SU(2)_charged: 192 edges (2 per fermion)
- SU(2)_neutral: 48 edges (0.5 per fermion, BB-type)
- SU(3): 384 edges (4 per fermion, CC-type)
Total: 720

Selection rules (from exp136):
- A-A = A-B = B-B = 0 (no gauge-gauge direct edges)
- Only C-C, A-C, B-C edges exist
"""

import numpy as np
from itertools import product, permutations
from collections import Counter, defaultdict
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix, lil_matrix

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-164: Gauge-Decomposed Propagator")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell
# ============================================================
print("\n--- Step 1: Build 600-cell ---")

verts_set = set()
for i in range(4):
    for s in [1, -1]:
        v = [0,0,0,0]; v[i] = s
        verts_set.add(tuple(v))
for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = [p for p in permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in product([1,-1], repeat=len(nz)):
        v = list(base)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

verts = np.array(sorted(verts_set))
N = len(verts)
print(f"  Vertices: {N}")

# Adjacency
dots = verts @ verts.T
adj = (np.abs(dots - PHI/2) < 0.01).astype(int)
np.fill_diagonal(adj, 0)

# ============================================================
# Step 2: Classify vertices
# ============================================================
print("\n--- Step 2: Vertex classification ---")

type_A, type_B, type_C = [], [], []
for i in range(N):
    c = verts[i]
    nz = np.sum(np.abs(c) > 0.01)
    if nz == 1 and np.isclose(np.max(np.abs(c)), 1.0):
        type_A.append(i)
    elif nz == 4 and np.allclose(np.abs(c), 0.5, atol=0.01):
        type_B.append(i)
    else:
        type_C.append(i)

print(f"  A: {len(type_A)}, B: {len(type_B)}, C: {len(type_C)}")
vertex_type = {}
for i in type_A: vertex_type[i] = 'A'
for i in type_B: vertex_type[i] = 'B'
for i in type_C: vertex_type[i] = 'C'

# ============================================================
# Step 3: Classify edges
# ============================================================
print("\n--- Step 3: Edge classification ---")

edges = []
edge_types = defaultdict(list)
for i in range(N):
    for j in range(i+1, N):
        if adj[i, j]:
            ti, tj = vertex_type[i], vertex_type[j]
            etype = ''.join(sorted([ti, tj]))
            edges.append((i, j, etype))
            edge_types[etype].append((i, j))

for etype, elist in sorted(edge_types.items()):
    print(f"  {etype}: {len(elist)} edges")

# Sub-classify CC edges into U(1) and SU(3) using shared neighbor structure
# From exp143: U(1) = "special CC" edge (degree 4, specific shared-neighbor pattern)
# SU(3) = remaining CC edges

# For each CC edge, count shared A-neighbors
cc_u1 = []
cc_su3 = []
for i, j in edge_types['CC']:
    # Shared neighbors
    shared = np.where((adj[i] > 0) & (adj[j] > 0))[0]
    shared_A = sum(1 for s in shared if vertex_type[s] == 'A')
    shared_B = sum(1 for s in shared if vertex_type[s] == 'B')
    shared_C = sum(1 for s in shared if vertex_type[s] == 'C')

    # U(1) edge: unique per fermion, has specific pattern
    # From exp143: special CC edge has shared_A=1, shared_B=0
    if shared_A == 1 and shared_B == 0:
        cc_u1.append((i, j))
    else:
        cc_su3.append((i, j))

print(f"\n  CC sub-classification:")
print(f"    U(1) (shared_A=1, shared_B=0): {len(cc_u1)}")
print(f"    SU(3) (other CC): {len(cc_su3)}")
print(f"    AC edges (SU(2/3)-related): {len(edge_types['AC'])}")
print(f"    BC edges (SU(2)-related): {len(edge_types['BC'])}")

# Full gauge classification
gauge_edges = {
    'U1': cc_u1,
    'SU2_AC': edge_types['AC'],
    'SU2_BC': edge_types.get('BC', []) + edge_types.get('BB', []),
    'SU3': cc_su3
}

for gtype, elist in gauge_edges.items():
    print(f"  {gtype}: {len(elist)} edges")

# ============================================================
# Step 4: Full Green's function
# ============================================================
print("\n--- Step 4: Full Green's function ---")

L = np.diag(adj.sum(axis=1).astype(float)) - adj.astype(float)
mu2 = 0.1  # mass regularization
G_full = np.linalg.inv(L + mu2 * np.eye(N))

# Fermion-fermion propagator statistics
CC_props = []
for i in type_C:
    for j in type_C:
        if i < j:
            CC_props.append(G_full[i, j])

CC_props = np.array(CC_props)
print(f"  G(C,C) stats: mean={np.mean(CC_props):.6f}, std={np.std(CC_props):.6f}")
print(f"  G(C,C) range: [{np.min(CC_props):.6f}, {np.max(CC_props):.6f}]")

# ============================================================
# Step 5: Channel-restricted propagators
# ============================================================
print("\n--- Step 5: Channel-restricted propagators ---")

# Build adjacency matrices for each gauge channel
def build_channel_adj(edge_list, n):
    """Build adjacency matrix from edge list."""
    A = np.zeros((n, n))
    for i, j in edge_list:
        A[i, j] = 1
        A[j, i] = 1
    return A

# Method: For each channel, compute the Green's function on the subgraph
# that only includes those edges

for channel_name, channel_edges in gauge_edges.items():
    A_ch = build_channel_adj(channel_edges, N)
    deg_ch = A_ch.sum(axis=1)
    n_connected = np.sum(deg_ch > 0)

    # Green's function on this channel
    L_ch = np.diag(deg_ch) - A_ch
    G_ch = np.linalg.inv(L_ch + mu2 * np.eye(N))

    # Fermion-fermion amplitude through this channel
    cc_amps = []
    for i in type_C[:20]:  # Sample
        for j in type_C[:20]:
            if i < j:
                cc_amps.append(G_ch[i, j])

    cc_amps = np.array(cc_amps)
    print(f"\n  Channel {channel_name}:")
    print(f"    Edges: {len(channel_edges)}, Connected vertices: {n_connected}")
    print(f"    G(C,C) mean: {np.mean(cc_amps):.6f}, std: {np.std(cc_amps):.6f}")
    print(f"    G(C,C) max |amplitude|: {np.max(np.abs(cc_amps)):.6f}")

# ============================================================
# Step 6: Effective coupling from propagator ratios
# ============================================================
print("\n--- Step 6: Effective coupling ratios ---")

# For each pair of fermions, compute the ratio of amplitudes through different channels
# This gives the RELATIVE coupling strength

# Build all channel Green's functions
G_channels = {}
for channel_name, channel_edges in gauge_edges.items():
    A_ch = build_channel_adj(channel_edges, N)
    L_ch = np.diag(A_ch.sum(axis=1)) - A_ch
    G_channels[channel_name] = np.linalg.inv(L_ch + mu2 * np.eye(N))

# Sample fermion pairs at different graph distances
from collections import deque

def bfs_dist(adj_matrix, source):
    n = len(adj_matrix)
    dist = [-1]*n; dist[source] = 0
    q = deque([source])
    while q:
        u = q.popleft()
        for v in range(n):
            if adj_matrix[u,v] and dist[v]==-1:
                dist[v] = dist[u]+1; q.append(v)
    return dist

# Pick a reference fermion
ref = type_C[0]
dists = bfs_dist(adj, ref)

print(f"\n  Reference fermion: vertex {ref}")
print(f"  Channel amplitudes by distance:")

for d in range(1, 6):
    targets = [j for j in type_C if dists[j] == d]
    if not targets:
        continue

    amps_by_channel = {}
    for ch_name, G_ch in G_channels.items():
        ch_amps = [G_ch[ref, j] for j in targets[:10]]
        amps_by_channel[ch_name] = np.mean(ch_amps)

    full_amp = np.mean([G_full[ref, j] for j in targets[:10]])

    print(f"\n  d={d} ({len(targets)} fermions):")
    print(f"    Full: {full_amp:.6f}")
    for ch_name, amp in amps_by_channel.items():
        ratio = amp / full_amp if abs(full_amp) > 1e-10 else 0
        print(f"    {ch_name}: {amp:.6f} (ratio={ratio:.3f})")

# ============================================================
# Step 7: Vertex amplitude (3-point function)
# ============================================================
print("\n--- Step 7: 3-point vertex amplitudes ---")

# For three fermions sharing a common gauge vertex (A or B),
# compute the effective 3-point coupling
# This is: V(i,j,k) = sum_{g in gauge} adj(g,i)*adj(g,j)*adj(g,k)

# For each gauge vertex, find its fermion neighbors
gauge_fermion_triples = defaultdict(list)
for g in type_A + type_B:
    g_type = vertex_type[g]
    fermion_nbrs = [v for v in range(N) if adj[g, v] and vertex_type[v] == 'C']

    # All triples of fermion neighbors
    n_triples = 0
    for a in range(len(fermion_nbrs)):
        for b in range(a+1, len(fermion_nbrs)):
            for c in range(b+1, len(fermion_nbrs)):
                n_triples += 1
    gauge_fermion_triples[g_type].append((g, len(fermion_nbrs), n_triples))

for g_type in ['A', 'B']:
    data = gauge_fermion_triples[g_type]
    nbr_counts = [d[1] for d in data]
    triple_counts = [d[2] for d in data]
    print(f"\n  {g_type}-vertex gauge bosons:")
    print(f"    Fermion neighbors: {Counter(nbr_counts)}")
    print(f"    3-fermion triples per gauge vertex: {Counter(triple_counts)}")

# Total 3-point vertices
total_A_triples = sum(d[2] for d in gauge_fermion_triples['A'])
total_B_triples = sum(d[2] for d in gauge_fermion_triples['B'])
print(f"\n  Total 3-point vertices: A-mediated={total_A_triples}, B-mediated={total_B_triples}")
print(f"  Ratio B/A: {total_B_triples/total_A_triples:.3f}")

# ============================================================
# Step 8: 4-point function (box diagram)
# ============================================================
print("\n--- Step 8: 4-point amplitudes ---")

# For 4 fermions forming a "box" through gauge exchange:
# i -- g1 -- j
# |          |
# g4         g2
# |          |
# l -- g3 -- k

# Count: how many 4-cycles go through exactly 2 gauge and 2 fermion vertices?
n_box = 0
n_box_types = Counter()

for g1 in type_A + type_B:
    g1_fnbrs = set(v for v in range(N) if adj[g1, v] and vertex_type[v] == 'C')
    for g2 in type_A + type_B:
        if g2 <= g1:
            continue
        g2_fnbrs = set(v for v in range(N) if adj[g2, v] and vertex_type[v] == 'C')

        # Shared fermion neighbors
        shared = g1_fnbrs & g2_fnbrs
        n_shared = len(shared)
        if n_shared >= 2:
            # Each pair of shared fermions gives a "box"
            box_count = n_shared * (n_shared - 1) // 2
            n_box += box_count
            key = (vertex_type[g1], vertex_type[g2])
            n_box_types[tuple(sorted(key))] += box_count

print(f"  4-point box diagrams (gauge-fermion-gauge-fermion):")
print(f"  Total: {n_box}")
for key, count in sorted(n_box_types.items()):
    print(f"    {key[0]}-{key[1]} mediated: {count}")

# ============================================================
# Step 9: Scattering matrix between fermion pairs
# ============================================================
print("\n--- Step 9: Scattering matrix ---")

# For nearest-neighbor fermion pairs, compute the "scattering amplitude"
# T(i->j) = G_full(i,j) - G_free(i,j) where G_free is diagonal

# More interesting: for each fermion, decompose its coupling to neighbors
# by gauge channel

for c_idx in type_C[:3]:
    c_nbrs = [v for v in range(N) if adj[c_idx, v]]
    c_nbrs_by_type = defaultdict(list)
    for n in c_nbrs:
        c_nbrs_by_type[vertex_type[n]].append(n)

    print(f"\n  Fermion {c_idx}: {len(c_nbrs)} neighbors")
    for vt, vlist in sorted(c_nbrs_by_type.items()):
        # For each neighbor, the "coupling strength" through each channel
        for v in vlist[:2]:
            amps = {}
            for ch_name, G_ch in G_channels.items():
                amps[ch_name] = G_ch[c_idx, v]
            full = G_full[c_idx, v]
            print(f"    -> {vt}-{v}: full={full:.5f}, " +
                  ", ".join(f"{k}={v:.5f}" for k, v in amps.items()))

# ============================================================
# Step 10: Confinement signal
# ============================================================
print("\n--- Step 10: Confinement analysis ---")

# The propagator changes sign at d=3 (from exp137)
# This means color-singlet states have different propagation than colored states

# Compute: for each channel, at what distance does the propagator change sign?
for ch_name, G_ch in G_channels.items():
    ref = type_C[0]
    sign_changes = []
    for d in range(1, 6):
        targets_d = [j for j in range(N) if dists[j] == d]
        if targets_d:
            mean_amp = np.mean([G_ch[ref, j] for j in targets_d[:20]])
            sign_changes.append((d, mean_amp))

    print(f"\n  {ch_name} propagator by distance:")
    for d, amp in sign_changes:
        sign = '+' if amp > 0 else '-'
        print(f"    d={d}: {amp:.6f} ({sign})")

# Full propagator
print(f"\n  Full propagator by distance:")
for d in range(1, 6):
    targets_d = [j for j in range(N) if dists[j] == d]
    if targets_d:
        mean_amp = np.mean([G_full[ref, j] for j in targets_d[:20]])
        sign = '+' if mean_amp > 0 else '-'
        print(f"    d={d}: {mean_amp:.6f} ({sign})")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key questions answered:
1. How do different gauge channels mediate fermion interactions?
2. What are the effective coupling ratios from graph structure?
3. How many 3-point and 4-point vertices exist?
4. At what distance does each channel show confinement?

STATUS: Check output above.
""")
