import numpy as np
from scipy import linalg
from itertools import permutations

# Constante
PHI = (1 + 5**0.5) / 2

def get_600_cell():
    vertices = []
    vertex_types = []
    vertex_subgroup = []

    # 1. Type A (16-cell) - Gauge
    for i in range(4):
        for s in [-1, 1]:
            v = [0, 0, 0, 0]
            v[i] = s
            vertices.append(v)
            vertex_types.append('A')
            vertex_subgroup.append(-1)

    # 2. Type B (tesseract) - Gauge
    for s in [[(i>>j)&1 for j in range(4)] for i in range(16)]:
        v = [0.5 if x==0 else -0.5 for x in s]
        vertices.append(v)
        vertex_types.append('B')
        vertex_subgroup.append(-1)

    # 3. Type C (fermions) - 96 vertices
    base = [PHI/2, 0.5, 1/(2*PHI), 0]
    perms = [p for p in permutations([0,1,2,3])]
    
    def get_parity(p):
        p = list(p)
        res = 0
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if p[i] > p[j]: res += 1
        return res % 2

    even_perms = [p for p in perms if get_parity(p) == 0]
    
    seen = set()
    for p in even_perms:
        coords = [base[i] for i in p]
        for signs in range(8):
            s = [1 if (signs >> i) & 1 else -1 for i in range(3)]
            v = [0, 0, 0, 0]
            si = 0
            zero_pos = -1
            for i in range(4):
                if coords[i] == 0:
                    v[i] = 0
                    zero_pos = i
                else:
                    v[i] = s[si] * coords[i]
                    si += 1
            
            v_tuple = tuple(np.round(v, 8))
            if v_tuple not in seen:
                seen.add(v_tuple)
                vertices.append(list(v))
                vertex_types.append('C')
                vertex_subgroup.append(zero_pos)

    vertices = np.array(vertices)
    N = len(vertices)
    
    # Adjacency
    dist = np.zeros((N, N))
    for i in range(N):
        dist[i] = np.linalg.norm(vertices - vertices[i], axis=1)
    
    edge = 1/PHI
    adj = (np.abs(dist - edge) < 1e-3).astype(float)
    return vertices, vertex_types, vertex_subgroup, adj

print("Generating 600-cell...")
verts, v_types, v_groups, adj = get_600_cell()
N = len(verts)

# Laplacian
deg = np.sum(adj, axis=1)
L = np.diag(deg) - adj

# Green's function (Moore-Penrose Pseudoinverse)
print("Computing Green's function (Pseudoinverse)...")
G = np.linalg.pinv(L)

# Grupare indici
idx_A = [i for i, t in enumerate(v_types) if t == 'A']
idx_B = [i for i, t in enumerate(v_types) if t == 'B']
idx_C = [i for i, t in enumerate(v_types) if t == 'C']
groups = {0:[], 1:[], 2:[], 3:[]}
for i in idx_C:
    groups[v_groups[i]].append(i)

# Presupunem axa 3 ca fiind "timpul"
# Atunci G3 este sectorul "sterile" (deoarece x3=0 în G3)
print("\nAssuming x3 as temporal axis:")
print("Visible Generations: G0, G1, G2")
print("Sterile Generation: G3 (x3=0)")

# Analiza propagatorului din G0
src = groups[0][0]
print(f"\nSource vertex: Index {src} (Type C, Group G0)")

avg_G = {}
for g in range(4):
    vals = [G[src, j] for j in groups[g]]
    avg_G[f"G{g}"] = np.mean(vals)
    
avg_G["Gauge_A"] = np.mean([G[src, j] for j in idx_A])
avg_G["Gauge_B"] = np.mean([G[src, j] for j in idx_B])

print("\nAverage Propagator Magnitude:")
for k, v in avg_G.items():
    ratio = v / avg_G["G0"]
    print(f"  To {k:<8}: {v:>10.6f} (Ratio: {ratio:.4f})")

# Distante grafice
def get_dists(adj, start):
    dists = [-1] * len(adj)
    dists[start] = 0
    q = [start]
    while q:
        u = q.pop(0)
        for v in np.where(adj[u] > 0)[0]:
            if dists[v] == -1:
                dists[v] = dists[u] + 1
                q.append(v)
    return dists

dists_from_src = get_dists(adj, src)
print("\nPropagator by Graph Distance:")
dist_vals = {}
for i in range(N):
    d = dists_from_src[i]
    if d not in dist_vals: dist_vals[d] = []
    dist_vals[d].append(G[src, i])

for d in sorted(dist_vals.keys()):
    avg = np.mean(dist_vals[d])
    print(f"  Distance {d}: {avg:>10.6f} (Count: {len(dist_vals[d])})")

leakage_sterile = avg_G["G3"]
print(f"\nSterile leakage (G0 -> G3): {leakage_sterile:.6f}")

eigs = np.linalg.eigvalsh(L)
print(f"\nSpectral Gap (lambda_1): {eigs[1]:.6f}")
print(f"Largest Eigenvalue: {eigs[-1]:.6f}")

print("\n600-cell is highly symmetric. Mass ratio 10^-7 is NOT directly here.")