
import numpy as np
from itertools import permutations
from collections import deque

PHI = (1 + 5**0.5) / 2

def build_600cell():
    phi = PHI
    vertices = []
    # A, B, C combined
    for i in range(4):
        for s in [1, -1]:
            vertices.append(tuple(np.round([s if j==i else 0 for j in range(4)], 10)))
    for s in [[(i>>j)&1 for j in range(4)] for i in range(16)]:
        vertices.append(tuple(np.round([0.5 if x==0 else -0.5 for x in s], 10)))
    base_vals = [phi/2, 0.5, 1/(2*phi), 0]
    for s in [[(i>>j)&1 for j in range(3)] for i in range(8)]:
        sign = [1 if x==0 else -1 for x in s]
        for p in permutations(range(4)):
            inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
            if inv % 2 == 0:
                v = [0]*4
                si = 0
                for i in range(4):
                    if base_vals[p[i]] == 0: v[i] = 0
                    else:
                        v[i] = sign[si] * base_vals[p[i]]
                        si += 1
                vertices.append(tuple(np.round(v, 10)))
    unique_v = list(set(vertices))
    verts = [np.array(v) for v in unique_v]
    adj = np.zeros((len(verts), len(verts)), dtype=int)
    for i in range(len(verts)):
        for j in range(i+1, len(verts)):
            if np.linalg.norm(verts[i] - verts[j]) < 1/PHI + 0.01:
                adj[i,j] = adj[j,i] = 1
    return verts, adj

def build_120cell(verts_600, adj_600):
    neighbors = [set(np.where(adj_600[i] > 0)[0]) for i in range(len(verts_600))]
    tetrahedra = []
    for i in range(len(verts_600)):
        for j in neighbors[i]:
            if j <= i: continue
            common_ij = neighbors[i] & neighbors[j]
            for k in common_ij:
                if k <= j: continue
                common_ijk = common_ij & neighbors[k]
                for l in common_ijk:
                    if l <= k: continue
                    tetrahedra.append((i, j, k, l))
    
    face_to_tet = {}
    for idx, tet in enumerate(tetrahedra):
        for skip in range(4):
            face = frozenset(tet[m] for m in range(4) if m != skip)
            face_to_tet.setdefault(face, []).append(idx)
    
    adj_120 = np.zeros((600, 600), dtype=int)
    for tets in face_to_tet.values():
        if len(tets) == 2:
            i, j = tets
            adj_120[i, j] = adj_120[j, i] = 1
    return adj_120

print("Building 600-cell...")
v6, a6 = build_600cell()
print("Building 120-cell...")
a120 = build_120cell(v6, a6)

# Laplacian 120-cell
L = np.diag(np.sum(a120, axis=1)) - a120
print("Computing 120-cell Green's function...")
G = np.linalg.pinv(L)

# Analyse suppression by distance
def bfs_dists(adj, start):
    dists = [-1] * 600
    dists[start] = 0
    q = deque([start])
    while q:
        u = q.popleft()
        for v in np.where(adj[u] > 0)[0]:
            if dists[v] == -1:
                dists[v] = dists[u] + 1
                q.append(v)
    return dists

dists = bfs_dists(a120, 0)
print("\n120-Cell Propagator by Distance (Source vertex 0):")
print(f"{'Dist':>5} {'Avg G':>15} {'Ratio to self':>15}")
self_G = G[0, 0]
for d_target in range(16):
    indices = [i for i, d in enumerate(dists) if d == d_target]
    if indices:
        avg = np.mean([G[0, i] for i in indices])
        print(f"{d_target:5d} {avg:15.8f} {avg/self_G:15.8f}")

# Max suppression
max_suppression = np.min(np.abs(G[0])) / self_G
print(f"\nMax relative suppression: {max_suppression:.8f}")

# Compare with m_nu / m_e ~ 10^-7
print(f"Target suppression (m_nu/m_e): ~1e-7")
if max_suppression < 1e-4:
    print("Promising suppression found!")
else:
    print("Suppression is not enough for neutrinos directly.")
