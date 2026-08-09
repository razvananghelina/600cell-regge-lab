import numpy as np
import scipy.linalg as la

# Constants
PHI = (1 + np.sqrt(5)) / 2

def generate_600_cell_vertices():
    # Same code as before
    vertices = []
    for i in range(16):
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
        vertices.append(signs * 0.5)
    for i in range(4):
        v = np.zeros(4); v[i] = 1.0
        vertices.append(v); vertices.append(-v)
    v_base = np.array([PHI, 1.0, 1.0/PHI, 0.0]) * 0.5
    from itertools import permutations
    perms = [p for p in list(permutations([0, 1, 2, 3])) 
             if sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j]) % 2 == 0]
    for p in perms:
        perm_vals = v_base[list(p)]
        for i in range(16):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
            vertices.append(perm_vals * signs)
    vertices = np.array(vertices)
    _, idx = np.unique(np.round(vertices, 8), axis=0, return_index=True)
    return vertices[idx]

def get_120_cell_vertices(verts_600):
    """
    The 120-cell is the dual of the 600-cell.
    Its vertices are the centers of the 600-cell's cells (tetrahedra).
    
    Algorithm:
    1. Find all tetrahedra in 600-cell.
       A tetrahedron is formed by 4 vertices that are all connected to each other.
       This is computationally expensive (Clique problem).
       
    Optimization:
    The 600-cell is a regular polytope. We can generate 120-cell vertices directly using coordinates.
    The 120-cell vertices are permutations of:
    - (0, 0, 2, 2)
    - (1, 1, 1, sqrt(5))
    - (phi^-2, phi, phi, phi)
    - (phi^-1, phi^-1, phi^-1, phi^2)
    ... scaled appropriately.
    
    BUT, easier method:
    The 600 vertices of the 120-cell lie on shells.
    Let's trust the "Center of Tetrahedra" method but implement it efficiently.
    Each vertex in 600-cell has 20 neighbors. These form the faces.
    Wait, finding cliques in a 120-node graph is fast enough.
    """
    
    print("Generating 120-cell vertices via Dual construction...")
    
    # 1. Build Adjacency
    dists = la.norm(verts_600[:, None, :] - verts_600[None, :, :], axis=2)
    flat_dists = dists[dists > 0.0001]
    edge_len = np.min(flat_dists)
    adj = (np.abs(dists - edge_len) < 0.001)
    
    # 2. Find Tetrahedra (Cliques of size 4)
    # Since we know the geometry, we can iterate:
    # For each vertex v1, iterate neighbors v2.
    #   Iterate common neighbors of v1, v2 -> v3.
    #     Iterate common neighbors of v1, v2, v3 -> v4.
    
    tetrahedra_centers = []
    
    # This might produce duplicates, we'll filter later.
    # To speed up, we only search for indices i < j < k < l
    
    num_verts = len(verts_600)
    
    # Pre-calculate neighbor lists
    neighbors = [set(np.where(adj[i])[0]) for i in range(num_verts)]
    
    found_tets = set()
    
    for i in range(num_verts):
        for j in neighbors[i]:
            if j <= i: continue
            
            common_ij = neighbors[i].intersection(neighbors[j])
            
            for k in common_ij:
                if k <= j: continue
                
                common_ijk = common_ij.intersection(neighbors[k])
                
                for l in common_ijk:
                    if l <= k: continue
                    
                    # Found a tetrahedron {i, j, k, l}
                    tet_indices = tuple(sorted((i, j, k, l)))
                    if tet_indices not in found_tets:
                        found_tets.add(tet_indices)
                        
                        # Calculate Center
                        pts = verts_600[list(tet_indices)]
                        center = np.mean(pts, axis=0)
                        
                        # Normalize to put on hypersphere (optional but good for regularity)
                        # center = center / np.linalg.norm(center) 
                        # Actually 120-cell vertices should be at specific radius.
                        # We'll just take raw centers first.
                        tetrahedra_centers.append(center)
                        
    print(f"Found {len(found_tets)} tetrahedra (Expected 600).")
    return np.array(tetrahedra_centers)

def analyze_strong_force():
    verts_600 = generate_600_cell_vertices()
    verts_120 = get_120_cell_vertices(verts_600)
    
    if len(verts_120) != 600:
        print(f"WARNING: Expected 600 vertices for 120-cell, found {len(verts_120)}")
    
    # Normalize vertices to radius 1 for consistency
    norms = np.linalg.norm(verts_120, axis=1)
    avg_norm = np.mean(norms)
    verts_120 = verts_120 / avg_norm
    
    print("\n--- Analyzing 120-cell (Dual) Graph ---")
    
    # Build Laplacian
    dists = la.norm(verts_120[:, None, :] - verts_120[None, :, :], axis=2)
    flat = dists[dists > 0.0001]
    if len(flat) == 0:
        print("Error: No edges found.")
        return
        
    edge_len = np.min(flat)
    print(f"Edge length (Dual): {edge_len}")
    
    adj = (np.abs(dists - edge_len) < 0.001).astype(float)
    deg = np.diag(np.sum(adj, axis=0))
    print(f"Average degree (Should be 4): {np.mean(np.diag(deg))}")
    
    lap = deg - adj
    d_inv_sqrt = np.linalg.inv(np.sqrt(deg))
    norm_lap = np.dot(d_inv_sqrt, np.dot(lap, d_inv_sqrt))
    
    evals = np.linalg.eigvalsh(norm_lap)
    
    # Analyze Spectrum
    lambda_1 = evals[1]
    print(f"Spectral Gap (Dual): {lambda_1}")
    
    # Calculate Trace(G^2) with mass = lambda_1
    # G = (L + m^2)^-1
    tr_g2 = np.sum(1.0 / (evals + lambda_1)**2)
    print(f"Trace(G^2): {tr_g2}")
    
    # Compare with Strong Force
    # alpha_s ~ 0.1179
    # 1/alpha_s ~ 8.48
    
    print(f"\n--- Checking Strong Force Candidates ---")
    print(f"Target 1/alpha_s: 8.48")
    print(f"Trace(G^2): {tr_g2}")
    
    ratio = tr_g2 / 8.48
    print(f"Ratio Trace/Target: {ratio}")
    
    # Check if Trace relates to 600 (number of vertices)
    print(f"Trace / 600: {tr_g2 / 600}")
    
    # Check other moments
    tr_g = np.sum(1.0 / (evals + lambda_1))
    print(f"Trace(G): {tr_g}")
    print(f"Trace(G) / 8.48: {tr_g / 8.48}")

if __name__ == "__main__":
    analyze_strong_force()
