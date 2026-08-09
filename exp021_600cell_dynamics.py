import numpy as np
from scipy.spatial.distance import pdist, squareform
import scipy.linalg as la

# Constants
PHI = (1 + np.sqrt(5)) / 2

def generate_600_cell_vertices():
    """
    Generates the 120 vertices of the 600-cell with edge length 1/PHI.
    Radius is 1.
    """
    vertices = []
    
    # 1. Permutations of (+-1, +-1, +-1, +-1) / 2  -> 16 vertices
    # Note: 0.5 * 2 = 1 (norm)
    for i in range(16):
        # binary representation to signs
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
        vertices.append(signs * 0.5)
        
    # 2. Permutations of (+-2, 0, 0, 0) / 2 -> (+-1, 0, 0, 0) -> 8 vertices
    # This creates the cross polytope vertices
    for i in range(4):
        v = np.zeros(4)
        v[i] = 1.0
        vertices.append(v)
        vertices.append(-v)
        
    # 3. Even permutations of (+-phi, +-1, +-1/phi, 0) / 2 -> 96 vertices
    # Base vector:
    v_base = np.array([PHI, 1.0, 1.0/PHI, 0.0]) * 0.5
    
    # We need to generate even permutations of coordinates and all sign combinations
    from itertools import permutations
    
    # Get unique permutations of indices
    perms = list(permutations([0, 1, 2, 3]))
    
    # Filter for even permutations
    # A permutation is even if it can be produced by an even number of swaps
    even_perms = []
    for p in perms:
        # Calculate number of swaps (bubble sort style)
        swaps = 0
        temp_p = list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if temp_p[i] > temp_p[j]:
                    swaps += 1
        if swaps % 2 == 0:
            even_perms.append(p)
            
    # Apply signs
    for p in even_perms:
        # Base permutation of values
        perm_vals = v_base[list(p)]
        
        # Iterate all 2^4 sign combinations
        for i in range(16):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
            candidate = perm_vals * signs
            
            # Check if already added (avoid duplicates if any, though geometrically distinct here)
            # Actually for this set, simpler logic:
            # Just generate them, we'll filter duplicates later if needed
            # For the specific (phi, 1, 1/phi, 0) set, all 96 are distinct
            vertices.append(candidate)

    # Convert to numpy array and unique to be safe
    vertices = np.array(vertices)
    
    # Remove duplicates (floating point tolerance)
    # Rounding for unique check
    _, idx = np.unique(np.round(vertices, 8), axis=0, return_index=True)
    vertices = vertices[idx]
    
    return vertices

def analyze_structure(vertices):
    print(f"Number of vertices found: {len(vertices)}")
    
    # Calculate all pairwise distances
    dists = squareform(pdist(vertices))
    
    # Find the edge length (shortest non-zero distance)
    # Filter 0 (self-distance)
    flat_dists = dists[dists > 0.0001]
    edge_length = np.min(flat_dists)
    
    print(f"Edge length (min dist): {edge_length:.6f}")
    print(f"In terms of Phi: 1/phi = {1/PHI:.6f}")
    
    # Build Adjacency Matrix
    # Connect if distance is approx edge_length
    adj_matrix = (np.abs(dists - edge_length) < 0.001).astype(float)
    
    degree = np.sum(adj_matrix, axis=0)
    avg_degree = np.mean(degree)
    print(f"Average Degree (neighbors per vertex): {avg_degree}")
    # Should be 12 for icosahedral vertex figure? No, 600-cell has 12 edges per vertex?
    # Let's verify. 20 tetrahedra meet at a vertex.
    # The vertex figure is an icosahedron. An icosahedron has 12 vertices.
    # So each vertex in 600-cell should have 12 neighbors.
    
    return adj_matrix

def compute_spectrum(adj_matrix):
    # Graph Laplacian: L = D - A
    degree_matrix = np.diag(np.sum(adj_matrix, axis=0))
    laplacian = degree_matrix - adj_matrix
    
    # Normalized Laplacian: L_norm = D^-1/2 * L * D^-1/2
    # This is better for spectral geometry
    d_inv_sqrt = np.linalg.inv(np.sqrt(degree_matrix))
    norm_laplacian = np.dot(d_inv_sqrt, np.dot(laplacian, d_inv_sqrt))
    
    eigenvalues = np.linalg.eigvalsh(norm_laplacian)
    
    # Sort
    eigenvalues = np.sort(eigenvalues)
    
    return eigenvalues

if __name__ == "__main__":
    verts = generate_600_cell_vertices()
    
    if len(verts) != 120:
        print(f"ERROR: Generated {len(verts)} vertices, expected 120.")
    else:
        adj = analyze_structure(verts)
        
        print("\n--- Calculating Laplacian Spectrum ---")
        evals = compute_spectrum(adj)
        
        print("Smallest 10 eigenvalues:")
        print(evals[:10])
        
        print("\nLargest 10 eigenvalues:")
        print(evals[-10:])
        
        # Spectral Gap
        gap = evals[1] - evals[0] # evals[0] is usually 0
        print(f"\nSpectral Gap (Mass Gap?): {gap}")
        
        # Search for alpha or phi ratios
        # Alpha is approx 1/137 ~ 0.00729
        
        print("\nSearching for physical ratios in spectrum...")
        # Check gap vs max
        ratio1 = gap / np.max(evals)
        print(f"Gap / Max_Eigenvalue: {ratio1}")
        
        # Check if any eigenvalue corresponds to known numbers
        # The eigenvalues of the graph Laplacian are related to the energies of free particles on the graph.
