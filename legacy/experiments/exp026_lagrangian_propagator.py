import numpy as np
import scipy.linalg as la

# Constants
PHI = (1 + np.sqrt(5)) / 2
PI = np.pi

def generate_600_cell_vertices():
    # ... (Same generation code as before for consistency) ...
    vertices = []
    # 1. Permutations of (+-1, +-1, +-1, +-1) / 2
    for i in range(16):
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
        vertices.append(signs * 0.5)
    # 2. Permutations of (+-2, 0, 0, 0) / 2
    for i in range(4):
        v = np.zeros(4)
        v[i] = 1.0
        vertices.append(v)
        vertices.append(-v)
    # 3. Even permutations of (+-phi, +-1, +-1/phi, 0) / 2
    v_base = np.array([PHI, 1.0, 1.0/PHI, 0.0]) * 0.5
    from itertools import permutations
    perms = list(permutations([0, 1, 2, 3]))
    even_perms = []
    for p in perms:
        swaps = 0
        temp_p = list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if temp_p[i] > temp_p[j]:
                    swaps += 1
        if swaps % 2 == 0:
            even_perms.append(p)
    for p in even_perms:
        perm_vals = v_base[list(p)]
        for i in range(16):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
            vertices.append(perm_vals * signs)
    vertices = np.array(vertices)
    _, idx = np.unique(np.round(vertices, 8), axis=0, return_index=True)
    return vertices[idx]

def build_laplacian(vertices):
    dists = la.norm(vertices[:, None, :] - vertices[None, :, :], axis=2)
    flat_dists = dists[dists > 0.0001]
    edge_length = np.min(flat_dists)
    adj_matrix = (np.abs(dists - edge_length) < 0.001).astype(float)
    degree_matrix = np.diag(np.sum(adj_matrix, axis=0))
    laplacian = degree_matrix - adj_matrix
    
    # Normalized Laplacian
    d_inv_sqrt = np.linalg.inv(np.sqrt(degree_matrix))
    norm_laplacian = np.dot(d_inv_sqrt, np.dot(laplacian, d_inv_sqrt))
    
    return norm_laplacian

def analyze_propagator():
    print("--- Constructing Scalar Field Theory on 600-cell ---")
    
    verts = generate_600_cell_vertices()
    L = build_laplacian(verts)
    
    # 1. Define Mass / Gap
    # We found lambda_1 = 1 / (2*phi^2)
    # Mass squared term m^2 corresponds to this energy scale
    
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_1 = eigenvalues[1] # First non-zero
    
    print(f"Spectral Gap (Mass^2 parameter): {lambda_1}")
    
    # 2. Construct Propagator G = (L + m^2)^-1
    # In QFT, Propagator is inverse of the kinetic operator (L) + mass term
    # Since L has a zero mode, we technically need the mass term to invert it properly,
    # or we compute the pseudo-inverse.
    # Here we use the physical mass gap we found.
    
    m_squared = lambda_1 # Using the gap itself as the mass term
    
    # Operator O = L + m^2 * I
    operator = L + m_squared * np.eye(len(L))
    
    # Propagator G
    G = np.linalg.inv(operator)
    
    print("\n--- Propagator Analysis (Green's Function) ---")
    
    # A. Trace of Propagator (Total Amplitude)
    # Sum of diagonal elements = Sum of 1 / (lambda_i + m^2)
    trace_G = np.trace(G)
    print(f"Trace(G) [Total Propagation]: {trace_G}")
    
    # Check against 20*phi^4
    target = 20 * PHI**4
    print(f"Target (20*phi^4): {target}")
    print(f"Ratio Trace/Target: {trace_G / target}")
    
    # B. Self-Interaction (Diagonal elements)
    # On a vertex-transitive graph like 600-cell, all diagonal elements should be equal
    diag_avg = np.mean(np.diag(G))
    print(f"Average Self-Interaction (Diagonal): {diag_avg}")
    
    # Check against 2*pi or alpha
    # Our formula had 2*pi*alpha.
    # Maybe G_diag ~ alpha?
    
    print(f"Inverse of Diagonal: {1/diag_avg}")
    
    # C. "Action" Value
    # S = Sum(G) over all elements?
    sum_all_G = np.sum(G)
    print(f"Sum of all Propagator entries: {sum_all_G}")
    
    # 3. The "Loop" Correction in QFT
    # One-loop correction is often Log(Det(Operator))
    
    log_det = np.sum(np.log(eigenvalues[1:] + m_squared)) # Exclude zero mode from det if pure L, but here we have mass
    log_det_full = np.linalg.slogdet(operator)[1]
    
    print(f"One-loop effective action (Log Det): {log_det_full}")
    print(f"Compare with 2*pi: {2*PI}")
    
    
    # 4. Hypothesis Check:
    # Does Trace(G) relate to 20*phi^4 directly?
    # We found earlier 1/alpha_bare = 5 * (1/lambda_1^2).
    # Trace(G) is sum(1/(lambda + m^2)).
    # If lambda_1 is dominant...
    
    # Let's calculate the "Bare Coupling" from this Lagrangian
    # beta_bare = Sum(G_ij squared)?
    
    frob_norm = np.linalg.norm(G, 'fro')**2
    print(f"Frobenius Norm Squared (Sum |G_ij|^2): {frob_norm}")
    
    print("\n--- Trying to recover 1/alpha ---")
    # Our derived formula: 5 * (1/lambda_1^2)
    # In terms of G, 1/(lambda+m)^2 is the eigenvalue of G^2.
    # So we are looking at Trace(G^2) or similar.
    
    trace_G2 = np.trace(np.dot(G, G))
    print(f"Trace(G^2) = Sum(1/(lambda_i + m^2)^2): {trace_G2}")
    
    # Remember our derivation: 5 * (1/lambda_1^2) = 137.08
    # Let's see what the full Trace(G^2) gives (including all high frequency modes)
    
    print(f"Target 137.08. Calculated Trace(G^2): {trace_G2}")
    
    # If Trace(G^2) is close to 137/5 = 27.4, then the "5" factor is indeed external (topological)
    # and the propagator handles the rest.
    
    print(f"Trace(G^2) * 5: {trace_G2 * 5}")

if __name__ == "__main__":
    analyze_propagator()
