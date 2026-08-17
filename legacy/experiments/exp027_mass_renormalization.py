import numpy as np
import scipy.linalg as la
from scipy.optimize import newton

# Constants
PHI = (1 + np.sqrt(5)) / 2
ALPHA_INV_TARGET = 137.035999084 # Experimental value
TERM_GEOMETRIC = 20 * PHI**4      # 137.082039...

def generate_600_cell_vertices():
    # Same generation logic
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

def get_eigenvalues():
    verts = generate_600_cell_vertices()
    dists = la.norm(verts[:, None, :] - verts[None, :, :], axis=2)
    edge_length = np.min(dists[dists > 0.0001])
    adj = (np.abs(dists - edge_length) < 0.001).astype(float)
    deg = np.diag(np.sum(adj, axis=0))
    lap = deg - adj
    # Normalized Laplacian
    d_inv_sqrt = np.linalg.inv(np.sqrt(deg))
    norm_lap = np.dot(d_inv_sqrt, np.dot(lap, d_inv_sqrt))
    evals = np.linalg.eigvalsh(norm_lap)
    return evals

def trace_g2(m_squared, eigenvalues):
    # Tr(G^2) = Sum(1 / (lambda_i + m^2)^2)
    return np.sum(1.0 / (eigenvalues + m_squared)**2)

def find_renormalized_mass():
    print("--- EXP-027: Mass Renormalization ---")
    
    # 1. Get Spectrum
    evals = get_eigenvalues()
    lambda_1 = evals[1] # Gap
    print(f"Base Spectral Gap (lambda_1): {lambda_1:.8f}")
    
    # Check current Trace with base mass
    current_trace = trace_g2(lambda_1, evals)
    print(f"Trace(G^2) with lambda_1: {current_trace:.4f}")
    print(f"Target Value (Geometric 20*phi^4): {TERM_GEOMETRIC:.4f}")
    
    # 2. Optimization Loop
    # We need to find m_sq such that trace_g2(m_sq) = TERM_GEOMETRIC
    # Function to find root of:
    func = lambda m_sq: trace_g2(m_sq, evals) - TERM_GEOMETRIC
    
    # Derivative (for Newton method): d/dm^2 [ (lam + m^2)^-2 ] = -2 * (lam + m^2)^-3
    dfunc = lambda m_sq: -2 * np.sum(1.0 / (evals + m_sq)**3)
    
    # Solve
    m_renorm_sq = newton(func, x0=lambda_1, fprime=dfunc)
    
    print(f"\n--- Renormalization Results ---")
    print(f"Renormalized Mass Squared (m*^2): {m_renorm_sq:.8f}")
    
    # 3. Analyze the Ratio (Z factor)
    Z_factor = m_renorm_sq / lambda_1
    print(f"Renormalization Factor Z (m*^2 / lambda_1): {Z_factor:.8f}")
    
    # 4. Search for Geometric Meaning of Z
    # Is Z related to alpha, phi, or curvature?
    
    print("\n--- Searching for Meaning of Z ---")
    print(f"Z = {Z_factor}")
    
    # Check simple fractions involving Phi
    print(f"phi: {PHI}")
    print(f"sqrt(phi): {np.sqrt(PHI)}")
    print(f"phi^2/2: {PHI**2/2}")
    
    # Check Curvature related
    # In 4D, curvature scaling often involves factors of 1.05... or 0.95...
    
    diff_from_1 = Z_factor - 1
    print(f"Deviation from 1 (Mass correction): {diff_from_1:.8f}")
    
    # Check relation to Alpha
    # alpha ~ 0.007
    print(f"Ratio Deviation/Alpha: {diff_from_1 / (1/137.036)}")
    
    # Check relation to 1/12 (curvature of sphere vs flat) or similar
    
    return m_renorm_sq, Z_factor

if __name__ == "__main__":
    find_renormalized_mass()
