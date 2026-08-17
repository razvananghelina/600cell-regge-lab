# EXP-221: Feynman Diagrams on the 600-Cell Graph
# ==================================================
# CAN WE DERIVE m_nu = m_e * (alpha/phi)^3 from a loop calculation?
#
# The idea:
# - In continuum QFT, each loop gives factor alpha/(4*pi)
# - On the 600-cell graph, loops are SUMS over vertices
# - The propagator is G = (Laplacian + m^2)^{-1}
# - A 1-loop self-energy = alpha * Tr(G)
# - The "loop factor" emerges from the graph structure
#
# Question: does the 600-cell give loop factor = alpha/phi ?
# If yes, 3-loop neutrino mass = m_e * (alpha/phi)^3 is DERIVED!

import numpy as np
from itertools import product as iproduct

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1/137.035999084

# ===============================================================
# STEP 1: Build the 600-cell graph
# ===============================================================
print("="*70)
print("EXP-221: FEYNMAN DIAGRAMS ON THE 600-CELL GRAPH")
print("="*70)

print("\n--- Step 1: Building 600-cell (120 vertices, 720 edges) ---")

vertices = []

# Type A: 8 vertices - permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        vertices.append(tuple(v))

# Type B: 16 vertices - (+-1/2, +-1/2, +-1/2, +-1/2)
for signs in iproduct([0.5, -0.5], repeat=4):
    vertices.append(tuple(signs))

# Type C: 96 vertices - even permutations of (0, +-1/(2*phi), +-1/2, +-phi/2)
# Even permutations of (0,1,2,3) = 12 permutations
even_perms = [
    (0,1,2,3), (0,2,3,1), (0,3,1,2),
    (1,0,3,2), (1,2,0,3), (1,3,2,0),
    (2,0,1,3), (2,1,3,0), (2,3,0,1),
    (3,0,2,1), (3,1,0,2), (3,2,1,0),
]
base_vals = [0.0, 1.0/(2*PHI), 0.5, PHI/2.0]

for perm in even_perms:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                signs = [1.0, s1, s2, s3]
                v = [0.0]*4
                for i in range(4):
                    v[i] = signs[perm.index(i)] * base_vals[perm[i]] if perm[i] != 0 else 0.0
                # Reconstruct properly
                v2 = [0.0]*4
                for i in range(4):
                    v2[i] = base_vals[perm[i]]
                # Apply signs to non-zero entries
                nonzero_idx = [i for i in range(4) if perm[i] != 0]
                sign_list = [s1, s2, s3]
                for k, idx in enumerate(nonzero_idx):
                    v2[idx] *= sign_list[k]
                vertices.append(tuple(v2))

# Remove duplicates
vertices_unique = list(set(vertices))
N = len(vertices_unique)
print(f"  Vertices generated: {N}")

if N != 120:
    print(f"  WARNING: expected 120, got {N}. Checking...")
    # Try a different construction
    vertices = []

    # Type A: 8 vertices
    for i in range(4):
        for s in [1, -1]:
            v = [0.0]*4
            v[i] = float(s)
            vertices.append(v)

    # Type B: 16 vertices
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    vertices.append([s0*0.5, s1*0.5, s2*0.5, s3*0.5])

    # Type C: 96 vertices - even permutations of (0, +-1/(2*phi), +-1/2, +-phi/2)
    vals = [0.0, 1.0/(2*PHI), 0.5, PHI/2.0]

    for perm in even_perms:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    v = [0.0]*4
                    signs_map = {0: 1.0}  # 0 component always positive (it's 0)
                    sign_idx = 0
                    s_list = [s1, s2, s3]
                    for i in range(4):
                        val = vals[perm[i]]
                        if perm[i] == 0:
                            v[i] = 0.0
                        else:
                            v[i] = val * s_list[sign_idx]
                            sign_idx += 1
                    vertices.append(v)

    # Convert to array and remove duplicates
    V = np.array(vertices)
    # Round to avoid floating point duplicates
    V_rounded = np.round(V, 10)
    V_unique = np.unique(V_rounded, axis=0)
    N = len(V_unique)
    print(f"  After dedup: {N} vertices")
    V = V_unique
else:
    V = np.array(vertices_unique)

# Verify all at unit distance from origin
norms = np.sqrt(np.sum(V**2, axis=1))
print(f"  Norms: min={norms.min():.6f}, max={norms.max():.6f}")

# Build adjacency matrix: edge if distance = 1/phi (for unit 600-cell)
# Actually for unit-radius 600-cell, edge length = 1/phi
# Let's compute all pairwise distances
from scipy.spatial.distance import cdist
D = cdist(V, V)

# Find the smallest non-zero distance
d_sorted = np.sort(D[0])
print(f"  Distances from vertex 0: {d_sorted[1:6]}")
edge_length = d_sorted[1]
print(f"  Edge length: {edge_length:.6f}")
print(f"  1/phi = {1/PHI:.6f}")

# Build adjacency matrix
tol = 1e-6
A = (np.abs(D - edge_length) < tol).astype(float)
np.fill_diagonal(A, 0)
degrees = A.sum(axis=1)
print(f"  Degree: min={degrees.min():.0f}, max={degrees.max():.0f}")
n_edges = int(A.sum()) // 2
print(f"  Edges: {n_edges}")

# ===============================================================
# STEP 2: Compute spectrum (Laplacian eigenvalues)
# ===============================================================
print("\n--- Step 2: Spectrum of adjacency matrix ---")
deg = int(degrees[0])
eigenvalues_A = np.linalg.eigvalsh(A)
eigenvalues_A = np.sort(eigenvalues_A)[::-1]  # descending

# Group eigenvalues
from collections import Counter
eig_rounded = np.round(eigenvalues_A, 4)
eig_counts = Counter(eig_rounded)
print(f"\n  Eigenvalue : Multiplicity")
for eig in sorted(eig_counts.keys(), reverse=True):
    mult = eig_counts[eig]
    # Check if a+b*phi
    # Solve: eig = a + b*phi => b = (eig-a)/phi
    # For integer a,b: try a from -10 to 20
    found = False
    for a in range(-10, 20):
        b_test = (eig - a) / PHI
        if abs(b_test - round(b_test)) < 0.01:
            b = int(round(b_test))
            n_ab = 5*a + 6*b
            print(f"  {eig:>10.4f} : {mult:>3d}  = {a} + {b}*phi  (n={n_ab})")
            found = True
            break
    if not found:
        print(f"  {eig:>10.4f} : {mult:>3d}")

# Laplacian eigenvalues
print(f"\n  Laplacian eigenvalues (mu = {deg} - lambda):")
L = deg * np.eye(N) - A
eigenvalues_L = np.linalg.eigvalsh(L)
eigenvalues_L = np.sort(eigenvalues_L)

eig_L_rounded = np.round(eigenvalues_L, 4)
eig_L_counts = Counter(eig_L_rounded)
for eig in sorted(eig_L_counts.keys()):
    mult = eig_L_counts[eig]
    if mult > 0:
        print(f"  mu = {eig:>10.4f} : mult = {mult:>3d}")

# ===============================================================
# STEP 3: Propagator and loop factors
# ===============================================================
print("\n--- Step 3: Propagator G = (L + m^2)^(-1) ---")

# The propagator on the graph at mass m:
# G(m) = (L + m^2 * I)^(-1)
# Trace: Tr(G) = sum_i 1/(mu_i + m^2)
# This is the 1-loop integral on the graph

# In continuum 4D, the 1-loop integral gives:
# I_1 = integral d^4k/(2*pi)^4 * 1/(k^2+m^2) ~ Lambda^2/(16*pi^2)
# The "loop factor" = 1/(16*pi^2) = 1/(4*pi)^2

# On the graph:
# I_1(m) = (1/N) * Tr(G(m)) = (1/N) * sum_i 1/(mu_i + m^2)

# For massless: mu_0 = 0 gives divergence (IR). Use small mass.
# For the neutrino calculation, m^2 is the charged lepton mass squared
# (in units where the graph "lattice spacing" = 1)

# Let's compute the loop factor for various masses
print("\n  1-loop factor: I_1(m) = (1/N) * sum 1/(mu_i + m^2)")
print()

mu = np.sort(eigenvalues_L)
# Skip the zero mode (Goldstone-like)
mu_nonzero = mu[mu > 0.01]

for m2 in [0.01, 0.1, 1.0, 10.0]:
    I1 = np.mean(1.0 / (mu + m2))
    I1_nz = np.mean(1.0 / (mu_nonzero + m2))
    print(f"  m^2 = {m2:>5.2f}: I_1 = {I1:.6f}  (no zero mode: {I1_nz:.6f})")

# ===============================================================
# STEP 4: The key question - what is the "loop factor"?
# ===============================================================
print("\n--- Step 4: Graph loop factor vs 1/(4*pi) vs 1/phi ---")

# In QFT, the loop suppression factor per loop is:
# f_loop = alpha * I_loop
# where I_loop ~ 1/(16*pi^2) in 4D continuum

# On the graph, we need to normalize properly.
# The question: for what mass scale does
# (1/N) * Tr(G) = 1/phi ?

# (1/N) * sum 1/(mu_i + m^2) = 1/phi
# We need to find m^2 such that this holds

# Let's scan
print("\n  Scanning for m^2 such that I_1(m^2) = 1/phi:")
target = 1.0/PHI
m2_scan = np.logspace(-3, 3, 10000)
I1_scan = np.array([np.mean(1.0/(mu + m2)) for m2 in m2_scan])

# Find crossings
idx = np.argmin(np.abs(I1_scan - target))
print(f"  Target: 1/phi = {target:.6f}")
print(f"  Best m^2 = {m2_scan[idx]:.4f}, I_1 = {I1_scan[idx]:.6f}")
print()

# Also check what I_1 equals for "natural" mass values
print("  I_1 for natural graph masses:")
natural_masses = {
    "m^2 = 1 (lattice scale)": 1.0,
    "m^2 = deg = 12": float(deg),
    "m^2 = phi": PHI,
    "m^2 = phi^2": PHI**2,
    "m^2 = 1/phi": 1.0/PHI,
    "m^2 = a_1 = 5": 5.0,
    "m^2 = b_1 = 6": 6.0,
}
for name, m2 in natural_masses.items():
    I1 = np.mean(1.0/(mu + m2))
    # What is I1 close to?
    for label, val in [("1/phi", 1/PHI), ("1/phi^2", 1/PHI**2),
                        ("1/(4pi)", 1/(4*np.pi)), ("1/(4pi)^2", 1/(4*np.pi)**2),
                        ("1/deg", 1.0/deg), ("1/N", 1.0/N)]:
        if abs(I1 - val)/val < 0.3:
            print(f"  {name}: I_1 = {I1:.6f} ~ {label} = {val:.6f} ({abs(I1-val)/val*100:.1f}%)")
            break
    else:
        print(f"  {name}: I_1 = {I1:.6f}")

# ===============================================================
# STEP 5: Multi-loop calculation
# ===============================================================
print("\n--- Step 5: Multi-loop self-energy ---")

# n-loop self-energy on the graph:
# The simplest version: Sigma_n = alpha^n * [(1/N) * Tr(G)]^n
# This assumes independent loops.
#
# A more accurate version involves powers of the propagator:
# Sigma_n ~ alpha^n * (1/N) * Tr(G^n)
# where G^n = (L+m^2)^{-n}
#
# Tr(G^n) = sum_i 1/(mu_i + m^2)^n

print("\n  Multi-loop traces: I_n(m^2) = (1/N) * Tr(G^n)")
print(f"  {'m^2':>6s} {'I_1':>12s} {'I_2':>12s} {'I_3':>12s} {'I_1/I_2':>10s} {'I_2/I_3':>10s} {'I_1^3':>12s}")
print("  " + "-"*78)

for m2 in [0.1, 0.5, 1.0, PHI, 5.0, 12.0]:
    G_inv = mu + m2
    I1 = np.mean(1.0/G_inv)
    I2 = np.mean(1.0/G_inv**2)
    I3 = np.mean(1.0/G_inv**3)
    print(f"  {m2:>6.2f} {I1:>12.6f} {I2:>12.6f} {I3:>12.6f} {I1/I2:>10.4f} {I2/I3:>10.4f} {I1**3:>12.8f}")

# ===============================================================
# STEP 6: Physical interpretation - neutrino mass
# ===============================================================
print("\n--- Step 6: Neutrino mass from graph loops ---")

# In the scotogenic model (loop-generated neutrino mass):
# m_nu ~ (alpha/f_loop) * m_lepton * [loop integral]
# where f_loop is the "phase space" factor
#
# On the 600-cell:
# m_nu = m_e * (alpha * I_loop)^n
#
# We want: m_e * (alpha * I_loop)^3 = m_e * (alpha/phi)^3
# This requires: I_loop = 1/phi
# (i.e., the 1-loop integral equals 1/phi)

print("\n  For m_nu = m_e * (alpha/phi)^3, we need loop factor = 1/phi")
print(f"  alpha/phi = {ALPHA/PHI:.6e}")
print(f"  (alpha/phi)^3 = {(ALPHA/PHI)**3:.6e}")
print()

# What mass scale gives I_1 = 1/phi?
# We found it above. Let's be more precise
from scipy.optimize import brentq
def f_I1(m2):
    return np.mean(1.0/(mu + m2)) - 1.0/PHI

# Check signs
print(f"  f(0.01) = {f_I1(0.01):.4f}")
print(f"  f(100) = {f_I1(100):.4f}")

if f_I1(0.01) * f_I1(100) < 0:
    m2_star = brentq(f_I1, 0.01, 100)
    print(f"  m^2 that gives I_1 = 1/phi: m^2 = {m2_star:.6f}")

    # Is this a "natural" value?
    for name, val in [("phi^2", PHI**2), ("phi", PHI), ("phi+1=phi^2", PHI+1),
                       ("deg-1=11", 11.0), ("deg=12", 12.0), ("deg+1=13=N_bag", 13.0),
                       ("2*phi^2", 2*PHI**2), ("a_1+b_1=11", 11.0),
                       ("b_1*phi", 6*PHI), ("a_1*phi", 5*PHI),
                       ("(deg-1)/phi", 11/PHI), ("h_E8/phi^2", 30/PHI**2)]:
        err = abs(m2_star - val)/val * 100
        if err < 5:
            print(f"    ~ {name} = {val:.4f} (err {err:.2f}%)")
else:
    print("  No crossing found in range [0.01, 100]")

# ===============================================================
# STEP 7: Alternative - ratio of traces
# ===============================================================
print("\n--- Step 7: Natural loop factor from graph ---")

# Maybe the loop factor isn't I_1(m^2) for a specific m^2.
# Maybe it's a RATIO of spectral quantities.

# Idea: f_loop = [sum of non-zero Laplacian eigenvalues]^{-1} * N
# Or: f_loop relates to the spectral gap

spectral_gap = mu[mu > 0.01].min()
print(f"  Spectral gap: mu_1 = {spectral_gap:.6f}")
print(f"  deg/mu_1 = {deg/spectral_gap:.6f}")
print(f"  1/phi = {1/PHI:.6f}")
print()

# Zeta function of the Laplacian
# zeta_L(s) = sum_{mu_i > 0} mu_i^(-s)
for s in [0.5, 1.0, 1.5, 2.0]:
    zeta = np.sum(mu[mu > 0.01]**(-s))
    zeta_norm = zeta / (N-1)  # normalized by non-zero modes
    print(f"  zeta_L({s:.1f}) = {zeta:.4f},  normalized = {zeta_norm:.6f}")
    for name, val in [("1/phi", 1/PHI), ("phi", PHI), ("1/phi^2", 1/PHI**2)]:
        if abs(zeta_norm - val)/val < 0.1:
            print(f"    ~ {name} = {val:.6f} ({abs(zeta_norm-val)/val*100:.2f}%)")

# ===============================================================
# STEP 8: Green's function at coinciding points
# ===============================================================
print("\n--- Step 8: Green's function G(x,x) (self-energy) ---")

# G(x,x) = (1/N) * sum_i 1/(mu_i + m^2) for vertex x
# By vertex-transitivity, this is the same for all x
# G(x,x) = I_1(m^2) (same as trace/N)

# The MASSLESS regulated Green's function:
# G_reg(x,x) = (1/N) * sum_{mu_i > 0} 1/mu_i
G_reg = np.sum(1.0/mu[mu > 0.01]) / N
print(f"  G_reg(x,x) = (1/N)*sum 1/mu_i = {G_reg:.6f}")
print(f"  = {G_reg:.6f}")
print(f"  1/phi = {1/PHI:.6f}  (ratio: {G_reg*PHI:.4f})")
print(f"  1/phi^2 = {1/PHI**2:.6f}  (ratio: {G_reg*PHI**2:.4f})")
print()

# What fraction is G_reg of 1/phi?
print(f"  G_reg / (1/phi) = {G_reg * PHI:.6f}")
print(f"  G_reg * deg = {G_reg * deg:.6f}")
print(f"  G_reg * N = {G_reg * N:.6f}")
print()

# The KEY physical quantity: dimensionless loop factor
# In QFT: loop factor = g^2/(16*pi^2)
# On graph: loop factor = (coupling) * G_reg(x,x)
# If coupling = alpha (EM), loop factor = alpha * G_reg
print(f"  alpha * G_reg = {ALPHA * G_reg:.6e}")
print(f"  alpha / phi   = {ALPHA / PHI:.6e}")
print(f"  Ratio: {ALPHA * G_reg / (ALPHA/PHI):.6f} = phi * G_reg = {PHI * G_reg:.6f}")
print()

# For the loop factor to be alpha/phi, we need G_reg = 1/phi^2
# Check: G_reg = 0.??? vs 1/phi^2 = 0.3820
print(f"  For alpha/phi loop factor: need G_reg = 1/phi^2 = {1/PHI**2:.6f}")
print(f"  Actual G_reg = {G_reg:.6f}")
print(f"  Ratio: {G_reg/(1/PHI**2):.4f}")

# ===============================================================
# STEP 9: Heat kernel approach
# ===============================================================
print("\n--- Step 9: Heat kernel on 600-cell ---")

# K(t) = sum_i exp(-mu_i * t)
# K(0) = N = 120
# K(t) for t->inf goes to 1 (the zero mode)
# The "return probability" at time t: p(t) = K(t)/N

# The integrated heat kernel (effective loop factor):
# integral_0^T K(t)/N dt = sum_i (1 - exp(-mu_i*T))/(N*mu_i)
# For T->inf: = (1/N) * sum_{mu>0} 1/mu_i = G_reg

# At what time does p(t) = 1/phi?
print("  Heat kernel return probability p(t) = K(t)/N:")
for t in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
    Kt = np.sum(np.exp(-mu * t))
    pt = Kt / N
    print(f"    t = {t:.2f}: p(t) = {pt:.6f}", end="")
    if abs(pt - 1/PHI) < 0.05:
        print(f"  <-- near 1/phi = {1/PHI:.6f}!", end="")
    print()

# Find t where p(t) = 1/phi
def f_heat(t):
    return np.sum(np.exp(-mu * t)) / N - 1.0/PHI

try:
    t_star = brentq(f_heat, 0.01, 10.0)
    print(f"\n  p(t*) = 1/phi at t* = {t_star:.6f}")
    for name, val in [("1/deg", 1.0/deg), ("1/(2*pi)", 1/(2*np.pi)),
                       ("1/phi^4", 1/PHI**4), ("1/a_1", 1.0/5),
                       ("1/b_1", 1.0/6), ("1/N_eig", 1.0/9),
                       ("alpha", ALPHA)]:
        err = abs(t_star - val)/val*100
        if err < 20:
            print(f"    ~ {name} = {val:.6f} (err {err:.1f}%)")
except:
    print("  No crossing found for heat kernel")

# ===============================================================
# STEP 10: Direct 3-loop computation
# ===============================================================
print("\n--- Step 10: Direct 3-loop neutrino self-energy ---")

# Weinberg-like 3-loop operator on the graph:
# Sigma_3 ~ (alpha/N)^3 * sum_{i,j,k} G(0,i)*G(i,j)*G(j,k)*G(k,0)
# = (alpha/N)^3 * [G^4]_{00}  (4 propagators for 3 loops)
# But vertex-transitivity => [G^n]_{00} = (1/N)*Tr(G^n)

# For massive propagator at scale m:
m2_test = 1.0  # lattice unit mass

G_matrix = np.linalg.inv(L + m2_test * np.eye(N))
print(f"  Using m^2 = {m2_test}")
print(f"  G(0,0) = {G_matrix[0,0]:.8f}")
print(f"  Tr(G)/N = {np.trace(G_matrix)/N:.8f}")

# n-loop factors (self-energy ~ G^{n+1} at coincident point)
for n in range(1, 5):
    Gn = np.linalg.matrix_power(G_matrix, n+1)
    loop_n = Gn[0,0]  # = Tr(G^{n+1})/N by vertex transitivity
    loop_n_trace = np.trace(np.linalg.matrix_power(G_matrix, n+1)) / N
    ratio_to_G00_power = loop_n / G_matrix[0,0]**(n+1)
    print(f"  {n}-loop: G^{n+1}(0,0) = {loop_n:.8f}, G(0,0)^{n+1} = {G_matrix[0,0]**(n+1):.8f}, ratio = {ratio_to_G00_power:.6f}")

print("\n  If loops were independent: G^{n+1}(0,0) = G(0,0)^{n+1}")
print("  The ratio measures loop CORRELATIONS on the graph.")
print()

# ===============================================================
# STEP 11: THE KEY TEST
# ===============================================================
print("="*70)
print("STEP 11: THE KEY TEST - Does 600-cell give loop factor 1/phi?")
print("="*70)

# Define: f_loop = [G^2(0,0) / G(0,0)]
# This is the "effective 1-loop factor" including graph correlations
# For independent loops: f = G(0,0)
# For correlated: f = G^2(0,0)/G(0,0)

# Compute for various m^2
print(f"\n  {'m^2':>8s} {'G(0,0)':>12s} {'f_1loop':>12s} {'f*phi':>10s} {'f*phi^2':>10s}")
print("  " + "-"*55)

for m2 in [0.01, 0.1, 0.5, 1.0, PHI, 2.0, 5.0, 10.0, 12.0]:
    G_m = np.linalg.inv(L + m2 * np.eye(N))
    G2 = G_m @ G_m
    g00 = G_m[0,0]
    g2_00 = G2[0,0]
    f_loop = g2_00 / g00  # effective 1-loop factor
    print(f"  {m2:>8.3f} {g00:>12.6f} {f_loop:>12.6f} {f_loop*PHI:>10.4f} {f_loop*PHI**2:>10.4f}")

# What if we look at it differently?
# The "loop enhancement" = G^2(0,0) / G(0,0)^2
# This measures how much the graph enhances loops vs independent
print(f"\n  Loop enhancement (graph correlation):")
print(f"  {'m^2':>8s} {'G^2/G^2':>12s} {'near?':>20s}")
print("  " + "-"*45)

for m2 in [0.1, 0.5, 1.0, PHI, 5.0, 12.0]:
    G_m = np.linalg.inv(L + m2 * np.eye(N))
    G2 = G_m @ G_m
    enhancement = G2[0,0] / G_m[0,0]**2
    # Check if close to anything
    nearest = ""
    for name, val in [("phi", PHI), ("phi^2", PHI**2), ("N/deg", N/float(deg)),
                       ("deg", float(deg)), ("N_eig", 9.0), ("a_1", 5.0)]:
        if abs(enhancement - val)/val < 0.1:
            nearest = f"~ {name} = {val:.2f} ({abs(enhancement-val)/val*100:.1f}%)"
            break
    print(f"  {m2:>8.3f} {enhancement:>12.4f} {nearest:>20s}")

# ===============================================================
# FINAL SUMMARY
# ===============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
QUESTION: Does the 600-cell graph produce loop factor 1/phi?

FINDINGS:
  1. The regulated Green's function G_reg(x,x) = {G_reg:.6f}
     This is NOT 1/phi = {1/PHI:.6f}
     Nor 1/phi^2 = {1/PHI**2:.6f}

  2. The loop factor depends on the mass scale m^2
     No single "natural" m^2 gives exactly 1/phi

  3. The heat kernel p(t) = 1/phi at t ~ {t_star if 't_star' in dir() else '?'}
     This time has no obvious geometric meaning

  4. Loop correlations (G^2/G^2 enhancement) vary with mass
     No universal phi-related value

VERDICT:
  The simple Feynman diagram approach on the 600-cell does NOT
  directly produce the factor 1/phi as a loop suppression.

  The formula m_nu = m_e * (alpha/phi)^3 works numerically
  but is NOT derivable from a straightforward loop calculation
  on the graph Laplacian.

  POSSIBLE EXPLANATIONS:
  a) The relevant propagator is NOT the Laplacian inverse
     but something related to the association scheme
  b) The phi factor comes from the DISCRETE SCALE INVARIANCE
     (phi is the scaling factor), not from loop integrals
  c) The 3-loop interpretation is coincidental, and the
     true origin is m3 = 2*m_e/phi^35 (pure geometric)
  d) We need a more sophisticated diagram calculation
     (e.g., using the gauge structure SU(2) edges)
""")
