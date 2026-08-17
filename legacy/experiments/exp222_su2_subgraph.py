# EXP-222: Neutrino Mass from SU(2) Subgraph + Spectral Identities
# ==================================================================
# The Feynman approach on the FULL Laplacian failed (EXP-221).
#
# NEW IDEA: Neutrinos couple ONLY to SU(2) (weak force).
# The relevant propagator should use the SU(2) SUBGRAPH, not the full graph!
#
# Plan:
# Part A: Build SU(2) subgraph (288 edges) and compute its spectrum
# Part B: Check if SU(2) propagator gives phi-related loop factors
# Part C: Spectral identities - can we derive 35 = h + a_1?
# Part D: Random walk / return probability approach
# Part E: Galois-broken DSI (Majorana mechanism)

import numpy as np
from itertools import product as iproduct
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1/137.035999084
ALPHA_S = 0.1179
SIN2TW = 6.0/26.0

m_e = 0.51099895e6
m3_exp = np.sqrt(2.455e-3)
m2_exp = np.sqrt(7.53e-5)

# Invariants
a_1 = 5; b_1 = 6; N_eig = 9; deg = 12
h_E8 = 30; dim_R4 = 4; N_bag = 13; N_gen = 3

print("="*70)
print("EXP-222: SU(2) SUBGRAPH + SPECTRAL IDENTITIES")
print("="*70)

# ===============================================================
# Build 600-cell
# ===============================================================
print("\n--- Building 600-cell ---")
vertices = []
for i in range(4):
    for s in [1, -1]:
        v = [0.0]*4; v[i] = float(s); vertices.append(v)
for s0 in [1,-1]:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                vertices.append([s0*0.5, s1*0.5, s2*0.5, s3*0.5])
even_perms = [
    (0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
    (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)]
vals = [0.0, 1.0/(2*PHI), 0.5, PHI/2.0]
for perm in even_perms:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                v = [0.0]*4; sign_idx = 0; s_list = [s1,s2,s3]
                for i in range(4):
                    if perm[i] == 0: v[i] = 0.0
                    else: v[i] = vals[perm[i]]*s_list[sign_idx]; sign_idx += 1
                vertices.append(v)
V = np.unique(np.round(np.array(vertices), 10), axis=0)
N = len(V)
print(f"  Vertices: {N}")

from scipy.spatial.distance import cdist
D_dist = cdist(V, V)
edge_length = np.sort(D_dist[0])[1]
tol = 1e-6
A_full = (np.abs(D_dist - edge_length) < tol).astype(float)
np.fill_diagonal(A_full, 0)
n_edges = int(A_full.sum()) // 2
print(f"  Edges: {n_edges}, degree: {int(A_full.sum(axis=1)[0])}")

# ===============================================================
# PART A: Classify edges into U(1), SU(2), SU(3)
# ===============================================================
print("\n" + "="*70)
print("PART A: EDGE CLASSIFICATION")
print("="*70)

# For each edge (i,j), count shared neighbors and their types
# Vertex types: A (8 with one coord +-1), B (16 with all +-1/2), C (96 golden)
vertex_types = []
for i in range(N):
    v = V[i]
    n_zero = np.sum(np.abs(v) < tol)
    n_half = np.sum(np.abs(np.abs(v) - 0.5) < tol)
    n_one = np.sum(np.abs(np.abs(v) - 1.0) < tol)
    if n_zero == 3 and n_one == 1:
        vertex_types.append('A')
    elif n_half == 4:
        vertex_types.append('B')
    else:
        vertex_types.append('C')

type_counts = Counter(vertex_types)
print(f"  Vertex types: A={type_counts['A']}, B={type_counts['B']}, C={type_counts['C']}")

# Classify each edge
edge_classes = {'U1': [], 'SU2': [], 'SU3': []}
edge_list = []
for i in range(N):
    for j in range(i+1, N):
        if A_full[i,j] > 0.5:
            # Count shared neighbors
            shared = np.where((A_full[i,:] > 0.5) & (A_full[j,:] > 0.5))[0]
            n_shared = len(shared)
            # Types of shared neighbors
            shared_A = sum(1 for s in shared if vertex_types[s] == 'A')
            shared_B = sum(1 for s in shared if vertex_types[s] == 'B')
            shared_C = sum(1 for s in shared if vertex_types[s] == 'C')

            edge_type = vertex_types[i] + vertex_types[j]

            # Classification based on shared neighbor profile
            if shared_A >= 1 and n_shared <= 2:
                edge_classes['U1'].append((i,j))
            elif shared_B >= 2 and shared_A == 0:
                edge_classes['SU2'].append((i,j))
            else:
                edge_classes['SU3'].append((i,j))
            edge_list.append((i,j, n_shared, shared_A, shared_B, shared_C))

print(f"  U(1) edges:  {len(edge_classes['U1'])}")
print(f"  SU(2) edges: {len(edge_classes['SU2'])}")
print(f"  SU(3) edges: {len(edge_classes['SU3'])}")
print(f"  Total:       {sum(len(v) for v in edge_classes.values())}")

# Build SU(2) adjacency matrix
A_SU2 = np.zeros((N, N))
for i,j in edge_classes['SU2']:
    A_SU2[i,j] = 1.0
    A_SU2[j,i] = 1.0
deg_SU2 = A_SU2.sum(axis=1)
print(f"\n  SU(2) subgraph degree: min={deg_SU2.min():.0f}, max={deg_SU2.max():.0f}, mean={deg_SU2.mean():.1f}")

# Build U(1) adjacency matrix
A_U1 = np.zeros((N, N))
for i,j in edge_classes['U1']:
    A_U1[i,j] = 1.0
    A_U1[j,i] = 1.0
deg_U1 = A_U1.sum(axis=1)
print(f"  U(1) subgraph degree: min={deg_U1.min():.0f}, max={deg_U1.max():.0f}")

# Build SU(3) adjacency matrix
A_SU3 = np.zeros((N, N))
for i,j in edge_classes['SU3']:
    A_SU3[i,j] = 1.0
    A_SU3[j,i] = 1.0
deg_SU3 = A_SU3.sum(axis=1)
print(f"  SU(3) subgraph degree: min={deg_SU3.min():.0f}, max={deg_SU3.max():.0f}")

# ===============================================================
# PART B: SU(2) subgraph spectrum
# ===============================================================
print("\n" + "="*70)
print("PART B: SU(2) SUBGRAPH SPECTRUM")
print("="*70)

eigs_SU2 = np.linalg.eigvalsh(A_SU2)
eigs_SU2_sorted = np.sort(eigs_SU2)[::-1]

eigs_r = np.round(eigs_SU2_sorted, 3)
eig_counts_SU2 = Counter(eigs_r)
print(f"\n  SU(2) adjacency eigenvalues:")
for eig in sorted(eig_counts_SU2.keys(), reverse=True):
    mult = eig_counts_SU2[eig]
    if mult > 0:
        # Check phi decomposition
        found = ""
        for a in range(-15, 15):
            b_test = (eig - a) / PHI
            if abs(b_test - round(b_test)) < 0.05:
                b = int(round(b_test))
                found = f"= {a}+{b}*phi"
                break
        print(f"    {eig:>10.3f} : mult {mult:>3d}  {found}")

# Laplacian of SU(2) subgraph
# Note: SU(2) subgraph is NOT regular (varying degree)
# Use normalized Laplacian or unnormalized
d_su2 = deg_SU2.mean()
L_SU2 = np.diag(deg_SU2) - A_SU2
eigs_L_SU2 = np.linalg.eigvalsh(L_SU2)
eigs_L_SU2_sorted = np.sort(eigs_L_SU2)

print(f"\n  SU(2) Laplacian eigenvalues (first 15):")
for i, e in enumerate(eigs_L_SU2_sorted[:15]):
    print(f"    mu_{i} = {e:.6f}")

print(f"\n  SU(2) spectral gap: {eigs_L_SU2_sorted[1]:.6f}")
print(f"  Full graph spectral gap: {12 - eigs_SU2_sorted[0]:.6f}... no, that's the SU2 eigenvalue.")

# ===============================================================
# PART B2: SU(2) Green's function and loop factors
# ===============================================================
print("\n--- SU(2) Green's function ---")

mu_SU2 = eigs_L_SU2_sorted
mu_SU2_nz = mu_SU2[mu_SU2 > 0.01]

G_reg_SU2 = np.sum(1.0/mu_SU2_nz) / N
print(f"  G_reg_SU2(x,x) = {G_reg_SU2:.6f}")
print(f"  G_reg_full(x,x) = 0.097910")
print(f"  Ratio: {G_reg_SU2/0.097910:.4f}")
print(f"  1/phi = {1/PHI:.6f}")
print(f"  1/phi^2 = {1/PHI**2:.6f}")
print()

# Check if any simple expression matches
for name, val in [("1/phi", 1/PHI), ("1/phi^2", 1/PHI**2), ("1/2", 0.5),
                   ("1/3", 1/3), ("1/4", 0.25), ("1/(2*phi)", 1/(2*PHI)),
                   ("sin2tW", SIN2TW), ("1/a_1", 0.2), ("1/b_1", 1/6.0),
                   ("1/deg", 1/12.0), ("alpha", ALPHA)]:
    err = abs(G_reg_SU2 - val)/val * 100
    if err < 30:
        print(f"  G_reg_SU2 ~ {name} = {val:.6f} (err {err:.1f}%)")

# Loop factors at various masses
print(f"\n  SU(2) loop factor I_1(m^2) = (1/N)*sum 1/(mu+m^2):")
for m2 in [0.1, 0.5, 1.0, PHI, 5.0]:
    I1 = np.mean(1.0/(mu_SU2 + m2))
    for name, val in [("1/phi", 1/PHI), ("1/phi^2", 1/PHI**2), ("1/(4pi)", 1/(4*np.pi))]:
        err = abs(I1 - val)/val * 100
        if err < 15:
            print(f"    m^2={m2:.2f}: I_1 = {I1:.6f} ~ {name} = {val:.6f} ({err:.1f}%)")

# ===============================================================
# PART C: Spectral identities for n=35
# ===============================================================
print("\n" + "="*70)
print("PART C: SPECTRAL IDENTITIES - WHY 35 = h + a_1?")
print("="*70)

# Full graph eigenvalues (adjacency)
eigs_full = np.linalg.eigvalsh(A_full)
eigs_full_sorted = np.sort(eigs_full)[::-1]
mu_full = 12.0 - eigs_full_sorted  # Laplacian eigenvalues

# Identity 1: Spectral zeta function
print("\n  Spectral zeta function zeta_L(s) = sum_{mu>0} mu^(-s):")
mu_nz = mu_full[mu_full > 0.01]
for s in [0.5, 1.0, 1.5, 2.0, 3.0]:
    zeta = np.sum(mu_nz**(-s))
    # Check against clean values
    for name, val in [("deg", 12.0), ("N", 120.0), ("h", 30.0), ("phi^k", None),
                       ("a_1*b_1", 30.0), ("deg-1/4", 11.75)]:
        if val is not None:
            if abs(zeta-val)/val < 0.05:
                print(f"    zeta_L({s}) = {zeta:.6f} ~ {name} = {val} ({abs(zeta-val)/val*100:.3f}%)")
                break
    else:
        print(f"    zeta_L({s}) = {zeta:.6f}")
        # Check powers of phi
        for k in range(-5, 10):
            if abs(zeta - PHI**k)/PHI**k < 0.02:
                print(f"      ~ phi^{k} = {PHI**k:.6f} ({abs(zeta-PHI**k)/PHI**k*100:.2f}%)")

# Identity 2: Trace identities
print("\n  Trace identities of adjacency matrix A:")
for p in range(1, 8):
    tr = np.trace(np.linalg.matrix_power(A_full.astype(int), p)) / N
    # Tr(A^p)/N = sum eigenvalue^p * mult / N
    print(f"    Tr(A^{p})/N = {tr:.4f}", end="")
    if p == 1:
        print(f"  (= 0, traceless)")
    elif p == 2:
        print(f"  (= deg = {deg})")
    elif p == 3:
        # Number of triangles per vertex
        print(f"  (= triangles per vertex * 2)")
    else:
        # Check for phi-related
        for name, val in [("phi^k", None), ("deg^(p/2)", deg**(p/2.0)),
                           ("N", 120.0), ("h*a_1", 150.0)]:
            if val is not None and abs(tr-val)/max(abs(val),1) < 0.1:
                print(f"  ~ {name}")
                break
        else:
            print()

# Identity 3: Can we get 35 from eigenvalue products/sums?
print("\n  Searching for 35 from spectral quantities:")
# Unique eigenvalues
unique_eigs = [12.0, 6*PHI, 4*PHI, 3.0, 0.0, -2.0, 4-4*PHI, -3.0, 6-6*PHI]
unique_mults = [1, 4, 9, 16, 25, 36, 9, 16, 4]

# Adjacency eigenvalues in a+b*phi form
ab_eigs = [(12,0), (0,6), (0,4), (3,0), (0,0), (-2,0), (4,-4), (-3,0), (6,-6)]

# Try sums/products of a,b values
a_vals = [ab[0] for ab in ab_eigs]
b_vals = [ab[1] for ab in ab_eigs]
print(f"  sum(a_i) = {sum(a_vals)}")
print(f"  sum(b_i) = {sum(b_vals)}")
print(f"  sum(|a_i|) = {sum(abs(a) for a in a_vals)}")
print(f"  sum(|b_i|) = {sum(abs(b) for b in b_vals)}")
print(f"  sum(a_i^2) = {sum(a**2 for a in a_vals)}")
print(f"  sum(a_i*mult_i) = {sum(a*m for a,m in zip(a_vals, unique_mults))}")
print(f"  sum(b_i*mult_i) = {sum(b*m for b,m in zip(b_vals, unique_mults))}")
print(f"  sum(n_i) where n=5a+6b = {sum(5*a+6*b for a,b in ab_eigs)}")
print(f"  sum(|n_i|) = {sum(abs(5*a+6*b) for a,b in ab_eigs)}")
print(f"  sum(n_i * mult_i) = {sum((5*a+6*b)*m for (a,b),m in zip(ab_eigs, unique_mults))}")
print(f"  h(E8) = {h_E8}, a_1 = {a_1}")
print(f"  h + a_1 = {h_E8 + a_1}")

# Try: sum of positive n_i values
pos_n = [(5*a+6*b) for a,b in ab_eigs if 5*a+6*b > 0]
print(f"  Positive n_i: {pos_n}")
print(f"  sum(positive n_i) = {sum(pos_n)}")

# n values
n_vals = [5*a+6*b for a,b in ab_eigs]
print(f"  All n_i: {n_vals}")
print(f"  sum(n_i > 0) = {sum(n for n in n_vals if n > 0)} = {sum(n for n in n_vals if n > 0)}")

# ===============================================================
# PART D: Random walk and return probability
# ===============================================================
print("\n" + "="*70)
print("PART D: RANDOM WALK ON 600-CELL")
print("="*70)

# Transition matrix T = A/deg (random walk on regular graph)
T = A_full / deg

# Return probability after n steps: p(n) = [T^n]_{00}
# By vertex-transitivity: p(n) = Tr(T^n)/N = (1/N)*sum (lambda_k/deg)^n * mult_k

print("\n  Return probability p(n) = (1/N)*sum(lambda_k/deg)^n * mult_k:")
print(f"  {'n':>4s} {'p(n)':>14s} {'p*N':>10s} {'p*phi^n':>14s} {'log(p)/log(phi)':>16s}")
print("  " + "-"*60)

for n in [1, 5, 10, 15, 20, 25, 30, 33, 34, 35, 36, 40, 50]:
    pn = sum(m * (lam/deg)**n for lam, m in zip(unique_eigs, unique_mults)) / N
    pn_real = pn.real if isinstance(pn, complex) else pn
    log_ratio = np.log(abs(pn_real))/np.log(PHI) if abs(pn_real) > 1e-30 else -999
    phi_n = pn_real * PHI**n if abs(pn_real) > 1e-30 else 0
    print(f"  {n:>4d} {pn_real:>14.6e} {pn_real*N:>10.4e} {phi_n:>14.6e} {log_ratio:>16.2f}")

# Check: p(n)*phi^n for n ~ 35?
print("\n  p(n) * phi^n near n=35:")
for n in range(30, 41):
    pn = sum(m * (lam/deg)**n for lam, m in zip(unique_eigs, unique_mults)) / N
    pn_real = pn.real if isinstance(pn, complex) else pn
    phi_n = pn_real * PHI**n
    m_pred = m_e * abs(pn_real)
    marker = ""
    if abs(m_pred - m3_exp)/m3_exp < 0.2:
        marker = " <-- NEAR m3!"
    print(f"    n={n}: p*phi^n = {phi_n:.6e}, m_e*p(n) = {m_pred:.4e} eV{marker}")

# ===============================================================
# PART E: Galois-broken DSI (Majorana mechanism)
# ===============================================================
print("\n" + "="*70)
print("PART E: GALOIS-BROKEN DSI (MAJORANA MECHANISM)")
print("="*70)

# For charged fermions: E = V(z) + V(z') = 0 (Galois invariant)
# V(z) = product(z - lambda_k) over eigenvalues
#
# For MAJORANA neutrinos: E = V(z) ONLY (no Galois partner!)
# The Galois symmetry is BROKEN because Majorana = self-conjugate
#
# Then: V(z) = 0 gives z = eigenvalue (9 solutions)
# The LIGHTEST Majorana state corresponds to the eigenvalue
# furthest from zero on the DSI lattice
#
# Which eigenvalue has the LARGEST n = 5a + 6b?
# From the list: n=60, 36, 24, 15, 0, -10, -4, -15, -6
# LARGEST positive: n=60 (eigenvalue 12, mult 1)
# Next: n=36 (eigenvalue 6*phi, mult 4)

print("\n  Eigenvalues ranked by DSI level n = 5a + 6b:")
for (a,b), m, lam in sorted(zip(ab_eigs, unique_mults, unique_eigs),
                              key=lambda x: -(5*x[0][0]+6*x[0][1])):
    n_dsi = 5*a + 6*b
    print(f"    n={n_dsi:>4d}: lambda={lam:>8.4f} (a,b)=({a},{b}) mult={m}")

print(f"\n  Maximum DSI level: n=60 (eigenvalue 12, the constant eigenvector)")
print(f"  This is the 'vacuum' state. Mass ~ m_e/phi^60 = {m_e/PHI**60:.2e} eV")
print(f"  Way too small! n=60 is too deep.")

# But what if we sum n over ALL eigenvalues?
n_sum = sum(5*a+6*b for a,b in ab_eigs)
n_sum_pos = sum(5*a+6*b for a,b in ab_eigs if 5*a+6*b > 0)
n_sum_weighted = sum((5*a+6*b)*m for (a,b),m in zip(ab_eigs, unique_mults))
print(f"\n  sum(n_i) = {n_sum}")
print(f"  sum(positive n_i) = {n_sum_pos}")
print(f"  sum(n_i * mult_i) / N = {n_sum_weighted/N:.4f}")

# INTERESTING: look at n_sum_pos
# n=60+36+24+15 = 135 = h(E8)*a_1 - 15 = 150-15? No, 135 = 9*15 = 27*5
# Hmm. Or 135 = 3*45 = 3*a_1*N_eig
print(f"  135 = 3 * 45 = N_gen * a_1 * N_eig = {N_gen*a_1*N_eig}")
print(f"  135 / 120 = {135/120}")

# ===============================================================
# PART F: The h(E8) + a_1 derivation attempt
# ===============================================================
print("\n" + "="*70)
print("PART F: ATTEMPTING TO DERIVE n_nu = h(E8) + a_1 = 35")
print("="*70)

# ARGUMENT:
# 1. Charged fermion spectrum: n in [0, n_max] where n_max = 5*4+6*1 = 26 (top)
# 2. The spectral range is 0 to n_max = 26
# 3. The FULL spectral range of the graph is 0 to n_top = 60 (eigenvalue 12)
# 4. The NEUTRINO, being the lightest massive particle, sits at:
#    n_nu = n_max + N_eig = 26 + 9 = 35 (??)
#    or: n_nu = n_max + (h - n_max) = h + something
#
# Let's check: n_top - n_nu = 60 - 35 = 25 = 5^2 = a_1^2 = mult of zero eigenvalue!

n_top = 60  # DSI level of the highest eigenvalue (constant eigenvector)
n_nu = 35   # target: h + a_1

print(f"\n  n_top (eigenvalue 12) = {n_top}")
print(f"  n_nu (neutrino)       = {n_nu}")
print(f"  n_top - n_nu          = {n_top - n_nu}")
print(f"  25 = a_1^2 = mult(lambda=0) = {a_1**2}")
print()
print(f"  n_nu = n_top - a_1^2 = 60 - 25 = 35  (**NEW!**)")
print()

# Is 60 = 5*12 = a_1 * deg? YES!
print(f"  n_top = 60 = a_1 * deg = 5 * 12 = {a_1 * deg}")
print(f"  So: n_nu = a_1 * deg - a_1^2 = a_1 * (deg - a_1) = 5 * 7 = 35")
print(f"  Also: n_nu = a_1 * (deg - a_1) = a_1 * (b_1 + 1) = 5 * 7 = 35")
print(f"  And: deg - a_1 = 12 - 5 = 7 = b_1 + 1")
print()

# Check: is this equivalent to h + a_1?
print(f"  h(E8) + a_1 = 30 + 5 = 35")
print(f"  a_1*(deg-a_1) = 5*7 = 35")
print(f"  SAME! Because h = a_1*b_1 = 30 and deg = a_1+b_1+1 = 12")
print(f"  So a_1*(deg-a_1) = a_1*(b_1+1) = a_1*b_1 + a_1 = h + a_1 = 35 QED!")
print()

# NEW INTERPRETATION:
print("  NEW INTERPRETATION:")
print(f"  n_top = (12+0)*a_1 = 60  (eigenvalue 12, DSI level of constant mode)")
print(f"  mult_0 = a_1^2 = 25      (multiplicity of zero eigenvalue)")
print(f"  n_nu = n_top - mult_0 = a_1*deg - a_1^2 = a_1*(deg - a_1)")
print(f"       = a_1*(b_1+1) = h + a_1 = 35")
print()
print("  PHYSICAL PICTURE:")
print("  The neutrino sits a_1^2 = 25 DSI levels BELOW the top of the spectrum.")
print("  25 = multiplicity of the zero eigenvalue = dimension of the 'vacuum' space.")
print("  The neutrino mass is suppressed by phi^25 relative to the maximum:")
print(f"    m_nu_max = m_e/phi^{n_top-25} = m_e/phi^35 (needs factor 2 for Majorana)")

# Check: does the factor 2 come from somewhere?
print(f"\n  Factor 2 = dim(psi_5) = real irrep of 2T (Majorana)")
print(f"  Or: 2 = mult of zero eigenvalue in SU(2) subgraph? Let me check...")

# Zero eigenvalues of SU(2) subgraph
n_zero_SU2 = np.sum(np.abs(eigs_L_SU2) < 0.01)
print(f"  Multiplicity of mu=0 in SU(2) Laplacian: {n_zero_SU2}")

# Connected components of SU(2) subgraph
# n_zero = number of connected components
print(f"  => SU(2) subgraph has {n_zero_SU2} connected component(s)")

# ===============================================================
# PART G: Putting it all together
# ===============================================================
print("\n" + "="*70)
print("PART G: DERIVATION ATTEMPT")
print("="*70)

print("""
PROPOSED DERIVATION of m_nu = 2 * m_e / phi^35:

Step 1: SPECTRAL RANGE
  The 600-cell adjacency matrix has maximum eigenvalue lambda_0 = 12 = deg.
  Its DSI level is n_top = 5*a + 6*b where (a,b) = (12, 0), giving n_top = 60.
  n_top = a_1 * deg = 5 * 12 = 60.

Step 2: VACUUM DEGENERACY
  The zero eigenvalue (lambda_4 = 0) has multiplicity a_1^2 = 25.
  This is the dimension of the "kernel" of the adjacency action.
  These 25 modes are the "vacuum" - they carry no spectral weight.

Step 3: NEUTRINO LEVEL
  The neutrino, being the lightest massive fermion, sits at the DSI level
  obtained by subtracting the vacuum degeneracy from the spectral top:

    n_nu = n_top - mult(0) = a_1*deg - a_1^2 = a_1*(deg - a_1) = 35

  Equivalently: n_nu = a_1*(b_1 + 1) = h(E8) + a_1 = 35.

Step 4: MAJORANA FACTOR
  The factor 2 comes from the real 2D irrep psi_5 of 2T:
  dim(psi_5) = 2 = the Majorana doubling.
  (Majorana neutrino = particle IS its own antiparticle,
   contributing 2 degrees of freedom to the mass.)

Step 5: MASS FORMULA
  m_nu3 = dim(psi_5) * m_e * phi^(-n_nu)
        = 2 * m_e * phi^(-(h(E8) + a_1))
        = 2 * m_e / phi^35
""")

m3_pred = 2 * m_e / PHI**35
m3_err = abs(m3_pred - m3_exp)/m3_exp*100
print(f"  RESULT: m3 = {m3_pred:.5f} eV (exp: {m3_exp:.5f}, err: {m3_err:.3f}%)")
print()

# STATUS assessment
print("  STATUS ASSESSMENT:")
print("  Step 1: DERIVED (n_top = a_1*deg from spectral structure)")
print("  Step 2: DERIVED (mult(0) = a_1^2 = 25)")
print("  Step 3: MOTIVATED but NOT DERIVED")
print("          WHY subtract mult(0) from n_top?")
print("          This is the weakest link. Needs physical justification.")
print("  Step 4: DERIVED (psi_5 = real 2D irrep, from 2I/2T branching)")
print("  Step 5: DERIVED (phi from DSI, m_e as reference)")
print()
print("  OVERALL: MOTIVATED PATTERN (sub-0.1%)")
print("  The formula is no longer 'numerological' - it uses")
print("  derived spectral quantities in a specific combination.")
print("  But Step 3 needs a rigorous argument for WHY n_nu = n_top - mult(0).")
print()

# Alternative check: is n_top - mult(0) = 35 a coincidence?
print("  COINCIDENCE CHECK:")
print(f"  n_top = 60 = a_1*deg (the ONLY way to decompose 12 in DSI: (12,0))")
print(f"  mult(0) = 25 = a_1^2 (from spectral theory)")
print(f"  60 - 25 = 35 = h+a_1 = a_1*(deg-a_1)")
print(f"  This follows from: a_1*deg - a_1^2 = a_1*(deg-a_1) = a_1*(b_1+1)")
print(f"  And h = a_1*b_1, so a_1*(b_1+1) = h + a_1. QED.")
print(f"  The identity is NOT a coincidence - it's algebraic.")
print(f"  The QUESTION is: WHY does the neutrino level = n_top - mult(0)?")
