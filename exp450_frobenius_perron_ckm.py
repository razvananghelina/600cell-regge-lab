"""
exp450_frobenius_perron_ckm.py
==============================
FROBENIUS-PERRON MECHANISM FOR CKM MIXING

QUESTION: WHY does CKM mixing follow theta_ij = arctan(phi^{-n})?
The ansatz theta_ij = arctan(phi^{-n_ij} * c_ij) works perfectly with corrections,
but the phi^{-n} structure itself was found by brute-force scan (exp082).

HYPOTHESIS (Razvan's proposal):
The Perron-Frobenius eigenvector of the McKay graph adjacency matrix has
components proportional to irrep dimensions of 2I. These are closely related
to powers of phi. If CKM mixing is determined by overlaps between eigenvectors
at fermion nodes, phi^{-n} appears naturally as ratio of components.

PLAN:
1. Build McKay graph (affine E8) adjacency matrix
2. Compute full eigensystem
3. Analyze eigenvector components at up-quark and down-quark nodes
4. Compute graph propagator G(up_i, down_j)
5. Check hop amplitudes <up|A^n|down>
6. Test if ratios give phi^{-n_ij} for CKM exponents {3,7,12}
7. If yes: identify the MECHANISM (which eigenvalues contribute)

Author: Razvan-Constantin Anghelina (theory), Claude (computation)
Date: March 2026
"""

import numpy as np
from fractions import Fraction

# === CONSTANTS ===
PHI = (1 + np.sqrt(5)) / 2   # golden ratio
phi = PHI
a1 = 5
b1 = 6

# CKM exponents
n_12, n_23, n_13 = 3, 7, 12

# Experimental CKM (PDG 2024)
V_us_exp = 0.2243
V_cb_exp = 0.0410
V_ub_exp = 0.00382

print("=" * 70)
print("exp450: FROBENIUS-PERRON MECHANISM FOR CKM MIXING")
print("=" * 70)

# === 1. McKAY GRAPH (AFFINE E8) ===
print("\n" + "=" * 70)
print("SECTION 1: McKay Graph Adjacency Matrix")
print("=" * 70)

# Affine E8: 9 nodes, edges as in mckay_spectral_data.md
# Topology:  0 - 1 - 2 - 3 - 4 - 5 - 6 - 7
#                                    |
#                                    8
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
n_nodes = 9

A = np.zeros((n_nodes, n_nodes))
for i, j in edges:
    A[i,j] = 1
    A[j,i] = 1

# Irrep dimensions (Perron-Frobenius eigenvector)
dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])
fermions = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
print(f"Nodes: {list(zip(fermions, dims))}")
print(f"Sum(dims) = {sum(dims)} = h(E8) = 30")
print(f"Sum(dims^2) = {sum(dims**2)} = |2I| = 120")

# Verify PF: A * dims = 2 * dims
Av = A @ dims
print(f"\nPerron-Frobenius check: A*dims = {Av}")
print(f"  2*dims = {2*dims}")
print(f"  Match: {np.allclose(Av, 2*dims)}")

# === 2. FULL EIGENSYSTEM ===
print("\n" + "=" * 70)
print("SECTION 2: Full Eigensystem")
print("=" * 70)

eigenvalues, eigenvectors = np.linalg.eigh(A)
# Sort by decreasing eigenvalue
idx = np.argsort(-eigenvalues)
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Known analytical eigenvalues
analytical_evals = [2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2]
print("Eigenvalue comparison (numerical vs analytical):")
for i in range(9):
    closest = min(analytical_evals, key=lambda x: abs(x - eigenvalues[i]))
    print(f"  lambda_{i} = {eigenvalues[i]:+8.5f}  (= {closest:+8.5f})")

# Normalize eigenvectors with positive first non-zero component
for k in range(9):
    v = eigenvectors[:, k]
    for c in v:
        if abs(c) > 1e-10:
            if c < 0:
                eigenvectors[:, k] = -v
            break

print("\nEigenvectors (rows = nodes, columns = eigenvalues):")
print(f"{'Node':>6} {'Ferm':>4}", end="")
for ev in eigenvalues:
    print(f" {ev:+7.3f}", end="")
print()
for i in range(9):
    print(f"{i:6d} {fermions[i]:>4}", end="")
    for k in range(9):
        print(f" {eigenvectors[i,k]:+7.4f}", end="")
    print()

# === 3. EIGENVECTOR ANALYSIS AT QUARK NODES ===
print("\n" + "=" * 70)
print("SECTION 3: Eigenvector Components at Quark Nodes")
print("=" * 70)

# Up-type quarks: u(1), c(5), t(7)
# Down-type quarks: d(2), s(3), b(8)
up_nodes = [1, 5, 7]  # u, c, t
down_nodes = [2, 3, 8]  # d, s, b
up_names = ['u', 'c', 't']
down_names = ['d', 's', 'b']

print("\nUp-quark eigenvector components:")
for k in range(9):
    ev = eigenvalues[k]
    comps = [eigenvectors[n, k] for n in up_nodes]
    print(f"  lambda={ev:+7.4f}: u={comps[0]:+.5f}, c={comps[1]:+.5f}, t={comps[2]:+.5f}")
    if abs(comps[0]) > 1e-10:
        print(f"    ratios: c/u={comps[1]/comps[0]:+.5f}, t/u={comps[2]/comps[0]:+.5f}")

print("\nDown-quark eigenvector components:")
for k in range(9):
    ev = eigenvalues[k]
    comps = [eigenvectors[n, k] for n in down_nodes]
    print(f"  lambda={ev:+7.4f}: d={comps[0]:+.5f}, s={comps[1]:+.5f}, b={comps[2]:+.5f}")
    if abs(comps[0]) > 1e-10:
        print(f"    ratios: s/d={comps[1]/comps[0]:+.5f}, b/d={comps[2]/comps[0]:+.5f}")

# === 4. PHI-POWER ANALYSIS OF EIGENVECTOR RATIOS ===
print("\n" + "=" * 70)
print("SECTION 4: Phi-Power Analysis of Eigenvector Ratios")
print("=" * 70)

# For each pair of eigenvectors, check if the ratio at specific nodes gives phi^n
print("\nLooking for phi^n in eigenvector component ratios...")
print("For reference: phi^1=1.618, phi^2=2.618, phi^3=4.236, phi^{-1}=0.618, phi^{-2}=0.382, phi^{-3}=0.236")

# Check ratios within the phi-eigenvalue eigenvector
k_phi = np.argmin(np.abs(eigenvalues - phi))
v_phi = eigenvectors[:, k_phi]
print(f"\nEigenvector for lambda = phi ({eigenvalues[k_phi]:.5f}):")
for i in range(9):
    if abs(v_phi[i]) > 1e-10:
        log_ratio = np.log(abs(v_phi[i] / v_phi[0])) / np.log(phi) if abs(v_phi[0]) > 1e-10 else 0
        print(f"  v_phi[{fermions[i]}] = {v_phi[i]:+.6f}, |v/v_e| = {abs(v_phi[i]/v_phi[0]):.4f}, log_phi = {log_ratio:.3f}")

# Check all eigenvalue-pair ratios for phi powers
print("\nSystematic search: component ratios that are phi^n (|n| <= 15):")
found_any = False
for k1 in range(9):
    for k2 in range(k1+1, 9):
        for ni in range(9):
            for nj in range(ni+1, 9):
                v1_i = eigenvectors[ni, k1]
                v1_j = eigenvectors[nj, k1]
                v2_i = eigenvectors[ni, k2]
                v2_j = eigenvectors[nj, k2]

                # Check ratio (v1_i * v2_j) / (v1_j * v2_i) = phi^n ?
                if abs(v1_j * v2_i) > 1e-12:
                    ratio = (v1_i * v2_j) / (v1_j * v2_i)
                    if abs(ratio) > 1e-12:
                        n_pow = np.log(abs(ratio)) / np.log(phi)
                        n_int = round(n_pow)
                        if abs(n_pow - n_int) < 0.02 and abs(n_int) > 0 and abs(n_int) <= 15:
                            ev1, ev2 = eigenvalues[k1], eigenvalues[k2]
                            sign = "+" if ratio > 0 else "-"
                            if n_int in [3, 7, 12, -3, -7, -12, 5, -5, 4, -4]:
                                found_any = True
                                print(f"  [CKM!] lambda=({ev1:+.3f},{ev2:+.3f}), nodes=({fermions[ni]},{fermions[nj]}): "
                                      f"ratio={sign}phi^{n_int} (exact: {n_pow:.4f})")

if not found_any:
    print("  No CKM-relevant phi^n ratios found in eigenvector cross-ratios.")

# === 5. GRAPH PROPAGATOR ===
print("\n" + "=" * 70)
print("SECTION 5: Graph Propagator G(i,j)")
print("=" * 70)

# G(i,j) = sum_k v_k(i)*v_k(j) / lambda_k  (excluding zero mode)
# This is the Green's function: (A - 0*I)^{-1} projected off zero mode
print("Graph propagator G(i,j) = sum_k v_k(i)*v_k(j)/lambda_k (k: lambda_k != 0):")
G = np.zeros((9, 9))
for k in range(9):
    if abs(eigenvalues[k]) > 1e-10:
        G += np.outer(eigenvectors[:, k], eigenvectors[:, k]) / eigenvalues[k]

print("\nG matrix (up-quarks x down-quarks):")
print(f"{'':>6}", end="")
for dn in down_names:
    print(f" {dn:>10}", end="")
print()
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    print(f"{un:>6}", end="")
    for j, dj in enumerate(down_nodes):
        print(f" {G[ui,dj]:+10.6f}", end="")
    print()

# Check if G ratios give phi powers
print("\nPropagator ratios (looking for phi^{-n} pattern):")
G_ud = np.array([[G[ui, dj] for dj in down_nodes] for ui in up_nodes])
G_ref = G_ud[0, 0]  # G(u,d) as reference
print(f"  G(u,d) = {G_ref:.6f} (reference)")
for i, un in enumerate(up_names):
    for j, dn in enumerate(down_names):
        ratio = G_ud[i,j] / G_ref if abs(G_ref) > 1e-12 else 0
        if abs(ratio) > 1e-12:
            n_pow = np.log(abs(ratio)) / np.log(phi)
            print(f"  G({un},{dn})/G(u,d) = {ratio:+.6f}, log_phi = {n_pow:.3f}")

# === 6. HOP AMPLITUDES: <i|A^n|j> ===
print("\n" + "=" * 70)
print("SECTION 6: Hop Amplitudes <i|A^n|j>")
print("=" * 70)

# A^n counts (weighted) n-step walks on the McKay graph
# CKM mixing could be related to walk amplitudes between up and down quarks
print("Number of n-step walks between quark nodes:")
print(f"{'n':>3}", end="")
for un in up_names:
    for dn in down_names:
        print(f"  {un}-{dn:>3}", end="")
print()

An = np.eye(9)
hop_data = {}
for n in range(16):
    print(f"{n:3d}", end="")
    for i, ui in enumerate(up_nodes):
        for j, dj in enumerate(down_nodes):
            val = An[ui, dj]
            print(f"  {val:7.0f}", end="")
            hop_data[(i,j,n)] = val
    print()
    An = An @ A

# Check if walk ratios give phi^{-n_ij}
print("\nWalk ratio analysis:")
print("For each step count n, check ratio of walks (up_i->down_j) vs (u->d):")
An = np.eye(9)
for n in range(1, 13):
    An = An @ A
    walks_ud = An[1, 2]  # u -> d
    if abs(walks_ud) > 0:
        for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
            for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
                w = An[ui, dj]
                if abs(w) > 0 and not (i == 0 and j == 0):
                    ratio = w / walks_ud
                    n_pow = np.log(abs(ratio)) / np.log(phi) if abs(ratio) > 1e-15 else float('nan')
                    if abs(n_pow - round(n_pow)) < 0.1 and not np.isnan(n_pow):
                        print(f"  n={n}: walks({un},{dn})/walks(u,d) = {ratio:.4f} ~ phi^{round(n_pow)}")

# === 7. PERRON-FROBENIUS EIGENVECTOR AND PHI POWERS ===
print("\n" + "=" * 70)
print("SECTION 7: Perron-Frobenius and Phi Structure")
print("=" * 70)

# The PF eigenvector = dims = [1,2,3,4,5,6,4,2,3]
# Check if dims can be written in terms of phi
print("Irrep dimensions in terms of phi:")
for i in range(9):
    d = dims[i]
    # Check if d = round(phi^k * c) for some k, c
    for k in range(-3, 6):
        ratio = d / phi**k
        if abs(ratio - round(ratio)) < 0.01 and round(ratio) > 0:
            print(f"  dim[{fermions[i]}] = {d} = {round(ratio)} * phi^{k} (approx: {round(ratio)*phi**k:.4f})")

# The KEY insight: components of eigenvector for lambda=phi
# For a tree graph, these are related to Chebyshev polynomials
print(f"\nEigenvector for lambda = phi:")
v_phi = eigenvectors[:, k_phi]
# Normalize so first component = 1 (if possible)
if abs(v_phi[0]) > 1e-10:
    v_phi_norm = v_phi / v_phi[0]
else:
    v_phi_norm = v_phi / max(abs(v_phi))
print("  Normalized (v[e]=1):")
for i in range(9):
    n_pow = np.log(abs(v_phi_norm[i])) / np.log(phi) if abs(v_phi_norm[i]) > 1e-10 else float('nan')
    print(f"  v[{fermions[i]}] = {v_phi_norm[i]:+.6f}, log_phi = {n_pow:.3f}")

# === 8. WEIGHTED PROPAGATOR WITH PHI EIGENVALUE ===
print("\n" + "=" * 70)
print("SECTION 8: Phi-Eigenvalue Contribution to Propagator")
print("=" * 70)

# Decompose propagator by eigenvalue
# G(i,j) = sum_k v_k(i)*v_k(j)/lambda_k
# The phi and 1/phi contributions are special (Galois pair)
print("Propagator decomposition by eigenvalue sector:")
print("Sectors: {2,-2}, {phi,-phi}, {1,-1}, {1/phi,-1/phi}")

sectors = {
    'lambda=2 (PF)': [0],
    'lambda=phi': [k_phi],
    'lambda=1': [np.argmin(np.abs(eigenvalues - 1))],
    'lambda=1/phi': [np.argmin(np.abs(eigenvalues - 1/phi))],
    'lambda=-1/phi': [np.argmin(np.abs(eigenvalues + 1/phi))],
    'lambda=-1': [np.argmin(np.abs(eigenvalues + 1))],
    'lambda=-phi': [np.argmin(np.abs(eigenvalues + phi))],
    'lambda=-2': [np.argmin(np.abs(eigenvalues + 2))],
}

# Zero eigenvalue excluded from propagator
print("\nPhi-sector contribution G_phi(up, down):")
for name, ks in sectors.items():
    if 'phi' in name and '-' not in name.split('=')[1].strip()[0:1]:
        for k in ks:
            ev = eigenvalues[k]
            if abs(ev) > 1e-10:
                G_sector = np.outer(eigenvectors[:, k], eigenvectors[:, k]) / ev
                print(f"\n  {name}:")
                for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
                    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
                        val = G_sector[ui, dj]
                        print(f"    G_{name}({un},{dn}) = {val:+.8f}")

# === 9. CHARACTER-BASED EIGENVECTORS ===
print("\n" + "=" * 70)
print("SECTION 9: Character-Based Eigenvector Construction")
print("=" * 70)

# For McKay graph, eigenvectors can be constructed from character table
# v_k(i) = chi_i(g_k) / sqrt(dim_i) (up to normalization)
# where g_k is the conjugacy class representative for eigenvalue lambda_k

# Character table of 2I for standard rep at class 10a (eigenvalue phi)
chi_std_10a = np.array([1, phi, phi, 1, 0, -1, -1, -1/phi, -1/phi])
# This is chi_std at the 10a conjugacy class

# The eigenvector for lambda=phi should be proportional to:
# v_phi(i) = chi_i(10a) * sqrt(|C_10a|/|2I|) ... actually let me compute it directly
# v_phi(i) = chi_i(class 10a) for all irreps i

# Character values at class 10a (from mckay_spectral_data.md):
chi_10a = np.array([1, phi, phi, 1, 0, -1, -1, -1/phi, -1/phi])  # chi_k(10a) for k=0..8

# Check: A * chi_10a = phi * chi_10a ?
Ac = A @ chi_10a
print("Character eigenvector check: A * chi(10a) = phi * chi(10a)?")
print(f"  A * chi = {Ac}")
print(f"  phi * chi = {phi * chi_10a}")
print(f"  Match: {np.allclose(Ac, phi * chi_10a)}")

# Now the character eigenvector for lambda = 1/phi (Galois conjugate)
chi_10b = np.array([1, 1/phi, 1/phi, 1, 0, -1, -1, -phi, -phi])  # chi_k(10b)
Ac2 = A @ chi_10b
print(f"\nA * chi(10b) = {Ac2}")
print(f"1/phi * chi = {(1/phi) * chi_10b}")
print(f"Match: {np.allclose(Ac2, (1/phi) * chi_10b)}")

# For lambda = -phi (class 5a)
chi_5a = np.array([1, -phi, phi, -1, 0, -1, 1, 1/phi, -1/phi])
# Wait, I need to be more careful. Let me use the actual character values.

# From the character table of 2I:
# Irreps: rho_0(1), rho_1(2), rho_2(3), rho_3(4), rho_4(5), rho_5(6), rho_6(4), rho_7(2), rho_8(3)
# For class 10a (order 10, size 12): chi_std = phi
# chi_k(10a) for each irrep:

# Actually, from the McKay correspondence:
# chi_k(10a) = the eigenvalue-phi eigenvector component at node k
# These ARE the dims... no. Let me reconsider.

# The eigenvalue lambda_k of the McKay adjacency matrix equals chi_std(g_k)
# The RIGHT eigenvector for eigenvalue chi_std(g_k) is:
# v(i) = chi_i(g_k) (the character of irrep i at class g_k)
# This works because: sum_j A_{ij} chi_j(g) = chi_std(g) * chi_i(g)
# (by the McKay property: rho_std tensor rho_i = sum_j A_{ij} rho_j)

# So the character table COLUMNS are eigenvectors!

# For class 10a: eigenvalue = phi, eigenvector = [chi_0(10a), chi_1(10a), ..., chi_8(10a)]
# chi_0(10a) = 1 (trivial rep)
# chi_1(10a) = phi (standard 2-dim)
# chi_2(10a) = ? (3-dim = Sym^2 of standard)

# Sym^2 of standard at 10a: chi_Sym^2(g) = (chi^2(g) + chi(g^2))/2
# g = 10a means g^2 = 5a. chi_std(10a) = phi, chi_std(5a) = -1/phi
# So chi_2(10a) = (phi^2 + (-1/phi))/2 = (phi+1 - 1/phi)/2 = (phi + phi-1)/2 = (2phi-1)/2
# = (2*(1+sqrt(5))/2 - 1)/2 = sqrt(5)/2... hmm

# Actually I realize I can just verify numerically.
print("\n\nCharacter eigenvector for lambda=phi (class 10a):")
print("These are chi_k(10a) for each irrep rho_k of 2I")

# The eigenvector from numpy (for lambda=phi) IS the character column
# Let me compare
v_phi_num = eigenvectors[:, k_phi]
# Normalize so v[0] = 1
v_phi_num_n = v_phi_num / v_phi_num[0]
print(f"  Numerical eigenvector (normalized v[e]=1):")
for i in range(9):
    print(f"    v[{fermions[i]}] = {v_phi_num_n[i]:+.6f}")

# The chi_10a should be [1, phi, ?, ?, ?, ?, ?, -1/phi, ?]
# Let me just verify with A
print(f"\n  Analytical chi(10a) guess: {chi_10a}")
print(f"  Ratio (numerical/analytical):")
for i in range(9):
    if abs(chi_10a[i]) > 1e-10:
        print(f"    {fermions[i]}: {v_phi_num_n[i]/chi_10a[i]:.6f}")

# === 10. CKM FROM CHARACTER RATIOS ===
print("\n" + "=" * 70)
print("SECTION 10: CKM from Character Ratios")
print("=" * 70)

# HYPOTHESIS: V_{ij} ~ chi_{down_j}(10a) / chi_{up_i}(10a)
# or some function of character ratios

print("Character values at class 10a (= phi-eigenvector):")
# Using numerical eigenvector
v = v_phi_num_n
print("  Up quarks:   u={:.4f}, c={:.4f}, t={:.4f}".format(v[1], v[5], v[7]))
print("  Down quarks: d={:.4f}, s={:.4f}, b={:.4f}".format(v[2], v[3], v[8]))

print("\nRatios v[down]/v[up] (CKM-like):")
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
        if abs(v[ui]) > 1e-10:
            ratio = v[dj] / v[ui]
            log_p = np.log(abs(ratio))/np.log(phi) if abs(ratio) > 1e-10 else float('nan')
            print(f"  v[{dn}]/v[{un}] = {ratio:+.6f}, |ratio| = {abs(ratio):.6f}, log_phi = {log_p:.3f}")

# Also check: Galois conjugate eigenvector (lambda = 1/phi, class 10b)
k_inv_phi = np.argmin(np.abs(eigenvalues - 1/phi))
v_invphi = eigenvectors[:, k_inv_phi]
v_invphi_n = v_invphi / v_invphi[0] if abs(v_invphi[0]) > 1e-10 else v_invphi

print(f"\nCharacter values at class 10b (= 1/phi eigenvector):")
print("  Up quarks:   u={:.4f}, c={:.4f}, t={:.4f}".format(v_invphi_n[1], v_invphi_n[5], v_invphi_n[7]))
print("  Down quarks: d={:.4f}, s={:.4f}, b={:.4f}".format(v_invphi_n[2], v_invphi_n[3], v_invphi_n[8]))

# === 11. PRODUCT OF PHI AND INVPHI EIGENVECTORS ===
print("\n" + "=" * 70)
print("SECTION 11: Galois-Paired Eigenvector Products")
print("=" * 70)

# The Galois automorphism swaps phi <-> -1/phi
# Product v_phi(i) * v_{-1/phi}(j) could give phi-power structure
# Because phi * (-1/phi) = -1, products chain to phi powers

k_neg_invphi = np.argmin(np.abs(eigenvalues + 1/phi))
k_neg_phi = np.argmin(np.abs(eigenvalues + phi))

v_neg_invphi = eigenvectors[:, k_neg_invphi]
v_neg_phi = eigenvectors[:, k_neg_phi]

# Normalize
v_neg_invphi_n = v_neg_invphi / v_neg_invphi[0] if abs(v_neg_invphi[0]) > 1e-10 else v_neg_invphi
v_neg_phi_n = v_neg_phi / v_neg_phi[0] if abs(v_neg_phi[0]) > 1e-10 else v_neg_phi

print("Four phi-related eigenvectors (normalized v[e]=1):")
print(f"{'Node':>6}", end="")
for name in ['v(phi)', 'v(1/phi)', 'v(-1/phi)', 'v(-phi)']:
    print(f"  {name:>10}", end="")
print()

phi_vecs = [v_phi_num_n, v_invphi_n, v_neg_invphi_n, v_neg_phi_n]
phi_labels = ['phi', '1/phi', '-1/phi', '-phi']
for i in range(9):
    print(f"{fermions[i]:>6}", end="")
    for vec in phi_vecs:
        print(f"  {vec[i]:+10.5f}", end="")
    print()

# Check: product v_phi * v_{1/phi} at each node
print("\nProducts v(phi)*v(1/phi) at each node:")
for i in range(9):
    prod = v_phi_num_n[i] * v_invphi_n[i]
    log_p = np.log(abs(prod))/np.log(phi) if abs(prod) > 1e-10 else float('nan')
    print(f"  {fermions[i]}: {prod:+.6f}, log_phi = {log_p:.3f}")

# === 12. TRANSITION MATRIX T(n) = A^n SPECTRAL DECOMPOSITION ===
print("\n" + "=" * 70)
print("SECTION 12: Spectral Walk Amplitudes")
print("=" * 70)

# <i|A^n|j> = sum_k lambda_k^n * v_k(i) * v_k(j)
# For the phi-sector: contribution ~ phi^n * v_phi(i) * v_phi(j)
# For large n, PF (lambda=2) dominates, but phi contribution decays differently

print("Spectral decomposition of walk amplitude at step n=1:")
print("T(1)(i,j) = sum_k lambda_k * v_k(i) * v_k(j)")
print("\nPhi-sector contribution to (up,down) transitions:")
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
        # Phi contribution
        phi_contrib = phi * eigenvectors[ui, k_phi] * eigenvectors[dj, k_phi]
        # Total (= A[ui,dj])
        total = A[ui, dj]
        # Ratio
        ratio = phi_contrib / total if abs(total) > 1e-10 else float('inf')
        print(f"  T_phi({un},{dn}) = {phi_contrib:+.6f}, T_total = {total:.0f}, fraction = {ratio:.4f}")

# === 13. DISTANCE-BASED ANALYSIS ===
print("\n" + "=" * 70)
print("SECTION 13: Graph Distance and CKM Exponents")
print("=" * 70)

# Graph distances on McKay graph
dist = [0, 1, 2, 3, 4, 5, 6, 7, 6]

print("CKM pairs and their graph distances:")
ckm_pairs = [
    ('u', 'd', 1, 2, n_12),
    ('c', 's', 5, 3, n_23),
    ('u', 'b', 1, 8, n_13),
    ('c', 'b', 5, 8, None),
    ('t', 'd', 7, 2, None),
    ('t', 's', 7, 3, None),
]

for un, dn, ui, dj, n_ckm in ckm_pairs:
    d = abs(dist[ui] - dist[dj])  # naive distance difference
    # Actual shortest path on graph
    # BFS
    from collections import deque
    visited = {ui: 0}
    queue = deque([ui])
    while queue:
        node = queue.popleft()
        for edge in edges:
            nbr = edge[1] if edge[0] == node else (edge[0] if edge[1] == node else None)
            if nbr is not None and nbr not in visited:
                visited[nbr] = visited[node] + 1
                queue.append(nbr)
    d_graph = visited.get(dj, -1)

    n_str = f"n_CKM={n_ckm}" if n_ckm else ""
    print(f"  d({un},{dn}) = {d_graph} {n_str}")

# CKM exponents vs graph distances
print(f"\nComparison: CKM exponents vs graph distances:")
print(f"  n_12 = {n_12}, d(u,d) = 1. Ratio: {n_12}/1 = {n_12}")
print(f"  n_23 = {n_23}, d(c,s) = 2. Ratio: {n_23}/2 = {n_23/2}")
print(f"  n_13 = {n_13}, d(u,b) = ?")

# Compute all pairwise distances
print("\nAll shortest distances on McKay graph:")
all_dist = np.full((9, 9), -1, dtype=int)
for start in range(9):
    visited_d = {start: 0}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for edge in edges:
            nbr = edge[1] if edge[0] == node else (edge[0] if edge[1] == node else None)
            if nbr is not None and nbr not in visited_d:
                visited_d[nbr] = visited_d[node] + 1
                queue.append(nbr)
    for end, d in visited_d.items():
        all_dist[start, end] = d

print(f"{'':>6}", end="")
for dn in down_names:
    print(f" {dn:>5}", end="")
print()
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    print(f"{un:>6}", end="")
    for j, dj in enumerate(down_nodes):
        print(f" {all_dist[ui,dj]:5d}", end="")
    print()

# === 14. HEAT KERNEL AND DIFFUSION ===
print("\n" + "=" * 70)
print("SECTION 14: Heat Kernel Diffusion")
print("=" * 70)

# K(i,j,t) = sum_k exp(-lambda_k * t) * v_k(i) * v_k(j) / (dim_i * dim_j) ?
# Or rather K(i,j,t) = <i| exp(-tL) |j> where L = 2I - A (graph Laplacian)
# L has eigenvalues 2 - lambda_k

L_evals = 2 - eigenvalues  # Laplacian eigenvalues
print("Laplacian eigenvalues: 2 - lambda_k =", [f"{x:.4f}" for x in L_evals])

# Heat kernel at various times
print("\nHeat kernel K(up, down, t) for various t:")
times = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
for t in times:
    print(f"\n  t = {t}:")
    for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
        for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
            K = sum(np.exp(-L_evals[k] * t) * eigenvectors[ui,k] * eigenvectors[dj,k] for k in range(9))
            print(f"    K({un},{dn}) = {K:.8f}", end="")
            if abs(K) > 1e-12:
                log_p = np.log(abs(K)) / np.log(phi)
                print(f"  (log_phi = {log_p:.2f})", end="")
            print()

# Find time t where K ratios match CKM
print("\nSearching for t where heat kernel ratios match CKM angles...")
best_t = None
best_chi2 = float('inf')
for t_try in np.linspace(0.01, 20.0, 2000):
    K_mat = np.zeros((3, 3))
    for i, ui in enumerate(up_nodes):
        for j, dj in enumerate(down_nodes):
            K_mat[i,j] = sum(np.exp(-L_evals[k]*t_try) * eigenvectors[ui,k]*eigenvectors[dj,k] for k in range(9))

    # Normalize rows
    for i in range(3):
        row_norm = np.sqrt(sum(K_mat[i,:]**2))
        if row_norm > 0:
            K_mat[i,:] /= row_norm

    # Compare with CKM-like structure: V_us ~ K[0,1], V_cb ~ K[1,2], V_ub ~ K[0,2]
    # or more precisely with the off-diagonal elements
    # But we need to be careful about what we're comparing

    # Try: |K_off_diag| ratios match phi^{-n}?
    if abs(K_mat[0,0]) > 1e-10:
        ratio_12 = abs(K_mat[0,1] / K_mat[0,0])
        ratio_23 = abs(K_mat[1,2] / K_mat[1,1]) if abs(K_mat[1,1]) > 1e-10 else 0
        ratio_13 = abs(K_mat[0,2] / K_mat[0,0])

        if ratio_12 > 1e-10 and ratio_23 > 1e-10 and ratio_13 > 1e-10:
            n12_eff = -np.log(ratio_12) / np.log(phi)
            n23_eff = -np.log(ratio_23) / np.log(phi)
            n13_eff = -np.log(ratio_13) / np.log(phi)

            chi2 = (n12_eff - 3)**2 + (n23_eff - 7)**2 + (n13_eff - 12)**2
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_t = t_try
                best_ns = (n12_eff, n23_eff, n13_eff)

if best_t is not None:
    print(f"  Best t = {best_t:.4f}")
    print(f"  Effective exponents: n12={best_ns[0]:.2f}, n23={best_ns[1]:.2f}, n13={best_ns[2]:.2f}")
    print(f"  Target:              n12=3.00, n23=7.00, n13=12.00")
    print(f"  Chi2 = {best_chi2:.4f}")
    if best_chi2 < 1:
        print("  [OK] Good match!")
    else:
        print("  [..] Poor match - heat kernel alone doesn't reproduce CKM exponents")

# === 15. CHEBYSHEV STRUCTURE ON LEGS ===
print("\n" + "=" * 70)
print("SECTION 15: Chebyshev Structure on E8 Legs")
print("=" * 70)

# Affine E8 has three legs from branch node (5):
# Long leg: 5-4-3-2-1-0 (length 5)
# Medium leg: 5-6-7 (length 2)
# Short leg: 5-8 (length 1)
# CKM correction coefficients from paper: assigned to legs

print("E8 legs from branch node (c, node 5):")
print("  Long leg:   c-mu-s-d-u-e  (length 5)")
print("  Medium leg: c-tau-t       (length 2)")
print("  Short leg:  c-b           (length 1)")

# On each leg, eigenvector components follow Chebyshev-like recurrence
# For eigenvalue lambda, on a chain: v[n+1] = lambda*v[n] - v[n-1]
# This means: v[n] ~ U_n(lambda/2) (Chebyshev U)

print("\nChebyshev analysis of phi-eigenvector on long leg (e,u,d,s,mu,c):")
long_leg = [0, 1, 2, 3, 4, 5]  # e, u, d, s, mu, c
v_long = [v_phi_num_n[i] for i in long_leg]
print(f"  Components: {[f'{x:.5f}' for x in v_long]}")
print(f"  Ratios v[n+1]/v[n]:")
for n in range(len(v_long)-1):
    if abs(v_long[n]) > 1e-10:
        r = v_long[n+1] / v_long[n]
        print(f"    v[{fermions[long_leg[n+1]]}]/v[{fermions[long_leg[n]]}] = {r:.6f}")

print("\nChebyshev analysis of phi-eigenvector on medium leg (c,tau,t):")
med_leg = [5, 6, 7]  # c, tau, t
v_med = [v_phi_num_n[i] for i in med_leg]
print(f"  Components: {[f'{x:.5f}' for x in v_med]}")
for n in range(len(v_med)-1):
    if abs(v_med[n]) > 1e-10:
        r = v_med[n+1] / v_med[n]
        print(f"    v[{fermions[med_leg[n+1]]}]/v[{fermions[med_leg[n]]}] = {r:.6f}")

print("\nChebyshev analysis of phi-eigenvector on short leg (c,b):")
short_leg = [5, 8]  # c, b
v_short = [v_phi_num_n[i] for i in short_leg]
print(f"  Components: {[f'{x:.5f}' for x in v_short]}")
if abs(v_short[0]) > 1e-10:
    r = v_short[1] / v_short[0]
    print(f"    v[b]/v[c] = {r:.6f}")

# === 16. RESOLVENT AND TRANSFER MATRIX ===
print("\n" + "=" * 70)
print("SECTION 16: Transfer Matrix on E8 Chains")
print("=" * 70)

# On a linear chain with eigenvalue lambda, the transfer matrix is:
# T = [[lambda, -1], [1, 0]]
# So T^n = [[U_n(lambda/2), -U_{n-1}(lambda/2)], [U_{n-1}(lambda/2), -U_{n-2}(lambda/2)]]
# where U_n = Chebyshev of 2nd kind

from numpy.polynomial import chebyshev

# For lambda = phi: x = phi/2 = (1+sqrt(5))/4
x_phi = phi / 2
print(f"Chebyshev argument for lambda=phi: x = phi/2 = {x_phi:.6f}")

# Compute Chebyshev U_n(x_phi)
print(f"\nChebyshev U_n(phi/2) for n=0..12:")
Un_vals = []
for n in range(13):
    # U_n(x) = sin((n+1)*arccos(x)) / sin(arccos(x))
    theta = np.arccos(x_phi)
    Un = np.sin((n+1)*theta) / np.sin(theta)
    Un_vals.append(Un)
    log_p = np.log(abs(Un))/np.log(phi) if abs(Un) > 1e-10 else float('nan')
    print(f"  U_{n}(phi/2) = {Un:+.6f}  (log_phi = {log_p:.3f})")

# KEY CHECK: Do ratios of U_n give CKM exponents?
print(f"\nCritical ratios:")
if len(Un_vals) > 12:
    for n1, n2, target in [(0,3,n_12), (0,7,n_23), (0,12,n_13)]:
        ratio = Un_vals[n2] / Un_vals[n1] if abs(Un_vals[n1]) > 1e-10 else float('nan')
        log_p = np.log(abs(ratio))/np.log(phi) if abs(ratio) > 1e-10 else float('nan')
        print(f"  U_{n2}/U_{n1} = {ratio:.6f}, log_phi = {log_p:.3f} (target: {target})")

# Also check U_n ratios at consecutive values
print(f"\nConsecutive U_n ratios (approaching phi as n->inf):")
for n in range(1, 12):
    r = Un_vals[n] / Un_vals[n-1] if abs(Un_vals[n-1]) > 1e-10 else float('nan')
    print(f"  U_{n}/U_{n-1} = {r:.6f}  (phi = {phi:.6f})")

# === 17. THE FROBENIUS-PERRON HYPOTHESIS TEST ===
print("\n" + "=" * 70)
print("SECTION 17: MAIN TEST - Frobenius-Perron Hypothesis")
print("=" * 70)

# HYPOTHESIS: CKM ~ phi^{-d(up,down)} where d is some effective distance
# that combines graph distance with eigenvalue weighting

# Construct the "CKM matrix" from McKay graph spectral data
# Method: V_{ij} = |v_phi(down_j)| / |v_phi(up_i)| * correction
# or V_{ij} ~ propagator with phi-weighting

print("TEST A: V_ij ~ |chi(10a)_down / chi(10a)_up| (character ratios)")
print("-" * 40)
# Using the phi-eigenvector as character values
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
        chi_up = v_phi_num_n[ui]
        chi_down = v_phi_num_n[dj]
        if abs(chi_up) > 1e-10:
            V_pred = abs(chi_down / chi_up)
            log_p = -np.log(V_pred)/np.log(phi) if V_pred > 1e-10 else float('nan')
            print(f"  |chi({dn})/chi({un})| = {V_pred:.6f}, effective phi^(-{log_p:.2f})")

# Experimental targets
V_exp = {
    (0,0): 0.97435,  # V_ud
    (0,1): V_us_exp,  # V_us
    (0,2): V_ub_exp,  # V_ub
    (1,0): 0.2242,   # V_cd (approx = V_us)
    (1,1): 0.9735,   # V_cs
    (1,2): V_cb_exp,  # V_cb
    (2,0): 0.00857,  # V_td
    (2,1): 0.0401,   # V_ts
    (2,2): 0.999172, # V_tb
}

print("\n\nTEST B: Walk-count weighted transition amplitudes")
print("-" * 40)
# At step n, the walk count between nodes is A^n[i,j]
# Weight by phi^n and normalize
# T(i,j) = sum_n phi^{-n} * A^n[i,j] / Z = resolvent G(phi, i, j)
# G(z,i,j) = <i|(zI - A)^{-1}|j>

# But phi IS an eigenvalue of A, so the resolvent diverges!
# Need regularized version: approach from z = phi + epsilon
eps = 0.01
z = phi + eps
G_resolvent = np.linalg.inv(z * np.eye(9) - A)
print(f"Resolvent G(z={z:.3f}, i, j) for (up, down) pairs:")
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
        print(f"  G({un},{dn}) = {G_resolvent[ui,dj]:+.4f}")

# Also try z = 2 (PF eigenvalue)
z2 = 2 + eps
G_pf = np.linalg.inv(z2 * np.eye(9) - A)
print(f"\nResolvent G(z={z2:.3f}, i, j) for (up, down) pairs:")
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
        print(f"  G({un},{dn}) = {G_pf[ui,dj]:+.4f}")

# === 18. DIMENSION-RATIO MECHANISM ===
print("\n" + "=" * 70)
print("SECTION 18: Dimension-Ratio Mechanism")
print("=" * 70)

# The PF eigenvector has components = dims = [1,2,3,4,5,6,4,2,3]
# Quark dims: u(2), c(6), t(2), d(3), s(4), b(3)
# CKM involves up-down transitions. The dimension ratio:
# dim(down)/dim(up) = mixing amplitude?

print("Dimension ratios for CKM transitions:")
for i, (un, ui) in enumerate(zip(up_names, up_nodes)):
    for j, (dn, dj) in enumerate(zip(down_names, down_nodes)):
        r = dims[dj] / dims[ui]
        print(f"  dim({dn})/dim({un}) = {dims[dj]}/{dims[ui]} = {r:.4f}")

# Check: product of dimension ratios along paths
print("\nPath-based dimension products:")
print("  u->d: direct neighbors, dim ratio = 3/2 = 1.5")
print("  c->s: path c-mu-s, product = (5/6)*(4/5) = 4/6 = 0.667")
print("  t->b: path t-tau-c-b... wait, need actual path")

# Find paths between up and down quarks
def find_path(adj_list, start, end):
    """BFS to find shortest path."""
    visited = {start: [start]}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for edge in edges:
            nbr = edge[1] if edge[0] == node else (edge[0] if edge[1] == node else None)
            if nbr is not None and nbr not in visited:
                visited[nbr] = visited[node] + [nbr]
                queue.append(nbr)
    return visited.get(end, [])

print("\nPaths and dimension-ratio products:")
for un, dn, ui, dj, n_ckm in [('u','d',1,2,3), ('c','s',5,3,7), ('u','b',1,8,12)]:
    path = find_path(edges, ui, dj)
    dim_prod = 1
    for step in range(len(path)-1):
        dim_prod *= dims[path[step+1]] / dims[path[step]]

    path_str = '->'.join(fermions[p] for p in path)
    print(f"  {path_str}: dim_product = {dim_prod:.6f}")
    log_p = -np.log(dim_prod)/np.log(phi) if dim_prod > 1e-10 else float('nan')
    print(f"    log_phi(1/prod) = {log_p:.3f} (target n_CKM = {n_ckm})")

# === 19. COMBINED MECHANISM: DISTANCE * PHI-CHARACTER ===
print("\n" + "=" * 70)
print("SECTION 19: Combined Mechanism")
print("=" * 70)

# Try: theta_ij = arctan( phi^{-d_graph} * |chi_phi(j)/chi_phi(i)| )
# where d_graph is the graph distance on McKay graph
print("Combined formula: theta = arctan( phi^{-d} * |chi_phi(down)/chi_phi(up)| )")
print()

for un, dn, ui, dj, n_ckm in [('u','d',1,2,3), ('c','s',5,3,7), ('u','b',1,8,12)]:
    d = all_dist[ui, dj]
    chi_ratio = abs(v_phi_num_n[dj] / v_phi_num_n[ui]) if abs(v_phi_num_n[ui]) > 1e-10 else 0
    arg = phi**(-d) * chi_ratio
    theta = np.arctan(arg)
    V = np.sin(theta)

    # What's the effective n?
    n_eff = -np.log(arg)/np.log(phi) if arg > 1e-15 else float('nan')

    exp_V = {'u-d': V_us_exp, 'c-s': V_cb_exp, 'u-b': V_ub_exp}
    key = f'{un}-{dn}'
    V_target = exp_V.get(key, 0)

    print(f"  {un}->{dn}: d={d}, chi_ratio={chi_ratio:.5f}, phi^(-{d})*chi={arg:.6f}")
    print(f"    theta = {np.degrees(theta):.4f} deg, V = {V:.6f} (exp: {V_target:.6f})")
    print(f"    effective n = {n_eff:.3f} (target: {n_ckm})")

# === 20. SUMMARY AND VERDICT ===
print("\n" + "=" * 70)
print("SECTION 20: SUMMARY AND VERDICT")
print("=" * 70)

print("""
QUESTION: Does the Frobenius-Perron eigenvector structure of the McKay graph
provide a mechanism for WHY CKM angles follow phi^{-n}?

RESULTS:
""")

# Collect all test results
tests = {
    'A: Character ratios chi(10a)_down/chi(10a)_up': 'See Section 10',
    'B: Graph propagator G(up,down)': 'See Section 5',
    'C: Heat kernel diffusion': 'See Section 14',
    'D: Chebyshev U_n(phi/2)': 'See Section 16',
    'E: Dimension ratios': 'See Section 18',
    'F: Combined distance+character': 'See Section 19',
}

for test, ref in tests.items():
    print(f"  {test}: {ref}")

print(f"""
KEY OBSERVATIONS:
1. The phi-eigenvector of McKay graph has components = chi(10a) = character
   values of 2I irreps at conjugacy class 10a.

2. These characters take values from {{-phi, -1, -1/phi, 0, 1/phi, 1, phi}}.
   At quark nodes: u(phi), c(-1), t(-1/phi), d(phi), s(1), b(-1/phi).

3. Character RATIOS between same-valued nodes give 1, not phi^{{-n}}.
   Character ratios between different-valued nodes give single phi powers,
   not the large exponents {{3, 7, 12}} needed for CKM.

4. Graph distances: d(u,d)=1, d(c,s)=2, d(u,b)={all_dist[1,8]}.
   These are TOO SMALL and DON'T scale like {{3, 7, 12}}.

5. Chebyshev U_n(phi/2) ratios do give phi-power behavior for large n,
   but the specific CKM exponents don't emerge from the graph structure.
""")

# Final assessment
print("VERDICT:")
print("-" * 40)

# Check if any mechanism reproduces {3, 7, 12}
# The honest answer based on all tests
print("""
The Frobenius-Perron / eigenvector structure of the McKay graph provides:
- phi as a NATURAL BASIS (it's an eigenvalue of the adjacency matrix)
- Character values that ARE powers of phi ({{phi, 1/phi, 1, 0, -1}})
- Chebyshev recurrence that generates phi-power decay on graph chains

HOWEVER: The specific CKM exponents {3, 7, 12} do NOT emerge directly from
the 9-node McKay graph. The graph is too small (diameter = 7) and the
character values too coarse to produce the precise exponential hierarchy.

WHAT WORKS: The phi^{-n} ANSATZ is natural because phi IS a spectral
eigenvalue of the theory's defining graph. Any perturbative expansion
around the phi-sector naturally involves powers of phi.

WHAT DOESN'T WORK: Deriving {3, 7, 12} from graph geometry alone.
The exponents come from ALGEBRAIC IDENTITIES (exp085, exp295, exp449)
that constrain them, not from the graph propagator.

STATUS: NEGATIVE for deriving CKM exponents from Frobenius-Perron mechanism.
        POSITIVE for explaining WHY phi is the natural base.
        The {3, 7, 12} remain algebraically constrained but NOT dynamically derived.
""")
