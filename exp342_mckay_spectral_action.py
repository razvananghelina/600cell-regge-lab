"""
exp342: Spectral Action on the McKay Graph
============================================
Key question: Can D_F on the McKay graph reproduce the mass formula n = 5a + 6b?

From aaaa.txt insight: the adjacency spectrum of the McKay graph is claimed to be
{-2, -phi, -1, -1/phi, 0, +1/phi, +1, +phi, +2} -- powers of the golden ratio!

Investigates:
1. Build affine E8 graph (standard McKay of 2I), verify spectrum
2. Build A5 tensor graph (McKay with 3-dim irrep), its spectrum
3. Define weighted D_F operators
4. Compute propagator from electron node to each fermion
5. Check if propagator reproduces mass exponents 5a + 6b
6. Galois automorphism of spectrum
7. Heat kernel analysis
8. Connection to mass corrections
"""
import numpy as np
from numpy.linalg import eigh, inv, eigvalsh

PHI = (1 + np.sqrt(5)) / 2
PHIp = (1 - np.sqrt(5)) / 2  # = -1/phi
A1 = 5
B1 = 6
N_GEN = 3
ALPHA = 1/137.035999084
ALPHA_S = 1/(2*PHI**3)
C_COEFF = 4/(A1**2 + 1)  # = 2/13
m_e = 0.51099895  # MeV

# PDG masses
pdg = {
    'e': 0.51099895, 'mu': 105.6583755, 'tau': 1776.86,
    'u': 2.16, 'c': 1270.0, 't': 172760.0,
    'd': 4.67, 's': 93.4, 'b': 4180.0,
}
# Quantum numbers
a_vals = {'e':0, 'mu':1, 'tau':1, 'u':3, 'c':2, 't':4, 'd':1, 's':1, 'b':-1}
b_vals = {'e':0, 'mu':1, 'tau':2, 'u':-2, 'c':1, 't':1, 'd':0, 's':1, 'b':4}
n_bare = {f: A1*a_vals[f] + B1*b_vals[f] for f in pdg}
n_eff = {f: np.log(pdg[f]/m_e)/np.log(PHI) for f in pdg}

# ================================================================
# PART 1: Affine E8 Dynkin Diagram (standard McKay graph of 2I)
# ================================================================
print("="*72)
print("PART 1: Affine E8 Adjacency Matrix - Spectrum")
print("="*72)

# Affine E8 diagram: 1-2-3-4-5-6-7-8 with branch 6-9
# Node indices 0-8, dimensions of 2I irreps
dims_2I = [1, 2, 3, 4, 5, 6, 4, 2, 3]
labels_2I = ['rho1', 'rho2', 'rho3', 'rho4', 'rho5', 'rho6', 'rho7', 'rho8', 'rho9']

# Adjacency matrix
A_E8 = np.zeros((9, 9))
edges_E8 = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
for i, j in edges_E8:
    A_E8[i,j] = 1
    A_E8[j,i] = 1

print(f"\nAffine E8 adjacency matrix (9x9):")
print(f"  Nodes: {labels_2I}")
print(f"  Dims:  {dims_2I}")
print(f"  Edges: {edges_E8}")

# Eigenvalues
evals_E8 = np.sort(eigvalsh(A_E8))[::-1]
print(f"\nEigenvalues: {[f'{v:.6f}' for v in evals_E8]}")

# Check against golden ratio values
golden_vals = sorted([0, 1/PHI, -1/PHI, 1, -1, PHI, -PHI, 2, -2], reverse=True)
print(f"\nGolden target: {[f'{v:.6f}' for v in golden_vals]}")
print(f"\nComparison:")
for i, (ev, gv) in enumerate(zip(evals_E8, golden_vals)):
    match = "YES" if abs(ev - gv) < 1e-10 else f"NO (diff={ev-gv:.6f})"
    print(f"  lambda_{i}: {ev:+.10f} vs {gv:+.10f}  {match}")

is_golden = all(abs(ev - gv) < 1e-10 for ev, gv in zip(evals_E8, golden_vals))
print(f"\n  SPECTRUM IS GOLDEN RATIO POWERS: {'YES!!!' if is_golden else 'NO'}")

# ================================================================
# PART 2: Algebraic structure of the spectrum
# ================================================================
print(f"\n{'='*72}")
print("PART 2: Algebraic Structure of Affine E8 Spectrum")
print("="*72)

if is_golden:
    print(f"\nSpectrum = {{0, +-1/phi, +-1, +-phi, +-2}}")
    print(f"  = {{0, +-phi^(-1), +-phi^0, +-phi^1, +-2}}")
    print(f"  Note: 2 is NOT a power of phi. But 2 = phi + 1/phi = phi + phi^(-1)")
    print(f"  Verification: phi + 1/phi = {PHI + 1/PHI:.6f} (should be sqrt(5)+... )")
    print(f"  Actually: phi + 1/phi = (1+sqrt(5))/2 + (sqrt(5)-1)/2 = sqrt(5) = {np.sqrt(5):.6f}")
    print(f"  So 2 != phi + 1/phi. But 2 = phi^2 - phi + 1 ... no.")
    print(f"  2 = 2*1 = 2*phi^0. It's just the degree (max eigenvalue of tree).")
    print(f"\n  Galois automorphism phi -> -1/phi maps:")
    for v in [2, PHI, 1, 1/PHI, 0]:
        # Under phi -> -1/phi:
        # phi -> -1/phi (i.e., phi -> phi' = -1/phi)
        # 1/phi -> -phi
        # 1 -> 1
        # 2 -> 2 (rational)
        # 0 -> 0
        if abs(v - PHI) < 1e-10:
            gv = PHIp  # = -1/phi
        elif abs(v - 1/PHI) < 1e-10:
            gv = -PHI
        elif abs(v - 1) < 1e-10:
            gv = 1
        elif abs(v - 2) < 1e-10:
            gv = 2
        elif abs(v) < 1e-10:
            gv = 0
        else:
            gv = v
        print(f"    {v:+.4f} -> {gv:+.4f}")

    print(f"\n  So Galois maps: {{0,+-1/phi,+-1,+-phi,+-2}} -> {{0,+-phi,+-1,+-1/phi,+-2}}")
    print(f"  The Galois conjugation PERMUTES the eigenvalues!")
    print(f"  It swaps phi <-> -1/phi and 1/phi <-> -phi")
    print(f"  The rational eigenvalues {{0, +-1, +-2}} are FIXED.")

    # Connection to 2*cos(k*pi/5)
    print(f"\n  Connection to pentagonal angles:")
    for k in range(6):
        v = 2*np.cos(k*np.pi/5)
        print(f"    2*cos({k}*pi/5) = 2*cos({k*36}deg) = {v:.6f}")
    print(f"\n  So eigenvalues = {{2*cos(k*pi/5) : k=0,1,2,3,4}} union {{0}}")
    print(f"  with multiplicities: 2 has mult 1, phi and 1/phi have mult 1 each,")
    print(f"  -phi and -1/phi have mult 1, -2 has mult 1, 0 has mult 1, +-1 have mult 1.")

# ================================================================
# PART 3: Eigenvectors and node localization
# ================================================================
print(f"\n{'='*72}")
print("PART 3: Eigenvectors - Which Nodes Do They Localize On?")
print("="*72)

evals_E8_full, evecs_E8 = eigh(A_E8)
# Sort by eigenvalue descending
idx = np.argsort(evals_E8_full)[::-1]
evals_E8_full = evals_E8_full[idx]
evecs_E8 = evecs_E8[:, idx]

print(f"\n{'Eigenvalue':>12} | Node weights (|v_i|^2)")
print("-" * 72)
for j in range(9):
    ev = evals_E8_full[j]
    vec = evecs_E8[:, j]
    weights = vec**2
    line = f"  {ev:+.6f}   | "
    for i in range(9):
        line += f"{weights[i]:.3f} "
    print(line)

# Which eigenvector is the "Perron" (all positive)?
print(f"\nPerron vector (eigenvalue 2):")
pv = evecs_E8[:, 0]
if pv[0] < 0:
    pv = -pv
for i in range(9):
    print(f"  rho_{i+1} (dim {dims_2I[i]}): v = {pv[i]:.6f}, v^2 = {pv[i]**2:.6f}")
print(f"  Note: Perron vector components proportional to dims? "
      f"Ratio: {[f'{pv[i]/pv[0]:.3f}' for i in range(9)]}")
print(f"  Dims:  {dims_2I}")
print(f"  Ratio matches dims: {all(abs(pv[i]/pv[0] - dims_2I[i]) < 0.01 for i in range(9))}")

# ================================================================
# PART 4: Graph Laplacian and Green's Function
# ================================================================
print(f"\n{'='*72}")
print("PART 4: Graph Laplacian and Green's Function")
print("="*72)

# Degree matrix
D_E8 = np.diag([sum(A_E8[i,:]) for i in range(9)])
L_E8 = D_E8 - A_E8  # Graph Laplacian

print(f"\nDegrees: {[int(D_E8[i,i]) for i in range(9)]}")
print(f"Laplacian eigenvalues: {sorted(eigvalsh(L_E8))}")

# Green's function: G = (L + epsilon*I)^{-1} (regularized)
# Or pseudoinverse: L^+ = sum_{lambda>0} (1/lambda) |v><v|
L_evals, L_evecs = eigh(L_E8)
print(f"\nLaplacian eigenvalues: {[f'{v:.6f}' for v in sorted(L_evals)]}")

# Pseudoinverse
L_pseudo = np.zeros((9,9))
for k in range(9):
    if abs(L_evals[k]) > 1e-10:
        L_pseudo += (1/L_evals[k]) * np.outer(L_evecs[:,k], L_evecs[:,k])

print(f"\nGreen's function G(0, j) from node 0 (electron):")
for j in range(9):
    gval = L_pseudo[0,j]
    # Effective distance = G(0,0) - 2*G(0,j) + G(j,j) = resistance distance
    R_eff = L_pseudo[0,0] - 2*L_pseudo[0,j] + L_pseudo[j,j]
    print(f"  G(e, rho_{j+1}): G={gval:+.4f}, R_eff={R_eff:.4f}, sqrt(R)={np.sqrt(max(R_eff,0)):.4f}")

# ================================================================
# PART 5: Propagator with mass parameter
# ================================================================
print(f"\n{'='*72}")
print("PART 5: Propagator K(e,f) = <e|(L + m^2)^{-1}|f>")
print("="*72)

# Try various m^2 values
for m2 in [0.01, 0.1, 1.0]:
    K = inv(L_E8 + m2*np.eye(9))
    print(f"\nm^2 = {m2}:")
    print(f"  {'Node':>6} {'K(e,f)':>10} {'-ln_phi(K)':>12} {'n_bare':>8} {'n_eff':>8}")
    K_e = K[0, :]  # from electron (node 0)
    for j in range(9):
        kval = K_e[j]
        if kval > 0:
            log_val = -np.log(kval) / np.log(PHI)
        else:
            log_val = float('nan')
        print(f"  rho_{j+1:1d}  {kval:>10.6f} {log_val:>12.4f}     ---      ---")

# ================================================================
# PART 6: Heat Kernel
# ================================================================
print(f"\n{'='*72}")
print("PART 6: Heat Kernel K(e,f;t) = <e|exp(-t*L)|f>")
print("="*72)

# Heat kernel: K(i,j;t) = sum_k exp(-t*lambda_k) * v_k(i) * v_k(j)
for t in [0.1, 0.5, 1.0, 2.0, PHI]:
    expL = np.zeros((9,9))
    for k in range(9):
        expL += np.exp(-t*L_evals[k]) * np.outer(L_evecs[:,k], L_evecs[:,k])

    print(f"\nt = {t:.4f}:")
    print(f"  {'Node':>6} {'K(e,f;t)':>12} {'-ln_phi(K)':>12}")
    for j in range(9):
        kval = expL[0,j]
        if kval > 0:
            log_val = -np.log(kval) / np.log(PHI)
        else:
            log_val = float('nan')
        print(f"  rho_{j+1:1d}  {kval:>12.8f} {log_val:>12.4f}")

# ================================================================
# PART 7: Weighted Adjacency with phi^5, phi^6
# ================================================================
print(f"\n{'='*72}")
print("PART 7: Weighted Adjacency Matrix")
print("="*72)

# Idea: edges on the "main arm" carry weight phi^5 (a-type),
# edges on the "branch" carry weight phi^6 (b-type)
# But which edges are which?

# In the affine E8: main arm = 0-1-2-3-4-5-6-7, branch = 5-8
# Try: main arm edges weight phi^5, branch edge weight phi^6
A_weighted_1 = np.zeros((9,9))
main_edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7)]
branch_edges = [(5,8)]
for i,j in main_edges:
    A_weighted_1[i,j] = PHI**5
    A_weighted_1[j,i] = PHI**5
for i,j in branch_edges:
    A_weighted_1[i,j] = PHI**6
    A_weighted_1[j,i] = PHI**6

print(f"\nWeighting 1: main arm = phi^5 = {PHI**5:.4f}, branch = phi^6 = {PHI**6:.4f}")
evals_w1 = sorted(eigvalsh(A_weighted_1), reverse=True)
print(f"  Eigenvalues: {[f'{v:.4f}' for v in evals_w1]}")

# Try: all edges weight 1, but transform to phi^n basis
# More physical: edge weight = phi^(|n_i - n_j|) where n is mass exponent
# But we don't know n for the intermediate nodes yet

# Try: edge weight = exp(-|n_i - n_j| * ln(phi)) = phi^(-|n_i-n_j|)
# For the mass formula on a weighted graph: m_f/m_e = product of edge weights along path

# ================================================================
# PART 8: Fermion-Node Mapping Search
# ================================================================
print(f"\n{'='*72}")
print("PART 8: Fermion-Node Mapping - Can Any Assignment Work?")
print("="*72)

fermions = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
n_targets = [n_bare[f] for f in fermions]

print(f"\nTarget exponents: {dict(zip(fermions, n_targets))}")
print(f"  e=0, mu=11, tau=17, u=3, c=16, t=26, d=5, s=11, b=19")
print(f"  Note: mu and s are DEGENERATE (both n=11)")

# Graph distances from node 0
from collections import deque
def bfs_distances(adj, start):
    n = adj.shape[0]
    dist = [-1]*n
    dist[start] = 0
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in range(n):
            if adj[u,v] > 0 and dist[v] == -1:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist

gdist = bfs_distances(A_E8, 0)
print(f"\nGraph distances from rho_1 (node 0): {gdist}")
print(f"  rho_1:0, rho_2:1, rho_3:2, rho_4:3, rho_5:4, rho_6:5, rho_7:6, rho_8:7, rho_9:6")

# The graph distances are 0,1,2,3,4,5,6,7,6
# The mass exponents are 0,3,5,11,11,16,17,19,26
# These are NOT proportional. So simple graph distance doesn't work.
# But maybe weighted distance with phi^5, phi^6 on edges?

# On affine E8, the shortest path from node 0 to node j:
# rho1->rho2->...->rhoj (along main arm) or through the branch
# Path lengths in unweighted graph are just gdist

# With weights: if main edge contributes w_m and branch contributes w_b,
# then the weighted distance to each node from 0 is:
# rho1: 0
# rho2: w_m
# rho3: 2*w_m
# rho4: 3*w_m
# rho5: 4*w_m
# rho6: 5*w_m
# rho7: 6*w_m
# rho8: 7*w_m
# rho9: 5*w_m + w_b

# We need these to match mass exponents for SOME assignment of fermions to nodes.
# Since rho1 = electron (n=0), the ordering along the main arm gives:
# 0, w_m, 2w_m, 3w_m, 4w_m, 5w_m, 6w_m, 7w_m
# and the branch gives 5w_m + w_b

# For w_m = A1 = 5: 0, 5, 10, 15, 20, 25, 30, 35 and 25+w_b
# For w_m = B1 = 6: 0, 6, 12, 18, 24, 30, 36, 42 and 30+w_b
# Neither matches {0,3,5,11,16,17,19,26} since the weighted distances are evenly spaced

print(f"\n  Simple weighted distance with uniform w_m CANNOT match:")
print(f"  The mass exponents are NOT evenly spaced.")
print(f"  We need NON-UNIFORM edge weights or a different mechanism (propagator, not distance).")

# ================================================================
# PART 9: The KEY insight - Spectral Propagator
# ================================================================
print(f"\n{'='*72}")
print("PART 9: Spectral Propagator on Affine E8")
print("="*72)

# The propagator on a graph is NOT just the shortest path distance.
# It's the FULL sum over all paths, weighted by the adjacency spectrum.
#
# For adjacency matrix A with eigenvalues {lambda_k} and eigenvectors {v_k}:
# (A^n)_{ij} = sum_k lambda_k^n * v_k(i) * v_k(j)
#
# The resolvent/propagator:
# G(z)_{ij} = <i| (z - A)^{-1} |j> = sum_k v_k(i)*v_k(j) / (z - lambda_k)
#
# For mass generation, consider:
# K_{ij} = <i| f(A) |j> for some function f
#
# If f(A) = A^n (powers of adjacency), then K_{ij} counts walks of length n.
# If f(A) = exp(t*A) (heat kernel on adjacency), diffusion interpretation.
# If f(A) = (z - A)^{-1} (resolvent), propagator interpretation.

# The mass exponent n_f should be related to:
# m_f/m_e = phi^{n_f} = (something involving the propagator from e to f)

# Key observation: if spectrum = {0, +-1/phi, +-1, +-phi, +-2},
# then phi APPEARS IN THE SPECTRUM. So functions of A naturally produce
# powers of phi!

# Test: (A^n)_{0j} for various n
print(f"\nPowers of adjacency matrix A^n, element (0,j):")
An = np.eye(9)
print(f"\n{'n':>3} | {'rho'+str(i+1):>8s}" for i in range(9))
header = f"{'n':>3} |"
for i in range(9):
    header += f" rho{i+1:1d}  "
print(header)
print("-" * 72)

for n in range(16):
    row = f"{n:>3} |"
    for j in range(9):
        row += f" {An[0,j]:>5.1f} "
    print(row)
    An = An @ A_E8

# ================================================================
# PART 10: Resolvent and phi-structure
# ================================================================
print(f"\n{'='*72}")
print("PART 10: Resolvent G(z) = (z*I - A)^{-1} at Special z Values")
print("="*72)

# Try z = phi^k for various k
for z_val, z_name in [(PHI**2, 'phi^2'), (PHI**3, 'phi^3'), (3.0, '3'),
                       (PHI**2 + 1, 'phi^2+1'), (2+PHI, '2+phi'),
                       (PHI + 1/PHI, 'phi+1/phi=sqrt5')]:
    try:
        G = inv(z_val * np.eye(9) - A_E8)
        print(f"\nz = {z_name} = {z_val:.4f}:")
        row = f"  G(0,j): "
        for j in range(9):
            row += f"{G[0,j]:>8.4f} "
        print(row)
        # Log_phi of |G(0,j)/G(0,0)|
        row = f"  -log_phi: "
        for j in range(9):
            ratio = abs(G[0,j])
            if ratio > 1e-15:
                log_val = -np.log(ratio / abs(G[0,0])) / np.log(PHI)
                row += f"{log_val:>8.3f} "
            else:
                row += f"     inf "
        print(row)
    except np.linalg.LinAlgError:
        print(f"\nz = {z_name}: singular (z is an eigenvalue)")

# ================================================================
# PART 11: The Dirac Operator Approach
# ================================================================
print(f"\n{'='*72}")
print("PART 11: Finite Dirac Operator D_F")
print("="*72)

# In NCG, the finite Dirac operator D_F has eigenvalues proportional to
# the Yukawa couplings (fermion masses after SSB).
# Our D_F should be built from the McKay graph.
#
# Approach: D_F = A (adjacency) or D_F = some weighted version
# Then the "mass" of fermion at node j is related to an eigenvalue or
# matrix element of D_F.
#
# For the standard NCG, the fermion mass matrix IS D_F restricted to
# the appropriate sector.
#
# Key idea: if we USE the eigenvectors of A as the fermion basis,
# then each fermion has eigenvalue lambda_k, and the mass is phi^{f(lambda_k)}.

# The eigenvectors of A are already computed (evecs_E8)
# Each eigenvalue is associated with a specific pattern of node weights.

# IMPORTANT: the spectrum has 9 eigenvalues and we have 9 fermions.
# Can we MAP eigenvalues to fermions such that lambda_k -> n_k?

print(f"\n9 eigenvalues of A: {[f'{v:.6f}' for v in evals_E8_full]}")
print(f"9 mass exponents:   {sorted(n_targets)}")
print(f"\nCan we find f such that f(lambda_k) = n_k for some assignment?")

# Sort both
sorted_evals = sorted(evals_E8_full, reverse=True)
sorted_n = sorted(n_targets, reverse=True)  # 26,19,17,16,11,11,5,3,0

print(f"\n  Sorted eigenvalues: {[f'{v:.4f}' for v in sorted_evals]}")
print(f"  Sorted exponents:   {sorted_n}")

# Try: f(lambda) = a * lambda + b (linear)
# Or: f(lambda) = a * lambda^2 + b (quadratic)
# Or: f(lambda) = a / (z - lambda) (resolvent-type)

# Linear: if lambda ranges from -2 to +2 and n from 0 to 26,
# then f(lambda) = 13*(lambda/2 + 1) = 6.5*(lambda + 2)
# This gives: f(2)=26, f(phi)=23.6, f(1)=19.5, f(1/phi)=17.0, f(0)=13,
#             f(-1/phi)=9.0, f(-1)=6.5, f(-phi)=2.4, f(-2)=0
lin_map = [(6.5*(v+2), v) for v in sorted_evals]
print(f"\n  Linear map f(lambda) = 6.5*(lambda+2):")
for n_lin, lam in lin_map:
    print(f"    lambda={lam:+.4f} -> n={n_lin:.1f}")

# Compare with targets
print(f"\n  Linear map values: {sorted([6.5*(v+2) for v in sorted_evals], reverse=True)}")
print(f"  Target values:     {sorted_n}")
# Check: f(2)=26 OK, f(phi)=23.56 vs 19 NO, f(0)=13 vs 11 NO
# Linear doesn't match

# Try: match largest to largest, find optimal linear
from itertools import permutations

# Brute force: try all 9!/2 permutations (since mu=s degenerate)
# Actually 9! = 362880, too many. Let's use sorted order and check
# if any monotonic mapping works.

# Since eigenvalues are symmetric about 0, and exponents are all >= 0,
# a monotonic map won't give both 0 and 26 from -2 and +2.
# Unless we use |lambda| or lambda^2.

# Try f(lambda) = a * lambda^2 + b
# lambda^2 values: 4, phi^2=2.618, 1, 1/phi^2=0.382, 0
# These are 5 distinct values for 9 eigenvalues (each +-pair gives same lambda^2)
# So lambda^2 can only distinguish 5 groups, but we need 8 (after mu/s degeneracy)
# This doesn't have enough resolution.

# Try: use eigenvector components at node 0 (electron node)
print(f"\n  Eigenvector component at node 0 (electron):")
for j in range(9):
    ev = evals_E8_full[j]
    v0 = evecs_E8[0, j]
    print(f"    lambda={ev:+.6f}, v_0 = {v0:+.8f}, v_0^2 = {v0**2:.8f}")

# ================================================================
# PART 12: The Cayley Graph on A5 Irreps
# ================================================================
print(f"\n{'='*72}")
print("PART 12: A5 Irrep Tensor Product Graph (McKay with 3-dim)")
print("="*72)

# A5 has 5 irreps: 1, 3, 3', 4, 5
# Tensor products with rho_3 (3-dim):
# 1 x 3 = 3
# 3 x 3 = 1 + 3 + 5
# 3' x 3 = 4 + 5
# 4 x 3 = 3' + 4 + 5
# 5 x 3 = 3 + 3' + 4 + 5

# Adjacency matrix (5x5), with self-loops
# Nodes: 0=1, 1=3, 2=3', 3=4, 4=5
A_A5 = np.array([
    [0, 1, 0, 0, 0],  # 1 x 3 = 3
    [1, 1, 0, 0, 1],  # 3 x 3 = 1 + 3 + 5
    [0, 0, 0, 1, 1],  # 3' x 3 = 4 + 5
    [0, 0, 1, 1, 1],  # 4 x 3 = 3' + 4 + 5
    [0, 1, 1, 1, 1],  # 5 x 3 = 3 + 3' + 4 + 5
], dtype=float)

labels_A5 = ['1', '3', "3'", '4', '5']
dims_A5 = [1, 3, 3, 4, 5]

print(f"\nA5 McKay graph (using 3-dim irrep):")
print(f"  Nodes: {labels_A5}, dims: {dims_A5}")
print(f"  Adjacency matrix:")
for i in range(5):
    print(f"    {labels_A5[i]:>3}: {A_A5[i,:]}")

evals_A5 = sorted(eigvalsh(A_A5), reverse=True)
print(f"\n  Eigenvalues: {[f'{v:.6f}' for v in evals_A5]}")

# The Cayley Laplacian eigenvalues (from character theory)
# lambda_i = dim(rho_3) - chi_3(class)/dim(rho_i) * dim(rho_3)
# Actually: for the tensor product graph, the eigenvalues are:
# chi_3(class_j) / dim(rho_i) for the appropriate formula
# Let me compute from character theory

# A5 character table
# Classes: 1a (1), 2a (15), 3a (20), 5a (12), 5b (12)
# Order:    1       2         3        5        5
chi = {
    '1':  [1,  1,  1,  1,  1],
    '3':  [3, -1,  0, PHI, PHIp],
    "3'": [3, -1,  0, PHIp, PHI],
    '4':  [4,  0,  1, -1, -1],
    '5':  [5,  1, -1,  0,  0],
}
class_sizes = [1, 15, 20, 12, 12]

# Cayley Laplacian eigenvalue for irrep rho_i:
# lambda_i = k - (k/d_i) * sum_g_in_S chi_i(g)
# where S is the generating set and k = |S|
# For tensor product with rho_3: the "generating set" is rho_3,
# but this is different from the Cayley graph.

# Actually, the adjacency matrix of the tensor product graph
# (A_A5)_{ij} = multiplicity of rho_j in rho_i tensor rho_3
# Its eigenvalues are chi_3(class_j) for j = 0,...,4
# by the general theory of McKay graphs.

# Wait, that's for the standard McKay graph where A_{ij} = n_{ij3} (the tensor
# product multiplicity). The eigenvalues of this matrix are:
# lambda_k = chi_3(g_k) / chi_1(g_k) = chi_3(g_k) for k-th conjugacy class
# But this isn't right either - the eigenvalues relate to the Fourier transform.

# For a commutative fusion ring (which A5 irreps form), the McKay matrix
# M_3 has eigenvalues = chi_3(class_j) evaluated at the "class" points.
# But there are 5 classes and 5 irreps, so 5 eigenvalues.
# The eigenvectors are rows of the character table (properly normalized).

print(f"\n  Character values of rho_3 on each conjugacy class:")
for k, (cls_name, chi_val) in enumerate(zip(['1a','2a','3a','5a','5b'], chi['3'])):
    print(f"    chi_3({cls_name}) = {chi_val:.6f}")

print(f"\n  If McKay graph eigenvalues = chi_3(class_j):")
chi3_vals = sorted(chi['3'], reverse=True)
print(f"    {{chi_3}} = {[f'{v:.6f}' for v in chi3_vals]}")
print(f"    = {{3, phi, 0, -1, phi'=-1/phi}}")
print(f"    A_A5 eigenvalues: {[f'{v:.6f}' for v in evals_A5]}")

# Check match
match_chi = sorted(chi['3'], reverse=True)
match_eig = sorted(evals_A5, reverse=True)
print(f"\n  Match: {all(abs(a-b) < 1e-6 for a,b in zip(match_chi, match_eig))}")

# ================================================================
# PART 13: Connecting the Two Graphs
# ================================================================
print(f"\n{'='*72}")
print("PART 13: Connection Between Affine E8 and A5 Spectra")
print("="*72)

print(f"\n  Affine E8 spectrum: {{0, +-1/phi, +-1, +-phi, +-2}}")
print(f"  = 9 values, symmetric about 0")
print(f"\n  A5 tensor graph spectrum: {{phi', -1, 0, phi, 3}}")
print(f"  = 5 values, NOT symmetric")
print(f"\n  RELATIONSHIP:")
print(f"  The A5 spectrum {{phi', -1, 0, phi, 3}} is related to")
print(f"  the E8 spectrum by: the 5 non-negative E8 eigenvalues")
print(f"  {{0, 1/phi, 1, phi, 2}} are HALF of the symmetric pairs,")
print(f"  and the A5 eigenvalues involve phi and phi' = -1/phi.")
print(f"\n  Both graphs 'know' about phi through their spectra!")

# ================================================================
# PART 14: Green's Function on A5 Graph
# ================================================================
print(f"\n{'='*72}")
print("PART 14: Propagator on A5 Tensor Graph")
print("="*72)

# The A5 graph is not a tree (has self-loops and cycles).
# Let's compute the Green's function.

L_A5 = np.diag(np.sum(A_A5, axis=1)) - A_A5
print(f"\nA5 Laplacian eigenvalues: {sorted(eigvalsh(L_A5))}")

# Resolvent at various z
for z_val, z_name in [(PHI, 'phi'), (PHI**2, 'phi^2'), (3+0.01, '3.01'),
                       (4.0, '4'), (A1+0.0, '5=a1'), (B1+0.0, '6=b1')]:
    try:
        G_A5 = inv(z_val * np.eye(5) - A_A5)
        print(f"\n  z={z_name}: G(1, j) = ", end="")
        for j in range(5):
            print(f"{G_A5[0,j]:+.4f} ", end="")
        print()
    except:
        print(f"\n  z={z_name}: SINGULAR")

# ================================================================
# PART 15: The Mass Formula from Spectral Data
# ================================================================
print(f"\n{'='*72}")
print("PART 15: Can the Mass Formula Emerge from Spectral Data?")
print("="*72)

# The mass formula is n_f = a1*a + b1*b where (a,b) are Z[phi] coordinates
# On the McKay graph, we need to connect this to spectral properties.
#
# Key observation: in Z[phi], z = a + b*phi has:
#   Trace(z) = z + z' = 2a + b (since phi + phi' = 1)
#   Norm(z) = z*z' = a^2 + ab - b^2
#
# The mass exponent n = a1*a + b1*b = 5a + 6b
# But also: n = a1*(a+b) + b = a1*Tr(z) + b*(1-a1+a1) ... hmm
# n = 5a + 6b = 5a + 5b + b = 5(a+b) + b = 5*Tr(z)/... no, Tr(z) = 2a+b
#
# n = 5a + 6b
# Tr(z) = 2a + b => a = (Tr - b*1) ... this isn't clean
#
# Alternative: n = a1*a + b1*b = a1*a + (a1+1)*b = a1*(a+b) + b
# Since a+b = Tr(z) - a:  n = a1*Tr(z) - a1*a + a1*b + b ... circular
#
# Direct: n = a1*a + b1*b is simply a linear functional on Z^2.
# In matrix form: n = [a1, b1] . [a, b]^T

# The spectral interpretation:
# If we think of the McKay graph as encoding the Z[phi] structure,
# then each node rho_k has some Z[phi] coordinate (a_k, b_k)
# and the mass exponent is the evaluation of the linear functional (a1, b1).

# For the affine E8 graph (9 nodes), we need to assign (a,b) to each node.
# Can the graph structure determine these assignments?

# The Perron vector (eigenvector for lambda_max = 2) has components
# proportional to the dimensions: [1,2,3,4,5,6,4,2,3]
# These are the Dynkin marks of affine E8.

# The DIMENSION is the norm of the representation in Z[phi]?
# For A5: dims = {1, 3, 3, 4, 5}
# In our framework: 3 = 3 (3-dim), 3' = 3 (3-dim), 4 = 4, 5 = 5 = a1
# And 1+3+3+4+5 = 16 = (a1-1)^2

print(f"\n  A5 irrep dimensions: {dims_A5}")
print(f"  Sum = {sum(dims_A5)} = (a1-1)^2 = {(A1-1)**2}")
print(f"  2I irrep dimensions: {dims_2I}")
print(f"  Sum = {sum(dims_2I)} = h(E8) = {A1*(A1+1)}")

# For 2I: sum of dims = 30 = h(E8) = Coxeter number
# For 2I: sum of dims^2 = 1+4+9+16+25+36+16+4+9 = 120 = |2I| = N
print(f"  Sum of dims^2 = {sum(d**2 for d in dims_2I)} = N = {A1*24}")

# ================================================================
# PART 16: Spectral Action Approach
# ================================================================
print(f"\n{'='*72}")
print("PART 16: Spectral Action Tr(f(D_F/Lambda))")
print("="*72)

# The spectral action is S = Tr(f(D/Lambda))
# For the finite Dirac operator D_F (our adjacency or Laplacian),
# the moments are:
# a_0 = Tr(1) = 9 (or dim H_F)
# a_2 = Tr(D_F^2)
# a_4 = Tr(D_F^4)
# etc.

# For D_F = A (adjacency):
print(f"\nD_F = A (adjacency of affine E8):")
for k in range(1, 7):
    Ak = np.linalg.matrix_power(A_E8, k)
    tr = np.trace(Ak)
    print(f"  Tr(A^{k}) = {tr:.0f}", end="")
    if k == 2:
        print(f"  = 2 * (number of edges) = 2 * {len(edges_E8)} = {2*len(edges_E8)}", end="")
    print()

# Using the eigenvalues:
print(f"\n  Cross-check with eigenvalues:")
for k in [2, 4, 6]:
    tr_from_evals = sum(ev**k for ev in evals_E8_full)
    print(f"  sum(lambda^{k}) = {tr_from_evals:.6f}")

# a_2 = Tr(D_F^2) in our framework should be = rank(E8) = 8?
print(f"\n  Tr(A^2) = {np.trace(A_E8 @ A_E8):.0f}")
print(f"  Compare: rank(E8) = 8, dim(G) = 12, h(E8) = 30")
print(f"  Tr(A^2) = 16 = 2 * 8 = 2 * rank(E8)!")
print(f"  Also: 16 = (a1-1)^2 = N(z_Planck) = fermions/generation")

# ================================================================
# PART 17: Dimension-Weighted Propagator
# ================================================================
print(f"\n{'='*72}")
print("PART 17: Dimension-Weighted Propagator")
print("="*72)

# Idea: weight each node by its dimension d_i.
# The propagator on a weighted graph:
# D_ij = A_ij * sqrt(d_i * d_j) (or some other weighting)
# This makes the operator "know" about the irrep structure.

# Try: D_F = diag(d) * A * diag(d) / 30 (normalized by |2I|/4 or similar)
d = np.array(dims_2I, dtype=float)
D_dim = np.diag(d) @ A_E8 @ np.diag(d)
print(f"\nD_dim = diag(d) * A * diag(d):")
evals_dim = sorted(eigvalsh(D_dim), reverse=True)
print(f"  Eigenvalues: {[f'{v:.4f}' for v in evals_dim]}")
print(f"  Tr(D_dim^2) = {np.trace(D_dim @ D_dim):.1f}")

# Try: D_F = A with each edge weight = sqrt(d_i * d_j)
D_sqrt = np.zeros((9,9))
for i,j in edges_E8:
    w = np.sqrt(d[i] * d[j])
    D_sqrt[i,j] = w
    D_sqrt[j,i] = w
print(f"\nD_sqrt (edge weight = sqrt(d_i * d_j)):")
evals_sqrt = sorted(eigvalsh(D_sqrt), reverse=True)
print(f"  Eigenvalues: {[f'{v:.4f}' for v in evals_sqrt]}")
print(f"  Tr(D_sqrt^2) = {np.trace(D_sqrt @ D_sqrt):.1f}")
print(f"  Compare: sum(d_i*d_j) over edges = {sum(d[i]*d[j] for i,j in edges_E8):.1f}")

# ================================================================
# PART 18: Transfer Matrix / Path Amplitudes
# ================================================================
print(f"\n{'='*72}")
print("PART 18: Transfer Matrix - Path Amplitudes on E8")
print("="*72)

# The key insight from aaaa.txt: "each step contributes phi^5 or phi^6"
# This means the TRANSFER MATRIX has entries phi^5 or phi^6.
# Define T_ij = phi^{w_ij} where w_ij depends on the edge type.
#
# For the mass formula n = 5a + 6b, we need the path amplitude
# <e|T^k|f> to give phi^{n_f} for the appropriate k (or sum over all k).
#
# If we define T = phi^5 * P_main + phi^6 * P_branch
# where P_main and P_branch are the adjacency matrices of the two edge types,
# this is just the weighted adjacency from Part 7.
#
# But we need a DIFFERENT weighting: not phi^5 on ALL main edges,
# but phi^5 on "a-type steps" and phi^6 on "b-type steps".
#
# The question is: which edges are a-type and which are b-type?
#
# In Z[phi], multiplying by phi adds (0,1) to (a,b):
#   z * phi = (a + b*phi) * phi = b + (a+b)*phi = (b, a+b)
# So (a,b) -> (b, a+b), i.e., da=b-a, db=a.
#
# The "a-direction" in Z[phi] is adding 1: (a,b) -> (a+1, b), dz = 1
# The "b-direction" is adding phi: (a,b) -> (a, b+1), dz = phi
#
# So the mass exponent change:
# a-step: dn = 5*1 + 6*0 = 5 (a1)
# b-step: dn = 5*0 + 6*1 = 6 (b1)
# phi-step: dn = 5*0 + 6*1 = 6 (same as b-step!)
# 1-step: dn = 5*1 + 6*0 = 5 (same as a-step!)

# So the McKay graph needs edges that correspond to "add 1" and "add phi"
# in the Z[phi] structure.

# In the TENSOR PRODUCT interpretation:
# Tensoring with rho_2 (2-dim, the natural rep of SU(2)) is one step.
# The 2-dim rep of 2I gives phi (its character on 5-cycles).
# So each tensor step "adds" phi to the Z[phi] coordinate? Not exactly...

# Let's think about it from a different angle.
# The CHARACTER of rho_2 on the 5-cycle class is phi.
# So rho_2 "carries" a factor of phi.
# If the mass exponent is determined by the character,
# then n_f = (some function of chi_f evaluated at the 5-cycle class)

# For A5 irreps:
# chi_1(5a) = 1,  chi_3(5a) = phi,  chi_3'(5a) = phi', chi_4(5a) = -1, chi_5(5a) = 0
# And the fermion exponents: need to match...

print(f"\n  A5 characters on 5-cycles and mass exponents:")
print(f"  {'Irrep':>5} {'dim':>4} {'chi(5a)':>8} {'chi(5b)':>8}")
for name in ['1', '3', "3'", '4', '5']:
    d = chi[name][0]
    c5a = chi[name][3]
    c5b = chi[name][4]
    print(f"  {name:>5} {d:>4} {c5a:>+8.4f} {c5b:>+8.4f}")

# ================================================================
# PART 19: The Fundamental Discovery - exp(A*t) Structure
# ================================================================
print(f"\n{'='*72}")
print("PART 19: exp(A*t) at Golden Ratio Times")
print("="*72)

# If the spectrum is {0, +-1/phi, +-1, +-phi, +-2}, then
# exp(A*t) has entries that are sums of exp(+-2t), exp(+-phi*t),
# exp(+-t), exp(+-t/phi), and exp(0) = 1.
#
# At t = ln(phi), we get factors phi^{+-2}, phi^{+-phi}, phi^{+-1},
# phi^{+-1/phi}, phi^0 = 1.
#
# The mass formula needs phi^n for integer n, but phi^{phi} is NOT
# phi^integer. So t = ln(phi) is not the right scale.
#
# At t = k*ln(phi) for integer k, we get phi^{+-2k}, etc.
# Still not matching 5a+6b pattern.
#
# BUT: if D_F is NOT just A but a FUNCTION of A, say
# D_F = a1*A + b1*A_branch or similar, then the eigenvalues change.

# Try: what function g(x) maps {0, 1/phi, 1, phi, 2} -> {0, 5, ?, ?, 26}?
# Or to any subset of {0,3,5,11,16,17,19,26}?

# Simplest: g(2) = 26, g(0) = 0 => linear g(x) = 13x
# g(phi) = 13*phi = 21.03... not integer
# g(1) = 13, g(1/phi) = 8.03... These are FIBONACCI NUMBERS!
# 13 = F_7, 8 = F_6, 21 = F_8, 0 = F_0

# WAIT! 13*phi = 13*(1+sqrt(5))/2 = 21.03... and F_8 = 21!
# And 13*1/phi = 13*(sqrt(5)-1)/2 = 8.03... and F_6 = 8!
# And 13*1 = 13 = F_7!
# And 13*2 = 26 and F? ... not Fibonacci

# Actually: 13 * {0, 1/phi, 1, phi, 2} ~ {0, 8.03, 13, 21.03, 26}
# The integer parts are {0, 8, 13, 21, 26}
# Hmm, 8 and 21 are Fibonacci, 13 is Fibonacci, but 26 = 2*13

# More precisely: C/(a1^2+1) * eigenvalue where C = a1^2+1 ...
# Or: (a1^2+1)/2 * eigenvalue = 13 * eigenvalue
# 13 * {0, 1/phi, 1, phi, 2} = {0, 13/phi, 13, 13*phi, 26}
# = {0, 8.034, 13, 21.034, 26}

print(f"\n  13 * eigenvalues (positive): ")
for v in [0, 1/PHI, 1, PHI, 2]:
    print(f"    13 * {v:.4f} = {13*v:.4f}")
print(f"  These approximate {{0, 8, 13, 21, 26}} = Fibonacci-like!")
print(f"  F_6=8, F_7=13, F_8=21")

# What if we scale by (a1^2+1)/2 = 13?
# The mass exponents are NOT {0, 8, 13, 21, 26} but {0, 3, 5, 11, 16, 17, 19, 26}
# The 9 mass exponents include duplicates (mu=s=11)

# The 5 UNIQUE A5 eigenvalues * 13 give {0, 8, 13, 21, 26}
# But our mass exponents on the main chain would be different

# For 2I (9 eigenvalues): 13 * {-2,-phi,-1,-1/phi,0,1/phi,1,phi,2}
# = {-26, -21.03, -13, -8.03, 0, 8.03, 13, 21.03, 26}
# Using floor/round: {-26, -21, -13, -8, 0, 8, 13, 21, 26}
# These are symmetric! Mass exponents are NOT symmetric.

# ================================================================
# PART 20: Key Test - Dimension-Weighted Spectral Function
# ================================================================
print(f"\n{'='*72}")
print("PART 20: Mass from Dimension-Weighted Characters")
print("="*72)

# For A5 irreps, the Cayley Laplacian eigenvalues are:
# lambda_rho = 12(1 - chi_rho(5a)/dim(rho))
# where 12 = vertex degree = 2(a1+1) = number of 5-cycles = |class 5a| + |class 5b|

# Wait, |5a| = 12 and |5b| = 12, total 24 five-cycles.
# But the Cayley graph on A5 with generating set = class 5a (12 elements):

print(f"\nCayley Laplacian eigenvalues (generating set = class 5a, |S|=12):")
print(f"  lambda_rho = 12 - 12*chi_rho(5a)/dim(rho)")
for name in ['1', '3', "3'", '4', '5']:
    d_i = chi[name][0]
    c5a = chi[name][3]
    lam = 12 - 12*c5a/d_i
    print(f"  rho={name:>3}: chi(5a)/d = {c5a:.4f}/{d_i} = {c5a/d_i:.6f}, "
          f"lambda = {lam:.6f}")
    # For 3: 12 - 12*phi/3 = 12 - 4*phi = 5.528
    # For 3': 12 - 12*phi'/3 = 12 + 4/phi = 12 + 4*(phi-1) = 8+4*phi = 14.472
    # For 4: 12 - 12*(-1)/4 = 12 + 3 = 15
    # For 5: 12 - 12*0/5 = 12
    # For 1: 12 - 12*1/1 = 0

# Now the question: can we get mass exponents from these Cayley eigenvalues?
# lambda_3 = 12 - 4*phi = 5.528, lambda_3' = 8 + 4*phi = 14.472
# Spectral gap = 4*sqrt(5) = 8.944 (confirmed in exp312)

# What if mass exponent n is related to Cayley eigenvalue lambda?
# For A5 irreps, we have 5 eigenvalues for 5 "types" of fermions.
# But we have 9 fermions assigned to these 5 types via the (a,b) quantum numbers.

# The A5 irreps decompose the 12-dim permutation rep as 12 = 1 + 3 + 3' + 5
# And the fermion assignments come from the McKay mass correspondence (exp323).

# Let me try: n_f = alpha_s_inv * lambda_f / 12 * 26 or similar scaling
alpha_s_inv = 2*PHI**3
for name in ['1', '3', "3'", '4', '5']:
    d_i = chi[name][0]
    c5a = chi[name][3]
    lam = 12 - 12*c5a/d_i
    # Scale by a1 or some framework constant
    scaled = lam * (A1**2 + 1) / 12  # = lambda * 26/12 = lambda * 13/6
    print(f"  rho={name:>3}: lambda={lam:.4f}, lambda*13/6 = {scaled:.4f}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*72}")
print("SUMMARY")
print("="*72)

print(f"""
KEY FINDINGS:

1. AFFINE E8 SPECTRUM VERIFICATION:
   The adjacency matrix of the affine E8 Dynkin diagram has eigenvalues
   {{0, +-1/phi, +-1, +-phi, +-2}} = CONFIRMED (or NOT) above.
   This is the McKay graph of 2I (binary icosahedral group).

2. A5 TENSOR GRAPH:
   The McKay graph of A5 with 3-dim irrep (5x5) has eigenvalues
   equal to chi_3(class_j) = {{phi', -1, 0, phi, 3}}.
   These are the CHARACTER VALUES of rho_3 on the 5 conjugacy classes.

3. PROPAGATOR RESULTS:
   [To be analyzed from the numerical results above]

4. GALOIS STRUCTURE:
   The Galois automorphism phi -> -1/phi acts on the E8 spectrum as:
   - Swaps phi <-> -1/phi (the two irrational eigenvalues)
   - Fixes the rational eigenvalues {{0, +-1, +-2}}
   This mirrors exactly how Galois acts on the mass formula!

5. FIBONACCI CONNECTION:
   13 * {{0, 1/phi, 1, phi, 2}} ~ {{0, 8, 13, 21, 26}}
   The first three are FIBONACCI NUMBERS (F_6, F_7, F_8)!
   26 = 2*13 = top mass exponent.

OPEN QUESTIONS:
- Which function f maps eigenvalues to mass exponents?
- How do the 9 fermions map to the 9 E8 nodes?
- Can the correction formula emerge from the spectral structure?
""")

print("Done.")
