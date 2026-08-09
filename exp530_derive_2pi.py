"""
EXP-530: Deriving 2*pi in the Alpha Equation
==============================================

STATUS: The alpha equation is  2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0.
  - 4*a1*phi^4: STRUCTURAL (spectral, KK identification)
  - 1: normalization (Vieta product alpha*alpha' = 1/(2*pi))
  - 2*pi: POSTULATED.  Paper claims Hopf holonomy (10 edges * pi/5 = 2*pi),
          but the IDENTIFICATION  delta(1/alpha) = -2*pi*alpha  is not derived.

GOAL: Can 2*pi emerge from a LATTICE SELF-ENERGY calculation on the 600-cell?

APPROACHES:
  A. Lattice propagator G = L^{-1} (pseudoinverse of Laplacian).
     Compute vacuum-polarization-like quantities: Tr(G^2), sum over edges, etc.
  B. Spectral zeta function at special values.
  C. Eta invariant of S^3/2I (spectral asymmetry).
  D. Hopf-fiber-restricted self-energy.
  E. Direct: does any graph-spectral object equal exactly 2*pi?
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, SQRT5, a1, b1, N, h, degree, d_ST, rank_E8,
                      alpha_em, inv_alpha, alpha_s, sin2_tW, N_gen,
                      SPECTRUM_600CELL, L0, L1, L2, L3, L4, L5, L6, L7, L8,
                      ALPHA_COEFF)
from commons.cell600 import build_600cell

print("=" * 70)
print("EXP-530: DERIVING 2*pi IN THE ALPHA EQUATION")
print("=" * 70)

# Framework:
# 1/alpha = C - 2*pi*alpha   where  C = 4*a1*phi^4
# => 2*pi*alpha^2 - C*alpha + 1 = 0
C = ALPHA_COEFF   # = 4*5*phi^4 = 137.082...
print(f"\n  Alpha equation: 2*pi*x^2 - {C:.4f}*x + 1 = 0")
print(f"  alpha = {alpha_em:.8f},  1/alpha = {inv_alpha:.6f}")
print(f"  Coefficient to derive: 2*pi = {2*np.pi:.6f}")

# ============================================================
# SECTION 1: BUILD 600-CELL AND LAPLACIAN
# ============================================================
print(f"\n{'='*70}")
print("SECTION 1: 600-CELL GRAPH AND LAPLACIAN")
print(f"{'='*70}")

vertices, A, L_mat = build_600cell()
n_vert = len(vertices)

# Extract edge list from adjacency matrix
edges = []
for i in range(n_vert):
    for j in range(i+1, n_vert):
        if A[i, j] > 0.5:
            edges.append((i, j))
n_edge = len(edges)
print(f"  Vertices: {n_vert},  Edges: {n_edge}")
print(f"  Vertex degree: {degree}")

# Verify spectrum
eigenvalues = np.linalg.eigvalsh(L_mat)
eigenvalues.sort()
print(f"  Laplacian eigenvalues: min={eigenvalues[0]:.4f}, max={eigenvalues[-1]:.4f}")
print(f"  Tr(L) = {np.trace(L_mat):.0f} = N*degree = {n_vert*degree}")

# Eigenvalue data from spectrum
eig_data = [(lam, mult, name) for lam, mult, name, _ in SPECTRUM_600CELL]
print(f"\n  Eigenvalue | Mult | Name")
for lam, mult, name in eig_data:
    print(f"  {lam:10.4f} | {mult:4d} | {name}")

# ============================================================
# SECTION 2: PROPAGATOR (PSEUDOINVERSE OF LAPLACIAN)
# ============================================================
print(f"\n{'='*70}")
print("SECTION 2: GRAPH PROPAGATOR G = L^+")
print(f"{'='*70}")

# Pseudoinverse: project out zero mode, then invert
# L = U @ diag(eigenvalues) @ U^T
evals, evecs = np.linalg.eigh(L_mat)

# Sort
idx_sort = np.argsort(evals)
evals = evals[idx_sort]
evecs = evecs[:, idx_sort]

# Build pseudoinverse
G_evals = np.zeros_like(evals)
for i in range(len(evals)):
    if evals[i] > 1e-10:
        G_evals[i] = 1.0 / evals[i]

G_mat = evecs @ np.diag(G_evals) @ evecs.T

# Check: L @ G should be projector onto non-zero eigenspace
LG = L_mat @ G_mat
print(f"  ||L @ G - (I - J/N)||_F = {np.linalg.norm(LG - (np.eye(n_vert) - np.ones((n_vert,n_vert))/n_vert)):.2e}")

# Propagator properties
G_diag = np.diag(G_mat)
print(f"\n  G(x,x) = {G_diag[0]:.8f}  (constant on vertex-transitive graph)")
print(f"  Tr(G) = {np.trace(G_mat):.8f}")

# Tr(G) = sum_{lambda!=0} mult/lambda
Tr_G_exact = sum(mult/lam for lam, mult, _ in eig_data if lam > 0)
print(f"  Tr(G) exact = {Tr_G_exact:.8f}")
print(f"  Tr(G)/N = {Tr_G_exact/N:.8f}")

# ============================================================
# SECTION 3: VACUUM POLARIZATION ON THE GRAPH
# ============================================================
print(f"\n{'='*70}")
print("SECTION 3: VACUUM POLARIZATION QUANTITIES")
print(f"{'='*70}")

# Several natural "self-energy" quantities:

# VP-1: Tr(G^2) = sum mult/lambda^2
Tr_G2 = np.trace(G_mat @ G_mat)
Tr_G2_exact = sum(mult/lam**2 for lam, mult, _ in eig_data if lam > 0)
print(f"  (a) Tr(G^2) = {Tr_G2_exact:.8f}")
print(f"      Tr(G^2) / N = {Tr_G2_exact/N:.8f}")

# VP-2: Sum over edges of G(i,j)^2
G_edge_sum = 0
for (i, j) in edges:
    G_edge_sum += G_mat[i, j]**2
print(f"\n  (b) Sum_edges G(i,j)^2 = {G_edge_sum:.8f}")
print(f"      / N_edges = {G_edge_sum/n_edge:.8f}")
print(f"      / N = {G_edge_sum/N:.8f}")

# VP-3: Sum over edges of [G(i,i)*G(j,j) - G(i,j)^2]  (bubble diagram)
G_bubble = 0
for (i, j) in edges:
    G_bubble += G_mat[i, i] * G_mat[j, j] - G_mat[i, j]**2
print(f"\n  (c) Bubble = Sum_edges [G(ii)*G(jj) - G(ij)^2] = {G_bubble:.8f}")
print(f"      / N_edges = {G_bubble/n_edge:.8f}")
print(f"      / N = {G_bubble/N:.8f}")

# VP-4: Average of G(i,j) over edges (mean-field propagator)
G_edge_avg = sum(G_mat[i, j] for (i, j) in edges) / n_edge
print(f"\n  (d) <G(i,j)>_edges = {G_edge_avg:.8f}")

# VP-5: G(i,j) for adjacent vertices (should be constant by vertex-transitivity)
G_adj_values = [G_mat[i, j] for (i, j) in edges]
G_adj_mean = np.mean(G_adj_values)
G_adj_std = np.std(G_adj_values)
print(f"\n  (e) G(i,j) for adjacent i~j: mean = {G_adj_mean:.8f}, std = {G_adj_std:.2e}")

# ============================================================
# SECTION 4: SEARCH FOR 2*PI
# ============================================================
print(f"\n{'='*70}")
print("SECTION 4: SEARCH FOR 2*pi IN SPECTRAL QUANTITIES")
print(f"{'='*70}")

target = 2 * np.pi
print(f"  Target: 2*pi = {target:.8f}\n")

# Compile all computed quantities
quantities = {
    "Tr(G)": Tr_G_exact,
    "Tr(G)/N": Tr_G_exact/N,
    "Tr(G^2)": Tr_G2_exact,
    "Tr(G^2)/N": Tr_G2_exact/N,
    "Sum_e G(ij)^2": G_edge_sum,
    "Sum_e G(ij)^2 / N_e": G_edge_sum/n_edge,
    "Bubble / N": G_bubble/N,
    "G(x,x)": G_diag[0],
    "G(x,y) [adj]": G_adj_mean,
    "N * G(x,x)": N * G_diag[0],
    "N * G(x,y)": N * G_adj_mean,
    "degree * G(x,x)": degree * G_diag[0],
    "N^2 * G(x,y)^2": N**2 * G_adj_mean**2,
    "Bubble": G_bubble,
    "Bubble/N_e": G_bubble/n_edge,
}

# Add spectral sums
spec_sums = {}
for p in [-2, -1, 0.5, 1, 2, 3]:
    val = sum(mult * lam**p for lam, mult, _ in eig_data if lam > 0)
    name = f"sum mult*lam^{p}"
    spec_sums[name] = val
    quantities[name] = val

# Add products and ratios of eigenvalues
quantities["L8/L1"] = L8/L1  # = phi^4
quantities["L6/L2"] = L6/L2  # = phi^2
quantities["L1*L8"] = L1*L8  # = 36 = b1^2
quantities["L2*L6"] = L2*L6  # = 80
quantities["(L1+L8)/2"] = (L1+L8)/2  # = 9
quantities["(L2+L6)/2"] = (L2+L6)/2  # = 10

# Galois norms
quantities["prod(lambda_broken)"] = L1*L2*L6*L8
quantities["prod/N^2"] = L1*L2*L6*L8/N**2

# Spectral zeta derived
quantities["zeta'(0) = -log det'(L)"] = -sum(mult*np.log(lam) for lam, mult, _ in eig_data if lam > 0)
quantities["det'(L)^{1/N}"] = np.exp(sum(mult*np.log(lam) for lam, mult, _ in eig_data if lam > 0) / (N-1))

# Test function related
quantities["N/(2pi)"] = N / (2*np.pi)
quantities["N^2/(2pi)"] = N**2 / (2*np.pi)
quantities["Tr(G) * 2pi / N"] = Tr_G_exact * 2*np.pi / N

# Now check which are close to 2*pi
print(f"  {'Quantity':40s} {'Value':>14s} {'Ratio to 2pi':>14s} {'Dev%':>10s}")
print(f"  {'-'*40} {'-'*14} {'-'*14} {'-'*10}")

close_matches = []
for name, val in sorted(quantities.items(), key=lambda x: abs(x[1]/target - 1) if x[1] != 0 else 999):
    if val != 0 and abs(val) > 1e-10:
        ratio = val / target
        dev = (ratio - 1) * 100
        if abs(dev) < 50:  # only show within 50%
            marker = " <---" if abs(dev) < 5 else ""
            print(f"  {name:40s} {val:14.6f} {ratio:14.6f} {dev:+10.4f}%{marker}")
            if abs(dev) < 5:
                close_matches.append((name, val, dev))

# ============================================================
# SECTION 5: SIMPLE ALGEBRAIC COMBINATIONS
# ============================================================
print(f"\n{'='*70}")
print("SECTION 5: ALGEBRAIC COMBINATIONS TARGETING 2*pi")
print(f"{'='*70}")

# Can 2*pi be expressed as a ratio/product of known spectral quantities?
# Key building blocks from the 600-cell:
# N=120, degree=12, a1=5, b1=6, phi, sqrt(5)

candidates = {
    "N*pi/a1!": N * np.pi / 120,  # = pi
    "degree * pi/a1!*2": degree * np.pi / 60,
    "N/a1!": N / 120,
    "h*pi/a1^2": h * np.pi / a1**2,
    "degree*pi/b1": degree * np.pi / b1,
    "2*degree*pi/degree": 2*np.pi,
    # Try from Tr(G)
    "Tr(G)^2 * pi": Tr_G_exact**2 * np.pi,
    "Tr(G) * Tr(G^2)": Tr_G_exact * Tr_G2_exact,
    "Tr(G) * pi / a1": Tr_G_exact * np.pi / a1,
    # From eigenvalues
    "L1+L2+L6+L8": L1+L2+L6+L8,  # = 38
    "pi*(L1+L2)": np.pi * (L1+L2),
    "sqrt(L1*L2*L6*L8)*pi/N": np.sqrt(L1*L2*L6*L8)*np.pi/N,
    # Volume of S^3/2I
    "Vol(S^3)/|2I|": 2*np.pi**2/N,
    "Vol(S^3)": 2*np.pi**2,
    "4*pi": 4*np.pi,
    # Hopf fiber
    "2*a1 * pi/a1": 2*a1 * np.pi / a1,  # = 2pi trivially
    "N_fiber_edges * theta_edge": 2*a1 * np.pi/a1,
    "N_fibers * Vol(fiber)": 72 * 2*np.pi/72,  # each fiber is C_10
    # Eta invariant candidates
    "pi^2/a1": np.pi**2/a1,
    "pi*sqrt(a1)": np.pi*np.sqrt(a1),
    # From determinant
    "exp(zeta'(0)/N)": np.exp(-sum(mult*np.log(lam) for lam, mult, _ in eig_data if lam > 0)/N),
}

print(f"\n  {'Expression':40s} {'Value':>12s} {'= 2pi?':>12s} {'Dev%':>10s}")
print(f"  {'-'*40} {'-'*12} {'-'*12} {'-'*10}")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]/target - 1)):
    ratio = val / target
    dev = (ratio - 1) * 100
    if abs(dev) < 20:
        marker = " <==" if abs(dev) < 0.01 else (" <-" if abs(dev) < 5 else "")
        print(f"  {name:40s} {val:12.6f} {ratio:12.6f} {dev:+10.4f}%{marker}")

# ============================================================
# SECTION 6: HOPF FIBER ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("SECTION 6: HOPF FIBER STRUCTURE")
print(f"{'='*70}")

# The 600-cell has 72 Hopf fibers, each a decagonal great circle (C_10).
# Each fiber has 10 vertices and 10 edges.
# Total: 72 * 10 = 720 = N_edges.  Each edge belongs to exactly 1 fiber.
# (Actually each edge might belong to 2 fibers? Let's check.)

# To find Hopf fibers: project S^3 vertices onto S^2 via Hopf map,
# then group vertices that map to the same point on S^2.

# Hopf map: (z1, z2) in S^3 subset C^2 -> [z1:z2] in CP^1 = S^2
# In coordinates: (x0,x1,x2,x3) -> (2(x0*x2+x1*x3), 2(x1*x2-x0*x3), x0^2+x1^2-x2^2-x3^2)

verts = np.array(vertices)
hopf_image = np.zeros((n_vert, 3))
for i in range(n_vert):
    x0, x1, x2, x3 = verts[i]
    hopf_image[i] = [
        2*(x0*x2 + x1*x3),
        2*(x1*x2 - x0*x3),
        x0**2 + x1**2 - x2**2 - x3**2
    ]

# Group vertices by their Hopf image (cluster by distance)
from collections import defaultdict
fiber_labels = -np.ones(n_vert, dtype=int)
n_fibers = 0
tolerance = 0.01

for i in range(n_vert):
    if fiber_labels[i] >= 0:
        continue
    fiber_labels[i] = n_fibers
    for j in range(i+1, n_vert):
        if fiber_labels[j] >= 0:
            continue
        dist = np.linalg.norm(hopf_image[i] - hopf_image[j])
        if dist < tolerance:
            fiber_labels[j] = n_fibers
    n_fibers += 1

fiber_sizes = defaultdict(int)
for fl in fiber_labels:
    fiber_sizes[fl] += 1

size_counts = defaultdict(int)
for s in fiber_sizes.values():
    size_counts[s] += 1

print(f"  Number of Hopf fibers: {n_fibers}")
print(f"  Fiber size distribution: {dict(size_counts)}")
print(f"  Total vertices: {sum(fiber_sizes.values())}")

# Check that fibers are cycles (C_10) by looking at adjacency within fibers
fiber_is_cycle = []
for f_id in range(n_fibers):
    f_verts = [i for i in range(n_vert) if fiber_labels[i] == f_id]
    if len(f_verts) != 10:
        fiber_is_cycle.append(False)
        continue
    # Count internal edges
    internal_edges = 0
    for (i, j) in edges:
        if fiber_labels[i] == f_id and fiber_labels[j] == f_id:
            internal_edges += 1
    fiber_is_cycle.append(internal_edges == 10)  # C_10 has 10 edges

n_cycles = sum(fiber_is_cycle)
print(f"  Fibers that are C_10 cycles: {n_cycles}/{n_fibers}")

# Propagator restricted to a single fiber
if n_cycles > 0:
    f_id = [i for i, c in enumerate(fiber_is_cycle) if c][0]
    f_verts = [v for v in range(n_vert) if fiber_labels[v] == f_id]

    # G restricted to fiber vertices
    G_fiber = G_mat[np.ix_(f_verts, f_verts)]

    # Internal edges of this fiber
    fiber_edges = [(i, j) for (i, j) in edges if fiber_labels[i] == f_id and fiber_labels[j] == f_id]

    # Self-energy on the fiber
    Pi_fiber = sum(G_mat[i, j]**2 for (i, j) in fiber_edges)
    Bubble_fiber = sum(G_mat[i, i]*G_mat[j, j] - G_mat[i, j]**2 for (i, j) in fiber_edges)

    print(f"\n  Single-fiber self-energy:")
    print(f"    Sum G(ij)^2 (fiber) = {Pi_fiber:.8f}")
    print(f"    Bubble (fiber)       = {Bubble_fiber:.8f}")
    print(f"    Pi_fiber * N_fibers  = {Pi_fiber * n_fibers:.8f}  (should relate to total?)")
    print(f"    Total / fiber: {G_edge_sum/n_fibers:.8f}  vs  {Pi_fiber:.8f}")

    # Phase accumulation along the fiber
    # The propagator G(i,j) for adjacent vertices on the fiber
    G_fiber_adj = [G_mat[i, j] for (i, j) in fiber_edges]
    print(f"    G(i,j) on fiber edges: {np.mean(G_fiber_adj):.8f} (mean)")

    # Holonomy: product of "transfer matrices" along the fiber
    # For a cycle C_10, the natural holonomy is the determinant of
    # the chain of propagators
    # A simpler measure: sum of G along the fiber
    G_fiber_sum = sum(G_fiber_adj)
    print(f"    Sum of G(i,j) along fiber = {G_fiber_sum:.8f}")
    print(f"    * N_fibers = {G_fiber_sum * n_fibers:.8f}")

# ============================================================
# SECTION 7: SPECTRAL ZETA AND ETA INVARIANT
# ============================================================
print(f"\n{'='*70}")
print("SECTION 7: SPECTRAL ZETA AND ETA INVARIANT")
print(f"{'='*70}")

# Spectral zeta of the Laplacian (excluding zero mode):
# zeta_L(s) = sum_{lambda>0} mult(lambda) * lambda^{-s}

# Special values:
zeta_values = {}
for s in [-2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]:
    val = sum(mult * lam**(-s) for lam, mult, _, _ in SPECTRUM_600CELL if lam > 0)
    zeta_values[s] = val

print(f"  Spectral zeta zeta_L(s):")
print(f"  {'s':>6s} {'zeta(s)':>14s} {'zeta/2pi':>14s}")
for s, val in sorted(zeta_values.items()):
    print(f"  {s:>6.1f} {val:>14.6f} {val/(2*np.pi):>14.6f}")

# Derivative at s=0 (regularized determinant)
# zeta'(0) = -sum mult*log(lambda)
zeta_prime_0 = -sum(mult * np.log(lam) for lam, mult, _, _ in SPECTRUM_600CELL if lam > 0)
log_det = -zeta_prime_0  # log det'(L)
print(f"\n  zeta'(0) = {zeta_prime_0:.6f}")
print(f"  log det'(L) = {log_det:.6f}")
print(f"  det'(L) = exp({log_det:.2f})")

# Eta invariant: needs DIRAC operator, not Laplacian.
# On S^3/2I, the Dirac eigenvalues are +/-(k+3/2) with multiplicities
# from McKay recursion.
# eta(s) = sum sign(lambda)*|lambda|^{-s}
# For a spherical space form S^3/Gamma:
# eta(S^3/Gamma) = (1/|Gamma|) * sum_{g!=1} csc^2(theta_g/2) * cot(theta_g/2) / 4
# where theta_g are the rotation angles of elements g in SU(2).

# Conjugacy classes of 2I:
# Class | Size | Order | theta (SU(2) angle)
# C0    |  1   |  1    | 0
# C1    |  1   |  2    | 2*pi (= -I, eigenvalues (-1,-1))
# C2    | 12   |  10   | pi/5
# C3    | 12   |  10   | 3*pi/5
# C4    | 20   |  6    | pi/3
# C5    | 20   |  6    | 2*pi/3
# C6    | 30   |  4    | pi/2
# C7    | 12   |  5    | 2*pi/5
# C8    | 12   |  5    | 4*pi/5

conj_classes = [
    # (size, order, theta, name)
    (1,   1,  0,            "identity"),
    (1,   2,  np.pi,        "-I"),       # theta = pi for SU(2) angle (eigenvalues -1,-1)
    (12, 10,  np.pi/5,      "C_10 (pi/5)"),
    (12, 10,  3*np.pi/5,    "C_10 (3pi/5)"),
    (20,  6,  np.pi/3,      "C_6 (pi/3)"),
    (20,  6,  2*np.pi/3,    "C_6 (2pi/3)"),
    (30,  4,  np.pi/2,      "C_4 (pi/2)"),
    (12,  5,  2*np.pi/5,    "C_5 (2pi/5)"),
    (12,  5,  4*np.pi/5,    "C_5 (4pi/5)"),
]

# Eta invariant formula (Donnelly / Bar):
# eta(S^3/Gamma) = (1/|Gamma|) * sum_{g!=e} size(C) * h(theta)
# where h(theta) = csc^2(theta/2) * cot(theta/2) / 4
# But wait, for SU(2), the element g has eigenvalues exp(+/- i*theta/2),
# and the formula uses theta as the full rotation angle.
# Let me use the Donnelly formula carefully.

# For an element g in SU(2) with eigenvalues exp(+/-i*theta/2):
# f(g) = (1/4) * cot(theta/2) * csc^2(theta/2)   [Donnelly 1978]
# eta(S^3/Gamma) = -(1/|Gamma|) * sum_{g!=1} f(g)

# But actually for the Poincare homology sphere S^3/2I,
# the Dirac operator has a specific eta. Let me compute:

eta_sum = 0
print(f"\n  Eta invariant of S^3/2I:")
print(f"  {'Class':15s} {'Size':>5s} {'theta':>8s} {'f(g)':>14s} {'size*f(g)':>14s}")

for size, order, theta, name in conj_classes:
    if theta < 1e-10 or abs(theta - np.pi) < 1e-10:
        # Identity: theta=0, f is infinite (skip, it's excluded)
        # -I: theta=pi, cot(pi/2)=0, so f=0
        if abs(theta - np.pi) < 1e-10:
            f_g = 0  # cot(pi/2) = 0
            eta_sum += size * f_g
            print(f"  {name:15s} {size:5d} {theta:8.4f} {f_g:14.6f} {size*f_g:14.6f}")
        else:
            print(f"  {name:15s} {size:5d} {theta:8.4f} {'(excluded)':>14s} {'---':>14s}")
        continue

    half_theta = theta / 2
    f_g = 0.25 * np.cos(half_theta) / np.sin(half_theta)**3  # cot * csc^2 / 4
    eta_sum += size * f_g
    print(f"  {name:15s} {size:5d} {theta:8.4f} {f_g:14.6f} {size*f_g:14.6f}")

eta_invariant = -eta_sum / N
print(f"\n  eta(S^3/2I) = -{eta_sum:.6f} / {N} = {eta_invariant:.8f}")
print(f"  eta / (2*pi) = {eta_invariant / (2*np.pi):.8f}")
print(f"  2*pi / eta = {2*np.pi / eta_invariant:.8f}" if abs(eta_invariant) > 1e-10 else "  eta = 0")
print(f"  pi * eta = {np.pi * eta_invariant:.8f}")

# ============================================================
# SECTION 8: SELF-ENERGY ON THE GRAPH (RIGOROUS)
# ============================================================
print(f"\n{'='*70}")
print("SECTION 8: LATTICE SELF-ENERGY (ONE-LOOP)")
print(f"{'='*70}")

# The one-loop vacuum polarization on the graph.
# In lattice gauge theory, the polarization tensor is:
#   Pi_e = sum_loops_through_e (product of propagators around loop)
#
# The simplest (scalar) self-energy is:
#   delta(1/g^2) = (g^2/N) * sum_{lambda!=0} f(lambda) / lambda
#
# where f(lambda) depends on the vertex structure.
# For standard lattice: f(lambda) = lambda itself (from the vertex factor).
# => delta(1/g^2) = (g^2/N) * sum mult*1 = g^2 * (N-1)/N ~ g^2
#    (this is the one-loop trace, trivial)
#
# The non-trivial structure comes from the HOPF FIBER restriction.
# The gauge field lives on the FIBER, so the self-energy is:
#   delta(1/g^2) = integral_around_fiber = phase_accumulation
#
# For a U(1) gauge field on a fiber C_10:
# The Wilson loop around the fiber = exp(i * sum_edges A_edge)
# At one loop: <W> ~ exp(-g^2 * Pi_fiber)
#
# The key identity from lattice QED (Wilson 1974):
# For a plaquette (elementary loop), the one-loop correction is:
#   delta(1/g^2) = (1/2) * sum_p [1 - cos(theta_p)] / N_plaq
#
# For the Hopf fiber (a single large loop):
#   delta(1/g^2) = (1/2) * [sum of angles around loop]
#
# If each edge carries angle pi/a1 = pi/5:
#   delta(1/g^2) = (1/2) * [10 * pi/5] = (1/2) * 2*pi = pi

# But we need delta(1/alpha) = -2*pi*alpha, which has the coupling alpha in it.
# This is like a "dressed" self-energy: Pi = 2*pi * alpha

# Let's be more careful. The general structure is:
# 1/alpha_ren = 1/alpha_bare - Pi(alpha)
# At one loop: Pi = c * alpha  where c is the geometric coefficient.
# We need c = 2*pi.

# FROM THE HOPF FIBER:
# The fiber C_10 has Laplacian eigenvalues: lambda_k = 2*(1-cos(2*pi*k/10)) for k=0,...,9
# The spectral gap: lambda_1 = 2*(1-cos(pi/5)) = 2*(1-phi/2) = 2-phi = 1/phi^2 = L1 (the 600-cell eigenvalue!)
# Total holonomy = 10 * pi/5 = 2*pi (from the lattice dispersion relation)

# One-loop self-energy on C_10:
# Pi_fiber = sum_{k=1}^{9} 1/lambda_k(C_10)
fiber_eigenvalues_C10 = [2*(1-np.cos(2*np.pi*k/10)) for k in range(10)]
Pi_C10 = sum(1/lam for lam in fiber_eigenvalues_C10 if lam > 1e-10)
print(f"  C_10 (Hopf fiber) eigenvalues: {[f'{l:.4f}' for l in sorted(fiber_eigenvalues_C10)]}")
print(f"  Sum 1/lambda_k (C_10, k>=1) = {Pi_C10:.8f}")
print(f"  Pi_C10 / (2*pi) = {Pi_C10 / (2*np.pi):.8f}")
print(f"  Pi_C10 / pi = {Pi_C10 / np.pi:.8f}")
print(f"  2*Pi_C10 = {2*Pi_C10:.8f}")
print(f"  Pi_C10 * pi/a1 = {Pi_C10 * np.pi/a1:.8f}")

# The C_10 zeta at s=1: sum 1/lambda = ?
# For C_n: lambda_k = 2*(1-cos(2*pi*k/n)), k=0,...,n-1
# sum_{k=1}^{n-1} 1/(2*(1-cos(2*pi*k/n))) = (n^2-1)/12
# For n=10: (100-1)/12 = 99/12 = 33/4 = 8.25
print(f"\n  Analytic: sum_C10 1/lambda = (n^2-1)/12 = {(10**2-1)/12:.4f}")
print(f"  Numeric: {Pi_C10:.4f}")
print(f"  Match: {'YES' if abs(Pi_C10 - 99/12) < 1e-6 else 'NO'}")

# So Pi_C10 = (n^2-1)/12 = (4*a1^2-1)/12 = 99/12 = 33/4
# For the 600-cell: n = 2*a1 = 10 (Hopf fiber has 2*a1 vertices)
# Pi = ((2*a1)^2-1)/12 = (4*a1^2-1)/12
Pi_analytic = (4*a1**2 - 1) / 12
print(f"  Pi_C10 = (4*a1^2-1)/12 = {Pi_analytic:.6f}")

# Now: can we get 2*pi from this?
# Pi_C10 = 33/4 = 8.25
# 2*pi = 6.283...
# Ratio: 2*pi / Pi_C10 = 0.7616...

# What if the self-energy involves the TOTAL over all fibers normalized?
# Total = 72 * Pi_C10 = 72 * 33/4 = 72 * 8.25 = 594
# 594 / N = 594/120 = 4.95 = 99/20
# 594 / (N*degree) = 594/1440 = 0.4125 = 33/80
# 594 / (2*pi) = 94.5...  nope

# What about the ANGULAR self-energy?
# Instead of 1/lambda, use the angle: theta_k = arccos(1-lambda_k/2)
# theta_k = 2*pi*k/n for C_n

angles_C10 = [2*np.pi*k/10 for k in range(10)]
angular_sum = sum(angles_C10[1:])  # exclude k=0
print(f"\n  Angular self-energy (Hopf fiber):")
print(f"    Sum of angles = {angular_sum:.6f}")
print(f"    = 2*pi * (0+1+2+...+9)/10 = 2*pi * 45/10 = 2*pi * 4.5 = {2*np.pi*4.5:.6f}")
print(f"    = 9*pi = {9*np.pi:.6f}")
print(f"    Ratio to 2*pi: {angular_sum/(2*np.pi):.4f}")

# Hmm. Sum of angles is 9*pi, not 2*pi.
# The 2*pi comes from the FIRST harmonic only: theta_1 = 2*pi/10 = pi/5
# Multiplied by the number of edges: 10 * pi/5 = 2*pi
# This is the HOLONOMY, not the self-energy sum.

holonomy = 10 * np.pi / a1  # 10 edges * pi/5 per edge
print(f"\n  Holonomy = 2*a1 * pi/a1 = {holonomy:.6f} = 2*pi? {abs(holonomy-2*np.pi)<1e-10}")

# ============================================================
# SECTION 9: CONNECTION BETWEEN SELF-ENERGY AND HOLONOMY
# ============================================================
print(f"\n{'='*70}")
print("SECTION 9: SELF-ENERGY vs HOLONOMY")
print(f"{'='*70}")

# The question: WHY does the holonomy (2*pi) enter the coupling equation?
#
# In Kaluza-Klein theory:
# The 4D gauge coupling g^2 = 1/(R * M_KK^{d_int})
# where R is the radius of the internal space.
# For a U(1) fiber with periodicity 2*pi*R:
# The Vieta product alpha*alpha' = 1/Vol(fiber) = 1/(2*pi*R)
# If R=1 (natural units): alpha*alpha' = 1/(2*pi)
#
# This is EXACTLY the Vieta product of the alpha equation:
# alpha * alpha' = 1/(2*pi)  [constant term / quadratic coefficient]

print(f"  Vieta product from alpha equation:")
print(f"    alpha * alpha' = 1/(2*pi) = {1/(2*np.pi):.8f}")
print(f"    Check: {alpha_em:.8f} * {1/(2*np.pi*alpha_em):.8f} = {1/(2*np.pi):.8f}")

# The KK interpretation:
# The fiber is the Hopf fiber C_10 with circumference = 10 edges * (pi/a1 per edge)
# In lattice units: C = 10 * lattice_spacing * pi/a1
# If lattice_spacing = 1/a1: C = 10/(a1) * pi/a1 = 10*pi/25 = 2*pi/5
# Hmm, that gives 2*pi/5, not 2*pi.
# If lattice_spacing = 1: C = 10 * pi/5 = 2*pi. YES!

# So the KK picture is:
# 1. The Hopf fiber is a circle of circumference 2*pi (in natural lattice units)
# 2. KK gives: the quadratic coefficient = Vol(fiber) = 2*pi
# 3. The linear coefficient = 4*a1*phi^4 comes from the spectral action
# 4. The constant = 1 is normalization

# Is this a DERIVATION?
# The argument is:
# (a) The Hopf fiber IS topologically S^1 (a circle)  — FACT
# (b) Its holonomy in the 600-cell IS 2*pi  — DERIVED (from spectral gap)
# (c) In KK theory, the coupling constant involves 1/Vol(fiber)  — STANDARD KK
# (d) Therefore the quadratic coefficient IS 2*pi  — DERIVED (from a+b+c)

# The WEAK point is (c): why should we USE the KK identification?
# This is the same as asking: why is alpha the coupling of a U(1) gauge
# field that lives on the Hopf fiber?

# Can we prove (c) from the spectral action alone?
# In Connes NCG: the spectral action gives
# S = Tr(f(D^2/Lambda^2))
# = f_0 * Lambda^4 * a_0 + f_2 * Lambda^2 * a_2 + f(0) * a_4 + ...
# The gauge coupling comes from the f(0) * a_4 term.
# a_4 contains the Yang-Mills action: (1/g^2) * F^2

# For the test function f with moments f_k = N^k/(2*pi) (from exp466):
# f(0) = ? (value of f at the origin)
# If f(x) = (N/2pi) * delta(x) [concentrated at 0], then f(0) = N/(2*pi)
# But f_0 = int f(v)/sqrt(v) dv = N^0/(2*pi) = 1/(2*pi)

# The gauge coupling: 1/g_i^2 = f(0) * C_i
# If the fiber structure gives C_EM = 1 (trivially, since U(1) has no
# non-abelian structure), then:
# 1/alpha = f(0) = ... involves 2*pi in the denominator

# BUT: this gives 1/alpha = something/(2*pi), not the QUADRATIC 2*pi*alpha^2.
# The quadratic form comes from SOLVING for alpha when you have both
# the bare coupling AND the one-loop correction.

print(f"\n  KK interpretation:")
print(f"    Vol(Hopf fiber) = 2*pi (from holonomy)")
print(f"    KK: alpha*alpha' = 1/Vol = 1/(2*pi)")
print(f"    This IS the Vieta product of the quadratic.")
print(f"    => 2*pi is the VOLUME of the Hopf fiber in the KK picture.")

# ============================================================
# SECTION 10: LATTICE SELF-ENERGY FROM REPRESENTATION THEORY
# ============================================================
print(f"\n{'='*70}")
print("SECTION 10: REPRESENTATION-THEORETIC SELF-ENERGY")
print(f"{'='*70}")

# On the 600-cell (Cayley graph of 2I), the propagator in the irrep basis:
# G_rho = 1/lambda_rho
# The "self-energy" (one-loop correction) at zero momentum:
# Pi(0) = (1/N) * sum_rho dim(rho)^2 / lambda_rho^2
#        = Tr(G^2) / N

Pi_rep = Tr_G2_exact / N
print(f"  Pi(0) = Tr(G^2)/N = {Pi_rep:.8f}")
print(f"  Pi(0) * N = Tr(G^2) = {Tr_G2_exact:.8f}")

# The self-energy correction: delta(1/alpha) = -Pi * alpha
# => 1/alpha = C - Pi*alpha
# => Pi*alpha^2 - C*alpha + 1 = 0  (comparing with 2*pi*alpha^2 - C*alpha + 1 = 0)
# => We need Pi = 2*pi iff Tr(G^2)/N = 2*pi

print(f"\n  For consistency with alpha equation:")
print(f"    Need: Pi(0) = 2*pi = {2*np.pi:.6f}")
print(f"    Have: Pi(0) = Tr(G^2)/N = {Pi_rep:.6f}")
print(f"    Ratio: {Pi_rep/(2*np.pi):.6f}")
print(f"    THIS DOES NOT MATCH.")

# What if the self-energy involves a different power or normalization?
# Try: (1/N) * sum dim^2 * lambda^{-p} for various p
print(f"\n  Tr(G^p)/N for various p:")
for p in [0.5, 1, 1.5, 2, 2.5, 3]:
    val = sum(mult * lam**(-p) for lam, mult, _, _ in SPECTRUM_600CELL if lam > 0) / N
    print(f"    p={p:.1f}: {val:.8f}   /2pi = {val/(2*np.pi):.6f}")

# What normalization of Tr(G^p) gives 2*pi?
for p in [0.5, 1, 1.5, 2, 2.5, 3]:
    val = sum(mult * lam**(-p) for lam, mult, _, _ in SPECTRUM_600CELL if lam > 0)
    need = 2*np.pi
    factor = need / val
    print(f"    For Tr(G^{p:.1f})*{factor:.6f} = 2*pi  (factor = {factor:.6f})")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  QUESTION: Can 2*pi be derived from a lattice self-energy on the 600-cell?

  RESULTS:

  1. SPECTRAL QUANTITIES: No combination of Tr(G^p), spectral zeta values,
     or graph-theoretic sums naturally produces 2*pi.
     - Tr(G^2)/N = {Pi_rep:.4f} (not 2*pi = {2*np.pi:.4f})
     - Pi_C10 = (4a1^2-1)/12 = {Pi_analytic} (not 2*pi either)

  2. ETA INVARIANT: eta(S^3/2I) = {eta_invariant:.6f}
     - Not directly 2*pi.

  3. HOPF HOLONOMY: 2*a1 * pi/a1 = 2*pi.  EXACT.
     - This IS the total phase around a Hopf fiber.
     - It IS derived from the spectral gap of C_10.

  4. KK IDENTIFICATION: alpha*alpha' = 1/(2*pi) = 1/Vol(fiber).
     - This identifies 2*pi as the VOLUME of the Hopf fiber.
     - The Vieta product of the quadratic IS 1/(2*pi).

  5. THE GAP: The derivation chain is
       spectral gap of C_10  -->  angle pi/a1 per edge  -->  holonomy 2*pi
       +  KK identification  -->  2*pi enters as fiber volume
     The weak point is the KK IDENTIFICATION itself.
     It is NOT derived from a lattice propagator calculation.

  HONEST ASSESSMENT:
    - The 2*pi is geometrically NATURAL (Hopf holonomy = fiber volume).
    - The identification with the coupling equation coefficient is
      MOTIVATED by Kaluza-Klein theory.
    - A direct lattice self-energy Pi = 2*pi has NOT been found.
    - The eta invariant does NOT produce 2*pi directly.

  STATUS: STRUCTURAL (KK-motivated geometric identification),
  NOT DERIVED (from first-principles lattice calculation).
""")

print("=" * 70)
print("EXP-530 COMPLETE")
print("=" * 70)
