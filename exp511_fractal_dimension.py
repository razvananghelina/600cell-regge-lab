"""
EXP-511: Spectral Dimension of Icosahedral Fractal Structures
===============================================================
QUESTION: Is the spectral dimension of an icosahedral quasicrystal
or fractal equal to phi = (1+sqrt(5))/2?

If yes: the CC formula Lambda = 4*sqrt(5) * |E_vac|^phi * exp(-N*L_1)
has a PHYSICAL justification (vacuum lives on a fractal of dimension phi).

APPROACH:
1. Build a Fibonacci chain (1D quasicrystal) and measure d_s
2. Build a Penrose-like tiling (2D) and measure d_s
3. Build a fractal from iterating the McKay/Galois map on the 600-cell
4. Compute spectral dimension from random walk return probability:
   P(t) ~ t^{-d_s/2}
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '.')
from commons import (build_600cell, PHI, PHI_CONJ, SQRT5, a1, b1, N,
                      alpha_em, ln_inv_alpha, degree)

print("=" * 70)
print("EXP-511: SPECTRAL DIMENSION OF ICOSAHEDRAL FRACTALS")
print("=" * 70)

# ============================================================
# PART 1: Fibonacci Chain (1D Quasicrystal)
# ============================================================
print(f"\n{'='*70}")
print("PART 1: FIBONACCI CHAIN")
print(f"{'='*70}")

# The Fibonacci chain is a 1D quasicrystal with two intervals:
# L (long, length phi) and S (short, length 1).
# Substitution rule: L -> LS, S -> L
# This generates a quasiperiodic sequence with phi-scaling.

def fibonacci_chain(n_generations):
    """Generate Fibonacci chain by substitution."""
    chain = [1]  # start with L
    for _ in range(n_generations):
        new_chain = []
        for c in chain:
            if c == 1:  # L -> LS
                new_chain.extend([1, 0])
            else:  # S -> L
                new_chain.extend([1])
        chain = new_chain
    return chain

# Build the chain and its adjacency graph
n_gen = 12
chain = fibonacci_chain(n_gen)
n_sites = len(chain)
print(f"\n  Fibonacci chain: {n_gen} generations, {n_sites} sites")

# Build adjacency: each site connects to its neighbors
# Hopping strengths: t_L = 1 (long bond), t_S = 1/phi (short bond)
# Actually, for spectral dimension measurement, use uniform hopping
A_fib = np.zeros((n_sites, n_sites))
for i in range(n_sites - 1):
    A_fib[i, i+1] = 1
    A_fib[i+1, i] = 1

# Laplacian
degrees_fib = A_fib.sum(axis=1)
L_fib = np.diag(degrees_fib) - A_fib

# Spectral dimension from heat kernel: K(t) = Tr(exp(-t*L)) / N
# P(t) = K(t)/N ~ t^{-d_s/2} for intermediate t
eigvals_fib = eigh(L_fib, eigvals_only=True)
eigvals_fib = np.sort(eigvals_fib)

# Compute heat trace for various t
t_values = np.logspace(-1, 3, 200)
P_t_fib = np.zeros_like(t_values)
for i, t in enumerate(t_values):
    P_t_fib[i] = np.sum(np.exp(-t * eigvals_fib)) / n_sites

# Spectral dimension from log-log slope: d_s = -2 * d(log P)/d(log t)
log_t = np.log(t_values)
log_P = np.log(P_t_fib)

# Use intermediate range to estimate d_s
# Avoid t too small (lattice effects) and t too large (finite size)
mask = (t_values > 1) & (t_values < n_sites / 10)
if np.sum(mask) > 10:
    coeffs = np.polyfit(log_t[mask], log_P[mask], 1)
    d_s_fib = -2 * coeffs[0]
    print(f"  Spectral dimension (from slope): d_s = {d_s_fib:.4f}")
    print(f"  phi = {PHI:.4f}")
    print(f"  Match: {abs(d_s_fib - PHI) < 0.1}")
else:
    d_s_fib = None
    print(f"  Insufficient data for slope estimate")

# Also compute running d_s(t)
print(f"\n  Running spectral dimension d_s(t):")
print(f"  {'t':>10} {'P(t)':>14} {'d_s(t)':>10}")
for i in range(1, len(t_values)-1):
    if t_values[i] > 0.5 and t_values[i] < n_sites / 5:
        ds_local = -2 * (log_P[i+1] - log_P[i-1]) / (log_t[i+1] - log_t[i-1])
        if i % 15 == 0:
            print(f"  {t_values[i]:10.2f} {P_t_fib[i]:14.6e} {ds_local:10.4f}")

# ============================================================
# PART 2: 600-Cell as Fractal Seed
# ============================================================
print(f"\n{'='*70}")
print("PART 2: 600-CELL SPECTRAL DIMENSION")
print(f"{'='*70}")

# The 600-cell is NOT a fractal, but we can measure its spectral dimension
verts, adj, lap = build_600cell()
eigvals_600 = np.sort(eigh(lap, eigvals_only=True))

# Heat trace
P_t_600 = np.zeros_like(t_values)
for i, t in enumerate(t_values):
    P_t_600[i] = np.sum(np.exp(-t * eigvals_600)) / N

log_P_600 = np.log(np.maximum(P_t_600, 1e-300))

print(f"\n  Running spectral dimension d_s(t) for 600-cell:")
print(f"  {'t':>10} {'P(t)':>14} {'d_s(t)':>10}")
ds_values_600 = []
for i in range(1, len(t_values)-1):
    if P_t_600[i] > 1e-200 and P_t_600[i] < 0.99:
        ds_local = -2 * (log_P_600[i+1] - log_P_600[i-1]) / (log_t[i+1] - log_t[i-1])
        ds_values_600.append((t_values[i], ds_local))
        if i % 10 == 0:
            print(f"  {t_values[i]:10.4f} {P_t_600[i]:14.6e} {ds_local:10.4f}")

# Find plateau value
if ds_values_600:
    t_arr = np.array([x[0] for x in ds_values_600])
    ds_arr = np.array([x[1] for x in ds_values_600])
    # Look for plateau in the middle range
    mid_mask = (t_arr > 0.1) & (t_arr < 2)
    if np.sum(mid_mask) > 3:
        ds_plateau = np.mean(ds_arr[mid_mask])
        print(f"\n  Plateau d_s (0.1 < t < 2): {ds_plateau:.4f}")
        print(f"  phi = {PHI:.4f}")
        print(f"  4 (topological dim) = 4")
    else:
        # Try wider range
        mid_mask = (t_arr > 0.05) & (t_arr < 5)
        if np.sum(mid_mask) > 3:
            ds_plateau = np.mean(ds_arr[mid_mask])
            print(f"\n  Plateau d_s (0.05 < t < 5): {ds_plateau:.4f}")

# ============================================================
# PART 3: Density of States Analysis
# ============================================================
print(f"\n{'='*70}")
print("PART 3: DENSITY OF STATES")
print(f"{'='*70}")

# For a d-dimensional system: rho(lambda) ~ lambda^{d/2 - 1} near lambda = 0
# This gives an alternative definition of spectral dimension.

# For the 600-cell: discrete spectrum, so we use the integrated DOS
# N(lambda) = number of eigenvalues <= lambda ~ lambda^{d_s/2}

# Eigenvalues of 600-cell Laplacian
print(f"\n  600-cell Laplacian eigenvalues (9 distinct):")
distinct_600 = []
for ev in eigvals_600:
    if not any(abs(ev - d) < 0.001 for d in distinct_600):
        distinct_600.append(ev)

for i, ev in enumerate(distinct_600):
    mult = sum(1 for x in eigvals_600 if abs(x - ev) < 0.001)
    cumul = sum(1 for x in eigvals_600 if x <= ev + 0.001)
    print(f"    lambda_{i} = {ev:8.4f}, mult = {mult:3d}, N(lambda) = {cumul:3d}")

# N(lambda) vs lambda for the nonzero eigenvalues
nonzero_evals = sorted([ev for ev in eigvals_600 if ev > 0.01])
lambdas = np.unique(np.round(nonzero_evals, 3))
N_lambda = [sum(1 for x in nonzero_evals if x <= lam + 0.001) for lam in lambdas]

# Fit N(lambda) ~ lambda^{d_s/2} in log-log
log_lam = np.log(lambdas)
log_N = np.log(np.array(N_lambda, dtype=float))
coeffs_dos = np.polyfit(log_lam, log_N, 1)
d_s_dos = 2 * coeffs_dos[0]

print(f"\n  N(lambda) ~ lambda^{coeffs_dos[0]:.4f}")
print(f"  => d_s = 2 * {coeffs_dos[0]:.4f} = {d_s_dos:.4f}")
print(f"  Compare: phi = {PHI:.4f}, 4.0, 3.0")

# ============================================================
# PART 4: Iterated Galois Map (Fractal Construction)
# ============================================================
print(f"\n{'='*70}")
print("PART 4: ITERATED GALOIS MAP")
print(f"{'='*70}")

# Construct a fractal by iterating the Galois map on the 600-cell.
# The Galois map sigma: phi -> phi' acts on the eigenvalues:
# L_k -> sigma(L_k) for Galois-broken eigenvalues.

# The "fractal" arises from the RENORMALIZATION of the heat kernel:
# At scale t: d_s(t) varies. The running d_s might approach phi.

# For the 600-cell, the spectral zeta function is:
# zeta(s) = sum_k lambda_k^{-s} (for lambda_k > 0)
# This encodes the spectral dimension via the pole structure.

print(f"\n  Spectral zeta function of 600-cell:")
# Compute zeta(s) for various s
s_values = np.linspace(0.5, 5, 20)
zeta_vals = []
nonzero = [ev for ev in eigvals_600 if ev > 0.001]
for s in s_values:
    z = sum(ev**(-s) for ev in nonzero)
    zeta_vals.append(z)
    if abs(s - 1) < 0.1 or abs(s - 2) < 0.1 or abs(s - PHI) < 0.1:
        print(f"    zeta({s:.3f}) = {z:.6f}")

# Check: does zeta have a pole at s = d_s/2?
# For a smooth d-dim manifold: zeta has a pole at s = d/2.
# For a fractal: the pole position gives d_s/2.

# Approximate: find where zeta(s) is largest / diverges
print(f"\n  Zeta values near potential poles:")
for i, (s, z) in enumerate(zip(s_values, zeta_vals)):
    print(f"    s = {s:.3f}: zeta = {z:.4f}")

# ============================================================
# PART 5: Walk Dimension and Fractal Dimension
# ============================================================
print(f"\n{'='*70}")
print("PART 5: FRACTAL ANALYSIS OF THE 600-CELL")
print(f"{'='*70}")

# For a fractal: d_s = 2*d_H / d_w
# where d_H = Hausdorff dimension, d_w = walk dimension

# For the 600-cell (a graph):
# d_H = topological dimension = 3 (S^3 is a 3-manifold)
# d_w = walk dimension = 2 (for a "normal" graph, random walk ~ sqrt(t))
# d_s = 2*3/2 = 3 (expected for S^3)

# But if the structure is FRACTAL (quasicrystal):
# d_w > 2 (anomalous diffusion, sub-diffusive)
# d_s < d_H (spectral dimension less than Hausdorff)

# Check: mean squared displacement of random walk on 600-cell
n_walks = 1000
n_steps = 500

# Random walk on the 600-cell graph
msd = np.zeros(n_steps)

for walk in range(n_walks):
    pos = 0  # start at vertex 0
    positions = [verts[pos]]

    for step in range(n_steps):
        # Find neighbors
        neighbors = np.where(adj[pos] > 0)[0]
        # Choose random neighbor
        pos = np.random.choice(neighbors)
        positions.append(verts[pos])

        # Distance from origin (on S^3: geodesic distance)
        dot = np.clip(np.dot(verts[0], verts[pos]), -1, 1)
        geo_dist = np.arccos(dot)
        msd[step] += geo_dist**2

msd /= n_walks

# Fit MSD ~ t^{2/d_w}
# Use log-log fit in intermediate range
t_walk = np.arange(1, n_steps + 1)
log_t_walk = np.log(t_walk)
log_msd = np.log(np.maximum(msd, 1e-10))

# MSD saturates on a finite graph, so use early-to-mid range
fit_mask = (t_walk > 2) & (t_walk < 50)
if np.sum(fit_mask) > 5:
    coeffs_msd = np.polyfit(log_t_walk[fit_mask], log_msd[fit_mask], 1)
    exponent = coeffs_msd[0]
    d_w = 2 / exponent
    print(f"\n  MSD ~ t^{exponent:.4f}")
    print(f"  Walk dimension d_w = 2/{exponent:.4f} = {d_w:.4f}")
    print(f"  Normal diffusion: d_w = 2.0, exponent = 1.0")
    print(f"  Sub-diffusive: d_w > 2, exponent < 1")

    # Spectral dimension from d_s = 2*d_H/d_w
    d_H = 3  # S^3
    d_s_from_walk = 2 * d_H / d_w
    print(f"\n  d_s = 2*d_H/d_w = 2*{d_H}/{d_w:.4f} = {d_s_from_walk:.4f}")
    print(f"  phi = {PHI:.4f}")

# ============================================================
# PART 6: The Key Test - Penrose Quasicrystal
# ============================================================
print(f"\n{'='*70}")
print("PART 6: PENROSE-TYPE ANALYSIS")
print(f"{'='*70}")

# The 600-cell can be seen as a PROJECTION from 8D (E8 lattice) to 4D.
# This projection is analogous to the projection that creates
# Penrose tilings (from 5D to 2D).

# In the cut-and-project method:
# E8 lattice in R^8 -> project to R^4 (physical) = 600-cell
#                    -> project to R^4 (perpendicular) = "dark" space

# The perpendicular-space image determines the "acceptance window"
# and creates the quasicrystalline order.

# For a Penrose tiling (2D projection from 5D):
# d_H = 2 (embedding dimension)
# Spectral dimension d_s depends on the vertex connectivity

# For the 600-cell viewed as a 4D quasicrystal:
# The SCALING of the quasicrystal is governed by phi.
# The inflation factor is phi^2 (area scaling) or phi (linear scaling).

# Key property: in a Penrose tiling, the number of vertices
# scales as phi^{2n} when you inflate by phi^n.
# The "spectral" scaling involves the ratio log(N)/log(R).

# For the 600-cell "inflated" by factor phi:
# The number of vertices scales as phi^d_H where d_H is the Hausdorff dim.
# For a 3-manifold (S^3): phi^3 scaling for volume.
# But for a QUASICRYSTAL: the scaling might differ.

print(f"""
  The 600-cell as a CUT-AND-PROJECT quasicrystal:
    Physical space: R^4 (H4 subspace of R^8)
    Perpendicular space: R^4 (Galois conjugate)
    Lattice: E8 root lattice in R^8
    Projection: H4 eigenspace of Coxeter element

  In this picture:
    Inflation symmetry: scaling by phi (from Coxeter eigenvalue)
    The quasicrystal has phi-scaling built in.

  For a quasicrystal with inflation factor tau = phi:
    If N(R) ~ R^d_H (number of points in ball of radius R):
      d_H = d (embedding dimension) for a quasicrystal
      d_H = 3 (for S^3 topology)

    The SPECTRAL dimension d_s depends on the connectivity.
    For a regular lattice: d_s = d_H = d.
    For a quasicrystal: d_s can differ due to critical states.
""")

# ============================================================
# PART 7: The Exact Spectral Dimension
# ============================================================
print(f"{'='*70}")
print("PART 7: EXACT SPECTRAL DIMENSION FROM EIGENVALUE SCALING")
print(f"{'='*70}")

# For the 600-cell: the 9 eigenvalues are in Z[phi].
# The Galois pairs: (12-6phi, 6+6phi) and (12-2sqrt5, 12+2sqrt5)
# The Galois map scales eigenvalues by phi in some sense:
# L_1 = 12-6*phi, L_8 = 6+6*phi
# L_8/L_1 = (6+6*phi)/(12-6*phi) = 6*(1+phi)/(6*(2-phi))
# = (1+phi)/(2-phi) = phi^2/phi'^2 = (phi/phi')^2

ratio_18 = (6 + 6*PHI) / (12 - 6*PHI)
print(f"\n  L_8/L_1 = {ratio_18:.6f}")
print(f"  phi^2/phi'^2 = {(PHI/abs(PHI_CONJ))**2:.6f}")
print(f"  = phi^4 = {PHI**4:.6f}")
print(f"  Match: {abs(ratio_18 - PHI**4) < 0.001}")

# L_8/L_1 = phi^4! This is the Galois scaling factor.
print(f"\n  GALOIS SCALING: L_8 = L_1 * phi^4")
print(f"    L_1 = {12-6*PHI:.6f}")
print(f"    L_1 * phi^4 = {(12-6*PHI)*PHI**4:.6f}")
print(f"    L_8 = {6+6*PHI:.6f}")
print(f"    Match: {abs((12-6*PHI)*PHI**4 - (6+6*PHI)) < 0.001}")

# For the second pair:
ratio_26 = (12 + 2*SQRT5) / (12 - 2*SQRT5)
print(f"\n  L_6/L_2 = {ratio_26:.6f}")
print(f"  = {ratio_26:.6f}")
# (12+2sqrt5)/(12-2sqrt5) = ?
# Let me compute: 12+2*2.236 = 16.472, 12-2*2.236 = 7.528
# 16.472/7.528 = 2.188. Not a clean phi power.
print(f"  phi^2 = {PHI**2:.6f}")
print(f"  phi^{np.log(ratio_26)/np.log(PHI):.4f} = {ratio_26:.6f}")

galois_exp_26 = np.log(ratio_26) / np.log(PHI)
print(f"  Galois scaling for pair 2: phi^{galois_exp_26:.4f}")

# ============================================================
# PART 8: SPECTRAL DIMENSION FROM GALOIS SCALING
# ============================================================
print(f"\n{'='*70}")
print("PART 8: SPECTRAL DIMENSION FROM GALOIS SCALING")
print(f"{'='*70}")

# The Galois map sigma acts on eigenvalues:
# L_1 -> L_8 = L_1 * phi^4
# L_2 -> L_6 = L_2 * phi^{1.62}

# In a fractal with spectral dimension d_s:
# The eigenvalue ratio under scaling by factor s is:
# lambda(sR) / lambda(R) = s^{2/d_s}

# For our Galois scaling (which IS a "reflection" that scales by phi):
# L_8/L_1 = phi^4 = phi^{2/d_s * ?}

# If the Galois "scale factor" is phi^2 (area scaling):
# phi^4 = (phi^2)^{2/d_s}
# 4 = 2 * 2/d_s
# d_s = 1. Not phi.

# If the Galois "scale factor" is phi (linear scaling):
# phi^4 = phi^{2/d_s}
# 4 = 2/d_s
# d_s = 1/2. Not phi.

# Hmm. What if d_s comes from a DIFFERENT scaling relation?

# In the Weyl law for fractals:
# N(lambda) ~ lambda^{d_s/2}
# From our DOS fit: d_s_dos was computed above

# The SCALING of the spectrum under Galois:
# ratio = phi^4 for pair 1, phi^{1.62} for pair 2
# If these represent scaling at TWO different levels of the fractal:
# Level 1 (fine): ratio phi^4, involving 4 modes
# Level 2 (coarse): ratio phi^{1.62}, involving 9 modes

# The mean scaling exponent:
# (4 * 4 + 1.62 * 9) / (4 + 9) = (16 + 14.58) / 13 = 2.35
# Not phi either.

# Actually: let me think about this differently.
# The SPECTRAL dimension is about how the heat kernel decays,
# not about individual eigenvalue ratios.

# For the 600-cell with its discrete spectrum:
# P(t) = (1/N) * sum_k exp(-t * lambda_k)
# At intermediate t, the dominant term is the smallest nonzero eigenvalue:
# P(t) ~ (1/N) * mult_1 * exp(-t * L_1) for t >> 1/L_2
# This gives d_s -> 0 (exponential decay, not power law)

# So the 600-cell ALONE doesn't have a well-defined spectral dimension
# (it's a finite graph). The spectral dimension concept requires
# either an INFINITE graph or a LIMIT of graphs.

# The FRACTAL would arise from iterating: 600-cell -> inflate by phi ->
# new 600-cell at larger scale -> iterate...
# The resulting infinite graph would have a well-defined d_s.

print(f"""
  KEY FINDING: L_8/L_1 = phi^4 (Galois scaling = phi^4)

  The Galois map scales the first eigenvalue pair by phi^4.
  This is the SPECTRAL signature of the icosahedral scaling.

  For a fractal generated by iterating this scaling:
    If scale factor s = phi and eigenvalue ratio = phi^4 = s^4:
    Then the spectral exponent is 4.
    And d_s = 2 * d_H / d_w, where the walk exponent d_w = 2*d_H/d_s.

  But we need to be more careful. The relation
    lambda(s*R) = s^{d_w} * lambda(R)
  for a fractal gives: phi^4 = phi^{d_w}, so d_w = 4.

  Then: d_s = 2*d_H/d_w = 2*d_H/4 = d_H/2.
  For d_H = 3 (S^3): d_s = 3/2 = 1.5
  For d_H = 4 (R^4): d_s = 4/2 = 2

  Neither gives phi directly.

  HOWEVER: if d_H is itself related to phi:
  For a quasicrystal, d_H can be NON-INTEGER.
  If d_H = 2*phi (a fractal with Hausdorff dim 2*phi = 3.236):
    d_s = 2*2*phi/4 = phi. BINGO!

  So: d_s = phi iff d_H = 2*phi AND d_w = 4.

  d_w = 4 comes from the Galois scaling L_8/L_1 = phi^4.
  d_H = 2*phi would need to come from the geometry of the
  quasicrystalline inflation.
""")

# ============================================================
# Check: is d_H = 2*phi for the 600-cell quasicrystal?
# ============================================================
print(f"{'='*70}")
print("d_H = 2*phi ?")
print(f"{'='*70}")

# The Hausdorff dimension of the 600-cell quasicrystal:
# In the cut-and-project method, the "quasicrystal" is the set of
# E8 lattice points that project into the acceptance window.
# The Hausdorff dimension of this set IN THE PHYSICAL SPACE is d_H = 4
# (it fills R^4 densely).

# But if we consider the 600-cell as a SUBSET of S^3 (not all of S^3),
# its "fractal dimension" on S^3 is:
# d_H = log(N) / log(1/angular_spacing)
# angular_spacing = angular distance to nearest neighbor

angular_spacing = np.arccos(PHI/2)  # nearest neighbor dot product = phi/2
d_H_est = np.log(N) / np.log(2*np.pi / angular_spacing)

print(f"\n  Angular spacing (nearest neighbor): {angular_spacing:.6f} rad = {np.degrees(angular_spacing):.2f} deg")
print(f"  Estimated d_H = log(N)/log(2*pi/spacing) = {d_H_est:.4f}")
print(f"  2*phi = {2*PHI:.4f}")
print(f"  3 (dim S^3) = 3")

# Alternative: use the scaling of vertex count under subdivision
# If we subdivide the 600-cell by factor phi:
# Number of new vertices = N * phi^d_H (for a d_H-dim object)
# The dual (120-cell) has 600 vertices, ratio = 5 = a1
# phi^d_H = a1 => d_H = log(a1)/log(phi)
d_H_from_dual = np.log(a1) / np.log(PHI)
print(f"\n  From dual: phi^d_H = a1 = {a1}")
print(f"  d_H = log({a1})/log(phi) = {d_H_from_dual:.4f}")
print(f"  2*phi = {2*PHI:.4f}")
print(f"  Match: {abs(d_H_from_dual - 2*PHI) < 0.1}")

# WOW: log(5)/log(phi) = 3.3399... and 2*phi = 3.2360...
# Close but not exact. Difference: 0.104

# What about log(a1)/log(phi) more precisely?
# log(5)/log(phi) = ln(5)/ln(phi) = 1.6094/0.4812 = 3.3437
# 2*phi = 1+sqrt(5) = 3.2361
# Difference: 0.108

# Not a match. But interestingly close.

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
FINDINGS:
=========

1. GALOIS SCALING: L_8/L_1 = phi^4 (EXACT)
   The Galois map scales the first eigenvalue pair by phi^4.
   This implies walk dimension d_w = 4.

2. SPECTRAL DIMENSION: d_s = phi requires d_H = 2*phi = {2*PHI:.4f}
   From d_s = 2*d_H/d_w: phi = 2*d_H/4, so d_H = 2*phi.

3. HAUSDORFF DIMENSION: d_H = log(a1)/log(phi) = {d_H_from_dual:.4f}
   Close to 2*phi = {2*PHI:.4f} but NOT exact (off by {abs(d_H_from_dual - 2*PHI):.4f}).

4. 600-CELL AS FINITE GRAPH: no well-defined spectral dimension
   (exponential decay, not power law). Need infinite fractal limit.

STATUS:
  - The Galois scaling phi^4 IS exact and derived.
  - d_s = phi requires d_H = 2*phi, which is CLOSE to log(5)/log(phi) but not exact.
  - The fractal picture is SUGGESTIVE but not proven.
  - The 600-cell alone is too small for a well-defined d_s.
  - An iterated fractal construction might give d_s = phi, but this
    requires defining the inflation rule precisely.

THE HONEST ANSWER: d_s = phi is PLAUSIBLE but NOT PROVEN.
The Galois scaling phi^4 is the strongest evidence.
""")

print("=" * 70)
print("EXP-511 COMPLETE")
print("=" * 70)
