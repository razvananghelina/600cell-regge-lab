"""
EXP-496: Coleman-Weinberg Effective Potential on 600-Cell
==========================================================
QUESTION: Does the one-loop effective potential on the 600-cell
generate a non-trivial Higgs VEV and mass hierarchy?

SETUP:
  - Fermion "masses" = eigenvalues of Delta_0 (scalar Laplacian)
  - Background Higgs field phi couples as mass term
  - D(phi)^2 on 0-forms = phi^2 * I + Delta_0
  - One-loop effective potential (Coleman-Weinberg):
    V(phi) = (1/64pi^2) * sum_k (phi^2 + lam_k)^2 * [ln((phi^2+lam_k)/Lambda^2) - 3/2]
  - On finite graph: no UV divergence, Lambda = natural cutoff = max eigenvalue

ALSO TEST:
  - Does the VEV relate to phi (golden ratio)?
  - Does the Higgs mass relate to known values?
  - Does the hierarchy m_H/m_W emerge?

STATUS: FIRST COMPUTATION (unexplored direction)
"""

import numpy as np
from scipy.spatial.distance import cdist

phi_gr = (1 + np.sqrt(5)) / 2  # golden ratio (renamed to avoid confusion with Higgs field)
a1 = 5
b1 = 6

print("=" * 70)
print("EXP-496: COLEMAN-WEINBERG ON 600-CELL")
print("=" * 70)

# === Build 600-cell and get spectrum ===
vertices = []
for i in range(4):
    for s in [1,-1]:
        v = [0.0]*4; v[i] = s; vertices.append(v)
for s0 in [1,-1]:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                vertices.append([s0/2, s1/2, s2/2, s3/2])
p, h, q = phi_gr/2, 0.5, 1/(2*phi_gr)
for t in [[p,h,q,0],[p,q,0,h],[p,0,h,q],[h,p,0,q],[h,q,p,0],[h,0,q,p],
          [q,p,h,0],[q,h,0,p],[q,0,p,h],[0,p,q,h],[0,h,p,q],[0,q,h,p]]:
    for s0 in [1,-1]:
        for s1 in [1,-1]:
            for s2 in [1,-1]:
                for s3 in [1,-1]:
                    v = [s0*t[0], s1*t[1], s2*t[2], s3*t[3]]
                    if np.linalg.norm(v) > 0.01: vertices.append(v)
unique = []
for v in vertices:
    va = np.array(v)
    if not any(np.allclose(va, u, atol=1e-10) for u in unique):
        unique.append(va)
verts = np.array(unique)
N = len(verts)
norms = np.linalg.norm(verts, axis=1)
verts = verts / norms[:, np.newaxis]

dist_matrix = cdist(verts, verts, 'euclidean')
np.fill_diagonal(dist_matrix, np.inf)
nn_dist = np.median(dist_matrix.min(axis=1))
adj = (np.abs(dist_matrix - nn_dist) < 0.02).astype(float)
np.fill_diagonal(adj, 0)

edges = []
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j] > 0.5: edges.append((i,j))
E = len(edges)

d0 = np.zeros((E, N))
for idx, (i,j) in enumerate(edges):
    d0[idx, i] = -1; d0[idx, j] = 1

Delta_0 = d0.T @ d0
evals_0 = np.sort(np.linalg.eigvalsh(Delta_0))

print(f"Scalar Laplacian spectrum ({N} eigenvalues):")
print(f"  Range: [{evals_0[0]:.4f}, {evals_0[-1]:.4f}]")
print(f"  Distinct: {len(set(round(e, 2) for e in evals_0))}")

# Group eigenvalues
unique_evals = []
for e in evals_0:
    if not unique_evals or abs(e - unique_evals[-1][0]) > 0.01:
        unique_evals.append([e, 1])
    else:
        unique_evals[-1][1] += 1

print(f"  Eigenvalues with multiplicities:")
for lam, mult in unique_evals:
    print(f"    lam = {lam:8.4f}, mult = {mult:3d}")

# ============================================================
# COLEMAN-WEINBERG POTENTIAL
# ============================================================
print("\n" + "=" * 70)
print("COLEMAN-WEINBERG EFFECTIVE POTENTIAL")
print("=" * 70)

# V_CW(h) = (1/64pi^2) * sum_k (h^2 + lam_k)^2 * [ln((h^2+lam_k)/Lambda^2) - 3/2]
# h = Higgs field VEV
# lam_k = eigenvalues of Delta_0
# Lambda^2 = cutoff = max eigenvalue

Lambda2 = evals_0[-1]  # b1 * phi^2

# Also test with other cutoffs
cutoffs = {
    'Lambda^2 = lam_max': evals_0[-1],
    'Lambda^2 = 12 (degree)': 12.0,
    'Lambda^2 = Tr(D0)/N': np.sum(evals_0) / N,
}

for cutoff_name, Lambda2 in cutoffs.items():
    print(f"\n--- Cutoff: {cutoff_name} = {Lambda2:.4f} ---")

    # Compute V(h) for h in [0, 5]
    h_values = np.linspace(0.001, 6, 2000)
    V_values = np.zeros_like(h_values)

    for i, h in enumerate(h_values):
        h2 = h**2
        V = 0.0
        for lam in evals_0:
            x = h2 + lam
            if x > 1e-30:
                V += x**2 * (np.log(x / Lambda2) - 1.5)
        V_values[i] = V / (64 * np.pi**2)

    # Find minimum
    min_idx = np.argmin(V_values)
    h_min = h_values[min_idx]
    V_min = V_values[min_idx]

    # Also check V'(h) = 0
    dV = np.gradient(V_values, h_values)
    d2V = np.gradient(dV, h_values)

    # Find zero crossings of dV (excluding h=0)
    zero_crossings = []
    for i in range(1, len(dV)-1):
        if dV[i-1] * dV[i] < 0 and h_values[i] > 0.1:
            h_zero = h_values[i-1] - dV[i-1] * (h_values[i] - h_values[i-1]) / (dV[i] - dV[i-1])
            zero_crossings.append(h_zero)

    print(f"  V minimum at h = {h_min:.4f}")
    print(f"  V''(h_min) = {d2V[min_idx]:.4f} (Higgs mass^2)")
    if zero_crossings:
        print(f"  V'=0 crossings at h = {[f'{z:.4f}' for z in zero_crossings]}")

    # Check if h_min is related to golden ratio
    print(f"  h_min / phi = {h_min / phi_gr:.4f}")
    print(f"  h_min / sqrt(phi) = {h_min / np.sqrt(phi_gr):.4f}")
    print(f"  h_min^2 = {h_min**2:.4f}")
    print(f"  h_min^2 / phi = {h_min**2 / phi_gr:.4f}")

    # Higgs mass from curvature at minimum
    if d2V[min_idx] > 0:
        m_H_lattice = np.sqrt(d2V[min_idx])
        print(f"  m_H (lattice units) = {m_H_lattice:.4f}")
        print(f"  m_H / h_min = {m_H_lattice / max(h_min, 1e-10):.4f}")

# ============================================================
# INTERPRETATION
# ============================================================
print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)

# The CW potential with 120 eigenvalues and cutoff Lambda^2 = lam_max
Lambda2 = evals_0[-1]
h_fine = np.linspace(0.001, 6, 10000)
V_fine = np.zeros_like(h_fine)

for i, h in enumerate(h_fine):
    h2 = h**2
    V = 0.0
    for lam in evals_0:
        x = h2 + lam
        if x > 1e-30:
            V += x**2 * (np.log(x / Lambda2) - 1.5)
    V_fine[i] = V / (64 * np.pi**2)

dV = np.gradient(V_fine, h_fine)
d2V = np.gradient(dV, h_fine)

min_idx = np.argmin(V_fine)
h_vev = h_fine[min_idx]

print(f"\nWith natural cutoff Lambda^2 = lam_max = b1*phi^2 = {Lambda2:.4f}:")
print(f"  Higgs VEV: h = {h_vev:.6f}")
print(f"  V(h_vev) = {V_fine[min_idx]:.6f}")
print(f"  V''(h_vev) = {d2V[min_idx]:.6f}")

# Ratios to theory constants
print(f"\n  Ratios to known constants:")
print(f"  h_vev = {h_vev:.4f}")
print(f"  phi = {phi_gr:.4f}")
print(f"  h_vev / phi = {h_vev / phi_gr:.4f}")
print(f"  h_vev^2 = {h_vev**2:.4f}")
print(f"  h_vev^2 - phi = {h_vev**2 - phi_gr:.4f}")
print(f"  h_vev / sqrt(lam_max) = {h_vev / np.sqrt(Lambda2):.4f}")

# Nearby eigenvalues: is h_vev near any special value?
for name, val in [('sqrt(2)', np.sqrt(2)), ('phi', phi_gr),
                   ('sqrt(phi)', np.sqrt(phi_gr)), ('1/phi', 1/phi_gr),
                   ('sqrt(a1)', np.sqrt(a1)), ('sqrt(b1)', np.sqrt(b1)),
                   ('sqrt(12-6phi)', np.sqrt(12-6*phi_gr)),
                   ('sqrt(7-4phi)', np.sqrt(7-4*phi_gr)),
                   ('phi^2-1', phi_gr**2 - 1)]:
    print(f"  |h_vev - {name}| = {abs(h_vev - val):.4f}")

# ============================================================
# ALTERNATIVE: Use ADJACENCY eigenvalues (not Laplacian)
# ============================================================
print("\n" + "=" * 70)
print("ALTERNATIVE: CW WITH ADJACENCY EIGENVALUES")
print("=" * 70)

evals_adj = 12 - evals_0  # adjacency eigenvalues

# In this version: "masses" are proportional to h * y_k where y_k = adjacency eigenvalue
# V(h) = sum_k (h * y_k)^4 * [ln((h*y_k)^2 / Lambda^2) - 3/2]

Lambda2_adj = (12)**2  # max eigenvalue squared

h_values = np.linspace(0.01, 5, 5000)
V_adj = np.zeros_like(h_values)

for i, h in enumerate(h_values):
    V = 0.0
    for y in evals_adj:
        m2 = (h * y)**2
        if m2 > 1e-30:
            V += m2**2 * (np.log(m2 / Lambda2_adj) - 1.5)
    V_adj[i] = V / (64 * np.pi**2)

min_idx_adj = np.argmin(V_adj)
h_vev_adj = h_values[min_idx_adj]

dV_adj = np.gradient(V_adj, h_values)
d2V_adj = np.gradient(dV_adj, h_values)

print(f"  Higgs VEV (adjacency): h = {h_vev_adj:.4f}")
print(f"  h_vev / phi = {h_vev_adj / phi_gr:.4f}")
print(f"  h_vev^2 = {h_vev_adj**2:.4f}")

# ============================================================
# THE REAL TEST: Does CW generate mass SPLITTING?
# ============================================================
print("\n" + "=" * 70)
print("KEY TEST: MASS SPLITTING FROM CW")
print("=" * 70)

print("""
The CW mechanism generates masses through dimensional transmutation.
On the 600-cell, the "bare" masses are the Laplacian eigenvalues
(9 distinct values with multiplicities 1,4,9,16,25,36,9,16,4).

Question: does the one-loop correction SPLIT these degenerate
eigenvalues into non-degenerate ones that match the fermion masses?
""")

# At the VEV h = h_vev, the effective masses are:
# m_eff_k = sqrt(h_vev^2 + lam_k)
# These are NOT split within each eigenspace (all modes with same lam_k
# get the same effective mass).

# For splitting, we need an ANISOTROPIC Higgs field (not uniform).
# A uniform Higgs shifts all eigenvalues equally -> no splitting.

print(f"At VEV h = {h_vev:.4f}:")
print(f"  Effective masses m_k = sqrt(h^2 + lam_k):")
for lam, mult in unique_evals:
    m_eff = np.sqrt(h_vev**2 + lam)
    print(f"    lam = {lam:8.4f} (mult {mult:3d}): m_eff = {m_eff:8.4f}")

print(f"\n  A UNIFORM Higgs field cannot split the degeneracies.")
print(f"  For mass splitting, need ANISOTROPIC Higgs (breaks 2I symmetry).")

# ============================================================
# ANISOTROPIC HIGGS: phi-eigenvector direction
# ============================================================
print("\n" + "=" * 70)
print("ANISOTROPIC HIGGS: phi-EIGENVECTOR DIRECTION")
print("=" * 70)

# Use the eigenvector of A with eigenvalue 6*phi_gr as the Higgs direction
evals_adj_full, evecs_adj = np.linalg.eigh(adj.astype(float))
idx_sort = np.argsort(evals_adj_full)[::-1]
evecs_adj = evecs_adj[:, idx_sort]

# The phi-eigenspace (eigenvalue 6*phi, multiplicity 4)
# These are modes 1-4 (after sorting by descending eigenvalue)
phi_modes = evecs_adj[:, 1:5]  # 4 modes with eigenvalue 6*phi

# Use the FIRST phi-mode as Higgs direction
higgs_dir = phi_modes[:, 0]  # (120,) vector

# Higgs field: H(v) = h * higgs_dir(v)
# Mass matrix: M_ij = h * higgs_dir(i) * delta_ij (diagonal mass on vertices)
# Modified Laplacian: Delta_0(h) = Delta_0 + h^2 * diag(higgs_dir^2)

print(f"Higgs direction: phi-eigenvector of adjacency matrix")
print(f"  Eigenvalue: 6*phi = {6*phi_gr:.4f}")
print(f"  ||higgs_dir|| = {np.linalg.norm(higgs_dir):.4f}")

h_values = np.linspace(0.0, 20, 500)
spectra = []

for h in h_values:
    Delta_h = Delta_0 + h**2 * np.diag(higgs_dir**2)
    evals_h = np.sort(np.linalg.eigvalsh(Delta_h))
    spectra.append(evals_h)

spectra = np.array(spectra)

# At h=0: 9 distinct eigenvalues (degenerate)
# At h>0: how many distinct eigenvalues?
print(f"\nNumber of distinct eigenvalues vs Higgs VEV h:")
for h_idx in [0, 50, 100, 200, 300, 499]:
    h = h_values[h_idx]
    ev = spectra[h_idx]
    n_distinct = len(set(round(e, 4) for e in ev))
    spread = ev[-1] - ev[0]
    print(f"  h = {h:6.2f}: {n_distinct:3d} distinct eigenvalues, "
          f"spread = [{ev[0]:.2f}, {ev[-1]:.2f}]")

# Check: does the phi-Higgs direction split the degeneracies?
print(f"\nDegeneracy splitting at h = {h_values[100]:.2f}:")
ev_h = spectra[100]
groups_h = []
for e in ev_h:
    if not groups_h or abs(e - groups_h[-1][0]) > 0.001:
        groups_h.append([e, 1])
    else:
        groups_h[-1][1] += 1

for lam_h, mult_h in groups_h:
    print(f"  lam = {lam_h:8.4f}, mult = {mult_h}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
COLEMAN-WEINBERG ON 600-CELL:
==============================

1. UNIFORM HIGGS:
   V_CW(h) has a minimum at h ~ {h_vev:.2f} (depends on cutoff).
   But a uniform Higgs cannot split eigenvalue degeneracies.
   All modes within the same irrep get the same effective mass.
   -> NO mass hierarchy from uniform CW.

2. ANISOTROPIC HIGGS (phi-eigenvector direction):
   Breaking 2I symmetry along the phi-eigenvector direction
   DOES split the degeneracies. The 9 degenerate eigenvalues
   become {len(groups_h)} distinct values at finite h.

3. THE FUNDAMENTAL ISSUE:
   The CW mechanism generates a VEV and splits degeneracies,
   but the PATTERN of splitting depends on the CHOICE of Higgs
   direction. Different directions give different splittings.
   There's no unique prediction without knowing the Higgs direction.

4. WHAT WOULD MAKE THIS WORK:
   If the 600-cell geometry SELECTS a unique Higgs direction
   (e.g., through the Galois kernel, or through minimization of
   the full anisotropic CW potential), then the splitting pattern
   would be a PREDICTION.

STATUS: The CW mechanism EXISTS on the 600-cell and CAN split
degeneracies, but without a principle to select the Higgs direction,
it doesn't predict the mass hierarchy.
""")

print("=" * 70)
print("EXP-496 COMPLETE")
print("=" * 70)
