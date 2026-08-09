"""
exp562: Discrete dispersion relation and light cone on the 600-cell.

With c^2 = a1 = 5 confirmed, the wave equation Box*psi = 0 defines a
discrete light cone. The joint spectrum of (L_fiber, L_cross) gives the
full dispersion relation. Modes ON the light cone (mu=0) are massless;
modes OFF it have discrete "mass" mu = L_cross - a1*L_fiber.

GOALS:
  1. Simultaneously diagonalize L_fiber and L_cross (they commute)
  2. Map the full joint spectrum (lambda_f, lambda_c)
  3. Identify the discrete light cone lambda_c = a1*lambda_f
  4. Compute the "mass spectrum" from deviations
  5. Check if mass splittings relate to known physics
  6. Analyze the dispersion relation: is c^2 = a1 exact at all momenta?
"""

import numpy as np
from numpy.linalg import eigh
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
TOL = 1e-8

verts, adj, lap = build_600cell()

# === Build Hopf fibration ===
def qmul(p, q):
    w = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
    x = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
    y = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
    z = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]
    return np.array([w, x, y, z])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    target_w = PHI / 2.0
    for i in range(N):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6):
                    ok = False
            if not ok: continue
            used = set()
            fibers = []
            subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10):
                idx = find_idx(pp, verts)
                subg.append(idx)
                pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si])
                    idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx)
                        used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
assert fibers is not None

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5:
                A_fiber[i,j] = 1.0
A_cross = adj - A_fiber
L_fiber = np.diag(np.sum(A_fiber, axis=1)) - A_fiber
L_cross = np.diag(np.sum(A_cross, axis=1)) - A_cross
Box = L_cross - a1 * L_fiber


# =====================================================================
# STEP 1: SIMULTANEOUS DIAGONALIZATION
# =====================================================================
print("=" * 70)
print("STEP 1: SIMULTANEOUS DIAGONALIZATION OF L_fiber AND L_cross")
print("=" * 70)

# Use irrational linear combination to break all degeneracies
# M = L_fiber + pi*L_cross has (generically) non-degenerate spectrum
M = L_fiber + np.pi * L_cross
evals_M, evecs_M = eigh(M)

# Read off L_fiber and L_cross eigenvalues for each joint eigenvector
lf_evals = np.zeros(N)
lc_evals = np.zeros(N)
for i in range(N):
    v = evecs_M[:, i]
    lf_evals[i] = v @ L_fiber @ v
    lc_evals[i] = v @ L_cross @ v

# Verify these are actual eigenvalues (not just expectation values)
# Check: L_fiber @ v = lf * v ?
residuals_f = []
residuals_c = []
for i in range(N):
    v = evecs_M[:, i]
    res_f = np.linalg.norm(L_fiber @ v - lf_evals[i] * v)
    res_c = np.linalg.norm(L_cross @ v - lc_evals[i] * v)
    residuals_f.append(res_f)
    residuals_c.append(res_c)

print(f"\n  Max residual L_fiber: {max(residuals_f):.2e}")
print(f"  Max residual L_cross: {max(residuals_c):.2e}")
print(f"  Simultaneous eigenvectors: {'YES' if max(residuals_f) < 1e-8 and max(residuals_c) < 1e-8 else 'APPROXIMATE'}")

# Box eigenvalue for each mode
mu_evals = lc_evals - a1 * lf_evals

# Also get adjacency eigenvalues
A_evals = np.zeros(N)
for i in range(N):
    v = evecs_M[:, i]
    A_evals[i] = v @ adj @ v


# =====================================================================
# STEP 2: JOINT SPECTRUM TABLE
# =====================================================================
print("\n" + "=" * 70)
print("STEP 2: JOINT SPECTRUM (lf, lc, mu)")
print("=" * 70)

# Group by (lf, lc) pairs
from collections import Counter
pair_counts = Counter()
pair_data = {}
for i in range(N):
    key = (round(lf_evals[i], 4), round(lc_evals[i], 4))
    pair_counts[key] += 1
    pair_data[key] = {
        'lf': lf_evals[i], 'lc': lc_evals[i],
        'mu': mu_evals[i], 'lam_A': A_evals[i],
        'l_full': lf_evals[i] + lc_evals[i]
    }

print(f"\n  {'lf':>10s} {'lc':>10s} {'l_full':>10s} {'mu=lc-5lf':>10s} {'lam_A':>10s} {'mult':>5s} {'light cone?':>12s}")
print(f"  {'-'*70}")

on_cone = 0
for (lf_key, lc_key), mult in sorted(pair_counts.items()):
    d = pair_data[(lf_key, lc_key)]
    cone = "MASSLESS" if abs(d['mu']) < 0.01 else ""
    print(f"  {d['lf']:10.4f} {d['lc']:10.4f} {d['l_full']:10.4f} {d['mu']:10.4f} {d['lam_A']:10.4f} {mult:5d} {cone:>12s}")
    if abs(d['mu']) < 0.01:
        on_cone += mult

print(f"\n  Modes on light cone (mu=0): {on_cone} (expect 9)")


# =====================================================================
# STEP 3: FIBER AND CROSS SPECTRA (CLEAN)
# =====================================================================
print("\n" + "=" * 70)
print("STEP 3: INDIVIDUAL SPECTRA")
print("=" * 70)

# Clean fiber eigenvalues
print(f"\n  L_fiber eigenvalues (from simultaneous diag):")
lf_distinct = Counter()
for lf in lf_evals:
    lf_distinct[round(lf, 6)] += 1
for val, mult in sorted(lf_distinct.items()):
    # Try to identify exact value
    exact = ""
    for name, expr in [("0", 0), ("2-phi", 2-PHI), ("3-phi", 3-PHI),
                         ("1+phi", 1+PHI), ("2+phi", 2+PHI), ("4", 4),
                         ("1.5", 1.5), ("2", 2), ("2.5", 2.5), ("3", 3)]:
        if abs(val - expr) < 0.001:
            exact = f" = {name}"
    print(f"    lf = {val:10.6f}  mult = {mult}{exact}")

# Clean cross eigenvalues
print(f"\n  L_cross eigenvalues (from simultaneous diag):")
lc_distinct = Counter()
for lc in lc_evals:
    lc_distinct[round(lc, 4)] += 1
for val, mult in sorted(lc_distinct.items()):
    print(f"    lc = {val:10.4f}  mult = {mult}")


# =====================================================================
# STEP 4: DISPERSION RELATION ANALYSIS
# =====================================================================
print("\n" + "=" * 70)
print("STEP 4: DISPERSION RELATION")
print("=" * 70)

# For each mode, compute the "effective c^2" = lc/lf (when lf > 0)
print(f"\n  Effective c^2(mode) = lc/lf for each mode:")
print(f"  {'lf':>10s} {'lc':>10s} {'c^2_eff':>10s} {'mu':>10s} {'mult':>5s}")
print(f"  {'-'*50}")

for (lf_key, lc_key), mult in sorted(pair_counts.items()):
    d = pair_data[(lf_key, lc_key)]
    if d['lf'] > 0.01:
        c2_eff = d['lc'] / d['lf']
        print(f"  {d['lf']:10.4f} {d['lc']:10.4f} {c2_eff:10.4f} {d['mu']:10.4f} {mult:5d}")

# For kernel modes: c^2_eff = a1 exactly
# For non-kernel modes: c^2_eff varies -- this is the "mass" effect


# =====================================================================
# STEP 5: MASS SPECTRUM FROM BOX
# =====================================================================
print("\n" + "=" * 70)
print("STEP 5: MASS SPECTRUM (Box eigenvalues)")
print("=" * 70)

mu_distinct = Counter()
for mu in mu_evals:
    mu_distinct[round(mu, 4)] += 1

print(f"\n  {'mu':>10s} {'mult':>5s} {'|mu|':>10s} {'interpretation':>20s}")
print(f"  {'-'*50}")

for val, mult in sorted(mu_distinct.items()):
    interp = ""
    if abs(val) < 0.01:
        interp = "MASSLESS (light cone)"
    elif val > 0:
        interp = "spacelike"
    else:
        interp = "timelike"

    # Try to identify exact value
    exact = ""
    for name, expr in [("0", 0),
                         ("b1-phi", b1-PHI), ("-(b1-phi)", -(b1-PHI)),
                         ("a1+sqrt5", a1+np.sqrt(5)), ("-(a1+sqrt5)", -(a1+np.sqrt(5))),
                         ("a1-sqrt5", a1-np.sqrt(5)),
                         ("b1+phi", b1+PHI), ("-(b1+phi)", -(b1+PHI)),
                         ("12", 12), ("-10", -10),
                         ("phi^2", PHI**2), ("-phi^2", -PHI**2),
                         ]:
        if abs(val - expr) < 0.01:
            exact = f" = {name}"

    print(f"  {val:10.4f} {mult:5d} {abs(val):10.4f} {interp:>20s}{exact}")

# Mass ratios
nonzero_mu = sorted(set(round(abs(m), 4) for m in mu_evals if abs(m) > 0.01))
print(f"\n  Distinct |mu| values (masses): {len(nonzero_mu)}")
if len(nonzero_mu) > 0:
    mu_min = nonzero_mu[0]
    print(f"  Smallest nonzero |mu| = {mu_min:.6f}")
    print(f"\n  Mass ratios (relative to smallest):")
    for m in nonzero_mu:
        ratio = m / mu_min
        print(f"    |mu| = {m:10.4f},  ratio = {ratio:10.4f}")


# =====================================================================
# STEP 6: CONTINUUM LIMIT CHECK
# =====================================================================
print("\n" + "=" * 70)
print("STEP 6: CONTINUUM DISPERSION RELATION")
print("=" * 70)

# On S^3, the scalar Laplacian has eigenvalues l(l+2) with mult (l+1)^2
# for l = 0, 1, 2, ...
# On the 600-cell, the Laplacian matches these for l = 0, 1, ..., 5
# (beyond l=5, lattice artifacts dominate)

# The continuum dispersion relation for massless particles on S^3:
#   omega^2 = c^2 * k^2 (flat space limit)
# On S^3 with radius R:
#   omega_l^2 = l(l+2)/R^2

# On the lattice with Hopf time direction:
#   omega^2 ~ lambda_fiber (temporal eigenvalue)
#   k^2 ~ lambda_cross (spatial eigenvalue)
#   c^2 = omega^2/k^2 for massless modes

# For each "angular momentum level" l, compute the average c^2
print(f"\n  Continuum S^3 eigenvalues l(l+2) for l=0,...,5:")
for l in range(6):
    expected = l * (l + 2)
    mult = (l + 1)**2
    print(f"    l={l}: l(l+2) = {expected:3d}, mult = {mult:3d}")

# Map 600-cell Laplacian eigenvalues to angular momentum
# Known mapping (from paper):
# l=0: lambda_L = 0 (mult 1)
# l=1: lambda_L = 12-6phi (mult 4)
# l=2: lambda_L = 10-2sqrt5 (mult 9)
# l=3: lambda_L = 9 (mult 16)
# l=4: lambda_L = 12 (mult 25)
# l=5: lambda_L = 14 (mult 36)
# Beyond: Galois-conjugate / lattice artifacts

l_map = {
    0: (0, 1),
    1: (12-6*PHI, 4),
    2: (10-2*np.sqrt(5), 9),
    3: (9, 16),
    4: (12, 25),
    5: (14, 36),
}

# For each angular momentum level, compute average (lf, lc) for modes
# at that L eigenvalue
print(f"\n  Per-level analysis: average fiber/cross eigenvalues")
print(f"  {'l':>3s} {'L_full':>10s} {'<lf>':>10s} {'<lc>':>10s} {'<mu>':>10s} {'c^2_eff':>10s} {'mult':>5s}")
print(f"  {'-'*55}")

for l, (L_val, expected_mult) in l_map.items():
    # Find modes at this L eigenvalue
    indices = [i for i in range(N) if abs(lf_evals[i] + lc_evals[i] - L_val) < 0.1]
    if len(indices) > 0:
        avg_lf = np.mean([lf_evals[i] for i in indices])
        avg_lc = np.mean([lc_evals[i] for i in indices])
        avg_mu = np.mean([mu_evals[i] for i in indices])
        c2_eff = avg_lc / avg_lf if avg_lf > 0.001 else float('inf')
        print(f"  {l:3d} {L_val:10.4f} {avg_lf:10.4f} {avg_lc:10.4f} "
              f"{avg_mu:10.4f} {c2_eff:10.4f} {len(indices):5d}")

# For the Galois-conjugate levels:
galois_levels = {
    '2G': (10+2*np.sqrt(5), 9),
    '1G': (6+6*PHI, 4),
    '3+': (15, 16),
}

print(f"\n  Galois-conjugate / lattice levels:")
print(f"  {'label':>5s} {'L_full':>10s} {'<lf>':>10s} {'<lc>':>10s} {'<mu>':>10s} {'c^2_eff':>10s} {'mult':>5s}")
print(f"  {'-'*60}")

for label, (L_val, expected_mult) in galois_levels.items():
    indices = [i for i in range(N) if abs(lf_evals[i] + lc_evals[i] - L_val) < 0.1]
    if len(indices) > 0:
        avg_lf = np.mean([lf_evals[i] for i in indices])
        avg_lc = np.mean([lc_evals[i] for i in indices])
        avg_mu = np.mean([mu_evals[i] for i in indices])
        c2_eff = avg_lc / avg_lf if avg_lf > 0.001 else float('inf')
        print(f"  {label:>5s} {L_val:10.4f} {avg_lf:10.4f} {avg_lc:10.4f} "
              f"{avg_mu:10.4f} {c2_eff:10.4f} {len(indices):5d}")


# =====================================================================
# STEP 7: LIGHT CONE STRUCTURE
# =====================================================================
print("\n" + "=" * 70)
print("STEP 7: LIGHT CONE = ker(Box)")
print("=" * 70)

kernel_modes = [i for i in range(N) if abs(mu_evals[i]) < 0.01]
print(f"\n  Light cone modes: {len(kernel_modes)}")

for i, idx in enumerate(kernel_modes):
    l_full = lf_evals[idx] + lc_evals[idx]
    A_val = A_evals[idx]
    # Identify irrep
    if abs(l_full) < 0.01:
        irrep = "rho_0 (l=0, trivial)"
    elif abs(l_full - (12-6*PHI)) < 0.01:
        irrep = "rho_1 (l=1, physical)"
    elif abs(l_full - (6+6*PHI)) < 0.01:
        irrep = "rho_8 (l=1, Galois conjugate)"
    else:
        irrep = "?"

    print(f"  Mode {i}: lf={lf_evals[idx]:.6f}, lc={lc_evals[idx]:.6f}, "
          f"L={l_full:.4f}, A={A_val:.4f}, {irrep}")

print(f"""
  LIGHT CONE DECOMPOSITION:
    rho_0 (1 mode):  l=0, constant mode (vacuum)
    rho_1 (4 modes): l=1, physical sector
    rho_8 (4 modes): l=1 Galois conjugate (dark sector)
    Total: 1 + 4 + 4 = 9 = N_eig

  Physical interpretation:
    4 massless modes at l=1 in 3+1 dimensions:
    In the continuum on S^3, l=1 has mult (1+1)^2 = 4.
    These are the 4 components of a massless vector field.

  On the light cone:  lc = {a1}*lf  (exact, c^2 = a1)
  Off the light cone: lc != {a1}*lf  (massive modes, discrete mass mu)
""")


# =====================================================================
# STEP 8: SUMMARY
# =====================================================================
print("=" * 70)
print("SUMMARY")
print("=" * 70)

n_massless = sum(1 for m in mu_evals if abs(m) < 0.01)
n_timelike = sum(1 for m in mu_evals if m < -0.01)
n_spacelike = sum(1 for m in mu_evals if m > 0.01)

print(f"\n  Total modes: {N}")
print(f"  Massless (on light cone): {n_massless}")
print(f"  Timelike (mu < 0): {n_timelike}")
print(f"  Spacelike (mu > 0): {n_spacelike}")
print(f"  Light cone speed: c^2 = a1 = {a1}")
print(f"  Distinct mass values: {len(nonzero_mu)}")
print(f"  Mass range: [{nonzero_mu[0]:.4f}, {nonzero_mu[-1]:.4f}]" if nonzero_mu else "")

print("\n" + "=" * 70)
print("EXP-562 COMPLETE")
print("=" * 70)
