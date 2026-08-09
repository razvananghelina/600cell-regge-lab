"""
exp563: Are W and Z bosons the first massive modes of Box?

The Box operator has 9 massless modes (light cone) and 111 massive modes.
The first two nonzero masses are:
  |mu1| = 6/phi - 2 = 4*phi - 4 ~ 1.708
  |mu2| = 5 - sqrt(5)            ~ 2.764
  ratio mu2/mu1 = phi (exact!)

QUESTION: Do these relate to m_W or m_Z through framework-derived factors?

The electroweak masses in the framework:
  m_W = m_e * phi^25 * R_W   (from Section 8)
  m_Z = m_W / cos(theta_W)
  m_Z/m_W = 1/cos(theta_W) = sqrt(26/20) = sqrt(13/10)

But Box eigenvalues are in LATTICE UNITS. To compare with physical masses,
we need the lattice-to-physical conversion, which involves m_e and phi^25.

The more useful comparison is DIMENSIONLESS RATIOS:
  mu2/mu1 = phi (we proved this)
  m_Z/m_W = sqrt(13/10) = 1.1402 (from sin^2(tW) = 6/26)

Also: m_W and m_Z are massive vector bosons (spin-1).
In the Box spectrum, which modes carry the right quantum numbers?

APPROACH:
  1. List all Box mass values with their A5 quantum numbers
  2. Identify which modes transform as gauge bosons
  3. Compare mass ratios with m_Z/m_W, m_H/m_W, etc.
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

verts, adj, lap = build_600cell()

# Build Hopf fibration
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
                subg.append(idx); pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si])
                    idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0
A_cross = adj - A_fiber
L_fiber = np.diag(np.sum(A_fiber, axis=1)) - A_fiber
L_cross = np.diag(np.sum(A_cross, axis=1)) - A_cross
Box = L_cross - a1 * L_fiber

# Simultaneous diagonalization
M = L_fiber + np.pi * L_cross
evals_M, evecs_M = eigh(M)
lf_all = np.array([evecs_M[:,i] @ L_fiber @ evecs_M[:,i] for i in range(N)])
lc_all = np.array([evecs_M[:,i] @ L_cross @ evecs_M[:,i] for i in range(N)])
mu_all = lc_all - a1 * lf_all
A_all = np.array([evecs_M[:,i] @ adj @ evecs_M[:,i] for i in range(N)])
L_all = lf_all + lc_all

# =====================================================================
# PART 1: Full mass spectrum with quantum numbers
# =====================================================================
print("=" * 70)
print("PART 1: MASS SPECTRUM WITH QUANTUM NUMBERS")
print("=" * 70)

# 600-cell Laplacian eigenvalues -> irrep labels
irrep_map = {
    0.0: ('rho_0', 1, 1, 'l=0'),
    12-6*PHI: ('rho_1', 2, 4, 'l=1'),
    10-2*np.sqrt(5): ('rho_2', 3, 9, 'l=2'),
    9.0: ('rho_3', 4, 16, 'l=3'),
    12.0: ('rho_4', 5, 25, 'l=4'),
    14.0: ('rho_5', 6, 36, 'l=5'),
    10+2*np.sqrt(5): ('rho_6', 3, 9, 'l=2G'),
    15.0: ('rho_7', 4, 16, 'l=3+'),
    6+6*PHI: ('rho_8', 2, 4, 'l=1G'),
}

print(f"\n  {'mu':>10s} {'|mu|':>8s} {'L_full':>8s} {'irrep':>8s} {'dim_irr':>8s} {'l':>5s} {'lf':>8s} {'lc':>8s} {'mult':>5s}")
print(f"  {'-'*75}")

from collections import Counter
seen = Counter()
entries = []
for i in range(N):
    key = (round(mu_all[i], 3), round(L_all[i], 3))
    seen[key] += 1

for key in sorted(seen.keys()):
    mu_val, L_val = key
    mult = seen[key]
    # Find irrep
    irrep_name = '?'
    dim_irr = '?'
    l_label = '?'
    for L_ref, (name, dim, mult_expected, ll) in irrep_map.items():
        if abs(L_val - L_ref) < 0.1:
            irrep_name = name
            dim_irr = dim
            l_label = ll
            break
    # Find (lf, lc) for this mode
    idx = None
    for i in range(N):
        if abs(mu_all[i] - mu_val) < 0.01 and abs(L_all[i] - L_val) < 0.01:
            idx = i; break
    lf = lf_all[idx] if idx is not None else 0
    lc = lc_all[idx] if idx is not None else 0

    cone = "*" if abs(mu_val) < 0.01 else " "
    print(f"  {mu_val:10.4f} {abs(mu_val):8.4f} {L_val:8.4f} {str(irrep_name):>8s} {str(dim_irr):>8s} {str(l_label):>5s} "
          f"{lf:8.4f} {lc:8.4f} {mult:5d} {cone}")
    entries.append({'mu': mu_val, 'L': L_val, 'irrep': irrep_name, 'mult': mult,
                    'lf': lf, 'lc': lc, 'dim_irr': dim_irr})


# =====================================================================
# PART 2: ELECTROWEAK MASS RATIOS
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: ELECTROWEAK MASS RATIOS")
print("=" * 70)

sin2_tW = b1 / (a1**2 + 1)  # = 6/26
cos2_tW = 1 - sin2_tW       # = 20/26
mZ_over_mW = 1 / np.sqrt(cos2_tW)  # = sqrt(26/20) = sqrt(13/10)

print(f"\n  Framework electroweak predictions:")
print(f"    sin^2(theta_W) = {sin2_tW:.6f} = 6/26")
print(f"    cos^2(theta_W) = {cos2_tW:.6f} = 20/26")
print(f"    m_Z / m_W = 1/cos(theta_W) = {mZ_over_mW:.6f} = sqrt(13/10)")

# Box mass values (absolute, sorted, distinct)
mu_abs_distinct = sorted(set(round(abs(mu_all[i]), 6) for i in range(N)))
mu_nonzero = [m for m in mu_abs_distinct if m > 0.01]

print(f"\n  Box mass spectrum (|mu|, distinct):")
for i, m in enumerate(mu_nonzero):
    print(f"    mu_{i+1} = {m:.6f}")

# Compare ratios
print(f"\n  --- RATIO COMPARISON ---")
print(f"  m_Z/m_W = {mZ_over_mW:.6f}")
print(f"  mu2/mu1 = {mu_nonzero[1]/mu_nonzero[0]:.6f} = phi = {PHI:.6f}")
print(f"  These are DIFFERENT: phi != sqrt(13/10)")
print(f"  phi = {PHI:.6f}, sqrt(13/10) = {mZ_over_mW:.6f}")

# Check if any pair of masses gives m_Z/m_W
print(f"\n  Searching for mass pair with ratio = sqrt(13/10) = {mZ_over_mW:.6f}:")
found_ew = False
for i, m1 in enumerate(mu_nonzero):
    for j, m2 in enumerate(mu_nonzero):
        if j > i:
            ratio = m2 / m1
            if abs(ratio - mZ_over_mW) < 0.01:
                print(f"    mu_{j+1}/mu_{i+1} = {m2:.4f}/{m1:.4f} = {ratio:.6f} ~ sqrt(13/10)!")
                found_ew = True

if not found_ew:
    print(f"    No exact match found.")

# Check individual masses against EW scale factors
# In the framework: m_Z = m_e * phi^n0 * ... where n0 = a1^2 = 25
# So the EW scale is m_e * phi^25 ~ 91 GeV
# The lattice mass mu is in units of the lattice spacing
# The conversion factor is: m_phys = mu * m_e * phi^(some power) / (lattice factor)

# Let's check dimensionless combinations:
print(f"\n  --- DIMENSIONLESS COMBINATIONS ---")

# m_H/m_W in framework
# m_H^2/m_W^2 = phi^2 - 16*alpha*phi (from paper)
alpha_phys = (4*a1*PHI**4 - np.sqrt((4*a1*PHI**4)**2 - 8*np.pi)) / (4*np.pi)
mH2_over_mW2 = PHI**2 - 16*alpha_phys*PHI
mH_over_mW = np.sqrt(mH2_over_mW2)
print(f"  m_H/m_W = {mH_over_mW:.6f} (framework)")
print(f"  m_Z/m_W = {mZ_over_mW:.6f} (framework)")
print(f"  m_H/m_Z = {mH_over_mW/mZ_over_mW:.6f} (framework)")

# Check: mu_i / mu_j for all pairs vs mH/mW, mZ/mW, mH/mZ
targets = {
    'm_Z/m_W': mZ_over_mW,
    'm_H/m_W': mH_over_mW,
    'm_H/m_Z': mH_over_mW / mZ_over_mW,
    'phi': PHI,
    '1/phi': 1/PHI,
    'phi^2': PHI**2,
    'sqrt(5)': np.sqrt(5),
    'a1-1': a1-1,
}

print(f"\n  Mass ratio scan (all pairs vs EW ratios):")
for label, target in targets.items():
    best_match = None
    best_err = 999
    for i, m1 in enumerate(mu_nonzero):
        for j, m2 in enumerate(mu_nonzero):
            if j != i:
                ratio = m2 / m1
                err = abs(ratio - target)
                if err < best_err:
                    best_err = err
                    best_match = (i+1, j+1, m1, m2, ratio)
    if best_match and best_err < 0.05:
        i, j, m1, m2, ratio = best_match
        print(f"    {label:>12s} = {target:.4f}: best mu_{j}/mu_{i} = {m2:.4f}/{m1:.4f} = {ratio:.4f} (err {best_err:.4f})")
    else:
        print(f"    {label:>12s} = {target:.4f}: no close match (best err {best_err:.4f})")


# =====================================================================
# PART 3: BOX MODES AS GAUGE BOSONS?
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: GAUGE BOSON INTERPRETATION")
print("=" * 70)

print(f"""
  The light cone (ker Box) has 9 modes: rho_0(1) + rho_1(4) + rho_8(4).

  In the SM: 12 gauge bosons = 8 gluons + W+ + W- + Z + photon.
  Before EWSB: all 12 are massless.
  After EWSB:  9 massless (8 gluons + photon) + 3 massive (W+, W-, Z).

  The Box kernel has dim = 9. This matches the POST-EWSB count!
  If we identify:
    rho_0 (dim 1) = U(1) photon
    rho_1 (dim 4) = 4 of the 8 SU(3) gluons (or spacetime vectors)
    rho_8 (dim 4) = 4 Galois-conjugate modes

  The 3 missing massless modes (W+, W-, Z) are the first massive modes.
  In the Box spectrum, the first massive modes have |mu1| = {mu_nonzero[0]:.4f} (mult 12).

  But we need mult = 3, not 12. And the decomposition should match
  the SU(2) adjoint (dim 3).
""")

# Which modes have the smallest nonzero |mu|?
for e in entries:
    if abs(abs(e['mu']) - mu_nonzero[0]) < 0.01:
        print(f"  First massive mode: mu={e['mu']:.4f}, L={e['L']:.4f}, "
              f"irrep={e['irrep']}, mult={e['mult']}")

print(f"\n  The first massive modes are at |mu| = {mu_nonzero[0]:.4f}")
print(f"  with multiplicity 12 in irrep rho_5 (l=5, dim_irr=6)")
print(f"  This does NOT match W/Z (need mult 3 in SU(2) adjoint).")

# What about the modes with lf closest to the kernel?
# Kernel has lf = 0 (rho_0), lf = 2-phi (rho_1), lf = 1+phi (rho_8)
# The "nearest" massive modes to the kernel on the light cone would be
# modes with lf close to these values but different lc
print(f"\n  Modes adjacent to the light cone:")
for e in entries:
    if abs(e['mu']) > 0.01 and abs(e['mu']) < 3.0:
        print(f"    mu={e['mu']:+8.4f}, L={e['L']:.4f}, irrep={e['irrep']}, "
              f"lf={e['lf']:.4f}, lc={e['lc']:.4f}, mult={e['mult']}")


# =====================================================================
# PART 4: SYMMETRY BREAKING PATTERN
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: MASS SPLITTINGS AND GOLDEN RATIO")
print("=" * 70)

# The mass spectrum has a golden ratio structure.
# Let's express all masses in terms of phi and integers.
print(f"\n  Algebraic mass spectrum (all in Z[phi]):")
mu_exact = [
    (0, "0", 9, "ker Box (massless)"),
    (6/PHI - 2, "6/phi - 2 = 4(phi-1)", 12, "rho_5, l=5"),
    (5 - np.sqrt(5), "5 - sqrt(5) = a1 - sqrt(a1)", 6, "rho_2, l=2"),
    (6/PHI, "6/phi = 6(phi-1)", 10, "rho_4, l=4"),
    (10 - 2*np.sqrt(5), "10 - 2*sqrt(5) = L_2", 3, "rho_2, l=2 (lf=0)"),
    (6*PHI - 3, "6*phi - 3", 16, "rho_3/rho_7"),
    (5 + np.sqrt(5), "5 + sqrt(5) = a1 + sqrt(a1)", 6, "rho_6, l=2G"),
    (6*PHI, "6*phi", 10, "rho_3, l=3 (timelike)"),
    (10, "10 = 2*a1", 12, "rho_5, l=5 (max timelike)"),
    (6*PHI + 2, "6*phi + 2", 12, "rho_5, l=5 (spacelike)"),
    (12, "12 = degree", 5, "rho_4, l=4 (lf=0)"),
    (10 + 2*np.sqrt(5), "10 + 2*sqrt(5) = L_6", 3, "rho_6, l=2G (lf=0)"),
]

for mu_val, expr, mult, desc in mu_exact:
    print(f"    |mu| = {mu_val:8.4f} = {expr:25s}  mult={mult:3d}  ({desc})")

# Golden ratio chain
print(f"\n  Mass ratios involving phi:")
m1 = 6/PHI - 2
m2 = 5 - np.sqrt(5)
m3 = 6/PHI
m5 = 5 + np.sqrt(5)
print(f"    m2/m1 = phi = {m2/m1:.6f} (exact)")
print(f"    m3/m1 = {m3/m1:.6f} = 6/(6-2*phi) = 3/(3-phi)")
print(f"    m5/m1 = {m5/m1:.6f} = (5+sqrt5)/(6/phi-2)")

# Rationalize m5/m1:
# (5+sqrt5)/(6/phi-2) = (5+sqrt5)/(-5+3*sqrt5)
# multiply by (-5-3sqrt5)/(-5-3sqrt5):
# num = (5+sqrt5)(-5-3sqrt5) = -25-15sqrt5-5sqrt5-3*5 = -40-20sqrt5 = -20(2+sqrt5)
# den = 25-45 = -20
# ratio = 2+sqrt5 = phi^3 (since phi^3 = 2phi+1 = 2+sqrt5... wait: phi^3 = phi*phi^2 = phi*(phi+1) = phi^2+phi = (phi+1)+phi = 2phi+1)
# 2phi+1 = 2*1.618+1 = 4.236 = 2+sqrt5 = 2+2.236 = 4.236. YES!
print(f"    m5/m1 = 2+sqrt(5) = phi^3 = {PHI**3:.6f} (exact)")
print(f"    (Proof: (5+sqrt5)/(6/phi-2) = 20(2+sqrt5)/20 = 2+sqrt5 = phi^3)")


# =====================================================================
# PART 5: SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
  1. The Box mass spectrum is FULLY ALGEBRAIC in Z[phi].
     12 distinct mass values, all determined by a1=5.

  2. Mass ratios involve phi:
     m2/m1 = phi (exact)
     m5/m1 = phi^3 (exact)

  3. The EW ratio m_Z/m_W = sqrt(13/10) does NOT appear as a Box mass ratio.
     The Box masses are LATTICE masses, not physical particle masses.

  4. ker(Box) = 9 matches the post-EWSB massless boson count (8 gluons + photon).
     But the first massive modes (mult 12) do NOT have the right quantum numbers
     for W/Z (which need mult 3 in SU(2) adjoint).

  5. The Box operator acts on SCALAR modes on the 600-cell vertices.
     W and Z are VECTOR bosons -- they would appear in the coexact spectrum,
     not in the scalar spectrum. The Box spectrum gives the scalar sector.

  VERDICT: The Box mass spectrum is beautiful (golden ratio structure),
  but it describes SCALAR mode masses, not gauge boson masses.
  W/Z should come from the VECTOR (coexact) sector, not from Box.
  The connection, if any, is indirect.
""")
