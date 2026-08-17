"""
exp574: One-loop effective potential from the full Dirac spectrum.

D^2 = Delta_0 + Delta_1 + Delta_2 + Delta_3 (Hodge Laplacians).
Only Delta_0 and Delta_1 are affected by the U(1) gauge twist:
  Delta_0(theta) = d0(theta)^* d0(theta)   [120 x 120]
  Delta_1(theta) = d0(theta) d0(theta)^* + d1^* d1  [720 x 720]

The one-loop effective potential:
  Gamma(theta) = -(N_gen/2) * [log det Delta_0(theta)/Delta_0(0)
                               + log det Delta_1(theta)/Delta_1(0)]

Combined with tree-level YM:
  S_eff(theta) = K_YM * theta^2 / g^2 + Gamma(theta)

Saddle point dS_eff/dtheta = 0 determines the self-consistent coupling.
"""

import numpy as np
from numpy.linalg import eigvalsh, norm
from collections import defaultdict
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
N_gen = 3

verts, adj, lap = build_600cell()

# === Hopf fibration ===
def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    for i in range(N):
        if abs(verts[i, 0] - PHI/2) < 1e-6:
            g = verts[i]; p = g.copy(); ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12:
                cycles = []
                for fib in fibers:
                    fib_set = set(fib)
                    cycle = [fib[0]]; visited = {fib[0]}
                    while len(cycle) < 10:
                        curr = cycle[-1]
                        for j in fib_set - visited:
                            if adj[curr, j] > 0.5: cycle.append(j); visited.add(j); break
                    cycles.append(cycle)
                return fibers, cycles
    return None, None

fibers, cycles = find_fibration()
fiber_forward = set()
for cyc in cycles:
    for k in range(10):
        fiber_forward.add((cyc[k], cyc[(k+1)%10]))

# === Simplicial complex ===
print("Building simplicial complex...")
adj_list = defaultdict(set)
edges = []; edge_to_idx = {}
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j] > 0.5:
            adj_list[i].add(j); adj_list[j].add(i)
            edge_to_idx[(i,j)] = len(edges); edges.append((i,j))
Ne = len(edges)

triangles = []
for i in range(N):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j: triangles.append((i,j,k))
Nf = len(triangles)

# d1 (untwisted)
d1 = np.zeros((Nf, Ne))
for f_idx, (i,j,k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i,j)]] = +1.0
    d1[f_idx, edge_to_idx[(j,k)]] = +1.0
    d1[f_idx, edge_to_idx[(i,k)]] = -1.0

# Coexact Laplacian (theta-independent)
C = d1.T @ d1  # 720 x 720

print(f"  V={N}, E={Ne}, F={Nf}")


def build_d0(theta):
    """Gauge-twisted d0: vertices -> edges."""
    d0 = np.zeros((Ne, N), dtype=complex)
    for e_idx, (i, j) in enumerate(edges):
        if (i,j) in fiber_forward:
            d0[e_idx, i] = -np.exp(1j * theta / 10)
        elif (j,i) in fiber_forward:
            d0[e_idx, i] = -np.exp(-1j * theta / 10)
        else:
            d0[e_idx, i] = -1.0
        d0[e_idx, j] = 1.0
    return d0


def compute_hodge_logdets(theta):
    """Compute log det' of Delta_0(theta) and Delta_1(theta)."""
    d0_th = build_d0(theta)
    d0_dag = d0_th.conj().T

    # Delta_0 = d0^dag d0  (120 x 120)
    Delta_0 = d0_dag @ d0_th
    evals_0 = np.sort(np.real(eigvalsh(Delta_0)))
    nonzero_0 = evals_0[evals_0 > 1e-8]
    logdet_0 = np.sum(np.log(nonzero_0))

    # Delta_1 = d0 d0^dag + C  (720 x 720)
    B_th = d0_th @ d0_dag  # exact part (720 x 720)
    Delta_1 = B_th + C
    evals_1 = np.sort(np.real(eigvalsh(Delta_1)))
    nonzero_1 = evals_1[evals_1 > 1e-8]
    logdet_1 = np.sum(np.log(nonzero_1))

    return logdet_0, logdet_1, len(nonzero_0), len(nonzero_1)


# =====================================================================
# STEP 1: Reference values at theta = 0
# =====================================================================
print("\n" + "=" * 70)
print("STEP 1: REFERENCE (theta = 0)")
print("=" * 70)

logdet0_ref, logdet1_ref, n0, n1 = compute_hodge_logdets(0.0)
print(f"  log det'(Delta_0(0)) = {logdet0_ref:.6f}  ({n0} nonzero eigenvalues)")
print(f"  log det'(Delta_1(0)) = {logdet1_ref:.6f}  ({n1} nonzero eigenvalues)")
print(f"  Total D^2 contribution = {logdet0_ref + logdet1_ref:.6f}")


# =====================================================================
# STEP 2: One-loop potential scan
# =====================================================================
print("\n" + "=" * 70)
print("STEP 2: ONE-LOOP POTENTIAL SCAN")
print("=" * 70)

n_scan = 50  # Reduced for speed (720x720 eigenvalues are expensive)
thetas = np.linspace(0.001, 0.3, n_scan)

Gamma_0 = np.zeros(n_scan)  # from Delta_0
Gamma_1 = np.zeros(n_scan)  # from Delta_1
Gamma_total = np.zeros(n_scan)

print(f"  Computing {n_scan} flux values (this takes a moment)...")
for idx, theta in enumerate(thetas):
    ld0, ld1, _, _ = compute_hodge_logdets(theta)
    Gamma_0[idx] = ld0 - logdet0_ref
    Gamma_1[idx] = ld1 - logdet1_ref
    Gamma_total[idx] = Gamma_0[idx] + Gamma_1[idx]
    if idx % 10 == 0:
        print(f"    theta = {theta:.4f}: Gamma_0 = {Gamma_0[idx]:+.6f}, "
              f"Gamma_1 = {Gamma_1[idx]:+.6f}, total = {Gamma_total[idx]:+.6f}")

print(f"\n  Summary:")
print(f"  {'theta':>8s} {'Gamma_0':>12s} {'Gamma_1':>12s} {'Gamma_total':>12s}")
print(f"  {'-'*48}")
for idx in [0, n_scan//4, n_scan//2, 3*n_scan//4, n_scan-1]:
    print(f"  {thetas[idx]:8.4f} {Gamma_0[idx]:+12.6f} {Gamma_1[idx]:+12.6f} "
          f"{Gamma_total[idx]:+12.6f}")


# =====================================================================
# STEP 3: Small-theta expansion
# =====================================================================
print("\n" + "=" * 70)
print("STEP 3: SMALL-THETA EXPANSION")
print("=" * 70)

# Fit Gamma(theta) = a * theta^2 + b * theta^4 at small theta
mask_small = thetas < 0.1
x = thetas[mask_small]**2
y0 = Gamma_0[mask_small]
y1 = Gamma_1[mask_small]
yt = Gamma_total[mask_small]

fit0 = np.polyfit(x, y0, 1)
fit1 = np.polyfit(x, y1, 1)
fitt = np.polyfit(x, yt, 1)

a0, b0 = fit0[1], fit0[0]  # Gamma_0 ~ a0*theta^2
a1_coeff, b1_coeff = fit1[1], fit1[0]  # Gamma_1 ~ a1*theta^2
at, bt = fitt[1], fitt[0]

print(f"  Gamma_0(theta) ~ {a0:.6f} * theta^2  (scalar sector)")
print(f"  Gamma_1(theta) ~ {a1_coeff:.6f} * theta^2  (vector sector)")
print(f"  Gamma_total(theta) ~ {at:.6f} * theta^2  (both sectors)")
print(f"  Gamma_total = Gamma_0 + Gamma_1: {abs(at - a0 - a1_coeff) < 0.01}")

# The fermionic one-loop correction for N_gen generations:
# delta(1/g^2) = -N_gen * at / 2
delta_coupling = -N_gen * at / 2
print(f"\n  Fermionic one-loop correction to 1/g^2:")
print(f"    delta(1/g^2) = -N_gen * Gamma_coeff / 2 = -{N_gen} * {at:.4f} / 2 = {delta_coupling:.6f}")

# Tree-level: S_YM = K_YM * theta^2 / g^2, so effective:
# K_YM / g^2_eff = K_YM / g^2_tree + delta_coupling
# For K_YM = 6:
K_YM = b1  # = 6
print(f"    K_YM = {K_YM}")
print(f"    Effective: K_YM/g^2_eff = K_YM/g^2_tree + delta_coupling")

# What is the effective coupling correction?
# 1/alpha_eff = 4*pi * K_YM/g^2_eff
# The correction to 1/alpha from the one-loop:
# delta(1/alpha) = 4*pi * delta(K_YM/g^2) = ... this needs more care

# Actually, the effective action at small theta is:
# S_eff = (K_YM/g^2 - N_gen*at/2) * theta^2
# = (K_YM/g^2 + delta) * theta^2
# where delta = -N_gen*at/2

# For self-consistency at theta = 2*pi*alpha:
# The total coefficient of theta^2 should be 1/(4*pi*alpha)...
# Actually this depends on the normalization conventions.

# Let me instead just check: what is at (the Dirac one-loop coefficient)?
print(f"\n  The key coefficient: Gamma_total/theta^2 = {at:.6f}")
print(f"  Compare with known quantities:")
print(f"    K_YM = {K_YM}")
print(f"    2*pi = {2*np.pi:.4f}")
print(f"    at * N_gen = {at * N_gen:.6f}")
print(f"    at * N_gen / 2 = {at * N_gen / 2:.6f}")

# The holonomy correction in the alpha equation is 2*pi*alpha.
# For the self-consistency:
# S_eff = (K_YM/g^2) * theta^2 - (N_gen/2) * at * theta^2
# = (K_YM/g^2 - N_gen*at/2) * theta^2

# Setting theta = 2*pi*alpha and g^2 = 4*pi*alpha:
# S_eff = (K_YM/(4*pi*alpha) - N_gen*at/2) * (2*pi*alpha)^2
# = (K_YM/(4*pi*alpha) - N_gen*at/2) * 4*pi^2*alpha^2
# = K_YM*pi*alpha - 2*N_gen*at*pi^2*alpha^2

# The alpha equation: 2*pi*alpha^2 - B*alpha + 1 = 0
# Needs: the coefficient of alpha^2 = 2*pi and the coefficient of alpha = B
# From the effective action: ???

# Let me try a DIFFERENT approach.
# The self-consistency equation is:
#   1/alpha = 1/alpha_0 - delta(alpha)
# where delta is the holonomy/one-loop correction.
# We need delta(alpha) = 2*pi*alpha for the alpha equation.
#
# From the one-loop: the correction to the INVERSE coupling is:
#   delta(1/alpha) = at_effective * (2*pi*alpha)^2 / (something)
# This requires more careful matching of conventions.

# Let me instead check: what FRACTION of the tree-level does the one-loop represent?
alpha_phys = (4*a1*PHI**4 - np.sqrt((4*a1*PHI**4)**2 - 8*np.pi)) / (4*np.pi)
theta_phys = 2 * np.pi * alpha_phys

# Compute V_eff at theta_phys
ld0_phys, ld1_phys, _, _ = compute_hodge_logdets(theta_phys)
Gamma_phys = (ld0_phys - logdet0_ref) + (ld1_phys - logdet1_ref)
S_YM_phys = 2400 * np.sin(theta_phys / 20)**2

print(f"\n  AT PHYSICAL FLUX theta = 2*pi*alpha = {theta_phys:.6f}:")
print(f"    S_YM(theta_phys) = {S_YM_phys:.8f}")
print(f"    Gamma_total(theta_phys) = {Gamma_phys:.8f}")
print(f"    Ratio Gamma/S_YM = {Gamma_phys/S_YM_phys:.6f}")
print(f"    N_gen * |Gamma|/S_YM = {N_gen * abs(Gamma_phys)/S_YM_phys:.6f}")


# =====================================================================
# STEP 4: Does the Dirac loop give 2*pi?
# =====================================================================
print("\n" + "=" * 70)
print("STEP 4: DOES THE DIRAC LOOP GIVE 2*pi?")
print("=" * 70)

# In the alpha equation: the holonomy correction is 2*pi*alpha.
# In the effective action language: at the self-consistent point,
# the one-loop correction to 1/alpha is -2*pi*alpha.
#
# So: delta(1/alpha) = -2*pi*alpha
# And: delta(1/alpha) is proportional to -N_gen * at * (something)
#
# The "something" involves the normalization between theta and alpha.
# If theta = 2*pi*alpha: delta(1/alpha) = ?? * N_gen * Gamma(2*pi*alpha)

# Let me check the RATIO more carefully.
# The tree-level: 1/alpha_0 = B = N*phi^4/K_YM = 137.08
# The correction: 2*pi*alpha = 2*pi/137 = 0.0459
# The effective: 1/alpha = B - 2*pi*alpha = 137.08 - 0.046 = 137.03

# The one-loop gives a correction to log det, not directly to 1/alpha.
# The Coleman-Weinberg potential at small field gives:
# V_CW ~ (coefficient) * theta^2 * [log(theta^2/mu^2) - constant]
# At the physical point: this gives a running coupling.

# For the LATTICE: there's no running (fixed lattice spacing).
# The one-loop correction is just Gamma(theta).
# The effective coupling: 1/g^2_eff = K_YM/g^2_tree + N_gen * at_coeff / 2...
# wait, the sign depends on statistics.

# For FERMIONS: the loop gives NEGATIVE contribution to the coupling
# (fermion loops SCREEN the charge, making alpha run to smaller values at low energy)
# delta(1/alpha)_fermion = + positive (ANTI-screening... wait, QED has screening)

# Actually in QED: fermion loops make alpha LARGER at higher energy (screening).
# At the lattice scale (UV): 1/alpha_lattice = 1/alpha_IR + (b0/2pi)*log(Lambda/mu)
# b0 > 0 for QED (fermion loops), so 1/alpha_lattice > 1/alpha_IR.
# Our tree-level = lattice scale. The physical = IR.
# 1/alpha_IR = 1/alpha_lattice - (b0/2pi)*log(Lambda/mu)

# In our framework: 1/alpha = 1/alpha_0 - 2*pi*alpha
# This looks like: 1/alpha = B - 2*pi*alpha
# Solving: alpha = (B - sqrt(B^2-8pi))/(4pi)

# The correction 2*pi*alpha comes from... what exactly?
# In the KK picture: it's the holonomy of the gauge field around the fiber.
# In the QFT picture: it would be a one-loop correction with specific coefficient.

# The coefficient at from the Dirac loop is at ~ {at:.6f}.
# Is N_gen * at related to 2*pi?
print(f"  One-loop coefficient: at = Gamma/theta^2 = {at:.6f}")
print(f"  N_gen * at = {N_gen * at:.6f}")
print(f"  2 * pi = {2*np.pi:.6f}")
print(f"  Ratio N_gen*at / (2*pi) = {N_gen * at / (2*np.pi):.6f}")
print(f"  Ratio at / (2*pi/N_gen) = {at / (2*np.pi/N_gen):.6f}")

# Check other combinations
print(f"\n  Other ratios:")
print(f"  at * N = {at * N:.6f}")
print(f"  at * N / (2*pi) = {at * N / (2*np.pi):.6f}")
print(f"  at * K_YM = {at * K_YM:.6f}")
print(f"  at * degree = {at * 12:.6f}")


print("\n" + "=" * 70)
print("EXP-574 COMPLETE")
print("=" * 70)
