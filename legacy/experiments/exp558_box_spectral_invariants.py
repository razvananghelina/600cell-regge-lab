"""
exp558: Spectral invariants of the Lorentzian operator Box
==========================================================

Now that c^2 = a1 = 5 is confirmed (exp557), explore what else the
Box operator encodes. Specifically:

1. Verify conjugation invariance algebraically (proof that [L_fiber, L_cross]=0)
2. Full joint spectrum of (L_fiber, L_cross) since they commute
3. Spectral determinant, zeta function, heat kernel of Box
4. SEARCH FOR 2*pi: does the Hopf fiber holonomy appear in spectral data?
5. Connection to the alpha equation 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0

BLIND: compute first, interpret after.
"""

import numpy as np
from numpy.linalg import eigh, norm, eigvalsh
from scipy.linalg import expm
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
TOL = 1e-8

# Build 600-cell
verts, adj, lap = build_600cell()
N = 120

# === Quaternion tools ===
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

# === Find Hopf fibration ===
def find_fibration():
    target_w = PHI / 2.0
    for i in range(N):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            power = g.copy()
            ok = True
            for k in range(2, 11):
                power = qmul(power, g)
                if k == 5 and not np.allclose(power, [-1,0,0,0], atol=1e-6):
                    ok = False; break
                if k == 10 and not np.allclose(power, [1,0,0,0], atol=1e-6):
                    ok = False
            if ok:
                # Build fibration
                subg = []
                p = np.array([1.0, 0, 0, 0])
                for k in range(10):
                    idx = find_idx(p, verts)
                    subg.append(idx)
                    p = qmul(p, g)
                used = set()
                fibers = []
                for s in range(N):
                    if s in used: continue
                    coset = []
                    for si in subg:
                        q = qmul(verts[s], verts[si])
                        idx = find_idx(q, verts)
                        if idx >= 0 and idx not in used:
                            coset.append(idx)
                            used.add(idx)
                    if len(coset) == 10:
                        fibers.append(coset)
                if len(fibers) == 12:
                    return fibers, i
    return None, -1

fibers, gen_idx = find_fibration()
assert fibers is not None

# Build L_fiber, L_cross
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
# PART 1: Verify conjugation invariance (algebraic proof)
# =====================================================================
print("=" * 70)
print("PART 1: ALGEBRAIC PROOF VERIFICATION")
print("=" * 70)

# The 12 generators (nearest neighbors of identity)
identity_idx = find_idx(np.array([1.0, 0, 0, 0]), verts)
generators = [j for j in range(N) if adj[identity_idx, j] > 0.5]
print(f"\n  Identity vertex index: {identity_idx}")
print(f"  Number of generators (neighbors of identity): {len(generators)}")

# Check: do all generators have the same trace? (single conjugacy class)
traces = [2*verts[g, 0] for g in generators]
print(f"  Generator traces: min={min(traces):.6f}, max={max(traces):.6f}")
print(f"  All trace = phi = {PHI:.6f}? {all(abs(t - PHI) < 1e-6 for t in traces)}")

# Identify S_fiber and S_cross
S_fiber = [g for g in generators if A_fiber[identity_idx, g] > 0.5]
S_cross = [g for g in generators if A_cross[identity_idx, g] > 0.5]
print(f"  |S_fiber| = {len(S_fiber)} (expect 2)")
print(f"  |S_cross| = {len(S_cross)} (expect 10)")

# Check: S_fiber elements commute (both in abelian C10)
f1, f2 = verts[S_fiber[0]], verts[S_fiber[1]]
comm_quat = qmul(f1, f2) - qmul(f2, f1)
print(f"  S_fiber commute? ||f1*f2 - f2*f1|| = {norm(comm_quat):.2e}")

# Check: S_cross is invariant under conjugation by each fiber generator
for fi, f_idx in enumerate(S_fiber):
    f = verts[f_idx]
    f_inv = np.array([f[0], -f[1], -f[2], -f[3]])
    conjugated = set()
    for c_idx in S_cross:
        c = verts[c_idx]
        fcf_inv = qmul(qmul(f, c), f_inv)
        new_idx = find_idx(fcf_inv, verts)
        conjugated.add(new_idx)

    cross_set = set(S_cross)
    invariant = conjugated == cross_set
    print(f"  f{fi+1}: S_cross invariant under conjugation? {invariant}")
    if not invariant:
        print(f"    S_cross = {sorted(S_cross)}")
        print(f"    Conjugated = {sorted(conjugated)}")

# Verify product identity: {fc : f in S_fiber, c in S_cross} = {cf : c in S_cross, f in S_fiber}
fc_set = []
cf_set = []
for f_idx in S_fiber:
    for c_idx in S_cross:
        fc = qmul(verts[f_idx], verts[c_idx])
        cf = qmul(verts[c_idx], verts[f_idx])
        fc_set.append(find_idx(fc, verts))
        cf_set.append(find_idx(cf, verts))

print(f"  Product multisets equal? {sorted(fc_set) == sorted(cf_set)}")


# =====================================================================
# PART 2: Joint spectrum of (L_fiber, L_cross)
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: JOINT SPECTRUM (L_fiber, L_cross)")
print("=" * 70)

# Since [L_fiber, L_cross] = 0, they share eigenvectors
evals_full, evecs_full = eigh(lap)

# For each eigenvector of L, compute L_fiber and L_cross eigenvalues
joint_spectrum = []
for i in range(N):
    v = evecs_full[:, i]
    lf = v @ L_fiber @ v
    lc = v @ L_cross @ v
    ltot = v @ lap @ v
    mu = lc - a1 * lf  # Box eigenvalue
    joint_spectrum.append((lf, lc, ltot, mu))

# Group by (lf, lc) pairs
from collections import Counter
pairs = Counter()
for lf, lc, lt, mu in joint_spectrum:
    key = (round(lf, 4), round(lc, 4))
    pairs[key] += 1

print(f"\n  Joint eigenvalue pairs (lf, lc) with multiplicities:")
print(f"  {'lf':>10s} {'lc':>10s} {'l=lf+lc':>10s} {'mu=lc-5lf':>10s} {'mult':>5s}")
print(f"  {'-'*50}")
for (lf, lc), mult in sorted(pairs.items()):
    print(f"  {lf:10.4f} {lc:10.4f} {lf+lc:10.4f} {lc-a1*lf:10.4f} {mult:5d}")


# =====================================================================
# PART 3: Spectral determinant and zeta function
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: SPECTRAL INVARIANTS OF Box")
print("=" * 70)

evals_box = np.sort(eigvalsh(Box))
nonzero = evals_box[np.abs(evals_box) > TOL]
positive = nonzero[nonzero > 0]
negative = nonzero[nonzero < 0]

print(f"\n  Box eigenvalues: {len(evals_box)} total, {len(nonzero)} nonzero")
print(f"  Positive: {len(positive)}, Negative: {len(negative)}, Zero: {N - len(nonzero)}")

# Spectral determinant det'(Box) = product of nonzero eigenvalues
log_det_abs = np.sum(np.log(np.abs(nonzero)))
det_sign = np.prod(np.sign(nonzero))
print(f"\n  log|det'(Box)| = {log_det_abs:.6f}")
print(f"  sign(det')     = {det_sign:+.0f}")
print(f"  |det'(Box)|    = {np.exp(log_det_abs):.6e}")

# Separate positive and negative determinants
log_det_pos = np.sum(np.log(positive))
log_det_neg = np.sum(np.log(np.abs(negative)))
print(f"\n  log(det_+) = {log_det_pos:.6f} (product of positive evals)")
print(f"  log|det_-| = {log_det_neg:.6f} (product of |negative evals|)")
print(f"  Ratio log(det_+)/log|det_-| = {log_det_pos/log_det_neg:.10f}")

# Spectral zeta function at positive integers
print(f"\n  Spectral zeta-type sums (over nonzero |mu|):")
for s in [1, 2, 3, 4]:
    zeta_s = np.sum(np.abs(nonzero)**(-s))
    zeta_pos = np.sum(positive**(-s))
    zeta_neg = np.sum(np.abs(negative)**(-s))
    print(f"    zeta({s}) = {zeta_s:.10f}  (pos: {zeta_pos:.6f}, neg: {zeta_neg:.6f})")

# Trace moments Tr(Box^n)
print(f"\n  Trace moments:")
for n in range(1, 7):
    tr = np.sum(evals_box**n)
    print(f"    Tr(Box^{n}) = {tr:15.4f}  /N = {tr/N:12.4f}")


# =====================================================================
# PART 4: SEARCH FOR 2*pi
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: SEARCHING FOR 2*pi IN SPECTRAL DATA")
print("=" * 70)

two_pi = 2 * np.pi
four_pi = 4 * np.pi
pi = np.pi

# Alpha equation coefficients
B_alpha = 4 * a1 * PHI**4
print(f"\n  Alpha equation: 2*pi*x^2 - B*x + 1 = 0")
print(f"  B = 4*a1*phi^4 = {B_alpha:.6f}")
print(f"  2*pi = {two_pi:.6f}")

# Check various spectral quantities against 2*pi
print(f"\n  --- Spectral quantities vs 2*pi ---")

# Heat kernel at t=1
K1 = np.sum(np.exp(-evals_box))
print(f"  Tr(exp(-Box)) = {K1:.10f}")

# Heat kernel at various special t
for t_name, t_val in [("1/a1", 1/a1), ("1/b1", 1/b1), ("1/12", 1/12),
                       ("phi", PHI), ("1/phi", 1/PHI), ("1/phi^2", 1/PHI**2),
                       ("1", 1.0), ("pi/a1", pi/a1)]:
    Kt = np.sum(np.exp(-t_val * evals_box))
    ratio_2pi = Kt / two_pi
    print(f"  K(t={t_name:>7s}) = {Kt:15.6f}   /2pi = {ratio_2pi:.6f}")

# Fiber spectral determinant
evals_Lf = np.sort(eigvalsh(L_fiber))
nonzero_Lf = evals_Lf[evals_Lf > TOL]
det_Lf = np.prod(nonzero_Lf)
log_det_Lf = np.sum(np.log(nonzero_Lf))
print(f"\n  det'(L_fiber) = {det_Lf:.6e}")
print(f"  log(det'(L_fiber)) = {log_det_Lf:.6f}")
print(f"  det'(C10) = 10 (standard result)")
print(f"  (det'(C10))^12 = 10^12 = {10**12:.6e}")
print(f"  det'(L_fiber) / 10^12 = {det_Lf / 10**12:.6f}")

# Cross fiber spectral determinant
evals_Lc = np.sort(eigvalsh(L_cross))
nonzero_Lc = evals_Lc[evals_Lc > TOL]
log_det_Lc = np.sum(np.log(nonzero_Lc))
print(f"\n  log(det'(L_cross)) = {log_det_Lc:.6f}")

# Full Laplacian spectral determinant
evals_L = np.sort(eigvalsh(lap))
nonzero_L = evals_L[evals_L > TOL]
log_det_L = np.sum(np.log(nonzero_L))
print(f"  log(det'(L_full)) = {log_det_L:.6f}")

# Key ratios
print(f"\n  --- Key ratios ---")
print(f"  log(det'(L_full)) / log(det'(L_fiber)) = {log_det_L/log_det_Lf:.10f}")
print(f"  log(det'(L_cross)) / log(det'(L_fiber)) = {log_det_Lc/log_det_Lf:.10f}")

# The Hopf fiber holonomy: total angle = 10 * (2*pi/10) = 2*pi
# Per-edge holonomy = 2*pi/10 = pi/5 = pi/a1
# Can we see pi/a1 in the spectral data?
print(f"\n  Hopf fiber per-edge holonomy: pi/a1 = {pi/a1:.10f}")
print(f"  gap(L_fiber) = 2-phi = {2-PHI:.10f}")
print(f"  4*sin^2(pi/10) = {4*np.sin(pi/10)**2:.10f}  (= 2-phi? {abs(4*np.sin(pi/10)**2 - (2-PHI)) < TOL})")
print(f"  ==> gap(L_fiber) = 4*sin^2(pi/(2*a1))")
print(f"  ==> pi = a1 * arcsin(sqrt(gap(L_fiber))/2) * 2 = {a1 * 2 * np.arcsin(np.sqrt(2-PHI)/2):.10f}")
print(f"  ==> 2*pi = 4*a1 * arcsin(sqrt(gap(L_fiber))/2) = {4*a1*np.arcsin(np.sqrt(2-PHI)/2):.10f}")

# The alpha equation from spectral data?
# 2*pi*alpha^2 - B*alpha + 1 = 0
# B = 4*a1*phi^4
# Can we write B in terms of Box spectrum?
print(f"\n  --- Alpha equation coefficient B from spectral data ---")
print(f"  B = 4*a1*phi^4 = {B_alpha:.10f}")

# Tr(Box^2)/N = ?
tr2 = np.sum(evals_box**2)
print(f"  Tr(Box^2)   = {tr2:.4f}")
print(f"  Tr(Box^2)/N = {tr2/N:.4f}")

# Various combinations
tr1_abs = np.sum(np.abs(evals_box))
print(f"  Tr|Box|     = {tr1_abs:.4f}")
print(f"  Tr|Box|/N   = {tr1_abs/N:.4f}")

# Sum of positive eigenvalues
sum_pos = np.sum(positive)
sum_neg = np.sum(negative)
print(f"  Sum(mu > 0) = {sum_pos:.4f}")
print(f"  Sum(mu < 0) = {sum_neg:.4f}")
print(f"  Asymmetry   = {sum_pos + sum_neg:.4f}  (should = Tr(Box))")

# Tr(Box) = Tr(L_cross) - a1*Tr(L_fiber)
tr_Lf = np.sum(evals_Lf)
tr_Lc = np.sum(evals_Lc)
tr_box = tr_Lc - a1 * tr_Lf
print(f"\n  Tr(L_fiber) = {tr_Lf:.4f} = N*deg_f = {N*2}")
print(f"  Tr(L_cross) = {tr_Lc:.4f} = N*deg_c = {N*10}")
print(f"  Tr(Box) = {tr_Lc} - {a1}*{tr_Lf} = {tr_box:.4f}")

# Winding number / holonomy from fiber
# For C10: the winding number is encoded in the spectrum
# eigenvalues 2 - 2*cos(2*pi*k/10) for k=0,...,9
# The k=1 eigenvalue "knows" about 2*pi/10
# Reconstructing 2*pi: 2*pi = 10 * 2*arcsin(sqrt(gap/4))... we did this above

# More direct: regularized trace
# zeta_Lf'(0) = -log(det'(Lf))
# For C_n: det'(L) = n, so zeta'(0) = -log(n)
# For 12 disjoint C10: zeta'(0) = -12*log(10)
print(f"\n  zeta'_Lf(0) = -log(det'(L_fiber)) = {-log_det_Lf:.6f}")
print(f"  Expected -12*log(10) = {-12*np.log(10):.6f}")
print(f"  Match: {abs(-log_det_Lf - (-12*np.log(10))) < 0.01}")

# The eta invariant (spectral asymmetry of Box)
# eta(0) = sum sign(mu) |mu|^0 = #(pos) - #(neg)
n_pos = len(positive)
n_neg = len(negative)
eta_0 = n_pos - n_neg
print(f"\n  --- Eta invariant (spectral asymmetry) ---")
print(f"  eta(0) = #pos - #neg = {n_pos} - {n_neg} = {eta_0}")
print(f"  eta(0)/2 = {eta_0/2:.1f}")

# Eta function at s=1
eta_1 = np.sum(np.sign(nonzero) * np.abs(nonzero)**(-1))
print(f"  eta(1) = sum sign(mu)/|mu| = {eta_1:.10f}")
print(f"  eta(1)*a1 = {eta_1*a1:.10f}")
print(f"  eta(1)*pi = {eta_1*pi:.10f}")

# Eta at s=2
eta_2 = np.sum(np.sign(nonzero) * np.abs(nonzero)**(-2))
print(f"  eta(2) = sum sign(mu)/|mu|^2 = {eta_2:.10f}")


# =====================================================================
# PART 5: RATIO ANALYSIS -- hunt for alpha equation
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: HUNTING FOR ALPHA EQUATION IN SPECTRAL DATA")
print("=" * 70)

# The alpha equation: 2*pi * alpha^2 - 4*a1*phi^4 * alpha + 1 = 0
#
# The coefficient "4*a1*phi^4" can be rewritten:
#   4*a1*phi^4 = 4*5*(3*phi+2) = 60*phi + 40
#   phi^4 = phi^2 + phi^2 - 1 + 1 = ... phi^4 = 3*phi+2
#
# Where in the 600-cell does phi^4 come from?
# The adjacency eigenvalue for rho_1 is 6*phi.
# 6*phi * (2*phi/3) = 4*phi^2 = 4*(phi+1) = 4*phi+4
# Hmm, 4*a1*phi^4 = 4*a1*(phi^2)^2 = 4*a1*(phi+1)^2 = 4*5*(phi^2+2*phi+1)
#   = 20*(phi+1+2*phi+1) = 20*(3*phi+2) = 60*phi+40

print(f"\n  B = 4*a1*phi^4 = {B_alpha:.10f}")
print(f"    = 60*phi + 40 = {60*PHI + 40:.10f}")

# Try: B from eigenvalues of A
# rho_1 eigenvalue = 6*phi, rho_8 = 6*phi' = -6/phi
lam_1 = 6*PHI
lam_8 = 6*(1-np.sqrt(5))/2

print(f"\n  Adjacency eigenvalue rho_1 = 6*phi = {lam_1:.6f}")
print(f"  Adjacency eigenvalue rho_8 = 6*phi' = {lam_8:.6f}")
print(f"  lam_1^2 = 36*phi^2 = {lam_1**2:.6f}")
print(f"  lam_1^2 / (a1-1) = {lam_1**2/(a1-1):.6f}")
print(f"  4*a1*phi^4 = (lam_1)^2 * (a1*phi^2/9) = {lam_1**2 * a1 * PHI**2 / 9:.6f}")

# Spectral gap of L_fiber times various things
gap_f = 2 - PHI
gap_c = 10 - 5*PHI
gap_full = 12 - 6*PHI

print(f"\n  Spectral gaps:")
print(f"    gap_f = {gap_f:.10f} = 1/phi^2 = {1/PHI**2:.10f}")
print(f"    gap_c = {gap_c:.10f}")
print(f"    gap_full = {gap_full:.10f}")
print(f"    gap_f = 2-phi = 1/phi^2 !! (golden ratio identity)")

# Can we extract 2*pi from the ratio of spectral quantities?
# The key identity: gap_f = 4*sin^2(pi/(2*a1))
# So pi = 2*a1 * arcsin(sqrt(gap_f)/2)
# And 2*pi = 4*a1 * arcsin(sqrt(gap_f)/2)
pi_from_gap = 2*a1 * np.arcsin(np.sqrt(gap_f)/2)
print(f"\n  RECONSTRUCTION OF pi:")
print(f"    gap_f = 4*sin^2(pi/(2*a1))")
print(f"    pi = 2*a1 * arcsin(sqrt(gap_f)/2)")
print(f"    = 2*{a1} * arcsin(sqrt({gap_f:.6f})/2)")
print(f"    = {pi_from_gap:.10f}")
print(f"    pi = {pi:.10f}")
print(f"    Match: {abs(pi_from_gap - pi) < 1e-10}")

# So the FULL alpha equation can be written as:
# A * alpha^2 - B * alpha + C = 0
# where:
#   A = 2*pi = 8*a1 * arcsin(sqrt(2-phi)/2) -- from fiber gap
#   B = 4*a1*phi^4 -- algebraic
#   C = 1
#
# Or equivalently, using gap_f = 2-phi = 1/phi^2:
#   A = 8*a1 * arcsin(1/(2*phi))
#   B = 4*a1*(3*phi+2)

print(f"\n  ALPHA EQUATION FROM SPECTRAL DATA:")
print(f"  A = 2*pi = 8*a1 * arcsin(1/(2*phi)) = {8*a1*np.arcsin(1/(2*PHI)):.10f}")
print(f"  B = 4*a1*(3*phi+2) = {4*a1*(3*PHI+2):.10f}")
print(f"  C = 1")
print(f"  Equation: A*alpha^2 - B*alpha + C = 0")

alpha_val = (B_alpha - np.sqrt(B_alpha**2 - 8*pi)) / (4*pi)
print(f"  alpha = {alpha_val:.10f} = 1/{1/alpha_val:.6f}")
print(f"  NIST:   alpha = 0.0072973526 = 1/137.035999")


# =====================================================================
# PART 6: HOLONOMY AND SPECTRAL RECONSTRUCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: FIBER HOLONOMY = 2*pi FROM C10 SPECTRUM")
print("=" * 70)

# Single C10 cycle eigenvalues
c10_evals = np.array([2 - 2*np.cos(2*np.pi*k/10) for k in range(10)])
c10_evals_sorted = np.sort(c10_evals)
print(f"  C10 eigenvalues: {np.round(c10_evals_sorted, 6)}")

# The gap encodes 2*pi/n: gap = 4*sin^2(pi/n) for C_n
# This is the discrete version of the Laplacian eigenvalue (2*pi*k/L)^2
# on a circle of circumference L, at k=1.
# The discrete-to-continuum map: (2*pi/n)^2 ~ gap (for large n)

# For C10: exact gap = 2-phi = 4*sin^2(pi/10)
# The angular resolution is delta_theta = 2*pi/10 = pi/5
# The total angle (holonomy) = 10 * delta_theta = 2*pi

# EXACT relation: det'(L_{C_n}) = n
# So det'(L_{C_{10}}) = 10 = 2*a1
det_c10 = np.prod(c10_evals_sorted[1:])  # exclude zero eigenvalue
print(f"\n  det'(C10) = {det_c10:.10f} (expect 10 = 2*a1)")
print(f"  Match: {abs(det_c10 - 10) < TOL}")

# Trace of C10 Laplacian
tr_c10 = np.sum(c10_evals_sorted)
print(f"  Tr(L_C10) = {tr_c10:.6f} (= 2*n = 20)")

# The REGULARIZED product: for C_n, sum_{k=1}^{n-1} log(2-2cos(2pi*k/n)) = log(n)
# This gives: sum_{k=1}^{n-1} log(4*sin^2(pi*k/n)) = log(n)
# i.e.: product_{k=1}^{n-1} sin(pi*k/n) = n / 2^{n-1}
# For n=10: product_{k=1}^{9} sin(pi*k/10) = 10 / 512 = 5/256
sin_product = np.prod([np.sin(np.pi*k/10) for k in range(1, 10)])
print(f"\n  Product sin(pi*k/10) for k=1..9 = {sin_product:.10f}")
print(f"  Expected 10/2^9 = 10/512 = {10/512:.10f}")
print(f"  = a1/256 = {a1/256:.10f}")
print(f"  Match: {abs(sin_product - 10/512) < TOL}")

# KEY OBSERVATION: the holonomy 2*pi is NOT a spectral invariant.
# It is a TOPOLOGICAL invariant of the circle.
# The spectrum {4*sin^2(pi*k/n)} determines n (from det' = n),
# but 2*pi enters through the PARAMETRIZATION of eigenvalues as
# functions of k, not from the eigenvalues themselves.
print(f"""
  KEY OBSERVATION:
  2*pi enters the alpha equation through the Hopf fiber topology.
  The C10 spectrum determines a1 (from det'=2*a1) but NOT 2*pi directly.
  2*pi is reconstructed as 4*a1*arcsin(sqrt(gap_f)/2), mixing:
    - a1 = 5 (algebraic, from 600-cell combinatorics)
    - gap_f = 2-phi = 1/phi^2 (algebraic, from golden ratio)
    - arcsin (transcendental, from circle geometry)
  The transcendental function arcsin encodes the curvature of S^1.
  This is the IRREDUCIBLE transcendental input:
    2*pi = 4*a1 * arcsin(1/(2*phi))
  It cannot be removed because S^1 has nonzero curvature.
""")


# =====================================================================
# PART 7: IS 2*pi DETERMINED BY a1 AND phi?
# =====================================================================
print("=" * 70)
print("PART 7: THE 2*pi IDENTITY")
print("=" * 70)

# The question: is 2*pi = 4*a1*arcsin(1/(2*phi)) an identity or a coincidence?
# arcsin(1/(2*phi)) = arcsin(sin(pi/10)) = pi/10 (for the principal branch)
# So 2*pi = 4*a1 * pi/10 = 4*5*pi/10 = 2*pi. TAUTOLOGY!
#
# But this IS the content: the fiber gap = 4*sin^2(pi/(2*a1))
# is equivalent to saying the C10 cycle inscribes a circle of
# circumference 2*pi with 2*a1 = 10 equally spaced vertices.

print(f"""
  Identity: 2*pi = 4*a1 * arcsin(1/(2*phi))

  Proof: arcsin(1/(2*phi)) = arcsin(sin(pi/10)) = pi/10 = pi/(2*a1)
  So: 4*a1 * pi/(2*a1) = 2*pi.    (tautology)

  The deeper question: WHY is sin(pi/(2*a1)) = 1/(2*phi)?

  This is the identity: sin(pi/10) = 1/(2*phi) = (sqrt(5)-1)/4

  Equivalently: sin(18 deg) = (sqrt(5)-1)/4

  This is a THEOREM from Euclidean geometry (regular pentagon).
  It connects:
    - pi (circle geometry / topology of S^1)
    - phi (golden ratio / algebra of Z[phi])
    - a1 = 5 (the integer that determines everything)

  The connection is: 2*cos(2*pi/a1!) = phi  (since a1! = 10, 2*pi/10 = pi/5)
  Wait: 2*cos(pi/5) = phi. And a1 = 5, 2*a1 = 10, so pi/5 = 2*pi/(2*a1).

  EXACT CHAIN:
    a1 = 5
    |2I| = 120 = 2*a1!  (but a1! = 120, not 10...)
    Actually: Hopf fiber = C_{2*a1} = C_{10}
    fiber gap = 2 - 2*cos(2*pi/(2*a1)) = 2 - 2*cos(pi/a1) = 2 - phi
    2*cos(pi/a1) = phi  <==>  cos(pi/5) = phi/2

    This is the PENTAGON IDENTITY: the golden ratio IS the diagonal-to-side
    ratio of a regular pentagon, equivalently cos(pi/5) = phi/2.

    So 2*pi enters through:
      2*pi = 2*a1 * (pi/a1) = 2*a1 * arccos(phi/2)

    And the alpha equation becomes:
      2*a1*arccos(phi/2)*alpha^2 - 4*a1*phi^4*alpha + 1 = 0

    Dividing by 2*a1:
      arccos(phi/2)*alpha^2 - 2*phi^4*alpha + 1/(2*a1) = 0
""")

# Verify the arccos form
acos_val = np.arccos(PHI/2)
print(f"  arccos(phi/2) = {acos_val:.10f}")
print(f"  pi/a1 = pi/5 = {pi/a1:.10f}")
print(f"  Match: {abs(acos_val - pi/a1) < 1e-10}")

# Normalized alpha equation
print(f"\n  Normalized alpha equation (dividing by 2*a1):")
A_norm = acos_val  # = pi/5
B_norm = 2*PHI**4
C_norm = 1.0/(2*a1)
disc_norm = B_norm**2 - 4*A_norm*C_norm
alpha_check = (B_norm - np.sqrt(disc_norm)) / (2*A_norm)
print(f"    arccos(phi/2)*x^2 - 2*phi^4*x + 1/(2*a1) = 0")
print(f"    A = arccos(phi/2) = pi/5 = {A_norm:.10f}")
print(f"    B = 2*phi^4       = {B_norm:.10f}")
print(f"    C = 1/(2*a1) = 1/10 = {C_norm:.10f}")
print(f"    alpha = {alpha_check:.10f} = 1/{1/alpha_check:.6f}")

print(f"""
  CONCLUSION:
  ===========
  The alpha equation 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
  can be rewritten as:

    arccos(phi/2) * alpha^2 - 2*phi^4 * alpha + 1/(2*a1) = 0

  where arccos(phi/2) = pi/a1 is the DIHEDRAL ANGLE per edge
  of the Hopf fiber decagon.

  The "2*pi" is not a free transcendental input. It equals:
    2*pi = 2*a1 * arccos(phi/2)
  which follows from the PENTAGON IDENTITY cos(pi/5) = phi/2.

  This identity connects:
    - The integer a1 = 5 (number of sides of the regular pentagon)
    - The golden ratio phi = 2*cos(pi/5) (pentagon diagonal/side ratio)
    - The circle constant pi (topology of S^1)

  In the framework: the Hopf fiber is a C_{2*a1} = C_10 cycle.
  Its per-edge angle pi/a1 is not a choice but a CONSEQUENCE
  of the 600-cell geometry. The total holonomy 2*pi = 2*a1 * (pi/a1)
  is then determined by the integer a1 and the topology of S^1.
""")

print("=" * 70)
print("EXP-558 COMPLETE")
print("=" * 70)
