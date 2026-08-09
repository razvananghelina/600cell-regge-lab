"""
exp463_newton_spectral.py
==========================
NEWTON'S CONSTANT FROM SPECTRAL DATA ON 600-CELL

GOAL: Derive G (or equivalently m_e/M_P = alpha^{4*phi^2}) from
spectral invariants of the Hodge Laplacians on the 600-cell.

APPROACH:
  1. Compute FULL D^2 spectrum (2640 eigenvalues)
  2. Extract spectral invariants: det', zeta, heat kernel traces
  3. Use Connes relation G = pi/(f_2*Lambda^2) with framework f_2
  4. Test KK reduction: G_4D from G_3D = 1/12 (DERIVED, exp424)
  5. Look for z_Planck = 4*phi^2 in spectral ratios

KEY IDENTITY (exp462): z_Planck = 2*D^2 - d_ST  (DERIVED)
where D^2 = a1 + sqrt(a1) (TQFT total quantum dim).

KNOWN: z_Planck DERIVED. m_e is the ONE free parameter.
QUESTION: Can the spectral action ALSO determine the absolute scale?

Author: Razvan-Constantin Anghelina
Date: March 2026
"""

import numpy as np
from itertools import permutations, product
from collections import defaultdict
from math import sqrt, pi, factorial, log, log10, exp

# =============================================================================
# CONSTANTS
# =============================================================================
a1 = 5
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
PHI_c = (1 - sqrt(a1)) / 2
N = factorial(a1)  # 120
degree = 2 * b1    # 12
N_gen = 3
N_eig = 9
rank_E8 = 8

# Couplings
sin2tW = b1 / (a1**2 + 1)  # 6/26
alpha_s = 1 / (2 * PHI**3)
alpha_em = (-(-4*a1*PHI**4) - sqrt((-4*a1*PHI**4)**2 - 4*2*pi*1)) / (2*2*pi)

# Seeley-DeWitt
A0 = 2*a1 + 1        # 11
A1 = 2*a1**2+2*a1+2  # 62
A2 = 6*a1**2+15*a1+8 # 233
c0_fw = 2*N*A0       # 2640
c1_fw = 2*N*A1       # 14880
c2_fw = 2*N*A2       # 55920

# Hierarchy
z_Planck = 4 * PHI**2  # = 2*D^2 - d_ST = 10.472
D_sq = a1 + sqrt(a1)
d_ST = 4

# Physical
m_e_GeV = 0.51099895e-3
M_P_GeV = 1.22089e19
ln_inv_alpha = -log(alpha_em)

print("=" * 72)
print("EXP-463: NEWTON'S CONSTANT FROM SPECTRAL DATA")
print("=" * 72)
print()
print(f"z_Planck = 4*phi^2 = 2*D^2-d_ST = {z_Planck:.6f}")
print(f"m_e/M_P = alpha^z_Planck = {alpha_em**z_Planck:.6e}")
print(f"Observed: {m_e_GeV/M_P_GeV:.6e}")
print()

# =============================================================================
# PART 1: BUILD 600-CELL + FULL D^2 SPECTRUM
# =============================================================================
print("=" * 72)
print("PART 1: BUILD 600-CELL AND COMPUTE SPECTRUM")
print("=" * 72)

# Build vertices (120)
verts_set = set()
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0]*4; v[i] = s; verts_set.add(tuple(v))
for s0 in [0.5,-0.5]:
    for s1 in [0.5,-0.5]:
        for s2 in [0.5,-0.5]:
            for s3 in [0.5,-0.5]:
                verts_set.add((s0,s1,s2,s3))
base = [PHI/2, 0.5, 1/(2*PHI), 0.0]
even_perms = [p for p in permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0]
for perm in even_perms:
    coords = [base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if abs(coords[i])>1e-12]
    for signs in product([1,-1], repeat=len(nz)):
        v = list(coords)
        for idx, s in zip(nz, signs): v[idx] *= s
        verts_set.add(tuple(round(x,10) for x in v))

verts = np.array(sorted(verts_set))
Nv = len(verts)
assert Nv == 120

# Build edges (720)
dots = verts @ verts.T
edge_thresh = PHI / 2
edges = [(i,j) for i in range(Nv) for j in range(i+1,Nv) if abs(dots[i,j]-edge_thresh)<0.001]
Ne = len(edges)
assert Ne == 720

edge_to_idx = {}
for idx, (i,j) in enumerate(edges):
    edge_to_idx[(i,j)] = idx; edge_to_idx[(j,i)] = idx

adj_list = defaultdict(set)
for i,j in edges: adj_list[i].add(j); adj_list[j].add(i)

# Build triangles (1200)
triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j: triangles.append((i,j,k))
Nf = len(triangles)
assert Nf == 1200

face_to_idx = {f: idx for idx, f in enumerate(triangles)}

# Build tetrahedra (600)
tetrahedra = []
for i in range(Nv):
    ni = adj_list[i]
    for j in ni:
        if j > i:
            cij = ni & adj_list[j]
            for k in cij:
                if k > j:
                    for l in cij & adj_list[k]:
                        if l > k: tetrahedra.append((i,j,k,l))
Nc = len(tetrahedra)
assert Nc == 600

print(f"  f-vector: ({Nv}, {Ne}, {Nf}, {Nc}), Euler={Nv-Ne+Nf-Nc}")

# Boundary operators
d0 = np.zeros((Ne, Nv))
for e_idx, (i,j) in enumerate(edges):
    d0[e_idx,i] = -1; d0[e_idx,j] = 1

d1 = np.zeros((Nf, Ne))
for f_idx, (i,j,k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i,j)]] = 1
    d1[f_idx, edge_to_idx[(j,k)]] = 1
    d1[f_idx, edge_to_idx[(i,k)]] = -1

d2 = np.zeros((Nc, Nf))
for c_idx, (i,j,k,l) in enumerate(tetrahedra):
    for face, sign in [((j,k,l),1),((i,k,l),-1),((i,j,l),1),((i,j,k),-1)]:
        f_idx = face_to_idx.get(face)
        if f_idx is not None: d2[c_idx, f_idx] = sign

assert np.max(np.abs(d1@d0)) < 1e-10
assert np.max(np.abs(d2@d1)) < 1e-10

# Hodge Laplacians
print("  Computing Hodge Laplacians and eigenvalues...")
Delta_0 = d0.T @ d0
Delta_1 = d0 @ d0.T + d1.T @ d1
Delta_2 = d1 @ d1.T + d2.T @ d2
Delta_3 = d2 @ d2.T

evals_0 = np.sort(np.linalg.eigvalsh(Delta_0))
evals_1 = np.sort(np.linalg.eigvalsh(Delta_1))
evals_2 = np.sort(np.linalg.eigvalsh(Delta_2))
evals_3 = np.sort(np.linalg.eigvalsh(Delta_3))

all_evals = np.sort(np.concatenate([evals_0, evals_1, evals_2, evals_3]))
dim_H = len(all_evals)
assert dim_H == 2640

# Zero modes (Betti numbers)
tol = 1e-8
b = [np.sum(ev < tol) for ev in [evals_0, evals_1, evals_2, evals_3]]
print(f"  Betti numbers: {b} = (1,0,0,1) for S^3")
assert b == [1,0,0,1]

# Seeley-DeWitt verification
c0 = dim_H
c1 = np.sum(all_evals)
c2 = np.sum(all_evals**2) / 2
print(f"\n  Seeley-DeWitt coefficients:")
print(f"    c0 = {c0} (fw: {c0_fw})")
print(f"    c1 = {c1:.2f} (fw: {c1_fw})")
print(f"    c2 = {c2:.2f} (fw: {c2_fw})")

# =============================================================================
# PART 2: SPECTRAL INVARIANTS PER FORM DEGREE
# =============================================================================
print()
print("=" * 72)
print("PART 2: SPECTRAL INVARIANTS PER FORM DEGREE")
print("=" * 72)
print()

for k, (evs, label) in enumerate([(evals_0,"0-forms (120)"),
                                     (evals_1,"1-forms (720)"),
                                     (evals_2,"2-forms (1200)"),
                                     (evals_3,"3-forms (600)")]):
    nz_evs = evs[evs > tol]
    n_zero = len(evs) - len(nz_evs)
    lam_min = nz_evs[0] if len(nz_evs) > 0 else 0
    lam_max = nz_evs[-1] if len(nz_evs) > 0 else 0
    tr = np.sum(nz_evs)
    tr2 = np.sum(nz_evs**2)
    log_det = np.sum(np.log(nz_evs)) if len(nz_evs) > 0 else 0

    print(f"  Delta_{k} [{label}]:")
    print(f"    Zero modes: {n_zero} (Betti b_{k} = {b[k]})")
    print(f"    Nonzero: {len(nz_evs)}")
    print(f"    Range: [{lam_min:.6f}, {lam_max:.6f}]")
    print(f"    Tr = {tr:.4f},  Tr^2 = {tr2:.4f}")
    print(f"    log det' = {log_det:.6f}")
    print()

# Full spectrum statistics
nz_all = all_evals[all_evals > tol]
log_det_full = np.sum(np.log(nz_all))
lam_max_full = nz_all[-1]
lam_min_full = nz_all[0]

print(f"  FULL D^2 spectrum (2638 nonzero):")
print(f"    Range: [{lam_min_full:.6f}, {lam_max_full:.6f}]")
print(f"    log det'(D^2) = {log_det_full:.6f}")
print(f"    det'(D^2) = exp({log_det_full:.2f})")

# =============================================================================
# PART 3: FRAMEWORK NUMBER SEARCH IN SPECTRAL DATA
# =============================================================================
print()
print("=" * 72)
print("PART 3: FRAMEWORK NUMBER SEARCH")
print("=" * 72)
print()

# Key quantities from the spectrum
log_dets = []
for k, evs in enumerate([evals_0, evals_1, evals_2, evals_3]):
    nz = evs[evs > tol]
    ld = np.sum(np.log(nz)) if len(nz) > 0 else 0
    log_dets.append(ld)
    print(f"  log det'(Delta_{k}) = {ld:.6f}")

# Analytic torsion
log_T = 0.5 * sum((-1)**k * k * log_dets[k] for k in range(4))
print(f"\n  Analytic torsion: log T = {log_T:.10f} (should be 0 for S^3)")

# Spectral ratios
print(f"\n  Spectral ratios:")
for i in range(4):
    for j in range(i+1,4):
        if log_dets[j] != 0 and log_dets[i] != 0:
            r = log_dets[i] / log_dets[j]
            print(f"    log det'(D_{i}) / log det'(D_{j}) = {r:.6f}")

# Maximum eigenvalue analysis
print(f"\n  Maximum eigenvalues:")
maxes = [evs[-1] for evs in [evals_0, evals_1, evals_2, evals_3]]
for k, mx in enumerate(maxes):
    print(f"    max(Delta_{k}) = {mx:.6f}")

print(f"\n  Lambda_max = max over all = {lam_max_full:.6f}")
print(f"  sqrt(Lambda_max) = {sqrt(lam_max_full):.6f}")

# Check framework numbers near Lambda_max
fw_nums = {
    '2*degree': 2*degree,
    'a1^2': a1**2,
    '4*a1': 4*a1,
    'N_gen*degree': N_gen*degree,
    'degree+a1': degree+a1,
    'b1*phi^2': b1*PHI**2,
    'N_eig*phi^2': N_eig*PHI**2,
}
print(f"\n  Framework numbers near Lambda_max = {lam_max_full:.4f}:")
for name, val in sorted(fw_nums.items(), key=lambda x: abs(x[1]-lam_max_full)):
    print(f"    {name:>15} = {val:10.4f}  (diff: {val-lam_max_full:+.4f})")

# =============================================================================
# PART 4: CONNES SPECTRAL ACTION FOR G
# =============================================================================
print()
print("=" * 72)
print("PART 4: CONNES SPECTRAL ACTION - G FROM f_2")
print("=" * 72)
print()

# In Connes: G = pi / (f_2 * Lambda^2 * normalization)
# For SM: 1/(16*pi*G) = f_2*Lambda^2 * c_2_Yuk / (96*pi^2)
# where c_2_Yuk = Tr(Y^dag*Y) involves Yukawa couplings
# In our framework: Tr(D_F^2) = 8 = rank(E8) (from McKay)
# So c_2_Yuk ~ 8 (framework value)

# Test: various f_2 ansatze
print("Testing f_2 ansatze for G = pi/(f_2*Lambda^2*c_2_norm):")
print()

# For each f_2 and Lambda, compute G and z_Planck_pred
# z_Planck_pred = -log(m_e/M_P)/log(1/alpha) if G is determined

# Natural cutoffs:
cutoffs = {
    'Lambda_max': lam_max_full,
    'Tr/dim': c1 / c0,   # average eigenvalue
    'degree': float(degree),
    'phi^4': PHI**4,
    'phi^2': PHI**2,
}

f2_candidates = {
    'N^2/(2pi)': N**2/(2*pi),
    'c0=dim(H)': float(c0),
    'c1/c0': c1/c0,
    'N*A0': float(N*A0),
    'N*A1': float(N*A1),
    'N^3/(2pi)': N**3/(2*pi),
    'N^4/(2pi)=f4': N**4/(2*pi),
}

print(f"  {'f_2':>15}  {'Lambda^2':>12}  {'f2*L2':>12}  {'G=pi/f2L2':>12}  {'z_pred':>8}")
for f2name, f2val in f2_candidates.items():
    for Lname, Lval in cutoffs.items():
        f2L2 = f2val * Lval
        G_pred = pi / f2L2
        # M_P = 1/sqrt(G), m_e/M_P = alpha^z
        # z = log(m_e/M_P) / log(alpha) = log(sqrt(G)*m_e) / log(alpha)
        # In natural units where the lattice = Planck scale:
        # If G_pred is in lattice units, M_P_lattice = 1/sqrt(G_pred)
        # m_e_lattice = M_P_lattice * alpha^z_Planck
        # The PREDICTED z_Planck from the spectral action:
        # G = pi/(f2*L2) should give z_Planck = log(1/sqrt(G_pred)) / log(1/alpha) ... no
        # Actually the RATIO f_4/f_2 determines the scalaron mass, not z_Planck directly.
        # z_Planck is already DERIVED algebraically. The question is whether f_2*Lambda^2
        # gives a self-consistent G.
        #
        # Just report the numbers:
        if f2L2 > 0:
            print(f"  {f2name:>15}  {Lname:>12}  {f2L2:12.2f}  {G_pred:12.6e}  ---")

# =============================================================================
# PART 5: KK REDUCTION FROM 3D TQFT
# =============================================================================
print()
print("=" * 72)
print("PART 5: KK REDUCTION FROM G_3D")
print("=" * 72)
print()

# From exp424: G_3D = 1/(4k) = 1/12 (DERIVED)
# KK: G_4D = G_3D * L_KK where L_KK is the compactification length
# From alpha*alpha' = 1/(2*pi): the KK circle has "volume" 1/(2*pi)
# From the spectral action: alpha = e^2/(4*pi), KK radius R links to alpha

G_3D = 1.0 / (4 * (a1-2))  # 1/12
print(f"  G_3D = 1/(4k) = 1/12 = {G_3D:.6f}  [DERIVED, exp424]")
print()

# Various KK radii
kk_candidates = {
    'R=1/(2pi)': 1/(2*pi),
    'R=alpha': alpha_em,
    'R=alpha_s': alpha_s,
    'R=1/N': 1.0/N,
    'R=1/(degree*phi)': 1/(degree*PHI),
    'R=alpha^2': alpha_em**2,
}

print(f"  {'KK radius':>20}  {'R':>12}  {'G_4D':>12}  {'M_P/M_latt':>12}  {'z_eff':>8}")
for name, R in kk_candidates.items():
    G_4D = G_3D * 2*pi*R
    M_P_ratio = 1/sqrt(G_4D)
    z_eff = log(M_P_ratio) / log(1/alpha_em) if M_P_ratio > 1 else -log(M_P_ratio)/log(1/alpha_em)
    print(f"  {name:>20}  {R:12.6e}  {G_4D:12.6e}  {M_P_ratio:12.4f}  {z_eff:8.4f}")

print()
print(f"  Target: z_Planck = {z_Planck:.4f}")
print(f"  Need G_4D = alpha^(2*z_Planck) = {alpha_em**(2*z_Planck):.6e}")
print()

# Solve: G_3D * 2*pi*R = alpha^(2*z_Planck)
# R = alpha^(2*z_Planck) / (2*pi*G_3D) = alpha^{20.944} / (pi/6)
R_needed = alpha_em**(2*z_Planck) / (2*pi*G_3D)
print(f"  Required R_KK = {R_needed:.6e}")
print(f"  = alpha^{log(R_needed)/log(alpha_em):.4f}")
z_R = log(R_needed)/log(alpha_em)
print(f"  So R_KK = alpha^{z_R:.4f}")
print()

# Is z_R a framework number?
print(f"  Framework numbers near z_R = {z_R:.4f}:")
fw_z = {
    '2*z_Planck': 2*z_Planck,
    '2*z_Planck-2': 2*z_Planck-2,
    '2*z_Planck+ln(6/pi)/lna': 2*z_Planck + log(6/pi)/log(alpha_em),
    'z_Planck+a1': z_Planck + a1,
    'z_Planck*2-ln(2pi*G3)/lna': 2*z_Planck - log(2*pi*G_3D)/log(alpha_em),
}
for name, val in sorted(fw_z.items(), key=lambda x: abs(x[1]-z_R)):
    print(f"    {name:>35} = {val:10.4f}  (diff: {val-z_R:+.6f})")

# =============================================================================
# PART 6: SPECTRAL ZETA AT SPECIAL VALUES
# =============================================================================
print()
print("=" * 72)
print("PART 6: SPECTRAL ZETA FUNCTION")
print("=" * 72)
print()

# zeta_D(s) = sum_{lambda>0} lambda^{-s}
# For D^2, using eigenvalues mu_i:
# zeta(s) = sum_{mu_i > 0} mu_i^{-s}

def spectral_zeta(evals, s):
    """Compute spectral zeta function at s."""
    nz = evals[evals > tol]
    return np.sum(nz**(-s))

# Per-form and total
print(f"  {'s':>6}  {'zeta_0':>12}  {'zeta_1':>12}  {'zeta_2':>12}  {'zeta_3':>12}  {'zeta_total':>12}")
for s in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, -0.5, -1.0]:
    z0 = spectral_zeta(evals_0, s)
    z1 = spectral_zeta(evals_1, s)
    z2 = spectral_zeta(evals_2, s)
    z3 = spectral_zeta(evals_3, s)
    zt = spectral_zeta(all_evals, s)
    print(f"  {s:6.1f}  {z0:12.4f}  {z1:12.4f}  {z2:12.4f}  {z3:12.4f}  {zt:12.4f}")

# Special: zeta(0) = regularized number of eigenvalues
# zeta'(0) = -log det' (for Ray-Singer torsion)
print()

# Look for z_Planck in ratios
print("  Searching for z_Planck = %.4f in spectral zeta ratios:" % z_Planck)
s_vals = [0.5, 1.0, 1.5, 2.0, 3.0]
for s1 in s_vals:
    for s2 in s_vals:
        if s1 != s2:
            zt1 = spectral_zeta(all_evals, s1)
            zt2 = spectral_zeta(all_evals, s2)
            ratio = zt1/zt2 if zt2 != 0 else 0
            if abs(ratio - z_Planck) < 0.5:
                print(f"    zeta({s1})/zeta({s2}) = {ratio:.6f}  (gap: {ratio-z_Planck:+.6f})")

# =============================================================================
# PART 7: HEAT KERNEL AT FRAMEWORK TIMES
# =============================================================================
print()
print("=" * 72)
print("PART 7: HEAT KERNEL AT FRAMEWORK TIMES")
print("=" * 72)
print()

def heat_kernel_trace(evals, t):
    return np.sum(np.exp(-t * evals))

# Natural times
times = {
    '1/phi^2': 1/PHI**2,
    '1/degree': 1.0/degree,
    '1/a1': 1.0/a1,
    '1/phi': 1/PHI,
    '1': 1.0,
    'phi': PHI,
    'phi^2': PHI**2,
    'a1': float(a1),
    'degree': float(degree),
    'a1^2': float(a1**2),
    'N_eig*phi': N_eig*PHI,
}

print(f"  {'time':>12}  {'K(t)':>14}  {'K(t)/c0':>10}  {'log K(t)':>12}")
for name, t in sorted(times.items(), key=lambda x: x[1]):
    K = heat_kernel_trace(all_evals, t)
    print(f"  {name:>12}  {K:14.6f}  {K/c0:10.6f}  {log(K):12.6f}")

# Special: K(t) at t where K = specific framework values
print()
print("  Time t* where K(t*) = framework number:")
# K(t) is monotonically decreasing from c0=2640 to 2 (zero modes)
# Find t where K = N, N/2, etc.
from scipy.optimize import brentq

def K_minus_target(t, target):
    return heat_kernel_trace(all_evals, t) - target

targets_K = {
    'N=120': 120.0,
    'c0/a1=528': c0/a1,
    'c0/degree=220': c0/degree,
    'a1^2=25': float(a1**2),
    'degree=12': float(degree),
}

for name, target_val in targets_K.items():
    try:
        t_star = brentq(lambda t: K_minus_target(t, target_val), 0.001, 100)
        print(f"    K(t*) = {name}: t* = {t_star:.6f}")
        # Check if t* is a framework number
        for fname, fval in [('1/phi', 1/PHI), ('1/a1', 1.0/a1), ('phi/degree', PHI/degree),
                            ('1/phi^2', 1/PHI**2), ('alpha_s', alpha_s), ('1/b1', 1.0/b1)]:
            if abs(t_star - fval) < 0.1*fval:
                print(f"         Near {fname} = {fval:.6f} (gap: {(t_star-fval)/fval*100:+.2f}%)")
    except:
        pass

# =============================================================================
# PART 8: RATIO log det' / ln(alpha)
# =============================================================================
print()
print("=" * 72)
print("PART 8: SPECTRAL DETERMINANT RATIOS")
print("=" * 72)
print()

# Per-form log det'
print("  Per-form log det' and ratios with ln(alpha):")
for k, ld in enumerate(log_dets):
    z = ld / log(alpha_em) if ld != 0 else 0
    print(f"    log det'(Delta_{k}) = {ld:12.6f}   / ln(alpha) = {z:12.6f}")

print()
# Combined determinants
log_det_01 = log_dets[0] + log_dets[1]
log_det_23 = log_dets[2] + log_dets[3]
log_det_even = log_dets[0] + log_dets[2]
log_det_odd = log_dets[1] + log_dets[3]

combos = {
    'D0+D1': log_det_01,
    'D2+D3': log_det_23,
    'even (D0+D2)': log_det_even,
    'odd (D1+D3)': log_det_odd,
    'D0-D3': log_dets[0] - log_dets[3],
    'D1-D2': log_dets[1] - log_dets[2],
    'full': log_det_full,
    '(D0*D3)/(D1*D2)': log_dets[0]+log_dets[3]-log_dets[1]-log_dets[2],
}

print(f"  {'Combination':>20}  {'log det':>12}  {'/ln(alpha)':>12}  {'near z_Planck?':>15}")
for name, val in combos.items():
    z = val / log(alpha_em) if val != 0 else 0
    near = abs(z - z_Planck) < 1.0 or abs(z - 2*z_Planck) < 1.0
    tag = f"z_Pl gap={z-z_Planck:+.4f}" if abs(z-z_Planck) < 2.0 else ""
    if abs(z - 2*z_Planck) < 2.0: tag = f"2*z_Pl gap={z-2*z_Planck:+.4f}"
    print(f"  {name:>20}  {val:12.4f}  {z:12.4f}  {tag:>15}")

# =============================================================================
# PART 9: SELF-CONSISTENCY: f_4 DETERMINES f_2?
# =============================================================================
print()
print("=" * 72)
print("PART 9: SELF-CONSISTENCY CHAIN")
print("=" * 72)
print()

# If f_4 = N^4/(2*pi) (exp408 PATTERN), the scalaron mass is:
# M_scal^2 = M_P^2 * (4*pi / (3*f_4)) * (c_2_F / c_0_F)
# In Connes' framework: M_scal = Lambda * sqrt(2*c_0 / (3*f_4 * pi))
# But this involves Lambda which is unknown.

# Alternative approach: f_4/f_2 ratio determines physics independently of Lambda.
# A_s = N_e^2 * f_2^2 / (72*pi^2 * f_4 * c_2_norm)
# If A_s is known (2.1e-9) and f_4 = N^4/(2*pi), can we get f_2?

As_obs = 2.1e-9
N_e = N / 2  # = 60 (framework: N/2 = |A5|)

# Connes normalization (simplified): A_s = N_e^2 / (12*pi^2 * f_4 * something_involving_f2)
# The exact relation depends on the full NCG model.
# Key ratio: f_4/f_2 = M_scal^2 * something / M_P^2
# M_scal ~ 2.9e13 GeV (from exp408)
# f_4/f_2 ~ (M_scal/M_P)^2 * numerical = (2.9e13/2.44e18)^2 ~ 1.4e-10

M_scal = 2.9e13  # GeV
M_P_red = 2.44e18  # reduced Planck mass

f4_over_f2_needed = 12 * As_obs * pi**2  # ~ 2.5e-7 (simplified Starobinsky)

f4 = N**4 / (2*pi)
f2_from_f4 = f4 / f4_over_f2_needed if f4_over_f2_needed > 0 else 0

print(f"  f_4 = N^4/(2*pi) = {f4:.4e}  [PATTERN, exp408]")
print(f"  f_4/f_2 needed ~ {f4_over_f2_needed:.4e}  (Starobinsky simplified)")
print(f"  => f_2 ~ {f2_from_f4:.4e}")
print()

# Is f2_from_f4 a framework number?
print(f"  Framework numbers near f_2 = {f2_from_f4:.4e}:")
f2_fw = {
    'N^8/(2pi)^2': N**8/(2*pi)**2,
    'N^6/(2pi)': N**6/(2*pi),
    'N^7/(2pi)^2': N**7/(2*pi)**2,
    'N^4*phi^4': N**4*PHI**4,
    'c0*f4': c0*f4,
}
for name, val in sorted(f2_fw.items(), key=lambda x: abs(log10(abs(x[1]))-log10(abs(f2_from_f4)))):
    ratio = val/f2_from_f4 if f2_from_f4 > 0 else 0
    print(f"    {name:>15} = {val:.4e}  (ratio: {ratio:.4f})")

# =============================================================================
# PART 10: SUMMARY
# =============================================================================
print()
print("=" * 72)
print("PART 10: SUMMARY")
print("=" * 72)
print()

print("WHAT IS DERIVED:")
print(f"  z_Planck = 4*phi^2 = 2*D^2 - d_ST = {z_Planck:.6f}  [exp413+exp462]")
print(f"  m_e/M_P = alpha^z_Planck  [DERIVED]")
print(f"  G_3D = 1/(4k) = 1/12  [DERIVED, exp424]")
print()

print("WHAT IS COMPUTED (new spectral data):")
print(f"  Full D^2 spectrum: 2640 eigenvalues on [{lam_min_full:.4f}, {lam_max_full:.4f}]")
print(f"  log det'(D^2) = {log_det_full:.4f}")
for k in range(4):
    print(f"  log det'(Delta_{k}) = {log_dets[k]:.4f}")
print()

print("WHAT REMAINS OPEN:")
print(f"  1. m_e (absolute value) = the ONE free parameter")
print(f"  2. f_2 (test function moment) - determines G via Connes formula")
print(f"  3. Lambda (cutoff) - determines the absolute scale")
print(f"  4. These are ALL equivalent: knowing any one determines the others")
print()
print("  STRUCTURAL CONCLUSION:")
print("  The framework determines ALL dimensionless ratios (masses, couplings,")
print("  hierarchy exponents). The absolute scale requires ONE dimensionful")
print("  input (m_e or M_P or Lambda). This is expected: any physical theory")
print("  needs at least one dimensionful parameter to set the scale.")
print()
print("  The test function problem (Connes 2008) is the NCG formulation of")
print("  this universal requirement. It is NOT specific to the 600-cell.")
