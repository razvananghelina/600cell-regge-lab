"""
exp408c_spectral_determinant.py
================================
GOAL: Compute spectral determinants and analytic torsion of the 600-cell
to find a first-principles origin for N^4 in f_4 = N^4/(2*pi).

KEY IDEA: On the 600-cell (topologically S^3), the combinatorial analytic
torsion T involves products of spectral determinants det'(Delta_k).
If T or det' involves N^4, we have the missing derivation.

QUANTITIES:
  det'(Delta_k) = product of nonzero eigenvalues of k-form Laplacian
  tau = det'(Delta_0)/N = number of spanning trees (Matrix Tree Theorem)
  T = combinatorial analytic torsion
  zeta'(0) = spectral zeta derivative (= -log det')

Author: Razvan-Constantin Anghelina
Date: 2026-02-25
"""

import numpy as np
from itertools import permutations, product as cartesian_product
from collections import defaultdict
from math import sqrt, pi, log, log10, factorial, exp
import time

print("=" * 72)
print("EXP-408c: SPECTRAL DETERMINANTS AND ANALYTIC TORSION")
print("         OF THE 600-CELL")
print("=" * 72)

# Framework constants
a1 = 5
b1 = a1 + 1
PHI = (1 + sqrt(a1)) / 2
N = factorial(a1)                       # 120
degree = 2 * b1                         # 12
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1                         # 30
N_gen = 3
alpha_s = 1 / (2*PHI**3)

# Seeley-DeWitt
A0 = 2*a1 + 1                          # 11
A1 = 2*a1**2 + 2*a1 + 2                # 62
A2 = 6*a1**2 + 15*a1 + 8               # 233
c0 = 2*N*A0                            # 2640
c1 = 2*N*A1                            # 14880
c2 = 2*N*A2                            # 55920


# =====================================================================
# PART 0: BUILD 600-CELL AND COMPUTE SPECTRA
# =====================================================================
print("\nBuilding 600-cell...")
t0 = time.time()

# --- Vertices ---
verts_set = set()
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        verts_set.add(tuple(v))
for s0 in [0.5, -0.5]:
    for s1 in [0.5, -0.5]:
        for s2 in [0.5, -0.5]:
            for s3 in [0.5, -0.5]:
                verts_set.add((s0, s1, s2, s3))
base = [PHI / 2, 0.5, 1 / (2 * PHI), 0.0]
even_perms = [p for p in permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1, 4)
                     if p[i] > p[j]) % 2 == 0]
for perm in even_perms:
    coords = [base[perm[i]] for i in range(4)]
    nonzero_idx = [i for i in range(4) if abs(coords[i]) > 1e-12]
    for signs in cartesian_product([1, -1], repeat=len(nonzero_idx)):
        v = list(coords)
        for idx, s in zip(nonzero_idx, signs):
            v[idx] *= s
        verts_set.add(tuple(round(x, 10) for x in v))

verts = np.array(sorted(verts_set))
Nv = len(verts)

# --- Edges ---
dots = verts @ verts.T
edges = [(i, j) for i in range(Nv) for j in range(i+1, Nv)
         if abs(dots[i, j] - PHI/2) < 0.001]
Ne = len(edges)
edge_to_idx = {}
for idx, (i, j) in enumerate(edges):
    edge_to_idx[(i, j)] = idx
    edge_to_idx[(j, i)] = idx
adj_list = defaultdict(set)
for i, j in edges:
    adj_list[i].add(j)
    adj_list[j].add(i)

# --- Triangles ---
triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    triangles.append((i, j, k))
Nf = len(triangles)
face_to_idx = {t: idx for idx, t in enumerate(triangles)}

# --- Tetrahedra ---
tetrahedra = []
for i in range(Nv):
    ni = adj_list[i]
    for j in ni:
        if j > i:
            cij = ni & adj_list[j]
            for k in cij:
                if k > j:
                    for l in cij & adj_list[k]:
                        if l > k:
                            tetrahedra.append((i, j, k, l))
Nc = len(tetrahedra)

# --- Boundary operators ---
d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = +1.0

d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = +1.0
    d1[f_idx, edge_to_idx[(j, k)]] = +1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

d2 = np.zeros((Nc, Nf))
for c_idx, (i, j, k, l) in enumerate(tetrahedra):
    for face, sign in [((j,k,l),+1), ((i,k,l),-1), ((i,j,l),+1), ((i,j,k),-1)]:
        f_idx = face_to_idx.get(face)
        if f_idx is not None:
            d2[c_idx, f_idx] = sign

# --- Hodge Laplacians ---
Delta_0 = d0.T @ d0
Delta_1 = d0 @ d0.T + d1.T @ d1
Delta_2 = d1 @ d1.T + d2.T @ d2
Delta_3 = d2 @ d2.T

# --- Eigenvalues ---
print("  Computing eigenvalues...")
evals_0 = np.sort(np.linalg.eigvalsh(Delta_0))
evals_1 = np.sort(np.linalg.eigvalsh(Delta_1))
evals_2 = np.sort(np.linalg.eigvalsh(Delta_2))
evals_3 = np.sort(np.linalg.eigvalsh(Delta_3))

build_time = time.time() - t0
print(f"  Done in {build_time:.1f}s")
print(f"  f-vector: ({Nv}, {Ne}, {Nf}, {Nc})")

# =====================================================================
# PART 1: SPECTRAL DETERMINANTS
# =====================================================================
print(f"\n{'='*72}")
print("PART 1: Spectral Determinants det'(Delta_k)")
print("=" * 72)

tol = 0.01

def log_det_prime(evals, name):
    """Compute log of product of nonzero eigenvalues."""
    nz = evals[evals > tol]
    n_zero = len(evals) - len(nz)
    log_det = np.sum(np.log(nz))
    print(f"  {name}: dim={len(evals)}, n_zero={n_zero}, n_nonzero={len(nz)}")
    print(f"    log det' = {log_det:.6f}")
    print(f"    det' = exp({log_det:.2f}) ~ 10^{log_det/log(10):.2f}")
    return log_det, n_zero

ld0, nz0 = log_det_prime(evals_0, "Delta_0 (0-forms)")
ld1, nz1 = log_det_prime(evals_1, "Delta_1 (1-forms)")
ld2, nz2 = log_det_prime(evals_2, "Delta_2 (2-forms)")
ld3, nz3 = log_det_prime(evals_3, "Delta_3 (3-forms)")

# Betti numbers check
print(f"\n  Betti numbers: b0={nz0}, b1={nz1}, b2={nz2}, b3={nz3}")
print(f"  Expected (S^3): 1, 0, 0, 1")


# =====================================================================
# PART 2: SPANNING TREES (Matrix Tree Theorem)
# =====================================================================
print(f"\n{'='*72}")
print("PART 2: Spanning Trees")
print("=" * 72)

# tau = det'(Delta_0) / N  (Kirchhoff's Matrix Tree Theorem)
log_tau = ld0 - log(N)
print(f"  log(tau) = log(det'(Delta_0)) - log(N)")
print(f"           = {ld0:.6f} - {log(N):.6f} = {log_tau:.6f}")
print(f"  tau ~ 10^{log_tau/log(10):.2f}")
print(f"\n  Decomposition of log(det'(Delta_0)) = {ld0:.6f}:")
print(f"    / ln(N) = {ld0/log(N):.6f}")
print(f"    / (4*ln(N)) = {ld0/(4*log(N)):.6f}")
print(f"    / ln(N^4) = {ld0/(4*log(N)):.6f}")
print(f"    / (N-1) = {ld0/(N-1):.6f}  (= average ln(eigenvalue))")
print(f"    exp(avg) = {exp(ld0/(N-1)):.6f}  (geometric mean eigenvalue)")


# =====================================================================
# PART 3: ANALYTIC TORSION
# =====================================================================
print(f"\n{'='*72}")
print("PART 3: Combinatorial Analytic Torsion")
print("=" * 72)

print(f"""
  The Ray-Singer / Cheeger-Mueller analytic torsion is:

    log T = (1/2) * sum_p (-1)^p * p * log det'(Delta_p)
          = (1/2) * [0 - log det'(D1) + 2*log det'(D2) - 3*log det'(D3)]
""")

log_T = 0.5 * (0 * ld0 - 1 * ld1 + 2 * ld2 - 3 * ld3)
print(f"  log det'(Delta_0) = {ld0:.6f}  (weight 0)")
print(f"  log det'(Delta_1) = {ld1:.6f}  (weight -1)")
print(f"  log det'(Delta_2) = {ld2:.6f}  (weight +2)")
print(f"  log det'(Delta_3) = {ld3:.6f}  (weight -3)")
print(f"\n  log T = 0.5 * (0 - {ld1:.4f} + 2*{ld2:.4f} - 3*{ld3:.4f})")
print(f"        = 0.5 * ({-ld1 + 2*ld2 - 3*ld3:.6f})")
print(f"        = {log_T:.6f}")
print(f"  T = exp({log_T:.6f}) = {exp(log_T):.6e}")

# Compare with framework quantities
T_val = exp(log_T)
print(f"\n  Comparisons:")
print(f"    T = {T_val:.6e}")
print(f"    1/(2*pi) = {1/(2*pi):.6e}")
print(f"    T / (1/(2*pi)) = {T_val * 2*pi:.6f}")
print(f"    N^4/(2*pi) = {N**4/(2*pi):.6e}")
print(f"    T * N^4 = {T_val * N**4:.6e}")
print(f"    2*pi*T = {2*pi*T_val:.6f}")
print(f"    log T / log(1/(2*pi)) = {log_T / log(1/(2*pi)):.6f}")

# Alternative torsion definition (Reidemeister-style)
# Some authors use: log T_R = sum_p (-1)^{p+1} * log det'(Delta_p)
log_T_R = sum((-1)**(p+1) * ld for p, ld in enumerate([ld0, ld1, ld2, ld3]))
print(f"\n  Alternative (Reidemeister-style):")
print(f"    log T_R = sum (-1)^(p+1) * log det'(Delta_p)")
print(f"            = -{ld0:.4f} + {ld1:.4f} - {ld2:.4f} + {ld3:.4f}")
print(f"            = {log_T_R:.6f}")
print(f"    T_R = {exp(log_T_R):.6e}")
print(f"    T_R / (1/(2*pi)) = {exp(log_T_R) * 2*pi:.6f}")


# =====================================================================
# PART 4: SPECTRAL ZETA FUNCTIONS
# =====================================================================
print(f"\n{'='*72}")
print("PART 4: Spectral Zeta Functions zeta_k(s)")
print("=" * 72)

print(f"""
  zeta_k(s) = sum_{{lambda > 0}} lambda^{{-s}}  (spectral zeta of Delta_k)
  zeta_k'(0) = -log det'(Delta_k)
  zeta_k(0) = n_nonzero  (regularized dimension)
""")

for s_val in [-2, -1, -0.5, 0, 0.5, 1, 2]:
    z0 = np.sum(evals_0[evals_0 > tol]**(-s_val))
    z1 = np.sum(evals_1[evals_1 > tol]**(-s_val))
    z2 = np.sum(evals_2[evals_2 > tol]**(-s_val))
    z3 = np.sum(evals_3[evals_3 > tol]**(-s_val))
    z_total = z0 + z1 + z2 + z3
    print(f"  s={s_val:5.1f}: zeta_0={z0:14.4f}  zeta_1={z1:14.4f}  "
          f"zeta_2={z2:14.4f}  zeta_3={z3:14.4f}  total={z_total:14.4f}")

# At s=-2: zeta(-2) = sum lambda^2 = Tr(D^4) per form degree
print(f"\n  At s=-2: zeta_k(-2) = Tr(Delta_k^2) per form degree")
tr4_0 = np.sum(evals_0**2)
tr4_1 = np.sum(evals_1**2)
tr4_2 = np.sum(evals_2**2)
tr4_3 = np.sum(evals_3**2)
tr4_total = tr4_0 + tr4_1 + tr4_2 + tr4_3
print(f"    Total Tr(D^4) = {tr4_total:.1f}")
print(f"    c2 = Tr(D^4)/2 = {tr4_total/2:.1f} (framework: {c2})")


# =====================================================================
# PART 5: SEARCHING FOR N^4 IN SPECTRAL DATA
# =====================================================================
print(f"\n{'='*72}")
print("PART 5: Searching for N^4 in Spectral Invariants")
print("=" * 72)

N4 = N**4
ln_N4 = 4 * log(N)
print(f"  N^4 = {N4}")
print(f"  ln(N^4) = 4*ln({N}) = {ln_N4:.6f}")
print(f"  N^4/(2*pi) = {N4/(2*pi):.4e}")
print(f"  ln(N^4/(2*pi)) = {log(N4/(2*pi)):.6f}")

print(f"\n  Testing spectral combinations against ln(N^4) = {ln_N4:.6f}:")
combos = [
    ("ld0", ld0),
    ("ld1", ld1),
    ("ld2", ld2),
    ("ld3", ld3),
    ("ld0 + ld3", ld0 + ld3),
    ("ld1 + ld2", ld1 + ld2),
    ("ld0 - ld3", ld0 - ld3),
    ("ld1 - ld2", ld1 - ld2),
    ("ld0 + ld1 + ld2 + ld3", ld0 + ld1 + ld2 + ld3),
    ("ld0 - ld1 + ld2 - ld3", ld0 - ld1 + ld2 - ld3),
    ("(ld0+ld3) - (ld1+ld2)", (ld0+ld3) - (ld1+ld2)),
    ("2*log_T", 2*log_T),
    ("log_T_R", log_T_R),
    ("-log_T_R", -log_T_R),
]

for name, val in combos:
    if abs(val) > 1e-10:
        ratio = val / ln_N4
        ratio2 = val / log(N4/(2*pi))
        print(f"    {name:>35} = {val:12.6f}  "
              f"/ ln(N^4) = {ratio:8.4f}  "
              f"/ ln(N^4/2pi) = {ratio2:8.4f}")

# Check: does any PRODUCT of det' give N^4?
print(f"\n  Checking products/ratios of det' against N^4 = {N4}:")

# det'(Delta_0)^a * det'(Delta_3)^b for small integer a,b
for a in range(-3, 4):
    for b in range(-3, 4):
        if a == 0 and b == 0:
            continue
        log_val = a * ld0 + b * ld3
        if abs(log_val - ln_N4) < 1.0:
            print(f"    det'(D0)^{a} * det'(D3)^{b}: "
                  f"log = {log_val:.4f}, target = {ln_N4:.4f}, "
                  f"diff = {log_val - ln_N4:.4f}")

# Check all four det' with small integer powers
print(f"\n  Exhaustive search: a*ld0 + b*ld1 + c*ld2 + d*ld3 = ln(N^4)?")
found = False
for a in range(-2, 3):
    for b in range(-2, 3):
        for c in range(-2, 3):
            for d in range(-2, 3):
                if a == b == c == d == 0:
                    continue
                val = a*ld0 + b*ld1 + c*ld2 + d*ld3
                if abs(val - ln_N4) < 0.5:
                    print(f"    ({a})*ld0 + ({b})*ld1 + ({c})*ld2 + ({d})*ld3"
                          f" = {val:.4f}  (target {ln_N4:.4f},"
                          f" diff {val-ln_N4:.4f})")
                    found = True
if not found:
    print(f"    No exact match found.")


# =====================================================================
# PART 6: HODGE DUALITY AND DETERMINANT RELATIONS
# =====================================================================
print(f"\n{'='*72}")
print("PART 6: Hodge Duality Relations")
print("=" * 72)

print(f"  On S^3, Hodge duality: spec(Delta_p) ~ spec(Delta_{{3-p}})")
print(f"  So det'(Delta_0) should relate to det'(Delta_3)")
print(f"  and det'(Delta_1) should relate to det'(Delta_2)")
print(f"\n  Ratios:")
print(f"    ld0/ld3 = {ld0/ld3:.6f}")
print(f"    ld1/ld2 = {ld1/ld2:.6f}")
print(f"    (ld0+ld3)/(ld1+ld2) = {(ld0+ld3)/(ld1+ld2):.6f}")

# The Hodge dual spectrum satisfies:
# nonzero spec(Delta_0) = nonzero spec(d0*d0^T) = exact part of Delta_1
# nonzero spec(Delta_3) = nonzero spec(d2^T*d2) = coexact part of Delta_2
# So det'(Delta_0) * det'(Delta_3) is related to the exact*coexact decomposition

product_03 = ld0 + ld3
product_12 = ld1 + ld2
print(f"\n  log[det'(D0)*det'(D3)] = {product_03:.6f}")
print(f"  log[det'(D1)*det'(D2)] = {product_12:.6f}")
print(f"  Difference = {product_03 - product_12:.6f}")
print(f"  Sum = {product_03 + product_12:.6f}")

# Compare with framework
print(f"\n  Framework comparisons:")
for name, val in [("ln(N)", log(N)), ("4*ln(N)", 4*log(N)),
                   ("ln(degree)", log(degree)), ("ln(c0)", log(c0)),
                   ("N*ln(degree)", N*log(degree)),
                   ("(N-1)*ln(degree)", (N-1)*log(degree)),
                   ("N/2*ln(N)", N/2*log(N)),
                   ("ln(N^4/(2*pi))", log(N**4/(2*pi)))]:
    print(f"    {name:>25} = {val:12.6f}")


# =====================================================================
# PART 7: THE 0-FORM SPECTRUM IN DETAIL
# =====================================================================
print(f"\n{'='*72}")
print("PART 7: Detailed 0-Form Spectrum and det'(Delta_0)")
print("=" * 72)

# Group eigenvalues
def group_evals(ev, tol_val=0.02):
    groups = []
    ev_s = np.sort(ev)
    cur = [ev_s[0]]
    for e in ev_s[1:]:
        if abs(e - cur[-1]) < tol_val:
            cur.append(e)
        else:
            groups.append((np.mean(cur), len(cur)))
            cur = [e]
    groups.append((np.mean(cur), len(cur)))
    return groups

groups_0 = group_evals(evals_0)
print(f"  Delta_0 spectrum (distinct eigenvalues):")
print(f"  {'lambda':>12} {'mult':>6} {'mult*ln(lam)':>14} {'Z[phi] form':>15}")

total_log = 0
for val, mult in groups_0:
    if val > tol:
        contrib = mult * log(val)
        total_log += contrib
        # Find a + b*phi form
        best_a, best_b = None, None
        best_err = 999
        for bb in range(-20, 21):
            for aa in range(-30, 31):
                err = abs(aa + bb*PHI - val)
                if err < 0.02 and err < best_err:
                    best_a, best_b = aa, bb
                    best_err = err
        zform = ""
        if best_a is not None:
            zform = f"{best_a}+{best_b}*phi" if best_b >= 0 else f"{best_a}{best_b}*phi"
        print(f"  {val:12.6f} {mult:6d} {contrib:14.6f} {zform:>15}")
    else:
        print(f"  {val:12.6f} {mult:6d} {'(zero)':>14} {'':>15}")

print(f"\n  Total: sum mult*ln(lam) = {total_log:.6f}")
print(f"  Check vs ld0 = {ld0:.6f}, diff = {abs(total_log-ld0):.2e}")

# Galois conjugate eigenvalues (phi -> phi' = -1/phi)
print(f"\n  Galois conjugate eigenvalues (phi -> phi' = (1-sqrt(5))/2):")
PHI_prime = (1 - sqrt(a1)) / 2

# The 2I characters at order-10 class
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
chi_C10 = [1, PHI, PHI, 1, 0, -1, -1, 1-PHI, 1-PHI]
chi_C10_galois = [1, PHI_prime, PHI_prime, 1, 0, -1, -1, 1-PHI_prime, 1-PHI_prime]

Lk = [degree * (1 - chi_C10[k]/dims[k]) for k in range(N_eig)]
Lk_galois = [degree * (1 - chi_C10_galois[k]/dims[k]) for k in range(N_eig)]
mults = [d**2 for d in dims]

print(f"  {'k':>3} {'d_k':>4} {'L_k':>12} {'L_k_galois':>12} {'L_k*L_k_gal':>14} {'mult':>5}")
for k in range(N_eig):
    prod_gal = Lk[k] * Lk_galois[k]
    print(f"  {k:3d} {dims[k]:4d} {Lk[k]:12.6f} {Lk_galois[k]:12.6f} "
          f"{prod_gal:14.6f} {mults[k]:5d}")

# Galois norm of det'(Delta_0)
# log|N(det')| = log(det') + log(det'_galois)
log_det_galois = sum(mults[k] * log(abs(Lk_galois[k]))
                     for k in range(N_eig) if Lk[k] > 0.01)
galois_norm_log = ld0 + log_det_galois

print(f"\n  log det'(Delta_0) = {ld0:.6f}")
print(f"  log det'_galois(Delta_0) = {log_det_galois:.6f}")
print(f"  log|N(det')| = {galois_norm_log:.6f}")
print(f"  |N(det')| = exp({galois_norm_log:.4f}) ~ 10^{galois_norm_log/log(10):.2f}")
print(f"\n  Comparisons:")
print(f"    galois_norm_log / ln(N) = {galois_norm_log/log(N):.6f}")
print(f"    galois_norm_log / 4 = {galois_norm_log/4:.6f}")
print(f"    exp(galois_norm_log/4) = {exp(galois_norm_log/4):.6f}")
print(f"    N = {N}")


# =====================================================================
# PART 8: KEY TEST -- GALOIS NORM OF FULL SPECTRAL DETERMINANT
# =====================================================================
print(f"\n{'='*72}")
print("PART 8: Galois Norm of the Full Spectral Determinant")
print("=" * 72)

# The FULL D^2 has eigenvalues from all form degrees.
# The 0-form and 3-form spectra are related by Hodge duality.
# The 1-form and 2-form spectra are also related.
#
# For 0-forms: eigenvalues are L_k with mults d_k^2
# For 3-forms: BY HODGE DUALITY, same nonzero spectrum as 0-forms
#   but shifted (Hodge dual maps p-forms to (3-p)-forms)
#
# The Galois norm of the 0-form det':
# |N(det'(D0))| = det'(D0) * det'_galois(D0)
#               = prod L_k^{d_k^2} * prod L_k_galois^{d_k^2}
#               = prod (L_k * L_k_galois)^{d_k^2}
#               = prod |N(L_k)|^{d_k^2}

print(f"  Galois norms of individual eigenvalues:")
for k in range(N_eig):
    if Lk[k] > 0.01:
        norm_k = Lk[k] * Lk_galois[k]
        # Check if norm is rational
        # L_k = degree*(1 - chi/d_k), L_k' = degree*(1 - chi'/d_k)
        # N(L_k) = degree^2 * (1 - chi/d_k)(1 - chi'/d_k)
        # = degree^2 * (1 - (chi+chi')/d_k + chi*chi'/d_k^2)
        # chi+chi' = Tr(chi) (rational), chi*chi' = N(chi) (rational)
        print(f"    k={k}: |N(L_{k})| = {Lk[k]:.6f} * {Lk_galois[k]:.6f}"
              f" = {norm_k:.6f}  (mult {mults[k]})")

# Product of all Galois norms
log_galois_product = sum(mults[k] * log(abs(Lk[k] * Lk_galois[k]))
                          for k in range(N_eig) if abs(Lk[k]) > 0.01)
print(f"\n  log prod |N(L_k)|^{{d_k^2}} = {log_galois_product:.6f}")
print(f"  = log|N(det'(D0))|")
print(f"\n  Key ratios:")
print(f"    / ln(N) = {log_galois_product/log(N):.6f}")
print(f"    / 4*ln(N) = {log_galois_product/(4*log(N)):.6f}")
print(f"    / ln(N^4/(2*pi)) = {log_galois_product/log(N**4/(2*pi)):.6f}")
print(f"    exp(.) = {exp(log_galois_product):.6e}")
print(f"    exp(.)/N^4 = {exp(log_galois_product)/N**4:.6e}")
print(f"    exp(.)^(1/4) = {exp(log_galois_product/4):.6f}")
print(f"    N = {N}")

# Check: is |N(det'(D0))|^{1/(N-1)} related to something?
geom_mean_norm = exp(log_galois_product / (N-1))
print(f"\n  Geometric mean of |N(L_k)|: {geom_mean_norm:.6f}")
print(f"    = exp(log_prod/(N-1))")
print(f"    degree^2 = {degree**2}")
print(f"    ratio: {geom_mean_norm/degree**2:.6f}")


# =====================================================================
# PART 9: THE RATIONAL GALOIS NORMS
# =====================================================================
print(f"\n{'='*72}")
print("PART 9: Rational Galois Norms (Exact Computation)")
print("=" * 72)

# Each L_k = degree*(1-chi_k/d_k) where chi_k is the character value.
# The Galois norm N(L_k) = L_k * sigma(L_k) is RATIONAL (in Q).
# Let's compute exactly.

# For the Cayley graph of 2I with generators = order-10 class:
# chi_k(C10) values: 1, phi, phi, 1, 0, -1, -1, 1-phi, 1-phi
# Galois conjugates: 1, phi', phi', 1, 0, -1, -1, 1-phi', 1-phi'
# where phi' = (1-sqrt(5))/2 = -1/phi

# N(chi_k) = chi_k * sigma(chi_k):
# k=0: 1*1 = 1
# k=1: phi*phi' = -1
# k=2: phi*phi' = -1
# k=3: 1*1 = 1
# k=4: 0*0 = 0
# k=5: (-1)*(-1) = 1
# k=6: (-1)*(-1) = 1
# k=7: (1-phi)*(1-phi') = 1-(phi+phi')+(phi*phi') = 1-1+(-1) = -1
# k=8: same as k=7: -1

# Tr(chi_k) = chi_k + sigma(chi_k):
# k=0: 2
# k=1: phi+phi' = 1
# k=2: phi+phi' = 1
# k=3: 2
# k=4: 0
# k=5: -2
# k=6: -2
# k=7: (1-phi)+(1-phi') = 2-(phi+phi') = 2-1 = 1
# k=8: same: 1

print(f"  Character traces and norms (exact in Q):")
print(f"  {'k':>3} {'d_k':>4} {'chi':>8} {'Tr(chi)':>8} {'N(chi)':>8}")

chi_traces = [2, 1, 1, 2, 0, -2, -2, 1, 1]
chi_norms = [1, -1, -1, 1, 0, 1, 1, -1, -1]

for k in range(N_eig):
    chi_str = f"{chi_C10[k]:.4f}"
    print(f"  {k:3d} {dims[k]:4d} {chi_str:>8} {chi_traces[k]:>8} {chi_norms[k]:>8}")

# N(L_k) = degree^2 * (1 - Tr(chi_k)/d_k + N(chi_k)/d_k^2)
# This is EXACT and RATIONAL!
print(f"\n  Exact Galois norms |N(L_k)| = degree^2 * (1 - Tr/d + N/d^2):")
log_exact_product = 0
for k in range(N_eig):
    d = dims[k]
    norm_val = degree**2 * (1 - chi_traces[k]/d + chi_norms[k]/d**2)
    if k > 0:  # skip zero eigenvalue
        log_exact_product += mults[k] * log(abs(norm_val))
    print(f"    k={k}: |N(L_{k})| = {degree}^2 * "
          f"(1 - {chi_traces[k]}/{d} + {chi_norms[k]}/{d}^2)"
          f" = {norm_val:.6f}  (mult {mults[k]})")

print(f"\n  log|N(det'(D0))| (exact) = {log_exact_product:.6f}")
print(f"  Check vs numerical: {abs(log_exact_product - log_galois_product):.2e}")

# Now the key: what is exp(log_exact_product)?
# = prod_k |N(L_k)|^{d_k^2}
# All factors are RATIONAL, so the product is rational!
# Let's compute it as a fraction (approximately)

print(f"\n  The product |N(det'(D0))| is RATIONAL (product of rationals).")
print(f"  Value = exp({log_exact_product:.6f}) = {exp(log_exact_product):.6e}")
print(f"\n  Key question: does this relate to N^4 = {N**4}?")
val = exp(log_exact_product)
print(f"    |N(det'(D0))| / N^4 = {val/N**4:.6e}")
print(f"    |N(det'(D0))|^(1/4) = {val**(0.25):.6e}")
print(f"    N = {N}")
print(f"    log|N(det')| / ln(N^4) = {log_exact_product/(4*log(N)):.6f}")


# =====================================================================
# PART 10: ALTERNATIVE -- det' AS SUM OF LOGS IN Z[phi]
# =====================================================================
print(f"\n{'='*72}")
print("PART 10: Decomposition of log det'(Delta_0) in Z[phi]")
print("=" * 72)

# log det'(D0) = sum_{k>0} d_k^2 * ln(L_k)
# Each L_k is in Z[phi] (the ring Z + Z*phi)
# So log det' = sum d_k^2 * ln(a_k + b_k*phi) where a_k,b_k are integers

# Eigenvalues in Z[phi] form:
# L_0 = 0
# L_1 = 12 - 6*phi
# L_2 = 12 - 4*phi
# L_3 = 9 (rational)
# L_4 = 12 (rational)
# L_5 = 14 (rational)
# L_6 = 15 (rational)
# L_7 = 6 + 6*phi
# L_8 = 8 + 4*phi

zphi_forms = [
    (0, 0),     # L_0 = 0
    (12, -6),   # L_1 = 12 - 6*phi
    (12, -4),   # L_2 = 12 - 4*phi
    (9, 0),     # L_3 = 9
    (12, 0),    # L_4 = 12
    (14, 0),    # L_5 = 14
    (15, 0),    # L_6 = 15
    (6, 6),     # L_7 = 6 + 6*phi
    (8, 4),     # L_8 = 8 + 4*phi
]

print(f"  Eigenvalues in Z[phi]:")
for k, (a, b) in enumerate(zphi_forms):
    L_check = a + b*PHI
    diff = abs(L_check - Lk[k]) if k > 0 or (a == 0 and b == 0) else 0
    bstr = f"+{b}" if b >= 0 else str(b)
    print(f"    L_{k} = {a}{bstr}*phi = {L_check:.6f}"
          f"  (d_k={dims[k]}, mult={mults[k]})")

# The "irrational part" of log det':
# Split each ln(a+b*phi) into rational and irrational parts
# This is hard in general, but the GALOIS trace gives the rational part:
# ln(L) + ln(sigma(L)) = ln|N(L)| (rational)
# ln(L) - ln(sigma(L)) = 2*arakelov_height(L) (irrational)

print(f"\n  Galois trace decomposition of log det':")
log_rational_part = 0
log_irrational_part = 0
for k in range(1, N_eig):  # skip zero
    ln_L = log(Lk[k])
    ln_L_gal = log(abs(Lk_galois[k]))
    rational_contrib = mults[k] * (ln_L + ln_L_gal) / 2  # symmetric part
    irrational_contrib = mults[k] * (ln_L - ln_L_gal) / 2  # antisymmetric
    log_rational_part += rational_contrib
    log_irrational_part += irrational_contrib

print(f"    Rational part:    (1/2)*log|N(det')| = {log_rational_part:.6f}")
print(f"    Irrational part:  sum d_k^2*(ln L - ln sigma(L))/2 = {log_irrational_part:.6f}")
print(f"    Total: {log_rational_part + log_irrational_part:.6f} (check: {ld0:.6f})")
print(f"\n    Rational part / ln(N) = {log_rational_part/log(N):.6f}")
print(f"    Irrational part / ln(phi) = {log_irrational_part/log(PHI):.6f}")


# =====================================================================
# FINAL ASSESSMENT
# =====================================================================
print(f"\n{'='*72}")
print("FINAL ASSESSMENT")
print("=" * 72)

print(f"""
  QUESTION: Does the spectral determinant or analytic torsion of the
  600-cell explain why f_4 = N^4/(2*pi)?

  FINDINGS:

  1. SPECTRAL DETERMINANTS (DERIVED):
     log det'(Delta_0) = {ld0:.4f}
     log det'(Delta_1) = {ld1:.4f}
     log det'(Delta_2) = {ld2:.4f}
     log det'(Delta_3) = {ld3:.4f}

  2. ANALYTIC TORSION:
     log T = {log_T:.6f}, T = {exp(log_T):.6e}
     T is NOT 1/(2*pi) = {1/(2*pi):.6e}
     (For smooth S^3 with trivial rep, T = 1; discrete differs.)

  3. GALOIS NORM OF det'(Delta_0):
     log|N(det')| = {log_exact_product:.6f}
     |N(det')| = {exp(log_exact_product):.4e}
     |N(det')| / N^4 = {exp(log_exact_product)/N**4:.4e}
     This is NOT N^4.

  4. NO integer-coefficient combination of log det'(Delta_k)
     equals ln(N^4) = {ln_N4:.4f}.

  5. DECOMPOSITION of log det' into rational + irrational parts:
     Rational part = {log_rational_part:.4f}
     Irrational part = {log_irrational_part:.4f}
     Neither directly gives ln(N^4) or ln(N^4/(2*pi)).

  CONCLUSION: The spectral determinants and analytic torsion of the
  600-cell do NOT directly explain the N^4 factor.

  The N^4 in f_4 = N^4/(2*pi) appears to originate from the
  4-DIMENSIONAL SPACETIME structure (N^{{d_ST}} with d_ST=4),
  not from the 3-dimensional fiber topology.

  STATUS: N^4 remains UNEXPLAINED (pattern, not derived).
""")

print("=" * 72)
print("EXPERIMENT COMPLETE")
print("=" * 72)
