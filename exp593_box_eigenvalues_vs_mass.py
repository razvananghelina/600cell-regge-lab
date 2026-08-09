#!/usr/bin/env python3
"""
EXP-593: BOX EIGENVALUES VS MASS EXPONENTS — Exista formula?
==============================================================

9 eigenvalues distincte Box <-> 9 fermioni pe 9 noduri McKay.

INTREBARE: Exista o relatie algebrica |mu_k| -> n_f?
Daca da -> masele fermionice vin DIRECT din spectrul Box.

Lantul: Box eigenvalue -> distanta McKay -> exponent masa
Verificam daca lantul e ALGEBRIC, nu doar monoton.
"""
import numpy as np
from numpy.linalg import eigh
import sys
sys.path.insert(0, ".")
from commons import build_600cell
from commons.constants import PHI, SQRT5, a1, b1, N

print("=" * 80)
print("EXP-593: BOX EIGENVALUES VS MASS EXPONENTS")
print("=" * 80)

# Build 600-cell + Hopf + Box
verts, adj, lap = build_600cell()

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

fibers = None
for i in range(N):
    if abs(verts[i, 0] - PHI/2) < 1e-6:
        g = verts[i]; p = g.copy(); ok = True
        for k in range(2, 11):
            p = qmul(p, g)
            if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
            if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
        if not ok: continue
        used = set(); fibers_tmp = []; subg = []
        pp = np.array([1.0,0,0,0])
        for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
        for s in range(N):
            if s in used: continue
            fib = []
            for si in subg:
                q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
            if len(fib) == 10: fibers_tmp.append(fib)
        if len(fibers_tmp) == 12:
            fibers = fibers_tmp
            break

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0

Box = b1 * A_fiber - adj

# Simultaneous diag to get eigenvalues per irrep
M_combo = lap + np.e * Box
evals_M, evecs_M = eigh(M_combo)
L_vals = np.array([evecs_M[:,i] @ lap @ evecs_M[:,i] for i in range(N)])
mu_vals = np.array([evecs_M[:,i] @ Box @ evecs_M[:,i] for i in range(N)])

# Known Laplacian eigenvalues per irrep
L_known = {
    0: 0.0, 1: 12-6*PHI, 2: 10-2*SQRT5, 3: 9.0,
    4: 12.0, 5: 14.0, 6: 10+2*SQRT5, 7: 15.0, 8: 6+6*PHI
}

irrep_dims = [1, 2, 3, 4, 5, 6, 3, 4, 2]

# =====================================================================
# PART 1: Box eigenvalues per irrep — detailed
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: BOX EIGENVALUES PER IRREP")
print("=" * 80)

box_data = {}
for rho in range(9):
    L_target = L_known[rho]
    d = irrep_dims[rho]
    indices = sorted(range(N), key=lambda i: abs(L_vals[i] - L_target))[:d**2]
    mus = sorted([mu_vals[i] for i in indices])

    # Distinct eigenvalues (within tolerance)
    distinct = []
    for m in mus:
        if not any(abs(m - d_val) < 0.01 for d_val in distinct):
            distinct.append(m)

    mean_abs = np.mean(np.abs(mus))
    box_data[rho] = {
        'dim': d, 'mult': d**2, 'mus': mus, 'distinct': distinct,
        'mean_abs': mean_abs, 'L': L_target
    }

print(f"\n  {'rho':>4s} {'dim':>4s} {'L':>8s} {'distinct Box evals':>45s} {'<|mu|>':>8s}")
print(f"  {'-'*75}")
for rho in range(9):
    d = box_data[rho]
    dist_str = ', '.join(f'{v:.4f}' for v in d['distinct'])
    print(f"  rho_{rho}  {d['dim']:>3d}  {d['L']:>8.4f}  {dist_str:>45s}  {d['mean_abs']:>8.4f}")

# =====================================================================
# PART 2: Fermion assignment — which node = which fermion
# =====================================================================
print("\n" + "=" * 80)
print("PART 2: FERMION -> NOD McKAY -> BOX EIGENVALUE")
print("=" * 80)

# From the paper's mass formula derivation (Theorem 13 + Theorem 5):
# Main chain walk with cumulative exponents:
# rho_0(n=0) -> rho_1(n=3) -> rho_2(n=5) -> rho_3(n=11) -> rho_4(n=11) -> rho_5(n=16)
# Branch: rho_6(n=17), rho_7(n=19), rho_8(n=26)

fermion_to_node = {
    'e':   0,   # electron,  n=0
    'u':   1,   # up,        n=3
    'd':   2,   # down,      n=5
    's':   3,   # strange,   n=11
    'mu':  4,   # muon,      n=11
    'c':   5,   # charm,     n=16
    'tau': 6,   # tau,       n=17
    'b':   7,   # bottom,    n=19
    't':   8,   # top,       n=26
}

fermion_n = {
    'e': 0, 'u': 3, 'd': 5, 's': 11, 'mu': 11,
    'c': 16, 'tau': 17, 'b': 19, 't': 26
}

print(f"\n  {'Fermion':>7s} {'rho':>5s} {'n_f':>4s} {'dim':>4s} {'<|mu|>':>8s} {'L':>8s}")
print(f"  {'-'*42}")

# Collect data for analysis
mu_values = []
n_values = []

for f in ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 'b', 't']:
    rho = fermion_to_node[f]
    n = fermion_n[f]
    d = box_data[rho]
    print(f"  {f:>7s}  rho_{rho}  {n:>3d}  {d['dim']:>3d}  {d['mean_abs']:>8.4f}  {d['L']:>8.4f}")
    mu_values.append(d['mean_abs'])
    n_values.append(n)

mu_arr = np.array(mu_values)
n_arr = np.array(n_values)

# =====================================================================
# PART 3: Correlation analysis
# =====================================================================
print("\n" + "=" * 80)
print("PART 3: CORELATIE |mu| vs n_f")
print("=" * 80)

# Pearson correlation
corr = np.corrcoef(mu_arr, n_arr)[0, 1]
print(f"\n  Pearson correlation: r = {corr:.6f}")

# Spearman (rank) correlation
from scipy.stats import spearmanr
rho_s, p_s = spearmanr(mu_arr, n_arr)
print(f"  Spearman rank correlation: rho = {rho_s:.6f}, p = {p_s:.6e}")

# Monotonic?
sorted_by_mu = sorted(zip(mu_arr, n_arr))
print(f"\n  Sorted by |mu|:")
for mu, n in sorted_by_mu:
    print(f"    |mu| = {mu:.4f}  ->  n = {n}")

monotonic = all(sorted_by_mu[i][1] <= sorted_by_mu[i+1][1]
                for i in range(len(sorted_by_mu)-1))
print(f"\n  Monotonic? {monotonic}")

# =====================================================================
# PART 4: Search for algebraic relation
# =====================================================================
print("\n" + "=" * 80)
print("PART 4: CAUTARE RELATIE ALGEBRICA |mu| -> n_f")
print("=" * 80)

# Try: n = a * |mu| + b (linear)
# Least squares
A_matrix = np.column_stack([mu_arr, np.ones(9)])
coeffs_lin, residuals_lin, _, _ = np.linalg.lstsq(A_matrix, n_arr, rcond=None)
n_pred_lin = A_matrix @ coeffs_lin
err_lin = np.sqrt(np.mean((n_pred_lin - n_arr)**2))
print(f"\n  Linear: n = {coeffs_lin[0]:.4f}*|mu| + {coeffs_lin[1]:.4f}")
print(f"  RMSE = {err_lin:.4f}")
for i, f in enumerate(['e','u','d','s','mu','c','tau','b','t']):
    print(f"    {f:>3s}: pred={n_pred_lin[i]:.1f}, actual={n_arr[i]}, err={n_pred_lin[i]-n_arr[i]:+.1f}")

# Try: n = a * |mu|^2 + b * |mu| + c (quadratic)
A_quad = np.column_stack([mu_arr**2, mu_arr, np.ones(9)])
coeffs_quad, _, _, _ = np.linalg.lstsq(A_quad, n_arr, rcond=None)
n_pred_quad = A_quad @ coeffs_quad
err_quad = np.sqrt(np.mean((n_pred_quad - n_arr)**2))
print(f"\n  Quadratic: n = {coeffs_quad[0]:.4f}*|mu|^2 + {coeffs_quad[1]:.4f}*|mu| + {coeffs_quad[2]:.4f}")
print(f"  RMSE = {err_quad:.4f}")

# Try: n = a * L + b (using Laplacian eigenvalue instead)
L_arr = np.array([box_data[fermion_to_node[f]]['L'] for f in
                   ['e','u','d','s','mu','c','tau','b','t']])
A_L = np.column_stack([L_arr, np.ones(9)])
coeffs_L, _, _, _ = np.linalg.lstsq(A_L, n_arr, rcond=None)
n_pred_L = A_L @ coeffs_L
err_L = np.sqrt(np.mean((n_pred_L - n_arr)**2))
print(f"\n  Linear in L: n = {coeffs_L[0]:.4f}*L + {coeffs_L[1]:.4f}")
print(f"  RMSE = {err_L:.4f}")

# Try: n = a * L^2 + b * L + c
A_L2 = np.column_stack([L_arr**2, L_arr, np.ones(9)])
coeffs_L2, _, _, _ = np.linalg.lstsq(A_L2, n_arr, rcond=None)
n_pred_L2 = A_L2 @ coeffs_L2
err_L2 = np.sqrt(np.mean((n_pred_L2 - n_arr)**2))
print(f"\n  Quadratic in L: n = {coeffs_L2[0]:.6f}*L^2 + {coeffs_L2[1]:.4f}*L + {coeffs_L2[2]:.4f}")
print(f"  RMSE = {err_L2:.4f}")

# Try: n = a * dim * something
dim_arr = np.array([irrep_dims[fermion_to_node[f]] for f in
                     ['e','u','d','s','mu','c','tau','b','t']])

# Try: combinations of dim and L
A_dL = np.column_stack([dim_arr * L_arr, dim_arr, L_arr, np.ones(9)])
coeffs_dL, _, _, _ = np.linalg.lstsq(A_dL, n_arr, rcond=None)
n_pred_dL = A_dL @ coeffs_dL
err_dL = np.sqrt(np.mean((n_pred_dL - n_arr)**2))
print(f"\n  Combined dim*L: n = {coeffs_dL[0]:.4f}*d*L + {coeffs_dL[1]:.4f}*d + {coeffs_dL[2]:.4f}*L + {coeffs_dL[3]:.4f}")
print(f"  RMSE = {err_dL:.4f}")
if err_dL < 0.5:
    print(f"  EXCELLENT FIT!")
    for i, f in enumerate(['e','u','d','s','mu','c','tau','b','t']):
        print(f"    {f:>3s}: pred={n_pred_dL[i]:.2f}, actual={n_arr[i]}")

# =====================================================================
# PART 5: Try EXACT relations using Z[phi] structure
# =====================================================================
print("\n" + "=" * 80)
print("PART 5: RELATII EXACTE IN Z[phi]")
print("=" * 80)

# The Box eigenvalues per irrep should be expressible in Z[phi]
# Let's compute them exactly

# Box = b1*A_fiber - A_full
# On irrep rho_k with adjacency eigenvalue lambda_k:
# Box has eigenvalue mu_k = b1*lambda_fiber_k - lambda_k
# where lambda_fiber_k is the adjacency eigenvalue of A_fiber on irrep rho_k

# The fiber adjacency A_fiber has its own eigenvalues per irrep
# We can compute these: project A_fiber to each irrep

# Actually, easier: Box = b1*A_fiber - A_full
# lambda_Box = b1*lambda_fiber - lambda_full
# But lambda_fiber and lambda_full are different for each irrep

# Let me compute the mean Box eigenvalue per irrep more carefully
# Actually, each irrep has d^2 Box eigenvalues, not just one

# For the MEAN: <mu>_rho = Tr(Box|_rho) / d^2
# But Tr(Box) = 0 per irrep? Let me check

print(f"\n  Box eigenvalue structure per irrep:")
for rho in range(9):
    d = box_data[rho]
    mus = d['mus']
    tr = sum(mus)
    print(f"  rho_{rho} (d={d['dim']}): Tr={tr:.4f}, distinct={[f'{v:.4f}' for v in d['distinct']]}")

# Key observation: Box has distinct eigenvalues per irrep
# Some irreps have only ONE distinct value (up to sign)
# rho_0: {0}
# rho_1: {0}
# rho_2: {-2.764, 5.528} = {-L2/2, L2} ≈ {-(a1-sqrt5), 2(a1-sqrt5)}
# rho_3: {-6.708, 6.708} = {-3*sqrt5, 3*sqrt5} ≈ ±sqrt(45) = ±3sqrt5
# rho_4: multiple values
# etc.

# Let me identify exact forms
print(f"\n  Identificare forme exacte:")
for rho in range(9):
    d = box_data[rho]
    for v in d['distinct']:
        # Try common forms
        candidates = []
        for expr_name, expr_val in [
            ("0", 0), ("a1", a1), ("-a1", -a1),
            ("sqrt5", SQRT5), ("-sqrt5", -SQRT5),
            ("3*sqrt5", 3*SQRT5), ("-3*sqrt5", -3*SQRT5),
            ("phi", PHI), ("-phi", -PHI),
            ("2*phi", 2*PHI), ("-2*phi", -2*PHI),
            ("b1*phi", b1*PHI), ("-b1*phi", -b1*PHI),
            ("b1/phi", b1/PHI), ("-b1/phi", -b1/PHI),
            ("a1-sqrt5", a1-SQRT5), ("a1+sqrt5", a1+SQRT5),
            ("-(a1-sqrt5)", -(a1-SQRT5)), ("-(a1+sqrt5)", -(a1+SQRT5)),
            ("2*(a1-sqrt5)", 2*(a1-SQRT5)), ("2*(a1+sqrt5)", 2*(a1+SQRT5)),
            ("L2", 10-2*SQRT5), ("-L2", -(10-2*SQRT5)),
            ("L6", 10+2*SQRT5), ("-L6", -(10+2*SQRT5)),
            ("L2/2", (10-2*SQRT5)/2), ("-L2/2", -(10-2*SQRT5)/2),
            ("L6/2", (10+2*SQRT5)/2), ("-L6/2", -(10+2*SQRT5)/2),
            ("b1-2*phi", b1-2*PHI), ("12-b1*phi", 12-b1*PHI),
            ("b1*phi-12", b1*PHI-12), ("2*b1-12", 2*b1-12),
            ("12", 12), ("-12", -12), ("10", 10), ("-10", -10),
            ("6*phi", 6*PHI), ("-6*phi", -6*PHI),
            ("6/phi", 6/PHI), ("-6/phi", -6/PHI),
            ("b1*(phi-2)", b1*(PHI-2)), ("b1*(2-phi)", b1*(2-PHI)),
        ]:
            if abs(v - expr_val) < 0.001:
                candidates.append(expr_name)
        if candidates:
            print(f"    rho_{rho}: mu = {v:>8.4f} = {candidates[0]}")
        else:
            print(f"    rho_{rho}: mu = {v:>8.4f} = ???")

# =====================================================================
# PART 6: The key test — |mu| vs n sorted
# =====================================================================
print("\n" + "=" * 80)
print("PART 6: |mu| SORTAT VS n_f SORTAT — EXISTA FORMULA?")
print("=" * 80)

# Sort fermions by <|mu|>
sorted_fermions = sorted(
    [(f, fermion_to_node[f], fermion_n[f], box_data[fermion_to_node[f]]['mean_abs'])
     for f in fermion_n],
    key=lambda x: x[3]
)

print(f"\n  {'Fermion':>7s} {'rho':>5s} {'n_f':>4s} {'<|mu|>':>8s} {'n/|mu|':>8s} {'n*|mu|':>8s}")
print(f"  {'-'*48}")
for f, rho, n, mu in sorted_fermions:
    ratio = n/mu if mu > 0.01 else 0
    prod = n * mu
    print(f"  {f:>7s}  rho_{rho}  {n:>3d}  {mu:>8.4f}  {ratio:>8.4f}  {prod:>8.4f}")

# The 3 kernel fermions (mu=0) have n=0,3,26... NOT all zero!
# So mu=0 does NOT mean n=0.
print(f"""
  PROBLEMA CRITICA:
  Kernel-ul Box (mu=0) contine rho_0(e,n=0), rho_1(u,n=3), rho_8(t,n=26).
  Deci |mu|=0 corespunde la n=0, n=3, SI n=26!

  NU exista functie |mu| -> n (nu e injectiva pe kernel).
  Trei fermioni cu masa Box ZERO au exponenti de masa DIFERITI.

  Concluzie: |mu| SINGUR nu determina n_f.
  Trebuie si informatia despre CARE irrep (dim, Galois type, etc.)
""")

# =====================================================================
# PART 7: Use (|mu|, dim) or (|mu|, Galois) as input
# =====================================================================
print("=" * 80)
print("PART 7: (|mu|, dim) SAU (|mu|, Galois type) -> n_f?")
print("=" * 80)

# Since |mu| alone fails, try |mu| + dim
print(f"\n  {'Fermion':>7s} {'rho':>5s} {'n_f':>4s} {'dim':>4s} {'<|mu|>':>8s} {'Galois':>8s}")
print(f"  {'-'*48}")
for f, rho, n, mu in sorted_fermions:
    gtype = 'fixed' if rho in [0,3,4,5] else 'broken'  # from earlier analysis
    print(f"  {f:>7s}  rho_{rho}  {n:>3d}  {irrep_dims[rho]:>3d}  {mu:>8.4f}  {gtype:>8s}")

# With (dim, |mu|, Galois type) as input, is the map injective?
# e:   (1, 0,   fixed)  -> 0
# u:   (2, 0,   broken) -> 3
# d:   (3, 3.7, broken) -> 5
# s:   (4, 6.7, broken) -> 11  (actually fixed for rho_3?)
# mu:  (5, 7.8, fixed)  -> 11
# c:   (6, 7.8, fixed)  -> 16
# tau: (3, 9.6, broken) -> 17
# b:   (4, 6.7, broken) -> 19
# t:   (2, 0,   broken) -> 26

# Still not injective: rho_3 and rho_7 both have dim=4, mu=6.7, broken
# But s(rho_3)=11 and b(rho_7)=19

# =====================================================================
# PART 8: VERDICT
# =====================================================================
print(f"\n" + "=" * 80)
print("V E R D I C T")
print("=" * 80)

print(f"""
  REZULTAT: |mu_k| (Box eigenvalue) NU determina unic n_f (mass exponent).

  PROBLEME:
  1. Kernel-ul Box (|mu|=0) contine 3 fermioni cu n=0, 3, 26
     -> functia nu e injectiva
  2. Galois pairs (rho_3, rho_7) au ACELASI |mu| dar n diferit (11 vs 19)
     -> nici (dim, |mu|) nu e suficient

  CE FUNCTIONEAZA:
  - Corelatia Pearson r = {corr:.3f} (buna dar nu perfecta)
  - Spearman rho = {rho_s:.3f} (rang, nu valor)
  - Tendinta monotona generala (cu exceptii pe kernel)

  CE NU FUNCTIONEAZA:
  - Formula algebrica |mu| -> n
  - Functie injectiva din spectrul Box in exponentii de masa

  CONCLUZIE:
  Spectrul Box CORELEAZA cu masele fermionice (structura generala),
  dar NU le DETERMINA unic. Informatia suplimentara (din McKay walk,
  Fibonacci, Galois cascade) e NECESARA.

  Derivarea maselor ramane pe lantul multi-step verificat in exp585:
    Fibonacci -> Generation Theorem -> Casimir offsets -> Sum rule

  Box contribuie INDIRECT (prin structura spectrala care defineste
  alpha, c^2, gravitatia) dar nu da masele fermionice DIRECT.
""")
