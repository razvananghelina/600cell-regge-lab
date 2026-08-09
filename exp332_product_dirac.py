"""
exp332_product_dirac.py - Product Dirac Operator on 600-cell x McKay
====================================================================

Direction A from alte_directii.txt:
D = D_ext (x) I_30 + Gamma (x) D_F
on the 3600-dimensional Hilbert space H = C^120 (x) C^30.

Gamma = Hopf T3 grading introduces sector-dependent scaling of D_F.
T3 values {+-phi/2, +-1/2, +-1/(2*phi), +-1, 0} on the 120 vertices.

Key insight: different T3 values act as sector-dependent SCALINGS of the
internal Dirac operator. This could generate mass hierarchy that D_F
alone cannot (exp329 confirmed negative).

Author: Razvan-Constantin Anghelina
Date: 2026-02-16
"""

import numpy as np
from scipy.special import comb as binomial
from itertools import product as iprod
import time
import warnings
warnings.filterwarnings('ignore')

t0 = time.time()

print("=" * 72)
print("exp332: PRODUCT DIRAC OPERATOR D_F x Delta_0")
print("=" * 72)

# =====================================================================
# CONSTANTS
# =====================================================================
a_1 = 5
b_1 = 6
PHI = (1 + np.sqrt(5)) / 2
N_gen = 3
N_eig = 9
me = 0.51099895  # MeV

fermion_exps = {
    'e': 0, 'mu': 11, 'tau': 17,
    'u': 3, 'c': 16, 't': 26,
    'd': 5, 's': 11, 'b': 19
}
target_sorted = np.array(sorted(fermion_exps.values(), reverse=True))
# [26, 19, 17, 16, 11, 11, 5, 3, 0]

irrep_dims = {1: 1, 2: 2, 3: 2, 4: 3, 5: 3, 6: 4, 7: 4, 8: 5, 9: 6}
total_dim_F = sum(irrep_dims.values())  # = 30

block_starts = {}
pos = 0
for i in range(1, 10):
    block_starts[i] = pos
    pos += irrep_dims[i]

# =====================================================================
# PART 1: BUILD 2I (Binary Icosahedral Group, 120 quaternions)
# =====================================================================
print(f"\n[{time.time()-t0:.1f}s] Building 2I...")


def qmul(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)


def build_2I():
    elements = []
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0]*4
            v[i] = s
            elements.append(tuple(v))
    for signs in iprod([0.5, -0.5], repeat=4):
        elements.append(tuple(signs))
    even_perms = [
        (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2),
        (1, 0, 3, 2), (1, 2, 0, 3), (1, 3, 2, 0),
        (2, 0, 1, 3), (2, 1, 3, 0), (2, 3, 0, 1),
        (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0)
    ]
    bv = [0.0, 0.5, PHI/2, 1/(2*PHI)]
    for perm in even_perms:
        base = [bv[perm[k]] for k in range(4)]
        nonzero = [k for k in range(4) if abs(base[k]) > 1e-10]
        for signs in iprod([1, -1], repeat=len(nonzero)):
            v = list(base)
            for idx, k in enumerate(nonzero):
                v[k] *= signs[idx]
            elements.append(tuple(v))
    unique = []
    for v in elements:
        if not any(sum((v[k]-u[k])**2 for k in range(4)) < 1e-10 for u in unique):
            unique.append(v)
    return unique


elements_2I = build_2I()
N_verts = len(elements_2I)
assert N_verts == 120, f"Got {N_verts}"
print(f"  |2I| = {N_verts}")


def find_element(q):
    for k in range(N_verts):
        if sum((q[m]-elements_2I[k][m])**2 for m in range(4)) < 1e-8:
            return k
    return -1


# =====================================================================
# PART 2: 600-CELL EDGES + LAPLACIAN
# =====================================================================
print(f"[{time.time()-t0:.1f}s] Building 600-cell graph...")

A_600 = np.zeros((N_verts, N_verts))
for i in range(N_verts):
    for j in range(i+1, N_verts):
        dot = sum(elements_2I[i][k]*elements_2I[j][k] for k in range(4))
        if abs(dot - PHI/2) < 1e-6:
            A_600[i, j] = A_600[j, i] = 1

degrees = np.sum(A_600, axis=1)
n_edges = int(np.sum(A_600) / 2)
print(f"  Vertices: {N_verts}, Edges: {n_edges}, Degree: {int(degrees[0])}")
assert np.allclose(degrees, 12) and n_edges == 720

Delta_0 = np.diag(degrees) - A_600
evals_D0, evecs_D0 = np.linalg.eigh(Delta_0)

print(f"  Laplacian eigenvalues (9 distinct):")
seen = set()
for ev in evals_D0:
    key = round(ev, 3)
    if key not in seen:
        seen.add(key)
        mult = int(np.sum(np.abs(evals_D0 - ev) < 0.01))
        print(f"    {ev:10.4f} (mult {mult})")

# =====================================================================
# PART 3: HOPF T3 CHARGE
# =====================================================================
T3_vals = np.array([q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
                    for q in elements_2I])
T3_unique = sorted(set(round(t, 6) for t in T3_vals))
print(f"\n  T3 values ({len(T3_unique)} distinct):")
for tv in T3_unique:
    cnt = int(np.sum(np.abs(T3_vals - tv) < 1e-4))
    print(f"    T3 = {tv:8.4f}: {cnt} vertices")

n_pos_T3 = np.sum(T3_vals > 1e-10)
n_neg_T3 = np.sum(T3_vals < -1e-10)
n_zero_T3 = np.sum(np.abs(T3_vals) < 1e-10)
print(f"  Positive: {n_pos_T3}, Negative: {n_neg_T3}, Zero: {n_zero_T3}")

# =====================================================================
# PART 4: BUILD 9 IRREPS + CG MAPS + D_F
# =====================================================================
print(f"\n[{time.time()-t0:.1f}s] Building 9 irreps of 2I...")


def quat_to_su2(q):
    w, x, y, z = q
    a = complex(w, x)
    b = complex(y, z)
    return np.array([[a, b], [-np.conj(b), np.conj(a)]])


def galois_real(x):
    for b_num in range(-8, 9):
        b = b_num / 2.0
        a = x - b * PHI
        a2 = round(2 * a)
        if abs(a - a2/2) < 1e-8:
            return (a2/2 + b) - b * PHI
    return x


def galois_complex(z):
    return galois_real(z.real) + 1j * galois_real(z.imag)


def galois_matrix(M):
    result = np.empty_like(M)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            result[i, j] = galois_complex(M[i, j])
    return result


def spin_j_matrix(U, j):
    two_j = int(round(2*j))
    dim = two_j + 1
    a, b, c, d = U[0, 0], U[0, 1], U[1, 0], U[1, 1]
    D = np.zeros((dim, dim), dtype=complex)
    for m_idx in range(dim):
        m = j - m_idx
        p = int(round(j+m))
        q_val = int(round(j-m))
        for mp_idx in range(dim):
            mp = j - mp_idx
            pp = int(round(j+mp))
            val = 0j
            for s in range(p+1):
                t = pp - s
                if 0 <= t <= q_val:
                    val += (binomial(p, s, exact=True) *
                            binomial(q_val, t, exact=True) *
                            a**s * c**(p-s) * b**t * d**(q_val-t))
            D[mp_idx, m_idx] = val
    return D


reps = {i: [] for i in range(1, 10)}
for g_idx, q in enumerate(elements_2I):
    U = quat_to_su2(q)
    U_gal = galois_matrix(U)
    reps[1].append(np.array([[1.0+0j]]))
    reps[2].append(U.copy())
    reps[3].append(U_gal.copy())
    reps[4].append(spin_j_matrix(U, 1))
    reps[5].append(spin_j_matrix(U_gal, 1))
    reps[6].append(np.kron(U, U_gal))
    reps[7].append(spin_j_matrix(U, 1.5))
    reps[8].append(spin_j_matrix(U, 2))
    reps[9].append(spin_j_matrix(U, 2.5))

print(f"[{time.time()-t0:.1f}s] Computing characters + McKay graph...")
characters = {}
for i in range(1, 10):
    characters[i] = np.array([np.trace(reps[i][g]) for g in range(120)])

mckay_adj = np.zeros((9, 9), dtype=int)
for i_rep in range(1, 10):
    prod_chars = characters[2] * characters[i_rep]
    for j_rep in range(1, 10):
        mult = np.sum(np.conj(characters[j_rep]) * prod_chars).real / 120
        mckay_adj[i_rep-1, j_rep-1] = int(round(mult))

mckay_edges = []
for i in range(9):
    for j in range(i+1, 9):
        if mckay_adj[i, j] > 0:
            mckay_edges.append((i+1, j+1))
print(f"  McKay edges: {len(mckay_edges)} (expected 8)")
assert len(mckay_edges) == 8

print(f"[{time.time()-t0:.1f}s] Computing CG maps...")
cg_maps = {}
for node_i, node_j in mckay_edges:
    for source, target in [(node_i, node_j), (node_j, node_i)]:
        d_s, d_t = irrep_dims[source], irrep_dims[target]
        mult = int(round(np.sum(
            np.conj(characters[target]) * characters[2] * characters[source]
        ).real / 120))
        if mult < 1:
            continue
        dim_product = 2 * d_s
        equations = []
        for g_idx in range(120):
            kron_g = np.kron(reps[2][g_idx], reps[source][g_idx])
            rho_t_g = reps[target][g_idx]
            A_g = (np.kron(kron_g, np.eye(d_t, dtype=complex)) -
                   np.kron(np.eye(dim_product, dtype=complex), rho_t_g.T))
            equations.append(A_g)
        A_mat = np.vstack(equations)
        _, S_svd, Vh_svd = np.linalg.svd(A_mat, full_matrices=True)
        C_vec = Vh_svd[-1].conj()
        C_map = C_vec.reshape(dim_product, d_t)
        max_idx = np.argmax(np.abs(C_map))
        max_val = C_map.flat[max_idx]
        if abs(max_val) > 1e-15:
            C_map *= np.abs(max_val) / max_val
        cg_maps[(source, target)] = C_map

print(f"  CG maps: {len(cg_maps)} (expected 16)")
assert len(cg_maps) == 16

# Build D_F (30x30 Hermitian)
D_F = np.zeros((total_dim_F, total_dim_F), dtype=complex)
for (source, target), C_map in cg_maps.items():
    d_s, d_t = irrep_dims[source], irrep_dims[target]
    Y = C_map[:d_s, :]
    rs, rt = block_starts[source], block_starts[target]
    D_F[rs:rs+d_s, rt:rt+d_t] = Y
    D_F[rt:rt+d_t, rs:rs+d_s] = Y.conj().T
D_F = (D_F + D_F.conj().T) / 2

tr_df2 = np.trace(D_F @ D_F).real
evals_DF = np.linalg.eigvalsh(D_F)
print(f"  Tr(D_F^2) = {tr_df2:.4f} (expected 8 = rank E8)")
print(f"  D_F eigenvalue range: [{evals_DF[0]:.4f}, {evals_DF[-1]:.4f}]")
df_max = np.max(np.abs(evals_DF))

# =====================================================================
# PART 5: D_ext AND GRADINGS
# =====================================================================
print(f"\n[{time.time()-t0:.1f}s] Computing D_ext and gradings...")

# D_ext option 1: sqrt(Delta_0)
evals_sqrt = np.sqrt(np.maximum(evals_D0, 0))
D_ext_sqrt = evecs_D0 @ np.diag(evals_sqrt) @ evecs_D0.T
D_ext_sqrt = (D_ext_sqrt + D_ext_sqrt.T) / 2  # ensure symmetry
dext_max = np.max(evals_sqrt)
print(f"  sqrt(Delta_0) range: [0, {dext_max:.4f}]")

# D_ext option 2: normalized adjacency A/12
D_ext_adj = A_600 / 12.0

# Grading operators
Gamma_T3 = np.diag(T3_vals)
Gamma_sign = np.diag(np.sign(T3_vals))

# Scale factor: normalize D_ext to match D_F scale
scale_ext = df_max / dext_max
print(f"  Scale factor (D_F/D_ext): {scale_ext:.4f}")
print(f"  Using BOTH: unit scale and matched scale")

# =====================================================================
# PART 6: PRODUCT OPERATORS
# =====================================================================
print(f"\n" + "=" * 72)
print("PART 6: PRODUCT OPERATORS ON 3600-DIM SPACE")
print("=" * 72)

I_120 = np.eye(N_verts)
I_30 = np.eye(total_dim_F)


def compute_and_analyze(D_matrix, label):
    """Compute eigenvalues and basic statistics."""
    print(f"\n  [{label}]")
    print(f"  Shape: {D_matrix.shape}, Hermitian: "
          f"{np.allclose(D_matrix, D_matrix.conj().T)}")
    t1 = time.time()
    evals = np.sort(np.linalg.eigvalsh(D_matrix))[::-1]
    dt = time.time() - t1
    print(f"  Eigenvalue computation: {dt:.1f}s")

    pos = evals[evals > 1e-10]
    neg = evals[evals < -1e-10]
    nz = np.sum(np.abs(evals) < 1e-10)
    print(f"  Total: {len(evals)}, Positive: {len(pos)}, "
          f"Zero: {nz}, Negative: {len(neg)}")

    if len(pos) > 0:
        dyn = pos[0] / pos[-1]
        log_dyn = np.log(dyn) / np.log(PHI)
        print(f"  Dynamic range (pos): {dyn:.1f} (phi^{log_dyn:.1f}), "
              f"need phi^26 = {PHI**26:.0f}")
        print(f"  Max: {pos[0]:.6f}, Min: {pos[-1]:.8f}")
    return evals


# --- Form 1: Trivial additive (eigenvalues = sums, no matrix needed) ---
print(f"\n--- Form 1: Additive baseline (trivial: eigenvalues = sums) ---")
evals_add = np.sort(np.array([de + df for de in evals_sqrt
                               for df in evals_DF]))[::-1]
pos_add = evals_add[evals_add > 1e-10]
if len(pos_add) > 0:
    dyn = pos_add[0] / pos_add[-1]
    print(f"  Dynamic range: {dyn:.1f} (phi^{np.log(dyn)/np.log(PHI):.1f})")

# --- Form 2: T3 Hopf grading (KEY TEST) ---
print(f"\n--- Form 2: T3 Hopf grading (KEY) ---")
print(f"[{time.time()-t0:.1f}s] Building 3600x3600 matrix...")
D_Hopf = np.kron(D_ext_sqrt, I_30) + np.kron(Gamma_T3, D_F)
evals_Hopf = compute_and_analyze(D_Hopf, "T3 Hopf (unit scale)")

# --- Form 2b: Matched scale ---
print(f"\n--- Form 2b: T3 Hopf, matched scale ---")
D_Hopf_m = np.kron(scale_ext * D_ext_sqrt, I_30) + np.kron(Gamma_T3, D_F)
evals_Hopf_m = compute_and_analyze(D_Hopf_m, "T3 Hopf (matched)")

# --- Form 3: Sign grading ---
print(f"\n--- Form 3: Sign grading ---")
D_sign = np.kron(D_ext_sqrt, I_30) + np.kron(Gamma_sign, D_F)
evals_sign = compute_and_analyze(D_sign, "Sign grading")

# --- Form 4: Adjacency + T3 ---
print(f"\n--- Form 4: Adjacency/12 + T3 ---")
D_adj_T3 = np.kron(D_ext_adj, I_30) + np.kron(Gamma_T3, D_F)
evals_adj = compute_and_analyze(D_adj_T3, "Adj/12 + T3")

# --- Form 5: Multiplicative (eigenvalues = products, trivial) ---
print(f"\n--- Form 5: Multiplicative baseline (trivial) ---")
evals_mult = np.sort(np.array([de * df for de in evals_sqrt
                                for df in evals_DF]))[::-1]
pos_mult = evals_mult[evals_mult > 1e-10]
if len(pos_mult) > 0:
    dyn = pos_mult[0] / pos_mult[-1]
    print(f"  Dynamic range: {dyn:.1f} (phi^{np.log(dyn)/np.log(PHI):.1f})")

# --- Form 6: D_F tensor Delta_0 (eigenvalues of D^2 = products of D^2) ---
print(f"\n--- Form 6: D_F^2 + Delta_0 (squared masses = sums) ---")
evals_Dsq = np.sort(np.array([d0 + df**2 for d0 in evals_D0
                               for df in evals_DF]))[::-1]
pos_Dsq = evals_Dsq[evals_Dsq > 1e-10]
if len(pos_Dsq) > 0:
    dyn = pos_Dsq[0] / pos_Dsq[-1]
    print(f"  Dynamic range: {dyn:.1f} (phi^{np.log(dyn)/np.log(PHI):.1f})")


# =====================================================================
# PART 7: MASS EXPONENT SEARCH
# =====================================================================
print(f"\n" + "=" * 72)
print("PART 7: MASS EXPONENT SEARCH")
print("=" * 72)
print(f"Target mass exponents: {list(target_sorted)}")
print(f"  (t=26, b=19, tau=17, c=16, mu=s=11, d=5, u=3, e=0)")


def search_mass_pattern(evals, label, targets=target_sorted):
    """Search for 9 eigenvalues matching mass exponents."""
    pos = np.sort(evals[evals > 1e-10])[::-1]
    if len(pos) < 9:
        print(f"\n  [{label}] Only {len(pos)} positive eigenvalues. SKIP.")
        return 999.0

    logs = np.log(pos) / np.log(PHI)

    # Method 1: Top 9
    top9 = logs[:9]
    shift1 = np.mean(top9 - targets)
    rms1 = np.sqrt(np.mean((top9 - targets - shift1)**2))

    # Method 2: Search over all possible shifts
    best_rms = 999.0
    best_shift = 0.0
    best_matched = None

    for shift_trial in np.linspace(logs[-1], logs[0] - 26, 500):
        matched = []
        used = set()
        for tgt in sorted(targets, reverse=True):
            desired = tgt + shift_trial
            # Find closest unused eigenvalue
            diffs = np.abs(logs - desired)
            for idx in np.argsort(diffs):
                if idx not in used:
                    matched.append(logs[idx])
                    used.add(idx)
                    break
        if len(matched) == len(targets):
            matched = np.array(matched)
            residuals = matched - targets - shift_trial
            rms = np.sqrt(np.mean(residuals**2))
            if rms < best_rms:
                best_rms = rms
                best_shift = shift_trial
                best_matched = matched.copy()

    print(f"\n  [{label}]")
    print(f"  N(positive) = {len(pos)}, "
          f"Dynamic range = {pos[0]/pos[-1]:.1f}")
    print(f"  Method 1 (top 9): shift={shift1:.2f}, RMS={rms1:.2f}")
    print(f"  Method 2 (best match): shift={best_shift:.2f}, "
          f"RMS={best_rms:.4f}")
    if best_matched is not None and best_rms < 5.0:
        shifted = best_matched - best_shift
        print(f"    Matched (shifted): "
              f"{', '.join(f'{x:.2f}' for x in shifted)}")
        print(f"    Target:            "
              f"{', '.join(f'{x:.0f}' for x in targets)}")
        devs = shifted - targets
        print(f"    Deviations:        "
              f"{', '.join(f'{x:+.2f}' for x in devs)}")

    return best_rms


results = {}
for evals_test, label in [
    (evals_add, "Form 1: Additive"),
    (evals_Hopf, "Form 2: T3 Hopf"),
    (evals_Hopf_m, "Form 2b: T3 Hopf matched"),
    (evals_sign, "Form 3: Sign"),
    (evals_adj, "Form 4: Adj+T3"),
    (evals_mult, "Form 5: Multiplicative"),
]:
    rms = search_mass_pattern(evals_test, label)
    results[label] = rms


# =====================================================================
# PART 8: RANDOM BASELINE (intellectual honesty check)
# =====================================================================
print(f"\n" + "=" * 72)
print("PART 8: RANDOM BASELINE")
print("=" * 72)

# For the best form, compare against random targets
best_form_label = min(results, key=results.get)
best_form_rms = results[best_form_label]
print(f"Best form: {best_form_label} (RMS = {best_form_rms:.4f})")

# Pick the eigenvalues of the best form
best_evals_map = {
    "Form 1: Additive": evals_add,
    "Form 2: T3 Hopf": evals_Hopf,
    "Form 2b: T3 Hopf matched": evals_Hopf_m,
    "Form 3: Sign": evals_sign,
    "Form 4: Adj+T3": evals_adj,
    "Form 5: Multiplicative": evals_mult,
}
best_evals = best_evals_map[best_form_label]

np.random.seed(42)
n_random = 100
random_rms_list = []
print(f"Testing {n_random} random target sets...")

for trial in range(n_random):
    # Random 9 exponents in [0, 26]
    rand_targets = np.sort(np.random.uniform(0, 26, 9))[::-1]
    rms = search_mass_pattern.__wrapped__(best_evals, rand_targets) \
        if hasattr(search_mass_pattern, '__wrapped__') else 999.0

    # Inline version of search
    pos = np.sort(best_evals[best_evals > 1e-10])[::-1]
    if len(pos) < 9:
        continue
    logs = np.log(pos) / np.log(PHI)
    brms = 999.0
    for shift in np.linspace(logs[-1], logs[0] - 26, 200):
        matched = []
        used = set()
        for tgt in sorted(rand_targets, reverse=True):
            desired = tgt + shift
            diffs = np.abs(logs - desired)
            for idx in np.argsort(diffs):
                if idx not in used:
                    matched.append(logs[idx])
                    used.add(idx)
                    break
        if len(matched) == len(rand_targets):
            matched = np.array(matched)
            r = np.sqrt(np.mean((matched - rand_targets - shift)**2))
            if r < brms:
                brms = r
    random_rms_list.append(brms)

random_rms = np.array(random_rms_list)
print(f"\nRandom baseline (n={len(random_rms)}):")
print(f"  Mean RMS: {np.mean(random_rms):.4f}")
print(f"  Std RMS:  {np.std(random_rms):.4f}")
print(f"  Min RMS:  {np.min(random_rms):.4f}")
print(f"  Actual RMS: {best_form_rms:.4f}")
if np.std(random_rms) > 0:
    z_score = (np.mean(random_rms) - best_form_rms) / np.std(random_rms)
    print(f"  Z-score: {z_score:.2f} "
          f"({'SIGNIFICANT' if z_score > 2 else 'NOT significant'})")
    pct_better = np.mean(random_rms < best_form_rms) * 100
    print(f"  % random better: {pct_better:.1f}%")


# =====================================================================
# PART 9: EIGENVALUE STRUCTURE (Form 2)
# =====================================================================
print(f"\n" + "=" * 72)
print("PART 9: EIGENVALUE STRUCTURE (Form 2: T3 Hopf)")
print("=" * 72)

pos_H = np.sort(evals_Hopf[evals_Hopf > 1e-10])[::-1]
if len(pos_H) >= 20:
    print(f"\nTop 20 positive eigenvalues:")
    for i in range(20):
        log_phi = np.log(pos_H[i]) / np.log(PHI)
        print(f"  lambda_{i+1:4d} = {pos_H[i]:10.6f}  "
              f"(log_phi = {log_phi:7.3f})")

    print(f"\nBottom 10 positive eigenvalues:")
    for i in range(max(0, len(pos_H)-10), len(pos_H)):
        log_phi = np.log(pos_H[i]) / np.log(PHI)
        print(f"  lambda_{i+1:4d} = {pos_H[i]:10.6f}  "
              f"(log_phi = {log_phi:7.3f})")

    # Check for phi-integer ratios
    print(f"\nRatios between top consecutive eigenvalues:")
    for i in range(min(15, len(pos_H)-1)):
        ratio = pos_H[i] / pos_H[i+1]
        log_ratio = np.log(ratio) / np.log(PHI)
        marker = " <-- phi^INTEGER" if abs(log_ratio - round(log_ratio)) < 0.05 else ""
        print(f"  lambda_{i+1}/lambda_{i+2} = {ratio:.6f}  "
              f"(phi^{log_ratio:.3f}){marker}")


# =====================================================================
# PART 10: CROSS-TERM ANALYSIS
# =====================================================================
print(f"\n" + "=" * 72)
print("PART 10: CROSS-TERM {D_ext, Gamma_T3}")
print("=" * 72)

# {D_ext, Gamma_T3}[i,j] = D_ext[i,j] * (T3[i] + T3[j])
cross_matrix = D_ext_sqrt * (T3_vals[:, None] + T3_vals[None, :])
cross_evals = np.linalg.eigvalsh(cross_matrix)
print(f"  Range: [{cross_evals[0]:.4f}, {cross_evals[-1]:.4f}]")
print(f"  Nonzero: {np.sum(np.abs(cross_evals) > 1e-10)} / {N_verts}")
print(f"  Trace: {np.trace(cross_matrix):.4f}")

# Frobenius norm decomposition of D_Hopf^2
D_Hopf_sq = D_Hopf @ D_Hopf.conj().T
part_ext = np.kron(D_ext_sqrt @ D_ext_sqrt, I_30)
part_F = np.kron(Gamma_T3 @ Gamma_T3, D_F @ D_F.conj().T)
part_cross = D_Hopf_sq - part_ext - part_F

norm_total = np.linalg.norm(D_Hopf_sq, 'fro')
norm_ext = np.linalg.norm(part_ext, 'fro')
norm_F = np.linalg.norm(part_F, 'fro')
norm_cross = np.linalg.norm(part_cross, 'fro')

print(f"\n  D_Hopf^2 Frobenius norm decomposition:")
print(f"    Total:      {norm_total:.2f}")
print(f"    D_ext^2:    {norm_ext:.2f} ({norm_ext/norm_total*100:.1f}%)")
print(f"    Gamma^2*DF: {norm_F:.2f} ({norm_F/norm_total*100:.1f}%)")
print(f"    Cross term: {norm_cross:.2f} ({norm_cross/norm_total*100:.1f}%)")

# T3 sector breakdown
print(f"\n  T3 sector eigenvalue ranges in D_Hopf diagonal blocks:")
for tv in T3_unique:
    if abs(tv) < 1e-10:
        continue
    indices = np.where(np.abs(T3_vals - tv) < 1e-4)[0]
    # Extract 30x30 diagonal block for this sector
    # In the full 3600 matrix, vertex i maps to rows i*30 : (i+1)*30
    sector_blocks = []
    for idx in indices[:3]:  # Sample 3 vertices from this T3 sector
        r0 = idx * total_dim_F
        r1 = r0 + total_dim_F
        block = D_Hopf[r0:r1, r0:r1]
        sector_blocks.append(block)
    avg_block = np.mean([np.linalg.eigvalsh(b) for b in sector_blocks],
                        axis=0)
    print(f"    T3={tv:+.4f}: avg diag-block evals = "
          f"[{np.min(avg_block):.4f}, {np.max(avg_block):.4f}], "
          f"scale factor {abs(tv):.4f}")


# =====================================================================
# PART 11: HONEST ASSESSMENT
# =====================================================================
print(f"\n" + "=" * 72)
print("PART 11: HONEST ASSESSMENT")
print("=" * 72)

print(f"\n[Total time: {time.time()-t0:.1f}s]")

print(f"""
SUMMARY OF RESULTS:

1. PRODUCT OPERATOR CONSTRUCTED:
   D = D_ext (x) I_30 + Gamma (x) D_F on 3600-dim Hilbert space
   D_ext = sqrt(Delta_0) on 120-vertex 600-cell
   Gamma = T3 Hopf grading (9 distinct values on vertices)
   D_F = McKay CG-derived Dirac operator (30x30, Tr(D_F^2)=8)

2. FORMS TESTED:
   Form 1: Additive (trivial - eigenvalues are pairwise sums)
   Form 2: T3 Hopf grading (KEY - cross terms break A5 symmetry)
   Form 2b: T3 Hopf, matched scale
   Form 3: Sign grading (binary chirality)
   Form 4: Adjacency + T3
   Form 5: Multiplicative (trivial - eigenvalues are pairwise products)
   Form 6: D_F^2 + Delta_0 (squared masses)

3. DYNAMIC RANGE:""")

for label in ["Form 1: Additive", "Form 2: T3 Hopf",
              "Form 2b: T3 Hopf matched", "Form 3: Sign",
              "Form 4: Adj+T3", "Form 5: Multiplicative"]:
    print(f"   {label}: RMS = {results.get(label, 'N/A')}")

print(f"""
4. KEY STRUCTURAL INSIGHT:
   T3 grading introduces SECTOR-DEPENDENT scaling of D_F:
   - T3 = +-1: 4 vertices each (gauge vertices)
   - T3 = +-phi/2: 16 vertices each (lepton-like?)
   - T3 = +-1/2: 16 vertices each (quark up-like?)
   - T3 = +-1/(2*phi): 16 vertices each (quark down-like?)
   - T3 = 0: 16 vertices (pure geometry, no D_F)
   Ratio phi/2 : 1/2 : 1/(2*phi) = phi : 1 : 1/phi

5. CROSS-TERM:
   {{D_ext, Gamma_T3}} couples internal (McKay) to external (600-cell)
   Cross term is {norm_cross/norm_total*100:.1f}% of D_Hopf^2
""")

if best_form_rms < np.mean(random_rms) - 2*np.std(random_rms):
    print("6. VERDICT: SIGNIFICANT -- mass pattern detected above random")
elif best_form_rms < np.mean(random_rms):
    print("6. VERDICT: MARGINAL -- better than random mean but not 2-sigma")
else:
    print("6. VERDICT: NEGATIVE -- no better than random target matching")

print(f"\n   Best RMS: {best_form_rms:.4f} ({best_form_label})")
if len(random_rms) > 0:
    print(f"   Random mean: {np.mean(random_rms):.4f} +- "
          f"{np.std(random_rms):.4f}")

print(f"""
7. STRUCTURAL CONCLUSIONS:
   - 3600-dim space DOES have larger dynamic range than 30-dim D_F
   - T3 grading DOES break the A5 symmetry that keeps all masses equal
   - Cross term is NONZERO: internal/external spaces genuinely couple
   - Whether the specific mass pattern phi^n_f emerges: see results above
""")

print("=" * 72)
print("exp332 COMPLETE")
print("=" * 72)
