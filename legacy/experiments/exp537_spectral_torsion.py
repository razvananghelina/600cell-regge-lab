"""
EXP-537: Spectral Torsion on the McKay Graph
==============================================

First computation of spectral torsion (Dabrowski-Mukhopadhyay-Pozar 2025)
on the 600-cell's finite spectral triple.

SETUP:
  - Hilbert space H_F = direct sum of V_{rho_k}, dim = 30
  - Algebra A = C^9 (one scalar per irrep)
  - Dirac operator D_F: off-diagonal blocks Y_{kl} at McKay edges
  - Yukawa norm: ||Y_{kl}||_F = 1/sqrt(2) per edge

COMPUTE:
  1. D_F as explicit 30x30 matrix
  2. Tr(D_F^2), Tr(D_F^4)
  3. 1-forms Omega^1_D(A) = span{a[D,b]}
  4. Spectral Ricci scalar R_D(a) = Tr(a * D^2)
  5. Spectral torsion T_D(u,v,w) = Tr(u*v*w*D)
  6. Einstein tensor G_D(u,v) = Tr(u*{v,D}*D)
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, a1, b1, N, rank_E8, IRREP_DIMS_2I, MCKAY_EDGES)

print("=" * 70)
print("EXP-537: SPECTRAL TORSION ON THE McKAY GRAPH")
print("=" * 70)

# ============================================================
# SECTION 1: BUILD THE FINITE DIRAC OPERATOR D_F
# ============================================================
print(f"\n{'='*70}")
print("SECTION 1: FINITE DIRAC OPERATOR D_F")
print(f"{'='*70}")

dims = IRREP_DIMS_2I  # [1, 2, 3, 4, 5, 6, 4, 2, 3]
n_irreps = len(dims)
dim_H = sum(dims)  # = 30

print(f"  Irrep dims: {dims}")
print(f"  dim(H_F) = {dim_H}")
print(f"  McKay edges: {MCKAY_EDGES}")

# Block offsets in the 30-dim Hilbert space
offsets = [0]
for d in dims:
    offsets.append(offsets[-1] + d)
print(f"  Block offsets: {offsets}")

# Build D_F with universal Yukawa coupling
# For each edge (k,l): Y_{kl} is a dims[k] x dims[l] matrix with ||Y||_F = 1/sqrt(2)
# Choice: proportional to the "democratic" matrix (all entries equal)
# Y_{kl} = c * ones(dk, dl) where c = 1/sqrt(2*dk*dl)

D_F = np.zeros((dim_H, dim_H))

for (k, l) in MCKAY_EDGES:
    dk, dl = dims[k], dims[l]
    # Democratic Yukawa: all CG coefficients equal
    c = 1.0 / np.sqrt(2 * dk * dl)
    Y_kl = c * np.ones((dk, dl))

    # Place in D_F (upper block)
    D_F[offsets[k]:offsets[k]+dk, offsets[l]:offsets[l]+dl] = Y_kl
    # Hermitian conjugate (lower block)
    D_F[offsets[l]:offsets[l]+dl, offsets[k]:offsets[k]+dk] = Y_kl.T

# Verify Hermiticity
print(f"\n  ||D_F - D_F^T|| = {np.linalg.norm(D_F - D_F.T):.2e} (should be 0)")

# Check Yukawa norms
print(f"\n  Yukawa norms per edge:")
for (k, l) in MCKAY_EDGES:
    dk, dl = dims[k], dims[l]
    Y = D_F[offsets[k]:offsets[k]+dk, offsets[l]:offsets[l]+dl]
    norm_sq = np.trace(Y.T @ Y)
    print(f"    ({k},{l}): ||Y||_F^2 = {norm_sq:.6f} (should be 0.5)")

# ============================================================
# SECTION 2: SPECTRAL TRACES
# ============================================================
print(f"\n{'='*70}")
print("SECTION 2: SPECTRAL TRACES")
print(f"{'='*70}")

D2 = D_F @ D_F
D4 = D2 @ D2

tr_D2 = np.trace(D2)
tr_D4 = np.trace(D4)

print(f"  Tr(D_F^2) = {tr_D2:.6f}")
print(f"  Expected: 2 * |edges| * ||Y||^2 = 2*8*0.5 = {2*8*0.5}")
print(f"  = rank(E8) = {rank_E8}? {'YES' if abs(tr_D2 - rank_E8) < 0.01 else 'NO'}")

print(f"\n  Tr(D_F^4) = {tr_D4:.6f}")

# Eigenvalues of D_F
evals_D = np.linalg.eigvalsh(D_F)
evals_D.sort()
nonzero_evals = evals_D[np.abs(evals_D) > 1e-10]
print(f"\n  Eigenvalues of D_F: {len(nonzero_evals)} nonzero out of {dim_H}")
print(f"  Nonzero eigenvalues: {[f'{e:.4f}' for e in sorted(nonzero_evals)]}")
print(f"  Zero eigenvalues: {dim_H - len(nonzero_evals)}")

# ============================================================
# SECTION 3: ALGEBRA AND 1-FORMS
# ============================================================
print(f"\n{'='*70}")
print("SECTION 3: ALGEBRA A = C^9 AND 1-FORMS")
print(f"{'='*70}")

# Algebra element e_k = projection onto rho_k block
def algebra_element(coeffs):
    """Build algebra element a = diag(a_0*I_{d0}, a_1*I_{d1}, ...)"""
    a = np.zeros((dim_H, dim_H))
    for k in range(n_irreps):
        a[offsets[k]:offsets[k]+dims[k], offsets[k]:offsets[k]+dims[k]] = coeffs[k] * np.eye(dims[k])
    return a

# Basis projections
E = [algebra_element([1 if j == k else 0 for j in range(n_irreps)]) for k in range(n_irreps)]

# Commutator [D, e_k]
def commutator(A, B):
    return A @ B - B @ A

# 1-form basis: omega_{kl} = E_k @ [D, E_l] for each edge (k,l)
# [D, E_l] has nonzero blocks at edges incident to l
print(f"  Computing 1-form basis...")

one_forms = {}
for (k, l) in MCKAY_EDGES:
    # omega_{k->l}: projection onto "outgoing from k to l"
    dE_l = commutator(D_F, E[l])
    omega_kl = E[k] @ dE_l  # a[D,b] with a=E_k, b=E_l

    # Also the reverse
    dE_k = commutator(D_F, E[k])
    omega_lk = E[l] @ dE_k

    one_forms[(k, l)] = omega_kl
    one_forms[(l, k)] = omega_lk

print(f"  Number of 1-form generators: {len(one_forms)}")

# Check: are they nonzero?
print(f"\n  1-form norms:")
for key, omega in sorted(one_forms.items()):
    norm = np.linalg.norm(omega, 'fro')
    print(f"    omega_{key}: ||omega||_F = {norm:.6f}")

# ============================================================
# SECTION 4: SPECTRAL RICCI SCALAR
# ============================================================
print(f"\n{'='*70}")
print("SECTION 4: SPECTRAL RICCI SCALAR R_D(a) = Tr(a * D^2)")
print(f"{'='*70}")

# R_D(a) = Tr(a * D_F^2) for algebra elements a
print(f"  R_D(e_k) = Tr(E_k * D^2) for each basis element:")
print(f"  {'k':>3s} {'dim':>4s} {'R_D(e_k)':>10s} {'R_D/dim':>10s}")

R_values = []
for k in range(n_irreps):
    R_k = np.trace(E[k] @ D2)
    R_values.append(R_k)
    print(f"  {k:3d} {dims[k]:4d} {R_k:10.6f} {R_k/dims[k]:10.6f}")

print(f"\n  Sum R_D(e_k) = {sum(R_values):.6f} = Tr(D^2) = {tr_D2:.6f}")
print(f"  R_D(1) = Tr(D^2) = {tr_D2:.6f}")

# Per-irrep Ricci scalar (normalized by dimension)
print(f"\n  Per-irrep Ricci (R_k/dim_k):")
for k in range(n_irreps):
    bars = '#' * int(R_values[k]/dims[k] * 20)
    print(f"    rho_{k} (dim {dims[k]}): {R_values[k]/dims[k]:.4f} {bars}")

# ============================================================
# SECTION 5: SPECTRAL TORSION
# ============================================================
print(f"\n{'='*70}")
print("SECTION 5: SPECTRAL TORSION T_D(u,v,w) = Tr(u*v*w*D)")
print(f"{'='*70}")

# Compute T_D for all triples of 1-form basis elements
# T_{(e1),(e2),(e3)} = Tr(omega_e1 * omega_e2 * omega_e3 * D_F)

edge_keys = sorted(one_forms.keys())
n_1forms = len(edge_keys)

print(f"  Computing torsion for {n_1forms}^3 = {n_1forms**3} triples...")

# First: check a few specific triples
nonzero_count = 0
max_torsion = 0
max_triple = None
torsion_values = []

for i, e1 in enumerate(edge_keys):
    for j, e2 in enumerate(edge_keys):
        for k_idx, e3 in enumerate(edge_keys):
            T = np.trace(one_forms[e1] @ one_forms[e2] @ one_forms[e3] @ D_F)
            if abs(T) > 1e-12:
                nonzero_count += 1
                torsion_values.append((e1, e2, e3, T))
                if abs(T) > max_torsion:
                    max_torsion = abs(T)
                    max_triple = (e1, e2, e3, T)

print(f"  Nonzero torsion values: {nonzero_count} / {n_1forms**3}")
print(f"  Max |T|: {max_torsion:.8f}")

if max_triple:
    print(f"  Max triple: T({max_triple[0]},{max_triple[1]},{max_triple[2]}) = {max_triple[3]:.8f}")

# Show top 10 nonzero torsions
if torsion_values:
    torsion_values.sort(key=lambda x: -abs(x[3]))
    print(f"\n  Top nonzero torsion values:")
    print(f"  {'(e1)':>8s} {'(e2)':>8s} {'(e3)':>8s} {'T_D':>14s}")
    for e1, e2, e3, T in torsion_values[:15]:
        print(f"  {str(e1):>8s} {str(e2):>8s} {str(e3):>8s} {T:+14.8f}")

# Total torsion (sum of all |T|)
total_torsion = sum(abs(T) for _, _, _, T in torsion_values)
print(f"\n  Total |torsion| = {total_torsion:.8f}")

# ============================================================
# SECTION 6: TORSION STRUCTURE
# ============================================================
print(f"\n{'='*70}")
print("SECTION 6: TORSION STRUCTURE AND ALGEBRAIC CONTENT")
print(f"{'='*70}")

if torsion_values:
    # Extract unique |T| values
    unique_T = sorted(set(round(abs(T), 10) for _, _, _, T in torsion_values))
    print(f"  Distinct |T| values: {len(unique_T)}")
    for t_val in unique_T[:10]:
        count = sum(1 for _, _, _, T in torsion_values if abs(abs(T) - t_val) < 1e-9)
        # Check if related to framework constants
        ratio_phi = t_val / (1/PHI) if t_val > 0 else 0
        ratio_alpha = t_val / (1/137.036) if t_val > 0 else 0
        print(f"    |T| = {t_val:.10f} (x{count}), "
              f"T*sqrt(2) = {t_val*np.sqrt(2):.6f}, "
              f"T*30 = {t_val*30:.6f}")

# ============================================================
# SECTION 7: EINSTEIN TENSOR
# ============================================================
print(f"\n{'='*70}")
print("SECTION 7: EINSTEIN TENSOR G_D(u,v)")
print(f"{'='*70}")

# G_D(u,v) = Tr(u * {v, D} * D) = Tr(u*v*D^2) + Tr(u*D*v*D)
# G'_D(u,v) = Tr(u*v*D^2)
# G''_D(u,v) = Tr(u*D*v*D)

# Compute for a few pairs of 1-forms
if one_forms:
    print(f"  G'_D and G''_D for edge pairs:")
    print(f"  {'(e1)':>8s} {'(e2)':>8s} {'Gprime':>12s} {'Gdouble':>12s} {'G_total':>12s}")

    for e1 in list(edge_keys)[:6]:
        for e2 in list(edge_keys)[:6]:
            u = one_forms[e1]
            v = one_forms[e2]
            G_prime = np.trace(u @ v @ D2)
            G_double = np.trace(u @ D_F @ v @ D_F)
            G_total = G_prime + G_double
            if abs(G_total) > 1e-12:
                print(f"  {str(e1):>8s} {str(e2):>8s} {G_prime:+12.6f} {G_double:+12.6f} {G_total:+12.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  FINITE SPECTRAL TRIPLE ON McKAY E8~:
    dim(H_F) = {dim_H}
    Tr(D_F^2) = {tr_D2:.4f} {'= rank(E8) = 8' if abs(tr_D2-8)<0.01 else ''}
    Tr(D_F^4) = {tr_D4:.4f}
    Nonzero eigenvalues of D_F: {len(nonzero_evals)}

  SPECTRAL RICCI SCALAR:
    R_D(1) = Tr(D^2) = {tr_D2:.4f}
    Per-irrep: varies from {min(R_values[k]/dims[k] for k in range(n_irreps)):.4f}
               to {max(R_values[k]/dims[k] for k in range(n_irreps)):.4f}

  SPECTRAL TORSION:
    Nonzero triples: {nonzero_count} / {n_1forms**3}
    Max |T|: {max_torsion:.8f}
    TORSION IS {'NONZERO' if nonzero_count > 0 else 'ZERO'}

  INTERPRETATION:
    The spectral torsion measures how the Yukawa couplings "twist"
    the internal geometry. Nonzero torsion means the geometry is
    NOT torsion-free -- the Yukawa couplings induce genuine curvature.
    This is the first computation of this quantity on the 600-cell
    finite spectral triple.
""")

print("=" * 70)
print("EXP-537 COMPLETE")
print("=" * 70)
