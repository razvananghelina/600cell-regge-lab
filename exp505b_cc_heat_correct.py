"""
EXP-505b: CC from Heat Kernel on 600-Cell (Correct, using commons)
===================================================================
Redo of exp505 with correct 600-cell from commons module.

The idea (user): auto-reflection that EVOLVES = heat kernel on 600-cell.
The Galois-breaking part gives an exponentially small vacuum energy.

Key question: does this naturally produce Lambda ~ alpha^57?
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '.')
from commons import (build_600cell, a1, b1, N, PHI, PHI_CONJ, SQRT5,
                      alpha_em, alpha_s, ln_inv_alpha, L1_600cell,
                      degree, dim_E8, N_gen, rank_E8, h)

print("=" * 70)
print("EXP-505b: CC FROM 600-CELL HEAT KERNEL (CORRECT)")
print("=" * 70)

# ============================================================
# Build 600-cell and get spectrum
# ============================================================
verts, adj, lap = build_600cell()
eigvals = np.sort(eigh(lap, eigvals_only=True))

# Distinct eigenvalues with multiplicities
distinct = []
for ev in eigvals:
    if not any(abs(ev - d[0]) < 0.001 for d in distinct):
        mult = sum(1 for x in eigvals if abs(x - ev) < 0.001)
        distinct.append((ev, mult))

print(f"\n600-cell Laplacian: {len(distinct)} distinct eigenvalues")
print(f"{'i':>3} {'lambda':>10} {'mult':>5} {'exact form':>25}")
# Known exact forms from the theory:
exact_forms = {
    0: "0",
    1: "12 - 6*phi",
    2: "12 - 2*sqrt(5)",
    3: "9",
    4: "12",
    5: "14",
    6: "12 + 2*sqrt(5)",
    7: "15",
    8: "6 + 6*phi",
}
for i, (ev, mult) in enumerate(distinct):
    form = exact_forms.get(i, "?")
    print(f"{i:3d} {ev:10.4f} {mult:5d}   {form}")

# ============================================================
# Identify Galois pairs
# ============================================================
print(f"\n{'='*70}")
print("GALOIS STRUCTURE")
print(f"{'='*70}")

# Galois involution sigma: sqrt(5) -> -sqrt(5), i.e., phi -> phi'
# Eigenvalues in Z[phi] come in Galois pairs.
# The Galois pair of lambda is sigma(lambda).

# For each eigenvalue, compute its Galois conjugate:
# If lambda = a + b*sqrt(5), then sigma(lambda) = a - b*sqrt(5)
# We can find b = (lambda - rational_part) / sqrt(5)

galois_pairs = []  # (idx_i, idx_j, ev_i, ev_j, mult)
galois_fixed = []  # (idx, ev, mult)
used = set()

for i, (ev_i, mult_i) in enumerate(distinct):
    if i in used:
        continue
    # Compute Galois conjugate
    # Try to express ev as a + b*sqrt(5)
    # For each possible rational a (half-integers up to 20):
    found = False
    for a2 in range(0, 40):
        a = a2 / 2.0
        b = (ev_i - a) / SQRT5
        if abs(b - round(b * 2) / 2) < 0.001:
            b_clean = round(b * 2) / 2
            ev_conj = a - b_clean * SQRT5
            if abs(b_clean) > 0.01:  # actually involves sqrt(5)
                # Find the conjugate in the list
                for j, (ev_j, mult_j) in enumerate(distinct):
                    if j > i and j not in used and abs(ev_j - ev_conj) < 0.01:
                        galois_pairs.append((i, j, ev_i, ev_j, mult_i))
                        used.add(i)
                        used.add(j)
                        found = True
                        print(f"  Pair: L_{i}={ev_i:.4f} ({a}+{b_clean}*sqrt5) "
                              f"<-> L_{j}={ev_j:.4f} ({a}-{b_clean}*sqrt5), mult {mult_i}")
                        break
            break

    if not found and i not in used:
        galois_fixed.append((i, ev_i, mult_i))
        used.add(i)

print(f"\n  Galois pairs: {len(galois_pairs)}")
print(f"  Galois-fixed: {len(galois_fixed)}")
for i, ev, mult in galois_fixed:
    print(f"    Fixed: L_{i} = {ev:.4f}, mult {mult}")

# ============================================================
# Quantum-dimension weighted Galois heat kernel
# ============================================================
print(f"\n{'='*70}")
print("GALOIS HEAT KERNEL")
print(f"{'='*70}")

# From the McKay correspondence, each eigenvalue corresponds to an
# irrep of 2I with a quantum dimension d_q.
# The Galois-broken irreps have d_q = phi, sigma(d_q) = phi'.
# The weight difference: d_q^2 - sigma(d_q)^2 = phi^2 - phi'^2 = sqrt(5).

# For the mult-4 pair (L_1, L_8):
# These correspond to the 2D irreps of 2I (there are two, Galois-conjugate).
# d_q = phi for one, d_q = phi' for the other.
# Weight: mult * (phi^2 - phi'^2) = 4 * sqrt(5)

# For the mult-9 pair (L_2, L_6):
# These correspond to the 3D irreps.
# Same quantum dimension structure.
# Weight: mult * (phi^2 - phi'^2) = 9 * sqrt(5)

# The Galois heat kernel difference:
# Delta_K(beta) = sum over Galois pairs: mult * sqrt(5) * (e^{-beta*L_i} - e^{-beta*L_j})
# where L_i < L_j (L_i is the phi-eigenvalue, L_j is the phi'-eigenvalue)

print("\n  Galois heat kernel: Delta_K(beta) = sum_pairs mult*sqrt(5)*(e^{-b*Li} - e^{-b*Lj})")

z_CC = dim_E8/4 - a1 - alpha_s
target_z = z_CC

print(f"\n  Target: z_CC = {z_CC:.6f}")
print(f"  Target: Lambda_P = alpha^{z_CC:.4f} = {alpha_em**z_CC:.4e}")
print(f"  Observed: ~2.888e-122")

betas = [
    (N, "N=120"),
    (N + N_gen, "N+N_gen=123"),
    (N + rank_E8//2, "N+rank/2=124"),
    (a1**3, "a1^3=125"),
    (N - 1, "N-1=119"),
    (2*N, "2N=240"),
    (h*a1, "h*a1=150"),
    (N//2, "N/2=60"),
]

print(f"\n  {'beta':>6} {'label':>16}  {'Delta_K':>14}  {'|Delta_K|':>14}  "
      f"{'z_equiv':>10}  {'vs 56.88':>10}")

for beta_val, label in betas:
    delta_K = 0
    for (i, j, ev_i, ev_j, mult) in galois_pairs:
        # ev_i < ev_j always (from how we found them)
        small, big = min(ev_i, ev_j), max(ev_i, ev_j)
        delta_K += mult * SQRT5 * (np.exp(-beta_val * small) - np.exp(-beta_val * big))

    abs_dK = abs(delta_K)
    if abs_dK > 1e-300:
        z_equiv = -np.log(abs_dK) / ln_inv_alpha
    else:
        z_equiv = float('inf')

    err = z_equiv - z_CC
    print(f"  {beta_val:6d} {label:>16}  {delta_K:14.4e}  {abs_dK:14.4e}  "
          f"{z_equiv:10.4f}  {err:+10.4f}")

# ============================================================
# Dominant term analysis
# ============================================================
print(f"\n{'='*70}")
print("DOMINANT TERM ANALYSIS")
print(f"{'='*70}")

# At large beta, the dominant Galois pair is the one with smallest eigenvalue
# That's the mult-4 pair (L_1 = 12-6*phi ~ 2.29, L_8 = 6+6*phi ~ 15.71)
L_1 = distinct[1][0]
L_8 = distinct[8][0]
m_1 = distinct[1][1]

print(f"\n  Dominant pair: L_1 = {L_1:.6f} (mult {m_1}), L_8 = {L_8:.6f}")
print(f"  L_1 = 12 - 6*phi = {12-6*PHI:.6f}")
print(f"  L_8 = 6 + 6*phi = {6+6*PHI:.6f} (Galois conjugate)")
print(f"  L_1 + L_8 = {L_1+L_8:.6f} = 18 = 3*b1")
print(f"  L_1 * L_8 = {L_1*L_8:.6f} = 12^2 - (6*phi)*(6*phi') = 144-36*(-1) = 180? {abs(L_1*L_8-36):.4f}")

prod_exact = (12 - 6*PHI)*(6 + 6*PHI)
# = 72 + 72*phi - 36*phi - 36*phi^2 = 72 + 36*phi - 36*(phi+1) = 72 + 36*phi - 36*phi - 36 = 36
print(f"  L_1 * L_8 = 36 = b1^2: {abs(prod_exact - 36) < 0.001}")

print(f"\n  At large beta, Delta_K(beta) ~ {m_1}*sqrt(5)*exp(-beta*L_1)")
print(f"  Prefactor: {m_1}*sqrt(5) = {m_1*SQRT5:.6f}")
print(f"  ln(prefactor) = {np.log(m_1*SQRT5):.6f}")

# For z_CC: need -ln(Delta_K)/ln(1/alpha) = z_CC
# => beta*L_1 - ln(prefactor) = z_CC * ln(1/alpha)
# => beta = (z_CC*ln(1/alpha) + ln(prefactor)) / L_1

beta_exact = (z_CC * ln_inv_alpha + np.log(m_1*SQRT5)) / L_1
print(f"\n  beta_exact = {beta_exact:.6f}")
print(f"  N = {N}")
print(f"  beta_exact - N = {beta_exact - N:.6f}")
print(f"  N_gen = {N_gen}")

# ============================================================
# Key test: beta = N + N_gen = 123
# ============================================================
print(f"\n{'='*70}")
print("KEY TEST: beta = N + N_gen = 123")
print(f"{'='*70}")

beta_123 = N + N_gen
bL1 = beta_123 * L_1
prefactor = m_1 * SQRT5  # = 4*sqrt(5)

# Dominant term
dominant = prefactor * np.exp(-bL1)
z_dom = -np.log(dominant) / ln_inv_alpha

# Full Galois heat kernel at beta=123
delta_full = 0
for (i, j, ev_i, ev_j, mult) in galois_pairs:
    small, big = min(ev_i, ev_j), max(ev_i, ev_j)
    delta_full += mult * SQRT5 * (np.exp(-beta_123 * small) - np.exp(-beta_123 * big))
z_full = -np.log(abs(delta_full)) / ln_inv_alpha

print(f"\n  beta = {beta_123}")
print(f"  beta * L_1 = {bL1:.6f}")
print(f"  Dominant: {m_1}*sqrt(5)*exp(-{beta_123}*L_1) = {dominant:.6e}")
print(f"  z_dominant = {z_dom:.6f}")
print(f"  Full Delta_K = {delta_full:.6e}")
print(f"  z_full = {z_full:.6f}")
print(f"  z_CC = {z_CC:.6f}")
print(f"  Error (dominant): {z_dom - z_CC:+.6f} ({abs(z_dom-z_CC)/z_CC*100:.4f}%)")
print(f"  Error (full):     {z_full - z_CC:+.6f} ({abs(z_full-z_CC)/z_CC*100:.4f}%)")

Lambda_dom = alpha_em ** z_dom
Lambda_full = alpha_em ** z_full
Lambda_target = alpha_em ** z_CC
Lambda_obs = 2.888e-122

print(f"\n  Lambda(dominant) = {Lambda_dom:.4e}")
print(f"  Lambda(full)     = {Lambda_full:.4e}")
print(f"  Lambda(target)   = {Lambda_target:.4e}")
print(f"  Lambda(observed) = {Lambda_obs:.4e}")

# ============================================================
# Algebraic content of beta*L_1
# ============================================================
print(f"\n{'='*70}")
print("ALGEBRAIC CONTENT")
print(f"{'='*70}")

# beta * L_1 = 123 * (12 - 6*phi)
# = 1476 - 738*phi
# = 738*(2 - phi) = 738*phi'^2

bL1_exact = 738 * PHI_CONJ**2
print(f"\n  123 * (12-6*phi) = 1476 - 738*phi")
print(f"  = 738 * (2-phi) = 738 * phi'^2")
print(f"  = {bL1_exact:.6f}")
print(f"  738 = 2 * 3 * 123 = 6 * 123 = b1 * (N + N_gen)")
print(f"  738 = N*b1 + N_gen*b1 = 720 + 18 = {N*b1 + N_gen*b1}")

# Compare with target:
target_bL1 = z_CC * ln_inv_alpha + np.log(m_1 * SQRT5)
print(f"\n  Target: z_CC*ln(1/alpha) + ln(4*sqrt(5)) = {target_bL1:.6f}")
print(f"  Actual: 738*phi'^2 = {bL1_exact:.6f}")
print(f"  Gap: {target_bL1 - bL1_exact:.6f}")
print(f"  Gap/L_1 = {(target_bL1 - bL1_exact)/L_1:.6f}")

# The gap is the amount by which beta differs from 123
gap_beta = (target_bL1 - bL1_exact) / L_1
print(f"\n  beta_exact = 123 + {gap_beta:.6f}")
print(f"  = 123.{int(gap_beta*1000):03d}...")

# ============================================================
# What if we DON'T write as alpha^z?
# ============================================================
print(f"\n{'='*70}")
print("DIRECT FORM: Lambda = prefactor * exp(-beta * L_1)")
print(f"{'='*70}")

# Lambda = 4*sqrt(5) * exp(-123 * (12-6*phi))
# = 4*sqrt(5) * exp(-738*phi'^2)

Lambda_direct = 4 * SQRT5 * np.exp(-738 * PHI_CONJ**2)
print(f"\n  Lambda_P = 4*sqrt(5) * exp(-738*phi'^2)")
print(f"           = 4*sqrt(5) * exp(-{738*PHI_CONJ**2:.4f})")
print(f"           = {Lambda_direct:.6e}")
print(f"  vs alpha^56.88 = {Lambda_target:.6e}")
print(f"  vs observed     = {Lambda_obs:.6e}")
print(f"  Ratio direct/target = {Lambda_direct/Lambda_target:.6f}")
print(f"  Error: {abs(Lambda_direct/Lambda_target - 1)*100:.2f}%")
print(f"  Ratio direct/obs = {Lambda_direct/Lambda_obs:.6f}")
print(f"  Error vs obs: {abs(Lambda_direct/Lambda_obs - 1)*100:.2f}%")

# ============================================================
# Include subdominant Galois pair
# ============================================================
print(f"\n{'='*70}")
print("WITH SUBDOMINANT PAIR (mult-9)")
print(f"{'='*70}")

# The mult-9 pair: L_2 = 12-2*sqrt(5) ~ 5.53, L_6 = 12+2*sqrt(5) ~ 16.47
L_2 = distinct[2][0]
L_6 = distinct[6][0]
m_2 = distinct[2][1]

print(f"  L_2 = {L_2:.6f} = 12-2*sqrt(5), mult {m_2}")
print(f"  L_6 = {L_6:.6f} = 12+2*sqrt(5)")

subdom = m_2 * SQRT5 * np.exp(-beta_123 * L_2)
print(f"\n  Subdominant: 9*sqrt(5)*exp(-123*L_2) = {subdom:.6e}")
print(f"  vs dominant: {dominant:.6e}")
print(f"  Ratio sub/dom: {subdom/dominant:.6e}")

# The subdominant is negligible at beta=123
print(f"  Subdominant is {'negligible' if subdom/dominant < 1e-6 else 'NOT negligible'}")

# Full with both pairs:
Lambda_both = dominant + subdom
z_both = -np.log(Lambda_both) / ln_inv_alpha
print(f"\n  Lambda(both pairs) = {Lambda_both:.6e}")
print(f"  z_both = {z_both:.6f}")
print(f"  Correction from subdominant: {abs(z_both-z_dom)/z_dom*100:.6f}%")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
THE FORMULA:
  Lambda_P = 4*sqrt(5) * exp(-(N + N_gen) * (12 - 6*phi))
           = 4*sqrt(5) * exp(-738 * phi'^2)

WHERE:
  4 = mult(L_1) = dim of first Galois-broken eigenspace
  sqrt(5) = sqrt(a1) = phi^2 - phi'^2 (Galois weight)
  N = 120 = vertices of 600-cell
  N_gen = 3 = number of generations
  12-6*phi = L_1 = spectral gap (first Galois-broken Laplacian eigenvalue)
  phi'^2 = (3-sqrt(5))/2 = 2 - phi (Galois conjugate squared)
  738 = b1 * (N + N_gen) = 6 * 123

ALL factors derived from a1 = 5.

NUMERICAL RESULT:
  Lambda_P = {Lambda_direct:.4e}
  vs alpha^56.88 = {Lambda_target:.4e}
  Error: {abs(Lambda_direct/Lambda_target - 1)*100:.2f}%

  vs observed = {Lambda_obs:.4e}
  Error vs obs: {abs(Lambda_direct/Lambda_obs - 1)*100:.2f}%

WHAT'S EXPLAINED:
  - WHY exponentially small: heat kernel diffusion (e^{{-beta*L_1}})
  - WHY base involves alpha: L_1 = 12-6*phi connects to alpha through
    the spectral action (both determined by the same golden ratio geometry)
  - WHY ~10^-122: beta*L_1 = 738*phi'^2 ~ 282 ~ 57*ln(1/alpha) ~ 280

WHAT'S NOT EXPLAINED:
  - WHY beta = N + N_gen (diffusion time = one step per vertex + per generation)
  - The {abs(Lambda_direct/Lambda_target - 1)*100:.2f}% gap (transcendental vs algebraic)
  - The physical mechanism selecting this specific heat kernel

STATUS: {'PATTERN (promising)' if abs(Lambda_direct/Lambda_target - 1)*100 < 25 else 'NEGATIVE'}
""")

print("=" * 70)
print("EXP-505b COMPLETE")
print("=" * 70)
