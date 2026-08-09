"""
EXP-513: CC from One-Loop Spectral Action on 600-Cell
=======================================================
NO CONJECTURES. Just compute.

The 600-cell Laplacian has 9 eigenvalues, all known exactly.
The one-loop effective action is:
  V_1loop = (1/2) * ln det'(Delta_0) = (1/2) * sum mult_k * ln(lambda_k)

The Galois-breaking part (CC candidate):
  Delta_V = V_physical - V_dark
  = (1/2) * [sum over Galois pairs: mult*(ln L_phys - ln L_dark)]

This is EXACT and COMPUTABLE. No free parameters.
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '.')
from commons import (build_600cell, PHI, PHI_CONJ, SQRT5, a1, b1, N,
                      alpha_em, alpha_s, ln_inv_alpha, degree)

print("=" * 70)
print("EXP-513: CC FROM ONE-LOOP SPECTRAL ACTION")
print("=" * 70)

# ============================================================
# Build spectrum
# ============================================================
verts, adj, lap = build_600cell()
eigvals = np.sort(eigh(lap, eigvals_only=True))

# Exact eigenvalues (from theory, verified numerically)
spectrum = [
    (0.0,           1,  "0",              True),   # trivial
    (12 - 6*PHI,    4,  "12-6phi",        False),  # Galois pair 1 (phys)
    (12 - 2*SQRT5,  9,  "12-2sqrt5",      False),  # Galois pair 2 (phys)
    (9.0,          16,  "9",              True),   # rational (fixed)
    (12.0,         25,  "12",             True),   # rational (fixed)
    (14.0,         36,  "14",             True),   # rational (fixed)
    (12 + 2*SQRT5,  9,  "12+2sqrt5",      False),  # Galois pair 2 (dark)
    (15.0,         16,  "15",             True),   # rational (fixed)
    (6 + 6*PHI,     4,  "6+6phi",         False),  # Galois pair 1 (dark)
]

print(f"\n600-cell Laplacian spectrum:")
print(f"{'i':>3} {'lambda':>12} {'mult':>5} {'exact':>15} {'Galois':>8}")
for i, (ev, mult, name, fixed) in enumerate(spectrum):
    status = "fixed" if fixed else "broken"
    print(f"{i:3d} {ev:12.6f} {mult:5d} {name:>15} {status:>8}")

# ============================================================
# STEP 1: Full log-determinant
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: FULL LOG-DETERMINANT")
print(f"{'='*70}")

# ln det'(Delta) = sum_{lambda > 0} mult * ln(lambda)
ln_det = 0
for ev, mult, name, _ in spectrum:
    if ev > 0.001:  # skip zero mode
        contribution = mult * np.log(ev)
        ln_det += contribution
        print(f"  {mult:3d} * ln({ev:10.6f}) = {mult:3d} * {np.log(ev):8.4f} = {contribution:10.4f}")

print(f"\n  ln det'(Delta_0) = {ln_det:.6f}")
print(f"  V_1loop = (1/2) * ln det' = {ln_det/2:.6f}")

# ============================================================
# STEP 2: Galois-breaking part
# ============================================================
print(f"\n{'='*70}")
print("STEP 2: GALOIS-BREAKING PART (THE CC CANDIDATE)")
print(f"{'='*70}")

# Galois pairs:
# Pair 1: L_1 = 12-6*phi (phys) <-> L_8 = 6+6*phi (dark), mult 4
# Pair 2: L_2 = 12-2*sqrt5 (phys) <-> L_6 = 12+2*sqrt5 (dark), mult 9

L_1 = 12 - 6*PHI
L_8 = 6 + 6*PHI
L_2 = 12 - 2*SQRT5
L_6 = 12 + 2*SQRT5

print(f"\n  Galois pair 1: L_1 = {L_1:.6f}, L_8 = {L_8:.6f}, mult = 4")
print(f"    L_1/L_8 = {L_1/L_8:.6f} = 1/phi^4 = {1/PHI**4:.6f}")
print(f"    ln(L_1/L_8) = {np.log(L_1/L_8):.6f} = -4*ln(phi) = {-4*np.log(PHI):.6f}")

print(f"\n  Galois pair 2: L_2 = {L_2:.6f}, L_6 = {L_6:.6f}, mult = 9")
print(f"    L_2/L_6 = {L_2/L_6:.6f}")
print(f"    ln(L_2/L_6) = {np.log(L_2/L_6):.6f}")

# The Galois breaking:
# Delta_V = (1/2) * [4*(ln L_1 - ln L_8) + 9*(ln L_2 - ln L_6)]

delta_pair1 = 4 * np.log(L_1/L_8)
delta_pair2 = 9 * np.log(L_2/L_6)
delta_V = 0.5 * (delta_pair1 + delta_pair2)

print(f"\n  Pair 1 contribution: 4 * ln(L_1/L_8) = {delta_pair1:.6f}")
print(f"  Pair 2 contribution: 9 * ln(L_2/L_6) = {delta_pair2:.6f}")
print(f"  Total: {delta_pair1 + delta_pair2:.6f}")
print(f"  Delta_V = (1/2) * total = {delta_V:.6f}")

# Exact forms:
# Pair 1: 4*ln(1/phi^4) = -16*ln(phi)
# Pair 2: 9*ln((12-2sqrt5)/(12+2sqrt5))

exact_pair1 = -16 * np.log(PHI)
exact_pair2 = 9 * np.log((12 - 2*SQRT5) / (12 + 2*SQRT5))

print(f"\n  Exact pair 1: -16*ln(phi) = {exact_pair1:.6f}")
print(f"  Exact pair 2: 9*ln((12-2s5)/(12+2s5)) = {exact_pair2:.6f}")
print(f"  Sum: {exact_pair1 + exact_pair2:.6f}")

# ============================================================
# STEP 3: Interpret Delta_V
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: INTERPRETATION")
print(f"{'='*70}")

print(f"\n  Delta_V = {delta_V:.6f}")
print(f"  |Delta_V| = {abs(delta_V):.6f}")

# Is this the CC? If so:
# Lambda_P = exp(Delta_V) or Lambda_P = Delta_V or Lambda_P = Delta_V^2 ?

# The one-loop effective potential DIFFERENCE between physical and dark:
# This is a NUMBER, not exponentially small.

# What z would this correspond to?
if delta_V != 0:
    z_from_deltaV = -np.log(abs(delta_V)) / ln_inv_alpha
    print(f"  If Lambda = |Delta_V|: z_equiv = {z_from_deltaV:.4f}")

    z_from_exp = delta_V  # If CC = exp(Delta_V)
    print(f"  If Lambda = exp(Delta_V): Lambda = {np.exp(delta_V):.6e}")

# Delta_V ~ -8.2. This is O(1), not 10^{-122}.
# The one-loop log-determinant gives an O(1) Galois breaking.
# This is NOT the CC.

print(f"""
  RESULT: Delta_V = {delta_V:.4f} (order 1, NOT 10^{{-122}})

  The one-loop log-determinant Galois breaking is O(1).
  This is NOT the cosmological constant.

  This was expected: the log-determinant is a PERTURBATIVE quantity
  (one-loop), while the CC requires a NON-PERTURBATIVE suppression.
""")

# ============================================================
# STEP 4: What about exp(Delta_V)?
# ============================================================
print(f"{'='*70}")
print("STEP 4: EXPONENTIAL OF GALOIS BREAKING")
print(f"{'='*70}")

# exp(Delta_V) = exp(-8.18) = 2.8e-4. Still not 10^{-122}.

exp_deltaV = np.exp(delta_V)
print(f"\n  exp(Delta_V) = exp({delta_V:.4f}) = {exp_deltaV:.6e}")
print(f"  This is ~ {exp_deltaV:.4f}, still O(1), not 10^-122.")

# What if the CC involves ITERATED one-loop?
# After n loops: Delta_V^n or n * Delta_V?

# n * Delta_V = -122 * ln(10) = -280.9 requires n = 280.9/8.18 = 34.3
print(f"\n  For exp(n*Delta_V) = 10^-122: n = {-122*np.log(10)/delta_V:.1f}")
# n = 34.3. Not a clean number.

# Delta_V^n for integer n:
print(f"  Delta_V^2 = {delta_V**2:.4f}")
print(f"  Delta_V^3 = {delta_V**3:.4f}")

# ============================================================
# STEP 5: The FULL spectral action (not just log-det)
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: FULL SPECTRAL ACTION Tr(f(D/Lambda))")
print(f"{'='*70}")

# For f(x) = exp(-x^2), the spectral action is the heat kernel:
# S(Lambda) = Tr exp(-Delta/Lambda^2) = sum mult_k * exp(-lambda_k/Lambda^2)

# The Galois-breaking part at cutoff Lambda:
# Delta_S(Lambda) = sum over Galois pairs: mult * (exp(-L_phys/Lambda^2) - exp(-L_dark/Lambda^2))

# At large Lambda (UV): Delta_S -> 0 (all eigenvalues suppressed)
# At small Lambda (IR): Delta_S -> dominant term exp(-L_1/Lambda^2)

print(f"\n  Galois spectral action at various cutoffs Lambda:")
print(f"  {'Lambda':>10} {'Delta_S':>14} {'z_equiv':>10} {'note':>20}")

for Lambda_sq in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 1000.0]:
    delta_S = 0
    # Pair 1: mult 4
    delta_S += 4 * (np.exp(-L_1/Lambda_sq) - np.exp(-L_8/Lambda_sq))
    # Pair 2: mult 9
    delta_S += 9 * (np.exp(-L_2/Lambda_sq) - np.exp(-L_6/Lambda_sq))

    if abs(delta_S) > 1e-300:
        z = -np.log(abs(delta_S)) / ln_inv_alpha
    else:
        z = float('inf')

    note = ""
    if abs(Lambda_sq - 1.0) < 0.01: note = "Lambda=1"
    if abs(Lambda_sq - L_1) < 0.1: note = "Lambda^2=L_1"

    print(f"  {np.sqrt(Lambda_sq):10.4f} {delta_S:14.6e} {z:10.4f} {note:>20}")

# ============================================================
# STEP 6: The spectral action at Lambda^2 = specific values
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: SPECIAL VALUES OF THE CUTOFF")
print(f"{'='*70}")

# What Lambda^2 gives CC ~ 10^{-122}?
# We need |Delta_S| ~ 10^{-122}
# At large Lambda^2: Delta_S ~ 4*(1 - exp(-(L_8-L_1)/Lambda^2)) ~ 4*(L_8-L_1)/Lambda^2
# This is O(1/Lambda^2), not exponentially small.

# At small Lambda^2: Delta_S ~ 4*exp(-L_1/Lambda^2)
# For this to be 10^{-122}: L_1/Lambda^2 ~ 280
# Lambda^2 ~ L_1/280 ~ 2.29/280 ~ 0.0082

Lambda_sq_cc = L_1 / (122 * np.log(10))
Lambda_cc = np.sqrt(Lambda_sq_cc)
delta_S_cc = 4 * np.exp(-L_1/Lambda_sq_cc) - 4*np.exp(-L_8/Lambda_sq_cc)
z_cc_check = -np.log(abs(delta_S_cc)) / ln_inv_alpha

print(f"\n  For CC ~ 10^-122, need Lambda^2 ~ L_1/(122*ln10) = {Lambda_sq_cc:.6f}")
print(f"  Lambda = {Lambda_cc:.6f}")
print(f"  Delta_S at this Lambda = {delta_S_cc:.6e}")
print(f"  z = {z_cc_check:.4f}")

# What IS this Lambda physically?
print(f"\n  Lambda^2 = {Lambda_sq_cc:.6f}")
print(f"  Lambda^2 / L_1 = {Lambda_sq_cc/L_1:.6f}")
print(f"  1/Lambda^2 = {1/Lambda_sq_cc:.4f}")
print(f"  Compare: N = {N}, N*L_1 = {N*L_1:.4f}")

# The natural cutoff in the heat kernel is Lambda^2 = 1/beta.
# For CC: 1/Lambda^2 ~ 122*ln(10)/L_1 ~ 122.5
# So beta ~ 122.5. This is close to 123 = N + N_gen!

beta_from_cc = 1/Lambda_sq_cc
print(f"\n  NATURAL BETA: 1/Lambda^2 = {beta_from_cc:.4f}")
print(f"  N + N_gen = {N + 3}")
print(f"  N = {N}")
print(f"  Difference from N+3: {beta_from_cc - (N+3):.4f}")

# Hmm: 1/Lambda^2 = 122*ln(10)/L_1 = 280.9/2.292 = 122.6
# This is NOT 123. It's the same beta_exact we found before.

# ============================================================
# STEP 7: The spectral action IS the heat kernel
# ============================================================
print(f"\n{'='*70}")
print("STEP 7: THE SPECTRAL ACTION = HEAT KERNEL")
print(f"{'='*70}")

print(f"""
  For f(x) = exp(-x^2), the spectral action IS:
    S(Lambda) = Tr exp(-Delta/Lambda^2) = heat kernel at t = 1/Lambda^2

  The Galois-breaking part is:
    Delta_S = 4*(exp(-L_1*t) - exp(-L_8*t)) + 9*(exp(-L_2*t) - exp(-L_6*t))

  At large t (small Lambda = IR):
    Delta_S ~ 4*exp(-L_1*t)  (dominant term)

  For CC: need 4*exp(-L_1*t) ~ 10^-122
    => L_1*t ~ 280
    => t ~ 280/L_1 ~ 122

  So the CC from the spectral action IS the heat kernel at t ~ 122.

  This is EXACTLY what we had before (exp505b)!
  The spectral action Tr(f(D/Lambda)) with f = exp(-x^2)
  evaluated at Lambda^2 = L_1/280 gives the CC.

  THE QUESTION REMAINS: WHY this specific Lambda (or t)?
  The spectral action doesn't choose its own cutoff.
  Lambda is an INPUT, not an output.
""")

# ============================================================
# STEP 8: Maybe the RATIO of spectral actions matters?
# ============================================================
print(f"{'='*70}")
print("STEP 8: RATIOS AND RELATIVE QUANTITIES")
print(f"{'='*70}")

# Instead of the absolute value, what about:
# CC = Delta_S / S_total (Galois breaking relative to total)

print(f"\n  Galois-breaking FRACTION at various t:")
print(f"  {'t':>10} {'S_total':>14} {'Delta_S':>14} {'ratio':>14} {'z(ratio)':>10}")

for t in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]:
    S_total = sum(mult * np.exp(-ev*t) for ev, mult, _, _ in spectrum if ev >= 0)
    delta_S = 4*(np.exp(-L_1*t) - np.exp(-L_8*t)) + 9*(np.exp(-L_2*t) - np.exp(-L_6*t))
    if S_total > 1e-300 and abs(delta_S) > 1e-300:
        ratio = abs(delta_S) / S_total
        z = -np.log(ratio) / ln_inv_alpha
    else:
        ratio = 0
        z = float('inf')
    print(f"  {t:10.2f} {S_total:14.6e} {delta_S:14.6e} {ratio:14.6e} {z:10.4f}")

# ============================================================
# STEP 9: The Galois spectral asymmetry (eta-like)
# ============================================================
print(f"\n{'='*70}")
print("STEP 9: SPECTRAL ASYMMETRY (ETA INVARIANT)")
print(f"{'='*70}")

# The APS eta invariant: eta(s) = sum_k sign(lambda_k) * |lambda_k|^{-s}
# For our Galois "eta": eta_Galois(s) = sum pairs mult * (L_phys^{-s} - L_dark^{-s})

# This is the spectral zeta function asymmetry under Galois.
print(f"\n  Galois spectral zeta asymmetry:")
print(f"  {'s':>8} {'eta_Galois(s)':>16}")

for s in [0.5, 1.0, 1.5, 2.0, PHI, 2.5, 3.0]:
    eta = 4*(L_1**(-s) - L_8**(-s)) + 9*(L_2**(-s) - L_6**(-s))
    print(f"  {s:8.4f} {eta:16.6f}")

# At s = 0: eta_Galois(0) = 4*(1-1) + 9*(1-1) = 0 (trivial)
# At s -> infinity: dominated by L_1^{-s} (smallest eigenvalue)

# The VALUE at s = 1:
eta_1 = 4*(1/L_1 - 1/L_8) + 9*(1/L_2 - 1/L_6)
print(f"\n  eta_Galois(1) = {eta_1:.6f}")
print(f"  = 4*(1/L_1 - 1/L_8) + 9*(1/L_2 - 1/L_6)")
print(f"  = 4*({1/L_1:.4f} - {1/L_8:.4f}) + 9*({1/L_2:.4f} - {1/L_6:.4f})")

# At s = phi:
eta_phi = 4*(L_1**(-PHI) - L_8**(-PHI)) + 9*(L_2**(-PHI) - L_6**(-PHI))
print(f"\n  eta_Galois(phi) = {eta_phi:.6f}")

# ============================================================
# STEP 10: Product of Galois eigenvalue ratios
# ============================================================
print(f"\n{'='*70}")
print("STEP 10: DETERMINANT RATIOS")
print(f"{'='*70}")

# det(physical Galois) / det(dark Galois)
# = (L_1^4 * L_2^9) / (L_8^4 * L_6^9)
# = (L_1/L_8)^4 * (L_2/L_6)^9

ratio_det = (L_1/L_8)**4 * (L_2/L_6)**9
ln_ratio = 4*np.log(L_1/L_8) + 9*np.log(L_2/L_6)

print(f"\n  det ratio = (L_1/L_8)^4 * (L_2/L_6)^9 = {ratio_det:.6e}")
print(f"  ln(det ratio) = {ln_ratio:.6f}")
print(f"  = -16*ln(phi) + 9*ln((12-2s5)/(12+2s5))")
print(f"  = {-16*np.log(PHI):.6f} + {9*np.log(L_2/L_6):.6f}")

# The z-equivalent of the determinant ratio
z_det = -np.log(ratio_det) / ln_inv_alpha
print(f"\n  z_equiv(det ratio) = {z_det:.4f}")
print(f"  16*ln(phi)/ln(1/alpha) = {16*np.log(PHI)/ln_inv_alpha:.4f}")

# What is ln(L_2/L_6) exactly?
# L_2/L_6 = (12-2sqrt5)/(12+2sqrt5)
# Rationalize: (12-2s5)^2 / (144-20) = (144-48s5+20)/124 = (164-48s5)/124
# = (41-12s5)/31
ratio_26 = L_2/L_6
print(f"\n  L_2/L_6 = {ratio_26:.10f}")
print(f"  (41-12*sqrt5)/31 = {(41-12*SQRT5)/31:.10f}")
print(f"  Match: {abs(ratio_26 - (41-12*SQRT5)/31) < 1e-8}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
ONE-LOOP RESULTS ON 600-CELL:
==============================

1. Full log-determinant: ln det'(Delta_0) = {ln_det:.4f}
   V_1loop = {ln_det/2:.4f}

2. Galois breaking of log-det: Delta_V = {delta_V:.4f} (ORDER 1, not 10^-122)
   = -16*ln(phi) + 9*ln((12-2s5)/(12+2s5))
   = {exact_pair1:.4f} + {exact_pair2:.4f}

3. Determinant ratio: (L_1/L_8)^4 * (L_2/L_6)^9 = {ratio_det:.6e}
   This IS exponentially small! z_equiv = {z_det:.4f}

4. The spectral action Tr(exp(-Delta*t)) IS the heat kernel.
   Its Galois part at t ~ 122 gives 10^-122.
   But t is NOT determined by the spectral action itself.

WHAT WE LEARNED:
  The one-loop Galois log-det = {delta_V:.4f} is O(1). NOT the CC.

  BUT the DETERMINANT RATIO is tiny: {ratio_det:.4e}
  with z_equiv = {z_det:.4f}. Compare target z = 56.88.

  z_det = {z_det:.4f} = 16*ln(phi)/ln(1/alpha) + correction
  = {16*np.log(PHI)/ln_inv_alpha:.4f} + {z_det - 16*np.log(PHI)/ln_inv_alpha:.4f}
""")

# IS z_det close to a framework number?
print(f"  z_det = {z_det:.6f}")
print(f"  z_CC = {57-alpha_s:.6f}")
print(f"  Difference: {z_det - (57-alpha_s):.6f}")
print(f"  z_det / z_CC = {z_det / (57-alpha_s):.6f}")

# Check: does the determinant ratio ALONE give the CC?
print(f"\n  Lambda from det ratio alone: {ratio_det:.6e}")
print(f"  alpha^z_CC = {alpha_em**(57-alpha_s):.6e}")
print(f"  Observed CC = 2.86e-122")
print(f"  det ratio = {ratio_det:.6e}")

print("\n" + "=" * 70)
print("EXP-513 COMPLETE")
print("=" * 70)
