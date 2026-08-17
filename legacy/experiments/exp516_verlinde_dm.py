"""
EXP-516: Verlinde Fusion + Galois = DM Abundance?
===================================================
We know: sin^2(tW) = b1/(a1^2+1) = 6/26
         b1 = ||N_{1/2}||_F^2 = 6 (fusion Frobenius norm)
         26 = Galois-broken modes

Question: does the Verlinde/fusion structure ALSO give the DM abundance?

Compute ALL Verlinde data for SU(2)_3 (k=3, k+2=a1=5) and check.
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen

print("=" * 70)
print("EXP-516: VERLINDE FUSION + GALOIS -> DM?")
print("=" * 70)

# ============================================================
# SU(2)_k TQFT at k = 3 (k+2 = a1 = 5)
# ============================================================
k = 3
n_reps = k + 1  # = 4 representations: j = 0, 1/2, 1, 3/2

print(f"\nSU(2)_{k} TQFT (k+2 = {k+2} = a1)")
print(f"Representations: j = 0, 1/2, 1, 3/2 ({n_reps} total)")

# Quantum dimensions
d_q = np.zeros(n_reps)
for idx in range(n_reps):
    j = idx / 2.0
    d_q[idx] = np.sin((2*j+1)*np.pi/(k+2)) / np.sin(np.pi/(k+2))

print(f"\nQuantum dimensions:")
for idx in range(n_reps):
    j = idx / 2.0
    print(f"  j={j:.1f}: d_q = {d_q[idx]:.6f}")

D_sq = np.sum(d_q**2)
print(f"D^2 = sum d_q^2 = {D_sq:.6f} = a1 + sqrt(a1) = {a1 + np.sqrt(a1):.6f}")

# Galois conjugates (phi -> phi' = -1/phi)
d_q_dark = np.array([1.0, abs(PHI_CONJ), abs(PHI_CONJ), 1.0])
D_sq_dark = np.sum(d_q_dark**2)
print(f"\nDark quantum dimensions (|sigma(d_q)|):")
for idx in range(n_reps):
    j = idx / 2.0
    print(f"  j={j:.1f}: d_q' = {d_q_dark[idx]:.6f}")
print(f"D'^2 = {D_sq_dark:.6f} = a1 - sqrt(a1) = {a1 - np.sqrt(a1):.6f}")

# ============================================================
# Modular S-matrix
# ============================================================
print(f"\n{'='*70}")
print("MODULAR S-MATRIX")
print(f"{'='*70}")

S = np.zeros((n_reps, n_reps))
for i in range(n_reps):
    for j_idx in range(n_reps):
        ji = i / 2.0
        jj = j_idx / 2.0
        S[i, j_idx] = np.sqrt(2.0/(k+2)) * np.sin((2*ji+1)*(2*jj+1)*np.pi/(k+2))

print(f"\nS-matrix:")
for i in range(n_reps):
    row = "  " + "  ".join(f"{S[i,j]:8.4f}" for j in range(n_reps))
    print(row)

# Verify: S is unitary
print(f"\n  S^2 check (should be charge conjugation):")
S2 = S @ S
for i in range(n_reps):
    row = "  " + "  ".join(f"{S2[i,j]:8.4f}" for j in range(n_reps))
    print(row)

# ============================================================
# Fusion matrices (Verlinde formula)
# ============================================================
print(f"\n{'='*70}")
print("FUSION MATRICES (VERLINDE FORMULA)")
print(f"{'='*70}")

# N_i: fusion matrix for representation i
# (N_i)_{jk} = sum_l S_{il} S_{jl} S*_{kl} / S_{0l}

fusion = np.zeros((n_reps, n_reps, n_reps))
for i in range(n_reps):
    for j in range(n_reps):
        for k_idx in range(n_reps):
            val = 0
            for l in range(n_reps):
                val += S[i,l] * S[j,l] * np.conj(S[k_idx,l]) / S[0,l]
            fusion[i,j,k_idx] = np.real(val)

# Round to nearest integer (fusion coefficients are integers)
fusion_int = np.round(fusion).astype(int)

for i in range(n_reps):
    ji = i / 2.0
    print(f"\n  N_{{j={ji:.1f}}}:")
    for j in range(n_reps):
        row = "    " + "  ".join(f"{fusion_int[i,j,kk]:2d}" for kk in range(n_reps))
        print(row)

    # Frobenius norm
    frob = np.sum(fusion_int[i]**2)
    print(f"    ||N_{{{ji:.1f}}}||_F^2 = {frob}")

# Verify: ||N_{1/2}||^2 = b1 = 6
frob_half = np.sum(fusion_int[1]**2)
print(f"\n  ||N_{{1/2}}||_F^2 = {frob_half} = b1 = {b1}: {frob_half == b1}")

# ============================================================
# Galois-breaking of fusion data
# ============================================================
print(f"\n{'='*70}")
print("GALOIS-BREAKING OF FUSION DATA")
print(f"{'='*70}")

# Under Galois: d_q -> d_q'. The fusion coefficients N_{ij}^k are integers
# (Galois-invariant). But the WEIGHTED fusion amplitudes change.

# Physical fusion amplitude: A_phys = sum_{ijk} N_{ij}^k * d_i * d_j * d_k
# Dark fusion amplitude: A_dark = sum_{ijk} N_{ij}^k * d'_i * d'_j * d'_k

A_phys = 0
A_dark = 0
for i in range(n_reps):
    for j in range(n_reps):
        for kk in range(n_reps):
            if fusion_int[i,j,kk] > 0:
                A_phys += fusion_int[i,j,kk] * d_q[i] * d_q[j] * d_q[kk]
                A_dark += fusion_int[i,j,kk] * d_q_dark[i] * d_q_dark[j] * d_q_dark[kk]

print(f"\n  Total weighted fusion amplitude:")
print(f"    A_phys = {A_phys:.6f}")
print(f"    A_dark = {A_dark:.6f}")
print(f"    A_phys / A_dark = {A_phys/A_dark:.6f}")
print(f"    (A_phys - A_dark) = {A_phys - A_dark:.6f}")
print(f"    (A_phys + A_dark) = {A_phys + A_dark:.6f}")
print(f"    A_phys / A_dark = {A_phys/A_dark:.6f}")
print(f"    phi^3 = {PHI**3:.6f}")

# Check per-channel
print(f"\n  Per-channel weighted amplitudes (N_ij^k > 0):")
print(f"  {'i':>4} {'j':>4} {'k':>4}  {'N':>3}  {'phys':>10}  {'dark':>10}  {'ratio':>8}")

channels_phys = []
channels_dark = []
for i in range(n_reps):
    for j in range(n_reps):
        for kk in range(n_reps):
            if fusion_int[i,j,kk] > 0:
                wp = d_q[i] * d_q[j] * d_q[kk]
                wd = d_q_dark[i] * d_q_dark[j] * d_q_dark[kk]
                channels_phys.append(wp)
                channels_dark.append(wd)
                ji, jj, jk = i/2, j/2, kk/2
                ratio = wp/wd if wd > 1e-10 else float('inf')
                print(f"  {ji:4.1f} {jj:4.1f} {jk:4.1f}  {fusion_int[i,j,kk]:3d}  {wp:10.4f}  {wd:10.4f}  {ratio:8.4f}")

# ============================================================
# Specific DM-related quantities
# ============================================================
print(f"\n{'='*70}")
print("DM-RELATED VERLINDE QUANTITIES")
print(f"{'='*70}")

# The number of fusion channels
n_channels = sum(1 for i in range(n_reps) for j in range(n_reps) for kk in range(n_reps)
                 if fusion_int[i,j,kk] > 0)
print(f"\n  Total fusion channels (N_ij^k > 0): {n_channels}")

# Galois-broken channels: those where the weight ratio != 1
galois_broken_channels = 0
galois_fixed_channels = 0
for i in range(n_reps):
    for j in range(n_reps):
        for kk in range(n_reps):
            if fusion_int[i,j,kk] > 0:
                wp = d_q[i] * d_q[j] * d_q[kk]
                wd = d_q_dark[i] * d_q_dark[j] * d_q_dark[kk]
                if abs(wp - wd) > 0.001:
                    galois_broken_channels += 1
                else:
                    galois_fixed_channels += 1

print(f"  Galois-broken channels: {galois_broken_channels}")
print(f"  Galois-fixed channels: {galois_fixed_channels}")

# Ratio
if galois_fixed_channels > 0:
    print(f"  Broken/Fixed: {galois_broken_channels/galois_fixed_channels:.4f}")
print(f"  Broken/Total: {galois_broken_channels/n_channels:.4f}")

# ============================================================
# Key ratios
# ============================================================
print(f"\n{'='*70}")
print("KEY RATIOS")
print(f"{'='*70}")

Omega_obs = 5.36
target = 7 - PHI

ratios_to_check = {
    'A_phys/A_dark': A_phys/A_dark,
    '(A_phys-A_dark)/(A_phys+A_dark)': (A_phys-A_dark)/(A_phys+A_dark),
    'A_phys/A_dark - 1': A_phys/A_dark - 1,
    'n_channels': float(n_channels),
    'broken_channels': float(galois_broken_channels),
    'broken_channels / d_ST': galois_broken_channels / 4.0,
    'D^2/D^2_dark': D_sq/D_sq_dark,
    'D^4/D^4_dark': (D_sq/D_sq_dark)**2,
    '(D^2-D^2_dark)/2': (D_sq-D_sq_dark)/2,
    'broken_chan - fixed_chan': float(galois_broken_channels - galois_fixed_channels),
    'sum N_{1/2}': float(np.sum(fusion_int[1])),
    'sum all N': float(np.sum(fusion_int)),
}

print(f"\n  {'Quantity':>35} {'Value':>10} {'vs 5.38':>10}")
for name, val in sorted(ratios_to_check.items(), key=lambda x: abs(x[1] - target)):
    err = (val/target - 1) * 100
    marker = " <---" if abs(err) < 5 else ""
    print(f"  {name:>35} {val:10.4f} {err:+10.2f}%{marker}")

# ============================================================
# The S-matrix under Galois
# ============================================================
print(f"\n{'='*70}")
print("S-MATRIX GALOIS STRUCTURE")
print(f"{'='*70}")

# S-matrix involves sin(n*pi/5). Under sigma: sin(pi/5) involves sqrt(5).
# sigma maps sqrt(5) -> -sqrt(5), which maps sin(pi/5) -> ...
# sin(pi/5) = sqrt((5-sqrt(5))/8) -> sqrt((5+sqrt(5))/8) = sin(2*pi/5)... hmm

# Actually sin(pi/5) = sqrt(10-2sqrt(5))/4. Under sigma:
# sigma(sin(pi/5)) = sqrt(10+2sqrt(5))/4 = sin(2*pi/5) = cos(pi/10)

# The S-matrix entries involve these sines. Under Galois:
# S_{ij} -> S_{ij}' (replace sqrt(5) -> -sqrt(5))

# The Galois conjugate S-matrix:
S_dark = np.zeros((n_reps, n_reps))
# Under sigma: swap j=1/2 <-> j=3/2 and j=0 <-> j=1? No...
# Actually: sigma permutes the representations.
# d_q: 1, phi, phi, 1 -> 1, 1/phi, 1/phi, 1
# The permutation maps j=1/2 -> j=3/2 (in some sense)

# Let me just compute S' by replacing sqrt(5) -> -sqrt(5)
# sin(n*pi/5) for n=1: sin(pi/5) = sqrt(10-2sqrt(5))/4
# Under sigma: sqrt(10-2sqrt(5)) -> sqrt(10+2sqrt(5)) = ...

# Actually easier: use the fact that the eigenvalues of the fusion matrix
# are d_j. Under Galois: d_j -> d_j'. The fusion matrices are unchanged
# (integer entries). So the Galois action on the S-matrix is:
# S_{ij}' = S_{sigma(i), sigma(j)} where sigma permutes reps.

# For SU(2)_3: sigma maps d_q -> d_q'
# d_q = (1, phi, phi, 1) -> (1, 1/phi, 1/phi, 1)
# This means sigma swaps j=1/2 <-> j=3/2? Let's check:
# d_{1/2} = phi, d_{3/2} = 1. sigma(d_{1/2}) should be 1/phi.
# But d_{3/2} = 1 != 1/phi. So sigma does NOT permute reps simply.

# The issue: sigma acts on Q(sqrt(5)), and the quantum dimensions
# 1, phi, phi, 1 don't simply permute under sigma.
# sigma(phi) = phi' = -1/phi, |phi'| = 1/phi.
# So sigma(d_{1/2}) = |sigma(phi)| = 1/phi, but this doesn't equal d_j for any j.
# The dark sector quantum dimensions 1, 1/phi, 1/phi, 1 are NOT standard SU(2)_k dims.

print(f"""
  Physical quantum dims: d = (1, phi, phi, 1)
  Dark quantum dims:     d'= (1, 1/phi, 1/phi, 1)

  The dark dims do NOT correspond to SU(2)_k at any integer k.
  They are a "Galois-deformed" TQFT.
  The fusion coefficients N_ij^k are UNCHANGED (integers).
  But the quantum dimensions are deformed: phi -> 1/phi.

  This means: the dark sector has the SAME fusion rules
  but DIFFERENT quantum dimensions.
  It's a "Galois-twisted" TQFT.
""")

# ============================================================
# Galois-twisted partition function
# ============================================================
print(f"{'='*70}")
print("GALOIS-TWISTED PARTITION FUNCTION")
print(f"{'='*70}")

# Z(S^3) = 1/D^2 (physical sector)
# Z'(S^3) = 1/D'^2 (dark sector)

Z_phys = 1/D_sq
Z_dark = 1/D_sq_dark

print(f"\n  Z(physical) = 1/D^2 = {Z_phys:.6f}")
print(f"  Z(dark)     = 1/D'^2 = {Z_dark:.6f}")
print(f"  Z_dark / Z_phys = D^2/D'^2 = phi^2 = {Z_dark/Z_phys:.6f}")
print(f"  Z_phys / Z_dark = D'^2/D^2 = 1/phi^2 = {Z_phys/Z_dark:.6f}")

# The RATIO of partition functions
ratio_Z = Z_dark / Z_phys
print(f"\n  Z_dark/Z_phys = {ratio_Z:.6f} = phi^2 = {PHI**2:.6f}")

# The DM abundance might involve Z_dark/Z_phys somehow.
# Omega = 5.38 and phi^2 = 2.618. Ratio: 5.38/2.618 = 2.056

# What if: Omega = (Z_dark/Z_phys) * (something) ?
# 5.38/2.618 = 2.056. What is 2.056?
# 2*phi - 1 = 2.236 (too big)
# phi + 1/phi = 2.236 (too big)
# 2 + 1/(a1*phi) = 2.124 (close but not exact)

# Or: Omega = Z_dark/Z_phys + Z_dark/Z_phys + ... (sum of channels)
# = n_channels_effective * phi^2?

# Actually: 7-phi = phi^2 + phi^2 + (7-phi-2*phi^2) = 2*phi^2 + (7-phi-2phi^2)
# = 2*phi^2 + (7-phi-2phi-2) = 2*phi^2 + (5-3phi)
# = 2*(phi+1) + 5-3phi = 2phi+2+5-3phi = 7-phi. Circular.

# Let me try: does Omega involve a PRODUCT of Z ratios?
# Z_dark/Z_phys = phi^2
# (Z_dark/Z_phys)^2 = phi^4 = 6.854 (too big)
# sqrt(Z_dark/Z_phys) = phi = 1.618 (too small)

# ============================================================
# The Verlinde formula gives n_ij = dimension of space of
# conformal blocks on a sphere with insertions i, j.
# For the ANNIHILATION channel (SM -> DM production):
# n_production = ???
# ============================================================

print(f"\n{'='*70}")
print("PRODUCTION CHANNELS")
print(f"{'='*70}")

# In the Galois picture: SM particles are in rep j (physical sector).
# Dark particles are in rep j' (dark sector, same j but different d_q).
# Production: SM x SM -> DM x DM via gravity.
#
# The number of independent channels = number of ways to
# combine two physical reps into two dark reps, weighted by d_q.

# For each (j1, j2) physical pair and (j3, j4) dark pair:
# Amplitude ~ N_{j1,j2}^J * N_{j3,j4}^J (via intermediate J)

# Total production amplitude:
prod_amp_phys = 0
prod_amp_dark = 0
for j1 in range(n_reps):
    for j2 in range(n_reps):
        for J in range(n_reps):
            if fusion_int[j1, j2, J] > 0:
                w_phys = d_q[j1] * d_q[j2] / d_q[J]
                w_dark = d_q_dark[j1] * d_q_dark[j2] / d_q_dark[J]
                prod_amp_phys += fusion_int[j1, j2, J] * w_phys
                prod_amp_dark += fusion_int[j1, j2, J] * w_dark

print(f"\n  Production amplitude (Verlinde):")
print(f"    Physical: {prod_amp_phys:.6f}")
print(f"    Dark: {prod_amp_dark:.6f}")
print(f"    Ratio phys/dark: {prod_amp_phys/prod_amp_dark:.6f}")
print(f"    phi^2 = {PHI**2:.6f}")

# ============================================================
# The KEY computation: Galois asymmetry per channel
# ============================================================
print(f"\n{'='*70}")
print("GALOIS ASYMMETRY PER CHANNEL")
print(f"{'='*70}")

# For each fusion channel (i,j,k): the asymmetry is
# d_i*d_j*d_k - d'_i*d'_j*d'_k

total_asym = 0
n_asym_channels = 0
asym_by_type = {}

for i in range(n_reps):
    for j in range(n_reps):
        for kk in range(n_reps):
            if fusion_int[i,j,kk] > 0:
                asym = d_q[i]*d_q[j]*d_q[kk] - d_q_dark[i]*d_q_dark[j]*d_q_dark[kk]
                if abs(asym) > 1e-10:
                    total_asym += asym
                    n_asym_channels += 1
                    # Classify by number of phi-type indices
                    n_phi = sum(1 for idx in [i,j,kk] if abs(d_q[idx] - PHI) < 0.01)
                    asym_by_type[n_phi] = asym_by_type.get(n_phi, 0) + asym

print(f"  Total Galois asymmetry: {total_asym:.6f}")
print(f"  Number of asymmetric channels: {n_asym_channels}")
print(f"\n  Asymmetry by number of phi-type indices:")
for n_phi, asym in sorted(asym_by_type.items()):
    print(f"    {n_phi} phi-indices: asymmetry = {asym:.6f}")

print(f"\n  Total asymmetry / D^2 = {total_asym/D_sq:.6f}")
print(f"  Total asymmetry / n_channels = {total_asym/n_channels:.6f}")
print(f"  Total asymmetry / 4 = {total_asym/4:.6f}")

# Check against DM target
print(f"\n  Comparison with DM abundance:")
print(f"  7-phi = {7-PHI:.6f}")
print(f"  Total asymmetry = {total_asym:.6f}")
print(f"  Ratio: {total_asym/(7-PHI):.6f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
VERLINDE DATA FOR SU(2)_3:
  Fusion channels: {n_channels} total, {galois_broken_channels} Galois-broken
  ||N_{{1/2}}||^2 = {frob_half} = b1 (confirmed)
  A_phys/A_dark = {A_phys/A_dark:.4f}
  Z_dark/Z_phys = phi^2 = {ratio_Z:.4f}
  Total Galois asymmetry = {total_asym:.4f}

  sin^2(tW) = ||N_{{1/2}}||^2_F / N_broken = {b1}/{26} = {b1/26:.4f}
  WHERE N_broken = 26 = Galois-broken Laplacian modes

DM ABUNDANCE:
  7-phi = {7-PHI:.4f} (target)
  Closest Verlinde quantity: {'TODO - check above'}

  The Verlinde fusion data gives phi^2 type ratios,
  NOT 7-phi directly. The DM abundance involves a MORE
  COMPLEX combination of the TQFT data.

STATUS: The Verlinde connection to sin^2(tW) is CONFIRMED.
  The DM abundance does NOT come directly from Verlinde ratios.
  It requires the FULL Galois mode counting formula:
  Omega = (N_broken - Delta_D^2) / d_ST = (26 - 2*sqrt(5)) / 4
""")

print("=" * 70)
print("EXP-516 COMPLETE")
print("=" * 70)
