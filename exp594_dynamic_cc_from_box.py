#!/usr/bin/env python3
"""
EXP-594: CONSTANTA COSMOLOGICA DINAMICA DIN MODURI BOX
=======================================================

IDEEA: Box are 9 eigenvalues discrete. Pe masura ce universul se expandeaza,
modurile "ingheata" (devin non-relativiste) unul cate unul.
Fiecare freeze-out schimba vacuum energy cu un pas discret.

DESI (2025, 4.2 sigma): w(z) traverseaza -1 la z ~ 0.5.
Dark energy se SLABESTE. Poate un mod Box a inghetat recent?

PLAN:
1. Calculam contributia fiecarui mod Box la vacuum energy
2. Vedem CAND (la ce z) fiecare mod ingheata
3. Verificam daca freeze-out-ul corespunde cu z ~ 0.5
4. Calculam w(z) efectiv
"""
import numpy as np
import sys
sys.path.insert(0, ".")
from commons.constants import PHI, SQRT5, a1, b1, N, N_gen, alpha_em

print("=" * 80)
print("EXP-594: CC DINAMIC DIN MODURI BOX DISCRETE")
print("=" * 80)

# =====================================================================
# PART 1: Box eigenvalues and their physical masses
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: EIGENVALUES BOX SI MASELE FIZICE ASOCIATE")
print("=" * 80)

# 9 distinct Box eigenvalue magnitudes (from exp591)
# Associated with 9 irreps of 2I
# Format: (irrep, dim, mult, |eigenvalue|, exact_form)
box_modes = [
    (0, 1,  1,   0.0,           "0 (kernel)"),
    (1, 2,  4,   0.0,           "0 (kernel)"),
    (2, 3,  9,   5.5279,        "a1-sqrt5 and 2*(a1-sqrt5)"),
    (3, 4, 16,   6.7082,        "3*sqrt5"),
    (4, 5, 25,   7.7666,        "mixed: b1*phi, b1/phi, 12"),
    (5, 6, 36,   7.8055,        "mixed: 10, ~1.71, ~11.71"),
    (6, 3,  9,   9.6481,        "a1+sqrt5 and 2*(a1+sqrt5)"),
    (7, 4, 16,   6.7082,        "3*sqrt5 (Galois of rho_3)"),
    (8, 2,  4,   0.0,           "0 (kernel)"),
]

print(f"\n  {'rho':>4s} {'dim':>4s} {'mult':>5s} {'<|mu|>':>8s} {'form'}")
print(f"  {'-'*50}")
for rho, d, mult, mu, form in box_modes:
    print(f"  rho_{rho}  {d:>3d}  {mult:>4d}  {mu:>8.4f}  {form}")

# Physical mass of each mode: m_mode = m_Planck * |mu| / Lambda_cutoff
# In Planck units, the Box eigenvalue IS the mass squared (up to normalization)
# The key: each mode contributes to vacuum energy proportional to m^4 (or m^2 * Lambda^2)

# In the spectral action framework:
# Lambda_CC = sum over modes of f(m_mode / Lambda_cutoff)
# When m_mode >> T (temperature), the mode is "frozen" and doesn't contribute to pressure
# When m_mode << T, the mode is relativistic and contributes

# =====================================================================
# PART 2: Spectral action vacuum energy decomposition
# =====================================================================
print("\n" + "=" * 80)
print("PART 2: SPECTRAL ACTION — VACUUM ENERGY PER MOD")
print("=" * 80)

# Seeley-DeWitt coefficients from paper:
# Tr(Box^0) = 120 = c_cosmological (proportional to CC)
# Tr(Box^2) = 7200 = c_EH (proportional to Einstein-Hilbert)
# Actual Seeley-DeWitt: c0=2640, c1=14880, c2=55920

# But for vacuum energy, the relevant quantity is the ZERO-POINT energy:
# E_vac = (1/2) * sum_modes omega_k
# For a discrete system: E_vac = (1/2) * sum |lambda_k|

# Each irrep rho_k contributes d_k^2 modes (multiplicity) with eigenvalue |mu_k|
# Total zero-point: E_vac = (1/2) * sum_k d_k^2 * <|mu_k|>

E_vac_total = 0.5 * sum(mult * mu for _, _, mult, mu, _ in box_modes)
print(f"\n  Zero-point energy: E_vac = (1/2) * sum d^2 * <|mu|>")
print(f"  = {E_vac_total:.4f} (in lattice units)")

print(f"\n  Contributie per mod:")
for rho, d, mult, mu, form in box_modes:
    contrib = 0.5 * mult * mu
    frac = contrib / E_vac_total * 100 if E_vac_total > 0 else 0
    print(f"    rho_{rho} (d={d}): {contrib:>8.2f} ({frac:>5.1f}%)")

# =====================================================================
# PART 3: Freeze-out redshifts
# =====================================================================
print("\n" + "=" * 80)
print("PART 3: FREEZE-OUT REDSHIFTS")
print("=" * 80)

# A mode "freezes out" when the cosmic temperature drops below its mass
# T(z) = T_0 * (1+z) where T_0 = 2.725 K = 2.35e-4 eV
# Mode mass (physical): m_phys = m_Planck * |mu| * (some scale factor)

# The question is: what sets the physical mass scale of Box modes?
# In the framework: m_Z = m_e * phi^25 sets the EW scale
# Box eigenvalues are dimensionless. Their physical mass depends on the cutoff.

# Approach 1: modes are at Planck scale -> ALL frozen, no dynamics
# Approach 2: modes are at EW scale -> some might be cosmologically relevant
# Approach 3: modes are at a DARK SECTOR scale

# Let's think about this differently.
# The CC in the framework comes from Lambda_P = alpha^(57-alpha_s)
# In Planck units: Lambda ~ 10^{-122}
# This is the OBSERVED CC, incredibly small.

# What if the CC is the RESIDUAL after cancellation between modes?
# Massive modes contribute +E_vac, and the near-cancellation leaves a tiny remnant.

print(f"""
  PROBLEMA DE SCALA:

  Box eigenvalues sunt O(1-10) in unitati de lattice (Planck).
  Vacuum energy E_vac ~ sum |mu| ~ 500 (Planck units).
  Observat: Lambda ~ 10^{{-122}} Planck units.

  Gap de 122 ordine de marime!

  Aceasta e PROBLEMA COSMOLOGICA STANDARD — nu e specifica framework-ului.
  Niciun framework discret nu rezolva asta prin simpla sumare de moduri.
""")

# =====================================================================
# PART 4: Alternative — CC from Galois cancellation
# =====================================================================
print("=" * 80)
print("PART 4: CC DIN ANULARE GALOIS")
print("=" * 80)

# Key idea: physical and dark (Galois conjugate) sectors contribute
# with OPPOSITE signs to vacuum energy.
# The CC is the DIFFERENCE between physical and dark contributions.

# Galois pairs:
# rho_1 (phys) <-> rho_8 (dark): both have |mu|=0
# rho_2 (phys) <-> rho_6 (dark): |mu|=3.685 vs 9.648
# rho_3 (phys) <-> rho_7 (dark): both |mu|=6.708
# Fixed: rho_0 (|mu|=0), rho_4 (|mu|=7.767), rho_5 (|mu|=7.806)

print(f"\n  Galois pairs si vacuum energy:")
print(f"  {'Pair':>12s} {'E_phys':>10s} {'E_dark':>10s} {'Delta':>10s}")
print(f"  {'-'*47}")

galois_pairs = [
    ("rho1-rho8", 4, 0.0, 4, 0.0),       # kernel pair
    ("rho2-rho6", 9, 3.6852, 9, 9.6481),  # broken pair
    ("rho3-rho7", 16, 6.7082, 16, 6.7082),  # broken pair (SAME |mu|!)
]
galois_fixed = [
    ("rho0", 1, 0.0),
    ("rho4", 25, 7.7666),
    ("rho5", 36, 7.8055),
]

total_phys = 0
total_dark = 0
for name, mult_p, mu_p, mult_d, mu_d in galois_pairs:
    E_p = 0.5 * mult_p * mu_p
    E_d = 0.5 * mult_d * mu_d
    delta = E_p - E_d
    total_phys += E_p
    total_dark += E_d
    print(f"  {name:>12s} {E_p:>10.4f} {E_d:>10.4f} {delta:>+10.4f}")

print(f"\n  Galois-fixed modes (contribute to both sectors equally):")
for name, mult, mu in galois_fixed:
    E = 0.5 * mult * mu
    print(f"    {name}: E = {E:.4f}")

print(f"\n  Total physical: {total_phys:.4f}")
print(f"  Total dark:     {total_dark:.4f}")
print(f"  Difference:     {total_phys - total_dark:+.4f}")
print(f"  Ratio |Delta/Total|: {abs(total_phys-total_dark)/(total_phys+total_dark)*100:.2f}%")

# The key: rho3 and rho7 have EXACTLY the same |mu| -> perfect cancellation!
# The residual comes entirely from the rho2-rho6 pair
delta_26 = 0.5 * 9 * (3.6852 - 9.6481)
print(f"\n  OBSERVATIE: rho3-rho7 se ANULEAZA EXACT (acelasi |mu|=3*sqrt5)")
print(f"  Residuul vine DOAR din rho2-rho6: Delta = {delta_26:.4f}")
print(f"  rho2: <|mu|> = a1-sqrt5 = {a1-SQRT5:.4f}")
print(f"  rho6: <|mu|> = a1+sqrt5 = {a1+SQRT5:.4f}")
print(f"  Diferenta |mu|: 2*sqrt5 = {2*SQRT5:.4f}")

# =====================================================================
# PART 5: Dynamic CC from Galois asymmetry evolution
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 5: CC DINAMIC DIN EVOLUTIA ASIMETRIEI GALOIS")
print("=" * 80)

# Idea: if the Galois asymmetry (difference between physical and dark
# sector contributions) depends on the cosmic scale, then CC is dynamic.

# The asymmetry comes from ONE pair: rho_2 vs rho_6
# rho_2 has eigenvalue L = 10-2*sqrt5 (physical Galois)
# rho_6 has eigenvalue L = 10+2*sqrt5 (dark Galois)

# Their Box eigenvalues:
# rho_2: mu = {-(a1-sqrt5), 2*(a1-sqrt5)} with mean |mu| = (a1-sqrt5)*(2+1)/3... no
# Actually from exp593:
# rho_2: distinct eigenvalues -2.764, 5.528
#   = -(a1-sqrt5), 2*(a1-sqrt5)
#   multiplicities: 6, 3 (dim^2 = 9, split 2:1 by Galois)
# rho_6: distinct eigenvalues -7.236, 14.472
#   = -(a1+sqrt5), 2*(a1+sqrt5)
#   multiplicities: 6, 3

# The vacuum energy per pair depends on which eigenvalues are "active"
# If only the POSITIVE eigenvalue contributes (modes with mu > 0 are "particles"):
E_rho2_pos = 0.5 * 3 * 2*(a1-SQRT5)  # 3 modes with mu = 2*(a1-sqrt5) = 5.528
E_rho6_pos = 0.5 * 3 * 2*(a1+SQRT5)  # 3 modes with mu = 2*(a1+sqrt5) = 14.472
E_rho2_neg = 0.5 * 6 * abs(-(a1-SQRT5))  # 6 modes with |mu| = a1-sqrt5 = 2.764
E_rho6_neg = 0.5 * 6 * abs(-(a1+SQRT5))  # 6 modes with |mu| = a1+sqrt5 = 7.236

print(f"\n  rho_2 (physical Galois):")
print(f"    Positive modes (3): mu = 2*(a1-sqrt5) = {2*(a1-SQRT5):.4f}, E = {E_rho2_pos:.4f}")
print(f"    Negative modes (6): |mu| = a1-sqrt5 = {a1-SQRT5:.4f}, E = {E_rho2_neg:.4f}")
print(f"  rho_6 (dark Galois):")
print(f"    Positive modes (3): mu = 2*(a1+sqrt5) = {2*(a1+SQRT5):.4f}, E = {E_rho6_pos:.4f}")
print(f"    Negative modes (6): |mu| = a1+sqrt5 = {a1+SQRT5:.4f}, E = {E_rho6_neg:.4f}")

# The Galois difference in vacuum energy:
Delta_pos = E_rho2_pos - E_rho6_pos  # from positive eigenvalue modes
Delta_neg = E_rho2_neg - E_rho6_neg  # from negative eigenvalue modes
Delta_total = Delta_pos + Delta_neg

print(f"\n  Delta (pos modes): {Delta_pos:+.4f}")
print(f"  Delta (neg modes): {Delta_neg:+.4f}")
print(f"  Delta (total):     {Delta_total:+.4f}")

# The ratio of Galois asymmetry to total
# This could be related to Lambda in Planck units
ratio = abs(Delta_total) / E_vac_total
print(f"\n  |Delta|/E_total = {ratio:.6f}")
print(f"  log10(ratio) = {np.log10(ratio):.2f}")

# =====================================================================
# PART 6: The spectral determinant approach
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 6: EXPONENTUL CC DIN SPECTRAL DATA")
print("=" * 80)

# From the paper: Lambda_P = alpha^(57-alpha_s) where 57 = (N-2*N_gen)/2
# The exponent 57 is derived 7 ways. But the BASE alpha^z is pattern.

# Can we get the exponent from Box spectral data?
# log(Lambda)/log(alpha) should be ~57 - alpha_s ~ 56.88

# Spectral determinant: log det'(L) / ln(alpha)
# This was tried before (near-miss: 56.97 vs 57)

# What about log det'(Box)?
# Box eigenvalues (nonzero, with multiplicities):
nonzero_box_evals = []
for rho, d, mult, mu_mean, form in box_modes:
    if mu_mean > 0.01:
        nonzero_box_evals.extend([mu_mean] * mult)

log_det_box = sum(np.log(abs(e)) for e in nonzero_box_evals if abs(e) > 0.01)
print(f"\n  log det'(|Box|) = sum log|mu_k| = {log_det_box:.6f}")
print(f"  Divide by ln(alpha): {log_det_box / np.log(alpha_em):.4f}")
print(f"  Target: 57 - alpha_s = {57 - 1/(2*PHI**3):.4f}")

# Try with the EXACT eigenvalues from exp593
exact_evals_with_mult = [
    # rho_2: 6 modes at a1-sqrt5, 3 modes at 2*(a1-sqrt5)
    *([a1-SQRT5]*6), *([2*(a1-SQRT5)]*3),
    # rho_3: 8 modes at 3*sqrt5 (both signs, but |mu|=3*sqrt5)
    *([3*SQRT5]*16),
    # rho_4: mixed
    *([6*PHI]*10), *([6/PHI]*10), *([12]*5),
    # rho_5: mixed (need exact values)
    *([10]*12), *([abs(12-6*PHI)]*12), *([abs(12-6/PHI)]*12),
    # rho_6: 6 modes at a1+sqrt5, 3 modes at 2*(a1+sqrt5)
    *([a1+SQRT5]*6), *([2*(a1+SQRT5)]*3),
    # rho_7: 8+8 modes at 3*sqrt5
    *([3*SQRT5]*16),
]

# Actually, let me just use the mean values with multiplicities
# This is approximate but gives the right ballpark
print(f"\n  Nota: calculul exact necesita eigenvalues individuale, nu medii.")
print(f"  Rezultatul de mai sus e APROXIMATIV.")

# =====================================================================
# PART 7: Key question — can DESI w(z) come from this?
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 7: POATE DESI w(z) SA VINA DIN FRAMEWORK?")
print("=" * 80)

print(f"""
  DESI DR2 (2025): w(z) traverseaza -1 la z ~ 0.5
  w0 ~ -0.8, wa ~ -0.7 (CPL parametrization)

  Ce ar trebui framework-ul sa produca:
  1. CC care scade in timp (dark energy se slabeste)
  2. Traversare w = -1 la z ~ 0.5 (acum ~5 Gyr)
  3. w < -1 in trecut (phantom), w > -1 acum (quintessence)

  CE POATE FRAMEWORK-UL:

  A) Anulare Galois partiala:
     rho_2 (phys) vs rho_6 (dark) da un reziduu.
     Daca raportul phys/dark depinde de scala -> CC dinamic.
     DAR: de ce ar depinde de scala? Eigenvalues sunt FIXE.

  B) Spectral action cu cutoff variabil:
     S = c0*f4*Lambda^4 + c1*f2*Lambda^2 + ...
     Daca Lambda(t) scade cu expansiunea -> CC scade.
     NATURAL: Lambda = 1/(scale factor) in unele modele.
     DA Lambda^4 scade prea REPEDE (ca a^{{-4}}).

  C) Modul SPECIFIC care ingheata:
     Daca un mod Box devine non-relativist la z~0.5,
     contributia lui la presiune se schimba (p -> 0).
     Asta schimba w de la -1 spre -1+delta.
     PROBLEMA: modurile Box sunt la scala PLANCK.
     Un mod Planckian "ingheata" la T ~ M_Planck,
     adica in primele 10^{{-43}} secunde. Nu la z=0.5.

  D) Sectorul dark cu masa mica:
     Daca particulele dark Galois au mase in eV range,
     ele ar putea deveni non-relativiste la z ~ 0.5.
     T(z=0.5) ~ 3.7 K ~ 3e-4 eV.
     Masa dark particle ~ eV -> freeze-out la z ~ eV/(2.7K*k_B) ~ 4000.
     Prea devreme. Masa ~ meV -> z ~ 1. Mai aproape!

     Neutrino m3 = 49.5 meV -> freeze-out la z ~ 200. Prea devreme.

     Exista un mod cu masa ~ 0.1 meV?
     Asta ar da freeze-out la z ~ 0.5!

  E) RUNNING al exponentului CC:
     Lambda(z) = alpha(mu(z))^(57-alpha_s(mu(z)))
     alpha si alpha_s ruleaza cu energia.
     La z=0.5: mu ~ few eV.
     Running de la mu=0 la mu=few eV e NEGLIJABIL.

  CONCLUZIE PRELIMINARA:
  Niciun mecanism simplu din framework nu da w(z) la z~0.5.
  Problema fundamentala: Box e la scala PLANCK, DESI e la scala COSMOLOGICA.
  Gap-ul de 60+ ordine de marime intre cele doua scale e neinlaturat.
""")

# =====================================================================
# PART 8: Ce AM PUTEA sa functioneze
# =====================================================================
print("=" * 80)
print("PART 8: DIRECTIE POSIBILA — CC CA REZIDUU SPECTRAL")
print("=" * 80)

print(f"""
  IDEEA CARE MERITA EXPLORATA:

  In loc sa derivam Lambda din alpha^z (pattern),
  derivam Lambda din REZIDUUL spectral action.

  Spectral action: S = Tr(f(Box/Lambda_cutoff^2))
  Expandat: S = c0*f4*L^4 + c1*f2*L^2 + c2*f0 + ...

  CC = c0*f4*L^4  (leading term, enorma)
  EH = c1*f2*L^2  (gravitatie)
  YM = c2*f0       (gauge)

  Problema clasica: c0*L^4 ~ M_Planck^4 ~ 10^{{76}} GeV^4
  Observat: Lambda ~ 10^{{-47}} GeV^4
  Gap: 10^{{123}} ordine.

  CE AR PUTEA FRAMEWORK-UL SA ADUCA:

  1. Galois ANULARE: sectorul fizic si cel dark contribuie cu semne OPUSE.
     Total CC = E_phys - E_dark ~ (delta_Galois) * M_Planck^4
     Daca delta_Galois ~ alpha^57...
     Asta ar DERIVA formula Lambda = alpha^(57-alpha_s)!

  2. Concret:
     E_phys = (1/2) sum_{k in physical} d_k^2 * |mu_k|
     E_dark = (1/2) sum_{k in dark} d_k^2 * |mu_k|
     CC proportional cu |E_phys - E_dark| / E_total

  3. Am calculat mai sus:
     Delta = {Delta_total:.4f} (din pair rho2-rho6)
     E_total = {E_vac_total:.4f}
     Ratio = {ratio:.6f}
     log10(ratio) = {np.log10(ratio):.2f}

     Dar avem nevoie de ratio ~ 10^{{-122}} in Planck units.
     {np.log10(ratio):.2f} e departe de -122.

  4. TOTUSI: calculul de mai sus e cu eigenvalue MEDII, nu exacte.
     Si nu include factorul de normalizare spectral (f4, f2, f0).
     Cu normalizarea corecta, ratio-ul ar putea fi diferit.

  VERDICT: Directia Galois anulare e PROMITATOARE conceptual
  dar necesita calcul mult mai detaliat. Nu e un experiment de o zi.
""")
