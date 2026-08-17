"""
EXP-272: DARK MATTER CANDIDATES FROM 600-CELL
==============================================
Explore what the 600-cell framework says about dark matter.
Key ideas:
1. Extra spectral modes beyond SM particles
2. The 601st DOF (conformal mode from gravity)
3. Sterile neutrinos from McKay graph
4. Topological defects on S^3
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3
ALPHA = 1/137.035999084
ALPHA_S = 1/(2*PHI**3)
m_e = 0.51099895  # MeV
m_P = 1.22093e22  # MeV (Planck mass)

print("="*72)
print("EXP-272: DARK MATTER CANDIDATES FROM 600-CELL")
print("="*72)

# === Part 1: Counting degrees of freedom ===
print("\n--- Part 1: Degrees of Freedom Budget ---")
# 600-cell simplicial complex:
V = 120; E = 720; F = 1200; C = 600
total = V + E + F + C  # = 2640
print(f"  Simplicial complex: V={V}, E={E}, F={F}, C={C}")
print(f"  Total DOF (Hilbert space): {total}")
print(f"  Betti numbers: (1, 0, 0, 1) -> 2 topological DOF")

# Gravity: 601 = 5*120+1 physical DOF
grav_DOF = a1*N + 1
print(f"  Graviton DOF: {grav_DOF} = a1*N + 1 = {a1}*{N}+1 (PRIME)")
print(f"  = 600 graviton modes + 1 conformal mode")

# SM particle count:
# 12 gauge bosons (1+3+8) -> polarizations
# 1 Higgs
# 3 generations * (2 quarks * 3 colors + 2 leptons) = 3*(6+2) = 24 Weyl fermions
# But with L+R: 48 Weyl or 24 Dirac
# Plus neutrinos: 3 or 6 (Majorana vs Dirac)
SM_bosons = 12 + 1  # gauge + Higgs
SM_fermions = 3 * (2*3 + 2) * 2  # gens * (quarks+leptons) * (L+R)
SM_total = SM_bosons + SM_fermions
print(f"\n  SM particles: {SM_bosons} bosons + {SM_fermions} fermion DOF = {SM_total}")

# McKay graph: 9 nodes * 3 generations = 27 particle types (without color)
# With color: 9*3 + quarks*3colors = ...
McKay_particles = N_eig * N_gen
print(f"  McKay graph: {N_eig} nodes * {N_gen} gens = {McKay_particles} types")

# Product spectral triple: 2640 * 9 = 23760
product_DOF = total * N_eig
print(f"  Product triple DOF: {total} * {N_eig} = {product_DOF}")

# === Part 2: Extra nodes on McKay graph ===
print("\n--- Part 2: McKay Graph Analysis ---")
# Affine E8 has 9 nodes. We assign 9 fermion types.
# The 9th node (affine extension) is special.
# In standard NCG: the 10th "node" would be right-handed neutrino
# But affine E8 has exactly 9 nodes - no room for extra particles!

# However: Coxeter labels sum to h=30
coxeter = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3])
names = ['e','u','d','s','c','b','t','tau','mu']
print(f"  Coxeter labels: {list(coxeter)}")
print(f"  Sum = {sum(coxeter)} = h(E8) = {h}")

# The Coxeter labels count multiplicity in the highest root
# Could "dark" particles correspond to HIGHER representations of E8?
# E8 has 248 generators, we use 9 nodes for fermions
# Remaining: 248 - 9 = 239 could be "dark sector"?

dim_E8 = 248
print(f"\n  dim(E8) = {dim_E8}")
print(f"  McKay nodes = {N_eig}")
print(f"  Remaining: {dim_E8 - N_eig} = {dim_E8 - N_eig}")
print(f"  Or: 248 = (N_eig-1)*(h+1) = 8*31")
print(f"  31 is prime, 8 = rank(E8)")

# === Part 3: The 5-partite structure ===
print("\n--- Part 3: 5-Partite DM Candidate ---")
# 600-cell is 5-partite: 5 classes of 24 vertices each
# 24 = |2T| (binary tetrahedral group)
# Each class could correspond to a "flavor" of dark matter

# SM uses 3 generations (from generation theorem)
# 5-partite gives 5 classes
# Visible: 3 generations, Dark: 2 extra classes?
print(f"  5-partite structure: 5 classes of 24 vertices")
print(f"  Visible generations: {N_gen}")
print(f"  Dark 'generations': {5 - N_gen} = 2")
print(f"  Ratio visible/total: {N_gen}/5 = {N_gen/5}")
print(f"  This gives Omega_DM/Omega_visible ~ 2/3 = {2.0/3:.3f}")
print(f"  Observed: Omega_DM/Omega_b = 5.36")
print(f"  NOT a match (off by factor ~8)")

# But: the 2 dark classes might have different multiplicities
# In 5-partite: each class has 24 vertices, 72 edges to each other class
# Symmetric structure -> each class has SAME coupling strength

# === Part 4: Spectral DM from coexact modes ===
print("\n--- Part 4: Spectral Dark Matter ---")
# The coexact (graviton) spectrum has a gap: lambda_min = 7 - 4*phi
gap_coexact = 7 - 4*PHI
gap_exact = 12 - 6*PHI
print(f"  Coexact gap: lambda = 7-4*phi = {gap_coexact:.6f}")
print(f"  Exact gap: lambda = 12-6*phi = {gap_exact:.6f}")
print(f"  Gap Galois: 7+4*phi = {7+4*PHI:.6f} (top eigenvalue pair)")
print(f"")
print(f"  N(coexact gap) = a1 = {a1}")
print(f"  N(exact gap) = b1^2 = {b1**2}")

# Mass of lightest graviton mode:
# m_graviton = sqrt(gap)/R where R = Hubble radius (or other macro scale)
R_Hubble = 4.4e26  # meters
hbar_c = 197.3e-15  # MeV*m
m_grav = np.sqrt(gap_coexact) * hbar_c / R_Hubble  # in MeV
print(f"\n  Lightest graviton mass:")
print(f"    m = sqrt({gap_coexact:.4f}) * hbar*c / R_Hubble")
print(f"    = {m_grav:.2e} MeV = {m_grav*1e6:.2e} eV")
print(f"    This is ~ 10^-34 eV (sub-LIGO)")

# But: at smaller R, the mass increases
# For R ~ 1/Lambda_QCD ~ 1 fm:
R_QCD = 1e-15  # meters
m_grav_QCD = np.sqrt(gap_coexact) * hbar_c / R_QCD
print(f"  At QCD scale R ~ 1 fm:")
print(f"    m = {m_grav_QCD:.2f} MeV")

# For DM: need m ~ keV to GeV range
# R ~ 10^-18 to 10^-15 m (sub-nuclear to nuclear)
for R_name, R_val in [('Hubble', 4.4e26), ('galaxy', 3e20), ('solar', 1.5e11),
                       ('Earth', 6.4e6), ('1 m', 1), ('1 fm', 1e-15), ('Planck', 1.6e-35)]:
    m_val = np.sqrt(gap_coexact) * hbar_c / R_val
    print(f"    R = {R_name:10s}: m = {m_val:.2e} MeV")

# === Part 5: Topological dark matter ===
print("\n--- Part 5: Topological DM ---")
# S^3 has pi_3(S^3) = Z -> topological defects (instantons)
# On 600-cell discretization: "instantons" = non-trivial maps S^3 -> S^3
# The instanton number is quantized

# The 600-cell fundamental group is trivial (it triangulates S^3)
# But: the HOLONOMY group of the 600-cell is binary icosahedral = 2I
# |2I| = 120 = N

# Dark matter as topological modes:
print(f"  pi_3(S^3) = Z (instanton number)")
print(f"  Binary icosahedral group: |2I| = {N}")
print(f"  On 600-cell: N = {N} vertices")
print(f"  Instanton-like configurations: maps preserving 600-cell structure")
print(f"  These would be MASSIVE, stable, weakly interacting")

# === Part 6: Z2 dark sector from Galois ===
print("\n--- Part 6: Galois Dark Sector ---")
# The Galois symmetry phi <-> phi' maps Z[phi] to itself
# Physical spectrum (phi sector) vs shadow spectrum (phi' sector)
# The phi' sector is NOT directly observable (Galois conjugate)

# If Galois symmetry is EXACT: the phi' sector contains "dark" copies
# of all particles, but with different (Galois-conjugate) quantum numbers

# n(phi) = 5a + 6b (physical mass exponent)
# n(phi') = 5a - b (Galois conjugate)
# For n > 0: physical particles
# For n' > 0: dark particles?

fermion_ab = {'e':(0,0), 'mu':(1,1), 'tau':(1,2), 'u':(3,-2),
              'c':(2,1), 't':(4,1), 'd':(1,0), 's':(1,1), 'b':(-1,4)}

print(f"  Galois shadow masses (n' = 5a - b):")
print(f"  {'Fermion':>8s}  {'(a,b)':>8s}  {'n':>4s}  {'n_prime':>7s}  {'m_phys':>10s}  {'m_dark':>10s}")
for name in ['e','mu','tau','u','d','s','c','t','b']:
    a, b = fermion_ab[name]
    n = a1*a + b1*b
    n_prime = a1*a - b
    m_phys = m_e * PHI**n  # MeV
    m_dark = m_e * PHI**n_prime if n_prime >= 0 else m_e / PHI**(-n_prime)
    print(f"  {name:>8s}  ({a:2d},{b:2d})  {n:4d}  {n_prime:7d}  {m_phys:10.3f}  {m_dark:10.6f}")

# The Galois dark sector: different masses but same interactions?
# No: interactions also involve phi -> the Galois conjugate coupling is different
# alpha_s' = 1/(2*phi'^3) = 1/(2*(-1/phi)^3) = -phi^3/2 < 0 (!)
# NEGATIVE coupling => dark sector is CONFINING in a different way

alpha_s_dark = 1/(2*(-1/PHI)**3)
print(f"\n  Galois dark coupling:")
print(f"    alpha_s' = 1/(2*phi'^3) = {alpha_s_dark:.4f} (NEGATIVE!)")
print(f"    This means: the dark sector has REPULSIVE strong force")
print(f"    -> No bound states -> dark 'quarks' are free!")
print(f"    -> This would be a GAS of weakly interacting particles")

# === Part 7: Dark matter mass prediction ===
print("\n--- Part 7: Mass Prediction ---")
# If DM = Galois partner of a known particle:
# Best candidate: a particle with n' in the right range

# WIMP miracle: m_DM ~ 100 GeV = 10^5 MeV
# Need phi^n' * m_e ~ 10^5
# n' ~ ln(10^5/0.511)/ln(phi) ~ 25.3
# n' = 25: m ~ m_e * phi^25 = 85.8 GeV (Z-like mass!)

n_target = 25  # same as n_Z
m_DM_25 = m_e * PHI**25 / 1000  # GeV
print(f"  WIMP-scale DM: n' = 25 -> m = m_e*phi^25 = {m_DM_25:.1f} GeV")
print(f"  This is the Z boson exponent!")
print(f"  But n'=25 requires (a,b) with 5a-b=25")
print(f"  Solutions: a=5,b=0 or a=4,b=-5 or a=3,b=-10...")
print(f"  a=5,b=0: n = 5*5+6*0 = 25 (physical n = n_Z!)")
print(f"  -> The Galois partner of the n=25 particle IS the n'=25 particle")
print(f"  -> SELF-CONJUGATE under Galois!")

# Check all fermions for self-conjugate:
print(f"\n  Galois self-conjugate fermions (n = n'):")
for name in ['e','mu','tau','u','d','s','c','t','b']:
    a, b = fermion_ab[name]
    n = a1*a + b1*b
    n_prime = a1*a - b
    if n == n_prime:
        print(f"    {name}: n = n' = {n} (b=0)")

# Down quark: a=1,b=0 -> n=5, n'=5. Self-conjugate!
# Electron: a=0,b=0 -> n=0, n'=0. Self-conjugate!
# These have b=0 (rational in Z[phi])

# For n'=25: not a known fermion, but could be a BOSON
# The Z boson has n_Z = 25 -> it IS self-conjugate
# m_Z = m_e * phi^25 * alpha(m_Z)/alpha(0) (with running)

# Alternative: dark matter mass from spectral gap
# If DM = lightest mode above the Higgs:
print(f"\n  Alternative DM masses:")
print(f"    From n=35 (neutrino complement): m = 2*m_e*phi^35 = heavy")
print(f"    n_Z + n_nu = 60 = N/2")
print(f"    DM at n = 60 - n_visible?")

# === Part 8: Relic abundance ===
print("\n--- Part 8: Relic Abundance Estimate ---")
# For WIMP with m ~ 100 GeV, sigma_ann ~ pi*alpha^2/m^2
# Omega_DM ~ 0.1 pb / <sigma*v>

# If DM couples with alpha_s' (Galois dark coupling):
# But alpha_s' is negative -> confusing
# Use |alpha_s'| = phi^3/2 = 2.118

print(f"  Galois coupling: |alpha_s'| = phi^3/2 = {PHI**3/2:.4f}")
print(f"  This is NOT perturbative (> 1)!")
print(f"  -> Galois dark sector is STRONGLY coupled")
print(f"  -> Cannot compute perturbatively")
print(f"  -> Lattice/non-perturbative methods needed")

# Instead: if DM has only gravitational interactions
# Then it's like a PLANCKION
# sigma ~ G_N^2 * m^2 ~ (m/m_P)^4 / m^2
# Way too small for thermal relic

# === Part 9: Observable signatures ===
print("\n--- Part 9: Predictions/Signatures ---")
print(f"  IF DM = Galois dark sector:")
print(f"    1. Mass spectrum is Galois-conjugate (n' = 5a-b)")
print(f"    2. Strong coupling is REPULSIVE (alpha_s' < 0)")
print(f"    3. No bound states -> free dark quarks")
print(f"    4. Interacts via gravity + possibly Higgs portal")
print(f"    5. Self-conjugate particles (b=0): e, d are same in both sectors")
print(f"")
print(f"  IF DM = extra 5-partite classes:")
print(f"    1. Same mass spectrum as visible matter")
print(f"    2. But different 'color' under 5-partite labeling")
print(f"    3. Ratio: 2/3 of visible (factor 8 too small)")
print(f"")
print(f"  IF DM = topological modes:")
print(f"    1. Instanton-like configurations on S^3")
print(f"    2. Mass ~ instanton action * Lambda")
print(f"    3. Stable by topological charge conservation")

# === Part 10: Summary ===
print("\n--- Part 10: Summary ---")
print(f"  DARK MATTER STATUS: SPECULATIVE")
print(f"  Three candidate mechanisms:")
print(f"    1. Galois dark sector (phi' copies of SM)")
print(f"    2. 5-partite extra classes (2 dark 'generations')")
print(f"    3. Topological modes on S^3")
print(f"  None fully worked out. The Galois mechanism is most")
print(f"  interesting (predicts mass spectrum) but has issues")
print(f"  (strong coupling is negative/non-perturbative).")
print(f"  The 5-partite gives wrong ratio.")
print(f"  OPEN PROBLEM: needs more work.")

print("\n" + "="*72)
print("EXP-272 COMPLETE")
print("="*72)
