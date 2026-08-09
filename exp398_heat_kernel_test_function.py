"""
EXP-398: Graph Heat Kernel as Canonical Test Function
======================================================

THE BIG PICTURE:
  Every unsolved problem (DM abundance, scalaron mass, m_e absolute value,
  Yukawa weights) involves the spectral action TEST FUNCTION f(x).
  In NCG, f is a free choice. But on a FINITE graph, the heat kernel
  IS canonical — no ambiguity.

IDEA: Use the Cayley graph heat kernel K(t,x) = exp(-t*L/Lambda^2)
evaluated at t=1 and Lambda^2=DEG=12 as the test function. This gives
fully determined moments f_0, f_2, f_4 from the 600-cell geometry alone.

QUESTIONS:
  1. What are the graph-determined moments f_0, f_2, f_4?
  2. Do they give physical predictions?
  3. Does f_4/f_2 determine the scalaron mass?
  4. Does the resulting T_reh give Omega_DM/Omega_b = 7-phi?

Author: Claude (exp398)
Date: 2026-02-25
"""

import math
import numpy as np

# ===== CONSTANTS =====
a1 = 5
b1 = 6
PHI = (1 + math.sqrt(5)) / 2
PHI_C = (1 - math.sqrt(5)) / 2
N = 120
N_eig = 9
DEG = 12
rank_E8 = 8
h_E8 = 30

m_e_eV = 0.51099895e6
M_Pl = 1.22e28  # eV

A_eq = 2 * math.pi
B_eq = -4 * a1 * PHI**4
disc = B_eq**2 - 4 * A_eq
ALPHA = (-B_eq - math.sqrt(disc)) / (2 * A_eq)
ALPHA_S = 1 / (2 * PHI**3)

print("=" * 72)
print("EXP-398: GRAPH HEAT KERNEL AS CANONICAL TEST FUNCTION")
print("=" * 72)

# =====================================================================
# PART 1: Compute ALL Cayley Laplacian eigenvalues
# =====================================================================
print("\n--- PART 1: Cayley Laplacian spectrum ---")

# Irrep data: (dim, chi_k at generator class)
dims = [1, 2, 3, 4, 5, 6, 4, 3, 2]
chi = [1, PHI, PHI, 1, 0, -1, -1, PHI_C, PHI_C]

# Adjacency eigenvalues: lambda_k = DEG * chi_k / dim_k
# Laplacian eigenvalues: L_k = DEG - lambda_k
adj_eig = []
lap_eig = []
for k in range(N_eig):
    lam = DEG * chi[k] / dims[k]
    L = DEG - lam
    adj_eig.append(lam)
    lap_eig.append(L)

print(f"\n  {'k':>2} {'dim':>4} {'chi_k':>8} {'lambda_k':>10} {'L_k':>10} {'mult':>6}")
print(f"  {'-'*44}")
for k in range(N_eig):
    mult = dims[k]**2
    print(f"  {k:>2} {dims[k]:>4} {chi[k]:>8.4f} {adj_eig[k]:>10.4f} {lap_eig[k]:>10.4f} {mult:>6}")

# Total check
total = sum(d**2 for d in dims)
print(f"\n  Total multiplicity: {total} = N = {N}")


# =====================================================================
# PART 2: Spectral action moments (Tr(L^n))
# =====================================================================
print("\n--- PART 2: Spectral moments Tr(L^n) ---")

def moment(n):
    """Tr(L^n) with multiplicities dim_k^2"""
    return sum(dims[k]**2 * lap_eig[k]**n for k in range(N_eig))

# Compute moments
for n in range(7):
    val = moment(n)
    val_over_N = val / N
    # Try to identify as framework number
    print(f"  Tr(L^{n}) = {val:>14.2f}   Tr(L^{n})/N = {val_over_N:>12.4f}")

print(f"\n  Key ratios:")
print(f"  Tr(L^0) = {moment(0):.0f} = N = {N}")
print(f"  Tr(L^1) = {moment(1):.2f} = N*DEG = {N*DEG}")
print(f"  Tr(L^2) = {moment(2):.2f}")
print(f"  Tr(L^1)/N = {moment(1)/N:.2f} = DEG = {DEG}")

# The key: Tr(L^2)/N
trL2_over_N = moment(2) / N
print(f"  Tr(L^2)/N = {trL2_over_N:.4f}")
print(f"  = DEG^2 + DEG = {DEG**2 + DEG} = {DEG*(DEG+1)}")
print(f"  = DEG*(DEG+1) = 12*13 = 156")
print(f"  Match: {abs(trL2_over_N - DEG*(DEG+1)) < 0.01}")

# Why? Tr(L^2) = Tr((DEG*I-A)^2) = N*DEG^2 - 2*DEG*Tr(A) + Tr(A^2)
# Tr(A) = 0 (always for non-trivial group Cayley graph)
# Tr(A^2) = N*DEG (each vertex has DEG neighbors, loops counted)
# So Tr(L^2) = N*(DEG^2 + DEG) = N*DEG*(DEG+1)
print(f"  THEOREM: Tr(L^2) = N*DEG*(DEG+1) for any regular Cayley graph")
print(f"           = {N}*{DEG}*{DEG+1} = {N*DEG*(DEG+1)}")


# =====================================================================
# PART 3: Heat kernel partition function
# =====================================================================
print("\n--- PART 3: Heat kernel Z(t) = Tr(exp(-t*L/DEG)) ---")

def Z(t):
    """Graph partition function at inverse temperature t"""
    return sum(dims[k]**2 * math.exp(-t * lap_eig[k] / DEG) for k in range(N_eig))

# Evaluate at various t
print(f"\n  {'t':>6} {'Z(t)':>12} {'Z(t)/N':>10}")
print(f"  {'-'*32}")
for t_val in [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    z = Z(t_val)
    print(f"  {t_val:>6.1f} {z:>12.4f} {z/N:>10.6f}")

Z1 = Z(1.0)
print(f"\n  Z(1) = {Z1:.6f}")
print(f"  Z(1)/N = {Z1/N:.6f}")

# Decompose Z(1) into positive/negative/zero sectors
Z_pos = sum(dims[k]**2 * math.exp(-lap_eig[k]/DEG) for k in range(N_eig) if adj_eig[k] > 0.01)
Z_zero = sum(dims[k]**2 * math.exp(-lap_eig[k]/DEG) for k in range(N_eig) if abs(adj_eig[k]) < 0.01)
Z_neg = sum(dims[k]**2 * math.exp(-lap_eig[k]/DEG) for k in range(N_eig) if adj_eig[k] < -0.01)

print(f"\n  Sector decomposition of Z(1):")
print(f"  Z_+ (lambda>0) = {Z_pos:.4f} (N_+ = {sum(dims[k]**2 for k in range(N_eig) if adj_eig[k]>0.01)})")
print(f"  Z_0 (lambda=0) = {Z_zero:.4f} (N_0 = {sum(dims[k]**2 for k in range(N_eig) if abs(adj_eig[k])<0.01)})")
print(f"  Z_- (lambda<0) = {Z_neg:.4f} (N_- = {sum(dims[k]**2 for k in range(N_eig) if adj_eig[k]<-0.01)})")
print(f"  Ratio Z_-/Z_+ = {Z_neg/Z_pos:.6f}")


# =====================================================================
# PART 4: NCG spectral action moments from heat kernel
# =====================================================================
print("\n--- PART 4: NCG moments from graph heat kernel ---")

# In NCG, the spectral action S = Tr(f(D^2/Lambda^2)) has Seeley-DeWitt expansion:
# S ~ f_0*a_0 + f_2*a_2*Lambda^{-2} + f_4*a_4*Lambda^{-4} + ...
# where f_k = integral_0^infty f(x)*x^{k/2-1} dx
#
# For heat kernel f(x) = exp(-x):
# f_0 = integral exp(-x)/x dx ... diverges! But on a finite graph, well-defined.
#
# On the graph, the "Seeley-DeWitt coefficients" are:
# a_0 = Tr(1) = N = 120
# a_2 = Tr(L) = N*DEG = 1440
# a_4 = Tr(L^2)/2 = N*DEG*(DEG+1)/2 = 9360

a_0 = moment(0)
a_2 = moment(1)
a_4 = moment(2) / 2  # convention: a_4 ~ Tr(L^2)/2

print(f"  Graph Seeley-DeWitt coefficients:")
print(f"    a_0 = Tr(1) = {a_0:.0f} = N")
print(f"    a_2 = Tr(L) = {a_2:.0f} = N*DEG = {N*DEG}")
print(f"    a_4 = Tr(L^2)/2 = {a_4:.0f} = N*DEG*(DEG+1)/2 = {N*DEG*(DEG+1)//2}")

# The PHYSICAL Seeley-DeWitt coefficients in the paper:
# c_k = 2*N*A_k where A_0=11, A_1=62, A_2=233
print(f"\n  Paper's Seeley-DeWitt (sec:cosmo): c_k = 2*N*A_k")
print(f"    c_0 = 2*N*A_0 = 2*120*11 = {2*N*11}")
print(f"    c_1 = 2*N*A_1 = 2*120*62 = {2*N*62}")
print(f"    c_2 = 2*N*A_2 = 2*120*233 = {2*N*233}")

# Ratios
print(f"\n  Ratios of graph moments:")
print(f"    a_2/a_0 = {a_2/a_0:.2f} = DEG = {DEG}")
print(f"    a_4/a_0 = {a_4/a_0:.2f} = DEG*(DEG+1)/2 = {DEG*(DEG+1)/2}")
print(f"    a_4/a_2 = {a_4/a_2:.2f} = (DEG+1)/2 = {(DEG+1)/2}")

# KEY: a_4/a_2 = (DEG+1)/2 = 13/2 = 6.5
print(f"\n  KEY RATIO: a_4/a_2 = (DEG+1)/2 = {(DEG+1)/2}")
print(f"  = (2*b1+1)/2 = {(2*b1+1)/2}")
print(f"  = (a1^2+1)/(a1-1) ... no, that's {(a1**2+1)/(a1-1)}")
print(f"  = b1 + 1/2")
print(f"  Note: 13 = (a1^2+1)/2 = Tr(7-phi)!")


# =====================================================================
# PART 5: Scalaron mass from graph heat kernel
# =====================================================================
print("\n--- PART 5: Scalaron mass from graph moments ---")

# In the spectral action NCG framework:
# M_scalaron^2 = (96*pi^2 * f_2) / (f_0 * c_2)
# where c_2 is the gravitational Seeley-DeWitt coefficient
# and f_0, f_2 are test function moments
#
# For the graph heat kernel, the effective ratio is:
# f_2/f_0 = a_2/a_0 = DEG = 12

# Using the paper's c_2 = 2*N*233 = 55920:
c_2_paper = 2 * N * 233
f2_over_f0 = DEG  # from graph heat kernel

# M^2 = 96*pi^2 * f_2 / (f_0 * c_2)
# = 96*pi^2 * DEG / c_2 (in units where Lambda^2 = Lambda_cutoff^2)
# But Lambda is the cutoff, which we need in physical units.

# In the framework: Lambda ~ M_Pl (the Planck scale)
# Then M_scalaron = M_Pl * sqrt(96*pi^2*DEG/c_2)

ratio = 96 * math.pi**2 * DEG / c_2_paper
M_scal_over_MPl = math.sqrt(ratio)
M_scal = M_scal_over_MPl * M_Pl

print(f"  M_scalaron^2/M_Pl^2 = 96*pi^2*DEG/c_2 = {ratio:.6f}")
print(f"  M_scalaron/M_Pl = {M_scal_over_MPl:.6f}")
print(f"  M_scalaron = {M_scal:.3e} eV = {M_scal/1e9:.3e} GeV")

# Is this a clean framework number?
z_M = math.log(M_scal / m_e_eV) / math.log(PHI)
print(f"  M_scalaron = m_e * phi^{z_M:.2f}")
print(f"  M_scalaron/m_e = {M_scal/m_e_eV:.3e}")

# Starobinsky reheating: T_reh ~ 0.3 * (M^3/M_Pl^2)^(1/2)
# (from perturbative decay of scalaron)
T_reh = 0.3 * math.sqrt(M_scal**3 / M_Pl**2)
print(f"\n  Starobinsky reheating:")
print(f"  T_reh ~ 0.3 * sqrt(M^3/M_Pl^2) = {T_reh:.3e} eV = {T_reh/1e9:.3e} GeV")
z_T = math.log(T_reh / m_e_eV) / math.log(PHI)
print(f"  T_reh = m_e * phi^{z_T:.2f}")


# =====================================================================
# PART 6: DM abundance from graph-determined T_reh
# =====================================================================
print("\n--- PART 6: DM abundance prediction ---")

# sum of Galois masses
n_f = {'e': 0, 'mu': 11, 'tau': 17, 'u': 3, 'c': 16, 't': 26, 'd': 5, 's': 11, 'b': 19}
sum_mfp = m_e_eV * sum(1/PHI**nf for nf in n_f.values())

# Gravitational freeze-in
N_species = 36  # 9 fermions * 2 spin * 2 (particle/anti)
zeta3 = 1.20206
g_star = 106.75
C_yield = 135 * zeta3 / (16 * math.pi**5 * g_star)

eta_B = 6.1e-10
m_proton = 938.272e6  # eV

Y_DM = N_species * C_yield * (T_reh / M_Pl)**3
Omega_ratio = sum_mfp * Y_DM / (m_proton * eta_B) * (2 * math.pi**2 * g_star / 45)
# Simplified: Omega_DM/Omega_b = sum(m_f') * N_sp * C * (T_reh/M_Pl)^3 / (m_p * eta_B)
# Need to be more careful with s/n_gamma factors

# Actually, the standard gravitational production yield for massless species:
# Y = n/s ~ 1.5 * (T_reh/M_Pl)^3 * (per species, spin 1/2, Boltzmann approx)
# More precisely: Y ~ T_reh^3 / (M_Pl^3 * g_star) * numerical_factor

# Let me use the order-of-magnitude estimate:
Y_DM_simple = (T_reh / M_Pl)**3
rho_DM_over_s = sum_mfp * N_species * Y_DM_simple
rho_b_over_s = m_proton * eta_B
omega_est = rho_DM_over_s / rho_b_over_s

print(f"  T_reh = {T_reh:.3e} eV")
print(f"  (T_reh/M_Pl)^3 = {(T_reh/M_Pl)**3:.3e}")
print(f"  sum(m_f') = {sum_mfp:.2f} eV")
print(f"  N_species = {N_species}")
print(f"")
print(f"  Order-of-magnitude estimate:")
print(f"  Omega_DM/Omega_b ~ sum(m_f')*N_sp*(T_reh/M_Pl)^3 / (m_p*eta_B)")
print(f"                   = {sum_mfp:.0f} * {N_species} * {(T_reh/M_Pl)**3:.2e} / ({m_proton:.2e} * {eta_B:.1e})")
print(f"                   = {omega_est:.4f}")
print(f"  Target (7-phi)   = {7-PHI:.4f}")
print(f"")
if omega_est > 0:
    print(f"  Ratio predicted/target = {omega_est/(7-PHI):.4f}")
    print(f"  Orders of magnitude off: {math.log10(omega_est/(7-PHI)):.1f}")


# =====================================================================
# PART 7: Alternative — cutoff at Lambda^2 = DEG
# =====================================================================
print("\n--- PART 7: Graph-natural scales ---")

# What if Lambda = DEG (not M_Pl)?
# Then the spectral action lives entirely on the graph.
# Physical scales emerge from the graph topology.

print(f"  Framework scale identifications:")
print(f"    DEG = 12 (graph degree)")
print(f"    N = 120 (group order)")
print(f"    h = 30 (Coxeter number)")
print(f"    a_4/a_2 = (DEG+1)/2 = 13/2")
print(f"")

# In the paper, the scalaron mass formula involves:
# M^2 = a_4/(a_2 * something)
# The ratio a_4/a_2 = 13/2 on the graph.
# 13/2 = (a1^2+1)/2 = half the dimension of the Weinberg denominator
# Note: sin^2(tW) = b1/(a1^2+1) = 6/26 = 6/(2*13)

print(f"  Connection to Weinberg angle:")
print(f"    a_4/a_2 = 13/2 = (a1^2+1)/(2*1)")
print(f"    sin^2(tW) = b1/(a1^2+1) = 6/26 = 3/13")
print(f"    Product: (a_4/a_2) * sin^2(tW) = 13/2 * 3/13 = 3/2")
print(f"    = N_gen/2 = {3/2}")

print(f"\n  INTERESTING: a_4/a_2 * sin^2(tW) = 3/2 = N_gen/2")
print(f"  This connects graph moments to gauge structure!")


# =====================================================================
# PART 8: Moment generating function
# =====================================================================
print("\n--- PART 8: Higher moments and patterns ---")

print(f"\n  {'n':>3} {'Tr(L^n)/N':>20} {'/(DEG^n)':>12} {'Factored':>30}")
for n in range(8):
    m = moment(n) / N
    ratio = m / DEG**n if DEG**n > 0 else 0
    # Try to identify the pattern
    # Tr(L^n)/N for regular graph with Tr(A^k) known
    print(f"  {n:>3} {m:>20.2f} {ratio:>12.6f}")

# The moments Tr(L^n)/N for a regular graph satisfy:
# n=0: 1
# n=1: DEG
# n=2: DEG^2 + DEG = DEG*(DEG+1)
# n=3: DEG^3 + 3*DEG^2 + DEG*t3 where t3 = Tr(A^3)/N = number of triangles per vertex
# For the 600-cell Cayley graph of 2I, t3 depends on the generating set.

# Actually Tr(A^3) counts the number of closed walks of length 3 (triangles).
# For 2I with generators = 12 elements of order ... need to compute.

print(f"\n  Tr(L^n)/(N*DEG^n) for n=0,...,6:")
for n in range(7):
    m = moment(n) / (N * DEG**n)
    print(f"  n={n}: {m:.8f}")


# =====================================================================
# PART 9: Connection to 7-phi?
# =====================================================================
print("\n--- PART 9: Does any moment ratio give 7-phi? ---")

target = 7 - PHI
print(f"  Target: 7-phi = {target:.6f}")
print()

# Try various ratios of moments
ratios = []
for n1 in range(0, 6):
    for n2 in range(n1+1, 7):
        m1 = moment(n1)
        m2 = moment(n2)
        if m1 > 0:
            r = m2 / m1
            r_over_DEG = r / DEG**(n2-n1)
            ratios.append((f"Tr(L^{n2})/Tr(L^{n1})", r, r_over_DEG))

# Also try with N, DEG, etc.
ratios.append(("a_4/(a_2*a_0^(1/2))", moment(2)/(2*moment(1)*math.sqrt(moment(0))), 0))
ratios.append(("Tr(L^3)/(N*DEG^2*phi)", moment(3)/(N*DEG**2*PHI), 0))

print(f"  {'Ratio':>30} {'Value':>12} {'/(DEG^diff)':>12}")
print(f"  {'-'*58}")
for name, val, norm in sorted(ratios, key=lambda x: abs(x[1] - target)):
    err = abs(val - target) / target * 100
    if err < 50:
        print(f"  {name:>30} {val:>12.4f} {norm:>12.6f}  ({err:.1f}%)")


# =====================================================================
# PART 10: HONEST ASSESSMENT
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
  FINDINGS:

  1. GRAPH MOMENTS are fully determined (no test function ambiguity):
     a_0 = N = 120
     a_2 = N*DEG = 1440
     a_4 = N*DEG*(DEG+1)/2 = 9360

  2. KEY RATIO: a_4/a_2 = (DEG+1)/2 = 13/2
     This connects to: Tr(7-phi) = 13 = a_4/a_2 * 2
     And: a_4/a_2 * sin^2(tW) = 3/2 = N_gen/2

  3. SCALARON MASS: M ~ M_Pl * sqrt(DEG/c_2) ~ {M_scal/1e9:.0e} GeV
     This is too high (close to M_Pl) because the graph moments
     don't encode the PHYSICAL normalization correctly.
     The graph gives TOPOLOGICAL moments, not PHYSICAL ones.

  4. DM ABUNDANCE: cannot be computed this way because T_reh
     depends on the physical normalization, not just ratios.

  CONCLUSION:
  The graph heat kernel gives well-defined RATIOS of moments,
  but the ABSOLUTE scale requires mapping to the continuum limit,
  which reintroduces the test function ambiguity.

  The ratio a_4/a_2 = 13/2 IS an interesting new identity connecting
  graph topology to the Weinberg angle. But it does NOT solve the
  test function problem.

  STATUS: The connection a_4/a_2 * sin^2(tW) = N_gen/2 is a
  NEW PATTERN worth recording, but the test function problem
  remains OPEN.
""")
