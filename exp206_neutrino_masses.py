#!/usr/bin/env python3
"""
EXP-206: NEUTRINO MASSES FROM 600-CELL STRUCTURE
==================================================
The biggest gap in the theory: neutrino masses and PMNS mixing.

What we know:
- Mass formula: m_f = m_e * phi^(a*d_eff + b*k_eff) works for 9 fermions
- Generation theorem: b in {0,1,2}, N_gen = 3
- CKM exponents: n(b1,b2) = 9*b2 - 5*b1 - 6
- Galois complementarity: leptons=shadow sector, quarks=real sector
- Neutrinos are MASSLESS in bare theory (no right-handed partner)
- Experimentally: m_nu ~ 0.01-0.1 eV (10^6 times lighter than electron)

APPROACHES:
A. Seesaw from Galois structure: m_nu ~ m_charged^2 / M_R
   where M_R comes from the dual (phi') sector
B. Norm-based: |N(z)| for non-unit z gives heavy particles
   Neutrinos might have |N(z)| >> 1 (heavy Majorana partner)
C. Schur complement: neutrinos as "integrated out" degrees of freedom
D. Perturbative: m_nu/m_e ~ phi^(-large n) for some geometric n
E. PMNS angles from the CKM formula with different coefficients

EXPERIMENTAL DATA:
  Dm^2_21 = 7.53e-5 eV^2 (solar)
  Dm^2_31 = 2.453e-3 eV^2 (atmospheric, normal ordering)
  theta_12 = 33.44 deg
  theta_23 = 49.0 deg (close to maximal!)
  theta_13 = 8.57 deg

STATUS: EXPLORATORY (all previous attempts NEGATIVE/PARTIAL)
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math

PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179
SIN2_TW = 0.23122
M_E = 0.51099895  # MeV
M_W = 80377.0  # MeV
M_Z = 91187.6  # MeV

# Neutrino data
DM2_SOL = 7.53e-5   # eV^2
DM2_ATM = 2.453e-3  # eV^2
TH12_NU = math.radians(33.44)
TH23_NU = math.radians(49.0)
TH13_NU = math.radians(8.57)

print("=" * 72)
print("EXP-206: NEUTRINO MASSES FROM 600-CELL STRUCTURE")
print("=" * 72)
print()

# ============================================================
# PART A: MASS SCALE FROM GALOIS SEESAW
# ============================================================
print("=" * 72)
print("PART A: GALOIS SEESAW MECHANISM")
print("=" * 72)
print()

# In the 600-cell: z = a + b*phi, z' = a + b*phi' = a + b*(1-phi) = (a+b) - b*phi
# The Galois conjugate maps phi -> phi' = -1/phi
# For lepton masses: corrections are in the phi' sector (Galois shadow)
# Idea: the Majorana mass comes from the Galois conjugate sector
# M_R ~ m_f * phi^(large number from phi' sector)

# The seesaw: m_nu = m_D^2 / M_R
# where m_D is the Dirac mass (same order as charged lepton)
# and M_R is the Majorana mass (very large)

# In our framework: m_D = m_e * phi^n (where n comes from z = a+b*phi)
# What about M_R? It should come from the Galois conjugate.
# z' = sigma(z). For units phi^k: sigma(phi^k) = (phi')^k = (-1/phi)^k

# Key idea: M_R ~ m_e * phi^(-n) (Galois conjugate reverses sign of exponent?)
# Then m_nu = (m_e * phi^n)^2 / (m_e * phi^(-n)) = m_e * phi^(3n)
# This gives LARGER mass, not smaller. WRONG direction.

# Alternative: M_R ~ M_Planck? No, we don't have Planck scale.
# Alternative: M_R ~ m_e * phi^(degree * something)?

# Let's think differently.
# In the 600-cell, the "dual" (120-cell) has 600 vertices.
# The mass ratio m_W/m_e ~ 1.57e5.
# phi^24 ~ 1.03e5. Close but not exact.
# phi^25 ~ 1.67e5. Also close.

# Key: the 600-cell has 720 edges.
# phi^(720/some_number)?

# Actually, let me try: what exponent n gives the right neutrino mass scale?
# m_nu ~ 0.05 eV = 5e-8 MeV
# m_nu/m_e = 5e-8 / 0.511 ~ 1e-7
# phi^(-n) = 1e-7 => n = 7 * ln(10) / ln(phi) = 7 * 2.3026 / 0.4812 = 33.5
# So n ~ 33-34

print("  Neutrino mass scale:")
print("  m_nu ~ 0.05 eV = 5e-8 MeV")
print("  m_nu/m_e ~ %.2e" % (5e-8 / M_E * 1e3))
print("  Need phi^(-n) ~ %.2e" % (5e-8 / M_E * 1e3))
n_needed = -math.log(5e-8 / (M_E * 1e-3)) / math.log(PHI)
print("  n = %.2f (need phi^-%.0f ~ %.2e)" % (n_needed, round(n_needed), PHI**(-round(n_needed))))
print()

# Interesting numbers near 33-34:
print("  Geometric numbers near n ~ 33:")
print("  30 = 5 * 6 = diam * b_1")
print("  33 = 3 * 11 = N_gen * (diam + b_1)")
print("  35 = 5 * 7 = diam * 7")
print("  36 = 6^2 = b_1^2 = eigenvalue multiplicity (k=5)")
print("  33 = n_12 + n_23 + n_13 + n_tot = 3+7+12+11...")
print("  22 = 3+7+12 (sum of CKM exponents)")
print()

# Test: phi^-33 / m_e for each generation
for name, n_guess in [("33=3*11", 33), ("34=2*17", 34), ("30=5*6", 30),
                       ("36=6^2", 36), ("22=CKM_sum", 22)]:
    m_nu = M_E * 1e3 * PHI**(-n_guess)  # in eV
    print("  phi^-%s: m_nu = %.4e eV (target ~0.05 eV)" % (name, m_nu))

print()

# ============================================================
# PART B: MASS RATIOS FROM NORM STRUCTURE
# ============================================================
print("=" * 72)
print("PART B: NEUTRINO MASS RATIOS")
print("=" * 72)
print()

# Even if we can't get the absolute scale, can we get the ratios?
# Dm^2_21/Dm^2_31 = 7.53e-5 / 2.453e-3 = 0.0307
# For normal ordering: m3 >> m2 > m1
# m2^2 - m1^2 = Dm^2_21
# m3^2 - m1^2 = Dm^2_31

# If m1 ~ 0, then m2 ~ sqrt(Dm^2_21) ~ 8.7 meV, m3 ~ sqrt(Dm^2_31) ~ 49.5 meV
m2_approx = math.sqrt(DM2_SOL) * 1000  # in meV
m3_approx = math.sqrt(DM2_ATM) * 1000  # in meV
print("  For m1 ~ 0 (normal ordering):")
print("  m2 ~ %.2f meV" % m2_approx)
print("  m3 ~ %.2f meV" % m3_approx)
print("  m3/m2 ~ %.3f" % (m3_approx / m2_approx))
print()

# Compare m3/m2 with phi powers
r32 = m3_approx / m2_approx
print("  m3/m2 = %.4f" % r32)
print("  phi^3 = %.4f" % PHI**3)
print("  phi^2 = %.4f" % PHI**2)
print("  phi^2.5 = %.4f" % PHI**2.5)
n_r = math.log(r32) / math.log(PHI)
print("  m3/m2 = phi^%.4f" % n_r)
print()

# The charged lepton ratios for comparison:
m_mu = 105.6584  # MeV
m_tau = 1776.86  # MeV
print("  Charged lepton ratios:")
print("  m_mu/m_e = %.2f = phi^%.4f" % (m_mu/M_E, math.log(m_mu/M_E)/math.log(PHI)))
print("  m_tau/m_e = %.2f = phi^%.4f" % (m_tau/M_E, math.log(m_tau/M_E)/math.log(PHI)))
print("  m_tau/m_mu = %.2f = phi^%.4f" % (m_tau/m_mu, math.log(m_tau/m_mu)/math.log(PHI)))
print()

# For charged leptons: (a,b) = e(0,0), mu(1,1), tau(1,2)
# n = 5a+6b: e=0, mu=11, tau=17
# m_mu/m_e ~ phi^11, m_tau/m_e ~ phi^17
# Ratio: m_tau/m_mu ~ phi^6 = phi^b_1

# For neutrinos: what (a,b) would give the right ratio?
# m3/m2 = phi^3.72 ~ phi^(some simple number)
# Actually 3.72 is not very clean.

# Dm^2 ratio is cleaner:
r_dm2 = DM2_SOL / DM2_ATM
print("  Dm^2_sol / Dm^2_atm = %.6f" % r_dm2)
print("  = phi^%.4f" % (math.log(r_dm2) / math.log(PHI)))
print("  = 1/phi^%.4f" % (-math.log(r_dm2) / math.log(PHI)))
print()

# Test: is Dm^2_sol/Dm^2_atm ~ 1/phi^7?
for n in range(1, 15):
    val = PHI**(-n)
    err = abs(val - r_dm2) / r_dm2 * 100
    if err < 50:
        print("  phi^-%d = %.6f (error %.1f%%)" % (n, val, err))

# Or some combination
print()
print("  Testing Dm^2 ratio against geometric expressions:")
tests = [
    ("1/(8*phi)", 1/(8*PHI)),
    ("1/(6*phi^2)", 1/(6*PHI**2)),
    ("alpha", ALPHA),
    ("alpha/phi", ALPHA/PHI),
    ("sin^2(tW)/phi^4", SIN2_TW/PHI**4),
    ("1/phi^7", PHI**(-7)),
    ("1/(2*phi^5)", 1/(2*PHI**5)),
    ("1/(phi^4*pi)", 1/(PHI**4 * math.pi)),
    ("alpha*phi^3", ALPHA*PHI**3),
    ("3/(2*phi^5)", 3/(2*PHI**5)),
    ("(alpha_s*phi)^2", (ALPHA_S*PHI)**2),
]
for name, val in tests:
    err = abs(val - r_dm2) / r_dm2 * 100
    print("  %s = %.6f (error %.2f%%)" % (name, val, err))


# ============================================================
# PART C: PMNS ANGLES FROM CKM FORMULA
# ============================================================
print()
print("=" * 72)
print("PART C: PMNS ANGLES")
print("=" * 72)
print()

# CKM formula: n(b1,b2) = 9*b2 - 5*b1 - 6
# CKM angles: theta_ij = arctan(phi^(-n(b_i,b_j)))
# PMNS angles are MUCH LARGER than CKM angles.
# theta_12(PMNS) ~ 33 deg vs theta_12(CKM) ~ 13 deg
# theta_23(PMNS) ~ 49 deg vs theta_23(CKM) ~ 2 deg
# theta_13(PMNS) ~ 8.5 deg vs theta_13(CKM) ~ 0.2 deg

# From exp191: CKM_n = 3 * PMNS_n for theta_12 and theta_13
# theta_12(CKM) uses phi^-3, theta_12(PMNS) uses phi^-1?
# theta_13(CKM) uses phi^-12, theta_13(PMNS) uses phi^-4?

print("  PMNS angles vs arctan(phi^-n):")
for n in range(0, 10):
    angle = math.degrees(math.atan(PHI**(-n)))
    print("  arctan(phi^-%d) = %.4f deg" % (n, angle))

print()
print("  Target PMNS angles:")
print("  theta_12 = %.2f deg" % math.degrees(TH12_NU))
print("  theta_23 = %.2f deg" % math.degrees(TH23_NU))
print("  theta_13 = %.2f deg" % math.degrees(TH13_NU))
print()

# The PMNS angles don't fit arctan(phi^-n) for any integer n!
# theta_12 = 33.44 deg: arctan(phi^-1) = 31.72 deg (diff 5%)
# theta_23 = 49.0 deg: close to 45 deg = arctan(1) = arctan(phi^0)
# theta_13 = 8.57 deg: arctan(phi^-3) = 13.28 deg (not close)

# Let's try: PMNS = CKM/3 rule (from exp191)
# CKM n: 3, 7, 12
# PMNS n: 1, 7/3, 4
# Only theta_12 works cleanly (3/3=1)

# Alternative: try arctan(phi^-n * correction)
print("  Trying arctan(phi^-n * (1+x)) for PMNS:")
for angle_name, angle_val in [("th12", TH12_NU), ("th23", TH23_NU), ("th13", TH13_NU)]:
    print("\n  %s = %.4f deg:" % (angle_name, math.degrees(angle_val)))
    tan_val = math.tan(angle_val)
    for n in range(-1, 8):
        bare = PHI**(-n)
        if bare > 0:
            corr = tan_val / bare - 1
            if abs(corr) < 2:
                print("    n=%d: correction = %.6f" % (n, corr))
                # Check if correction matches known constants
                for cname, cval in [("0", 0), ("alpha", ALPHA), ("-alpha", -ALPHA),
                                     ("alpha_s/pi", ALPHA_S/math.pi),
                                     ("-2alpha_s/3pi", -2*ALPHA_S/(3*math.pi)),
                                     ("1/phi", 1/PHI), ("-1/phi", -1/PHI),
                                     ("sin2tW", SIN2_TW), ("-sin2tW", -SIN2_TW),
                                     ("1/phi^2", 1/PHI**2), ("-1/phi^2", -1/PHI**2),
                                     ("phi-1", PHI-1), ("2-phi", 2-PHI),
                                     ("1/2", 0.5), ("-1/2", -0.5),
                                     ("1/3", 1/3.0), ("2/3", 2/3.0),
                                     ("1/5", 0.2), ("1/sqrt(5)", 1/math.sqrt(5)),
                                     ("phi/5", PHI/5),
                                     ("-1/5", -0.2)]:
                    if abs(corr - cval) < 0.05 and cname != "0":
                        print("      ~= %s (%.6f, diff %.4f)" % (cname, cval, corr-cval))


# ============================================================
# PART D: TRIBIMAXIMAL AND PHI-DEPENDENT MIXING
# ============================================================
print()
print("=" * 72)
print("PART D: TRIBIMAXIMAL AND GEOMETRIC MIXING MATRICES")
print("=" * 72)
print()

# Tribimaximal mixing: sin^2(th12)=1/3, sin^2(th23)=1/2, th13=0
# This is close but not exact for PMNS
print("  Tribimaximal prediction vs experiment:")
print("  sin^2(th12): 1/3 = 0.3333 vs %.4f (exp)" % math.sin(TH12_NU)**2)
print("  sin^2(th23): 1/2 = 0.5000 vs %.4f (exp)" % math.sin(TH23_NU)**2)
print("  sin^2(th13): 0 vs %.4f (exp)" % math.sin(TH13_NU)**2)
print()

# Can we get these from 600-cell geometry?
# 1/3 = 1/N_gen. This is natural.
# 1/2 = maximal mixing.

# Golden ratio mixing: sin^2(th12) = 1/phi^2 = 1/(phi+1) = phi-1 ~ 0.382?
# No, 1/phi^2 = 0.382. But experimental = 0.304. Not great.

# Actually: sin(th12) = 1/sqrt(3) for tribimaximal
# Try: sin(th12) = 1/sqrt(phi+2) = 1/sqrt(3.618)?
sin2_test = 1/(PHI+2)
print("  sin^2(th12) = 1/(phi+2) = %.4f vs %.4f (exp)" % (sin2_test, math.sin(TH12_NU)**2))
print("  diff: %.2f%%" % (100*abs(sin2_test - math.sin(TH12_NU)**2)/math.sin(TH12_NU)**2))

# Or: sin^2(th12) = 2/(phi+5)?
sin2_test2 = 2/(PHI+5)
print("  sin^2(th12) = 2/(phi+5) = %.4f vs %.4f" % (sin2_test2, math.sin(TH12_NU)**2))
print("  diff: %.2f%%" % (100*abs(sin2_test2 - math.sin(TH12_NU)**2)/math.sin(TH12_NU)**2))

# Or from the Schur complement: sin^2 = 6/26 was for Weinberg
# Try: sin^2(th12_PMNS) = 8/26?
sin2_test3 = 8/26
print("  sin^2(th12) = 8/26 = %.4f vs %.4f" % (sin2_test3, math.sin(TH12_NU)**2))
print("  diff: %.2f%%" % (100*abs(sin2_test3 - math.sin(TH12_NU)**2)/math.sin(TH12_NU)**2))

# 8/26 = 4/13 = 0.3077. Experimental = 0.3040. That's 1.2% off!
print("  ==> sin^2(th12_PMNS) = 4/13 is 1.2%% off!")
print()

# What about th23?
# sin^2(th23) ~ 0.569
print("  sin^2(th23_PMNS) = %.4f" % math.sin(TH23_NU)**2)
# 0.569 ~ close to something?
for num in range(1, 30):
    for den in range(num+1, 50):
        val = num/den
        err = abs(val - math.sin(TH23_NU)**2) / math.sin(TH23_NU)**2 * 100
        if err < 1.0:
            print("  %d/%d = %.4f (%.2f%%)" % (num, den, val, err))

# sin^2(th13) ~ 0.0222
print("\n  sin^2(th13_PMNS) = %.4f" % math.sin(TH13_NU)**2)
for name, val in [("alpha/pi", ALPHA/math.pi), ("1/(6*phi^2)", 1/(6*PHI**2)),
                   ("alpha*3", 3*ALPHA), ("sin2tW/10", SIN2_TW/10),
                   ("1/45", 1/45.0), ("1/(6*b_1+1)", 1/37.0),
                   ("alpha_s/(a_1+b_1-2)", ALPHA_S/9),
                   ("1/(5*9)", 1/45.0),
                   ("phi^-8", PHI**(-8)),
                   ("1/(2*phi^5)", 1/(2*PHI**5))]:
    err = abs(val - math.sin(TH13_NU)**2) / math.sin(TH13_NU)**2 * 100
    if err < 30:
        print("  %s = %.6f (%.2f%%)" % (name, val, err))

# ============================================================
# PART E: INTERSECTION-NUMBER APPROACH TO PMNS
# ============================================================
print()
print("=" * 72)
print("PART E: NON-DISTANCE-REGULAR STRUCTURE AND PMNS")
print("=" * 72)
print()

# The 600-cell is NOT distance-regular (exp205 finding!)
# At distance d=2: some vertices see b=3, others b=6
# At distance d=3: some see b=2, others b=6
# This non-uniformity might encode PMNS mixing!

# Let's count the "types" at each distance
from itertools import permutations as iperms

def build_600cell():
    verts = set()
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.add(tuple(v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))
    half = 0.5
    phi_half = PHI / 2.0
    iphi_half = 1.0 / (2.0 * PHI)
    base = [0.0, half, phi_half, iphi_half]
    even_perms = []
    for p in iperms(range(4)):
        inv = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if p[a] > p[b]:
                    inv += 1
        if inv % 2 == 0:
            even_perms.append(p)
    for p in even_perms:
        cu = [base[p[0]], base[p[1]], base[p[2]], base[p[3]]]
        nz = [i for i in range(4) if abs(cu[i]) > 1e-12]
        for sb in range(1 << len(nz)):
            v = list(cu)
            for k, pos in enumerate(nz):
                if (sb >> k) & 1:
                    v[pos] = -v[pos]
            verts.add(tuple(round(x, 12) for x in v))
    return sorted(verts)

vertices = build_600cell()
N = len(vertices)
varray = np.array(vertices)
dot_mat = varray @ varray.T

# All possible inner products give distance classes
# Inner products between vertices of 600-cell
inner_prods = set()
for i in range(N):
    for j in range(i+1, N):
        ip = round(dot_mat[i,j], 6)
        inner_prods.add(ip)

inner_prods_sorted = sorted(inner_prods, reverse=True)
print("  All inner products between 600-cell vertices:")
for ip in inner_prods_sorted:
    count = 0
    for j in range(1, N):
        if abs(dot_mat[0, j] - ip) < 1e-4:
            count += 1
    if count > 0:
        print("  ip = %8.6f: %d vertices" % (ip, count))

# The non-distance-regularity means different vertices at same graph distance
# can have different inner products!
# Let's check: at graph distance 2, what are the inner products?
adj_bool = np.abs(dot_mat - PHI / 2.0) < 1e-8
np.fill_diagonal(adj_bool, False)

# BFS distances from vertex 0
dist0 = [-1] * N
dist0[0] = 0
queue = [0]
while queue:
    v = queue.pop(0)
    for j in range(N):
        if adj_bool[v, j] and dist0[j] == -1:
            dist0[j] = dist0[v] + 1
            queue.append(j)

print("\n  Inner products at each graph distance from vertex 0:")
for d in range(6):
    verts_d = [j for j in range(N) if dist0[j] == d and j != 0]
    if not verts_d:
        continue
    ips_d = [round(dot_mat[0, j], 6) for j in verts_d]
    ip_counts = {}
    for ip in ips_d:
        found = False
        for key in ip_counts:
            if abs(ip - key) < 1e-4:
                ip_counts[key] += 1
                found = True
                break
        if not found:
            ip_counts[ip] = 1

    print("  d=%d (%d vertices):" % (d, len(verts_d)))
    for ip, cnt in sorted(ip_counts.items(), reverse=True):
        # What is this inner product?
        # ip = cos(angle) on S^3
        # Check if a + b*phi form
        desc = ""
        for a in range(-5, 5):
            for b in range(-5, 5):
                test = (a + b * PHI) / 2.0
                if abs(test - ip) < 1e-4:
                    desc = "(%d+%d*phi)/2" % (a, b)
                    break
            if desc:
                break
        if not desc:
            for a in range(-5, 5):
                for b in range(-5, 5):
                    test = a + b * PHI
                    if abs(test - ip) < 1e-4:
                        desc = "%d+%d*phi" % (a, b)
                        break
                if desc:
                    break
        print("    ip = %.6f (%s): %d vertices" % (ip, desc if desc else "?", cnt))

# The "association scheme" classes
# Since 600-cell is NOT distance-regular, the association scheme has more classes
# than the graph distance would suggest. Let's count.
print("\n  Association scheme classes (by inner product):")
all_ips = {}
for j in range(1, N):
    ip = round(dot_mat[0, j], 6)
    found = False
    for key in all_ips:
        if abs(ip - key) < 1e-4:
            all_ips[key] += 1
            found = True
            break
    if not found:
        all_ips[ip] = 1

n_classes = len(all_ips)
print("  Number of classes: %d" % n_classes)
print("  (Distance-regular would need %d = diam)" % 5)
print("  (We need %d for 9 eigenvalues)" % 8)  # n_classes = n_eigenvalues - 1
print()

# The additional classes (beyond distance-regular) encode the mixing structure!
# At d=2: TWO types of vertices. At d=3: TWO types.
# This gives 2 + 2 = 4 extra classes beyond distance-regular.
# Total: 5 (distance-regular) + 4 (extra) = 9 classes? Let's check.
# Actually, we have 1(d=0) + 1(d=1) + 2(d=2) + 2(d=3) + 2(d=4) + 1(d=5) = 9
# But that's 9 total partitions including d=0 self.
# Association scheme classes = 9 - 1 = 8 (exclude self). Eigenvalues = 8 + 1 = 9.
print("  MATCH: 8 association scheme classes + 1 (identity) = 9 eigenvalues!")


# ============================================================
# PART F: SCHUR COMPLEMENT FOR NEUTRINOS
# ============================================================
print()
print("=" * 72)
print("PART F: SCHUR COMPLEMENT - NEUTRINO AS INTEGRATED OUT")
print("=" * 72)
print()

# Idea: In the Schur complement L_eff (24x24), integrating out C,
# we got eigenvalues 0(x1), 7(x4), 450/43(x9), 126/11(x8), 12(x2).
# The multiplicities (1, 4, 9, 8, 2) could relate to particle content:
# 1 = singlet (photon?)
# 4 = ?
# 9 = 3 families x 3 colors?
# 8 = gluons?
# 2 = W+, W-?

# If this interpretation is correct, where are neutrinos?
# They might be in the KERNEL (eigenvalue 0, mult 1)
# Or they might be encoded in the eigenvalue SPLITTINGS

print("  Schur complement L_eff eigenvalue multiplicities:")
print("  0 (x1), 7 (x4), 450/43 (x9), 126/11 (x8), 12 (x2)")
print()
print("  Particle interpretation attempt:")
print("  1 = photon (massless, singlet)")
print("  4 = ??? (4 components of something)")
print("  9 = 3x3 (3 families x 3 colors) or 3^2 (adjoint of SU(3)?)")
print("  8 = 8 gluons (adjoint of SU(3))")
print("  2 = W+, W- (charged weak bosons)")
print()
print("  Problem: neutrinos don't fit cleanly.")
print("  1+4+9+8+2 = 24 gauge d.o.f. Neutrinos are fermionic.")
print()

# Actually, the Schur complement integrates out fermions (C vertices)
# to get an effective gauge theory. Neutrino masses would appear
# as modifications to the effective Laplacian when solitons are present.

# Let me try: the RATIO of Schur eigenvalues
print("  Ratios of Schur eigenvalues:")
ev_schur = [0, 7, 450/43, 126/11, 12]
for i in range(1, len(ev_schur)):
    for j in range(i+1, len(ev_schur)):
        r = ev_schur[j] / ev_schur[i]
        print("    ev[%d]/ev[%d] = %.6f/%.6f = %.6f" % (
            j, i, ev_schur[j], ev_schur[i], r))

# ============================================================
# PART G: SEESAW FROM Z AND Z'
# ============================================================
print()
print("=" * 72)
print("PART G: SEESAW FROM GALOIS STRUCTURE")
print("=" * 72)
print()

# For each generation, z = 1 + b*phi (unit in Z[phi])
# Galois conjugate: z' = 1 + b*phi' = 1 + b*(1-phi) = (1+b) - b*phi
# Dirac mass: m_D propto phi^n where n = 5a + 6b (charged lepton exponent)
# For leptons: (a,b) = e(0,0), mu(1,1), tau(1,2)
# n_e = 0, n_mu = 11, n_tau = 17

# Seesaw: m_nu = m_D^2 / M_R
# If M_R comes from Galois conjugate: M_R propto phi^(n')
# where n' involves a and b through the conjugate

# For z = a + b*phi: z' = a + b*phi' = a + b - b*phi = (a+b) - b*phi
# So (a, b) -> (a+b, -b) under Galois
# n -> n' = 5(a+b) + 6(-b) = 5a + 5b - 6b = 5a - b

print("  Galois conjugate of exponents:")
print("  n = 5a + 6b -> n' = 5a - b")
print()

leptons = [("e", 0, 0), ("mu", 1, 1), ("tau", 1, 2)]
for name, a, b in leptons:
    n = 5*a + 6*b
    np_conj = 5*a - b
    print("  %s: (a,b)=(%d,%d), n=%d, n'=%d" % (name, a, b, n, np_conj))
    print("    m_D ~ phi^%d, M_R ~ phi^%d" % (n, np_conj))
    if np_conj != 0:
        n_seesaw = 2*n - np_conj
        print("    m_nu ~ m_D^2/M_R ~ phi^(2*%d - %d) = phi^%d" % (n, np_conj, n_seesaw))
        m_nu_pred = M_E * 1e3 * PHI**(n_seesaw)  # in eV if m_e in eV
        # m_e = 0.511 MeV = 511000 eV
        m_nu_pred_eV = 511000 * PHI**(n_seesaw)
        print("    m_nu = m_e * phi^%d = %.4e eV" % (n_seesaw, m_nu_pred_eV))
    else:
        print("    n' = 0, no seesaw")

print()
print("  Problem: n_seesaw is POSITIVE (large mass), not negative (small mass)")
print("  The Galois seesaw gives the WRONG DIRECTION.")
print()

# Try the INVERSE: M_R ~ phi^(-n')
print("  Alternative: M_R ~ phi^(-n') (negative Galois exponent):")
for name, a, b in leptons:
    n = 5*a + 6*b
    np_conj = 5*a - b
    n_seesaw = 2*n + np_conj  # m_D^2/M_R where M_R ~ phi^(-n')
    print("  %s: n_seesaw = 2*%d + %d = %d" % (name, n, np_conj, n_seesaw))
    m_nu_pred_eV = 511000 * PHI**(n_seesaw)
    print("    m_nu = m_e * phi^%d = %.4e eV" % (n_seesaw, m_nu_pred_eV))

print()
print("  Still wrong direction (masses too large).")
print()

# The fundamental problem: seesaw needs M_R >> m_D
# In our framework, all exponents are moderate (0-17 for charged leptons)
# To get m_nu/m_e ~ 10^-7, need phi^(-33) or so
# This requires n_nu ~ -33, which doesn't come from moderate (a,b) values

# Let's try: what (a,b) gives n = -33?
# 5a + 6b = -33
# Simplest: a = -3, b = -3 gives -15-18 = -33
# Or: a = 3, b = -8 gives 15-48 = -33
# Or: a = -9, b = 2 gives -45+12 = -33
# None of these correspond to units (|N(z)| = 1)

# Check: z = a + b*phi, N(z) = a^2 + ab - b^2
print("  Looking for (a,b) with n = 5a+6b ~ -33 and |N(z)|=1:")
for a in range(-15, 15):
    for b in range(-15, 15):
        n = 5*a + 6*b
        norm = a*a + a*b - b*b
        if abs(n + 33) <= 2 and abs(abs(norm) - 1) < 0.01:
            print("    (a,b) = (%d,%d): n=%d, N(z)=%d" % (a, b, n, norm))

print()

# ============================================================
# CONCLUSIONS
# ============================================================
print()
print("=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print()

print("1. MASS SCALE: phi^-33 ~ 7e-8 close to m_nu/m_e ratio")
print("   33 = 3 * 11 = N_gen * (a_1 + b_1)")
print("   But no (a,b) with |N(z)|=1 gives n ~ -33.")
print("   STATUS: PATTERN (numerological)")
print()

print("2. MASS RATIOS: m3/m2 ~ phi^3.7 (not a clean phi power)")
print("   Dm^2_sol/Dm^2_atm = 0.0307 ~ 1/phi^7.2 (not clean)")
print("   STATUS: NEGATIVE (no clean formula found)")
print()

print("3. PMNS ANGLES:")
print("   sin^2(th12) ~ 4/13 (1.2%% off) -- interesting!")
print("   sin^2(th23) ~ maximal (~ 1/2)")
print("   sin^2(th13) ~ 0.022 (no clean match)")
print("   STATUS: PATTERN for th12 only")
print()

print("4. NON-DISTANCE-REGULARITY:")
print("   600-cell has 8 association classes (not 5)")
print("   9 eigenvalues = 8 classes + 1 (identity)")
print("   The 'extra' classes at d=2,3,4 encode finer structure")
print("   that could relate to lepton mixing")
print("   STATUS: DERIVED (mathematical fact), SPECULATIV (physics)")
print()

print("5. GALOIS SEESAW: all attempts give WRONG DIRECTION")
print("   (masses too large, not too small)")
print("   STATUS: NEGATIVE")
print()

print("6. OVERALL: Neutrino masses remain the biggest gap.")
print("   The framework naturally gives MASSLESS neutrinos.")
print("   Getting tiny masses requires a mechanism outside")
print("   the basic (a,b) quantum number structure.")
print()

print("EXP-206 COMPLETE")
