"""
EXP-118: Proton Mass from 600-Cell Framework
=============================================

The proton mass ratio m_p/m_e = 1836.15 is one of the most precisely
known constants. Can it be derived from 600-cell geometry?

APPROACHES:
A. Direct phi-power: m_p/m_e = phi^N for some N?
B. Via LAMBDA_QCD: alpha_s -> LAMBDA_QCD -> m_p (dimensional transmutation)
C. Composite formula: m_p = f(m_e, alpha, alpha_s, phi, geometric constants)
D. Spectral: m_p from eigenvalue ratios or combinatorial factors

KNOWN: alpha_s = 1/(2*phi^3), alpha = 1/137.036
The proton mass comes from QCD confinement, not the Higgs mechanism.
99% of the proton mass is from the strong force binding energy.

Author: Claude (exp118)
Date: 2026-02-08
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1 / 137.036
ALPHA_S = 1 / (2 * PHI**3)

# Experimental values
m_e = 0.51099895e-3   # GeV
m_p = 0.93827208       # GeV
m_p_over_m_e = m_p / m_e  # = 1836.15

print("=" * 80)
print("EXP-118: Proton Mass from 600-Cell Framework")
print("=" * 80)
print()
print(f"  m_p/m_e = {m_p_over_m_e:.4f}")
print(f"  ln(m_p/m_e) / ln(phi) = {np.log(m_p_over_m_e)/np.log(PHI):.6f}")
print()

# ============================================================
# SECTION 1: Direct phi power
# ============================================================

print("-" * 80)
print("SECTION 1: Direct phi powers")
print("-" * 80)
print()

for n in range(14, 18):
    val = PHI**n
    err = abs(val - m_p_over_m_e) / m_p_over_m_e * 100
    print(f"  phi^{n} = {val:.2f}, error = {err:.1f}%")

print()
print(f"  m_p/m_e = phi^{np.log(m_p_over_m_e)/np.log(PHI):.4f}")
print(f"  Closest integer: phi^16 = {PHI**16:.2f} (error {abs(PHI**16 - m_p_over_m_e)/m_p_over_m_e*100:.1f}%)")
print(f"  NOT a clean phi power.")
print()

# ============================================================
# SECTION 2: Dimensional transmutation (LAMBDA_QCD)
# ============================================================

print("-" * 80)
print("SECTION 2: LAMBDA_QCD from alpha_s")
print("-" * 80)
print()

# The QCD scale LAMBDA_QCD is related to alpha_s at scale mu by:
# LAMBDA_QCD = mu * exp(-2*pi / (b_0 * alpha_s(mu)))
# where b_0 = (33 - 2*N_f) / (12*pi) for SU(3) with N_f flavors

# At mu = M_Z = 91.19 GeV, alpha_s = 0.1179
M_Z = 91.1876  # GeV

for N_f in [3, 4, 5, 6]:
    b_0 = (33 - 2 * N_f) / (12 * np.pi)
    Lambda_QCD = M_Z * np.exp(-1 / (2 * b_0 * ALPHA_S))
    print(f"  N_f={N_f}: b_0 = {b_0:.4f}, Lambda_QCD = {Lambda_QCD*1000:.1f} MeV")

print()
# Standard value: Lambda_QCD ~ 200-350 MeV (depends on scheme and N_f)
# Using our alpha_s = 1/(2*phi^3):
N_f = 5  # at M_Z scale
b_0 = (33 - 2 * N_f) / (12 * np.pi)
Lambda_QCD = M_Z * np.exp(-1 / (2 * b_0 * ALPHA_S))
Lambda_QCD_exp = 0.210  # GeV (MSbar, N_f=5)

print(f"  Our prediction: Lambda_QCD(N_f=5) = {Lambda_QCD*1000:.1f} MeV")
print(f"  Experimental: ~210 MeV (MSbar)")
print(f"  Error: {abs(Lambda_QCD - Lambda_QCD_exp)/Lambda_QCD_exp*100:.1f}%")
print()

# Proton mass from LAMBDA_QCD:
# Lattice QCD gives m_p ~ 3.2 * Lambda_QCD (approximately)
# More precisely: m_p ~ Lambda_QCD * exp(C/alpha_s) for some constant
# But the simplest relation: m_p ~ (2-4) * Lambda_QCD

ratio = m_p / Lambda_QCD
print(f"  m_p / Lambda_QCD(ours) = {ratio:.2f}")
print(f"  In terms of phi: {np.log(ratio)/np.log(PHI):.4f}")
print()

# ============================================================
# SECTION 3: Composite formulas
# ============================================================

print("-" * 80)
print("SECTION 3: Composite formulas for m_p/m_e")
print("-" * 80)
print()

target = m_p_over_m_e  # 1836.15

# Test various combinations
candidates = []

# A) m_p/m_e = phi^n / alpha^k
for n in range(10, 20):
    for k_num in range(1, 10):
        for k_den in range(1, 5):
            k = k_num / k_den
            val = PHI**n / ALPHA**k
            if 0.001 < val < 1e8:
                err = abs(val - target) / target * 100
                if err < 1:
                    candidates.append((err, f"phi^{n} / alpha^({k_num}/{k_den})", val))

# B) m_p/m_e = alpha_s^a * phi^n
for a in range(-3, 4):
    for n in range(-5, 30):
        val = ALPHA_S**a * PHI**n
        if 100 < val < 10000:
            err = abs(val - target) / target * 100
            if err < 2:
                candidates.append((err, f"alpha_s^{a} * phi^{n}", val))

# C) m_p/m_e related to spectral constants
# 20*phi^4 = coefficient in alpha equation
coeff_alpha = 20 * PHI**4
for k in range(-2, 8):
    for n in range(-5, 20):
        val = coeff_alpha**k * PHI**n
        if 100 < val < 10000:
            err = abs(val - target) / target * 100
            if err < 2:
                candidates.append((err, f"(20*phi^4)^{k} * phi^{n}", val))

# D) Products of 600-cell constants
# 120 (vertices), 720 (edges), 1200 (faces), 600 (cells)
for a, name_a, val_a in [(120, "V", 120), (720, "E", 720), (1200, "F", 1200), (600, "C", 600)]:
    for n in range(-5, 10):
        val = val_a * PHI**n
        err = abs(val - target) / target * 100
        if err < 2:
            candidates.append((err, f"{name_a}({a}) * phi^{n}", val))

    # Ratios
    for b, name_b, val_b in [(120, "V", 120), (720, "E", 720)]:
        if val_a != val_b:
            for n in range(-5, 10):
                val = (val_a / val_b) * PHI**n
                if 100 < val < 10000:
                    err = abs(val - target) / target * 100
                    if err < 2:
                        candidates.append((err, f"({name_a}/{name_b}) * phi^{n}", val))

# E) m_p/m_e = (some integer)^2 * phi^n
for base in [3, 5, 6, 10, 12, 14, 26, 42]:
    for n in range(-5, 15):
        val = base**2 * PHI**n
        if 100 < val < 10000:
            err = abs(val - target) / target * 100
            if err < 2:
                candidates.append((err, f"{base}^2 * phi^{n}", val))

# F) via alpha and phi together
# m_p/m_e = phi^a * (1/alpha)^b
for a in range(-5, 25):
    for b_num in range(1, 20):
        for b_den in range(1, 8):
            b = b_num / b_den
            val = PHI**a * (1/ALPHA)**b
            if abs(val - target) / target < 0.01:
                err = abs(val - target) / target * 100
                candidates.append((err, f"phi^{a} * (1/alpha)^({b_num}/{b_den})", val))

# Sort by error
candidates.sort()

print("  Best composite formulas (error < 2%):")
print(f"  {'Error':>8s} {'Formula':>40s} {'Value':>12s}")
print(f"  {'-'*64}")
seen_formulas = set()
count = 0
for err, formula, val in candidates:
    if formula not in seen_formulas and count < 20:
        seen_formulas.add(formula)
        print(f"  {err:7.3f}% {formula:>40s} {val:12.4f}")
        count += 1

print()

# ============================================================
# SECTION 4: The 3-quark structure
# ============================================================

print("-" * 80)
print("SECTION 4: Proton as 3-quark composite (uud)")
print("-" * 80)
print()

# Proton = uud. Current quark masses:
m_u = 2.16e-3  # GeV
m_d = 4.67e-3  # GeV

m_quarks = 2 * m_u + m_d
print(f"  Current quark masses: 2*m_u + m_d = {m_quarks*1000:.2f} MeV")
print(f"  This is only {m_quarks/m_p*100:.1f}% of the proton mass!")
print(f"  The remaining {(1 - m_quarks/m_p)*100:.1f}% comes from QCD binding energy.")
print()

# In our framework:
# m_u = m_e * phi^3 * C_up(d_eff, k_eff)
# Can we express m_p in terms of the binding energy?

# The proton mass in terms of alpha_s:
# m_p ~ Lambda_QCD * (1 + c1*alpha_s + c2*alpha_s^2 + ...)
# Or: m_p ~ 3 * constituent_quark_mass
# Constituent quark mass ~ Lambda_QCD ~ 300 MeV -> m_p ~ 900 MeV ~ 940 MeV

# More interesting: m_p/Lambda_QCD in terms of phi?
print(f"  m_p / Lambda_QCD = {m_p/Lambda_QCD:.4f}")
print()

# Try: m_p = Lambda_QCD * phi^n
for n in range(0, 10):
    val = Lambda_QCD * PHI**n
    err = abs(val - m_p) / m_p * 100
    if err < 20:
        print(f"    Lambda_QCD * phi^{n} = {val:.4f} GeV (error {err:.1f}%)")

print()

# ============================================================
# SECTION 5: Koide-like formula for proton
# ============================================================

print("-" * 80)
print("SECTION 5: Algebraic relations")
print("-" * 80)
print()

# Test: m_p/m_e in terms of phi and small integers
# Key identity: phi^2 = phi + 1
# Try: m_p/m_e = A * phi^B + C for small integers

print("  Testing m_p/m_e = A * phi^B + C:")
best = []
for A in range(-50, 51):
    if A == 0: continue
    for B in range(0, 20):
        for C in range(-50, 51):
            val = A * PHI**B + C
            if val > 0:
                err = abs(val - target) / target * 100
                if err < 0.5:
                    best.append((err, A, B, C, val))

best.sort()
for err, A, B, C, val in best[:10]:
    print(f"    {A}*phi^{B} + {C} = {val:.4f} (error {err:.3f}%)")

print()

# Test: m_p/m_e = phi^a * N where N is a "nice" number
print("  Testing m_p/m_e = phi^a * N:")
best2 = []
for a_num in range(0, 40):
    for a_den in range(1, 5):
        a = a_num / a_den
        N = target / PHI**a
        N_int = round(N)
        if N_int > 0 and abs(N - N_int) / N < 0.005:
            err = abs(N_int * PHI**a - target) / target * 100
            best2.append((err, a_num, a_den, N_int))

best2.sort()
for err, a_num, a_den, N in best2[:10]:
    a_str = f"{a_num}/{a_den}" if a_den > 1 else str(a_num)
    print(f"    phi^({a_str}) * {N} = {N * PHI**(a_num/a_den):.4f} (error {err:.3f}%)")

print()

# ============================================================
# SECTION 6: 600-cell combinatorial formulas
# ============================================================

print("-" * 80)
print("SECTION 6: Combinatorial formulas from 600-cell")
print("-" * 80)
print()

# 600-cell numbers: V=120, E=720, F=1200, C=600, |H4|=14400
# Also: 5 (diameter), 6 (decagons), 12 (degree), 26 (phi-sector)

target = m_p_over_m_e

# Test products and ratios of 600-cell constants
constants_600 = {
    'V': 120, 'E': 720, 'F': 1200, 'C': 600,
    'deg': 12, 'diam': 5, 'dec': 6, 'phi_sec': 26,
    '|H4|': 14400, 'a1': 5, 'b1': 6
}

print("  Testing A * B / C for 600-cell constants:")
best3 = []
items = list(constants_600.items())
for i, (n1, v1) in enumerate(items):
    for j, (n2, v2) in enumerate(items):
        for k, (n3, v3) in enumerate(items):
            if v3 == 0: continue
            val = v1 * v2 / v3
            if 100 < val < 10000:
                err = abs(val - target) / target * 100
                if err < 5:
                    best3.append((err, f"{n1}*{n2}/{n3}", val))

            # Also with phi multiplier
            for p in range(-3, 4):
                val2 = v1 * v2 / v3 * PHI**p
                if 100 < val2 < 10000:
                    err2 = abs(val2 - target) / target * 100
                    if err2 < 2:
                        best3.append((err2, f"{n1}*{n2}/{n3}*phi^{p}", val2))

best3.sort()
seen = set()
for err, formula, val in best3[:15]:
    if formula not in seen:
        seen.add(formula)
        print(f"    {formula:30s} = {val:10.2f} (error {err:.2f}%)")

print()

# ============================================================
# SECTION 7: Dimensional transmutation chain
# ============================================================

print("-" * 80)
print("SECTION 7: Full derivation chain alpha_s -> Lambda_QCD -> m_p")
print("-" * 80)
print()

# The only rigorous way: alpha_s determines Lambda_QCD, which determines m_p
# m_p / Lambda_QCD ~ 4.5 (from lattice QCD)
# Can we get this factor from geometry?

# Our alpha_s = 1/(2*phi^3) = 0.1180
# Lambda_QCD = M_Z * exp(-1/(2*b_0*alpha_s))

# Then: m_p = Lambda_QCD * R where R ~ 4.5 is the "QCD ratio"
# The question: is R expressible in terms of phi?

R = m_p / Lambda_QCD
print(f"  m_p / Lambda_QCD = {R:.4f}")
print(f"  In phi: phi^{np.log(R)/np.log(PHI):.4f}")
print()

# Test: is R close to some phi expression?
for a in range(0, 8):
    for b in range(-3, 4):
        val = PHI**a + b
        if val > 0:
            err = abs(val - R) / R * 100
            if err < 5:
                print(f"    phi^{a} + {b} = {val:.4f} (error {err:.1f}%)")

# Also test simple fractions times phi
for num in range(1, 20):
    for den in range(1, 10):
        for p in range(-3, 6):
            val = (num/den) * PHI**p
            if abs(val - R) / R < 0.02:
                err = abs(val - R) / R * 100
                print(f"    ({num}/{den})*phi^{p} = {val:.4f} (error {err:.2f}%)")

print()

# ============================================================
# SECTION 8: Summary
# ============================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print(f"  m_p/m_e = {m_p_over_m_e:.4f}")
print(f"  = phi^{np.log(m_p_over_m_e)/np.log(PHI):.4f} (not a clean power)")
print()

if candidates:
    best_err, best_formula, best_val = candidates[0]
    print(f"  Best composite formula: {best_formula}")
    print(f"  Value: {best_val:.4f}, error: {best_err:.3f}%")
    print()

print(f"  Dimensional transmutation chain:")
print(f"    alpha_s = 1/(2*phi^3) = {ALPHA_S:.4f}")
print(f"    Lambda_QCD = M_Z * exp(-1/(2*b_0*alpha_s)) = {Lambda_QCD*1000:.1f} MeV")
print(f"    m_p / Lambda_QCD = {R:.2f}")
print()

print("  CLASSIFICATION: The proton mass is a QCD COMPOSITE quantity.")
print("  It cannot be simply expressed as phi^n * m_e.")
print("  The chain alpha_s -> Lambda_QCD -> m_p is the correct approach,")
print("  but the last step (Lambda_QCD -> m_p) requires non-perturbative QCD")
print("  (lattice calculations) that cannot be reduced to simple algebra.")
print()
