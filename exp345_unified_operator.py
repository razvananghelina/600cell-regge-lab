"""
exp345: Unified Mass Correction Operator - Systematic Search
=============================================================
Goal: Can the 7 mass correction cases (Eq. normlog) be unified into fewer formulas?

Current 7 cases:
  1. e: delta=0 (z=0)
  2. mu,tau: delta = c_ell * sign(z') * |z'|^(3/4)
  3. u: delta=0 (|N|=1, T3=+1/2)
  4. c,t: delta = +C * ln|N|  (T3=+1/2, |N|>1)
  5. d: delta = -2/a1  (T3=-1/2, b=0, rational)
  6. s: delta = -3C/phi^2  (T3=-1/2, |N|=1, unit)
  7. b: delta = -C*ln|N|/phi  (T3=-1/2, |N|>1)

Key insight to test: cases 4,7 unify via phi^(2T3-1).
"""
import numpy as np

A1 = 5
B1 = 6
N_GEN = 3
PHI = (1 + np.sqrt(5)) / 2
PHIp = (1 - np.sqrt(5)) / 2
C = 4 / (A1**2 + 1)  # 2/13
ALPHA = 7.2973525693e-3
C_ELL = PHI / (2 * A1)  # phi/10

m_e = 0.51099895  # MeV
pdg = dict(e=0.51099895, mu=105.6583755, tau=1776.86,
           u=2.16, c=1270.0, t=172760.0,
           d=4.67, s=93.4, b=4180.0)

fermions = ["e","mu","tau","u","c","t","d","s","b"]
a_vals = dict(e=0,mu=1,tau=1,u=3,c=2,t=4,d=1,s=1,b=-1)
b_vals = dict(e=0,mu=1,tau=2,u=-2,c=1,t=1,d=0,s=1,b=4)
T3_vals = dict(e=0,mu=0,tau=0,u=0.5,c=0.5,t=0.5,d=-0.5,s=-0.5,b=-0.5)

# Compute algebraic quantities
z, zp, Nz, n_bare = {}, {}, {}, {}
for f in fermions:
    a, b = a_vals[f], b_vals[f]
    z[f] = a + b * PHI
    zp[f] = a + b * PHIp
    Nz[f] = a**2 + a*b - b**2
    n_bare[f] = 5*a + 6*b

# Experimental delta
n_eff = {f: np.log(pdg[f]/m_e)/np.log(PHI) for f in fermions}
delta_exp = {f: n_eff[f] - n_bare[f] for f in fermions}

# Unit decomposition: z = sign * phi^k for units
def phi_rank(f):
    """For units |N|=1, find k such that |z| = phi^k"""
    if z[f] == 0: return None  # electron
    if abs(Nz[f]) > 1: return None  # not a unit
    if abs(z[f]) < 1e-10: return 0
    k = np.log(abs(z[f])) / np.log(PHI)
    k_round = round(k)
    if abs(k - k_round) < 0.001:
        return k_round
    return None

def mass_pred(f, delta):
    return m_e * PHI**(n_bare[f] + delta)

def pct_err(f, delta):
    return (mass_pred(f, delta) - pdg[f]) / pdg[f] * 100

def rms(errs):
    return np.sqrt(np.mean(np.array(errs)**2))

SEP = "=" * 78

# ============================================================
# TABLE 0: Algebraic structure
# ============================================================
print(SEP)
print("EXP345: UNIFIED MASS CORRECTION OPERATOR")
print(SEP)
print()
print("TABLE 0: Algebraic structure of fermions in Z[phi]")
print("-" * 78)
print("f     (a ,b )       z      z'     N(z)  |N|    T3  phi_k  type")
print("-" * 78)
for f in fermions:
    a, b = a_vals[f], b_vals[f]
    pk = phi_rank(f)
    if pk is not None:
        pks = "%+d" % pk
    else:
        pks = " -"
    atype = "ZERO" if z[f]==0 else ("unit" if abs(Nz[f])<=1 else "prime")
    print("%-5s (%2d,%2d) %+7.3f %+8.4f %+8d %4d %+5.1f  %4s   %s" %
          (f, a, b, z[f], zp[f], Nz[f], abs(Nz[f]), T3_vals[f], pks, atype))

# Check z_b = phi * sigma(z_t)
print()
print("KEY ALGEBRAIC RELATION:")
print("  z_b = phi * sigma(z_t) = phi * (4+phi') = phi*(%.4f) = %.4f" %
      (zp["t"], PHI * zp["t"]))
print("  Actual z_b = %.4f  CHECK: %s" %
      (z["b"], "MATCH" if abs(PHI*zp["t"] - z["b"]) < 1e-10 else "FAIL"))

# ============================================================
# PARTIAL UNIFICATION: Up/Down quarks with |N|>1
# ============================================================
print()
print(SEP)
print("RESULT 1: Partial unification of quarks with |N| > 1")
print(SEP)
print()
print("  delta_q = 2*T3 * C * ln|N(z)| * phi^(T3 - 1/2)")
print()
print("  T3 = +1/2 (up):   +C * ln|N| * phi^0   = +C * ln|N|      [c, t]")
print("  T3 = -1/2 (down): -C * ln|N| * phi^(-1) = -C * ln|N|/phi  [b]")
print()

for f in ["c", "t", "b"]:
    T3 = T3_vals[f]
    delta_unified = 2*T3 * C * np.log(abs(Nz[f])) * PHI**(T3-0.5)
    err = pct_err(f, delta_unified)
    print("  %s: delta_unified = %+.4f, delta_known = %+.4f, mass err = %+.3f%%" %
          (f, delta_unified, delta_exp[f], err))

# Also check u (|N|=1, T3=+1/2): formula gives 0 automatically
delta_u = np.sign(0.5) * C * np.log(max(abs(Nz["u"]), 1e-100)) * PHI**(2*0.5-1)
# ln(1) = 0, so delta = 0
print("  u: |N|=1, ln|N|=0, delta_unified = 0 (automatic). mass err = %+.3f%%" %
      pct_err("u", 0))

# ============================================================
# THE STUBBORN CASES: d, s (unit down-quarks)
# ============================================================
print()
print(SEP)
print("RESULT 2: Unit down-quark residuals (the hard problem)")
print(SEP)
print()
print("d: z=1 (rational, phi^0), delta = -2/a1 = %.4f" % (-2/A1))
print("s: z=phi^2 (unit), delta = -3C/phi^2 = %.4f" % (-3*C/PHI**2))
print()

# Test various unified formulas for d and s
print("Can d and s be written as f(phi_rank, T3)?")
print()

# d: k=0, delta=-0.4
# s: k=2, delta=-0.1763
k_d, k_s = 0, 2
delta_d, delta_s = -2/A1, -3*C/PHI**2

# Test: delta = -A * phi^(-B*k) for down units
# d: -A * 1 = -0.4, so A = 0.4 = 2/a1
# s: -(2/a1) * phi^(-B*2) = -0.1763
# phi^(-2B) = 0.1763/0.4 = 0.4408
# -2B * ln(phi) = ln(0.4408)
B_test = -np.log(0.4408) / (2 * np.log(PHI))
print("  Test delta = -(2/a1) * phi^(-B*k):")
print("  B = %.4f (need clean fraction)" % B_test)
# B = 0.8512, not clean

# Test: delta = -A * |z'|^(3/4) (same exponent as leptons!)
# d: z'=1, |z'|^(3/4) = 1. -A*1 = -0.4, A = 0.4
# s: z'=1/phi^2, |z'|^(3/4) = phi^(-3/2) = 0.4835. -0.4*0.4835 = -0.1934. Need -0.1763.
delta_s_test34 = -2/A1 * abs(zp["s"])**(3.0/4)
print()
print("  Test delta = -(2/a1) * |z'|^(3/4) (lepton exponent):")
print("  d: pred=%.4f, need=%.4f (MATCH)" % (-2/A1 * 1.0, delta_d))
print("  s: pred=%.4f, need=%.4f (err=%.1f%%)" %
      (delta_s_test34, delta_s, (delta_s_test34-delta_s)/delta_s*100))

# Test: delta = -A * |z'|^k for what k?
# d: -A*1 = -0.4 always (regardless of k)
# s: -0.4 * (1/phi^2)^k = -0.1763 => (0.382)^k = 0.4408
k_opt = np.log(0.4408) / np.log(abs(zp["s"]))
print()
print("  Optimal exponent for quark units: k = %.4f" % k_opt)
print("  Compare: lepton k = 3/4 = 0.7500, quark k_opt = %.4f" % k_opt)
print("  Ratio: k_opt/k_lep = %.4f" % (k_opt/0.75))

# Test k = 5/6
k56 = 5/6
delta_s_56 = -2/A1 * abs(zp["s"])**k56
print()
print("  Test k = 5/6 = 0.8333:")
print("  s: pred=%.4f, need=%.4f (err=%.1f%%)" %
      (delta_s_56, delta_s, (delta_s_56-delta_s)/delta_s*100))

# Test k = 1 (linear in z')
delta_s_1 = -2/A1 * abs(zp["s"])
print("  Test k = 1:")
print("  s: pred=%.4f, need=%.4f (err=%.1f%%)" %
      (delta_s_1, delta_s, (delta_s_1-delta_s)/delta_s*100))

# ============================================================
# CANDIDATE UNIFIED FORMULAS
# ============================================================
print()
print(SEP)
print("RESULT 3: Testing unified formula candidates")
print(SEP)

def delta_known(f):
    """The 7-case formula from the paper"""
    T3 = T3_vals[f]
    a, b = a_vals[f], b_vals[f]
    absN = abs(Nz[f])
    if T3 == 0:
        if z[f] == 0: return 0.0
        return C_ELL * np.sign(zp[f]) * abs(zp[f])**(3.0/4)
    elif T3 > 0:
        return C * np.log(absN) if absN > 1 else 0.0
    else:
        if b == 0: return -2.0/A1
        elif absN > 1: return -C * np.log(absN) / PHI
        else: return -3*C/PHI**2

def test_formula(name, delta_fn):
    errs = []
    print()
    print("  --- %s ---" % name)
    for f in fermions:
        d = delta_fn(f)
        e = pct_err(f, d)
        errs.append(e)
        print("  %-5s: delta=%+7.4f  m_pred=%10.3f  m_pdg=%10.3f  err=%+7.3f%%" %
              (f, d, mass_pred(f, d), pdg[f], e))
    quarks = [e for f,e in zip(fermions,errs) if T3_vals[f]!=0]
    leptons_e = [e for f,e in zip(fermions,errs) if T3_vals[f]==0]
    print("  RMS: all=%.4f%%, quarks=%.4f%%, leptons=%.4f%%" %
          (rms(errs), rms(quarks), rms(leptons_e)))
    return rms(errs)

# Formula A: Original 7-case (baseline)
rms_A = test_formula("Formula A: Original 7-case (baseline)", delta_known)

# Formula B: 3-case unified (lose d, s accuracy)
def delta_B(f):
    T3 = T3_vals[f]
    absN = abs(Nz[f])
    if T3 == 0:
        if z[f] == 0: return 0.0
        return C_ELL * np.sign(zp[f]) * abs(zp[f])**(3.0/4)
    else:
        if absN <= 1: return 0.0  # gives 0 for u, d, s
        return 2*T3 * C * np.log(absN) * PHI**(T3-0.5)

rms_B = test_formula("Formula B: 3-case (leptons + unified quarks, units->0)", delta_B)

# Formula C: Use |z'|^(3/4) for ALL sectors, different coefficient per T3
def delta_C(f):
    T3 = T3_vals[f]
    if z[f] == 0: return 0.0
    if abs(zp[f]) < 1e-10: return 0.0
    if T3 == 0:
        return C_ELL * np.sign(zp[f]) * abs(zp[f])**(3.0/4)
    elif T3 > 0:
        # Need to match: +C*ln|N| for c,t. What coeff with |z'|^(3/4)?
        return C * np.sign(zp[f]) * abs(zp[f])**(3.0/4) * PHI**(T3-0.5)
    else:
        return -2/A1 * np.sign(zp[f]) * abs(zp[f])**(3.0/4)

rms_C = test_formula("Formula C: |z'|^(3/4) for ALL, T3-dependent coeff", delta_C)

# Formula D: ln|N| for quarks, residual for units via |z'|^(3/4)
def delta_D(f):
    T3 = T3_vals[f]
    absN = abs(Nz[f])
    if T3 == 0:
        if z[f] == 0: return 0.0
        return C_ELL * np.sign(zp[f]) * abs(zp[f])**(3.0/4)
    else:
        # Main term (works for |N|>1, gives 0 for |N|=1)
        main = 0.0
        if absN > 1:
            main = 2*T3 * C * np.log(absN) * PHI**(T3-0.5)
        # Residual for units (T3<0 only, since u gets 0 correctly)
        resid = 0.0
        if absN <= 1 and T3 < 0:
            resid = -2/A1 * abs(zp[f])**(3.0/4) * np.sign(zp[f])
        return main + resid

rms_D = test_formula("Formula D: ln|N| quarks + |z'|^(3/4) residual for unit-downs", delta_D)

# Formula E: Optimal power-law for down-quark units
k_opt_val = k_opt  # ~0.851
def delta_E(f):
    T3 = T3_vals[f]
    absN = abs(Nz[f])
    if T3 == 0:
        if z[f] == 0: return 0.0
        return C_ELL * np.sign(zp[f]) * abs(zp[f])**(3.0/4)
    else:
        if absN > 1:
            return 2*T3 * C * np.log(absN) * PHI**(T3-0.5)
        elif T3 > 0:
            return 0.0
        else:
            return -2/A1 * np.sign(zp[f]) * abs(zp[f])**k_opt_val

rms_E = test_formula("Formula E: Optimal k=%.4f for unit-down residual" % k_opt_val, delta_E)

# Formula F: Try c_ell = 2/a1 universally (same coeff for leptons and unit-downs)
# with k=3/4 everywhere
def delta_F(f):
    T3 = T3_vals[f]
    absN = abs(Nz[f])
    if z[f] == 0: return 0.0
    if T3 == 0:
        return (2/A1) * np.sign(zp[f]) * abs(zp[f])**(3.0/4)
    else:
        if absN > 1:
            return 2*T3 * C * np.log(absN) * PHI**(T3-0.5)
        elif T3 > 0:
            return 0.0
        else:
            return -(2/A1) * abs(zp[f])**(3.0/4) * np.sign(zp[f])

rms_F = test_formula("Formula F: c_ell=2/a1=0.4 for leptons too, k=3/4", delta_F)

# ============================================================
# THE DEEP STRUCTURE: N = z * z' decomposition
# ============================================================
print()
print(SEP)
print("RESULT 4: N(z) = z * z' structure")
print(SEP)
print()
print("For non-units: ln|N| = ln|z| + ln|z'|")
print("For units: ln|N| = 0, but ln|z| = -ln|z'| != 0")
print()
print("%-5s  ln|z|    ln|z'|   ln|N|   delta_exp" % "f")
print("-" * 50)
for f in fermions:
    if z[f] == 0:
        print("%-5s  -inf     -inf     -inf    %+.4f" % (f, delta_exp[f]))
    else:
        lz = np.log(abs(z[f]))
        lzp = np.log(abs(zp[f]))
        lN = np.log(max(abs(Nz[f]), 1e-100))
        print("%-5s  %+7.4f  %+7.4f  %+7.4f  %+.4f" % (f, lz, lzp, lN, delta_exp[f]))

# ============================================================
# FORMULA G: The "Galois Archimedean" formula
# ============================================================
# Idea: delta = coeff(T3) * ln|z'| for quarks, c_ell * sign(z')*|z'|^(3/4) for leptons
print()
print(SEP)
print("RESULT 5: Galois-archimedean formula: delta_q ~ ln|z'|")
print(SEP)
print()
print("Test: delta_quark = A(T3) * ln|z'|")
print()
for f in ["u","c","t","d","s","b"]:
    lzp_val = np.log(abs(zp[f]))
    ratio = delta_exp[f] / lzp_val if abs(lzp_val) > 0.01 else float('inf')
    print("  %s: ln|z'|=%+7.4f, delta_exp=%+7.4f, ratio=%+7.4f" %
          (f, lzp_val, delta_exp[f], ratio))

# ============================================================
# BEST PARTIAL UNIFICATION
# ============================================================
print()
print(SEP)
print("SUMMARY: Achievable unification")
print(SEP)
print()
print("The 7 cases REDUCE to 4 structural types:")
print()
print("TYPE 1 (trivial): delta = 0")
print("  Applies to: e (z=0), u (|N|=1, T3>0)")
print("  Condition: z=0 OR (unit AND T3>0)")
print()
print("TYPE 2 (Galois power): delta = c_ell * sign(z') * |z'|^(3/4)")
print("  Applies to: mu, tau")
print("  Condition: T3 = 0, z != 0")
print()
print("TYPE 3 (norm-log): delta = 2*T3 * C * ln|N(z)| * phi^(T3-1/2)")
print("  Applies to: c, t, b")
print("  Condition: |N(z)| > 1 (non-unit quarks)")
print("  This is NEW: unifies up/down with |N|>1!")
print()
print("TYPE 4 (unit-down residual): special")
print("  Applies to: d (-2/a1), s (-3C/phi^2)")
print("  Condition: |N(z)| = 1, T3 < 0")
print("  d: -2/a1 = -(a1-N_gen)/a1 [rational sector, b=0]")
print("  s: -N_gen*C/phi^2 = -2*sin^2(tW)/phi^2 [unit sector]")
print()

# Can Type 4 be absorbed?
print("TYPE 4 alternatives tested:")
print("  delta = -(2/a1)*|z'|^(3/4): d EXACT, s err=9.7%% [k=3/4 from leptons]")
print("  delta = -(2/a1)*|z'|^(%.3f): d EXACT, s EXACT [optimal but ugly]" % k_opt_val)
print("  delta = -(2/a1)*|z'|:        d EXACT, s err=13.3%% [linear]")
print()
print("CONCLUSION: 7 cases -> 3 formulas + 1 residual")
print("  The unified quark formula delta = 2*T3*C*ln|N|*phi^(T3-1/2)")
print("  covers c,t,b,u (4 of 6 quarks). The d,s residuals resist")
print("  unification due to their unit-norm algebraic nature.")
print()

# ============================================================
# BONUS: The spectral interpretation
# ============================================================
print(SEP)
print("BONUS: Connection to spectral gap")
print(SEP)
print()
print("The Cayley Laplacian eigenvalue ratio:")
print("  lambda(3')/lambda(3) = (8+4*phi)/(12-4*phi) = %.4f" %
      ((8+4*PHI)/(12-4*PHI)))
print("  phi^2 = %.4f" % PHI**2)
print("  MATCH: lambda(3')/lambda(3) = phi^2 (EXACT)" )
print()
print("This is the SAME phi^2 that appears in L(3')/L(3) = phi^2")
print("(Cayley eigenvalue product giving alpha, Section 5 of paper)")
print()
print("The spectral gap 4*sqrt(5) = %.4f separates the 3 and 3' sectors." %
      (4*np.sqrt(5)))
print("The mass corrections operate WITHIN each sector (norm-log)")
print("while the bare hierarchy (phi^n) operates BETWEEN sectors.")

# Final RMS comparison
print()
print(SEP)
print("RMS COMPARISON")
print(SEP)
print("  A (original 7-case):      %.4f%%" % rms_A)
print("  B (3-case, units->0):     %.4f%%" % rms_B)
print("  C (|z'|^3/4 all):         %.4f%%" % rms_C)
print("  D (ln|N| + |z'|^3/4 res): %.4f%%" % rms_D)
print("  E (optimal k unit-down):  %.4f%%" % rms_E)
print("  F (c_ell=2/a1):           %.4f%%" % rms_F)
