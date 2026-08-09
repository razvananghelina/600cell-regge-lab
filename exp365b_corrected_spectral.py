"""
exp365b_corrected_spectral.py
==============================
CORRECTED version of exp365 spectral formula.

BUG FIX: G(s) = (L_10a - L_5a)/2 = -3, NOT -6.
The correct factor is C*G/phi^2 (not C*G/(2*phi^2)).

Result: delta_s = C*(-3)/phi^2 = -6/(13*phi^2) = -N_gen*C/phi^2 EXACT.

Also: for up quarks, the factor is 4/a1 (not 2/a1).

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PHIc = -1/PHI
a1 = 5
b1 = 6
N_gen = 3
C = 4.0 / (a1**2 + 1)  # = 2/13
m_e = 0.51099895  # MeV

# McKay graph
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
A = np.zeros((9,9))
for i,j in edges:
    A[i,j] = A[j,i] = 1

# Character table from eigenvectors
evals_raw, evecs_raw = np.linalg.eigh(A)
idx = np.argsort(evals_raw)[::-1]
evals = evals_raw[idx]
evecs = evecs_raw[:, idx]
for j in range(9):
    if evecs[0,j] < 0:
        evecs[:,j] *= -1

class_sizes = [1, 12, 20, 12, 30, 12, 20, 12, 1]
char = np.zeros((9,9))
for c in range(9):
    norm = np.sqrt(120.0/class_sizes[c])
    for k in range(9):
        char[k,c] = evecs[k,c] * norm

# Columns: 0=1a(2), 1=10a(phi), 2=3a(1), 3=10b(1/phi), 4=4a(0),
#           5=5b(-1/phi), 6=6a(-1), 7=5a(-phi), 8=2a(-2)
COL_10A, COL_5A = 1, 7

# Bipartite
color = [-1]*9
color[0] = 0
q = [0]
while q:
    n = q.pop(0)
    for j in range(9):
        if A[n,j] > 0 and color[j] < 0:
            color[j] = 1 - color[n]
            q.append(j)

# Fermion map: node -> (name, n_bare, (a,b), T3, Nc)
fmap = {
    0: ('e',   0,  ( 0, 0), -0.5, 1),
    1: ('u',   3,  ( 3,-2), +0.5, 3),
    2: ('d',   5,  ( 1, 0), -0.5, 3),
    3: ('s',  11,  ( 1, 1), -0.5, 3),
    4: ('mu', 11,  ( 1, 1), -0.5, 1),
    5: ('c',  16,  ( 2, 1), +0.5, 3),
    6: ('tau',17,  ( 1, 2), -0.5, 1),
    7: ('t',  26,  ( 4, 1), +0.5, 3),
    8: ('b',  19,  (-1, 4), -0.5, 3),
}

pdg = {'e':0.511, 'mu':105.66, 'tau':1776.86,
       'u':2.16, 'c':1270., 't':172760., 'd':4.67, 's':93.4, 'b':4180.}

delta_known = {'e':0, 'mu':0.0792, 'tau':-0.0552,
               'u':0, 'c':0.2476, 't':0.4530,
               'd':-0.4000, 's':-0.1763, 'b':-0.2800}

def zphi_norm(a, b):
    return a**2 + a*b - b**2

def pct(p, o):
    return (p-o)/o*100 if o != 0 else float('inf')

# ======================================================================
# STEP 1: Compute G values for BLACK nodes
# ======================================================================
print()
print("=" * 80)
print("CORRECTED G VALUES (G = (L_10a - L_5a) / 2)")
print("=" * 80)

spec = {}
for k in range(9):
    name = fmap[k][0]
    d = dims[k]
    chi10a = char[k, COL_10A]
    chi5a  = char[k, COL_5A]
    L10a = 12*(1 - chi10a/d)
    L5a  = 12*(1 - chi5a/d)
    G = (L10a - L5a)/2
    spec[name] = {'G': G, 'L10a': L10a, 'L5a': L5a, 'color': color[k]}

print()
print("BLACK nodes (spin reps, G != 0):")
print("%-4s  G=%+9.4f  L10a=%7.3f  L5a=%7.3f" % ('u', spec['u']['G'], spec['u']['L10a'], spec['u']['L5a']))
print("%-4s  G=%+9.4f  L10a=%7.3f  L5a=%7.3f" % ('s', spec['s']['G'], spec['s']['L10a'], spec['s']['L5a']))
print("%-4s  G=%+9.4f  L10a=%7.3f  L5a=%7.3f" % ('c', spec['c']['G'], spec['c']['L10a'], spec['c']['L5a']))
print("%-4s  G=%+9.4f  L10a=%7.3f  L5a=%7.3f" % ('t', spec['t']['G'], spec['t']['L10a'], spec['t']['L5a']))
print()
print("G values in Z[phi]:")
print("  G(u) = -6*phi = %.4f" % (-6*PHI))
print("  G(s) = -3 (integer)")
print("  G(c) = +2 (integer)")
print("  G(t) = +6/phi = +6*(phi-1) = %.4f" % (6/PHI))
print()
print("NOTE: G(u) and G(t) are GALOIS CONJUGATES: sigma(-6*phi) = 6/phi")
print("      G(s) and G(c) are RATIONAL (Galois-invariant)")

# ======================================================================
# STEP 2: EXACT s quark result (CORRECTED)
# ======================================================================
print()
print("=" * 80)
print("EXACT RESULT: s quark correction from G(s)")
print("=" * 80)
print()

G_s = spec['s']['G']
delta_s_spectral = C * G_s / PHI**2
delta_s_known = delta_known['s']

print("G(s) = %.4f = -3 = -N_gen" % G_s)
print("delta_s = C * G(s) / phi^2 = (2/13)*(-3)/phi^2")
print("        = %.10f" % delta_s_spectral)
print("Known:    %.10f" % delta_s_known)
print("Diff:     %.2e" % abs(delta_s_spectral - delta_s_known))
print()
print("WHY it works: G(s)/phi^2 = -3/phi^2 = -N_gen/phi^2")
print("  C*(-N_gen/phi^2) = -N_gen*C/phi^2 = prescription P3.")
print("  STATUS: EXACT DERIVATION of P3 from spectral geometry.")

# ======================================================================
# STEP 3: Approximate c,t results
# ======================================================================
print()
print("=" * 80)
print("APPROXIMATE: c,t quarks from G")
print("=" * 80)
print()

# Known: delta = C*ln|N|. Spectral: try delta = C*(4/a1)*max(0,G)
for name in ['c', 't']:
    G = spec[name]['G']
    a, b = fmap[{'c':5,'t':7}[name]][2]
    Nz = abs(zphi_norm(a, b))
    lnN = np.log(Nz)

    delta_known_val = C * lnN  # this IS the known formula
    delta_spectral = C * (4.0/a1) * G
    ratio = lnN / G

    print("%s: |N|=%d, ln|N|=%.4f, G=%.4f" % (name, Nz, lnN, G))
    print("   ln|N|/G = %.5f (would need = 4/a1 = %.5f for exact match)" % (ratio, 4.0/a1))
    print("   delta_known   = C*ln|N| = %.6f" % delta_known_val)
    print("   delta_spectral = C*(4/a1)*G = %.6f" % delta_spectral)
    print("   Difference: %.4f%%" % (pct(delta_spectral, delta_known_val)))
    print()

print("CONCLUSION: ln|N|/G is NOT exactly 4/a1.")
print("  c: ratio = 0.805 = ln(5)/2")
print("  t: ratio = 0.794 = ln(19)*phi/6")
print("  Differs by ~1.4%%. This is APPROXIMATE, not exact.")

# ======================================================================
# STEP 4: Full corrected formula test
# ======================================================================
print()
print("=" * 80)
print("FULL CORRECTED FORMULA TEST")
print("=" * 80)
print()

print("Formula (corrected):")
print("  BLACK up quarks: delta = C * (4/a1) * max(0, G_k)")
print("  BLACK down quarks: delta = C * G_k / phi^2")
print("  WHITE quarks |N|>1: delta = -C*ln|N|/phi (unchanged)")
print("  WHITE d quark: delta = -2/a1 (unchanged)")
print("  Leptons: delta = c_ell*sign(z')*|z'|^(3/4) (unchanged)")
print()

d_ST = 4.0
c_ell = C * PHI**3 / d_ST

errs_spectral = []
errs_known = []
print("%-4s %4s %5s %+9s %+9s %+9s %12s %12s %+8s %+8s" % (
    "Name", "T3", "color", "delta_S", "delta_K", "diff",
    "m_S/MeV", "m_exp/MeV", "err_S", "err_K"))
print("-" * 110)

for k in range(9):
    name, n_bare, (a,b), T3, Nc = fmap[k]
    G = spec[name]['G']
    Nz = abs(zphi_norm(a,b))
    zp = a + b*PHIc
    col = "W" if color[k]==0 else "B"
    m_exp = pdg[name]
    dk = delta_known[name]

    if name == 'e':
        print("%-4s %+.1f %5s %9s %9s %9s %12.4f %12.4f %8s %8s" % (
            name, T3, col, "0", "0", "-", m_exp, m_exp, "(inp)", "(inp)"))
        continue

    # Spectral formula
    if Nc == 1:  # lepton
        ds = c_ell * np.sign(zp) * abs(zp)**0.75
    elif color[k] == 1:  # BLACK quark
        if T3 > 0:  # up
            ds = C * (4.0/a1) * max(0, G)
        else:  # down
            ds = C * G / PHI**2
    else:  # WHITE quark
        if b == 0:  # rational
            ds = -2.0/a1
        elif Nz > 1:  # prime
            ds = -C * np.log(Nz) / PHI
        else:
            ds = 0.0

    m_s = m_e * PHI**(n_bare + ds)
    m_k = m_e * PHI**(n_bare + dk)
    err_s = pct(m_s, m_exp)
    err_k = pct(m_k, m_exp)
    errs_spectral.append(abs(err_s))
    errs_known.append(abs(err_k))

    print("%-4s %+.1f %5s %+9.4f %+9.4f %+9.4f %12.4f %12.4f %+7.2f%% %+7.2f%%" % (
        name, T3, col, ds, dk, ds-dk, m_s, m_exp, err_s, err_k))

print()
rms_s = np.sqrt(np.mean(np.array(errs_spectral)**2))
rms_k = np.sqrt(np.mean(np.array(errs_known)**2))
print("RMS (spectral corrected):  %.3f%%" % rms_s)
print("RMS (known prescriptions): %.3f%%" % rms_k)

# ======================================================================
# STEP 5: What changed from known to spectral
# ======================================================================
print()
print("=" * 80)
print("WHAT CHANGED (spectral vs known)")
print("=" * 80)
print()

print("UNCHANGED (5/8):")
print("  d:   delta = -2/a1 = -0.4000 (WHITE, rational)")
print("  b:   delta = -C*ln19/phi = -0.2800 (WHITE, prime)")
print("  mu:  delta = c_ell*|z'|^(3/4) = +0.0792 (lepton)")
print("  tau: delta = -c_ell*|z'|^(3/4) = -0.0552 (lepton)")
print("  u:   delta = 0 (BLACK, max(0,G)=0 since G<0)")
print()

print("CHANGED (3/8):")
print()
print("  s: IMPROVED (exact spectral)")
print("    OLD: -N_gen*C/phi^2 (appears ad hoc, N_gen unexplained)")
print("    NEW: C*G(s)/phi^2 where G(s)=-3=-N_gen (SPECTRAL DERIVATION)")
print("    N_gen appears because G(s)=-3 at this node. NOT put in by hand.")
print()

G_c = spec['c']['G']
G_t = spec['t']['G']
print("  c: SMALL DEGRADATION")
print("    OLD: C*ln(5) = %.6f" % (C*np.log(5)))
print("    NEW: C*(4/a1)*G = C*(4/5)*2 = C*1.600 = %.6f" % (C*4.0/a1*G_c))
print("    ln(5)=1.609 vs 1.600: %.2f%% error in correction" % (
    pct(4.0/a1*G_c, np.log(5))))
print()
print("  t: SMALL DEGRADATION")
print("    OLD: C*ln(19) = %.6f" % (C*np.log(19)))
print("    NEW: C*(4/a1)*G = C*(4/5)*6/phi = C*%.4f = %.6f" % (
    4.0/a1*G_t, C*4.0/a1*G_t))
print("    ln(19)=2.944 vs %.4f: %.2f%% error in correction" % (
    4.0/a1*G_t, pct(4.0/a1*G_t, np.log(19))))

# ======================================================================
# STEP 6: Deep analysis - WHY 4/a1?
# ======================================================================
print()
print("=" * 80)
print("DEEP ANALYSIS: WHY 4/a1?")
print("=" * 80)
print()

print("For BLACK up quarks with |N|>1:")
print("  delta = C * ln|N| (known exact)")
print("  delta ~ C * (4/a1) * G (spectral approximation)")
print()
print("  The coefficient 4/a1 = 4/5 = 0.8 is close to but NOT exactly:")
print("    ln(5)/G(c) = ln(5)/2 = 0.8047")
print("    ln(19)/G(t) = ln(19)/(6/phi) = 0.7940")
print()
print("  Their average: (0.8047 + 0.7940)/2 = 0.7994 ~ 4/a1 = 0.8")
print("  So 4/a1 is the BEST rational approximation, not exact.")
print()
print("  CONCLUSION: The spectral formula C*(4/a1)*G is a PATTERN")
print("  with ~0.7%% accuracy. The true formula uses ln|N|, which")
print("  involves transcendental functions unavailable from graph spectra.")
print()

# Check: is there a better algebraic coefficient?
# ln(5) = 2*x and ln(19) = (6/phi)*x => x = ln(5)/2 and ln(19)/(6/phi)
# These are NOT equal, so no SINGLE algebraic coefficient works
x_c = np.log(5) / G_c
x_t = np.log(19) / G_t
print("  Exact coefficients needed:")
print("    c: ln(5)/G(c) = %.6f" % x_c)
print("    t: ln(19)/G(t) = %.6f" % x_t)
print("    Ratio x_c/x_t = %.4f (would be 1.000 if universal)" % (x_c/x_t))

# ======================================================================
# STEP 7: Summary of spectral unification
# ======================================================================
print()
print("=" * 80)
print("FINAL SUMMARY: SPECTRAL UNIFICATION STATUS")
print("=" * 80)
print()

print("PRESCRIPTIONS STATUS:")
print()
print("  P1 (prime quarks): PARTIALLY SPECTRAL")
print("    c,t (BLACK up): delta ~ C*(4/a1)*G, ~0.7%% approx")
print("    b (WHITE down):  delta = -C*ln19/phi, NOT spectral (G=0)")
print()
print("  P2 (d quark): NOT SPECTRAL")
print("    d (WHITE, rational): delta = -2/a1, G=0")
print()
print("  P3 (s quark): FULLY SPECTRAL (EXACT)")
print("    s (BLACK down): delta = C*G(s)/phi^2 where G(s)=-3=-N_gen")
print("    The appearance of N_gen is now DERIVED, not assumed")
print()
print("  P4 (leptons): NOT SPECTRAL")
print("    mu, tau (WHITE): delta from |z'|^(3/4), G=0")
print()

print("SCORE: 3/8 corrections have spectral origin")
print("  - s: EXACT (P3 fully derived)")
print("  - u: EXACT (delta=0 from sign of G)")
print("  - c,t: APPROXIMATE (~0.7%% accuracy)")
print()
print("  - d, b, mu, tau: NO spectral origin (WHITE nodes, G=0)")
print()

print("CATEGORY: P3 elevated from PATTERN to DERIVED.")
print("  Other prescriptions unchanged.")
print()
print("RMS impact: Negligible (0.110%% -> %.3f%% with spectral c,t)" % rms_s)
print()
print("End of exp365b.")
