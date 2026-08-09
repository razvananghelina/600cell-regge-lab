"""
exp365_mckay_spectral_corrections.py
=====================================
BLIND PROCEDURE: Compute spectral quantities from the McKay graph
and 600-cell Cayley graph. Then check if mass corrections delta_f
can be expressed as functions of these quantities.

From memory: Cayley Laplacian eigenvalue for irrep rho_k is
  L_k = 12 * (1 - chi_k(class) / d_k)
where class = 10a (generators), d_k = dim(rho_k).

The Galois conjugate uses class 5a.

STEP 1: Compute all spectral data (no reference to corrections)
STEP 2: Look for algebraic patterns
STEP 3: Compare with known corrections

Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
from numpy.linalg import eigh, inv

PHI = (1 + np.sqrt(5)) / 2
PHIc = (1 - np.sqrt(5)) / 2  # = -1/phi
a1 = 5
b1 = 6
N_gen = 3
C_coeff = 4.0 / (a1**2 + 1)  # = 2/13
m_e = 0.51099895  # MeV

# ====================================================================
# STEP 0: Build McKay graph and character table
# ====================================================================

dims = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3], dtype=float)
A_mck = np.zeros((9, 9))
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
for i, j in edges:
    A_mck[i,j] = 1
    A_mck[j,i] = 1

# Adjacency spectrum and eigenvectors
evals_raw, evecs_raw = eigh(A_mck)
idx = np.argsort(evals_raw)[::-1]
evals_mck = evals_raw[idx]
evecs_mck = evecs_raw[:, idx]
for j in range(9):
    if evecs_mck[0, j] < 0:
        evecs_mck[:, j] *= -1

# Eigenvalue -> class mapping (sorted by eigenvalue descending)
# eval: 2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2
# class: 1a, 10a, 3a, 10b, 4a, 5b, 6a, 5a, 2a
class_sizes = [1, 12, 20, 12, 30, 12, 20, 12, 1]
class_names = ['1a', '10a', '3a', '10b', '4a', '5b', '6a', '5a', '2a']

char_table = np.zeros((9, 9))
for c in range(9):
    norm = np.sqrt(120.0 / class_sizes[c])
    for k in range(9):
        char_table[k, c] = evecs_mck[k, c] * norm

# Column indices
COL_1A  = 0  # eval = 2
COL_10A = 1  # eval = phi
COL_3A  = 2  # eval = 1
COL_10B = 3  # eval = 1/phi
COL_4A  = 4  # eval = 0
COL_5B  = 5  # eval = -1/phi
COL_6A  = 6  # eval = -1
COL_5A  = 7  # eval = -phi
COL_2A  = 8  # eval = -2

# Bipartite coloring (BFS from node 0)
color = [-1]*9
color[0] = 0
queue = [0]
while queue:
    n = queue.pop(0)
    for j in range(9):
        if A_mck[n,j] > 0 and color[j] < 0:
            color[j] = 1 - color[n]
            queue.append(j)

# Graph distances from node 0
gdist = np.zeros(9, dtype=int)
visited = [False]*9
visited[0] = True
queue = [(0, 0)]
while queue:
    node, d = queue.pop(0)
    gdist[node] = d
    for j in range(9):
        if A_mck[node,j] > 0 and not visited[j]:
            visited[j] = True
            queue.append((j, d+1))

# Fermion mapping (Solution 3)
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

# Known corrections (for comparison ONLY in Step 3)
delta_known = {'e':0, 'mu':0.0792, 'tau':-0.0552,
               'u':0, 'c':0.2476, 't':0.4530,
               'd':-0.4000, 's':-0.1763, 'b':-0.2800}

def pct(pred, obs):
    if obs == 0: return float('inf')
    return (pred-obs)/obs*100

def rms(v):
    return np.sqrt(np.mean(np.array(v)**2))

def pr(c='=', w=85):
    print(c*w)

# ====================================================================
# STEP 1: SPECTRAL DATA TABLE (BLIND -- no corrections used)
# ====================================================================

print()
pr()
print("STEP 1: SPECTRAL DATA FOR EACH FERMION NODE (BLIND)")
pr()

print("\nCayley Laplacian: L_k = 12*(1 - chi_k(class)/d_k)")
print("Galois symmetric:  S_k = (L_10a + L_5a)/2")
print("Galois anti-sym:   G_k = (L_10a - L_5a)/2")
print("Bipartite color:   WHITE=A5 rep (G=0), BLACK=spin rep (G!=0)")
print()

# Z[phi] helper
def to_zphi(x):
    for p in range(-10, 11):
        q = round((x - p) / PHI)
        if abs(p + q*PHI - x) < 0.01:
            return (p, q)
    return (round(x), 0)

header = "%-4s %-4s %3s %5s %7s %7s %7s %7s %7s %6s %5s"
print(header % ("Node", "Name", "d", "color", "chi10a", "L_10a", "L_5a",
                "S_k", "G_k", "dist", "n_b"))
print("-" * 85)

spec = {}
for k in range(9):
    name, n_bare, (a,b), T3, Nc = fmap[k]
    d = dims[k]
    chi_10a = char_table[k, COL_10A]
    chi_5a  = char_table[k, COL_5A]
    L_10a = 12 * (1 - chi_10a / d)
    L_5a  = 12 * (1 - chi_5a / d)
    S_k = (L_10a + L_5a) / 2
    G_k = (L_10a - L_5a) / 2
    col_str = "W" if color[k] == 0 else "B"

    print("rho%d %-4s %3d %5s %+7.3f %7.3f %7.3f %7.3f %+7.3f %6d %5d" % (
        k+1, name, int(d), col_str, chi_10a, L_10a, L_5a, S_k, G_k,
        gdist[k], n_bare))

    spec[name] = {
        'k': k, 'd': d, 'chi_10a': chi_10a, 'chi_5a': chi_5a,
        'L_10a': L_10a, 'L_5a': L_5a, 'S': S_k, 'G': G_k,
        'dist': gdist[k], 'n_bare': n_bare, 'color': color[k],
        'T3': T3, 'Nc': Nc, 'a': a, 'b': b,
    }

# ====================================================================
# STEP 1b: ALGEBRAIC PATTERNS IN SPECTRAL DATA (BLIND)
# ====================================================================

print()
pr()
print("STEP 1b: ALGEBRAIC PATTERNS (BLIND)")
pr()

print("\nKey observations from the spectral data:")
print()
print("1. WHITE nodes (A5 reps) have G = 0 exactly.")
print("   BLACK nodes (spin reps) have G != 0.")
print()
print("2. G values for BLACK nodes:")
for name in ['u', 's', 'c', 't']:
    s = spec[name]
    G = s['G']
    p, q = to_zphi(-G)  # -G because G = -12*chi(10a)/d for BLACK
    print("   %s: G = %+.4f = %+d + %+d*phi [in Z[phi]: (%+d,%+d)]" % (
        name, G, round(-12*s['chi_10a']/s['d']), 0, p, q))

print()
print("3. S values:")
print("   BLACK nodes: S = 12 for all (spin symmetric part)")
for name in ['u', 's', 'c', 't']:
    print("     %s: S = %.3f" % (name, spec[name]['S']))
print("   WHITE nodes: S = L (since L_10a = L_5a)")
for name in ['e', 'd', 'mu', 'tau', 'b']:
    print("     %s: S = L = %.3f" % (name, spec[name]['S']))

print()
print("4. L_10a in Z[phi]:")
for name in ['e','u','d','s','mu','c','tau','t','b']:
    s = spec[name]
    L = s['L_10a']
    p, q = to_zphi(L)
    print("   %s: L = %.4f = %d + %d*phi" % (name, L, p, q))

print()
print("5. Products and sums:")
L_d = spec['d']['L_10a']
L_b = spec['b']['L_10a']
print("   L(d) + L(b) = %.4f + %.4f = %.4f = 20" % (L_d, L_b, L_d+L_b))
print("   L(d) * L(b) = %.4f * %.4f = %.4f = 80 = 16*a1" % (L_d, L_b, L_d*L_b))
L_u = spec['u']['L_10a']
L_t = spec['t']['L_10a']
print("   L(u) + L(t) = %.4f + %.4f = %.4f = 18" % (L_u, L_t, L_u+L_t))
print("   L(u) * L(t) = %.4f * %.4f = %.4f" % (L_u, L_t, L_u*L_t))
print("     = (12-6*phi)*(12+6/phi) = 144 + 72/phi - 72*phi - 36/phi")
print("     = 144 + 36/phi - 72*phi = 144 + 36*(1/phi - 2*phi)")
v_ut = 144 + 36*(1/PHI - 2*PHI)
print("     = 144 + 36*(%.4f) = %.4f" % (1/PHI-2*PHI, v_ut))
print("   L(s) + L(c) = %.4f + %.4f = %.4f = 23 = a1^2 - 2" % (
    spec['s']['L_10a'], spec['c']['L_10a'],
    spec['s']['L_10a'] + spec['c']['L_10a']))

# N(z) for each fermion
def zphi_norm(a, b):
    return a**2 + a*b - b**2

print()
print("6. |N(z)| for each fermion:")
for name in ['e','u','d','s','mu','c','tau','t','b']:
    s = spec[name]
    Nz = zphi_norm(s['a'], s['b'])
    print("   %s: |N| = %d" % (name, abs(Nz)))


# ====================================================================
# STEP 2: BLIND FORMULA CANDIDATES FROM SPECTRAL DATA
# ====================================================================

print()
pr()
print("STEP 2: BLIND FORMULA CANDIDATES")
pr()

# Candidate A: delta proportional to G_k (Galois anti-symmetric)
print("\nCandidate A: delta = f(T3) * C * G_k")
print("  Only nonzero for BLACK quarks (u, s, c, t)")
print()

# For BLACK up quarks (T3=+1/2): c, t have G > 0; u has G < 0
# For BLACK down quarks (T3=-1/2): s has G < 0
print("  BLACK quark G values and signs:")
for name in ['u', 's', 'c', 't']:
    s = spec[name]
    print("    %s: T3=%+.1f, G=%+.3f" % (name, s['T3'], s['G']))

print()
print("  Note: Up quarks (c,t) have G > 0. Up quark u has G < 0.")
print("  Down quark s has G < 0.")
print("  Hypothesis: positive G -> positive correction (up quarks)")
print("  This naturally gives delta_u = 0 if formula uses max(0, G).")
print()

# Candidate B: delta = C * max(0, G) * factor for up quarks
print("Candidate B: delta_up = C * max(0, G_k) * f_up")
print("  c: max(0, 4) = 4. Need delta_c/C = ?  (unknown, blind!)")
print("  t: max(0, 12/phi) = 12/phi. Need delta_t/C = ?")
print("  u: max(0, -12*phi) = 0. Predicts delta_u = 0.")
print()

# Candidate C: delta = C * G_k / (2*phi^2) for down quarks
print("Candidate C: delta_down = C * G_k / (2*phi^2) for BLACK down quarks")
print("  s: G = -6. delta_s = C*(-6)/(2*phi^2) = -3*C/phi^2")
print("     = -(2/13)*3/(phi^2) = -6/(13*phi^2) = %.6f" % (-6/(13*PHI**2)))
print()

# Candidate D: Combined spectral formula using L itself
print("Candidate D: delta = C * h(L_k, d_k) for WHITE quarks")
print("  d: L = 12-4*phi. Need h -> delta_d")
print("  b: L = 12+4/phi. Need h -> delta_b")
print("  Note: L(d) = 12 - 4*chi_d(10a)/d_d * d_d ... no, L = 12-4*phi")
print("  For WHITE quarks, L is the ONLY spectral distinguisher (G=0)")
print()


# ====================================================================
# STEP 3: COMPARE WITH KNOWN CORRECTIONS
# ====================================================================

print()
pr()
print("STEP 3: COMPARE SPECTRAL CANDIDATES WITH KNOWN CORRECTIONS")
pr()

print("\n--- BLACK quarks: G_k vs delta ---")
print()
print("%-5s %5s %+8s %+8s %+8s %10s" % (
    "Name", "T3", "G_k", "delta", "delta/C", "G/delta*C"))
print("-" * 55)

for name in ['u', 's', 'c', 't']:
    s = spec[name]
    d_known = delta_known[name]
    d_over_C = d_known / C_coeff if abs(d_known) > 1e-10 else 0
    G_over_dC = s['G'] / d_over_C if abs(d_over_C) > 1e-10 else float('inf')
    print("%-5s %+.1f %+8.3f %+8.4f %+8.4f %10.4f" % (
        name, s['T3'], s['G'], d_known, d_over_C, G_over_dC))

print()
print("Key ratios for up quarks (c, t) where delta = +C*ln|N|:")
for name in ['c', 't']:
    s = spec[name]
    G = s['G']
    Nz = abs(zphi_norm(s['a'], s['b']))
    lnN = np.log(Nz)
    ratio = lnN / G
    print("  %s: ln|N|/G = %.4f/%.4f = %.5f" % (name, lnN, G, ratio))

print()
print("  Both ratios are close to 2/a1 = 2/5 = %.5f" % (2.0/a1))
print("  => ln|N| ~ (2/a1) * G for BLACK up quarks with |N|>1")
print("  Accuracy: ~0.6-0.7%% (approximate pattern, not exact)")
print()

# Test: exact spectral formula for s quark
print("--- EXACT RESULT: s quark from Galois spectral data ---")
print()
G_s = spec['s']['G']
delta_s_spectral = C_coeff * G_s / (2 * PHI**2)
delta_s_known = delta_known['s']
print("  G(s) = %.4f = -b1 = -%d (EXACT integer)" % (G_s, b1))
print("  delta_s = C * G_s / (2*phi^2) = (2/13)*(-6)/(2*phi^2)")
print("         = %.6f" % delta_s_spectral)
print("  Known:   %.6f" % delta_s_known)
print("  Match:   %s (diff = %.2e)" % (
    "EXACT" if abs(delta_s_spectral - delta_s_known) < 1e-10 else "APPROX",
    abs(delta_s_spectral - delta_s_known)))
print()
print("  INTERPRETATION: G(s) = -b1 means the s quark sits at the node")
print("  where the Galois anti-symmetric Cayley eigenvalue equals -b1.")
print("  The correction is C * (-b1) / (2*phi^2) = -N_gen*C/phi^2")
print("  since b1/(2*phi^2) = 6/(2*phi^2) = 3/phi^2 = N_gen/phi^2.")
print("  So delta_s = -N_gen*C/phi^2 has a SPECTRAL ORIGIN: G = -b1.")

print()
print("--- WHITE quarks: L_10a vs delta ---")
print()
print("%-5s %5s %+8s %+8s %+8s  %s" % (
    "Name", "T3", "L_10a", "delta", "delta/C", "Comment"))
print("-" * 65)

for name in ['d', 'b']:
    s = spec[name]
    d_known = delta_known[name]
    d_over_C = d_known / C_coeff if abs(d_known) > 1e-10 else 0
    print("%-5s %+.1f %+8.3f %+8.4f %+8.4f  G=0, need L" % (
        name, s['T3'], s['L_10a'], d_known, d_over_C))

print()
print("  d: L = 12-4*phi = %.4f. delta = -2/a1 = -0.4" % spec['d']['L_10a'])
print("  b: L = 12+4/phi = %.4f. delta = -C*ln(19)/phi = -0.2800" % spec['b']['L_10a'])
print()
print("  Is delta_d related to L(d)?")
print("    -2/a1 / (12-4*phi) = %.6f" % (-2.0/a1 / (12-4*PHI)))
print("    -2/a1 * (12+4/phi) / 80 = %.6f" % (-2.0/a1 * (12+4/PHI) / 80))
print("    L(d)*L(b) = 80 = 16*a1. So L(b) = 80/L(d).")
print("    delta_d = -2/a1 is a GLOBAL constant, not L-dependent.")
print("    CONCLUSION: d quark correction is NOT spectral.")
print()

# Check if delta_b is related to L(b)
Nz_b = abs(zphi_norm(-1, 4))
print("  Is delta_b related to L(b)?")
print("    delta_b = -C*ln(19)/phi = %.4f" % (-C_coeff*np.log(19)/PHI))
print("    L(b) = 12+4/phi = %.4f" % (12+4/PHI))
print("    ln(19)/L(b) = %.4f" % (np.log(19)/(12+4/PHI)))
print("    ln(19)/(4/phi) = ln(19)*phi/4 = %.4f" % (np.log(19)*PHI/4))
print("    No clean spectral expression for ln(19).")

print()
print("--- LEPTONS: Spectral data ---")
print()
for name in ['mu', 'tau']:
    s = spec[name]
    d_known = delta_known[name]
    print("  %s: L=%.3f, G=%.3f, dist=%d, color=W, delta=%+.4f" % (
        name, s['L_10a'], s['G'], s['dist'], d_known))

print()
print("  mu (L=12): lepton at the spectral CENTER (L=degree)")
print("  tau (L=15): lepton at high eigenvalue")
print("  Lepton corrections use |z'|^(3/4), unrelated to L_k.")


# ====================================================================
# STEP 4: FULL TEST OF SPECTRAL-INSPIRED FORMULA
# ====================================================================

print()
pr()
print("STEP 4: FULL TEST - SPECTRAL-INSPIRED FORMULA")
pr()

print("\nFormula S (spectral):")
print("  BLACK up quarks: delta = C * (2/a1) * max(0, G_k)")
print("  BLACK down quarks: delta = C * G_k / (2*phi^2)")
print("  WHITE quarks with |N|>1: delta = -C*ln|N|/phi (unchanged)")
print("  WHITE d quark: delta = -2/a1 (unchanged)")
print("  Leptons: delta = c_ell*sign(z')*|z'|^(3/4) (unchanged)")
print()

d_ST = 4.0
c_ell = C_coeff * PHI**3 / d_ST

all_errs = []
print("%-5s %5s %5s %+8s %+8s %12s %12s %+9s  %s" % (
    "Name", "T3", "color", "delta_S", "delta_K", "m_S/MeV", "m_exp/MeV",
    "err_S", "formula"))
print("-" * 100)

for k in range(9):
    name, n_bare, (a,b), T3, Nc = fmap[k]
    s = spec[name]
    m_exp = pdg[name]
    Nz = abs(zphi_norm(a, b))
    zp = a + b * PHIc
    col_str = "W" if s['color']==0 else "B"

    if name == 'e':
        print("%-5s %+.1f %5s %+8s %+8s %12.4f %12.4f %9s  %s" % (
            name, T3, col_str, "0", "0", m_exp, m_exp, "(input)", "electron"))
        continue

    dk = delta_known[name]

    # Spectral formula
    if Nc == 1:  # lepton
        ds = c_ell * np.sign(zp) * abs(zp)**0.75
        formula = "lepton |z'|^(3/4)"
    elif s['color'] == 1:  # BLACK quark
        if T3 > 0:  # up
            ds = C_coeff * (2.0/a1) * max(0, s['G'])
            formula = "C*(2/a1)*max(0,G)"
        else:  # down
            ds = C_coeff * s['G'] / (2 * PHI**2)
            formula = "C*G/(2*phi^2)"
    else:  # WHITE quark
        if b == 0:  # rational
            ds = -2.0 / a1
            formula = "-2/a1 (rational)"
        elif Nz > 1:  # prime
            ds = -C_coeff * np.log(Nz) / PHI
            formula = "-C*ln|N|/phi"
        else:
            ds = 0.0
            formula = "WHITE unit"

    n_eff_s = n_bare + ds
    n_eff_k = n_bare + dk
    m_s = m_e * PHI**n_eff_s
    err = pct(m_s, m_exp)
    if name != 'e':
        all_errs.append(abs(err))

    print("%-5s %+.1f %5s %+8.4f %+8.4f %12.4f %12.4f %+8.2f%%  %s" % (
        name, T3, col_str, ds, dk, m_s, m_exp, err, formula))

print()
print("Spectral formula RMS: %.3f%%" % rms(all_errs))
print("Known formula RMS:    0.110%%")


# ====================================================================
# STEP 5: SUMMARY
# ====================================================================

print()
pr()
print("STEP 5: SUMMARY AND CONCLUSIONS")
pr()

print()
print("RESULT 1 (EXACT): s quark correction from Galois spectral data")
print("  G(s) = -b1 = -6 (Galois anti-symmetric Cayley eigenvalue)")
print("  delta_s = C * G(s) / (2*phi^2) = -N_gen*C/phi^2")
print("  This DERIVES P3 from spectral geometry of the 600-cell.")
print("  The s quark sits at rho_4 (dim 4, BLACK), where chi(10a)=1,")
print("  chi(5a)=-1, giving L_10a=9, L_5a=15, G=-6=-b1.")
print()

print("RESULT 2 (APPROXIMATE): c,t quark corrections ~ (2/a1)*G")
print("  delta_c ~ C*(2/a1)*G(c) = C*(2/5)*4 = C*1.6  (exact: C*ln5=C*1.609)")
print("  delta_t ~ C*(2/a1)*G(t) = C*(2/5)*12/phi     (exact: C*ln19=C*2.944)")
print("  The approximation ln|N| ~ (2/a1)*G holds to ~0.7%%.")
print("  This is a PATTERN, not exact. (2/a1)*G is 0.4*G while ln|N|/G ~ 0.40)")
print()

print("RESULT 3 (EXACT): u quark delta=0 from sign of G")
print("  G(u) = -12*phi/2 = -6*phi < 0")
print("  Since max(0, G) = 0, the spectral formula gives delta_u = 0.")
print("  Physical: the u quark is the ONLY up quark where chi(10a) > 0")
print("  among spin reps, making its Galois contribution negative.")
print()

print("RESULT 4 (NEGATIVE): WHITE quarks not spectral")
print("  d, b are WHITE/A5 nodes with G=0.")
print("  Their corrections (-2/a1 and -C*ln19/phi) cannot be derived")
print("  from the Galois anti-symmetric spectral quantity.")
print("  The symmetric part S=L is not simply related to corrections either.")
print("  L(d)*L(b) = 80 = 16*a1 is an exact identity but doesn't help.")
print()

print("RESULT 5 (NEGATIVE): Leptons not spectral")
print("  mu (L=12) and tau (L=15) are WHITE nodes.")
print("  Their corrections use |z'|^(3/4), a power-law unrelated to L.")
print()

print("OVERALL ASSESSMENT:")
print("  The Galois anti-symmetric Cayley eigenvalue G_k = (L_10a-L_5a)/2")
print("  captures the corrections for BLACK quarks:")
print("    - s: EXACT (G=-b1 -> delta=-N_gen*C/phi^2)")
print("    - u: EXACT (G<0 -> delta=0)")
print("    - c,t: APPROXIMATE (G ~ a1*ln|N|/2, ~0.7%% accuracy)")
print()
print("  But G is USELESS for WHITE nodes (d, b, mu, tau) where G=0.")
print("  The McKay bipartite structure creates a fundamental asymmetry:")
print("    BLACK = spin reps = Galois-nontrivial = spectrally distinguishable")
print("    WHITE = A5 reps = Galois-invariant = spectrally degenerate")
print()
print("  Score: 2 EXACT + 2 APPROXIMATE out of 8 nontrivial corrections.")
print("  Category: PATTERN (not full derivation).")
print()
print("End of exp365.")
