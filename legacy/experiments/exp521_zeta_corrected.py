"""
EXP-521: Spectral Zeta with CORRECT Eigenvalues
=================================================
All eigenvalues from commons/constants.py (verified numerically).

Key corrections from exp520:
  Pair B: (10-2sqrt5, 10+2sqrt5) NOT (12-2sqrt5, 12+2sqrt5)
  L2*L6 = 80 = d_ST^2 * a1  (NOT 124)

THREE zeta functions to explore:
  1. zeta_broken(s) = 4*L1^{-s} + 9*L2^{-s} + 9*L6^{-s} + 4*L8^{-s}
  2. zeta_asym(s) = 4*(L1^{-s} - L8^{-s}) + 9*(L2^{-s} - L6^{-s})
  3. zeta_norm(s) = 4*36^{-s} + 9*80^{-s}  (Galois norms)
"""

import numpy as np
import sys
sys.path.insert(0, '.')
from commons import (PHI, PHI_CONJ, SQRT5, a1, b1, N, N_gen, degree, d_ST,
                      alpha_em, sin2_tW, inv_alpha, Omega_DM,
                      L0, L1, L2, L3, L4, L5, L6, L7, L8,
                      GALOIS_NORM_A, GALOIS_NORM_B,
                      D2_PHYS, D2_DARK, d_dark,
                      SPECTRUM_600CELL, N_GALOIS_BROKEN)

print("=" * 70)
print("EXP-521: SPECTRAL ZETA WITH CORRECT EIGENVALUES")
print("=" * 70)

# Quick verification
print(f"\n  Eigenvalues: L1={L1:.6f}, L2={L2:.6f}, L6={L6:.6f}, L8={L8:.6f}")
print(f"  Norms: L1*L8={GALOIS_NORM_A:.0f}=b1^2, L2*L6={GALOIS_NORM_B:.0f}=d_ST^2*a1")
print(f"  L8/L1={L8/L1:.6f} = phi^4={PHI**4:.6f}")

# ============================================================
# DEFINE THE THREE ZETA FUNCTIONS
# ============================================================

def zeta_total(s):
    return sum(m * ev**(-s) for ev, m, _, _ in SPECTRUM_600CELL if ev > 0.01)

def zeta_broken(s):
    return 4*L1**(-s) + 9*L2**(-s) + 9*L6**(-s) + 4*L8**(-s)

def zeta_fixed(s):
    return 16*L3**(-s) + 25*L4**(-s) + 36*L5**(-s) + 16*L7**(-s)

def zeta_asym(s):
    """Physical minus dark within each Galois pair."""
    return 4*(L1**(-s) - L8**(-s)) + 9*(L2**(-s) - L6**(-s))

def zeta_norm(s):
    """Galois norm zeta: uses products L1*L8=36, L2*L6=80."""
    return 4 * 36.0**(-s) + 9 * 80.0**(-s)

# ============================================================
# PART 1: VALUES AT CLEAN s
# ============================================================
print(f"\n{'='*70}")
print("PART 1: ZETA AT CLEAN FRAMEWORK s VALUES")
print(f"{'='*70}")

clean_s = [
    ("1/a1",     1/a1),
    ("1/d_ST",   1/d_ST),
    ("1/N_gen",  1/N_gen),
    ("1/b1",     1/b1),
    ("1/phi",    1/PHI),
    ("1/2",      0.5),
    ("1/phi^2",  1/PHI**2),
    ("phi-1",    PHI-1),   # = 1/phi
    ("ln(phi)",  np.log(PHI)),
    ("1",        1.0),
    ("phi",      PHI),
    ("phi^2",    PHI**2),
    ("2",        2.0),
    ("phi^3",    PHI**3),
    ("3",        3.0),
    ("pi",       np.pi),
    ("d_ST",     float(d_ST)),
    ("a1",       float(a1)),
    ("b1",       float(b1)),
]

targets = {
    '1/alpha': inv_alpha,
    'alpha': alpha_em,
    '4a1*phi^4': 4*a1*PHI**4,
    'sin2tW': sin2_tW,
    '7-phi': 7-PHI,
    'b1': float(b1),
    'a1': float(a1),
    'degree': float(degree),
    'phi': PHI,
    'phi^2': PHI**2,
    'phi^3': PHI**3,
    'phi^4': PHI**4,
    '2sqrt5': 2*SQRT5,
    'sqrt5': SQRT5,
    '2': 2.0,
    'pi': np.pi,
    'N_brk': 26.0,
    'N': 120.0,
}

print(f"\n  {'s_name':>10} {'s':>8} | {'z_total':>10} {'z_brk':>10} {'z_asym':>10} {'z_norm':>10} {'z_fix':>10}")
print(f"  {'':>10} {'':>8} | {'':>10} {'':>10} {'':>10} {'':>10} {'':>10}")

for sname, sval in clean_s:
    zt = zeta_total(sval)
    zb = zeta_broken(sval)
    za = zeta_asym(sval)
    zn = zeta_norm(sval)
    zf = zeta_fixed(sval)
    print(f"  {sname:>10} {sval:>8.4f} | {zt:>10.4f} {zb:>10.4f} {za:>10.4f} {zn:>10.6f} {zf:>10.4f}")

# ============================================================
# PART 2: SCAN FOR EXACT MATCHES
# ============================================================
print(f"\n{'='*70}")
print("PART 2: HIGH-PRECISION SCAN FOR FRAMEWORK CONSTANTS")
print(f"{'='*70}")

# For each zeta variant, find s* where it hits each target
best_matches = []
for s_test in np.arange(0.001, 10.001, 0.0005):
    funcs = {
        'z_total': zeta_total(s_test),
        'z_brk': zeta_broken(s_test),
        'z_asym': zeta_asym(s_test),
        'z_norm': zeta_norm(s_test),
        '1/z_asym': 1/zeta_asym(s_test) if abs(zeta_asym(s_test)) > 1e-12 else 0,
    }
    for fname, fval in funcs.items():
        if abs(fval) < 1e-12:
            continue
        for tname, tval in targets.items():
            if tval < 1e-6:
                continue
            err = abs(fval/tval - 1)
            if err < 0.0005:  # 0.05% precision
                best_matches.append((err, s_test, fname, fval, tname, tval))

best_matches.sort()
seen = set()
print(f"\n  {'s*':>10} {'function':>10} {'value':>12} {'target':>12} {'name':>10} {'err':>10}")
print(f"  {'-'*10} {'-'*10} {'-'*12} {'-'*12} {'-'*10} {'-'*10}")
for err, s, fn, fv, tn, tv in best_matches:
    key = (fn, tn)
    if key in seen:
        continue
    seen.add(key)
    # Check if s is a clean expression
    s_id = ""
    for sname, sval in clean_s:
        if abs(s - sval) < 0.002:
            s_id = f" = {sname}"
            break
    print(f"  {s:>10.4f} {fn:>10} {fv:>12.6f} {tv:>12.6f} {tn:>10} {err:>10.6f}{s_id}")

# ============================================================
# PART 3: GALOIS NORM ZETA AT FINE GRID
# ============================================================
print(f"\n{'='*70}")
print("PART 3: GALOIS NORM ZETA - ANALYTIC STRUCTURE")
print(f"{'='*70}")

# zeta_norm(s) = 4*36^{-s} + 9*80^{-s}
# = 4*b1^{-2s} + 9*(d_ST^2*a1)^{-s}
# = 4*6^{-2s} + 9*80^{-s}

# At s=0: zeta_norm(0) = 4+9 = 13 = half of N_broken
print(f"  zeta_norm(0) = 4+9 = {zeta_norm(0):.0f} = N_broken/2 = {N_GALOIS_BROKEN//2}")

# At s -> inf: zeta_norm -> 4*36^{-s} (dominated by smaller norm)
# At s -> -inf: zeta_norm -> 9*80^{-s} (dominated by larger norm)

# Find s* where zeta_norm = framework constants
print(f"\n  Solving zeta_norm(s*) = target for each framework constant:")
print(f"  {'target':>12} {'value':>10} {'s*':>12} {'s* clean?':>30}")

for tname in ['degree', 'phi^4', 'b1', '7-phi', 'a1', 'phi^3', 'pi',
              'phi^2', 'sqrt5', '2', 'phi', 'sin2tW', 'alpha']:
    tval = targets[tname]
    # Binary search for s*
    s_lo, s_hi = -0.5, 5.0
    if tval > 13:  # larger than zeta_norm(0)
        s_lo, s_hi = -2.0, 0.0
    for _ in range(100):
        s_mid = (s_lo + s_hi) / 2
        if zeta_norm(s_mid) > tval:
            s_lo = s_mid
        else:
            s_hi = s_mid
    s_star = (s_lo + s_hi) / 2

    # Check if s* is a clean expression
    s_clean = ""
    for sname, sval in clean_s:
        if abs(s_star - sval) < 0.005:
            s_clean = f"~ {sname} = {sval:.6f}"
            break
    # Also check simple fractions
    for num in range(1, 20):
        for den in range(1, 20):
            if abs(s_star - num/den) < 0.003:
                s_clean = f"~ {num}/{den} = {num/den:.6f}"
                break

    print(f"  {tname:>12} {tval:>10.4f} {s_star:>12.6f} {s_clean:>30}")

# ============================================================
# PART 4: ASYMMETRY ZETA DEEP ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("PART 4: ASYMMETRY ZETA - THE KEY FUNCTION")
print(f"{'='*70}")

# zeta_asym(s) = 4*(L1^{-s} - L8^{-s}) + 9*(L2^{-s} - L6^{-s})
# Since L8/L1 = phi^4 and L2*L6 = 80:
# L1^{-s} - L8^{-s} = L1^{-s}*(1 - (L8/L1)^{-s}) = L1^{-s}*(1 - phi^{-4s})

# So zeta_asym(s) = 4*L1^{-s}*(1-phi^{-4s}) + 9*L2^{-s}*(1-(L6/L2)^{-s})

# At s=1: asym = 4*(1/L1-1/L8) + 9*(1/L2-1/L6)
print(f"  zeta_asym(s) = 4*L1^{{-s}}*(1-phi^{{-4s}}) + 9*L2^{{-s}}*(1-(L6/L2)^{{-s}})")
print(f"  L6/L2 = {L6/L2:.6f}")

# The factor (1-phi^{-4s}) is interesting:
# At s=1/4: 1-phi^{-1} = 1-1/phi = 1-(phi-1) = 2-phi = 1/phi^2 = d_dark^2
# At s=1/2: 1-phi^{-2} = 1-1/phi^2 = 1-(2-phi) = phi-1 = 1/phi = d_dark

print(f"\n  The phi-factor 1-phi^{{-4s}}:")
for sname, sval in [("1/4", 0.25), ("1/2", 0.5), ("1", 1.0), ("phi/4", PHI/4),
                     ("1/phi", 1/PHI), ("phi", PHI)]:
    pf = 1 - PHI**(-4*sval)
    print(f"    s={sname:>6}: 1-phi^{{-4s}} = 1-phi^{{{-4*sval:>6.3f}}} = {pf:.6f}")

# At s=1/4: 1-phi^{-1} = 2-phi = 0.382 = d_dark^2
# At s=1/2: 1-phi^{-2} = phi-1 = 1/phi = d_dark

print(f"\n  KEY: at s=1/4, the phi-factor = 1-1/phi = 2-phi = d_dark^2 = {(2-PHI):.6f}")
print(f"       at s=1/2, the phi-factor = 1-1/phi^2 = phi-1 = d_dark = {(PHI-1):.6f}")

# Now compute zeta_asym at these special points
za_quarter = zeta_asym(0.25)
za_half = zeta_asym(0.5)

print(f"\n  zeta_asym(1/4) = {za_quarter:.6f}")
print(f"  zeta_asym(1/2) = {za_half:.6f}")

# Check: at s=1/2, asym = 4*L1^{-1/2}*(1/phi) + 9*L2^{-1/2}*(1-(L6/L2)^{-1/2})
L6_over_L2 = L6/L2
factor_B_half = 1 - L6_over_L2**(-0.5)
print(f"\n  At s=1/2:")
print(f"    Pair A contrib: 4*L1^{{-1/2}}*d_dark = 4/{np.sqrt(L1):.4f} * {d_dark:.4f} = {4/np.sqrt(L1)*d_dark:.6f}")
print(f"    Pair B factor: 1-(L6/L2)^{{-1/2}} = {factor_B_half:.6f}")
print(f"    Pair B contrib: 9*L2^{{-1/2}}*{factor_B_half:.4f} = {9/np.sqrt(L2)*factor_B_half:.6f}")

# ============================================================
# PART 5: PRODUCT DETERMINANT APPROACH
# ============================================================
print(f"\n{'='*70}")
print("PART 5: DETERMINANT AND LOG-DETERMINANT")
print(f"{'='*70}")

# Broken det = L1^4 * L2^9 * L6^9 * L8^4
# = (L1*L8)^4 * (L2*L6)^9
# = 36^4 * 80^9
# = b1^8 * (d_ST^2*a1)^9

log_det_broken = 4*np.log(L1) + 9*np.log(L2) + 9*np.log(L6) + 4*np.log(L8)
log_det_norms = 4*np.log(36) + 9*np.log(80)

print(f"  det_broken = L1^4 * L2^9 * L6^9 * L8^4")
print(f"             = 36^4 * 80^9")
print(f"             = b1^8 * (d_ST^2*a1)^9")
print(f"  log(det_broken) = {log_det_broken:.6f}")
print(f"  via norms:        {log_det_norms:.6f}")
print(f"  Match: {abs(log_det_broken - log_det_norms) < 1e-10}")

# Galois-conjugate det: sigma swaps L1<->L8, L2<->L6
# So det_dark = L8^4 * L6^9 * L2^9 * L1^4 = SAME as det_broken!
# The determinant is Galois-INVARIANT.
print(f"\n  det_dark = L8^4*L6^9*L2^9*L1^4 = det_broken (Galois-invariant)")

# But we can split WITHIN each pair:
# det_phys_half = L1^4 * L2^9  (physical members only)
# det_dark_half = L8^4 * L6^9  (dark members only)
log_det_phys = 4*np.log(L1) + 9*np.log(L2)
log_det_dark = 4*np.log(L8) + 9*np.log(L6)

ratio = np.exp(log_det_dark - log_det_phys)
print(f"\n  det_phys_half = L1^4 * L2^9,  log = {log_det_phys:.6f}")
print(f"  det_dark_half = L8^4 * L6^9,  log = {log_det_dark:.6f}")
print(f"  det_dark/det_phys = {ratio:.6f}")
print(f"  = (L8/L1)^4 * (L6/L2)^9 = phi^16 * (L6/L2)^9")
print(f"  = {PHI**16:.4f} * {(L6/L2)**9:.4f} = {PHI**16 * (L6/L2)**9:.4f}")
print(f"  log(ratio) = {log_det_dark-log_det_phys:.6f}")
print(f"  = 16*ln(phi) + 9*ln(L6/L2)")
print(f"  = {16*np.log(PHI):.6f} + {9*np.log(L6/L2):.6f}")

# Check if log(ratio) relates to anything
lr = log_det_dark - log_det_phys
print(f"\n  log(det_dark/det_phys) = {lr:.6f}")
print(f"    / (4*pi) = {lr/(4*np.pi):.6f}")
print(f"    / phi = {lr/PHI:.6f}")
print(f"    / degree = {lr/degree:.6f}")
print(f"    / a1 = {lr/a1:.6f}")
print(f"    / (2*pi) = {lr/(2*np.pi):.6f}")
print(f"    vs ln(1/alpha) = {np.log(inv_alpha):.6f}")
print(f"    vs 4*pi/alpha = {4*np.pi/alpha_em:.6f}")

# ============================================================
# PART 6: THE s* VALUES - ARE THEY CLEAN?
# ============================================================
print(f"\n{'='*70}")
print("PART 6: SEARCHING FOR CLEAN s* EXPRESSIONS")
print(f"{'='*70}")

# For the most interesting matches, check if s* is a clean expression
# of framework constants

# Specifically for zeta_asym(s*) = phi^2:
# We found s* ≈ 0.398 and s* ≈ 0.665. Are these clean?

# s ≈ 0.398: check if it's 2/a1 = 0.4, or 1/phi^2 = 0.382, or ln(phi) = 0.481
# or 2/(a1+something)
print(f"  Looking for s* where zeta_asym(s*) = phi^2 = {PHI**2:.6f}:")

# Fine search around the two solutions
for region_lo, region_hi in [(0.35, 0.45), (0.62, 0.70)]:
    s_lo, s_hi = region_lo, region_hi
    for _ in range(200):
        s_mid = (s_lo + s_hi) / 2
        if zeta_asym(s_mid) > PHI**2:
            s_lo = s_mid
        else:
            s_hi = s_mid
    s_star = (s_lo + s_hi) / 2
    za_check = zeta_asym(s_star)
    print(f"    s* = {s_star:.10f},  zeta_asym = {za_check:.10f}")

    # Check against framework expressions
    candidates = [
        ("2/a1", 2/a1), ("1/phi^2", 1/PHI**2), ("ln(phi)", np.log(PHI)),
        ("1/SQRT5", 1/SQRT5), ("2/b1", 2/b1), ("1/(phi*SQRT5)", 1/(PHI*SQRT5)),
        ("phi/d_ST", PHI/d_ST), ("1/pi", 1/np.pi), ("2/(a1+1)", 2/(a1+1)),
        ("phi/(2*a1)", PHI/(2*a1)), ("1/N_gen", 1/N_gen), ("(phi-1)/phi", (PHI-1)/PHI),
        ("phi/a1", PHI/a1), ("2/degree", 2/degree), ("(sqrt5-1)/d_ST", (SQRT5-1)/d_ST),
        ("d_dark/phi", d_dark/PHI), ("2/(phi*a1)", 2/(PHI*a1)),
        ("ln(2)/SQRT5", np.log(2)/SQRT5), ("phi^2/d_ST", PHI**2/d_ST),
        ("(a1-d_ST)/(a1*phi)", (a1-d_ST)/(a1*PHI)),
        ("2*ln(phi)", 2*np.log(PHI)), ("b1/degree", b1/degree),
        ("1/(phi+1)", 1/(PHI+1)), ("phi/(phi+a1)", PHI/(PHI+a1)),
    ]
    for cname, cval in sorted(candidates, key=lambda x: abs(x[1]-s_star)):
        diff = abs(cval - s_star)
        if diff < 0.01:
            print(f"      {cname:>20} = {cval:.10f}  diff = {diff:.2e}")

# Now check for zeta_asym(s*) = sin^2(tW) = 6/26
print(f"\n  Looking for s* where zeta_asym(s*) = sin^2(tW) = {sin2_tW:.6f}:")
s_lo, s_hi = 3.0, 4.0
for _ in range(200):
    s_mid = (s_lo + s_hi) / 2
    if zeta_asym(s_mid) > sin2_tW:
        s_lo = s_mid
    else:
        s_hi = s_mid
s_star = (s_lo + s_hi) / 2
za_check = zeta_asym(s_star)
print(f"    s* = {s_star:.10f},  zeta_asym = {za_check:.10f}")
for cname, cval in [("phi^2+phi", PHI**2+PHI), ("a1/phi", a1/PHI), ("phi^3", PHI**3),
                     ("2*phi", 2*PHI), ("b1/phi", b1/PHI), ("N_gen+1/phi", N_gen+1/PHI),
                     ("pi+1/phi", np.pi+d_dark), ("a1-phi", a1-PHI), ("pi+ln(phi)", np.pi+np.log(PHI)),
                     ("degree/phi^2", degree/PHI**2), ("7/2", 3.5), ("d_ST-1/phi", d_ST-1/PHI),
                     ("ln(b1^2)", np.log(36)), ("ln(80)/2", np.log(80)/2)]:
    diff = abs(cval - s_star)
    if diff < 0.05:
        print(f"      {cname:>20} = {cval:.6f}  diff = {diff:.2e}")

# And for zeta_norm(s*) = sin^2(tW)
print(f"\n  Looking for s* where zeta_norm(s*) = sin^2(tW):")
s_lo, s_hi = 0.5, 1.5
for _ in range(200):
    s_mid = (s_lo + s_hi) / 2
    if zeta_norm(s_mid) > sin2_tW:
        s_lo = s_mid
    else:
        s_hi = s_mid
s_star = (s_lo + s_hi) / 2
zn_check = zeta_norm(s_star)
print(f"    s* = {s_star:.10f},  zeta_norm = {zn_check:.10f}")
for cname, cval in [("1", 1.0), ("phi/phi", 1.0), ("ln(phi^2)", 2*np.log(PHI)),
                     ("a1/a1", 1.0), ("1/phi+1/phi^2", 1/PHI+1/PHI**2),
                     ("phi^2-phi", PHI**2-PHI), ("d_dark+d_dark^2", d_dark+d_dark**2)]:
    diff = abs(cval - s_star)
    if diff < 0.05:
        print(f"      {cname:>25} = {cval:.6f}  diff = {diff:.2e}")

# Check zeta_norm at s=1 (the simplest value)
print(f"\n  zeta_norm(1) = 4/36 + 9/80 = {zeta_norm(1):.10f}")
print(f"    = 1/9 + 9/80 = {1/9+9/80:.10f}")
print(f"    = (80+81)/(9*80) = 161/720 = {161/720:.10f}")
print(f"    = 161/720")
print(f"    161 = ? (prime)")
print(f"    720 = N*b1 = {N}*{b1} = {N*b1}")

# ============================================================
# PART 7: zeta_asym AT s=1 - ANALYTIC
# ============================================================
print(f"\n{'='*70}")
print("PART 7: ANALYTIC EXPRESSIONS AT s=1")
print(f"{'='*70}")

# zeta_asym(1) = 4*(1/L1 - 1/L8) + 9*(1/L2 - 1/L6)
# = 4*(L8-L1)/(L1*L8) + 9*(L6-L2)/(L2*L6)
# = 4*(L8-L1)/36 + 9*(L6-L2)/80

diff_A = L8 - L1  # = 6+6phi - 12+6phi = 12phi-6 = 6(2phi-1) = 6sqrt5
diff_B = L6 - L2  # = 10+2s5 - 10+2s5 = 4sqrt5

za1 = 4*diff_A/36 + 9*diff_B/80
print(f"  zeta_asym(1) = 4*(L8-L1)/(L1*L8) + 9*(L6-L2)/(L2*L6)")
print(f"    L8-L1 = {diff_A:.6f} = 6*sqrt5 = {6*SQRT5:.6f}: {abs(diff_A-6*SQRT5)<1e-10}")
print(f"    L6-L2 = {diff_B:.6f} = 4*sqrt5 = {4*SQRT5:.6f}: {abs(diff_B-4*SQRT5)<1e-10}")
print(f"  = 4*6*sqrt5/36 + 9*4*sqrt5/80")
print(f"  = (24*sqrt5)/36 + (36*sqrt5)/80")
print(f"  = (2*sqrt5)/3 + (9*sqrt5)/20")
print(f"  = sqrt5*(2/3 + 9/20)")
print(f"  = sqrt5*(40+27)/60")
print(f"  = 67*sqrt5/60")

exact_za1 = 67*SQRT5/60
print(f"  = 67*sqrt5/60 = {exact_za1:.10f}")
print(f"  Numerical:      {zeta_asym(1):.10f}")
print(f"  Match: {abs(exact_za1 - zeta_asym(1)) < 1e-10}")

# Now what is 67/60 ?
# 67 is prime. 60 = 5! / 2 = a1! / 2
print(f"\n  67 is prime")
print(f"  60 = {60} = a1! / 2 = N / 2")
print(f"  So: zeta_asym(1) = 67*sqrt(a1) / (N/2) = 2*67*sqrt(a1)/N")

# zeta_broken(1) = 4/L1 + 9/L2 + 9/L6 + 4/L8
zb1 = 4/L1 + 9/L2 + 9/L6 + 4/L8
# = 4*(L1+L8)/(L1*L8) + 9*(L2+L6)/(L2*L6)
# = 4*18/36 + 9*20/80
# = 2 + 9/4 = 2 + 2.25 = 4.25 = 17/4
exact_zb1 = 4*18/36 + 9*20/80
print(f"\n  zeta_broken(1) = 4/L1+9/L2+9/L6+4/L8")
print(f"    = 4*(L1+L8)/(L1*L8) + 9*(L2+L6)/(L2*L6)")
print(f"    = 4*18/36 + 9*20/80")
print(f"    = 2 + 9/4 = 17/4 = {17/4}")
print(f"    Numerical: {zb1:.10f}")
print(f"    Match: {abs(zb1-17/4)<1e-10}")

# zeta_total(1) = zeta_broken(1) + zeta_fixed(1)
zf1 = 16/9 + 25/12 + 36/14 + 16/15
print(f"\n  zeta_fixed(1) = 16/9+25/12+36/14+16/15 = {zf1:.10f}")
# Common denom: LCM(9,12,14,15) = 1260
# = 2240/1260 + 2625/1260 + 3240/1260 + 1344/1260
# = 9449/1260
from fractions import Fraction
zf1_exact = Fraction(16,9)+Fraction(25,12)+Fraction(36,14)+Fraction(16,15)
print(f"    = {zf1_exact} = {float(zf1_exact):.10f}")

zt1_exact = Fraction(17,4) + zf1_exact
print(f"\n  zeta_total(1) = 17/4 + {zf1_exact} = {zt1_exact} = {float(zt1_exact):.10f}")
print(f"    Numerical: {zeta_total(1):.10f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
  EXACT RESULTS AT s=1:
    zeta_broken(1) = 17/4 = {17/4}  (RATIONAL)
    zeta_asym(1)   = 67*sqrt5/60    (IN Q(sqrt5))
    zeta_total(1)  = {zt1_exact}

  GALOIS NORMS (EXACT):
    L1*L8 = 36 = b1^2
    L2*L6 = 80 = d_ST^2 * a1
    L8/L1 = phi^4

  EIGENVALUE DIFFERENCES (EXACT):
    L8-L1 = 6*sqrt5
    L6-L2 = 4*sqrt5

  zeta_norm(1) = 161/720 = 161/(N*b1)

  SCAN RESULTS:
    Multiple framework constants appear in the spectral zeta,
    but the s* values require further analysis for clean expressions.
    The ANALYTIC structure at s=1 is rich and tractable.
""")

print("=" * 70)
print("EXP-521 COMPLETE")
print("=" * 70)
