"""
EXP-241: UNIFIED MASS CORRECTIONS - NON-UNITS + GALOIS EXPONENT
================================================================

Key insight from exp099: sector corrections (d_eff, k_eff) already reduce
charm, top, bottom to sub-2%. The REAL remaining problem is down (10%)
and strange (7%), which are UNITS.

This experiment:
1. Confirms non-units are handled by sector corrections
2. Introduces n' = Galois conjugate exponent as new ingredient
3. Searches for a formula using n' that works for ALL quarks
4. Separates DERIVED vs FITTING honestly

GALOIS CONJUGATE EXPONENT:
  If z = a + b*phi has DSI level n = 5a + 6b,
  then z' = (a+b) - b*phi has DSI level n' = 5(a+b) + 6(-b) = 5a - b.

  So n' = 5a - b (just flip the sign of b and add 5b to a).
  For units: n + n' = 10a + 5b = 5*(2a+b). n - n' = 7b = (a_1+2)*b.
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-241: UNIFIED MASS CORRECTIONS WITH GALOIS EXPONENT")
print("=" * 70)

a_1 = 5
b_1 = 6
phi = PHI
alpha = ALPHA
alpha_s = 1 / (2 * phi**3)
sin2tW = 6.0 / 26.0
m_e = 0.51099895  # MeV

# PDG 2024 masses
masses_exp = {
    'electron': 0.51099895,
    'muon': 105.6583755,
    'tau': 1776.86,
    'up': 2.16,
    'down': 4.67,
    'strange': 93.4,
    'charm': 1270.0,
    'top': 172760.0,
    'bottom': 4180.0,
}

# Quantum numbers
qn = {
    'electron': (0, 0),
    'muon': (1, 1),
    'tau': (1, 2),
    'up': (3, -2),
    'down': (1, 0),
    'strange': (1, 1),
    'charm': (2, 1),
    'top': (4, 1),
    'bottom': (-1, 4),
}

sectors = {
    'electron': 'lepton', 'muon': 'lepton', 'tau': 'lepton',
    'up': 'up', 'charm': 'up', 'top': 'up',
    'down': 'down', 'strange': 'down', 'bottom': 'down',
}

# =====================================================================
# PART A: Compute all exponents and Galois conjugates
# =====================================================================
print("\n" + "=" * 70)
print("PART A: DSI EXPONENTS AND GALOIS CONJUGATES")
print("=" * 70)

print(f"\n{'Fermion':<10} {'(a,b)':>7} {'n':>4} {'n_prime':>7} {'N(z)':>5} "
      f"{'n-n_prime':>9} {'n+n_prime':>9} {'n_eff':>8} {'delta_n':>8}")
print("-" * 80)

data = {}
for name in ['electron', 'muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    a, b = qn[name]
    n = 5*a + 6*b
    n_prime = 5*a - b  # Galois conjugate exponent
    norm = a**2 + a*b - b**2  # Z[phi] norm

    m_exp = masses_exp[name]
    if m_exp > m_e:
        n_eff = np.log(m_exp / m_e) / np.log(phi)
    else:
        n_eff = 0.0
    delta_n = n_eff - n

    data[name] = {
        'a': a, 'b': b, 'n': n, 'n_prime': n_prime, 'norm': norm,
        'n_eff': n_eff, 'delta_n': delta_n, 'sector': sectors[name],
        'm_exp': m_exp, 'm_geom': m_e * phi**n
    }

    is_unit = abs(abs(norm) - 1) < 0.01
    unit_str = "U" if is_unit else f"{norm}"

    print(f"{name:<10} ({a:+d},{b:+d}) {n:4d} {n_prime:7d} {unit_str:>5s} "
          f"{n-n_prime:9d} {n+n_prime:9d} {n_eff:8.3f} {delta_n:+8.4f}")

# =====================================================================
# PART B: Verify sector corrections from exp099
# =====================================================================
print("\n" + "=" * 70)
print("PART B: SECTOR CORRECTIONS FROM EXP-099")
print("=" * 70)

# Sector-dependent effective d and k
sector_corrections = {
    'lepton': {'dd': sin2tW, 'dk': -1/phi**4},
    'up': {'dd': 2*alpha_s/np.pi, 'dk': alpha_s},
    'down': {'dd': -1.0/a_1, 'dk': -alpha_s},
}

print(f"\nSector corrections:")
for sec, corr in sector_corrections.items():
    d_eff = 5 + corr['dd']
    k_eff = 6 + corr['dk']
    print(f"  {sec:7s}: dd = {corr['dd']:+.4f}, dk = {corr['dk']:+.4f} "
          f"-> d_eff = {d_eff:.4f}, k_eff = {k_eff:.4f}")

print(f"\n{'Fermion':<10} {'n_bare':>6} {'n_sector':>8} {'m_sector':>10} "
      f"{'m_exp':>10} {'err_bare':>9} {'err_sector':>11}")
print("-" * 75)

for name in ['electron', 'muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    a, b = d['a'], d['b']
    sec = d['sector']
    corr = sector_corrections[sec]

    d_eff = 5 + corr['dd']
    k_eff = 6 + corr['dk']
    n_sector = a * d_eff + b * k_eff
    m_sector = m_e * phi**n_sector

    err_bare = (d['m_geom'] - d['m_exp']) / d['m_exp'] * 100
    err_sector = (m_sector - d['m_exp']) / d['m_exp'] * 100

    d['n_sector'] = n_sector
    d['m_sector'] = m_sector
    d['err_sector'] = err_sector
    d['err_bare'] = err_bare

    marker = " <-- SUB-2%" if abs(err_sector) < 2 else ""
    print(f"{name:<10} {d['n']:6d} {n_sector:8.3f} {m_sector:10.1f} "
          f"{d['m_exp']:10.1f} {err_bare:+8.1f}% {err_sector:+10.2f}%{marker}")

print(f"\nSUMMARY: Sector corrections solve non-units (c,t,b) to sub-2%!")
print(f"Remaining problems: down ({data['down']['err_sector']:+.1f}%), "
      f"strange ({data['strange']['err_sector']:+.1f}%)")

# =====================================================================
# PART C: Galois exponent n' as correction ingredient
# =====================================================================
print("\n" + "=" * 70)
print("PART C: GALOIS EXPONENT n' ANALYSIS")
print("=" * 70)

print(f"""
KEY OBSERVATION: n - n' = (a_1 + 2) * b = 7*b

This is EXACT. The difference between phi and phi' sectors is
proportional to the k-quantum number b.

For units: n_eff - n is small (< 0.4).
For non-units: n_eff - n can be large (up to 0.5).

HYPOTHESIS: The correction involves n' / b_1^2:
  n_eff = n + n' / b_1^2 = n + n'/36
""")

print(f"Test: n_eff = n + n'/36:")
print(f"{'Fermion':<10} {'n':>4} {'n_prime':>7} {'n+n_prime/36':>12} "
      f"{'n_eff_exp':>9} {'error':>8}")
print("-" * 55)

for name in ['charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + d['n_prime'] / 36.0
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    print(f"{name:<10} {d['n']:4d} {d['n_prime']:7d} {n_pred:12.3f} "
          f"{d['n_eff']:9.3f} {err:+7.2f}%")

# Now test for ALL fermions
print(f"\nExtended test for ALL fermions (n + n'/36):")
print(f"{'Fermion':<10} {'n+n_prime/36':>12} {'m_pred':>10} {'m_exp':>10} {'error':>8}")
print("-" * 55)

for name in ['electron', 'muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + d['n_prime'] / 36.0
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    marker = " <-- SUB-2%" if abs(err) < 2 else ""
    print(f"{name:<10} {n_pred:12.3f} {m_pred:10.1f} {d['m_exp']:10.1f} "
          f"{err:+7.2f}%{marker}")

# =====================================================================
# PART D: Optimized coefficient c in n + c*n'
# =====================================================================
print("\n" + "=" * 70)
print("PART D: OPTIMIZE COEFFICIENT c IN n + c*n'")
print("=" * 70)

# Scan c to minimize total error for non-units
best_c = 0
best_err = 1e10
results_scan = []

for c_int in range(0, 200):
    c = c_int / 1000.0
    total_sq = 0
    for name in ['charm', 'top', 'bottom']:
        d = data[name]
        n_pred = d['n'] + c * d['n_prime']
        m_pred = m_e * phi**n_pred
        err = (m_pred - d['m_exp']) / d['m_exp'] * 100
        total_sq += err**2
    results_scan.append((c, total_sq))
    if total_sq < best_err:
        best_err = total_sq
        best_c = c

print(f"Best c for non-units only: c = {best_c:.4f}")
print(f"  1/b_1^2 = 1/36 = {1/36:.4f}")
print(f"  alpha_s/(2*pi) = {alpha_s/(2*np.pi):.4f}")
print(f"  alpha_s/a_1 = {alpha_s/a_1:.4f}")
print(f"  1/(a_1*N_eig) = {1/(a_1*9):.4f}")

print(f"\nWith c = {best_c:.4f}:")
for name in ['charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + best_c * d['n_prime']
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    print(f"  {name}: n_eff = {n_pred:.3f}, m = {m_pred:.1f}, err = {err:+.2f}%")

# Also optimize for ALL quarks (not just non-units)
print(f"\nOptimizing c for ALL 6 quarks:")
best_c_all = 0
best_err_all = 1e10

for c_int in range(-100, 200):
    c = c_int / 1000.0
    total_sq = 0
    for name in ['up', 'down', 'strange', 'charm', 'top', 'bottom']:
        d = data[name]
        n_pred = d['n'] + c * d['n_prime']
        m_pred = m_e * phi**n_pred
        err = (m_pred - d['m_exp']) / d['m_exp'] * 100
        total_sq += err**2
    if total_sq < best_err_all:
        best_err_all = total_sq
        best_c_all = c

print(f"Best c for all quarks: c = {best_c_all:.4f}")
print(f"\nWith c = {best_c_all:.4f}:")
for name in ['up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + best_c_all * d['n_prime']
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    print(f"  {name}: n_eff = {n_pred:.3f}, m = {m_pred:.1f}, err = {err:+.2f}%")

# =====================================================================
# PART E: Two-parameter fit: n_eff = n + c1*n' + c2*N(z)
# =====================================================================
print("\n" + "=" * 70)
print("PART E: TWO-PARAMETER FORMULA n + c1*n' + c2*N(z)")
print("=" * 70)

# For 6 quarks, minimize sum of squares
best_params = (0, 0)
best_err_2p = 1e10

for c1_int in range(-100, 200):
    c1 = c1_int / 1000.0
    for c2_int in range(-100, 100):
        c2 = c2_int / 1000.0
        total_sq = 0
        for name in ['up', 'down', 'strange', 'charm', 'top', 'bottom']:
            d = data[name]
            n_pred = d['n'] + c1 * d['n_prime'] + c2 * d['norm']
            m_pred = m_e * phi**n_pred
            err = (m_pred - d['m_exp']) / d['m_exp'] * 100
            total_sq += err**2
        if total_sq < best_err_2p:
            best_err_2p = total_sq
            best_params = (c1, c2)

c1_best, c2_best = best_params
print(f"Best (c1, c2) = ({c1_best:.4f}, {c2_best:.4f})")
print(f"  c1 ~ 1/{1/c1_best:.1f}" if c1_best > 0.001 else "  c1 ~ 0")
if abs(c2_best) > 0.001:
    print(f"  c2 ~ 1/{1/c2_best:.1f}" if c2_best != 0 else "  c2 ~ 0")

print(f"\nWith c1={c1_best:.4f}, c2={c2_best:.4f}:")
for name in ['up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + c1_best * d['n_prime'] + c2_best * d['norm']
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    marker = " <-- SUB-2%" if abs(err) < 2 else ""
    print(f"  {name}: n_eff = {n_pred:.3f}, m = {m_pred:.1f}, err = {err:+.2f}%{marker}")

# =====================================================================
# PART F: Sector + Galois combined
# =====================================================================
print("\n" + "=" * 70)
print("PART F: SECTOR CORRECTIONS + GALOIS REFINEMENT")
print("=" * 70)

print("""
Can we combine the sector corrections (which work for non-units)
with a Galois refinement (which might fix down/strange)?

Formula: n_eff = a*(5+dd) + b*(6+dk) + c*n'

This adds a Galois term c*n' to the sector-corrected formula.
""")

# For each sector, optimize c
for sec_name in ['lepton', 'up', 'down']:
    fermion_names = [n for n in data if data[n]['sector'] == sec_name]

    corr = sector_corrections[sec_name]
    d_eff = 5 + corr['dd']
    k_eff = 6 + corr['dk']

    best_c_sec = 0
    best_err_sec = 1e10

    for c_int in range(-200, 200):
        c = c_int / 10000.0
        total_sq = 0
        for name in fermion_names:
            d = data[name]
            n_pred = d['a'] * d_eff + d['b'] * k_eff + c * d['n_prime']
            m_pred = m_e * phi**n_pred
            err = ((m_pred - d['m_exp']) / d['m_exp'] * 100)
            total_sq += err**2
        if total_sq < best_err_sec:
            best_err_sec = total_sq
            best_c_sec = c

    print(f"\n{sec_name.upper()} sector: best c = {best_c_sec:.5f}")
    for name in fermion_names:
        d = data[name]
        n_pred = d['a'] * d_eff + d['b'] * k_eff + best_c_sec * d['n_prime']
        m_pred = m_e * phi**n_pred
        err = (m_pred - d['m_exp']) / d['m_exp'] * 100
        print(f"  {name}: err = {err:+.2f}%")

# =====================================================================
# PART G: The Galois embedding formula
# =====================================================================
print("\n" + "=" * 70)
print("PART G: GALOIS EMBEDDING MASS FORMULA")
print("=" * 70)

print("""
PHYSICAL IDEA: The mass is determined by the Z[phi] element z = a+b*phi
through BOTH embeddings:

  z  -> R (phi-sector):  |z| = |a + b*phi|
  z' -> R (phi'-sector): |z'| = |a + b*phi'|

For UNITS: |z|*|z'| = 1, so knowing |z| suffices.
For NON-UNITS: |z|*|z'| = |N(z)| != 1, second embedding matters.

FORMULA: m = m_e * phi^n * |N(z)|^c

For units: |N| = 1, correction = 1 (no change).
""")

# Test: m = m_e * phi^n * |N(z)|^c
print(f"Testing m = m_e * phi^n * |N(z)|^c for non-units:")

# What c is needed for each?
for name in ['charm', 'top', 'bottom']:
    d = data[name]
    needed_factor = d['m_exp'] / d['m_geom']
    abs_norm = abs(d['norm'])
    if abs_norm > 1:
        c_needed = np.log(needed_factor) / np.log(abs_norm)
        print(f"  {name}: |N|={abs_norm}, factor={needed_factor:.4f}, c_needed={c_needed:.5f}")

# The sign of N matters - try signed version
print(f"\nSigned version: m = m_e * phi^n * |N(z)|^(sign(N)*c)")
for name in ['charm', 'top', 'bottom']:
    d = data[name]
    sgn = 1 if d['norm'] > 0 else -1
    needed_factor = d['m_exp'] / d['m_geom']
    abs_norm = abs(d['norm'])
    c_needed = np.log(needed_factor) / (sgn * np.log(abs_norm))
    print(f"  {name}: N={d['norm']:+d}, sign={sgn:+d}, c_needed={c_needed:.5f}")

# =====================================================================
# PART H: n' as the KEY - the Galois mass formula
# =====================================================================
print("\n" + "=" * 70)
print("PART H: n' AS THE KEY INGREDIENT")
print("=" * 70)

print("""
Let's look at what n' represents physically:

  n  = 5a + 6b  = DSI level in phi-sector
  n' = 5a - b   = DSI level in phi'-sector

For the 600-cell, the phi and phi' sectors contribute differently.
The mass in the phi-sector is phi^n, in the phi'-sector would be phi'^{n'}.

Since phi' < 0, we use |phi'|^{n'} = phi^{-n'} (since |phi'| = 1/phi).

The GALOIS-WEIGHTED mass would be:
  m = m_e * phi^n * (phi'/phi)^{c*n'}

But phi'/phi = -1/phi^2, so |phi'/phi|^{n'} = phi^{-2*n'}.

A simpler version: m = m_e * phi^{n + c*n'} for small c.

Let me check what c value matches the n' formula for EACH fermion.
""")

print(f"Required c = delta_n / n' for each fermion:")
print(f"{'Fermion':<10} {'delta_n':>9} {'n_prime':>7} {'c_req':>8} {'N(z)':>5}")
print("-" * 45)

for name in ['muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    if d['n_prime'] != 0:
        c_req = d['delta_n'] / d['n_prime']
    else:
        c_req = float('inf')
    print(f"{name:<10} {d['delta_n']:+9.4f} {d['n_prime']:7d} {c_req:+8.5f} {d['norm']:5d}")

# =====================================================================
# PART I: Sector-specific c values for n'
# =====================================================================
print("\n" + "=" * 70)
print("PART I: SECTOR-SPECIFIC GALOIS COEFFICIENTS")
print("=" * 70)

print("""
If c depends on SECTOR (lepton/up/down), what are the best c values?
""")

for sec_name in ['lepton', 'up', 'down']:
    fermion_names = [n for n in data if data[n]['sector'] == sec_name and n != 'electron']

    if not fermion_names:
        continue

    best_c_sec = 0
    best_err_sec = 1e10

    for c_int in range(-500, 500):
        c = c_int / 10000.0
        total_sq = 0
        for name in fermion_names:
            d = data[name]
            n_pred = d['n'] + c * d['n_prime']
            m_pred = m_e * phi**n_pred
            err = (m_pred - d['m_exp']) / d['m_exp'] * 100
            total_sq += err**2
        if total_sq < best_err_sec:
            best_err_sec = total_sq
            best_c_sec = c

    print(f"\n{sec_name.upper()} sector: best c = {best_c_sec:.5f}")

    # Check if c matches any framework constant
    if abs(best_c_sec) > 1e-5:
        print(f"  Closest matches:")
        candidates = [
            ('alpha', alpha),
            ('alpha_s', alpha_s),
            ('alpha_s/pi', alpha_s/np.pi),
            ('alpha/(2*pi)', alpha/(2*np.pi)),
            ('1/b_1^2', 1/b_1**2),
            ('1/(a_1*N_eig)', 1/(a_1*9)),
            ('alpha_s/a_1', alpha_s/a_1),
            ('sin2tW/a_1^2', sin2tW/a_1**2),
            ('-alpha_s/a_1', -alpha_s/a_1),
            ('-1/b_1^2', -1/b_1**2),
        ]
        for label, val in candidates:
            if abs(val) > 1e-6:
                ratio = best_c_sec / val
                if 0.5 < abs(ratio) < 2.0:
                    print(f"    {label} = {val:.5f} (ratio = {ratio:.3f})")

    for name in fermion_names:
        d = data[name]
        n_pred = d['n'] + best_c_sec * d['n_prime']
        m_pred = m_e * phi**n_pred
        err = (m_pred - d['m_exp']) / d['m_exp'] * 100
        print(f"  {name}: err = {err:+.2f}%")

# =====================================================================
# PART J: The UNIFIED formula attempt
# =====================================================================
print("\n" + "=" * 70)
print("PART J: UNIFIED FORMULA - ALL 9 FERMIONS")
print("=" * 70)

print("""
ATTEMPT: A single formula for all 9 fermions:
  m = m_e * phi^{n + c_sector * n'}

where c_sector depends on the sector (lepton/up/down).

Can we DERIVE c_sector from the framework?
""")

# Try a UNIVERSAL formula: n_eff = n + c_universal * n' * alpha_s
print(f"Test: n_eff = n + alpha_s * n' / (a_1 + 1)")
c_test = alpha_s / (a_1 + 1)
print(f"  c = alpha_s / b_1 = {c_test:.5f}")

for name in ['electron', 'muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + c_test * d['n_prime']
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    marker = " <-- SUB-2%" if abs(err) < 2 else ""
    print(f"  {name:<10}: n_eff = {n_pred:.3f}, err = {err:+.2f}%{marker}")

# Try: c = alpha_s / (2*pi)
print(f"\nTest: c = alpha_s / (2*pi) = {alpha_s/(2*np.pi):.5f}")
c_test2 = alpha_s / (2*np.pi)

for name in ['electron', 'muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    n_pred = d['n'] + c_test2 * d['n_prime']
    m_pred = m_e * phi**n_pred
    err = (m_pred - d['m_exp']) / d['m_exp'] * 100
    marker = " <-- SUB-2%" if abs(err) < 2 else ""
    print(f"  {name:<10}: n_eff = {n_pred:.3f}, err = {err:+.2f}%{marker}")

# =====================================================================
# PART K: QCD running from Galois scale
# =====================================================================
print("\n" + "=" * 70)
print("PART K: QCD RUNNING FROM GALOIS SCALE")
print("=" * 70)

print("""
PHYSICAL MECHANISM: The phi-sector mass m_e*phi^n is the mass at the
phi-sector scale. QCD running from this scale to the physical scale
gives a correction proportional to n' (the phi'-sector exponent).

The running factor is:
  m_phys / m_geom = [alpha_s(m_phys) / alpha_s(m_geom)]^{gamma_0/(2*beta_0)}

For non-units, the geometric and physical scales are DIFFERENT.
For units, they coincide (approximately).
""")

# Alpha_s running
alpha_s_mZ = 0.1179
M_Z = 91187.6

def alpha_s_run(mu, n_f):
    beta_0 = 11 - 2*n_f/3.0
    return alpha_s_mZ / (1 + beta_0/(2*np.pi) * alpha_s_mZ * np.log(mu/M_Z))

def running_ratio(mu1, mu2, n_f):
    gamma_0 = 8.0
    beta_0 = 11 - 2*n_f/3.0
    return (alpha_s_run(mu1, n_f) / alpha_s_run(mu2, n_f))**(gamma_0/(2*beta_0))

# Test: run from phi^n-related scale to physical scale
# The "Galois scale" could be m_e * phi^{|n'|}
print(f"Running from Galois scale m_e * phi^|n'| to physical scale:")
for name in ['charm', 'top', 'bottom']:
    d = data[name]
    galois_scale = m_e * phi**abs(d['n_prime'])  # MeV
    phys_scale = d['m_exp']
    n_f = 4 if name == 'charm' else (6 if name == 'top' else 5)

    ratio = running_ratio(phys_scale, galois_scale, n_f)
    m_pred = d['m_geom'] * ratio
    err = (m_pred - d['m_exp'])/d['m_exp']*100

    print(f"  {name}: galois_scale = {galois_scale:.0f} MeV, "
          f"run_ratio = {ratio:.4f}, m = {m_pred:.0f}, err = {err:+.1f}%")

# =====================================================================
# PART L: Summary table
# =====================================================================
print("\n" + "=" * 70)
print("PART L: COMPREHENSIVE COMPARISON")
print("=" * 70)

print(f"\n{'Method':<25} {'charm':>8} {'top':>8} {'bottom':>8} {'down':>8} {'strange':>8}")
print("-" * 70)

# Method 1: Bare phi^n
errs = {}
for name in ['charm', 'top', 'bottom', 'down', 'strange']:
    d = data[name]
    errs[name] = d['err_bare']
print(f"{'Bare phi^n':<25} {errs['charm']:+7.1f}% {errs['top']:+7.1f}% "
      f"{errs['bottom']:+7.1f}% {errs['down']:+7.1f}% {errs['strange']:+7.1f}%")

# Method 2: Sector corrections
for name in ['charm', 'top', 'bottom', 'down', 'strange']:
    d = data[name]
    errs[name] = d['err_sector']
print(f"{'Sector (exp099)':<25} {errs['charm']:+7.1f}% {errs['top']:+7.1f}% "
      f"{errs['bottom']:+7.1f}% {errs['down']:+7.1f}% {errs['strange']:+7.1f}%")

# Method 3: n + c*n' with best c for all quarks
for name in ['charm', 'top', 'bottom', 'down', 'strange']:
    d = data[name]
    n_pred = d['n'] + best_c_all * d['n_prime']
    m_pred = m_e * phi**n_pred
    errs[name] = (m_pred - d['m_exp'])/d['m_exp']*100
print(f"{'n + c*n_prime (c=' + f'{best_c_all:.3f})':<25} {errs['charm']:+7.1f}% {errs['top']:+7.1f}% "
      f"{errs['bottom']:+7.1f}% {errs['down']:+7.1f}% {errs['strange']:+7.1f}%")

# Method 4: n + c1*n' + c2*N(z)
for name in ['charm', 'top', 'bottom', 'down', 'strange']:
    d = data[name]
    n_pred = d['n'] + c1_best * d['n_prime'] + c2_best * d['norm']
    m_pred = m_e * phi**n_pred
    errs[name] = (m_pred - d['m_exp'])/d['m_exp']*100
label4 = f'n+c1*n_prime+c2*N'
print(f"{label4:<25} {errs['charm']:+7.1f}% {errs['top']:+7.1f}% "
      f"{errs['bottom']:+7.1f}% {errs['down']:+7.1f}% {errs['strange']:+7.1f}%")

# =====================================================================
# PART M: Key identities involving n'
# =====================================================================
print("\n" + "=" * 70)
print("PART M: GALOIS EXPONENT IDENTITIES")
print("=" * 70)

print(f"\nExact identities (from n = 5a+6b, n' = 5a-b):")
print(f"  n - n' = 7b = (a_1+2)*b   [EXACT for all fermions]")
print(f"  n + n' = 10a + 5b = 5*(2a+b)")
print(f"  n*n' = 25a^2 + 25ab - 6b^2 = 25*N(z) + 19*b^2")

print(f"\nFermion Galois products n*n':")
for name in ['muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    prod = d['n'] * d['n_prime']
    check = 25*d['norm'] + 19*d['b']**2
    print(f"  {name:<10}: n*n' = {prod:5d} = 25*{d['norm']:+3d} + 19*{d['b']**2:2d} = {check:5d}")

print(f"\nCross-sector sums n + n':")
for name in ['muon', 'tau', 'up', 'down', 'strange', 'charm', 'top', 'bottom']:
    d = data[name]
    s = d['n'] + d['n_prime']
    print(f"  {name:<10}: n+n' = {s:3d} = 5*{s//5}", end="")
    if s % 5 == 0:
        val = s // 5
        # Check what 2a+b is
        print(f" = 5*(2*{d['a']}+{d['b']}) = 5*{2*d['a']+d['b']}")
    else:
        print()

# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY AND CONCLUSIONS")
print("=" * 70)

print("""
MAIN FINDINGS:
==============

1. SECTOR CORRECTIONS (exp099) ALREADY SOLVE NON-UNITS:
   - charm: 1.0% (was -11%)
   - top: -1.8% (was -20%)
   - bottom: +0.3% (was +14%)
   All three non-unit fermions are SUB-2% with sector corrections!

2. The REAL remaining problem is DOWN (10%) and STRANGE (7%),
   which are UNITS in Z[phi]. The non-unit problem was a RED HERRING.

3. GALOIS EXPONENT n' = 5a - b provides a new structural ingredient:
   - n - n' = (a_1+2)*b  (EXACT, algebraic identity)
   - For non-units, n' encodes the phi'-sector contribution

4. A universal Galois correction n_eff = n + c*n' with c ~ 1/b_1^2
   improves non-units but does NOT solve down/strange.

5. Down quark (n=5, n'=5): UNIQUELY has n = n' (Galois self-conjugate!)
   This means it sits ON the fixed point of the Galois action.
   Its mass correction cannot come from the n' mechanism.

6. Strange quark shares (a,b) = (1,1) with muon. The 7% error is the
   QCD correction for strange at 2 GeV relative to the muon mass.

HONEST STATUS:
  - Non-units (c,t,b): SOLVED by sector corrections (exp099)
  - Down quark: OPEN (10%, within 1-sigma of PDG uncertainty)
  - Strange quark: OPEN (7%, QCD running partially explains)
  - Formula m_f = m_e * phi^{a*d_eff + b*k_eff}: 7/9 sub-2%
""")

print("\n" + "=" * 70)
print("EXP-241 COMPLETE")
print("=" * 70)
