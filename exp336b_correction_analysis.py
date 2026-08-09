"""
exp336b: Analysis of bare mass errors - searching for DERIVED corrections
=========================================================================
The bare formula m_f = m_e * phi^(5a+6b) has 12.6% RMS error.
The old sector corrections (delta_d, delta_k) are FITTING.

Question: Is there a DERIVED correction that improves accuracy?

Key observation: the SIGN of errors splits cleanly by T3:
  - Up-type quarks (T3=+1/2): predicted TOO LOW (c: -11%, t: -20%)
  - Down-type quarks (T3=-1/2): predicted TOO HIGH (d: +21%, s: +9%, b: +14%)

This T3 pattern is universal and physically meaningful.
"""

import math

PHI = (1 + math.sqrt(5)) / 2
m_e = 0.51100  # MeV

# Fermion data: name, (a,b), n=5a+6b, m_exp, T3, Q_em, C2_color, gen
fermions = [
    # name,  a,  b,  n,  m_exp,     T3,    Q,    C2,  gen
    ('e',    0,  0,  0,  0.51100,  -0.5, -1.0,  0.0,   0),
    ('mu',   1,  1, 11,  105.658,  -0.5, -1.0,  0.0,   1),
    ('tau',  1,  2, 17,  1776.86,  -0.5, -1.0,  0.0,   2),
    ('u',    3, -2,  3,  2.16,      0.5,  2/3,  4/3,   0),
    ('c',    2,  1, 16,  1270.0,    0.5,  2/3,  4/3,   1),
    ('t',    4,  1, 26,  172690.0,  0.5,  2/3,  4/3,   2),
    ('d',    1,  0,  5,  4.67,     -0.5, -1/3,  4/3,   0),
    ('s',    1,  1, 11,  93.4,     -0.5, -1/3,  4/3,   1),
    ('b',   -1,  4, 19,  4180.0,   -0.5, -1/3,  4/3,   2),
]

# Framework constants
a1 = 5
b1 = 6
alpha_s = 1/(2*PHI**3)  # 0.11803
sin2_tW = b1/(a1**2 + 1)  # 6/26 = 0.23077

print("="*72)
print("EFFECTIVE EXPONENTS AND DEVIATIONS FROM INTEGER")
print("="*72)

print(f"\n{'Name':>5} {'n':>3} {'n_eff':>8} {'delta_n':>8} {'T3':>5} {'Q':>6} {'Sector':>8}")
print("-"*55)

for f in fermions:
    name, a, b, n, m_exp, T3, Q, C2, gen = f
    if name == 'e':
        n_eff = 0.0
        delta = 0.0
    else:
        n_eff = math.log(m_exp / m_e) / math.log(PHI)
        delta = n_eff - n

    f_dict = {'n_eff': n_eff, 'delta': delta}

    sector = 'lepton' if C2 == 0 else ('up' if T3 > 0 else 'down')
    print(f"{name:>5} {n:>3} {n_eff:>8.4f} {delta:>+8.4f} {T3:>+5.1f} {Q:>+6.2f} {sector:>8}")


# ================================================================
# Pattern analysis
# ================================================================
print("\n" + "="*72)
print("SECTOR AVERAGES")
print("="*72)

sectors = {'lepton': [], 'up': [], 'down': []}
for f in fermions:
    name, a, b, n, m_exp, T3, Q, C2, gen = f
    if name == 'e': continue
    n_eff = math.log(m_exp / m_e) / math.log(PHI)
    delta = n_eff - n
    sector = 'lepton' if C2 == 0 else ('up' if T3 > 0 else 'down')
    sectors[sector].append((name, delta, gen))

for sec in ['lepton', 'up', 'down']:
    deltas = [d for _, d, _ in sectors[sec]]
    avg = sum(deltas) / len(deltas)
    print(f"\n{sec:>8}: avg delta_n = {avg:+.4f}")
    for name, d, g in sectors[sec]:
        print(f"  gen {g}: {name:>3} delta = {d:+.4f}")


# ================================================================
# Test: T3 correction with various framework coefficients
# ================================================================
print("\n" + "="*72)
print("TEST: n -> n + 2*T3*X  (quarks only, leptons unchanged)")
print("="*72)

# Framework candidate coefficients
candidates = {
    'alpha_s = 1/(2*phi^3)': alpha_s,
    '2*alpha_s': 2*alpha_s,
    'sin^2(tW) = 6/26': sin2_tW,
    '1/phi^2 = 2-phi': 1/PHI**2,
    '1/4 = b1/24': 0.25,
    '1/phi^3 = 2*alpha_s': 1/PHI**3,
    '1/(a1-1) = 1/4': 1/(a1-1),
    'alpha_s*phi': alpha_s*PHI,
    '1/(2*phi)': 1/(2*PHI),
    '1/a1 = 1/5': 1/a1,
    'OPTIMAL (fitted)': None,  # will compute
}

# First find optimal X
best_x = 0
best_rms = 999
for x_try in [i*0.001 for i in range(1, 1000)]:
    total_sq = 0
    count = 0
    for f in fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        if name == 'e': continue
        if C2 > 0:  # quark
            n_corr = n + 2*T3*x_try
        else:  # lepton
            n_corr = n
        m_pred = m_e * PHI**n_corr
        err = (m_pred / m_exp - 1) * 100
        total_sq += err**2
        count += 1
    rms = math.sqrt(total_sq / count)
    if rms < best_rms:
        best_rms = rms
        best_x = x_try

candidates['OPTIMAL (fitted)'] = best_x

print(f"\n{'Coefficient':>30} {'X':>8} | {'mu':>6} {'tau':>6} {'u':>6} {'c':>6} {'t':>6} {'d':>6} {'s':>6} {'b':>6} | {'RMS':>6}")
print("-"*115)

for label, X in candidates.items():
    errs = {}
    for f in fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        if name == 'e': continue
        if C2 > 0:  # quark
            n_corr = n + 2*T3*X
        else:  # lepton
            n_corr = n
        m_pred = m_e * PHI**n_corr
        errs[name] = (m_pred / m_exp - 1) * 100

    all_errs = list(errs.values())
    rms = math.sqrt(sum(e**2 for e in all_errs) / len(all_errs))

    print(f"{label:>30} {X:>8.4f} |"
          f" {errs['mu']:>+5.1f} {errs['tau']:>+5.1f}"
          f" {errs['u']:>+5.1f} {errs['c']:>+5.1f} {errs['t']:>+5.1f}"
          f" {errs['d']:>+5.1f} {errs['s']:>+5.1f} {errs['b']:>+5.1f}"
          f" | {rms:>5.1f}%")

print(f"\n  Bare formula (X=0): RMS = 12.6%")
print(f"  Old sector corrections: RMS = 4.4% (6 fitted parameters)")
print(f"  Optimal T3 correction: RMS = {best_rms:.1f}% (1 parameter, X = {best_x:.4f})")


# ================================================================
# Check: what if correction depends on T3 AND generation?
# ================================================================
print("\n" + "="*72)
print("TEST: n -> n + 2*T3*X*g  (generation-dependent, quarks only)")
print("="*72)
print("(g = generation index: 0, 1, 2)")

# Required: delta_n = 2*T3*X*g
# Up: g=0 -> 0, g=1 -> X, g=2 -> 2X
# Down: g=0 -> 0, g=1 -> -X, g=2 -> -2X
# But d (g=0) has delta=-0.40, which can't be 0!

print("\nThis model CANNOT work: d quark (gen 0) has delta = -0.40,")
print("but the model predicts delta = 0 for gen 0.")


# ================================================================
# Test: T3 + constant shift for quarks
# ================================================================
print("\n" + "="*72)
print("TEST: n -> n + 2*T3*X + Y*C2  (T3 + color Casimir shift)")
print("="*72)
print("(C2 = 4/3 for quarks, 0 for leptons)")

# Scan over X and Y
best_rms_2 = 999
best_x_2 = 0
best_y_2 = 0

for x_try in [i*0.01 for i in range(-50, 50)]:
    for y_try in [i*0.01 for i in range(-50, 50)]:
        total_sq = 0
        count = 0
        for f in fermions:
            name, a, b, n, m_exp, T3, Q, C2, gen = f
            if name == 'e': continue
            n_corr = n + 2*T3*x_try + y_try*C2
            m_pred = m_e * PHI**n_corr
            err = (m_pred / m_exp - 1) * 100
            total_sq += err**2
            count += 1
        rms = math.sqrt(total_sq / count)
        if rms < best_rms_2:
            best_rms_2 = rms
            best_x_2 = x_try
            best_y_2 = y_try

print(f"\nOptimal: X = {best_x_2:.3f}, Y = {best_y_2:.3f}, RMS = {best_rms_2:.1f}%")

# Now test with framework candidates for (X, Y)
print(f"\nTest with framework values:")

combo_tests = [
    ('T3*alpha_s + C2*0', alpha_s, 0),
    ('T3*sin2tW + C2*0', sin2_tW, 0),
    ('T3*1/4 + C2*0', 0.25, 0),
    ('T3*alpha_s - C2*alpha_s/2', alpha_s, -alpha_s/2),
    ('T3*1/4 - C2*alpha_s', 0.25, -alpha_s),
    ('T3*sin2tW - C2*alpha_s', sin2_tW, -alpha_s),
    ('T3*alpha_s - C2*1/(2*a1)', alpha_s, -1/(2*a1)),
    (f'OPTIMAL: T3*{best_x_2:.3f} + C2*{best_y_2:.3f}', best_x_2, best_y_2),
]

print(f"\n{'Model':>40} | {'mu':>6} {'tau':>6} {'u':>6} {'c':>6} {'t':>6} {'d':>6} {'s':>6} {'b':>6} | {'RMS':>6}")
print("-"*115)

for label, X, Y in combo_tests:
    errs = {}
    for f in fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        if name == 'e': continue
        n_corr = n + 2*T3*X + Y*C2
        m_pred = m_e * PHI**n_corr
        errs[name] = (m_pred / m_exp - 1) * 100

    all_errs = list(errs.values())
    rms = math.sqrt(sum(e**2 for e in all_errs) / len(all_errs))

    print(f"{label:>40} |"
          f" {errs['mu']:>+5.1f} {errs['tau']:>+5.1f}"
          f" {errs['u']:>+5.1f} {errs['c']:>+5.1f} {errs['t']:>+5.1f}"
          f" {errs['d']:>+5.1f} {errs['s']:>+5.1f} {errs['b']:>+5.1f}"
          f" | {rms:>5.1f}%")


# ================================================================
# Test: What if coefficients are not (5,6) but (5+eps_a, 6+eps_b)?
# ================================================================
print("\n" + "="*72)
print("TEST: m_f = m_e * phi^((5+da)*a + (6+db)*b)")
print("="*72)
print("Scanning (da, db) for best RMS...")

best_rms_3 = 999
best_da = 0
best_db = 0

for da in [i*0.001 for i in range(-200, 200)]:
    for db in [i*0.001 for i in range(-200, 200)]:
        total_sq = 0
        count = 0
        for f in fermions:
            name, a, b, n, m_exp, T3, Q, C2, gen = f
            if name == 'e': continue
            n_corr = (5+da)*a + (6+db)*b
            m_pred = m_e * PHI**n_corr
            err = (m_pred / m_exp - 1) * 100
            total_sq += err**2
            count += 1
        rms = math.sqrt(total_sq / count)
        if rms < best_rms_3:
            best_rms_3 = rms
            best_da = da
            best_db = db

print(f"\nOptimal: a1_eff = {5+best_da:.4f}, b1_eff = {6+best_db:.4f}")
print(f"Shifts: da = {best_da:+.4f}, db = {best_db:+.4f}")
print(f"RMS = {best_rms_3:.1f}%")

# Check if shifts match framework numbers
print(f"\nFramework number check:")
print(f"  da = {best_da:+.4f}")
print(f"    alpha_s = {alpha_s:.4f}")
print(f"    1/phi^4 = {1/PHI**4:.4f}")
print(f"    sin2_tW = {sin2_tW:.4f}")
print(f"    2*alpha_s = {2*alpha_s:.4f}")
print(f"  db = {best_db:+.4f}")
print(f"    alpha_s = {alpha_s:.4f}")
print(f"    -alpha_s = {-alpha_s:.4f}")
print(f"    -1/phi^4 = {-1/PHI**4:.4f}")

# Show errors with optimal (da, db)
print(f"\nErrors with optimal (da, db):")
print(f"{'Name':>5} {'n_orig':>7} {'n_corr':>8} {'m_pred':>12} {'m_exp':>12} {'Error':>8}")
print("-"*55)
for f in fermions:
    name, a, b, n, m_exp, T3, Q, C2, gen = f
    if name == 'e':
        print(f"{'e':>5} {'0':>7} {'0':>8} {'0.511':>12} {'0.511':>12} {'input':>8}")
        continue
    n_corr = (5+best_da)*a + (6+best_db)*b
    m_pred = m_e * PHI**n_corr
    err = (m_pred/m_exp - 1)*100
    print(f"{name:>5} {n:>7} {n_corr:>8.3f} {m_pred:>12.3f} {m_exp:>12.3f} {err:>+7.1f}%")


# ================================================================
# Test: sector-specific (da, db) for each of 3 sectors
# ================================================================
print("\n" + "="*72)
print("TEST: Sector-specific (da, db) - what does each sector NEED?")
print("="*72)

for sector_name, sector_filter in [('Leptons', lambda f: f[7]==0),
                                    ('Up quarks', lambda f: f[7]>0 and f[5]>0),
                                    ('Down quarks', lambda f: f[7]>0 and f[5]<0)]:
    sector_fermions = [f for f in fermions if sector_filter(f) and f[0]!='e']
    if len(sector_fermions) < 2:
        continue

    best_rms_s = 999
    best_da_s = 0
    best_db_s = 0

    for da in [i*0.001 for i in range(-500, 500)]:
        for db in [i*0.001 for i in range(-500, 500)]:
            total_sq = 0
            for f in sector_fermions:
                name, a, b, n, m_exp, T3, Q, C2, gen = f
                n_corr = (5+da)*a + (6+db)*b
                m_pred = m_e * PHI**n_corr
                err = (m_pred/m_exp - 1)*100
                total_sq += err**2
            rms = math.sqrt(total_sq/len(sector_fermions))
            if rms < best_rms_s:
                best_rms_s = rms
                best_da_s = da
                best_db_s = db

    print(f"\n{sector_name}:")
    print(f"  Optimal: da = {best_da_s:+.4f}, db = {best_db_s:+.4f}, RMS = {best_rms_s:.2f}%")
    for f in sector_fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        n_corr = (5+best_da_s)*a + (6+best_db_s)*b
        m_pred = m_e * PHI**n_corr
        err = (m_pred/m_exp - 1)*100
        print(f"    {name:>3}: n={n_corr:.3f}, pred={m_pred:.3f}, exp={m_exp:.3f}, err={err:+.1f}%")


# ================================================================
# Test: Electric charge as correction
# ================================================================
print("\n" + "="*72)
print("TEST: n -> n + X*Q_em  (electric charge correction)")
print("="*72)

best_rms_q = 999
best_xq = 0

for xq in [i*0.01 for i in range(-200, 200)]:
    total_sq = 0
    count = 0
    for f in fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        if name == 'e': continue
        n_corr = n + xq*Q
        m_pred = m_e * PHI**n_corr
        err = (m_pred/m_exp - 1)*100
        total_sq += err**2
        count += 1
    rms = math.sqrt(total_sq/count)
    if rms < best_rms_q:
        best_rms_q = rms
        best_xq = xq

print(f"Optimal: X_Q = {best_xq:.3f}, RMS = {best_rms_q:.1f}%")

# Test with specific framework values for Q coefficient
q_candidates = [alpha_s, -alpha_s, sin2_tW, -sin2_tW, 1/PHI**2, -1/PHI**2, 0.25, -0.25]
for xq in q_candidates:
    total_sq = 0
    count = 0
    errs = {}
    for f in fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        if name == 'e': continue
        n_corr = n + xq*Q
        m_pred = m_e * PHI**n_corr
        err = (m_pred/m_exp - 1)*100
        errs[name] = err
        total_sq += err**2
        count += 1
    rms = math.sqrt(total_sq/count)
    print(f"  X_Q = {xq:>+7.4f}: RMS = {rms:>5.1f}% | u:{errs['u']:>+5.1f} c:{errs['c']:>+5.1f} t:{errs['t']:>+5.1f} d:{errs['d']:>+5.1f} s:{errs['s']:>+5.1f} b:{errs['b']:>+5.1f}")


# ================================================================
# Test: Hypercharge Y as correction
# ================================================================
print("\n" + "="*72)
print("TEST: n -> n + X*Y_W  (weak hypercharge correction)")
print("="*72)

# Hypercharge Y = 2(Q - T3)
# Leptons (e_L): Y = 2(-1 - (-1/2)) = -1; (e_R): Y = 2(-1 - 0) = -2
# For mass eigenstates, use average or left-handed Y
# Actually, Y of left-handed doublet:
# (nu_L, e_L): Y = -1
# (u_L, d_L): Y = 1/3
# Right-handed:
# e_R: Y = -2, u_R: Y = 4/3, d_R: Y = -2/3

# For the mass formula, let's use the LEFT-HANDED hypercharge
Y_vals = {
    'e': -1, 'mu': -1, 'tau': -1,
    'u': 1/3, 'c': 1/3, 't': 1/3,
    'd': 1/3, 's': 1/3, 'b': 1/3,
}

best_rms_y = 999
best_xy = 0

for xy in [i*0.01 for i in range(-200, 200)]:
    total_sq = 0
    count = 0
    for f in fermions:
        name, a, b, n, m_exp, T3, Q, C2, gen = f
        if name == 'e': continue
        Y = Y_vals[name]
        n_corr = n + xy*Y
        m_pred = m_e * PHI**n_corr
        err = (m_pred/m_exp - 1)*100
        total_sq += err**2
        count += 1
    rms = math.sqrt(total_sq/count)
    if rms < best_rms_y:
        best_rms_y = rms
        best_xy = xy

print(f"Optimal: X_Y = {best_xy:.3f}, RMS = {best_rms_y:.1f}%")
print(f"(Y_L: leptons = -1, quarks = 1/3)")


# ================================================================
# COMPREHENSIVE: All single-quantum-number corrections
# ================================================================
print("\n" + "="*72)
print("SUMMARY: Best RMS with single-parameter corrections")
print("="*72)

print(f"""
  Bare formula (no correction):           RMS = 12.6%
  Old sector corrections (6 params):      RMS = 4.4%

  Single-parameter derived corrections:
    n -> n + X*2*T3 (quarks only):        RMS = {best_rms:.1f}%  (X = {best_x:.4f})
    n -> n + X*2*T3 + Y*C2:              RMS = {best_rms_2:.1f}%  (X={best_x_2:.3f}, Y={best_y_2:.3f})
    n -> (5+da)*a + (6+db)*b:            RMS = {best_rms_3:.1f}%  (da={best_da:+.4f}, db={best_db:+.4f})
    n -> n + X*Q_em:                      RMS = {best_rms_q:.1f}%  (X = {best_xq:.3f})
    n -> n + X*Y_L:                       RMS = {best_rms_y:.1f}%  (X = {best_xy:.3f})
""")

# ================================================================
# Test: The STRUCTURAL argument - integer exponent precision limit
# ================================================================
print("="*72)
print("STRUCTURAL: What is the BEST possible RMS with integer exponents?")
print("="*72)
print("\nFor any formula m_f = m_e * phi^(integer), the error comes from")
print("rounding n_eff to the nearest integer. The residual is:")

total_sq = 0
count = 0
for f in fermions:
    name, a, b, n, m_exp, T3, Q, C2, gen = f
    if name == 'e': continue
    n_eff = math.log(m_exp/m_e)/math.log(PHI)
    n_int = round(n_eff)
    m_best = m_e * PHI**n_int
    err = (m_best/m_exp - 1)*100
    delta = n_eff - n_int
    total_sq += err**2
    count += 1
    print(f"  {name:>3}: n_eff = {n_eff:.4f}, n_best = {n_int}, delta = {delta:+.4f}, err = {err:+.1f}%")

rms_best_int = math.sqrt(total_sq/count)
print(f"\nBest possible RMS with ANY integer exponents: {rms_best_int:.1f}%")
print(f"Current formula RMS: 12.6%")
print(f"\nGap = {12.6 - rms_best_int:.1f}% = room for improvement from WRONG integers")
print("(Some fermions may have the wrong integer exponent in the current formula)")

# Check which fermions have the wrong integer
print("\n  Fermions with possibly WRONG integer exponent:")
for f in fermions:
    name, a, b, n, m_exp, T3, Q, C2, gen = f
    if name == 'e': continue
    n_eff = math.log(m_exp/m_e)/math.log(PHI)
    n_best = round(n_eff)
    if n_best != n:
        m_better = m_e * PHI**n_best
        err_old = (m_e*PHI**n/m_exp - 1)*100
        err_new = (m_better/m_exp - 1)*100
        print(f"  {name}: current n={n} (err={err_old:+.1f}%), better n={n_best} (err={err_new:+.1f}%)")
        print(f"         Would need (a,b) giving 5a+6b = {n_best}")
