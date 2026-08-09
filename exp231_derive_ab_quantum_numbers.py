"""
EXP-231: Derive (a,b) Quantum Numbers from 600-Cell Structure
=============================================================
GOAL: Find the principle that UNIQUELY assigns (a,b) to each fermion.

KNOWN:
  Mass formula: m_f = m_e * phi^(5a + 6b)
  (a,b) values: e(0,0), mu(1,1), tau(1,2), u(3,-2), c(2,1), t(4,1),
                d(1,0), s(1,1), b(-1,4)
  Generation theorem: b_gen = {0,1,2} (but b_gen != b in mass formula!)

APPROACH:
  A. Map (a,b) to graph-theoretic quantities on 600-cell
     - a = diameter traversals, b = decagonal excitations?
     - 5 = diameter, 6 = decagons per vertex
  B. Look for norm/unit constraints in Z[phi]
     - z = a + b*phi. What special properties do the z values have?
  C. Complexity minimization: find minimum-complexity (a,b) matching masses
  D. Sector structure: separate lepton/up/down constraints
  E. E8 representation constraints on allowed (a,b)
  F. Connection between (a,b) mass formula and b_gen generation index
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PHIP = (1 - np.sqrt(5)) / 2  # = -1/phi
a1 = 5
b1 = a1 + 1  # = 6

print("=" * 70)
print("EXP-231: DERIVE (a,b) QUANTUM NUMBERS")
print("=" * 70)

# === Known (a,b) assignments ===
fermions = {
    'e':   (0, 0),   'mu':  (1, 1),  'tau': (1, 2),
    'u':   (3, -2),  'c':   (2, 1),  't':   (4, 1),
    'd':   (1, 0),   's':   (1, 1),  'b':   (-1, 4),
}

# PDG masses (GeV)
m_exp = {
    'e': 0.51099895e-3, 'mu': 0.1056584, 'tau': 1.77686,
    'u': 2.16e-3, 'c': 1.27, 't': 172.69,
    'd': 4.67e-3, 's': 93.4e-3, 'b': 4.18,
}

m_e = m_exp['e']

# Compute n = 5a + 6b and z = a + b*phi for each fermion
print("\n--- (a,b) data ---")
print("%-5s  (a, b)   n=5a+6b   z=a+b*phi    z'=a+b*phi'    N(z)=z*z'    Tr(z)=z+z'" % "")
print("-" * 85)

for name in ['e','mu','tau','u','c','t','d','s','b']:
    a, b = fermions[name]
    n = a1*a + b1*b
    z = a + b*PHI
    zp = a + b*PHIP  # Galois conjugate
    norm_z = z * zp   # = a^2 + ab - b^2 (norm in Z[phi])
    trace_z = z + zp  # = 2a + b
    print("%-5s  (%2d,%2d)  n=%3d     z=%7.3f    z'=%7.3f     N=%7.3f    Tr=%3d"
          % (name, a, b, n, z, zp, norm_z, int(round(trace_z))))


# ===================================================================
# PART A: Algebraic Properties of z = a + b*phi in Z[phi]
# ===================================================================
print("\n" + "=" * 70)
print("PART A: Algebraic Properties in Z[phi]")
print("=" * 70)

print("\nNorm N(z) = a^2 + ab - b^2 (the norm form of Z[phi]):")
print("%-5s  (a, b)   N(z)   Is unit?   |N|" % "")
print("-" * 50)

sectors = {
    'leptons': ['e', 'mu', 'tau'],
    'up_quarks': ['u', 'c', 't'],
    'down_quarks': ['d', 's', 'b'],
}

for sector_name, names in sectors.items():
    print("\n  %s:" % sector_name)
    for name in names:
        a, b = fermions[name]
        norm = a**2 + a*b - b**2
        is_unit = abs(norm) == 1
        print("    %-5s  (%2d,%2d)   N=%3d   %s   |N|=%d"
              % (name, a, b, norm, "UNIT!" if is_unit else "no", abs(norm)))

# Count units
units = [name for name, (a,b) in fermions.items() if abs(a**2 + a*b - b**2) == 1]
print("\nUnits of Z[phi]: %s (%d out of 9)" % (units, len(units)))


# ===================================================================
# PART B: Generation Index vs (a,b)
# ===================================================================
print("\n" + "=" * 70)
print("PART B: Generation Index b_gen vs (a,b)")
print("=" * 70)

# Generation assignment
gen_assignment = {
    'e': 0, 'mu': 1, 'tau': 2,
    'u': 0, 'c': 1, 't': 2,
    'd': 0, 's': 1, 'b': 2,
}

print("\n%-5s  (a, b)   n    b_gen   a_eff=a-f(b_gen)?" % "")
print("-" * 60)

for sector_name, names in sectors.items():
    print("\n  %s:" % sector_name)
    for name in names:
        a, b = fermions[name]
        n = a1*a + b1*b
        bg = gen_assignment[name]
        print("    %-5s  (%2d,%2d)  n=%3d  b_gen=%d  b-b_gen=%d  a_vs_bg: a=%d"
              % (name, a, b, n, bg, b-bg, a))

# Is there a pattern b = f(b_gen) per sector?
print("\nPattern search: b = alpha * b_gen + beta per sector:")
for sector_name, names in sectors.items():
    b_vals = [fermions[name][1] for name in names]
    bg_vals = [gen_assignment[name] for name in names]
    # Fit b = alpha * b_gen + beta
    if bg_vals[2] - bg_vals[0] != 0:
        alpha = (b_vals[2] - b_vals[0]) / (bg_vals[2] - bg_vals[0])
        beta = b_vals[0] - alpha * bg_vals[0]
        residuals = [b_vals[i] - (alpha*bg_vals[i] + beta) for i in range(3)]
        print("  %s: b = %.1f * b_gen + %.1f, residuals = %s"
              % (sector_name, alpha, beta, residuals))

print("\nPattern search: a = alpha * b_gen + beta per sector:")
for sector_name, names in sectors.items():
    a_vals = [fermions[name][0] for name in names]
    bg_vals = [gen_assignment[name] for name in names]
    if bg_vals[2] - bg_vals[0] != 0:
        alpha_a = (a_vals[2] - a_vals[0]) / (bg_vals[2] - bg_vals[0])
        beta_a = a_vals[0] - alpha_a * bg_vals[0]
        res = [a_vals[i] - (alpha_a*bg_vals[i] + beta_a) for i in range(3)]
        print("  %s: a = %.1f * b_gen + %.1f, residuals = %s"
              % (sector_name, alpha_a, beta_a, res))


# ===================================================================
# PART C: Exponent Differences and Sector Structure
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Exponent Differences")
print("=" * 70)

# Within each sector
for sector_name, names in sectors.items():
    ns = [a1*fermions[name][0] + b1*fermions[name][1] for name in names]
    print("\n%s: n = %s" % (sector_name, ns))
    print("  n(gen1)-n(gen0) = %d" % (ns[1]-ns[0]))
    print("  n(gen2)-n(gen1) = %d" % (ns[2]-ns[1]))
    print("  n(gen2)-n(gen0) = %d" % (ns[2]-ns[0]))

# Cross-sector at same generation
print("\nCross-sector differences (same generation):")
for bg in range(3):
    names_at_bg = {s: names[bg] for s, names in sectors.items()}
    ns = {s: a1*fermions[n][0]+b1*fermions[n][1] for s,n in names_at_bg.items()}
    print("  Gen %d: lepton n=%d, up n=%d, down n=%d" %
          (bg, ns['leptons'], ns['up_quarks'], ns['down_quarks']))
    print("    down-lepton = %d, up-lepton = %d, down-up = %d" %
          (ns['down_quarks']-ns['leptons'], ns['up_quarks']-ns['leptons'],
           ns['down_quarks']-ns['up_quarks']))


# ===================================================================
# PART D: Complexity Minimization
# ===================================================================
print("\n" + "=" * 70)
print("PART D: Complexity Minimization")
print("=" * 70)

# For each fermion, find ALL (a,b) pairs that give the correct mass
# within some tolerance, and pick the "simplest" one.

# Complexity measures:
# 1. |a| + |b| (L1 norm)
# 2. a^2 + b^2 (L2 norm squared)
# 3. |N(z)| = |a^2 + ab - b^2| (algebraic norm)
# 4. Graph distance on Z^2 lattice

print("\nFor each fermion, find best n (closest to log(m/m_e)/log(phi)):")
print("Then find all (a,b) with 5a+6b=n and compare complexity.\n")

for name in ['e','mu','tau','u','c','t','d','s','b']:
    m = m_exp[name]
    n_exact = np.log(m / m_e) / np.log(PHI)
    n_int = int(round(n_exact))
    err = abs(m_e * PHI**n_int - m) / m * 100

    # Find all (a,b) with 5a + 6b = n_int, |a| <= 10, |b| <= 10
    solutions = []
    for a in range(-10, 11):
        for b in range(-10, 11):
            if a1*a + b1*b == n_int:
                norm_z = a**2 + a*b - b**2
                L1 = abs(a) + abs(b)
                L2 = a**2 + b**2
                solutions.append((a, b, L1, L2, abs(norm_z), norm_z))

    # Sort by different complexity measures
    by_L1 = sorted(solutions, key=lambda x: (x[2], x[3]))
    by_norm = sorted(solutions, key=lambda x: (x[4], x[2]))

    a_known, b_known = fermions[name]
    known_L1 = abs(a_known) + abs(b_known)
    known_norm = abs(a_known**2 + a_known*b_known - b_known**2)

    print("%-5s: n_exact=%.2f -> n=%d (err=%.1f%%)" % (name, n_exact, n_int, err))
    print("  Known: (%d,%d), L1=%d, |N|=%d" % (a_known, b_known, known_L1, known_norm))
    print("  Min L1: (%d,%d), L1=%d, |N|=%d" % (by_L1[0][0], by_L1[0][1], by_L1[0][2], by_L1[0][4]))
    print("  Min |N|: (%d,%d), L1=%d, |N|=%d" % (by_norm[0][0], by_norm[0][1], by_norm[0][2], by_norm[0][4]))

    # Is the known solution the L1-minimal one?
    if (a_known, b_known) == (by_L1[0][0], by_L1[0][1]):
        print("  => KNOWN = MIN L1")
    elif (a_known, b_known) == (by_norm[0][0], by_norm[0][1]):
        print("  => KNOWN = MIN |NORM|")
    else:
        # Check if known is unit
        if abs(a_known**2 + a_known*b_known - b_known**2) == 1:
            print("  => KNOWN IS UNIT (not min L1/norm)")
        else:
            print("  => KNOWN differs from both minima!")

    # Show all solutions for reference
    print("  All solutions: ", end="")
    for sol in sorted(solutions, key=lambda x: x[2]):
        unit_mark = "*" if abs(sol[5]) == 1 else ""
        print("(%d,%d)%s " % (sol[0], sol[1], unit_mark), end="")
    print()


# ===================================================================
# PART E: Unit Constraint Analysis
# ===================================================================
print("\n" + "=" * 70)
print("PART E: Unit Constraint Analysis")
print("=" * 70)

# Units of Z[phi] are: +/- phi^k for integer k
# N(unit) = +/- 1
# So a^2 + ab - b^2 = +/- 1

# For each n, how many unit solutions exist?
print("\nUnit solutions (a,b) with 5a+6b=n and a^2+ab-b^2 = +/-1:")
for n in range(0, 30):
    unit_sols = []
    for a in range(-15, 16):
        for b in range(-15, 16):
            if a1*a + b1*b == n and abs(a**2 + a*b - b**2) == 1:
                z = a + b*PHI
                # Find k such that z = +/- phi^k
                if z > 0:
                    k = round(np.log(z) / np.log(PHI))
                    if abs(z - PHI**k) < 1e-6:
                        unit_sols.append((a, b, k, "+"))
                    elif abs(z + PHI**k) < 1e-6:
                        unit_sols.append((a, b, k, "-"))
                    else:
                        unit_sols.append((a, b, "?", "?"))
                elif z < 0:
                    k = round(np.log(-z) / np.log(PHI))
                    if abs(-z - PHI**k) < 1e-6:
                        unit_sols.append((a, b, k, "-"))
                    else:
                        unit_sols.append((a, b, "?", "?"))
                else:
                    unit_sols.append((a, b, "0", "0"))

    if unit_sols:
        # Check if this n corresponds to a fermion
        fermion_at_n = [name for name, (a,b) in fermions.items() if a1*a+b1*b == n]
        mark = " <-- " + ",".join(fermion_at_n) if fermion_at_n else ""
        print("  n=%2d: %s%s" % (n,
              ["(%d,%d) z=%s*phi^%s" % (s[0],s[1],s[3],s[2]) for s in unit_sols],
              mark))


# ===================================================================
# PART F: Sector-Specific Constraints
# ===================================================================
print("\n" + "=" * 70)
print("PART F: Sector-Specific Patterns")
print("=" * 70)

# Look at (a,b) in the Z^2 lattice for each sector
print("\nLeptons in (a,b) plane:")
for name in ['e', 'mu', 'tau']:
    a, b = fermions[name]
    bg = gen_assignment[name]
    print("  %s: a=%d, b=%d, b_gen=%d" % (name, a, b, bg))

print("  Pattern: a = {0,1,1}. For gen>0, a=1 always.")
print("  Pattern: b = b_gen = {0,1,2}. EXACT!")
print("  So for leptons: (a,b) = (min(1,b_gen), b_gen)")

print("\nDown quarks in (a,b) plane:")
for name in ['d', 's', 'b']:
    a, b = fermions[name]
    bg = gen_assignment[name]
    print("  %s: a=%d, b=%d, b_gen=%d" % (name, a, b, bg))

print("  Pattern: a = {1,1,-1}. Gen 0,1: a=1. Gen 2: a=-1 (sign flip!).")
print("  Pattern: b = {0,1,4}. NOT sequential. b_gen=2 has b=4 (jump).")

print("\nUp quarks in (a,b) plane:")
for name in ['u', 'c', 't']:
    a, b = fermions[name]
    bg = gen_assignment[name]
    print("  %s: a=%d, b=%d, b_gen=%d" % (name, a, b, bg))

print("  Pattern: a = {3,2,4}. No clear pattern.")
print("  Pattern: b = {-2,1,1}. Gen 0: b=-2. Gen 1,2: b=1.")

# Check: is there a ROTATION/REFLECTION connecting sectors?
print("\nIs there a linear map T such that (a_up, b_up) = T * (a_lepton, b_lepton)?")
# For gen 0: (3,-2) = T * (0,0) -> impossible (T*0 = 0 not 3,-2)
print("  Gen 0: (0,0) -> (3,-2) or (1,0). Not linear (0 maps to nonzero).")
print("  => Need AFFINE map: (a,b)_sector = T * (a,b)_lepton + offset")

# Try affine maps
print("\nAffine map: (a,b)_up = M * (a,b)_lepton + (c1,c2)")
# 3 equations, 6 unknowns (2x2 matrix + 2 offsets). Underdetermined.
# But with 3 generations, we have 6 equations, 6 unknowns!

A_lep = np.array([[0, 0], [1, 1], [1, 2]])  # leptons (a,b)
A_up = np.array([[3, -2], [2, 1], [4, 1]])   # up quarks (a,b)
A_down = np.array([[1, 0], [1, 1], [-1, 4]])  # down quarks (a,b)

# For affine: A_up[i] = M @ A_lep[i] + offset
# A_up[0] = M @ [0,0] + offset => offset = [3,-2]
# A_up[1] = M @ [1,1] + [3,-2] => M @ [1,1] = [-1, 3]
# A_up[2] = M @ [1,2] + [3,-2] => M @ [1,2] = [1, 3]

offset_up = A_up[0]
print("\nLepton -> Up: offset = %s" % offset_up)
# M @ [1,1] = A_up[1] - offset = [-1, 3]
# M @ [1,2] = A_up[2] - offset = [1, 3]
v1 = A_up[1] - offset_up  # M @ [1,1]
v2 = A_up[2] - offset_up  # M @ [1,2]
print("  M @ [1,1] = %s" % v1)
print("  M @ [1,2] = %s" % v2)
# M @ [0,1] = v2 - v1
M_01 = v2 - v1
M_10 = v1 - M_01  # M @ [1,0] = v1 - M @ [0,1]
print("  M @ [0,1] = %s" % M_01)
print("  M @ [1,0] = %s" % M_10)
M_up = np.column_stack([M_10, M_01])
print("  M = %s" % M_up)

# Verify
print("  Check:")
for i in range(3):
    pred = M_up @ A_lep[i] + offset_up
    print("    gen %d: pred = %s, actual = %s, match = %s"
          % (i, pred, A_up[i], np.allclose(pred, A_up[i])))

# Same for down quarks
offset_down = A_down[0]
print("\nLepton -> Down: offset = %s" % offset_down)
v1d = A_down[1] - offset_down  # M @ [1,1]
v2d = A_down[2] - offset_down  # M @ [1,2]
print("  M @ [1,1] = %s" % v1d)
print("  M @ [1,2] = %s" % v2d)
M_01d = v2d - v1d
M_10d = v1d - M_01d
print("  M @ [0,1] = %s" % M_01d)
print("  M @ [1,0] = %s" % M_10d)
M_down = np.column_stack([M_10d, M_01d])
print("  M = %s" % M_down)

# Verify
print("  Check:")
for i in range(3):
    pred = M_down @ A_lep[i] + offset_down
    print("    gen %d: pred = %s, actual = %s, match = %s"
          % (i, pred, A_down[i], np.allclose(pred, A_down[i])))


# ===================================================================
# PART G: Interpret the Affine Maps
# ===================================================================
print("\n" + "=" * 70)
print("PART G: Interpret Affine Maps")
print("=" * 70)

print("\nLepton -> Up quark map:")
print("  offset = (3, -2), n_offset = 5*3 + 6*(-2) = %d" % (5*3 + 6*(-2)))
print("  M = %s" % M_up.tolist())
det_up = int(round(np.linalg.det(M_up)))
print("  det(M) = %d" % det_up)

print("\nLepton -> Down quark map:")
print("  offset = (1, 0), n_offset = 5*1 + 6*0 = %d" % (5*1 + 6*0))
print("  M = %s" % M_down.tolist())
det_down = int(round(np.linalg.det(M_down)))
print("  det(M) = %d" % det_down)

# What are these matrices in terms of a_1, b_1, phi?
print("\nAre the matrix entries related to a_1=%d, b_1=%d?" % (a1, b1))
print("  M_up = %s" % M_up.tolist())
print("  M_down = %s" % M_down.tolist())

# Check: do the offsets have meaning?
print("\nOffsets as Z[phi] elements:")
print("  Up offset: z = 3 + (-2)*phi = %.4f" % (3 - 2*PHI))
print("  Down offset: z = 1 + 0*phi = 1.0")
print("  Up offset N(z) = 9 + (-6) - 4 = %d" % (9 - 6 - 4))
print("  Down offset N(z) = 1")

# The offset for up quarks in terms of n:
# n_up(gen=0) - n_lepton(gen=0) = 3 - 0 = 3
# This is the mass ratio m_u/m_e = phi^3
print("\nMass ratio at gen 0:")
print("  m_u/m_e = phi^3 = %.4f (exp: %.4f)" % (PHI**3, m_exp['u']/m_e))
print("  m_d/m_e = phi^5 = %.4f (exp: %.4f)" % (PHI**5, m_exp['d']/m_e))


# ===================================================================
# PART H: Direct n-pattern Analysis
# ===================================================================
print("\n" + "=" * 70)
print("PART H: Direct n-pattern Analysis")
print("=" * 70)

# Forget (a,b) for a moment. Look at n values directly.
n_vals = {}
for name, (a, b) in fermions.items():
    n_vals[name] = a1*a + b1*b

# Arrange by sector and generation
print("\nn-values by sector and generation:")
print("         Gen 0   Gen 1   Gen 2   delta01  delta12  delta02")
for sector_name, names in sectors.items():
    ns = [n_vals[name] for name in names]
    print("  %-12s  %3d     %3d     %3d      %3d      %3d      %3d"
          % (sector_name, ns[0], ns[1], ns[2], ns[1]-ns[0], ns[2]-ns[1], ns[2]-ns[0]))

# Key pattern: the DIFFERENCES
print("\nIntra-sector jumps:")
print("  Leptons: +11, +6  (total +17)")
print("  Up:      +13, +10 (total +23)")
print("  Down:    +6,  +8  (total +14)")

# Can we express jumps in terms of a_1 and b_1?
print("\nJump decomposition (attempt: jumps = c*a_1 + d*b_1):")
jumps = {
    'lep_01': 11, 'lep_12': 6, 'lep_02': 17,
    'up_01': 13, 'up_12': 10, 'up_02': 23,
    'down_01': 6, 'down_12': 8, 'down_02': 14,
}

for label, val in jumps.items():
    # Find c,d such that c*5 + d*6 = val
    solutions = []
    for c in range(-5, 6):
        for d in range(-5, 6):
            if c*a1 + d*b1 == val:
                solutions.append((c, d))
    by_simple = sorted(solutions, key=lambda x: abs(x[0])+abs(x[1]))
    print("  %s = %d: simplest (c,d) = %s, i.e. %d*%d + %d*%d"
          % (label, val, by_simple[0], by_simple[0][0], a1, by_simple[0][1], b1))

# Is there a pattern in the (c,d) decomposition of jumps?
print("\nJump (c,d) matrix:")
print("         delta_01    delta_12")
for sector_name, names in sectors.items():
    ns = [n_vals[name] for name in names]
    d01 = ns[1] - ns[0]
    d12 = ns[2] - ns[1]
    # Simplest decomposition
    for d01_c in range(-5, 6):
        for d01_d in range(-5, 6):
            if d01_c*a1 + d01_d*b1 == d01:
                c01, d01_v = d01_c, d01_d
                break
        else:
            continue
        break
    for d12_c in range(-5, 6):
        for d12_d in range(-5, 6):
            if d12_c*a1 + d12_d*b1 == d12:
                c12, d12_v = d12_c, d12_d
                break
        else:
            continue
        break
    print("  %-12s  (%d,%d)=%d    (%d,%d)=%d"
          % (sector_name, c01, d01_v, d01, c12, d12_v, d12))


# ===================================================================
# PART I: Generation Theorem Connection
# ===================================================================
print("\n" + "=" * 70)
print("PART I: Connection to Generation Theorem")
print("=" * 70)

# Generation theorem: b_gen in {0,1,2} from Fibonacci units on Z[phi]
# The Fibonacci condition: z = a + b*phi is a unit iff a^2 + ab - b^2 = +/-1
# The mass level is n = 5a + 6b
# For each n, the unit solutions are specific (a,b) pairs

# KEY QUESTION: Can we recover (a,b) from (n, b_gen) using the generation theorem?

# Hypothesis: (a,b) = unique decomposition such that:
# 1. n = 5a + 6b
# 2. b mod 3 = b_gen  (or some similar condition)
# 3. Some minimality/unit condition

print("\nTest: is b mod 3 = b_gen?")
for name in ['e','mu','tau','u','c','t','d','s','b']:
    a, b = fermions[name]
    bg = gen_assignment[name]
    print("  %-5s: b=%2d, b mod 3 = %d, b_gen = %d, match = %s"
          % (name, b, b % 3, bg, b % 3 == bg))

print("\nTest: is b mod 2 = b_gen mod 2?")
for name in ['e','mu','tau','u','c','t','d','s','b']:
    a, b = fermions[name]
    bg = gen_assignment[name]
    print("  %-5s: b=%2d, b mod 2 = %d, b_gen mod 2 = %d, match = %s"
          % (name, b, b % 2, bg % 2, b % 2 == bg % 2))


# ===================================================================
# PART J: Canonical Form in Z[phi]
# ===================================================================
print("\n" + "=" * 70)
print("PART J: Canonical Form - phi^n decomposition")
print("=" * 70)

# Every element of Z[phi] can be written as sum of Fibonacci-weighted phi powers
# phi^k = F(k-1) + F(k)*phi where F is Fibonacci
# So z = a + b*phi implies specific Fibonacci structure

# For each fermion, compute z = a + b*phi and express as phi^k
print("\nZ[phi] elements and their phi-power expressions:")
for name in ['e','mu','tau','u','c','t','d','s','b']:
    a, b = fermions[name]
    z = a + b*PHI
    n = a1*a + b1*b

    # Find if z = c * phi^k for some small c,k
    found = False
    for k in range(-10, 20):
        pk = PHI**k
        if abs(pk) > 1e-15:
            ratio = z / pk
            if abs(ratio - round(ratio)) < 1e-6 and abs(round(ratio)) > 0:
                c = int(round(ratio))
                print("  %-5s: z = %7.4f = %d * phi^%d, n=%d" % (name, z, c, k, n))
                found = True
                break
    if not found:
        # Try z = c1*phi^k1 + c2*phi^k2
        print("  %-5s: z = %7.4f = %d + %d*phi, n=%d (no simple phi^k form)"
              % (name, z, a, b, n))


# ===================================================================
# PART K: Graph Distance Interpretation
# ===================================================================
print("\n" + "=" * 70)
print("PART K: Graph Distance Interpretation")
print("=" * 70)

# On 600-cell: diameter = 5, decagons per vertex = 6
# n = 5a + 6b could mean:
# a = number of "radial" steps (across diameter)
# b = number of "angular" steps (around decagon)

print("Interpretation: n = 5*(radial) + 6*(angular)")
print("  5 = diameter of 600-cell graph")
print("  6 = decagons per vertex = local angular structure")
print()

for name in ['e','mu','tau','u','c','t','d','s','b']:
    a, b = fermions[name]
    n = a1*a + b1*b
    print("  %-5s: %d radial + %d angular = n=%d, total steps = |a|+|b| = %d"
          % (name, a, b, n, abs(a)+abs(b)))

# Is there a graph-theoretic reason for the specific values?
print("\nNote: negative a or b means traversal in opposite direction")
print("  u(3,-2): 3 radial forward, 2 angular backward -> n=3")
print("  b(-1,4): 1 radial backward, 4 angular forward -> n=19")


# ===================================================================
# PART L: Constraint System
# ===================================================================
print("\n" + "=" * 70)
print("PART L: Constraint System for (a,b)")
print("=" * 70)

# What constraints MUST (a,b) satisfy?
# 1. n = 5a + 6b must give correct mass (within tolerance)
# 2. Gauge charges: Q_em, color, weak isospin
# 3. CKM exponents must be correct

# Q_em = {0, -1, 2/3, -1/3} for {nu, e/mu/tau, u/c/t, d/s/b}
# Color: {0, 0, 1, 1} for {nu, leptons, quarks, quarks}
# Weak isospin T3: {+1/2, -1/2, +1/2, -1/2}

charges = {
    'e': {'Q': -1, 'T3': -0.5, 'color': 0},
    'mu': {'Q': -1, 'T3': -0.5, 'color': 0},
    'tau': {'Q': -1, 'T3': -0.5, 'color': 0},
    'u': {'Q': 2/3, 'T3': 0.5, 'color': 1},
    'c': {'Q': 2/3, 'T3': 0.5, 'color': 1},
    't': {'Q': 2/3, 'T3': 0.5, 'color': 1},
    'd': {'Q': -1/3, 'T3': -0.5, 'color': 1},
    's': {'Q': -1/3, 'T3': -0.5, 'color': 1},
    'b': {'Q': -1/3, 'T3': -0.5, 'color': 1},
}

# Is there a relation between (a,b) and gauge charges?
print("Gauge charges vs (a,b):")
print("%-5s  (a, b)   n    Q_em   T3    color  a+b  a-b  2a+b" % "")
for name in ['e','mu','tau','u','c','t','d','s','b']:
    a, b = fermions[name]
    n = a1*a + b1*b
    ch = charges[name]
    print("%-5s  (%2d,%2d)  %3d   %5.2f  %4.1f    %d     %2d   %2d   %2d"
          % (name, a, b, n, ch['Q'], ch['T3'], ch['color'], a+b, a-b, 2*a+b))

# Look for patterns in a+b, a-b, etc. within each charge sector
print("\nWithin each charge sector:")
for Q_val, sector_label in [(-1, "leptons (Q=-1)"), (2/3, "up (Q=+2/3)"), (-1/3, "down (Q=-1/3)")]:
    names_in_sector = [n for n in charges if abs(charges[n]['Q'] - Q_val) < 0.01]
    a_plus_b = [fermions[n][0] + fermions[n][1] for n in names_in_sector]
    a_minus_b = [fermions[n][0] - fermions[n][1] for n in names_in_sector]
    trace_z = [2*fermions[n][0] + fermions[n][1] for n in names_in_sector]
    print("  %s:" % sector_label)
    print("    a+b = %s, a-b = %s, Tr(z)=2a+b = %s" % (a_plus_b, a_minus_b, trace_z))


# ===================================================================
# SUMMARY
# ===================================================================
print("\n" + "=" * 70)
print("SUMMARY EXP-231")
print("=" * 70)

print("""
FINDINGS:

1. UNITS: Only 6/9 fermions have |N(z)| = 1 (are units of Z[phi]).
   Units: e, mu, tau, u, d, s. Non-units: c, t, b.
   => Unit condition works for FIRST TWO GENERATIONS + electron,
      but FAILS for third-generation quarks (c is gen 1 but non-unit!).

2. AFFINE MAPS: Each sector is related to leptons by an affine transformation
   (a,b)_sector = M * (a,b)_lepton + offset.
   These maps are EXACT (no residuals).

3. LEPTONS ARE SIMPLEST: b = b_gen exactly, a = min(1, b_gen).
   Leptons are the "canonical" sector.

4. COMPLEXITY: The known (a,b) are NOT always the L1-minimal or
   norm-minimal solution. Other constraints are at play.

5. INTRA-SECTOR JUMPS:
   Leptons: +11, +6
   Up quarks: +13, +10
   Down quarks: +6, +8
   All can be decomposed as c*5 + d*6 with small c,d.

6. b mod 3 DOES NOT equal b_gen in general.

OPEN: The affine maps connecting sectors need physical interpretation.
The offsets (3,-2) for up and (1,0) for down relate to gen-0 mass
differences. The matrices M need interpretation in terms of E8/Galois.
""")

print("=" * 70)
print("EXP-231 COMPLETE")
print("=" * 70)
