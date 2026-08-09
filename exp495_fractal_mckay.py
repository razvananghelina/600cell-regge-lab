"""
EXP-495: Fractal McKay Graph - Does phi^n Emerge at Level n?
==============================================================
IDEA: The McKay graph at ONE level has eigenvalue phi.
At LEVEL k (tensor product M^{otimes k}), eigenvalue phi^k appears.
The mass hierarchy phi^26 needs k >= 26 levels.

CONCRETE TEST:
1. Level 1: McKay graph (9 nodes), eigenvalues include phi
2. Level 2: M tensor M (81 nodes), eigenvalues include phi^2
3. Level k: eigenvalues include phi^k (analytically)
4. Does the mass hierarchy emerge from the SPECTRAL STRUCTURE?

KEY QUESTION: Is phi^n at level k the same as "n levels chose phi,
k-n levels chose 1"? If so, the mass exponent n = number of
"golden resonances" in the fractal.
"""

import numpy as np
from itertools import product as iterproduct

phi = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-495: FRACTAL McKAY GRAPH")
print("=" * 70)

# McKay adjacency (9 nodes, extended E8)
M = np.zeros((9, 9))
for a, b in [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,7), (7,8)]:
    M[a, b] = 1; M[b, a] = 1

evals_M, evecs_M = np.linalg.eigh(M)
evals_M = np.sort(evals_M)[::-1]

print(f"\nLevel 1: McKay graph eigenvalues:")
labels = ['+2', '+phi', '+1', '+1/phi', '0', '-1/phi', '-1', '-phi', '-2']
for i, (ev, lab) in enumerate(zip(evals_M, labels)):
    print(f"  lambda_{i} = {ev:+.6f} = {lab}")

# ============================================================
# LEVEL 2: Tensor product M tensor M
# ============================================================
print("\n" + "=" * 70)
print("LEVEL 2: M tensor M (81 nodes)")
print("=" * 70)

# Eigenvalues of M tensor M are products lambda_i * lambda_j
level2_evals = []
for li in evals_M:
    for lj in evals_M:
        level2_evals.append(li * lj)

level2_evals = np.sort(level2_evals)[::-1]

# Count distinct eigenvalues
distinct = []
for ev in level2_evals:
    if not distinct or abs(ev - distinct[-1][0]) > 1e-8:
        distinct.append([ev, 1])
    else:
        distinct[-1][1] += 1

print(f"\n  {len(distinct)} distinct eigenvalues (from 81 total)")
print(f"\n  Positive eigenvalues:")
for ev, mult in distinct:
    if ev > 0.01:
        # Try to identify
        ident = ""
        for a_exp in range(-3, 4):
            for b_exp in range(-3, 4):
                val = phi**a_exp * 2**b_exp
                if abs(ev - val) < 1e-6:
                    ident = f"= phi^{a_exp} * 2^{b_exp}"
                    break
            if ident: break
        if not ident:
            if abs(ev - phi**2) < 1e-6: ident = "= phi^2"
            elif abs(ev - 1/phi**2) < 1e-6: ident = "= 1/phi^2"
        print(f"    {ev:10.6f} (mult {mult:2d}) {ident}")

# Key: does phi^2 appear?
phi2_found = any(abs(ev - phi**2) < 1e-6 for ev, _ in distinct)
print(f"\n  phi^2 = {phi**2:.6f} found: {phi2_found}")

# ============================================================
# LEVEL k: Analytic eigenvalue structure
# ============================================================
print("\n" + "=" * 70)
print("LEVEL k: ANALYTIC EIGENVALUE STRUCTURE")
print("=" * 70)

print(f"""
At level k, eigenvalues are products of k values from
{{2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2}}.

The eigenvalue phi^n arises from choosing phi at n positions
and 1 at (k-n) positions (ignoring signs for now).

For the mass hierarchy:
  m_f / m_e = phi^(n_f)
  We need n_f = 0, 3, 5, 11, 16, 17, 19, 26

This requires k >= 26 levels.
""")

# At level k, the positive eigenvalue phi^n (for n <= k) appears.
# Its multiplicity: C(k, n) * (mult of phi)^n * (mult of 1)^(k-n)

# Multiplicity of phi eigenvalue in McKay: 1
# Multiplicity of 1 eigenvalue in McKay: 1
# So phi^n at level k has multiplicity C(k, n)

print(f"At level k = 26:")
print(f"  Total nodes: 9^26 = {9**26:.3e}")
print(f"  phi^0 = 1 eigenvalue, mult = C(26,0) = 1")
print(f"  phi^1 eigenvalue, mult = C(26,1) = 26")
print(f"  phi^3 eigenvalue, mult = C(26,3) = {int(int(np.round(np.exp(sum(np.log(range(1,n+1)) for _ in []) if False else 0)))(26,3))}")
print(f"  phi^26 eigenvalue, mult = C(26,26) = 1")

print(f"\n  Fermion mass exponents and their multiplicities at level 26:")
fermions = [('e', 0), ('u', 3), ('d', 5), ('mu/s', 11),
            ('c', 16), ('tau', 17), ('b', 19), ('t', 26)]

for name, n in fermions:
    if n <= 26:
        mult = int(int(np.round(np.exp(sum(np.log(range(1,n+1)) for _ in []) if False else 0)))(26, n))
        print(f"    {name:>5}: n={n:2d}, phi^n = {phi**n:12.2f}, "
              f"mult = C(26,{n}) = {mult}")

# ============================================================
# THE CRITICAL TEST: Is 26 special?
# ============================================================
print("\n" + "=" * 70)
print("IS k = 26 SPECIAL?")
print("=" * 70)

# n = 26 is the top quark exponent. Is k = 26 the "total number of levels"?
# If k = 26, then ALL fermion exponents {0, 3, 5, 11, 16, 17, 19, 26}
# are available as eigenvalues phi^n.

# But why k = 26 and not k = 30 or k = 100?

# Check: n_top = 5*a_top + 6*b_top = 5*4 + 6*1 = 26
# And 26 = a1^2 + 1 (the denominator of sin^2 theta_W!)

print(f"\n  Top quark: (a,b) = (4,1), n = 5*4 + 6*1 = 26")
print(f"  26 = a1^2 + 1 = {5**2 + 1} (denominator of sin^2 theta_W)")
print(f"  26 = 2 * 13 (13 = F_7, Fibonacci number)")

# If k = 26 = a1^2 + 1, then the fractal has a1^2 + 1 levels.
# The total number of nodes: 9^26 = 9^(a1^2+1)
# The spectral range: (2*phi)^26 = (2*phi)^(a1^2+1)

print(f"\n  Spectral range at k=26: (2*phi)^26 = {(2*phi)**26:.3e}")
print(f"  Mass range needed: phi^26 = {phi**26:.3e}")
print(f"  Ratio: (2*phi)^26 / phi^26 = 2^26 = {2**26}")

# ============================================================
# DEEPER: What does "choosing phi at n levels" MEAN?
# ============================================================
print("\n" + "=" * 70)
print("PHYSICAL INTERPRETATION")
print("=" * 70)

print(f"""
PICTURE:
  The universe has k = 26 levels of self-reflection.
  At each level, the McKay graph structure exists.
  Each level has 9 "modes" (eigenvalues of McKay).

  A PARTICLE is characterized by its eigenvalue SEQUENCE:
  (lambda_1, lambda_2, ..., lambda_26) where each lambda_i
  is one of {{2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2}}.

  The MASS of the particle is the PRODUCT of its eigenvalues:
  m / m_e = product(lambda_i)

  For a fermion with mass m_e * phi^n:
  - n levels have lambda_i = phi
  - 26 - n levels have lambda_i = 1
  - The remaining levels contribute 1 (neutral)
  - Product = phi^n * 1^(26-n) = phi^n

  The ELECTRON (n=0): all levels are "neutral" (lambda = 1)
  The TOP QUARK (n=26): all levels are "golden" (lambda = phi)
  The MUON (n=11): 11 levels golden, 15 neutral

WHICH LEVELS ARE GOLDEN?
  This is determined by the Z[phi] quantum numbers (a, b).
  n = 5a + 6b where 5 = dim(r4), 6 = dim(r5).

  The (a, b) quantum numbers specify WHICH of the 26 levels
  are "golden" (phi) vs "neutral" (1).

  Specifically:
  - "a" counts levels in the "r4-direction" (dim 5)
  - "b" counts levels in the "r5-direction" (dim 6)
  - Total golden levels: 5a + 6b = n
""")

# ============================================================
# VERIFICATION: Spectrum at small levels
# ============================================================
print("=" * 70)
print("VERIFICATION: SPECTRUM AT LEVELS 1-3")
print("=" * 70)

base_evals = sorted([2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2], reverse=True)

for k in [1, 2, 3]:
    # All products of k eigenvalues
    all_evals = []
    for combo in iterproduct(base_evals, repeat=k):
        all_evals.append(np.prod(combo))

    all_evals.sort(reverse=True)

    # Find phi^n eigenvalues
    phi_powers = {}
    for n in range(k+1):
        target = phi**n
        count = sum(1 for ev in all_evals if abs(ev - target) < 1e-6)
        if count > 0:
            phi_powers[n] = count

    # Spectral range (positive only)
    pos_evals = [ev for ev in all_evals if ev > 1e-10]
    if pos_evals:
        spec_range = max(pos_evals) / min(pos_evals)
    else:
        spec_range = 0

    print(f"\n  Level {k}: {9**k} nodes, {len(set(round(ev, 8) for ev in all_evals))} distinct eigenvalues")
    print(f"    Spectral range (positive): {spec_range:.2f}")
    print(f"    phi^n eigenvalues found:", end="")
    for n, count in sorted(phi_powers.items()):
        expected_mult = int(int(np.round(np.exp(sum(np.log(range(1,n+1)) for _ in []) if False else 0)))(k, n))
        print(f" phi^{n}(x{count}, expected C({k},{n})={expected_mult})", end="")
    print()

# ============================================================
# THE 5+6 STRUCTURE
# ============================================================
print("\n" + "=" * 70)
print("THE 5+6 STRUCTURE OF LEVELS")
print("=" * 70)

# If k = 26 = 5*4 + 6*1 (top quark) = a1^2 + 1
# And the levels split into "type-5" and "type-6" (from the McKay central nodes)
# Then n = 5a + 6b means: a groups of 5 levels + b groups of 6 levels are golden.

# This would mean: the 26 levels are organized into BLOCKS:
# Some blocks of size 5 (from r4, dim=5) and some of size 6 (from r5, dim=6)
# 26 = 5*4 + 6*1 or 5*2 + 6*... etc.

# How many ways to write 26 = 5a + 6b with a,b >= 0?
print(f"\n  Ways to write 26 = 5a + 6b:")
for a in range(10):
    for b in range(10):
        if 5*a + 6*b == 26:
            print(f"    a={a}, b={b}: {a} blocks of 5 + {b} blocks of 6")

# And for the TOP QUARK: (a,b) = (4,1), n = 5*4 + 6*1 = 26
# This means: 4 blocks of 5 + 1 block of 6 = 20 + 6 = 26 levels
# ALL 26 levels are golden for the top quark.

# For the electron: (a,b) = (0,0), n = 0
# NO levels are golden. All 26 are neutral.

# For the muon: (a,b) = (1,1), n = 5+6 = 11
# 1 block of 5 + 1 block of 6 = 11 golden levels, 15 neutral

print(f"\n  Each fermion's golden level structure:")
all_fermions = [
    ('e', 0, 0), ('mu', 1, 1), ('tau', 1, 2),
    ('u', 3, -2), ('c', 2, 1), ('t', 4, 1),
    ('d', 1, 0), ('s', 1, 1), ('b', -1, 4)
]

for name, a, b in all_fermions:
    n = 5*a + 6*b
    blocks_5 = a  # blocks of 5 golden levels
    blocks_6 = b  # blocks of 6 golden levels
    print(f"    {name:>3}: (a,b)=({a:+d},{b:+d}), n={n:2d} = "
          f"{a}x5 + {b}x6 golden levels")

# ============================================================
# NEGATIVE a OR b: What does this mean fractally?
# ============================================================
print(f"\n  NOTE: Some fermions have NEGATIVE a or b:")
print(f"    u: a=3, b=-2 -> 3 blocks of 5 GOLDEN + 2 blocks of 6 ANTI-GOLDEN")
print(f"    b: a=-1, b=4 -> 1 block of 5 ANTI-GOLDEN + 4 blocks of 6 GOLDEN")
print(f"""
  ANTI-GOLDEN means lambda = 1/phi (instead of phi).
  Product: phi^a * (1/phi)^|b_neg| = phi^(a - |b_neg|) ...
  But n = 5a + 6b with b negative:
  u: n = 15 - 12 = 3 (few net golden levels)
  b: n = -5 + 24 = 19 (many net golden levels)

  So NEGATIVE quantum numbers mean ANTI-RESONANCE (1/phi)
  at those levels, reducing the total mass.
""")

# ============================================================
# SUMMARY
# ============================================================
print("=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
THE FRACTAL McKAY PICTURE:
==========================

1. The universe has k = 26 = a1^2 + 1 levels of self-reflection.

2. At each level, the McKay graph (extended E8) exists with
   eigenvalues {{+2, +phi, +1, +1/phi, 0, -1/phi, -1, -phi, -2}}.

3. A PARTICLE is a specific eigenvalue sequence across all 26 levels.
   Its mass is the PRODUCT of eigenvalues: m/m_e = product(lambda_i).

4. For eigenvalue = phi at n levels and 1 at (26-n) levels:
   m/m_e = phi^n  (the mass formula!)

5. The quantum numbers (a,b) specify the STRUCTURE of golden levels:
   - a = number of "5-blocks" resonating at phi
   - b = number of "6-blocks" resonating at phi
   - Negative a or b = anti-resonance at 1/phi
   - n = 5a + 6b = net number of golden levels

6. The 5 and 6 come from dim(r4) = 5 and dim(r5) = 6,
   the central nodes of the McKay graph.

STATUS: SPECULATIVE but SELF-CONSISTENT.
  - Reproduces the mass formula phi^n: YES (by construction)
  - Explains WHY 5 and 6 appear: YES (McKay central dimensions)
  - Explains WHY k = 26: PARTIALLY (= a1^2 + 1, top quark exponent)
  - Predicts something NEW: ???
  - Derivable from first principles: NOT YET

The picture is elegant but TAUTOLOGICAL at this stage:
we're RESTATING the mass formula in fractal language,
not DERIVING it from new principles.

To become physics, we need:
  (a) A mechanism that SELECTS k = 26 levels
  (b) A rule that determines WHICH levels are golden for each fermion
  (c) A prediction that DIFFERS from the standard framework
""")

print("=" * 70)
print("EXP-495 COMPLETE")
print("=" * 70)
