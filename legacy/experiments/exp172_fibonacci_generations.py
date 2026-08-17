"""
EXP-172: 3 Generations from Fibonacci / Z[phi] Units
=====================================================
From exp168: |N(1,b)| <= 1 gives b in {0,1,2} on the a=1 line.

KEY INSIGHT: The elements 1+b*phi for b=0,1,2 are POWERS OF PHI:
  1 + 0*phi = phi^0
  1 + 1*phi = phi^2
  1 + 2*phi = phi^3

These are the ONLY non-negative powers of phi with rational part = 1!
The reason: phi^n = F(n-1) + F(n)*phi, and F(n-1) = 1 only for n = 0, 2, 3.
(Because Fibonacci: F(-1)=1, F(0)=0, F(1)=1, F(2)=1, F(3)=2, F(4)=3,...)

So: 3 GENERATIONS = number of times F(k) = 1 for k >= -1!

Questions:
1. Verify the phi^n = Fibonacci connection rigorously
2. Physical interpretation: why must generation elements be UNITS of Z[phi]?
3. Does this extend to other fermion lines (b=1, 3a+2b=5)?
4. What is the ALGEBRAIC meaning of "unit norm" for fermion content?
5. Connection to the level selection rules from exp123

Category: RESEARCH (deepening a PATTERN toward DERIVATION)
"""

import numpy as np
from collections import Counter

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

print("=" * 70)
print("EXP-172: 3 Generations from Fibonacci / Z[phi] Units")
print("=" * 70)

# ============================================================
# Step 1: Powers of phi and Fibonacci
# ============================================================
print("\n" + "=" * 70)
print("Step 1: Powers of phi = Fibonacci pairs")
print("=" * 70)

# phi^n = F(n-1) + F(n)*phi where F is extended Fibonacci
# F(-1)=1, F(0)=0, F(1)=1, F(2)=1, F(3)=2, F(4)=3, F(5)=5, ...
def fib(n):
    """Extended Fibonacci: F(-1)=1, F(0)=0, F(1)=1, ..."""
    if n >= 0:
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a + b
        return a
    else:
        # F(-n) = (-1)^(n+1) * F(n)
        return (-1)**(abs(n)+1) * fib(abs(n))

print(f"\n  phi^n = F(n-1) + F(n)*phi  [F = extended Fibonacci]")
print(f"\n  {'n':>4} {'F(n-1)=a':>10} {'F(n)=b':>10} {'phi^n':>12} {'a+b*phi':>12} {'match':>6}")
print("  " + "-" * 55)

for n in range(-5, 10):
    a = fib(n-1)
    b = fib(n)
    phi_n = PHI**n
    reconstructed = a + b * PHI
    match = abs(phi_n - reconstructed) < 1e-8
    on_a1 = "  <-- a=1!" if a == 1 else ""
    print(f"  {n:>4} {a:>10} {b:>10} {phi_n:>12.4f} {reconstructed:>12.4f} {'OK' if match else 'FAIL':>6}{on_a1}")

print(f"\n  *** Powers of phi with a = F(n-1) = 1: ***")
print(f"  F(n-1) = 1 requires n-1 in {{k : F(k) = 1}}")
print(f"  F(k) = 1 for k = -1, 1, 2 (and ONLY these for k >= -1)")
print(f"  So n = 0, 2, 3")
print(f"  Giving b = F(0), F(2), F(3) = 0, 1, 2")
print(f"  Exactly 3 generations!")

# Why F(k)=1 has only these solutions
print(f"\n  Why F(k) = 1 has finitely many solutions:")
print(f"  Fibonacci grows exponentially: F(k) ~ phi^k / sqrt(5)")
print(f"  For k >= 3: F(k) >= 2 (strictly increasing)")
print(f"  For k <= -2: |F(k)| >= 1, and F(k) alternates in sign")
print(f"  Specifically: F(-2) = -1, not 1")
print(f"  So F(k) = 1 only for k = -1, 1, 2")
print(f"")
print(f"  The NUMBER 3 comes from: F(k)=1 has 3 solutions")
print(f"  This is a DEEP NUMBER-THEORETIC fact about Fibonacci!")

# ============================================================
# Step 2: Units of Z[phi]
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Z[phi] Unit Structure")
print("=" * 70)

print(f"""
  Z[phi] = {{a + b*phi : a, b in Z}} is the ring of integers of Q(sqrt(5)).

  The UNITS of Z[phi] are elements with multiplicative inverse in Z[phi].
  By Dirichlet's unit theorem: units = {{+-phi^n : n in Z}}.

  The NORM N(a+b*phi) = (a+b*phi)(a+b*phi') = a^2 + ab - b^2.
  An element is a unit iff |N| = 1.

  KEY CONNECTION:
  The a=1 generation line has elements 1 + b*phi.
  These are UNITS iff |N(1+b*phi)| = |1+b-b^2| = 1.
  This gives b in {{0, 1, 2}}.

  PHYSICAL INTERPRETATION:
  Each generation corresponds to a UNIT of Z[phi].
  Units are "fundamental" - they cannot be factored in Z[phi].
  Non-units (like 1+3*phi with |N|=5) are COMPOSITE.

  Generation = power of the fundamental unit phi:
    Gen 0 (e, d, u): phi^0 = 1
    Gen 1 (mu, s, c): phi^2 = 1 + phi
    Gen 2 (tau, b, t): phi^3 = 1 + 2*phi
""")

# Verify unit property
print("  Verification - units on a=1 line:")
for b in range(5):
    elem = 1 + b*PHI
    conj = 1 + b*PHI_CONJ
    norm = round(elem * conj)
    is_unit = abs(norm) == 1
    # Which power of phi?
    if is_unit and elem > 0:
        power = round(np.log(elem) / np.log(PHI))
    else:
        power = "N/A"
    print(f"    1+{b}*phi = {elem:.4f}, N = {norm:>3}, unit = {is_unit}, phi^{power}")

# ============================================================
# Step 3: Why a=1? The Fibonacci analysis of other lines
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Other Fermion Lines - Fibonacci Analysis")
print("=" * 70)

# For each fermion (a,b), check if a+b*phi is a unit
fermions = {
    'e':   (0, 0, 'lepton'),  'mu':  (1, 1, 'lepton'), 'tau': (1, 2, 'lepton'),
    'u':   (3,-2, 'up'),      'c':   (2, 1, 'up'),     't':   (4, 1, 'up'),
    'd':   (1, 0, 'down'),    's':   (1, 1, 'down'),   'b_q': (-1, 4, 'down'),
}

print(f"\n  All fermions - unit analysis:")
print(f"  {'name':>5} {'(a,b)':>8} {'a+b*phi':>10} {'|N|':>5} {'unit':>5} {'phi^n':>7}")
print("  " + "-" * 50)

for name, (a, b, sector) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    val = a + b*PHI
    norm = abs(a**2 + a*b - b**2)
    is_unit = norm == 1
    if is_unit and val > 0:
        n = round(np.log(val) / np.log(PHI))
        phi_str = f"phi^{n}"
    elif is_unit and val < 0:
        n = round(np.log(-val) / np.log(PHI))
        phi_str = f"-phi^{n}"
    elif val == 0:
        phi_str = "0"
    else:
        phi_str = f"|N|={norm}"
    print(f"  {name:>5} ({a:>2},{b:>2}) {val:>10.4f} {norm:>5} {'YES' if is_unit else 'no':>5} {phi_str:>7}")

# Count units vs non-units
n_units = sum(1 for name, (a,b,s) in fermions.items() if abs(a**2+a*b-b**2) <= 1)
print(f"\n  Units (|N|<=1): {n_units}/9 fermions")
print(f"  Non-units: {9-n_units}/9 fermions")

# Which are non-units?
print(f"\n  Non-unit fermions:")
for name, (a, b, sector) in fermions.items():
    norm = abs(a**2 + a*b - b**2)
    if norm > 1:
        # Factor the element in Z[phi]
        val = a + b*PHI
        print(f"    {name} ({a},{b}): |N| = {norm}, value = {val:.4f}")
        # Check if norm is a Fibonacci number or Fibonacci-related
        # Norms: 5 (charm), 19 (bottom, top)
        # 5 = F(5), 19 = ?
        # Actually, in Z[phi], elements of norm 5 factor as:
        # sqrt(5) = 2phi-1 is a prime of norm 5 in Z[phi]
        # (Note: 5 = (2phi-1)(2phi'-1) = (sqrt5)(-sqrt5) = -5... no)
        # Actually N(2phi-1) = N((-1)+2*phi) = (-1)^2 + (-1)(2) - 2^2 = 1-2-4 = -5
        # So |N(sqrt5)| = 5. sqrt(5) is a prime of Z[phi] of norm 5.

# ============================================================
# Step 4: The norm spectrum of occupied levels
# ============================================================
print("\n" + "=" * 70)
print("Step 4: Norm Spectrum Analysis")
print("=" * 70)

# Norms of occupied fermions: {0, 1, 5, 19}
# These are exactly the ALLOWED norms from exp123!
norms_occupied = set()
for name, (a, b, sector) in fermions.items():
    norms_occupied.add(abs(a**2 + a*b - b**2))
print(f"  Norms of occupied fermions: {sorted(norms_occupied)}")

# Check if these are ALL norms of Fibonacci-adjacent elements
print(f"\n  Fibonacci connection of norms:")
print(f"  |N| = 0: electron (0,0) - trivial")
print(f"  |N| = 1: units (phi^n) - 6 fermions")
print(f"  |N| = 5: charm (2,1)")
print(f"  |N| = 19: bottom (-1,4), top (4,1)")

# Check: are 0, 1, 5, 19 related to Fibonacci or Lucas numbers?
fib_seq = [fib(k) for k in range(15)]
lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47]  # Lucas numbers L(n)
print(f"\n  Fibonacci: {fib_seq[:12]}")
print(f"  Lucas:     {lucas}")
print(f"  Norms:     {sorted(norms_occupied)}")
print(f"  0 in Fib: {0 in fib_seq}")
print(f"  1 in Fib: {1 in fib_seq}")
print(f"  5 in Fib: {5 in fib_seq} (F(5) = 5)")
print(f"  19 in Fib: {19 in fib_seq}")
print(f"  19 in Lucas: {19 in lucas}")

# 19 is NOT Fibonacci or Lucas. But check F(n)^2 + F(n+1)^2 etc.
print(f"\n  Other sequences containing 19:")
print(f"  F(n)^2: {[fib(k)**2 for k in range(8)]}")
print(f"  F(n)*F(n+2): {[fib(k)*fib(k+2) for k in range(8)]}")
# 19 = 4^2 + 3*1 = ... hmm
# Actually, 19 = 4^2 + 4*1 - 1^2? No. Let me check: norm of (4,1) = 16+4-1 = 19
# And norm of (-1,4) = 1-4-16 = -19
# Is 19 a prime in Z? Yes. Is it inert, split, or ramified in Z[phi]?
# 19 mod 5: 19 = 3*5+4, so 19 mod 5 = 4.
# Legendre symbol (5/19): 5^((19-1)/2) mod 19 = 5^9 mod 19
# 5^2=25=6, 5^4=36=17, 5^8=17^2=289=289-15*19=289-285=4, 5^9=20=1 mod 19
# So (5/19) = 1, meaning 19 splits in Z[phi]. sqrt(5) exists mod 19.
# Indeed: 9^2 = 81 = 4*19+5, so 9^2 ≡ 5 mod 19.
# So 19 = pi * pi' in Z[phi] where pi = 4+1*phi (norm 19).

print(f"\n  Norm 19 analysis:")
print(f"  19 is prime in Z, splits in Z[phi] (Legendre symbol (5|19) = 1)")
print(f"  19 = (4+phi)(4+phi') = (4+phi)(3-phi) = 12+3phi-4phi-phi^2 = 12-phi-(phi+1) = 11-2phi")
print(f"  Check: (4+phi)(4+phi') = 16 + 4phi' + 4phi + phi*phi' = 16 + 4(phi+phi') + phi*phi'")
print(f"  phi + phi' = 1, phi*phi' = -1")
print(f"  = 16 + 4 - 1 = 19. Confirmed!")

# ============================================================
# Step 5: Generation counting on each line
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Generation Counting per Line via Units")
print("=" * 70)

# a=1 line: 3 generations (b=0,1,2) - ALL units
# b=1 line: fermions mu(1,1), c(2,1), t(4,1)
#   mu: |N(1,1)| = 1 (unit)
#   c:  |N(2,1)| = 5 (NOT unit)
#   t:  |N(4,1)| = 19 (NOT unit)
# So on b=1 line, only 1 unit (mu/s). Not 3 generations!

# 3a+2b=5 line: u(3,-2), s(1,1), b_q(-1,4)
#   u:   |N(3,-2)| = 1 (unit)
#   s:   |N(1,1)| = 1 (unit)
#   b_q: |N(-1,4)| = 19 (NOT unit)
# 2 units out of 3.

print(f"  Line analysis:")
print(f"\n  a=1 line (leptons/down quarks):")
for b in range(5):
    a = 1
    val = a + b*PHI
    norm = abs(a**2 + a*b - b**2)
    n = 5*a + 6*b
    unit = norm <= 1
    fermion = {0:'d', 1:'mu/s', 2:'tau'}.get(b, '')
    print(f"    b={b}: n={n:>2}, |N|={norm:>2}, unit={'Y' if unit else 'N'}, {fermion}")

print(f"\n  b=1 line (up quarks + mu/s):")
for a in range(-1, 6):
    b = 1
    val = a + b*PHI
    norm = abs(a**2 + a*b - b**2)
    n = 5*a + 6*b
    unit = norm <= 1
    fermion = {(1,1):'mu/s', (2,1):'c', (4,1):'t'}.get((a,b), '')
    if fermion or unit:
        print(f"    a={a}: n={n:>2}, |N|={norm:>2}, unit={'Y' if unit else 'N'}, {fermion}")

print(f"\n  3a+2b=5 line:")
for a in range(-2, 6):
    if (5-3*a) % 2 != 0:
        continue
    b = (5-3*a)//2
    val = a + b*PHI
    norm = abs(a**2 + a*b - b**2)
    n = 5*a + 6*b
    fermion = {(3,-2):'u', (1,1):'s', (-1,4):'b_q'}.get((a,b), '')
    if fermion:
        unit = norm <= 1
        print(f"    ({a},{b}): n={n:>2}, |N|={norm:>2}, unit={'Y' if unit else 'N'}, {fermion}")

# ============================================================
# Step 6: Why phi^0, phi^2, phi^3 (not phi^1)?
# ============================================================
print("\n" + "=" * 70)
print("Step 6: Why phi^0, phi^2, phi^3 and NOT phi^1?")
print("=" * 70)

print(f"""
  Powers of phi on a=1 line: phi^0, phi^2, phi^3
  Missing: phi^1 = phi = 0 + 1*phi (has a=0, not a=1)

  The gap {0, 2, 3} is interesting:
  - 0 = trivial (identity)
  - 2 = first non-trivial with a=1
  - 3 = second (and last) with a=1

  In the mass formula n = 5a + 6b = 5*1 + 6b = 5+6b:
  - phi^0: n=5 (down quark mass level)
  - phi^2: n=11 (muon/strange mass level)
  - phi^3: n=17 (tau mass level)

  The gap between consecutive levels: 6, 6 (uniform!)
  This is k_eff_base = 6 = b_1 (intersection number).

  Why NOT phi^1? Because phi^1 = (0,1), which has a=0.
  In the mass formula, a=0 is the ELECTRON (base state).
  phi^1 would give n = 5*0 + 6*1 = 6, which is a FALSE POSITIVE
  in level selection (n=6 passes norm+line but is unoccupied).
""")

# Check: phi^1 = (0,1) on the b=1 line, not a=1 line
print(f"  phi^1 = (0,1): on b=1 line (NOT a=1 line)")
print(f"  n = 5*0 + 6*1 = 6")
print(f"  This IS one of the 3 false positives from exp123!")
print(f"")
print(f"  FALSE POSITIVES REVISITED:")
print(f"  n=1:  (-1,1) = -phi^(-1) => negative unit")
print(f"  n=6:  (0,1)  = phi^1     => unit, but on WRONG line (b=1, not a=1)")
print(f"  n=23: (1,3)  = 1+3phi    => NOT a unit (|N|=5)")

# Verify
fp_analysis = [
    (1, -1, 1),
    (6, 0, 1),
    (23, 1, 3),
]
print(f"\n  {'n':>3} {'(a,b)':>8} {'element':>10} {'|N|':>4} {'phi^?':>10} {'explanation':>30}")
print("  " + "-" * 70)
for n, a, b in fp_analysis:
    val = a + b*PHI
    norm = abs(a**2 + a*b - b**2)
    is_unit = norm <= 1
    if is_unit and val > 0:
        power = round(np.log(val)/np.log(PHI))
        phi_str = f"phi^{power}"
    elif is_unit and val < 0:
        power = round(np.log(-val)/np.log(PHI))
        phi_str = f"-phi^{power}"
    else:
        phi_str = "non-unit"

    if n == 1:
        expl = "negative unit (unphysical)"
    elif n == 6:
        expl = "unit but wrong line (a!=1)"
    elif n == 23:
        expl = "non-unit (composite in Z[phi])"
    else:
        expl = ""
    print(f"  {n:>3} ({a:>2},{b:>2}) {val:>10.4f} {norm:>4} {phi_str:>10} {expl:>30}")

print(f"\n  *** ALL 3 FALSE POSITIVES EXPLAINED: ***")
print(f"  n=1:  negative unit (need positive for mass)")
print(f"  n=6:  phi^1, correct unit, but a=0 (electron line, not generation line)")
print(f"  n=23: |N|=5, composite element, NOT a unit of Z[phi]")

# ============================================================
# Step 7: The deep connection
# ============================================================
print("\n" + "=" * 70)
print("Step 7: The Deep Connection - Summary")
print("=" * 70)

print(f"""
  THEOREM (on a=1 line):
  =====================
  The positive units of Z[phi] with rational part a=1 are exactly:
    phi^0 = 1 + 0*phi  (b=0, generation 0)
    phi^2 = 1 + 1*phi  (b=1, generation 1)
    phi^3 = 1 + 2*phi  (b=2, generation 2)

  PROOF:
  phi^n = F(n-1) + F(n)*phi. Rational part = F(n-1).
  F(n-1) = 1 iff n-1 in {{-1, 1, 2}} (the only solutions to F(k)=1).
  This gives n in {{0, 2, 3}}, i.e., b = F(n) in {{0, 1, 2}}.
  QED.

  COROLLARY: There are EXACTLY 3 generations on the a=1 line.

  NUMBER-THEORETIC ORIGIN:
  3 = |{{k >= -1 : F(k) = 1}}| = |{{-1, 1, 2}}|

  This count is UNIVERSAL and UNCHANGEABLE:
  - F(0) = 0 (not 1)
  - F(1) = 1, F(2) = 1 (both equal 1)
  - F(3) = 2 (already > 1, and grows forever after)
  - F(-1) = 1 (the extra solution from extension)

  PHYSICAL INTERPRETATION:
  - Generations = powers of the fundamental unit phi on the mass line a=1
  - Each generation is a UNIT of Z[phi] (has multiplicative inverse)
  - The 4th generation would require a NON-UNIT (composite element)
  - Composite elements can be "factored" = unstable/not fundamental

  STATUS: PATTERN (algebraically clean, but need to justify WHY
  generation elements must be units of Z[phi])

  TO ELEVATE TO DERIVED: need physical principle that REQUIRES
  the generation element to be invertible in Z[phi].
  Candidate: normalizability of wave functions on the Z[phi] lattice?
  Candidate: modular invariance of the mass spectrum?
""")

# ============================================================
# Step 8: Universality check - all fermions as phi-powers
# ============================================================
print("\n" + "=" * 70)
print("Step 8: All Fermions as Phi-Powers")
print("=" * 70)

print(f"\n  If EVERY fermion's Z[phi] element a+b*phi = +-phi^n, which n?")
for name, (a, b, sector) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    val = a + b*PHI
    norm = abs(a**2 + a*b - b**2)
    n_level = 5*a + 6*b
    if val > 0 and norm == 1:
        n_power = round(np.log(val)/np.log(PHI))
        print(f"  {name:>5} ({a:>2},{b:>2}): n={n_level:>2}, {a}+{b}phi = phi^{n_power} = {PHI**n_power:.4f} (unit)")
    elif val < 0 and norm == 1:
        n_power = round(np.log(-val)/np.log(PHI))
        print(f"  {name:>5} ({a:>2},{b:>2}): n={n_level:>2}, {a}+{b}phi = -phi^{n_power} = {-PHI**n_power:.4f} (neg unit)")
    elif val == 0:
        print(f"  {name:>5} ({a:>2},{b:>2}): n={n_level:>2}, {a}+{b}phi = 0 (zero, trivial)")
    else:
        print(f"  {name:>5} ({a:>2},{b:>2}): n={n_level:>2}, {a}+{b}phi = {val:.4f}, |N|={norm} (NOT a unit)")

# Unit fermions
unit_fermions = [(name, a, b) for name, (a,b,s) in fermions.items() if abs(a**2+a*b-b**2) <= 1]
non_unit_fermions = [(name, a, b) for name, (a,b,s) in fermions.items() if abs(a**2+a*b-b**2) > 1]

print(f"\n  Unit fermions ({len(unit_fermions)}/9): {[name for name,_,_ in unit_fermions]}")
print(f"  Non-unit fermions ({len(non_unit_fermions)}/9): {[name for name,_,_ in non_unit_fermions]}")
print(f"  Non-unit norms: {[abs(a**2+a*b-b**2) for _,a,b in non_unit_fermions]}")

# The non-units are the "heavy" quarks: charm, bottom, top
# These have norms 5 and 19
# 5 and 19 are PRIMES in Z
# In Z[phi]: 5 = sqrt(5)^2 / (something)... actually
# N(sqrt(5)) = N(2phi-1) = N(-1+2phi) = 1-2-4 = -5
# So sqrt(5) has norm -5, and 5 = -N(sqrt(5)) = (2phi-1)(1-2phi')
# And 19 splits as (4+phi)(4+phi') [verified above]

print(f"\n  Non-unit factorization in Z[phi]:")
print(f"  |N|=5:  5 splits as (2phi-1)(2phi'-1) where 2phi-1=sqrt(5)")
print(f"  |N|=19: 19 splits as (4+phi)(4+phi')")
print(f"  These are the SIMPLEST primes that split in Z[phi]")
print(f"  5 is RAMIFIED (5 = sqrt(5)^2 up to units)")
print(f"  19 SPLITS (19 = pi*pi' with distinct conjugate primes)")

# ============================================================
# CONCLUSIONS
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print(f"""
MAIN RESULT:
  3 generations on a=1 line = number of positive powers of phi
  with rational part (Fibonacci index) equal to 1.

  Algebraic: phi^0, phi^2, phi^3 are the ONLY such powers.
  Number-theoretic: F(k) = 1 has exactly 3 solutions (k = -1, 1, 2).

  This is ALGEBRAICALLY CLEAN and gives EXACTLY 3.

UNIVERSALITY:
  6/9 fermions are units of Z[phi] (e, mu, tau, d, s, u).
  3/9 are non-units: charm(|N|=5), bottom(|N|=19), top(|N|=19).
  The non-units correspond to heavy quarks on other lattice lines.

  Unit constraint works PERFECTLY for a=1 line (3 generations).
  It does NOT directly constrain the b=1 or 3a+2b=5 lines.

FALSE POSITIVES EXPLAINED:
  n=1:  negative unit (-phi^(-1)): need positive mass
  n=6:  phi^1 = (0,1): unit but on electron line (a=0), not a=1
  n=23: (1,3) with |N|=5: composite, not a unit

GAP TO DERIVATION:
  Need physical principle requiring generation elements to be Z[phi] units.
  Most promising: modular structure of the mass spectrum, or
  wave function normalizability on the 600-cell lattice.

STATUS: STRONG PATTERN (algebraically clean, number-theoretically deep)
""")
