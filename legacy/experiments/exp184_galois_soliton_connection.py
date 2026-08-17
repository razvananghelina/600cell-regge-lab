"""
EXP-184: Galois Balance and Soliton Stability

Key discoveries from exp183:
  - Units on a=1 line have |z/z'| = phi^k (EXACT integer powers)
  - z = 1+b*phi = phi^{0,2,3} for b={0,1,2} (powers of phi!)
  - Non-units DON'T have phi-power Galois ratios

Goals:
  1. Show z = phi^k explicitly and connect to Fibonacci
  2. Compute the "shadow level" n_T under Galois conjugation
  3. Find a physical principle why ONLY units are stable solitons
  4. Connect everything into a single derivation chain
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PHI_C = (1 - np.sqrt(5)) / 2  # Galois conjugate

# ============================================================
# PART 1: Elements on a=1 line are powers of phi
# ============================================================
print("=" * 70)
print("PART 1: POWERS OF PHI ON THE a=1 LINE")
print("=" * 70)

# phi^k = F(k-1) + F(k)*phi where F is Fibonacci
# We define F(-1)=1, F(0)=0, F(1)=1, F(2)=1, F(3)=2, ...
def fib(n):
    """Extended Fibonacci: F(-1)=1, F(0)=0, F(1)=1, ..."""
    if n == -1: return 1
    if n == 0: return 0
    if n == 1: return 1
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b

print("\nphi^k = F(k-1) + F(k)*phi, with (a,b) = (F(k-1), F(k)):")
print(f"  {'k':>3} {'phi^k':>10} {'F(k-1)=a':>10} {'F(k)=b':>8} {'a=1?':>6} {'n=5a+6b':>8}")
print("  " + "-" * 52)
for k in range(0, 12):
    pk = PHI ** k
    a_k = fib(k-1)
    b_k = fib(k)
    on_a1 = "<<< YES" if a_k == 1 else ""
    n = 5*a_k + 6*b_k
    print(f"  {k:3d} {pk:10.4f} {a_k:10d} {b_k:8d} {on_a1:>7} {n:8d}")

print()
print("  phi^k lies on the a=1 line when F(k-1) = 1.")
print("  F(k-1) = 1 for k-1 in {-1, 1, 2} (since F(-1)=F(1)=F(2)=1)")
print("  => k in {0, 2, 3}")
print("  => b = F(k) = F(0), F(2), F(3) = 0, 1, 2")
print()
print("  THIS IS THE FIBONACCI CONDITION: F(m)=1 has exactly 3 solutions")
print("  for m >= -1. After m=2, Fibonacci grows: F(3)=2, F(4)=3, ...")
print("  So F(m)=1 NEVER returns. EXACTLY 3 powers of phi lie on a=1.")

# ============================================================
# PART 2: Galois Shadow Levels
# ============================================================
print()
print("=" * 70)
print("PART 2: GALOIS SHADOW LEVELS")
print("=" * 70)

print("""
For z = a + b*phi, the Galois conjugate is z' = a + b*phi'.
Since phi' = 1 - phi, we get z' = (a+b) - b*phi = a' + b'*phi
where a' = a+b, b' = -b.

The shadow level: n' = 5*a' + 6*b' = 5(a+b) + 6(-b) = 5a + 5b - 6b = 5a - b.
On the a=1 line: n' = 5 - b.
""")

print("Shadow levels on a=1 line:")
print(f"  {'b':>3} {'n_S':>5} {'n_T':>5} {'n_S+n_T':>8} {'n_S-n_T':>8} {'n_S*n_T':>8} {'|N|':>5}")
for b in range(0, 6):
    n_S = 5 + 6*b
    n_T = 5 - b
    N = 1 + b - b**2
    print(f"  {b:3d} {n_S:5d} {n_T:5d} {n_S+n_T:8d} {n_S-n_T:8d} {n_S*n_T:8d} {abs(N):5d}")

print()
print("  n_S + n_T = 10 + 5b (always)")
print("  n_S - n_T = 7b (always)")
print("  n_S * n_T = (5+6b)(5-b)")
print()
print("  Shadow becomes NEGATIVE at b >= 6: n_T = 5-b < 0")
print("  But the interesting constraint is earlier...")

# ============================================================
# PART 3: DSI Level Compatibility
# ============================================================
print()
print("=" * 70)
print("PART 3: DSI LEVEL COMPATIBILITY IN BOTH SECTORS")
print("=" * 70)

print("""
The DSI potential V = A*sin^2(pi*ln|psi|/ln(phi)) has minima at ALL
integer levels n. But the COMPLEXITY MINIMIZATION says that only
levels with (a,b) = argmin|a|+|b| s.t. 5a+6b=n are realized.

For the shadow level n_T = 5-b, what is the minimal (a',b')?
""")

# For each b, find minimal complexity of shadow level
print("  Shadow level complexity:")
print(f"  {'b':>3} {'n_T':>5} {'best (a,b)':>12} {'|a|+|b|':>8} {'on a=1?':>8}")
for b in range(0, 6):
    n_T = 5 - b
    # Find argmin |a|+|b| s.t. 5a+6b=n_T
    best = None
    best_cost = 999
    for a2 in range(-20, 21):
        for b2 in range(-20, 21):
            if 5*a2 + 6*b2 == n_T:
                cost = abs(a2) + abs(b2)
                if cost < best_cost:
                    best_cost = cost
                    best = (a2, b2)
    on_a1 = "yes" if best and best[0] == 1 else "no"
    print(f"  {b:3d} {n_T:5d} {str(best):>12} {best_cost:8d} {on_a1:>8}")

# ============================================================
# PART 4: The Cross-Sector Energy
# ============================================================
print()
print("=" * 70)
print("PART 4: CROSS-SECTOR ENERGY ARGUMENT")
print("=" * 70)

print("""
In the E8 lattice, S and T sectors are coupled through the
Galois automorphism sigma: phi -> phi'.

Hypothesis: A soliton at level z = 1+b*phi in sector S has a
"shadow" at z' = (1+b) - b*phi in sector T. The TOTAL energy
of the soliton has two contributions:

  E_total = E_soliton(S) + E_coupling(S,T)

where E_coupling depends on the norm N(z) = z * z'.

For UNITS (|N|=1): E_coupling is minimal (cross terms cancel).
For NON-UNITS (|N|>1): E_coupling adds |N|-1 extra energy.
""")

# The soliton energy is E_flip = 3 (uniform). But could there be
# a Galois correction?
#
# In the star graph S_3, the perturbation eigenvalues are {0, -1, -4}.
# The TOTAL perturbation energy is Tr(DeltaL) = -3-1-1-1 = -6.
# This splits as: 0 (bulk) + 2*(-1) (shell) + (-4) (core) = -6.
#
# For the Galois partner: the shadow soliton in T has the SAME
# local structure (E8 is symmetric under S<->T modulo scaling).
# The cross-sector coupling would add terms proportional to N(z).

# Can we show that |N|>1 destabilizes the soliton?
# The star eigenvalue at the core is -(k+1) = -4.
# If the Galois coupling adds |N|-1 to the core energy:
# E_core = 4 + (|N|-1) = 3 + |N|
# For units: E_core = 4 (standard)
# For |N|=5: E_core = 8 (costs double!)
# For |N|=19: E_core = 22 (very costly)

print("  If E_core = (k+1) + (|N|-1) = k + |N|:")
print(f"  {'b':>3} {'|N|':>5} {'E_core':>8} {'E_flip=3?':>10}")
for b in range(0, 5):
    N_abs = abs(1 + b - b**2)
    E_core = 3 + N_abs  # k=3
    stable = "STABLE" if N_abs <= 1 else "COSTLY"
    print(f"  {b:3d} {N_abs:5d} {E_core:8d} {stable:>10}")

print()
print("  For units: E_core = 4 = standard soliton energy")
print("  For b=3: E_core = 8 = double cost -> energetically forbidden")
print("  SPECULATIVE: needs derivation of E_coupling = |N|-1")

# ============================================================
# PART 5: The Key Identity z = phi^k
# ============================================================
print()
print("=" * 70)
print("PART 5: THE POWER-OF-PHI STABILITY CRITERION")
print("=" * 70)

print("""
OBSERVATION: On the a=1 line, the ALLOWED generations correspond
to elements that are POWERS OF PHI:

  b=0: z = 1       = phi^0    (ALLOWED - electron/down)
  b=1: z = 1+phi   = phi^2    (ALLOWED - muon/strange)
  b=2: z = 1+2*phi = phi^3    (ALLOWED - tau/bottom)
  b=3: z = 1+3*phi = ???      (FORBIDDEN)

Is 1+3*phi a power of phi? NO.
""")

val = 1 + 3*PHI
k_test = np.log(val) / np.log(PHI)
print(f"  1 + 3*phi = {val:.6f}")
print(f"  log_phi(1+3*phi) = {k_test:.6f}")
print(f"  Nearest integer: {round(k_test)}")
print(f"  phi^{round(k_test)} = {PHI**round(k_test):.6f}")
print(f"  Difference: {abs(val - PHI**round(k_test)):.6f}")
print(f"  1+3*phi is NOT a power of phi (off by {abs(val - PHI**round(k_test)):.4f})")
print()

# Check: which integer combinations a+b*phi are powers of phi?
print("  Powers of phi and their (a,b) decomposition:")
print(f"  {'k':>3} {'phi^k':>10} {'a':>4} {'b':>4} {'Unit?':>6} {'Line':>8}")
for k in range(-3, 10):
    pk = PHI ** k
    a_k = fib(k-1)
    b_k = fib(k)
    N_k = a_k**2 + a_k*b_k - b_k**2
    # N(phi^k) = (-1)^k, so always unit
    unit = "YES" if abs(N_k) <= 1 else "no"
    line = f"a={a_k}"
    print(f"  {k:3d} {pk:10.4f} {a_k:4d} {b_k:4d} {unit:>6} {line:>8}")

print()
print("  ALL powers of phi are units (|N(phi^k)| = |(-1)^k| = 1). ALWAYS.")
print("  Conversely: NOT all units are powers of phi (e.g., 2+phi is a unit")
print("  but 2+phi = phi^2 + 1 = phi^2 + phi^0... it IS phi^3 + phi^(-1))")
print()

# Actually, in Z[phi], the ONLY units are +/- phi^k.
# So z is a unit iff z = +/- phi^k for some integer k.
print("  THEOREM (Dirichlet unit theorem for Z[phi]):")
print("  The units of Z[phi] are EXACTLY {+/- phi^k : k in Z}.")
print("  So 'z is a unit' is equivalent to 'z is a power of phi (up to sign)'.")
print()

# ============================================================
# PART 6: The Complete Argument
# ============================================================
print("=" * 70)
print("PART 6: THE COMPLETE ARGUMENT")
print("=" * 70)

print("""
PROPOSITION: Fermion generations on the a=1 line correspond to
elements z = 1+b*phi that are UNITS in Z[phi], i.e., powers of phi.

ALGEBRAIC PROOF that this gives exactly 3:

  1. By Dirichlet's unit theorem, units of Z[phi] = {+/- phi^k : k in Z}.

  2. phi^k = F(k-1) + F(k)*phi where F is the Fibonacci sequence.

  3. The a=1 line requires F(k-1) = 1.

  4. F(m) = 1 has solutions m in {-1, 1, 2} (for m >= -1).
     For m >= 3: F(m) >= 2 (Fibonacci is strictly increasing for m >= 1).
     This gives k in {0, 2, 3}.

  5. Therefore b = F(k) = F(0), F(2), F(3) = 0, 1, 2.

  QED: Exactly 3 generations on the a=1 line.

WHY THE NUMBER 3 IS IMMUTABLE:
  - It equals the number of Fibonacci values equal to 1
  - F(m) = 1 for m = -1, 1, 2 (and only these)
  - This is a number-theoretic fact, independent of any physical parameter
  - Even if the 600-cell were different, as long as the mass formula
    uses phi (the golden ratio), you get exactly 3 generations on
    any "a = 1" type line.
""")

# ============================================================
# PART 7: Why must fermions be units?
# ============================================================
print("=" * 70)
print("PART 7: WHY MUST FERMIONS BE UNITS? (Candidate arguments)")
print("=" * 70)

print("""
The remaining question: WHY must z be a unit of Z[phi]?

ARGUMENT A: DSI Consistency (STRONGEST)
  The DSI potential V = A*sin^2(pi*ln|psi|/ln(phi)) has minima at phi^n.
  A fermion at level n is a soliton sitting at a minimum.
  In the E8 picture, the soliton exists in BOTH Galois sectors.

  In sector S: mass = m_e * phi^n_S where n_S = 5 + 6b
  In sector T: mass = m_e * phi^n_T where n_T = 5 - b (Galois shadow)

  The soliton is described by z = 1+b*phi in S and z' = (1+b)-b*phi in T.
  The product z*z' = N(z) = norm.

  For the COMBINED soliton (S x T) to sit at a DSI minimum,
  we need BOTH n_S and n_T to be valid DSI levels.

  Condition: z and z' must BOTH be "DSI-compatible", meaning
  they can be expressed as products of the DSI scaling factor phi.

  z is DSI-compatible iff z = +/- phi^k iff z is a UNIT.

  For non-units: z = phi^k * (irreducible factor). The irreducible
  factor is NOT a power of phi, so z does NOT sit at a pure DSI minimum.
""")

# Verify: for each generation, show z and z' are both powers of phi
print("  Verification: z and z' as powers of phi:")
print(f"  {'b':>3} {'z':>12} {'z_as_phi^k':>12} {'z_bar':>12} {'zbar_as_phi^k':>14}")
for b in range(0, 5):
    z = 1 + b * PHI
    zbar = 1 + b * PHI_C

    # Find k such that phi^k = z (or -z)
    k_z = None
    k_zbar = None
    for k in range(-10, 15):
        if abs(PHI**k - z) < 1e-8:
            k_z = k
        if abs(PHI**k - (-z)) < 1e-8:
            k_z = f"-phi^{k}"
        if abs(PHI**k - zbar) < 1e-8:
            k_zbar = k
        if abs(PHI**k - (-zbar)) < 1e-8:
            k_zbar = f"-phi^{k}"

    z_str = f"phi^{k_z}" if isinstance(k_z, int) else str(k_z) if k_z else "NOT phi^k"
    zbar_str = f"phi^{k_zbar}" if isinstance(k_zbar, int) else str(k_zbar) if k_zbar else "NOT phi^k"
    marker = "   <<<" if isinstance(k_z, int) or (k_z and "phi" in str(k_z)) else ""
    print(f"  {b:3d} {z:12.6f} {z_str:>12} {zbar:12.6f} {zbar_str:>14}{marker}")

print()

# ============================================================
# PART 8: Extend to ALL fermions
# ============================================================
print("=" * 70)
print("PART 8: ALL FERMIONS - UNIT ANALYSIS")
print("=" * 70)

fermions = [
    ('electron', 0, 0, 0),
    ('muon',     1, 1, 11),
    ('tau',      1, 2, 17),
    ('up',       3,-2, 3),
    ('charm',    2, 1, 16),
    ('top',      4, 1, 26),
    ('down',     1, 0, 5),
    ('strange',  1, 1, 11),
    ('bottom',  -1, 4, 19),
]

print(f"\n  {'Fermion':>10} {'(a,b)':>7} {'z=a+b*phi':>10} {'|N|':>5} {'Unit':>5} {'z=phi^k?':>12}")
print("  " + "-" * 55)
for name, a, b, n in fermions:
    z = a + b * PHI
    N = a**2 + a*b - b**2
    is_unit = abs(N) <= 1

    # Check if z = +/- phi^k
    phi_k = None
    for k in range(-10, 30):
        if abs(PHI**k - abs(z)) < 1e-6:
            phi_k = k
            break

    k_str = f"phi^{phi_k}" if phi_k is not None else "no"
    if z < 0 and phi_k is not None:
        k_str = f"-phi^{phi_k}"

    print(f"  {name:>10} ({a:2d},{b:2d}) {z:10.4f} {abs(N):5d} {'YES' if is_unit else 'no':>5} {k_str:>12}")

print()

# PART 9: The generation selection on OTHER lines
print("=" * 70)
print("PART 9: GENERATION COUNTING ON OTHER LINES")
print("=" * 70)

print("""
On the a=1 line: 3 generations (b=0,1,2), all units.
What about other lines?

Line b=1: charm(2,1), muon/strange(1,1), top(4,1)
  z = a + phi for a = 1, 2, 4.
  N(a,1) = a^2 + a - 1.
""")

print("Line b=1 (fixed b, varying a):")
print(f"  {'a':>3} {'n':>5} {'z':>10} {'N':>6} {'|N|':>5} {'Fermion':>12}")
for a in range(-2, 8):
    n = 5*a + 6
    z = a + PHI
    N = a**2 + a - 1
    name = {1: 'mu/s', 2: 'charm', 4: 'top'}.get(a, '')
    print(f"  {a:3d} {n:5d} {z:10.4f} {N:6d} {abs(N):5d} {name:>12}")

print()
print("Line a=1 summary: units at b={0,1,2}, first non-unit at b=3 (|N|=5)")
print("Line b=1 summary: unit at a=1 (|N|=1), non-units elsewhere")
print("  a=0: |N|=1 (unit, but n=6 is EMPTY)")
print("  a=1: |N|=1 (mu/s)")
print("  a=2: |N|=5 (charm, NON-UNIT)")
print("  a=4: |N|=19 (top, NON-UNIT)")
print()
print("  On b=1: only a={0,1} are units. a=0 gives n=6 (empty).")
print("  So only 1 occupied unit on this line (muon/strange).")

# ============================================================
# PART 10: FINAL SYNTHESIS
# ============================================================
print()
print("=" * 70)
print("FINAL SYNTHESIS: THE GENERATION THEOREM")
print("=" * 70)

print("""
THEOREM: The number of fermion generations on the a=1 line of
the DSI ladder is EXACTLY 3.

PROOF (combining algebraic and topological arguments):

ROUTE 1 (Algebraic - Fibonacci):
  (i)   Mass levels on a=1: n = 5 + 6b, element z = 1 + b*phi.
  (ii)  z is a unit of Z[phi] iff z = +-phi^k (Dirichlet).
  (iii) phi^k = F(k-1) + F(k)*phi lies on a=1 iff F(k-1) = 1.
  (iv)  F(m) = 1 for m = {-1, 1, 2} only (Fibonacci grows for m >= 3).
  (v)   => k = {0, 2, 3} => b = {0, 1, 2} => 3 generations.

ROUTE 2 (Topological - Star Graph):
  (i)   600-cell has proper 5-coloring (all edges inter-coset).
  (ii)  Soliton support = star graph S_3 (leaves not connected).
  (iii) Star with k >= 2 has exactly 3 distinct eigenvalues.
  (iv)  => 3 spectral sub-levels: bulk, shell, core.

ROUTE 3 (Dimensional - 4-to-3):
  (i)   96 C-vertices split into 4 groups of 24 by zero coordinate.
  (ii)  Temporal axis selection removes 1 group.
  (iii) => 3 observable groups.

CONNECTION (Routes 1 <-> 2):
  - Both give 3 because of phi's minimal polynomial x^2-x-1.
  - Route 1: Fibonacci condition F(m)=1, controlled by phi.
  - Route 2: k = 12/4 = 3 neighbors per coset, controlled by the
    600-cell geometry which is built from phi (edge angle = pi/5).
  - The identity a_1 = 5 = disc(Z[phi]) = |2I/2T| unifies them:
    the graph structure and the algebraic norm are both controlled
    by the discriminant of phi's number field.

WHAT REMAINS PATTERN:
  - WHY must generations be units? The DSI consistency argument
    (both Galois sectors must sit at phi-power minima) is
    physically motivated but not rigorously derived.
  - The mapping sub-level <-> generation index b is qualitative
    (localization matches mass hierarchy) but not quantitative
    (shift ratios are not phi^6).

CLASSIFICATION:
  N_gen = 3 on a=1 line: DERIVED (both routes independent)
  N_gen = 3 from soliton splitting: DERIVED (star graph theorem)
  Unit condition => 3 generations: DERIVED (Fibonacci/Dirichlet)
  Physical necessity of unit condition: PATTERN (best candidate:
    DSI consistency in both E8 Galois sectors)
  Mapping sub-levels to specific masses: PATTERN
""")

# Scorecard
print("SCORECARD UPDATE:")
print("  DERIVED:")
print("    - 3 = Fibonacci solutions F(m)=1 (exact, number-theoretic)")
print("    - 3 = star graph eigenvalues (exact, graph-theoretic)")
print("    - 3 = 12/4 = neighbors per coset (exact, combinatorial)")
print("    - 5 = disc = a_1 = cosets (exact, algebraic)")
print("    - norm gap 1 -> 5 specific to phi (verified for other fields)")
print("    - units of Z[phi] = +-phi^k (Dirichlet unit theorem)")
print("    - phi^k on a=1 iff F(k-1)=1 (Fibonacci, 3 solutions)")
print("    - star topology forced by proper 5-coloring (no intra-coset edges)")
print()
print("  PATTERN:")
print("    - generation = unit (physical principle not derived)")
print("    - sub-level = generation (qualitative localization match)")
print("    - DSI consistency in both Galois sectors (motivated but unproven)")
print()
print("  NEGATIVE:")
print("    - shift ratios do NOT give mass ratios phi^6")
print("    - defect amplitudes do NOT match mass hierarchy quantitatively")
