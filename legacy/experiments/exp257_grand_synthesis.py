"""
EXP-257: Grand Synthesis - E8 Exponents, Coxeter Pairs, and Mass Structure
===========================================================================
Ties together all discoveries from exp253-256 into a unified picture.
Explores: down quark {6,14}, unused E8 exponents, mass corrections,
and the complete Coxeter pair assignment.
"""
import math

PHI = (1 + math.sqrt(5)) / 2
a_1 = 5
b_1 = 6
h = 30
N_gen = 3
rank_E8 = 8

e8_exponents = [1, 7, 11, 13, 17, 19, 23, 29]
coxeter_pairs = [(1,29), (7,23), (11,19), (13,17)]

# Fermion data
fermions = [
    ('e',   0,  'lep', 0), ('mu',  11, 'lep', 1), ('tau', 17, 'lep', 2),
    ('u',   3,  'up',  0), ('c',   16, 'up',  1), ('t',   26, 'up',  2),
    ('d',   5,  'dn',  0), ('s',   11, 'dn',  1), ('b',   19, 'dn',  2),
]

sector_bases = {'lep': 0, 'up': 3, 'dn': 5}
gen_base = [0, 11, 17]

m_e = 0.51099895
# Experimental masses (MeV)
m_exp = {'e': 0.511, 'mu': 105.658, 'tau': 1776.86,
         'u': 2.16, 'c': 1270, 't': 172500,
         'd': 4.67, 's': 93.4, 'b': 4180}

print("=" * 70)
print("EXP-257: GRAND SYNTHESIS")
print("=" * 70)

# =====================================================================
# PART 1: Down quark shifts {6, 14} = 2 * {3, 7}
# =====================================================================
print()
print("PART 1: Down quark generation shifts")
print("-" * 50)
print()

# Down n-values: {5, 11, 19}
# Shifts from sector base (5): {0, 6, 14}
print("  Down quark n-values: {5, 11, 19}")
print("  Shifts from base (a_1=5): {0, 6, 14}")
print()

# Key observation: 6 = 2*3 and 14 = 2*7
print("  6 = 2 * 3   where 3 = N_gen = sector_base(up)")
print("  14 = 2 * 7   where 7 = E8 exponent!")
print()
print("  So down_shift = 2 * {0, 3, 7}")
print("  {3, 7}: sum = 10 = 2*a_1")
print("  {3, 7}: product = 21 = 3*7")
print()

# Is 7 special? 7 is the 2nd E8 exponent
# 7 = a_1 + 2 = a_1 + F(3)
# 7 is the PARTNER of 23: 7 + 23 = 30 = h
# 23 is used by up quarks (n_top - 3 = 23)
# So 7 = h - (n_top - base_up) = h - 23
print("  7 = a_1 + 2 = a_1 + F(3)")
print("  7 + 23 = 30 = h (Coxeter partner of 23)")
print("  23 = n_top - base_up (up quark uses 23)")
print("  So 7 is the COXETER PARTNER of the top quark's generation shift!")
print()

# And 3 = N_gen = base_up
# 3 is NOT an E8 exponent (E8 exps are coprime to h=30, gcd(3,30)=3)
# But 3 = a_1 - F(3) = base offset for up quarks
print("  3 = N_gen = base_up = a_1 - F(3)")
print("  3 is NOT an E8 exponent (gcd(3,30) = 3, not coprime)")
print("  But 3 + 11 = 14 (down g=2 shift = 2*7)")
print()

# COMPLETE PICTURE:
print("  COMPLETE DECOMPOSITION:")
print("  Sector    Shifts (from 0)    Shifts (from base)")
print("  Lepton    {0, 11, 17}        {0, 11, 17}")
print("  Up        {3, 16, 26}        {0, 13, 23}")
print("  Down      {5, 11, 19}        {0,  6, 14}")
print()
print("  From base, factored:")
print("  Lepton: {0, 11, 17} = {0, E8, E8}")
print("  Up:     {0, 13, 23} = {0, E8, E8}")
print("  Down:   {0, 2*3, 2*7} = 2*{0, N_gen, E8}")
print()
print("  Down quarks have a FACTOR OF 2 in their generation shifts!")
print("  This factor 2 = F(3) = Fibonacci third.")

# =====================================================================
# PART 2: Complete Coxeter pair assignment
# =====================================================================
print()
print("=" * 70)
print("PART 2: Complete Coxeter pair assignment")
print("=" * 70)
print()

print("  4 Coxeter pairs of E8: (1,29), (7,23), (11,19), (13,17)")
print()
print("  ASSIGNMENT (via generation shifts from sector base):")
print()

# Pair (11, 19):
print("  Pair (11, 19): h-complementary")
print("    11 -> lepton g=1 (muon), also down g=1 (strange)")
print("    19 -> direct: n_bottom = 19 (down g=2 before base subtraction)")
print("    19 - base_dn = 14 = 2*7 (down g=2 shift)")
print("    SERVES: leptons g=1, down quarks g=1+g=2")
print()

# Pair (13, 17):
print("  Pair (13, 17): h-complementary")
print("    17 -> lepton g=2 (tau)")
print("    13 -> up g=1 (charm): n_c - base_up = 16 - 3 = 13")
print("    SERVES: leptons g=2, up quarks g=1")
print("    >>> THIS PAIR CONNECTS TAU AND CHARM <<<")
print()

# Pair (7, 23):
print("  Pair (7, 23): h-complementary")
print("    23 -> up g=2 (top): n_t - base_up = 26 - 3 = 23")
print("    7  -> down: 2*7 = 14 = down g=2 shift")
print("    SERVES: up quarks g=2, down quarks g=2 (via doubling)")
print("    >>> THIS PAIR CONNECTS TOP AND BOTTOM <<<")
print()

# Pair (1, 29):
print("  Pair (1, 29): h-complementary")
print("    1 -> electron? n_e = 0 (off by 1)")
print("    29 -> ? No direct fermion use")
print("    BUT: 29 = h - 1 = a_1*b_1 - 1")
print("    AND: 1 + 29 = 30 = h (trivially)")
print()

# Check: does pair (1,29) have ANY role?
print("  Pair (1,29) - the 'vacuum pair':")
print("    n_e = 0 = 1 - 1 (ground state, one below trivial exponent)")
print("    base_up = 3 = 1 + F(3) (one above trivial + Fibonacci)")
print("    base_dn = 5 = 1 + 4 = 1 + (a_1-1)")
print("    29 = h - 1 = highest non-trivial E8 exponent")
print()

# Is 29 hidden somewhere?
for name, n, sec, g in fermions:
    if abs(n - 29) <= 3:
        print("    %s: n=%d, dist to 29 = %d" % (name, n, abs(n-29)))
print("    n_top = 26 = 29 - 3 = 29 - N_gen")
print("    So: pair (1,29) brackets the ENTIRE mass spectrum: 0..26")
print("    And (29-1)/2 = 14 = down g=2 shift!")

# =====================================================================
# PART 3: Coxeter pair - sector assignment table
# =====================================================================
print()
print("=" * 70)
print("PART 3: Coxeter pair - sector assignment table")
print("=" * 70)
print()

print("  Each Coxeter pair (m, h-m) contributes to fermion sectors:")
print()
print("  Pair       Element  Sector     Generation  Role")
print("  " + "-" * 55)
print("  (11, 19)   11       Lepton     g=1         muon shift")
print("             11       Down       g=1         strange (direct)")
print("             19       Down       g=2         bottom (direct)")
print("  (13, 17)   17       Lepton     g=2         tau shift")
print("             13       Up         g=1         charm shift")
print("  (7, 23)    23       Up         g=2         top shift")
print("             7        Down       g=2         2*7=14 (doubled)")
print("  (1, 29)    1        ---        ---         vacuum/ground")
print("             29       ---        ---         h-1 boundary")
print()

# Count usage:
print("  Usage count per pair:")
print("    (11,19): 3 fermion connections (muon, strange, bottom)")
print("    (13,17): 2 fermion connections (tau, charm)")
print("    (7,23):  2 fermion connections (top, bottom via 2*7)")
print("    (1,29):  0 direct, brackets spectrum")
print()
print("  Total: 7 connections for 6 non-ground-state fermions")
print("  (muon/strange share n=11, so count as 1 E8 connection)")

# =====================================================================
# PART 4: The doubling phenomenon for down quarks
# =====================================================================
print()
print("=" * 70)
print("PART 4: Why down quarks have factor 2")
print("=" * 70)
print()

# down_shift = 2 * {0, 3, 7}
# But we know delta_down(g) = F(3) + N_gen*sigma(g) - g*(g+1)
# = {5, 0, 2}

# And down n from base = n_base(g) + delta_down(g)
# = {0+5, 11+0, 17+2} = {5, 11, 19}
# Shifts from 5: {0, 6, 14}

# Actually, from n=0: {5, 11, 19}
# 5 = a_1
# 11 = a_1 + b_1 = a_1 + (a_1+1) = 2*a_1+1
# 19 = a_1^2 - b_1 = a_1^2 - a_1 - 1

# Alternative: from LEPTON n-values:
# n_down - n_lepton per generation:
print("  n_down - n_lepton (sector offset per generation):")
for g in range(3):
    n_l = gen_base[g]
    n_d = gen_base[g] + [5,0,2][g]
    diff = n_d - n_l
    print("    g=%d: n_down=%d - n_lep=%d = %d = delta_down(%d)" % (g, n_d, n_l, diff, g))

print()
print("  delta_down = {5, 0, 2}")
print("  These come from: F(3) + N_gen*sigma(g) - g*(g+1)")
print("  sigma = (01)(2): swaps generations 0 and 1")
print()

# The down-quark generation shifts (from base) are:
# 0, 6, 14
# = 0, 2*3, 2*7
# But this is DERIVED from n_base + delta_down - base_down
# = (0+5-5), (11+0-5), (17+2-5) = (0, 6, 14)

# Can we see the factor of 2 from the derivation?
# shift(g=1) = 11 + 0 - 5 = 6 = 2*3 = 2*N_gen
# shift(g=2) = 17 + 2 - 5 = 14 = 2*7 = 2*(a_1+2)

# shift(g=1) = n_base(1) + delta_down(1) - delta_down(0)
# = (a_1+b_1) + 0 - a_1 = b_1 = 6 = 2*N_gen (since b_1 = 2*N_gen for a_1=5)
print("  shift(g=1) = b_1 = 6 = 2*N_gen (since b_1 = 2*N_gen iff a_1=5)")
print("  Verify: b_1 = a_1+1 = 6, N_gen = 3, 2*3 = 6. EXACT for a_1=5.")
print()

# shift(g=2) = n_base(2) + delta_down(2) - delta_down(0)
# = (a_1+2*b_1) + 2 - a_1 = 2*b_1 + 2 = 2*(b_1+1) = 2*7
# b_1+1 = 7 = E8 exponent!
print("  shift(g=2) = 2*b_1 + 2 = 2*(b_1+1) = 2*7")
print("  b_1 + 1 = a_1 + 2 = 7 = E8 EXPONENT!")
print()
print("  So: down_shift(g=2) = 2*(b_1+1) = 2*(a_1+2)")
print("  And a_1+2 = 7 is an E8 exponent iff a_1+2 is coprime to h=a_1*b_1")
print("  gcd(7, 30) = 1. YES!")
print()
print("  >>> DOWN QUARK SHIFTS ARE 2*(E8-related quantities) <<<")
print("  shift(1) = 2*N_gen = 2*3 = b_1")
print("  shift(2) = 2*(a_1+2) = 2*7 (7 = E8 exponent)")

# =====================================================================
# PART 5: Why the factor 2? From Galois structure
# =====================================================================
print()
print("=" * 70)
print("PART 5: Origin of the factor 2")
print("=" * 70)
print()

# The factor 2 appears because:
# n_base(g) + delta_down(g) - base_down
# delta_down(g) = F(3) + N_gen*sigma(g) - g*(g+1)
# base_down = delta_down(0) = F(3) + N_gen*sigma(0) - 0 = F(3) + N_gen = 2 + 3 = 5

# So shift = n_base(g) + delta_down(g) - delta_down(0)
# = n_base(g) + N_gen*(sigma(g) - sigma(0)) - g*(g+1)
# For g=1: = (a_1+b_1) + N_gen*(0-1) - 2 = 11 - 3 - 2 = 6
# For g=2: = (a_1+2*b_1) + N_gen*(2-1) - 6 = 17 + 3 - 6 = 14

print("  shift(g) = n_base(g) + N_gen*(sigma(g)-1) - g*(g+1)")
print("  g=1: (a_1+b_1) + N_gen*(0-1) - 2 = 11 - 3 - 2 = 6")
print("  g=2: (a_1+2b_1) + N_gen*(2-1) - 6 = 17 + 3 - 6 = 14")
print()

# Factor 2: can we rewrite?
# shift(1) = a_1 + b_1 - N_gen - 2 = 5 + 6 - 3 - 2 = 6 = 2*3
# shift(2) = a_1 + 2b_1 + N_gen - 6 = 5 + 12 + 3 - 6 = 14 = 2*7

# For a_1 = 5:
# shift(1) = 2*(a_1+b_1-N_gen-2)/2... no, just 6 = 2*3
# 6 is EVEN because a_1+b_1-N_gen-2 = 5+6-3-2 = 6 is even.
# Always even? a_1+b_1 = 2a_1+1 (odd), N_gen = 3 (odd), 2 (even)
# odd - odd - even = even. YES, always even!

print("  shift(1) = (2a_1+1) - N_gen - F(3) = (2a_1+1) - 3 - 2 = 2a_1-4 = 2(a_1-2)")
print("  For a_1=5: 2*(5-2) = 2*3 = 6. CHECK!")
print()
print("  shift(2) = (a_1+2b_1) + N_gen - g*(g+1)")
print("  = a_1 + 2(a_1+1) + 3 - 6 = 3a_1 + 2 + 3 - 6 = 3a_1 - 1")
print("  For a_1=5: 15-1 = 14 = 2*7. CHECK!")
print()
print("  ALWAYS EVEN?")
print("  shift(1) = 2(a_1-2): always even for integer a_1. YES!")
print("  shift(2) = 3a_1-1: even iff a_1 is odd. a_1=5 is odd. YES!")
print()
print("  So the factor 2 in down_shift(1) is ALGEBRAIC (always there)")
print("  And for a_1=5: shift(1)/2 = 3 = N_gen, shift(2)/2 = 7 = E8 exp")

# =====================================================================
# PART 6: Mass corrections from E8 exponent structure
# =====================================================================
print()
print("=" * 70)
print("PART 6: Can E8 exponents explain mass errors?")
print("=" * 70)
print()

print("  Bare formula: m = m_e * phi^n")
print("  Hypothesis: correction depends on distance to nearest E8 exponent")
print()

# For each fermion: distance to nearest E8 exponent
print("  %-5s  n  nearest_E8  dist  m_bare(MeV)  m_exp(MeV)  err%%  sign" % "name")
for name, n, sec, g in fermions:
    dists = [(abs(n-e), e) for e in e8_exponents]
    dists.sort()
    nearest_e8 = dists[0][1]
    dist = n - nearest_e8  # signed distance
    m_bare = m_e * PHI**n
    m_ex = m_exp[name]
    err = (m_bare - m_ex) / m_ex * 100
    print("  %-5s %2d  %10d  %+3d  %11.2f  %10.2f  %+6.1f  %s" %
          (name, n, nearest_e8, dist, m_bare, m_ex, err,
           "E8" if dist == 0 else ""))

print()
print("  Fermions ON E8 exponents: mu, tau, s, b")
print("  Errors for these: mu(-3.8%), tau(+2.7%), s(+8.9%), b(+14.3%)")
print("  NOT systematically better! Being on an E8 exponent doesn't reduce error.")
print()

# =====================================================================
# PART 7: Mass ratios and E8 exponent differences
# =====================================================================
print()
print("=" * 70)
print("PART 7: Mass ratios from generation shifts")
print("=" * 70)
print()

print("  Generation mass ratios (bare formula):")
print("  m(g+1)/m(g) = phi^(n(g+1)-n(g)) = phi^spacing")
print()

for sec_name in ['lep', 'up', 'dn']:
    sec_fermions = [(name,n,g) for name,n,s,g in fermions if s == sec_name]
    print("  %s sector:" % sec_name.upper())
    for i in range(len(sec_fermions)-1):
        name1, n1, g1 = sec_fermions[i]
        name2, n2, g2 = sec_fermions[i+1]
        spacing = n2 - n1
        ratio_pred = PHI**spacing
        ratio_exp = m_exp[name2] / m_exp[name1]
        err = (ratio_pred - ratio_exp) / ratio_exp * 100
        print("    %s->%s: spacing=%d, phi^%d=%.2f, m_exp_ratio=%.2f, err=%+.1f%%" %
              (name1, name2, spacing, spacing, ratio_pred, ratio_exp, err))
    print()

# =====================================================================
# PART 8: The 600-cell spectrum as mass operator
# =====================================================================
print()
print("=" * 70)
print("PART 8: Spectral interpretation")
print("=" * 70)
print()

# 9 eigenvalues of 600-cell correspond to 9 McKay nodes
# Can we map them to 9 fermions (e,mu,tau,u,c,t,d,s,b)?
# Problem: mu and s both have n=11

# The McKay nodes have dims {1,2,3,4,5,6,4,2,3}
# And Coxeter labels {1,2,3,4,5,6,4,2,3}

# Natural: 9 nodes <-> 9 fermions
# But which mapping?

# Approach: sort by Coxeter label (= irrep dimension)
print("  McKay nodes sorted by Coxeter label:")
nodes = [
    (0, 1, 12.0,       "trivial"),
    (1, 2, 6*PHI,      "dim-2 (phi)"),
    (7, 2, 6-6*PHI,    "dim-2 (phi')"),
    (2, 3, 4*PHI,      "dim-3 (phi)"),
    (8, 3, 4-4*PHI,    "dim-3 (phi')"),
    (3, 4, 3.0,        "dim-4 (rational)"),
    (6, 4, -3.0,       "dim-4 (rational)"),
    (4, 5, 0.0,        "dim-5 (zero)"),
    (5, 6, -2.0,       "dim-6 (rational)"),
]

print("  Node  dim  lambda          type")
for node, dim, lam, typ in nodes:
    print("  %4d  %3d  %8.4f  %s" % (node, dim, lam, typ))

print()
print("  9 nodes, 9 fermions. Possible bijection?")
print()

# Attempt: by generation and sector
# Generation 0 (lightest): e(0), u(3), d(5) -> use smallest dims?
# Generation 1 (middle): mu(11), c(16), s(11) -> middle dims
# Generation 2 (heaviest): tau(17), t(26), b(19) -> largest dims

# By dimension:
# dim 1: only one -> must be unique particle
# dim 2: two nodes (Galois pair) -> two particles
# dim 3: two nodes (Galois pair) -> two particles
# dim 4: two nodes (same rational?) -> two particles
# dim 5: only one -> must be unique particle
# dim 6: only one -> must be unique particle

# So: {1} + {2,2} + {3,3} + {4,4} + {5} + {6} = 1+2+2+2+1+1 = 9

# Natural assignment respecting Galois pairs:
# dim 1 (trivial): electron (n=0, ground state, unique)
# dim 6 (largest): top? (n=26, heaviest, unique)
# dim 5 (second): tau? (n=17, but up sector...)

# Actually, let me try by n = 5a + 6b and dim:
# dim 1: n=0 (e), dim^2-1=0 MATCH
# dim 2: n=3 (u), dim^2-1=3 MATCH!
# dim 3: n=8 (rank), dim^2-1=8 (not a fermion n)
# dim 4: n=15, dim^2-1=15 (not a fermion n)
# dim 5: n=24, dim^2-1=24 (not a fermion n)
# dim 6: n=35, dim^2-1=35 = n_nu (neutrino!)

print("  dim^2-1 mapping:")
for dim in sorted(set(d for _,d,_,_ in nodes)):
    val = dim**2 - 1
    matches = [f[0] for f in fermions if f[1] == val]
    print("    dim=%d: dim^2-1=%d  %s" % (dim, val, "= n(%s)" % ",".join(matches) if matches else ""))

# =====================================================================
# PART 9: Grand identity table
# =====================================================================
print()
print("=" * 70)
print("PART 9: Grand identity table")
print("=" * 70)
print()

print("  ALL identities involving n-values that hold for a_1=5:")
print()

# 1. Column sums
print("  1. Column sum g=0->1 = h = 30")
print("     3*(a_1+b_1) + N_gen*(sigma(1)-sigma(0)) = 6*a_1 = a_1*(a_1+1)")
print("     True iff a_1 = 5.")
print()

# 2. Column sum g=1->2
print("  2. Column sum g=1->2 = |2T| = 24")
print("     3*b_1 + N_gen*(sigma(2)-sigma(1)) = 3*a_1 + 9 = 24")
print("     True iff a_1 = 5.")
print()

# 3. det(N)
print("  3. det(N) = 4 = a_1-1")
print("     -a_1^3 + 3*a_1^2 + 10*a_1 + 4 = 4")
print("     (a_1-5)(a_1+2) = 0. True iff a_1 = 5.")
print()

# 4. Trace
print("  4. Tr(N) = 35 = b_1^2 - 1 = |eta|")
print("     n_e + n_c + n_b = 0 + 16 + 19 = 35")
print()

# 5. Column g=0 = rank
print("  5. Column g=0 = 0 + 3 + 5 = 8 = rank(E8)")
print("     base_lep + base_up + base_dn = 0 + (a_1-2) + a_1 = 2*a_1-2 = 8")
print("     True iff a_1 = 5.")
print()

# 6. h * |2T| = edges
print("  6. h * |2T| = 30 * 24 = 720 = edges of 600-cell")
print("     Follows from a_1! = 4*a_1*(a_1+1). True iff a_1 = 5.")
print()

# 7. Row sums
print("  7. Row sums = {28, 45, 35}")
print("     28 = h - 2, 45 = a_1*N_eig, 35 = |eta|")
print()

# 8. Char poly
print("  8. Char poly: x^3 - 35*x^2 - 100*x - 4")
print("     = x^3 - |eta|*x^2 - (a_1-1)*a_1^2*x - (a_1-1)")
print("     p(-1) = 60 = N/2")
print()

# 9. 2x2 minors
print("  9. Spacing matrix 2x2 minors: 2^a_1, 4*13, 4*11")
print("     Sum = 128 = 2^(a_1+2)")
print()

# 10. E8 exponent usage
print(" 10. E8 exponents used: {11,13,17,19,23} (5/8)")
print("     Sum = 83. Unused: {1,7,29}, sum = 37 (prime)")
print("     Used sum + unused sum = 120 = N")
all_e8_sum = sum(e8_exponents)
used = [11, 13, 17, 19, 23]
unused = [1, 7, 29]
print("     Total: %d + %d = %d = N = a_1!" % (sum(used), sum(unused), all_e8_sum))
print()

# 11. n_u + n_d = rank
print(" 11. Sum rules: n_u+n_d=8(rank), n_tau+n_b=36(b_1^2),")
print("     n_t+n_b=45(a_1*N_eig), n_tau+n_u=20(4*a_1)")
print()

# =====================================================================
# PART 10: The UNUSED E8 exponents {1, 7, 29}
# =====================================================================
print()
print("=" * 70)
print("PART 10: The unused E8 exponents {1, 7, 29}")
print("=" * 70)
print()

print("  Sum of used E8 exps: 11+13+17+19+23 = %d" % sum(used))
print("  Sum of unused E8 exps: 1+7+29 = %d" % sum(unused))
print("  Sum of ALL E8 exps: %d = N = a_1! = 120" % all_e8_sum)
print()
print("  This is REMARKABLE: sum of all E8 exponents = N = |600-cell vertices|")
print()

# Verify: E8 exponents are {m : 1 <= m <= h-1, gcd(m,h) = 1}
# For h=30: coprime to 30 means coprime to 2, 3, 5
# i.e., odd, not divisible by 3, not divisible by 5
# = {1,7,11,13,17,19,23,29}
# Sum = 120. This is Euler's totient: sum of m coprime to h.
# Actually: sum = h * phi(h) / 2 where phi = Euler totient
# phi(30) = 30 * (1-1/2) * (1-1/3) * (1-1/5) = 30 * 1/2 * 2/3 * 4/5 = 8
# So sum = 30 * 8 / 2 = 120. YES!

print("  WHY sum = 120?")
print("  E8 exponents = {m : 1 <= m <= h-1, gcd(m,h) = 1}")
print("  Sum = h * phi(h) / 2 where phi = Euler totient")
print("  phi(30) = 30 * (1-1/2)(1-1/3)(1-1/5) = 8 = rank(E8)")
print("  Sum = 30 * 8 / 2 = 120 = N")
print()
print("  >>> rank(E8) = phi(h(E8)) = Euler totient of Coxeter number <<<")
print("  >>> sum(E8 exponents) = N = a_1! <<<")
print("  Both are STANDARD FACTS in Coxeter theory!")
print()

# So the unused exponents satisfy:
# sum(unused) = N - sum(used) = 120 - 83 = 37
# 37 is prime. And 37 = h + 7 = h + E8_exp(2)
# Also: 37 = a_1^2 + a_1 + 7 = 25 + 5 + 7. Not clean.

print("  Unused: {1, 7, 29}")
print("  1 is the trivial exponent (identity of Coxeter group)")
print("  29 = h-1 is its Coxeter partner")
print("  7 is the partner of 23 (which IS used)")
print()
print("  INTERPRETATION:")
print("  (1, 29): ground state pair. Defines the RANGE of mass spectrum.")
print("    n_min = 0 (slightly below 1), n_max = 26 (slightly below 29)")
print("    The mass spectrum fits INSIDE the (1,29) pair.")
print()
print("  7: appears as half of down_shift(g=2) = 2*7 = 14")
print("    So 7 IS used, but through DOUBLING.")
print("    While 23 (its partner) is used DIRECTLY.")
print("    Asymmetry: up quarks use E8 exps directly,")
print("    down quarks use them doubled.")

# =====================================================================
# PART 11: Summary and interpretation
# =====================================================================
print()
print("=" * 70)
print("GRAND SYNTHESIS SUMMARY")
print("=" * 70)
print()

print("THE MASS EXPONENT STRUCTURE (complete picture):")
print()
print("1. THREE BUILDING BLOCKS:")
print("   a) Generation base {0, 11, 17} from Fibonacci + Z[phi] units")
print("   b) Sector offsets (0, 3, 5) summing to rank(E8) = 8")
print("   c) Casimir correction g*(g+1) from generation SU(2)")
print()

print("2. E8 EXPONENT ASSIGNMENT:")
print("   - 4 Coxeter pairs serve different sector-generation combinations")
print("   - Pair (11,19): muon/strange + bottom (lepton-down link)")
print("   - Pair (13,17): tau + charm (lepton-up link)")
print("   - Pair (7,23): top + bottom via doubling (up-down link)")
print("   - Pair (1,29): vacuum pair, brackets the spectrum")
print()

print("3. COLUMN SUM RULES (THEOREM):")
print("   - g=0->1 spacings sum to h(E8) = 30")
print("   - g=1->2 spacings sum to |2T| = 24")
print("   - Product: h*|2T| = 720 = edges (from Diophantine)")
print("   - Both are algebraic identities for a_1 = 5")
print()

print("4. MATRIX IDENTITIES (THEOREM):")
print("   - det(N) = 4 = a_1-1 (identity for a_1=5)")
print("   - Tr(N) = 35 = |eta| = spectral asymmetry")
print("   - Char poly coefficients = framework constants")
print("   - p(-1) = 60 = N/2")
print()

print("5. DOWN QUARK DOUBLING:")
print("   - Down shifts = 2*{N_gen, E8_exp(7)} = {6, 14}")
print("   - Factor 2 is algebraic (shift is always even)")
print("   - Connects down quarks to up quarks via Coxeter partners")
print()

print("6. OPEN QUESTIONS:")
print("   - WHY does pair (1,29) not contribute directly?")
print("   - Can E8 exponents explain mass errors (3-21%)?")
print("     Answer: NO (being on E8 exp doesn't improve accuracy)")
print("   - The mass formula m=m_e*phi^n remains BARE prediction")
print("   - Corrections, if any, come from dynamics (QCD, QED running)")
print()

print("CATEGORY: ALL column sum rules and det(N) are DERIVAT.")
print("Coxeter pair assignment is PATTERN (derived structure, not mechanism).")
print("Mass corrections remain OPEN.")
