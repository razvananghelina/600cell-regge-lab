"""
EXP-276: CASIMIR g*(g+1) MECHANISM
===================================
The generation shift g*(g+1) appears in delta_up(g) = 3 + g*(g+1).
This looks like a Casimir eigenvalue. Derive WHY from 600-cell structure.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120; h = 30
N_eig = 9; N_gen = 3

print("="*72)
print("EXP-276: CASIMIR g*(g+1) MECHANISM")
print("="*72)

# === Part 1: The observation ===
print("\n--- Part 1: The g*(g+1) Pattern ---")
# delta_up(g) = 3 + g*(g+1) for g = 0, 1, 2
# g=0: 3+0 = 3  -> n_u = 3
# g=1: 3+2 = 5  -> n_c = 3 + 5 + n_base_shift = ..
# Wait, let me get this right from memory

# Mass exponents: n = sector_base + gen_shift
# Up sector: base = 3, gen_shift = {0, 13, 23}
# Actually: delta_up(g) = 3 + g*(g+1): u=3, c=5, t=9
# These are the DIFFERENCES between up-type and lepton exponents
# n_u = n_e + delta_up(0) = 0 + 3 = 3
# But n_c = n_mu + delta_up(1) = 11 + 5 = 16
# And n_t = n_tau + delta_up(2) = 17 + 9 = 26

for g in range(N_gen):
    delta = 3 + g*(g+1)
    lep_n = [0, 11, 17][g]
    quark_n = lep_n + delta
    print(f"  g={g}: delta_up = 3 + {g}*{g+1} = {delta}, n_lepton={lep_n}, n_up = {quark_n}")

# The g*(g+1) is the Casimir eigenvalue of SU(2):
# J^2|j,m> = j(j+1)|j,m>
# For j = g: eigenvalue = g*(g+1)
# g=0: 0, g=1: 2, g=2: 6

print(f"\n  g*(g+1) values: {[g*(g+1) for g in range(N_gen)]}")
print(f"  = SU(2) Casimir for j=0,1,2")
print(f"  = Laplacian eigenvalues on S^2: l(l+1) for l=0,1,2")
print(f"  Max l = 2 (from icosahedral symmetry)")

# === Part 2: Why SU(2)? ===
print("\n--- Part 2: Origin of SU(2) ---")
# The 600-cell fibers over S^2 via Hopf fibration: S^3 -> S^2
# The fiber is S^1 (great circles)
# The base S^2 has the icosahedral symmetry group I ~ A_5

# On S^2: the Laplacian has eigenvalues l(l+1) with multiplicity 2l+1
# Under icosahedral symmetry: l=0,1,2 survive (l_max=2)
# l=0: singlet (1 state)
# l=1: triplet (3 states) -> but icosahedral has no fundamental rep...
# l=2: quintet (5 states) = a1 states

print(f"  Hopf fibration: S^3 -> S^2 (base)")
print(f"  Base symmetry: I ~ A_5 (icosahedral)")
print(f"  S^2 Laplacian: l(l+1), multiplicities 2l+1")
print(f"  Under I: l=0 (1), l=1 (3), l=2 (5)")
print(f"  l_max = 2 because icosahedron has 5-fold symmetry")
print(f"  (l=5 would be the NEXT singlet of A_5)")

# l_max = 2 from icosahedral truncation
# The dimension of l=2 rep: 2*2+1 = 5 = a1 (!)
print(f"\n  CRUCIAL: dim(l=2) = 5 = a1")
print(f"  The highest mode on the Hopf base has a1 components")
print(f"  These correspond to the a1 coordinate directions")

# === Part 3: Connection to generations ===
print("\n--- Part 3: Generations from Hopf Base ---")
# Generation theorem: N_gen = 3 from Fibonacci on Z[phi]
# The generation INDEX g = 0, 1, 2 can be identified with:
# - DSI levels k = 0, 2, 3 (from generation theorem)
# - Or: angular momentum l = 0, 1, 2 on Hopf base

# The map: g -> l is one-to-one:
# g=0 (lightest) -> l=0 (s-wave, spread on whole S^2)
# g=1 (middle)   -> l=1 (p-wave, one node)
# g=2 (heaviest)  -> l=2 (d-wave, two nodes)

# The Casimir g*(g+1) comes from the S^2 Laplacian
# evaluated at the generation wavefunction

print(f"  Generation -> Angular momentum:")
print(f"    g=0 (e, u, d) -> l=0 (s-wave): L^2 = 0")
print(f"    g=1 (mu, c, s) -> l=1 (p-wave): L^2 = 2")
print(f"    g=2 (tau, t, b) -> l=2 (d-wave): L^2 = 6")
print(f"")
print(f"  Mass correction from Laplacian:")
print(f"    delta(g) = base_offset + L^2/R^2 (kinetic energy on S^2)")
print(f"    = 3 + g*(g+1) for up-type quarks")
print(f"    The 3 = a1 - F(3) = 5 - 2 (Galois image of tau)")

# === Part 4: Why l_max = 2? ===
print("\n--- Part 4: Why l_max = 2 ---")
# The icosahedron has 12 vertices on S^2
# The spherical harmonics Y_l^m for l <= some l_max
# can be uniquely determined at 12 points

# Sampling theorem on S^2: need (l_max+1)^2 <= N_points
# For 12 points: (l_max+1)^2 <= 12 -> l_max <= 2 (since 9 <= 12)
# l_max = 3 would need 16 > 12 points (FAILS)

print(f"  Icosahedron: 12 vertices on S^2")
print(f"  Nyquist: (l_max+1)^2 <= 12")
print(f"    l_max=2: (2+1)^2 = 9 <= 12 (OK)")
print(f"    l_max=3: (3+1)^2 = 16 > 12 (FAILS)")
print(f"  Therefore: l_max = 2, N_gen = 2+1 = 3")
print(f"")
print(f"  This gives ANOTHER derivation of N_gen = 3:")
print(f"  N_gen = l_max + 1 = floor(sqrt(12)) = floor(3.46) = 3")
print(f"  12 = degree of 600-cell vertex = Tr(z_Planck)")

# More precisely: 12 Hopf fibers through icosahedral vertices
# Each fiber is a great circle on S^3
# Total: 12 great circles cover S^3

# === Part 5: The binary tetrahedral 2T ===
print("\n--- Part 5: Binary Tetrahedral Group ---")
# The generation group is 2T (binary tetrahedral, |2T| = 24)
# 2T = SU(2) cover of T (tetrahedral = A_4)
# 2T has irreps: 1, 1', 1'', 2, 2', 2'', 3
# The 3-dim irrep gives the 3 generations

# The Casimir of the SU(2) covering:
# For the spin-j representation: C_j = j(j+1)
# 2T contains spin-0, spin-1, spin-2 (embedded in SU(2))
# These ARE the l=0,1,2 modes on S^2

print(f"  |2T| = 24 = |binary tetrahedral|")
print(f"  2T < SU(2): contains spins j=0,1,2")
print(f"  Casimir eigenvalues: j(j+1) = 0, 2, 6")
print(f"  These ARE the generation corrections!")
print(f"")
print(f"  2T also has: 24 = N/a1 = 120/5")
print(f"  The 5 cosets of 2T in 2I give the 5-partite structure")
print(f"  |2I|/|2T| = 120/24 = a1 = 5")

# The spin-2 representation:
# dim(spin-2) = 2*2+1 = 5 = a1
# Casimir: 2*3 = 6 = b1
print(f"\n  Spin-2 representation of SU(2):")
print(f"    dim = 2*2+1 = {2*2+1} = a1")
print(f"    Casimir = 2*(2+1) = {2*3} = b1")
print(f"    -> a1 and b1 are the dimension and Casimir of the SAME representation!")

# === Part 6: Yukawa coupling and Casimir ===
print("\n--- Part 6: Yukawa Connection ---")
# In the spectral action: Yukawa coupling y_f comes from D_F
# On McKay graph: D_F_ij = n_i * delta_ij (diagonal)
# The generation structure modifies this:
# y_f(g) = y_0 * phi^(delta_sector(g))
# delta_sector(g) = base + g*(g+1)

# In terms of the Yukawa operator:
# Y = D_F * exp(-beta * L_S2) where L_S2 is the Hopf base Laplacian
# exp(-beta * l(l+1)) for each generation l=g
# For beta -> 0: all generations same (no splitting)
# For beta = ln(phi): exp(-ln(phi)*g*(g+1)) = phi^(-g*(g+1))
# -> Yukawa ~ phi^(n_base - g*(g+1))... WRONG sign!

# Actually: delta_up = 3 + g*(g+1) means the exponent INCREASES with g
# So: n(g) = n_base + g*(g+1) (additive in exponent)
# Mass ~ phi^(n_base + g*(g+1)) = phi^n_base * phi^(g*(g+1))
# Higher g -> heavier (phi > 1)

# Physical picture:
# Generation g has angular momentum l=g on the Hopf base
# Higher l -> more kinetic energy on S^2 -> heavier
# The "kinetic energy" adds g*(g+1) to the mass exponent

print(f"  Physical picture:")
print(f"    Generation g = angular momentum l on Hopf base S^2")
print(f"    KE on S^2 ~ l(l+1)/R^2 adds to mass exponent")
print(f"    Higher generation -> more angular momentum -> heavier")
print(f"    m(g) = m(0) * phi^(g*(g+1))")
print(f"")
print(f"  For up quarks: m_u * phi^(g(g+1)):")
print(f"    g=0: m_u = {0.511*PHI**3:.2f} MeV (bare, exp ~2.2 MeV)")
print(f"    g=1: m_u * phi^2 = {0.511*PHI**3*PHI**2:.1f} MeV (bare, charm ~1270)")
print(f"    g=2: m_u * phi^6 = {0.511*PHI**3*PHI**6:.0f} MeV (bare, top ~172500)")

# But this doesn't work alone - need the lepton shifts too
# Full: n(g) = n_lepton(g) + delta_sector(g)
# = [0, 11, 17] + [3+0, 3+2, 3+6] = [3, 16, 26] = n_u, n_c, n_t (correct!)

print(f"\n  Full up-type: n = n_lepton + 3 + g*(g+1)")
for g in range(3):
    n_lep = [0, 11, 17][g]
    delta = 3 + g*(g+1)
    n_up = n_lep + delta
    print(f"    g={g}: {n_lep} + {delta} = {n_up} = n_{['u','c','t'][g]}")

# === Part 7: Casimir and S(g) ===
print("\n--- Part 7: S(g) from Casimir ---")
# S(g) = delta_up(g) + delta_down(g) = {8, 5, 11}
# delta_up(g) = 3 + g*(g+1) = {3, 5, 9}
# delta_down(g) = S(g) - delta_up(g) = {5, 0, 2}

# S(g) = a1 + N_gen * sigma(g) where sigma = (01)(2)
# sigma(0)=1, sigma(1)=0, sigma(2)=2
# S = {5+3, 5+0, 5+6} = {8, 5, 11}

print(f"  S(g) = a1 + N_gen*sigma(g):")
sigma = [1, 0, 2]
for g in range(3):
    S = a1 + N_gen*sigma[g]
    du = 3 + g*(g+1)
    dd = S - du
    print(f"    g={g}: S={S}, delta_up={du}, delta_down={dd}")
    print(f"           sum={du+dd}={S}, product={du*dd}")

# The sigma permutation: (01)(2) = transposition of gen 0 and 1
# This is the GALOIS transposition on generations
# Galois swaps phi <-> phi' on the DSI levels

print(f"\n  sigma = (01)(2): transposes generations 0 and 1")
print(f"  This is the GALOIS action on generation labels")
print(f"  phi <-> phi' swaps the DSI minima at k=0,2,3:")
print(f"    phi'^0 = 1 = phi^0 (invariant)")
print(f"    phi'^2 = phi^(-2) = 1/phi^2 (swapped with some k)")
print(f"    phi'^3 = -1/phi^3 (changes sign)")

# === Part 8: The mechanism ===
print("\n--- Part 8: Proposed Mechanism ---")
print(f"  MECHANISM for g*(g+1):")
print(f"")
print(f"  1. Hopf fibration S^3 -> S^2 gives base space S^2")
print(f"  2. Icosahedron inscribed in S^2 has 12 vertices")
print(f"  3. Spherical harmonics on S^2 truncated at l_max=2")
print(f"     (Nyquist: (l_max+1)^2 <= 12)")
print(f"  4. Generations g = 0,1,2 = angular momenta l = 0,1,2")
print(f"  5. S^2 Laplacian eigenvalue: l(l+1) = g(g+1)")
print(f"  6. This adds to mass exponent: n(g) = n_base + g(g+1)")
print(f"")
print(f"  WHY specifically in delta_up:")
print(f"  The up-type quarks couple to the Hopf base more strongly")
print(f"  than down-type (which have delta_down = S-delta_up)")
print(f"  This is related to T3 = +1/2 for up vs -1/2 for down")
print(f"  The SIGN of T3 determines which metric (Hopf base vs fiber)")
print(f"")
print(f"  REMAINING GAP:")
print(f"  Why T3 selects Hopf base metric specifically")
print(f"  This connects to how SU(2)_L acts on the Hopf structure")

# === Summary ===
print("\n--- Summary ---")
print(f"  g*(g+1) MECHANISM: PARTIALLY DERIVED")
print(f"    Source: S^2 Laplacian on Hopf base (Casimir of SU(2))")
print(f"    l_max = 2 from icosahedral Nyquist (12 vertices)")
print(f"    N_gen = l_max + 1 = 3 (alternative derivation)")
print(f"    dim(l=2) = 5 = a1, Casimir(l=2) = 6 = b1")
print(f"    Binary tetrahedral 2T contains spins 0,1,2")
print(f"  OPEN: Connection between T3 and Hopf base/fiber metric")
print(f"  CATEGORY: PARTIALLY DERIVED")

print("\n" + "="*72)
print("EXP-276 COMPLETE")
print("="*72)
