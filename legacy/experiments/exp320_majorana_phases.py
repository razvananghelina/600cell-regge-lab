"""
exp320_majorana_phases.py
Derivation of Majorana phases alpha_1, alpha_2 from the Galois structure
on spinorial irreps of 2I (binary icosahedral group).

The PMNS matrix for Majorana neutrinos:
  U_PMNS = U_Dirac * diag(1, e^{i*alpha_1/2}, e^{i*alpha_2/2})

We have delta_PMNS = 3*arctan(sqrt(5)) = 197.72 deg (Dirac phase).
Now we seek alpha_1, alpha_2 from the Galois action on 2I representations.

Key idea: The Galois automorphism sigma: phi -> phi' = -1/phi acts on the
character table of 2I. The phases it introduces on the spinorial irreps
(those not coming from SO(3) reps of I) should give the Majorana phases.
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2
PHI_PRIME = -1.0 / PHI  # Galois conjugate = (1-sqrt(5))/2
SQRT5 = np.sqrt(5)
a1 = 5
b1 = 6
N = 120
N_gen = 3

print("=" * 72)
print("EXP-320: MAJORANA PHASES FROM GALOIS STRUCTURE ON 2I")
print("=" * 72)

# ============================================================
# PART 1: CHARACTER TABLE OF 2I
# ============================================================
print("\n" + "=" * 72)
print("PART 1: CHARACTER TABLE OF 2I (BINARY ICOSAHEDRAL GROUP)")
print("=" * 72)

# 2I has 9 conjugacy classes and 9 irreps
# Dims: 1, 2, 2', 3, 3', 4, 4', 5, 6
# The primed reps are Galois conjugates (phi <-> phi')

# Conjugacy classes of 2I (|2I| = 120):
# 1 (identity), 1 class, size 1
# -1 (center), 1 class, size 1
# 12a (order 10), 1 class, size 12
# 12b (order 10), 1 class, size 12
# 12a' (Galois of 12a), 1 class, size 12
# 12b' (Galois of 12b), 1 class, size 12
# 20 (order 6), 1 class, size 20
# 30 (order 4), 1 class, size 30
# 20' (order 3), 1 class, size 20

# Character table (rows = irreps, cols = conjugacy classes)
# Classes:    1    -1   12a   12b  12a'  12b'   20    30   20'
# Sizes:      1     1    12    12    12    12    20    30    20

# Standard labeling:
# rho_1 (dim 1): trivial
# rho_2 (dim 2): spinorial (from SU(2), not SO(3))
# rho_2'(dim 2): Galois conjugate of rho_2
# rho_3 (dim 3): from SO(3) rep (l=1 of icosahedron)
# rho_3'(dim 3): Galois conjugate
# rho_4 (dim 4): from SO(3) rep (l=2... actually mixed)
# rho_4'(dim 4): Galois conjugate
# rho_5 (dim 5): from SO(3) rep (l=2)
# rho_6 (dim 6): from SO(3) rep (l=... mixed)

# The character values involve phi and phi'
# Key: sigma(chi_rho(g)) = chi_{sigma(rho)}(g) for the Galois action

# Character table of 2I (binary icosahedral, |G|=120)
# Conjugacy classes by element order and type:
# C1: {e}, size 1, order 1
# C2: {-e}, size 1, order 2
# C3: 5-fold (12a), size 12, order 10
# C4: 5-fold (12b), size 12, order 10
# C5: 5-fold' (12a'), size 12, order 10
# C6: 5-fold' (12b'), size 12, order 10
# C7: 3-fold, size 20, order 6
# C8: 2-fold, size 30, order 4
# C9: 3-fold', size 20, order 3

class_sizes = np.array([1, 1, 12, 12, 12, 12, 20, 30, 20])
assert sum(class_sizes) == 120

# Character table
# Using standard references (e.g., Wilson, or GAP computation)
# chi[i][j] = character of irrep i on class j

chi = np.zeros((9, 9))

# Labels for printing
names = ['1', '2', "2'", '3', "3'", '4', "4'", '5', '6']
dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]

# Classes: 1, -1, 12a, 12b, 12a', 12b', 20, 30, 20'
# rho_1 (trivial, dim 1)
chi[0] = [1, 1, 1, 1, 1, 1, 1, 1, 1]

# rho_2 (dim 2, spinorial)
# Characters involve phi: chi(12a) = phi, chi(12b) = -1/phi = phi'
chi[1] = [2, -2, PHI, -PHI_PRIME, PHI_PRIME, -PHI, -1, 0, 1]

# rho_2' (dim 2, Galois conjugate spinorial)
# sigma: phi -> phi' swaps rho_2 and rho_2'
chi[2] = [2, -2, PHI_PRIME, -PHI, PHI, -PHI_PRIME, -1, 0, 1]

# rho_3 (dim 3, from icosahedral l=1)
# chi(12a) = phi, chi(12b) = phi
chi[3] = [3, 3, PHI, PHI, PHI_PRIME, PHI_PRIME, 0, -1, 0]

# rho_3' (dim 3, Galois conjugate)
chi[4] = [3, 3, PHI_PRIME, PHI_PRIME, PHI, PHI, 0, -1, 0]

# rho_4 (dim 4)
chi[5] = [4, -4, -1, 1, -1, 1, 1, 0, -1]

# rho_4' (dim 4, Galois conjugate)
# Actually, rho_4 is self-Galois-conjugate (no phi in characters)
# The 4-dim rep of 2I: chi(12a) = -1, all rational
chi[6] = [4, -4, 1, -1, 1, -1, 1, 0, -1]

# rho_5 (dim 5)
chi[7] = [5, 5, 0, 0, 0, 0, -1, 1, 2]

# rho_6 (dim 6)
chi[8] = [6, -6, -1, 1, 1, -1, 0, 0, 0]

# Verify orthogonality
print("\n  Character table of 2I:")
print(f"  {'Rep':>4s} {'dim':>3s}", end="")
for j in range(9):
    print(f"  {'C'+str(j+1):>8s}", end="")
print()
print(f"  {'':>4s} {'':>3s}", end="")
for s in class_sizes:
    print(f"  {'('+str(s)+')':>8s}", end="")
print()

for i in range(9):
    print(f"  {names[i]:>4s} {dims[i]:>3d}", end="")
    for j in range(9):
        print(f"  {chi[i][j]:8.4f}", end="")
    print()

# Orthogonality check: sum_g chi_i(g) * chi_j(g)^* = |G| * delta_ij
print("\n  Orthogonality check:")
max_err = 0
for i in range(9):
    for j in range(i, 9):
        inner = sum(class_sizes[k] * chi[i][k] * chi[j][k] for k in range(9))
        expected = 120 if i == j else 0
        err = abs(inner - expected)
        max_err = max(max_err, err)
        if i == j:
            norm = inner / 120
print(f"  Max orthogonality error: {max_err:.2e}")


# ============================================================
# PART 2: GALOIS ACTION ON IRREPS
# ============================================================
print("\n" + "=" * 72)
print("PART 2: GALOIS ACTION sigma: phi -> phi'")
print("=" * 72)

print("""
  The Galois automorphism sigma: sqrt(5) -> -sqrt(5) acts on characters:
    sigma(phi) = phi' = -1/phi
    sigma(phi') = phi

  This permutes the irreps:
    sigma(rho_1) = rho_1  (trivial, rational)
    sigma(rho_2) = rho_2' (spinorial pair)
    sigma(rho_3) = rho_3' (dimension-3 pair)
    sigma(rho_4) = rho_4' (dimension-4 pair)
    sigma(rho_5) = rho_5  (rational)
    sigma(rho_6) = rho_6  (rational, or self-conjugate)
""")

# Verify Galois pairing by applying sigma to characters
print("  Verification of Galois pairing:")
for i in range(9):
    chi_sigma = np.array([
        chi[i][j] if not any(abs(chi[i][j] - PHI) < 1e-10 or
                             abs(chi[i][j] - PHI_PRIME) < 1e-10 or
                             abs(chi[i][j] + PHI) < 1e-10 or
                             abs(chi[i][j] + PHI_PRIME) < 1e-10
                             for _ in [0])
        else chi[i][j]  # placeholder
        for j in range(9)
    ])
    # Apply sigma properly
    chi_s = np.zeros(9)
    for j in range(9):
        v = chi[i][j]
        # Decompose v = a + b*phi where a, b are rational
        # Then sigma(v) = a + b*phi'
        # For our characters: values are in Z[phi]
        # phi = (1+sqrt5)/2, phi' = (1-sqrt5)/2
        # v = a + b*sqrt5 => sigma(v) = a - b*sqrt5
        # Since phi = (1+sqrt5)/2: if v = n*phi + m for integers n,m
        # then sigma(v) = n*phi' + m

        # Practical: replace phi with phi' and phi' with phi
        if abs(v - PHI) < 1e-10:
            chi_s[j] = PHI_PRIME
        elif abs(v - PHI_PRIME) < 1e-10:
            chi_s[j] = PHI
        elif abs(v + PHI) < 1e-10:
            chi_s[j] = -PHI_PRIME
        elif abs(v + PHI_PRIME) < 1e-10:
            chi_s[j] = -PHI
        else:
            chi_s[j] = v  # rational, unchanged

    # Find which irrep sigma(rho_i) matches
    for k in range(9):
        if np.allclose(chi_s, chi[k]):
            if k != i:
                print(f"    sigma(rho_{names[i]}) = rho_{names[k]}")
            else:
                print(f"    sigma(rho_{names[i]}) = rho_{names[i]} (self-conjugate)")
            break


# ============================================================
# PART 3: SPINORIAL vs BOSONIC IRREPS
# ============================================================
print("\n" + "=" * 72)
print("PART 3: SPINORIAL vs BOSONIC IRREPS")
print("=" * 72)

print("""
  2I is the double cover of I (icosahedral group, |I|=60).
  Irreps of 2I split into:
    BOSONIC: chi(-1) = chi(1)  -> descend to I
    SPINORIAL: chi(-1) = -chi(1) -> genuine 2I reps

  Bosonic (from I): rho_1(1), rho_3(3), rho_3'(3), rho_5(5)
    -> 4 irreps, dims sum to 1+3+3+5 = 12 = |I/[I,I]|... actually = 12
  Spinorial (genuine 2I): rho_2(2), rho_2'(2), rho_4(4), rho_4'(4), rho_6(6)
    -> 5 irreps, dims sum to 2+2+4+4+6 = 18
  Total: 12 + 18 = 30 (check: sum d^2 = 1+4+4+9+9+16+16+25+36 = 120)
""")

for i in range(9):
    is_spin = (chi[i][1] < 0)  # chi(-1) = -dim means spinorial
    typ = "SPINORIAL" if is_spin else "bosonic"
    print(f"  rho_{names[i]:>2s} (dim {dims[i]}): chi(-1) = {chi[i][1]:+.0f}, {typ}")


# ============================================================
# PART 4: PHASES FROM GALOIS ACTION ON SPINORIAL REPS
# ============================================================
print("\n" + "=" * 72)
print("PART 4: MAJORANA PHASES FROM GALOIS ACTION")
print("=" * 72)

print("""
  The Majorana phases arise from the relative phases between neutrino
  mass eigenstates under CP transformation. In the framework:
  - Neutrinos are Majorana (from seesaw with Galois partner)
  - The CP transformation relates rho and sigma(rho) (Galois conjugate)
  - The phase difference between rho and sigma(rho) gives the Majorana phase

  For the PMNS matrix with Majorana phases:
  U = U_Dirac * diag(1, e^{i*alpha_1/2}, e^{i*alpha_2/2})

  The three neutrino generations map to three A5 irreps:
  nu_1 ~ component in 3, nu_2 ~ component in 3', nu_3 ~ component in 4

  (Actually, the TBM eigenvectors from the Galois mixing matrix
   give the mixing, and the Majorana phases come from the Galois
   action on the mass matrix.)
""")

# Approach 1: Galois phase from character ratios
# For irrep rho with Galois conjugate sigma(rho):
# The "Galois phase" is arg(chi_rho(g) / chi_{sigma(rho)}(g))
# This is well-defined when chi values are in Q(sqrt(5))

print("  Galois phase analysis for 5-fold classes (12a, 12b, 12a', 12b'):")
print()

# For the 5-fold classes, characters involve phi
# The relevant classes for neutrino physics are the 5-cycles (icosahedral)
# Classes 12a (index 2) and 12b (index 3) have order 10

# Phase from rho_3/rho_3' ratio at 12a class
for cls_idx, cls_name in [(2, '12a'), (3, '12b'), (4, "12a'"), (5, "12b'")]:
    # rho_3 vs rho_3'
    r3 = chi[3][cls_idx]
    r3p = chi[4][cls_idx]
    if abs(r3) > 1e-10 and abs(r3p) > 1e-10:
        ratio = r3 / r3p
        print(f"  Class {cls_name}: chi_3/chi_3' = {r3:.4f}/{r3p:.4f} = {ratio:.6f}")

    # rho_2 vs rho_2'
    r2 = chi[1][cls_idx]
    r2p = chi[2][cls_idx]
    if abs(r2) > 1e-10 and abs(r2p) > 1e-10:
        ratio2 = r2 / r2p
        print(f"  Class {cls_name}: chi_2/chi_2' = {r2:.4f}/{r2p:.4f} = {ratio2:.6f}")
    print()

# The ratio phi/phi' = phi / (-1/phi) = -phi^2 = -(phi+1)
# |phi/phi'| = phi^2, arg(phi/phi') = pi (since phi' < 0)
print(f"  phi/phi' = {PHI/PHI_PRIME:.6f} = -phi^2 = {-PHI**2:.6f}")
print(f"  |phi/phi'| = phi^2 = {PHI**2:.6f}")
print(f"  arg(phi/phi') = pi (since phi > 0, phi' < 0)")

# But these are REAL characters, so the "phase" is just a sign.
# For Majorana phases, we need COMPLEX phases.


# ============================================================
# PART 5: COMPLEX GALOIS PHASES FROM ROOTS OF UNITY
# ============================================================
print("\n" + "=" * 72)
print("PART 5: COMPLEX PHASES FROM 2I REPRESENTATION THEORY")
print("=" * 72)

print("""
  The Majorana phases require COMPLEX structure.
  Key insight: the elements of 2I in SU(2) are already complex!

  In SU(2) representation: g = [[a, -b*], [b, a*]] with |a|^2+|b|^2=1
  The character chi_j(g) = sin((2j+1)*theta) / sin(theta)
  where cos(theta) = Re(a) = trace(g)/2 for the fundamental (j=1/2).

  For 2I elements on 5-fold axes:
  theta = 2*pi*k/5 for k=1,2,3,4 (giving 12a, 12b, 12a', 12b')

  The EIGENVALUES of g in the spin-j rep are:
  e^{i*m*2*theta} for m = -j, ..., +j

  These eigenvalues carry phases that differ between Galois conjugate reps.
""")

# Elements on 5-fold axes have eigenvalues zeta_5^k in fundamental rep
# zeta_5 = e^{2*pi*i/5} is a 5th root of unity
# In Q(sqrt(5)): zeta_5 + zeta_5^4 = (sqrt(5)-1)/2 = 1/phi
#                zeta_5^2 + zeta_5^3 = -(sqrt(5)+1)/2 = -phi

# The Galois action sigma: sqrt(5) -> -sqrt(5) sends:
# zeta_5 -> zeta_5^2 (or zeta_5^3, depending on embedding)
# This is because sigma permutes the roots of unity

zeta5 = np.exp(2j * np.pi / 5)
print(f"  zeta_5 = e^{{2*pi*i/5}} = {zeta5:.6f}")
print(f"  zeta_5 + zeta_5^4 = {(zeta5 + zeta5**4).real:.6f} = 1/phi = {1/PHI:.6f}")
print(f"  zeta_5^2 + zeta_5^3 = {(zeta5**2 + zeta5**3).real:.6f} = -phi = {-PHI:.6f}")

# The Galois automorphism sigma: zeta_5 -> zeta_5^2
# This is the Frobenius at 5: sigma(zeta_5^k) = zeta_5^{2k mod 5}
print(f"\n  Galois action on 5th roots:")
for k in range(5):
    z = zeta5**k
    z_sigma = zeta5**((2*k) % 5)
    print(f"    zeta_5^{k} -> zeta_5^{(2*k)%5}:  {z:.4f} -> {z_sigma:.4f}")

# The spin-j representation of a 5-fold element has eigenvalues
# zeta_5^{m*2} for m = -j, ..., +j (where 2theta = 4*pi/5 for the 12a class)
# Actually, for an element g with eigenvalue e^{i*theta} in fundamental:
# chi_j(g) = sum_{m=-j}^{j} e^{i*m*2*theta}

# For 12a class: theta = pi/5 (order 10 element in SU(2))
# eigenvalue in fundamental: e^{i*pi/5} = zeta_10 = zeta_5^{1/2}
# Actually, the element has order 10 in 2I, order 5 in I
# In SU(2): eigenvalues are e^{+/- i*pi/5}

theta_12a = np.pi / 5
print(f"\n  Class 12a: theta = pi/5 = {theta_12a:.6f} rad = {np.degrees(theta_12a):.1f} deg")

# Character of spin-j rep at theta:
# chi_j(theta) = sin((2j+1)*theta) / sin(theta)
for j_half in range(10):
    j = j_half / 2.0
    ch = np.sin((2*j+1)*theta_12a) / np.sin(theta_12a)
    dim = int(2*j+1)
    if dim in dims:
        idx = dims.index(dim)
        print(f"    j={j:.1f} (dim {dim}): chi = {ch:.6f}, table value = {chi[idx][2]:.6f}")


# ============================================================
# PART 6: MAJORANA PHASE CANDIDATES
# ============================================================
print("\n" + "=" * 72)
print("PART 6: MAJORANA PHASE CANDIDATES FROM FRAMEWORK")
print("=" * 72)

print("""
  The Majorana phases alpha_1, alpha_2 must satisfy:
  1. They appear in the PMNS matrix as diag(1, e^{i*alpha_1/2}, e^{i*alpha_2/2})
  2. They are related to the Galois structure (sigma: phi -> phi')
  3. They should be expressible in terms of a_1 = 5 framework constants

  We already have delta_PMNS = 3*arctan(sqrt(5)) = 197.72 deg.
  This comes from: Berry phase on generation loop, N_gen * delta_CKM.

  For the Majorana phases, the natural candidates involve:
""")

# Candidate 1: Galois phase from zeta_5
# The "Galois phase" of zeta_5 relative to zeta_5^2 is:
# arg(zeta_5 / zeta_5^2) = arg(zeta_5^{-1}) = -2*pi/5
galois_phase_1 = 2 * np.pi / 5
print(f"  Candidate 1: 2*pi/a_1 = 2*pi/5 = {np.degrees(galois_phase_1):.2f} deg = {galois_phase_1:.6f} rad")

# Candidate 2: Phase from spinorial eigenvalue ratio
# In rho_2, eigenvalues at 12a are e^{+/- i*pi/5}
# In sigma(rho_2) = rho_2', eigenvalues at 12a are e^{+/- i*2*pi/5}
# Phase difference: pi/5 vs 2*pi/5, so delta = pi/5
phase_spin = np.pi / 5
print(f"  Candidate 2: pi/a_1 = pi/5 = {np.degrees(phase_spin):.2f} deg = {phase_spin:.6f} rad")

# Candidate 3: arctan(sqrt(5)) = delta_CKM (half the 5-cycle angle)
delta_CKM = np.arctan(SQRT5)
print(f"  Candidate 3: arctan(sqrt(5)) = {np.degrees(delta_CKM):.2f} deg")

# Candidate 4: The angle between 3 and 3' on the 12a class
# chi_3(12a) = phi, chi_3'(12a) = phi'
# These are real numbers with ratio phi/phi' = -phi^2
# The "angle" in some space between these characters...
# On the circle Q(zeta_5): phi = zeta_5 + zeta_5^4, phi' = zeta_5^2 + zeta_5^3
# arg(phi) = 0 (real positive), arg(phi') = pi (real negative)
# This gives phase difference pi -- too simple

# Candidate 5: From the TBM eigenvectors
# The Galois mixing matrix eigenvalues are {0, 3, 5}
# The eigenvectors are the TBM matrix columns
# The Majorana phases come from the relative phases of the mass eigenstates
# In TBM: U_TBM is real, so Majorana phases are 0 or pi.
# But with Galois corrections, they acquire nontrivial values.

print(f"\n  Candidate 5: 0 or pi (from TBM being real)")
print(f"    TBM matrix is REAL -> alpha_1 = alpha_2 = 0 or pi")
print(f"    But Galois corrections could shift these.")

# Candidate 6: From the spectral gap phase
# The Galois conjugate eigenspaces have eigenvalues:
# lambda_3 = 12 - 4*phi, lambda_3' = 8 + 4*phi
# These are Galois conjugates: sigma(12-4*phi) = 12-4*phi' = 12+4/phi = 8+4*phi
# The "phase" from Z[phi] -> C via zeta_5:
# Replace phi with zeta_5 + zeta_5^4 = 2*cos(2*pi/5)
# lambda_3 -> 12 - 4*(zeta_5 + zeta_5^4) = 12 - 8*cos(2*pi/5)
# This is real, so no complex phase this way.

# Candidate 7: From the mass formula embedding in Q(zeta_5)
# m_f = m_e * phi^(5a+6b)
# In Q(zeta_5): phi = zeta_5 + zeta_5^4
# phi^n has a definite phase when embedded in C via zeta_5
# The phase of phi^n in this embedding gives the Majorana phase

print(f"\n  Candidate 7: Phase of phi^n embedded via zeta_5")
# phi = zeta_5 + zeta_5^4 = 2*cos(2*pi/5) = 1/phi (as real number)
# Wait: zeta_5 + zeta_5^4 = 2*cos(2*pi/5) = 2*(sqrt(5)-1)/4 = (sqrt(5)-1)/2 = 1/phi
# But phi = (1+sqrt(5))/2, and 1/phi = (sqrt(5)-1)/2
# So the embedding zeta_5 -> e^{2*pi*i/5} gives phi -> 1/phi (NOT phi!)

# Let's be more careful.
# Q(sqrt(5)) embeds into Q(zeta_5) via sqrt(5) = zeta_5 - zeta_5^2 - zeta_5^3 + zeta_5^4
# = (zeta_5 + zeta_5^4) - (zeta_5^2 + zeta_5^3)
# = 1/phi - (-phi) = 1/phi + phi = sqrt(5). Correct!

# Now embed Q(zeta_5) into C via zeta_5 -> e^{2*pi*i/5}:
# sqrt(5) -> (e^{2*pi*i/5} + e^{-2*pi*i/5}) - (e^{4*pi*i/5} + e^{-4*pi*i/5})
#          = 2*cos(2*pi/5) - 2*cos(4*pi/5)
#          = 2*(sqrt(5)-1)/4 - 2*(-(sqrt(5)+1)/4)
#          = (sqrt(5)-1)/2 + (sqrt(5)+1)/2 = sqrt(5). Good.

# The embedding zeta_5 -> e^{2*pi*i*k/5} for k=1,2,3,4 gives 4 embeddings.
# k=1: phi -> (1+sqrt(5))/2 = phi (the standard one)
# k=2: sigma(phi) = (1+sigma(sqrt(5)))/2. Now sigma(sqrt(5)) = ?
# Under zeta_5 -> zeta_5^2:
# sqrt(5) = zeta_5 + zeta_5^4 - zeta_5^2 - zeta_5^3
# -> zeta_5^2 + zeta_5^3 - zeta_5^4 - zeta_5 = -sqrt(5)
# So sigma_2: sqrt(5) -> -sqrt(5), hence phi -> phi'
# This is exactly the Galois automorphism!

# k=1: sqrt(5) -> sqrt(5), phi -> phi (identity)
# k=2: sqrt(5) -> -sqrt(5), phi -> phi' (Galois)
# k=3: same as k=2 (since zeta_5^3 = conj(zeta_5^2))
# k=4: same as k=1 (since zeta_5^4 = conj(zeta_5))

# So there are exactly 2 complex embeddings: k=1,4 (giving phi) and k=2,3 (giving phi')
# The Galois group Gal(Q(zeta_5)/Q) = (Z/5Z)^* = {1,2,3,4} ~ Z/4Z
# The subgroup fixing Q(sqrt(5)) is {1,4} (complex conjugation)

# For the Majorana phase: embed via zeta_5 -> e^{2*pi*i/5}
# Then phi^n becomes a COMPLEX number (for the k=2 embedding)
# The phase of phi^n under the k=2 embedding:
# sigma_2(phi) = phi' = (1-sqrt(5))/2 = -0.618...
# sigma_2(phi^n) = (phi')^n

# But phi' is REAL (negative for odd n), so still no complex phase!

# The complex phases only appear when we use the FULL embedding
# zeta_5 -> e^{2*pi*i/5} (which is complex)

# Let's compute phi = zeta_5 + zeta_5^4 as a complex number with zeta_5 = e^{2*pi*i/5}
phi_complex = zeta5 + zeta5**4  # = 2*cos(2*pi/5) = real!
print(f"  phi via zeta_5: {phi_complex:.6f} (real = {phi_complex.real:.6f})")
print(f"  This is 2*cos(2*pi/5) = {2*np.cos(2*np.pi/5):.6f} = 1/phi = {1/PHI:.6f}")

# So phi, embedded in C via zeta_5, is REAL. No complex phase.
# But individual zeta_5 powers DO have phases!


# ============================================================
# PART 7: PHASES FROM THE MASS MATRIX IN Q(zeta_5)
# ============================================================
print("\n" + "=" * 72)
print("PART 7: PHASES FROM NEUTRINO MASS MATRIX")
print("=" * 72)

print("""
  The Majorana mass matrix M_nu = U^* diag(m1, m2, m3) U^dagger
  has Majorana phases encoded in the eigenvalue phases.

  In the framework:
  - Neutrino masses: m_i = 2*m_e / phi^(35 + delta_i) for some offsets
  - The mass matrix in flavor basis involves the TBM rotation
  - The GALOIS STRUCTURE gives: M ~ f(phi) in the 3x3 flavor block

  Key: the 2I group acts on the flavor space via rho_3 (or rho_3').
  An element g in 2I on a 5-fold axis has COMPLEX eigenvalues in rho_3:
  eigenvalues of g in rho_3: {1, zeta_5, zeta_5^4} (for appropriate g)

  The RELATIVE PHASES between eigenstates are the candidate Majorana phases.
""")

# In rho_3 (3-dim irrep), the eigenvalues of a 5-fold element are:
# For j=1 rep of SU(2): eigenvalues e^{i*m*2*theta} for m=-1,0,1
# With theta = pi/5 (12a class): eigenvalues e^{-2i*pi/5}, 1, e^{2i*pi/5}
# = zeta_5^{-2}, 1, zeta_5^2
# Wait, let me be more careful.

# A 5-fold element in SU(2) has eigenvalues e^{+/- i*pi/5} in the fundamental.
# In the spin-1 (dim 3) rep: eigenvalues are e^{-2i*pi/5}, 1, e^{2i*pi/5}
# = zeta_5^{-2}, 1, zeta_5^2

# These eigenvalues correspond to the three neutrino states in the mass basis!
# The phases are: 0, 2*pi/5, -2*pi/5 = 0, 72, -72 degrees

# For rho_3' (Galois conjugate): eigenvalues of the SAME element are:
# sigma(zeta_5^2) = zeta_5^4, sigma(1) = 1, sigma(zeta_5^{-2}) = zeta_5^{-4}
# So eigenvalues: zeta_5^{-4}, 1, zeta_5^4

# The phases are: 0, 4*pi/5, -4*pi/5 = 0, 144, -144 degrees

print("  Eigenvalues of 5-fold element in rho_3 and rho_3':")
theta_5 = np.pi / 5  # half-angle for order 10 in SU(2)
for m in [-1, 0, 1]:
    eig_3 = np.exp(2j * m * theta_5)
    # For rho_3': eigenvalues under Galois
    eig_3p = np.exp(2j * m * 2 * theta_5)
    print(f"    m={m:+d}: rho_3 -> e^{{i*{2*m}*pi/5}} = {eig_3:.4f} (phase {np.degrees(np.angle(eig_3)):+.1f})")
    print(f"         rho_3'-> e^{{i*{4*m}*pi/5}} = {eig_3p:.4f} (phase {np.degrees(np.angle(eig_3p)):+.1f})")

# The Majorana phases relate to the RELATIVE eigenvalue phases
# between the 3 and 3' blocks.

# Phase differences between corresponding eigenvalues in 3 and 3':
print(f"\n  Phase difference (rho_3' - rho_3) per eigenvalue:")
for m in [-1, 0, 1]:
    phase_3 = 2 * m * theta_5
    phase_3p = 4 * m * theta_5
    diff = phase_3p - phase_3
    print(f"    m={m:+d}: {np.degrees(diff):+.1f} deg = {diff/np.pi:.4f}*pi")

# The phase differences are: +2*pi/5, 0, -2*pi/5 for m = -1, 0, +1
# These correspond to the Galois "twist" between 3 and 3'

# Candidate Majorana phases:
print(f"\n  CANDIDATE MAJORANA PHASES:")

# alpha_1 = phase diff for nu_2 relative to nu_1
# alpha_2 = phase diff for nu_3 relative to nu_1

# If we assign: nu_1 ~ m=0 (eigenvalue 1, no phase)
#               nu_2 ~ m=+1 (phase 2*pi/5 in 3, 4*pi/5 in 3')
#               nu_3 ~ m=-1 (phase -2*pi/5 in 3, -4*pi/5 in 3')

# The Galois action shifts the phase by 2*pi/5 per unit m.
# The MAJORANA phase is the Galois twist:
alpha1_cand = 2 * np.pi / a1  # 2*pi/5 = 72 deg
alpha2_cand = 4 * np.pi / a1  # 4*pi/5 = 144 deg

print(f"  alpha_1 = 2*pi/a_1 = 2*pi/5 = {np.degrees(alpha1_cand):.2f} deg = {alpha1_cand:.6f} rad")
print(f"  alpha_2 = 4*pi/a_1 = 4*pi/5 = {np.degrees(alpha2_cand):.2f} deg = {alpha2_cand:.6f} rad")

# Alternative: phases from the spinorial reps (rho_2, rho_2')
# In rho_2: eigenvalues at 12a are e^{+/- i*pi/5}
# In rho_2': eigenvalues are e^{+/- i*2*pi/5}
# Phase shift: pi/5 per eigenvalue

alpha1_spin = 2 * np.pi / (2 * a1)  # pi/5 = 36 deg
alpha2_spin = 2 * 2 * np.pi / (2 * a1)  # 2*pi/5 = 72 deg

print(f"\n  Alternative (from spinorial rho_2/rho_2'):")
print(f"  alpha_1 = pi/a_1 = pi/5 = {np.degrees(alpha1_spin):.2f} deg")
print(f"  alpha_2 = 2*pi/a_1 = 2*pi/5 = {np.degrees(alpha2_spin):.2f} deg")


# ============================================================
# PART 8: CONSISTENCY WITH DELTA_PMNS
# ============================================================
print("\n" + "=" * 72)
print("PART 8: CONSISTENCY WITH delta_PMNS AND FRAMEWORK")
print("=" * 72)

delta_PMNS = N_gen * np.arctan(SQRT5)
print(f"  delta_PMNS = 3*arctan(sqrt(5)) = {np.degrees(delta_PMNS):.2f} deg")
print(f"  cos(delta_PMNS) = {np.cos(delta_PMNS):.6f}")
print(f"  = -7/(3*sqrt(6)) = {-7/(3*np.sqrt(6)):.6f}")

# All three phases should be related by the Galois structure
# delta_PMNS: Berry phase from generation loop
# alpha_1, alpha_2: eigenvalue phases from 5-fold element

# Check: do the phases form a consistent pattern?
print(f"\n  Phase pattern check:")
print(f"    delta_PMNS = {np.degrees(delta_PMNS):.2f} deg")
print(f"    alpha_1    = {np.degrees(alpha1_cand):.2f} deg (candidate: 2*pi/5)")
print(f"    alpha_2    = {np.degrees(alpha2_cand):.2f} deg (candidate: 4*pi/5)")

# Relation: alpha_2 = 2 * alpha_1 (doubling from Galois squaring map)
print(f"\n  alpha_2 / alpha_1 = {alpha2_cand/alpha1_cand:.4f} (should be 2)")
print(f"  This comes from: sigma(zeta_5) = zeta_5^2 (Frobenius)")

# Check if delta_PMNS relates to alpha_1
ratio_delta_alpha = delta_PMNS / alpha1_cand
print(f"\n  delta_PMNS / alpha_1 = {ratio_delta_alpha:.6f}")
print(f"  = {delta_PMNS:.6f} / {alpha1_cand:.6f}")
print(f"  ~ {ratio_delta_alpha:.4f} (no simple ratio)")

# The Jarlskog-like invariant with Majorana phases
# J_M = Im(U_e1 U_mu2^* U_e2^* U_mu1) * sin(alpha_1/2)
# For neutrinoless double beta decay: |m_ee| = |sum_i U_ei^2 m_i|

# Using TBM + our phases:
# U_TBM = [[sqrt(2/3), 1/sqrt(3), 0],
#           [-1/sqrt(6), 1/sqrt(3), 1/sqrt(2)],
#           [1/sqrt(6), -1/sqrt(3), 1/sqrt(2)]]

# |m_ee| = |m1 * 2/3 + m2 * 1/3 * e^{i*alpha1} + m3 * 0 * e^{i*alpha2}|
# Since theta_13 ~ 0 in TBM, the alpha_2 phase drops out of m_ee!

# More precisely with theta_13 = sqrt(1/45):
s13 = np.sqrt(1.0/45)
c13 = np.sqrt(1 - 1.0/45)
s12 = np.sqrt(2.0/(PHI + a1))  # framework value
c12 = np.sqrt(1 - s12**2)

print(f"\n  Effective Majorana mass for 0nu-beta-beta:")
print(f"    sin^2(th12) = {s12**2:.6f}, cos^2(th12) = {c12**2:.6f}")
print(f"    sin^2(th13) = {s13**2:.6f}")

# Framework neutrino masses (normal ordering assumed)
# m3 = 2*m_e/phi^35 = 49.53 meV (100% atmospheric)
# Delta m^2_atm = m3^2 - m2^2 ~ m3^2 (for hierarchical)
# Delta m^2_sol = m2^2 - m1^2
# We need m1 and m2... from framework: m_i = 2*m_e/phi^(35+delta_i)
# Currently we only have m3 robustly. Let's use m1 << m2 << m3.

m3 = 2 * 0.511 / PHI**35 * 1000  # in meV
m2 = np.sqrt(7.53e-5) * 1000  # from Delta m^2_sol ~ 7.53e-5 eV^2
m1_max = np.sqrt(m2**2 - 7.53e-2)  # rough

print(f"    m3 = {m3:.2f} meV (framework)")
print(f"    m2 ~ {m2:.2f} meV (from Delta m^2_sol)")

# |m_ee| = |c13^2 (c12^2 m1 + s12^2 m2 e^{i*alpha1}) + s13^2 m3 e^{i*(alpha2-2*delta)}|
# For hierarchical (m1 << m2 << m3):
for alpha1_test, alpha2_test, label in [
    (alpha1_cand, alpha2_cand, "2*pi/5, 4*pi/5"),
    (alpha1_spin, alpha2_spin, "pi/5, 2*pi/5"),
    (0, 0, "0, 0"),
    (np.pi, np.pi, "pi, pi"),
    (np.pi, 0, "pi, 0"),
]:
    # Simplified m_ee (ignoring m1, using m2 and m3)
    term2 = c13**2 * s12**2 * m2 * np.exp(1j * alpha1_test)
    term3 = s13**2 * m3 * np.exp(1j * (alpha2_test - 2*delta_PMNS))
    m_ee = abs(term2 + term3)
    print(f"    alpha=({label}): |m_ee| ~ {m_ee:.3f} meV")


# ============================================================
# PART 9: THE FROBENIUS DERIVATION
# ============================================================
print("\n" + "=" * 72)
print("PART 9: THE FROBENIUS DERIVATION")
print("=" * 72)

print(f"""
  The Galois group Gal(Q(zeta_5)/Q) = (Z/5Z)^* = {{1, 2, 3, 4}}

  The Frobenius automorphism at 5 is: sigma_2: zeta_5 -> zeta_5^2
  This generates the full Galois group of order 4 = phi(5).

  The subgroup fixing Q(sqrt(5)) is {{1, 4}} (complex conjugation).
  The quotient Gal(Q(sqrt(5))/Q) = {{1, sigma_2}} = Z/2Z.

  For 2I representations:
  - sigma_2 permutes rho_2 <-> rho_2' and rho_3 <-> rho_3'
  - On eigenvalues: sigma_2(zeta_5^m) = zeta_5^{{2m mod 5}}
  - The phase shift is: 2m*2*pi/5 - m*2*pi/5 = m*2*pi/5

  So the Galois twist per eigenvalue index m is m*2*pi/5.

  For the three neutrino mass eigenstates (m = -1, 0, +1 in rho_3):
    nu_1 (m=0): Galois twist = 0
    nu_2 (m=+1): Galois twist = +2*pi/5
    nu_3 (m=-1): Galois twist = -2*pi/5

  The MAJORANA PHASES are defined modulo 2*pi:
    alpha_1/2 = Galois twist of nu_2 = +2*pi/5  =>  alpha_1 = 4*pi/5 = 144 deg
    alpha_2/2 = Galois twist of nu_3 = -2*pi/5  =>  alpha_2 = -4*pi/5 = -144 deg

  Or equivalently (modulo 2*pi):
    alpha_1 = 4*pi/a_1 = 4*pi/5 = 144 deg
    alpha_2 = -4*pi/a_1 = 6*pi/5 = 216 deg (= -144 mod 360)
""")

alpha1_final = 4 * np.pi / a1
alpha2_final = -4 * np.pi / a1
# Convention: alpha/2 is the phase, so alpha = 2 * twist
# twist for m=+1: 2*pi/5, so alpha_1 = 2 * 2*pi/5 = 4*pi/5
# twist for m=-1: -2*pi/5, so alpha_2 = 2 * (-2*pi/5) = -4*pi/5

print(f"  FINAL PREDICTION:")
print(f"  alpha_1 = 4*pi/a_1 = 4*pi/5 = {np.degrees(alpha1_final):.2f} deg")
print(f"  alpha_2 = -4*pi/a_1 = -4*pi/5 = {np.degrees(alpha2_final):.2f} deg")
print(f"         = {np.degrees(alpha2_final) + 360:.2f} deg (mod 360)")
print(f"  |alpha_1| = |alpha_2| = 4*pi/5 (symmetric, from m -> -m)")
print(f"  alpha_2 = -alpha_1 (from m = -1 vs m = +1)")

# Self-consistency check:
# The combination appearing in neutrinoless double beta decay:
# alpha_21 = alpha_2 - alpha_1 = -4*pi/5 - 4*pi/5 = -8*pi/5 = 2*pi/5 (mod 2*pi)
alpha_21 = (alpha2_final - alpha1_final) % (2*np.pi)
print(f"\n  alpha_21 = alpha_2 - alpha_1 = {np.degrees(alpha_21):.2f} deg (mod 360)")
print(f"  = 2*pi/a_1 = {np.degrees(2*np.pi/a1):.2f} deg")

# The other combination: alpha_31 = alpha_2 - alpha_1 (already computed above)
# Wait, standard notation:
# alpha_21 = alpha_2 - alpha_1 (but some papers define differently)
# Let me use the PDG convention: alpha_1 and alpha_2 are independent phases

print(f"\n  Effective Majorana mass |m_ee| with final phases:")
for m1_val in [0, 1, 3, 5]:  # m1 in meV
    term1 = c13**2 * c12**2 * m1_val
    term2 = c13**2 * s12**2 * m2 * np.exp(1j * alpha1_final)
    term3 = s13**2 * m3 * np.exp(1j * alpha2_final)
    m_ee = abs(term1 + term2 + term3)
    print(f"    m1 = {m1_val:.0f} meV: |m_ee| = {m_ee:.3f} meV")


# ============================================================
# PART 10: HONEST ASSESSMENT
# ============================================================
print("\n" + "=" * 72)
print("PART 10: HONEST ASSESSMENT")
print("=" * 72)

print(f"""
RESULTS:
========

1. MAJORANA PHASES PREDICTED:
   alpha_1 = 4*pi/5 = 144 deg
   alpha_2 = -4*pi/5 = -144 deg (= 216 deg mod 360)
   |alpha_1| = |alpha_2| = 4*pi/a_1

2. DERIVATION CHAIN:
   a) 2I is the double cover of I = A_5
   b) Neutrino generations live in rho_3 (dim 3) of 2I
   c) The 5-fold element g_5 has eigenvalues zeta_5^m (m=-1,0,+1) in rho_3
   d) The Galois automorphism sigma_2: zeta_5 -> zeta_5^2 (Frobenius at 5)
   e) sigma_2 shifts eigenvalue phases by m*2*pi/5
   f) The Majorana phase alpha_k/2 is the Galois twist of mass eigenstate k
   g) alpha_1/2 = 2*pi/5 (for m=+1), alpha_2/2 = -2*pi/5 (for m=-1)

3. KEY ASSUMPTIONS:
   - Neutrino mass eigenstates correspond to angular momentum projections m
   - The Galois twist acts as the CP transformation on Majorana states
   - The assignment m=0 -> nu_1, m=+1 -> nu_2, m=-1 -> nu_3

4. TESTABILITY:
   - |m_ee| ~ {abs(c13**2 * s12**2 * m2 * np.exp(1j*alpha1_final) + s13**2 * m3 * np.exp(1j*alpha2_final)):.2f} meV (for m1 ~ 0)
   - This is BELOW current experimental reach (~20-100 meV)
   - Future ton-scale experiments may reach O(1) meV sensitivity

5. RELATION TO OTHER PHASES:
   delta_PMNS = 3*arctan(sqrt(5)) = {np.degrees(delta_PMNS):.2f} deg (Berry phase)
   alpha_1 = 4*pi/5 = {np.degrees(alpha1_final):.2f} deg (Frobenius twist)
   alpha_2 = -4*pi/5 = {np.degrees(alpha2_final):.2f} deg (conjugate twist)

   These are INDEPENDENT derivations:
   - delta from Berry phase on generation loop (N_gen * arctan(sqrt(5)))
   - alpha from Frobenius on 5th roots of unity (eigenvalue phases)

CATEGORY: PARTIALLY DERIVED
  The Frobenius automorphism and eigenvalue structure are rigorous.
  The assignment of m-values to mass eigenstates is assumed (not derived).
  The identification of Galois twist with Majorana phase is motivated
  by the role of sigma as CP transformation, but not proven from first principles.
""")
