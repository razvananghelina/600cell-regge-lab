#!/usr/bin/env python3
"""Exact audit of the Perron-normalized Fibonacci fusion dynamics."""

import sympy as sy


passed = 0
tests = 0


def check(name, condition, detail=""):
    global passed, tests
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"       {detail}")


sqrt5 = sy.sqrt(5)
phi = (1 + sqrt5) / 2
sigma_phi = (1 - sqrt5) / 2

# Basis (1, X), with X tensor X = 1 + X.  N acts by fusion with X.
N = sy.Matrix([[0, 1], [1, 1]])
d = sy.Matrix([1, phi])

check("Fibonacci fusion Perron vector is exact",
      sy.simplify(N * d - phi * d) == sy.zeros(2, 1))

# Doob/Perron normalization: P[a,b] = N[a,b] d[b]/(phi d[a]).
P = sy.simplify(sy.diag(*[1 / x for x in d]) * N * sy.diag(*d) / phi)
check("P is entrywise positive where fusion is allowed",
      all(x >= 0 for x in P), str(P))
check("P is exactly stochastic", P * sy.ones(2, 1) == sy.ones(2, 1))

D2 = sy.simplify(sum(x**2 for x in d))
pi = sy.Matrix([[sy.simplify(x**2 / D2) for x in d]])
check("quantum-dimension squares give the stationary law",
      sy.simplify(pi * P - pi) == sy.zeros(1, 2), str(pi))
check("the dynamics obeys detailed balance",
      sy.simplify(pi[0, 0] * P[0, 1] - pi[0, 1] * P[1, 0]) == 0)

eigenvalues = set(P.eigenvals())
galois_ratio = sy.simplify(sigma_phi / phi)
check("the nonstationary mode is the normalized Galois conjugate",
      eigenvalues == {sy.Integer(1), galois_ratio},
      f"spectrum={eigenvalues}, sigma(phi)/phi={galois_ratio}")
check("the Galois mode is alternating and contractive",
      sy.simplify(galois_ratio + phi**-2) == 0
      and abs(float(galois_ratio)) < 1)

# Applying the field automorphism entrywise preserves row sums but destroys
# positivity.  It is therefore a signed conjugate dynamics, not a second
# stochastic world.
P_sigma = P.applyfunc(lambda x: sy.simplify(x.xreplace({sqrt5: -sqrt5})))
check("Galois-conjugate transition matrix retains total weight one",
      P_sigma * sy.ones(2, 1) == sy.ones(2, 1), str(P_sigma))
check("Galois-conjugate transition matrix is not probabilistic",
      any(x < 0 for x in P_sigma) and any(x > 1 for x in P_sigma))

# Attack the uniqueness claim.  Every rank-two based ring
# X^2 = 1 + m X has the same construction.  Fibonacci is m=1, not selected
# by Perron normalization alone.
m = sy.symbols("m", positive=True, integer=True)
delta = sy.sqrt(m**2 + 4)
d_m = (m + delta) / 2
sigma_d_m = (m - delta) / 2
P_m = sy.Matrix([[0, 1], [d_m**-2, m / d_m]])
check("the rank-two family is stochastic symbolically",
      sy.simplify(P_m * sy.ones(2, 1) - sy.ones(2, 1)) == sy.zeros(2, 1))
check("the conjugate-ratio mode persists for every rank-two m",
      sy.simplify(P_m.det() - sigma_d_m / d_m) == 0)

# For integer m >= 1, d_m strictly increases, so |sigma(d_m)/d_m|=d_m^-2
# strictly decreases.  A finite exact audit is only a sanity check; the
# derivative d'(m)=(1+m/sqrt(m^2+4))/2 > 0 is the proof.
derivative = sy.diff((m + sy.sqrt(m**2 + 4)) / 2, m)
check("Fibonacci maximizes conjugate-mode persistence in the m>=1 family",
      derivative.is_positive
      and all(sy.simplify(
          ((k + sy.sqrt(k*k + 4))/2)**-2
          - ((k + 1 + sy.sqrt((k + 1)**2 + 4))/2)**-2) > 0
              for k in range(1, 20)),
      "conditional on the rank-two family X^2=1+mX")

# The independently known C10 cycle has the same gap phi^-2.  Attack the
# naive identification: its gap eigenspace is a real plane, not a canonical
# real line.  Equivariance forces a phase-bearing two-dimensional completion.
S = sy.zeros(10)
for i in range(10):
    S[i, (i + 1) % 10] = 1
L = 2 * sy.eye(10) - S - S.T
gap = sy.simplify(phi**-2)
gap_space = (L - gap * sy.eye(10)).nullspace()
check("the C10 fiber gap equals the Fibonacci decay magnitude",
      sy.simplify(gap + galois_ratio) == 0)
check("the C10 gap eigenspace is exactly a real plane",
      len(gap_space) == 2)
check("the gap plane contains no C10-invariant line with trivial action",
      sy.Matrix.vstack(L - gap * sy.eye(10), S - sy.eye(10)).nullspace() == [])
phase_polynomial_holds = all(
    sy.simplify((S**2 - phi * S + sy.eye(10)) * vector)
    == sy.zeros(10, 1)
    for vector in gap_space
)
check("fiber rotation on the gap plane has phase polynomial t^2-phi*t+1",
      phase_polynomial_holds,
      "eigenphases exp(+- i*pi/5); minimal real completion has dimension 2")

# The phase lift is cyclotomic and unique up to orientation.  Taking the
# field norm Q(sqrt(5)) -> Q of its quadratic polynomial gives Phi_10.
z = sy.symbols("z")
phase_poly = z**2 - phi*z + 1
phase_poly_sigma = z**2 - sigma_phi*z + 1
cyclotomic_10 = sy.cyclotomic_poly(10, z)
check("the phase lift has exact cyclotomic norm Phi_10",
      sy.expand(phase_poly * phase_poly_sigma - cyclotomic_10) == 0,
      f"Norm(z^2-phi*z+1)={cyclotomic_10}")
phase_roots = sy.solve(phase_poly, z)
check("the unitary lift is unique up to orientation",
      len(phase_roots) == 2
      and sy.simplify(phase_roots[0] * phase_roots[1] - 1) == 0
      and all(sy.simplify(root * sy.conjugate(root) - 1) == 0
              for root in phase_roots),
      "roots are exp(+-i*pi/5)")

# Hostile categorical boundary.  A strong tensor functor to ordinary finite
# dimensional representations would send X to an integer-dimensional V with
# V tensor V = 1 + V, hence n^2=1+n, which has no nonnegative integer root.
integer_dimension_solutions = [n for n in range(100) if n*n == 1+n]
check("no ordinary finite-dimensional fiber functor can realize Fibonacci",
      integer_dimension_solutions == [])

# Any group grading assigns degree g to X.  Since both 1 and X occur in X^2,
# it requires g^2=e and g^2=g, hence g=e.  Verify this exhaustively for C10.
grading_solutions = [
    g for g in range(10) if (2*g) % 10 == 0 and (2*g) % 10 == g
]
check("every C10 grading of the Fibonacci fusion rule is trivial",
      grading_solutions == [0])

# Exact induction of the phase character C10 -> 2I.  Character-table class
# order: 1A,2A,4A,6A,3A,10A,5A,5B,10B.  The irreducibles are generated from
# the two Galois-conjugate defining spinors by the SU(2) recurrence.
phip = sigma_phi
x = (2, -2, 0, 1, -1, phi, -phi, phi-1, phip)
xp = (2, -2, 0, 1, -1, phip, phi-1, -phi, phi)


def symmetric_power_character(power, values):
    if power == 0:
        return tuple(sy.Integer(1) for _ in values)
    if power == 1:
        return values
    previous2 = tuple(sy.Integer(1) for _ in values)
    previous1 = values
    for _ in range(2, power + 1):
        current = tuple(sy.simplify(v*a-b)
                        for v, a, b in zip(values, previous1, previous2))
        previous2, previous1 = previous1, current
    return previous1


irreps_2i = (
    symmetric_power_character(0, x),
    symmetric_power_character(1, x),
    symmetric_power_character(1, xp),
    symmetric_power_character(2, x),
    symmetric_power_character(2, xp),
    tuple(sy.expand(a*b) for a, b in zip(x, xp)),
    symmetric_power_character(3, x),
    symmetric_power_character(4, x),
    symmetric_power_character(5, x),
)
irrep_dims = (1, 2, 2, 3, 3, 4, 4, 5, 6)

# Powers of an order-10 element in class 10A, certified by the defining
# traces 2*cos(k*pi/5).
c10_power_classes = (0, 5, 7, 8, 6, 1, 6, 8, 7, 5)
check("C10 power classes reproduce all defining-spinor traces",
      all(sy.simplify(x[class_index] - 2*sy.cos(sy.pi*k/5)) == 0
          for k, class_index in enumerate(c10_power_classes)))

zeta10 = sy.exp(sy.I * sy.pi / 5)


def induced_multiplicities(harmonic):
    answer = []
    for character in irreps_2i:
        inner = sum(character[class_index] * zeta10**(-harmonic*k)
                    for k, class_index in enumerate(c10_power_classes)) / 10
        answer.append(sy.simplify(sy.expand_complex(inner)))
    return tuple(answer)


induced_1 = induced_multiplicities(1)
induced_3 = induced_multiplicities(3)
expected_1 = (0, 1, 0, 0, 0, 0, 1, 0, 1)
expected_3 = (0, 0, 1, 0, 0, 0, 1, 0, 1)
check("Ind_C10^2I chi_1 is exactly 2+4+6",
      induced_1 == expected_1
      and sum(a*b for a, b in zip(induced_1, irrep_dims)) == 12,
      f"multiplicities={induced_1}")
check("Galois phase chi_3 induces exactly 2-prime+4+6",
      induced_3 == expected_3
      and sum(a*b for a, b in zip(induced_3, irrep_dims)) == 12,
      f"multiplicities={induced_3}")

common = tuple(min(a, b) for a, b in zip(induced_1, induced_3))
amalgamated_union = tuple(max(a, b) for a, b in zip(induced_1, induced_3))
odd_mckay = (0, 1, 1, 0, 0, 0, 1, 0, 1)
check("the two phase characters have common irreducible content 4+6",
      common == (0, 0, 0, 0, 0, 0, 1, 0, 1)
      and sum(a*b for a, b in zip(common, irrep_dims)) == 10)
check("their multiplicity-wise union is the 14D odd McKay half",
      amalgamated_union == odd_mckay
      and sum(a*b for a, b in zip(amalgamated_union, irrep_dims)) == 14,
      "2 + 2-prime + 4 + 6")
check("the induced phase modules are genuinely spinorial",
      sum(mult*character[1]
          for mult, character in zip(induced_1, irreps_2i)) == -12
      and sum(mult*character[1]
              for mult, character in zip(induced_3, irreps_2i)) == -12,
      "central -1 acts with character -12")

# Gluing-phase audit.  Since 4 and 6 occur once, Schur's lemma gives
# End_2I(4+6) = C + C and unitary automorphisms U(1)^2.  This is an ambiguity
# before imposing McKay adjacency, not yet a gauge group.
check("the common 4+6 content has a two-complex-dimensional commutant",
      sum(mult**2 for mult in common) == 2,
      "Aut_2I(4+6) unitary = U(1)^2")

class_sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)
mckay_edges = []
for i in range(9):
    for j in range(i + 1, 9):
        multiplicity = sy.simplify(sum(
            class_sizes[c] * x[c] * irreps_2i[i][c] * irreps_2i[j][c]
            for c in range(9)
        ) / 120)
        if multiplicity:
            mckay_edges.append((i, j, multiplicity))
expected_mckay_edges = [
    (0, 1, 1), (1, 3, 1), (2, 5, 1), (3, 6, 1),
    (4, 8, 1), (5, 8, 1), (6, 7, 1), (7, 8, 1),
]
check("exact character products reproduce the affine-E8 McKay tree",
      mckay_edges == expected_mckay_edges)

# A diagonal node phase commutes with adjacency iff its endpoint phases agree
# on every edge.  The incidence matrix has rank 8, hence only the global
# scalar phase survives; modulo that physically trivial scalar, dimension 0.
incidence = sy.zeros(len(mckay_edges), 9)
for row, (i, j, _) in enumerate(mckay_edges):
    incidence[row, i] = 1
    incidence[row, j] = -1
check("McKay compatibility kills every relative gluing phase",
      incidence.rank() == 8 and len(incidence.nullspace()) == 1,
      "only global U(1) commutes with D; effective relative torus is trivial")

# Canonical spectral-carrier screen on M1 (+ chirality) plus M3 (- chirality).
# Schur intertwiners exist only on the common 4+6 summands.  Consequently an
# equivariant off-diagonal D has rank at most 2*(4+6)=20 on H of dimension 24.
hom_multiplicity = sum(a*b for a, b in zip(induced_1, induced_3))
common_dimension = sum(min(a, b)*d0
                       for a, b, d0 in zip(induced_1, induced_3, irrep_dims))
check("equivariant Galois-sheet Dirac has two Schur channels",
      hom_multiplicity == 2 and common_dimension == 10)
check("every equivariant Dirac on the 12+12 carrier has endpoint kernel",
      2*(12-common_dimension) == 4,
      "max rank 20 on H=24; the 2 and 2-prime endpoints cannot couple")

# Canonical central node algebra on irreducible types (2,2',4,6).  Its
# graded multiplicities are (+2,-2,0,0), and J permutes 2<->2' while fixing
# 4 and 6.  Therefore Cap_ij = s_i delta_{i,sigma(j)} has rank only two.
graded_dimensions = sy.Matrix([2, -2, 0, 0])
galois_nodes = (1, 0, 2, 3)
cap_phase_carrier = sy.zeros(4)
for i in range(4):
    for j in range(4):
        if i == galois_nodes[j]:
            cap_phase_carrier[i, j] = graded_dimensions[i]
check("canonical C^4 node algebra fails Poincare duality",
      cap_phase_carrier.rank() == 2 and cap_phase_carrier.det() == 0,
      f"Cap={cap_phase_carrier.tolist()}")
check("endpoint projectors make the canonical carrier disconnected",
      12-common_dimension > 0,
      "the uncoupled 2 and 2-prime projectors commute with every equivariant D")

# Full Hopf-harmonic support theorem.  Orientation identifies q with -q at
# the character-support level, so q=0,...,5 suffice.  Even harmonics induce
# precisely the even McKay nodes; odd harmonics induce precisely the odd
# spinorial nodes.  This is support union, not direct-sum multiplicity.
all_induced = tuple(induced_multiplicities(q) for q in range(10))
expected_induced = (
    (1, 0, 0, 1, 1, 0, 0, 1, 0),
    (0, 1, 0, 0, 0, 0, 1, 0, 1),
    (0, 0, 0, 1, 0, 1, 0, 1, 0),
    (0, 0, 1, 0, 0, 0, 1, 0, 1),
    (0, 0, 0, 0, 1, 1, 0, 1, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 2),
    (0, 0, 0, 0, 1, 1, 0, 1, 0),
    (0, 0, 1, 0, 0, 0, 1, 0, 1),
    (0, 0, 0, 1, 0, 1, 0, 1, 0),
    (0, 1, 0, 0, 0, 0, 1, 0, 1),
)
check("all ten C10 harmonic inductions are exact",
      all_induced == expected_induced
      and all(sum(a*b for a, b in zip(row, irrep_dims)) == 12
              for row in all_induced))


def support_union(rows):
    rows = tuple(rows)
    return tuple(int(any(row[i] for row in rows)) for i in range(9))


even_harmonic_support = support_union(all_induced[q] for q in (0, 2, 4))
odd_harmonic_support = support_union(all_induced[q] for q in (1, 3, 5))
even_mckay = (1, 0, 0, 1, 1, 1, 0, 1, 0)
check("even C10 harmonics reproduce exactly the 16D even McKay support",
      even_harmonic_support == even_mckay
      and sum(a*b for a, b in zip(even_harmonic_support, irrep_dims)) == 16)
check("odd C10 harmonics reproduce exactly the 14D odd McKay support",
      odd_harmonic_support == odd_mckay
      and sum(a*b for a, b in zip(odd_harmonic_support, irrep_dims)) == 14)
check("Hopf harmonic parity equals central spinorial parity",
      all(all(character[1] == ((-1)**q)*character[0]
              for mult, character in zip(all_induced[q], irreps_2i) if mult)
          for q in range(10)))
full_support = tuple(max(a, b)
                     for a, b in zip(even_harmonic_support, odd_harmonic_support))
check("even plus odd harmonic supports recover all nine McKay nodes",
      full_support == (1,)*9
      and sum(a*b for a, b in zip(full_support, irrep_dims)) == 30,
      "16+14=30 with no fitted node assignment")

# Canonical categorical realization: summing all five harmonics of a fixed
# parity retains their true multiplicities.  The result is exactly the two
# central ideals e_+ C[2I] and e_- C[2I] in the regular representation.
even_regular_mult = tuple(sum(all_induced[q][i] for q in (0, 2, 4, 6, 8))
                          for i in range(9))
odd_regular_mult = tuple(sum(all_induced[q][i] for q in (1, 3, 5, 7, 9))
                         for i in range(9))
expected_even_regular = tuple(d0 if sign else 0
                              for d0, sign in zip(irrep_dims, even_mckay))
expected_odd_regular = tuple(d0 if sign else 0
                             for d0, sign in zip(irrep_dims, odd_mckay))
check("even harmonic sum is the canonical central ideal e_+ C[2I]",
      even_regular_mult == expected_even_regular
      and sum(m*d0 for m, d0 in zip(even_regular_mult, irrep_dims)) == 60,
      f"multiplicities={even_regular_mult}")
check("odd harmonic sum is the canonical central ideal e_- C[2I]",
      odd_regular_mult == expected_odd_regular
      and sum(m*d0 for m, d0 in zip(odd_regular_mult, irrep_dims)) == 60,
      f"multiplicities={odd_regular_mult}")
check("the two parity ideals reconstruct the 120D regular representation",
      tuple(a+b for a, b in zip(even_regular_mult, odd_regular_mult))
      == irrep_dims)

# Hostile dynamical boundary: the two central ideals have disjoint irreducible
# support, so their equivariant Hom space is zero.  Tensoring by the defining
# spinor connects their simple spectra, but that is a functor, not an
# equivariant odd endomorphism of the regular module.
equivariant_odd_hom_dimension = sum(a*b for a, b
                                    in zip(even_regular_mult, odd_regular_mult))
check("no nonzero 2I-equivariant odd Dirac connects the central ideals",
      equivariant_odd_hom_dimension == 0,
      "McKay adjacency is tensor-by-2 functorial data, not such an endomorphism")

# The correct equivariant object is an odd correspondence: tensoring by the
# defining spinor.  On multiplicity vectors it is the McKay adjacency.  It
# sends either regular parity ideal to two copies of the opposite ideal.
mckay_adj = sy.zeros(9)
for i, j, multiplicity in mckay_edges:
    mckay_adj[i, j] = multiplicity
    mckay_adj[j, i] = multiplicity
even_mult_vector = sy.Matrix(even_regular_mult)
odd_mult_vector = sy.Matrix(odd_regular_mult)
check("tensor-by-2 sends the even regular ideal to twice the odd ideal",
      mckay_adj * even_mult_vector == 2 * odd_mult_vector)
check("tensor-by-2 sends the odd regular ideal to twice the even ideal",
      mckay_adj * odd_mult_vector == 2 * even_mult_vector)

# Perron/Doob normalization of the McKay correspondence gives a canonical
# stochastic dynamics on simple types.
dimension_vector = sy.Matrix(irrep_dims)
check("McKay dimensions form the exact Perron eigenvector",
      mckay_adj * dimension_vector == 2 * dimension_vector)
P_mckay = sy.zeros(9)
for i in range(9):
    for j in range(9):
        P_mckay[i, j] = sy.Rational(1, 2) * mckay_adj[i, j] \
            * sy.Rational(irrep_dims[j], irrep_dims[i])
check("Perron-normalized McKay correspondence is stochastic",
      P_mckay * sy.ones(9, 1) == sy.ones(9, 1))
plancherel = sy.Matrix([[sy.Rational(d0*d0, 120) for d0 in irrep_dims]])
check("Plancherel weights are the exact stationary law",
      plancherel * P_mckay == plancherel
      and sum(plancherel) == 1)
node_parity = tuple(1 if x0 else -1 for x0 in even_mckay)
check("every McKay transition flips Hopf/central parity",
      all(node_parity[i] == -node_parity[j]
          for i, j, _ in mckay_edges))
check("stationary weight is balanced exactly 1/2 plus and 1/2 minus",
      sum(plancherel[0, i] for i in range(9) if node_parity[i] == 1)
      == sy.Rational(1, 2)
      and sum(plancherel[0, i] for i in range(9) if node_parity[i] == -1)
      == sy.Rational(1, 2))

# Spectral closure back to the Fibonacci seed.  P_McKay is similar to A/2;
# its spectrum is the defining 2I character divided by two.  The ratio of
# the Galois-conjugate golden modes is exactly the original two-state memory
# eigenvalue sigma(phi)/phi = -phi^-2.
t = sy.symbols("t")
expected_mckay_charpoly = sy.expand(
    t * (t**2 - 1) * (t**2 - sy.Rational(1, 4))
    * (t**4 - sy.Rational(3, 4)*t**2 + sy.Rational(1, 16))
)
check("McKay Markov characteristic polynomial factors exactly",
      sy.expand(P_mckay.charpoly(t).as_expr() - expected_mckay_charpoly) == 0,
      f"charpoly={sy.factor(expected_mckay_charpoly)}")
expected_spectrum = {
    sy.Integer(1), sy.Integer(-1), sy.Integer(0),
    sy.Rational(1, 2), sy.Rational(-1, 2),
    sy.simplify(phi/2), sy.simplify(-phi/2),
    sy.simplify(1/(2*phi)), sy.simplify(-1/(2*phi)),
}
check("McKay Markov spectrum is the normalized defining character",
      set(P_mckay.eigenvals()) == expected_spectrum)
check("Fibonacci memory is the projective Galois ratio of McKay modes",
      sy.simplify((sigma_phi/2)/(phi/2) - galois_ratio) == 0
      and sy.simplify((sigma_phi/2)/(phi/2) + phi**-2) == 0,
      "(sigma(phi)/2)/(phi/2)=sigma(phi)/phi=-phi^-2")

# Two steps return to the same chirality.  The nontrivial slow eigenvalue in
# either parity block is phi^2/4, giving exact two-step spectral gap.
P2 = sy.simplify(P_mckay**2)
even_indices = [i for i, sign in enumerate(node_parity) if sign == 1]
odd_indices = [i for i, sign in enumerate(node_parity) if sign == -1]
P2_even = P2.extract(even_indices, even_indices)
P2_odd = P2.extract(odd_indices, odd_indices)
slow_mode = sy.simplify(phi**2 / 4)
fast_mode = sy.simplify(1/(4*phi**2))
u = sy.symbols("u")
expected_even_p2_charpoly = sy.expand(
    u*(u-1)*(u-sy.Rational(1, 4))*(u-slow_mode)*(u-fast_mode))
expected_odd_p2_charpoly = sy.expand(
    (u-1)*(u-sy.Rational(1, 4))*(u-slow_mode)*(u-fast_mode))
check("two-step parity chains have the predicted exact spectra",
      sy.simplify(P2_even.charpoly(u).as_expr()
                  - expected_even_p2_charpoly) == 0
      and sy.simplify(P2_odd.charpoly(u).as_expr()
                      - expected_odd_p2_charpoly) == 0)
check("same-chirality spectral gap is exact",
      sy.simplify(1 - slow_mode - (3-phi)/4) == 0,
      "gap=1-phi^2/4=(3-phi)/4")

print("-" * 72)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT: one positive stationary sheet plus one alternating Galois mode;")
print("         canonical given Fibonacci, but not a derivation of Fibonacci.")
if passed != tests:
    raise SystemExit(1)
