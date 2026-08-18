# Fibonacci non-binary dynamics: one process, two algebraic modes

## Result

The foundational rule

`X tensor X = 1 direct-sum X`

has fusion matrix

`N_X = [[0,1],[1,1]]`.

Let `d=(1,phi)` be its positive Perron vector.  The canonical Perron/Doob
normalization

`P_ab = (N_X)_ab d_b / (phi d_a)`

is

`P = [[0,1],[phi^-2,phi^-1]]`.

It is an exact stochastic matrix.  Its stationary distribution is
proportional to the squared quantum dimensions,

`pi = (1,phi^2)/(1+phi^2)`,

and it satisfies detailed balance.  **DERIVED conditional on Axiom S01.**

## The non-binary lead

The two eigenvalues of the single stochastic process are

`1` and `sigma(phi)/phi = -phi^-2`,

where `sigma(phi)=1-phi=-phi^-1` is the nontrivial Galois conjugate.  Thus the
physical Perron embedding and its algebraic conjugate occur simultaneously:

- the Perron mode is the positive stationary component;
- the conjugate mode is an alternating memory mode whose magnitude decays as
  `phi^(-2n)`.

No extra Hilbert-space sheet or trivial chirality double is inserted.
**DERIVED algebraic statement.**  Reading it ontologically as “reality is
non-binary” is **STRUCTURAL interpretation**, not a theorem of physics.

Applying `sigma` entrywise to `P` preserves row sums but produces a negative
entry and an entry greater than one.  Therefore the conjugate matrix is a
signed dynamics, not a second probability space.  Positivity selects one
embedding as probabilistic while retaining the other internally as a decay
mode.  **DERIVED.**

## Hostile boundary

This construction does not derive Fibonacci.  For every rank-two based rule

`X tensor X = 1 direct-sum m X`, `m>=1`,

Perron normalization gives a stochastic matrix with nontrivial eigenvalue

`sigma(d_m)/d_m = -d_m^-2`,

where `d_m=(m+sqrt(m^2+4))/2`.  Hence “one stationary mode plus one Galois
memory mode” is a property of the entire rank-two family.  **DERIVED
negative.**

Within that family, `d_m` grows strictly with `m`, so `m=1` uniquely maximizes
the persistence `|sigma(d_m)/d_m|=d_m^-2` among productive `m>=1` rules.
Thus Fibonacci is the maximally persistent non-pointed rank-two
self-reference.  **DERIVED conditional on having already restricted to this
rank-two family.**  Promoting maximal persistence to a law of nature would be
a new **STRUCTURAL axiom**, not a derivation.

## Concrete research lead

The next falsifiable bridge is to compare this intrinsic decay mode with the
independently obtained `phi^-2` spectral modes of the Hopf/fiber and geometric
operators.  Equality of eigenvalues alone is only a **PATTERN**.  A real
bridge requires an explicit intertwiner taking the fusion-chain conjugate
eigenvector to a geometrically defined `phi^-2` eigenspace while respecting
the relevant symmetry action.

For the `C10` cycle this first attempt has an exact answer.  Its Laplacian gap
is `phi^-2`, but the corresponding real eigenspace has dimension two.  It has
no nonzero vector fixed by the `C10` rotation, so the one-dimensional real
Fibonacci decay mode cannot map to it equivariantly if the source carries the
trivial `C10` action.  **DERIVED no-go for the naive intertwiner.**

The minimal equivariant completion is nevertheless informative.  On that
real gap plane the fiber rotation obeys

`S^2 - phi S + 1 = 0`,

and therefore has eigenphases `exp(+-i pi/5)`.  Equivalently, the required
completion of the real decay coordinate is one complex phase-bearing line.
**DERIVED.**  Interpreting this as an emergence of quantum phase from the
geometric completion of non-binary self-reference is a **STRUCTURAL lead**.
It gives only discrete `C10` phase at present, not a derived continuum
`U(1)` gauge symmetry.  A canonical functor producing this completion from
the Fibonacci category remains **OPEN**.

## Functor no-go and the surviving phase lift

The obvious functorial interpretation is actually impossible.  A strong
tensor functor from Fibonacci to ordinary finite-dimensional vector spaces
or `Rep(C10)` would send `X` to a representation of integer dimension `n`
obeying

`n^2 = 1+n`.

There is no nonnegative integer solution.  Likewise, assigning a group
degree `g` to `X` requires simultaneously `g^2=e` (because `1` occurs in
`X^2`) and `g^2=g` (because `X` occurs there), hence `g=e`.  Every `C10`
grading is trivial.  **DERIVED NO-GO:** the desired bridge is neither an
ordinary fiber functor nor a nontrivial group grading of the bare Fibonacci
fusion category.

What survives is not a tensor functor but a unique unitary phase lift.  The
real Perron value satisfies

`phi = z + z^-1`

precisely when

`z^2 - phi z + 1 = 0`.

Its two roots are `exp(+-i pi/5)`, exchanged by orientation reversal.  The
field norm of this polynomial is

`(z^2-phi z+1)(z^2-sigma(phi) z+1)`
`= z^4-z^3+z^2-z+1 = Phi_10(z)`.

Therefore the `C10` phase is the exact cyclotomic lift of the real Fibonacci
dimension, unique up to the two orientations.  **DERIVED.**  This is weaker
than a categorical functor but stronger than a numerical match.

The proposed foundational reading is: a real observable `phi` is the
orientation-blind trace of a conjugate pair of phases `z,z^-1`; geometry
restores the orientation information suppressed by the real trace.
**STRUCTURAL interpretation.**  Whether this trace/lift operation is the
correct physical meaning of “non-binary reality” remains **OPEN**.

## Exact induction to the binary icosahedral McKay sector

The phase lift supplies the fundamental character

`chi_1(g)=exp(i pi/5)`

of a Hopf subgroup `C10 < 2I`.  Frobenius reciprocity and the exact `2I`
character table give

`Ind_C10^2I(chi_1) = 2 direct-sum 4 direct-sum 6`.

The dimensions sum to 12, as required by the subgroup index.  Every `C10`
subgroup is conjugate, so the isomorphism class of the induced module does
not depend on which of the six Hopf carriers is used.  Reversing orientation
(`chi_1 -> chi_-1`) gives the same complex `2I` character.  **DERIVED.**

The arithmetic Galois operation sends the primitive harmonic `1` to `3` and
exchanges the two defining spinors.  Exact induction gives

`Ind(chi_3) = 2-prime direct-sum 4 direct-sum 6`.

Thus the two Galois phase characters have common irreducible content
`4 direct-sum 6`, of total dimension 10.  Their multiplicity-wise union in
the representation semiring is

`2 direct-sum 2-prime direct-sum 4 direct-sum 6`,

which is exactly the 14-dimensional odd half of the affine-`E8` McKay node
module already derived independently in the repository.  This equality is
**DERIVED at the level of representation isomorphism classes.**  It does not
canonically identify the two concrete copies of `4 direct-sum 6`; an actual
amalgamation requires intertwiners and remains **STRUCTURAL/OPEN**.  The central
element `-1` has character `-12` on each induced module, so these are
genuinely spinorial rather than `A5`-factorized modules.  **DERIVED.**

This is an exact bridge

`Fibonacci Perron value -> cyclotomic C10 phase -> induced 2I module`

and not merely equality of dimensions.  However, it constructs a
representation/module, not a finite spectral triple, a Standard-Model
matter representation, or a dynamical law.  Those upgrades remain **OPEN**.
Calling the shared `4+6` sector physical matter is currently **PATTERN**.

## Gluing-phase audit: no gauge `U(1)` yet

Because the common irreducibles `4` and `6` each occur with multiplicity one,
Schur's lemma gives

`End_2I(4 direct-sum 6) = C direct-sum C`,

whose unitary automorphism group is `U(1)^2`.  Before coupling to the rest of
the McKay graph there is therefore one relative phase after quotienting a
common scalar.  **DERIVED representation-theoretic ambiguity.**  Calling
this relative phase a gauge field would be premature.

The exact affine-`E8` McKay adjacency is a connected tree.  Requiring a
diagonal node-phase transformation to commute with its adjacency forces the
phases to agree across every edge.  Only the global scalar `U(1)` remains;
after quotienting that scalar, the effective relative torus is trivial.
**DERIVED NO-GO:** the unfluctuated McKay operator does not turn the gluing
ambiguity into a nontrivial gauge `U(1)`.

A nontrivial gauge field would require a specified represented algebra,
opposite action and inner fluctuations rather than the adjacency commutant
alone.  Whether the induced phase modules support such a finite spectral
triple is **OPEN**.

## Canonical spectral-carrier screen

The most direct nontrivial double is

`H+ = Ind(chi_1) = 2+4+6`,

`H- = Ind(chi_3) = 2-prime+4+6`,

with Galois exchange as the candidate real structure.  Unlike a trivial
chirality doubling, its two sheets are inequivalent.  Nevertheless Schur's
lemma permits equivariant off-diagonal Dirac blocks only between the common
`4` and `6` summands.  Every such `D` has rank at most 20 on the
24-dimensional carrier and leaves the `2` and `2-prime` endpoints in a
four-dimensional kernel.  **DERIVED.**  A kernel alone is not fatal, but the
endpoint projectors commute with every such `D`, so the canonical node
geometry is disconnected.  **DERIVED negative.**

For the canonical central node algebra `C^4` on types `(2,2-prime,4,6)`, the
graded dimension vector is `(2,-2,0,0)`.  Galois exchanges the first two
types and fixes the common ones.  The resulting intersection form has rank
two rather than four and determinant zero.  **DERIVED NO-GO:** this natural
24-dimensional phase-induced double fails Poincare duality and is not yet a
viable matter carrier.

Changing the algebra by assigning the common summands to coarser blocks
could alter the result, but no canonical assignment has been derived; doing
so now would be **STRUCTURAL/FITTED**.  Thus the new bridge presently ends at
an exact spinorial McKay module, not a manifold-like finite spectral triple.

## Hopf-harmonic support theorem: the exact `16+14` split

Inducing only the fundamental phase understates the structure.  Exhaustive
Frobenius-reciprocity calculation for all ten characters `chi_q` of `C10`
shows that every induced module has dimension 12 and that orientation pairs
`q` with `-q` at the level of irreducible support.

Take support union, meaning that an irreducible is retained once whenever it
appears, rather than adding repeated multiplicities.  Then

- the even harmonics `q=0,2,4` have support
  `1+3+3-prime+4+5`, of total dimension 16;
- the odd harmonics `q=1,3,5` have support
  `2+2-prime+4-spinorial+6`, of total dimension 14.

These are exactly the two bipartite halves of the affine-`E8` McKay node
module.  Their union contains all nine irreducibles once and has dimension
`16+14=30`.  No node assignment was fitted.  Moreover `-1 in C10` acts on
`chi_q` as `(-1)^q`, so Fourier parity agrees exactly with integer versus
half-integer spin parity.  **DERIVED exact support theorem.**

This supplies a new rigorous bridge:

`Hopf C10 Fourier parity -> McKay chirality gamma_F`.

The hostile boundary is important.  Support union is an idempotent operation
on representation types, not the ordinary direct sum of induced modules.
The three even inductions and three odd inductions each have total dimension
36 before repeated irreducibles are identified.  Why physical state content
should use support union rather than multiplicity is **OPEN/STRUCTURAL**.
Therefore the theorem derives the node inventory and its chirality grading,
but not yet its Hilbert multiplicities or dynamics.

### Canonical multiplicity completion

There is, however, a canonical categorical object behind the support
statement.  Sum all five even harmonics with their actual induction
multiplicities.  The result is

`direct-sum_(q even) Ind(chi_q) = e_+ C[2I]`,

where `e_+=(1+z_c)/2` and `z_c=-1` is the central involution.  Likewise

`direct-sum_(q odd) Ind(chi_q) = e_- C[2I]`.

Each ideal has dimension 60.  In the even ideal every integer-spin irrep
occurs with multiplicity equal to its dimension; in the odd ideal every
spinorial irrep does likewise.  Their direct sum is the full 120-dimensional
regular representation.  **DERIVED.**

Thus support union is not needed to construct the Hilbert modules: it is
simply the operation of forgetting multiplicities in the simple spectra of
the two canonical central ideals.  The `16` and `14` count dimensions of the
simple-spectrum representatives, whereas the canonical Hilbert dimensions
are balanced `60+60`.  **DERIVED clarification.**

This does not yet supply dynamics.  The two central ideals have disjoint
irreducible support, so

`Hom_2I(e_+ C[2I], e_- C[2I]) = 0`.

No nonzero `2I`-equivariant odd Dirac connects them.  McKay adjacency does
connect their simple spectra, but as the tensor-by-defining-spinor functor,
not as an equivariant endomorphism of the regular module.  **DERIVED NO-GO.**
A physical odd `D` must therefore be functorial/categorical, symmetry
covariant rather than invariant, or involve additional bimodule structure.
Choosing among these possibilities is **OPEN**.

### Canonical relational dynamics

Although there is no equivariant odd endomorphism between the two central
ideals, tensoring by the defining spinor `2` is a canonical odd
correspondence.  Exact character arithmetic gives

`2 tensor e_+ C[2I] = 2 copies of e_- C[2I]`,

`2 tensor e_- C[2I] = 2 copies of e_+ C[2I]`.

On the simple spectrum this correspondence is precisely the affine-`E8`
McKay adjacency matrix `A`.  The dimension vector

`d=(1,2,2,3,3,4,4,5,6)`

satisfies `A d=2d`.  Its Perron normalization

`P_ij=A_ij d_j/(2d_i)`

is therefore an exact stochastic matrix.  Its stationary distribution is
the Plancherel law

`pi_i=d_i^2/120`.

Every transition changes central/Hopf parity, while the stationary weight is
exactly one half on each sector.  **DERIVED.**  This gives a parameter-free
categorical dynamics in which chirality is alternated by the fundamental
interaction rather than selected as an independent binary state.

Calling this Markov correspondence physical time evolution or a Dirac
operator would be incorrect: it is stochastic functorial data on simple
representation types.  Promoting it to quantum amplitudes, a spectral
correspondence, or an unbounded Kasparov cycle is **STRUCTURAL/OPEN**.

### Spectral closure back to Fibonacci

The normalized McKay dynamics closes exactly back onto the original
two-state Fibonacci process.  Its characteristic polynomial is

`t (t^2-1) (t^2-1/4) (t^4-(3/4)t^2+1/16)`,

and its nine eigenvalues are the values of the defining `2I` character
divided by two:

`1,-1,0,+-1/2,+-phi/2,+-1/(2phi)`.

In particular, the two Galois-conjugate golden modes are `phi/2` and
`sigma(phi)/2`.  Their projective ratio is

`(sigma(phi)/2)/(phi/2)=sigma(phi)/phi=-phi^-2`,

exactly the nonstationary memory eigenvalue of the original Fibonacci
Perron process.  **DERIVED spectral closure.**

Thus the chain now closes algebraically:

`Fibonacci memory -> C10 cyclotomic phase -> 2I induction -> McKay dynamics`

`-> projective Galois ratio -> Fibonacci memory`.

After two steps the process returns to the same chirality.  Each parity
block has slow nontrivial eigenvalue `phi^2/4`, hence exact spectral gap

`1-phi^2/4=(3-phi)/4`.

These are exact dynamics of the categorical Markov correspondence.
Identifying the projective eigenmode ratio with a measured physical decay,
mass ratio, or coupling would be **PATTERN** until an observable map is
constructed.

The exact audit is
`reproducible/verify_fibonacci_nonbinary_dynamics.py`.
