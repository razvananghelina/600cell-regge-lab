# Derived bimodule audit: Krajewski and commutant no-go

Date: 2026-07-22 (session with Claude)

## Decision

The McKay node-space route is closed under the explicit derived compatibility
condition below.  There is no compatible Krajewski multiplicity matrix on
either the 30-dimensional node space `W` or the 60-dimensional Galois double
`W + W^sigma`.  The obstruction occurs before the Dirac-legality test: the
allowed restricted bimodule blocks do not contain all nine `2I` irreducibles.

Independently, the largest algebra supplied by the doubled node module is

`End_2I(W + W^sigma) = product_(rho in Irr(2I)) M2(C)`.

It cannot contain `M3(C)`.  Imposing the available grading and the conditional
real structure only makes this commutant smaller.  Thus the node space cannot
carry the requested algebra as a subalgebra of its derived equivariant
endomorphisms.

The 240-dimensional bidirected edge/Hom module has a different outcome: its
commutant is

`M16(C) + M6(C) + M16(C) + M22(C)`

and admits abstract faithful representations of `C + H + M3(C)`.  This is a
**STRUCTURAL possibility**, not a derived finite triple: no edge grading,
real structure, first-order Dirac, or canonical embedding selects it.

All finite counts and algebraic assertions in this note are checked by
`reproducible/verify_bimodule_krajewski.py`.

**2026-07-27 scope correction (DERIVED).**  The later free-arena Q8
counterexample shows that order-zero opposite actions may mix a regular
factor with multiplicity space.  It does not alter the calculations below,
which concern the explicitly stated diagonal node restriction and its
commutant.  This note is not a no-go for arbitrary multiplicity-mixing
antiunitaries on `C[2I] tensor C^m`.

## 1. Compatibility condition

Let the five complex Krajewski labels be

`(1, 1bar, 2, 3, 3bar)`

with dimensions `(1,1,2,3,3)`.  This keeps conjugate orientations distinct at
the algebra level even though the finite subgroup cannot always distinguish
them.

The compatibility condition used in the exhaustive search is:

1. the weak defining module restricts along the derived binary-icosahedral
   embedding as `rho_2`;
2. the color defining module restricts along the derived real color embedding
   as `rho_5 = 3'`;
3. both scalar orientations restrict as `rho_1`;
4. a Krajewski block `(i,j) = n_i tensor overline(m_j)` is a `2I` submodule
   whose restriction is exactly `r_i tensor r_j*`;
5. on the Galois sheet the entire restriction is outer-twisted, exchanging
   `rho_2 <-> rho_3` and `rho_4 <-> rho_5`.

Items 1--3 are **DERIVED** from the existing weak and color embeddings.
Item 4 is **STRUCTURAL and necessary for this route**: it says the same
diagonal `2I` action organizes both the proposed algebra bimodule and the
McKay module.  Without it, an arbitrary vector-space reblocking would not be
a derivation of the algebra action from `2I`.  Item 5 is the strongest
Galois-compatible version currently available; it gives the doubled route
every block allowed by either sheet.

This result does not claim that no unrelated, externally assigned
representation of the Standard-Model algebra can exist on a 30- or
60-dimensional vector space.  It proves that it is not obtained by this
derived diagonal restriction.

## 2. Exhaustive Krajewski search

### Dimension-only search

For a nonnegative `5 x 5` multiplicity matrix `mu`, its Hilbert dimension is

`sum_(i,j) mu_ij dim(n_i) dim(m_j)`.

The generating-function coefficient calculation exhaustively gives:

- **DERIVED:** 188,908,396 matrices of total dimension 30;
- **DERIVED:** 1,362,811,872,984 matrices of total dimension 60.

Consequently, dimension arithmetic alone does not kill the construction.  In
particular, the presence of a five-dimensional McKay node is not by itself a
no-go: a reblocking could cut across nodes if the representation content
allowed it.

### Exact restricted block types

Because `3'` is real on `A5`, conjugating the color label does not change its
finite-subgroup restriction.  The only fusion types are

`1 tensor r = r`,

`2 tensor 2 = rho_1 + rho_4`,

`2 tensor rho_5 = rho_9`,

`rho_5 tensor rho_5 = rho_1 + rho_5 + rho_8`.

Therefore every untwisted compatible block is supported on

`{rho_1, rho_2, rho_4, rho_5, rho_8, rho_9}`.

The required McKay module

`W = rho_1 + rho_2 + rho_3 + rho_4 + rho_5 + rho_6 + rho_7 + rho_8 + rho_9`

also contains `rho_3`, `rho_6`, and `rho_7`, none of which can occur in any
allowed block.  The exact bounded dynamic program over all 25 matrix entries
therefore finds zero compatible matrices.  **DERIVED negative.**

Adding all 25 outer-twisted block types on the second sheet enlarges the total
support to

`{rho_1, rho_2, rho_3, rho_4, rho_5, rho_8, rho_9}`.

It still omits `rho_6` and `rho_7`, while `W + W^sigma` requires multiplicity
two of every irrep.  The exhaustive 50-entry two-sheet search again finds
zero compatible matrices.  **DERIVED negative.**

### Dirac legality

The McKay graph is the connected nine-node, eight-edge affine-E8 tree.
However, no compatible block decomposition reaches the stage at which an
edge can be tested for a shared left or right Krajewski index.  Thus:

- **DERIVED negative:** the killing constraint is `2I` representation
  compatibility, not dimension arithmetic;
- **DERIVED negative:** neither McKay adjacency nor a restricted branch
  Dirac can become a legal Krajewski Dirac on a nonexistent compatible
  diagram;
- **OPEN:** Dirac legality for a different matter space or a different,
  independently derived algebra action.

This is deliberately not reported as a vacuous success of the first-order
condition.

## 3. The maximal node algebra actually supported

### Schur commutants

Every irrep occurs once in `W`, hence

`End_2I(W) = C^9`.  **DERIVED.**

After Galois doubling, every irrep occurs twice (the second occurrence has a
twisted node label), hence

`End_2I(W + W^sigma) = M2(C)^9`,

of complex dimension 36.  **DERIVED.**

With opposite sheet chirality, gamma has one `+` and one `-` line in each
two-dimensional multiplicity space.  Its commutant inside the preceding
algebra is `C^18`.  **DERIVED conditional on the doubled grading.**

For the sheet-swap antiunitary `J = S K` with `J^2=+1`, the `J`-commuting
real algebra is `M2(R)^9`.  Requiring both `J` and gamma leaves

`{diag(a,conjugate(a)) : a in C}^9`,

which is `C^9` as a real algebra.  **DERIVED conditional on the chosen
anti-linear isometric sheet intertwiners.**  Their geometric selection remains
**STRUCTURAL/OPEN**, as in `galois_doubling_triple.md`.

### Fundamental `M3` obstruction

Every nonzero complex representation of the simple algebra `M3(C)` has
dimension at least three.  Thus there is no nonzero algebra homomorphism
`M3(C) -> M2(C)`, and no unital embedding into a product of `M2(C)` factors.

Therefore:

- **DERIVED negative:** the maximal equivariant doubled-node algebra does not
  contain `M3(C)`;
- **DERIVED negative:** its gamma-, `J`-, or Dirac-selected subalgebras cannot
  contain `M3(C)` either, since they are subalgebras of the same maximal
  algebra;
- **DERIVED negative:** `C + H + M3(C)` is not the `J`-real,
  gamma-even part of the node commutant.

This conclusion does not require choosing a node-level lift of adjacency.
That is important because the McKay edges encode tensor-product incidence;
they are not by themselves specified `2I`-equivariant linear maps between
inequivalent summands of the 30-dimensional direct sum.  A commutant with a
particular full block Dirac would require extra Clebsch--Gordan maps and is
presently **OPEN**.  Any such commutant would nevertheless be smaller and
could not evade the `M3` obstruction.

The maximal algebra genuinely supported by the derived node representation
data is therefore `M2(C)^9`; with the conditional KO6 data imposed it reduces
as described above.  It is a multiplicity algebra, not the Standard-Model
finite algebra.

## 4. Secondary edge/Hom check

The existing exact decomposition is

`E_Hom = 16 rho_2 + 6 rho_3 + 16 rho_7 + 22 rho_9`.

Schur's lemma gives

`End_2I(E_Hom) = M16(C) + M6(C) + M16(C) + M22(C)`.  **DERIVED.**

Unlike the node commutant, this algebra is large enough.  For example:

- represent `H` on a multiplicity-16 factor by eight copies of its defining
  complex two-dimensional representation;
- represent `M3(C)` on the multiplicity-6 factor by two defining copies;
- represent `C` faithfully on either remaining factor.

Together with zero action of each direct summand on the other allocated
sectors, this is a faithful unital representation of the real algebra
`C + H + M3(C)` on the edge/Hom multiplicity space.  **STRUCTURAL existence.**

It is not yet derived or distinguished:

- the allocation among large multiplicity factors is not selected by the
  current discrete geometry;
- no derived edge/Hom `gamma` or KO6 antiunitary has been constructed;
- no first-order edge/Hom Dirac has been constructed;
- the real `3` versus `bar(3)` restriction and hypercharge orientation remain
  unresolved.

Thus the edge/Hom space passes only the algebra-size gate.  It is the honest
remaining direction after the node no-go, not a completed matter bimodule.

## 5. Status ledger

### Strengthened

- **DERIVED:** the Krajewski search is finite and exhaustive, including every
  nonnegative multiplicity matrix at dimensions 30 and 60.
- **DERIVED negative:** both compatible solution counts are exactly zero.
- **DERIVED:** the precise killing constraint is missing `2I` fusion support,
  not the five-dimensional node or total dimension.
- **DERIVED negative:** the maximal doubled node commutant cannot contain
  color `M3(C)`.
- **DERIVED:** the edge/Hom commutant has the exact four matrix factors above.

### Downgraded or delimited

- A Krajewski diagram on the node space is not merely unknown under the stated
  compatibility; it is ruled out.
- The doubled KO6 signs do not imply a Standard-Model algebra: their
  simultaneous commutant is still only a product of scalar real algebras.
- The edge/Hom embedding is **STRUCTURAL**, not **DERIVED** or distinguished.

### Open

- a derived edge/Hom Krajewski organization;
- an edge/Hom grading, antiunitary, and legal first-order Dirac;
- a geometric selector for a particular `C + H + M3(C)` embedding;
- color orientation and the generation-blind commutant generator;
- hypercharge, which Route C anomalies can constrain only after that
  bimodule and generator exist.
