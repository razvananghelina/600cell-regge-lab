# Oriented barycentric chamber double

Date: 2026-07-28

## Decision

The 120 complete barycentric flags of the icosahedron produce the first
repository-derived carrier supporting an explicit finite real-even spectral
triple with every audited finite-geometry gate:

- real-even with KO6 signs `(+,+,-)`;
- orientable in metric dimension zero;
- Poincare-dual, with a unimodular rank-120 intersection form;
- equipped with a derived nonzero odd chamber-adjacency operator;
- equipped with nonzero one-forms;
- first order;
- connectedness.

The algebra is an explicitly certified partition algebra `C^36`.  Its
selection is **STRUCTURAL/FITTED**, not derived: a deterministic combinatorial
search selected one 90-bit vacuum among symmetry-breaking choices.  Thus a
finite spectral triple is constructed, but a derived matter model and any
Standard-Model interpretation are not.

All claims are checked locally by
`reproducible/verify_oriented_chamber_double.py`.

## 1. Derived double

A chamber is a complete flag

`vertex < edge < face`.

There are `20*3*2=120` chambers.  The rotation group `A5` acts freely and
has exactly two chamber orbits of size 60.  These are the two orientations.
A geometric reflection is a fixed-point-free involution exchanging them.

Let `gamma` be `+1` on one orbit and `-1` on the other, and let `J` be
reflection followed by coefficient conjugation.  Then

`J^2=1`, `J gamma=-gamma J`.

This is not a trivial inserted copy: the sheets are the two derived chamber
orbits.  **DERIVED.**

## 2. Dirac and Poincare pairing

Two barycentric chambers are adjacent when their flags differ in exactly
one entry.  The resulting chamber graph is 3-regular.  Every adjacency
crosses orientation, hence its adjacency operator `D` satisfies

`D gamma=-gamma D`.

Reflection preserves adjacency, hence

`JD=DJ`.

For the chamber-function algebra `A=C^120`, minimal projections are chamber
projectors.  The opposite projector of chamber `j` is the projector of its
reflected chamber.  Therefore

`cap_ij=gamma_i delta_(i,Jj)`.

The intersection form is an antisymmetric signed permutation matrix.  It
has rank 120 and determinant of absolute value one.  **DERIVED.**

Orientability also holds: `gamma` is itself a represented chamber function,
so `gamma=pi(gamma)Jpi(1)J^-1` is a Hochschild 0-cycle.  **DERIVED.**

## 3. Remaining obstruction

Both the represented and opposite chamber-function algebras are diagonal,
so order zero holds.  Full independent chamber projectors fail first order
against chamber adjacency.  The verifier prints an exact minimal-projector
witness.  One-forms are nonzero.  **DERIVED negative.**

Thus this carrier passes precisely the axioms killed by trivial doubling,
but it does not yet pass the bilinear Dirac gate.

### Exact partition-algebra search

For an algebra of functions constant on a partition `P` of the chambers,
first order is equivalent to the following edge condition:

> for every chamber edge `{x,y}`, either `x` and `y` lie in one `P` block,
> or `Jx` and `Jy` lie in one `P` block.

This turns the first-order problem into a finite combinatorial problem.
**DERIVED.**

The 180 chamber edges form three `A5` orbits of size 60, one for each
changed flag coordinate.  Central inversion `J` preserves every orbit.
Exhausting all `2^3` `A5`-invariant contraction choices gives only one
first-order survivor: the scalar algebra.  Its fluctuations vanish.
**DERIVED no-go for `A5`-invariant partition algebras.**

After allowing symmetry breaking, a deterministic search plus local rank
optimization found an explicit partition with

`A=C^32`,

for which:

- first order holds exhaustively;
- inner one-forms are nonzero;
- the `32 x 32` intersection form has determinant `1`.

This is an **EXISTS/PATTERN** witness, not a derived vacuum: its 90-bit
partition choice is search-selected.

It still fails orientability.  Sixteen of its 92 nonempty intersections
between an algebra block and an opposite-algebra block contain both
orientations.  A subsequent 250-step local optimization did not find a
complete survivor.  This is evidence, not a no-go.  **OPEN.**

That first witness was then used as the starting point for an orientability
search.  It found a second explicit partition with

`A=C^36`

for which:

- all first-order partition constraints hold;
- all common algebra/opposite-algebra intersections are orientation-pure,
  so orientability holds;
- the first complete witness had exact rank 36 and determinant 4;
- inner one-forms are nonzero;
- the quotient graph is connected, hence `[D,a]=0` only for scalars;
- the representation is faithful and unital.

A subsequent search restricted entirely to already legal, orientable,
full-rank partitions found an integral witness with

`det(cap)=+1`.

Its 90-bit certificate is printed as `INTEGRAL_SURVIVOR_BITS` and registered
in the verifier.  Hence the pairing is unimodular over `Z`, not merely
nondegenerate over `Q`.  Every property is rechecked from the certificate,
independently of the optimization score.  **DERIVED verification of an
EXISTS/STRUCTURAL construction.**

The integral partition has trivial stabilizer in `A5`; its rotation orbit
therefore contains exactly 60 equivalent partitions.  Its opposite `J(P)`
does not belong to that orbit, so the full icosahedral orbit has size 120.
The selected vacuum breaks the full `A5` rotation symmetry and comes with a
distinct opposite orbit.  **DERIVED; particle/antiparticle language would
be PATTERN.**

For this exact witness:

- `dim_C Omega_D^1=128`, counted by directed quotient-edge blocks;
- the unitary algebra is `U(1)^36`;
- the kernel of `u -> uJuJ^-1` is the global `U(1)`;
- the effective gauge torus therefore has Lie dimension 35.

These are **DERIVED algebraic outputs**.  They are not Standard-Model gauge
fields: the group is a large abelian torus.

## Local non-uniqueness and a first selection functional

The integral witness is not isolated.  Exhausting every partition at
Hamming distance one or two from it gives:

- 17 legal orientable full-rank partitions;
- 15 integral (`|det cap|=1`) partitions;
- none of those 15 lies in the original witness's `A5` orbit;
- all 15 have algebra `C^36`;
- their one-form dimensions are `126` (one case), `128` (eight cases), and
  `130` (six cases).

Thus determinant one and algebra dimension 36 do not select a unique
vacuum.  **DERIVED local non-uniqueness.**

Within this completely enumerated radius-two neighborhood, the functional

`A -> dim_C Omega_D^1(A)`

has a unique minimum, `126`.  The corresponding integral certificate is
printed as `MINIMAL_OMEGA_SURVIVOR_BITS`; it also has trivial `A5`
stabilizer.  **DERIVED local minimum; PATTERN as a proposed global
variational selection principle.**

The later exact `C5`-symmetric contraction witness has integral Poincare
duality and `dim_C Omega_D^1=110`.  Therefore 126 is **REFUTED** as the
global minimum; it remains only the certified local minimum in the
radius-two neighborhood described above.  The actual global minimum and any
physical principle requiring its minimization remain **OPEN**.

Two global verification attempts have now been made:

- exact union-find branch-and-bound reached 10,000,000 nodes and 2,220,453
  completed orientable partitions without exhausting the tree;
- a global SMT query asking for `dim Omega <=124` returned `unknown` after a
  600-second timeout.

Those two older searches neither found a lower integral witness nor produced
an exhaustion or UNSAT certificate.  The later `C5` search did find the
exact 110 counterexample.  Consequently 126 is only a certified local
minimum, and 110 is the current certified global upper bound.  See
`chamber_global_minimum_audit.md`.  **DERIVED counterexample; global minimum
OPEN.**

Numerically, `126` is exactly the number of roots of `E7`.  The search did
not target this number: it minimized the one-form dimension among the local
integral survivors.  This equality is recorded as a **PATTERN only**.  The
actual gauge group of the witness is the abelian torus described above, not
`E7`, so no exceptional gauge symmetry is claimed.

There is a sharper exact count, but also a decisive negative.  The minimal
quotient graph has 36 vertices and 63 undirected edges.  Independently,
`E6` has 36 positive roots, while `E7` has 63 positive roots; under the
standard `E6` root subsystem the latter split in cardinality as

`63 = 36 + 27`.

These root counts and the quotient counts are **DERIVED**.  Identifying the
36 quotient vertices with positive `E6` roots, the 63 quotient edges with
positive `E7` roots, or the excess 27 with the fundamental `27` of `E6`
remains a **PATTERN**: no canonical bijection or compatible action has been
constructed.  In fact the most direct proposed identification is
**REFUTED**: the Hasse graph of the positive-root poset of `E6` has 60
edges, whereas the quotient graph has 63, so the two graphs are not
isomorphic.

The exact new boundary is therefore:

- the full algebra: orientable and PD, but not first order;
- the fitted `C^32` algebra: first order, PD and fluctuating, but not
  orientable;
- the fitted `C^36` integral algebra: every audited finite spectral-triple
  axiom passes, with `det(cap)=1`;
- the `A5`-invariant first-order algebra: only scalars.

## 4. Galois correction

The orientation swap is not the nontrivial outer automorphism of `A5`.
The full geometric icosahedral symmetry group is

`A5 x C2`,

with the second factor generated by central inversion.  Every improper
geometric symmetry is central inversion times a rotation, so its
conjugation on `A5` is inner.  Consequently:

- **REFUTED:** geometric chamber reflection realizes the Galois outer
  automorphism;
- **DERIVED:** geometric reflection exchanges the orientation sheets;
- **OPEN:** any separate arithmetic action relating the chamber carrier to
  the Galois character permutation.

## 5. Innovation boundary

The chamber set is naturally the regular set of the full order-120
icosahedral Coxeter symmetry.  Its 3-regular adjacency is the chamber/Coxeter
Cayley graph, and orientation is the determinant character.  This suggests
that the correct algebra may be a chamber groupoid, Hecke algebra, or
one-sided convolution algebra rather than all chamber functions.

The next exact problem is:

> classify geometrically canonical subalgebras of `C^120`, or canonical
> convolution/groupoid actions on the same chamber carrier, for which first
> order holds while the rank-120 pairing and nonzero one-forms survive.

No algebra reduction is selected yet.  **OPEN.**

## Status ledger

### DERIVED

- two oriented free `A5` chamber orbits;
- geometric reflection `J`;
- chamber orientation grading;
- 3-regular odd chamber Dirac;
- KO6 signs;
- orientability;
- unimodular rank-120 Poincare pairing;
- order zero and nonzero one-forms;
- failure of first order for `C^120`.

### STRUCTURAL

- the symmetry-breaking `C^36` partition algebra;
- promoting that existence witness to a physical internal algebra.

### PATTERN

- interpreting chamber orientation as particle/antiparticle chirality.

### OPEN

- a geometric/dynamical rule selecting the certified `C^36` partition, or
  another member of its solution set, without search;
- relation to arithmetic Galois doubling;
- any particle, generation, gauge, or Yukawa interpretation.
