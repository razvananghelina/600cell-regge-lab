# The `(2,3,5)` orbifold incidence route: corrected audit and closure

Date: 2026-08-09

## Executive verdict

The preregistered `N=3` result is **refuted**.  Nontrivial twisted incidence
operators do not vanish: there are exactly 62 canonical incidence-map adjoint
pairs, and there is one genuine signed cellular complex

`F0 -> E2 -> V0`.

After correcting the enumeration, the scoped matter route still reaches its
kill boundary.  None of the 62 kernels or cokernels, and not the sole middle
cohomology, is `M15`, `M16`, `3*M15`, or `3*M16`.  The hit rate is

`0/63 canonical objects = 0/125 operator-dependent output slots`.

This is a **DERIVED NEGATIVE** for complex line modules induced from the
derived `C10`, `C4`, and `C6` stabilizers, with maps supported on the single
actual incidence double coset.  It is not a theorem against every conceivable
higher-rank orbifold bundle.

The endomorphism-algebra, real-form, KO6, finite-triple, central-lift, and
anomaly stages were not run: the corrected matter-output gate failed first.

## Provenance ledger

The shared repository contains two original protocol commits:

- original Step 1: `36bd6825b6b64edc6c95c1ac67a21a80693318fc`;
- original Step 3: `d2b02eb5fc1a42428bc06ff67b0cef4c6d2da295`.

They do have the requested parent ordering, but they do not certify the
corrected result.  Step 1 enumerated only three untwisted maps because of a
twisted-kernel bug.  Step 3 compared that incomplete list and, in the same
commit, accidentally added the corrected 62-map data.  The corrected list
therefore did not precede all target comparison in Git history.

The corrective sequence is:

- corrected target-independent census:
  `42c4c709d82005c1372a72a273775e8791c3b4da`;
- corrected comparison:
  `decdff66d41d2ef359eba11a1e20323d6b0e8de5`.

This second sequence freezes code in the requested order, but it cannot
retroactively restore blindness: the target and the earlier comparison were
already known.

### Why local Git preregistration is not proof of blindness

Commit ordering proves tree ancestry, not when a human inspected a target.
Local history can also be rebased, amended, or recreated unless its hash is
published outside the committer's control.  Here the target was already in
the prompt and prior note, and a concurrent process committed files while the
independent audit was running.  Therefore Git ordering is useful provenance,
but not external proof that no target knowledge influenced the list.

A stronger protocol would:

1. put the enumeration in a clean repository that contains no target module;
2. publish the signed enumeration commit hash to an append-only remote log or
   independent timestamp service;
3. have a different party release the target after that publication;
4. run the comparison in CI from the published tree without allowing history
   rewriting.

This provenance defect would invalidate a claimed positive.  It does not
create the corrected negative, which follows from an exhaustive list and a
target-independent central-parity obstruction.

## Hypotheses and acceptance test

All claims below assume:

1. `G=2I` acts on the icosahedral vertex, edge, and face sets;
2. their binary stabilizers are `C10`, `C4`, and `C6`;
3. the fibers are one-dimensional complex characters of those cyclic groups;
4. a primitive canonical map is the unique-up-to-scale equivariant kernel
   supported on the single actual cell-incidence orbit;
5. the reverse map is its Hermitian adjoint;
6. a short complex is a normalized `F -> E -> V` incidence pair whose
   composition vanishes exactly;
7. arbitrary sums of double-coset kernels, arbitrary Schur coefficients, and
   target-selected direct sums are forbidden.

The corrected acceptance test is operator-dependent.  A route advances only
if `ker`, `coker`, or the middle cohomology equals `M15`, `M16`, or three
identical generation modules.  A virtual index does not count because

`[ker T]-[coker T]=[source]-[target]`

for every linear map.  Likewise, an Euler character is endpoint data; the
middle cohomology is the signal.

## 1. Independent geometric reconstruction

The exact quaternion action over `Q(phi)` gives twelve vertex rays.  Their
normalized mutual dot products are

`-1, -sqrt(5)/5, +sqrt(5)/5`.

Taking maximal-dot-product pairs gives 30 edges, and taking graph triangles
gives 20 faces.  This rebuilds the combinatorial icosahedron without rotation
eigenvectors.

The three pair spaces decompose as follows:

| pair | double cosets | size of each group double coset | pair-orbit incidence census |
|---|---:|---:|---|
| vertex--edge | 6 | 20 | one orbit has 60; five have 0 |
| vertex--face | 4 | 30 | one orbit has 60; three have 0 |
| edge--face | 10 | 12 | one orbit has 60; nine have 0 |

Every conjugate-stabilizer intersection is exactly the center `{+1,-1}`.
Thus incidence really is one pure double coset in each case.  **DERIVED.**

The reported integers `(0,3,6)` are not geometric invariants.  With the same
right-coset ordering but exact oriented quaternion axes, the labels are
`(0,1,8)`.  Flipping only the order-six axis to its antipode reproduces
`(0,3,6)`.  A rotation eigensolver is free to make precisely that sign flip.
The label changes; the one-pure-incidence-orbit statement does not.

The right-coset invariant `H*x*y^-1*K`, rather than `H*x^-1*y*K`, is confirmed.

## 2. Why `N=3` is false

The old coefficient for a relative element `d` summed over all `(h,k)` in
`H x K` while testing whether `h*d*k` remained in `H d K`.  It always does.
The sum therefore factorized into full stabilizer-character sums and killed
every nontrivial character.

The Mackey orbit is `(H x K)/L_d`, not `H x K`, where

`L_d = H intersect d K d^-1`.

A covariant line-valued kernel exists exactly when the two characters agree
on `L_d`, and then the orbit-supported Hom summand is one-dimensional.  Here
`L_d=C2` on every cross-cell double coset, so existence is exactly matching
central parity.  The corrected count, made before any corrected comparison,
is

- edge to vertex: `4*10/2 = 20`;
- face to vertex: `6*10/2 = 30`;
- face to edge: `6*4/2 = 12`.

Hence **62**, not 3, nonzero incidence-map adjoint pairs exist.  **DERIVED
refutation.**

## 3. Exact Hom and rank census

For the twenty induced modules, the exact Gram matrix of equivariant Hom
dimensions has diagonal histogram

`3:8, 4:2, 7:4, 8:2, 15:2, 16:2`

and ordered off-diagonal histogram

`0:200, 2:32, 3:8, 4:60, 6:48, 7:4, 10:24, 14:2, 15:2`.

There are zero one-dimensional full Hom spaces.  Geometry nevertheless
selects one one-dimensional Mackey summand by requiring support only on the
incidence orbit.  This defends the chosen canonicity criterion; it does not
select a stabilizer character.

For every one of the 62 maps, the exact source/target multiplicities bound
each irreducible channel rank from above.  Reduction in both finite fields
`F_601` and `F_1801`, using a primitive 60th root, attains that upper bound.
This sandwiches every characteristic-zero rank exactly.  The complete table
of source, target, full Hom dimension, double-coset support, rank, kernel,
cokernel, and virtual index is in
`reproducible/orbifold_incidence_preregistered.json`.

The old three untwisted ranks `12,12,20` are confirmed.  Their apparent
`4/4'` discrepancy is only an irrep-order convention: after translating to
the repository order, the untwisted outputs agree.  The missing 59 twisted
maps are the fatal error in the old census.

Kernel dimensions across the corrected list are

`0,1,8,9,10,18,19`,

and cokernel dimensions are

`0,1,2,10,11`.

## 4. Short complexes

All sixty parity-compatible normalized pairs `F -> E -> V` were multiplied
exactly in `Z[z]/Phi_60(z)`.  Exactly one composition vanishes:

`F0 -> E2 -> V0`.

The `E2` character is the edge-orientation sign.  This is the ordinary
oriented cellular complex of the icosahedral sphere, with

`H2=rho1, H1=0, H0=rho1`.

Thus the old assertion that no three-term complex exists is **DERIVED FALSE**.
Its unsigned `E0` incidence matrices naturally had nonzero composition; they
were not the cellular boundary maps.

There are 62 map/adjoint pairs plus one complex/adjoint pair, so the
look-elsewhere count is `N=63`.  This is a moderate finite family, not a
single rigid operator.  Its 63 endpoint indices have 28 distinct characters.

## 5. Target provenance and the rank-eleven warning

The target is independently derived from the diagonal homomorphism

`g -> (rho5(g),rho2(g)) in SU(3) x SU(2)`.

Exact character multiplication gives `rho5 tensor rho2=rho9`, hence

`M16|2I = 2 rho1 + rho2 + 2 rho5 + rho9`

and

`M15|2I = rho1 + rho2 + 2 rho5 + rho9`.

Their dimensions are 16 and 15.  The central-parity vector of the nine irreps
is

`(+,-,-,+,+,+,-,+,-)`,

so both targets mix central parities; `M16` has central trace zero.
**DERIVED.**

The `20 x 9` all-induction matrix has rank nine.  Since the displayed virtual
identity is an integral solution, all virtual solutions form an affine
integral lattice of rank

`20-9=11`.

Thus the old virtual identity has zero selection weight.  This is stronger
than a generic Artin-induction warning.  **DERIVED.**

## 6. Corrected comparison and kill boundary

The relevant outputs are the 124 separate kernels/cokernels and the sole
middle cohomology, not the 63 virtual endpoint indices.

| comparison module | canonical objects hit | output slots hit |
|---|---:|---:|
| `M15` | 0/63 | 0/125 |
| `M16` | 0/63 | 0/125 |
| `3*M15` | 0/63 | 0/125 |
| `3*M16` | 0/63 | 0/125 |

There are two independent exact obstructions:

1. actual output dimensions are only `0,1,2,8,9,10,11,18,19`, never
   `15,16,45,48`;
2. every kernel, cokernel, and middle cohomology belongs to one central parity,
   while every comparison module mixes the two parities.

Therefore no canonical primitive incidence operator or admitted short
complex produces the required matter character or chirality.  **DERIVED
NEGATIVE: KILL BOUNDARY REACHED for the stated line-incidence route.**

The virtual-index comparison is also `0/63`, but is recorded only to complete
the preregistered multiset.  It is not evidence about the operators.

## 7. Canonicity and higher-rank scope

“Supported on the incidence coset alone” is a defensible sharp definition for
a primitive incidence map: it is invariant, mechanically testable, and leaves
only overall scale.  Allowing non-incidence double-coset weights immediately
restores a multi-parameter Hom space and is **STRUCTURAL**, not canonical.

The remaining weakness is one level earlier: geometry does not select one of
the 62 character pairs.  The exhaustive count prevents hiding that
look-elsewhere freedom, and the zero result makes it harmless here.

Over the complex numbers every representation of a cyclic stabilizer splits
into line characters.  A generic higher-rank equivariant bundle therefore
adds multiplicities and matrix-valued intertwiners; without an independently
derived tangent, spinor, or other geometric bundle, it recreates the forbidden
Schur-coefficient fitting freedom.  It is worth pursuing only if the bundle
and its fiber representation are fixed before a matter comparison.  Generic
higher rank is not a surviving continuation of this route.

## Complete status ledger

- **DERIVED:** exact combinatorial icosahedron and pure incidence double cosets.
- **DERIVED:** old `(0,3,6)` labels are one axis-sign convention; labels are
  not invariant.
- **DERIVED:** Hom Gram census, zero full Hom spaces of dimension one, and
  relation-lattice rank eleven.
- **DERIVED NEGATIVE:** the old `N=3` enumeration and its twisted-vanishing
  claim.
- **DERIVED:** 62 nonzero canonical incidence-map adjoint pairs; exact kernels
  and cokernels.
- **DERIVED NEGATIVE:** the old “no short complex” claim.
- **DERIVED:** exactly one short complex, with `H1=0`.
- **DERIVED:** the diagonal `M15/M16` comparison characters.
- **DERIVED NEGATIVE:** `0/63` canonical-object hits and `0/125` output-slot
  hits; the scoped route is closed.
- **STRUCTURAL:** incidence-only support as canonicity; exclusion of arbitrary
  orbit sums and target-selected direct sums.
- **PATTERN:** none promoted.  Virtual identities are explicitly assigned zero
  evidential weight.
- **OPEN OUTSIDE SCOPE:** an independently selected higher-rank geometric
  bundle.  Arbitrary higher rank is forbidden fitting.
- **NOT REACHED AFTER KILL:** index endomorphism algebra, selected real form of
  `C+H+M3(C)`, KO6, orientability, Poincare duality, order zero, first order,
  connectedness, central lift, and generation-blind anomaly equations.

Exact verifiers:

- `reproducible/verify_incidence_operator_enumeration.py`;
- `reproducible/verify_orbifold_incidence_route.py`.

No PDF build was performed.
