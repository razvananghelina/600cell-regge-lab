# Edge-space Krajewski audit: exact no-go for the canonical two-sided route

Date: 2026-07-22

## Decision

The bipartite edge observation is correct and important, but it does not
produce the Standard-Model finite algebra.  For the full two-sided finite
group, every arrow is a distinct irreducible outer product

`Hom(V_s,V_t) = V_s* boxtimes V_t`

of `2I x 2I`.  Hence

`End_(2I x 2I)(E) = C^16`.

**DERIVED (negative):** this commutant contains neither `H` nor `M3(C)`.
The earlier large commutant `M16+M6+M16+M22` belongs to the *diagonal* `2I`
restriction.  Passing to the genuinely two-sided action destroys, rather than
selects, those multiplicities.

More decisively, the canonical endpoint bimodule and orientation grading have
no nonzero first-order odd Dirac.  Krajewski-legal blocks share a source or a
target.  On a bipartite graph such blocks have the same orientation and the
same grading.  Odd blocks reverse orientation, and then share neither the
left nor right endpoint.  Therefore

`first order + D gamma = -gamma D  ==>  D=0`.

**DERIVED no-go:** the edge-space construction dies at the combined
first-order/evenness gate.  No commutant generator `Y`, anomaly forcing, or
`M15/M16 x N_gen` census is licensed.

All finite statements are checked by `reproducible/verify_edge_matter.py`.

**2026-07-27 scope correction (DERIVED).**  The zero-Dirac theorem below is
for the canonical endpoint bimodule and orientation grading.  The free-arena
Q8 counterexample proves that it cannot be extrapolated to arbitrary
multiplicity-mixing opposite actions on another Hilbert space.

## 1. Exact edge census and two-sided commutant

Use the McKay-chain convention in the question:

| node | irrep | dimension | parity |
|---|---:|---:|---|
| `rho0` | `1` | 1 | integer |
| `rho1` | `2` | 2 | spinor |
| `rho2` | `3` | 3 | integer |
| `rho3` | `4s` | 4 | spinor |
| `rho4` | `5` | 5 | integer |
| `rho5` | `6` | 6 | spinor |
| `rho6` | `4` | 4 | integer |
| `rho7` | `2'` | 2 | spinor |
| `rho8` | `3'` | 3 | integer |

The 16 oriented blocks are:

| arrow | `2I x 2I` irrep | dimension |
|---|---|---:|
| `rho0 -> rho1` | `rho0* boxtimes rho1` | 2 |
| `rho1 -> rho0` | `rho1* boxtimes rho0` | 2 |
| `rho1 -> rho2` | `rho1* boxtimes rho2` | 6 |
| `rho2 -> rho1` | `rho2* boxtimes rho1` | 6 |
| `rho2 -> rho3` | `rho2* boxtimes rho3` | 12 |
| `rho3 -> rho2` | `rho3* boxtimes rho2` | 12 |
| `rho3 -> rho4` | `rho3* boxtimes rho4` | 20 |
| `rho4 -> rho3` | `rho4* boxtimes rho3` | 20 |
| `rho4 -> rho5` | `rho4* boxtimes rho5` | 30 |
| `rho5 -> rho4` | `rho5* boxtimes rho4` | 30 |
| `rho5 -> rho6` | `rho5* boxtimes rho6` | 24 |
| `rho6 -> rho5` | `rho6* boxtimes rho5` | 24 |
| `rho6 -> rho7` | `rho6* boxtimes rho7` | 8 |
| `rho7 -> rho6` | `rho7* boxtimes rho6` | 8 |
| `rho5 -> rho8` | `rho5* boxtimes rho8` | 18 |
| `rho8 -> rho5` | `rho8* boxtimes rho5` | 18 |

Their dimensions sum to 240.  Each outer product is irreducible for the direct
product, and ordered pairs of irreducible labels classify these products.
All 16 pairs are different.  Schur's lemma therefore proves the stated
`C^16` commutant.  **DERIVED.**

### What the seed embeddings do—and do not do

`rho0=1`, `rho1=2`, `rho8=3'`, and the Galois twin `rho2=3` supply the stated
ambient scalar, weak, and two color embeddings.  Their images act *on* the
relevant representation factors.  They are not elements of the commutant of
those actions.  Asking their commutant to contain the same noncommutative
algebra reverses the representation/commutant roles.  **DERIVED scope
correction.**

The only arrows whose two endpoints both have one of the available ambient
actions are the two orientations of `0--1`, `1--2`, `2--3`, and `5--8`, where
the weak action on `rho3` is `Sym^3(2)` and the last edge uses the factorization
below.  Thus there are 8, not a generation census of
uniform color/weak blocks.  **DERIVED.**  The number has no SM interpretation.

If one makes the **STRUCTURAL** minimal ambient-extension choice—`rho3` as the
weak quartet and `rho5` as the weak sextet—the eight defined arrows have
product-group types

| unoriented edge | forward / reverse types | multiplicities |
|---|---|---:|
| `0--1` | `(1,2)`, `(1,2)` | 2 |
| `1--2` | `(bar3,2)`, `(3,2)` for the Galois-color embedding | 1, 1 |
| `2--3` | `(bar3,4)`, `(3,4)` | 1, 1 |
| `5--8` | `(bar3',6)`, `(3',6)` | 1, 1 |

Thus the conditional commutant on this 76-dimensional defined subspace is
`M2(C) + C^6`; its only multiplicity larger than one is the two orientations
of the colorless weak doublet.  **DERIVED conditional on the displayed
extensions.**  The other 164 dimensions have no supplied color action, so
there is no commutant of a globally defined derived SM subgroup on all `E`.
In particular, neither `M3` nor `H` appears as the requested commutant.

There is an additional ambiguity at `rho5=6`: the same finite `2I` irrep can
be extended either as the weak sextet `Sym^5(2)` or, using
`6=2 tensor 3'`, as the `(3',2)` product module.  These ambient extensions are
inequivalent and the finite character does not select between them.
**DERIVED ambiguity / OPEN selector.**  The table uses the parity-motivated
weak-sextet choice; choosing the product extension changes the ambient block
decomposition and confirms that it is not canonical.

## 2. Exact tensor factorizations and their limits

Exact character inner products give

`2 tensor 3' = 6`,

`2 tensor 3 = 2 + 4s`,

`2 tensor 4 = 6 + 2'`.

Also `4s=Sym^3(2)`, while `2'` is the Galois-twin defining spinor.
**DERIVED.**

Consequences:

- `rho5=6` can carry the external `SU(2) x SU(3)` action associated with
  `rho1 tensor rho8`, up to the multiplicity-one intertwiner.  The isomorphism
  class is derived; a normalized CG map is not selected.  **STRUCTURAL choice.**
- `rho3=4s` carries the weak action `Sym^3(2)`, but it is only a summand of
  `2 tensor 3`; it is not a color-times-weak module.  Projecting to that
  summand breaks the external product action.  **DERIVED (negative).**
- `rho7=2'` occurs only together with `6` in `2 tensor 4`; it is not by itself
  a color-times-weak module.  **DERIVED (negative).**
- `rho4=5` does carry `Sym^4(2)` as an `SU(2)` module, but as an integer node
  it has no supplied color action; `rho6=4` likewise has no supplied color
  action.  A finite `A5` action alone does not extend them to arbitrary
  `SU(3)` actions.  **OPEN if new structure is added; unavailable here.**

The Galois pair `3,3'` therefore gives two inequivalent real `A5` embeddings,
not color fundamental versus antifundamental: both `3` and `bar3` restrict
identically along either real embedding.  **DERIVED (negative), consistent
with `galois_doubling_triple.md`.**  The proposed orientation resolution does
not occur.

## 3. Canonical `gamma` and `J`

Let `gamma=+1` on `Hom(even,odd)` and `-1` on `Hom(odd,even)`.  Give every
block its Hilbert--Schmidt inner product and define

`J(phi:s->t) = phi^dagger:t->s`.

Adjoint is anti-linear, isometric, and involutive.  It reverses orientation.
Thus

`J^2=+1`, `J gamma=-gamma J`.  **DERIVED.**

For any self-adjoint `D` satisfying `J D=D J`, the remaining KO6 sign is
`JD=DJ`.  In particular `D=0` obeys it, so the exact sign table is

`(J^2, JD/DJ, Jgamma/gammaJ)=(+,+,-)`.

**DERIVED but degenerate:** the signs exist; the Dirac gate below admits only
zero for the canonical endpoint algebra.  Composing `J` with Galois requires
chosen isometric intertwiners and adds no escape.  **STRUCTURAL variant.**

## 4. Dirac candidates and the killing constraint

The canonical endpoint algebra is the vertex algebra, acting on an arrow
`(s,t)` from its left/source and right/target indices.  The first-order
Krajewski rule permits a `D` block only between arrows sharing a left index or
sharing a right index.

Because the graph is bipartite:

- arrows with a common source have the same orientation;
- arrows with a common target have the same orientation;
- arrows of opposite orientation share neither corresponding index.

Therefore every first-order-legal block commutes with `gamma`, while every
odd block is illegal.  Their intersection is zero.  **DERIVED decisive
no-go.**

The proposed canonical candidates do not evade it:

1. Composing two edge arrows produces a length-two path.  Its endpoints have
   the same parity (or it backtracks into `End(V_s)`), so it is never another
   edge block.  Projection back to `E` is exactly zero.  **DERIVED.**
2. A fixed-edge multiplication is not an endomorphism of `E` before that
   zero projection; it maps degree-one paths to degree two.  **DERIVED
   negative.**
3. The preprojective moment map `sum_e [phi_e,phi_ebar]` is quadratic and
   vertex-`End`-valued.  It is a constraint on quiver representations, not a
   linear Dirac operator on `E`.  **DERIVED scope correction.**
4. Arbitrary CG maps could manufacture other linear blocks, but with the full
   endpoint algebra they still face the legality/oddness theorem.  Shrinking
   or replacing the algebra could change the test, but no derived SM
   subalgebra has been selected.  **OPEN outside this route.**

## 5. Hypercharge and matter census

There is no surviving nonzero first-order odd `D` and no selected
`C+H+M3(C)` representation.  Hence there is no derived commutant `u(1)`
generator to require to be doublet-constant and generation-blind.  Route C
cannot be applied.  **DERIVED scope limit.**

Accordingly, the exact comparison with `M15/M16 x N_gen` is: **undefined as a
module comparison, not a numerical near-match**.  The 240-dimensional total
alone is `16 x 15`, but the 16 inequivalent arrow blocks have dimensions

`2,2,6,6,12,12,20,20,30,30,24,24,8,8,18,18`,

not uniform 15- or 16-dimensional generations.  Reading `240=16*15` as a
matter census would be forcing.  **PATTERN rejected.**

## Status ledger

### Strengthened

- **DERIVED:** complete 16-arrow census and dimension 240.
- **DERIVED:** full two-sided commutant is exactly `C^16`.
- **DERIVED:** exact tensor factorizations and their extension limits.
- **DERIVED:** canonical adjoint reversal gives the KO6 `J`/`gamma` signs.
- **DERIVED no-go:** endpoint first order plus orientation oddness forces
  `D=0`.
- **DERIVED:** path multiplication projects to zero on the edge space; the
  preprojective moment map is not a Dirac operator.

### Downgraded

- The diagonal commutant's abstract `M3` and `H` embeddings are not selected
  by restoring the two-sided action; the two-sided commutant is smaller.
- `3` versus `3'` is Galois/outer twisting, not `3` versus `bar3`.
- `240=16*15` is dimension arithmetic only.

### Open

- a different, independently derived quotient/subalgebra of the endpoint
  action for which a nonzero odd first-order Dirac survives;
- a selected associative `C+H+M3(C)` representation;
- color orientation, generation blocks, and a generation-blind `Y`.

Those require new structure.  They are not hidden in the canonical edge/Hom
space described here.
