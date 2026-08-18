# Preregistered protocol: canonical one-step time-slab carrier

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- NO LORENTZIAN DYNAMICS CLAIMED**

## 1. Question and complete scope

Let `K` be exactly the fixed simplicial boundary of the 600-cell, with
`f(K)=(120,720,1200,600)` and its full `H4` action. The previous audit gave a
canonical cotangent arena for its 720 edge lengths but no selected kinetic
Hamiltonian.

This audit asks the prior carrier question:

> Does `K`, functorially and without a chosen vertex ordering, supply a
> canonical four-dimensional one-step carrier with two spatial boundary
> slices; and does any already certified 115,200-state local walk coincide
> with the chamber incidence of that carrier?

The audit is exhaustive only in the following explicitly frozen elementary
class:

1. the cone `C K`;
2. the cellular product cylinder `K x I`;
3. a simplicial triangulation of `K x I` using only its product vertices;
4. the barycentric subdivision `sd(K x I)`;
5. one-layer elementary 4D Regge/Pachner moves on the initial slice;
6. the committed robust eight-component walk under a projection-preserving
   identification `(spatial chamber,component) -> (same spatial chamber,
   product flag label)`.

It is **not** an exhaustive no-go against arbitrary cobordisms, extra
vertices, spin foams, causal sets or future non-factorized identifications.

## 2. External structural benchmark

Canonical simplicial gravity uses the action/Hamilton principal function as
the generator of discrete evolution, with 4-simplex gluings inducing Pachner
moves on the hypersurface. This is external methodology, not a theorem of the
600-cell model:

- Dittrich--Hoehn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>;
- Hoehn, *Canonical linearized Regge Calculus*,
  <https://arxiv.org/abs/1411.5672>;
- Bahr--Dittrich, *(Broken) Gauge Symmetries and Constraints in Regge
  Calculus*, <https://arxiv.org/abs/0905.1670>.

The last source is load-bearing for the warning that curved finite Regge
backgrounds need not possess exact continuum-like first-class constraints.
No external result may be relabelled as derived by this theory.

## 3. Frozen combinatorial checks

The verifier must rebuild `K` and check the following without a physical
target comparison.

### Cone and cylinder

1. `C K` has the cone `f`-vector and Euler characteristic one, with exactly
   one `K` boundary. It is a cap/filling, not a two-slice evolution carrier.
2. The product CW complex `K x I` has two boundary copies of `K`, 600
   tetrahedral-prism four-cells, its exact product-cell counts and Euler
   characteristic zero.
3. Its 120 vertical edges form one `H4` orbit. This supplies a natural
   vertex-lapse carrier but does not fix its metric value or make it a
   multiplier.

### Simplicial canonicity

4. The stabilizer of one tetrahedron acts as the full `S4` on its vertices.
5. On every square `edge x I`, an endpoint transposition exchanges the two
   diagonals. Since a vertex-preserving simplicial triangulation must select
   exactly one, no such triangulation can be invariant under full `H4`.
6. Barycentric subdivision is functorial and evades that diagonal choice by
   adding all cell barycentres. The verifier must count its maximal
   four-flags as

```text
600 prisms * 24 spatial flags * 2 interval endpoints * 4 shuffles.
```

The total must be recorded before comparing it with any existing carrier.

### Elementary Pachner layer

7. Count legal initial `1-4`, `2-3`, `3-2` and `4-1` moves using the exact
   local incidence criteria.
8. Check the `H4` orbit count for each nonempty move class.
9. Build the conflict graph of legal `2-3` moves, where moves conflict when
   their two-tetrahedron supports overlap. If the class is one transitive
   orbit and the full orbit is not independent, no nonempty `H4`-invariant
   parallel `2-3` layer exists. This conclusion is scoped to one parallel
   layer; orbit-complete schedules and quantum sums remain open.

## 4. Frozen 115,200-state comparison

The barycentric chambers of one tetrahedral prism have the canonical label

```text
(epsilon,r),  epsilon in {0,1},  r in {0,1,2,3},
```

where `epsilon` is the starting endpoint of `I` and `r` is the position at
which the temporal promotion is interleaved among the three spatial
promotions. Thus a projection-preserving identification of the robust walk's
eight components with product flags has exactly `8!` possible bijections.

For each of all `8!` bijections, the verifier must test the literal committed
robust stages and macro map against the exact chamber adjacency of
`sd(K x I)`:

- local product-flag adjacency changes `r` by one, or flips `epsilon` only at
  `r=0`;
- an adjacency crossing spatial chamber colour `i` preserves `(epsilon,r)`
  and is available exactly when `r != i`;
- identity/idling is permitted, but every nontrivial transition must be one
  chamber edge.

Record the number of bijections passing each stage separately, all stages,
and the macro map. Do not inspect a preferred labeling and do not fit a
component permutation.

Before registration, an exploratory count noticed the equality
`14,400*8=600*192=115,200` and the fact that the robust cross-chamber rules
also change component labels. Therefore neither the cardinality equality nor
an anticipated factorized failure is blind evidence. The exhaustive table is
the certified content.

## 5. Decision boundary

- **DERIVED CANONICAL CW SLAB:** `K x I` passes the exact two-boundary and
  symmetry checks.
- **DERIVED SIMPLICIAL CHOICE OBSTRUCTION:** no full-`H4`, no-new-vertex
  simplicial triangulation exists, while `sd(K x I)` is canonical only after
  adding barycentric vertices.
- **DERIVED EXISTING-WALK SLAB BRIDGE:** at least one projection-preserving
  bijection makes all three robust stages and their macro incidence-local on
  `sd(K x I)`, and the passing bijection is selected uniquely up to interval
  reversal.
- **DERIVED CURRENT SLAB-DYNAMICS GAP:** the canonical barycentric carrier
  exists but no such bijection passes; equal state counts then carry no
  incidence evidence.
- **OPEN/INCOMPLETE:** any action, Lorentzian metric, lapse constraint or
  exact finite permutation cannot be certified.

Even a positive carrier result does not select a Regge action or prove a
Hamiltonian constraint. A topological cylinder is not yet physical time.

Only the targeted verifier and a static registry check may run. No full suite
and no PDF build.
