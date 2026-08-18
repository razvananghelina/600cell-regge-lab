# Conservation exposes a selection/locality conflict

Date: 2026-08-11

Preregistration commit: `21cbc05`

Targeted verifier:
`reproducible/verify_whitney_projection_selection.py`

Targeted result: **8/8 PASS**.  The verifier is registered exactly once.  No
candidate spectrum and no phenomenological target were used.  The full suite
was not run, by explicit user request.

## Headline

The two natural local broken-FEEC projections cannot be distinguished by
preserving the complete topological sector of the closed three-dimensional
carrier.  Requiring preservation of every conforming Whitney moment does
select a unique projection, but that projection is exactly nonlocal.

> **DERIVED NEGATIVE FOR TOPOLOGICAL-MOMENT UNIQUENESS ON THE FROZEN
> TOWER:** both counting and diagonal-Whitney recovery preserve every
> nonvacuous harmonic moment on `Esd_k(sd boundary Delta^4)`, for
> `k=1,2,4`.

> **DERIVED SELECTION/LOCALITY CONFLICT:** preservation of all conforming
> moments uniquely selects the full-metric orthogonal projection.  Its
> recovery has the exactly nonzero coefficient `243/7480` outside strict
> occurrence support.

This closes the most immediate proposed repair of the projection ambiguity.
Topology is robust, but it contains too little information to choose the
positive-energy dynamics.  The complete Whitney metric contains enough
information, but its inverse propagates globally.

## The algebraic selection theorem

Let `V` carry a positive definite mass `M`, let `W=im(J)`, and suppose
`P^2=P` with `im(P)=W`.  The full-moment requirement is

\[
 \langle Pv,w\rangle_M=\langle v,w\rangle_M
 \quad\text{for every }v\in V,\ w\in W.
\]

It says exactly that `v-Pv` is orthogonal to `W`.  Since `Pv` is in `W`, the
uniqueness of orthogonal decomposition gives

\[
 P=P^A=J(J^TMJ)^{-1}J^TM.
\]

Conversely this projector satisfies the moment identity.  Thus full-moment
preservation is not one heuristic among many: under the stated hypotheses it
is equivalent to the unique full-metric orthogonal projection.  **DERIVED.**

For a local recovery `L`, the verifier evaluates the equivalent defect

\[
 \Delta=L^T(J^TMJ)-MJ.
\]

Already in degree zero, both local candidates have nonzero defects on every
frozen level:

| `k` | global/local dimensions | candidate | nonzero defect entries | maximum absolute defect | relative defect |
|---:|---:|---|---:|---:|---:|
| 1 | `30 / 480` | counting | 5,280 | 0.00462962963 | 0.416666667 |
| 1 | `30 / 480` | diagonal | 5,280 | 0.00462962963 | 0.416666667 |
| 2 | `180 / 3,840` | counting | 50,400 | 0.000578703704 | 0.416666667 |
| 2 | `180 / 3,840` | diagonal | 50,400 | 0.000578703704 | 0.416666667 |
| 4 | `1,320 / 30,720` | counting | 422,880 | 0.0000723379630 | 0.416666667 |
| 4 | `1,320 / 30,720` | diagonal | 422,880 | 0.0000723379630 | 0.416666667 |

The absolute entries shrink with tetrahedron volume, so their raw decrease
is not evidence for convergence.  After normalization by the largest entry
of `MJ`, the observed defect is the same `0.416666667` at all three levels.
This repeated decimal is a **DERIVED NUMERICAL CONTROL**, not promoted here
to an exact all-level value.

The nonlocality conclusion does not rest on a floating threshold.  The new
verifier independently reconstructs the rational `k=1,p=0` mass system.  In
row 5 of the unique orthogonal recovery, local copy 52 belongs to global
vertex 3, so strict occurrence locality requires the coefficient to vanish.
Instead it is exactly

\[
 (L^A)_{5,52}={243\over7480}\ne0.
\]

This reproduces the independent certificate in the earlier canonicity audit.
By the uniqueness theorem, no strict-occurrence-local projector can satisfy
all conforming moment identities.

## Why harmonic conservation cannot select here

Every carrier is a subdivision of a closed triangulated `S^3`, hence its
conforming harmonic dimensions are exactly

\[
 (b_0,b_1,b_2,b_3)=(1,0,0,1).
\]

This exhausts the test rather than sampling harmonic fields:

- `p=0`: the harmonic line is the global constant cochain.  Every top
  tetrahedron has the same exact volume at a fixed level:
  `1/9`, `1/72`, and `1/576` for `k=1,2,4`.  Equal counting is therefore the
  same volume weighting as diagonal scalar Whitney recovery.  Both preserve
  the constant pairing exactly.
- `p=1,2`: there is no harmonic subspace, so the requirement is vacuous.
- `p=3`: every top simplex has exactly one local copy.  Both projections are
  the identity on the entire broken top-form space and hence preserve its
  harmonic line automatically.

The numerical constant-mode residuals range from `2.60e-18` to `6.10e-20`,
but the conclusion uses the exact volume and weight identities, not those
roundoff values.

Therefore harmonic conservation is too coarse for this task for a structural
reason: all nontrivial candidate freedom lives in degrees one and two, while
the closed 3-sphere has no harmonic modes in precisely those degrees.

## Hodge and Poincare do not presently repair the gap

Homological Poincare duality is already reflected in
`(1,0,0,1)`.  Since both complexes have the same cohomology, this level of
duality cannot distinguish their projections.

A metric Hodge-star condition could be stronger, but the present construction
does not derive a dual Hilbert carrier, its metric, or a primal-to-dual star.
The unequal primal dimensions also preclude interpreting the star as a square
endomorphism pairing the existing primal layers.  A barycentric or other dual
cellulation can be introduced, but the associated metric star and its
compatibility with recovery would be additional data that must be
preregistered and justified independently.

> **OPEN:** a derived primal-dual metric Hodge star might impose a new local
> compatibility condition.  This audit does not prove that every such
> construction fails.

It does prove that invoking Poincare duality only at the level currently
present cannot select between the two candidates.

## Attack on the physical framing

Calling the full-moment identity a conservation law is stronger rhetoric than
the derivation supports.  Mathematically it is the variational
characterization of orthogonal projection.  No time translation, Noether
symmetry, or microscopic unitary evolution has derived it.  It is therefore
a **STRUCTURAL projection axiom**, despite the exact theorem that follows if
one assumes it.

Likewise, preserving harmonic pairings protects global topological charges;
it does not determine propagation, dispersion or the positive spectrum.
Passing that gate is necessary for a topologically faithful discretization,
not evidence for physical dynamics.

Combined with the preceding audits, the current situation is a sharp
trilemma:

1. strict local naturality admits at least two distinct projections;
2. complete metric orthogonality selects one projection but makes it global;
3. topological conservation remains local but cannot distinguish the
   projections.

An auxiliary local KKT system can represent the global solve, as established
earlier, but its descriptor metric is singular and its exact constraints are
second class.  It therefore does not yet provide the missing local unitary
tick.

## Status ledger

- **DERIVED:** full conforming-moment preservation uniquely characterizes the
  `M`-orthogonal projection.
- **DERIVED:** both local candidates fail that full condition already in
  degree zero on all three frozen controls.
- **DERIVED:** the selected orthogonal recovery is nonlocal by an exact
  rational coefficient.
- **DERIVED:** both local candidates preserve all harmonic moments available
  on the frozen `S^3` tower.
- **DERIVED NEGATIVE:** topology/harmonic conservation does not select the
  local dynamics.
- **DERIVED NUMERICAL:** the normalized degree-zero full-moment defect is
  `0.416666667` at `k=1,2,4`.
- **STRUCTURAL NEGATIVE:** existing geometric and conservation data do not
  select a unique bounded-star projection.
- **OPEN:** an independently derived local variational or primal-dual axiom.
- **NOT CLAIMED:** a no-go theorem for every discretization, or any derivation
  of Lorentzian time, causal speed, inertia, mass or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_projection_selection.py
```

Expected result: `8/8`.
