# Verdict: the embedded self-stress carrier is not dynamically closed

Date: 2026-08-18

## Status ledger

| statement | status | evidence |
|---|---|---|
| The regular embedded 600-cell rigidity map has rank `470`. | **KNOWN / CONTROL** | Whiteley's convex-polytopal rigidity theorem; reproduced in both schedules. |
| Its tangent-vertex restriction has rank `354`. | **KNOWN / CONTROL** | Only the six ambient rotations remain in the tangent kernel. |
| The conformal and tangent images intersect in dimension `4`. | **KNOWN / CONTROL** | `120 + 354 - 470 = 4`; reproduced numerically with no open singular values. |
| `ker R^T` has dimension `250`. | **KNOWN / CONTROL** | Framework self-stress count `720 - 470`; it is not a graviton count. |
| The frozen binary-tetrahedral sectors restore the full `470/250` split. | **DERIVED COMPUTATIONAL** | All seven minimal-sector ranks are resolved and restore the global theorem count in both schedules. |
| The self-stress carrier is preserved by the centered Regge recurrence. | **DERIVED COMPUTATIONAL NEGATIVE** | Refuted: 408 of 448 required cross/leakage quantities are resolved nonzero. |
| This carrier is a physical transverse-traceless/York sector. | **OPEN** | No action-derived constraint quotient or DeWitt-orthogonal York projector selects it. |
| The calculation is externally novel. | **OPEN** | A targeted literature search found no identical construction, but this is not a novelty proof. |

## Frozen question and provenance

The complete hypothesis was:

1. the carrier is the `720` logarithmic squared-edge variations of the fixed,
   regular, embedded 600-cell;
2. `R` is the literal Euclidean rigidity differential of its `120` vertices in
   `R^4`;
3. the candidate complementary carrier is exactly the Euclidean edge-orthogonal
   self-stress space `S = ker R^T`;
4. the dynamics is the already committed centered dust-Regge recurrence, tested
   in every minimal binary-tetrahedral sector, both parity schedules and all four
   frozen derivative variants;
5. no post-result rotation, fitted coefficient, continuum spectrum, polarization
   count or desired wave speed is loaded.

Commit ordering:

| stage | commit |
|---|---|
| prior-art and framing gate | `a318d6e` |
| preregistered protocol, before dynamic inspection | `f9b692e` |
| registered verifier | `b9d4de2` |
| first certified artifact | `2cb3015` |

The result artifact is
`reproducible/gravity_600cell_dust_rigidity_york.json`, with SHA-256

```text
251851c08f81ba2f0d41c2d0da428ab11f1ba918b9cb59e0a1e347143c883981
```

Two completed targeted runs produced this identical byte sequence.  The second
run ended with `9/9 PASS` and
`RIGIDITY_YORK_DECOUPLING_REFUTED`.  No full-suite run was performed.

## Geometry controls

For an edge `uv`, the normalized rigidity differential is

```text
(R z)_uv = 2 (x_u-x_v) dot (z_u-z_v) / ell^2.
```

Radial displacements `z_v = sigma_v x_v` reproduce the previously certified
conformal map exactly, `R J = C`.  The two schedule matrices are related by the
literal edge permutation: their operator residual is `1.335e-15`, while the
common-edge squared-length spread is `5.602e-11` within the frozen arithmetic
envelope.

Both schedules resolve

```text
rank R                = 470
rank(R P_tangent)     = 354
rank(C, R P_tangent)  = 470
dim(im C intersect im D) = 4
dim ker R^T           = 250
```

The high-precision symmetry bases have maximum residual `1.5452e-98`.
Minimal-sector ranks and self-stresses are identical for both schedules:

| irrep dimension | rank of `R_d` | minimal self-stress | weighted self-stress |
|---:|---:|---:|---:|
| 3 | 59 | 31 | 93 |
| 2 | 40 | 20 | 40 |
| 2 | 40 | 20 | 40 |
| 2 | 38 | 22 | 44 |
| 1 | 20 | 10 | 10 |
| 1 | 20 | 10 | 10 |
| 1 (trivial) | 17 | 13 | 13 |
| **restored total** | **470** | — | **250** |

These numbers validate the carrier construction.  They are not evidence for
gravity because the global counts follow from rigidity theory before the Regge
operators are loaded.

## Primary dynamic result

For `M,N,V`, both bilinear cross blocks between `im R` and `S` were tested.  For
the normalized recurrence matrices `Gamma=M^-1 N` and `Omega=M^-1 V`, leakage
out of `S` was tested.  The `448` preregistered classifications are:

| operator | resolved nonzero | open | zero-consistent |
|---|---:|---:|---:|
| `M` | 112 | 0 | 0 |
| `N` | 96 | 16 | 0 |
| `V` | 96 | 16 | 0 |
| `Gamma` | 48 | 0 | 8 |
| `Omega` | 56 | 0 | 0 |
| **total** | **408** | **32** | **8** |

The weakest resolved nonzero `Omega` leakage is still `295.43` times its full
error envelope.  Every `M` cross block is resolved nonzero by at least
`3.61e6` error units.  The eight zero-consistent `Gamma` entries occur only in
the trivial sector and cannot rescue global closure.  All `56` schedule
comparisons are `SCHEDULE_ROBUST`; their largest distance/error ratio is
`3.30e-10`.

Therefore:

> **DERIVED COMPUTATIONAL / STRUCTURAL NEGATIVE.**  On the fixed coarse
> 600-cell and the committed centered dust-Regge recurrence, the canonical
> embedded framework self-stress carrier is not dynamically invariant.  The
> proposed route from framework self-stress to a York/tensor sector is closed.

The `32` open quantities do not weaken this refutation: certification of
closure required every quantity to be zero-consistent, while hundreds are
separated from zero.  This is a universal-condition failure, not a majority
vote.

## Attack on the framing

Calling `ker R^T` a "York sector" without qualification is too strong.  In
continuum general relativity the York split is selected by the spatial metric,
its covariant derivative, the DeWitt supermetric and the Hamiltonian and
momentum constraints.  Here `ker R^T` is selected by the Euclidean bar-and-joint
embedding.  Framework self-stress means a left null vector of an edge-length
differential; it is not matter stress-energy and it is not automatically a
transverse-traceless metric perturbation.

The negative result is consequently narrow but clean.  It rules out this
particular geometry-selected proxy.  It does **not** prove that the Regge model
has no tensor perturbations, no gravitons or no continuum limit.

There is also no exact coarse action-derived gauge quotient available to repair
the result: the already certified boundary Legendre cross block is invertible.
Discarding the `354` tangent directions by hand would therefore manufacture a
constraint absent from the fixed curved action.

## Post-result prior-art check

The technical terms exposed by the calculation were searched again after the
outcome: `Regge rigidity matrix`, `Regge equilibrium self-stress`, `self-stress
invariant subspace`, and `curved Regge pseudo-constraints`.

- Bahr and Dittrich show that curved Regge solutions generally lose exact
  discrete diffeomorphism gauge symmetries and replace constraints by
  pseudo-constraints: <https://arxiv.org/abs/0905.1670>.
- Hoehn's lattice-graviton construction uses action-derived vertex-displacement
  constraints in linearized Regge calculus around a **flat** background, not an
  embedded framework self-stress complement on a curved dust `S^3`:
  <https://arxiv.org/abs/1411.5672>.
- The exact pre/post-constraint framework for simplicial evolution is developed
  by Dittrich and Hoehn: <https://arxiv.org/abs/1108.1974>.

No primary source located in this targeted search performs the same 600-cell
self-stress/dynamic-closure comparison.  That absence is only a search result;
external novelty remains **OPEN** pending a dedicated review.

## What remains physically load-bearing

The next legitimate question is not another imposed quotient on this fixed
carrier.  It is whether an **action-derived** weak/pseudo-constraint subspace
emerges and approaches the continuum vertex-displacement constraints under a
declared refinement family.  Such a test must examine the small singular
values and vectors of the refined Legendre/Hessian blocks, preregister their
scaling law, and compare them to geometric tangent displacements without
forcing equality.

Until that refinement result exists:

- exact coarse gauge symmetry: **DERIVED NEGATIVE**;
- canonical conformal minority-sign carrier: **DERIVED COMPUTATIONAL**;
- framework self-stress as a closed tensor carrier: **DERIVED NEGATIVE**;
- physical two-polarization graviton sector: **OPEN**;
- wave equation and an effective limiting speed: **OPEN**.
