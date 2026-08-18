# Result: the full generalized cotangent lift is not transported

Date: 2026-08-18

## Status ledger

| Statement | Status |
|---|---|
| The generalized rank-15 configuration fiber changes between the two middle slices, conditional on the frozen four derivative schedules | **DERIVED COMPUTATIONAL** (upstream) |
| The full phase lift `E direct-sum E*` is transported by the action tangent | **DERIVED COMPUTATIONAL REFUTATION** |
| A separate Cholesky/principal-angle implementation reproduces the refutation | **STRUCTURAL INDEPENDENT CORROBORATION** |
| A nonzero constraint-reduced transported intersection exists | **OPEN** |
| That intersection is a graph or Lagrangian subspace selected by the action | **OPEN** |
| A propagator, dispersion relation, particle inertia, mass or limiting speed follows | **OPEN / NOT COMPUTED** |
| External novelty | **OPEN** |

## 1. Frozen construction

For every parity, symmetry sector `4,5` and all four derivative schedules, the
frozen phase fibers were

```text
F_old     = E_old     direct-sum E_old*,
F_shifted = E_shifted direct-sum E_shifted*,
Q_t       = diag(P_t, conjugate(P_t)).
```

The unique second-slab tangent `T_2=[A B;C D]` was reconstructed from the
Regge--dust action.  The preregistered condition was

```text
(I-Q_shifted) T_2 Q_old = 0.
```

No fitted graph, alignment or momentum relation was used.

## 2. Execution chronology and repair

- prior-art/framing commit: `60aafe1`;
- preregistered protocol commit: `3794418`;
- initially registered verifier commit: `63cfa6b`;
- permanent adversarial-replication rule: `ef7f1fe`;
- documented precision-context repair: `b47f5c9`;
- accepted high-precision artifact commit: `86f6b00`.

The first execution at `63cfa6b` correctly terminated `6/7` with
`GENERALIZED_PHASE_TRANSPORT_CONTROL_FAILED`.  Importing the frozen residual
verifier after raising the ambient `mpmath` precision changed the last digits
of its legacy `1e-70` error floor and therefore changed its JSON hash.  The
distances and classifications did not change, but that run is diagnostic only.

The failure and cause were recorded before rerunning.  The repair makes the
legacy 15-digit initialization context explicit during replay, restores the
caller's 100-digit context afterwards and creates the new verifier's own
`1e-70` only after fixing its precision.  No scientific equation, threshold or
outcome branch changed.

After a standalone `10/10` reconstruction restored the exact frozen residual
hash

```text
3244185127aecf7c9a44261cced0be521c9dc42bf8e44f909d8a0ce10a96eadf,
```

the phase verifier ran twice at `7/7` with byte-identical artifact SHA-256

```text
45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743.
```

No full suite was run.

## 3. High-precision result

All 16 exact second-slab tangent balls passed the determinant, boundary-map,
principal-identity and symplecticity controls.  The complete census gave

```text
A blocks: 16/16 LEAKAGE_NONZERO_RESOLVED
B blocks: 16/16 LEAKAGE_NONZERO_RESOLVED
C blocks: 16/16 LEAKAGE_NONZERO_RESOLVED
D blocks: 16/16 LEAKAGE_NONZERO_RESOLVED
full maps: 16/16 LEAKAGE_NONZERO_RESOLVED
```

The block residual norms span approximately

```text
3.4512e-5 ... 3.8759e1,
```

against complete errors

```text
1.7389e-54 ... 2.3004e-47.
```

Even the least separated block lies `3.76e45` error units from zero.  The full
phase residual is approximately `38.75939`, at `1.68e48...2.02e48` error
units.  The result is therefore not a threshold or finite-precision marginal.

Accepted mechanical outcome:

```text
GENERALIZED_PHASE_TRANSPORT_REFUTED
```

## 4. Adversarial independent replication

The project-wide rule added in commit `ef7f1fe` requires more than rerunning
one implementation.  A separate post-result audit was therefore fixed and
executed:

- independence/prior-art commit: `f071366`;
- adversarial protocol commit: `8f0e7f6`;
- registered verifier commit: `0534f27`;
- accepted adversarial artifact commit: `24d0ee9`.

It did not use the residual-certified projectors, the four-block leakage
formula or the new tangent reconstruction as its decisive route.  It instead

1. rebuilt the generalized bases from the earlier direct `M,V` midpoints;
2. used an explicitly coded Cholesky whitening and `numpy.linalg.eigh`;
3. loaded the earlier committed two-step tangent archive;
4. QR-orthonormalized `T_2 F_old` and measured direct sine-of-principal-angle
   leakage into `F_shifted`'s orthogonal complement.

Both audit executions passed `8/8` and produced byte-identical SHA-256

```text
c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880.
```

All 16 canonical cells were `INDEPENDENT_NONCLOSING`.  Their ranges were

```text
maximum principal sine       0.9813993142476888 ... 0.9813993142476927
relative least-squares leak  0.7006578008167527 ... 0.7006578008167836
```

The synthetic exact-transport control gave `7.9e-16...1.1e-15`; the known
one-direction leakage control gave `1` to rounding.  Basis rephasing/reversal
changed the result by at most `1.9e-15`, and the direct principal-sine and
projector-distance calculations agreed within `3.4e-16`.  The alternative
same-representation and swapped-dual convention stresses also failed in all
16 cells, though they are not part of the decisive definition.

Accepted adversarial outcome:

```text
ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED
```

This float64 audit is not a second exact certificate.  It is the mechanically
independent falsification attempt required before accepting the high-precision
certificate.

## 5. What the negative does and does not mean

**DERIVED COMPUTATIONAL.** The complete 30-dimensional cotangent lift of the
generalized configuration fiber is not an invariant time-dependent phase
bundle under the canonical second tick.  Therefore it cannot itself be used as
the reduced propagating carrier.

It does **not** follow that the Regge model has no dynamics, no perturbations or
no wave-like sector.  Full cotangent invariance was deliberately a sufficient
and strong condition.  The action may propagate a smaller subspace selected by
pre/post constraints or a Lagrangian relation.

The post-result literature search makes that distinction load-bearing.
Marsden and West derive symplectic evolution from discrete variational
principles (<https://doi.org/10.1017/S096249290100006X>).  Dittrich and Hoehn
show that singular discrete systems evolve between pre- and post-constraint
surfaces by a pre-symplectic map, and that the propagating degrees of freedom
depend on both endpoints (<https://doi.org/10.1063/1.4818895>,
<https://arxiv.org/abs/1303.4294>).  Generalized Hamilton--Jacobi treatments
describe the evolution graph as a Lagrangian relation between endpoint phase
spaces (<https://doi.org/10.3934/jgm.2022010>).

These sources do not predict our result.  They show why failure of the full
lift logically points to a constrained relation rather than closing the whole
dynamics programme.

## 6. Next honest gate

The next canonical object is fixed without fitting:

```text
K_old = F_old intersection T_2^-1(F_shifted).
```

The next mission must compute its exact rank in every cell before inspecting a
desired dimension.  If `K_old=0` throughout, the present generalized-fiber
route is closed.  If it is nonzero, only then test whether it is symplectic,
isotropic/Lagrangian, a graph over configuration variables, stable across the
four derivative schedules and transported consistently over the next slab.

No graph may be selected by minimizing leakage after seeing this negative.
Mass, inertia, dispersion and a limiting speed remain downstream and **OPEN**.

