# Result: the generalized negative-mode recurrence closes

Date: 2026-08-18

## Headline

The kinetic--stiffness generalized negative fiber is common across the old and
shifted centered slices and is preserved, within the frozen conservative error
bounds, by both normalized recurrence operators `Gamma` and `Omega`.

The preregistered outcome is

```text
GENERALIZED_MODE_RECURRENCE_CLOSURE_CERTIFIED.
```

This resolves the preceding `B,D` phase-leakage negative in the intended
direction: the Euclidean stiffness projector was not the dynamical mode
projector; the Hermitian definite pair `(A,B)` selects the closed
configuration-mode fiber without fitting.

## Provenance ledger

| stage | commit |
|---|---|
| generalized-mode prior-art/framing gate | `52ec90a` |
| frozen recurrence-closure protocol | `dfd833e` |
| registered verifier before execution | `c8a8544` |
| deterministic artifact | `ae36c24` |

The verifier source has SHA-256

```text
0a84c8ec4fab1c9626d5e4c711f89c6f9638cf37c15ef6b0050d6b66dfdde6c1
```

and its twice-reproduced byte-identical artifact has SHA-256

```text
53e046e2020a97fc992559546ce3d45479c0c0de7ce2e01322b09998ba85cf80.
```

Both runs reported `7/7` checks passed.  Only targeted verifiers were run; the
full suite was not run.

## Complete result

For both centered slices, both parities, sectors `4,5` and all four derivative
schedules:

```text
32/32 Hermitian-definite pencil projectors resolved as 15 negative + 10 positive,
16/16 old/shifted projectors GENERALIZED_COMMON_FIBER_RESOLVED,
64/64 Gamma/Omega leakages LEAKAGE_ZERO_CONSISTENT.
```

Broken down by time and operator:

| slice | operator | zero-consistent | nonzero/open |
|---|---|---:|---:|
| old | `Gamma` | 16 | 0 |
| old | `Omega` | 16 | 0 |
| shifted | `Gamma` | 16 | 0 |
| shifted | `Omega` | 16 | 0 |

The midpoint leakage norms are about `8.4e-8 ... 2.1e-7`.  The deliberately
conservative propagated leakage bounds are about `0.005 ... 0.014`, so the
residuals lie near `1.5e-5 ... 1.7e-5` error units.  The old/shifted projector
distance is about `2.05e-5` of its combined error bound.

## Numerical caveat

The generalized-projector error is conservative:

```text
old      about 0.0034,
shifted  about 0.0159.
```

Therefore “certified” means that closure is robustly consistent with zero
under the preregistered finite-family propagation.  It is not a formal
analytic interval theorem for the exact derivative, and the bound is too wide
to quantify a tiny nonzero rotation.  A direct high-precision generalized
projector reconstruction is the natural independent tightening.

## Why this object is more dynamical

On the action-selected shape carrier,

```text
A=-V_S,  B=-M_S>0,  A v=lambda B v.
```

This is the standard mass/kinetic--stiffness mode problem.  Its negative
subspace is an invariant spectral subspace of the normalized stiffness part
`Omega`; the nontrivial result is that the drift part `Gamma` also preserves
it on both slices.  Thus the full second-order centered recurrence can be
restricted to the same rank-`15` generalized fiber in each of sectors `4,5`.

Post-result literature confirms the correct next mathematical object.  The
standard quadratic mode problem has the form

```text
(lambda^2 M + lambda D + K)x=0,
```

with mass, damping/drift and stiffness matrices; see Tisseur and Meerbergen,
[*The Quadratic Eigenvalue Problem*](https://doi.org/10.1137/S0036144500381988),
and the [Netlib eigenvalue templates](https://netlib.org/utk/people/JackDongarra/etemplates/node19.html).
For Regge gravity, Höhn's [*Canonical linearized Regge Calculus: counting
lattice gravitons with Pachner moves*](https://arxiv.org/abs/1411.5672)
shows that identifying propagating lattice-gravity degrees of freedom also
requires gauge/constraint control.  None of these sources establishes the
present finite 600-cell result.

## Status ledger

- **DERIVED COMPUTATIONAL, conditional on the frozen derivative family:** the
  generalized rank-`15` fiber in each target sector is common across two
  centered slices and closed under `Gamma,Omega`.
- **DERIVED COMPUTATIONAL NEGATIVE:** the Euclidean stiffness fiber's naive
  phase lift is not closed; the kinetic metric is load-bearing.
- **STRUCTURAL:** the closed generalized fiber is the correct finite
  mass/kinetic--stiffness modal object.
- **OPEN:** a direct high-precision generalized-projector confirmation,
  longer-time persistence, constraint/gauge quotient and refinement.
- **OPEN:** the reduced quadratic/discrete eigenfrequencies, stability,
  dispersion, effective speed and any mass gap.
- **OPEN:** identification with gravitons, particle inertia or Standard-Model
  mass.

## Next load-bearing calculation

Restrict the complete centered recurrence to the now selected rank-`30`
generalized carrier and solve its quadratic/discrete characteristic problem
without fitting.  Report whether its roots are oscillatory, growing/decaying
or unit-modulus, and only then compare their spatial symmetry labels with the
600-cell Laplacian spectrum.  A mass or limiting speed claim is forbidden
until that comparison and a refinement test exist.
