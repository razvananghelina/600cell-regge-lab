# Preserved first result: static Lorentz-covariance control failure

Date: 2026-08-19

## Frozen provenance

| item | value |
|---|---|
| prior-art commit | `403059c` |
| protocol commit | `f2dc6d3` |
| registered-verifier commit | `7a7e411` |
| verifier SHA-256 | `308d97cc0b057d3ac79cbc8a4706a63fe1b3d76792d84d8576a84df7d7d63514` |
| artifact SHA-256 | `3ac5cce9db2b2f828e0ced2114f301f761dd9371847b712ad47119709396cf7d` |

The first registered execution returned

```text
12/13
RELATIVE_POINCARE_CONTROL_FAILED.
```

This result is preserved rather than overwritten.

## What passed

Exact rational arithmetic gave

```text
det(T) = -256 tau (lambda-1)^3.
```

The direct fixed-bottom kernel and the Poincare-parameter kernel agreed
exactly.  At the three preregistered representatives the untransformed
classification was

```text
(lambda,tau)   rank(T)  pure translations  Lorentz projection
(1,5)             1            3                   3
(2,5)             4            0                   6
(3,11)            4            0                   6.
```

Origin-shift and paired-vertex-permutation controls passed.  The rational
Lorentz-boost control passed at both expanding representatives and failed at
the static representative.

## Why no scientific conclusion is accepted yet

The verifier's covariance tuple includes the separate ranks of the
projections onto a preselected rotation basis and a preselected boost basis.
That split depends on the observer's timelike direction.  A Lorentz boost
conjugates the static spatial-rotation stabilizer into a mixed rotation/boost
subalgebra, so those two separate ranks need not be invariant even when the
total Lorentz subspace transforms exactly.

This is a suspected classifier error, not yet a repaired result.  The next
gate must preserve this failed artifact and test separately:

1. the exact intertwining equations for the displacement and constraint
   matrices;
2. invariant ranks: total Lorentz projection and pure-translation kernel;
3. covariance of the static Lorentz image by conjugation, rather than
   equality of its coordinate rotation/boost split;
4. the transformed observer-normal stabilizer interpretation.

Until that correction is preregistered and independently checked, both the
expanding Lorentz chart and the static stabilizer interpretation remain
**OPEN** despite their passing raw ranks.
