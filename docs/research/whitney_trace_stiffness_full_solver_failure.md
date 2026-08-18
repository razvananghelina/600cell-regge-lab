# Recorded solver failure: complete trace-stiffness preregistered LOBPCG

Date: 2026-08-11

Protocol commit: `7ed6d49`

## Status

The preregistered generalized block-LOBPCG calculation returned **10/11**.
This is a numerical-method failure, not a physical negative and not an
accepted spectral result.

The dense control calibration passed all six degrees/levels, with maximum
relative value error below (8\times10^{-15}).  Five of the six complete
base/refined degree blocks also passed the frozen (10^{-7}) recomputed Ritz
gate.

The refined degree-two block failed because its largest returned eigenpair
had relative Ritz residual

\[
 6.56\times10^{-6},
\]

above the frozen gate.  LOBPCG reached its preregistered 2,000-iteration
ceiling for that largest block.

## Consequence

The provisional complete ratios

```text
(4.3793, 3.2485, 3.8092)
```

and their provisional spread `1.3481` are **NOT ACCEPTED** at this commit.
They are recorded only so the failed attempt cannot disappear from
provenance.

## Non-target-driven repair available

Degree two has an exact structural simplification absent in degrees zero and
one: each triangle occurrence graph has one row and rank one, so the complete
row image is the entire row carrier.  The generalized quotient can therefore
be transformed directly into an ordinary symmetric positive sparse matrix

\[
 H^{1/2} R M^{-1}R^*H^{1/2}
\]

of dimension 28,800.  A symmetric Lanczos solver can certify its extremal
blocks without changing the eigenvalue definition or the (10^{-7}) Ritz
gate.

Using that solver is a post-failure method correction.  It must be committed
before recomputation and cannot be described as part of the original blind
protocol.

## Ledger

- **DERIVED:** dense-control quotient representation is correct.
- **DERIVED NUMERICAL:** five complete degree/level blocks pass.
- **FAILED GATE:** refined degree-two maximum Ritz residual.
- **OPEN:** corrected degree-two extremal certificate.
- **NOT A RESULT:** provisional complete ratios and balance verdict.
