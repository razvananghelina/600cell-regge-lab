# Protocol: exact rational global lift of canonical boundary data

Date: 2026-08-19

Freeze this target-disclosed protocol before constructing the first exact
rational pivot system or seeing any lift coefficient.

## Frozen inputs

- prior-art/framing gate SHA-256
  `5a1a695a4526a70320f9316771f015d3c14bffdf3bca3d115a52fa8e9fe73d27`;
- modular data-space verifier SHA-256
  `21022ab2014f5c95bd3e2f06bf1137f713533df465d152a40c069f107dd153ae`;
- its first frozen artifact SHA-256
  `3db0b9ce8c90cba9de3fbbff818129388d79a98e0483a0ca3ae53b2e4d271434`;
- complete admissibility source and artifact SHA-256
  `4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e`
  and
  `fa45c80739ca0dda4f82c9da98a4b22f4d8a18c182a40696a2a22d1d26ec89a1`;
- refuted local-lift source and artifact SHA-256
  `3adb80448e19fd99f0b8ec205497f11325b0a6f7a72c9a2785f9b65778707750`
  and
  `4065950aaac4180ec1cdd0b82f7a8bc403b2969c50d26cf14cc28592085cb2c5`.

## Complete hypotheses and target

For the exact rational complete face equations

```text
F f + E e + S s = 0,
```

use 3600 cell-flex variables, 720 upper spatial squared-edge data, and 120
strut squared-length data.  Define `U` independently from the sorted
600-cell graph by

```text
(U sigma)_{uv} = 8 lambda (sigma_u + sigma_v).
```

The disclosed candidate is the 240-column rational data block

```text
B = [E U, S].
```

The target is an exact rational matrix `L` satisfying

```text
F L + B = 0
```

on every original face equation.  No coefficient may be selected, fitted, or
discarded after a residual is seen.

Run the complete test for:

1. both baseline rational representatives `(lambda,tau)=(2,5),(3,11)`;
2. both representatives with the alternate exact local right inverse;
3. reversed face orientation at `(2,5)`;
4. an odd canonical relabelling at `(2,5)`;
5. reversed metric sign at `(2,5)`.

The lift coordinates may differ between right-inverse graphs and conventions;
existence, uniqueness in the declared flex coordinates, and zero physical
residual must agree.

## Mechanically independent exact method

Do not call the modular-rank routine and do not infer rational inclusion from
the two finite fields.

1. Reconstruct `F` and `B` directly from the exact face-equation dictionaries.
2. Perform deterministic sparse Gaussian elimination over `sympy.Rational`,
   choosing the least remaining flex column as pivot.  Determine pivots from
   the flex block, while carrying all 240 right-hand-side columns.
3. Require 3600 flex pivots.  Derive `L` by exact reverse substitution through
   those pivot equations.
4. The solution is determined using the pivot equations only.  Independently
   substitute it into all original equations, including every non-pivot row,
   and require every one of the 240 coefficients of every residual to vanish
   exactly.
5. Canonically serialize every nonzero coefficient of `L` as
   `(flex,data,numerator,denominator)` and record its SHA-256, nonzero count,
   support census, and maximum numerator/denominator bit lengths.  The source
   must regenerate the complete matrix; the artifact need not duplicate all
   coefficients.

Because `F` has exact full column rank, a successful lift is unique in the
declared cell-flex coordinates.  This uniqueness is kinematic and is not a
claim of physical uniqueness.

## Negative construction attack

Let `{u0,v0}` be the lexicographically first sorted edge.  Replace only its
incidence row `(u0:1,v0:1)` by `(u0:1)`, retaining the same factor
`8 lambda`, to form `U_bad` and `B_bad=[E U_bad,S]`.

Require exact rational ranks

```text
rank(U)=120, rank(U_bad)=120, rank([U U_bad])=121.
```

Reduce `B_bad` through the same rational flex pivots.  At least one exact
nonzero zero-flex obstruction must remain in every construction.  Record the
first obstruction.  If both candidate and corrupted block lift, declare a
control failure because the already frozen edge-compatible space has
dimension 120 modulo both primes.

## Required controls

1. Every frozen input hash and upstream outcome reproduces byte-for-byte.
2. The exact 600-cell has `f=(120,720,1200,600)`, its graph is connected and
   contains a triangle, and the three incidence ranks are `120,120,121`.
3. Every local geometry, transition, inverse, and face control remains true.
4. Every exact flex elimination has precisely 3600 pivots.
5. Candidate consistency is decided by exact zero rows, not a numerical
   tolerance or a modular proxy.
6. A positive lift passes direct substitution on all original rows.
7. The corrupted image leaves an exact nonzero obstruction.
8. Record lift support rather than assuming it is local.  Whether the support
   is local or global is descriptive, not an acceptance criterion.
9. Reproduce, but do not overwrite or reinterpret, the frozen negative status
   of the earlier local lift.

## Outcome hierarchy

- `RATIONAL_DATA_LIFT_CONTROL_FAILED`: provenance, combinatorics, exact
  pivot, direct-residual, corruption, or upstream controls fail.
- `RATIONAL_DATA_LIFT_DISAGREEMENT_OPEN`: legitimate constructions disagree
  on existence of the candidate lift.
- `RATIONAL_VERTEX_STRUT_DATA_LIFT_REFUTED`: all controlled constructions
  agree that an exact candidate obstruction remains.
- `RATIONAL_VERTEX_STRUT_DATA_LIFT_DERIVED`: every construction has 3600
  exact pivots, the candidate has an exact unique lift with zero direct
  residual on every equation, and the corrupted image is rejected.

The positive branch establishes only the rational first-order kinematic
carrier and its cell-flex response.  The action Hessian, symplectic form,
constraint/gauge split, physical lapse, tick, tensor modes, `c`, `G`, and
Planck units remain **OPEN**.

## Reproducibility discipline

Register and commit the verifier before its first execution.  Run only this
targeted verifier.  Freeze its first artifact before constructing or testing
the boundary action/Hessian.

