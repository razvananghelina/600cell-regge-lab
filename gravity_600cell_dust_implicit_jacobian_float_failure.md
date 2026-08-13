# Recorded first-run outcome: dust implicit Jacobian unresolved

Date: 2026-08-13

Prior-art commit: `31717a8`

Frozen protocol commit: `41acf7b`

Implementation commit: `6621b6a`

Targeted first run: **13/15 implementation checks passed**.  The full suite
was not run.

This file records the failure before any precision correction or altered step
is evaluated.

## 1. Frozen gate that failed

Both schedule parities failed exactly the same preregistered convergence gate:

```text
relative change between h=5e-4 and h=2.5e-4
    even = 4.411880364631418e-4
    odd  = 4.411880563683628e-4

frozen maximum = 3e-4.
```

Every other implementation control passed:

- all 390 displaced points per parity remained Lorentzian and off all branch
  boundaries;
- relative imaginary contamination was at most `1.002e-12`;
- Richardson antisymmetry was at most `1.655e-10`;
- Richardson boundary cross-reciprocity error was at most `5.133e-10`;
- all four independent 60-decimal action curvatures passed, including the
  weakest computed direction, with normalized errors at most `1.062e-6`.

The output is retained in
`reproducible/gravity_600cell_dust_implicit_jacobian.json` and mechanically
labels both parities `NUMERICALLY_UNRESOLVED`.

## 2. Preregistered spectral outcome

For both parities the Richardson spectrum had:

```text
relative rank 1e-7 : 30
relative rank 1e-9 : 30
relative rank 1e-11: 34
boundary rank      : 30 at all three thresholds
combined rank      : 35 at all three thresholds.
```

The thirty large internal singular values lie between approximately `67.15`
and `2402.66`.  They are cleanly separated from five small values:

```text
even: 4.59886e-8, 4.59862e-8, 4.59859e-8, 4.59853e-8, 6.55e-11
odd : 4.59885e-8, 4.59862e-8, 4.59855e-8, 4.56477e-8, 1.51e-10.
```

Because the global Richardson difference has spectral norm `0.3533`, the
smallest singular value is not separated from the frozen empirical envelope.
The preregistered full-rank acceptance boundary is therefore not met.

## 3. Post-result diagnostic, not preregistered evidence

The following inspection used only the already written matrices and is
labelled **PATTERN** until independently frozen and checked.

- The leading `30 x 30` staircase-diagonal block is well conditioned:
  `s_min=67.1531`, condition approximately `35.78`.
- All five small eigenvectors have pole-subspace norm equal to one to the
  printed precision.  Their subspace has all five principal cosines equal to
  one against the five-dimensional pole coordinate space.
- The weakest vector has pole components close to the normalized collective
  vector `(-1,-1,-1,-1,-1)/sqrt(5)` and a staircase-diagonal norm only
  `3.31e-5`.
- Under step halving, its eigenvalue scales approximately as `h^2`:

```text
even: 1.095e-8 -> 2.735e-9 -> 6.367e-10,
```

  while the other four approach approximately `4.6e-8`.
- Eliminating the well-conditioned 30-dimensional block gives a `5 x 5`
  Schur complement with the same small spectrum.

This is exactly the qualitative location expected for lapse/Bianchi or
pseudo-constraint directions, but that interpretation is not yet derived.
The binary64 full matrix cannot decide whether the four `4.6e-8` values are
real broken-gauge stiffnesses and whether the collective mode is exactly
zero.

## 4. Decision

Do not relax the failed `3e-4` threshold and do not call the matrix full rank.
The next correction must leave the 30-dimensional regular block alone and
resolve only the canonically exposed five-dimensional pole Schur complement
using the independent arbitrary-precision action.

Preregister deterministic Schur-lifted pole directions, reconstruct their
quadratic form from action differences, and test whether the four relative
modes remain nonzero while the collective mode tends to zero.  This is a
precision correction to an unresolved result, not a new target search.
