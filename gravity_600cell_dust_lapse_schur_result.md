# Five-dimensional lapse Schur correction: result and framing failure

Date: 2026-08-13

Prior-art update: `0882934`

Frozen correction protocol: `5c0372a`

Implementation commit: `52772c4`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_lapse_schur.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_lapse_schur.json`

Targeted run: **13/13 implementation checks passed**.  The full suite was not
run.

## 1. Mechanical preregistered outcome

Both parities receive the frozen label

```text
FIVE_STIFF.
```

This label is reported honestly and is not silently replaced below.  It was
triggered because all five computed eigenvalues exceeded `100*epsilon_5`,
where `epsilon_5` measured only action-step extrapolation error.

The arbitrary-precision reconstruction itself is exceptionally stable:

| quantity | even | odd |
|---|---:|---:|
| action precision | 80 decimals | 80 decimals |
| largest action imaginary part | `3.37e-77` | `2.91e-76` |
| `epsilon_5` | `1.185e-23` | `1.561e-23` |
| parity Schur difference | \- | `7.34e-18` normalized |

The regular `30 x 30` block has rank 30 at all three thresholds and condition
`35.78` in both parities.  All 180 displaced geometries remain Lorentzian and
off branch boundaries.

## 2. Reconstructed five-dimensional spectrum

The collective values are

```text
even: 2.1458130482e-17
odd : 2.0882154548e-17.
```

The four relative phase-lapse eigenvalues are

```text
even:
  4.604967055079e-8
  4.604967055170e-8
  4.604967055637e-8
  4.604967055742e-8

odd:
  4.604967055134e-8
  4.604967055353e-8
  4.604967055556e-8
  4.604967056203e-8.
```

The reconstructed matrix is, to the displayed scale, the
permutation-symmetric projector form

```text
S_5 approximately a * [I-(1/5) 11^T],
a = 4.604967055...e-8.
```

Its diagonal entries are approximately `3.683973644...e-8` and off-diagonal
entries `-9.20993410...e-9`.  Thus the `4+1` split is not a fragile sorting of
five unrelated small numbers.

## 3. Independent lapse-family control

At each frozen point

```text
rho = rho0 exp(eta),
q   = l0^2-rho,
eta in {-1e-3,0,+1e-3},
```

all 35 analytic per-edge residuals pass `1e-7`.  Maximum values were

```text
even: 7.05e-11, 5.24e-12, 1.00e-10
odd : 3.68e-9,  5.26e-10, 3.49e-9.
```

The normalized analytic tangent

```text
(-rho/q on the 30 diagonals, 1 on the five poles)
```

overlaps the Schur-lifted collective direction by `1.0` to printed precision.
It overlaps the weakest original eigenvector by `0.99999999967` (even) and
`0.99999319` (odd; the original double matrix did not resolve the direction
inside the tiny cluster).

## 4. Attack on the correction protocol

**DERIVED POST-RESULT DIAGNOSTIC:** `epsilon_5` is not the full error budget.
The Schur lifts were formed from the earlier binary64 matrix.  The normalized
lifted collective vector differs from the analytic lapse tangent by

```text
even: norm 5.6191e-10
odd : norm 5.4533e-10.
```

Evaluating the recorded quadratic form on precisely that vector difference
gives

```text
even: 2.1458201e-17
odd : 2.0882261e-17,
```

which accounts for the entire computed collective eigenvalue.  The protocol
compared that value only with the much smaller action-extrapolation error
`epsilon_5`, so `FIVE_STIFF` is mechanically correct under its definition but
does not distinguish physical collective stiffness from lift uncertainty.

This is not grounds to subtract the value post hoc.  The exact analytic lapse
path must be preregistered as the evaluated direction and checked directly.

## 5. Status ledger

| Claim | Status |
|---|---|
| Thirty staircase-diagonal directions are regular | **DERIVED COMPUTATIONAL** |
| Four relative phase-lapse curvatures are approximately `4.604967055e-8` | **DERIVED COMPUTATIONAL** |
| The five-dimensional matrix has a permutation-projector `4+1` pattern | **DERIVED COMPUTATIONAL** |
| Frozen Schur protocol outcome is `FIVE_STIFF` | **DERIVED BY PREREGISTERED RULE** |
| That label establishes physical collective stiffness | **REFUTED AS AN INFERENCE: error model omitted lift uncertainty** |
| The collective direction is exactly null | **OPEN; strongly supported, not yet directly action-tested on the exact path** |
| Lapse/pseudo-constraint mechanism is new | **REFUTED BY PRIOR ART** |
| Exact `4+1` realization on this carrier has external novelty | **OPEN** |

## 6. Next valid correction

Freeze the exact one-parameter path `rho=rho0 exp(t)`,
`q=l0^2-rho` itself and evaluate the complete 80-decimal action on that path.
Do not obtain the collective direction from a double-precision inverse.  The
four relative curvatures need no further search; retain their present frozen
values and error envelope.

Only if the exact-path action is constant and the full residual remains zero
away from the base may the final scoped label become
`ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF`.
