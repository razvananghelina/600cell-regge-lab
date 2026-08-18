# Post-result preregistration: falsifiability controls for conformal/shape closure

Date: 2026-08-18

Primary prior-art commit: `42313fb`  
Primary protocol commit: `4039244`  
Registered primary implementation: `3edfa46`  
Preserved first artifact commit: `9bab7a3`

## Why this addendum is necessary

The first frozen execution returned

```text
CONFORMAL_SHAPE_DYNAMICS_DECOUPLED
224/224 primary residuals ZERO_CONSISTENT
112/112 schedule comparisons SCHEDULE_ROBUST.
```

This result is known while writing the present addendum.  Therefore the
following controls are **post-result robustness tests**, not part of the
original blind evidence.  The first artifact and verifier are preserved by
the commits above.

The primary zero is scientifically useful only if the same matrices and error
classifier can resolve mixing for declared alternative splits of the same
dimensions.  Otherwise `ZERO_CONSISTENT` could merely mean that the error band
is too broad or that every symmetry-compatible split is invariant.

## Frozen negative controls

For every schedule, sector, derivative variant and each `X` in
`{Gamma,Omega}`, use the same dimension

```text
r=5d,   n-r=25d.
```

No control is a physical candidate and none changes the primary conformal
carrier.

### 1. Euclidean spectral split of `H`

Let `U_+` and `U_-` be the positive and negative Euclidean eigenspaces of

```text
H=(M+M*)/2.
```

Their dimensions are already frozen as `r` and `n-r`.  With `g` the gap
between the smallest positive and largest negative eigenvalues, define

```text
eta_spec = 2 epsilon_H/(g-2 epsilon_H)
         + 1000 eps_machine n.
```

The control residuals are

```text
(I-U_+U_+*) X U_+,
(I-U_-U_-*) X U_-.
```

This split was computed and shown distinct from the conformal image before
the primary dynamic test.  It is therefore not chosen because of the new
residuals.

### 2. Fixed dense Fourier coordinate split

Let

```text
F[j,k] = exp(2*pi*i*j*k/n)/sqrt(n),
```

in the already frozen minimal-sector coordinate order.  Use its first `r`
columns and remaining `n-r` columns as a deterministic dense split.  Record

```text
eta_F = ||F*F-I||_2 + 1000 eps_machine n
```

and test invariance of both factors.  This split is deliberately a numerical
power control, not a geometric proposal.

For either control basis `Z`, use

```text
R(Z,X)=(I-ZZ*) X Z,
epsilon_control = epsilon_X
                + 2 eta_Z (||X||_2+epsilon_X)
                + 1000 eps_machine n max(1,||X||_2).
```

Apply the unchanged `10/100` zero/open/nonzero classification.

## Complete power ledger and frozen verdict

There are

```text
2 schedules * 7 sectors * 4 variants * 2 operators = 112
```

power cells and four residuals per cell, hence `448` negative-control
classifications.  A cell is:

- `POWER_HIT` if at least one of its four controls is `NONZERO_RESOLVED`;
- `POWER_OPEN` if none is resolved nonzero and at least one is open;
- `POWER_ZERO` if all four are zero-consistent.

The original positive result is upgraded to

```text
CONFORMAL_SHAPE_DYNAMICS_DECOUPLED_POWER_CERTIFIED
```

only if all `112` cells are `POWER_HIT` while all `224` primary residuals
remain zero-consistent.  If at least one cell lacks a hit, the final positive
verdict is instead

```text
CONFORMAL_SHAPE_DECOUPLING_POWER_OPEN.
```

A primary resolved nonzero continues to give
`CONFORMAL_SHAPE_MIXING_REFUTED`, regardless of the controls.

## Interpretation firewall

Passing these controls would show that the primary zero is selective at the
declared resolution.  It would not prove that the conformal/shape split is
the only invariant split, that the shape factor is transverse-traceless, or
that it contains physical gravitons.  Failure would weaken the evidence; it
would not manufacture a negative physical result from an arbitrary control.
