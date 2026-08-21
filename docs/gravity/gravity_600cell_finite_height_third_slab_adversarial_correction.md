# Frozen correction: dual outer brackets and negative-height direct checks

Date: 2026-08-21.

Adversarial protocol commit: `2f1c6f0`.

Adversarial verifier registration commit: `36643b1`.

Preserved first adversarial artifact commit: `8bb6b39`.

First adversarial artifact SHA-256:

```text
b8150a5223b4a3eb01102b8220212ee9da1ce37a6ae78403e6365f1558b9d6ca
```

## 1. Observed failure

The first equal-`mu` run returned `5/9` and

```text
THIRD_SLAB_EXTENDIBILITY_ADVERSARIAL_OPEN.
```

It reproduced the primary negative-side roots, but on branch B printed the
same positive root `q=31.279...` twice, missed the primary exterior root
`q=99.627...`, and consequently reported physical counts `A:0, B:0`.

It also rejected the negative-height algebraic roots with large direct-action
momentum residuals.

## 2. Outer-bracket implementation error

For an interval `(q_stationary,+infinity)`, the helper searching for a point
with the asymptotic sign began at `q=10` without requiring

```text
q > q_stationary.
```

On branch B, `q_stationary=47.637...`; the helper accepted `q=10`, reversed
the intended bracket and bisected back to the already-counted interior root
`q=31.279...`.  This is directly visible in the duplicated artifact rows.

### Sole allowed bracket correction

- An infinite-tail point must start beyond the adjacent finite stationary
  point, with magnitude at least `2*abs(q_stationary)`.
- A one-sided-zero point must start strictly between zero and the adjacent
  stationary point, with magnitude at most `abs(q_stationary)/2`.
- The stationary points, one-sided limit signs, root tolerances and monotone
  interval count remain unchanged.

## 3. Negative-height domain error in the protocol

The constraint-first elimination is algebraic in signed `h` and correctly
contains negative-height reverse roots.  The complete action, however, is
parameterised by

```text
rho=h^2,
sqrt(rho)>0.
```

Substituting a negative reconstructed `h` into `rho=h^2` and then dividing
the direct lapse equation by that negative `h` does not represent the same
positive-oriented action branch.  It is therefore invalid to demand direct
full-action residuals for negative-height algebraic controls.

The primary protocol already handled this correctly: it required direct
action and junction checks only for physical roots.

### Sole allowed domain correction

- Require both reduced equations `C=P=0` for every algebraic root, including
  negative height or negative endpoint scale.
- Require the direct full-action equations and junction only for roots with
  `h>0` and `1+h*q>0`.
- Do not change any algebraic root, physical inequality or count.

## 4. Rerun boundary

Add this correction commit to the artifact provenance, change only the two
mechanisms above, and rerun only the targeted adversarial third-slab
verifier.  Any remaining disagreement is substantive and stops the result.
