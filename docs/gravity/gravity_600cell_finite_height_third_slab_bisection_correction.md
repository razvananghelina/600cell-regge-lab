# Frozen correction: separate topology sign margins from root convergence

Date: 2026-08-21.

Protocol commit: `0f31fe8`.

Verifier registration commit: `beac60d`.

Preserved first artifact commit: `564dca9`.

First artifact SHA-256:

```text
8e7f716984ee8779a51a3678f8a69b1dad0405e0a1520d2a20cbd4b44e696806
```

## Observed failure

The first run returned `6/8` and

```text
THIRD_SLAB_EXTENDIBILITY_OPEN.
```

It completed the monotone all-real census and provisionally found physical
counts

```text
A: 0,
B: 1,
```

but the direct full-action residual for the physical B root was exactly at
approximately `2e-90`, just above the preregistered strict `1e-90` bound.

## Diagnosed implementation error

The bisection routine used

```text
SIGN_TOL=1e-90
```

both to certify that interval and tail signs were safely separated from zero
and to choose the next bisection half.  Values with magnitude below
`SIGN_TOL` were therefore treated as sign zero inside the root iteration,
even though the requested root convergence was

```text
ROOT_TOL=1e-115.
```

The loop consequently converged to the boundary of the `1e-90` sign band,
not to the zero of the elimination function.  Every printed reduced
constraint residual was correspondingly `+-2e-90`.  This is an algorithmic
precision mismatch, not evidence against the recurrence or the provisional
root count.

## Sole allowed correction

Keep the `SIGN_TOL` logic unchanged for:

- critical-level separation;
- stationary-point nonzero certificates;
- infinite-tail signs;
- deciding whether a monotone interval contains a root.

Inside `bisect_root` only:

1. require the already-certified endpoint product to be strictly negative;
2. test exact zero directly;
3. stop only at `abs(f)<ROOT_TOL` or interval width below `ROOT_TOL`;
4. choose the next half from the raw product `f_left*f_middle`, not from the
   topology-margin sign helper.

No seed, root interval, state recurrence, branch label, physical inequality,
tail formula, tolerance or outcome rule may change.  Rerun only the targeted
third-slab verifier.
