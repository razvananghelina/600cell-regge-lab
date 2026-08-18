# Second-tick canonical-target homotopy: fixed-grid result

Date: 2026-08-16

## Provenance

- prior-art gate: `e760462`;
- frozen protocol: `a2564d5`;
- implementation and registration before evaluation: `f520498`;
- artifact SHA-256:
  `882f2ebbe9146a90b573eca6f9b3c0b7080004825264761c4f2a077cfe0906b3`.

Only the targeted verifier was run.  It returned **7/7**.  The full suite was
not run.

## Mechanical verdict

```text
SECOND_TICK_HOMOTOPY_NEWTON_OPEN
```

The fixed 32-step, zero-order-predictor homotopy does not reach `lambda=1`.
Both schedule parities pass the complete gates through `lambda=4/32=0.125`
and fail the six-iteration Newton boundary at `lambda=5/32=0.15625`.

This is a second algorithmic negative, not proof that the second-tick equations
have no root.

## Certified path segment

At the five accepted nodes `lambda=0,...,4/32`, the absolute upper-scale log,
absolute lapse-square log and weak singular value are

```text
lambda      b                    r                 s_min
0           -3.1160596e-6       -3.5592531e-6     4.2445618e-9
1/32        -3.0186805e-6       -0.1290814870      3.4974065e-9
2/32        -2.9213014e-6       -0.2670683715      2.8435162e-9
3/32        -2.8239221e-6       -0.4152857978      2.2766732e-9
4/32        -2.7265426e-6       -0.5753731572      1.7906599e-9
```

Every accepted node satisfies all 35 internal equations, all 30 interpolated
momentum equations, the Lorentzian branch and the calibrated rank-two gate.
The two parities agree to the stored precision.

At `lambda=5/32`, after the frozen six corrections, both parities stop at

```text
(b,r)=(-2.6710256603e-6,-0.7206140179),
F=(-4.1021338e-11,-1.8614943e-8).
```

The node is not accepted and has no endpoint rank certificate.

## Interpretation ledger

- **DERIVED COMPUTATIONAL NEGATIVE:** the preregistered fixed-grid homotopy
  fails at node 5 under its exact corrector limit.
- **DERIVED COMPUTATIONAL:** a connected rank-two solution segment exists at
  least to `lambda=0.125` and is schedule independent.
- **PATTERN:** the lapse falls rapidly and the weak singular value decreases
  along that segment.  This is consistent with approach to a fold or
  pseudo-constraint boundary, but does not prove either.
- **OPEN:** whether a root exists at the desired `lambda=1`, whether the
  connected branch folds before it, and whether disconnected branches exist.
- **NOT DERIVED:** a second physical tick, continued contraction, a collapse
  time or emergent time.

## Why solver tuning now stops

The direct Newton attempt and the fixed homotopy have both failed by
algorithmic limits.  Increasing iterations or refining lambda after seeing
the failures would create a solver look-elsewhere sequence without resolving
the physical question.

The next calculation should instead remove the target from the search:
parameterize the complete internally stationary homothetic curve at the fixed
second lower geometry and compute its canonical pre-momentum envelope.  Its
fold/maximal momentum is defined without the desired target.  Only after that
envelope is committed should it be compared with the inherited momentum.

If the target lies above a certified global envelope on the connected
Lorentzian branch, the second canonical homothetic tick is impossible on that
branch.  If it lies inside, a bracketed scalar solve is licensed without
further Newton fitting.
