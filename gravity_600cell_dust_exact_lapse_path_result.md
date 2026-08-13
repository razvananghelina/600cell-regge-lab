# Exact collective lapse path: one null plus four pseudo-constraints

Date: 2026-08-13

Prior-art update: `0882934`

Upstream five-direction result: `dc927a5`

Frozen exact-path protocol: `515a509`

Implementation commit: `80e24fa`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_exact_lapse_path.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_exact_lapse_path.json`

Targeted run: **11/11 implementation checks passed**.  The full suite was not
run.

## Headline

> **DERIVED COMPUTATIONAL:** on both ordered-schedule parities, the complete
> one-slab dust action has a stationary one-parameter collective lapse family
> `rho=tau^2 exp(t)`, `q=l0^2-rho`.  The collective tangent is action-null,
> while four relative phase-lapse directions retain the resolved curvature
> `4.604967055...e-8`.

The scoped structure is therefore

```text
30 regular staircase-diagonal directions
+4 small but nonzero relative lapse pseudo-constraints
+1 collective lapse null direction.
```

## 1. Direct 100-decimal path result

The action was evaluated, without a numerical Hessian lift, at

```text
t = -0.1, -0.03, -0.01, -0.003, -0.001,
     0,
     0.001, 0.003, 0.01, 0.03, 0.1.
```

For every point,

```text
rho(t) = tau^2 exp(t),
q(t)   = l0^2-rho(t).
```

The maximum normalized change in the complete 100-decimal action was

```text
even: 8.250e-97
odd : 9.668e-94,
```

against the preregistered `1e-50` threshold.  Maximum imaginary action parts
were `2.30e-97` and `4.61e-97`.

The extrapolated exact-path second derivatives were

```text
even: 3.898e-91 + 3.24e-94 i
odd : 1.923e-88 - 3.84e-92 i,
```

against the frozen `1e-40` null threshold.  The corresponding first
derivatives are also compatible with zero (`3.70e-95` and `1.90e-91` in
absolute real part).

## 2. Full stationarity, not merely tangential constancy

At all eleven points the verifier also evaluated all 35 analytic internal
equations, rather than only the path derivative.  Maximum absolute per-edge
residuals over the whole path were

```text
even: 2.381e-10
odd : 7.791e-9,
```

well below the frozen `1e-7` threshold.  All representative simplices remain
Lorentzian, with minimum Gram moduli `9.414e-5` and `4.707e-5`, and minimum
angle-argument modulus `0.9953`.

Thus this is computational evidence for a stationary family, not merely an
accidental zero directional derivative at one point.

## 3. Four relative modes retained without reselection

The verifier did not rediagonalize or choose four favorable directions.  It
loaded the already committed 80-decimal relative eigenvalues:

```text
even: 4.604967055079e-8 ... 4.604967055742e-8
odd : 4.604967055134e-8 ... 4.604967056203e-8.
```

Every value exceeds `100*epsilon_5` by many orders of magnitude.  Their
mechanism is **KNOWN** from Regge pseudo-constraint literature; this
calculation supplies the explicit five-phase 600-cell realization.

## 4. Reconciliation with the previous `FIVE_STIFF` label

The previous preregistered rule called the approximate Schur matrix
`FIVE_STIFF`.  That historical outcome remains recorded.  Its collective
`2.1e-17` eigenvalue was produced entirely by a `5.5e-10` error in the
double-derived Schur lift, an error not included in its `epsilon_5`.

The present protocol did not subtract that eigenvalue.  It replaced the
approximate lift by the exact published path specified before evaluation.
The much stronger direct-action result resolves the earlier framing failure
in favor of one collective null direction.

## 5. Physical meaning

- **DERIVED COMPUTATIONAL:** `tau` is not selected by the one-slab equations
  on this time-symmetric dust family.  It is collective lapse/gauge data.
- **DERIVED COMPUTATIONAL:** after fixing that gauge, the invariant internal
  Hessian has 34 nonzero directions: thirty regular and four very soft
  relative pseudo-constraint directions.
- **STRUCTURAL:** the five-stage schedule breaks independent phase-lapse
  freedom down to one collective freedom; the four relative combinations
  acquire small stiffness.
- **NOT DERIVED:** a physical duration, speed of light, Planck time or
  preferred tick.
- **NOT DERIVED:** a continuum Hamiltonian constraint algebra or a graviton.
- **NOT TESTED:** the unreduced 840-edge carrier and multi-slab propagation.
- **NOT NEW IN GENERAL:** lapse nullity and curvature-induced
  pseudo-constraints are established Regge phenomena.
- **OPEN:** external novelty of this explicit `30+4+1` 600-cell computation.

This costs an attractive but false interpretation: the input
`tau=0.0102` reproduced from the paper is a coordinate/lapse choice, not a
fundamental tick produced by the geometry.

## 6. Correct next step

The full `35 x 35` implicit-function theorem was the wrong formulation because
one exact gauge direction must be quotiented out.  The next target-free audit
is the gauge-fixed quotient:

1. use the exact normalized lapse tangent `w`;
2. construct a deterministic orthonormal basis of `w^perp`;
3. project the internal Hessian to `34 x 34` and test its regularity;
4. verify the boundary-response block is orthogonal to the redundant gauge
   equation;
5. compute the unique 34-coordinate internal response to all thirty allowed
   final-boundary directions, modulo collective lapse.

If those gates pass, the one-slab action supports local nonhomogeneous
boundary evolution in this invariant sector, with one gauge parameter.  That
would still be a linear local result; displaced nonlinear roots and multiple
ticks would remain separate preregistered missions.

## 7. Post-result coordinate correction

An initial audit after the frozen result incorrectly projected the published
raw per-edge residual `(partial S/partial x)/24` onto a Hessian of the
logarithmic equations `x*(partial S/partial x)/24`.  It omitted the factors
`x`, most dramatically `rho approximately 1e-4` on the pole equations.  The
resulting claimed odd correction `8.675e-3` is **RETRACTED**.

Using coordinate-consistent logarithmic residuals gives at `t=0`:

```text
                         even            odd
norm(log residual)       2.576e-13       2.729e-12
predicted correction     2.521e-10       9.019e-7
```

Both correction proxies are below the subsequently frozen `1e-5` weak-scale
tolerance.  Therefore the printed base remains adequate for the linear
matrix audit.  What is still not established by the old `1e-7` raw-gradient
gate is a high-accuracy stationary-family theorem across the complete
collective interval.  That narrower question requires smaller-step
arbitrary-precision logarithmic derivatives.
