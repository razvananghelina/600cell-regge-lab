# Preregistration: connected canonical continuation of the 600-cell dust slab

Date: 2026-08-16

Prior-art gate: `52a6d50`

Upstream canonical-rank result: `715b6ad`

Status: **frozen before evaluating any nonlinear state away from the
published canonical datum**.

## 1. Frozen equation, targets and carrier

For both derived schedule parities retain the complete action, branch and
coordinate order of the canonical-rank verifier.  Hold `q_old` at the
published regular boundary and solve

```text
F(y;lambda) = (g_internal[35], p_pre(y)-p(lambda)[30]) = 0,
y = (log x[35], log q_new[30]),
p_pre(y) = -g_old(y),
p(lambda) = p_pre,pub + lambda*(P p_post,pub-p_pre,pub).
```

Load the two endpoint momenta and the derived orbit permutation `P` from the
corrected two-slab artifact.  No momentum coefficient is optimized.
`lambda` is a numerical homotopy only; it is not physical time.

The only physical candidate endpoint in this protocol is `lambda=1`.
`lambda=0` is a reproduction control.

## 2. Symmetry-fixed predictor

For invariant data, the locally unique solution furnished by the certified
rank-65 Jacobian lies in the regular fixed subspace.  Parameterize it by

```text
r = log(common staircase-diagonal magnitude),
t = log(common positive pole magnitude),
b = log(common final-boundary square).
```

The reduced residual is the mean of each of the three full equation types:

```text
R = (mean g_staircase[30], mean g_pole[5],
     mean(p_pre-p(lambda))[30]).
```

Every accepted reduced state must be expanded back to 65 coordinates and
pass all full residual and within-type spread gates below.  The reduction is
not used to claim absence of nonsymmetric roots after branch loss.

The committed base inverse gives the design diagnostic

```text
Delta(r,t,b) at lambda=1 approximately
(+5.4091038624e-5, -4, 0).
```

This value is printed only as an upstream linear comparison.  It is not a
nonlinear initial guess for `lambda=1`.

## 3. Arbitrary-precision evaluator and Jacobian

Use 100 decimal digits and the analytic Schlaefli gradient independently
implemented in the canonical-rank verifier.  At every evaluated state require:

- all representative 4-simplices have exactly one timelike Gram direction;
- every leading principal minor is nonzero;
- minimum complex angle-argument modulus exceeds `1e-6`;
- maximum imaginary contamination in action and gradient is below `1e-70`;
- every positive-magnitude coordinate remains positive.

Differentiate the three reduced residuals in `(r,t,b)` by calibrated central
differences.  Freeze the same pairs as the rank census:

```text
operational: (1e-20, 1e-15),
validation : (3e-20, 3e-15).
```

The smaller step is primary and primary minus shadow is its error proxy.
Require every operational/validation cross difference to be at most ten
times the sum of the two proxy magnitudes plus `1e-70`.  Define the combined
spectral Jacobian error exactly as in the rank census.  A Newton inversion is
allowed only when the smallest operational singular value exceeds 100 times
that error.  Report all three singular values at every accepted point.

## 4. Reproduction control

Use `lambda=0`.  Define the exact reduced collective-lapse tangent

```text
w3 = (-rho/(l0^2-rho), 1, 0)
```

and normalize it in the inherited multiplicity metric

```text
30 dr^2 + 5 dt^2 + 30 db^2.
```

Start at the published logarithms plus `1e-3*w3_normalized`.  This seed and
amplitude are frozen before evaluation.  Apply the Newton corrector of
section 5.  The control passes only if:

1. the solver reaches an operational residual infinity norm below `1e-50`;
2. a disjoint validation Jacobian passes its calibration;
3. the full 65-residual infinity norm is below `1e-40` and the maximum
   within-type spread is below `1e-50`;
4. the recovered logarithms differ from the published point by less than
   `1e-10` in infinity norm;
5. all branch gates pass.

Stop before forward continuation if either parity fails.

## 5. Frozen Newton corrector

At a fixed `lambda`, use the operational primary Jacobian and solve

```text
J delta = -R.
```

Test damping values `alpha=1,1/2,...,2^-20` in that order.  Accept the first
branch-valid trial satisfying the fixed Armijo condition

```text
||R_trial||_infinity <= (1-alpha/4)*||R||_infinity.
```

Recompute the Jacobian after every accepted step.  Stop successfully at
`||R||_infinity < 1e-50`; otherwise stop after 30 accepted Newton steps, no
accepted damping, a Jacobian error-band failure or a branch failure.  No
Broyden update, SciPy optimizer, alternate seed, random restart or tolerance
change is allowed.

## 6. Connected continuation rule

After reproduction, attempt the fixed coarse targets

```text
lambda_k = k/64,  k=1,...,64,
```

in order.  At each accepted point compute the tangent predictor from

```text
J dy/dlambda = (0,0,mean(P p_post-p_pre)).
```

Use the predicted point as the sole Newton seed for the next target.

If a coarse target fails, bisect only the interval between the last accepted
`lambda` and that failed target.  Always continue from the last accepted root
with its tangent predictor.  Perform at most 20 such bisections.  Do not skip
past a failed interval and do not search another branch.  This defines the
connected branch operationally without choosing roots after inspection.

At every accepted root require the section-3 calibration and branch gates,
full 65-residual infinity norm below `1e-40`, and within-type spread below
`1e-50`.

## 7. Endpoint diagnostics

At every root record:

```text
rho/rho0,
(slant_square+rho-l0^2)/l0^2,
Delta log(q_new/l0^2),
Delta spatial length scale = 0.5 Delta log(q_new/l0^2),
all reduced Jacobian singular values,
all 65 residuals and branch margins.
```

At `lambda=1`, classify the spatial endpoint as expanding, contracting or
zero-consistent using the frozen band:

```text
expanding:   Delta length log > +1e-12,
contracting: Delta length log < -1e-12,
static:      absolute value <= 1e-12.
```

This classification is kinematic.  It does not identify physical elapsed
time or a Friedmann continuum limit.

## 8. Mechanical outcome hierarchy

Assign outcomes in this order for each parity:

1. `CANONICAL_CONTINUATION_CONTROL_FAILED` if the evaluator or reproduction
   control fails;
2. one of the three `CANONICAL_FORWARD_ROOT_*` outcomes if the connected
   branch reaches `lambda=1` and passes every endpoint gate;
3. `CONNECTED_BRANCH_APPROACHES_ZERO_LAPSE` only if, after all 20 bisections,
   the last accepted point has

```text
0.49 <= lambda < 0.5,
rho/rho0 < 1e-12,
abs(Delta log(q_new/l0^2)) < 1e-8,
abs((slant_square+rho-l0^2)/l0^2) < 1e-10;
```

4. `CANONICAL_CONTINUATION_BRANCH_TERMINATED` if the unresolved boundary of
   the connected path is a Lorentzian/positivity branch failure and item 3
   does not apply;
5. `CANONICAL_CONTINUATION_NUMERICALLY_OPEN` for any remaining unresolved
   Jacobian, damping or iteration failure.

Item 4 is a computational statement about the frozen connected-continuation
procedure.  It is not a theorem that no other root exists.

## 9. Acceptance and claim boundary

The first-frame acceptance boundary is that **both** parity carriers reach
`lambda=1`, pass all controls and agree on the kinematic scale change within
the larger of the combined final Newton correction proxy and the already
frozen absolute kinematic band `1e-12` in logarithmic length scale.  A reached
root may be expanding, contracting or static; none is relabelled to obtain a
preferred answer.

Any other clean outcome blocks the claim of a canonically selected next
frame.  An unrelated root is not a rescue.

Even acceptance establishes only a symmetric reduced Regge/Friedmann step.
Perturbative shape propagation, refinement, a physical clock, inertia,
inflation, `c`, Planck units and particle masses remain OPEN or NOT TESTED.
