# Preregistration: cubic generic-velocity formal-integrability gate

Date: 2026-08-21

Prior-art gate commit: `70e7ca2`.

Status: frozen before evaluating `C2`, `P2`, any root in `c`, or any
exceptional cubic velocity.

## 1. Frozen provenance, carrier and physical state

Use the certified homogeneous cellular 600-cell Regge-plus-conserved-dust
action on the positive Lorentzian branch with zero cosmological constant.
Freeze the accepted primary and adversarial next-order artifacts:

```text
gravity_600cell_generic_velocity_next_order.json
SHA-256 4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18

gravity_600cell_generic_velocity_next_order_adversarial.json
SHA-256 3ab16e6d19b527590b3dce6e8b3caa093efb6cc504a2a7824362ffc529a83a05
```

Set `L_minus=1` and keep exactly the same physical incoming state

```text
M=mu(v),
p0=p(v),
```

for real

```text
v!=0,
K(v^2)!=0.
```

The already-classified pair `v=+-v_star` is excluded because no quadratic
endpoint jet exists there.  Freeze

```text
a(v)=-B(v^2)/K(v^2)
```

and vary only the new coefficient `c` in

```text
L_plus=exp(v h+a(v)h^2+c h^3),
rho=h^2,
h>0.
```

No `h`-dependent mass, momentum, velocity or lower-order endpoint
coefficient is permitted.

## 2. Registered coefficients and extraction route

Let

```text
F     =rho partial_rho S,
p_pre =-L_minus partial_(L_minus)S/2.
```

The accepted branch cancels the constant and linear terms of `2F/h` and
`p_pre-p0`.  Define the next coefficients exactly by

```text
C2(v,c)=lim_(h->0+) [2F/h]/h^2,
P2(v,c)=lim_(h->0+) [p_pre-p0]/h^2.
```

The primary route will differentiate the already-derived scaled exact
expressions with variables `(lm,lp,q,w,tau)`.  For

```text
lp=exp(vh+ah^2+ch^3),
q=(lp-1)/h,
w=600 sqrt(3)(1-lp^2)/h,
tau=h,
```

use the exact path data

```text
base:
  lm=1, lp=1, q=v, w=-1200 sqrt(3)v, tau=0

coefficient of h:
  lm=0
  lp=v
  q=a+v^2/2
  w=-1200 sqrt(3)(a+v^2)
  tau=1

coefficient of h^2:
  lm=0
  lp=a+v^2/2
  q=c+va+v^3/6
  w=-1200 sqrt(3)(c+2va+2v^3/3)
  tau=0.
```

For a scaled expression `f`, its registered `h^2` coefficient is

```text
sum_i f_i y2_i + (1/2) sum_(i,j) f_ij y1_i y1_j
```

at the base point.  Derive this coefficient before substituting
`a=-B/K`; then substitute the frozen `a` and normalize exact positive
radicals.  If normalization exposes a new composite radical, preserve the
failed run first and freeze a complete positive-factorization inventory in a
new commit before changing the implementation.

## 3. Controls before the cubic census

Before interpreting `C2` or `P2`, the verifier must reconstruct exactly:

1. the accepted `C1(v,a)` and `P1(v,a)`;
2. their common root `a=-B/K`;
3. the exact cross identity proving equality of the two quadratic roots;
4. vanishing of both lower-order residuals after the frozen substitution.

Two hostile changes must be detected:

1. `M -> mu(v)+1/10` gives the already-certified leading lapse defect
   `-4*pi/5`;
2. `a -> -B/K+1/10` leaves a nonzero first-correction residual throughout
   the registered domain `v!=0`, `K!=0`.

The hostile tests validate state and order sensitivity; they are not part of
the root census.

## 4. Complete exact root census

Treat `C2` and `P2` as exact functions of `c`.  Before numerical sampling,
record:

```text
degrees in c,
leading coefficients,
numerators and denominators,
gcd or resultant,
all generic common real roots,
all real degree-drop loci,
all denominator loci,
all inverse-function branch loci,
parity under v -> -v.
```

Every candidate introduced by clearing denominators or squaring must be
substituted back into the original unsquared coefficients.  The intervals

```text
0<|v|<v_star,
|v|>v_star
```

must be covered separately.  A sign scan, a finite grid or a collection of
sampled velocities cannot establish completeness.  If the transcendental
zero set cannot be certified over the full registered domain, the global
outcome is `GENERIC_CUBIC_OPEN`.

## 5. Arbitrary-precision coefficient controls

Only after the symbolic coefficients and census exist, use 100-decimal
unexpanded evaluations at

```text
v in {1/2,3/2,3},
c in {0,1/11},
h in {1/1000,1/2000,1/4000,1/8000}.
```

Compare

```text
[2F/h]/h^2,
[p_pre-p0]/h^2
```

with `C2(v,c)` and `P2(v,c)`.  For each quantity require either all errors
below `1e-70`, or strictly decreasing resolved errors whose three successive
base-two orders all lie in `[0.8,1.2]`.  Mixed resolved/unresolved sequences
are failures.  These points validate extraction only and may not choose or
remove a branch.

## 6. Interpretation and outcome hierarchy

### `GENERIC_DURATION_FREE_TO_CUBIC_ORDER`

Use **DERIVED EXACT / STRUCTURAL** only if both equations possess the same
real `c(v)` on the complete registered domain, with every exceptional locus
classified.  This means that an arbitrary-small-duration endpoint jet is
formally integrable through cubic endpoint order.  It does not select `h`
and is not a tick.

### `GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION`

Use **DERIVED NEGATIVE, scoped** if the exact equations have no common `c`
generically and every exceptional locus is completely classified.  This is
a local pseudo-constraint/discretization obstruction, not an isolated
positive duration.

### `GENERIC_CUBIC_EXCEPTIONAL_STRATA`

Use **STRUCTURAL / OPEN selection** if common roots survive only on a finite
or lower-dimensional set that is completely classified.  No exceptional
velocity is a limiting speed or tick without a separate propagation and
refinement theorem.

### `GENERIC_CUBIC_OPEN`

Use **OPEN** for incomplete real-root classification, unresolved limits,
failed controls, or disagreement between exact and numerical routes.

This gate tests one slab only.  Composition is a separate, later gate and
will be attempted only after one-slab cubic formal integrability is settled.
No outcome here derives an absolute unit, seconds, `c`, `G`, Planck time or
external novelty.  Only the targeted verifier will be run; the full suite is
out of scope.
