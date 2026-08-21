# Preregistration: generic-velocity next-order lapse and composition gate

Date: 2026-08-21

Prior-art gate commit: `62d92ab`.

Status: frozen before evaluating any next-order coefficient or common root.

## 1. Frozen exact input and state

Use the complete cellular action reconstructed from the certified homothetic
frustum artifact.  Freeze its SHA-256 and the hashes of both accepted leading
generic-velocity artifacts in the verifier before evaluating a series.

For symbolic real `v!=0`, define exactly

```text
r(v)       =sqrt(v^2+4),
theta(v)   =acos((v^2+2)/(2(v^2+3))),
eta(v)     =asinh(v/sqrt(8(v^2+3))),
epsilon_v  =2*pi-5*theta(v),
M          =mu(v)=180 epsilon_v/[pi r(v)],
p0         =p(v)=180 v epsilon_v/r(v)-600 sqrt(3) eta(v).
```

Set `L_minus=1`.  The same literal `M` and `p0` are used for a coarse slab,
the first fine slab and every comparison at the same `v`.  Do not introduce
`M(h)`, `p0(h)`, a fitted endpoint coefficient or a fitted velocity.

Let

```text
h=sqrt(rho)>0
```

be the proper slab height.  It is the candidate relational interval, not a
coordinate label.

## 2. One-slab next-order jet

Use the complete analytic endpoint ansatz

```text
L_plus=exp(v h+a h^2).
```

The omitted `O(h^3)` endpoint jet cannot affect the registered first
correction.  From the unexpanded action define

```text
F       =rho partial_rho S,
p_pre   =-L_minus partial_(L_minus)S/2,

C1(v,a)=lim_(h->0+) [2F/h]/h,
P1(v,a)=lim_(h->0+) [p_pre-p0]/h.
```

First verify that the zeroth-order residuals vanish exactly by the accepted
leading theorem.  Then derive `C1` and `P1` symbolically without substituting
a numerical velocity.

## 3. Complete common-root census

Determine all real common roots in `a` of

```text
C1(v,a)=0,
P1(v,a)=0
```

for real `v!=0`.  Record:

```text
degrees in a,
leading coefficients,
gcd or resultant,
all generic real roots,
all exceptional real nonzero velocities,
all excluded denominator or inverse-function branch loci.
```

No sampled velocity may establish a global root count.  Squared equations
must be substituted back into the unsquared residuals.  If the exceptional
set cannot be classified exactly, the global verdict is `OPEN` even when all
numerical controls agree.

Interpret the census separately from composition:

- a common `a(v)` gives a local arbitrary-small-duration family through this
  order and therefore does **not** select a tick;
- no common `a` is a local fixed-state discretization obstruction, not an
  isolated positive tick;
- isolated finite positive roots are outside this asymptotic census and
  require a later bounded exact-finite protocol.

## 4. Turning-point and hostile controls

The generic theorem excludes `v=0`, but a separately evaluated static ansatz

```text
L_plus=exp(A h^2),
M=90 epsilon/pi,
p0=180 epsilon h
```

must reproduce the already-certified coarse turning-point factors, up to one
explicit common nonzero normalization per equation:

```text
A(D A+4 epsilon),
D A+4 epsilon,

epsilon=2*pi-5*acos(1/3),
D=5sqrt(2)/3-epsilon.
```

This is a control of the series machinery, not part of the generic branch
census.

Two hostile state changes must fail before the next-order census is accepted:

1. replace `M` by `mu(v)+1/10`; the zeroth-order lapse residual must be
   exactly `-4*pi/5`;
2. replace `p0` by `p(v)+1/10`; the zeroth-order momentum residual must be
   exactly `-1/10` with the registered sign convention.

## 5. Numerical coefficient controls

Only after the exact expressions exist, use 100 decimals at

```text
v in {1/3,4/5,3/2},
a in {0,1/7},
h in {1/400,1/800}.
```

Compare the direct unexpanded quotients

```text
[2F/h]/h,
[p_pre-p0]/h
```

with `C1(v,a)` and `P1(v,a)`.  For

```text
error(h)=abs(direct-exact)/max(1,abs(exact)),
```

require either both errors below `1e-70`, or decreasing resolved errors with

```text
0.8 <= log2[error(1/400)/error(1/800)] <= 1.2.
```

Mixed resolved/unresolved pairs are `OPEN`.  These points validate the
coefficient extraction and never select a branch.

## 6. Conditional two-half-slab composition

Execute this section only if the one-slab census gives a unique common
analytic root `a(v)` on a completely classified real domain.

Define

```text
L_mid    =exp(v h/2+a(v) h^2/4),
L_coarse =exp(v h+a(v) h^2),
L_fine   =exp(v h+b h^2),
rho_c    =h^2,
rho_1=rho_2=h^2/4.
```

For the second fine slab derive the first nonzero coefficients of

```text
F_2=0,
G_mid=P_plus(1,L_mid,rho_1;M)
      +P_minus(L_mid,L_fine,rho_2;M)=0,
```

where `P_boundary=L_boundary partial_(L_boundary)S/2`.  Enumerate every real
common root `b(v)` exactly.  A valid two-slab history must also retain the
already-passed first-slab lapse and initial-momentum coefficients.

For every surviving `b(v)`, record separately:

```text
endpoint defect:       b(v)-a(v),
final momentum defect: lim_(h->0+) [P_plus,fine-P_plus,coarse]/h,
action defect:         lim_(h->0+) [S_1+S_2-S_coarse]/h^2.
```

All three must vanish for `NEXT_ORDER_COMPOSITION`.  A nonzero defect is a
variational-integrator/refinement error; it is not evidence for a fundamental
tick.

## 7. Outcome hierarchy

### `GENERIC_DURATION_FREE_TO_NEXT_ORDER`

Use **DERIVED EXACT / STRUCTURAL** if one common one-slab root exists on the
complete registered real domain.  State explicitly that no relational tick
is selected to this order.  Add `_COMPOSITIONAL` only if the conditional
two-slab endpoint, momentum and action defects all vanish.

### `GENERIC_FIXED_STATE_NEXT_ORDER_OBSTRUCTION`

Use **DERIVED NEGATIVE, scoped** if the exact one-slab equations have no
common root for generic `v` and every exceptional nonzero velocity has been
classified.  Do not call this a tick.  The next allowed calculation is a
separately preregistered exact finite positive-root census.

### `GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES`

Use **STRUCTURAL / OPEN selection** if only a finite or lower-dimensional
nonzero-velocity set survives.  Report its exact set before any physical
interpretation.

### `GENERIC_NEXT_ORDER_OPEN`

Use **OPEN** for unresolved limits, transcendental sign/domain questions,
incomplete exceptional sets, failed controls or incompatible symbolic and
numerical routes.

No outcome derives an absolute unit, `c`, `G`, Planck time or external
novelty.  Only the new targeted verifier will be run; the full suite remains
out of scope.
