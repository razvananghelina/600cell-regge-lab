# Preregistration: generic-velocity leading relational map

Date: 2026-08-21

Prior-art gate commit: `1fcab34`.

Status: frozen before taking any generic-velocity action limit or solving any
mass--velocity constraint.

## 1. Frozen definitions

Use the complete exact cellular action stated in the prior-art note.  Set the
incoming scale to one and retain symbolic

```text
s>0, e>0, v real, mu>0,
L_minus=1,
L_plus=exp(s*v*e),
rho=s^2*e^2.
```

Here `s=1` is a coarse nominal interval and `s=1/2` is a half interval with
the endpoint displacement scaled so that the physical leading velocity `v`
is unchanged.  No numerical `mu` or `v` may enter the derivation.

## 2. Leading principal function

Construct directly

```text
L0(s,v,mu)=lim_(e->0+) S(1,exp(s*v*e),s^2*e^2;mu)/(s*e).
```

The limit must be real on the stated branch.  Prove or refute exact
independence from `s` before evaluating `s=1` or `1/2`.  Record every
trigonometric and inverse-hyperbolic branch used.

## 3. Fixed-endpoint lapse constraint

Let `tau=s*e` and regard the leading endpoint displacement as
`Delta=tau*v`.  A lapse variation holds `Delta`, not `v`, fixed.  Therefore
derive the leading constraint in two independent ways:

1. from the full action,

   ```text
   C_direct(s,v,mu)=lim_(e->0+) 2 F/(s*e),
   F=rho*partial S/partial rho;
   ```

2. from the leading principal function,

   ```text
   C_HJ(v,mu)=L0-v*partial L0/partial v.
   ```

Require exact equality and exact independence from `s`.  Differentiating
`L0` at fixed `v` without the `-v L0'` term is a control failure.

## 4. Incoming momentum

From the complete action define

```text
p_pre=-Pminus,
Pminus=(L_minus/2)*partial S/partial L_minus.
```

Compute

```text
p0(s,v,mu)=lim_(e->0+) p_pre.
```

Require exact independence from `s`.  Independently derive the same result
from the endpoint derivative of the leading principal function, including
the dependence of the dimensionless mass ratio on the incoming scale when
the physical mass is held fixed.  If this chain-rule convention cannot be
closed exactly, report `OPEN` rather than silently treating `mu` as fixed
under a scale derivative.

## 5. Mass--velocity branch census

The constraint is expected to be affine in `mu`, but this is not assumed.
Solve `C(v,mu)=0` exactly for all real positive-mass branches and record

```text
N_mu(v),
mu_j(v),
p_j(v).
```

No branch may be dropped by choosing a numerical velocity.  For every branch
that is real near `v=0`, derive exactly:

```text
mu_j(0),
mu_j'(0),
mu_j''(0),
p_j(0),
p_j'(0).
```

Only after these expressions are frozen, compare `mu_j(0)` with the known
static normalization `(90/pi)epsilon`.  Determine the sign of the first
nonzero mass correction without a floating fit.

## 6. Interval-factor composition gate

A branch is `LEADING_SAME_STATE_COMPATIBLE` only if all of the following hold
exactly:

- `L0`, `C` and `p_pre` are independent of `s`;
- the same `(v,mu_j(v),p_j(v))` satisfies the coarse `s=1` and fine `s=1/2`
  leading equations;
- the fine endpoint displacement is exactly half the coarse displacement;
- two identical leading fine displacements add to the coarse displacement.

Report the hit fraction over all positive-mass branches.  This is a leading
kinematic gate only; no next-order endpoint or momentum comparison is made.

## 7. Independent controls

- The `v->0` limit must recover the exact static action, mass and zero
  physical leading momentum.
- Time reversal must send `v->-v`, leave every mass branch invariant and
  reverse the physical momentum.
- At the disclosed post-symbolic control values

  ```text
  v in {1/5,1/2,1}, s in {1,1/2}, e in {1/200,1/400},
  ```

  direct 100-decimal evaluations of `S/(s e)`, `2F/(s e)` and `p_pre` must
  converge to their exact limits.  These values are controls, not branch
  selectors.

## 8. Outcomes

### `GENERIC_VELOCITY_LEADING_REPARAMETRIZATION`

Report **DERIVED EXACT / STRUCTURAL** only if the exact `s` independence,
constraint, momentum, complete positive-mass branch census, time reversal and
all direct controls pass, and at least one branch is uniquely same-state
compatible near `v=0`.

### `GENERIC_VELOCITY_LEADING_INTERVAL_DEPENDENCE`

Report **DERIVED NEGATIVE, scoped** if a certified nonzero `s` dependence or
coarse/fine leading mismatch survives with all branch conventions resolved.

### `GENERIC_VELOCITY_LEADING_NONUNIQUE`

Report **STRUCTURAL / OPEN selection** if multiple positive-mass branches
pass and the action supplies no selection among them.

### `GENERIC_VELOCITY_LEADING_OPEN`

Use **OPEN** for unresolved symbolic limits, branch domains, chain-rule
conventions or numerical controls.

## 9. Interpretation boundary

A positive result is expected continuum reparametrization kinematics and
does not refute the turning-point `O(e^2)` half-step no-go.  It licenses only
the next-order generic-velocity composition test.  It does not derive a
fundamental or absolute tick, Einstein convergence, a limiting speed, Planck
units, anisotropic propagation or external novelty.

Only the new targeted verifier will be run.  The full suite will not be run.

