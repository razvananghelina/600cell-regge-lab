# Adversarial protocol: cubic fixed-state obstruction

Date: 2026-08-21

Primary result commit: `08bdde5`.

Primary artifact SHA-256:
`1d35b46cd4db20df0af3ed3e6b5de676d69753cf5059e0eb607d1eec949b9103`.

Status: frozen after the primary result, before implementing or executing the
adversarial derivation.

## 1. Claim under attack

For the same certified action, fixed incoming state and domain

```text
L_minus=1,
M=mu(v),
p0=p(v),
v!=0,
K(v^2)!=0,
a=-B/K,
L_plus=exp(vh+ah^2+ch^3),
rho=h^2,
```

the primary route claims that the cubic lapse and incoming-momentum
coefficients have no common real `c`.  Its decisive expression is

```text
Delta(v)=129600 epsilon(v)^2/(v^2+4)>0.
```

The adversarial route must rederive `Delta`; reading it from the primary JSON
does not count.

## 2. Mechanically different route

Differentiate the complete unscaled action first with respect to `rho` and
`L_minus`.  Only afterward substitute

```text
L_minus=1,
L_plus=1+tau q,
rho=tau^2.
```

This produces exact two-variable functions `Fhat(tau,q)` and
`Phat(tau,q)`.  Use the univariate endpoint quotient

```text
q(tau)=v+q1 tau+q2 tau^2+O(tau^3),
q1=a+v^2/2,
q2=c+va+v^3/6.
```

Derive its second Taylor coefficient directly as

```text
(1/2)Fhat_tt + q1 Fhat_tq + (1/2)q1^2 Fhat_qq + q2 Fhat_q
```

and likewise for `Phat`, at `(tau,q)=(0,v)`.  This two-variable
derivative-first route must not import the primary path Hessian or its
`(lm,lp,q,w,tau)` coefficient formulas.

Read and compare the primary formulas only after both adversarial
coefficients and their cross-resultant exist.

## 3. Exact acceptance and falsification gates

The adversarial verifier must establish independently:

1. both zeroth residuals vanish for the frozen state;
2. both first coefficients recover the same `a=-B/K` branch;
3. both second coefficients are affine in `c`;
4. their slopes recursively equal the corresponding first-order slopes;
5. their independently formed cross-resultant is exactly
   `129600 epsilon(v)^2/(v^2+4)`;
6. the two full coefficients match the frozen primary artifact only after
   items 1--5 exist.

Use a different exact positivity proof from the primary route.  For
`x=v^2>=0`, prove

```text
z(x)=(x+2)/(2(x+3))>=1/3>cos(2*pi/5).
```

Since `acos` is strictly decreasing on its real principal branch,

```text
acos(z(x))<2*pi/5,
epsilon(x)>0.
```

Thus the cross-resultant has no real zero.  Because the lapse slope is a
nonzero multiple of `v*K` on the registered domain, there is no common `c`.

Hostile controls:

1. `M -> mu(v)+1/13` must give the exact leading lapse defect `-8*pi/13`;
2. replacing `q2` by `c` must fail the exact exponential-jet identity by
   `va+v^3/6` before any residual is evaluated;
3. substituting the lapse root into the momentum coefficient must leave a
   nonzero expression, not an exact zero.

## 4. Independent numerical controls

After the exact derivation and comparison, use 110-decimal direct evaluations
at

```text
v in {-7/5,2/3,5/2},
c in {-2/7,1/5},
h in {1/1200,1/2400,1/4800,1/9600}.
```

Use the full unexpanded action derivatives.  Require either all normalized
errors below `1e-75`, or strict decrease with all three successive base-two
orders in `[0.8,1.2]`.  No sample may determine the global root census.

## 5. Outcome

Only if every exact and numerical gate passes, assign

```text
GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION_ADVERSARIALLY_CORROBORATED.
```

Any disagreement leaves the theorem **OPEN** and becomes the headline.  A
confirmed obstruction is a **DERIVED NEGATIVE, scoped** for the fixed finite
homogeneous action.  It is not a derived tick, speed, unit of time or theorem
about a refined/perfect action.
