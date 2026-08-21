# Adversarial protocol: next-order exceptional-velocity theorem

Date: 2026-08-21

Primary implementation commit: `98acd61`.

Frozen primary artifact:

```text
reproducible/gravity_600cell_generic_velocity_next_order.json
SHA-256 4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18
```

Status: frozen after the primary verifier returned
`GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES` (`13/13`), before constructing or
running an adversarial verifier.

## 1. Claim under attack

For the exact fixed incoming state `(L_minus=1,M=mu(v),p0=p(v))` and endpoint

```text
L_plus=exp(v h+a h^2),
rho=h^2,
v real and nonzero,
```

the first nonzero lapse and momentum residuals are affine in `a`.  Their
generic roots agree.  There is exactly one `x_star>0` such that the common
degree drops at `v=+-sqrt(x_star)`; at those two velocities the lapse constant
is nonzero and no branch exists.  Away from them, duration remains free to
this order.  The primary conditional composition jet has zero endpoint,
final-momentum and action defects.

## 2. Mechanically independent derivative-first route

Do not reuse the primary scaled-action differentiation or its `path_first`
function.

1. Reconstruct the complete unexpanded action.
2. Differentiate it first in `(rho,L_minus,L_plus)` to obtain exact `F`,
   `P_minus` and `P_plus`.
3. Only after differentiation introduce independent positive `tau` and real
   `q` through

   ```text
   L_minus=1,
   L_plus=1+tau*q,
   rho=tau^2.
   ```

4. Reduce the resulting derivative expressions at fixed `(tau,q)` and take
   their partial derivatives at `(tau,q)=(0,v)`.
5. Insert only the endpoint tangent

   ```text
   q(tau)=v+(a+v^2/2)tau+O(tau^2)
   ```

   to obtain the direct total coefficients.  This reverses the decisive order
   of operations used by the primary verifier.
6. Prove the eight positive-radical identities frozen at `71b8312` and use no
   other branch rewrite.
7. Only after the direct coefficients and root census exist, compare them to
   primary artifact hash
   `4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18`.

If exact derivative-first reduction cannot be completed, the theorem remains
unaccepted; numerical agreement alone is insufficient.

## 3. Exceptional-set proof

Independently recover exact functions `K`, `B` and a nonzero prefactor such
that

```text
coefficient_a(C1)=prefactor*K,
constant(C1)=prefactor*B,
constant(C1)*coefficient_a(P1)
 -constant(P1)*coefficient_a(C1)=0.
```

Repeat the monotonicity certificate for

```text
H(x)=(x+3)*sqrt((3x+8)/(x+4))*epsilon(x),
K=sqrt(x+4)*(10-H),
```

including exact derivative identities, positive endpoint inequalities and
the infinite limit.  Use the new numerical bracket

```text
5 < x_star < 6
```

only as a control.  Require agreement with the primary value to 80 decimal
digits.

## 4. New direct controls

Use 100 decimals and points not used by the primary verifier:

```text
v in {-6/5,2/5,11/5},
a in {-1/5,2/9},
h in {1/500,1/1000}.
```

For both unexpanded direct residual quotients, require the same resolved
first-order interval `[0.8,1.2]`, with the same `1e-70` exact-at-precision
rule.  Record every error and order.

At generic control velocities, evaluate the primary common root `a(v)` and
the predicted coarse/mid/final endpoints.  Require the first-slab lapse and
pre-momentum, second-slab lapse, seam momentum and coarse-versus-fine final
momentum/action residuals to converge to their registered zero jets.  These
are independent numerical controls of composition; the exact composition
claim remains separately labelled as primary-route unless rederived exactly.

Use the scaled residuals

```text
4*F_1/h^2,
(p_pre,1-p0)/h,
4*F_2/h^2,
G_mid/h,
(P_plus,2-P_plus,coarse)/h,
(S_1+S_2-S_coarse)/h^2.
```

At `h in {1/500,1/1000}`, each must either be below `1e-70` at both
resolutions or decrease with a halving order in `[0.8,1.2]`.  A mixed
resolved/unresolved pair is disagreement.  No absolute smallness threshold
may substitute for this convergence test.

## 5. Hostile controls

- Change the fixed mass by `+1/10`; the zeroth lapse defect must be
  `-4*pi/5`.
- Change the fixed momentum by `+1/10`; the zeroth momentum defect must be
  `-1/10`.
- Replace `K` by `K+1/10`; its root must move outside an 80-digit agreement
  with `x_star`.

## 6. Outcomes

- `GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED` only
  if every exact coefficient, exceptional-set, hostile and new precision gate
  passes.
- `PRIMARY_GENERIC_NEXT_ORDER_RESULT_REFUTED` if a certified nonzero exact
  difference exists.
- `GENERIC_NEXT_ORDER_ADVERSARIAL_DISAGREEMENT` for any unresolved branch,
  domain, precision or composition control.

No outcome derives a tick.  The corroborated primary interpretation, if it
survives, is: duration is free on the generic velocity domain and the two
isolated velocities are fixed-state discretization obstructions, not selected
finite durations.  Absolute time, `c`, `G`, Planck units and external novelty
remain outside this result.  Only the targeted adversarial verifier will be
run.
