# Frozen adversarial protocol: direct two-slab nonuniqueness audit

Date: 2026-08-21.

Primary protocol commit: `386466e`.

Registered primary implementation: `c7bffc3`.

Preserved symbolic-timeout correction: `1fe35f4`.

Primary correction implementation: `ac6527a`.

Primary artifact commit: `0f832c6`.

Primary artifact SHA-256:
`d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447`.

Status: frozen before the adversarial verifier exists or is executed.

## 1. Independence boundary

The primary composition proof eliminated height to the scalar equation

```text
E2(q)=4*pi[mu(q)-m1]+q[p(q)-pi1]
```

and counted its roots through the monotone branches of `p`.  Reusing that
elimination is only reproducibility.

The adversarial verifier must:

1. redifferentiate the complete action;
2. solve the original two equations directly in `(h,q)`;
3. compute `p_post` from the complete derivative;
4. construct the unnormalised two-slab action and verify its shared-slice
   derivative directly;
5. compare with the primary artifact only after both second roots exist.

No scalar determinant or `E2` root finder may be used in the decisive direct
solve.

## 2. Frozen counterexample state and seeds

Use only the preregistered admitted state

```text
v=3/2.
```

At each of 80, 120 and 180 decimal digits, solve the first full-action
constraint and incoming momentum equations from the rational seed

```text
(h1,q1)=(1/5,10).
```

Compute `(L1,m1,pi1)` independently.  Then solve the second full-action
equations from the two disjoint rational seeds

```text
branch A: (h2,q2)=(7,1/50),
branch B: (h2,q2)=(1/14,31).
```

Also solve from

```text
reverse control: (h2,q2)=(-1/14,10),
```

using the signed affine extension only as a control; it must have `h2<0` and
is not a physical slab.

The two physical roots must remain distinct under independent seed changes
of `+-5%` in both coordinates.  A seed basin is not evidence unless all
original residuals and Jacobians pass.

## 3. Acceptance tests

For each of branches A and B require:

```text
h2>0,
1+h2*q2>0,
abs(C2)<1e-90,
abs(P2)<1e-90,
abs(det d(C2,P2)/d(h2,q2))>1e-20,
abs((L1/2)*d(S1+S2)/dL1)<1e-90.
```

Require the two roots to differ by more than `1` in `q2` and to produce
different next scales.  Reconstruct them independently at all three
precisions and require agreement to at least 60 decimal digits.

## 4. Convention attacks

- `p_post/L1` instead of `p_post/L1^2` must move both roots or destroy the
  matching residual.
- `-p_post/L1^2` must fail the shared-slice sign.
- Replacing conserved `M/L1` by `mu(q1)` must change the second equations.
- The state-curve closure diagnostic is not used anywhere in the direct
  solve.

## 5. Outcome boundary

If both independent physical roots survive, emit

```text
FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED.
```

This is **DERIVED NEGATIVE for unique deterministic composition under the
frozen action**, and **STRUCTURAL / OPEN selection** for evolution: a further
principle might select one branch, but none is presently derived.

Any disagreement returns the result to `FINITE_HEIGHT_TWO_SLAB_OPEN`.
Nothing here derives a fundamental tick or excludes other canonical state
spaces, boundary conditions or improved actions.
