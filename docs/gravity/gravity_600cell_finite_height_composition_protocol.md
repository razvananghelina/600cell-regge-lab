# Preregistration: canonical composition of the finite-height update

Date: 2026-08-21.

Prior-art gate commit: `c9eb996`.

Input theorem commit: `2890c7c`.

Status: frozen before evaluating any outgoing momentum or solving any second
slab.

## 1. Frozen action and conventions

Reconstruct the complete homogeneous cellular action

```text
S(L_minus,L_plus,rho;M)
```

and define

```text
F=rho*partial S/partial rho,
p_pre=-(L_minus/2)*partial S/partial L_minus,
p_post=(L_plus/2)*partial S/partial L_plus.
```

Use only the positive principal branches already certified by the finite
classification.  No sign convention may be changed after a root is seen.

## 2. Exact composition certificates

Before numerical evaluation, prove exactly:

```text
S(lambda*Lm,lambda*Lp,lambda^2*rho;lambda*M)
  =lambda^2*S(Lm,Lp,rho;M),

p_post(lambda data)=lambda^2*p_post(data).
```

For two slabs, differentiate

```text
S_total=S(L0,L1,rho1;M)+S(L1,L2,rho2;M)
```

at the shared scale and prove

```text
(L1/2)*partial S_total/partial L1
  =p_post,1-p_pre,2.
```

This fixes both the sign and the `L1^2` normalization of the second incoming
momentum.

## 3. First-slab image

For every first state in

```text
v in (v_A,v_star) union (v_star,v_C),
```

let `(h1(v),q1(v))` be the already-proved unique physical root and set

```text
l1(v)=1+h1*q1,
m1(v)=mu(v)/l1,
pi1(v)=p_post(1,l1,h1^2;mu(v))/l1^2.
```

Derive `p_post` from the full action before simplifying it on shell.  Record
the exact or certified sign structure of `m1`, `pi1` and the outgoing
Jacobian.  If a complete all-real classification cannot be proved on the
first run, keep it **OPEN**.

## 4. State-curve closure diagnostic

Separately ask whether there is a real `w` satisfying both

```text
mu(w)=m1(v),
p(w)=pi1(v).
```

This diagnostic must report the number of solutions, but it is not the
canonical composition gate.  Failure here may only be labelled failure of
the special one-parameter family to close.

## 5. General second-slab gate

Reconstruct the second normalized slab directly from `(m1,pi1)`.  Prove the
general affine residuals and elimination formula

```text
C2=8*pi[mu(q2)-m1]+4*pi*h2*q2*mu(q2),
P2=p(q2)-pi1-2*pi*h2*mu(q2),

E2=4*pi[mu(q2)-m1]+q2[p(q2)-pi1].
```

Enumerate every real root, including `q2=0` and all zero-slope cases, then
reconstruct

```text
h2=[p(q2)-pi1]/[2*pi*mu(q2)],
l2/l1=1+h2*q2=2*m1/mu(q2)-1.
```

Accept only roots with

```text
h2>0,
l2/l1>0,
det partial(C2,P2)/partial(h2,q2) != 0.
```

Multiple physical roots are not a selected evolution.

## 6. Frozen representative controls

After the exact objects are printed, evaluate exactly these three first
states and no substituted alternatives:

```text
v=3/2    expanding outer-q first branch,
v=3      expanding inner-q first branch,
v=20     contracting negative-q first branch.
```

For each:

1. solve the first full-action equations directly;
2. compute the unnormalised `p_post`;
3. solve the second full-action equations directly in `(h2,q2)`;
4. verify the shared derivative of `S_total` below `1e-90`;
5. require each lapse and momentum residual below `1e-90`;
6. enumerate all competing real second roots through the monotone intervals
   of `p`, not through a finite plotting box.

Numerical representatives may falsify a global statement but may not prove
one.

## 7. Hostile controls

- Replace `p_post/L1^2` by `p_post/L1`; this wrong scale convention must
  change the second equations whenever `L1!=1`.
- Reverse the post-momentum sign; it must fail the shared-slice derivative.
- Reset the conserved mass to `mu(q1)` instead of retaining `M/L1`; this is a
  changed physical state and must be detected.
- A state-curve miss must not be counted as a general composition miss.

## 8. Outcome hierarchy

### `FINITE_HEIGHT_TWO_SLAB_SELF_CLOSED`

Use **DERIVED EXACT, two-slab scoped** only if the outgoing data lie on the
special state curve and the unique general second slab agrees with that
curve-selected slab.

### `FINITE_HEIGHT_TWO_SLAB_GENERAL_COMPOSITION`

Use **DERIVED EXACT / STRUCTURAL, two-slab scoped** if the state curve does not
close but the exact outgoing data select one isolated physical second slab.

### `FINITE_HEIGHT_TWO_SLAB_NONUNIQUE`

Use **STRUCTURAL / OPEN selection** if more than one physical second root
exists for any admitted first state and no canonical selector is derived.

### `FINITE_HEIGHT_TWO_SLAB_NO_GO`

Use **DERIVED NEGATIVE, scoped** only after a complete first-state and
all-real second-root proof shows that no physical second slab exists.

### `FINITE_HEIGHT_TWO_SLAB_OPEN`

Use **OPEN** for incomplete outgoing-image, tail, exceptional-stratum or
root-count proofs regardless of representative agreement.

No outcome derives a universal tick, an absolute time unit, `c`, `G` or
Planck time.  Only the targeted verifier will be run.
