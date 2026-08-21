# Prior-art gate: composition of the selected finite-height update

Date: 2026-08-21.

Input theorem commit: `2890c7c`.

Primary finite-height artifact SHA-256:
`9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03`.

Adversarial finite-height artifact SHA-256:
`da8d60e95b5196beaf93ea234fbf9dfb93e3d5e6bd00fb0a85ed2ef4ba388996`.

Status: prior-art and framing gate before any outgoing-momentum evaluation or
two-slab root solve.

## 1. Exact object and hypotheses

Use the same homogeneous cellular 600-cell Regge-plus-conserved-dust action
at zero cosmological constant,

```text
S(L0,L1,rho;M),
rho=h^2,
h>0,
L0>0,
L1>0.
```

For the first slab fix `L0=1` and an incoming state on the already-classified
one-parameter family

```text
M=mu(v),
p0=p(v),
v in (v_A,v_star) union (v_star,v_C).
```

Let its unique physical root be `(h1,q1)` with `L1=1+h1*q1`.  Define the
canonical data directly from the complete action:

```text
p_pre  =-(L0/2)*partial S/partial L0,
p_post = (L1/2)*partial S/partial L1.
```

The action is homogeneous under

```text
(L0,L1,h,M)->lambda*(L0,L1,h,M),
S->lambda^2*S.
```

Therefore normalize the second slab by `L1`, while holding the physical dust
mass fixed:

```text
m1=M/L1,
pi1=p_post/L1^2.
```

These, not the original `(mu(v),p(v))`, are the incoming canonical data of
the second normalized slab.

## 2. Canonical composition condition

For a second slab with normalized unknowns `(h2,q2)` and
`L2/L1=1+h2*q2`, require

```text
F2=0,
p_pre,2=pi1,
h2>0,
1+h2*q2>0.
```

In unnormalised variables, `p_post,1=p_pre,2` is exactly the discrete
Euler--Lagrange equation at the shared slice.  No momentum, height or endpoint
may be adjusted after this condition is imposed.

The known affine reduction for arbitrary normalized incoming `(m,pi)` is

```text
C=8*pi[mu(q)-m]+4*pi*h*q*mu(q),
P=p(q)-pi-2*pi*h*mu(q).
```

Hence the second-slab elimination equation is expected to be

```text
E_general(m,pi;q)=4*pi[mu(q)-m]+q[p(q)-pi]=0,
```

with

```text
h=[p(q)-pi]/[2*pi*mu(q)],
L_next/L_current=2*m/mu(q)-1.
```

These formulas must be rederived from the complete two-slab action before
they are used.

## 3. Framing correction

There are two different questions:

1. **State-curve closure:** does `(m1,pi1)` equal `(mu(w),p(w))` for some
   real `w`?
2. **Canonical composition:** does the general second-slab problem above
   have an isolated physical root?

The first is a useful structural diagnostic but is not a necessary condition
for the second.  The family `(mu(v),p(v))` was selected as a leading
same-state input family; it has not been proved to be the complete canonical
constraint surface.  Failure of state-curve closure may not be cited as a
failure of evolution.

Conversely, one successful second slab is not yet an indefinitely iterable
flow.  It is only the first composition gate.

## 4. Primary prior art

- [Dittrich--Hoehn, *Canonical simplicial gravity*, arXiv:1108.1974](https://arxiv.org/abs/1108.1974)
  formulate the discrete action as Hamilton's principal function generating
  canonical evolution and explain how later consistency conditions can fix
  data that were initially free.
- [Dittrich--Hoehn, *From covariant to canonical formulations of discrete
  gravity*, arXiv:0912.1817](https://arxiv.org/abs/0912.1817) derive pre/post
  constraints and pseudo-constraints from broken discrete gauge symmetry.
- [Di Bartolo--Gambini--Porto--Pullin, arXiv:gr-qc/0405131](https://arxiv.org/abs/gr-qc/0405131)
  show that discrete consistency conditions can fix continuum Lagrange
  multipliers and generate canonical transformations.
- [De Felice--Fabri, arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093)
  and [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077) implement
  multi-step dust-filled 600-cell evolution and identify its causal stopping
  point.

Thus momentum matching, dynamically fixed discrete lapse data and multi-step
600-cell evolution are **KNOWN** mechanisms.  A successful composition is a
consistency check of the present frozen action, not external novelty by
itself.

## 5. KNOWN / CONTROL / OPEN

### KNOWN

- The first slab has one physical root on the exact registered state set.
- A discrete action supplies pre/post momenta and the shared-slice equation is
  momentum matching.
- Global scale covariance forbids deriving an absolute duration here.

### CONTROL

- Reprove the degree-two action homogeneity and the scaling of `p_post`.
- Obtain `p_post` from the complete action, not from a guessed velocity.
- Verify the shared-slice derivative of `S1+S2` equals the post/pre momentum
  mismatch with the complete sign convention.
- Treat every zero-slope and `q=0` case before division.
- Separate state-curve closure from general canonical composition.
- Require positive height, positive endpoint and a nonzero two-equation
  Jacobian.
- Include direct full-action residuals and a mass-conservation hostile
  control.

### OPEN

- The exact image `(m1(v),pi1(v))` of the first physical branch.
- Whether that image lies on the special state curve.
- Whether every, some or no first-slab state admits a unique physical second
  slab.
- Global root counts of the general second-slab elimination equation.
- Third and later composition, nonhomogeneous stability and refinement.

## 6. Acceptance hierarchy

Use **DERIVED NEGATIVE, scoped** if a complete all-real proof shows that no
first physical state admits a second physical slab.  Use **STRUCTURAL / OPEN
selection** if multiple second roots occur without a derived selector.  Use
**DERIVED EXACT, two-slab scoped** only if the outgoing canonical data and a
unique second physical root are both certified without fitting.

No outcome here is a fundamental tick.  External novelty remains **OPEN**.
Only targeted verifiers may be run.
