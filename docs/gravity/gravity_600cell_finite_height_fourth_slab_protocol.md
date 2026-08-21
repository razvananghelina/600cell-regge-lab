# Prior-art gate and protocol: fourth slab of the surviving history

Date: 2026-08-21.

Accepted third-slab artifacts:

```text
reproducible/gravity_600cell_finite_height_third_slab.json
SHA-256 6b0e92d031aa891fdc3e1b2045c35bd135a955bb1374c92f015dcd5727d3d8fc

reproducible/gravity_600cell_finite_height_third_slab_adversarial.json
SHA-256 df689f5360ace94d2212e1d71c799ed4e8019457d2702e989bf045ea566abda8
```

Input theorem commit: `9204871`.

Status: frozen before evaluating the third slab's outgoing momentum or
solving any fourth slab.

## 1. Prior-art framing

Repeated discrete pre/post momentum matching is standard variational-
integrator and canonical-Regge methodology; see Marsden--West,
<https://doi.org/10.1017/S096249290100006X>, and Dittrich--Hoehn,
<https://arxiv.org/abs/1108.1974>.  Later simplicial consistency conditions
can remove earlier freedom, but no cited theorem guarantees indefinite
extendibility or global uniqueness for the present nonlinear action.

Published dust-filled 600-cell Regge evolutions also stop after finitely many
steps at a causal boundary; see De Felice--Fabri,
<https://arxiv.org/abs/gr-qc/0009093> and
<https://arxiv.org/abs/gr-qc/0106077>.  Hence a fourth-slab success or stop is
an internal consistency datum, not external novelty by itself.

## 2. Complete frozen history

Use the same homogeneous tetrahedral-frustum 600-cell action at zero
cosmological constant, conserved global dust, positive proper height,
positive endpoint scale and committed pre/post momentum convention.

Fix the history without a new choice:

1. incoming state `v=3/2`;
2. its unique physical first slab;
3. second branch B, selected only because branch A has no physical third
   continuation;
4. the unique physical third slab of branch B.

Reconstruct all three slabs from the complete action at working precision.
The stored decimal roots are comparison data, not substitutes for the direct
solve.

## 3. Fourth incoming state

For the third slab `(h3,q3)` with scale ratio `r3=1+h3*q3`, derive

```text
m3=m2/r3,
pi3=p_post,3/r3^2
   =p(q3)+2*pi*h3*mu(q3)/r3.
```

Require direct recurrence agreement below `1e-110`.  Wrong `r3` power,
reversed sign and mass reset remain hostile controls.

## 4. Complete primary root census

Classify

```text
E4(q)=4*pi[mu(q)-m3]+q[p(q)-pi3]
```

on the entire real line through

```text
E4'(q)=p(q)-pi3.
```

Enumerate all stationary points from every monotone interval of `p`, certify
all finite endpoint signs and both infinite tails, and treat `q=0` directly.
No finite plotting box may support a root count.

For each root reconstruct

```text
h4=[p(q)-pi3]/[2*pi*mu(q)],
r4=1+h4*q.
```

Require reduced equations for every algebraic root.  Require complete-action
and shared-slice residuals below `1e-90` for every physical root, defined by

```text
h4>0,
r4>0.
```

## 5. Frozen outcomes

### `SURVIVING_HISTORY_HAS_UNIQUE_FOURTH_SLAB`

Exactly one physical fourth slab.  Label **DERIVED COMPUTATIONAL,
four-slab scoped / STRUCTURAL**.  Do not infer indefinite evolution.

### `SURVIVING_HISTORY_BRANCHES_AT_FOURTH_SLAB`

More than one physical fourth slab.  The finite-horizon uniqueness at three
slabs does not persist.

### `SURVIVING_HISTORY_STOPS_BEFORE_FOURTH_SLAB`

No physical fourth slab.  The future-extendibility selector only delayed a
finite stopping point.

### `FOURTH_SLAB_EXTENDIBILITY_OPEN`

Incomplete recurrence, tail, exceptional-stratum, provenance or direct-
action gate.

No outcome derives a fundamental tick, infinite history, `c`, `G` or a
Planck scale.  Only targeted verifiers may be run.
