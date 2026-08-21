# Prior-art gate and protocol: third-slab future extendibility

Date: 2026-08-21.

Input composition artifact:

```text
reproducible/gravity_600cell_finite_height_composition.json
SHA-256 d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447
```

Input selector theorem commit: `d778e78`.

Status: frozen before evaluating either second branch's outgoing momentum or
solving any third slab.

## 1. Prior-art framing

- Marsden and West, [*Discrete mechanics and variational
  integrators*](https://doi.org/10.1017/S096249290100006X), formulate a
  multi-step discrete trajectory by repeated pre/post momentum matching.  A
  regular discrete Legendre transform gives local evolution; global
  injectivity is stronger.
- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), explain that data left free by a
  simplicial move may be fixed by consistency conditions in later moves.
- Dittrich and Hoehn, [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817), derive preservation conditions
  and pseudo-constraints for later discrete steps.
- De Felice and Fabri, [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093)
  and [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077), evolve
  dust-filled 600-cell geometries through multiple steps and encounter a
  later causal stopping point.

Therefore testing later consistency and extendibility is **KNOWN
methodology**.  A branch distinction at the third slab would be an internal
property of the present frozen action, not by itself a new selection
principle.

## 2. Complete hypotheses and frozen inputs

Use the same homogeneous tetrahedral-frustum 600-cell action at zero
cosmological constant, conserved global dust, positive proper height,
positive endpoint scale and committed pre/post momentum convention.

Fix only the first state

```text
v=3/2
```

and the two already accepted physical second slabs.  Label them before the
calculation by increasing second-slab `q`:

```text
A = lower-q physical second root,
B = higher-q physical second root.
```

Neither branch may be dropped because it is near-null, has a large/small
height or looks unlike a desired cosmology.  The selector audit proved both
are future oriented, causal and locally regular.

## 3. Exact outgoing-state recurrence

For arbitrary normalized incoming state `(m,pi)` and a physical slab root
`(h,q)`, let

```text
r=1+h*q.
```

Derive from the complete action, before evaluating A or B,

```text
m_next=m/r,
pi_next=p_post(1,r,h^2;m)/r^2.
```

Use exact endpoint reversal and action homogeneity as an independent
simplification gate.  The reversed normalized slab has

```text
h_reverse=h/r,
q_reverse=-q,
m_reverse=m/r,
p_pre,reverse=-pi_next.
```

Only if the complete action confirms it may the closed recurrence

```text
pi_next=p(q)+2*pi*h*mu(q)/r
```

be used.  A sign or power of `r` may not be repaired after seeing roots.

## 4. Complete all-real third-root census

For each independently computed outgoing state `(m2,pi2)`, solve

```text
E3(q)=4*pi[mu(q)-m2]+q[p(q)-pi2]=0
```

on the complete real line.  Do not scan a finite plotting box.  Use

```text
E3'(q)=p(q)-pi2
```

and enumerate every stationary point from all four monotone intervals of
`p`:

```text
(-infinity,-v_star), (-v_star,0),
(0,v_star),          (v_star,+infinity).
```

Treat explicitly:

- `q=0`;
- equality with either asymptotic value of `p`;
- roots at stationary points;
- both infinite tails;
- every zero-slope or zero-height case before division.

For every real root reconstruct

```text
h3=[p(q)-pi2]/[2*pi*mu(q)],
r3=1+h3*q=2*m2/mu(q)-1.
```

Accept it as physical only if

```text
h3>0,
r3>0,
det partial(C3,P3)/partial(h3,q)>0.
```

The last determinant must equal the already-proved positive expression, but
direct full-action residuals and the unnormalised shared-slice momentum must
also pass below `1e-90`.

## 5. Hostile controls

- Use `p_post/r` instead of `p_post/r^2`; the third-root state must change.
- Reverse the post-momentum sign; direct junction matching must fail.
- Reset conserved mass to `mu(q)`; this must change the state.
- A finite root box may be printed diagnostically but may not certify the
  all-real count.

## 6. Frozen outcome hierarchy

### `BOTH_SECOND_BRANCHES_EXTEND`

Both A and B have at least one physical third slab.  Branching survives the
first future-extendibility test.  Report each multiplicity; do not call the
evolution deterministic.

### `ONE_SECOND_BRANCH_EXTENDS_UNIQUELY`

Exactly one of A or B has exactly one physical third slab and the other has
none.  Label this **DERIVED COMPUTATIONAL, three-slab scoped / STRUCTURAL
selection candidate**.  It is not yet a physical selection principle.

### `ONE_SECOND_BRANCH_EXTENDS_NONUNIQUELY`

Only one branch extends, but it has multiple physical third slabs.  Global
extendibility prunes one branch while leaving evolution multivalued.

### `NEITHER_SECOND_BRANCH_EXTENDS`

Neither branch has a physical third slab.  The present finite solution is a
two-slab boundary-value relation, not an iterated evolution at this state.

### `THIRD_SLAB_EXTENDIBILITY_OPEN`

Use **OPEN** for an incomplete tail, stationary-point, recurrence,
provenance or direct-action gate.

No outcome derives a fundamental tick, a branch probability, `c`, `G` or a
Planck scale.  Only targeted verifiers may be run.
