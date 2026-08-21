# Result: one two-slab branch has a unique physical third extension

Date: 2026-08-21.

## Provenance

```text
third-slab prior-art gate and protocol            0f31fe8
primary verifier registered before first run      beac60d
first primary precision failure preserved         564dca9
primary bisection correction frozen               9c40419
corrected primary implementation                  cceb0cc
accepted primary artifact                         a0abf0a
dual adversarial protocol                         2f1c6f0
adversarial verifier registered before first run  36643b1
first dual failure preserved                      8bb6b39
dual bracket/domain correction frozen             e1898f6
corrected dual implementation                     5439781
accepted adversarial artifact                     8d41256
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_finite_height_third_slab.py
  8/8 PASS

reproducible/verify_gravity_600cell_finite_height_third_slab_adversarial.py
  9/9 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_finite_height_third_slab.json
SHA-256 6b0e92d031aa891fdc3e1b2045c35bd135a955bb1374c92f015dcd5727d3d8fc

reproducible/gravity_600cell_finite_height_third_slab_adversarial.json
SHA-256 df689f5360ace94d2212e1d71c799ed4e8019457d2702e989bf045ea566abda8
```

Both targeted verifiers were rerun together and reproduced these hashes byte
for byte.  No full-suite run was performed.  Static registry audit: `411/411`
distinct registrations, two deliberate exclusions, zero duplicates, zero
unregistered verifiers and zero missing files.

## Headline

```text
ONE_SECOND_BRANCH_EXTENDS_UNIQUELY_ADVERSARIALLY_CORROBORATED
```

> **DERIVED COMPUTATIONAL, THREE-SLAB SCOPED / STRUCTURAL SELECTION
> CANDIDATE / ADVERSARIALLY CORROBORATED:** of the two physical second slabs
> from the frozen first state `v=3/2`, branch A has no physical third-slab
> continuation, while branch B has exactly one.  Requiring one further
> physical slab therefore separates the two branches at this finite horizon.

This does not prove that future extendibility is a fundamental physical
selection axiom, nor that branch B extends indefinitely.  It is the first
derived structural condition in the present homogeneous action that
distinguishes the two previously admissible continuations.

## Complete hypotheses

Use the fixed homogeneous tetrahedral-frustum 600-cell action at zero
cosmological constant, conserved global dust, the committed canonical
pre/post momentum convention, positive proper height and positive endpoint
scale.  Fix the admitted first state `v=3/2` and retain both physical second
slabs exactly as classified by the two-slab theorem.

The result does not cover another initial state, nonhomogeneous edge data,
another action or carrier, refinement, a cosmological constant, a quantum
boundary state or an added matter clock.

## Exact outgoing-state recurrence

For a normalized incoming state `(m,pi)` and any physical slab `(h,q)`, put

```text
r=1+h*q.
```

Action homogeneity and endpoint reversal give

```text
m_next=m/r,
pi_next=p_post(1,r,h^2;m)/r^2
       =p(q)+2*pi*h*mu(q)/r.
```

The verifier first reconstructed the complete action and proved its endpoint
symmetry.  It then confirmed the closed recurrence directly on both second
branches with errors below `4e-149`.  Wrong powers of `r`, reversed momentum
sign and mass reset all change both outgoing states.

## Primary complete real-line proof

For each outgoing state `(m2,pi2)`, the primary route classified

```text
E(q)=4*pi[mu(q)-m2]+q[p(q)-pi2].
```

Its decisive exact identity is

```text
E'(q)=p(q)-pi2.
```

All stationary points were enumerated from the four monotone intervals of
`p`, and every interval between them was certified by endpoint and infinite-
tail signs.  No finite plotting box was used.  The point `q=0` and all tail
coefficients were evaluated separately.

## Adversarial dual proof

The independent route solved the constraint first for `q!=0`:

```text
h_C(q)=2[m2-mu(q)]/[q*mu(q)],
R(q)=p(q)-pi2+4*pi[mu(q)-m2]/q.
```

It did not use `E` or `E'=p-pi2`.  Instead it derived

```text
R'(q)=4*pi[m2-mu(q)]/q^2.
```

Thus the adversarial monotone partition is controlled by every equal-`mu`
point.  It handled the two one-sided limits at zero, both infinite tails and
the excluded `q=0` point, then reconstructed every root from the complete
action at 100 and 160 decimal digits.  The censuses nested beyond 60 digits
before the primary artifact was read.

## Complete frozen census

### Branch A

Outgoing state:

```text
m2  = 3.43466651778226612123913325384932...,
pi2 = 141.732061579487411233398103800390....
```

| `q3` | `h3` | `L3/L2` | Physical failure |
|---:|---:|---:|---|
| `-2.36620856...` | `0.62644994...` | `-0.48231121...` | negative endpoint scale |
| `0.02125943...` | `-6.31126685...` | `0.86582604...` | negative height |

Therefore branch A has two real algebraic roots and zero physical third
slabs.

### Branch B

Outgoing state:

```text
m2  = 1.25686888697869468141790509851156...,
pi2 = -171.794841746295772087674469515661....
```

| `q3` | `h3` | `L3/L2` | Physical failure |
|---:|---:|---:|---|
| `-0.16186931...` | `8.35800604...` | `-0.35290466...` | negative endpoint scale |
| `31.27922362...` | `-0.02184079...` | `0.31683711...` | negative height |
| `99.62760169...` | `0.02184095...` | `3.17596147...` | none |

The last root is isolated, has positive local determinant and passes both
complete-action residuals and the shared-slice junction below `9e-116`.
Therefore branch B has three real algebraic roots and exactly one physical
third slab.

## Preserved failures and corrections

The first primary run found the same provisional counts but stopped at
`6/8`: the bisection routine used the `1e-90` topology sign margin inside a
root iteration whose requested tolerance was `1e-115`.  It consequently
stopped at residuals `+-2e-90`.  The failure was committed before the raw-sign
bisection correction was frozen.

The first adversarial run stopped at `5/9`.  Its `+infinity` helper chose
`q=10` for an interval beginning at `q=47.637...`, duplicated an interior
root and missed the exterior root.  It also incorrectly demanded a positive-
orientation full-action check for algebraic roots with negative `h`, despite
the action using `rho=h^2` and `sqrt(rho)>0`.  Both errors were visible in the
preserved artifact.  The frozen correction required tail points to lie beyond
their adjacent stationary point and reserved direct-action checks for
physical roots, while retaining reduced-equation checks for every algebraic
root.  No root interval, physical condition or outcome rule changed.

## Status ledger

| Claim | Status |
|---|---|
| The outgoing normalized state recurrence is correct | **DERIVED / direct-action checked** |
| Branch A has no physical third slab at `v=3/2` | **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED** |
| Branch B has exactly one physical third slab | **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED** |
| One-step future extendibility separates A and B | **DERIVED STRUCTURAL, three-slab scoped** |
| Future extendibility is a fundamental selection axiom | **OPEN** |
| The B history extends to a fourth slab | **OPEN** |
| The B history extends indefinitely | **OPEN** |
| The resulting finite-horizon history is stable | **OPEN** |
| The result persists for other admitted first states | **OPEN** |
| The result persists nonhomogeneously or under refinement | **OPEN** |
| A deterministic fundamental tick is derived | **NO** |
| External novelty of this coefficient-level branch census | **OPEN** |

## Interpretation and next gate

At two slabs the action was a globally multivalued canonical relation.  At
three slabs, demanding one more physical future removes the slow A branch and
leaves one B continuation.  This is more informative than choosing a branch
by size, causality margin or local regularity: it follows from the next
canonical consistency problem itself.

The next non-fitted test is the fourth slab of the unique surviving history.
If it has no physical continuation, the apparent selector only delays the
stopping point.  If it has one, the finite history extends uniquely by one
more step.  If it branches again, finite-horizon uniqueness is lost.  No
claim about a dynamical law should outrun that census.
