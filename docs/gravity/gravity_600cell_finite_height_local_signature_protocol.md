# Protocol: local stability of the `v=3/2` branch signature

Date: 2026-08-22.

Prior-art gate commit: `5bace00`.

Status: frozen before the local verifier exists or any Arb certificate is
evaluated.

## 1. Claim and evidential boundary

Under the complete hypotheses below, certify that there exists an
**unspecified** real `epsilon>0` such that every incoming state on the analytic
curve

```text
(m,pi)=(mu(v),p(v)),  |v-3/2|<epsilon,
```

has the same ordered physical branch tree through slab four as the accepted
representative at `v=3/2`:

```text
slab 1: one nontrivial physical child;
slab 2: two physical children A<B in increasing q;
branch A at slab 3: no physical child (DEAD);
branch B at slab 3: one physical child;
branch B at slab 4: one physical child whose state enters D.
```

Here

```text
D={(m,x): 0<m<=2/5, x=m*q>=125}.
```

This is a local topological-stability theorem. It does not derive `v=3/2`,
select an explicit radius, classify the complete incoming domain, or promote
complete extendibility to a local physical law.

## 2. Complete hypotheses and equations

Use only:

- the fixed homogeneous tetrahedral-frustum 600-cell action;
- zero cosmological constant;
- conserved global dust;
- the committed pre/post canonical-momentum convention;
- positive proper height and positive endpoint scale;
- the accepted invariant theorem on `D`;
- the already derived real-analytic scalar functions

  ```text
  epsilon(q)=2*pi-5*acos[(q^2+2)/(2*(q^2+3))],

  mu(q)=180*epsilon(q)/(pi*sqrt(q^2+4)),

  p(q)=180*q*epsilon(q)/sqrt(q^2+4)
       -600*sqrt(3)*asinh[q/sqrt(8*(q^2+3))].
  ```

At every normalized canonical state `(m,pi)`, use

```text
E(m,pi,q)=4*pi*(mu(q)-m)+q*(p(q)-pi),

h=[p(q)-pi]/[2*pi*mu(q)],
r=2*m/mu(q)-1=1+h*q,

m_plus=m/r,
pi_plus=p(q)+[p(q)-pi]/r.
```

The physical gates are exactly `h>0` and `r>0`.

## 3. Symbolic analytic skeleton

Before numerical certification, verify symbolically from the displayed
functions that

```text
partial_q E = p(q)-pi.
```

On the initial curve verify identically

```text
E(mu(v),p(v),v)=0,
partial_q E(mu(v),p(v),v)=0.
```

The diagonal solution is therefore a persistent zero-height stationary root,
not a physical child. Certify `p'(3/2)!=0`, so this tangency is nondegenerate
and persists as exactly one double root locally. Every other root used in the
tree must be simple.

The real functions are analytic on a neighbourhood of every certified finite
argument. No statement about complex analyticity is required.

## 4. Complete real-root certificate at `v=3/2`

Use python-flint Arb at 192 decimal digits or more. Every printed sign must be
an outward-rounded ball excluding zero. Decimal `mpmath` values may seed
brackets but are not certificates.

At each visited state, certify the complete real-root census by the accepted
stationary partition

```text
E_q=p(q)-pi.
```

The certificate must include:

1. rigorous comparison of `pi` with `p_star`, `p_infinity` and `0`, fixing
   the complete number of stationary points on both axes;
2. disjoint Arb enclosures for every stationary point, with `p'(q)` excluding
   zero;
3. `E` at every stationary point, with zero excluded except for the exact
   initial diagonal tangency;
4. nonzero right and left tail coefficients;
5. nonzero `E(0)`;
6. one disjoint interval enclosure for every real root implied by the ordered
   sign partition;
7. `E_q` excluding zero at every non-diagonal root.

Initial rational brackets may be taken from the committed discovery artifact,
but every bracket must be certified internally and the complete count must
come from stationary signs and analytic tails, not from the bracket list.

The required ordered all-real/physical counts are frozen as

```text
initial state:       2 all-real including q=v, 1 nontrivial physical;
second state:        3 all-real, 2 physical;
branch-A third:      2 all-real, 0 physical;
branch-B third:      3 all-real, 1 physical;
branch-B fourth:     3 all-real, 1 physical.
```

If another real root exists, a listed root is not unique, or a stationary/tail
sign is unresolved, the outcome is `OPEN`.

## 5. Strict gates and recursive state

For every real root, not only physical ones, serialize Arb enclosures for

```text
h,
r,
mu(q)-2*m,
E_q.
```

Certify the expected signs without tolerance. Recursively compute each
accepted outgoing state using Arb enclosures and prove that the branch ordering
is disjoint.

For the branch-B slab-four entry, certify strictly

```text
0<m<2/5,
q>0,
m*q>125.
```

Equality is not accepted. Once these inequalities hold, invoke only the
already accepted invariant theorem for later continuation.

## 6. Local-constancy lemma

The verifier and result note must state the complete logical implication:

- the nondegenerate initial diagonal tangency persists identically;
- every other root continues uniquely by the real-analytic implicit-function
  theorem;
- strict stationary, origin and tail signs prevent an unseen real root from
  appearing locally or escaping from infinity;
- strict `h`, `r`, branch-ordering and `D` inequalities preserve every
  physical/terminal label;
- the outgoing-state formulas are analytic because every accepted `r`
  excludes zero.

These facts imply existence of some `epsilon>0`. They do not compute a maximal
or physically meaningful radius. Do not perform a radius search, grid scan,
node fraction or interval expansion.

## 7. Frozen controls

After the primary certificate is assembled in memory:

- compare its ordered root midpoints with the committed discovery tree only
  as a reproducibility control;
- redifferentiate the complete slab action and verify constraint, pre-momentum
  and post-momentum residuals on every physical edge of the tree;
- require the wrong post-momentum scale and reversed sign to fail at least one
  recursive-state enclosure;
- replace `m*q>125` by the false hostile inequality `m*q>126` at the accepted
  entry and require rejection.

Controls cannot repair a failed primary sign certificate.

## 8. Outcomes

### `LOCAL_SIGNATURE_PRIMARY_CERTIFIED`

All symbolic identities, complete Arb root censuses, strict gates, recursive
states and controls pass. This establishes the local theorem only as a
**primary certificate**. Rule 4 still requires a separately preregistered,
mechanically different adversarial replication before consolidation.

### `LOCAL_SIGNATURE_OPEN`

Use for any interval dependency loss, overlapping enclosure, stationary/tail
ambiguity, zero gate, root-count disagreement or failed control. Preserve the
artifact and do not introduce an explicit fitted radius as a repair.

## 9. Interpretation boundary

A corroborated positive result would refute only the claim that the accepted
history is literally confined to one isolated incoming point. Because no
radius is selected and the `D` margins are small, it would not justify the
words `generic`, `large basin`, `physical selection` or `new local physics`.

