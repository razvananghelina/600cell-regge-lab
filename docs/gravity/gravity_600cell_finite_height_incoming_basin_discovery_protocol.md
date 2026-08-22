# Protocol: target-free incoming-state branch skeleton

Date: 2026-08-22.

Prior-art gate commit: `edd5db0`.

Status: frozen before the discovery verifier exists or any new incoming-state
node is evaluated.

## 1. Purpose and evidential boundary

Construct a deterministic candidate skeleton for the complete branch diagram
over the already derived one-slab physical incoming domain

```text
I=(v_A,v_star) union (v_star,v_C).
```

This stage is **DISCOVERY ONLY**.  A finite node census cannot prove a
continuum branch count, exclude a narrow bifurcation or establish a basin.
Its only admissible output is a preregistered list of candidate signature
cells and intrinsic equations for a later rigorous interval proof.

The representative results at `v=3/2`, `v=3` and `v=20` are already known and
cannot be made blind.  They must not choose the grid, depth, brackets or
signature rules.  Compare them only after the candidate skeleton has been
serialized.

## 2. Complete hypotheses and exact state update

Use only:

- the fixed homogeneous tetrahedral-frustum 600-cell action;
- zero cosmological constant;
- conserved global dust;
- the committed pre/post canonical-momentum convention;
- positive proper height and positive endpoint scale;
- the accepted exact one-slab thresholds and invariant half-strip.

Define

```text
epsilon(q)=2*pi-5*acos[(q^2+2)/(2*(q^2+3))],

mu(q)=180*epsilon(q)/(pi*sqrt(q^2+4)),

p(q)=180*q*epsilon(q)/sqrt(q^2+4)
     -600*sqrt(3)*asinh[q/sqrt(8*(q^2+3))].
```

For a normalized canonical state `(m,pi)`, every nonzero candidate slab root
satisfies

```text
E(m,pi,q)=4*pi*(mu(q)-m)+q*(p(q)-pi)=0,

h=[p(q)-pi]/[2*pi*mu(q)],
r=1+h*q=2*m/mu(q)-1,

m_plus=m/r,
pi_plus=p(q)+[p(q)-pi]/r.
```

The physical gates are exactly

```text
h>0,
r>0.
```

At the initial state use

```text
m0=mu(v),
pi0=p(v),
```

and exclude the diagonal zero-height root `q=v`.  Its unique nontrivial
physical first root is already certified by the one-slab theorem.

## 3. Complete all-real root census at each node

No finite `q` box or unconstrained polynomial root finder is permitted.  Use
one of the already certified complete scalar partitions:

```text
E_q=p(q)-pi
```

with stationary points given by all solutions of `p(q)=pi`, or

```text
R(q)=p(q)-pi+4*pi*(mu(q)-m)/q,
R'(q)=4*pi*(m-mu(q))/q^2
```

with `q=0` checked separately.  Use analytic tails at both infinities and
certified bisection between consecutive stationary points.  Every real root
must be recorded before the physical gates are applied.

Near a multiple or endpoint root, return `UNRESOLVED`; do not assign a count
through a tolerance.

## 4. Frozen discovery nodes

Split the incoming domain only at the intrinsic accepted thresholds:

```text
(v_A,v_star),
(v_star,v_M),
(v_M,v_C).
```

On each open component evaluate exactly `1024` Gauss--Chebyshev nodes

```text
v_j=(a+b)/2+(b-a)/2*cos[(2*j+1)*pi/(2*1024)],
j=0,...,1023,
```

sorted increasingly and computed at 100 decimal digits with 40 guard digits.
The irrational node locations are deterministic diagnostics, not interval
endpoints for the later proof.

Do not add endpoint probes, random seeds, adaptive nodes or a denser fallback.
If any root census is unresolved, preserve it and return discovery `OPEN`.

## 5. Frozen branch depth and terminal labels

From every incoming node construct the complete physical branch tree through
slab four:

```text
q1, q2, q3, q4.
```

At each slab, sort all physical children by increasing `q`.  A branch is
terminally labelled:

- `DEAD` if it has no physical child before or at slab four;
- `ENTERED_D` as soon as a current root satisfies

  ```text
  0<m<=2/5,
  x=m*q>=125;
  ```

- `LIVE_OUTSIDE_D_AT_DEPTH_4` if it remains physical at `q4` without entering
  `D`;
- `UNRESOLVED` if any multiple-root, endpoint, precision or census ambiguity
  occurs.

Once `ENTERED_D` is assigned, do not iterate further: the rigorous invariant
theorem supplies its unique later continuation.

Depth four is frozen because it is the first depth at which the already known
`v=3/2` branch enters `D`.  This is disclosed post-hoc motivation, not a
derived universal horizon.  If any branch is `LIVE_OUTSIDE_D_AT_DEPTH_4`, the
complete basin remains **OPEN**; do not increase the depth in this mission.

## 6. Signature and candidate change cells

For each incoming node serialize the ordered rooted-tree signature containing:

- all-real and physical child counts at every visited state;
- the ordered physical child `q,h,r,m_plus,pi_plus,x` values;
- terminal labels and first entry depth;
- all nonzero Jacobian and noncritical-level margins used numerically.

Compress only consecutive incoming nodes with identical combinatorial tree
signature and terminal-label multiset.  Every boundary between unequal
signatures is a candidate change cell bounded by its two adjacent frozen
nodes.

Within constant-signature runs, also create candidate cells whenever one of
the following monitored quantities changes sign between adjacent nodes:

```text
h,
r,
m-2/5,
m*q-125,
E at a stationary point,
p(q_stationary)-pi,
mu(q)-2*m.
```

Do not interpret linear interpolation inside a cell as a threshold.

## 7. Intrinsic equations for the later proof

Each candidate cell must be assigned one or more equation types, with all
variables and the complete branch path included:

1. branch birth/merger:

   ```text
   E(m,pi,q)=0,
   p(q)=pi;
   ```

2. zero height / state-curve contact:

   ```text
   E=0,
   p(q)=pi,
   mu(q)=m;
   ```

3. zero endpoint:

   ```text
   E=0,
   mu(q)=2*m;
   ```

4. invariant entry:

   ```text
   m=2/5
   ```

   or

   ```text
   m*q=125;
   ```

5. inherited incoming endpoint or puncture:

   ```text
   v=v_A, v_star, v_M, or v_C.
   ```

The later rigorous protocol must freeze the exact candidate list and a
complete interval-cover strategy before solving or comparing any threshold.
No candidate may be deleted merely because its numerical margin looks small.

## 8. Frozen controls

Only after the discovery artifact has been written in memory, check:

- `v=3/2` has two physical second branches, with one dying and the other
  entering `D` at the already accepted depth;
- `v=3` has two physical second branches;
- `v=20` has no physical second successor;
- reversed post-momentum sign, `p_post/r` instead of `p_post/r^2`, and dust
  mass reset change the branch signature;
- the exact action residuals reproduce one node per distinct signature at
  140 decimal digits.

These are controls, not evidence for completeness.

## 9. Outcomes

### `INCOMING_BASIN_CANDIDATE_SKELETON_FROZEN`

Every frozen node has a complete scalar root census, the exact-action and
hostile controls pass, and a finite candidate-cell list is serialized.  This
outcome permits preregistration of the later interval proof but proves no
continuum basin.

### `INCOMING_BASIN_DISCOVERY_OPEN`

Use for any unresolved node, resource limit, incomplete tail partition or
control disagreement.  Preserve all resolved rows and do not change node
count, precision or depth.

No discovery outcome supports a fraction-of-nodes claim, probabilistic basin
measure or physical selection rule.

## 10. Scope warning

Even a later successful continuum theorem over `I` classifies only the
special one-parameter incoming curve `(mu(v),p(v))`.  It does not classify the
full two-dimensional canonical state space `(m,pi)`, introduce local degrees
of freedom or derive an absolute tick.
