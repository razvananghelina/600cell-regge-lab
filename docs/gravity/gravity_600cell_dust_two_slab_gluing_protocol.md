# Preregistration: two-slab 600-cell dust gluing and sign control

Date: 2026-08-16

Prior-art gate: `8c45290`

Pre-protocol framing correction: `620461d`

Consecutive-schedule convention correction: `6c4a377`

Status: **frozen before evaluating a new two-slab action, shared-slice
derivative or orbitwise dust momentum relation**.

This protocol is a structural control.  It contains no deformed evolution
target and performs no root search.

## 1. Exact question and hypotheses

For each of the already derived ordered-schedule parity representatives
`sigma in {even, odd}`, use two consecutive copies of the same five-phase
staircase schedule.  Let

```text
S1 = S_sigma(q0, x1, q1),
S2 = S_sigma(q1, x2, q2),
Ssum = S1 + S2,
```

where each one-slab action is the certified complete Lorentzian Regge
curvature action plus the published De Felice--Fabri dust action.  The three
spatial layers contain vertices

```text
layer 0:   0,...,119,
layer 1: 120,...,239,
layer 2: 240,...,359.
```

The first slab uses layers `0,1`; the second uses layers `1,2`.  Both use the
same phase function on logical vertex labels modulo 120.

Freeze the published time-symmetric control on both factors:

```text
q0 = q1 = q2 = l0^2 on every spatial edge,
x1 = x2 = (d^2 on 30 diagonal orbits, tau^2 on 5 pole orbits),
tau = 0.0102,
d^2 = l0^2 - tau^2,
M = (90/pi) (2*pi-5*acos(1/3)) l0.
```

These values are external controls.  No number may be optimized.

## 2. Geometric construction and orbit map

Rebuild the 600-cell, its 600 tetrahedra, the five independent cells of 24
vertices, both schedule representatives and both 2400-simplex staircase
slabs from the same exact code path certified upstream.

The combined complex must contain exactly:

```text
360 vertices in three labelled layers,
4800 four-simplices,
one shared 600-cell layer,
```

with no duplicated four-simplex.

The order-24 stabilizer acts on logical labels in all three layers.  Derive
the first-final to second-old edge map only by

```text
(u+120,v+120) in first final layer
    ->
(u,v) in a local copy of the second old layer,
```

equivalently by common logical endpoints `{u,v}`.  Map whole orbit sets and
require each of the 30 first-final orbits to land in exactly one second-old
orbit.  The result must be a bijection.  Sorting momentum values or searching
over `30!` permutations is forbidden.

Independently derive the old/new layer-reversal map within one slab and record
its induced 30-orbit permutation.  It is an audit of time-reversal signs, not
the schedule used for the second forward slab.

## 3. Direct combined action

Construct a direct two-slab Regge evaluator on the union of all 4800
four-simplices.  Its hinge constants are:

```text
pi    on triangles in outer layers 0 and 2,
2*pi  on every other triangle, including layer 1.
```

The direct dust action is the sum of the two five-pole world-line terms.  No
dust term is attached to spatial boundary edges.

At the published control require:

1. the direct two-slab action and `Ssum` agree;
2. every simplex in both evaluations has one timelike Gram direction;
3. all angle-argument and minimum-Gram branch gates inherited from the
   published control pass;
4. imaginary contamination is below `1e-70` at 100 decimal digits.

Also check action gluing at all sixty one-coordinate shared-boundary audit
points

```text
q1_j -> q1_j exp(+1e-6),
q1_j -> q1_j exp(-1e-6),       j=1,...,30.
```

This is an exhaustive orbitwise gluing audit, not a fitted sample.  At every
point require direct-versus-summed action agreement within the inherited
upstream `5e-8` relative action class and require identical branch counts.

## 4. Arbitrary-precision boundary derivatives

Use 100 decimal digits for variables, actions and contractions.  For each of
the 30 old and 30 new boundary orbit squares of a one-slab factor compute

```text
m_j(h) = [S(q_j exp(h))-S(q_j exp(-h))]/(48 h).
```

The division by 24 reports the per-edge logarithmic orbit derivative.  Define

```text
p_pre  = -m_old,
p_post = +m_new.
```

Use the already calibrated disjoint derivative pairs without alteration:

```text
operational: (1e-20, 1e-15),
validation:  (3e-20, 3e-15).
```

In each pair the smaller-step row is primary and the signed primary-minus-
shadow difference is the stability proxy.  Require operational/validation
agreement componentwise within ten times the sum of their proxies, with
arithmetic floor `1e-60`.  Require all branch audits and maximum imaginary
contamination below `1e-70`.

Use the geometrically derived reversal permutation `R` and preregister the
time-symmetric one-slab sign relation

```text
p_pre + R p_post = 0.
```

It must pass componentwise inside the same ten-proxy band.  No equality of
sorted multisets counts.

## 5. Shared-slice derivative and negative sign control

Compute the direct two-slab derivative with respect to every shared orbit by
the same operational and validation pairs:

```text
r_shared,j = (1/24) partial S_direct / partial log(q1_j).
```

Using the consecutive-slab orbit bijection `P`, require componentwise

```text
r_shared = p_post(first) - P p_pre(second).
```

At the repeated published sandwich the two factor records are identical in
their own old/final orbit coordinates, but the consecutive-slab gluing map
`P` is not assumed to equal the inverse of the time-reversal audit map `R`.
Consequently no `2 p_post` simplification is allowed.  The preregistered
negative-control vector is exactly

```text
cusp = p_post(first) - P p_pre(second).
```

It must be assembled from the two independently evaluated one-slab boundary
derivatives and the incidence-derived map `P`, before the direct two-slab
derivative is inspected.  Require both:

1. `r_shared` agrees with `cusp` inside ten times the combined
   direct/component stability proxies;
2. `norm(cusp) > 100 max(norm(error_proxy),1e-60)` and the independently
   evaluated `r_shared` passes the same nonzero gate.

Thus a zero shared residual is a **failure**, not a positive result.  It would
signal an incorrect boundary sign, orbit map or gluing convention.

## 6. Frozen attempts and outcome hierarchy

There are exactly two attempts:

```text
even repeated schedule,
odd repeated schedule.
```

No phase-order search, alternative orbit permutation, branch choice,
deformation, optimizer or tolerance adjustment is allowed.

The only outcomes are:

- `TWO_SLAB_GLUING_CONTROL_PASSED` if both parities pass every geometry,
  branch, direct-action, derivative, orbitwise sign and nonzero-cusp gate;
- `TWO_SLAB_GEOMETRY_OR_ORBIT_MAP_FAILED` if the combined complex or either
  derived permutation fails;
- `TWO_SLAB_ACTION_GLUING_FAILED` if direct and summed actions disagree;
- `TWO_SLAB_DERIVATIVE_CONTROL_FAILED` if either derivative pair fails its
  accuracy or branch gates;
- `TWO_SLAB_MOMENTUM_SIGN_CONTROL_FAILED` if an orbitwise time-reversal or
  shared-derivative identity fails.

Stop at the first applicable failure class and record all evaluations already
performed.  Do not repair a failure in the same experiment.

## 7. Acceptance, kill and claim boundary

**Acceptance boundary:** both parity representatives return
`TWO_SLAB_GLUING_CONTROL_PASSED`, including the nonzero repeated-sandwich
shared residual.  This is **DERIVED CONTROL** only.

**Kill boundary:** failure of geometry, action gluing or orbitwise momentum
signs blocks the present canonical-evolution route.  Do not search for a next
slab until a separately preregistered correction identifies the exact
convention error or the construction is abandoned.

No outcome here derives an evolving next frame.  In particular it does not
derive:

- a second-slab solution;
- expansion or contraction;
- a tick duration;
- lapse selection;
- a causal speed;
- a continuum limit;
- `c`, Planck time, Planck mass or particle masses.

If and only if this control passes, the next separate protocol may test the
canonical map

```text
given q_old and p_target,
solve internal equations and p_pre=p_target for x and q_new,
```

first with `p_target=p_pre(published)` as a reproduction control, and only
then with `p_target=p_post(published)` as the first candidate forward
continuation.  The known collective lapse/gauge direction must be handled
explicitly and must not be mistaken for a physical time scale.
