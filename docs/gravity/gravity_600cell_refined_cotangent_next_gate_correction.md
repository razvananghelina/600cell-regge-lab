# Correction: temporal reflection does not select the missing refined momentum

Date: 2026-08-21

## Disclosed framing error

The first consolidated cotangent-lift note proposed a
"time-reflection-symmetric two-slab boundary-value problem" as an example of
an action-derived way to select refined initial data.  That example was too
strong as written.

The repository had already derived, on the published regular dust sandwich,

```text
p_pre + R p_post = 0
```

under the geometric time-reversal map `R`, and for a shared boundary

```text
dS_total/dq_shared = p_post - P p_pre.
```

The individual time-symmetric sandwich momenta are nonzero:

```text
p_pre,j  = -0.00090810444890653157621005256593...,
p_post,j = +0.00090810444890653157621005256593....
```

Thus temporal reflection supplies a relation between pre- and post-momenta,
not an absolute zero covector and not a unique inverse of the coarse/fine
pullback.  A slab and its correctly reversed partner can obey the reflection
relation while retaining the five refined cotangent directions found in the
new theorem.  Treating time symmetry as `p_refined=0` would silently choose a
boundary canonical convention and would not solve the underdetermination.

The established source is
`gravity_600cell_dust_two_slab_gluing_result.md`, whose action and momentum
identities were already independently certified.  No new numerical result is
claimed here.

## Corrected next gate

The immediate action-derived seed remains the refined analogue of the
published equal-boundary stationary sandwich itself.  Its internal equations
must possess an on-shell `H4` fill before its six boundary momenta can be used
as refined canonical data.

Consequently the paused nested `6+3+1` census is now load-bearing:

- a certified finite root would supply a refined stationary fill whose
  boundary derivatives are selected by the refined action;
- a certified no-root-on-the-followed-branch result would close that specific
  equal-boundary seed and force the programme toward a genuine perfect-action
  coarse graining or a newly derived matter/time discretization;
- a branch-unresolved result would leave the seed **OPEN**, not license a
  pseudoinverse lift.

The earlier complete execution was user-witnessed with no accepted root but
lost its artifact.  The deterministic checkpoint contains `2/12` classes.
At the user's explicit direction the rerun is deferred until exclusion of
this branch is required for a concrete decision; the witnessed result is not
promoted in the meantime.

## Status

- **DERIVED EXACT:** the five-dimensional lift freedom remains unchanged.
- **REFUTED AS A SHORTCUT:** temporal reflection alone removes that freedom.
- **OPEN / DEFERRED:** existence of an equal-boundary refined on-shell seed.
- **OPEN:** action-selected refined dynamics, tensor modes, dispersion, `c`,
  `G` and Planck units.
