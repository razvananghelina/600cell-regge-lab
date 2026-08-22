# Protocol: bounded no-go for an exact finite-height constraint quotient

Date: 2026-08-22

Status: preregistered before the verifier source, registry entry, result
artifact, or consolidated verdict exists.

Prior-art gate commit: `cbbf7b6`.

## Frozen question

Under every hypothesis in
`gravity_600cell_finite_height_constraint_quotient_prior_art.md`, do the first
and second accepted finite-height slabs admit a nontrivial local boundary
phase-space quotient selected by exact pre/post constraints of the frozen
Regge-plus-conserved-dust action?

The claim is restricted to Legendre-degeneracy constraints and exact gauge
directions of the accepted local canonical relation.  It is not a claim about
future propagated constraints, arbitrary symmetry reduction, refinement, a
different action, or an enlarged matter carrier.

## Frozen inputs and hashes

The verifier shall read only these scientific input artifacts:

| input | SHA-256 |
|---|---|
| first full-boundary adversarial tangent JSON | `ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7` |
| second full-boundary adversarial tangent JSON | `1355f8cf339d18c1cf2855ecb1228e97e868d73f7a1ef739e4c11ce9521fcd4b` |
| internal-kernel canonical reconciliation JSON | `81ec0379247023451e82ab42f5beb026ee2d1b083aa5e2553e42b894554266f6` |
| old full anisotropic Legendre-rank control JSON | `7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226` |
| prior-art gate note | `53de3f6262df6e3a0cdff916d11d435fc462289b5db080bb845dcb425bb270c5` |

No tangent eigenvalue, desired mode count, continuum dispersion law, limiting
speed, `G`, Planck scale or particle target may be parsed.

## Two proof routes

The verifier must implement two logically distinct finite-dimensional routes.

### Route A: rank and image

For every accepted first/second slab realization:

1. require strict certified regularity of the complete pre-Legendre matrix,
   including `sigma_min > frozen gap gate` at every stored Richardson level;
2. require the complete boundary tangent to carry its accepted canonical
   label;
3. require all four two-step products to carry their accepted canonical
   labels;
4. infer, by the implicit-function theorem, a local update defined on an open
   subset of the complete `1440`-dimensional incoming phase space;
5. use the invertibility of a symplectic tangent to infer an open image in the
   complete outgoing phase space.

The resulting local pre- and post-constraint codimensions are both zero.

### Route B: constraint-covector contradiction

Suppose an irreducible local post-constraint `C` vanished on the image of a
full-rank tangent `T`.  Differentiating `C(T(z))=0` gives

```text
dC * DT = 0.
```

Since `DT` is invertible, `dC=0`, contradicting irreducibility.  The same
argument applied to the inverse tangent excludes an irreducible local
pre-constraint.  The verifier shall reproduce this implication with exact
rational linear algebra on deterministic matrices, independently of the
rank/image bookkeeping in Route A.

Agreement of the two routes is required.  A mismatch leaves the verdict
`OPEN`.

## Frozen controls

### Positive regular control

Use an exact rational quadratic generating function with one old, one
internal and one new coordinate.  Its combined pre-Legendre matrix must be
invertible, its exact canonical tangent must be symplectic, and both proof
routes must return constraint codimension zero.

### Singular negative control

Change the exact mixed coefficients so that the combined pre-Legendre matrix
is singular.  The exact pre/post image must then have positive codimension,
and the no-go classifier must refuse to fire.

### Small-nonzero hostile control

Use an exact rational mixed coefficient `10^-100`.  Exact arithmetic must keep
the system regular even though a deliberately hostile `10^-12` numerical
cutoff calls it zero.  The threshold-derived quotient must be explicitly
rejected.

### Future-singular control

Compose a regular first move with a deliberately singular later move.  The
later pre-constraint must pull back to a nontrivial condition on the earlier
phase data.  This control must pass so that the verifier cannot overstate two
regular slabs as an infinite-history theorem.

### Repository weak-mode control

The old complete anisotropic artifact must retain `1560/1560` regularity and
its pseudo-constraint interpretation must remain non-gauge/open.  It is a
control only; it shall not be used as a premise for the finite-height rank.

## Frozen checks

The registered verifier shall report exactly these checks:

1. all five frozen input hashes agree;
2. the verifier is registered exactly once and the registry contains no
   duplicates;
3. the first adversarial artifact retains its accepted `22/22` outcome;
4. both first-slab parities are strictly regular and canonical above their
   frozen rank gates;
5. the second adversarial artifact retains its accepted `28/28` outcome and
   its target firewall;
6. all six first/second stage realizations are strictly regular and canonical
   above their frozen rank gates;
7. all four two-step products are canonical and schedule robust;
8. the internal lapse tangent is still removed by fixed incoming momentum;
9. the exact regular positive control passes both proof routes;
10. the exact singular negative control is detected;
11. the exact `10^-100` hostile control refutes threshold quotienting;
12. the future-singular control propagates a constraint backward and bounds
    the claim;
13. the old weak-mode artifact remains a regular non-gauge control;
14. Routes A and B agree on zero local constraint codimension for every
    accepted finite-height map and on the bounded outcome.

The result artifact shall contain all tested margins, exact control matrices,
ranks, codimensions, firewall flags, scope exclusions and reopening
conditions.  It must not contain a physical mode spectrum.

## Frozen outcomes

### Acceptance / bounded no-go

Emit

```text
FINITE_HEIGHT_TWO_STEP_EXACT_CONSTRAINT_QUOTIENT_BOUNDED_NO_GO
```

only if all 14 checks pass.  The licensed interpretation is:

> On the frozen two-slab branch-B region, the action supplies a regular local
> symplectic map on the complete boundary phase carrier.  It supplies no
> positive-codimension pre/post-constraint surface and therefore no
> nontrivial exact local gauge quotient selected by Legendre degeneracy.

### Positive constraint signal

If an accepted slab has a certified singular pre-Legendre matrix or a
nontrivial exact covector annihilating the map image, emit

```text
FINITE_HEIGHT_EXACT_CONSTRAINT_SIGNAL
```

and keep the quotient route open.  This would be the more interesting result.

### Otherwise

Any failed provenance, control, rank margin, canonicality, route agreement or
scope firewall yields

```text
FINITE_HEIGHT_CONSTRAINT_QUOTIENT_OPEN
```

with no no-go claim.

## Reopening conditions

The bounded no-go, if obtained, may be reopened only by at least one new
derived premise:

1. a later accepted slab with an exact singular Legendre transform whose
   constraint propagates to the present data;
2. an exact continuous gauge or momentum-map symmetry of an enlarged action;
3. a declared refinement family in which a geometrically frozen weak carrier
   converges to an exact kernel under a preregistered scaling law;
4. a different/perfect action with exact discrete diffeomorphism symmetry.

Weak but nonzero singular directions at one resolution are not a reopening
condition.

## Execution policy

Register the verifier once in `reproducible/run_all.py`, execute only this
targeted verifier, preserve any first failure, and do not run the full suite.
Because the scientific input ranks already have mechanically different dense
replications, the new adversarial requirement is the independent covector
proof plus the singular, tiny-nonzero and future-singular controls above.
