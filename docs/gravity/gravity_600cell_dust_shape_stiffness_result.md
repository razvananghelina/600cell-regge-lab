# Result: the autonomous shape carrier has resolved negative stiffness

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL NEGATIVE.**  The complete action-relative
`600`-position shape carrier is dynamically autonomous and has a resolved
definite kinetic form, but its Hermitian generalized stiffness is not
nonnegative.  On one complete carrier the blind inertia census is

```text
250 positive resolved
 30 negative resolved
108 zero-consistent
212 open
-------------------
600 full-multiplicity shape coordinates per schedule/variant decomposition.
```

The `30` certified negative directions occur only in two one-dimensional
binary-tetrahedral sectors, with inertia `15 negative + 10 positive` in each.
The same pattern occurs in both staircase schedules and all four derivative
variants.

Therefore the shortcut

```text
autonomous shape carrier  =>  positive stable wave carrier
```

is **REFUTED**.  This is not yet a ghost theorem and does not kill Regge
gravity: no pseudo-constraint/gauge quotient has determined whether the
negative directions survive as physical perturbations.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art gate | `e37d80c` |
| original blind protocol | `4c29645` |
| preregistered sign-consistency clarification | `c719a87` |
| registered verifier before first execution | `6efa0b9` |
| pre-spectrum missing-import repair | `f22ce7b` |
| disclosed post-result count correction | `f3babf4` |
| post-result audit fields added | `936df3d` |
| first representative-count artifact | `443b9b0` |
| disclosed irrep-multiplicity correction | `0b112c3` |
| multiplicity-aware verifier | `f0ea0eb` |
| final certified artifact | `d34260d` |

Targeted verifier:

```text
reproducible/verify_gravity_600cell_dust_shape_stiffness.py
```

Certified artifact:

```text
reproducible/gravity_600cell_dust_shape_stiffness.json
SHA-256 03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868
```

Two complete scientific executions returned byte-identical artifacts and

```text
12/12 PASS
SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED.
```

The first attempted execution stopped before reconstructing any spectral
carrier because the new wrapper omitted its `scipy.sparse` import.  The repair
was committed without changing the protocol, and no scientific number had
been seen at that point.

Only the mission-specific verifier was run.  The full suite was not run.

## Disclosed count and representation-multiplicity corrections

The original protocol wrote

```text
2 schedules * 4 variants * sum(25 d),  sum(d)=12,
```

and initially displayed `4,800` per object.  Commit `f3babf4` first corrected
the direct evaluation to

```text
2 * 4 * 25 * sum(d) = 2,400 directly diagonalized values per object.
```

That correction was incomplete as a count of the full carrier.  Each minimal
irrep block of dimension `d` occurs with representation multiplicity `d` in
the regular carrier.  Therefore

```text
representative block dimension = 25 sum(d)   = 300,
full shape-carrier dimension   = 25 sum(d^2) = 600,

direct values per object across eight audits = 2,400,
full-multiplicity values per object           = 4,800.
```

Commit `0b112c3` discloses this second post-result correction.  The verifier
now records both ledgers explicitly.  The certified negative sectors both
have `d=1`, so the `30` negative directions per full carrier, thresholds,
outcome branch and scientific conclusion are unchanged.

## Complete hypotheses

The result is conditional on all of the following:

1. the fixed regular 600-cell and its literal `720` logarithmic
   signed-squared spatial edge coordinates;
2. the first two accepted nonstationary fixed-mass dust-Regge slabs;
3. the adjacent-slice edge identification used by the centered Jacobi
   recurrence;
4. all seven frozen minimal binary-tetrahedral sectors, both schedules and all
   four derivative variants;
5. the canonical vertex-conformal incidence map
   `(C sigma)_uv=sigma_u+sigma_v`;
6. the action-selected shape complement

   ```text
   S_H = ker(C* H_M),  H_M=(M+M*)/2;
   ```

7. the committed component midpoint/radius enclosures and the frozen
   `10/100` classification bands;
8. no independent dust perturbations, constraint quotient, proper-time unit,
   continuum refinement or tensor-harmonic target.

It is not a theorem for later ticks, nonlinear perturbations, arbitrary Regge
triangulations or a continuum limit.

## Exact finite objects

For the orthonormal shape basis `W`, the verifier formed

```text
M_S     = W* [(M+M*)/2] W,
V_S     = W* [(V+V*)/2] W,
B       = -M_S,
A       = -V_S,
Omega_S = W* (M^-1 V) W.
```

The Hermitian definite pencil

```text
A x = lambda B x
```

has the same generalized eigenvalues as `V_S x=lambda M_S x`.  The overall
action sign cancels, so the negative generalized signs are not a sign-convention
artifact.

Because `B` is positive definite, Sylvester inertia gives

```text
inertia(B^-1/2 A B^-1/2) = inertia(A).
```

Thus signs were classified from the better-conditioned Hermitian form `A`,
not by forcing small generalized eigenvalues away from zero.

## Kinetic and carrier controls

All `56/56` shape kinetic restrictions are
`POSITIVE_DEFINITE_RESOLVED`.  Their minimum-eigenvalue margins span

```text
5.86e5 ... 1.30e7 error units.
```

All conformal and shape carrier ranks are resolved, the direct sums remain
regular, the geometry retains exact `120/720/1200` vertex/edge/triangle
counts, and the seven high-precision symmetry bases have residual below
`1.55e-98`.

Changing the overall action sign swaps the raw sign of both restricted forms
but neither the generalized pencil nor the conclusion.

## Complete sign census

Across schedules and derivative variants, the Hermitian-pencil counts are

| label | 2,400 direct representatives | 4,800 full multiplicity | one full carrier |
|---|---:|---:|---:|
| positive resolved | 1,008 | 2,000 | 250 |
| negative resolved | 240 | 240 | 30 |
| zero-consistent | 288 | 864 | 108 |
| open | 864 | 1,696 | 212 |

The negative-form eigenvalues have margins

```text
135.87 ... 264.33 error units,
```

strictly beyond the preregistered `100`-unit resolved threshold.  The weakest
negative midpoint is about `-1.9583e-4` against a complete restricted-form
error about `1.4413e-6`.

Every certified negative occurs in sector `4` or `5`, both of irrep dimension
one:

```text
sector 4: 15 negative, 10 positive,
sector 5: 15 negative, 10 positive.
```

All other possible negative-looking entries remain honestly `OPEN` or
`ZERO_CONSISTENT`.  Consequently `30` is a certified lower count on this
finite carrier, not a claim that the remaining `212` open directions are
nonnegative.

## Actual recurrence block versus Hermitian pencil

The normalized shape block has real midpoints in

```text
-7.3153e-6 ... 1.44721e-4.
```

All `2,400/2,400` directly diagonalized eigenvalues, equivalently all
`4,800/4,800` full-multiplicity entries, are `REAL_CONSISTENT` under their
conservative Bauer--Fike envelopes.  Their full-multiplicity sign census is
weaker:

| label | count |
|---|---:|
| positive resolved | 244 |
| zero-consistent | 3,404 |
| open | 1,152 |
| negative resolved | 0 |

This does not contradict the Hermitian inertia.  Normalizing by `M` amplifies
the stored solve uncertainty, and the actual negative midpoints do not cross
the `100`-error sign threshold.  Inertia of `A` supplies the stronger sign
certificate while definiteness of `B` protects that sign.

The independent action-compatibility residual

```text
M_S Omega_S - V_S
```

is `ZERO_CONSISTENT` in all `56/56` cells, at only
`0.00347 ... 0.11993` error units.  Hence the two constructions are compatible
within the frozen enclosure; the result is not produced by silently replacing
the recurrence with an unrelated symmetric matrix.

## Schedule and implementation robustness

All `56/56` ordered even/odd comparisons are `SCHEDULE_ROBUST`.  Their largest
distance/error ratio is `3.11e-9`.  Operational/shadow and
validation-primary/validation-shadow variants give the same inertia pattern.

The artifact was regenerated twice after exposing all inertia margins and was
byte-identical.  Registry inspection gives `276` entries, `276` unique names,
no missing verifier and no unregistered verifier outside the two declared
exclusions.

## What the negative means physically

- **DERIVED COMPUTATIONAL:** the action-selected `600`-position shape carrier
  is autonomous for the finite centered recurrence.
- **DERIVED COMPUTATIONAL:** its kinetic restriction is definite and far from
  singular.
- **DERIVED COMPUTATIONAL NEGATIVE:** its stiffness has at least `30`
  resolved negative generalized directions.
- **STRUCTURAL:** these are finite dust-Regge shape directions on a changing
  background.
- **OPEN:** whether the two one-dimensional symmetry sectors are scalar,
  longitudinal, dust-coupled, pseudo-constraint, or physical curvature
  modes.
- **OPEN:** whether the full recurrence including `Gamma` grows or remains
  bounded on these directions.

A negative stiffness in a stationary undamped oscillator would indicate an
exponential mode.  This background is neither stationary nor reduced: it has
a resolved first-difference term and unresolved gravitational constraints.
Calling the `30` directions ghosts or instabilities now would overstate the
calculation.

Conversely, ignoring them and fitting only the positive entries to an `S^3`
Laplacian would be selection bias.  The entire carrier must first be reduced
by an action-derived constraint test.

## Post-result prior-art check

The second search used the observed terms `mixed shape stiffness inertia`,
`negative Regge FLRW Hessian modes`, `one-dimensional symmetry sectors` and
`curved-background pseudo-constraints`.

- Bahr--Dittrich's warning remains directly relevant: curvature breaks exact
  discrete gauge symmetry and produces pseudo-constraints:
  <https://arxiv.org/abs/0905.1670>.
- Hoehn's lattice-graviton count requires a gauge/curvature separation before
  physical interpretation: <https://arxiv.org/abs/1411.5672>.
- Rostworowski's continuum FLRW reduction separates two gravitational master
  scalars from a matter transport scalar only after solving constraints:
  <https://arxiv.org/abs/1902.05090>.
- De Felice--Fabri find nontrivial limitations in generalized 600-cell
  evolution, but do not publish this shape-pencil inertia:
  <https://arxiv.org/abs/gr-qc/0106077>.
- Liu--Williams discuss local/global variation and finite-resolution failure
  in closed Regge FLRW models, not this spectrum:
  <https://arxiv.org/abs/1501.07614>.

No located primary source identifies the exact two-sector `15+15` negative
carrier found here.  Search failure is not proof of novelty; external novelty
remains **OPEN**.

## Next load-bearing gate

Do not compare the `250` positive entries with a desired spatial spectrum yet.
An inventory performed after this result gives an important framing
correction.  The already certified complete pre-Legendre Jacobian has

```text
rank 1560/1560, error-consistent nullity 0.
```

Therefore this fixed finite system has no exact pre/post-constraint kernel
from which to construct a canonical quotient.  Its `120` weak directions are
resolved nonzero pseudo-constraint candidates.  Quotienting them by a chosen
singular-value threshold would be fitted and is forbidden.  A true reduced
phase space could reappear only after an independently derived exact symmetry,
matter extension or refinement limit.

The next calculation on the fixed carrier should instead ask:

1. are the two `15`-dimensional negative-stiffness eigenspaces invariant under
   both restricted `Gamma` and `Omega`, or do they mix with the other ten
   directions in their symmetry sectors;
2. what are the eigenvalues and singular amplification of the complete
   doubled recurrence in both full `25`-dimensional sectors, including
   `Gamma`, not merely zero-frequency stiffness;
3. does the negative carrier intersect the expanding/contracting invariant
   spaces selected independently by that recurrence;
4. only after this, derive a spatial curvature/tensor intertwiner and compare
   dispersion under refinement.

Resolved growing multipliers would establish local hyperbolicity of these
finite sectors, not automatically a continuum ghost.  Absence of growth would
show that the drift/mixing term cures the naive stiffness reading.  Either
outcome is meaningful and requires no fitted quotient.
