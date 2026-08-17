# Result: the strong tangent modes are not the geometric lapse subspace

Date: 2026-08-17

## Provenance

```text
prior-art gate                                  24d2ce6
confirmatory protocol                           e50a0ea
registered implementation                       614d083
preserved first OPEN run (12/14)                2eb8b40
valid-OPEN implementation clarification          8c5dae2
corrected implementation                         8d964b8
passing artifact                                 816473e
```

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py`.

Artifact:

```text
reproducible/gravity_600cell_dust_hyperbolic_lapse_alignment.json
SHA-256 a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff
```

Only this targeted verifier and its direct 43-control geometry import were
run.  The full suite was not run.  Final result:

```text
14/14 PASS
HYPERBOLIC_EXTREME_SUBSPACE_OPEN
```

The first run's nonzero exit was procedural: it correctly assigned the valid
OPEN outcome but counted the failed physical gap predicates as verifier
failures.  The correction changed only pass semantics.  It did not change a
gap, subspace, angle, uncertainty, label or outcome.

## What was compared

For both schedule parities and all seven minimal `2T` sectors, the verifier
reconstructed

```text
Y = J^-1 R : boundary phase -> (internal,new) slab variations
```

under all four high-precision derivative estimates.  It then fixed, before
reading eigenvectors:

- `C_weak`: the `5d`-dimensional canonical Schur lift of the five pole-edge
  coordinates;
- `G`: the independently derived `5d` geometric vertex-lapse space;
- `E_plus`: the `5d` most expanding tangent invariant directions;
- `E_minus`: the `5d` most contracting directions.

The spectral spaces were transported through `Y` and compared separately
with `C_weak` and `G`.  Combining candidates or choosing the better branch
sector by sector was forbidden.

All Flint determinants and response ranks passed.  The old
canonical-versus-geometric lapse distances were independently reproduced
with maximum error `2.64e-11`, well inside the frozen `2e-8` control.

## Complete look-elsewhere ledger

There were exactly 56 comparisons:

```text
2 parities x 7 sectors x 2 branches x 2 candidates.
```

The fixed hit fractions were

| comparison | IDENTIFIED | OPEN | SEPARATED |
|---|---:|---:|---:|
| expanding vs canonical weak lift | 0/14 | 10/14 | 4/14 |
| expanding vs geometric lapse | 0/14 | 0/14 | 14/14 |
| contracting vs canonical weak lift | 0/14 | 0/14 | 14/14 |
| contracting vs geometric lapse | 0/14 | 0/14 | 14/14 |

Total: `0 IDENTIFIED`, `10 NUMERICALLY_OPEN`, `46 SEPARATED`.

## Clean negative result

**DERIVED COMPUTATIONAL NEGATIVE, under the frozen logarithmic edge norm:**
neither extreme branch is the geometric vertex-lapse subspace.  All 28
geometric comparisons are resolved `SEPARATED`.

The distances are small in absolute terms but far above their calibrated
uncertainties:

```text
expanding vs geometric       about 9.30e-6 ... 1.05e-5,
contracting vs geometric     about 4.55e-5 ... 5.11e-5,
distance uncertainties       about 1.6e-8 ... 6.4e-8.
```

Thus it would now be false to dismiss the 119 strong reciprocal pairs merely
as geometric lapse directions.  The equality of dimensions
`119+1=120` was a useful clue, but not a subspace theorem.

## What survives of the weak-sector explanation

**PATTERN / OPEN:** the expanding branch is extraordinarily close to the
*canonical algebraic* Schur lift:

```text
nontrivial-sector projector distances  1.64e-6 ... 3.30e-6,
angles                                  9.37e-5 ... 1.89e-4 degrees.
```

Ten comparisons are numerically open and four are resolved separated; none
is identified.  By contrast, the contracting branch lies farther from the
canonical lift, at roughly `5.4e-5`--`5.9e-5`, and is separated in every
sector.

This asymmetric pattern supports a limited statement: the weak Schur block
probably controls the large one-step amplification through the inverse
Legendre solve.  It does **not** show that the amplified modes are gauge or
that they can be quotient out.  The canonical weak lift and geometric lapse
space were already distinct at about `1.2e-5`, so replacing one by the other
would repeat the very ambiguity this test was designed to expose.

## Why the combined outcome is OPEN

All nontrivial sectors have a clean extreme/non-extreme modulus gap:

```text
minimum nontrivial gap >= 8.51.
```

The trivial homogeneous sector does not.  Its fifth candidate pair is the
near-unit pair

```text
0.99390373, 1.00613367,
```

so the extreme boundary gap is only `1.006`, below the preregistered `>2`
gate.  Consequently the complete `120/120` identification is numerically
open even before considering principal angles.  The 119 strong modes remain
well selected; the homogeneous 120th pair does not.

## Physical judgement

- **DERIVED:** the full one-step canonical map exists and contains 119 strong
  expanding plus 119 reciprocal contracting directions.
- **DERIVED NEGATIVE:** those strong invariant spaces are not exactly the
  frozen geometric vertex-lapse space.
- **PATTERN:** the expanding space is close to the canonical weak Schur lift,
  indicating algebraic amplification by the almost-singular Legendre block.
- **OPEN:** whether the modes change gauge-invariant Regge curvature.
- **OPEN:** whether they are constraint-violation modes, physical
  scalar/vector/tensor perturbations, or mixtures.
- **OPEN:** multi-tick growth; a one-step eigenvalue about 47 is not yet a
  Lyapunov exponent on an evolving background.

Therefore the large eigenvalues must not be advertised as gravitational-wave
instabilities, but they also cannot be removed as lapse artifacts.  The next
logical gate is no longer another length-space projection: it is the derived
linear response of four-dimensional deficit angles (and then intrinsic
boundary curvature) on these same invariant subspaces.  Exact or suppressed
curvature response would support a pseudo-gauge interpretation; a resolved
large curvature response would make the modes physically load-bearing.

## Post-result prior-art reconciliation

The post-result search strengthens that choice.  Hoehn's
[canonical linearized Regge analysis](https://arxiv.org/abs/1411.5672)
identifies lattice gravitons as propagating **curvature** degrees of freedom,
not by edge-vector overlap alone.  Dittrich and Hoehn's
[canonical framework](https://arxiv.org/abs/0912.1817) explains
background-dependent pseudo-constraints but does not equate every weak
Legendre direction with gauge.  In continuum perturbation theory,
[curvature-based gauge-invariant methods](https://arxiv.org/abs/gr-qc/9801071)
likewise separate physical modes from constraint-violation subsystems through
curvature observables.

No located primary source supplies this exact 600-cell comparison.  External
novelty remains **OPEN**; the calculation is a project-specific diagnostic,
not a new general theorem.
