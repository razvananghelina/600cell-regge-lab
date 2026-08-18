# Result: complete two-step canonical tangent cocycle

Date: 2026-08-18

## Provenance

```text
prior-art gate                         d51fbbf
blind protocol                         4ec04b1
registered implementation              950ac48
preserved first-run control failure    b3cc397
protocol correction                    03cb878
corrected verifier                     97d5a6f
passing blind artifact                 25d64a4
```

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_two_step_full_tangent.py`.

Artifacts:

```text
reproducible/gravity_600cell_dust_two_step_full_tangent.json
SHA-256 f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc

reproducible/gravity_600cell_dust_two_step_full_tangent.npz
SHA-256 ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d
```

No target spectrum, continuum harmonic degeneracy or desired count was used
to build either tangent.  Only the targeted verifier was run; the full suite
was not run.  The corrected targeted run was repeated in `319.76 s` and
returned `16/16 PASS` with the same JSON and NPZ SHA-256 values.

## The transparent first-run failure

The first execution stopped at `15/16`, after all mathematical controls had
passed, because the protocol incorrectly counted the arrays in an archive
containing two maps.  The correct cardinality is

```text
2 schedules * 7 sectors * 4 variants * 2 maps * 4 fields = 448,
```

not `224`.  The failed JSON and the numeric archive were preserved before the
correction.  The protocol and verifier correction changed only this expected
cardinality; no operator, threshold, state or spectrum changed.  In
particular, the failed and passing runs have the same NPZ SHA-256 printed
above.

## Exact object and background controls

The calculation reconstructs the complete action-generated maps

```text
T_1 : T Phase(boundary 0) -> T Phase(boundary 1),
T_2 : T Phase(boundary 1) -> T Phase(boundary 2),
C_21 = T_2 T_1.
```

Each phase carrier has `720+720=1440` real position/momentum dimensions.  The
two slabs are the already accepted, canonically sewn, fixed-total-mass dust
solutions.  The maximum complete thirty-component seam residual is

```text
6.0468e-44,       inherited bound 3.6514e-21.
```

The second slab uses its committed unequal old/new scales and lapse, with the
same conserved dust mass.  No mass is recomputed from the later scale.

For both independently derived order-24 schedules:

- the full 600-cell boundary carrier and literal boundary identification
  reconstruct;
- all Lorentzian branch and reality controls pass;
- all `28` second-slab pre-Legendre determinant balls exclude zero;
- all seven blocks of `T_2` are canonical within propagated ball errors;
- all seven directly multiplied blocks of `C_21` are canonical within
  propagated ball errors;
- all fourteen singular-spectrum and all fourteen eigen-spectrum schedule
  comparisons receive `SCHEDULE_ROBUST`.

**DERIVED COMPUTATIONAL:** this certifies a finite, regular, symplectic,
non-autonomous two-step tangent cocycle on the complete declared carrier.

## Blind census

The first and second one-step maps have the same blind branch count:

| map | contracting | expanding | resolved off-unit | unit-consistent | open |
|---|---:|---:|---:|---:|---:|
| `T_2`, even | 119 | 119 | 238 | 1200 | 2 |
| `T_2`, odd  | 119 | 119 | 238 | 1200 | 2 |

For the two-step product the previously unresolved pair is resolved:

| map | contracting | expanding | resolved off-unit | unit-consistent | open |
|---|---:|---:|---:|---:|---:|
| `C_21`, even | 120 | 120 | 240 | 1120 | 80 |
| `C_21`, odd  | 120 | 120 | 240 | 1118 | 82 |

The `80` versus `82` difference is only a conservative uncertainty-label
difference for near-unit, highly non-normal eigenvalues.  It is not a
schedule-dependent spectrum: all preregistered schedule comparisons are
robust.

The resolved product moduli range from approximately

```text
0.00510794326067 to 195.773513715.
```

The count is restored with the correct representation multiplicities; it is
not a count of minimal-block eigenvalues.

**DERIVED COMPUTATIONAL:** the two-step map has exactly `120` resolved
contracting and `120` reciprocal expanding eigenvalues under the frozen
classification.  The one-step `119+1` split does not persist.

## What happened to the homogeneous pair

A post-result restriction of the committed trivial-sector midpoint to the
literal uniform `(q,p)` plane gives

```text
T_2:   -1.0123048012,       -0.9878447665,
C_21: -13.9284502399,       -0.07179549647.
```

The two product values are reciprocal within the certified symplectic error.
The uniform-plane leakage of the binary midpoint is below `9.6e-11` in both
schedules.  This restriction was not part of the blind verdict and is
therefore labelled **POST-RESULT DIAGNOSTIC**, not a new preregistered hit.
It identifies which pair accounts for the count change, but no physical mode
name is assigned to it.

The result sharply refutes the earlier interpretation of the unresolved
one-step pair as an exact neutral time/lapse/gauge direction.  It was a close
reciprocal pair of a non-stationary one-step map, not a protected line.

## A two-slab curvature-observability consequence

The preceding certified response calculation established

```text
rank F_1 = rank F_2 = 1439,
ker F_1 = K_1,   ker F_2 = K_2,
T_1 K_1 != K_2.
```

Define the two-slab observation map

```text
O(x) = (F_1 x, F_2 T_1 x).
```

Then

```text
ker O = K_1 intersect T_1^{-1} K_2 = {0}.
```

The last equality follows because both kernels are one-dimensional and their
transport separation was resolved by `3.42e13` uncertainty units.

**DERIVED COMPUTATIONAL:** the stacked two-slab internal-curvature response
is injective on the full `1440`-dimensional initial phase tangent.  Every
nonzero initial tangent changes the chosen Regge internal-curvature data in
at least one of the first two slabs.

This is not a count of physical gravitons.  It says only that the particular
instantaneous curvature-blind line is not an exact multi-step gauge orbit.
Constraint covectors, matter perturbations and a physical quotient have not
yet been constructed.

## Interpretation firewall

- **DERIVED COMPUTATIONAL:** `T_2` and `C_21` exist, are regular and
  symplectic on the full finite carrier.
- **DERIVED COMPUTATIONAL:** the blind two-step branch census is `120+120`,
  not `119+1`.
- **DERIVED NEGATIVE:** the missing one-step `+1` is not a dynamically
  protected neutral time/lapse/gauge direction.
- **DERIVED COMPUTATIONAL:** two successive curvature observations have no
  common blind tangent direction.
- **STRUCTURAL:** singular values refer to the declared Euclidean norm on
  logarithmic positions and momenta, not a derived physical energy norm.
- **OPEN:** which resolved directions are physical scalar, vector or tensor
  perturbations, pseudo-constraint artifacts, or frozen-matter artifacts.
- **OPEN:** nonlinear anisotropic propagation, refinement convergence,
  continuum dispersion, a limiting speed and Planck units.
- **OPEN:** external novelty.

In particular, “off the unit circle” is not yet synonymous with “physical
instability”.  The map is non-normal, the background evolves, and the dust
perturbations have not been included as independent canonical variables.

## Post-result prior-art check

The post-result search used the more precise terms `non-autonomous canonical
Regge cocycle`, `curved-background Hessian`, `pseudo-constraints`,
`curvature perturbations` and `600-cell dust stability`.

- [Bahr--Dittrich](https://arxiv.org/abs/0905.1670) show that curved Regge
  solutions generically lack exact discrete gauge symmetries and replace
  constraints by pseudo-constraints.
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) derive the corresponding
  background-dependent canonical structure.
- [Hoehn](https://arxiv.org/abs/1411.5672) identifies propagating lattice
  gravitons only after a gauge and curvature analysis on flat backgrounds.
- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) establish action-generated
  composable canonical simplicial evolution.
- [De Felice--Fabri](https://arxiv.org/abs/gr-qc/0009093) evolve the homogeneous
  dust 600-cell through multiple steps but do not print the present complete
  tangent cocycle.

These sources make the absence of an exact curved-background gauge line
plausible, but they do not prove the finite result above.  No located primary
source contains the present complete `1440`-dimensional two-slab spectrum or
the two-slab curvature-observability statement.  A search cannot establish
novelty, so external novelty remains **OPEN**.

## Correct next physical gate

Do not interpret the `120` pairs by their eigenvalue count alone.  The
one-step geometric-lapse comparison and the internal-curvature response have
already been done: the strong subspace is neither the geometric lapse space
nor curvature-blind.  Repeating another projector comparison would add no
new information.

The next action-selected object is instead the three-slice Jacobi equation

```text
K_- delta q_0 + K_0 delta q_1 + K_+ delta q_2 = 0,
```

obtained by eliminating each slab's internal variables and linearizing the
canonical seam.  It is the position-space recurrence equivalent to the
two-step canonical cocycle and the first legitimate candidate for a discrete
wave operator.  It must first be constructed and certified blindly.  Only
after that commit may its spatial part be compared with an independently
derived 600-cell operator.  A limiting speed would require a stable
temporal/spatial factorization and refinement; it cannot be read from the
present eigenvalue count.
