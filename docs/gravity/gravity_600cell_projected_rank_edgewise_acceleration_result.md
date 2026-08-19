# Result: the canonical carrier passes the homogeneous Regge calibration

Date: 2026-08-19

## Headline

For the complete direct cellular Lorentzian Regge action with one homogeneous
scale, one global lapse and one conserved total dust mass, the canonical
projected carriers give

```text
carrier                              a
P(sd K_600)                    -0.505139799211
P(Esd_2(sd K_600))             -0.501337041736
closed-dust FLRW target        -0.500000000000
```

where

```text
log(R_1/R_0)=a*(tau/R_0)^2+O((tau/R_0)^4).
```

The absolute error decreases from `0.00513979921` to `0.00133704174`, a factor
of `3.84416`, in the first canonical refinement.  The coefficients were
committed before the FLRW target was loaded.

The result is **DERIVED NUMERICAL CONTROL, ADVERSARIALLY CORROBORATED**.  It
shows that the non-arbitrary carrier passes the known homogeneous continuum
calibration.  It is not a new law of gravity and does not yet contain local
gravitational degrees of freedom.

## Provenance ledger

| stage | commit |
|---|---|
| canonical-carrier consolidated result | `8248691` |
| acceleration prior-art gate | `1921519` |
| blind primary protocol | `ebe3889` |
| registered blind verifier | `b8c319f` |
| frozen blind artifact | `9469e33` |
| registered target comparison | `ac92931` |
| frozen target comparison | `f65dfc8` |
| adversarial protocol | `3b7bd6c` |
| registered first audit | `f6a74f8` |
| preserved first audit failure | `391f758` |
| preregistered audit correction | `1c3318d` |
| corrected audit implementation | `914df3f` |
| frozen corroborating artifact | `13c0856` |

The blind artifact is

```text
reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json
SHA-256 2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2.
```

It explicitly records

```text
continuum_target_loaded = false
projected_red_coefficients_loaded = false.
```

The corrected adversarial artifact is

```text
reproducible/gravity_600cell_projected_rank_edgewise_acceleration_adversarial_corrected.json
SHA-256 9b995e7f331697a00cc939869773b75eb039ac3cc6c921e20b6dbeb2fed30409.
```

Two corrected audit executions produced this identical digest.  This second
execution checks reproducibility; the mechanical independence is supplied by
the different orbit-compressed construction, not by rerunning one program.

## 1. Primary blind calculation

The two carriers were independently reconstructed as

```text
P(sd K_600):
  f=(2640,17040,28800,14400)
  V_bar=19.147932918312847

P(Esd_2(sd K_600)):
  f=(19680,134880,230400,115200)
  V_bar=19.583480465413963.
```

Every face has incidence two, both Euler characteristics vanish, all
tetrahedral volumes are positive and all sampled frustum metrics have
Lorentzian inertia `(3,1)`.

The direct action includes:

- every timelike edge-trapezoid hinge;
- both complete spatial-boundary terms;
- no averaged Schlaefli symbol or averaged angle;
- one total dust term `-8*pi*M*sqrt(rho)`.

At each level the mass was selected before the dynamic coefficient by the
static global lapse constraint at unit volume radius.  Seam and lapse
equations reconstructed the same coefficient within their preregistered
tolerances.  The blind verifier passed `10/10`.

Only after commit `9469e33` did the comparison verifier load `a_FLRW=-1/2`.
It passed `5/5` and assigned

```text
CANONICAL_CARRIER_HOMOGENEOUS_CALIBRATION_PASS.
```

## 2. Adversarial reconstruction

The independent audit did not evaluate the final coefficient by the primary
full-array derivative route.  It:

- found the 600 source tetrahedra with NetworkX maximal cliques;
- rebuilt the rank-selected direct eight-child split;
- classified unique edges, faces and tetrahedra by intrinsic squared-distance
  signatures up to vertex permutation;
- evaluated one representative per intrinsic class with its independently
  counted multiplicity;
- used real five-point logarithmic derivatives with different `eta` values,
  steps and coefficient sentinels.

At 13-digit signatures the fine carrier compresses from 115,200 tetrahedra to
168 numerical intrinsic classes.  Its three held-out action errors are at
most `1.15e-10`; the base errors are at most `5.16e-12`.

The deliberately wrong control adds one copy to the largest local
tetrahedron class while leaving the global census fixed.  It changes the
held-out action by about `2.02e-2` relative, so the audit is capable of
detecting an incorrect multiplicity.

The audit coefficients are

```text
P(sd K_600)                 -0.505137095442
P(Esd_2(sd K_600))          -0.501336162358.
```

Their differences from the full primary route are respectively
`2.70e-6` and `8.79e-7`, within the frozen `2e-5` gate.  The already certified
regular-600-cell positive control is reproduced within `1.44e-5`.  The
corrected audit passed `11/11` and assigned

```text
ADVERSARIAL_CANONICAL_ACCELERATION_CORROBORATED.
```

## 3. The preserved failed audit

The first audit result is not erased.  It is preserved at commit `391f758`
with artifact digest

```text
e41a7f57be6995cb39b8fdb89fb981263473be9bbf8982ae0bb84f65d07fc8f6.
```

It failed `2/11` gates for two distinct numerical reasons.

First, the static *total* action vanishes by the selected dust cancellation,
so its relative error used an approximately zero denominator and returned an
artificial value near one.  The correction compares the nonzero static
gravitational action.

Second, 11-digit intrinsic signatures produced held-out errors up to
`2.78e-8`, slightly above the `2e-8` gate.  A target-free precision census
showed that 13 digits reduce the error to approximately `1e-10`.  The
correction selected 13 digits and tightened the within-class residual gate;
it did not loosen the acceptance threshold or change any coefficient grid.

## 4. Interpretation

- **DERIVED NUMERICAL CONTROL:** the direct homogeneous action is well defined
  on both finite canonical carriers and the first refinement moves strictly
  toward closed-dust FLRW.
- **ADVERSARIALLY CORROBORATED:** a weighted intrinsic-shape route with a
  different derivative stencil reproduces both coefficients.
- **PATTERN:** the one-step error-reduction factor `3.844` is compatible with
  a second-order trend.
- **OPEN:** an asymptotic order or infinite-refinement theorem.  One step
  cannot establish either.
- **STRUCTURAL:** the declared radial projection and round `S3` background.
- **OPEN:** local dust, local lapse constraints, the physical quadratic
  fluctuation spectrum and propagation.

The result removes one genuine ambiguity from the earlier positive
calibration: it no longer depends on choosing one of the `3^600`
projected-red central diagonals.  But the homogeneous ansatz itself supplies
only the global scale mode.  Agreement with Friedmann is therefore a required
continuum control, not evidence for new gravitational degrees of freedom.

No claim is established about gravitons, a limiting speed, a physical tick,
inertia, Planck units, inflation or particle masses.

## 5. Literature reconciliation after the result

The broad result is firmly within known Regge cosmology.  Barrett, Galassi,
Miller, Sorkin, Tuckey and Williams already illustrated a local implicit
evolution scheme with homogeneous dust on a 600-cell Friedmann cosmology:

- *A Parallelizable Implicit Evolution Scheme for Regge Calculus*,
  arXiv:`gr-qc/9411008`.

Tsuda and Fujiwara project subdivisions of the 600-cell to form geodesic
4-domes, but replace their cumbersome direct irregular dynamics by a
pseudo-regular angularly averaged model:

- arXiv:`2011.04120`, DOI `10.1093/ptep/ptab079`, Section 6.

Edgewise subdivisions are independently used as regular refinement sequences
in a Regge/BF setting:

- M. Kisielowski, arXiv:`1704.00998`, DOI
  `10.1007/s00023-018-0747-6`, refinement discussion around Figure 4.

Thus neither 600-cell Friedmann dynamics, projected geodesic domes nor
edgewise refinement is new.  The search did not locate this exact complete
direct irregular coefficient on the projected rank-edgewise carrier, but
absence from a search is not a novelty proof.  External novelty remains
**OPEN**.

The literature also warns that repeated pure barycentric subdivision can
produce increasingly poor simplex shapes.  The present object is one initial
barycentric chamber complex followed by rank-edgewise refinement and radial
projection, not repeated barycentric subdivision.  Its first shape gate
passes, but an all-level shape theorem remains **OPEN**.

## 6. Next discriminating question

The next mission must not compute another homogeneous number.  It must first
select a local conserved dust discretization and local lapse equations on the
canonical carrier, then form the constraint-reduced quadratic Regge action
for boundary-edge perturbations.

Only a stable low-mode generalized spectrum under refinement could test the
proposed discrete wave equation and an effective propagation speed.  Without
that constraint and matter step, a large Hessian would mix physical,
pseudo-gauge and discretization modes and would not answer a well-posed
question.

Only the three mission-specific verifiers and static registry/documentation
checks were run.  No full-suite claim is made.
