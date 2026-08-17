# Result: the full 600-cell dust slab generates a canonical 1440-dimensional tangent map

Date: 2026-08-17

## Provenance ledger

```text
prior-art gate                                      5dfc2a7
blind protocol                                      bc114bf
verifier registration                               35567be
preserved first failure (17/19)                     6fbdc08
disclosed SVD-floor correction                      bf88d8e
corrected verifier                                  1574d6e
passing blind artifact                              360253d
```

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py`.

Machine-readable artifacts:

```text
reproducible/gravity_600cell_dust_full_boundary_tangent.json
SHA-256 4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5

reproducible/gravity_600cell_dust_full_boundary_tangent.npz
SHA-256 816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b
```

Only the targeted verifier and its direct 43-control geometry import were
run.  The full suite was not run, by instruction.  Final targeted result:

```text
19/19 PASS
FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED
```

## Exact result

**DERIVED COMPUTATIONAL.**  The accepted non-static Regge--dust slab defines
a regular canonical self-map on the complete boundary phase carrier

```text
(delta q_old, delta p_pre) -> (delta q_new, delta p_post),
dimension 720 + 720 = 1440.
```

This is the first result in this route that does not restrict boundary
perturbations to the 30-orbit invariant quotient.

The result uses seven minimal blocks for the free binary-tetrahedral action,

```text
irrep dimensions d                         3,2,2,2,1,1,1
minimal tangent sizes                    180,120,120,120,60,60,60
real dimensions after multiplicity       540,240,240,240,60,60,60
total                                                   1440.
```

For both schedule parities:

- the literal shift `(u,v) -> (u+120,v+120)` fixes all 720 old/final edge
  identifications; the orbit map is the identity in the frozen ordering and
  the group coordinate is unchanged;
- all four high-precision Hessian estimates have 3,305 nonzero representative
  convolution entries;
- all 28 Flint pre-Legendre determinant balls exclude zero;
- the seven minimal canonical spectra reproduce the prior full-rank artifact
  with maximum normalized errors `1.70e-13` (even) and `8.99e-14` (odd);
- the trivial block reproduces the previously committed `60 x 60` tangent
  with relative error below `5.76e-43`;
- every complete projected Hessian is Hermitian inside its calibrated error;
- all seven tangent blocks satisfy the correct complex-sector identity
  `T* Omega T = Omega` inside calibration.

The symplectic defect norms are only `1e-12` of their already conservative
allowed bounds.  Thus existence and canonicity of the finite tangent map are
not marginal numerical claims.

## The failed run and why it mattered

The first run was preserved at `6fbdc08` with `17/19`.  It passed the primary
symplectic and determinant controls but failed the secondary reciprocal-SVD
test.  The tangent matrices have condition numbers

```text
7.05e10 <= condition_2(T) <= 2.66e11.
```

The original proxy allowed an error of order machine epsilon in a product of
extreme singular values.  A backward-stable SVD only controls each singular
value to absolute scale `eps_machine*sigma_max`, so the extreme product has
floor `eps_machine*sigma_max^2`, equivalently
`eps_machine*condition_2(T)`.  The observed residuals `2.6e-7`--`2.2e-6`
were below that correct bound.

The correction was disclosed and committed before rerunning.  Crucially, the
compressed archive of every tangent midpoint and Flint radius has the same
SHA-256 in the failed and passing runs.  No geometric or spectral number was
changed; only the invalid binary64 post-processing floor changed.

## Blind spectrum

**DERIVED COMPUTATIONAL, but not yet a physical stability spectrum.**  In
both schedule parities there are exactly

```text
119 resolved expanding eigenvalues,
119 reciprocal resolved contracting eigenvalues,
238 resolved off-unit eigenvalues in total.
```

The resolved moduli occupy

```text
contracting: 0.0210862 ... 0.117468
expanding:    8.51297  ... 47.4244.
```

The remaining eigenvalues are unit-consistent or numerically open under the
conservative non-normal eigensolver gate:

```text
even: 1184 unit-consistent, 18 open
odd:  1150 unit-consistent, 52 open.
```

Those different open counts do not represent different spectra.  The two
operational spectra match at `1.1e-10`--`4.7e-10`; the label changes because
the eigenvector condition number changes the conservative uncertainty.  In
particular, the trivial homogeneous block contains one unresolved reciprocal
pair

```text
0.99390373, 1.00613367.
```

No claim of exact unit modulus is made for an `OPEN` pair.

## Schedule comparison

**DERIVED COMPUTATIONAL:** all seven blocks receive the preregistered
`SCHEDULE_ROBUST` spectral label.

The eigenvalue comparison alone is weak because non-normality gives empirical
uncertainties `1.1e-4`--`8.0e-4`.  The independent ordered singular spectra
are more informative: even/odd distances are
`5.82e-11`--`2.33e-10`, versus uncertainties
`1.18e-9`--`2.29e-9`.  Hence schedule robustness does not rest only on the
ill-conditioned eigenvectors.

This is robustness of spectral invariants.  The even and odd matrices are
not equal in the present orbit bases, and no full `1440 x 1440` carrier
conjugacy was proved here.  That stronger statement remains **OPEN**.

The result also reconciles an old apparent conflict.  The first reduced
tangent artifact labelled its spectra schedule-dependent because its
eigenvalue uncertainty omitted non-normal conditioning.  A later exhaustive
reduced conjugacy audit already proved the reduced maps directly conjugate.
The present full result agrees with that stronger control rather than with
the obsolete sensitivity label.

## The 120-versus-119 observation

**PATTERN, POST-RESULT.**  The previously certified weak pre-Legendre Schur
sector has dimension 120.  The blind tangent spectrum has 119 strong
hyperbolic reciprocal pairs plus one unresolved near-unit reciprocal pair in
the homogeneous sector:

```text
119 + 1 = 120.
```

The count holds sector by sector.  Each minimal sector of irrep dimension
`d` contains `5d` candidate pairs: all `5d` are strongly off-unit except in
the trivial sector, where four are strong and the collective pair is open.
This is too structured to ignore, but it was noticed after the spectrum and
therefore is not evidence yet that the hyperbolic modes *are* lapse or
pseudo-constraint modes.

A falsifiable follow-up is now clear: construct the geometry-selected
`5d` pole/lapse response inside every sector before looking at eigenvectors,
form its canonical phase completion, and compare it with the `10d`
stable-plus-unstable spectral subspace.  Exact identification would demote
the 119 large pairs from candidate gravitational instabilities to the weak
lapse/pseudo-constraint sector.  Separation would make them a genuine
physical warning.  No spatial Laplacian or desired wave spectrum should be
loaded before this distinction is settled.

## Physical status ledger

- **DERIVED COMPUTATIONAL:** a full, regular, action-generated, symplectic
  one-slab tangent exists on all 1,440 boundary phase variables.
- **DERIVED COMPUTATIONAL:** 119 strong expanding/contracting reciprocal
  pairs occur in both schedule spectra.
- **DERIVED COMPUTATIONAL:** the spectral invariants pass the frozen
  schedule-robustness criterion.
- **PATTERN:** the `119+1=120` match with the weak vertex-lapse dimension.
- **OPEN:** whether the strong pairs are lapse/pseudo-constraint artifacts,
  physical scalar/vector/tensor perturbations, or a mixture.
- **OPEN:** a gauge-invariant curvature-mode quotient on this curved dust
  background.
- **OPEN:** a second dynamically solved anisotropic slab, so no multi-tick
  propagation or long-time stability exists yet.
- **OPEN:** refinement convergence and comparison with continuous `S^3`
  perturbations.
- **OPEN:** dispersion, a physical clock, a limiting speed, Planck units and
  particle masses.

The present result is therefore real progress toward dynamics, but not yet a
graviton spectrum and not a derivation of `c`.

## Post-result prior-art check

A second primary-source search used the newly relevant terms `hyperbolic
modes`, `curved-background Regge Hessian`, `pseudo-constraints`, `lapse
modes` and `symplectic transfer map`.  The located primary works continue to
support the established interpretation:

- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) explain that curvature
  breaks exact discrete gauge symmetry and produces background-dependent
  pseudo-constraints;
- [Hoehn](https://arxiv.org/abs/1411.5672) identifies lattice gravitons only
  after a gauge analysis around flat backgrounds;
- [Dittrich--Freidel--Speziale](https://arxiv.org/abs/0707.4513) connect
  Hessian null modes with remnant diffeomorphisms and require a continuum
  limit for graviton interpretation.

No located source supplies the present full dust-slab spectrum or proves the
`119+1` identification.  A search is not a novelty proof; external novelty
remains **OPEN**.
