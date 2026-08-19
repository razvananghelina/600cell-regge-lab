# Consolidated result: pseudo-longitudinal non-invariance persists for two ticks

Date: 2026-08-19

## Verdict

**DERIVED COMPUTATIONAL / STRUCTURAL:** on both accepted centered backgrounds
of the fixed 600-cell recurrence, the 15-dimensional tangential
vertex-displacement image in each of sectors 4 and 5 is not invariant under
the generalized action stiffness.  The shifted result holds for both
parities and four derivative schedules: 16 of 16 cells in two mechanically
different constructions.

This is temporal persistence on two finite, nearly homothetic backgrounds.
It is not a curvature limit, refinement result, conserved quantity, physical
instability, propagating mode or continuum gauge statement.

## Provenance ledger

| stage | commit | literal outcome |
|---|---|---|
| prior-art/framing gate | `d1d7193` | shifted residual still OPEN |
| primary protocol | `7f7e06e` | target-disclosed temporal test |
| primary registered verifier | `9c29111` | no shifted result yet |
| primary artifact | `93d775e` | `PERSISTS`, `10/10` |
| direct adversarial protocol | `1bfe9e9` | archive-free SVD route |
| direct registered verifier | `3ea94a0` | no direct result yet |
| direct artifact | `f0f0625` | literal `OPEN`, `18/18` |
| disclosed classifier analysis | `b0549a9` | redundant rank threshold identified |
| logical correction protocol | `5ba7830` | target-disclosed exact-rank correction |
| correction registered verifier | `04d5411` | original OPEN preserved |
| correction artifact | `27e103d` | residual confirmation, `7/7` |

Artifact hashes:

```text
primary shifted
0480f5d49d24e0f5d8e4e95f0cf62b7d0d9242459ed2b8f6d8e835ecd6e103a7

direct adversarial OPEN
9e9f7253fd10422f3534914fae020857162862123fd4eae889e3570083552179

rank-classifier correction
42ef59f9c31a3bc16c78c0964b0450f5e62288fe6585ad460814f104face2eb3
```

## Independent numerical evidence

The primary route loaded the shifted centered binary64 archive, rebuilt the
carrier and used column-pivoted QR.  The adversarial route loaded no centered
`M,V` archive: it reconstructed slabs 2 and 3 from local 4-simplex Hessians
at 100 decimal digits, rebuilt exact golden-ratio vertices and all 720 edges,
and used full SVD.  It opened the primary residual artifact only after its
direct 16-cell census.

Both routes returned

```text
rho_span  NONZERO_RESOLVED  16/16
rho_comm  NONZERO_RESOLVED  16/16
inertia   15 negative + 10 positive per cell
```

For `even_sector4/operational_primary`:

```text
current rho_span          0.0777626757830970
shifted direct rho_span   0.0777626677060581
shifted/current ratio     0.999999896132189

current rho_comm          0.00901883801272048
shifted direct rho_comm   0.00901889333050163
shifted/current ratio     1.00000613358185
```

The direct values agree with the shifted archived route to about two parts in
`10^8`.  The near equality between ticks is a **PATTERN**, not a conservation
law: the residuals are normalized and the recurrence is close to homothetic.

## Preserved negative about the classifier

The direct protocol conjunctively required an auxiliary augmented numerical
rank greater than 15.  Its conservative threshold was about `4.25805e-5`,
above the sixteenth augmented singular value `2.96331e-5`, so the literal
direct outcome remained `OPEN`.  This is preserved.

The later correction did not lower that threshold.  It used the exact lemma

```text
rank([X,Y]) = rank(X) + rank((I-P_X)Y)
```

with `X=BL`, `Y=AL`.  A resolved nonzero span residual is already the witness
that `Y` leaves `im X`; asking a second threshold to establish the same fact
was redundant.  The correction therefore confirms persistence of
non-invariance while retaining `PRESERVED OPEN` for the original numerical
rank conjunction.

## Scientific status

- **DERIVED COMPUTATIONAL / STRUCTURAL:** the 30 negative modes remain
  pseudo-longitudinal rather than exact gauge on two consecutive accepted
  backgrounds.  Quotienting them away is not authorized.
- **PATTERN:** their dimensionless leakage is almost unchanged between the
  two backgrounds.
- **OPEN:** whether the leakage tends to zero under controlled flattening or
  projected refinement.
- **OPEN:** whether the negative subsystem is a physical instability,
  pseudo-constraint artefact or coarse-lattice artefact.
- **OPEN:** tensor quotient, wave equation, polarizations, dispersion and
  effective speed.

## Next decision

Another tick at the same resolution has low evidential value.  The decisive
test now requires a genuine parameter that changes geometry:

```text
rho_span(curvature, resolution),
rho_comm(curvature, resolution).
```

The repository's projected red-refinement tower currently certifies only the
homogeneous scale/lapse dynamics.  It does not contain the full inhomogeneous
boundary-edge Hessian or a canonical prolongation for these 30 modes.  Thus a
true refinement test is feasible but is a new sparse/symmetry-reduced
construction, not a reuse of the existing homogeneous verifier.  Until that
construction is preregistered and completed, continuum recovery remains
OPEN.

