# Result: a canonical nested tangential displacement carrier exists

Date: 2026-08-22

## Verdict

> **DERIVED COMPUTATIONAL / STRUCTURAL INFRASTRUCTURE:** the selected radial
> edgewise refinement defines a canonical injective prolongation from all
> spatial tangential vertex displacements of `K0=P(sd K_600)` to those of
> `K1=P(Esd_2(sd K_600))`.

The map is fixed by normalized midpoint geometry, is independent of parent
chamber and temporal staircase, is covariant under both determinant classes
of `O(4)`, and has an exact old-vertex left inverse.  Its rank is therefore
`7920` without a numerical rank threshold.

This is a positive feasibility result.  It is not a restored Regge
constraint, a canonical momentum lift or a physical graviton carrier.

## Provenance

| stage | commit |
|---|---|
| prior-art and repository audit | `79b612b` |
| target-free protocol | `9005af2` |
| verifier registered before first execution | `077164e` |
| preserved `13/14` first failure | `acbf74d` |
| precision correction protocol | `db6845a` |
| corrected implementation | `74438af` |

Accepted verifier and artifact:

```text
reproducible/verify_gravity_600cell_nested_vertex_displacement.py
SHA-256 2b4351258e0e2f0f1a6fccb61abd6d8f11e574e1975ed28f07fb34d69088250d

reproducible/gravity_600cell_nested_vertex_displacement.json
SHA-256 ef7f565fd54487885e9de459d97fbe61af8dcc6e9add30768ae2fe94ca7d4250
```

The corrected targeted verifier passed `14/14` twice and produced a
byte-identical artifact.  No action Hessian, spectrum or full suite was run.

## Two carrier constructions

Route A traversed the eight exact `r=2` edgewise children of every ranked
barycentric chamber and merged `460800` fine-vertex occurrences by their
unordered repeated-vertex keys.

Route B did not use the edgewise-facet traversal.  It enumerated the edges of
`K0` directly and formed

```text
{(i,i): i in V(K0)} union E(K0).
```

The routes agreed exactly and reconstructed

```text
K0 f-vector = (2640,17040,28800,14400),
K1 f-vector = (19680,134880,230400,115200),
fine vertices = 2640 retained + 17040 projected midpoints.
```

Every fine triangular face has incidence two.  The minimum denominator
`||x_i+x_j||` is `1.962455...`, so no normalized midpoint is singular.

## Canonical derivative

For `y=P(x_i+x_j)` the prolongation is

```text
U = (I-y y^T)(u_i+u_j)/||x_i+x_j||.
```

All parent occurrences agree, with maximum coordinate discrepancy zero and
maximum tangent discrepancy `7.006e-18`.  Coarse and fine tangency residuals
are respectively below `8.24e-18` and `5.21e-18`.

Every old key `(i,i)` occurs exactly once and its derivative is `u_i`.
Restricting the fine field to the old keys therefore gives

```text
R_old T = I.
```

This is an exact index-level left inverse, so the rank on the direct sum of
the `2640` three-dimensional tangent spaces is

```text
rank(T) = 3*2640 = 7920.
```

No claim is made here about the rank after converting this carrier to edge
lengths.

## Independent differential controls

Centered differentiation of the complete nonlinear normalized construction
gave maximum fine-vertex component errors between `2.09e-14` and
`8.12e-14`.  On all `134880` fine squared chord lengths, the errors were
between `1.57e-13` and `7.25e-13`.

The first binary64 execution had edge errors of order `10^-9` which grew as
the step decreased.  It correctly failed the preregistered convergence rule
despite being far below the absolute tolerance.  That `13/14` artifact is
preserved.  The correction changed only the finite-difference control path to
extended precision, retaining both steps, formulas, carriers, probes and
thresholds.

Two fixed orthogonal transformations, with determinants `+1` and `-1`, gave
maximum covariance residual `3.331e-16`.  A corrupted midpoint weight was
separated from the canonical point by `0.0527397...`, and deleting one old
key destroyed the left-inverse census.

## Temporal schedules and the phase-space warning

At both levels all 24 temporal schedules have the same spatial boundary
carrier, so the prolongation reads no schedule.  The verifier simultaneously
retains the accepted fact that the 24 internal temporal edge sets are
distinct and have empty total intersection.  The result therefore does not
select or average a temporal staircase.

More importantly, configuration injectivity does not reverse the direction
of cotangent pullback.  The already certified homogeneous coarse-to-fine
momentum lift still has a five-dimensional affine ambiguity.  The present
result must not be used to import a coarse momentum into `K1` by transpose,
pseudoinverse, minimum norm or another undeclared choice.

## Scientific status

- **DERIVED COMPUTATIONAL:** the two finite carriers, key census, parent
  independence, tangency, exact left inverse and complete edge derivative.
- **STRUCTURAL:** radial projection onto the declared round `S^3` and the
  interpretation of these tangential fields as shift-like displacement
  candidates.
- **CLOSED:** the claim that no canonical configuration carrier can be
  matched across `K0` and `K1`.
- **STILL OPEN:** normal/lapse displacement transport and a dynamically
  selected cotangent relation.
- **STILL OPEN:** an on-shell `K1` finite-height background with the same
  action and conserved-matter prescription as `K0`.
- **STILL OPEN:** convergence of pseudo-constraint couplings or singular
  values, an exact constraint quotient, gravitons, `c`, `G`, Planck units and
  particle physics.

## Next admissible gate

Do not compute a cross-resolution Hessian spectrum yet.  First construct and
preregister a matched `K1` on-shell finite-height seed using the same radial
edgewise carrier, Regge-plus-dust action, physical radius and conserved-matter
normalization as the accepted `K0` seed.  Its internal stationarity and
temporal-schedule handling must be certified before evaluating the Hessian on
the nested carrier.

If no such matched seed exists, constraint-restoration under this refinement
family is not a well-posed comparison and the route closes.  If it exists,
the nested carrier fixed here can be used target-blind; eigenvector matching
after spectral inspection remains forbidden.
