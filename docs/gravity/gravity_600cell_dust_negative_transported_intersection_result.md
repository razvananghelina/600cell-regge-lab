# Result: the transported negative-stiffness intersection is zero

Date: 2026-08-18

## Headline

The source-certified negative-stiffness phase fibers have zero transported
intersection under the canonical second-slab tangent in every frozen cell:

```text
16/16 cells: rank((I-Q_1) T_2 Q_0) = 30,
16/16 cells: dim(F^-_0 intersection T_2^(-1)(F^-_1)) = 0.
```

The exact primary result is reproducible and the initially inconclusive
binary64 audit has been resolved by a separate fixed-input calculation at
`100` and `140` decimal digits.  The consolidated classification is
**DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED**.

This closes the Euclidean negative-fiber phase route relative to the frozen
carrier Hilbert metric.  It does not close full Regge dynamics or establish
anything about physical pre/post constraint surfaces not yet derived from the
action.

## Complete provenance ledger

| stage | commit |
|---|---|
| primary prior-art and canonicity gate | `65419f6` |
| primary rank protocol | `2927fb9` |
| registered primary verifier | `e5c168d` |
| reproducible primary artifact | `5cc576e` |
| adversarial independence gate | `8d532e8` |
| adversarial binary64 protocol | `411df3f` |
| registered adversarial verifier | `80d9854` |
| preserved adversarial OPEN artifact | `899a98c` |
| roundoff-resolution framing gate | `cf64f34` |
| exact-dyadic framing correction | `ba42d40` |
| roundoff-resolution protocol | `b54f80a` |
| registered roundoff resolver | `f032c18` |
| reproducible resolution artifact | `d9e2bc4` |

No new intersection singular value was computed before the primary protocol
commit.  The adversarial target was necessarily disclosed post-result, but
its algorithm and thresholds were committed before its first execution.  The
roundoff resolver was likewise committed only after preserving the first
audit's OPEN outcome, and did not modify that audit's frozen threshold.

## Result A: source-certified rank

The primary verifier reconstructed, from the existing Flint source balls:

- all `32` old/shifted ordinary negative-stiffness projectors;
- the separated `15/10` stiffness split in every case;
- all `16` exact second-slab tangent balls;
- all `60` singular values of the full residual in each cell.

Every one of the `32` projectors overlaps its independently committed binary
control.  Their source-certified projector errors lie between approximately
`3.66e-55` and `3.71e-54`; the weakest sign separation exceeds approximately
`6.00e51` error units.

The complete singular census is

```text
480 SINGULAR_NONZERO_RESOLVED,
480 SINGULAR_ZERO_CONSISTENT.
```

The lower `30` singular values are the structural zeros forced by the exact
rank-`30` right factor `Q_0`.  The upper `30` are all nonzero-resolved.  The
smallest resolved value is approximately

```text
3.9917386e-7,
```

while the complete source-and-projector residual error is at most about
`1.09e-48`.  The weakest nonzero singular value is therefore about
`3.67e41` complete error units from zero.  Weyl's inequality plus the exact
upper bound proves rank exactly `30`, not merely a numerical rank.

Two complete executions reported `10/10` and wrote byte-identical artifact

```text
c490431bdaeae3026692cd358f60d0b47ef5d63aa59217e400daac807ed21be0.
```

## Result B: the first independent audit was honestly inconclusive

The mechanically different binary64 audit used the earlier binary projectors,
fresh SciPy spectral bases, the earlier tangent archive, a complete-QR target
complement and the square leakage

```text
L = W_1,perp^* T_2 W_0.
```

All controls and convention stresses passed, and its midpoint singular values
were numerically compatible with the primary spectrum.  However, its frozen
roundoff envelope was

```text
e = 1000 eps_machine 60 max(1,||T_2||,||L||)
  = approximately 3.54e-6,
```

larger than the weakest singular directions.  It classified, per cell,

```text
20 nonzero-resolved, 4 open, 6 zero-consistent,
```

and correctly returned

```text
ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN.
```

This was not overwritten or reinterpreted as a passing audit.  Two executions
reported `11/11` with byte-identical artifact

```text
d5074507326bb981ad7573bd562c1aa9f0af4e1eb6b6924e3ac959a5fa1d3340.
```

The lesson is methodological: an intentionally conservative error floor can
be too broad to falsify a small singular value even when independent numeric
implementations agree on that value.

## Result C: fixed-input roundoff resolution

The resolution test treated every stored float of the independent binary
projectors and tangent midpoint as its exact dyadic rational, reconstructed
their spectral fibers at both `100` and `140` decimal digits, and used the
explicit bottom spectral space as target complement.

For all `16` cells and both precisions:

- all `30` square-leakage singular values were nonzero-resolved;
- direct SVD agreed with the Gram spectrum;
- `1/||L^-1||` independently recovered the smallest singular value;
- full-intersection and zero-intersection controls discriminated;
- reversal/rephasing preserved the spectrum;
- the complete `100`- and `140`-digit spectra agreed within the preregistered
  residual bound.

The complete census is

```text
precision 100: 480 FIXED_SINGULAR_NONZERO_RESOLVED,
precision 140: 480 FIXED_SINGULAR_NONZERO_RESOLVED.
```

The weakest fixed-input value exceeds its residual bound by approximately
`5.01e62`.  The largest observed difference between the two precisions is
about `1.18e-97`, while the smallest permitted stability error is about
`1.59e-69`.

Two executions reported `11/11` and wrote byte-identical artifact

```text
83b5cca5ba1b4cc81c5e2e3be2a9405837df23e24b667690f0035156683dfc0e.
```

Thus the first audit's open labels arose from its global float64 envelope, not
from a rank loss in its fixed binary data.  This resolution does not claim
that binary midpoints enclose the source; that role remains with Result A.

## Why the full kernel vanishes despite the earlier block pattern

The earlier blockwise closure test found `A,C` leakage zero-consistent under
very conservative binary projector errors, while `B,D` leakage was
nonzero-resolved.  That pattern did not prove a surviving half-dimensional
kernel.  The complete leakage combines all four blocks.  At source-certified
precision, their coupling has full column rank `30`; even the weakest coupled
direction is nonzero.

Equal configuration fibers and a persistent negative inertia therefore do
not imply a common propagating phase fiber.

## Canonicity and physical scope

An indefinite stiffness form intrinsically fixes its inertia but not a unique
maximal negative subspace without a positive metric.  Here the negative
spectral fiber is canonical only relative to the positive Hermitian product
inherited from the frozen orthonormal binary-action carrier.  No coefficient,
alignment, schedule or target was fitted, but this extra metric remains
**STRUCTURAL** rather than an action-derived pre/post constraint condition.

Consequently:

- **DERIVED COMPUTATIONAL:** the frozen stiffness form has a separated
  `15/10` inertia in all old/shifted cells;
- **STRUCTURAL:** the particular negative spectral fiber is selected relative
  to the frozen carrier Hilbert metric;
- **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED:** its transported
  intersection is zero in all `16/16` cells;
- **CLOSED:** this negative-spectral phase-carrier route;
- **CLOSED:** the earlier generalized kinetic--stiffness phase-carrier route;
- **OPEN:** action-derived pre/post constraint surfaces, a constraint-reduced
  Jacobi carrier, full anisotropic perturbation propagation, refinement and
  continuum limits;
- **NOT COMPUTED:** graph or Lagrangian structure, because the present
  intersection is zero;
- **NOT ESTABLISHED:** graviton modes, dispersion, inertia, particle mass or a
  limiting speed.

The full anisotropic Legendre map remains regular and the full canonical
tangent still exists.  The result says that two proposed spectral subbundles
do not carry nonzero data across this tick.  It does not say that the full
Regge system has no dynamics.

## Post-result literature check

The refined search used the terms transported spectral intersection,
discrete Jacobi map, negative eigenspace and pre/post constraints.  It found no
primary source studying this exact `600`-cell dust construction.

The closest relevant framework remains:

- Dittrich and Hoehn, *Canonical simplicial gravity*, DOI
  `10.1088/0264-9381/29/11/115009`, arXiv:`1108.1974`;
- Dittrich and Hoehn, *Constraint analysis for variational discrete systems*,
  DOI `10.1063/1.4818895`, arXiv:`1303.4294`;
- Dittrich and Hoehn, *From covariant to canonical formulations of discrete
  gravity*, DOI `10.1088/0264-9381/27/15/155001`, arXiv:`0912.1817`.

Those works support the caution rather than a physics claim: propagating data
in a singular discrete theory are selected by action-derived pre/post
constraint surfaces, and curved Regge backgrounds need not retain flat-space
gauge directions.  External novelty remains **OPEN** pending expert review.

Only the three targeted verifiers described above were run.  The full suite
was deliberately not run.
