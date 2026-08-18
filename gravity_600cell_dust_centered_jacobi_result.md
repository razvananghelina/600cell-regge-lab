# Result: centered Regge stencil and a blind `1:5` kinetic-signature pattern

Date: 2026-08-18

## Provenance

```text
prior-art gate                 0b69556
blind protocol                 20cf463
registered implementation     18d3574
blind passing artifact         a2d0ce0
```

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_centered_jacobi.py`.

Artifacts:

```text
reproducible/gravity_600cell_dust_centered_jacobi.json
SHA-256 fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56

reproducible/gravity_600cell_dust_centered_jacobi.npz
SHA-256 1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef
```

The first targeted run completed in `51.27 s` with

```text
7/7 PASS
CENTERED_JACOBI_CERTIFIED.
```

A second complete targeted run took `53.32 s`, again returned `7/7 PASS`,
and reproduced both artifacts byte for byte with the same SHA-256 hashes.

No spatial operator, desired inertia, continuum harmonic spectrum or speed
was loaded.  The full suite was not run.

## Certified centered equation

For every schedule, sector and derivative variant, the committed Jacobi
operator was decomposed as

```text
M = (K_- + K_+)/2,
N = (K_+ - K_-)/2,
V =  K_- + K_0 + K_+,
```

giving exactly

```text
M (delta q_2 - 2 delta q_1 + delta q_0)
+ N (delta q_2 - delta q_0)
+ V delta q_1 = 0.
```

All `56` Flint determinant balls of `M` exclude zero.  All original and
normalized recurrence identities hold entrywise in ball arithmetic.  Hence

```text
Gamma = M^-1 N,
Omega = M^-1 V
```

exist throughout the complete declared carrier.

**DERIVED COMPUTATIONAL:** the first two Regge dust slabs select a regular
finite “second difference + first difference + stiffness” equation on all
`720` boundary edge-metric perturbations.  No coefficient was fitted.

## Near-Hermiticity under the literal slice identification

The operational midpoint adjoint-defect ratios over all sectors are

```text
M : 8.85e-15 ... 1.71e-14,
N : 7.26e-14 ... 1.51e-13,
V : 1.09e-13 ... 2.23e-13.
```

This was not assumed by the protocol.  The earlier future/past asymmetry of
about `0.22--0.24` is therefore predominantly a difference between two
nearly Hermitian coefficients, not a large anti-Hermitian contamination.

**DERIVED COMPUTATIONAL:** in the declared logarithmic edge coordinates, the
centered coefficients are Hermitian to the calibrated construction floor.
An analytic equality and a coordinate-free superspace transport remain
**OPEN**.

## Blind inertia

Every minimal sector of irrep dimension `d` has exactly

```text
5d  resolved eigenvalues of one sign,
25d resolved eigenvalues of the other sign,
0 zero-consistent,
0 open.
```

After restoring representation multiplicity, each schedule separately has

```text
120 positive, 600 negative, 0 zero, 0 open
```

for the Hermitian part of `M`.  The JSON aggregate prints twice these values
because it records both independent schedules; there are not `1440` physical
position variables in one schedule.

The sector extrema span approximately

```text
-554.04 <= lambda(H_M) <= 164.86,
```

with inertia errors between `9.59e-8` and `2.31e-6`.  The sign decisions are
not marginal.

**DERIVED COMPUTATIONAL:** the centered kinetic coefficient is indefinite.
Multiplying the entire action by `-1` exchanges the two signs but cannot
remove the `120:600 = 1:5` split.

This is not automatically a ghost theorem.  General relativity itself uses
an indefinite DeWitt supermetric, and the present `M` connects discretely
identified time fibres before a constraint quotient.  Positivity of all raw
edge variables was never a correct continuum requirement.

## A striking post-result pattern

At each point of a three-dimensional continuum geometry, a symmetric metric
perturbation has six components: one trace/conformal direction and five
traceless directions.  The DeWitt supermetric has the corresponding `1:5`
indefinite signature, up to overall sign convention.

The blind result

```text
720 edge variables = 120 * 6,
inertia             = 120 * (1 sign) + 120 * (5 opposite signs)
```

therefore has precisely the DeWitt count.

**PATTERN, POST-RESULT:** the centered Regge coefficient has the global count
expected from one conformal plus five traceless metric components per 600-cell
vertex.

This is not yet an identification.  The simplicial supermetric literature
shows that discrete signatures can change and can contain additional
physical negative directions.  A dimension match cannot prove that the
`120`-dimensional sign subspace is conformal.

There is, however, a canonical falsifier.  Vertex scaling defines the
unfitted discrete-conformal map

```text
(C sigma)_(uv) = sigma_u + sigma_v
```

from the `120` vertex scalars to the `720` edge logarithms.  The 600-cell
contains triangles, so this unsigned incidence has rank `120`.  Comparing
`im C` with the `120`-dimensional inertia subspace can prove or refute the
DeWitt reading without choosing Schur coefficients.

## Generalized drift and stiffness

The operational eigenvalues are numerically real to the displayed arithmetic
floor.  For `Omega` the real parts lie in

```text
-7.0635e-6 ... 0.8809942,
```

and the largest imaginary midpoint is below `2.1e-14`.  All `720` eigenvalues
per schedule are `REAL_CONSISTENT` under the preregistered Bauer--Fike error;
none is resolved complex or open.

This does not certify nonnegativity.  The eigenvalue errors range from about
`1.75e-7` to `4.05e-5`, so the small negative values are not uniformly
resolved away from zero.  The singular condition numbers of `Omega` are
large, about

```text
3.86e5 ... 1.37e6,
```

although its eigenvector condition numbers are only about `3.45--23.53`.

For `Gamma`, the eigenvalue real parts lie between approximately
`-0.66073` and `-2.35e-6`.  No Hubble-friction interpretation is made: the
small end is within some eigenvalue errors and no proper-time normalization
has been supplied.

- **DERIVED COMPUTATIONAL:** the finite normalized operators exist and have
  schedule-stable, reality-consistent spectra.
- **STRUCTURAL:** these are coefficients per declared tick and edge
  coordinate.
- **OPEN:** real analytic spectrum, sign of the near-zero stiffness modes,
  continuum wave frequencies and physical damping.

## Schedule robustness

All primary comparisons pass:

```text
14/14 singular spectra SCHEDULE_ROBUST,
7/7 Omega eigen spectra SCHEDULE_ROBUST.
```

The largest distance/error ratio is below `6.53e-9` for singular spectra and
below `3.60e-9` for the secondary eigen spectra.  Thus neither the inertia nor
the displayed normalized spectrum is being selected by one staircase parity.

## Prior-art reconciliation

The continuum `1:5` DeWitt signature is standard.  A useful primary discrete
warning is [Hartle--Miller--Williams](https://arxiv.org/abs/gr-qc/9609028):
the Lund--Regge metric on squared edge-length configuration space always has
at least one timelike direction, but its complete signature can depend on the
triangulation and point, and additional negative directions need not be
gauge.  [Williams' review](https://arxiv.org/abs/gr-qc/9702006) likewise
reports degeneracy and signature change in simplicial superspace.

Therefore the count found here is consistent with continuum superspace but
is not guaranteed by generic Regge calculus.  No located source computes the
present two-slab centered coefficient or its exact `120:600` split.  External
novelty remains **OPEN**.

## Physical status

- **DERIVED:** regular centered finite stencil on all edge perturbations.
- **DERIVED:** nearly Hermitian `M,N,V` and real-consistent `Omega` spectrum.
- **DERIVED:** resolved inertia `120:600` in each schedule.
- **PATTERN:** exact match to the continuum DeWitt `1:5` local count.
- **OPEN:** whether the 120-dimensional sign subspace is the canonical
  vertex-conformal image.
- **OPEN:** constraint quotient and identification of the five traceless
  directions with two propagating tensor polarizations.
- **OPEN:** refinement, physical time, dispersion, `c` and Planck units.

## Next falsifiable gate

Construct the unsigned vertex--edge incidence map directly from the literal
600-cell and project it into the same seven symmetry sectors before reading
any inertia eigenvector.  Compare its `5d`-dimensional image with the
precommitted `5d` sign subspace in every sector and both schedules.

- Equality would upgrade the `1:5` count to a finite DeWitt-type conformal
  decomposition.
- Resolved separation would demote it to a numerical coincidence.
- Partial overlap would show that the raw kinetic sign mixes conformal and
  shape variables, blocking the simple continuum interpretation.

No continuum harmonic or desired wave speed belongs in that test.
