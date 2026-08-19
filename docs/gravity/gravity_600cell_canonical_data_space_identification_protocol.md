# Protocol: identify the modular 120+120 compatible data space

Date: 2026-08-19

Freeze this target-disclosed protocol before constructing `E U` or computing
its first global rank.

## Frozen inputs

- complete compatibility source and target-blind artifact from the canonical
  data admissibility mission;
- projection artifact committed as `9b97775`, SHA-256
  `f011ef9848a6139408a9f8495a12e0d8e0050e04f39aa5c00ca88c02dde26beb`;
- projection protocol deviation committed as `a3fe2b9`;
- prior-art and framing gate
  `gravity_600cell_canonical_data_space_identification_prior_art.md`.

Use the same representatives `(lambda,tau)=(2,5),(3,11)`, primes 1000003
and 1000033, baseline and alternate right-inverse graphs, reversed faces, odd
canonical relabelling, and reversed metric sign.

## Target disclosed before calculation

The proposed compatible boundary-data space is

\[
\operatorname{im}U\oplus\mathbb Q^{120}_{\rm strut},
\]

where `U` is the unsigned 600-cell vertex-edge incidence map.  On edge
`{u,v}`, its physical squared-length coefficient is
`8 lambda (sigma_u + sigma_v)`.  No coefficient may be fitted after a rank is
seen.

This target concerns data coordinates only.  It does not claim that the
refuted local 3600-cell-flex lift is correct.

## Exact construction

For every complete equation row:

1. retain its 3600 `F` coefficients;
2. compose every upper-edge coefficient with `U`, producing 120 vertex-scale
   columns and the exact factor `8 lambda`;
3. compute `rank_p([F E U])`;
4. independently retain the 120 strut columns and compute `rank_p([F S])`.

The candidate edge image is compatible iff `rank_p([F E U])=3600`.  The
strut ambient space is compatible iff `rank_p([F S])=3600`.  Verify
`rank(U)=120` directly from the exact graph and by the triangle/connectedness
argument in the prior-art note.

## Negative construction attack

Let `{u0,v0}` be the lexicographically first sorted 600-cell edge.  Construct
`U_bad` by replacing only that row's entries `(u0:1,v0:1)` with `(u0:1)`.
Before global inclusion is evaluated, require

```text
rank(U) = 120
rank(U_bad) = 120
rank([U U_bad]) = 121
```

over the rationals and both primes.  Compute `rank_p([F E U_bad])`.  If the
candidate is included and the corrupted distinct image is also included, the
frozen 120-dimensional `K_E` census is contradicted and the run is a control
failure, not a second positive.

## Required controls

1. All input hashes and the frozen `(240,120,120,120,120)` projection tuple
   reproduce byte-for-byte.
2. Every local/face construction control remains true.
3. `rank_p(F)=3600` and `rank_p([F E])=4200` reproduce for each construction.
4. The exact incidence, corrupted-incidence, and joint ranks above hold.
5. Candidate inclusion decisions agree across primes, representatives, exact
   right-inverse graphs, reversed faces, odd relabelling, and metric sign.
6. The corrupted image is rejected whenever the candidate image is accepted.
7. Record all ranks even for a negative outcome.

## Outcome hierarchy

- `CANONICAL_DATA_SPACE_CONTROL_FAILED`: provenance, construction, frozen
  dimensions, incidence ranks, or a conditional corruption control fails.
- `CANONICAL_DATA_SPACE_DISAGREEMENT_OPEN`: legitimate constructions or
  primes disagree on candidate inclusion.
- `CANONICAL_DATA_VERTEX_SCALE_IMAGE_REFUTED`: the complete agreeing tests
  reject `im(U)`.
- `CANONICAL_DATA_MODULAR_VERTEX_SCALE_PLUS_ARBITRARY_STRUTS`: every complete
  test includes `im(U)` and the whole strut space, while rejecting `U_bad`.

The positive branch proves equality only over the two frozen finite fields,
using inclusion plus the already frozen dimension.  Rational equality, the
global cell-flex lift, action dynamics, lapse selection, tensor propagation,
`c`, `G`, and Planck units remain **OPEN**.

## Reproducibility and acceptance

Register the verifier before first execution and run only it.  Freeze the
first artifact before attempting a rational lift.  Because the preceding
projection mission has a recorded blindness deviation, no material carrier
claim is accepted without a later mechanically independent exact
implementation that does not reuse this modular rank composition.
