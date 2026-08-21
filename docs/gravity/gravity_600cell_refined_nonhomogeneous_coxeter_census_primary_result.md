# Primary result: all schedule kernels reduce to the duration line

Date: 2026-08-21

Status: **32/32 primary checks passed.**  Schedule 0 remains adversarially
corroborated; the eleven newly computed representatives remain **PRIMARY
ONLY** until a mechanically different replication is complete.

## Frozen provenance

- finite-group prior-art gate: `295e90d`;
- accepted schedule-0 result: `66e47a7`;
- target-free all-schedule protocol and structural manifest: `2eef0e1`;
- registered implementation before first execution: `b7bc50b`;
- result artifact:

  ```text
  reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_census.json
  SHA-256 5a50b8179fc272d75c0811dcf34fb8c6a464e564729956bdd6baaa82a6f058b6
  ```

Only the active verifier was run.  The full suite and the failed old sparse
factorization census were not rerun.  The static registry audit reports 388
registered and 388 distinct verifier names, with no duplicate, unregistered
or missing entry; two declared experiments remain deliberately skipped.

## Complete scope and hypotheses

The calculation covers all 24 labelled staircase orders through twelve
explicit time-reversal pairs on the fixed carrier `P(sd K_600)`.  Every
matrix uses:

```text
spatial f-vector           (2640, 17040, 28800, 14400)
complete slab edges        53760
fixed boundary edges       34080
internal directions        19680
input duration             tau0 = 0.0102
vertex masses              m_v = K_v/(8*pi)
coordinates                logarithmic signed squared edge lengths
```

For each representative `s`, the tested bordered matrix is

```text
K_s = [[C_s, n_s],
       [n_s^T, 0]],
```

where `n_s` is the analytic product-duration tangent.  All 24 complete source
matrices reproduce the CSR hashes frozen before the run.  Each forward matrix
is congruent to its explicitly assembled reverse matrix inside its declared
operator-error envelope.

The old artifact had outcome `LOCAL_EXTENSION_INVALID`; only its structural
schedule map and CSR hashes were extracted into the preregistered manifest.
No old eigenvalue or kernel field was loaded by this verifier.

## Exhaustive construction

For every representative, the fixed left Coxeter element has order 30 and
acts on the `19,680` internal edges as exactly

```text
656 cycles x 30 edges/cycle.
```

All sixteen independent cyclic sectors `k=0,...,15` were diagonalized in
full.  Conjugate weights give the exact bordered dimension `19,681` for every
schedule.  Trace, Frobenius/Parseval, explicit `k=1` versus `k=29`, source
digest, stationarity, reality, tangent-nullity, time reversal and covariance
checks all pass.  Omitting sector `k=15` is mechanically detected by the
weighted-dimension control.  A preregistered `1e-4` diagonal corruption is
detected for every representative.

Schedule 0 reproduces its frozen complete primary spectrum exactly at stored
binary64 precision: maximum difference `0`, under a `4.87324e-10` comparison
gate.

## Complete numerical census

The table reports the least zero-exclusion margin over all `19,681`
eigenvalues of each bordered representative.

| pair | schedules | order | weakest sector | eigenvalue | gate | margin |
|---:|:---:|:---:|---:|---:|---:|---:|
| 0 | 0/23 | 0123 | 1 | `1.455649e-9` | `2.424503e-10` | 6.0039 |
| 1 | 1/17 | 0132 | 1 | `1.455649e-9` | `2.424713e-10` | 6.0034 |
| 2 | 2/21 | 0213 | 1 | `1.455649e-9` | `2.532741e-10` | 5.7473 |
| 3 | 3/11 | 0231 | 1 | `1.455649e-9` | `2.835968e-10` | **5.1328** |
| 4 | 4/15 | 0312 | 1 | `1.455649e-9` | `2.728661e-10` | 5.3347 |
| 5 | 5/9 | 0321 | 1 | `1.455649e-9` | `2.835130e-10` | 5.1343 |
| 6 | 6/22 | 1023 | 1 | `1.455649e-9` | `2.425006e-10` | 6.0027 |
| 7 | 7/16 | 1032 | 1 | `1.455649e-9` | `2.424721e-10` | 6.0034 |
| 8 | 8/19 | 1203 | 1 | `1.455649e-9` | `2.532307e-10` | 5.7483 |
| 9 | 10/13 | 1302 | 1 | `1.455649e-9` | `2.730496e-10` | 5.3311 |
| 10 | 12/20 | 2013 | 1 | `1.455649e-9` | `2.428178e-10` | 5.9948 |
| 11 | 14/18 | 2103 | 1 | `1.455649e-9` | `2.426707e-10` | 5.9985 |

No representative contains a zero-compatible eigenvalue.  The global least
margin is `5.132811`, not a near-threshold rounding decision.

The almost common weakest eigenvalue and sector are a post-hoc **PATTERN**,
not an invariant established by this calculation.  The complete sorted
spectra are not identical: their maximum differences from schedule 0 reach
approximately `3.33e-5`.  No universality or physical interpretation is
claimed from the repeated weakest value.

## Verdict and limits

- **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED:** under the complete
  fixed hypothesis list, `ker(C_0) = span(n_0)`.
- **DERIVED COMPUTATIONAL, PRIMARY ONLY:** the same statement holds for each
  of the other eleven representatives within the printed error model.
- **DERIVED COMPUTATIONAL NEGATIVE, PRIMARY ONLY:** the exhaustive complete
  internal carrier contains no additional nonhomogeneous zero mode in any of
  the 24 labelled schedules, using time-reversal congruence.
- **OPEN:** mechanically different adversarial replication of the eleven new
  spectra.
- **OPEN:** an exact symbolic certificate and external novelty of this exact
  finite instance.
- **OPEN / NOT TESTED:** boundary-to-boundary propagation, gravitons, a wave
  equation, `c`, `G`, or Planck units.

This still does **not** produce a physical tick.  The value `tau0 = 0.0102`
is an input, and `n_s` is precisely the surviving product-duration direction.
The result says that all transverse internal directions are stiff at this
configuration.  It does not select a nonzero displacement along `n_s`, a
duration between boundaries, or a law of evolution.

## Next admissible step

Replicate the eleven new representatives using the already accepted
mechanically different construction: ambient `R4` reflections, geometrically
induced edge permutations, explicit sparse Fourier isometries `Q_k`, sparse
products `Q_k^* C_s Q_k`, and a different dense Hermitian solver.  The
adversarial protocol must be committed before that implementation is run.
Only after this replication passes should the internal-kernel question be
treated as closed and the separate boundary-evolution/tick-selection problem
resume.
