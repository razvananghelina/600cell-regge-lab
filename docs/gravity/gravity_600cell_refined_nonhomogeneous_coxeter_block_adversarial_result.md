# Adversarial result: the schedule-0 internal kernel is one-dimensional

Date: 2026-08-21

Status: **16/16 adversarial checks passed.**  Under the fixed hypotheses and
the declared numerical error model, the primary schedule-0 kernel result is
adversarially corroborated.

## Frozen provenance

- prior-art gate: `295e90d`;
- primary target-free protocol: `acfe795`;
- registered primary implementation: `69d2293`;
- frozen primary result: `fd75a55`;
- adversarial protocol: `b54f590`;
- registered adversarial implementation: `fa1ffc3`;
- adversarial artifact:

  ```text
  reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_blocks_adversarial.json
  SHA-256 3104422c0d2418061aad7fdfb7927a7ca6123239c2197c7b2779f2613103baf5
  ```

Only the active verifier was run.  The full suite and the failed old sparse
twelve-pair census were not rerun.

## Complete hypotheses

The claim concerns only the lexicographic schedule 0 on the fixed refined
carrier `P(sd K_600)`, with:

```text
19,680 internal slab edges,
tau0 = 0.0102,
m_v = K_v/(8*pi),
logarithmic signed-squared-edge coordinates,
fixed spatial boundary data,
the frozen binary64 local Hessian C0,
the analytic normalized product-duration tangent n0.
```

The tested bordered matrix is

```text
K0 = [[C0, n0],
      [n0^T, 0]].
```

For a symmetric matrix satisfying `C0 n0 = 0`, nonsingularity of `K0` is
equivalent to `ker(C0) = span(n0)`.  The computation certifies the premise and
nonsingularity only inside its printed floating-point error envelopes; it is
not a symbolic theorem about the transcendental Regge-action entries.

## Mechanically different reconstruction

The primary code found the order-30 action as a colour-preserving chamber
automorphism and accumulated Fourier blocks entry by entry.  Neither decisive
step was reused here.

The adversarial verifier instead reconstructed the four ambient `H4`
reflections from a geometric base flag in `R4`, induced their Coxeter element
on all coarse cells and then on internal slab edges, and built every block as

```text
B_k = Q_k^* (C0 Q_k)
```

with an explicit sparse Fourier isometry `Q_k`.  It diagonalized all sixteen
independent blocks using SciPy's divide-and-conquer Hermitian driver and
checked lifted eigenvector residuals against the original sparse matrix.
The primary spectrum was loaded only after the complete adversarial spectrum
had been constructed.

The geometric action again consists of exactly

```text
656 cycles x 30 edges/cycle.
```

The bordered weighted dimension is therefore `19,681`.  Forward and reverse
Coxeter words both reproduce the carrier action; the reverse word induces the
inverse edge permutation.

## Numerical result and controls

All `19,681` eigenvalues are separated from zero.  The least favourable block
is `k=1`:

```text
minimum |lambda| = 1.4556485274722962e-9,
zero-exclusion gate = 2.4565813337806738e-10,
margin ratio = 5.9255.
```

The largest lifted residual is `1.55547e-14`.  The complete sorted adversarial
and primary spectra differ by at most `7.10543e-15`, under their combined
comparison gate `4.90534e-10`.  Weighted trace and Frobenius/Parseval checks
agree exactly at the stored precision.

The preregistered negative controls also discriminate:

- one reflection has order 2 rather than 30;
- a `1e-4` diagonal corruption raises the covariance defect to
  `9.66667e-5`;
- the wrong Fourier-phase convention misses the same-sector block by
  `1.96368`, far above its `2.51879e-10` gate.

## Verdict

- **DERIVED COMPUTATIONAL / ADVERSARIALLY CORROBORATED:** for schedule 0 and
  the complete `19,680`-direction internal carrier, `K0` is numerically
  nonsingular in two mechanically different exhaustive calculations.
- **DERIVED COMPUTATIONAL NEGATIVE:** no additional nonhomogeneous internal
  zero mode exists for this schedule within the declared error model.
- **DERIVED COMPUTATIONAL:** consequently `ker(C0) = span(n0)` under the full
  hypothesis list above.
- **OPEN:** the other eleven time-reversal schedule representatives.  Nothing
  here permits extrapolation from schedule 0 to them.
- **OPEN:** an exact symbolic certificate and the external novelty of this
  particular finite instance.
- **OPEN / NOT TESTED:** boundary propagation, a graviton dispersion law,
  `c`, `G`, or Planck units.

Most importantly, this result does **not** derive a physical tick.  The number
`tau0 = 0.0102` remains an input.  The surviving vector `n0` is the known
product-duration tangent: showing that it is the only internal degeneracy
says that all transverse internal directions are stiff, not that the
dynamics selects a nonzero duration.  A derived tick would require a
boundary-to-boundary evolution problem whose equations select `tau` rather
than leave this line free.

## Next admissible question

The immediate mathematical completion is to reuse the already fixed
geometric Coxeter decomposition for the remaining eleven representative
schedules and certify each complete bordered spectrum.  That extension must
be preregistered before assembling or inspecting any new spectrum.  Only if
all twelve representatives have the same one-dimensional kernel is it
defensible to move to the separate physical question of whether boundary
data or matter conservation selects motion along that line.
