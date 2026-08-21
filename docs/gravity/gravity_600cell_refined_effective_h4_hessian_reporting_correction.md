# Reporting correction: singular effective-Hessian gate

Date: 2026-08-21

The first result and byte-identical artifact are preserved in commit
`57da48b`.  All 24 internal blocks were singular, so the frozen protocol
forbade every Schur complement.  Three downstream controls were consequently
not executed, but generic boolean fall-through printed them as `PASS` with
zero samples.

This is not a scientific disagreement, but the labels could be mistaken for
evidence that an effective matrix had been tested.  Make only these reporting
changes:

1. when the internal singular gate fires, explicitly label Schur,
   directional, time-reversal/class and corruption branches as skipped by
   that gate;
2. require their record lists and effective-matrix lists to be empty;
3. serialize unavailable time-reversal, class-count, directional and
   corruption quantities as `null`, not `false` or zero;
4. print `NOT_COMPUTED_INTERNAL_SINGULAR` for the unavailable effective
   results.

Do not alter input hashes, finite-difference steps, precision, envelopes,
eigenvalues, singularity criterion or outcome hierarchy.  Rerun only this
targeted verifier twice and require byte identity.

