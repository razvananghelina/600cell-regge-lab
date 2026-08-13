# Disclosed search-result and implementation correction

Date: 2026-08-13

Root-search protocol commit: `5b687a3`

The first complete sequential run has now been observed:

```text
validated roots        0
even terminations      4 Jacobian-boundary, 2 iteration-limit
odd terminations       4 Jacobian-boundary, 2 iteration-limit
targeted checks        14/14
terminal ||r||_2       2.50579 ... 2.85780.
```

Inspection showed that every generic `jacobian_failure` occurred because one
centered point crossed the frozen lower logarithmic box boundary `y=-6`.
These eight searches must therefore be relabeled `boundary_contact`, as the
original protocol explicitly requires.  They are inconclusive, not evidence
of an interior minimum or root obstruction.

Before rerunning, freeze two implementation-only changes:

1. classify a centered Jacobian failure caused by `y +/- 2e-5` crossing the
   artificial box as `boundary_contact`;
2. evaluate the same 70 centered Jacobian points concurrently using the local
   eight-worker fork pool.  Candidate damping and backtracking remain serial,
   so the first accepted candidate and the entire mathematical trajectory are
   unchanged;
3. add an explicit consistency check that all 12 frozen searches terminate in
   a protocol-listed state;
4. retain the observed first-run counts and terminal residual interval in the
   JSON and final note.

No start, bound, damping value, backtracking factor, tolerance, causal gate or
iteration limit changes.  No root-search outcome has been used to add a new
start.  The full rerun must reproduce zero roots and the disclosed terminal
data to relative `2e-8`; otherwise the parallel implementation is rejected.
