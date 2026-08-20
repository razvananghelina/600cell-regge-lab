# Result: no finite refined H4 stationary root was found

Date: 2026-08-20

Prior-art commit: `29162db`  
Protocol commit: `722fb3c`  
First implementation commit: `36ac45b`  
Accepted deterministic implementation commit: `1dfdfd4`

## Complete hypotheses

Use `K0=P(sd K_600)`, exact rank-derived geometry, equal fixed lower and upper
spatial boundaries, the corrected complex Lorentzian Regge action including
boundary terms, the previously selected total dust mass, conditional `P1`
weights and one representative of each of the 12 certified
schedule/time-reversal classes.

Search for zeros of all ten total-orbit internal log equations in six positive
cross-diagonal squared lengths and four positive lapse squares. The main box
is

```text
-0.35 <= y_cross <= +0.35,
-8 <= y_rho <= +2.
```

Use six frozen seeds per class and four independently initialized lower-lapse
boxes per class. This gives 72 main attempts and 48 boundary-ladder attempts.
No continuum target, physical constant or effective boundary operator enters
the search.

## Disclosed execution history

The first run stopped at `10/11` before every solve because its displaced
time-reversal anchors were outside an incorrectly defined branch gate. The
smaller-anchor correction also stopped at `10/11`; diagnosis showed that the
gate wrongly required every individual Lorentzian hinge curvature to be real.
The complete action and equations were real to `1e-77`. Both failures are
retained in commits `d2c018a` and `72c1749`.

After restoring the correct complex Lorentzian reality convention, a literal
50-decimal `3-point` run was interrupted after about 20 CPU-minutes without
finishing one six-seed class. No endpoint was recorded. A preregistered
complex128 evaluator then reproduced the high-precision action at 72 anchors
before it was allowed inside the unchanged solver.

The first complete run passed `12/12`. A final serialization defect was then
disclosed: wall-clock fields contradicted the byte-identical-artifact
requirement. Removing only those non-scientific fields and adding an exact
comparison with the first complete artifact produced two `13/13` runs with
identical SHA-256

```text
e945dc54a0768b00358aca6bef9e9a105ab3d0080d22dd83dfd140b038adf14d
```

No full suite was run.

## Evaluator and symmetry controls

Across all 72 pre-solve anchors, the largest fast/high-precision discrepancies
were

```text
relative action       2.165e-9
relative gradient     8.942e-11
minimum argument      7.855e-13
```

The largest binary64 angle-identity residual was `3.791e-12`, physical
imaginary contamination was `1.879e-11`, and fast time-reversal difference
was `2.799e-9`, all inside their preregistered gates. The independent
80-decimal time-reversal difference was `4.877e-75`.

The same solver exactly recovered its synthetic interior root and correctly
rejected a synthetic no-root control.

## DERIVED COMPUTATIONAL BOUNDED NEGATIVE

The hit counts are

```text
validated finite positive roots: 0
classes with a root:              0/12
eligible main endpoints:          0/72
accepted attempts:                0/120
```

The smallest main preconditioned residual was `1.2091e-5`, over two orders of
magnitude above the `1e-7` refinement gate. Its raw ten-equation norm was not
near zero; the minimum raw main norm over the census was `0.1355`. Optimizer
success flags occurred for 34/72 main attempts, illustrating why optimizer
status was not used as a scientific zero criterion.

The 48 ladder attempts likewise produced no eligible endpoint. Their best
preconditioned residual was `2.6470e-5`; only 9/48 optimizer flags were
successful.

The frozen outcome is

```text
REFINED_H4_NO_FINITE_ROOT_FOUND_OTHER
```

## Refuted simple interpretation

Zero of the 12 classes met the preregistered
`ZERO_LAPSE_BOUNDARY_PATTERN` criterion. The ladder endpoints did not
systematically activate their lower lapse bounds, and their residuals were
not monotone as the bound moved from `-4` to `-16`.

Therefore the local Newton proposal `y_rho approximately -2` was only a
**PATTERN**. The completed search refutes the stronger description that the
frozen trust-region trajectories simply track the zero-lapse boundary.

## Scope and honest verdict

**DERIVED COMPUTATIONAL BOUNDED NEGATIVE:** within the declared ten-dimensional
box, seeds, 12 schedule classes and solver, no finite stationary fill was
found. The inherited static background cannot be repaired by this complete
bounded search and cannot license an effective boundary Hessian.

**OPEN:** global nonexistence. The search is not interval Newton, Krawczyk,
degree theory or an exhaustive algebraic certificate. Its many interior
nonzero endpoints may indicate poor conditioning or merit-function stationary
points rather than a mathematical obstruction. No claim stronger than the
frozen search space could have been falsified by this verifier.

This result has reproducibility and a mechanically distinct fast/high-
precision equation cross-check, but it does not yet have an independent
root-exclusion method. Under Rule 4 the global physical conclusion remains
unaccepted until the reduction below attacks it differently.

## Next adversarial test

At a fixed common lapse coordinate `t`, solve only the six cross-diagonal
equations. Then evaluate the four remaining rank-lapse equations on that
six-equation solution branch. This is a `6+4` Lyapunov--Schmidt/Schur-style
reduction, not another ten-dimensional multistart search.

It must report:

- existence, uniqueness and branch validity of the six-equation solution;
- the four reduced residuals and their independent combinations;
- signs, zeros and rank over a preregistered continuous lapse interval;
- all 12 schedule classes and time reversal;
- whether the reduced residual has a finite zero or only a degenerate limit.

Only a certified finite common zero advances the static route. A sign or
interval exclusion can close it mathematically. An unresolved reduced branch
leaves the route **OPEN** and motivates moving to unequal boundaries, not
fitting dust weights.

Nothing here derives a physical tick, `c`, `G`, Planck scales, graviton or
particle mass.
