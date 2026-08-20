# Prior-art gate: positive stationary roots of the refined H4 slab

Date: 2026-08-20

Status: completed after the internal-Jacobian census and before evaluating
the ten equations away from the inherited fill.

## Exact object and complete hypotheses

For `K0=P(sd K_600)`, use the exact rank-derived spatial geometry, both equal
fixed boundaries, the corrected complex Lorentzian Regge action with boundary
terms, the previously selected total dust mass, the conditional `P1` weights
and one representative of each of the 12 certified schedule/time-reversal
classes.

For each representative, solve the ten internal equations in the logarithms
of six positive cross-diagonal squared lengths and four positive lapse
squares. The question is whether a finite, nondegenerate root exists on the
same Lorentzian angle branch as the inherited fill. The supplied
`tau0=0.0102` is an initial scale only; finding a root near it would not derive
an absolute unit of time.

No boundary evolution, effective Hessian, mode spectrum, continuum target or
physical constant is part of this mission.

## KNOWN from primary literature

- Internal Regge equations must be imposed before the slab action can be
  treated as Hamilton's principal function for boundary evolution. The
  discrete equations can fix data that would be lapse/shift gauge in the
  continuum: Dittrich and Hoehn,
  [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974).
- Curved Regge discretizations generically replace exact continuum gauge
  constraints by background-dependent pseudo-constraints: Bahr and Dittrich,
  [(Broken) Gauge Symmetries and Constraints in Regge Calculus](https://arxiv.org/abs/0905.1670),
  and Dittrich and Hoehn,
  [From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817).
- Newton damping, trust regions and residual merit functions are standard
  numerical globalization devices. Inexact Newton convergence does not turn
  solver failure into root nonexistence: Dembo, Eisenstat and Steihaug,
  [Inexact Newton Methods](https://doi.org/10.1137/0719025), and Pawlowski et
  al., [Inexact Newton Dogleg Methods](https://doi.org/10.1137/050632166).
- A rigorous exclusion over a specified box requires interval Newton,
  Krawczyk, degree or comparable exhaustive machinery. A finite multistart
  search supplies no such theorem.

The literature therefore predicts neither a finite root nor a zero-lapse
endpoint for this carrier. A focused post-Jacobian search using the terms
`full-rank Regge internal Hessian`, `zero lapse`, `tent move`,
`pseudo-constraint` and `stationary simplicial slab` located no primary source
that computes this ten-variable projected 600-cell system. Search absence is
not novelty evidence; external novelty remains **OPEN**.

## CONTROL from the repository

- The inherited fill is off shell in exactly the four rank-lapse equations.
- Its ten-by-ten Jacobian has certified rank ten and inertia `(8,0,2)` for all
  schedules.
- The 24 schedules form 12 matrix classes paired by time reversal.
- The initial residual is common to all schedules.
- The unapplied local Newton correction sends all four log lapse-square
  variables approximately to `-2` and moves cross variables by at most
  `0.0112`.

## Framing attack

A positive numerical endpoint is not automatically a physical clock. It is
first a stationary internal fill of one fixed triangulation with an inherited
dimensionful initialization. It advances the construction only if it is
branch-valid, independently validated and compatible across all schedule
classes.

Conversely, convergence toward lapse zero is not by itself a no-root theorem.
Log coordinates never reach zero at a finite iterate, so an optimizer can
make the residual small merely by suppressing the complete action with the
shrinking time thickness. Such endpoints must be classified as degenerate,
not as roots.

The action has an indefinite local Hessian. Least-squares minimization can
therefore terminate at a nonzero stationary point of the merit function.
Optimizer success flags are never scientific acceptance gates; the original
ten equations and branch data decide.

## OPEN

- existence of a finite positive common-branch root in any class;
- equality of roots between different schedule classes;
- whether every accepted trajectory instead approaches lapse zero;
- global nonexistence outside a frozen bounded search domain;
- any relation between a finite lapse and a derived physical tick.

## Next admissible calculation

Preregister a deterministic bounded multistart solve with frozen coordinates,
seeds, residual scaling, trust-region settings, positive/branch gates and an
independent high-precision final residual. Report search-space coverage and
distinguish `ROOT`, `DEGENERATE_LIMIT`, `UNRESOLVED` and control failure.

Only an interval or other exhaustive certificate may later upgrade a bounded
negative search to a mathematical no-root result.
