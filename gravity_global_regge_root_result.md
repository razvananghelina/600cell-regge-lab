# No stationary root from the frozen regular-boundary searches

Date: 2026-08-13

Root-search protocol commit: `5b687a3`

Disclosed correction commits: `82ab7b8`, `6594398`

Certified evaluator commit: `d9fe159`

Registered verifier:
`reproducible/verify_gravity_global_regge_roots.py`

Machine-readable search history:
`reproducible/gravity_global_regge_roots.json`

Final targeted run: **16/16 passed**.  The full repository suite was not run.

## Headline

> **DERIVED SEARCH NEGATIVE, STRICTLY SCOPED:** with both 600-cell boundary
> metrics fixed regular, no stationary point was found from any of the six
> preregistered starts in either of the two 35-variable schedule-invariant
> Lorentzian sectors.

This is not a no-root theorem for either invariant box, and it says nothing
conclusive about the full 840-variable space.  Seven trajectories reach the
artificial lower variable bound, one reaches a causal/angle-branch boundary,
and four exhaust the fixed 80 accepted iterations.  None has terminal
residual remotely close to the `1e-10` acceptance scale.

The result narrows the next question: a global vacuum tick should not be sought
as a return between two independently fixed identical spatial metrics.  The
final boundary metric and its canonical momentum must enter as boundary data
of an evolution problem.

## 1. Fixed problem

For each ordered phase parity, the search used exactly:

- 30 independent staircase-diagonal squared-length orbits;
- five independent positive pole-magnitude orbits;
- old and final 600-cell boundary squared edge lengths fixed to one;
- the certified zero-volume Lorentzian Regge action;
- the corrected complex-angle branch;
- all 35 per-edge internal equations, without component weighting;
- logarithmic variables with the frozen box `exp(-6)<x_i<exp(6)`.

The same six deterministic starts were used for both parities.  No random
restart or start chosen after seeing a residual was added.

The Levenberg--Marquardt rule, damping list, backtracking factors, causal gates,
Jacobian step and 80-iteration limit were all committed before the first
iteration.

## 2. Complete outcome table

The residual below is the Euclidean norm of the 35 per-edge gradients after
dividing out the common orbit size 24.

| parity | start | termination | accepted steps | terminal `||r||_2` | smallest terminal variable |
|---|---:|---|---:|---:|---:|
| even | S0 | artificial-box contact | 25 | 2.50596515 | pole `0.0024787540` |
| even | S1 | 80-step limit | 80 | 2.72121041 | diagonal `0.01591335` |
| even | S2 | artificial-box contact | 17 | 2.50579274 | diagonal `0.0024787552` |
| even | S3 | 80-step limit | 80 | 2.71601642 | diagonal `0.02428634` |
| even | S4 | causal/branch Jacobian failure | 64 | 2.57311420 | pole `0.08142793` |
| even | S5 | artificial-box contact | 24 | 2.60343554 | diagonal `0.0024787704` |
| odd | S0 | artificial-box contact | 14 | 2.68321548 | diagonal `0.0024787647` |
| odd | S1 | artificial-box contact | 18 | 2.85027779 | diagonal `0.0024787977` |
| odd | S2 | artificial-box contact | 19 | 2.58959499 | diagonal `0.0024787925` |
| odd | S3 | 80-step limit | 80 | 2.85780221 | diagonal `0.05864107` |
| odd | S4 | 80-step limit | 80 | 2.60563937 | diagonal `0.01691197` |
| odd | S5 | artificial-box contact | 26 | 2.59247754 | diagonal `0.0024787996` |

For reference,

```text
exp(-6) = 0.00247875217667...
```

so the seven box-contact values are mechanically identifiable, not subjective
labels.  The exceptional `even/S4` stays far from that box boundary but has
minimum angle-argument modulus `0.00149484`; its centered Jacobian attempts to
leave the admitted causal/branch region.

The initial norms range from `3.09219` to `3.21861`.  The terminal range is

```text
2.50579274 <= ||r||_2 <= 2.85780221.
```

Thus none is a near-miss against the required `3e-10` norm.

## 3. Reproducibility and corrections

The first sequential run returned `14/14`, zero roots and eight generically
named Jacobian failures.  The solver was then parallelized only across the 70
independent centered points used to build each Jacobian; candidate selection
remained serial.

The first parallel rerun reproduced all twelve terminal residuals exactly but
returned `15/16`, because a correction note had incorrectly classified all
eight Jacobian failures as artificial-box contacts.  Inspection showed the
`even/S4` causal/branch case above.  That correction was committed before the
final rerun.

The final targeted run reproduces every disclosed terminal residual with

```text
maximum relative difference = 0.
```

and obtains the exact status census

```text
even: 3 box contacts + 1 branch failure + 2 iteration limits
odd:  4 box contacts +                    2 iteration limits.
```

This provenance is intentionally retained rather than presenting the final
classification as if it had been obvious initially.

## 4. What the negative means

- **DERIVED:** all twelve frozen starts are initially Lorentzian and lie on
  the certified real action branch.
- **DERIVED SEARCH NEGATIVE:** none yields a validated stationary root under
  the frozen algorithm.
- **OPEN:** roots elsewhere in the same 35-dimensional causal box.
- **OPEN:** roots outside the artificial box.
- **OPEN:** roots without the order-24 symmetry, in all 840 variables.
- **OPEN:** roots for nonregular final boundary metrics.
- **NOT CLAIMED:** nonexistence of global Lorentzian evolution.

In particular, the seven artificial-box contacts are inconclusive by the
preregistered rule.  They cannot be cited as evidence that a minimum occurs at
zero length; the search was stopped before allowing degeneration.

## 5. Attack on the framing

Fixing both boundary metrics to the same regular 600-cell asks for a special
two-boundary return solution.  A dynamical law normally does not select both
initial and final configuration data.  In the canonical tent-move picture,
the action is a generating function: boundary momenta accompany the boundary
lengths, while only internal-edge variations must vanish.

Therefore failure of the equal-boundary root search does not mean that the
global carrier lacks dynamics.  It means the strongest static-return ansatz
has no support from the frozen local searches.

This interpretation is consistent with the earlier local result: the
one-orbit regular Lorentzian tent has no stationary pole, while an asymmetric
final star can have one.  The global calculation now says that merely making
the 35 internal orbits asymmetric, while keeping both complete boundaries
regular, did not recover a nearby root.

## 6. Status ledger

| Claim | Status |
|---|---|
| The 35-variable root evaluator is the certified full-action restriction | **DERIVED COMPUTATIONAL** upstream |
| A regular equal-boundary stationary slab exists | **DERIVED NEGATIVE** at the exact regular point |
| One of six nearby starts finds an equal-boundary root | **DERIVED SEARCH NEGATIVE** |
| The two phase parities are distinguished by root existence | **DERIVED NEGATIVE for this search:** neither produced one |
| No equal-boundary root exists in the whole causal box | **OPEN** |
| No invariant global vacuum evolution exists | **OPEN** |
| No unconstrained 840-variable evolution exists | **OPEN** |
| Final boundary data should be fixed regular | **REFUTED AS A GENERAL DYNAMICAL REQUIREMENT** |

## 7. Next correct calculation

Keep the initial regular 600-cell fixed, but allow the final spatial boundary
edges to vary in all orbits permitted by the same ordered-schedule stabilizer.
Do not set their action derivatives to zero: record them as post-momenta.  Then
ask whether the 35 internal equations define a canonical relation between old
and new boundary phase space.

The target-independent order is:

1. enumerate final-boundary edge orbits and the full old/new momentum orbits;
2. derive the rectangular internal Jacobian with respect to internal and final
   boundary variables;
3. determine its exact/numerical rank before choosing boundary perturbations;
4. preregister a continuation problem based on that rank;
5. accept dynamics only when all 840 individual internal gradients vanish and
   both boundary momentum vectors are real;
6. only afterward ask whether either phase parity is selected.

This tests evolution rather than an overconstrained static return.
