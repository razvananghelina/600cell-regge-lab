# Second refinement improves degree balance but refutes an exact repeated factor

Date: 2026-08-11

Original preregistration: `702fa5b`  
Recorded matrix-free solver failure: `6f21619`  
Explicit quotient correction frozen before rerun: `e3e4294`

Targeted verifier:
`reproducible/verify_whitney_trace_stiffness_second_refinement.py`

Corrected targeted result: **13/13 PASS**.  The verifier is registered.  The
full suite was not run by explicit user request.

## Headline

The first-step stiffness factors do not repeat exactly under the second
barycentric refinement:

\[
 R^{01}=(11.6371,,6.81065,,3.80924),
\]

\[
 R^{12}=(7.27607,,5.41028,,5.44882).
\]

Thus:

> **DERIVED NUMERICAL NEGATIVE:** there is no exact repeated scalar
> refinement factor on this three-level control.

However, the three form degrees align much better at the second step:

\[
 \operatorname{spread}(R^{01})=3.05498,
 \qquad
 \operatorname{spread}(R^{12})=1.34486.
\]

Therefore:

> **PATTERN:** the exact trace energy flows toward a common degree scale over
> the two tested steps.

This is not yet a fixed point, exponent, or physical renormalization law.

## Complete three-level geometry

Every cell is retained:

| level | f-vector | tetrahedra | duplicated dimension | element metric types | face metric types |
|---:|---|---:|---:|---:|---:|
| 0 | `(5,10,10,5)` | 5 | 75 | 1 | 1 |
| 1 | `(30,150,240,120)` | 120 | 1,800 | 1 | 4 |
| 2 | `(540,3420,5760,2880)` | 2,880 | 43,200 | 24 | 57 |

The level-two tetrahedra are not replaced by regular or averaged elements.
All 24 exact metric types and all 57 exact face types are used in their actual
locations.

Every shared face has the same exact induced metric from both parents.  All
occurrence graphs are connected, and the exact constraint ranks are:

\[
 (15,20,10),
 \quad
 (450,570,240),
 \quad
 (10980,13860,5760).
\]

Thus every quotient kernel is precisely conformity at all three levels.

## Solver certification and recorded failure

The first preregistered matrix-free run returned `10/11`.  Its level-two
degree-zero and degree-one extremal residuals were respectively
(4.32\times10^{-4}) and (8.39\times10^{-4}), so its provisional ratios
were rejected and preserved in commit `6f21619`.

Commit `e3e4294` froze an explicit sparse construction of

\[
 A=V^*HRM^{-1}R^*HV,
 \qquad
 B=V^*HV
\]

and generalized symmetric Lanczos, without changing the operator or the
(10^{-7}) Ritz gate.

The corrected calculation gives:

- dense calibration maximum relative error:
  (5.33\times10^{-15});
- dense calibration maximum residual:
  (1.51\times10^{-12});
- maximum residual over all 18 spectral edges:
  (1.01\times10^{-11});
- maximum pre-symmetrization sparse transpose residual:
  (1.10\times10^{-13});
- accepted level-one matrix-free cross-check error:
  (3.78\times10^{-15}).

## Positive gaps and local Dirac norms

| level | worst local Dirac norm (a_r) | (g_{r,0}) | (g_{r,1}) | (g_{r,2}) |
|---:|---:|---:|---:|---:|
| 0 | 3.872983 | 8.660254 | 5.196152 | 3.464102 |
| 1 | 11.936870 | 2.293663 | 2.351464 | 2.802834 |
| 2 | 63.308026 | 1.671865 | 2.305084 | 2.728120 |

All nine quotient gaps remain positive.  Finite trace stiffness therefore
survives the second finite refinement algebraically.

The corresponding scales (a_r/g_{r,p}) are:

| level | degree 0 | degree 1 | degree 2 |
|---:|---:|---:|---:|
| 0 | 0.447214 | 0.745356 | 1.118034 |
| 1 | 5.204283 | 5.076356 | 4.258856 |
| 2 | 37.866715 | 27.464519 | 23.205736 |

## What the two steps say

The componentwise changes of the refinement factors are

\[
 \frac{R^{12}}{R^{01}}
 =
 (0.625246,,0.794386,,1.430423).
\]

They are far from one, so exact repetition is decisively absent.  Yet degrees
one and two have almost coincident second-step factors:

\[
 5.41028\quad\text{and}\quad5.44882.
\]

Degree zero remains the bottleneck at (7.27607), but it moved toward the
other two relative to the first step.  This is the precise content of the
positive pattern.

## Hostile interpretation: mesh degeneration may be the signal

The worst local Dirac norm grows by factors

\[
 3.08\quad\text{then}\quad5.30,
\]

while the number of ordered element metric types grows

\[
 1\to1\to24.
\]

Thus the apparent flow is not purely a change of physical resolution.  It
also includes worsening and proliferating simplex shapes under iterated
barycentric subdivision.  The level-two element norms range from about
32.31 to 63.31.

This attacks the encouraging interpretation directly: a clustered stiffness
factor could be compensating mesh anisotropy rather than revealing a physical
renormalization group.

Consequently a brute-force third barycentric step is not yet the best next
experiment.  First one must measure shape quality and compare against a
shape-regular refinement rule.  Otherwise increasing factors may merely
track degenerating numerical elements.

## Comparison with the complete 600-cell first step

The complete 600-cell first-step ratios were

\[
 (4.3793,3.2485,3.8092),
\]

whereas the control's first and second steps differ substantially.  Therefore
there is no evidence yet for a topology-independent universal scalar.

## Physical verdict

The exact face metric remains the correct local construction found so far:
it preserves positive quotient gaps and increasingly aligns the form degrees.
But the fourth ingredient needed for physics is still missing: a
shape-regular continuum/refinement prescription selected by the theory.

Nothing here derives the absolute (kappa), a finite causal speed, Lorentzian
time, mass, or a Planck scale.

## Status ledger

- **DERIVED:** complete three-level carriers and exact metric-type census.
- **DERIVED:** exact face agreement, kernels, ranks, and positive gaps.
- **DERIVED NUMERICAL:** all 18 quotient spectral edges pass.
- **DERIVED NUMERICAL NEGATIVE:** exact refinement factor does not repeat.
- **PATTERN:** degree spread improves from `3.05498` to `1.34486`.
- **PATTERN:** degree-one and degree-two second-step factors nearly coincide.
- **STRUCTURAL WARNING:** worst-element Dirac growth accelerates and element
  types proliferate.
- **OPEN:** shape-quality audit and shape-regular refinement comparison.
- **OPEN:** complete second-refined 600-cell.
- **OPEN:** asymptotic law, absolute stiffness, chirality, and causal time.
- **NOT CLAIMED:** mass, inertia, (c), (hbar), Newton's (G), or a Planck
  scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_trace_stiffness_second_refinement.py
```

Expected corrected result: `13/13`.
