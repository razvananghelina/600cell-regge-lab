# Prior-art gate: all-index weak-lapse recurrence of the cellular 600-cell map

Date: 2026-08-17

Status: written after the independently committed four-step derivation
`76a09ab` and comparison `ae5c3bc`, before constructing or running an
all-index recurrence verifier.

## 1. Exact object and hypotheses

Use the same fixed homogeneous cellular Regge action as in
`gravity_600cell_cellular_weak_lapse_protocol.md`, with:

- one regular 600-cell `S^3` carrier on every boundary;
- one conserved total dust mass
  `M/L0=(90/pi)[2*pi-5*acos(1/3)]`;
- positive Lorentzian lapse variable `rho_n`;
- homothetic endpoints only;
- the lapse equation `F_n=0` for each slab;
- the additive canonical seam equation `G_n=0` at each shared boundary;
- initial static slab `L_-1=L_0=L0`, `rho_0=e^2 L0^2`;
- the real contracting branch continuous at `e=0`.

Write, with `x=e^2`,

```text
log(L_n/L0)       = A_n x+B_n x^2+O(x^3),
log(rho_n/(e^2L0^2)) = R_n x+O(x^2).
```

The new question is whether the coefficient equations imply a closed
recurrence and solution for every formal integer `n>=1`, not merely the four
indices already evaluated.

## 2. Primary prior art

- Collins and Williams introduced regular 5-, 16- and 600-cell Regge models
  of closed Friedmann dynamics:
  <https://doi.org/10.1103/PhysRevD.7.965>.
- De Felice and Fabri evolved dust on the 600-cell with a Sorkin algorithm:
  <https://arxiv.org/abs/gr-qc/0009093> and
  <https://arxiv.org/abs/gr-qc/0106077>.
- Gentle and Miller describe an implicit local Regge evolution scheme and a
  preliminary 600-cell Friedmann application:
  <https://arxiv.org/abs/gr-qc/9411008>.
- Dittrich and Hoehn formulate additive simplex actions as generating
  functions for canonical discrete evolution:
  <https://arxiv.org/abs/1108.1974>.
- Liu and Williams derive global and local Collins--Williams equations and
  analyze their continuum approximation and resolution dependence:
  <https://doi.org/10.1103/PhysRevD.93.024032>.
- Tsuda and Fujiwara give nonlinear regular 4-polytopal frustum recurrence
  equations and their continuum-time Friedmann limit:
  <https://arxiv.org/abs/2011.04120>.

These sources establish the broad recurrence mechanism.  Searches for
`600-cell dust weak lapse expansion`, `Collins-Williams 600-cell time
symmetry recurrence`, and `weak lapse Regge cosmology` did not locate the
specific coefficient formulas below.  A negative search is not evidence of
novelty.

## 3. KNOWN / CONTROL / OPEN

### KNOWN

- The closed cellular action is exactly subdivision independent on the
  homogeneous family.
- The first four coefficient systems have leading rank one and next rank
  two, with determinant `16200 epsilon_3^2`.
- For `n<=4`, exact action differentiation produced
  triangular `A_n`, square `R_n`, and odd post-momentum coefficients.
- The resulting jet explains all committed finite-lapse four-tick data.
- Quadratic turning-point motion is the expected local form of closed
  Friedmann dynamics.

### CONTROL

- Derive the coefficient equations with symbolic generic boundary data;
  do not infer a recurrence by interpolating the four stored outputs.
- Require an exact algebraic induction step, including the nuisance
  coefficient `B_n` needed to determine `R_n`.
- Require the rank and determinant claims to hold symbolically for generic
  integer index `n`.
- Verify the closed form at additional indices not used to discover it.
- Keep the formal all-index statement separate from uniform convergence of
  the asymptotic series as `n` grows.

### OPEN

- Whether the four-step formulas obey the same recurrence for all `n`.
- A derived closed form for `B_n`.
- The radius of validity in the combined variables `n` and `e`.
- Spatial-refinement convergence and anisotropic dynamics.
- External novelty of the exact coefficient recurrence.

## 4. Proposed difference and evidential value

The narrow possible addition is an exact induction theorem for the
weak-lapse coefficient jet of this particular conserved-dust cellular action.
Even if proved, it is an internal mathematical closure of a known Regge
minisuperspace mechanism, not a new gravitational theory.

The physically important next test remains spatial refinement of the
acceleration coefficient.  An all-index theorem only prevents us from
mistaking a four-term coincidence for the local recurrence.

## 5. Procedural disclosure

Before this note was committed, an exploratory interpolation was applied to
the already committed `B_-1,...,B_4` values.  It suggested a degree-four
polynomial.  This occurred after the literature search but before the formal
prior-art-gate commit.  It is labelled **PATTERN** and is inadmissible as
evidence for the theorem.  The preregistered verifier must derive the generic
coefficient equations from the frozen action and must be able to refute that
interpolant.

**External novelty: OPEN.**

