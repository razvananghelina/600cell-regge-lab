# Prior-art gate: H4-invariant stationary fill of the refined slab

Date: 2026-08-20

Status: completed before evaluating any refined internal Regge residual.

## 1. Exact question and complete hypotheses

Use the fixed projected barycentric carrier

```text
K0 = P(sd K_600),   f=(2640,17040,28800,14400),
```

its rank colouring, all 24 standard colour-ordered staircase triangulations
of `K0 x I`, the corrected Lorentzian Regge action with boundary terms, and
the already selected static scale, supplied lapse `tau0=0.0102`, total dust
mass and conditional local `P1` dust distribution.

Restrict only by the full spatial `H4` symmetry retained by every schedule.
This gives six old spatial edge types, six new spatial edge types, six
cross-diagonal types and four vertical-edge types.  The question is whether
the induced flat static fill is stationary in **each** of the ten internal
types, not merely under their previously tested common homothetic variation.

This gate does not solve for a new fill, compute a Hessian or spectrum, choose
a schedule, compare with a continuum mode or infer a physical tick.

## 2. KNOWN from primary sources

- In canonical simplicial gravity, the action's Hamilton principal function
  generates the boundary evolution only after the bulk equations associated
  with internal variables are imposed.  Free data can be fixed by later
  constraints rather than being gauge automatically: Dittrich and Hoehn,
  [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974).
- Discrete constraints are exact on flat linearized backgrounds but become
  background-dependent pseudo-constraints beyond that regime: Dittrich and
  Hoehn,
  [From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817).
- Four-dimensional Regge theory is not generically triangulation independent.
  Its linearized Hessians under Pachner moves have special structure but do
  not supply a general equality theorem for the present curved closed slab:
  Dittrich and Steinhaus,
  [Path integral measure and triangulation independence in discrete gravity](https://arxiv.org/abs/1110.6866),
  and Dittrich, Kaminski and Steinhaus,
  [Discretization independence implies non-locality in 4D discrete quantum gravity](https://arxiv.org/abs/1404.5288).
- Propagating curvature degrees of freedom can be identified in canonical
  linearized Regge calculus, but only after separating vertex-displacement
  gauge directions and imposing the nontrivial equations of motion: Hoehn,
  [Canonical linearized Regge Calculus: counting lattice gravitons with
  Pachner moves](https://arxiv.org/abs/1411.5672).

These results make internal stationarity a necessary gate.  None proves that
the 24 refined 600-cell fills agree or that the inherited static fill solves
the enlarged equations.

## 3. Repository distinction

The accepted projected acceleration used a direct **homothetic cellular**
action with only `s_minus`, `s_plus` and a common lapse square.  Its static
lapse equation therefore certifies only one collective internal variation.
The conditional `P1` construction distributes the selected mass locally and
equivariantly, but it explicitly did not construct independent vertex-lapse
dynamics.

Consequently neither accepted result implies ten separate internal equations.
Calling the static fill a full refined solution before this census would be a
minisuperspace-to-full-system inference and is forbidden.

## 4. Framing attack

The full `H4` restriction is a **necessary falsification sector**, not the
physical phase space.  Passing it cannot establish schedule independence or
gravitons.  Failing it is nevertheless decisive against using the inherited
fill as the background of an effective Hessian: a Schur complement evaluated
away from the internal equations is not the action-generated canonical map.

The ten-variable bound is exact only on this invariant sector.  It says
nothing about the many non-invariant internal variables, which remain for a
later test even if this gate passes.

## 5. Novelty status

The variational principle and triangulation-dependence warning are **KNOWN**.
The exact all-24 `H4` residual census on the projected barycentric 600-cell is
not located in the cited primary literature, but a targeted search cannot
prove absence.  External novelty is **OPEN**.

## 6. Next admissible calculation

Preregister a target-blind, high-precision computation of all ten internal
log-length residuals for every schedule, with an independent finite-difference
control and an exact rank-type reconstruction of the projected chamber
geometry.  Only after its frozen outcome may a stationary-fill solve or an
effective Hessian be attempted.
