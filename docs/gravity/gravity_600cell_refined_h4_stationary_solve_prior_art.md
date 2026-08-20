# Prior-art gate: local well-posedness of the refined H4 stationary solve

Date: 2026-08-20

Status: completed before constructing any ten-by-ten internal Jacobian.

## Exact question and complete hypotheses

On `K0=P(sd K_600)`, keep both equal spatial boundaries, the selected total
dust mass and the conditional local `P1` weights fixed.  For each of the 24
rank-ordered staircase slabs, use all six cross-diagonal and four positive
vertical lapse-square coordinates.  Evaluate the derivative of the ten
internal Regge equations with respect to their ten logarithmic coordinates at
the already certified induced static fill.

The immediate question is only whether the local stationary equations are
full rank, gauge-degenerate or schedule dependent.  No Newton step is applied
in this gate, and the point is already known to be off shell.

## KNOWN from primary literature

The source audit in
`gravity_600cell_refined_h4_stationary_fill_prior_art.md` applies unchanged:

- an action-generated canonical map requires the internal Regge equations,
  and its Hessian/constraints can have null directions associated with
  discrete gauge symmetry: Dittrich and Hoehn,
  [From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817)
  and [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974);
- exact vertex-displacement symmetry is special to flat linearized
  backgrounds, while nonlinear curved discretizations generically produce
  pseudo-constraints rather than an automatically free lapse;
- four-dimensional Regge Hessians need not be triangulation independent:
  Dittrich and Steinhaus,
  [Path integral measure and triangulation independence in discrete gravity](https://arxiv.org/abs/1110.6866),
  and Dittrich, Kaminski and Steinhaus,
  [Discretization independence implies non-locality in 4D discrete quantum gravity](https://arxiv.org/abs/1404.5288).

No cited theorem fixes the rank or equality of the present 24 matrices.

## Framing attack

The Jacobian of the internal gradient is an internal block of the action
Hessian.  At an off-shell point it is a legitimate Newton-conditioning
diagnostic, but it is **not** an effective boundary Hessian, constraint
algebra, graviton operator or physical spectrum.

A small singular value cannot be called gauge by threshold alone.  It must be
separated from the numerical error envelope and its vector must be compared
with the exact tangent of the induced common-lapse family.  Conversely, full
rank at the off-shell point does not prove that a nearby root exists or that
the root is unique globally.

Failure of a local Newton search cannot prove nonexistence.  Only a certified
root is positive evidence; a failed bounded solver remains **OPEN** unless a
separate exclusion certificate is constructed.

## Repository consequence

The old collective lapse cancellation does not establish a Hessian null
direction.  Along the induced family, cross diagonals and vertical edges move
together and the off-shell gradient itself can vary.  The rank must therefore
be measured before imposing a clock gauge or dropping an equation.

External novelty of this exact all-24 matrix census is **OPEN**.  The
variational and triangulation warnings are **KNOWN**.

## Next admissible calculation

Preregister a high-precision all-schedule internal Jacobian census with
step-doubling, Hessian-symmetry, directional-second-derivative and
time-reversal controls.  Its frozen outcome alone selects whether the later
stationary solve is square, gauge-fixed or unresolved.
