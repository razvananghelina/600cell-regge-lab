# Prior-art gate: projection structure of compatible 600-cell boundary data

Date: 2026-08-19

## Exact object and hypotheses

Let `F`, `E`, and `S` be the rational sparse blocks of the complete
variable-face flat-frustum compatibility matrix already frozen by
`verify_gravity_600cell_canonical_data_admissibility.py`:

```text
F: 3600 cell-flex columns
E:  720 upper-spatial-edge squared-length columns
S:  120 strut squared-length columns
```

The lower regular 600-cell is fixed.  The two nonstatic rational
representatives are `(lambda,tau)=(2,5)` and `(3,11)`.  The immediate object
is the finite-field compatible-data space

\[
K_p=\{(e,s)\in\mathbb F_p^{720}\oplus\mathbb F_p^{120}:
\exists z\in\mathbb F_p^{3600},\;Fz+Ee+Ss=0\}.
\]

This mission asks only for the dimensions of `K_p`, its intersections with
`s=0` and `e=0`, and its projections onto edge and strut data.  It does not
yet construct an action Hessian, a symplectic form, a physical constraint
surface, propagating modes, a tick, or a speed.

## KNOWN

Dittrich and Hoehn derive canonical dynamics from a discrete action and show
that linearized Regge calculus on flat backgrounds has vertex-displacement
constraints; nonlinear corrections can turn them into pseudo-constraints.
See [arXiv:0912.1817](https://arxiv.org/abs/0912.1817), especially the flat
linearized analysis and conclusion.

Their general simplicial formalism allows moves to introduce free data that
may later be fixed by constraints.  This is the correct general warning
against identifying raw boundary coordinates with freely admissible Cauchy
data.  See [arXiv:1108.1974](https://arxiv.org/abs/1108.1974), abstract and
the pre/post-constraint construction.

The pre/post-constraint and reduced-phase-space framework for variational
discrete systems is developed in
[arXiv:1303.4294](https://arxiv.org/abs/1303.4294).  It also stresses that
propagating observables can depend on the chosen initial and final moves.

For linearized four-dimensional Regge calculus, Hoehn identifies
vertex-displacement gauge variables and gauge-invariant lattice-graviton
degrees of freedom under Pachner moves in
[arXiv:1411.5672](https://arxiv.org/abs/1411.5672).  Therefore a future
interpretation of any nonhomogeneous modes must separate gauge directions
from curvature modes rather than infer gravitons from dimension alone.

The use of the 600-cell in Regge evolution is established prior art.  A
generalized 600-cell evolution with more free variables and a
causality-breaking endpoint was studied by De Felice and Fabri in
[arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077).  Homogeneous
600-cell dust evolution and its stopping point were studied in
[arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093).  These works mean
that neither “Regge + 600-cell” nor a finite admissibility obstruction is by
itself new.

## CONTROL

Pure linear algebra gives the projection census without constructing a basis.
Because the frozen calculation has `rank(F)=3600`, `F` has no kernel.  Hence,
for each prime,

\[
\begin{aligned}
\dim K_p &=4440-\operatorname{rank}[F\ E\ S],\\
\dim(K_p\cap\{s=0\})&=4320-\operatorname{rank}[F\ E],\\
\dim(K_p\cap\{e=0\})&=3720-\operatorname{rank}[F\ S],\\
\dim\pi_E(K_p)&=\dim K_p-\dim(K_p\cap\{e=0\}),\\
\dim\pi_S(K_p)&=\dim K_p-\dim(K_p\cap\{s=0\}).
\end{aligned}
\]

These are quotient/image identities, not a physical interpretation.  The two
homogeneous rational scale/lapse tangents are known controls and must remain
compatible.

## OPEN and proposed difference

The search above did not identify a primary source computing this exact
`F/E/S` projection census for the schedule-free, complete variable-face,
flat-frustum 600-cell carrier.  A web search cannot prove absence; external
novelty remains **OPEN**.

The proposed difference is narrow: a target-blind, convention-attacked census
of how the already frozen 240 modular modes project onto spatial-edge and
strut data.  It is not a claim of new gravity or of physical degrees of
freedom.  Any rational equality, geometric carrier, gauge interpretation, or
propagation claim requires a later exact construction and mechanically
independent replication.
