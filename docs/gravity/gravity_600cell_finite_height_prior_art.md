# Prior-art gate: exact finite-height fixed-state census

Date: 2026-08-21

Status: written after the adversarially corroborated cubic obstruction and
before evaluating the finite-height elimination function or any positive
root.

## 1. Exact object and hypotheses

Use the same certified homogeneous cellular 600-cell
Regge-plus-conserved-dust action, positive Lorentzian branch and zero
cosmological constant.  Normalize `L_minus=1`.  For every real state
parameter

```text
v in R,
```

freeze exactly

```text
M=mu(v),
p0=p(v).
```

Do not impose the failed Taylor endpoint.  Solve the complete unexpanded
equations

```text
F(1,L_plus,h^2;mu(v))=0,
p_pre(1,L_plus,h^2;mu(v))=p(v)
```

for

```text
h>0,
L_plus>0.
```

The question is whether the same incoming canonical state selects no
positive solution, an isolated positive solution, or a positive-dimensional
family.

Unlike the preceding Taylor census, the exact finite problem must include
`v=0` and `K(v^2)=0`.  The state functions `mu(v)` and `p(v)` are regular
there; only the lower-order division by `v*K` failed.  Excluding those strata
would make the finite census incapable of falsifying an isolated-root claim.

## 2. Framing correction: a bounded box is not enough

There is no geometric upper bound on `h` in the present action, and endpoint
positivity alone does not compactify `(L_plus,h)`.  Therefore a numerical
search in an arbitrary rectangle cannot prove absence and cannot close the
tick route.

Introduce instead the exact dimensionless endpoint quotient

```text
q=(L_plus-1)/h,
L_plus=1+h q.
```

For real `q`, all square-root and inverse-function arguments in the positive
Lorentzian action remain real.  Physical endpoint positivity becomes

```text
1+h q>0.
```

After differentiating the complete action first and substituting
`rho=h^2`, both exact residuals are affine in `h`:

```text
C(v;q,h)=C0(v;q)+h Ch(q),
P(v;q,h)=P0(v;q)-p(v)+h Ph(q).
```

This affineness follows structurally from `L_plus=1+hq` and
`(1-L_plus^2)/h=-2q-hq^2`; no finite root has yet been evaluated.

The complete problem can therefore be reduced without a height cutoff to
the one-variable elimination equation

```text
D(v,q)=C0(v;q) Ph(q)-[P0(v;q)-p(v)] Ch(q)=0,
```

together with exact treatment of all zero-slope cases and the reconstructed
height

```text
h=-C0/Ch=-[P0-p(v)]/Ph>0,
1+hq>0.
```

A later protocol must classify `D(v,q)` on the full real `q` line and the
complete registered `v` domain.  A grid is only a diagnostic.

## 3. Primary prior art

- [De Felice--Fabri, arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093)
  evolve a dust-filled 600-cell with a Sorkin algorithm.  They report that
  earlier authors also found a stopping point and interpret their stop as a
  causality-breaking spacetime singularity.
- [De Felice--Fabri, arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077)
  generalize the 600-cell evolution by allowing more variables and again
  study the causality-breaking endpoint and matter implementation.
- [Liu--Williams, arXiv:1510.05771](https://arxiv.org/abs/1510.05771)
  review Collins--Williams Regge cosmology for closed FLRW and lattice
  universes.
- [Liu--Williams, arXiv:1502.03000](https://arxiv.org/abs/1502.03000)
  study closed Regge lattice universes, including stability regions and
  perturbed-mass evolution.
- [Bahr--Dittrich, arXiv:0905.1670](https://arxiv.org/abs/0905.1670)
  explain why lapse-dependent pseudo-constraints arise when discrete
  diffeomorphism symmetry is broken.

These sources establish that 600-cell dust evolution, stopping points and
discrete lapse determination are **KNOWN mechanisms**.  A positive finite
root here would not by itself be a new physical tick.

The located papers do not use the exact certified cellular frustum action,
the same global conserved-dust term and the same fixed canonical-state
two-equation census printed above.  That narrower comparison remains
**OPEN** until their detailed conventions are audited; search alone does not
prove novelty.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The absolute time scale is forbidden by global classical scale covariance.
- The generic same-state branch exists formally through quadratic endpoint
  order and fails at cubic order.
- Finite Regge evolution may fix lapse-like data through pseudo-constraints.
- Published 600-cell dust evolutions can stop at a causality-breaking point.

### CONTROL

- Recover the exact `h=0,q=v` boundary solution for the frozen state.
- Reproduce the positive cubic cross-resultant from the local expansion of
  the finite elimination equation.
- Enumerate all simultaneous zero-slope cases before dividing by `Ch` or
  `Ph`.
- Require `h>0` and `1+hq>0` after reconstruction.
- Keep positive- and negative-velocity components separate and evaluate the
  `K=0` and `v=0` strata explicitly rather than excluding them.

### OPEN

- The complete real zero set of `D(v,q)` away from the trivial `h=0`
  boundary.
- Whether any positive finite solution is isolated for a fixed state.
- Whether a finite root is connected to the published 600-cell causality
  stop or is specific to the present action.
- Composition, stability, carrier refinement and action improvement.
- Any conversion of a dimensionless height to seconds or Planck time.
- External novelty of the exact elimination formula.

## 5. Acceptance boundary

An isolated positive `h/L_minus` is only a **candidate relational interval**.
It becomes evidence for a tick only after the same branch composes, is stable
under perturbations and survives carrier/action refinement without a fitted
cutoff.

No positive solutions on the complete domain would be a **DERIVED NEGATIVE**
for the finite homogeneous tick route of this action.  A continuum of
solutions would be **STRUCTURAL reparametrization**, not tick selection.  An
incomplete transcendental root classification remains **OPEN** regardless of
how dense a numerical scan is.
