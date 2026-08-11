# Preregistration: can a finite pure constraint remove Whitney assembly leakage?

Date: 2026-08-11

## Fixed question

The exact element-local Whitney operator on the 9,000-dimensional duplicated
tetrahedron carrier leaks out of the 2,640-dimensional conforming cochain
subspace.  Metric projection repairs it, but projection is not a reversible
microscopic step.

Test the most general additive **pure-constraint penalty**, not a selected
matrix ansatz.

Let

\[
\mathcal H_{\rm loc}=\bigoplus_{p=0}^3\mathcal H_p^{\rm loc},
\qquad
J=\bigoplus_{p=0}^3J_p,
\qquad
\mathcal C=\operatorname{im}J,
\]

and let (D_{\rm loc}=d_{\rm loc}+\delta_{\rm loc}) be exactly the operator
already fixed by protocol commit `162ce61`.  Let (P) be **any** linear
projection onto \(\mathcal C\).

A pure constraint is any linear operator (K) satisfying

\[
KJ=0.
\]

This includes every positive or self-adjoint penalty whose kernel contains
the conforming states, but the theorem below does not need positivity,
self-adjointness, locality or a chosen equality graph.  For an arbitrary
finite scalar κ define

\[
H_\kappa=D_{\rm loc}+\kappa K.
\]

The question is generator-level invariance:

\[
H_\kappa\mathcal C\subseteq\mathcal C.
\]

## Frozen universal identity

Before any new calculation, the hypotheses imply algebraically

\[
(I-P)H_\kappa J
=(I-P)D_{\rm loc}J+\kappa(I-P)KJ
=(I-P)D_{\rm loc}J.
\]

Thus a pure penalty cannot change the instantaneous quotient leakage at any
finite κ.  The only empirical input still requiring certification is the
rank of the fixed leakage map.

The prior verifier used exact equality-difference maps (Q_p) with

\[
\ker Q_p=\operatorname{im}J_p.
\]

It reported exact full-column ranks for

\[
Q_p\delta_p^{\rm loc}J_{p+1},\qquad p=0,1,2,
\]

namely `(720,1200,600)`.  Recompute these matrices rather than trusting the
stored report.  Also verify (Q_pJ_p=0) and that every (Q_p) has rank
`dim(H_p^loc)-dim(C_p)`, establishing its kernel exactly.

Because the upward maps preserve conformity and the three downward leakage
blocks land in distinct form degrees, the full quotient leakage must then
have rank

\[
720+1200+600=2520
\]

on the 2,520-dimensional positive-degree conforming sector.  The only
conforming inputs with no leakage are the 120 degree-zero inputs.

## Decision boundary

- **DERIVED FINITE PURE-PENALTY NO-GO:** all equality-difference kernels and
  the three full leakage ranks are certified exactly.  Then no finite scalar
  and no operator (K) satisfying (KJ=0) makes the conforming subspace
  invariant under (D_{\rm loc}+\kappa K).
- **REFUTED:** any claimed difference kernel is larger than the conforming
  subspace, or any recomputed leakage rank differs from the frozen values.

No coefficient scan is permitted or needed: the identity covers the whole
linear class at once.

## What this does *not* establish

State every boundary explicitly:

1. It does not prove that a physical infinity exists.  It says only that an
   exact invariant generator cannot be obtained at finite penalty strength
   while the penalty vanishes on all physical states.
2. It does not prove that the singular κ→∞ limit produces the desired
   global Whitney operator.  That additionally requires a selected
   self-adjoint constraint whose kernel projection is the metric projector
   and a controlled finite-dimensional low-energy limit.
3. It does not exclude a finite **non-pure** correction (C) with (CJ\ne0).
   Such a correction must act on conforming states and cancel the quotient
   leakage exactly.
4. It does not exclude stroboscopic return at one selected tick,
   time-dependent controls, a product of local reflections, ancillas,
   discontinuous-Galerkin fluxes, or an exact polynomial implementation of a
   projector/reflection.
5. It does not select a penalty scale, a physical time, a mass, inertia or
   the speed of light.

The targeted verifier alone will be run.  The full suite is excluded by the
user's current instruction.
