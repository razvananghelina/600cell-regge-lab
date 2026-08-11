# Finite pure penalties cannot repair Whitney assembly leakage

Date: 2026-08-11

Preregistration commit: `4ce5d3e`

Targeted verifier:
`reproducible/verify_whitney_element_local_assembly.py`

Targeted result: **16/16 PASS**.  The full suite was not run, by explicit
request.

## Complete hypotheses

Let \(D_{\rm loc}=d_{\rm loc}+\delta_{\rm loc}\) be the exact element-local
Whitney Kähler--Dirac operator on the 9,000-dimensional duplicated carrier,
let (J) inject the 2,640 global cochains as equal oriented copies, and let
(P) be any projection with range \(\operatorname{im}J\).

The result applies to every finite scalar κ and every linear operator (K)
such that

\[
KJ=0.
\]

No locality, positivity or self-adjointness assumption is needed.  Therefore
the claim automatically covers the smaller physical class of local positive
self-adjoint constraint Hamiltonians whose kernel contains the conforming
states.

## Universal no-go

For

\[
H_\kappa=D_{\rm loc}+\kappa K
\]

one has identically

\[
(I-P)H_\kappa J
=(I-P)D_{\rm loc}J+\kappa(I-P)KJ
=(I-P)D_{\rm loc}J.
\]

Thus a pure penalty does not merely fail for the tested equality Laplacian:
it cannot change the generator-level leakage at all.  This is a theorem for
the whole class, not evidence from a coefficient search.

The targeted verifier independently reconstructs the complete 600-cell and
the exact rational Whitney blocks.  Its equality-difference maps have exact
ranks

\[
(2280,2880,1200,0)
\]

and kernels of dimensions

\[
(120,720,1200,600),
\]

exactly the four conforming cochain spaces.  Both the structural-rank upper
bounds and finite-field lower bounds coincide, and every product (Q_pJ_p)
is exactly zero.

The recomputed downward quotient-leakage ranks are

\[
(720,1200,600).
\]

They land in different form degrees, so the complete quotient map has rank

\[
720+1200+600=2520
\]

on a positive-degree conforming domain of dimension 2,520.  It is therefore
injective there.  The complete conforming kernel has dimension 120 and
consists precisely of degree-zero inputs, whose upward differential already
preserves conformity.

> **DERIVED FINITE PURE-PENALTY NO-GO:** no finite additive penalty that
> vanishes on every conforming state can make the conforming subspace
> invariant under the element-local Whitney generator.

## Does this force an infinity?

No.  It gives a conditional and much more precise answer.

Suppose (D) is self-adjoint in the chosen metric, (P) is the corresponding
orthogonal conforming projector, and (Q=I-P).  The finite correction

\[
C_*=-QDP-PDQ
\]

is self-adjoint and gives

\[
D+C_*=PDP+QDQ,
\]

which preserves both sectors exactly.  Among self-adjoint corrections with
no diagonal blocks, its off-diagonal action is forced; it is not fitted.
Moreover,

\[
C_*P=-QDP\ne0,
\]

so it evades the no-go precisely by acting nontrivially on physical states.

This proves that an infinite parameter is **not mathematically necessary**
for exact block invariance.  But (C_*) is not yet a local construction: it
is written using the desired metric projector (P), whose Whitney formula
contains the inverse assembled mass.  It repackages the missing operation
rather than deriving a microscopic implementation of it.

If one insists instead that a constraint vanish on physical states, exact
finite-κ invariance is impossible.  A large-gap or quantum-Zeno construction
can only approach an effective projected generator through the singular
limit κ→∞ (or an equivalent time-scale separation).  Whether such a limit
exists for a geometrically selected local constraint, and whether its kernel
projection is the required Whitney metric projector, remain separate
questions.

## Scope attacks

This result does **not** exclude:

- an exact finite product of local reflections;
- stroboscopic cancellation at one geometrically selected tick;
- time-dependent controls or ancillas;
- a finite non-pure coupling whose action on conforming states is selected by
  local geometry;
- discontinuous-Galerkin flux carriers;
- an exact polynomial implementation of a projector or reflection.

Nor does it establish a physical infinite quantity.  A limit is a feature of
one enforcement strategy, not an ontology.

## Status ledger

- **DERIVED:** every equality-difference kernel is exactly the conforming
  subspace.
- **DERIVED:** the quotient leakage rank is `2520/2520` on positive-degree
  conforming cochains.
- **DERIVED NEGATIVE:** every finite additive pure penalty fails at generator
  level, independently of its coefficients.
- **DERIVED:** a finite abstract non-pure block-cancelling correction exists.
- **STRUCTURAL:** that correction uses the global metric projector and is not
  yet a local microscopic rule.
- **STRUCTURAL:** interpreting κ→∞ as a low-energy or Zeno mechanism.
- **OPEN:** a coefficient-free finite local realization of the metric glue.
- **OPEN:** a selected constraint scale and controlled refinement limit.
- **NOT CLAIMED:** mass, inertia, physical time, Lorentz invariance or the
  measured speed of light.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_element_local_assembly.py
```

Expected result: `16/16`.
