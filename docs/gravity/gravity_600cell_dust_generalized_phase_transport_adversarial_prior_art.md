# Prior-art and independence gate: adversarial phase-transport replication

Date: 2026-08-18

## Exact purpose

The preregistered high-precision calculation has already reported
`GENERALIZED_PHASE_TRANSPORT_REFUTED`.  This is therefore not a blind discovery
test and cannot repair provenance after the result is known.  Its narrower
purpose is to try to falsify an implementation error by recomputing the
decisive subspace statement through a mechanically different path.

The frozen claim under attack is

```text
T_2 (E_old direct-sum E_old*) is not contained in
    (E_shifted direct-sum E_shifted*).
```

The audit must not use the high-precision residual verifier's projectors,
four-block leakage formula or tangent reconstruction as its decisive numeric
path.

## Primary numerical literature

**KNOWN.** Knyazev and Argentati formulate principal angles between subspaces
and SVD-based algorithms for their cosines and sines, including the warning
that cosine-only calculations lose accuracy for small angles:
<https://doi.org/10.1137/S1064827500377332>.

**KNOWN.** Knyazev, Jujunashvili and Argentati relate squared cosines of
subspace angles to Ritz values and characterize the corresponding projectors:
<https://doi.org/10.1016/j.jfa.2010.05.018> and
<https://arxiv.org/abs/0705.1023>.

**KNOWN.** The LAPACK Users' Guide reduces the Hermitian-definite problem
`A z = lambda B z`, with `B=L L^H`, to the ordinary Hermitian problem
`(L^-1 A L^-H)y=lambda y`:
<https://www.netlib.org/lapack/lug/node54.html>.

These are standard validation tools.  They supply no evidence that the present
600-cell subspaces should close or fail to close.  External novelty remains
**OPEN**.

## Mechanically distinct path

The accepted calculation used:

1. residual-certified `mpmath` generalized projectors;
2. a newly reconstructed Flint tangent ball;
3. the projector residual `(I-Q1) T Q0`, both whole and by four blocks.

This audit will instead use:

1. the earlier direct-precision `M,V` midpoints, before the residual
   refinement;
2. an independently coded null-space reduction followed by explicit Cholesky
   whitening and `numpy.linalg.eigh`, rather than `scipy.linalg.eigh(A,B)`;
3. the earlier committed two-step tangent archive, which predates the phase
   calculation, rather than rebuilding the tangent through its helper;
4. QR bases for `T F_old`, direct sine-of-principal-angle leakage into the
   orthogonal complement of `F_shifted`, and a least-squares image residual.

The common dependency is the same derived Regge action and carrier.  The audit
is independent at the decisive spectral solver, tangent source and subspace
comparison, but is not an independent physical model.

## Required adversarial controls

- A synthetic map `W1 W0^H` must transport the source basis exactly into the
  target basis.
- A synthetic map replacing one target direction by a target-complement
  direction must be detected as nonclosing.
- Rephasing and reversing every source/target basis must leave all subspace
  measures invariant.
- The canonical cotangent lift `diag(U,conjugate(U))` is decisive.  The
  noncanonical same-representation and swapped-dual lifts are reported only as
  convention stress tests and cannot replace the frozen object.

## Epistemic limit

This float64 path cannot replace the residual-certified error proof.  A match
is **STRUCTURAL INDEPENDENT CORROBORATION** of that proof.  A disagreement is
more important: under project rule 4 the consolidated verdict returns to
**OPEN** until resolved.

