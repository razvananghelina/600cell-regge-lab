# Protocol: broken-FEEC reconciliation of the finite-stiffness mechanism

Date: 2026-08-11

This protocol is committed before constructing the new three-dimensional
CONGA matrices or evaluating their kernels.  It uses no phenomenological
target and does not select a projection by spectral performance.

## Prior-art correction

The circle construction is not a new numerical principle.  It is the
one-dimensional lowest-order instance of the projection-based broken finite
element exterior calculus (broken-FEEC/CONGA) framework of Martin Campos Pinto
and Yaman Güçlü, *Broken-FEEC discretizations and Hodge Laplace problems*,
Mathematics of Computation, DOI
[`10.1090/mcom/4085`](https://doi.org/10.1090/mcom/4085), preprint
[`arXiv:2109.02553v3`](https://arxiv.org/abs/2109.02553).

That work defines a broken differential by `d_h=d_pw P`, its Hilbert adjoint
using the broken metric, and a symmetric stabilization.  Its kernel theorem
states, under the complete Hilbert-complex hypotheses, that any positive
stabilization recovers precisely the harmonic fields of the conforming
subcomplex.

The present calculation does not claim to rediscover that theorem.  It tests
whether its hypotheses and algebraic identities hold on this project's
specific duplicated Whitney carrier, and whether the projection ambiguity
found in commit `add2d2b` changes the protected kernel.

## Frozen carrier and candidates

Use the closed controls

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,
\]

with the exact piecewise-Euclidean Whitney masses and oriented simplex bases.
Let `J_p` be the conforming copy injection and use both preregistered local
left inverses:

1. `L^C_p`: equal counting recovery;
2. `L^D_p`: diagonal-Whitney (mass-lumped) recovery.

For `X in {C,D}`, set

\[
 P^X_p=J_pL^X_p.
\]

Both candidates are included.  There is no post-result choice and no fitted
mixture.

## Frozen broken-FEEC operators

Let `D_pw,p` be the block-diagonal local simplicial coboundary and let
`M_loc,p` be the block-diagonal exact local Whitney mass.  Define

\[
 d^X_{h,p}=D_{\rm pw,p}P^X_p
\]

and its exact broken-metric adjoint in coefficient form

\[
 (d^X_{h,p})^*
 =M_{\rm loc,p}^{-1}(d^X_{h,p})^T M_{\rm loc,p+1}.
\]

For degree `p` and the frozen positive value `alpha=1`, use the symmetric weak
Hodge-Laplace pencil

\[
 \begin{split}
 K^X_p={}&(d^X_{h,p})^TM_{p+1}d^X_{h,p}\\
 &+((d^X_{h,p-1})^*)^TM_{p-1}(d^X_{h,p-1})^*\\
 &+(I-P^X_p)^TM_p(I-P^X_p),
 \end{split}
\]

against `M_p`.  Missing endpoint terms are zero.

`alpha=1` is not a physical parameter claim.  Positivity makes the exact
kernel independent of its magnitude; it is one numerical control only.

## Frozen exact and numerical gates

For both candidates, both levels and all degrees:

1. verify `P_p^2=P_p`, `P_pJ_p=J_p`, and `L_pJ_p=I`;
2. verify `D_pw,p J_p=J_{p+1}d_p` with the global simplicial coboundary;
3. verify `P_{p+1}d_{h,p}=d_{h,p}` and
   `d_{h,p+1}d_{h,p}=0`;
4. verify broken-metric adjointness with maximum residual below `1e-11`;
5. verify symmetry and positive semidefiniteness of every weak pencil;
6. on the complete `k=1` control, compute the lowest eight generalized
   eigenvalues without a target fit and require nullities by degree
   `(1,0,0,1)` at absolute threshold `1e-9`;
7. at `k=2`, certify the same kernel algebraically from positivity,
   `ker(I-P)=im(J)`, and the Betti numbers of the subdivided 3-sphere;
8. record whether the two candidates give identical or different broken
   differentials and positive spectra.  No agreement is required for a PASS.

The exact sparse identities are evaluated over floating realizations of
rational matrices, with nonzero structural patterns independently checked.
Residuals below `1e-11` are labelled **DERIVED NUMERICAL**.  The projection
identities also follow algebraically from `LJ=I`, and nilpotency follows from
the verified conforming intertwiner; those implications are labelled
**DERIVED** with every hypothesis stated.

## Frozen interpretation

- If both candidates have the conforming harmonic kernel, report a
  **DERIVED APPLICATION OF KNOWN THEOREM**: topology is robust to this
  projection ambiguity, although the finite operator remains unselected.
- If either candidate has another kernel, report a **DERIVED REFUTATION OF
  THE ASSUMPTION MAPPING**: at least one hypothesis of the broken-FEEC theorem
  was mapped incorrectly.
- If their positive spectra differ, report **STRUCTURAL AMBIGUITY**, not a
  winner.
- A physical Kähler--Dirac evolution is **not** established: the stabilized
  Hodge Laplacian is positive and topology-correct, but its stabilization is
  not the square of a newly derived odd first-order operator.

## Scope exclusions

- no claim of mathematical novelty for broken-FEEC/CONGA;
- no selection of `P^C` over `P^D`;
- no continuum fit or extrapolation from two levels;
- no derivation of Lorentzian time, causal speed, inertia, mass or Planck
  units;
- no full suite run, by explicit user request.

