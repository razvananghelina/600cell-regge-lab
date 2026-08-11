# Preregistration: complete 600-cell exact trace stiffness

Date: 2026-08-11

## Question

The exact Whitney trace-jump energy improved degree balance on the complete
boundary-of-4-simplex control.  This test asks whether that improvement
survives on the actual complete 600-cell and its full first barycentric
subdivision.

The operator is frozen by commits `b9a4104` and `a92c911`.  No new weight,
coefficient, or candidate is introduced.

## Complete carriers

Use exactly:

- base f-vector `(120,720,1200,600)`, duplicated dimension 9,000;
- refined f-vector `(2640,17040,28800,14400)`, duplicated dimension 216,000;
- exact regular-parent and barycentric-child Whitney element masses;
- every shared-triangle trace jump in degrees (p=0,1,2);
- exact face Whitney Gram blocks (M_{F,p}).

The positive weak penalty is

\[
 B_p=R_p^*H_pR_p,
 \qquad
 H_p=\operatorname{diag}_F(M_{F,p}).
\]

The candidate count remains (N=1).

## Quotient spectrum without zero-mode contamination

The complete refined matrices have exact conforming nullities too large for a
direct smallest-eigenvalue query.  Let (V_p) be an orthonormal basis of
(operatorname{im}R_p), constructed as the direct sum of incidence-image
bases of every connected occurrence graph.

The positive generalized spectrum of

\[
 B_pv=\lambda M_pv
\]

equals the spectrum of the strictly positive reduced pencil

\[
 \left(
 V_p^*H_pR_pM_p^{-1}R_p^*H_pV_p
 \right)x
 =
 \lambda
 \left(V_p^*H_pV_p\right)x.
\]

Both reduced matrices are applied as sparse linear operators.  No assembled
mass inverse, dense full matrix, or fitted deflation shift is formed.

## Frozen iterative protocol

Use SciPy block LOBPCG with:

- deterministic seed `60020260811 + 100*level + degree`;
- block size 5;
- `largest=False` for the five smallest quotient eigenpairs;
- `largest=True` for the five largest quotient eigenpairs;
- tolerance (10^{-9});
- maximum 2,000 iterations;
- no post-result preconditioner or tolerance change.

The reported gap and maximum are the outer eigenvalues of those blocks.  A
pair is certified only if its directly recomputed generalized relative Ritz
residual is below (10^{-7}).

This is a numerical extremal certificate, not an exact algebraic spectrum.

## Mandatory dense calibration

Before the complete calculation, apply the identical quotient construction
and iterative solver to the base and refined boundary-of-4-simplex control.
It must reproduce all six dense exact-trace gaps and maxima stored by
`verify_whitney_trace_stiffness.py`:

\[
g_{0,p}=(8.66025403784,5.19615242271,3.46410161514),
\]

\[
g_{1,p}=(2.29366260722,2.35146405514,2.80283446025),
\]

with relative value error below (5\times10^{-7}) and relative Ritz
residual below (10^{-7}).  The six maxima are read from that committed JSON
certificate rather than copied by hand.

Failure of this calibration invalidates the solver but is not evidence
against the physical operator.

## Frozen complete outputs

For each complete level and degree record:

1. row count, exact rank, redundancy, and connectedness;
2. face metric types and exact parent-face agreement;
3. smallest five and largest five quotient eigenvalues;
4. their maximum recomputed Ritz residual;
5. positive gap (g_{r,p}) and maximum (q_{r,p});
6. (s_{r,p}=a_r/g_{r,p});
7. first-step ratios (R_p=s_{1,p}/s_{0,p});
8. ratio spread (max R_p/\min R_p).

The already committed unweighted complete certificate
`whitney_stiffness_refinement.json`, protocol `03e0abc`, is the paired
comparison because it uses these exact same two complete carriers and local
element metrics.

## Decision protocol

- If all (R_p=1) within (5\times10^{-7}), label exact trace stiffness
  **DERIVED NUMERICAL: first-step compatible on the complete 600-cell**.
- Otherwise label exact compatibility **DERIVED NUMERICAL NEGATIVE** and
  report every ratio.
- If the complete trace ratio spread is smaller than the paired unweighted
  spread, label **PATTERN: improved complete degree balance**.
- If it is equal or larger, the proposed geometric repair reaches its kill
  boundary.
- Even an improvement does not derive a flow or exponent from two levels.
- The overall (kappa) normalization remains OPEN in all outcomes.

## Scope attacks fixed in advance

The exact trace form is unique only within the stated (L^2)-jump contract.
Other discontinuous-Galerkin fluxes, auxiliary face fields, or chiral
dilations are different theories and are not covered.

The positive penalty remains degree-preserving and breaks Kähler--Dirac
oddness at finite (kappa).  Improved stiffness scaling would not solve that
gate or the singular exact-(kappa\to\infty) limit.

## Status before execution

- **DERIVED INPUT:** exact trace operator and uniform (h^{-1}) scaling.
- **DERIVED INPUT:** dense control spectra and paired unweighted complete
  spectra.
- **OPEN:** calibration of the quotient LOBPCG method.
- **OPEN:** complete trace gaps and degree ratios.
- **OPEN:** repeated refinement, overall stiffness, chirality, and causal
  dynamics.
- **NOT CLAIMED:** physical time, mass, inertia, (c), (hbar), Newton's
  (G), or a Planck scale.
