# Protocol: edgewise Whitney continuum-dynamics gate

Date: 2026-08-12

This protocol is committed before constructing or evaluating the exact
edgewise Whitney inclusions below.  It uses no particle, mass, coupling,
clock, `a1=5`, Planck or Standard-Model target.

## Question

The consistent Whitney metric gives exact Galerkin induction but no strict
finite-level propagation cone.  The former barycentric tower also failed the
shape-regularity hypothesis needed by standard finite-element exterior
calculus (FEEC) convergence results.

The canonical tower derived in protocol/result commits `58fa9fc` and
`0eddf27`,

\[
 K_n=\operatorname{Esd}_{2^n}(\operatorname{sd}K),
\]

is nested, conforming and uniformly shape regular.  The present gate asks:

1. does the exact all-degree Whitney--Kaehler--Dirac induction survive on
   this tower, rather than only on barycentric subdivision?
2. do the now-proved mesh properties place the family inside the hypotheses
   of established FEEC Hodge--Laplacian spectral convergence?
3. if so, what causal statement follows for the *continuum* Hodge--Dirac
   evolution, and what still does not follow physically?

## Fixed geometry and exact finite control

Use the rank-ordered barycentric orthoscheme of a regular tetrahedron.  Build
`Esd_1` and `Esd_2` independently from the Edelsbrunner--Grayson color-scheme
definition, using rational barycentric coordinates.  Do not import assembled
mass or inclusion matrices from an earlier verifier.

For every degree `p=0,1,2,3`:

1. assemble the exact consistent Whitney mass `M_(k,p)` from its defining
   affine integral;
2. construct the inclusion `P_p` by integrating each coarse Whitney form on
   each fine `p`-simplex, equivalently by the exact barycentric determinant;
3. verify exactly

\[
 d_fP_p=P_{p+1}d_c,
 \qquad
 P_p^*M_{f,p}P_p=M_{c,p};
\]

4. assemble the weak Kaehler--Dirac form `A=M D` and verify

\[
 P^*A_fP=A_c,
 \qquad
 \gamma_fP=P\gamma_c.
\]

Record, without an acceptance target, the exact strong-adjoint leakage ranks
and the exact row-sum-lumped isometry residual ranks in every degree.

Local affine naturality, face conformity and the edgewise composition theorem
are the all-level/global step.  The finite calculation is a defining identity
control, not an extrapolation of spectral numbers.

## Frozen theorem-hypothesis audit

The verifier must load the prior rank-edgewise certificate and require:

- exact face conformity;
- exact nesting `Esd_4` over `Esd_2`;
- exact `S4` equivariance;
- a level-independent set of three normalized shape classes;
- no direct-edgewise fixed variant without the one-time rank layer.

Together with the direct calculation, record the following implication chain:

\[
 \text{finite shape set}
 \Longrightarrow \text{uniform fullness/shape regularity},
\]

\[
 h_n=h_0,2^{-n}\longrightarrow0,
\]

\[
 \text{Whitney subcomplex + bounded commuting projections}
 \Longrightarrow \text{stable Hodge--Laplacian/eigenpair convergence}.
\]

The last implication is a cited mathematical theorem, not a finite verifier
claim.  The relevant primary references are:

- Edelsbrunner and Grayson, *Edgewise Subdivision of a Simplex*, DOI
  `10.1007/s004540010063` (composition and finite shape classes);
- Arnold, Falk and Winther, *Finite element exterior calculus: from Hodge
  theory to numerical stability*, arXiv:`0906.4325` (subcomplex, bounded
  cochain projection, Hodge--Laplacian and eigenvalue approximation);
- Christiansen, *Finite element systems of differential forms*,
  arXiv:`1006.4779` (cellular complexes, polyhedral grids and eigenpair
  approximation);
- Christiansen, *Stability of Hodge decompositions...*, arXiv:`1007.1120`
  (discrete Poincare/Rellich compactness on compact manifolds).

### Scope trap to test explicitly

The continuum carrier used by the exact masses is the fixed compact
piecewise-flat 600-cell boundary, not the round smooth 3-sphere.  If the
cellular/polyhedral FEEC theorems do not cover the glued piecewise-flat
Hilbert complex with its codimension-two Regge singularities, the continuum
claim must remain **OPEN ANALYTIC GAP**.  A smooth-manifold theorem may not be
silently substituted, and radial projection may not be introduced.

## Frozen continuum causal control

On the eight-dimensional exterior algebra of a Euclidean 3-dimensional
cotangent space, construct exterior multiplication `epsilon(xi)` and
contraction `iota(xi)` exactly.  For the Hodge--Dirac principal symbol

\[
 \sigma_D(\xi)=i\bigl(\epsilon(\xi)-\iota(\xi)\bigr),
\]

verify

\[
 \sigma_D(\xi)^2=|\xi|^2I.
\]

Therefore the continuum equation

\[
 i\partial_t\psi=(cD+\mu\gamma)\psi
\]

has characteristic speed `|c|`; the zeroth-order mass term does not change
the principal symbol.  Finite propagation is an established theorem for a
self-adjoint Dirac-type operator on the relevant complete continuum carrier.
This is conditional on the continuum gate above and on supplied `c,t`.

Retain the calibrated consistent-Whitney circle dispersion as a hostile
ultraviolet control.  Verify symbolically that its group velocity tends to
`c` at fixed physical momentum as `h -> 0`, while the already known
cutoff-scale value `sqrt(2)c` remains.  Passing the continuum gate must not be
reported as a strict causal cone at any finite mesh.

## Decision boundaries

### Finite acceptance

All exact inclusion, metric, weak-Dirac and grading identities pass.  Label
this **DERIVED EDGEWISE GALERKIN-INDUCTIVE DYNAMICS**.

### Continuum acceptance

The finite acceptance passes and every hypothesis of the cited
cellular/polyhedral FEEC result applies to the fixed piecewise-flat carrier.
Label spectral convergence **STRUCTURAL (ESTABLISHED THEOREM APPLIED)**, not
a new theorem proved by the verifier.

Then the principal-symbol identity and the standard finite-propagation
theorem give **DERIVED CONDITIONAL CONTINUUM CAUSALITY**: a finite speed for
the emergent continuum equation with externally supplied `t` and `c`.

### Analytic-gap verdict

If the Regge singular carrier is outside the cited theorem's stated scope,
retain finite acceptance but label continuum convergence and causality
**OPEN ANALYTIC GAP**.

### Refutation

Any exact inclusion, isometry, compression, conformity, nesting or uniform
shape gate fails.

No passing outcome derives a fourth dimension, Lorentzian time, the numerical
value of `c`, mass, `hbar`, Newton's `G` or Planck units.
