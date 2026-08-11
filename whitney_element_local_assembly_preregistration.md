# Preregistration: can element-local Whitney dynamics preserve assembly?

Date: 2026-08-11

## Question frozen before computation

The signed Grover--Szegedy tick is local, unitary and vertex-isotropic, but
its unweighted discriminant is not the consistent Whitney Kähler--Dirac
operator under tetrahedral refinement.  The global strong Whitney adjoint
contains an inverse mass matrix and loses one-step locality.

Test the minimal carrier enlargement that could separate these roles:

> Duplicate every simplex coefficient once for every incident tetrahedron,
> apply the exact metric operator independently inside each tetrahedron, and
> ask whether the conforming assembled subspace is invariant.

No dispersion, particle, speed, mass or phenomenological target enters this
test.

## Frozen carriers

For each of the 600 tetrahedra use all 15 nonempty local faces:

\[
4+6+4+1=15.
\]

The discontinuous element carrier therefore has dimension

\[
600\times15=9000.
\]

For degree `p`, let `J_p` restrict one global oriented cochain coefficient to
all of its tetrahedron-local copies, including the orientation sign.  The
conforming subspace is `im(J_p)`.

Before dynamics, verify the complete occurrence counts:

- each vertex occurs in 20 tetrahedra;
- each edge occurs in 5 tetrahedra;
- each triangle occurs in 2 tetrahedra;
- each tetrahedron occurs once.

Thus the duplicated dimensions must be `(2400,3600,2400,600)`, while the
global dimensions remain `(120,720,1200,600)`.

## Frozen local metric operator

Use the exact rational Whitney mass matrices of a regular affine
tetrahedron.  Every 600-cell facet is congruent, and uniform rescaling does
not affect the zero/nonzero invariance question.

Let

\[
M_p^{\rm loc}=\bigoplus_{t=1}^{600}M_{p,t},
\]

and let `d_p^loc` be the direct sum of the usual signed tetrahedron
coboundary.  Define the element-local codifferential

\[
\delta_p^{\rm loc}
=(M_p^{\rm loc})^{-1}(d_p^{\rm loc})^TM_{p+1}^{\rm loc}.
\]

All local matrices must be generated from the defining Whitney-form
integral, not copied as unexplained constants.

## Exact assembly tests

1. Verify the upward operator preserves conformity exactly:

   \[
   d_p^{\rm loc}J_p=J_{p+1}d_p.
   \]

2. For every downward degree `p=0,1,2`, test whether there exists any global
   matrix `X_p` such that

   \[
   \delta_p^{\rm loc}J_{p+1}=J_pX_p.
   \]

   This is an exact subspace-membership question.  Compute the unique
   equality-average candidate

   \[
   X_p=(J_p^TJ_p)^{-1}J_p^T
       \delta_p^{\rm loc}J_{p+1}
   \]

   and record the residual, its nonzero count, rank, and an exact rational
   witness whenever it fails.  Do not infer a general theorem from one
   witness alone; evaluate all three complete matrices.

3. Assemble

   \[
   M_p=J_p^TM_p^{\rm loc}J_p.
   \]

   Verify the weak identity

   \[
   J_p^TM_p^{\rm loc}\delta_p^{\rm loc}J_{p+1}
   =d_p^TM_{p+1}.
   \]

   This identity is the coefficient-level statement that the
   `M_loc`-orthogonal projection back to `im(J_p)` reproduces the global
   Whitney adjoint

   \[
   \delta_p=M_p^{-1}d_p^TM_{p+1}.
   \]

4. Record the support of every local codifferential block and the Hasse
   distance of any new within-element coupling.

## Numerical/exact standard

The local mass and codifferential blocks are exact SymPy rational matrices.
The global restriction and incidence matrices have integer entries.  Clear
all local denominators before the support and assembly comparisons so the
zero tests are exact integers.

Numerical singular values may be used only to report residual ranks after an
exact nonzero certificate exists.  A tolerance-only failure is not accepted.

## Decision boundaries

- **DERIVED LOCAL METRIC INVARIANCE:** all three downward maps preserve every
  conforming subspace.  This would provide the missing element-local metric
  generator without an assembly projection.
- **DERIVED ASSEMBLY-LEAKAGE NO-GO:** at least one complete downward map has
  exact nonzero leakage.  Report the hit fraction over all three degrees.
- **DERIVED PROJECTED FACTORIZATION:** independently, all three weak assembly
  identities hold exactly.  Then the correct metric operator is local
  element evolution followed by metric projection, and the unresolved
  locality problem is isolated in that projection.

## Hostile scope boundary

Failure closes only the direct-sum element operator as an invariant dilation.
It does not close multi-step consensus/gluing walks, discontinuous Galerkin
fluxes, larger ancillas or approximate local solvers for the metric
projection.

Conversely, the algebraic projected factorization is not yet a unitary tick:
an instantaneous orthogonal projection is not reversible, and its exact
reflection may be nonlocal because it contains the assembled inverse mass.
Any continuation must implement the glue reversibly rather than silently
project after each step.
