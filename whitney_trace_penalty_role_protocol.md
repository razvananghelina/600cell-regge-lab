# Protocol: what the Whitney trace penalty measures

Date: 2026-08-11

This protocol is frozen before the new kernel, valence, and scaling audit.  It
attacks the interpretation of the recently computed trace-stiffness spectra.

## Complete hypotheses

For each mesh `K_h` in the canonical rank-edgewise tower, fix a form degree
`p=0,1,2` and let

- `V_h^disc` be the direct sum of local Whitney `p`-cochains over all
  tetrahedra;
- `V_h^conf` be the globally assembled Whitney cochains;
- `J_h : V_h^conf -> V_h^disc` be the copy injection;
- `R_h` be the complete shared-face copy-difference map;
- `M_h>0` be the block-local exact Whitney mass;
- `H_h>0` be the exact Whitney trace mass on the shared faces;
- `B_h=R_h^* H_h R_h` be the trace-jump penalty.

No claim below applies to a different discontinuous Galerkin form containing
volume derivatives, consistency fluxes, or upwind terms.  It applies to this
pure jump penalty and to the already tested pencil

\[
 W_h+\kappa B_h.
\]

## Claim A: zero physical compression

The exact assembly identity should be

\[
 R_hJ_h=0,
 \qquad
 J_h^*B_hJ_h=0.
\]

If `H_h` is positive and every occurrence graph is connected, then

\[
 \ker B_h=\ker R_h=\operatorname{im}J_h.
\]

Therefore the penalty has zero quadratic form on every conforming cochain,
not only on constants or harmonic modes.  The audit must exhibit a nonconstant
conforming scalar cochain with zero jump but nonzero simplicial coboundary.

**Interpretation fixed in advance:** if these identities pass, the pure trace
penalty has zero tangent/principal action on the physical conforming sector.
Its positive spectrum measures copy mismatch, not propagation of an assembled
Whitney field.

## Claim B: uniform separator scaling

The same negative does not make the penalty useless.  Under all of the
following hypotheses:

1. a conforming, quasi-uniform tetrahedral tower with `h -> 0`;
2. a level-independent finite set of normalized tetrahedron and face shapes;
3. a level-independent bound on the number of tetrahedra incident to each
   simplex;
4. exact Whitney element and trace masses;

affine scaling and finite-dimensional norm equivalence imply

\[
 c_p h^{3-2p}I\le M_{T,p}\le C_p h^{3-2p}I,
\]

\[
 c'_p h^{2-2p}I\le H_{F,p}\le C'_p h^{2-2p}I.
\]

The nonzero singular values of the bounded occurrence-graph incidences are
then bounded away from zero.  Consequently the smallest positive generalized
penalty eigenvalue and the local Kähler--Dirac norm obey

\[
 c h^{-1}\le g_h\le C h^{-1},
 \qquad
 a_h\le C_D h^{-1},
\]

and hence `a_h/g_h` is uniformly bounded.

This is a theorem conditional on all four stated hypotheses.  Finite-level
numerics may control the construction but may not be cited as proof of the
all-level bound.

## Exact and numerical gates

For the complete closed control at edgewise resolutions `k=1,2,4`:

1. construct `J_h` and `R_h` independently as sparse integer matrices;
2. check `R_h J_h=0` exactly;
3. check every occurrence graph is connected and
   `rank R_h = dim V_h^disc - dim V_h^conf` combinatorially;
4. record the maximum occurrence count in every degree;
5. exhibit the nonconstant conforming zero-jump/nonzero-coboundary witness;
6. reuse only the accepted `9/9` rank-edgewise stiffness certificate and
   record `h a`, `h g`, `a/g`, and the separation threshold `2a/g` with
   `h=1/k`;
7. confirm the exact affine mass powers from the independently certified
   uniform-dilation control.

No value of `kappa` is selected.  The levelwise `2a/g` values are necessary
threshold diagnostics for the earlier sufficient separation bound, not
measurements of a physical coupling.

## Acceptance and kill boundaries

- **DERIVED NEGATIVE FOR DYNAMICS** if Claim A passes: the pure penalty cannot
  be the missing propagation operator on assembled fields.
- **DERIVED STRUCTURAL POSITIVE FOR SEPARATION** if the four hypotheses of
  Claim B are established for the rank-edgewise tower: a finite dimensionless
  `kappa` can in principle keep the conforming and mismatch sectors uniformly
  separated under refinement.
- **KILL FOR UNIFORM SEPARATION** if normalized shapes or occurrence sizes are
  unbounded, or if the exact Whitney/trace scale powers differ by degree.

Even a positive separation theorem does not select `kappa`, prove spectral
convergence of the coupled pencil, provide a tensor-product causal system, or
derive time, mass, `c`, `hbar`, `G`, or a Planck scale.

