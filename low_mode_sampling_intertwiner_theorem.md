# Exact low-mode sampling/intertwining theorem

Date: 2026-07-22

## Theorem (DERIVED)

Let `B` be the 12-point Hopf base obtained from the decagonal fibers of the
600-cell. Let `L_B` be the unweighted graph Laplacian of its adjacency graph.
The graph is the icosahedron. Choose the regular icosahedral embedding
`B -> S^2` supplied by the first three-dimensional eigenspace of `L_B`, and
let

`E : H_0 direct_sum H_1 direct_sum H_2 -> R^B`

be evaluation of real spherical harmonics at those 12 points. Then:

1. `E` is injective and has rank `1+3+5=9`.
2. Its image is exactly the sum of the `L_B` eigenspaces with eigenvalues
   `0`, `5-sqrt(5)`, and `6`.
3. Define

   `C_B = ((3-sqrt(5))/(4 sqrt(5))) L_B^2`
   `      + ((5 sqrt(5)-9)/(2 sqrt(5))) L_B`.

   On `H_0 + H_1 + H_2`, evaluation exactly intertwines the round-sphere
   scalar Casimir and this discrete operator:

   `C_B E = E (-Delta_S2)`.

4. With `R=(E^T E)^{-1}E^T`, one has `R E=I`; moreover `E R` is exactly the
   spectral projector onto `1+3+5` and `R` annihilates the remaining
   three-dimensional alias eigenspace.

The verifier is `reproducible/verify_low_mode_sampling_intertwiner.py`.

## Scope

- **DERIVED:** exact finite representation and reconstruction of scalar round
  `S^2` modes through `l=2`, including an exact Laplacian/Casimir intertwiner.
- **OPEN:** any refinement limit or convergence theorem beyond this fixed
  12-point sampling set.
- **OPEN:** vector-valued gauge potentials, local parallel transport, a Lie
  bracket, and a nonabelian continuum connection.

This theorem uses a route distinct from the known failed operator-compression
and signed-lift constructions: it is an exact scalar sampling theorem on a
band-limited subspace, not a claim about local transport.
