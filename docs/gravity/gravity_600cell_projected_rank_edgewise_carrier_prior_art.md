# Prior-art gate: canonical projected rank-edgewise 600-cell carrier

Date: 2026-08-19

Status: written before constructing the new full-carrier `Esd_2(sd K)` mesh
or evaluating any of its geometric diagnostics.  External novelty remains
**OPEN**.

## 1. Exact proposed object and hypotheses

Let `K` be the boundary complex of the regular 600-cell in its already
declared unit-quaternion realization in `S^3 subset R^4`.  Define

```text
K_n = P(Esd_(2^n)(sd K)),       P(x)=x/||x||,       n=0,1,...
```

as follows.

1. Apply barycentric subdivision `sd` exactly once.  A maximal chamber is a
   complete flag of nonempty faces of dimensions `(0,1,2,3)`.
2. Order its four barycentric vertices by those face dimensions.  This order
   is preserved by every automorphism of `K`; no global vertex identifier is
   used.
3. Apply the Edelsbrunner--Grayson edgewise subdivision at resolution `2^n`
   with that derived order.
4. Realize every barycentric/edgewise vertex by its affine barycentric
   combination in `R^4`, then radially normalize it to the unit three-sphere.
   Shared abstract vertices are merged before projection.
5. Use straight Euclidean chords between the projected vertices when a
   piecewise-flat spatial simplex is needed later.

The proposed first calculation constructs `K_0=sd K` and
`K_1=Esd_2(sd K)`, verifies their global topology, shape and `H4`
equivariance, and records their volume and chordal approximation diagnostics.
It does not yet evaluate a Lorentzian action or a fluctuation spectrum.

## 2. Why the projected-red carrier is not sufficient

The already tested projected `1 -> 8` red refinement divides the central
octahedron of every regular parent tetrahedron by one of three equal
opposite-midpoint diagonals.  This gives `3^600` possible first-level global
assignments.

The rotational stabilizer of a tetrahedron is `A4`.  Its action on the three
opposite-edge pairings factors through `A4/V4 = C3` and is transitive.
Therefore no diagonal is fixed by the cell stabilizer.  Since the 600
tetrahedra form one `H4` orbit, the standard fixed-point criterion for an
equivariant section implies:

> no assignment of one central-octahedron diagonal per tetrahedron is
> `H4`-equivariant.

This is the same obstruction independently encountered by the repository's
direct edgewise audit on an unranked tetrahedron.  The four red-refinement
tie rules used in the homogeneous acceleration calculation are legitimate
regulator probes, and their close agreement is evidence that the homogeneous
scalar coefficient is insensitive to those four choices.  They do not select
a unique anisotropic Hessian.  Computing such a spectrum before replacing
the regulator would make its eigenvectors choice-dependent.

## 3. KNOWN primary results

Edelsbrunner and Grayson construct edgewise subdivision of a simplex and
prove equal child volumes, face compatibility, composition in resolution and
a dimension-dependent finite number of shape classes:

- H. Edelsbrunner and D. R. Grayson, *Edgewise Subdivision of a Simplex*,
  DOI `10.1007/s004540010063`.

The repository has independently reconstructed the color schemes and proved
that applying `Esd_(2^n)` after one rank-ordered barycentric subdivision is
functorial, nested and uniformly shape regular in the affine
piecewise-Euclidean realization.  It has also reconstructed the negative
`A4` orbit of the three direct midpoint subdivisions.

Projected subdivisions of the 600-cell are established geodesic 4-dome
prior art.  Tsuda and Fujiwara subdivide 600-cell tetrahedra and octahedra,
project new vertices to the circumsphere, and use a pseudo-regular angular
average to obtain a closed-FLRW continuum limit:

- R. Tsuda and T. Fujiwara, *Oscillating 4-Polytopal Universe in Regge
  Calculus*, arXiv:`2011.04120`, DOI `10.1093/ptep/ptab079`, especially
  Section 6.

That construction establishes the broad projected-geodesic-dome programme.
It does not use the present rank-selected simplicial tower or the complete
direct irregular Regge action.

The physical interpretation of later linearized variables must respect the
scope of canonical Regge results.  On a flat background, Hoehn identifies
vertex-displacement gauge variables and lattice gravitons under Pachner
moves:

- P. A. Hoehn, *Canonical linearized Regge Calculus: counting lattice
  gravitons with Pachner moves*, arXiv:`1411.5672`, DOI
  `10.1007/JHEP05(2015)036`.

The present closed curved dust background is not flat, so those mode labels
cannot be copied without a new action-derived constraint audit.

## 4. CONTROL

- The source 600-cell must reproduce f-vector `(120,720,1200,600)` and every
  source vertex must lie on the unit sphere within its declared decimal
  precision.
- `sd K` must reproduce the already certified projected carrier
  `(2640,17040,28800,14400)`.
- The local rank-ordered `Esd_2(sd Delta^3)` complex must contain `192`
  tetrahedra and be invariant under all `24` permutations of the parent
  tetrahedron.
- The explicit 600-cell rotational cell stabilizer must induce all `12`
  elements of `A4` and act without a fixed central-octahedron diagonal.
- Radial projection must commute with the left/right binary-icosahedral
  rotations and quaternion conjugation within source-coordinate error.
- Shared-face incidence, Euler characteristic, positive volume and vertex
  merging must be checked globally rather than inferred from local counts.

## 5. OPEN before calculation

- the full f-vector and element-quality range of projected
  `Esd_2(sd K)`;
- whether maximum chord length and element-centroid radial sag decrease at
  the first edgewise step;
- whether radial projection preserves an adequate finite-level shape bound;
- whether total chordal volume moves toward `vol(S^3)=2 pi^2`;
- a uniform all-level shape theorem after nonlinear radial projection;
- the direct homogeneous Regge acceleration on this carrier;
- local lapse equations, a refined local dust quadrature, the constraint
  quotient and an anisotropic generalized spectrum;
- external novelty of this exact carrier/action combination.

## 6. Attack on the framing

The construction is canonical only after the homogeneous round-sphere
background is declared.  The map `P(x)=x/||x||` uses the distinguished centre
and ambient Euclidean norm.  Both are `H4`-invariant and already present in
the unit-quaternion realization, but choosing a round `S^3` target is a
**STRUCTURAL** background hypothesis, not an equation of motion.

Pure affine subdivision without `P` would be canonical and uniformly shape
regular, but it would only subdivide the same piecewise-flat 600-cell
geometry.  It can refine finite-element fields on that fixed polyhedron; it
cannot by itself improve the polyhedron toward a round sphere.  Radial
projection is therefore load-bearing rather than cosmetic.

Conversely, passing this carrier gate will not produce gravity.  It only
removes a regulator ambiguity before a later action calculation.  A positive
homogeneous Friedmann calibration would still be known Regge convergence,
not new physics.  The first physically discriminating result remains a
constraint-reduced inhomogeneous quadratic action stable under refinement.

## 7. Route decision

Proceed with the projected rank-edgewise carrier before adding Brown--Kuchar
dust fields.  The latter contain new proper-time, label and momentum variables
not selected by the current point-particle action.  Refining the already
declared geometry therefore introduces less external physical structure.
This decision does not claim that dynamical dust is unnecessary; it defers
that larger model extension until the geometric regulator is controlled.
