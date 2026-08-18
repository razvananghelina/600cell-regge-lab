# Prior-art gate: geometric conjugacy of the two dynamic tangent maps

Date: 2026-08-17

## Exact object, carrier and hypotheses

The input is the already committed pair of calibrated `60 x 60` canonical
tangent maps about the first non-static fixed-mass dust tick.  Each map acts on
the order-24 invariant boundary phase space with 30 logarithmic edge-orbit
coordinates and their 30 canonical momenta.  The even and odd slabs use the
same 600-cell, five binary-tetrahedral cover cells and order-24 stabilizer, but
the odd schedule exchanges the first two phase cells.

The question is whether their calibration-indistinguishable characteristic
spectra are explained by a map supplied by the finite carrier:

1. an `H4` action that maps the complete even two-layer slab to the complete
   odd slab, either preserving or reversing its two boundary layers; or
2. failing that, a spatial `H4` action that maps the 30-orbit boundary
   partition to itself and whose canonical or anti-canonical phase lift
   intertwines the two tangent maps.

The complete-slab maps are geometrically canonical candidates.  Boundary-only
maps are retained separately as **STRUCTURAL** candidates because they need
not map the four-dimensional triangulation or action.

No arbitrary similarity matrix, eigenvector matrix, Schur-block coefficient,
continuum eigenvalue, speed or desired degeneracy is admissible.

## Primary literature

- Dittrich and Hoehn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>, derive action-generated pre/post canonical
  evolution under simplex gluing and Pachner moves.  This establishes how a
  genuine carrier isomorphism must act on boundary data, but does not imply
  that two curved triangulations give conjugate maps.
- Hoehn, *Canonical linearized Regge Calculus: counting lattice gravitons with
  Pachner moves*, <https://arxiv.org/abs/1411.5672>, derives gauge and curvature
  degrees of freedom around flat backgrounds.  The present background is
  curved and dust-filled, so its flat-background symmetry conclusions cannot
  be imported.
- Dittrich and Steinhaus, *Path integral measure and triangulation independence
  in discrete gravity*, <https://arxiv.org/abs/1110.6866>, show that
  triangulation independence is special rather than automatic, particularly
  when passing from 3D to 4D Regge calculus.
- Dittrich, Kaminski and Steinhaus, *Discretization independence implies
  non-locality in 4D discrete quantum gravity*,
  <https://arxiv.org/abs/1404.5288>, exhibit explicit 4D restrictions on local
  triangulation independence.  Hence shared spectra must be demonstrated, not
  inferred from the topology of the slabs.
- Watkins, *Product Eigenvalue Problems*,
  <https://doi.org/10.1137/S0036144504443110>, supplies a general setting in
  which cyclically reordered matrix products can share spectra.  It does not
  show that the present maps possess such a factorization.

No located primary source constructs these two order-24 dust tangent maps or a
conjugacy between their schedule parities.  External novelty remains **OPEN**;
a search cannot establish absence.

## Existing controls

- **DERIVED:** the consecutive old-to-final boundary orbit map is the identity
  in each schedule's independently sorted quotient coordinates.
- **DERIVED:** each schedule has 24 phase-reversing `H4` actions, all inducing
  one quotient reversal permutation.
- **DERIVED:** the two static complete-action Hessian singular spectra differ
  at resolved order `1e-2`; the two schedules are not merely the same stored
  Hessian under coordinate permutation.
- **DERIVED COMPUTATIONAL:** the dynamic canonical shape spectra agree at 160
  digits within the calibrated uncertainty, while the maps are nonnormal.
- **OPEN:** an exact carrier relation explaining that agreement.

## Framing attack

Any two diagonalizable matrices with the same spectrum are similar over the
complex numbers.  Solving for an unrestricted intertwiner after observing
isospectrality would therefore be content-free.  Even a successful boundary
permutation is weaker than a slab isomorphism: it can relate the endpoint
coordinates without relating the bulk triangulations that generated the
maps.

The next calculation must consequently enumerate the finite geometric family
before reading a single intertwining residual, report its exact size and hit
fraction, and maintain separate complete-slab and boundary-only verdicts.

## Status

- **KNOWN:** action-generated canonical covariance under an actual simplicial
  isomorphism.
- **CONTROL:** finite `H4` carrier, two committed schedules, boundary orbit
  partitions, reversal maps and calibrated tangent matrices.
- **OPEN:** whether any complete-slab cross-parity isomorphism exists.
- **OPEN:** whether any finite carrier-derived lift intertwines the maps.
- **OPEN:** whether the numerical isospectrality has an exact algebraic or
  geometric explanation.
