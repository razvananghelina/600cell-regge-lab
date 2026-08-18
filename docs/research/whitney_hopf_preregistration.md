# STEP 1 preregistration: blind Whitney/Hopf kinetic enumeration

Date: 2026-08-10

## Protocol declaration

This file and the accompanying JSON are committed **before** comparing any
observed Whitney spectral ratio with a bootstrap integer or a proposed speed.
No target comparison is performed in this step.

The committed artifacts are:

- `reproducible/verify_whitney_hopf_blind_enumeration.py`;
- `reproducible/whitney_hopf_blind_enumeration.json`;
- this protocol note.

## Complete hypotheses

1. The carrier is the unrefined piecewise-Euclidean boundary complex of the
   regular 600-cell.
2. The scalar Hilbert space is the lowest-order Whitney zero-form space.
3. `M_0` is the exact consistent `L2` mass assembled from all 600 Euclidean
   tetrahedral facets.
4. `K_0=d^T M_1 d` is the exact scalar Whitney stiffness.
5. A discrete Hopf fibration is one of the six partitions into twelve
   order-ten cosets already derived in the repository.
6. A coarse edge is `fiber` precisely when its endpoints lie in the same
   Hopf fiber; every other 600-cell edge is `cross`.
7. Since the assembled regular-facet stiffness has one common edge weight,
   its fiber/cross split is fixed by this edge support.  Component diagonals
   are fixed by zero row sums.
8. Generalized spectra are computed for `(K_fiber,M_0)` and
   `(K_cross,M_0)`.  Raw stiffness spectra are also recorded so that the
   influence of the mass matrix cannot be hidden.
9. Every distinct eigenvalue, its multiplicity, both first positive gaps,
   their ratio, kernel multiplicities, edge counts and the per-tetrahedron
   fiber-edge histogram are written to JSON for all six fibrations.

No level coefficient is varied and no Schur block is fitted.

## Scope boundary fixed before comparison

The fiber/cross split is canonical on the coarse 600-cell edge set.  A
canonical continuation of that binary split to the new barycentric edges has
not yet been defined.  The enumeration therefore makes no refinement-stability
claim.  The per-tetrahedron fiber-edge histogram is preregistered because it
constrains whether such a continuation can be unique.

## Evidence labels before comparison

- **DERIVED:** the Whitney matrices and all six Hopf partitions are fixed by
  the stated geometry.
- **DERIVED:** the complete raw and generalized spectral multisets can be
  enumerated without a physical target.
- **OPEN:** whether any observed ratio agrees with an independently derived
  bootstrap datum.
- **OPEN:** whether the split extends canonically under refinement.
- **NOT CLAIMED:** a speed of light, Lorentzian time or physical units.
