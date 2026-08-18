# Prior-art and framing gate: transported negative-stiffness intersection

Date: 2026-08-18

Status: **TARGET-DISCLOSED, NO TRANSPORTED-INTERSECTION SPECTRUM COMPUTED.**

## Exact question

On each of the already frozen `16` cells (two parities, sectors `4,5`, four
derivative schedules), let `E^-_0,E^-_1` be the old and shifted rank-`15`
negative spectral subspaces of the Euclidean shape-stiffness form `-V_S`.
Let

```text
F^-_t = E^-_t direct-sum conjugate(E^-_t)
```

inside the `60`-dimensional boundary phase carrier, and let `T_2` be the
already certified second-slab tangent.  The question is the exact dimension
of

```text
K^- = F^-_0 intersection T_2^(-1)(F^-_1).
```

Equivalently, for the orthogonal phase projectors `Q_0,Q_1`, compute the rank
of

```text
R^- = (I-Q_1) T_2 Q_0.
```

No desired rank, graph or physical interpretation is used to choose the
operator.

## Primary literature checked before the protocol

1. Dittrich and Hoehn, *Canonical simplicial gravity*, Class. Quantum Grav.
   **29** (2012) 115009, DOI
   `10.1088/0264-9381/29/11/115009`, arXiv:`1108.1974`.  The discrete action is
   the generating function for canonical evolution; in singular systems the
   map is defined on a pre-constraint surface, has a post-constraint surface
   as image, and preserves the canonical two-form only after restriction.
2. Dittrich and Hoehn, *Constraint analysis for variational discrete
   systems*, J. Math. Phys. **54** (2013) 093505, DOI
   `10.1063/1.4818895`, arXiv:`1303.4294`.  Propagating observables and reduced
   phase space depend on both the initial and final evolution steps; equality
   of a carrier dimension at one time is not a propagation theorem.
3. Dittrich and Hoehn, *From covariant to canonical formulations of discrete
   gravity*, Class. Quantum Grav. **27** (2010) 155001, DOI
   `10.1088/0264-9381/27/15/155001`, arXiv:`0912.1817`.  Exact gauge symmetry
   is special to linearized flat Regge backgrounds; curvature generically
   produces background-dependent pseudo-constraints.

None of these papers studies the binary-icosahedral `600`-cell dust carrier or
the particular subspace `K^-`.  External novelty remains **OPEN** pending a
dedicated review.

## Framing attack

### What is intrinsic

The inertia `(15 negative, 10 positive)` of the nondegenerate Hermitian form
`-V_S` on the shape quotient is invariant under congruence.  The shape
carrier itself is selected by the action through the `M`-orthogonal removal
of the derived conformal image.  These are **DERIVED COMPUTATIONAL** facts of
the frozen model.

### What requires extra structure

An indefinite form alone does not select a unique maximal negative subspace.
The spectral space `E^-` is selected only after using the positive Hermitian
inner product inherited from the orthonormal binary-action carrier.  That
inner product is geometry-fixed and was committed upstream, so it introduces
no adjustable coefficient, alignment or target search.  It is nevertheless
additional to the Regge quadratic form and symplectic evolution.

Therefore:

- the full census below is canonical **relative to the frozen carrier
  Hilbert metric**;
- a nonzero `K^-` would be **STRUCTURAL**, not yet a physical propagating
  mode or a pre/post constraint surface in the sense of the cited canonical
  formalism;
- the generalized kinetic--stiffness fiber was dynamically better motivated,
  and its transported intersection has already been certified zero;
- a zero `K^-` would close this remaining Euclidean negative-fiber phase
  route, but not full Regge dynamics.

This qualification is load-bearing.  Calling `E^-` "action-selected" without
also naming the carrier metric would overstate canonicity.

## KNOWN / CONTROL / OPEN

### KNOWN

- The old and shifted Euclidean negative projectors have rank `15` in all
  `32` schedule-local cases at binary precision, with a separated `15/10`
  gap.
- Their configuration subspaces agree within the existing conservative
  enclosures in all `16` matched cells.
- The unrestricted lift `F^-_0` is not carried wholly into `F^-_1`: the
  `B,D` leakage blocks are nonzero-resolved.  This does not determine the
  kernel of the full leakage operator.
- Exact Flint source balls for `M,V` and exact second-slab tangent balls have
  already been reconstructed independently of the binary archive in the
  accepted generalized-fiber chain.

### CONTROL

- Reconstruct `E^-_t` anew from those high-precision source balls, using the
  frozen carrier metric and no generalized-pencil eigenvectors.
- Require all reconstructed projectors to overlap the earlier binary
  projectors within their committed conservative errors.
- Require a certified rank-`15` spectral split and propagate source,
  shape-carrier and eigenspace errors.
- Use all `16` cells and all `60` singular values per cell.  The exact right
  factor `Q_0` gives `rank(R^-) <= 30`.

### OPEN

- Whether `rank(R^-)=30` (zero transported intersection), is smaller, or
  cannot be resolved.
- If a nonzero intersection is certified, whether it is a graph over
  configuration space, Lagrangian, symplectically transported or physical.
- Constraint reduction, long-time persistence, refinement, continuum
  behavior, graviton interpretation, inertia, mass and limiting speed.

## No-go conditions before implementation

The route must remain **OPEN** rather than claim a positive dimension if the
calculation resolves fewer than `30` nonzero singular values without proving
an independent exact rank upper bound.  "Zero-consistent" singular values are
not exact zeros.  Likewise, a convenient `15+15` visual split is not evidence
for a `15`-dimensional intersection.

No Procrustes alignment, basis optimization, schedule deletion, tolerance
tuning after inspection or fitted Lagrangian graph is admissible.
