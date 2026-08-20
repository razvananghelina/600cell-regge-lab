# Prior-art gate: feasibility of a refined unrestricted canonical map

Date: 2026-08-20

Status: completed before constructing or diagonalizing any refined action
Hessian.

## 1. Exact object and complete hypotheses

Use only the two certified spatial carriers

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)),
P(x)=x/||x||,
```

their certified proper four-colourings, and the standard colour-consecutive
staircase triangulations of `Ki x I`.  A temporal schedule is one of the 24
linear orders of the four colour classes.  The time orientation is fixed, so
an order and its reverse remain distinct labelled slabs.

For a chosen schedule, the proposed unrestricted Regge carrier consists of
all old-boundary, new-boundary and slab-interior edge squares.  The desired
eventual object is the action-generated boundary canonical/Jacobi map after
eliminating all slab-interior edge variables at a homogeneous Lorentzian
background with the already certified local `P1` dust weights.

This gate does **not** choose one schedule, average schedules, compare a
spectrum with a continuum target, or define a coarse-to-fine mode matching.
It asks whether the finite problem is sufficiently selected and sufficiently
small to justify constructing that operator.

## 2. KNOWN

The relevant constructions and obstructions are standard in broad form.

- Joswig and Witte construct simplicial products from linear vertex orders
  and show that the order matters; colour-consecutive orders give balanced
  products but are not selected by the product itself:
  <https://arxiv.org/abs/math/0508180>.
- Dittrich and Hoehn obtain canonical simplicial evolution from the action as
  Hamilton's principal function:
  <https://arxiv.org/abs/1108.1974>.
- Four-dimensional Regge dynamics is not generally triangulation independent.
  Linearized actions and measures under Pachner moves are analyzed in
  <https://arxiv.org/abs/1110.6866> and
  <https://arxiv.org/abs/1404.5288>.
- Improved or perfect actions require dynamical coarse graining; a common
  refinement alone does not produce one:
  <https://arxiv.org/abs/0907.4323>.

Therefore a schedule-independent refined canonical map is not supplied by a
general theorem.

## 3. Repository controls

The following results are already certified and are not recomputed as new
physics:

- both projected spatial carriers are closed, choice-free and `H4`
  equivariant;
- each has a proper four-colouring unique up to global colour renaming;
- all 24 colour orders give conforming slabs and all 24 are compatible with
  the tested spatial symmetry;
- no colour order is selected by the spatial carrier plus time orientation;
- the direct cellular homogeneous action is schedule independent because all
  schedules triangulate one flat homothetic Lorentzian frustum;
- this flat-frustum theorem does not extend automatically to arbitrary
  anisotropic edge perturbations;
- the nodal `P1` dust weights are unique only conditional on the structural
  `P1` matter ansatz.

## 4. Framing attack

The phrase "refine the unrestricted canonical map" hides three logically
separate choices.

1. **Temporal carrier.**  An unrestricted simplicial Hessian is attached to
   one of 24 distinct internal diagonal sets.  Homogeneous action equality
   does not prove equality of their anisotropic quadratic actions.
2. **Elimination.**  Comparing raw Hessians is meaningless because their
   internal edge labels differ.  The candidate invariant object is the
   effective boundary canonical map after eliminating the internal equations.
3. **Inter-level transport.**  Even if each level supplies a unique boundary
   map, their phase spaces have different dimensions.  A prolongation or
   restriction between coarse and fine metric edge variations must be derived
   before matching modes.  The coordinate radial projection used to construct
   the background is not automatically an intrinsic phase-space map.

A full dense matrix is also not required by the physics.  If it is too large,
a matrix-free low-mode computation may remain possible.  Consequently memory
size can close a dense implementation, but not the scientific route.

## 5. CONTROL / OPEN / FORBIDDEN

- **CONTROL:** reconstruct the two spatial f-vectors and the 24 conforming
  schedule count from frozen sources.
- **OPEN:** exact total, internal and boundary edge counts; distinct schedule
  internal-edge sets; direct dense storage sizes; local sparse-support upper
  bounds.
- **OPEN:** equality of all 24 effective quadratic boundary maps.
- **OPEN:** a geometry-selected coarse/fine phase-space transport.
- **FORBIDDEN:** choosing the numerically most convenient order, averaging
  schedules, fitting a transport from eigenvectors, or reading physical modes
  before these choices are settled.

## 6. Next gate

Preregister an exact combinatorial feasibility census.  It may license only a
later all-schedule quadratic boundary-covariance test or return a structural
blocker.  It cannot establish gravitons, a dispersion relation, an effective
speed, a physical tick, `G` or Planck scales.

External novelty remains **OPEN**.  A targeted source search cannot prove it.
