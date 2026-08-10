# Preregistration: shape-controlled smooth Hopf red refinement

Date: 2026-08-10

## Scope

The first projected-barycentric refinement passed its continuum-mode gates,
but one level cannot establish a convergent tower and repeated barycentric
subdivision has a known shape-degeneration obstruction.  This protocol fixes
a different refinement before any spectra at its new levels are inspected.

The question is purely numerical geometry:

> Does a standard projected tetrahedral `1 -> 8` refinement preserve element
> quality through two levels while the canonical smooth-Hopf modes converge
> toward the exact round-`S3` spectrum?

No bootstrap integer, speed, mass, or target involving `a_1` is admitted.

## Refinement rule fixed before execution

Starting from the 600 tetrahedra of the 600-cell boundary:

1. create one globally shared midpoint for every mesh edge;
2. radially project that midpoint to the unit `S3`;
3. split each parent into four corner tetrahedra plus the central octahedron;
4. the octahedron has three possible opposite-midpoint diagonals; choose the
   shortest Euclidean chord, breaking an exact tie lexicographically by the
   global endpoint pair;
5. triangulate the octahedron into four tetrahedra around that diagonal.

The choice is local geometric mesh quality, not spectral information.  It
does not inspect Hopf modes or continuum eigenvalues.

The expected combinatorial sizes are fixed independently:

```text
level      vertices      tetrahedra
0               120               600
1               840              4800
2              6480             38400
```

Every triangular face must occur twice, since the carrier is a closed
three-manifold.

## Shape gate fixed before execution

For a chordal tetrahedron with Euclidean three-volume `V` and six squared edge
lengths, use the mean-ratio quality

`q = 12*(3V)^(2/3) / sum_edges(length^2)`.

A regular tetrahedron has `q=1`; a collapsing tetrahedron has `q->0`.

The preregistered first-two-level gate is:

- maximum chord length decreases strictly at each level;
- minimum `q` at levels 1 and 2 is at least `0.5`;
- `q_min(level2) >= 0.8*q_min(level1)`.

This is a finite two-level shape gate, not a proof of uniform shape regularity
for every future level.

## Hopf operator and canonical calibration modes

Use exactly the smooth orthogonal construction already frozen in commit
`5f78826`:

`X(q)=q*u`, `P_V=X_t tensor X_t`, `P_H=P_tangent-P_V`.

The assembled identity `K_V+K_H=K_full` must hold on every level.  The
projected Hopf vector must have norm greater than `0.98` before normalization
in every element.

The comparison spaces and exact continuum values are unchanged:

- four ambient coordinate modes: `(lambda_V,lambda_H,lambda_full)=(1,2,3)`;
- three Hopf-base coordinate pullbacks: `(0,8,8)`.

For every component and every calibration space, the maximum absolute target
error must decrease strictly from level 0 to 1 and again from level 1 to 2.
The combined low spectrum must have exactly one resolved zero at all levels.

Low separated spectra are recorded blindly.  No multiplicity clustering or
near-zero threshold beyond the explicit calibration spaces is an acceptance
criterion.

## Interpretation boundary

Passing would establish a robust two-level numerical approximation of the
round Hopf geometry.  It would not prove an infinite convergence theorem, and
it would not select an anisotropy coefficient `r` in `K_H+r K_V`.

Failure would kill this projected red-refinement implementation.  It would not
refute the analytic Hopf spectrum.

