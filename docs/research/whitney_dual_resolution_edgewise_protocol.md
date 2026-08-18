# Preregistration: dual-resolution locality on the rank-edgewise tower

Date: 2026-08-12

This protocol is committed before enumerating any new vertex or edge link on
the `k=2,4` controls.  The preceding dual-resolution result used the base
600-cell and one barycentric subdivision; its observed degree growth cannot
be transferred silently to the distinct rank-edgewise tower.

No spectrum, target value or phenomenological quantity is used.

## Frozen carrier

Use exactly the already-preregistered shape-regular controls

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,4,
\]

with the canonical flag rank order and exact vertex identifications supplied
by `whitney_trace_refinement_tools.py`.

The complete dual-cell constraint resolution is not re-selected here.  Its
incidence definition and exactness theorem are inherited from protocol
`c5f9bee` and result commit `799966f`.  This audit asks whether its incidence
degrees behave differently on the actual refinement mechanism.

## Complete hypotheses to transfer exactness

For every level verify combinatorially:

1. every triangle has exactly two parent tetrahedra;
2. for every edge, its incident tetrahedra and triangles form one connected
   cycle;
3. for every vertex, its link is connected, every link edge has two link
   faces, and its Euler characteristic is exactly two;
4. every tetrahedron occurrence graph for a fixed simplex is connected.

These are the finite link hypotheses used by the dual-cell proof:

- a triangle has a dual interval;
- an edge has a dual disk whose boundary is its link circle;
- a vertex has a dual 3-ball whose boundary is its link 2-sphere.

Passing them transfers the exact signed resolution to each frozen control.
It does not prove an arbitrary-mesh theorem.

## Frozen locality observables

Record complete histograms and maxima of:

- `a_0(k)`: tetrahedron occurrences at one vertex, the maximum row degree of
  the quotient `A_0`;
- `a_1(k)`: tetrahedron occurrences at one edge, simultaneously the maximum
  row degree of `A_1` and the maximum boundary length in `R_0,2` and
  `R_1,2`;
- `r_3(k)`: edges incident at one vertex, the maximum boundary length in
  `R_0,3`;
- the already-proved neighbour constraint node maxima `(3,2,1,0)`.

All quantities are exact integer incidences.  No support tolerance is used.

## Frozen flow gate

The coarse `k=1` level may have exceptional symmetry.  The decision is based
on both transitions and is fixed before seeing the values:

- **PATTERN TOWARD BOUNDED DUAL LOCALITY** if every observable is
  nonincreasing from `k=2` to `k=4`; equality or decrease is accepted, but
  three levels do not prove a uniform all-level bound;
- **PATTERN NEGATIVE FOR BOUNDED DUAL LOCALITY** if any observable increases
  from `k=2` to `k=4`;
- separately record every increase from `k=1` to `k=2`; it cannot be hidden
  by the primary gate;
- an exact all-level boundedness claim requires a combinatorial proof for
  general `k` and remains **OPEN** regardless of the finite outcome.

No candidate map or subfamily may be dropped because its degree is larger.

## Attack on the framing

Bounded incidence of this kinematic resolution would only show that the
canonical redundant constraint/ghost hierarchy can be stored and applied by
a bounded local stencil.  It would not make the original constraints first
class, derive a BRST charge with physical cohomology, or remove the Gram
inverse from gauge-invariant observables.

Conversely, degree growth on three levels would be a negative pattern rather
than a theorem of unbounded growth.  The hierarchy could still be subdivided
or represented by a different local resolution, but any such construction
must be specified independently of a desired spectrum.

## Outputs and exclusions

The registered verifier will write
`reproducible/whitney_dual_resolution_edgewise.json` and a result note will
record the protocol commit, every histogram, the exact link gates and the
finite-flow label.

Excluded:

- no operator spectrum;
- no metric, Hamiltonian or physical gauge claim;
- no continuum extrapolation or fitted bound;
- no time, causality, inertia, mass or Planck units;
- no full-suite run, by explicit user request.
