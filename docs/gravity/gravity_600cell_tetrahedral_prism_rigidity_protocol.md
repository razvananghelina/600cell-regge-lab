# Protocol: constrained rigidity of a tetrahedral time prism

Date: 2026-08-19

Prior-art commit: `dd181a3`.
Pre-protocol framing-correction commit: `018cb5c`.

## Scope and complete hypotheses

The local cell is the combinatorial four-polytope `Delta_3 x I`.  Its input
metric data are the 16 squared lengths of the six bottom edges, six top edges
and four corresponding struts.  In addition, being this polytope means that
each of its six lateral four-cycles is a planar quadrilateral.  These
planarity constraints are part of the hypothesis, not fitted measurements.

The test is local and exact.  It does not construct a 600-cell slab, choose a
staircase, evaluate a Regge action or compare with continuum gravity.

## Exact control family

Work first in rational Euclidean coordinates in `R^4`.  The bottom
tetrahedron is

```text
b0=( 1,0,0,0),  b1=(-1,0,0,0),
b2=( 0,1,0,0),  b3=( 0,0,1,0).
```

For scale ratio `q` and translation `t`, set

```text
ui = q*bi + t.
```

The preregistered scale ratios are

```text
q in {1, 9/10, 11/10, 2},       t=(0,0,0,2).
```

All six lateral faces are planar trapezoids.  Because all `bi` have the same
Euclidean spatial norm, all four struts also have a common length for each
`q`.  The same coordinate family will be tested with the Lorentzian bilinear
form `diag(1,1,1,-1)`.

To attack symmetry dependence, a second exact control is obtained by a
rational non-affine projective transformation of the standard prism.  A
projective map preserves lines, planes and the face lattice but removes the
parallel/homothetic symmetry.  Its explicit integer matrix, translation and
rational denominator covector will be frozen in the verifier.

## Constraint Jacobians

For eight labelled vertices in four dimensions, form:

1. `J_edge`, the Jacobian of the 16 natural squared lengths;
2. `J_plane`, the Jacobian of all four `3x3` minors of the three difference
   vectors on each of the six lateral quadrilaterals;
3. `J_poly`, the vertical concatenation of the two.

Although 24 minors are written, their rank rather than their row count is the
quantity recorded.  Euclidean rigid motions contribute ten kernel
dimensions.  The quotient infinitesimal-flex count is

```text
32-rank(J_poly)-10.
```

For Lorentzian edge lengths, replace the Euclidean quadratic form by
`diag(1,1,1,-1)` and retain the same affine planarity equations.

The scratch values that forced the framing correction are disclosed: rank 19
at `q=1` and rank 22 at tested `q!=1`.  They are not counted as evidence.  The
registered verifier must reconstruct every symbolic Jacobian from scratch.

## Finite nonuniqueness attack at equal scale

Infinitesimal rank loss at a symmetric point need not integrate to a finite
flex.  Therefore the verifier will construct explicit pairs.

At `q=1`, translating the entire top tetrahedron by any vector `t` preserves
top and bottom edge lengths and makes every lateral face a parallelogram.  In
Euclidean signature compare

```text
t=(1,2,3,4),       t'=(4,2,3,1),       |t|^2=|t'|^2=30.
```

In Lorentzian signature compare

```text
t=(1,2,3,4),       t'=(0,1,1,2),
<t,t>=<t',t'>=-2.
```

For each pair, verify exactly:

1. all 16 natural squared lengths agree;
2. all lateral quadrilaterals are planar;
3. the absolute four-volumes differ.

Different volume proves that the two labelled cells are not related by an
ambient isometry, so this is finite metric nonuniqueness rather than only an
infinitesimal mechanism.

## Diagonal controls and relation to the 24 schedules

Each of the six lateral quadrilaterals has two cross-diagonals.  Enumerate all
`2^6=64` ways of choosing one per face.  On the equal-scale rational control:

1. record the exact edge-rigidity rank after the six diagonal lengths are
   added;
2. delete each added diagonal in turn and record the rank;
3. identify the 24 transitive choices corresponding to staircase orders and
   compare them with the 40 cyclic choices.

This diagnoses whether the number six is a rigidity count or specifically a
staircase phenomenon.  No cyclic choice is accepted as a triangulation.

## Frozen interpretation

- Edge-only rank 16 with six graph flexes is a control, not the polytopal
  verdict.
- Rank 22 for `q!=1` establishes only local infinitesimal determination on
  those exact controls.
- A finite same-length/different-volume pair at `q=1` proves genuine
  schedule-free underdetermination at the static equal-scale cell.
- The three equal-scale flexes are to be compared with tangential translation
  of the top tetrahedron.  Their physical interpretation as discrete shift
  variables is **STRUCTURAL/OPEN**, not automatically gauge.
- If six diagonals restore rank, this shows what extra metric data a
  simplicial schedule supplies.  It does not select one schedule.

## Decision boundary

The edge-count no-go is refuted if planarity raises the generic unequal-scale
controls to rank 22.  The schedule-free static cell nevertheless fails if the
finite equal-scale pairs have identical declared lengths and different
volumes.

In that mixed outcome, no Hessian based only on the 16 natural lengths may be
called canonical at `q=1`.  The next admissible question is whether a shift
variable or a derived normal-evolution condition supplies the missing three
data consistently across the full carrier.

## Artifact

The registered verifier will be

```text
reproducible/verify_gravity_600cell_tetrahedral_prism_rigidity.py
```

and will write

```text
reproducible/gravity_600cell_tetrahedral_prism_rigidity.json.
```

Only this verifier and static guards are required.  The full suite is
excluded by user instruction.

