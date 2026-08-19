# Prior-art gate: rigidity of a schedule-free tetrahedral time prism

Date: 2026-08-19

## Exact object and question

Let `T=Delta_3 x I` be one four-dimensional time cell with eight vertices:
four on a bottom tetrahedron and four corresponding vertices on a top
tetrahedron.  Its natural one-skeleton contains

```text
6 bottom edges + 6 top edges + 4 struts = 16 edges.
```

Assume:

1. the cell is a flat nondegenerate four-dimensional Euclidean or Lorentzian
   polytope;
2. its metric variables are only the squared lengths of those 16 natural
   edges;
3. adjacent cells shape-match on every shared three-face;
4. no cross-diagonal, normal, area-angle datum, embedding coordinate or
   staircase order is supplied.

The question is whether these data determine the cell volume and four-
dimensional dihedral angles up to isometry, as required by a local polytopal
Regge action.  If not, how many independent metric constraints are missing,
and do the six diagonals introduced by a staircase triangulation supply
exactly them?

This is a local rigidity question.  If one cell is underdetermined, summing
over the complete 600-cell carrier cannot make its local dihedral angles a
function of the declared edge lengths without additional variables.

## Known rigidity criterion

For a framework with `n` vertices in `d` dimensions, the rigidity matrix is
the Jacobian of the squared edge-length map.  A nondegenerate infinitesimally
rigid framework has rank

```text
d*n - d*(d+1)/2,
```

where the subtracted term is the dimension of translations and rotations.
This is the standard Asimow--Roth framework:

- L. Asimow and B. Roth, *The Rigidity of Graphs, II*, Journal of
  Mathematical Analysis and Applications 68 (1979), 171--190, DOI
  `10.1016/0022-247X(79)90108-2`.

For `n=8,d=4`, full infinitesimal rigidity requires rank 22.  A 16-row
rigidity matrix has rank at most 16, so the natural prism graph has at least
six non-isometric infinitesimal motions.  This necessary-count negative is
already analytic; computation will test that the bound is sharp on an exact
flat frustum and that the staircase diagonals remove it.

The same constraint count applies in nondegenerate pseudo-Euclidean
signature.  Generic global rigidity has also been developed explicitly for
pseudo-Euclidean spaces, but the present negative needs only the row-count
bound: changing the nondegenerate bilinear form cannot turn 16 scalar length
constraints into rank 22.

## Regge and hyperfrustum literature

Standard length Regge calculus uses simplices, for which edge lengths fix the
flat metric.  The directly relevant cosmology literature also exposes the
extra-data issue.

Tsuda and Fujiwara use four-dimensional frusta as homogeneous building
blocks in the Collins--Williams formalism.  Their regular ansatz assumes
parallel regular top and bottom cells and equal struts; polygon and trapezoid
areas, dihedral angles and volumes then have explicit formulae:

- R. Tsuda and T. Fujiwara, *Oscillating 4-Polytopal Universe in Regge
  Calculus*, arXiv:`2011.04120`, DOI `10.1093/ptep/ptab079`.

That is a symmetry-reduced construction, not a proof that 16 arbitrary local
edge lengths determine a general flat prism.

Dittrich, Gielen and Schander study tetrahedral frusta and state that flatness
plus isosceles trapezoidal faces fixes their homogeneous building block.
They then subdivide it into four 4-simplices by adding six face diagonals;
the subdivision permits an ordinary Regge action.  Setting the diagonal
lengths to their flat-frustum values produces a generalized frustum action:

- B. Dittrich, S. Gielen and S. Schander, *Lorentzian quantum cosmology goes
  simplicial*, arXiv:`2109.00875v3`, Section V.A, especially pp. 23--25.

Donà, Fanizza, Sarno and Speziale discuss Regge actions for general 4D
polytopes.  Flat polytope data require closure, shape matching and
flat-embedding conditions; more general area-angle data allow conformal
shape mismatch and are not ordinary Regge geometry.  They explicitly
distinguish polytope data from edge-length Regge data and note where
subdivision is needed to represent curvature:

- P. Donà, M. Fanizza, G. Sarno and S. Speziale, *SU(2) graph invariants,
  Regge actions and polytopes*, arXiv:`1708.01727`, Section 5.

Thus a schedule-free polytopal action is possible only if enough flat-polytope
data are supplied.  Whether the present theory derives those extra data is
the question; their mere availability is not selection.

## Repository relation

The preceding mission found exactly 24 distinct staircase schedules.  Each
schedule chooses one diagonal on each of the six quadrilateral faces
`edge x I`.  The numerical coincidence

```text
22 required constraints - 16 natural edges = 6 missing constraints
```

was noticed only after that result.  It is therefore **PATTERN** before the
registered rigidity test, not evidence.

The old homogeneous cellular-frustum action remains a valid symmetry-reduced
control.  It supplies parallelism, common scale and a common lapse by ansatz.
It does not falsify a six-flex local no-go.

## Classification before execution

- **KNOWN:** rigidity-matrix rank criterion; standard simplicial Regge data;
  homogeneous hyperfrustum actions with extra symmetry/flatness assumptions.
- **DERIVED ANALYTIC:** 16 natural edge constraints cannot reach the required
  rank 22; at least six infinitesimal non-isometric motions remain under the
  stated length-only hypotheses.
- **CONTROL:** an exact rational homothetic frustum must realize a flat
  nondegenerate `Delta_3 x I` cell.
- **OPEN:** whether its natural 16 rows have full row rank 16; whether every
  transitive six-diagonal staircase completion has exact rank 22; whether
  nontransitive diagonal choices behave differently.
- **OPEN:** a theory-derived set of normals, areas or embedding constraints
  that could replace the six diagonals without fitting.

## Decision boundary

If the natural prism has six flexes and each staircase completion removes
them, a schedule-free edge-length polytopal action is not locally determined.
The temporal route remains **STRUCTURAL** unless new polytope variables and
their own selection law are derived.

If the exact rank does not show this, the count-only framing was insufficient
and the result stays **OPEN**.  No local Regge Hessian may be evaluated from
underdetermined dihedral angles.

