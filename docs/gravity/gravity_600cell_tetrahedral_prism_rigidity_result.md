# Result: the missing temporal datum is a three-component shift

Date: 2026-08-19

## Headline

The naive claim that a schedule-free tetrahedral prism must have six metric
flexes was wrong: it counted only its 16 natural edges and forgot that a
four-polytope also requires its six lateral quadrilaterals to be planar.  That
error was found and committed before the protocol.

With the complete hypothesis, the exact result is sharper:

```text
unequal-scale or generic prism:  locally infinitesimally determined;
equal-scale prism:               exactly 3 non-isometric modes.
```

The three modes translate the whole top tetrahedron tangentially relative to
the bottom one.  Finite Euclidean and Lorentzian examples have identical 16
natural squared lengths but different four-volumes.  Kinematically these are
precisely shift-like data.  Whether they are gauge in the discrete theory is
**OPEN**, not assumed.

This mixed result is **DERIVED EXACT, ADVERSARIALLY CORROBORATED**.  It blocks
a length-only static Hessian, but it opens a better next route: add or derive
the shift/normal data on the full carrier instead of choosing one of 24
staircase triangulations.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art gate | `dd181a3` |
| pre-protocol framing correction | `018cb5c` |
| primary protocol | `360d2a5` |
| registered exact-minor verifier | `8825c76` |
| frozen primary artifact | `f5d7c9b` |
| adversarial protocol | `e424c98` |
| registered projective-coordinate audit | `68c4175` |
| frozen adversarial artifact | `a26af8b` |

The primary artifact was reproduced byte-for-byte with SHA-256

```text
ce9eb1917dd647c6dd8155a0f9646a72dc7734c0310f763ec31e070403230db8.
```

The adversarial artifact was reproduced byte-for-byte with SHA-256

```text
511c9d5f6747357e42e8299aa38f6cadea70525d86ae8d35b4d2f105ed199689.
```

## 1. The framing error and its correction

A framework of eight vertices in four dimensions needs rigidity rank

```text
4*8 - 4*5/2 = 22.
```

The tetrahedral-prism graph has only

```text
6 bottom + 6 top + 4 struts = 16 edges,
```

so its edge-rigidity matrix cannot have rank above 16.  It was tempting to
call the difference of six a no-go and identify it with the six staircase
diagonals.

That argument applies to a graph framework, not yet to a polytope.  A
realization of `Delta_3 x I` has six planar quadrilateral side faces.  The
Jacobian of their determinant constraints has exact rank eight on every
registered control.  Combined with edge lengths, those constraints can raise
the rank to 22 without any diagonal.

This was caught before preregistration.  The git history therefore preserves
both the bad initial framing and its correction rather than silently editing
the hypothesis after the result.

## 2. Exact constrained ranks

The bottom tetrahedron was fixed in rational coordinates, and the top was

```text
u_i = q b_i + (0,0,0,2)
```

for `q=1,9/10,11/10,2`.  All bottom vertices have equal spatial norm, so the
four struts have equal length at each `q`.  Every lateral face is planar.

In both Euclidean signature and Lorentzian signature `diag(1,1,1,-1)`:

| scale ratio `q` | edge rank | planarity rank | combined rank | non-isometric infinitesimal modes |
|---:|---:|---:|---:|---:|
| `1` | 16 | 8 | 19 | 3 |
| `9/10` | 16 | 8 | 22 | 0 |
| `11/10` | 16 | 8 | 22 | 0 |
| `2` | 16 | 8 | 22 | 0 |

A rational non-affine projective image of the standard prism also has
combined rank 22 in both signatures.  Thus the positive unequal-scale result
is not an artifact of parallel regular coordinates.

The conclusion is local and infinitesimal.  It does not assert a global
single-valued inverse for arbitrary length data or global compatibility on
the complete 600-cell carrier.

## 3. The equal-scale modes are finite, not numerical zero modes

After pinning the bottom tetrahedron, the equal-scale constraint Jacobian has
top-vertex nullity three.  Every basis mode gives the same velocity to all
four top vertices: it translates the top rigidly in a direction tangent to
the fixed-strut-length sphere or hyperboloid.

This mechanism integrates exactly.  At `q=1`, let the top be a translated
copy of the bottom.  Top and bottom edges are independent of the translation,
all side faces are parallelograms, and the four struts all have the squared
norm of the translation vector.

The registered finite pairs are:

```text
Euclidean:
    t =(1,2,3,4),  t'=(4,2,3,1),  |t|^2=|t'|^2=30
    volumes: 4/3 versus 1/3

Lorentzian:
    t =(1,2,3,4),  t'=(0,1,1,2),
    <t,t>=<t',t'>=-2
    absolute volumes: 4/3 versus 2/3.
```

All 16 natural squared lengths agree exactly in each pair.  Ten of the twelve
cross-diagonal squared lengths differ, and the four-volumes differ.  Hence
the cells are genuinely non-isometric; their dihedral/volume data cannot be
functions of the 16 natural lengths alone.

This is **DERIVED EXACT NEGATIVE** for a schedule-free equal-scale length-only
cell.

## 4. Independent projective-coordinate audit

The audit used no face-planarity minors.  Instead, it parameterized the local
realization space by the 24 infinitesimal directions of `PGL(5)` acting on
homogeneous vertex coordinates.  This is justified by the standard fact that
products of simplices are projectively unique; `Delta_3 x Delta_1` is such a
product.

The rank of the 16-length map on that projective chart is:

```text
q=1:                  11, kernel 13 = 10 isometries + 3 modes;
q=9/10,11/10,2:       14, kernel 10 = 10 isometries only.
```

Euclidean and Lorentzian ranks agree exactly.  The audit also differentiated
three explicit tangential translations and found zero derivative for all 16
lengths, then independently reconstructed the finite pairs and the missing
cross-diagonal data.

This mechanically different route corroborates the complete mixed result.

## 5. What the six diagonals do—and do not do

At the equal-scale cell, every one of the `2^6=64` choices of one diagonal on
each lateral quadrilateral raises the edge-rigidity rank to 22 in both
signatures.  Removing any one of the six added diagonals drops the rank to 21
in all `64*6=384` tests.

Thus six diagonal lengths are an exact isostatic completion of the natural
edge data at this cell.  But rigidity does not distinguish:

```text
24 transitive choices that come from staircase orders,
40 cyclic choices that are not staircase triangulations.
```

All 64 have rank 22.  Therefore rigidity explains why six additional metric
numbers suffice; it does not select a temporal triangulation and does not
repair the preceding 24-schedule no-go.

## 6. Physical reading

In a `3+1` decomposition, a lapse gives normal separation while a
three-component shift gives tangential displacement between slices.  The
mathematics here has exactly that kinematics:

- fixed strut norm fixes a proper-time-like magnitude;
- the three undetermined directions move the top tangentially;
- changing the scale away from one makes the four struts geometrically
  distinct enough to determine the relative placement locally.

This match is **STRUCTURAL**, not yet a derivation of ADM gauge symmetry.  A
discretization can break diffeomorphism symmetry, and a shift-like variable
need not be a null gauge direction after cells are glued and matter is added.
The full-carrier constraint rank must decide that.

For the homogeneous `H4` background, a nonzero invariant tangent vector is
not expected, so zero shift is compatible with symmetry.  That does not
select local shifts once symmetry is perturbed.

## 7. Literature reconciliation

The post-result search located a precise justification for the adversarial
parameterization: all products of simplices are projectively unique.

- J. Gouveia, A. Macchia, R. R. Thomas and A. Wiebe, *Projectively unique
  polytopes and toric slack ideals*, arXiv:`1808.01692`, DOI
  `10.1016/j.jpaa.2019.106229`, Section 4.

The result also sharpens the relation to frustum cosmology.  Homogeneous
Collins--Williams/Tsuda--Fujiwara blocks impose parallel regular cells and
equal struts, thereby fixing a shift-free sector by ansatz:

- R. Tsuda and T. Fujiwara, arXiv:`2011.04120`, DOI
  `10.1093/ptep/ptab079`.

Dittrich--Gielen--Schander add six diagonals to define a simplicial Regge
action and set them to flat-frustum values in their generalized action.  The
present rank census explains locally why six is the isostatic number, but
does not make their diagonal choice canonical here:

- B. Dittrich, S. Gielen and S. Schander, arXiv:`2109.00875v3`, Section V.A.

General polytopal Regge data require flat embedding and shape matching, in
agreement with the fact that edge lengths alone become singular on the
equal-scale stratum:

- P. Donà, M. Fanizza, G. Sarno and S. Speziale, arXiv:`1708.01727`,
  Section 5.

No novelty is claimed for rigidity matrices, projective uniqueness,
hyperfrusta or ADM shift.  The exact application to this carrier and its
physical novelty remain **OPEN** pending external review.

## 8. Status and next step

- **DERIVED EXACT, ADVERSARIALLY CORROBORATED:** planar-face constraints make
  the tested unequal-scale/generic prisms locally infinitesimally determined.
- **DERIVED EXACT NEGATIVE:** equal-scale natural edge lengths admit a finite
  three-parameter non-isometric family.
- **STRUCTURAL/OPEN:** identify those three parameters with discrete shift
  and set them to zero as gauge.
- **DERIVED NEGATIVE:** rigidity supplies no selection among the 24
  staircase schedules.
- **OPEN:** glue shift/normal variables across the complete carrier, compute
  the constraint rank, and determine whether their modes are gauge,
  constrained or physical.
- **NOT YET JUSTIFIED:** a canonical local Regge Hessian at the static point.

So the programme is not blocked, but the next variable was misidentified.
Before asking for a local lapse Hessian, the theory must account for shift.

Only the two mission verifiers and static guards are run.  The full suite is
not run by explicit user instruction.

