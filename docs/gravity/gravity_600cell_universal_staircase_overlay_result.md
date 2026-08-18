# Universal staircase overlay: exact carrier result

Date: 2026-08-17

Only the targeted exact combinatorial verifier was run.  No gravity action and
no full suite were run.

## 1. Provenance

- prior-art gate: `4fb243b`
- preregistered protocol: `4dc2feb`
- verifier registered before enumeration: `b824cab`
- verifier:
  `reproducible/verify_gravity_600cell_universal_staircase_overlay.py`
- result artifact:
  `reproducible/gravity_600cell_universal_staircase_overlay.json`
- result SHA-256:
  `0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc`

The verifier passed `12/12` preregistered controls with exact rational Z3
arithmetic.  It performed 1,116 satisfiability checks in a pruned exhaustive
traversal of all `2^14=16,384` possible sign words, then rebuilt and rechecked
every feasible full word in a fresh solver instance.

## 2. Exact local result

Inside `Delta^3 x I`, the 14 walls

```text
t = sum_(i in A) lambda_i,
empty != A != {0,1,2,3},
```

cut the strict prism interior into exactly

```text
148 full-dimensional convex chambers.
```

For every chamber and every one of the 24 vertex-order staircase
triangulations, exactly one of the four staircase four-simplices contains the
chamber.  There were zero assignment failures among `148*24=3,552` exact
chamber/order tests.  For every order, the four simplex counts are

```text
(19, 55, 55, 19).
```

Consequently, the closed cell arrangement is the common polyhedral refinement
of all 24 staircases.  Its barycentric subdivision is a schedule-free common
simplicial refinement.  This is **DERIVED COMBINATORIAL**.

## 3. Symmetry and gluing

The chamber set is invariant under all 48 elements of `S4 x C2`, with `C2`
the time reflection `t -> 1-t`.  The 148 chambers form 14 orbits:

```text
orbit size  2: 2 orbits
orbit size  4: 3 orbits
orbit size  8: 3 orbits
orbit size 12: 3 orbits
orbit size 24: 3 orbits
```

Their sizes sum to 148 and each divides 48.  On every face `lambda_j=0`, the
14 labelled walls restrict to the same six internal subset walls of a
triangular prism; each internal wall occurs twice, while `t=0` and `t=1`
occur once.  Thus the tetrahedral-prism overlays agree on shared faces and
glue functorially over the 600-cell boundary.

The resulting global polyhedral cylinder has

```text
600*148 = 88,800 full-dimensional cells.
```

This count uses only disjoint tetrahedron-prism interiors.  The complete lower
face vector and the number of four-simplices after barycentric subdivision are
still **OPEN**.

## 4. Independent identification of the number 148

After the calculation, the sign patterns were recognized as a classical
object.  Define

```text
f(A)=1  iff  sum_(i in A) lambda_i > t.
```

Every chamber gives a nonconstant positive threshold Boolean function on four
variables, and every such function admits a strict positive realization after
normalization.  The classical count is 150 positive threshold functions when
the two constants are included; excluding them gives exactly 148.  This is an
independent literature cross-check of the solver census.

Therefore:

- the number 148 is **KNOWN**, not a discovery;
- the exact correspondence provides an independent audit against a software
  counting error;
- the use of this arrangement as the universal local staircase overlay is
  **STRUCTURAL**, with external novelty **OPEN**.

Primary threshold-function enumerations include S. Muroga, I. Toda and
M. Kondo, *Majority decision functions of up to six variables*, Mathematics
of Computation 16 (1962), 459--472; the classical tables are consolidated in
S. Muroga, *Threshold Logic and Its Applications* (1971).  The modern sequence
and source trail are indexed by
[OEIS A002078](https://oeis.org/A002078).

## 5. Hostile interpretation

The result removes one precise ambiguity: no local even/odd or vertex-order
schedule needs to be selected at the carrier level.  It does **not** show that
the corresponding Regge actions agree, and it does not select a unique action
on the refinement.

The cost is substantial.  One staircase has four four-simplices per prism;
the overlay has 148 polyhedral four-cells, a factor of 37 before any
barycentric triangulation.  A direct fine action may therefore be expensive,
and integrating out the new variables requires a measure and boundary map.
Those are precisely the nontrivial ingredients of a perfect-action
construction; averaging the old even and odd actions would not supply them.

Canonicity is also relative, not absolute: the overlay is canonical after
choosing the complete class of staircase triangulations as the objects that
must be refined.  Geometry and `S4 x C2` symmetry justify treating all 24
equally, but no physical principle yet proves that a universal staircase
refinement is the correct microscopic carrier.  This limitation is
**STRUCTURAL**, not a failed control.

## 6. Status ledger

| Claim | Status |
|---|---|
| Exact local chamber count is 148 | **KNOWN / independently reproduced** |
| 148 equals nonconstant positive four-variable threshold functions | **KNOWN** |
| Every chamber refines every one of the 24 staircases | **DERIVED COMBINATORIAL** |
| Full local `S4 x C2` invariance | **DERIVED COMBINATORIAL** |
| Face restrictions glue over the 600-cell | **DERIVED COMBINATORIAL** |
| Global top-cell count is 88,800 | **DERIVED COMBINATORIAL** |
| Complete refined global f-vector | **OPEN** |
| Nondegenerate Lorentzian simplicial realization | **OPEN** |
| Fine Regge--dust action and integration measure | **OPEN** |
| Coarse/perfect action independent of schedule | **OPEN** |
| A selected nonlinear physical tick | **OPEN** |

## 7. Consequence and next falsifier

We now have a mathematically valid, schedule-free common carrier.  The next
load-bearing question is no longer whether a common refinement exists, but
whether dynamics can be transferred to it without a new arbitrary choice.

Before evaluating any gravity action, the next preregistered test should
enumerate the overlay's full face poset and ask whether there is a unique
`S4 x C2`-invariant local subdivision/coarse-graining map compatible with all
24 staircase boundary inclusions.  If multiple invariant measures or maps
remain, the carrier alone has merely relocated the fitting freedom.  If a
unique one exists, it becomes the first defensible input for an improved
Regge--dust action.

