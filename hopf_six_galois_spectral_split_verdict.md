# Fibre incidence canonically separates the golden `M12` pair

Date: 2026-08-11

Protocol commit: `43c6dd3`.

Registered verifier:
`reproducible/verify_hopf_six_galois_spectral_split.py`.
Targeted exact result: `20/20`.

No Hessian, Standard-Model module, mass, coupling or later target was used.

## Headline

The earlier six-fibration audit established the real algebra

```text
B_R=M6(R)+M6(R)+M12(R)+M12(R)
```

and showed that imposing exact golden-Galois descent makes every KO6
intersection form degenerate.  It then stated that the geometry supplied no
reason to distinguish the two golden `M12` nodes.

The no-go under rational descent remains correct, but the final statement was
too strong.  The actual 600-cell fibre incidence supplies a canonical
self-adjoint central operator whose real spectral projectors distinguish the
two nodes exactly.

**DERIVED CORRECTION:** golden conjugation is not a symmetry of the
incidence-decorated fibre.  It maps fibre edges to non-edge distance-two
chords.

This reopens node selection over the ordered real numbers.  It does not yet
select a Krajewski support or a nondegenerate Poincare pairing.

## The incidence-selected central element

Fix any derived Hopf fibration and let its effective `A5` stabilizer be the
already certified `D5`.  Each of its twelve fibres is a decagon in the
600-cell graph.  If two vertices are adjacent along such a fibre, their
relative binary element induces one of an inverse pair of order-five
rotations in `D5`:

```text
{r,r^-1}.
```

Reversing the edge exchanges the two elements, so the unoriented incidence
relation defines

```text
u_edge=r+r^-1 in R[D5].
```

There is no generator orientation or Schur coefficient in this definition.
The verifier checks every internal edge:

```text
six qH fibrations:  6*120 = 720 undirected fibre edges,
six Hq fibrations:  6*120 = 720 undirected fibre edges.
```

For all 1,440 edges, the relative actions give exactly one inverse pair.  The
other pair `{r^2,r^-2}` occurs exactly on graph-distance-two chords and never
on edges.  The selected pairs also transport correctly under every one of
the 60 `A5` rotations and all six possible base labels.  Hence the result is
independent of fibre, vertex, edge orientation, handed family and enumeration
order.  **DERIVED.**

## Exact spectrum and projectors

In the ten-dimensional left regular representation of `D5`, exact integer
arithmetic gives

```text
charpoly(u_edge)
  =(lambda-2)^2 (lambda^2+lambda-1)^4.
```

Writing `phi=(1+sqrt(5))/2`, the spectrum is

```text
2       multiplicity 2,
phi-1   multiplicity 4,
-phi    multiplicity 4.
```

The two four-dimensional regular isotypic spaces are precisely the two real
`M2(R)` Wedderburn blocks of `R[D5]`.  After the canonical six-point Morita
amplification they become the two `M12(R)` blocks of `B_R`.

The spectral projectors `P_+` and `P_-` are exact symmetric central
idempotents of ranks four.  Golden conjugation exchanges them.  Their sum is
rational, while neither projector separately descends to `Q[D5]`.

The sign is not a fitted threshold: `phi-1>0` and `-phi<0`, so zero separates
the two spectra.  More invariantly, the complete spectral decomposition
labels both nodes even if no one node is privileged.  Replacing adjacency by
its negative would exchange the words positive and negative but would not
erase the separation.

## Why the old Galois constraint is not geometric

The abstract `D5` automorphism

```text
tau: r -> r^2
```

implements the golden character conjugation.  The verifier proves exactly

```text
tau(u_edge)=r^2+r^-2=u_chord,
tau(P_-)=P_+.
```

But `u_chord` is supported on distance-two non-edges, not fibre edges.
Therefore `tau` is an automorphism of the undecorated stabilizer and its
character table, but not of the stabilizer together with the derived
600-cell incidence.

The previous Pfaffian calculation remains a valid theorem with its full
hypothesis:

> If the four-node real intersection form is required to descend to the
> unsplit rational algebra, equivalently to preserve or reverse grading under
> `tau`, then its Pfaffian is zero.

What is no longer justified is silently elevating that arithmetic descent to
a symmetry axiom after the real incidence operator has distinguished edges
from chords.  Real finite spectral triples are built over `R`; ordered-real
functional calculus is legitimate.  Whether the theory should nevertheless
impose rational descent is now an extra hypothesis, not a geometric
necessity.

## All four nodes are geometrically labelled

The unique five-element reflection class of `D5` supplies a second canonical
central operator

```text
v_ref=sum_(s reflection) s.
```

The joint exact spectrum is

```text
real Wedderburn node       (u_edge,v_ref)
trivial                    (2, 5)
reflection sign            (2,-5)
positive golden doublet    (phi-1,0)
negative golden doublet    (-phi,0).
```

Thus the decorated geometry labels all four simple real summands without an
arbitrary same-dimension convention.  **DERIVED NODE SEPARATION.**

## Hostile framing audit

This is useful, but it is not yet the finite Standard-Model geometry.

1. A central separator labels nodes; it does not choose oriented arrows or
   their multiplicities in a Krajewski diagram.
2. Consequently it does not choose any of the six independent entries of a
   generic four-node antisymmetric KO6 intersection form.
3. The old canonical faithful doubles still fail order zero or
   orientability.  This calculation does not repair those carriers.
4. Selecting the negative block alone is unnecessary for the theorem and
   would require explaining the physical use of a spectral cut.  The robust
   result is the canonical joint labelling of all blocks.
5. The operator is derived from the already selected Hopf-fibration sector.
   It does not solve the separate OPEN problem of why one handed fibration
   sector is dynamically realized; the split itself was nevertheless checked
   in both hands.
6. The protocol was preregistered after the candidate was noticed, so the
   result is disclosed rather than blind.  Its evidential force comes from
   exhaustive incidence and exact algebra, not from a target hit.

## Status ledger

- **DERIVED:** all 1,440 internal edges across both handed six-fibration
  families select one inverse rotation class; all distance-two chords select
  the other.
- **DERIVED:** the selection is equivariant under all 60 `A5` rotations.
- **DERIVED:** `u_edge` has exact spectrum
  `2, phi-1, -phi` with multiplicities `2,4,4`.
- **DERIVED:** its two golden spectral projectors are the two real `M2`
  blocks and hence the two crossed-product `M12` blocks.
- **DERIVED:** `tau:r->r^2` realizes golden conjugation but maps edges to
  chords, so it is not an incidence symmetry.
- **DERIVED:** `(u_edge,v_ref)` labels all four real Wedderburn nodes.
- **DERIVED, scoped negative:** rational/Galois-compatible KO6 Poincare
  duality remains impossible.
- **STRUCTURAL OPENING:** the split real four-node arena no longer lacks a
  canonical node labelling.
- **OPEN:** a geometry-selected Krajewski support, grading, real structure,
  nondegenerate intersection form, first-order Dirac operator and nonzero
  inner one-forms.
- **NO TARGET COMPARISON:** no Hessian or particle module was inspected.

## Next admissible gate

The next calculation must use the now-labelled central projectors to build
all bimodule correspondences directly induced by fibre edge and chord
incidence, before looking at any desired intersection matrix.  Their
antisymmetrized multiplicity matrices must then be computed exactly.

The route advances only if that incidence construction itself selects a
nondegenerate KO6 form.  Merely choosing one of the many possible four-node
skew matrices after seeing that its Pfaffian is nonzero would reintroduce the
same fitting problem this audit was designed to remove.

## Subsequent carrier result

The preregistered continuation is recorded in
`hopf_six_spectral_krajewski_verdict.md`.  The full off-diagonal enveloping
bimodule passes order zero, KO6 signs, metric-zero orientability and
unimodular Poincare duality for all eight legitimate signed lexicographic
readings of the joint node spectrum.  This is robust structural existence,
not yet a selected Dirac operator or finite spectral triple.
