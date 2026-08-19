# Independence gate: projected rank-edgewise carrier audit

Date: 2026-08-19

Primary artifact commit: `cb4fc24`

Primary artifact SHA-256:

```text
b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84
```

Status: written before the independent direct-split carrier is constructed.
This is an adversarial audit of a project result, not a new external novelty
claim.

## Exact claim under attack

The primary verifier used an Edelsbrunner--Grayson color-scheme enumeration
to construct

```text
P(Esd_2(sd K_600))
```

with f-vector `(19680,134880,230400,115200)`, full declared equivariance and
healthy finite geometry.  The decisive risk is that the color-scheme
implementation, its global vertex merging or its Gram-volume calculation
silently built a different complex from the intended rank-selected
`1 -> 8` refinement.

## Mechanically different reconstruction

For an ordered chamber `(0,1,2,3)`, name its edge midpoints `m_ij`.  Construct
the eight children directly as four corner tetrahedra plus the four tetrahedra
around the rank-selected central diagonal `m_02--m_13`.  This is the standard
direct description of the `k=2` edgewise split and contains no color words,
weak compositions or abacus enumeration.

Further independent choices are fixed:

- enumerate the 600 maximal tetrahedra with NetworkX maximal cliques rather
  than the primary nested common-neighbour loops;
- assemble `sd K` independently from face-containment flags;
- merge old vertices and global edge midpoints through explicit unordered
  pairs;
- compute tetrahedral volumes from the Cayley--Menger determinant rather than
  a Gram determinant;
- compare the resulting counts and geometry only with the frozen primary
  JSON, never with mutable in-memory primary objects.

## Controls

1. On one abstract tetrahedron the direct rank-selected construction must
   give 192 distinct children after all 24 barycentric flags and must be
   invariant under all `S4` permutations.
2. Applying one fixed direct central diagonal to the unranked parent must fail
   `A4` invariance.  This is the deliberately non-equivariant red control.
3. A regular Euclidean tetrahedron must have identical Gram and
   Cayley--Menger volumes at machine precision, while the actual audit uses
   only Cayley--Menger for its decisive geometry.
4. The frozen artifact hash and passing outcome must match exactly.

## Scope and framing

Agreement corroborates the finite carrier, its topology and reported
geometry.  It does not independently reconstruct the full quaternionic `H4`
action used by the primary verifier.  The exact local `S4` invariance supplies
an independent canonicity mechanism, while the primary full action remains
the global realization check.

No literature search can establish novelty for this project-specific
replication.  The primary prior-art gate already records the relevant
edgewise-subdivision and geodesic-dome literature.  No Lorentzian action,
Friedmann target or particle quantity enters this audit.
