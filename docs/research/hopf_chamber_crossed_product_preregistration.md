# Preregistration: the canonical six-axis crossed product on oriented chambers

Date: 2026-08-11

## Disclosed construction and risk

The static orientation audit in commit `ef2a920` closes every
`A5`-equivariant section of the oriented Hopf-axis cover.  It also points to
a different canonical object already in the repository: each of the two
oriented chamber sheets is a free 60-point `A5` orbit.

Every icosahedral chamber is a flag

```text
vertex < edge < face.
```

Its vertex determines one of the six antipodal fivefold axes.  Therefore
each chamber sheet has an exact equivariant projection

```text
p : A5/1 -> A5/D5.
```

This supplies the standard covariant representation of the already derived
crossed product

```text
B_R = R(A5/D5) crossed_product A5
    = M6(R)+M6(R)+M12(R)+M12(R)
```

on `R^60`, with axis functions acting diagonally through `p` and `A5`
acting by chamber permutations.  Doubling over the two chamber orientations
gives a canonical representation on the existing `R^120` carrier, which
already has geometric `D`, `gamma` and reflection `J`.

The following likely obstruction was noticed before implementation and is
disclosed.  Central inversion preserves each antipodal axis and commutes with
all rotations.  Hence geometric `J` may normalize the represented
noncommutative algebra rather than place its opposite in the commutant.  If
so, order zero fails.  The calculation is not blind.

## Frozen exact geometry

Reconstruct, without importing the previous gate outcomes:

1. the 12 vertices, 30 edges, 20 faces and 120 complete flags of the
   icosahedron;
2. the exact 60-element orientation-preserving permutation group `A5`;
3. its two free chamber orbits of size 60;
4. central inversion, its exchange of the sheets and its commutation with
   `A5`;
5. the six antipodal vertex axes and the chamber-to-axis projection;
6. the 3-regular chamber adjacency `D` and orientation grading `gamma`.

All incidence and permutation statements must be exact/combinatorial after
the coordinate carrier is constructed.

## Frozen crossed-product representation

On each 60-state sheet define exactly

```text
pi(delta_x) e_c = 1_(p(c)=x) e_c,
pi(u_g) e_c     = e_(g c).
```

Verify covariance on all `6*60` generator pairs.  Build the span of all
`pi(delta_x u_g)` and compute its exact dimension.  Faithfulness requires
dimension 360, the dimension of the crossed product.

Independently compute its commutant and its real Wedderburn multiplicities.
The expected decomposition `(1,1,2,2)` for simple modules of dimensions
`(6,6,12,12)` is disclosed as a consistency candidate, not assumed.

Use the same geometrically transported representation on both chamber
sheets.  No base chamber, right-regular relabelling or fitted unitary
conjugation is allowed.

## Frozen real-triple gates

For central inversion followed by conjugation, and for every one of the 60
improper geometric reflections obtained by composing it with a rotation,
compute exactly:

- `J^2` and `J gamma`;
- `JD=DJ`;
- whether `J pi(B_R) J^-1` equals, normalizes or commutes with `pi(B_R)`;
- exhaustive order zero on an exact spanning basis;
- exhaustive first order if order zero survives;
- rank of `a -> [D,pi(a)]` and the algebra-commutant dimension;
- whether inner one-forms are nonzero.

Dimension or source naming is not evidence for a lift.  The exact covariant
action above is the lift being tested.

## Acceptance and kill boundaries

- **DERIVED CHAMBER LIFT:** the natural representation is faithful and at
  least one geometric `J` passes order zero, first order, connectedness and
  nonzero forms with the derived chamber `D`.
- **DERIVED CANONICAL-LIFT NO-GO:** the representation is faithful but every
  geometric `J` fails order zero or first order.
- **CONSTRUCTION FAILURE:** the purported 360-dimensional crossed-product
  image is not faithful; then the free-flag framing was wrong.

A basepoint-dependent identification of one sheet with a right regular
torsor is excluded.  Such an identification could enlarge the commutant, but
the free orbit has 60 equally legitimate base chambers; choosing one after
the failure would be a new symmetry-breaking datum.

No matter character, mass, coupling or Standard-Model target will be used.
Only a targeted verifier will be run; the full suite remains excluded by
user instruction.
