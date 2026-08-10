# Exact Hopf symmetry-selector verdict

Date: 2026-08-10

## Provenance

The algorithm, scope and falsifiers were frozen before implementation in
commit `35b39c7` (`Preregister exact Hopf symmetry-selector test -- NO
output`).  The calculation is implemented by
`reproducible/verify_hopf_symmetry_selector.py`; its machine-readable exact
output is `reproducible/hopf_symmetry_selector.json`.

The load-bearing calculation uses rational-pair arithmetic in
`Q(sqrt(5))`.  No numerical tolerance enters any group, orbit, projector,
rank or isotropy conclusion.

## Exact group census

The independently rebuilt quaternion coordinates close exactly as a group of
120 unit elements.  The element-order distribution is

```text
order          1   2   3   4   5   6  10
multiplicity   1   1  20  30  24  20  24
```

The 24 order-ten elements generate exactly six distinct cyclic subgroups
`C10`; each subgroup has four order-ten generators.  All four generators of a
fixed subgroup determine the same unoriented imaginary line and hence the
same exact rank-one projector.

## Framing correction: six axes, twelve fibrations

For each of the six subgroups `H`, both coset constructions give a partition
of the 120 vertices into twelve ten-point fibres:

- six distinct partitions of the form `qH`;
- six distinct partitions of the form `Hq`;
- zero overlap between the two partition families.

Thus the old phrase "all six Hopf fibrations" was incomplete.  There are six
in one handed coset class and twelve after including its mirror.  Exact
quaternion conjugation bijects the two six-element families.  This is a
**DERIVED correction** to the scope of the earlier numerical enumeration; it
does not invalidate spectra computed for its explicitly chosen handed class.

## The conjugation orbit and coefficient freedom

Conjugation by `2I` induces 60 distinct permutations of the six `C10`
subgroups.  The action is transitive.  The integer constraint matrix for an
invariant coefficient vector `(c_1,...,c_6)` has rank five, hence its fixed
space is exactly one-dimensional:

`c_1=c_2=...=c_6`.

Consequently an unbroken symmetry-invariant linear combination of the six
vertical projectors is forced to use equal weights.  There is no second
invariant coefficient that could retain a preferred quadratic direction.

## Exact equiangular tight frame

For the six rank-one projectors `P_i`, the exact Gram data are

```text
Tr(P_i P_i) = 1,
Tr(P_i P_j) = 1/5  for every i != j,
sum_i P_i   = 2 I_3.
```

At every `q in S^3`, left or right unit-quaternion multiplication carries this
identity to the tangent space.  Both handed fields therefore obey

```text
sum_i P_i(q)       = 2 P_T(q),
(1/6) sum_i P_i(q) = (1/3) P_T(q).
```

The verifier checks the lifted identities exactly at all 120 group vertices;
the everywhere-continuum extension follows algebraically because unit
quaternion multiplication is orthogonal and maps the imaginary subspace onto
`T_q S^3`.

For a one-fibration anisotropic kinetic tensor

`K_i = P_T + (r-1) P_i`,

the unique symmetry-invariant average is therefore

`K_average = ((r+2)/3) P_T`.

Only an overall scalar remains.  The anisotropy itself is erased for every
value of `r`.

## What the exact five does and does not mean

The calculation gives a new exact occurrence of five:

`1 / Tr(P_i P_j) = 5` for distinct fibration axes.

Equivalently, each axis has five equiangular alternatives.  This is a
**DERIVED icosahedral frame invariant**.  Identifying it with the repository's
bootstrap symbol `a1=5` is only **STRUCTURAL** unless a separate physical map
is derived.  In particular it does not imply a propagation speed, a mass, a
coupling, Lorentzian signature or the coefficient `r=5`.

## Verdict and next opening

- **DERIVED POSITIVE:** the complete order-ten-subgroup family is finite,
  exact and rigid; the two chiral classes and their mirror relation are now
  classified.
- **DERIVED NEGATIVE:** unbroken full 600-cell symmetry selects neither one
  Hopf fibration nor a quadratic anisotropic kinetic tensor in their linear
  span.
- **DERIVED:** symmetry restoration produces the round metric, not a fourth
  direction and not the old factor five.
- **OPEN:** a chosen fibration requires a symmetry-breaking order parameter,
  boundary condition or dynamics.

This negative result points to a sharp next question.  The six axes form an
icosahedral tight frame, so their quadratic moment is isotropic.  One must
determine the first higher-order invariant moment that is not constant on the
unit two-sphere.  If it has isolated minima or maxima on the six axes, it can
support six degenerate symmetry-broken vacua.  Geometry would then provide
the allowed vacua but not, without a derived sign and coefficient, the
dynamics that chooses one.
