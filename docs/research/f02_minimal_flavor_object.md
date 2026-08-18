# F02 Closure: Minimal Flavor Object

## Exact claim

Define the minimal flavor object that can be built from the admissible
post-gauge-collapse input set without importing gauge structure.

## Formal setting

Arithmetic scaffold over `Z[phi]` with:

- a three-slot unit structure;
- an integer suppression lattice;
- an optional conditional exponent assignment.

## Inputs used

- [f01_allowed_flavor_inputs.md](D:\infinity\ToE\science\f01_allowed_flavor_inputs.md)
- [one_integer_paper_exact_core.tex](D:\infinity\ToE\science\one_integer_paper_exact_core.tex)
- [flavor_first_program.md](D:\infinity\ToE\science\flavor_first_program.md)

## Output status

- `Conditional flavor statement`

## The minimal object

The minimal admissible flavor object is the quadruple
\[
\mathcal{F}_{\mathrm{min}} = (S,L,w,N)
\]
with:

1. family-slot set
   \[
   S=\{0,1,2\},
   \]
   coming from the exact theorem of three stable unit sectors on the chiral
   line;

2. charge lattice
   \[
   L=\mathbb{Z}^2,
   \]
   whose points are written `(a,b)`;

3. weight map
   \[
   w:L\to\mathbb{Z},\qquad w(a,b)=5a+6b;
   \]

4. arithmetic norm
   \[
   N(a+b\phi)=a^2+ab-b^2.
   \]

This object is exact up to here.

## Conditional flavor reading

The minimal conditional flavor reading is:

- `S` gives three candidate family slots;
- `w(a,b)` gives a candidate suppression exponent;
- `N(a+b\phi)` gives a candidate secondary arithmetic discriminator;
- if one accepts the conditional exponent set
  \[
  \{0,3,5,11,11,16,17,19,26\},
  \]
  then the exact `(a,b)` theorem supplies a unique lattice placement of the
  nine listed exponents.

## Why this is nontrivial

This object is not empty formalism. It already contains:

- a rigid three-slot structure;
- a nontrivial integer grading;
- an exact arithmetic norm;
- a conditional uniqueness theorem for placing the known exponent set.

What it does not yet contain is any claim about Yukawa operators, gauge
quantum numbers, or electroweak representations.

## Failure criterion

`F02` would fail only if every candidate flavor object required forbidden
gauge-derived input. That does not happen here: the object above uses only the
admissible arithmetic and McKay-based ingredients.

## Decision

`F02` is closed as `Conditional flavor statement`.
