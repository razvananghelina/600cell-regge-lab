# F03 Closure: Exponent Lattice as Suppression-Charge Scaffold

## Exact claim

Decide whether the exponent lattice
\[
n=5a+6b
\]
can be read as a genuine suppression-charge structure without importing gauge
or electroweak data.

## Formal setting

Arithmetic lattice `L = Z^2` with weight map
\[
w(a,b)=5a+6b
\]
and norm
\[
N(a+b\phi)=a^2+ab-b^2.
\]

## Inputs used

- [f02_minimal_flavor_object.md](D:\infinity\ToE\science\f02_minimal_flavor_object.md)
- [f01_allowed_flavor_inputs.md](D:\infinity\ToE\science\f01_allowed_flavor_inputs.md)
- [one_integer_paper_exact_core.tex](D:\infinity\ToE\science\one_integer_paper_exact_core.tex)

## Output status

- `Conditional flavor statement`

## Exact arithmetic content

The map
\[
w:L\to Z,\qquad w(a,b)=5a+6b
\]
is a surjective group homomorphism because `gcd(5,6)=1`.

Its kernel is the rank-1 lattice
\[
\ker w = Z\cdot (6,-5),
\]
since
\[
5a+6b=0 \iff (a,b)=t(6,-5),\quad t\in Z.
\]

Therefore `w` defines an exact one-dimensional integer grading on the
2-dimensional lattice `L`.

## Conditional suppression-charge reading

This exact grading can be read as a suppression-charge scaffold in the
following precise and limited sense:

1. differences in exponent are additive:
   \[
   \Delta n = w(\Delta a,\Delta b);
   \]
2. if one introduces a phenomenological small parameter `epsilon`, then a
   candidate hierarchy ansatz has the generic form
   \[
   \text{weight} \sim \epsilon^{\,w(a,b)};
   \]
3. if one accepts the conditional exponent set and the exact `(a,b)` theorem,
   the corresponding suppression charges are uniquely placed on the lattice.

So the exact lattice data support a flavor-style suppression language.

## What this does not prove

The map `w` alone does not provide:

1. Yukawa operators;
2. gauge quantum numbers;
3. electroweak representations;
4. texture zeros or full mixing matrices;
5. separation of states with equal exponent, such as the repeated value `11`.

That last point is important: `w` is only a one-dimensional grading. Distinct
states with the same exponent require additional data, such as the full
`(a,b)` placement, the family-slot structure, or further conditional flavor
dictionary.

## Decision

`F03` is closed as `Conditional flavor statement`:

- exact theorem-level content:
  `w(a,b)=5a+6b` is a surjective integer grading with kernel `Z(6,-5)`;
- conditional flavor content:
  this grading may be used as a suppression-charge scaffold;
- forbidden overread:
  it is not yet a full FN/Yukawa model.
