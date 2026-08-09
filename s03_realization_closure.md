# S03 Realization Closure

Scop: sa inchidem `S03` in forma maximala adevarata.

Problema initiala era formulata prea tare:

- `a_1 = 5 fixes the binary icosahedral / 600-cell realization`.

In aceasta forma, afirmatia nu este demonstrata.

Trebuie inlocuita cu doua rezultate mai precise:

1. un `derived lemma` de selectie intr-o clasa bine definita;
2. un `scoped no-go` impotriva pretentiei de selectie absoluta din `a_1=5`
   singur.

## Exact Input

Din `S02` avem:

\[
a_1 = 5,
\qquad
\phi = \frac{1+\sqrt{5}}{2},
\qquad
\mathbb{Q}(\sqrt{5}).
\]

## Derived Lemma

> **Lemma S03.**
> Assume one seeks a realization of the bootstrap output inside the class of
> regular convex 4-polytopes, and require that the realization intrinsically
> carry:
> 1. the golden-ratio field `\mathbb{Q}(\sqrt{5})`;
> 2. the `H_4` / McKay-`E_8` package;
> 3. vertex set cardinality `|V| = |2I| = 120`;
> 4. local degree `12`.
>
> Then the unique regular convex 4-polytope satisfying all four requirements is
> the 600-cell.

## Proof

Among the six regular convex 4-polytopes:

1. only the `H_4` dual pair
   \[
   \{5,3,3\},\ \{3,3,5\}
   \]
   intrinsically carries the golden-ratio field `\mathbb{Q}(\sqrt{5})`;
2. this same pair is the only one lying in the `H_4` / McKay-`E_8` class;
3. among that pair, the 120-cell has
   \[
   |V|=600,\ \deg(v)=4,
   \]
   while the 600-cell has
   \[
   |V|=120,\ \deg(v)=12.
   \]

Therefore the simultaneous requirements

\[
\mathbb{Q}(\sqrt{5}),\quad H_4/E_8,\quad |V|=120,\quad \deg(v)=12
\]

select the 600-cell uniquely.

This is verified independently by

- [reproducible/verify_polytope_selection_intrinsic.py](D:\infinity\ToE\science\reproducible\verify_polytope_selection_intrinsic.py)

which passes on all six regular convex 4-polytopes.

## Scoped No-Go

> **No-go.**
> The arithmetic datum `a_1 = 5` by itself does not uniquely force the
> 600-cell realization inside the class of regular convex 4-polytopes.

## Reason

The same bootstrap field `\mathbb{Q}(\sqrt{5})` is intrinsically shared by both
members of the `H_4` dual pair:

- 120-cell;
- 600-cell.

Therefore the bootstrap arithmetic alone narrows the realization class to
`H_4`, but does not break the dual-pair ambiguity.

The ambiguity is broken only after adding the exact intrinsic selectors:

- `|V| = |2I| = 120`;
- local degree `12`.

## Closure Decision

`S03` is closed in the following honest form:

- `Derived lemma`:
  within the regular-convex 4-polytope realization class carrying the
  bootstrap field and the `H_4/E_8` package, the 600-cell is uniquely selected
  by vertex count and local degree.
- `Scoped no-go`:
  `a_1 = 5` alone does not absolutely force the 600-cell.

## What Is Not Claimed

This closure does not claim:

- that every conceivable realization of `a_1 = 5` must be the 600-cell;
- that pure arithmetic alone singles out the 600-cell with no auxiliary
  structural conditions.

It claims only the strongest exact statement presently supported.
