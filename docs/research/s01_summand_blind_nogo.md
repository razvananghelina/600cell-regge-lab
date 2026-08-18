# S01 Summand-Blind No-Go

Scop: sa vedem daca putem formula un `no-go` mai general decat simplele
contraexemple izolate.

Ideea este urmatoarea:

- poate exista un principiu mai slab decat `no-branching`;
- dar daca acel principiu nu vede explicit ce simple netriviale distincte apar
  in `X \otimes X`, atunci probabil nu poate exclude ramificarea.

Acest document formalizeaza exact aceasta idee.

## Informal Definition

Vom numi un principiu `summand-blind at tensor-square level` daca el foloseste
doar proprietati de tip:

- `X` genereaza categoria;
- `X` este self-dual;
- `X` este non-invertibil;
- `1` apare in `X \otimes X`;
- subcategoria pointed este triviala;
- alte proprietati globale care nu disting ce simplu netrivial nou apare in
  descompunerea lui `X \otimes X`.

Adica:

- principiul nu inspecteaza identitatea summand-urilor netriviale distincte de
  `X` din `X \otimes X`.

## Why This Class Matters

Toate candidatele naturale incercate pana acum au fost de acest tip:

1. `productive = non-invertible`;
2. `self-return`;
3. `trivial pointed part`;
4. `one-generator`;
5. `self-duality`.

Dar contraexemplele au aratat ca niciun astfel de pachet nu ajunge.

## Core Counterexample Pattern

Categoria Fibonacci:

\[
\tau \otimes \tau \cong 1 \oplus \tau.
\]

Categoria `\mathrm{Rep}(A_5)` cu generator `X = 3`:

\[
3 \otimes 3 \cong 1 \oplus 3 \oplus 5.
\]

Pentru `X = 3` in `\mathrm{Rep}(A_5)` avem simultan:

1. `X` este simplu;
2. `X` este self-dual;
3. `X` este non-invertibil;
4. `1` apare in `X \otimes X`;
5. subcategoria pointed este triviala;
6. `X` genereaza categoria.

Singura diferenta relevanta fata de Fibonacci la nivelul tensorului patrat este:

- apare un simplu netrivial nou `5`.

## Theorem-Style No-Go

> **Scoped no-go theorem.**
> No tensor-square principle that is blind to the presence of extra nontrivial
> simple summands in `X \otimes X` can force Fibonacci from the data
> `one-generator + self-dual + non-invertible + unit return + trivial pointed
> subcategory`.

## Proof Idea

Fibonacci satisfies all five displayed properties.

So does `\mathrm{Rep}(A_5)` with `X = 3`.

But Fibonacci has

\[
X \otimes X \cong 1 \oplus X,
\]

while `\mathrm{Rep}(A_5)` has

\[
X \otimes X \cong 1 \oplus X \oplus Y
\]

with `Y = 5 \not\cong X`.

Therefore any principle that does not distinguish these two situations at the
level of extra nontrivial summands cannot force the Fibonacci rule.

## Consequence

This is stronger than the earlier isolated counterexamples.

It says:

- the missing ingredient is not merely
  - non-invertibility,
  - self-duality,
  - trivial pointed part,
  - generation,
  - or unit return;
- the missing ingredient must be explicitly `summand-sensitive`.

In plain language:

- any successful foundational axiom must see branching.

## Structural Corollary

Any realistic candidate principle for `S01` must do one of the following:

1. directly forbid new simple isoclasses in `X \otimes X`;
2. imply a theorem that forbids them;
3. be strong enough to be essentially equivalent to `no-branching`.

If it does none of these, it cannot close `S01`.

## Strategic Interpretation

Acesta este aproape verdictul final pentru `S01`.

Nu demonstreaza inca in mod absolut ca `no-branching` este axioma ireductibila.

Dar demonstreaza ceva foarte apropiat:

- orice principiu viabil trebuie sa fie `summand-sensitive`;
- adica trebuie sa controleze explicit aparitia lui `Y \neq X` in `X \otimes X`.

Prin urmare, distanta dintre un astfel de principiu si axioma `no-branching`
este acum foarte mica.

## Current Best Reading

Cea mai onesta concluzie de lucru este:

- fie gasim un principiu `summand-sensitive` cu adevarat independent;
- fie acceptam `no-branching` ca axioma ireductibila.

In acest moment, a doua varianta pare mai probabila.
