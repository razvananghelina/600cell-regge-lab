# S01 Trivial-Pointed Counterexample

Scop: sa testam propunerea:

- `one generator + self-dual + trivial pointed subcategory`

ar putea forta Fibonacci.

Verdict:

- nu este suficient.

## Candidate Principle Under Test

Ideea era:

- contraexemplele Ising si `Rep(S_3)` esueaza deoarece apar obiecte extra
  invertibile;
- poate daca cerem ca subcategoria pointed sa fie triviala,
  atunci `X \otimes X` nu mai poate produce obiecte noi.

Aceasta ipoteza este naturala, dar trebuie testata.

## Counterexample: \texorpdfstring{$\mathrm{Rep}(A_5)$}{Rep(A5)}

Categoria `\mathrm{Rep}(A_5)` are cinci obiecte simple:

\[
\mathbf{1},\ 3,\ 3',\ 4,\ 5.
\]

Proprietati relevante:

1. subcategoria pointed este triviala:
   - singurul obiect invertibil este `\mathbf{1}`;
2. obiectul `3` este simplu si auto-dual;
3. `3` genereaza categoria.

Tensor products relevante:

\[
3 \otimes 3 \cong \mathbf{1} \oplus 3 \oplus 5,
\]

\[
3 \otimes 5 \cong 3 \oplus 3' \oplus 4 \oplus 5.
\]

Prin urmare:

- `3` este self-dual;
- pointed part este triviala;
- categoria este generata de un singur simplu netrivial `3`;
- dar
  \[
  3 \otimes 3
  \]
  introduce un nou simplu netrivial `5`.

Deci nu avem forma Fibonacci.

## Consequence

Propunerea:

- `one generator + self-dual + trivial pointed subcategory`

nu forteaza rang 2 si nu forteaza Fibonacci.

## Exact No-Go Statement

> **No-go.**
> Triviality of the pointed subcategory is not enough to exclude new simple
> summands in `X \otimes X`.

Counterexample:

- `\mathrm{Rep}(A_5)` with generator `X = 3`.

## What This Means

Directia

- `no global symmetries` / `trivial pointed part`

ramane o constrangere buna, dar este prea slaba ca principiu unic.

Ea poate fi pastrata eventual ca axiomă auxiliara, dar nu inchide `S01`.

## Verification Note

Decompozitiile de mai sus au fost verificate direct din tabelul de caractere al
lui `A_5`, folosind produsul scalar pe clase:

- `3 \otimes 3 = 1 + 3 + 5`
- `3 \otimes 5 = 3 + 3' + 4 + 5`

Asadar `3` genereaza intr-adevar intreaga categorie.
