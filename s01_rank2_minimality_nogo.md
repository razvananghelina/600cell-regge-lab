# S01 Rank-2 Minimality No-Go

Scop: sa testam daca expresii precum

- `minimal productive self-reference`
- `prima deformare netriviala`
- `minimal closure`

pot forteza in mod real `rank 2`, fara sa ascunda deja concluzia.

## Intrebarea

Putem obtine

\[
X \otimes X \cong \mathbf{1} \oplus nX
\]

cu un argument de minimalitate care sa fie substantial mai slab decat
presupunerea directa a inchiderii de rang 2?

## Observatia Centrala

Ca sa ajungi la regula de mai sus, trebuie sa excluzi explicit posibilitatea ca

\[
X \otimes X
\]

sa contina inca un simplu netrivial distinct

\[
Y \not\cong X.
\]

Dar exact aceasta excludere este continutul matematic al trecerii la rang 2.

## Lemma 1

In orice categorie semisimpla monoidala generata de un obiect simplu `X`,
descompunerea lui

\[
X \otimes X
\]

are forma

\[
X \otimes X \cong a\,\mathbf{1} \oplus b\,X \oplus \bigoplus_i c_i Y_i,
\]

unde `Y_i` sunt simple netriviale distincte de `X`.

Concluzie:

- `rank 2` este exact afirmatia ca toti termenii `Y_i` lipsesc.

## Lemma 2

Orice principiu de forma:

- `self-reference introduces the minimal number of new nontrivial simple types`

nu produce singur o teoremă pana cand nu se specifica exact:

1. clasa peste care se ia minimul;
2. functia de cost;
3. de ce acea functie de cost este structurala si nu arbitrara.

Fara acestea, cuvantul `minimal` nu are continut theorematic.

## Lemma 3

Daca functia de cost este aleasa sa fie

- `number of nontrivial simple isoclasses appearing in X ⊗ X`,

atunci minimul pozitiv egal cu 1 este exact o reformulare a presupunerii:

- `X ⊗ X introduces no second distinct nontrivial simple isoclass`.

Dar aceasta este tocmai presupunerea care forteaza rangul 2.

Concluzie:

- criteriul de minimalitate nu mai este substantial mai slab decat concluzia.

## No-Go Statement

> **No-go statement.**
> Any derivation of rank-2 closure from a bare slogan such as
> `minimal productive self-reference`
> is empty unless `minimal` is formalized by an independent structural
> criterion. If `minimal` is formalized directly as
> `no second distinct nontrivial simple isoclass appears in X \otimes X`,
> then the argument does not derive rank 2; it assumes it in disguised form.

## Ce Este Totusi Salvat

Acest no-go nu ucide programul.

El spune doar:

- nu putem folosi cuvantul `minimal` in mod liber;
- trebuie gasit un criteriu independent, care sa forteze lipsa lui `Y`.

Programul rămâne posibil doar daca identificam un principiu de alt tip, de
exemplu:

1. un principiu universal de auto-similaritate;
2. o constrangere de rigiditate / universalitate a generatorului;
3. un no-go pentru branching al auto-compozitiei intr-o clasa bine definita.

## Verdict

Rezultat curent pentru `S01`:

- nu avem inca theorem care sa derive `rank 2`;
- dar avem un no-go util:
  - `minimal closure` in forma naiva nu este o demonstratie;
  - fie este vag, fie este doar concluzia reformulata.

Acesta este progres real, pentru ca elimina o directie inselatoare.
