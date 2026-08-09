# Rank-2 Self-Reference Theorem Program

Scop: sa reformulam ideea

- `auto-referinta productiva forteaza Fibonacci`

intr-o forma matematica exacta.

## Executive Result

Exista deja un nucleu theorem foarte puternic:

in limbaj de categorie de fuziune / fusion ring de rang 2,

- forma cuadratica nu este arbitrara;
- coeficientul lui `1` nu este arbitrar;
- singurul caz ne-pointed categorificabil este exact Fibonacci.

Mai precis, daca `C` este o categorie de fuziune de rang 2 cu obiecte simple

\[
\mathbf{1}, X,
\]

atunci:

1. produsul lui `X` cu el insusi are in mod necesar forma
   \[
   X \otimes X \cong \mathbf{1} \oplus n X
   \]
   sau, in cazul pointed,
   \[
   X \otimes X \cong \mathbf{1};
   \]
2. daca categoria este ne-pointed, atunci Ostrik arata ca singura posibilitate
   categorificabila este
   \[
   X \otimes X \cong \mathbf{1} \oplus X.
   \]

Deci:

- `Fibonacci` nu trebuie vazut ca un ansatz liber in clasa rank-2 non-pointed;
- el este deja teorema de clasificare in acea clasa.

## Sursa Primara

- Viktor Ostrik, *Fusion categories of rank 2*, arXiv:math/0203255
  https://arxiv.org/abs/math/0203255

Punctele relevante din sursa:

- rangul 2 inseamna doua clase de simple `1, X`;
- regulile de fuziune sunt determinate de un numar `n` din
  \[
  X^2 = \mathbf{1} + nX;
  \]
- theorem principal:
  exista doar 4 categorii de rang 2:
  doua pentru `n=0` si doua pentru `n=1`;
  pentru `n \ge 2` nu exista categorificare.

## Ce Este Deja Demonstrabil

### Lemma 1: De ce relatia este de grad 2

Intr-o categorie de fuziune de rang 2, baza Grothendieck este

\[
\{\mathbf{1}, X\}.
\]

Prin urmare,

\[
X \otimes X
\]

nu are unde sa se inchida decat in span-ul acestei baze:

\[
X \otimes X \cong a\,\mathbf{1} \oplus b\,X
\]

cu coeficienti intregi nenegativi.

Concluzie:

- `gradul 2` nu vine din estetica;
- vine din faptul ca prima auto-compozitie netriviala a unui unic generator
  intr-o structura de rang 2 trebuie sa se inchida quadratic.

### Lemma 2: De ce coeficientul lui `1` este exact 1

Intr-o categorie rigida, pentru orice simplu `X`,

\[
\dim \operatorname{Hom}(\mathbf{1}, X \otimes X^*)
=
\dim \operatorname{Hom}(X, X)
= 1.
\]

In rang 2, daca nu suntem in cazul pointed, atunci `X^* = X`.

Deci multiplicity-ul unitatii in `X \otimes X` este exact 1:

\[
X \otimes X \cong \mathbf{1} \oplus n X.
\]

Concluzie:

- termenul `+1` nu este un artificiu;
- el este fortat de rigiditate + simplitatea lui `X`.

### Lemma 3: De ce cazul pointed este separat

Daca `X^* = \mathbf{1}`, atunci `X` este invertibil si

\[
X \otimes X \cong \mathbf{1}.
\]

Acesta este cazul trivial / pointed.

Deci exista o bifurcatie structurala reala:

1. pointed:
   \[
   X^2 = 1;
   \]
2. ne-pointed:
   \[
   X^2 = 1 + nX.
   \]

### Theorem 1: Singurul caz ne-pointed categorificabil este Fibonacci

Prin clasificarea lui Ostrik:

- pentru fusion rings de rang 2, cazurile categorificabile sunt doar `n=0`
  si `n=1`;
- `n=0` este pointed;
- `n=1` este cazul Fibonacci / Yang-Lee.

Prin urmare:

> **Theorem.**
> In a rank-2 fusion category, if the category is non-pointed, then the fusion
> rule is necessarily
> \[
> X \otimes X \cong \mathbf{1} \oplus X.
> \]

Acesta este exact seed-ul Fibonacci.

## Ce Inseamna Pentru Programul Tau

Programul tau poate fi reformulat astfel:

Nu trebuie sa demonstram direct

- `existenta -> x^2 = x + 1`

ci ceva mai precis si mai realist:

- `un principiu structural minimal forteaza o categorie de fuziune de rang 2 sa fie ne-pointed`.

Odata obtinut acest pas, Fibonacci vine gratis din clasificare.

Aceasta este o schimbare strategica majora.

## Noul Lant Corect

In loc de:

- `existenta -> auto-referinta -> x^2 = x + 1`

putem tinti:

- `principiu de auto-referinta productiva`
  ->
  `structura monoidala rigida cu un unic generator netrivial`
  ->
  `rang 2`
  ->
  `ne-pointed`
  ->
  `Fibonacci prin clasificare`
  ->
  `phi`
  ->
  restul lantului.

Acest lant este mult mai serios matematic.

## Unde Este Acum Adevarata Dificultate

Nu in:

- `rank 2 non-pointed => Fibonacci`

aceasta parte este practic rezolvata.

Dificultatea reala este sa justificam dintr-un principiu mai primitiv:

1. de ce structura trebuie sa aiba un singur generator netrivial;
2. de ce inchiderea trebuie sa aiba rang 2;
3. de ce trebuie exclus cazul pointed.

Acestea sunt acum cele trei probleme fundamentale.

## Cea Mai Buna Formulare a Problemei

### Problem A

Gaseste un principiu structural care forteaza:

\[
\mathrm{Irr}(C)=\{\mathbf{1}, X\}.
\]

Adica:

- exact un obiect simplu netrivial.

### Problem B

Gaseste un principiu structural care exclude:

\[
X \otimes X = \mathbf{1}.
\]

Adica:

- auto-referinta nu este trivial-invertibila;
- ea trebuie sa fie productiva.

### Problem C

Arata ca principiul folosit este mai slab decat a presupune direct regula
Fibonacci.

Altfel nu am castigat nimic.

## Cea Mai Buna Teorema-Tinta

O versiune buna a theorem-ului tinta ar fi:

> **Target Theorem.**
> Let `C` be a rigid semisimple monoidal category generated by one nontrivial
> simple object `X`. Assume:
> 1. the closure of tensor powers of `X` has rank 2;
> 2. `X` is not invertible;
> 3. `C` is a fusion category over an algebraically closed field of
>    characteristic 0.
>
> Then
> \[
> X \otimes X \cong \mathbf{1} \oplus X.
> \]

Aceasta theorem nu este noua ca rezultat final, dar este exact forma corecta
pentru programul nostru.

## Ce Ramane De Demonstrat Cu Adevarat

Ca sa impingem mai jos nivelul de axiome, trebuie sa derivam:

- `rank 2`
- `one nontrivial generator`
- `non-invertibility`

dintr-un principiu de `auto-referinta productiva`.

Asta este partea cu adevarat noua si grea.

## Concluzie Strategica

Programul de lucru se schimba astfel:

1. Nu mai incercam sa derivam direct coeficientii Fibonacci.
2. Folosim theorem-ul de clasificare existent ca bloc final.
3. Atacam doar ce este cu adevarat deschis:
   - de ce apare exact clasa `rank 2 non-pointed`.

Acesta este, probabil, cel mai eficient mod de a transforma ideea
`zero-postulate` intr-un program matematic serios.
