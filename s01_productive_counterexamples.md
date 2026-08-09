# S01 Productive Self-Reference Counterexamples

Scop: sa testam daca axiome candidate aparent rezonabile sunt suficiente pentru
Fibonacci sau macar pentru rang 2.

## Intrebarea 1

Este suficient:

1. un singur generator netrivial `X`;
2. `X` ne-invertibil;

ca sa forteze inchiderea de rang 2?

Raspuns:

- nu.

## Counterexample 1: Ising

Categoria Ising are trei obiecte simple:

\[
\mathbf{1}, \psi, \sigma
\]

cu reguli de fuziune standard

\[
\psi \otimes \psi \cong \mathbf{1}, \qquad
\psi \otimes \sigma \cong \sigma, \qquad
\sigma \otimes \sigma \cong \mathbf{1} \oplus \psi.
\]

Proprietati:

1. `\sigma` este simplu si ne-invertibil;
2. categoria este generata de `\sigma`;
3. dar
   \[
   \sigma \otimes \sigma \not\cong \mathbf{1} \oplus n \sigma.
   \]

Concluzie:

- `one generator + non-invertible` nu forteaza forma Fibonacci;
- nici macar nu forteaza ca `X` sa reapara in `X \otimes X`.

Aceasta elimina ideea:

- `productive = non-invertible` este suficient.

## Intrebarea 2

Poate este suficient sa adaugam si:

3. `X` reapare in propria auto-compozitie.

Adica:

\[
X \subset X \otimes X.
\]

Raspuns:

- tot nu.

## Counterexample 2: \texorpdfstring{$\mathrm{Rep}(S_3)$}{Rep(S3)}

In categoria de reprezentari finite-dimensionale a lui `S_3`, exista trei
simple:

\[
\mathbf{1}, \varepsilon, \rho
\]

unde:

- `\mathbf{1}` este triviala;
- `\varepsilon` este reprezentarea semn;
- `\rho` este reprezentarea standard de dimensiune 2.

Regula de tensor relevantă este:

\[
\rho \otimes \rho \cong \mathbf{1} \oplus \varepsilon \oplus \rho.
\]

Proprietati:

1. `\rho` este simplu;
2. `\rho` este ne-invertibil;
3. `\rho` genereaza categoria;
4. `\rho` reapare in `\rho \otimes \rho`.

Dar:

\[
\rho \otimes \rho
\]

contine un simplu suplimentar distinct, anume `\varepsilon`.

Concluzie:

- `one generator + non-invertible + self-return`
  nu forteaza rang 2;
- nici nu forteaza forma
  \[
  X \otimes X \cong \mathbf{1} \oplus nX.
  \]

## No-Go Statement 1

> **No-go.**
> The slogan `productive self-reference = non-invertible self-composition`
> is too weak to force Fibonacci or rank-2 closure.

Counterexample:

- Ising.

## No-Go Statement 2

> **No-go.**
> Even the stronger package
> `one generator + non-invertible + self-return`
> is still too weak to force rank 2.

Counterexample:

- `\mathrm{Rep}(S_3)` with generator `\rho`.

## Ce Ramane Posibil

Aceste contraexemple nu distrug programul.

Ele spun doar ca, pentru a ajunge la Fibonacci, trebuie o axioma suplimentara
care exclude aparitia unei noi clase simple distincte in `X \otimes X`.

Deci nu este suficient:

- nici `productive`;
- nici `productive + self-return`.

Trebuie ceva mai tare, dar inca legitim.

## Miezul Exact al Problemei

Acum stim mai precis ce trebuie sa forteze un principiu bun:

\[
X \otimes X
\]

trebuie sa nu contina niciun simplu netrivial distinct de `X`.

Adica, principiul cautat trebuie sa excluda exact termenii de tip:

\[
Y \not\cong X.
\]

Acesta este punctul matematic precis ramas deschis.

## Surse

- Ising fusion rules are standard; see e.g. the overview statement in
  *Categorical Fermionic Actions and Minimal Modular Extensions*:
  https://sigma-journal.com/2025/085/sigma25-085.pdf
- General fusion-category background:
  Etingof--Nikshych--Ostrik, *On fusion categories*,
  https://annals.math.princeton.edu/2005/162-2/p01
