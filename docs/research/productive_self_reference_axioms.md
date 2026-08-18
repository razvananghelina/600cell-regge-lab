# Productive Self-Reference Axioms

Scop: sa formulam un set de axiome candidate pentru programul

- `self-reference -> rank-2 non-pointed fusion structure -> Fibonacci`

si sa evaluam onest care dintre ele sunt:

- legitime;
- prea tari;
- aproape echivalente cu concluzia.

## Punctul de Pornire

Din [rank2_self_reference_theorem.md](D:\infinity\ToE\science\rank2_self_reference_theorem.md)
rezulta deja:

- daca ajungem la o categorie de fuziune de rang 2, ne-pointed,
  atunci regula Fibonacci urmeaza din clasificarea lui Ostrik.

Deci problema reala este:

- ce axiome minimale forteaza exact clasa `rank 2 + non-pointed`?

## Cadrul Formal Minim

Lucram intr-o categorie monoidala rigida, semisimpla, liniara, peste un corp
algebric inchis de caracteristica 0.

Obiectivul este sa justificam aparitia unui singur obiect simplu netrivial

\[
X
\]

care joaca rolul de `self-reference object`.

## Axiome Candidate

### A0. Unitate

Exista un obiect unitate

\[
\mathbf{1}.
\]

Status:

- inevitabil;
- nu este problematica.

### A1. Unic generator netrivial

Categoria este generata monoidal de un singur obiect simplu netrivial `X`.

Interpretare:

- `X` este singurul mod primar in care structura se aplica ei insesi.

Ce castigi:

- programul devine one-generator.

Risc:

- este deja o axiomă puternica;
- nu este echivalentă cu Fibonacci, dar este un input real.

Verdict:

- acceptabila ca axiomă candidat;
- dar trebuie justificata ulterior, nu ascunsa.

### A2. Auto-referinta interna

Prima compozitie netriviala este

\[
X \otimes X.
\]

Interpretare:

- auto-referinta este modelata prin tensorarea lui `X` cu el insusi.

Ce castigi:

- gradul 2 apare natural.

Verdict:

- buna;
- aceasta explica de ce relatia cautata este cuadratica, nu cubica.

### A3. Minimal closure of self-composition

Descompunerea lui

\[
X \otimes X
\]

nu introduce mai mult de un singur nou tip simplu netrivial.

Versiune tare:

- toate componentele netriviale din `X \otimes X` sunt izomorfe cu `X`.

Atunci rezulta imediat

\[
X \otimes X \cong a\,\mathbf{1} \oplus b\,X.
\]

Observatie:

- aceasta este exact axiomă de tip `rank-2 closure`.

Ce castigi:

- foarte mult;
- practic reduci problema la un fusion ring de rang 2.

Risc:

- este aproape toata dificultatea.

Verdict:

- buna ca teoremă-țintă;
- prea tare ca axiomă finală daca nu este motivată independent.

### A4. Productivitate

Auto-referinta nu este invertibila / triviala:

\[
X \otimes X \not\cong \mathbf{1}.
\]

Interpretare:

- self-reference trebuie sa produca structura noua, nu doar sa se anuleze
  inapoi in unitate.

Ce castigi:

- excluzi cazul pointed.

Verdict:

- esentiala;
- probabil cea mai defensabila formulare a cuvantului `productiva`.

### A5. Simplitate interna a unitatii

Multiplicitatea lui `\mathbf{1}` in `X \otimes X` este una.

Status real:

- in cadrul rigid, aceasta nu trebuie pusa ca axiomă separata;
- ea urmeaza din simplitatea lui `X`:

\[
\dim \mathrm{Hom}(\mathbf{1}, X \otimes X^*)
=
\dim \mathrm{Hom}(X, X)
= 1.
\]

Verdict:

- lemma, nu axiomă.

### A6. Auto-similaritate netriviala

Orice componenta netriviala aparuta in `X \otimes X` este din nou `X`.

Aceasta este forma explicita a auto-similaritatii.

Rezultat:

\[
X \otimes X \cong \mathbf{1} \oplus nX.
\]

Verdict:

- buna ca formulare matematica;
- dar aproape echivalenta cu presupunerea de rang 2.

## Ce Rezulta Imediat

Daca luam:

- `A1` unic generator netrivial,
- `A2` auto-referinta interna,
- `A4` productivitate,
- `A6` auto-similaritate netriviala,

atunci obtinem:

\[
X \otimes X \cong \mathbf{1} \oplus nX
\]

cu `n \ge 1`.

Apoi clasificarea lui Ostrik da:

\[
n = 1.
\]

Deci Fibonacci urmeaza.

## Problema de Onestitate

Trebuie spus clar:

- daca pui `A6` ca axiomă finală, ai împins deja aproape toata concluzia
  în premise.

Mai precis:

- `A6` este mult mai slabă decât Fibonacci exact,
- dar este încă foarte puternică.

Prin urmare, nu aceasta trebuie sa fie forma finală a programului.

## Variante de Program

### Varianta 1: Programul scurt

Axiome:

- one-generator;
- productive;
- rank-2 closure.

Avantaj:

- rapid;
- curat;
- foloseste direct clasificarea.

Dezavantaj:

- `rank-2 closure` este încă foarte aproape de concluzie.

Verdict:

- bun pentru un theorem intermediar;
- insuficient ca fundament ultim.

### Varianta 2: Programul serios

Trebuie sa derivam `rank-2 closure` din ceva mai primitiv.

Cel mai plauzibil principiu este:

- `minimal productive self-reference creates exactly one nontrivial similarity class`.

Asta trebuie formalizat.

O formulare posibila:

> dintre toate descompunerile lui `X \otimes X` compatibile cu rigiditatea,
> unitatea si productivitatea, alegem pe cea cu numar minim de clase simple
> netriviale distincte.

Daca minimul pozitiv este 1, atunci rangul 2 urmeaza.

Problema:

- trebuie demonstrat ca acest criteriu de minimalitate este structural,
  nu estetic.

### Varianta 3: Programul prin defect fata de idempotenta

Pornim de la cazul trivial

\[
X^2 = X
\]

ca model de auto-identitate pura.

Cerem:

- sa nu fie idempotent;
- sa difere minimal de idempotenta;
- cresterea sa fie exact o unitate structurala noua.

Matematic, asta ar trebui sa devina ceva de tip:

- defect minim in conul coeficientilor pozitivi ai semiringului.

Aceasta varianta este conceptually apropiata de discutia ta initiala.

Problema:

- in forma actuala nu este inca formalizata.

Verdict:

- promițătoare conceptual;
- nedefinită suficient matematic.

## Cea Mai Buna Țintă Imediată

Cea mai buna teoremă intermediară de atacat acum este:

> **Intermediate Target.**
> Let `C` be a rigid semisimple monoidal category generated by one nontrivial
> simple object `X`. Assume:
> 1. `X` is not invertible;
> 2. the tensor square `X \otimes X` contains no nontrivial simple isoclass
>    other than `X`.
>
> Then
> \[
> X \otimes X \cong \mathbf{1} \oplus nX,
> \]
> and if `C` is a fusion category over an algebraically closed field of
> characteristic 0, then necessarily `n=1`.

Aceasta este buna pentru ca:

- separa clar partea structurala proprie de partea de clasificare;
- spune exact unde mai este problema deschisa.

## Ce Trebuie Demonstrat Cu Adevarat

Nu:

- `why n = 1?`

Ci:

- `why does productive self-reference forbid the appearance of a second distinct nontrivial simple isoclass in X \otimes X?`

Asta este miezul nou.

## Criteriu de Progres Real

Vom considera ca facem progres numai daca reusim una dintre urmatoarele:

1. Demonstram un principiu care forteaza `rank-2 closure`.
2. Demonstram un no-go theorem pentru aparitia a doua sau mai multe
   clase simple netriviale in `X \otimes X` sub ipoteze minimale.
3. Aratam ca orice astfel de principiu este inevitabil echivalent cu
   o presupunere aproape-Fibonacci, caz in care programul trebuie
   reformulat.

## Verdict Curent

Ce este deja solid:

- `rank 2 + non-pointed => Fibonacci`.

Ce este inca deschis:

- cum coboram natural la `rank 2 + non-pointed`.

Deci focusul corect este:

- nu coeficientii Fibonacci;
- ci teorema de selectie a clasei `rank-2 productive self-reference`.
