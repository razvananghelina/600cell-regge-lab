# S04 Spectrum Closure

Scop: sa inchidem `S04` in forma maxima justificata.

## Exact Claim

Pentru graful 600-cell construit ca graf Cayley al lui `2I`, Laplacianul scalar

\[
\Delta_0 = 12I - A
\]

are exact noua valori proprii distincte, toate in `\mathbb{Z}[\phi]`, cu
multiplicitati

\[
1,4,9,16,25,36,9,16,4.
\]

## Status

- `Computational fact`

dar cu substructura:

- valorile proprii si multiplicitatile exacte: `Computational fact`;
- proprietatea de patrate perfecte: `Derived representation-theoretic consequence`;
- localizarea lui `3` si `3'`: `Representation-theoretic consequence` verified by character computation.

## What Is Exactly Verified

Scriptul

- [reproducible/verify_spectrum_600cell.py](D:\infinity\ToE\science\reproducible\verify_spectrum_600cell.py)

verifica direct:

1. constructia celor 120 varfuri ale lui `2I`;
2. graful 12-regulat cu 720 muchii;
3. existenta exact a 9 valori proprii distincte;
4. faptul ca toate valorile proprii sunt in `\mathbb{Z}[\phi]`;
5. multiplicitatile exacte;
6. perechile Galois;
7. localizarea ireps-urilor `3` si `3'` in eigenspatiile de multiplicitate 9.

Rularea efectiva a scriptului trece integral.

## Exact Spectrum

Valorile proprii ale lui `\Delta_0` sunt:

\[
0,
\quad
12-6\phi,
\quad
12-4\phi,
\quad
9,
\quad
12,
\quad
14,
\quad
8+4\phi,
\quad
15,
\quad
6+6\phi.
\]

cu multiplicitati:

\[
1,4,9,16,25,36,9,16,4.
\]

## Classification of the Additional Features

Nu toate proprietatile frumoase ale spectrului au acelasi statut logic.

### A. Exact finite spectral data

Acestea sunt pur si simplu rezultate de calcul finit verificat:

- cele 9 valori proprii;
- multiplicitatile lor exacte;
- faptul ca valorile stau in `\mathbb{Z}[\phi]`.

### B. Why the multiplicities are perfect squares

Aceasta nu este tratata aici ca o coincidenta misterioasa.

Pentru un graf Cayley al unui grup finit, operatorul de adiacenta este
convolutie pe reprezentarea regulata stanga. Cand generating set-ul este o
clasa de conjugare, operatorul apartine centrului algebrei de grup, deci prin
Schur actioneaza scalar pe fiecare bloc izotipic ireductibil.

In reprezentarea regulata a lui `2I`, fiecare irrep `\rho` de dimensiune
`d_\rho` apare cu multiplicitate `d_\rho`, deci blocul sau are dimensiune

\[
d_\rho^2.
\]

Prin urmare, cand o valoare proprie corespunde unui singur bloc ireductibil,
multiplicitatea ei este exact un patrat perfect.

In cazul de fata, cele noua multiplicitati sunt

\[
1^2,2^2,3^2,4^2,5^2,6^2,3^2,4^2,2^2,
\]

adică pătratele dimensiunilor celor nouă ireps ale lui `2I`.

Deci aceasta proprietate este:

- nu doar observatie;
- ci consecinta structurala a decompunerii reprezentarii regulate.

### C. What exactly means `3` and `3'` localize

Formularea corecta nu este „se localizeaza” in sens vag.

Afirmația exactă este:

- după factorizarea acțiunii prin `A_5`,
- eigenspațiul lui `\Delta_0` pentru valoarea proprie `12-4\phi`, de
  dimensiune 9, se descompune ca trei copii ale irrep-ului `3` al lui `A_5`;
- eigenspațiul pentru valoarea proprie `8+4\phi`, tot de dimensiune 9, se
  descompune ca trei copii ale irrep-ului `3'`.

Aceasta nu este o simpla inevitabilitate dimensionala.

Motiv:

- un spatiu invariant de dimensiune 9 ar putea, a priori, sa se descompuna si
  altfel sub `A_5`;
- scriptul verifica explicit prin caractere că pentru primul eigenspațiu
  apare `3` de trei ori și `3'` deloc, iar pentru al doilea exact invers.

Deci este o consecință reprezentațională verificată, nu doar o observație vagă.

## Additional Exact Features

Rezulta de asemenea:

1. toate multiplicitatile sunt patrate perfecte;
2. suma multiplicitatilor este
   \[
   120 = |2I|;
   \]
3. valorile proprii Galois-conjugate au multiplicitati egale;
4. eigenspatiul pentru
   \[
   12-4\phi
   \]
   contine exact 3 copii ale irrep-ului `3` al lui `A_5`;
5. eigenspatiul pentru
   \[
   8+4\phi
   \]
   contine exact 3 copii ale irrep-ului `3'` al lui `A_5`.

## Closure Decision

`S04` este inchis ca:

- `Computational fact`

in sens tare:

- constructia finita,
- valorile proprii,
- multiplicitatile,
- si structura `\mathbb{Z}[\phi]`

sunt toate verificate explicit.

## What Is Not Claimed

Nu se pretinde aici:

- ca spectrul este dedus conceptual dintr-un theorem analitic inchis;
- ca fiecare proprietate spectrala are deja o interpretare fizica unica.

Se pretinde doar:

- exactitatea calculului finit si a structurii aritmetice asociate.
