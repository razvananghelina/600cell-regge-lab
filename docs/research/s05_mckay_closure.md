# S05 McKay Closure

Scop: sa inchidem `S05` strict pe nucleul standard si exact.

## Exact Claim

Pentru grupul binar icosaedral `2I`:

1. graful McKay al reprezentarii definitorii de dimensiune 2 este diagrama
   Dynkin afină `\widetilde{E}_8`;
2. cele 9 ireps au dimensiuni
   \[
   1,2,3,4,5,6,3,4,2.
   \]

## Status

- `Theorem` for the McKay correspondence statement;
- `Computational confirmation` for the explicit local realization used here.

## Standard Theorem Layer

Acesta este rezultat standard din clasificarea subgrupurilor finite ale lui
`SU(2)` și corespondența McKay:

- finite subgroups of `SU(2)` are matched with affine ADE diagrams;
- binary icosahedral group `2I` corresponds to affine `E_8`.

Prin urmare:

\[
2I \longleftrightarrow \widetilde{E}_8.
\]

Aceasta parte nu este o descoperire a framework-ului; este matematică standard.

## Local Explicit Confirmation

Scriptul

- [reproducible/verify_mckay_chirality.py](D:\infinity\ToE\science\reproducible\verify_mckay_chirality.py)

verifica explicit:

1. tabelul de caractere al lui `2I`;
2. ortogonalitatea caracterelor;
3. faptul că
   \[
   \sum_i d_i^2 = 120 = |2I|;
   \]
4. faptul că dimensiunile ireps sunt
   \[
   1,2,2,3,3,4,4,5,6
   \]
   în ordonarea internă a scriptului;
5. faptul că graful McKay obținut din produsul cu reprezentarea de dimensiune 2
   are 9 noduri, 8 muchii, este arbore, și are lungimile de brațe
   \[
   (1,2,5),
   \]
   adică exact forma afină `E_8`.

În paper, ordonarea dimensiunilor este scrisă ca

\[
1,2,3,4,5,6,3,4,2,
\]

adică aceeași multime de dimensiuni, doar în ordonarea McKay folosită acolo.

## Closure Decision

`S05` este închis ca:

- `Theorem`:
  `2I` has affine `E_8` McKay graph.
- `Computational confirmation`:
  the explicit character-table construction used in this workspace reproduces
  the expected graph and irrep dimensions.

## What Is Not Claimed

Nu se pretinde aici:

- că deja am derivat conținut fizic din `E_8`;
- că orice proprietate ulterioară din scriptul `verify_mckay_chirality.py`
  intră automat în exact-core.

În exact-core intră doar:

- identificarea standard `2I -> \widetilde{E}_8`;
- datele exacte ale ireps-urilor și ale grafului.
