# Protocol adversarial: suportul geometric al eliminării prism-shift

Data: 2026-08-19

Protocolul principal şi artefactul său au fost îngheţate în commit-urile
`e411f93` şi `78fa42c`.  Acest audit nu recalculează acţiunea Regge sau
complementul Schur.  Atacă ipoteza geometrică rămasă: că extensia polytopală
a celor 120 de variaţii de strut este un graf numai în coordonatele interne
puternice, pentru ambele triangulări staircase.

## 1. Reconstrucţie independentă

Pornind numai de la coordonatele 600-cell din `commons/cell600.py`:

1. reconstruieşte produsul grupului binar icosaedral prin înmulţire de
   cuaternioni şi nearest-neighbour unic;
2. reconstruieşte cele cinci clase drepte ale subgrupului binar tetraedral;
3. obţine colorările `even` şi `odd`, ultima prin schimbarea primelor două
   clase;
4. reconstruieşte cele 600 tetraedre ca 4-cliques ale grafului de vârfuri;
5. construieşte fiecare slab staircase direct din ordinea culorilor, fără a
   importa constructorul Regge folosit de verificatoarele îngheţate.

Pentru fiecare paritate se cer exact:

```text
720 muchii jos + 720 sus + 120 poli + 720 diagonale = 2280,
2400 4-simplexe.
```

Fiecare muchie spaţială `{i,j}` trebuie să producă exact o diagonală, de la
vârful de jos cu culoare mai mare la copia de sus a vârfului cu culoare mai
mică.  Nicio diagonală sau pol nu poate fi muchie de frontieră.

## 2. Formula geometrică exactă

Într-un tetraedru regulat centrat, notează cu `q != 1` raportul de scară,
cu `r_i` variaţia pătratului strutului `i`, iar cu `D_ij` variaţia pătratului
diagonalei orientate de la `B_i` la `T_j`.  Din geometria afină trebuie
rederivat exact

```text
r_j-r_i = 2(q-1)(b_j-b_i).s,
D_ij = [q r_j-r_i]/(q-1).
```

Pentru coordonatele logaritmice `z_i` ale magnitudinii polului,
`r_i=-rho*z_i`, graful în coordonatele logaritmice ale diagonalelor este

```text
delta log D_ij = rho*(z_i-q*z_j)/[(q-1)*D_ij].
```

Se cer două controale exacte, raţionale şi nesimetrice, calculate direct din
coordonate, nu din formula ţintă.  Formula cu orientarea inversată trebuie să
eşueze pe cel puţin un control.

## 3. Rang şi suport

Construieşte matricea rară `G` cu 840 rânduri interne (720 diagonale urmate de
120 poli) şi 120 coloane `z`.  Pentru două valori raţionale distincte ale lui
`q` şi coeficienţi nenuli:

- `rank(G)=120` peste `F_101` şi `F_1000003`;
- restricţia la hiperplanul relativ are rang `119`;
- proiecţia colectivă reproduce coeficientul geometric îngheţat
  `-rho/D` pe toate diagonalele şi `1` pe poli;
- toate cele 720 rânduri diagonale aparţin coordonatelor interne non-pole,
  iar cele 120 rânduri pol sunt exact coordonatele slabe;
- suportul pe cele 1440 coordonate de frontieră este exact zero.

Nu se testează o valoare dorită pentru spectru, viteză sau masă.

## 4. Control negativ şi verdict

Controalele negative elimină diagonalele aferente unui vârf şi, separat,
înlocuiesc orientarea geometrică prin cea inversată.  Cel puţin unul dintre
rang, identitatea directă şi suportul complet trebuie să eşueze.

Verdictul

```text
RELATIVE_SHIFT_ELIMINATION_GEOMETRICALLY_CORROBORATED
```

este permis numai dacă toate porţile trec.  Altfel verdictul principal este
redeschis ca

```text
RELATIVE_SHIFT_ELIMINATION_GEOMETRY_OPEN.
```

## 5. Limita interpretării

Un rezultat pozitiv confirmă numai afirmaţia locală, omogenă, la geometria
veche şi impulsul de intrare fixate.  Nu dovedeşte că răspunsul indus pe
frontieră este zero şi nu identifică modurile tensoriale fizice.  Clasificarea
„auxiliar/pseudo-constrângere” rămâne **STRUCTURAL**; propagarea tensorială,
dispersia şi viteza limită rămân **OPEN**.

