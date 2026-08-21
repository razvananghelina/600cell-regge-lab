# Prior-art gate: blocuri Coxeter pentru Hessianul intern rafinat

Data: 2026-08-21

Status: completat după eșecul diagnostic `2a7064a` și înainte de construirea
permutării sau a vreunui bloc spectral.

## 1. Obiectul exact

Se păstrează neschimbate carrier-ul, acțiunea, masele, coordonatele și clasa
temporală lexicografică din protocoalele anterioare.  Obiectul este matricea
real-simetrică internă

```text
C0 in R^(19680 x 19680)
```

și bordarea ei cu tangenta analitică de durată.  Nu se schimbă matricea și nu
se caută o altă soluție.  Se caută doar o bază unitară exhaustivă în care
aceeași matrice se descompune în blocuri.

Carrier-ul barycentric are acțiunea exactă a grupului Coxeter `H4`.  Se alege
elementul Coxeter obținut din produsul celor patru reflexii colorate în
ordinea fixă `(0,1,2,3)`.  Acesta trebuie să aibă ordin exact `30`.  Acțiunea
stângă pe chambers se reconstruiește drept automorfismul colorat care trimite
chamber-ul de bază în imaginea lui Coxeter și comută cu cele patru adiacențe
colorate.  Din ea se induce mecanic o permutare a celulelor barycentrice, a
celor două straturi și a celor `19680` muchii interne.

Dacă `C0` comută cu această permutare, reprezentarea grupului ciclic `C30`
descompune exact spațiul complexificat în cele 30 de caractere Fourier.
Fiecare ciclu de lungime `L|30` contribuie o direcție la caracterul `k` dacă

```text
k*L = 0 mod 30.
```

Aceasta este o schimbare de bază unitară, nu o selecție de moduri.  Suma
dimensiunilor tuturor sectoarelor trebuie să fie exact `19680`; coordonata de
bordare intră numai în sectorul invariant `k=0`.

## 2. KNOWN din literatura primară

- Descompunerea bloc a operatorilor invarianți sub un grup finit este
  reprezentare standard, nu o idee fizică nouă.  Un tratament explicit al
  reducerii operatorilor/matricelor prin reprezentări de grup este F.
  Vallentin, [Symmetry in semidefinite
  programs](https://arxiv.org/abs/0706.4233).
- Pentru matrici geometrice de rigiditate, aceeași idee este demonstrată
  explicit de B. Schulze,
  [Block-diagonalized rigidity matrices of symmetric frameworks and
  applications](https://arxiv.org/abs/0906.3377).
- Simetriile și modurile nule ale Hessienilor Regge pe fundaluri plate sau
  curbate sunt cunoscute; interpretarea lor rămâne cea din B. Bahr și B.
  Dittrich, [(Broken) Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670), și P. A. Hoehn,
  [Canonical linearized Regge Calculus](https://arxiv.org/abs/1411.5672).
- Evoluția omogenă a 600-cell-ului este cunoscută, inclusiv oprirea cauzală:
  A. De Felice și E. Fabri,
  [Singularities of the closed RW metric in Regge
  Calculus](https://arxiv.org/abs/gr-qc/0106077).

Aceste surse justifică metoda și avertismentele.  Căutarea nu a găsit
descompunerea Coxeter a Hessianului intern complet pentru
`P(sd K_600) x I`, cu masele selectate aici.  Absența din căutare nu dovedește
noutatea; aceasta rămâne **OPEN**.

## 3. Ce ar fi redundant și ce nu

Este redundant să „descoperim” că simetria bloc-diagonalizează un operator
invariant.  Conținutul calculului nostru este numai:

1. dacă matricea concretă construită de noi este exact covariantă sub
   elementul Coxeter indus combinatoric;
2. dacă transformarea ciclică acoperă toate cele `19680` direcții;
3. nulitatea exhaustivă a fiecărui bloc al acestei matrice concrete.

Nu se va prezenta metoda Fourier drept fizică nouă.  Un eventual rezultat
nou poate fi doar spectrul/nucleul instanței exacte, cu noutatea externă
etichetată **OPEN** până la review dedicat.

## 4. Framing attack

Un singur element Coxeter nu separă în general toate ireprezentările `H4` și
nu explică multiplicitățile fizic.  El este însă suficient pentru un test de
nulitate dacă și numai dacă baza Fourier este completă și `C0` comută cu
permutarea: uniunea spectrelor celor 30 de blocuri este atunci spectrul
complet, indiferent că unele blocuri mai conțin mai multe ireprezentări.

Riscurile care trebuie controlate înainte de spectru sunt:

- confundarea acțiunii drepte a generatorilor de chamber cu acțiunea stângă
  geometrică;
- o permutare care nu este bine definită pe reziduurile/celulele de rang fix;
- o convenție Fourier greșită;
- tăierea ciclurilor de lungime proprie `L<30`;
- declararea exhaustivității din câteva valori proprii ale fiecărui bloc.

De aceea protocolul ulterior trebuie să reconstruiască matricea din toate
blocurile sau să verifice Parseval și toate dimensiunile, să diagonalizeze
dense **fiecare** bloc și să testeze direct toate perechile proprii relevante.

## 5. Statut înainte de calcul

- **KNOWN:** teoria reducerii prin simetrie finită și formalismul Hessianului
  Regge.
- **CONTROL:** relațiile Coxeter colorate `(3,3,5)` și carrier-ul barycentric
  sunt deja certificate în repository, dar acțiunea pe muchiile slab-ului va
  fi reconstruită independent în noul verificator.
- **OPEN:** echivarianța exactă a `C0`, dimensiunile sectoarelor și nulitatea
  lor.
- **NOT TESTED:** identificarea drept graviton, dispersia, `c`, `G`, timpul sau
  scara Planck.
