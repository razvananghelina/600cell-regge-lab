# Protocol: certificat Coxeter exhaustiv pentru clasa internă 0

Data: 2026-08-21

Status: preregistrat după prior-art gate `295e90d`, înainte de asamblarea lui
`C0` sau diagonalizarea vreunui bloc.

## 1. Ipoteze și obiect

Se folosește exact matricea internă `C0` a schedule-ului lexicografic
`(0,1,2,3)` din protocolul `fdf6f89`, cu:

```text
P(sd K_600), 19680 muchii interne,
tau0 = 0.0102,
m_v = K_v/(8*pi),
coordonate logaritmice ale lungimilor pătrate semnate.
```

Se folosește aceeași tangentă normalizată `n0` și aceeași matrice bordată

```text
K0 = [[C0,n0],[n0^T,0]].
```

Nu se schimbă acțiunea, nu se simetrizează după `H4` înainte de a măsura
defectul de covarianță și nu se încarcă nicio țintă de mod, multiplicitate,
continuum sau fizică.

## 2. Acțiunea Coxeter preregistrată

Din complexul grosier se reconstruiește independent graful celor `14400`
chambers.  Adiacența de culoare `r` unește cele două flags care coincid la
celelalte trei ranguri.  Cele patru aplicații trebuie să fie involuții și să
aibă relațiile Coxeter `(3,3,5)`.

Se fixează chamber-ul de bază cu indice `0`.  Imaginea sa este obținută prin
aplicarea culorilor în ordinea exactă `(0,1,2,3)`.  Acțiunea **stângă** este
unicul automorfism colorat construit prin propagarea condiției

```text
L(c_r x) = c_r L(x)
```

din acea imagine.  Aceasta evită confundarea produsului drept al
generatorilor cu simetria geometrică stângă.

Calculul structural efectuat înainte de spectru și înghețat aici dă:

```text
ordin pe chambers = 30,
ordin pe celulele barycentrice = 30,
ordin pe muchiile interne = 30,
permutarea muchiilor = 656 cicluri de lungime 30,
puncte fixe = 0.
```

Verificatorul trebuie să rederive toate aceste valori.  Fiecare caracter
`k=0,...,29` al lui `C30` are astfel dimensiunea exactă `656`, iar suma
dimensiunilor este `19680`.

## 3. Baza Fourier și exhaustivitatea

Pentru un ciclu ordonat

```text
i_j = L^j(i_0),  j=0,...,29,
```

coloana normalizată a sectorului `k` are coeficienți

```text
q_k(i_j) = exp(-2*pi*i*k*j/30)/sqrt(30).
```

Cele `30*656` coloane trebuie să fie o bază unitară completă.  Nu se păstrează
doar sectoarele mici: se construiesc toate.  Deoarece matricea este reală,
sectoarele `k` și `30-k` trebuie să aibă spectre conjugate/identice.  Se
diagonalizează exhaustiv blocurile `k=0,...,15`; blocurile `1,...,14` au
pondere spectrală doi, iar `0` și `15` pondere unu.  Dimensiunea ponderată a
matricei bordate trebuie să fie

```text
657 + 14*(2*656) + 656 = 19681.
```

Coordonata de bordare și proiecția lui `n0` apar numai în blocul `k=0`.
Proiecțiile tangentei în toate celelalte sectoare trebuie să fie zero în
anvelopa de înmulțire.

## 4. Covarianța numerică și operatorul efectiv testat

Asamblarea `binary64` poate rupe covarianța exactă prin ordinea sumelor.  Nu
se ignoră acest lucru.  Pentru toate puterile `L^r`, `r=0,...,29`, se măsoară

```text
delta_r = ||C0 - P_r^T C0 P_r||_inf.
```

Media de grup

```text
Cbar = (1/30) sum_r P_r^T C0 P_r
```

este exact bloc-diagonală.  Fără a o construi dens, se folosește limita

```text
||C0-Cbar||_2 <= (1/30) sum_r delta_r,
```

care se adaugă explicit anvelopei spectrale.  Covarianța este acceptată doar
dacă această limită este cel mult de 100 de ori anvelopa operatorului local.

Blocurile Fourier sunt blocurile diagonale ale lui `C0` și deci exact
blocurile lui `Cbar`.  Pentru fiecare intrare se acumulează și suma modulelor,
și numărul real de termeni; limita Higham corespunzătoare dă eroarea de
asamblare a blocului.  Norma spectrală a acestei erori este mărginită prin
norma maximă de rând a matricei de anvelope.

## 5. Diagonalizare și criteriul de zero

Fiecare dintre cele 16 blocuri independente, de dimensiune `656` sau `657`,
este diagonalizat complet cu un solver Hermitian dens.  Pentru fiecare
pereche proprie se calculează direct reziduul `2`-norm.  O valoare este
separată de zero numai dacă

```text
|lambda| > 100 * (
    operator_error
  + symmetry_averaging_bound
  + block_roundoff_bound
  + eigenpair_residual
).
```

Nu se folosește reziduul unei rezolvări LU și nu se ajustează numărul de
valori cerute: toate valorile fiecărui bloc sunt calculate.

Verdictul este:

- dacă toate blocurile bordate sunt separate de zero, `K0` este
  numeric certificat nesingular în anvelopa declarată, deci
  `ker(C0)=span(n0)` pentru această clasă și această instanță numerică;
- dacă orice bloc are o valoare compatibilă cu zero, se reconstruiește
  vectorul complet, se raportează caracterul ciclic `k` și reziduul direct,
  iar nucleul suplimentar rămâne candidat până la precizie independentă;
- dacă echivarianța, unitaritatea, dimensiunile sau controalele eșuează,
  verdictul este `COXETER_BLOCK_CONSTRUCTION_INVALID` și nu se interpretează
  spectrul.

Acesta este un certificat numeric cu anvelopă, nu o demonstrație simbolică
asupra numerelor transcendente din acțiune.

## 6. Controale pozitive și negative

Înainte de verdict trebuie să treacă:

1. digestul CSR înghețat al lui `C0`;
2. relațiile Coxeter, ordinul 30, cele 656 cicluri complete și dimensiunea
   Fourier totală;
3. covarianța în anvelopă și identitatea urmelor dintre matrice și blocuri;
4. egalitatea spectrelor blocurilor conjugate, calculată explicit pentru cel
   puțin un reprezentant nenul, nu doar presupusă;
5. reproducerea primelor opt valori Ritz vechi numai ca control post-hoc al
   aceleiași matrice, fără a le folosi la alegerea blocurilor;
6. o corupție preregistrată `+1e-4` în prima intrare diagonală trebuie să
   depășească poarta de covarianță;
7. încercarea de a folosi produsul **drept** al culorilor ca acțiune pe
   celulele de rang trebuie să eșueze testul de well-definedness sau
   covarianță; acesta este controlul convenției.

## 7. Scope

Un rezultat pozitiv certifică numai clasa temporală `0`.  Nu este permisă
extrapolarea la celelalte unsprezece clase și nici interpretarea automată a
valorilor nenule drept gravitoni.  Extinderea la toate clasele se
preregistrează separat numai dacă această reducere trece.

Nu se rulează suita completă și nu se rulează vechiul recensământ rar.
