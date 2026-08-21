# Protocol diagnostic: primul Hessian intern neomogen rafinat

Data: 2026-08-21

Status: preregistrat după artefactul eșuat `8f7113c`, înaintea oricărui
recalcul spectral.

## 1. Scop și limită

Prima execuție completă a verificatorului
`verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py` este
înghețată cu verdictul corect

```text
LOCAL_EXTENSION_INVALID
```

deoarece două controale au eșuat.  Aceasta nu este nici dovadă de nucleu
suplimentar, nici dovadă de nesingularitate bordată.  Protocolul de față
diagnostichează cele două controale și scara reziduurilor pe **o singură**
clasă temporală.  Nu repetă recensământul de două ore, nu încarcă moduri
spațiale sau ținte de continuum și nu produce un verdict fizic.

Carrier-ul, acțiunea, masele și coordonatele sunt exact cele din protocolul
`fdf6f89`.  Clasa aleasă este ordinea lexicografică

```text
(0,1,2,3), schedule_index = 0,
```

aleasă înainte de spectru.  Celelalte unsprezece clase nu sunt factorizate.

## 2. Diagnosticul pullback-ului `10 x 10`

Pentru clasa 0 se construiesc trei matrici fără a schimba formula acțiunii:

1. `P64`: asamblarea rară `binary64`, urmată de pullback-ul pe cele zece
   orbite, exact ca în prima execuție;
2. `PLD`: suma directă a tuturor incidențelor locale direct în cele `100`
   intrări de orbită, cu stencils la `140` zecimale convertite în
   `complex long double`; aceasta nu construiește matricea rară și nu face
   înmulțirea cu indicatorii de orbită;
3. `PAGG`: Hessianul acțiunii agregate, diferențiat independent la `100` și
   `140` zecimale, exact ca în controlul înghețat.

Se raportează fiecare matrice, toate cele trei diferențe entrywise, poziția și
semnul maximului, precum și anvelopele forward.  Pentru `PLD` anvelopa include
diferența stencil-urilor între precizii și limita Higham calculată din numărul
real de termeni adunați; pentru `P64` se păstrează anvelopa asamblării
înghețate; pentru `PAGG` se păstrează anvelopa diferențierii agregate.  O
egalitate este acceptată numai dacă

```text
max |X-Y| <= 100 * (err_X + err_Y).
```

Ierarhia de diagnostic este fixată astfel:

- dacă `PLD=PAGG`, dar `P64` diferă, problema este acumularea/anvelopa
  `binary64`;
- dacă `P64=PLD`, dar `PAGG` diferă, problema este controlul agregat prin
  diferențe finite;
- dacă `P64=PAGG`, dar `PLD` diferă, suma directă nouă este greșită;
- dacă nu există o singură pereche compatibilă sau toate trei diferă,
  formula rămâne **OPEN** și nicio corecție nu este acceptată;
- dacă toate trei coincid, eșecul vechi a fost doar o anvelopă prea strânsă,
  dar aceasta se afirmă numai dacă noua anvelopă este derivată, nu mărită
  până la rezultat.

Digestul CSR al matricei reconstruite trebuie să coincidă exact cu digestul
clasei 0 din artefactul înghețat; altfel diagnosticul se oprește.

## 3. Controlul de corupție

Regula veche alegea prima incidență cu rând și coloană interne.  Artefactul
înghețat arată că gradientul ariei acelei incidențe este zero pentru indicii
forward `6, 7, 12, 14` și nenul pentru `0, 1, 2, 3, 4, 5, 8, 10`.
Diagnosticul trebuie mai întâi să reproducă exact acel set de eșecuri; aceasta
este controlul negativ cunoscut.

Regula corectată este declarată înainte de execuție:

> Alege lexicografic prima incidență cu rând și coloană interne pentru care
> gradientul ariei calculat la 140 zecimale are modul strict mai mare decât
> `1e-20`.

Pragul separă un coeficient geometric de ordin unitate de zero numeric; nu
este ales dintr-un spectru.  Corupția rămâne exact `1e-4` într-un singur
coeficient al derivatei unghiului.  Regula se verifică pe toate cele
12 reprezentante de inversare temporală, fără asamblări sau factorizări
complete.  Ea trece numai dacă schimbarea prezisă depășește de 100 de ori
anvelopa locală înghețată în fiecare caz.

## 4. Diagnosticul spectral pe clasa 0

Se factorizează o singură dată matricea bordată

```text
K0 = [[C0,n0],[n0^T,0]].
```

Pentru aceleași opt membre drepte deterministe din prima execuție se
raportează separat:

```text
r_raw = ||Kx-b||_inf / ||b||_inf,
eta    = max_i |Kx-b|_i / (|K||x|+|b|)_i.
```

`r_raw` este o eroare forward raportată la un membru drept de ordin unitate;
nu este, prin definiție, o incertitudine de valoare proprie.  `eta` este
eroarea backward componentwise relevantă rezolvării liniare.  Se fac exact
trei corecții iterative cu aceeași factorizare și se publică întregul șir;
îmbunătățirea nu este presupusă.

O rezolvare este etichetată numeric backward-stable numai dacă

```text
eta <= 100 * gamma_q,
gamma_q = q*u/(1-q*u),
```

unde `q` este numărul maxim de termeni nenuli într-un rând al matricei
bordate și `u` este unit roundoff `binary64`.  Acest test nu este transformat
într-o limită asupra spectrului.

Cu aceeași factorizare se cer, înainte de rezultat, cele 32 de valori Ritz
cele mai apropiate de zero, în două rulări cu toleranțe `1e-10` și `1e-12` și
vectori inițiali determinist diferiți.  Pentru fiecare se verifică direct

```text
||K v - lambda v||_2.
```

O valoare Ritz observată este separată de zero față de eroarea operatorului
numai dacă

```text
|lambda| > 100 * (operator_error + Ritz_residual).
```

Chiar dacă toate cele 32 trec, această execuție nu certifică faptul că ARPACK
nu a omis un mod nul.  Prin urmare rezultatul maxim admis aici este:

- reziduul brut vechi este sau nu o categorie de eroare nepotrivită;
- clusterul moale observat este sau nu rezolvat față de eroarea de asamblare;
- ce certificat suplimentar de exhaustivitate ar fi necesar pentru nucleu.

Nu se acceptă verdictul `ker(C_s)=span(n_s)` din acest diagnostic.

## 5. Controale și rezultat permis

- **CONTROL pozitiv:** digestul matricei și cele opt valori Ritz vechi ale
  clasei 0 trebuie reproduse în anvelopele declarate;
- **CONTROL negativ:** regula veche de corupție trebuie să reproducă cele
  patru selecții nule din artefact;
- schimbarea preciziei, a normei sau a numărului de vectori Ritz nu poate
  schimba retroactiv protocolul complet `fdf6f89`;
- dacă oricare control cunoscut nu este reprodus, verdictul este
  `DIAGNOSTIC_INVALID`;
- altfel diagnosticul poate localiza o eroare de implementare sau poate lăsa
  problema `OPEN`, dar nu poate declara fizică nouă.

Nu se rulează suita completă.  Se rulează numai verificatorul diagnostic și
auditul static de înregistrare/duplicate al `run_all.py`.
