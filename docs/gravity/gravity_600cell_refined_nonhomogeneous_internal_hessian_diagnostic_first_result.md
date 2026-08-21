# Primul rezultat al diagnosticului Hessianului intern neomogen

Data: 2026-08-21

Status: **12/13 CONTROL FAILED; diagnosticul nu este acceptat integral.**

## Proveniență înghețată

- protocol diagnostic înainte de execuție: `cdd0f69`;
- implementare înregistrată înainte de execuție: `0118b41`;
- artefact:

  ```text
  reproducible/gravity_600cell_refined_nonhomogeneous_internal_hessian_diagnostic.json
  SHA-256 ca0729becf027aa0b4181fe86abf2daede6ff769738a16683eab42100bc35e80
  ```

A fost asamblată și factorizată numai clasa lexicografică `0`.  Nu s-a rulat
suita completă și nu s-a repetat recensământul celor douăsprezece clase.

## Ce a fost reprodus

- digestul CSR al matricei clasei 0 coincide exact cu prima execuție;
- diferența înghețată `P64-PAGG` este reprodusă exact:

  ```text
  max = 6.293703336268663e-10 la intrarea (1,1);
  fracție din poarta veche = 2.6368154.
  ```

- regula veche de corupție selectează coeficient zero exact la indicii
  forward `6,7,12,14`;
- regula preregistrată „prima incidență cu |gradient arie| > 1e-20” produce
  o corupție detectabilă în toate cele 12 reprezentante;
- toate cele opt rezolvări LU inițiale au eroare backward componentwise între
  aproximativ `1.25e-12` și `8.19e-12`, sub poarta `2.18e-10`;
- o corecție iterativă le aduce tipic la ordinul `1e-16`, deși reziduul brut
  raportat doar la membrul drept poate rămâne `1e-9`--`1e-7`;
- primele opt valori Ritz din artefactul vechi sunt reproduse.

Aceste rezultate susțin diagnosticul matematic că reziduul brut al unei
rezolvări aproape singulare nu este o incertitudine de valoare proprie.
Acesta nu trebuie înmulțit cu norma matricei și adăugat reziduului Ritz.

## Framing attack asupra diagnosticului pullback

Clasificatorul preregistrat a returnat mecanic

```text
AGGREGATE_FINITE_DIFFERENCE_CONTROL_LOCALIZED
```

din pattern-ul de compatibilitate `(P64=PLD, P64!=PAGG, PLD!=PAGG)`.  Nu
acceptăm această etichetă drept cauză stabilită.  Compatibilitatea unor
intervale cu anvelope foarte diferite nu este tranzitivă, iar valorile
centrale spun contrariul la intrarea dominantă `(1,1)`:

```text
|P64-PAGG| = 6.2937e-10,
|P64-PLD|  = 6.3036e-10,
|PLD-PAGG| = 9.87e-13 aproximativ.
```

`PLD-PAGG` eșuează la o altă intrare mică, `(2,9)`, cu fracția de poartă
`1.7095`, deși diferența absolută este numai `2.91e-16`.  Prin urmare:

- **PATTERN:** valorile centrale ale sumei directe `long double` urmăresc
  controlul agregat la discrepanța dominantă și sugerează o problemă de
  acumulare/pullback `binary64`;
- **OPEN:** sursa exactă a discrepanței entrywise; ierarhia relațională din
  protocol nu o poate decide;
- eticheta mecanică din JSON rămâne păstrată ca rezultat al protocolului, dar
  nu este promovată în ledger.

Aceasta este o eroare de framing a diagnosticului nostru, nu o modificare a
acțiunii fizice.

## Controlul spectral care a eșuat

Cele două rulări de câte 32 de vectori Ritz au reziduuri directe foarte mici
și fiecare valoare observată este separată de zero față de eroarea de
asamblare.  Totuși listele nu coincid entrywise.  Gruparea valorilor centrale
arată:

```text
rulare 1: multiplicități observate 4, 9, 13, 6;
rulare 2: multiplicități observate 4, 9, 16, 3.
```

Primele două clustere coincid; diferența apare când fereastra de 32 taie
clusterele al treilea și al patrulea.  Este plauzibil că ARPACK a returnat
baze diferite dintr-un multiplet degenerat, dar protocolul cerea egalitate a
listelor complete.  Așadar acesta este un **PATTERN**, nu o scuză
retroactivă, iar verdictul rămâne

```text
DIAGNOSTIC_INVALID.
```

## Verdict și continuare admisă

- **DERIVED COMPUTATIONAL:** matricea veche este reproductibilă; controlul de
  corupție avea o eroare de selecție și regula corectată o repară; rezolvările
  LU sunt backward-stable după metrica potrivită.
- **PATTERN:** clusterele moi apar cu multiplicități `4,9,16,...`, compatibile
  cu multiplete de simetrie, dar acest lucru nu este încă demonstrat.
- **OPEN:** nucleul complet al lui `C_s`.
- **NOT TESTED:** gravitoni, propagare, viteză, `G` sau scară Planck.

O nouă cerere ARPACK cu o fereastră mai mare ar putea evita tăierea
multipletului, dar tot nu ar certifica exhaustivitatea.  Continuarea utilă nu
este repetarea oarbă a diagonalizării, ci o descompunere exactă după o
simetrie finită a carrier-ului (de exemplu elementul Coxeter de ordin 30),
urmată de diagonalizarea exhaustivă a tuturor blocurilor.  Aceasta necesită un
protocol separat și nu este executată în această notă.
