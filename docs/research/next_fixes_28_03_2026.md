# Next Fixes Prompt

Folosește acest workspace:
`D:\infinity\ToE\science`

Context minim pe care trebuie să-l presupui ca deja stabilit:

- `one_integer_paper.tex` a fost reparat parțial pe linie `theory-first`.
- Unicitatea pentru asignarea `(a,b)` nu mai este doar bounded search:
  există acum un lanț constructiv bazat pe:
  McKay main chain, minimal edge lifts, branch identification, și
  relația Galois `z_b = phi * sigma(z_t)`.
- Sectorul spectral action are acum un nucleu discret exact verificat:
  `c0 = 2640`, `c1 = 14880`, `c2 = 55920`,
  tripla redusă `(11,62,233)`,
  și identitatea
  `2 A1^2 + 1 = 3 A0 A2`.
- Sectorul gravity / scalar response a fost coborât corect la teoremă discretă exactă:
  `gamma_disc = 1`,
  adică `B^+ d0 = d0 Delta0^+`, `h = d0 Phi`, `d1 h = 0`, `d0^+ h = Phi`.
- Paper-ul a fost deja curățat în concluzii / sumar astfel încât:
  unicitatea constructivă este `Derived`,
  prefactorii de gauge `(8/15, 1/3, 2/15)` sunt încă deschiși,
  iar interpretarea continuum PPN nu este revendicată.

Fișiere cheie de citit mai întâi:

- `theory_fix_progress.md`
- `one_integer_paper.tex`
- `spectral_action_discrete_theorem.md`
- `discrete_scalar_response_theorem.md`
- `minimal_edge_lift_theorem.md`
- `reproducible/README.md`
- `reproducible/run_all.py`

Verificatoare deja relevante:

- `reproducible/verify_global_uniqueness_constructive.py`
- `reproducible/verify_branch_identification.py`
- `reproducible/verify_minimal_edge_lifts.py`
- `reproducible/verify_spectral_action.py`
- `reproducible/verify_discrete_scalar_response.py`

Stare de reproducibilitate cunoscută:

- `python reproducible\\run_all.py` trecea ultima dată cu `20/20 PASS`.

Problemele reale rămase, în ordinea recomandată:

1. Spectral action gauge prefactors
   Țintă:
   decide riguros dacă prefactorii
   `(8/15, 1/3, 2/15)`
   pot fi derivați structural din datele discrete, sau trebuie coborâți definitiv
   la `pattern / normalization ansatz`.
   Cerință:
   nu lăsa formulări ambigue în paper.
   Rezultat dorit:
   ori un verifier nou + text întărit, ori downgrade explicit în paper și supplementary.

2. Curățare completă a wording-ului stale din paper și supplement
   Țintă:
   găsește toate locurile unde textul încă afirmă mai mult decât este demonstrat,
   sau încă vorbește în termenii versiunii vechi (`bounded search`, continuum claims,
   spectral action prea tare, PPN prea tare).
   Caută în special:
   `one_integer_paper.tex`
   și
   `one_integer_supplementary.tex`.
   Rezultat dorit:
   formulări consistente peste tot între theorem / derived / structural / pattern / open.

3. Alpha derivation audit
   Țintă:
   separă exact ce parte din derivarea lui `alpha` este teoremă,
   ce parte este identificare geometrică,
   și ce parte este ansatz de normalizare.
   Rezultat dorit:
   o notă scurtă de status + eventual patch în paper ca să nu supra-claim.

4. Gravity wording hardening
   Țintă:
   verifică dacă mai există formulări care sugerează derivarea completă a lui
   `gamma_PPN = 1` sau a unei gravitații 4D neliniare.
   Rezultat dorit:
   toate afirmațiile despre gravity să fie strict compatibile cu:
   teorema discretă exactă + polarizations + open nonlinear completion.

5. Supplement synchronization
   Țintă:
   adu `one_integer_supplementary.tex` la aceeași stare conceptuală cu main paper.
   În special:
   spectral action,
   uniqueness,
   gravity,
   honesty/status paragraphs.

Reguli de lucru:

- Nu reintroduce claim-uri continuum fără verifier sau derivare clară.
- Când o afirmație este doar analogie, spune explicit `pattern`, `structural`, sau `motivates`.
- Când o afirmație este exactă doar pe stratul discret, spune explicit `exact discrete theorem`.
- Nu slăbi ce este deja demonstrat constructiv pe uniqueness.
- Păstrează separarea:
  `exact discrete result` vs `continuum interpretation`.

Ce vreau de la tine în această sesiune:

1. Citește starea actuală din fișierele de mai sus.
2. Identifică cele mai periculoase contradicții sau overclaims rămase.
3. Repare direct textele, nu te opri la analiză.
4. Dacă faci o afirmație teoretică nouă, încearcă să adaugi și un verifier în `reproducible/`.
5. Actualizează `theory_fix_progress.md` la final cu ce ai schimbat.
6. Dacă poți, rulează `python reproducible\\run_all.py`.
7. Dacă PDF-ul nu compilează, spune clar de ce; ultima problemă cunoscută a fost
   MiKTeX local neconfigurat complet.

Formatul ideal al rezultatului:

- o listă scurtă cu ce ai schimbat;
- ce claims au fost întărite;
- ce claims au fost downgradate;
- ce a rămas încă deschis;
- dacă testele / verificatoarele au trecut.

Prioritatea absolută:

Nu căuta să faci paper-ul mai impresionant.
Fă-l intern coerent și defensabil matematic.
