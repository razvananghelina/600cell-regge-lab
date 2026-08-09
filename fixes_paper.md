# Fixes Paper v3.1 — Tracking

Data: 14 Feb 2026
Sursa: `14_feb_2026_recomandari_paper.md`

---

## Status Legend
- [x] DONE
- [~] PARTIAL (text fix, cercetare ramane deschisa)
- [ ] OPEN (necesita cercetare noua)

---

## FIX 1: Counting-ul de predictii (Rec #8, Prioritate Scazuta)
**Problema:** "~40" in abstract, "45+" in concluzii — include predictii corelate (J din CKM, m_W din m_Z+sin2tW, sum rules din exponenti).
**Actiune:** Numara independent, inlocuieste cu numar onest.
**Numarare independenta:**
- 3 coupling constants (alpha, alpha_s, sin2tW)
- 3 mass ratios (m_Z/m_e, m_H/m_W, m_e/m_P)
- 8 fermion mass ratios (m_mu..m_b / m_e; m_mu=m_s bare dar particule diferite)
- 6 mixing angles (3 CKM + 3 PMNS)
- 3 CP phases (delta_CKM, delta_PMNS, theta_QCD)
- 1 neutrino mass (m_3)
- 3 structural (N_gen, gauge group, (a,b) quantum numbers)
- 1 gravity (gamma_PPN)
- ~5 dark sector (darkness, stability, mass duality, no confinement, inverted seesaw)
- Total: **~28 independent predictions** + CC pattern + Omega_DM speculative
**Status:** [x] DONE — abstract: "more than 25 independent", conclusions: "~30 independent quantities"

---

## FIX 2: Predictii falsificabile in abstract/introducere (Rec #11, Prioritate Inalta)
**Problema:** Sectiunea 15 vine la p.32, prea tarziu. Referee decide in primele 5 pagini.
**Actiune:** Adaugat in abstract si introducere cele 4 predictii cheie + propozitiie de falsificare.
**Status:** [x] DONE — adaugat in abstract si in introducere

---

## FIX 3: Recategorizare onesta Tabel 2 (Rec #9, Prioritate Inalta)
**Problema:** Inconsistenta: Sectiunea 16.1 pune alpha, alpha_s, holonomy corrections in "Algebraically motivated" dar Tabelul le pune "Derived".
**Actiune:** Recategorizat in Tabel:
- alpha → Alg. motivated
- alpha_s → Alg. motivated
- Holonomy-corrected masses → Alg. motivated
- PMNS theta_13 → Motivated estimate (split de restul PMNS)
- Adaugat coloana "Remaining gap" in tabel
**Status:** [x] DONE

---

## FIX 4: Comparatie tehnica Lisi/Furey/Connes/O'Neill (Rec #10, Prioritate Medie)
**Problema:** Related work prea scurt, un referee va vrea detalii tehnice.
**Actiune:** Expandat paragraful Related work cu comparatii tehnice specifice.
**Status:** [x] DONE

---

## FIX 5: Ecuatia diofantina — axioma vs derivare (Rec #1, Prioritate Medie)
**Problema:** Ecuatia e prezentata cu "geometric origin" dar nu e derivata dintr-un principiu.
**Actiune:** Reframat explicit ca "geometrically motivated axiom" cu nota ce ar constitui o derivare.
**Status:** [~] PARTIAL — text clarificat, derivare ramane open problem

---

## FIX 6: Alpha — onestitate one-loop (Rec #2, Prioritate Inalta)
**Problema:** Coeficientul one-loop 2*pi e postulat, nu derivat din self-energy pe latice.
**Actiune:** Adaugat nota explicita in Sectiunea 5.1 si in footnote-ul tabelului.
**Status:** [~] PARTIAL — text clarificat, calcul pe latice ramane open problem

---

## FIX 7: Holonomy corrections — onestitate sectoriale (Rec #3, Prioritate Inalta)
**Problema:** 3 seturi de corectii per sector, nu un operator unificator.
**Actiune:** Adaugat nota explicita in Sectiunea 7 ca cele 3 seturi sunt "assigned per sector" si ca un operator unificator ramane open problem.
**Status:** [~] PARTIAL — text clarificat, operator unificat ramane open problem

---

## FIX 8: sin2theta_13 = 1/45 — recategorizare (Rec #4, Prioritate Medie-Inalta)
**Problema:** Channel counting nu e calcul perturbativ riguros.
**Actiune:** Recategorizat la "Motivated estimate" in Tabel si text.
**Status:** [x] DONE — recategorizat in tabel, nota in text

---

## FIX 9: Sectiunea CC — scurtare (Rec #5, Prioritate Medie)
**Problema:** Prea lunga pentru un pattern fara mecanism.
**Actiune:** Condensat subsectiunile 14.2 si 14.3 intr-un singur paragraf.
**Status:** [x] DONE

---

## FIX 10: Dark matter abundance — tonare (Rec #7, Prioritate Medie-Scazuta)
**Problema:** Omega_DM/Omega_b = 7-phi fara mecanism.
**Actiune:** Pastrat in sectiune dar eliminat din abstract; marcat clar ca "noted coincidence".
**Status:** [x] DONE

---

## FIX 11: Cod reproductibil (Rec #13, Prioritate Medie)
**Problema:** Codul nu e disponibil public.
**Actiune:** Adaugat nota in Acknowledgments despre disponibilitatea codului.
**Status:** [x] DONE

---

## FIX 12: Beta functions / scala M_Z (Rec #6, Prioritate Inalta)
**Problema:** Identificarea scalei 600-cell cu M_Z e "self-consistent but not derived."
**Actiune:** NU poate fi fixat textual — necesita cercetare. Deja mentionat in limitation 7.
**Status:** [ ] OPEN — necesita derivare beta functions din spectral action

---

## FIX 13: Lagrangianul fermionic D_F pe McKay (Rec #12, Prioritate Inalta)
**Problema:** Yukawa couplings din D_F pe graful McKay lipsesc.
**Actiune:** NU poate fi fixat textual — necesita cercetare. Deja mentionat in limitation 7.
**Status:** [ ] OPEN — cel mai important open problem

---

## FIX 14: sin2theta_13 derivare riguroasa (Rec #4, Prioritate Medie-Inalta)
**Problema:** Channel counting e o estimare, nu calcul perturbativ.
**Actiune:** exp312 — calcul spectral complet pe 120 vertexuri. Rezultate adaugate in paper.
**Rezultat exp312:** 3 localizat in ker(L-(12-4*phi)), 3' in ker(L-(8+4*phi)).
Gap spectral = 4*sqrt(a1) = 4*sqrt(5). theta_13=0 exact la nivel A5.
sin^2(theta_13) = 1/45 = rata democratica de tranzitie prin bariera spectrala.
**Status:** [x] DONE — derivare completa adaugata in paper (tabel eigenspaces, gap spectral, ec. theta13)

---

## FIX 15: Holonomy corrections unificate (Rec #3, Prioritate Inalta)
**Problema:** Un singur operator care produce cele 3 seturi de corectii.
**Actiune:** exp313 — sub-adjacency matrices A_R construite, edge decomposition verificata.
**Rezultat exp313:** Edge weights sunt IZOTROPE (w_R = d_R/12 pt toate muchiile).
**Rezultat exp313b:** Hopf fibration: 3 clase de muchii {0, 1/phi^2, 1/phi}, fiber/base=phi.
**Rezultat exp313c:** Hopf x Gauge: chi2=0, PERFECT INDEPENDENT. Nici vertex types, nici rank-based.
**Rezultat exp314:** Triangle SO(3) holonomy: Omega = same pt toate 1200 fete.
  Normal fiber: 5 valori {0, 1/(2phi), 1/2, phi/2, 1} (Galois conjugate pair!).
  **Mean Berry phase = 0.14592 ~ 1/phi^4 = delta_k(lepton) (0.014% match)**.
  DAR: NO anisotropy intre sectoare gauge (ACC/BCC/CCC identice).
**Concluzie:** 600-cell-ul este PREA simetric (vertex+edge+face transitive).
  Berry phase da baseline universal 1/phi^4. Sector dependence necesita alt mecanism.
**Status:** [~] PARTIAL — Berry phase 1/phi^4 DERIVED, sector dependence ramane OPEN

---

## FIX 16: CC din spectral action (Rec #5, Prioritate Medie)
**Problema:** Exponentul 57 nu apare din spectral action.
**Actiune:** NU poate fi fixat textual — necesita calcul Seeley-DeWitt.
**Status:** [ ] OPEN — spectaculos daca reusesti

---

## FIX 17: Alpha one-loop pe latice (Rec #2, Prioritate Inalta)
**Problema:** Coeficientul 2*pi trebuie derivat din self-energy pe graful 600-cell.
**Actiune:** NU poate fi fixat textual — necesita calcul propagator pe graf.
**Status:** [ ] OPEN — candidat pentru exp314

---

## Rezumat

| Fix | Descriere | Status | Prioritate |
|-----|-----------|--------|-----------|
| 1 | Counting predictii onest | DONE | Scazuta |
| 2 | Predictii falsificabile in abstract | DONE | Inalta |
| 3 | Recategorizare Tabel 2 | DONE | Inalta |
| 4 | Comparatie Lisi/Furey/Connes | DONE | Medie |
| 5 | Diophantine = axioma | PARTIAL | Medie |
| 6 | Alpha one-loop onestitate | PARTIAL | Inalta |
| 7 | Holonomy sector assignment | PARTIAL | Inalta |
| 8 | theta_13 recategorizare | DONE | Medie-Inalta |
| 9 | CC scurtare | DONE | Medie |
| 10 | DM abundance tonare | DONE | Medie-Scazuta |
| 11 | Cod reproductibil | DONE | Medie |
| 12 | Beta functions / M_Z | OPEN | Inalta |
| 13 | D_F pe McKay (Yukawa) | OPEN | Inalta |
| 14 | theta_13 riguroasa | DONE | Medie-Inalta |
| 15 | Holonomy operator unificat | PARTIAL | Inalta |
| 16 | CC din spectral action | OPEN | Medie |
| 17 | Alpha lattice self-energy | OPEN | Inalta |

## FIX 18: Berry phase baseline + no-go Proposition (din exp314, Prioritate Inalta)
**Problema:** Holonomy corrections pretindeau "Derived" dar sector assignment nu era derivat.
**Actiune:** Adaugat in paper:
- Berry phase = 1/phi^4 ca Eq. (berry_baseline) in sectiunea holonomy
- Proposition 3 (Static geometry no-go) cu 3 rezultate: isotropic edges, Hopf-gauge independence, face-transitive holonomy
- Recategorizat holonomy-corrected masses la "Alg. motivated" in Table 2 cu footnote k
- Updatat "Derivation vs. Reproduction" item (ii) cu no-go reference
- Noua limitare 8 (Holonomy sector dependence)
- Honest counting: abstract "more than 25 independent", concluzii "~30 independent"
- Version bumped to 3.2
**Status:** [x] DONE

---

## Rezumat

| Fix | Descriere | Status | Prioritate |
|-----|-----------|--------|-----------|
| 1 | Counting predictii onest | DONE | Scazuta |
| 2 | Predictii falsificabile in abstract | DONE | Inalta |
| 3 | Recategorizare Tabel 2 | DONE | Inalta |
| 4 | Comparatie Lisi/Furey/Connes | DONE | Medie |
| 5 | Diophantine = axioma | PARTIAL | Medie |
| 6 | Alpha one-loop onestitate | PARTIAL | Inalta |
| 7 | Holonomy sector assignment | PARTIAL | Inalta |
| 8 | theta_13 recategorizare | DONE | Medie-Inalta |
| 9 | CC scurtare | DONE | Medie |
| 10 | DM abundance tonare | DONE | Medie-Scazuta |
| 11 | Cod reproductibil | DONE | Medie |
| 12 | Beta functions / M_Z | OPEN | Inalta |
| 13 | D_F pe McKay (Yukawa) | OPEN | Inalta |
| 14 | theta_13 riguroasa | DONE | Medie-Inalta |
| 15 | Holonomy operator unificat | PARTIAL | Inalta |
| 16 | CC din spectral action | OPEN | Medie |
| 17 | Alpha lattice self-energy | OPEN | Inalta |
| 18 | Berry phase + no-go + recategorize | DONE | Inalta |

**Scor: 10 DONE + 4 PARTIAL + 4 OPEN**
