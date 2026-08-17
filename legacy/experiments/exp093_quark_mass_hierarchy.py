"""
EXP-093: Ierarhia Maselor Quarkilor din Geometria 600-Cell
==========================================================

OBIECTIV:
=========
Investigam daca mecanismul holonomic care explica masele leptonice
(phi^11, phi^6, phi^17) se aplica si la quarks.

OBSERVATIE CHEIE din exp092:
============================
Strange quark are n ≈ 11, ACELASI ca muonul!
Aceasta nu poate fi coincidenta.

IPOTEZE DE TESTAT:
==================
1. Quarkii urmeaza acelasi pattern phi^n ca leptonii?
2. Diferenta vine din factorul de culoare (3)?
3. Corectiile QCD modifica exponentii?
4. Exista o structura geometrica pentru quarks pe 600-cell?

MASE QUARKS (la scala ~2 GeV, MS-bar):
======================================
up:     2.16 +/- 0.49 MeV
down:   4.67 +/- 0.48 MeV
strange: 93.4 +/- 8.6 MeV
charm:  1.27 +/- 0.02 GeV
bottom: 4.18 +/- 0.03 GeV
top:    172.69 +/- 0.30 GeV (masa pol)
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-093: IERARHIA MASELOR QUARKILOR")
print("=" * 70)

phi = PHI
log_phi = np.log(phi)

# ============================================================
# PASUL 1: DATELE EXPERIMENTALE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 1: MASELE QUARKILOR (PDG 2024)")
print("-" * 70)

# Mase in MeV (scala ~2 GeV pentru quarks usori, masa pol pentru top)
quarks = {
    'up':      {'mass': 2.16, 'err': 0.49, 'charge': 2/3},
    'down':    {'mass': 4.67, 'err': 0.48, 'charge': -1/3},
    'strange': {'mass': 93.4, 'err': 8.6, 'charge': -1/3},
    'charm':   {'mass': 1270, 'err': 20, 'charge': 2/3},
    'bottom':  {'mass': 4180, 'err': 30, 'charge': -1/3},
    'top':     {'mass': 172690, 'err': 300, 'charge': 2/3},
}

# Mase leptoni pentru comparatie
leptons = {
    'electron': {'mass': 0.511, 'charge': -1},
    'muon':     {'mass': 105.66, 'charge': -1},
    'tau':      {'mass': 1776.8, 'charge': -1},
}

m_e = leptons['electron']['mass']

print("Quarks:")
for name, data in quarks.items():
    print(f"  {name:8s}: {data['mass']:10.2f} MeV (Q = {data['charge']:+.2f})")

print("\nLeptons (pentru comparatie):")
for name, data in leptons.items():
    print(f"  {name:8s}: {data['mass']:10.2f} MeV")

# ============================================================
# PASUL 2: EXPONENTII phi^n PENTRU QUARKS
# ============================================================
print("\n" + "-" * 70)
print("PASUL 2: EXPONENTII phi^n PENTRU QUARKS")
print("-" * 70)

print("""
Formula: m = m_e * phi^n
Deci: n = log(m/m_e) / log(phi)
""")

print(f"{'Quark':<10} {'m (MeV)':<12} {'n_exact':<10} {'n_round':<8} {'m_pred':<12} {'Err %':<10}")
print("-" * 70)

quark_exponents = {}
for name, data in quarks.items():
    m = data['mass']
    n_exact = np.log(m / m_e) / log_phi
    n_round = round(n_exact)
    m_pred = m_e * phi**n_round
    err = abs(m_pred - m) / m * 100
    quark_exponents[name] = {'n_exact': n_exact, 'n_round': n_round, 'err': err}
    print(f"{name:<10} {m:<12.2f} {n_exact:<10.2f} {n_round:<8} {m_pred:<12.2f} {err:<10.1f}")

# ============================================================
# PASUL 3: COMPARATIE CU LEPTONII
# ============================================================
print("\n" + "-" * 70)
print("PASUL 3: COMPARATIE CU LEPTONII")
print("-" * 70)

print(f"{'Particula':<12} {'Tip':<8} {'n':<6} {'Observatie':<30}")
print("-" * 60)
print(f"{'electron':<12} {'lepton':<8} {0:<6} {'Baza':<30}")
print(f"{'muon':<12} {'lepton':<8} {11:<6} {'5 + 6 (diametru + decagoane)':<30}")
print(f"{'tau':<12} {'lepton':<8} {17:<6} {'11 + 6 (muon + decagoane)':<30}")
print()
print(f"{'up':<12} {'quark':<8} {quark_exponents['up']['n_round']:<6} {'dim(Im(H)) = 3 ?':<30}")
print(f"{'down':<12} {'quark':<8} {quark_exponents['down']['n_round']:<6} {'diametru = 5 ?':<30}")
print(f"{'strange':<12} {'quark':<8} {quark_exponents['strange']['n_round']:<6} {'ACELASI ca muon!':<30}")
print(f"{'charm':<12} {'quark':<8} {quark_exponents['charm']['n_round']:<6} {'~tau (17) - 1?':<30}")
print(f"{'bottom':<12} {'quark':<8} {quark_exponents['bottom']['n_round']:<6} {'17 + 2?':<30}")
print(f"{'top':<12} {'quark':<8} {quark_exponents['top']['n_round']:<6} {'Exceptional':<30}")

# ============================================================
# PASUL 4: STRUCTURA PE GENERATII
# ============================================================
print("\n" + "-" * 70)
print("PASUL 4: STRUCTURA PE GENERATII")
print("-" * 70)

print("""
IPOTEZA: Quarkii si leptonii din aceeasi generatie au exponenti legati.

Generatia 1: electron (0), up (3), down (5)
Generatia 2: muon (11), charm (16), strange (11)
Generatia 3: tau (17), top (26), bottom (19)
""")

generations = [
    {'name': 'Gen 1', 'lepton': ('electron', 0), 'up_type': ('up', 3), 'down_type': ('down', 5)},
    {'name': 'Gen 2', 'lepton': ('muon', 11), 'up_type': ('charm', 16), 'down_type': ('strange', 11)},
    {'name': 'Gen 3', 'lepton': ('tau', 17), 'up_type': ('top', 26), 'down_type': ('bottom', 19)},
]

print(f"{'Generatie':<12} {'Lepton (n)':<15} {'Up-type (n)':<15} {'Down-type (n)':<15}")
print("-" * 60)
for gen in generations:
    l_name, l_n = gen['lepton']
    u_name, u_n = gen['up_type']
    d_name, d_n = gen['down_type']
    print(f"{gen['name']:<12} {l_name} ({l_n}){'':<5} {u_name} ({u_n}){'':<5} {d_name} ({d_n})")

print("""
OBSERVATII:
1. Strange (n=11) = Muon (n=11) - IDENTIC!
2. Down (n=5) = Diametru 600-cell
3. Up (n=3) = dim(Im(H)) sau 3 culori?
""")

# ============================================================
# PASUL 5: DIFERENTE INTRE GENERATII
# ============================================================
print("\n" + "-" * 70)
print("PASUL 5: DIFERENTE INTRE GENERATII")
print("-" * 70)

print("Diferente de exponenti intre generatii consecutive:")
print()

# Leptoni
print("LEPTONI:")
print(f"  muon - electron: 11 - 0 = 11 = 5 + 6")
print(f"  tau - muon: 17 - 11 = 6 (decagoane)")

# Up-type quarks
print("\nUP-TYPE QUARKS:")
print(f"  charm - up: 16 - 3 = 13")
print(f"  top - charm: 26 - 16 = 10 (decagon complet!)")

# Down-type quarks
print("\nDOWN-TYPE QUARKS:")
print(f"  strange - down: 11 - 5 = 6 (decagoane!)")
print(f"  bottom - strange: 19 - 11 = 8 (varfuri 16-cell!)")

# ============================================================
# PASUL 6: IPOTEZA CULORILOR
# ============================================================
print("\n" + "-" * 70)
print("PASUL 6: IPOTEZA CULORILOR (FACTOR 3)")
print("-" * 70)

print("""
IPOTEZA: Quarkii au 3 stari de culoare, deci masa efectiva
ar putea fi modificata cu un factor legat de 3.

Daca m_quark = m_lepton * phi^(delta_n), unde delta_n vine din culoare?
""")

# Comparatie strange-muon
m_strange = quarks['strange']['mass']
m_muon = leptons['muon']['mass']
ratio_strange_muon = m_strange / m_muon

print(f"m_strange / m_muon = {ratio_strange_muon:.4f}")
print(f"Daca ar fi identice geometric: raport = 1")
print(f"Raport actual ~ {ratio_strange_muon:.2f} ~ 1 (aproape!)")

# Ce factor ar trebui?
delta_n = np.log(ratio_strange_muon) / log_phi
print(f"Delta_n pentru strange vs muon: {delta_n:.3f}")
print(f"Deci strange ~ muon * phi^{delta_n:.2f}")

# Verificare cu 3 culori
print(f"\nDaca includem 3 culori:")
print(f"  phi^0.5 = {phi**0.5:.4f} ~ sqrt(phi)")
print(f"  3^(1/3) = {3**(1/3):.4f}")
print(f"  log(3)/log(phi) = {np.log(3)/log_phi:.4f} ~ 2.28")

# ============================================================
# PASUL 7: PATTERN PENTRU UP-TYPE VS DOWN-TYPE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 7: PATTERN UP-TYPE VS DOWN-TYPE")
print("-" * 70)

print("""
Observatie: Up-type si down-type au sarcini diferite (+2/3 vs -1/3).
Poate exponentii reflecta aceasta diferenta?
""")

print("Diferenta n(up-type) - n(down-type) per generatie:")
for gen in generations:
    u_name, u_n = gen['up_type']
    d_name, d_n = gen['down_type']
    diff = u_n - d_n
    print(f"  Gen {generations.index(gen)+1}: {u_n} - {d_n} = {diff}")

print("""
Diferentele: -2, 5, 7
Nu exista un pattern clar... sau poate:
  -2 = ? (negativ!)
  5 = diametru
  7 = dim(Im(O)) octonioni
""")

# ============================================================
# PASUL 8: RAPOARTE INTRA-GENERATIE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 8: RAPOARTE INTRA-GENERATIE")
print("-" * 70)

print("Rapoarte de masa in cadrul fiecarei generatii:")
print()

for i, gen in enumerate(generations):
    l_name, l_n = gen['lepton']
    u_name, u_n = gen['up_type']
    d_name, d_n = gen['down_type']

    m_l = leptons[l_name]['mass'] if l_name in leptons else m_e * phi**l_n
    m_u = quarks[u_name]['mass']
    m_d = quarks[d_name]['mass']

    print(f"Generatia {i+1}:")
    print(f"  m_{u_name}/m_{l_name} = {m_u/m_l:.2f} ~ phi^{u_n - l_n} = phi^{u_n - l_n} = {phi**(u_n - l_n):.2f}")
    print(f"  m_{d_name}/m_{l_name} = {m_d/m_l:.2f} ~ phi^{d_n - l_n} = phi^{d_n - l_n} = {phi**(d_n - l_n):.2f}")
    print(f"  m_{u_name}/m_{d_name} = {m_u/m_d:.2f} ~ phi^{u_n - d_n} = phi^{u_n - d_n} = {phi**(u_n - d_n):.2f}")
    print()

# ============================================================
# PASUL 9: CONEXIUNEA CU GEOMETRIA 600-CELL
# ============================================================
print("\n" + "-" * 70)
print("PASUL 9: CONEXIUNEA CU GEOMETRIA 600-CELL")
print("-" * 70)

print("""
NUMERELE CHEIE din 600-cell:
  120 = varfuri
  720 = muchii
  1200 = fete (triunghiuri)
  600 = celule (tetraedre)
  12 = grad (vecini/vertex)
  5 = diametru
  6 = decagoane/vertex
  20 = tetraedre/vertex
  8 = varfuri 16-cell
  16 = varfuri tesseract
  96 = varfuri fermionice

EXPONENTII QUARKS SI SEMNIFICATIA LOR:
""")

exponent_meanings = {
    3: "dim(Im(H)) = 3 (quaternioni) sau 3 culori",
    5: "Diametru graf 600-cell",
    11: "5 + 6 (diametru + decagoane) = MUON!",
    16: "Varfuri tesseract (8-cell)",
    19: "16 + 3? sau 11 + 8?",
    26: "2 * 13? sau 20 + 6?",
}

print(f"{'n':<6} {'Quark':<10} {'Interpretare posibila':<40}")
print("-" * 60)
for name, data in quark_exponents.items():
    n = data['n_round']
    meaning = exponent_meanings.get(n, "?")
    print(f"{n:<6} {name:<10} {meaning:<40}")

# ============================================================
# PASUL 10: FORMULA UNIFICATA PROPUSA
# ============================================================
print("\n" + "-" * 70)
print("PASUL 10: FORMULA UNIFICATA PROPUSA")
print("-" * 70)

print("""
IPOTEZA: Toti fermionii urmeaza formula

  m = m_e * phi^n

Unde n e determinat de:
  1. TIPUL particulei (lepton vs quark)
  2. GENERATIA (1, 2, 3)
  3. SARCINA (up-type vs down-type pentru quarks)

STRUCTURA PROPUSA:

LEPTONI:
  n_l(gen) = 0, 11, 17  (diferente: 11, 6)

DOWN-TYPE QUARKS:
  n_d(gen) = 5, 11, 19  (diferente: 6, 8)

UP-TYPE QUARKS:
  n_u(gen) = 3, 16, 26  (diferente: 13, 10)
""")

# Verificare
print("\nVERIFICARE NUMERICA:")
print(f"{'Particula':<12} {'n propus':<10} {'m_calc':<12} {'m_exp':<12} {'Err %':<10}")
print("-" * 60)

proposed = [
    ('electron', 0, 0.511),
    ('muon', 11, 105.66),
    ('tau', 17, 1776.8),
    ('up', 3, 2.16),
    ('down', 5, 4.67),
    ('strange', 11, 93.4),
    ('charm', 16, 1270),
    ('bottom', 19, 4180),
    ('top', 26, 172690),
]

for name, n, m_exp in proposed:
    m_calc = m_e * phi**n
    err = abs(m_calc - m_exp) / m_exp * 100
    print(f"{name:<12} {n:<10} {m_calc:<12.2f} {m_exp:<12.2f} {err:<10.1f}")

# ============================================================
# PASUL 11: OBSERVATIA CHEIE - STRANGE = MUON
# ============================================================
print("\n" + "-" * 70)
print("PASUL 11: OBSERVATIA CHEIE - STRANGE = MUON")
print("-" * 70)

print("""
CEA MAI IMPORTANTA OBSERVATIE:
==============================
Strange quark si muon au ACELASI exponent n = 11!

Aceasta sugereaza ca:
1. Ambele sunt "prima excitatie majora" deasupra bazei
2. Exponentul 11 = 5 + 6 e UNIVERSAL pentru a doua generatie
3. Diferenta de masa (93.4 vs 105.66 MeV) vine din alte efecte
   (QCD running, culoare, etc.)

RAPORTUL:
  m_strange / m_muon = 93.4 / 105.66 = 0.884

Aceasta diferenta de ~12% ar putea veni din:
  - Corectii QCD (running mass)
  - Factorul de culoare (medie pe 3 stari)
  - Scala diferita de renormalizare
""")

# ============================================================
# PASUL 12: PATTERN-URI EMERGENTE
# ============================================================
print("\n" + "-" * 70)
print("PASUL 12: PATTERN-URI EMERGENTE")
print("-" * 70)

print("""
DIFERENTE INTRE GENERATII CONSECUTIVE:

Leptoni:     11, 6       (total 17)
Down-type:   6, 8        (total 14)
Up-type:     13, 10      (total 23)

OBSERVATII:
1. Leptonii au cea mai simpla structura (5+6, apoi 6)
2. Down-type: 6 (decagoane), 8 (16-cell)
3. Up-type: 13 (?), 10 (decagon complet)

SUMA DIFERENTELOR:
  Leptoni: 11 + 6 = 17
  Down: 6 + 8 = 14
  Up: 13 + 10 = 23

  Total: 17 + 14 + 23 = 54 = ?

VERIFICARE: 54 = 2 * 27 = 2 * 3^3
            sau 54 = 6 * 9 = decagoane * ?
""")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZII EXP-093")
print("=" * 70)

print("""
REZULTATE PRINCIPALE:
=====================

1. STRANGE = MUON (n = 11):
   - Cea mai importanta descoperire
   - Acelasi exponent pentru a doua generatie
   - Diferenta de masa (~12%) vine din efecte secundare

2. DOWN = DIAMETRU (n = 5):
   - Down quark are exponentul egal cu diametrul 600-cell
   - Sugereaza ca down e "excitarea minima"

3. UP = QUATERNIONI (n = 3):
   - Up quark are n = 3 = dim(Im(H))
   - Posibil legatura cu structura algebrica

4. STRUCTURA PE GENERATII:
   Gen 1: e(0), u(3), d(5)   - baza + cuaternioni + diametru
   Gen 2: mu(11), c(16), s(11) - strange=muon, charm=tesseract
   Gen 3: tau(17), t(26), b(19) - extensii ale gen 2

5. PATTERN DIFERENTE:
   - Leptoni: +11, +6
   - Down-type: +6, +8
   - Up-type: +13, +10

   Numerele 6, 8, 10 apar frecvent (decagoane, 16-cell, decagon)

STATUS:
=======
PATTERN PARTIAL - exponentii au sens geometric partial
DERIVARE INCOMPLETA - nu intelegem inca diferentele up vs down

URMATORUL PAS:
==============
Investigarea rolului CULORII (SU(3)) in modificarea exponentilor.
Posibil: n_quark = n_lepton + f(culoare)
""")

# ============================================================
# BONUS: GRAFIC CONCEPTUAL
# ============================================================
print("\n" + "-" * 70)
print("BONUS: STRUCTURA IERARHICA")
print("-" * 70)

print("""
Exponent n:

  26 |                                              top
  25 |
  24 |
  23 |
  22 |
  21 |
  20 |
  19 |                                    bottom
  18 |
  17 |                          tau
  16 |                                    charm
  15 |
  14 |
  13 |
  12 |
  11 |              muon        strange
  10 |
   9 |
   8 |
   7 |
   6 |
   5 |                          down
   4 |
   3 |                          up
   2 |
   1 |
   0 |  electron
     +------------------------------------------
        Lepton      Down-type   Up-type

OBSERVATIE: Strange si muon sunt pe aceeasi linie orizontala!
""")
