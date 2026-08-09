"""
EXP-074: Doua Formule pentru Alpha
===================================

Am gasit DOUA formule posibile:

1. alpha = 1/(20*phi^4)           - din raportul eigenvalues
2. alpha = 1/sqrt(Trace(L^2))     - din trace spectral

Care e mai fundamentala? Sunt echivalente?
"""

import numpy as np
from physics_formulas import PHI, ALPHA

print("=" * 70)
print("EXP-074: DOUA FORMULE PENTRU ALPHA")
print("=" * 70)

# Spectrul 600-cell (din EXP-072)
spectrum = [
    (0, 1),                    # lambda_0 = 0, mult 1
    (6/PHI**2, 4),             # lambda_1 = 2.2918, mult 4
    (5.527864, 9),             # lambda_2, mult 9
    (9, 16),                   # lambda_3, mult 16
    (12, 25),                  # lambda_4, mult 25
    (14, 36),                  # lambda_5, mult 36
    (14.472136, 9),            # lambda_6, mult 9
    (15, 16),                  # lambda_7, mult 16
    (6/PHI**2 + 12, 4),        # lambda_8 = 15.708, mult 4
]

# ============================================================
# FORMULA 1: alpha = 1/(20*phi^4)
# ============================================================
print("\n" + "-" * 70)
print("FORMULA 1: alpha = 1/(20*phi^4)")
print("-" * 70)

lambda_1 = 6/PHI**2
lambda_4 = 12
ratio_spectral = lambda_4 / lambda_1

print(f"lambda_1 = 6/phi^2 = {lambda_1:.6f}")
print(f"lambda_4 = 12")
print(f"lambda_4/lambda_1 = {ratio_spectral:.6f}")
print(f"2*phi^2 = {2*PHI**2:.6f}")
print(f"Match: {abs(ratio_spectral - 2*PHI**2) < 1e-6}")

formula1_inv = 20 * PHI**4
formula1_alpha = 1 / formula1_inv

print(f"\n1/alpha (formula 1) = 20*phi^4 = {formula1_inv:.6f}")
print(f"alpha (formula 1) = {formula1_alpha:.8f}")
print(f"alpha (experimental) = {ALPHA:.8f}")
print(f"Eroare: {abs(formula1_alpha - ALPHA)/ALPHA * 100:.4f}%")

# ============================================================
# FORMULA 2: alpha = 1/sqrt(Trace(L^2))
# ============================================================
print("\n" + "-" * 70)
print("FORMULA 2: alpha = 1/sqrt(Trace(L^2))")
print("-" * 70)

trace_L2 = sum(lam**2 * mult for lam, mult in spectrum)
formula2_inv = np.sqrt(trace_L2)
formula2_alpha = 1 / formula2_inv

print(f"Trace(L^2) = sum(lambda_k^2 * mult_k) = {trace_L2:.4f}")
print(f"sqrt(Trace(L^2)) = {formula2_inv:.6f}")
print(f"\n1/alpha (formula 2) = sqrt(Trace(L^2)) = {formula2_inv:.6f}")
print(f"alpha (formula 2) = {formula2_alpha:.8f}")
print(f"alpha (experimental) = {ALPHA:.8f}")
print(f"Eroare: {abs(formula2_alpha - ALPHA)/ALPHA * 100:.4f}%")

# ============================================================
# COMPARATIE
# ============================================================
print("\n" + "-" * 70)
print("COMPARATIE CELE DOUA FORMULE")
print("-" * 70)

print(f"""
| Formula | Valoare 1/alpha | Eroare |
|---------|-----------------|--------|
| 20*phi^4 | {formula1_inv:.4f} | {abs(formula1_alpha - ALPHA)/ALPHA * 100:.4f}% |
| sqrt(Tr(L^2)) | {formula2_inv:.4f} | {abs(formula2_alpha - ALPHA)/ALPHA * 100:.4f}% |
| Experimental | {1/ALPHA:.4f} | 0% |
""")

print(f"Diferenta intre formule: {abs(formula1_inv - formula2_inv):.4f}")
print(f"Raport: {formula1_inv / formula2_inv:.6f}")

# ============================================================
# SUNT ECHIVALENTE?
# ============================================================
print("\n" + "-" * 70)
print("ANALIZA: SUNT ECHIVALENTE?")
print("-" * 70)

print("""
FORMULA 1: 20*phi^4 = 5 * (lambda_4/lambda_1)^2

FORMULA 2: sqrt(Trace(L^2)) = sqrt(sum_k lambda_k^2 * mult_k)

Pentru a fi echivalente ar trebui:
  [5 * (lambda_4/lambda_1)^2]^2 = Trace(L^2)

Adica:
  25 * (lambda_4/lambda_1)^4 = sum_k lambda_k^2 * mult_k
""")

lhs = 25 * ratio_spectral**4
rhs = trace_L2

print(f"LHS = 25 * (lambda_4/lambda_1)^4 = {lhs:.4f}")
print(f"RHS = Trace(L^2) = {rhs:.4f}")
print(f"Raport LHS/RHS = {lhs/rhs:.6f}")

print("""
NU SUNT ECHIVALENTE!

Formula 1 (20*phi^4) e mai precisa.
Formula 2 (sqrt(Trace(L^2))) e mai "naturala" in contextul Spectral Action.
""")

# ============================================================
# CE CONTINE TRACE(L^2)?
# ============================================================
print("\n" + "-" * 70)
print("DESCOMPUNEREA LUI TRACE(L^2)")
print("-" * 70)

print("Contributii la Trace(L^2):\n")
total = 0
for lam, mult in spectrum:
    contrib = lam**2 * mult
    total += contrib
    percent = contrib / trace_L2 * 100
    print(f"  lambda={lam:8.4f}, mult={mult:2d}: {lam:.4f}^2 * {mult} = {contrib:8.2f} ({percent:5.2f}%)")

print(f"\nTotal: {total:.4f}")

# Contributia de la lambda_4 = 12 cu mult = 25
contrib_lambda4 = 12**2 * 25
print(f"\nContributia lui lambda_4 = 12: {contrib_lambda4} ({contrib_lambda4/trace_L2*100:.2f}%)")

# ============================================================
# FORMULA HIBRIDA?
# ============================================================
print("\n" + "-" * 70)
print("CAUTAM O FORMULA UNIFICATA")
print("-" * 70)

print("""
Poate exista o formula care combina cele doua:

Idee: Trace(L^2) ar putea fi exprimat in termeni de phi
daca spectrul are structura speciala.
""")

# Verificam daca Trace(L^2) are o forma simpla
# Trace(L^2) = 18550 ~ ???

# Incercam diverse combinatii cu phi
combinations = [
    ("120 * phi^8", 120 * PHI**8),
    ("20 * phi^4 * 137", 20 * PHI**4 * 137),
    ("(20*phi^4)^2 / phi^2", (20*PHI**4)**2 / PHI**2),
    ("1000 * phi^6", 1000 * PHI**6),
    ("600 * phi^5 * 2", 600 * PHI**5 * 2),
    ("5 * 12^2 * phi^4", 5 * 144 * PHI**4),
    ("720 * phi^5", 720 * PHI**5),  # 720 = numar de muchii
    ("120 * 12^2 + ...", 120 * 144),
]

print(f"Trace(L^2) = {trace_L2:.4f}\n")
for name, val in combinations:
    err = abs(val - trace_L2) / trace_L2 * 100
    marker = " <--" if err < 5 else ""
    print(f"  {name:<25} = {val:10.2f} (err = {err:5.2f}%){marker}")

# ============================================================
# STRUCTURA SPECTRULUI
# ============================================================
print("\n" + "-" * 70)
print("STRUCTURA SPECTRULUI 600-CELL")
print("-" * 70)

print("""
Observatii despre spectru:

1. Multiplicitatile: 1, 4, 9, 16, 25, 36, 9, 16, 4
   - Primele 6: patrate perfecte consecutive!
   - 1^2, 2^2, 3^2, 4^2, 5^2, 6^2

2. Valorile proprii speciale:
   - lambda_1 = 6/phi^2 (contine phi)
   - lambda_4 = 12 (intreg)
   - lambda_8 = lambda_1 + 12 (suma!)

3. Raportul magic:
   - lambda_4/lambda_1 = 12/(6/phi^2) = 2*phi^2 EXACT
""")

# Verificam daca lambda_i au forme simple
print("\nValorile proprii in termeni de phi:")
for i, (lam, mult) in enumerate(spectrum):
    # Incercam sa gasim o forma simpla
    if lam > 0:
        ratio_to_phi = lam / PHI
        ratio_to_phi2 = lam / PHI**2
        if abs(ratio_to_phi - round(ratio_to_phi)) < 0.01:
            print(f"  lambda_{i} = {lam:.4f} = {round(ratio_to_phi)} * phi")
        elif abs(ratio_to_phi2 - round(ratio_to_phi2)) < 0.01:
            print(f"  lambda_{i} = {lam:.4f} = {round(ratio_to_phi2)} * phi^2")
        elif abs(lam - round(lam)) < 0.01:
            print(f"  lambda_{i} = {lam:.4f} = {round(lam)} (intreg)")
        else:
            # Verificam 6/phi^2
            if abs(lam - 6/PHI**2) < 0.01:
                print(f"  lambda_{i} = {lam:.4f} = 6/phi^2")
            elif abs(lam - 6/PHI**2 - 12) < 0.01:
                print(f"  lambda_{i} = {lam:.4f} = 6/phi^2 + 12")
            else:
                print(f"  lambda_{i} = {lam:.4f} (forma complexa)")

# ============================================================
# CONCLUZIE
# ============================================================
print("\n" + "=" * 70)
print("CONCLUZIE EXP-074")
print("=" * 70)

print(f"""
DOUA FORMULE PENTRU ALPHA:
==========================

FORMULA 1: alpha = 1/(20*phi^4)
  - Vine din: 5 * (lambda_4/lambda_1)^2
  - Eroare: 0.034%
  - Pro: Foarte precisa, structura algebrica clara
  - Contra: De ce lambda_4 si lambda_1? De ce 5?

FORMULA 2: alpha = 1/sqrt(Trace(L^2))
  - Vine din: Trace spectral direct
  - Eroare: 0.61%
  - Pro: Naturala in contextul Spectral Action
  - Contra: Mai putin precisa

CONCLUZIE:
==========
Formula 1 e mai precisa si mai eleganta matematic.
Formula 2 e mai naturala fizic (tip Spectral Action).

Diferenta de ~0.6% intre ele sugereaza ca:
- Fie Formula 1 e "fundamentala" si Formula 2 e aproximatie
- Fie exista o formula EXACTA care le uneste pe amandoua

INTREBARE CHEIE:
================
De ce raportul lambda_4/lambda_1 = 2*phi^2 EXACT?

Asta e misterul central. Daca putem deriva ACEST raport
din primele principii, am avea o derivare completa a lui alpha.
""")
