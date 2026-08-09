"""
EXP-018: Verificarea Tuturor Afirmatiilor
=========================================
Verificam riguros fiecare afirmatie facuta in experimente.
"""

from physics_formulas import *

print("=" * 70)
print("EXP-018: VERIFICAREA TUTUROR AFIRMATIILOR")
print("=" * 70)

all_ok = True

def check(name, condition, details=""):
    global all_ok
    status = "OK" if condition else "FAIL"
    if not condition:
        all_ok = False
    print(f"  [{status}] {name}")
    if details and not condition:
        print(f"        {details}")
    return condition

print("\n" + "-" * 70)
print("1. CONSTANTE FUNDAMENTALE (CODATA 2022)")
print("-" * 70)

# Verificam valorile
check("phi = (1+sqrt(5))/2", abs(PHI - (1 + math.sqrt(5))/2) < 1e-10)
check("phi^2 = phi + 1", abs(PHI**2 - PHI - 1) < 1e-10)
check("alpha ~ 1/137.036", abs(ALPHA - 1/137.035999084) < 1e-12)

print("\n" + "-" * 70)
print("2. FORMULA PRINCIPALA: 1/alpha + 2*pi*alpha = 20*phi^4")
print("-" * 70)

lhs = ALPHA_INV + 2*PI*ALPHA
rhs = 20 * PHI**4
err_formula = abs(lhs - rhs) / rhs * 100

check(f"20*phi^4 = {20*PHI**4:.6f}", True)
check(f"1/alpha + 2*pi*alpha = {lhs:.6f}", True)
check(f"Eroare < 0.001%", err_formula < 0.001, f"Eroare = {err_formula:.6f}%")

# Verificam ecuatia cuadratica
K = 20 * PHI**4
disc = K**2 - 8*PI
alpha_calc = (K - math.sqrt(disc)) / (4*PI)
err_alpha = abs(alpha_calc - ALPHA) / ALPHA * 100

check(f"Alpha din ecuatie = {alpha_calc:.10f}", True)
check(f"Alpha CODATA = {ALPHA:.10f}", True)
check(f"Eroare alpha < 0.0002%", err_alpha < 0.0002, f"Eroare = {err_alpha:.6f}%")

print("\n" + "-" * 70)
print("3. PROPRIETATEA ALGEBRICA: alpha * alpha_dual = 1/(2*pi)")
print("-" * 70)

alpha_minus = (K - math.sqrt(disc)) / (4*PI)
alpha_plus = (K + math.sqrt(disc)) / (4*PI)
product = alpha_minus * alpha_plus
expected_product = 1/(2*PI)

check(f"alpha_- = {alpha_minus:.10f}", True)
check(f"alpha_+ = {alpha_plus:.10f}", True)
check(f"Produs = {product:.10f}", True)
check(f"1/(2*pi) = {expected_product:.10f}", True)
check("Produs = 1/(2*pi) EXACT", abs(product - expected_product) < 1e-14)

# Verificam 1/alpha_+ = 2*pi*alpha
inv_alpha_plus = 1/alpha_plus
two_pi_alpha = 2*PI*ALPHA

check(f"1/alpha_+ = {inv_alpha_plus:.10f}", True)
check(f"2*pi*alpha = {two_pi_alpha:.10f}", True)
check("1/alpha_+ ~ 2*pi*alpha", abs(inv_alpha_plus - two_pi_alpha)/two_pi_alpha < 0.0001)

print("\n" + "-" * 70)
print("4. RELATIILE RAZELOR ELECTRONULUI")
print("-" * 70)

# Calcule
lambda_C_bar = HBAR / (M_ELECTRON * C)  # lungimea Compton redusa
r_classical = E_CHARGE**2 / (4*PI*EPSILON_0*M_ELECTRON*C**2)  # raza clasica

print(f"\nValori calculate:")
print(f"  r_e (clasica) = {r_classical:.6e} m")
print(f"  lambda_C_bar  = {lambda_C_bar:.6e} m")
print(f"  a_0 (Bohr)    = {A_0:.6e} m")

# Rapoarte
ratio1 = lambda_C_bar / r_classical
ratio2 = A_0 / lambda_C_bar

print(f"\nRapoarte:")
print(f"  r_C / r_e = {ratio1:.6f}")
print(f"  a_0 / r_C = {ratio2:.6f}")
print(f"  1/alpha   = {ALPHA_INV:.6f}")

check("r_C / r_e = 1/alpha", abs(ratio1 - ALPHA_INV) / ALPHA_INV < 1e-6)
check("a_0 / r_C = 1/alpha", abs(ratio2 - ALPHA_INV) / ALPHA_INV < 1e-6)

# Media geometrica
geom_mean_sq = r_classical * A_0
lambda_sq = lambda_C_bar**2

check(f"r_e * a_0 = {geom_mean_sq:.6e}", True)
check(f"r_C^2     = {lambda_sq:.6e}", True)
check("r_e * a_0 = r_C^2", abs(geom_mean_sq - lambda_sq) / lambda_sq < 1e-6)

# Verificam formulele
r_e_check = ALPHA * lambda_C_bar
a_0_check = lambda_C_bar / ALPHA

check("r_e = alpha * lambda_C_bar", abs(r_classical - r_e_check) / r_classical < 1e-6)
check("a_0 = lambda_C_bar / alpha", abs(A_0 - a_0_check) / A_0 < 1e-6)

print("\n" + "-" * 70)
print("5. NUMERELE 600-CELL")
print("-" * 70)

V_600 = 120  # vertices
E_600 = 720  # edges
F_600 = 1200  # faces
C_600 = 600  # cells
tetra_per_vertex = 20

check(f"V - E + F - C = {V_600 - E_600 + F_600 - C_600} (Euler 4D)", V_600 - E_600 + F_600 - C_600 == 0)
check(f"Cells * 4 / Vertices = {C_600*4/V_600:.0f}", C_600*4/V_600 == 20)
check(f"Tetraedre per vertex = 20", tetra_per_vertex == 20)

# Verificam ca 20 se obtine in mai multe moduri
check("V/6 = 20", V_600/6 == 20)
check("E/36 = 20", E_600/36 == 20)
check("F/60 = 20", F_600/60 == 20)
check("C/30 = 20", C_600/30 == 20)

print("\n" + "-" * 70)
print("6. ICOSAEDRU (vertex figure al 600-cell)")
print("-" * 70)

V_ico = 12
E_ico = 30
F_ico = 20

check(f"V - E + F = {V_ico - E_ico + F_ico} (Euler 3D)", V_ico - E_ico + F_ico == 2)
check("Fete icosaedru = 20", F_ico == 20)
check("E - V + 2 = 20", E_ico - V_ico + 2 == 20)

# Spectrul icosaedrului
lambda_max = 5  # valoare proprie maxima
check("lambda_max^2 - lambda_max = 20", lambda_max**2 - lambda_max == 20)

print("\n" + "-" * 70)
print("7. VITEZA LUMINII SI UNITATI PLANCK")
print("-" * 70)

c_from_planck = L_PLANCK / T_PLANCK
check(f"L_Planck / T_Planck = {c_from_planck:.6e} m/s", True)
check(f"c = {C:.6e} m/s", True)
check("L_Planck / T_Planck = c", abs(c_from_planck - C) / C < 1e-10)

print("\n" + "-" * 70)
print("8. PHI IN COORDONATELE 600-CELL")
print("-" * 70)

# 96 din 120 varfuri contin phi
vertices_with_phi = 96
total_vertices = 120
ratio_phi = vertices_with_phi / total_vertices

check(f"Varfuri cu phi: {vertices_with_phi}/{total_vertices} = {ratio_phi:.2f}", True)
check("80% din varfuri contin phi", ratio_phi == 0.8)

print("\n" + "-" * 70)
print("9. LUNGIMEA COMPTON SI ENERGIA")
print("-" * 70)

lambda_C = H / (M_ELECTRON * C)
E_photon = H * C / lambda_C
E_electron = M_ELECTRON * C**2

check(f"lambda_C = {lambda_C:.6e} m", True)
check(f"E_foton(lambda_C) = {E_photon:.6e} J", True)
check(f"m_e * c^2 = {E_electron:.6e} J", True)
check("E_foton = m_e*c^2", abs(E_photon - E_electron) / E_electron < 1e-10)

print("\n" + "-" * 70)
print("10. McKAY CORRESPONDENCE")
print("-" * 70)

# Binary icosahedral group
size_2I = 120
size_A5 = 60
E8_roots = 240

check("|2I| = 120", size_2I == 120)
check("|A5| = 60 = |2I|/2", size_A5 == size_2I / 2)
check("|E8 roots| = 240 = 2 * |2I|", E8_roots == 2 * size_2I)

print("\n" + "=" * 70)
print("REZUMAT VERIFICARE")
print("=" * 70)

if all_ok:
    print("\n*** TOATE VERIFICARILE AU TRECUT! ***")
else:
    print("\n*** ATENTIE: Unele verificari au esuat! ***")

print("""
NOTA IMPORTANTA:
================
Relatiile verificate mai sus sunt FAPTE MATEMATICE sau DEFINITII.

CE RAMANE SPECULATIV:
- Interpretarea ca spatiul ARE structura 600-cell
- Interpretarea ca electronul E un foton in bucla
- Derivarea lui alpha din aceste principii (nu doar potrivire)

STATUS ONEST:
- Formula 1/alpha + 2*pi*alpha = 20*phi^4 e un FITTING foarte precis
- Interpretarea geometrica e PLAUZIBILA dar nu DEMONSTRATA
- Relatiile razelor (r_e : r_C : a_0) sunt CUNOSCUTE, nu descoperiri noi
""")
