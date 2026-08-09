"""
EXP-003: Investigarea relatiei diferenta ~ alpha * 2*pi
=======================================================
Din EXP-002 am observat ca diferenta / alpha / (2*pi) ~ 1.004
"""

from physics_formulas import *

print("=" * 60)
print("EXP-003: DIFERENTA ~ ALPHA * 2*PI ?")
print("=" * 60)

# Valorile
val_icosaedru = 20 * PHI**4
val_alpha_inv = ALPHA_INV
diferenta = val_icosaedru - val_alpha_inv

print(f"\n20*phi^4 = {val_icosaedru:.10f}")
print(f"1/alpha  = {val_alpha_inv:.10f}")
print(f"diferenta = {diferenta:.10f}")

print("\n" + "-" * 60)
print("TEST: diferenta vs alpha * 2*pi")
print("-" * 60)

alpha_2pi = ALPHA * 2 * PI
print(f"\nalpha * 2*pi = {alpha_2pi:.10f}")
print(f"diferenta    = {diferenta:.10f}")
print(f"Raport       = {diferenta / alpha_2pi:.10f}")
print(f"Eroare       = {(diferenta - alpha_2pi) / alpha_2pi * 100:.4f}%")

print("\n" + "-" * 60)
print("Daca ar fi exact, am avea:")
print("-" * 60)

# 1/alpha = 20*phi^4 - alpha*2*pi
# Deci: 1/alpha + alpha*2*pi = 20*phi^4 ?

val_test = val_alpha_inv + alpha_2pi
print(f"\n1/alpha + alpha*2*pi = {val_test:.10f}")
print(f"20*phi^4             = {val_icosaedru:.10f}")
print(f"Diferenta            = {val_icosaedru - val_test:.10f}")
print(f"Eroare               = {(val_icosaedru - val_test) / val_icosaedru * 100:.6f}%")

print("\n" + "-" * 60)
print("Reformulare algebrica")
print("-" * 60)

# Daca 1/alpha = 20*phi^4 - 2*pi*alpha, atunci:
# 1/alpha + 2*pi*alpha = 20*phi^4
# Inmultim cu alpha:
# 1 + 2*pi*alpha^2 = 20*phi^4 * alpha
#
# Sau, din 1/alpha = 20*phi^4 - 2*pi*alpha:
# 1 = alpha * (20*phi^4 - 2*pi*alpha)
# 1 = 20*phi^4*alpha - 2*pi*alpha^2

print("\nDaca presupunem: 1/alpha = 20*phi^4 - 2*pi*alpha")
print("Atunci: 1 = 20*phi^4*alpha - 2*pi*alpha^2")
print()

lhs = 1
rhs = 20 * PHI**4 * ALPHA - 2 * PI * ALPHA**2
print(f"LHS (1)                        = {lhs}")
print(f"RHS (20*phi^4*alpha - 2pi*alpha^2) = {rhs:.10f}")
print(f"Eroare = {abs(lhs - rhs) / lhs * 100:.4f}%")

print("\n" + "-" * 60)
print("Alta abordare: ecuatie pentru alpha")
print("-" * 60)

# Daca 1/alpha + 2*pi*alpha = 20*phi^4
# Fie x = alpha, K = 20*phi^4
# 1/x + 2*pi*x = K
# 1 + 2*pi*x^2 = K*x
# 2*pi*x^2 - K*x + 1 = 0
# x = (K ± sqrt(K^2 - 8*pi)) / (4*pi)

K = 20 * PHI**4
discriminant = K**2 - 8 * PI
print(f"\nK = 20*phi^4 = {K:.10f}")
print(f"Discriminant = K^2 - 8*pi = {discriminant:.10f}")
print(f"sqrt(discriminant) = {math.sqrt(discriminant):.10f}")

alpha_plus = (K + math.sqrt(discriminant)) / (4 * PI)
alpha_minus = (K - math.sqrt(discriminant)) / (4 * PI)

print(f"\nSolutii ale ecuatiei 2*pi*x^2 - K*x + 1 = 0:")
print(f"alpha_+ = {alpha_plus:.10f}  =>  1/alpha_+ = {1/alpha_plus:.6f}")
print(f"alpha_- = {alpha_minus:.10f}  =>  1/alpha_- = {1/alpha_minus:.6f}")

print(f"\nalpha real = {ALPHA:.10f}")
print(f"Eroare alpha_- vs real = {(alpha_minus - ALPHA) / ALPHA * 100:.4f}%")

print("\n" + "-" * 60)
print("CONCLUZIE")
print("-" * 60)

print(f"""
Relatia 1/alpha = 20*phi^4 - 2*pi*alpha NU e exacta.
- Eroare: ~0.4%

Dar daca am PRESUPUNE aceasta relatie, am obtine alpha cu eroare de {(alpha_minus - ALPHA) / ALPHA * 100:.4f}%.

Aceasta nu e o derivare - e un fitting cu forma presupusa!
Fara o justificare FIZICA pentru de ce ar aparea 2*pi*alpha ca corectie,
ramane doar o coincidenta numerica.

STATUS: COINCIDENTA, nu derivare
""")
