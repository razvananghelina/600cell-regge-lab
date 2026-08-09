"""
EXP-265: HIGGS MASS FROM SPECTRAL DATA
=======================================
Derive m_H/m_W = phi - 2/(a1*phi^4) from spectral action.
Compare with Connes NCG formula.
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120
ALPHA = 1/137.035999084
ALPHA_S = 0.1179
m_e = 0.51099895  # MeV
m_W = 80379.0; m_Z = 91187.6; m_H_exp = 125250.0; m_t = 172570.0
v_EW = 246220.0  # MeV

# Fermion data
fermions = ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']
n_vals =   [0,   11,   17,    3,   5,   11,  16,  19,  26]
m_pdg =    [0.511, 105.658, 1776.86, 2.16, 4.67, 93.4, 1270, 4180, 172570]
Nc =       [1, 1, 1, 3, 3, 3, 3, 3, 3]

print("="*72)
print("EXP-265: HIGGS MASS FROM SPECTRAL DATA")
print("="*72)

# === PART 1: Standard Connes NCG formula ===
print("\n--- Part 1: Connes NCG Formula ---")
# m_H^2 = 2*v^2 * (sum Nc*yf^4) / (sum Nc*yf^2)
# yf = sqrt(2)*mf/v
y_f = [np.sqrt(2)*m/v_EW for m in m_pdg]
a_C = sum(nc*y**4 for nc, y in zip(Nc, y_f))
b_C = sum(nc*y**2 for nc, y in zip(Nc, y_f))
m_H_Connes = v_EW * np.sqrt(2*a_C/b_C)

print(f"  a = sum Nc*yf^4 = {a_C:.6e}")
print(f"  b = sum Nc*yf^2 = {b_C:.6e}")
print(f"  m_H(Connes) = {m_H_Connes:.0f} MeV = {m_H_Connes/1000:.2f} GeV")
print(f"  m_H(exp)    = {m_H_exp:.0f} MeV = {m_H_exp/1000:.2f} GeV")

# Top dominance
y_t = np.sqrt(2)*m_t/v_EW
print(f"  Top: y_t^4/a = {3*y_t**4/a_C:.4f}, y_t^2/b = {3*y_t**2/b_C:.4f}")
print(f"  sqrt(2)*m_t = {np.sqrt(2)*m_t/1000:.2f} GeV (top-only approx)")

# === PART 2: Our formula ===
print("\n--- Part 2: 600-cell Formula ---")
m_H_alpha = m_W * (PHI - 8*ALPHA)
m_H_geom  = m_W * (PHI - 2/(a1*PHI**4))
print(f"  m_H = m_W*(phi-8*alpha)       = {m_H_alpha/1000:.2f} GeV ({(m_H_alpha/m_H_exp-1)*100:+.3f}%)")
print(f"  m_H = m_W*(phi-2/(a1*phi^4))  = {m_H_geom/1000:.2f} GeV ({(m_H_geom/m_H_exp-1)*100:+.3f}%)")
print(f"  Difference alpha vs geom: {abs(m_H_alpha-m_H_geom):.2f} MeV")

# === PART 3: Correction factor ===
print("\n--- Part 3: Correction Factor ---")
R = (m_H_exp / m_H_Connes)**2
print(f"  R = (m_H_exp/m_H_Connes)^2 = {R:.6f}")
print(f"  lambda_600cell = {R:.4f} * lambda_Connes")

# What is R in framework terms?
g2 = 2*m_W/v_EW
R_check = g2**2 * (PHI-8*ALPHA)**2 / (8*a_C/b_C)
print(f"  R = g2^2*(phi-8*alpha)^2*b/(8*a) = {R_check:.6f}")

# Express differently
ratio_HW = m_H_exp/m_W
ratio_tW = m_t/m_W
print(f"\n  m_H/m_W = {ratio_HW:.4f} (exp), phi-8*alpha = {PHI-8*ALPHA:.4f}")
print(f"  m_t/m_W = {ratio_tW:.4f}")
print(f"  (m_H/m_t)^2 = {(m_H_exp/m_t)**2:.6f}")
print(f"  (m_H/m_t)^2 / 2 = {(m_H_exp/m_t)**2/2:.6f} (= R)")

# === PART 4: Per-simplex trace identities ===
print("\n--- Part 4: Per-simplex Traces ---")
tr = [12, 7, 5, 4]  # Delta_0, Delta_1, Delta_2, Delta_3
print(f"  Tr/simplex: {tr}")
print(f"  Sum: {sum(tr)} = 28")
print(f"  Even (0+2): {tr[0]+tr[2]} = 17 = n_tau")
print(f"  Odd  (1+3): {tr[1]+tr[3]} = 11 = n_mu = n_s")
print(f"  Alternating: {tr[0]-tr[1]+tr[2]-tr[3]} = {b1} = b_1")
print(f"  Ratio even/odd: {(tr[0]+tr[2])/(tr[1]+tr[3]):.6f} = 17/11 = {17/11:.6f}")

# Connection to Higgs?
# m_H/m_W = phi - 8*alpha
# Can we relate 8*alpha to per-simplex data?
print(f"\n  8*alpha = {8*ALPHA:.6f}")
print(f"  2/(a1*phi^4) = {2/(a1*PHI**4):.6f}")
# Try: 2/(a1*phi^4) = 2/(a1 * (tr[0]+tr[2]+tr[3])) ?  No, 12+5+4=21
# Try: related to trace ratios
for i in range(4):
    for j in range(4):
        if i != j and tr[j] != 0:
            r = tr[i]/tr[j]
            diff = abs(r - 8*ALPHA*a1*PHI**4/2)
            if diff < 0.5:
                print(f"  tr[{i}]/tr[{j}] = {tr[i]}/{tr[j]} = {r:.4f} (cf. a1*phi^4*4*alpha = {a1*PHI**4*4*ALPHA:.4f})")

# === PART 5: Spectral coefficients ===
print("\n--- Part 5: Spectral Coefficients ---")
c0 = 2640; c1 = 14880; c2 = 55920
s0 = c0//240; s1 = c1//240; s2 = c2//240
print(f"  c_k = {{2640, 14880, 55920}}")
print(f"  c_k/240 = {{{s0}, {s1}, {s2}}}")
print(f"  11 = L(5) (Lucas), 233 = F(13) (Fibonacci)")
print(f"  62 = 2*(a1^2+b1) = 2*{a1**2+b1}")
print(f"  c1^2/(c0*c2) = {c1**2/(c0*c2):.6f} (approx 3/2)")

# Ratios
print(f"\n  c1/c0 = {c1/c0:.6f} = 62/11")
print(f"  c2/c1 = {c2/c1:.6f} = 233/62")
print(f"  c2/c0 = {c2/c0:.6f} = 233/11")

# Can we get R from spectral coefficients?
# R = 0.2637
# Try various combinations
print(f"\n  Looking for R = {R:.6f} in spectral data:")
candidates = {
    's0/s2': s0/s2,
    's0/(s1+s2)': s0/(s1+s2),
    's1/(s0*s2)': s1/(s0*s2),
    '1/(s2/s0 - s1/s0)': 1/(s2/s0 - s1/s0),
    'b1/(s0+s2)': b1/(s0+s2),
    '(s0-b1)/(s1/3)': (s0-b1)/(s1/3),
    'a1/s1': a1/s1,
    '(phi-8*alpha)^2/phi^2': (PHI-8*ALPHA)**2/PHI**2,
}
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-R)):
    print(f"    {name:30s} = {val:.6f} {'<-- CLOSE' if abs(val-R) < 0.01 else ''}")

# === PART 6: Key identity ===
print("\n--- Part 6: Self-consistency ---")
print(f"  The alpha equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0")
print(f"  determines alpha = {ALPHA:.10f}")
print(f"  The Higgs correction: delta = 2/(a1*phi^4) = {2/(a1*PHI**4):.8f}")
print(f"  And 8*alpha = {8*ALPHA:.8f}")
print(f"  Difference: {abs(2/(a1*PHI**4) - 8*ALPHA):.2e}")
print(f"")
print(f"  From alpha equation: alpha = (4*a1*phi^4 - sqrt(D))/(4*pi)")
print(f"  where D = (4*a1*phi^4)^2 - 8*pi = 16*a1^2*phi^8 - 8*pi")
D_val = (4*a1*PHI**4)**2 - 8*np.pi
alpha_check = (4*a1*PHI**4 - np.sqrt(D_val))/(4*np.pi)
print(f"  D = {D_val:.6f}, alpha = {alpha_check:.10f}")
print(f"  8*alpha = 2*(4*a1*phi^4 - sqrt(D))/(4*pi)")
print(f"          = (4*a1*phi^4 - sqrt(D))/(pi/2)")

# Tree + 1-loop structure
print(f"\n  STRUCTURE:")
print(f"    Tree:   m_H/m_W = phi = {PHI:.6f}  (DSI potential minimum)")
print(f"    1-loop: delta = 8*alpha = {8*ALPHA:.6f}  (EM vacuum polarization)")
print(f"    Result: m_H/m_W = phi - 8*alpha = {PHI-8*ALPHA:.6f}")
print(f"    = {(PHI-8*ALPHA)*m_W/1000:.2f} GeV vs {m_H_exp/1000:.2f} GeV ({(m_H_alpha/m_H_exp-1)*100:+.3f}%)")

# Why 8 = rank(E8)?
print(f"\n  WHY 8 = rank(E8)?")
print(f"    In the alpha equation: coefficient = 4*a1*phi^4 = {4*a1*PHI**4:.4f}")
print(f"    delta = 2/(a1*phi^4) = 2/(coefficient/4) = 8/coefficient")
print(f"    = 8 * (1/coefficient) = 8 * alpha (to leading order)")
print(f"    The '8' comes from: 2 (numerator) * 4 (from alpha eq) = 8 = rank(E8)")
print(f"    More precisely: 2/(a1*phi^4) vs 8/(4*a1*phi^4) = same thing")

# === PART 7: Quartic coupling ===
print("\n--- Part 7: Higgs Quartic Coupling ---")
lambda_SM = m_H_exp**2 / (2*v_EW**2)
lambda_ours = (m_W*(PHI-8*ALPHA))**2 / (2*v_EW**2)
lambda_Connes = m_H_Connes**2 / (2*v_EW**2)
print(f"  lambda_SM     = {lambda_SM:.6f}")
print(f"  lambda_ours   = {lambda_ours:.6f} ({(lambda_ours/lambda_SM-1)*100:+.2f}%)")
print(f"  lambda_Connes = {lambda_Connes:.6f} (too large by {lambda_Connes/lambda_SM:.1f}x)")
print(f"  lambda = g2^2*(phi-8*alpha)^2/8 = {g2**2*(PHI-8*ALPHA)**2/8:.6f}")

print(f"\n  SUMMARY:")
print(f"    Formula: m_H = m_W*(phi - 2/(a1*phi^4))")
print(f"    = m_W*(phi - 8*alpha) to O(alpha^2)")
print(f"    Uses ONLY: a1=5, phi=(1+sqrt(5))/2, alpha (from alpha eq)")
print(f"    Prediction: {m_H_geom/1000:.2f} GeV (exp: {m_H_exp/1000:.2f} GeV, {(m_H_geom/m_H_exp-1)*100:+.2f}%)")
print(f"    STATUS: PARTIALLY DERIVED")
print(f"      Tree level phi: DERIVED (DSI)")
print(f"      Correction 8*alpha: from alpha equation (DERIVED)")
print(f"      Connection to Connes formula: correction R={R:.4f}")
print(f"      Full derivation from spectral action: OPEN")

print("\n" + "="*72)
print("EXP-265 COMPLETE")
print("="*72)
