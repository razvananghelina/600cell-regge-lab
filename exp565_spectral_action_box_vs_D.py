"""
exp565: Spectral action Tr f(Box/Lambda) vs Tr f(D/Lambda).

Box has a fully algebraic spectrum with 12 distinct masses in Z[phi].
D is the simplicial Dirac operator on 2640-dimensional space.
Compare their Seeley-DeWitt coefficients and look for connections to alpha.

KEY IDEA: The ratio c_k(D)/a_k(Box) at each order k might encode alpha
or other coupling constants.

Also: Box^2 is positive semidefinite (proper heat kernel exists).
Compute zeta_Box(s), det'(Box), and spectral invariants.
"""

import numpy as np
import math
from numpy.linalg import eigvalsh, eigh
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
degree = 12

verts, adj, lap = build_600cell()

# Build Hopf fibration and Box
def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

def find_fibration():
    for i in range(N):
        if abs(verts[i, 0] - PHI/2) < 1e-6:
            g = verts[i]; p = g.copy(); ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
                if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
            if not ok: continue
            used = set(); fibers = []; subg = []
            pp = np.array([1.0,0,0,0])
            for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
            for s in range(N):
                if s in used: continue
                fib = []
                for si in subg:
                    q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                    if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
                if len(fib) == 10: fibers.append(fib)
            if len(fibers) == 12: return fibers
    return None

fibers = find_fibration()
A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0

L_fiber = np.diag(np.sum(A_fiber, axis=1)) - A_fiber
L_cross = np.diag(np.sum(adj - A_fiber, axis=1)) - (adj - A_fiber)
Box = L_cross - a1 * L_fiber

evals_box = np.sort(eigvalsh(Box))
nonzero_box = evals_box[np.abs(evals_box) > 1e-6]
positive_box = nonzero_box[nonzero_box > 0]
negative_box = nonzero_box[nonzero_box < 0]


# =====================================================================
# PART 1: SEELEY-DEWITT COEFFICIENTS COMPARISON
# =====================================================================
print("=" * 70)
print("PART 1: SEELEY-DEWITT COEFFICIENTS")
print("=" * 70)

# Box trace moments (already known)
tr = {}
for n in range(0, 9):
    tr[n] = np.sum(evals_box**n)

# D coefficients from paper (Section 9.5)
# c_0 = dim(H) = 2640, c_1 = Tr(D^2) = 14880, c_2 = Tr(D^4)/2 = 55920
# But actually: c_k in Seeley-DeWitt: S ~ sum c_k Lambda^{d-2k}
# c_0 = Tr(1) = dim, c_1 = Tr(D^2), c_2 = (1/2)Tr(D^4)
c_D = {0: 2640, 1: 14880, 2: 55920}

# Box: a_k = Tr(Box^{2k}) / (2k)! ... no, the standard Seeley-DeWitt is:
# Tr(e^{-t*Box^2}) ~ sum_{k>=0} a_k(Box^2) * t^{k-d/2}
# For the operator Box^2 (positive semidefinite):
evals_box2 = evals_box**2

a_Box = {}
for n in range(7):
    a_Box[n] = np.sum(evals_box2**n)  # = Tr(Box^{2n})

print(f"\n  D operator (simplicial, dim 2640):")
print(f"    c_0 = Tr(1)   = {c_D[0]}")
print(f"    c_1 = Tr(D^2) = {c_D[1]}")
print(f"    c_2 = Tr(D^4)/2 = {c_D[2]}")
print(f"    c_0/240 = {c_D[0]/240} = L_5 (Lucas)")
print(f"    c_1/240 = {c_D[1]/240}")
print(f"    c_2/240 = {c_D[2]/240} = F_13 (Fibonacci)")

print(f"\n  Box operator (vertex, dim 120):")
print(f"    Tr(1)     = {N}")
print(f"    Tr(Box)   = {tr[1]:.0f}")
print(f"    Tr(Box^2) = {tr[2]:.0f}")
print(f"    Tr(Box^3) = {tr[3]:.0f}")
print(f"    Tr(Box^4) = {tr[4]:.0f}")
print(f"    Tr(Box^5) = {tr[5]:.0f}")
print(f"    Tr(Box^6) = {tr[6]:.0f}")

print(f"\n  Box^2 operator (positive semidefinite):")
for n in range(5):
    print(f"    Tr(Box^{2*n}) = {a_Box[n]:.0f}")


# =====================================================================
# PART 2: RATIOS D vs BOX
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: RATIOS c_k(D) / Tr(Box^{2k})")
print("=" * 70)

# c_0(D) / Tr(Box^0) = 2640/120 = 22
ratio_0 = c_D[0] / N
print(f"  c_0(D) / N = {ratio_0:.4f} = {c_D[0]}/{N}")
print(f"    = 22 = dim(simplicial)/dim(vertex)")
print(f"    = 1 + 6 + 10 + 5 (vertex contribution per vertex)")

# c_1(D) / Tr(Box^2) = 14880/7200
ratio_1 = c_D[1] / tr[2]
print(f"\n  c_1(D) / Tr(Box^2) = {ratio_1:.6f}")
print(f"    = {c_D[1]}/{tr[2]:.0f}")
print(f"    = 14880/7200 = 62/30 = 31/15")

# c_2(D) / Tr(Box^4) = 55920/756000
ratio_2 = c_D[2] / tr[4]
print(f"\n  c_2(D) / Tr(Box^4) = {ratio_2:.6f}")
print(f"    = {c_D[2]}/{tr[4]:.0f}")

# c_2(D) / (Tr(Box^4)/2) = 55920/(756000/2) = 55920/378000
ratio_2b = c_D[2] / (tr[4]/2)
print(f"  c_2(D) / (Tr(Box^4)/2) = {ratio_2b:.6f}")

# Check: c_1(D) / (N * degree * a1) = 14880/7200 = 31/15
# And 31 = ? Not obviously related to framework.
# 15 = 3*5 = N_gen * a1

# Try: c_1(D) / (N * degree) = 14880 / (120*12) = 14880/1440 = 31/3
# And 31 = a1^2 + b1 = 25+6 = 31. YES!

print(f"\n  --- Algebraic identities ---")
print(f"  c_1(D) / (N * degree) = {c_D[1]/(N*degree):.4f}")
print(f"    = {c_D[1]}/{N*degree} = 31/3")
print(f"    31 = a1^2 + b1 = {a1**2 + b1}")

# So c_1(D) = N * degree * (a1^2+b1)/3
predicted_c1 = N * degree * (a1**2 + b1) / 3
print(f"  Predicted c_1 = N*deg*(a1^2+b1)/3 = {predicted_c1:.0f}")
print(f"  Actual c_1    = {c_D[1]}")
print(f"  Match: {abs(predicted_c1 - c_D[1]) < 0.01}")

# c_0(D) = N * 22. And 22 = 2*L_5 = 2*11. Also 22 = 2(a1^2-a1+1) = 2*21. No, 21 != 11.
# 22 = 2*11 = 2*L_5 where L_5 is the 5th Lucas number.
# Also: dim(H) = V+E+F+C = 120+720+1200+600 = 2640 = 120*22.
# 22 = 1 + 6 + 10 + 5. These are: 1, degree/2, 10, a1.
# Or: f-vector per vertex = (1, 12/2, 20/2, 10/2)... hmm.

# c_2(D) = 55920. c_2/N = 466. c_2/240 = 233 = F_13.
# 233/30 = 7.767... not clean.


# =====================================================================
# PART 3: SPECTRAL ZETA FUNCTION OF BOX
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: SPECTRAL ZETA FUNCTION OF Box")
print("=" * 70)

# zeta_Box(s) = sum_{mu != 0} |mu|^{-s}
# For integer s:
for s in [1, 2, 3, 4]:
    zeta_s = np.sum(np.abs(nonzero_box)**(-s))
    print(f"  zeta_Box({s}) = {zeta_s:.10f}")

# Separate positive and negative:
print(f"\n  Separated by sign:")
for s in [1, 2]:
    z_pos = np.sum(positive_box**(-s))
    z_neg = np.sum(np.abs(negative_box)**(-s))
    print(f"  zeta_+({s}) = {z_pos:.10f},  zeta_-({s}) = {z_neg:.10f}")
    print(f"    ratio z+/z- = {z_pos/z_neg:.10f}")

# Spectral determinant
log_det = np.sum(np.log(np.abs(nonzero_box)))
log_det_pos = np.sum(np.log(positive_box))
log_det_neg = np.sum(np.log(np.abs(negative_box)))
print(f"\n  log|det'(Box)| = {log_det:.6f}")
print(f"  log(det_+)     = {log_det_pos:.6f}")
print(f"  log|det_-|     = {log_det_neg:.6f}")

# zeta'(0) = -log|det'|
print(f"  zeta'(0) = -log|det'| = {-log_det:.6f}")

# Check if det' has a clean form
det_pos = np.exp(log_det_pos)
det_neg = np.exp(log_det_neg)
print(f"\n  det_+ = {det_pos:.6e}")
print(f"  |det_-| = {det_neg:.6e}")
print(f"  |det'| = det_+ * |det_-| = {det_pos * det_neg:.6e}")


# =====================================================================
# PART 4: EXACT ALGEBRAIC SPECTRUM AND CLOSED-FORM EXPRESSIONS
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: EXACT ALGEBRAIC DETERMINANT")
print("=" * 70)

# The Box eigenvalues are all in Z[phi]:
# mu values with multiplicities (from exp562):
mu_spectrum = [
    (-10, 12), (-(6*PHI), 10), (-(5+np.sqrt(5)), 6), (-(6*PHI-3), 8),
    (-(6*PHI-3), 8), (-(5-np.sqrt(5)), 6), (-(6/PHI-2), 12),
    (0, 9),
    (6/PHI, 10), (10-2*np.sqrt(5), 3), (6*PHI-3, 8), (6*PHI-3, 8),
    (6*PHI+2, 12), (12, 5), (10+2*np.sqrt(5), 3),
]

# Cleaner: the distinct |mu| with exact expressions
exact_masses = [
    (0, 9, "0"),
    (6/PHI - 2, 12, "6/phi - 2 = 4phi-4"),
    (5-np.sqrt(5), 6, "5-sqrt5 = a1-sqrt(a1)"),
    (6/PHI, 10, "6/phi = 6(phi-1)"),
    (10-2*np.sqrt(5), 3, "10-2sqrt5 = L2"),
    (6*PHI-3, 16+16, "6phi-3"),  # appears twice (positive + negative)
    (5+np.sqrt(5), 6, "5+sqrt5 = a1+sqrt(a1)"),
    (6*PHI, 10, "6phi"),
    (10, 12, "2a1"),
    (6*PHI+2, 12, "6phi+2"),
    (12, 5, "degree"),
    (10+2*np.sqrt(5), 3, "10+2sqrt5 = L6"),
]

# Compute det'(Box) exactly from algebraic eigenvalues
# det' = product of all nonzero eigenvalues
# Use the exact pairing from the joint spectrum

# From the full spectrum, the nonzero eigenvalues with signs and multiplicities:
signed_spectrum = [
    (-10, 12), (-(6*PHI), 10), (-(5+np.sqrt(5)), 6),
    (-(6*PHI-3), 16), (-(5-np.sqrt(5)), 6), (-(6/PHI-2), 12),
    (6/PHI, 10), (10-2*np.sqrt(5), 3), (6*PHI-3, 16),
    (6*PHI+2, 12), (12, 5), (10+2*np.sqrt(5), 3),
]

# Wait, I need to be more careful. Let me use the actual computed spectrum.
# Group nonzero eigenvalues by value
from collections import Counter
mu_counter = Counter()
for mu in nonzero_box:
    mu_counter[round(mu, 4)] += 1

print(f"  Nonzero Box eigenvalues (signed, with multiplicities):")
log_det_check = 0
for val, mult in sorted(mu_counter.items()):
    log_det_check += mult * np.log(abs(val))
    print(f"    mu = {val:+10.4f}  mult = {mult:3d}  log|mu|*mult = {mult*np.log(abs(val)):8.4f}")

print(f"\n  Sum = log|det'| = {log_det_check:.6f}")
print(f"  Direct computation: {log_det:.6f}")
print(f"  Match: {abs(log_det_check - log_det) < 0.01}")


# =====================================================================
# PART 5: RATIOS WITH alpha, 2*pi, phi
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: SEARCHING FOR alpha IN SPECTRAL INVARIANTS")
print("=" * 70)

alpha_phys = (4*a1*PHI**4 - np.sqrt((4*a1*PHI**4)**2 - 8*np.pi)) / (4*np.pi)
B_alpha = 4*a1*PHI**4
two_pi = 2*np.pi

print(f"  alpha = {alpha_phys:.10f} = 1/{1/alpha_phys:.4f}")
print(f"  B = 4*a1*phi^4 = {B_alpha:.6f}")
print(f"  2*pi = {two_pi:.6f}")

# Compute various spectral quantities and check against alpha
print(f"\n  --- Spectral quantity / alpha comparisons ---")

# Tr(Box^2) / N / B
r1 = tr[2] / N / B_alpha
print(f"  Tr(Box^2)/(N*B) = 7200/(120*{B_alpha:.2f}) = {r1:.6f}")

# Tr(Box^2) / (N * degree) = 60/12 = 5 = a1
print(f"  Tr(Box^2)/(N*deg) = {tr[2]/(N*degree):.4f} = a1")

# Tr(Box^3) / Tr(Box^2) = N^2/(N*60) = N/60 = 2
print(f"  Tr(Box^3)/Tr(Box^2) = {tr[3]/tr[2]:.4f} = N/60 = 2 = deg/b1")

# Tr(Box^4) / Tr(Box^2)^2 = 756000/7200^2
print(f"  Tr(Box^4)/Tr(Box^2)^2 = {tr[4]/tr[2]**2:.6f}")
print(f"    = 756000/51840000 = 7/480 = 7/(4N)")

# 7/(4N) = 7/480. Is 7 = a1+2 = degree/2+1?
print(f"    7 = a1+2? YES. 7/(4N) = (a1+2)/(4*a1!)")

# Kurtosis-like ratio
print(f"\n  Tr(Box^4)*Tr(Box^0) / Tr(Box^2)^2 = {tr[4]*N/tr[2]**2:.6f}")
print(f"    = 756000*120/7200^2 = {756000*120/7200**2:.4f}")
print(f"    = 7*N/(4N) * ... hmm = 7/4 * (N/Tr2) * N = ...")

# Heat kernel at special t values
print(f"\n  --- Heat kernel K(t) = Tr(exp(-t*Box^2)) ---")
evals_box2_sorted = np.sort(evals_box**2)
for t_name, t_val in [("1/(2*pi)", 1/(2*np.pi)), ("alpha", alpha_phys),
                       ("1/B", 1/B_alpha), ("1/a1", 1/a1),
                       ("1/N", 1/N), ("phi/N", PHI/N)]:
    K = np.sum(np.exp(-t_val * evals_box2_sorted))
    print(f"  K(t={t_name:>10s} = {t_val:.6f}) = {K:.6f}   K/N = {K/N:.6f}")

# Spectral action with cutoff Lambda
print(f"\n  --- Spectral action S(Lambda) = count(|mu| < Lambda) ---")
for L_name, L_val in [("a1", a1), ("b1", b1), ("degree", degree),
                       ("2*a1", 2*a1), ("B/a1", B_alpha/a1)]:
    count = np.sum(np.abs(evals_box) < L_val)
    count_nz = np.sum((np.abs(evals_box) > 1e-6) & (np.abs(evals_box) < L_val))
    print(f"  |mu| < {L_name:>8s} ({L_val:8.3f}): {count:3d} modes ({count_nz:3d} massive)")


# =====================================================================
# PART 6: THE GENERATING FUNCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: CHARACTERISTIC POLYNOMIAL OF Box")
print("=" * 70)

# The characteristic polynomial det(Box - lambda*I) encodes ALL information.
# Since Box has algebraic eigenvalues, this polynomial has integer coefficients
# (or coefficients in Z[phi]).

# Compute characteristic polynomial coefficients via Newton's identities
# p_k = Tr(Box^k) are the power sums
# The elementary symmetric polynomials e_k give det(Box-lambda) = sum (-1)^k e_k lambda^{N-k}

# For an alternative approach: use numpy
coeffs = np.polynomial.polynomial.polyfromroots(evals_box)
# This gives coefficients from low to high degree
# The constant term = product of eigenvalues = det(Box)

print(f"  Characteristic polynomial degree: {len(coeffs)-1}")
print(f"  det(Box) = constant term = {coeffs[0]:.6f}")
print(f"  (Should be 0 since Box has kernel)")

# The polynomial has a factor lambda^9 (from the 9 zero eigenvalues)
# Dividing out: det(Box)/lambda^9 evaluated at lambda=0 = det'(Box)
# = product of nonzero eigenvalues

# The reduced polynomial (after removing lambda^9):
reduced_roots = nonzero_box
reduced_coeffs = np.polynomial.polynomial.polyfromroots(reduced_roots)

print(f"\n  Reduced polynomial (111 nonzero roots):")
print(f"  Constant term (= det') = {reduced_coeffs[0]:.6e}")
print(f"  Top coefficient = {reduced_coeffs[-1]:.6f}")

# Check: does det' have a clean form?
# det' = product of all nonzero mu
# Using exact values:
# product = (-10)^12 * (-6phi)^10 * (-(5+sqrt5))^6 * (-(6phi-3))^16 * ...
# This is very large but might simplify.

# Let's compute log|det'| using exact expressions:
print(f"\n  Exact log|det'| composition:")
terms = [
    (10, 12, "|-10|^12"),
    (6*PHI, 10, "|6phi|^10"),
    (5+np.sqrt(5), 6, "|(5+sqrt5)|^6"),
    (6*PHI-3, 16, "|6phi-3|^16"),
    (5-np.sqrt(5), 6, "|(5-sqrt5)|^6"),
    (6/PHI-2, 12, "|6/phi-2|^12"),
    (6/PHI, 10, "|6/phi|^10"),
    (10-2*np.sqrt(5), 3, "|(10-2sqrt5)|^3"),
    (6*PHI+2, 12, "|(6phi+2)|^12"),
    (12, 5, "|12|^5"),
    (10+2*np.sqrt(5), 3, "|(10+2sqrt5)|^3"),
]

total_log = 0
for val, mult, name in terms:
    contribution = mult * np.log(val)
    total_log += contribution
    print(f"    {name:>20s} : {contribution:8.4f}")

print(f"    {'TOTAL':>20s} : {total_log:8.4f}")
print(f"    Direct:                {log_det:.4f}")
print(f"    Match: {abs(total_log - log_det) < 0.01}")

# Some key products:
# (10-2sqrt5)(10+2sqrt5) = 100-20 = 80
# (5-sqrt5)(5+sqrt5) = 25-5 = 20
# (6phi-3)(6/phi-2) = (6phi-3)(6phi-8) = 36phi^2-48phi-18phi+24 = 36(phi+1)-66phi+24
#   = 36phi+36-66phi+24 = 60-30phi = 30(2-phi)
# Check: (6*1.618-3)*(6/1.618-2) = 6.708*1.708 = 11.454. And 30*(2-1.618) = 30*0.382 = 11.46. YES!

print(f"\n  Galois-paired products:")
print(f"    (10-2sqrt5)(10+2sqrt5) = 80 = 16*a1")
print(f"    (5-sqrt5)(5+sqrt5) = 20 = 4*a1")
print(f"    (6phi-3)(6/phi-2) = 30*(2-phi) = 30/phi^2 = {30/PHI**2:.6f}")
print(f"    Check: {(6*PHI-3)*(6/PHI-2):.6f}")
print(f"    6phi * 6/phi = 36")


# =====================================================================
# PART 7: THE RATIO c_1(D) / Tr(Box^2) AND alpha
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: D/Box RATIOS AND COUPLING CONSTANTS")
print("=" * 70)

# c_1(D)/Tr(Box^2) = 14880/7200 = 31/15
# c_2(D)/Tr(Box^4) = 55920/756000 = 233/3150
# Simplify: 233/3150 = 233/(3150). gcd(233,3150)? 233 is prime. So 233/3150.

print(f"  Ratio c_1(D)/Tr(Box^2) = {c_D[1]/tr[2]:.10f}")
print(f"    = 14880/7200 = 31/15")
print(f"    31 = a1^2 + b1,  15 = N_gen * a1")

print(f"\n  Ratio c_2(D)/Tr(Box^4) = {c_D[2]/tr[4]:.10f}")
print(f"    = 55920/756000 = 233/3150")
print(f"    233 = F_13 (Fibonacci)")
print(f"    3150 = 2*3^2*5^2*7 = 2*a1^2*9*7/5 ... not clean")

# Perhaps the better comparison is per-vertex:
c1_per_v = c_D[1] / N   # = 124 = dim(E8)/2
tr2_per_v = tr[2] / N   # = 60 = deg*a1

print(f"\n  Per vertex:")
print(f"    c_1(D)/N = {c1_per_v:.0f} = dim(E8)/2")
print(f"    Tr(Box^2)/N = {tr2_per_v:.0f} = deg*a1")
print(f"    Ratio = {c1_per_v/tr2_per_v:.6f} = 124/60 = 31/15")

# 31/15 = (a1^2+b1)/(3*a1) = (25+6)/(15) = 31/15
# Can we relate this to alpha?
# 1/alpha ~ 137 = B = 4*a1*phi^4
# 31/15 = 2.0667... Not directly alpha.

# But: c_1(D) / Tr(Box^2) * something = 1/alpha?
# 137.08 / (31/15) = 137.08 * 15/31 = 66.33... not clean.

# Try: c_2(D) * alpha = ?
print(f"\n  c_2(D) * alpha = {c_D[2] * alpha_phys:.6f}")
print(f"    = 55920 * 0.00730 = {55920 * alpha_phys:.2f}")
print(f"    = 55920 / 137.036 = {55920/137.036:.2f}")
# 55920/137 = 408.0. And 408 = 8*51 = 8*3*17. Not obvious.

# c_2(D) / (N*B) = 55920/(120*137.08) = 55920/16449.6 = 3.400
print(f"  c_2(D) / (N*B) = {c_D[2]/(N*B_alpha):.6f}")
# 3.40 ~ ? Not clean.

# Tr(Box^2) * alpha = 7200 * 0.00730 = 52.54
print(f"  Tr(Box^2) * alpha = {tr[2]*alpha_phys:.6f}")
# 52.5 ~ ? Not clean.

# Try: Tr(Box^2) / (2*pi) = 7200/6.283 = 1146.0
print(f"\n  Tr(Box^2) / (2*pi) = {tr[2]/two_pi:.4f}")
# 1146 ~ 1/alpha * ??? 1146/137 = 8.36 ~ ??? Not clean.

# N * B * alpha = 120 * 137.08 * 0.00730 = 120 * 1 = 120
print(f"  N * B * alpha = {N * B_alpha * alpha_phys:.6f}")
print(f"    = N * 4*a1*phi^4 * alpha = N * 1/alpha * alpha = ... wait")
print(f"    B*alpha = {B_alpha * alpha_phys:.6f}")
# B*alpha = 4*a1*phi^4 * alpha ~ 137*0.0073 = 1.000... ALMOST!
# Actually: from the alpha equation, B*alpha = 1 + 2*pi*alpha^2
# B*alpha - 1 = 2*pi*alpha^2 ~ 0.000335
print(f"    B*alpha - 1 = {B_alpha*alpha_phys - 1:.10f}")
print(f"    2*pi*alpha^2 = {two_pi*alpha_phys**2:.10f}")
print(f"    Match (Vieta): {abs(B_alpha*alpha_phys - 1 - two_pi*alpha_phys**2) < 1e-10}")

# So Tr(Box^2) * alpha / N = 60 * alpha
# And 60 * alpha = 60/137.036 = 0.4378
# And: Tr(Box^2)/(N*B) = 60/137.08 = 60*alpha (approximately)
# Actually: Tr(Box^2)/(N*B) = 60/(4*a1*phi^4) = 60/(20*phi^4)
#   = 3/phi^4 = 3/(3phi+2) = 3phi^{-4}
print(f"\n  Tr(Box^2)/(N*B) = 3/phi^4 = {3/PHI**4:.10f}")
print(f"  Check: {tr[2]/(N*B_alpha):.10f}")
print(f"  Match: {abs(tr[2]/(N*B_alpha) - 3/PHI**4) < 1e-8}")

# So: Tr(Box^2) = N * B * 3/phi^4 = N * 4*a1*phi^4 * 3/phi^4 = 12*N*a1 = 12*120*5 = 7200. TAUTOLOGY.
# Because B = 4*a1*phi^4 and Tr(Box^2)/N = 12*a1 = deg*a1 = 60.
# 60/(4*a1*phi^4) = 60/(20*phi^4) = 3/phi^4. Circular.

print(f"\n  CONCLUSION: Tr(Box^2)/N = deg*a1 = 60 is a spectral identity,")
print(f"  and B = 4*a1*phi^4 = 1/alpha_0. Their ratio 3/phi^4 is algebraic")
print(f"  but does not independently produce alpha.")

# Final attempt: the spectral action ratio
# S_D/S_Box = c_2(D) / (Tr(Box^4)/2) = 55920/378000 = 233/1575
# 233/1575 = F_13 / (3^2 * 5^2 * 7) = F_13 / (9*175)
# Not clean.

print(f"\n  S_D(k=2) / S_Box(k=2) = c_2(D)/(Tr(Box^4)/2) = {c_D[2]/(tr[4]/2):.6f}")
print(f"    = 55920/378000 = {55920/378000:.6f} = 233/1575")

print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)
print(f"""
  The spectral action Tr f(Box/Lambda) has Seeley-DeWitt coefficients
  that are simple algebraic functions of a1, phi, and degree.
  The D coefficients involve additional structure (E8 dimension, Fibonacci).

  KEY FINDING: c_1(D)/N = 124 = dim(E8)/2 while Tr(Box^2)/N = 60 = deg*a1.
  Their ratio 31/15 = (a1^2+b1)/(3*a1) is clean but does not yield alpha.

  The alpha equation CANNOT be extracted from the ratio of D and Box
  spectral actions. The relationship between alpha and the spectral
  data remains the structural identity (not variational):
    2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0

  The 2*pi comes from topology (Hopf fiber holonomy), the 4*a1*phi^4
  from algebraic spectral data. These are DIFFERENT sectors and their
  combination into the alpha equation requires the KK identification,
  which this analysis does not reproduce.
""")

print("=" * 70)
print("EXP-565 COMPLETE")
print("=" * 70)
