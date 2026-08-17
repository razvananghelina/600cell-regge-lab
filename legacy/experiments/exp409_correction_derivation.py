"""
exp409_correction_derivation.py
================================
GOAL: Derive the correction factor a1/(a1-alpha_s) in
      f_4 = N^{d_ST}*alpha*alpha' * a1/(a1-alpha_s)

The correction is 1/(1 - alpha_s/a1) ~ 1.02418, which closes the
A_s gap from +3.2% to +0.8%.

APPROACHES:
  1. Algebraic structure and decomposition
  2. One-loop spectral action correction (fermionic determinant)
  3. Ihara zeta function of 600-cell graph
  4. Heat kernel / partition function at specific temperatures
  5. Spectral zeta residues
  6. RG running interpretation
  7. Connection to det'(Delta_0) prime factorization

BLIND APPROACH: Compute everything, identify the derivation.

Author: Razvan-Constantin Anghelina
Date: 2026-02-25
"""

import numpy as np
from math import sqrt, pi, log, log10, factorial, exp, cos, sin
from itertools import permutations, product as cartesian_product
from collections import defaultdict
import time

print("=" * 72)
print("EXP-409: DERIVATION OF THE CORRECTION a1/(a1 - alpha_s)")
print("=" * 72)

# Framework constants
a1 = 5
b1 = a1 + 1                            # 6
PHI = (1 + sqrt(a1)) / 2               # golden ratio
PHI_prime = (1 - sqrt(a1)) / 2         # Galois conjugate
N = factorial(a1)                       # 120
degree = 2 * b1                         # 12
N_eig = 9
rank_E8 = 8
h_E8 = a1 * b1                         # 30
N_gen = 3
d_ST = 4
alpha_s = 1 / (2*PHI**3)               # 0.11803
alpha_em = 1/137.035999  # approximate for numerical checks

# Seeley-DeWitt
A0 = 2*a1 + 1                          # 11
A1 = 2*a1**2 + 2*a1 + 2                # 62
A2 = 6*a1**2 + 15*a1 + 8               # 233
c0 = 2*N*A0                            # 2640
c1 = 2*N*A1                            # 14880
c2 = 2*N*A2                            # 55920

# 2I irrep data
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
mults = [d**2 for d in dims]
chi_phi = [1, PHI, PHI, 1, 0, -1, -1, 1-PHI, 1-PHI]

# Laplacian eigenvalues
Lk = [degree * (1 - chi_phi[k]/dims[k]) for k in range(N_eig)]

# Adjacency eigenvalues (mu_k = degree - L_k)
mu_k = [degree - Lk[k] for k in range(N_eig)]

# The correction
correction = a1 / (a1 - alpha_s)
epsilon = alpha_s / a1  # the small parameter


# =====================================================================
# PART 1: ALGEBRAIC STRUCTURE
# =====================================================================
print(f"\n{'='*72}")
print("PART 1: Algebraic Structure of the Correction")
print("=" * 72)

print(f"""
  The correction factor:
    C = a1/(a1 - alpha_s) = {a1}/({a1} - {alpha_s:.6f}) = {correction:.8f}

  Equivalent forms:
    C = 1/(1 - alpha_s/a1) = 1/(1 - epsilon)
    epsilon = alpha_s/a1 = 1/(2*a1*phi^3) = {epsilon:.8f}

  Algebraic:
    alpha_s = 1/(2*phi^3)
    epsilon = 1/(2*a1*phi^3) = 1/(10*phi^3)
    2*a1*phi^3 = {2*a1*PHI**3:.6f}
    10*phi^3 = {10*PHI**3:.6f}
""")

# Express in terms of framework quantities
print(f"  Key: 2*a1*phi^3 = a1/alpha_s = {a1/alpha_s:.6f}")
print(f"       = 10*(2+sqrt(5)) = {10*(2+sqrt(5)):.6f}")
print(f"       = 10*phi^3 = {10*PHI**3:.6f}")

# Perturbative expansion
print(f"\n  Perturbative expansion: C = sum (alpha_s/a1)^n")
for n in range(6):
    term = epsilon**n
    partial = sum(epsilon**k for k in range(n+1))
    print(f"    n={n}: (alpha_s/a1)^{n} = {term:.8e}, "
          f"partial sum = {partial:.8f} (exact: {correction:.8f})")

# Alternative algebraic forms
print(f"\n  Alternative forms:")
print(f"    C = 2*a1*phi^3/(2*a1*phi^3 - 1)")
print(f"      = {2*a1*PHI**3:.4f}/({2*a1*PHI**3:.4f} - 1)")
print(f"      = {2*a1*PHI**3:.4f}/{2*a1*PHI**3-1:.4f}")
print(f"      = {2*a1*PHI**3/(2*a1*PHI**3-1):.8f}")

# In terms of N/b1 (since N/b1 = 20 and 2*a1*phi^3 = 10*phi^3)
print(f"\n    N/b1 = {N/b1} = {N//b1}")
print(f"    a1/alpha_s = 2*a1*phi^3 = 10*phi^3")
print(f"    = (N/b1) * phi^3 / 2 = {N/b1 * PHI**3 / 2:.6f}")
print(f"    Wait: N/b1 = 20, 20*phi^3/2 = 10*phi^3 = {10*PHI**3:.6f} =/= a1/alpha_s")
print(f"    Actually: 2*a1*phi^3 = 10*phi^3. And N/(2*b1) = 10.")
print(f"    So: a1/alpha_s = (N/(2*b1)) * phi^3 = (N*phi^3)/(2*b1)")


# =====================================================================
# PART 2: ONE-LOOP SPECTRAL ACTION INTERPRETATION
# =====================================================================
print(f"\n{'='*72}")
print("PART 2: One-Loop Spectral Action Correction")
print("=" * 72)

print(f"""
  In the spectral action, the BOSONIC part is:
    S_bos = Tr(f(D^2/Lambda^2)) = f_4*Lambda^4*a_0 + f_2*Lambda^2*a_2 + f(0)*a_4

  The ONE-LOOP FERMIONIC correction adds:
    S_1-loop = -(1/2)*log det(D_A^2/Lambda^2) = (1/2)*zeta'_D(0) + ...

  This RENORMALIZES the moments f_k.

  For f_4 (the cosmological constant moment):
    f_4^eff = f_4 * (1 + delta_1-loop)

  KEY QUESTION: Does delta_1-loop = alpha_s/a1?
""")

# The fermionic correction involves Tr(D^{-2s})|_{s=0} type quantities
# Let's compute relevant spectral sums

# Spectral sums of Delta_0
nz_evals = [(Lk[k], mults[k]) for k in range(N_eig) if Lk[k] > 0.01]

# zeta(s) = sum mult * lambda^{-s}
def zeta_0(s):
    return sum(m * lam**(-s) for lam, m in nz_evals)

print(f"  Spectral zeta function of Delta_0:")
print(f"  zeta(1) = Tr(Delta_0^{{-1}}) = {zeta_0(1):.6f}")
print(f"  zeta(2) = Tr(Delta_0^{{-2}}) = {zeta_0(2):.6f}")

# The one-loop correction to the cosmological constant involves
# the spectral asymmetry and the eta invariant.
# For a GRAPH Laplacian (positive semidefinite), eta = 0.
# The correction comes from the HEAT KERNEL at finite cutoff.

# Specifically, the effective f_4 includes:
# f_4^{eff} = f_4 + (1/(4*pi^2)) * integral_0^{Lambda^2} [Tr(exp(-t*D^2))-c_0] dt/t * ...
# This is the STANDARD one-loop correction in spectral action NCG.

# Let's compute the "one-loop" type correction from the spectral data
# The key ratio is: Tr(D^{-2}) / c_0 = zeta(1) / c_0
ratio_1loop = zeta_0(1) / c0
print(f"\n  Candidate one-loop corrections:")
print(f"    zeta(1)/c_0 = {zeta_0(1)}/{c0} = {ratio_1loop:.8f}")
print(f"    alpha_s/a1 = {epsilon:.8f}")
print(f"    Ratio: {ratio_1loop/epsilon:.6f}")

# Try other combinations
print(f"\n    zeta(1)/N = {zeta_0(1)/N:.8f}")
print(f"    zeta(1)/(N-1) = {zeta_0(1)/(N-1):.8f}")
print(f"    zeta(1)/degree = {zeta_0(1)/degree:.8f}")
print(f"    zeta(2)/zeta(1) = {zeta_0(2)/zeta_0(1):.8f}")
print(f"    1/zeta(1) = {1/zeta_0(1):.8f}")
print(f"    a1/zeta(1) = {a1/zeta_0(1):.8f}")

# Check: is zeta(1) related to epsilon or correction?
print(f"\n    zeta(1) = {zeta_0(1):.6f}")
print(f"    a1*degree = {a1*degree} = 60")
print(f"    N/2 = {N/2} = 60")
print(f"    zeta(1) / (N/2) = {zeta_0(1)/(N/2):.8f}")

# Key spectral sums
print(f"\n  More spectral sums of Delta_0:")
for s in [0.5, 1, 1.5, 2, 3]:
    val = zeta_0(s)
    print(f"    zeta({s:.1f}) = {val:.6f}")

# What about the FULL Dirac spectrum (all form degrees)?
# From exp408c: the full D^2 has eigenvalues from Delta_0, Delta_1, Delta_2, Delta_3
# But we only have Delta_0 eigenvalues analytically.
# However, c_1 = 2*N*A1 gives Tr(D^2) = c_1 for the full Dirac.
# Tr(D^{-2}) for the FULL Dirac is harder.

# The INTERNAL Dirac D_F has Tr(D_F^2) = c_0 * <lambda> where <lambda> is
# related to the McKay graph eigenvalues.


# =====================================================================
# PART 3: IHARA ZETA FUNCTION
# =====================================================================
print(f"\n{'='*72}")
print("PART 3: Ihara Zeta Function of the 600-cell Graph")
print("=" * 72)

print(f"""
  The Ihara zeta function of a graph G:
    Z_G(u) = prod_{{primes}} (1 - u^{{length(p)}})^{{-1}}

  For a (q+1)-regular graph (q = degree-1 = {degree-1}):
    Z_G(u)^{{-1}} = (1-u^2)^{{E-V}} * det(I - u*A + q*u^2*I)

  For 600-cell (V={N}, E=720, q={degree-1}):
    Z_G(u)^{{-1}} = (1-u^2)^600 * prod_k (1 - u*mu_k + {degree-1}*u^2)^{{d_k^2}}

  Convergence: |u| < 1/(q) = 1/{degree-1} = {1/(degree-1):.6f}
  epsilon = alpha_s/a1 = {epsilon:.6f} < {1/(degree-1):.6f}: CONVERGES!
""")

# Adjacency eigenvalues
print(f"  Adjacency eigenvalues mu_k = degree - L_k:")
for k in range(N_eig):
    print(f"    mu_{k} = {mu_k[k]:.6f} (mult {mults[k]})")

# Compute Ihara zeta at u = epsilon = alpha_s/a1
def log_ihara_inv(u):
    """log of 1/Z_G(u)"""
    # (1-u^2)^600 part
    part1 = 600 * log(abs(1 - u**2))
    # det part: prod (1 - u*mu_k + 11*u^2)^{mult_k}
    part2 = 0
    for k in range(N_eig):
        val = 1 - u*mu_k[k] + (degree-1)*u**2
        if val > 0:
            part2 += mults[k] * log(val)
        else:
            part2 += mults[k] * log(abs(val))  # would need complex
    return part1 + part2

def ihara_det_factor(u):
    """det(I - u*A + q*u^2*I) as product over eigenvalues"""
    result = 1.0
    for k in range(N_eig):
        val = 1 - u*mu_k[k] + (degree-1)*u**2
        result *= val**mults[k]
    return result

# Evaluate at key values
print(f"\n  Ihara zeta evaluations:")
test_u = [
    ("alpha_s/a1 = epsilon", epsilon),
    ("alpha_s", alpha_s),
    ("1/N", 1/N),
    ("1/degree", 1/degree),
    ("1/(2*degree)", 1/(2*degree)),
    ("1/h_E8", 1/h_E8),
    ("epsilon/2", epsilon/2),
]

for name, u in test_u:
    if abs(u) < 1/(degree-1):
        det_val = ihara_det_factor(u)
        topological = (1-u**2)**600
        Z_inv = topological * det_val
        if abs(Z_inv) > 1e-300:
            log_Z = -log(abs(Z_inv))
            print(f"  u = {name:>25} = {u:.8f}")
            print(f"    det factor = {det_val:.6e}")
            print(f"    (1-u^2)^600 = {topological:.6e}")
            print(f"    Z_G(u) = exp({log_Z:.4f}) = {exp(log_Z):.6e}")
            # Compare with correction
            print(f"    log Z_G / log(correction) = {log_Z/log(correction):.6f}")
        else:
            print(f"  u = {name}: Z_G^(-1) too small")

# KEY TEST: Z_G(epsilon) vs correction
print(f"\n  KEY: Does Ihara zeta at u=epsilon give the correction?")
det_eps = ihara_det_factor(epsilon)
topo_eps = (1-epsilon**2)**600
Z_inv_eps = topo_eps * det_eps
log_Z_eps = -log(abs(Z_inv_eps))
Z_eps = exp(log_Z_eps)
print(f"    Z_G(epsilon) = {Z_eps:.8f}")
print(f"    correction = {correction:.8f}")
print(f"    Z_G(epsilon) / correction = {Z_eps/correction:.8f}")
print(f"    Z_G(epsilon)^(1/N) = {Z_eps**(1/N):.8f}")

# What about just the det factor?
print(f"\n    det(I - epsilon*A + 11*epsilon^2*I) = {det_eps:.8e}")
print(f"    det^(1/N) = {det_eps**(1/N):.8f}")
print(f"    det^(1/(N-1)) = {det_eps**(1/(N-1)):.8f}")


# =====================================================================
# PART 4: HEAT KERNEL AT SPECIFIC TEMPERATURES
# =====================================================================
print(f"\n{'='*72}")
print("PART 4: Heat Kernel and Partition Function")
print("=" * 72)

def Z_heat(beta):
    """Partition function Z(beta) = Tr(exp(-beta*Delta_0))"""
    return sum(mults[k] * exp(-beta * Lk[k]) for k in range(N_eig))

# Z(0) = N = 120
print(f"  Z(0) = {Z_heat(0):.1f} = N")
print(f"  Z(inf) -> {mults[0]} (zero mode multiplicity)")

# The RETURN PROBABILITY at the identity:
# P(t) = Z(t)/N = probability of returning to start at time t
# For a random walk on the Cayley graph.

print(f"\n  Looking for beta where Z(beta)/Z(0) = 1 - epsilon = {1-epsilon:.8f}")
print(f"  or equivalently Z(beta) = N*(1-epsilon) = {N*(1-epsilon):.6f}")

# For small beta: Z(beta) ~ N - beta*Tr(D^2) + ...
# Z(beta)/N ~ 1 - beta*(c_1/c_0)/N * something...
# Actually Z(beta) = sum d_k^2 * exp(-beta*L_k) for Delta_0
# So Z(beta)/N = (1/N) * sum d_k^2 * exp(-beta*L_k)
# For small beta: Z(beta)/N ~ 1 - beta*<L> where <L> = sum d_k^2*L_k/N

avg_L = sum(mults[k] * Lk[k] for k in range(N_eig)) / N
print(f"\n  <L> = sum d_k^2 * L_k / N = {avg_L:.6f}")
print(f"  = degree * (1 - sum d_k^2*chi_k/(N*d_k)) / 1")

# Actually: <L> = sum d_k^2 * L_k / N = degree - (degree/N)*sum d_k*chi_k
sum_dk_chi = sum(dims[k] * chi_phi[k] for k in range(N_eig))
print(f"  sum d_k * chi_k(phi-class) = {sum_dk_chi:.6f}")
print(f"  <L> = degree - degree*sum(d_k*chi_k)/N = {degree} - {degree}*{sum_dk_chi:.4f}/{N}")
print(f"       = {degree - degree*sum_dk_chi/N:.6f}")

# For Z/N = 1-epsilon: beta = epsilon/<L>
beta_target = epsilon / avg_L
print(f"\n  beta for Z/N = 1-epsilon: beta ~ epsilon/<L> = {beta_target:.8f}")
print(f"  Z(beta_target)/N = {Z_heat(beta_target)/N:.8f}")
print(f"  Target: 1 - epsilon = {1-epsilon:.8f}")

# More precise: find beta where Z(beta)/N = (a1-alpha_s)/a1
from_target = 1 - epsilon  # = (a1-alpha_s)/a1
beta_lo, beta_hi = 0.0, 0.1
for _ in range(200):
    beta_mid = (beta_lo + beta_hi) / 2
    if Z_heat(beta_mid)/N > from_target:
        beta_lo = beta_mid
    else:
        beta_hi = beta_mid
beta_exact = (beta_lo + beta_hi) / 2
print(f"\n  EXACT: Z(beta*)/N = (a1-alpha_s)/a1 when beta* = {beta_exact:.10f}")
print(f"  Verification: Z({beta_exact:.8f})/N = {Z_heat(beta_exact)/N:.10f}")

# Is beta* a framework quantity?
print(f"\n  Comparing beta* = {beta_exact:.8f} with framework:")
print(f"    epsilon/<L> = {epsilon/avg_L:.8f} (linear approx)")
print(f"    epsilon/degree = {epsilon/degree:.8f}")
print(f"    1/(N*degree) = {1/(N*degree):.8f}")
print(f"    alpha_s/(a1*degree) = {alpha_s/(a1*degree):.8f}")
print(f"    1/(a1*N) = {1/(a1*N):.8f}")
print(f"    1/(c0) = {1/c0:.8f}")
print(f"    1/(c1) = {1/c1:.8f}")
print(f"    epsilon/(2*b1) = {epsilon/(2*b1):.8f}")

# Ratio to simple quantities
print(f"\n  beta* * <L> = {beta_exact * avg_L:.8f} (should be ~ epsilon = {epsilon:.8f})")
print(f"  beta* * degree = {beta_exact * degree:.8f}")
print(f"  beta* * c0 = {beta_exact * c0:.8f}")
print(f"  beta* * N = {beta_exact * N:.8f}")


# =====================================================================
# PART 5: THE REGULARIZED TRACE AND MULTIPLICATIVE ANOMALY
# =====================================================================
print(f"\n{'='*72}")
print("PART 5: Regularized Traces and the Multiplicative Anomaly")
print("=" * 72)

print(f"""
  The multiplicative anomaly (Wodzicki residue) for two operators A, B:
    log det(AB) = log det(A) + log det(B) + anomaly(A,B)

  For the spectral action, the relevant quantity is:
    Tr(f(D^2/Lambda^2)) vs Tr(f(D^2)) / Lambda^{d_ST}

  The "finite-size correction" to f_4 comes from the discrete nature
  of the spectrum. On a lattice with N^{{d_ST}} points:

    f_4^{{lattice}} = f_4^{{continuum}} * (1 + O(1/N))

  The leading finite-size correction for a gauge theory is:
    delta f_4 / f_4 = c * g^2 / (4*pi)^2

  where g is the gauge coupling and c is a coefficient.
""")

# In our framework, g^2 = 4*pi*alpha_s for the strong coupling
g_s_sq = 4 * pi * alpha_s
print(f"  g_s^2 = 4*pi*alpha_s = {g_s_sq:.6f}")
print(f"  g_s^2/(4*pi)^2 = alpha_s/(4*pi) = {alpha_s/(4*pi):.8f}")
print(f"  alpha_s/a1 = {epsilon:.8f}")
print(f"  Ratio: (alpha_s/a1) / (alpha_s/(4*pi)) = {epsilon/(alpha_s/(4*pi)):.6f}")
print(f"       = 4*pi/a1 = {4*pi/a1:.6f}")

# So if the one-loop correction is alpha_s/(4*pi) * something,
# we need "something" = 4*pi/a1 to get alpha_s/a1.
# 4*pi/a1 = 4*pi/5 = 2.513...

# Actually, the standard one-loop correction in lattice gauge theory is:
# delta f_4 / f_4 = -b_0 * alpha_s / (2*pi)
# where b_0 is the one-loop beta function coefficient.

# For SU(3): b_0 = 11 - 2/3*n_f = 11 - 2/3*6 = 7
# For SU(2): b_0 = 22/3 - 2/3*n_f
# For our framework: b_0 might be different.

# But our correction is POSITIVE (f_4 increases), so it's not the usual
# asymptotic freedom correction. Let's think differently.

# In the spectral action, the correction to f_4 comes from the
# SPECTRAL GEOMETRY of the internal space. The relevant quantity is:
# f_4^{eff}/f_4 = 1 + Tr(D_F^{-2}) * (something) + ...

# What is Tr(D_F^{-2}) for the internal Dirac operator?
# The internal Dirac D_F has eigenvalues related to the McKay graph.
# If we interpret D_F^2 as Delta_0 on the 600-cell:
# Tr(Delta_0^{-1}) = zeta(1) = sum d_k^2 / L_k (k>0)

z1 = zeta_0(1)
print(f"\n  Tr(Delta_0^{{-1}}) = zeta(1) = {z1:.6f}")
print(f"  zeta(1) / degree = {z1/degree:.6f}")
print(f"  zeta(1) * alpha_s = {z1*alpha_s:.6f}")
print(f"  zeta(1) * alpha_s / a1 = {z1*alpha_s/a1:.8f}")
print(f"  zeta(1) / (N-1) = {z1/(N-1):.8f} (avg 1/lambda)")
print(f"  degree * zeta(1)/(N-1) = {degree*z1/(N-1):.6f}")

# HYPOTHESIS: The correction involves the HEAT KERNEL TRACE at the
# coupling scale, normalized by the volume.

# The spectral action gives (schematically):
# f_4^{eff} = f_4 * [1 + alpha_s * F(spectral data)]
# where F is some functional of the spectrum.

# If F = 1/a1, we get the correction.

# Check: in the Connes-Chamseddine framework, the one-loop effective action
# modifies the spectral action through:
# S_eff = Tr(f(D^2/Lambda^2)) + (1/2)*Tr(log(D_F^2))
# The second term is the one-loop fermionic correction.

# log det'(D_F^2) = sum d_k^2 * log(L_k) = log det'(Delta_0)
# This was computed in exp408c/d:
log_det_prime = sum(mults[k] * log(Lk[k]) for k in range(1, N_eig))
print(f"\n  log det'(Delta_0) = {log_det_prime:.6f}")
print(f"  log det' / N = {log_det_prime/N:.6f}")
print(f"  log det' / c0 = {log_det_prime/c0:.6f}")
print(f"  exp(log det'/c0) = {exp(log_det_prime/c0):.8f}")
print(f"  correction = {correction:.8f}")


# =====================================================================
# PART 6: RG RUNNING INTERPRETATION
# =====================================================================
print(f"\n{'='*72}")
print("PART 6: Renormalization Group Running")
print("=" * 72)

print(f"""
  HYPOTHESIS: f_4 runs from the bare value N^4/(2*pi) at the cutoff
  to an effective value at the physical scale.

  In the RG framework, the running of f_4 is:
    d(ln f_4)/d(ln mu) = gamma_f4 = some anomalous dimension

  If gamma_f4 = alpha_s/a1, then over one "e-folding" of running:
    f_4(mu) = f_4(Lambda) * exp(alpha_s/a1 * ln(Lambda/mu))

  For Lambda/mu ~ e (one e-folding):
    f_4(mu) = f_4(Lambda) * exp(alpha_s/a1) ~ f_4(Lambda) * (1 + alpha_s/a1)

  But we need the RESUMMED form: f_4 * a1/(a1-alpha_s) = f_4/(1-alpha_s/a1)

  This corresponds to running from mu=0 to mu=Lambda with:
    f_4^{{-1}}(mu) = f_4^{{-1}}(Lambda) - alpha_s/(a1*f_4(Lambda)) * (...)

  More precisely, the EXACT RG equation:
    f_4(mu)^{{-1}} = f_4(Lambda)^{{-1}} * (1 - alpha_s/a1)

  gives:
    f_4(mu) = f_4(Lambda) / (1 - alpha_s/a1) = f_4 * a1/(a1-alpha_s)
""")

# The beta function interpretation:
# In the SM, the one-loop beta functions are:
# b_1 = 41/6, b_2 = -19/6, b_3 = -7

# In our framework:
# beta_i/alpha_i propto b_i (the one-loop coefficient)
# The RATIO of beta functions: b_3:b_2:b_1 = 8:5:2 in our framework

# The correction involves alpha_s, suggesting it's from the STRONG sector.
# The coefficient 1/a1 = 1/5 could be:
# - 1/a1 from the normalization of the Cayley graph
# - Related to the Dynkin index of the representation

print(f"  One-loop beta function coefficients (framework):")
print(f"    b_3 : b_2 : b_1 = 8 : 5 : 2")
print(f"    Sum = 15 = a1*N_gen")
print(f"    b_3 = 8 = rank(E8)")
print(f"    b_3/N_gen = 8/3 = {8/3:.4f}")
print(f"    b_3/(a1*b1) = 8/30 = {8/30:.6f}")

# Could the correction involve b_3?
# alpha_s * b_3 / (something) = alpha_s/a1?
# b_3 / something = 1/a1 => something = a1*b_3 = 5*8 = 40
# 40 = N/3 = N/N_gen. Interesting!
print(f"\n  If correction = alpha_s * b_3 / (a1*b_3):")
print(f"    = alpha_s/a1 (trivially)")
print(f"  But b_3/(a1*b_3) = 1/a1 is tautological.")

# Let's try: correction from EACH gauge sector
# The total one-loop correction would be:
# delta = sum_i alpha_i * b_i * c_i
# where c_i is some coefficient.

# For the correction to be alpha_s/a1:
# We need only the strong sector contribution (alpha_s),
# with coefficient b_3_eff = 1/a1 = 1/5.

# Actually, 1/a1 = 1/5. And the Dynkin index of the fundamental
# representation of SU(3) is T(fund) = 1/2.
# But for 2I: the "Dynkin index" analog... hmm.

# Let me think about this differently.
# The 600-cell has degree 12. The spectrum of the normalized Laplacian
# L_norm = Delta/degree has eigenvalues in [0, 2].
# The "gap" is lambda_1/degree = L_1/degree = 2.292/12 = 0.191
print(f"\n  Spectral gap: lambda_1/degree = {Lk[1]/degree:.6f}")
print(f"  alpha_s = {alpha_s:.6f}")
print(f"  Ratio: gap/alpha_s = {Lk[1]/degree/alpha_s:.6f}")

# Very interesting: lambda_1/degree ~ 1.6*alpha_s
# Not exact. Let me check lambda_1/(2*degree):
print(f"  lambda_1/(2*degree) = {Lk[1]/(2*degree):.6f}")
print(f"  alpha_s = {alpha_s:.6f}")
print(f"  Ratio: {Lk[1]/(2*degree)/alpha_s:.6f}")

# lambda_1 = 12 - 6*phi = 12 - 6*phi
# alpha_s = 1/(2*phi^3)
# lambda_1 * alpha_s = (12-6*phi)/(2*phi^3)
print(f"\n  lambda_1 * alpha_s = {Lk[1]*alpha_s:.8f}")
print(f"  lambda_1 * alpha_s * a1 = {Lk[1]*alpha_s*a1:.8f}")


# =====================================================================
# PART 7: THE EFFECTIVE ACTION APPROACH
# =====================================================================
print(f"\n{'='*72}")
print("PART 7: Effective Action -- Direct Computation")
print("=" * 72)

print(f"""
  The most direct approach: compute the ONE-LOOP effective potential
  on the 600-cell and extract the correction to f_4.

  The effective potential at one loop:
    V_eff = V_tree + V_1-loop

  where V_1-loop = (1/(64*pi^2)) * Tr[M^4 * (ln(M^2/mu^2) - c)]

  For the spectral action, M^2 = D_F^2 (the internal Dirac squared).
  The eigenvalues of D_F^2 are the L_k with multiplicities d_k^2.

  V_1-loop = (1/(64*pi^2)) * sum_k d_k^2 * L_k^2 * (ln(L_k/mu^2) - 3/2)

  The ratio V_1-loop / V_tree at the cosmological constant level
  determines the correction to f_4.
""")

# Tree-level CC contribution (from f_4 * Lambda^4 * a_0):
# V_tree ~ f_4 * a_0 = f_4 * c_0

# One-loop CC contribution:
# V_1-loop ~ (1/(64*pi^2)) * Tr(D_F^4) * (ln terms)
# Tr(D_F^4) = sum d_k^2 * L_k^2 = Tr(Delta_0^2) = ...

tr_L2 = sum(mults[k] * Lk[k]**2 for k in range(N_eig))
tr_L4 = sum(mults[k] * Lk[k]**4 for k in range(N_eig))

print(f"  Tr(Delta_0^2) = sum d_k^2 * L_k^2 = {tr_L2:.2f}")
print(f"  Tr(Delta_0^4) = sum d_k^2 * L_k^4 = {tr_L4:.2f}")

# The Coleman-Weinberg one-loop correction
# At mu^2 = <L^2> (the natural scale):
avg_L2 = tr_L2 / (N-1)  # average over nonzero modes
mu2 = avg_L2
print(f"  <L^2> (avg over nonzero) = {avg_L2:.6f}")
print(f"  <L> = {avg_L:.6f}")

# V_1-loop / V_tree at the CC level
# V_tree = f_4 * c_0 (schematic)
# V_1-loop = (1/(64*pi^2)) * sum d_k^2 * L_k^2 * (ln(L_k/mu) - 3/2)

V_1loop_sum = sum(mults[k] * Lk[k]**2 * (log(Lk[k]/mu2) - 1.5)
                  for k in range(1, N_eig))
print(f"\n  V_1-loop sum = {V_1loop_sum:.6f}")
print(f"  / (64*pi^2) = {V_1loop_sum/(64*pi**2):.6f}")

# Normalize by V_tree = f_4 * c_0
# If f_4 = N^4/(2*pi), V_tree ~ N^4/(2*pi) * c_0
V_tree = N**4 / (2*pi) * c0
ratio_V = V_1loop_sum / (64*pi**2) / V_tree
print(f"  V_tree = f_4 * c_0 = {V_tree:.4e}")
print(f"  V_1-loop/V_tree = {ratio_V:.8e}")
print(f"  alpha_s/a1 = {epsilon:.8e}")
print(f"  Ratio: V_1-loop/V_tree / (alpha_s/a1) = {ratio_V/epsilon:.6f}")


# =====================================================================
# PART 8: THE CRITICAL INSIGHT -- LATTICE PLAQUETTE CORRECTION
# =====================================================================
print(f"\n{'='*72}")
print("PART 8: Lattice Plaquette Correction")
print("=" * 72)

print(f"""
  In Wilson's lattice gauge theory, the partition function is:
    Z = prod_links int dU * exp(beta * sum_plaq Re Tr U_p / N_c)

  At STRONG COUPLING (beta = 0): Z = |G|^{{#links}} = N^d per vertex
  At WEAK COUPLING (beta >> 1): continuum gauge theory emerges

  The LEADING CORRECTION to the strong-coupling result is:
    Z = N^d * (1 + beta * <plaquette> + ...)

  For our framework:
    beta = 2*N_c / g^2 = 1/(2*alpha_s)  [for "SU(N_c)"]

  QUESTION: What plays the role of the plaquette expectation value?

  A PLAQUETTE is a closed path of length 4 on the lattice.
  On the Cayley graph of 2I, the number of 4-cycles passing through
  a vertex is related to Tr(A^4)/N.
""")

# Number of closed walks of length 4 from a vertex
# = (A^4)_{ii} = sum_k d_k^2 * mu_k^4 / N (by character theory)
walks_4 = sum(mults[k] * mu_k[k]**4 for k in range(N_eig)) / N
print(f"  Walks of length 4 per vertex: (1/N)*Tr(A^4) = {walks_4:.2f}")
print(f"  Walks of length 2 per vertex: (1/N)*Tr(A^2) = {sum(mults[k]*mu_k[k]**2 for k in range(N_eig))/N:.2f}")
print(f"  degree = {degree}")
print(f"  degree^2 = {degree**2}")

# The "plaquette" in our theory could be:
# <plaq> = (1/N) * Tr(A^4) / degree^4 (normalized)
plaq_norm = walks_4 / degree**4
print(f"\n  Normalized plaquette: <plaq> = walks_4/degree^4 = {plaq_norm:.8f}")
print(f"  beta * <plaq> = (1/(2*alpha_s)) * {plaq_norm:.6f} = {plaq_norm/(2*alpha_s):.8f}")
print(f"  alpha_s/a1 = {epsilon:.8f}")
print(f"  Ratio: {plaq_norm/(2*alpha_s)/epsilon:.6f}")

# Try: strong coupling expansion with beta = 1/(something)
# Z_1vertex / N^d = 1 + correction
# correction = alpha_s/a1 => what beta gives this?

# Actually in the strong coupling expansion:
# Z = N^d * exp(d*beta/N + ...) for U(1) on a d-dim lattice
# For a general group: Z = N^d * (1 + d*(d-1)/2 * (1/N^2) * ...)

# The SIMPLEST model: the correction comes from the fact that NOT
# all gauge configurations contribute equally. The "effective" number
# of configurations is N^d * (1 - alpha_s/a1) instead of N^d,
# or equivalently the EFFECTIVE f_4 is N^d/(1-alpha_s/a1) * alpha*alpha'.


# =====================================================================
# PART 9: COMBINATORIAL IDENTITY SEARCH
# =====================================================================
print(f"\n{'='*72}")
print("PART 9: Systematic Search for alpha_s/a1 in Spectral Data")
print("=" * 72)

target = epsilon  # = alpha_s/a1

# Search among ratios of spectral invariants
spectral_quantities = {
    'zeta(1)': zeta_0(1),
    'zeta(2)': zeta_0(2),
    'Tr(L^2)': tr_L2,
    'Tr(L^4)': tr_L4,
    '<L>': avg_L,
    '<L^2>': avg_L2,
    'degree': float(degree),
    'N': float(N),
    'N-1': float(N-1),
    'c0': float(c0),
    'c1': float(c1),
    'c2': float(c2),
    'A0': float(A0),
    'A1': float(A1),
    'A2': float(A2),
    'h_E8': float(h_E8),
    'rank_E8': float(rank_E8),
    'a1': float(a1),
    'b1': float(b1),
    'N_gen': float(N_gen),
    'phi': PHI,
    'phi^2': PHI**2,
    'phi^3': PHI**3,
    'phi^4': PHI**4,
}

print(f"  Target: alpha_s/a1 = {target:.8f}")
print(f"\n  Searching for ratios A/B = alpha_s/a1:")
matches = []
for name_a, val_a in spectral_quantities.items():
    for name_b, val_b in spectral_quantities.items():
        if name_a != name_b and abs(val_b) > 1e-10:
            ratio = val_a / val_b
            if abs(ratio - target) / target < 0.01:  # within 1%
                rel_err = (ratio - target) / target * 100
                matches.append((rel_err, name_a, name_b, ratio))

matches.sort(key=lambda x: abs(x[0]))
for rel_err, na, nb, ratio in matches[:15]:
    print(f"    {na}/{nb} = {ratio:.8f} ({rel_err:+.4f}%)")

# Also search for alpha_s/a1 as a simple combination
print(f"\n  Searching for a*X + b*Y = alpha_s/a1 (a,b in {{-2..2}}/N-type):")
simple_vals = {
    '1/N': 1/N,
    '1/c0': 1/c0,
    '1/c1': 1/c1,
    '1/c2': 1/c2,
    '1/(N*degree)': 1/(N*degree),
    'zeta(1)/N^2': zeta_0(1)/N**2,
    'zeta(2)/N': zeta_0(2)/N,
    '<L>/N^2': avg_L/N**2,
    'A0/c0': A0/c0,
}

for name, val in simple_vals.items():
    ratio = target / val
    print(f"    alpha_s/a1 / ({name}) = {ratio:.6f}")
    # Check if ratio is close to a simple integer or fraction
    for num in range(1, 20):
        for den in range(1, 20):
            if abs(ratio - num/den) < 0.02:
                print(f"      ~ {num}/{den} = {num/den:.6f} (err {(ratio-num/den)/ratio*100:.3f}%)")


# =====================================================================
# PART 10: THE KEY IDENTITY
# =====================================================================
print(f"\n{'='*72}")
print("PART 10: Summary and Key Identity")
print("=" * 72)

# Let me check one more thing: the AVERAGE of 1/L_k weighted by d_k
# (not d_k^2)
avg_inv_L_dim = sum(dims[k] / Lk[k] for k in range(1, N_eig))
avg_inv_L_dim2 = sum(dims[k]**2 / Lk[k] for k in range(1, N_eig))
print(f"  sum d_k / L_k (k>0) = {avg_inv_L_dim:.6f}")
print(f"  sum d_k^2 / L_k (k>0) = zeta(1) = {avg_inv_L_dim2:.6f}")
print(f"  sum 1/L_k (k>0) = {sum(1/Lk[k] for k in range(1, N_eig)):.6f}")

# Character-weighted trace
# Tr_reg(Delta^{-1}) = zeta(1) = sum d_k^2/L_k
# Tr_fund(Delta^{-1}) = sum d_k/L_k (if "fund" means dimension-weighted)

# Now: 1/(2*phi^3) = alpha_s
# And: (1/a1) * alpha_s = epsilon
# So: 1/(2*a1*phi^3) = 1/(10*phi^3)

inv_10phi3 = 1/(10*PHI**3)
print(f"\n  1/(10*phi^3) = {inv_10phi3:.8f} = alpha_s/a1")

# Check: 10*phi^3 = 10*(2+sqrt(5))
# phi^3 = phi^2 * phi = (phi+1)*phi = phi^2+phi = (phi+1)+phi = 2*phi+1
# So 10*phi^3 = 10*(2*phi+1) = 20*phi+10
val_10phi3 = 20*PHI + 10
print(f"  10*phi^3 = 20*phi + 10 = {val_10phi3:.6f}")
print(f"  = 10*(2*phi+1) = 10*phi^3")

# And a1/(a1-alpha_s) = (2*a1*phi^3)/(2*a1*phi^3-1)
# = (20*phi+10)/(20*phi+9)
denom = 20*PHI + 9
print(f"\n  a1/(a1-alpha_s) = (20*phi+10)/(20*phi+9)")
print(f"  = {(20*PHI+10)/denom:.8f}")
print(f"  Check: {correction:.8f}")

# Is 20*phi+9 = 20*phi+(2*a1-1)?
print(f"\n  20*phi + 9 = 2*a1*phi^3 - 1")
print(f"  = 2*a1*(2*phi+1) - 1 = 4*a1*phi + 2*a1 - 1 = {4*a1*PHI+2*a1-1:.6f}")
print(f"  Check: {2*a1*PHI**3-1:.6f}")

# FACTORED FORM:
# 2*a1*phi^3 - 1 = 10*phi^3 - 1
# = 10*(2+sqrt(5)) - 1 = 19 + 10*sqrt(5)
print(f"\n  2*a1*phi^3 - 1 = 19 + 10*sqrt(5) = {19+10*sqrt(5):.6f}")

# Is 19+10*sqrt(5) related to something in the 600-cell?
# (3+sqrt(5))^2 = 14+6*sqrt(5). No.
# (4+sqrt(5))^2 = 21+8*sqrt(5). No.
# 5*(3+2*sqrt(5)) = 15+10*sqrt(5). Close but not 19.
# 4+10*sqrt(5)+15 = 19+10*sqrt(5). Hmm, 4 = d_ST, 15 = degree+3?

# phi^5 = 5*phi+3 = 11.09...
# 2*phi^5 = 10*phi+6 = 22.18...
# Not matching.

# But we know: a1/(a1-alpha_s) = 1 + alpha_s/a1 + (alpha_s/a1)^2 + ...
# This is a GEOMETRIC SERIES. The physical interpretation is:
# Each power of alpha_s/a1 represents one more "loop" of the
# strong interaction correction.

# FINAL ASSESSMENT
print(f"""
{'='*72}
FINAL ASSESSMENT
{'='*72}

  The correction a1/(a1-alpha_s) = 1/(1-alpha_s/a1) = {correction:.6f}

  ALGEBRAIC STRUCTURE:
    = (20*phi+10)/(20*phi+9) = 2*a1*phi^3/(2*a1*phi^3-1)
    = a1/alpha_s / (a1/alpha_s - 1)

  INTERPRETATION CANDIDATES:

  1. ONE-LOOP RG RUNNING: f_4^{{eff}} = f_4 / (1 - alpha_s/a1)
     This is the resummation of the perturbative series
     sum_n (alpha_s/a1)^n. Status: STRUCTURAL (no explicit derivation
     of why the anomalous dimension is exactly alpha_s/a1).

  2. LATTICE PLAQUETTE CORRECTION: The effective number of gauge
     configurations at a vertex is not N^d but N^d/(1-alpha_s/a1).
     Status: HEURISTIC.

  3. SPECTRAL ZETA REGULARIZATION: None of the spectral invariants
     of Delta_0 directly gives alpha_s/a1. The correction does NOT
     come from the 600-cell graph spectrum.

  4. ALGEBRAIC: a1/(a1-alpha_s) = a1*2*phi^3/(2*a1*phi^3-1).
     The numerator is a1/alpha_s = 10*phi^3 (= D^2 = total quantum
     dimension squared of Fibonacci TQFT... wait!)
""")

# WAIT -- check this!
D2_TQFT = a1 + sqrt(a1)  # = 5 + sqrt(5)
D2_TQFT_prime = a1 - sqrt(a1)  # = 5 - sqrt(5)
D2_product = D2_TQFT * D2_TQFT_prime  # = a1*(a1-1) = 20

print(f"  TQFT check:")
print(f"    D^2 = a1 + sqrt(a1) = {D2_TQFT:.6f}")
print(f"    D'^2 = a1 - sqrt(a1) = {D2_TQFT_prime:.6f}")
print(f"    D^2 * D'^2 = {D2_product:.6f} = a1*(a1-1) = {a1*(a1-1)}")
print(f"    N/b1 = {N/b1} = 20")

# a1/alpha_s = 2*a1*phi^3 = 10*phi^3
# D^2 = 5+sqrt(5) = 5+2.236 = 7.236
# 10*phi^3 = 42.36
# These are not the same.

# But: D^2/2 = (a1+sqrt(a1))/2 = phi^2*... hmm
# phi^2 = phi+1 = (3+sqrt(5))/2 = 2.618
# 10*phi^3 = 10*phi*phi^2 = 10*phi*(phi+1) = 10*(phi^2+phi) = 10*(2*phi+1)

# Let me try: a1/(a1-alpha_s) in terms of TQFT quantities
# D^2 = a1+sqrt(a1), so sqrt(a1) = D^2-a1
# phi = (1+sqrt(a1))/2 = (1+D^2-a1)/2 = (D^2-4)/2 (since a1=5)
# phi^3 = ((D^2-4)/2)^3 ... getting complicated

# Let me try a COMPLETELY DIFFERENT angle:
# What if the correction comes from the NUMBER OF E-FOLDS?

N_e = N / 2  # = 60 (number of e-folds, DERIVED)
# The CMB amplitude A_s depends on N_e.
# If the exact N_e is not 60 but 60*(1-alpha_s/a1):

N_e_eff = N_e * (1 - epsilon)
print(f"\n  E-fold correction:")
print(f"    N_e = N/2 = {N_e}")
print(f"    N_e_eff = N_e*(1-alpha_s/a1) = {N_e_eff:.4f}")
print(f"    N_e_eff = 60*(1-{epsilon:.6f}) = {N_e_eff:.4f}")
print(f"    From exp408: N_e for exact match = 59.06")
print(f"    N_e_eff = {N_e_eff:.4f}")
print(f"    Difference from 59.06: {N_e_eff - 59.06:.4f}")

# Hmm, 60*(1-0.02361) = 58.58, not 59.06.
# So the N_e correction is NOT simply (1-alpha_s/a1).

# What if the correction enters through A_s differently?
# A_s propto 1/f_4 (from exp408 analysis)
# So f_4 * a1/(a1-alpha_s) means A_s * (a1-alpha_s)/a1
# A_s_corr = A_s_bare * (1-alpha_s/a1) = A_s_bare * 0.9764

A_s_bare = 2.167e-9  # from exp408
A_s_corr = A_s_bare * (1 - epsilon)
A_s_obs = 2.1e-9
print(f"\n  A_s correction:")
print(f"    A_s(bare) = {A_s_bare:.3e}")
print(f"    A_s(corr) = A_s * (1-alpha_s/a1) = {A_s_corr:.3e}")
print(f"    A_s(obs) = {A_s_obs:.3e}")
print(f"    (corr-obs)/obs = {(A_s_corr-A_s_obs)/A_s_obs*100:.2f}%")

# Even better: check if A_s propto N_e^2 / f_4
# If N_e -> N_e_eff and f_4 -> f_4_eff:
# A_s propto N_e^2/f_4
# If ONLY f_4 gets corrected (N_e stays at 60):
# A_s_eff = A_s_bare * (a1-alpha_s)/a1
print(f"\n  If only f_4 is corrected:")
print(f"    A_s = {A_s_bare*(1-epsilon):.4e} ({(A_s_bare*(1-epsilon)/A_s_obs-1)*100:+.2f}%)")
print(f"  If A_s propto N_e^2/f_4 and N_e also shifts:")
for delta_Ne in [0.0, -0.5, -0.94, -1.0, -1.42]:
    Ne_trial = 60 + delta_Ne
    A_s_trial = A_s_bare * (Ne_trial/60)**2 * (1 - epsilon)
    print(f"    N_e={Ne_trial:.2f}: A_s = {A_s_trial:.4e} ({(A_s_trial/A_s_obs-1)*100:+.2f}%)")


# =====================================================================
# PART 11: THE PLANCHEREL-LATTICE CONNECTION
# =====================================================================
print(f"\n{'='*72}")
print("PART 11: Plancherel Measure and Gauge Volume")
print("=" * 72)

print(f"""
  The Plancherel formula for 2I:
    sum d_k^2 / |G| = 1

  For d_ST copies (one per lattice direction):
    (sum d_k^2 / |G|)^{{d_ST}} = 1^{{d_ST}} = 1

  But the DISCRETE sum differs from the INTEGRAL (Weyl integration):
    For SU(2): integral |chi_j(g)|^2 dg = 1 (each irrep, Haar normalized)
    For 2I:    (1/|G|) sum_g |chi_j(g)|^2 = 1 (each irrep, counting measure)

  The lattice-to-continuum transition involves:
    |G|^{{d_ST}} (counting) -> Vol(G)^{{d_ST}} (Haar) * normalization

  Vol(SU(2)) = 2*pi^2 (volume of S^3)
  |2I| = 120

  Ratio per direction: |2I|/Vol(SU(2)) = 120/(2*pi^2) = {N/(2*pi**2):.6f}

  For d_ST = 4 directions:
    (|2I|/Vol(SU(2)))^4 = {(N/(2*pi**2))**4:.4e}
""")

# The correction factor: when going from DISCRETE to CONTINUUM,
# the effective counting changes.
# Could the correction 1/(1-alpha_s/a1) come from this transition?

# The NORMALIZED character sum (Plancherel with phi-class characters):
planch_phi = sum(dims[k] * chi_phi[k] for k in range(N_eig))
print(f"  sum d_k * chi_k(phi) = {planch_phi:.6f}")
print(f"  This should be 0 for non-identity class by orthogonality")

# The REPRESENTATION ZETA function:
# Z_rep(s) = sum_k d_k^{-s} (sum over irreps)
print(f"\n  Representation zeta function:")
for s in [0, 1, 2, 3, 4]:
    Z_rep = sum(d**(-s) for d in dims)
    print(f"    Z_rep({s}) = sum d_k^(-{s}) = {Z_rep:.6f}")

# Check: Z_rep(0) = 9 = N_eig (number of irreps)
# Z_rep(-2) = sum d_k^2 = N = 120
# Z_rep(-4) = sum d_k^4 = 2628

# Is alpha_s/a1 related to Z_rep at some value?
# 1/Z_rep(-2) = 1/N = 0.00833
# 1/Z_rep(-4) = 1/2628 = 0.000381
# Z_rep(2) = sum 1/d_k^2 = 1+1/4+1/9+1/16+1/25+1/36+1/16+1/4+1/9
Z_rep_2 = sum(1/d**2 for d in dims)
print(f"\n    Z_rep(2) = sum 1/d_k^2 = {Z_rep_2:.6f}")
print(f"    Z_rep(2)/N = {Z_rep_2/N:.8f}")
print(f"    alpha_s/a1 = {epsilon:.8f}")
print(f"    Z_rep(2)/(a1*N) = {Z_rep_2/(a1*N):.8f}")

# HMMMM Z_rep(2) = 1 + 1/4 + 1/9 + 1/16 + 1/25 + 1/36 + 1/16 + 1/4 + 1/9
# = 1 + 0.25 + 0.1111 + 0.0625 + 0.04 + 0.02778 + 0.0625 + 0.25 + 0.1111
# = 1.9149...

# Interesting: Z_rep(2) = 1.915 ~ ln(N^4)/(4*ln(N)) * ... no
# 4*ln(N) = 19.15. And Z_rep(2)*10 = 19.15!
print(f"    Z_rep(2) * 10 = {Z_rep_2*10:.6f}")
print(f"    4*ln(N) = {4*log(N):.6f}")
print(f"    !!! Z_rep(2) * 10 = 4*ln(N) to within...")
print(f"    Difference: {abs(Z_rep_2*10 - 4*log(N)):.6f}")
# Not exact. 19.149 vs 19.150. Close but probably coincidental.

print(f"\n{'='*72}")
print("EXPERIMENT COMPLETE")
print("=" * 72)
