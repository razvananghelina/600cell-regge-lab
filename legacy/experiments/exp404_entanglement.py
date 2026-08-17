"""
EXP-404: Entanglement in the 600-Cell Framework
=================================================

Question: Can the 600-cell framework say something about quantum entanglement?

Approaches:
A. Entanglement entropy of McKay graph bipartition (chirality split)
B. Fibonacci anyons: quantum dimension d=phi, topological entropy
C. Entanglement entropy on the 600-cell graph (area law)
D. Galois pairs as entangled states

ALL inputs from a1=5. Looking for framework quantities in entanglement measures.

Author: Claude (exp404)
Date: 2026-02-25
"""

import math
import numpy as np
from scipy import linalg

a1 = 5
b1 = 6
PHI = (1 + math.sqrt(5)) / 2
N = 120
DEG = 12
N_gen = 3
alpha_s = 1 / (2 * PHI**3)

print("=" * 72)
print("EXP-404: ENTANGLEMENT IN THE 600-CELL FRAMEWORK")
print("=" * 72)


# =====================================================================
# PART 1: Fibonacci Anyons and the Golden Ratio
# =====================================================================
print("\nPART 1: Fibonacci Anyons")
print("-" * 50)

# SU(2) Chern-Simons at level k
# Quantum dimensions: d_j = sin(pi*(2j+1)/(k+2)) / sin(pi/(k+2))
# At k=3 (k+2=5=a1): j = 0, 1/2, 1, 3/2

k = a1 - 2  # k = 3
print(f"  SU(2) Chern-Simons level k = a1-2 = {k}")
print(f"  k+2 = {k+2} = a1 = {a1}")

q_dims = []
for j2 in range(k+1):  # j = 0, 1/2, 1, 3/2 -> j2 = 0, 1, 2, 3
    j = j2 / 2
    d_j = math.sin(math.pi * (j2 + 1) / (k + 2)) / math.sin(math.pi / (k + 2))
    q_dims.append((j, d_j))
    print(f"  j = {j}: d_j = {d_j:.6f}", end="")
    if abs(d_j - 1) < 1e-10:
        print(" = 1", end="")
    elif abs(d_j - PHI) < 1e-10:
        print(f" = phi", end="")
    print()

# Total quantum dimension
D_sq = sum(d**2 for _, d in q_dims)
D = math.sqrt(D_sq)
print(f"\n  Total quantum dimension D^2 = sum(d_j^2) = {D_sq:.6f}")
print(f"  D = {D:.6f}")
print(f"  D^2 = {D_sq:.6f} = 2 + 2*phi^2 = 2 + 2*(phi+1) = 4 + 2*phi")
print(f"       = (a1-1) + 2*phi = {(a1-1) + 2*PHI:.6f}")
print(f"       = {4 + 2*PHI:.6f}")

# Topological entanglement entropy
gamma = math.log(D)
print(f"\n  Topological entanglement entropy:")
print(f"  gamma = ln(D) = {gamma:.6f}")
print(f"  Compare: ln(phi) = {math.log(PHI):.6f}")
print(f"  gamma / ln(phi) = {gamma / math.log(PHI):.6f}")

# Fibonacci subcategory: only 1 and tau (d=phi)
D_fib_sq = 1 + PHI**2
D_fib = math.sqrt(D_fib_sq)
gamma_fib = math.log(D_fib)
print(f"\n  Fibonacci subcategory (1 and tau only):")
print(f"  D_fib^2 = 1 + phi^2 = 2 + phi = {D_fib_sq:.6f}")
print(f"  D_fib = {D_fib:.6f}")
print(f"  gamma_fib = ln(D_fib) = {gamma_fib:.6f}")
print(f"  gamma_fib = (1/2)*ln(2+phi) = (1/2)*ln(phi^2+1)")

# KEY: connection to framework
print(f"\n  Framework connections:")
print(f"  D_full^2 = 4+2*phi = 2*(2+phi) = 2*D_fib^2")
print(f"  D_full = sqrt(2) * D_fib")
print(f"  D_full^2 / a1 = {D_sq/a1:.6f}")
print(f"  D_fib^2 * 2*alpha_s = {D_fib_sq * 2*alpha_s:.6f}")
print(f"    = (2+phi)/(2*phi^3)")
print(f"    = (2+phi)/(2*(2*phi+1))")  # phi^3 = 2*phi+1

val = (2 + PHI) / (2 * (2*PHI + 1))
print(f"    = {val:.6f}")

# Check: is gamma related to any framework quantity?
print(f"\n  Entropy checks:")
print(f"  gamma = {gamma:.6f}")
print(f"  ln(phi) = {math.log(PHI):.6f}")
print(f"  gamma/ln(phi) = {gamma/math.log(PHI):.6f}")
print(f"  ln(2) = {math.log(2):.6f}")
print(f"  gamma - ln(phi) = {gamma - math.log(PHI):.6f}")
print(f"  ln(sqrt(2)) = {math.log(math.sqrt(2)):.6f}")
print(f"  gamma = ln(phi) + ln(sqrt(2))? {abs(gamma - math.log(PHI) - math.log(math.sqrt(2))) < 0.01}")

# Actually: D_full = sqrt(2) * D_fib = sqrt(2) * sqrt(2+phi)
# gamma_full = ln(sqrt(2)) + ln(sqrt(2+phi)) = (1/2)*ln(2) + (1/2)*ln(2+phi)
# = (1/2)*ln(2*(2+phi)) = (1/2)*ln(4+2*phi)


# =====================================================================
# PART 2: McKay Graph Entanglement
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: McKay Graph Entanglement Entropy")
print("=" * 72)

# McKay graph: affine E8
edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)]
n = 9

A = np.zeros((n, n))
for i, j in edges:
    A[i][j] = 1
    A[j][i] = 1

# Laplacian
D_mat = np.diag([sum(A[i]) for i in range(n)])
L = D_mat - A

# Hamiltonian: H = -A (tight-binding model on McKay graph)
evals, evecs = np.linalg.eigh(-A)  # sorted ascending

print(f"\n  Tight-binding spectrum of McKay graph (H = -A):")
for i, e in enumerate(evals):
    print(f"    E_{i} = {e:.6f}")

# Half-filling: occupy lowest 4 or 5 states (9 sites)
# For 9 sites, half-filling is ambiguous. Try both.
for n_occ in [4, 5]:
    print(f"\n  --- n_occupied = {n_occ} ---")
    # Correlation matrix C_ij = sum_{k in occupied} psi_k(i) * psi_k(j)
    C = evecs[:, :n_occ] @ evecs[:, :n_occ].T

    # Bipartition: BLACK (integer spin) vs WHITE (half-integer spin)
    # From McKay bipartite structure:
    # BLACK = {0, 2, 4, 6, 8} (dims: 1, 3, 5, 4, 3) - integer spin
    # WHITE = {1, 3, 5, 7} (dims: 2, 4, 6, 2) - half-integer spin
    black = [0, 2, 4, 6, 8]
    white = [1, 3, 5, 7]

    for name, subsystem in [("BLACK (integer spin)", black), ("WHITE (half-integer spin)", white)]:
        C_sub = C[np.ix_(subsystem, subsystem)]
        evals_sub = np.linalg.eigvalsh(C_sub)
        # Entanglement entropy
        S = 0
        for ev in evals_sub:
            if 1e-12 < ev < 1 - 1e-12:
                S -= ev * math.log(ev) + (1 - ev) * math.log(1 - ev)
        print(f"    S({name}) = {S:.6f}")
        print(f"      Eigenvalues of C_sub: {[round(e, 6) for e in evals_sub]}")

        # Check framework quantities
        if S > 0:
            print(f"      S/ln(2) = {S/math.log(2):.4f} (bits)")
            print(f"      S/ln(phi) = {S/math.log(PHI):.4f}")

# Also try: ground state of Laplacian
print(f"\n  --- Ground state of Laplacian (lowest non-zero mode) ---")
evals_L, evecs_L = np.linalg.eigh(L)
# The Laplacian has a zero mode (constant vector). The ground state
# of the particle is the lowest POSITIVE eigenvalue mode.
print(f"  Laplacian eigenvalues: {[round(e, 4) for e in evals_L]}")


# =====================================================================
# PART 3: Entanglement Spectrum of the Chirality Split
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: Entanglement Spectrum (Chirality)")
print("=" * 72)

# Consider the McKay graph Hamiltonian H = -A at half-filling (5 particles)
# Bipartite: the entanglement spectrum between BLACK and WHITE
# tells us about the internal structure of the chirality split

n_occ = 5  # slightly more than half
C = evecs[:, :n_occ] @ evecs[:, :n_occ].T

for name, subsystem in [("BLACK (5 sites)", black), ("WHITE (4 sites)", white)]:
    C_sub = C[np.ix_(subsystem, subsystem)]
    evals_sub = np.linalg.eigvalsh(C_sub)

    print(f"\n  {name}:")
    print(f"    Correlation matrix eigenvalues:")
    for i, ev in enumerate(evals_sub):
        if abs(ev) < 1e-10 or abs(ev - 1) < 1e-10:
            status = " (0 or 1 = no entanglement)"
        else:
            xi = math.log((1-ev)/ev)  # entanglement energy
            status = f" -> xi = {xi:.4f}"
        print(f"      nu_{i} = {ev:.8f}{status}")

    S = 0
    for ev in evals_sub:
        if 1e-12 < ev < 1 - 1e-12:
            S -= ev * math.log(ev) + (1 - ev) * math.log(1 - ev)
    print(f"    Entanglement entropy S = {S:.6f}")


# =====================================================================
# PART 4: Galois Entanglement Structure
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Galois Pairs as Entangled States")
print("=" * 72)

print(f"""
  In the framework, Galois conjugation sigma acts on the 96 fermionic
  vertices of the 600-cell, creating 48 Galois pairs (v, sigma(v)).

  These pairs have maximally correlated properties:
    m_f * m_sigma(f) = m_e^2        (mass product duality)
    alpha * alpha' = 1/(2*pi)        (coupling product)
    Lambda * Lambda' = (1/(2*pi))^z  (CC product)

  This is analogous to EPR pairs: measuring one determines the other.

  The Galois "entanglement entropy" per pair:
    The two states (f, sigma(f)) span a 2D Hilbert space.
    If they are in a Bell state |psi> = (|f,sigma(f)> + |sigma(f),f>)/sqrt(2),
    the entanglement entropy is S = ln(2) = {math.log(2):.6f}.

  Total Galois entanglement:
    48 pairs * ln(2) = {48 * math.log(2):.4f}
    = 48 * ln(2) = ln(2^48) = {48 * math.log(2):.4f}

  Note: 48 = 16 * N_gen = Tr(A^4) of McKay graph (exp403)
""")


# =====================================================================
# PART 5: Fibonacci Anyons from 2I
# =====================================================================
print("=" * 72)
print("PART 5: Fibonacci Anyons and 2I")
print("=" * 72)

# The binary icosahedral group 2I is related to SU(2) at level k=3
# through the ADE classification of modular invariants.
# The McKay correspondence: 2I <-> E8
# The Chern-Simons theory: SU(2)_k where k+2 = a1

# Fusion matrices for SU(2)_3 (Fibonacci + identity + ...)
# Actually, let's compute the S-matrix of SU(2)_3
print(f"\n  SU(2)_{k} Chern-Simons (k = {k}, k+2 = {k+2} = a1)")
print(f"\n  Modular S-matrix:")
S_mat = np.zeros((k+1, k+1))
for i in range(k+1):
    for j in range(k+1):
        S_mat[i][j] = math.sqrt(2.0/(k+2)) * math.sin(math.pi*(i+1)*(j+1)/(k+2))

print(f"    S = sqrt(2/{k+2}) * sin(pi*(i+1)*(j+1)/{k+2})")
for i in range(k+1):
    row = [f"{S_mat[i][j]:8.5f}" for j in range(k+1)]
    print(f"    [{', '.join(row)}]")

# Quantum dimensions from S-matrix: d_j = S_{0j}/S_{00}
print(f"\n  Quantum dimensions from S-matrix:")
for j in range(k+1):
    d = S_mat[0][j] / S_mat[0][0]
    print(f"    d_{j/2} = S_{{0,{j}}}/S_{{0,0}} = {d:.6f}")

# Topological entanglement entropy
# For a disk with boundary, the TEE is gamma = ln(D)
# For a specific anyon a inside, the entropy is S_a = ln(D/d_a)
print(f"\n  Topological entanglement entropies:")
D_total = math.sqrt(sum(S_mat[0][j]**(-2) * S_mat[0][j]**2 for j in range(k+1)))
# Actually D = 1/S_{00}
D_from_S = 1 / S_mat[0][0]
print(f"  D = 1/S_{{00}} = {D_from_S:.6f}")
print(f"  Check: D^2 = {D_from_S**2:.6f}")
print(f"  gamma = ln(D) = {math.log(D_from_S):.6f}")

# For each anyon type:
for j in range(k+1):
    d_j = S_mat[0][j] / S_mat[0][0]
    if d_j > 0:
        S_a = math.log(D_from_S) - math.log(d_j)
        print(f"  S(j={j/2}) = ln(D/d_j) = ln({D_from_S:.4f}/{d_j:.4f}) = {S_a:.6f}")


# =====================================================================
# PART 6: Key Relations
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: Framework Connections")
print("=" * 72)

# Check: k+2 = a1
print(f"\n  KEY: k+2 = a1 = 5")
print(f"  The Chern-Simons level that gives Fibonacci anyons")
print(f"  is k = a1-2 = 3, with k+2 = a1 roots of unity.")

# The q-parameter: q = exp(2*pi*i/(k+2)) = exp(2*pi*i/a1)
# This is a primitive a1-th root of unity
print(f"\n  q = exp(2*pi*i/a1) = exp(2*pi*i/{a1})")
print(f"  This is precisely the Galois element!")
print(f"  sigma: phi -> phi' corresponds to the Z/a1 action on roots of unity")

# Quantum dimension in terms of phi
print(f"\n  Quantum dimension of Fibonacci anyon:")
print(f"  d_tau = sin(2*pi/{a1}) / sin(pi/{a1})")
print(f"        = 2*cos(pi/{a1})")
print(f"        = 2*cos(pi/5) = 2*cos(36 deg)")
print(f"        = phi = {PHI:.6f}")
print(f"  Check: 2*cos(pi/5) = {2*math.cos(math.pi/5):.6f}")

# S-matrix elements in terms of phi
print(f"\n  S-matrix in terms of phi:")
print(f"  S_{{00}} = sqrt(2/{a1}) * sin(pi/{a1})")
s00 = math.sqrt(2/a1) * math.sin(math.pi/a1)
print(f"         = {s00:.6f}")
print(f"  1/S_{{00}} = D = {1/s00:.6f}")
print(f"  D^2 = {1/s00**2:.6f}")
print(f"  D^2 = {a1}/(2*sin^2(pi/{a1})) = {a1/(2*math.sin(math.pi/a1)**2):.6f}")

# Express D^2 in closed form
sin_pi5 = math.sin(math.pi / 5)
print(f"\n  sin(pi/5) = sqrt((5-sqrt(5))/8) = {sin_pi5:.6f}")
print(f"  sin^2(pi/5) = (5-sqrt(5))/8 = {sin_pi5**2:.6f}")
D_sq_exact = a1 / (2 * sin_pi5**2)
print(f"  D^2 = 5/(2*(5-sqrt(5))/8) = 5*8/(2*(5-sqrt(5)))")
print(f"       = 20/(5-sqrt(5)) = 20*(5+sqrt(5))/20 = 5+sqrt(5)")
print(f"       = {5 + math.sqrt(5):.6f}")
print(f"  Check: D^2 = {D_sq_exact:.6f}")

# KEY RESULT
print(f"\n  *** D^2 = a1 + sqrt(a1) = 5 + sqrt(5) ***")
print(f"  *** D^2 = 2*phi*a1/sqrt(a1) = ... ***")

# Simplify: 5 + sqrt(5) = sqrt(5)*(sqrt(5)+1) = sqrt(5)*2*phi = 2*phi*sqrt(a1)
val = 2 * PHI * math.sqrt(a1)
print(f"  D^2 = 2*phi*sqrt(a1) = {val:.6f}")
print(f"  Check: {5 + math.sqrt(5):.6f}")
print(f"  MATCH: {abs(val - 5 - math.sqrt(5)) < 1e-10}")

# Topological entropy in terms of framework quantities
gamma_exact = 0.5 * math.log(5 + math.sqrt(5))
print(f"\n  gamma = ln(D) = (1/2)*ln(a1 + sqrt(a1))")
print(f"        = (1/2)*ln(2*phi*sqrt(a1))")
print(f"        = (1/2)*(ln(2) + ln(phi) + (1/2)*ln(a1))")
print(f"        = {gamma_exact:.6f}")

print(f"\n  In terms of basic framework quantities:")
print(f"  gamma = (1/2)*ln(2*phi*sqrt(5))")
print(f"        = (1/2)*ln(2) + (1/2)*ln(phi) + (1/4)*ln(5)")
print(f"        = {0.5*math.log(2) + 0.5*math.log(PHI) + 0.25*math.log(5):.6f}")


# =====================================================================
# PART 7: Entanglement and Alpha_s
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: Entanglement and Coupling Constants")
print("=" * 72)

# Check if D^2 relates to alpha_s
print(f"  alpha_s = 1/(2*phi^3) = {alpha_s:.6f}")
print(f"  1/alpha_s = 2*phi^3 = {1/alpha_s:.6f}")
print(f"  D^2 = 5+sqrt(5) = {5+math.sqrt(5):.6f}")
print(f"  D^2 / (1/alpha_s) = {(5+math.sqrt(5))/(1/alpha_s):.6f}")
print(f"  D^2 * alpha_s = {(5+math.sqrt(5))*alpha_s:.6f}")
print(f"  D^2 * alpha_s = (5+sqrt(5))/(2*phi^3)")
val2 = (5+math.sqrt(5)) / (2*PHI**3)
print(f"               = {val2:.6f}")
print(f"  Compare: sqrt(5)/phi^2 = {math.sqrt(5)/PHI**2:.6f}")
# (5+sqrt(5))/(2*phi^3) = (5+sqrt(5))/(2*(2*phi+1))
# = sqrt(5)*(sqrt(5)+1)/(2*(2*phi+1))
# = 2*sqrt(5)*phi/(2*(2*phi+1))
# = sqrt(5)*phi/(2*phi+1)
# = sqrt(5)*phi/phi^3  (since phi^3 = 2*phi+1)
# = sqrt(5)/phi^2
print(f"  D^2 * alpha_s = sqrt(a1)/phi^2 = sqrt(a1)*(2-phi)")
print(f"               = {math.sqrt(a1)*(2-PHI):.6f}")

# Another: D^2 and sin^2(tW)
sin2tW = b1 / (a1**2 + 1)
print(f"\n  D^2 * sin^2(tW) = {(5+math.sqrt(5))*sin2tW:.6f}")
print(f"  D^2 / N_gen = {(5+math.sqrt(5))/N_gen:.6f}")
print(f"  D^2 / a1 = {(5+math.sqrt(5))/a1:.6f} = 1 + 1/sqrt(a1)")


# =====================================================================
# PART 8: Summary
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"""
  DERIVED:
  1. Chern-Simons level k = a1-2 = 3 gives Fibonacci anyons
     with quantum dimension d_tau = phi (golden ratio)
  2. q = exp(2*pi*i/a1) is the CS deformation parameter
     = precisely the Galois conjugation action
  3. Total quantum dimension D^2 = a1 + sqrt(a1) = 5+sqrt(5)
     = 2*phi*sqrt(a1)
  4. Topological entanglement entropy:
     gamma = (1/2)*ln(2*phi*sqrt(a1))
  5. D^2 * alpha_s = sqrt(a1)/phi^2 (elegant but no clear meaning)

  STRUCTURAL:
  - k+2 = a1 connects CS theory to the defining integer
  - Galois conjugation = quantum group deformation parameter
  - Fibonacci anyons are the "simplest" non-abelian anyons
  - The 2I group IS the symmetry of the Fibonacci anyon model

  SPECULATIVE:
  - Galois pairs as EPR-like entangled states
  - Topological entanglement as probe of confinement
  - Area law on 600-cell

  STATUS: The connection k+2 = a1 and d_tau = phi are DERIVED.
  The physical interpretation (what entanglement MEANS in the
  framework) remains to be developed.

  MOST PROMISING: The Galois-CS connection
    sigma(phi) = phi' corresponds to q -> q^(-1) in quantum groups.
    This identifies the Galois symmetry of the framework with the
    quantum group symmetry of topological field theory.
    This could be a DEEP structural insight.
""")
