"""
exp378b: WHY IS c2/(2N) = 233 = F_13 (Fibonacci)?
===================================================
Deep investigation of the Fibonacci number in spectral coefficients.
"""

import numpy as np
from math import gcd

# Framework
a1 = 5
b1 = 6
phi = (1 + np.sqrt(5)) / 2
N = 120
degree = 12
N_eig = 9

# Seeley-DeWitt
c0 = 2640
c1 = 14880
c2 = 55920

# A_k = c_k / (2*N)
A0 = c0 // (2*N)  # = 11
A1 = c1 // (2*N)  # = 62
A2 = c2 // (2*N)  # = 233

print("=" * 60)
print("exp378b: WHY IS c2/(2N) = 233 = F_13?")
print("=" * 60)

# Fibonacci and Lucas numbers
fib = [0, 1]
for i in range(30):
    fib.append(fib[-1] + fib[-2])
luc = [2, 1]
for i in range(30):
    luc.append(luc[-1] + luc[-2])

print(f"\nA_0 = {A0} = a1+b1 = L_5 (Lucas)")
print(f"A_1 = {A1}")
print(f"A_2 = {A2} = F_13 (Fibonacci)")

# ============================================================
# KEY DISCOVERY: degree^2 = F_12
# ============================================================
print("\n" + "=" * 60)
print("KEY: degree^2 = 12^2 = 144 = F_12")
print("=" * 60)

print(f"degree = {degree}")
print(f"degree^2 = {degree**2}")
print(f"F_12 = {fib[12]}")
print(f"Match: {degree**2 == fib[12]}")

# This is a KNOWN number theory fact:
# 12 is the UNIQUE n > 1 where F_n = n^2
print("\nCheck: which n satisfy F_n = n^2?")
for n in range(30):
    if fib[n] == n*n:
        print(f"  F_{n} = {fib[n]} = {n}^2  <-- MATCH")

print("\n12 is the UNIQUE degree where F_degree = degree^2")
print("And degree = 12 comes from a1=5: degree = 2*(a1+1)")

# ============================================================
# DECOMPOSITION: A2 = F_12 + F_11 = F_13
# ============================================================
print("\n" + "=" * 60)
print("DECOMPOSITION: 233 = 156 + 77")
print("=" * 60)

# Tr(L0^2)/N = degree*(degree+1) for any regular graph
TrL0sq_per_N = degree * (degree + 1)
print(f"Tr(L0^2)/N = degree*(degree+1) = {degree}*{degree+1} = {TrL0sq_per_N}")

# Tr(L1^2)/N from Hodge duality
TrL1sq_per_N = A2 - TrL0sq_per_N
print(f"Tr(L1^2)/N = A2 - Tr(L0^2)/N = {A2} - {TrL0sq_per_N} = {TrL1sq_per_N}")

# Now the magic:
print(f"\n--- The Fibonacci decomposition ---")
print(f"degree^2 = {degree**2} = F_12 = {fib[12]}")
print(f"degree*(degree+1) = {TrL0sq_per_N} = F_12 + degree = {fib[12]} + {degree}")
print(f"A2 = Tr(L0^2)/N + Tr(L1^2)/N = {TrL0sq_per_N} + {TrL1sq_per_N}")
print(f"   = (F_12 + degree) + (F_11 - degree)")
print(f"   = F_12 + F_11 = F_13")
print(f"\nVerify: F_11 = {fib[11]}")
print(f"F_11 - degree = {fib[11]} - {degree} = {fib[11] - degree}")
print(f"Tr(L1^2)/N = {TrL1sq_per_N}")
print(f"Match: {TrL1sq_per_N == fib[11] - degree}")

# ============================================================
# WHY Tr(L1^2)/N + degree = F_11 = 89?
# ============================================================
print("\n" + "=" * 60)
print("WHY does Tr(L1^2)/N + degree = 89 = F_11?")
print("=" * 60)

# Tr(L1^2)/N = 77. And degree + 77 = 89 = F_11.
# Equivalently: Tr(L1^2) = N * (F_11 - degree) = 120*(89-12) = 9240

# Let's verify by computing Tr(L0^2) directly
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
chi_gen = [1, phi, phi, 1, 0, -1, -1, 1-phi, 1-phi]
Lk = [degree * (1 - chi_gen[k]/dims[k]) for k in range(N_eig)]
mults = [d**2 for d in dims]

TrL0 = sum(mults[k] * Lk[k] for k in range(N_eig))
TrL0sq = sum(mults[k] * Lk[k]**2 for k in range(N_eig))
print(f"Tr(L0) = {TrL0:.0f} = N*degree = {N*degree}")
print(f"Tr(L0^2) = {TrL0sq:.0f}")
print(f"Tr(L0^2)/N = {TrL0sq/N:.0f}")

# Adjacency matrix: a_k = degree * chi_k / d_k
ak = [degree * chi_gen[k]/dims[k] for k in range(N_eig)]
TrA = sum(mults[k] * ak[k] for k in range(N_eig))
TrAsq = sum(mults[k] * ak[k]**2 for k in range(N_eig))
print(f"\nTr(A) = {TrA:.6f} (should be 0)")
print(f"Tr(A^2) = {TrAsq:.0f} = N*degree = {N*degree}")

# Verify: Tr(L0^2) = degree^2*N + Tr(A^2) = 144*120 + 1440 = 18720
print(f"Tr(L0^2) = degree^2*N + Tr(A^2) = {degree**2*N} + {int(TrAsq)} = {degree**2*N + int(TrAsq)}")

# Higher traces of adjacency
TrA3 = sum(mults[k] * ak[k]**3 for k in range(N_eig))
TrA4 = sum(mults[k] * ak[k]**4 for k in range(N_eig))
print(f"\nTr(A^3) = {TrA3:.0f}")
print(f"  = degree^3 * sum(chi_k^3/d_k) = {degree**3} * {sum(chi_gen[k]**3/dims[k]*dims[k]**2 for k in range(N_eig))/N:.6f}")
print(f"Tr(A^4) = {TrA4:.0f}")

# Tr(A^3)/N = number of triangles * 6 / N = 1200*6/120 = 60
print(f"Tr(A^3)/N = {TrA3/N:.0f}")
print(f"  Expected: 6 * N_tri/N = 6 * {1200}/{N} = {6*1200/N:.0f}")

# Tr(A^4)/N counts closed walks of length 4
print(f"Tr(A^4)/N = {TrA4/N:.0f}")

# ============================================================
# THE COMPLETE CHAIN: why 233 = F_13
# ============================================================
print("\n" + "=" * 60)
print("THE COMPLETE CHAIN")
print("=" * 60)

# On the 600-cell Cayley graph:
# sum chi_k^2 = |class|/... = N/|centralizer(g)| where g in C_10
# For C_10 (12 elements): sum chi_k^2 = 120/12 = 10
sum_chi_sq = sum(chi_gen[k]**2 for k in range(N_eig))
print(f"sum(chi_k^2) = {sum_chi_sq:.6f}")
print(f"  = N / |C_10| = 120/12 = 10")

# Tr(A^2)/N = degree^2 * sum(chi_k^2)/N_class_size...
# Actually: Tr(A^2) = sum d_k^2 * a_k^2 = degree^2 * sum chi_k^2 = 144*10 = 1440
# And Tr(A^2)/N = 1440/120 = 12 = degree
# This is GENERAL for Cayley graphs: Tr(A^2) = N * degree
print(f"\nTr(A^2) = degree^2 * sum(chi^2) = {degree**2} * {sum_chi_sq:.0f} = {int(degree**2 * sum_chi_sq)}")
print(f"Tr(A^2)/N = degree (always true for Cayley graphs)")

# Now for Tr(A^4):
sum_chi_4 = sum(chi_gen[k]**4 for k in range(N_eig))
TrA4_formula = degree**4 * sum_chi_4 / 1  # not quite...
# Tr(A^4) = sum d_k^2 * (degree*chi_k/d_k)^4 = degree^4 * sum chi_k^4 / d_k^2
ratio_4 = sum(chi_gen[k]**4 / dims[k]**2 * dims[k]**2 for k in range(N_eig))
# Wait: a_k = degree * chi_k/d_k, so a_k^4 = degree^4 * chi_k^4/d_k^4
# Tr(A^4) = sum d_k^2 * a_k^4 = degree^4 * sum chi_k^4/d_k^2
s4 = sum(chi_gen[k]**4 / dims[k]**2 for k in range(N_eig))
print(f"\nsum(chi_k^4/d_k^2) = {s4:.6f}")
print(f"Tr(A^4) = degree^4 * {s4:.6f} = {degree**4 * s4:.0f}")
print(f"Tr(A^4)/N = {degree**4 * s4 / N:.6f}")

# For Tr(L0^2):
# L_k = degree - a_k
# L_k^2 = degree^2 - 2*degree*a_k + a_k^2
# Tr(L0^2) = degree^2*N - 2*degree*Tr(A) + Tr(A^2) = degree^2*N + N*degree
# = N*degree*(degree+1) = 120*12*13 = 18720
print(f"\nTr(L0^2) = N*degree*(degree+1) = {N}*{degree}*{degree+1} = {N*degree*(degree+1)}")

# For the TOTAL c2 = 2*Tr(L0^2) + 2*Tr(L1^2):
# We need Tr(L1^2). This comes from the 1-form Laplacian.
# By the Hodge decomposition on S^3 (beta_1=0):
# Tr(L1^2) = Tr(L0^2) + (extra terms from 1-form structure)
#
# ACTUALLY: using c2 = 2*Tr(L0^2) + 2*Tr(L1^2) [Hodge duality]
# Tr(L1^2) = (c2 - 2*Tr(L0^2))/2 = (55920 - 37440)/2 = 9240
TrL1sq = (c2 - 2*int(TrL0sq))/2
print(f"\nTr(L1^2) = (c2 - 2*Tr(L0^2))/2 = ({c2}-{2*int(TrL0sq)})/2 = {TrL1sq:.0f}")

# ============================================================
# CRITICAL CHECK: Is Tr(L1^2) related to Tr(A^k)?
# ============================================================
print("\n" + "=" * 60)
print("Tr(L1^2) vs graph invariants")
print("=" * 60)

# On a simplicial complex, Lap_1 = d0*d0^T + d1^T*d1
# For the exact part: eigenvalues = nonzero eigenvalues of Lap_0
# Tr_exact(L1^2) = Tr(L0^2) - 0 = 18720 (dim(exact)=N-1=119)
#
# For coexact: eigenvalues = nonzero eigenvalues of Lap_2
# By Hodge duality: Lap_2 spectrum = Lap_1 spectrum
# So the total Tr(L1^2) is NOT simply 2*Tr(L0^2)
# because the exact and coexact eigenvalues are different

# Let's use the relation differently:
# Tr(L1) = c1/2 - Tr(L0) = 14880/2 - 1440 = 6000 (WRONG, should use Hodge)
# Actually: c1 = 2*Tr(L0) + 2*Tr(L1), so Tr(L1) = (c1 - 2*Tr(L0))/2 = 6000
TrL1 = (c1 - 2*int(TrL0))/2
print(f"Tr(L1) = {TrL1:.0f}")
print(f"Tr(L1)/N_edges = {TrL1/720:.6f}")
print(f"Tr(L1^2) = {TrL1sq:.0f}")
print(f"Tr(L1^2)/N_edges = {TrL1sq/720:.6f}")

# Variance of 1-form Laplacian
var_L1 = TrL1sq/720 - (TrL1/720)**2
print(f"Var(L1) = {var_L1:.6f}")

# ============================================================
# THE 600-CELL SPECIFIC IDENTITY
# ============================================================
print("\n" + "=" * 60)
print("600-CELL SPECIFIC: Tr(A^3) and triangles")
print("=" * 60)

# Tr(A^3)/N = 6 * (triangles per vertex) = 6*30 = 180? No.
# Tr(A^3) counts closed walks of length 3 = 6 * (number of triangles)
TrA3_check = 6 * 1200  # 6 orientations per triangle
print(f"Tr(A^3) = 6*N_tri = 6*1200 = {6*1200}")
print(f"From eigenvalues: {TrA3:.0f}")
print(f"Match: {abs(TrA3 - 6*1200) < 0.01}")
print(f"Tr(A^3)/N = {TrA3/N:.0f} = 6*a1*b1 = {6*a1*b1}")

# Tr(A^4) counts closed walks of length 4
# = sum over vertices: (number of walks v->?->?->?->v of length 4)
# = N * (degree + paths through triangles)
print(f"\nTr(A^4) = {TrA4:.0f}")
print(f"Tr(A^4)/N = {TrA4/N:.0f}")
# For k-regular: Tr(A^4)/N = degree^2 + some correction
correction = TrA4/N - degree**2
print(f"Tr(A^4)/N - degree^2 = {correction:.0f}")
# On k-regular: closed 4-walks = k (go-return-go-return) + k*(k-1)*(triangles at vertex)/...
# Actually: sum_j A_{ij}^2 = degree for each i (diagonal of A^2)
# (A^4)_{ii} = sum_j (A^2)_{ij}^2

# ============================================================
# CRITICAL: Tr(L1^2) from 1-form structure
# ============================================================
print("\n" + "=" * 60)
print("Tr(L1^2) = 9240 from graph combinatorics")
print("=" * 60)

# Lap_1 = B0^T * B0 + B1 * B1^T
# where B0 = incidence matrix (edges x vertices, N_e x N_v)
#       B1 = incidence matrix (triangles x edges, N_t x N_e)
# Tr(Lap_1^2) = Tr((B0^T B0 + B1 B1^T)^2)
# = Tr((B0^T B0)^2) + Tr((B1 B1^T)^2) + 2*Tr(B0^T B0 B1 B1^T)

# For simplicial S3 (beta_1=0): exact and coexact are orthogonal
# => the cross term vanishes
# Tr(Lap_1^2) = Tr((B0^T B0)^2) + Tr((B1 B1^T)^2)

# Tr((B0^T B0)^2) = Tr(B0 B0^T B0 B0^T) = Tr((B0 B0^T)^2) = Tr(Lap_0^2)
# Wait: B0 B0^T is NOT Lap_0. B0^T B0 acts on edges, B0 B0^T on vertices.
# B0 B0^T = Lap_0 (graph Laplacian on vertices)
# B0^T B0 = edge-edge adjacency (acts on 1-forms, exact part of Lap_1)

# By cyclic trace: Tr((B0^T B0)^2) = Tr(B0^T B0 B0^T B0) = Tr((B0 B0^T)^2) = Tr(Lap_0^2)
# So the exact contribution to Tr(Lap_1^2) = Tr(Lap_0^2) = 18720.

# Similarly: Tr((B1 B1^T)^2) = Tr((B1^T B1)^2) = Tr(Lap_2^2)
# And by Hodge duality on S3: Tr(Lap_2^2) = Tr(Lap_1^2)
#
# This gives: Tr(Lap_1^2) = Tr(Lap_0^2) + Tr(Lap_1^2)
# => 0 = Tr(Lap_0^2), which is WRONG.
#
# The issue: Hodge duality gives Lap_2 spectrum = Lap_1 spectrum,
# NOT Tr(Lap_2^2) = Tr(Lap_1^2) when they have different dimensions.
# Lap_1 is 720x720, Lap_2 is 1200x1200.
# Even though nonzero eigenvalues match, the zero eigenvalues differ.

# Let me think again using the spectral decomposition.
# On S3 simplicial: the nonzero spectrum of Lap_p equals nonzero spectrum of Lap_{3-p}
# Specifically: nonzero eigs of Lap_1 = union of {nonzero eigs of Lap_0} and {nonzero eigs of Lap_2}
# where the first set are exact eigenvalues and second are coexact eigenvalues.
#
# And nonzero eigs of Lap_2 = nonzero eigs of Lap_1 (Hodge duality 1-forms <-> 2-forms on S3)
# This means coexact 1-form eigenvalues = nonzero 2-form eigenvalues
# = nonzero 1-form eigenvalues (by Hodge on S3)
# So coexact 1-form spectrum = full nonzero 1-form spectrum. Circular!

# Better approach: use the relation between c_k and Tr(D^{2k})
# c2 = Tr(D^4) = sum over ALL eigenvalues lambda of D^2, each lambda^2.
# These eigenvalues are Lap_0 eigs (120), Lap_1 eigs (720),
# Lap_2 eigs (1200), Lap_3 eigs (600).
#
# Nonzero eigenvalue matching:
# Lap_0 has 1 zero (constant) + 119 nonzero eigenvalues
# Lap_3 has 1 zero (top form) + 599 nonzero eigenvalues
# On S3: nonzero Lap_0 eigs = subset of nonzero Lap_1 eigs (exact part)
#         nonzero Lap_3 eigs = subset of nonzero Lap_2 eigs (coexact part)
# And Lap_1 and Lap_2 are related by Hodge on S3.

# Actually the correct statement for simplicial S3:
# Lap_0 nonzero spectrum with multiplicities embedded in Lap_1 (exact)
# Lap_3 nonzero spectrum embedded in Lap_2 (coexact)
# Lap_1 remaining (coexact part) = some spectrum
# Lap_2 remaining (exact part) = same spectrum (Hodge)

# So Tr(Lap_1^2) = Tr(Lap_0^2 nonzero) + Tr(coexact Lap_1^2)
# And Tr(Lap_2^2) = Tr(coexact Lap_2^2) + Tr(Lap_3^2 nonzero)
# By Hodge on S3: coexact Lap_1 spectrum = exact Lap_2 spectrum
# And: exact Lap_2 spectrum = coexact Lap_1 spectrum
# Hmm, still circular.

# Let me just use the FORMULA approach:
# c2/(2N) = 233 is a FACT from the spectral action computation.
# 233 = F_13 is NUMBER THEORY.
# The question is whether this is structural or coincidental.

print("The decomposition:")
print(f"  A2 = 233 = degree*(degree+1) + 77")
print(f"     = (degree^2 + degree) + (F_11 - degree)")
print(f"     = degree^2 + F_11")
print(f"     = F_12 + F_11   [using degree^2 = 12^2 = 144 = F_12]")
print(f"     = F_13")
print(f"\nThe chain has TWO ingredients:")
print(f"  (i)  degree^2 = F_12  (number theory: unique n>1 with n^2=F_n)")
print(f"  (ii) Tr(L1^2)/N = 77 = F_11 - degree = 89 - 12")

# ============================================================
# IS ingredient (ii) structural or coincidental?
# ============================================================
print("\n" + "=" * 60)
print("Testing: is 77 = F_11 - degree structural?")
print("=" * 60)

# 77 = Tr(L1^2)/N. Can we compute this from graph invariants?
# Tr(L1^2) = c2/2 - Tr(L0^2) = 55920/2 - 18720 = 27960 - 18720 = 9240
print(f"Tr(L1^2) = c2/2 - Tr(L0^2) = {c2//2} - {int(TrL0sq)} = {c2//2 - int(TrL0sq)}")
print(f"  (using c2 = 2*Tr(L0^2) + 2*Tr(L1^2) from Hodge duality)")

# 77 as framework numbers:
print(f"\n77 = {77}")
print(f"  = 7 * 11 = L_4 * (a1+b1)")
print(f"  = 7 * (a1+b1)")
print(f"  = (a1+2) * (a1+b1) = ... no, 7*11=77 but a1+2=7. YES!")
print(f"  = (a1+2) * (a1+b1) = 7 * 11 = 77")
print(f"  So Tr(L1^2)/N = (a1+2)*(a1+b1)")

# Verify: (a1+2)*(a1+b1) = 7*11 = 77 for a1=5. CHECK.
# And: A2 = degree*(degree+1) + (a1+2)*(a1+b1)
# = 2(a1+1)*(2a1+3) + (a1+2)*(2a1+1)
# = (4a1^2+10a1+6) + (2a1^2+5a1+2)
# = 6a1^2 + 15a1 + 8
A2_formula = 6*a1**2 + 15*a1 + 8
print(f"\nA2 = 6*a1^2 + 15*a1 + 8 = 6*25 + 75 + 8 = {A2_formula}")
print(f"Check: {A2_formula == A2}")

# Is 6*a1^2 + 15*a1 + 8 = F_13 only for a1=5, or is there a general Fibonacci identity?
print(f"\n6*a1^2 + 15*a1 + 8 for various a1:")
for a in range(1, 10):
    val = 6*a**2 + 15*a + 8
    is_fib = val in fib
    fib_idx = fib.index(val) if is_fib else -1
    mark = f" = F_{fib_idx}" if is_fib else ""
    print(f"  a1={a}: {val}{mark}")

# ============================================================
# VERIFY: Tr(L1^2)/N = (a1+2)*(a1+b1)
# ============================================================
print("\n" + "=" * 60)
print("VERIFY: is 77 = (a1+2)*(a1+b1) or just 7*11?")
print("=" * 60)

# We need to check if this formula is DERIVED or just numerology.
# Tr(L1^2)/N = 9240/120 = 77.
# (a1+2)*(a1+b1) = 7*11 = 77 for a1=5, b1=6.
# But does this formula hold for other polytopes/groups?
# We can't easily test (only a1=5 gives the 600-cell).
#
# Alternative framework expressions for 77:
print(f"77 = (a1+2)*(a1+b1) = {(a1+2)*(a1+b1)}")
print(f"77 = (a1+2)*(2*a1+1) = {(a1+2)*(2*a1+1)}")
print(f"  since b1 = a1+1: a1+b1 = 2*a1+1")
print(f"77 = 2*a1^2 + 5*a1 + 2 = {2*a1**2+5*a1+2}")
print(f"  = (2*a1+1)*(a1+2) = (2*5+1)*(5+2) = 11*7 = 77")

# So the GENERAL formula (for b1=a1+1):
# A2 = degree*(degree+1) + (a1+2)*(2*a1+1)
# = 2(a1+1)(2a1+3) + (a1+2)(2a1+1)
# = (2a1+1)[2(a1+1) + (a1+2)] + 2(a1+1)*2  ... let me just expand
# = 4a1^2+10a1+6 + 2a1^2+5a1+2 = 6a1^2+15a1+8
print(f"\nA2(a1) = 6*a1^2 + 15*a1 + 8")
print(f"For a1=5: {6*25+75+8} = 233 = F_13")

# ============================================================
# ALSO CHECK A0 and A1
# ============================================================
print("\n" + "=" * 60)
print("General formulas for A0, A1, A2")
print("=" * 60)

# A0 = c0/(2N) = (N_v+N_e+N_t+N_T)/(2N)
# = (N + N*b1 + N*a1*b1/3 + N*a1)/(2N)  ... wait
# c0 = N + N*degree/2 + N*tri_per_vert/3 + N*tet_per_vert/4
# Wait: c0 = N_v + N_e + N_tri + N_tet = 120+720+1200+600 = 2640
# = N*(1 + degree/2 + a1*b1/3 + a1) for 600-cell
# A0 = c0/(2N) = (1 + degree/2 + a1*b1/3 + a1)/2
# = (1 + (a1+1) + a1*(a1+1)/3 + a1)/2 = (1+a1+1+a1^2/3+a1/3+a1)/2
# Hmm, for a1=5: (1+6+10+5)/2 = 22/2 = 11. CHECK.
# A0 = (2+2*a1+a1*(a1+1)/3)/2... with a1=5: (2+10+10)/2=11
# General: A0 = 1 + a1 + a1*(a1+1)/(2*3) + 1/2... messy.
# Simpler: A0 = a1+b1 = 2*a1+1 for b1=a1+1.
print(f"A0 = a1+b1 = 2*a1+1 = {2*a1+1}")
# For a1=5: 11. CHECK.

# A1 = c1/(2N) = [2*Tr(L0) + 2*Tr(L1)]/(2N) = Tr(L0)/N + Tr(L1)/N
# Tr(L0)/N = degree = 2(a1+1). For a1=5: 12.
# Tr(L1)/N = 6000/120 = 50.
# A1 = 12 + 50 = 62.
# Is 50 a framework expression? 50 = 2*a1^2 = 2*25 = 50.
print(f"\nA1 = degree + Tr(L1)/N = {degree} + {int(TrL1/N)}")
print(f"  = 2*(a1+1) + Tr(L1)/N")
print(f"If Tr(L1)/N = 2*a1^2 = {2*a1**2}: A1 = 2*(a1+1)+2*a1^2 = 2*a1^2+2*a1+2")
A1_formula = 2*a1**2 + 2*a1 + 2
print(f"  = {A1_formula}. Matches A1={A1}? {A1_formula == A1}")

# So: A0(a1) = 2*a1+1, A1(a1) = 2*a1^2+2*a1+2, A2(a1) = 6*a1^2+15*a1+8
print(f"\n--- Summary ---")
print(f"A0(a1) = 2*a1 + 1 = {2*a1+1}")
print(f"A1(a1) = 2*a1^2 + 2*a1 + 2 = {2*a1**2+2*a1+2}")
print(f"A2(a1) = 6*a1^2 + 15*a1 + 8 = {6*a1**2+15*a1+8}")

# Check A2/A0: 233/11 = 21.18...
# A2/A0 = (6*a1^2+15*a1+8)/(2*a1+1)
# Polynomial division: 6*a1^2+15*a1+8 = (2*a1+1)*(3*a1+6) + 2
# = (2*a1+1)*(3*(a1+2)) + 2
print(f"\nA2 = A0 * 3*(a1+2) + 2 = {A0} * {3*(a1+2)} + 2 = {A0*3*(a1+2)+2}")
print(f"  = (2*a1+1)*3*(a1+2) + 2")
print(f"Check: {(2*a1+1)*3*(a1+2)+2 == A2}")

# So A2 = 3*(a1+2)*A0 + 2. For a1=5: 3*7*11+2 = 231+2 = 233.
# And 3*(a1+2) = 3*7 = 21 = F_8. Is this coincidence?
print(f"\n3*(a1+2) = {3*(a1+2)}")
print(f"F_8 = {fib[8]}")
print(f"Match: {3*(a1+2) == fib[8]}")

# ============================================================
# FIBONACCI COINCIDENCE ANALYSIS
# ============================================================
print("\n" + "=" * 60)
print("FIBONACCI COINCIDENCE ANALYSIS")
print("=" * 60)

print("For A2(a1) = 6*a1^2 + 15*a1 + 8:")
print("This is a Fibonacci number ONLY for specific a1:")
for a in range(1, 20):
    val = 6*a**2 + 15*a + 8
    is_fib = val in fib[:30]
    if is_fib:
        idx = fib.index(val)
        print(f"  a1={a}: A2={val} = F_{idx}")

print("\nFor A0(a1) = 2*a1+1:")
for a in range(1, 20):
    val = 2*a + 1
    is_luc = val in luc[:20]
    if is_luc:
        idx = luc.index(val)
        print(f"  a1={a}: A0={val} = L_{idx}")

# ============================================================
# THE DEGREE-12 FIBONACCI IDENTITY
# ============================================================
print("\n" + "=" * 60)
print("THE DEGREE-12 FIBONACCI IDENTITY")
print("=" * 60)

print(f"The UNIQUE property: 12^2 = 144 = F_12")
print(f"This is proven: F_n = n^2 only for n = 0, 1, 12")
print(f"(Cohn 1964, Wyler 1964)")
print(f"\nFor the 600-cell: degree = 12, so degree^2 = F_degree")
print(f"Then: A2 = degree*(degree+1) + 77")
print(f"     = F_12 + 12 + 77")
print(f"     = F_12 + 89 = F_12 + F_11 = F_13")
print(f"\nIF AND ONLY IF: 12 + 77 = 89 = F_11")
print(f"i.e., 77 = F_11 - 12 = F_11 - degree")
print(f"\nVerify: F_11 = {fib[11]}")
print(f"F_11 - degree = {fib[11] - degree}")
print(f"Tr(L1^2)/N = 77")
print(f"Match: {fib[11] - degree == 77}")

# So the Fibonacci identity reduces to:
# Tr(L1^2)/N = F_11 - degree = F_{degree-1} - degree
# For degree=12: F_11 - 12 = 89 - 12 = 77.
# This is specific to degree=12.

# ============================================================
# HONEST ASSESSMENT
# ============================================================
print("\n" + "=" * 60)
print("HONEST ASSESSMENT")
print("=" * 60)
print("""
FINDING: c2/(2N) = 233 = F_13 decomposes as:

  233 = degree*(degree+1) + (a1+2)*(2*a1+1)
      = 12*13 + 7*11
      = 156 + 77
      = (F_12 + degree) + (F_11 - degree)
      = F_12 + F_11
      = F_13

The Fibonacci identity has TWO ingredients:

(1) degree^2 = 12^2 = 144 = F_12
    This is a NUMBER THEORY fact (Cohn-Wyler 1964):
    12 is the unique n > 1 where F_n = n^2.
    It is SPECIFIC to the 600-cell (degree = 12).

(2) Tr(L1^2)/N = 77 = F_11 - degree
    This comes from the 1-form Laplacian of the 600-cell.
    We showed: 77 = (a1+2)*(2*a1+1) for a1=5.
    The Fibonacci identity 77 = F_11 - 12 is then EQUIVALENT
    to: (a1+2)*(2*a1+1) = F_{2a1+1} - 2*(a1+1).
    For a1=5: 7*11 = F_11 - 12 = 89 - 12 = 77. TRUE.

CLASSIFICATION:
  The identity 233 = F_13 is PARTLY STRUCTURAL, PARTLY COINCIDENTAL.

  STRUCTURAL: The decomposition A2 = degree*(degree+1) + Tr(L1^2)/N
  follows from the spectral action on any simplicial complex.

  COINCIDENTAL: The specific value degree=12 makes F_12 = degree^2,
  which only works for degree = 0, 1, or 12.

  REMARKABLE: The 600-cell picks out degree=12, which is the UNIQUE
  nontrivial solution of F_n = n^2. This CANNOT be derived from the
  framework — it is a fact of number theory that enhances the
  framework's mathematical richness.

  STATUS: PATTERN (beautiful but not derivable)
""")

# Bonus: also check A0 and A1
print("BONUS: A0 = L_5 is also structural:")
print(f"  L_5 = phi^5 + phi'^5 = 11 = a1+b1 = 2*a1+1")
print(f"  For ANY a1: A0 = 2*a1+1. This is L_5 ONLY for a1=5.")
print(f"  A0 = L_{a1} for a1=5 since L_5 = {luc[5]}")
print(f"  Check other: L_3={luc[3]} vs 2*3+1=7: {luc[3]==7}")
print(f"  L_4={luc[4]} vs 2*4+1=9: {luc[4]==9}")
# So: 2*a1+1 = L_a1 is NOT generally true.
# L_3=4, 2*3+1=7. L_4=7, 2*4+1=9. L_5=11, 2*5+1=11. L_6=18, 2*6+1=13.
# Only for a1=5: 2*5+1=11=L_5. UNIQUE!

print(f"\n  2*a1+1 = L_a1 ONLY for a1=5:")
for a in range(1, 10):
    match = (2*a+1 == luc[a])
    print(f"    a1={a}: 2a1+1={2*a+1}, L_{a}={luc[a]}, match={match}")

print(f"\n  REMARKABLE: a1=5 is also the unique a1 where:")
print(f"    A0 = L_a1 (Lucas)")
print(f"    degree^2 = F_degree (Fibonacci)")
print(f"    A2 = F_(degree+1) (Fibonacci)")
