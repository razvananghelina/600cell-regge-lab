#!/usr/bin/env python3
"""
EXP-209: MECHANISTIC DERIVATION OF CKM EXPONENTS
====================================================
We know n(b1,b2) = 9*b2 - 5*b1 - 6 gives {3,7,12} exactly.
But WHY does the CKM angle scale as phi^(-n)?

HYPOTHESIS: The CKM matrix comes from the overlap of generation
eigenstates on the 600-cell graph. The transition amplitude between
generations b1 and b2 is determined by the spectral overlap:

  <b1|V|b2> = sum_k f(lambda_k) * <b1|k><k|b2>

where |k> are the eigenstates of the 600-cell adjacency matrix
and V is the gauge interaction vertex.

The key question: do the generation eigenstates phi^(b*k) have an
overlap that decays as phi^(-n(b1,b2))?

APPROACHES:
A. Spectral overlap: define generation states from the eigenvalue structure
B. Green's function: G(z1, z2) on the 600-cell spectrum
C. Random walk: transition probability from gen-b1 to gen-b2
D. Perturbation theory: CKM as off-diagonal mass matrix element

STATUS: EXPLORATORY (attempting to upgrade from PATTERN to DERIVED)
Author: Razvan-Constantin Anghelina
Date: February 2026
"""

import numpy as np
import math
from itertools import permutations as iperms

PHI = (1 + math.sqrt(5)) / 2
ALPHA = 7.2973525693e-3
ALPHA_S = 0.1179

print("=" * 72)
print("EXP-209: MECHANISTIC DERIVATION OF CKM EXPONENTS")
print("=" * 72)
print()

# Build 600-cell
def build_600cell():
    verts = set()
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.add(tuple(v))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    verts.add((s0, s1, s2, s3))
    half = 0.5
    phi_half = PHI / 2.0
    iphi_half = 1.0 / (2.0 * PHI)
    base = [0.0, half, phi_half, iphi_half]
    even_perms = []
    for p in iperms(range(4)):
        inv = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if p[a] > p[b]:
                    inv += 1
        if inv % 2 == 0:
            even_perms.append(p)
    for p in even_perms:
        cu = [base[p[0]], base[p[1]], base[p[2]], base[p[3]]]
        nz = [i for i in range(4) if abs(cu[i]) > 1e-12]
        for sb in range(1 << len(nz)):
            v = list(cu)
            for k, pos in enumerate(nz):
                if (sb >> k) & 1:
                    v[pos] = -v[pos]
            verts.add(tuple(round(x, 12) for x in v))
    return sorted(verts)

vertices = build_600cell()
N = len(vertices)
varray = np.array(vertices)
dot_mat = varray @ varray.T
adj = (np.abs(dot_mat - PHI / 2.0) < 1e-8).astype(np.float64)
np.fill_diagonal(adj, 0)

eigvals_adj, eigvecs_adj = np.linalg.eigh(adj)
# Sort descending
idx_sort = np.argsort(eigvals_adj)[::-1]
eigvals_adj = eigvals_adj[idx_sort]
eigvecs_adj = eigvecs_adj[:, idx_sort]

print("  600-cell: %d vertices" % N)
print()

# ============================================================
# PART A: GENERATION STATES ON THE SPECTRUM
# ============================================================
print("=" * 72)
print("PART A: GENERATION STATES FROM EIGENVALUE STRUCTURE")
print("=" * 72)
print()

# The 9 eigenvalue clusters with their (a,b) parameters:
# k=0: lambda=12=12+0*phi, mult=1, a=12, b=0
# k=1: lambda=6*phi~9.71, mult=4, a=0, b=6
# k=2: lambda=4*phi~6.47, mult=9, a=0, b=4
# k=3: lambda=3, mult=16, a=3, b=0
# k=4: lambda=0, mult=25, a=0, b=0
# k=5: lambda=-2, mult=36, a=-2, b=0
# k=6: lambda=4-4*phi~-2.47, mult=9, a=4, b=-4
# k=7: lambda=-3, mult=16, a=-3, b=0
# k=8: lambda=6-6*phi~-3.71, mult=4, a=6, b=-6

# Generation quantum number: b in {0, 1, 2}
# The eigenvalues with b != 0 are: k=1(b=6), k=2(b=4), k=6(b=-4), k=8(b=-6)
# The eigenvalues with b = 0 are: k=0(a=12), k=3(a=3), k=4(a=0), k=5(a=-2), k=7(a=-3)

# Group eigenvalues by cluster
clusters = []
tol = 0.01
for i in range(N):
    ev = eigvals_adj[i]
    found = False
    for cl in clusters:
        if abs(ev - cl[0]) < tol:
            cl.append(i)  # store index
            found = True
            break
    if not found:
        clusters.append([i])

print("  Eigenvalue clusters:")
for ci, cl in enumerate(clusters):
    ev = eigvals_adj[cl[0]]
    # (a,b) decomposition
    for a in range(-10, 15):
        b_test = (ev - a) / PHI
        if abs(b_test - round(b_test)) < 0.01:
            b = int(round(b_test))
            break
    n = 5*a + 6*b
    print("  k=%d: lambda=%.4f, (a,b)=(%d,%d), n=%d, mult=%d" % (
        ci, ev, a, b, n, len(cl)))

# The key idea: each eigenvalue has n = 5a + 6b.
# A "generation-b_gen state" is the projection onto eigenspaces
# where the eigenvalue parameter b matches.
# But b is the spectral parameter, not the generation number!

# From the generation theorem:
# Generation 1: b_gen=0, z=1, eigenvalues with b=0 in (a,b)
# Generation 2: b_gen=1, z=1+phi, eigenvalues with b=1?
# No, this doesn't work directly because the spectral b values are {0,4,6,-4,-6}

# Alternative: the generation states are LINEAR COMBINATIONS of eigenstates
# weighted by phi^(b_gen * something)

# Let me try: define generation-b_gen projector as
# P_{b_gen} = sum_k w_k(b_gen) * |k><k|
# where w_k(b_gen) depends on b_gen and the eigenvalue parameters

# For the CKM matrix: V_{ij} ~ <P_i | gauge | P_j>
# If the "gauge" operator is the adjacency matrix itself, then
# V_{ij} = sum_k lambda_k * <i|k><k|j> (spectral decomposition)

# Try the simplest: V_{b1,b2} = (1/N) * sum_k lambda_k^p * exp(i*(b2-b1)*theta_k)
# where theta_k is some phase associated with eigenvalue k

# Actually, let me try a more physical approach.
# The adjacency matrix A has spectral decomposition A = sum_k lambda_k P_k
# where P_k is the projector onto eigenspace k.

# Define: Phi(b_gen, k) = response of generation b_gen to eigenvalue k
# Phi(b, k) = phi^(b * b_k) where b_k is the spectral b-parameter

print()
print("  Defining generation-eigenvalue coupling:")
print("  Phi(b_gen, k) = phi^(b_gen * b_k)")
print()

# Spectral b-values
b_spectral = [0, 6, 4, 0, 0, 0, -4, 0, -6]
mults = [1, 4, 9, 16, 25, 36, 9, 16, 4]
lambdas_unique = [eigvals_adj[cl[0]] for cl in clusters]

print("  Phi values:")
print("  %8s" % "b_gen" + "".join("  k=%d" % k for k in range(9)))
for b_gen in [0, 1, 2]:
    vals = [PHI**(b_gen * b_spectral[k]) for k in range(9)]
    line = "  %8d" % b_gen + "".join("%8.3f" % v for v in vals)
    print(line)

# Now: CKM element V_{b1, b2} = sum_k Phi(b1, k) * Phi(b2, k) * mult_k * f(lambda_k)
# where f(lambda_k) is some weighting function

# The simplest: f = 1/N (uniform)
# V_{b1,b2} = (1/N) * sum_k phi^((b1+b2)*b_k) * mult_k
print()
print("  Overlap O(b1,b2) = sum_k phi^((b1+b2)*b_k) * m_k / N:")

for b1 in range(3):
    for b2 in range(b1, 3):
        overlap = sum(PHI**((b1+b2)*b_spectral[k]) * mults[k] for k in range(9)) / N
        print("  O(%d,%d) = %.6f" % (b1, b2, overlap))

# That doesn't give small numbers for off-diagonal. Let me try difference:
print()
print("  Overlap O(b1,b2) = sum_k phi^((b2-b1)*b_k) * m_k / N:")
for b1 in range(3):
    for b2 in range(b1, 3):
        db = b2 - b1
        overlap = sum(PHI**(db*b_spectral[k]) * mults[k] for k in range(9)) / N
        print("  O(%d,%d) = %.6f" % (b1, b2, overlap))

# ============================================================
# PART B: GREEN'S FUNCTION APPROACH
# ============================================================
print()
print("=" * 72)
print("PART B: GREEN'S FUNCTION ON SPECTRUM")
print("=" * 72)
print()

# Define Green's function: G(z1, z2; t) = sum_k P_k(z1,z2) * exp(-lambda_k * t)
# where z_gen = 1 + b_gen*phi is the generation z-value
# and P_k(z1, z2) is the spectral projection weighted by z

# This is abstract. Let me try: the HEAT KERNEL on the 600-cell graph
# H(v, w; t) = <v| exp(-t*L) |w> = sum_k exp(-t*(12-lambda_k)) * psi_k(v)*psi_k(w)

# For a vertex-transitive graph, if we average over starting points:
# H_avg(d; t) = sum over w at distance d of H(v,w;t) / (count at d)
# This gives a function of distance d and "time" t

# The CKM angle could be related to the off-diagonal elements of H
# at some characteristic time t* (e.g., t* = 1/12 or t* = 1/(gap))

# Let's compute the heat kernel averaged over association scheme classes
L = np.diag(adj.sum(axis=1)) - adj  # Laplacian
eigvals_L = 12 - eigvals_adj  # Laplacian eigenvalues

# Association scheme classes (8 classes by inner product)
ip_classes = {}
for j in range(1, N):
    ip = round(dot_mat[0, j], 4)
    found = False
    for key in ip_classes:
        if abs(ip - key) < 0.001:
            ip_classes[key].append(j)
            found = True
            break
    if not found:
        ip_classes[ip] = [j]

print("  Association scheme classes (from vertex 0):")
for ip in sorted(ip_classes.keys(), reverse=True):
    verts_in = ip_classes[ip]
    print("  ip=%.4f: %d vertices" % (ip, len(verts_in)))

# For each class, compute the heat kernel at various times
print()
print("  Heat kernel H(class; t) for t = 1/12, 1/7, 1/5:")
times = [1.0/12, 1.0/7, 1.0/5, 1.0/3, 1.0]

exp_tL = {}
for t in times:
    # exp(-t*L) = V * diag(exp(-t*lambda)) * V^T
    exp_evals = np.exp(-t * eigvals_L)
    H_mat = eigvecs_adj @ np.diag(exp_evals) @ eigvecs_adj.T
    exp_tL[t] = H_mat

print("  %10s" % "ip" + "".join("  t=%.3f" % t for t in times))
for ip in sorted(ip_classes.keys(), reverse=True):
    verts_in = ip_classes[ip]
    vals = []
    for t in times:
        H = exp_tL[t]
        h_avg = np.mean([H[0, j] for j in verts_in])
        vals.append(h_avg)
    line = "  %10.4f" % ip + "".join("%10.4e" % v for v in vals)
    print(line)

# ============================================================
# PART C: TRANSITION AMPLITUDES BETWEEN GENERATIONS
# ============================================================
print()
print("=" * 72)
print("PART C: TRANSITION AMPLITUDES BETWEEN GENERATIONS")
print("=" * 72)
print()

# Here's the key idea:
# Each generation b_gen = {0, 1, 2} corresponds to a z-value z = 1 + b*phi
# In the spectrum, eigenvalue k has (a_k, b_k) with n_k = 5*a_k + 6*b_k
# The generation "response function" to eigenvalue k is:
# R(b_gen, k) = phi^(n_k * b_gen / (a_1+b_1-2))
# = phi^(n_k * b_gen / 9)
# (dividing by N_eig = 9 to normalize)

# Or more simply: R(b_gen, k) = exp(2*pi*i*b_gen*k/9) (mod 9 Fourier)
# Since there are 9 eigenvalue clusters and 3 generations, and 9 = 3*3

# Try: V_{b1,b2} = sum_k R(b1,k)* R(b2,k) * w_k
# where w_k is some weight from the eigenvalue

# For the Fourier approach:
# V_{b1,b2} = (1/9) * sum_{k=0}^{8} exp(2*pi*i*(b2-b1)*k/9) * w_k

print("  Fourier approach: V(db) = (1/9) * sum_k exp(2*pi*i*db*k/9) * w_k")
print()

for db_name, db in [("db=0 (diagonal)", 0), ("db=1 (1-2 or 2-3)", 1), ("db=2 (1-3)", 2)]:
    for weight_name, weights in [
        ("uniform", [1]*9),
        ("lambda_k", lambdas_unique),
        ("mult_k", mults),
        ("lambda*mult", [l*m for l,m in zip(lambdas_unique, mults)]),
        ("|lambda_k|", [abs(l) for l in lambdas_unique]),
    ]:
        V = sum(np.exp(2j*np.pi*db*k/9) * weights[k] for k in range(9)) / 9
        if abs(V) > 1e-10:
            mag = abs(V)
            # Compare with phi^-n
            if mag > 0 and mag < 1:
                n_eff = -math.log(mag)/math.log(PHI)
            elif mag >= 1:
                n_eff = math.log(mag)/math.log(PHI)
            else:
                n_eff = 999
            pass  # store for later
        else:
            mag = abs(V)
            n_eff = 999

        if abs(mag) > 1e-12:
            print("  %s w=%s: |V| = %.6f = phi^(%.2f)" % (
                db_name[:4], weight_name[:10], mag,
                math.log(mag)/math.log(PHI) if mag > 0 else 0))

print()

# ============================================================
# PART D: SPECTRAL ZETA FUNCTION APPROACH
# ============================================================
print()
print("=" * 72)
print("PART D: SPECTRAL ZETA FUNCTION")
print("=" * 72)
print()

# The spectral zeta function: zeta_G(s) = sum_k mult_k * lambda_k^(-s)
# (sum over non-zero eigenvalues)

# Idea: the CKM angle is related to the spectral zeta function
# evaluated at s = n(b1,b2) somehow

# For each non-zero eigenvalue:
print("  Non-zero eigenvalues and powers:")
nz_lambdas = [(l, m) for l, m in zip(lambdas_unique, mults) if abs(l) > 0.01]

for s in [3, 7, 12]:
    zeta_s = sum(m * abs(l)**(-s) for l, m in nz_lambdas)
    zeta_s_pos = sum(m * l**(-s) for l, m in nz_lambdas if l > 0.01)
    print("  zeta(%d) = %.6e (positive only: %.6e)" % (s, zeta_s, zeta_s_pos))

print()

# Also: the Ihara zeta function of the graph
# And: transition probability after n steps
# P(v -> w in n steps) = (A^n)_{vw} / degree^n

# For n = 3, 7, 12: what is the average transition probability
# to vertices at various distances?
print("  Transition probabilities (A^n / 12^n):")
A_norm = adj / 12  # normalized adjacency (transition matrix)

for n_step in [3, 7, 12]:
    An = np.linalg.matrix_power(A_norm, n_step)
    print("  n = %d steps:" % n_step)
    for ip in sorted(ip_classes.keys(), reverse=True):
        verts_in = ip_classes[ip]
        p_avg = np.mean([An[0, j] for j in verts_in])
        if abs(p_avg) > 1e-15:
            print("    ip=%.4f: P = %.6e" % (ip, p_avg))

print()

# ============================================================
# PART E: THE n = 5a+6b CONNECTION
# ============================================================
print()
print("=" * 72)
print("PART E: DIRECT SPECTRAL n-VALUE ANALYSIS")
print("=" * 72)
print()

# Each eigenvalue has n_k = 5*a_k + 6*b_k
# The CKM exponent is n(b1,b2) = 9*b2 - 5*b1 - 6
# Can we express this as a function of the n_k values?

n_values = []
for ci, cl in enumerate(clusters):
    ev = eigvals_adj[cl[0]]
    for a in range(-10, 15):
        b_test = (ev - a) / PHI
        if abs(b_test - round(b_test)) < 0.01:
            b = int(round(b_test))
            break
    n_k = 5*a + 6*b
    n_values.append((ci, ev, a, b, n_k, len(cl)))

print("  Eigenvalue n-values:")
for ci, ev, a, b, n_k, mult in n_values:
    print("  k=%d: lambda=%.4f, (a,b)=(%d,%d), n=%d, mult=%d" % (
        ci, ev, a, b, n_k, mult))

# The n-values: [60, 36, 24, 15, 0, -10, -4, -15, -6]
nvals = [nv[4] for nv in n_values]
print("\n  n-values:", nvals)

# Can we construct {3, 7, 12} from these?
# 36 - 24 = 12! (adjacent pair)
# 15 - 0 = 15? No.
# 24 - 15 = 9 = N_eig. Close to but not 3 or 7.
# 0 - (-4) = 4. Not 3 or 7.
# 60 - 36 = 24. Not 3, 7, or 12.
# -6 - (-10) = 4.
# 36 - 24 = 12 YES!
# -4 - (-10) = 6 = b_1
# 0 - (-4) = 4
# 15 - 0 = 15

# Differences between ALL pairs:
print("\n  All pairwise differences of n-values:")
diffs_set = set()
for i, ni in enumerate(nvals):
    for j, nj in enumerate(nvals):
        d = ni - nj
        if d > 0:
            diffs_set.add(d)
print("  Positive differences:", sorted(diffs_set))
print("  Contains 3:", 3 in diffs_set)
print("  Contains 7:", 7 in diffs_set)
print("  Contains 12:", 12 in diffs_set)

# Find which pairs give 3, 7, 12
for target in [3, 7, 12]:
    print("  n-diff = %d:" % target)
    for i, ni in enumerate(nvals):
        for j, nj in enumerate(nvals):
            if ni - nj == target:
                print("    n_%d - n_%d = %d - %d = %d" % (i, j, ni, nj, target))

print()

# Interesting: 12 = n_1 - n_2 = 36 - 24
# But 3 and 7 don't appear as n-differences!
# Instead: 3 = sum of smaller pieces?
# 3 = -(-4) + (-10) + 0 - (-6) + ... no clean way

# Let's try a different approach:
# The CKM exponent n(b1,b2) = 9*b2 - 5*b1 - 6
# = N_eig * b2 - diam * b1 - b_1
# This can be rewritten as:
# n = 9*(b2 - b1) + 4*b1 - 6
# For b2 - b1 = 1: n = 9 + 4*b1 - 6 = 3 + 4*b1
#   b1=0: n=3, b1=1: n=7
# For b2 - b1 = 2: n = 18 + 4*b1 - 6 = 12 + 4*b1
#   b1=0: n=12

print("  FACTORED FORM: n = 9*db + 4*b1 - 6 where db = b2 - b1")
print("  For db=1: n = 3 + 4*b1 (b1=0 => 3, b1=1 => 7)")
print("  For db=2: n = 12 + 4*b1 (b1=0 => 12)")
print()

# The number 4 = 2^2 = ... or 4 = a_1 - 1
# 9 = N_eig
# 6 = b_1
# 3 = 9 - 6 = N_eig - b_1

# So: n = N_eig * db + (a_1-1)*b1 - b_1
# = N_eig*db + 4*b1 - 6
# n(0,1) = 9 - 6 = 3 = N_eig - b_1
# n(1,2) = 9 + 4 - 6 = 7 = N_eig + (a_1-1) - b_1 = 9 + 4 - 6
# n(0,2) = 18 - 6 = 12 = 2*N_eig - b_1 = degree

print("  n = N_eig * db + (a_1-1) * b1 - b_1")
print("  = 9*db + 4*b1 - 6")
print()
print("  Remarkable: n(0,2) = 2*N_eig - b_1 = 18-6 = 12 = degree!")
print("  The 1-3 mixing exponent equals the degree of the graph!")
print()

# ============================================================
# PART F: RANDOM WALK AND PHI^-n DECAY
# ============================================================
print()
print("=" * 72)
print("PART F: RANDOM WALK PHI^-n DECAY")
print("=" * 72)
print()

# On a regular graph, the random walk transition probability
# P(v->w in t steps) ~ sum_k exp(-t*gap_k) * f_k(d(v,w))
# For large t, this is dominated by the spectral gap.

# The spectral gap of the 600-cell adjacency: lambda_1/lambda_0 = 9.708/12 = 0.809
# This equals phi/2! (phi/2 = 0.809)

spectral_ratio = eigvals_adj[1] / eigvals_adj[0]
print("  lambda_1/lambda_0 = %.6f" % spectral_ratio)
print("  phi/2 = %.6f" % (PHI/2))
print("  Match: %.4f%%" % (100*abs(spectral_ratio - PHI/2)/(PHI/2)))
print()

# Wait: lambda_1 = 6*phi and lambda_0 = 12
# lambda_1/lambda_0 = 6*phi/12 = phi/2 EXACTLY!
print("  EXACT: lambda_1/lambda_0 = 6*phi/12 = phi/2")
print()

# After t steps of random walk:
# P(t) ~ (lambda_1/lambda_0)^t = (phi/2)^t
# For phi/2 = phi^(1-log2/logphi):
# (phi/2)^t = phi^(t*(1-log2/logphi))
# 1 - log2/logphi = 1 - 0.6931/0.4812 = 1 - 1.440 = -0.440
# So (phi/2)^t = phi^(-0.440*t)

# For the CKM: if the "time" t is related to the generation quantum numbers,
# then we need t such that 0.440*t = n(b1,b2)
# n=3: t = 3/0.440 = 6.82
# n=7: t = 7/0.440 = 15.91
# n=12: t = 12/0.440 = 27.27

# Not clean integers. This approach doesn't work directly.

# But: if the decay involves phi directly (not phi/2), then
# P ~ phi^(-n) requires a different mechanism.

# Key observation: on the 600-cell, the random walk has a GOLDEN RATIO
# decay because the spectral ratio is phi/2.
# The factor of 2 is the "entropic" contribution from the degree-12 graph.

# What if the CKM involves the RESOLVENT rather than the heat kernel?
# G(z) = (z - A)^{-1} has poles at the eigenvalues.
# If z = phi (the golden ratio itself), then:
# G(phi) = sum_k 1/(phi - lambda_k) * P_k

# At the eigenvalue 6*phi (lambda_1): phi - 6*phi = -5*phi
# At lambda_2 = 4*phi: phi - 4*phi = -3*phi
# At lambda_0 = 12: phi - 12 = phi - 12 = -10.382

# The resolvent at z = phi^m for various m:
print("  Resolvent trace at z = phi^m:")
for m in range(0, 10):
    z = PHI**m
    trace_G = sum(mult / (z - lam) for lam, mult in zip(lambdas_unique, mults) if abs(z-lam) > 0.01)
    if abs(trace_G) > 0.01:
        print("  z = phi^%d = %.4f: Tr(G) = %.6f = phi^(%.3f)" % (
            m, z, trace_G, math.log(abs(trace_G))/math.log(PHI)))

print()

# ============================================================
# CONCLUSIONS
# ============================================================
print()
print("=" * 72)
print("CONCLUSIONS")
print("=" * 72)
print()

print("1. SPECTRAL RATIO: lambda_1/lambda_0 = phi/2 EXACTLY")
print("   This gives golden-ratio decay in random walks.")
print("   STATUS: DERIVED")
print()

print("2. FOURIER APPROACH: DFT with 9 eigenvalues and 3 generations")
print("   gives overlap integrals, but they don't match phi^-n directly.")
print("   STATUS: NEGATIVE (this specific approach)")
print()

print("3. FACTORED FORM: n = N_eig*db + (a_1-1)*b1 - b_1")
print("   = 9*(b2-b1) + 4*b1 - 6")
print("   The base step is 9 per generation gap, shifted by 4*b1.")
print("   n(0,2) = 12 = degree (maximal mixing exponent = degree)")
print("   STATUS: PATTERN (reformulation, not derivation)")
print()

print("4. The CKM exponents {3,7,12} cannot be extracted from")
print("   simple spectral operations (DFT, heat kernel, resolvent)")
print("   on the 600-cell. The formula n = 9*b2 - 5*b1 - 6 remains")
print("   at PATTERN level. A true derivation likely requires")
print("   understanding the DYNAMICS of generation mixing on the")
print("   soliton background, which is beyond pure graph theory.")
print()

print("5. KEY INSIGHT: n(0,2) = 12 = degree of the graph")
print("   This means the 1-3 mixing is maximally suppressed")
print("   (exponent equals the coordination number).")
print("   And n = 9*db + 4*b1 - 6 = N_eig*db + (a_1-1)*b1 - b_1")
print("   has ALL coefficients from the 600-cell geometry.")
print()

print("EXP-209 COMPLETE")
