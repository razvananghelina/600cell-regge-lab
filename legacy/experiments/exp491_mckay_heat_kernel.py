"""
EXP-491: Heat Kernel on McKay Graph -> Mass Hierarchy?
=======================================================
CONTEXT: D_F eigenvalues -> masses FAILED (exp307-362).
Wilson lines WORK but are not operator-spectral.
Heat kernel / Green's function: UNEXPLORED direction.

IDEA: The heat kernel K(t, i, j) = sum_k exp(-t * lambda_k^2) * v_k(i) * v_k(j)
on the McKay graph (9 nodes, extended E8) propagates from the
vacuum (r0) to each node. At some critical time t*, the propagator
ratios K(t*, 0, j) / K(t*, 0, 0) might reproduce the mass hierarchy.

WHY THIS COULD WORK:
  - McKay eigenvalues: {2, phi, 1, 1/phi, 0, -1/phi, -1, -phi, -2}
  - exp(-t*phi^2) / exp(-t*1) = exp(-t*(phi^2 - 1)) = exp(-t*phi)
    (since phi^2 = phi + 1, so phi^2 - 1 = phi)
  - This naturally generates phi-dependent ratios!
  - Wilson line weights come from PATH SUMS, heat kernel from DIFFUSION
  - Both are "transport" phenomena on the McKay graph

ALSO TEST: Green's function G = (D_F^2)^{-1} or resolvent (D_F^2 - z)^{-1}

STATUS: EXPLORATORY (first attempt at unexplored direction)
"""

import numpy as np
from scipy.linalg import expm

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6

print("=" * 70)
print("EXP-491: HEAT KERNEL ON McKAY GRAPH")
print("=" * 70)

# ============================================================
# Build McKay adjacency matrix
# ============================================================
# Ordering: r0(1), r1(2), r2(3), r3(4), r4(5), r5(6), r6(3), r7(4), r8(2)
dims = np.array([1, 2, 3, 4, 5, 6, 3, 4, 2])
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,7), (7,8)]
n_nodes = 9

M = np.zeros((n_nodes, n_nodes))
for a, b in edges:
    M[a, b] = 1; M[b, a] = 1

# McKay Laplacian
L = np.diag(M.sum(axis=1)) - M

# Eigendecomposition
evals_M, evecs_M = np.linalg.eigh(M)
idx = np.argsort(evals_M)[::-1]
evals_M = evals_M[idx]
evecs_M = evecs_M[:, idx]

evals_L, evecs_L = np.linalg.eigh(L)
idx_L = np.argsort(evals_L)
evals_L = evals_L[idx_L]
evecs_L = evecs_L[:, idx_L]

print(f"\nMcKay adjacency eigenvalues: {evals_M.round(6)}")
print(f"McKay Laplacian eigenvalues: {evals_L.round(6)}")

# Known fermion data
fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_n = np.array([0, 11, 17, 3, 16, 26, 5, 11, 19])
fermion_mass = 0.511 * phi**fermion_n  # MeV

# Wilson line edge weights on main chain (from exp323)
wl_weights_main = [3, 2, 6, 0, 5]  # edges r0-r1, r1-r2, ..., r4-r5
# Branch: r5-r6 weight and r5-r7, r7-r8 weights unknown from this notation
# Full path exponents:
# e -> r0: n=0
# d -> r0-r1-r2: n=5 (path sum = 3+2)
# u -> r0-r1: n=3 (path sum = 3)
# s,mu -> r0-r1-r2-r3: n=11 (path sum = 3+2+6)
# c -> r0-r1-r2-r3-r4-r5: n=16 (path sum = 3+2+6+0+5)
# tau -> r0-...-r5-r6: n=17 (path = 16+1? or different)
# b -> path to branch: n=19
# t -> path to end: n=26

# Actually, the mapping node -> fermion isn't established.
# Let's just compute the heat kernel and see what ratios it gives.

# ============================================================
# SECTION 1: Heat kernel on McKay adjacency
# ============================================================
print("\n" + "=" * 70)
print("SECTION 1: HEAT KERNEL K(t) = exp(-t * M^2)")
print("=" * 70)

# K(t) = exp(-t * M^2) since M is the "Dirac" operator, M^2 is the "Laplacian"
# K(t, i, j) = sum_k exp(-t * lambda_k^2) * v_k(i) * v_k(j)

# Eigenvalues of M^2: {4, phi^2, 1, 1/phi^2, 0, 1/phi^2, 1, phi^2, 4}
# = {0, 1/phi^2, 1/phi^2, 1, 1, phi^2, phi^2, 4, 4}
# = {0, phi^2-2*phi+1, ..., phi+1, ..., 4}
# Wait: 1/phi^2 = (phi-1)^2 = phi^2 - 2*phi + 1 ... no.
# 1/phi = phi - 1, so 1/phi^2 = (phi-1)^2 = phi^2 - 2*phi + 1 = (phi+1) - 2*phi + 1 = 2 - phi
# And phi^2 = phi + 1

M2_evals = evals_M**2
print(f"\nM^2 eigenvalues: {np.sort(M2_evals).round(6)}")
print(f"  = {{0, {2-phi:.4f}(=2-phi), {2-phi:.4f}, 1, 1, {phi+1:.4f}(=phi+1), {phi+1:.4f}, 4, 4}}")

# Compute K(t, 0, j) for various t
print(f"\nHeat kernel K(t, r0, rj) for various t:")
print(f"{'t':>8}", end="")
for j in range(n_nodes):
    print(f"  r{j}({dims[j]})", end="")
print(f"  {'ratio max/min':>14}")

M2 = M @ M
interesting_t = []

for t in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0,
          1/phi, phi, 1/(phi**2), phi**2, np.log(phi),
          1/(2*phi), 2*phi]:
    Kt = expm(-t * M2)
    K0j = Kt[0, :]

    nonzero = K0j[K0j > 1e-15]
    if len(nonzero) > 1:
        ratio = np.max(nonzero) / np.min(nonzero)
    else:
        ratio = 0

    print(f"{t:8.4f}", end="")
    for j in range(n_nodes):
        print(f"  {K0j[j]:8.5f}", end="")
    print(f"  {ratio:14.2f}")

    interesting_t.append((t, K0j.copy(), ratio))

# ============================================================
# SECTION 2: Log-ratios vs mass exponents
# ============================================================
print("\n" + "=" * 70)
print("SECTION 2: LOG-RATIOS OF HEAT KERNEL")
print("=" * 70)

print("""
If K(t, 0, j) ~ phi^{-n_j}, then -log_phi(K(t,0,j)/K(t,0,0)) = n_j.
We need to find t such that these log-ratios match mass exponents.
""")

# For each t, compute -log_phi(K(t,0,j)/K(t,0,0)) and compare with sorted n values
print(f"{'t':>8}  {'log-ratios -log_phi(K(t,0,j)/K(t,0,0))':>50}")

best_t = None
best_score = float('inf')

for t_val in np.concatenate([np.linspace(0.01, 10, 200),
                              [1/phi, phi, phi**2, 1/(phi**2), np.log(phi)]]):
    Kt = expm(-t_val * M2)
    K0j = Kt[0, :]

    # Compute log ratios
    log_ratios = []
    for j in range(n_nodes):
        if K0j[j] > 1e-20 and K0j[0] > 1e-20:
            ratio = K0j[j] / K0j[0]
            if ratio > 0:
                lr = -np.log(ratio) / np.log(phi)
            else:
                lr = float('inf')
        else:
            lr = float('inf')
        log_ratios.append(lr)

    log_ratios = np.array(log_ratios)

    # Sort and compare with sorted fermion exponents
    finite_lr = sorted([lr for lr in log_ratios if np.isfinite(lr)])
    sorted_n = sorted(fermion_n)

    # Try to find best match by minimizing sum of squared differences
    if len(finite_lr) >= 9:
        # Try all permutations... too many. Just compare sorted.
        score = sum((finite_lr[i] - sorted_n[i])**2 for i in range(min(len(finite_lr), 9)))
        if score < best_score:
            best_score = score
            best_t = t_val

if best_t is not None:
    print(f"\nBest t = {best_t:.6f} (score = {best_score:.4f})")
    Kt = expm(-best_t * M2)
    K0j = Kt[0, :]

    log_ratios = []
    for j in range(n_nodes):
        if K0j[j] > 1e-20 and K0j[0] > 1e-20:
            ratio = K0j[j] / K0j[0]
            lr = -np.log(max(ratio, 1e-100)) / np.log(phi)
        else:
            lr = float('inf')
        log_ratios.append(lr)

    log_ratios_sorted = sorted([(lr, j) for j, lr in enumerate(log_ratios) if np.isfinite(lr)])

    print(f"\n  {'Node':>6} {'K(t,0,j)':>12} {'-log_phi ratio':>15} {'Nearest n':>10} {'delta':>8}")
    for lr, j in log_ratios_sorted:
        nearest_n = min(fermion_n, key=lambda x: abs(x - lr))
        delta = lr - nearest_n
        print(f"  r{j}({dims[j]}) {K0j[j]:12.6e} {lr:15.4f} {nearest_n:10d} {delta:8.4f}")

# ============================================================
# SECTION 3: Heat kernel on McKay LAPLACIAN
# ============================================================
print("\n" + "=" * 70)
print("SECTION 3: HEAT KERNEL K(t) = exp(-t * L)")
print("=" * 70)

print("Using the graph Laplacian L = D - M instead of M^2.")
print(f"Laplacian eigenvalues: {evals_L.round(4)}")

best_t_L = None
best_score_L = float('inf')

for t_val in np.concatenate([np.linspace(0.01, 20, 400),
                              np.linspace(0.001, 0.01, 20),
                              [1/phi, phi, phi**2, np.log(phi), 5*np.log(phi)]]):
    Kt = expm(-t_val * L)
    K0j = Kt[0, :]

    log_ratios = []
    for j in range(n_nodes):
        if K0j[j] > 1e-30 and K0j[0] > 1e-30:
            ratio = K0j[j] / K0j[0]
            if ratio > 0:
                lr = -np.log(ratio) / np.log(phi)
            else:
                lr = float('inf')
        else:
            lr = float('inf')
        log_ratios.append(lr)

    finite_lr = sorted([lr for lr in log_ratios if np.isfinite(lr)])
    sorted_n = sorted(fermion_n)

    if len(finite_lr) >= 9:
        score = sum((finite_lr[i] - sorted_n[i])**2 for i in range(min(len(finite_lr), 9)))
        if score < best_score_L:
            best_score_L = score
            best_t_L = t_val

print(f"\nBest t = {best_t_L:.6f} (score = {best_score_L:.4f})")
Kt = expm(-best_t_L * L)
K0j = Kt[0, :]

log_ratios = []
for j in range(n_nodes):
    if K0j[j] > 1e-30 and K0j[0] > 1e-30:
        ratio = K0j[j] / K0j[0]
        lr = -np.log(max(ratio, 1e-100)) / np.log(phi)
    else:
        lr = float('inf')
    log_ratios.append(lr)

log_ratios_sorted = sorted([(lr, j) for j, lr in enumerate(log_ratios) if np.isfinite(lr)])

print(f"\n  {'Node':>6} {'K(t,0,j)':>12} {'-log_phi ratio':>15} {'Nearest n':>10} {'delta':>8}")
for lr, j in log_ratios_sorted:
    nearest_n = min(fermion_n, key=lambda x: abs(x - lr))
    delta = lr - nearest_n
    print(f"  r{j}({dims[j]}) {K0j[j]:12.6e} {lr:15.4f} {nearest_n:10d} {delta:8.4f}")

# ============================================================
# SECTION 4: WEIGHTED McKay graph
# ============================================================
print("\n" + "=" * 70)
print("SECTION 4: WEIGHTED McKAY GRAPH (Wilson line weights)")
print("=" * 70)

print("""
The Wilson line weights on the main chain are [3, 2, 6, 0, 5].
What if we use these as EDGE WEIGHTS in the Laplacian?
L_w[i,j] = -w_{ij} for adjacent i,j, L_w[i,i] = sum_j w_{ij}
""")

# We need the full 8 edge weights. From exp323:
# Main chain r0-r1-r2-r3-r4-r5: weights [3, 2, 6, 0, 5]
# Branches from r5: r5-r6 and r5-r7-r8
# The mass exponents by node assignment determine these.
# But we don't KNOW the node->fermion mapping!

# Let's TRY: if the graph distance from r0 gives the sorted masses,
# then weights are the DIFFERENCES of consecutive exponents.
# Sorted exponents: 0, 3, 5, 11, 11, 16, 17, 19, 26
# Differences: 3, 2, 6, 0, 5, 1, 2, 7

# But the graph has only 8 edges and 9 nodes (tree).
# Sorted by graph distance:
# d=0: r0 (1 node)
# d=1: r1 (1 node)
# d=2: r2 (1 node)
# d=3: r3 (1 node)
# d=4: r4 (1 node)
# d=5: r5 (1 node)
# d=6: r6, r7 (2 nodes - the BRANCH!)
# d=7: r8 (1 node)

# If we assign fermions in order of increasing mass exponent
# to nodes in order of increasing graph distance:
# n=0 -> r0 (d=0): e
# n=3 -> r1 (d=1): u
# n=5 -> r2 (d=2): d
# n=11 -> r3 (d=3): mu or s (both n=11!)
# n=11 -> ??? (need another d=3 node, but there isn't one)

# Problem: two fermions (mu, s) have the SAME exponent n=11.
# The McKay graph has no two nodes at the same distance from r0
# except r6 and r7 (both at d=6).

# Alternative: assign by something other than distance.
# But the branching at r5 gives two branches:
# Branch A: r5-r6 (length 1)
# Branch B: r5-r7-r8 (length 2)

# Possible assignment using the branch:
# Main chain r0->r5 for 6 fermions
# Branch A (r6) for 1 fermion
# Branch B (r7, r8) for 2 fermions
# Total: 6 + 1 + 2 = 9. Works!

# Let's try: main chain = {e, u, d, mu/s, c(?), ...}
# The two n=11 fermions (mu, s) are in DIFFERENT sectors (lepton, down).
# Perhaps one goes on main chain, the other on a branch?

# Actually, let me just test with KNOWN Wilson line weights on main chain
# and see if heat kernel on the weighted graph matches better.

# Test with main chain weights [3, 2, 6, 0, 5] and guess branch weights
for w_56, w_57, w_78 in [(1, 3, 7), (1, 2, 7), (3, 3, 7), (1, 3, 10), (2, 1, 8)]:
    M_w = np.zeros((9, 9))
    main_weights = [3, 2, 6, 0, 5]  # edges 0-1, 1-2, 2-3, 3-4, 4-5
    for idx, (a, b) in enumerate(edges[:5]):
        M_w[a, b] = main_weights[idx]
        M_w[b, a] = main_weights[idx]
    # Branches
    M_w[5, 6] = w_56; M_w[6, 5] = w_56  # r5-r6
    M_w[5, 7] = w_57; M_w[7, 5] = w_57  # r5-r7
    M_w[7, 8] = w_78; M_w[8, 7] = w_78  # r7-r8

    L_w = np.diag(M_w.sum(axis=1)) - M_w

    # Find best t
    best_t_w = None
    best_score_w = float('inf')

    for t_val in np.linspace(0.01, 5, 500):
        Kt = expm(-t_val * L_w)
        K0j = Kt[0, :]

        log_ratios = []
        for j in range(n_nodes):
            if K0j[j] > 1e-30 and K0j[0] > 1e-30:
                ratio = K0j[j] / K0j[0]
                lr = -np.log(max(ratio, 1e-100)) / np.log(phi)
            else:
                lr = float('inf')
            log_ratios.append(lr)

        finite_lr = sorted([lr for lr in log_ratios if np.isfinite(lr)])
        sorted_n = sorted(fermion_n)

        if len(finite_lr) >= 9:
            score = sum((finite_lr[i] - sorted_n[i])**2 for i in range(9))
            if score < best_score_w:
                best_score_w = score
                best_t_w = t_val

    # Show result
    Kt = expm(-best_t_w * L_w)
    K0j = Kt[0, :]
    log_ratios = [-np.log(max(K0j[j]/K0j[0], 1e-100))/np.log(phi)
                  if K0j[j] > 1e-30 else float('inf') for j in range(9)]

    print(f"\n  Weights [{','.join(str(w) for w in main_weights)}|{w_56},{w_57},{w_78}], "
          f"best t={best_t_w:.3f}, score={best_score_w:.2f}")
    for j in range(9):
        if np.isfinite(log_ratios[j]):
            print(f"    r{j}({dims[j]}): -log_phi = {log_ratios[j]:.2f}", end="")
            if abs(log_ratios[j]) < 0.01:
                print(f"  [= 0 = n_e]")
            else:
                nearest = min(fermion_n, key=lambda x: abs(x - log_ratios[j]))
                print(f"  [nearest n = {nearest}]")

# ============================================================
# SECTION 5: Green's function approach
# ============================================================
print("\n" + "=" * 70)
print("SECTION 5: GREEN'S FUNCTION G = L^{-1}")
print("=" * 70)

print("The Green's function G = L^+ (pseudoinverse of Laplacian)")
print("gives the 'effective resistance' between nodes.")

# Pseudoinverse of L
# L has one zero eigenvalue (constant vector), so use pseudoinverse
L_pinv = np.linalg.pinv(L)

print(f"\nGreen's function G(r0, rj) = L^+(0, j):")
for j in range(n_nodes):
    print(f"  G(r0, r{j}) = {L_pinv[0, j]:.6f}")

# Effective resistance between r0 and rj (on a tree = graph distance)
print(f"\nEffective resistance R(r0, rj) = G(0,0) + G(j,j) - 2*G(0,j):")
for j in range(n_nodes):
    R = L_pinv[0,0] + L_pinv[j,j] - 2*L_pinv[0,j]
    print(f"  R(r0, r{j}) = {R:.6f}  (graph distance = {sum(1 for a,b in edges if min(a,b) < j and max(a,b) <= j)})")
    # Note: for a tree, effective resistance = graph distance (unweighted)

# Log of Green's function
print(f"\n-log_phi(G(r0,rj) / G(r0,r0)):")
for j in range(n_nodes):
    if L_pinv[0, j] > 0 and L_pinv[0, 0] > 0:
        lr = -np.log(L_pinv[0, j] / L_pinv[0, 0]) / np.log(phi)
        print(f"  r{j}: {lr:.4f}")
    elif L_pinv[0, j] < 0:
        print(f"  r{j}: (negative, = {L_pinv[0, j]:.6f})")

# ============================================================
# SECTION 6: t = infinity limit (stationary distribution)
# ============================================================
print("\n" + "=" * 70)
print("SECTION 6: LARGE-t BEHAVIOR AND CRITICAL TIME")
print("=" * 70)

# At large t, exp(-t*L) -> projection onto zero eigenspace = uniform (1/9) for connected graph
# At t = 0, exp(-t*L) = I (delta function)
# The interesting regime is intermediate t where the propagator
# has spread but not yet equilibrated.

# The key: at what t does the propagator ratio K(t,0,j)/K(t,0,0)
# transition from order 1 to exponentially small?
# This is determined by the spectral gap of L.

print(f"Spectral gap of McKay Laplacian: {evals_L[1]:.6f}")
print(f"Inverse gap (relaxation time): {1/evals_L[1]:.4f}")
print(f"phi*relaxation time: {phi/evals_L[1]:.4f}")

# Detailed scan around relaxation time
print(f"\nDetailed scan: log-ratios at t near relaxation time")
t_relax = 1/evals_L[1]

for t_val in [t_relax * f for f in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]]:
    Kt = expm(-t_val * L)
    K0j = Kt[0, :]
    log_ratios = [-np.log(max(K0j[j]/K0j[0], 1e-100))/np.log(phi)
                  for j in range(9)]
    max_lr = max(lr for lr in log_ratios if np.isfinite(lr))
    print(f"  t = {t_val:8.3f} ({t_val/t_relax:.1f} * t_relax): "
          f"max log-ratio = {max_lr:.2f}, "
          f"range needed: 26 (for top quark)")

# ============================================================
# SECTION 7: SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SECTION 7: SUMMARY")
print("=" * 70)

print(f"""
RESULTS:
========

1. UNWEIGHTED HEAT KERNEL (M^2 or L):
   The heat kernel on the unweighted McKay graph produces log-ratios
   that are MONOTONIC in graph distance but have the WRONG values.
   Max dynamic range: limited by spectral gap ~ {evals_L[1]:.3f}
   Needed: log-ratio up to 26 (top quark).
   At achievable t: max log-ratio ~ {max_lr:.1f}

   VERDICT: Unweighted heat kernel CANNOT reproduce the mass hierarchy.
   The graph is too small (9 nodes) and too symmetric.

2. WEIGHTED HEAT KERNEL (with Wilson line weights):
   Adding edge weights helps but the mapping node->fermion is unknown.
   Best scores are still large (poor fit).

   VERDICT: CIRCULAR - if we know the weights, we already know the masses.
   The heat kernel doesn't add information beyond the weights themselves.

3. GREEN'S FUNCTION:
   G = L^+ gives effective resistances = graph distances.
   On a tree, this is just the path length.
   Log-ratios don't match mass exponents.

   VERDICT: Too simple. Graph distances {0,1,2,3,4,5,6,6,7} don't
   match mass exponents {{0,3,5,11,11,16,17,19,26}}.

4. THE FUNDAMENTAL PROBLEM:
   The McKay graph has only 9 nodes and 8 edges.
   The mass exponents range from 0 to 26.
   Any diffusion process on 9 nodes equilibrates quickly,
   making it impossible to generate ratios spanning phi^26 ~ 271,000.

   The Wilson line approach WORKS because it uses CUMULATIVE path sums
   with specific edge weights, not diffusion. The heat kernel tries to
   extract the same information from spectral decomposition, but the
   9-node graph doesn't have enough spectral resolution.

CONCLUSION:
===========
The heat kernel on the 9-node McKay graph CANNOT reproduce the fermion
mass hierarchy. The graph is too small. This direction appears CLOSED.

The Wilson line approach (path integrals with edge weights) remains the
correct formulation for masses. The open problem is not the transport
mechanism (Wilson line vs heat kernel vs D_F eigenvalues) but the
DERIVATION of the edge weights [3, 2, 6, 0, 5, ?, ?, ?] from
first principles.
""")

print("=" * 70)
print("EXP-491 COMPLETE")
print("=" * 70)
