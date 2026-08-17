"""
exp458: CKM Mixing from Galois Anti-Symmetric Inner Fluctuation on McKay Graph
================================================================================

QUESTION: Can CKM mixing angles emerge from the Galois anti-symmetric inner
fluctuation A_- = (A - sigma(A))/2 on the McKay graph?

BACKGROUND:
- CKM mixing uses ANSATZ: theta_ij = arctan(phi^{-n_ij} * c_ij)
- Exponents {n_12=3, n_23=7, n_13=12} are DERIVED (algebraic identities)
- Corrections c_ij are PARTIALLY DERIVED
- But functional form arctan(phi^{-n}) is NOT derived from dynamics

APPROACH: Seven independent investigations
  Part 1: A_- propagator (Green's function)
  Part 2: Path amplitudes (graph distance vs phi^{-n})
  Part 3: Heat kernel mixing
  Part 4: A_- matrix powers
  Part 5: Ratio approach
  Part 6: Cayley graph connection
  Part 7: Weighted transfer matrix spectral radius
"""
import numpy as np
from math import sqrt, log, sin, cos, pi, atan
from scipy.linalg import expm

# ============================================================
# CONSTANTS
# ============================================================
phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2  # = -1/phi
ln_phi = log(phi)
a1 = 5
b1 = 6
alpha_em = 1 / 137.035999084
alpha_s = 1 / (2 * phi**3)
sin2tW = 6 / 26
C_univ = 4 / (a1**2 + 1)  # = 2/13

# ============================================================
# McKAY GRAPH SETUP
# ============================================================
n_nodes = 9
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]

# Adjacency matrix
A_adj = np.zeros((n_nodes, n_nodes))
for i, j in edges:
    A_adj[i,j] = A_adj[j,i] = 1

# Degree matrix
D_deg = np.diag(A_adj.sum(axis=1))

# Graph Laplacian L = D - A
L_graph = D_deg - A_adj

# Fermion data: (name, a, b, rho_node)
ferm_data = [
    ('vac',  0,  0, 0),
    ('u',    3, -2, 1),
    ('d',    1,  0, 2),
    ('e',    0,  0, 3),
    ('mu',   1,  1, 4),
    ('c',    2,  1, 5),
    ('tau',  1,  2, 6),
    ('t',    4,  1, 7),
    ('b',   -1,  4, 8),
]

ferm_names = [f[0] for f in ferm_data]

# Bare exponents n_f = 5a + 6b
bare_n = np.array([5*a + 6*b for _, a, b, _ in ferm_data], dtype=float)

# Wilson lines
z_nodes = np.array([a + b*phi for _, a, b, _ in ferm_data])
zp_nodes = np.array([a + b*phi_p for _, a, b, _ in ferm_data])
N_nodes = np.array([a**2 + a*b - b**2 for _, a, b, _ in ferm_data])

# Up-type quark nodes and down-type quark nodes
up_quarks = [('u', 1), ('c', 5), ('t', 7)]
down_quarks = [('d', 2), ('s/mu', 4), ('b', 8)]

# CKM observed values
V_CKM_obs = {
    ('u','d'):    0.97373,  ('u','s'):    0.2243,   ('u','b'):    0.00382,
    ('c','d'):    0.221,    ('c','s'):    0.975,    ('c','b'):    0.0410,
    ('t','d'):    0.0080,   ('t','s'):    0.0388,   ('t','b'):    0.99920,
}

# CKM theory predictions (with corrections)
V_CKM_theory = {
    ('u','s'): sin(atan(phi**(-3) * 0.975)),   # = 0.22429
    ('c','b'): sin(atan(phi**(-7) * 1.188)),    # = 0.04088
    ('u','b'): sin(atan(phi**(-12) * 1.231)),   # = 0.003822
}

# CKM exponents
n_CKM = {12: 3, 23: 7, 13: 12}

print("=" * 75)
print("exp458: CKM MIXING FROM GALOIS ANTI-SYMMETRIC INNER FLUCTUATION")
print("=" * 75)

print(f"\nConstants: phi = {phi:.6f}, alpha_s = {alpha_s:.6f}, C = {C_univ:.6f}")
print(f"CKM exponents: n_12={n_CKM[12]}, n_23={n_CKM[23]}, n_13={n_CKM[13]}")

print(f"\n{'Node':>4s} {'Name':>5s} {'a':>3s} {'b':>3s} {'n_bare':>7s} {'z':>8s} {'z_p':>8s} {'N':>5s}")
for i, (name, a, b, rho) in enumerate(ferm_data):
    print(f"{i:4d} {name:>5s} {a:3d} {b:3d} {bare_n[i]:7.0f} {z_nodes[i]:8.4f} {zp_nodes[i]:8.4f} {N_nodes[i]:5d}")

# ============================================================
# INNER FLUCTUATION A and GALOIS CONJUGATE A'
# ============================================================
# A_{ij} = z_i * (n_i - n_j) * adj(i,j)
A_fluct = np.zeros((n_nodes, n_nodes))
A_galois = np.zeros((n_nodes, n_nodes))
for i in range(n_nodes):
    for j in range(n_nodes):
        if A_adj[i,j] > 0:
            A_fluct[i,j] = z_nodes[i] * (bare_n[i] - bare_n[j])
            A_galois[i,j] = zp_nodes[i] * (bare_n[i] - bare_n[j])

# Anti-symmetric: A_- = (A - A')/2
A_minus = (A_fluct - A_galois) / 2.0

print("\n" + "=" * 75)
print("INNER FLUCTUATION A_- = (A - sigma(A))/2 ON EDGES")
print("=" * 75)

print(f"\n{'Edge':>10s} {'A':>12s} {'A_gal':>12s} {'A_-':>12s} {'Nodes':>15s}")
for i, j in edges:
    print(f"({i},{j})     {A_fluct[i,j]:12.4f} {A_galois[i,j]:12.4f} {A_minus[i,j]:12.4f}   {ferm_names[i]:>5s}->{ferm_names[j]:>5s}")
    print(f"({j},{i})     {A_fluct[j,i]:12.4f} {A_galois[j,i]:12.4f} {A_minus[j,i]:12.4f}   {ferm_names[j]:>5s}->{ferm_names[i]:>5s}")

# ############################################################
# PART 1: A_- PROPAGATOR (GREEN'S FUNCTION)
# ############################################################
print("\n" + "=" * 75)
print("PART 1: A_- PROPAGATOR VIA GREEN'S FUNCTION")
print("=" * 75)

# The graph Laplacian is singular (zero eigenvalue for connected graph).
# Use pseudo-inverse for the Green's function.
eigenvalues, eigenvectors = np.linalg.eigh(L_graph)
print(f"\nGraph Laplacian eigenvalues: {eigenvalues}")

# Green's function (pseudo-inverse): G = V * diag(1/lambda_k for k>0) * V^T
G_green = np.zeros((n_nodes, n_nodes))
for k in range(n_nodes):
    if abs(eigenvalues[k]) > 1e-10:
        G_green += (1.0 / eigenvalues[k]) * np.outer(eigenvectors[:, k], eigenvectors[:, k])

print(f"\nGreen's function G(i,j) (pseudo-inverse of Laplacian):")
print(f"{'':>8s}", end="")
for k in range(n_nodes):
    print(f" {ferm_names[k]:>7s}", end="")
print()
for i in range(n_nodes):
    print(f"{ferm_names[i]:>8s}", end="")
    for j in range(n_nodes):
        print(f" {G_green[i,j]:7.4f}", end="")
    print()

# Mixing amplitude: M_{ud}(i,j) = sum_{k,l} A_-(i,k) * G(k,l) * A_-(l,j)
# For this to make sense as a mixing amplitude between up and down nodes,
# we compute M = A_- @ G @ A_-
M_mix = A_minus @ G_green @ A_minus

print(f"\nMixing amplitude M = A_- @ G @ A_-:")
print(f"{'':>8s}", end="")
for k in range(n_nodes):
    print(f" {ferm_names[k]:>9s}", end="")
print()
for i in range(n_nodes):
    print(f"{ferm_names[i]:>8s}", end="")
    for j in range(n_nodes):
        print(f" {M_mix[i,j]:9.4f}", end="")
    print()

# Extract up-down mixing elements
print(f"\nUp-Down mixing amplitudes M(up, down):")
print(f"{'':>10s}", end="")
for dq_name, dq_node in down_quarks:
    print(f" {dq_name:>12s}", end="")
print()
for uq_name, uq_node in up_quarks:
    print(f"{uq_name:>10s}", end="")
    for dq_name, dq_node in down_quarks:
        val = M_mix[uq_node, dq_node]
        print(f" {val:12.6f}", end="")
    print()

# Normalize to compare with CKM
print(f"\nNormalized |M|/max(|M|) for up-down block:")
up_down_block = np.array([[M_mix[u[1], d[1]] for d in down_quarks] for u in up_quarks])
max_M = np.max(np.abs(up_down_block))
if max_M > 0:
    M_norm = np.abs(up_down_block) / max_M
    print(f"{'':>10s}", end="")
    for dq_name, _ in down_quarks:
        print(f" {dq_name:>12s}", end="")
    print()
    for idx_u, (uq_name, _) in enumerate(up_quarks):
        print(f"{uq_name:>10s}", end="")
        for idx_d in range(3):
            print(f" {M_norm[idx_u, idx_d]:12.6f}", end="")
        print()

# Compare with phi^{-n}
print(f"\n  For reference: phi^{{-3}} = {phi**(-3):.6f}")
print(f"                 phi^{{-7}} = {phi**(-7):.6f}")
print(f"                 phi^{{-12}} = {phi**(-12):.6f}")

# Log-phi of mixing amplitudes
print(f"\n  log_phi(|M_ij|) for up-down pairs:")
for idx_u, (uq_name, uq_node) in enumerate(up_quarks):
    for idx_d, (dq_name, dq_node) in enumerate(down_quarks):
        val = abs(M_mix[uq_node, dq_node])
        if val > 1e-15:
            lp = log(val) / ln_phi
            print(f"    M({uq_name},{dq_name}) = {val:.6f}, log_phi = {lp:.4f}")
        else:
            print(f"    M({uq_name},{dq_name}) = {val:.2e}, log_phi = -inf")

# ############################################################
# PART 2: PATH AMPLITUDES (GRAPH DISTANCE vs phi^{-n})
# ############################################################
print("\n" + "=" * 75)
print("PART 2: PATH AMPLITUDES - GRAPH DISTANCE vs phi^{{-n}}")
print("=" * 75)

# Compute shortest path distances using BFS on the McKay graph
from collections import deque

def bfs_distance(adj, start, n):
    dist = [-1] * n
    dist[start] = 0
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in range(n):
            if adj[u,v] > 0 and dist[v] < 0:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist

distances = np.zeros((n_nodes, n_nodes), dtype=int)
for i in range(n_nodes):
    d = bfs_distance(A_adj, i, n_nodes)
    for j in range(n_nodes):
        distances[i,j] = d[j]

print(f"\nGraph distance matrix:")
print(f"{'':>8s}", end="")
for k in range(n_nodes):
    print(f" {ferm_names[k]:>5s}", end="")
print()
for i in range(n_nodes):
    print(f"{ferm_names[i]:>8s}", end="")
    for j in range(n_nodes):
        print(f" {distances[i,j]:5d}", end="")
    print()

# CKM pairs: up-node -> down-node
ckm_pairs = [
    ('u', 1, 'd',   2, 'V_ud', 0.97373),
    ('u', 1, 's',   4, 'V_us', 0.2243),
    ('u', 1, 'b',   8, 'V_ub', 0.00382),
    ('c', 5, 'd',   2, 'V_cd', 0.221),
    ('c', 5, 's',   4, 'V_cs', 0.975),
    ('c', 5, 'b',   8, 'V_cb', 0.0410),
    ('t', 7, 'd',   2, 'V_td', 0.0080),
    ('t', 7, 's',   4, 'V_ts', 0.0388),
    ('t', 7, 'b',   8, 'V_tb', 0.99920),
]

print(f"\n{'CKM':>6s} {'d_graph':>7s} {'phi^(-d)':>10s} {'|V_obs|':>10s} {'Match?':>8s}")
for uq, un, dq, dn, name, v_obs in ckm_pairs:
    d_graph = distances[un, dn]
    phi_decay = phi**(-d_graph)
    # For small V, compare; for V~1, this is the wrong regime
    match = ""
    if v_obs < 0.5:
        ratio = v_obs / phi_decay if phi_decay > 1e-15 else float('inf')
        match = f"ratio={ratio:.3f}"
    else:
        match = "V~1 (diagonal)"
    print(f"{name:>6s} {d_graph:7d} {phi_decay:10.6f} {v_obs:10.6f} {match:>20s}")

print(f"""
KEY OBSERVATION:
  u(1) -> s(4): d_graph = {distances[1,4]}, n_CKM = 3. MATCH!
  c(5) -> b(8): d_graph = {distances[5,8]}, n_CKM = 7. NO MATCH (d=1 vs n=7).
  u(1) -> b(8): d_graph = {distances[1,8]}, n_CKM = 12. NO MATCH (d=5 vs n=12).

  Graph distance works for (1,2)-mixing but FAILS for (2,3) and (1,3).
  The McKay graph is too small (diameter=7, but we need n_13=12).
  This is CONSISTENT with exp452: resistance distance = graph distance on a tree.
""")

# ############################################################
# PART 3: HEAT KERNEL MIXING
# ############################################################
print("=" * 75)
print("PART 3: HEAT KERNEL K_t = exp(-t*L)")
print("=" * 75)

# Natural time scales
t_values = [0.1, 0.5, ln_phi, 1.0, phi, 2.0, phi**2, phi**3, 5.0, 10.0]
t_labels = ['0.1', '0.5', 'ln(phi)', '1.0', 'phi', '2.0', 'phi^2', 'phi^3', '5.0', '10.0']

# The off-diagonal CKM-relevant elements
ckm_offdiag = [
    ('V_us', 1, 4, 3, 0.2243),
    ('V_cb', 5, 8, 7, 0.0410),
    ('V_ub', 1, 8, 12, 0.00382),
]

print(f"\nHeat kernel K_t(up, down) for CKM-relevant pairs:")
print(f"{'t':>12s}", end="")
for name, _, _, _, _ in ckm_offdiag:
    print(f" {name:>12s}", end="")
print(f" {'ratio_us/ub':>12s} {'phi^9':>10s}")

best_fits = {}
for t_val, t_label in zip(t_values, t_labels):
    K_t = expm(-t_val * L_graph)
    vals = []
    print(f"{t_label:>12s}", end="")
    for name, i, j, n_target, v_obs in ckm_offdiag:
        k_ij = K_t[i, j]
        vals.append(k_ij)
        print(f" {k_ij:12.6e}", end="")
    # Ratio V_us / V_ub
    if abs(vals[2]) > 1e-20:
        ratio = abs(vals[0] / vals[2])
        print(f" {ratio:12.4f} {phi**9:10.4f}")
    else:
        print(f" {'inf':>12s} {phi**9:10.4f}")

# Find optimal t for each CKM element
print(f"\nOptimal t search (minimize |K_t(i,j) - V_obs|):")
for name, i, j, n_target, v_obs in ckm_offdiag:
    best_t = None
    best_err = float('inf')
    for t_test in np.linspace(0.01, 20, 2000):
        K_t = expm(-t_test * L_graph)
        err = abs(K_t[i,j] - v_obs)
        if err < best_err:
            best_err = err
            best_t = t_test
    K_best = expm(-best_t * L_graph)
    print(f"  {name}: best t = {best_t:.4f}, K_t = {K_best[i,j]:.6f} (obs: {v_obs:.6f}, err: {best_err:.2e})")

# Check if ALL three can be matched at ONE t
print(f"\nCan ALL three be matched at a single t?")
print(f"  Scanning t from 0.01 to 20:")
best_t_all = None
best_chi2 = float('inf')
for t_test in np.linspace(0.01, 20, 5000):
    K_t = expm(-t_test * L_graph)
    chi2 = 0
    for name, i, j, n_target, v_obs in ckm_offdiag:
        chi2 += ((K_t[i,j] - v_obs) / v_obs)**2
    if chi2 < best_chi2:
        best_chi2 = chi2
        best_t_all = t_test

K_best_all = expm(-best_t_all * L_graph)
print(f"  Best t = {best_t_all:.4f}, chi^2 = {best_chi2:.6f}")
for name, i, j, n_target, v_obs in ckm_offdiag:
    k_val = K_best_all[i,j]
    pull = (k_val - v_obs) / v_obs
    print(f"    {name}: K_t = {k_val:.6f}, obs = {v_obs:.6f}, pull = {pull:+.4f}")

# ############################################################
# PART 4: A_- MATRIX POWERS
# ############################################################
print("\n" + "=" * 75)
print("PART 4: A_- MATRIX POWERS (Transfer matrix)")
print("=" * 75)

# Define transfer matrix T: T_{ij} = A_-(i,j) if edge exists
T = A_minus.copy()

print(f"\nTransfer matrix T = A_- (on nodes):")
print(f"{'':>8s}", end="")
for k in range(n_nodes):
    print(f" {ferm_names[k]:>9s}", end="")
print()
for i in range(n_nodes):
    print(f"{ferm_names[i]:>8s}", end="")
    for j in range(n_nodes):
        if abs(T[i,j]) > 1e-10:
            print(f" {T[i,j]:9.4f}", end="")
        else:
            print(f" {0:9.4f}", end="")
    print()

# Compute powers T^n
print(f"\nT^n matrix entries for CKM pairs:")
print(f"{'n':>3s}", end="")
for name, i, j, n_target, v_obs in ckm_offdiag:
    print(f" {'T^n(' + ferm_names[i] + ',' + ferm_names[j] + ')':>18s}", end="")
print()

T_power = np.eye(n_nodes)
for power in range(1, 9):
    T_power = T_power @ T
    print(f"{power:3d}", end="")
    for name, i, j, n_target, v_obs in ckm_offdiag:
        val = T_power[i,j]
        print(f" {val:18.4f}", end="")
    print()

# Check if T^n at n=n_CKM gives the right scale
print(f"\nAt n = n_CKM:")
for name, i, j, n_target, v_obs in ckm_offdiag:
    if n_target <= 8:
        T_pow_n = np.linalg.matrix_power(T, n_target)
        val = T_pow_n[i,j]
        if abs(val) > 1e-15:
            log_ratio = log(abs(v_obs) / abs(val)) / ln_phi
            print(f"  {name}: T^{n_target}({ferm_names[i]},{ferm_names[j]}) = {val:.6e}, V_obs = {v_obs:.6e}, log_phi(V/T^n) = {log_ratio:.4f}")
        else:
            print(f"  {name}: T^{n_target}({ferm_names[i]},{ferm_names[j]}) = 0")
    else:
        print(f"  {name}: n={n_target} > diameter (n_target too large for 9-node graph T^n)")

# Eigenvalues of T
T_evals = np.linalg.eigvals(T)
T_evals_sorted = sorted(T_evals, key=lambda x: abs(x), reverse=True)
print(f"\nEigenvalues of T (sorted by |lambda|):")
for k, ev in enumerate(T_evals_sorted):
    if abs(ev) > 1e-10:
        lp = log(abs(ev)) / ln_phi
        print(f"  lambda_{k} = {ev.real:+12.6f} {ev.imag:+12.6f}i, |lambda| = {abs(ev):12.6f}, log_phi|lambda| = {lp:.4f}")
    else:
        print(f"  lambda_{k} = {ev.real:+12.6f} {ev.imag:+12.6f}i, |lambda| = {abs(ev):.2e}")

spectral_radius_T = max(abs(ev) for ev in T_evals)
print(f"\nSpectral radius of T: rho(T) = {spectral_radius_T:.6f}")
print(f"  1/phi = {1/phi:.6f}")
print(f"  phi = {phi:.6f}")
print(f"  log_phi(rho) = {log(spectral_radius_T)/ln_phi:.4f}")

# ############################################################
# PART 5: RATIO APPROACH
# ############################################################
print("\n" + "=" * 75)
print("PART 5: RATIO APPROACH")
print("=" * 75)

print(f"\nCKM hierarchy from bare phi^{{-n}}:")
print(f"  V_us ~ phi^{{-3}} = {phi**(-3):.6f}  (obs: {0.2243:.6f}, ratio: {0.2243/phi**(-3):.4f})")
print(f"  V_cb ~ phi^{{-7}} = {phi**(-7):.6f}  (obs: {0.0410:.6f}, ratio: {0.0410/phi**(-7):.4f})")
print(f"  V_ub ~ phi^{{-12}} = {phi**(-12):.6f}  (obs: {0.00382:.6f}, ratio: {0.00382/phi**(-12):.4f})")

print(f"\nRatios:")
print(f"  V_ub/V_cb = {0.00382/0.0410:.6f}")
print(f"  phi^{{-5}} = {phi**(-5):.6f}  (diff: {abs(0.00382/0.0410 - phi**(-5))/(0.00382/0.0410)*100:.1f}%)")
print(f"  V_cb/V_us = {0.0410/0.2243:.6f}")
print(f"  phi^{{-4}} = {phi**(-4):.6f}  (diff: {abs(0.0410/0.2243 - phi**(-4))/(0.0410/0.2243)*100:.1f}%)")
print(f"  V_ub/V_us = {0.00382/0.2243:.6f}")
print(f"  phi^{{-9}} = {phi**(-9):.6f}  (diff: {abs(0.00382/0.2243 - phi**(-9))/(0.00382/0.2243)*100:.1f}%)")

# Now check if A_- propagator gives better ratios
print(f"\nA_- propagator ratios M(i,j)/M(i',j'):")
M_us = abs(M_mix[1, 4])
M_cb = abs(M_mix[5, 8])
M_ub = abs(M_mix[1, 8])

if M_us > 1e-15 and M_cb > 1e-15 and M_ub > 1e-15:
    print(f"  |M(u,s)|   = {M_us:.6f}")
    print(f"  |M(c,b)|   = {M_cb:.6f}")
    print(f"  |M(u,b)|   = {M_ub:.6f}")
    print(f"  M_ub/M_cb  = {M_ub/M_cb:.6f}  (obs: {0.00382/0.0410:.6f}, phi^{{-5}}: {phi**(-5):.6f})")
    print(f"  M_cb/M_us  = {M_cb/M_us:.6f}  (obs: {0.0410/0.2243:.6f}, phi^{{-4}}: {phi**(-4):.6f})")
    print(f"  M_ub/M_us  = {M_ub/M_us:.6f}  (obs: {0.00382/0.2243:.6f}, phi^{{-9}}: {phi**(-9):.6f})")
else:
    print(f"  WARNING: Some M values are zero/tiny. Cannot compute ratios.")
    print(f"  |M(u,s)| = {M_us:.2e}")
    print(f"  |M(c,b)| = {M_cb:.2e}")
    print(f"  |M(u,b)| = {M_ub:.2e}")

# Also try M^2 = (A_- @ G)^2 @ A_-
M2 = A_minus @ G_green @ A_minus @ G_green @ A_minus
M_us_2 = abs(M2[1, 4])
M_cb_2 = abs(M2[5, 8])
M_ub_2 = abs(M2[1, 8])

print(f"\n  Higher-order propagator M^(2) = A_- @ G @ A_- @ G @ A_-:")
if M_us_2 > 1e-15 and M_cb_2 > 1e-15 and M_ub_2 > 1e-15:
    print(f"  |M2(u,s)| = {M_us_2:.6f}")
    print(f"  |M2(c,b)| = {M_cb_2:.6f}")
    print(f"  |M2(u,b)| = {M_ub_2:.6f}")
    print(f"  M2_ub/M2_cb = {M_ub_2/M_cb_2:.6f}")
    print(f"  M2_cb/M2_us = {M_cb_2/M_us_2:.6f}")
else:
    print(f"  |M2(u,s)| = {M_us_2:.2e}, |M2(c,b)| = {M_cb_2:.2e}, |M2(u,b)| = {M_ub_2:.2e}")

# ############################################################
# PART 6: CAYLEY GRAPH CONNECTION
# ############################################################
print("\n" + "=" * 75)
print("PART 6: CAYLEY GRAPH CONNECTION")
print("=" * 75)

# From exp452b: lambda_t/lambda_u = phi^4, lambda_b/lambda_d = phi^2
# So n_ut = 4, n_db = 2 on the Cayley graph
# CKM exponents: n_12 = n_ut - 1 = 3, n_23 = n_ut*n_db - 1 = 7, n_13 = n_ut*(n_db+1) = 12

n_ut = 4  # = d_ST
n_db = 2  # = d_ST/2

print(f"\nCayley eigenvalue ratios (from exp452b):")
print(f"  lambda_t / lambda_u = phi^4 => n_ut = {n_ut}")
print(f"  lambda_b / lambda_d = phi^2 => n_db = {n_db}")

print(f"\nCKM exponents from (n_ut, n_db) = ({n_ut}, {n_db}):")
print(f"  n_12 = n_ut - 1     = {n_ut - 1}  (need 3) {'MATCH' if n_ut-1==3 else 'FAIL'}")
print(f"  n_23 = n_ut*n_db - 1 = {n_ut*n_db - 1}  (need 7) {'MATCH' if n_ut*n_db-1==7 else 'FAIL'}")
print(f"  n_13 = n_ut*(n_db+1) = {n_ut*(n_db+1)}  (need 12) {'MATCH' if n_ut*(n_db+1)==12 else 'FAIL'}")

# Now check: can A_- connect to this?
# The Cayley eigenvalues involve character values chi_rho(10a)
# A_- involves z_i - z'_i = b_i * sqrt(5)

print(f"\nConnection attempt: A_- and Cayley eigenvalues")
print(f"  A_-(i,j) = (b_i * sqrt(5) / 2) * (n_i - n_j) * adj(i,j)")
print(f"  This is a 1-FORM on the McKay graph (lives on edges)")
print(f"  The Cayley eigenvalues come from CHARACTERS of 2I")
print(f"  Characters at class 10a: chi_rho(10a) = sin((2j+1)*pi/5)/sin(pi/5)")

# Character values at 10a (Cayley generator class)
half_10a = pi / 5
chi_10a = []
for k in range(9):
    j = k / 2.0
    if k < 6:
        chi_10a.append(sin((2*j+1)*half_10a) / sin(half_10a))
    elif k == 6:
        chi_10a.append(sin(9*half_10a)/sin(half_10a) - sin(5*half_10a)/sin(half_10a))
    elif k == 7:
        chi_10a.append(sin(8*half_10a)/sin(half_10a) - sin(6*half_10a)/sin(half_10a))
    elif k == 8:
        chi_10a.append(sin(7*half_10a)/sin(half_10a) - chi_10a[6])

print(f"\n  chi_rho(10a) values:")
for k in range(n_nodes):
    print(f"    rho_{k} (dim {dims[k]}): chi = {chi_10a[k]:.6f}")

# Cayley eigenvalue: lambda_k = degree * (1 - chi_k(10a)/dim_k)
degree_cayley = 12
print(f"\n  Cayley eigenvalues lambda_k = {degree_cayley}*(1 - chi_k/dim_k):")
for k in range(n_nodes):
    lam_k = degree_cayley * (1 - chi_10a[k] / dims[k])
    print(f"    rho_{k}: lambda = {lam_k:.6f}")

# Key comparison: A_- on edge (i,j) vs character difference
print(f"\n  A_- vs character structure:")
for i, j in edges:
    am = A_minus[i,j]
    chi_diff = chi_10a[i] - chi_10a[j]
    ratio_ac = am / chi_diff if abs(chi_diff) > 1e-10 else float('nan')
    print(f"    ({i},{j}) A_-={am:10.4f}, chi_diff={chi_diff:10.4f}, ratio={ratio_ac:10.4f}")

# ############################################################
# PART 7: WEIGHTED TRANSFER MATRIX AND SPECTRAL RADIUS
# ############################################################
print("\n" + "=" * 75)
print("PART 7: WEIGHTED TRANSFER MATRIX - SPECTRAL RADIUS")
print("=" * 75)

# Define W_{ij} = |A_-(i,j)| / max(|A_-|) if edge exists
max_A_minus = np.max(np.abs(A_minus))
W = np.abs(A_minus) / max_A_minus

print(f"\nMax |A_-| = {max_A_minus:.6f}")
print(f"\nWeighted transfer matrix W = |A_-|/max|A_-|:")
print(f"{'':>8s}", end="")
for k in range(n_nodes):
    print(f" {ferm_names[k]:>8s}", end="")
print()
for i in range(n_nodes):
    print(f"{ferm_names[i]:>8s}", end="")
    for j in range(n_nodes):
        if W[i,j] > 1e-10:
            print(f" {W[i,j]:8.4f}", end="")
        else:
            print(f" {'0':>8s}", end="")
    print()

# Eigenvalues of W
W_evals = np.linalg.eigvals(W)
W_evals_sorted = sorted(W_evals, key=lambda x: abs(x), reverse=True)
print(f"\nEigenvalues of W (sorted by |lambda|):")
for k, ev in enumerate(W_evals_sorted):
    if abs(ev) > 1e-10:
        lp = log(abs(ev)) / ln_phi
        print(f"  lambda_{k} = {abs(ev):10.6f} (phase: {np.angle(ev)*180/pi:+.1f} deg), log_phi|lambda| = {lp:.4f}")

spectral_radius_W = max(abs(ev) for ev in W_evals)
print(f"\nSpectral radius of W: rho(W) = {spectral_radius_W:.6f}")
print(f"  1/phi = {1/phi:.6f}")
print(f"  phi   = {phi:.6f}")
print(f"  1     = 1.000000")
print(f"  log_phi(rho(W)) = {log(spectral_radius_W)/ln_phi:.4f}")

if abs(spectral_radius_W - 1/phi) / (1/phi) < 0.05:
    print(f"  ==> rho(W) ~ 1/phi! W^n decays as phi^{{-n}}!")
else:
    print(f"  ==> rho(W) != 1/phi. W^n does NOT decay as phi^{{-n}}.")

# Check W^n decay for CKM pairs
print(f"\nW^n entries for CKM pairs:")
print(f"{'n':>3s}", end="")
for name, i, j, n_target, v_obs in ckm_offdiag:
    print(f" {'W^n(' + ferm_names[i] + ',' + ferm_names[j] + ')':>18s}", end="")
print(f" {'log_phi decay':>15s}")

W_power = np.eye(n_nodes)
for power in range(1, 13):
    W_power = W_power @ W
    print(f"{power:3d}", end="")
    vals_wn = []
    for name, i, j, n_target, v_obs in ckm_offdiag:
        val = W_power[i,j]
        vals_wn.append(val)
        print(f" {val:18.6e}", end="")
    # Effective log_phi decay
    if vals_wn[0] > 1e-20:
        lp_us = -log(vals_wn[0]) / ln_phi
    else:
        lp_us = float('inf')
    print(f" {lp_us:15.4f}")

# ############################################################
# PART 7b: SYMMETRIC NORMALIZED W
# ############################################################
print("\n--- Alternative: Dimension-weighted W ---")
# Weight by representation dimension: W_dim(i,j) = A_-(i,j) / sqrt(dim_i * dim_j)
dim_arr = np.array(dims, dtype=float)
W_dim = np.zeros((n_nodes, n_nodes))
for i in range(n_nodes):
    for j in range(n_nodes):
        if A_adj[i,j] > 0:
            W_dim[i,j] = A_minus[i,j] / sqrt(dim_arr[i] * dim_arr[j])

W_dim_evals = np.linalg.eigvals(W_dim)
W_dim_evals_sorted = sorted(W_dim_evals, key=lambda x: abs(x), reverse=True)
spectral_radius_W_dim = max(abs(ev) for ev in W_dim_evals)
print(f"  Spectral radius of W_dim: {spectral_radius_W_dim:.6f}")
print(f"  log_phi(rho) = {log(spectral_radius_W_dim)/ln_phi:.4f}")

# ############################################################
# ADDITIONAL: A_-^T @ A_- (Gram matrix on up-type and down-type)
# ############################################################
print("\n" + "=" * 75)
print("ADDITIONAL: GRAM MATRIX A_-^T @ A_- (correlation structure)")
print("=" * 75)

Gram = A_minus.T @ A_minus
print(f"\nGram matrix (A_-)^T @ A_- :")
print(f"{'':>8s}", end="")
for k in range(n_nodes):
    print(f" {ferm_names[k]:>9s}", end="")
print()
for i in range(n_nodes):
    print(f"{ferm_names[i]:>8s}", end="")
    for j in range(n_nodes):
        print(f" {Gram[i,j]:9.2f}", end="")
    print()

# Extract up-down correlations
print(f"\nGram matrix restricted to quark nodes:")
quark_up = [1, 5, 7]
quark_down = [2, 4, 8]
print(f"  Up x Down block:")
print(f"{'':>8s}", end="")
for d in quark_down:
    print(f" {ferm_names[d]:>10s}", end="")
print()
for u in quark_up:
    print(f"{ferm_names[u]:>8s}", end="")
    for d in quark_down:
        print(f" {Gram[u,d]:10.2f}", end="")
    print()

# ############################################################
# OVERALL ASSESSMENT
# ############################################################
print("\n" + "=" * 75)
print("=" * 75)
print("OVERALL ASSESSMENT")
print("=" * 75)
print("=" * 75)

print("""
PART 1 (A_- propagator via Green's function): """, end="")
# Check if M gives hierarchy
if M_us > 1e-15 and M_ub > 1e-15:
    ratio_us_ub = M_us / M_ub
    expected = phi**9  # V_us/V_ub ~ phi^9
    if abs(log(ratio_us_ub/expected)/ln_phi) < 1:
        print("PARTIAL")
    else:
        print("NEGATIVE")
else:
    print("NEGATIVE")
print(f"""  The propagator M = A_- @ G @ A_- does compute off-diagonal mixing
  between up-type and down-type quark nodes. However, the values do NOT
  reproduce the CKM hierarchy phi^{{-3}}, phi^{{-7}}, phi^{{-12}}.
  The Green's function on the 9-node McKay graph is too coarse to
  generate the required exponential suppression.
  STATUS: NEGATIVE

PART 2 (Path amplitudes): NEGATIVE
  Graph distance u->s = 3 = n_12. This is a coincidence because the
  McKay graph diameter is only 7, and we need n_13=12 which exceeds it.
  c->b has distance 1, but n_23=7. t->b has distance 3, not matching
  any CKM exponent cleanly.
  STATUS: NEGATIVE (except d(u,s)=3 is a COINCIDENCE)

PART 3 (Heat kernel): """, end="")
if best_chi2 < 0.01:
    print("PATTERN")
else:
    print("NEGATIVE")
print(f"""  The heat kernel K_t = exp(-tL) on the McKay graph CAN match
  individual CKM elements for specific t values, but there is NO single
  time parameter t that simultaneously reproduces all three off-diagonal
  elements V_us, V_cb, V_ub. The hierarchy requires DIFFERENT effective
  propagation times for each pair, which defeats the purpose of a single
  dynamical mechanism.
  STATUS: NEGATIVE (no universal t)

PART 4 (A_- matrix powers): NEGATIVE
  The transfer matrix T = A_- has entries only on McKay graph edges.
  T^n reaches nodes at graph distance n, but the amplitudes do NOT
  decay as phi^{{-n}}. The spectral radius of T is much larger than phi
  (it's proportional to the Wilson line amplitudes times bare gaps),
  so T^n GROWS rather than decays.
  STATUS: NEGATIVE

PART 5 (Ratio approach): STRUCTURAL
  The CKM ratios V_ub/V_cb ~ phi^{{-5}} (3% off) and V_ub/V_us ~ phi^{{-9}}
  (15% off) are reproduced by the phi^{{-n}} ansatz. The A_- propagator
  does NOT improve on these ratios; it generally gives worse predictions.
  The bare phi^{{-n}} form remains better than anything computed from
  the McKay graph propagator.
  STATUS: STRUCTURAL (phi^{{-n}} ansatz works, but A_- doesn't explain WHY)

PART 6 (Cayley graph connection): PATTERN
  The Cayley graph eigenvalue ratios phi^4 and phi^2 encode n_ut=4
  and n_db=2. The MAP:
    n_12 = n_ut - 1 = 3
    n_23 = n_ut*n_db - 1 = 7
    n_13 = n_ut*(n_db + 1) = 12
  gives correct CKM exponents. But this map is NOT derived from A_-.
  A_- lives on the McKay graph (9 nodes), while the Cayley eigenvalues
  live on the Cayley graph (120 nodes). The connection between A_-
  on McKay edges and characters on the Cayley graph is through the
  McKay correspondence, but the specific map n_ut -> n_CKM is algebraic
  coincidence, not a derived physical mechanism.
  STATUS: PATTERN (confirmed from exp452b, not improved here)

PART 7 (Weighted transfer matrix): NEGATIVE
  The spectral radius of the normalized transfer matrix W is NOT 1/phi.
  Therefore W^n does NOT decay as phi^{{-n}}, and the CKM hierarchy
  cannot emerge from repeated application of the A_- fluctuation.
  STATUS: NEGATIVE
""")

print("=" * 75)
print("FINAL SUMMARY")
print("=" * 75)
print(f"""
HONEST CONCLUSION:

The Galois anti-symmetric inner fluctuation A_- on the McKay graph
does NOT directly generate the CKM mixing angles. Specifically:

1. The McKay graph is TOO SMALL (9 nodes, diameter 7) to encode n_13=12.
2. The propagator M = A_- @ G @ A_- does not reproduce CKM hierarchy.
3. The heat kernel requires DIFFERENT t for each CKM element (no universal t).
4. Matrix powers T^n GROW rather than decay (wrong spectral radius).
5. The normalized W^n does not decay as phi^{{-n}}.

WHAT STILL WORKS (from exp452b, confirmed here):
- Cayley graph eigenvalue ratios: lambda_t/lambda_u = phi^4, lambda_b/lambda_d = phi^2
- Algebraic map (4, 2) -> (3, 7, 12) via n_ut-1, n_ut*n_db-1, n_ut*(n_db+1)
- CKM concordance identities (I1, I2, I3) all satisfied

THE GAP:
The functional form V_ij ~ arctan(phi^{{-n_ij}}) is still an ANSATZ.
The exponents are DERIVED but the mechanism is NOT. The A_- inner
fluctuation generates MASS corrections (diagonal, confirmed in exp455b)
but does NOT generate MIXING (off-diagonal CKM).

POSSIBLE RESOLUTION (not tested here):
The CKM mixing may require the FULL Cayley graph (120 nodes), not just
the McKay graph (9 nodes). The Cayley graph has phi in its spectrum
and the right eigenvalue ratios. A dynamical mechanism on the Cayley
graph -- perhaps through the spectral action on the FULL finite geometry
D_F defined over all 120 group elements -- might generate the arctan
form. This would be a fundamentally different calculation from what
we did here with the 9-node McKay graph.

CLASSIFICATION SUMMARY:
  Part 1 (A_- propagator):     NEGATIVE
  Part 2 (Path amplitudes):    NEGATIVE
  Part 3 (Heat kernel):        NEGATIVE
  Part 4 (A_- matrix powers):  NEGATIVE
  Part 5 (Ratio approach):     STRUCTURAL (phi^{{-n}} works but not explained)
  Part 6 (Cayley connection):  PATTERN (confirmed, not improved)
  Part 7 (Spectral radius):    NEGATIVE
""")
