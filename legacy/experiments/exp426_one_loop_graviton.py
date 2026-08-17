"""
exp426_one_loop_graviton.py
===========================
One-loop graviton self-energy on the 600-cell.

Key insight: since all 600 tets are CONGRUENT (regular, edge=1/phi),
the per-tet derivative matrices M = d(theta)/dl and N = d^2(theta)/(dl dl)
are UNIVERSAL. Compute once, assemble 600 times.

Regge calculus identities:
  dS/dl_a = epsilon_a  (Schlaefli)
  d^2S/(dl_a dl_b) = -sum_T d(theta_a^T)/dl_b
  d^3S/(dl_a dl_b dl_c) = -sum_T d^2(theta_a^T)/(dl_b dl_c)
"""

import numpy as np
import itertools
from collections import defaultdict
import time

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 1.0 / 137.035999084

print("=" * 72)
print("EXP-426: ONE-LOOP GRAVITON SELF-ENERGY ON THE 600-CELL")
print("=" * 72)

# ============================================================
# PART 1: INFRASTRUCTURE
# ============================================================
print("\n" + "=" * 72)
print("PART 1: INFRASTRUCTURE")
print("=" * 72)
t0 = time.time()

# Build 600-cell
verts_set = set()
for i in range(4):
    for s in [1, -1]:
        v = [0,0,0,0]; v[i] = s
        verts_set.add(tuple(v))
for signs in itertools.product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = [p for p in itertools.permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in itertools.product([1,-1], repeat=len(nz)):
        v = list(base)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

verts = np.array(sorted(verts_set))
N_v = len(verts)
assert N_v == 120

dot_threshold = PHI / 2
dots = verts @ verts.T
edges = []
edge_to_idx = {}
adj_list = defaultdict(set)
for i in range(N_v):
    for j in range(i+1, N_v):
        if abs(dots[i,j] - dot_threshold) < 0.001:
            idx = len(edges)
            edges.append((i,j))
            edge_to_idx[(i,j)] = idx
            edge_to_idx[(j,i)] = idx
            adj_list[i].add(j)
            adj_list[j].add(i)
N_e = len(edges)
assert N_e == 720

edge_set = {(min(i,j), max(i,j)) for (i,j) in edges}
cells = set()
for (i,j) in edges:
    common = adj_list[i] & adj_list[j]
    for k in common:
        for l in common:
            if l > k and (min(k,l), max(k,l)) in edge_set:
                cells.add(tuple(sorted([i,j,k,l])))
cells = sorted(cells)
N_c = len(cells)
assert N_c == 600

triangles = []
for i in range(N_v):
    for j in adj_list[i]:
        if j > i:
            for k in (adj_list[i] & adj_list[j]):
                if k > j:
                    triangles.append((i,j,k))
N_f = len(triangles)
assert N_f == 1200

edge_to_cells = [[] for _ in range(N_e)]
for ci, cell in enumerate(cells):
    for pair in itertools.combinations(cell, 2):
        ei = edge_to_idx.get((pair[0], pair[1]))
        if ei is not None:
            edge_to_cells[ei].append(ci)

# Tet-to-edges map
tet_edges_list = []
for ci, cell in enumerate(cells):
    te = [edge_to_idx[(p[0],p[1])] for p in itertools.combinations(cell, 2)]
    tet_edges_list.append(te)

# Boundary operators
d0 = np.zeros((N_e, N_v))
for e_idx, (i,j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = 1.0

d1 = np.zeros((N_f, N_e))
for f_idx, (i,j,k) in enumerate(triangles):
    e_ij = edge_to_idx.get((i,j))
    e_ik = edge_to_idx.get((i,k))
    e_jk = edge_to_idx.get((j,k))
    if e_ij is not None and e_ik is not None and e_jk is not None:
        d1[f_idx, e_ij] = 1.0
        d1[f_idx, e_jk] = 1.0
        d1[f_idx, e_ik] = -1.0

assert np.max(np.abs(d1 @ d0)) < 1e-10

# Coexact Laplacian and Green's function
L_coexact = d1.T @ d1
evals_co, evecs_co = np.linalg.eigh(L_coexact)
coexact_mask = evals_co > 0.01
n_coexact = np.sum(coexact_mask)
V_co = evecs_co[:, coexact_mask]
lambda_co = evals_co[coexact_mask]
P_coexact = V_co @ V_co.T
G_matrix = V_co @ np.diag(1.0 / lambda_co) @ V_co.T

# Exact projector
L_exact = d0 @ d0.T
evals_ex, evecs_ex = np.linalg.eigh(L_exact)
V_ex = evecs_ex[:, evals_ex > 0.01]
P_exact = V_ex @ V_ex.T

mass_gap = lambda_co.min()
print(f"  V={N_v}, E={N_e}, F={N_f}, C={N_c}")
print(f"  Coexact DOF: {n_coexact}, Mass gap: {mass_gap:.6f} = 7-4*phi = {7-4*PHI:.6f}")
print(f"  Infrastructure: {time.time()-t0:.1f}s")


# ============================================================
# PART 2: PER-TET DIHEDRAL ANGLE DERIVATIVES (UNIVERSAL)
# ============================================================
print("\n" + "=" * 72)
print("PART 2: UNIVERSAL TET DERIVATIVE MATRICES M, N")
print("=" * 72)
t1 = time.time()

TET_EDGE_PAIRS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

def tet_reconstruct_3d(L):
    """Reconstruct 3D tet vertices from 6 edge lengths.
    L = [l01, l02, l03, l12, l13, l23]"""
    l01, l02, l03, l12, l13, l23 = L
    v0 = np.array([0.0, 0.0, 0.0])
    v1 = np.array([l01, 0.0, 0.0])
    x2 = (l01**2 + l02**2 - l12**2) / (2*l01)
    y2 = np.sqrt(max(l02**2 - x2**2, 0))
    v2 = np.array([x2, y2, 0.0])
    x3 = (l01**2 + l03**2 - l13**2) / (2*l01)
    dot23 = (l02**2 + l03**2 - l23**2) / 2
    y3 = (dot23 - x2*x3) / y2 if y2 > 1e-15 else 0.0
    z3 = np.sqrt(max(l03**2 - x3**2 - y3**2, 0))
    v3 = np.array([x3, y3, z3])
    return np.array([v0, v1, v2, v3])

def dihedral_3d(va, vb, vc, vd):
    """Dihedral angle at edge (a,b) with opposite vertices c,d."""
    e = vb - va
    e_hat = e / np.linalg.norm(e)
    ac = vc - va; ad = vd - va
    ac_p = ac - np.dot(ac, e_hat)*e_hat
    ad_p = ad - np.dot(ad, e_hat)*e_hat
    n1, n2 = np.linalg.norm(ac_p), np.linalg.norm(ad_p)
    if n1 < 1e-15 or n2 < 1e-15:
        return 0.0
    return np.arccos(np.clip(np.dot(ac_p, ad_p)/(n1*n2), -1, 1))

def tet_dihedral_angles(L):
    """6 dihedral angles from 6 edge lengths."""
    v = tet_reconstruct_3d(L)
    angles = np.zeros(6)
    for idx, (i,j) in enumerate(TET_EDGE_PAIRS):
        others = [x for x in range(4) if x != i and x != j]
        angles[idx] = dihedral_3d(v[i], v[j], v[others[0]], v[others[1]])
    return angles

# Regular tet: all edges = 1/phi
l_edge = 1.0 / PHI
L0 = np.array([l_edge]*6)
theta0 = tet_dihedral_angles(L0)
print(f"  Regular tet dihedral: {theta0[0]:.6f} rad = {np.degrees(theta0[0]):.4f} deg")
print(f"  arccos(1/3) = {np.arccos(1/3):.6f} rad")
print(f"  All equal: {np.allclose(theta0, theta0[0])}")

# Deficit angle
epsilon0 = 2*np.pi - 5*theta0[0]
print(f"  Deficit angle: {epsilon0:.6f} rad = {np.degrees(epsilon0):.4f} deg")

# --- M matrix: M[i,j] = d(theta_i)/dl_j ---
h_fd = 1e-7
theta_plus = np.zeros((6, 6))  # theta_plus[j] = theta vector when l_j perturbed +h
theta_minus = np.zeros((6, 6))
for j in range(6):
    Lp = L0.copy(); Lp[j] += h_fd
    Lm = L0.copy(); Lm[j] -= h_fd
    theta_plus[j] = tet_dihedral_angles(Lp)
    theta_minus[j] = tet_dihedral_angles(Lm)

M = np.zeros((6, 6))
for i in range(6):
    for j in range(6):
        M[i,j] = (theta_plus[j,i] - theta_minus[j,i]) / (2*h_fd)

print(f"\n  M matrix (d theta_i / dl_j) for regular tet:")
# Check structure: diagonal, adjacent, opposite
diag_vals = [M[i,i] for i in range(6)]
adj_vals = []
opp_vals = []
for i in range(6):
    vi, vj = TET_EDGE_PAIRS[i]
    for j in range(6):
        if j == i: continue
        vk, vl = TET_EDGE_PAIRS[j]
        shared = {vi,vj} & {vk,vl}
        if shared:
            adj_vals.append(M[i,j])
        else:
            opp_vals.append(M[i,j])
print(f"  M_diag = {np.mean(diag_vals):.8f} (all: {np.std(diag_vals):.2e})")
print(f"  M_adj  = {np.mean(adj_vals):.8f} (all: {np.std(adj_vals):.2e})")
print(f"  M_opp  = {np.mean(opp_vals):.8f} (all: {np.std(opp_vals):.2e})")
print(f"  M symmetry: |M - M^T| = {np.max(np.abs(M - M.T)):.2e}")

# Verify Schlaefli: sum_i l_i * M[i,j] = 0 for all j
schlaefli_check = L0 @ M
print(f"  Schlaefli: |sum l*M[:,j]| = {np.max(np.abs(schlaefli_check)):.2e}")

# --- N tensor: N[i,j,k] = d^2(theta_i)/(dl_j dl_k) ---
print(f"\n  Computing N tensor (d^2 theta / dl dl)...")
N = np.zeros((6, 6, 6))
for j in range(6):
    # Diagonal: d^2 theta_i / dl_j^2
    N[:, j, j] = (theta_plus[j] - 2*theta0 + theta_minus[j]) / h_fd**2
    # Off-diagonal: d^2 theta_i / (dl_j dl_k)
    for k in range(j+1, 6):
        Lpp = L0.copy(); Lpp[j] += h_fd; Lpp[k] += h_fd
        Lpm = L0.copy(); Lpm[j] += h_fd; Lpm[k] -= h_fd
        Lmp = L0.copy(); Lmp[j] -= h_fd; Lmp[k] += h_fd
        Lmm = L0.copy(); Lmm[j] -= h_fd; Lmm[k] -= h_fd
        tpp = tet_dihedral_angles(Lpp)
        tpm = tet_dihedral_angles(Lpm)
        tmp_ = tet_dihedral_angles(Lmp)
        tmm = tet_dihedral_angles(Lmm)
        N[:, j, k] = (tpp - tpm - tmp_ + tmm) / (4*h_fd**2)
        N[:, k, j] = N[:, j, k]

# Verify N symmetry in last two indices (should be exact by construction)
for i in range(6):
    print(f"  N[{i}] symmetry: |N[i]-N[i]^T| = {np.max(np.abs(N[i]-N[i].T)):.2e}")

print(f"  N range: [{N.min():.4f}, {N.max():.4f}]")
print(f"  Per-tet computation: {time.time()-t1:.1f}s")


# ============================================================
# PART 3: ASSEMBLE GLOBAL HESSIAN
# ============================================================
print("\n" + "=" * 72)
print("PART 3: ASSEMBLE REGGE HESSIAN H_edge (720x720)")
print("=" * 72)
t2 = time.time()

# H[ea, eb] = -sum_{T containing ea,eb} M[loc_ea, loc_eb]
H_edge = np.zeros((N_e, N_e))
for ci in range(N_c):
    te = tet_edges_list[ci]
    for i_loc in range(6):
        for j_loc in range(6):
            H_edge[te[i_loc], te[j_loc]] -= M[i_loc, j_loc]

H_edge = 0.5 * (H_edge + H_edge.T)  # enforce symmetry
print(f"  H_edge symmetry: |H - H^T| = {np.max(np.abs(H_edge - H_edge.T)):.2e}")

evals_H = np.linalg.eigvalsh(H_edge)
n_pos = np.sum(evals_H > 0.01)
n_neg = np.sum(evals_H < -0.01)
n_zero_H = np.sum(np.abs(evals_H) < 0.01)
print(f"  Spectrum: {n_pos} positive, {n_neg} negative, {n_zero_H} zero")
print(f"  Range: [{evals_H[0]:.4f}, {evals_H[-1]:.4f}]")

# Project to coexact
H_co = P_coexact @ H_edge @ P_coexact
evals_Hco = np.linalg.eigvalsh(H_co)
evals_Hco_nz = evals_Hco[np.abs(evals_Hco) > 0.01]
norm_H_co = np.linalg.norm(H_co)
print(f"  H_co: {len(evals_Hco_nz)} nonzero eigenvalues, ||H_co|| = {norm_H_co:.6f}")
if len(evals_Hco_nz) > 0:
    print(f"  Range: [{evals_Hco_nz.min():.6f}, {evals_Hco_nz.max():.6f}]")

# Compare with L_coexact
ratio_HL = np.linalg.norm(H_co) / np.linalg.norm(L_coexact) if np.linalg.norm(L_coexact) > 0 else 0
print(f"  ||H_co|| / ||L_coexact|| = {ratio_HL:.6f}")
print(f"  Hessian assembly: {time.time()-t2:.1f}s")


# ============================================================
# PART 4: ASSEMBLE CUBIC VERTEX
# ============================================================
print("\n" + "=" * 72)
print("PART 4: CUBIC VERTEX V_eee")
print("=" * 72)
t3 = time.time()

# V[ea,eb,ec] = d^3S/(dl_ea dl_eb dl_ec) = -sum_T N[loc_ea, loc_eb, loc_ec]
# Store as dict of sorted triples
V_cubic = defaultdict(float)
for ci in range(N_c):
    te = tet_edges_list[ci]
    g2l = {te[k]: k for k in range(6)}
    # All unique sorted triples from this tet's 6 edges
    triples_seen = set()
    for a in te:
        for b in te:
            for c in te:
                triple = tuple(sorted([a, b, c]))
                if triple not in triples_seen:
                    triples_seen.add(triple)
                    la, lb, lc = g2l[triple[0]], g2l[triple[1]], g2l[triple[2]]
                    V_cubic[triple] -= N[la, lb, lc]

n_entries = len(V_cubic)
V_vals = np.array(list(V_cubic.values()))
V_nonzero = V_vals[np.abs(V_vals) > 1e-6]
print(f"  Unique sorted triples: {n_entries}")
print(f"  Nonzero (|V| > 1e-6): {len(V_nonzero)}")
print(f"  V range: [{V_vals.min():.6f}, {V_vals.max():.6f}]")
print(f"  |V| RMS: {np.sqrt(np.mean(V_vals**2)):.6e}")

# Verify symmetry: V[ea,eb,ec] should equal V computed with different theta index
# Check a few random triples
print(f"  Symmetry check (V[a,b,c] via different theta indices):")
n_check = 0
max_asym = 0.0
for triple, val in list(V_cubic.items())[:20]:
    ea, eb, ec = triple
    # The stored value uses triple[0] as theta index.
    # Check: compute V[eb,ea,ec] = -sum_T N[loc_eb, loc_ea, loc_ec]
    val2 = 0.0
    for ci in range(N_c):
        te = tet_edges_list[ci]
        te_set = set(te)
        if ea in te_set and eb in te_set and ec in te_set:
            g2l = {te[k]: k for k in range(6)}
            val2 -= N[g2l[eb], g2l[ea], g2l[ec]]
    asym = abs(val - val2)
    if abs(val) > 1e-10:
        max_asym = max(max_asym, asym / abs(val))
    n_check += 1
print(f"  Max relative asymmetry over {n_check} triples: {max_asym:.2e}")
print(f"  Cubic vertex: {time.time()-t3:.1f}s")


# ============================================================
# PART 5: ONE-LOOP SELF-ENERGY
# ============================================================
print("\n" + "=" * 72)
print("PART 5: ONE-LOOP SELF-ENERGY Pi (bubble diagram)")
print("=" * 72)
t4 = time.time()

# Pi[ea,eb] = sum_{e3,e4} V[ea,e3,e4] * G[e3,e3'] * G[e4,e4'] * V[eb,e3',e4']
#           = sum entries_a v_a * sum entries_b v_b * G[e3a,e3b] * G[e4a,e4b]

# Build per-edge V lookup: ea -> [(e3, e4, val), ...]
# V is symmetric in all 3 indices, so for each sorted triple (a,b,c):
# it contributes to V_per_edge[a] as (b,c,val), V_per_edge[b] as (a,c,val), etc.
edge_V = defaultdict(list)
for (e1, e2, e3), val in V_cubic.items():
    if abs(val) < 1e-14:
        continue
    for perm in set(itertools.permutations([e1, e2, e3])):
        ea, eb, ec = perm
        edge_V[ea].append((eb, ec, val))

# Deduplicate per-edge entries
V_per_edge = {}
for ea, entries in edge_V.items():
    seen = set()
    deduped = []
    for (eb, ec, v) in entries:
        key = (eb, ec)
        if key not in seen:
            seen.add(key)
            deduped.append((eb, ec, v))
    V_per_edge[ea] = deduped

n_active = len(V_per_edge)
avg_ent = np.mean([len(v) for v in V_per_edge.values()])
print(f"  Active edges: {n_active}, avg entries per edge: {avg_ent:.1f}")

# Compute Pi via M_ea approach:
# For each ea: M_ea = sum v * outer(G[e3,:], G[e4,:])
# Then Pi[ea,eb] = sum_{(e3',e4',vb) in V_eb} M_ea[e3',e4'] * vb
Pi = np.zeros((N_e, N_e))
active_sorted = sorted(V_per_edge.keys())

print(f"  Computing Pi ({n_active} active edges)...")
t_pi = time.time()
last_report = time.time()

for ea_idx, ea in enumerate(active_sorted):
    if time.time() - last_report > 10:
        elapsed = time.time() - t_pi
        frac = (ea_idx+1)/n_active
        print(f"    {ea_idx+1}/{n_active} ({100*frac:.0f}%), "
              f"elapsed {elapsed:.0f}s, ETA {elapsed/frac - elapsed:.0f}s")
        last_report = time.time()

    entries_a = V_per_edge[ea]
    # Build M_ea = sum v * outer(G[e3], G[e4])
    M_ea = np.zeros((N_e, N_e))
    for (e3, e4, v) in entries_a:
        M_ea += v * np.outer(G_matrix[e3], G_matrix[e4])

    # Contract with V_eb for all eb
    for eb in active_sorted:
        val = 0.0
        for (e3p, e4p, vb) in V_per_edge[eb]:
            val += M_ea[e3p, e4p] * vb
        Pi[ea, eb] += val

Pi = 0.5 * (Pi + Pi.T)
t_pi_done = time.time()

print(f"  Pi computed in {t_pi_done - t_pi:.1f}s")
print(f"  Pi symmetry: |Pi - Pi^T| = {np.max(np.abs(Pi - Pi.T)):.2e}")
print(f"  ||Pi|| (Frobenius): {np.linalg.norm(Pi):.6e}")
print(f"  Tr(Pi): {np.trace(Pi):.6e}")
print(f"  Pi range: [{Pi.min():.6e}, {Pi.max():.6e}]")


# ============================================================
# PART 6: FINITENESS ANALYSIS
# ============================================================
print("\n" + "=" * 72)
print("PART 6: FINITENESS ANALYSIS")
print("=" * 72)

Pi_co = P_coexact @ Pi @ P_coexact
norm_Pi_co = np.linalg.norm(Pi_co)
Pi_co_spectral = V_co.T @ Pi @ V_co  # 601x601
evals_Pi = np.linalg.eigvalsh(Pi_co_spectral)

print(f"  ||Pi_co|| = {norm_Pi_co:.6e}")
print(f"  Tr(Pi_co) = {np.trace(Pi_co):.6e}")
print(f"  Pi_co eigenvalue range: [{evals_Pi.min():.6e}, {evals_Pi.max():.6e}]")
print(f"  Positive: {np.sum(evals_Pi > 1e-15)}, Negative: {np.sum(evals_Pi < -1e-15)}")

alpha_loop = norm_Pi_co / norm_H_co if norm_H_co > 1e-15 else float('inf')

print(f"\n  FINITENESS: Every entry of Pi is a FINITE sum of finite products.")
print(f"  In continuum QG: same diagram diverges as Lambda^2. Here: FINITE.")
print(f"\n  Loop expansion parameter:")
print(f"    alpha_loop = ||Pi_co|| / ||H_co|| = {alpha_loop:.6e}")
print(f"  Comparison:")
print(f"    alpha_em   = {ALPHA:.6e}")
print(f"    alpha_s    = {1/(2*PHI**3):.6e}")
print(f"    1/N        = {1/120:.6e}")
print(f"    1/degree   = {1/12:.6e}")
print(f"    phi^{{-2}}   = {PHI**-2:.6e}")

Tr_Pi_co = np.trace(Pi_co)
Tr_G_Pi = np.trace(G_matrix @ Pi_co)
print(f"\n  Vacuum energy contributions:")
print(f"    Tr(Pi_co) = {Tr_Pi_co:.6e}")
print(f"    Tr(G @ Pi_co) = {Tr_G_Pi:.6e}")


# ============================================================
# PART 7: GRAVITON MASS PROTECTION
# ============================================================
print("\n" + "=" * 72)
print("PART 7: GRAVITON MASS PROTECTION")
print("=" * 72)

L_co_spectral = V_co.T @ L_coexact @ V_co
L_eff_spectral = L_co_spectral + Pi_co_spectral
evals_eff = np.sort(np.linalg.eigvalsh(L_eff_spectral))

tree_gap = lambda_co.min()
loop_gap = evals_eff[0]
delta_gap = loop_gap - tree_gap

print(f"  Tree-level mass gap: {tree_gap:.6f}")
print(f"  One-loop mass gap:  {loop_gap:.6f}")
print(f"  Shift: {delta_gap:.6e}")
print(f"  Relative shift: {abs(delta_gap)/tree_gap:.6e}")

if abs(delta_gap)/tree_gap < 0.01:
    gap_status = "PROTECTED (< 1%)"
elif abs(delta_gap)/tree_gap < 0.1:
    gap_status = "Moderately shifted"
else:
    gap_status = "Significantly shifted"
print(f"  --> Mass gap: {gap_status}")

print(f"\n  Lowest eigenvalues (tree vs one-loop):")
print(f"  {'tree':>12s} {'one-loop':>12s} {'shift':>12s} {'rel':>12s}")
lco_sorted = np.sort(lambda_co)
for k in range(min(15, len(lco_sorted))):
    tv, lv = lco_sorted[k], evals_eff[k]
    sh = lv - tv
    print(f"  {tv:12.6f} {lv:12.6f} {sh:12.6e} {abs(sh)/tv:12.6e}")

# Ward identity
Pi_exact_proj = P_exact @ Pi @ P_exact
norm_Pi_exact = np.linalg.norm(Pi_exact_proj)
norm_Pi = np.linalg.norm(Pi)
ward_ratio = norm_Pi_exact / norm_Pi if norm_Pi > 1e-15 else 0
Pi_cross = P_exact @ Pi @ P_coexact
norm_cross = np.linalg.norm(Pi_cross)

print(f"\n  Ward identity:")
print(f"    ||P_exact @ Pi @ P_exact|| / ||Pi|| = {ward_ratio:.6e}")
print(f"    ||P_exact @ Pi @ P_coexact|| = {norm_cross:.6e}")
ward_status = "SATISFIED" if ward_ratio < 1e-3 else f"BROKEN ({ward_ratio:.2e})"
print(f"    --> {ward_status}")


# ============================================================
# PART 8: FRAMEWORK CONNECTIONS
# ============================================================
print("\n" + "=" * 72)
print("PART 8: FRAMEWORK CONNECTIONS")
print("=" * 72)

print(f"\n  alpha_loop = {alpha_loop:.6e}")
test_vals = {
    'alpha': ALPHA, 'alpha^2': ALPHA**2,
    '1/(2*phi^3)': 1/(2*PHI**3), 'alpha*phi': ALPHA*PHI,
    '1/N': 1/120, '1/N^2': 1/14400,
    '1/degree': 1/12, 'phi^(-4)': PHI**(-4),
    'phi^(-5)': PHI**(-5), '1/(8*pi^2)': 1/(8*np.pi**2),
    '1/(16*pi^2)': 1/(16*np.pi**2),
}
print(f"  {'expression':>16s} {'value':>14s} {'ratio':>14s}")
for name, val in sorted(test_vals.items(),
        key=lambda x: abs(np.log(max(abs(x[1]),1e-30)/max(abs(alpha_loop),1e-30)))):
    if abs(alpha_loop) > 1e-30 and abs(val) > 1e-30:
        print(f"  {name:>16s} {val:14.6e} {alpha_loop/val:14.6f}")

G_0 = 1.0/12
Tr_G_L = np.trace(G_matrix @ L_coexact)
if abs(Tr_G_L) > 1e-10:
    delta_G = Tr_G_Pi / Tr_G_L
    print(f"\n  Newton's constant correction:")
    print(f"    G_0 = 1/12 = {G_0:.6f}")
    print(f"    delta_G = Tr(G*Pi)/Tr(G*L) = {delta_G:.6e}")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
total_time = time.time() - t0

print(f"""
600-CELL ONE-LOOP GRAVITON SELF-ENERGY
======================================

INFRASTRUCTURE:
  600-cell: V=120, E=720, F=1200, C=600
  Coexact DOF: {n_coexact}, Mass gap: {tree_gap:.6f} = 7-4*phi

REGGE HESSIAN (tet-local, Schlaefli-based):
  M_diag = {np.mean(diag_vals):.8f}, M_adj = {np.mean(adj_vals):.8f}, M_opp = {np.mean(opp_vals):.8f}
  H_edge: {n_pos} pos, {n_neg} neg, {n_zero_H} zero eigenvalues
  ||H_co|| = {norm_H_co:.6e}

CUBIC VERTEX:
  {n_entries} unique triples, {len(V_nonzero)} nonzero, |V| RMS = {np.sqrt(np.mean(V_vals**2)):.6e}

ONE-LOOP SELF-ENERGY:
  ||Pi|| = {norm_Pi:.6e}, ||Pi_co|| = {norm_Pi_co:.6e}
  alpha_loop = ||Pi_co||/||H_co|| = {alpha_loop:.6e}

FINITENESS: TRIVIALLY TRUE (finite sums on finite graph)
  Continuum: Lambda^2 divergence. Here: graph = UV regulator.

MASS PROTECTION: {gap_status}
  Tree gap: {tree_gap:.6f}, One-loop gap: {loop_gap:.6f}
  Shift: {delta_gap:.6e} (relative: {abs(delta_gap)/tree_gap:.6e})

WARD IDENTITY: {ward_status}

TOTAL TIME: {total_time:.1f}s
""")
