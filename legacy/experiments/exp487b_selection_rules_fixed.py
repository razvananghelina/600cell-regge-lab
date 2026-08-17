"""
EXP-487b: Selection Rules FIXED - Using d0*psi as gauge modes
==============================================================
FIX: exp487 used arbitrary eigenvectors of C for exact modes.
The CORRECT gauge modes are d0*psi_alpha, the images of vertex
eigenmodes under the boundary operator. These inherit the 2I
representation structure from the vertex modes.

KEY FORMULA:
  V(alpha, beta, gamma) = sum_e [psi_a(j)-psi_a(i)] * psi_b(i) * psi_g(j)

  where alpha = gauge mode (d0*psi_alpha), beta,gamma = fermion modes,
  and (i,j) = edge endpoints.

All three indices are vertex eigenmodes -> the selection rules are
determined by the Clebsch-Gordan coefficients of 2I.
"""

import numpy as np
from collections import defaultdict
from scipy.spatial.distance import cdist

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6

print("=" * 70)
print("EXP-487b: SELECTION RULES (FIXED GAUGE MODE BASIS)")
print("=" * 70)

# === Build 600-cell (compact) ===
vertices = []
for i in range(4):
    for s in [1,-1]:
        v = [0.0]*4; v[i] = s; vertices.append(v)
for s0 in [1,-1]:
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            for s3 in [1,-1]:
                vertices.append([s0/2, s1/2, s2/2, s3/2])
p, h, q = phi/2, 0.5, 1/(2*phi)
for t in [[p,h,q,0],[p,q,0,h],[p,0,h,q],[h,p,0,q],[h,q,p,0],[h,0,q,p],
          [q,p,h,0],[q,h,0,p],[q,0,p,h],[0,p,q,h],[0,h,p,q],[0,q,h,p]]:
    for s0 in [1,-1]:
        for s1 in [1,-1]:
            for s2 in [1,-1]:
                for s3 in [1,-1]:
                    v = [s0*t[0], s1*t[1], s2*t[2], s3*t[3]]
                    if np.linalg.norm(v) > 0.01: vertices.append(v)
unique = []
for v in vertices:
    va = np.array(v)
    if not any(np.allclose(va, u, atol=1e-10) for u in unique):
        unique.append(va)
vertices = np.array(unique)
N = len(vertices)
norms = np.linalg.norm(vertices, axis=1)
vertices = vertices / norms[:, np.newaxis]
assert N == 120

dist_matrix = cdist(vertices, vertices, 'euclidean')
np.fill_diagonal(dist_matrix, np.inf)
nn_dist = np.median(dist_matrix.min(axis=1))
adj = (np.abs(dist_matrix - nn_dist) < 0.02).astype(int)
np.fill_diagonal(adj, 0)

edges = []
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j]:
            edges.append((i,j))
E = len(edges)

from scipy.sparse import lil_matrix, csr_matrix
d0 = lil_matrix((E, N))
for idx, (i,j) in enumerate(edges):
    d0[idx, i] = -1; d0[idx, j] = 1
d0_dense = d0.toarray()

print(f"  N={N}, E={E}")

# === Eigendecomposition of Delta_0 ===
Delta_0 = d0_dense.T @ d0_dense  # (120 x 120)
evals_0, evecs_0 = np.linalg.eigh(Delta_0)

# Adjacency eigenvalues
adj_matrix = adj.astype(float)
evals_adj = 12 - evals_0  # Delta_0 = 12*I - A

# Group vertex modes by representation
rep_groups = []
current_eval = None
current_group = []
for i in range(N):
    if current_eval is None or abs(evals_0[i] - current_eval) > 0.01:
        if current_group:
            dim = int(round(np.sqrt(len(current_group))))
            rep_groups.append((current_eval, current_group, dim))
        current_eval = evals_0[i]
        current_group = [i]
    else:
        current_group.append(i)
if current_group:
    dim = int(round(np.sqrt(len(current_group))))
    rep_groups.append((current_eval, current_group, dim))

n_reps = len(rep_groups)
print(f"\n  9 representations of 2I:")
print(f"  {'Rep':>4} {'dim':>4} {'lam_adj':>10} {'lam_Delta':>10} {'mult':>6}")
for r, (lam, modes, dim) in enumerate(rep_groups):
    print(f"  {r:4d} {dim:4d} {12-lam:10.4f} {lam:10.4f} {len(modes):6d}")

# ============================================================
# KEY COMPUTATION: Gauge modes = d0 * vertex eigenmodes
# ============================================================
print("\n" + "=" * 70)
print("GAUGE MODES = d0 * vertex eigenmodes")
print("=" * 70)

# G[e, alpha] = (d0 @ psi_alpha)(e) = psi_alpha(j) - psi_alpha(i)
# Shape: (720, 120)
G = d0_dense @ evecs_0  # (720, 120)

# Check: G[:, 0] should be zero (constant mode)
print(f"  ||d0 * psi_0|| = {np.linalg.norm(G[:, 0]):.2e} (should be ~0, trivial mode)")

# Gauge modes per rep:
print(f"  Gauge modes per representation:")
for r, (lam, modes, dim) in enumerate(rep_groups):
    norms = [np.linalg.norm(G[:, m]) for m in modes]
    print(f"    Rep {r} (dim={dim}): {len(modes)} modes, "
          f"||d0*psi||: {np.mean(norms):.4f} (should be sqrt(lam)={np.sqrt(lam):.4f})")

# ============================================================
# PAIR CREATION VERTEX: V(alpha, beta, gamma)
# ============================================================
print("\n" + "=" * 70)
print("PAIR CREATION VERTEX V(gauge, fermion1, fermion2)")
print("=" * 70)

edge_v1 = np.array([edges[e][0] for e in range(E)])
edge_v2 = np.array([edges[e][1] for e in range(E)])

# W[e, beta, gamma] = psi_beta(v1(e)) * psi_gamma(v2(e))
psi_v1 = evecs_0[edge_v1, :]  # (720, 120)
psi_v2 = evecs_0[edge_v2, :]  # (720, 120)

# W_flat[e, beta*120 + gamma] = psi_v1[e, beta] * psi_v2[e, gamma]
W_flat = (psi_v1[:, :, np.newaxis] * psi_v2[:, np.newaxis, :]).reshape(E, N*N)

# V[alpha, beta, gamma] = G[:, alpha]^T @ W_flat -> shape (120, 120*120)
print("Computing full vertex tensor (120 x 120 x 120)...")
V_full = G.T @ W_flat  # (120, 14400) -> reshape to (120, 120, 120)
V_full = V_full.reshape(N, N, N)

# The gauge mode alpha=0 (trivial rep) should give V=0
print(f"  V[trivial, :, :] max = {np.max(np.abs(V_full[0])):.2e} (should be ~0)")

# ============================================================
# REPRESENTATION-RESOLVED VERTEX
# ============================================================
print("\n" + "=" * 70)
print("REPRESENTATION-RESOLVED |V|^2")
print("=" * 70)

V_rep = np.zeros((n_reps, n_reps, n_reps))

for r_g, (lam_g, modes_g, dim_g) in enumerate(rep_groups):
    for r_f1, (lam_f1, modes_f1, dim_f1) in enumerate(rep_groups):
        for r_f2, (lam_f2, modes_f2, dim_f2) in enumerate(rep_groups):
            total = 0.0
            for a in modes_g:
                for b in modes_f1:
                    for c in modes_f2:
                        total += V_full[a, b, c]**2
            V_rep[r_g, r_f1, r_f2] = total

# Normalize by number of modes for comparison
V_rep_norm = np.zeros((n_reps, n_reps, n_reps))
for r_g in range(n_reps):
    for r_f1 in range(n_reps):
        for r_f2 in range(n_reps):
            n_modes = len(rep_groups[r_g][1]) * len(rep_groups[r_f1][1]) * len(rep_groups[r_f2][1])
            if n_modes > 0:
                V_rep_norm[r_g, r_f1, r_f2] = V_rep[r_g, r_f1, r_f2] / n_modes

# ============================================================
# SELECTION RULE MATRIX
# ============================================================
print("\n  Selection rule matrix (sum |V|^2 per rep triple):")
print(f"  Nonzero entries (|V|^2 > 1e-10):\n")

threshold = 1e-10

for r_g in range(n_reps):
    dim_g = rep_groups[r_g][2]
    if r_g == 0:
        continue  # trivial rep gives zero gauge mode

    total_nonzero = sum(1 for r1 in range(n_reps)
                       for r2 in range(n_reps)
                       if V_rep[r_g, r1, r2] > threshold)
    total_pairs = n_reps * n_reps

    print(f"  Gauge rep {r_g} (dim={dim_g}, lam_adj={12-rep_groups[r_g][0]:.2f}):")
    print(f"  {'':>8}", end="")
    for r2 in range(n_reps):
        print(f" r{r2}({rep_groups[r2][2]})", end="")
    print()

    for r1 in range(n_reps):
        print(f"  r{r1}({rep_groups[r1][2]:>1})", end=" ")
        for r2 in range(n_reps):
            val = V_rep[r_g, r1, r2]
            if val > threshold:
                # Use log scale for readability
                mag = np.log10(val)
                if mag > -1:
                    print(f" {val:6.3f}", end="")
                else:
                    print(f"  {val:5.0e}", end="")
            else:
                print(f"      .", end="")
        print()
    print(f"    [{total_nonzero}/{total_pairs} nonzero]\n")

# ============================================================
# McKAY GRAPH EXTRACTION
# ============================================================
print("=" * 70)
print("McKAY GRAPH FROM PAIR CREATION")
print("=" * 70)

# The McKay graph M[k,j] encodes: rho_1 x rho_k contains rho_j
# In our vertex: gauge mode from rep k creates fermion pair.
# The gauge mode is d0*psi_k, which is in rep rho_1 x rho_k.
# So V(gauge_k, ferm_j, ferm_?) != 0 means rho_j x rho_? contains rho_1 x rho_k.

# For the McKay graph with the fundamental rep (r=1, dim=2):
# Check which fermion pairs (j1, j2) are created by gauge rep 1:
print(f"\nGauge rep 1 (fundamental, dim=2, lam_adj={12-rep_groups[1][0]:.4f}):")
print("  Fermion pairs with nonzero V:")
for r1 in range(n_reps):
    for r2 in range(r1, n_reps):
        val = V_rep[1, r1, r2]
        if val > threshold:
            print(f"    ({r1},{r2}) = (dim {rep_groups[r1][2]}, dim {rep_groups[r2][2]})  |V|^2 = {val:.6f}")

# Actually, for the McKay graph, I need to check:
# When the gauge mode comes from rep k (i.e., d0 * psi_k),
# which single-fermion reps does it connect to?
# V(k, j1, j2) != 0 means rho_1 x rho_k has a component in rho_j1 x rho_j2.
# The McKay adjacency is about which reps appear in rho_1 x rho_k.

# DIRECT APPROACH: The McKay adjacency matrix for rho_1 can be read from
# the adjacency matrix spectrum. M[k,j] = 1 iff rho_j appears in rho_1 x rho_k.
# This is equivalent to: A restricted to rho_k eigenspace has overlap with rho_j eigenspace.

# Compute: P_j * A * P_k, where P_j projects onto eigenspace j
print("\n--- McKay graph from adjacency matrix ---")

projectors = {}
for r, (lam, modes, dim) in enumerate(rep_groups):
    # Projector onto eigenspace r
    P = np.zeros((N, N))
    for m in modes:
        P += np.outer(evecs_0[:, m], evecs_0[:, m])
    projectors[r] = P

# McKay matrix: M[k,j] = Tr(P_j A P_k) / (dim_k * eigenvalue)
# Or simpler: M[k,j] = ||P_j A P_k||_F^2 / (dim_k^2 * dim_j^2) ... normalized somehow

# Actually, the cleanest way: for each rep k, compute A * P_k and project onto each rep j
M_mckay = np.zeros((n_reps, n_reps))
for k in range(n_reps):
    for j in range(n_reps):
        AP_k = adj_matrix @ projectors[k]
        val = np.trace(projectors[j] @ AP_k)
        # Normalize: M[k,j] should give the multiplicity of rho_j in rho_1 x rho_k
        # Tr(P_j A P_k) = sum_{m in j, n in k} <psi_m, A psi_n>
        # = sum_{m in j, n in k} lambda_k * delta_{mn} if j == k
        # In general: Tr(P_j A P_k) = lambda_k * Tr(P_j P_k) = lambda_k * dim_k * delta_{jk}
        # No, that's only for the DIAGONAL. The OFF-DIAGONAL gives the McKay structure.
        M_mckay[k, j] = val

# The diagonal should give lambda_k * dim_k^2
print(f"  McKay matrix (raw Tr(P_j A P_k)):")
print(f"  {'':>5}", end="")
for j in range(n_reps):
    print(f"  r{j}({rep_groups[j][2]})", end="")
print()
for k in range(n_reps):
    print(f"  r{k}({rep_groups[k][2]})", end="")
    for j in range(n_reps):
        val = M_mckay[k, j]
        if abs(val) > 0.01:
            print(f"  {val:6.1f}", end="")
        else:
            print(f"      .", end="")
    print()

# Normalize to get multiplicity: M_norm[k,j] = M[k,j] / (dim_k * dim_j)
# Wait, better: for the McKay graph, we expect M[k,j] to be 0 or 1 (for simply-laced)
# Let's normalize: the adjacency eigenvalue lambda_k = sum_j M[k,j] * dim_j
# And Tr(P_j A P_k) = M[k,j] * dim_j * dim_k (if I understand correctly)

print(f"\n  Normalized McKay adjacency M[k,j] = Tr(P_j A P_k) / (dim_k * dim_j):")
print(f"  {'':>5}", end="")
for j in range(n_reps):
    print(f" r{j}({rep_groups[j][2]})", end="")
print()

M_norm = np.zeros((n_reps, n_reps))
for k in range(n_reps):
    dim_k = rep_groups[k][2]
    for j in range(n_reps):
        dim_j = rep_groups[j][2]
        M_norm[k, j] = M_mckay[k, j] / (dim_k**2 * dim_j)
        # This should give the EIGENVALUE of the McKay adjacency matrix?

for k in range(n_reps):
    print(f"  r{k}({rep_groups[k][2]})", end="")
    for j in range(n_reps):
        val = M_norm[k, j]
        if abs(val) > 0.01:
            print(f" {val:6.2f}", end="")
        else:
            print(f"      .", end="")
    print()

# SIMPLEST METHOD: The McKay graph adjacency is just whether
# |Tr(P_j A P_k)| > 0 for j != k.
print(f"\n  McKay graph edges (|Tr(P_j A P_k)| > 0.01, j != k):")
mckay_edges = []
for k in range(n_reps):
    for j in range(k+1, n_reps):
        if abs(M_mckay[k, j]) > 0.01:
            dim_k = rep_groups[k][2]
            dim_j = rep_groups[j][2]
            print(f"    r{k}(dim={dim_k}) -- r{j}(dim={dim_j})  "
                  f"Tr={M_mckay[k,j]:.2f}")
            mckay_edges.append((k, j))

print(f"\n  Total edges: {len(mckay_edges)} (expected 8 for extended E8)")
degrees_mckay = [0] * n_reps
for k, j in mckay_edges:
    degrees_mckay[k] += 1
    degrees_mckay[j] += 1
print(f"  Degrees: {degrees_mckay}")
print(f"  Branching node(s): {[r for r in range(n_reps) if degrees_mckay[r] == 3]}")

# ============================================================
# SELECTION RULES FROM V: WHICH PAIRS ARE FORBIDDEN?
# ============================================================
print("\n" + "=" * 70)
print("PAIR CREATION SELECTION RULES")
print("=" * 70)

# For each gauge rep, list ALLOWED and FORBIDDEN fermion pairs
for r_g in range(1, n_reps):  # skip trivial
    dim_g = rep_groups[r_g][2]
    lam_g = 12 - rep_groups[r_g][0]

    allowed = []
    forbidden = []
    for r1 in range(n_reps):
        for r2 in range(r1, n_reps):
            val = V_rep[r_g, r1, r2]
            if val > threshold:
                allowed.append((r1, r2, val))
            else:
                forbidden.append((r1, r2))

    print(f"\n  Gauge rep {r_g} (dim={dim_g}, lam_adj={lam_g:.2f}):")
    print(f"    ALLOWED ({len(allowed)} channels):")
    for r1, r2, val in sorted(allowed, key=lambda x: -x[2])[:10]:
        print(f"      r{r1}(d{rep_groups[r1][2]}) x r{r2}(d{rep_groups[r2][2]})  |V|^2={val:.4f}")
    if len(allowed) > 10:
        print(f"      ... and {len(allowed)-10} more")
    print(f"    FORBIDDEN ({len(forbidden)} channels):")
    for r1, r2 in forbidden[:5]:
        print(f"      r{r1}(d{rep_groups[r1][2]}) x r{r2}(d{rep_groups[r2][2]})")
    if len(forbidden) > 5:
        print(f"      ... and {len(forbidden)-5} more")

# ============================================================
# CROSS-CHECK: Does V respect the McKay graph?
# ============================================================
print("\n" + "=" * 70)
print("CROSS-CHECK: V vs McKay GRAPH")
print("=" * 70)

# The gauge mode from rep k is d0*psi_k, which lives in rho_1 x rho_k.
# The McKay graph tells us which reps appear in rho_1 x rho_k.
# The pair creation V(k, j1, j2) should be nonzero iff
# some rep in rho_1 x rho_k also appears in rho_j1 x rho_j2.

# First: compute which reps are in rho_1 x rho_k (from McKay graph)
print("\n  rho_1 x rho_k decomposition (from McKay edges):")
for k in range(n_reps):
    neighbors = [j for j in range(n_reps)
                 if (min(k,j), max(k,j)) in [(a,b) for a,b in mckay_edges]
                 or k == j and abs(M_mckay[k,k]) > 0.01]
    # Actually, recompute from M_mckay directly
    components = [j for j in range(n_reps) if abs(M_mckay[k,j]) > 0.01]
    dims = [rep_groups[j][2] for j in components]
    total = sum(d for d in dims)  # should equal dim(rho_1)*dim(rho_k) = 2*dim_k
    dim_k = rep_groups[k][2]
    expected = 12  # This is lambda_k... wait no
    # rho_1 x rho_k has dimension 2*dim_k
    print(f"    rho_1 x rho_{k}(d{dim_k}) = "
          + " + ".join(f"rho_{j}(d{rep_groups[j][2]})" for j in components)
          + f"  [dim check: {total} vs 2*{dim_k}={2*dim_k}]")

# ============================================================
# ANTISYMMETRIC VERTEX (FERMION PAIRS)
# ============================================================
print("\n" + "=" * 70)
print("ANTISYMMETRIC VERTEX (FERMION PAIR CREATION)")
print("=" * 70)

# V_anti(alpha, beta, gamma) = V(alpha, beta, gamma) - V(alpha, gamma, beta)
V_anti = V_full - V_full.transpose(0, 2, 1)

# Check: V_anti for same rep
print("\n  Same-rep fermion pairs (V_anti should be ~ 0 if within same rep and mode):")
for r in range(n_reps):
    modes = rep_groups[r][1]
    total = 0
    count = 0
    for b in modes:
        for c in modes:
            if b < c:
                for a in range(1, N):  # all non-trivial gauge modes
                    total += V_anti[a, b, c]**2
                    count += 1
    if count > 0:
        print(f"    Rep {r} (dim={rep_groups[r][2]}): "
              f"mean |V_anti|^2 = {total/count:.6f} ({count} pairs)")

# Antisymmetric vertex by representation
V_anti_rep = np.zeros((n_reps, n_reps, n_reps))
for r_g in range(n_reps):
    for r_f1 in range(n_reps):
        for r_f2 in range(n_reps):
            total = 0.0
            for a in rep_groups[r_g][1]:
                for b in rep_groups[r_f1][1]:
                    for c in rep_groups[r_f2][1]:
                        total += V_anti[a, b, c]**2
            V_anti_rep[r_g, r_f1, r_f2] = total

print(f"\n  Top antisymmetric channels (gauge -> fermion pair):")
channels = []
for r_g in range(1, n_reps):
    for r1 in range(n_reps):
        for r2 in range(r1+1, n_reps):
            val = V_anti_rep[r_g, r1, r2]
            if val > threshold:
                channels.append((val, r_g, r1, r2))

channels.sort(reverse=True)
print(f"  {'|V_anti|^2':>12} {'gauge':>8} {'ferm1':>8} {'ferm2':>8}")
for val, rg, r1, r2 in channels[:20]:
    print(f"  {val:12.4f}  r{rg}(d{rep_groups[rg][2]})  "
          f"r{r1}(d{rep_groups[r1][2]})  r{r2}(d{rep_groups[r2][2]})")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

n_allowed_full = sum(1 for rg in range(1, n_reps)
                     for r1 in range(n_reps)
                     for r2 in range(r1, n_reps)
                     if V_rep[rg, r1, r2] > threshold)
n_total_full = sum(1 for rg in range(1, n_reps)
                   for r1 in range(n_reps)
                   for r2 in range(r1, n_reps))

n_forbidden = n_total_full - n_allowed_full

print(f"\n  Total (gauge, ferm1, ferm2) triples: {n_total_full}")
print(f"  Allowed: {n_allowed_full}")
print(f"  Forbidden: {n_forbidden}")
print(f"  Selection ratio: {n_allowed_full/n_total_full:.3f}")

print(f"\n  McKay graph: {len(mckay_edges)} edges")
print(f"  Is tree (|E|=|V|-1): {len(mckay_edges) == n_reps - 1}")
print(f"  Degree sequence: {sorted(degrees_mckay, reverse=True)}")

if n_forbidden > 0:
    print(f"\n  RESULT: {n_forbidden} channels are FORBIDDEN by the geometry.")
    print(f"  The pair creation vertex has NON-TRIVIAL selection rules")
    print(f"  derived from the boundary operator d0 and 2I representation theory.")
else:
    print(f"\n  RESULT: ALL channels are allowed.")
    print(f"  The vertex V(alpha, beta, gamma) is nonzero for all rep triples.")
    print(f"  Selection rules must come from ANTISYMMETRY or ENERGY CONSERVATION,")
    print(f"  not from representation theory alone.")

if n_forbidden == 0:
    # Check antisymmetric separately
    n_anti_forbidden = sum(1 for rg in range(1, n_reps)
                          for r1 in range(n_reps)
                          for r2 in range(r1+1, n_reps)
                          if V_anti_rep[rg, r1, r2] < threshold)
    n_anti_total = sum(1 for rg in range(1, n_reps)
                       for r1 in range(n_reps)
                       for r2 in range(r1+1, n_reps))

    print(f"\n  Antisymmetric vertex (fermion pairs):")
    print(f"    Total: {n_anti_total}, Allowed: {n_anti_total - n_anti_forbidden}, "
          f"Forbidden: {n_anti_forbidden}")

    if n_anti_forbidden > 0:
        print(f"\n    FORBIDDEN antisymmetric channels:")
        for rg in range(1, n_reps):
            for r1 in range(n_reps):
                for r2 in range(r1+1, n_reps):
                    if V_anti_rep[rg, r1, r2] < threshold:
                        print(f"      gauge r{rg}(d{rep_groups[rg][2]}) -> "
                              f"r{r1}(d{rep_groups[r1][2]}) ^ r{r2}(d{rep_groups[r2][2]})")

print(f"""
STATUS:
  DERIVED:  Pair creation vertex from d0 (boundary operator)
  COMPUTED: Full 120x120x120 vertex tensor
  COMPUTED: Representation-resolved vertex (9x9x9)
  COMPUTED: McKay graph from adjacency projectors

  The selection rules, antisymmetry structure, and representation
  decomposition are ALL consequences of the 600-cell geometry.
""")

print("=" * 70)
print("EXP-487b COMPLETE")
print("=" * 70)
