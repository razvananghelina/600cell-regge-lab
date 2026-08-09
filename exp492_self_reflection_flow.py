"""
EXP-492: Self-Reflection Spectral Flow on the 600-Cell
========================================================
INTUITION: A point "looking inward" through successive self-reflections.
Each level of self-reflection creates a new particle.
phi = eigenvalue of self-reflection. phi^n = n levels deep.

CONCRETE TESTS:
1. Ball expansion: start at identity e in 2I, expand radius n.
   At each n, which irreps of 2I become "visible" (have support in the ball)?
   Does the ORDER of appearance match the mass hierarchy?

2. Shell spectral content: what is the spectral decomposition of each
   concentric shell around the identity?

3. Self-reflection as McKay walk: M^n * e_0 projected onto phi-eigenvector
   gives phi^n. This IS the mass formula.

4. Cumulative spectral weight: as we include more shells, does the
   "spectral weight" at each irrep grow like phi^{something}?

STATUS: EXPLORATORY (testing the self-reflection hypothesis)
"""

import numpy as np
from scipy.spatial.distance import cdist
from collections import defaultdict, deque

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6

print("=" * 70)
print("EXP-492: SELF-REFLECTION SPECTRAL FLOW")
print("=" * 70)

# ============================================================
# Build 600-cell
# ============================================================
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
verts = np.array(unique)
N = len(verts)
norms = np.linalg.norm(verts, axis=1)
verts = verts / norms[:, np.newaxis]
assert N == 120

dist_matrix = cdist(verts, verts, 'euclidean')
np.fill_diagonal(dist_matrix, np.inf)
nn_dist = np.median(dist_matrix.min(axis=1))
adj = (np.abs(dist_matrix - nn_dist) < 0.02).astype(float)
np.fill_diagonal(adj, 0)

# Identity vertex
v0 = 0  # (1,0,0,0)

# Full adjacency spectrum
evals_adj, evecs_adj = np.linalg.eigh(adj)

# Group into 9 irreps by eigenvalue
rep_groups = []
i = 0
sorted_idx = np.argsort(evals_adj)[::-1]
evals_sorted = evals_adj[sorted_idx]
evecs_sorted = evecs_adj[:, sorted_idx]

i = 0
while i < N:
    lam = evals_sorted[i]
    group = [i]
    while i+1 < N and abs(evals_sorted[i+1] - lam) < 0.01:
        i += 1
        group.append(i)
    dim = int(round(np.sqrt(len(group))))
    rep_groups.append((lam, group, dim))
    i += 1

rep_names = [f"r{i}(d{rep_groups[i][2]},l={rep_groups[i][0]:.1f})" for i in range(9)]
print(f"\n9 irreps: {', '.join(rep_names)}")

# Build projection operators for each irrep
projectors = []
for r, (lam, modes, dim) in enumerate(rep_groups):
    P = np.zeros((N, N))
    for m in modes:
        P += np.outer(evecs_sorted[:, m], evecs_sorted[:, m])
    projectors.append(P)

# ============================================================
# SECTION 1: Distance shells from identity
# ============================================================
print("\n" + "=" * 70)
print("SECTION 1: DISTANCE SHELLS FROM IDENTITY")
print("=" * 70)

# BFS to find distances from v0
graph_dist = {v0: 0}
queue = deque([v0])
while queue:
    node = queue.popleft()
    for j in range(N):
        if adj[node, j] > 0.5 and j not in graph_dist:
            graph_dist[j] = graph_dist[node] + 1
            queue.append(j)

max_dist = max(graph_dist.values())
shells = defaultdict(list)
for v, d in graph_dist.items():
    shells[d].append(v)

print(f"\nShell structure from identity vertex {v0}:")
print(f"  Diameter of 600-cell graph: {max_dist}")
print(f"  {'Shell':>6} {'Size':>6} {'Cumul':>6}")
cumul = 0
for d in range(max_dist + 1):
    size = len(shells[d])
    cumul += size
    print(f"  {d:6d} {size:6d} {cumul:6d}")

# ============================================================
# SECTION 2: Spectral content of each shell
# ============================================================
print("\n" + "=" * 70)
print("SECTION 2: SPECTRAL CONTENT OF EACH SHELL")
print("=" * 70)

print("""
For each shell (vertices at distance d from identity),
compute the "spectral weight" in each irrep:
  w(d, r) = sum_{v in shell(d)} P_r[v, v]
where P_r is the projector onto irrep r.

This measures how much of each irrep is "present" in shell d.
""")

print(f"{'Shell':>6}", end="")
for r in range(9):
    print(f"  r{r}({rep_groups[r][2]})", end="")
print(f"  {'total':>6}")

shell_weights = np.zeros((max_dist + 1, 9))

for d in range(max_dist + 1):
    verts_in_shell = shells[d]
    print(f"  d={d:3d}", end="")
    total = 0
    for r in range(9):
        weight = sum(projectors[r][v, v] for v in verts_in_shell)
        shell_weights[d, r] = weight
        total += weight
        print(f"  {weight:6.2f}", end="")
    print(f"  {total:6.1f}")

# ============================================================
# SECTION 3: Cumulative spectral weight (ball expansion)
# ============================================================
print("\n" + "=" * 70)
print("SECTION 3: CUMULATIVE SPECTRAL WEIGHT (BALL EXPANSION)")
print("=" * 70)

print(f"""
As the "ball of self-reflection" expands from radius 0 to max,
which irreps are "discovered" first?
""")

print(f"{'Radius':>7}", end="")
for r in range(9):
    print(f"  r{r}({rep_groups[r][2]})", end="")
print(f"  {'total':>6}")

cumul_weights = np.zeros((max_dist + 1, 9))
first_appearance = {}  # irrep -> first radius where weight > threshold

threshold_appearance = 0.01

for d in range(max_dist + 1):
    if d == 0:
        cumul_weights[d] = shell_weights[d]
    else:
        cumul_weights[d] = cumul_weights[d-1] + shell_weights[d]

    print(f"  R<={d:3d}", end="")
    total = 0
    for r in range(9):
        w = cumul_weights[d, r]
        total += w
        if r not in first_appearance and w > threshold_appearance:
            first_appearance[r] = d
        print(f"  {w:6.2f}", end="")
    print(f"  {total:6.1f}")

print(f"\nFirst appearance of each irrep (cumul weight > {threshold_appearance}):")
for r in range(9):
    d = first_appearance.get(r, -1)
    dim = rep_groups[r][2]
    lam = rep_groups[r][0]
    print(f"  r{r}(dim={dim}, lam={lam:+.2f}): first at radius {d}")

# Order of appearance
appearance_order = sorted(first_appearance.items(), key=lambda x: (x[1], x[0]))
print(f"\nOrder of appearance: {' -> '.join(f'r{r}(d{rep_groups[r][2]})@R={d}' for r, d in appearance_order)}")

# ============================================================
# SECTION 4: SELF-REFLECTION OPERATOR
# ============================================================
print("\n" + "=" * 70)
print("SECTION 4: SELF-REFLECTION = ADJACENCY AT IDENTITY")
print("=" * 70)

print("""
The "self-reflection" at vertex v0 is: looking at the neighbors.
The adjacency operator A restricted to the local neighborhood
has the same eigenvalues as the full A, but different weights.

Key: A^n[v0, v] = number of walks of length n from v0 to v.
This IS the n-th level of self-reflection: how many ways
can you reach vertex v in exactly n steps from v0?
""")

# Compute A^n[v0, :] for various n
An = np.eye(N)
print(f"{'n':>4} {'||A^n[v0,:]||':>14} {'A^n[v0,v0]':>12} {'sum':>10} {'log_phi(sum)':>13}")

walk_from_v0 = []
for n in range(20):
    row = An[v0, :]
    norm_row = np.linalg.norm(row)
    total = np.sum(row)
    diag = An[v0, v0]
    log_total = np.log(total) / np.log(phi) if total > 0 else 0

    walk_from_v0.append(row.copy())

    print(f"  {n:3d} {norm_row:14.2f} {diag:12.2f} {total:10.0f} {log_total:13.4f}")
    An = An @ adj

# ============================================================
# SECTION 5: Projection of walks onto irrep eigenvectors
# ============================================================
print("\n" + "=" * 70)
print("SECTION 5: WALK PROJECTION ONTO IRREPS")
print("=" * 70)

print("""
A^n[v0,:] projected onto each irrep's eigenspace gives:
  w_r(n) = sum_{modes m in rep r} (A^n[v0,:] . evec_m)^2

For irrep r with eigenvalue lambda_r:
  w_r(n) should be proportional to lambda_r^{2n}
  (since A^n has eigenvalue lambda_r^n)

The phi-irrep (r1, lambda=6*phi) gives weight ~ (6*phi)^{2n}
""")

print(f"{'n':>4}", end="")
for r in range(9):
    print(f"  r{r}({rep_groups[r][2]})", end="")
print()

for n in range(min(8, len(walk_from_v0))):
    row = walk_from_v0[n]
    print(f"  {n:3d}", end="")
    for r in range(9):
        # Project onto eigenspace r
        weight = 0
        for m in rep_groups[r][1]:
            coeff = np.dot(row, evecs_sorted[:, m])
            weight += coeff**2
        # Normalize by eigenvalue^(2n) to see the structure
        print(f"  {weight:6.1f}", end="")
    print()

# ============================================================
# SECTION 6: THE KEY TEST - Does phi^n emerge?
# ============================================================
print("\n" + "=" * 70)
print("SECTION 6: DOES phi^n EMERGE FROM SELF-REFLECTION?")
print("=" * 70)

# The walk A^n from v0, restricted to v0 itself, gives A^n[v0,v0]
# = number of closed walks of length n.
# In terms of eigenvalues: A^n[v0,v0] = sum_k lambda_k^n * |v_k(v0)|^2

# The dominant eigenvalue is 12 (trivial rep), so A^n[v0,v0] ~ 12^n / 120

# But the RATIO of the r1 component to the r0 component:
# (6*phi)^n * |v_r1(v0)|^2 / (12^n * |v_r0(v0)|^2)
# = (phi/2)^n * (4/120) / (1/120) = (phi/2)^n * 4

# So the phi-component is SUPPRESSED by (phi/2)^n relative to the trivial.
# This goes to ZERO, not to phi^n!

# To get phi^n, we need to ISOLATE the phi-eigenvalue component.
# On the McKay graph (not the 600-cell), the adjacency eigenvalue IS phi.
# On the 600-cell, it's 6*phi (scaled by the degree 12).

# McKay walk analysis
print("\n--- McKay graph walk analysis ---")
M_mckay = np.zeros((9, 9))
for a, b in [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (5,7), (7,8)]:
    M_mckay[a, b] = 1; M_mckay[b, a] = 1

evals_mk, evecs_mk = np.linalg.eigh(M_mckay)
idx_mk = np.argsort(evals_mk)[::-1]
evals_mk = evals_mk[idx_mk]
evecs_mk = evecs_mk[:, idx_mk]

# Find phi-eigenvector index
phi_idx = np.argmin(np.abs(evals_mk - phi))
v_phi = evecs_mk[:, phi_idx]

print(f"McKay eigenvalue: {evals_mk[phi_idx]:.6f} = phi")
print(f"phi-eigenvector: {v_phi.round(4)}")
print()

# Walk M^n * e_0 projected onto v_phi
print(f"{'n':>4} {'(M^n*e0).v_phi':>15} {'phi^n*v_phi[0]':>15} {'ratio':>10}")
Mn = np.eye(9)
e0 = np.zeros(9); e0[0] = 1

for n in range(20):
    walk = Mn @ e0
    proj = np.dot(walk, v_phi)
    expected = phi**n * v_phi[0]
    ratio = proj / expected if abs(expected) > 1e-15 else 0

    print(f"  {n:3d} {proj:15.6f} {expected:15.6f} {ratio:10.6f}")
    Mn = Mn @ M_mckay

# ============================================================
# SECTION 7: SPECTRAL UNFOLDING - mass from walk
# ============================================================
print("\n" + "=" * 70)
print("SECTION 7: MASS FROM SELF-REFLECTION WALK")
print("=" * 70)

print("""
KEY FORMULA: (M^n * e_0) . v_phi = phi^n * v_phi[0]

This means: the phi-component of the n-step walk from r0 IS phi^n.
The mass formula m_f = m_e * phi^(n_f) says:
  "The mass of fermion f is the phi-projection of n_f steps of self-reflection."

So the mass exponent n_f = number of self-reflection steps to reach fermion f.

But WHERE on the McKay graph is each fermion? That's the missing link.
Let's check: does M^n * e_0 at NODE j give phi^(n_j)?
""")

# For each node j, at which n does (M^n * e_0)[j] first become significant?
# And what is the phi-projection of (M^n * e_0)?

print(f"Walk M^n * e_0 at each node:")
print(f"{'n':>4}", end="")
for j in range(9):
    print(f"  r{j}({rep_groups[j][2]})", end="")
print()

Mn = np.eye(9)
for n in range(10):
    walk = Mn @ e0
    print(f"  {n:3d}", end="")
    for j in range(9):
        val = walk[j]
        if abs(val) < 0.5:
            print(f"  {val:6.1f}", end="")
        else:
            print(f"  {val:6.0f}", end="")
    print()
    Mn = Mn @ M_mckay

# The walk at node j at step n = (M^n)[0,j] = number of paths of length n from 0 to j
# Using spectral decomposition: (M^n)[0,j] = sum_k lambda_k^n * v_k(0) * v_k(j)

# For large n, dominated by lambda_max = 2:
# (M^n)[0,j] ~ 2^n * v_Perron(0) * v_Perron(j) = 2^n * dim(0)/(sum dims) * dim(j)/(sum dims)

# The phi-component: phi^n * v_phi(0) * v_phi(j)
# This is ALWAYS present but overwhelmed by the 2^n term.

# To ISOLATE the phi-component, subtract the Perron-Frobenius term:
print(f"\nphi-component isolated: (M^n)[0,j] - 2^n * d_0 * d_j / sum(d^2)")
print(f"{'n':>4}", end="")
for j in range(9):
    print(f"  r{j}({rep_groups[j][2]})", end="")
print()

dim_vec = np.array([r[2] for r in rep_groups], dtype=float)
dim_norm2 = np.sum(dim_vec**2)

Mn = np.eye(9)
for n in range(10):
    walk = Mn @ e0
    print(f"  {n:3d}", end="")
    for j in range(9):
        perron = 2**n * dim_vec[0] * dim_vec[j] / dim_norm2
        residual = walk[j] - perron
        print(f"  {residual:6.2f}", end="")
    print()
    Mn = Mn @ M_mckay

# ============================================================
# SECTION 8: SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SECTION 8: SYNTHESIS")
print("=" * 70)

print(f"""
RESULTS:
========

1. BALL EXPANSION on 600-cell:
   All 9 irreps appear within radius 2 (fast convergence).
   This is because the 600-cell is highly connected (degree 12).
   There is NO sequential emergence of irreps - they all appear
   almost simultaneously. The "self-reflection" picture doesn't
   create particles ONE BY ONE on the 600-cell.

2. SHELL SPECTRAL CONTENT:
   The spectral weight in each shell is determined by the
   eigenvector values at the vertices in that shell.
   The trivial rep (r0) has weight 1/120 per vertex (uniform).
   Other reps have varying weights per shell.

3. McKAY WALK - KEY RESULT:
   (M^n * e_0) . v_phi = phi^n * v_phi[0]  EXACTLY.

   This is the spectral decomposition theorem applied to the
   McKay adjacency. It's MATHEMATICALLY EXACT and holds for all n.

   The mass formula m_f = m_e * phi^(n_f) IS the statement:
   "The phi-projection of n_f steps of McKay self-reflection
   gives the mass ratio phi^(n_f)."

4. BUT: the walk M^n * e_0 at a specific NODE j does NOT give phi^(n_j).
   It gives 2^n * dim_j/120 + phi^n * v_phi(0)*v_phi(j) + ...
   The phi-component is always present but mixed with other eigenvalues.

5. THE ISOLATION PROBLEM:
   To get PURE phi^n, you need to project onto the phi-eigenvector.
   This is a GLOBAL operation (needs all nodes), not a LOCAL one.
   The "self-reflection" at a SINGLE point doesn't isolate phi^n.

CONCLUSION:
===========
The mathematical identity (M^n * e_0) . v_phi = phi^n * v_phi[0]
is CORRECT and EXACT. It connects n steps of McKay self-reflection
to the mass ratio phi^n. But:

(a) It works on the McKay graph (9 nodes), not on the 600-cell (120 nodes)
(b) It requires projection onto the phi-eigenvector (a global operation)
(c) It doesn't explain WHY each fermion has its specific exponent n_f

The "self-reflection" intuition is VALID as a metaphor for the McKay walk,
but the mass hierarchy (the specific values of n_f) still comes from
the Z[phi] selection rules and Wilson line structure, not from the walk itself.

WHAT WOULD MAKE THIS DEEPER:
If we could show that the phi-projection is NATURAL (not arbitrary) -
for example, if the 600-cell geometry SELECTS the phi-eigenvector
as the unique physically relevant direction - then the mass formula
would truly be "n levels of self-reflection give mass phi^n."

The Galois kernel theorem DOES select rho_1 (which has eigenvalue 6*phi
on the 600-cell and phi on the McKay graph) as special. This is suggestive
but not yet a complete argument.
""")

print("=" * 70)
print("EXP-492 COMPLETE")
print("=" * 70)
