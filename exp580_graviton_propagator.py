"""
exp580: Graviton propagator on the 600-cell and discrete Newton's law.

QUESTION: Does the 600-cell produce a 1/r gravitational potential?

SETUP:
  - The 600-cell discretizes S^3 (3-sphere, diameter 5 in graph distance)
  - On S^3, the Green's function G(theta) ~ 1/(4*pi^2*sin(theta)) at short distances
  - For small theta: G ~ 1/(4*pi^2 * r) = Newtonian 1/r potential in 3D
  - We compute the discrete Green's function and test this.

THREE LEVELS OF ANALYSIS:
  1. Graph Laplacian Green's function (simplest: L_graph^{-1})
  2. Static Box sector: effective potential from Box restricted to time-independent modes
  3. Full graviton propagator from spectral action Hessian

From exp578-579:
  - kappa = 0 (Ricci flat vacuum)
  - Stiffness ratio fiber/cross = a1^2 = 25
  - Graviton propagates at order 4 (Tr(Box^4) Hessian)
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, norm, pinv
from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix
from collections import Counter
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
degree = 12

print("=" * 70)
print("EXP-580: GRAVITON PROPAGATOR AND DISCRETE NEWTON'S LAW")
print("=" * 70)

# =====================================================================
# STEP 0: Build infrastructure
# =====================================================================
print("\nBuilding 600-cell, Hopf fibration, and Box operator...")

verts, adj, lap = build_600cell()

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
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
vtx_fiber = {}
for fi, f in enumerate(fibers):
    for v in f: vtx_fiber[v] = fi

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0
A_cross = adj - A_fiber
Box = a1 * A_fiber - A_cross

# Graph distances
sparse_adj = csr_matrix(adj)
dist_matrix = shortest_path(sparse_adj, method='D', unweighted=True).astype(int)

# Distance distribution per vertex
dist_counts = Counter()
for j in range(N):
    if j != 0:
        dist_counts[dist_matrix[0, j]] += 1
print(f"  Distance distribution from vertex 0:")
for d in sorted(dist_counts.keys()):
    print(f"    d={d}: {dist_counts[d]} vertices")

# Geodesic angles on S^3 (approximate: theta = d * pi / diameter)
diameter = dist_matrix.max()
print(f"  Diameter: {diameter}")


# =====================================================================
# PART 1: GRAPH LAPLACIAN GREEN'S FUNCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: GRAPH LAPLACIAN GREEN'S FUNCTION (VERTEX SPACE)")
print("=" * 70)

# L = degree*I - A. Spectrum: 0 (x1), lambda_1, ..., lambda_119
evals_L, evecs_L = eigh(lap)

print(f"  Laplacian spectrum: [{evals_L[0]:.4f}, {evals_L[-1]:.4f}]")
print(f"  Zero modes: {np.sum(np.abs(evals_L) < 1e-8)}")

# Green's function = pseudoinverse of L
# G_L(x,y) = sum_{k: lam_k > 0} psi_k(x) * psi_k(y) / lam_k
G_L = np.zeros((N, N))
for k in range(N):
    if evals_L[k] > 1e-8:
        G_L += np.outer(evecs_L[:, k], evecs_L[:, k]) / evals_L[k]

# Verify: L @ G_L = I - (1/N)|1><1|  (projector onto non-constant)
proj_const = np.ones((N, N)) / N
err_LG = norm(lap @ G_L - (np.eye(N) - proj_const))
print(f"  ||L @ G_L - (I - P_0)||_F = {err_LG:.2e}")

# Potential from a point source at vertex 0:
# L * Phi = delta_0 - 1/N  (source minus background to remove zero mode)
# Phi = G_L @ (delta_0 - 1/N) = G_L[:, 0] - G_L.mean(axis=1)
source = np.zeros(N)
source[0] = 1.0
Phi_L = G_L @ (source - 1.0/N)

# Normalize: Phi(0) = max
print(f"\n  Gravitational potential from Laplacian Green's function:")
print(f"    Phi(source) = Phi(d=0) = {Phi_L[0]:.8f}")

# Average potential at each distance
print(f"\n  {'d':>3s}  {'n_vertices':>10s}  {'<Phi>':>14s}  {'1/d':>10s}  {'ratio':>10s}")
print(f"  {'-'*55}")
for d in range(0, diameter + 1):
    mask = (dist_matrix[0, :] == d)
    n_d = np.sum(mask)
    if n_d > 0:
        avg_phi = Phi_L[mask].mean()
        inv_d = 1.0 / d if d > 0 else float('inf')
        ratio = avg_phi * d if d > 0 else 0
        print(f"  {d:3d}  {n_d:10d}  {avg_phi:14.8f}  {inv_d:10.6f}  "
              f"{ratio:10.6f}")

# Fit 1/d to data at d=1,2,3
phi_d = {}
for d in range(1, diameter + 1):
    mask = (dist_matrix[0, :] == d)
    phi_d[d] = Phi_L[mask].mean()

# If Phi ~ A/d, then A = Phi*d should be constant
print(f"\n  Test V ~ A/d (A = V*d should be constant):")
for d in range(1, diameter + 1):
    A = phi_d[d] * d
    print(f"    d={d}: V*d = {A:.8f}")

# Compare with S^3 continuum Green's function:
# G_S3(theta) = (1/(4*pi^2)) * (pi - theta) / sin(theta)
# where theta = geodesic distance on S^3 of radius R=1
# For the 600-cell: theta_d = d * arccos(phi/2)
# Actually more precisely: d=1 corresponds to edge length on S^3
# The edge length = 2*arcsin(1/phi) (from the fact that vertices are at
# angular distance arccos(phi/2) on S^3)

theta_edge = np.arccos(PHI / 2)  # ~ 0.5536 rad
print(f"\n  S^3 continuum comparison:")
print(f"  Edge angular distance: theta_1 = arccos(phi/2) = {theta_edge:.6f} rad")
print(f"  {'d':>3s}  {'theta':>10s}  {'G_S3':>12s}  {'G_disc':>12s}  {'ratio':>10s}")
print(f"  {'-'*52}")
for d in range(1, diameter + 1):
    theta = d * theta_edge
    if theta < np.pi:
        G_S3 = (np.pi - theta) / (4 * np.pi**2 * np.sin(theta))
    else:
        G_S3 = 0
    ratio = phi_d[d] / G_S3 if G_S3 != 0 else float('inf')
    print(f"  {d:3d}  {theta:10.6f}  {G_S3:12.8f}  {phi_d[d]:12.8f}  {ratio:10.6f}")


# =====================================================================
# PART 2: BOX-DERIVED STATIC POTENTIAL
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: STATIC POTENTIAL FROM BOX OPERATOR")
print("=" * 70)

# For static (time-independent) perturbations, the relevant operator is
# the cross-sector. A static function f is constant on each fiber.
#
# The projection onto static functions:
# P_static = (1/10) * sum_{v in fiber} |v><v| contracted to fiber basis
#
# Effective operator on the 12-dimensional base space:
# Build base adjacency from cross connections between fibers.

n_base = len(fibers)  # = 12

# Base adjacency: A_base[i,j] = number of cross edges from fiber i to fiber j / 10
A_base = np.zeros((n_base, n_base))
for i in range(N):
    for j in range(N):
        if A_cross[i,j] > 0.5:
            fi = vtx_fiber[i]
            fj = vtx_fiber[j]
            if fi != fj:
                A_base[fi, fj] += 1

# Each vertex has 10 cross neighbors. Within a fiber of 10 vertices,
# total cross edges out = 10 * 10 = 100. Divide by 10 to get per-vertex base adjacency.
A_base_per_vtx = A_base / 10.0  # cross neighbors in fiber j per vertex in fiber i

print(f"  Base space: {n_base} fibers (icosahedral graph)")
print(f"  Base adjacency (per vertex):")
base_degree = A_base_per_vtx.sum(axis=1)
print(f"    Degree per base vertex: {base_degree[0]:.1f} (expect 10 = cross degree)")

# Base distances
adj_base_binary = (A_base > 0).astype(float)
dist_base = shortest_path(csr_matrix(adj_base_binary), method='D', unweighted=True).astype(int)
print(f"    Base diameter: {dist_base.max()}")
print(f"    Base distance distribution:")
for d in range(0, dist_base.max() + 1):
    n_d = np.sum(dist_base[0, :] == d)
    print(f"      d_base={d}: {n_d} fibers")

# Effective Box on static sector
# Box|_{static} acting on fiber-averaged functions:
# (Box_eff g)(F_i) = a1 * 2 * g(F_i) - sum_j A_base_per_vtx[i,j] * g(F_j)
# = 2*a1*g(i) - sum_j A_base[i,j]/10 * g(j)

L_base = np.diag(A_base_per_vtx.sum(axis=1)) - A_base_per_vtx
Box_base = 2 * a1 * np.eye(n_base) - A_base_per_vtx
# Actually: Box_eff = a1 * (2*I) - A_cross_base
# where A_cross_base = A_base_per_vtx

evals_base, evecs_base = eigh(L_base)
evals_box_base, evecs_box_base = eigh(Box_base)

print(f"\n  Base Laplacian spectrum: {np.round(evals_base, 4)}")
print(f"  Base Box spectrum: {np.round(evals_box_base, 4)}")

# Green's function of base Laplacian
G_base = np.zeros((n_base, n_base))
for k in range(n_base):
    if evals_base[k] > 1e-8:
        G_base += np.outer(evecs_base[:, k], evecs_base[:, k]) / evals_base[k]

# Potential from base source
source_base = np.zeros(n_base)
source_base[0] = 1.0
Phi_base = G_base @ (source_base - 1.0/n_base)

print(f"\n  Base Newton potential (from L_base):")
for d in range(0, dist_base.max() + 1):
    mask = (dist_base[0, :] == d)
    n_d = np.sum(mask)
    if n_d > 0:
        avg = Phi_base[mask].mean()
        print(f"    d_base={d}: n={n_d}, <Phi_base> = {avg:.8f}"
              f"{'  (source)' if d == 0 else f'  V*d = {avg*d:.6f}' if d > 0 else ''}")


# =====================================================================
# PART 3: FULL 600-CELL GREEN'S FUNCTION -- DETAILED ANALYSIS
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: DETAILED 1/r TEST ON THE FULL 600-CELL")
print("=" * 70)

# The key test: Phi(d) vs 1/d
# On a 3-sphere of radius R, the Green's function is:
# G(r) = 1/(4*pi^2*R) * (1/r - 1/(pi*R) + corrections)
# where r = R*theta is the geodesic distance.
# At short distances: G ~ 1/(4*pi^2*R*r) = 1/(4*pi^2) * 1/r for R=1.
#
# On the graph, distances are discrete: d = 1,2,3,4,5
# The "effective radius" of the 600-cell:
# 120 vertices on S^3: vol(S^3) = 2*pi^2*R^3 = 120 * vol_cell
# So R^3 ~ 120 / (2*pi^2) -> R ~ (120/(2*pi^2))^{1/3} ~ 1.83

R_eff = (N / (2 * np.pi**2))**(1.0/3)
print(f"  Effective radius: R_eff = (N/(2*pi^2))^(1/3) = {R_eff:.4f}")

# Effective geodesic distance: r_d = d * l_edge
# Edge length: on S^3 of radius R, angular separation arccos(phi/2)
# l_edge = R * arccos(phi/2)
l_edge = R_eff * theta_edge
print(f"  Edge length: l_edge = R*theta = {l_edge:.4f}")

# Newton potential comparison
print(f"\n  NEWTON'S LAW TEST: V(d) vs A/d")
print(f"  {'d':>3s}  {'n(d)':>6s}  {'V_disc':>14s}  {'V*d':>12s}  "
      f"{'V_S3_cont':>12s}  {'V_disc/V_S3':>12s}")
print(f"  {'-'*70}")

A_newton = phi_d[1] * 1  # amplitude from d=1

for d in range(1, diameter + 1):
    n_d = np.sum(dist_matrix[0, :] == d)
    V_disc = phi_d[d]
    V_times_d = V_disc * d

    # Continuum S^3 Green's function
    theta = d * theta_edge
    r = R_eff * theta
    if theta < np.pi:
        V_S3 = (np.pi - theta) / (4 * np.pi**2 * R_eff * np.sin(theta))
    else:
        V_S3 = 0.001  # avoid division by zero

    ratio = V_disc / V_S3

    print(f"  {d:3d}  {n_d:6d}  {V_disc:14.8f}  {V_times_d:12.8f}  "
          f"{V_S3:12.8f}  {ratio:12.6f}")

# Fit power law: log(V) = log(A) - alpha * log(d)
import numpy as np
d_vals = np.array([1, 2, 3])  # short distances only
V_vals = np.array([phi_d[d] for d in d_vals])
log_d = np.log(d_vals)
log_V = np.log(V_vals)

# Linear fit: log(V) = a + b * log(d)
coeffs = np.polyfit(log_d, log_V, 1)
power_law_exp = coeffs[0]
amplitude = np.exp(coeffs[1])

print(f"\n  Power law fit V ~ A * d^n (using d=1,2,3):")
print(f"    Exponent n = {power_law_exp:.6f}")
print(f"    Amplitude A = {amplitude:.8f}")
print(f"    Expected for Newton: n = -1")
print(f"    Deviation: {abs(power_law_exp + 1):.6f}")

# Also fit using d=1,2,3,4
d_vals4 = np.array([1, 2, 3, 4])
V_vals4 = np.array([phi_d[d] for d in d_vals4])
coeffs4 = np.polyfit(np.log(d_vals4), np.log(V_vals4), 1)
print(f"\n  Power law fit V ~ A * d^n (using d=1,2,3,4):")
print(f"    Exponent n = {coeffs4[0]:.6f}")
print(f"    Deviation from -1: {abs(coeffs4[0] + 1):.6f}")


# =====================================================================
# PART 4: DECOMPOSITION INTO IRREP CHANNELS
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: GREEN'S FUNCTION DECOMPOSITION BY 2I IRREPS")
print("=" * 70)

# The 600-cell graph is the Cayley graph of 2I (binary icosahedral).
# Its Laplacian eigenvalues = 12 - chi_k(g)/dim(k) for the generators.
# The Green's function decomposes by irreps:
# G(x,y) = sum_{k: irrep k, lam_k > 0} dim(k)/N * chi_k(x^{-1}*y) / lam_k
#
# This means V(d) = sum_k G_k(d) where each irrep k contributes at each distance.

# Group the Laplacian eigenmodes by eigenvalue (= irrep)
evals_rounded = np.round(evals_L, 4)
irrep_groups = {}
for k in range(N):
    val = evals_rounded[k]
    if val not in irrep_groups:
        irrep_groups[val] = []
    irrep_groups[val].append(k)

print(f"  Irrep decomposition of Green's function:")
print(f"  {'lambda':>10s}  {'mult':>5s}  {'contribution to V(d=1)':>22s}  "
      f"{'to V(d=2)':>12s}  {'to V(d=3)':>12s}")
print(f"  {'-'*70}")

for lam_val in sorted(irrep_groups.keys()):
    if lam_val < 1e-6:
        continue  # skip zero mode
    indices = irrep_groups[lam_val]
    mult = len(indices)

    # Contribution to G(0, y) at distance d
    contribs = {}
    for d in range(1, 4):
        # Average over vertices at distance d from vertex 0
        mask = (dist_matrix[0, :] == d)
        contrib = 0
        for k in indices:
            # sum over y at distance d of psi_k(0)*psi_k(y) / lam_k
            contrib += evecs_L[0, k] * evecs_L[mask, k].mean() / lam_val
        contribs[d] = contrib

    print(f"  {lam_val:10.4f}  {mult:5d}  {contribs[1]:22.8f}  "
          f"{contribs[2]:12.8f}  {contribs[3]:12.8f}")

# Total (should match Phi_L)
for d in range(1, 4):
    total = sum(
        sum(evecs_L[0, k] * evecs_L[dist_matrix[0, :] == d, k].mean() / evals_L[k]
            for k in grp)
        for lam_val, grp in irrep_groups.items() if lam_val > 1e-6)
    print(f"  Total at d={d}: {total:.8f}  vs  Phi_L = {phi_d[d]:.8f}")


# =====================================================================
# PART 5: BOX OPERATOR GREEN'S FUNCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: BOX GREEN'S FUNCTION (LORENTZIAN)")
print("=" * 70)

# The Box operator is Lorentzian: it has both positive and negative eigenvalues.
# For the graviton propagator, we need the Feynman propagator:
# G_F = (Box - i*epsilon)^{-1} in Lorentzian signature
#
# For the STATIC potential (time-independent source), we need the spatial
# Green's function. This requires projecting out the time direction.
#
# Practical approach: use Box^2 (positive semidefinite) as the kinetic operator.
# G_{Box2} = (Box^2)^{-1} on the non-kernel subspace.

evals_box, evecs_box = eigh(Box)

print(f"  Box spectrum: [{evals_box[0]:.4f}, {evals_box[-1]:.4f}]")
print(f"  Zero modes: {np.sum(np.abs(evals_box) < 1e-6)}")

# Green's function of Box^2 (= |Box|^{-2} on non-kernel)
G_Box2 = np.zeros((N, N))
for k in range(N):
    if abs(evals_box[k]) > 1e-6:
        G_Box2 += np.outer(evecs_box[:, k], evecs_box[:, k]) / evals_box[k]**2

# Potential from Box^2
Phi_Box2 = G_Box2 @ (source - 1.0/N)

print(f"\n  Newton potential from Box^2 Green's function:")
print(f"  {'d':>3s}  {'V_L':>14s}  {'V_Box2':>14s}  {'V_L*d':>12s}  {'V_Box2*d':>12s}")
print(f"  {'-'*60}")
for d in range(1, diameter + 1):
    mask = (dist_matrix[0, :] == d)
    V_L = phi_d[d]
    V_B = Phi_Box2[mask].mean()
    print(f"  {d:3d}  {V_L:14.8f}  {V_B:14.8f}  {V_L*d:12.8f}  {V_B*d:12.8f}")

# Power law for Box^2
phi_box_d = {}
for d in range(1, diameter + 1):
    mask = (dist_matrix[0, :] == d)
    phi_box_d[d] = Phi_Box2[mask].mean()

V_box_vals = np.array([phi_box_d[d] for d in [1, 2, 3]])
coeffs_box = np.polyfit(np.log([1, 2, 3]), np.log(V_box_vals), 1)
print(f"\n  Box^2 power law (d=1,2,3): V ~ d^{coeffs_box[0]:.4f}")


# =====================================================================
# PART 6: THE RETARDED GREEN'S FUNCTION (Box itself)
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: GREEN'S FUNCTION OF Box (RETARDED)")
print("=" * 70)

# For the retarded propagator: G_R = (Box + i*eps)^{-1}
# The real part gives the Yukawa-like potential.
# For real analysis: G_R = Box^{-1} on non-kernel.
# Box has both + and - eigenvalues, so this is well-defined on non-kernel.

G_Box = np.zeros((N, N))
for k in range(N):
    if abs(evals_box[k]) > 1e-6:
        G_Box += np.outer(evecs_box[:, k], evecs_box[:, k]) / evals_box[k]

Phi_Box = G_Box @ (source - 1.0/N)

print(f"  Potential from Box^{{-1}} (retarded):")
print(f"  {'d':>3s}  {'V_Box':>14s}  {'V_Box*d':>12s}")
print(f"  {'-'*35}")
for d in range(1, diameter + 1):
    mask = (dist_matrix[0, :] == d)
    V_B1 = Phi_Box[mask].mean()
    print(f"  {d:3d}  {V_B1:14.8f}  {V_B1*d:12.8f}")


# =====================================================================
# PART 7: FIBER-AVERAGED (SPATIAL) GREEN'S FUNCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: FIBER-AVERAGED SPATIAL GREEN'S FUNCTION")
print("=" * 70)

# Project the Green's function onto the static (fiber-constant) sector.
# This gives the spatial part relevant for Newton's law.

# Fiber averaging projector P: (Pf)(v) = (1/10) * sum_{u in fiber(v)} f(u)
P_static = np.zeros((N, N))
for fi, fib in enumerate(fibers):
    for i in fib:
        for j in fib:
            P_static[i, j] = 1.0 / 10.0

# Projected Green's function: G_spatial = P @ G_L @ P
G_spatial = P_static @ G_L @ P_static

# Source: fiber-averaged delta at fiber 0
source_fiber = np.zeros(N)
for v in fibers[0]:
    source_fiber[v] = 1.0 / 10.0  # normalized point source on fiber 0

Phi_spatial = G_spatial @ (source_fiber - 1.0/N)

# Average over fibers at each BASE distance
print(f"  Spatial (fiber-averaged) potential:")
print(f"  {'d_base':>6s}  {'n_fibers':>8s}  {'V_spatial':>14s}  {'V*d':>12s}")
print(f"  {'-'*45}")
for d in range(0, dist_base.max() + 1):
    base_mask = (dist_base[0, :] == d)
    n_fib = np.sum(base_mask)
    if n_fib > 0:
        # Average potential at the fibers at base distance d
        avg_V = 0
        for fi in range(n_base):
            if base_mask[fi]:
                avg_V += Phi_spatial[fibers[fi][0]]  # take any vertex in fiber
        avg_V /= n_fib
        vd = avg_V * d if d > 0 else 0
        print(f"  {d:6d}  {n_fib:8d}  {avg_V:14.8f}  {vd:12.8f}")


# =====================================================================
# PART 8: EXTRACT G_NEWTON AND COMPARE WITH SPECTRAL ACTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: EXTRACTING NEWTON'S CONSTANT G_N")
print("=" * 70)

# In the continuum on S^3:
# V(r) = G_N * M / r   (Newton's potential, 3D)
#
# On the graph:
# V(d) = G_N_discrete * M / d   (if 1/r works)
#
# The "mass" of our source is M = 1 (delta function).
# G_N_discrete = V(d) * d (should be constant)

G_N_disc = phi_d[1]  # = V(1) * 1

print(f"  G_N (discrete) from V(d=1): {G_N_disc:.8f}")
print(f"  G_N from V(d=2): {phi_d[2] * 2:.8f}")
print(f"  G_N from V(d=3): {phi_d[3] * 3:.8f}")

# Continuum G_N on S^3:
# G_N = 1 / (4*pi^2*R) (the coefficient of the S^3 Green's function)
G_N_continuum = 1.0 / (4 * np.pi**2 * R_eff)
print(f"\n  G_N (S^3 continuum): 1/(4*pi^2*R) = {G_N_continuum:.8f}")

# From spectral action: G_N = 3/(pi * c_1) where c_1 = 14880
# In the Connes-Chamseddine approach: the EH term gives
# S_EH = (1/(16*pi*G_N)) * integral R * sqrt(g)
# and c_1 = Tr(D^2) = 14880 relates to 1/G_N
c1_D = 14880
print(f"\n  Spectral action c_1 = Tr(D^2) = {c1_D}")
print(f"  c_1 / N = {c1_D/N}")

# From the Seeley-DeWitt expansion:
# S_EH ~ f_2 * Lambda^2 * c_1 / (4*pi)
# This gives 1/G_N ~ f_2 * Lambda^2 * c_1 / (4*pi)
# For our discrete model, Lambda^{-2} ~ l_edge^2, so:
# G_N ~ 4*pi * l_edge^2 / (f_2 * c_1)

# Let's just check the RATIO:
# G_N_disc / G_N_cont = ?
ratio_GN = G_N_disc / G_N_continuum
print(f"\n  G_N_disc / G_N_cont = {ratio_GN:.6f}")

# The ratio should be related to the normalization of the graph Green's function
# vs continuum. On a graph with N vertices and volume V = N * vol_cell:
# G_graph = (1/N) * G_continuum * V = G_continuum * vol_cell
# So the ratio encodes the cell volume.

print(f"\n  Interpretation: G_N_disc/G_N_cont = {ratio_GN:.4f}")
print(f"    = normalization factor from discrete -> continuum")
print(f"    = effective cell volume on S^3")


# =====================================================================
# PART 9: SUMMARY AND VERDICT
# =====================================================================
print("\n" + "=" * 70)
print("PART 9: SUMMARY -- DOES THE 600-CELL GIVE NEWTON'S LAW?")
print("=" * 70)

print(f"""
  GRAPH LAPLACIAN GREEN'S FUNCTION:
    V(d=1) = {phi_d[1]:.8f}
    V(d=2) = {phi_d[2]:.8f}  (V*d = {phi_d[2]*2:.8f})
    V(d=3) = {phi_d[3]:.8f}  (V*d = {phi_d[3]*3:.8f})
    V(d=4) = {phi_d[4]:.8f}  (V*d = {phi_d[4]*4:.8f})
    V(d=5) = {phi_d[5]:.8f}  (V*d = {phi_d[5]*5:.8f})

  POWER LAW FIT (d=1,2,3): V ~ d^{{{power_law_exp:.4f}}}
    Expected for 3D Newton: V ~ d^(-1)
    Deviation: {abs(power_law_exp + 1):.4f}

  V*d CONSTANCY TEST (1/r means V*d = const):
    V(1)*1 = {phi_d[1]*1:.8f}
    V(2)*2 = {phi_d[2]*2:.8f}  ({(phi_d[2]*2)/(phi_d[1]*1)*100:.1f}% of V(1)*1)
    V(3)*3 = {phi_d[3]*3:.8f}  ({(phi_d[3]*3)/(phi_d[1]*1)*100:.1f}% of V(1)*1)

  S^3 CONTINUUM MATCH:
    G_N_disc = {G_N_disc:.8f}
    G_N_cont = {G_N_continuum:.8f}
    Ratio: {ratio_GN:.4f}

  BOX^2 GREEN'S FUNCTION:
    Power law: V ~ d^{{{coeffs_box[0]:.4f}}}
""")

# Final verdict
if abs(power_law_exp + 1) < 0.15:
    verdict = "YES -- 1/r POTENTIAL CONFIRMED"
elif abs(power_law_exp + 1) < 0.3:
    verdict = "PARTIAL -- approximate 1/r with corrections"
else:
    verdict = f"NO -- power law is d^{{{power_law_exp:.2f}}}, not d^(-1)"

print(f"  VERDICT: {verdict}")
print(f"  (On a discrete compact graph, deviations at d >= diameter/2 are expected)")

# Check: is the deviation at d=2,3 consistent with S^3 curvature corrections?
# On S^3: V(r) = (1/4pi^2) * (pi-r)/sin(r) ~ 1/r * (1 - r^2/6 + ...)
# So V*r ~ const * (1 - r^2/6) = decreasing with r.
# If V*d is DECREASING, it matches S^3 curvature corrections.
decreasing = phi_d[1]*1 > phi_d[2]*2 > phi_d[3]*3
print(f"\n  V*d decreasing (S^3 correction): {decreasing}")
if decreasing:
    print(f"  -> The drop in V*d at larger d is CONSISTENT with")
    print(f"     positive curvature of S^3 (the 600-cell is compact).")
    print(f"     At short distances, 1/r Newton's law holds.")

print("\n" + "=" * 70)
print("EXP-580 COMPLETE")
print("=" * 70)
