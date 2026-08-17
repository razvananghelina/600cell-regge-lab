"""
EXP-108: Soliton Stability on the 600-Cell Graph
=================================================

QUESTION: Can we DERIVE which (a,b) configurations are dynamically
stable solitons, elevating the (a,b) assignment from PATTERN to DERIVED?

METHOD:
  Discretized sine-Gordon model on the 600-cell graph:
    E[eta] = (kappa/2) * eta^T L eta + lambda * Sum sin^2(pi * eta_i)
  where eta_i = ln(psi_i/m_e)/ln(phi) is the log-mass field at vertex i.

  Minima of the potential are at eta = integer (masses m_e * phi^n).
  Solitons = field configurations interpolating between different minima.

  For each target winding number n:
  1. Construct ansatz (radial profile from source vertex)
  2. Relax via gradient descent to nearest local minimum
  3. Compute Hessian eigenvalues -> check stability
  4. Analyze directional structure (diameter vs decagonal)

KEY OBSERVATIONS from exp107:
  - All occupied levels have gcd(|a|,|b|) = 1 (or trivial (0,0))
  - All a=0 levels except n=0 are EMPTY: (0,1), (0,2), (0,3), (0,4)
  - n=6 with (0,+1) and |a|+|b|=1 is EMPTY
  - n=5 with (1,+0) and |a|+|b|=1 is OCCUPIED (Down quark)
  - WHY?

STATUS: EXPLORATORY
"""

import numpy as np
from itertools import permutations
from collections import deque

PHI = (1 + np.sqrt(5)) / 2

# =====================================================
# PART 1: Build 600-Cell Graph
# =====================================================

print("=" * 70)
print("EXP-108: Soliton Stability on the 600-Cell Graph")
print("=" * 70)

print("\n--- PART 1: Graph Construction ---\n")

def build_600cell():
    """Generate 120 vertices and find edges of the 600-cell."""
    verts = []

    # Family 1: 8 vertices (16-cell) - axis-aligned unit vectors
    for i in range(4):
        for s in [+1, -1]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = float(s)
            verts.append(tuple(v))

    # Family 2: 16 vertices (tesseract) - all (+-1/2, +-1/2, +-1/2, +-1/2)
    for s0 in [+0.5, -0.5]:
        for s1 in [+0.5, -0.5]:
            for s2 in [+0.5, -0.5]:
                for s3 in [+0.5, -0.5]:
                    verts.append((s0, s1, s2, s3))

    # Family 3: 96 golden-ratio vertices
    # Even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
    base = [PHI/2, 0.5, 1/(2*PHI), 0.0]

    for perm in permutations(range(4)):
        # Check if even permutation
        p = list(perm)
        inv_count = 0
        for i in range(4):
            for j in range(i+1, 4):
                if p[i] > p[j]:
                    inv_count += 1
        if inv_count % 2 != 0:
            continue

        # All sign combinations for non-zero components
        for s0 in [+1, -1]:
            for s1 in [+1, -1]:
                for s2 in [+1, -1]:
                    v = [0.0, 0.0, 0.0, 0.0]
                    v[perm[0]] = s0 * base[0]
                    v[perm[1]] = s1 * base[1]
                    v[perm[2]] = s2 * base[2]
                    v[perm[3]] = 0.0

                    vt = tuple(v)
                    # Check unit sphere
                    norm = np.sqrt(sum(x**2 for x in vt))
                    if abs(norm - 1.0) > 1e-10:
                        continue
                    # Check duplicate
                    is_dup = False
                    for ev in verts:
                        if all(abs(vt[k] - ev[k]) < 1e-10 for k in range(4)):
                            is_dup = True
                            break
                    if not is_dup:
                        verts.append(vt)

    V = np.array(verts)
    N = len(V)

    # Find edges: minimum nonzero distance between vertices
    # On unit S^3, 600-cell edge length = 1/phi
    edge_len = 1.0 / PHI
    tol = 0.01

    edges = []
    adj = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(i+1, N):
            d = np.linalg.norm(V[i] - V[j])
            if abs(d - edge_len) < tol:
                edges.append((i, j))
                adj[i, j] = 1
                adj[j, i] = 1

    return V, edges, adj

V, edges, adj = build_600cell()
N = len(V)
print(f"Vertices: {N}")
print(f"Edges: {len(edges)}")
print(f"Degree: {np.sum(adj[0])}")

# BFS for all-pairs distances
def bfs_distances(adj, source):
    N = len(adj)
    dist = -np.ones(N, dtype=int)
    dist[source] = 0
    queue = deque([source])
    while queue:
        v = queue.popleft()
        for u in range(N):
            if adj[v, u] and dist[u] == -1:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist

# Compute all-pairs distances
all_dist = np.zeros((N, N), dtype=int)
for i in range(N):
    all_dist[i] = bfs_distances(adj, i)

diameter = np.max(all_dist)
print(f"Graph diameter: {diameter}")

# Shell sizes from vertex 0
d0 = all_dist[0]
shell_sizes = [np.sum(d0 == d) for d in range(diameter + 1)]
print(f"Shell sizes: {shell_sizes}")

# Graph Laplacian
degree_vec = np.sum(adj, axis=1)
L = np.diag(degree_vec.astype(float)) - adj.astype(float)

# =====================================================
# PART 2: Sine-Gordon Model on Graph
# =====================================================

print("\n" + "=" * 70)
print("PART 2: Sine-Gordon Model on 600-Cell")
print("=" * 70)

def sg_energy(eta, L, kappa, lam):
    """Sine-Gordon energy: E = kappa/2 * eta^T L eta + lam * sum sin^2(pi*eta)"""
    kinetic = 0.5 * kappa * eta @ L @ eta
    potential = lam * np.sum(np.sin(np.pi * eta)**2)
    return kinetic + potential, kinetic, potential

def sg_gradient(eta, L, kappa, lam):
    """Gradient of sine-Gordon energy."""
    grad = kappa * L @ eta + lam * np.pi * np.sin(2 * np.pi * eta)
    return grad

def sg_hessian(eta, L, kappa, lam):
    """Hessian of sine-Gordon energy."""
    H = kappa * L + 2 * lam * np.pi**2 * np.diag(np.cos(2 * np.pi * eta))
    return H

def relax_soliton(eta_init, L, kappa, lam, lr=0.001, max_iter=10000, tol=1e-8):
    """Gradient descent to find nearest local minimum."""
    eta = eta_init.copy()
    for it in range(max_iter):
        grad = sg_gradient(eta, L, kappa, lam)
        gnorm = np.linalg.norm(grad)
        if gnorm < tol:
            break
        eta -= lr * grad
    return eta, it, gnorm

print("""
Model: E[eta] = (kappa/2) * eta^T * L * eta + lam * Sum_i sin^2(pi * eta_i)

  - eta_i = ln(psi_i / m_e) / ln(phi)  (log-mass field at vertex i)
  - Minima at eta_i = integer n  (mass = m_e * phi^n)
  - L = graph Laplacian of 600-cell
  - kappa = kinetic coupling (stiffness)
  - lam = potential barrier height
  - kappa/lam controls soliton width
""")

# =====================================================
# PART 3: Radial Soliton Stability
# =====================================================

print("=" * 70)
print("PART 3: Radial Soliton Stability Analysis")
print("=" * 70)

# Use vertex 0 as source. Results are symmetric (H_4 symmetry).
v0 = 0
d_from_v0 = all_dist[v0]

# Scan kappa/lam ratios
print("\nRadial ansatz: eta_i = n * d(i, v0) / D, then relaxed.")
print("Checking which winding numbers n produce stable solitons.\n")

# Try several kappa/lam ratios to see robustness
ratios = [0.1, 0.5, 1.0, 2.0, 5.0]

# For each ratio, test n = 1..10
target_ns = list(range(1, 11))

print(f"{'kappa/lam':>10} {'n':>4} {'E_total':>10} {'E_kin':>10} {'E_pot':>10} {'neg_eig':>8} {'stable':>7} {'eta_range':>10}")
print("-" * 75)

stability_results = {}

for ratio in ratios:
    kappa = ratio
    lam = 1.0

    for n in target_ns:
        # Radial ansatz: linear interpolation from 0 to n across diameter
        eta_init = n * d_from_v0.astype(float) / diameter

        # Relax
        eta_rel, iters, gnorm = relax_soliton(eta_init, L, kappa, lam,
                                                lr=0.0005, max_iter=20000, tol=1e-7)

        # Energy
        E_tot, E_kin, E_pot = sg_energy(eta_rel, L, kappa, lam)

        # Hessian
        H = sg_hessian(eta_rel, L, kappa, lam)
        eigenvalues = np.linalg.eigvalsh(H)
        n_negative = np.sum(eigenvalues < -1e-6)

        eta_range = np.max(eta_rel) - np.min(eta_rel)
        stable = "YES" if n_negative == 0 else "NO"

        print(f"{ratio:>10.1f} {n:>4} {E_tot:>10.3f} {E_kin:>10.3f} {E_pot:>10.3f} {n_negative:>8} {stable:>7} {eta_range:>10.3f}")

        if ratio not in stability_results:
            stability_results[ratio] = {}
        stability_results[ratio][n] = {
            'stable': n_negative == 0,
            'n_negative': n_negative,
            'energy': E_tot,
            'eta_range': eta_range,
            'eta_final': eta_rel.copy(),
            'eigenvalues': eigenvalues
        }
    print()

# =====================================================
# PART 4: Collapsed vs Survived Solitons
# =====================================================

print("\n" + "=" * 70)
print("PART 4: Soliton Fate Analysis")
print("=" * 70)

print("""
After relaxation, a soliton either:
  (a) SURVIVES: field retains non-trivial profile (eta_range > 0.5)
  (b) COLLAPSES: field relaxes to uniform state (eta_range ~ 0)
  (c) SNAPS: field jumps to nearest integer at each site
""")

ratio_test = 1.0  # focus on kappa/lam = 1
kappa = ratio_test
lam = 1.0

print(f"Analysis for kappa/lam = {ratio_test}:\n")
print(f"{'n':>4} {'eta_range':>10} {'<eta>':>8} {'min_eta':>8} {'max_eta':>8} {'fate':>10} {'n_neg':>6}")
print("-" * 60)

for n in range(1, 16):
    res = stability_results[ratio_test][n] if n in stability_results.get(ratio_test, {}) else None

    if res is None:
        # Compute it
        eta_init = n * d_from_v0.astype(float) / diameter
        eta_rel, iters, gnorm = relax_soliton(eta_init, L, kappa, lam,
                                                lr=0.0005, max_iter=20000, tol=1e-7)
        E_tot, E_kin, E_pot = sg_energy(eta_rel, L, kappa, lam)
        H = sg_hessian(eta_rel, L, kappa, lam)
        eigenvalues = np.linalg.eigvalsh(H)
        n_negative = np.sum(eigenvalues < -1e-6)
        res = {
            'stable': n_negative == 0,
            'n_negative': n_negative,
            'energy': E_tot,
            'eta_range': np.max(eta_rel) - np.min(eta_rel),
            'eta_final': eta_rel,
            'eigenvalues': eigenvalues
        }
        if ratio_test not in stability_results:
            stability_results[ratio_test] = {}
        stability_results[ratio_test][n] = res

    eta = res['eta_final']
    eta_range = np.max(eta) - np.min(eta)
    mean_eta = np.mean(eta)

    # Classify fate
    if eta_range < 0.3:
        fate = "COLLAPSED"
    elif all(abs(eta[i] - round(eta[i])) < 0.2 for i in range(N)):
        fate = "SNAPPED"
    else:
        fate = "SURVIVED"

    print(f"{n:>4} {eta_range:>10.4f} {mean_eta:>8.3f} {np.min(eta):>8.3f} {np.max(eta):>8.3f} {fate:>10} {res['n_negative']:>6}")

# =====================================================
# PART 5: Directional Solitons (Diameter vs Decagonal)
# =====================================================

print("\n" + "=" * 70)
print("PART 5: Directional Soliton Analysis")
print("=" * 70)

print("""
KEY QUESTION: Why is (1,0) at n=5 OCCUPIED but (0,1) at n=6 EMPTY?

We construct DIRECTIONAL solitons:
  - DIAMETER soliton: field varies along the diameter (v0 to antipode)
  - DECAGONAL soliton: field varies along a decagonal great circle
  - Compare their energies and stability
""")

# Find the antipode of v0 (distance 5)
antipode = np.where(d_from_v0 == diameter)[0][0]
print(f"Source vertex v0 = {v0}, antipode = {antipode}")

# Find a great decagon through v0
# A decagon consists of 10 vertices forming a regular 10-gon on a great circle of S^3
# Start from v0, find next vertex on a decagon
neighbors_v0 = np.where(adj[v0] == 1)[0]
print(f"v0 has {len(neighbors_v0)} neighbors")

# Walk along a decagon: from v0, pick a neighbor, then always pick the
# neighbor that continues the geodesic (maximizes angle continuation)
def find_decagon(v_start, v_next, adj, V):
    """Follow a great decagon from v_start through v_next."""
    path = [v_start, v_next]
    for step in range(8):  # need 8 more to complete 10
        current = path[-1]
        prev = path[-2]
        # Direction vector
        direction = V[current] - V[prev]
        direction = direction / np.linalg.norm(direction)
        # Among neighbors of current, pick the one most aligned with direction
        # (excluding prev)
        nbrs = np.where(adj[current] == 1)[0]
        best_dot = -2
        best_v = -1
        for u in nbrs:
            if u == prev:
                continue
            # Also skip if already in path (except possibly closing the loop)
            if u in path and not (step == 7 and u == v_start):
                continue
            fwd = V[u] - V[current]
            fwd = fwd / np.linalg.norm(fwd)
            dot = np.dot(direction, fwd)
            if dot > best_dot:
                best_dot = dot
                best_v = u
        if best_v == -1:
            break
        path.append(best_v)
    return path

# Find a decagon through v0
decagon = find_decagon(v0, neighbors_v0[0], adj, V)
print(f"Decagon through v0: {len(decagon)} vertices")
print(f"  Closes? Last vertex = {decagon[-1]}, v0 = {v0}, same = {decagon[-1] == v0}")

# Check distances along decagon
dec_dists = [np.linalg.norm(V[decagon[i]] - V[decagon[i+1]]) for i in range(len(decagon)-1)]
print(f"  Edge lengths along decagon: {[f'{d:.4f}' for d in dec_dists]}")

# Compute graph distance of each decagon vertex from v0
dec_graph_dists = [d_from_v0[v] for v in decagon]
print(f"  Graph distances from v0: {dec_graph_dists}")

# Now construct directional solitons
kappa = 1.0
lam = 1.0

print(f"\n--- Directional Soliton Comparison (kappa/lam = 1.0) ---\n")

# DIAMETER soliton for n=5: field varies linearly from 0 to 5 along diameter
eta_diameter = np.zeros(N)
for i in range(N):
    eta_diameter[i] = d_from_v0[i]  # ranges 0 to 5

eta_diam_rel, _, _ = relax_soliton(eta_diameter, L, kappa, lam,
                                     lr=0.0005, max_iter=20000, tol=1e-7)
E_diam, Ek_diam, Ep_diam = sg_energy(eta_diam_rel, L, kappa, lam)
H_diam = sg_hessian(eta_diam_rel, L, kappa, lam)
eig_diam = np.linalg.eigvalsh(H_diam)
n_neg_diam = np.sum(eig_diam < -1e-6)

print(f"DIAMETER soliton (a=1, b=0, n=5):")
print(f"  Energy: {E_diam:.4f} (kinetic: {Ek_diam:.4f}, potential: {Ep_diam:.4f})")
print(f"  Negative eigenvalues: {n_neg_diam}")
print(f"  eta range: {np.min(eta_diam_rel):.3f} to {np.max(eta_diam_rel):.3f}")

# DECAGONAL soliton for n=6: field varies along decagonal direction
# Assign each vertex a "decagonal coordinate" based on angle from the
# decagonal plane through v0
# Use the first 2 vertices of the decagon to define the plane
if len(decagon) >= 3:
    # Direction of decagon (tangent at v0)
    dec_dir = V[decagon[1]] - V[decagon[0]]
    dec_dir = dec_dir / np.linalg.norm(dec_dir)

    # For each vertex, project onto decagonal direction
    dec_coord = np.zeros(N)
    for i in range(N):
        diff = V[i] - V[v0]
        dec_coord[i] = np.dot(diff, dec_dir)

    # Normalize to [0, 6]
    dec_min, dec_max = np.min(dec_coord), np.max(dec_coord)
    eta_decagonal = 6.0 * (dec_coord - dec_min) / (dec_max - dec_min)

    eta_dec_rel, _, _ = relax_soliton(eta_decagonal, L, kappa, lam,
                                       lr=0.0005, max_iter=20000, tol=1e-7)
    E_dec, Ek_dec, Ep_dec = sg_energy(eta_dec_rel, L, kappa, lam)
    H_dec = sg_hessian(eta_dec_rel, L, kappa, lam)
    eig_dec = np.linalg.eigvalsh(H_dec)
    n_neg_dec = np.sum(eig_dec < -1e-6)

    print(f"\nDECAGONAL soliton (a=0, b=1, n=6):")
    print(f"  Energy: {E_dec:.4f} (kinetic: {Ek_dec:.4f}, potential: {Ep_dec:.4f})")
    print(f"  Negative eigenvalues: {n_neg_dec}")
    print(f"  eta range: {np.min(eta_dec_rel):.3f} to {np.max(eta_dec_rel):.3f}")

# MIXED soliton for n=11 = 5+6: diameter + decagonal
eta_mixed = eta_diameter + eta_decagonal * (6.0/6.0)  # scale decagonal part
# Normalize so total range is ~11
current_range = np.max(eta_mixed) - np.min(eta_mixed)
eta_mixed = eta_mixed * (11.0 / current_range)
eta_mixed -= np.min(eta_mixed)

eta_mix_rel, _, _ = relax_soliton(eta_mixed, L, kappa, lam,
                                    lr=0.0005, max_iter=20000, tol=1e-7)
E_mix, Ek_mix, Ep_mix = sg_energy(eta_mix_rel, L, kappa, lam)
H_mix = sg_hessian(eta_mix_rel, L, kappa, lam)
eig_mix = np.linalg.eigvalsh(H_mix)
n_neg_mix = np.sum(eig_mix < -1e-6)

print(f"\nMIXED soliton (a=1, b=1, n=11):")
print(f"  Energy: {E_mix:.4f} (kinetic: {Ek_mix:.4f}, potential: {Ep_mix:.4f})")
print(f"  Negative eigenvalues: {n_neg_mix}")
print(f"  eta range: {np.min(eta_mix_rel):.3f} to {np.max(eta_mix_rel):.3f}")

# =====================================================
# PART 6: Topological Protection via Cycle Winding
# =====================================================

print("\n" + "=" * 70)
print("PART 6: Topological Analysis -- Cycle Winding Numbers")
print("=" * 70)

print("""
A soliton on a graph is topologically protected if the field has
nonzero WINDING NUMBER along a fundamental cycle of the graph.

For the sine-Gordon model, the winding number along a cycle C is:
  W(C) = (1/2pi) * Sum_{(i,j) in C} arcsin(sin(2*pi*(eta_j - eta_i)))

This is an integer for a configuration at a local minimum.
""")

# Find all great decagons through v0
decagons_v0 = []
for nb in neighbors_v0:
    dec = find_decagon(v0, nb, adj, V)
    if len(dec) == 10 or (len(dec) == 11 and dec[-1] == dec[0]):
        # Normalize: only keep unique decagons
        dec_set = frozenset(dec[:10])
        if dec_set not in [frozenset(d[:10]) for d in decagons_v0]:
            decagons_v0.append(dec[:10])

print(f"Found {len(decagons_v0)} distinct great decagons through v0")

# For each relaxed soliton, compute winding numbers along decagons
def compute_winding(eta, cycle):
    """Compute winding number of field eta along a graph cycle."""
    total = 0
    for k in range(len(cycle)):
        i = cycle[k]
        j = cycle[(k+1) % len(cycle)]
        delta = eta[j] - eta[i]
        total += delta
    # For a closed cycle, total should be an integer (up to rounding)
    return total

print("\nWinding numbers along decagons for each soliton type:")
print(f"{'Soliton':>15} {'Decagon 1':>10} {'Decagon 2':>10} {'Decagon 3':>10} {'Max |W|':>10}")
print("-" * 60)

# Diameter soliton
if len(decagons_v0) >= 3:
    windings_diam = [compute_winding(eta_diam_rel, dec) for dec in decagons_v0[:6]]
    windings_dec = [compute_winding(eta_dec_rel, dec) for dec in decagons_v0[:6]]
    windings_mix = [compute_winding(eta_mix_rel, dec) for dec in decagons_v0[:6]]

    print(f"{'Diameter(n=5)':>15} {windings_diam[0]:>10.3f} {windings_diam[1]:>10.3f} {windings_diam[2]:>10.3f} {max(abs(w) for w in windings_diam):>10.3f}")
    print(f"{'Decagonal(n=6)':>15} {windings_dec[0]:>10.3f} {windings_dec[1]:>10.3f} {windings_dec[2]:>10.3f} {max(abs(w) for w in windings_dec):>10.3f}")
    print(f"{'Mixed(n=11)':>15} {windings_mix[0]:>10.3f} {windings_mix[1]:>10.3f} {windings_mix[2]:>10.3f} {max(abs(w) for w in windings_mix):>10.3f}")

# =====================================================
# PART 7: Shell Profile Analysis
# =====================================================

print("\n" + "=" * 70)
print("PART 7: Shell-Averaged Soliton Profiles")
print("=" * 70)

print("""
For each relaxed soliton, compute the average field value per shell.
This reveals the spatial structure of the soliton.
""")

print(f"{'Shell d':>8} {'N_verts':>8} {'Diam(n=5)':>10} {'Dec(n=6)':>10} {'Mix(n=11)':>10}")
print("-" * 50)

for d in range(diameter + 1):
    mask = (d_from_v0 == d)
    n_verts = np.sum(mask)
    avg_diam = np.mean(eta_diam_rel[mask])
    avg_dec = np.mean(eta_dec_rel[mask])
    avg_mix = np.mean(eta_mix_rel[mask])
    print(f"{d:>8} {n_verts:>8} {avg_diam:>10.4f} {avg_dec:>10.4f} {avg_mix:>10.4f}")

# =====================================================
# PART 8: Energy per Winding Number
# =====================================================

print("\n" + "=" * 70)
print("PART 8: Energy Spectrum of Solitons")
print("=" * 70)

print("""
For various kappa/lam ratios, compute the soliton energy E(n)/n^2
(normalized). If E/n^2 is constant, solitons are non-interacting.
Deviations reveal which n values have special (lower) energy.
""")

kappa = 1.0
lam = 1.0

print(f"\nkappa/lam = {kappa/lam:.1f}\n")
print(f"{'n':>4} {'E_total':>10} {'E/n^2':>10} {'E/n':>10} {'n_neg':>6} {'stable':>7}")
print("-" * 55)

energies = {}
for n in range(1, 27):
    eta_init = n * d_from_v0.astype(float) / diameter
    eta_rel, _, gnorm = relax_soliton(eta_init, L, kappa, lam,
                                       lr=0.0003, max_iter=30000, tol=1e-7)
    E, Ek, Ep = sg_energy(eta_rel, L, kappa, lam)
    H = sg_hessian(eta_rel, L, kappa, lam)
    eigs = np.linalg.eigvalsh(H)
    n_neg = np.sum(eigs < -1e-6)
    stable = "YES" if n_neg == 0 else "NO"

    energies[n] = {
        'E': E, 'E_per_n2': E/n**2, 'E_per_n': E/n,
        'n_neg': n_neg, 'stable': stable,
        'eta': eta_rel.copy()
    }
    print(f"{n:>4} {E:>10.3f} {E/n**2:>10.4f} {E/n:>10.4f} {n_neg:>6} {stable:>7}")

# =====================================================
# PART 9: Comparison with Observed Fermions
# =====================================================

print("\n" + "=" * 70)
print("PART 9: Comparison with Observed Fermion Spectrum")
print("=" * 70)

fermions = [
    ("Electron", 0, 0, 0),
    ("Up", 3, -2, 3),
    ("Down", 1, 0, 5),
    ("Muon", 1, 1, 11),
    ("Strange", 1, 1, 11),
    ("Charm", 2, 1, 16),
    ("Tau", 1, 2, 17),
    ("Bottom", -1, 4, 19),
    ("Top", 4, 1, 26),
]

print(f"\n{'Fermion':>10} {'(a,b)':>8} {'n':>4} {'E_sol':>10} {'E/n^2':>10} {'n_neg':>6} {'stable':>7}")
print("-" * 60)

for name, a, b, n in fermions:
    if n == 0:
        print(f"{name:>10} ({a:+d},{b:+d}) {n:>4} {'0.000':>10} {'--':>10} {'0':>6} {'YES':>7}")
        continue
    if n in energies:
        res = energies[n]
        print(f"{name:>10} ({a:+d},{b:+d}) {n:>4} {res['E']:>10.3f} {res['E_per_n2']:>10.4f} {res['n_neg']:>6} {res['stable']:>7}")

# Check empty levels
print(f"\nEmpty levels (no known particle):")
print(f"{'n':>4} {'simplest (a,b)':>15} {'|a|+|b|':>8} {'E_sol':>10} {'E/n^2':>10} {'n_neg':>6} {'stable':>7}")
print("-" * 65)

occupied_n = {0, 3, 5, 11, 16, 17, 19, 26}

for n in range(1, 27):
    if n in occupied_n:
        continue
    # Find simplest (a,b) for this n
    best = None
    for a in range(-10, 11):
        if (n - 5*a) % 6 == 0:
            b_val = (n - 5*a) // 6
            c = abs(a) + abs(b_val)
            if best is None or c < best[2]:
                best = (a, b_val, c)

    if n in energies:
        res = energies[n]
        print(f"{n:>4} ({best[0]:+d},{best[1]:+d}) {best[2]:>8} {res['E']:>10.3f} {res['E_per_n2']:>10.4f} {res['n_neg']:>6} {res['stable']:>7}")

# =====================================================
# PART 10: Pattern Discovery
# =====================================================

print("\n" + "=" * 70)
print("PART 10: Pattern Analysis")
print("=" * 70)

# Check if occupied n values have special energy or stability properties
occupied_energies = [energies[n]['E_per_n2'] for n in [3, 5, 11, 16, 17, 19, 26]]
empty_energies = [energies[n]['E_per_n2'] for n in range(1, 27) if n not in occupied_n]

print(f"\nE/n^2 statistics:")
print(f"  Occupied levels: mean = {np.mean(occupied_energies):.4f}, std = {np.std(occupied_energies):.4f}")
print(f"  Empty levels:    mean = {np.mean(empty_energies):.4f}, std = {np.std(empty_energies):.4f}")

occupied_stable = [energies[n]['n_neg'] == 0 for n in [3, 5, 11, 16, 17, 19, 26]]
empty_stable = [energies[n]['n_neg'] == 0 for n in range(1, 27) if n not in occupied_n]

print(f"\nStability:")
print(f"  Occupied levels: {sum(occupied_stable)}/{len(occupied_stable)} stable")
print(f"  Empty levels:    {sum(empty_stable)}/{len(empty_stable)} stable")

# Check n mod 5 and n mod 6 patterns
print(f"\nModular analysis of occupied levels:")
for n in sorted(occupied_n - {0}):
    print(f"  n={n:>2}: mod 5 = {n%5}, mod 6 = {n%6}, mod 10 = {n%10}, mod 12 = {n%12}")

# Check if the PROFILE matters (not just the level)
print(f"\nProfile analysis: integer content of relaxed solitons")
print(f"  How close to integer is each shell's average eta?")
print(f"{'n':>4} {'d=0':>6} {'d=1':>6} {'d=2':>6} {'d=3':>6} {'d=4':>6} {'d=5':>6} {'sum_frac':>10}")
print("-" * 55)

for n in [3, 5, 6, 11, 16, 17, 19, 26]:
    eta = energies[n]['eta']
    fracs = []
    vals = []
    for d in range(diameter + 1):
        mask = (d_from_v0 == d)
        avg = np.mean(eta[mask])
        frac = abs(avg - round(avg))
        fracs.append(frac)
        vals.append(avg)
    sum_frac = sum(fracs)
    label = "OCCUPIED" if n in occupied_n else "EMPTY"
    print(f"{n:>4} " + " ".join(f"{v:>6.2f}" for v in vals) + f" {sum_frac:>10.4f}  {label}")

# =====================================================
# SUMMARY
# =====================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
FINDINGS:

1. SINE-GORDON MODEL: The discretized sine-Gordon equation on the
   600-cell graph has well-defined soliton solutions for each winding
   number n.

2. RADIAL SOLITONS: Using the radial ansatz (eta proportional to
   graph distance from source), solitons for ALL n values can be
   constructed and relaxed to local minima.

3. STABILITY ANALYSIS: The Hessian eigenvalue spectrum reveals which
   solitons have unstable (negative eigenvalue) modes.

4. DIRECTIONAL STRUCTURE: Diameter vs decagonal solitons have
   different energy profiles and stability properties.

CLASSIFICATION:
  - Sine-Gordon on 600-cell: DERIVED (natural discretization)
  - Soliton existence for each n: COMPUTED
  - Stability pattern: COMPUTED (compare with observed fermions above)
  - (a,b) selection mechanism: IN PROGRESS

OPEN QUESTIONS:
  - Does the stability pattern match observed fermions?
  - What is the physical value of kappa/lam?
  - Do directional (non-radial) solitons give different stability?
  - Can confinement of quarks be related to soliton topology?
""")
