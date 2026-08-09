"""
EXP-165: Topological Solitons from 5-Coset Structure
=====================================================
The 5 inscribed 24-cells give a natural Z_5 coloring of the 600-cell.
This is like a Potts model with q=5 states.

Excitations = "flips" of vertices between cosets.
Solitons = topological defects (domain walls, kinks, vortices).

Key structure:
- 600-cell is complete 5-partite (0 intra-coset edges, 72 between each pair)
- Each fermion vertex has exactly 3 neighbors per OTHER coset (from 12 total)
- Domain wall energy = number of "frustrated" edges

Questions:
1. What are the elementary excitations (single-vertex flips)?
2. Do stable multi-vertex solitons exist?
3. Is there a topological charge (winding number)?
4. Do soliton masses relate to fermion masses?
"""

import numpy as np
from itertools import product, permutations
from collections import Counter, defaultdict

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-165: Topological Solitons from 5-Coset Structure")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell
# ============================================================
print("\n--- Step 1: Build 600-cell ---")

verts_set = set()
for i in range(4):
    for s in [1, -1]:
        v = [0,0,0,0]; v[i] = s
        verts_set.add(tuple(v))
for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = [p for p in permutations(range(4))
              if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]
for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in product([1,-1], repeat=len(nz)):
        v = list(base)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

verts = np.array(sorted(verts_set))
N = len(verts)
dots = verts @ verts.T
adj = (np.abs(dots - PHI/2) < 0.01).astype(int)
np.fill_diagonal(adj, 0)
print(f"  N={N}, edges={adj.sum()//2}")

# ============================================================
# Step 2: Find 5 inscribed 24-cells
# ============================================================
print("\n--- Step 2: 5 inscribed 24-cells ---")

# Quaternion multiplication
def qmul(q1, q2):
    a,b,c,d = q1; e,f,g,h = q2
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

# Type A + B = coset 0 (binary tetrahedral group 2T)
coset0_idx = []
for i in range(N):
    c = verts[i]
    nz = np.sum(np.abs(c) > 0.01)
    if nz == 1 and np.isclose(np.max(np.abs(c)), 1.0):
        coset0_idx.append(i)  # Type A
    elif nz == 4 and np.allclose(np.abs(c), 0.5, atol=0.01):
        coset0_idx.append(i)  # Type B

cosets = [set(coset0_idx)]
used = set(coset0_idx)

# Find other cosets by quaternion multiplication
type_C = [i for i in range(N) if i not in used]
for rep_idx in type_C:
    if rep_idx in used:
        continue
    rep = verts[rep_idx]
    new_coset = set()
    for v_idx in coset0_idx:
        prod = qmul(rep, verts[v_idx])
        dists_sq = np.sum((verts - prod)**2, axis=1)
        closest = np.argmin(dists_sq)
        if dists_sq[closest] < 0.001:
            new_coset.add(closest)
    if len(new_coset) == 24 and not new_coset & used:
        cosets.append(new_coset)
        used |= new_coset
    if len(cosets) == 5:
        break

print(f"  Found {len(cosets)} cosets")
for ci, c in enumerate(cosets):
    print(f"    Coset {ci}: {len(c)} vertices")

# Assign coset to each vertex
vertex_coset = np.zeros(N, dtype=int)
for ci, c in enumerate(cosets):
    for v in c:
        vertex_coset[v] = ci

# Verify: 0 intra-coset edges
intra = 0
for i in range(N):
    for j in range(i+1, N):
        if adj[i,j] and vertex_coset[i] == vertex_coset[j]:
            intra += 1
print(f"  Intra-coset edges: {intra} (should be 0)")

# Cross-coset neighbor count
print(f"\n  Cross-coset neighbor pattern:")
for ci in range(5):
    ci_verts = [v for v in range(N) if vertex_coset[v] == ci]
    nbr_counts = Counter()
    for v in ci_verts:
        for n in range(N):
            if adj[v, n]:
                nbr_counts[vertex_coset[n]] += 1
    # Average per vertex
    for cj in range(5):
        if ci != cj:
            avg = nbr_counts[cj] / len(ci_verts)
            if ci == 0:
                print(f"    Coset 0 -> Coset {cj}: {avg:.1f} neighbors/vertex")

# ============================================================
# Step 3: Potts model Hamiltonian
# ============================================================
print("\n--- Step 3: Potts model ---")

# H = -J * sum_{<i,j>} delta(sigma_i, sigma_j)
# Ground state: all sigma_i = vertex_coset[i] (natural coloring)
# Energy = -J * (number of same-coset neighbor pairs) = -J * 0 = 0
# (because 0 intra-coset edges!)

# Actually, the ANTI-ferromagnetic Potts model is more natural:
# H = +J * sum_{<i,j>} delta(sigma_i, sigma_j)
# Ground state minimizes same-coset neighbor pairs = 0 (the natural coloring IS the ground state!)

# Energy of a configuration sigma:
def potts_energy(sigma, adj_matrix):
    """Energy = number of same-color neighbor pairs."""
    E = 0
    for i in range(len(sigma)):
        for j in range(i+1, len(sigma)):
            if adj_matrix[i,j] and sigma[i] == sigma[j]:
                E += 1
    return E

# Ground state energy
E0 = potts_energy(vertex_coset, adj)
print(f"  Ground state energy: {E0} (should be 0)")

# ============================================================
# Step 4: Elementary excitations (single-vertex flip)
# ============================================================
print("\n--- Step 4: Single-vertex flip excitations ---")

# Flip vertex v from coset c to coset c': energy cost = number of c'-neighbors
# Because flipping to c' creates new same-color pairs with c'-neighbors

flip_energies = defaultdict(list)

for v in range(N):
    c_orig = vertex_coset[v]
    nbrs = np.where(adj[v] > 0)[0]
    nbr_cosets = Counter(vertex_coset[n] for n in nbrs)

    for c_new in range(5):
        if c_new == c_orig:
            continue
        # Energy cost = number of neighbors in coset c_new
        delta_E = nbr_cosets.get(c_new, 0)
        flip_energies[(c_orig, c_new)].append(delta_E)

print(f"  Flip energy statistics (coset_from -> coset_to):")
all_flip_E = []
for (c_from, c_to), energies in sorted(flip_energies.items()):
    mean_E = np.mean(energies)
    all_flip_E.extend(energies)
    if c_from == 0:
        print(f"    {c_from} -> {c_to}: mean={mean_E:.1f}, "
              f"values={Counter(energies).most_common(5)}")

all_flip_E = np.array(all_flip_E)
print(f"\n  Overall flip energies: {Counter(all_flip_E).most_common()}")
print(f"  Mean flip energy: {np.mean(all_flip_E):.2f}")
print(f"  This is the 'mass' of a single quasiparticle excitation")

# ============================================================
# Step 5: Two-vertex excitations (pair flips)
# ============================================================
print("\n--- Step 5: Pair excitations ---")

# Flip two adjacent vertices to each other's coset (swap)
# This creates a "dipole" excitation
swap_energies = []

for i in range(N):
    for j in range(i+1, N):
        if not adj[i,j]:
            continue

        ci, cj = vertex_coset[i], vertex_coset[j]
        if ci == cj:
            continue  # Same coset, no swap possible

        # Energy before swap: 0 (ground state)
        # After swap: i gets coset cj, j gets coset ci
        # New same-color pairs:
        # - i's neighbors in coset cj (excluding j)
        # - j's neighbors in coset ci (excluding i)
        # But we lose nothing (ground state had 0 same-color pairs)

        nbrs_i = np.where(adj[i] > 0)[0]
        nbrs_j = np.where(adj[j] > 0)[0]

        delta_E_i = sum(1 for n in nbrs_i if n != j and vertex_coset[n] == cj)
        delta_E_j = sum(1 for n in nbrs_j if n != i and vertex_coset[n] == ci)

        # Also: i and j now have same coset? No - they SWAPPED, so i=cj, j=ci, still different
        total_delta_E = delta_E_i + delta_E_j
        swap_energies.append(total_delta_E)

swap_E = np.array(swap_energies)
print(f"  Pair swap energies: {Counter(swap_E).most_common()}")
print(f"  Mean: {np.mean(swap_E):.2f}, Min: {np.min(swap_E)}, Max: {np.max(swap_E)}")

# ============================================================
# Step 6: Ring excitations (cyclic permutation along a path)
# ============================================================
print("\n--- Step 6: Ring/loop excitations ---")

# Find shortest cycles in the graph
# Triangle: 3 vertices, each in different coset
# Cyclic permutation: a->b->c->a (shift cosets along the triangle)

triangle_energies = []
for i in range(N):
    nbrs_i = set(np.where(adj[i] > 0)[0])
    for j in nbrs_i:
        if j <= i:
            continue
        common = nbrs_i & set(np.where(adj[j] > 0)[0])
        for k in common:
            if k <= j:
                continue
            # Triangle (i, j, k)
            ci, cj, ck = vertex_coset[i], vertex_coset[j], vertex_coset[k]
            if len({ci, cj, ck}) < 3:
                continue  # Need all different cosets for cyclic permutation

            # Cyclic: i->cj, j->ck, k->ci
            sigma_new = vertex_coset.copy()
            sigma_new[i] = cj
            sigma_new[j] = ck
            sigma_new[k] = ci

            # Energy change: count new same-color pairs minus old
            delta_E = 0
            for v in [i, j, k]:
                for n in np.where(adj[v] > 0)[0]:
                    if n in {i, j, k}:
                        continue
                    # Old: different (ground state). New: check
                    if sigma_new[v] == vertex_coset[n]:
                        delta_E += 1

            # Also check within the triangle
            for v1, v2 in [(i,j), (j,k), (i,k)]:
                if sigma_new[v1] == sigma_new[v2]:
                    delta_E += 1

            triangle_energies.append(delta_E)

tri_E = np.array(triangle_energies) if triangle_energies else np.array([])
if len(tri_E) > 0:
    print(f"  Triangle cyclic permutation energies:")
    print(f"    Distribution: {Counter(tri_E).most_common()}")
    print(f"    Mean: {np.mean(tri_E):.2f}, Min: {np.min(tri_E)}")
    print(f"    Zero-energy triangles: {np.sum(tri_E == 0)}/{len(tri_E)}")
else:
    print(f"  No valid triangle permutations found")

# ============================================================
# Step 7: Topological charge
# ============================================================
print("\n--- Step 7: Topological structure ---")

# The coset structure defines a map: 600-cell -> Z_5
# On the boundary of each triangle, we can define a "winding number"
# W(triangle) = (sigma_j - sigma_i + sigma_k - sigma_j + sigma_i - sigma_k) mod 5
# This is trivially 0 for any configuration.

# Better: define a Z_5 gauge field on edges:
# A(i,j) = (sigma_j - sigma_i) mod 5
# Flux through triangle = (A(i,j) + A(j,k) + A(k,i)) mod 5

# For ground state:
flux_counts = Counter()
n_triangles = 0
for i in range(N):
    nbrs_i = set(np.where(adj[i] > 0)[0])
    for j in nbrs_i:
        if j <= i:
            continue
        for k in (nbrs_i & set(np.where(adj[j] > 0)[0])):
            if k <= j:
                continue
            ci, cj, ck = vertex_coset[i], vertex_coset[j], vertex_coset[k]
            A_ij = (cj - ci) % 5
            A_jk = (ck - cj) % 5
            A_ki = (ci - ck) % 5
            flux = (A_ij + A_jk + A_ki) % 5
            flux_counts[flux] += 1
            n_triangles += 1

print(f"  Z_5 flux through {n_triangles} triangles:")
for f, count in sorted(flux_counts.items()):
    print(f"    Flux={f}: {count} ({100*count/n_triangles:.1f}%)")

# ============================================================
# Step 8: Ordered vs disordered phases
# ============================================================
print("\n--- Step 8: Phase structure ---")

# At T=0: ground state = natural coloring (E=0)
# At T=infinity: random coloring
# Q: is there a phase transition?

# Random coloring energy
n_samples = 1000
random_energies = []
for _ in range(n_samples):
    sigma_rand = np.random.randint(0, 5, size=N)
    E = 0
    for i in range(N):
        for j in range(i+1, N):
            if adj[i,j] and sigma_rand[i] == sigma_rand[j]:
                E += 1
    random_energies.append(E)

random_E = np.array(random_energies)
print(f"  Random coloring energy: mean={np.mean(random_E):.1f}, std={np.std(random_E):.1f}")
print(f"  Expected: 720 * 1/5 = {720/5:.0f} (random same-color probability 1/5)")
print(f"  Ground state: {E0}")
print(f"  Energy gap: {np.mean(random_E) - E0:.1f}")

# ============================================================
# Step 9: Single-flip spectrum (mass spectrum)
# ============================================================
print("\n--- Step 9: Excitation mass spectrum ---")

# Each single flip has energy = number of neighbors in target coset
# This gives a discrete "mass spectrum"

# For each vertex, compute its flip energy to each other coset
mass_spectrum = defaultdict(list)
for v in range(N):
    nbrs = np.where(adj[v] > 0)[0]
    nbr_cosets = Counter(vertex_coset[n] for n in nbrs)
    c_orig = vertex_coset[v]

    for c_new in range(5):
        if c_new == c_orig:
            continue
        E = nbr_cosets.get(c_new, 0)
        mass_spectrum[E].append((v, c_orig, c_new))

print(f"  Mass spectrum of single-flip excitations:")
total_flips = sum(len(v) for v in mass_spectrum.values())
for E in sorted(mass_spectrum.keys()):
    count = len(mass_spectrum[E])
    print(f"    E={E}: {count} excitations ({100*count/total_flips:.1f}%)")

# Check if mass values relate to phi
print(f"\n  Mass values as phi-expressions:")
for E in sorted(mass_spectrum.keys()):
    for a in range(-5, 6):
        for b in range(-5, 6):
            if abs(a + b * PHI - E) < 0.01:
                print(f"    E={E} = {a} + {b}*phi")

# ============================================================
# Step 10: Stability analysis
# ============================================================
print("\n--- Step 10: Stability of excitations ---")

# A single-flip excitation with energy E: is it stable or can it lower energy?
# Check: can the flipped vertex move to a LOWER energy coset?
# Also: does the excitation "attract" further flips (instability)?

# For a vertex flipped from c0 to c1:
# Its neighbors now have a same-color neighbor (if they're in c1)
# These neighbors can lower their energy by flipping AWAY from c1
# This is a "cascade" effect

for E_val in sorted(mass_spectrum.keys()):
    example = mass_spectrum[E_val][0]
    v, c_from, c_to = example

    # After flipping v from c_from to c_to:
    # Check each neighbor in c_to: can they lower energy by flipping?
    nbrs = np.where(adj[v] > 0)[0]
    nbrs_in_target = [n for n in nbrs if vertex_coset[n] == c_to]

    cascade_possible = 0
    for n in nbrs_in_target:
        # n is now "frustrated" (has same-color neighbor v)
        # Can n flip to any other coset with lower energy?
        n_nbrs = np.where(adj[n] > 0)[0]
        n_nbr_cosets = Counter(vertex_coset[nn] for nn in n_nbrs if nn != v)
        n_nbr_cosets[c_to] += 1  # v is now in c_to

        current_E = n_nbr_cosets.get(c_to, 0)  # same-color neighbors after v's flip
        for c_alt in range(5):
            if c_alt == c_to:
                continue
            alt_E = n_nbr_cosets.get(c_alt, 0)
            if alt_E < current_E:
                cascade_possible += 1
                break

    print(f"  E={E_val}: {len(nbrs_in_target)} frustrated neighbors, "
          f"{cascade_possible} can cascade")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key findings:
1. Antiferromagnetic q=5 Potts ground state = natural coset coloring (E=0)
2. Single-flip excitations have discrete mass spectrum
3. Pair swaps and triangle permutations create higher excitations
4. Z_5 gauge flux through triangles characterizes topological sectors
5. Cascade stability determines if excitations are bound or unbound

STATUS: Check output above
""")
