"""
EXP-166: Soliton Bound States on 600-Cell
==========================================
From exp165: E_flip=3 for single, E_swap=4 for pair.
Binding formula E=3n-2e for swaps (each internal swap-edge saves 2).

Questions:
1. What are the densest subgraphs? (max edges e for given vertices n)
2. How does the energy depend on the excitation type (swap vs cyclic)?
3. Are there zero-energy or near-zero solitons?
4. Which bound states could correspond to hadrons?
5. Connection E_single=3 to M_vacancy=3 and dim(S^3)?
"""

import numpy as np
from itertools import product, permutations, combinations
from collections import Counter, defaultdict

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-166: Soliton Bound States")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell + coset structure
# ============================================================

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

# Quaternion mult
def qmul(q1, q2):
    a,b,c,d = q1; e,f,g,h = q2
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

# Build 5 cosets
coset0 = []
for i in range(N):
    c = verts[i]
    nz = np.sum(np.abs(c) > 0.01)
    if (nz == 1 and np.isclose(np.max(np.abs(c)), 1.0)) or \
       (nz == 4 and np.allclose(np.abs(c), 0.5, atol=0.01)):
        coset0.append(i)

cosets = [set(coset0)]
used = set(coset0)
type_C = [i for i in range(N) if i not in used]
for rep_idx in type_C:
    if rep_idx in used: continue
    new_coset = set()
    for v_idx in coset0:
        prod = qmul(verts[rep_idx], verts[v_idx])
        closest = np.argmin(np.sum((verts - prod)**2, axis=1))
        if np.sum((verts[closest] - prod)**2) < 0.001:
            new_coset.add(closest)
    if len(new_coset) == 24 and not new_coset & used:
        cosets.append(new_coset)
        used |= new_coset
    if len(cosets) == 5: break

vertex_coset = np.zeros(N, dtype=int)
for ci, c in enumerate(cosets):
    for v in c: vertex_coset[v] = ci

print(f"  N={N}, 5 cosets of 24, edges={adj.sum()//2}")

# ============================================================
# Step 2: General excitation energy calculator
# ============================================================

def excitation_energy(flips, vertex_coset_arr, adj_matrix):
    """
    flips: dict {vertex_index: new_coset}
    Returns number of frustrated (same-color neighbor) edges.
    """
    sigma = vertex_coset_arr.copy()
    for v, c in flips.items():
        sigma[v] = c

    E = 0
    for v in flips:
        for n in range(len(sigma)):
            if n > v and adj_matrix[v, n] and sigma[v] == sigma[n]:
                E += 1
            elif n < v and n not in flips and adj_matrix[v, n] and sigma[v] == sigma[n]:
                E += 1
    # Also count edges between non-flipped neighbors that got frustrated
    # Actually, only flipped vertices create frustration (ground state has E=0)
    # So only edges involving at least one flipped vertex matter

    # Recount properly
    E = 0
    for v in flips:
        for n in np.where(adj_matrix[v] > 0)[0]:
            if sigma[v] == sigma[n]:
                E += 1
    # Each edge counted once per flipped endpoint. Edges between two flipped: counted twice.
    # Need to subtract double-counted edges between flipped pairs
    flipped_set = set(flips.keys())
    double_counted = 0
    for v1 in flipped_set:
        for v2 in flipped_set:
            if v1 < v2 and adj_matrix[v1, v2] and sigma[v1] == sigma[v2]:
                double_counted += 1

    return (E - double_counted)

# ============================================================
# Step 3: Verify known results
# ============================================================
print("\n--- Step 3: Verify known results ---")

# Single flip
v0 = list(cosets[0])[0]
for target_c in range(1, 5):
    E = excitation_energy({v0: target_c}, vertex_coset, adj)
    print(f"  Single flip {v0}(c0->c{target_c}): E={E}")

# Pair swap
v0 = list(cosets[0])[0]
nbrs = np.where(adj[v0] > 0)[0]
for n in nbrs[:3]:
    cn = vertex_coset[n]
    if cn == 0: continue
    E = excitation_energy({v0: cn, n: 0}, vertex_coset, adj)
    print(f"  Pair swap {v0}(c0)-{n}(c{cn}): E={E}")
    break

# Triangle cyclic
for i in range(N):
    nbrs_i = set(np.where(adj[i] > 0)[0])
    for j in nbrs_i:
        if j <= i: continue
        for k in (nbrs_i & set(np.where(adj[j] > 0)[0])):
            if k <= j: continue
            ci, cj, ck = vertex_coset[i], vertex_coset[j], vertex_coset[k]
            if len({ci, cj, ck}) == 3:
                E = excitation_energy({i: cj, j: ck, k: ci}, vertex_coset, adj)
                print(f"  Triangle cyclic ({i},{j},{k}) cosets ({ci},{cj},{ck}): E={E}")
                break
        else: continue
        break
    else: continue
    break

# ============================================================
# Step 4: Tetrahedron excitations
# ============================================================
print("\n--- Step 4: Tetrahedron excitations ---")

# Find tetrahedra (4-cliques) with vertices in different cosets
tet_count = 0
tet_energies = Counter()
tet_coset_patterns = Counter()

for i in range(N):
    ni = set(np.where(adj[i] > 0)[0])
    for j in ni:
        if j <= i: continue
        nij = ni & set(np.where(adj[j] > 0)[0])
        for k in nij:
            if k <= j: continue
            nijk = nij & set(np.where(adj[k] > 0)[0])
            for l in nijk:
                if l <= k: continue
                coset_set = {vertex_coset[v] for v in [i,j,k,l]}
                n_cosets = len(coset_set)
                tet_coset_patterns[n_cosets] += 1

                if n_cosets >= 3 and tet_count < 20:
                    ci, cj, ck, cl = [vertex_coset[v] for v in [i,j,k,l]]

                    # Try all permutations of target cosets
                    best_E = 999
                    best_perm = None
                    coset_list = [ci, cj, ck, cl]
                    for perm in permutations(range(len(coset_list))):
                        targets = [coset_list[perm[p]] for p in range(4)]
                        # Each vertex gets the coset of the permuted position
                        if targets == coset_list: continue  # identity
                        flips = {}
                        for idx_v, v in enumerate([i,j,k,l]):
                            if targets[idx_v] != coset_list[idx_v]:
                                flips[v] = targets[idx_v]
                        if flips:
                            E = excitation_energy(flips, vertex_coset, adj)
                            if E < best_E:
                                best_E = E
                                best_perm = perm

                    if best_E < 999:
                        tet_energies[best_E] += 1
                        if tet_count < 5:
                            print(f"  Tet ({i},{j},{k},{l}) cosets=({ci},{cj},{ck},{cl}): best E={best_E}")

                tet_count += 1

print(f"\n  Total tetrahedra: {tet_count}")
print(f"  Coset patterns: {dict(sorted(tet_coset_patterns.items()))}")
print(f"  Best energies (first 20 tets): {dict(sorted(tet_energies.items()))}")

# ============================================================
# Step 5: Systematic dense subgraph search
# ============================================================
print("\n--- Step 5: Dense subgraph solitons ---")

# For n=2,3,4,5: find subgraphs with maximum edges
# and compute their minimum excitation energy

# Build neighbor lists
nbr_list = [set(np.where(adj[i] > 0)[0]) for i in range(N)]

# n=2: all edges (e=1)
print(f"  n=2: e=1 (edge), E=4 (all identical)")

# n=3: triangles (e=3)
print(f"  n=3: e=3 (triangle), E=6 cyclic, check swap...")

# For a triangle, try all possible flip patterns
sample_tri = None
for i in range(N):
    for j in nbr_list[i]:
        if j <= i: continue
        for k in (nbr_list[i] & nbr_list[j]):
            if k <= j: continue
            sample_tri = (i, j, k)
            break
        if sample_tri: break
    if sample_tri: break

if sample_tri:
    i, j, k = sample_tri
    ci, cj, ck = vertex_coset[i], vertex_coset[j], vertex_coset[k]
    print(f"\n  Sample triangle ({i},{j},{k}), cosets ({ci},{cj},{ck}):")

    # Try all possible target assignments
    min_E = 999
    for ti in range(5):
        for tj in range(5):
            for tk in range(5):
                flips = {}
                if ti != ci: flips[i] = ti
                if tj != cj: flips[j] = tj
                if tk != ck: flips[k] = tk
                if not flips: continue
                E = excitation_energy(flips, vertex_coset, adj)
                if E < min_E:
                    min_E = E
                    print(f"    Targets ({ti},{tj},{tk}): E={E} [n_flips={len(flips)}]")
                    min_E = E
    print(f"  Min energy triangle excitation: {min_E}")

# n=4: look for 4-cliques and compute minimum energy
print(f"\n  n=4 (tetrahedra):")
# Already computed above

# n=5: pentagons or other 5-vertex subgraphs
print(f"\n  n=5: searching for dense 5-vertex subgraphs...")

# Find 5-vertex subgraphs with many edges
best_5_edges = 0
best_5_verts = None
# Sample: take a vertex, its neighbor, and their common neighbors
v0 = 0
for v1 in nbr_list[v0]:
    common = nbr_list[v0] & nbr_list[v1]
    for v2 in common:
        if v2 <= v1: continue
        for v3 in (common & nbr_list[v2]):
            if v3 <= v2: continue
            for v4 in (common & nbr_list[v2] & nbr_list[v3]):
                if v4 <= v3: continue
                subset = [v0, v1, v2, v3, v4]
                e = sum(1 for a in range(5) for b in range(a+1,5) if adj[subset[a], subset[b]])
                if e > best_5_edges:
                    best_5_edges = e
                    best_5_verts = subset

if best_5_verts:
    print(f"  Best 5-vertex subgraph: {best_5_verts}, edges={best_5_edges}")
    subset_cosets = [vertex_coset[v] for v in best_5_verts]
    print(f"    Cosets: {subset_cosets}")

# ============================================================
# Step 6: Binding energy analysis
# ============================================================
print("\n--- Step 6: Binding energy formula ---")

# For a general excitation with n flipped vertices and e internal edges:
# Savings per edge depends on the flip pattern
#
# Case 1: SWAP (i->cj, j->ci): saves 2 per edge (both directions match)
# Case 2: CYCLIC (i->cj, j->ck, k->ci): saves 1 per edge (one direction matches)
# Case 3: RANDOM: saves 0-2 per edge depending on matches

# General formula: E = 3n - sum_edges(savings)
# where savings(e) = (1 if j_orig == i_target) + (1 if i_orig == j_target)

# For uniform 3-nbrs-per-coset: each vertex has exactly 3 external neighbors in any coset
# So the key variable is how many of the internal edges are "saving" edges

print(f"  E = 3n - S, where S = sum of savings per internal edge")
print(f"  Savings per edge: 0 (unmatched), 1 (one-way match), 2 (swap)")
print(f"  Maximum savings: S_max = 2e (all swaps)")
print(f"  Minimum energy: E_min = 3n - 2e")

# Check: what's the maximum e for n vertices in the 600-cell?
# 600-cell is 12-regular, but clique number is limited
# Find clique number
max_clique = 0
for i in range(N):
    for j in nbr_list[i]:
        if j <= i: continue
        common_ij = nbr_list[i] & nbr_list[j]
        for k in common_ij:
            if k <= j: continue
            common_ijk = common_ij & nbr_list[k]
            for l in common_ijk:
                if l <= k: continue
                common_ijkl = common_ijk & nbr_list[l]
                if len(common_ijkl) > 0:
                    max_clique = max(max_clique, 5)
                else:
                    max_clique = max(max_clique, 4)

print(f"\n  Max clique size in 600-cell: {max_clique}")

# For each n, what's the theoretical minimum E?
print(f"\n  Theoretical minimum E (assuming all-swap is possible):")
for n in range(1, 7):
    # Max edges in a subgraph of 600-cell with n vertices
    if n == 1: e_max = 0
    elif n == 2: e_max = 1
    elif n == 3: e_max = 3  # triangle
    elif n == 4: e_max = 6  # tetrahedron (K4)
    elif n == 5: e_max = best_5_edges if best_5_verts else 9
    elif n == 6: e_max = 12  # estimate

    E_swap = 3*n - 2*e_max
    E_cyclic = 3*n - e_max  # if only one-way matches
    print(f"    n={n}: e_max={e_max}, E_swap={max(0,E_swap)}, E_cyclic={max(0,E_cyclic)}")

# ============================================================
# Step 7: Can we achieve E=0?
# ============================================================
print("\n--- Step 7: Zero-energy solitons? ---")

# E=0 requires 3n = 2e (for all-swap), i.e. e = 3n/2
# n=2: e=3 -> impossible (only 1 edge between 2 vertices)
# n=4: e=6 -> K4 (tetrahedron), possible!
# n=6: e=9 -> need 9 edges among 6 vertices

# But for tetrahedron: can we do all-swap? Need to pair up 4 vertices into 2 swapping pairs
# With vertices in cosets c0,c1,c2,c3: swap (c0,c1) and (c2,c3)
# The saving per edge:
# - (v0,v1): swap, saves 2
# - (v2,v3): swap, saves 2
# - (v0,v2): v0->c1, v2->c3. v0 was c0≠c3, v2 was c2≠c1. Saves 0.
# Total savings: 4, not 12. E = 12-4 = 8.

# For ALL edges to save 2: every edge must be a swap edge
# This means: for every edge (i,j), i's target = j's original AND j's target = i's original
# This is ONLY possible for perfect matchings of the flip graph

# Alternatively: can we find a vertex set where a cyclic permutation gives E=0?
# E_cyclic = 3n - e (one-way saves). E=0 requires e=3n.
# But max edges in n-vertex graph = n(n-1)/2. Need n(n-1)/2 >= 3n, so n >= 7.
# And all edges must be one-way saves. Very restrictive.

# Conclusion: E=0 soliton likely impossible for small n
print(f"  For all-swap: E=0 needs e=3n/2. For K4: 3*4/2=6=K4 edges.")
print(f"  But K4 swap only saves 4 (not 12), because cross-pair edges save 0.")

# Actually verify with the calculator
print(f"\n  Verifying K4 all-swap configurations:")
# Find a tetrahedron with 4 different cosets
for i in range(N):
    found = False
    for j in nbr_list[i]:
        if j <= i: continue
        for k in (nbr_list[i] & nbr_list[j]):
            if k <= j: continue
            for l in (nbr_list[i] & nbr_list[j] & nbr_list[k]):
                if l <= k: continue
                cs = [vertex_coset[v] for v in [i,j,k,l]]
                if len(set(cs)) == 4:
                    print(f"  Tet ({i},{j},{k},{l}), cosets {cs}")

                    # Try all 24 permutations
                    best = 999
                    for perm in permutations(cs):
                        if list(perm) == cs: continue
                        flips = {}
                        for idx, v in enumerate([i,j,k,l]):
                            if perm[idx] != cs[idx]:
                                flips[v] = perm[idx]
                        E = excitation_energy(flips, vertex_coset, adj)
                        if E < best:
                            best = E
                            print(f"    Perm {cs}->{list(perm)}: E={E} (n_flips={len(flips)})")

                    print(f"    Best E = {best}")
                    found = True
                    break
            if found: break
        if found: break
    if found: break

# ============================================================
# Step 8: Physical interpretation
# ============================================================
print("\n--- Step 8: Physical interpretation ---")

# A single flip has E=3. This is the "bare mass" of a quasiparticle.
# Bound states have E < n*3 = lower mass than n free particles.
# The binding energy b = 3n - E > 0 means the bound state is stable.

# In hadron physics:
# - Mesons = quark-antiquark (n=2): E = 4, binding = 2
# - Baryons = 3 quarks (n=3): need 3-vertex configuration
#   - Triangle: E=6, binding = 3
#   - Chain (n=3, e=2): E = ?

# Check chain (3 vertices, 2 edges, no triangle)
print(f"\n  Chain excitation (3 vertices, 2 edges):")
for i in range(N):
    for j in nbr_list[i]:
        if j <= i: continue
        for k in nbr_list[j]:
            if k <= i or k == i: continue
            if adj[i, k]: continue  # We want a CHAIN, not triangle

            ci, cj, ck = vertex_coset[i], vertex_coset[j], vertex_coset[k]
            if len({ci, cj, ck}) < 3: continue

            # Try cyclic: i->cj, j->ck, k->ci
            E_cyc = excitation_energy({i: cj, j: ck, k: ci}, vertex_coset, adj)
            # Try swap + single: i<->j swap, k separate
            E_swap = excitation_energy({i: cj, j: ci, k: ci}, vertex_coset, adj)

            print(f"  Chain ({i},{j},{k}), cosets ({ci},{cj},{ck}):")
            print(f"    Cyclic: E={E_cyc}")
            print(f"    Swap i-j + flip k: E={E_swap}")
            break
        else: continue
        break
    else: continue
    break

# ============================================================
# Step 9: Excitation spectrum summary
# ============================================================
print("\n--- Step 9: Excitation spectrum ---")

print(f"""
  Excitation spectrum on 600-cell (antiferromagnetic q=5 Potts):

  n=1: E=3 (single flip, 'bare quasiparticle')
  n=2: E=4 (pair swap, 'meson-like', binding=2)
  n=3: E=6 (triangle cyclic, binding=3)

  General: E = 3n - S, where S counts matched internal edges

  KEY: The number 3 = neighbors per other coset
       = dim(S^3) = M_vacancy(exp109)

  Binding fraction: b/E_free = (3n-E)/(3n)
  n=2: 2/6 = 33%
  n=3: 3/9 = 33% (same!)

  Mass ratios:
  E(pair)/E(single) = 4/3 = 1.333
  E(triangle)/E(single) = 6/3 = 2.000
  E(triangle)/E(pair) = 6/4 = 1.500
""")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
