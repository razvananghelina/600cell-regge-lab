"""
EXP-418: Full Coexact Spectrum of Hodge-1 Laplacian on 600-cell
================================================================
Goal: Compute the COMPLETE coexact spectrum of Delta_1 on the 600-cell
and determine if a spin-2 graviton mode (multiplicity 9) exists.

Key question (from quantum_gravity.txt):
- The coexact sector has dim 601. What are ALL multiplicities?
- Is there a mult-9 eigenspace? (= symmetric traceless 2-tensor = graviton)
- If not, graviton is COMPOSITE, not fundamental on 600-cell.

Method:
- Build 600-cell simplicial complex (120V, 720E, 1200F, 600C)
- Compute boundary operators d0 (V->E), d1 (E->F)
- Full dense diagonalization of d1^T*d1 (720x720)
- Nonzero eigenvalues = coexact spectrum (601 values)
- Group by value, report ALL multiplicities

Under H4 (order 14400): irrep dims include 1,4,5,6,8,9,10,16,18,20,24,25,30,36,40,48
Graviton = Sym^2_0(standard rep) = 9-dim irrep of H4
"""

import numpy as np
from itertools import combinations, product, permutations
from collections import defaultdict
import sys

PHI = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXP-418: Full Coexact Spectrum - Graviton Search")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell vertices
# ============================================================
print("\n--- Step 1: Build 600-cell ---")

verts_set = set()

# Type A: axis vertices (8)
for i in range(4):
    for s in [1, -1]:
        v = [0,0,0,0]; v[i] = s
        verts_set.add(tuple(v))

# Type B: half-integer vertices (16)
for signs in product([0.5, -0.5], repeat=4):
    verts_set.add(tuple(signs))

# Type C: golden ratio vertices (96)
vals_base = [PHI/2, 0.5, 1/(2*PHI), 0]
even_perms = []
for p in permutations(range(4)):
    inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
    if inv % 2 == 0:
        even_perms.append(p)

for perm in even_perms:
    base = [vals_base[perm[i]] for i in range(4)]
    nz = [i for i in range(4) if base[i] != 0]
    for signs in product([1,-1], repeat=len(nz)):
        v = list(base)
        for idx, s in zip(nz, signs):
            v[idx] *= s
        verts_set.add(tuple(np.round(v, 10)))

verts = np.array(sorted(verts_set))
N_v = len(verts)
print(f"  Vertices: {N_v}")
assert N_v == 120, f"Expected 120 vertices, got {N_v}"

# ============================================================
# Step 2: Build edges (720)
# ============================================================
print("  Building edges...")
dot_threshold = PHI / 2
dots = verts @ verts.T

edges = []
edge_to_idx = {}
adj_list = defaultdict(set)
for i in range(N_v):
    for j in range(i+1, N_v):
        if abs(dots[i,j] - dot_threshold) < 0.001:
            idx = len(edges)
            edges.append((i, j))
            edge_to_idx[(i, j)] = idx
            edge_to_idx[(j, i)] = idx
            adj_list[i].add(j)
            adj_list[j].add(i)

N_e = len(edges)
print(f"  Edges: {N_e}")
assert N_e == 720, f"Expected 720 edges, got {N_e}"

# Verify regularity
degrees = [len(adj_list[i]) for i in range(N_v)]
print(f"  Degree: {min(degrees)}-{max(degrees)} (should be 12)")

# ============================================================
# Step 3: Build faces (1200 triangles)
# ============================================================
print("  Building faces...")
triangles = []
face_to_idx = {}
for i in range(N_v):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    idx = len(triangles)
                    triangles.append((i, j, k))
                    face_to_idx[(i, j, k)] = idx

N_f = len(triangles)
print(f"  Faces: {N_f}")
assert N_f == 1200, f"Expected 1200 faces, got {N_f}"

# ============================================================
# Step 4: Boundary operators (dense for full diag)
# ============================================================
print("\n--- Step 2: Boundary operators ---")

# d0: N_e x N_v
print("  Building d0 (120 -> 720)...")
d0 = np.zeros((N_e, N_v))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = 1.0

# d1: N_f x N_e
print("  Building d1 (720 -> 1200)...")
d1 = np.zeros((N_f, N_e))
for f_idx, (i, j, k) in enumerate(triangles):
    e_ij = edge_to_idx.get((i, j))
    e_ik = edge_to_idx.get((i, k))
    e_jk = edge_to_idx.get((j, k))
    if e_ij is not None and e_ik is not None and e_jk is not None:
        d1[f_idx, e_ij] = 1.0
        d1[f_idx, e_jk] = 1.0
        d1[f_idx, e_ik] = -1.0

# Verify chain complex
check = d1 @ d0
print(f"  |d1*d0| = {np.max(np.abs(check)):.2e} (should be 0)")

# ============================================================
# Step 5: Full spectrum of Delta_0 (graph Laplacian)
# ============================================================
print("\n--- Step 3: Graph Laplacian spectrum (= exact part of Delta_1) ---")

Delta_0 = d0.T @ d0  # 120 x 120
evals_0 = np.linalg.eigvalsh(Delta_0)
evals_0 = np.sort(evals_0)

# Group eigenvalues
def group_eigenvalues(evals, tol=0.001):
    """Group eigenvalues by value, return (value, multiplicity) pairs."""
    groups = []
    i = 0
    while i < len(evals):
        val = evals[i]
        count = 1
        while i + count < len(evals) and abs(evals[i + count] - val) < tol:
            count += 1
        groups.append((np.mean(evals[i:i+count]), count))
        i += count
    return groups

groups_0 = group_eigenvalues(evals_0)
print(f"  Delta_0 spectrum ({N_v} eigenvalues):")
for val, mult in groups_0:
    # Try a + b*phi form
    phi_str = ""
    for a in range(-30, 31):
        for b in range(-30, 31):
            if abs(a + b * PHI - val) < 0.01:
                phi_str = f"  = {a}+{b}*phi" if b != 0 else f"  = {a}"
                break
        if phi_str:
            break
    print(f"    lambda = {val:10.5f}  mult = {mult:3d}{phi_str}")

print(f"\n  Nonzero eigenvalues: {sum(m for v,m in groups_0 if v > 0.01)}")
print(f"  Zero eigenvalues: {sum(m for v,m in groups_0 if v < 0.01)} (= b0 = 1)")

# ============================================================
# Step 6: Coexact Laplacian = d1^T * d1
# ============================================================
print("\n--- Step 4: Coexact Laplacian (d1^T * d1) ---")

L_coexact = d1.T @ d1  # 720 x 720
print(f"  Shape: {L_coexact.shape}")
print(f"  Computing full eigendecomposition (720 x 720)...")

evals_coexact_full, evecs_coexact = np.linalg.eigh(L_coexact)
evals_coexact_full = np.sort(evals_coexact_full)

# Separate zero and nonzero
n_zero = np.sum(np.abs(evals_coexact_full) < 0.01)
n_nonzero = N_e - n_zero
print(f"  Zero eigenvalues: {n_zero} (= dim(exact) + b1 = {N_v - 1} + 0 = {N_v - 1})")
print(f"  Nonzero eigenvalues: {n_nonzero} (= dim(coexact) = 601?)")

coexact_evals = evals_coexact_full[np.abs(evals_coexact_full) > 0.01]
print(f"  Coexact eigenvalue range: [{coexact_evals[0]:.5f}, {coexact_evals[-1]:.5f}]")

# ============================================================
# Step 7: FULL coexact spectrum with multiplicities
# ============================================================
print("\n--- Step 5: FULL COEXACT SPECTRUM ---")
print("=" * 70)

groups_coexact = group_eigenvalues(coexact_evals)
print(f"  {len(groups_coexact)} distinct eigenvalues, {n_nonzero} total (should be 601)")
print()

total_check = 0
has_mult_9 = False
mult_list = []

for idx, (val, mult) in enumerate(groups_coexact):
    total_check += mult
    mult_list.append(mult)

    # Try a + b*phi form
    phi_str = ""
    for a in range(-50, 51):
        for b in range(-50, 51):
            if abs(a + b * PHI - val) < 0.01:
                phi_str = f"{a}+{b}*phi" if b != 0 else f"{a}"
                break
        if phi_str:
            break

    marker = ""
    if mult == 9:
        marker = "  <-- GRAVITON CANDIDATE (mult=9)!"
        has_mult_9 = True
    elif mult == 5:
        marker = "  <-- (mult=5, could be vector)"
    elif mult == 1:
        marker = "  <-- (scalar)"

    print(f"  [{idx+1:2d}] lambda = {val:10.5f}  mult = {mult:3d}  [{phi_str}]{marker}")

print(f"\n  Total eigenvalues: {total_check} (should be {n_nonzero})")
print(f"  Distinct eigenvalues: {len(groups_coexact)}")

# ============================================================
# Step 8: Multiplicity analysis
# ============================================================
print("\n--- Step 6: Multiplicity Analysis ---")
print("=" * 70)

from collections import Counter
mult_counts = Counter(mult_list)
print("  Multiplicity distribution:")
for m in sorted(mult_counts.keys()):
    count = mult_counts[m]
    print(f"    mult={m:3d}: appears {count} time(s)")

print(f"\n  H4 irrep dimensions: 1, 4, 5, 6, 8, 9, 10, 16, 18, 20, 24, 25, 30, 36, 40, 48")
print(f"  Graviton (Sym^2_0 of standard rep) = dim 9")

if has_mult_9:
    print("\n  *** RESULT: Multiplicity 9 FOUND in coexact spectrum! ***")
    print("  *** This is the GRAVITON CANDIDATE ***")
    # Find which eigenvalue
    for val, mult in groups_coexact:
        if mult == 9:
            print(f"  *** Eigenvalue: {val:.6f} ***")
else:
    print("\n  *** RESULT: NO multiplicity 9 in coexact spectrum ***")
    print("  *** Graviton does NOT appear as fundamental mode ***")
    print("  *** Implication: graviton is COMPOSITE on 600-cell ***")

# Check if 9 appears as sum of smaller irreps in same eigenspace
print("\n  Check: could 9 arise as sum of irreps at same eigenvalue?")
for val, mult in groups_coexact:
    if mult > 9:
        # Check if 9 divides into it or is contained
        remainder = mult
        n9 = remainder // 9
        if n9 > 0 and remainder % 9 == 0:
            print(f"    lambda={val:.4f}, mult={mult} = {n9}x9")
        # Check if mult = 9 + something
        for d in [1, 4, 5, 6, 8, 10, 16]:
            if mult - d == 9:
                print(f"    lambda={val:.4f}, mult={mult} = 9 + {d}")
            if mult - 2*d == 9:
                print(f"    lambda={val:.4f}, mult={mult} = 9 + 2x{d}")

# ============================================================
# Step 9: Exact spectrum (for comparison)
# ============================================================
print("\n--- Step 7: Exact spectrum (from Delta_0) ---")

L_exact = d0 @ d0.T  # 720 x 720
evals_exact_full = np.linalg.eigvalsh(L_exact)
evals_exact_full = np.sort(evals_exact_full)

exact_evals = evals_exact_full[np.abs(evals_exact_full) > 0.01]
groups_exact = group_eigenvalues(exact_evals)

print(f"  {len(groups_exact)} distinct exact eigenvalues, {len(exact_evals)} total")
print(f"  (Should match nonzero spectrum of Delta_0 = graph Laplacian)")
for val, mult in groups_exact:
    phi_str = ""
    for a in range(-30, 31):
        for b in range(-30, 31):
            if abs(a + b * PHI - val) < 0.01:
                phi_str = f"  [{a}+{b}*phi]"
                break
        if phi_str:
            break
    print(f"    lambda = {val:10.5f}  mult = {mult:3d}{phi_str}")

# ============================================================
# Step 10: Consistency check
# ============================================================
print("\n--- Step 8: Consistency ---")

Delta_1 = d0 @ d0.T + d1.T @ d1  # Full Hodge-1 Laplacian
evals_1_full = np.linalg.eigvalsh(Delta_1)
evals_1_full = np.sort(evals_1_full)

# Delta_1 = exact + coexact (orthogonal)
# Full spectrum should be union of exact and coexact spectra
print(f"  Delta_1: {N_e} eigenvalues")
print(f"  Range: [{evals_1_full[0]:.5f}, {evals_1_full[-1]:.5f}]")
print(f"  Zero modes: {np.sum(np.abs(evals_1_full) < 0.01)} (= b1 = 0)")

# Full Delta_1 spectrum with multiplicities
groups_1 = group_eigenvalues(evals_1_full[evals_1_full > 0.01])
print(f"\n  Full Delta_1 nonzero spectrum ({len(groups_1)} distinct):")
for val, mult in groups_1:
    # Identify if exact, coexact, or both
    in_exact = any(abs(val - ve) < 0.02 for ve, me in groups_exact)
    in_coexact = any(abs(val - vc) < 0.02 for vc, mc in groups_coexact)
    src = ""
    if in_exact and in_coexact:
        src = " [BOTH exact+coexact]"
    elif in_exact:
        src = " [exact]"
    elif in_coexact:
        src = " [coexact]"

    phi_str = ""
    for a in range(-50, 51):
        for b in range(-50, 51):
            if abs(a + b * PHI - val) < 0.01:
                phi_str = f"  = {a}+{b}*phi" if b != 0 else f"  = {a}"
                break
        if phi_str:
            break

    print(f"    lambda = {val:10.5f}  mult = {mult:3d}{phi_str}{src}")

# ============================================================
# Step 11: Key physical quantities
# ============================================================
print("\n--- Step 9: Physical Interpretation ---")
print("=" * 70)

# Mass gap
print(f"  Coexact mass gap (lowest nonzero): {coexact_evals[0]:.6f}")
gap_val = coexact_evals[0]
for a in range(-20, 21):
    for b in range(-20, 21):
        if abs(a + b * PHI - gap_val) < 0.01:
            print(f"    = {a} + {b}*phi = {a + b*PHI:.6f}")

# Total DOF count
print(f"\n  DOF count:")
print(f"    Exact (gauge): {len(exact_evals)}")
print(f"    Coexact (physical): {len(coexact_evals)}")
print(f"    Harmonic (topological): 0")
print(f"    Total: {len(exact_evals) + len(coexact_evals)} (= {N_e})")

# In Regge calculus: edge length perturbations = 1-cochains
# Gauge = vertex displacements = exact 1-forms
# Physical DOF = coexact 1-forms
print(f"\n  Regge calculus interpretation:")
print(f"    Total metric DOF (edge lengths): {N_e}")
print(f"    Gauge DOF (diffeomorphisms): {N_v - 1} = {N_v}-1")
print(f"    Physical DOF: {N_e - N_v + 1} = 601")
print(f"    Per vertex: 601/120 = {601/120:.2f}")
print(f"    Continuum 4D graviton: 2 DOF per point")
print(f"    Continuum 3D on S^3: tower of TT modes")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
600-cell simplicial complex: V=120, E=720, F=1200, C=600
Euler char: 0 (S^3). Betti: (1,0,0,1).

Hodge-1 Laplacian Delta_1 (720 x 720):
  Exact sector: {len(exact_evals)} eigenvalues ({len(groups_exact)} distinct)
  Coexact sector: {len(coexact_evals)} eigenvalues ({len(groups_coexact)} distinct)
  Harmonic: 0 (b1=0)

Coexact multiplicities: {sorted(mult_list)}

Multiplicity 9 found: {has_mult_9}
""")

if has_mult_9:
    print("CONCLUSION: Graviton spin-2 mode EXISTS as fundamental coexact mode.")
    print("STATUS: DERIVED")
else:
    print("CONCLUSION: NO spin-2 graviton in fundamental coexact spectrum.")
    print("The graviton on 600-cell must be COMPOSITE (emergent from lower modes).")
    print("STATUS: NEGATIVE (but informative)")
    print()
    print("This means: discrete gravity on 600-cell is NOT like lattice QCD")
    print("(where gluons are fundamental links). Instead, the graviton emerges")
    print("only in the continuum limit from combinations of edge modes.")
