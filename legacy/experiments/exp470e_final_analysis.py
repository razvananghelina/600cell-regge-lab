"""
exp470e: FINAL ANALYSIS — 40 decompositions, Grassmannian angles
=================================================================
KEY FINDINGS FROM exp470d:
- 40 distinct 600-cell decompositions of E8
- Overlaps: 60, 72, 84, 96 (multiples of 12)
- Grassmannian singular values: 1/sqrt(5) and sqrt(2/5)
- All decompositions isomorphic

NOW: Deep structural analysis, framework connections
STATUS: EXPLORATORY
"""

import numpy as np
from itertools import product as iprod, permutations
from collections import Counter
import random

phi = (1 + np.sqrt(5)) / 2
a1 = 5

# ============================================================
# Setup (copied from exp470d)
# ============================================================
E8_list = []
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8); v[i] = si; v[j] = sj
                E8_list.append(v)
for signs in iprod([1, -1], repeat=8):
    if sum(1 for s in signs if s == -1) % 2 == 0:
        E8_list.append(np.array(signs) / 2.0)
E8 = np.array(E8_list)

simple = np.zeros((8, 8))
simple[0] = np.array([1, -1, -1, -1, -1, -1, -1, 1]) / 2.0
simple[1] = np.array([1, 1, 0, 0, 0, 0, 0, 0])
simple[2] = np.array([-1, 1, 0, 0, 0, 0, 0, 0])
simple[3] = np.array([0, -1, 1, 0, 0, 0, 0, 0])
simple[4] = np.array([0, 0, -1, 1, 0, 0, 0, 0])
simple[5] = np.array([0, 0, 0, -1, 1, 0, 0, 0])
simple[6] = np.array([0, 0, 0, 0, -1, 1, 0, 0])
simple[7] = np.array([0, 0, 0, 0, 0, -1, 1, 0])

def make_refl(alpha):
    return np.eye(8) - 2 * np.outer(alpha, alpha) / (alpha @ alpha)

def build_real_basis(eigvecs, indices):
    raw = []
    used = set()
    for j in indices:
        if j in used: continue
        v = eigvecs[:, j]
        for k in indices:
            if k > j and k not in used:
                if np.allclose(v, np.conj(eigvecs[:, k]), atol=1e-10):
                    raw.append(np.real(v))
                    raw.append(np.imag(v))
                    used.add(j); used.add(k)
                    break
        else:
            if j not in used:
                raw.append(np.real(v))
                used.add(j)
    B = np.array(raw).T
    Q, R = np.linalg.qr(B)
    return Q[:, :4]

def get_full_decomposition(perm):
    """Get decomposition AND projector from ordering."""
    C = np.eye(8)
    for i in perm:
        C = make_refl(simple[i]) @ C
    if not np.allclose(np.linalg.matrix_power(C, 30), np.eye(8)):
        return None
    eigvals, eigvecs = np.linalg.eig(C)
    args = np.angle(eigvals)
    h4_idx = []
    for i in range(8):
        m = round(args[i] * 30 / (2*np.pi))
        if m <= 0: m += 30
        if m in {1, 11, 19, 29}:
            h4_idx.append(i)
    if len(h4_idx) != 4: return None
    P = build_real_basis(eigvecs, h4_idx)
    if P.shape != (8, 4): return None
    proj = E8 @ P
    norms2 = np.sum(proj**2, axis=1)
    inner = frozenset(i for i in range(240) if norms2[i] < 1.0)
    outer = frozenset(i for i in range(240) if norms2[i] >= 1.0)
    if len(inner) != 120 or len(outer) != 120: return None
    return (inner, outer, P)

# ============================================================
# Collect ALL decompositions with projectors
# ============================================================
print("=" * 70)
print("COLLECTING ALL DECOMPOSITIONS AND PROJECTORS")
print("=" * 70)

decomp_data = {}  # partition -> list of (perm, projector)
n = 0
for perm in permutations(range(8)):
    result = get_full_decomposition(perm)
    n += 1
    if result is not None:
        inner, outer, P = result
        key = (inner, outer) if min(inner) < min(outer) else (outer, inner)
        if key not in decomp_data:
            decomp_data[key] = []
        decomp_data[key].append((perm, P))
    if n % 10000 == 0:
        print(f"  Processed {n}/40320...")

print(f"\nTotal distinct partitions: {len(decomp_data)}")

# ============================================================
# ANALYSIS 1: The number 40
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 1: THE NUMBER 40")
print(f"{'='*70}")

print(f"""
Found: 40 distinct 600-cell decompositions of E8.

Framework connections to test:
  40 = 8 * 5 = rank(E8) * a1        = {8 * a1}
  40 = 120 / 3 = N / N_gen           = {120 // 3}
  40 = (a1^2 + 1) + (a1^2 + 1)/2 + 1 = {26 + 13 + 1}  NO
  40 = 2 * (a1^2 - a1) = 2*20        = {2*(a1**2 - a1)} YES
  40 = dim(2-forms on R^8) - C(8,2) + 4 = 28-28+40 NO
  40 = |A5| / 3 = 60/3               = {60//3}  NO, 60/3=20
  40 = h(E8) + dim(H4) = 30 + 10     = 40 INTERESTING!

  h(E8) = 30 = Coxeter number
  dim of what? The Lie algebra of H4 would be dim(so(4)) = 6... no.
  The H4 ROOT SYSTEM has 120 roots. But dim(Gr(4,8)) = 16, not 10.

  Actually: 40 = C(a1,2) * 4 = 10*4  = {10*4}
  40 = C(a1,2) * (a1-1)              = {10*(a1-1)}

  Most intriguing: 40 = rank(E8) * a1 = 8 * 5
  Or: 40 = 2 * 4 * a1 = 2 * dim(R^4) * a1
""")

# ============================================================
# ANALYSIS 2: Overlap structure
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: OVERLAP GRAPH ON 40 DECOMPOSITIONS")
print(f"{'='*70}")

keys = list(decomp_data.keys())
n_decomp = len(keys)

# Compute all overlaps
overlap_matrix = np.zeros((n_decomp, n_decomp), dtype=int)
for i in range(n_decomp):
    for j in range(i+1, n_decomp):
        inner_i = keys[i][0]
        inner_j = keys[j][0]
        overlap = len(inner_i & inner_j)
        # Also check inner_i vs outer_j
        overlap_alt = len(inner_i & keys[j][1])
        # Use the larger overlap (one of them is > 60, the other < 60)
        ov = max(overlap, overlap_alt)
        overlap_matrix[i, j] = ov
        overlap_matrix[j, i] = ov

overlap_vals = Counter()
for i in range(n_decomp):
    for j in range(i+1, n_decomp):
        overlap_vals[overlap_matrix[i,j]] += 1

print(f"Overlap distribution between 40 decompositions:")
for ov, cnt in sorted(overlap_vals.items()):
    diff = 120 - ov
    print(f"  overlap = {ov}/120 (differ by {diff} = {diff//12}*12): {cnt} pairs")

total_pairs = n_decomp * (n_decomp - 1) // 2
print(f"Total pairs: {total_pairs}")

# Check: overlap = 60 means DISJOINT inner sets (60 = half of 120)
# Actually: overlap = 60 could mean one inner = other's outer
for i in range(n_decomp):
    for j in range(i+1, n_decomp):
        if overlap_matrix[i,j] == 60:
            # Check if inner_i = outer_j exactly
            if keys[i][0] == keys[j][1]:
                print(f"  Decomp {i} and {j}: inner/outer SWAP (exact)")

# ============================================================
# ANALYSIS 3: Grassmannian principal angles
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: GRASSMANNIAN ANGLES BETWEEN 4D SUBSPACES")
print(f"{'='*70}")

# For each distinct decomposition, pick one representative projector
rep_projectors = []
rep_keys = []
for key in keys:
    _, P = decomp_data[key][0]
    rep_projectors.append(P)
    rep_keys.append(key)

# Compute principal angles between all pairs
angle_data = {}
for i in range(min(40, n_decomp)):
    for j in range(i+1, min(40, n_decomp)):
        M = rep_projectors[i].T @ rep_projectors[j]
        svd = np.linalg.svd(M, compute_uv=False)
        svd_sorted = tuple(np.sort(svd)[::-1])
        angles = tuple(np.round(np.arccos(np.clip(svd, -1, 1)) * 180 / np.pi, 2))

        # Round for counting
        svd_key = tuple(np.round(svd_sorted, 4))
        if svd_key not in angle_data:
            angle_data[svd_key] = 0
        angle_data[svd_key] += 1

print(f"Distinct singular value patterns between {min(40, n_decomp)} projectors:")
for svd_key, cnt in sorted(angle_data.items(), key=lambda x: -x[1]):
    angles = [round(np.arccos(min(s, 1.0)) * 180 / np.pi, 2) for s in svd_key]
    print(f"  SVs = {svd_key}: {cnt} pairs  (angles = {angles} deg)")

# Check special values
print(f"\nFramework values check:")
print(f"  1/sqrt(5) = 1/sqrt(a1) = {1/np.sqrt(5):.6f}")
print(f"  sqrt(2/5) = sqrt(2/a1) = {np.sqrt(2/5):.6f} = c_bare (mass formula!)")
print(f"  arccos(1/sqrt(5)) = {np.arccos(1/np.sqrt(5))*180/np.pi:.2f} deg")
print(f"  arccos(sqrt(2/5)) = {np.arccos(np.sqrt(2/5))*180/np.pi:.2f} deg")

# ============================================================
# ANALYSIS 4: How many distinct 4D subspaces?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: DISTINCT 4D SUBSPACES (Grassmannian points)")
print(f"{'='*70}")

# Two projectors give the same subspace iff their SVs are all 1
# Count distinct subspaces
subspace_classes = {}
for i in range(n_decomp):
    found = False
    for cls_idx, cls_rep in subspace_classes.items():
        M = rep_projectors[i].T @ rep_projectors[cls_rep]
        svd = np.linalg.svd(M, compute_uv=False)
        if np.allclose(svd, 1.0, atol=1e-6):
            found = True
            break
    if not found:
        subspace_classes[i] = i

print(f"Distinct 4D subspaces among {n_decomp} decompositions: {len(subspace_classes)}")

# Actually, each DECOMPOSITION can have multiple Coxeter elements (orderings)
# giving different subspaces. Let me check within each decomposition.
for idx, key in enumerate(keys[:5]):
    projectors_for_key = [P for _, P in decomp_data[key]]
    n_proj = len(projectors_for_key)

    # Count distinct subspaces among these
    sub_classes = {}
    for pi in range(n_proj):
        found = False
        for ci, cr in sub_classes.items():
            M = projectors_for_key[pi].T @ projectors_for_key[cr]
            svd = np.linalg.svd(M, compute_uv=False)
            if np.allclose(svd, 1.0, atol=1e-6):
                found = True
                break
        if not found:
            sub_classes[pi] = pi

    print(f"  Decomp {idx}: {n_proj} orderings -> {len(sub_classes)} distinct subspaces")

# ============================================================
# ANALYSIS 5: Graph structure of the 40 decompositions
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: GRAPH ON 40 DECOMPOSITIONS")
print(f"{'='*70}")

# Build adjacency at different overlap levels
for threshold in [96, 84, 72]:
    adj = (overlap_matrix >= threshold).astype(int)
    np.fill_diagonal(adj, 0)
    degrees = np.sum(adj, axis=1)
    print(f"Adjacency at overlap >= {threshold}: degree distribution = {dict(Counter(degrees))}")

# Check: is the graph of 40 decompositions regular?
# At overlap = 96: adjacent decompositions differ by only 24 roots
adj_96 = (overlap_matrix == 96).astype(int)
np.fill_diagonal(adj_96, 0)
deg_96 = np.sum(adj_96, axis=1)
print(f"\nGraph at overlap = 96 (differ by 24 = 2*12 roots):")
print(f"  Degree distribution: {dict(Counter(deg_96))}")

adj_84 = (overlap_matrix == 84).astype(int)
np.fill_diagonal(adj_84, 0)
deg_84 = np.sum(adj_84, axis=1)
print(f"Graph at overlap = 84 (differ by 36 = 3*12 roots):")
print(f"  Degree distribution: {dict(Counter(deg_84))}")

adj_72 = (overlap_matrix == 72).astype(int)
np.fill_diagonal(adj_72, 0)
deg_72 = np.sum(adj_72, axis=1)
print(f"Graph at overlap = 72 (differ by 48 = 4*12 roots):")
print(f"  Degree distribution: {dict(Counter(deg_72))}")

adj_60 = (overlap_matrix == 60).astype(int)
np.fill_diagonal(adj_60, 0)
deg_60 = np.sum(adj_60, axis=1)
print(f"Graph at overlap = 60 (differ by 60 = 5*12 roots):")
print(f"  Degree distribution: {dict(Counter(deg_60))}")

# Total degree check
for i in range(n_decomp):
    total = sum(overlap_matrix[i, j] > 0 for j in range(n_decomp) if j != i)
assert total == n_decomp - 1  # all connected since overlap > 0

# ============================================================
# ANALYSIS 6: The 24-root exchange
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: STRUCTURE OF ROOT EXCHANGES")
print(f"{'='*70}")

# When two decompositions differ by 24 roots (overlap = 96):
# which 24 roots switch between inner and outer?
for i in range(n_decomp):
    for j in range(i+1, n_decomp):
        if overlap_matrix[i,j] == 96:
            inner_i = keys[i][0]
            inner_j = keys[j][0]

            # If overlap of inner_i with inner_j is 96, then 24 switch
            ov_ii = len(inner_i & inner_j)
            ov_io = len(inner_i & keys[j][1])

            if ov_ii == 96:
                switched_out = inner_i - inner_j  # 24 roots that left inner
                switched_in = inner_j - inner_i   # 24 roots that entered inner
            else:
                continue

            # What kind of E8 roots are these?
            out_list = sorted(switched_out)
            in_list = sorted(switched_in)

            # Check if they form a root subsystem
            out_roots = E8[out_list]
            in_roots = E8[in_list]

            IP_out = out_roots @ out_roots.T
            ip_out_vals = sorted(set(np.round(IP_out[np.triu_indices(24, k=1)], 1)))
            print(f"  24 switched-out roots: inner products = {ip_out_vals}")

            IP_in = in_roots @ in_roots.T
            ip_in_vals = sorted(set(np.round(IP_in[np.triu_indices(24, k=1)], 1)))
            print(f"  24 switched-in roots: inner products = {ip_in_vals}")

            # Check: are the 24 roots a 24-cell?
            # 24-cell has 24 vertices, all unit quaternions with integer coords
            # Inner products of 24-cell: -1 (12), -1/2 (32), 0 (18), 1/2 (32), 1 (12)...
            # Actually 24-cell has degree 8.

            # Count: number of roots that are antipodal pairs in the 24
            n_anti = 0
            for k in out_list:
                if any(np.allclose(E8[k], -E8[l]) for l in out_list):
                    n_anti += 1
            print(f"  Antipodal pairs in 24 switched roots: {n_anti}/24")

            # Adjacency in E8 (IP=1) among the 24
            adj_out = (np.abs(IP_out - 1) < 0.01).astype(int)
            np.fill_diagonal(adj_out, 0)
            deg_out = Counter(np.sum(adj_out, axis=1))
            print(f"  E8-adjacency degree in 24 switched: {dict(deg_out)}")

            # Just analyze first example
            break
    else:
        continue
    break

# ============================================================
# ANALYSIS 7: Spectral properties
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: SPECTRAL ANALYSIS OF DECOMPOSITION SPACE")
print(f"{'='*70}")

# Build the "exchange graph": nodes = 40 decompositions,
# edges weighted by overlap
W = overlap_matrix.astype(float)
np.fill_diagonal(W, 0)

# Normalized Laplacian
D = np.diag(W.sum(axis=1))
L = D - W

eigvals_L = np.sort(np.linalg.eigvalsh(L))
print(f"Laplacian eigenvalues of exchange graph (40 nodes):")
for k, ev in enumerate(eigvals_L):
    print(f"  lambda_{k} = {ev:.4f}")

# Check algebraic connectivity
print(f"\nAlgebraic connectivity (lambda_1) = {eigvals_L[1]:.4f}")
print(f"Spectral gap = {eigvals_L[1] / eigvals_L[-1]:.6f}")

# ============================================================
# FINAL COMPREHENSIVE SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("=" * 70)
print("COMPREHENSIVE SUMMARY: exp470 QUANTUM GRAVITY FROM E8")
print("=" * 70)
print(f"{'='*70}")

print(f"""
TASK 1: Number of 600-cell decompositions of E8
================================================
ANSWER: 40 DISTINCT DECOMPOSITIONS

The 240 E8 roots can be partitioned into two concentric 600-cells
("inner" at radius n1, "outer" at radius n2) in exactly 40 ways.

  n1^2 = 1 - 1/sqrt(a1)  (inner 600-cell)
  n2^2 = 1 + 1/sqrt(a1)  (outer 600-cell)
  n2/n1 = phi             (golden ratio!)
  n1^2 * n2^2 = 4/a1      (product involves a1 = 5)

  The number 40 = rank(E8) * a1 = 8 * 5

TASK 2: Graph structure comparison
===================================
ANSWER: ALL 40 DECOMPOSITIONS ARE ISOMORPHIC

- Inner and outer 600-cells have IDENTICAL E8-adjacency spectra
- Same inner product distribution within each set
- Same Gram matrix eigenvalues
- The two 600-cells within each decomposition are also isomorphic

The decompositions differ only in WHICH 120 roots are "inner" vs "outer".

TASK 3: Grassmannian structure
===============================
ANSWER: RICH STRUCTURE WITH FRAMEWORK CONSTANTS

The 40 decompositions correspond to different 4D subspaces in Gr(4,8).
Principal angles between subspaces involve:

  cos(theta) = 1/sqrt(a1)      = 1/sqrt(5)     [= {1/np.sqrt(5):.6f}]
  cos(theta) = sqrt(2/a1)      = sqrt(2/5)     [= {np.sqrt(2/5):.6f}]

  The SECOND value sqrt(2/a1) is EXACTLY c_bare from the mass formula!
  c_eff = sqrt(2/a1) * (1 - b1*alpha^2)

Overlap structure:
  Decompositions differ by 24, 36, 48, or 60 roots (multiples of 12).
  The 12-fold periodicity matches the vertex degree of the 600-cell.

TASK 4: Path integral potential
================================
STATUS: The finite number (40) of decompositions makes a discrete
"path integral" Z = sum_{{D}} exp(-S[D]) well-defined in principle.

  Z = sum over 40 decompositions D
  S[D] = spectral action of the projected 600-cell

Each decomposition gives the same spectral action (by isomorphism),
so Z = 40 * exp(-S_0). This is TRIVIAL -- no dynamics.

UNLESS the action depends on the 4D subspace (not just the partition).
The different 4D subspaces DO have different geometric properties
(different angles, different orientations in R^8).

TASK 5: Connection to Regge calculus
======================================
STATUS: The 40 decompositions could parameterize 40 discrete
gravitational configurations. But since all are isomorphic,
they give the same Regge curvature. The gravitational content
is in the ANGLES between subspaces, not in the 600-cell itself.

VERDICT:
========
The idea of quantum gravity from E8 rotations produces a FINITE
and SPECIFIC structure: 40 decompositions forming a graph with
overlaps in multiples of 12 and Grassmannian angles involving
sqrt(2/a1) = c_bare.

However, the gravitational interpretation is LIMITED because
all decompositions are isomorphic. The different 4D subspaces
provide a 40-element discrete geometry, but this is far from
a continuum graviton.

STATUS: INTERESTING STRUCTURAL RESULT (40 = rank(E8) * a1),
        NOT a viable quantum gravity mechanism.
        The Grassmannian angle sqrt(2/a1) = c_bare is intriguing
        but may be coincidental.
""")
