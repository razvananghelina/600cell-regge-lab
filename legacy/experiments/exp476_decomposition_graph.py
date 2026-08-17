"""
exp476: The Decomposition Graph -- Spectral Analysis of the 40 E8/H4 splittings
================================================================================

WHAT WE HAVE: 40 decompositions of 240 E8 roots into 2x 600-cell.
Complete graph K_40 with edges 4-colored by overlap m in {5,6,7,8}.
Pair counts: p5=176, p6=296, p7=184, p8=124.

WHAT WE DON'T KNOW:
  1. The spectral properties of the 4 adjacency matrices A_m
  2. Whether this is an association scheme (A_i @ A_j = sum c_{ij}^k A_k)
  3. The automorphism group
  4. Whether framework numbers appear in the spectra
  5. Whether there's structure beyond CKM

PLAN:
  Part 1: Build the 40 decompositions and the overlap matrix
  Part 2: Extract 4 adjacency matrices, compute spectra
  Part 3: Test association scheme axioms (Bose-Mesner algebra)
  Part 4: Hunt for framework numbers in eigenvalues
  Part 5: Analyze graph-theoretic properties (cliques, chromatic, etc.)
  Part 6: Connection to known combinatorial objects

Author: Razvan-Constantin Anghelina
Date: 2026-03-28
STATUS: EXPLORATORY
"""

import numpy as np
from itertools import permutations
from math import cos, sin, pi, sqrt, gcd, factorial
from collections import Counter
import sys

phi = (1 + sqrt(5)) / 2
phi_p = (1 - sqrt(5)) / 2
a1 = 5
b1 = 6
h = a1 * b1  # 30
N = 120
rank_E8 = 8
degree_600 = 12  # vertex degree of 600-cell

# ============================================================
# PART 1: BUILD ALL 40 DECOMPOSITIONS
# ============================================================
print("=" * 70)
print("PART 1: BUILDING 40 E8 -> 2x H4 DECOMPOSITIONS")
print("=" * 70)

# E8 roots
E8_roots = []
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8); v[i] = si; v[j] = sj
                E8_roots.append(v)
for bits in range(256):
    signs = [(-1)**((bits >> k) & 1) for k in range(8)]
    n_minus = sum(1 for s in signs if s < 0)
    if n_minus % 2 == 0:
        v = np.array(signs) / 2.0
        E8_roots.append(v)
E8 = np.array(E8_roots)
assert len(E8) == 240, f"Expected 240 E8 roots, got {len(E8)}"

# Simple roots (Bourbaki)
simple = np.zeros((8, 8))
simple[0] = np.array([1, -1, 0, 0, 0, 0, 0, 0])
simple[1] = np.array([0, 1, -1, 0, 0, 0, 0, 0])
simple[2] = np.array([0, 0, 1, -1, 0, 0, 0, 0])
simple[3] = np.array([0, 0, 0, 1, -1, 0, 0, 0])
simple[4] = np.array([0, 0, 0, 0, 1, -1, 0, 0])
simple[5] = np.array([0, 0, 0, 0, 0, 1, -1, 0])
simple[6] = np.array([0, 0, 0, 0, 0, 1, 1, 0])
simple[7] = np.array([-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5])

reflections = [np.eye(8) - 2 * np.outer(simple[i], simple[i]) / np.dot(simple[i], simple[i])
               for i in range(8)]

def make_coxeter(perm):
    C = np.eye(8)
    for i in perm:
        C = reflections[i] @ C
    return C

def get_H4_projector(C):
    if not np.allclose(np.linalg.matrix_power(C, 30), np.eye(8), atol=1e-8):
        return None
    eigvals, eigvecs = np.linalg.eig(C)
    args = np.angle(eigvals)
    exponents = []
    for i in range(8):
        m = round(args[i] * 30 / (2*pi)) % 30
        if m == 0: m = 30
        exponents.append((m, i))
    h4_idx = [i for m, i in exponents if m in {1, 11, 19, 29}]
    if len(h4_idx) != 4:
        return None
    v_raw = eigvecs[:, h4_idx]
    basis_real = []
    used = set()
    for k in range(4):
        if k in used: continue
        v = v_raw[:, k]
        if np.max(np.abs(np.imag(v))) > 0.01:
            basis_real.append(np.real(v))
            basis_real.append(np.imag(v))
            for k2 in range(k+1, 4):
                if k2 not in used and np.allclose(v_raw[:, k2], np.conj(v)):
                    used.add(k2)
                    break
        else:
            basis_real.append(np.real(v))
    if len(basis_real) != 4:
        return None
    B = np.array(basis_real).T
    Q, _ = np.linalg.qr(B)
    P = Q[:, :4]
    return P

# Find all 40 distinct decompositions
print("Searching for 40 distinct H4 subspaces...")
projectors = []
np.random.seed(42)
perms_to_try = list(permutations(range(8)))
np.random.shuffle(perms_to_try)

for perm in perms_to_try:
    C = make_coxeter(perm)
    P = get_H4_projector(C)
    if P is None:
        continue
    is_new = True
    for P_old in projectors:
        M = P_old.T @ P
        svs = np.linalg.svd(M, compute_uv=False)
        if np.allclose(svs, np.ones(4), atol=1e-6):
            is_new = False
            break
    if is_new:
        projectors.append(P)
    if len(projectors) >= 40:
        break

n_found = len(projectors)
print(f"Found {n_found}/40 H4 subspaces")
assert n_found == 40, f"Need all 40, got {n_found}"

# ============================================================
# PART 2: OVERLAP MATRIX AND ADJACENCY MATRICES
# ============================================================
print(f"\n{'='*70}")
print("PART 2: OVERLAP MATRIX AND 4 ADJACENCY MATRICES")
print("="*70)

# Compute full overlap matrix
overlap_matrix = np.zeros((40, 40), dtype=int)
for i in range(40):
    for j in range(i+1, 40):
        M = projectors[i].T @ projectors[j]
        svs = np.linalg.svd(M, compute_uv=False)
        ov = round(h * np.sum(svs**2))
        overlap_matrix[i, j] = ov
        overlap_matrix[j, i] = ov
    overlap_matrix[i, i] = N  # self-overlap = 120

# Verify overlap values
off_diag = []
for i in range(40):
    for j in range(i+1, 40):
        off_diag.append(overlap_matrix[i, j])

counts = Counter(off_diag)
print(f"\nOverlap distribution (off-diagonal):")
for ov in sorted(counts.keys()):
    m = ov // 12
    print(f"  overlap={ov} (m={m}): {counts[ov]} pairs")
print(f"  Total: {sum(counts.values())} = C(40,2) = {40*39//2}")

# Build adjacency matrices for each m-level
A = {}
for m in [5, 6, 7, 8]:
    A[m] = np.zeros((40, 40), dtype=int)
    for i in range(40):
        for j in range(40):
            if overlap_matrix[i, j] == 12 * m and i != j:
                A[m][i, j] = 1

# Verify: sum of A[m] = J - I (complete graph minus diagonal)
A_sum = sum(A[m] for m in [5,6,7,8])
J_minus_I = np.ones((40, 40), dtype=int) - np.eye(40, dtype=int)
assert np.array_equal(A_sum, J_minus_I), "Adjacency matrices don't sum to J-I!"
print("\nVerification: A5 + A6 + A7 + A8 = J - I: PASS")

# Degree sequences
print(f"\nDegree sequences (valencies per node):")
for m in [5, 6, 7, 8]:
    degrees = np.sum(A[m], axis=1)
    unique_degs = sorted(set(degrees))
    print(f"  A_{m}: degrees = {unique_degs}, count = {Counter(degrees.tolist())}")

# ============================================================
# PART 3: SPECTRA OF ADJACENCY MATRICES
# ============================================================
print(f"\n{'='*70}")
print("PART 3: SPECTRAL ANALYSIS")
print("="*70)

all_spectra = {}
for m in [5, 6, 7, 8]:
    eigvals = np.sort(np.linalg.eigvalsh(A[m].astype(float)))[::-1]
    all_spectra[m] = eigvals

    # Find distinct eigenvalues
    distinct = []
    for ev in eigvals:
        if not distinct or abs(ev - distinct[-1][0]) > 1e-6:
            distinct.append([ev, 1])
        else:
            distinct[-1][1] += 1

    print(f"\nSpectrum of A_{m}:")
    print(f"  {'eigenvalue':>12} {'mult':>5} {'as fraction':>20}")
    print(f"  {'-'*40}")
    for ev, mult in distinct:
        # Try to identify as framework number
        frac_str = ""
        if abs(ev - round(ev)) < 1e-6:
            frac_str = f"= {round(ev)}"
        elif abs(ev * a1 - round(ev * a1)) < 1e-4:
            frac_str = f"~ {round(ev*a1)}/{a1}"
        elif abs(ev * b1 - round(ev * b1)) < 1e-4:
            frac_str = f"~ {round(ev*b1)}/{b1}"
        elif abs(ev / phi - round(ev / phi)) < 1e-4:
            frac_str = f"~ {round(ev/phi)}*phi"
        elif abs(ev * phi - round(ev * phi)) < 1e-4:
            frac_str = f"~ {round(ev*phi)}/phi"
        print(f"  {ev:12.6f} {mult:5d} {frac_str:>20}")

    print(f"  Trace = {np.sum(eigvals):.6f} (should be 0)")

# ============================================================
# PART 4: ASSOCIATION SCHEME TEST
# ============================================================
print(f"\n{'='*70}")
print("PART 4: ASSOCIATION SCHEME TEST (Bose-Mesner algebra)")
print("="*70)

# For an association scheme: A_i @ A_j = sum_k p_{ij}^k A_k
# This means the span of {I, A_5, A_6, A_7, A_8} is closed under matrix multiplication

I40 = np.eye(40)
basis = [I40, A[5].astype(float), A[6].astype(float), A[7].astype(float), A[8].astype(float)]
basis_names = ['I', 'A5', 'A6', 'A7', 'A8']

# Compute all products A_i @ A_j and check if they're in the span
print("\nProducts A_i @ A_j -- checking if in span{I, A5, A6, A7, A8}:")
is_scheme = True
structure_constants = {}

for i_idx in range(1, 5):  # A5, A6, A7, A8
    for j_idx in range(i_idx, 5):
        prod = basis[i_idx] @ basis[j_idx]

        # Try to decompose: prod = c0*I + c1*A5 + c2*A6 + c3*A7 + c4*A8
        # Use the fact that these are 0/1 matrices with disjoint support
        # The product may have arbitrary integer entries

        # Method: for each pair (r,s), prod[r,s] should equal sum_k c_k * basis[k][r,s]
        # Since basis matrices have disjoint 0/1 patterns (except overlaps in product),
        # we can read off coefficients from specific positions

        # Better: least squares
        n_basis = len(basis)
        B_mat = np.column_stack([b.flatten() for b in basis])
        coeffs, residual, _, _ = np.linalg.lstsq(B_mat, prod.flatten(), rcond=None)

        # Check quality of fit
        reconstructed = sum(c * b for c, b in zip(coeffs, basis))
        max_err = np.max(np.abs(prod - reconstructed))

        m_i = [5,6,7,8][i_idx-1]
        m_j = [5,6,7,8][j_idx-1]

        # Round coefficients
        coeffs_r = np.round(coeffs).astype(int)
        reconstructed_r = sum(c * b for c, b in zip(coeffs_r, basis))
        max_err_r = np.max(np.abs(prod - reconstructed_r))

        status = "OK" if max_err_r < 0.5 else "FAIL"
        if max_err_r >= 0.5:
            is_scheme = False

        coeff_str = " + ".join(f"{c}*{n}" for c, n in zip(coeffs_r, basis_names) if c != 0)
        print(f"  A{m_i} @ A{m_j} = {coeff_str}  (max_err={max_err_r:.4f}) [{status}]")
        structure_constants[(m_i, m_j)] = dict(zip(basis_names, coeffs_r))

if is_scheme:
    print(f"\n*** RESULT: The 4-colored graph IS an association scheme! ***")
else:
    print(f"\n*** RESULT: NOT an association scheme (products not in span). ***")
    print(f"Checking if it's a COHERENT CONFIGURATION (weaker condition)...")

    # For coherent configuration: products are non-negative integer combos
    # but the graph need not be vertex-transitive
    # Actually let's check if graph is vertex-transitive by looking at degree sequence
    all_regular = True
    for m in [5,6,7,8]:
        degs = np.sum(A[m], axis=1)
        if len(set(degs)) > 1:
            all_regular = False
            break
    print(f"  All color classes regular? {all_regular}")

# ============================================================
# PART 5: FRAMEWORK NUMBER HUNT IN SPECTRA
# ============================================================
print(f"\n{'='*70}")
print("PART 5: FRAMEWORK NUMBER HUNT")
print("="*70)

# Collect all distinct eigenvalues
all_eigs = set()
for m in [5,6,7,8]:
    for ev in all_spectra[m]:
        # Round to avoid duplicates
        ev_r = round(ev, 6)
        all_eigs.add(ev_r)

all_eigs = sorted(all_eigs)
print(f"\nAll distinct eigenvalues across A5..A8: {len(all_eigs)}")

# Framework numbers to test against
framework = {
    'a1': 5, 'b1': 6, 'N': 120, 'h': 30, 'rank': 8, 'degree': 12,
    'phi': phi, 'phi^2': phi**2, 'phi^3': phi**3,
    'n_12': 3, 'n_23': 7, 'n_13': 12, 'n_ut': 4, 'n_db': 2,
    'N_gen': 3, 'N_eig': 9, 'alpha_s_inv': 2*phi**3,
    'sqrt(a1)': sqrt(5), 'a1^2': 25, 'a1+1': 6, 'a1-1': 4,
    '2*a1': 10, '(a1^2+1)': 26, 'b1/a1': 6/5,
    'sin2tW_inv': 26/6,
}

print(f"\nChecking eigenvalues against framework numbers:")
for ev in all_eigs:
    if abs(ev) < 1e-6:
        continue
    matches = []
    for name, val in framework.items():
        if val == 0: continue
        ratio = ev / val
        if abs(ratio - round(ratio)) < 0.01 and abs(round(ratio)) <= 20:
            matches.append(f"{round(ratio)}*{name}")
        if abs(1/ratio - round(1/ratio)) < 0.01 and abs(round(1/ratio)) <= 20:
            matches.append(f"{name}/{round(1/ratio)}")
    if matches:
        print(f"  {ev:10.4f}: {', '.join(matches[:5])}")

# ============================================================
# PART 6: COMBINED SPECTRA AND OVERLAP OPERATOR
# ============================================================
print(f"\n{'='*70}")
print("PART 6: OVERLAP OPERATOR SPECTRUM")
print("="*70)

# The overlap operator: O = sum_m m * A_m (weighted by overlap level)
O = sum(m * A[m].astype(float) for m in [5,6,7,8])
O_eigs = np.sort(np.linalg.eigvalsh(O))[::-1]

distinct_O = []
for ev in O_eigs:
    if not distinct_O or abs(ev - distinct_O[-1][0]) > 1e-4:
        distinct_O.append([ev, 1])
    else:
        distinct_O[-1][1] += 1

print(f"\nOverlap operator O = 5*A5 + 6*A6 + 7*A7 + 8*A8:")
print(f"  {'eigenvalue':>12} {'mult':>5}")
print(f"  {'-'*20}")
for ev, mult in distinct_O:
    # Check if eigenvalue/40 or /39 is a framework number
    frac_str = ""
    for denom in [1, a1, b1, N, h, rank_E8, 39, 40]:
        ratio = ev / denom
        if abs(ratio - round(ratio, 2)) < 0.01:
            if abs(ratio - round(ratio)) < 0.01:
                frac_str += f" = {round(ratio)}*{denom}"
    print(f"  {ev:12.6f} {mult:5d} {frac_str}")

# Also try: exchanged operator E = sum_m (10-m)*A_m (counting exchanged roots / 12)
E_op = sum((2*a1 - m) * A[m].astype(float) for m in [5,6,7,8])
E_eigs = np.sort(np.linalg.eigvalsh(E_op))[::-1]

distinct_E = []
for ev in E_eigs:
    if not distinct_E or abs(ev - distinct_E[-1][0]) > 1e-4:
        distinct_E.append([ev, 1])
    else:
        distinct_E[-1][1] += 1

print(f"\nExchanged operator E = 5*A5 + 4*A6 + 3*A7 + 2*A8:")
print(f"  {'eigenvalue':>12} {'mult':>5}")
print(f"  {'-'*20}")
for ev, mult in distinct_E:
    print(f"  {ev:12.6f} {mult:5d}")

# ============================================================
# PART 7: GRAPH STRUCTURE
# ============================================================
print(f"\n{'='*70}")
print("PART 7: GRAPH-THEORETIC PROPERTIES")
print("="*70)

# For each color, find clique and independence numbers (approximation)
for m in [5,6,7,8]:
    adj = A[m]
    n_edges = np.sum(adj) // 2
    avg_deg = np.mean(np.sum(adj, axis=1))
    max_deg = np.max(np.sum(adj, axis=1))
    min_deg = np.min(np.sum(adj, axis=1))

    # Find a maximal clique (greedy)
    nodes = list(range(40))
    clique = [nodes[0]]
    for n in nodes[1:]:
        if all(adj[n, c] == 1 for c in clique):
            clique.append(n)

    # Complement: independent set
    complement = 1 - adj - np.eye(40)
    ind_set = [0]
    for n in range(1, 40):
        if all(adj[n, c] == 0 for c in ind_set) and n not in ind_set:
            ind_set.append(n)

    print(f"\nA_{m} (overlap = {12*m}):")
    print(f"  Edges: {n_edges}, Avg degree: {avg_deg:.1f}, Range: [{min_deg}, {max_deg}]")
    print(f"  Greedy clique: {len(clique)}, Greedy independent set: {len(ind_set)}")

# ============================================================
# PART 8: COMMON EIGENSPACES
# ============================================================
print(f"\n{'='*70}")
print("PART 8: SIMULTANEOUS DIAGONALIZATION")
print("="*70)

# If it's an association scheme, all A_m share eigenspaces
# Check if A_m commute
print("\nCommutativity check [A_i, A_j]:")
for i in [5,6,7,8]:
    for j in range(i+1, 9):
        if j > 8: break
        comm = A[i] @ A[j] - A[j] @ A[i]
        max_comm = np.max(np.abs(comm))
        print(f"  [A{i}, A{j}] max = {max_comm:.6f}  {'COMMUTE' if max_comm < 0.5 else 'NONCOMMUTATIVE'}")

# If they commute, simultaneously diagonalize
# Use the overlap operator as primary
O_eigvals, O_eigvecs = np.linalg.eigh(O.astype(float))
idx = np.argsort(O_eigvals)[::-1]
O_eigvals = O_eigvals[idx]
O_eigvecs = O_eigvecs[:, idx]

# Project each A_m into the eigenspaces of O
print(f"\nA_m eigenvalues in each O-eigenspace:")
print(f"  {'O-eigval':>10} {'mult':>5} | {'A5':>8} {'A6':>8} {'A7':>8} {'A8':>8} | {'sum':>8}")
print(f"  {'-'*65}")

# Group O eigenvalues by multiplicity
O_groups = []
i = 0
while i < 40:
    ev = O_eigvals[i]
    j = i
    while j < 40 and abs(O_eigvals[j] - ev) < 0.01:
        j += 1
    O_groups.append((ev, i, j))
    i = j

for ev, start, end in O_groups:
    mult = end - start
    V = O_eigvecs[:, start:end]

    am_vals = []
    for m in [5,6,7,8]:
        # A_m restricted to this eigenspace
        Am_proj = V.T @ A[m].astype(float) @ V
        # If diagonal (scalar multiple of identity), read off eigenvalue
        am_eig = np.trace(Am_proj) / mult
        am_vals.append(am_eig)

    total = sum(am_vals)
    print(f"  {ev:10.4f} {mult:5d} | {am_vals[0]:8.3f} {am_vals[1]:8.3f} {am_vals[2]:8.3f} {am_vals[3]:8.3f} | {total:8.3f}")

# ============================================================
# PART 9: PAIR COUNT ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("PART 9: DEEP PAIR COUNT ANALYSIS")
print("="*70)

p = {5: 176, 6: 296, 7: 184, 8: 124}
total = sum(p.values())

print(f"\nPair counts and ratios:")
for m in [5,6,7,8]:
    print(f"  p_{m} = {p[m]:4d}  ({100*p[m]/total:.1f}%)  p_{m}/4 = {p[m]/4:.1f}")

print(f"\nPair count differences:")
print(f"  p6 - p5 = {p[6]-p[5]} = {p[6]-p[5]}")
print(f"  p6 - p7 = {p[6]-p[7]} = {p[6]-p[7]}")
print(f"  p7 - p8 = {p[7]-p[8]} = {p[7]-p[8]}")
print(f"  p5 - p8 = {p[5]-p[8]} = {p[5]-p[8]}")

print(f"\nFramework decompositions:")
print(f"  p5 = {p[5]} = 8*22 = rank*22")
print(f"  p6 = {p[6]} = 8*37 = rank*37")
print(f"  p7 = {p[7]} = 8*23 = rank*23")
print(f"  p8 = {p[8]} = 4*31 = n_ut*31")
print(f"  p5+p8 = {p[5]+p[8]} = 12*{(p[5]+p[8])//12} = degree*25 = degree*a1^2")
print(f"  p6+p7 = {p[6]+p[7]} = 12*{(p[6]+p[7])//12} = degree*40 = degree*rank*a1")
print(f"  p5+p7 = {p[5]+p[7]} = {p[5]+p[7]}")
print(f"  p6+p8 = {p[6]+p[8]} = {p[6]+p[8]}")

# Weighted sums
print(f"\nWeighted sums:")
S1 = sum(m * p[m] for m in [5,6,7,8])
S2 = sum(m**2 * p[m] for m in [5,6,7,8])
S_exch = sum((2*a1-m) * p[m] for m in [5,6,7,8])
print(f"  sum m*p_m = {S1}")
print(f"  sum m*p_m / C(40,2) = {S1/780:.6f} (average m)")
print(f"  sum (2a1-m)*p_m = {S_exch} (weighted by multipliers)")
print(f"  sum m^2*p_m = {S2}")

# Check: average m * 780 should equal S1
avg_m = S1 / 780
print(f"\n  Average overlap level: {avg_m:.6f}")
print(f"  = {avg_m:.6f}")
for name, val in [('b1', 6), ('h/a1', 6), ('(a1+b1+n_23+rank)/4', (5+6+7+8)/4)]:
    print(f"    vs {name} = {val}: diff = {abs(avg_m - val):.6f}")

# ============================================================
# PART 10: SUMMARY AND HONEST ASSESSMENT
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY AND HONEST ASSESSMENT")
print("="*70)
