"""
EXP-502: Algebraic Proof that n96 = dim(E8)/2
================================================
GOAL: Prove that exactly 124 = dim(E8)/2 pairs of Coxeter H4
subspaces share a 3D subspace with cos^2 = 1/a1 in the 4th direction.

APPROACH:
The 40 decompositions come from Coxeter elements c = s_{p(1)}...s_{p(8)}.
Two Coxeter elements give the SAME H4 subspace iff they're in the same
orbit of some stabilizer. Two give overlap 96 iff their H4 subspaces
share a 3D intersection.

KEY FACTS TO USE:
- 40320 permutations / 40 decompositions = 1008 orderings per decomposition
- 124 pairs at overlap 96, each with unique 3D subspace
- Profile: (1,1,1,1/a1) for ALL 124 pairs
- Degree distribution: {5:4, 6:24, 7:12}

STRATEGY:
1. Characterize WHEN two Coxeter elements give overlap 96
2. Count using the Weyl group orbit structure
3. Look for a formula in terms of dim(E8), rank, h, a1
"""

import numpy as np
from itertools import permutations, product as iprod
from collections import Counter, defaultdict

phi = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
h = 30
rank = 8

print("=" * 70)
print("EXP-502: ALGEBRAIC PROOF THAT n96 = dim(E8)/2")
print("=" * 70)

# === Build E8 (verified code) ===
E8_roots = []
for i in range(8):
    for j in range(i+1, 8):
        for si in [1, -1]:
            for sj in [1, -1]:
                v = np.zeros(8); v[i] = si; v[j] = sj
                E8_roots.append(v)
for signs in iprod([1, -1], repeat=8):
    if sum(1 for s in signs if s == -1) % 2 == 0:
        E8_roots.append(np.array(signs) / 2.0)
E8 = np.array(E8_roots)

simple = np.zeros((8, 8))
for i in range(6):
    simple[i, i] = 1; simple[i, i+1] = -1
simple[6] = np.array([0,0,0,0,0,1,1,0])
simple[7] = np.array([-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,0.5])

def make_refl(a):
    return np.eye(8) - np.outer(a, a)

# === Map each permutation to its decomposition index ===
print("Building permutation -> decomposition mapping...")

def get_decomp_key(perm):
    C = np.eye(8)
    for i in perm:
        C = make_refl(simple[i]) @ C
    if not np.allclose(np.linalg.matrix_power(C, 30), np.eye(8), atol=1e-8):
        return None
    eigvals, eigvecs = np.linalg.eig(C)
    args = np.angle(eigvals)
    h4_idx = []
    for i in range(8):
        m = round(args[i] * 30 / (2*np.pi))
        if m <= 0: m += 30
        if m in {1, 11, 19, 29}:
            h4_idx.append(i)
    if len(h4_idx) != 4:
        return None
    # Build real basis
    raw = []
    used = set()
    for j in h4_idx:
        if j in used: continue
        v = eigvecs[:, j]
        for k in h4_idx:
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
    P = Q[:, :4]
    if P.shape != (8, 4):
        return None
    proj = E8 @ P
    norms2 = np.sum(proj**2, axis=1)
    inner = frozenset(i for i in range(240) if norms2[i] < 1.0)
    if len(inner) != 120:
        return None
    outer = frozenset(range(240)) - inner
    key = (inner, outer) if min(inner) < min(outer) else (outer, inner)
    return key, P

# Collect ALL decompositions with their permutations
decomp_perms = defaultdict(list)  # decomp_key -> list of permutations
decomp_projectors = {}  # decomp_key -> projector

for perm in permutations(range(8)):
    result = get_decomp_key(perm)
    if result is not None:
        key, P = result
        decomp_perms[key].append(perm)
        if key not in decomp_projectors:
            decomp_projectors[key] = P

n_decomp = len(decomp_perms)
print(f"Found {n_decomp} decompositions")
assert n_decomp == 40

keys = list(decomp_perms.keys())

# === Analyze the Coxeter elements themselves ===
print(f"\n{'='*70}")
print("STEP 1: COXETER ELEMENT STRUCTURE")
print(f"{'='*70}")

# For each decomposition, compute the Coxeter element matrix
coxeter_matrices = []
for key in keys:
    perm = decomp_perms[key][0]  # take first permutation
    C = np.eye(8)
    for i in perm:
        C = make_refl(simple[i]) @ C
    coxeter_matrices.append(C)

# How many DISTINCT Coxeter matrices are there?
distinct_cox = []
cox_to_decomp = {}
for idx, C in enumerate(coxeter_matrices):
    found = False
    for d_idx, C_ref in enumerate(distinct_cox):
        if np.allclose(C, C_ref, atol=1e-8):
            cox_to_decomp.setdefault(d_idx, []).append(idx)
            found = True
            break
    if not found:
        cox_to_decomp[len(distinct_cox)] = [idx]
        distinct_cox.append(C)

print(f"Distinct Coxeter matrices: {len(distinct_cox)}")
print(f"Decompositions per Coxeter matrix: {[len(v) for v in cox_to_decomp.values()]}")

# Key: multiple decompositions can share the SAME Coxeter element
# (different orderings of the same reflections can give the same matrix)
# OR different Coxeter elements can give the same H4 subspace
# (the H4 subspace depends on the eigenspaces, which are shared by conjugate elements)

# === Analyze overlap 96 pairs in terms of Coxeter elements ===
print(f"\n{'='*70}")
print("STEP 2: OVERLAP 96 IN TERMS OF COXETER ELEMENTS")
print(f"{'='*70}")

# Find all 124 pairs at overlap 96
pairs_96 = []
for i in range(40):
    for j in range(i+1, 40):
        inner_i = keys[i][0]
        ov = max(len(inner_i & keys[j][0]), len(inner_i & keys[j][1]))
        if ov == 96:
            pairs_96.append((i, j))

print(f"Pairs at overlap 96: {len(pairs_96)}")

# For each pair: how are their Coxeter elements related?
print(f"\nRelationship between Coxeter elements in overlap-96 pairs:")

# Check: C_i * C_j^{-1} or other relations
relation_types = Counter()
for i, j in pairs_96[:20]:  # sample first 20
    Ci = coxeter_matrices[i]
    Cj = coxeter_matrices[j]

    # Product
    prod = Ci @ np.linalg.inv(Cj)
    prod_order = None
    P = np.eye(8)
    for n in range(1, 61):
        P = P @ prod
        if np.allclose(P, np.eye(8), atol=1e-6):
            prod_order = n
            break

    # Trace of product (= character)
    tr_prod = np.trace(prod)

    # Commutator
    comm = Ci @ Cj - Cj @ Ci
    comm_norm = np.linalg.norm(comm)

    relation_types[prod_order] = relation_types.get(prod_order, 0) + 1

print(f"Orders of C_i * C_j^(-1) for overlap-96 pairs: {dict(relation_types)}")

# === What PERMUTATION relates two decompositions at overlap 96? ===
print(f"\n{'='*70}")
print("STEP 3: PERMUTATION STRUCTURE OF OVERLAP-96 PAIRS")
print(f"{'='*70}")

# For each pair (i,j) at overlap 96:
# Take one permutation from decomposition i and one from j.
# What is the "distance" between these permutations?
# (e.g., number of transpositions, cycle structure of p_i * p_j^{-1})

def perm_distance(p1, p2):
    """Compute the cycle structure of p1 * p2^{-1} in S_8."""
    # p2^{-1}
    inv_p2 = [0]*8
    for i, v in enumerate(p2):
        inv_p2[v] = i
    # composition p1 * p2^{-1}
    comp = tuple(p1[inv_p2[i]] for i in range(8))
    # cycle structure
    visited = set()
    cycles = []
    for i in range(8):
        if i in visited:
            continue
        cycle = []
        j = i
        while j not in visited:
            visited.add(j)
            cycle.append(j)
            j = comp[j]
        if len(cycle) > 0:
            cycles.append(len(cycle))
    return tuple(sorted(cycles, reverse=True))

cycle_structures = Counter()
for i, j in pairs_96:
    p_i = decomp_perms[keys[i]][0]
    p_j = decomp_perms[keys[j]][0]
    cs = perm_distance(p_i, p_j)
    cycle_structures[cs] += 1

print(f"Cycle structures of p_i * p_j^(-1) for overlap-96 pairs:")
for cs, count in sorted(cycle_structures.items(), key=lambda x: -x[1]):
    print(f"  {cs}: {count} pairs")

# === The H4 projectors: what's their algebraic structure? ===
print(f"\n{'='*70}")
print("STEP 4: PROJECTOR ALGEBRA")
print(f"{'='*70}")

# Each H4 projector P_i is a rank-4 projection matrix in R^8.
# P_i = Q_i @ Q_i^T where Q_i is 8x4 orthonormal.
# The overlap is Tr(P_i @ P_j).

projector_mats = []
for idx in range(40):
    Q = decomp_projectors[keys[idx]]
    P = Q @ Q.T  # 8x8 rank-4 projector
    projector_mats.append(P)

# Verify: Tr(P_i) = 4 for all i
traces = [np.trace(P) for P in projector_mats]
print(f"Traces of projectors: {sorted(set(np.round(traces, 6)))}")

# Tr(P_i @ P_j) at overlap 96
print(f"\nTr(P_i @ P_j) for overlap-96 pairs:")
tr_values = []
for i, j in pairs_96[:10]:
    tr = np.trace(projector_mats[i] @ projector_mats[j])
    tr_values.append(tr)
    if len(tr_values) <= 5:
        print(f"  Pair ({i},{j}): Tr(Pi*Pj) = {tr:.6f}")

print(f"  All values: {sorted(set(np.round(tr_values, 6)))}")
print(f"  Expected: 1+1+1+1/a1 = {3 + 1/a1:.6f}")

# Verify for ALL 124 pairs
all_tr_96 = []
for i, j in pairs_96:
    tr = np.trace(projector_mats[i] @ projector_mats[j])
    all_tr_96.append(tr)

print(f"  All 124 values: min={min(all_tr_96):.6f}, max={max(all_tr_96):.6f}")
print(f"  All equal to 3+1/a1 = 3.2: {all(abs(t - 3.2) < 0.001 for t in all_tr_96)}")

# === The KEY algebraic quantity: Tr(P_i @ P_j) ===
print(f"\n{'='*70}")
print("STEP 5: Tr(P_i @ P_j) FOR ALL 780 PAIRS")
print(f"{'='*70}")

# Tr(P_i @ P_j) = sum cos^2(theta_k) for the 4 principal angles
# This is the "Hilbert-Schmidt inner product" of the projectors.
# It takes values 2m/a1 for m = 5,6,7,8 (from the overlap quantization).

all_traces = defaultdict(int)
for i in range(40):
    for j in range(i+1, 40):
        tr = np.trace(projector_mats[i] @ projector_mats[j])
        tr_rounded = round(tr * a1) / a1  # round to multiples of 1/a1
        all_traces[tr_rounded] += 1

print(f"Tr(P_i @ P_j) distribution:")
for tr, count in sorted(all_traces.items()):
    m = tr * a1 / 2
    ov = int(h * tr)
    print(f"  Tr = {tr:.4f} = 2*{m:.0f}/a1, overlap = {ov}, count = {count}")

# === Sum of ALL Tr(P_i @ P_j) ===
print(f"\n{'='*70}")
print("STEP 6: TOTAL TRACE SUM")
print(f"{'='*70}")

total_trace = 0
for i in range(40):
    for j in range(i+1, 40):
        total_trace += np.trace(projector_mats[i] @ projector_mats[j])

print(f"Sum of Tr(P_i @ P_j) over all 780 pairs: {total_trace:.4f}")

# Also: sum of P_i over all 40 decompositions
P_sum = sum(projector_mats)
print(f"\nSum of all 40 projectors:")
print(f"  Trace: {np.trace(P_sum):.4f} = 40 * 4 = {40*4}")
print(f"  Is proportional to identity? ")
eigenvals_Psum = np.sort(np.linalg.eigvalsh(P_sum))[::-1]
print(f"  Eigenvalues: {eigenvals_Psum.round(4)}")
is_proportional = np.allclose(eigenvals_Psum, eigenvals_Psum[0], atol=0.01)
print(f"  Proportional to I: {is_proportional}")
if is_proportional:
    print(f"  P_sum = {eigenvals_Psum[0]:.4f} * I_8")
    print(f"  {eigenvals_Psum[0]:.4f} = 40*4/8 = {40*4/8}")
    print(f"  = n_decomp * sub_dim / rank = {40*4/8}")

# If P_sum = (40*4/8) * I = 20 * I, then:
# sum_{i<j} Tr(P_i P_j) = (1/2)(Tr(P_sum^2) - sum_i Tr(P_i^2))
# = (1/2)(Tr((20I)^2) - 40*Tr(P_i^2))
# = (1/2)(400*8 - 40*4) = (1/2)(3200 - 160) = 1520

print(f"\n  Verification: sum Tr(PiPj) from P_sum:")
sum_sq = np.trace(P_sum @ P_sum)
sum_diag = sum(np.trace(P @ P) for P in projector_mats)
computed = (sum_sq - sum_diag) / 2
print(f"  (Tr(P_sum^2) - sum Tr(Pi^2)) / 2 = ({sum_sq:.4f} - {sum_diag:.4f}) / 2 = {computed:.4f}")
print(f"  Direct sum: {total_trace:.4f}")
print(f"  Match: {abs(computed - total_trace) < 0.01}")

# === The identity we want to prove ===
print(f"\n{'='*70}")
print("STEP 7: THE IDENTITY")
print(f"{'='*70}")

# We know:
# P_sum = (n_decomp * sub_dim / rank) * I = 20 * I
# This means: the 40 H4 subspaces form a "tight frame" in R^8.
# (They're uniformly distributed on the Grassmannian in some sense.)

# From the frame condition:
# sum_i P_i = 20 * I
# sum_i Tr(P_i^2) = 40 * 4 = 160 (since each P_i has rank 4)
# sum_{i<j} Tr(P_i P_j) = (Tr((20I)^2) - 160) / 2 = (3200 - 160) / 2 = 1520

# Now: 1520 = 176*(10/5) + 296*(12/5) + 184*(14/5) + 124*(16/5)
# = 176*2 + 296*2.4 + 184*2.8 + 124*3.2
# = 352 + 710.4 + 515.2 + 396.8 = wait let me recompute

check = 176*2.0 + 296*2.4 + 184*2.8 + 124*3.2
print(f"\n  Weighted sum check: {check}")
print(f"  Frame identity: {(3200-160)/2}")
print(f"  Match: {abs(check - 1520) < 0.01}")

# So: 176*2 + 296*2.4 + 184*2.8 + 124*3.2 = 1520
# In terms of m (sum cos^2 = 2m/a1):
# 176*(10/5) + 296*(12/5) + 184*(14/5) + 124*(16/5)
# = (1/5) * (176*10 + 296*12 + 184*14 + 124*16)
# = (1/5) * (1760 + 3552 + 2576 + 1984) = (1/5) * 9872 = ???

raw = 176*10 + 296*12 + 184*14 + 124*16
print(f"\n  176*10 + 296*12 + 184*14 + 124*16 = {raw}")
print(f"  /5 = {raw/5}")
print(f"  Should equal 1520*5/... hmm. 1520 = {raw}/5? {raw/5 == 1520}")

# So the constraint is:
# sum_m n(m) * 2m/a1 = 1520
# i.e., sum_m n(m) * m = 1520 * a1 / 2 = 1520 * 5 / 2 = 3800

sum_nm = 176*5 + 296*6 + 184*7 + 124*8
print(f"\n  sum n(m)*m = {sum_nm}")
print(f"  1520*a1/2 = {1520*a1/2}")

# And sum_m n(m) = 780
sum_n = 176 + 296 + 184 + 124
print(f"  sum n(m) = {sum_n}")

# Two equations, four unknowns: need more constraints.
# Third constraint: sum_m n(m) * m^2 = ???

# Tr(P_sum^2 * P_sum) = Tr((20I)^3) = 20^3 * 8 = 64000? No.
# Actually: sum_{i,j,k distinct} Tr(P_i P_j P_k) involves triple products.

# Let's compute sum_m n(m) * (2m/a1)^2:
sum_tr_sq = 0
for i in range(40):
    for j in range(i+1, 40):
        tr = np.trace(projector_mats[i] @ projector_mats[j])
        sum_tr_sq += tr**2

print(f"\n  sum Tr(PiPj)^2 = {sum_tr_sq:.4f}")

check2 = 176*2.0**2 + 296*2.4**2 + 184*2.8**2 + 124*3.2**2
print(f"  176*4 + 296*5.76 + 184*7.84 + 124*10.24 = {check2:.4f}")
print(f"  Match: {abs(sum_tr_sq - check2) < 0.01}")

# So we have:
# sum n(m) = 780
# sum n(m) * (2m/a1) = 1520
# sum n(m) * (2m/a1)^2 = {sum_tr_sq}

# From these 3 constraints on 4 unknowns, we can determine 3 of the 4.
# The 4th is fixed by... the tight frame condition? Or something else?

print(f"\n{'='*70}")
print("STEP 8: CONSTRAINT SYSTEM")
print(f"{'='*70}")

# Let x = (n5, n6, n7, n8) = (176, 296, 184, 124)
# Constraints:
# (1) n5 + n6 + n7 + n8 = C(40,2) = 780
# (2) n5*10 + n6*12 + n7*14 + n8*16 = 1520*5 = 7600
#     (from tight frame: sum Tr(PiPj) = 1520, and Tr = 2m/5)
# (3) sum n_m * (2m/5)^2 = sum_tr_sq
# (4) ???

# Simplify constraint (2): 10n5 + 12n6 + 14n7 + 16n8 = 7600
# From (1): n5 = 780 - n6 - n7 - n8
# Sub into (2): 10(780-n6-n7-n8) + 12n6 + 14n7 + 16n8 = 7600
#               7800 - 10n6 - 10n7 - 10n8 + 12n6 + 14n7 + 16n8 = 7600
#               2n6 + 4n7 + 6n8 = -200
# Hmm, that gives negative. Let me recheck.

# Constraint (2): sum n(m) * (2m/a1) = 1520
# 176*(10/5) + 296*(12/5) + 184*(14/5) + 124*(16/5) = 1520
# 352 + 710.4 + 515.2 + 396.8 = 1974.4 ???

actual_sum = 176*2.0 + 296*2.4 + 184*2.8 + 124*3.2
print(f"\n  Actual sum Tr(PiPj) = {actual_sum:.4f}")
print(f"  Frame prediction = {(3200-160)/2:.4f}")
# These should match!

# Wait: 176*2 + 296*2.4 + 184*2.8 + 124*3.2
# = 352 + 710.4 + 515.2 + 396.8 = 1974.4
# But frame gives 1520. MISMATCH!

print(f"\n  MISMATCH: sum = {actual_sum:.1f} but frame gives {(3200-160)/2:.1f}")
print(f"  Direct computation: {total_trace:.4f}")

# So sum Tr(PiPj) = 1974.4 (from overlap distribution)
# But P_sum might NOT be proportional to I!

# Recheck
print(f"\n  P_sum eigenvalues: {eigenvals_Psum.round(4)}")
print(f"  Is P_sum = 20*I? {is_proportional}")

# If P_sum != 20*I, then the frame is NOT tight, and the constraint is different.

if not is_proportional:
    print(f"\n  P_sum is NOT proportional to identity!")
    print(f"  The 40 H4 subspaces do NOT form a tight frame.")
    print(f"  This means the 'frame constraint' doesn't apply directly.")
    print(f"  The constraint system has fewer equations than assumed.")

# === Even without tight frame, what constraints DO we have? ===
print(f"\n{'='*70}")
print("STEP 9: ACTUAL CONSTRAINTS")
print(f"{'='*70}")

# The ONLY algebraic constraints are:
# (1) sum n(m) = C(40, 2) = 780 (total pairs)
# (2) overlap values are in {60, 72, 84, 96} (from Grassmannian quantization)
# (3) P_sum has specific eigenvalues (not proportional to I)

# Without additional constraints, the system is UNDERDETERMINED.
# We have 4 unknowns and only 1 equation from pair counting.

# The Grassmannian quantization (fact that overlaps are in {60,72,84,96})
# is a DEEP constraint from the Weyl group, but it doesn't determine the counts.

# CONCLUSION: A purely ALGEBRAIC proof of n_96 = 124 requires
# understanding the Weyl group orbit structure at a level beyond
# simple counting arguments.

print(f"""
  P_sum eigenvalues: {eigenvals_Psum.round(4)}
  P_sum is {'NOT ' if not is_proportional else ''}proportional to identity.

  Constraints:
  (1) n5 + n6 + n7 + n8 = 780 (pair counting)
  (2) Overlaps in {{60, 72, 84, 96}} (Grassmannian quantization)
  (3) P_sum eigenspectrum: {eigenvals_Psum.round(2)}

  These do NOT uniquely determine n8 = 124.
  A proof requires Weyl group orbit analysis.
""")

# === What DO the P_sum eigenvalues tell us? ===
print(f"{'='*70}")
print("STEP 10: P_SUM EIGENSTRUCTURE")
print(f"{'='*70}")

# P_sum = sum of 40 rank-4 projectors in R^8
# Tr(P_sum) = 40*4 = 160
# If P_sum were 20*I, eigenvalues would be (20,20,20,20,20,20,20,20)
# Actual eigenvalues tell us about the non-uniformity

print(f"P_sum = sum of 40 projectors:")
print(f"  Eigenvalues: {eigenvals_Psum.round(6)}")
print(f"  Trace: {sum(eigenvals_Psum):.1f} = 40*4 = 160")
print(f"  If uniform: all eigenvalues = 160/8 = 20")

# Eigenvectors of P_sum: these are the "preferred directions" of the configuration
_, evecs_Psum = np.linalg.eigh(P_sum)
print(f"\n  The 40 H4 subspaces have a PREFERRED frame:")
for i in range(8):
    print(f"    Direction {i}: eigenvalue {eigenvals_Psum[7-i]:.4f}")

# Deviation from uniform
dev = eigenvals_Psum - 20
print(f"\n  Deviations from 20: {dev.round(4)}")
print(f"  Sum of deviations: {sum(dev):.6f} (should be 0)")
print(f"  Sum of squared deviations: {sum(dev**2):.4f}")

# === SUMMARY ===
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"""
ALGEBRAIC PROOF STATUS:
========================

1. The identity n_96 = 124 = dim(E8)/2 is NUMERICALLY EXACT.

2. The tight frame assumption (P_sum = 20*I) is {'VALID' if is_proportional else 'INVALID'}.
   P_sum eigenvalues: {eigenvals_Psum.round(2)}

3. Simple counting constraints (pair totals, overlap quantization)
   are INSUFFICIENT to determine n_96 uniquely.

4. A proof requires deeper Weyl group analysis:
   - Count W(E8) orbits on pairs of H4 subspaces
   - Or use the specific structure of Coxeter elements
   - Or find an independent algebraic identity

5. The fact that n_96 = dim(E8)/2 = (|roots| + rank)/2 = 124
   WHILE being the count of codimension-1 Grassmannian neighbors
   suggests a connection to the Lie algebra structure,
   but the MECHANISM is not yet identified.

WHAT A PROOF WOULD NEED:
  Show that the number of pairs of Coxeter-related H4 subspaces
  in Gr(4,8) with Tr(P_i P_j) = 16/5 equals exactly
  (|roots(E8)| + rank(E8)) / 2.

  This is likely a theorem in the representation theory of
  the Weyl group W(E8) = Aut(E8 root system).
""")

print("=" * 70)
print("EXP-502 COMPLETE")
print("=" * 70)
