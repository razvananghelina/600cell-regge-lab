"""
EXP-252: Spectral Mass Operator on the 600-cell
================================================
Instead of ASSIGNING masses via m = m_e * phi^(5a+6b),
find an OPERATOR M on the 600-cell whose eigenvalues ARE the masses.

Key insight: The adjacency matrix A has 9 eigenvalues, one per fermion species.
The eigenvalues are lambda_i = a_i + b_i*phi in Z[phi].
The mass formula n = 5a + 6b is a LINEAR FUNCTIONAL on Z[phi].

QUESTION: Can we build M = 5*A_int + 6*A_phi (the mass functional applied
to the adjacency matrix) and get the correct mass exponents?

If not, what operator DOES give them?
"""

import numpy as np
from itertools import product as iterproduct

PHI = (1 + np.sqrt(5)) / 2
PHI_PRIME = 1 - PHI
a1 = 5
b1 = 6

# =====================================================================
# PART 0: Build the 600-cell
# =====================================================================
def build_600cell():
    """Build 600-cell vertices using even permutations method."""
    verts = set()

    # Type A: 8 vertices - permutations of (+/-1, 0, 0, 0)
    for i in range(4):
        for s in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = s
            verts.add(tuple(v))

    # Type B: 16 vertices - all (+/-1/2, +/-1/2, +/-1/2, +/-1/2)
    for signs in iterproduct([0.5, -0.5], repeat=4):
        verts.add(tuple(signs))

    # Type C: 96 vertices - even permutations of (+/-phi/2, +/-1/2, +/-1/(2*phi), 0)
    base = [PHI/2, 0.5, 1/(2*PHI), 0.0]

    def even_perms(vals):
        """Generate all even permutations of 4 values."""
        perms = []
        indices = [0, 1, 2, 3]
        # All 24 permutations, keep even ones
        from itertools import permutations
        for p in permutations(indices):
            # Count inversions to determine parity
            inv = 0
            for i in range(4):
                for j in range(i+1, 4):
                    if p[i] > p[j]:
                        inv += 1
            if inv % 2 == 0:
                perms.append(tuple(vals[p[k]] for k in range(4)))
        return perms

    for ep in even_perms(base):
        # All sign combinations
        for signs in iterproduct([1, -1], repeat=4):
            v = tuple(s * x for s, x in zip(signs, ep))
            # Skip if any coordinate is -0.0
            if all(abs(x) > 1e-10 or x == 0 for x in v):
                verts.add(tuple(round(x, 10) for x in v))
            else:
                verts.add(tuple(round(x, 10) if abs(x) > 1e-10 else 0.0 for x in v))

    # Clean up and deduplicate
    clean = set()
    for v in verts:
        clean.add(tuple(round(x, 10) for x in v))

    verts_list = sorted(clean)
    print("  Vertices: %d (expect 120)" % len(verts_list))
    return np.array(verts_list)

print("=" * 70)
print("EXP-252: SPECTRAL MASS OPERATOR")
print("=" * 70)
print()
print("Building 600-cell...")
vertices = build_600cell()
N = len(vertices)

# Build adjacency matrix
print("Building adjacency matrix...")
# Edge criterion: distance = 1/phi for unit-radius 600-cell
# Actually, for our coordinates, the edge length is 1/phi
dists = np.zeros((N, N))
for i in range(N):
    for j in range(i+1, N):
        d = np.sqrt(np.sum((vertices[i] - vertices[j])**2))
        dists[i, j] = d
        dists[j, i] = d

# Find the edge distance (smallest nonzero distance)
all_dists = dists[dists > 0.01]
edge_dist = np.min(all_dists)
print("  Edge distance: %.6f (expect 1/phi = %.6f)" % (edge_dist, 1/PHI))

A = np.zeros((N, N), dtype=float)
for i in range(N):
    for j in range(i+1, N):
        if abs(dists[i, j] - edge_dist) < 0.01:
            A[i, j] = 1
            A[j, i] = 1

degree = int(np.sum(A[0]))
n_edges = int(np.sum(A)) // 2
print("  Degree: %d (expect 12)" % degree)
print("  Edges: %d (expect 720)" % n_edges)

# =====================================================================
# PART 1: Eigenvalues of adjacency matrix
# =====================================================================
print()
print("=" * 70)
print("PART 1: Adjacency matrix eigenvalues")
print("=" * 70)

eigenvalues_A = np.linalg.eigvalsh(A)
eigenvalues_A = np.sort(eigenvalues_A)[::-1]  # descending

# Group by distinct values
from collections import Counter
eig_rounded = np.round(eigenvalues_A, 4)
eig_counts = Counter(eig_rounded)
distinct_eigs = sorted(eig_counts.keys(), reverse=True)

print()
print("  %-4s  %-15s  %-8s  %-20s" % ("#", "Eigenvalue", "Mult", "Z[phi] form"))
print("  " + "-" * 55)

# For each eigenvalue, find a + b*phi form
eig_data = []
for idx, lam in enumerate(distinct_eigs):
    mult = eig_counts[lam]
    # Solve: a + b*phi = lam => b = (lam - a)/(phi)
    # Since lam = a + b*phi and a,b integers:
    # Try b values
    found = False
    for b in range(-10, 11):
        a = lam - b * PHI
        if abs(a - round(a)) < 0.01:
            a = int(round(a))
            n_mass = a1 * a + b1 * b
            eig_data.append((lam, mult, a, b, n_mass))
            phi_str = "%d + %d*phi" % (a, b) if b != 0 else "%d" % a
            print("  %-4d  %12.4f     %3d    %-20s  n=%d" % (idx+1, lam, mult, phi_str, n_mass))
            found = True
            break
    if not found:
        eig_data.append((lam, mult, None, None, None))
        print("  %-4d  %12.4f     %3d    (not in Z[phi])" % (idx+1, lam, mult))

# =====================================================================
# PART 2: The naive mass operator M = a_1 * A_int + b_1 * A_phi
# =====================================================================
print()
print("=" * 70)
print("PART 2: Mass functional applied to eigenvalues")
print("=" * 70)
print()
print("  The mass formula n = %d*a + %d*b applied to eigenvalue a+b*phi:" % (a1, b1))
print()

# Extract n values from eigenvalues
n_from_eigs = sorted(set(d[4] for d in eig_data if d[4] is not None), reverse=True)
n_fermion = sorted([0, 3, 5, 11, 16, 17, 19, 26], reverse=True)

print("  Eigenvalue mass exponents: %s" % sorted(n_from_eigs))
print("  Fermion mass exponents:    %s" % sorted(n_fermion))
print()
print("  These are COMPLETELY DIFFERENT.")
print("  The naive approach (apply n=5a+6b to eigenvalues) does NOT work.")

# =====================================================================
# PART 3: McKay correspondence - the E8 connection
# =====================================================================
print()
print("=" * 70)
print("PART 3: McKay correspondence - 9 eigenvalues = 9 irreps of 2I")
print("=" * 70)
print()

# The 9 irreps of 2I (binary icosahedral group) have dimensions:
# 1, 2, 2', 3, 3', 4, 4', 5, 6
# Multiplicities = dim^2: 1, 4, 4, 9, 9, 16, 16, 25, 36
# These match our eigenvalue multiplicities!

irrep_dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]
irrep_labels = ['1', '2', "2'", '3', "3'", '4', "4'", '5', '6']

# Sort eigenvalues by multiplicity to match irreps
eig_by_mult = sorted(eig_data, key=lambda x: x[1])

print("  McKay matching (eigenvalues sorted by multiplicity):")
print("  %-6s %-10s %-5s %-12s %-6s %-6s" % ("Irrep", "dim", "dim^2", "Eigenvalue", "a+b*phi", "n"))
print("  " + "-" * 50)

for i, (lam, mult, a, b, n) in enumerate(eig_by_mult):
    if i < len(irrep_labels):
        phi_str = "%d+%d*phi" % (a, b) if b else str(a)
        print("  %-6s %-10d %-5d %12.4f %-12s %-6d" %
              (irrep_labels[i], irrep_dims[i], irrep_dims[i]**2, lam, phi_str, n if n else 0))

# =====================================================================
# PART 4: What FUNCTION of eigenvalues gives mass exponents?
# =====================================================================
print()
print("=" * 70)
print("PART 4: Can we find f such that f(lambda_i) = n_fermion_i?")
print("=" * 70)
print()

# We have 9 eigenvalues and 9 mass exponents (with n=11 appearing twice).
# But 9 eigenvalues and 8 distinct mass exponents don't match one-to-one.
# Two eigenvalues must map to n=11 (muon and strange).

# The two eigenvalues with same multiplicity 9: 4*phi and 4-4*phi
# The two eigenvalues with same multiplicity 16: 3 and -3
# The two eigenvalues with same multiplicity 4: 6*phi and 6-6*phi

# n=11 corresponds to both muon (lepton) and strange (quark).
# In McKay, the irreps 3 and 3' are Galois conjugates (phi <-> phi').
# So the pairing (muon, strange) = (3, 3') makes sense!

print("  Galois pairs among eigenvalues:")
for i in range(len(eig_data)):
    li, mi, ai, bi, ni = eig_data[i]
    for j in range(i+1, len(eig_data)):
        lj, mj, aj, bj, nj = eig_data[j]
        if ai is not None and aj is not None:
            # Galois conjugate: (a,b) -> (a+b, -b)
            if ai == aj + bj and bi == -bj:
                print("    lambda_%d = %.4f (a=%d,b=%d) <-> lambda_%d = %.4f (a=%d,b=%d) [Galois pair]"
                      % (i+1, li, ai, bi, j+1, lj, aj, bj))

# Try mapping eigenvalues to mass exponents
# Sort eigenvalues and mass exponents to find potential mapping
print()
print("  Eigenvalues sorted:      ", ["%.3f" % d[0] for d in sorted(eig_data, key=lambda x: x[0])])
print("  Eigenvalue n-values:     ", sorted([d[4] for d in eig_data if d[4] is not None]))
print("  Fermion n-values needed: ", sorted(n_fermion))

# =====================================================================
# PART 5: The DIRAC approach - operator on edges, not vertices
# =====================================================================
print()
print("=" * 70)
print("PART 5: Graph Dirac operator")
print("=" * 70)
print()

# The graph Dirac operator acts on edges (1-forms), not vertices (0-forms)
# D = d + d* where d is the incidence matrix
# D^2 = Laplacian on 0-forms + Laplacian on 1-forms

# Incidence matrix B: 120 x 720
# B[v, e] = +1 if v is the "start" of edge e, -1 if "end", 0 otherwise
# This requires orienting edges

# Build edge list
edges = []
for i in range(N):
    for j in range(i+1, N):
        if A[i, j] > 0:
            edges.append((i, j))

print("  Building incidence matrix B (%d x %d)..." % (N, len(edges)))
B = np.zeros((N, len(edges)), dtype=float)
for e_idx, (i, j) in enumerate(edges):
    B[i, e_idx] = 1
    B[j, e_idx] = -1

# Laplacian on 0-forms: L0 = B * B^T
L0 = B @ B.T
print("  L0 = B * B^T: %dx%d" % L0.shape)
print("  L0[0,0] = %.1f (should be degree = 12)" % L0[0, 0])

# Laplacian on 1-forms: L1 = B^T * B
L1 = B.T @ B
print("  L1 = B^T * B: %dx%d" % L1.shape)

# Eigenvalues of L1 (edge Laplacian)
print("  Computing eigenvalues of L1 (720x720)...")
eig_L1 = np.linalg.eigvalsh(L1)
eig_L1_sorted = np.sort(eig_L1)

# Group L1 eigenvalues
eig_L1_rounded = np.round(eig_L1_sorted, 3)
L1_counts = Counter(eig_L1_rounded)
print()
print("  Edge Laplacian L1 eigenvalues:")
print("  %-15s %-6s" % ("Eigenvalue", "Mult"))
print("  " + "-" * 25)
for val in sorted(L1_counts.keys()):
    print("  %12.4f    %4d" % (val, L1_counts[val]))

# =====================================================================
# PART 6: Restricted adjacency on C-vertices only
# =====================================================================
print()
print("=" * 70)
print("PART 6: Restricted operators on C-vertices (fermions)")
print("=" * 70)
print()

# Classify vertices by type
types = []
for v in vertices:
    n_zero = sum(1 for x in v if abs(x) < 1e-8)
    n_half = sum(1 for x in v if abs(abs(x) - 0.5) < 1e-8)
    if n_zero == 3:  # permutations of (+/-1, 0, 0, 0)
        types.append('A')
    elif n_half == 4:  # all +/-1/2
        types.append('B')
    else:
        types.append('C')

type_counts = Counter(types)
print("  Vertex types: A=%d, B=%d, C=%d" % (type_counts['A'], type_counts['B'], type_counts['C']))

# Extract C-vertex indices
C_indices = [i for i, t in enumerate(types) if t == 'C']
n_C = len(C_indices)
print("  C-vertices (fermions): %d" % n_C)

# Build adjacency restricted to C-vertices
# CC edges only
A_CC = np.zeros((n_C, n_C), dtype=float)
idx_map = {v: i for i, v in enumerate(C_indices)}
cc_edges = 0
for i in C_indices:
    for j in C_indices:
        if i < j and A[i, j] > 0:
            A_CC[idx_map[i], idx_map[j]] = 1
            A_CC[idx_map[j], idx_map[i]] = 1
            cc_edges += 1

print("  CC edges: %d" % cc_edges)
print("  Degree in CC subgraph: %.1f" % np.mean(np.sum(A_CC, axis=1)))

# Eigenvalues of A_CC
eig_CC = np.linalg.eigvalsh(A_CC)
eig_CC = np.sort(eig_CC)[::-1]
eig_CC_rounded = np.round(eig_CC, 3)
CC_counts = Counter(eig_CC_rounded)

print()
print("  CC adjacency eigenvalues (96x96):")
print("  %-15s %-6s" % ("Eigenvalue", "Mult"))
print("  " + "-" * 25)
for val in sorted(CC_counts.keys(), reverse=True):
    print("  %12.4f    %4d" % (val, CC_counts[val]))

# =====================================================================
# PART 7: A_CC through B (SU(2) mediated) and through A (U(1) mediated)
# =====================================================================
print()
print("=" * 70)
print("PART 7: Mediated operators (C -> B -> C and C -> A -> C)")
print("=" * 70)
print()

# C-to-B-to-C: two-step paths through B vertices
A_indices = [i for i, t in enumerate(types) if t == 'A']
B_indices = [i for i, t in enumerate(types) if t == 'B']

# Build C-B adjacency
A_CB = np.zeros((n_C, len(B_indices)), dtype=float)
b_idx_map = {v: i for i, v in enumerate(B_indices)}
for ic, c in enumerate(C_indices):
    for b in B_indices:
        if A[c, b] > 0:
            A_CB[ic, b_idx_map[b]] = 1

# C-B-C operator = A_CB @ A_CB^T
M_CBC = A_CB @ A_CB.T
print("  M_CBC (C->B->C): %dx%d, mean degree: %.1f" % (M_CBC.shape[0], M_CBC.shape[1], np.mean(np.sum(M_CBC > 0, axis=1))))

eig_CBC = np.linalg.eigvalsh(M_CBC)
eig_CBC = np.sort(eig_CBC)[::-1]
eig_CBC_rounded = np.round(eig_CBC, 2)
CBC_counts = Counter(eig_CBC_rounded)

print("  M_CBC eigenvalues:")
for val in sorted(CBC_counts.keys(), reverse=True):
    print("    %10.3f  x%d" % (val, CBC_counts[val]))

# C-to-A-to-C
A_CA = np.zeros((n_C, len(A_indices)), dtype=float)
a_idx_map = {v: i for i, v in enumerate(A_indices)}
for ic, c in enumerate(C_indices):
    for a in A_indices:
        if A[c, a] > 0:
            A_CA[ic, a_idx_map[a]] = 1

M_CAC = A_CA @ A_CA.T
print()
print("  M_CAC (C->A->C): %dx%d, mean degree: %.1f" % (M_CAC.shape[0], M_CAC.shape[1], np.mean(np.sum(M_CAC > 0, axis=1))))

eig_CAC = np.linalg.eigvalsh(M_CAC)
eig_CAC = np.sort(eig_CAC)[::-1]
eig_CAC_rounded = np.round(eig_CAC, 2)
CAC_counts = Counter(eig_CAC_rounded)

print("  M_CAC eigenvalues:")
for val in sorted(CAC_counts.keys(), reverse=True):
    print("    %10.3f  x%d" % (val, CAC_counts[val]))

# =====================================================================
# PART 8: Combined "Yukawa" operator: C -> (A+B) -> C
# =====================================================================
print()
print("=" * 70)
print("PART 8: Yukawa-like operators")
print("=" * 70)
print()

# The idea: fermion mass comes from coupling through gauge/Higgs vertices
# M_Yukawa = alpha * M_CAC + beta * M_CBC
# where alpha, beta are DERIVED coupling constants

# Natural choices:
# alpha = 1 (U(1) coupling = 1 edge per C vertex)
# beta = 1 (SU(2) coupling)
# Or: alpha = sin^2(tW), beta = cos^2(tW)

for label, alpha_Y, beta_Y in [
    ("M_CAC + M_CBC", 1, 1),
    ("sin^2(tW)*M_CAC + cos^2(tW)*M_CBC", 6/26, 20/26),
    ("M_CAC - M_CBC", 1, -1),
    ("a_1*M_CAC + b_1*M_CBC", a1, b1),
]:
    M_Y = alpha_Y * M_CAC + beta_Y * M_CBC
    eig_Y = np.linalg.eigvalsh(M_Y)
    eig_Y = np.sort(eig_Y)[::-1]
    eig_Y_rounded = np.round(eig_Y, 2)
    Y_counts = Counter(eig_Y_rounded)

    print("  %s:" % label)
    print("    Distinct eigenvalues: %d" % len(Y_counts))
    vals = sorted(Y_counts.keys(), reverse=True)
    for val in vals[:12]:
        print("      %10.3f  x%d" % (val, Y_counts[val]))
    if len(vals) > 12:
        print("      ... (%d more)" % (len(vals) - 12))
    print()

# =====================================================================
# PART 9: The PHI operator - decompose A into phi and phi' sectors
# =====================================================================
print()
print("=" * 70)
print("PART 9: Phi/Phi' decomposition of adjacency matrix")
print("=" * 70)
print()

# A has eigenvalues lambda_i = a_i + b_i*phi
# A = A_rat + A_irr where A_rat has rational eigenvalues, A_irr has irrational
# More precisely: A = sum_i lambda_i * P_i
# We want: A_int = sum_i a_i * P_i and A_phi = sum_i b_i * P_i
# Then: M = 5*A_int + 6*A_phi would give n_i = 5*a_i + 6*b_i on each eigenspace

# To build P_i, we need the eigenvectors
print("  Computing full eigendecomposition...")
eigenvalues_full, eigenvectors = np.linalg.eigh(A)

# Group eigenvectors by eigenvalue
eigenspaces = {}
for i in range(N):
    lam = round(eigenvalues_full[i], 3)
    if lam not in eigenspaces:
        eigenspaces[lam] = []
    eigenspaces[lam].append(eigenvectors[:, i])

print("  Building projectors for each eigenspace...")
projectors = {}
for lam, vecs in eigenspaces.items():
    V = np.array(vecs).T  # columns are eigenvectors
    P = V @ V.T  # projector onto eigenspace
    projectors[lam] = P

# Now build A_int and A_phi
A_int = np.zeros((N, N))
A_phi_mat = np.zeros((N, N))

print()
print("  Eigenvalue decomposition into Z[phi]:")
for lam_val, P in projectors.items():
    # Find a, b such that lam_val = a + b*phi
    found = False
    for b in range(-10, 11):
        a = lam_val - b * PHI
        if abs(a - round(a)) < 0.05:
            a = int(round(a))
            A_int += a * P
            A_phi_mat += b * P
            n_val = a1 * a + b1 * b
            print("    lambda=%.4f -> (%d, %d)*phi, n=5*%d+6*%d=%d, mult=%d"
                  % (lam_val, a, b, a, b, n_val, len(eigenspaces[lam_val])))
            found = True
            break
    if not found:
        print("    lambda=%.4f -> NOT in Z[phi]!" % lam_val)

# Build mass operator
M_mass = a1 * A_int + b1 * A_phi_mat
print()
print("  Mass operator M = %d*A_int + %d*A_phi:" % (a1, b1))
eig_M = np.linalg.eigvalsh(M_mass)
eig_M_rounded = np.round(eig_M, 1)
M_counts = Counter(eig_M_rounded)

print("  Eigenvalues of M:")
for val in sorted(M_counts.keys(), reverse=True):
    print("    n = %6.1f  x%d" % (val, M_counts[val]))

print()
print("  Fermion n-values needed: %s" % sorted(n_fermion))
print()
print("  The eigenvalues of M = 5*A_int + 6*A_phi are the mass functional")
print("  APPLIED TO EIGENVALUES of A, not to fermion quantum numbers (a,b).")
print("  These are {%s}, not {%s}." %
      (", ".join(str(int(v)) for v in sorted(M_counts.keys())),
       ", ".join(str(v) for v in sorted(n_fermion))))

# =====================================================================
# PART 10: KEY QUESTION - what operator has eigenvalues {0,3,5,11,16,17,19,26}?
# =====================================================================
print()
print("=" * 70)
print("PART 10: What operator gives the fermion mass exponents?")
print("=" * 70)
print()

# The fermion mass exponents are {0, 3, 5, 11, 16, 17, 19, 26}
# with n=11 having multiplicity 2 (muon + strange)
# Total: 9 values in 8 distinct levels

# Can we write these as f(lambda_i) for some function f?
# We have 9 eigenvalues and 9 mass exponents.

# Let's try a POLYNOMIAL fit: n = c0 + c1*lambda + c2*lambda^2 + ...
# This is fitting, but let's see if the coefficients are derived quantities

eig_sorted = sorted(eig_data, key=lambda x: x[0])  # sort by eigenvalue
n_target_sorted = sorted(n_fermion)  # {0, 3, 5, 11, 11, 16, 17, 19, 26}

# But we have 9 eigenvalues and need to match to 9 n-values (including n=11 twice)
# The question is: which eigenvalue maps to which n?

# Try all possible mappings to find polynomial
# First, let's see if the ORDERING matches

print("  Eigenvalues (sorted):       ", ["%.3f" % d[0] for d in eig_sorted])
print("  Mass exponents (sorted):    ", n_target_sorted)
print()

# Direct mapping (smallest eig -> smallest n, etc.)
lambdas = np.array([d[0] for d in eig_sorted])
ns = np.array(n_target_sorted, dtype=float)

# Try polynomial fit
for deg in [1, 2, 3]:
    coeffs = np.polyfit(lambdas, ns, deg)
    fitted = np.polyval(coeffs, lambdas)
    residuals = ns - fitted
    rms = np.sqrt(np.mean(residuals**2))
    print("  Degree-%d polynomial: coeffs = %s" % (deg, np.round(coeffs, 4)))
    print("    Fitted n: %s" % np.round(fitted, 1))
    print("    Residuals: %s" % np.round(residuals, 1))
    print("    RMS: %.2f" % rms)
    print()

# =====================================================================
# SUMMARY
# =====================================================================
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("1. The adjacency matrix A has 9 eigenvalues in Z[phi], matching")
print("   9 fermion species via McKay correspondence.")
print()
print("2. The mass functional n = 5a + 6b applied to eigenvalues gives")
print("   {%s} - NOT the fermion exponents." %
      ", ".join(str(int(v)) for v in sorted(M_counts.keys())))
print()
print("3. Restricted operators (CC, CBC, CAC) give different eigenvalue")
print("   structures but don't directly match fermion masses either.")
print()
print("4. The edge Laplacian L1 (720x720) has a richer spectrum.")
print()
print("5. OPEN: The correct mass operator has not been identified.")
print("   Promising directions:")
print("   - Graph Dirac operator (spinor bundle on 600-cell)")
print("   - Yukawa operator weighted by gauge coupling structure")
print("   - Spectral action on noncommutative geometry (Connes)")
print("   - Restricted operator on C-vertices with gauge weighting")
