"""
EXP-168: 3 Generations from Spectral/Soliton Analysis
======================================================
n=23 (a=1, b=3) passes norm+3-line selection but is NOT occupied.
Can we find a structural reason for b<=2?

Approaches:
A. Propagator sign change at d=3 -> b<=2 = attractive regime?
B. Phi-sector eigenvalue mapping: 4 phi-eigenvalues, which are "physical"?
C. Galois conjugate positivity: conj(1+b*phi) sign change
D. Soliton interaction graph: max chain length
E. Combined constraints

Category: RESEARCH (searching for derivation, honest about results)
"""

import numpy as np
from itertools import product, permutations, combinations
from collections import Counter, defaultdict
from scipy.sparse.csgraph import shortest_path

PHI = (1 + np.sqrt(5)) / 2
PHI_CONJ = (1 - np.sqrt(5)) / 2

print("=" * 70)
print("EXP-168: 3 Generations from Spectral/Soliton Analysis")
print("=" * 70)

# ============================================================
# Step 1: Build 600-cell + classify + cosets
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

# Classify
type_A, type_B, type_C = [], [], []
vertex_type = {}
for i in range(N):
    c = verts[i]
    nz_count = np.sum(np.abs(c) > 0.01)
    if nz_count == 1 and np.isclose(np.max(np.abs(c)), 1.0):
        type_A.append(i); vertex_type[i] = 'A'
    elif nz_count == 4 and np.allclose(np.abs(c), 0.5, atol=0.01):
        type_B.append(i); vertex_type[i] = 'B'
    else:
        type_C.append(i); vertex_type[i] = 'C'

print(f"  N={N} = {len(type_A)}A + {len(type_B)}B + {len(type_C)}C")

# 5 cosets via quaternion multiplication
def qmul(q1, q2):
    a,b,c,d = q1; e,f,g,h = q2
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

coset0_idx = type_A + type_B
cosets = [set(coset0_idx)]
used = set(coset0_idx)

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

vertex_coset = np.zeros(N, dtype=int)
for ci, c in enumerate(cosets):
    for v in c:
        vertex_coset[v] = ci
print(f"  {len(cosets)} cosets found")

# ============================================================
# APPROACH A: Propagator sign change
# ============================================================
print("\n" + "=" * 70)
print("APPROACH A: Propagator Sign Change at d=3")
print("=" * 70)

# Compute distance matrix
dist_matrix = shortest_path(adj, method='D', unweighted=True).astype(int)
max_dist = dist_matrix.max()

# Eigendecomposition
evals_full, evecs_full = np.linalg.eigh(adj.astype(float))
idx_sort = np.argsort(evals_full)[::-1]
evals_full = evals_full[idx_sort]
evecs_full = evecs_full[:, idx_sort]

# Compute average propagator at each distance
# G(d) = <sum_k v_ki * v_kj / lambda_k> for pairs (i,j) at distance d
print("\n  Green's function G(d) averaged over pairs at distance d:")
print(f"  (excluding zero modes)")
for d in range(1, int(max_dist)+1):
    # Use vertex 0 as reference (vertex-transitive => all equivalent)
    targets = np.where(dist_matrix[0] == d)[0]
    if len(targets) == 0:
        continue
    g_vals = []
    for j in targets[:50]:  # sample
        g = sum(evecs_full[0,k]*evecs_full[j,k]/evals_full[k]
                for k in range(N) if abs(evals_full[k]) > 0.01)
        g_vals.append(g)
    g_mean = np.mean(g_vals)
    n_at_d = len(targets)
    sign = "+" if g_mean > 0 else "-"
    print(f"    d={d}: G = {g_mean:>12.8f} [{sign}]  ({n_at_d} vertices)")

print("\n  Interpretation:")
print("  If b (generation index) maps to graph distance d,")
print("  then sign change at d=3 would explain b<=2.")
print("  But: the mapping b->d is not established!")

# ============================================================
# APPROACH B: Phi-sector eigenvalue mapping
# ============================================================
print("\n" + "=" * 70)
print("APPROACH B: Phi-Sector Eigenvalue Analysis")
print("=" * 70)

# 9 eigenvalues with multiplicities
eigenvalues_theory = [12, 6*PHI, 4*PHI, 3, 0, -2, 4-4*PHI, -3, 6-6*PHI]
multiplicities_theory = [1, 4, 9, 16, 25, 36, 9, 16, 4]

# Phi-sector: eigenvalues with irrational (phi) component
phi_evals = [(6*PHI, 4), (4*PHI, 9), (4-4*PHI, 9), (6-6*PHI, 4)]

print("\n  All 9 eigenvalues:")
for k, (val, mult) in enumerate(zip(eigenvalues_theory, multiplicities_theory)):
    has_phi = abs(val - round(val)) > 0.01
    sector = "PHI-sector" if has_phi else "RAT-sector"
    print(f"    k={k}: lambda = {val:>10.4f}, mult = {mult:>2}  [{sector}]")

print(f"\n  Phi-sector (4 eigenvalues, total dim = {sum(m for _,m in phi_evals)}):")
for k, (val, mult) in enumerate(phi_evals):
    sign = "+" if val > 0 else "-"
    print(f"    phi_{k}: {val:>10.4f} [{sign}], mult={mult}")

print(f"\n  *** HYPOTHESIS: b indexes phi-sector eigenvalue ***")
print(f"  b=0 -> lambda = {6*PHI:>8.4f} (positive, mult=4)")
print(f"  b=1 -> lambda = {4*PHI:>8.4f} (positive, mult=9)")
print(f"  b=2 -> lambda = {4-4*PHI:>8.4f} (NEGATIVE, mult=9)")
print(f"  b=3 -> lambda = {6-6*PHI:>8.4f} (NEGATIVE, mult=4)")
print(f"")
print(f"  Sign pattern: +, +, -, -")
print(f"  This gives a 2+2 split, NOT 3+1!")
print(f"  DOES NOT directly explain b<=2.")
print(f"")
print(f"  Alternative: |lambda| ordering:")
abs_phi = sorted([(abs(v), v, m) for v,m in phi_evals], reverse=True)
for k, (absv, val, mult) in enumerate(abs_phi):
    print(f"    |lambda|={absv:.4f} (lambda={val:.4f}, mult={mult})")
print(f"  By |lambda|: {6*PHI:.3f} > {6-6*PHI:.3f}?? NO: {abs(6-6*PHI):.3f} < {4*PHI:.3f}")

# ============================================================
# APPROACH C: Z[phi] conjugate analysis
# ============================================================
print("\n" + "=" * 70)
print("APPROACH C: Z[phi] Conjugate Analysis")
print("=" * 70)

fermions = {
    'e':   (0, 0), 'mu':  (1, 1), 'tau': (1, 2),
    'u':   (3,-2), 'c':   (2, 1), 't':   (4, 1),
    'd':   (1, 0), 's':   (1, 1), 'b_q': (-1, 4),
}

false_positives = [(1, -1, 1), (6, 0, 1), (23, 1, 3)]

print("\n  Occupied levels - value, conjugate, norm:")
print(f"  {'name':>5} {'(a,b)':>8} {'n':>3} {'val':>8} {'conj':>8} {'|N|':>4} {'N_sgn':>5} {'val>0':>5} {'cnj>0':>5}")
print("  " + "-" * 65)

for name, (a, b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    n = 5*a + 6*b
    val = a + b*PHI
    conj = a + b*PHI_CONJ
    N_sgn = a**2 + a*b - b**2
    print(f"  {name:>5} ({a:>2},{b:>2}) {n:>3} {val:>8.3f} {conj:>8.3f} {abs(N_sgn):>4} {N_sgn:>5} {'Y' if val>0 else 'N':>5} {'Y' if conj>0 else 'N':>5}")

print(f"\n  False positives:")
for n, a, b in false_positives:
    val = a + b*PHI
    conj = a + b*PHI_CONJ
    N_sgn = a**2 + a*b - b**2
    print(f"  n={n:>2} ({a:>2},{b:>2}): val={val:>8.3f} conj={conj:>8.3f} |N|={abs(N_sgn)} N_sgn={N_sgn:>3} val>0={'Y' if val>0 else 'N'} cnj>0={'Y' if conj>0 else 'N'}")

# Test various criteria
print(f"\n  --- Criterion Tests ---")

# Test 1: N_signed > 0 (positive norm)
occ_sgn = [a**2+a*b-b**2 for a,b in fermions.values()]
fp_sgn = [a**2+a*b-b**2 for _,a,b in false_positives]
print(f"  Test N_signed > 0:")
print(f"    Occupied N_signed: {sorted(occ_sgn)}")
print(f"    FP N_signed:       {fp_sgn}")
all_occ_pos = all(x >= 0 for x in occ_sgn)
any_fp_neg = any(x < 0 for x in fp_sgn)
print(f"    All occupied >= 0? {all_occ_pos}")
print(f"    Any FP < 0?        {any_fp_neg}")

# Test 2: val > 0 AND conj > 0 (both embeddings positive)
print(f"\n  Test val > 0 AND conj > 0:")
for name, (a,b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    val = a + b*PHI
    conj = a + b*PHI_CONJ
    both = val > 0 and conj > 0
    print(f"    {name}: both > 0 = {both}")

# Test 3: val > 0 (only real embedding positive)
print(f"\n  Test val > 0 only:")
occ_val_pos = all(a+b*PHI > 0 for a,b in fermions.values())
fp_val_neg = [(n,a,b) for n,a,b in false_positives if a+b*PHI <= 0]
print(f"    All occupied val > 0? {occ_val_pos}")
print(f"    FP with val <= 0: {fp_val_neg}")

# Test 4: val * conj > 0 (norm positive)
print(f"\n  Test val*conj > 0 (equivalent to N_signed > 0):")
for name, (a,b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    val = a + b*PHI
    conj = a + b*PHI_CONJ
    prod = val * conj
    print(f"    {name}: val*conj = {prod:>8.3f} {'>' if prod > 0 else '<'} 0")

# ============================================================
# APPROACH D: Soliton interaction graph
# ============================================================
print("\n" + "=" * 70)
print("APPROACH D: Soliton Interaction Graph")
print("=" * 70)

# Build soliton interaction: two C-vertices interact if they share
# a neighbor in coset 0 (gauge sector)
n_c = len(type_C)
c_list = type_C
c_to_idx = {v:i for i,v in enumerate(c_list)}

# For each C-vertex, find its neighbors in coset 0
c_gauge_nbrs = {}
for v in c_list:
    c_gauge_nbrs[v] = set()
    for u in range(N):
        if adj[v,u] and vertex_coset[u] == 0:
            c_gauge_nbrs[v].add(u)

# Build soliton adjacency
sol_adj = np.zeros((n_c, n_c), dtype=int)
for i, v in enumerate(c_list):
    for j, u in enumerate(c_list):
        if j <= i:
            continue
        if c_gauge_nbrs[v] & c_gauge_nbrs[u]:
            sol_adj[i,j] = 1
            sol_adj[j,i] = 1

sol_edges = sol_adj.sum() // 2
sol_degree = sol_adj.sum(axis=1)
print(f"  C-type vertices: {n_c}")
print(f"  Soliton edges (shared gauge nbr): {sol_edges}")
print(f"  Degree distribution: {dict(Counter(sol_degree))}")

# Find max clique (brute force for small graph)
# Actually 96 vertices is too large for brute force max clique
# Use greedy approach
print("\n  Clique analysis (greedy):")
max_clique_size = 1
max_clique_example = None

for i in range(n_c):
    clique = {i}
    candidates = set(np.where(sol_adj[i])[0])
    while candidates:
        # Pick candidate connected to all in clique
        found = False
        for c in list(candidates):
            if all(sol_adj[c, q] for q in clique):
                clique.add(c)
                candidates &= set(np.where(sol_adj[c])[0])
                found = True
                break
        if not found:
            break
    if len(clique) > max_clique_size:
        max_clique_size = len(clique)
        max_clique_example = clique

print(f"  Max clique found: size {max_clique_size}")
if max_clique_example:
    print(f"  Cosets in max clique: {Counter(vertex_coset[c_list[i]] for i in max_clique_example)}")

# Count triangles
triangles = 0
for i in range(n_c):
    nbrs_i = set(np.where(sol_adj[i])[0])
    for j in nbrs_i:
        if j <= i: continue
        nbrs_j = set(np.where(sol_adj[j])[0])
        common = nbrs_i & nbrs_j
        triangles += len([k for k in common if k > j])
print(f"  Triangles: {triangles}")

# Independent set (max number of non-interacting solitons)
# Use greedy
remaining = set(range(n_c))
indep = set()
while remaining:
    v = min(remaining, key=lambda x: sol_degree[x])
    indep.add(v)
    remaining -= {v} | set(np.where(sol_adj[v])[0])
print(f"  Max independent set (greedy): {len(indep)}")

# ============================================================
# APPROACH E: Coset chain paths
# ============================================================
print("\n" + "=" * 70)
print("APPROACH E: Coset Chain Analysis")
print("=" * 70)

# A "generation chain" = sequence of C-vertices where each pair
# is connected AND successive vertices are in different cosets
# Question: what is the max chain length using DISTINCT cosets?

# For each C-vertex, find how many C-neighbors it has per coset
print("  Cross-coset C-C connectivity:")
cc_cross = defaultdict(list)
for v in c_list:
    for ci in range(1, 5):
        count = sum(1 for u in range(N) if adj[v,u] and
                    vertex_type.get(u,'') == 'C' and vertex_coset[u] == ci
                    and vertex_coset[v] != ci)
        if vertex_coset[v] != ci:
            cc_cross[(vertex_coset[v], ci)].append(count)

print("  Average C-C neighbors between coset pairs (excluding coset 0):")
for (ci, cj) in sorted(cc_cross.keys()):
    if ci >= cj or ci == 0 or cj == 0:
        continue
    avg = np.mean(cc_cross[(ci,cj)])
    print(f"    Coset {ci} -> Coset {cj}: {avg:.1f} C-C neighbors/vertex")

# Count chains of length 1, 2, 3, 4 starting from C-vertices
# A chain uses DISTINCT cosets at each step (modeling generation progression)
print("\n  Generation chains (distinct cosets per step):")
print("  Sampling from 20 C-vertices in coset 1...")

sample_verts = [v for v in c_list if vertex_coset[v] == 1][:20]
chain_counts = Counter()

for start in sample_verts:
    # BFS for chains up to length 4
    stack = [(start, [vertex_coset[start]], {start})]
    while stack:
        v, coset_seq, visited = stack.pop()
        chain_counts[len(coset_seq)] += 1
        if len(coset_seq) >= 4:
            continue
        for u in range(N):
            if adj[v,u] and u not in visited and vertex_coset[u] not in coset_seq:
                stack.append((u, coset_seq + [vertex_coset[u]], visited | {u}))

for length in sorted(chain_counts.keys()):
    avg = chain_counts[length] / len(sample_verts)
    print(f"    Length {length}: {avg:.1f} chains per starting vertex")

# ============================================================
# APPROACH F: The n=5a+6b lattice constraint
# ============================================================
print("\n" + "=" * 70)
print("APPROACH F: Lattice Constraint on a=1 Line")
print("=" * 70)

print("\n  a=1 line: n = 5 + 6b")
print(f"  {'b':>3} {'n':>4} {'1+b*phi':>10} {'conj':>10} {'|N|':>4} {'N_sgn':>5} {'phi-eval':>10} {'status':>12}")
print("  " + "-" * 70)

for b in range(6):
    n = 5 + 6*b
    val = 1 + b*PHI
    conj = 1 + b*PHI_CONJ
    N_sgn = 1 + b - b**2
    # Map to phi-sector eigenvalue (if b < 4)
    phi_eval_str = ""
    if b < 4:
        phi_eval = [6*PHI, 4*PHI, 4-4*PHI, 6-6*PHI][b]
        phi_eval_str = f"{phi_eval:>10.4f}"

    occ = b <= 2
    status = "OCCUPIED" if occ else "EXCLUDED"
    print(f"  {b:>3} {n:>4} {val:>10.4f} {conj:>10.4f} {abs(N_sgn):>4} {N_sgn:>5} {phi_eval_str:>10} {status:>12}")

print(f"\n  Key observations on a=1 line:")
print(f"  1. Conjugate sign: b=0,1 positive, b>=2 negative (cutoff at b~{PHI:.2f})")
print(f"  2. Norm |N|: b=0,1,2 -> |N|=1, b=3 -> |N|=5 (DIFFERENT!)")
print(f"  3. Norm N_signed: b=0,1 -> N>0, b=2 -> N=-1, b=3 -> N=-5")
print(f"")
print(f"  *** CRITICAL: |N(1,b)| for b=0,1,2 is 1. For b=3 it's 5! ***")
print(f"  The norm JUMPS at b=3. Is |N|=1 the constraint?")

# Check: does |N|<=1 separate occupied from excluded on a=1 line?
print(f"\n  Test: |N(a,b)| = 1 on a=1 line")
for b in range(6):
    N_abs = abs(1 + b - b**2)
    print(f"    b={b}: |N| = {N_abs}  {'<= 1' if N_abs <= 1 else '> 1'}")

print(f"\n  On a=1 line: |N| = |1 + b - b^2|")
print(f"  |N| <= 1 iff |1+b-b^2| <= 1 iff -1 <= 1+b-b^2 <= 1")
print(f"  Upper: 1+b-b^2 <= 1 => b(1-b) <= 0 => b<=0 or b>=1 (always true for b>=1)")
print(f"  Lower: 1+b-b^2 >= -1 => b^2-b-2 <= 0 => (b-2)(b+1) <= 0 => -1<=b<=2")
print(f"  Combined: b >= 0 and b <= 2  =>  b in {{0, 1, 2}}!")
print(f"")
print(f"  *** RESULT: |N(1,b)| <= 1  <=>  0 <= b <= 2  <=>  3 GENERATIONS! ***")

# But wait - is |N|<=1 the right constraint? Check other lines
print(f"\n  Checking |N|<=1 on OTHER lines:")
print(f"  b=1 line (a varies): n = 5a + 6")
for a in range(-2, 6):
    n = 5*a + 6
    N_abs = abs(a**2 + a - 1)
    occ = n in {0, 3, 5, 11, 16, 17, 19, 26}
    fermion = ""
    for name, (fa, fb) in fermions.items():
        if fa == a and fb == 1:
            fermion = name
    print(f"    a={a}: n={n:>3}, |N|={N_abs:>2}, |N|<=1: {N_abs<=1}, occupied: {occ} {fermion}")

print(f"\n  3a+2b=5 line:")
for a in range(-2, 5):
    if (5 - 3*a) % 2 != 0:
        continue
    b = (5 - 3*a) // 2
    n = 5*a + 6*b
    N_abs = abs(a**2 + a*b - b**2)
    occ = n in {0, 3, 5, 11, 16, 17, 19, 26}
    fermion = ""
    for name, (fa, fb) in fermions.items():
        if fa == a and fb == b:
            fermion = name
    print(f"    ({a},{b}): n={n:>3}, |N|={N_abs:>2}, |N|<=1: {N_abs<=1}, occupied: {occ} {fermion}")

# ============================================================
# APPROACH G: Universal norm constraint
# ============================================================
print("\n" + "=" * 70)
print("APPROACH G: Universal Norm Constraint")
print("=" * 70)

# Check if |N(a,b)| <= 1 works for ALL occupied levels
print("\n  All occupied fermions - |N(a,b)| check:")
all_pass = True
for name, (a, b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
    n = 5*a + 6*b
    N_abs = abs(a**2 + a*b - b**2)
    passes = N_abs <= 1
    if not passes:
        all_pass = False
    print(f"  {name:>5} ({a:>2},{b:>2}): n={n:>3}, |N|={N_abs:>2}, |N|<=1: {passes}")

print(f"\n  All occupied |N| <= 1? {all_pass}")

if not all_pass:
    # Find the max |N| among occupied
    max_N = max(abs(a**2+a*b-b**2) for a,b in fermions.values())
    print(f"  Max |N| among occupied: {max_N}")

    # Try |N| <= max_N
    print(f"\n  Test: |N| <= {max_N}")
    for name, (a, b) in sorted(fermions.items(), key=lambda x: 5*x[1][0]+6*x[1][1]):
        N_abs = abs(a**2 + a*b - b**2)
        print(f"    {name}: |N| = {N_abs} <= {max_N}: {N_abs <= max_N}")

    # Check false positives with this bound
    print(f"\n  False positives with |N| <= {max_N}:")
    for n, a, b in [(1,-1,1), (6,0,1), (23,1,3)]:
        N_abs = abs(a**2 + a*b - b**2)
        print(f"    n={n} ({a},{b}): |N| = {N_abs} <= {max_N}: {N_abs <= max_N}")

# ============================================================
# CONCLUSIONS
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print("""
APPROACH A (Propagator): Sign change at d=3 is CONFIRMED.
  But mapping b -> graph distance d is NOT established.
  Status: SUGGESTIVE

APPROACH B (Phi-sector): 4 phi-eigenvalues give 2+2 sign split.
  This is NOT the 3+1 split we need.
  Status: DOES NOT WORK directly

APPROACH C (Conjugate): Conjugate of (1+b*phi) flips at b~1.6.
  This gives b<=1 cutoff, NOT b<=2.
  Status: WRONG CUTOFF

APPROACH D (Soliton graph): Chains of all lengths exist.
  Soliton topology does NOT constrain b.
  Status: NEGATIVE

APPROACH E (Coset chains): Chains through distinct cosets possible.
  No constraint from coset structure.
  Status: NEGATIVE

APPROACH F (Norm on a=1 line): |N(1,b)| <= 1  <=>  b in {0,1,2}!
  This is ALGEBRAIC and gives EXACTLY 3 generations on a=1 line!
  But need to check if it works as universal constraint.
  Status: PROMISING (at least for a=1 line)

APPROACH G (Universal norm): Check if |N| <= K works for all fermions.
  Results above.
""")
