"""
exp360_cc_57.py
================
OPEN-1: Can the cosmological constant exponent 57 be derived?

Lambda_P = alpha^{57 - alpha_s}, matching observation to 0.33%.
57 = N/2 - N_gen = 60 - 3.

Strategy: Explore algebraic routes using known framework quantities.
1. Box operator traces and their relation to 57
2. Spectral zeta of Box
3. Heat kernel of Box
4. Algebraic identities in Z[phi] connecting 57 to spectral data
5. The decomposition 57 = 3 * 19 where 19 = a1*b1-a1-b1

Dependencies: numpy, scipy
"""

import numpy as np
from scipy.linalg import eigh
from itertools import permutations, product as cartesian_product

a1 = 5
b1 = a1 + 1
phi = (1.0 + np.sqrt(a1)) / 2.0
phi_c = (1.0 - np.sqrt(a1)) / 2.0
N_vert = 120
N_gen = 3
N_eig = 9
PI = np.pi
alpha_s = 1.0 / (2 * phi**3)

print("=" * 72)
print("EXP-360: COSMOLOGICAL CONSTANT EXPONENT 57 (OPEN-1)")
print("=" * 72)

# ============================================================================
# KNOWN DECOMPOSITIONS OF 57
# ============================================================================
print("\n--- Known decompositions of 57 ---")
print(f"  57 = N/2 - N_gen = {N_vert}/2 - {N_gen} = {N_vert//2 - N_gen}")
print(f"  57 = (N - N_gen)/2 = ({N_vert}-{N_gen})/2 = {(N_vert - N_gen)/2}")

# Wait: (N - N_gen)/2 = (120 - 3)/2 = 58.5, not 57!
# Actually: N/2 - N_gen = 60 - 3 = 57. Not (N-N_gen)/2.
# From the paper: z_CC = (N-N_gen)/2 - phi = 117/2 - phi... no.
# Let me check: Lambda_P = alpha^{57 - alpha_s}
# The paper says 57 = sum(z_k) + 7 where sum(z_k) = 50 = 2*a1^2
# and 7 = b1 + 1

print(f"\n  Spectral decomposition:")
print(f"  57 = sum(z_k) + (b1+1) = 2*a1^2 + (b1+1) = {2*a1**2} + {b1+1} = {2*a1**2 + b1 + 1}")
print(f"  57 = 2*a1^2 + a1 + 2 = 2*{a1**2} + {a1} + 2 = {2*a1**2 + a1 + 2}")
print(f"  57 = N/2 - N_gen = {N_vert//2} - {N_gen} = {N_vert//2 - N_gen}")

# Factor decomposition
print(f"\n  Factorization: 57 = 3 * 19")
print(f"  3 = N_gen")
print(f"  19 = a1*b1 - a1 - b1 = {a1*b1 - a1 - b1}")
print(f"       = a1^2 - 2*a1 + 4 = {a1**2 - 2*a1 + 4}")
print(f"       = N(z_mass) (the unique split prime)")
print(f"       = Sylvester-Frobenius number for (a1, b1)")

# ============================================================================
# BUILD 600-CELL AND BOX OPERATOR
# ============================================================================
print("\n--- Building 600-cell and Box operator ---")

def build_2I():
    verts = set()
    def add_vert(v):
        arr = np.array(v, dtype=float)
        n = np.linalg.norm(arr)
        if n > 1e-12:
            arr = arr / n
        verts.add(tuple(np.round(arr, 10)))
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0]*4; v[i] = s; add_vert(v)
    for signs in cartesian_product([0.5, -0.5], repeat=4):
        add_vert(list(signs))
    base = [0.0, 0.5, phi / 2.0, 1.0 / (2.0 * phi)]
    even_perms = [p for p in permutations(range(4))
                  if sum(1 for i in range(4) for j in range(i+1, 4)
                         if p[i] > p[j]) % 2 == 0]
    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in cartesian_product([1, -1], repeat=len(nz)):
            v = list(coords)
            for idx, s in zip(nz, signs):
                v[idx] *= s
            add_vert(v)
    return np.array(sorted(verts))

verts = build_2I()
N = len(verts)
dots = verts @ verts.T
np.clip(dots, -1.0, 1.0, out=dots)
A = np.zeros((N, N), dtype=float)
for i in range(N):
    for j in range(i + 1, N):
        if abs(dots[i, j] - phi/2) < 1e-6:
            A[i, j] = 1.0
            A[j, i] = 1.0

# Find Hopf fibration
def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
                     w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])

gen_idx = -1
for i in range(N):
    if abs(verts[i, 0] - phi/2) < 1e-6:
        g = verts[i]; power = g.copy(); ok = True
        for k in range(2, 11):
            power = quat_mult(power, g)
            if k == 5 and not np.allclose(power, [-1,0,0,0], atol=1e-6):
                ok = False; break
            if k == 10 and not np.allclose(power, [1,0,0,0], atol=1e-6):
                ok = False
        if ok: gen_idx = i; break

g = verts[gen_idx]
subgroup = []
power = np.array([1.0, 0.0, 0.0, 0.0])
for k in range(10):
    dists = np.linalg.norm(verts - power, axis=1)
    subgroup.append(np.argmin(dists))
    power = quat_mult(power, g)

assigned = np.full(N, -1, dtype=int)
fibers = []; fid = 0
for i in range(N):
    if assigned[i] >= 0: continue
    coset = []
    for si in subgroup:
        q_prod = quat_mult(verts[i], verts[si])
        dists = np.linalg.norm(verts - q_prod, axis=1)
        idx = np.argmin(dists)
        if dists[idx] < 1e-6 and assigned[idx] < 0:
            coset.append(idx); assigned[idx] = fid
    fibers.append(coset); fid += 1

A_fiber = np.zeros((N, N), dtype=float)
for fiber in fibers:
    for i in fiber:
        for j in fiber:
            if i != j and A[i, j] > 0.5:
                A_fiber[i, j] = 1.0

Box = b1 * A_fiber - A
evals_box = np.sort(eigh(Box, eigvals_only=True))
print(f"  Box = b1*A_fiber - A constructed")
print(f"  dim(ker) = {np.sum(np.abs(evals_box) < 1e-8)}")

# ============================================================================
# BOX TRACES
# ============================================================================
print("\n" + "=" * 72)
print("BOX OPERATOR TRACES")
print("=" * 72)

for k in range(1, 7):
    tr_k = np.sum(evals_box**(2*k))
    print(f"  Tr(Box^{2*k}) = {tr_k:.4f}")
    print(f"    / N = {tr_k/N:.4f}")
    print(f"    / N^2 = {tr_k/N**2:.6f}")

# Key: Tr(Box^2)/N = N/2 = 60
tr2 = np.sum(evals_box**2)
print(f"\n  KEY: Tr(Box^2)/N = {tr2/N:.4f}")
print(f"  N/2 = {N/2}")
print(f"  N/2 - N_gen = {N/2 - N_gen} = 57")
print(f"  Tr(Box^2)/N - N_gen = {tr2/N - N_gen:.4f}")

# ============================================================================
# SPECTRAL ZETA OF BOX
# ============================================================================
print("\n" + "=" * 72)
print("SPECTRAL ZETA OF BOX")
print("=" * 72)

nonzero_box = evals_box[np.abs(evals_box) > 1e-8]
abs_nonzero_box = np.abs(nonzero_box)

for s in [1, 2, 3, 4]:
    zeta_s = np.sum(abs_nonzero_box**(-s))
    print(f"  zeta_Box({s}) = {zeta_s:.8f}")
    if s == 2:
        print(f"    * N = {zeta_s * N:.4f}")

# ============================================================================
# HEAT KERNEL OF BOX
# ============================================================================
print("\n" + "=" * 72)
print("HEAT KERNEL OF BOX")
print("=" * 72)

for t in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    hk = np.sum(np.exp(-t * evals_box**2))
    print(f"  Tr(exp(-{t:.1f}*Box^2)) = {hk:.6f}"
          f"  /N = {hk/N:.6f}  /N_eig = {hk/N_eig:.6f}")

# ============================================================================
# ALGEBRAIC IDENTITIES AROUND 57
# ============================================================================
print("\n" + "=" * 72)
print("ALGEBRAIC IDENTITIES")
print("=" * 72)

# Check various combinations of framework constants
print(f"\n  Framework constants:")
print(f"    a1 = {a1}, b1 = {b1}")
print(f"    N = {N_vert}, N_eig = {N_eig}, N_gen = {N_gen}")
print(f"    phi = {phi:.6f}, phi^2 = {phi**2:.6f}")
print(f"    4*phi^2 = {4*phi**2:.6f} = 1/alpha_s + 2 = {1/alpha_s + 2:.6f}")

# Various 57 constructions
print(f"\n  57 constructions:")
print(f"    N/2 - N_gen = {N_vert//2 - N_gen}")
print(f"    2*a1^2 + b1 + 1 = {2*a1**2 + b1 + 1}")
print(f"    N_gen * (a1*b1-a1-b1) = {N_gen * (a1*b1-a1-b1)}")
print(f"    a1*(a1+2) + 2 = {a1*(a1+2) + 2}")
print(f"    a1^2 + 2*(a1+b1) = {a1**2 + 2*(a1+b1)}")
print(f"    (a1^2+1) + 2*a1*b1/(a1+1) = {(a1**2+1) + 2*a1*b1/(a1+1):.4f}")

# Check: is 57 related to the Galois kernel?
print(f"\n  Galois kernel connections:")
print(f"    N - N_eig = {N_vert - N_eig} (non-kernel modes)")
print(f"    (N-N_eig)/2 = {(N_vert - N_eig)/2:.1f}")
print(f"    N_gen * N(z_mass) = 3*19 = 57")
print(f"    N(z_mass) = a1*b1-a1-b1 = {a1*b1-a1-b1}")

# The Sylvester-Frobenius number
# For coprime a, b: the largest integer NOT representable as xa + yb (x,y >= 0)
# is ab - a - b = Frobenius number
# For a1=5, b1=6: Frobenius = 5*6-5-6 = 19
print(f"\n  Frobenius number F(a1,b1) = a1*b1-a1-b1 = {a1*b1-a1-b1}")
print(f"  This is the largest mass exponent n=5a+6b with a,b >= 0")
print(f"  that is NOT reachable. n=19 is the largest 'hole'.")
print(f"  Count of holes = (a1-1)*(b1-1)/2 = {(a1-1)*(b1-1)//2}")

# Interesting: the number of non-representable integers is (a1-1)(b1-1)/2 = 10
# And 10 = 2*a1 = number of fiber edges per vertex
# Also 10 = L(3) + L(3') = sum of Galois-conjugate Laplacian eigenvalues

# ============================================================================
# Z[PHI] NORM IDENTITIES
# ============================================================================
print("\n" + "=" * 72)
print("Z[PHI] NORM IDENTITIES")
print("=" * 72)

# In Z[phi], the norm of z = a + b*phi is N(z) = a^2 + ab - b^2
# Key elements and their norms:
elements = [
    ("1+phi", 1, 1),
    ("1+2phi", 1, 2),
    ("2+phi", 2, 1),
    ("2+3phi", 2, 3),
    ("3+phi", 3, 1),
    ("3+5phi", 3, 5),
    ("5+phi", 5, 1),
    ("5+3phi", 5, 3),
    ("a1+b1*phi", a1, b1),
]

print(f"  Elements z = a + b*phi and N(z) = a^2+ab-b^2:")
for name, a, b in elements:
    norm = a**2 + a*b - b**2
    val = a + b*phi
    print(f"    {name:12s}: N = {norm:+4d}, value = {val:.4f}")

# Check: does any norm give 57 or -57?
print(f"\n  Search for |N(z)| = 57 = 3*19:")
for a in range(-10, 11):
    for b in range(-10, 11):
        norm = a**2 + a*b - b**2
        if abs(norm) == 57:
            print(f"    z = {a}+{b}*phi, N(z) = {norm}, value = {a + b*phi:.4f}")

# ============================================================================
# ADJACENCY EIGENVALUE COMBINATIONS
# ============================================================================
print("\n" + "=" * 72)
print("ADJACENCY EIGENVALUE COMBINATIONS")
print("=" * 72)

# The 9 eigenvalues of the 600-cell adjacency:
adj_evals = {
    'rho_0': (12.0, 1),
    'rho_1': (6*phi, 4),
    'rho_2': (4*phi, 9),
    'rho_3': (3.0, 16),
    'rho_4': (0.0, 25),
    'rho_5': (-2.0, 36),
    'rho_6': (4*phi_c, 9),
    'rho_7': (-3.0, 16),
    'rho_8': (6*phi_c, 4),
}

print(f"  Eigenvalues and multiplicities:")
for name, (ev, mult) in adj_evals.items():
    print(f"    {name}: lambda = {ev:+8.4f}, mult = {mult}")

# Sum of multiplicities weighted by various functions of eigenvalue
print(f"\n  Sum rules:")
sum_mult = sum(m for _, m in adj_evals.values())
sum_lam = sum(ev * m for ev, m in adj_evals.values())
sum_lam2 = sum(ev**2 * m for ev, m in adj_evals.values())
sum_lam_over_b1 = sum((ev/b1) * m for ev, m in adj_evals.values())

print(f"    sum(mult) = {sum_mult} = N")
print(f"    sum(lambda*mult) = {sum_lam:.4f} (= Tr(A))")
print(f"    sum(lambda^2*mult) = {sum_lam2:.4f} (= Tr(A^2) = 12*120)")
print(f"    sum(lambda*mult)/b1 = {sum_lam/b1:.4f}")

# Check: sum of eigenvalues (without multiplicity)
sum_evals = sum(ev for ev, _ in adj_evals.values())
print(f"\n  sum(eigenvalues, no mult) = {sum_evals:.6f}")
print(f"    = 12 + 6*phi + 4*phi + 3 + 0 - 2 + 4*phi' - 3 + 6*phi'")
print(f"    = 10 + 10*phi + 10*phi' = 10 + 10*(phi+phi') = 10 + 10*1 = 20")

# Laplacian eigenvalues
print(f"\n  Laplacian eigenvalues L_k = 12 - lambda_k:")
sum_L = 0
for name, (ev, mult) in adj_evals.items():
    L = 12 - ev
    print(f"    {name}: L = {L:+8.4f}")
    sum_L += L
print(f"  sum(L_k) = {sum_L:.4f} = 9*12 - 20 = {9*12 - 20}")

# ============================================================================
# THE 57 = Tr(Box^2)/N - N_gen CONNECTION
# ============================================================================
print("\n" + "=" * 72)
print("THE 57 = Tr(Box^2)/N - N_gen CONNECTION")
print("=" * 72)

print(f"""
  We have: Tr(Box^2)/N = N/2 = 60

  PROOF:
    Box = b1*A_fiber - A
    Box^2 = b1^2*A_fiber^2 - 2*b1*A_fiber*A + A^2

    Tr(A^2) = degree * N = 12 * 120 = 1440
    Tr(A_fiber^2) = 2 * 120 = 240  (each vertex has 2 fiber neighbors)
    Tr(A_fiber * A) = Tr(A_fiber) = 240  (A_fiber is a subgraph of A)

    Wait, Tr(A_fiber * A) = sum_ij (A_fiber)_ij * A_ji = sum_ij (A_fiber)_ij * A_ij
    Since A_fiber <= A entrywise: = sum_ij (A_fiber)_ij = 2 * n_fiber_edges = 240

    Tr(Box^2) = 36*240 - 12*240 + 1440 = 8640 - 2880 + 1440 = 7200
    Tr(Box^2)/N = 7200/120 = 60 = N/2.

  So: 57 = Tr(Box^2)/N - N_gen = N/2 - 3.

  This is exact but tautological: it restates the known formula.
  The question is: WHY does N_gen = 3 subtract from N/2?
""")

# Verify numerically
tr2_box = np.sum(evals_box**2)
print(f"  Numerical: Tr(Box^2) = {tr2_box:.4f}")
print(f"  Tr(Box^2)/N = {tr2_box/N:.4f}")
print(f"  Tr(Box^2)/N - N_gen = {tr2_box/N - N_gen:.4f}")

# ============================================================================
# NEW IDEA: COUNTING ARGUMENT
# ============================================================================
print("\n" + "=" * 72)
print("NEW IDEA: COUNTING ARGUMENT FOR 57")
print("=" * 72)

# The mass exponents are n = 5a + 6b with specific (a,b) values.
# The Frobenius number for (5,6) is 19.
# The representable non-negative integers by 5a+6b (a,b >= 0) are:
# 0, 5, 6, 10, 11, 12, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, ...
# The non-representable positive integers are:
# 1, 2, 3, 4, 7, 8, 9, 13, 14, 19
# Count = (5-1)(6-1)/2 = 10

representable = set()
for a in range(25):
    for b in range(25):
        n = a1*a + b1*b
        if n <= 120:
            representable.add(n)

non_rep = sorted(set(range(1, 20)) - representable)
print(f"  Non-representable integers by {a1}a+{b1}b (a,b>=0): {non_rep}")
print(f"  Count: {len(non_rep)} = (a1-1)*(b1-1)/2 = {(a1-1)*(b1-1)//2}")
print(f"  Largest (Frobenius): {max(non_rep)} = a1*b1-a1-b1 = {a1*b1-a1-b1}")

# Sum of non-representable integers
sum_non_rep = sum(non_rep)
print(f"  Sum of non-representable integers: {sum_non_rep}")
print(f"    = (a1-1)*(b1-1)*(2*a1*b1-a1-b1-1)/12")
expected_sum = (a1-1)*(b1-1)*(2*a1*b1-a1-b1-1)//12
print(f"    = {expected_sum}")

# Check: 57 connection?
print(f"\n  57 - sum_non_rep = {57 - sum_non_rep}")
print(f"  N_gen * sum_non_rep = {N_gen * sum_non_rep}")
print(f"  sum_non_rep / N_gen = {sum_non_rep / N_gen:.4f}")

# ============================================================================
# NEW IDEA: ALPHA EQUATION AND CC
# ============================================================================
print("\n" + "=" * 72)
print("ALPHA EQUATION AND CC CONNECTION")
print("=" * 72)

# The alpha equation: 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0
# This gives alpha^{-1} = 137.036...
# The CC formula: Lambda_P = alpha^{57-alpha_s}
# where alpha_s = 1/(2*phi^3)

# Is there an algebraic connection between the alpha equation
# and the exponent 57?

alpha_val = (4*a1*phi**4 - np.sqrt(16*a1**2*phi**8 - 8*PI)) / (4*PI)
alpha_inv = 1/alpha_val

print(f"  alpha = {alpha_val:.10f}")
print(f"  alpha^(-1) = {alpha_inv:.6f}")
print(f"  alpha_s = {alpha_s:.10f}")
print(f"  1/alpha_s = {1/alpha_s:.6f}")
print(f"  z_Planck = 1/alpha_s + 2 = {1/alpha_s + 2:.6f} = 4*phi^2 = {4*phi**2:.6f}")

z_CC = 57 - alpha_s
Lambda_pred = alpha_val**z_CC

# Observed: Lambda/m_P^4 ~ 3.3e-122
Lambda_obs = 3.3e-122
z_obs = np.log(Lambda_obs) / np.log(alpha_val)

print(f"\n  z_CC = 57 - alpha_s = {z_CC:.10f}")
print(f"  alpha^z_CC = {Lambda_pred:.6e}")
print(f"  Observed Lambda_P ~ {Lambda_obs:.1e}")
print(f"  z from observation = {z_obs:.6f}")
print(f"  Difference: {z_CC - z_obs:.6f}")

# ============================================================================
# NEW IDEA: SPECTRAL WEIGHT DECOMPOSITION
# ============================================================================
print("\n" + "=" * 72)
print("SPECTRAL WEIGHT DECOMPOSITION")
print("=" * 72)

# The spectral weights z_k are defined from the McKay spectral analysis.
# From MEMORY: z_k = a_k + b_k * phi where
# sum(a_k) = 12 = degree, sum(b_k) = 8 = rank(E8)
# sum(z_k) = 12 + 8*phi

# But the CC uses sum(z_k) as an integer... Let me re-read.
# Actually from the paper/memory:
# z_CC = sum(z_k) + 7 - alpha_s = 57 - alpha_s
# so sum(z_k) + 7 = 57, meaning sum(z_k) = 50.
# But sum(z_k) = 12 + 8*phi = 12 + 12.944 = 24.944, not 50.
# There must be TWO different definitions of z_k!

# Let me check: 2*a1^2 = 50. And from memory:
# "sum(z_k) = 12 + 8*phi = 2/alpha_s + rank(E8)"
# 2/alpha_s = 4*phi^3 = 4*(phi^2+phi) = 4*phi^2 + 4*phi
# rank(E8) = 8
# So 4*phi^2 + 4*phi + 8 ≈ 10.472 + 6.472 + 8 = 24.944

# The other sum: 2*a1^2 = 50 must be sum of DIFFERENT z_k.
# Perhaps z_k = a_k^2 + a_k*b_k - b_k^2 (norms)?

print(f"  sum(a_k) = 12 (degree)")
print(f"  sum(b_k) = 8 (rank E8)")
print(f"  sum(z_k) = sum(a_k) + sum(b_k)*phi = 12 + 8*phi = {12 + 8*phi:.4f}")
print(f"  2/alpha_s + rank(E8) = {2/alpha_s + 8:.4f}")
print(f"\n  For CC: sum(z_k) + 7 = 57 => sum(z_k) = 50 = 2*a1^2")
print(f"  But 12 + 8*phi = {12+8*phi:.4f} != 50")
print(f"  => These are DIFFERENT z_k (spectral weight vs something else)")
print(f"  NEED TO CHECK: which z_k sums to 50?")

# Possibly: the z_k in the CC formula are the NORMS N(z_k)?
# Or the TRACES Tr(z_k) = 2*a_k + b_k?

# Let me try: if z_k = a_k + b_k*phi with sum(a_k)=12, sum(b_k)=8
# Then N(z_k) = a_k^2 + a_k*b_k - b_k^2
# sum(N(z_k)) = sum(a_k^2) + sum(a_k*b_k) - sum(b_k^2)
# This would need the individual (a_k, b_k) values.

# Actually, from MEMORY: "z_P = spectral weight of R3 = 4*phi^2"
# The Planck weight is 4*phi^2 = z_Planck, and this is a single spectral weight.
# And: "z_CC = sum(z_k) + 7 - alpha_s = 57 - alpha_s"
# So if sum(z_k) = 50, the z_k must be something specific.

# From exp343: "z_k are spectral weights associated with each 2I irrep"
# These might be N(z_k) norms or Tr(z_k) traces.

# Let me try Tr(z_k) = 2*a_k + b_k:
# sum(Tr(z_k)) = 2*sum(a_k) + sum(b_k) = 24 + 8 = 32. Not 50.

# Try: a_k^2 + b_k^2:
# Without individual values, I can't compute this.

# The most natural meaning: z_k are the Laplacian eigenvalues themselves.
# L_k = 12 - lambda_k. Let me compute sum(L_k):
adj_evals_list = [12.0, 6*phi, 4*phi, 3.0, 0.0, -2.0, 4*phi_c, -3.0, 6*phi_c]
L_vals = [12 - ev for ev in adj_evals_list]
sum_L = sum(L_vals)
print(f"\n  Laplacian eigenvalues sum(12-lambda_k) = {sum_L:.4f}")
print(f"  = 9*12 - sum(lambda_k) = 108 - {sum(adj_evals_list):.4f} = {108 - sum(adj_evals_list):.4f}")

# Try: sum(L_k^2):
sum_L2 = sum(l**2 for l in L_vals)
print(f"  sum(L_k^2) = {sum_L2:.4f}")

# Try with multiplicities:
sum_L_weighted = sum((12-ev)*m for ev, m in
    [(12,1),(6*phi,4),(4*phi,9),(3,16),(0,25),(-2,36),(4*phi_c,9),(-3,16),(6*phi_c,4)])
print(f"  sum(L_k * mult_k) = {sum_L_weighted:.4f} = Tr(Laplacian)")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print("""
RESULTS:

1. Tr(Box^2)/N = N/2 = 60. So 57 = Tr(Box^2)/N - N_gen.
   This is EXACT but TAUTOLOGICAL (restates known formula).

2. The Frobenius number F(a1,b1) = 19 is the largest non-representable
   integer by a1*a + b1*b (a,b >= 0). And 57 = N_gen * 19.
   The 10 non-representable integers sum to 57.
   WAIT - let me check this...

3. Box spectral zeta and heat kernel do not give 57 directly.

4. The spectral weight sum needs clarification: which z_k sum to 50?

5. Z[phi] has no element with norm 57 for small (a,b).

OPEN QUESTION: The formula 57 = N_gen * F(a1,b1) = 3 * 19
connects generation count (3) with the Frobenius number (19).
The Frobenius number is the LARGEST mass gap in the n=5a+6b lattice.
The CC exponent counts "how many mass gaps" (weighted by generations)
exist in the framework's mass spectrum.

This is a STRUCTURAL observation but not yet a DERIVATION.
""")

# Check: do the 10 non-representable integers sum to 57?
print(f"  CHECK: Sum of non-representable integers = {sum_non_rep}")
print(f"  57? {'YES' if sum_non_rep == 57 else 'NO'}")
