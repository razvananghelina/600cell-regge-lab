"""
exp103: Two focused questions about Box_edge spectrum.

Q1: Do the first massive modes encode W/Z?
    - Multiplicity 3? Mass ratio cos(theta_W)?

Q2: Trace moments Tr(Box_edge^n) — framework identities?
    - Are they all integers? (Box_edge is integer matrix => YES trivially)
    - Does Tr(Box_1^3) = something clean?
    - Bootstrap analog: Tr(Box_0^3) = N^2 on vertices. What about edges?
"""

import numpy as np
from numpy.linalg import eigh, norm
from collections import defaultdict, Counter
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
N = 120
SQRT5 = np.sqrt(5)

# =====================================================================
# BUILD (compact, reuse from exp100)
# =====================================================================
print("Building...")
verts, adj, lap = build_600cell()
Nv = 120

adj_list = defaultdict(set)
edges = []; edge_to_idx = {}
for i in range(Nv):
    for j in range(i+1, Nv):
        if adj[i, j] > 0.5:
            adj_list[i].add(j); adj_list[j].add(i)
            edge_to_idx[(i, j)] = len(edges); edges.append((i, j))
Ne = 720

def qmul(p, q):
    return np.array([p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1-tol else -1

# Hopf fibration
fibers = None
for i in range(Nv):
    if abs(verts[i, 0] - PHI/2) >= 1e-6: continue
    g = verts[i]; p = g.copy(); ok = True
    for k in range(2, 11):
        p = qmul(p, g)
        if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok = False; break
        if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok = False
    if not ok: continue
    used = set(); fibs = []; subg = []
    pp = np.array([1.0,0,0,0])
    for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
    for s in range(Nv):
        if s in used: continue
        fib = []
        for si in subg:
            idx = find_idx(qmul(verts[s], verts[si]), verts)
            if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
        if len(fib) == 10: fibs.append(fib)
    if len(fibs) == 12: fibers = fibs; break

vtf = {}
for fi, f in enumerate(fibers):
    for v in f: vtf[v] = fi

is_fiber = np.zeros(Ne, dtype=bool)
for e, (i, j) in enumerate(edges):
    if vtf[i] == vtf[j]: is_fiber[e] = True

vte = defaultdict(list)
for e, (i, j) in enumerate(edges): vte[i].append(e); vte[j].append(e)

A_line = np.zeros((Ne, Ne))
for v in range(Nv):
    inc = vte[v]
    for a in range(len(inc)):
        for b in range(a+1, len(inc)):
            A_line[inc[a], inc[b]] = 1.0; A_line[inc[b], inc[a]] = 1.0

A_fib_line = np.zeros((Ne, Ne))
for fi, fiber in enumerate(fibers):
    fe = []
    for k in range(10):
        u, v = fiber[k], fiber[(k+1)%10]
        fe.append(edge_to_idx[(min(u,v), max(u,v))])
    for k in range(10):
        A_fib_line[fe[k], fe[(k+1)%10]] = 1.0; A_fib_line[fe[(k+1)%10], fe[k]] = 1.0

Lf = np.diag(np.sum(A_fib_line, axis=1)) - A_fib_line
Lc = np.diag(np.sum(A_line - A_fib_line, axis=1)) - (A_line - A_fib_line)
Box1 = Lc - a1 * Lf

# Vertex Box
Afv = adj * np.array([[float(vtf[i]==vtf[j]) for j in range(Nv)] for i in range(Nv)])
Lfv = np.diag(Afv.sum(1)) - Afv; Lcv = np.diag((adj-Afv).sum(1)) - (adj-Afv)
Box0 = Lcv - a1 * Lfv

evals1 = np.sort(np.linalg.eigvalsh(Box1))
evals0 = np.sort(np.linalg.eigvalsh(Box0))

print("Done.\n")


# =====================================================================
# Q1: MASSIVE MODES — W/Z SEARCH
# =====================================================================
print("=" * 72)
print("Q1: MASSIVE MODES NEAR KERNEL")
print("=" * 72)

# Distinct eigenvalues with multiplicities
ev1_dist = Counter(np.round(evals1, 4))
ev0_dist = Counter(np.round(evals0, 4))

print("\n  First 10 eigenvalues of Box_edge (above kernel):")
print(f"  {'lambda':>10s} {'mult':>5s} {'irrep':>20s} {'sqrt(lambda)':>12s}")
print(f"  {'-'*55}")

# From exp100d: irrep assignments
irrep_map = {
    0.0: "ker: rho_0+2*rho_5",
    1.502: "2*rho_4(5)", 1.6454: "2*rho_6(3)", 1.6698: "2*rho_7(4)", 1.8737: "2*rho_1(2)",
    3.7175: "2*rho_2(3)", 3.9733: "2*rho_3(4)",
    5.5279: "1*rho_2(3)", 5.682: "2*rho_5(6)",
}

for val, mult in sorted(ev1_dist.items()):
    if val > 6: break
    irr = ""
    for k, v in irrep_map.items():
        if abs(val - k) < 0.01: irr = v; break
    print(f"  {val:10.4f} {mult:5d} {irr:>20s} {np.sqrt(abs(val)):12.6f}")

# Check: any multiplicity 3 or 4 near the kernel?
print(f"\n  Multiplicities in first band (lambda < 2):")
band1 = [(v, m) for v, m in sorted(ev1_dist.items()) if 0.01 < v < 2]
for v, m in band1:
    print(f"    {v:10.4f} x {m}")
total_band1 = sum(m for _, m in band1)
print(f"  Total modes in band 1: {total_band1}")

# W/Z mass ratio test
cos2_tW = 20.0/26.0  # sin^2(tW) = 6/26 => cos^2(tW) = 20/26
print(f"\n  Expected m_W^2/m_Z^2 = cos^2(theta_W) = 20/26 = {cos2_tW:.6f}")
print(f"  Expected m_Z^2/m_W^2 = 26/20 = {26/20:.6f}")

# Check all ratios in first band
print(f"\n  Eigenvalue ratios in first band:")
band_vals = [v for v, _ in band1]
for i in range(len(band_vals)):
    for j in range(i+1, len(band_vals)):
        r = band_vals[j] / band_vals[i]
        match = ""
        if abs(r - 26/20) < 0.05: match = " *** = 26/20 ***"
        if abs(r - 20/26) < 0.05: match = " *** = 20/26 ***"
        if abs(r - PHI) < 0.05: match = " ~= phi"
        if abs(r - 1/PHI) < 0.05: match = " ~= 1/phi"
        mi, mj = ev1_dist[band_vals[i]], ev1_dist[band_vals[j]]
        if match or abs(r - 1) > 0.01:
            print(f"    {band_vals[j]:.4f}/{band_vals[i]:.4f} = {r:.6f} (mults {mi},{mj}){match}")

# EW sector: the SU(2)xU(1) has dim 4 = 3+1
# Before symmetry breaking: 4 massless
# After: 3 massive (W+,W-,Z) + 1 massless (photon)
# Do any irreps give 3+1?
print(f"\n  Looking for 3+1 = 4 pattern in first band:")
print(f"  rho_1(2): dim 4 via 2*rho_1 => 4 modes at lambda = 1.8737")
print(f"  rho_6(3): dim 6 via 2*rho_6 => 6 modes at lambda = 1.6454")
print(f"  rho_7(4): dim 8 via 2*rho_7 => 8 modes at lambda = 1.6698")
print(f"  rho_4(5): dim 10 via 2*rho_4 => 10 modes at lambda = 1.5020")
print(f"  NONE has multiplicity exactly 3.")


# Alternative: look at VERTEX Box massive modes for comparison
print(f"\n  Vertex Box_0 first massive modes:")
for val, mult in sorted(ev0_dist.items()):
    if abs(val) > 0.01:
        print(f"    {val:10.4f} x {mult}")
        if val > 6: break


# =====================================================================
# Q2: TRACE MOMENTS Tr(Box_edge^n)
# =====================================================================
print("\n" + "=" * 72)
print("Q2: TRACE MOMENTS Tr(Box_edge^n)")
print("=" * 72)

# Box1 is integer => all Tr(Box1^n) are integers
Box1_int = np.round(Box1).astype(np.int64)

# Compute Tr(Box^n) using eigenvalues (faster for high n)
# But also verify with matrix powers for small n
print("\n  Computing trace moments (using matrix powers for verification)...")

M = np.eye(Ne)
moments_edge = []
moments_vert = []

M0 = np.eye(Nv)
for n in range(1, 11):
    M = M @ Box1_int
    tr_edge = int(np.trace(M))
    moments_edge.append(tr_edge)

    M0 = M0 @ Box0
    tr_vert = round(np.trace(M0))
    moments_vert.append(tr_vert)

print(f"\n  {'n':>3s} {'Tr(Box_1^n)':>18s} {'Tr(Box_0^n)':>15s} {'ratio':>12s} {'factored':>30s}")
print(f"  {'-'*85}")

for n in range(10):
    tr_e = moments_edge[n]
    tr_v = moments_vert[n]
    ratio = tr_e / tr_v if tr_v != 0 else float('inf')

    # Try to factor the edge moment
    factored = ""
    # Check: Tr(Box1^n) = c * N^k for some c, k?
    for k in range(1, 6):
        if N**k != 0 and tr_e % (N**k) == 0:
            c = tr_e // (N**k)
            factored = f"{c} * {N}^{k}"
            break
    # Check: = c * Ne^k?
    if not factored:
        for k in range(1, 4):
            if Ne**k != 0 and tr_e % (Ne**k) == 0:
                c = tr_e // (Ne**k)
                factored = f"{c} * {Ne}^{k}"
                break

    print(f"  {n+1:3d} {tr_e:18d} {tr_v:15d} {ratio:12.4f} {factored:>30s}")

# Key checks
print(f"\n  KEY CHECKS:")
print(f"  Tr(Box_0^1) = {moments_vert[0]} (expect 0: traceless)")
print(f"  Tr(Box_1^1) = {moments_edge[0]} (expect 14400 = N^2)")
print(f"  Tr(Box_0^2) = {moments_vert[1]} (= 7200 = 60*N)")
print(f"  Tr(Box_1^2) = {moments_edge[1]} (= 324000 = 45*7200 = 9*a1*Tr(Box_0^2))")
print(f"  Tr(Box_0^3) = {moments_vert[2]} (= 14400 = N^2 -- BOOTSTRAP!)")
print(f"  Tr(Box_1^3) = {moments_edge[2]}")

# Analyze Tr(Box_1^3)
tr3 = moments_edge[2]
print(f"\n  ANALYZING Tr(Box_1^3) = {tr3}:")
print(f"    / N^2 = {tr3 / N**2:.4f}")
print(f"    / N^3 = {tr3 / N**3:.4f}")
print(f"    / Ne^2 = {tr3 / Ne**2:.4f}")
print(f"    / Ne*N = {tr3 / (Ne*N):.4f}")
print(f"    / Tr(Box_0^3) = {tr3 / moments_vert[2]:.4f}")

# Factor tr3
print(f"    Factorization: ", end="")
n = abs(tr3)
factors = []
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    while n % p == 0:
        factors.append(p)
        n //= p
if n > 1: factors.append(n)
print(f"  {tr3} = {' * '.join(str(f) for f in factors)}")

# Framework constant identification
print(f"\n  Tr(Box_1^3) / Tr(Box_0^3) = {moments_edge[2] / moments_vert[2]:.6f}")
r31 = moments_edge[2] / moments_vert[2]
# Check if ratio is a framework constant
for name, val in [("N/2", 60), ("a1*N/2", 300), ("Ne/2", 360),
                   ("a1^2*a1", 125), ("9*a1^2", 225),
                   ("(Ne/N)^3", (720/120)**3), ("521.5", 521.5),
                   ("a1*Ne/N", 30), ("b1*Ne/N", 36)]:
    if abs(r31 - val) < 0.5:
        print(f"    *** matches {name} = {val} ***")


# =====================================================================
# Q2b: GALOIS INVARIANCE OF MOMENTS
# =====================================================================
print("\n" + "=" * 72)
print("Q2b: GALOIS INVARIANCE")
print("=" * 72)

# All Tr(Box1^n) are integers (Box1 is integer matrix).
# But are the INDIVIDUAL eigenvalue moments Galois-invariant?
# i.e., do Galois-conjugate eigenvalues appear in pairs?

# The eigenvalues come in pairs related by phi -> -1/phi?
# For vertex Box: eigenvalues are in Z[phi], Galois pairs give same moments.
# For edge Box: eigenvalues are NOT all in Z[phi] (only 6/30 were).

# Check: for each eigenvalue lambda, is there a "Galois conjugate" lambda'
# such that lambda + lambda' is rational and lambda * lambda' is rational?
print("\n  Eigenvalue Galois pairing analysis:")
distinct = sorted(set(np.round(evals1, 6)))
paired = set()
galois_pairs = []

for i, lam in enumerate(distinct):
    if i in paired: continue
    # Try to find a Galois partner: lambda' such that lambda + lambda' in Z
    # and lambda * lambda' in Z
    found_partner = False
    for j, lam2 in enumerate(distinct):
        if j in paired or j == i: continue
        s = lam + lam2
        p = lam * lam2
        if abs(s - round(s)) < 0.01 and abs(p - round(p)) < 0.01:
            mult_i = ev1_dist[round(lam, 4)]
            mult_j = ev1_dist[round(lam2, 4)]
            if mult_i == mult_j:  # Galois pairs should have same multiplicity
                galois_pairs.append((lam, lam2, int(round(s)), int(round(p)), mult_i))
                paired.add(i); paired.add(j)
                found_partner = True
                break
    if not found_partner:
        # Self-paired (rational eigenvalue)
        mult = ev1_dist[round(lam, 4)]
        if abs(lam - round(lam)) < 0.01:
            galois_pairs.append((lam, lam, int(round(2*lam)), int(round(lam**2)), mult))
            paired.add(i)

print(f"  {'lambda':>10s} {'lambda_G':>10s} {'sum':>6s} {'prod':>8s} {'mult':>5s}")
print(f"  {'-'*45}")
for lam, lam2, s, p, m in galois_pairs:
    tag = "(self)" if abs(lam - lam2) < 0.01 else ""
    print(f"  {lam:10.4f} {lam2:10.4f} {s:6d} {p:8d} {m:5d} {tag}")

n_paired = len([1 for l1, l2, _, _, _ in galois_pairs if abs(l1-l2) > 0.01])
n_self = len([1 for l1, l2, _, _, _ in galois_pairs if abs(l1-l2) < 0.01])
n_unpaired = len(distinct) - 2*n_paired - n_self
print(f"\n  Galois pairs: {n_paired}, self-paired (rational): {n_self}, unpaired: {n_unpaired}")


# =====================================================================
# Q2c: BOOTSTRAP IDENTITY SEARCH
# =====================================================================
print("\n" + "=" * 72)
print("Q2c: BOOTSTRAP IDENTITY SEARCH")
print("=" * 72)

# On vertices: Tr(Box_0^3) = N^2 = 14400, and this is equivalent to a1=5.
# On edges: what identity does Tr(Box_1^n) satisfy?

# Systematic search: Tr(Box_1^n) = f(a1, b1, N, Ne, ...)?
print("\n  Searching for moment identities...")

# Known: Tr(Box_1^1) = N^2 = 14400
# Known: Tr(Box_1^2) = 9*a1 * Tr(Box_0^2) = 45 * 7200 = 324000

# Tr(Box_1^n) for n=1..10
for n in range(10):
    tr = moments_edge[n]
    # Try: tr = c * N^a * Ne^b * a1^c ...
    # Simple: tr / (N^2) and tr / (Ne^2)
    r_N2 = tr / N**2
    r_Ne2 = tr / Ne**2
    r_N2Ne = tr / (N**2 * Ne)
    print(f"\n  n={n+1}: Tr = {tr}")
    print(f"    /N^2 = {r_N2:.2f}, /Ne^2 = {r_Ne2:.4f}, /N^3 = {tr/N**3:.4f}")

    # Check if it's a product of small framework constants
    for label, val in [
        ("N^2", N**2), ("N^3", N**3), ("N^2*Ne", N**2*Ne),
        ("Ne^2*N", Ne**2*N), ("Ne^3", Ne**3),
        ("N^2*a1", N**2*a1), ("N^2*b1", N**2*b1),
        ("N^2*degree", N**2*12), ("N^2*a1^2", N**2*a1**2),
        ("9*a1*Tr(Box0^n)", 9*a1*moments_vert[n] if n < 10 else 0),
    ]:
        if val != 0 and abs(tr / val - round(tr / val)) < 0.01 and abs(tr/val) < 10000:
            print(f"    = {int(round(tr/val))} * {label}")


# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"""
  Q1 (W/Z from massive modes):
  - First band has multiplicities 10, 6, 8, 4 (total 28)
  - NO multiplicity-3 mode found
  - NO eigenvalue ratio matches cos^2(theta_W) = 20/26
  - W/Z mass splitting NOT visible in Box_edge spectrum
  - STATUS: NEGATIVE for direct W/Z identification

  Q2 (Trace moments):
  - All Tr(Box_1^n) are integers (trivially, integer matrix)
  - Tr(Box_1^1) = {moments_edge[0]} = N^2
  - Tr(Box_1^2) = {moments_edge[1]} = 9*a1 * Tr(Box_0^2)
  - Tr(Box_1^3) = {moments_edge[2]}
  - Tr(Box_1^3)/Tr(Box_0^3) = {moments_edge[2]/moments_vert[2]:.4f}
""")

print("=" * 72)
print("EXP103 COMPLETE")
print("=" * 72)
