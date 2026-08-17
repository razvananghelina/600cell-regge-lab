"""
exp557: Verify c² = a₁ = 5 from 600-cell Hopf fibration spectral gaps.

BLIND PROCEDURE: compute spectral gaps first, compare with a₁ and b₁ after.

Steps:
  1. Build 600-cell (120 vertices, degree 12, 720 edges)
  2. Construct Hopf fibrations via 2I quaternion multiplication (order-10 elements)
  3. Build L_fiber and L_cross, verify L_full = L_fiber + L_cross
  4. Compute spectral gaps of L_full, L_fiber, L_cross
  5. Check [L_fiber, L_cross] commutativity
  6. Gap eigenvector analysis
  7. Verify □ = L_cross - a₁·L_fiber
  8. Check □ kernel dimension (expect 9 from Theorem 9)
  9. Repeat for all 27 Hopf fibrations
"""

import numpy as np
from numpy.linalg import eigh, norm
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
TOL = 1e-8

# ─── Quaternion arithmetic ───────────────────────────────────────────

def qmul(p, q):
    """Hamilton product of two quaternions (w, x, y, z)."""
    w = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
    x = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
    y = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
    z = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]
    return np.array([w, x, y, z])


def qconj(q):
    """Quaternion conjugate."""
    return np.array([q[0], -q[1], -q[2], -q[3]])


def find_vertex_index(v, verts, tol=1e-6):
    """Find index of vertex v in the vertex list."""
    dots = verts @ v
    idx = np.argmax(dots)
    if dots[idx] > 1 - tol:
        return idx
    # Also check antipodal (unit quaternions: q and -q are distinct in 2I)
    return -1


# ─── Step 1: Build 600-cell ─────────────────────────────────────────

print("=" * 70)
print("STEP 1: Build 600-cell")
print("=" * 70)

verts, adj, lap = build_600cell()
n_verts = len(verts)
n_edges = int(np.sum(adj) / 2)
degree = int(np.sum(adj[0]))

print(f"  Vertices: {n_verts} (expect 120)")
print(f"  Edges:    {n_edges} (expect 720)")
print(f"  Degree:   {degree} (expect 12)")

# Verify spectrum
evals_full = np.sort(eigh(lap)[0])
distinct_full = []
for e in evals_full:
    if not distinct_full or abs(e - distinct_full[-1][0]) > 1e-6:
        distinct_full.append([e, 1])
    else:
        distinct_full[-1][1] += 1

print(f"  Distinct eigenvalues: {len(distinct_full)} (expect 9)")
print(f"  Spectrum:")
for val, mult in distinct_full:
    print(f"    λ = {val:12.6f}  mult = {mult}")

gap_full = distinct_full[1][0]
gap_full_exact = 12 - 6*PHI
print(f"\n  Spectral gap λ₁(L_full) = {gap_full:.10f}")
print(f"  Expected 12 - 6φ        = {gap_full_exact:.10f}")
print(f"  Match: {abs(gap_full - gap_full_exact) < TOL}")


# ─── Step 2: Construct Hopf fibrations ──────────────────────────────

print("\n" + "=" * 70)
print("STEP 2: Find Hopf fibrations via order-10 elements of 2I")
print("=" * 70)

# Find all elements of order 10 in 2I.
# An element g has order 10 if g^10 = ±identity and g^5 ≠ ±identity.
# For unit quaternion: order 10 iff trace(g) = 2cos(π/5) = φ or 2cos(3π/5) = -1/φ
# (since g^10 = id means angle = 2πk/10, so cos(πk/5) for k=1,...,4,
#  but order exactly 10 means k=1 or k=3: cos(π/5)=φ/2, cos(3π/5)=-1/(2φ))

def qpow(q, n):
    """Compute q^n by repeated multiplication."""
    result = np.array([1.0, 0.0, 0.0, 0.0])
    for _ in range(n):
        result = qmul(result, q)
    return result

# Identity element
identity = np.array([1.0, 0.0, 0.0, 0.0])

def is_close_to_vertex(q, v, tol=1e-6):
    """Check if quaternion q is close to vertex v (or -v since they're both in 2I)."""
    return abs(np.dot(q, v)) > 1 - tol

# Find elements of order 10
order10_elements = []
for i, v in enumerate(verts):
    trace = 2 * v[0]  # trace of quaternion = 2*w
    # Order 10: trace = 2cos(π/5) = φ or trace = 2cos(3π/5) = -(φ-1)
    if abs(abs(trace) - PHI) < 1e-6 or abs(abs(trace) - (PHI - 1)) < 1e-6:
        # Verify: g^10 = +1 and g^5 != +1 (g^5 = -1 is OK for order 10)
        g5 = qpow(v, 5)
        g10 = qpow(v, 10)
        is_order10 = (np.dot(g10, identity) > 1 - 1e-6 and
                      np.dot(g5, identity) < 1 - 1e-6)
        if is_order10:
            order10_elements.append(i)

print(f"  Elements of order 10 in 2I: {len(order10_elements)}")

# Build Hopf fibrations: a fiber through vertex q₀ using generator g
# is {q₀, q₀·g, q₀·g², ..., q₀·g⁹}
# Two generators g and g' give the same fibration iff they generate the same
# cyclic subgroup (g' = g^k with gcd(k,10)=1).

def build_fibration(gen_idx):
    """Build a fibration using the given generator.
    Returns list of 12 fibers, each a sorted tuple of 10 vertex indices."""
    g = verts[gen_idx]
    used = set()
    fibers = []

    for start_idx in range(n_verts):
        if start_idx in used:
            continue
        # Generate fiber: {v_start · g^k : k=0,...,9}
        fiber = []
        q = verts[start_idx].copy()
        for k in range(10):
            idx = find_vertex_index(q, verts)
            if idx < 0:
                break
            fiber.append(idx)
            q = qmul(q, g)

        if len(fiber) == 10 and len(set(fiber)) == 10:
            fiber_set = frozenset(fiber)
            if not fiber_set & used:
                fibers.append(tuple(sorted(fiber)))
                used |= fiber_set

    return fibers if len(fibers) == 12 and len(used) == 120 else None


# Also build fibrations using RIGHT multiplication: {g^k · q₀}
def build_fibration_right(gen_idx):
    """Build a fibration using RIGHT action: fiber = {g^k · v : k=0,...,9}."""
    g = verts[gen_idx]
    used = set()
    fibers = []

    for start_idx in range(n_verts):
        if start_idx in used:
            continue
        fiber = []
        q = verts[start_idx].copy()
        for k in range(10):
            idx = find_vertex_index(q, verts)
            if idx < 0:
                break
            fiber.append(idx)
            q = qmul(g, q)  # LEFT multiply by g (right action on cosets)

        if len(fiber) == 10 and len(set(fiber)) == 10:
            fiber_set = frozenset(fiber)
            if not fiber_set & used:
                fibers.append(tuple(sorted(fiber)))
                used |= fiber_set

    return fibers if len(fibers) == 12 and len(used) == 120 else None


# Find all distinct fibrations using left, right, and conjugated generators
all_fibrations = []
fibration_signatures = set()

# Try all order-10 generators with left and right actions
for idx in order10_elements:
    for builder in [build_fibration, build_fibration_right]:
        fibs = builder(idx)
        if fibs is not None:
            sig = frozenset(frozenset(f) for f in fibs)
            if sig not in fibration_signatures:
                fibration_signatures.add(sig)
                all_fibrations.append(fibs)

# Also try conjugated generators: h·g·h^{-1} for all h in 2I
for idx_g in order10_elements:
    g = verts[idx_g]
    for idx_h in range(n_verts):
        h = verts[idx_h]
        h_inv = qconj(h)
        g_conj = qmul(qmul(h, g), h_inv)
        # Find this conjugate in the vertex list
        conj_idx = find_vertex_index(g_conj, verts)
        if conj_idx >= 0:
            for builder in [build_fibration, build_fibration_right]:
                fibs = builder(conj_idx)
                if fibs is not None:
                    sig = frozenset(frozenset(f) for f in fibs)
                    if sig not in fibration_signatures:
                        fibration_signatures.add(sig)
                        all_fibrations.append(fibs)

print(f"  Distinct Hopf fibrations found: {len(all_fibrations)} (expect 27)")

# Verify first fibration
if all_fibrations:
    fibs0 = all_fibrations[0]
    all_verts_in_fibers = set()
    for f in fibs0:
        all_verts_in_fibers |= set(f)
    print(f"  First fibration: {len(fibs0)} fibers × {len(fibs0[0])} vertices = {len(all_verts_in_fibers)} total")


# ─── Step 3: Build L_fiber and L_cross ──────────────────────────────

print("\n" + "=" * 70)
print("STEP 3: Build L_fiber, L_cross for first Hopf fibration")
print("=" * 70)


def build_fiber_matrices(fibration, adj_mat, lap_mat):
    """Build fiber and cross-fiber adjacency and Laplacian matrices."""
    n = len(adj_mat)
    adj_fiber = np.zeros((n, n))

    # Fiber edges: edges between vertices in the same fiber
    for fiber in fibration:
        fiber_set = set(fiber)
        for i in fiber_set:
            for j in fiber_set:
                if i != j and adj_mat[i, j] > 0.5:
                    adj_fiber[i, j] = 1.0

    # Cross-fiber edges: all other edges
    adj_cross = adj_mat - adj_fiber

    # Verify: no negative entries
    assert np.all(adj_cross >= -TOL), "Negative entries in adj_cross!"
    adj_cross = np.clip(adj_cross, 0, 1)

    # Degrees
    deg_fiber = np.sum(adj_fiber, axis=1)
    deg_cross = np.sum(adj_cross, axis=1)

    # Laplacians
    lap_fiber = np.diag(deg_fiber) - adj_fiber
    lap_cross = np.diag(deg_cross) - adj_cross

    return adj_fiber, adj_cross, lap_fiber, lap_cross


fibs0 = all_fibrations[0]
adj_fiber, adj_cross, lap_fiber, lap_cross = build_fiber_matrices(fibs0, adj, lap)

# Verify degrees
deg_fiber = np.sum(adj_fiber, axis=1)
deg_cross = np.sum(adj_cross, axis=1)
print(f"  Fiber degree:  min={deg_fiber.min():.0f}, max={deg_fiber.max():.0f} (expect 2)")
print(f"  Cross degree:  min={deg_cross.min():.0f}, max={deg_cross.max():.0f} (expect 10)")

# Verify L_full = L_fiber + L_cross
lap_sum = lap_fiber + lap_cross
diff = norm(lap - lap_sum)
print(f"  ||L_full - (L_fiber + L_cross)|| = {diff:.2e} (expect 0)")
print(f"  Decomposition exact: {diff < TOL}")


# ─── Step 4: Compute spectral gaps (BLIND) ──────────────────────────

print("\n" + "=" * 70)
print("STEP 4: Spectral gaps (BLIND computation)")
print("=" * 70)

evals_fiber, evecs_fiber = eigh(lap_fiber)
evals_cross, evecs_cross = eigh(lap_cross)
evals_full_sorted, evecs_full = eigh(lap)

# Spectral gaps (smallest nonzero eigenvalue)
def spectral_gap(evals):
    """Return smallest eigenvalue > TOL."""
    pos = evals[evals > TOL]
    return pos[0] if len(pos) > 0 else None

gap_fiber = spectral_gap(evals_fiber)
gap_cross = spectral_gap(evals_cross)
gap_full_check = spectral_gap(evals_full_sorted)

print(f"  λ₁(L_full)  = {gap_full_check:.10f}")
print(f"  λ₁(L_fiber) = {gap_fiber:.10f}")
print(f"  λ₁(L_cross) = {gap_cross:.10f}")

print(f"\n  --- BLIND RATIOS (computed before comparing with a₁, b₁) ---")
ratio_full_fiber = gap_full_check / gap_fiber
ratio_cross_fiber = gap_cross / gap_fiber
print(f"  λ₁(L_full)  / λ₁(L_fiber) = {ratio_full_fiber:.10f}")
print(f"  λ₁(L_cross) / λ₁(L_fiber) = {ratio_cross_fiber:.10f}")

print(f"\n  --- NOW COMPARE ---")
print(f"  b₁ = {b1},  ratio_full/fiber  = {ratio_full_fiber:.10f},  match: {abs(ratio_full_fiber - b1) < TOL}")
print(f"  a₁ = {a1},  ratio_cross/fiber = {ratio_cross_fiber:.10f},  match: {abs(ratio_cross_fiber - a1) < TOL}")

# Expected exact values
print(f"\n  Expected values:")
print(f"    λ₁(L_fiber) = 2 - φ       = {2 - PHI:.10f}")
print(f"    λ₁(L_cross) = 10 - 5φ     = {10 - 5*PHI:.10f}")
print(f"    λ₁(L_full)  = 12 - 6φ     = {12 - 6*PHI:.10f}")
print(f"  Fiber gap match (2-φ):  {abs(gap_fiber - (2-PHI)) < TOL}")
print(f"  Cross gap match (10-5φ): {abs(gap_cross - (10-5*PHI)) < TOL}")

# Full fiber spectrum
print(f"\n  Full L_fiber spectrum (distinct eigenvalues):")
distinct_fiber = []
for e in evals_fiber:
    if not distinct_fiber or abs(e - distinct_fiber[-1][0]) > 1e-6:
        distinct_fiber.append([e, 1])
    else:
        distinct_fiber[-1][1] += 1
for val, mult in distinct_fiber:
    print(f"    λ = {val:12.6f}  mult = {mult}")

print(f"\n  Full L_cross spectrum (distinct eigenvalues):")
distinct_cross = []
for e in evals_cross:
    if not distinct_cross or abs(e - distinct_cross[-1][0]) > 1e-6:
        distinct_cross.append([e, 1])
    else:
        distinct_cross[-1][1] += 1
for val, mult in distinct_cross:
    print(f"    λ = {val:12.6f}  mult = {mult}")


# ─── Step 5: Commutativity ──────────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 5: Check [L_fiber, L_cross] commutativity")
print("=" * 70)

commutator = lap_fiber @ lap_cross - lap_cross @ lap_fiber
comm_norm = norm(commutator)
print(f"  ||[L_fiber, L_cross]||_F = {comm_norm:.10f}")
print(f"  Commute: {comm_norm < TOL}")

if comm_norm > TOL:
    comm_max = np.max(np.abs(commutator))
    print(f"  Max entry of [L_fiber, L_cross] = {comm_max:.10f}")


# ─── Step 6: Gap eigenvector analysis ───────────────────────────────

print("\n" + "=" * 70)
print("STEP 6: Gap eigenvector analysis")
print("=" * 70)

# Find gap eigenvectors of L_full (mult = 4)
gap_indices_full = np.where(np.abs(evals_full_sorted - gap_full_check) < 1e-6)[0]
print(f"  L_full gap multiplicity: {len(gap_indices_full)}")

for i, idx in enumerate(gap_indices_full):
    v = evecs_full[:, idx]
    expectation_fiber = v @ lap_fiber @ v
    expectation_cross = v @ lap_cross @ v
    print(f"  v_{i}: <v|L_fiber|v> = {expectation_fiber:.10f},  "
          f"<v|L_cross|v> = {expectation_cross:.10f},  "
          f"sum = {expectation_fiber + expectation_cross:.10f}")

print(f"\n  Expected: <v|L_fiber|v> = 2 - φ  = {2 - PHI:.10f}")
print(f"  Expected: <v|L_cross|v> = 10 - 5φ = {10 - 5*PHI:.10f}")


# ─── Step 7: Verify □ = L_cross - a₁·L_fiber ───────────────────────

print("\n" + "=" * 70)
print("STEP 7: Verify □ = L_cross - a₁·L_fiber")
print("=" * 70)

# From paper: □ = b₁·A_fiber - A
# Where A = adj, A_fiber = adj restricted to fiber edges
# □ = b₁·A_fiber - A
box_from_def = b1 * adj_fiber - adj
# Also: □ = L_cross - a₁·L_fiber (algebraic identity)
box_from_lap = lap_cross - a1 * lap_fiber

diff_box = norm(box_from_def - box_from_lap)
print(f"  ||□_def - □_lap|| = {diff_box:.2e}")
print(f"  Identity □ = L_cross - a₁·L_fiber holds: {diff_box < TOL}")

# Alternative derivation to verify:
# □ = b₁·A_fiber - A
#   = b₁·A_fiber - (A_fiber + A_cross)
#   = (b₁ - 1)·A_fiber - A_cross
#   = a₁·A_fiber - A_cross
# And L_cross - a₁·L_fiber
#   = (d_cross·I - A_cross) - a₁·(d_fiber·I - A_fiber)
#   = (10·I - A_cross) - a₁·(2·I - A_fiber)
#   = 10·I - A_cross - 10·I + 5·A_fiber  [since a₁=5, d_fiber=2]
#   = 5·A_fiber - A_cross
#   = a₁·A_fiber - A_cross ✓
# So □ = a₁·A_fiber - A_cross = L_cross - a₁·L_fiber ✓

# Verify the algebraic identity step by step
box_alt = a1 * adj_fiber - adj_cross
diff_alt = norm(box_from_def - box_alt)
print(f"  Also: □ = a₁·A_fiber - A_cross, diff = {diff_alt:.2e}")


# ─── Step 8: □ kernel dimension ─────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 8: □ kernel dimension (Theorem 9)")
print("=" * 70)

evals_box = np.sort(np.linalg.eigvalsh(box_from_lap))
kernel_dim = np.sum(np.abs(evals_box) < 1e-6)
print(f"  dim(ker □) = {kernel_dim} (expect 9 from Theorem 9)")

# Show full □ spectrum
distinct_box = []
for e in evals_box:
    if not distinct_box or abs(e - distinct_box[-1][0]) > 1e-6:
        distinct_box.append([e, 1])
    else:
        distinct_box[-1][1] += 1

print(f"  □ spectrum ({len(distinct_box)} distinct eigenvalues):")
for val, mult in distinct_box:
    marker = " <-- KERNEL" if abs(val) < 1e-6 else ""
    print(f"    μ = {val:12.6f}  mult = {mult}{marker}")


# ─── Step 9: All 27 fibrations ──────────────────────────────────────

print("\n" + "=" * 70)
print(f"STEP 9: Repeat for all {len(all_fibrations)} Hopf fibrations")
print("=" * 70)

results = []
for i, fibs in enumerate(all_fibrations):
    af, ac, lf, lc = build_fiber_matrices(fibs, adj, lap)

    ef = np.sort(eigh(lf)[0])
    ec = np.sort(eigh(lc)[0])

    gf = spectral_gap(ef)
    gc = spectral_gap(ec)
    ratio = gc / gf if gf > TOL else float('nan')

    box = lc - a1 * lf
    ebox = np.sort(np.linalg.eigvalsh(box))
    kdim = np.sum(np.abs(ebox) < 1e-6)

    comm = norm(lf @ lc - lc @ lf)

    results.append({
        'idx': i,
        'gap_fiber': gf,
        'gap_cross': gc,
        'ratio': ratio,
        'kernel_dim': kdim,
        'commutator_norm': comm,
    })

    match_str = "YES" if abs(ratio - a1) < TOL else "NO"
    print(f"  Fibration {i+1:2d}: λ₁_fiber={gf:.6f}, λ₁_cross={gc:.6f}, "
          f"ratio={ratio:.6f}, c²=a₁? {match_str}, "
          f"ker□={kdim}, ||[Lf,Lc]||={comm:.2e}")


# ─── Summary ────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

n_match = sum(1 for r in results if abs(r['ratio'] - a1) < TOL)
n_kernel9 = sum(1 for r in results if r['kernel_dim'] == 9)
n_commute = sum(1 for r in results if r['commutator_norm'] < TOL)

print(f"  Total fibrations found:        {len(all_fibrations)}")
print(f"  c² = a₁ = 5 (exact):          {n_match}/{len(results)}")
print(f"  ker□ = 9:                      {n_kernel9}/{len(results)}")
print(f"  [L_fiber, L_cross] = 0:        {n_commute}/{len(results)}")
print()

if n_match == len(results) and len(results) > 0:
    print("  *** RESULT: c² = a₁ = 5 CONFIRMED for ALL fibrations ***")
    print(f"  The speed of light c = √a₁ = √5 = {np.sqrt(5):.10f} in lattice units")
    print(f"  Note: √5 = 2φ - 1 = {2*PHI - 1:.10f}")
elif n_match > 0:
    print(f"  PARTIAL: c² = a₁ for {n_match} of {len(results)} fibrations")
else:
    print("  NEGATIVE: c² ≠ a₁ for any fibration")

if n_match == 0 and len(results) > 0:
    # Report what the ratio actually is
    ratios = [r['ratio'] for r in results]
    print(f"\n  Actual ratios: min={min(ratios):.6f}, max={max(ratios):.6f}")
    print(f"  Mean ratio = {np.mean(ratios):.6f}")
