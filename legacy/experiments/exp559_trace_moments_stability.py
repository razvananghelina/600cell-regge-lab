"""
exp559: Verify trace moments of Box are fibration-independent.

From exp558:
  Tr(Box^1)/N = 0
  Tr(Box^2)/N = 60 = degree * a1
  Tr(Box^3)   = N^2 = 14400
  Tr(Box^5)   = (2*a1)! = 10! = 3628800

Are these exact theorems or fibration-dependent numerical accidents?
Check across all fibrations.
"""

import numpy as np
from numpy.linalg import eigvalsh
import sys
sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
b1 = 6
TOL = 1e-4  # relaxed for higher powers

verts, adj, lap = build_600cell()
N = 120

def qmul(p, q):
    w = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
    x = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
    y = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
    z = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]
    return np.array([w, x, y, z])

def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

# Find ALL fibrations (left + right actions of all order-10 elements)
def find_all_fibrations():
    fibrations = []
    sigs = set()
    order10 = []
    for i in range(N):
        g = verts[i]
        p = g.copy()
        ok = True
        for k in range(2, 11):
            p = qmul(p, g)
            if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6):
                ok = False; break
            if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6):
                ok = False
        if ok: order10.append(i)

    for gi in order10:
        g = verts[gi]
        for side in ['left', 'right']:
            used = set()
            fibs = []
            for s in range(N):
                if s in used: continue
                fib = []
                q = verts[s].copy()
                for k in range(10):
                    idx = find_idx(q, verts)
                    if idx < 0 or idx in used: break
                    fib.append(idx)
                    used.add(idx)
                    q = qmul(q, g) if side == 'left' else qmul(g, q)
                if len(fib) == 10:
                    fibs.append(tuple(sorted(fib)))
            if len(fibs) == 12 and len(used) == N:
                sig = frozenset(frozenset(f) for f in fibs)
                if sig not in sigs:
                    sigs.add(sig)
                    fibrations.append(fibs)
    return fibrations

print("Finding all Hopf fibrations...")
all_fibs = find_all_fibrations()
print(f"Found {len(all_fibs)} distinct fibrations\n")

# For each fibration, compute Box and its trace moments
print(f"{'Fib':>3s} {'Tr1':>10s} {'Tr2':>10s} {'Tr3':>12s} {'Tr4':>12s} {'Tr5':>14s} {'Tr6':>14s} {'ker':>4s}")
print("-" * 90)

expected = {
    1: 0,
    2: 7200,       # N * 60 = N * degree * a1
    3: 14400,       # N^2
    4: None,        # unknown
    5: 3628800,     # (2*a1)! = 10!
    6: None,        # unknown
}

all_match = True
trace_values = {n: [] for n in range(1, 7)}

for fi, fibs in enumerate(all_fibs):
    A_f = np.zeros((N, N))
    for fib in fibs:
        for i in fib:
            for j in fib:
                if i != j and adj[i,j] > 0.5:
                    A_f[i,j] = 1.0

    Box = b1 * A_f - adj
    evals = eigvalsh(Box)
    ker = int(np.sum(np.abs(evals) < 1e-6))

    traces = {}
    for n in range(1, 7):
        traces[n] = np.sum(evals**n)
        trace_values[n].append(traces[n])

    print(f"{fi+1:3d} {traces[1]:10.1f} {traces[2]:10.1f} {traces[3]:12.1f} "
          f"{traces[4]:12.1f} {traces[5]:14.1f} {traces[6]:14.1f} {ker:4d}")

    for n, exp_val in expected.items():
        if exp_val is not None and abs(traces[n] - exp_val) > TOL:
            all_match = False

print()
print("=" * 70)
print("STABILITY ANALYSIS")
print("=" * 70)

for n in range(1, 7):
    vals = trace_values[n]
    mn, mx = min(vals), max(vals)
    spread = mx - mn
    mean = np.mean(vals)
    stable = spread < TOL
    print(f"  Tr(Box^{n}): mean = {mean:14.2f},  spread = {spread:.2e},  "
          f"stable: {'YES' if stable else 'NO'}")

# Check specific identities
print()
mean2 = np.mean(trace_values[2])
mean3 = np.mean(trace_values[3])
mean5 = np.mean(trace_values[5])

print(f"  Tr(Box^2) = {mean2:.1f} = N * degree * a1 = {N*12*a1}? "
      f"{abs(mean2 - N*12*a1) < TOL}")
print(f"  Tr(Box^3) = {mean3:.1f} = N^2 = {N**2}? "
      f"{abs(mean3 - N**2) < TOL}")
import math
fact_2a1 = math.factorial(2*a1)
print(f"  Tr(Box^5) = {mean5:.1f} = (2*a1)! = {fact_2a1}? "
      f"{abs(mean5 - fact_2a1) < TOL}")

# Also check: Tr(Box^2)/N = degree*(degree-2) = 12*10 = 120?
# No: Tr(Box^2)/N = 60, not 120.
# Tr(Box^2)/N = 60 = degree*a1 = 12*5
# Or: Tr(Box^2) = N*degree*a1 = 120*12*5 = 7200

# Alternative: Tr(Box^2) via adjacency matrix identity
# Box = b1*A_f - A, so Box^2 = b1^2*Af^2 - 2*b1*Af*A + A^2
# Tr(A^2) = 2*edges = 1440
# Tr(Af^2) = 2*fiber_edges = 240
# Tr(Af*A) = sum_ij Af_ij * A_ij = 2*fiber_edges = 240 (since Af <= A)
# Tr(Box^2) = 36*240 - 12*240 + 1440 = 8640 - 2880 + 1440 = 7200 CHECK
print(f"\n  Algebraic verification of Tr(Box^2):")
print(f"    b1^2 * Tr(Af^2) = {b1**2} * 240 = {b1**2 * 240}")
print(f"    2*b1 * Tr(Af*A) = {2*b1} * 240 = {2*b1 * 240}")
print(f"    Tr(A^2)          = 1440")
print(f"    Total: {b1**2*240} - {2*b1*240} + 1440 = {b1**2*240 - 2*b1*240 + 1440}")
print(f"    = 7200 = N*degree*a1. THEOREM (fibration-independent).")

# For Tr(Box^3): need Tr(Af^3)
# Af is disjoint C10 cycles. C10 has no triangles (girth 10 > 3).
# So Tr(Af^3) = 0.
# Tr(Box^3) = b1^3*Tr(Af^3) - 3*b1^2*Tr(Af^2*A) + 3*b1*Tr(Af*A^2) - Tr(A^3)
# Tr(Af^3) = 0 (no fiber triangles)
# Tr(A^3) = 6 * (#triangles in 600-cell)
# Need: Tr(Af^2*A) and Tr(Af*A^2) -- these might be fibration-dependent!

# Check empirically
fibs0 = all_fibs[0]
A_f0 = np.zeros((N, N))
for fib in fibs0:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5:
                A_f0[i,j] = 1.0

tr_Af3 = np.trace(A_f0 @ A_f0 @ A_f0)
tr_Af2A = np.trace(A_f0 @ A_f0 @ adj)
tr_AfA2 = np.trace(A_f0 @ adj @ adj)
tr_A3 = np.trace(adj @ adj @ adj)

print(f"\n  Components of Tr(Box^3) [fibration 1]:")
print(f"    Tr(Af^3) = {tr_Af3:.1f} (expect 0, no fiber triangles)")
print(f"    Tr(Af^2*A) = {tr_Af2A:.1f}")
print(f"    Tr(Af*A^2) = {tr_AfA2:.1f}")
print(f"    Tr(A^3) = {tr_A3:.1f} = 6 * (#triangles) => triangles = {tr_A3/6:.0f}")

# Check Tr(Af^2*A) for another fibration
if len(all_fibs) > 1:
    fibs1 = all_fibs[1]
    A_f1 = np.zeros((N, N))
    for fib in fibs1:
        for i in fib:
            for j in fib:
                if i != j and adj[i,j] > 0.5:
                    A_f1[i,j] = 1.0

    tr_Af2A_1 = np.trace(A_f1 @ A_f1 @ adj)
    tr_AfA2_1 = np.trace(A_f1 @ adj @ adj)
    print(f"\n  Components of Tr(Box^3) [fibration 2]:")
    print(f"    Tr(Af^2*A) = {tr_Af2A_1:.1f} (was {tr_Af2A:.1f})")
    print(f"    Tr(Af*A^2) = {tr_AfA2_1:.1f} (was {tr_AfA2:.1f})")
    fib_indep = (abs(tr_Af2A - tr_Af2A_1) < TOL and
                 abs(tr_AfA2 - tr_AfA2_1) < TOL)
    print(f"    Fibration-independent? {fib_indep}")

# Algebraic prediction for Tr(Box^3) = N^2
predicted = b1**3 * tr_Af3 - 3*b1**2 * tr_Af2A + 3*b1 * tr_AfA2 - tr_A3
print(f"\n  Predicted Tr(Box^3) = {predicted:.1f}")
print(f"  Expected N^2 = {N**2}")
print(f"  Match: {abs(predicted - N**2) < TOL}")

print("\n" + "=" * 70)
if all_match:
    print("ALL TRACE IDENTITIES CONFIRMED ACROSS ALL FIBRATIONS")
else:
    print("SOME TRACE IDENTITIES FAILED")
print("=" * 70)
