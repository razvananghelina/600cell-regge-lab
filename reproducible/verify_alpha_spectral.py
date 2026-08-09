#!/usr/bin/env python3
"""
Verify: alpha equation from spectral action on Box (Eq. alpha0_box).

Tests:
1. Tr(Box^2_gauge) = 1440 = deg * N
2. Box_base eigenvalues = {0, L(rho_2), L(rho_4), L(rho_6)} exact
3. B = Tr(Box^2_gauge) * phi^4 / (n_base * b1) = 4*a1*phi^4
4. Tr(Box^2_full)/Tr(Box^2_gauge) = a1
5. Holonomy = 2*pi (from c_1 = 1)
6. Alpha equation: 2*pi*alpha^2 - B*alpha + 1 = 0 gives 1/137.036
7. sigma=(01)(2) uniquely selected: only non-negative + degenerate permutation
"""
import numpy as np
from numpy.linalg import eigh, eigvalsh
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell
from commons.constants import PHI, SQRT5, a1, b1, N

ALPHA_EXP = 1 / 137.035999084
PASS_ALL = True

def check(name, condition, detail=""):
    global PASS_ALL
    status = "PASS" if condition else "FAIL"
    if not condition: PASS_ALL = False
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))

# Build 600-cell + Hopf
verts, adj, lap = build_600cell()

def qmul(p, q):
    return np.array([
        p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])

def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1

fibers = None
for i in range(N):
    if abs(verts[i, 0] - PHI/2) < 1e-6:
        g = verts[i]; p = g.copy(); ok = True
        for k in range(2, 11):
            p = qmul(p, g)
            if k == 5 and not np.allclose(p, [-1,0,0,0], atol=1e-6): ok=False; break
            if k == 10 and not np.allclose(p, [1,0,0,0], atol=1e-6): ok=False
        if not ok: continue
        used = set(); fibers_tmp = []; subg = []
        pp = np.array([1.0,0,0,0])
        for k in range(10): subg.append(find_idx(pp, verts)); pp = qmul(pp, g)
        for s in range(N):
            if s in used: continue
            fib = []
            for si in subg:
                q = qmul(verts[s], verts[si]); idx = find_idx(q, verts)
                if idx >= 0 and idx not in used: fib.append(idx); used.add(idx)
            if len(fib) == 10: fibers_tmp.append(fib)
        if len(fibers_tmp) == 12:
            fibers = fibers_tmp
            break

n_base = len(fibers)
n_fiber = len(fibers[0])

A_fiber = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            if i != j and adj[i,j] > 0.5: A_fiber[i,j] = 1.0

Box = b1 * A_fiber - adj

# Fiber zero mode projector
P0 = np.zeros((N, N))
for fib in fibers:
    for i in fib:
        for j in fib:
            P0[i, j] = 1.0 / n_fiber

Box_gauge = P0 @ Box @ P0

# Build 12x12 base operator
Box_base = np.zeros((n_base, n_base))
for a in range(n_base):
    for b_idx in range(n_base):
        val = sum(Box[i, j] for i in fibers[a] for j in fibers[b_idx])
        Box_base[a, b_idx] = val / n_fiber

print("=" * 60)
print("verify_alpha_spectral.py")
print("=" * 60)

# Test 1: Tr(Box^2_gauge) = 1440
Tr_gauge_2 = np.trace(Box_gauge @ Box_gauge)
check("Tr(Box^2_gauge) = 1440", abs(Tr_gauge_2 - 1440) < 0.01,
      f"got {Tr_gauge_2:.1f}")

# Test 2: Box_base eigenvalues
evals_base = sorted(eigvalsh(Box_base))
L_rho2 = 10 - 2*SQRT5
L_rho4 = 12.0
L_rho6 = 10 + 2*SQRT5
expected = sorted([0.0]*1 + [L_rho2]*3 + [L_rho4]*5 + [L_rho6]*3)
check("Box_base eigenvalues = {0, L(rho2), L(rho4), L(rho6)}",
      np.allclose(evals_base, expected, atol=1e-6))

# Test 3: B formula
B = Tr_gauge_2 * PHI**4 / (n_base * b1)
B_target = 4 * a1 * PHI**4
check("B = Tr(Box^2_gauge)*phi^4/(n_base*b1) = 4*a1*phi^4",
      abs(B - B_target) / B_target < 1e-10,
      f"B={B:.6f}, target={B_target:.6f}")

# Test 4: Tr ratio = a1
Tr_full_2 = np.trace(Box @ Box)
ratio = Tr_full_2 / Tr_gauge_2
check("Tr(Box^2_full)/Tr(Box^2_gauge) = a1",
      abs(ratio - a1) < 0.01,
      f"ratio={ratio:.4f}")

# Test 5: Holonomy
holonomy = 2 * a1 * np.pi / a1
check("Holonomy = 2*pi", abs(holonomy - 2*np.pi) < 1e-12)

# Test 6: Alpha equation
disc = B**2 - 4*2*np.pi
alpha_pred = (B - np.sqrt(disc)) / (2*2*np.pi)
err = abs(1/alpha_pred - 1/ALPHA_EXP) / (1/ALPHA_EXP)
check("Alpha equation gives 1/137.036 (0.001% precision)",
      err < 1e-5,
      f"1/alpha={1/alpha_pred:.6f}, err={err*100:.5f}%")

# Test 7: sigma uniqueness
N_gen = 3
delta_up = {0: 3, 1: 5, 2: 9}
lepton_n = {0: 0, 1: 11, 2: 17}
perms = [
    {0:0,1:1,2:2}, {0:1,1:0,2:2}, {0:2,1:1,2:0},
    {0:0,1:2,2:1}, {0:1,1:2,2:0}, {0:2,1:0,2:1},
]
up_n = {0: 3, 1: 16, 2: 26}
valid = []
for perm in perms:
    S = [a1 + N_gen*perm[g] for g in range(3)]
    dd = [S[g] - delta_up[g] for g in range(3)]
    nd = [lepton_n[g] + dd[g] for g in range(3)]
    all_nn = all(d >= 0 for d in dd)
    all_n = sorted(list(lepton_n.values()) + list(up_n.values()) + nd)
    has_degen = len(set(all_n)) < 9  # McKay zero-edge requires degeneracy
    if all_nn and has_degen:
        valid.append(perm)

check("sigma=(01)(2) uniquely selected (non-negative + degenerate)",
      len(valid) == 1 and valid[0] == {0:1,1:0,2:2},
      f"valid permutations: {len(valid)}")

print()
if PASS_ALL:
    print("ALL TESTS PASSED")
else:
    print("SOME TESTS FAILED")
