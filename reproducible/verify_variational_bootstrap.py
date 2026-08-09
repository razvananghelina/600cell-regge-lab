"""
verify_variational_bootstrap.py — Theorem: variational bootstrap.

Verifies that Tr((c*A_fiber - A)^3) = 3*N_f*(c - 2) is linear in c,
and the unique solution of Tr(Box(c)^3) = N^2 is c = b1 = 6.

All computation in exact integer arithmetic (int64).
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5; b1 = 6; N = 120
PASS = 0; tests_run = 0; tests_pass = 0

def check(name, cond, detail=""):
    global tests_run, tests_pass, PASS
    tests_run += 1
    if cond: tests_pass += 1; print(f"  [PASS] {name}")
    else: PASS = 1; print(f"  [FAIL] {name}")
    if detail: print(f"         {detail}")

# Build
verts, adj, lap = build_600cell()
A = np.round(adj).astype(np.int64)

def qmul(p, q):
    return np.array([p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
        p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
        p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
        p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])
def find_idx(v, vs, tol=1e-6):
    dots = vs @ v; idx = np.argmax(dots)
    return idx if dots[idx] > 1-tol else -1

fibers = None
for i in range(N):
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
    for s in range(N):
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

Af = np.zeros((N, N), dtype=np.int64)
for i in range(N):
    for j in range(N):
        if A[i, j] == 1 and vtf[i] == vtf[j]:
            Af[i, j] = 1

# Exact traces (int64)
Af3 = Af @ Af @ Af
A3 = A @ A @ A
Af2A = (Af @ Af) @ A
AfA2 = Af @ (A @ A)
AfAAf = Af @ A @ Af

t_Af3 = int(np.trace(Af3))
t_Af2A = int(np.trace(Af2A))
t_AfAAf = int(np.trace(AfAAf))
t_AfA2 = int(np.trace(AfA2))
t_A2Af = int(np.trace((A @ A) @ Af))
t_AAfA = int(np.trace(A @ Af @ A))
t_A3 = int(np.trace(A3))

Nf = 1200
h = a1 * b1  # 30

print("=" * 60)
print("VARIATIONAL BOOTSTRAP VERIFICATION")
print("=" * 60)

# Polynomial coefficients
alpha_3 = t_Af3
alpha_2 = -(t_Af2A + t_AfAAf + int(np.trace(A @ Af @ Af)))
alpha_1 = t_AfA2 + t_AAfA + t_A2Af
alpha_0 = -t_A3

check("alpha_3 = Tr(A_f^3) = 0", alpha_3 == 0, f"got {alpha_3}")
check("alpha_2 = 0 (no 2-fiber-edge triangles)", alpha_2 == 0, f"got {alpha_2}")
check("alpha_1 = 3*N_f = 3600", alpha_1 == 3 * Nf, f"got {alpha_1}, expect {3*Nf}")
check("alpha_0 = -6*N_f = -7200", alpha_0 == -6 * Nf, f"got {alpha_0}, expect {-6*Nf}")
check("Polynomial is LINEAR (c^3=c^2=0)", alpha_3 == 0 and alpha_2 == 0)

# Unique solution
numerator = N * N - alpha_0  # N^2 + 6*N_f = 14400 + 7200 = 21600
denominator = alpha_1        # 3600
check("N^2 - alpha_0 divisible by alpha_1", numerator % denominator == 0,
      f"{numerator} / {denominator} = {numerator // denominator}")
c_sol = numerator // denominator
check("Unique solution c = b1 = 6", c_sol == b1, f"c = {c_sol}")

# Structural identities
check("c = N/h + 2 = 120/30 + 2 = 6", N // h + 2 == b1)
check("c = (a1-1) + 2 = 6 (arithmetic identity)", (a1 - 1) + 2 == b1)
check("c = a1 + 1 = b1", a1 + 1 == b1)
check("3*N_f/N = h = Coxeter number", 3 * Nf == N * h)

# Direct verification
Box = (b1 * Af - A).astype(np.int64)
tr3 = int(np.trace(Box @ Box @ Box))
check("Tr(Box(b1)^3) = N^2 = 14400 (exact)", tr3 == N * N, f"Tr = {tr3}")

# Triangle classification
from collections import defaultdict
adj_list = defaultdict(set)
for i in range(N):
    for j in range(N):
        if A[i,j]: adj_list[i].add(j)
T = {0: 0, 1: 0, 2: 0, 3: 0}
for i in range(N):
    for j in adj_list[i]:
        if j > i:
            for k in adj_list[i] & adj_list[j]:
                if k > j:
                    nf = int(Af[i,j]) + int(Af[j,k]) + int(Af[i,k])
                    T[nf] += 1
check("T_0 = T_1 = 600", T[0] == 600 and T[1] == 600, f"T_0={T[0]}, T_1={T[1]}")
check("T_2 = T_3 = 0 (Hopf geometry)", T[2] == 0 and T[3] == 0, f"T_2={T[2]}, T_3={T[3]}")
check("Total triangles N_f = 1200", T[0]+T[1]+T[2]+T[3] == Nf)

print(f"\n{'='*60}")
print(f"RESULTS: {tests_pass}/{tests_run} tests passed")
print(f"{'='*60}")
sys.exit(PASS)
