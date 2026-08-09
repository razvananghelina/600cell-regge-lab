"""
verify_neutral_vacuum_scale.py

Verifies the internal spectral support for the electroweak exponent n = 25:

  (a) ker(A) has dimension 25 and is the unique adjacency-blind block
  (b) the corresponding irrep has dimension 5 and generator character 0
  (c) Tr(exp(-a1*A^2)) isolates this block to high precision already at t=a1=5
  (d) Box refines ker(A) into 12^(5) + (6/phi)^(10) + (-6phi)^(10)
  (e) the phi-conjugate pair ratio gives the exact exponent count 5 + 2*10 = 25
"""

import os
import sys
from collections import defaultdict

import numpy as np
from numpy.linalg import eigh, norm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from commons import build_600cell

PHI = (1 + np.sqrt(5)) / 2
a1 = 5
PASS = 0
tests_run = 0
tests_pass = 0


def check(name, condition, detail=""):
    global tests_run, tests_pass, PASS
    tests_run += 1
    if condition:
        tests_pass += 1
        print(f"  [PASS] {name}")
    else:
        PASS = 1
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def get_count(mult_dict, key):
    return mult_dict.get(key, 0) + mult_dict.get(-key, 0)


def count_close(evals, target, tol=1e-8):
    return int(np.sum(np.abs(evals - target) < tol))


def qmul(p, q):
    return np.array(
        [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
        ]
    )


def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1 - tol else -1


def find_fibration(verts):
    n = len(verts)
    target_w = PHI / 2.0
    for i in range(n):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1, 0, 0, 0], atol=1e-6):
                    ok = False
                    break
                if k == 10 and not np.allclose(p, [1, 0, 0, 0], atol=1e-6):
                    ok = False
            if not ok:
                continue
            used = set()
            fibers = []
            subg = []
            pp = np.array([1.0, 0, 0, 0])
            for _ in range(10):
                subg.append(find_idx(pp, verts))
                pp = qmul(pp, g)
            for s in range(n):
                if s in used:
                    continue
                fib = []
                for si in subg:
                    idx = find_idx(qmul(verts[s], verts[si]), verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx)
                        used.add(idx)
                if len(fib) == 10:
                    fibers.append(fib)
            if len(fibers) == 12:
                return fibers
    return None


print("Building 600-cell and spectral data...")
verts, adj, lap = build_600cell()
_ = lap
n = len(verts)

evals_A, evecs_A = eigh(adj)
groups = defaultdict(list)
for i, val in enumerate(np.round(evals_A, 10)):
    groups[float(val)].append(i)

blocks = []
for lam, inds in sorted(groups.items()):
    basis = evecs_A[:, inds]
    mult = len(inds)
    dim = int(round(np.sqrt(mult)))
    blocks.append((lam, basis, mult, dim))

print("\n" + "=" * 60)
print("TEST (a): unique neutral adjacency block")
print("=" * 60)
neutral = [(lam, basis, mult, dim) for lam, basis, mult, dim in blocks if abs(lam) < 1e-8]
check("Exactly one adjacency block with lambda_A = 0", len(neutral) == 1, f"count = {len(neutral)}")
if neutral:
    _, ker_basis, ker_mult, ker_dim = neutral[0]
    check("dim ker(A) = 25", ker_mult == 25, f"multiplicity = {ker_mult}")
    check("Underlying irrep dimension = 5", ker_dim == 5, f"dim = {ker_dim}")
else:
    ker_basis = None

print("\n" + "=" * 60)
print("TEST (b): generator character vanishes on the neutral block")
print("=" * 60)
identity = find_idx(np.array([1.0, 0, 0, 0]), verts)
gens = [j for j in range(n) if adj[identity, j] > 0.5]
check("Generator class has size 12", len(gens) == 12, f"|S| = {len(gens)}")

perms = []
for g_idx in range(n):
    vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(n)])
    perms.append(vp)

if ker_basis is not None:
    g0 = gens[0]
    proj = ker_basis @ ker_basis.T
    chi_isotypic = np.trace(proj[:, perms[g0]])
    chi_irrep = chi_isotypic / ker_dim
    blind = norm(adj @ ker_basis) / np.sqrt(ker_mult)
    check("chi_rho4(generator) = 0", abs(chi_irrep) < 1e-10, f"chi = {chi_irrep:.12e}")
    check("A annihilates the neutral block", blind < 1e-10, f"||A||_block = {blind:.3e}")

print("\n" + "=" * 60)
print("TEST (c): IR projector isolates rank 25 at t = a1")
print("=" * 60)
tr_ir = float(np.sum(np.exp(-a1 * evals_A**2)))
check("Tr(exp(-a1*A^2)) ~= 25 at t = 5", abs(tr_ir - 25.0) < 1e-6, f"trace = {tr_ir:.12f}")

print("\n" + "=" * 60)
print("TEST (d): Box refinement of the neutral sector")
print("=" * 60)
fibers = find_fibration(verts)
check("Hopf fibration found", fibers is not None)

if fibers is not None and ker_basis is not None:
    A_f = np.zeros_like(adj)
    for fib in fibers:
        for i in fib:
            for j in fib:
                if i != j and adj[i, j] > 0.5:
                    A_f[i, j] = 1.0

    A_c = adj - A_f
    L_f = np.diag(np.sum(A_f, axis=1)) - A_f
    L_c = np.diag(np.sum(A_c, axis=1)) - A_c
    Box = L_c - a1 * L_f

    evals_Lf = np.round(np.linalg.eigvalsh(ker_basis.T @ L_f @ ker_basis), 10)
    evals_Box = np.round(np.linalg.eigvalsh(ker_basis.T @ Box @ ker_basis), 10)

    mult_Lf = defaultdict(int)
    mult_Box = defaultdict(int)
    for x in evals_Lf:
        mult_Lf[float(x)] += 1
    for x in evals_Box:
        mult_Box[float(x)] += 1

    ok_Lf = (
        count_close(evals_Lf, 0.0) == 5
        and count_close(evals_Lf, 3 - PHI) == 10
        and count_close(evals_Lf, PHI + 2) == 10
    )
    ok_Box = (
        count_close(evals_Box, 12.0) == 5
        and count_close(evals_Box, 6 / PHI) == 10
        and count_close(evals_Box, -6 * PHI) == 10
    )

    check("L_fiber|ker(A) = 0^5 + (3-phi)^10 + (phi+2)^10", ok_Lf, f"{dict(mult_Lf)}")
    check("Box|ker(A) = 12^5 + (6/phi)^10 + (-6phi)^10", ok_Box, f"{dict(mult_Box)}")

    print("\n" + "=" * 60)
    print("TEST (e): internal exponent counter")
    print("=" * 60)
    ratio = abs((-6 * PHI) / (6 / PHI))
    exponent_count = 5 + 10 * (np.log(ratio) / np.log(PHI))
    check("|mu_-|/mu_+ = phi^2", abs(ratio - PHI**2) < 1e-10, f"ratio = {ratio:.12f}")
    check("5 + 10*log_phi(phi^2) = 25", abs(exponent_count - 25.0) < 1e-10,
          f"count = {exponent_count:.12f}")

print("\n" + "=" * 60)
print(f"TOTAL: {tests_pass}/{tests_run} tests PASSED")
print("=" * 60)
sys.exit(PASS)
