"""
verify_hopf_fibration_invariants.py
===================================

Focused verification of the discrete Hopf-fibration invariants used in S06.

This script verifies exactly the following four propositions on the full class
of discrete Hopf fibrations of the 600-cell arising from left cosets of
order-10 subgroups of 2I:

  P1. There are exactly 6 distinct such fibrations.
  P2. For every fibration F, the unique nontrivial kernel coefficient is c = 6.
  P3. For every F, ker(Box_F(6)) has dimension 9 and equals
      E_A(12) + E_A(6*phi) + E_A(6*phi').
  P4. For every F, lambda_1(L_cross(F)) / lambda_1(L_fiber(F)) = 5.

No physical interpretation is tested here.
"""

import numpy as np
from scipy.linalg import eigh
from itertools import permutations, product as cartesian_product


a1 = 5
b1 = 6
phi = (1.0 + np.sqrt(a1)) / 2.0
phi_conj = (1.0 - np.sqrt(a1)) / 2.0
TOL = 1e-8
EDGE_TOL = 1e-6

N_PASS = 0
N_FAIL = 0


def check(condition, label, detail=""):
    global N_PASS, N_FAIL
    if condition:
        N_PASS += 1
        print(f"[PASS] {label}")
    else:
        N_FAIL += 1
        print(f"[FAIL] {label}")
    if detail:
        print(f"       {detail}")


def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


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
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            add_vert(v)

    for signs in cartesian_product([0.5, -0.5], repeat=4):
        add_vert(list(signs))

    base = [0.0, 0.5, phi / 2.0, 1.0 / (2.0 * phi)]
    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nz_indices = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in cartesian_product([1, -1], repeat=len(nz_indices)):
            v = list(coords)
            for idx, s in zip(nz_indices, signs):
                v[idx] *= s
            add_vert(v)

    return np.array(sorted(verts))


def find_vertex_index(verts, q, tol=1e-6):
    dists = np.linalg.norm(verts - q, axis=1)
    idx = np.argmin(dists)
    if dists[idx] < tol:
        return idx
    return -1


def build_adjacency(verts):
    n = len(verts)
    dots = verts @ verts.T
    np.clip(dots, -1.0, 1.0, out=dots)
    A = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(dots[i, j] - phi / 2.0) < EDGE_TOL:
                A[i, j] = 1.0
                A[j, i] = 1.0
    return A


def find_all_hopf_fibrations(verts):
    n = len(verts)
    fibrations = []
    seen = set()

    order10_indices = []
    for i in range(n):
        g = verts[i]
        power = g.copy()
        found = False
        for k in range(2, 121):
            power = quat_mult(power, g)
            if np.allclose(power, [1, 0, 0, 0], atol=1e-6):
                if k == 10:
                    order10_indices.append(i)
                found = True
                break
        if not found:
            continue

    for gen_idx in order10_indices:
        g = verts[gen_idx]
        subgroup_indices = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        valid = True
        for _ in range(10):
            idx = find_vertex_index(verts, power)
            if idx < 0:
                valid = False
                break
            subgroup_indices.append(idx)
            power = quat_mult(power, g)
        if not valid or len(set(subgroup_indices)) != 10:
            continue

        assigned = np.full(n, -1, dtype=int)
        fib_list = []
        fid = 0
        for i in range(n):
            if assigned[i] >= 0:
                continue
            coset = []
            for si in subgroup_indices:
                q_prod = quat_mult(verts[i], verts[si])
                idx = find_vertex_index(verts, q_prod)
                if idx >= 0 and assigned[idx] < 0:
                    coset.append(idx)
                    assigned[idx] = fid
            if len(coset) != 10:
                valid = False
                break
            fib_list.append(tuple(sorted(coset)))
            fid += 1

        if not valid or fid != 12:
            continue

        signature = tuple(sorted(fib_list))
        if signature not in seen:
            seen.add(signature)
            fibrations.append(fib_list)

    return fibrations


def build_fiber_adjacency(A, fib_list):
    n = A.shape[0]
    Af = np.zeros((n, n), dtype=float)
    for fiber in fib_list:
        for i in fiber:
            for j in fiber:
                if i != j and A[i, j] > 0.5:
                    Af[i, j] = 1.0
    return Af


def grouped_eigenspace_basis(A):
    eigvals_A, eigvecs_A = eigh(A)
    tol_group = 0.1
    eigval_groups = {}
    for i, ev in enumerate(eigvals_A):
        found = False
        for key in eigval_groups:
            if abs(ev - key) < tol_group:
                eigval_groups[key].append(i)
                found = True
                break
        if not found:
            eigval_groups[ev] = [i]

    target_vals = [12.0, 6 * phi, 6 * phi_conj]
    target_indices = []
    for target in target_vals:
        for key in eigval_groups:
            if abs(key - target) < tol_group:
                target_indices.extend(eigval_groups[key])
                break
    return eigvecs_A[:, target_indices]


def main():
    print("=" * 72)
    print("VERIFY HOPF FIBRATION INVARIANTS")
    print("=" * 72)

    verts = build_2I()
    A = build_adjacency(verts)
    n = len(verts)
    target_space = grouped_eigenspace_basis(A)

    fibrations = find_all_hopf_fibrations(verts)
    n_found = len(fibrations)
    check(n_found == 6, "P1: exactly 6 distinct Hopf fibrations",
          f"Found {n_found}")

    all_unique_c6 = True
    all_kernel_dim9 = True
    all_kernel_sector = True
    all_gap_ratio5 = True

    print()
    for idx, fib_list in enumerate(fibrations, start=1):
        Af = build_fiber_adjacency(A, fib_list)
        Ac = A - Af
        Lf = 2.0 * np.eye(n) - Af
        Lc = 10.0 * np.eye(n) - Ac

        # P2: unique nontrivial kernel coefficient c = 6
        good_cs = []
        for c in range(1, 15):
            evals = eigh(c * Af - A, eigvals_only=True)
            kdim = int(np.sum(np.abs(evals) < TOL))
            if kdim > 0 and c != 1:
                good_cs.append(c)
        if good_cs != [6]:
            all_unique_c6 = False

        # P3: kernel dimension 9 and exact spectral sector
        evals_full, evecs_full = eigh(b1 * Af - A)
        kernel_mask = np.abs(evals_full) < TOL
        kernel_vecs = evecs_full[:, kernel_mask]
        if kernel_vecs.shape[1] != 9:
            all_kernel_dim9 = False
        proj = target_space.T @ kernel_vecs
        svals = np.linalg.svd(proj, compute_uv=False)
        if not (kernel_vecs.shape[1] == 9 and np.min(svals) > 0.99):
            all_kernel_sector = False

        # P4: gap ratio = 5
        evals_Lf = np.sort(eigh(Lf, eigvals_only=True))
        evals_Lc = np.sort(eigh(Lc, eigvals_only=True))
        gap_f = evals_Lf[evals_Lf > TOL][0]
        gap_c = evals_Lc[evals_Lc > TOL][0]
        ratio = gap_c / gap_f
        if abs(ratio - a1) > 1e-8:
            all_gap_ratio5 = False

        print(
            f"Fibration {idx}: c_list={good_cs}, "
            f"ker_dim={kernel_vecs.shape[1]}, "
            f"min_proj_sv={np.min(svals):.12f}, "
            f"gap_ratio={ratio:.12f}"
        )

    print()
    check(all_unique_c6, "P2: unique nontrivial kernel coefficient is c = 6")
    check(all_kernel_dim9, "P3a: dim ker(Box_F(6)) = 9 for all 6 fibrations")
    check(
        all_kernel_sector,
        "P3b: ker(Box_F(6)) = E_A(12) + E_A(6phi) + E_A(6phi') for all 6 fibrations",
    )
    check(all_gap_ratio5, "P4: lambda_1(L_cross)/lambda_1(L_fiber) = 5 for all 6 fibrations")

    print()
    print(f"TOTAL: {N_PASS} passed, {N_FAIL} failed")
    return 0 if N_FAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
