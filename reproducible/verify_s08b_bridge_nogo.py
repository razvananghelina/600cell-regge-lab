"""
verify_s08b_bridge_nogo.py
==========================

Focused verification of the S08b bridge obstruction.

This script checks, on the full six-fibration class:

  R1. The central element -1 acts on ker(Box_1) with fixed subspace of
      dimension exactly 1.
  R2. Therefore any quotient-compatible A5-equivariant map from the
      12-dimensional fiber module into ker(Box_1) can only land in the
      trivial 1-dimensional sector.
  R3. The canonical fiber-edge lift from R^12 to edge space has kernel
      projection of rank 1, not 12.

No claim is made here about arbitrary non-natural maps.
"""

import numpy as np
from scipy.linalg import eigh

from verify_s08_edge_fibration_uniformity import (
    a1,
    TOL,
    build_2I,
    build_adjacency,
    build_edges,
    build_line_graph_adj,
    build_group_actions,
    find_all_hopf_fibrations,
    build_fiber_edge_data,
)


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


def matrix_rank_tol(M, tol=1e-8):
    svals = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(svals > tol))


def main():
    print("=" * 72)
    print("VERIFY S08b BRIDGE NO-GO")
    print("=" * 72)

    verts = build_2I()
    A = build_adjacency(verts)
    edges, edge_to_idx = build_edges(A)
    A_line = build_line_graph_adj(len(verts), edges)
    all_vperm, all_eperm = build_group_actions(verts, edges, edge_to_idx)
    fibrations = find_all_hopf_fibrations(verts)

    minus1 = None
    for g in range(len(verts)):
        if np.allclose(verts[g], [-1, 0, 0, 0], atol=1e-6):
            minus1 = g
            break
    assert minus1 is not None

    all_fixed_dim1 = True
    all_lift_proj_rank1 = True
    all_lift_nonker = True

    print()
    for idx, fib_list in enumerate(fibrations, start=1):
        is_fiber, A_fiber_line = build_fiber_edge_data(fib_list, A, edges, edge_to_idx)
        A_cross_line = A_line - A_fiber_line
        L_fiber = np.diag(np.sum(A_fiber_line, axis=1)) - A_fiber_line
        L_cross = np.diag(np.sum(A_cross_line, axis=1)) - A_cross_line
        Box1 = L_cross - a1 * L_fiber

        evals_box, evecs_box = eigh(Box1)
        ker_mask = np.abs(evals_box) < TOL
        ker_basis = evecs_box[:, ker_mask]

        minus1_on_ker = ker_basis.T @ ker_basis[all_eperm[minus1], :]
        evals_m1 = np.sort(eigh(minus1_on_ker, eigvals_only=True))
        fixed_dim = int(np.sum(np.abs(evals_m1 - 1.0) < 1e-8))
        if fixed_dim != 1:
            all_fixed_dim1 = False

        # Canonical lift: one basis vector per fiber, constant on the 10 fiber
        # edges belonging to that fiber.
        lift = np.zeros((len(edges), 12), dtype=float)
        for fi, fiber in enumerate(fib_list):
            edge_ids = []
            for i in fiber:
                for j in fiber:
                    if i < j and A[i, j] > 0.5:
                        edge_ids.append(edge_to_idx[(i, j)])
            assert len(edge_ids) == 10
            lift[edge_ids, fi] = 1.0 / np.sqrt(10.0)

        box_norm = np.linalg.norm(Box1 @ lift, ord="fro")
        if box_norm < 1e-8:
            all_lift_nonker = False

        proj_to_ker = ker_basis @ (ker_basis.T @ lift)
        proj_rank = matrix_rank_tol(proj_to_ker, tol=1e-8)
        if proj_rank != 1:
            all_lift_proj_rank1 = False

        print(
            f"Fibration {idx}: "
            f"fixed_dim={fixed_dim}, "
            f"lift_Box_norm={box_norm:.8f}, "
            f"ker_proj_rank={proj_rank}"
        )

    print()
    check(
        all_fixed_dim1,
        "R1: (-1)-fixed subspace of ker(Box_1) has dimension 1 for all 6 fibrations",
    )
    check(
        all_lift_nonker,
        "R3a: the canonical fiber-edge lift is not itself contained in ker(Box_1)",
    )
    check(
        all_lift_proj_rank1,
        "R3b: the canonical fiber-edge lift has kernel projection rank 1 for all 6 fibrations",
    )

    print()
    print("Interpretation:")
    print("  Any quotient-compatible A5-equivariant map into ker(Box_1)")
    print("  must land in the (+1)-eigenspace of -1, hence at most in the")
    print("  unique trivial 1-dimensional sector.")
    print()
    print(f"TOTAL: {N_PASS} passed, {N_FAIL} failed")
    return 0 if N_FAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
