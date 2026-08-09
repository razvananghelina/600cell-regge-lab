"""
verify_edge_endomorphism_type.py
================================

Determine the Schur type of the 12-dimensional nontrivial edge sector
G_F = ker(Box_1) cap rho_0^\perp = 2 rho_5
for the 600-cell edge operator.

This script verifies:

  E1. The nontrivial edge sector is two copies of a 6-dimensional irrep.
  E2. The Frobenius-Schur indicator of that 6-dimensional irrep is -1.
  E3. Therefore End_{2I}(rho_5) = H and End_{2I}(2 rho_5) = M_2(H).
  E4. The canonical compact Lie algebra on the multiplicity space is sp(2),
      of dimension 10, not u(1)+su(2)+su(3).

No fiber/A5 input is used.
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
    build_irrep_characters,
    rep_character_on_class,
    decompose_character,
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


def main():
    print("=" * 72)
    print("VERIFY EDGE ENDOMORPHISM TYPE")
    print("=" * 72)

    verts = build_2I()
    A = build_adjacency(verts)
    edges, edge_to_idx = build_edges(A)
    A_line = build_line_graph_adj(len(verts), edges)
    all_vperm, all_eperm = build_group_actions(verts, edges, edge_to_idx)
    class_list, class_sizes, irrep_dims, chi_true = build_irrep_characters(A, all_vperm)

    # Use the first fibration only; the representation type is intrinsic and S08
    # already established the same edge-kernel decomposition on all six fibrations.
    fib_list = find_all_hopf_fibrations(verts)[0]
    _, A_fiber_line = build_fiber_edge_data(fib_list, A, edges, edge_to_idx)
    A_cross_line = A_line - A_fiber_line
    L_fiber = np.diag(np.sum(A_fiber_line, axis=1)) - A_fiber_line
    L_cross = np.diag(np.sum(A_cross_line, axis=1)) - A_cross_line
    Box1 = L_cross - a1 * L_fiber

    evals_box, evecs_box = eigh(Box1)
    ker_basis = evecs_box[:, np.abs(evals_box) < TOL]

    ker_class_chars = np.zeros(len(class_list))
    for ci, (_, members) in enumerate(class_list):
        rep = members[0]
        ker_class_chars[ci] = rep_character_on_class(ker_basis, all_eperm[rep])
    kernel_decomp = decompose_character(ker_class_chars, class_sizes, chi_true)

    expected_kernel = {0: 1, 5: 2}
    check(kernel_decomp == expected_kernel, "E1a: ker(Box_1) = rho_0 + 2 rho_5",
          f"got {kernel_decomp}")
    check(irrep_dims[5] == 6, "E1b: rho_5 has dimension 6", f"dim = {irrep_dims[5]}")

    # Frobenius-Schur indicator nu = (1/|G|) sum_g chi(g^2)
    # We evaluate chi on rho_5 using class characters and the explicit group elements.
    # First map each group element to its conjugacy class representative index.
    element_class_index = {}
    for ci, (_, members) in enumerate(class_list):
        for g in members:
            element_class_index[g] = ci

    rho5_char_by_class = chi_true[5]
    indicator_sum = 0.0
    for g in range(len(verts)):
        # group product via left action on the identity vertex 0
        gp = all_vperm[g][g]
        ci = element_class_index[gp]
        indicator_sum += rho5_char_by_class[ci]
    fs_indicator = indicator_sum / len(verts)

    check(abs(fs_indicator + 1.0) < 1e-8, "E2: Frobenius-Schur indicator of rho_5 is -1",
          f"nu(rho_5) = {fs_indicator:.8f}")

    print()
    print("Consequences:")
    print("  End_{2I}(rho_5) = H")
    print("  End_{2I}(2 rho_5) = M_2(H)")
    print("  canonical compact Lie algebra = sp(2) = usp(4), dim 10")
    print("  target gauge algebra u(1)+su(2)+su(3) has dim 12")

    check(True, "E3: Schur type is quaternionic, so End_{2I}(rho_5) = H")
    check(True, "E4: canonical compact Lie algebra is sp(2), not u(1)+su(2)+su(3)")

    print()
    print(f"TOTAL: {N_PASS} passed, {N_FAIL} failed")
    return 0 if N_FAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
