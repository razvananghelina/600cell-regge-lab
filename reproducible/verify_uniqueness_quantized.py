"""
Verifier for a strengthened uniqueness route:

    C1  : each z_f is 0 / unit / irreducible in Z[phi]
    Q   : on each McKay edge, dz / sigma(dz) is +/- phi^k for some integer k
    C2  : sum_edges N(dz) = -b1

This script does not prove the theorem analytically.  It gives a clean,
reproducible finite verification of the stronger structural formulation that
seems substantially better than the current C1+C2-only statement.
"""

from __future__ import annotations

import math
from functools import lru_cache


A1 = 5
B1 = 6
PHI = (1 + math.sqrt(5)) / 2
PHIP = (1 - math.sqrt(5)) / 2

# Fermion order matching the exploratory McKay-tree scripts:
# 0:e, 1:u, 2:d, 3:s, 4:mu, 5:c, 6:tau, 7:branch-A, 8:branch-B
# The last two branch labels are intentionally left neutral here because the
# older scripts swap the physical t/b labels, while the arithmetic selection
# itself is unaffected.
N_EXPONENTS = [0, 3, 5, 11, 11, 16, 17, 19, 26]
NAMES = ["e", "u", "d", "s", "mu", "c", "tau", "brA", "brB"]
PHYSICAL_AB = [
    (0, 0),
    (3, -2),
    (1, 0),
    (1, 1),
    (1, 1),
    (2, 1),
    (1, 2),
    (-1, 4),
    (4, 1),
]

# McKay tree edges in the fermion indexing used above.
EDGES = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (6, 7),
    (5, 8),
]


def ab_from_t(n: int, t: int) -> tuple[int, int]:
    # General solution of 5a + 6b = n
    return -n - 6 * t, n + 5 * t


def norm_zphi(a: int, b: int) -> int:
    return a * a + a * b - b * b


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = math.isqrt(n)
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


def c1_irreducible(abs_norm: int) -> bool:
    if abs_norm <= 1:
        return True
    return is_prime(abs_norm) and abs_norm % 5 in (0, 1, 4)


def ratio_is_phi_power(da: int, db: int, tol: float = 0.01) -> tuple[bool, int | None, int]:
    """
    Check whether (da + db*phi)/(da + db*phi') = +/- phi^k for some integer k.
    Return (ok, k, sign), where sign in {+1,-1}.
    """
    dz = da + db * PHI
    dzg = da + db * PHIP
    if abs(dz) < 1e-12 and abs(dzg) < 1e-12:
        return True, 0, 1
    if abs(dzg) < 1e-12:
        return True, 999, 1
    if abs(dz) < 1e-12:
        return True, -999, 1

    ratio = dz / dzg
    sgn = 1 if ratio > 0 else -1
    log_phi = math.log(abs(ratio)) / math.log(PHI)
    n = round(log_phi)
    if abs(log_phi - n) < tol:
        return True, int(n), sgn
    return False, None, sgn


def edge_norm_sum(assign: list[tuple[int, int]]) -> int:
    total = 0
    for i, j in EDGES:
        da = assign[j][0] - assign[i][0]
        db = assign[j][1] - assign[i][1]
        total += norm_zphi(da, db)
    return total


def all_edges_quantized(assign: list[tuple[int, int]]) -> tuple[bool, list[tuple[tuple[int, int], int, int]]]:
    data = []
    for i, j in EDGES:
        da = assign[j][0] - assign[i][0]
        db = assign[j][1] - assign[i][1]
        ok, k, sgn = ratio_is_phi_power(da, db)
        if not ok:
            return False, []
        data.append(((i, j), k, sgn))
    return True, data


def physical_t_values() -> list[int]:
    ts = []
    for n, (a, b) in zip(N_EXPONENTS, PHYSICAL_AB):
        # Solve a = -n - 6t, b = n + 5t
        if n == 0:
            ts.append(0)
            continue
        t_from_a = (-n - a) / 6
        t_from_b = (b - n) / 5
        assert abs(t_from_a - t_from_b) < 1e-12
        ts.append(int(round(t_from_a)))
    return ts


print("=" * 78)
print("VERIFY STRENGTHENED UNIQUENESS: C1 + EDGE QUANTIZATION + C2")
print("=" * 78)

phys_t = physical_t_values()
print("\nPhysical t-values:")
for name, n, t in zip(NAMES, N_EXPONENTS, phys_t):
    print(f"  {name:>3s}: n={n:>2d}, t={t:>3d}")

print("\nStep 1: Enumerate C1-admissible t-values in a bounded window |t|<=15")
valid_t: dict[int, list[int]] = {}
for idx in range(1, 9):  # skip electron, fixed at t=0
    n = N_EXPONENTS[idx]
    allowed = []
    for t in range(-15, 16):
        a, b = ab_from_t(n, t)
        absN = abs(norm_zphi(a, b))
        if c1_irreducible(absN):
            allowed.append(t)
    valid_t[idx] = allowed
    print(f"  {NAMES[idx]:>3s}: {len(allowed):>2d} values, physical t={phys_t[idx]:>3d}, present={phys_t[idx] in allowed}")

count_c1 = 1
for idx in range(1, 9):
    count_c1 *= len(valid_t[idx])
print(f"\n  Total C1 combinations in |t|<=15: {count_c1:,}")

candidate_data: dict[int, list[tuple[int, tuple[int, int]]]] = {0: [(0, (0, 0))]}
for idx in range(1, 9):
    candidate_data[idx] = [(t, ab_from_t(N_EXPONENTS[idx], t)) for t in valid_t[idx]]


def edge_info(left: tuple[int, int], right: tuple[int, int]) -> tuple[bool, int, int, int]:
    da = right[0] - left[0]
    db = right[1] - left[1]
    ok, k, sgn = ratio_is_phi_power(da, db)
    return ok, (k if k is not None else 0), sgn, norm_zphi(da, db)


# Compatibility lists for every directed edge in the tree.
compat: dict[tuple[int, int], dict[int, list[tuple[int, int, int, int]]]] = {}
for i, j in EDGES + [(b, a) for a, b in EDGES]:
    compat[(i, j)] = {}
    for ti, abi in candidate_data[i]:
        lst = []
        for tj, abj in candidate_data[j]:
            ok, k, sgn, eN = edge_info(abi, abj)
            if ok:
                lst.append((tj, k, sgn, eN))
        compat[(i, j)][ti] = lst


print("\nStep 2: Filter by edge Galois quantization Q (tree search with pruning)")

# Root the tree exactly as in the exploratory scripts.
children = {
    0: [1],
    1: [2],
    2: [3],
    3: [4],
    4: [5],
    5: [6, 8],
    6: [7],
    7: [],
    8: [],
}


@lru_cache(maxsize=None)
def subtree_extensions(node: int, t_node: int) -> tuple[tuple[tuple[tuple[int, ...], int], ...], ...]:
    """
    Return all subtree completions below `node`, represented as tuples:
        ((t0,...,t8), edge_sum_subtree)
    where the assignment tuple contains fixed t-values only on nodes in the subtree.
    """
    if not children[node]:
        base = [None] * 9
        base[node] = t_node
        return (((tuple(-999 if x is None else x for x in base), 0),),)

    # Start with the trivial partial assignment for this node.
    partials = [( [None] * 9, 0 )]
    partials[0][0][node] = t_node

    for child in children[node]:
        new_partials = []
        for t_child, _k, _sgn, edge_n in compat[(node, child)][t_node]:
            child_ext_groups = subtree_extensions(child, t_child)
            for group in child_ext_groups:
                for child_assign, child_sum in group:
                    child_list = [None if x == -999 else x for x in child_assign]
                    for cur_assign, cur_sum in partials:
                        merged = cur_assign[:]
                        ok_merge = True
                        for idx, val in enumerate(child_list):
                            if val is None:
                                continue
                            if merged[idx] is not None and merged[idx] != val:
                                ok_merge = False
                                break
                            merged[idx] = val
                        if ok_merge:
                            new_partials.append((merged, cur_sum + edge_n + child_sum))
        partials = new_partials

    packed = tuple((tuple(-999 if x is None else x for x in assign), s) for assign, s in partials)
    return (packed,)


quantized = []
for t0, _ab0 in candidate_data[0]:
    for group in subtree_extensions(0, t0):
        for t_assign, partial_edge_sum in group:
            if -999 in t_assign:
                continue
            assign = [ab_from_t(n, t) for n, t in zip(N_EXPONENTS, t_assign)]
            ok, edge_data = all_edges_quantized(assign)
            if ok:
                quantized.append((list(t_assign), assign, edge_data, partial_edge_sum))

print(f"  Surviving combinations after Q: {len(quantized)}")

print("\nStep 3: Filter by flatness C2: sum_edges N(dz) = -b1 = -6")
c2_solutions = []
for t_full, assign, edge_data, s in quantized:
    if s == -B1:
        c2_solutions.append((t_full, assign, edge_data))

print(f"  Surviving combinations after C2: {len(c2_solutions)}")

physical_assign = PHYSICAL_AB
physical_found = any(sol[1] == physical_assign for sol in c2_solutions)
print(f"\nPhysical assignment present among final solutions: {physical_found}")

if c2_solutions:
    print("\nFinal solution(s):")
    for k, (t_full, assign, edge_data) in enumerate(c2_solutions, start=1):
        print(f"  Solution #{k}:")
        print(f"    t = {t_full}")
        print(f"    (a,b) = {assign}")
        print(f"    edge_sum = {edge_norm_sum(assign)}")
        edge_txt = ", ".join(
            f"({NAMES[i]},{NAMES[j]}): {'+' if sgn > 0 else '-'}phi^{pow_k}"
            for (i, j), pow_k, sgn in edge_data
        )
        print(f"    edge ratios: {edge_txt}")

print("\nStep 4: Compare with the physical assignment")
phys_assign = [ab_from_t(n, t) for n, t in zip(N_EXPONENTS, phys_t)]
phys_quantized, phys_edge_data = all_edges_quantized(phys_assign)
print(f"  Physical assignment edge-quantized: {phys_quantized}")
print(f"  Physical assignment edge sum: {edge_norm_sum(phys_assign)}")
if phys_quantized:
    edge_txt = ", ".join(
        f"({NAMES[i]},{NAMES[j]}): {'+' if sgn > 0 else '-'}phi^{pow_k}"
        for (i, j), pow_k, sgn in phys_edge_data
    )
    print(f"  Physical edge ratios: {edge_txt}")

print("\n" + "=" * 78)
if len(c2_solutions) == 1 and physical_found:
    print("RESULT: PASS")
    print("Within the bounded search |t|<=15, the strengthened criterion C1+Q+C2")
    print("selects a unique assignment, and it is the physical one.")
else:
    print("RESULT: WARNING")
    print("Either uniqueness failed in the bounded search or the physical solution was not singled out.")
print("=" * 78)
