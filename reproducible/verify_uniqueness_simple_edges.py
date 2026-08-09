"""
Verifier for an even stronger uniqueness principle than C1+Q+C2.

Constraint set:
  C1: each node value z_f is 0 / unit / irreducible in Z[phi]
  S : each McKay-edge difference dz lies in
        {0} union { +/- phi^r } union { +/- 2 phi^r }, r in Z

Empirical result in the bounded domain |t|<=15:
  C1 + S already singles out the physical assignment uniquely.

This suggests that the true primitive principle may live on edge differences
themselves, with Q (Galois-ratio quantization) following as a consequence.
"""

from __future__ import annotations

import math
from functools import lru_cache


PHI = (1 + math.sqrt(5)) / 2

NAMES = ["e", "u", "d", "s", "mu", "c", "tau", "brA", "brB"]
N_EXPONENTS = [0, 3, 5, 11, 11, 16, 17, 19, 26]
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
EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (5, 8)]


def ab_from_t(n: int, t: int) -> tuple[int, int]:
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


def simple_edge_class(da: int, db: int) -> tuple[bool, tuple[int, int] | None]:
    """
    Check if dz = da + db*phi lies in {0} U {+/- phi^r} U {+/- 2 phi^r}.
    Return (ok, (coeff, power)) where coeff is in {0, +/-1, +/-2}.
    """
    dz = da + db * PHI
    if abs(dz) < 1e-12:
        return True, (0, 0)

    for coeff in (1, 2):
        for power in range(-12, 13):
            for sign in (-1, 1):
                if abs(dz - sign * coeff * (PHI**power)) < 1e-9:
                    return True, (sign * coeff, power)
    return False, None


def physical_t_values() -> list[int]:
    ts = []
    for n, (a, b) in zip(N_EXPONENTS, PHYSICAL_AB):
        if n == 0:
            ts.append(0)
            continue
        t_from_a = (-n - a) / 6
        t_from_b = (b - n) / 5
        assert abs(t_from_a - t_from_b) < 1e-12
        ts.append(int(round(t_from_a)))
    return ts


valid_t: dict[int, list[int]] = {0: [0]}
for idx in range(1, 9):
    vals = []
    n = N_EXPONENTS[idx]
    for t in range(-15, 16):
        a, b = ab_from_t(n, t)
        if c1_irreducible(abs(norm_zphi(a, b))):
            vals.append(t)
    valid_t[idx] = vals


candidate_data = {idx: [(t, ab_from_t(N_EXPONENTS[idx], t)) for t in valid_t[idx]] for idx in range(9)}

compat: dict[tuple[int, int], dict[int, list[tuple[int, int, tuple[int, int]]]]] = {}
for i, j in EDGES + [(b, a) for a, b in EDGES]:
    compat[(i, j)] = {}
    for ti, abi in candidate_data[i]:
        lst = []
        for tj, abj in candidate_data[j]:
            da = abj[0] - abi[0]
            db = abj[1] - abi[1]
            ok, tag = simple_edge_class(da, db)
            if ok:
                lst.append((tj, norm_zphi(da, db), tag))
        compat[(i, j)][ti] = lst


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
def subtree_extensions(node: int, t_node: int) -> tuple[tuple[int | None, ...], ...]:
    if not children[node]:
        base = [None] * 9
        base[node] = t_node
        return (tuple(base),)

    partials = [[None] * 9]
    partials[0][node] = t_node

    for child in children[node]:
        new_partials = []
        for t_child, _edge_norm, _tag in compat[(node, child)][t_node]:
            for child_assign in subtree_extensions(child, t_child):
                child_list = list(child_assign)
                for cur_assign in partials:
                    merged = cur_assign[:]
                    ok = True
                    for idx, val in enumerate(child_list):
                        if val is None:
                            continue
                        if merged[idx] is not None and merged[idx] != val:
                            ok = False
                            break
                        merged[idx] = val
                    if ok:
                        new_partials.append(merged)
        partials = new_partials

    return tuple(tuple(x) for x in partials)


solutions: list[tuple[tuple[int, ...], list[tuple[int, int]]]] = []
for t0 in valid_t[0]:
    for t_assign in subtree_extensions(0, t0):
        if any(v is None for v in t_assign):
            continue
        assign = [ab_from_t(n, t) for n, t in zip(N_EXPONENTS, t_assign)]
        solutions.append((tuple(int(x) for x in t_assign), assign))


seen = set()
unique_solutions = []
for t_assign, assign in solutions:
    if t_assign in seen:
        continue
    seen.add(t_assign)
    unique_solutions.append((t_assign, assign))


print("=" * 78)
print("VERIFY UNIQUENESS VIA SIMPLE EDGE WILSON LINES")
print("=" * 78)

phys_t = physical_t_values()
print("\nPhysical t-values:")
for name, n, t in zip(NAMES, N_EXPONENTS, phys_t):
    print(f"  {name:>3s}: n={n:>2d}, t={t:>3d}")

count_c1 = 1
for idx in range(1, 9):
    count_c1 *= len(valid_t[idx])
print(f"\nC1 count in |t|<=15: {count_c1:,}")
print(f"Simple-edge count in |t|<=15: {len(unique_solutions)}")

physical_found = any(assign == PHYSICAL_AB for _, assign in unique_solutions)
print(f"Physical assignment selected: {physical_found}")

for k, (t_assign, assign) in enumerate(unique_solutions, start=1):
    print(f"\nSolution #{k}")
    print(f"  t = {t_assign}")
    print(f"  (a,b) = {assign}")
    edge_desc = []
    for i, j in EDGES:
        da = assign[j][0] - assign[i][0]
        db = assign[j][1] - assign[i][1]
        ok, tag = simple_edge_class(da, db)
        assert ok
        coeff, power = tag
        edge_desc.append(f"({NAMES[i]},{NAMES[j]}): {coeff}*phi^{power}")
    print("  edges:")
    for item in edge_desc:
        print(f"    {item}")

print("\n" + "=" * 78)
if len(unique_solutions) == 1 and physical_found:
    print("RESULT: PASS")
    print("Within the bounded domain |t|<=15, C1 plus the simple-edge condition")
    print("dz in {0, +/-phi^r, +/-2phi^r} already selects the physical assignment uniquely.")
else:
    print("RESULT: WARNING")
    print("The simple-edge condition did not uniquely recover the physical assignment.")
print("=" * 78)
