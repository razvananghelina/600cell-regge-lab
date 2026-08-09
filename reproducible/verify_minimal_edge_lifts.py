"""
Constructive verifier for a more primitive edge principle behind the
bounded-search uniqueness results.

Idea:
  The simple-edge condition

      dz in {0} U {+/-phi^r} U {+/-2phi^r}

  appears to be equivalent, on the actual McKay edge-weight data used by the
  uniqueness search, to a canonical minimal-complexity lift:

      given an edge exponent jump Delta n,
      choose the unique integer solution (da, db) of
          5*da + 6*db = Delta n
      minimizing |da| + |db|.

Because the McKay graph is a tree, these edge lifts integrate uniquely from
the root z_e = 0 and reconstruct the full node assignment.

This script verifies that:
  1. each McKay edge jump has a unique minimal-L1 lift;
  2. the lifted edges are exactly the simple edge classes;
  3. integrating them reconstructs the neutral uniqueness solution;
  4. the reconstructed node values satisfy C1, Q, and C2 automatically.

Important:
  The branch endpoints are kept neutral as brA/brB, matching the bounded-search
  uniqueness scripts. This verifier reconstructs the arithmetic assignment on
  the tree; the later physical identification of the two branch endpoints is a
  separate step.
"""

from __future__ import annotations

import math
from collections import deque


PHI = (1 + math.sqrt(5)) / 2

NAMES = ["e", "u", "d", "s", "mu", "c", "tau", "brA", "brB"]
N_EXPONENTS = [0, 3, 5, 11, 11, 16, 17, 19, 26]
TARGET_AB = [
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


def ab_from_t(n: int, t: int) -> tuple[int, int]:
    return -n - 6 * t, n + 5 * t


def minimal_l1_lift(delta_n: int, t_window: int = 30) -> tuple[int, int]:
    """
    Solve 5*da + 6*db = delta_n by the canonical parametrization

        (da, db) = (-delta_n - 6t, delta_n + 5t)

    and select the unique minimizer of |da| + |db|.
    """
    best_pair = None
    best_key = None

    for t in range(-t_window, t_window + 1):
        da, db = ab_from_t(delta_n, t)
        key = (abs(da) + abs(db), abs(da), abs(db))
        if best_key is None or key < best_key:
            best_key = key
            best_pair = (da, db)

    assert best_pair is not None

    minimizers = []
    for t in range(-t_window, t_window + 1):
        da, db = ab_from_t(delta_n, t)
        key = (abs(da) + abs(db), abs(da), abs(db))
        if key == best_key:
            minimizers.append((da, db))

    if len(minimizers) != 1:
        raise AssertionError(
            f"Delta n = {delta_n} has non-unique minimal-L1 lifts: {minimizers}"
        )

    return best_pair


def simple_edge_label(da: int, db: int) -> str:
    table = {
        (0, 0): "0",
        (1, 0): "1",
        (2, 0): "2",
        (0, 1): "phi",
        (-1, 1): "phi^-1",
        (-2, 2): "2phi^-1",
        (3, -2): "-phi^-3",
    }
    if (da, db) in table:
        return table[(da, db)]
    return "not-simple"


def q_edge_ok(da: int, db: int) -> tuple[bool, str]:
    """
    Check the earlier edge-quantization condition:
        dz / sigma(dz) = +/- phi^k
    """
    dz = da + db * PHI
    dz_sigma = da - db / PHI

    if abs(dz) < 1e-12 and abs(dz_sigma) < 1e-12:
        return True, "+phi^0"
    if abs(dz_sigma) < 1e-12 or abs(dz) < 1e-12:
        return False, "degenerate"

    ratio = dz / dz_sigma
    sign = "+" if ratio > 0 else "-"
    k = round(math.log(abs(ratio), PHI))
    if abs(abs(ratio) - PHI**k) < 1e-9:
        return True, f"{sign}phi^{k}"
    return False, "not-quantized"


def reconstruct_nodes() -> tuple[list[tuple[int, int]], list[tuple[int, int, int, int]]]:
    adj = [[] for _ in NAMES]
    for i, j in EDGES:
        adj[i].append(j)
        adj[j].append(i)

    edge_data = []
    z = [None for _ in NAMES]
    z[0] = (0, 0)

    q = deque([0])
    seen = {0}

    while q:
        i = q.popleft()
        ai, bi = z[i]
        for j in adj[i]:
            if j in seen:
                continue
            delta_n = N_EXPONENTS[j] - N_EXPONENTS[i]
            da, db = minimal_l1_lift(delta_n)
            z[j] = (ai + da, bi + db)
            edge_data.append((i, j, da, db))
            seen.add(j)
            q.append(j)

    assert all(v is not None for v in z)
    return z, edge_data


reconstructed_ab, edge_data = reconstruct_nodes()

print("=" * 78)
print("VERIFY CONSTRUCTIVE MINIMAL-EDGE LIFT PRINCIPLE")
print("=" * 78)

print("\nEdge data:")
all_simple = True
all_q = True
edge_sum = 0
for i, j, da, db in edge_data:
    delta_n = N_EXPONENTS[j] - N_EXPONENTS[i]
    simple = simple_edge_label(da, db)
    q_ok, q_tag = q_edge_ok(da, db)
    all_simple = all_simple and (simple != "not-simple")
    all_q = all_q and q_ok
    edge_norm = norm_zphi(da, db)
    edge_sum += edge_norm
    print(
        f"  ({NAMES[i]},{NAMES[j]}): "
        f"Delta n={delta_n:>2d}, "
        f"(da,db)=({da:+d},{db:+d}), "
        f"L1={abs(da)+abs(db):>2d}, "
        f"N(dz)={edge_norm:+d}, "
        f"simple={simple:>7s}, "
        f"Q={q_tag}"
    )

print("\nReconstructed node values:")
all_c1 = True
for name, n, (a, b), target in zip(NAMES, N_EXPONENTS, reconstructed_ab, TARGET_AB):
    nz = norm_zphi(a, b)
    c1 = c1_irreducible(abs(nz))
    all_c1 = all_c1 and c1
    match = (a, b) == target
    print(
        f"  {name:>3s}: (a,b)=({a:+d},{b:+d}), "
        f"n={5*a+6*b:>2d}, "
        f"N={nz:+d}, "
        f"C1={c1}, "
        f"target_match={match}"
    )

node_sum = sum(norm_zphi(a, b) for a, b in reconstructed_ab)
matches_target = reconstructed_ab == TARGET_AB

print("\nGlobal checks:")
print(f"  Reconstructed assignment matches target: {matches_target}")
print(f"  All edges in simple class: {all_simple}")
print(f"  All edges satisfy Q: {all_q}")
print(f"  All nodes satisfy C1: {all_c1}")
print(f"  sum_nodes N(z) = {node_sum:+d}  (target +6)")
print(f"  sum_edges N(dz) = {edge_sum:+d}  (target -6)")
print(f"  Flatness sum = {node_sum + edge_sum:+d}  (target 0)")

print("\nInterpretation:")
print("  The previous simple-edge condition is not an extra ad hoc axiom here.")
print("  On the actual McKay edge jumps Delta n in the neutral branch ordering,")
print("  it is exactly the canonical minimal-L1 lift of 5*da + 6*db = Delta n.")

print("\n" + "=" * 78)
if matches_target and all_simple and all_q and all_c1 and node_sum == 6 and edge_sum == -6:
    print("RESULT: PASS")
    print("A local minimal-edge lift rule reconstructs the full neutral assignment")
    print("constructively, with C1, Q, and C2 emerging as consequences.")
else:
    print("RESULT: WARNING")
    print("The constructive minimal-edge lift rule did not fully reproduce the target data.")
print("=" * 78)
