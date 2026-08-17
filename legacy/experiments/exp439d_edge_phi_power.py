"""
exp439d_edge_phi_power.py
==========================
Is the "all edge ratios are phi-powers" property AUTOMATIC from (C1)+(C2),
or is it an INDEPENDENT constraint?

Test: among all Z[phi] assignments satisfying (C1) alone (irreducibility),
how many ALSO have the phi-power edge ratio property?

If very few => it's an independent constraint (new theorem!)
If all/most => it's a consequence of (C1)+(C2) (not new)
"""
import numpy as np
from itertools import product as cartprod

phi = (1 + np.sqrt(5)) / 2
phi_p = -1 / phi
a1 = 5

# McKay graph
names = ['e', 'u', 'd', 's', 'mu', 'c', 'tau', 't', 'b']
n_bare = [0, 3, 5, 11, 11, 16, 17, 19, 26]
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (5,8)]

# For each n, find ALL (a,b) with 5a+6b=n and |N(z)| in {0,1,p prime split/ramified}
# i.e., satisfying (C1)
def norm_zphi(a, b):
    return a**2 + a*b - b**2

def is_irreducible_norm(N):
    """Check if |N| is 0, 1, or a prime that splits/ramifies in Z[phi]."""
    absN = abs(N)
    if absN <= 1:
        return True
    # Check if absN is prime
    if absN < 2:
        return False
    for p in range(2, int(absN**0.5) + 1):
        if absN % p == 0:
            return False
    # It's prime. Check if it splits or ramifies in Z[phi]:
    # p ramifies iff p = 5 (= a1)
    # p splits iff p ≡ ±1 (mod 5)
    # p is inert iff p ≡ ±2 (mod 5)
    if absN == 5:
        return True  # ramified
    if absN % 5 == 1 or absN % 5 == 4:
        return True  # splits (±1 mod 5)
    return False  # inert => NOT irreducible in Z[phi]

# For each bare exponent n, find all valid (a,b) in range |t| <= 15
# where (a,b) = (-n-6t, n+5t)
def find_valid_ab(n, t_range=15):
    valid = []
    for t in range(-t_range, t_range + 1):
        a = -n - 6*t
        b = n + 5*t
        assert 5*a + 6*b == n
        N = norm_zphi(a, b)
        if is_irreducible_norm(N):
            valid.append((a, b, N))
    return valid

print("=" * 80)
print("PART 1: Valid (a,b) assignments per fermion (satisfying C1)")
print("=" * 80)

valid_per_node = []
for i in range(9):
    v = find_valid_ab(n_bare[i])
    valid_per_node.append(v)
    print(f"  {names[i]:3s} (n={n_bare[i]:2d}): {len(v)} valid assignments")
    for a, b, N in v[:5]:
        print(f"      (a={a:+3d}, b={b:+3d}), N={N:+5d}, |N|={abs(N)}")
    if len(v) > 5:
        print(f"      ... and {len(v)-5} more")

total = 1
for v in valid_per_node:
    total *= len(v)
print(f"\nTotal combinations satisfying C1: {total:,}")

# ============================================================
print("\n" + "=" * 80)
print("PART 2: Testing phi-power property on edges")
print("=" * 80)

def edge_is_phi_power(a_i, b_i, a_j, b_j, tol=0.01):
    """Check if dz/dz' is a power of phi (or degenerate)."""
    da, db = a_j - a_i, b_j - b_i
    dz = da + db * phi
    dzg = da + db * phi_p

    if abs(dz) < 1e-12 and abs(dzg) < 1e-12:
        return True, 0, "degenerate"  # both zero
    if abs(dzg) < 1e-12:
        return True, float('inf'), "dz'=0"
    if abs(dz) < 1e-12:
        return True, 0, "dz=0"

    ratio = dz / dzg
    log_phi = np.log(abs(ratio)) / np.log(phi)
    n = round(log_phi)

    if abs(log_phi - n) < tol:
        return True, n, f"phi^{n}"
    return False, log_phi, f"~phi^{log_phi:.3f}"

# Test the ACTUAL assignment
print("\nActual assignment (the unique solution of C1+C2):")
actual_ab = [(0,0), (3,-2), (1,0), (1,1), (1,1), (2,1), (1,2), (4,1), (-1,4)]
all_phi = True
for i, j in edges:
    ai, bi = actual_ab[i]
    aj, bj = actual_ab[j]
    is_pp, n, label = edge_is_phi_power(ai, bi, aj, bj)
    all_phi = all_phi and is_pp
    print(f"  ({names[i]:>3s},{names[j]:>3s}): {label:>10s}  {is_pp}")
print(f"  ALL phi-power: {all_phi}")

# ============================================================
print("\n" + "=" * 80)
print("PART 3: How many C1-valid assignments have ALL edges phi-power?")
print("=" * 80)

# This is a large search space. Let's use the tree structure to prune.
# McKay graph is a tree: e-u-d-s-mu-c with branches c-tau-t and c-b
#
# Strategy: enumerate assignments along the chain, checking edge
# phi-power at each step (pruning early).

count_c1 = 0
count_phi = 0
count_c2 = 0
count_both = 0

# The edge sum for C2: sum over edges of N(z_i - z_j) = -6
TARGET_EDGE_SUM = -6  # = -b_1

# We need to enumerate all 9-tuples from valid_per_node
# Total is ~500M, too many to brute force.
# Use tree decomposition: process along the tree, prune early.

# Tree structure:
# Chain: 0 -> 1 -> 2 -> 3 -> 4 -> 5
# Branch from 5: 5 -> 6 -> 7
# Branch from 5: 5 -> 8

# Process: fix nodes 0-5 along chain, then branch to 6,7 and 8.

print("\nSearching (pruning on chain, then branches)...")
print("This tests ALL ~500M combinations but prunes aggressively.")

phi_power_assignments = []

for a0, b0, N0 in valid_per_node[0]:  # e
    for a1_v, b1_v, N1 in valid_per_node[1]:  # u
        # Check edge (0,1)
        ok01, _, _ = edge_is_phi_power(a0, b0, a1_v, b1_v)

        for a2, b2, N2 in valid_per_node[2]:  # d
            ok12, _, _ = edge_is_phi_power(a1_v, b1_v, a2, b2)

            for a3, b3, N3 in valid_per_node[3]:  # s
                ok23, _, _ = edge_is_phi_power(a2, b2, a3, b3)

                for a4, b4, N4 in valid_per_node[4]:  # mu
                    ok34, _, _ = edge_is_phi_power(a3, b3, a4, b4)

                    for a5, b5, N5 in valid_per_node[5]:  # c
                        ok45, _, _ = edge_is_phi_power(a4, b4, a5, b5)

                        # Check chain phi-power
                        chain_phi = ok01 and ok12 and ok23 and ok34 and ok45

                        for a6, b6, N6 in valid_per_node[6]:  # tau
                            ok56, _, _ = edge_is_phi_power(a5, b5, a6, b6)

                            for a7, b7, N7 in valid_per_node[7]:  # t
                                ok67, _, _ = edge_is_phi_power(a6, b6, a7, b7)

                                for a8, b8, N8 in valid_per_node[8]:  # b
                                    ok58, _, _ = edge_is_phi_power(a5, b5, a8, b8)

                                    count_c1 += 1

                                    all_pp = (chain_phi and ok56 and ok67 and ok58)

                                    # Check C2: edge sum = -6
                                    assignment = [(a0,b0),(a1_v,b1_v),(a2,b2),(a3,b3),
                                                  (a4,b4),(a5,b5),(a6,b6),(a7,b7),(a8,b8)]

                                    edge_sum = 0
                                    for ii, jj in edges:
                                        ai, bi = assignment[ii]
                                        aj, bj = assignment[jj]
                                        da, db = aj-ai, bj-bi
                                        edge_sum += norm_zphi(da, db)

                                    has_c2 = (edge_sum == TARGET_EDGE_SUM)

                                    if all_pp:
                                        count_phi += 1
                                    if has_c2:
                                        count_c2 += 1
                                    if all_pp and has_c2:
                                        count_both += 1
                                        phi_power_assignments.append(assignment)

print(f"\nResults:")
print(f"  Total C1-valid assignments:           {count_c1:>12,}")
print(f"  With ALL edges phi-power:             {count_phi:>12,} "
      f"({100*count_phi/count_c1:.4f}%)")
print(f"  With C2 (edge sum = -6):              {count_c2:>12,} "
      f"({100*count_c2/count_c1:.4f}%)")
print(f"  With BOTH phi-power AND C2:           {count_both:>12,}")

if count_phi > 0:
    print(f"\n  P(C2 | phi-power) = {count_both}/{count_phi} = "
          f"{count_both/count_phi:.4f}")
if count_c2 > 0:
    print(f"  P(phi-power | C2) = {count_both}/{count_c2} = "
          f"{count_both/count_c2:.4f}")

# ============================================================
print("\n" + "=" * 80)
print("PART 4: The phi-power assignments")
print("=" * 80)

if phi_power_assignments:
    for idx, assign in enumerate(phi_power_assignments[:20]):
        print(f"\n  Assignment {idx+1}:")
        for i in range(9):
            a, b = assign[i]
            Nz = norm_zphi(a, b)
            z_val = a + b*phi
            print(f"    {names[i]:3s}: (a={a:+3d}, b={b:+3d}), N={Nz:+5d}, z={z_val:+8.4f}")

        # Edge analysis
        edge_sum = 0
        print("    Edges:")
        for ii, jj in edges:
            ai, bi = assign[ii]
            aj, bj = assign[jj]
            da, db = aj-ai, bj-bi
            N_edge = norm_zphi(da, db)
            edge_sum += N_edge
            _, n, label = edge_is_phi_power(ai, bi, aj, bj)
            print(f"      ({names[ii]:>3s},{names[jj]:>3s}): dz/dz'={label:>8s}, N(dz)={N_edge:+4d}")
        print(f"    Edge sum = {edge_sum} (target: {TARGET_EDGE_SUM})")
else:
    print("  No assignments found with both properties!")

# ============================================================
print("\n" + "=" * 80)
print("PART 5: Phi-power WITHOUT C2 - what do they look like?")
print("=" * 80)

# Show a few phi-power assignments that DON'T satisfy C2
count_shown = 0
for a0, b0, N0 in valid_per_node[0]:
    if count_shown >= 5: break
    for a1_v, b1_v, N1 in valid_per_node[1]:
        if count_shown >= 5: break
        ok01, _, _ = edge_is_phi_power(a0, b0, a1_v, b1_v)
        if not ok01: continue
        for a2, b2, N2 in valid_per_node[2]:
            if count_shown >= 5: break
            ok12, _, _ = edge_is_phi_power(a1_v, b1_v, a2, b2)
            if not ok12: continue
            for a3, b3, N3 in valid_per_node[3]:
                if count_shown >= 5: break
                ok23, _, _ = edge_is_phi_power(a2, b2, a3, b3)
                if not ok23: continue
                for a4, b4, N4 in valid_per_node[4]:
                    if count_shown >= 5: break
                    ok34, _, _ = edge_is_phi_power(a3, b3, a4, b4)
                    if not ok34: continue
                    for a5, b5, N5 in valid_per_node[5]:
                        if count_shown >= 5: break
                        ok45, _, _ = edge_is_phi_power(a4, b4, a5, b5)
                        if not ok45: continue
                        for a6, b6, N6 in valid_per_node[6]:
                            if count_shown >= 5: break
                            ok56, _, _ = edge_is_phi_power(a5, b5, a6, b6)
                            if not ok56: continue
                            for a7, b7, N7 in valid_per_node[7]:
                                if count_shown >= 5: break
                                ok67, _, _ = edge_is_phi_power(a6, b6, a7, b7)
                                if not ok67: continue
                                for a8, b8, N8 in valid_per_node[8]:
                                    if count_shown >= 5: break
                                    ok58, _, _ = edge_is_phi_power(a5, b5, a8, b8)
                                    if not ok58: continue

                                    assign = [(a0,b0),(a1_v,b1_v),(a2,b2),(a3,b3),
                                              (a4,b4),(a5,b5),(a6,b6),(a7,b7),(a8,b8)]
                                    edge_sum = sum(norm_zphi(assign[jj][0]-assign[ii][0],
                                                             assign[jj][1]-assign[ii][1])
                                                   for ii, jj in edges)

                                    if count_shown < 5:
                                        c2_ok = "YES" if edge_sum == -6 else "NO"
                                        norms = [norm_zphi(a,b) for a,b in assign]
                                        print(f"  [{count_shown+1}] norms={norms}, "
                                              f"edge_sum={edge_sum:+d}, C2={c2_ok}")
                                        count_shown += 1

# ============================================================
print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"""
Out of {count_c1:,} assignments satisfying C1 (irreducibility):
  - {count_phi:,} have ALL edges with phi-power Galois ratio ({100*count_phi/count_c1:.2f}%)
  - {count_c2:,} satisfy C2 (edge flatness sum = -6) ({100*count_c2/count_c1:.2f}%)
  - {count_both} satisfy BOTH

If phi-power is RARE (< 1%), it's a strong independent constraint.
If phi-power is COMMON (> 10%), it's likely a consequence.

The key question: does phi-power + C1 IMPLY C2?
  P(C2 | phi-power) = {count_both}/{count_phi if count_phi > 0 else '?'}
""")
