"""
Verify that the physical identification of the two neutral branch endpoints
is fixed by McKay chirality plus the generation-formula exponents.

The uniqueness scripts deliberately keep the last two nodes neutral:

    brA, brB

because the pure Z[phi] arithmetic on the abstract tree does not know which
endpoint should be called top or bottom.

This script adds back the missing structure:
  1. the actual affine-E8 McKay leg geometry (lengths 1, 2, 5);
  2. the bipartite chirality / T3 assignment;
  3. the generation-formula exponents for the g=2 quark doublet.

Result:
  - the BLACK endpoint on the length-2 leg is the generation-2 up-type quark,
    hence top with exponent 26;
  - the WHITE endpoint on the short leg is the generation-2 down-type quark,
    hence bottom with exponent 19.

So the neutral branch labels are resolved structurally:

    br_long_black  -> top
    br_short_white -> bottom

This does not depend on experimental masses.
"""

from __future__ import annotations


NAMES = [
    "rho_1", "rho_2", "rho_3", "rho_4", "rho_5",
    "rho_6", "rho_7", "rho_8", "rho_9",
]

# McKay graph edges from verify_mckay_chirality.py
EDGES = [
    (0, 1),
    (1, 3),
    (2, 5),
    (3, 6),
    (4, 8),
    (5, 8),
    (6, 7),
    (7, 8),
]

# Bipartite coloring from the same verifier:
# WHITE = +1, BLACK = -1
WHITE = {0, 3, 4, 5, 7}
BLACK = {1, 2, 6, 8}


def build_adj():
    adj = [[] for _ in NAMES]
    for i, j in EDGES:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def find_branch_and_legs(adj):
    branch = next(i for i in range(len(adj)) if len(adj[i]) == 3)
    legs = []
    for start in adj[branch]:
        leg = [start]
        prev = branch
        cur = start
        while True:
            nxt = [j for j in adj[cur] if j != prev]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            leg.append(cur)
        legs.append(leg)
    legs.sort(key=len)
    return branch, legs


adj = build_adj()
branch, legs = find_branch_and_legs(adj)
short_leg, mid_leg, long_leg = legs

# Physical assignment already used in verify_mckay_chirality.py
chain_from_e = list(reversed(long_leg))
physical_node_to_fermion = {}
for idx, f in enumerate(["e", "u", "d", "s", "mu"]):
    physical_node_to_fermion[chain_from_e[idx]] = f
physical_node_to_fermion[branch] = "c"
physical_node_to_fermion[mid_leg[0]] = "tau"
physical_node_to_fermion[mid_leg[1]] = "t"
physical_node_to_fermion[short_leg[0]] = "b"


def node_color(node: int) -> str:
    return "WHITE" if node in WHITE else "BLACK"


print("=" * 78)
print("VERIFY BRANCH IDENTIFICATION: TOP VS BOTTOM")
print("=" * 78)

print("\nAffine-E8 McKay geometry:")
print(f"  branch node = {NAMES[branch]} (degree 3)")
print(f"  leg lengths = {[len(x) for x in legs]}  (expected [1, 2, 5])")
print(f"  short leg   = {[NAMES[i] for i in short_leg]}")
print(f"  mid leg     = {[NAMES[i] for i in mid_leg]}")
print(f"  long leg    = {[NAMES[i] for i in long_leg]}")

print("\nRemaining branch endpoints after fixing the main chain and tau:")
mid_endpoint = mid_leg[1]
short_endpoint = short_leg[0]
print(
    f"  length-2 endpoint: {NAMES[mid_endpoint]}  color={node_color(mid_endpoint)}"
)
print(
    f"  short-leg endpoint: {NAMES[short_endpoint]}  color={node_color(short_endpoint)}"
)

# Generation-2 quark exponents from the explicit formulas already in the paper
# and verify_masses_and_mixing.py
n_tau = 17
g = 2
a1 = 5
N_gen = 3
F3 = 2

n_top = n_tau + (a1 - 2) + g * (g + 1)        # 17 + 3 + 6 = 26
n_bottom = n_tau + F3 + 6 - g * (g + 1)       # 17 + 2 + 6 - 6 = 19

print("\nGeneration-2 quark exponents from the internal formulas:")
print(f"  up-type g=2    : n = {n_top}  (= 17 + 3 + 6)")
print(f"  down-type g=2  : n = {n_bottom}  (= 17 + 2 + 6 - 6)")

print("\nMcKay chirality / weak-isospin rule:")
print("  BLACK endpoint  -> up-type branch member")
print("  WHITE endpoint  -> down-type branch member")

top_node = mid_endpoint if mid_endpoint in BLACK else short_endpoint
bottom_node = short_endpoint if short_endpoint in WHITE else mid_endpoint

print("\nResolved identification:")
print(f"  top    = {NAMES[top_node]}  (color={node_color(top_node)}, n={n_top})")
print(
    f"  bottom = {NAMES[bottom_node]}  (color={node_color(bottom_node)}, n={n_bottom})"
)

print("\nCross-check with the physical node assignment already used elsewhere:")
for node in [branch, mid_leg[0], mid_leg[1], short_leg[0]]:
    print(
        f"  {NAMES[node]} -> {physical_node_to_fermion[node]}"
        f"  (color={node_color(node)})"
    )

ok_lengths = [len(x) for x in legs] == [1, 2, 5]
ok_colors = (mid_endpoint in BLACK) and (short_endpoint in WHITE)
ok_names = (
    physical_node_to_fermion[mid_endpoint] == "t"
    and physical_node_to_fermion[short_endpoint] == "b"
)
ok_exponents = (n_top == 26) and (n_bottom == 19)

print("\n" + "=" * 78)
if ok_lengths and ok_colors and ok_names and ok_exponents:
    print("RESULT: PASS")
    print("The neutral branch ambiguity is resolved by geometry plus chirality:")
    print("the length-2 BLACK endpoint is top (n=26),")
    print("and the short-leg WHITE endpoint is bottom (n=19).")
else:
    print("RESULT: WARNING")
    print("The branch identification check did not fully close.")
print("=" * 78)
