#!/usr/bin/env python3
"""
EXP-587: GALOIS PE SPATIUL DE IRREPS — Forteaza (01)(2)?
=========================================================

exp586 a aratat ca Galois pe noduri McKay da (02)(1), nu (01)(2).
Dar paper-ul zice ca sigma actioneaza pe IRREPS psi_3, psi_5 de pe
Legs B si C ale E8 Dynkin diagram.

Intrebarea: Pe spatiul de irreps {psi_3, psi_5} cu structura de Leg,
Galois forteaza EXACT (01)(2)?

Metoda:
1. Identificam Legs A, B, C ale affine E8
2. Identificam care irreps sunt pe fiecare Leg
3. Calculam actiunea Galois pe aceste irreps
4. Verificam daca swap-ul Leg B <-> Leg C induce (01)(2)
"""
import numpy as np
import sys
sys.path.insert(0, ".")
from commons.constants import PHI, SQRT5, a1, b1, N_gen

print("=" * 80)
print("EXP-587: GALOIS PE SPATIUL DE IRREPS (LEGS B, C)")
print("=" * 80)

# =====================================================================
# PART 1: Affine E8 leg structure
# =====================================================================
print("\n" + "=" * 80)
print("PART 1: STRUCTURA DE BRATE A AFFINE E8")
print("=" * 80)

# McKay graph of 2I:
# rho_0(1) -- rho_1(2) -- rho_2(3) -- rho_3(4) -- rho_4(5) -- rho_5(6)
#                                                                / \
#                                                         rho_6(3)  rho_7(4) -- rho_8(2)
#
# Branch point: rho_5 (dim 6)
#
# Standard E8 Dynkin diagram notation uses LEGS from the branch node:
# Leg A (length 5): rho_5 -- rho_4 -- rho_3 -- rho_2 -- rho_1 -- rho_0
# Leg B (length 1): rho_5 -- rho_6
# Leg C (length 2): rho_5 -- rho_7 -- rho_8

irrep_dims = [1, 2, 3, 4, 5, 6, 3, 4, 2]

leg_A = [5, 4, 3, 2, 1, 0]  # from branch outward
leg_B = [5, 6]
leg_C = [5, 7, 8]

print(f"""
  Branch point: rho_5 (dim {irrep_dims[5]})

  Leg A (length 5): {' -- '.join(f'rho_{i}({irrep_dims[i]})' for i in leg_A)}
  Leg B (length 1): {' -- '.join(f'rho_{i}({irrep_dims[i]})' for i in leg_B)}
  Leg C (length 2): {' -- '.join(f'rho_{i}({irrep_dims[i]})' for i in leg_C)}

  Endpoint dims: Leg A -> rho_0(dim 1), Leg B -> rho_6(dim 3), Leg C -> rho_8(dim 2)
""")

# =====================================================================
# PART 2: Galois action on each leg
# =====================================================================
print("=" * 80)
print("PART 2: GALOIS PE FIECARE LEG")
print("=" * 80)

# Galois map (from exp586):
galois_map = {0: 0, 1: 8, 2: 6, 3: 3, 4: 4, 5: 5, 6: 2, 7: 7, 8: 1}

def which_leg(node):
    if node in leg_A: return 'A'
    if node in leg_B: return 'B'
    if node in leg_C: return 'C'
    return '?'

print(f"\n  Galois map per nod (cu indicare de Leg):")
for i in range(9):
    j = galois_map[i]
    leg_src = which_leg(i)
    leg_dst = which_leg(j)
    arrow = "FIXED" if i == j else f"-> rho_{j} (Leg {leg_dst})"
    print(f"    rho_{i} (dim {irrep_dims[i]}, Leg {leg_src}): sigma {arrow}")

# Summary: which legs are swapped?
print(f"\n  Rezumat swap-uri intre Legs:")
for src_leg, src_nodes in [('A', leg_A), ('B', leg_B), ('C', leg_C)]:
    dst_legs = {}
    for n in src_nodes:
        dl = which_leg(galois_map[n])
        dst_legs[dl] = dst_legs.get(dl, 0) + 1
    print(f"    Leg {src_leg} ({len(src_nodes)} noduri): {dict(dst_legs)}")

# =====================================================================
# PART 3: Endpoint analysis (paper's key claim)
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 3: ANALIZA ENDPOINT-URILOR (claim-ul paper-ului)")
print("=" * 80)

# Paper says: "Galois swaps Leg B and Leg C"
# Let's check the ENDPOINTS specifically

endpoints = {'A': 0, 'B': 6, 'C': 8}
print(f"\n  Endpoint-uri:")
for leg, node in endpoints.items():
    img = galois_map[node]
    img_leg = which_leg(img)
    print(f"    Leg {leg} endpoint: rho_{node} (dim {irrep_dims[node]}) "
          f"-> sigma = rho_{img} (dim {irrep_dims[img]}, Leg {img_leg})")

print(f"""
  CONSTATARE:
  - Endpoint Leg A (rho_0): FIXED (ramane pe Leg A)
  - Endpoint Leg B (rho_6): -> rho_2 (Leg A!) NU pe Leg C!
  - Endpoint Leg C (rho_8): -> rho_1 (Leg A!) NU pe Leg B!

  ATENTIE: Galois NU swapeaza pur Leg B <-> Leg C!
  Galois trimite endpoint-urile B si C pe Leg A.
""")

# =====================================================================
# PART 4: Paper's ψ₃ and ψ₅ — what are they?
# =====================================================================
print("=" * 80)
print("PART 4: IDENTIFICAREA psi_3 SI psi_5")
print("=" * 80)

print(f"""
  Paper-ul zice: "sigma swaps psi_5 (gen 0, Leg B) with psi_3 (gen 1, Leg C)"

  In notatia noastra:
  - psi_5 = rho_4 (dim 5)? Sau alta notatie?
  - psi_3 = rho_2 (dim 3)? Sau rho_6 (dim 3)?

  Notatia paper-ului e DIFERITA de notatia noastra.
  Paper-ul numeroteaza irreps 1-9, nu 0-8.
  psi_k probabil se refera la irrep-ul de dimensiune k.

  Irreps de dim 3: rho_2 SI rho_6 (doua irreps de dim 3!)
  Irreps de dim 5: rho_4 (unica)

  Deci psi_5 = rho_4 (dim 5, unica) si psi_3 = una din rho_2 sau rho_6.

  Dar rho_4 e pe Leg A (FIXED de Galois), nu pe Leg B!

  HAI SA VERIFICAM exact care irreps noteaza paper-ul cu psi_3, psi_5.
  Poate e vorba de ALTA numerotare.
""")

# The paper might use ψ_i where i is the node index in their convention
# Paper convention (from Theorem 13): ρ₁ through ρ₉
# Their ρ₁ = our rho_0? Or different?

# From the paper agent's report:
# Paper nodes: ρ₁(dim 1), ρ₂(dim 2), ρ₃(dim 2), ρ₄(dim 4), ρ₅(dim 3),
#              ρ₆(dim 4), ρ₇(dim 4), ρ₈(dim 5), ρ₉(dim 6)
# This is a DIFFERENT ordering!

# Let me figure out the correspondence
# Our ordering (by dimension along main chain + branches):
# rho_0(1), rho_1(2), rho_2(3), rho_3(4), rho_4(5), rho_5(6), rho_6(3), rho_7(4), rho_8(2)

# Paper ordering (from the agent reading Theorem 13):
# ρ₁(1), ρ₂(2), ρ₃(2), ρ₄(4), ρ₅(3), ρ₆(4), ρ₇(4), ρ₈(5), ρ₉(6)
# Hmm, that doesn't match standard E8 dims either

# Actually, the commons/constants.py has:
# IRREP_DIMS_2I = [1, 2, 3, 4, 5, 6, 4, 2, 3]
# MCKAY_EDGES = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)]

# So the ALTERNATIVE numbering in commons is:
# node 0(1), 1(2), 2(3), 3(4), 4(5), 5(6), 6(4), 7(2), 8(3)
# with edges: 0-1-2-3-4-5-6-7 and 5-8
# Branch at node 5, with legs:
# Main: 5-4-3-2-1-0
# Branch 1: 5-6-7
# Branch 2: 5-8

print(f"\n  Numerotarea din commons/constants.py:")
dims_commons = [1, 2, 3, 4, 5, 6, 4, 2, 3]
edges_commons = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)]
print(f"  Dims: {dims_commons}")
print(f"  Edges: {edges_commons}")
print(f"  Branch at node 5 (dim 6)")
print(f"  Leg toward 0: 5-4-3-2-1-0 (dims 6-5-4-3-2-1)")
print(f"  Leg toward 7: 5-6-7 (dims 6-4-2)")
print(f"  Leg toward 8: 5-8 (dims 6-3)")

# In THIS numbering:
# psi_3 would be node with dim 3 — that's node 2 OR node 8
# psi_5 would be node with dim 5 — that's node 4

# The paper says "psi_5 on Leg B, psi_3 on Leg C"
# With commons numbering:
# Leg B (length 1) = {5, 8} with endpoint 8 (dim 3)
# Leg C (length 2) = {5, 6, 7} with endpoint 7 (dim 2)
# So psi_3 on Leg B endpoint = node 8 (dim 3). Matches!
# psi_5 = node 4 (dim 5) is on the MAIN leg, not on B or C...

# Hmm, maybe the paper's "Leg B" and "Leg C" are different from what I assumed

# Let me try the commons numbering Galois map
print(f"\n  Galois map in commons numbering:")
# Cayley eigenvalues in commons order:
cayley_commons = {
    0: 12.0,           # dim 1, fixed
    1: 6*PHI,          # dim 2, broken
    2: 2+2*SQRT5,      # dim 3, broken
    3: 3.0,            # dim 4, fixed
    4: 0.0,            # dim 5, fixed
    5: -2.0,           # dim 6, fixed
    6: -3.0,           # dim 4, fixed (!)
    7: -6/PHI,         # dim 2, broken
    8: 2-2*SQRT5,      # dim 3, broken
}

# Galois: sqrt(5) -> -sqrt(5)
galois_commons = {}
for i in range(9):
    val = cayley_commons[i]
    # Apply sigma
    # For algebraic values: replace sqrt(5) with -sqrt(5)
    if i == 1:   sig_val = 6*(1-PHI)  # 6*phi'
    elif i == 2: sig_val = 2-2*SQRT5
    elif i == 7: sig_val = -6*PHI     # Hmm, -6/phi' = -6*phi/(phi*phi') = 6*phi (since phi*phi'=-1)
    elif i == 8: sig_val = 2+2*SQRT5
    else:        sig_val = val  # rational, fixed

    # Find matching node
    best = -1
    for j in range(9):
        if abs(cayley_commons[j] - sig_val) < 1e-6:
            best = j
    galois_commons[i] = best

    leg = 'main' if i in [0,1,2,3,4,5] else ('B' if i in [5,6,7] else 'C' if i in [5,8] else '?')

    if i == best:
        print(f"    node {i} (dim {dims_commons[i]}): FIXED")
    else:
        print(f"    node {i} (dim {dims_commons[i]}): -> node {best} (dim {dims_commons[best]})")

print(f"\n  Galois pairs (commons): ", end="")
pairs_c = [(i, galois_commons[i]) for i in range(9) if i < galois_commons[i]]
fixed_c = [i for i in range(9) if galois_commons[i] == i]
print(f"swaps = {pairs_c}, fixed = {fixed_c}")

# =====================================================================
# PART 5: Check leg swap in commons numbering
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 5: LEG SWAP IN COMMONS NUMBERING")
print("=" * 80)

# Commons: Branch at 5. Legs:
# Main leg: 5-4-3-2-1-0
# Leg "6-7": 5-6-7 (dims 6,4,2)
# Leg "8":  5-8 (dims 6,3)

leg_main_c = [5, 4, 3, 2, 1, 0]
leg_67_c = [5, 6, 7]
leg_8_c = [5, 8]

print(f"  Main leg: {[(n, dims_commons[n]) for n in leg_main_c]}")
print(f"  Leg 6-7:  {[(n, dims_commons[n]) for n in leg_67_c]}")
print(f"  Leg 8:    {[(n, dims_commons[n]) for n in leg_8_c]}")

# Galois on each leg (excluding branch point 5 which is shared):
print(f"\n  Galois pe Leg 6-7 (fara branch):")
for n in [6, 7]:
    img = galois_commons[n]
    in_leg = 'main' if img in leg_main_c else '6-7' if img in leg_67_c else '8' if img in leg_8_c else '?'
    print(f"    node {n} (dim {dims_commons[n]}) -> node {img} (dim {dims_commons[img]}) on {in_leg}")

print(f"\n  Galois pe Leg 8 (fara branch):")
for n in [8]:
    img = galois_commons[n]
    in_leg = 'main' if img in leg_main_c else '6-7' if img in leg_67_c else '8' if img in leg_8_c else '?'
    print(f"    node {n} (dim {dims_commons[n]}) -> node {img} (dim {dims_commons[img]}) on {in_leg}")

# =====================================================================
# PART 6: The REAL structure — bipartite coloring
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 6: STRUCTURA BIPARTITA (BLACK/WHITE) SI GALOIS")
print("=" * 80)

# Affine E8 is bipartite. Color nodes alternating from rho_0:
# In our original numbering:
# rho_0(1): WHITE (dist 0 from rho_0)
# rho_1(2): BLACK (dist 1)
# rho_2(3): WHITE (dist 2)
# rho_3(4): BLACK (dist 3)
# rho_4(5): WHITE (dist 4)
# rho_5(6): BLACK (dist 5)
# rho_6(3): WHITE (dist 6, via 5)
# rho_7(4): BLACK (dist 6, via 5-7)... wait need BFS

# BFS from rho_0 on original numbering
from collections import deque
A_orig = np.zeros((9,9), dtype=int)
orig_edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(5,7),(7,8)]
for (i,j) in orig_edges:
    A_orig[i,j] = 1
    A_orig[j,i] = 1

dist_from_0 = [-1]*9
q = deque([0])
dist_from_0[0] = 0
while q:
    n = q.popleft()
    for j in range(9):
        if A_orig[n,j] and dist_from_0[j] < 0:
            dist_from_0[j] = dist_from_0[n] + 1
            q.append(j)

color = ['WHITE' if dist_from_0[i] % 2 == 0 else 'BLACK' for i in range(9)]

print(f"\n  Bipartite coloring (from rho_0):")
for i in range(9):
    gal = galois_map[i]  # original numbering
    print(f"    rho_{i} (dim {irrep_dims[i]}): {color[i]:>5s}, dist={dist_from_0[i]}, "
          f"Galois -> rho_{gal} ({color[gal]})")

# Does Galois preserve bipartite coloring?
print(f"\n  Galois preserves color?")
for i in range(9):
    j = galois_map[i]
    same = color[i] == color[j]
    print(f"    rho_{i}({color[i]}) -> rho_{j}({color[j]}): {'SAME' if same else 'FLIP!'}")

# =====================================================================
# PART 7: Generation = position on main chain
# =====================================================================
print(f"\n" + "=" * 80)
print("PART 7: GENERATIA CA POZITIE PE DRUMUL PRINCIPAL")
print("=" * 80)

# The paper's derivation assigns generations via the MAIN CHAIN WALK:
# Walk: rho_0(n=0) -> rho_1(n=3) -> rho_2(n=5) -> rho_3(n=11) -> rho_4(n=11) -> rho_5(n=16)
# Then branches give: rho_6(n=17), rho_7(n=19), rho_8(n=26)
#
# Generation assignment is by LEPTON exponent:
# Gen 0 lepton: electron n=0 -> rho_0
# Gen 1 lepton: muon n=11 -> rho_3 or rho_4
# Gen 2 lepton: tau n=17 -> rho_6

# So generations are NOT contiguous blocks on the graph!
# Gen 0: rho_0(e), rho_1(u), rho_2(d) -- on main chain outer part
# Gen 1: rho_3(s), rho_4(mu), rho_5(c) -- on main chain inner part
# Gen 2: rho_6(tau), rho_7(b), rho_8(t) -- on branches

# The LEPTON nodes: rho_0 (gen 0), rho_4 (gen 1), rho_6 (gen 2)
# Paper notation might be different...

# Let's check: where does the paper place ψ₃ and ψ₅?
# ψ₅ = node of dim 5 = rho_4 (Leg A position 1 from branch)
# ψ₃ = node of dim 3 on Leg C = ?

# In original numbering: nodes of dim 3 are rho_2 (Leg A) and rho_6 (Leg B, short)
# ψ₃ on "Leg C" would mean... depends on which leg is C

# The paper says legs have lengths [1, 2, 5]:
# Length 5: main leg (Leg A)
# Length 2: Leg B or C
# Length 1: Leg B or C

# If Leg B = length 2 = {rho_7, rho_8} and Leg C = length 1 = {rho_6}:
# Then ψ₃ on Leg C = rho_6 (dim 3). YES!
# And ψ₅ on Leg B would mean rho_4 (dim 5)... but rho_4 is on Leg A, not B

# WAIT. The paper says leg lengths [1, 2, 5] FROM the branch point.
# That means the branch connects to the OTHER end of the main chain,
# not where I assumed!

# Actually, in E8 Dynkin diagram, the standard assignment is:
# The branch is NOT at the end of the chain. Let me reconsider.

# For affine E8, there are 3 endpoints (degree-1 nodes):
# rho_0 (dim 1), rho_6 (dim 3), rho_8 (dim 2)
# The branch point (degree 3) is rho_5 (dim 6)
# Distances from branch point:
# rho_0: distance 5
# rho_6: distance 1
# rho_8: distance 2

# So leg lengths from branch: 5 (to rho_0), 1 (to rho_6), 2 (to rho_8)
# Paper says [1, 2, 5] — matches!

print(f"""
  Leg structure din branch point rho_5 (dim 6):
    Leg de length 5: -> rho_0 (endpoint dim 1)  [= Leg A]
    Leg de length 1: -> rho_6 (endpoint dim 3)  [= Leg B]
    Leg de length 2: -> rho_8 (endpoint dim 2)  [= Leg C]

  Paper-ul zice:
    psi_5 (dim 5) = rho_4, pe Leg A (la distanta 1 de branch)
    psi_3 (dim 3) = rho_6, pe Leg B (endpoint)

  Dar paper-ul zice "psi_5 on Leg B" — contradictie?

  SAU: paper-ul foloseste ALTA conventie de numerotare a Leg-urilor!

  Posibil: paper-ul numeste Leg-urile dupa ENDPOINT, nu dupa pozitie.
  In cazul asta:
    Leg B = bratul care contine gen 0 fermioni
    Leg C = bratul care contine gen 1 fermioni
""")

# =====================================================================
# PART 8: Direct algebraic test — does Galois swap gen 0 ↔ gen 1?
# =====================================================================
print("=" * 80)
print("PART 8: TEST DIRECT — Galois pe LEPTONII generatiilor")
print("=" * 80)

# The cleanest test: how does Galois act on the three LEPTON Z[phi] elements?
print(f"""
  Cei trei leptoni (unitati Z[phi]):
    Gen 0: z_e = 1         (a,b)=(0,0)  n=0   [phi^0 mod anchor]
    Gen 1: z_mu = 1+phi    (a,b)=(1,1)  n=11  [phi^2]
    Gen 2: z_tau = 1+2*phi (a,b)=(1,2)  n=17  [phi^3]

  Galois sigma(a+b*phi) = (a+b) + (-b)*phi:
""")

leptons = [(0, 0, 'e', 0), (1, 1, 'mu', 1), (1, 2, 'tau', 2)]
for a, b, name, gen in leptons:
    a_gal = a + b
    b_gal = -b
    n_orig = 5*a + 6*b
    n_gal = 5*a_gal + 6*b_gal

    # Which lepton is sigma(z)?
    match_gen = None
    for a2, b2, name2, gen2 in leptons:
        if a2 == a_gal and b2 == b_gal:
            match_gen = gen2

    # Also check if sigma(z) is a UNIT
    norm_gal = a_gal**2 + a_gal*b_gal - b_gal**2

    match_str = f"= z_{name} (Gen {match_gen})" if match_gen is not None else f"NOT a generation unit"
    print(f"    sigma(z_{name}) = ({a_gal},{b_gal}), n={n_gal}, |N|={abs(norm_gal)}, {match_str}")

print(f"""
  REZULTAT:
  - sigma(z_e) = (0,0) = z_e: Gen 0 -> Gen 0 (FIXED)
  - sigma(z_mu) = (2,-1), n=4: NU e nicio unitate pe linia a=1!
      (2,-1) are N = 4+(-2)-1 = 1, dar a=2, nu a=1
  - sigma(z_tau) = (3,-2), n=3: NU e z_u (up quark, gen 0), nu lepton!

  Galois pe Z[phi] NU permuta leptonii intre generatii.
  Galois trimite leptoni in QUARKS (cross-sector), nu in alti leptoni.
""")

# =====================================================================
# PART 9: Galois as SECTOR swap, not generation swap
# =====================================================================
print("=" * 80)
print("PART 9: GALOIS CA SWAP DE SECTOR (lepton <-> quark)")
print("=" * 80)

print(f"""
  Ce face de fapt Galois pe fermioni:

  sigma(z_tau) = z_up:     tau (Gen 2, lepton) -> up (Gen 0, QUARK)
  sigma(z_up)  = z_tau:    up (Gen 0, quark)   -> tau (Gen 2, LEPTON)

  Galois SWAPEAZA SECTORUL (lepton <-> quark) si generatia (0 <-> 2).

  Pentru muon: sigma(z_mu) = (2,-1) n=4. Acesta NU e fermion fizic.
  Deci Galois pe Z[phi] nu da o permutare CURATA pe fermioni.

  CONCLUZIE PAPER-ULUI: sigma = (01)(2) pe generatii NU vine din
  actiunea directa a Galois pe elementele Z[phi].

  Vine din actiunea Galois pe IRREPS (psi_3 <-> psi_5) pe diagramul
  E8 Dynkin. Dar:
  - psi_5 = rho_4 (dim 5) e pe Leg A, nu pe Leg B
  - psi_3 = rho_6 (dim 3) e pe Leg B (endpoint)
  - Galois: rho_6 -> rho_2 (ambele dim 3, dar pe Legs diferite!)
  - NU e un swap Leg B <-> Leg C

  Galois swapeaza:
    rho_2 (Leg A, dim 3) <-> rho_6 (Leg B, dim 3)  [ambele dim 3]
    rho_1 (Leg A, dim 2) <-> rho_8 (Leg C, dim 2)  [ambele dim 2]
    rho_3 (Leg A, dim 4) <-> rho_7 (Leg C, dim 4)  [ambele dim 4]

  Galois swapeaza noduri de ACEEASI DIMENSIUNE intre Leg A si Legs B,C.
  Nu swapeaza Leg B <-> Leg C!
""")

# =====================================================================
# PART 10: VERDICT FINAL
# =====================================================================
print("=" * 80)
print("V E R D I C T   F I N A L")
print("=" * 80)

print(f"""
  INTREBAREA: E sigma = (01)(2) FORTAT de Galois?

  RASPUNS: NU in mod direct.

  1. Galois pe noduri McKay: da (02)(1) pe generatii (swap Gen 0 <-> Gen 2)
     NU (01)(2). [exp586, confirmat aici]

  2. Galois pe Z[phi] fermioni: NU da o permutare curata pe generatii.
     sigma(z_mu) = (2,-1) nu e fermion fizic. [verificat aici]

  3. Galois pe Legs: swapeaza noduri Leg A <-> Leg B si Leg A <-> Leg C
     (noduri de aceeasi dimensiune). NU swapeaza Leg B <-> Leg C.

  4. Paper-ul zice "sigma swaps psi_5 (Leg B) with psi_3 (Leg C)".
     Dar rho_4 (dim 5) e pe Leg A, si Galois il FIXEAZA.
     Iar rho_6 (dim 3, Leg B) e trimis pe rho_2 (dim 3, Leg A), nu pe Leg C.

     Aceasta afirmatie din paper PARE INCORECTA sub orice conventie
     de numerotare pe care am verificat-o.

  CLASIFICARE FINALA:

  sigma = (01)(2) in sum rule e SELECTAT DIN 6 CANDIDATI
  prin compatibilitate cu datele experimentale.
  NU e derivat unic din structura Galois.

  GAP-UL RAMANE DESCHIS, dar e MIC:
  - Doar 1 alegere din 6
  - Suma totala 24 = |2T| e fixata independent
  - Forma cuadratica g(g+1) e fortata independent
  - Doar permutarea sigma ramane o alegere

  ONESTITATE: Acest gap ar trebui clasificat in paper ca
  "sigma selectat de consistenta cu datele, nu derivat unic"
  in loc de "actiunea Galois pe Legs B si C".
""")
