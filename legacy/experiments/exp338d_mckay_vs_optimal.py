"""
exp338d: Compare McKay assignment (exp323) with VEV optimal assignment
======================================================================
Two different numbering schemes exist:
- exp338 (Cayley): rho_0..8 with dims [1,2,3,4,5,6,4,2,3]
- exp307/323 (McKay): rho_1..9 with dims [1,2,2',3,3',4,4',5,6]

The McKay graph (affine E8):
  rho_3(2')
    |
  rho_6(4)
    |
  rho_1(1)--rho_2(2)--rho_4(3)--rho_7(4')--rho_8(5)--rho_9(6)--rho_5(3')

Mapping McKay -> Cayley:
  mckay rho_1(dim 1)  = cayley rho_0(dim 1)
  mckay rho_2(dim 2)  = cayley rho_1(dim 2)
  mckay rho_3(dim 2') = cayley rho_7(dim 2)  [Galois conjugate]
  mckay rho_4(dim 3)  = cayley rho_2(dim 3)
  mckay rho_5(dim 3') = cayley rho_8(dim 3)  [Galois conjugate]
  mckay rho_6(dim 4)  = cayley rho_3(dim 4)
  mckay rho_7(dim 4') = cayley rho_6(dim 4)  [Galois conjugate]
  mckay rho_8(dim 5)  = cayley rho_4(dim 5)
  mckay rho_9(dim 6)  = cayley rho_5(dim 6)
"""
import numpy as np
from itertools import permutations

PHI = (1 + np.sqrt(5)) / 2
ALPHA_S = 1 / (2 * PHI**3)

# Cayley numbering
dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
char_C1 = [1, PHI, PHI, 1, 0, -1, -1, -1/PHI, -1/PHI]
cayley_eigs = np.array([12 * (1 - char_C1[k] / dims[k]) for k in range(9)])

even_modes = [2, 4, 6, 8]
lam = np.array([cayley_eigs[m] for m in even_modes])
d_m = np.array([dims[m] for m in even_modes])

CG_self = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 1, 0, 1],
    [1, 0, 1, 0, 2, 0, 2, 0, 1],
    [1, 0, 2, 0, 3, 0, 2, 0, 2],
    [1, 0, 1, 0, 1, 0, 1, 0, 1],
    [1, 0, 0, 0, 0, 0, 0, 0, 1],
    [1, 0, 0, 0, 1, 0, 0, 0, 1],
])

fermion_names = ['e', 'mu', 'tau', 'u', 'c', 't', 'd', 's', 'b']
fermion_T3 = np.array([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5])
fermion_Nc = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
delta_exp = np.array([0.000, 0.080, -0.055, -0.004, 0.247, 0.456, -0.402, -0.177, -0.278])
mass_exp = np.array([0, 11, 17, 3, 16, 26, 5, 11, 19])  # n = 5a+6b

print("="*72)
print("exp338d: McKay Assignment vs VEV Optimal Assignment")
print("="*72)

# ================================================================
# McKay numbering -> Cayley numbering
# ================================================================
# McKay: rho_1..9 with dims 1,2,2',3,3',4,4',5,6
# Cayley: rho_0..8 with dims 1,2,3,4,5,6,4,2,3
mckay_to_cayley = {
    1: 0,  # dim 1
    2: 1,  # dim 2
    3: 7,  # dim 2' (Galois conj)
    4: 2,  # dim 3
    5: 8,  # dim 3' (Galois conj)
    6: 3,  # dim 4
    7: 6,  # dim 4' (Galois conj)
    8: 4,  # dim 5
    9: 5,  # dim 6
}

print("\nNumbering correspondence:")
print(f"{'McKay':>8} {'dim':>4} {'Cayley':>8} {'dim':>4} {'Galois?':>8}")
print("-"*40)
for mk, cy in sorted(mckay_to_cayley.items()):
    gal = ""
    if mk in [3, 5, 7]:
        conj_mk = {3: 2, 5: 4, 7: 6}[mk]
        gal = f"conj({conj_mk})"
    print(f"rho_{mk:>2}  {dims[cy]:>4}  rho_{cy:>2}  {dims[cy]:>4}  {gal:>8}")

# ================================================================
# McKay Solution 3 assignment (from exp323)
# ================================================================
# In McKay numbering:
# rho_1 -> e, rho_2 -> u, rho_3 -> t, rho_4 -> d
# rho_5 -> b, rho_6 -> tau, rho_7 -> mu, rho_8 -> s, rho_9 -> c
#
# Translating to Cayley numbering:
# e -> 0, u -> 1, t -> 7, d -> 2, b -> 8, tau -> 3, mu -> 6, s -> 4, c -> 5

mckay_sol3 = {
    'e': 0, 'mu': 6, 'tau': 3, 'u': 1, 'c': 5, 't': 7, 'd': 2, 's': 4, 'b': 8
}
assign_mckay = [mckay_sol3[f] for f in fermion_names]

# SCAN 5 optimal (free 4 params)
assign_optimal = [0, 3, 5, 4, 2, 8, 7, 1, 6]

print("\n" + "="*72)
print("ASSIGNMENT COMPARISON")
print("="*72)
print(f"\n{'Fermion':>8} {'n':>4} {'McKay':>8} {'Optimal':>8} {'Match':>6}")
print("-"*40)
for f in range(9):
    match = "YES" if assign_mckay[f] == assign_optimal[f] else ""
    print(f"{fermion_names[f]:>8} {mass_exp[f]:>4} rho_{assign_mckay[f]:>2} rho_{assign_optimal[f]:>2} {match:>6}")

# ================================================================
# Compute RMS for both assignments with free amplitudes (lstsq)
# ================================================================
def compute_optimal_rms(assignment, name=""):
    """Given assignment, compute optimal c_m by lstsq and return RMS."""
    M = np.zeros((9, 4))
    for f in range(9):
        k = assignment[f]
        for idx, m_idx in enumerate(even_modes):
            M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
    result = np.linalg.lstsq(M, delta_exp, rcond=None)
    c_opt = result[0]
    pred = M @ c_opt
    rms = np.sqrt(np.mean((delta_exp - pred)**2))
    return rms, c_opt, pred

rms_mckay, c_mckay, pred_mckay = compute_optimal_rms(assign_mckay, "McKay Sol3")
rms_opt, c_opt, pred_opt = compute_optimal_rms(assign_optimal, "Optimal")

print(f"\n{'='*72}")
print(f"FREE AMPLITUDES (4 params, lstsq)")
print(f"{'='*72}")

print(f"\nMcKay Solution 3: RMS = {rms_mckay:.6f}")
print(f"  c_m = [{', '.join(f'{x:.6f}' for x in c_mckay)}]")
print(f"  {'Fermion':>6} {'pred':>8} {'exp':>8} {'err':>8}")
for f in range(9):
    print(f"  {fermion_names[f]:>6} {pred_mckay[f]:>+8.4f} {delta_exp[f]:>+8.4f} {pred_mckay[f]-delta_exp[f]:>+8.4f}")

print(f"\nOptimal (SCAN 5): RMS = {rms_opt:.6f}")
print(f"  c_m = [{', '.join(f'{x:.6f}' for x in c_opt)}]")
print(f"  {'Fermion':>6} {'pred':>8} {'exp':>8} {'err':>8}")
for f in range(9):
    print(f"  {fermion_names[f]:>6} {pred_opt[f]:>+8.4f} {delta_exp[f]:>+8.4f} {pred_opt[f]-delta_exp[f]:>+8.4f}")

# ================================================================
# Check ALL 4 McKay solutions
# ================================================================
print(f"\n{'='*72}")
print("ALL McKay Solutions (from exp323)")
print("="*72)

# The 4 solutions differ only in {tau, b, t} placement on the branch
# Main chain is FIXED: e->rho_1, u->rho_2, d->rho_4, mu/s->rho_7/8, c->rho_9
# That gives in Cayley: e->0, u->1, d->2, mu->6 or 4, s->4 or 6, c->5
# Wait - mu and s are degenerate (both n=11), so mu<->s swap doesn't matter

# Main chain fixed:
# mckay: rho_1->e, rho_2->u, rho_4->d, rho_7/rho_8->{mu,s}, rho_9->c
# In Cayley: e->0, u->1, d->2, {mu,s}->{6,4}, c->5

# Branch varies: {tau, b, t} on {rho_3, rho_5, rho_6}
# In Cayley: {rho_7, rho_8, rho_3}
# Note: rho_6(dim 4) is already taken by mu, so the branch from rho_9:
# Actually, let me re-examine the McKay graph more carefully.

# McKay graph edges (from mckay_yukawa.md):
# rho_1(1)--rho_2(2)--rho_4(3)--rho_7(4')--rho_8(5)--rho_9(6)--rho_5(3')
#                        |
#                      rho_6(4)
#                        |
#                      rho_3(2')

# So the branch at rho_4 (NOT rho_9!) leads to rho_6->rho_3
# And the main chain goes: rho_1-rho_2-rho_4-rho_7-rho_8-rho_9-rho_5
# Additional branch: rho_4-rho_6-rho_3

# Wait, this is confusing. Let me re-read more carefully.
# From mckay_yukawa.md:
# rho_1(1)--rho_2(2)--rho_4(3)--rho_7(4')--rho_8(5)--rho_9(6)--rho_5(3')
#                        |
#                      rho_6(4)
#                        |
#                      rho_3(2')
# So the branch is at rho_4 (dim 3), NOT at rho_9 (dim 6)!

# But MEMORY.md says "branch at rho_9 (dim 6=b1)". Let me look at it from dims:
# Affine E8: dims 1,2,3,4,5,6,4,2,3 with sum=30
# The node of dim 6 has degree 3 in the Dynkin diagram
# So branch at dim 6 node.

# In the graph above, rho_9(6) connects to: rho_8(5) and rho_5(3')
# That's only degree 2 for rho_9. But in affine E8, the dim-6 node has degree 3.

# Let me reconsider. The affine E8 Dynkin diagram:
#  1 - 2 - 3 - 4 - 5 - 6 - 4 - 2
#                      |
#                      3
# The dim-6 node connects to: dim-5, dim-4, dim-3 (branch)
# So the graph should be:
# rho_0(1) - rho_1(2) - rho_2(3) - rho_3(4) - rho_4(5) - rho_5(6) - rho_6(4) - rho_7(2)
#                                                           |
#                                                         rho_8(3)

# In Cayley numbering! This matches our CG_self matrix.

# Now in McKay numbering from mckay_yukawa.md, the graph has branch at rho_4:
# rho_1(1)--rho_2(2)--rho_4(3)--rho_7(4')--rho_8(5)--rho_9(6)--rho_5(3')
#                        |
#                      rho_6(4)
#                        |
#                      rho_3(2')

# So rho_9(6) connects to rho_8(5) and rho_5(3') only - that's degree 2
# But rho_4(3) connects to rho_2(2), rho_7(4'), and rho_6(4) - that's degree 3
# Wait no, that means the BRANCH is at rho_4 (dim 3)?!

# That can't be right for affine E8 where the branch should be at the dim-6 node.

# I think there's a numbering confusion. Let me look at it from pure E8 structure.
# The affine E8 Dynkin diagram (label = dim of 2I irrep):
#
#  1 -- 2 -- 3 -- 4 -- 5 -- 6 -- 4' -- 2'
#                            |
#                            3'
#
# The branch is at dim 6. That means:
# dim-6 connects to: dim-5, dim-4', and dim-3'
# This is 3 neighbors.

# Now mapping to McKay numbering from the graph in mckay_yukawa.md,
# rho_9 has dim 6 but only 2 neighbors (rho_8 and rho_5).
# Unless the graph in mckay_yukawa.md is WRONG.

# Actually wait, re-reading more carefully:
# "Branch at rho_9 (dim 6=b1), legs (1,2,5)"
# Legs of length 1, 2, 5.
# So rho_9(6) has THREE neighbors:
# Leg of length 5: rho_9-rho_8-rho_7-rho_4-rho_2-rho_1 (dims 6,5,4',3,2,1)
# Leg of length 2: rho_9-rho_5-??
# Hmm, but rho_5 has dim 3' and is an endpoint in the graph shown.

# I think the graph from mckay_yukawa.md is drawn misleadingly. Let me reconstruct it:
# The 8 edges of affine E8 connect:
# 1-2, 2-3, 3-4, 4-5, 5-6, 6-4', 6-3', 4'-2'
#
# In McKay numbering:
# rho_1(1)-rho_2(2)-rho_4(3)-rho_7(4')-rho_8(5)-rho_9(6)
#                                                     |--rho_5(3')
#                                                     |--rho_6(4)--rho_3(2')
#
# Wait, that gives rho_9 connecting to rho_8, rho_5, and rho_6 => degree 3 for rho_9

# So the CORRECT graph is:
# rho_1(1)--rho_2(2)--rho_4(3)--rho_7(4')--rho_8(5)--rho_9(6)
#                                                      |    |
#                                                 rho_5(3') |
#                                                      rho_6(4)--rho_3(2')

# rho_9(6) connects to: rho_8(5), rho_5(3'), rho_6(4) => degree 3. YES!
# And rho_6(4) connects to: rho_9(6), rho_3(2') => degree 2

# The graph in mckay_yukawa.md was drawn incorrectly (had the branch at rho_4).
# The CORRECT branch is at rho_9(6) as stated in MEMORY.md.

# OK so the correct McKay graph in Cayley numbering:
# rho_0(1) -- rho_1(2) -- rho_2(3) -- rho_3(4) -- rho_4(5) -- rho_5(6)
#                                                                |  |  \
#                                                          rho_8(3) |   rho_6(4)--rho_7(2)
# Wait, this doesn't quite work either. Let me just use dims.

# Affine E8 dims with branch at 6:
# Chain: 1--2--3--4--5--6
# Branch at 6: 6--3', 6--4'--2'
# Legs from 6: length 5 (to 1), length 1 (to 3'), length 2 (to 2')

# In Cayley numbering (dims [1,2,3,4,5,6,4,2,3]):
# 1=rho_0, 2=rho_1, 3=rho_2, 4=rho_3, 5=rho_4, 6=rho_5
# 3'=rho_8, 4'=rho_6, 2'=rho_7
# Edges: 0-1, 1-2, 2-3, 3-4, 4-5, 5-8, 5-6, 6-7

# In McKay numbering:
# 1=rho_1, 2=rho_2, 3=rho_4, 4=rho_7, 5=rho_8, 6=rho_9
# 3'=rho_5, 4'=rho_6, 2'=rho_3

# So the McKay graph edges are:
# rho_1-rho_2, rho_2-rho_4, rho_4-rho_7, rho_7-rho_8, rho_8-rho_9
# rho_9-rho_5, rho_9-rho_6, rho_6-rho_3

# This is consistent with branch at rho_9 (dim 6).
# And the 3 legs from rho_9 are:
#  - rho_9--rho_8--rho_7--rho_4--rho_2--rho_1 (length 5)
#  - rho_9--rho_5 (length 1)
#  - rho_9--rho_6--rho_3 (length 2)

# Good. Now I can write the McKay graph in Cayley numbering.

# McKay graph adjacency (in Cayley numbering):
mckay_edges_cayley = [
    (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),  # main chain
    (5, 8),  # branch leg 1 (length 1)
    (5, 6), (6, 7),  # branch leg 2 (length 2)
]

print("\nMcKay graph edges (Cayley numbering):")
for i, j in mckay_edges_cayley:
    print(f"  rho_{i}(dim {dims[i]}) -- rho_{j}(dim {dims[j]})")

# Graph distance matrix
mckay_dist = np.full((9, 9), 99)
for i in range(9):
    mckay_dist[i, i] = 0
for i, j in mckay_edges_cayley:
    mckay_dist[i, j] = 1
    mckay_dist[j, i] = 1
# Floyd-Warshall
for k in range(9):
    for i in range(9):
        for j in range(9):
            if mckay_dist[i, k] + mckay_dist[k, j] < mckay_dist[i, j]:
                mckay_dist[i, j] = mckay_dist[i, k] + mckay_dist[k, j]

print("\nMcKay graph distance matrix (Cayley numbering):")
print("      ", "  ".join(f"r{k}" for k in range(9)))
for i in range(9):
    print(f"  r{i}: " + "  ".join(f"{mckay_dist[i,j]:2d}" for j in range(9)))

# ================================================================
# Test exp323 Solution 3 and all 4 solutions
# ================================================================
print(f"\n{'='*72}")
print("McKay Solutions with weighted graph distances")
print("="*72)

# Solution 3 edge weights: [3, 2, 6, 0, 5, 3, 1, 9]
# These are for edges in McKay order:
# rho_1-rho_2: 3, rho_2-rho_4: 2, rho_4-rho_7: 6, rho_7-rho_8: 0
# rho_8-rho_9: 5, rho_9-rho_5: 3, rho_9-rho_6: 1, rho_6-rho_3: 9

# In Cayley numbering:
# rho_0-rho_1: 3, rho_1-rho_2: 2, rho_2-rho_3: 6, rho_3-rho_4: 0
# rho_4-rho_5: 5, rho_5-rho_8: 3, rho_5-rho_6: 1, rho_6-rho_7: 9

sol3_weights = {(0,1): 3, (1,2): 2, (2,3): 6, (3,4): 0,
                (4,5): 5, (5,8): 3, (5,6): 1, (6,7): 9}

# Build weighted distance matrix
wdist = np.full((9, 9), 999.0)
for i in range(9):
    wdist[i, i] = 0
for (i, j), w in sol3_weights.items():
    wdist[i, j] = w
    wdist[j, i] = w
for k in range(9):
    for i in range(9):
        for j in range(9):
            if wdist[i, k] + wdist[k, j] < wdist[i, j]:
                wdist[i, j] = wdist[i, k] + wdist[k, j]

print("\nSolution 3 weighted distances from rho_0 (Cayley):")
for k in range(9):
    print(f"  d(rho_0, rho_{k}) = {wdist[0, k]:.0f} (dim={dims[k]})")

# Which fermion has n = d(rho_0, rho_k)?
print("\nMcKay Sol3: fermion assignment by weighted distance = mass exponent n:")
for k in range(9):
    d = int(wdist[0, k])
    matches = [fermion_names[f] for f in range(9) if mass_exp[f] == d]
    print(f"  rho_{k} (dim {dims[k]}): d=n={d} -> {matches if matches else 'NO MATCH'}")

# ================================================================
# FULL RMS COMPARISON
# ================================================================
print(f"\n{'='*72}")
print("RMS COMPARISON: McKay vs Optimal Assignments")
print("="*72)

# All possible McKay solutions differ in {tau, b, t} on branch
# Main chain (FIXED for all solutions):
# rho_0 -> e(0), rho_1 -> u(3), rho_2 -> d(5),
# rho_3/rho_4 degenerate -> {mu(11), s(11)}
# rho_5 is the branch center

# Actually, the main chain distances:
# d(0,0)=0 -> e(n=0)
# d(0,1)=3 -> u(n=3)
# d(0,2)=5 -> d(n=5)
# d(0,3)=11 -> mu or s (n=11, degenerate)
# d(0,4)=11 -> mu or s (n=11, degenerate)
# d(0,5)=16 -> c(n=16)
# Branch:
# d(0,8)=19 -> b(n=19) [via rho_5-rho_8, w=3, total=16+3=19]
# d(0,6)=17 -> tau(n=17) [via rho_5-rho_6, w=1, total=16+1=17]
# d(0,7)=26 -> t(n=26) [via rho_5-rho_6-rho_7, w=1+9=10, total=16+10=26]

# So McKay Sol3 (in Cayley numbering):
print("\nMcKay Sol3 assignment (Cayley numbering):")
print("  e->0, u->1, d->2, mu->3(or 4), s->4(or 3), c->5, tau->6, t->7, b->8")

# mu and s are degenerate (both n=11), so the assignment mu->rho_3/rho_4 is ambiguous
# Let's test both permutations:
# Option A: mu->3, s->4
assign_mckay_A = [0, 3, 6, 1, 5, 7, 2, 4, 8]  # e,mu,tau,u,c,t,d,s,b
# Option B: mu->4, s->3
assign_mckay_B = [0, 4, 6, 1, 5, 7, 2, 3, 8]

# Also, the 3 other solutions from exp323 swap {tau, b, t} on {rho_6, rho_7, rho_8}
# But the branch weights would change. Let me just try all 6 permutations of {tau, b, t}
# on the 3 branch nodes {rho_6, rho_7, rho_8}:
branch_nodes = [6, 7, 8]
branch_fermions_idx = [2, 5, 8]  # tau=2, t=5, b=8

# Base: e->0, mu->3, u->1, c->5, d->2, s->4
base_assign = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}  # fermion_idx -> cayley_rho

print(f"\n{'Assignment':>40} {'RMS(4p)':>10} {'RMS(quark)':>12}")
print("-"*70)

all_tested = []

for tau_node, t_node, b_node in permutations(branch_nodes):
    # Also try mu<->s swap on {3, 4}
    for mu_node, s_node in [(3, 4), (4, 3)]:
        assign = [0] * 9
        assign[0] = 0  # e
        assign[1] = mu_node  # mu
        assign[2] = tau_node  # tau
        assign[3] = 1  # u
        assign[4] = 5  # c
        assign[5] = t_node  # t
        assign[6] = 2  # d
        assign[7] = s_node  # s
        assign[8] = b_node  # b

        rms, c, pred = compute_optimal_rms(assign)
        # Quark-only RMS
        quark_idx = [3, 4, 5, 6, 7, 8]
        rms_q = np.sqrt(np.mean((delta_exp[quark_idx] - pred[quark_idx])**2))

        label = f"mu->{mu_node},s->{s_node},tau->{tau_node},t->{t_node},b->{b_node}"
        all_tested.append((rms, rms_q, assign[:], label, c, pred))
        print(f"{label:>40} {rms:>10.4f} {rms_q:>12.4f}")

# Sort by total RMS
all_tested.sort(key=lambda x: x[0])
print(f"\nBest McKay-based assignment:")
rms, rms_q, assign, label, c, pred = all_tested[0]
print(f"  {label}")
print(f"  Total RMS = {rms:.6f}, Quark RMS = {rms_q:.6f}")
print(f"  c_m = [{', '.join(f'{x:.6f}' for x in c)}]")
print(f"  {'Fermion':>6} {'pred':>8} {'exp':>8} {'err':>8}")
for f in range(9):
    print(f"  {fermion_names[f]:>6} {pred[f]:>+8.4f} {delta_exp[f]:>+8.4f} {pred[f]-delta_exp[f]:>+8.4f}")

# ================================================================
# Compare with SCAN 5 optimal
# ================================================================
print(f"\n{'='*72}")
print("COMPARISON: Best McKay vs Optimal")
print("="*72)

rms_mckay_best = all_tested[0][0]
rms_mckay_q = all_tested[0][1]

print(f"\n{'':>30} {'Total RMS':>12} {'Quark RMS':>12}")
print("-"*60)
print(f"{'McKay best':>30} {rms_mckay_best:>12.4f} {rms_mckay_q:>12.4f}")

rms_o, c_o, pred_o = compute_optimal_rms(assign_optimal)
quark_idx = [3, 4, 5, 6, 7, 8]
rms_oq = np.sqrt(np.mean((delta_exp[quark_idx] - pred_o[quark_idx])**2))
print(f"{'SCAN5 optimal':>30} {rms_o:>12.4f} {rms_oq:>12.4f}")
print(f"{'Bare formula':>30} {0.247:>12.4f} {'':>12}")

# ================================================================
# Key question: WHY is the optimal assignment different?
# ================================================================
print(f"\n{'='*72}")
print("STRUCTURAL ANALYSIS: Why assignments differ")
print("="*72)

# The CG_self matrix determines which even modes contribute to each irrep
print("\nCG_self for even modes (rows = Cayley irreps, cols = even modes 2,4,6,8):")
print(f"{'rho_k':>6} {'dim':>4}  m=2  m=4  m=6  m=8  CG_pattern")
print("-"*55)
for k in range(9):
    cg_vec = [CG_self[k, m] for m in even_modes]
    pattern = tuple(cg_vec)
    print(f"rho_{k:>2} {dims[k]:>4}  {cg_vec[0]:>3}  {cg_vec[1]:>3}  {cg_vec[2]:>3}  {cg_vec[3]:>3}  {pattern}")

# Identify UNIQUE CG patterns
patterns = {}
for k in range(9):
    pat = tuple(CG_self[k, m] for m in even_modes)
    if pat not in patterns:
        patterns[pat] = []
    patterns[pat].append(k)

print(f"\nDistinct CG patterns: {len(patterns)}")
for pat, nodes in patterns.items():
    node_str = ", ".join(f"rho_{k}(d={dims[k]})" for k in nodes)
    print(f"  {pat} -> {node_str}")

# For the VEV mechanism, fermions with SAME CG pattern get SAME correction
# So if two fermions need DIFFERENT corrections but have same CG pattern,
# they CANNOT be at the same pattern group.
print("\nImplication: fermions at same-pattern nodes get IDENTICAL corrections.")
print("McKay places rho_3 & rho_6 together (pattern (1,1,1,1))")
print("  -> McKay: tau(rho_3) and t/b at rho_6/7/8")
print("  -> If tau and t need different corrections, this is a problem.")

# ================================================================
# What if we CONSTRAIN to McKay-compatible but allow branch freedom?
# ================================================================
print(f"\n{'='*72}")
print("McKay-COMPATIBLE assignment with 2-param function forms")
print("="*72)

# For each of the 12 McKay-type assignments, try 2-param forms
alt_sign = np.array([-1, 1, -1, 1])
forms = {
    'A*alt+B': np.column_stack([alt_sign, np.ones(4)]),
    'A*lam+B': np.column_stack([lam, np.ones(4)]),
    'A/lam+B/d': np.column_stack([1.0/lam, 1.0/d_m]),
}

for form_name, F in forms.items():
    print(f"\n--- {form_name} ---")
    best_mck = (999, None, None)
    for rms, rms_q, assign, label, c, pred in all_tested:
        M = np.zeros((9, 4))
        for f in range(9):
            k = assign[f]
            for idx, m_idx in enumerate(even_modes):
                M[f, idx] = fermion_T3[f] * fermion_Nc[f] * CG_self[k, m_idx] * dims[m_idx]
        MF = M @ F
        result = np.linalg.lstsq(MF, delta_exp, rcond=None)
        params = result[0]
        pred2 = MF @ params
        rms2 = np.sqrt(np.mean((delta_exp - pred2)**2))
        if rms2 < best_mck[0]:
            best_mck = (rms2, label, params)
    print(f"  Best: RMS={best_mck[0]:.4f}, {best_mck[1]}, params={best_mck[2]}")

# ================================================================
# CRITICAL: Check the "effective CG vector" for each fermion
# ================================================================
print(f"\n{'='*72}")
print("EFFECTIVE CORRECTION VECTOR per fermion")
print("="*72)

print("\nFor McKay assignment (best):")
assign = all_tested[0][2]
c_best = all_tested[0][4]
print(f"  Assignment: {assign}")
print(f"  Amplitudes: {np.round(c_best, 4)}")
print(f"\n  {'Fermion':>6} {'T3':>5} {'Nc':>3} {'rho':>4} {'CG(2,4,6,8)':>20} {'eff_coupling':>14}")
for f in range(9):
    k = assign[f]
    cg = [CG_self[k, m]*dims[m] for m in even_modes]
    eff = sum(cg[i]*c_best[i] for i in range(4))
    print(f"  {fermion_names[f]:>6} {fermion_T3[f]:>+5.1f} {fermion_Nc[f]:>3} {k:>4}  {str([round(x,1) for x in cg]):>20} {eff:>+14.4f}")

print("\nFor SCAN5 optimal:")
c_opt = all_tested[0][4]  # recompute
_, c_opt_scan5, _ = compute_optimal_rms(assign_optimal)
print(f"  Assignment: {assign_optimal}")
print(f"  Amplitudes: {np.round(c_opt_scan5, 4)}")
print(f"\n  {'Fermion':>6} {'T3':>5} {'Nc':>3} {'rho':>4} {'CG(2,4,6,8)':>20} {'eff_coupling':>14}")
for f in range(9):
    k = assign_optimal[f]
    cg = [CG_self[k, m]*dims[m] for m in even_modes]
    eff = sum(cg[i]*c_opt_scan5[i] for i in range(4))
    print(f"  {fermion_names[f]:>6} {fermion_T3[f]:>+5.1f} {fermion_Nc[f]:>3} {k:>4}  {str([round(x,1) for x in cg]):>20} {eff:>+14.4f}")

# ================================================================
# GALOIS STRUCTURE: Can we connect assignments via Galois?
# ================================================================
print(f"\n{'='*72}")
print("GALOIS STRUCTURE of both assignments")
print("="*72)

# Galois conjugation maps:
# In Cayley: rho_1(2) <-> rho_7(2), rho_2(3) <-> rho_8(3), rho_3(4) <-> rho_6(4)
galois_pairs = [(1, 7), (2, 8), (3, 6)]

print("\nGalois pairs (Cayley numbering):")
for a, b in galois_pairs:
    print(f"  rho_{a}(dim {dims[a]}) <-> rho_{b}(dim {dims[b]})")

print(f"\nGalois consistency check:")
print(f"  McKay best:")
assign = all_tested[0][2]
for a, b in galois_pairs:
    fa = [fermion_names[f] for f in range(9) if assign[f] == a]
    fb = [fermion_names[f] for f in range(9) if assign[f] == b]
    print(f"    rho_{a}({fa}) <-> rho_{b}({fb})")

print(f"  SCAN5 optimal:")
for a, b in galois_pairs:
    fa = [fermion_names[f] for f in range(9) if assign_optimal[f] == a]
    fb = [fermion_names[f] for f in range(9) if assign_optimal[f] == b]
    print(f"    rho_{a}({fa}) <-> rho_{b}({fb})")

# Check if Galois pairs have conjugate mass exponents
print(f"\n  Mass exponent Galois pairs (n1 + n2 = ?):")
for a, b in galois_pairs:
    for label, asgn in [("McKay", all_tested[0][2]), ("Optimal", assign_optimal)]:
        fa = [f for f in range(9) if asgn[f] == a]
        fb = [f for f in range(9) if asgn[f] == b]
        if fa and fb:
            na = mass_exp[fa[0]]
            nb = mass_exp[fb[0]]
            print(f"    {label}: rho_{a}({fermion_names[fa[0]]},n={na}) + rho_{b}({fermion_names[fb[0]]},n={nb}) = {na+nb}")

print("\n\nDone.")
