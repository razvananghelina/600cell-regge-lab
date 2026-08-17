"""
EXP-158: Can alpha = 3/2 be justified from 600-cell geometry?

STG with alpha=3/2 gives beta=1 (GR). Can we derive 3/2 from the 600-cell?

Known: alpha=1/2 from 360/720 angular edges.
Question: What gives the factor 3?

Candidates:
1. 600-cell/120-cell duality factor = 3 (degree, diameter, eigenvalues all ratio 3)
2. Dimension of S^3 = 3
3. Number of spatial dimensions = 3
4. Dihedral angle 2*pi/5 has 5 tetrahedra per edge, and 5/... something
5. Vertex figure icosahedron has F/E = 20/30 = 2/3
6. Some spectral ratio
"""

import numpy as np
from itertools import combinations

print("=" * 70)
print("EXP-158: Alpha = 3/2 from 600-Cell Geometry")
print("=" * 70)

# =====================================================================
# Generate 600-cell
# =====================================================================
phi = (1 + np.sqrt(5)) / 2
phi_inv = 1 / phi

vertices = []
# 8 vertices: permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0, 0, 0, 0]
        v[i] = s
        vertices.append(v)

# 16 vertices: (+-1/2, +-1/2, +-1/2, +-1/2)
for s0 in [1, -1]:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                vertices.append([s0/2, s1/2, s2/2, s3/2])

# 96 vertices: even permutations of (+-phi/2, +-1/2, +-phi_inv/2, 0)
vals = [phi/2, 0.5, phi_inv/2, 0]
even_perms = []
for p in [[0,1,2,3],[0,2,3,1],[0,3,1,2],
          [1,0,3,2],[1,2,0,3],[1,3,2,0],
          [2,0,1,3],[2,1,3,0],[2,3,0,1],
          [3,0,2,1],[3,1,0,2],[3,2,1,0]]:
    even_perms.append(p)

for perm in even_perms:
    base = [vals[perm[i]] for i in range(4)]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = [0, 0, 0, 0]
                for i in range(4):
                    if base[i] == 0:
                        v[i] = 0
                    else:
                        signs = [s0, s1, s2]
                        nonzero_idx = sum(1 for j in range(i) if base[j] != 0)
                        if nonzero_idx < 3:
                            v[i] = base[i] * signs[nonzero_idx]
                        else:
                            v[i] = base[i]
                vertices.append(v)

# Deduplicate
verts = []
for v in vertices:
    v = np.array(v)
    is_dup = False
    for w in verts:
        if np.linalg.norm(v - w) < 1e-10:
            is_dup = True
            break
    if not is_dup:
        verts.append(v)

verts = np.array(verts)
print(f"\nVertices: {len(verts)}")

# Build edges
dots = verts @ verts.T
edge_dot = np.cos(np.pi/5)  # = phi/2 ~ 0.809
edges = []
adj = [[] for _ in range(len(verts))]
for i in range(len(verts)):
    for j in range(i+1, len(verts)):
        if abs(dots[i, j] - edge_dot) < 1e-6:
            edges.append((i, j))
            adj[i].append(j)
            adj[j].append(i)

print(f"Edges: {len(edges)}")
print(f"Degree: {len(adj[0])}")

# =====================================================================
# PART 1: The factor 3 from 120-cell duality
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: Factor 3 from 120-cell Duality")
print("=" * 70)

print("""
120-cell / 600-cell ratios (all exactly 3):
  degree:    12 / 4  = 3
  diameter:  15 / 5  = 3
  eigenvals: 27 / 9  = 3

If alpha_600 = 1/2 (from 600-cell edge fraction),
then alpha_120 = 3 * alpha_600 = 3/2?

Interpretation: The 600-cell gives the "angular fraction" 1/2.
The 120-cell (dual, representing cells/volume) amplifies this by 3.
The physical metric uses the VOLUME element, not the edge element.
In 3+1 dimensions, the volume element has 3 spatial factors.

alpha_physical = (spatial_dim) * alpha_edge = 3 * 1/2 = 3/2
""")

# =====================================================================
# PART 2: Dimension counting
# =====================================================================
print("=" * 70)
print("PART 2: Dimension Counting Arguments")
print("=" * 70)

# The STG metric: g_ang = C^(2*alpha) * r^2
# In d spatial dimensions, the solid angle element is r^(d-1) dOmega
# The metric on a sphere of radius r in curved space:
# ds^2_sphere = C^(2*alpha) * r^2 * dOmega^2
# The area of the sphere: A = C^(2*alpha*(d-1)) * Omega_{d-1} * r^(d-1)
# Wait, that's not right. For a 2-sphere in 3+1d:
# ds^2_sphere = C^(2*alpha) * r^2 * (dtheta^2 + sin^2(theta) dphi^2)
# Area = integral C^(2*alpha) r^2 sin(theta) dtheta dphi = 4*pi * C^(2*alpha) * r^2

# The areal radius: R_areal = C^alpha * r
# In d spatial dimensions: R_areal = C^(alpha * (d-1)/(d-1)) * r = C^alpha * r
# No, the areal radius doesn't depend on d in this way.

# Different argument: the 600-cell lives on S^3 (3-dimensional sphere).
# The metric on S^3 has 3 angular dimensions.
# Each dimension contributes alpha_edge = 1/2.
# Total: 3 * 1/2 = 3/2?

# In the STG metric, g_ang = C^(2*alpha) * r^2
# This multiplies the ENTIRE angular part (dtheta^2 + sin^2 dphi^2)
# If each angular dimension independently has factor C^(2*1/2) = C,
# then two angular dimensions give C^2, i.e., alpha=1 for 3+1d.
# Three angular dimensions (S^3) give C^3, i.e., alpha=3/2.

print("""
Argument: Each spatial dimension of S^3 contributes alpha=1/2.

In 3+1 dimensions (our spacetime):
  - 2 angular dimensions on each 2-sphere
  - Each contributes C^(2*1/2) = C to the metric
  - Total angular: C^2 => alpha_eff = 1

In the 600-cell (living on S^3):
  - 3 angular dimensions on S^3
  - Each contributes C^(2*1/2) = C
  - Total angular: C^3 => alpha_eff = 3/2

This gives beta = 5/2 - 3/2 = 1 = GR!

BUT: This argument assumes a specific dimensional reduction S^3 -> R^3.
The number of angular dimensions in the TARGET space is 2 (not 3).
So should we use 2 * 1/2 = 1 instead?
""")

# =====================================================================
# PART 3: Check alpha=1 (2 angular dims)
# =====================================================================
print("=" * 70)
print("PART 3: Alpha Values and Their Predictions")
print("=" * 70)

from sympy import Rational, pi as spi, sqrt as ssqrt, Symbol, simplify, series

alpha_sym = Symbol('alpha')
beta_formula = Rational(5, 2) - alpha_sym

for alpha_val, label, justification in [
    (Rational(1, 2), "1/2", "edge fraction 360/720"),
    (1, "1", "2 angular dims x 1/2 (conformal metric)"),
    (Rational(3, 2), "3/2", "3 dims of S^3 x 1/2, OR 120-cell factor"),
    (2, "2", "4 dims of R^4 x 1/2"),
]:
    beta_val = Rational(5, 2) - alpha_val
    gamma_val = 1
    prec_factor = (2 + 2*gamma_val - beta_val) / 3
    print(f"\nalpha = {label} ({justification}):")
    print(f"  beta = {beta_val}, gamma = 1")
    print(f"  Precession factor = {prec_factor} (GR = 1)")
    print(f"  Deflection factor = {(1+gamma_val)/2} (GR = 1)")
    if beta_val == 1:
        print(f"  *** MATCHES GR! ***")
    else:
        print(f"  Excluded by Solar System tests")

# =====================================================================
# PART 4: Angular fraction on the 600-cell - more careful analysis
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: Angular Fractions - Detailed Analysis")
print("=" * 70)

# BFS from vertex 0
from collections import deque

def bfs_shells(adj_list, start):
    n = len(adj_list)
    dist = [-1] * n
    dist[start] = 0
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for w in adj_list[v]:
            if dist[w] == -1:
                dist[w] = dist[v] + 1
                queue.append(w)
    return dist

# Classify edges as angular (same shell) or radial (different shell)
dist_0 = bfs_shells(adj, 0)
shells = {}
for i, d in enumerate(dist_0):
    if d not in shells:
        shells[d] = []
    shells[d].append(i)

print(f"\nShell sizes from vertex 0: {[len(shells[d]) for d in sorted(shells)]}")

angular = 0
radial = 0
angular_per_shell = {}
radial_per_shell = {}

for i, j in edges:
    di, dj = dist_0[i], dist_0[j]
    if di == dj:
        angular += 1
        angular_per_shell[di] = angular_per_shell.get(di, 0) + 1
    else:
        radial += 1
        d_min = min(di, dj)
        radial_per_shell[d_min] = radial_per_shell.get(d_min, 0) + 1

print(f"Angular edges: {angular}, Radial edges: {radial}")
print(f"Fraction angular: {angular}/{angular+radial} = {angular/(angular+radial):.6f}")

print(f"\nPer-shell angular edges:")
for d in sorted(angular_per_shell):
    total_d = angular_per_shell.get(d, 0) + radial_per_shell.get(d, 0) + radial_per_shell.get(d-1, 0) if d > 0 else 0
    print(f"  Shell {d}: {angular_per_shell.get(d, 0)} angular")

# Now compute angular fraction for FACES (triangles)
print("\n--- Triangular faces ---")
# Find all triangles
triangles = []
for i in range(len(verts)):
    for j in adj[i]:
        if j > i:
            for k in adj[j]:
                if k > j and k in adj[i]:
                    triangles.append((i, j, k))

print(f"Total triangles: {len(triangles)}")

# Classify triangles: how many vertices in each shell?
tri_types = {}
for i, j, k in triangles:
    di, dj, dk = dist_0[i], dist_0[j], dist_0[k]
    key = tuple(sorted([di, dj, dk]))
    tri_types[key] = tri_types.get(key, 0) + 1

angular_tri = 0
mixed_tri = 0
for key, count in sorted(tri_types.items()):
    if key[0] == key[1] == key[2]:
        angular_tri += count
        label = "ANGULAR"
    else:
        mixed_tri += count
        label = "mixed"
    print(f"  ({key[0]},{key[1]},{key[2]}): {count} [{label}]")

print(f"\nAngular triangles: {angular_tri}")
print(f"Mixed triangles: {mixed_tri}")
print(f"Fraction angular: {angular_tri}/{angular_tri+mixed_tri} = {angular_tri/(angular_tri+mixed_tri):.6f}")

# Classify tetrahedra
print("\n--- Tetrahedra (cells) ---")
# Find all tetrahedra
tetrahedra = []
for i in range(len(verts)):
    for j in adj[i]:
        if j > i:
            common_ij = set(adj[i]) & set(adj[j])
            for k in common_ij:
                if k > j:
                    common_ijk = common_ij & set(adj[k])
                    for l in common_ijk:
                        if l > k:
                            tetrahedra.append((i, j, k, l))

print(f"Total tetrahedra: {len(tetrahedra)}")

tet_types = {}
angular_tet = 0
for tet in tetrahedra:
    dists = tuple(sorted([dist_0[v] for v in tet]))
    tet_types[dists] = tet_types.get(dists, 0) + 1
    if dists[0] == dists[-1]:
        angular_tet += 1

for key, count in sorted(tet_types.items()):
    label = "ANGULAR" if key[0] == key[-1] else ""
    print(f"  {key}: {count} {label}")

print(f"\nAngular tetrahedra: {angular_tet}")
print(f"Total tetrahedra: {len(tetrahedra)}")
print(f"Fraction angular tet: {angular_tet/(len(tetrahedra)):.6f}")

# =====================================================================
# PART 5: The hierarchy of angular fractions
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: Hierarchy of Angular Fractions")
print("=" * 70)

frac_edge = angular / len(edges)
frac_tri = angular_tri / len(triangles)
frac_tet = angular_tet / len(tetrahedra)

print(f"Angular fraction (edges, dim=1):    {angular}/{len(edges)} = {frac_edge:.6f}")
print(f"Angular fraction (faces, dim=2):    {angular_tri}/{len(triangles)} = {frac_tri:.6f}")
print(f"Angular fraction (cells, dim=3):    {angular_tet}/{len(tetrahedra)} = {frac_tet:.6f}")

# Check if (1/2)^k pattern
print(f"\nRatio tri/edge: {frac_tri/frac_edge:.6f}")
print(f"Ratio tet/tri: {frac_tet/frac_tri:.6f}")
print(f"Ratio tet/edge: {frac_tet/frac_edge:.6f}")

# Key insight: if alpha = sum of angular fractions by dimension...
alpha_sum = frac_edge + frac_tri + frac_tet
alpha_prod = frac_edge * frac_tri * frac_tet
print(f"\nSum of angular fractions: {alpha_sum:.6f}")
print(f"Product: {alpha_prod:.6f}")
print(f"1 + frac_edge: {1 + frac_edge:.6f}")

# Alternative: alpha from the AVERAGE angular fraction
# In d-dim simplex, the k-dimensional angular fraction might give alpha_k
# Total alpha = integral over dimensions?

# =====================================================================
# PART 6: Topological argument for 3/2
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: Topological Arguments for alpha = 3/2")
print("=" * 70)

# The 600-cell lives on S^3. The metric on S^3 in terms of
# Hopf coordinates is:
# ds^2 = dpsi^2 + sin^2(psi)(dtheta^2 + sin^2(theta) dphi^2)
# = dpsi^2 + sin^2(psi) dOmega_2^2
#
# In the STG framework, C(r) modifies the radial direction.
# The angular part of the metric in 3+1D spacetime has 2 angular dimensions.
# But the 600-cell is on S^3, which has 3 "angular" dimensions.
#
# Key question: does the metric exponent alpha count:
# (a) angular dimensions of the TARGET manifold (2 for S^2 in 3+1D) => alpha = 1
# (b) angular dimensions of the SOURCE manifold S^3 (3 for S^3) => alpha = 3/2
# (c) edge fraction directly (1/2) => alpha = 1/2

# The Euler characteristic tells us something:
chi_S3 = 0  # chi(S^3) = 0
chi_600 = 120 - 720 + 1200 - 600  # V - E + F - C = 0
print(f"chi(600-cell) = {chi_600} (= chi(S^3) = 0)")

# f-vector
V, E, F, C_cells = 120, 720, 1200, 600
print(f"f-vector: V={V}, E={E}, F={F}, C={C_cells}")
print(f"E/V = {E/V} = degree/2 = 6")
print(f"F/E = {F/E} = {F/E:.6f}")
print(f"C/F = {C_cells/F} = {C_cells/F:.6f}")
print(f"F/V = {F/V} = 10")
print(f"C/V = {C_cells/V} = 5")

# Dehn-Sommerville relations for 3-sphere simplicial complex:
# f_0 - f_1 + f_2 - f_3 = 0 (chi = 0)
# 2*f_1 = 3*f_2 - ... no, for simplicial:
# Each triangle has 3 edges => 3*F = 2*(edges per triangle counted with sharing)
# Actually: each edge is in exactly 5 triangles (from 5 tet per edge, each tet has 3 triangles thru edge... no)
# Let me compute directly.

# Edges per triangle (each triangle has 3 edges)
# Sum = 3*F = 3*1200 = 3600
# Each edge is in how many triangles?
edge_in_tri = 3 * F / E
print(f"\nAvg triangles per edge: 3*F/E = {edge_in_tri}")

# Triangles per tetrahedron (each tet has 4 triangles)
# 4*C = 4*600 = 2400
tri_in_tet = 4 * C_cells / F
print(f"Avg tetrahedra per triangle: 4*C/F = {tri_in_tet}")

# Edges per tetrahedron (each tet has 6 edges)
# 6*C = 6*600 = 3600
tet_per_edge = 6 * C_cells / E
print(f"Avg tetrahedra per edge: 6*C/E = {tet_per_edge}")
print(f"(= dihedral 2*pi / (2*pi/5) = 5)")

# =====================================================================
# PART 7: Direct metric derivation attempt
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: Direct Metric from Regge-like Counting")
print("=" * 70)

# In Regge calculus, the Einstein-Hilbert action is:
# S = sum_hinges (deficit_angle * area_of_hinge)
# In 3+1D, hinges are 2D triangles.
# In 3D (spatial), hinges are 1D edges.
#
# For the 600-cell (3-dimensional simplicial complex on S^3):
# S_Regge = sum_edges (deficit_angle_e * length_e)
#
# The dihedral angle of the 600-cell is 2*pi/5 = 72 degrees.
# With 5 tetrahedra around each edge, total angle = 5 * 72 = 360 = 2*pi.
# Deficit angle = 2*pi - 5*(2*pi/5) = 0.
#
# This means the 600-cell IS FLAT (as it should be - it tiles S^3).

dihedral_600 = 2 * np.pi / 5
tets_per_edge = 5
total_angle = tets_per_edge * dihedral_600
deficit = 2 * np.pi - total_angle

print(f"Dihedral angle: {np.degrees(dihedral_600):.2f} deg = 2*pi/5")
print(f"Tetrahedra per edge: {tets_per_edge}")
print(f"Total angle around edge: {np.degrees(total_angle):.2f} deg")
print(f"Deficit angle: {deficit:.10f} (should be 0 for S^3 tiling)")

print("""
The 600-cell has ZERO deficit angle at every edge.
This means it perfectly tiles S^3 - it IS flat (in the S^3 sense).
There is no intrinsic curvature FROM the triangulation itself.

The curvature of S^3 is "stored" in the EMBEDDING, not in deficit angles.
This is the key distinction: Regge calculus on the 600-cell gives
FLAT space (S^3 with its intrinsic curvature, not extra curvature).

For gravity, we need to DEFORM the 600-cell (change edge lengths
or remove cells) to create deficit angles.
""")

# =====================================================================
# PART 8: Summary and Assessment
# =====================================================================
print("=" * 70)
print("PART 8: Summary")
print("=" * 70)

print(f"""
ANGULAR FRACTION HIERARCHY:
  Edges (1-simplices):  {angular}/{len(edges)} = {frac_edge:.6f}  (topological invariant)
  Faces (2-simplices):  {angular_tri}/{len(triangles)} = {frac_tri:.6f}
  Cells (3-simplices):  {angular_tet}/{len(tetrahedra)} = {frac_tet:.6f}

ALPHA CANDIDATES AND THEIR PPN PREDICTIONS:
  alpha = 1/2 (edge fraction):    beta = 2,   precession = 2/3 GR  [EXCLUDED]
  alpha = 1   (conformal / 2-dim): beta = 3/2, precession = 5/6 GR  [EXCLUDED]
  alpha = 3/2 (3-dim / 120-cell): beta = 1,   precession = GR      [MATCHES!]
  alpha = 2   (4-dim):            beta = 1/2,  precession = 7/6 GR  [EXCLUDED]

ARGUMENTS FOR alpha = 3/2:
  [A1] 120-cell/600-cell duality factor 3: alpha = 3 * (1/2) = 3/2
  [A2] S^3 has 3 dimensions, each contributes 1/2: alpha = 3 * (1/2) = 3/2
  [A3] This is the UNIQUE value giving beta = 1 (GR-compatible)

ARGUMENTS AGAINST alpha = 3/2:
  [C1] The edge fraction 360/720 = 1/2 is a TOPOLOGICAL INVARIANT.
       The number 3/2 is not a topological invariant of the 600-cell.
  [C2] The factor 3 = dim(S^3) is external input, not from the graph.
  [C3] No rigorous derivation connects 600-cell to 3/2.

ASSESSMENT:
  alpha = 3/2 giving beta = 1 is TANTALIZING but currently SPECULATIVE.
  The 120-cell duality factor 3 and S^3 dimensionality both suggest 3/2,
  but neither constitutes a DERIVATION from the geometry alone.

  Status: SPECULATIVE (multiple hints, no proof)
""")
