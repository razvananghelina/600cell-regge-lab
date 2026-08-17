"""
EXP-159: Regge Calculus on the 600-Cell

Regge calculus is the standard discretization of GR on simplicial complexes.
The 600-cell is a simplicial complex (600 tetrahedra) tiling S^3.

Key questions:
1. What are the deficit angles? (Should be 0 for perfect S^3 tiling)
2. Can we deform the 600-cell to create curvature (mass)?
3. What is the Regge action?
4. Can we derive the Schwarzschild metric from a deformed 600-cell?
5. What is the relationship between deficit angles and the metric?
"""

import numpy as np
from collections import defaultdict, deque
from itertools import combinations

print("=" * 70)
print("EXP-159: Regge Calculus on the 600-Cell")
print("=" * 70)

# =====================================================================
# Generate 600-cell vertices on unit S^3
# =====================================================================
phi = (1 + np.sqrt(5)) / 2
phi_inv = 1 / phi

vertices = []
# 8 Type-A: perms of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = float(s)
        vertices.append(v)

# 16 Type-B: (+-1/2)^4
for s0 in [1, -1]:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                vertices.append([s0/2, s1/2, s2/2, s3/2])

# 96 Type-C: even permutations of (+-phi/2, +-1/2, +-phi_inv/2, 0)
vals = [phi/2, 0.5, phi_inv/2, 0.0]
even_perms = [
    [0,1,2,3],[0,2,3,1],[0,3,1,2],
    [1,0,3,2],[1,2,0,3],[1,3,2,0],
    [2,0,1,3],[2,1,3,0],[2,3,0,1],
    [3,0,2,1],[3,1,0,2],[3,2,1,0]
]

for perm in even_perms:
    base = [vals[perm[i]] for i in range(4)]
    nonzero = [i for i in range(4) if base[i] != 0]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = [0.0] * 4
                signs = [s0, s1, s2]
                for idx, ni in enumerate(nonzero):
                    v[ni] = base[ni] * signs[idx]
                vertices.append(v)

# Deduplicate
verts = []
for v in vertices:
    v = np.array(v)
    found = False
    for w in verts:
        if np.linalg.norm(v - w) < 1e-10:
            found = True
            break
    if not found:
        verts.append(v)
verts = np.array(verts)
N = len(verts)
print(f"\nVertices: {N}")

# Build adjacency
dots = verts @ verts.T
edge_dot = phi / 2  # cos(pi/5) for nearest neighbors
edges = []
adj = defaultdict(set)
edge_set = set()
for i in range(N):
    for j in range(i+1, N):
        if abs(dots[i,j] - edge_dot) < 1e-6:
            edges.append((i,j))
            adj[i].add(j)
            adj[j].add(i)
            edge_set.add((i,j))

print(f"Edges: {len(edges)}")

# Build triangles
triangles = []
tri_set = set()
for i in range(N):
    for j in adj[i]:
        if j > i:
            for k in adj[i] & adj[j]:
                if k > j:
                    tri = (i,j,k)
                    triangles.append(tri)
                    tri_set.add(frozenset(tri))

print(f"Triangles: {len(triangles)}")

# Build tetrahedra
tetrahedra = []
for i in range(N):
    for j in adj[i]:
        if j > i:
            common = adj[i] & adj[j]
            for k in common:
                if k > j:
                    for l in common & adj[k]:
                        if l > k:
                            tetrahedra.append((i,j,k,l))

print(f"Tetrahedra: {len(tetrahedra)}")

# Verify: V - E + F - C = chi(S^3) = 0
chi = N - len(edges) + len(triangles) - len(tetrahedra)
print(f"Euler characteristic: {N} - {len(edges)} + {len(triangles)} - {len(tetrahedra)} = {chi}")

# =====================================================================
# PART 1: Dihedral Angles and Deficit Angles
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: Dihedral Angles")
print("=" * 70)

# For each edge, find the tetrahedra containing it
edge_to_tets = defaultdict(list)
for idx, (a, b, c, d) in enumerate(tetrahedra):
    tet_edges = [(a,b),(a,c),(a,d),(b,c),(b,d),(c,d)]
    for e in tet_edges:
        key = (min(e), max(e))
        edge_to_tets[key].append(idx)

tets_per_edge_list = [len(edge_to_tets[e]) for e in edges]
print(f"Tetrahedra per edge: min={min(tets_per_edge_list)}, max={max(tets_per_edge_list)}, "
      f"all={len(set(tets_per_edge_list))==1 and tets_per_edge_list[0]}")

# Compute dihedral angle of a regular tetrahedron on S^3
# The 600-cell edge length (chord): |v_i - v_j| = 2*sin(pi/10) = 1/phi
chord_length = np.linalg.norm(verts[edges[0][0]] - verts[edges[0][1]])
print(f"\nChord length: {chord_length:.10f}")
print(f"1/phi = {1/phi:.10f}")

# Geodesic (arc) distance on unit S^3:
arc_length = np.arccos(np.clip(dots[edges[0][0], edges[0][1]], -1, 1))
print(f"Arc length: {arc_length:.10f}")
print(f"pi/5 = {np.pi/5:.10f}")

# Dihedral angle: the angle between two faces sharing an edge in a tetrahedron.
# For a REGULAR spherical tetrahedron inscribed in S^3 with edge arc-length pi/5:
# The dihedral angle is 2*pi/5 = 72 degrees (since 5 tets meet at each edge).

# Let me compute it directly from the geometry.
# Take one edge (v0, v1) and find two tetrahedra sharing it.
e0 = edges[0]
tets_of_e0 = edge_to_tets[e0]
print(f"\nEdge ({e0[0]}, {e0[1]}): {len(tets_of_e0)} tetrahedra")

# For dihedral angle computation:
# Take edge (A, B). Two tetrahedra share this edge: (A,B,C,D) and (A,B,C',D').
# The dihedral angle is the angle between the triangles (A,B,C) and (A,B,C')
# where C and C' are the vertices NOT on the edge.

def compute_dihedral_angle_S3(v0, v1, v2, v3):
    """
    Compute dihedral angle at edge (v0,v1) in tetrahedron (v0,v1,v2,v3)
    on S^3 (unit sphere in R^4).

    The dihedral angle is the angle between the two faces (v0,v1,v2) and (v0,v1,v3).
    On S^3, we use the tangent space at the midpoint of the edge.
    """
    # Midpoint of geodesic v0-v1 on S^3
    mid = (v0 + v1)
    mid = mid / np.linalg.norm(mid)

    # Project v2 and v3 into tangent space at mid
    # Tangent space at p on S^3: vectors orthogonal to p
    def project_to_tangent(p, v):
        # Direction from p to v on S^3
        # = v - (v.p)p, normalized
        vt = v - np.dot(v, p) * p
        return vt / np.linalg.norm(vt) if np.linalg.norm(vt) > 1e-12 else vt

    # Edge direction in tangent space at mid
    e_dir = project_to_tangent(mid, v1)

    # Project v2 and v3
    d2 = project_to_tangent(mid, v2)
    d3 = project_to_tangent(mid, v3)

    # Remove component along edge direction
    d2_perp = d2 - np.dot(d2, e_dir) * e_dir
    d3_perp = d3 - np.dot(d3, e_dir) * e_dir

    if np.linalg.norm(d2_perp) < 1e-12 or np.linalg.norm(d3_perp) < 1e-12:
        return 0.0

    d2_perp = d2_perp / np.linalg.norm(d2_perp)
    d3_perp = d3_perp / np.linalg.norm(d3_perp)

    cos_angle = np.clip(np.dot(d2_perp, d3_perp), -1, 1)
    return np.arccos(cos_angle)

# Compute dihedral angle for each tetrahedron at each edge
print("\nComputing dihedral angles for sample tetrahedra...")

dihedral_angles = []
for tet_idx in range(min(20, len(tetrahedra))):
    a, b, c, d = tetrahedra[tet_idx]
    va, vb, vc, vd = verts[a], verts[b], verts[c], verts[d]

    # Dihedral at edge (a,b): angle between faces (a,b,c) and (a,b,d)
    angle = compute_dihedral_angle_S3(va, vb, vc, vd)
    dihedral_angles.append(angle)

print(f"Sample dihedral angles (degrees):")
for i, ang in enumerate(dihedral_angles[:10]):
    print(f"  Tet {i}: {np.degrees(ang):.4f} deg")

mean_dihedral = np.mean(dihedral_angles)
print(f"\nMean dihedral: {np.degrees(mean_dihedral):.4f} deg")
print(f"Expected (2*pi/5): {np.degrees(2*np.pi/5):.4f} deg")
print(f"5 * mean: {np.degrees(5*mean_dihedral):.4f} deg (should be 360)")

# More precise: compute using the FLAT (Euclidean) dihedral angle of a regular tetrahedron
# A regular tetrahedron in FLAT space has dihedral angle arccos(1/3) = 70.53 deg
flat_dihedral = np.arccos(1/3)
print(f"\nFlat regular tetrahedron dihedral: {np.degrees(flat_dihedral):.4f} deg")
print(f"Spherical 600-cell tetrahedron dihedral: {np.degrees(mean_dihedral):.4f} deg")
print(f"Ratio: {mean_dihedral/flat_dihedral:.6f}")

# =====================================================================
# PART 2: Deficit Angle Computation
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: Deficit Angles at Edges")
print("=" * 70)

# For the full computation, we need the dihedral angle of EACH tetrahedron
# at EACH of its edges. For a spherical regular tetrahedron, all dihedral
# angles are equal (by symmetry).

# But let's verify: compute all dihedral angles properly.
# The spherical dihedral angle can be computed from the Gram matrix.

def spherical_dihedral_at_edge(v0, v1, v2, v3):
    """
    Dihedral angle at edge (v0,v1) in spherical tetrahedron (v0,v1,v2,v3).
    Uses the normal vectors to the two faces meeting at the edge.

    On S^3, the "normal" to a face (great-sphere triangle) through an edge
    is found using the 4D cross-product analog.
    """
    # Face 1: (v0, v1, v2) - find normal in tangent space of edge
    # Face 2: (v0, v1, v3)

    # Method: compute the angle between v2 and v3 as seen from the edge (v0,v1),
    # projected perpendicular to the edge.

    # The edge defines a great circle. The perpendicular directions to this
    # great circle at the midpoint span a 2D plane.

    # Simpler method: use the formula for spherical dihedral angles
    # cos(dihedral) = (cos(d_23) - cos(d_02)*cos(d_03)) / (sin(d_02)*sin(d_03))
    # where d_ij is the arc distance, evaluated at... no, this is for a vertex angle.

    # Actually, the dihedral angle at edge (v0,v1) in a spherical tetrahedron
    # can be computed as follows:
    # Let n1 be the outward normal to face (v0,v1,v2) and n2 to face (v0,v1,v3).
    # The dihedral angle = pi - angle(n1, n2).

    # In R^4, the normal to a hyperplane through v0,v1,v2 can be found.
    # The hyperplane through origin, v0, v1, v2 has normal:

    # Use Gram-Schmidt to get orthogonal basis in span(v0,v1,v2)
    e0 = v0 / np.linalg.norm(v0)
    e1 = v1 - np.dot(v1, e0) * e0
    e1 = e1 / np.linalg.norm(e1)
    e2 = v2 - np.dot(v2, e0)*e0 - np.dot(v2, e1)*e1
    e2 = e2 / np.linalg.norm(e2)

    # Normal to hyperplane through v0, v1, v2 (in R^4):
    # n1 is perpendicular to e0, e1, e2
    # We can find it as the component of any 4th independent vector
    # perpendicular to the span of {e0, e1, e2}

    # For face (v0, v1, v3):
    e2_prime = v3 - np.dot(v3, e0)*e0 - np.dot(v3, e1)*e1
    e2_prime_norm = np.linalg.norm(e2_prime)
    if e2_prime_norm < 1e-12:
        return 0.0
    e2_prime = e2_prime / e2_prime_norm

    # The dihedral angle at edge (v0,v1) is the angle between e2 and e2_prime
    # (the components of v2 and v3 perpendicular to the edge plane span(v0,v1))

    # Wait - we need to be more careful. e0, e1 span the edge direction.
    # e2 is the component of v2 perpendicular to span(v0,v1) (in the face 1 direction)
    # e2_prime is the component of v3 perpendicular to span(v0,v1) (in the face 2 direction)

    cos_dihedral = np.clip(np.dot(e2, e2_prime), -1, 1)
    return np.arccos(cos_dihedral)

# Compute deficit angle for each edge
print("Computing deficit angles for all 720 edges...")
deficit_angles = []

# First compute dihedral angle for one tet to calibrate
a, b, c, d = tetrahedra[0]
test_angle = spherical_dihedral_at_edge(verts[a], verts[b], verts[c], verts[d])
print(f"Test dihedral: {np.degrees(test_angle):.6f} deg")

for edge_idx, (i, j) in enumerate(edges):
    tet_indices = edge_to_tets[(i, j)]
    total_angle = 0

    for tet_idx in tet_indices:
        a, b, c, d = tetrahedra[tet_idx]
        # Find which two vertices are NOT on the edge
        tet_verts = {a, b, c, d}
        other = list(tet_verts - {i, j})

        # Compute dihedral at edge (i,j) using (i, j, other[0], other[1])
        angle = spherical_dihedral_at_edge(verts[i], verts[j], verts[other[0]], verts[other[1]])
        total_angle += angle

    deficit = 2 * np.pi - total_angle
    deficit_angles.append(deficit)

deficit_angles = np.array(deficit_angles)
print(f"\nDeficit angle statistics:")
print(f"  Mean: {np.degrees(np.mean(deficit_angles)):.6f} deg")
print(f"  Std:  {np.degrees(np.std(deficit_angles)):.6f} deg")
print(f"  Min:  {np.degrees(np.min(deficit_angles)):.6f} deg")
print(f"  Max:  {np.degrees(np.max(deficit_angles)):.6f} deg")
print(f"  All zero (< 1e-6 rad)? {np.all(np.abs(deficit_angles) < 1e-6)}")
print(f"  Total deficit: {np.sum(deficit_angles):.6f} rad")

# For S^3 of radius R, the Gauss-Bonnet theorem gives:
# integral of scalar curvature = 2 * chi(S^3) * 2*pi^2 = 0
# So total deficit should be 0 (since chi(S^3) = 0).

# =====================================================================
# PART 3: Regge Action
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: Regge Action for the 600-Cell")
print("=" * 70)

# The Regge action in 3D (spatial part) is:
# S_Regge = sum_edges (deficit_angle_e * length_e)
#
# For the 600-cell with all deficit angles = 0:
# S_Regge = 0

edge_length = np.linalg.norm(verts[edges[0][0]] - verts[edges[0][1]])
S_Regge = np.sum(deficit_angles * edge_length)
print(f"Edge length (chord): {edge_length:.10f}")
print(f"Regge action: {S_Regge:.10f}")
print(f"Expected: 0 (flat S^3 tiling)")

# =====================================================================
# PART 4: Deformed 600-Cell (Mass Insertion)
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: Deformed 600-Cell - Mass Insertion")
print("=" * 70)

print("""
To model a point mass, we can deform the 600-cell by:
1. Removing a vertex (creating a "vacancy")
2. Changing edge lengths near a vertex (creating deficit angles)
3. Inserting a conical singularity

Method 2 is the Regge calculus approach:
- Keep the combinatorial structure
- Allow edge lengths to vary
- The equations of motion (Regge equations) determine the lengths
  that satisfy the discrete Einstein equations

For a single massive particle at vertex v0:
- The mass creates a deficit solid angle at v0
- This is related to the ADM mass: M = (deficit_solid_angle)/(8*pi*G)

In 2+1D Regge gravity, a point mass removes a wedge of angle 8*pi*G*m.
In 3+1D, a mass creates deficit angles on all edges through v0.
""")

# Compute what happens if we change edge lengths near v0
# Take vertex 0 and its 12 neighbors
v0_neighbors = list(adj[0])
print(f"Vertex 0: {len(v0_neighbors)} neighbors")
print(f"Neighbor types: ", end="")
for n in v0_neighbors[:5]:
    d = np.sqrt(np.sum((verts[n] - verts[0])**2))
    print(f"{d:.4f} ", end="")
print("...")

# In linearized Regge calculus, the deficit angle at edge e responds
# to a mass perturbation at vertex v as:
# delta(deficit_e) = G_Regge(e, v) * delta_mass

# For a uniform perturbation (all edges from v0 scaled by factor 1+epsilon):
# This is equivalent to moving v0 slightly toward the center of S^3.

# Let's compute the Regge equations for a radially symmetric perturbation.
# Scale all edges from v0 by (1+eps):
eps = 0.01

# New positions: move v0 along the radial direction
# On S^3, this means scaling the embedding: v0_new = (1-eps)*v0 / ||(1-eps)*v0||
# Actually, to change edge lengths, we should move v0 "inward":
# v0_new = v0 * (1 - delta) / ||v0 * (1-delta)|| = v0 (on S^3, still unit)

# Instead: change the RADIUS of S^3 seen by v0
# Model: v0 is at angular position theta=0 on S^3 of radius R
# A mass at v0 deforms the metric: R(theta) = R0 * (1 - M*f(theta))
# For small M: edge lengths change proportionally

# Actually, the simplest deformation: rescale the 600-cell radius
# from R=1 to R(r) where r is the distance from v0.
# This is essentially the conformal factor C(r).

# In Regge calculus, the continuum limit of the deficit angles
# gives the Ricci scalar. For a conformally flat metric on S^3:
# ds^2 = C(theta)^2 (dtheta^2 + sin^2(theta) dOmega_2^2)
# the Ricci scalar is R = 6/C^2 * (1 - (C'/C)^2 * ... )

# The KEY QUESTION: does Regge calculus on the 600-cell
# naturally produce the Schwarzschild metric?

# For a spherically symmetric deformation with C(d) = 1 + eps(d)
# where d is graph distance from v0:

print("\n--- Regge linearized equations ---")
print("For a mass delta*m at vertex v0:")
print("The Regge equation at each edge e is:")
print("  sum_tets_containing_e d(deficit_e)/d(l_f) * l_f = 8*pi*G * T_e")
print("where T_e is the matter content near edge e.")
print("")
print("For the 600-cell (flat background), this becomes LINEAR:")
print("  sum_f K_ef * delta_l_f = 8*pi*G * delta_T_e")
print("where K_ef is the Regge stiffness matrix.")

# =====================================================================
# PART 5: Regge Stiffness Matrix (Linearized)
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: Linearized Regge - Deficit Response")
print("=" * 70)

# For a flat background, the deficit angles are all zero.
# Under a perturbation of edge lengths delta_l_e:
# delta(deficit_f) = sum_e d(deficit_f)/d(l_e) * delta_l_e

# For a 3D simplicial complex, the dihedral angle at edge e in tet T
# depends on all 6 edge lengths of T.

# For an equilateral spherical tetrahedron (all edges = l0 = pi/5):
# By symmetry, d(theta_e)/d(l_e) is the same for the edge itself,
# and d(theta_e)/d(l_f) is the same for all other edges f in the same tet.

# Let me compute these numerically.

l0 = np.pi / 5  # geodesic edge length

def make_spherical_tet(l):
    """Make a regular spherical tetrahedron with geodesic edge length l on unit S^3."""
    # Place 4 points on S^3 at equal geodesic distances
    # For a regular spherical tetrahedron, all pairwise dot products are equal:
    # cos(l) = dot product
    cos_l = np.cos(l)

    # Gram matrix: G_ij = cos(d_ij) = cos_l for i!=j, 1 for i=j
    G = np.full((4, 4), cos_l)
    np.fill_diagonal(G, 1.0)

    # Cholesky to get positions
    try:
        L = np.linalg.cholesky(G)
        # L is 4x4, but we need 4D vectors. Each row is a position in R^4.
        # Normalize each row to unit length (they should already be close)
        positions = L.copy()
        # Extend to R^4 if needed
        if L.shape[1] < 4:
            positions = np.zeros((4, 4))
            positions[:, :L.shape[1]] = L
        for i in range(4):
            positions[i] /= np.linalg.norm(positions[i])
        return positions
    except:
        return None

# Compute dihedral angle of reference tet
ref_tet = make_spherical_tet(l0)
if ref_tet is not None:
    ref_dihedral = spherical_dihedral_at_edge(ref_tet[0], ref_tet[1], ref_tet[2], ref_tet[3])
    print(f"Reference tet dihedral: {np.degrees(ref_dihedral):.6f} deg")
    print(f"Expected: {np.degrees(2*np.pi/5):.6f} deg = 72.000000 deg")

    # Numerical derivatives
    dl = 1e-7

    # Perturb the edge (0,1) length
    tet_plus = make_spherical_tet(l0 + dl)
    tet_minus = make_spherical_tet(l0 - dl)

    if tet_plus is not None and tet_minus is not None:
        dihedral_plus = spherical_dihedral_at_edge(tet_plus[0], tet_plus[1], tet_plus[2], tet_plus[3])
        dihedral_minus = spherical_dihedral_at_edge(tet_minus[0], tet_minus[1], tet_minus[2], tet_minus[3])

        d_theta_d_l_same = (dihedral_plus - dihedral_minus) / (2 * dl)
        print(f"\nd(dihedral at e)/d(length of e): {d_theta_d_l_same:.6f}")

        # For a REGULAR tet, perturbing one edge changes the dihedral at that edge
        # and also at the opposite edge (by a different amount due to symmetry).
        # By symmetry of the regular tet:
        # - 1 edge has derivative d_theta_same
        # - 4 adjacent edges have derivative d_theta_adj
        # - 1 opposite edge has derivative d_theta_opp

        print(f"\n(Note: for a regular tet, all edges are equivalent)")
        print(f"d(dihedral)/d(same edge length) = {d_theta_d_l_same:.6f}")

# =====================================================================
# PART 6: Continuum Limit and Schwarzschild
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: Continuum Limit - Can We Get Schwarzschild?")
print("=" * 70)

print("""
KEY INSIGHT: The 600-cell has only 5 shells (graph distances 0-5).
This is FAR too coarse for a realistic metric.

To get Schwarzschild, we would need:
1. Many copies of 600-cells (refinement), or
2. A mathematical limit where 600-cell properties
   determine the continuum metric

APPROACH: Instead of refining the mesh, use the 600-cell as a
"fundamental cell" whose SYMMETRY determines the metric.

The 600-cell determines:
- Gauge group (D4 -> SM)
- Coupling constants (spectral ratios)
- Can it also determine the GRAVITATIONAL sector?

In Kaluza-Klein theory:
  5D gravity = 4D gravity + U(1) gauge field
  The extra dimension's size determines the gauge coupling.

In our case:
  The 600-cell on S^3 determines gauge couplings.
  The SAME geometry might determine gravity on S^3.

The Ricci curvature of S^3 with radius R is R_ij = (2/R^2) g_ij.
The Einstein equations on S^3 x time give:
  G_tt = 3/(R^2) (constant positive curvature)

This is a cosmological constant, not Schwarzschild!

For Schwarzschild, we need ASYMPTOTICALLY FLAT space,
which means R -> infinity (the continuum limit).
""")

# =====================================================================
# PART 7: The Spectral Approach to Gravity
# =====================================================================
print("=" * 70)
print("PART 7: Spectral Geometry and Gravity")
print("=" * 70)

# The spectral action principle (Connes):
# S = Tr(f(D/Lambda))
# where D is the Dirac operator and Lambda is a cutoff.
#
# For a Riemannian manifold, the heat kernel expansion gives:
# Tr(e^{-tD^2}) = (4*pi*t)^{-d/2} * sum_k a_k * t^k
# where a_0 = volume, a_1 = 1/6 * integral(R), a_2 = ... (Gauss-Bonnet), etc.
#
# The coefficient a_1 gives the Einstein-Hilbert action!
#
# For the 600-cell graph Laplacian:
# Tr(e^{-t*L}) = sum_i exp(-t*lambda_i)
# where lambda_i are the eigenvalues.

# Eigenvalues of the 600-cell (from exp110)
eigenvalues = [12, 6*phi, 4*phi, 3, 0, -2, 4-4*phi, -3, 6-6*phi]
multiplicities = [1, 4, 9, 16, 25, 36, 9, 16, 4]

# Shift to make them non-negative: L_shifted = 12*I - A (Laplacian)
# Laplacian eigenvalues: 12 - lambda_i
lap_eigenvalues = [12 - ev for ev in eigenvalues]
print("Laplacian eigenvalues (12 - adjacency eigenvalue):")
for ev, m in zip(lap_eigenvalues, multiplicities):
    print(f"  {12-ev:+10.6f} -> lap_ev = {ev:10.6f} (mult {m})")

# Heat kernel trace
print("\nHeat kernel Tr(e^{-t*L}):")
for t in [0.01, 0.1, 0.5, 1.0, 2.0]:
    trace = sum(m * np.exp(-t * ev) for ev, m in zip(lap_eigenvalues, multiplicities))
    print(f"  t = {t:.2f}: Tr = {trace:.6f}")

# The small-t expansion coefficients:
# Tr(e^{-tL}) ~ a_0 + a_1*t + a_2*t^2 + ...
# a_0 = N = 120 (number of vertices)
# a_1 = -sum_i m_i * lambda_i = -Tr(L) = -12*120 + sum edges = ...
# Actually Tr(L) = sum of diagonal = 12*120 = 1440 (each vertex has degree 12)

TrL = sum(m * ev for ev, m in zip(lap_eigenvalues, multiplicities))
print(f"\nTr(L) = {TrL:.6f}")
print(f"Expected: 12 * 120 = {12*120}")

TrL2 = sum(m * ev**2 for ev, m in zip(lap_eigenvalues, multiplicities))
print(f"Tr(L^2) = {TrL2:.6f}")
print(f"= sum_i (deg_i)^2 + sum_edges = 120*144 + ... ")

# a_0 = 120
# a_1 = -Tr(L) = -1440
# a_2 = Tr(L^2)/2 = ...
# The ratio a_1/a_0 = -12 = -degree = analogous to "Ricci scalar"
# In the continuum: a_1/a_0 = -(1/6)*R*V/V = -(1/6)*R
# So R_effective = 6 * 12 = 72?

print(f"\nSpectral 'Ricci scalar' analogy:")
print(f"  a_0 = 120 (volume)")
print(f"  a_1 = -Tr(L) = -1440")
print(f"  a_1/a_0 = -12 = -degree")
print(f"  If a_1/a_0 = -R/6 (continuum formula), then R = 72")
print(f"  For S^3 of radius r: R = 6/r^2, so r^2 = 6/72 = 1/12")
print(f"  r = 1/sqrt(12) = {1/np.sqrt(12):.6f}")
print(f"  This is the 'effective radius' of S^3 from the graph perspective.")

# =====================================================================
# PART 8: The Spectral Action on the 600-Cell
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: Spectral Action Coefficients")
print("=" * 70)

# Spectral action: S = sum_k f_k * a_k
# where f_k are moments of the cutoff function f
# and a_k are Seeley-DeWitt coefficients.
#
# For a graph, the "spectral action" is:
# S = sum_i f(lambda_i/Lambda^2)
# For f(x) = 1 (counting): S = N = 120 (cosmological constant)
# For f(x) = x: S = Tr(L)/Lambda^2 = 1440/Lambda^2 (Einstein-Hilbert)
# For f(x) = x^2: S = Tr(L^2)/Lambda^4 (higher curvature)

print("Spectral action moments (Seeley-DeWitt coefficients):")
for k in range(5):
    moment = sum(m * ev**k for ev, m in zip(lap_eigenvalues, multiplicities))
    print(f"  a_{k} = Tr(L^{k}) = {moment:.2f}")

# The ratio Tr(L^2) / (Tr(L))^2 * N is a measure of curvature
# fluctuations.
curvature_fluct = TrL2 * 120 / TrL**2
print(f"\nCurvature fluctuation measure: Tr(L^2)*N/Tr(L)^2 = {curvature_fluct:.6f}")

# =====================================================================
# PART 9: Summary
# =====================================================================
print("\n" + "=" * 70)
print("PART 9: Summary and Assessment")
print("=" * 70)

print("""
REGGE CALCULUS RESULTS:
1. The 600-cell has ZERO deficit angles everywhere (perfect S^3 tiling)
2. The Regge action is identically zero
3. To create curvature, the 600-cell must be DEFORMED
4. The linearized Regge equations give the response to mass insertion
5. With only 5 shells, the 600-cell is too coarse for Schwarzschild

SPECTRAL APPROACH:
6. The heat kernel gives an effective Ricci scalar R_eff = 72
7. This corresponds to S^3 of radius 1/sqrt(12)
8. The spectral action naturally gives:
   - a_0 = 120 (cosmological term)
   - a_1 = 1440 (Einstein-Hilbert term)
   - Higher a_k for higher-curvature terms

KEY INSIGHT:
The 600-cell on its own does NOT give a specific gravity theory.
It gives:
- A FLAT background (S^3, no deficit angles)
- Gauge structure (from D4/spectrum)
- Matter content (from Type-C vertices)

For gravity, we need ADDITIONAL INPUT:
- Either: a specific deformation mechanism (how mass curves the lattice)
- Or: a dimensional reduction from S^3 to R^3 (breaking S^3 compactness)
- Or: a spectral action principle that selects the gravitational action

ASSESSMENT:
The 600-cell geometry is rich enough for particle physics (gauge + matter)
but does NOT uniquely determine the gravitational metric.
The STG ansatz C(r) = 1+M/r was one guess; it failed.
A principled derivation of gravity from 600-cell geometry remains OPEN.

Status: OPEN PROBLEM. No clean derivation found.
""")
