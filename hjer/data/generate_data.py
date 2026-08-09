"""
Generate JSON data for 600-cell and 120-cell polytopes.
Run once: python generate_data.py
Outputs: 600cell.json, 120cell.json
"""

import numpy as np
import json
from itertools import permutations
from collections import deque

PHI = (1 + np.sqrt(5)) / 2


def generate_600cell_vertices():
    """
    Generate the 120 vertices of the 600-cell on the unit 3-sphere.
    Uses the proven algorithm from exp028_strong_force_dual.py:
    - 16 vertices: all sign combinations of (+-0.5, +-0.5, +-0.5, +-0.5)
    - 8 vertices: permutations of (+-1, 0, 0, 0)
    - 96 vertices: even permutations of (phi/2, 0.5, 1/(2*phi), 0) with all signs
    Then deduplicate via rounding.
    """
    vertices = []

    # Family 2 first (16 vertices): (+-0.5, +-0.5, +-0.5, +-0.5)
    for i in range(16):
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
        vertices.append(signs * 0.5)

    # Family 1 (8 vertices): permutations of (+-1, 0, 0, 0)
    for i in range(4):
        v = np.zeros(4)
        v[i] = 1.0
        vertices.append(v)
        vertices.append(-v)

    # Family 3 (96 vertices): even permutations of (phi/2, 0.5, 1/(2*phi), 0)
    # with all 16 sign combinations
    v_base = np.array([PHI, 1.0, 1.0 / PHI, 0.0]) * 0.5

    # Get even permutations (even number of inversions)
    even_perms = []
    for p in permutations([0, 1, 2, 3]):
        inv_count = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv_count % 2 == 0:
            even_perms.append(p)

    for p in even_perms:
        perm_vals = v_base[list(p)]
        for i in range(16):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
            vertices.append(perm_vals * signs)

    vertices = np.array(vertices)
    # Deduplicate using rounding
    _, idx = np.unique(np.round(vertices, 8), axis=0, return_index=True)
    return vertices[idx]


def find_edges(vertices):
    """Find edges (pairs at minimum distance) and build adjacency."""
    n = len(vertices)
    # Compute all pairwise squared distances
    dists_sq = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sum((vertices[i] - vertices[j]) ** 2)
            dists_sq[i, j] = d
            dists_sq[j, i] = d

    # Find minimum edge length
    upper = dists_sq[np.triu_indices(n, k=1)]
    min_dist_sq = np.min(upper)

    # Select edges within 10% tolerance
    edges = []
    adjacency = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if dists_sq[i, j] < min_dist_sq * 1.1:
                edges.append([i, j])
                adjacency[i].append(j)
                adjacency[j].append(i)

    return edges, adjacency


def bfs_distances(adjacency, source):
    """BFS from source, return distance array."""
    n = len(adjacency)
    dist = [-1] * n
    dist[source] = 0
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in adjacency[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


def compute_graph_diameter(adjacency):
    """Compute graph diameter via BFS from all vertices."""
    n = len(adjacency)
    diameter = 0
    for i in range(n):
        d = bfs_distances(adjacency, i)
        diameter = max(diameter, max(d))
    return diameter


def find_tetrahedra(adjacency):
    """Find all tetrahedra (4-cliques) in the graph."""
    n = len(adjacency)
    neighbor_sets = [set(adjacency[i]) for i in range(n)]
    tetrahedra = []

    for i in range(n):
        for j in neighbor_sets[i]:
            if j <= i:
                continue
            common_ij = neighbor_sets[i].intersection(neighbor_sets[j])
            for k in common_ij:
                if k <= j:
                    continue
                common_ijk = common_ij.intersection(neighbor_sets[k])
                for l in common_ijk:
                    if l <= k:
                        continue
                    tetrahedra.append([i, j, k, l])

    return tetrahedra


def build_120cell(vertices_600, tetrahedra):
    """Build 120-cell as dual: vertices are tetrahedra centers."""
    centers = []
    for tet in tetrahedra:
        pts = vertices_600[tet]
        center = np.mean(pts, axis=0)
        centers.append(center)

    centers = np.array(centers)
    # Normalize to unit sphere
    norms = np.linalg.norm(centers, axis=1, keepdims=True)
    centers = centers / norms

    return centers


def save_polytope_json(filename, vertices4d, edges, adjacency, degree, diameter, name):
    """Save polytope data as JSON."""
    data = {
        "name": name,
        "vertices4D": vertices4d.tolist(),
        "edges": edges,
        "adjacency": adjacency,
        "numVertices": len(vertices4d),
        "numEdges": len(edges),
        "degree": degree,
        "graphDiameter": diameter,
    }
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  Saved {filename}: {len(vertices4d)} vertices, {len(edges)} edges, "
          f"degree {degree}, diameter {diameter}")


def main():
    print("=" * 60)
    print("Generating polytope data for HJER engine")
    print("=" * 60)

    # --- 600-cell ---
    print("\n[1] Generating 600-cell vertices...")
    verts_600 = generate_600cell_vertices()
    print(f"  Vertices: {len(verts_600)} (expected 120)")

    print("  Finding edges...")
    edges_600, adj_600 = find_edges(verts_600)
    degree_600 = len(adj_600[0])
    print(f"  Edges: {len(edges_600)} (expected 720), Degree: {degree_600} (expected 12)")

    print("  Computing graph diameter...")
    diam_600 = compute_graph_diameter(adj_600)
    print(f"  Diameter: {diam_600} (expected 5)")

    save_polytope_json("600cell.json", verts_600, edges_600, adj_600,
                       degree_600, diam_600, "600-cell")

    # --- 120-cell (dual) ---
    print("\n[2] Finding tetrahedra for 120-cell dual...")
    tetrahedra = find_tetrahedra(adj_600)
    print(f"  Tetrahedra found: {len(tetrahedra)} (expected 600)")

    print("  Building 120-cell vertices from tetrahedra centers...")
    verts_120 = build_120cell(verts_600, tetrahedra)
    print(f"  Vertices: {len(verts_120)} (expected 600)")

    print("  Finding edges...")
    edges_120, adj_120 = find_edges(verts_120)
    degree_120 = len(adj_120[0]) if len(adj_120) > 0 else 0
    print(f"  Edges: {len(edges_120)} (expected 1200), Degree: {degree_120} (expected 4)")

    print("  Computing graph diameter (this may take a moment)...")
    diam_120 = compute_graph_diameter(adj_120)
    print(f"  Diameter: {diam_120}")

    save_polytope_json("120cell.json", verts_120, edges_120, adj_120,
                       degree_120, diam_120, "120-cell")

    print("\nDone! JSON files ready for HJER engine.")


if __name__ == "__main__":
    main()
