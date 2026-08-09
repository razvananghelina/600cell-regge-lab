"""
Vizualizare 600-cell: Spatiul Teoriei
=====================================
Proiectie stereografica din 4D in 3D
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

PHI = (1 + np.sqrt(5)) / 2

def generate_600cell_vertices():
    """
    Genereaza cele 120 de varfuri ale 600-cell.
    Coordonate 4D normalizate.
    """
    vertices = []

    # 1. 8 varfuri: permutatii de (±1, 0, 0, 0)
    for i in range(4):
        for sign in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = sign
            vertices.append(v)

    # 2. 16 varfuri: (±0.5, ±0.5, ±0.5, ±0.5)
    for s1 in [0.5, -0.5]:
        for s2 in [0.5, -0.5]:
            for s3 in [0.5, -0.5]:
                for s4 in [0.5, -0.5]:
                    vertices.append([s1, s2, s3, s4])

    # 3. 96 varfuri: permutatii pare de (±phi/2, ±0.5, ±1/(2*phi), 0)
    coords = [PHI/2, 0.5, 1/(2*PHI), 0]

    # Toate permutatiile
    from itertools import permutations
    for perm in permutations(range(4)):
        base = [coords[perm[0]], coords[perm[1]], coords[perm[2]], coords[perm[3]]]
        # Toate combinatiile de semne (pentru coordonatele nenule)
        for s0 in ([1, -1] if base[0] != 0 else [1]):
            for s1 in ([1, -1] if base[1] != 0 else [1]):
                for s2 in ([1, -1] if base[2] != 0 else [1]):
                    for s3 in ([1, -1] if base[3] != 0 else [1]):
                        v = [s0*base[0], s1*base[1], s2*base[2], s3*base[3]]
                        # Verificam paritatea (numar par de semne negative)
                        neg_count = sum(1 for x in [s0, s1, s2, s3] if x == -1)
                        if neg_count % 2 == 0:
                            vertices.append(v)

    # Eliminam duplicate si normalizam
    unique = []
    for v in vertices:
        v_norm = np.array(v)
        norm = np.linalg.norm(v_norm)
        if norm > 0:
            v_norm = v_norm / norm
        is_dup = False
        for u in unique:
            if np.allclose(v_norm, u, atol=1e-6):
                is_dup = True
                break
        if not is_dup:
            unique.append(v_norm)

    return np.array(unique)

def stereographic_projection(v4d, pole_index=0):
    """
    Proiectie stereografica din 4D in 3D.
    Proiecteaza de la polul nord (0,0,0,1).
    """
    w = v4d[3]  # Coordonata a 4-a
    if abs(w - 1) < 1e-6:
        return None  # La pol
    scale = 1 / (1 - w)
    return np.array([v4d[0] * scale, v4d[1] * scale, v4d[2] * scale])

def find_edges(vertices, edge_length_sq=None):
    """
    Gaseste muchiile (varfuri la distanta minima).
    """
    n = len(vertices)

    # Calculam toate distantele
    distances = []
    for i in range(n):
        for j in range(i+1, n):
            d = np.sum((vertices[i] - vertices[j])**2)
            distances.append((d, i, j))

    # Sortam si gasim lungimea minima a muchiei
    distances.sort()
    if edge_length_sq is None:
        edge_length_sq = distances[0][0]

    # Selectam muchiile cu aceasta lungime (cu toleranta)
    edges = []
    for d, i, j in distances:
        if d < edge_length_sq * 1.1:  # 10% toleranta
            edges.append((i, j))
        else:
            break

    return edges

def create_icosahedron():
    """Creeaza un icosaedru pentru referinta."""
    verts = []
    # 12 varfuri
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            verts.append([0, s1, s2*PHI])
            verts.append([s1, s2*PHI, 0])
            verts.append([s2*PHI, 0, s1])

    verts = np.array(verts)
    # Normalizam
    verts = verts / np.linalg.norm(verts[0])
    return verts

def main():
    print("Generez varfurile 600-cell...")
    vertices_4d = generate_600cell_vertices()
    print(f"Numar de varfuri 4D: {len(vertices_4d)}")

    # Proiectam in 3D
    print("Proiectez in 3D...")
    vertices_3d = []
    valid_indices = []
    for i, v in enumerate(vertices_4d):
        p = stereographic_projection(v)
        if p is not None and np.linalg.norm(p) < 10:  # Limitam raza
            vertices_3d.append(p)
            valid_indices.append(i)

    vertices_3d = np.array(vertices_3d)
    print(f"Varfuri proiectate: {len(vertices_3d)}")

    # Gasim muchiile in 4D
    print("Calculez muchiile...")
    edges_4d = find_edges(vertices_4d)
    print(f"Muchii gasite: {len(edges_4d)}")

    # Mapam muchiile la indicii valizi
    index_map = {old: new for new, old in enumerate(valid_indices)}
    edges_3d = []
    for i, j in edges_4d:
        if i in index_map and j in index_map:
            edges_3d.append((index_map[i], index_map[j]))

    print(f"Muchii vizibile: {len(edges_3d)}")

    # ===============================
    # VIZUALIZARE
    # ===============================
    fig = plt.figure(figsize=(16, 12))

    # ---- Plot 1: 600-cell complet ----
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.set_title('600-cell: Proiectie Stereografica', fontsize=12, fontweight='bold')

    # Desenam muchiile
    for i, j in edges_3d:
        p1, p2 = vertices_3d[i], vertices_3d[j]
        dist = np.linalg.norm(p1)  # Distanta de la centru
        alpha = max(0.1, 1 - dist/5)  # Mai transparente cele departate
        ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                'b-', alpha=alpha, linewidth=0.3)

    # Varfuri
    ax1.scatter(vertices_3d[:, 0], vertices_3d[:, 1], vertices_3d[:, 2],
               c='gold', s=10, alpha=0.7, edgecolors='orange')

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_box_aspect([1,1,1])

    # ---- Plot 2: Straturile concentrice ----
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    ax2.set_title('Straturi Concentrice (distanta de la centru)', fontsize=12, fontweight='bold')

    # Coloram varfurile dupa distanta de la centru
    distances = np.linalg.norm(vertices_3d, axis=1)
    scatter = ax2.scatter(vertices_3d[:, 0], vertices_3d[:, 1], vertices_3d[:, 2],
                         c=distances, cmap='plasma', s=30, alpha=0.8)
    plt.colorbar(scatter, ax=ax2, label='Distanta', shrink=0.5)

    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_box_aspect([1,1,1])

    # ---- Plot 3: Structura icosaedrica centrala ----
    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    ax3.set_title('Nucleul Icosaedric (20 tetraedre)', fontsize=12, fontweight='bold')

    # Selectam varfurile cele mai apropiate de centru
    central_mask = distances < np.percentile(distances, 15)
    central_verts = vertices_3d[central_mask]

    # Desenam conexiunile centrale
    for i, j in edges_3d:
        if central_mask[i] and central_mask[j]:
            p1, p2 = vertices_3d[i], vertices_3d[j]
            ax3.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                    'r-', linewidth=1.5, alpha=0.8)

    ax3.scatter(central_verts[:, 0], central_verts[:, 1], central_verts[:, 2],
               c='red', s=80, alpha=1, edgecolors='darkred', linewidths=2)

    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    ax3.set_box_aspect([1,1,1])

    # ---- Plot 4: Geodezice (cercuri mari) ----
    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    ax4.set_title('Geodezice pe S^3 (Fibre Hopf)', fontsize=12, fontweight='bold')

    # Desenam cateva cercuri mari
    theta = np.linspace(0, 2*np.pi, 100)

    colors = ['red', 'green', 'blue', 'orange', 'purple', 'cyan']
    for k, color in enumerate(colors):
        # Cercuri in diferite planuri
        angle = k * np.pi / 6
        x = np.cos(theta) * np.cos(angle)
        y = np.cos(theta) * np.sin(angle)
        z = np.sin(theta)
        ax4.plot(x, y, z, color=color, linewidth=2, alpha=0.7,
                label=f'Fibra {k+1}')

    # Adaugam varfurile ca noduri
    # Normalizam pe sfera unitate
    verts_unit = vertices_3d / (np.linalg.norm(vertices_3d, axis=1, keepdims=True) + 1e-6)
    verts_unit = verts_unit[distances < 3]  # Doar cele apropiate

    ax4.scatter(verts_unit[:, 0], verts_unit[:, 1], verts_unit[:, 2],
               c='gold', s=15, alpha=0.5)

    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z')
    ax4.set_box_aspect([1,1,1])
    ax4.legend(loc='upper left', fontsize=8)

    plt.tight_layout()
    plt.savefig('600cell_visualization.png', dpi=150, bbox_inches='tight')
    print("\nSalvat: 600cell_visualization.png")
    plt.show()

    # ===============================
    # FIGURA 2: DETALII
    # ===============================
    fig2, axes = plt.subplots(1, 3, figsize=(18, 6))

    # ---- Sectiune 2D prin centru ----
    ax = axes[0]
    ax.set_title('Sectiune 2D (plan XY)', fontsize=12, fontweight='bold')

    # Selectam varfurile aproape de planul z=0
    z_slice = np.abs(vertices_3d[:, 2]) < 0.5
    slice_verts = vertices_3d[z_slice]

    ax.scatter(slice_verts[:, 0], slice_verts[:, 1], c='blue', s=50, alpha=0.7)

    # Desenam conexiunile din aceasta sectiune
    for i, j in edges_3d:
        if z_slice[i] and z_slice[j]:
            p1, p2 = vertices_3d[i], vertices_3d[j]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'b-', alpha=0.3, linewidth=0.5)

    # Desenam un decagon (simetrie 10-fold)
    angles = np.linspace(0, 2*np.pi, 11)
    r = 1.5
    ax.plot(r*np.cos(angles), r*np.sin(angles), 'r--', linewidth=2, label='Decagon')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.axis('equal')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # ---- Histograma distantelor ----
    ax = axes[1]
    ax.set_title('Distributia Distantelor (straturi)', fontsize=12, fontweight='bold')

    ax.hist(distances, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(x=1, color='red', linestyle='--', label='r=1 (inner shell)')
    ax.axvline(x=PHI, color='green', linestyle='--', label=f'r=phi={PHI:.2f}')
    ax.axvline(x=PHI**2, color='orange', linestyle='--', label=f'r=phi^2={PHI**2:.2f}')

    ax.set_xlabel('Distanta de la centru')
    ax.set_ylabel('Numar de varfuri')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # ---- Schema conceptuala ----
    ax = axes[2]
    ax.set_title('Schema: Electronul ca Soliton', fontsize=12, fontweight='bold')

    # Desenam o grila hexagonala simplificata
    for i in range(-3, 4):
        for j in range(-3, 4):
            x = i + 0.5*j
            y = j * np.sqrt(3)/2
            circle = plt.Circle((x, y), 0.15, color='lightblue', alpha=0.5)
            ax.add_patch(circle)

    # Zona de densitate ridicata = "electronul"
    electron = plt.Circle((0, 0), 0.8, color='red', alpha=0.3, label='Electron (soliton)')
    ax.add_patch(electron)
    ax.scatter([0], [0], c='red', s=200, zorder=5)

    # Linii de camp
    for angle in np.linspace(0, 2*np.pi, 13)[:-1]:
        r = np.linspace(0.3, 2.5, 20)
        x = r * np.cos(angle)
        y = r * np.sin(angle)
        ax.plot(x, y, 'orange', alpha=0.4, linewidth=1)

    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig('600cell_details.png', dpi=150, bbox_inches='tight')
    print("Salvat: 600cell_details.png")
    plt.show()

if __name__ == "__main__":
    main()
