"""
Vizualizare Tetraedre din 600-cell
==================================
Desenam tetraedrele 3D cu fete colorate si unghiuri corecte.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import combinations

PHI = (1 + np.sqrt(5)) / 2

def create_icosahedron_vertices():
    """
    Creeaza cele 12 varfuri ale unui icosaedru.
    Acestea formeaza structura centrala a 600-cell.
    """
    vertices = []

    # Varfurile icosaedrului in coordonate carteziene
    # Trei dreptunghiuri aurii perpendiculare
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            # Dreptunghi in planul YZ
            vertices.append([0, s1 * 1, s2 * PHI])
            # Dreptunghi in planul XZ
            vertices.append([s1 * 1, s2 * PHI, 0])
            # Dreptunghi in planul XY
            vertices.append([s1 * PHI, 0, s2 * 1])

    vertices = np.array(vertices)
    # Normalizam la raza 1
    r = np.sqrt(1 + PHI**2)
    vertices = vertices / r
    return vertices

def find_icosahedron_faces(vertices, tolerance=0.01):
    """
    Gaseste cele 20 de fete triunghiulare ale icosaedrului.
    """
    n = len(vertices)

    # Calculam toate distantele si gasim lungimea muchiei minime
    all_distances = []
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(vertices[i] - vertices[j])
            all_distances.append(d)

    all_distances.sort()
    edge_length = all_distances[0]  # Cea mai scurta distanta = muchia
    print(f"Lungime muchie icosaedru: {edge_length:.4f}")

    # Gasim toate combinatiile de 3 varfuri care formeaza triunghiuri echilaterale
    faces = []
    for combo in combinations(range(n), 3):
        i, j, k = combo
        d1 = np.linalg.norm(vertices[i] - vertices[j])
        d2 = np.linalg.norm(vertices[j] - vertices[k])
        d3 = np.linalg.norm(vertices[k] - vertices[i])

        # Toate laturile trebuie sa fie egale cu lungimea muchiei
        if (abs(d1 - edge_length) < tolerance and
            abs(d2 - edge_length) < tolerance and
            abs(d3 - edge_length) < tolerance):
            faces.append(combo)

    return faces

def create_tetrahedra_from_center(vertices, faces):
    """
    Creeaza 20 de tetraedre de la centru (0,0,0) la fiecare fata.
    Acestea reprezinta cele 20 de tetraedre in jurul unui vertex al 600-cell.
    """
    tetrahedra = []
    center = np.array([0, 0, 0])

    for face in faces:
        i, j, k = face
        tetra = [center, vertices[i], vertices[j], vertices[k]]
        tetrahedra.append(tetra)

    return tetrahedra

def get_tetrahedron_faces(tetra):
    """
    Returneaza cele 4 fete triunghiulare ale unui tetraedru.
    """
    return [
        [tetra[0], tetra[1], tetra[2]],  # Fata 1
        [tetra[0], tetra[1], tetra[3]],  # Fata 2
        [tetra[0], tetra[2], tetra[3]],  # Fata 3
        [tetra[1], tetra[2], tetra[3]],  # Fata 4 (baza)
    ]

def main():
    print("=" * 60)
    print("VIZUALIZARE: 20 TETRAEDRE DIN 600-CELL")
    print("=" * 60)

    # Cream icosaedrul
    vertices = create_icosahedron_vertices()
    print(f"Varfuri icosaedru: {len(vertices)}")

    # Gasim fetele
    faces = find_icosahedron_faces(vertices)
    print(f"Fete icosaedru: {len(faces)}")

    # Cream tetraedrele
    tetrahedra = create_tetrahedra_from_center(vertices, faces)
    print(f"Tetraedre: {len(tetrahedra)}")

    # Calculam unghiurile
    tetra = tetrahedra[0]
    v0, v1, v2, v3 = [np.array(v) for v in tetra]

    # Vectorii de la centru
    e1 = v1 - v0
    e2 = v2 - v0
    e3 = v3 - v0

    # Unghiul diedru
    angle12 = np.arccos(np.dot(e1, e2) / (np.linalg.norm(e1) * np.linalg.norm(e2)))
    print(f"\nUnghiul intre muchii: {np.degrees(angle12):.2f} grade")
    print(f"Unghiul teoretic (arccos(1/phi^2)): {np.degrees(np.arccos(1/PHI**2)):.2f} grade")

    # ===============================
    # VIZUALIZARE
    # ===============================
    fig = plt.figure(figsize=(18, 12))

    # ---- Plot 1: Toate cele 20 de tetraedre ----
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.set_title('20 Tetraedre in jurul unui Vertex\n(Configuratia "Pixel" din 600-cell)',
                  fontsize=11, fontweight='bold')

    # Paleta de culori
    cmap = plt.cm.hsv
    n_tetra = len(tetrahedra)
    colors = [cmap(i/n_tetra) for i in range(n_tetra)]

    for idx, tetra in enumerate(tetrahedra):
        faces_3d = get_tetrahedron_faces(tetra)

        # Desenam fetele cu transparenta
        poly = Poly3DCollection(faces_3d, alpha=0.4,
                               facecolor=colors[idx],
                               edgecolor='black',
                               linewidths=0.5)
        ax1.add_collection3d(poly)

    # Centrul
    ax1.scatter([0], [0], [0], c='red', s=200, zorder=10, label='Centru (vertex)')

    # Varfurile icosaedrului
    ax1.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2],
               c='gold', s=100, zorder=5, edgecolors='orange',
               label='Varfuri vecine')

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_box_aspect([1, 1, 1])
    ax1.legend(loc='upper left')

    # ---- Plot 2: Un singur tetraedru cu dimensiuni ----
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    ax2.set_title('Un Singur Tetraedru\n(Unitate fundamentala de spatiu)',
                  fontsize=11, fontweight='bold')

    tetra = tetrahedra[0]
    faces_3d = get_tetrahedron_faces(tetra)

    poly = Poly3DCollection(faces_3d, alpha=0.6,
                           facecolor='steelblue',
                           edgecolor='darkblue',
                           linewidths=2)
    ax2.add_collection3d(poly)

    # Etichete varfuri
    labels = ['O (centru)', 'A', 'B', 'C']
    for i, (v, label) in enumerate(zip(tetra, labels)):
        ax2.scatter([v[0]], [v[1]], [v[2]], c='red' if i==0 else 'gold',
                   s=150, zorder=10)
        ax2.text(v[0]*1.15, v[1]*1.15, v[2]*1.15, label, fontsize=10,
                fontweight='bold')

    # Aratam lungimea muchiei
    v0, v1 = np.array(tetra[0]), np.array(tetra[1])
    edge_len = np.linalg.norm(v1 - v0)
    mid = (v0 + v1) / 2
    ax2.text(mid[0], mid[1], mid[2], f'L={edge_len:.3f}', fontsize=9, color='red')

    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_box_aspect([1, 1, 1])

    # ---- Plot 3: Vedere de sus (proiectie 2D) ----
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_title('Vedere de Sus (Proiectie Z)\nSimetrie Icosaedrica',
                  fontsize=11, fontweight='bold')

    # Proiectam pe planul XY
    for idx, tetra in enumerate(tetrahedra):
        # Desenam triunghiul bazei
        base = [tetra[1], tetra[2], tetra[3], tetra[1]]  # Inchidem triunghiul
        xs = [v[0] for v in base]
        ys = [v[1] for v in base]
        ax3.fill(xs, ys, alpha=0.3, color=colors[idx])
        ax3.plot(xs, ys, 'k-', linewidth=0.5)

    # Varfuri
    ax3.scatter(vertices[:, 0], vertices[:, 1], c='gold', s=100,
               edgecolors='orange', zorder=5)
    ax3.scatter([0], [0], c='red', s=200, zorder=10)

    # Cercuri pentru referinta
    theta = np.linspace(0, 2*np.pi, 100)
    ax3.plot(np.cos(theta), np.sin(theta), 'b--', alpha=0.3, label='r=1')

    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.axis('equal')
    ax3.grid(True, alpha=0.3)

    # ---- Plot 4: Unghiuri si proportii ----
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.set_title('Proportiile Aurului in Tetraedru',
                  fontsize=11, fontweight='bold')
    ax4.axis('off')

    info_text = f"""
    TETRAEDRUL DIN 600-CELL
    =======================

    Varfuri: 4
    Fete: 4 (triunghiuri echilaterale)
    Muchii: 6

    DIMENSIUNI:
    -----------
    Lungime muchie (baza): {np.linalg.norm(np.array(tetra[1])-np.array(tetra[2])):.4f}
    Lungime muchie (la centru): {edge_len:.4f}

    UNGHIURI:
    ---------
    Unghi intre muchii la centru: {np.degrees(angle12):.2f}°
    Unghi diedru tetraedru regulat: 70.53°

    NUMERE CHEIE:
    -------------
    phi = {PHI:.6f}
    phi^2 = {PHI**2:.6f}
    1/phi = {1/PHI:.6f}

    arccos(1/phi^2) = {np.degrees(np.arccos(1/PHI**2)):.2f}°

    SEMNIFICATIE:
    -------------
    Cele 20 de tetraedre in jurul fiecarui
    vertex formeaza structura fundamentala
    a spatiului in teoria 600-cell.

    Fiecare "pixel" de spatiu = 20 tetraedre
    Numarul 20 apare in formula:
    1/alpha = 20 * phi^4 - 2*pi*alpha
    """

    ax4.text(0.05, 0.95, info_text, transform=ax4.transAxes,
            fontsize=10, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('tetrahedra_600cell.png', dpi=150, bbox_inches='tight')
    print("\nSalvat: tetrahedra_600cell.png")
    plt.show()

    # ===============================
    # FIGURA 2: ICOSAEDRU + DODECAEDRU
    # ===============================
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7),
                                    subplot_kw={'projection': '3d'})

    # ---- Icosaedru cu fete ----
    ax1.set_title('Icosaedru (20 fete)\nDual al Dodecaedrului', fontsize=11, fontweight='bold')

    for idx, face in enumerate(faces):
        i, j, k = face
        triangle = [vertices[i], vertices[j], vertices[k]]
        poly = Poly3DCollection([triangle], alpha=0.6,
                               facecolor=colors[idx],
                               edgecolor='black',
                               linewidths=1)
        ax1.add_collection3d(poly)

    ax1.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2],
               c='gold', s=100, edgecolors='black', zorder=5)

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_box_aspect([1, 1, 1])

    # ---- Dodecaedru (dual) ----
    ax2.set_title('Dodecaedru (12 fete pentagonale)\nVertexurile sunt centrele fetelor icosaedrului',
                  fontsize=11, fontweight='bold')

    # Centrele fetelor icosaedrului = varfurile dodecaedrului
    dodeca_vertices = []
    for face in faces:
        i, j, k = face
        center = (vertices[i] + vertices[j] + vertices[k]) / 3
        center = center / np.linalg.norm(center) * 0.9  # Normalizam
        dodeca_vertices.append(center)

    dodeca_vertices = np.array(dodeca_vertices)

    # Desenam varfurile
    ax2.scatter(dodeca_vertices[:, 0], dodeca_vertices[:, 1], dodeca_vertices[:, 2],
               c='purple', s=100, edgecolors='black', zorder=5)

    # Conectam varfurile apropiate
    for i in range(len(dodeca_vertices)):
        for j in range(i+1, len(dodeca_vertices)):
            d = np.linalg.norm(dodeca_vertices[i] - dodeca_vertices[j])
            if d < 0.6:  # Muchii scurte
                ax2.plot([dodeca_vertices[i, 0], dodeca_vertices[j, 0]],
                        [dodeca_vertices[i, 1], dodeca_vertices[j, 1]],
                        [dodeca_vertices[i, 2], dodeca_vertices[j, 2]],
                        'purple', linewidth=2, alpha=0.7)

    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_box_aspect([1, 1, 1])

    plt.tight_layout()
    plt.savefig('icosahedron_dodecahedron.png', dpi=150, bbox_inches='tight')
    print("Salvat: icosahedron_dodecahedron.png")
    plt.show()

if __name__ == "__main__":
    main()
