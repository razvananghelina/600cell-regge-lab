import numpy as np

# Constants
PHI = (1 + np.sqrt(5)) / 2

def generate_600_cell_vertices():
    """
    Generates the 120 vertices of the 600-cell with edge length 1/PHI.
    Radius is 1.
    """
    vertices = []
    
    # 1. Permutations of (+-1, +-1, +-1, +-1) / 2
    for i in range(16):
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
        vertices.append(signs * 0.5)
        
    # 2. Permutations of (+-2, 0, 0, 0) / 2
    for i in range(4):
        v = np.zeros(4)
        v[i] = 1.0
        vertices.append(v)
        vertices.append(-v)
        
    # 3. Even permutations of (+-phi, +-1, +-1/phi, 0) / 2
    v_base = np.array([PHI, 1.0, 1.0/PHI, 0.0]) * 0.5
    from itertools import permutations
    perms = list(permutations([0, 1, 2, 3]))
    even_perms = []
    for p in perms:
        swaps = 0
        temp_p = list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if temp_p[i] > temp_p[j]:
                    swaps += 1
        if swaps % 2 == 0:
            even_perms.append(p)
            
    for p in even_perms:
        perm_vals = v_base[list(p)]
        for i in range(16):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(4)])
            vertices.append(perm_vals * signs)

    vertices = np.array(vertices)
    # Remove duplicates
    _, idx = np.unique(np.round(vertices, 8), axis=0, return_index=True)
    vertices = vertices[idx]
    return vertices

def analyze_loops(vertices):
    print(f"Total vertices: {len(vertices)}")
    
    # 1. Establish Edge Length (Chord)
    # Pick vertex 0 as reference
    v0 = vertices[0]
    
    # Distances to all other vertices
    dists = np.linalg.norm(vertices - v0, axis=1)
    
    # Filter 0
    dists = dists[dists > 0.00001]
    
    # Min distance = Edge Length
    chord_length = np.min(dists)
    print(f"Chord Length (Edge): {chord_length:.8f}")
    print(f"Check against 1/phi: {1/PHI:.8f}")
    
    # 2. Calculate Arc Length (Theta)
    # For a unit sphere (Radius = 1), Chord = 2 * sin(Theta/2)
    # Theta = 2 * arcsin(Chord / 2)
    
    radius = 1.0 # By construction
    theta_rad = 2 * np.arcsin(chord_length / (2 * radius))
    theta_deg = np.degrees(theta_rad)
    
    print(f"\n--- Geometric Angle of one Step ---")
    print(f"Theta (radians): {theta_rad:.8f}")
    print(f"Theta (degrees): {theta_deg:.8f}")
    
    # 3. Find the Great Circle (Geodesic Loop)
    # How many steps to return to origin along a straight line?
    # A full circle is 2*pi.
    # Number of steps = 2*pi / theta
    
    steps_float = (2 * np.pi) / theta_rad
    print(f"\n--- Loop Calculation ---")
    print(f"Theoretical steps for full circle (2*pi / theta): {steps_float:.8f}")
    
    steps_int = int(round(steps_float))
    deviation = abs(steps_float - steps_int)
    
    print(f"Integer Steps: {steps_int}")
    print(f"Deviation from integer: {deviation:.8f}")
    
    # 4. Verify Total Phase
    total_phase = steps_int * theta_rad
    print(f"\nTotal Loop Phase: {total_phase:.8f}")
    print(f"Target (2*pi):    {2 * np.pi:.8f}")
    
    if deviation < 1e-6:
        print("\nCONCLUSION: The geodesic loop is EXACTLY composed of integer steps.")
        print(f"The particle travels a discrete path of {steps_int} steps to complete a 2*pi rotation.")
    else:
        print("\nCONCLUSION: The loop does not close perfectly on integer vertices.")

if __name__ == "__main__":
    verts = generate_600_cell_vertices()
    analyze_loops(verts)
