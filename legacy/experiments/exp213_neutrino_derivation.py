import numpy as np
from scipy import linalg

PHI = (1 + 5**0.5) / 2

def build_icosian_e8():
    """Builds the 240 roots of E8 from S (600-cell) and T (phi' * S)."""
    vertices = []
    # Type A
    for i in range(4):
        for s in [1, -1]:
            v = [0,0,0,0]; v[i] = s
            vertices.append(v)
    # Type B
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    vertices.append([s0, s1, s2, s3])
    # Type C
    import itertools
    base = [PHI/2, 0.5, 1/(2*PHI), 0]
    for p in itertools.permutations(range(4)):
        inv = 0
        for i in range(4):
            for j in range(i+1, 4):
                if p[i] > p[j]: inv += 1
        if inv % 2 == 0:
            for s in itertools.product([1, -1], repeat=3):
                v = [0,0,0,0]
                si = 0
                for i in range(4):
                    if base[p[i]] == 0: v[i] = 0
                    else:
                        v[i] = s[si] * base[p[i]]
                        si += 1
                vertices.append(v)
    
    unique_S = []
    for v in vertices:
        v = np.round(v, 8)
        if not any(np.allclose(v, u) for u in unique_S):
            unique_S.append(v)
    
    S = np.array(unique_S)
    T = -1/PHI * S
    return S, T

print("Building E8 sectors S and T...")
S, T = build_icosian_e8()

n_tau = 17
n_nu_base = 2 * n_tau
print(f"Base Neutrino exponent (2 * n_tau): {n_nu_base}")

# Icosahedron neighbors
phi = PHI
ico_verts = []
for s in [1, -1]:
    for t in [phi, -phi]:
        ico_verts.append([0, s, t])
        ico_verts.append([t, 0, s])
        ico_verts.append([s, t, 0])
ico_verts = np.array(ico_verts) / np.linalg.norm([0, 1, phi])

ico_adj = np.zeros((12, 12))
for i in range(12):
    for j in range(i+1, 12):
        if np.allclose(np.linalg.norm(ico_verts[i]-ico_verts[j]), np.linalg.norm(ico_verts[0]-ico_verts[1]), atol=0.1):
            ico_adj[i,j] = ico_adj[j,i] = 1

L_bag = np.zeros((13, 13))
L_bag[0, 1:] = L_bag[1:, 0] = -1
for i in range(12):
    for j in range(12):
        if ico_adj[i,j]: L_bag[i+1, j+1] = -1
for i in range(13):
    L_bag[i,i] = -np.sum(L_bag[i, :])

eigs = np.linalg.eigvalsh(L_bag)
print("\nSoliton Bag Eigenvalues:")
for e in eigs:
    found = False
    for a in range(-10, 15):
        for b in range(-10, 10):
            if np.isclose(e, a + b*PHI, atol=1e-5):
                print(f"  {e:.4f} = {a} + {b}*phi")
                found = True
                break
        if found: break
    if not found: print(f"  {e:.4f} (irrational)")

lambda_600 = 12 - 6*phi if (12 - 6*phi) > 0 else 2.2918
lambda_core = eigs[-1]
lambda_shell = eigs[-2]

print(f"\nBulk 600-cell Lambda: {lambda_600:.4f}")
shift_core = np.log(lambda_core/lambda_600)/np.log(PHI)
shift_shell = np.log(lambda_shell/lambda_600)/np.log(PHI)

print(f"Shift Core: {shift_core:.4f} (~ 3.5 or 4)")
print(f"Shift Shell: {shift_shell:.4f} (~ 3)")