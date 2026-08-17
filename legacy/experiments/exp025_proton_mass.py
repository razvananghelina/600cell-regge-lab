import numpy as np

# Constants
PHI = (1 + np.sqrt(5)) / 2
PI = np.pi
ALPHA = 0.00729735256
PROTON_ELECTRON_MASS_RATIO = 1836.15267343

def search_mass_ratio():
    print(f"Target Mass Ratio (mp/me): {PROTON_ELECTRON_MASS_RATIO}")
    
    # 1. The "Simon Plouffe" coincidence (often cited in numerology, let's verify)
    val1 = 6 * PI**5
    print(f"6 * pi^5: {val1} (Error: {abs(val1 - PROTON_ELECTRON_MASS_RATIO)})")
    
    # 2. Relations involving 600-cell numbers
    # We have 600 cells.
    # We have 120 vertices.
    # We have 20 faces per vertex.
    
    # Hypothesis: Mass is related to volume.
    # Proton = The whole 600-cell boundary? (600 tetrahedra)
    # Electron = Some minimal excitation?
    
    # Let's try 600 * something
    # 1836 / 600 = 3.06025
    # Is 3.06 related to Phi or Pi?
    # PI = 3.14
    # sqrt(10) = 3.16
    # 600 * PI = 1884 (too high)
    
    # 3. Relation with Alpha
    # Classic scaling: mp/me ~ 1/alpha ? No, 137 vs 1836.
    # Is it related to (20*phi^4)?
    # 137.08 * X = 1836
    # X = 13.39
    
    # 4. Looking for "Phi" relations
    # 1836.15
    # phi^15 = 1364 (too low)
    # phi^16 = 2207 (too high)
    
    # 5. The "Crystal" Hypothesis
    # Crystal volume = V_boundary = 50 * sqrt(2) * a^3 (from exp024)
    # with a = 1/phi
    
    a = 1/PHI
    v_boundary = 50 * np.sqrt(2) * a**3
    print(f"Boundary Volume (V3): {v_boundary}")
    
    # Ratio of V3 / Volume_Tetrahedron = 600.
    
    # What if the electron is related to the edge length?
    # Or surface area?
    
    # Check: 10 * (20*phi^4) * (4/3) ?
    # Check: (20*phi^4)^1.5 ?
    # 137^1.5 = 1603 (low)
    
    # Let's try to find 1836 using ONLY {20, phi, pi, 5, 2}
    
    # 6. Combinations
    term_alpha_inv = 20 * PHI**4
    
    # Maybe mp/me = 4/3 * term_alpha_inv * 10?
    # 1.33 * 137 * 10 = 1826 (Close!)
    
    val2 = (4/3) * term_alpha_inv * 10
    print(f"(4/3) * (20*phi^4) * 10: {val2}")
    
    # Why 4/3? Volume of sphere factor?
    # Why 10? The Decagon loop we found in EXP-023!
    
    # Refined formula:
    # mp/me = (4/3) * 10 * (1/alpha_bare)
    #       = 40/3 * 20 * phi^4
    #       = 800/3 * phi^4
    
    val3 = (800/3) * PHI**4
    print(f"800/3 * phi^4: {val3} (Error: {abs(val3 - PROTON_ELECTRON_MASS_RATIO)})")
    
    # Let's add the pi correction?
    # 1/alpha_phys = 20*phi^4 - 2*pi*alpha
    # Maybe use 1/alpha_phys instead?
    
    val4 = (40/3) * (1/ALPHA)
    print(f"40/3 * (1/alpha_phys): {val4}")
    
    # It's 1827.
    # We are missing about 9 units.
    
    # Is there a phi power that gives 9?
    # phi^4 = 6.85
    # phi^5 = 11.09
    
    # Try: 6 * pi^5 (Simon Plouffe again)
    # It is remarkably close.
    # 6 * pi^5 = 1836.118
    # Error is 0.034
    
    # Can we derive 6*pi^5 from 600-cell?
    # 600-cell has 120 vertices.
    # S^3 surface is 2*pi^2.
    # 120 * (pi^3 / 20) ?
    
    print("\n--- Search Conclusion ---")
    print("Most promising candidates:")
    print(f"1. 6 * pi^5 = {val1}")
    print(f"2. (40/3) * (1/alpha) = {val4} (Too low by ~0.5%)")

if __name__ == "__main__":
    search_mass_ratio()
