import numpy as np

PHI = (1 + 5**0.5) / 2

# Experimental data (PDG 2024)
dm2_atm = 2.45e-3 # eV^2
dm2_sol = 7.50e-5 # eV^2

# Assume normal hierarchy: m1 < m2 < m3
m_e = 0.511e6 # eV

print("Neutrino Mass Scale Analysis")
print("-" * 30)

# Hypothesis 1: m = m_e * phi^-n
print("\nHypothesis 1: m = m_e * phi^-n")
for n in range(30, 45):
    m = m_e * PHI**(-n)
    print(f"n={n}: m={m:.6f} eV")

# Hypothesis 2: m_nu = m_e^3 / m_tau^2
m_tau = 1777e6 # eV
m_h1 = (m_e**3) / (m_tau**2)
print(f"\nHypothesis 2 (Double Seesaw): m = m_e^3 / m_tau^2 = {m_h1:.6f} eV")
print(f"Corresponding n = {-np.log(m_h1/m_e)/np.log(PHI):.2f} (Target 34)")

# Let's find n1, n2, n3 to match dm2
print("\nSearch for (n1, n2, n3) to match Delta m^2:")
best_err = 1e9
best_ns = (0,0,0)

for n3 in np.arange(32, 36, 0.1):
    for n2 in np.arange(n3+0.1, n3+10, 0.1):
        for n1 in np.arange(n2+0.1, n2+15, 0.1):
            m3 = m_e * PHI**(-n3)
            m2 = m_e * PHI**(-n2)
            m1 = m_e * PHI**(-n1)
            
            calc_dm2_atm = m3**2 - m2**2
            calc_dm2_sol = m2**2 - m1**2
            
            err = (calc_dm2_atm/dm2_atm - 1)**2 + (calc_dm2_sol/dm2_sol - 1)**2
            if err < best_err:
                best_err = err
                best_ns = (n1, n2, n3)

print(f"Best fit n: {best_ns}")
m1, m2, m3 = m_e * PHI**(-best_ns[0]), m_e * PHI**(-best_ns[1]), m_e * PHI**(-best_ns[2])
print(f"m3={m3:.6f} eV, m2={m2:.6f} eV, m1={m1:.6f} eV")
print(f"Delta m^2 atm calc: {m3**2-m2**2:.2e} (Target {dm2_atm:.2e})")
print(f"Delta m^2 sol calc: {m2**2-m1**2:.2e} (Target {dm2_sol:.2e})")

# Integer or half-integer check
print("\nHypothesis: m_nu = m_lepton * alpha^2")
alpha = 1/137.036
for name, ml in [("e", 0.511e6), ("mu", 105.66e6), ("tau", 1777e6)]:
    m_nu = ml * alpha**2
    print(f"nu_{name}: {m_nu:.6f} eV")

print("\nHypothesis: m_nu = m_lepton * alpha^3 * 2*pi")
for name, ml in [("e", 0.511e6), ("mu", 105.66e6), ("tau", 1777e6)]:
    m_nu = ml * alpha**3 * 2*np.pi
    print(f"nu_{name}: {m_nu:.6f} eV")