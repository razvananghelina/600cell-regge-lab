"""
FORMULE FUNDAMENTALE ALE FIZICII
================================
Toate unitatile sunt SI (metru, kilogram, secunda, amper, kelvin, mol, candela)
Numele variabilelor sunt consistente in tot fisierul.

Conventia de nume:
- Constante: MAJUSCULE (C, H, HBAR, G, etc.)
- Variabile: minuscule (m, r, t, v, etc.)
- Functii: snake_case (energia_cinetica, forta_gravitationala, etc.)
"""

import math
from typing import Tuple

# =============================================================================
# PARTEA 1: CONSTANTE FUNDAMENTALE (CODATA 2022)
# =============================================================================

# Viteza luminii in vid [m/s]
C = 299_792_458

# Constanta Planck [J·s]
H = 6.626_070_15e-34

# Constanta Planck redusa (h-bar) [J·s]
HBAR = H / (2 * math.pi)  # 1.054571817e-34

# Constanta gravitationala [m³/(kg·s²)]
G = 6.674_30e-11

# Sarcina elementara [C]
E_CHARGE = 1.602_176_634e-19

# Masa electronului [kg]
M_ELECTRON = 9.109_383_7015e-31

# Masa protonului [kg]
M_PROTON = 1.672_621_923_69e-27

# Masa neutronului [kg]
M_NEUTRON = 1.674_927_498_04e-27

# Constanta Boltzmann [J/K]
K_B = 1.380_649e-23

# Numarul lui Avogadro [1/mol]
N_A = 6.022_140_76e23

# Permitivitatea vidului [F/m]
EPSILON_0 = 8.854_187_8128e-12

# Permeabilitatea vidului [H/m]
MU_0 = 1.256_637_062_12e-6

# Constanta de structura fina (adimensionala)
ALPHA = 7.297_352_5693e-3  # ≈ 1/137.035999084

# Inversa constantei de structura fina
ALPHA_INV = 1 / ALPHA  # 137.035999084

# Constanta Rydberg [1/m]
R_INF = 10_973_731.568_160

# Raza Bohr [m]
A_0 = 5.291_772_109_03e-11

# Magnetonul Bohr [J/T]
MU_B = 9.274_010_0783e-24

# Masa Planck [kg]
M_PLANCK = math.sqrt(HBAR * C / G)  # 2.176e-8 kg

# Lungimea Planck [m]
L_PLANCK = math.sqrt(HBAR * G / C**3)  # 1.616e-35 m

# Timpul Planck [s]
T_PLANCK = math.sqrt(HBAR * G / C**5)  # 5.391e-44 s

# Energia Planck [J]
E_PLANCK = math.sqrt(HBAR * C**5 / G)  # 1.956e9 J

# Temperatura Planck [K]
TEMP_PLANCK = E_PLANCK / K_B  # 1.417e32 K

# Numarul de aur (phi)
PHI = (1 + math.sqrt(5)) / 2  # 1.6180339887...

# Pi
PI = math.pi


# =============================================================================
# PARTEA 2: MECANICA CLASICA
# =============================================================================

# 1. Legea a doua a lui Newton: F = m * a
def forta_newton(m: float, a: float) -> float:
    """
    F = m * a
    [N] = [kg] * [m/s²]
    """
    return m * a


# 2. Energia cinetica: E_k = (1/2) * m * v²
def energia_cinetica(m: float, v: float) -> float:
    """
    E_k = (1/2) * m * v²
    [J] = [kg] * [m/s]²
    """
    return 0.5 * m * v**2


# 3. Energia potentiala gravitationala (camp uniform): E_p = m * g * h
def energia_potentiala_uniform(m: float, g: float, h: float) -> float:
    """
    E_p = m * g * h
    [J] = [kg] * [m/s²] * [m]
    """
    return m * g * h


# 4. Energia potentiala gravitationala (Newton): U = -G * M * m / r
def energia_potentiala_gravitationala(M: float, m: float, r: float) -> float:
    """
    U = -G * M * m / r
    [J] = [m³/(kg·s²)] * [kg] * [kg] / [m]
    """
    return -G * M * m / r


# 5. Forta gravitationala: F = G * M * m / r²
def forta_gravitationala(M: float, m: float, r: float) -> float:
    """
    F = G * M * m / r²
    [N] = [m³/(kg·s²)] * [kg] * [kg] / [m²]
    """
    return G * M * m / r**2


# 6. Impuls: p = m * v
def impuls(m: float, v: float) -> float:
    """
    p = m * v
    [kg·m/s] = [kg] * [m/s]
    """
    return m * v


# 7. Moment cinetic: L = r * p (scalar, perpendicular)
def moment_cinetic(r: float, p: float) -> float:
    """
    L = r * p
    [kg·m²/s] = [m] * [kg·m/s]
    """
    return r * p


# 8. Lucru mecanic: W = F * d (forte colineare)
def lucru_mecanic(F: float, d: float) -> float:
    """
    W = F * d
    [J] = [N] * [m]
    """
    return F * d


# 9. Putere: P = W / t = F * v
def putere(F: float, v: float) -> float:
    """
    P = F * v
    [W] = [N] * [m/s]
    """
    return F * v


# 10. Forta elastica (Hooke): F = -k * x
def forta_elastica(k: float, x: float) -> float:
    """
    F = -k * x
    [N] = [N/m] * [m]
    """
    return -k * x


# 11. Perioada oscilatorului armonic: T = 2*pi*sqrt(m/k)
def perioada_oscilator(m: float, k: float) -> float:
    """
    T = 2*pi*sqrt(m/k)
    [s] = sqrt([kg] / [N/m]) = sqrt([kg] / [kg/s²])
    """
    return 2 * PI * math.sqrt(m / k)


# 12. Frecventa: f = 1/T
def frecventa(T: float) -> float:
    """
    f = 1/T
    [Hz] = 1/[s]
    """
    return 1 / T


# 13. Frecventa unghiulara: omega = 2*pi*f
def frecventa_unghiulara(f: float) -> float:
    """
    omega = 2*pi*f
    [rad/s] = [1/s]
    """
    return 2 * PI * f


# 14. Viteza pe orbita circulara: v = sqrt(G*M/r)
def viteza_orbitala(M: float, r: float) -> float:
    """
    v = sqrt(G*M/r)
    [m/s] = sqrt([m³/(kg·s²)] * [kg] / [m])
    """
    return math.sqrt(G * M / r)


# 15. Perioada orbitala (Kepler 3): T = 2*pi*sqrt(r³/(G*M))
def perioada_orbitala(r: float, M: float) -> float:
    """
    T = 2*pi*sqrt(r³/(G*M))
    [s] = sqrt([m³] / ([m³/(kg·s²)] * [kg]))
    """
    return 2 * PI * math.sqrt(r**3 / (G * M))


# =============================================================================
# PARTEA 3: ELECTROMAGNETISM
# =============================================================================

# 16. Legea lui Coulomb: F = k_e * q1 * q2 / r²
def forta_coulomb(q1: float, q2: float, r: float) -> float:
    """
    F = (1/(4*pi*epsilon_0)) * q1 * q2 / r²
    [N] = [N·m²/C²] * [C] * [C] / [m²]
    """
    k_e = 1 / (4 * PI * EPSILON_0)
    return k_e * q1 * q2 / r**2


# 17. Camp electric (sarcina punctiforma): E = k_e * q / r²
def camp_electric_punct(q: float, r: float) -> float:
    """
    E = (1/(4*pi*epsilon_0)) * q / r²
    [N/C] = [V/m]
    """
    k_e = 1 / (4 * PI * EPSILON_0)
    return k_e * q / r**2


# 18. Potential electric: V = k_e * q / r
def potential_electric(q: float, r: float) -> float:
    """
    V = (1/(4*pi*epsilon_0)) * q / r
    [V] = [J/C]
    """
    k_e = 1 / (4 * PI * EPSILON_0)
    return k_e * q / r


# 19. Energia potentiala electrica: U = k_e * q1 * q2 / r
def energia_potentiala_electrica(q1: float, q2: float, r: float) -> float:
    """
    U = (1/(4*pi*epsilon_0)) * q1 * q2 / r
    [J]
    """
    k_e = 1 / (4 * PI * EPSILON_0)
    return k_e * q1 * q2 / r


# 20. Legea lui Ohm: V = I * R
def tensiune_ohm(I: float, R: float) -> float:
    """
    V = I * R
    [V] = [A] * [Ohm]
    """
    return I * R


# 21. Putere electrica: P = I * V = I² * R = V² / R
def putere_electrica(I: float, V: float) -> float:
    """
    P = I * V
    [W] = [A] * [V]
    """
    return I * V


# 22. Capacitanta: C = Q / V
def capacitanta(Q: float, V: float) -> float:
    """
    C = Q / V
    [F] = [C] / [V]
    """
    return Q / V


# 23. Energia condensatorului: U = (1/2) * C * V²
def energia_condensator(C: float, V: float) -> float:
    """
    U = (1/2) * C * V²
    [J] = [F] * [V²]
    """
    return 0.5 * C * V**2


# 24. Inductanta bobinei: energie U = (1/2) * L * I²
def energia_bobina(L: float, I: float) -> float:
    """
    U = (1/2) * L * I²
    [J] = [H] * [A²]
    """
    return 0.5 * L * I**2


# 25. Forta Lorentz: F = q * (E + v x B), scalar pentru v perp B
def forta_lorentz_scalar(q: float, v: float, B: float) -> float:
    """
    F = q * v * B (pentru v perpendicular pe B)
    [N] = [C] * [m/s] * [T]
    """
    return q * v * B


# 26. Camp magnetic (fir infinit): B = mu_0 * I / (2*pi*r)
def camp_magnetic_fir(I: float, r: float) -> float:
    """
    B = mu_0 * I / (2*pi*r)
    [T] = [H/m] * [A] / [m]
    """
    return MU_0 * I / (2 * PI * r)


# 27. Flux magnetic: Phi = B * A (camp uniform, perpendicular)
def flux_magnetic(B: float, A: float) -> float:
    """
    Phi = B * A
    [Wb] = [T] * [m²]
    """
    return B * A


# 28. Legea inductiei Faraday: EMF = -d(Phi)/dt
def emf_faraday(dPhi: float, dt: float) -> float:
    """
    EMF = -dPhi/dt
    [V] = [Wb] / [s]
    """
    return -dPhi / dt


# 29. Viteza undelor EM: c = 1/sqrt(epsilon_0 * mu_0)
def viteza_em() -> float:
    """
    c = 1/sqrt(epsilon_0 * mu_0)
    [m/s]
    """
    return 1 / math.sqrt(EPSILON_0 * MU_0)


# 30. Densitatea energiei campului electric: u_E = (1/2) * epsilon_0 * E²
def densitate_energie_electric(E: float) -> float:
    """
    u_E = (1/2) * epsilon_0 * E²
    [J/m³]
    """
    return 0.5 * EPSILON_0 * E**2


# 31. Densitatea energiei campului magnetic: u_B = (1/2) * B² / mu_0
def densitate_energie_magnetic(B: float) -> float:
    """
    u_B = (1/2) * B² / mu_0
    [J/m³]
    """
    return 0.5 * B**2 / MU_0


# =============================================================================
# PARTEA 4: TERMODINAMICA
# =============================================================================

# 32. Ecuatia gazului ideal: P*V = n*R*T
R_GAZ = 8.314_462_618  # Constanta gazelor [J/(mol·K)]

def presiune_gaz_ideal(n: float, T: float, V: float) -> float:
    """
    P = n*R*T/V
    [Pa] = [mol] * [J/(mol·K)] * [K] / [m³]
    """
    return n * R_GAZ * T / V


# 33. Energia interna gaz ideal monoatomic: U = (3/2) * n * R * T
def energia_interna_monoatomic(n: float, T: float) -> float:
    """
    U = (3/2) * n * R * T
    [J]
    """
    return 1.5 * n * R_GAZ * T


# 34. Entropia (definitie Clausius): dS = dQ_rev / T
def variatie_entropie(dQ: float, T: float) -> float:
    """
    dS = dQ/T (reversibil)
    [J/K]
    """
    return dQ / T


# 35. Entropia Boltzmann: S = k_B * ln(W)
def entropia_boltzmann(W: float) -> float:
    """
    S = k_B * ln(W)
    [J/K] = [J/K] * [adimensional]
    W = numarul de microstari
    """
    return K_B * math.log(W)


# 36. Energia libera Helmholtz: F = U - T*S
def energia_helmholtz(U: float, T: float, S: float) -> float:
    """
    F = U - T*S
    [J] = [J] - [K] * [J/K]
    """
    return U - T * S


# 37. Energia libera Gibbs: G = H - T*S = U + P*V - T*S
def energia_gibbs(U: float, P: float, V: float, T: float, S: float) -> float:
    """
    G = U + P*V - T*S
    [J]
    """
    return U + P * V - T * S


# 38. Eficienta Carnot: eta = 1 - T_cold/T_hot
def eficienta_carnot(T_cold: float, T_hot: float) -> float:
    """
    eta = 1 - T_cold/T_hot
    [adimensional]
    """
    return 1 - T_cold / T_hot


# 39. Legea Stefan-Boltzmann: P = sigma * A * T⁴
SIGMA_SB = 5.670_374_419e-8  # [W/(m²·K⁴)]

def putere_radiatie_corp_negru(A: float, T: float) -> float:
    """
    P = sigma * A * T⁴
    [W] = [W/(m²·K⁴)] * [m²] * [K⁴]
    """
    return SIGMA_SB * A * T**4


# 40. Legea lui Wien: lambda_max * T = b
WIEN_B = 2.897_771_955e-3  # [m·K]

def lungime_unda_maxima(T: float) -> float:
    """
    lambda_max = b / T
    [m] = [m·K] / [K]
    """
    return WIEN_B / T


# =============================================================================
# PARTEA 5: RELATIVITATE SPECIALA
# =============================================================================

# 41. Factorul Lorentz: gamma = 1/sqrt(1 - v²/c²)
def gamma_lorentz(v: float) -> float:
    """
    gamma = 1/sqrt(1 - v²/c²)
    [adimensional]
    """
    beta = v / C
    return 1 / math.sqrt(1 - beta**2)


# 42. Dilatarea timpului: delta_t = gamma * delta_t_0
def dilatare_timp(delta_t_0: float, v: float) -> float:
    """
    delta_t = gamma * delta_t_0
    [s]
    delta_t_0 = timpul propriu
    """
    return gamma_lorentz(v) * delta_t_0


# 43. Contractia lungimii: L = L_0 / gamma
def contractie_lungime(L_0: float, v: float) -> float:
    """
    L = L_0 / gamma
    [m]
    L_0 = lungimea proprie
    """
    return L_0 / gamma_lorentz(v)


# 44. Energia relativista totala: E = gamma * m * c²
def energie_relativista_totala(m: float, v: float) -> float:
    """
    E = gamma * m * c²
    [J] = [kg] * [m/s]²
    """
    return gamma_lorentz(v) * m * C**2


# 45. Energia de repaus: E_0 = m * c²
def energie_repaus(m: float) -> float:
    """
    E_0 = m * c²
    [J] = [kg] * [m/s]²
    """
    return m * C**2


# 46. Impuls relativist: p = gamma * m * v
def impuls_relativist(m: float, v: float) -> float:
    """
    p = gamma * m * v
    [kg·m/s]
    """
    return gamma_lorentz(v) * m * v


# 47. Relatie energie-impuls: E² = (pc)² + (mc²)²
def energie_din_impuls(p: float, m: float) -> float:
    """
    E = sqrt((pc)² + (mc²)²)
    [J]
    """
    return math.sqrt((p * C)**2 + (m * C**2)**2)


# 48. Compunerea vitezelor: u = (v + w) / (1 + v*w/c²)
def compunere_viteze(v: float, w: float) -> float:
    """
    u = (v + w) / (1 + v*w/c²)
    [m/s]
    """
    return (v + w) / (1 + v * w / C**2)


# 49. Efectul Doppler relativist: f_obs = f_src * sqrt((1+beta)/(1-beta))
def doppler_relativist(f_src: float, v: float, approaching: bool = True) -> float:
    """
    f_obs = f_src * sqrt((1+beta)/(1-beta))  pentru apropiere
    f_obs = f_src * sqrt((1-beta)/(1+beta))  pentru deparare
    [Hz]
    """
    beta = v / C
    if approaching:
        return f_src * math.sqrt((1 + beta) / (1 - beta))
    else:
        return f_src * math.sqrt((1 - beta) / (1 + beta))


# 50. Interval invariant: ds² = c²dt² - dx² - dy² - dz²
def interval_spatiotemporal(dt: float, dx: float, dy: float, dz: float) -> float:
    """
    ds² = c²dt² - dx² - dy² - dz²
    [m²]
    """
    return (C * dt)**2 - dx**2 - dy**2 - dz**2


# =============================================================================
# PARTEA 6: RELATIVITATE GENERALA
# =============================================================================

# 51. Raza Schwarzschild: r_s = 2*G*M/c²
def raza_schwarzschild(M: float) -> float:
    """
    r_s = 2*G*M/c²
    [m] = [m³/(kg·s²)] * [kg] / [m/s]²
    """
    return 2 * G * M / C**2


# 52. Dilatarea gravitationala a timpului: dt_far = dt_local * sqrt(1 - r_s/r)
def dilatare_gravitationala(dt_local: float, M: float, r: float) -> float:
    """
    dt_far = dt_local / sqrt(1 - r_s/r)
    [s]
    """
    r_s = raza_schwarzschild(M)
    return dt_local / math.sqrt(1 - r_s / r)


# 53. Redshift gravitational: z = 1/sqrt(1 - r_s/r) - 1
def redshift_gravitational(M: float, r: float) -> float:
    """
    z = 1/sqrt(1 - r_s/r) - 1
    [adimensional]
    """
    r_s = raza_schwarzschild(M)
    return 1 / math.sqrt(1 - r_s / r) - 1


# 54. Temperatura Hawking: T = hbar*c³ / (8*pi*G*M*k_B)
def temperatura_hawking(M: float) -> float:
    """
    T = hbar*c³ / (8*pi*G*M*k_B)
    [K]
    """
    return HBAR * C**3 / (8 * PI * G * M * K_B)


# 55. Entropia Bekenstein-Hawking: S = A*k_B*c³ / (4*G*hbar)
def entropia_gaura_neagra(A: float) -> float:
    """
    S = A*k_B*c³ / (4*G*hbar)
    [J/K]
    A = aria orizontului de evenimente
    """
    return A * K_B * C**3 / (4 * G * HBAR)


# 56. Precesia Schwarzschild (per orbita): delta_phi = 6*pi*G*M / (a*c²*(1-e²))
def precesie_schwarzschild(M: float, a: float, e: float) -> float:
    """
    delta_phi = 6*pi*G*M / (a*c²*(1-e²))
    [rad per orbita]
    a = semiaxa mare, e = excentricitate
    """
    return 6 * PI * G * M / (a * C**2 * (1 - e**2))


# =============================================================================
# PARTEA 7: MECANICA CUANTICA
# =============================================================================

# 57. Energia fotonului: E = h*f = hbar*omega
def energia_foton(f: float) -> float:
    """
    E = h * f
    [J] = [J·s] * [Hz]
    """
    return H * f


# 58. Impulsul fotonului: p = h/lambda = E/c
def impuls_foton(wavelength: float) -> float:
    """
    p = h / lambda
    [kg·m/s] = [J·s] / [m]
    """
    return H / wavelength


# 59. Lungimea de unda de Broglie: lambda = h/p
def lungime_unda_de_broglie(p: float) -> float:
    """
    lambda = h / p
    [m] = [J·s] / [kg·m/s]
    """
    return H / p


# 60. Principiul incertitudinii Heisenberg: delta_x * delta_p >= hbar/2
def incertitudine_minima_pozitie(delta_p: float) -> float:
    """
    delta_x >= hbar / (2 * delta_p)
    [m]
    """
    return HBAR / (2 * delta_p)


# 61. Incertitudine energie-timp: delta_E * delta_t >= hbar/2
def incertitudine_minima_energie(delta_t: float) -> float:
    """
    delta_E >= hbar / (2 * delta_t)
    [J]
    """
    return HBAR / (2 * delta_t)


# 62. Energia nivelurilor hidrogenului: E_n = -13.6 eV / n²
E_RYDBERG_eV = 13.605_693_122_994  # [eV]
E_RYDBERG_J = E_RYDBERG_eV * E_CHARGE  # [J]

def energia_hidrogen(n: int) -> float:
    """
    E_n = -E_Rydberg / n²
    [J]
    n = numarul cuantic principal
    """
    return -E_RYDBERG_J / n**2


# 63. Raza orbitei Bohr: r_n = n² * a_0
def raza_bohr(n: int) -> float:
    """
    r_n = n² * a_0
    [m]
    """
    return n**2 * A_0


# 64. Frecventa tranzitiei: f = (E_i - E_f) / h
def frecventa_tranzitie(E_i: float, E_f: float) -> float:
    """
    f = (E_i - E_f) / h
    [Hz]
    """
    return (E_i - E_f) / H


# 65. Ecuatia Schrodinger (independenta de timp): H*psi = E*psi
# Operatorul Hamiltonian pentru particula libera: H = -hbar²/(2m) * d²/dx²
# (Aceasta e o ecuatie, nu o functie computabila direct)

# 66. Energia oscilatorului armonic cuantic: E_n = hbar*omega*(n + 1/2)
def energia_oscilator_cuantic(omega: float, n: int) -> float:
    """
    E_n = hbar * omega * (n + 1/2)
    [J]
    n = 0, 1, 2, ... (nivelul cuantic)
    """
    return HBAR * omega * (n + 0.5)


# 67. Efectul fotoelectric: E_k = h*f - W (functia de lucru)
def energia_cinetica_fotoelectron(f: float, W: float) -> float:
    """
    E_k = h*f - W
    [J]
    W = functia de lucru [J]
    """
    return H * f - W


# 68. Efectul Compton: delta_lambda = (h/(m_e*c)) * (1 - cos(theta))
COMPTON_WAVELENGTH = H / (M_ELECTRON * C)  # 2.426e-12 m

def compton_shift(theta: float) -> float:
    """
    delta_lambda = lambda_C * (1 - cos(theta))
    [m]
    theta = unghiul de imprastiere [rad]
    """
    return COMPTON_WAVELENGTH * (1 - math.cos(theta))


# 69. Densitatea de stari (3D, particule libere): g(E) = V*(2m)^(3/2) / (4*pi²*hbar³) * sqrt(E)
def densitate_stari_3d(V: float, m: float, E: float) -> float:
    """
    g(E) = V*(2m)^(3/2) / (4*pi²*hbar³) * sqrt(E)
    [1/J] = numar de stari per unitate de energie
    """
    return V * (2 * m)**1.5 / (4 * PI**2 * HBAR**3) * math.sqrt(E)


# 70. Distributia Fermi-Dirac: f(E) = 1 / (exp((E-mu)/(k_B*T)) + 1)
def distributie_fermi_dirac(E: float, mu: float, T: float) -> float:
    """
    f(E) = 1 / (exp((E-mu)/(k_B*T)) + 1)
    [adimensional] = probabilitatea de ocupare
    mu = potentialul chimic
    """
    x = (E - mu) / (K_B * T)
    return 1 / (math.exp(x) + 1)


# 71. Distributia Bose-Einstein: n(E) = 1 / (exp((E-mu)/(k_B*T)) - 1)
def distributie_bose_einstein(E: float, mu: float, T: float) -> float:
    """
    n(E) = 1 / (exp((E-mu)/(k_B*T)) - 1)
    [adimensional] = numarul mediu de ocupare
    mu = potentialul chimic (mu < E pentru bosoni)
    """
    x = (E - mu) / (K_B * T)
    return 1 / (math.exp(x) - 1)


# =============================================================================
# PARTEA 8: QED SI MODELUL STANDARD
# =============================================================================

# 72. Constanta de cuplaj QED: alpha = e² / (4*pi*epsilon_0*hbar*c)
def calcul_alpha() -> float:
    """
    alpha = e² / (4*pi*epsilon_0*hbar*c)
    [adimensional]
    """
    return E_CHARGE**2 / (4 * PI * EPSILON_0 * HBAR * C)


# 73. Lungimea Compton redusa: lambda_C_bar = hbar / (m*c)
def lungime_compton_redusa(m: float) -> float:
    """
    lambda_C_bar = hbar / (m*c)
    [m]
    """
    return HBAR / (m * C)


# 74. Raza clasica a electronului: r_e = alpha * a_0 = e²/(4*pi*epsilon_0*m_e*c²)
R_ELECTRON_CLASIC = E_CHARGE**2 / (4 * PI * EPSILON_0 * M_ELECTRON * C**2)

# 75. Relatia intre constante: a_0 = lambda_C_bar / alpha = hbar / (m_e * c * alpha)
def verifica_relatie_bohr_compton():
    """
    Verifica: a_0 = lambda_C_bar / alpha
    """
    lambda_c_bar = HBAR / (M_ELECTRON * C)
    return A_0, lambda_c_bar / ALPHA


# 76. Masa W boson [kg]
M_W_BOSON = 80.377e9 * E_CHARGE / C**2  # 80.377 GeV/c²

# 77. Masa Z boson [kg]
M_Z_BOSON = 91.1876e9 * E_CHARGE / C**2  # 91.1876 GeV/c²

# 78. Masa Higgs [kg]
M_HIGGS = 125.25e9 * E_CHARGE / C**2  # 125.25 GeV/c²

# 79. Constanta de cuplaj tare (la scala M_Z): alpha_s ≈ 0.1179
ALPHA_S = 0.1179

# 80. Constanta Fermi (interactie slaba): G_F / (hbar*c)³
G_FERMI = 1.166_378_7e-5  # [GeV^-2] in unitati naturale
G_FERMI_SI = G_FERMI * (1e9 * E_CHARGE)**(-2) * (HBAR * C)**3  # [J·m³]

# 81. Unghiul Weinberg: sin²(theta_W) ≈ 0.23122
SIN2_THETA_W = 0.231_22

# 82. Masa electron in eV
M_ELECTRON_eV = M_ELECTRON * C**2 / E_CHARGE  # 0.511 MeV

# 83. Masa muon [kg]
M_MUON = 105.6583755e6 * E_CHARGE / C**2  # 105.66 MeV/c²

# 84. Masa tau [kg]
M_TAU = 1776.86e6 * E_CHARGE / C**2  # 1776.86 MeV/c²


# =============================================================================
# PARTEA 9: COSMOLOGIE
# =============================================================================

# 85. Legea lui Hubble: v = H_0 * d
H_0_SI = 70 * 1000 / (3.086e22)  # 70 km/s/Mpc in [1/s]

def viteza_hubble(d: float) -> float:
    """
    v = H_0 * d
    [m/s] = [1/s] * [m]
    """
    return H_0_SI * d


# 86. Redshift cosmologic: z = (lambda_obs - lambda_emit) / lambda_emit
def redshift(lambda_obs: float, lambda_emit: float) -> float:
    """
    z = (lambda_obs - lambda_emit) / lambda_emit
    [adimensional]
    """
    return (lambda_obs - lambda_emit) / lambda_emit


# 87. Factor de scala din redshift: a = 1 / (1 + z)
def factor_scala(z: float) -> float:
    """
    a = 1 / (1 + z)
    [adimensional]
    a = factor de scala (a=1 astazi)
    """
    return 1 / (1 + z)


# 88. Densitatea critica: rho_c = 3*H²/(8*pi*G)
def densitate_critica(H: float) -> float:
    """
    rho_c = 3*H² / (8*pi*G)
    [kg/m³]
    H = constanta Hubble [1/s]
    """
    return 3 * H**2 / (8 * PI * G)


# 89. Parametrul de densitate: Omega = rho / rho_c
def parametru_densitate(rho: float, H: float) -> float:
    """
    Omega = rho / rho_c
    [adimensional]
    """
    rho_c = densitate_critica(H)
    return rho / rho_c


# 90. Ecuatia Friedmann (simplificata, plat): (da/dt)² = (8*pi*G/3) * rho * a²
# (Aceasta e o ecuatie diferentiala)

# 91. Temperatura CMB: T_CMB = 2.7255 K (astazi)
T_CMB = 2.7255  # [K]

# 92. Temperatura CMB in trecut: T(z) = T_0 * (1 + z)
def temperatura_cmb_z(z: float) -> float:
    """
    T(z) = T_CMB * (1 + z)
    [K]
    """
    return T_CMB * (1 + z)


# 93. Densitatea de energie a radiatiei: rho_r = a_r * T⁴ / c²
A_RADIATION = 4 * SIGMA_SB / C  # constanta de radiatie

def densitate_radiatie(T: float) -> float:
    """
    rho_r = a_r * T⁴ / c²
    [kg/m³]
    """
    return A_RADIATION * T**4 / C**2


# 94. Varsta universului (aproximare, materie dominanta): t = 2/(3*H)
def varsta_univers_materie(H: float) -> float:
    """
    t = 2 / (3*H)
    [s]
    """
    return 2 / (3 * H)


# =============================================================================
# PARTEA 10: RELATII IMPORTANTE SI CONSTANTE DERIVATE
# =============================================================================

# 95. Raza Bohr din constante fundamentale: a_0 = hbar / (m_e * c * alpha)
def calcul_raza_bohr() -> float:
    """
    a_0 = hbar / (m_e * c * alpha)
    [m]
    """
    return HBAR / (M_ELECTRON * C * ALPHA)


# 96. Energia Rydberg: E_R = m_e * c² * alpha² / 2
def calcul_energie_rydberg() -> float:
    """
    E_R = m_e * c² * alpha² / 2
    [J]
    """
    return M_ELECTRON * C**2 * ALPHA**2 / 2


# 97. Magnetonul Bohr: mu_B = e * hbar / (2 * m_e)
def calcul_magneton_bohr() -> float:
    """
    mu_B = e * hbar / (2 * m_e)
    [J/T]
    """
    return E_CHARGE * HBAR / (2 * M_ELECTRON)


# 98. Raportul maselor proton/electron
RAPORT_PROTON_ELECTRON = M_PROTON / M_ELECTRON  # ≈ 1836.15

# 99. Constanta de structura fina din alte constante
def alpha_din_constante() -> float:
    """
    alpha = e² / (4*pi*epsilon_0*hbar*c)
          = e² * mu_0 * c / (4*pi*hbar)
          = e² / (2*epsilon_0*h*c)
    [adimensional]
    """
    return E_CHARGE**2 / (4 * PI * EPSILON_0 * HBAR * C)


# 100. Relatie icosaedru-phi: 20 * phi^4 (pentru investigatia ta)
def icosaedru_phi() -> float:
    """
    20 * phi^4
    [adimensional]
    phi = (1 + sqrt(5)) / 2
    """
    return 20 * PHI**4


# =============================================================================
# VERIFICARI SI TESTE
# =============================================================================

def verifica_constante():
    """Verifica consistenta unor relatii fundamentale."""

    print("=== VERIFICARI CONSTANTE ===\n")

    # Viteza luminii din epsilon_0 si mu_0
    c_calculat = viteza_em()
    print(f"c din epsilon_0, mu_0: {c_calculat:.6e} m/s")
    print(f"c definit:             {C:.6e} m/s")
    print(f"Diferenta relativa:    {abs(c_calculat - C)/C:.2e}\n")

    # Alpha calculat
    alpha_calc = calcul_alpha()
    print(f"alpha calculat:  {alpha_calc:.10e}")
    print(f"alpha CODATA:    {ALPHA:.10e}")
    print(f"1/alpha:         {1/ALPHA:.6f}\n")

    # Raza Bohr
    a0_calc = calcul_raza_bohr()
    print(f"a_0 calculat:    {a0_calc:.6e} m")
    print(f"a_0 CODATA:      {A_0:.6e} m\n")

    # Energia Rydberg
    E_R = calcul_energie_rydberg()
    print(f"E_Rydberg calculat: {E_R:.6e} J")
    print(f"E_Rydberg in eV:    {E_R/E_CHARGE:.6f} eV\n")

    # Icosaedru
    ico = icosaedru_phi()
    print(f"20 * phi^4:      {ico:.6f}")
    print(f"1/alpha:         {ALPHA_INV:.6f}")
    print(f"Diferenta:       {abs(ico - ALPHA_INV):.6f}")
    print(f"Diferenta (%):   {abs(ico - ALPHA_INV)/ALPHA_INV * 100:.4f}%\n")

    # Mase Planck
    print(f"Masa Planck:     {M_PLANCK:.6e} kg")
    print(f"Lungime Planck:  {L_PLANCK:.6e} m")
    print(f"Timp Planck:     {T_PLANCK:.6e} s")
    print(f"Energie Planck:  {E_PLANCK:.6e} J")
    print(f"Temp Planck:     {TEMP_PLANCK:.6e} K\n")


def lista_constante():
    """Afiseaza toate constantele definite."""

    constante = {
        "C (viteza luminii)": (C, "m/s"),
        "H (Planck)": (H, "J·s"),
        "HBAR (Planck redus)": (HBAR, "J·s"),
        "G (gravitationala)": (G, "m³/(kg·s²)"),
        "E_CHARGE (sarcina elem.)": (E_CHARGE, "C"),
        "M_ELECTRON": (M_ELECTRON, "kg"),
        "M_PROTON": (M_PROTON, "kg"),
        "M_NEUTRON": (M_NEUTRON, "kg"),
        "K_B (Boltzmann)": (K_B, "J/K"),
        "N_A (Avogadro)": (N_A, "1/mol"),
        "EPSILON_0": (EPSILON_0, "F/m"),
        "MU_0": (MU_0, "H/m"),
        "ALPHA": (ALPHA, "adimensional"),
        "ALPHA_INV (1/alpha)": (ALPHA_INV, "adimensional"),
        "A_0 (Bohr)": (A_0, "m"),
        "PHI (nr. aur)": (PHI, "adimensional"),
        "M_PLANCK": (M_PLANCK, "kg"),
        "L_PLANCK": (L_PLANCK, "m"),
        "T_PLANCK": (T_PLANCK, "s"),
        "M_W_BOSON": (M_W_BOSON, "kg"),
        "M_Z_BOSON": (M_Z_BOSON, "kg"),
        "M_HIGGS": (M_HIGGS, "kg"),
        "20*PHI^4": (20 * PHI**4, "adimensional"),
    }

    print("=== CONSTANTE FUNDAMENTALE ===\n")
    for nume, (val, unit) in constante.items():
        print(f"{nume:25} = {val:.10e} [{unit}]")


if __name__ == "__main__":
    print("FORMULE FUNDAMENTALE ALE FIZICII")
    print("=" * 50)
    print()

    lista_constante()
    print()
    verifica_constante()
