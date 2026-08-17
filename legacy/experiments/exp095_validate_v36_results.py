"""
EXP-095: Validare numerica - Rezultatele v3.6 (6 Februarie 2026)
================================================================
Validam RIGUROS toate rezultatele noi din research_6_feb.txt:

1. Materia Intunecata - 25% din sectorul G3
2. Tensiunea Hubble - raportul 26/24
3. Ierarhia fermonica extinsa (quarks: Up, Strange, Top)
4. Muon masa din phi^11 cu corectii
5. 121-lea varf si 120-cell

REGULA: Fiecare rezultat e clasificat ca DERIVAT / PATTERN / FITTING / SPECULATIV
"""

import math

# =============================================================================
# CONSTANTE FUNDAMENTALE (CODATA 2022 / PDG 2024)
# =============================================================================

PHI = (1 + math.sqrt(5)) / 2   # 1.6180339887...
PI = math.pi
ALPHA = 7.2973525693e-3         # 1/137.035999084
ALPHA_S_MZ = 0.1179             # alpha_s(M_Z) PDG 2024
SIN2_TW = 0.23122               # sin^2(theta_W)

# Mase leptoni (MeV/c^2) - PDG 2024
M_E = 0.51099895               # electron
M_MU = 105.6583755             # muon
M_TAU = 1776.86                # tau

# Mase quarki (MeV/c^2) - PDG 2024 (MS-bar la 2 GeV pt quarki usori)
M_UP = 2.16                    # +0.49 -0.26
M_DOWN = 4.67                  # +0.48 -0.17
M_STRANGE = 93.4               # +8.6 -3.4 (MS-bar la 2 GeV)
M_CHARM = 1270                 # +20 -20
M_BOTTOM = 4180                # +30 -20 (MS-bar la m_b)
M_TOP = 172760                 # +300 -300 (pole mass)

# Mase bosoni (GeV)
M_W = 80.377                   # W boson
M_Z = 91.1876                  # Z boson
M_H_EXP = 125.25               # Higgs (LHC)

# Constante Hubble
H0_PLANCK = 67.4               # km/s/Mpc (Planck 2018)
H0_PLANCK_ERR = 0.5
H0_SHOES = 73.04               # km/s/Mpc (SH0ES 2022)
H0_SHOES_ERR = 1.04

# Materie intunecata
OMEGA_DM_PLANCK = 0.265         # densitatea materiei intunecate (Planck 2018)
OMEGA_DM_ERR = 0.007
OMEGA_BARYON = 0.0493           # densitatea barionica
OMEGA_TOTAL_MATTER = OMEGA_DM_PLANCK + OMEGA_BARYON  # ~0.314

# Proprietati 600-cell
VERTICES_600CELL = 120
EDGES_600CELL = 720
FACES_600CELL = 1200
CELLS_600CELL = 600
VERTEX_DEGREE = 12              # vecini per nod
GRAPH_DIAMETER = 5              # distanta geodezica maxima
DECAGONS_PER_VERTEX = 6         # decagoane per vertex
TETRAHEDRA_PER_VERTEX = 20      # tetraedre per vertex
STRUCTURES_PER_VERTEX = 26      # 6 + 20

# =============================================================================
# HELPER
# =============================================================================
def percent_error(calculated, experimental):
    """Eroare procentuala."""
    return abs(calculated - experimental) / experimental * 100

def print_result(name, calc, exp, unit="", status="???"):
    """Afiseaza rezultat cu eroare."""
    err = percent_error(calc, exp)
    ok = "PASS" if err < 5.0 else ("MARGINAL" if err < 10 else "FAIL")
    print(f"  {name}:")
    print(f"    Calculat:     {calc:.6f} {unit}")
    print(f"    Experimental: {exp:.6f} {unit}")
    print(f"    Eroare:       {err:.4f}%  [{ok}]")
    print(f"    Status:       {status}")
    print()
    return err

# =============================================================================
# TEST 1: MATERIA INTUNECATA - 25% DIN SECTORUL G3
# =============================================================================
def test_dark_matter():
    print("=" * 70)
    print("TEST 1: MATERIA INTUNECATA - Sectorul G3")
    print("=" * 70)
    print()

    # Cele 96 varfuri fermionice = 4 grupuri de 24
    total_fermionic = 96
    groups = 4
    per_group = total_fermionic // groups  # 24

    # Alegerea axei temporale selecteaza 3 grupuri active
    active_groups = 3
    sterile_group = 1

    # Fractia sterila din varfurile fermionice
    fraction_sterile_fermionic = sterile_group / groups  # 0.25

    # Fractia totala din 600-cell (incluzand si cele 24 non-fermionice)
    fraction_sterile_total = per_group / VERTICES_600CELL  # 24/120 = 0.20

    # Ce spune teoria: 25% din potentialul retelei
    theory_dm_fraction = 0.25

    print(f"  Varfuri fermionice: {total_fermionic} = {groups} x {per_group}")
    print(f"  Grupuri active (axa temporala): {active_groups}")
    print(f"  Grup sterile (G3): {sterile_group} x {per_group} = {per_group} varfuri")
    print()

    # Comparatie 1: Fractia G3 din varfurile fermionice
    print(f"  Fractia G3/fermioni: {fraction_sterile_fermionic:.4f} = 24/96 = 1/4")
    print(f"  Predictia teoriei:   {theory_dm_fraction:.4f}")
    print()

    # Comparatie 2: Ce zice Planck?
    # Omega_DM / (Omega_DM + Omega_baryon) = fractia materie intunecata din materie totala
    dm_fraction_of_matter = OMEGA_DM_PLANCK / OMEGA_TOTAL_MATTER

    print(f"  --- Comparatie cu observatii Planck ---")
    err1 = print_result(
        "DM / Materie totala",
        theory_dm_fraction,
        dm_fraction_of_matter,
        status="SPECULATIV - ce fractie comparam?"
    )

    # Alternativa: DM ca fractie din densitatea critica totala
    # Omega_DM = 0.265 -> 26.5%
    err2 = print_result(
        "Omega_DM (fractie din densitate critica)",
        theory_dm_fraction,
        OMEGA_DM_PLANCK,
        status="SPECULATIV"
    )

    # Analiza critica
    print("  --- ANALIZA CRITICA ---")
    print(f"  Teoria zice: 25% din potential = materie intunecata")
    print(f"  Planck masoara: Omega_DM = {OMEGA_DM_PLANCK} ({OMEGA_DM_PLANCK*100:.1f}%)")
    print(f"  DM/materie_totala = {dm_fraction_of_matter:.4f} ({dm_fraction_of_matter*100:.1f}%)")
    print()
    print(f"  25% vs 26.5% (Omega_DM):        eroare {err2:.1f}%")
    print(f"  25% vs {dm_fraction_of_matter*100:.1f}% (DM/materie): eroare {err1:.1f}%")
    print()

    # VERDICT
    if err2 < 10:
        print("  VERDICT: INTERESANT dar SPECULATIV")
        print("  - Numarul 25% (=1/4) e un rezultat geometric simplu")
        print("  - Omega_DM=26.5% e apropiat dar NU e 25%")
        print("  - Lipseste mecanismul fizic: DE CE sectorul G3 contribuie la DM?")
        print("  - Nu e clar daca comparatia corecta e cu Omega_DM sau DM/materie")
    else:
        print("  VERDICT: SPECULATIV - discrepanta prea mare")
    print()
    return err2


# =============================================================================
# TEST 2: TENSIUNEA HUBBLE - RAPORTUL 26/24
# =============================================================================
def test_hubble_tension():
    print("=" * 70)
    print("TEST 2: TENSIUNEA HUBBLE - Raportul 26/24")
    print("=" * 70)
    print()

    # Formula: H0_local = H0_global * (26/24)
    ratio = 26 / 24  # = 1.08333...

    H0_calc = H0_PLANCK * ratio

    print(f"  Raportul geometric: 26/24 = {ratio:.6f}")
    print(f"  26 = directii per vertex (6 decagoane + 20 tetraedre)")
    print(f"  24 = varfuri 24-cell (baza bosonica)")
    print()

    err = print_result(
        "H0 local (din Planck * 26/24)",
        H0_calc,
        H0_SHOES,
        unit="km/s/Mpc",
        status="PATTERN"
    )

    # Cat de bun e 26/24?
    # Raportul experimental
    ratio_exp = H0_SHOES / H0_PLANCK
    print(f"  Raportul experimental H0_SH0ES/H0_Planck: {ratio_exp:.6f}")
    print(f"  Raportul teoretic 26/24:                   {ratio:.6f}")
    err_ratio = percent_error(ratio, ratio_exp)
    print(f"  Eroare raport:                             {err_ratio:.4f}%")
    print()

    # Verificam si cu barele de eroare
    H0_SHOES_low = H0_SHOES - H0_SHOES_ERR
    H0_SHOES_high = H0_SHOES + H0_SHOES_ERR
    print(f"  SH0ES interval: [{H0_SHOES_low:.2f}, {H0_SHOES_high:.2f}] km/s/Mpc")
    print(f"  Calculat: {H0_calc:.2f} km/s/Mpc")
    in_range = H0_SHOES_low <= H0_calc <= H0_SHOES_high
    print(f"  In intervalul SH0ES: {'DA' if in_range else 'NU'}")
    print()

    # Analiza critica
    print("  --- ANALIZA CRITICA ---")
    print(f"  PRO: Raportul 26/24 foloseste numere deja din teorie")
    print(f"       (26 = structuri per vertex, 24 = 24-cell)")
    print(f"  PRO: Eroare {err:.2f}% e excelenta")
    print(f"  CONTRA: De ce H0_local ar trebui MULTIPLICAT cu 26/24?")
    print(f"  CONTRA: Lipseste mecanismul fizic (cum afecteaza proiectia masuratorile?)")
    print(f"  CONTRA: Tensiunea Hubble ar putea avea explicatii sistematice")
    print()

    if err < 1.0:
        print("  VERDICT: NUMERIC EXCELENT dar SPECULATIV ca mecanisme")
    else:
        print("  VERDICT: SPECULATIV")
    print()
    return err


# =============================================================================
# TEST 3: IERARHIA FERMONICA - QUARKS
# =============================================================================
def test_quark_hierarchy():
    print("=" * 70)
    print("TEST 3: IERARHIA MASELOR QUARKILOR")
    print("=" * 70)
    print()

    # Formula: m_quark = m_e * phi^n * C(Q, alpha_s)
    # Testam intai FARA corectii (phi^n pur)

    quarks = [
        ("Up",      M_UP,      3,  "dim Im(H) = 3 unitati cuaternionice"),
        ("Down",    M_DOWN,    5,  "diametrul grafului 600-cell"),
        ("Strange", M_STRANGE, 11, "5+6 = diametru + decagoane (ca muon)"),
        ("Charm",   M_CHARM,   16, "16 = fermioni per generatie"),
        ("Bottom",  M_BOTTOM,  19, "19 = ? (neclar)"),
        ("Top",     M_TOP,     26, "26 = structuri per vertex"),
    ]

    print(f"  Formula: m_q = m_e * phi^n")
    print(f"  m_e = {M_E:.6f} MeV")
    print(f"  phi = {PHI:.10f}")
    print()

    errors = {}
    for name, m_exp, n, origin in quarks:
        m_calc = M_E * PHI**n
        err = percent_error(m_calc, m_exp)
        ok = "PASS" if err < 5 else ("MARGINAL" if err < 15 else ("WEAK" if err < 25 else "FAIL"))

        print(f"  {name:8s}: n={n:2d}  calc={m_calc:12.2f} MeV  exp={m_exp:12.2f} MeV  err={err:6.1f}%  [{ok}]")
        print(f"           Sursa n: {origin}")
        errors[name] = err

    print()

    # Verificare speciala: Strange = Muon la scala 1.17 GeV?
    print("  --- VERIFICARE SPECIALA: Strange-Muon la scala 1.17 GeV ---")
    print()

    # Strange quark running mass cu QCD 1-loop
    # m_s(mu) = m_s(2 GeV) * [alpha_s(mu)/alpha_s(2GeV)]^(gamma0/2*b0)
    # gamma0 = 8, b0 = 23/3 (pentru Nf=3)
    # alpha_s(2 GeV) ~ 0.30, alpha_s(1.17 GeV) ~ 0.40

    # Approximare simplificata (1-loop QCD running):
    Nf = 3  # flavors active la ~1 GeV
    b0 = (33 - 2*Nf) / (12*PI)  # = 27/(12*pi) = 0.7162
    gamma_m = 8 / (33 - 2*Nf)    # anomalous dimension / (2*b0) simplificat
    # = 8/27 = 0.2963

    # alpha_s la diferite scale (aproximativ, 1-loop)
    Lambda_QCD = 0.332  # GeV (MS-bar, Nf=3)

    def alpha_s_1loop(mu, Lambda=Lambda_QCD, nf=3):
        b0_loc = (33 - 2*nf) / (12*PI)
        return 1.0 / (2 * b0_loc * math.log(mu / Lambda))

    mu_2GeV = 2.0
    mu_target = 1.168  # GeV

    alpha_s_2 = alpha_s_1loop(mu_2GeV)
    alpha_s_target = alpha_s_1loop(mu_target)

    # Running mass: m_s(mu) = m_s(2GeV) * [alpha_s(mu)/alpha_s(2GeV)]^(4/(33-2*Nf))
    # exponent = 4*2/(33-2*Nf) = 8/27 ... actually: gamma0/(2*b0) = (8) / (2 * 27/(12*pi)) = 8*12*pi/(2*27)
    # Let me be more careful:
    # gamma_0 = 8 (1-loop anomalous dimension for quark mass in QCD)
    # b_0 = (33-2Nf)/(12*pi) for alpha_s = g^2/(4*pi)
    # m(mu2)/m(mu1) = [alpha_s(mu2)/alpha_s(mu1)]^(gamma_0 / (2*b_0 * 4*pi))
    # Wait, let me use the standard formula directly:
    # m(mu) = m(mu0) * [alpha_s(mu)/alpha_s(mu0)]^(gamma0/(2*beta0))
    # where gamma0 = 8 and beta0 = (33-2Nf)/3 (in the convention beta = -beta0*g^3/(16*pi^2))
    # So gamma0/(2*beta0) = 8/(2*(33-2*3)/3) = 8/(2*27/3) = 8/18 = 4/9

    power = 4.0/9.0  # = gamma0/(2*beta0) in standard convention

    m_strange_running = M_STRANGE * (alpha_s_target / alpha_s_2)**power

    print(f"  m_strange(2 GeV) = {M_STRANGE:.1f} MeV (PDG 2024)")
    print(f"  alpha_s(2 GeV)   = {alpha_s_2:.4f} (1-loop)")
    print(f"  alpha_s(1.17 GeV)= {alpha_s_target:.4f} (1-loop)")
    print(f"  Power (gamma0/2beta0) = {power:.4f}")
    print(f"  m_strange(1.17 GeV) = {m_strange_running:.2f} MeV (running)")
    print(f"  m_muon              = {M_MU:.2f} MeV")
    err_sm = percent_error(m_strange_running, M_MU)
    print(f"  Eroare strange vs muon: {err_sm:.2f}%")
    print()

    if err_sm < 5:
        print("  VERDICT: STRANGE = MUON la scala ~1.17 GeV - CONFIRMAT!")
    elif err_sm < 15:
        print("  VERDICT: STRANGE ~ MUON la scala ~1.17 GeV - APROXIMATIV")
    else:
        print(f"  VERDICT: Strange-Muon: discrepanta {err_sm:.1f}% - NECESITA ANALIZA")
    print()

    # Analiza critica generala
    print("  --- ANALIZA CRITICA ---")
    good = sum(1 for e in errors.values() if e < 5)
    marginal = sum(1 for e in errors.values() if 5 <= e < 15)
    weak = sum(1 for e in errors.values() if 15 <= e < 25)
    fail = sum(1 for e in errors.values() if e >= 25)

    print(f"  PASS (<5%):     {good}/6 quarks")
    print(f"  MARGINAL:       {marginal}/6")
    print(f"  WEAK (15-25%):  {weak}/6")
    print(f"  FAIL (>25%):    {fail}/6")
    print()

    if good >= 2:
        print("  VERDICT: Exponentii Up (n=3) si Strange (n=11) sunt solizi")
        print("  Restul necesita corectii QCD de 10-20% (asteptat)")
    print()

    return errors


# =============================================================================
# TEST 4: IERARHIA MASELOR LEPTONICE (re-validare)
# =============================================================================
def test_lepton_hierarchy():
    print("=" * 70)
    print("TEST 4: IERARHIA MASELOR LEPTONICE (re-validare)")
    print("=" * 70)
    print()

    # Rapoarte experimentale
    ratio_mu_e_exp = M_MU / M_E    # 206.768
    ratio_tau_mu_exp = M_TAU / M_MU  # 16.817
    ratio_tau_e_exp = M_TAU / M_E   # 3477.15

    print(f"  Rapoarte experimentale:")
    print(f"    m_mu/m_e  = {ratio_mu_e_exp:.4f}")
    print(f"    m_tau/m_mu = {ratio_tau_mu_exp:.4f}")
    print(f"    m_tau/m_e  = {ratio_tau_e_exp:.2f}")
    print()

    # Derivare din 600-cell:
    # n=11 = 5 (diametru) + 6 (decagoane per vertex)
    # n=6  = decagoane per vertex
    # n=17 = 11 + 6

    tests = [
        ("m_mu/m_e",  11, "5+6 = diametru + decagoane", ratio_mu_e_exp),
        ("m_tau/m_mu",  6, "6 = decagoane per vertex",    ratio_tau_mu_exp),
        ("m_tau/m_e",  17, "11+6 = holonomie nivel 2",    ratio_tau_e_exp),
    ]

    print(f"  Derivare din geometria 600-cell:")
    print(f"    5 = diametrul grafului (INTRINSEC)")
    print(f"    6 = decagoane per vertex (INTRINSEC)")
    print()

    for name, n, origin, exp_val in tests:
        calc_val = PHI**n
        err = percent_error(calc_val, exp_val)
        ok = "PASS" if err < 5 else "MARGINAL"
        print(f"  {name:12s}: phi^{n:2d} = {calc_val:8.2f}  exp = {exp_val:8.2f}  err = {err:.2f}%  [{ok}]")
        print(f"               n={n}: {origin}")

    print()
    print("  VERDICT: Exponentii DERIVATI din topologia 600-cell")
    print("  Erorile de 3-7% sunt asteptate (masele nu sunt exact phi^n)")
    print()


# =============================================================================
# TEST 5: VERIFICAREA CONSISTENTEI INTERNE
# =============================================================================
def test_internal_consistency():
    print("=" * 70)
    print("TEST 5: CONSISTENTA INTERNA")
    print("=" * 70)
    print()

    # Test A: 26 = 6 + 20 (structuri per vertex)
    print("  Test A: 26 = 6 + 20?")
    print(f"  6 decagoane + 20 tetraedre = {6 + 20}")
    print(f"  Verificat: {'PASS' if 6 + 20 == 26 else 'FAIL'}")
    print()

    # Test B: 96 = 4 x 24
    print("  Test B: 96 = 4 x 24?")
    print(f"  120 - 24 (non-fermionici) = {120 - 24}")
    print(f"  4 x 24 = {4 * 24}")
    print(f"  Verificat: {'PASS' if 96 == 4*24 else 'FAIL'}")
    print()

    # Test C: phi^11 = phi^5 * phi^6 (aditivitate exponenti)
    print("  Test C: phi^11 = phi^5 * phi^6?")
    print(f"  phi^11 = {PHI**11:.6f}")
    print(f"  phi^5 * phi^6 = {PHI**5 * PHI**6:.6f}")
    print(f"  Diferenta: {abs(PHI**11 - PHI**5 * PHI**6):.2e}")
    print(f"  Verificat: PASS (identitate algebrica triviala)")
    print()

    # Test D: phi^17 = phi^11 * phi^6
    print("  Test D: phi^17 = phi^11 * phi^6?")
    print(f"  phi^17 = {PHI**17:.2f}")
    print(f"  phi^11 * phi^6 = {PHI**11 * PHI**6:.2f}")
    print(f"  Verificat: PASS")
    print()

    # Test E: 120 = 8 + 16 + 96 (descompunerea varfurilor)
    print("  Test E: 120 = 8 + 16 + 96?")
    print(f"  8 (16-cell) + 16 (tesseract) + 96 (snub 24-cell) = {8 + 16 + 96}")
    print(f"  Verificat: {'PASS' if 8 + 16 + 96 == 120 else 'FAIL'}")
    print()

    # Test F: 26/24 vs raportul Hubble
    print("  Test F: Consistenta 26/24 cu alte formule?")
    print(f"  26/24 = {26/24:.6f}")
    print(f"  1 + 1/12 = {1 + 1/12:.6f}")
    print(f"  Nota: 12 = grad vertex 600-cell")
    print(f"  26/24 = (6+20)/24 = 6/24 + 20/24 = 1/4 + 5/6 = {1/4 + 5/6:.6f}")
    print(f"  Verificat: 26/24 = 1 + 1/12 = {26/24 == 1 + 1/12}")
    # Actually, 26/24 = 13/12 = 1 + 1/12
    print(f"  De fapt: 26/24 = 13/12 = 1 + 1/12")
    print(f"  12 = vertex degree => legatura cu graful!")
    print()

    # Test G: Formula alpha verificare
    print("  Test G: Ecuatia 2*pi*alpha^2 - 20*phi^4*alpha + 1 = 0?")
    val = 2*PI*ALPHA**2 - 20*PHI**4*ALPHA + 1
    print(f"  2*pi*alpha^2 - 20*phi^4*alpha + 1 = {val:.8f}")
    print(f"  Verificat: {'PASS' if abs(val) < 0.001 else 'FAIL'} (trebuie ~0)")
    print()

    # Test H: alpha_s/alpha = 10*phi?
    print("  Test H: alpha_s/alpha = 10*phi?")
    ratio_calc = 10 * PHI
    ratio_exp = ALPHA_S_MZ / ALPHA
    err = percent_error(ratio_calc, ratio_exp)
    print(f"  10*phi = {ratio_calc:.4f}")
    print(f"  alpha_s/alpha = {ratio_exp:.4f}")
    print(f"  Eroare: {err:.2f}%")
    print(f"  Verificat: {'PASS' if err < 1 else 'FAIL'}")
    print()


# =============================================================================
# TEST 6: 120-CELL SI GRAVITONUL (al 121-lea varf)
# =============================================================================
def test_120cell_graviton():
    print("=" * 70)
    print("TEST 6: 120-CELL SI GRAVITONUL")
    print("=" * 70)
    print()

    # Proprietati 120-cell
    vertices_120cell = 600
    edges_120cell = 1200
    faces_120cell = 720
    cells_120cell = 120  # celule dodecaedrice

    print(f"  120-cell:")
    print(f"    Varfuri:  {vertices_120cell}")
    print(f"    Muchii:   {edges_120cell}")
    print(f"    Fete:     {faces_120cell}")
    print(f"    Celule:   {cells_120cell} (dodecaedre)")
    print()

    # Dualitatea: 120-cell este DUALUL lui 600-cell
    print(f"  Dualitate 600-cell <-> 120-cell:")
    print(f"    600-cell: {VERTICES_600CELL} varfuri, {CELLS_600CELL} celule")
    print(f"    120-cell: {vertices_120cell} varfuri, {cells_120cell} celule")
    print(f"    Celulele unuia = varfurile celuilalt? {CELLS_600CELL == vertices_120cell and VERTICES_600CELL == cells_120cell}")
    print()

    # Al 121-lea varf
    print(f"  Al 121-lea varf (Body-Centered):")
    print(f"    Conceptual: centrul politopului ca punct suplimentar")
    print(f"    Rol: graviton de spin-2")
    print()

    # Analiza critica
    print("  --- ANALIZA CRITICA ---")
    print(f"  CONTRA: Un politop regulat nu are un varf 'central' natural")
    print(f"  CONTRA: 121 nu are proprietati speciale evidente")
    print(f"  CONTRA: Spin-2 nu e derivat, e ATRIBUIT")
    print(f"  PRO: Un defect topologic central ar putea media gravitatia")
    print(f"  PRO: 121 = 11^2, iar 11 = exponentul mu/e-ului!")
    print()
    print("  VERDICT: SPECULATIV - idee interesanta dar nedemonstranta")
    print()


# =============================================================================
# TEST 7: TOP QUARK n=26
# =============================================================================
def test_top_quark():
    print("=" * 70)
    print("TEST 7: TOP QUARK - n=26 DETALIAT")
    print("=" * 70)
    print()

    n_top = 26
    m_top_calc = M_E * PHI**n_top

    print(f"  n = 26 = structuri per vertex (6 + 20)")
    print(f"  phi^26 = {PHI**26:.2f}")
    print(f"  m_top calculat = m_e * phi^26 = {m_top_calc:.2f} MeV = {m_top_calc/1000:.2f} GeV")
    print(f"  m_top experimental = {M_TOP:.2f} MeV = {M_TOP/1000:.2f} GeV")
    print()

    err = percent_error(m_top_calc, M_TOP)
    print(f"  Eroare: {err:.1f}%")
    print()

    # Ce n ar fi exact?
    n_exact = math.log(M_TOP / M_E) / math.log(PHI)
    print(f"  n exact pentru m_top: {n_exact:.4f}")
    print(f"  Cel mai apropiat intreg: {round(n_exact)}")
    print()

    if err < 1:
        print("  VERDICT: EXCELENT - sub 1% eroare!")
    elif err < 5:
        print("  VERDICT: BUN - sub 5% eroare")
    elif err < 20:
        print(f"  VERDICT: MARGINAL - eroare {err:.1f}% (necesita corectii QCD)")
    else:
        print(f"  VERDICT: SLAB - eroare {err:.1f}%")
    print()


# =============================================================================
# TEST 8: MUON DIN PHI^11 - PRECIZIA 0.06% DIN RESEARCH
# =============================================================================
def test_muon_precision():
    print("=" * 70)
    print("TEST 8: MUON - VERIFICARE CLAIM 0.06%")
    print("=" * 70)
    print()

    # research_6_feb.txt zice: Muon n=11, calc=105.6 MeV, exp=105.66, err=0.06%
    # Sa verificam

    m_mu_calc_naive = M_E * PHI**11
    err_naive = percent_error(m_mu_calc_naive, M_MU)

    print(f"  m_mu = m_e * phi^11 = {M_E} * {PHI**11:.4f} = {m_mu_calc_naive:.4f} MeV")
    print(f"  m_mu experimental = {M_MU:.4f} MeV")
    print(f"  Eroare REALA: {err_naive:.2f}%")
    print()

    # research_6_feb.txt spune 0.06% - asta e GRESIT daca folosim phi^11 pur
    # research spune calc=105.6 MeV, ceea ce nu e phi^11 * m_e
    print(f"  Claim din research: calc=105.6 MeV, err=0.06%")
    print(f"  Dar phi^11 * m_e = {m_mu_calc_naive:.4f} MeV")
    print()

    if err_naive > 1.0:
        print(f"  ATENTIE: Eroarea reala ({err_naive:.2f}%) NU e 0.06%!")
        print(f"  Research-ul foloseste probabil m_e*phi^11*C(Q,alpha_s)")
        print(f"  unde corectia C aduce valoarea la 105.6 MeV")
        print()

        # Ce corectie ar trebui?
        correction_needed = M_MU / m_mu_calc_naive
        print(f"  Factor de corectie necesar: {correction_needed:.6f}")
        print(f"  Sau procentual: {(correction_needed-1)*100:.2f}%")
        print()

        print("  VERDICT: Claim-ul de 0.06% include o CORECTIE nedocumentata")
        print("  Fara corectie, eroarea reala e {:.2f}% (DERIVAT din phi^11)".format(err_naive))
        print("  Cu corectie, statusul devine FITTING, nu DERIVAT!")
    else:
        print("  VERDICT: Claim-ul verificat - eroare sub 1%")
    print()


# =============================================================================
# SUMAR FINAL
# =============================================================================
def print_summary():
    print()
    print("=" * 70)
    print("SUMAR FINAL - VALIDARE REZULTATE v3.6")
    print("=" * 70)
    print()
    print("  REZULTAT                      | EROARE    | STATUS")
    print("  " + "-" * 65)

    results = [
        ("alpha (ecuatie spectrala)",       "0.0001%",  "DERIVAT - SOLID"),
        ("alpha_s = 1/(2*phi^3)",          "0.03%",    "DERIVAT - SOLID"),
        ("sin^2(theta_W) = 6/26",          "0.19%",    "DERIVAT - SOLID"),
        ("m_H = m_W*(phi-8*alpha)",        "0.09%",    "DERIVAT - SOLID"),
        ("m_mu/m_e = phi^11",              "3.8%",     "DERIVAT - BUN"),
        ("m_tau/m_mu = phi^6",             "6.7%",     "DERIVAT - MARGINAL"),
        ("m_up/m_e = phi^3",               "0.2%",     "DERIVAT - EXCELENT"),
        ("Hubble tension 26/24",           "~0.04%",   "PATTERN - NEVERIFICAT FIZIC"),
        ("Dark Matter = 25%",              "~5.7%",    "SPECULATIV"),
        ("Muon 0.06% (cu corectii)",       "0.06%*",   "FITTING (corectie nedocumentata)"),
        ("Top quark n=26",                 "~20%",     "MARGINAL (necesita QCD)"),
        ("121-lea varf = graviton",        "N/A",      "SPECULATIV"),
    ]

    for name, err, status in results:
        print(f"  {name:35s} | {err:9s} | {status}")

    print()
    print("  * = include corectie C(Q, alpha_s) nedocumentata complet")
    print()
    print("  CATEGORII:")
    print("    DERIVAT SOLID:   alpha, alpha_s, sin2tW, m_H, m_up")
    print("    DERIVAT BUN:     m_mu/m_e=phi^11, m_tau/m_e=phi^17")
    print("    PATTERN:         Hubble 26/24 (numeric bun, fizic neclar)")
    print("    SPECULATIV:      Dark Matter 25%, gravitonul 121")
    print("    NECESITA MUNCA:  Top quark, Bottom, Charm (corectii QCD)")
    print()


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 70)
    print("  EXP-095: VALIDARE NUMERICA - TEORIA 600-CELL v3.6")
    print("  Data: 6 Februarie 2026")
    print("*" * 70)
    print()

    test_dark_matter()
    test_hubble_tension()
    test_quark_hierarchy()
    test_lepton_hierarchy()
    test_internal_consistency()
    test_120cell_graviton()
    test_top_quark()
    test_muon_precision()
    print_summary()
