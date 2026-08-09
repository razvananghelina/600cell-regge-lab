"""
EXP-193: DERIVATION OF THE ELECTROWEAK VEV FROM 600-CELL GEOMETRY
==================================================================

The electroweak vacuum expectation value (vev = 246.22 GeV) sets the
overall energy scale of the electroweak sector. It is currently NOT
derived in the 600-cell theory. This experiment systematically searches
for geometric expressions that reproduce vev.

Parts:
  A. Dimensional analysis (ratios with known masses)
  B. Volume / geometric approach (600-cell invariants)
  C. Spectral structure (eigenvalues of adjacency matrix)
  D. Fermi constant / weak coupling route
  E. Higgs quartic coupling
  F. Number theoretic approach
  G. From the unified mass formula
  H. Systematic phi-power scan
  I. Combined alpha-phi expressions
  J. Summary and best candidates

Author: Theory of Everything - 600-cell project
Date: 2026-02-10
"""

import math
import sys

# =============================================================================
# CONSTANTS
# =============================================================================
PHI = (1 + math.sqrt(5)) / 2
PHI_PRIME = 1 - PHI  # = -1/phi
PI = math.pi
ALPHA = 7.2973525693e-3
ALPHA_INV = 1.0 / ALPHA
ALPHA_S = 0.1179  # = 1/(2*phi^3) in theory
SIN2_TW = 0.23122
COS2_TW = 1 - SIN2_TW

# Masses in GeV
M_E = 0.51099895e-3      # electron
M_MU = 0.1056583755      # muon
M_TAU = 1.77686           # tau
M_W = 80.377              # W boson
M_Z = 91.1876             # Z boson
M_H = 125.25              # Higgs
M_PLANCK = 1.22089e19     # Planck mass
VEV_EXP = 246.22          # experimental vev

# Fermi constant
G_F = 1.1663788e-5        # GeV^-2

# 600-cell geometric quantities
N_VERT = 120
N_EDGE = 720
N_FACE = 1200
N_CELL = 600

# Eigenvalues of 600-cell adjacency matrix
# {12, 6*phi, 4*phi, 3, 0, -2, 4-4*phi, -3, 6-6*phi}
EIGENVALS = [
    12,
    6*PHI,
    4*PHI,
    3,
    0,
    -2,
    4 - 4*PHI,
    -3,
    6 - 6*PHI,
]

# Multiplicities: {1, 4, 9, 16, 25, 36, 9, 16, 4}
MULTS = [1, 4, 9, 16, 25, 36, 9, 16, 4]

# Eigenvalue names for readability
EVAL_NAMES = [
    "12", "6*phi", "4*phi", "3", "0", "-2", "4-4phi", "-3", "6-6phi"
]

def pct_err(predicted, experimental):
    """Percent error."""
    return abs(predicted - experimental) / abs(experimental) * 100

def log_phi(x):
    """Log base phi."""
    if x <= 0:
        return float('nan')
    return math.log(x) / math.log(PHI)


# =============================================================================
# PART A: DIMENSIONAL ANALYSIS
# =============================================================================
def part_a():
    print("=" * 72)
    print("PART A: DIMENSIONAL ANALYSIS - Ratios with known masses")
    print("=" * 72)
    print()

    ratios = [
        ("vev / m_e",    VEV_EXP / M_E),
        ("vev / m_mu",   VEV_EXP / M_MU),
        ("vev / m_tau",  VEV_EXP / M_TAU),
        ("vev / m_W",    VEV_EXP / M_W),
        ("vev / m_Z",    VEV_EXP / M_Z),
        ("vev / m_H",    VEV_EXP / M_H),
        ("vev / m_P",    VEV_EXP / M_PLANCK),
        ("m_W / m_e",    M_W / M_E),
        ("m_Z / m_e",    M_Z / M_E),
        ("m_H / m_e",    M_H / M_E),
    ]

    print("  %-18s  %18s  %12s  %s" % ("Ratio", "Value", "log_phi", "Nearest int/half"))
    print("  " + "-" * 70)
    for name, val in ratios:
        lp = log_phi(val)
        nearest_int = round(lp)
        nearest_half = round(lp * 2) / 2.0
        phi_int = PHI ** nearest_int
        phi_half = PHI ** nearest_half
        err_int = pct_err(phi_int, val)
        err_half = pct_err(phi_half, val)
        print("  %-18s  %18.6f  %12.4f  n=%5.1f (%.2f%%), n=%5.1f (%.2f%%)" %
              (name, val, lp, nearest_int, err_int, nearest_half, err_half))

    print()
    # Special: vev/m_e
    vev_me = VEV_EXP / M_E
    lp = log_phi(vev_me)
    print("  KEY: vev/m_e = %.2f, log_phi = %.6f" % (vev_me, lp))
    print("  Check 5a+6b form:")
    found = []
    for a in range(-10, 30):
        for b in range(-10, 30):
            n = 5*a + 6*b
            if abs(n - lp) < 0.5:
                pred = M_E * PHI**n
                err = pct_err(pred, VEV_EXP)
                if err < 10:
                    found.append((a, b, n, pred, err))
    found.sort(key=lambda x: x[4])
    for a, b, n, pred, err in found[:10]:
        print("    (a=%2d, b=%2d) -> n=%3d, phi^n*m_e = %.4f GeV, err = %.2f%%" %
              (a, b, n, pred, err))

    # Special: vev/m_W
    vev_mw = VEV_EXP / M_W
    lp_mw = log_phi(vev_mw)
    print()
    print("  KEY: vev/m_W = %.6f, log_phi = %.6f" % (vev_mw, lp_mw))
    # Check simple fractions and phi expressions
    tests = [
        ("2/cos(tW)", 2 / math.sqrt(COS2_TW)),
        ("2/sqrt(6/26)", 2 / math.sqrt(6.0/26)),
        ("2/(g'/2) where g' = e/cos", 2.0 * math.sqrt(COS2_TW) / SIN2_TW),
        ("phi^(3/2)", PHI**1.5),
        ("sqrt(pi)", math.sqrt(PI)),
        ("phi + 1/phi", PHI + 1/PHI),
        ("phi^2", PHI**2),
        ("2*phi/sqrt(3)", 2*PHI/math.sqrt(3)),
        ("1/sin(tW)^2", 1/SIN2_TW),
    ]

    # Actually vev = 2*m_W/g_w and g_w = e/sin(tW) so
    # vev = 2*m_W*sin(tW)/e = 2*m_W*sin(tW)/(sqrt(4*pi*alpha))
    # Let's check this numerically
    g_w = math.sqrt(4 * PI * ALPHA) / math.sqrt(SIN2_TW)
    vev_check = 2 * M_W / g_w
    print("  Check SM formula: vev = 2*m_W/g_w = %.4f GeV (err %.4f%%)" %
          (vev_check, pct_err(vev_check, VEV_EXP)))

    print()
    print("  vev/m_W simple expressions:")
    for name, val in tests:
        err = pct_err(val, vev_mw)
        print("    %-30s = %10.6f  err = %8.4f%%" % (name, val, err))

    # Key relation: vev = 2*m_W/g and m_W = vev*g/2
    # g = e/sin(tW), so vev/m_W = 2*sin(tW)/e = 2*sin(tW)/sqrt(4*pi*alpha)
    # This is a SM identity, not a derivation.
    print()
    print("  SM identity: vev/m_W = 2*sin(tW)/e = 2*sqrt(sin2_tW)/sqrt(4*pi*alpha)")
    sm_ratio = 2 * math.sqrt(SIN2_TW) / math.sqrt(4 * PI * ALPHA)
    print("    = %.6f vs actual ratio %.6f (err %.4f%%)" %
          (sm_ratio, vev_mw, pct_err(sm_ratio, vev_mw)))
    print("  With sin^2(tW) = 6/26:")
    sm_ratio_600 = 2 * math.sqrt(6.0/26) / math.sqrt(4 * PI * ALPHA)
    print("    = %.6f (err vs actual ratio: %.4f%%)" %
          (sm_ratio_600, pct_err(sm_ratio_600, vev_mw)))


# =============================================================================
# PART B: VOLUME / GEOMETRIC APPROACH
# =============================================================================
def part_b():
    print()
    print("=" * 72)
    print("PART B: VOLUME AND GEOMETRIC APPROACH")
    print("=" * 72)
    print()

    # Volume of unit 600-cell (inscribed in S^3 of radius 1)
    # The 600-cell has 600 regular tetrahedra
    # Edge length of 600-cell inscribed in unit S^3: a = 1/phi
    a_edge = 1.0 / PHI
    print("  600-cell edge length (unit S^3): a = 1/phi = %.6f" % a_edge)

    # Volume of regular tetrahedron with edge a: V = a^3/(6*sqrt(2))
    V_tet = a_edge**3 / (6 * math.sqrt(2))
    V_600 = 600 * V_tet  # This is approximate for flat tetrahedra
    print("  Volume of one flat tet (edge 1/phi): %.6f" % V_tet)
    print("  600 * V_tet (flat approx): %.6f" % V_600)

    # Volume of S^3 with radius R: V = 2*pi^2*R^3
    V_S3 = 2 * PI**2
    print("  Volume of unit S^3: 2*pi^2 = %.6f" % V_S3)

    # Ratio
    print("  600*V_tet / V_S3 = %.6f" % (V_600 / V_S3))

    # The actual 3-volume of 600-cell in R^4 with circumradius 1:
    # V = 50*phi^2 / sqrt(2) (known result)
    # Actually: 600-cell 3-volume = (50/sqrt(2)) * phi^2 for circumradius 1
    V_600_exact = 50 * PHI**2 / math.sqrt(2)
    print("  Exact 600-cell 3-content: 50*phi^2/sqrt(2) = %.6f" % V_600_exact)

    print()
    print("  Testing vev = m_P * f(geometry):")
    # vev/m_P ~ 2.02e-17
    vev_mp = VEV_EXP / M_PLANCK
    print("  vev/m_P = %.6e" % vev_mp)

    # From known: m_e/m_P = alpha^(4*phi^2)
    # And: vev/m_e ~ phi^27 (from part A)
    # So: vev/m_P = vev/m_e * m_e/m_P ~ phi^27 * alpha^(4*phi^2)
    lp_vev_me = log_phi(VEV_EXP / M_E)
    print("  log_phi(vev/m_e) = %.6f" % lp_vev_me)
    print("  So vev ~ m_e * phi^%.2f" % lp_vev_me)
    print("  And m_e = m_P * alpha^(4*phi^2)")
    print("  => vev = m_P * alpha^(4*phi^2) * phi^%.2f" % lp_vev_me)
    val = M_PLANCK * ALPHA**(4*PHI**2) * PHI**round(lp_vev_me)
    print("  vev = m_P * alpha^(4*phi^2) * phi^%d = %.4f GeV (err %.2f%%)" %
          (round(lp_vev_me), val, pct_err(val, VEV_EXP)))

    # Test geometric combinations
    print()
    print("  Testing vev = m_e * (geometric quantity)^k * phi^j:")
    geom_quants = {
        "120 (vertices)": 120,
        "720 (edges)": 720,
        "1200 (faces)": 1200,
        "600 (cells)": 600,
        "120*720": 120*720,
        "sqrt(120)": math.sqrt(120),
        "sqrt(720)": math.sqrt(720),
        "120^(1/3)": 120**(1.0/3),
        "720/120": 6.0,
        "1200/120": 10.0,
    }

    target = VEV_EXP / M_E
    for name, gval in geom_quants.items():
        if gval > 0:
            ratio = target / gval
            lp = log_phi(ratio)
            nearest = round(lp)
            pred = M_E * gval * PHI**nearest
            err = pct_err(pred, VEV_EXP)
            print("    vev = m_e * %-16s * phi^%3d = %10.4f GeV  (err %7.3f%%)" %
                  (name, nearest, pred, err))


# =============================================================================
# PART C: SPECTRAL STRUCTURE
# =============================================================================
def part_c():
    print()
    print("=" * 72)
    print("PART C: SPECTRAL STRUCTURE OF 600-CELL")
    print("=" * 72)
    print()

    print("  Eigenvalues and multiplicities:")
    print("  %-12s  %12s  %6s  %12s" % ("Name", "Value", "Mult", "lam^2 * mult"))
    print("  " + "-" * 50)
    total_trace2 = 0
    for i in range(9):
        val = EIGENVALS[i]
        m = MULTS[i]
        contrib = val**2 * m
        total_trace2 += contrib
        print("  %-12s  %12.6f  %6d  %12.4f" % (EVAL_NAMES[i], val, m, contrib))

    print("  " + "-" * 50)
    print("  Tr(A^2) = sum(lam^2 * mult) = %.4f" % total_trace2)
    print("  This should equal 2*|edges| = 2*720 = 1440? No, = sum of degree^2")
    # For regular graph with degree d: Tr(A^2) = N*d where N=vertices, d=degree
    # Actually Tr(A^2) = 2*|edges| = 2*720 = 1440? No.
    # Tr(A^2)_ij = sum of (A^2)_ii = number of walks of length 2 from i to i = degree(i)
    # For a regular graph with degree 12: Tr(A^2) = 120 * 12 = 1440
    print("  Check: N*degree = 120*12 = %d" % (120 * 12))

    # Product of nonzero eigenvalues
    nonzero_evals = [(EIGENVALS[i], MULTS[i]) for i in range(9) if abs(EIGENVALS[i]) > 1e-10]
    log_prod = 0
    for val, m in nonzero_evals:
        log_prod += m * math.log(abs(val))
    prod_nonzero = math.exp(log_prod)
    print()
    print("  Product of |nonzero eigenvalues|^mult = exp(sum m*ln|lam|)")
    print("    = exp(%.4f) = %.6e" % (log_prod, prod_nonzero))

    # Spectral zeta: Z(s) = sum mult_i * |lam_i|^(-s) for nonzero
    print()
    print("  Spectral zeta Z(s) = sum mult_i * |lam_i|^(-s) [nonzero only]:")
    for s in [1, 2, 3, 0.5, -1, -2]:
        Z = 0
        for val, m in nonzero_evals:
            Z += m * abs(val)**(-s)
        print("    Z(%.1f) = %.6f" % (s, Z))
        # Check if any ratio gives vev
        for scale_name, scale in [("m_e", M_E), ("m_W", M_W)]:
            pred = scale * Z
            if pred > 1 and pred < 1e6:
                err = pct_err(pred, VEV_EXP)
                print("      %s * Z = %.4f GeV (err %.2f%%)" % (scale_name, pred, err))

    # Sum of eigenvalues (trace) = 0 for non-trivial
    trace = sum(EIGENVALS[i] * MULTS[i] for i in range(9))
    print()
    print("  Trace (sum lam*mult) = %.4f (should be degree*N if including self? No, = sum of eigenvalues)" % trace)
    print("  For a graph, Tr(A) = 0 always. Check: trace should = 12*1 + ... = %.4f" % trace)

    # Characteristic polynomial value at special points
    print()
    print("  Testing spectral quantities as vev expressions:")

    # Sum of positive eigenvalues * their multiplicities
    pos_sum = sum(EIGENVALS[i] * MULTS[i] for i in range(9) if EIGENVALS[i] > 0)
    neg_sum = sum(abs(EIGENVALS[i]) * MULTS[i] for i in range(9) if EIGENVALS[i] < 0)
    print("    Sum of (positive lam * mult) = %.4f" % pos_sum)
    print("    Sum of (|negative lam| * mult) = %.4f" % neg_sum)
    print("    Ratio pos/neg = %.6f" % (pos_sum / neg_sum))

    # Phi content
    # 6*phi, 4*phi, 4-4*phi, 6-6*phi are the phi-dependent eigenvalues
    phi_part = 6*PHI*4 + 4*PHI*9  # mult-weighted phi part from positive phi eigenvalues
    print("    Phi-weighted positive: 6*phi*4 + 4*phi*9 = %.4f" % phi_part)

    # Total spectral weight (Frobenius norm of adjacency matrix)
    frob = math.sqrt(total_trace2)
    print("    Frobenius norm = sqrt(Tr(A^2)) = %.6f" % frob)
    print("    = sqrt(1440) = %.6f" % math.sqrt(1440))

    # Try: vev = m_e * frob * phi^k
    for k in range(20, 30):
        pred = M_E * frob * PHI**k
        err = pct_err(pred, VEV_EXP)
        if err < 5:
            print("    vev = m_e * sqrt(1440) * phi^%d = %.4f GeV (err %.2f%%)" % (k, pred, err))

    # Sum of squared eigenvalues (without multiplicity)
    sum_sq_raw = sum(e**2 for e in EIGENVALS)
    print("    Sum of lam^2 (no mult) = %.4f" % sum_sq_raw)

    # Spectral gap
    gap = min(abs(e) for e in EIGENVALS if abs(e) > 1e-10)
    print("    Spectral gap (smallest |nonzero lam|) = %.6f" % gap)


# =============================================================================
# PART D: FERMI CONSTANT / WEAK COUPLING
# =============================================================================
def part_d():
    print()
    print("=" * 72)
    print("PART D: FERMI CONSTANT AND WEAK COUPLING ROUTE")
    print("=" * 72)
    print()

    # Standard derivation: vev = 1/sqrt(sqrt(2)*G_F)
    vev_from_gf = 1.0 / math.sqrt(math.sqrt(2) * G_F)
    print("  vev = 1/sqrt(sqrt(2)*G_F) = %.4f GeV (err %.4f%%)" %
          (vev_from_gf, pct_err(vev_from_gf, VEV_EXP)))

    # Can we express G_F in terms of alpha, phi, m_e?
    # G_F ~ 1.166e-5 GeV^-2
    # G_F * m_e^2 = 1.166e-5 * (0.511e-3)^2 = 3.045e-12 (dimensionless * GeV^0)
    gf_me2 = G_F * M_E**2
    print()
    print("  G_F * m_e^2 = %.6e (dimensionless)" % gf_me2)
    print("  log_phi(G_F * m_e^2) = %.6f" % log_phi(gf_me2))

    # Test: G_F = alpha^k / m_e^2
    la = math.log(ALPHA)
    target_gf = math.log(G_F * M_E**2)
    k_alpha = target_gf / la
    print("  If G_F = alpha^k / m_e^2, then k = ln(G_F*m_e^2)/ln(alpha) = %.6f" % k_alpha)
    pred_gf = ALPHA**k_alpha / M_E**2
    print("    Check: alpha^%.4f / m_e^2 = %.6e (exact by construction)" % (k_alpha, pred_gf))

    # Can k be a simple expression in phi?
    print("    k = %.6f" % k_alpha)
    print("    phi = %.6f, phi^2 = %.6f, 2*phi = %.6f" % (PHI, PHI**2, 2*PHI))
    print("    k/phi = %.6f, k/phi^2 = %.6f" % (k_alpha/PHI, k_alpha/PHI**2))
    # Check some simple forms
    candidates_k = [
        ("2*phi + 1", 2*PHI + 1),
        ("phi^2 + 1", PHI**2 + 1),
        ("2*phi^2 - 1", 2*PHI**2 - 1),
        ("phi^3/2", PHI**3/2),
        ("5/2", 2.5),
        ("phi^2", PHI**2),
        ("3", 3.0),
        ("2+phi", 2+PHI),
        ("phi^2 + phi/2", PHI**2 + PHI/2),
        ("1 + 2/phi", 1 + 2/PHI),
        ("phi + 1 + 1/phi", PHI + 1 + 1/PHI),
    ]
    print()
    print("    Matching k = %.6f to simple phi expressions:" % k_alpha)
    for name, val in candidates_k:
        err = pct_err(val, k_alpha)
        pred_gf2 = ALPHA**val / M_E**2
        pred_vev = 1.0 / math.sqrt(math.sqrt(2) * pred_gf2)
        err_vev = pct_err(pred_vev, VEV_EXP)
        if err < 20:
            print("      %-22s = %.6f (err on k: %.2f%%, vev: %.4f GeV, err: %.2f%%)" %
                  (name, val, err, pred_vev, err_vev))

    # Key approach: use already-derived relations
    # m_W = vev * g / 2 and g = e / sin(tW)
    # e = sqrt(4*pi*alpha)
    # So m_W = vev * sqrt(4*pi*alpha) / (2*sin(tW))
    # => vev = 2 * m_W * sin(tW) / sqrt(4*pi*alpha)
    # With sin^2(tW) = 6/26:
    print()
    print("  KEY ROUTE: vev = 2 * m_W * sin(tW) / e")
    print("  Using sin^2(tW) = 6/26 (DERIVED) and alpha from 600-cell equation:")
    sin_tW = math.sqrt(6.0/26)
    e_em = math.sqrt(4 * PI * ALPHA)
    vev_route1 = 2 * M_W * sin_tW / e_em
    print("    vev = 2 * %.3f * %.6f / %.6f = %.4f GeV (err %.4f%%)" %
          (M_W, sin_tW, e_em, vev_route1, pct_err(vev_route1, VEV_EXP)))
    print("    This uses m_W as INPUT => not a full derivation.")

    # Can we derive m_W from m_e?
    # m_W/m_e ~ 80377/0.511 = 157293
    mw_me = M_W / M_E
    lp_mw = log_phi(mw_me)
    print()
    print("  m_W / m_e = %.2f, log_phi = %.6f" % (mw_me, lp_mw))
    # 5a+6b candidates
    print("  Checking 5a+6b form for m_W/m_e:")
    found = []
    for a in range(-5, 20):
        for b in range(-5, 20):
            n = 5*a + 6*b
            if abs(n - lp_mw) < 0.5:
                pred = M_E * PHI**n
                err = pct_err(pred, M_W)
                if err < 15:
                    found.append((a, b, n, pred, err))
    found.sort(key=lambda x: x[4])
    for a, b, n, pred, err in found[:5]:
        print("    (a=%2d, b=%2d) -> n=%3d, m_e*phi^n = %.4f GeV (err %.2f%%)" %
              (a, b, n, pred, err))

    # Now try the full chain: m_W from m_e, then vev from m_W
    print()
    print("  FULL CHAIN ATTEMPT: m_e -> m_W -> vev")
    if found:
        a_best, b_best, n_best, mw_pred, mw_err = found[0]
        vev_chain = 2 * mw_pred * sin_tW / e_em
        print("    m_W = m_e * phi^%d = %.4f GeV (err %.2f%%)" %
              (n_best, mw_pred, mw_err))
        print("    vev = 2*m_W*sqrt(6/26)/sqrt(4*pi*alpha) = %.4f GeV (err %.2f%%)" %
              (vev_chain, pct_err(vev_chain, VEV_EXP)))

    # Alternative: express vev directly as alpha^a * m_e^b * phi^c
    # vev = m_e * phi^N * alpha^M for some N, M
    # log(vev/m_e) = N*log(phi) + M*log(alpha)
    # We need two conditions to solve for N and M.
    # But we can scan:
    print()
    print("  Scanning vev = m_e * alpha^M * phi^N:")
    log_target = math.log(VEV_EXP / M_E)
    log_a = math.log(ALPHA)
    log_p = math.log(PHI)
    best_results = []
    for M_num in range(-30, 31):  # M in units of 1/10
        for M_den in [1, 2, 3, 4, 5, 6, 10]:
            M = M_num / float(M_den)
            if abs(M) > 5:
                continue
            # N = (log_target - M*log_a) / log_p
            N = (log_target - M * log_a) / log_p
            N_round = round(N)
            pred = M_E * ALPHA**M * PHI**N_round
            err = pct_err(pred, VEV_EXP)
            if err < 2:
                best_results.append((M, N_round, pred, err))

    best_results.sort(key=lambda x: x[3])
    if best_results:
        print("    Best results (err < 2%%):")
        seen = set()
        count = 0
        for M, N, pred, err in best_results:
            key = (round(M*10), N)
            if key not in seen:
                seen.add(key)
                print("      alpha^(%.3f) * phi^%d * m_e = %.4f GeV (err %.4f%%)" %
                      (M, N, pred, err))
                count += 1
                if count >= 10:
                    break
    else:
        print("    No results with err < 2%%. Expanding search...")
        # Relax to 5%
        for M_num in range(-50, 51):
            for M_den in [1, 2, 3, 4, 5, 6, 10]:
                M = M_num / float(M_den)
                if abs(M) > 5:
                    continue
                N = (log_target - M * log_a) / log_p
                N_round = round(N)
                pred = M_E * ALPHA**M * PHI**N_round
                err = pct_err(pred, VEV_EXP)
                if err < 5:
                    best_results.append((M, N_round, pred, err))
        best_results.sort(key=lambda x: x[3])
        seen = set()
        count = 0
        for M, N, pred, err in best_results[:10]:
            key = (round(M*10), N)
            if key not in seen:
                seen.add(key)
                print("      alpha^(%.3f) * phi^%d * m_e = %.4f GeV (err %.4f%%)" %
                      (M, N, pred, err))
                count += 1


# =============================================================================
# PART E: HIGGS QUARTIC COUPLING
# =============================================================================
def part_e():
    print()
    print("=" * 72)
    print("PART E: HIGGS QUARTIC COUPLING LAMBDA")
    print("=" * 72)
    print()

    lam = M_H**2 / (2 * VEV_EXP**2)
    print("  lambda = m_H^2 / (2*vev^2) = %.6f" % lam)
    print("  alpha_s = %.6f" % ALPHA_S)
    print("  lambda / alpha_s = %.6f" % (lam / ALPHA_S))
    print("  Difference: %.2f%%" % pct_err(lam, ALPHA_S))
    print()

    # Test various phi expressions
    tests = [
        ("alpha_s * phi/2", ALPHA_S * PHI / 2),
        ("alpha_s * phi^2/4", ALPHA_S * PHI**2 / 4),
        ("alpha_s / (2*cos(pi/5))", ALPHA_S / (2 * math.cos(PI/5))),
        ("1/(2*phi^3)", 1.0 / (2 * PHI**3)),  # = alpha_s_theory
        ("1/(2*phi^3) * phi", PHI / (2 * PHI**3)),
        ("phi^(-4)", PHI**(-4)),
        ("phi^(-4) * (1 + alpha)", PHI**(-4) * (1 + ALPHA)),
        ("alpha_s * (1 + alpha_s)", ALPHA_S * (1 + ALPHA_S)),
        ("ALPHA * 2*phi^3", ALPHA * 2 * PHI**3),  # = alpha * alpha_s/alpha_theory = alpha * 10*phi
        ("ALPHA * 10 * phi * phi", ALPHA * 10 * PHI * PHI),
        ("pi * alpha / (2*sin^2_tW)", PI * ALPHA / (2 * SIN2_TW)),
        ("alpha / (2*alpha_s)", ALPHA / (2 * ALPHA_S)),
        ("phi^2/(4*pi)", PHI**2 / (4*PI)),
        ("1/8", 0.125),
        ("phi/12", PHI/12),
        ("phi/(4*pi)", PHI/(4*PI)),
        ("sqrt(alpha_s/pi)", math.sqrt(ALPHA_S/PI)),
    ]

    print("  Testing lambda = f(alpha, phi, alpha_s, sin^2_tW):")
    print("  %-35s  %10s  %10s" % ("Expression", "Value", "Error %"))
    print("  " + "-" * 60)
    for name, val in tests:
        err = pct_err(val, lam)
        marker = " <---" if err < 5 else ""
        print("  %-35s  %10.6f  %10.4f%%%s" % (name, val, err, marker))

    # Deeper: m_H = m_W*(phi - 8*alpha) is derived.
    # lambda = m_H^2/(2*vev^2) = [m_W*(phi-8*alpha)]^2 / (2*vev^2)
    # = m_W^2*(phi-8*alpha)^2 / (2*vev^2)
    # But vev = 2*m_W/g, so vev^2 = 4*m_W^2/g^2
    # lambda = m_W^2*(phi-8a)^2 / (2*4*m_W^2/g^2) = g^2*(phi-8a)^2 / 8
    # g = e/sin_tW = sqrt(4*pi*alpha)/sin_tW
    # g^2 = 4*pi*alpha/sin^2_tW
    # lambda = 4*pi*alpha*(phi-8a)^2 / (8*sin^2_tW)
    #        = pi*alpha*(phi-8a)^2 / (2*sin^2_tW)
    phi_8a = PHI - 8*ALPHA
    lam_derived = PI * ALPHA * phi_8a**2 / (2 * SIN2_TW)
    print()
    print("  DERIVED FORMULA: lambda = pi*alpha*(phi-8*alpha)^2 / (2*sin^2_tW)")
    print("    = pi * %.6f * %.6f / (2 * %.5f)" % (ALPHA, phi_8a**2, SIN2_TW))
    print("    = %.6f (err %.4f%%)" % (lam_derived, pct_err(lam_derived, lam)))

    # With sin^2_tW = 6/26:
    lam_600 = PI * ALPHA * phi_8a**2 / (2 * 6.0/26)
    print("    With sin^2_tW = 6/26: lambda = %.6f (err %.4f%%)" % (lam_600, pct_err(lam_600, lam)))
    print()
    print("  NOTE: This is a SM identity (not a new derivation), since")
    print("  lambda = g^2*(phi-8*alpha)^2/8 follows from m_H = m_W*(phi-8*alpha).")


# =============================================================================
# PART F: NUMBER THEORETIC APPROACH
# =============================================================================
def part_f():
    print()
    print("=" * 72)
    print("PART F: NUMBER THEORETIC APPROACH")
    print("=" * 72)
    print()

    print("  vev = %.2f GeV" % VEV_EXP)
    print("  246 = 2 * 3 * 41")
    print("  m_e = %.6f MeV" % (M_E * 1000))
    vev_me = VEV_EXP / M_E
    print("  vev/m_e = %.2f" % vev_me)

    # 600-cell number combinations
    print()
    print("  Testing vev = m_e * N_geom * phi^k:")
    geom_nums = {
        "120": 120,
        "720": 720,
        "1200": 1200,
        "600": 600,
        "120^2": 14400,
        "720/6": 120,
        "12 (degree)": 12,
        "12^2": 144,
        "12^3": 1728,
        "24 (24-cell verts)": 24,
        "48 (24-cell edges)": 48,
        "96 (24-cell faces)": 96,
        "8 (24-cell cells)": 8,
        "240 (E8 roots)": 240,
        "120*12": 1440,
        "120*12/2": 720,
        "2*pi^2 (V_S3)": 2*PI**2,
    }

    results = []
    for name, gval in geom_nums.items():
        if gval > 0:
            ratio = vev_me / gval
            lp = log_phi(ratio)
            # Try nearest integer and half-integer
            for n in [round(lp), round(lp*2)/2.0]:
                pred = M_E * gval * PHI**n
                err = pct_err(pred, VEV_EXP)
                if err < 5:
                    results.append((name, gval, n, pred, err))

    results.sort(key=lambda x: x[4])
    for name, gval, n, pred, err in results[:15]:
        n_str = "%d" % n if n == int(n) else "%.1f" % n
        print("    m_e * %-20s * phi^%-5s = %10.4f GeV  (err %7.3f%%)" %
              (name, n_str, pred, err))

    # Test: is vev/m_e close to a product of 600-cell numbers times phi power?
    print()
    print("  Fibonacci/Lucas connections:")
    # Fibonacci: 1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,...
    # Lucas:     2,1,3,4,7,11,18,29,47,76,123,199,322,521,843,...
    fibs = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
            1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393,
            196418, 317811, 514229, 832040]
    lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521, 843,
             1364, 2207, 3571, 5778, 9349, 15127, 24476, 39603, 64079, 103682,
             167761, 271443, 439204, 710647]

    # phi^n = F(n)*phi + F(n-1) for n>=1
    # So vev/m_e ~ phi^27.15
    # phi^27 = F(27)*phi + F(26) = 196418*phi + 121393
    phi27 = 196418 * PHI + 121393
    print("  phi^27 = F(27)*phi + F(26) = 196418*phi + 121393 = %.2f" % phi27)
    print("  m_e * phi^27 = %.4f GeV" % (M_E * phi27))
    print("  vev/m_e = %.2f, phi^27 = %.2f (ratio %.6f)" %
          (vev_me, phi27, vev_me / phi27))

    # What if vev/m_e = F(n) or L(n)?
    print()
    for n, fn in enumerate(fibs):
        pred = M_E * fn
        if 100 < pred < 500:
            err = pct_err(pred, VEV_EXP)
            print("    m_e * F(%d) = m_e * %d = %.4f GeV (err %.2f%%)" % (n, fn, pred, err))
    for n, ln in enumerate(lucas):
        pred = M_E * ln
        if 100 < pred < 500:
            err = pct_err(pred, VEV_EXP)
            print("    m_e * L(%d) = m_e * %d = %.4f GeV (err %.2f%%)" % (n, ln, pred, err))


# =============================================================================
# PART G: FROM THE UNIFIED MASS FORMULA
# =============================================================================
def part_g():
    print()
    print("=" * 72)
    print("PART G: FROM THE UNIFIED MASS FORMULA")
    print("=" * 72)
    print()

    # m_f = m_e * phi^(5a + 6b) at leading order
    # vev = m_e * phi^N where N = log_phi(vev/m_e)
    N_vev = log_phi(VEV_EXP / M_E)
    print("  log_phi(vev/m_e) = %.6f" % N_vev)
    print("  Nearest integer: %d" % round(N_vev))
    print()

    # Find all (a,b) with 5a + 6b = round(N_vev)
    N_int = round(N_vev)
    print("  Solutions to 5a + 6b = %d:" % N_int)
    solutions = []
    for a in range(-20, 30):
        for b in range(-20, 30):
            if 5*a + 6*b == N_int:
                solutions.append((a, b))
    for a, b in solutions[:15]:
        # Compute norm
        norm = a**2 + a*b - b**2
        print("    (a=%3d, b=%3d)  norm = a^2+ab-b^2 = %d" % (a, b, norm))

    pred_27 = M_E * PHI**27
    err_27 = pct_err(pred_27, VEV_EXP)
    print()
    print("  m_e * phi^27 = %.4f GeV (err %.4f%%)" % (pred_27, err_27))

    # What about phi^(5*a + 6*b + delta) where delta is a correction?
    # vev/m_e = phi^27.15... so delta ~ 0.15
    delta = N_vev - 27
    print("  Fractional part: delta = %.6f" % delta)

    # Can delta be expressed simply?
    print("  delta/alpha = %.4f" % (delta / ALPHA))
    print("  delta*phi = %.6f" % (delta * PHI))
    print("  delta*pi = %.6f" % (delta * PI))
    print("  delta*137 = %.4f" % (delta * 137))
    print("  1/delta = %.4f" % (1/delta if abs(delta) > 1e-10 else float('inf')))

    # Check: is delta ~ alpha*k for simple k?
    for k_num in range(1, 50):
        for k_den in [1, 2, 3, 4, 5, 6, PI]:
            candidate = ALPHA * k_num / k_den
            if abs(candidate - delta) / abs(delta) < 0.05:
                # compute vev
                N_test = 27 + candidate
                pred = M_E * PHI**N_test
                err = pct_err(pred, VEV_EXP)
                if isinstance(k_den, float):
                    print("    delta ~ alpha * %d / pi = %.6f => vev = %.4f GeV (err %.4f%%)" %
                          (k_num, candidate, pred, err))
                else:
                    print("    delta ~ alpha * %d/%d = %.6f => vev = %.4f GeV (err %.4f%%)" %
                          (k_num, k_den, candidate, pred, err))

    # What if vev corresponds to a BOSON slot in the mass formula?
    # Compare with W, Z, H positions:
    print()
    print("  Comparison with boson positions:")
    for name, mass in [("W", M_W), ("Z", M_Z), ("H", M_H), ("vev", VEV_EXP)]:
        lp = log_phi(mass / M_E)
        print("    %-4s: m/m_e = %.2f, log_phi = %.4f, 5a+6b candidates:" % (name, mass/M_E, lp))
        n_near = round(lp)
        for a in range(-5, 15):
            for b in range(-5, 15):
                if 5*a + 6*b == n_near:
                    norm = a**2 + a*b - b**2
                    pred = M_E * PHI**n_near
                    err = pct_err(pred, mass)
                    if abs(a) < 8 and abs(b) < 8:
                        print("      (a=%2d, b=%2d) norm=%3d  pred=%.2f GeV  err=%.2f%%" %
                              (a, b, norm, pred, err))
                        break


# =============================================================================
# PART H: SYSTEMATIC PHI-POWER SCAN (fine-grained)
# =============================================================================
def part_h():
    print()
    print("=" * 72)
    print("PART H: SYSTEMATIC SCAN - vev = m_e * phi^(p/q)")
    print("=" * 72)
    print()

    target = VEV_EXP / M_E
    lp = log_phi(target)
    print("  Exact log_phi(vev/m_e) = %.10f" % lp)
    print()

    results = []
    # Scan p/q for q up to 12
    for q in range(1, 13):
        for p in range(20*q, 35*q):
            frac = p / float(q)
            pred = PHI**frac
            err = pct_err(pred, target)
            if err < 0.5:
                results.append((p, q, frac, M_E * pred, err))

    results.sort(key=lambda x: x[4])
    print("  Best phi^(p/q) approximations (err < 0.5%%):")
    print("  %-8s  %10s  %12s  %10s" % ("p/q", "Value", "vev (GeV)", "Error %"))
    print("  " + "-" * 50)
    seen_frac = set()
    for p, q, frac, vev_pred, err in results[:20]:
        # Reduce fraction
        from math import gcd
        g = gcd(p, q)
        p_red, q_red = p // g, q // g
        key = (p_red, q_red)
        if key not in seen_frac:
            seen_frac.add(key)
            if q_red == 1:
                label = "%d" % p_red
            else:
                label = "%d/%d" % (p_red, q_red)
            print("  %-8s  %10.2f  %12.4f  %10.4f%%" % (label, frac, vev_pred, err))


# =============================================================================
# PART I: COMBINED ALPHA-PHI EXPRESSIONS
# =============================================================================
def part_i():
    print()
    print("=" * 72)
    print("PART I: COMBINED ALPHA-PHI EXPRESSIONS FOR VEV")
    print("=" * 72)
    print()

    # The most promising approach: combine known derived quantities
    # We know: alpha (DERIVED), sin^2(tW) = 6/26 (DERIVED), m_e/m_P = alpha^(4*phi^2) (DERIVED)
    # We need one more relation to get vev.

    # Approach 1: vev from the Planck mass
    # vev/m_P = alpha^k. What is k?
    k_vev = math.log(VEV_EXP / M_PLANCK) / math.log(ALPHA)
    print("  vev/m_P = %.6e" % (VEV_EXP / M_PLANCK))
    print("  If vev = m_P * alpha^k, then k = %.6f" % k_vev)
    print("  Compare with m_e/m_P: k_e = 4*phi^2 = %.6f" % (4*PHI**2))
    print("  Difference: k_vev - k_e = %.6f" % (k_vev - 4*PHI**2))
    print("  Ratio k_vev/k_e = %.6f" % (k_vev / (4*PHI**2)))
    print()

    # Try: vev = m_P * alpha^(a + b*phi + c*phi^2)
    print("  Scanning vev = m_P * alpha^(a + b*phi) for integer a,b:")
    best = []
    log_target = math.log(VEV_EXP / M_PLANCK) / math.log(ALPHA)
    for a in range(-10, 20):
        for b in range(-10, 20):
            k_test = a + b * PHI
            pred_vev = M_PLANCK * ALPHA**k_test
            if pred_vev > 0 and pred_vev < 1e6:
                err = pct_err(pred_vev, VEV_EXP)
                if err < 20:
                    best.append((a, b, k_test, pred_vev, err))
    best.sort(key=lambda x: x[4])
    for a, b, k, pred, err in best[:10]:
        print("    alpha^(%d + %d*phi) = alpha^(%.4f), vev = %.4f GeV (err %.4f%%)" %
              (a, b, k, pred, err))

    # Try a + b*phi^2
    print()
    print("  Scanning vev = m_P * alpha^(a + b*phi^2) for integer a,b:")
    best2 = []
    for a in range(-10, 20):
        for b in range(-10, 20):
            k_test = a + b * PHI**2
            pred_vev = M_PLANCK * ALPHA**k_test
            if pred_vev > 0 and pred_vev < 1e6:
                err = pct_err(pred_vev, VEV_EXP)
                if err < 20:
                    best2.append((a, b, k_test, pred_vev, err))
    best2.sort(key=lambda x: x[4])
    for a, b, k, pred, err in best2[:10]:
        print("    alpha^(%d + %d*phi^2) = alpha^(%.4f), vev = %.4f GeV (err %.4f%%)" %
              (a, b, k, pred, err))

    # Key test: vev = m_P * alpha^(4*phi^2) * phi^27
    # = m_e * phi^27 (since m_e = m_P * alpha^(4*phi^2))
    pred_me_phi27 = M_PLANCK * ALPHA**(4*PHI**2) * PHI**27
    err_27 = pct_err(pred_me_phi27, VEV_EXP)
    print()
    print("  KEY: vev = m_P * alpha^(4*phi^2) * phi^27 = m_e * phi^27")
    print("    = %.4f GeV (err %.4f%%)" % (pred_me_phi27, err_27))

    # With the fractional correction:
    exact_N = log_phi(VEV_EXP / M_E)
    pred_exact = M_PLANCK * ALPHA**(4*PHI**2) * PHI**exact_N
    print("    Exact exponent: phi^(%.6f), vev = %.4f GeV" % (exact_N, pred_exact))

    # Approach 2: from G_F directly
    # G_F = pi*alpha / (sqrt(2)*m_W^2*sin^2_tW)  -- tree level
    # If m_W = m_e * phi^N_W, then G_F can be expressed purely
    print()
    print("  Approach: G_F from alpha and m_W")
    gf_tree = PI * ALPHA / (math.sqrt(2) * M_W**2 * SIN2_TW)
    print("    G_F(tree) = pi*alpha/(sqrt(2)*m_W^2*sin^2_tW) = %.6e GeV^-2" % gf_tree)
    print("    G_F(exp) = %.6e GeV^-2" % G_F)
    print("    Ratio = %.6f (off by %.2f%%)" % (gf_tree/G_F, pct_err(gf_tree, G_F)))
    print("    NOTE: Tree-level G_F differs because G_F = g^2/(4*sqrt(2)*m_W^2)")
    gf_correct = 1.0 / (math.sqrt(2) * VEV_EXP**2)
    print("    Correct: G_F = 1/(sqrt(2)*vev^2) = %.6e" % gf_correct)

    # Approach 3: Look for the simplest expression
    print()
    print("  SIMPLEST EXPRESSIONS (combining all known derived quantities):")
    print()

    # vev in terms of alpha, phi, m_e or m_P
    expressions = [
        # (name, value)
        ("m_e * phi^27", M_E * PHI**27),
        ("m_e * phi^(27+alpha*20)", M_E * PHI**(27 + ALPHA*20)),
        ("m_e * phi^(27+alpha*21)", M_E * PHI**(27 + ALPHA*21)),
        ("m_e * phi^(27+alpha*22)", M_E * PHI**(27 + ALPHA*22)),
        ("m_e * phi^(5*5+6*0+alpha*20)", M_E * PHI**(25 + ALPHA*20)),
        ("m_e * phi^(5*3+6*2)", M_E * PHI**(5*3+6*2)),  # n=27
        ("m_P * alpha^(4*phi^2+1/phi)", M_PLANCK * ALPHA**(4*PHI**2 + 1/PHI)),
        ("m_P * alpha^(4*phi^2) * phi^27", M_PLANCK * ALPHA**(4*PHI**2) * PHI**27),
        ("m_W / sin(tW)^2 * sqrt(pi*alpha)", M_W / SIN2_TW * math.sqrt(PI * ALPHA)),
        ("2*m_W*sqrt(26/6)/sqrt(4*pi*alpha)", 2*M_W*math.sqrt(26.0/6)/math.sqrt(4*PI*ALPHA)),
        ("m_e * 120 * phi^17", M_E * 120 * PHI**17),
        ("m_e * 720 * phi^14", M_E * 720 * PHI**14),
        ("m_e * 12 * phi^22", M_E * 12 * PHI**22),
        ("m_W * 2/sqrt(4*pi*alpha*6/26)", M_W * 2 / math.sqrt(4*PI*ALPHA*6.0/26)),
        ("sqrt(m_e * m_P * phi^27)", math.sqrt(M_E * M_PLANCK * PHI**27)),
    ]

    print("  %-45s  %12s  %10s" % ("Expression", "Value (GeV)", "Error %"))
    print("  " + "-" * 72)
    for name, val in expressions:
        err = pct_err(val, VEV_EXP)
        marker = " ***" if err < 1 else (" ** " if err < 5 else "")
        print("  %-45s  %12.4f  %10.4f%%%s" % (name, val, err, marker))

    # The most interesting: vev = 2*m_W / g where g from derived sin^2_tW
    # This is SM identity but with derived sin^2_tW
    # To DERIVE vev we need to derive m_W independently
    print()
    print("  CRITICAL QUESTION: Can m_W be derived from m_e?")
    print("  m_W/m_e = %.2f, log_phi = %.6f" % (M_W/M_E, log_phi(M_W/M_E)))
    print("  Best (a,b): n=25 -> (a=5,b=0) or (a=-1,b=5) or (a=3,b=1+2/3)...")
    for n in [24, 25, 26]:
        pred = M_E * PHI**n
        err = pct_err(pred, M_W)
        print("    m_e * phi^%d = %.4f GeV (err %.2f%%)" % (n, pred, err))


# =============================================================================
# PART J: SUMMARY AND BEST CANDIDATES
# =============================================================================
def part_j():
    print()
    print("=" * 72)
    print("PART J: GRAND SUMMARY - BEST CANDIDATES FOR VEV DERIVATION")
    print("=" * 72)
    print()

    candidates = []

    # 1. vev = m_e * phi^27
    val1 = M_E * PHI**27
    candidates.append(("m_e * phi^27", val1, pct_err(val1, VEV_EXP),
                       "PATTERN", "From unified mass formula 5a+6b with (a=3,b=2)"))

    # 2. vev = m_P * alpha^(4*phi^2) * phi^27
    val2 = M_PLANCK * ALPHA**(4*PHI**2) * PHI**27
    candidates.append(("m_P * alpha^(4*phi^2) * phi^27", val2, pct_err(val2, VEV_EXP),
                       "PATTERN", "Equivalent to #1 via m_e = m_P*alpha^(4*phi^2)"))

    # 3. vev = 2*m_W/g with g from sin^2(tW)=6/26
    sin_tW = math.sqrt(6.0/26)
    e_em = math.sqrt(4*PI*ALPHA)
    g_w = e_em / sin_tW
    val3 = 2 * M_W / g_w
    candidates.append(("2*m_W*sin(tW)/e [sin^2=6/26]", val3, pct_err(val3, VEV_EXP),
                       "SM IDENTITY", "Requires m_W as input"))

    # 4. Using m_H formula: vev = m_H / sqrt(2*lambda) and m_H = m_W*(phi-8*alpha)
    # But lambda is not independently derived
    lam = M_H**2 / (2*VEV_EXP**2)
    val4 = M_H / math.sqrt(2*lam)  # tautological
    candidates.append(("m_H/sqrt(2*lambda)", val4, pct_err(val4, VEV_EXP),
                       "TAUTOLOGY", "lambda not independently derived"))

    # 5. Phi^(27 + 20*alpha) correction
    val5 = M_E * PHI**(27 + 20*ALPHA)
    candidates.append(("m_e * phi^(27+20*alpha)", val5, pct_err(val5, VEV_EXP),
                       "FITTING", "20 is a guess; could be edges/vertex"))

    # 6. m_e * 12 * phi^22
    val6 = M_E * 12 * PHI**22
    candidates.append(("m_e * 12 * phi^22", val6, pct_err(val6, VEV_EXP),
                       "FITTING", "12 = vertex degree, 22 arbitrary"))

    # 7. m_e * 120 * phi^17
    val7 = M_E * 120 * PHI**17
    candidates.append(("m_e * 120 * phi^17", val7, pct_err(val7, VEV_EXP),
                       "FITTING", "120 = vertices, 17 = tau exponent"))

    # 8. Try: vev = m_e * phi^(5^2 + alpha*phi*26)
    # 5^2 = 25, alpha*phi*26 ~ 0.0073*1.618*26 = 0.307
    val8 = M_E * PHI**(25 + ALPHA*PHI*26)
    candidates.append(("m_e * phi^(25+alpha*phi*26)", val8, pct_err(val8, VEV_EXP),
                       "SPECULATIVE", "Uses 5^2 and dim-of-phi-sector=26"))

    # 9. Search for vev = m_e * phi^N * alpha^M with special M values
    # We found that vev/m_e = phi^27.15..., and the fraction 0.15 ~ 20*alpha
    # 20 = number of tetrahedra per vertex in 600-cell!
    val9 = M_E * PHI**(27 + 20*ALPHA)
    err9 = pct_err(val9, VEV_EXP)
    # Let's compute what correction 'c' makes it exact
    exact_N = log_phi(VEV_EXP / M_E)
    c_exact = (exact_N - 27) / ALPHA
    print("  IMPORTANT: vev = m_e * phi^(27 + c*alpha)")
    print("  Exact c = (%.6f - 27) / alpha = %.6f" % (exact_N, c_exact))
    print()

    # What is c_exact close to?
    print("  c = %.6f. Close to:" % c_exact)
    c_tests = [
        ("20", 20),
        ("12+8", 20),
        ("phi^5", PHI**5),
        ("phi^4", PHI**4),
        ("3*phi^3", 3*PHI**3),
        ("2*phi^4", 2*PHI**4),
        ("4*phi^2", 4*PHI**2),
        ("phi^4+phi", PHI**4+PHI),
        ("26/phi", 26/PHI),
        ("20+phi", 20+PHI),
        ("20+1", 21),
        ("22", 22),
        ("7*pi", 7*PI),
        ("5^2/phi", 25/PHI),
        ("24/phi", 24/PHI),
    ]
    for name, val in c_tests:
        pred_N = 27 + val * ALPHA
        pred_vev = M_E * PHI**(pred_N)
        err = pct_err(pred_vev, VEV_EXP)
        print("    c = %-12s = %8.4f => vev = %.4f GeV (err %.4f%%)" %
              (name, val, pred_vev, err))

    print()
    print("  =" * 36)
    print("  RANKED CANDIDATES:")
    print("  =" * 36)
    print()

    candidates.sort(key=lambda x: x[2])
    print("  %-50s  %12s  %10s  %s" % ("Formula", "Value (GeV)", "Error %", "Status"))
    print("  " + "-" * 95)
    for name, val, err, status, note in candidates:
        print("  %-50s  %12.4f  %10.4f%%  %s" % (name, val, err, status))
        print("    Note: %s" % note)

    print()
    print("=" * 72)
    print("OVERALL ASSESSMENT")
    print("=" * 72)
    print()
    print("  STATUS: vev is NOT DERIVED from the 600-cell framework.")
    print()
    print("  BEST CANDIDATE: vev = m_e * phi^27 (err %.2f%%)" % pct_err(M_E * PHI**27, VEV_EXP))
    print("    - 27 = 5*3 + 6*2 in the unified mass formula")
    print("    - Quantum numbers (a=3, b=2) give norm = 9+6-4 = 11")
    print("    - This is PATTERN level, not DERIVED")
    print("    - The %.2f%% error is significant" % pct_err(M_E * PHI**27, VEV_EXP))
    print()
    print("  THE CORE PROBLEM: We cannot derive the overall SCALE.")
    print("    - The 600-cell gives us RATIOS (alpha, sin^2_tW, mass ratios)")
    print("    - But the absolute scale (vev, m_W, m_e in GeV) needs")
    print("      EITHER m_e as input OR a way to derive m_e/m_P beyond the")
    print("      pattern alpha^(4*phi^2).")
    print()
    print("  WHAT WOULD CONSTITUTE A DERIVATION:")
    print("    1. Derive m_W from geometric first principles (not phi^N fits)")
    print("    2. Then vev = 2*m_W/g with g from sin^2(tW)=6/26 (DERIVED)")
    print("    3. OR: Derive G_F geometrically, then vev = 1/sqrt(sqrt(2)*G_F)")
    print("    4. OR: Find a geometric principle fixing the ratio vev/m_P")
    print()
    print("  INTERESTING OBSERVATION:")
    exact_c = (log_phi(VEV_EXP / M_E) - 27) / ALPHA
    print("    vev = m_e * phi^(27 + c*alpha) where c = %.2f" % exact_c)
    print("    c is close to 20 (tetrahedra per vertex) -> err = %.2f%%" %
          pct_err(M_E * PHI**(27 + 20*ALPHA), VEV_EXP))
    print("    c is close to 4*phi^2 = %.2f -> err = %.2f%%" %
          (4*PHI**2, pct_err(M_E * PHI**(27 + 4*PHI**2*ALPHA), VEV_EXP)))
    print("    But these are FITS, not derivations.")
    print()
    print("  CATEGORY: SPECULATIV / PATTERN (not DERIVAT)")


# =============================================================================
# PART K: DEEP DIVE INTO BEST CANDIDATES
# =============================================================================
def part_k():
    print()
    print("=" * 72)
    print("PART K: DEEP DIVE INTO BEST CANDIDATES")
    print("=" * 72)
    print()

    # =========================================================================
    # CANDIDATE 1: vev = m_e * phi^(27 + 26*alpha)
    # =========================================================================
    print("  CANDIDATE 1: vev = m_e * phi^(27 + 26*alpha)")
    print("  " + "-" * 50)
    exact_N = log_phi(VEV_EXP / M_E)
    c_exact = (exact_N - 27) / ALPHA
    print("  Exact c = %.6f" % c_exact)
    print("  26 appears as: sin^2(tW) = 6/26, phi-sector dimension = 26")
    print("  26 = 20 (SU(2) directions) + 6 (U(1) directions) per vertex")
    print()

    val_26 = M_E * PHI**(27 + 26*ALPHA)
    err_26 = pct_err(val_26, VEV_EXP)
    print("  vev = m_e * phi^(27 + 26*alpha) = %.6f GeV" % val_26)
    print("  Error: %.4f%%" % err_26)
    print()

    # Test nearby integer values of c
    print("  Scanning integer c near 26:")
    for c in range(24, 30):
        val = M_E * PHI**(27 + c*ALPHA)
        err = pct_err(val, VEV_EXP)
        marker = " <--- BEST" if c == 26 else ""
        print("    c = %2d: vev = %.6f GeV, err = %.4f%%%s" % (c, val, err, marker))

    # c = 26 has a nice connection to the theory!
    # Write it differently: vev = m_e * phi^(27 + 26*alpha)
    # = m_e * phi^27 * phi^(26*alpha)
    # phi^(26*alpha) = phi^(26*7.297e-3) = phi^(0.18972)
    correction = PHI**(26*ALPHA)
    print()
    print("  Decomposition: vev = m_e * phi^27 * phi^(26*alpha)")
    print("    phi^(26*alpha) = phi^(%.6f) = %.6f" % (26*ALPHA, correction))
    print("    This is a %.4f%% correction to phi^27" % ((correction - 1)*100))
    print()

    # Connection: 26 = 1/sin^2(tW) at tree level (denominator)
    # So: vev = m_e * phi^(27 + alpha/sin^2_tW_600)
    # where sin^2_tW_600 = 6/26 => 1/sin^2 = 26/6
    # Wait, 26*alpha is not 1/sin^2. Let me be precise:
    # c = 26 exactly. And 26 appears in sin^2(tW) = 6/26.
    # Alternative reading: delta = alpha * 26 = alpha * (total structural dim per vertex)
    print("  PHYSICAL INTERPRETATION (SPECULATIVE):")
    print("    27 = base exponent (from 5*3 + 6*2 in mass formula)")
    print("    26*alpha = electromagnetic correction from all 26 structures/vertex")
    print("    This would be a 1-loop EM shift of the mass formula.")
    print("    26 = dim of phi-sector in 600-cell = denominator of sin^2(tW)")
    print()

    # Alternative: c = 1/(alpha*phi^2) ... no that's huge
    # c = ALPHA_INV / phi^2?
    # ALPHA_INV = 137.036, phi^2 = 2.618
    # 137.036 / 2.618 = 52.36 ... no
    # c = ALPHA_INV / 5 = 27.407? Close-ish
    val_ai5 = M_E * PHI**(27 + ALPHA_INV/5 * ALPHA)
    err_ai5 = pct_err(val_ai5, VEV_EXP)
    print("  Alternative: c = 1/(5*alpha)*alpha = 1/5 per vertex?")
    print("    vev = m_e * phi^(27+1/5) = %.4f GeV (err %.4f%%)" %
          (M_E * PHI**(27.2), pct_err(M_E * PHI**(27.2), VEV_EXP)))

    # =========================================================================
    # CANDIDATE 2: vev = m_e * 2*pi^2 * phi^21
    # =========================================================================
    print()
    print("  CANDIDATE 2: vev = m_e * 2*pi^2 * phi^21")
    print("  " + "-" * 50)
    val_2pi2 = M_E * 2 * PI**2 * PHI**21
    err_2pi2 = pct_err(val_2pi2, VEV_EXP)
    print("  2*pi^2 = Volume of unit S^3 = %.6f" % (2*PI**2))
    print("  21 = 3*7 or 5*3 + 6*1")
    print("  vev = m_e * V(S^3) * phi^21 = %.6f GeV" % val_2pi2)
    print("  Error: %.4f%%" % err_2pi2)
    print()
    print("  INTERPRETATION: The vev involves the volume of the 3-sphere")
    print("  on which the 600-cell is inscribed, times phi^21.")
    print("  21 = 5*3 + 6*1 => (a=3, b=1) in the mass formula")
    print("  norm(3+1*phi) = 9 + 3 - 1 = 11 (same as vev candidate (3,2)!)")
    # Actually: norm(a,b) = a^2 + ab - b^2
    norm_31 = 9 + 3 - 1
    norm_32 = 9 + 6 - 4
    print("  norm(a=3,b=1) = %d, norm(a=3,b=2) = %d" % (norm_31, norm_32))
    print("  Both give norm = 11! Coincidence?")
    print()
    # Check: is 2*pi^2 * phi^21 = phi^27 * f?
    ratio_cands = 2*PI**2 * PHI**21 / PHI**27
    print("  2*pi^2 * phi^21 / phi^27 = 2*pi^2 / phi^6 = %.6f" % ratio_cands)
    print("  2*pi^2 = %.6f, phi^6 = %.6f" % (2*PI**2, PHI**6))
    print("  Ratio = %.6f" % (2*PI**2 / PHI**6))
    print("  This is (2*pi^2/phi^6 - 1) = %.4f%% correction" %
          ((2*PI**2/PHI**6 - 1)*100))

    # =========================================================================
    # CANDIDATE 3: vev = m_e * phi^(27 + alpha * ALPHA_INV/phi^2)
    # =========================================================================
    # ALPHA_INV/phi^2 = 137.036/2.618 = 52.36
    # 52.36*alpha = 0.382 => phi^(27.382) ... too high
    # Not useful.

    # =========================================================================
    # CANDIDATE 4: G_F derivation
    # =========================================================================
    print()
    print("  CANDIDATE 3: G_F = alpha^k / m_e^2 approach")
    print("  " + "-" * 50)
    # k = 5.389430 from Part D
    # What if k = (phi^3 + phi) / phi = phi^2 + 1 = phi + 1 + 1 = phi^2 + 1?
    # phi^2 + 1 = 3.618. Not close.
    # k = 5 + alpha_s? 5 + 0.1179 = 5.1179. Not close.
    # k = 5 + 1/(2*phi) = 5.309. Closer.
    # k = 5 + 6*alpha/alpha_s? No, that's tiny.
    k_gf = math.log(G_F * M_E**2) / math.log(ALPHA)
    print("  k = ln(G_F*m_e^2)/ln(alpha) = %.6f" % k_gf)
    print()

    # Try rational multiples of phi
    best_k = []
    for p in range(1, 100):
        for q in range(1, 20):
            for sign in [1, -1]:
                for offset in range(-10, 15):
                    k_try = offset + sign * p * PHI / q
                    if abs(k_try - k_gf) / abs(k_gf) < 0.001:
                        pred_gf = ALPHA**k_try / M_E**2
                        pred_vev = 1.0 / math.sqrt(math.sqrt(2) * pred_gf)
                        err = pct_err(pred_vev, VEV_EXP)
                        best_k.append(("%d + %d*%d*phi/%d" %
                                       (offset, sign, p, q),
                                       k_try, pred_vev, err))

    best_k.sort(key=lambda x: x[3])
    print("  Best k = a + b*phi/c representations (err < 1%%):")
    seen = set()
    count = 0
    for name, k, vev_pred, err in best_k:
        if err < 1 and name not in seen:
            seen.add(name)
            print("    k = %-25s = %.6f, vev = %.4f GeV (err %.4f%%)" %
                  (name, k, vev_pred, err))
            count += 1
            if count >= 10:
                break

    # =========================================================================
    # CANDIDATE 5: vev = m_e * phi^(136/5)
    # =========================================================================
    print()
    print("  CANDIDATE 4: vev = m_e * phi^(136/5)")
    print("  " + "-" * 50)
    val_136_5 = M_E * PHI**(136.0/5)
    err_136_5 = pct_err(val_136_5, VEV_EXP)
    print("  136/5 = 27.2, phi^(136/5) = %.2f" % PHI**(136.0/5))
    print("  vev = %.6f GeV, err = %.4f%%" % (val_136_5, err_136_5))
    print("  136 = 8 * 17. 17 = tau lepton exponent!")
    print("  5 = diameter of 600-cell graph")
    print("  So: vev = m_e * phi^(8*17/5)")
    print("  = m_e * phi^(8 * n_tau / d_graph)")
    print("  LEVEL: PATTERN (numerology)")
    print()

    # =========================================================================
    # BONUS: Exact c decomposition
    # =========================================================================
    print()
    print("  EXACT DECOMPOSITION of c = %.6f:" % c_exact)
    print("  " + "-" * 50)
    # c = 26.384
    # c = 26 + 0.384
    # 0.384 ~ 1/phi^2 = 0.382 (err 0.5%)
    remainder = c_exact - 26
    print("  c = 26 + %.6f" % remainder)
    print("  1/phi^2 = %.6f (err %.2f%%)" %
          (1/PHI**2, pct_err(1/PHI**2, remainder)))
    print("  phi - 1 = 1/phi = %.6f (err %.2f%%)" %
          (1/PHI, pct_err(1/PHI, remainder)))
    print("  So c ~ 26 + 1/phi^2 = 26 + 2 - phi = 28 - phi")
    val_28mphi = M_E * PHI**(27 + (28 - PHI)*ALPHA)
    err_28mphi = pct_err(val_28mphi, VEV_EXP)
    print("  vev = m_e * phi^(27 + (28-phi)*alpha)")
    print("    = m_e * phi^(%.6f) = %.6f GeV (err %.4f%%)" %
          (27 + (28-PHI)*ALPHA, val_28mphi, err_28mphi))

    # Or: c = 26 + phi^(-2) where phi^(-2) = 2 - phi
    val_c_full = M_E * PHI**(27 + (26 + PHI**(-2))*ALPHA)
    err_c_full = pct_err(val_c_full, VEV_EXP)
    print()
    print("  REFINED: vev = m_e * phi^(27 + (26 + phi^{-2})*alpha)")
    print("    = m_e * phi^(%.8f)" % (27 + (26 + PHI**(-2))*ALPHA))
    print("    = %.6f GeV (err %.4f%%)" % (val_c_full, err_c_full))

    # Check: 26 + phi^{-2} = 26 + 2 - phi = 28 - phi
    print("  Note: 26 + phi^{-2} = 26 + (2-phi) = 28 - phi = %.6f" %
          (28 - PHI))

    # Very precise version: c = 26 + 0.384 = 26.384
    # What about c = 20 + 6/phi^2 = 20 + 6*(2-phi) = 20+12-6phi = 32-6phi
    val_test = 20 + 6*(2-PHI)
    print()
    print("  Test: c = 20 + 6*(2-phi) = 32-6*phi = %.6f" % val_test)
    val_t = M_E * PHI**(27 + val_test*ALPHA)
    print("  vev = %.6f GeV (err %.4f%%)" % (val_t, pct_err(val_t, VEV_EXP)))

    # What about c = 6/sin^2(tW) = 6/(6/26) = 26? That's just 26.
    # Or c = 1/sin^2(tW) = 26/6 = 4.333? No, that's small.
    # c = (1 + 1/sin^2(tW))*6 = (1+26/6)*6 = 32. Hmm.

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================
    print()
    print("  " + "=" * 60)
    print("  FINAL RANKING OF CANDIDATES")
    print("  " + "=" * 60)
    print()

    finals = [
        ("m_e * phi^(27+(26+phi^{-2})*alpha)", val_c_full, err_c_full,
         "FITTING", "c splits as 26+phi^{-2}, both geometric"),
        ("m_e * phi^(27+26*alpha)", val_26, err_26,
         "PATTERN", "26 = structural dim = denom(sin^2_tW)"),
        ("m_e * 2*pi^2 * phi^21", val_2pi2, err_2pi2,
         "PATTERN", "V(S^3)*phi^21, 21=5*3+6*1"),
        ("m_e * phi^(136/5)", val_136_5, err_136_5,
         "PATTERN", "136=8*17, 5=graph diameter"),
        ("m_e * phi^27", M_E * PHI**27, pct_err(M_E*PHI**27, VEV_EXP),
         "PATTERN", "Leading order, 27=5*3+6*2"),
    ]

    finals.sort(key=lambda x: x[2])
    print("  %-45s  %10s  %8s  %s" %
          ("Formula", "vev (GeV)", "Err %", "Level"))
    print("  " + "-" * 90)
    for name, val, err, level, note in finals:
        print("  %-45s  %10.4f  %8.4f  %s" % (name, val, err, level))
        print("    %s" % note)

    print()
    print("  HONEST ASSESSMENT:")
    print("    - The BEST fit is m_e*phi^(27+(26+phi^{-2})*alpha) at %.4f%%" % err_c_full)
    print("    - But this has 2 free parameters (26 and phi^{-2})")
    print("    - The SIMPLEST non-trivial is m_e*phi^(27+26*alpha) at %.4f%%" % err_26)
    print("    - 26 has clear geometric meaning (structures per vertex)")
    print("    - None of these constitute a DERIVATION")
    print("    - vev remains an OPEN PROBLEM in the 600-cell theory")


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print("EXP-193: DERIVATION OF THE ELECTROWEAK VEV")
    print("=" * 72)
    print()
    print("Target: vev = %.2f GeV" % VEV_EXP)
    print("PHI = %.10f" % PHI)
    print("ALPHA = %.10e (1/%.6f)" % (ALPHA, 1/ALPHA))
    print("sin^2(tW) = %.5f (theory: 6/26 = %.5f)" % (SIN2_TW, 6.0/26))
    print()

    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    part_g()
    part_h()
    part_i()
    part_j()
    part_k()

    print()
    print("EXP-193 COMPLETE.")
