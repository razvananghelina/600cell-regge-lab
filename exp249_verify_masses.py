"""
EXP-249: Verify ALL mass predictions in the paper
===================================================
Simple, direct computation: m_f = m_e * phi^n
No tricks, no corrections. Just the formula as stated.
"""

PHI = (1 + 5**0.5) / 2
m_e = 0.51099895  # MeV, CODATA 2022

print("="*70)
print("EXP-249: MASS FORMULA VERIFICATION")
print("="*70)
print(f"phi = {PHI:.10f}")
print(f"m_e = {m_e} MeV")
print()

# Compute phi^n for reference
for n in range(0, 30):
    print(f"phi^{n:2d} = {PHI**n:15.4f}")

print()
print("="*70)
print("FERMION MASS PREDICTIONS: m_f = m_e * phi^n")
print("="*70)

fermions = [
    ("Electron", 0, 0, 0, 0.51099895, "input"),
    ("Muon",     1, 1, 1, 105.6584, "0.3%?"),
    ("Tau",      2, 1, 2, 1776.86, "0.3%?"),
    ("Up",       0, 3, -2, 2.16, "0.3%?"),
    ("Charm",    1, 2, 1, 1270.0, "1.0%?"),
    ("Top",      2, 4, 1, 172760.0, "1.8%?"),
    ("Down",     0, 1, 0, 4.67, "10%?"),
    ("Strange",  1, 1, 1, 93.4, "6.6%?"),
    ("Bottom",   2, -1, 4, 4180.0, "0.3%?"),
]

print(f"{'Fermion':10s} {'g':>2s} {'(a,b)':>8s} {'n':>4s} {'phi^n':>12s} "
      f"{'m_pred':>12s} {'m_exp':>12s} {'Error':>8s} {'Paper claims':>14s}")
print("-"*90)

for name, g, a, b, m_exp, paper_claim in fermions:
    n = 5*a + 6*b
    phi_n = PHI**n
    m_pred = m_e * phi_n
    if m_exp > 0 and name != "Electron":
        err = abs(m_pred - m_exp) / m_exp * 100
        print(f"{name:10s} {g:2d} ({a:2d},{b:2d}) {n:4d} {phi_n:12.4f} "
              f"{m_pred:12.4f} {m_exp:12.4f} {err:7.2f}% {paper_claim:>14s}")
    else:
        print(f"{name:10s} {g:2d} ({a:2d},{b:2d}) {n:4d} {phi_n:12.4f} "
              f"{m_pred:12.4f} {m_exp:12.4f} {'input':>8s} {paper_claim:>14s}")

print()
print("="*70)
print("WHAT PAPER TABLE CLAIMS vs WHAT FORMULA GIVES")
print("="*70)

paper_predictions = {
    "Muon": 106.0,
    "Tau": 1772.0,
    "Up": 2.15,
    "Charm": 1283.0,
    "Top": 169600.0,
    "Down": 5.15,
    "Strange": 87.3,
    "Bottom": 4191.0,
}

for name, g, a, b, m_exp, _ in fermions:
    if name == "Electron":
        continue
    n = 5*a + 6*b
    m_formula = m_e * PHI**n
    m_paper = paper_predictions[name]
    disc = abs(m_formula - m_paper) / m_paper * 100
    print(f"{name:10s}: n={n:2d}, formula gives {m_formula:10.2f} MeV, "
          f"paper says {m_paper:10.2f} MeV, discrepancy {disc:.1f}%")

# Check: does ANY consistent m_e give the paper's values?
print()
print("="*70)
print("WHAT m_e WOULD GIVE THE PAPER'S VALUES?")
print("="*70)
for name, g, a, b, m_exp, _ in fermions:
    if name == "Electron":
        continue
    n = 5*a + 6*b
    m_paper = paper_predictions[name]
    m_e_needed = m_paper / PHI**n
    print(f"{name:10s}: n={n:2d}, paper_pred={m_paper:.1f} => "
          f"needs m_e = {m_e_needed:.6f} MeV")

# Check phi^n via Fibonacci
print()
print("="*70)
print("PHI^N VIA FIBONACCI (cross-check)")
print("="*70)

def fib(n):
    if n <= 0: return 0
    if n == 1: return 1
    a, b = 0, 1
    for _ in range(n-1):
        a, b = b, a+b
    return b

for n in [3, 5, 11, 16, 17, 19, 26]:
    # phi^n = F(n-1) + F(n)*phi
    fn_minus_1 = fib(n-1) if n > 0 else 1
    fn = fib(n)
    phi_n_fib = fn_minus_1 + fn * PHI
    phi_n_direct = PHI**n
    print(f"phi^{n:2d} = F({n-1})+F({n})*phi = {fn_minus_1}+{fn}*phi = "
          f"{phi_n_fib:.6f} (direct: {phi_n_direct:.6f}, match: {abs(phi_n_fib-phi_n_direct)<1e-6})")

# Now check what experiments actually found
print()
print("="*70)
print("CHECKING: ARE THE PAPER ERRORS COMPARING TO RUNNING MASSES?")
print("="*70)

# Running masses at mu=2 GeV (MS-bar) vs pole masses
running_masses_2GeV = {
    "Up": 2.16,      # MS-bar at 2 GeV
    "Down": 4.67,     # MS-bar at 2 GeV
    "Strange": 93.4,  # MS-bar at 2 GeV
    "Charm": 1270.0,  # MS-bar at m_c
    "Bottom": 4180.0, # MS-bar at m_b
    "Top": 172760.0,  # pole mass (sometimes 162.5 GeV MS-bar)
    "Muon": 105.6584, # pole mass
    "Tau": 1776.86,   # pole mass
}

print("If comparing to PDG running masses:")
for name, g, a, b, m_exp, _ in fermions:
    if name == "Electron":
        continue
    n = 5*a + 6*b
    m_pred = m_e * PHI**n
    m_run = running_masses_2GeV[name]
    err = abs(m_pred - m_run) / m_run * 100
    print(f"  {name:10s}: pred={m_pred:.2f}, exp(running)={m_run:.2f}, error={err:.1f}%")
