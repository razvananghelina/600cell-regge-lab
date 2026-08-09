import math

a1 = 5
b1 = 6
phi = (1 + math.sqrt(5)) / 2
phi_prime = (1 - math.sqrt(5)) / 2
A_coef = 2 * math.pi
B_coef = -4 * a1 * phi**4
C_coef = 1
disc = B_coef**2 - 4 * A_coef * C_coef
alpha_fw = (-B_coef - math.sqrt(disc)) / (2 * A_coef)
C_val = 4.0 / (a1**2 + 1)
m_e = 0.51099895
m_mu_exp = 105.6583755
m_tau_exp = 1776.86
sigma_mu = 0.0000023
sigma_tau = 0.12
z_prime_mu = 1 + 1 * phi_prime
z_prime_tau = 1 + 2 * phi_prime
SEP = "=" * 80
print(SEP)
print("           LEPTON COEFFICIENT DISCRIMINATION TEST")
print(SEP)
print()

print("--- Framework constants ---")
print("  a1 =", a1)
print("  b1 =", b1)
print("  phi = {:.10f}".format(phi))
print("  phi_prime = {:.10f}".format(phi_prime))
print("  alpha_fw = {:.10f}  (= 1/{:.4f})".format(alpha_fw, 1.0/alpha_fw))
print("  C = 4/26 = {:.10f}".format(C_val))
print()

print("--- Galois conjugates ---")
print("  z_prime_mu  = {:.10f}".format(z_prime_mu))
print("  z_prime_tau = {:.10f}".format(z_prime_tau))
print("  |z_prime_mu|  = {:.10f}".format(abs(z_prime_mu)))
print("  |z_prime_tau| = {:.10f}".format(abs(z_prime_tau)))
sign_mu = 1 if z_prime_mu > 0 else -1
sign_tau = 1 if z_prime_tau > 0 else -1
print("  sign(z_prime_mu)  =", sign_mu)
print("  sign(z_prime_tau) =", sign_tau)
print()

print("--- Experimental masses ---")
print("  m_e   = {:.8f} MeV".format(m_e))
print("  m_mu  = {:.7f} +/- {:.7f} MeV".format(m_mu_exp, sigma_mu))
print("  m_tau = {:.2f} +/- {:.2f} MeV".format(m_tau_exp, sigma_tau))
print()

m_mu_bare = m_e * phi**11
m_tau_bare = m_e * phi**17
err_mu_bare = (m_mu_bare - m_mu_exp) / m_mu_exp * 100
err_tau_bare = (m_tau_bare - m_tau_exp) / m_tau_exp * 100
print("--- Bare predictions (delta=0) ---")
print("  m_mu_bare  = {:.4f} MeV  (error: {:+.4f}%)".format(m_mu_bare, err_mu_bare))
print("  m_tau_bare = {:.2f} MeV  (error: {:+.4f}%)".format(m_tau_bare, err_tau_bare))
print()
cands = [
    ("phi/(2*a1)=phi/10", phi / (2 * a1)),
    ("C+alpha=4/26+alpha", C_val + alpha_fw),
    ("1/b1=1/6", 1.0 / b1),
    ("(a1-1)/a1^2=4/25", (a1 - 1.0) / a1**2),
]

print("--- Candidate values ---")
for nm, val in cands:
    print("  {:<30s} = {:.10f}".format(nm, val))
print()

def compute_delta(c_ell, z_p):
    sgn = 1 if z_p >= 0 else -1
    return c_ell * sgn * abs(z_p)**0.75

def compute_predictions(c_ell):
    dmu = compute_delta(c_ell, z_prime_mu)
    dtau = compute_delta(c_ell, z_prime_tau)
    mmu = m_e * phi**(11 + dmu)
    mtau = m_e * phi**(17 + dtau)
    return mmu, mtau, dmu, dtau

print(SEP)
print("DETAILED RESULTS FOR EACH CANDIDATE")
print(SEP)

results = []
for nm, c_ell in cands:
    mmu, mtau, dmu, dtau = compute_predictions(c_ell)
    emu = (mmu - m_mu_exp) / m_mu_exp * 100
    etau = (mtau - m_tau_exp) / m_tau_exp * 100
    c2mu = ((mmu - m_mu_exp) / sigma_mu)**2
    c2tau = ((mtau - m_tau_exp) / sigma_tau)**2
    c2tot = c2mu + c2tau
    rms = math.sqrt((emu**2 + etau**2) / 2)
    results.append(dict(name=nm, c_ell=c_ell, delta_mu=dmu, delta_tau=dtau,
        m_mu=mmu, m_tau=mtau, err_mu=emu, err_tau=etau,
        chi2_mu=c2mu, chi2_tau=c2tau, chi2=c2tot, rms=rms))
    print()
    print("  Candidate:", nm)
    print("  c_ell = {:.10f}".format(c_ell))
    print("  delta_mu  = {:+.8f}".format(dmu))
    print("  delta_tau = {:+.8f}".format(dtau))
    print("  m_mu_pred  = {:.7f} MeV  (exp: {:.7f})  err: {:+.6f}%".format(mmu, m_mu_exp, emu))
    print("  m_tau_pred = {:.4f} MeV  (exp: {:.2f})  err: {:+.6f}%".format(mtau, m_tau_exp, etau))
    print("  chi2_mu  = {:.2e}".format(c2mu))
    print("  chi2_tau = {:.4f}".format(c2tau))
    print("  chi2_tot = {:.4f}".format(c2tot))
    print("  RMS error = {:.4f}%".format(rms))
print()
print(SEP)
print("SUMMARY COMPARISON TABLE  (sorted by RMS)")
print(SEP)
print()
hdr = "{:<24s} | {:>10s} | {:>10s} | {:>10s} | {:>10s} | {:>12s} | {:>4s}".format(
    "Candidate", "c_ell", "err_mu(%)", "err_tau(%)", "RMS(%)", "chi2_total", "Rank")
print(hdr)
print("-" * len(hdr))
rs = sorted(results, key=lambda r: r["rms"])
for rank, r in enumerate(rs, 1):
    r["rank"] = rank
for r in rs:
    short = r["name"]
    if len(short) > 24:
        short = short[:21] + "..."
    print("{:<24s} | {:10.8f} | {:+10.6f} | {:+10.6f} | {:10.6f} | {:12.4f} | {:4d}".format(
        short, r["c_ell"], r["err_mu"], r["err_tau"], r["rms"], r["chi2"], r["rank"]))
print()
print(SEP)
print("DETAILED CHI-SQUARED BREAKDOWN")
print(SEP)
print()
hdr2 = "{:<24s} | {:>12s} | {:>12s} | {:>12s}".format("Candidate", "chi2_mu", "chi2_tau", "chi2_total")
print(hdr2)
print("-" * len(hdr2))
for r in rs:
    short = r["name"]
    if len(short) > 24:
        short = short[:21] + "..."
    print("{:<24s} | {:12.2e} | {:12.4f} | {:12.4f}".format(
        short, r["chi2_mu"], r["chi2_tau"], r["chi2"]))

print()
print(SEP)
print("ANALYSIS AND VERDICT")
print(SEP)
print()
winner = rs[0]
runner = rs[1]
print("RANKING:")
for r in rs:
    print("  #{:d}  {:<30s}  RMS={:.4f}%  chi2={:.4f}".format(r["rank"], r["name"], r["rms"], r["chi2"]))
print()
print("BEST candidate:", winner["name"])
print("  c_ell = {:.10f}".format(winner["c_ell"]))
print("  RMS error = {:.4f}%".format(winner["rms"]))
print("  chi2 = {:.4f}".format(winner["chi2"]))
print()
print("RUNNER-UP:", runner["name"])
print("  c_ell = {:.10f}".format(runner["c_ell"]))
print("  RMS error = {:.4f}%".format(runner["rms"]))
print("  chi2 = {:.4f}".format(runner["chi2"]))
print()
if runner["chi2"] > 0 and winner["chi2"] > 0:
    ratio = runner["chi2"] / winner["chi2"]
    print("Chi2 ratio (runner/winner) = {:.2f}".format(ratio))
    if ratio > 10:
        print("  => STRONG discrimination: winner is clearly preferred")
    elif ratio > 3:
        print("  => MODERATE discrimination: winner is preferred")
    else:
        print("  => WEAK discrimination: candidates are comparable")
print()
print("--- Sub-percent accuracy check ---")
for r in rs:
    short = r["name"]
    mu_ok = abs(r["err_mu"]) < 1.0
    tau_ok = abs(r["err_tau"]) < 1.0
    both = "YES" if (mu_ok and tau_ok) else "NO"
    mstr = "Y" if mu_ok else "N"
    tstr = "Y" if tau_ok else "N"
    print("  {:<28s} |err_mu|<1%: {:<3s}  |err_tau|<1%: {:<3s}  Both: {}".format(short, mstr, tstr, both))
print()

print("--- Reverse-engineering exact c_ell ---")
n_mu_exact = math.log(m_mu_exp / m_e) / math.log(phi)
n_tau_exact = math.log(m_tau_exp / m_e) / math.log(phi)
print("  Exact exponent for muon:  {:.10f}  (bare: 11)".format(n_mu_exact))
print("  Exact exponent for tau:   {:.10f}  (bare: 17)".format(n_tau_exact))
delta_mu_exact = n_mu_exact - 11
delta_tau_exact = n_tau_exact - 17
print("  Required delta_mu = {:+.10f}".format(delta_mu_exact))
print("  Required delta_tau = {:+.10f}".format(delta_tau_exact))
factor_mu = (1 if z_prime_mu >= 0 else -1) * abs(z_prime_mu)**0.75
factor_tau = (1 if z_prime_tau >= 0 else -1) * abs(z_prime_tau)**0.75
c_from_mu = delta_mu_exact / factor_mu
c_from_tau = delta_tau_exact / factor_tau
print()
print("  c_ell needed for exact m_mu:  {:.10f}".format(c_from_mu))
print("  c_ell needed for exact m_tau: {:.10f}".format(c_from_tau))
if c_from_mu != 0:
    print("  Ratio c_tau/c_mu: {:.6f}".format(c_from_tau / c_from_mu))
print("  (If ratio ~ 1.000, a single c_ell fits both perfectly)")
print()
c_ideal = (c_from_mu + c_from_tau) / 2
print("  Average ideal c_ell: {:.10f}".format(c_ideal))
print()
print("  Closeness of each candidate to ideal average:")
for r in rs:
    short = r["name"]
    dist = abs(r["c_ell"] - c_ideal) / c_ideal * 100
    print("    {:<28s} c={:.10f}  dist={:.4f}%".format(short, r["c_ell"], dist))
print()
print(SEP)
print("CONCLUSION")
print(SEP)
print()
print("The coefficient", winner["name"], "provides the best fit")
print("with RMS={:.4f}% across both leptons.".format(winner["rms"]))
print()
avg = (c_from_mu + c_from_tau) / 2
print("Consistency: c_from_mu vs c_from_tau differ by {:.4f}%,".format(abs(c_from_mu - c_from_tau) / avg * 100))
print("indicating how well a SINGLE coefficient can describe both leptons.")
print()
print("--- Sigma-level for each candidate ---")
for r in rs:
    short = r["name"]
    stau = abs(r["m_tau"] - m_tau_exp) / sigma_tau
    smu = abs(r["m_mu"] - m_mu_exp) / sigma_mu
    print("  {:<24s}  tau: {:7.2f} sigma   mu: {:10.0f} sigma".format(short, stau, smu))
