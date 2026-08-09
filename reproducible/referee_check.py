"""
Adversarial referee-style checks for mathematically sensitive claims in
one_integer_paper.tex.

This script is intentionally different from the main reproducibility suite:
it verifies some exact arithmetic claims, but it also flags proof gaps and
overclaims as WARNING / FAIL when the current paper text is stronger than the
mathematics implemented here can justify.

It does not attempt to settle the physics. It only checks:
  - exact algebraic identities that are straightforward to verify;
  - whether some cited mathematical justifications are actually correct;
  - whether some claimed uniqueness proofs are globally established or only
    bounded-search evidence;
  - where external input or continuum assumptions remain.

Expected use:
    python referee_check.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass


A1 = 5
B1 = A1 + 1
PHI = (1 + math.sqrt(5)) / 2
PHIP = (1 - math.sqrt(5)) / 2


@dataclass
class Result:
    level: str
    name: str
    detail: str


results: list[Result] = []


def add(level: str, name: str, detail: str) -> None:
    results.append(Result(level=level, name=name, detail=detail))
    print(f"[{level:<7}] {name}")
    print(f"          {detail}")


def pass_(name: str, detail: str) -> None:
    add("PASS", name, detail)


def warn(name: str, detail: str) -> None:
    add("WARNING", name, detail)


def fail(name: str, detail: str) -> None:
    add("FAIL", name, detail)


def fib(n: int) -> int:
    if n == 0:
        return 0
    if n > 0:
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a + b
        return a
    # F(-n) = (-1)^(n+1) F(n)
    return (-1) ** ((-n) + 1) * fib(-n)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = math.isqrt(n)
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


def norm_from_nt(n: int, t: int) -> int:
    # From the paper's parameterization:
    # (a,b)=(-n-6t, n+5t),  N(z)=-(n^2+9nt+19t^2)
    return -(n * n + 9 * n * t + 19 * t * t)


def admissible_c1(norm_abs: int) -> bool:
    # Theorem C1 in the paper allows 0, 1, or split/ramified primes p ≡ 0, ±1 mod 5
    return norm_abs in (0, 1) or (is_prime(norm_abs) and norm_abs % 5 in (0, 1, 4))


def section(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


section("1. Exact Arithmetic Claims")

# Unique integer
solutions = [a for a in range(1, 20) if math.factorial(a) == 4 * a * (a + 1)]
if solutions == [5]:
    pass_(
        "Diophantine uniqueness",
        f"Exact scan on a=1..19 gives the unique solution {solutions[0]} to a! = 4a(a+1).",
    )
else:
    fail("Diophantine uniqueness", f"Unexpected solutions found: {solutions}")

# Bootstrap uniqueness
bootstrap_hits = []
for a in range(1, 26):
    lhs = math.cos(2 * math.pi / a)
    rhs = (math.sqrt(a) - 1) / 4
    if abs(lhs - rhs) < 1e-12:
        bootstrap_hits.append(a)

if bootstrap_hits == [5]:
    pass_(
        "Bootstrap identity on finite domain",
        "For integers 1..25, cos(2pi/a) = (sqrt(a)-1)/4 holds only at a=5; "
        "for a>=26 the RHS exceeds 1, so no further solutions exist.",
    )
else:
    fail("Bootstrap identity", f"Unexpected hits: {bootstrap_hits}")

# Generation arithmetic core
generation_ks = [k for k in range(-8, 9) if fib(k - 1) == 1]
generation_bs = [fib(k) for k in generation_ks]
if generation_ks == [0, 2, 3] and generation_bs == [0, 1, 2]:
    pass_(
        "Generation arithmetic lemma",
        "On the line a=1, the unit solutions phi^k have k={0,2,3}, hence b={0,1,2}. "
        "This arithmetic step is correct.",
    )
else:
    fail(
        "Generation arithmetic lemma",
        f"Found k={generation_ks}, b={generation_bs}, expected k=[0,2,3], b=[0,1,2].",
    )


section("2. Coupling and Laplacian Checks")

# Alpha quadratic
A = 2 * math.pi
B = -4 * A1 * PHI**4
C = 1
disc = B * B - 4 * A * C
alpha = (-B - math.sqrt(disc)) / (2 * A)
alpha_inv = 1 / alpha
if abs(alpha_inv - 137.035999084) / 137.035999084 < 2e-6:
    pass_(
        "Alpha quadratic reproduces the claimed number",
        f"Solving 2*pi*alpha^2 - 4*a1*phi^4*alpha + 1 = 0 gives 1/alpha = {alpha_inv:.9f}.",
    )
else:
    fail("Alpha quadratic", f"Unexpected value 1/alpha = {alpha_inv:.9f}.")

# Icosahedron eigenvalue products
distinct_product = (5 - math.sqrt(5)) * 6 * (5 + math.sqrt(5))
full_product = (5 - math.sqrt(5)) ** 3 * 6**5 * (5 + math.sqrt(5)) ** 3
tree_count = round(full_product / 12)

if abs(distinct_product - 120) < 1e-10:
    pass_(
        "Distinct nonzero Laplacian eigenvalue product",
        "For the icosahedron, the distinct nonzero Laplacian eigenvalues satisfy "
        "(5-sqrt(5))*6*(5+sqrt(5)) = 120 exactly.",
    )
else:
    fail("Distinct eigenvalue product", f"Computed value {distinct_product}.")

if tree_count == 5_184_000:
    pass_(
        "Kirchhoff separated correctly from the distinct-eigenvalue identity",
        "For the icosahedron, tau = (1/12) * product(lambda_i with multiplicity) = 5,184,000. "
        "Therefore Kirchhoff must be kept separate from the algebraic identity "
        "(5-sqrt(5))*6*(5+sqrt(5)) = 120 for the distinct nonzero eigenvalues.",
    )
else:
    warn(
        "Kirchhoff sanity check",
        f"Computed tau={tree_count}; regardless of normalization details, Kirchhoff concerns the full spectrum "
        "with multiplicities, not the distinct-eigenvalue product.",
    )


section("3. Uniqueness-Proof Stress Test")

# Search for C1-admissible solutions beyond |t| <= 15
n_values = sorted({0, 3, 5, 11, 16, 17, 19, 26})
counterexamples: list[str] = []
for n in n_values:
    for t in range(16, 201):
        norm_abs = abs(norm_from_nt(n, t))
        if admissible_c1(norm_abs):
            counterexamples.append(f"n={n}, t={t}, |N|={norm_abs}")
            break

if counterexamples:
    warn(
        "Bounded search does not prove global uniqueness",
        "The current proof searches only |t|<=15, but C1-admissible values exist beyond that range, e.g. "
        + "; ".join(counterexamples[:4])
        + ". Therefore the present theorem still needs an analytic global bound or a complete infinite-family exclusion.",
    )
else:
    pass_(
        "No C1-admissible values found beyond |t|<=15 in a wider scan",
        "This does not prove the theorem, but the stress test did not find easy counterexamples.",
    )


section("4. Claim-Strength and Assumption Checks")

warn(
    "Generation theorem still uses physical assumptions",
    "The arithmetic core is correct, but the step 'stable state iff E(z)=0' and the choice of the "
    "tight-binding / nearest-neighbor dynamics are additional modeling assumptions, not purely forced consequences "
    "of algebraic number theory alone.",
)

warn(
    "Spectral action to full SM bosonic Lagrangian",
    "Finite checks can verify simplicial counts, Hodge decomposition, and Seeley-DeWitt-like coefficients, "
    "but they do not by themselves establish the continuum Chamseddine-Connes theorem for this discrete S^3 setup. "
    "This should be stated as structural reproduction, not strict identity with the full SM Lagrangian.",
)

warn(
    "PPN gamma = 1 claim is stronger than the discrete identity proved",
    "The exact graph identity h = d_0 Phi is meaningful, but it is weaker than a full continuum PPN derivation. "
    "A safer phrasing is that the discrete model has no independent static scalar mode in the tested sector.",
)

warn(
    "Zero-free-parameters language needs precision",
    "Dimensionless internal formulas may be parameter-free, but the strongest m_Z comparison uses the measured "
    "running ratio alpha(m_Z)/alpha(0), and some best-fit Higgs numbers use experimental m_W. "
    "These are external inputs, even if explicitly acknowledged in the text.",
)


section("5. Summary")

counts = {"PASS": 0, "WARNING": 0, "FAIL": 0}
for item in results:
    counts[item.level] += 1

print(f"PASS    : {counts['PASS']}")
print(f"WARNING : {counts['WARNING']}")
print(f"FAIL    : {counts['FAIL']}")

print()
print("Interpretation:")
print("- PASS: exact arithmetic / algebraic check succeeded.")
print("- WARNING: the mathematics may be internally consistent, but the proof in the paper is incomplete or the claim is stronger than what was shown.")
print("- FAIL: a cited mathematical justification is incorrect as stated.")

if counts["FAIL"] == 0:
    print()
    print("No outright mathematical misstatement was detected by this audit.")
else:
    print()
    print("At least one claim should be corrected in the paper before calling the argument rigorous.")
