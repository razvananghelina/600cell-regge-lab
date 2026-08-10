#!/usr/bin/env python3
"""Heat-kernel semiclassical-state gate for the HD-600 route.

This verifier tests whether the SU(2) heat-kernel semigroup removes the
vertical refinement freedom found in verify_hd600_connection_space.py.  No
particle target or fitted spectrum is used.
"""

import numpy as np
import sympy as sp


tests = passed = 0


def check(name, ok, detail=""):
    global tests, passed
    tests += 1
    ok = bool(ok)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("HD-600 GATE 3: HEAT-KERNEL STATES AND REFINEMENT WIDTHS")
print("=" * 78)

# In the SU(2) irrep of spin j=n/2, the bi-invariant Laplacian has Casimir
# c_n=n(n+2)/4.  The normalized heat kernel has scalar Fourier coefficient
# exp(-c_n tau).  Convolution multiplies Fourier coefficients.
s, t = sp.symbols("s t", positive=True)
casimirs = [sp.Rational(n*(n+2), 4) for n in range(21)]
semigroup_residuals = [
    sp.simplify(sp.exp(-casimir*s)*sp.exp(-casimir*t)
                - sp.exp(-casimir*(s+t)))
    for casimir in casimirs
]
check("SU(2) heat-kernel convolution semigroup holds in 21 irreps",
      all(value == 0 for value in semigroup_residuals),
      "Fourier coefficients exp[-j(j+1)t] multiply exactly")
check("heat-kernel normalization is fixed by the trivial irrep",
      casimirs[0] == 0 and sp.exp(-casimirs[0]*t) == 1,
      "integral K_t dg=1")

# Check the convention for a state centered at two background links.  With
# right noise U1=h1*x1 and U2=h2*x2, their product can be recentered at h2*h1;
# central heat-kernel noise is invariant under the conjugation that appears.
def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z,
                     w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X,
                     w*Z+x*Y-y*X+z*W))


def qconj(q):
    return np.array((q[0], -q[1], -q[2], -q[3]))


def random_unit(rng):
    q = rng.normal(size=4)
    return q/np.linalg.norm(q)


rng = np.random.default_rng(235600)
centering_residual = 0.0
for _ in range(100):
    h1, h2, noise1, noise2 = (random_unit(rng) for _ in range(4))
    direct = qmul(qmul(h2, noise2), qmul(h1, noise1))
    conjugated_noise2 = qmul(qconj(h1), qmul(noise2, h1))
    recentered = qmul(qmul(h2, h1),
                      qmul(conjugated_noise2, noise1))
    centering_residual = max(centering_residual,
                             float(np.linalg.norm(direct-recentered)))
check("two fine states recenter on the composed background holonomy",
      centering_residual < 2e-15,
      f"100 deterministic quaternion trials, max residual={centering_residual:.3e}")

# For independent equal fine noises, convolution gives tau_coarse=2*tau_fine.
tau_coarse, tau_fine = sp.symbols("tau_coarse tau_fine", positive=True)
equal_split_solution = sp.solve(
    sp.Eq(2*tau_fine, tau_coarse), tau_fine)[0]
check("independent equal subedges force tau_fine=tau_coarse/2",
      equal_split_solution == tau_coarse/2)
level, tau_zero = sp.symbols("level tau_zero", integer=True, nonnegative=True)
check("repeated binary refinement leaves one arbitrary global width",
      sp.simplify((tau_zero/2**(level+1))*2-tau_zero/2**level) == 0,
      "tau_level=tau_0/2^level; tau_0 remains free")

# Attack the hidden independence assumption in the tangent/Gaussian limit.
# The most general exchange-symmetric covariance of two fine noises is
# C=[[a,b],[b,a]].  Fixing the product/sum marginal only fixes a+b.  The
# difference-mode variance nu=a-b remains arbitrary and positive.
sigma2, nu = sp.symbols("sigma2 nu", positive=True)
a = sigma2/4+nu/2
b = sigma2/4-nu/2
covariance = sp.Matrix(((a, b), (b, a)))
sum_mode = sp.Matrix((1, 1))
difference_mode = sp.Matrix((1, -1))
sum_variance = sp.simplify((sum_mode.T*covariance*sum_mode)[0])
difference_variance = sp.simplify(
    (difference_mode.T*covariance*difference_mode)[0])
check("coarse heat-kernel marginal fixes only the sum-mode variance",
      sum_variance == sigma2,
      f"Var(x1+x2)={sum_variance}")
check("the vertical correlation width remains a free positive parameter",
      difference_variance == 2*nu,
      f"Var(x1-x2)={difference_variance}")

# Independence is b=0, hence nu=sigma2/2.  It is one point in the allowed
# family, not a consequence of the coarse marginal or exchange symmetry.
independent_nu = sp.solve(sp.Eq(b, 0), nu)[0]
check("factorized fine heat kernels select one correlation only by assumption",
      independent_nu == sigma2/2,
      "b=0 picks nu=sigma2/2")
covariance_examples = [covariance.subs({sigma2: 2, nu: value})
                       for value in (sp.Rational(1, 2), 1, 2)]
check("three inequivalent positive correlated states have the same coarse width",
      all(matrix.det() > 0 for matrix in covariance_examples) and
      len({tuple(matrix) for matrix in covariance_examples}) == 3 and
      all((sum_mode.T*matrix*sum_mode)[0] == 2
          for matrix in covariance_examples),
      "nu=(1/2,1,2), all Var(sum)=2")

print("\n" + "=" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("=" * 78)
print("SU2_HEAT_KERNEL_SEMIGROUP=DERIVED")
print("INDEPENDENT_EQUAL_EDGE_WIDTH_RECURSION=DERIVED")
print("ABSOLUTE_HEAT_WIDTH=OPEN_ONE_GLOBAL_PARAMETER")
print("FINE_EDGE_INDEPENDENCE=STRUCTURAL_NOT_DERIVED")
print("CORRELATED_SEMICLASSICAL_STATES=DERIVED_VERTICAL_FREEDOM")
print("HEAT_KERNEL_STATE_SELECTS_DIRAC_SCALE=DERIVED_NEGATIVE")
print("SM_TARGET_COMPARISON=NOT_PERFORMED")

if passed != tests:
    raise SystemExit(1)
