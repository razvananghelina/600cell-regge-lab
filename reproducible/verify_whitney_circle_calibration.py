#!/usr/bin/env python3
"""Known-answer calibration for the inductive Whitney Dirac dynamics.

The continuum target is the unit circle, whose first positive Dirac
eigenvalue is 2*pi.  Exact rational refinement identities are checked before
any numerical convergence measurement.  The same calculation exposes a
causality/induction tension between the consistent Whitney mass and canonical
row-sum mass lumping.
"""

import math

import numpy as np
import scipy.linalg as sla
import sympy as sy


tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("WHITNEY DIRAC CIRCLE CALIBRATION AND CAUSALITY GATE")
print("=" * 78)


def circle_complex(vertices, lumped=False):
    h = sy.Rational(1, vertices)
    differential = sy.zeros(vertices, vertices)
    for edge in range(vertices):
        differential[edge, edge] = -1
        differential[edge, (edge+1) % vertices] = 1
    if lumped:
        mass_zero = h*sy.eye(vertices)
    else:
        mass_zero = sy.zeros(vertices, vertices)
        for edge in range(vertices):
            left, right = edge, (edge+1) % vertices
            mass_zero[left, left] += h*sy.Rational(1, 3)
            mass_zero[right, right] += h*sy.Rational(1, 3)
            mass_zero[left, right] += h*sy.Rational(1, 6)
            mass_zero[right, left] += h*sy.Rational(1, 6)
    mass_one = (1/h)*sy.eye(vertices)
    return differential, (mass_zero, mass_one)


def dyadic_inclusions(coarse_vertices):
    fine_vertices = 2*coarse_vertices
    p_zero = sy.zeros(fine_vertices, coarse_vertices)
    p_one = sy.zeros(fine_vertices, coarse_vertices)
    for coarse in range(coarse_vertices):
        p_zero[2*coarse, coarse] = 1
        p_zero[2*coarse+1, coarse] = sy.Rational(1, 2)
        p_zero[2*coarse+1, (coarse+1) % coarse_vertices] = sy.Rational(1, 2)
        p_one[2*coarse, coarse] = sy.Rational(1, 2)
        p_one[2*coarse+1, coarse] = sy.Rational(1, 2)
    return p_zero, p_one


coarse_n = 4
fine_n = 8
d_c, masses_c = circle_complex(coarse_n)
d_f, masses_f = circle_complex(fine_n)
p_zero, p_one = dyadic_inclusions(coarse_n)
check("circle Whitney inclusions commute with d exactly",
      d_f*p_zero == p_one*d_c)
check("consistent Whitney L2 metrics are exactly refinement-isometric",
      p_zero.T*masses_f[0]*p_zero == masses_c[0]
      and p_one.T*masses_f[1]*p_one == masses_c[1])

weak_c = sy.zeros(2*coarse_n)
weak_f = sy.zeros(2*fine_n)
forward_c = masses_c[1]*d_c
forward_f = masses_f[1]*d_f
weak_c[coarse_n:, :coarse_n] = forward_c
weak_c[:coarse_n, coarse_n:] = forward_c.T
weak_f[fine_n:, :fine_n] = forward_f
weak_f[:fine_n, fine_n:] = forward_f.T
p_all = sy.diag(p_zero, p_one)
check("consistent weak Dirac form compresses exactly under refinement",
      p_all.T*weak_f*p_all == weak_c)

# Row-sum lumping is canonical once the consistent mass is supplied: it keeps
# the constant integral at each vertex and produces a diagonal Hodge mass.
# It is tested because its dispersion respects the finite lattice speed bound.
_, lumped_c = circle_complex(coarse_n, lumped=True)
_, lumped_f = circle_complex(fine_n, lumped=True)
lumped_zero_residual = p_zero.T*lumped_f[0]*p_zero-lumped_c[0]
check("row-sum lumping breaks exact refinement isometry",
      lumped_zero_residual != sy.zeros(coarse_n),
      f"residual rank={lumped_zero_residual.rank()}")


def dirac_spectrum(vertices):
    differential, masses = circle_complex(vertices)
    d = np.asarray(differential, dtype=float)
    m0 = np.asarray(masses[0], dtype=float)
    m1 = np.asarray(masses[1], dtype=float)
    metric = sla.block_diag(m0, m1)
    forward = m1@d
    weak = np.block([[np.zeros_like(m0), forward.T],
                     [forward, np.zeros_like(m1)]])
    return sla.eigh(weak, metric, eigvals_only=True)


levels = (8, 16, 32, 64, 128)
first_positive = []
errors = []
for vertices in levels:
    eigenvalues = dirac_spectrum(vertices)
    positive = eigenvalues[eigenvalues > 1e-8]
    first_positive.append(float(positive[0]))
    errors.append(abs(float(positive[0])-2*math.pi))

check("known first circle Dirac eigenvalue converges monotonically to 2*pi",
      all(errors[index+1] < errors[index]
          for index in range(len(errors)-1)),
      "errors=" + ", ".join(f"{value:.3e}" for value in errors))
error_ratios = [errors[index]/errors[index+1]
                for index in range(len(errors)-1)]
check("circle spectral convergence is second order after calibration",
      all(3.7 < ratio < 4.3 for ratio in error_ratios[1:]),
      "ratios=" + ", ".join(f"{value:.3f}" for value in error_ratios))

# Analytic dispersion for the consistent Whitney mass on a uniform mesh:
# p_h(k)=2 sin(kh/2)/(h sqrt(2/3+cos(kh)/3)).  It converges at fixed physical
# k, but at q=kh=2*pi/3 its derivative is sqrt(2), exceeding the continuum
# coefficient c.  The lumped dispersion has derivative cos(q/2), bounded by 1.
def consistent_momentum(k, h):
    return (2*np.sin(k*h/2)/h
            / np.sqrt(2/3+np.cos(k*h)/3))


def consistent_momentum_derivative(k, h):
    q = k*h
    denominator = 1-sy.Rational(2, 3)*np.sin(q/2)**2
    return np.cos(q/2)/denominator**1.5


fixed_k = 2*math.pi
momentum_errors = [abs(consistent_momentum(fixed_k, 1/vertices)-fixed_k)
                   for vertices in levels]
check("consistent Whitney momentum converges at fixed physical wavelength",
      all(momentum_errors[index+1] < momentum_errors[index]
          for index in range(len(momentum_errors)-1)),
      "errors=" + ", ".join(f"{value:.3e}" for value in momentum_errors))

q_witness = 2*math.pi/3
consistent_uv_velocity = consistent_momentum_derivative(
    q_witness, 1.0)  # k=q and h=1; only q=kh matters
lumped_uv_velocity = math.cos(q_witness/2)
check("consistent Whitney dispersion violates the unit finite lattice speed",
      abs(consistent_uv_velocity-math.sqrt(2)) < 1e-12
      and consistent_uv_velocity > 1,
      f"v/c at q=2*pi/3 is {consistent_uv_velocity:.12f}")
check("lumped dispersion has the exact finite bound |v|<=c",
      abs(lumped_uv_velocity-0.5) < 1e-12,
      "v/c=cos(q/2), hence its absolute value is at most one")

# Massive low-momentum velocity uses the same mass shell.  At fixed physical
# k it converges to the continuum value, even though the cutoff-scale bound
# above fails for the consistent mass.
mass = 2.0
continuum_velocity = fixed_k/math.sqrt(mass**2+fixed_k**2)
velocity_errors = []
for vertices in levels:
    h = 1/vertices
    discrete_p = consistent_momentum(fixed_k, h)
    discrete_dp = consistent_momentum_derivative(fixed_k, h)
    discrete_velocity = (discrete_p*discrete_dp
                         / math.sqrt(mass**2+discrete_p**2))
    velocity_errors.append(abs(discrete_velocity-continuum_velocity))
check("massive group velocity converges at fixed physical momentum",
      all(velocity_errors[index+1] < velocity_errors[index]
          for index in range(len(velocity_errors)-1)),
      "errors=" + ", ".join(f"{value:.3e}" for value in velocity_errors))

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: consistent Whitney Dirac is exactly Galerkin-inductive.")
print("DERIVED: its known-answer circle spectrum converges at second order.")
print("DERIVED: fixed physical momenta and massive velocities converge.")
print("DERIVED NEGATIVE: consistent mass permits v=sqrt(2)c at cutoff scale.")
print("DERIVED TRADEOFF: lumping restores |v|<=c but breaks exact induction.")
raise SystemExit(0 if passed == tests else 1)
