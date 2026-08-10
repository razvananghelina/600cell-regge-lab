#!/usr/bin/env python3
"""Exact finite-dimensional audit of the canonical-state shortcut.

The candidate shortcut is: choose a state, obtain its representation by GNS,
and read time from its modular flow.  This script checks the decisive finite
facts.  A tracial/Haar state has trivial modular flow.  A faithful nontracial
state on M2 has a modular frequency fixed by an arbitrary density-eigenvalue
ratio.  Its faithful GNS image is still M2 with a four-dimensional commutant;
the state cannot manufacture M4 or remove the multiplicity label.
"""

import sympy as sp


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
print("MISSING-LINK AUDIT: CAN A CANONICAL STATE SUPPLY REPRESENTATION AND TIME?")
print("=" * 78)

I2 = sp.eye(2)
units2 = []
for i in range(2):
    for j in range(2):
        unit = sp.zeros(2, 2)
        unit[i, j] = 1
        units2.append(unit)
E11, E12, E21, E22 = units2


def state(rho, observable):
    return sp.trace(rho * observable)


# Normalized trace: the finite analogue of the invariant Haar state on a
# matrix fibre.  Check traciality on a full algebra basis, not just examples.
rho_trace = I2 / 2
tracial = all(
    state(rho_trace, a * b) == state(rho_trace, b * a)
    for a in units2 for b in units2
)
check("normalized M2 state is tracial on the full matrix-unit basis", tracial)

logrho_trace = sp.log(sp.Rational(1, 2)) * I2
trace_generators = [logrho_trace * a - a * logrho_trace for a in units2]
check("tracial state has zero modular generator", all(x == sp.zeros(2) for x in trace_generators))

# Every state on a commutative algebra is tracial.  Diagonal matrices give a
# faithful finite model of the multiplication algebra of a finite connection
# sample.  Its density commutes with every observable, for arbitrary weights.
p = sp.symbols("p", positive=True)
rho_diag = sp.diag(p, 1 - p)
diagonal_basis = (E11, E22)
check("every faithful state on the commutative diagonal algebra is tracial",
      all(state(rho_diag, a*b) == state(rho_diag, b*a)
          for a in diagonal_basis for b in diagonal_basis))
check("commutative connection observables have trivial modular action",
      all(rho_diag*a == a*rho_diag for a in diagonal_basis))

# A nontracial faithful state creates a modular flow, but its frequency is the
# logarithm of the freely chosen density-eigenvalue ratio.
nontracial_defect = sp.simplify(
    state(rho_diag, E12*E21) - state(rho_diag, E21*E12))
check("M2 state is tracial iff p=1/2",
      nontracial_defect == 2*p - 1 and
      sp.solve(sp.Eq(nontracial_defect, 0), p) == [sp.Rational(1, 2)])

q = sp.symbols("q", positive=True)
logrho = sp.diag(sp.log(p), sp.log(q))
commutator_E12 = sp.simplify(logrho*E12 - E12*logrho)
check("E12 modular frequency is log(p/q)",
      sp.simplify(commutator_E12 - (sp.log(p)-sp.log(q))*E12)
      == sp.zeros(2))

omega_2 = sp.log(sp.Rational(2, 3) / sp.Rational(1, 3))
omega_3 = sp.log(sp.Rational(3, 4) / sp.Rational(1, 4))
check("two faithful states give inequivalent modular frequencies",
      omega_2 == sp.log(2) and omega_3 == sp.log(3) and omega_2 != omega_3,
      "p=2/3 gives log(2); p=3/4 gives log(3)")

omega = sp.symbols("omega", real=True)
p_from_omega = sp.exp(omega) / (1 + sp.exp(omega))
q_from_omega = 1 / (1 + sp.exp(omega))
check("any positive modular ratio can be inserted through the state",
      sp.simplify(p_from_omega / q_from_omega) == sp.exp(omega))

# Faithful GNS of M2 is its left-regular action on the four-dimensional
# Hilbert-Schmidt space.  It is A tensor I2, hence has image dimension four
# and a four-dimensional right-action commutant, not the full M4 algebra.
left_basis = [sp.kronecker_product(a, I2) for a in units2]
left_span = sp.Matrix.hstack(*(a.reshape(16, 1) for a in left_basis)).rank()
check("faithful M2 GNS image has complex dimension 4", left_span == 4)

ambient_units = []
for i in range(4):
    for j in range(4):
        unit = sp.zeros(4, 4)
        unit[i, j] = 1
        ambient_units.append(unit)
commutator_columns = []
for x in ambient_units:
    commutator_columns.append(sp.Matrix.vstack(*(
        (x*a-a*x).reshape(16, 1) for a in left_basis
    )))
commutator_map = sp.Matrix.hstack(*commutator_columns)
commutant_dim = 16 - commutator_map.rank()
check("faithful M2 GNS commutant has complex dimension 4", commutant_dim == 4)
check("a state on M2 cannot turn its image into M4",
      left_span == 4 and 4 < 16,
      "linear image dimension <= dim(M2)=4, whereas dim(M4)=16")

# The two copy-labelled vectors are exactly indistinguishable by A tensor I.
e0 = sp.Matrix((1, 0))
copy0 = sp.kronecker_product(e0, sp.Matrix((1, 0)))
copy1 = sp.kronecker_product(e0, sp.Matrix((0, 1)))
expectations_equal = all(
    (copy0.T*a*copy0)[0] == (copy1.T*a*copy1)[0]
    for a in left_basis
)
check("faithful GNS still cannot separate the multiplicity copy label",
      expectations_equal)

# Adding finite geometry symmetries in the most canonical way does not create
# thermal time either.  As a complete finite control, form C(Z3) crossed with
# the translation action of Z3.  Its nine basis elements delta_x U_g span M3,
# while the invariant coefficient-of-identity state remains the normalized
# matrix trace.  Verify traciality on all 81 ordered basis pairs.
def crossed_multiply(left, right):
    """Multiply basis labels (x,g) for delta_x U_g in C(Z3) rt Z3."""
    x, g = left
    y, h = right
    # (g.delta_y)(z)=delta_y(z-g), so support is y+g.
    if x != (y + g) % 3:
        return None
    return (x, (g + h) % 3)


def crossed_trace(label):
    if label is None:
        return sp.Integer(0)
    x, g = label
    return sp.Rational(1, 3) if g == 0 else sp.Integer(0)


crossed_basis = [(x, g) for x in range(3) for g in range(3)]
crossed_tracial = all(
    crossed_trace(crossed_multiply(a, b))
    == crossed_trace(crossed_multiply(b, a))
    for a in crossed_basis for b in crossed_basis
)
check("invariant state on the finite symmetry crossed product is tracial",
      crossed_tracial)

# Regular covariant representation: delta_x is |x><x| and U_g translates.
crossed_matrices = []
for x, g in crossed_basis:
    delta = sp.zeros(3, 3)
    delta[x, x] = 1
    shift = sp.zeros(3, 3)
    for source in range(3):
        shift[(source+g) % 3, source] = 1
    crossed_matrices.append(delta*shift)
crossed_span = sp.Matrix.hstack(
    *(matrix.reshape(9, 1) for matrix in crossed_matrices)
).rank()
check("finite translation crossed product is the full M3 algebra",
      crossed_span == 9)
check("its canonical invariant state still has zero modular flow",
      crossed_tracial and crossed_span == 9,
      "noncommutativity alone is insufficient; the state must be nontracial")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: Haar/tracial and commutative states have trivial modular flow.")
print("DERIVED: nontracial modular frequency is an arbitrary density ratio.")
print("DERIVED NEGATIVE: GNS of M2 remains M2 and does not select M4.")
print("DERIVED NEGATIVE: invariant finite crossed-product trace is modularly trivial.")
print("OPEN: a selected nontracial state on a selected noncommutative algebra.")
raise SystemExit(0 if passed == tests else 1)
