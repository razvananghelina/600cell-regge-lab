# Regular-lapse identity: interval-assumption control failure

Date: 2026-08-16

Frozen protocol: `cf492b9`

Exact-rational correction: `70d4477`

Second artifact SHA-256:
`70e7c7545c20e44a1aacce643c9d40c229ac986feeeedf49905e66936fbc5e7b`

Status: **SYMBOLIC ASSUMPTION FAILURE; SCIENTIFIC OUTCOME NOT CLOSED**.

After replacing the accidental binary half by an exact rational, the second
targeted run returned `12/13`.  Exact Gram inertia, local volume, angle
products, anchored curvature sums, action, internal gradients and boundary
momenta all passed.  All twelve 100-digit controls again passed.

The sole failure was the no-log-cut classifier.  The symbol was declared
positive, but SymPy was not told the additional theorem domain `u<1/2`.
Calls to `re()` and `im()` therefore introduced `Abs` and `atan2` around
radicals such as `sqrt(1-2*u)` and could not simplify their real/imaginary
parts globally.

The already frozen branch proof uses exact sign polynomials on
`0<u<1/2`.  The allowed implementation correction is to compare each of the
eleven angle arguments with its explicit radical normal form and use those
same exact polynomial certificates to prove that its real scale factor is
positive.  No formula, domain, test point, tolerance or outcome rule changes.

