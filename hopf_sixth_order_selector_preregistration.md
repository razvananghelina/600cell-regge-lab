# Preregistration: first canonical higher-moment Hopf selector

Date: 2026-08-10

## Provenance caveat

This is not a blind prediction.  Before this file was written, an exploratory
symbolic calculation had already shown the candidate identities

```text
S2 = 2 r^2,
S4 = (6/5) r^4,
critical S6 values = {26/25, 4/5, 34/45} on r=1.
```

It also suggested that the three values occur on the icosahedral fivefold,
twofold and threefold axes.  The purpose of this preregistration is to freeze
an exhaustive proof protocol capable of refuting that suggestion.  It must
not later be described as a blind hit.

## Fixed construction

Use the standard twelve exact icosahedron vertices

```text
(0, +/-1, +/-phi),
(+/-1, +/-phi, 0),
(+/-phi, 0, +/-1),
phi=(1+sqrt(5))/2.
```

Antipodal pairs give the six unoriented `C10`/Hopf axes derived independently
in `verify_hopf_symmetry_selector.py`.  For their exact rank-one projectors
`P_i`, define the canonical even moment polynomials

`S_(2m)(x) = sum_i (x^T P_i x)^m`.

No coefficient other than the equal weight forced by the transitive
six-axis action may be inserted.

## Exact gates

1. Reconstruct the icosahedron combinatorially: twelve vertices, thirty
   maximal-dot-product edges and twenty triangular faces.
2. Deduplicate the vertex, edge-midpoint and face-centre lines.  They must give
   respectively `6`, `15` and `10` unoriented axes.  This independently
   realizes the expected `C10`, `C4` and `C6` symmetry-axis orbits.
3. Derive `S2`, `S4` and `S6` in `Q(sqrt(5))[x,y,z]`.  Check exact radiality of
   `S2` and `S4`, and exact non-radiality of `S6`.
4. On `x^2+y^2+z^2=1`, form the Lagrange ideal

   `grad(S6)-2 lambda (x,y,z)=0`.

   Compute an exact Gröbner basis over `Q(sqrt(5))`.  The ideal must be
   zero-dimensional.  Extract and record its monic elimination polynomial in
   `lambda`; do not infer exhaustivity from sampled symmetry axes.
5. For each real root `lambda`, specialize the ideal, compute a lexicographic
   Gröbner basis and count its standard monomials.  This is the quotient-ring
   dimension, hence the number of complex solutions counted with
   multiplicity.
6. Independently verify that every signed point on the corresponding
   combinatorial axis orbit is a real critical point with that value.  Only if
   the number of distinct exhibited points equals the quotient dimension may
   the orbit be called exhaustive and reduced.
7. Use Euler homogeneity, `lambda=3 S6` on the unit sphere, and compactness to
   classify the global minimum and maximum.  An intermediate critical orbit
   must not be called a saddle unless its constrained Hessian is also checked.

## Decision boundary

- If `S2` or `S4` is non-radial, the earlier statement that quadratic and
  quartic moments cannot select a fibration is refuted.
- If the Lagrange ideal is not zero-dimensional or the elimination polynomial
  has additional real roots, the exploratory three-orbit picture is
  incomplete.
- If the exhibited real orbit count is below the quotient dimension, generic
  or complex critical points remain and exhaustivity is unproved.
- If the six Hopf axes are exactly the global maxima of `S6`, then the
  potential `-g S6`, with `g>0`, has six unoriented degenerate Hopf minima.
  This is a **DERIVED mathematical symmetry-breaking mechanism**, conditional
  on that sign.
- Neither the sign `g>0`, its magnitude, a kinetic anisotropy `r`, nor a
  physical identification of the order parameter is supplied by the finite
  geometry.  Until derived from an existing action, all of those remain
  **OPEN**.  The verifier may not promote the conditional potential to a
  physical prediction.

The exact appearance of `1/5` in the projector Gram matrix and `6/5` in the
quartic moment may be compared with the symbol `a1=5` only as **STRUCTURAL**;
no coupling or speed formula follows from these moment identities.
