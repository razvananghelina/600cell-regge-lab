# The unique copy-symmetric central flux repairs the circle obstruction

Date: 2026-08-11

Preregistration commit: `2b83f31`

Targeted verifier:
`reproducible/verify_whitney_central_flux_circle.py`

Targeted result: **13/13 PASS**.  The verifier is registered exactly once.
The full suite was not run by explicit user request.

## Headline

The pure-penalty failure on the circle is repaired by one face-local flux with
no fitted coefficient.

> **DERIVED POSITIVE ON S1:** equal averaging of the two copies at every
> vertex, followed by the assembled incidence block, preserves both Betti
> zero modes at every tested finite stiffness.

> **DERIVED NUMERICAL:** for four fixed values of `kappa`, the physical low
> branch converges to the known unit-circle Kähler--Dirac eigenvalue `2*pi`
> with unit limiting velocity.

The mismatch branch remains separated at exactly `12 kappa/h`.

This rescues finite stiffness as a numerical continuum mechanism in one
dimension.  It does not yet select three-dimensional dynamics or a physical
value of `kappa`.

## Why the coefficient is unique

At global vertex `n`, the two incident local copies are `a_n` and
`b_{n-1}`.  Let a vertex-local linear recovery map use weights `(w_L,w_R)`.

Copy-exchange symmetry requires

\[
 w_L=w_R.
\]

Reproducing a conforming value requires

\[
 w_L+w_R=1.
\]

Therefore, uniquely,

\[
 w_L=w_R=\frac12,
\]

and

\[
 (L_0x)_n=\frac{a_n+b_{n-1}}2.
\]

The verifier solves these two conditions symbolically; no numerical search is
involved.

## Flux construction

Let `J_0` copy global vertex values, `R_0` subtract the two copies, and let

\[
 F_h=M_{1,h}d_h
\]

be the assembled Whitney forward block.  Define

\[
 \widetilde F_h=F_hL_0,
\]

and

\[
 \widetilde W_h=
 \begin{pmatrix}
 0&\widetilde F_h^*\\
 \widetilde F_h&0
 \end{pmatrix}.
\]

The finite pencil is

\[
 \widetilde W_h+\kappa R_0^*R_0.
\]

Exact identities give

\[
 L_0J_0=I,
 \qquad
 R_0J_0=0,
 \qquad
 \widetilde F_hJ_0=F_h.
\]

Thus the weak conforming compression is exactly the known Whitney block.

The construction has endpoint-star support only.  An edge couples to the two
copies at each of its endpoints and no farther.  Before adding the positive
penalty, `tilde W_h` is exactly odd under form parity and Hermitian.

## Topology gate

For

\[
 N=8,16,
 \qquad
 \kappa\in\{1/2,1,2,4\},
\]

the complete generalized matrices give exactly two numerical zero modes:

1. the constant scalar;
2. the constant 1-form.

Both are also inserted explicitly and annihilated directly by the full
pencil.

The original no-flux operator, retained without modification as a negative
control, still has exactly one zero.  The repair therefore comes from the
consistency flux, not from a changed tolerance or reclassification.

At zero momentum the third branch is exactly

\[
 z_{\rm mismatch}=\frac{12\kappa}{h}.
\]

The two physical Betti modes remain at zero while this mismatch mode moves to
the cutoff.

## Known-answer spectral convergence

The first positive physical branch was evaluated at

\[
 N=8,16,32,64,128,256
\]

for every frozen `kappa`.  The target `2*pi` is known from the continuum unit
circle and was not used in the construction.

| `kappa` | value at `N=8` | value at `N=256` | final absolute error | initial/final error |
|---:|---:|---:|---:|---:|
| 0.5 | 6.352671 | 6.283340 | 1.551e-4 | 447.96 |
| 1 | 6.402651 | 6.283342 | 1.564e-4 | 763.78 |
| 2 | 6.424943 | 6.283342 | 1.571e-4 | 902.56 |
| 4 | 6.435490 | 6.283343 | 1.574e-4 | 967.72 |

All four error sequences decrease monotonically and improve by far more than
the preregistered factor-16 gate.

The final dimensionless velocities `z/(2*pi)` are

\[
 1.00002469,quad
 1.00002489,quad
 1.00002500,quad
 1.00002505.
\]

They approach one, not four different `kappa`-dependent speeds.

The consecutive error ratios approach four, consistent with the calibrated
second-order Whitney limit.  This is a diagnostic pattern; no exponent was
fitted or used as an acceptance gate.

For `N=8,16`, full-matrix first-positive eigenvalues agree with the selected
Bloch branch to maximum relative error `1.31e-14`.

## Attack on the positive framing

The construction starts from the already derived assembled Whitney forward
block `F_h`.  Therefore it does **not** independently select Kähler--Dirac
dynamics.  It supplies a canonical local lift of dynamics already chosen by
the Whitney complex.

The positive result means:

- the local duplicated carrier need not lose topology at finite stiffness;
- a consistency flux can replace the singular `kappa h -> infinity` route;
- the flux coefficient on the regular circle is fixed by copy symmetry.

It does not mean:

- that Kähler--Dirac itself has been derived from more primitive physics;
- that the same averaging is uniquely selected on a three-dimensional mesh
  with inequivalent incident tetrahedra;
- that the complete finite-`kappa` pencil is chiral—the penalty remains even;
- that a physical value of `kappa` or a time unit has appeared.

This distinction prevents a successful lifting construction from being
misreported as new dynamics.

## The three-dimensional canonicity problem

For a global simplex with several local copies, the direct analogue is an
averaging left inverse

\[
 L_pJ_p=I.
\]

Equal counting weights give

\[
 L_p=(J_p^*J_p)^{-1}J_p^*.
\]

This is combinatorially canonical, but three-dimensional incident copies need
not be metrically equivalent.  Other possibilities include exact
mass-weighted or face-trace-weighted averages.  The mass-orthogonal choice

\[
 (J_p^*M_pJ_p)^{-1}J_p^*M_p
\]

reintroduces the assembled inverse and is generally nonlocal.

Therefore the circle result authorizes a three-dimensional audit but does not
prejudge it.  The next acceptance boundary must ask whether geometry and
adjointness uniquely select a bounded-star `L_p`; simply choosing equal
weights would be too weak if the incident cells are inequivalent.

## Physical status

- **DERIVED:** unique symmetric weights `(1/2,1/2)` on the circle.
- **DERIVED:** exact locality, Hermiticity, odd flux, and conforming
  compression.
- **DERIVED:** both Betti modes survive every tested finite `kappa`.
- **DERIVED NUMERICAL:** all fixed-`kappa` low branches converge to `2*pi`.
- **DERIVED:** the mismatch branch is exactly `12 kappa/h`.
- **STRUCTURAL POSITIVE:** finite stiffness can be continuum-consistent when
  accompanied by the selected central flux.
- **STRUCTURAL LIMITATION:** the assembled Whitney operator is an input to
  the lift, not an output.
- **OPEN:** uniqueness and bounded-star locality of a three-dimensional flux.
- **OPEN:** finite-cutoff chirality, Lorentzian time, causal speed, mass, and
  Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_central_flux_circle.py
```

Expected result: `13/13`.

