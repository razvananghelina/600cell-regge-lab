# Protocol: unique copy-symmetric central flux on the circle

Date: 2026-08-11

This protocol is frozen before evaluating the new flux spectrum.  It tests
the minimal non-fitted repair of the finite-stiffness circle no-go.

## Fixed carrier and maps

Use the same unit circle with `N` equal edges, `h=1/N`, exact local Whitney
metric, and duplicated coefficients `(a_n,b_n,e_n)` as in protocol
`d499ca0`.

Let

- `J_0` copy a global vertex value to the two incident edge copies;
- `R_0` take their difference;
- `L_0` average them at each global vertex;
- `J_1=L_1=I` on edge 1-cochains.

Thus

\[
 (L_0x)_n=\frac{a_n+b_{n-1}}2,
 \qquad L_0J_0=I.
\]

The coefficient `1/2` is not fitted.  Among vertex-local linear maps invariant
under exchanging the two incident copies, the left-inverse condition forces
equal weights `w,w` with `2w=1`.

Let the assembled exact Whitney forward block be

\[
 F_h=M_{1,h}d_h.
\]

Define the central-flux forward block and odd weak operator by

\[
 \widetilde F_h=F_hL_0,
 \qquad
 \widetilde W_h=
 \begin{pmatrix}0&\widetilde F_h^*\\
 \widetilde F_h&0\end{pmatrix}.
\]

The tested pencil is

\[
 (\widetilde W_h+\kappa R_0^*R_0)v=zM_hv,
 \qquad \kappa>0.
\]

No other flux, upwind coefficient, mass lumping, or phenomenological target
enters the construction.

## Exact gates before continuum comparison

For `N=8,16` and `kappa=1/2,1,2,4`:

1. `L_0J_0=I` exactly;
2. `R_0J_0=0` exactly;
3. `tilde F_h J_0=F_h` exactly, so conforming weak compression is the
   assembled Whitney block;
4. `tilde W_h` is Hermitian and odd before the even penalty is added;
5. every matrix entry couples either copies at one vertex or an edge to its
   two endpoint stars; no longer-range coupling is allowed;
6. the complete generalized pencil has exactly two zero modes for every
   tested finite `kappa`;
7. those zeros are explicitly the constant scalar and constant 1-form.

The original local-without-flux pencil is retained as a negative control and
must still have only one zero.

## Known-answer spectral comparison

For fixed

\[
 \kappa\in\{1/2,1,2,4\}
\]

and

\[
 N\in\{8,16,32,64,128,256\},
\]

compute the `q=2*pi/N` Bloch block.  The smallest positive eigenvalue is
compared with the known first positive unit-circle Kähler--Dirac eigenvalue

\[
 2\pi.
\]

No convergence exponent is fitted.  Record every error and consecutive error
ratio.  The preregistered acceptance requires:

- the error tends monotonically downward over the frozen grid for every
  fixed `kappa`;
- the final error is smaller than the initial error by at least a factor 16;
- the dimensionless low-mode velocity `z/(2*pi)` tends toward one rather than
  a `kappa`-dependent constant.

The factor-16 gate is a coarse falsification threshold fixed before data, not
an estimate of order.

For `N=8,16`, the complete full-matrix first positive eigenvalue must agree
with the Bloch value within `1e-10` relative error.

## Stiff branch and chirality

At `q=0`, the mismatch branch is predicted to have exact eigenvalue

\[
 z_{\rm stiff}=\frac{12\kappa}{h},
\]

while both physical Betti modes remain zero.  Verify this exactly/numerically.

The flux operator `tilde W_h` is odd under form parity.  The positive penalty
is even, so the complete finite-`kappa` pencil is not odd.  Passing the
continuum test does not erase that finite-cutoff chirality defect.

## Acceptance and kill boundaries

**Acceptance on `S1`:** all exact gates pass and all four fixed-`kappa`
low branches converge to `2*pi` while preserving both Betti zeros.

This would establish only that a unique copy-symmetric consistency flux
repairs the known one-dimensional obstruction.  It would authorize, but not
prove, a three-dimensional analogue.

**Kill:** loss of either zero mode, a `kappa`-dependent limiting velocity,
nonlocal support, or failure of the frozen convergence gate.

No passing outcome selects `kappa`, a physical time unit, causal speed, mass,
or Planck scale.

