# Positive multiplier kinetics immediately changes the Whitney spectrum

Date: 2026-08-11

Preregistration commit: `c34c659`

Targeted verifier:
`reproducible/verify_whitney_positive_metric_completion.py`

Targeted result: **8/8 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Complete result

Start with the exact constrained descriptor

\[
K=\begin{pmatrix}A_{\mathrm{loc}}&C^*\\C&0\end{pmatrix},
\qquad
B_0=\begin{pmatrix}M_{\mathrm{loc}}&0\\0&0\end{pmatrix}.
\]

Give the independent multipliers an arbitrary positive metric (N>0):

\[
B_\varepsilon=
\begin{pmatrix}M_{\mathrm{loc}}&0\\0&\varepsilon N\end{pmatrix},
\qquad\varepsilon>0.
\]

This turns the descriptor into an ordinary positive-metric generalized
Hamiltonian, but it cannot preserve the exact Whitney spectrum.

> **DERIVED POSITIVE-METRIC NO-GO FOR THE MINIMAL COMPLETION:** every positive
> multiplier metric makes the second-class sector dynamical, adds (2r)
> finite spectral slots relative to the physical descriptor, and generically
> moves the original nonzero modes immediately.

The claim covers the block-fixed family above for every (N>0).  It does not
cover new Hamiltonian blocks or a different physical embedding.

## Exact dimension obstruction

For local dimension (n) and independent constraint rank (r), the singular
descriptor has (n-r) physical finite eigenvalues.  A positive metric on the
full ​((n+r))-dimensional carrier has (n+r) finite eigenvalues.  It therefore
adds exactly

\[
(n+r)-(n-r)=2r
\]

finite slots.

| level | (n) | (r) | physical descriptor count | positive-metric count | extra finite slots |
|---|---:|---:|---:|---:|---:|
| base | 9,000 | 6,360 | 2,640 | 15,360 | 12,720 |
| first barycentric | 216,000 | 153,120 | 62,880 | 369,120 | 306,240 |

The extra count equals the real second-class constraint count.  Calling those
modes unphysical would require an additional exact projection or gauge
principle, neither of which exists in this completion.

## First-order eigenvalue obstruction

For a simple finite descriptor eigenpair (y=(u,\lambda)), symmetry gives

\[
z'(0)=-z\,
\frac{\lambda^*N\lambda}{u^*M_{\mathrm{loc}}u}.
\]

Since (N>0), every nonzero mode with ​(lambda\ne0) shifts.  In a degenerate
eigenspace the splitting matrix is

\[
-zL^*NL,
\]

where the columns of (L) are its multiplier components.  The whole
eigenspace remains fixed only if (L=0).

## Independent control

On the five-tetrahedron boundary of a 4-simplex:

\[
n=75,qquad r=45,qquad n-r=30,qquad n+r=120.
\]

The descriptor has six nonzero eigenvalue clusters containing 28 branches,
plus two harmonic zero modes.  With the preregistered (N=I) control:

- nonzero eigenspace shift hit fraction: **6/6**;
- nonzero branch shift hit fraction: **28/28**;
- multiplier-equation residual: (2.20\times10^{-15});
- perturbative/finite-difference normalized residual at
  ​(arepsilon=10^{-6}): (1.36\times10^{-5});
- additional finite modes: exactly **90**.

The nonzero descriptor eigenvalues and their first-order derivative ranges
are:

| eigenvalue | multiplicity | derivative range |
|---:|---:|---:|
| -2.73861 | 4 | 0.867 to 1.141 |
| -2.58199 | 6 | 0.277 to 0.822 |
| -1.93649 | 4 | 0.161 to 0.635 |
| +1.93649 | 4 | -0.635 to -0.161 |
| +2.58199 | 6 | -0.822 to -0.277 |
| +2.73861 | 4 | -1.141 to -0.867 |

At ​(arepsilon=10^{-6}), the 90 extra modes have absolute eigenvalues from
approximately 1,378.5 to 6,697.6.  Their large scale is the expected singular
decoupling behaviour, not a derived physical mass.

## Why this does not derive the Planck scale

One could choose a tiny ​(arepsilon) and call the resulting heavy multiplier
modes a cutoff or Planck sector.  That would be fitting.  Neither the local
Whitney geometry nor the constraint kernel selects ​(arepsilon), and changing
it changes both the heavy scale and the physical eigenvalue shifts.

Only the singular limit

\[
\varepsilon\to0
\]

recovers the exact descriptor spectrum and sends the additional modes to
infinite frequency.  Thus the apparent Planck-like separation is exactly the
unselected infinity the programme was trying to avoid.

## What is now closed and open

- **DERIVED:** every positive multiplier metric adds (2r) finite slots.
- **DERIVED:** the simple-mode derivative formula for every (N>0).
- **DERIVED NUMERICAL CONTROL:** all 28 nonzero small-control branches shift.
- **DERIVED NEGATIVE:** the block-fixed positive-metric completion cannot be
  the exact Whitney tick.
- **DERIVED NEGATIVE:** a small multiplier kinetic coefficient does not derive
  a Planck scale; it inserts one.
- **OPEN:** a new first-class extension with additional Hamiltonian blocks.
- **OPEN:** a geometry-selected embedded physical subspace in a larger
  positive-metric carrier.
- **OPEN:** approximate finite-stiffness dynamics if exact Whitney spectra are
  relaxed.
- **NOT CLAIMED:** physical mass, time, (c), ​(hbar) or (G).

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_positive_metric_completion.py
```

Expected result: `8/8`.
