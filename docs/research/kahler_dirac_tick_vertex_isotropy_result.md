# The theory's own local tick has an isotropic vertex-centred causal front

Date: 2026-08-11  
Preregistration commit: `0015738`

## Result

Starting from the canonical embedded vertex cochain

\[
\psi_0=Ae_0,
\]

the signed Kähler--Dirac Grover--Szegedy tick is tangent-isotropic at every
one of the eight preregistered nonzero times.

> **DERIVED VERTEX-CENTRED ISOTROPIC TICK:** the hit fraction is `8/8`; the
> largest measured drift is `4.36e-17`, every covariance eigenvalue ratio is
> one to displayed precision, and the largest normalized traceless residual
> is `1.63e-16`.

The support also respects an exact angular speed bound

\[
c_{\rm angle}\leq\frac{\pi}{10}
\quad\text{radians per Hasse micro-tick},
\]

and its outer front numerically saturates that bound through tick eight at
the stored-coordinate precision.

This is the cleanest spacetime-kinematic result obtained in the present
route.  It is still dimensionless finite-complex kinematics, not yet a
Lorentzian continuum theory or the measured speed of light.

The targeted verifier passes `17/17` in about one second.  No full suite was
run.

## Why this succeeds where the flag-chamber walk failed

A complete flag chamber is a generic point of the Coxeter complex and has
trivial H4 stabilizer.  Its four neighbour directions are not a regular
tetrahedral frame, so the earlier four-direction scalar and fixed-coin tests
were strongly anisotropic.

A 600-cell vertex is different.  Its 12 graph neighbours form an
icosahedral vertex figure.  With equal weights their tangent shell has

\[
\mu=0,
\qquad
\Sigma=\frac13 I_3
\]

to numerical residual `5.66e-17`.  Doubling the probability of one frozen
direction makes the same estimator reject isotropy, so the positive is not
an estimator artefact.

The initial state `A e_0` is the signed equal-incidence embedding selected by
the normalized Kähler--Dirac construction.  The walk is H4-equivariant, so
the icosahedral stabilizer of the vertex preserves the evolved probability
distribution.  Its irreducible tangent action explains the zero vector mean
and scalar covariance.  The numerical calculation verifies that this
symmetry reasoning is realized by the actual signed incidence code.

Thus the two outcomes are not contradictory:

- **generic flag localization:** anisotropic;
- **canonical vertex-cochain localization:** exactly isotropic at second
  moment.

## Complete eight-tick record

| tick | occupied simplex tails | drift norm | covariance ratio | residual | RMS radius |
|---:|---:|---:|---:|---:|---:|
| 1 | 12 | `3.75e-17` | 1 | `2.66e-17` | 0.314159 |
| 2 | 43 | `4.34e-17` | 1 | `9.47e-17` | 0.294218 |
| 3 | 114 | `4.36e-17` | 1 | `9.14e-17` | 0.330446 |
| 4 | 255 | `1.74e-17` | 1 | `7.60e-17` | 0.439733 |
| 5 | 336 | `3.91e-17` | 1 | `1.63e-16` | 0.399634 |
| 6 | 687 | `2.70e-17` | 1 | `1.28e-16` | 0.421587 |
| 7 | 588 | `3.85e-17` | 1 | `1.49e-16` | 0.479310 |
| 8 | 1109 | `1.53e-17` | 1 | `5.48e-17` | 0.428743 |

All nine norms, including the initial state, remain one.  The form-degree
parity flips exactly at every tick, and no occupied simplex violates the
radius-`n` Hasse cone.

The nonmonotone RMS radius is important: the packet bulk is not simply a
ballistically translating shell.  The clean velocity statement concerns
the outer causal front, not the mean or RMS packet speed.

## The angular causal bound

Adjacent 600-cell vertices have spherical separation

\[
\arccos(\varphi/2)=\frac{\pi}{5}.
\]

The centre of an incident edge is halfway along that geodesic, so a
vertex--edge Hasse step has exact length

\[
\ell_{01}=\frac{\pi}{10}.
\]

The three numerically reconstructed rank-pair lengths are

\[
(\ell_{01},\ell_{12},\ell_{23})
=
(0.3141592654,\ 0.1887105308,\ 0.1354592602).
\]

Therefore `pi/10` is the largest possible distance crossed by one local
micro-tick.  Locality plus the spherical triangle inequality gives the exact
front bound

\[
r_{\max}(n)\leq n\frac{\pi}{10}.
\]

For `n=1,...,8`, the computed outer support satisfies

\[
r_{\max}(n)=n\frac{\pi}{10}
\]

within maximum stored-coordinate residual `1.34e-10`.  This saturation is
**DERIVED NUMERICAL**; the one-step upper bound itself is exact.

If the unit three-sphere later acquires physical radius `R_*` and one
micro-tick acquires duration `tau_*`, the bound becomes

\[
c_{\max}=\frac{\pi R_*}{10\tau_*}.
\]

Nothing currently selects `R_*` or `tau_*`.  Setting this expression equal
to the measured `c` would define their ratio, not predict it.

## It really uses the theory's operator sector

The verifier independently reconstructs

\[
A^*SA=Q^{-1/2}DQ^{-1/2}
\]

and both invariant-plane identities

\[
UA=SA,
\]

\[
USA=2SAT-A.
\]

Hence the chosen state and all its iterates remain in

\[
\operatorname{span}(A\mathcal H,SA\mathcal H),
\]

the 5,280-dimensional sector spectrally tied to the original 2,640-component
Kähler--Dirac operator.  The result is not produced by leaking into the
9,600 walk-only complement.

## What is and is not physical yet

What is now established:

- exact unitary reversible evolution;
- strict locality and a finite angular causal cone;
- a canonical own-operator initial state;
- exact second-moment isotropy about every equivalent 600-cell vertex;
- numerical saturation of the causal front through the frozen horizon.

What remains missing:

- a refinement limit preserving this cone;
- a long-wavelength quasienergy relation with Dirac/Lorentz form;
- proof that the relevant excitations occupy the selected invariant sector;
- one physical conversion from radius and tick to metres and seconds;
- a mass term and an explanation of inertia;
- gravitational normalization and therefore Planck units.

The three Hasse rank transitions also have unequal geometric lengths.  A
single graph tick is therefore not one common microscopic ruler step, even
though `pi/10` is a universal maximum.  This distinction prevents the
present bound from being advertised as a complete derivation of `c`.

## Status ledger

- **DERIVED:** exact unitary Kähler--Dirac local tick and invariant-sector
  identities.
- **DERIVED CONTROL:** the 12-direction icosahedral shell is isotropic, while
  a frozen one-direction perturbation is rejected.
- **DERIVED NUMERICAL:** vertex-centred isotropy hits `8/8`.
- **STRUCTURAL EXPLANATION:** the A5 vertex stabilizer forces the observed
  first two tangent moments.
- **DERIVED:** exact angular causal upper speed `pi/10` per micro-tick.
- **DERIVED NUMERICAL:** the front saturates that bound through tick eight.
- **OPEN:** Dirac dispersion and controlled refinement.
- **OPEN:** physical `R_*`, `tau_*`, mass, inertia and Planck normalization.
- **NOT CLAIMED:** observed SI `c` or Lorentz invariance.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_kahler_dirac_tick_vertex_isotropy.py
```

Expected result: `17/17`.
