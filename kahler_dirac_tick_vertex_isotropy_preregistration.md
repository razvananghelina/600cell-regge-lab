# Preregistration: vertex-centred isotropy of the theory's own local tick

Date: 2026-08-11

## Why this test replaces the chamber-local question

The four-colour flag chamber is a useful walk carrier, but a chamber has
trivial H4 stabilizer.  Exact rotational isotropy at one generic chamber is
therefore not enforced by the geometry.  By contrast, every 600-cell vertex
has an icosahedral tangent shell of 12 neighbours and all 120 vertices are
equivalent.

Before adding another coin, test the local unitary already derived from the
theory's own signed Kähler--Dirac operator:

\[
U=SC
\]

on the 14,880 directed Hasse incidences from commit `3783c52`.

No new coin, weight, phase or scale is allowed.

## Frozen initial state

Use vertex index 0 only, fixed before the computation.  Let `A` be the
certified signed normalized incidence embedding for which

\[
A^*SA=Q^{-1/2}DQ^{-1/2}.
\]

Start from

\[
\psi_0=Ae_0.
\]

This is the canonical embedding of a localized vertex cochain, not one
selected outgoing arc and not an optimized internal spin state.  It lies in
the previously certified 5,280-dimensional invariant own-operator sector.

Evolve exactly for the fixed times

\[
n=0,1,\ldots,8.
\]

The horizon eight is frozen before inspecting any moment result.  No time is
discarded because it looks bad.

## Frozen geometric observable

Assign every simplex its normalized spherical barycentre in `S^3`.  At each
tick, sum `|psi(x,y)|^2` over all outgoing arcs with the same tail simplex
`x`.  This gives a probability distribution on the 2,640 simplex centres.

Map every occupied centre to the tangent space at the initial 600-cell
vertex using the spherical logarithm.  Record at every tick:

- probability normalization;
- support size and form-degree probability distribution;
- mean displacement and its norm;
- centred tangent covariance eigenvalues;
- eigenvalue ratio and normalized traceless residual;
- mean, root-mean-square and maximum radial geodesic distance.

For every `n>=1` with nonzero covariance, the frozen isotropy gate is

\[
\|\mu_n\|<10^{-11},\qquad
|\lambda_{\max}/\lambda_{\min}-1|<10^{-10},
\]

and normalized traceless residual below `10^-10`.

## Calibration

Before using the evolved distribution:

1. verify that the 12 geodesic directions from vertex 0 to its graph
   neighbours have equal lengths, zero mean and covariance proportional to
   `I_3`;
2. perturb one of their equal probabilities by a frozen factor of two and
   renormalize; the same estimator must reject isotropy.

This distinguishes a symmetry result from an estimator that always returns
one.

## Exact structural checks

1. Reconstruct the oriented cochain complex and verify `d^2=0` over the
   integers.
2. Reconstruct the signed Grover reflections from their integer numerators
   and verify each block squares exactly.
3. Verify `C`, `S` and `U` preserve norm at every recorded tick.
4. Verify the state remains inside
   `span(A H, S A H)` by checking the exact walk/discriminant intertwining
   relation already used to define that sector; do not infer this merely
   from norm preservation.
5. Verify the strict Hasse cone: every occupied tail at tick `n` is at Hasse
   distance at most `n` from the initial vertex.

## Decision boundaries

- **DERIVED VERTEX-CENTRED ISOTROPIC TICK:** all eight nonzero-time moment
  audits pass.
- **DERIVED NEGATIVE:** at least one preregistered tick fails; report the full
  hit fraction and do not keep only a favourable time.

## Hostile scope boundary

A positive result would show exact second-moment isotropy for a canonical
vertex-localized state on the fixed 600-cell.  It would explain why the
generic flag-chamber anisotropy test was too local for a vertex observer.

It would not establish:

- a continuum or refinement limit;
- a Dirac rather than diffusion-like dispersion;
- one common physical distance per Hasse micro-tick;
- Lorentz invariance, physical `c`, mass, inertia or Planck units;
- that arbitrary initial arc states are isotropic.

The vertex stabilizer may make the positive result symmetry-forced.  That is
still useful, but it is a kinematic theorem rather than dynamical selection.
