# Preregistration: leading angular multipole of the Kähler--Dirac tick

Date: 2026-08-11

## Status and provenance

This is a **confirmatory, not blind**, test.  A colleague reported that a
quantity denoted (A_6) is much smaller at even than at odd ticks.  The
reported normalization and exact definition are not present in the
repository, so that numerical range is not yet reproducible.  No value from
the existing tick evolution has been inspected under the definitions below
before this protocol is committed.

The geometric input, walk, initial state and ticks (n=0,\ldots,8) are fixed
by protocol commit `0015738` and verifier
`reproducible/verify_kahler_dirac_tick_vertex_isotropy.py`:

- the signed Grover--Szegedy Kähler--Dirac tick on directed Hasse arcs;
- initial state (A e_0) at vertex 0;
- probabilities aggregated at the 5,280 simplices;
- the Riemannian logarithm at vertex 0 as the tangent coordinate.

No Standard-Model or experimental target enters this calculation.

## Coordinate-free definition

For each occupied simplex away from the base point, write its tangent
coordinate as (x_i=r_i u_i), with (r_i>0) and (u_i\in S^2), and let
(p_i) be its probability.  For any nonnegative weights (q_i), define

\[
 A_\ell(q)^2 =
 \frac{\sum_{i,j} q_i q_j P_\ell(u_i\mathbin{\cdot}u_j)}
      {(\sum_i q_i)^2}.
\]

This is the rotation-invariant norm of the spherical-harmonic coefficient at
degree \(\ell\), normalized to one for a point mass.  Tiny negative squared
values in ([-10^{-12},0)) are rounded to zero; a value below (-10^{-12})
fails the verifier.

Three natural weightings are frozen because the colleague's normalization is
unknown and normalization dependence is itself a possible artefact:

1. **conditional angular:** (q_i=p_i), normalized after deleting the base
   point;
2. **unconditional angular:** the same numerator but divided by one rather
   than by the moving mass, so stay-put probability can suppress it;
3. **solid-harmonic radial:** (q_i=p_i r_i^\ell), normalized by
   \(\sum_i p_i r_i^\ell\).

The conditional quantity measures angular shape; the unconditional quantity
also measures how much probability has left the origin; the radial quantity
tests cancellation of the degree-ℓ solid harmonic across shells.  None may
be selected after seeing which one gives the strongest parity effect.

## Frozen computations and controls

For every tick (1,\ldots,8), compute all three quantities for
(\ell=1,\ldots,12).  Degree six is the preregistered primary observable;
the other degrees are diagnostics, not 12 independent searches.

Independently compute the multiplicity of the trivial representation of the
icosahedral rotation group in the (SO(3)) harmonic of degree \(\ell):

\[
 m_\ell=\frac1{60}\left[(2\ell+1)+15\chi_\ell(\pi)
 +20\chi_\ell(2\pi/3)+12\chi_\ell(2\pi/5)
 +12\chi_\ell(4\pi/5)\right],
\]

where
(\chi_\ell(\theta)=\sin((\ell+\tfrac12)\theta)/\sin(\theta/2)).
The control must show no nonconstant invariant for
(1\leq\ell\leq5) and the first one at \(\ell=6\).  The uniform 12-neighbour
icosahedral shell must therefore have (A_1,\ldots,A_5<10^{-10}) and
(A_6>10^{-6}).  A deliberately distorted shell must generate at least one
lower multipole above (10^{-3}).

Because the entire evolved distribution and each radius are invariant under
the vertex stabilizer, all three variants must have
(A_1,\ldots,A_5<10^{-9}) at every nonzero tick.  This is an implementation
gate, not the sought result.

## Frozen parity diagnostics

Print the eight (A_6) values for each weighting, plus

\[
 R_{\rm sep}=\frac{\max_{n\in\{2,4,6,8\}} A_6(n)}
                    {\min_{n\in\{1,3,5,7\}} A_6(n)}.
\]

Classify each weighting before any physical interpretation:

- **DERIVED parity separation** if (R_{\rm sep}<0.1);
- **PATTERN only** if the even values are all below the odd values but
  (R_{\rm sep}\geq0.1), or if separation occurs only for some weightings;
- **DERIVED NEGATIVE** for that weighting if the even and odd ranges overlap.

Even a robust separation is not by itself Lorentz invariance or a continuum
limit.  Nonmonotonicity within the even subsequence, dependence on weighting,
and the support size at each tick must be reported.  A claim of physical
anisotropy suppression remains **OPEN** until it survives longer times and a
refinement family selected independently of this observation.

Only the targeted verifier will be run; the full suite is deliberately out of
scope at the user's request.
