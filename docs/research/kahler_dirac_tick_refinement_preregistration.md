# Preregistration: refinement gate for the local Kähler--Dirac tick

Date: 2026-08-11

## Prior result

Commit `3783c52`, evaluated after target-blind protocol commit `ff6b2ce`,
constructs a signed Grover--Szegedy walk on the oriented Hasse incidences of
the fixed 600-cell cochain complex.  It is exactly unitary and propagates at
most one Hasse edge per micro-tick.

That result is only a finite combinatorial kinematics.  This protocol asks
whether it is compatible with geometric refinement and with the metric
Whitney Kähler--Dirac operator already accepted by the theory.

No value involving `A1=5`, a particle target, a coupling, a physical mass,
the Planck scale or experimental data will be examined.

## Part A: known-answer dyadic circle

Use a unit circle divided into \(N\) oriented edges and its dyadic refinement
to \(2N\) edges.  The Hasse graph has \(2N\) vertices and is degree two.

For Fourier mode \(m\), derive rather than fit:

1. the positive singular value of the normalized signed incidence;
2. the walk quasienergy measured relative to the central phase
   \(\pi/2\);
3. the relation between one coarse physical interval and two fine
   micro-ticks.

The acceptance condition is exact for every resolved mode below Nyquist:

\[
\varepsilon_N(m)=2\varepsilon_{2N}(m).
\]

The corresponding physical frequency, with Hasse-edge length
\(1/(2N)\), must be \(2\pi m\), so the calibrated dimensionless speed is
exactly one.  No peak, slope regression or adjustable normalization is
allowed.

Failure on this control kills the refinement interpretation immediately.

## Part B: first barycentric tetrahedron

Use the same regular reference tetrahedron and exact barycentric flag child
as `verify_whitney_kahler_induction.py`.  Recompute over the rationals:

1. the four face areas of one flag child;
2. its exact Whitney 2-form and 3-form mass matrices;
3. the metric codifferential

   \[
   \delta_2=M_2^{-1}d_2^T M_3;
   \]

4. the top-to-face vector used by the uniform signed Grover coin.

The comparison is qualitative and exact, not fitted:

- if the Whitney codifferential is proportional to the uniform signed
  incidence vector, the local tick survives this metric gate;
- if it is not proportional, the unweighted tick is **not** a unitary lift of
  the accepted metric Whitney dynamics, even though it remains a valid walk
  on the abstract Hasse graph.

Also record whether the consistent Whitney mass or its inverse is diagonal.
Non-diagonal metric data explain why replacing the uniform coin by local
scalar weights is not automatically equivalent to the accepted operator.

## Decision boundary

- **DERIVED REFINEMENT-COMPATIBLE TICK:** the circle scaling passes and the
  tetrahedral metric direction agrees up to a geometry-fixed scalar.
- **DERIVED KINEMATIC ONLY:** the circle passes but the tetrahedral Whitney
  direction is not proportional.  Then the exact cone belongs to the
  combinatorial graph metric, not yet to the accepted piecewise-Euclidean
  spectral geometry.
- **KILL:** the known-answer circle scaling fails.

A tetrahedral failure does not prove that every possible metric-aware local
unitary dilation is impossible.  It kills only the unweighted construction
committed in `3783c52`.  A weighted replacement would need a new blind
protocol and a proof of canonicity; arbitrary edge amplitudes are forbidden.

