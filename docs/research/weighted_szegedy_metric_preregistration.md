# Preregistration: can arbitrary incidence weights repair the metric tick?

Date: 2026-08-11

## Fixed prior facts

Commit `cf7aefa` established two exact facts:

1. the unweighted signed Grover--Szegedy walk refines exactly on the
   homogeneous circle;
2. on one barycentric flag tetrahedron, the consistent Whitney metric
   codifferentials in degrees 0 and 1 contain 10 and 12 nonzero coefficient
   entries outside ordinary face/coface incidence.

This protocol tests the most obvious rescue: retain the same directed Hasse
arcs but choose arbitrary local amplitudes.

No physical target, spectrum, mass, coupling, `A1=5` or Planck quantity will
be inspected.

## Complete class under test

Let \(K\) be the full simplex lattice of one barycentric flag tetrahedron.
The coefficient space has one orthonormal state \(|x\rangle\) for every
simplex.  The walk carrier has exactly one state \(|x,y\rangle\) for every
directed codimension-one incidence \(x\sim y\), and no other arcs or ancillas.

For every simplex permit an arbitrary normalized complex vector

\[
|s_x\rangle=\sum_{y\sim x}a_{xy}|x,y\rangle,
\qquad
\sum_{y\sim x}|a_{xy}|^2=1.
\]

The coefficients may be unequal, complex or zero.  Permit arbitrary
unit-modulus reversal phases in

\[
S|x,y\rangle=\omega_{xy}|y,x\rangle,
\qquad
S^2=I.
\]

Define

\[
C=2\sum_x|s_x\rangle\langle s_x|-I,
\qquad U=SC,
\qquad A|x\rangle=|s_x\rangle.
\]

This includes the previously committed signed uniform walk and every
rank-one weighted Grover/Szegedy coin on the same incidence carrier.  A
nonorthogonal Whitney encoding, extra local copies, added arcs or a multi-step
effective walk is outside this class and is not silently ruled out.

## Universal support theorem to test

The discriminant is

\[
T=A^*SA.
\]

Enumerate its symbolic support before substituting any amplitudes.  Because
the outgoing arc sets of distinct simplices are disjoint and the shift only
reverses an existing arc, the proposed theorem is

\[
T_{xy}=0\quad\text{whenever }x\not\sim y,
\]

for every permitted choice of \(a_{xy}\) and \(\omega_{xy}\).

Diagonal basis gauges and a global scalar do not alter this zero pattern.

## Metric comparison

Independently integrate the exact Whitney mass matrices of the same flag
tetrahedron and calculate

\[
\delta_k=M_k^{-1}d_k^TM_{k+1},\qquad k=0,1,2.
\]

Record every nonzero entry of \(\delta_k\) outside the support of \(d_k^T\).

## Decision boundary

- **DERIVED NO-GO FOR WEIGHTED INCIDENCE COINS:** at least one exact Whitney
  entry lies off incidence while the universal discriminant support theorem
  holds.  Then no adjustment of weights or phases on the existing arcs can
  repair the metric mismatch.
- **RESCUE REMAINS OPEN IN THIS CLASS:** all metric adjoints have incidence
  support, or the general discriminant can connect a nonincident pair.

This is not a no-go for all local unitary dynamics.  A negative result says
precisely that the next construction must enlarge the carrier, introduce
intermediate local steps, or use a nonorthogonal metric encoding.  Each such
extension needs its own canonicity and refinement test.

