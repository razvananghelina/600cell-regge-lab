# Preregistration: a local unitary tick from the signed Kähler--Dirac incidence

Date: 2026-08-11

## Motivation

The repository has a derived spatial Kähler--Dirac operator

\[
D=d+d^*
\]

on the 2640 oriented cochains of the 600-cell boundary.  Its accepted Whitney
dynamics is continuous-time and metric-unitary, but the time parameter and
the coefficient called \(c\) are external, and the strong generator need not
have finite support.  Older exploratory scripts that iterate \(D\) or
normalized adjacency are not unitary evolutions.

This protocol tests one sharply defined replacement: the signed
Grover--Szegedy walk canonically associated with the oriented Hasse incidence
graph of the same cochain complex.  The construction is standard quantum-walk
mathematics; finding it here is not by itself a new physical theory.

No value involving \(A_1=5\), \(\alpha\), a particle mass, the number of
generations, the Planck scale or experimental data will be inspected.

## Frozen carrier and operator

Reconstruct exactly:

1. the 600-cell boundary f-vector `(120,720,1200,600)`;
2. all signed coboundaries over the integers;
3. the signed Hasse adjacency \(D=d+d^T\);
4. the Hasse degree \(q_x\) of every oriented simplex \(x\).

The walk Hilbert space has one basis state \(|x,y\rangle\) for every directed
codimension-one incidence \(x\sim y\).  No self-loops or fitted weights are
added.

For each simplex, define the normalized incidence state

\[
|s_x\rangle=\frac1{\sqrt{q_x}}
\sum_{y\sim x}\eta_{xy}|x,y\rangle,
\]

where the arc phases are fixed by form degree:

\[
\eta_{xy}=D_{xy}\quad\text{if }\dim y=\dim x+1,
\qquad
\eta_{xy}=1\quad\text{if }\dim y=\dim x-1.
\]

Thus \(\eta_{xy}\eta_{yx}=D_{xy}\).  Let

\[
C=2\sum_x|s_x\rangle\langle s_x|-I,
\qquad
S|x,y\rangle=|y,x\rangle,
\qquad
U=SC.
\]

Changing simplex orientations is allowed only as a diagonal gauge change.
`CS` versus `SC` is a shift conjugacy, time reversal gives \(U^{-1}\), and a
global phase is not counted as a distinct physical operator.  No coin angle
is permitted.

## Tests frozen before spectral evaluation

The verifier must check:

1. **Exact complex:** the integer coboundaries square to zero and reproduce
   the certified Kähler--Dirac carrier.
2. **Exact unitarity:** every local coin block obeys \(C_x^2=I\), the shift
   obeys \(S^2=I\), hence \(U^*U=I\).  This must be certified algebraically,
   not inferred from a rounded norm.
3. **Strict locality:** every nonzero transition of one tick changes the tail
   simplex by exactly one Hasse incidence.  After \(n\) ticks, support is
   bounded by Hasse distance \(n\).
4. **Grading:** one micro-tick exchanges even and odd form degree exactly.
5. **Own-operator spectral map:** with

   \[
   Q=\operatorname{diag}(q_x),\qquad
   T=Q^{-1/2}D Q^{-1/2},
   \]

   the discriminant of the walk must equal \(T\).  Every
   \(T f=\lambda f\), \(|\lambda|<1\), must yield walk phases

   \[
   e^{\pm i\arccos\lambda}.
   \]

   The two topological zero modes of \(D\) must remain identifiable under
   this map.
6. **Known-answer calibration:** on an oriented cycle, the same construction
   must give ballistic left/right propagation at exactly one Hasse edge per
   micro-tick.  Peaks, fitted velocities and continuum regressions are not
   allowed.
7. **Nontriviality:** a localized state must evolve while preserving norm.

All reported equalities that involve square roots may use a numerical audit
only after the corresponding algebraic identity has been proved.

## Decision boundaries

- **DERIVED LOCAL UNITARY LIFT:** all seven tests pass.  This establishes a
  canonical dimensionless micro-tick and an exact finite propagation cone on
  an enlarged incidence-arc carrier.
- **KILL:** unitarity requires an adjustable angle, transitions leave the
  incidence neighbourhood, or the discriminant is not the normalized signed
  Kähler--Dirac operator.
- **OPEN EVEN AFTER PASSING:** uniqueness among all possible local unitary
  lifts, refinement compatibility on the 600-cell tower, Lorentzian
  reconstruction, selection of an internal mass operator, and conversion of
  the dimensionless edge/tick units to metres/seconds/kilograms.

In particular, a successful finite-speed walk does **not** derive

\[
t_P=\sqrt{\frac{\hbar G}{c^5}},\qquad
m_P=\sqrt{\frac{\hbar c}{G}}.
\]

Those require a dimensionful normalization involving \(\hbar\), \(G\) and
\(c\), none of which is selected by this finite combinatorial construction.

## Framing attack

The signed incidence data select this Grover--Szegedy lift within the stated
reflection construction, but they do not prove that Nature must quantize the
incidence operator by this functor.  Therefore even a perfect pass is a
mathematical bridge, not yet a selected fundamental dynamics.  The next
physical gate would be refinement stability of its low-quasienergy cone.

