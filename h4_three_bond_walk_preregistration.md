# Preregistration: the three Coxeter bonds as a connected walk schedule

Date: 2026-08-11

## Motivation fixed before computation

Commit `9de718e` closed the literal two-colour robust transplant: arbitrary
local coins cannot connect its 1,440 decagonal orbit blocks.  The same audit
showed that colours 0 and 1 connect the quotient of those blocks.

The labelled (H_4) Coxeter diagram is the path

\[
0\xleftrightarrow{3}1\xleftrightarrow{3}2
\xleftrightarrow{5}3.
\]

It contains exactly three adjacent bonds.  Before testing any continuum
target, use all three as the only candidate translation directions:

\[
B_{01},\qquad B_{12},\qquad B_{23}.
\]

## Frozen construction

For an ordered colour pair (a<b), define (T_{ab}) by applying the
corrected robust macro template from commit `69447b9`, substituting
(s_a,s_b) for the paper's (s_2,s_3).  On the four active components,

\[
(T_{ab}\phi)(k)=
\begin{pmatrix}
\phi_2(k)\\
\phi_3(s_a k)\\
\phi_0(s_b k)\\
\phi_1(k)
\end{pmatrix}.
\]

The bond orientation (a<b) is fixed by flag rank.  Reversing individual
bonds is not allowed in this protocol.

Between translations, use only the nonzero support of the paper's fixed
directional coin

\[
\widehat C=I_2\otimes C,
\qquad
C=e^{i\pi/3}R_{\sigma_z}(\pi/2)R_{\sigma_x}(\pi/2).
\]

Both (2\times2) blocks of this coin are dense.  Connectivity depends only
on that exact zero pattern; no angle is optimized.

## Complete look-elsewhere set

There are (3!=6) orders of the three distinct bonds.  Enumerate all six.
Do not inspect one order, declare success and ignore the others.

The rank-forward schedule

\[
(01)\to(12)\to(23)
\]

is designated before execution.  Cyclic rotations correspond to a change of
the origin of the three-phase clock; nevertheless all six are retained in
the reported count.  Reverse order is not assumed equivalent and must be
tested.

Thus the preregistered attempt count is

\[
N=6.
\]

## Exact tests

1. Reconstruct the four-coloured 14,400-chamber (H_4) graph.
2. Check each (T_{01},T_{12},T_{23}) is a permutation on the 57,600 active
   states and crosses at most one chamber edge.
3. Record the single-bond chamber component multisets; no equality is
   assumed in advance beyond the Coxeter relations.
4. For every one of the six schedules, construct the exact nonzero-support
   graph of the periodic three-phase evolution

   \[
   \widehat C T_{b_3}\;\widehat C T_{b_2}\;
   \widehat C T_{b_1}.
   \]

5. Compute weak and strong connected components of the time-expanded graph
   `(phase, chamber, active component)` with 172,800 vertices.
6. Report the hit fraction: how many of the six schedules are strongly
   connected, and separately whether the designated rank-forward schedule
   is strongly connected.
7. Record the number of genuine schedule classes after quotienting only by
   cyclic phase-origin rotations.  Do not quotient by reversal unless an
   explicit conjugacy is constructed.

## Decision boundaries

- **DERIVED CONNECTED SCAFFOLD:** the designated schedule is strongly
  connected on the time-expanded active carrier.
- **DERIVED NEGATIVE:** it has more than one strong component.

Even a positive result is **STRUCTURAL**, not a Dirac derivation.  If all six
schedules pass, connectedness has no selecting power and that must be said
plainly.  If only a subset passes, the exact hit fraction is the evidence.

## Hostile framing boundary

Calling the three Coxeter bonds "spatial axes" is not yet derived.  They are
three canonical combinatorial directions in flag space, but one has Coxeter
order five while the other two have order three.  This anisotropy may prevent
an isotropic continuum limit.  The present test asks only whether the
geometry-led schedule repairs global propagation without free labels.

No Standard-Model target, mass, speed, Planck scale or desired dispersion is
used.
