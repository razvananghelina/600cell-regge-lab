# All three Coxeter bonds repair support connectivity

Date: 2026-08-11  
Preregistration commit: `9bf2fba`

## Result

Replacing the single published colour pair by the three consecutive bonds of
the (H_4) diagram repairs the exact support disconnection.

> **DERIVED CONNECTED SUPPORT SCAFFOLD:** the preregistered rank-forward
> schedule
>
> \[
> (01)\longrightarrow(12)\longrightarrow(23)
> \]
>
> is strongly connected on the 172,800-node three-phase active carrier.

Every one of the six possible orders is also strongly connected:

\[
\text{strong-connectivity hits}=6/6.
\]

Therefore connectedness validates the use of all three bonds but has no
power to select their temporal order.

The targeted verifier passes `10/10` in about 1.6 seconds.  No full suite was
run.

## Construction

The three rank-oriented bonds are

\[
B_{01},\qquad B_{12},\qquad B_{23}.
\]

For each (a<b), the corrected robust translation template acts on the four
active components by

\[
(T_{ab}\phi)(k)=
\begin{pmatrix}
\phi_2(k)\\
\phi_3(s_a k)\\
\phi_0(s_b k)\\
\phi_1(k)
\end{pmatrix}.
\]

All three maps are exact permutations of the 57,600 active states and cross
at most one chamber facet.  Between them, the verifier uses the exact
nonzero pattern of the published fixed coin

\[
\widehat C=I_2\otimes C,
\]

whose two (2\times2) blocks are dense.  No angle or coefficient was
optimized.

## Complete attempt count

The three distinct bonds have (3!=6) temporal orders.  All were enumerated
before reading the outcome.  Cyclic rotations of a schedule merely move the
origin of the three-phase clock, leaving two classes of three schedules each.
No reversal equivalence was assumed.

For every schedule, the exact directed graph on

\[
(\text{phase},\text{chamber},\text{component})
\]

has one strong component and one weak component.  Thus the hit fractions are

\[
6/6\quad\text{and}\quad6/6.
\]

This is not a look-elsewhere success: there was no rare hit.  It is the
stronger but less selective statement that any order using all three bonds
removes coordinate-support blocks.

## What changed relative to the literal walk

A single bond sees only one rank-two residue of the chamber geometry:

\[
\begin{array}{c|c}
\text{bond}&\text{chamber orbits}\\ \hline
(01)&2400\times6\\
(12)&2400\times6\\
(23)&1440\times10
\end{array}
\]

The literal robust walk used only the final line and was trapped in 1,440
decagons.  The three-bond schedule uses every colour (s_0,s_1,s_2,s_3), so
successive coin mixing permits transitions between all chamber-component
states.

## Hostile interpretation audit

This is only **support connectivity**.  It proves there is no invariant
partition spanned by subsets of coordinate basis states.  It does not prove:

- absence of nonlocal invariant subspaces;
- absence of exact destructive interference in powers of the Floquet map;
- an isotropic continuum limit;
- convergence to the Dirac equation;
- equality with the Whitney Kähler--Dirac operator.

Moreover, the three bonds are not microscopically equivalent: two generate
hexagons and the golden (H_4) bond generates decagons.  Calling them the
three physical tangent axes is therefore **STRUCTURAL**, not DERIVED.  The
fixed rank order is canonical as combinatorics, but physics still has to show
that its anisotropy washes out under refinement.

## Status ledger

- **DERIVED:** all three bond translations are local permutations.
- **DERIVED:** the designated rank-forward schedule is strongly
  support-connected.
- **DERIVED:** all six orders pass; connectivity selects none of them.
- **DERIVED:** the microscopic orbit data are (6,6,10), not isotropic.
- **STRUCTURAL:** three Coxeter bonds are a canonical combinatorial analogue
  of three directions.
- **OPEN:** exact amplitude-level invariant sectors and cancellations.
- **OPEN:** a target-blind isotropy/refinement gate.
- **OPEN:** relation to the theory's Whitney/Kähler--Dirac operator.
- **NOT CLAIMED:** a Dirac particle, mass, inertia, (c), or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_h4_three_bond_walk.py
```

Expected result: `10/10`.
