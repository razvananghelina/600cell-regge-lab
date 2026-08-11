# The direct tetrahedral Dirac-walk transplant fails on barycentric chambers

Date: 2026-08-11  
Preregistration commit: `e731871`

## Result

The first barycentric subdivision of the 600-cell supplies a remarkably
clean ordered chamber carrier, but the economical four-amplitude shift of
Nzongani *et al.* does not transplant to its natural handedness.

> **DERIVED POSITIVE, CARRIER:** the 14,400 complete flags form the exact
> four-coloured (H_4) chamber graph.  It is connected and has a canonical
> balanced orientation bipartition (7,200+7,200), up to global exchange.

> **DERIVED NEGATIVE, DIRECT DYNAMICS:** with that handedness, the published
> non-robust black shift is two-to-one rather than a permutation.  Therefore
> it is not unitary, for either global orientation convention.

The targeted verifier passes `12/12` in about 0.6 seconds.  No full suite was
run.

## Exact chamber geometry

Starting from the certified 600-cell complex

\[
(f_0,f_1,f_2,f_3)=(120,720,1200,600),
\]

the verifier enumerates all 24 flags inside each tetrahedron:

\[
v\subset e\subset f\subset t,
\qquad 600\times24=14{,}400.
\]

Changing one rank defines four fixed-point-free involutions
(s_0,s_1,s_2,s_3).  Their product orders are exactly

\[
\begin{pmatrix}
1&3&2&2\\
3&1&3&2\\
2&3&1&5\\
2&2&5&1
\end{pmatrix},
\]

the Coxeter relations of (H_4).  Thus the geometry really does provide
four ordered neighbours, without coordinates or fitted labels.

The chamber graph is connected and bipartite.  Its two sheets contain 7,200
chambers each and are exchanged across every colour.  This is the intrinsic
left/right datum available from chamber orientation.

## Why the published shift fails

The non-robust walk in [*Dirac quantum walk on tetrahedra*](https://arxiv.org/abs/2404.09840)
uses a handedness-dependent causal stage (S_B) followed by a local stage
(S_G).  The exact source count shows:

- components 0 and 2 are covered once iff the handedness label is invariant
  under (s_3);
- components 1 and 3 are covered once iff it is invariant under (s_2).

Hence (S_B) is a permutation iff handedness is constant on every orbit of

\[
\langle s_2,s_3\rangle.
\]

There are exactly 1,440 such orbits, each of size ten.  Natural chamber
orientation does the opposite: both (s_2) and (s_3) exchange its two
sheets.  On the 57,600-dimensional four-amplitude carrier the result is
exact:

\[
28{,}800\ \text{distinct inputs used twice},
\qquad
28{,}800\ \text{inputs never used}.
\]

The local grey stage remains a permutation, but its product with a
non-bijective black stage cannot be unitary.

## Can the failure be repaired by relabelling?

Algebraically, yes; canonically, no result has selected one.

Any binary value chosen independently on each of the 1,440 causal orbits
makes the necessary invariance possible, giving (2^{1440}) labelings.
Balanced labelings alone number

\[
\binom{1440}{720}.
\]

The full (H_4) chamber action is transitive, so a label invariant under the
complete geometry is constant, not a nontrivial left/right split.  The
paper's flat cube decomposition supplies a different chirality pattern in
which its chosen causal links preserve handedness.  Importing such a pattern
into the 600-cell without another selector would be fitting.

## The important surviving route

Appendix B of the same paper gives a robust construction that:

1. doubles the amplitudes;
2. uses three swap substeps;
3. is designed to remain unitary when links are absent or irregular.

This is not an ad hoc response invented after our failure.  It is a
preregisterable external construction, and it lands exactly in the two
escape classes left open by our earlier theorem:

- enlarged carrier / ancillas;
- incidence depth at least three.

Therefore the next route is sharply defined: test the published robust
three-swap shift on the same (H_4) chamber graph.  Passing that gate would
establish a local unitary scaffold, not a Dirac continuum limit and not the
Whitney metric.

## Status ledger

- **DERIVED:** canonical ordered four-neighbour (H_4) carrier.
- **DERIVED:** balanced chamber-orientation sheets (7,200+7,200).
- **DERIVED NEGATIVE:** the four-amplitude non-robust shift is not unitary on
  those sheets.
- **DERIVED:** repairing it by a free causal-orbit labeling introduces 1,440
  binary choices.
- **STRUCTURAL:** the robust doubled walk is a highly relevant external
  candidate because its design matches our independently derived escape
  boundary.
- **OPEN:** robust-shift transplant.
- **OPEN:** selection of spin coins by 600-cell geometry.
- **OPEN:** refinement convergence to the Dirac equation on (S^3).
- **NOT CLAIMED:** mass, inertia, SI light speed or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_tetrahedral_dirac_walk_bridge.py
```

Expected result: `12/12`.
