# The robust three-swap scaffold works, with a published formula discrepancy

Date: 2026-08-11  
Preregistration commit: `c80f448`

## Result

The three literal robust shift stages in Appendix B of Nzongani *et al.*
define an exact local permutation on the first-barycentric (H_4) chamber
carrier.

> **DERIVED STRUCTURAL BRIDGE:** on 14,400 chambers with eight components
> each, all three stages (S_0,S_1,S_2) and their synchronous global product
> are permutations of 115,200 states.  The construction needs no chamber
> handedness label and crosses at most one dual edge per macro tick.

> **DERIVED FORMULA DISCREPANCY:** Eq. (40) in the final Physical Review A
> article is not the global product of Eqs. (39a)--(39c).  The printed column
> is non-bijective; the literal operator product is bijective.

The targeted verifier passes `13/13` in about 0.7 seconds.  No full suite was
run.

## Literal global product

Let (n_2,n_3) be the two intrinsic rank-colour involutions used by the
paper's robust shift.  Composing the three published stages synchronously on
every chamber gives

\[
(S_2S_1S_0\phi)(k)=
\begin{pmatrix}
\phi_2(k)\\
\phi_3(n_2(k))\\
\phi_0(n_3(k))\\
\phi_1(k)\\
\phi_4(k)\\
\phi_5(k)\\
\phi_6(k)\\
\phi_7(k)
\end{pmatrix}.
\]

This is a permutation because (n_2) and (n_3) are involutions.  Output
components 1 and 2 cross one chamber facet; all others remain in the same
chamber.

Components 4--7 return exactly to their original state after a macro tick,
but the individual stages mix the active and ancillary sectors.  Thus the
ancillas mediate the stricter three-swap implementation even though the
macro operator preserves the active/ancilla split.

## Exact discrepancy in the published article

The final article is [Phys. Rev. A **110**, 042418
(2024)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.110.042418),
not merely an unreviewed preprint.  Its Eq. (40) prints components 1 and 2 as

\[
\phi_6(n_2(k)),\qquad \phi_5(n_3(k)),
\]

rather than the literal-composition values

\[
\phi_3(n_2(k)),\qquad \phi_0(n_3(k)).
\]

The difference affects exactly

\[
2\times14{,}400=28{,}800
\]

output positions.  Treated as a global map, the printed Eq. (40) uses only
86,400 distinct inputs out of 115,200: 28,800 inputs occur twice and 28,800
do not occur.  It therefore cannot be unitary.

The likely algebraic issue is that the displayed expansion applies (S_0)
to the observed chamber but not to the neighbouring amplitudes subsequently
read by (S_1).  A synchronous global local operation acts on those
neighbours as well.  This is an inference from the equations, not a claim
about the authors' intent.

The stage definitions are mutually consistent and yield the permutation
above.  They are therefore retained as the coherent robust construction;
the printed expanded column is not used as evidence.

## What this does and does not rescue

This result is genuinely useful:

- it supplies an external, already published three-substep unitary scaffold;
- it works on the canonical four-coloured (H_4) chamber graph;
- it avoids the 1,440 free handedness labels that killed the economical
  transplant;
- its use of ancillas and multiple stages is qualitatively consistent with
  the escape routes left open by our separate weighted-incidence no-go.

But it is not yet our physical evolution:

- its carrier is chamber-facet spinor data, not the 2,640 Kähler--Dirac
  cochains;
- it has not been related to the exact Whitney metric;
- the Pauli-direction coin was chosen in the paper to recover a target Dirac
  equation, not selected from the 600-cell axioms;
- the published continuum proof uses a flat periodic orthoscheme lattice;
- only colours 2 and 3 occur in this translation stage, so global propagation
  on (H_4) still needs a separate connectedness audit.

## Status ledger

- **DERIVED:** all three literal robust stages are exact permutations.
- **DERIVED:** the 115,200-state macro shift is unitary and strictly local on
  the chamber dual graph.
- **DERIVED:** ancillas are used at intermediate stages and return after the
  macro step.
- **DERIVED DISCREPANCY:** final published Eq. (40) is neither the literal
  product nor a permutation.
- **STRUCTURAL:** a robust local unitary scaffold now exists on the geometry.
- **OPEN:** whether local coins plus this shift connect the whole chamber
  carrier rather than ten-chamber causal orbits.
- **OPEN:** a geometry-selected three-axis schedule.
- **OPEN:** Whitney/Kähler--Dirac relation and refinement convergence.
- **NOT CLAIMED:** a new particle law, mass derivation, (c), (G),
  (hbar), Planck time or Planck mass.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_tetrahedral_dirac_walk_robust.py
```

Expected result: `13/13`.
