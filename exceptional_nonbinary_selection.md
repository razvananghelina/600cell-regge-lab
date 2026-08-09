# Exceptional non-binary selection: `E6/E7/E8` control

## Scope

This note compares only the three exceptional binary polyhedral groups

`2T <-> affine E6`, `2O <-> affine E7`, `2I <-> affine E8`.

It does not classify the binary-dihedral family, all finite subgroups of
`SU(2)`, or all fusion categories.  Every uniqueness statement below is
restricted to this exceptional ladder.

## Exact phase ladder

The traces of the fundamental maximal cyclic phases are

`2 cos(pi/3)=1`, `2 cos(pi/4)=sqrt(2)`, `2 cos(pi/5)=phi`.

Their phase polynomials give `Phi_6`, and by quadratic field norm `Phi_8`
and `Phi_10`, respectively.  **DERIVED.**

These dimensions have the familiar minimal fusion readings:

- `d=1`: pointed/invertible, represented by `x^2=1`;
- `d=sqrt(2)`: Ising branching, `sigma^2=1+psi`, with a new simple `psi`;
- `d=phi`: Fibonacci non-branching, `tau^2=1+tau`.

Consequently only the `2I/E8` entry is compatible with all four structural
contents of S01 under this dictionary: noninvertibility, unit return, self
return and no new simple type.  The dimension equations are **DERIVED**;
identifying a group-element phase trace with the quantum dimension of the
foundational fusion object is **STRUCTURAL**, not derived from group
representation theory.  The resulting selection is therefore
**DERIVED conditional on that dictionary and on the displayed three cases.**

## Galois-memory selection

The corresponding normalized conjugate ratios are

`1`, `-1`, `sigma(phi)/phi=-phi^-2`.

The first has no nontrivial arithmetic conjugate and no decay.  The second
is a pure period-two oscillation.  Only the third is nonzero and strictly
contractive.  It is exactly the nonstationary eigenvalue of the canonical
Fibonacci Perron process.  **DERIVED conditional selection within the
exceptional ladder.**

Thus two independent screens agree:

1. S01 non-branching productive self-reference;
2. nonzero strictly contractive Galois memory.

Both select `2I/E8` over `2T/E6` and `2O/E7`.  This is a genuine comparative
control, not a global uniqueness theorem.

## Binary-dihedral counterexample

The restriction to the exceptional ladder is load-bearing.  The
binary-dihedral group `Dic_5`, of order 20 and affine-`D7` McKay type, also
has a maximal cyclic `C10`.  Hence it has exactly the same fundamental phase
trace `phi`, the same cyclotomic polynomial `Phi_10`, and the same normalized
Galois ratio `-phi^-2`.

Therefore phase trace, S01 compatibility and strict Galois contraction do
not select `2I/E8` among all binary polyhedral groups.  `Dic_5/D7` is an
explicit **DERIVED counterexample** to that stronger claim.  An additional
selector must distinguish the full group/McKay topology, not merely its
maximal cyclic fiber.

## Honest boundary

- **DERIVED:** every displayed trigonometric, cyclotomic, fusion-dimension and
  Galois-ratio identity.
- **DERIVED CONDITIONAL:** uniqueness inside the three exceptional cases.
- **STRUCTURAL:** treating maximal cyclic phase trace as the input read by
  the foundational fusion seed.
- **DERIVED NEGATIVE:** `Dic_5/D7` defeats global phase-based uniqueness.
- **OPEN:** a non-phase principle explaining why the exceptional ladder, or
  specifically the affine-`E8` topology, is selected.
- **NOT CLAIMED:** selection of the Standard Model, a spectral triple, or a
  physical decay rate.

Exact verifier:
`reproducible/verify_exceptional_nonbinary_selection.py`.
