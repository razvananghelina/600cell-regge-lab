# The affine Hessian selector has no real-triple realization on the fixed ten-state carrier

Date: 2026-08-11

Initial protocol commit: `a0d7dd8`.

Protocol correction commit: `cae7233`.

Registered verifier: `reproducible/verify_hopf_hessian_ten_state_triple.py`.
Targeted exact result: `21/21`.

## Complete scope

Fix

```text
H=W_+ direct_sum W_-,       dim_C W_+=dim_C W_-=5,
gamma=diag(+I_5,-I_5),      U_g=diag(rho_5(g),rho_5(g)),
D_X=[[0,B_X],[B_X^*,0]],    B_X=bI_5+c Hhat_X, c!=0.
```

Here `rho_5` is the derived irreducible real five-dimensional `A5`
representation and `X` ranges over the complete admitted five-dimensional
Hopf--Box field module.  Ask for a unital complex finite-dimensional
`*`-algebra with faithful representation such that:

1. it commutes with `gamma`;
2. its represented image is normalized by the fixed `A5` action;
3. a KO6 sheet-exchanging real structure satisfies order zero;
4. the complete affine family `D_X` satisfies first order;
5. represented inner one-forms are nonzero;
6. connectedness holds.

Orientability and Poincare duality are later gates.  Since every algebra type
fails one of the cheaper necessary gates, this audit makes no claim that
those later axioms were tested.

## Correction made before the gate result

The initial expected sheet-image list `C,M5(C)` was false.  A diagonal
algebra is normalized by monomial, not merely permutation, matrices.  The
exact group calculation finds

```text
A4/V4 = C3
```

and constructs the two nontrivial character projectors in
`Q(omega)`, `omega^2+omega+1=0`.  For each character, its `A5` orbit consists
of exactly five projectors which are:

- rank-one idempotents;
- mutually orthogonal and sum to `I_5`;
- self-adjoint in the inherited physical Gram metric on `W_5`;
- permuted by all 60 exact `A5` matrices.

Thus there are two distinct normalized monomial `C^5` systems.  This was
preregistered as a correction in `cae7233` before their first-order support
was evaluated.

## Exact sheet-algebra classification

The exact commutant of the 60 matrices on `W_5` has complex dimension one,
so `W_5` is irreducible.  Let `r` be the number of primitive central supports
of a normalized sheet algebra.

- For `2<=r<=4`, a nontrivial centre action would inject the simple group
  `A5` into `S_r`, which is impossible.  A trivial action would make each
  support an invariant subspace, also impossible by irreducibility.
- For `r=5`, transitivity and total dimension five make every support
  one-dimensional.  The trivial `A4` stabilizer character would give the
  reducible permutation module `1+4`; the two nontrivial characters give
  exactly the two systems constructed above.
- For `r=1`, the image is `M_k(C) tensor I_m` with `km=5`.  Since five is
  prime, it is either `C I_5` or `M5(C)`.

Therefore the exhaustive normalized sheet-image list is

```text
C, C^5, M5(C).
```

This classification uses the `A5`-normalization hypothesis.  Dropping that
hypothesis would admit continuously many fitted embeddings and is outside
the stated result.

## The decisive `C^5` first-order calculation

For a faithful five-point algebra on each sheet, order zero makes the left
and opposite characters into two label maps.  The first-order double
commutator implies that an entry of the off-diagonal Dirac block may be
nonzero only when its endpoints share their left label or share their
opposite label.  After arbitrary relative relabelling, these two conditions
are two permutation matchings.  Hence every allowed support contains at most

```text
5+5=10 of the 25 matrix positions.
```

The verifier separately evaluates the five variable directions and the
constant affine baseline:

```text
I_5, Hhat_(X_1), ..., Hhat_(X_5)
```

between every ordered pair of the two character systems.  The exact supports
are

| plus system | minus system | variable Hessian | span including `I` | two-rook cover |
|---|---:|---:|---:|---:|
| `chi` | `chi` | 20/25 | 25/25 | no |
| `chi` | `chibar` | 25/25 | 25/25 | no |
| `chibar` | `chi` | 25/25 | 25/25 | no |
| `chibar` | `chibar` | 20/25 | 25/25 | no |

The calculation is exact in `Q(omega)`: there is no numerical tolerance.
Even on the baseline-free branch `b=0`, the support has at least 20 entries,
twice the first-order ceiling.  For `b!=0`, the affine span occupies all 25.
Phases and relative permutations cannot turn a nonzero position into zero.
Thus no normalized `C^5` embedding and no admissible KO6 relabelling can
contain the complete variable Hessian family for any value of `b`.

## Exhaustion of joint algebra types

Each sheet activates either one scalar node, five scalar nodes, or one
`M5(C)` node.  Faithfulness says that the union of nodes seen on the two
sheets is the complete algebra.

- If two five-point node sets overlap, their overlap is invariant under the
  transitive `A5` action.  It is therefore empty or all five.  These are the
  split case and the `C^5` graph case respectively.
- A scalar node cannot be partially shared with a transitive five-node
  system: it would define an `A5`-fixed point.  Different simple block sizes
  cannot be shared.
- Two `M5(C)` nodes are either the same graph block or two split blocks.

Consequently no intermediate subdirect product was omitted.  The necessary
gates are:

| joint category | exact obstruction |
|---|---|
| graph scalar `C` | every represented one-form is zero |
| graph `C^5` | `20/25` or `25/25` variable support violates the ten-position first-order ceiling |
| any split node sets | a sheet central idempotent gives `[[D,e],e^o]` with the complete nonzero `B_X` block |
| graph `M5(C)` | order zero asks a full `M5` to commute with a full `M5` on one sheet |

## Verdict and framing limits

**DERIVED TEN-STATE NO-GO.**  Under the complete hypotheses above, the fixed
carrier `W_+ direct_sum W_-` cannot realize the affine Hessian selector as a
real finite spectral triple with nonzero one-forms.  The difficult omitted
`C^5` case exists, but fails first order decisively: at least `20` occupied
positions versus an allowed maximum of `10`.

This is not a no-go for the spectral selector itself.  The exact fourth
moment and its Hopf-orbit selection remain valid as a **STRUCTURAL ADVANCE**.
What fails is their realization on the minimal ten-state carrier.

This is also not a theorem about larger Krajewski diagrams.  Extra bimodule
multiplicities can change the first-order support.  Such a completion must be
selected independently; adding spectator states solely to route around the
`at least 20 versus 10` obstruction would reintroduce fitting.

Freezing `X` at one already-selected configuration is likewise outside the
hypotheses.  The result concerns the complete field family, as preregistered,
not an algebra retrofitted to one vacuum.

## Status ledger

- **DERIVED CORRECTION:** `W_5` has two normalized monomial `C^5` algebras.
- **DERIVED:** the exhaustive normalized sheet images are `C,C^5,M5(C)`.
- **DERIVED:** the four character-system pairs have exact variable-Hessian
  supports `20,25,25,20`; adding the constant affine direction gives `25/25`
  in every pair.
- **DERIVED NEGATIVE:** no `C^5` first-order rook support contains the affine
  family.
- **DERIVED NEGATIVE:** split centres force the nonzero selector block to
  vanish; a full-matrix graph fails order zero; the scalar graph has zero
  one-forms.
- **DERIVED TEN-STATE NO-GO:** every joint type fails a necessary gate.
- **STRUCTURAL ADVANCE RETAINED:** the affine fourth moment still selects a
  signed Hopf orbit.
- **OPEN:** a geometry-selected larger bimodule/Krajewski completion.

## Next honest boundary

Do not search arbitrary larger matrices.  First derive the smallest
`A5`-equivariant bimodule multiplicity pattern for which a first-order rook
support can contain the five-dimensional Hessian family.  Only after that
blind carrier enumeration should the physical selector be compared with the
survivors.
