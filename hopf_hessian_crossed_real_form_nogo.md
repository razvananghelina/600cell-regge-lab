# The canonical real crossed product fails KO6 Poincare duality

Date: 2026-08-11

Protocol commit: `5eafb29`.

Registered verifier:
`reproducible/verify_hopf_hessian_crossed_real_form.py`.
Targeted exact result: `12/12`.

No Hessian or selector target is used.

## Complete statement

Let

```text
B_R=R(P) crossed_product A5
```

be the canonical real transformation-group algebra of the derived five-point
action.  Under KO6, metric-dimension-zero orientability and Poincare duality
over `K0(B_R)`, no finite real spectral triple over `B_R` exists.

The result covers arbitrary bimodule multiplicities and arbitrary Dirac
operators because it is an algebra-type parity obstruction.

## Exact real Wedderburn derivation

The point stabilizer is the exact `A4` subgroup.  Its conjugacy classes have
sizes and element orders

```text
sizes  =1,3,4,4,
orders =1,2,3,3.
```

The four complex irreducible characters have degrees `1,1,1,3`.  Constructed
exactly in `Q(omega)`, they are orthonormal and have Frobenius--Schur
indicators

```text
nu(1)=1,
nu(chi)=0,
nu(chibar)=0,
nu(3)=1.
```

Thus the trivial and three-dimensional representations are of real type,
while `chi,chibar` form one complex-type real irreducible pair.  Consequently

```text
R[A4] ~= R + C + M3(R)
```

and

```text
B_R ~= M5(R) + M5(C) + M15(R).
```

The real dimensions are

```text
25+50+225=300,
```

matching `|P| |A5|=5*60` exactly.

Complexification gives

```text
B_R tensor_R C
 ~= M5(C)+M5(C)+M5(C)+M15(C).
```

The middle real `M5(C)` block splits into the two complex blocks induced from
`chi` and `chibar`.  Therefore their conjugation is not an optional symmetry
added after the fact: it is how the canonical real algebra complexifies.

## KO6 parity obstruction

For minimal `K0` projections `p_i`, the intersection form is

```text
cap_ij=Tr(gamma pi(p_i) J pi(p_j) J^-1).
```

Conjugating the trace by the antiunitary `J`, using order zero and

```text
J gamma J^-1=-gamma,
```

gives

```text
cap_ji=-cap_ij.
```

The canonical real algebra has three simple summands, hence

```text
rank K0(B_R)=3.
```

Its generic KO6 intersection form is therefore

```text
[ 0  x  y]
[-x  0  z]
[-y -z  0],
```

whose determinant is identically zero and whose generic rank is two.
It can never implement a nondegenerate pairing on a rank-three `K0` group.

## Verdict

**DERIVED FULL-ARENA POINCARE NO-GO, under the canonical-real-form
hypothesis.**  No multiplicity matrix, real structure realization or Dirac
operator can repair the odd-rank antisymmetric form.

This is the real-algebra explanation of the previous four-complex-node
Pfaffian result.  In the complexification, reality exchanges `chi` and
`chibar`; imposing either grading behavior forced the four-node Pfaffian to
vanish.  Over the real form, the same obstruction is simply the impossibility
of a nondegenerate antisymmetric form on three generators.

## Scope and framing

This theorem does not say that every complex spectral triple over

```text
M5(C)^3+M15(C)
```

is impossible.  Treating the conjugate `chi,chibar` blocks as independent
complex algebra summands yields four `K0` nodes and evades the odd-rank
argument.  But that changes the canonical real transformation-group algebra
into its complexification and discards its real descent condition.  The
blind census found 256 abstract necessary-gate supports in precisely that
larger complex arena, with no geometric selection among them.

Other genuine exits are:

- changing KO dimension so the transpose sign changes;
- dropping Poincare duality;
- choosing a different real algebra independently from the geometry.

Each is a changed hypothesis, not a survivor of the proved real-form arena.

## Status ledger

- **DERIVED:** exact FS indicators are `1,0,0,1`.
- **DERIVED:** `R[A4]=R+C+M3(R)`.
- **DERIVED:** `B_R=M5(R)+M5(C)+M15(R)`, real dimension 300.
- **DERIVED:** its complexification is `M5^3+M15`, with conjugation pairing
  the two nontrivial-character blocks.
- **DERIVED:** `rank K0(B_R)=3` and the KO6 intersection form is
  antisymmetric.
- **DERIVED FULL-ARENA NO-GO:** strict Poincare duality is impossible for the
  canonical real crossed product, arbitrary multiplicities and `D` included.
- **STRUCTURAL EXISTENCE / SELECTION NEGATIVE:** the independent four-node
  complexification has many abstract supports but no selected one.
- **OPEN:** a different geometry-derived real algebra with even `K0` rank.

## Programme consequence

The real-form audit removes the last ambiguity in the conjugation argument.
The canonical five-point crossed product does not merely fail on its natural
carriers: its real algebra type is globally incompatible with strict KO6
Poincare duality.  Continuing inside its complexified four-node support space
would require first explaining why the real geometry has lost its canonical
real descent.
