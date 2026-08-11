# The existing action misses the selector, but its minimal affine extension works structurally

Date: 2026-08-11

Protocol commit: `94d8176`.

Verifier: `reproducible/verify_hopf_hessian_action_origin.py`.  Targeted
result: `16/16`.

## Verdict

Two claims must remain separate.

1. **DERIVED NEGATIVE for the existing certified action:** the 2640-state
   Kähler--Dirac certificate contains no label-Hessian block or coupling.  Its
   heat functional is written in `D^2`; inserting a baseline-free Hessian
   block produces only even powers and cannot generate
   `Tr(Hhat_X^3)`.
2. **STRUCTURAL ACTION ADVANCE:** on the same physical five-dimensional label
   module, the complete fourth moment of the unique equivariant affine
   baseline selects a signed Hopf orbit for every nonzero magnitude ratio.
   The desired positive orbit occurs for one of the two relative signs.

The second result is a new operator construction, not something that was
already hidden inside the first.

## What the certified spectral-action file contains

The authoritative file `reproducible/verify_spectral_action.py` constructs the
600-cell cochain carrier of dimension 2640 and certifies

```text
Tr(I), Tr(D^2), (1/2)Tr(D^4),
Tr exp(-tD^2)=c0-c1 t+c2 t^2+O(t^3).
```

It contains no `Hhat_X`, `H_X`, `C_box` or six-label coupling.  This is a
scoped source audit, supplemented by the algebraic parity result:

```text
Tr f(Hhat_(-X)^2)=Tr f(Hhat_X^2).
```

Likewise, a grading-odd baseline-free double has symmetric spectrum and all
odd full traces vanish.  Thus neither direct insertion supplies the required
cubic.

An unspecified coupling between the 2640-state carrier and the new label
module could evade this statement, but that coupling would be new input.

## Exhaustive constant baseline on the physical label module

The exact simultaneous commutant of all 60 `A5` matrices on `W_5=1^perp` has
dimension one:

```text
End_A5(W_5)=R I.
```

Therefore every `A5`-equivariant constant baseline on this fixed carrier is

```text
B_(b,c)(X)=b I_5+c Hhat_X.
```

This is exhaustive for constant affine baselines on `W_5`; it is not a
classification of arbitrary larger carriers or nonlinear operators.

The exact Euclidean realization uses the six-label projector

```text
P=I_6-J_6/6
```

and `B=bP+cH_X`.  Its grading-odd double is

```text
D_(b,c)(X)=[[0,B],[B,0]].
```

It is a self-adjoint 12-dimensional matrix with two constant zero modes, so
its physical restriction has dimension ten.

## Complete fourth moment

Writing `p_k=Tr(Hhat_X^k)`, exact coefficientwise expansion gives

```text
Tr(D_(b,c)(X)^4)
 =10b^4+12b^2c^2 p2+8bc^3 p3+2c^4 p4.
```

Nothing is truncated.  On `q=7200`:

- `p2` is constant;
- `p3=-23328 C_box` has its unique global minima at `+Box_i` and unique
  maxima at `-Box_i`;
- `p4=2e2^2-4e4` has its global minima at all twelve `+/-Box_i`.

Consequently the complete branch classification is

| parameter branch | exact global minima |
|---|---|
| `bc^3>0` | six `+Box_i` |
| `bc^3<0` | six `-Box_i` |
| `b=0, c!=0` | twelve `+/-Box_i` |
| `c=0` | the entire normalized sphere |

For either nonzero sign, the locus is independent of `|b/c|`.  The quartic
term reinforces rather than defeats the cubic selection, so no magnitude or
heat-scale tuning is hidden here.

## Look-elsewhere and sign

There are exactly two nonzero relative-sign branches.  Both select a signed
Hopf orbit; exactly one selects the desired positive orbit:

```text
desired +Box hit fraction = 1/2.
```

Equivariance fixes the baseline line but does not fix the sign of `bc^3`.
Calling the positive branch predicted would therefore overstate the result.
The old convention for `H_X` fixes what the two signs mean; it does not supply
a physical coupling sign.

## Missing spectral-triple gates

The ten-state physical double presently supplies only:

- self-adjointness;
- an even grading with an odd operator;
- exact `A5` covariance;
- a complete fourth-moment selector.

It does not yet supply a represented nontrivial algebra, a real structure
`J`, order zero, first order, orientability, Poincare duality or nonzero inner
one-forms.  In particular, using the scalar algebra would make the algebraic
axioms trivial but would give zero one-forms and no gauge content.

The fixed ten-state carrier has now been audited in
`hopf_hessian_ten_state_triple_verdict.md`.  Its normalized sheet algebras are
exactly `C,C^5,M5(C)`.  The previously omitted monomial `C^5` case exists, but
the variable Hessian family occupies `20` or `25` matrix positions (the span
with `I` occupies all `25`) while first order permits at most a union of two
five-entry matchings.  The other types fail zero forms, first order, or order
zero.  Thus the ten-state real-triple realization is a **DERIVED NEGATIVE**,
not still open.

The arbitrary-multiplicity continuation is now closed for the same
commutative algebra: `hopf_hessian_c5_bimodule_nogo.md` proves that `A5`, KO6
and metric-zero orientability forbid every nonzero `C^5` bimodule before a
Dirac operator is chosen.

The canonical noncommutative continuation has also been audited.  The
five-point crossed product has derived type `M5+M5+M5+M15`, but its natural
ten-state branches fail faithfulness plus order zero or first order; its
canonical faithful doubles fail order zero or orientability.

**OPEN:** derive a physical bimodule independently of the desired selector.
Retrofitting a proper Krajewski support after seeing the obstruction would
reintroduce the fitting problem.

The exhaustive blind central-support census now sharpens this warning:
`256/729` abstract supports pass necessary gates and 24 attain the minimum
dimension 300, but none is compatible with exact `chi<->chibar` conjugation.
The crossed product supplies existence without selection.

## Status ledger

- **DERIVED NEGATIVE:** the existing certified Kähler--Dirac action contains
  no label-Hessian coupling.
- **DERIVED NEGATIVE:** baseline-free `D^2` parity cannot generate the odd
  selector cubic.
- **DERIVED:** `dim End_A5(W_5)=1`; the constant affine baseline is `bI`.
- **DERIVED:** the full fourth moment has the displayed four-term expansion.
- **DERIVED:** every `bc^3>0` magnitude ratio selects exactly the six
  `+Box_i`; every negative relative sign selects `-Box_i`.
- **PATTERN/CHOICE:** the desired relative-sign hit is `1/2` and is not fixed
  by current physics.
- **STRUCTURAL ADVANCE:** a minimal grading-odd operator turns the exact
  Hessian cubic into a robust complete even spectral moment.
- **DERIVED TEN-STATE NO-GO:** no normalized algebra on the fixed physical
  double realizes the complete affine family with nonzero one-forms.
- **DERIVED FULL-ARENA C5 NO-GO:** arbitrary bimodule multiplicity cannot
  restore KO6 metric-zero orientability.
- **STRUCTURAL ALGEBRA ADVANCE:** the canonical crossed product is
  `M5+M5+M5+M15`.
- **DERIVED CANONICAL-CARRIER NEGATIVE:** its natural and functorial carriers
  fail before Hessian comparison.
- **PATTERN/SELECTION NEGATIVE:** hundreds of abstract proper supports survive
  necessary gates, while exact character conjugation forbids every
  nondegenerate one.
- **OPEN PHYSICS:** no noncommutative all-gate finite spectral triple or
  coupling to the certified 2640-state geometry has been derived.

## Next boundary

The ten-state question, the arbitrary-multiplicity commutative `C^5`
continuation and the canonical crossed-product carriers are all closed
negatively.  The geometry now selects a noncommutative algebra but not a
physical bimodule over it.  Until it does, the selector remains an effective
order-parameter construction rather than a sector of the theory's spectral
geometry.
