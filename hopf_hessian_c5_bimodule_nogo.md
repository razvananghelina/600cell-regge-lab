# Multiplicity cannot rescue the commutative `C^5` Hessian arena

Date: 2026-08-11

Protocol commit: `1bed02d`.

Registered verifier:
`reproducible/verify_hopf_hessian_c5_bimodule_nogo.py`.
Targeted exact result: `11/11`.

## Complete theorem

Let `A=C^5` have a unital faithful representation on a nonzero
finite-dimensional complex Hilbert space.  Assume:

1. order zero;
2. a grading commuting with `A`;
3. KO6 reality `J gamma=-gamma J`;
4. metric-dimension-zero orientability;
5. compatibility of `A`, `J` and `gamma` with the exact derived `A5` action
   on the five primitive algebra characters.

Then no such bimodule exists, for any finite multiplicities.

No Dirac operator, first-order equation, connectedness condition or Hessian
target enters this result.

## Exact five-point action

The verifier reconstructs the 60-element quotient action from the binary
geometry.  An involution has centralizer `V4`; its normalizer is an
index-five `A4`.  Conjugation on the five `A4` subgroups gives the exact
faithful action on the primitive projectors of either monomial `C^5` system.

Its ordered-pair orbit census is

| orbit | size | invariant under `(i,j)->(j,i)` | sees all left/right labels |
|---|---:|---:|---:|
| diagonal | 5 | yes | yes |
| distinct pairs | 20 | yes | yes |

Thus there are only three nonempty invariant support types: diagonal,
off-diagonal, or their union.  Each would already be faithful on both sides.

## Why orientability kills every support

Order zero decomposes every finite carrier as

```text
H=direct_sum_(i,j) H_(i,j).
```

On `H_(i,j)`, every represented Hochschild zero-cycle is a scalar times the
identity of the whole multiplicity space.  Metric-zero orientability
therefore forces

```text
gamma|H_(i,j)=epsilon_(i,j) I,   epsilon_(i,j) in {+1,-1}.
```

Equivariance makes `epsilon` constant on each `A5` orbit, while KO6 reality
requires

```text
epsilon_(j,i)=-epsilon_(i,j).
```

The diagonal orbit gives `epsilon_(i,i)=-epsilon_(i,i)`.  The off-diagonal
orbit contains both orders of every pair, so equivariance makes their signs
equal while KO6 makes them opposite.  The exact sign census gives

```text
nonempty invariant supports tested = 3
orientation-sign solutions        = 0, 0, 0.
```

Increasing a multiplicity cannot split the grading signs inside one block:
zero-cycles remain scalar on that multiplicity space.  This is why the
result covers arbitrary finite multiplicities rather than a bounded search.

## Verdict

**DERIVED FULL-ARENA NO-GO.**  A larger Krajewski carrier does not rescue the
commutative `C^5` route under the stated KO6, metric-zero-orientable and
`A5`-equivariant hypotheses.  The obstruction occurs before a Dirac operator
is chosen, so comparing larger supports with the Hessian would be pointless.

This strengthens the ten-state result `09142e8`:

- the ten-state carrier fails already at zero forms, order zero or first
  order, depending on its sheet algebra;
- every larger `C^5` bimodule fails orientability, independently of its
  multiplicities and independently of the Hessian.

The structural spectral selector remains mathematically valid.  What is
closed is its realization through the canonical commutative five-point
algebra under the current finite-triple axioms.

## Scope limits

The theorem does not cover:

- an independently derived noncommutative algebra;
- explicit breaking of `A5` before orientability;
- an `A5` action incompatible with `J`;
- `J gamma=+gamma J`;
- positive metric-dimension Hochschild orientability.

Each option changes a stated hypothesis.  None is a survivor hidden inside
the arbitrary-multiplicity census.

## Status ledger

- **DERIVED:** the five-point action is reconstructed from the exact group,
  not assumed as a standard permutation representation.
- **DERIVED:** its ordered-pair orbits have sizes `5` and `20`, both stable
  under reversal.
- **DERIVED:** all three nonempty invariant support types have zero KO6
  orientation-sign solutions.
- **DERIVED:** extra multiplicity cannot change a zero-cycle from scalar to
  sign-splitting on a fixed bimodule block.
- **DERIVED FULL-ARENA NO-GO:** no nonzero faithful `A5`-equivariant KO6
  metric-zero-orientable `C^5` bimodule exists.
- **OPEN:** a geometry-selected noncommutative algebra and carrier.
- **STRUCTURAL ADVANCE RETAINED:** the exact affine fourth-moment Hopf
  selection, without a licensed finite-triple realization.

## Next honest boundary

The next route cannot be “add more copies of `C^5`.”  It must first select a
noncommutative algebra from the geometry, blind to the Hessian support, and
only then classify its real bimodules.  The canonical candidate is the
transformation-group algebra of the five-point `A5` action; whether its
natural representation supplies a useful finite triple is still **OPEN**.
