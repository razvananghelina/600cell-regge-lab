# Preregistration: can the ten-state Hessian selector be a real finite triple?

Date: 2026-08-11

## Complete hypotheses

Fix the physical carrier and symmetry

```text
H = W_+ direct_sum W_-,
dim_C W_+=dim_C W_-=5,
gamma=diag(+I_5,-I_5),
U_g=diag(rho_5(g),rho_5(g)),
```

where `rho_5` is the already derived irreducible real five-dimensional
`A5` representation.  Fix the nonzero grading-odd selector family

```text
D_X=[[0,B_X],[B_X^*,0]],
B_X=bI+cHhat_X,
```

on a nonconstant branch `c!=0`.

Ask whether there exists a finite-dimensional unital complex `*`-algebra `A`
with faithful representation `pi` such that:

1. `[gamma,pi(A)]=0`;
2. the derived `A5` action normalizes `pi(A)`;
3. a KO6 real structure satisfies `Jgamma=-gamma J` and order zero;
4. `D_X` satisfies first order for every admitted `X`;
5. represented inner one-forms are nonzero;
6. connectedness holds.

Orientability and Poincare duality are later gates only if a candidate passes
the cheaper necessary conditions above.  Failure before them is not promoted
to a theorem under hypotheses that were not tested.

This protocol concerns the fixed ten-state carrier.  Adding spectator
bimodules or higher multiplicities is a different arena.

## STEP 1: classify algebra images, not partitions

For each orientation sheet, classify every unital `*`-subalgebra of `M5(C)`
normalized by the irreducible `A5` action.

Use the finite-dimensional structure theorem:

```text
image ~= direct_sum_alpha (M_(k_alpha)(C) tensor I_(m_alpha)),
sum_alpha k_alpha m_alpha=5.
```

The group permutes equal central supports.  Check all possible centre-orbit
sizes and use the minimal faithful permutation degree of `A5`.  In
particular, a transitive five-support action would make the sheet a
five-point permutation module containing a fixed vector, which `W_5` does
not.

Decision target:

```text
the only normalized sheet images are C I_5 and M5(C).
```

If this claim fails, enlarge the subsequent exact type list before testing
the selector.  Do not silently assume it.

## STEP 2: enumerate joint faithful types

Given the two sheet images, enumerate the joint faithful representation types,
including distinct kernels on the two sheets.  The expected list to test, not
assume, is:

```text
C,
C+C,
M5,
M5+C,
M5+M5,
```

with sheet placements and graph-of-automorphism cases identified.

## STEP 3: exact gate order

For each surviving type test, in this order:

1. order zero with an arbitrary KO6 sheet-exchanging antiunitary `J`;
2. first order for the full nonzero off-diagonal block `B_X`;
3. nonzero represented one-forms;
4. connectedness.

Use central idempotents wherever possible.  If `e` distinguishes the two
sheet summands, compute exactly

```text
[[D_X,pi(e)], J pi(e) J^-1].
```

Do not infer a general first-order no-go from a search that included only
diagonal matrix representatives.

## Decision boundary

- **Finite-triple advance:** a noncommutative algebra type passes order zero,
  first order, nonzero forms and connectedness with the complete selector
  block.  Then continue to `J^2`, `JD`, orientability and Poincare duality.
- **Abelian structural witness:** only a commutative type survives every cheap
  gate.  State plainly that it supplies no non-abelian gauge sector.
- **Ten-state kill:** scalar `C` has zero forms, all split-support types force
  `B_X=0` by first order, and full-matrix types fail order zero.  This closes
  the fixed carrier but not larger bimodule completions.

## Framing attack

`dim End_A5(W_5)=1` classifies equivariant *operators*, not automatically all
normalized subalgebras; the central-support argument is required separately.
Likewise, failure on ten states cannot rule out embedding the same selector
block into a larger Krajewski diagram.  Any final no-go must state both scope
limits explicitly.
