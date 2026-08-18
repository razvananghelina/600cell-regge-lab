# Preregistration: transformation-groupoid algebra of the six Hopf labels

Date: 2026-08-10

## Motivation and fixed arena

The commutative label algebra `C(F)=C^6` is closed negative under the current
`A5`-equivariant KO6 metric-zero axioms.  The next algebra is not chosen from
a Wedderburn list: it is the full transformation-groupoid algebra determined
by the already derived effective action

```text
B = C(F) crossed_product A5,
|F|=6.
```

The central `+/-1` in `2I` acts trivially on the unoriented fibration labels,
so this audit uses the effective group `A5`.  A crossed product by the
ineffective binary action would be a different, larger arena and is not
silently identified with this one.

No selector polynomial or Standard-Model target is to be examined before the
algebra type, centre and natural representation are fixed.

## Frozen exact tests

1. Reconstruct the 60-element action on the six fibrations from quaternionic
   conjugation.  Verify transitivity and compute the stabilizer `H` of one
   label, including its element-order census.
2. Use the transitive-groupoid theorem

   ```text
   C(G/H) crossed_product G ~= M_[G:H](C[H])
   ```

   and independently verify dimensions.
3. Derive the complex Wedderburn type of `C[H]` from the stabilizer character
   degrees.  Record the number of simple summands and the full crossed-product
   block sizes before any selector comparison.
4. Construct the natural covariant representation on `H_F=C^6`:

   ```text
   delta_i -> E_ii,
   u_g -> permutation matrix P_g.
   ```

   Compute the exact span rank of all `E_ii P_g`, its commutant and the kernel
   dimension of the 360-dimensional algebra representation.
5. Test whether the six label projections are central.  Compute their
   conjugacy/equivalence orbit under the group unitaries and the centre
   dimension of the full crossed product.
6. Audit the diagonal conditional expectation.  Distinguish:

   - the canonical regular-trace expectation taking the identity-group
     coefficient;
   - arbitrary `C(F)`-bimodular projections, which may retain stabilizer-loop
     data.

   Determine whether the crossed-product algebra itself uniquely forces the
   diagonal response without specifying a trace/state.
7. Only after Steps 1--6, locate `Phi(X)` and `D_aux` in the natural
   representation and state whether their label diagonality is central,
   superselected or merely a chosen conditional expectation.

## Decision boundary

- **Superselection advance:** the six label projections become central or a
  unique intrinsic expectation forces every admissible response into
  `C(F)` while retaining the full algebra's dynamics.
- **Kill for diagonal derivation:** the natural representation is `M6(C)`,
  the labels are mutually equivalent noncentral projections, or diagonal
  projection requires an additional state/trace choice.
- **Structural algebra advance:** even if superselection fails, a completely
  derived noncommutative Wedderburn type may be useful.  It must not be called
  the Standard-Model algebra merely because it is noncommutative.

The crossed product includes all symmetry transitions by construction.  If it
destroys rather than derives label superselection, that is the expected honest
negative, not a reason to delete its off-diagonal blocks.
