# Preregistration: can existing dynamics induce Hopf-label superselection?

Date: 2026-08-10

## Question and complete hypotheses

Let `F` be the six unoriented Hopf fibrations already derived from the
600-cell, let `H_F=C^6`, and let

```text
C(F) = span{E_00,...,E_55} subset M6(C)
```

be the label-diagonal algebra.  Use no new label algebra, weights, bath,
Hamiltonian, regulator or fitted linear combination.

The only two exact averaging/dephasing mechanisms already fixed by the
construction are:

1. the effective `A5` permutation action on `H_F`;
2. the complete label Hessian of the already audited extended cubic,

   ```text
   H_X(i,j)=3 Tr(X(Box_i Box_j+Box_j Box_i)),
   X in W=span_R{Box_i}, Tr(X^2)=7200.
   ```

This audit asks whether either mechanism makes the six rank-one projections
exact superselection sectors.  Approximate decoherence, a chosen environment
and the diagonal auxiliary operator `D_aux` are out of scope: the first two
are new inputs, while using `D_aux` to derive the diagonality assumed in its
own construction would be circular.

## Frozen tests

### A. Symmetry twirl

Construct the exact Reynolds map

```text
E_A5(M)=(1/60) sum_g P_g M P_g^-1.
```

Compute its superoperator rank and image.  Check whether its fixed algebra is
`C(F)`, contains the six individual `E_ii`, or instead is the commutant of the
permutation representation `1+5`.  Also apply the simultaneous group average
to the linear order parameter/Hessian family and record whether it preserves
any nonzero `X` direction.

### B. Hessian-generated dephasing

For a self-adjoint generator `H_X`, the infinite-time Heisenberg average is
the conditional expectation onto `{H_X}'`.  Exact label superselection
requires every `E_ii` to commute with `H_X`; equivalently `H_X` must be
diagonal in the fibration basis.

Compute the exact linear map from `X` to the off-diagonal entries of `H_X`.
Its kernel on `W` decides whether any nonzero normalized `X` can preserve all
six labels.  At every six desired `Box_i`, additionally compute:

- the commutators `[H_Box_i,E_jj]`;
- the connectivity of the nonzero off-diagonal graph;
- `C(F) intersect {H_Box_i}'`;
- the full fixed-algebra dimension from the exact eigenvalue
  multiplicities.

The last two quantities distinguish “some degeneracy remains” from “the six
label sectors remain.”

### C. Honest interpretation

Do not call a group or time average a derivation of locality merely because it
is a canonical conditional expectation.  The image algebra itself must be
the six-label diagonal algebra, or must contain its six minimal projections
as fixed observables.

## Decision boundary

- **Advance:** one of the two already-defined parameter-free mechanisms fixes
  all six `E_ii` and suppresses every inter-label matrix unit, without using
  `D_aux` or adding a state, bath or regulator.
- **Kill for existing-dynamics superselection:** the symmetry fixed algebra
  is smaller/different from `C(F)`, and the full Hessian has nonzero
  off-diagonal entries for every normalized `X` (in particular at all six
  `Box_i`).
- **Structural only:** a canonical averaging map exists but selects a
  different algebra, or a noncanonical bath/state could be chosen to obtain
  `C(F)`.

## Framing limitation

This protocol can close superselection from the dynamics already present in
the repository.  It cannot prove that no future, independently derived local
environment or dynamical field can generate superselection.  Conversely,
inventing such an environment after this test would be a new physical
hypothesis, not a repair inside the audited construction.
