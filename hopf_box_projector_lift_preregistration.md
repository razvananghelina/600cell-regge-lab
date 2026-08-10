# Preregistration: lift the 5D projector simplex into the six `Box_F` operators

Date: 2026-08-10

This protocol is committed before computing the operator-space Gram matrix or
any moment coefficient.  The proposed affine construction was recognized
before registration, so this is not blind target discovery.  The numerical
outcomes and, in particular, whether the cubic coefficient vanishes and what
sign it has, have not been inspected.

## Fixed input

Use only the six discrete Hopf fibrations already enumerated by
`verify_hopf_fibration_invariants.py`.  For each fibration `F_i`, use its
integer fibre-edge adjacency `A_f,i`, the fixed 600-cell adjacency `A`, and
the already derived coefficient `b1=6`:

```text
Box_i = 6 A_f,i - A.
```

No weights or alternative operators may be searched.

Define the affine centre and centered operators

```text
Box_bar = (1/6) sum_i Box_i = sum_i A_f,i - A,
B_i     = Box_i - Box_bar = 6 A_f,i - sum_j A_f,j.
```

Use the trace/Hilbert--Schmidt pairing.  Compare these `B_i` with the six
centered projectors `T_i=P_i-I/3` only through the unique linear map fixed by
`T_i -> B_i`; an overall rescaling is to be recorded, not adjusted.

## Tests and falsifiers

1. Compute the exact integer Gram matrix `Tr(B_i B_j)`, its rank, and the
   sum of the `B_i`.  If it is not a regular 5-simplex, the canonical linear
   lift of the projector order parameter fails.
2. Parameterize the zero-sum barycentric space without fitting and compute
   the homogeneous cubic part of

   ```text
   Tr((Box_bar + X)^4),  X in span{B_i}.
   ```

   Compare the complete polynomial, coefficient by coefficient, with the
   projector cubic `C3(Q)=sum_i Tr(QT_i)^3` under `T_i -> B_i`.
3. Record whether the fourth-moment cubic is zero.  If nonzero and
   proportional to `C3`, record the exact proportionality and its sign.
4. Also compute the cubic part of `Tr((Box_bar+X)^3)` and the full values of
   the third and fourth moments at all six vertices as consistency checks.

The route is killed at this stage if the operator simplex has rank other than
five or if the fourth-moment cubic vanishes.  A nonzero coefficient advances
only an operator-level structural route: `Box` is not the certified
Kahler--Dirac spectral-action operator, so no physical action claim follows
without an additional licensed bridge.

## Interpretation boundary

Even a successful result does not prove that arbitrary affine combinations
of the six `Box_i` are admissible Dirac operators, that the fourth moment is a
physical potential, or that its overall spectral coefficient has the needed
sign.  Those are separate gates.
