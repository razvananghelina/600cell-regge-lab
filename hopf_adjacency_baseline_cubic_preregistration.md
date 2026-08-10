# Preregistration: canonical adjacency-baseline cubic audit

Date: 2026-08-10

## Provenance and scope

The desired equal-weight projector cubic is already known.  This is therefore
not a blind discovery protocol.  Its purpose is to prevent changing the
operator family or the admissible trace words after seeing coefficients.

The previous audit closed the zero-baseline family `D(X)=X`: its third moment
does not select the six Hopf operators and its fourth moment is even.  The
next cheapest operator already present on the same 120-vertex carrier is the
600-cell adjacency `A`.  At a derived Hopf vertex,

```text
A + Box_i = 6 A_f,i.
```

This audit tests that fixed nonzero baseline.  It does not assume that
`A+X` is already a licensed fluctuated Dirac operator; that physical status
must be judged separately even if the polynomial test succeeds.

## Fixed data

Use only

```text
W = span_R{Box_i},                 dim_R W = 5,
Box_i = 6 A_f,i-A,                sum_i Box_i = 0,
D_A(X) = A+X.
```

No matrix insertion besides `A`, no weights on the six fibrations, and no
free linear combination of trace words may be introduced.

The only homogeneous cubic supplied by the fourth single-trace moment of
this fixed affine family is

```text
G_A(X) = 4 Tr(A X^3),
```

because cyclicity makes the four words with one `A` and three `X` equal.
For comparison, the zero-baseline third-moment cubic is fixed as

```text
G_0(X) = Tr(X^3).
```

## Frozen exact tests

1. Reconstruct the six fibrations, `A`, `Box_i` and the five-dimensional
   basis `E_a=Box_a-Box_5` from the registered integer geometry.
2. Independently compute the character of the induced `A5` action on `W` and
   the multiplicity of the trivial representation in `Sym^3(W*)`.  This
   gives the complete dimension of symmetry-allowed cubic freedom before
   comparing any trace word with the desired selector.
3. Compute `G_0` and `G_A` coefficientwise in exact integer arithmetic and
   determine their span rank.
4. Only then construct the equal-weight analysis cubic intrinsic to the
   already derived operator simplex,

   ```text
   C_box(X) = sum_i Tr(X Box_i)^3.
   ```

   Compare it coefficientwise with `G_A` and with the span of `(G_0,G_A)`.
5. If `G_A` is proportional to `C_box`, record the exact proportionality and
   determine whether a positive coefficient of `Tr((A+X)^4)` favours
   `+Box_i` or `-Box_i` on the fixed quadratic sphere.  Do not infer the
   global minima of the full fourth moment from its cubic part alone.
6. Independently test the complete fixed functional

   ```text
   S_4(X)=Tr((A+X)^4)
   ```

   on `Tr(X^2)=7200`: stationarity at all `+/-Box_i`, exact values there, and
   an exhaustive or exact counterexample test for global selection.  A
   nonstationary desired vertex is already a kill.

## Decision boundary

- **Algebraic advance:** `G_A` is a nonzero multiple of `C_box`, so the
  canonical adjacency baseline supplies the correct cubic tensor without
  fitted Schur coefficients.
- **Dynamical advance:** in addition, the complete fixed fourth moment has
  exactly the six `+Box_i` as its selected extrema with the required sign.
- **Algebraic kill:** `C_box` is not proportional to `G_A`.
- **Dynamical kill:** the correct cubic occurs but its sign is wrong, the
  desired vertices are nonstationary for the full `S_4`, or exact additional
  minimizers/stronger extrema exist.

Even a dynamical advance is only **STRUCTURAL** until `A+X` is realized as a
licensed fluctuated self-adjoint Dirac operator satisfying the spectral-triple
gates.  A match in this 120-vertex trace family cannot silently promote it to
the already certified 2640-dimensional Kahler--Dirac spectral action.
