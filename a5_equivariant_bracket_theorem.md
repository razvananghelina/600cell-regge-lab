# A5-equivariant Lie brackets on `3' + 5`

Date: 2026-07-22

## Equivariant-map space

Let `W=3'+5` over the reals, with the repository convention
`chi_3'(5A)=phi'` and `chi_3'(5B)=phi`. Character calculation gives

`Lambda^2 W = 2(3) + 3(3') + 2(4) + 5`.

Therefore:

- **DERIVED:** `dim Hom_A5(Lambda^2 W,3')=3`;
- **DERIVED:** `dim Hom_A5(Lambda^2 W,5)=1`;
- **DERIVED:** the full space of equivariant antisymmetric maps
  `Lambda^2 W -> W` has dimension `4`.

Use the explicit model `W=so(3)+Sym_0(3)`, with elements `(A,S)` and
`a=vee(A)`. Every equivariant antisymmetric map is uniquely

`[(A,S),(B,T)]_vector`
` = a0 (a cross b) + b0 (T a-S b) + d0 vee([S,T])`,

`[(A,S),(B,T)]_tensor`
` = c0 ([A,T]-[B,S])`.

The four displayed tensors are independent and exhaust the character count.

## Jacobi variety

The Jacobi identity is equivalent to

`c0(a0-c0)=0`, `a0 d0-b0^2-c0 d0=0`, `b0 c0=0`, `b0 d0=0`,
and `b0(2a0-3c0)=0`.

Over the reals this reduces exactly to `b0=0`, and either `c0=d0=0` or
`a0=c0`. The complete classification is:

- `c0=d0=0`, `a0=0`: abelian `R^8`;
- `c0=d0=0`, `a0!=0`: `so(3) direct_sum R^5`, with central `R^5`;
- `a0=c0!=0`, `a0 d0<0`: compact simple `su(3)`;
- `a0=c0!=0`, `a0 d0>0`: split simple `sl(3,R)`;
- `a0=c0!=0`, `d0=0`: `so(3) semidirect Sym_0(3)`;
- `a0=c0=0`, `d0!=0`: two-step nilpotent, with central `3'` and
  `[5,5]=3'`.

All statements above are **DERIVED** and verified in
`reproducible/verify_a5_equivariant_brackets.py`.

## Rigidity and residual family

An `A5`-equivariant change of basis acts by nonzero real scalars `r` on `3'`
and `s` on `5`:

`a0 -> r a0`, `b0 -> s b0`, `c0 -> r c0`,
`d0 -> (s^2/r)d0`.

Consequently:

- **DERIVED:** all points with `a0 d0<0` are equivariantly isomorphic to
  `(1,0,1,-1)`, the compact `su(3)` bracket;
- **DERIVED:** all points with `a0 d0>0` are equivariantly isomorphic to
  `(1,0,1,1)`, the split `sl(3,R)` bracket;
- **DERIVED:** the apparent continuous family through a compact point is only
  a block-rescaling orbit. After quotienting by equivariant isomorphism and
  overall scale, the compact-simple class has no continuous modulus.

Thus `su(3)` is the unique compact-simple equivariant bracket class, but
`A5` equivariance plus Jacobi alone does not select compact simplicity from
the other branches. Physical selection of the compact-simple class remains
**STRUCTURAL** unless compatibility with the canonical discrete metric is
required. The next section makes that condition exact.

## Canonical edge metric and compactness

The 12-dimensional gauge kernel is exactly the lift of one real amplitude on
each Hopf fiber to the alternating edge mode on that fiber's even ten-cycle.
It is not the constant lift. Hence its edge/Hodge Gram matrix is `10 I_12` in
fiber-amplitude coordinates. The arbitrary initial sign on each cycle is an
orthogonal basis convention and cannot change this metric. For the
explicit fivefold-axis orbit `n` on the icosahedral base, the color sampling
maps are

`v -> (n dot v)_n` on `3'`, and `S -> (n^T S n)_n` on `5`.

Their exact icosahedral moment identities are

`sum_n (n dot v)(n dot w) = 2 <hat(v),hat(w)>_F`,

`sum_n (n^T S n)(n^T T n) = (8/5) <S,T>_F`.

Therefore, conditional on identifying the abstract `3'+5` summand with these
explicit sampling maps:

- **DERIVED for the displayed sampling embedding:** the edge-amplitude inner
  product restricts to
  `20 <,>_F` on `3'` and `16 <,>_F` on `5`;
- **DERIVED:** the geometric Schur-scalar ratio is `lambda_3'/lambda_5=5/4`.

This is not the equal-block Frobenius representative used above for
`(1,0,1,-1)`. For the general bracket `(a0,b0,c0,d0)`, ad-invariance of the
actual edge metric is exactly

`b0=0`, `16 c0+20 d0=0`.

Intersecting these equations with the Jacobi variety leaves exactly:

- the abelian bracket;
- `so(3) direct_sum R^5` with central `R^5`;
- the compact-simple family `(a0,0,a0,-4a0/5)`, `a0!=0`, isomorphic to
  `su(3)`.

Thus:

- **DERIVED:** the split, semidirect, and nilpotent branches cannot preserve
  the canonical positive-definite edge metric;
- **DERIVED:** the geometric metric is compatible with the compact-simple
  branch, but selects its rescaled representative `d0/a0=-4/5` rather than
  `-1` in the displayed Frobenius coordinates;
- **DERIVED (conditional):** among brackets preserving the canonical metric,
  trivial center (equivalently, faithful adjoint in this finite list) uniquely
  selects `su(3)` up to overall bracket scale;
- **STRUCTURAL:** the minimal remaining color axiom is that the bracket makes
  the color sector a center-free metric Lie algebra for the already-derived
  edge inner product. The bracket classification is derived.  The verifier
  does not reconstruct the 720-edge gauge kernel and compare its actual
  isotypic intertwiners with these sampling maps, so uniqueness/canonicity of
  this embedding is not certified.  Required metric compatibility and
  center-freeness are not consequences of the current operators alone.

## Explicit `3'` embedding

Take `A=diag(-1,-1,1)` and an order-three rotation `B` about an axis whose
squared third component is `(3-sqrt(5))/6`. Then

`A^2=B^3=(AB)^5=1`, `Tr(AB)=phi'`,

and the generated rotation group has order 60. This is the repository's
`3'` convention. Embedded as real unitary matrices in `SU(3)`, its adjoint
character is `|chi_3'|^2-1=chi_3'+chi_5`.

- **DERIVED:** `ad(su(3))|A5=3'+5`, not `3+5`, in the requested convention.
- **DERIVED:** the Galois-conjugate `3` embedding analogously gives `3+5`.

There is no convention failure in the edge-kernel assignment.

## Killing normalization

At the compact point write `X=A+iS`. The bracket is the matrix commutator and
`B_su(3)(X,Y)=6 Re Tr(XY)`. Therefore

`-B_su(3)=6(<A,B>_F+<S,T>_F)`.

- **DERIVED:** the Killing form fixes the relative normalization of the `3'`
  and `5` blocks to the common defining-representation Frobenius trace. There
  is no remaining relative scale inside the color `8`, apart from one overall
  normalization of the simple factor.
- **OPEN:** relative normalization between the separate `su(3)`, `su(2)`, and
  `u(1)` factors. The direct sum permits independent invariant forms on its
  simple factors, while the abelian factor has no Killing normalization.
- **OPEN:** a common finite matter representation whose Dynkin indices compare
  the simple-factor traces, and a fixed `U(1)` charge unit.

Thus the bracket obstruction is closed, but the gauge prefactors
`(8/15,1/3,2/15)` remain **PATTERN**.

## Companion sectors

- **DERIVED:** `dim Hom_A5(Lambda^2 3,3)=1`; its nonzero bracket is the cross
  product and gives compact `su(2)=so(3)` up to scale.
- **DERIVED:** `Lambda^2(1)=0`; the `u(1)` sector is necessarily abelian and
  its norm/charge unit is not fixed by its bracket.
