# Preregistration: do spectral axioms derive fibration-label superselection?

Date: 2026-08-10

## Exact question and scope

Let `A_F=C^6` be the canonical algebra of complex functions on the six derived
Hopf fibrations.  Does order zero plus first order force every admissible
operator response to be diagonal in the fibration labels while retaining
connectedness and nonzero one-forms?

The primary claim tested is deliberately scoped to order zero, first order,
connectedness, nonzero forms, `A5` equivariance and the stated KO signs.  A
counterexample that fails orientability does **not** refute a theorem with
orientability added; its failure must be recorded separately.

## Arena A: minimal label representation

Use

```text
H_min=C^6,
pi(a)=diag(a_0,...,a_5),
J=complex conjugation.
```

For an arbitrary self-adjoint `D`, derive the complete matrix-element form of

```text
[[D,pi(a)],J pi(b) J^-1]=0.
```

Test whether this forces `D` diagonal.  If it does, also test connectedness

```text
{a in A_F:[D,pi(a)]=0}=C*1
```

and the represented one-form space.  Locality that is obtained only by making
all commutators zero is not a physical advance.

Apply the same audit to the doubled selector carrier with
`pi(a)=diag(a,a)`: its `Phi(X)` acts label-diagonally and must not be called a
fluctuated Dirac if its one-forms vanish.

## Arena B: canonical pair-groupoid bimodule

Use the complete ordered-pair carrier

```text
H_pair=C^6 tensor C^6,
L(a)|i,j>=a_i|i,j>,
R(b)|i,j>=b_j|i,j>,
J_pair|i,j>=|j,i> with complex conjugation.
```

Define the parameter-free rook-graph operator `D_pair` connecting two pairs
iff exactly one coordinate changes.  It is the sum of the two complete-graph
adjacencies on the factors; no edge weights are fitted.

Freeze the following tests:

1. order zero `[L(a),R(b)]=0` on all minimal projections;
2. first order `[[D_pair,L(a)],R(b)]=0` on all minimal projections;
3. self-adjointness and `J_pair D_pair=D_pair J_pair`;
4. diagonal `A5` equivariance under the already derived six-label action;
5. connectedness for the left algebra and nonzero represented one-forms.

Then form the standard odd double

```text
H=H_pair+H_pair,
gamma=diag(+I,-I),
D=[[0,D_pair],[D_pair,0]],
J=[[0,J_pair],[J_pair,0]] K.
```

Test the KO6 signs `(J^2,JD,Jgamma)=(+1,+1,-1)` and preserve all previous
gates.

## Orientability boundary

Independently determine whether metric-dimension-zero orientability can hold
for the odd double with the stated representation.  The span of

```text
pi(a) J pi(b) J^-1
```

must be computed exactly.  If it acts identically on the two grading sheets,
it cannot contain `gamma`; record orientability failure rather than hiding it.

Also test whether a single-copy KO6 grading on ordered distinct pairs can be
both `A5` invariant and satisfy `Jgamma=-gamma J`.  A transitive ordered-pair
orbit containing the swap would obstruct such a sign assignment.

## Decision boundary

- **Superselection advance:** the axioms force the diagonal label algebra in
  an arena that remains connected, fluctuating and passes the stated KO and
  orientability gates.
- **Scoped refutation:** the pair-groupoid witness passes order zero, first
  order, connectedness, nonzero forms, KO6 and `A5` equivariance while
  retaining off-diagonal label channels.  Then those gates alone do not imply
  superselection.
- **Full-gate open:** if that witness fails orientability, no conclusion may be
  claimed under orientability.  State whether the obstruction kills this
  arena or actually proves a wider no-go.

The aim is not to find any representation that gives the desired diagonal.
It is to learn whether the diagonal is forced without sacrificing the very
dynamics it was meant to license.
