# The canonical `C^6` label algebra does not derive viable superselection

Date: 2026-08-10

Protocol commit: `39e35b5`.  Registered verifier:
`reproducible/verify_hopf_label_superselection.py`.  Targeted result: `15/15`.

## Complete question

Let

```text
A_F=C^6
```

be the canonical algebra of functions on the six derived Hopf fibrations.
Can order zero and first order force the label-diagonal response used by
`Phi`, while preserving connectedness, nonzero one-forms, `A5` equivariance,
KO6 and metric-dimension-zero orientability?

The answer is no in two complementary senses:

1. the minimal representation gets diagonality only by killing all dynamics;
2. the canonical pair bimodule retains off-diagonal dynamics through order
   zero, first order, connectedness and KO6, while adding orientability kills
   the entire `A5`-equivariant `C^6` arena rather than selecting a viable
   diagonal sector.

## Minimal representation: locality without dynamics

Use

```text
H_min=C^6,
pi(a)=diag(a_0,...,a_5),
J=complex conjugation.
```

For a general matrix element `D_pq`, first order gives factors of the form

```text
D_pq (a_q-a_p)(b_q-b_p).
```

Testing the complete 36-matrix-unit basis against all 36 pairs of minimal
projections leaves exactly

```text
6 diagonal units,
0 off-diagonal units.
```

Thus first order forces `D` diagonal.  But then

```text
[D,A_F]=0,
Omega_D^1(A_F)=0,
{a in A_F:[D,a]=0}=A_F != C*1.
```

Connectedness fails maximally and all represented inner one-forms vanish.

The same applies to the doubled selector representation

```text
pi(a)=diag(a,a),
D_aux=[[I,Phi],[Phi,I]],
```

because every block is label diagonal.  It commutes with the whole `C^6`.

**DERIVED NEGATIVE:** minimal first-order locality cannot license `D_aux` as
a fluctuated connected Dirac operator.  It obtains the desired diagonal only
by removing the dynamics the diagonal was supposed to explain.

## Canonical pair bimodule: first order permits off-diagonal dynamics

Use the full pair-groupoid carrier

```text
H_pair=C^6 tensor C^6,
L(a)|i,j>=a_i|i,j>,
R(b)|i,j>=b_j|i,j>,
J_pair|i,j>=|j,i> K.
```

The parameter-free rook operator is

```text
D_pair=(A_K6 tensor I)+(I tensor A_K6).
```

It changes exactly one label coordinate at a time.  Exact matrix tests give:

```text
[L(a),R(b)]                         =0,
[[D_pair,L(a)],R(b)]                =0,
D_pair*=D_pair,
J_pair D_pair                       =D_pair J_pair.
```

The complete first-order legal support consists of

```text
396 matrix units total,
36 diagonal,
360 genuinely off-diagonal.
```

`D_pair` uses all 360 legal off-diagonal edges.  It is equivariant under all
60 elements of the derived `A5` action.

For the left algebra,

```text
rank(a -> [D_pair,L(a)])=5,
dim kernel=1,
dim_C Omega_D^1(A_F)=30.
```

Therefore it is connected and fluctuating.

The standard odd double

```text
D=[[0,D_pair],[D_pair,0]],
gamma=diag(I,-I),
J=[[0,J_pair],[J_pair,0]] K
```

has exactly

```text
J^2=+1,
JD=+DJ,
Jgamma=-gamma J,
gamma D=-D gamma.
```

Order zero and first order survive the double.

**SCOPED REFUTATION:** order zero, first order, connectedness, nonzero forms,
KO6 and `A5` equivariance do not imply fibration-label superselection.  A
canonical counterexample has 360 off-diagonal channels.

## Why this witness is not yet an all-gate counterexample

Metric-dimension-zero orientability asks for

```text
gamma in span{pi(a)Jpi(b)J^-1}.
```

On the odd double, every such zero-cycle acts identically on the two grading
sheets, while `gamma` has opposite signs.  Exact diagonal-vector ranks are

```text
rank zero-cycle span             =36,
rank after adjoining gamma       =37.
```

Hence this particular pair witness is not orientable.

## General orientability no-go for the whole arena

The failure is not repaired by changing multiplicities.

Any finite `C^6` bimodule decomposes as

```text
H = direct_sum_{i,j} H_ij,
```

where the left and opposite algebras act on `H_ij` by the scalar pair
`(a_i,b_j)`.  A metric-dimension-zero orienting cycle therefore acts on every
`H_ij` as

```text
epsilon_ij I_{m_ij},       epsilon_ij in {+1,-1}.
```

KO6 reality requires

```text
epsilon_ji=-epsilon_ij.
```

The derived `A5` action has exactly two orbits on the 36 ordered label pairs:

```text
diagonal pairs                         orbit size 6,
all ordered distinct pairs             orbit size 30.
```

Both are obstructed:

- on a diagonal block, KO6 demands
  `epsilon_ii=-epsilon_ii`;
- the off-diagonal orbit contains both `(i,j)` and `(j,i)`.  `A5` invariance
  demands equal signs, while KO6 demands opposite signs.

Extra multiplicity does not help because every zero-cycle remains scalar on
the multiplicity space of a fixed bimodule block.

Therefore:

> **DERIVED NO-GO.** Under the complete hypotheses
> `A=C^6`, `A5`-equivariant bimodule, KO6 and metric-dimension-zero
> orientability, no nonzero finite real bimodule exists, even with arbitrary
> multiplicities.

In particular, there is no faithful connected all-gate triple in this arena.
Adding orientability does not derive the desired diagonal; it deletes the
whole arena.

## Consequence for the Hopf selector

The three outcomes are now exhaustive for this algebraic route:

1. **Minimal representation:** diagonal, but disconnected and zero forms.
2. **Canonical pair representation:** connected, fluctuating and first-order,
   but retains all off-diagonal label dynamics and fails orientability.
3. **Full KO6 plus metric-zero orientability:** no nonzero `A5`-equivariant
   `C^6` bimodule at all.

Thus the canonical function algebra on the six fibrations cannot turn the
Gram recognizer into a licensed physical selector under the current axioms.

## Scope boundaries

The no-go depends on every stated hypothesis.  It does not cover:

- breaking `A5` before imposing orientability;
- a different KO sign with `Jgamma=+gamma J`;
- positive metric-dimension Hochschild orientability;
- a noncommutative or crossed-product label algebra;
- abandoning the six-point label carrier.

These are genuine changes of arena, not loopholes inside the proved result.

## Status ledger

- **DERIVED NEGATIVE:** minimal first order forces diagonal `D` and therefore
  zero one-forms plus failed connectedness.
- **SCOPED REFUTATION:** the pair-groupoid witness retains 360 off-diagonal
  channels while passing order zero, first order, connectedness, nonzero
  forms, KO6 and `A5` equivariance.
- **DERIVED NEGATIVE:** that witness fails metric-zero orientability.
- **DERIVED NO-GO:** no nonzero `A5`-equivariant KO6 metric-zero orientable
  `C^6` bimodule exists, arbitrary multiplicities included.
- **DERIVED NEGATIVE:** `C^6` does not license the diagonal selector under the
  current complete axiom set.
- **OPEN:** a separately motivated noncommutative label algebra or different
  KO/metric dimension.

## Next admissible route

If the label programme is continued, the most canonical changed arena is the
transformation-groupoid/crossed-product algebra

```text
C(F) crossed_product A5,
```

because it contains both label functions and the already derived symmetry
transitions.  Its Wedderburn type and orientability must be computed before
asking whether it contains the selector.  Choosing a smaller noncommutative
algebra after inspecting the desired blocks would reintroduce fitting.
