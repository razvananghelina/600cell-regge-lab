# Does the current certified action generate the Hopf selector?

Date: 2026-08-10

## Verdict

**DERIVED NEGATIVE for the current certified construction.**  The repository
now contains a canonical sixth-order icosahedral potential that would select
the six Hopf axes with one sign, but none of the presently certified actions
simultaneously supplies

1. an `A5`-equivariant selected order parameter;
2. a nonzero fluctuated selector coefficient;
3. the required sign.

This is not a no-go theorem against future dynamics.  It is an exact audit of
the operators and admissible finite triples currently in the repository.

Protocol: commit `e18c00b`.  Registered verifier:
`reproducible/verify_hopf_selector_action_gate.py`.  Exact structured result:
`reproducible/hopf_selector_action_gate.json`.

## What the finite spectral-action file actually contains

`verify_spectral_action.py` certifies

```text
c0 = Tr(I),
c1 = Tr(D^2),
c2 = (1/2) Tr(D^4)
```

for the fixed unfluctuated Kahler--Dirac operator.  It contains no fluctuated
`D_A`, no order parameter and no `D_A^6` moment.  Its coefficients are finite
Taylor moments, not a specified continuum spectral action.

If `D_A` depends linearly on the vector field `n`, entries of `D_A^p` have
polynomial degree at most `p`.  Consequently the certified `D_A^2/D_A^4`
ceiling cannot produce a nonzero homogeneous degree-six anisotropy in `n`.
This scoped degree bound is **DERIVED** and independent of matrix size.

The subsequent exact projector audit found a necessary correction to the
unqualified version of that statement.  A Hopf fibration is an unoriented
axis, whose natural tensor variable is

`Q(n)=n n^T-(n.n)I/3 in Sym^2_0(R^3)`.

The same angular selector is cubic in `Q`:

`C3(Q(n))=S6(n)-(34/45)(n.n)^3`.

Therefore a fourth moment can in principle contain the selector if the
licensed fluctuation is linear in `Q`.  A sixth moment is required only under
the original linear-in-`n` hypothesis; it is not a general lower bound.

Formally, a heat trace would contain

`Tr exp(-t D_A^2) = Tr I - t Tr D_A^2 + t^2 Tr D_A^4/2
                    - t^3 Tr D_A^6/6 + ...`.

The alternating sign makes a sixth moment an interesting possible route, but
the repository has not derived its anisotropic part, the scale `t`, or a
valid order parameter.  The formal Taylor sign is not itself a selector.

The cheaper projector-valued continuation was subsequently tested on the six
theory-defined wave operators `Box_F=6A_f-A`.  Their span is an exact regular
5-simplex, but the six fibre-edge sets partition all 720 edges, implying

`sum_F Box_F=0`.

Thus the affine baseline is zero and the cubic part
`4 Tr(Box_bar X^3)` of the fourth moment vanishes identically.  The fourth
moment is even under `X -> -X` and cannot select the six positive fibration
vertices over their negatives.  This is a **DERIVED NEGATIVE** for the
canonical `Box_F` realization of the fourth-moment shortcut.

The nonzero third moment `Tr(X^3)` is not proportional to the equal-weight
projector cubic.  Its complete constrained critical set remains a separate
open gate; equality `Tr(Box_F^3)=N^2` at the desired six vertices is not by
itself a selection theorem.

## The canonical free-cell algebra does not repair the gap

On the free 600-cell cochain arena:

- the canonical left `C[2I]` placement has identically zero inner one-forms;
- the right placement has nonzero one-forms, but every enumerated real
  structure fails at least one of the `JD`, order-zero or first-order gates.

Thus that arena has no licensed fluctuated real spectral action from which a
Hopf order parameter can currently be extracted.  These are previously
**DERIVED** results, not rerun or enlarged by this audit.

## What the valid chamber counterexample contributes

The all-gate chamber witness `A=M2(C)+C^3` is a genuine finite real-triple
witness with nonzero represented one-forms.  The new exact calculation
isolates its only noncommutative-to-scalar Dirac block:

```text
block shape                       4 x 25
nonzero fixed-D incidences            12
rank(D_block)                          4
complex one-form block dimension      4
```

The represented block is exactly

`{(C tensor I_2) D_block : C in M2(C)}`.

Under left `SU(2)`, `M2(C)` is its two columns, hence two copies of the
fundamental doublet.  Its commuting multiplicity algebra is the full right
`M2(C)`, of complex dimension four.  Therefore neither column is selected by
the algebra or by `D`.  A Bloch/Hopf vector can be built from a chosen
doublet, but choosing that doublet introduces new multiplicity data.

This is **DERIVED module structure**.  Calling either column the physical
Higgs is **STRUCTURAL**.

## The decisive symmetry test

The verifier independently reconstructs the 60 rotations of `A5` on the 120
oriented chambers and applies them to the committed B1 cell colouring.  Its
four central-support capacities are distinct:

`(4,25,12,19)`.

The exact result is

```text
A5 stabilizer order = 1,
colouring orbit size = 60.
```

Hence the valid B1 representation preserves no nontrivial icosahedral
rotation.  It cannot serve as an `A5`-equivariant origin for a potential whose
six vacua are supposed to be selected only after spontaneous breaking.  The
embedding has already broken all of `A5` by construction.  Using it would
replace a six-vacuum derivation by a 60-attempt structural choice.

This does not weaken the B1 refutation: that witness was designed to disprove
a universal commutativity theorem, for which existence is sufficient.  It
does prevent promoting the same witness to a geometry-selected physical
algebra.

## Status ledger

- **DERIVED:** the current finite spectral-action certificate contains only
  fixed `D^0,D^2,D^4` moments.
- **DERIVED NEGATIVE:** these moments cannot generate homogeneous degree-six
  anisotropy from a fluctuation linear in the vector `n`.
- **DERIVED CORRECTION:** the selector is cubic in the canonical
  five-real-dimensional projector variable `Q`, so the degree argument does
  not exclude it from a fourth moment linear in `Q`.
- **DERIVED:** the six canonical `Box_F` operators realize that regular 5D
  simplex up to a fixed squared scale `10800`.
- **DERIVED NEGATIVE:** their affine centre is zero, so their fourth moment
  has no cubic term and cannot distinguish `+Box_F` from `-Box_F`.
- **DERIVED NEGATIVE:** the existing constraints `Tr(X^2)=7200` and
  `Tr(X^3)=14400` contain an exact extra regular point and a local real
  three-dimensional continuum; they do not select the six `Box_F`.
- **DERIVED:** the valid chamber noncommutative fluctuation block is
  `M2(C)`, i.e. two left fundamental doublets.
- **DERIVED NEGATIVE:** the B1 central supports have trivial `A5` stabilizer
  and a 60-element orbit.
- **STRUCTURAL:** identifying one of its two columns with a Higgs doublet or
  its Bloch vector with the Hopf direction.
- **OPEN:** an `A5`-equivariant selected finite triple with nonzero forms.
- **OPEN:** a licensed `A5`-equivariant `Q` fluctuation whose action produces
  the distinct equal-weight projector cubic, with derived coefficient and
  sign.
- **OPEN:** alternatively, a derived `D_A^6` in a vector realization.
- **OPEN:** locking the internal triplet to one handed geometric Hopf sector.

## Next admissible continuation

There are three honest routes:

1. find an `A5`-equivariant all-gate finite triple with a five-real-dimensional
   projector channel and compute the cubic restriction of its fourth moment;
2. find a vector-valued all-gate fluctuation and compute its complete sixth
   moment; or
3. introduce a geometric order parameter directly and derive its effective
   action by integrating out an already-defined operator, while proving that
   the coupling and regulator were fixed before examining the six-axis
   anisotropy.

Merely appending `-g S6` is not a derivation.  The mathematical availability
of the potential is real; its presence in the theory remains **OPEN**.
