# Spectral-action sign on the certified round--Regge path

Date: 2026-08-12

Protocol commit: `e59d1eb`

Geometric input commit: `c97349c`

Registered verifier: `reproducible/verify_round_regge_spectral_action_sign.py`

Machine-readable result:
`reproducible/round_regge_spectral_action_sign.json`

## Headline

The sign objection is closed, but the complete physical-action gate is not:

> **DERIVED CONDITIONAL SIGN.** For every standard even, nonnegative,
> nonzero spectral cutoff, the asymptotic coefficient multiplying the
> certified de Rham `A2` is strictly positive. Therefore the spectral action
> does not reverse the round preference of `A2` at leading shape-dependent
> order.

> **OPEN COMPLETE FINITE ACTION.** Positivity does not select the cutoff
> moment, the absolute scale, or the higher/singular remainder. It therefore
> does not prove that the complete action at finite cutoff has the round
> endpoint as its unique minimum.

The targeted verifier passes `15/15`. No full-suite run was performed.

## 1. Complete hypotheses

Let `D_u` be the complete ordinary de Rham operator with the already frozen
transmittal/conic domain and put `P_u=D_u^2`. Assume

```text
Tr exp(-t P_u)
 ~ (4*pi*t)^(-3/2)[A0(u)+t*A2(u)+...].
```

The metrics have already been normalized to equal volume, so `A0` is
constant. The standard bosonic spectral action is

```text
S_chi(Lambda,u)=Tr chi(D_u/Lambda),
```

where `chi` is even, nonnegative, nonzero, and sufficiently regular and
decaying. Positivity is part of the spectral-action principle, not something
derived from the 600-cell. The original formulation explicitly uses a
positive cutoff function: [Chamseddine--Connes, *The Spectral Action
Principle*](https://arxiv.org/abs/hep-th/9606001).

## 2. Exact Mellin sign

Set `F(x)=chi(sqrt(x))`. Mellin functional calculus gives in dimension three

```text
S_chi(Lambda,u)
 ~ (4*pi)^(-3/2)[
      Lambda^3 C0(chi) A0(u)
    + Lambda   C2(chi) A2(u)
    + ...],

C0(chi)=4/sqrt(pi) integral_0^infinity chi(v)v^2 dv,
C2(chi)=2/sqrt(pi) integral_0^infinity chi(v) dv.
```

For a nonzero nonnegative cutoff,

```text
C0(chi)>0,  C2(chi)>0.
```

The continuous certificate proves

```text
Delta A2(u)=A2(u)-A2(1)>0,  u<1.
```

Hence the leading shape-dependent action difference is

```text
Delta S
 ~ (4*pi)^(-3/2) Lambda C2(chi) Delta A2 > 0.
```

At every fixed `u<1`, round therefore wins asymptotically, provided the
stated expansion and remainder apply. The sign is not fitted after seeing
the path result.

## 3. Why the weight and finite action remain free

The exact positive family

```text
chi_a(v)=exp(-a v^2),  a>0,
```

has

```text
C0=a^(-3/2),
C2=a^(-1/2),
C2/C0=a.
```

Thus positivity permits a continuum of relative weights. For example,
`a=1` gives `C2=1`, while `a=4` gives `C2=1/2`. This is the same unselected
dimensionless heat-time freedom previously exposed by the scale audit.

Dropping positivity would allow `-exp(-v^2)` and would reverse the sign, but
that function is outside the standard bosonic spectral-action hypothesis.
Accordingly the positive sign is universal spectral-action structure; it is
not a special consequence of `a1=5`.

Positivity also does not fix finite-cutoff ordering. The preregistered exact
control spectra

```text
X={0,1,10},  Y={0,2,3}
```

reverse their positive heat-trace ordering between
`exp(-t)=100/101` and `1/2`. These are deliberately generic control spectra,
not spectra assigned to the round--Regge path. They disprove the inference
“positive cutoff implies cutoff-independent minimizer”.

The path verifier computes no `A4`, higher singular coefficient, eigenvalue
family, or uniform finite-cutoff remainder. Consequently it cannot establish
the complete finite action merely by having certified `A2`.

## 4. Status ledger

| Claim | Status |
|---|---|
| Mellin coefficient of `A2` is `Lambda C2/(4*pi)^(3/2)` | **DERIVED** |
| Every standard positive cutoff has `C2>0` | **DERIVED FROM THE POSITIVITY AXIOM** |
| Leading spectral-action term preserves the certified round-path ordering | **DERIVED CONDITIONAL SIGN** |
| This positive sign is special to the 600-cell or `a1=5` | **REFUTED** |
| Positivity fixes `C2/C0`, `chi`, or `Lambda` | **REFUTED** |
| Positive heat traces have cutoff-independent ordering | **REFUTED IN GENERAL** |
| At each fixed `u<1`, sufficiently high cutoff prefers round | **STRUCTURAL / ASYMPTOTIC, conditional on the remainder** |
| Complete finite-cutoff action uniquely minimizes at round on the full path | **OPEN** |
| Newton's constant or a Planck scale follows from the sign | **REFUTED** |

## 5. Consequence

The earlier phrase “the physical sign still has to be supplied” was too
coarse. Under the standard spectral-action principle, positivity supplies
the sign and it is the favorable one. What remains genuinely missing is the
selected weight and dominance of the complete action.

The next hostile gate was subsequently completed in
`round_a2_transverse_hessian_result.md`. Non-gauge conformal `l=2` modes do
have the opposite sign, so the round point is a saddle of this `A2` on the
full smooth metric space even though it remains minimizing in the certified
Hopf and affine Regge sectors. The later finite audit also finds 150 negative
directions at the equilateral Regge point. Therefore the positive asymptotic
weight preserves, rather than repairs, both smooth and finite `A2`
instabilities; the complete finite-cutoff action remains open.
