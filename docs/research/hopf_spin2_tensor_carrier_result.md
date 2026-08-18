# The Hopf geometry contains a genuine round-`S^3` spin-two tensor carrier

Date: 2026-08-12  
Protocol commit: `4004c25`  
Registered verifier: `reproducible/verify_hopf_spin2_tensor_carrier.py`  
Machine-readable result: `reproducible/hopf_spin2_tensor_carrier.json`

## Headline

**DERIVED ROUND-`S^3` KINEMATIC POSITIVE.**  The repository already contained
the ingredients of the correct local tensor type for a gravitational field,
but they had not been assembled in that role.

For the six certified fivefold-axis projectors

```text
T_i = P_i - I/3,
```

left or right quaternion multiplication lifts the `T_i` to smooth symmetric
tracefree tangent tensors on the unit round three-sphere.  At every point they
span the complete five-dimensional fibre

```text
Sym^2_0(T*S^3).
```

With constant coefficients in either handed invariant frame, every such
tensor is transverse as well as tracefree and obeys

```text
nabla* nabla H = 6 H.
```

The left- and right-invariant homogeneous spaces each have dimension five and
intersect only at zero.  Thus the two handed Hopf families supply two exact
five-dimensional homogeneous TT eigenspaces.

This is a real kinematic advance.  It is not yet a graviton, a gravitational
action, or emergent general relativity.

## Provenance and scope

Protocol commit `4004c25` is post-recognition: the TT expectation and Casimir
value were understood before the verifier was written.  The commit freezes
the formulas and falsifiers but is not evidence of a blind discovery.

The complete metric hypothesis is the **unit round** `S^3=SU(2)`.  The result
does not silently transfer to the distinct fixed piecewise-flat Regge metric
of `whitney_regge_continuum_transfer_result.md`.

No mass, coupling, speed, Newton constant, Planck scale, or measured target
was used.

## 1. Exact tensor frame

The six exact projectors satisfy

```text
sum_i T_i = 0,
Tr(T_i^2) = 2/3,
Tr(T_i T_j) = -2/15  for i != j,
rank span{T_i} = 5.
```

Their frame operator on `Sym^2_0(R^3)` is

```text
sum_i |T_i><T_i| = (4/5) I_5.
```

Consequently every symmetric tracefree `3 x 3` tensor has the exact,
choice-free reconstruction

```text
H = (5/4) sum_i Tr(H T_i) T_i.
```

The six axes are therefore an overcomplete tight frame, not six adjustable
Schur coefficients.

## 2. Global quaternionic lift

For a unit quaternion `q` and `v in Im(H)`, define

```text
L_q(v)=qv,       R_q(v)=vq.
```

Exact polynomial calculation gives

```text
L_q^T q = R_q^T q = 0,
L_q^T L_q = R_q^T R_q = |q|^2 I_3.
```

Thus both maps are orthonormal identifications of `R^3` with `T_qS^3` when
`|q|=1`.  The lifts

```text
T_i^L(q)=L_q T_i L_q^T,
T_i^R(q)=R_q T_i R_q^T
```

are symmetric, tangent and tracefree for every `q`.  Because the tangent map
is an isometry and the `T_i` span the identity fibre, either handed family
spans all of `Sym^2_0(T_q^*S^3)` at every point.

This corrects a potentially damaging framing error: the Hopf order parameter
is not merely an internal abstract `A5` module.  On the round quaternionic
carrier it has a canonical realization as the symmetric-tracefree tangent
tensor bundle.

## 3. Transversality and connection spectrum

Use a left-invariant unit frame with

```text
[e_a,e_b] = 2 epsilon_abc e_c,
nabla_(e_a)e_b = epsilon_abc e_c.
```

Let `Gamma_a` be the corresponding skew `so(3)` matrices.  A tensor with
constant left-frame components satisfies

```text
nabla_a H = [Gamma_a,H].
```

For symmetric `H`, its divergence is

```text
(div H)_b = sum_a [Gamma_a,H]_(ab) = 0,
```

because the connection tensor is antisymmetric while `H` is symmetric.
Tracefreeness is already exact.  Therefore every constant-frame element of
the five-space is TT.

The connection Laplacian is the spin-two Casimir:

```text
nabla* nabla H
  = -sum_a [Gamma_a,[Gamma_a,H]]
  = 6 H.
```

The verifier proves this coefficientwise for a generic symbolic
symmetric-tracefree matrix, then checks all six Hopf tensors.  Reversing the
connection sign in the right-invariant frame leaves both divergence and the
double commutator unchanged.

For a round sphere of radius `R`, dimensional rescaling would give `6/R^2`.
The theory has not selected a physical `R`, so this is a scaling statement,
not a mass or frequency prediction.

## 4. The two handed spaces

If a tracefree tensor field were constant in both the left and right frames,
its identity-frame tensor would be invariant under the adjoint rotation group.
The verifier constructs exact order-three and order-five icosahedral rotations,
checks that they generate all 60 elements of `A5`, and solves the invariance
equations on `Sym^2_0(R^3)`.  Their rank is five, so the invariant subspace is
zero.

Hence

```text
dim H_TT^L = 5,
dim H_TT^R = 5,
H_TT^L intersect H_TT^R = {0}.
```

Under `Spin(4)=SU(2)_L x SU(2)_R`, these are the two handed homogeneous
spin-two modules `(0,2)` and `(2,0)`, up to the convention that names left and
right.  Calling the two modules physical helicities would require Lorentzian
dynamics and a gauge quotient and is not claimed.

## 5. Why the old coexact argument was insufficient

Under proper local tangent rotations,

```text
Lambda^*(R^3)
 = Lambda^0 + Lambda^1 + Lambda^2 + Lambda^3
 = 2 V_0 + 2 V_1.
```

The geometric fibre of the Kähler--Dirac/de Rham field therefore contains
scalars and vectors, but no intrinsic `V_2`.  A scalar or vector field can
have an orbital spherical harmonic with `l=2`; that does not make its local
tensor type spin two.  This is why the earlier coexact/scalar multiplicity
ratio could not establish a graviton.

By contrast,

```text
Sym^2(V_1) = V_0 + V_2,
Sym^2_0(V_1) = V_2.
```

The Hopf projector fields realize precisely this missing tracefree symmetric
piece.  They are a separate geometric tensor construction, not a hidden
component of the exterior-algebra fibre.

## 6. What this does and does not establish

**DERIVED:** a canonical six-element tight frame for
`Sym^2_0(T^*S^3)` on the round quaternionic carrier.

**DERIVED:** five homogeneous left-TT and five homogeneous right-TT fields,
each with connection-Laplacian eigenvalue six and with zero intersection.

**DERIVED CORRECTION:** the repository does possess the correct intrinsic
spin-two tensor carrier.  The statement that only cochain 1-forms were
available was too narrow.

**STRUCTURAL ADVANCE:** these fields are credible kinematic seeds for an
emergent gravitational sector.

**OPEN:** the existing certified Kähler--Dirac action has no selected coupling
to this five-component label/tensor sector.  The prior action-origin audit
explicitly records that absence.

**OPEN:** no variable metric, quadratic gravitational action, curvature shift
or Lichnerowicz operator has been selected.  Importing the desired Einstein or
Regge action would establish compatibility, not emergence from the theory.

**OPEN:** no diffeomorphism constraint, Lorentzian propagation, massless
quotient, universal stress-tensor coupling, nonlinear closure, Newton constant
or Planck scale follows.

**SCOPE BOUNDARY:** the smooth round metric and the fixed Regge continuum are
uniformly equivalent for the previously studied de Rham problem, but they are
not the same metric.  A tensor-action result on one cannot be assigned to the
other without a new transfer theorem.

## 7. Next strict gate

The next admissible question is no longer “does a spin-two carrier exist?”;
it does.  The next question is:

> Does the theory's already-selected metric dependence of the Kähler--Dirac
> or Whitney operator induce, without inserted coefficients, a nondegenerate
> quadratic action on the Hopf `Sym^2_0` tensor fields, and does that action
> extend naturally under the canonical refinement tower?

The action family and normalization must be frozen before inspecting its
spectrum.  A positive finite Hessian would still be stiffness; propagation
requires the later Lorentzian and source-coupling gates.  Failure to derive
the coupling would leave the tensor carrier kinematically real but physically
decoupled.

## Status ledger

| Claim | Status |
|---|---|
| Six centered projectors span `Sym^2_0(R^3)` | **DERIVED** |
| Left/right lifts span `Sym^2_0(T_q^*S^3)` at every `q` | **DERIVED** |
| Homogeneous lifted tensors are TT | **DERIVED** |
| Connection-Laplacian eigenvalue is `6` on unit round `S^3` | **DERIVED** |
| Left/right homogeneous intersection is zero | **DERIVED** |
| Exterior/Kähler--Dirac fibre contains intrinsic spin two | **REFUTED** |
| Hopf tensor construction supplies intrinsic spin two | **DERIVED** |
| These modes are physical gravitons | **OPEN** |
| Existing spectral action couples to them | **DERIVED NEGATIVE as currently defined** |
| Result transfers unchanged to the fixed Regge metric | **OPEN / not established** |
| Newton/Planck scale is selected | **OPEN / unsupported** |

## Reproduction

Run only the targeted verifier:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_hopf_spin2_tensor_carrier.py
```

Expected result: `20/20 checks passed`.  The full suite is not run under the
user's current instruction.
