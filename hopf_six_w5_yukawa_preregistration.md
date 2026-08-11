# Preregistration: can the derived `W5` order parameter supply the missing Dirac tensor?

Date: 2026-08-11

## Disclosed motivation and prior information

Commit `22717d8` proves that every odd first-order Dirac operator on the
936-state carrier that is natural as an `A5`-invariant scalar has a
non-scalar algebra commutant.  Static invariant incidence is therefore
exhausted on this carrier.

The repository independently contains a real five-dimensional `A5` order
parameter

```text
W5 = span{Box_i-Box_j} = 1^perp subset R^6.
```

Its six centered Hopf projectors form a regular 5-simplex.  The complete
label-Hessian affine fourth moment has two sign branches and, structurally,
selects one of the two signed six-point orbits for every nonzero magnitude
ratio.  This is disclosed prior work, not a blind discovery.  The action is
not yet a licensed sector of the finite spectral triple, and the relative
sign is not derived.

The candidate continuation is a covariant Yukawa/Dirac field rather than an
invariant scalar:

```text
T_ij : W5 -> Hom_R(V_i,V_j),
T_ij(g X) = rho_j(g) T_ij(X) rho_i(g)^-1.
```

At a nonzero vacuum `X`, `T_ij(X)` may see irreducible sectors that every
invariant intertwiner misses.  No coupling multiplicity or connectedness
result has been computed before this protocol.

## Phase 1: target-blind covariant-coupling census

Using the exact real `A5` character table and

```text
V0=1+5,
V1=3+3',
V2=3+4+5,
V3=3'+4+5,
W5=5,
```

compute for all six unordered node pairs

```text
m_ij = dim Hom_A5(W5, Hom_R(V_i,V_j))
     = <chi_5, chi_i chi_j>_A5.
```

Before evaluating a single Hopf vacuum or a connectedness rank, record and
commit:

1. the complete `4 x 4` symmetric multiplicity matrix;
2. the multiset of the six off-diagonal multiplicities;
3. for each of the eight already derived spectral readings, the three legal
   central links and their coupling multiplicities;
4. the number of readings for which every required link has multiplicity
   exactly one;
5. the number with a zero link and the number requiring a choice in a
   higher-dimensional coupling space.

Multiplicity zero kills a link.  Multiplicity one selects a covariant tensor
line up to scale.  Multiplicity greater than one is **STRUCTURAL freedom**;
no favorable line may be selected after seeing connectedness.

This Phase-1 file must be committed with an explicit statement that no
vacuum matrix, Dirac commutant or physical target has been inspected.

## Phase 2: exact vacuum evaluation, only after Phase 1 is committed

If at least one reading has multiplicity-one covariant lines on every legal
link, construct them independently in an exact model of the actual
six-fibration permutation action.  The construction must solve the covariance
equations for the full group; a character multiplicity alone does not supply
Clebsch--Gordan coefficients.

For every eligible reading and all twelve signed simplex points `+/-Box_i`:

1. evaluate every normalized covariant tensor `T_ij(X)`;
2. record its exact rank and singular-value polynomial;
3. insert it into all first-order-compatible `T tensor I` and `I tensor T`
   cell positions, with adjoints forced by self-adjointness and cell
   transpose forced by `J`;
4. check oddness, `JD=DJ`, exact first order and nonzero inner one-forms;
5. compute exactly the dimension of
   `{a in B_R : [D_X,pi(a)]=0}`;
6. verify that all six points in one orbit give conjugate results, rather
   than counting six symmetry copies as six attempts.

Independent nonzero scalar magnitudes on different links are not selected by
representation theory.  Connectedness is insensitive to such magnitudes, so
the maximal-support commutant may be used as a necessary gate.  A positive
spectral or mass claim may not use arbitrarily chosen relative magnitudes.

If no reading has three unique lines, Phase 2 may still compute the common
commutant of the **entire** covariant coupling span.  If even that maximal
span is nonconnected, the whole linear `W5` route is killed.  If it is
connected, the result is only a structural opening until geometry selects a
line in every multiplicity-greater-than-one space.

## Acceptance and kill boundaries

- **DERIVED COVARIANT BRIDGE:** at least one preregistered reading has unique
  covariant tensor lines on all links and every nonzero simplex vacuum yields
  a connected, real, odd, first-order `D_X` with nonzero forms.
- **STRUCTURAL OPENING:** the full covariant span can connect the carrier, but
  at least one required tensor line or relative coupling remains unselected.
- **DERIVED LINEAR-FIELD NO-GO:** even the complete legal `W5`-covariant span
  has a non-scalar algebra commutant for all eight readings.

Even the first outcome would not yet be a complete physical vacuum theorem:
the existing affine Hessian action is a structural extension, not a derived
inner fluctuation of the 936-state triple, and its relative sign has a `1/2`
look-elsewhere ambiguity.  These provenance limitations must remain in the
verdict.

No matter character, mass, coupling, Hessian target value or Standard-Model
comparison is allowed in Phase 1.  Only targeted verifiers will be run; the
full suite remains excluded by user instruction.
