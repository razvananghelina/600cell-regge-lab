# Preregistration: full-Hessian spectral route to Hopf selection

Date: 2026-08-10

## Fixed arena and framing

Retain every off-diagonal channel of the already audited full label Hessian

```text
H_X(i,j)=3 Tr(X(Box_i Box_j+Box_j Box_i)),
X in W=span_R{Box_i},
Tr(X^2)=7200.
```

The constant label vector is a universal zero mode.  Let `Hhat_X` be the
restriction of `H_X` to the canonical physical subspace `1^perp`.  No basis
of `1^perp` may enter a final invariant.

This is a possible bypass of label superselection, not a derivation of it.
Even a positive result remains **STRUCTURAL** until the extended cubic and
its label fluctuation field are licensed as physical dynamics.

## Complete primitive spectral list

For the five-dimensional operator `Hhat_X`, freeze the characteristic
polynomial convention

```text
det(lambda I-Hhat_X)
  = lambda^5-e1(X)lambda^4+e2(X)lambda^3
    -e3(X)lambda^2+e4(X)lambda-e5(X).
```

The exact number of primitive spectral coefficients is therefore

```text
N=5.
```

Also compute the power sums `s_k=Tr(Hhat_X^k)`, `1<=k<=5`, solely as an
independent Newton-identity check.  They are not five additional attempts.
All higher single traces and `det(Hhat_X^2)=e5^2` are polynomial functions of
the same five coefficients and must not inflate `N`.

An arbitrary linear combination, heat trace at a selected `t`, cutoff
function, or fitted polynomial in the `e_k` is forbidden.  Without a derived
regulator those form an infinite look-elsewhere family.

## Two-stage provenance protocol

### STEP 1: blind spectral enumeration

Using exact integer/rational arithmetic:

1. construct `Hhat_X` on five independent coordinates of `W`;
2. record all five `e_k(X)` and all five `s_k(X)` coefficientwise;
3. verify Newton identities and `A5` invariance;
4. record which coefficients vanish, which are functions only of the norm,
   and the number of algebraically distinct nonconstant primitive spectral
   characters on the normalized sphere;
5. write the complete enumeration to a machine-readable file.

Do **not** evaluate at `Box_i`, compare with `Tr(X^3)`, the equal-weight
projector cubic `C_box`, the six desired fibrations, or any signed target.
Commit this enumeration with an explicit no-target-comparison message.

### STEP 2: target comparison only after the STEP 1 commit

Then, and only then:

1. express each nonconstant primitive coefficient in the already known
   invariant basis, including `Tr(X^3)` and `C_box(X)` where degrees permit;
2. test exact constrained stationarity and tangent Hessian at all signed
   `Box_i`;
3. test whether an individual primitive coefficient, without an added
   coefficient or sign chosen after inspection, has exactly the desired
   global extremal/kernel locus;
4. report the hit count as a fraction of `N=5` and separately among the
   distinct nonconstant normalized-sphere characters;
5. audit the Gaussian interpretation: a bosonic Euclidean Hessian integral
   requires positive definiteness; a determinant of an indefinite Hessian is
   not silently called a convergent bosonic effective action.

## Acceptance and kill boundaries

- **Structural advance:** an individual primitive spectral coefficient of
  the complete `Hhat_X`, with every off-diagonal channel retained, uniquely
  selects the signed Hopf simplex or, together with the already fixed positive
  cubic condition, exactly the six `+Box_i`.
- **Pattern:** a desired orbit is merely stationary, one hit appears among
  several attempts/competitors, or a sign/regulator must be selected after
  comparison.
- **Kill for this bypass:** the primitive spectral data reduce to invariants
  already known not to select, the desired orbit is not extremal, exact
  competitors remain, or the only apparent selector needs a fitted spectral
  function.

## Attack on the protocol itself

Characteristic coefficients are a complete finite census of spectral data,
but they do not make every function of the spectrum a separate canonical
physical action.  The enumeration closes parameter-free primitive
single-operator evidence only.  Preregistration prevents adaptation after
this point; it cannot erase knowledge from the earlier Hessian calculation,
including the already known spectrum at one `Box_i`.  Any positive result
must therefore be labelled with that provenance limitation.
