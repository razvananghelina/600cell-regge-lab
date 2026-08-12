# Preregistered protocol: can the theory select the Lorentzian tent volume term?

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- PRELIMINARY NON-SELECTION DISCLOSED**

## 1. Exact question

The registered symmetric Lorentzian tent audit found no stationary vacuum
pole at zero volume coefficient. On the static branch, a generic action

```text
S_tent = sum_h A_h epsilon_h - lambda sum_sigma V_sigma
```

would require a positive dimensionless coefficient

```text
ell=lambda a^2
```

that depends on the desired proper-time ratio `rho=tau^2/a^2`.

This audit asks a narrower question before any nonsymmetric search:

> Does the already certified theory, without consulting a desired `rho`,
> select the sign and one unique value of `ell` on the same Lorentzian
> four-dimensional tent carrier?

Only constructions committed before this protocol may count as selectors.
No measured cosmological constant, Planck scale, desired tick, Euclidean
golden ratio or numerical coincidence may be used.

## 2. Acceptance standard

An existing construction selects the coefficient only if it supplies every
one of the following:

1. an action on the Lorentzian four-dimensional tent carrier, or a proved
   transfer from its own carrier to that action;
2. both curvature and four-volume terms in one convention;
3. their relative sign and normalization, modulo one irrelevant common
   nonzero action multiplier;
4. a unique dimensionless `ell=lambda a^2`;
5. if a spectral action is used, a selected cutoff function and a selected
   cutoff-to-edge ratio;
6. if the source is a spatial action, a derived `3+1`/lapse/extrinsic-curvature
   transfer rather than merely matching the algebraic form of two terms.

Equality with one value required by the tent equation is not evidence unless
items 1--6 were fixed first. Choosing a cutoff after selecting `rho` is
fitting.

## 3. Frozen candidate census

The verifier must inspect the authoritative certificates for every current
candidate class relevant to this question.

### A. Fixed finite Kähler--Dirac moments

The fixed 2640-state operator has

```text
c0=2640, c1=14880, c2=55920,
r=c1/(2c0)=31/11.
```

The dimension audit says these are Taylor moments of a finite matrix, not
continuum Seeley--DeWitt coefficients. The metric-selector certificate says
that under scaling `r` is an inverse-length-squared unit, but

```text
t=alpha/r,       alpha>0
```

leaves a free dimensionless heat time. The audit must verify whether any
committed result fixes `alpha` or maps these moments to a four-volume term.

### B. Smooth/conical three-dimensional de Rham spectral action

The certified continuum ordinary de Rham coefficients on the spatial
three-sphere are

```text
A0=8 Vol_3,
A2=-(2/3) integral_S3 R_3 dVol_3.
```

For a standard positive even cutoff,

```text
S_chi ~ (4*pi)^(-3/2)
        [Lambda^3 C0 A0 + Lambda C2 A2 + ...],

C0>0, C2>0.
```

After one common rescaling that normalizes the spatial curvature coefficient
to `+1`, derive rather than assume

```text
S_normalized = integral R_3 dVol_3 - lambda_3 Vol_3 + ...,
lambda_3 = 12 Lambda^2 C0/C2.
```

For the fully admissible Gaussian family

```text
chi_s(v)=exp(-s v^2),       s>0,
C0=s^(-3/2), C2=s^(-1/2),
```

the dimensionless spatial ratio is therefore

```text
lambda_3 a^2 = 12 (Lambda a)^2/s.
```

This family is a mandatory hostile control. It tests whether positivity fixes
only the favorable relative sign while leaving every positive magnitude
available.

The audit must not silently identify this three-dimensional potential with
the four-dimensional Lorentzian Regge action. It must check whether the
repository contains a derived extrinsic-curvature/lapse transfer.

### C. Four-dimensional carriers and other registered gravity operators

Inspect the certificates for:

- the canonical CW slab `K x I`;
- the Lorentzian tent result;
- the finite `Box` spectral/stiffness constructions;
- the registered electromagnetic `alpha` equation.

An operator that gives a static stiffness, a gauge coupling, a carrier count
or a dimensionless spectral invariant does not select `ell` unless it also
supplies the action bridge in items 1--6.

Historical exploratory scripts that explicitly tune a heat time or compare
to a cosmological target are not certified selectors. Their existence must
be reported, not promoted.

## 4. Frozen mathematical checks

The verifier must establish:

1. the exact `A0,A2` multipliers and the normalized factor `12`;
2. the Gaussian moment identities;
3. `lambda_3 a^2` ranges over all positive real values as the effective
   cutoff `(Lambda a)/sqrt(s)` varies;
4. the exact finite ratio `31/11` rescales and retains a free positive heat
   parameter;
5. the registered 4D slab certificate leaves Lorentzian edge data and a
   Regge action open;
6. the registered finite incidence operator has no variable metric or
   selected cutoff;
7. the `alpha` and legacy `Box` gravity verifiers contain no definition of a
   tent/slab Regge volume coefficient or cutoff-to-edge map;
8. no candidate in the frozen census passes all six acceptance fields.

Only after checks 1--8, compare with the already committed tent requirement.
For a desired positive `ell(rho)`, the optimistic three-dimensional
identification can always be fitted by

```text
(Lambda a)/sqrt(s) = sqrt[ell(rho)/12].
```

Record this as an anti-selection result, not a hit.

## 5. Decision boundary

- **DERIVED CURRENT VOLUME-SELECTION ABSENCE:** the complete frozen candidate
  census has no construction satisfying all six acceptance fields, and the
  strongest spectral candidate admits a continuous family covering every
  positive `ell`.
- **SELECTED VOLUME TERM:** an already committed target-independent chain is
  found that fixes the same-carrier 4D action, sign and unique `ell`.
- **STRUCTURAL SIGN ONLY:** positivity fixes a favorable sign but the carrier
  transfer or magnitude remains free.
- **OPEN/INCOMPLETE:** the candidate census or action conventions cannot be
  reconciled.

A negative closes only the attempted spectral-volume repair from existing
data. It does not prove that no future axiom, matter sector, refinement limit
or 4D operator can generate a cosmological term.

Only the new targeted verifier and the static registry guard may run. No full
suite and no PDF build.
