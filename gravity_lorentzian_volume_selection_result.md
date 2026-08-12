# The current theory does not select the Lorentzian tent volume coefficient

Date: 2026-08-12

Preregistered protocol commit: `dc7a6ab`

Registered verifier:
`reproducible/verify_gravity_lorentzian_volume_selection.py`

Machine-readable result:
`reproducible/gravity_lorentzian_volume_selection.json`

## Headline

> **DERIVED CURRENT VOLUME-SELECTION ABSENCE.** None of the six existing
> certified candidate classes supplies an action/transfer, sign and unique
> dimensionless coefficient `ell=lambda a^2` on the same Lorentzian tent
> carrier. The strongest spectral candidate has a continuous cutoff freedom
> that realizes every positive `ell`; matching the tent equation would be
> fitting.

> **STRUCTURAL SIGN ONLY.** The ordinary positive three-dimensional de Rham
> spectral action can be normalized to the form
>
> ```text
> integral R_3 - lambda_3 Vol_3,
> lambda_3=12 Lambda^2 C0/C2>0.
> ```
>
> Thus its relative sign is compatible with the positive coefficient needed
> by the static Lorentzian tent. It does not select its magnitude and is not a
> derived four-dimensional Lorentzian transfer.

The targeted verifier passes `24/24`. No full suite was run.

## 1. What selection required

A coefficient counted as selected only if one already committed construction
supplied all six fields:

1. an action on the same four-dimensional Lorentzian carrier, or a proved
   transfer to it;
2. curvature and four-volume in one action;
3. their relative sign and normalization;
4. one unique `ell=lambda a^2`;
5. a selected cutoff function and cutoff-to-edge scale when spectral input is
   used;
6. a native 4D construction or a derived `3+1` transfer including the missing
   temporal structure.

No desired pole time, Euclidean golden ratio, measured cosmological constant
or Planck scale was admissible input.

The candidate census was frozen before the final comparison:

| candidate class | fields passed | decisive gap |
|---|---:|---|
| generic Lorentzian tent extension | `4/6` | sign/value of `lambda` not selected |
| spatial de Rham spectral action | `0/6` | 3D carrier, free cutoff, no 3+1 transfer |
| fixed finite incidence moments | `0/6` | fixed matrix, no variable metric or cutoff |
| canonical CW slab | `2/6` | carrier exists, action does not |
| `Box` edge stiffness | `1/6` | static stiffness, no volume/action dictionary |
| `Box` electromagnetic-alpha equation | `1/6` | gauge normalization, no gravitational bridge |

Fully passing candidates: `0/6`.

The counts are an audit aid, not a claim that all future action types have
been enumerated. The result is explicitly about the present certified
repository.

## 2. Why the finite moments do not supply `lambda`

The fixed incidence Kähler--Dirac operator has the exact invariants

```text
(c0,c1,c2)=(2640,14880,55920),
r=c1/(2c0)=31/11.
```

These are Taylor moments of a finite heat trace:

```text
Tr exp(-tD^2)=c0-c1 t+c2 t^2+... .
```

They are not continuum Seeley--DeWitt coefficients and therefore cannot be
renamed as cosmological, Einstein--Hilbert and gauge coefficients. This was
already certified by the dimension reconciliation.

The ratio `r` is useful as an internal inverse-length-squared unit. Under
`g -> c^2 g`, it scales as `r -> r/c^2`. A scale-covariant heat time is

```text
t=beta/r,
```

but every `beta>0` is equally covariant. Hence `31/11` does not choose a heat
time or a physical cutoff. Setting `beta=1` because it looks simple would be
a normalization choice, not a derivation.

For a Gaussian finite spectral action, every positive eigenmode obeys

```text
d/dLambda exp(-lambda_i/Lambda^2)>0.
```

The exact finite heat trace therefore cannot select an interior cutoff by
extremizing itself either.

## 3. The strongest optimistic spectral calculation

The certified smooth ordinary de Rham coefficients on the spatial
three-sphere are

```text
A0=8 Vol_3,
A2=-(2/3) integral R_3 dVol_3.
```

For a positive even cutoff,

```text
S_chi ~ (4*pi)^(-3/2)
        [Lambda^3 C0 A0+Lambda C2 A2+...],
C0>0, C2>0.
```

Multiplying the whole classical functional by
`-3/(2 Lambda C2)` normalizes its curvature coefficient to `+1` and gives

```text
S_normalized
 = integral R_3 dVol_3
   -12 Lambda^2(C0/C2) Vol_3+... .
```

Thus

```text
lambda_3=12 Lambda^2 C0/C2>0.
```

Only the relative sign is robust. An overall action sign is irrelevant to
the isolated classical stationarity equation but would matter once matter
and path-integral weighting are included; no such complete coupling is
claimed here.

## 4. Exact anti-selection family

The admissible positive Gaussian cutoffs

```text
chi_s(v)=exp(-s v^2),       s>0,
```

have exact Mellin moments

```text
C0=s^(-3/2),       C2=s^(-1/2).
```

Consequently

```text
lambda_3 a^2=12 (Lambda a)^2/s
             =12 mu^2,
mu=(Lambda a)/sqrt(s)>0.
```

This map is onto all positive real numbers. Given any desired coefficient,

```text
mu=sqrt(ell/12)
```

reproduces it exactly. The cutoff-shape parameter and scale reduce to one
effective continuous freedom in this Gaussian family, but one freedom is
already enough to destroy selection.

Even importing the finite unit `r=31/11` leaves

```text
ell=12r/beta=372/(11 beta),       beta>0,
```

which again ranges from zero to infinity. The internal spectral number does
not close the missing dimensionless parameter.

## 5. Comparison with the tent only after the census

The already committed static Lorentzian tent controls require

| pole regime | required `ell` | fitted `mu=sqrt(ell/12)` |
|---|---:|---:|
| `rho -> 0+` | `1.307292211049...` | `0.330062142009...` |
| `rho=1` | `3.214668192665...` | `0.517579961026...` |
| `rho -> infinity` | `8.706236948324...` | `0.851774468797...` |

All three are reproduced merely by choosing the free effective cutoff. So a
numerical hit at any pole time would carry exactly zero evidence for that
time. The map would have been able to hit every positive target.

This is the same methodological failure seen earlier in the orbifold
induction identities: a flexible family that contains the answer does not
select the answer.

## 6. The 3D-to-4D category boundary

There is an even earlier obstruction. The displayed `A0,A2` live on a smooth
Riemannian `S3`. The tent action is a four-dimensional Lorentzian Regge
functional whose curvature resides on triangular spacetime hinges.

The repository currently has:

- a genuine spatial de Rham curvature response;
- a canonical four-dimensional cellular cylinder `K x I`;
- a local Lorentzian tent geometry;

but it has no single operator/action combining all three. In particular it
has no derived lapse/extrinsic-curvature term or proof that the spatial
`A0/A2` ratio transfers to the spacetime Regge ratio.

Therefore the formula for `lambda_3` is an intentionally optimistic hostile
control. It proves non-selection even after granting a transfer that has not
itself been established.

## 7. Historical-source audit

Two old exploratory files could otherwise create confusion:

- `exp261_spectral_action.py` calls the finite spectral action the full
  bosonic Lagrangian. It is unregistered and is superseded by the finite-
  moment/Seeley audit.
- `exp513_cc_spectral_action.py` tunes a heat time against `10^-122`. It is
  also unregistered and explicitly concludes that the cutoff is an input,
  not an output.

Their numerical calculations are not silently discarded, but neither is a
certificate for a tent volume coefficient. The registered electromagnetic
`alpha` verifier and legacy `Box` gravity verifier likewise contain no
tent/slab volume or cutoff-to-edge map.

## 8. Status ledger

| Claim | Status |
|---|---|
| Finite moment triple and `31/11` are exact | **DERIVED** |
| Those moments are continuum cosmological/EH coefficients | **REJECTED** |
| `31/11` fixes the dimensionless heat time | **REFUTED** |
| Spatial de Rham `A0/A2` gives a positive normalized `lambda_3` | **DERIVED CONTINUUM / STRUCTURAL FOR TENT** |
| Positive cutoff fixes the magnitude | **REFUTED** |
| Gaussian cutoff family covers every positive `ell` | **DERIVED** |
| Existing 3D action transfers to the 4D Lorentzian tent | **OPEN / ABSENT** |
| Canonical 4D slab already has a spectral/Regge action | **REFUTED AS CURRENTLY BUILT** |
| Existing theory selects the tent volume coefficient | **DERIVED CURRENT ABSENCE** |
| Matching a desired tent root by choosing cutoff is evidence | **REFUTED; FITTING** |
| No future construction can generate a cosmological term | **NOT CLAIMED** |
| `G`, Planck units or a physical cosmological constant follows | **OPEN** |

## 9. Consequence and next gate

The volume-term rescue cannot be used with current data. There are now two
honest options:

1. construct a genuinely four-dimensional Lorentzian operator/action that
   selects its own scale and coefficient; or
2. keep `lambda=0` and test whether target-independently selected
   nonsymmetric tent data or a different causal carrier can satisfy the
   Regge equation.

The second route was executed under preregistration in
`gravity_lorentzian_asymmetric_tent_result.md`. The per-simplex angle bound
does **not** extend to arbitrary final edge lengths: a strictly admissible
asymmetric zero-volume pole exists. Its boundary data are target-found and
not selected, so this repairs existence but not the physical-clock gate.

## 10. Reproduction history

The preliminary non-selection formula was disclosed in protocol commit
`dc7a6ab`. The first targeted implementation passed every mathematical and
candidate-census check but reported `21/23` because two Markdown source guards
required exact phrases across line breaks. They were replaced by conjunctions
of independent semantic fragments; no formula, candidate field or verdict
changed. A further hostile check proved that the Gaussian finite spectral
action is monotone in its cutoff and cannot select an interior value by
extremization. The final targeted verifier passes `24/24`.

No full suite and no PDF build were run.
