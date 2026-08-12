# Preregistered protocol: transverse conformal Hessian of round `S^3`

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- NO TRANSVERSE-HESSIAN RESULT**

This protocol is written after the homogeneous Hopf Hessian and the
round--Regge path results are known. It therefore cannot preregister those
inputs. It freezes a new hostile question before its exact conformal
calculation is encoded or run:

> Does the favorable ordinary de Rham `A2` Hessian persist from homogeneous
> trace-free directions to non-gauge conformal directions of the smooth round
> three-sphere?

No measured physical target, fitted coefficient or finite-cutoff spectrum is
used.

## 1. Complete scope and hypotheses

1. The carrier for this audit is the **smooth unit round three-sphere**, not
   the finite fixed 600-cell edge-length space.
2. `D=d+d*` is the ordinary, ungraded de Rham operator on the complete exterior
   algebra of a closed smooth three-manifold.
3. Its integrated coefficient is the already derived convention

   ```text
   A2(g)=-(2/3) integral R_g dVol_g.
   ```

4. Scale is removed by the normalized functional

   ```text
   A2_hat(g)=Vol(g)^(-1/3) A2(g).
   ```

   On an exactly equal-volume path this differs from `A2` only by a positive
   constant, so it has the same Hessian sign.
5. Conformal metrics are written

   ```text
   g_epsilon=(1+epsilon*f)^4 g_round,
   ```

   for sufficiently small real `epsilon`, with `integral f=0`.
6. The positive Laplacian convention is

   ```text
   integral |grad f|^2=lambda integral f^2.
   ```

   Scalar harmonics on unit `S^3` have `lambda_l=l(l+2)`.
7. Scale (`l=0`) and diffeomorphism-generated conformal (`l=1`) directions
   are excluded from the physical transverse test. The first tested mode is
   the explicit `l=2` harmonic

   ```text
   f=x1^2-x2^2,  lambda=8.
   ```

8. The previously certified homogeneous trace-free result

   ```text
   delta^2 A2_hat(H,H)>0 for H nonzero
   ```

   is an input used only to decide whether opposite signs make the full smooth
   Hessian indefinite.

This audit does **not** claim that every smooth conformal direction is an
edge-length direction of the fixed 600-cell. Extending the conclusion to a
particular discrete Regge configuration space would require a separate
discretization/refinement argument.

## 2. Frozen exact calculation

For `u=1+epsilon*f`, the three-dimensional conformal scalar-curvature formula
gives the normalized Einstein--Hilbert/Yamabe functional

```text
Y(u)=
  [8 integral |grad u|^2 + 6 integral u^2]
  / [integral u^6]^(1/3).
```

Write `V=Vol(S^3)` and `F=integral f^2`. For a mean-zero eigenfunction the
calculation must expand numerator and denominator independently through
`epsilon^2` and test the exact identity

```text
Y(1+epsilon*f)
=6 V^(2/3)
 + 8 (lambda-3) V^(-1/3) F epsilon^2
 + O(epsilon^3).
```

Therefore the frozen predicted coefficient and Hessian are

```text
[epsilon^2] A2_hat
  =-(16/3)(lambda-3)V^(-1/3)F,

delta^2 A2_hat(f,f)
  =-(32/3)(lambda-3)V^(-1/3)F.
```

The verifier must not merely assume the displayed result. It must derive it
symbolically from the two expansions and independently check the explicit
`l=2` harmonic data on `S^3`:

```text
V=2*pi^2,
integral f=0,
integral f^2=V/6,
integral |grad f|^2=8 integral f^2.
```

## 3. Hostile controls

The verifier must distinguish all three sectors:

- `lambda=0`: the displayed mean-zero formula does not apply; this is scale.
- `lambda=3`: the conformal Hessian vanishes; this is the Möbius/diffeomorphism
  zero sector and cannot establish stability.
- `lambda=8`: the first non-gauge scalar harmonic; its sign decides the hostile
  test.

It must also check that the sign reversal comes from the already fixed
ordinary de Rham multiplier `-2/3`. Replacing that multiplier by `+1` is a
formal control and must reverse the conclusion; it is not an admissible rescue
of the present theory.

## 4. Preregistered decision boundary

- **DERIVED SMOOTH SADDLE:** the exact `l=2` conformal Hessian is negative,
  the prior nonzero homogeneous trace-free Hessian is positive, and `l=2` is
  outside scale/diffeomorphism zero modes.
- **REFUTED HOSTILE PREDICTION:** the exact conformal Hessian is nonnegative in
  the `l=2` direction under the frozen conventions.
- **OPEN:** the conformal expansion, harmonic data, or compatibility with the
  ordinary de Rham coefficient cannot be established exactly.

If the saddle boundary is met, it refutes only the promotion of the
homogeneous/one-path `A2` selector to a minimum over all smooth metrics. It
does not refute the certified Hopf theorem or the certified affine Regge path.
It also means that a positive-cutoff `A2` term alone cannot be the complete
Euclidean gravitational vacuum selector; higher spectral terms, a contour or
constraint, or different dynamics would become load-bearing and remain
**OPEN**.

## 5. Reproducibility boundary

The eventual verifier must be registered in `reproducible/run_all.py`, but
only that targeted verifier and a non-executing static registry check will be
run in this mission. No full-suite run and no PDF build are authorized.
