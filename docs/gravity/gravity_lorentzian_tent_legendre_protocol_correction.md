# Protocol correction: causal-type sign in the complex dihedral cosine

Date: 2026-08-12

Parent protocol commit: `9b10ed4`

Status at this correction:
**POST-FAILURE CORRECTION -- PRELIMINARY FULL-RANK RESULT DISCLOSED**

## 1. What failed

The parent protocol froze

```text
C_ab=-N_ab/[sqrt(+N_aa)sqrt(+N_bb)]
```

as the complex dihedral cosine for every hinge. That unqualified formula is
wrong when the two adjacent facet normals have different causal type.

The first implementation passed the internal-angle and reality controls but
failed two load-bearing tests:

```text
full-action versus Schlaefli gradient error: 8.641,
relative Hessian antisymmetry:                0.3811.
```

The internal pole derivative still matched because timelike internal hinges
have same-type facet normals and do not encounter the missing sign.

## 2. Source-forced correction

Use equation A.38 of Borissova--Dittrich,
[arXiv:2303.07367](https://arxiv.org/abs/2303.07367):

```text
cos(theta_ab^+)
 = 16 [partial V_sigma/partial s_ab]
   /[sqrt(+V_face_a)sqrt(+V_face_b)]
```

for `d=4`, with signed squared volumes. Retain their A.39 formula for sine:

```text
sin(theta_ab^+)
 =-(4/3) sqrt(+V_h)sqrt(+V_sigma)
   /[sqrt(+V_face_a)sqrt(+V_face_b)].
```

At same-causal-type facets, A.38 agrees with the normalized-normal formula.
At opposite-causal-type facets it has the opposite sign. This is fixed by
the signed-volume branch and is not a coefficient choice.

With A.38, before this correction was committed, the exploratory rerun gave

```text
12/12 checks,
full-action/Schlaefli gradient error 2.81e-9,
Hessian antisymmetry                  6.27e-11,
mixed rank                           12,
s_min/s_max                          7.896e-3.
```

Therefore the corrected rank result is disclosed and is not blind.

## 3. Frozen independent confirmation

Before accepting the rank, extend the verifier with checks that do not reuse
the Schlaefli-gradient Hessian as their only evidence:

1. At the witness and two deterministic causal perturbations, differentiate
   all ten complex angles of every four-simplex by centered differences and
   certify the per-simplex complex Schlaefli sum

   ```text
   sum_(h in sigma) sqrt(+V_h) d theta_(sigma,h)^+ = 0
   ```

   to absolute tolerance `2e-7`.
2. Explicitly verify on all 200 simplex-hinge incidences that A.38 agrees
   with the normalized-normal cosine for same-type facet normals and is its
   negative for opposite-type normals.
3. Compute the raw entries needed for the on-shell mixed block directly from
   centered second differences of the **full complex action**, not from the
   Schlaefli gradient:

   ```text
   S_pq, S_p,rho, S_rho,q, S_rho,rho.
   ```

4. Use step `h=2e-5` and require the resulting Schur-complement mixed block
   to agree with the gradient-derived block at the same step to relative
   Frobenius error below `2e-4`.
5. Reapply the original relative rank threshold `1e-7`. Require direct-action
   rank `12` and every singular value to agree with the gradient route to
   relative `5e-3`.

Failure of any confirmation makes the result **INCONCLUSIVE**. No further
angle-sign modification is allowed under this correction; another failure
requires a new disclosed correction.

Only the targeted verifier and static registry audit may run. No full suite
and no PDF build.
