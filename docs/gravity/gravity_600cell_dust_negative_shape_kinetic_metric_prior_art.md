# Prior-art gate: kinetic-metric root certificate

Date: 2026-08-18

## Exact object, carrier, and hypotheses

The preceding target-disclosed calculation transferred the root count
`15/0/15` through the literal Euclidean coefficient balls in every one of the
`16` negative-shape recurrence cells.  Its smallest sampled
signal-to-uncertainty ratio was approximately `8.17`; hence all cells failed
the repository's separate `100 x error` convention.

This mission asks one narrower question: does the kinetic form already selected
by the same Regge action provide a canonical induced norm in which the
coefficient-ball Rouche margin is resolved at `100x`?

For each inherited negative carrier `E_-`, let

```text
B_- = E_-^* (-M_S) E_-,
```

where `M_S` is the Hermitian centered kinetic block on the shape carrier.  The
upstream blind census certifies `-M_S` positive definite before any root target
was used.  If `B_-` remains positive-resolved, define its unique positive
square root

```text
S = B_-^(1/2)
```

and the induced kinetic norm

```text
||x||_B = ||S x||_2,
||X||_B = ||S X S^(-1)||_2.
```

The polynomial is not changed, only represented in the canonically inherited
norm:

```text
Q_B(z) = S Q(z) S^(-1).
```

The Euclidean source balls must be transported rather than silently reused:

```text
||S DeltaQ(z) S^(-1)||_2
 <= kappa_2(S)
    [epsilon_Gamma |z^2-1| + epsilon_Omega].
```

No diagonal balancing, eigenvector matrix, optimized similarity, selected
angular weight or post-result metric is admissible.  A scalar multiple of
`B_-` is immaterial and is not a parameter.

## Primary literature

**KNOWN.**  Tisseur and Higham develop structured pseudospectra and stability
radii for polynomial eigenvalue problems, emphasizing that the perturbation
structure and its weights are part of the problem rather than a cosmetic
choice:
<https://doi.org/10.1137/S0895479800371451>.

**KNOWN.**  Tisseur and Meerbergen survey quadratic eigenvalue problems with
mass, damping and stiffness matrices, their linearizations, conditioning and
the need to preserve physical matrix structure:
<https://doi.org/10.1137/S0036144500381988>.

**KNOWN.**  Marsden and West derive discrete evolution from a discrete
variational principle and show why the action-selected discrete symplectic and
metric structures, rather than an arbitrary Euclidean coordinate norm, are the
natural geometric data of a variational integrator:
<https://doi.org/10.1017/S096249290100006X>.

These sources justify auditing a physically inherited norm and transporting
its perturbations.  They do not say that such a norm must improve a Rouche
margin, and they do not study this `600`-cell Regge recurrence.  External
novelty remains **OPEN**.

## Structural caveats

- The positive kinetic form is selected upstream, but the `15`-dimensional
  negative carrier is only invariant-consistent under the calibrated source
  balls.  The same autonomy caveat remains.
- A similarity leaves the exact roots unchanged.  Any change in the sufficient
  margin measures the sharpness of the norm bound, not new dynamics.
- The factor `kappa_2(S)` is mandatory.  Omitting it would make almost any
  nonnormal polynomial look artificially robust.
- The frozen polynomial is not star-palindromic; no reciprocal-pair theorem is
  restored by changing norm.
- Even a `100x` success would establish only local frozen hyperbolicity.  It
  would not establish later-slab persistence or nonlinear growth.

## KNOWN, CONTROL, OPEN, proposed difference

- **DERIVED UPSTREAM:** all `16` Euclidean literal covers transfer `15/0/15`.
- **OPEN UPSTREAM:** every Euclidean `100x` cover fails.
- **KNOWN CONTROL:** a unitary change of basis has `kappa(S)=1` and must leave
  the Euclidean certificate numerically unchanged.
- **KNOWN CONTROL:** scalar rescaling of `B_-` must leave `Q_B`, its transported
  error bound and every verdict unchanged.
- **OPEN:** positivity margin and condition number of `B_-` on the selected
  negative carriers.
- **OPEN:** whether the canonical kinetic similarity improves or worsens the
  complete continuous-contour lower bound.

The proposed difference is therefore not another root search.  It is a blind
comparison of two fully specified induced norms on the same already-disclosed
polynomial, with all coefficient errors transported before any result is read.

