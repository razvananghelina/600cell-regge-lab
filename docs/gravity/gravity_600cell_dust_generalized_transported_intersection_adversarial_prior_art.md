# Independence gate: adversarial replication of the zero intersection

Date: 2026-08-18

## Purpose and epistemic status

The residual-certified calculation has already reported
`TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL`.  This post-result audit cannot
make that discovery blind.  It can falsify an implementation error by testing
the same rank statement through a previously constructed, mechanically
different numeric route.

For orthonormal phase bases `W0,W1`, let `W1_perp` span the target orthogonal
complement.  The intersection is the kernel of the square leakage matrix

```text
L = W1_perp^H T_2 W0 : C^30 -> C^30.
```

Thus full rank of `L` is equivalent to zero transported intersection.

## Primary numerical basis

The SVD/principal-angle characterization and perturbation of subspace
intersections are standard; see Knyazev and Argentati
<https://doi.org/10.1137/S1064827500377332> and Knyazev, Jujunashvili and
Argentati <https://doi.org/10.1016/j.jfa.2010.05.018>.

The physical caveat remains the discrete-constraint analysis of Dittrich and
Hoehn: pre/post constraint surfaces, rather than an arbitrarily imposed full
cotangent bundle, select propagating data
(<https://doi.org/10.1063/1.4818895>,
<https://arxiv.org/abs/1303.4294>).

These sources define no 600-cell rank and supply no target value. External
novelty is **OPEN**.

## Independent path

The exact calculation used residual-certified 100-digit projectors, freshly
reconstructed Flint tangent balls and the 60-by-60 residual spectrum.

The audit will instead replay the already accepted adversarial phase verifier,
which independently constructs:

- float64 generalized bases through null-space reduction and explicit
  Cholesky whitening;
- the earlier committed two-step tangent archive;
- canonical phase bases without reading numeric values from the exact
  intersection artifact.

It will then compute the full singular spectrum of `L` directly.  This audit
shares the Regge action and carrier but not the exact projector, tangent source
or decisive rank matrix.

## Required controls and limitation

- `T_pos=W1 W0^H` must give rank-zero leakage and a 30-dimensional
  intersection.
- `T_neg=W1_perp W0^H` must give rank-30 leakage and zero intersection.
- reversing and rephasing both bases must preserve every singular value to the
  preregistered roundoff floor.

The audit is float64 and is not a second exact certificate. Agreement is
**STRUCTURAL INDEPENDENT CORROBORATION**. Any cell not robustly full rank sends
the consolidated result back to **OPEN** under project rule 4.

