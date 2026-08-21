# Protocol: homogeneous coarse-to-refined cotangent lift

Date: 2026-08-21

Prior-art gate commit: `3e188e7`.

This protocol is frozen before constructing either pullback matrix or its
kernel.

## 1. Frozen inputs

Use and hash exactly:

```text
reproducible/gravity_600cell_refined_canonical_map_feasibility.json
  ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e
docs/gravity/gravity_600cell_refined_canonical_map_feasibility_result.md
  222739a82680c35c337b127c455fc0a8a7c24c2bd4a2f6c8f9953eeb3251e681
docs/gravity/gravity_600cell_refined_homogeneous_cotangent_lift_prior_art.md
  fa391b36501f08dfee9a0bd588f651e772cebc38643acd4c852fc95bb8cd6f21
```

Require the feasibility artifact's accepted outcome, `K0` f-vector, six
colour-pair populations and their sum before evaluating the target rank.

## 2. Exact canonical pullbacks

Order the six rank-pair types lexicographically:

```text
(01,02,03,12,13,23).
```

For logarithmic squared-edge coordinates, homothetic scale has tangent

```text
h=(2,2,2,2,2,2).
```

If `P_i` is the total canonical momentum of orbit `i`, construct over the
rationals

```text
R_total = [2,2,2,2,2,2],       p_s = R_total P.
```

If `mu_i` is a common per-edge momentum and `N_i` is its edge population,
construct

```text
D=diag(N_i),
R_edge=R_total D,               p_s = R_edge mu.
```

Compute exact ranks and nullspaces without a numerical tolerance.  Verify
`D` is invertible and that multiplication by `D` maps the per-edge equation
to the orbit-total equation.

## 3. Controls

1. Reproduce exactly

   ```text
   N=(1440,3600,2400,3600,3600,2400), sum N=17040.
   ```

2. Reconstruct the six positive unit squared edge lengths from the exact
   golden-ratio rank formula used by the refined geometry.
3. In a synthetic one-orbit refinement, require rank one, nullity zero and
   the unique lift `P=p_s/2`.
4. Replace logarithmic squared length by logarithmic length, hence replace
   `h=2*1` by `h=1`; rank and nullity must be unchanged.
5. Reverse the six orbit labels and require unchanged rank and nullity.
6. Set one synthetic population to zero and require the convention map `D`
   to become singular.  This negative control ensures that convention
   equivalence is not asserted without positive orbit populations.

## 4. Frozen outputs and interpretation

Before interpretation, write:

- both exact pullback rows;
- rank and nullity in both conventions;
- an exact basis for each kernel;
- the number of free parameters in the affine solution for fixed nonzero
  coarse momentum;
- all control results.

No preferred particular solution is allowed.  In particular, the verifier
must not compute a pseudoinverse, minimum-norm lift, population-weighted lift,
edge-length-weighted lift or action-metric lift.

Use the first applicable outcome:

1. `REFINED_HOMOGENEOUS_COTANGENT_CONTROL_FAILED` if an input or control
   fails;
2. `REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNIQUE` if both actual pullbacks have
   nullity zero;
3. `REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNDERDETERMINED` if both have rank
   one, nullity five, convention equivalence holds and the one-orbit control
   is unique;
4. `REFINED_HOMOGENEOUS_COTANGENT_LIFT_OPEN` otherwise.

Outcome 3 is **STRUCTURAL**, not a no-go for an action-selected or
perfect-action transport.  It forbids only the claim that geometry, `H4` and
the symplectic pairing already select the refined momentum.

## 5. Reproduction boundary

Create one registered verifier and deterministic JSON, run only that verifier
twice, and run the static registry/coverage guard.  Do not run the full suite,
the nested `H4` root census, a refined slab solve, a Hessian or a spectrum.
