# The honest global Regge reduction has 35 variables and remembers phase parity

Date: 2026-08-12

Upstream five-phase schedule commit: `d439b07`

Orbit-reduction protocol commit: `06a1c6a`

Disclosed correction commits: `08b638c`, `4a28c41`, `074a82b`

Registered verifier:
`reproducible/verify_gravity_global_regge_orbits.py`

Machine-readable result:
`reproducible/gravity_global_regge_orbits.json`

Final targeted run: **43/43 passed**.  No full suite was run.

## Headline

> **DERIVED COMPUTATIONAL:** after fixing an ordered five-phase slab and its
> two regular 600-cell boundaries, the exact schedule stabilizer reduces the
> 840 internal Regge edge variables to 35 honest symmetry orbits: 30
> staircase-diagonal orbits and five tent-pole orbits.  A 100-simplex-orbit
> action reproduces the unreduced 2400-simplex action and all derivatives at
> three Lorentzian controls.

The regular metric is again nonstationary.  More unexpectedly, the two
`H4`-inequivalent phase-order parities have different Hessian spectra at this
same regular metric.  Thus the current discrete Regge action remembers the
phase order at quadratic level.

No stationary root was searched in this audit.

## 1. Exact orbit reduction

For either phase parity, the pointwise stabilizer of all five phase cells has
order 24.  It acts freely on every audited internal simplex layer:

```text
internal edges       840 =  35 x 24
  diagonals          720 =  30 x 24
  poles              120 =   5 x 24
internal triangles  3840 = 160 x 24
internal tetrahedra 5400 = 225 x 24
four-simplices      2400 = 100 x 24.
```

For each of the ten unordered pairs of phases there are three distinct
diagonal orbits.  Consequently a model with just one diagonal length per
phase pair would have 15 variables,

```text
10 diagonal classes + 5 poles,
```

but would impose twenty equalities not supplied by an automorphism.  Such a
model would be a fitted/truncated ansatz.  The smallest invariant system is
35-dimensional.

## 2. Full action versus reduced action

The verifier uses the corrected Borissova--Dittrich plus complex-angle branch
already certified locally.  It compares the 100-orbit evaluator with a fresh
sum over all 2400 four-simplices at three deterministic controls:

```text
R0  all diagonals 1; all pole rho=1/4
R1  small positive orbit-dependent perturbations
R2  small negative diagonal and positive pole perturbations.
```

Across both phase parities and all controls, the worst errors are:

```text
full/reduced action relative error       1.11e-14
full/reduced restricted-gradient error   1.42e-14
triangle-orbit curvature error           2.82e-15
direct full-action derivative error      2.37e-6
per-simplex Schlaefli residual            3.81e-10.
```

Every control has 100/100 reduced representatives and 2400/2400 unreduced
simplexes with exactly one negative Gram eigenvalue.  The largest real-branch
residual is below `7e-13`, and the smallest angle-argument modulus is `0.8567`.

Therefore the reduction is not manufacturing a smaller action: it is the
exact restriction of the complete action to the group-fixed subspace.

## 3. The regular global metric is not a vacuum tick

At `rho=1/4`, the local regular pole deficit is

```text
epsilon = 0.196545606092 radians.
```

Every one of the five pole-orbit gradients is

```text
30.8897138013,
```

up to roundoff, for both phase parities.  Hence all 120 individual pole
equations are nonzero.  Global assembly does not cure the local symmetric
Lorentzian no-go.

This is **DERIVED NEGATIVE** for regular boundary/interior data.  It is not a
no-go for asymmetric internal lengths.

## 4. Linearization is nonsingular

At regular data, centered differentiation of the 35 restricted gradients
gives a real symmetric Hessian for each parity.  At relative thresholds
`1e-7`, `1e-9`, and `1e-11`, both ranks are stably

```text
rank = 35.
```

Both signed Hessians have inertia

```text
29 positive, 6 negative, 0 zero.
```

Their conditioning controls are:

```text
even s_min/s_max = 0.0180807
odd  s_min/s_max = 0.0156534.
```

Thus an invariant Newton route is locally well posed near the regular point.
Since the point is not stationary, this Hessian is not called a physical
propagator or stability operator.

## 5. Phase parity survives at quadratic order

The first completed run was `38/39`: the sole failing assertion incorrectly
expected the two parity Hessian spectra to agree.  Their observed relative
separation was `0.0162335`.  This was recorded before changing the assertion.

A first direct-trace control at step `2e-5` then gave `39/41`: the real traces
agreed, but cancellation amplified action imaginary residuals to order
`1e-2`.  That failure was also recorded.  The final preregistered control uses
the originally proposed step `5e-4` and obtains:

```text
even direct trace   746.777097589
even Hessian trace  746.776929995

odd direct trace    746.777090768
odd Hessian trace   746.776929994

relative trace errors  2.27e-7 and 2.17e-7.
```

Changing the Hessian step from `2e-5` to `5e-4` changes each spectrum by only
about `1.1e-6`.  The recomputed parity separation remains

```text
relative maximum singular-spectrum difference = 0.0162326.
```

Some coarse invariants coincide:

```text
trace             746.776929995  (both)
Frobenius norm    208.042818224  (both)
signed inertia    (29+,6-)       (both).
```

But the spectra and log absolute determinants do not:

```text
even log|det H| = 108.722759763
odd  log|det H| = 108.653292930.
```

Therefore the difference is subtler than an overall scale or a changed mode
count, but it is not numerical noise at the frozen precision.

- **DERIVED COMPUTATIONAL:** the two ordered phase parities have different
  quadratic Regge actions.
- **DERIVED NEGATIVE:** present discrete Regge dynamics does not erase the
  phase-order choice at this finite triangulation.
- **OPEN:** whether refinement restores refoliation invariance.
- **OPEN:** whether either parity is selected dynamically.
- **NOT CLAIMED:** fermionic chirality, CP violation or an arrow of time.

## 6. What this does and does not establish

| Claim | Status |
|---|---|
| Exact 35-variable invariant reduction | **DERIVED COMPUTATIONAL** |
| Reduced evaluator equals complete action/gradient | **DERIVED COMPUTATIONAL** |
| Regular global metric is stationary | **DERIVED NEGATIVE** |
| Regular reduced Hessian has rank 35 | **DERIVED COMPUTATIONAL** |
| Even and odd phase orders are dynamically equivalent at quadratic order | **DERIVED NEGATIVE** |
| One parity is physically preferred | **OPEN** |
| A stationary point exists in either 35-variable sector | **NOT SEARCHED** |
| A stationary point exists in the full 840-variable system | **OPEN** |
| Five phases define a duration or speed | **OPEN** |

## 7. Next decisive calculation

Preregister a root search separately for both phase parities.  The search must:

1. use all 35 variables, with no extra orbit equalities;
2. start from regular data and from deterministic causal perturbations fixed
   before outcomes are known;
3. use analytic Schlaefli gradients and their verified Hessian;
4. reject every iterate that loses Lorentzian simplex inertia or crosses the
   angle branch cut;
5. require all 35 restricted derivatives and all 840 reconstructed individual
   edge derivatives to vanish;
6. validate any root using the unreduced 2400-simplex evaluator and an
   independent higher-precision calculation;
7. label failure only as a no-root result in the chosen invariant sector, not
   as a full 840-variable no-go.

If roots exist in only one parity, that would supply a dynamical selector.  If
both exist, their on-shell actions and boundary Legendre data must decide
whether they are physically distinct.  If neither exists, the five-phase slab
remains a valid Lorentzian carrier but not a vacuum evolution.
