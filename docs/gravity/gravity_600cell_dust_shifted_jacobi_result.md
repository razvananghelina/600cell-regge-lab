# Result: the shifted three-slice Jacobi operator is certified

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, BLIND.**  The second and third accepted dust slabs
define a regular, action-derived position recurrence on all `720` edge-metric
perturbations:

```text
K_-^(2) delta q_1 + K_0^(2) delta q_2
                         + K_+^(2) delta q_3 = 0.
```

All `112` boundary-twist determinant balls exclude zero, every reconstructed
principal-function identity holds entrywise, and the recurrence reproduces
the separately committed product `T_3T_2` in both implicit and solved forms.
All fourteen schedule comparisons are robust.  The preregistered outcome is

```text
SHIFTED_JACOBI_CERTIFIED.
```

The theory now has two consecutive complete Jacobi stencils on an evolving
background.  No earlier conformal/shape or negative-mode artifact entered
this construction.

## Provenance and reproduction

| stage | commit |
|---|---|
| common primary-literature and framing gate | `920ce5d` |
| certified third-slab tangent input | `d666315` |
| blind shifted-Jacobi protocol | `9ee1f3f` |
| verifier registered before first execution | `76e67dc` |

Verifier:

```text
reproducible/verify_gravity_600cell_dust_shifted_jacobi.py
```

Artifacts:

```text
reproducible/gravity_600cell_dust_shifted_jacobi.json
SHA-256 63b37b6000146d5d53dbbc01da5c9aba9a5e3373b8bc3830a404ef0f681ecf31

reproducible/gravity_600cell_dust_shifted_jacobi.npz
SHA-256 d2f507c4a2fa11c5d7a808849c199a986278516f422cf43654f6de153ab170d0
```

Two targeted runs returned `8/8` and reproduced both artifacts byte for byte.
The full suite was not run.

## What was certified

The committed tangent balls were partitioned as

```text
T_i = [A_i B_i; C_i D_i],   i=2,3.
```

Their regular twists reconstruct the quadratic principal functions, from
which the verifier forms

```text
K_-^(2) = S_2,10,
K_0^(2) = S_2,11 + S_3,00,
K_+^(2) = S_3,01,

P_2 = -(K_+^(2))^-1 K_0^(2),
Q_2 = -(K_+^(2))^-1 K_-^(2).
```

The four derivative variants, both schedule parities and all seven minimal
sectors are retained.  The `560` midpoint/radius arrays are exactly the five
matrix families above over every declared cell.

The boundary-twist midpoint diagnostics span

```text
minimum singular value  0.00180166 ... 0.00184123,
condition number        13.0887    ... 67.0916.
```

Thus the position recurrence is not produced by a marginal twist inversion.

The largest midpoint/radius Frobenius ledgers are

```text
principal identities   8.83e-11 midpoint, 2.94e-4 radius,
product identities     5.78e-10 midpoint, 2.06e-4 radius.
```

All residual entries contain zero in the propagated Flint balls.  The
enclosures are conservative rather than marginal.

## Blind time-asymmetry diagnostic

The preregistered coefficient diagnostic is

```text
||K_+^(2) - (K_-^(2))^*||_F
--------------------------------
max(||K_+^(2)||_F,||K_-^(2)||_F).
```

It lies between

```text
0.0338021 and 0.0363930
```

in every sector and schedule.  This is a target-free output: the two later
slabs are not time-translation invariant, but their shifted stencil is much
closer to adjoint symmetry than an order-one mismatch.

**POST-RESULT COMPARISON:** the earlier stencil reported `0.22098...0.23738`.
The decrease was not an acceptance target and no physical damping conclusion
is drawn from two points.  It may reflect background evolution, coordinate
normalization or both.

## Hostile interpretation audit

1. Every regular pair of variational maps admits such a recurrence.  Its
   existence is necessary structure, not proof of Einstein gravity.
2. The coefficient matrices act on raw edge-metric perturbations before an
   exact constraint quotient or matter extension.
3. A smaller adjoint-asymmetry diagnostic is not a continuum limit and is not
   evidence of isotropization without later ticks or refinement.
4. No spectrum of centered `M,N,V` was read and no old subspace was loaded.
   Therefore persistence of the previous `30` negative modes remains open.
5. The recurrence supplies no physical unit of time, frequency or limiting
   speed.

## Status ledger

- **DERIVED COMPUTATIONAL:** two consecutive complete Jacobi stencils now
  exist on the same literal `720`-position carrier.
- **DERIVED COMPUTATIONAL:** the shifted stencil is regular, variationally
  consistent and schedule robust.
- **STRUCTURAL:** its blind time-asymmetry diagnostic is `0.034...0.036`.
- **POST-RESULT PATTERN:** this is substantially smaller than in the first
  stencil; two time positions do not establish a trend.
- **OPEN:** persistence of the action-selected conformal/shape split and the
  negative-stiffness subsystem.
- **OPEN:** physical quotient, tensor polarizations, finite-time growth,
  refinement, dispersion, causal speed and external novelty.

## Next load-bearing gate

The operator is now safely committed before its target comparison.  The next
verifier may load both Jacobi artifacts and must, with a frozen outcome tree:

1. construct the shifted centered `M_2,N_2,V_2`;
2. reconstruct the vertex-conformal and kinetic-orthogonal shape carriers
   independently at slice 2;
3. classify the shifted shape stiffness before consulting the old negative
   sectors;
4. only then compare negative-sector identities and subspace transport;
5. if the carrier persists, form the genuine non-autonomous product of the
   two consecutive restricted companions.

Eigenvalues of that two-update product will remain finite-time diagnostics,
not asymptotic Lyapunov exponents.

