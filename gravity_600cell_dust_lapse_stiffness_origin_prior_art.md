# Prior-art gate: origin of the scalar lapse stiffness

Date: 2026-08-17  
This note follows the observed full Schur result and precedes its source
decomposition.  Any scalarity/cancellation claim is therefore post-result and
must remain labelled **PATTERN** unless independently derived.

## Exact object

The accepted full pole Schur audit found, in every irreducible `2T` sector and
both schedule parities, a matrix numerically consistent with

```text
S_total = alpha I_120,
alpha = -4.24456181727093e-9.
```

The new question is not whether this number is nonzero; that is already
certified.  It is which part comes from the gravitational Regge action and
which part comes from the uniform dust action.

For a formal dust multiplier `mu`, keep the accepted geometry fixed and write

```text
J(mu) = J_gravity + mu J_dust.
```

The dust action is a sum of 120 independent pole terms,

```text
S_dust = -(8 pi M/120) sum_p sqrt(rho_p).
```

In logarithmic pole variables its Hessian is

```text
h_dust I_120,
h_dust = -(2 pi M/120) sqrt(rho).
```

It has no strong-strong or strong-pole entries.  Therefore, for the frozen
strong partition, Schur algebra gives the exact affine identity

```text
S(mu) = S_gravity + mu h_dust I_120.
```

This derivation does not depend on the observed value of `alpha`.  What remains
empirical is whether `S_gravity` is scalar and how close the physical
normalization `mu=1` lies to the zero of the affine family.

## Primary prior art

### KNOWN

- Lapse variation imposing gravitational constraints, and its discrete
  alteration into pseudo-constraints when Regge gauge symmetry is broken, is
  standard in [Bahr--Dittrich](https://arxiv.org/abs/0905.1670) and
  [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817).
- Lorentzian simplicial cosmology with dust shells and lapse integration is
  treated by Dittrich, Gielen and Schander in [Lorentzian quantum cosmology
  goes simplicial](https://arxiv.org/abs/2109.00875).
- De Felice and Fabri introduced matter into the evolving dust 600-cell and
  studied the symmetric dynamics and causal endpoint:
  [2000](https://arxiv.org/abs/gr-qc/0009093),
  [2001](https://arxiv.org/abs/gr-qc/0106077).

### CONTROL

- The mass `M` and the accepted non-static state were fixed before the full
  120-dimensional calculation.
- The full pole Schur operator is regular with all 120 singular directions
  resolved by at least `4.08e12` times the frozen nonzero boundary.
- Its scalar appearance was noticed only after that rank result.

### OPEN

- Whether the gravitational Schur term is structurally proportional to the
  identity, rather than merely consistent with it at this one background.
- Whether the near cancellation is a discrete Hamiltonian-constraint
  mechanism, a consequence of the imposed mass normalization, or an accidental
  feature of this carrier.
- How the residual scales with refinement, curvature, dust density or tick
  size.

## Search result

The focused pre-decomposition search found the general constraint and dust
cosmology mechanisms, but no primary source reporting this exact scalar
gravity--dust Schur cancellation for the 600-cell.  This does not prove
external novelty.

The proposed calculation can establish an exact algebraic decomposition of
the already-computed matrix and a calibrated numerical cancellation ratio.  It
cannot by itself establish a continuum Hamiltonian constraint or a physical
clock.
