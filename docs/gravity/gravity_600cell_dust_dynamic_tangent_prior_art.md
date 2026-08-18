# Prior-art gate: linearized canonical map about the first dynamic dust tick

Date: 2026-08-16

## Exact object and hypotheses

Use the first accepted non-static fixed-mass slab, not the published
time-symmetric slab.  For each independently derived order-24 schedule parity,
the old boundary has thirty equal squared lengths, the final boundary has the
accepted contracted scale, and all 35 internal orbit variables are the
accepted homothetic solution.  The complete Lorentzian Regge curvature plus
De Felice--Fabri dust action and certified angle branch are unchanged.

In the 95 logarithmic coordinates

```text
z = (o[30], x[35], n[30]),
g = (1/24) partial S/partial z,
```

compute the calibrated complete Hessian `K=partial g/partial z`.  Eliminate
the 35 internal variations and solve the linearized pre-Legendre equations to
obtain the canonical tangent map

```text
T : (delta o, delta p_pre) -> (delta n, delta p_post),
```

where `p_pre=-g_o`, `p_post=g_n`.  Apply the already derived final-to-next-old
orbit permutation to both output coordinates and output momenta before any
one-step self-map spectrum is formed.

The phase variables are canonical logarithmic pairs: changing from squared
lengths to their logarithms requires the conjugate covectors `g`, and the
symplectic form is

```text
Omega = sum_i d(delta p_i) wedge d(delta z_i).
```

The object is the **60 x 60 tangent map in the order-24 invariant quotient**.
It is not the tangent dynamics on all 720 boundary edges.

## Primary prior art

- Dittrich and Hoehn derive action-generated pre/post canonical evolution and
  background-dependent data in simplicial gravity:
  <https://arxiv.org/abs/1108.1974> and
  <https://arxiv.org/abs/0912.1817>.
- Hoehn constructs canonical linearized Regge evolution and counts
  gauge-invariant lattice gravitons around flat backgrounds:
  <https://arxiv.org/abs/1411.5672>.
- Bahr and Dittrich show that curvature generally breaks exact Regge gauge
  symmetry into pseudo-constraints: <https://arxiv.org/abs/0905.1670>.
- Barrett et al. formulate implicit Sorkin evolution and demonstrate it on a
  homogeneous dust-filled 600-cell:
  <https://arxiv.org/abs/gr-qc/9411008>.
- De Felice and Fabri evolve the dust 600-cell with five ordered classes:
  <https://arxiv.org/abs/gr-qc/0009093> and
  <https://arxiv.org/abs/gr-qc/0106077>.
- Liu and Williams study a different inhomogeneous closed lattice cosmology:
  <https://arxiv.org/abs/1502.03000>.

Thus canonical linearized Regge dynamics, discrete gravitons in the flat
regime, inhomogeneous lattice cosmology and 600-cell dust evolution are all
**KNOWN**.  No located primary source prints the present dynamic-slab
`60 x 60` quotient tangent map.  External novelty is **OPEN**; search is not
proof.

## KNOWN / CONTROL / OPEN

- **KNOWN:** a regular type-I discrete action generates a symplectic tangent
  map after internal variables are eliminated.
- **KNOWN:** compact continuous `S^3` already has discrete harmonic spectra,
  gaps and symmetry degeneracies.  Those features alone cannot be claimed as
  discrete new physics.
- **CONTROL:** the first non-static homothetic slab is committed with all 35
  internal equations and all 30 pre-momentum seam components passing.
- **CONTROL:** the published static slab's complete `65 x 65` pre-Legendre
  Jacobian is resolved full rank, but the dynamic slab's rank remains to be
  recomputed rather than assumed.
- **CONTROL:** the tangent map must satisfy calibrated symplecticity and
  reciprocal spectral pairing independently in both schedule parities.
- **OPEN:** dynamic-slab canonical rank and conditioning.
- **OPEN:** whether the homogeneous scale phase plane is invariant under the
  complete tangent map or mixes with the 58-dimensional zero-sum shape phase
  space.
- **OPEN:** schedule-parity dependence of the target-independent spectrum.
- **OPEN:** stability or amplification of shape perturbations.
- **OPEN:** identification of gauge-invariant curvature/tensor modes,
  refinement and continuum propagation.

## Framing attack

The order-24 quotient contains one scale configuration direction and 29
zero-sum shape directions.  It does not carry the full `H4` representation on
720 edges.  Therefore degeneracies inside this quotient cannot be advertised
as the complete 600-cell graviton spectrum, and absent irreducible sectors
cannot be reconstructed from it.

A one-step eigenvalue is also not yet a physical wave frequency.  The
background is changing, the matter perturbations are frozen, and no proper
time or physical norm has been derived.  The defensible first questions are:

1. does the complete action produce a resolved symplectic tangent map;
2. is scale/shape separation dynamically invariant;
3. are target-independent spectral data robust between schedule parities;
4. are any amplifications resolved beyond derivative uncertainty?

Only after committing that blind census may it be compared with continuous
`S^3` tensor harmonics.  No value of `c`, desired degeneracy, continuum
frequency or Planck scale may enter the enumeration.

Only the new targeted verifier will be run.  The full suite will not be run.
