# Vector sampling/intertwining theorem and normalization no-go

Date: 2026-07-22

## Choice of discrete target

The target is the 30-dimensional space of oriented edge cochains on the
icosahedral Hopf base.

- **DERIVED:** geodesic edge integration is the discrete de Rham map
  `I(alpha)_ij = integral_[i,j] alpha` and is canonical once the oriented
  spherical icosahedron is fixed.
- **DERIVED:** it commutes with the exterior derivative on scalar functions:
  `I(dY)=d_0(EY)`.
- Vertex-tangent samples are not used because converting their 24 components
  into a 1-cochain requires an additional endpoint/midpoint projection and
  weighting convention. Thus they do not address the normalization question
  without already adding extra structure.

## Band-limited vector theorem

Let `V_l^E={dY:Y in H_l}` and `V_l^C={*dY:Y in H_l}` for `l=1,2` on the unit
round sphere. Let `S` be geodesic edge integration into icosahedral edge
cochains. Then:

1. **DERIVED:** `S` is injective on
   `V_1^E + V_2^E + V_1^C + V_2^C`, of dimension
   `3+5+3+5=16`.
2. **DERIVED:** the first two summands are exact discrete cochains and the
   latter two are coexact; the exact and coexact sampled spaces are orthogonal
   for the unweighted edge inner product.
3. **DERIVED:** for the combinatorial edge Hodge Laplacian
   `Delta_1=d_0 d_0^T+d_1^T d_1`, the four sampled blocks have eigenvalues

   - exact `l=1`: `5-sqrt(5)`;
   - exact `l=2`: `6`;
   - coexact `l=1`: `3-sqrt(5)`;
   - coexact `l=2`: `2`.

4. Define `P(x)` by

   `P(x)=(42-24sqrt(5))+(-24+22sqrt(5))x`
   `     +(3-6sqrt(5))x^2+(sqrt(5)/2)x^3`.

   **DERIVED:** on the 16-dimensional sampled space,

   `P(Delta_1) S = S Delta_Hodge`,

   where the continuum Hodge eigenvalue is `l(l+1)`, hence `2` on `l=1` and
   `6` on `l=2` for both Hodge types.
5. **DERIVED:** the Moore--Penrose reconstruction is a left inverse on these
   16 modes. The 14-dimensional alias complement splits into 3 exact and 11
   coexact dimensions.

The verifier is `reproducible/verify_vector_sampling_intertwiner.py`.

## Inner-product normalization no-go

For the unweighted edge norm, each irreducible block is conformal to its
continuum `L^2` norm, but with a different factor:

- exact `l=1`: `(4+phi^{-4})/pi`;
- exact `l=2`: `3/pi`;
- coexact `l=1`: `15 theta^2/(4pi)`, where
  `theta=arccos(1/sqrt(5))` is the spherical edge length;
- coexact `l=2`: `(4+phi^{-4})/pi`.

These are numerically approximately
`(1.31968033, 0.95492966, 1.46316505, 1.31968033)`.

- **DERIVED:** the factors are unequal, so one overall edge weight cannot make
  all four blocks isometric.
- **DERIVED:** all 30 icosahedral edges form one symmetry orbit. Therefore any
  local diagonal `A_5`-invariant edge metric has only one weight and cannot
  repair the mismatch.
- **STRUCTURAL:** one may force an isometry by assigning separate spectral
  weights to the four blocks. That metric is mode-dependent/nonlocal and its
  weights are additional choices; its extension to the 14-dimensional alias
  complement remains arbitrary.
- **OPEN:** a geometric principle selecting those spectral weights and their
  alias-sector extension.

Consequently this theorem does **not** supply the isometry axiom required in
`alpha_phi4_missing_axiom.md`. It does not canonically normalize the physical
`U(1)` field, so it neither derives nor contradicts
`1/alpha_0=20 phi^4`. That identification remains **STRUCTURAL**.

## Consequence for the gauge prefactors

The 12-dimensional gauge kernel is a vertical-amplitude space on the 12 Hopf
fibers, while the vector sampling target is a base-edge space. The canonical
map from amplitudes to base edges is `d_0`, which has rank 11:

- **DERIVED:** it kills the constant `1` sector;
- **DERIVED:** its low exact image contains only the `3+5` sectors;
- **DERIVED:** the remaining `3'` is the three-dimensional exact alias sector.

Thus the vector bridge is not an injective trace-preserving map of the full
`1+3+3'+5` gauge kernel and supplies no trace on `1+3+8`. In particular it
cannot determine `(8/15,1/3,2/15)`.

- Gauge prefactors: **PATTERN**.
- A gauge-kernel-to-vector-field map preserving the constant mode, together
  with a derived bracket on `3'+5` and a common matter trace: **OPEN**.
