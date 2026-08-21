# Protocol: local curvature mass identity on the refined static slab

Date: 2026-08-21

Prior-art/disclosure commit: `f12f56c`.

This protocol is frozen before writing or executing the registered verifier.
The exploratory match is already disclosed in the prior-art gate, so a
positive outcome is confirmatory rather than blind discovery evidence.

## 1. Frozen inputs

Require these exact files and SHA-256 values:

```text
reproducible/gravity_600cell_refined_h4_stationary_fill.json
  283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_canonical_map_feasibility.json
  ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e
reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py
  36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require the stationary artifact outcome
`REFINED_H4_INDUCED_FILL_OFF_SHELL`, `12/12` tests, 24 schedules, exactly
96 nonzero vertical entries, zero nonzero cross entries and the supplied
`tau0=0.0102`.

## 2. Independent spatial-curvature reconstruction

At 100 decimal digits reconstruct the six squared chordal rank distances from

```text
phi=(1+sqrt(5))/2,  c=phi/2,
N_k=k(1+(k-1)c),
u_a.u_b=a(1+(b-1)c)/sqrt(N_a N_b).
```

Reconstruct the chamber Gram matrix and all six spatial dihedral angles.  For
rank pair `r<s`, use the committed edge population `n_rs` and certify the
tetrahedron incidence number

```text
i_rs=14400/n_rs
```

is an integer.  Define

```text
C_rs=n_rs*s0*l_rs*(2*pi-i_rs*theta_rs),
K_r=(1/2) sum_(s != r) C_min(r,s),max(r,s).
```

Require all six deficits, all six `C_rs` and all four `K_r` to be positive,
and require `sum_r K_r` to reproduce the committed spatial Regge curvature
within `1e-68` absolute error.

This path must not read any lapse residual while constructing `K_r`.

## 3. Action-gradient reconstruction

For each committed schedule and rank, read only the certified total vertical
log residual `G_total,r`.  The frozen `P1` action assigned rank mass `M/4`, so
its analytic dust derivative is

```text
G_P1_dust,r=-4*pi*(M/4)*tau0=-pi*M*tau0.
```

Recover the gravitational derivative by

```text
G_grav,r=G_total,r+pi*M*tau0.
```

Require, independently for all `24*4=96` entries,

```text
abs(G_grav,r-tau0*K_r/2)<1e-68.                 (1)
```

Require all schedules to give one mass vector inside the same envelope.

## 4. Uniqueness and matter diagnostics

For arbitrary conserved rank masses `mu_r`, the four dust derivatives have
Jacobian

```text
-4*pi*tau0*I_4.
```

Certify symbolically that it has rank four for `tau0>0`.  Equation (1) then
selects uniquely

```text
mu_r=K_r/(8*pi).
```

Require every selected mass positive and

```text
sum_r mu_r=M
```

within `1e-68`.  Print before interpretation:

- the four `K_r/K` fractions;
- total and per-vertex rank masses;
- density ratios relative to the `P1` dual-volume masses, `4*K_r/K`;
- their maximum/minimum ratio.

The existing `P1` vector `(1,1,1,1)/4`, the original-vertices-only vector
`(1,0,0,0)` and the tempting Bernstein/binomial pattern `(1,3,3,1)/8` are
diagnostic alternatives, not targets.  Require all three to leave at least
one rank equation nonzero by more than `1e-6` in mass fraction.  In
particular, numerical proximity to `(1,3,3,1)/8` must not be promoted to an
exact identity.

## 5. Corruption and scope controls

Add `1e-10` to one recovered total lapse residual and require the identity
gate to fail.  Replace one edge population by an incompatible value and
require the integral-incidence/topology gate to fail.  These controls must be
executed, not merely described.

The verifier must record that it does not execute a root search, nested
census, Hessian, spectrum, continuum target, particle target or physical
constant extraction.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED` if provenance, geometry,
   topology, rank, corruption or scope gates fail.
2. `REFINED_LOCAL_CURVATURE_MASS_IDENTITY_REFUTED` if a valid reconstructed
   action entry violates (1).
3. `REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC` if every gate
   passes.

Outcome 3 means **DERIVED COMPUTATIONAL / STRUCTURAL, post-hoc confirmed**:
the fixed product geometry selects one curvature-matched rank-mass vector.
It does not prove homogeneous dust, refinement convergence or a physical
Hamiltonian-constraint algebra.  No outcome derives a tick, `c`, `G`, Planck
units or particle masses.

## 7. Deliverables and replication gate

Write a registered verifier and deterministic JSON artifact; run only that
verifier twice and require byte-identical output.  Then build a mechanically
different adversarial verifier which reconstructs actual spatial vertex-edge
incidences and evaluates a dust-free action derivative without recovering it
by subtracting the old `P1` term.  Do not accept the result before that gate.
Run the static registry audit only; do not run the full suite.
