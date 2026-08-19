# Prior-art and framing gate: classify the 240 canonical-data directions

Date: 2026-08-19

## Frozen observation before interpretation

The target-blind admissibility census is frozen at commit `e3a77fe`.  At both
predeclared nonstatic representatives and both finite-field primes it found

```text
rank(C_z)=3600,
rank(C_aug)=4200,
modular nullity(C_aug)=4440-4200=240.
```

Two rational homothetic directions were constructed, so the rational
admissible-data dimension is known only to lie in

```text
2 <= dim_Q A <= 240.
```

The equality `240=2*120` was seen only after the target-blind artifact was
committed.  It is therefore a **PATTERN**, not yet an identification.

## Target-disclosed candidate

For every global 600-cell vertex `v`, introduce two rational scalar values

```text
sigma_v  local radial/scale displacement,
nu_v     local normal/lapse displacement.
```

On a tetrahedron with canonical local vertex `p_i`, test the upper-vertex
displacement

```text
delta q_i=sigma_v p_i+nu_v n.
```

At the regular homothetic background its natural data would be

```text
delta l_uv^2 = 8 lambda (sigma_u+sigma_v),
delta s_v^2  = 6(lambda-1)sigma_v-2 tau nu_v.
```

These formulas must be derived again from the local length Jacobian.  Their
dimension match with the modular nullity is not evidence by itself.

## What the literature establishes

Vertex-based conformal edge variations and vertex-displacement variables are
standard discrete constructions, while connection/shape-matching data are
constrained before they reduce to length Regge geometry:

- Luo, *Combinatorial Yamabe flow on surfaces*,
  DOI [`10.4310/CAG.2004.v12.n3.a2`](https://doi.org/10.4310/CAG.2004.v12.n3.a2),
  and the later discrete-conformal literature use vertex scalars whose edge
  variations involve endpoint sums.  This establishes the broad pattern,
  not the present Lorentzian frustum carrier.
- Hoehn, *Canonical linearized Regge Calculus*,
  arXiv:[1411.5672](https://arxiv.org/abs/1411.5672), separates vertex
  lapse/shift variables from curvature-carrying lattice gravitons on flat
  backgrounds.
- Bahr and Dittrich, *Broken Gauge Symmetries and Constraints in Regge
  Calculus*, arXiv:[0905.1670](https://arxiv.org/abs/0905.1670), explain why
  curvature can lift those exact gauge directions into pseudo-constraints.
- Dittrich and Ryan, arXiv:[0807.2806](https://arxiv.org/abs/0807.2806), and
  Anza and Speziale, arXiv:[1409.0836](https://arxiv.org/abs/1409.0836),
  delimit the larger connection/area-angle space and its shape-matching
  reduction.

No located source computes the present 600-cell `4200/4440` augmented rank
or identifies its rational kernel.  Search absence does not prove novelty;
external novelty remains **OPEN**.

## Exact new question

Does the complete rational kernel of the already frozen augmented
compatibility matrix equal the image of the 240-column vertex carrier

```text
(sigma,nu) -> (cell flexes, upper-edge data, strut data)?
```

The claim requires all of the following, not merely dimension equality:

1. derive the local data from the exact Jacobian for arbitrary vertex
   values;
2. decompose every local displacement through the frozen right inverse and
   six-flex kernel;
3. verify every complete rational face equation exactly;
4. prove the 840-by-240 data projection has rational rank 240;
5. combine this rational lower bound with the already frozen modular upper
   bound 240.

## Classification before calculation

- **KNOWN:** endpoint-sum conformal variations and vertex lapse variables are
  standard broad structures.
- **PATTERN:** the observed modular nullity equals `2V`.
- **OPEN:** inclusion of the proposed 240-dimensional rational carrier.
- **OPEN:** exhaustion of the rational kernel.
- **NOT CLAIMED:** a tensor/graviton sector, wave speed, physical tick or
  value of `G`.

## Consequences fixed in advance

If exact inclusion and exhaustion hold, the schedule-free flat-frustum
tangent carrier is precisely two scalar functions on vertices: local scale
and local lapse.  It contains inhomogeneity but no independent upper-edge
shape/tensor direction.  That would close the linear gravitational-wave route
under the stated flat-cell hypotheses, while leaving open theories with
additional first-order connection/area-angle variables or an independently
selected simplicial schedule.

If inclusion fails, the `2V` reading is refuted.  If inclusion holds but does
not exhaust an independently stronger rational upper bound, the carrier
remains **OPEN**.  No outcome is promoted from a count alone.
