# Prior-art and framing gate: can the prism shift become a dynamical field?

Date: 2026-08-19

Status: **completed before an unequal-scale extension verifier or action
calculation**.

## 1. Question and complete hypotheses

The preceding exact result gives a 119-dimensional shape-matched potential
family on an equal-scale 600-cell slab and the restricted action Hessian

```text
H_phi = [2*pi-5*acos(1/3)]/(L*sqrt(rho)) * Delta_0.
```

Before interpreting that spatial Hessian as one half of a wave operator, ask
whether the same variables exist on a genuinely dynamical homogeneous slab.

The local test uses:

1. a centered regular tetrahedron with vertices `b_i` in Euclidean three-
   space, equal radii and affine span three;
2. bottom vertices `B_i=b_i` and top vertices

   ```text
   T_i=q*b_i+s+N*n,
   ```

   where `q>0` is the top/bottom scale ratio, `s` is tangential and `n` is a
   unit timelike normal;
3. planar lateral quadrilaterals, which this homothetic translated
   realization satisfies;
4. one common timelike squared strut length
   `(T_i-B_i)^2=-rho`, `rho>0`, for all four corresponding vertex pairs;
5. no cosmological, action or desired propagation target in the kinematic
   test.

This is the local building block of the homogeneous common-lapse carrier.
It is not a theorem about arbitrary nonhomothetic boundaries, four
independent strut lengths, a general ADM shift or the already triangulated
2,280-edge action carrier.

## 2. Framing attack before calculation

For every vertex,

```text
ell_i^2 = |(q-1)b_i+s|^2-N^2.
```

Because all `|b_i|` are equal, subtraction gives

```text
ell_i^2-ell_j^2
  =2(q-1) s.(b_i-b_j).
```

The three independent tetrahedral edge directions span the tangential
space.  Therefore common struts appear to impose

```text
(q-1)s=0.
```

If this elementary derivation survives an exact coordinate and rank audit,
the admissible set is not a smooth `(q,phi)` configuration space.  It is the
union of two branches:

```text
equal scale q=1:       s arbitrary,
unequal scale q!=1:    s=0.
```

Their intersection at `(q,s)=(1,0)` is singular.  Its linearized tangent
space contains both the scale and shift directions because the constraint is
bilinear, but a generic mixed tangent does not integrate even to second
order.  A Hessian evaluated on the large linear tangent space would not by
itself define a two-variable dynamics.

This observation attacks the final paragraph of
`gravity_600cell_prism_shift_action_result.md`: merely "releasing equal
scale" may destroy the shift carrier instead of producing its temporal
kinetic block.

## 3. What changes with nonuniform struts

If the four strut squares are allowed to differ, the same identity reads

```text
ell_i^2-ell_j^2=2(q-1)s.(b_i-b_j).
```

At `q!=1`, three independent strut differences determine the three
components of `s`.  Thus the shift is no longer an independent invisible
cell datum; it is encoded in inhomogeneous lapse/strut data.  The inverse map
is singular as `q` approaches one.

If confirmed, the correct dynamical arena is therefore the complete internal
edge carrier, or an explicitly derived nonuniform-strut reduction.  It is
not a freely propagated copy of the equal-scale vertex potential on every
tick.

## 4. Prior art

Homogeneous Collins--Williams and polytopal cosmologies use parallel or
homothetic regular spatial cells with common struts and thereby select a
shift-free sector by ansatz:

- R. Tsuda and T. Fujiwara, *Oscillating 4-Polytopal Universe in Regge
  Calculus*, [arXiv:2011.04120](https://arxiv.org/abs/2011.04120), DOI
  `10.1093/ptep/ptab079`;
- R. G. Liu and R. M. Williams, *Regge calculus models of the closed vacuum
  Lambda-FLRW universe*,
  [arXiv:1501.07614](https://arxiv.org/abs/1501.07614);
- B. Dittrich, S. Gielen and S. Schander, *Lorentzian quantum cosmology goes
  simplicial*, [arXiv:2109.00875](https://arxiv.org/abs/2109.00875), DOI
  `10.1088/1361-6382/ac42ad`.

Canonical linearized Regge calculus distinguishes lapse/shift-like vertex
displacements from curvature modes, and warns that the available variables
and constraints depend on the move and background:

- P. A. Hoehn, *Canonical linearized Regge Calculus: counting lattice
  gravitons with Pachner moves*,
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672);
- B. Bahr and B. Dittrich, *(Broken) Gauge Symmetries and Constraints in
  Regge Calculus*, [arXiv:0905.1670](https://arxiv.org/abs/0905.1670).

The repository's exact prism-rigidity theorem already found three
non-isometric translations at equal scale and none at the tested unequal
scales.  The proposed calculation sharpens that finite rank census into the
explicit branch equation and its global consequence.

The focused search used `Regge hyperfrustum shift`, `unequal scale struts`,
`Collins-Williams shift lapse` and `600-cell unequal strut evolution`.  It
located the standard homogeneous and canonical mechanisms above, but no
source stating this exact regular-tetrahedron branch equation or its
600-cell potential consequence.  Search absence is not proof; external
novelty is **OPEN**.

## 5. Frozen questions for the protocol

1. Does an exact coordinate calculation give
   `(ell_i^2-ell_0^2)=2(q-1)s.(b_i-b_0)` for all three independent rows?
2. Does the `3 x 3` edge-direction matrix have nonzero determinant?
3. Is the common-strut solution ideal exactly the ideal generated by the
   three components of `(q-1)s`, up to an invertible linear change?
4. Does the Jacobian at `(q,s)=(1,0)` vanish in the mixed constraints while
   the quadratic tangent cone rejects a generic simultaneous scale/shift
   direction?
5. On the complete connected 600-cell, does `q!=1` force all local shifts
   and hence all potential differences to zero?
6. With nonuniform struts, do three independent length differences recover
   `s` uniquely for `q!=1`, with inverse proportional to `1/(q-1)`?

## 6. Decision boundary

If all six points pass, report an **EXACT DYNAMICAL-EXTENSION OBSTRUCTION**:
the equal-scale `phi` Hessian is a spatial stiffness on a singular static
branch, not yet the spatial half of a propagated scalar field.  Correct the
previous next-step wording and move to the complete internal-edge canonical
carrier.

If a nonzero `s` survives at `q!=1` with common struts, the obstruction is
refuted and a genuine unequal-scale action `S(q_-,q_+,rho,phi)` becomes the
next calculation.

If only nonuniform struts rescue `s`, record that result without calling the
resulting derived variable an independent shift field.

No value of `c`, dispersion relation or continuum spectrum is loaded in this
gate.
