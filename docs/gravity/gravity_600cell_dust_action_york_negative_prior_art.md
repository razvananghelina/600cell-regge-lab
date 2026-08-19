# Prior-art and framing gate: action-weighted longitudinal test of the negative shape modes

Date: 2026-08-19

## Question

On the fixed regular embedded 600-cell, do the two already selected
rank-`15` negative-stiffness shape spaces coincide with the tangential
vertex-displacement images after the conformal component is removed by the
action-derived kinetic form?

This is a target-disclosed explanatory test.  The negative spaces and their
dimensions are known before this gate.  It is not a blind discovery protocol.

## Frozen input facts

The repository already certifies:

1. the edge carrier has dimension `720`;
2. the vertex-conformal map `C` has rank `120`;
3. the centered kinetic form `H` is nondegenerate, positive on `im C` and
   negative on the action-selected shape factor

   ```text
   S_H = ker(C* H),              dim S_H=600;
   ```

4. the embedded tangential rigidity differential `D` has rank `354`;
5. `dim(im C intersection im D)=4`;
6. the Euclidean framework self-stress `ker R*` is not dynamically closed;
7. two one-dimensional binary-tetrahedral sectors, numbered `4` and `5` in
   the frozen order, each contain a resolved `15`-dimensional negative
   stiffness space and a `10`-dimensional positive space;
8. those negative ranks persist at the next centered recurrence, but their
   naive and generalized phase lifts do not propagate as closed fibers.

The equality of `15` with the available tangential-shape dimension in the two
selected sectors is already visible.  Dimension equality alone is forbidden
as evidence.

## Canonical candidate and complete hypotheses

Let

```text
P_S = I-C(C* H C)^(-1)C* H
```

be the `H`-orthogonal projection along `im C` onto `S_H`.  Let `D` be the
literal tangent-vertex rigidity differential of the unit-quaternion embedding,
not the full ambient rigidity matrix.  Define

```text
L_H = im(P_S D) subset S_H.
```

This construction uses no spectral target.  It is fixed by:

- the declared spherical 600-cell embedding;
- the literal tangent projectors `I-x_v x_v*`;
- the action-derived centered kinetic bilinear form `H`;
- the exact conformal incidence `C`.

The expected global dimension follows before any stiffness comparison:

```text
dim L_H = rank D-dim(im C intersection im D)=354-4=350.
```

Its `H`-orthogonal complement inside `S_H` has dimension `250`.  These counts
parallel the rigidity/self-stress counts but the subspaces need not equal the
Euclidean `im R` and `ker R*` carriers already refuted.

## Prior art and the scope trap

Hoehn's canonical linearized Regge analysis identifies vertex displacements
as gauge/lapse-shift directions and lattice gravitons as curvature degrees of
freedom around **flat** backgrounds:

- P. A. Hoehn, *Canonical linearized Regge Calculus: counting lattice
  gravitons with Pachner moves*, arXiv:`1411.5672`.

The present dust background is curved and its complete finite Legendre block
is regular.  Bahr--Dittrich and Dittrich--Hoehn explain that curved discrete
backgrounds generically break the exact vertex-displacement symmetry and
replace constraints by background-dependent pseudo-constraints:

- B. Bahr and B. Dittrich, *Broken Gauge Symmetries and Constraints in Regge
  Calculus*, arXiv:`0905.1670`;
- B. Dittrich and P. A. Hoehn, *From covariant to canonical formulations of
  discrete gravity*, arXiv:`0912.1817`.

Hartle--Miller--Williams further show that signs of the Lund--Regge
supermetric do not by themselves separate gauge from physical directions and
can depend on the triangulation:

- J. B. Hartle, W. A. Miller and R. M. Williams, *Signature of the
  Simplicial Supermetric*, arXiv:`gr-qc/9609028`.

Therefore even an exact equality

```text
negative stiffness space = L_H
```

would establish an action-weighted **longitudinal structural identity**, not
an exact gauge theorem.  It would be evidence against interpreting these
thirty coarse modes as tensor gravitons, but it would not authorize quotienting
them from the curved finite dynamics.

Conversely, failure of equality would show only that the negative carrier is
mixed relative to this declared action/embedding split.  It would not make it
physical.

## Required hostile controls

A valid calculation must:

1. enumerate `dim L_H` and its complement in all seven minimal sectors before
   comparing the selected two;
2. verify directly that `P_S C=0`, `C*H P_S=0` and `P_S^2=P_S` within the
   inherited error model;
3. compare projectors/principal angles, not dimensions;
4. test the stiffness cross block between `L_H` and its `H`-orthogonal
   complement;
5. test the sign of both restricted stiffness forms;
6. use both schedules and all four derivative variants;
7. include a same-dimension rotated control which must fail the equality;
8. load no continuum harmonic, polarization, speed or mass target.

The two possible scientific outcomes are:

- `NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_RESOLVED` if subspace equality,
  cross-block closure and the `15 negative / 10 positive` longitudinal split
  all resolve;
- `NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_REFUTED_OR_OPEN` otherwise.

External novelty is **OPEN**.  A literature search found no primary source
performing this exact 600-cell comparison, but search absence is not proof.

