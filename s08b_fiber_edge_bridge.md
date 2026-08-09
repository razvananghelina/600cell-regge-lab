# S08b Bridge Program: From Fiber Module to Edge Space

## Exact claim

Construct or exclude a natural `A_5`-equivariant map
\[
\Psi:\mathbb{R}^{12}\longrightarrow \mathbb{R}^{720}
\]
from the 12-dimensional fiber permutation module to edge space such that:

1. the image is nonzero;
2. the image is contained in `\ker(\Box_1)`;
3. ideally, the image equals the 12-dimensional nontrivial part of
   `\ker(\Box_1)`.

## Why this step is needed

`S08` closes two separate exact facts:

- the edge-kernel decomposition
  \[
  \ker(\Box_1)=\rho_0\oplus 2\rho_5;
  \]
- the fiber permutation decomposition
  \[
  \mathbb{R}^{12}\cong \mathbf{1}\oplus \mathbf{3}\oplus \mathbf{3'}\oplus \mathbf{5}.
  \]

These live in different ambient spaces and are not automatically the same
module. Therefore the Lie-algebra candidate theorem in `S09` cannot yet be used
as part of the main derivation chain.

## Formal setting

- source: functions on the 12 Hopf fibers;
- target: functions on the 720 edges of the 600-cell;
- symmetry: left action of `2I` on fibers and on edges, factored to `A_5` on
  the source;
- operator constraint: target image should lie in `\ker(\Box_1)`.

## Natural candidate maps to test

1. constant-on-fiber-edge lift:
   assign to each fiber label a uniform function on the 10 edges of that fiber;
2. cross-incidence lift:
   lift a fiber-label function to cross edges according to their endpoint fiber
   labels;
3. incidence-averaged lift:
   average fiber values over the two endpoint fibers of each edge;
4. any other map built functorially from the fiber partition and edge
   incidence, without arbitrary basis choices.

## Success criteria

A candidate bridge is acceptable only if:

1. it is defined canonically from the fibration data;
2. it is `A_5`-equivariant;
3. its image is nonzero and lies in `\ker(\Box_1)`;
4. if it is claimed to identify the 12-dimensional spaces, its image must have
   dimension 12.

## Failure criterion

If every natural `A_5`-equivariant candidate map has zero image inside
`\ker(\Box_1)`, then the bridge fails and the fiber permutation module cannot
be used as the source of the `S09` Lie-algebra candidate inside the main chain.

## Decision

`S08b` is now the blocking step before `S09`.
