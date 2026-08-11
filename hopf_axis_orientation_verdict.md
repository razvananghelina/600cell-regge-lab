# Static binary and chamber data do not orient the Hopf axes

Date: 2026-08-11

Protocol commit: `1c8694b`.

Registered verifier:
`reproducible/verify_hopf_axis_orientation.py`.
Targeted exact result: `15/15`.

No Dirac matrix, matter character, mass, coupling or Standard-Model target
was evaluated.  The full suite was not run, by explicit user instruction.

## Headline

The orientation needed to escape the residual `C2` is not already hidden in
the binary handedness, the chamber grading or the fibre-edge incidence.

An oriented fivefold axis is the exact two-point homogeneous cover

```text
A5/C5 -> A5/D5.
```

This cover has **no `A5`-equivariant section**.  A reflection in the axis
stabilizer `D5` preserves the underlying axis and all unoriented incidence,
but swaps the two cover points and sends

```text
r -> r^-1
```

for every rotation in its `C5` subgroup.  Consequently any static
`A5`-natural construction must treat the two fibre directions equally.

This is a **DERIVED STATIC-ORIENTATION NO-GO** for all three preregistered
candidates.

## Existing candidates

### Fibre incidence

The earlier golden split correctly found the canonical inverse pair

```text
{r,r^-1}.
```

The present calculation shows why it cannot be sharpened to one member.  The
normalizer reflection is an actual incidence-preserving stabilizer, not an
enumeration convention.  Selecting `r` rather than `r^-1` breaks this exact
symmetry.

### `qH` versus `Hq`

The handed fibration double is

```text
(A5/D5) x {qH,Hq}.
```

There are exactly two equivariant maps from the six-axis set into it: choose
the `qH` copy or choose the `Hq` copy.  Both retain the same unoriented
`A5/D5` label.  Handedness distinguishes left from right cosets; it does not
choose a point in `A5/C5` and therefore does not distinguish `r` from
`r^-1`.

### Oriented chambers

Each chamber-orientation sheet is a free `A5` torsor `A5/1`.  Neither sheet
contains a point fixed by the axis stabilizer `D5`, so there is no
equivariant map from `A5/D5` into either sheet and no canonical chamber over
an axis.

Choosing a chamber would orient an entire flag, but it is a 60-fold
symmetry-breaking choice.  The existence of two chamber chiralities does not
select one of those 60 chambers.

## Complete global orientation census

A simultaneous fibre orientation chooses one of two points over each of six
axes.  All

```text
2^6=64
```

assignments were enumerated under the exact `A5` action.  Their complete
decomposition is

```text
orbit size 20, stabilizer C3:  2 orbits, 40 assignments,
orbit size 12, stabilizer C5:  2 orbits, 24 assignments.
```

Hence

```text
A5-invariant assignments     0/64,
trivial-stabilizer assignments 0/64,
free A5 orbits                  0.
```

This corrects a tempting but false shortcut: even an independently supplied
six-bit global orientation field does not by itself form a complete frame.
Every one of its possible vacua retains either `C3` or `C5`.

Choosing only the orientation of the selected Hopf axis gives one of the
twelve points `A5/C5` and retains stabilizer `C5` exactly.  It removes the
reflection part of `D5`, but still fails the connectedness requirement for a
covariant Dirac operator.

## Homogeneous-set proof

For a transitive `G`-set `G/H`, an equivariant map to a `G`-set `Y` is fixed
by the image of the base point, and that image must be fixed by `H`.

The verifier enumerates every possible base image and reconstructs every map
directly:

```text
target                              D5-fixed points / equivariant maps
A5/C5 oriented axes                                  0
(A5/D5) x {qH,Hq}                                    2
(A5/1) x {chamber +,-}                               0
```

The two middle maps select only the handed copy and project identically back
to the unoriented axis.  Thus none is a section of the oriented-axis cover.

## Hostile framing audit

1. Ambient orientation of three-space does not orient a line.  An axis needs
   a direction or an additional transverse datum.
2. A Hopf fibration has an oriented-circle description only after a `U(1)`
   generator is chosen.  Replacing the generator by its inverse preserves
   the unoriented fibre partition and is implemented here by `D5`.
3. `qH/Hq` is chirality of multiplication, not orientation along the cyclic
   subgroup.
4. Chamber orientation distinguishes the two flag orbits, not one chamber
   or one endpoint over every Hopf axis.
5. The 64-assignment census rules out the idea that an arbitrary global sign
   convention could accidentally break all residual symmetry.

## Status ledger

- **DERIVED:** `A5/C5 -> A5/D5` is an equivariant two-point cover.
- **DERIVED:** a `D5` reflection swaps its fibres and implements
  `r <-> r^-1`.
- **DERIVED NEGATIVE:** the cover has zero `A5`-equivariant sections.
- **DERIVED:** `qH/Hq` supplies exactly two handed-copy maps and zero axis
  orientations.
- **DERIVED NEGATIVE:** neither free chamber sheet supplies an equivariant
  axis section.
- **DERIVED:** the 64 global orientation assignments form two size-20 and
  two size-12 orbits, with stabilizers `C3` and `C5`.
- **DERIVED NEGATIVE:** zero assignments have trivial stabilizer.
- **DERIVED:** one oriented Hopf axis retains `C5`.
- **OPEN:** a second correlated order parameter whose common stabilizer with
  an oriented Hopf axis is trivial.
- **NO DIRAC TARGET:** no operator was selected or tested.

## Programme boundary

The static six-fibration geometry supplies neither a canonical direction
along the fibres nor enough orientation data to make a complete frame.

The smallest remaining possibility must combine two independently selected
pieces, for example an oriented Hopf point (`C5` stabilizer) and a threefold
axis (`D3` stabilizer), whose common stabilizer is necessarily trivial.  The
complete correlated-pair orbit and look-elsewhere census must precede any
Dirac construction, and the current action still supplies neither datum
simultaneously.
