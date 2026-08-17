# Cellular weak-lapse jet versus the four committed ticks

Date: 2026-08-17

Prior-art commit: `b77856a`  
Protocol commit: `71d10b4`  
Blind Stage-A result commit: `76a09ab`  
Stage-B verifier registration commit: `6f76f21`

Frozen Stage-B artifact:

```text
reproducible/gravity_600cell_cellular_weak_lapse_comparison.json
SHA-256 bc9dccddb5f4f07bf9904e53eccfc9936aa5bbb5d24d76d35e09d7e679bece58
```

## Verdict

**DERIVED:** `CELLULAR_JET_EXPLAINS_FOUR_TICKS` (`6/6` targeted
checks).

The coefficient jet derived from the closed cellular action and committed in
`76a09ab` explains every disclosed leading observable of the previously
committed four-tick numerical trajectory.  No state was re-solved in Stage B.

## Comparison ledger

### Normalized law

The blind exact prediction was

```text
u_n/u_1       = n,
a_n/u_1       = n(n+1)/2,
v_n/v_1       = 2n-1,
r_n/v_1       = n^2,
p_out,n/k     = 2n+1,
```

for `n=1,2,3,4`.  Stage B made 21 disclosed comparisons: the eleven stored
`n<=3` Richardson intercepts, and five `n=4` intercepts for each independently
constructed even/odd schedule.  Every blind value lies inside its frozen
uncertainty band.  The largest absolute intercept error is

```text
1.07728659518e-11
```

and occurs within the already committed fourth-tick comparison.

### Absolute coefficients

The committed carrier constants give

```text
e0^2 = (tau0/L0)^2
     = 0.0000135225767934642581923492147619670993883...
```

The action therefore predicts, without an adjustable scale coefficient,

```text
u1/lambda^2 = A1 e0^2
             = -3.1160599219482742840...e-6,

v1/lambda^2 = R1 e0^2
             = -3.5592558568436707966...e-6.
```

The Richardson discrepancies from the old numerical solves are respectively

```text
2.44138964237e-21  within band 3.66208415556e-19,
2.97485413615e-21  within band 4.46199046278e-19.
```

These absolute checks carry more information than recognizing triangular or
square integer sequences.

### Direct finite-state equations

All `3 lambdas x 4 slabs x 2 parities = 24` committed states were inserted
directly into the cellular action derivatives, without Newton correction.

```text
max |cellular lapse residual F| = 7.85687026814e-30,
max |cellular seam residual G|  = 6.89457992311e-24.
```

The even and odd stored states agree exactly at their serialized precision.
The comparison uses the closed cellular formula and no staircase carrier.

## What has actually been established

**DERIVED:** the numerical four-tick law was not a separate phenomenon.  It
is the finite-lapse realization of the unique leading canonical jet of the
homogeneous cellular Regge action.

**DERIVED:** close to the time-symmetric static slab, the fixed 600-cell dust
model has nonzero quadratic scale evolution.  In the contracting time
orientation used here, the logarithmic position is triangular, its increment
is linear in the tick index, and the lapse correction is quadratic.
Time reversal gives the corresponding expanding orientation.

**STRUCTURAL:** this is a local discrete Friedmann law around a turning point
for one regular finite carrier.  It is genuine gravitational dynamics in the
minisuperspace sector, but it is not yet a test of general gravitational
degrees of freedom.

**DERIVED CONTROL:** under the preregistered volume-radius map and half-step
convention, the magnitude of the first discrete acceleration coefficient is
`1.078979468...` times the closed continuum Friedmann value.  The fixed
600-cell is therefore about `7.90%` high in this comparison.

## Framing attack and provenance limit

The commit ordering proves a restricted but useful statement: the analytic
artifact was frozen before the numerical artifacts were parsed by Stage B,
and the Stage-B verifier has no coefficients to tune.  It does **not** prove
cognitive blindness.  The integer pattern was already known to the authors,
appears in the prior-art note, and could in principle have influenced choices
of variables, ansatz or branch.

This limitation does not undo the algebraic result: the action was fixed, the
canonical equations have ranks `1` and `2`, their determinant is nonzero, and
both absolute coefficients and all finite-state residuals agree.  But the
honest provenance label is:

```text
STRUCTURAL artifact-level blindness, not a double-blind prediction.
```

A stronger prospective protocol would give the frozen action and canonical
boundary datum to an independent implementation while withholding and
renaming all tick observables, then reveal their hashes only after the exact
jet artifact is signed.  That cannot be manufactured retroactively.

## Relation to prior work

The broad result is **KNOWN**, not a claim of a new theory of gravity:

- Collins and Williams introduced regular 5-, 16- and 600-cell Regge models
  of closed Friedmann dynamics: <https://doi.org/10.1103/PhysRevD.7.965>.
- De Felice and Fabri explicitly evolved dust on the 600-cell with a Sorkin
  scheme: <https://arxiv.org/abs/gr-qc/0009093> and
  <https://arxiv.org/abs/gr-qc/0106077>.
- Dittrich and Hoehn formulate additive one-step Regge actions as generating
  functions for canonical evolution: <https://arxiv.org/abs/1108.1974>.
- Tsuda and Fujiwara derive regular 4-polytopal frustum recurrences and their
  continuum-time Friedmann limit: <https://arxiv.org/abs/2011.04120>.
- Liu and Williams show that finite regular carriers reproduce qualitative
  FLRW evolution but exhibit resolution-dependent discrepancies:
  <https://doi.org/10.1103/PhysRevD.93.024032>.

The post-result searches did not locate the particular exact conserved-dust
weak-lapse coefficient formulas or the artifact-level reconciliation above.
That absence is not proof of novelty.  **External novelty remains OPEN.**

## What this does not establish

**OPEN:** an all-`n` theorem.  Four coefficients support the displayed
closed form, but the registered derivation stops at `n=4`.

**OPEN:** spatial-refinement convergence.  The `7.90%` coefficient excess
could be the expected fixed-resolution error; it cannot be interpreted until
at least one refined carrier is computed.

**OPEN:** anisotropic stability and tensor modes.  Homogeneous evolution has
one scale degree of freedom and cannot demonstrate gravitational waves.

**OPEN:** a physical fundamental tick.  Since the changes scale as
`lambda^2` while proper time scales as `lambda`, the weak-lapse limit gives a
finite acceleration; it does not select an absolute duration.

**OPEN:** a limiting speed, Newton's dimensional constant, Planck time,
Planck mass or particle masses.

## Next discriminating test

The cheap mathematical closure is to prove or refute the formulas for
arbitrary `n` directly from the coefficient recursion.  The next physically
load-bearing test is spatial refinement: compute the same coefficient on a
specified refined `S^3` carrier and ask whether the ratio `1.078979468...`
moves toward `1` under a preregistered radius and dust normalization.  Only
that can distinguish a convergent approximation to Einstein--Friedmann
dynamics from a special fixed-600-cell coincidence.

