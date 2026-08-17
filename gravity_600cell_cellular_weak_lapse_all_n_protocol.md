# Preregistration: all-index cellular weak-lapse coefficient theorem

Date: 2026-08-17

Prior-art gate commit: `bd1e7df`

Status: frozen before constructing or running the generic-index verifier.

## 1. Inputs and scope

The verifier may use only the closed cellular action and conventions frozen
in the Stage-A artifact

```text
reproducible/gravity_600cell_cellular_weak_lapse_blind.json
SHA-256 6d39e9a4594d9c9ead102f94cf9115d8474132ecce511fe7359826dcc73b9de0
```

and may reconstruct its symbolic derivatives.  The finite tick artifacts are
not inputs.  This is no longer a prediction-blind calculation: the candidate
closed forms are disclosed below and must be treated as hypotheses to be
falsified by a generic-index identity.

The statement concerns the formal coefficient jet at each fixed integer
`n>=1`.  It does not claim uniform convergence when `n` grows with `1/e`.

## 2. Frozen definitions

Put

```text
epsilon = 2*pi-5*acos(1/3),
q       = 5*sqrt(2)-3*epsilon,
x       = e^2,
t_n     = n(n+1).
```

Use exactly the action, lapse equation and additive seam equation from
`gravity_600cell_cellular_weak_lapse_protocol.md`.  Set

```text
A_-1=A_0=B_-1=B_0=R_0=0.
```

The candidate formulas are

```text
A_n = -6*epsilon*t_n/q,

R_n = -10*epsilon*(7*sqrt(2)*epsilon+60)*n^2/q^3,

B_n = -epsilon^2*t_n/q^4 *
      [(108*epsilon^2-395*sqrt(2)*epsilon+300)*t_n
       -54*epsilon^2+145*sqrt(2)*epsilon-600],

p_out,n/k = 2*n+1.
```

`A_n`, `R_n` and the momentum pattern are suggested by the exact four-step
result.  `B_n` was found by the procedurally disclosed interpolation in the
prior-art note.  None counts as evidence here.

## 3. Generic induction test

Reconstruct the universal action jet through the orders already certified in
Stage A, but keep `n` symbolic.  For a generic step:

1. substitute the candidate formulas at `n-2` and `n-1`;
2. retain independent unknowns `(A_n,B_n,R_n)` at the new endpoint;
3. extract the first nonzero lapse and seam coefficients;
4. prove that the leading seam equation is affine in `A_n` and uniquely
   selects the displayed `A_n`; the leading lapse equation may contain an
   additional algebraic factor, so require it to vanish on the seam root and
   require the joint local Jacobian there to have rank one;
5. substitute that value and prove that the next equations are affine in
   `(B_n,R_n)`, with exact rank two and determinant
   `16200*epsilon^2`;
6. prove that their unique Cramer solution is exactly the displayed
   `(B_n,R_n)`;
7. prove that the leading outgoing momentum is `(2*n+1)k`.

The verifier must reduce every residual exactly in
`Q(sqrt(2))(epsilon,n)`.  Checking many numerical integers, polynomial
interpolation, or simplification after substituting `n=1..4` is not a proof.

### 3.1 Correction frozen before implementation

Inspection of the already committed Stage-A equations after the first
protocol commit showed that the lapse equation is generally nonlinear in
`A_n` (already at `n=1` it contains `A_1` times the factor selected by the
seam equation).  Therefore the original phrase “the leading equations are
affine” was false.  No all-index verifier had yet been constructed or run.
The corrected test above is the actual property needed for a unique common
branch and is strictly capable of refuting the candidate.

## 4. Domain and branch

Verify analytically:

```text
epsilon>0,
q>0,
det=16200*epsilon^2>0.
```

Permitted elementary controls include

```text
cos(2*pi/5)<1/3  => acos(1/3)<2*pi/5 => epsilon>0,
acos(1/3)>pi/3   => epsilon<pi/3,
pi<5*sqrt(2)     => q>0.
```

Consequently `A_n<0` for every `n>=1`, selecting the contracting formal
orientation.  The time-reversed branch is not separately enumerated.

## 5. Independent controls

After the generic identities are established:

- specialize them at `n=1..4` and require exact equality with the frozen
  blind coefficients;
- evaluate the candidate at the out-of-training indices `n=5,7,11` and at
  `e=1/200` in the full unexpanded cellular equations;
- require the lapse and seam residuals to have the same first-omitted orders
  (`e^7` and `e^5`) under halving to `e=1/400` and `e=1/800`;
- report expression operation counts and runtime so computational
  intractability is not disguised as a negative theorem.

The finite-index controls cannot rescue a failed generic identity.

## 6. Outcomes

### `CELLULAR_WEAK_LAPSE_ALL_N_PROVED`

Report **DERIVED** only if every generic rank, determinant, Cramer,
coefficient and momentum identity vanishes exactly and all independent
controls pass.

### `CELLULAR_WEAK_LAPSE_ALL_N_REFUTED`

Report **DERIVED NEGATIVE** if any generic identity has a certified nonzero
residual.  Print the first residual and whether `A`, `B`, `R` or momentum
fails.

### `CELLULAR_WEAK_LAPSE_ALL_N_OPEN`

Report **OPEN** if the generic algebra cannot be decided within the registered
implementation.  A timeout is not evidence for or against the theorem.

## 7. Interpretation boundary

A positive result proves a formal all-index recurrence for the homogeneous
weak-lapse coefficient jet.  It does not prove:

- convergence for finite `e` at arbitrarily late index;
- continuum Einstein dynamics;
- a physical absolute tick;
- anisotropic or propagating gravitational degrees of freedom;
- external novelty.
