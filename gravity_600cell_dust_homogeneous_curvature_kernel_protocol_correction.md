# Protocol implementation correction: kernel-line coordinate and Schur controls

Date: 2026-08-17

Original protocol: `53dc168`.
Angular-resolution clarification: `c2cbcd3`.
Registered implementation: `21d0451`.
Preserved first run: `c048518`.

The first targeted run ended

```text
11/14 PASS
HOMOGENEOUS_CURVATURE_KERNEL_CONTROL_FAILED
```

before an admissible physical outcome.  Two implementation controls, not two
candidate comparisons, failed.

## 1. Direct-eigenvector versus Schur-plane control

The protocol required this discrepancy to be calibrated by tangent-ball,
conditioning, spectral separation and binary eigensolver error.  The code
instead inserted an undocumented absolute cutoff `distance < 1e-6`.  The
observed distances were `8.52e-6` and `1.83e-5`, so that accidental cutoff
failed even though the near-`-1` selection gap was about 159.

The correction removes the undocumented cutoff and implements the declared
calibration:

```text
epsilon_plane =
    (tangent_ball_radius
     + eps_binary*dimension*||T||_2)
    * condition(eigenvectors)
    / selected_to_unselected_spectral_separation.
```

The outer acceptance comparison already supplies the factor 10.  Including
another factor 10 inside the error would count it twice.  The control passes
only when the direct/Schur plane distance is at most `10*epsilon_plane` and
`epsilon_plane < 1e-2`.  The already-preregistered
`1e-2` angular-resolution cap remains in force.  This is necessarily a
post-failure implementation correction, but it restores the formula stated
before the run rather than choosing a threshold from the desired candidate
label.

## 2. Literal common boundary ordering

The code required the two internally sorted 30-orbit lists to be equal entry
by entry.  They are not: the schedule-dependent sort moves the eleven orbit
types at even positions 3--11 to odd positions 12--20.  A target-independent
literal edge-set comparison gives the unique type map

```text
even -> odd:
[0,1,2,12,13,14,15,16,17,18,19,20,3,4,5,6,7,8,9,10,11,21,22,23,24,25,26,27,28,29].
```

Within every matched orbit the 24 group-coordinate edge labels agree exactly;
there is one common identity group permutation.  The correction derives this
permutation solely from equality of literal 24-edge sets, requires uniqueness,
and applies it to both position and momentum halves before the cross-schedule
line comparison.  No kernel component or candidate distance enters the map.

## What is unchanged

- all frozen input hashes;
- the reconstructed `rank 59 / nullity 1` maps;
- the ten candidates and 20-attempt ledger;
- every candidate distance and uncertainty formula;
- the `1e-2` angular-resolution cap;
- the outcome hierarchy;
- the ban on fitted schedule conjugacies.

The failed JSON remains in git history at `c048518`.  The corrected run must
still be reported as a second run, not silently substituted for the first.
