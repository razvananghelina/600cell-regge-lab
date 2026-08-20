# First execution failure: nonlinear time-reversal anchor control

Date: 2026-08-20

Protocol commit: `722fb3c`  
Implementation commit: `36ac45b`

The first targeted execution stopped before every action root search. It
passed `10/11` checks and wrote artifact SHA-256

```text
ac759653317bc3e85e1e126885692c64462350468ac12710d0656bf9e402bf25
```

The failed check was the preregistered nonlinear time-reversal reduction. At
least one of the two displaced control anchors was not classified on the
same admissible Lorentzian branch as the inherited fill, so the implementation
assigned an infinite comparison rather than comparing complex action rows.

No main-box or boundary-ladder attempt ran:

```text
main attempts   = 0/72
ladder attempts = 0/48
```

Therefore `REFINED_H4_STATIONARY_ROOT_CONTROL_FAILED` is an infrastructure /
control result and contains no evidence for or against a finite stationary
root.

The narrow repair must first print the finite/branch diagnostics for every
predeclared anchor and schedule pair. It may either replace branch-invalid
anchors by smaller predetermined displacements or compare identical invalid
diagnostics, but the choice and new anchors must be committed before the root
search is rerun. The 120 scientific attempts, bounds, seeds, solver and root
gates may not change.
