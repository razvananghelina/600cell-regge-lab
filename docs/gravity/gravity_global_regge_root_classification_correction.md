# Correction to the boundary-contact count

Date: 2026-08-13

Prior correction commit: `82ab7b8`

The first parallel rerun reproduced every terminal residual exactly but
finished `15/16` because the prior correction overclassified one termination.
The exact reproduced statuses are:

```text
even: 3 artificial-box contacts, 1 causal/angle Jacobian failure,
      2 iteration limits
odd:  4 artificial-box contacts, 2 iteration limits.
```

The exceptional path is `even/S4`.  Its terminal minimum variable is
`0.0814279303`, far above `exp(-6)`, while its minimum angle-argument modulus
has fallen to `0.00149484`; a centered Jacobian point fails a causal/branch
gate.  It must remain `jacobian_failure`, not `boundary_contact`.

Freeze only this expected-status correction.  The second rerun must again
reproduce all twelve terminal residuals to relative `2e-8` and the exact
status counts above.  No mathematical search parameter changes.
