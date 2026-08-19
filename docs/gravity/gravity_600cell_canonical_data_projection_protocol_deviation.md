# Canonical-data projection protocol deviation

Date: 2026-08-19

This deviation is recorded while the first projection verifier is still
reconstructing its upstream baseline matrices, before it has printed a partial
projection rank or written its JSON artifact.

## Deviation

The preregistered protocol forbids comparison of the unknown projection
dimensions with 119, 120, or combinations suggested by prior carriers until
the first JSON is committed.  During the running calculation, the assistant
used `119+121` as an example of precisely the post-hoc pattern the protocol was
designed to prevent.  That sentence itself discloses a possible target before
the artifact and violates the strict target-blindness rule.

## What remains externally fixed

The prior-art gate (`faf5b08`), protocol (`36ddebd`), complete verifier source,
registry entry, and every outcome/control criterion (`d259477`) were committed
before this disclosure.  No projection value had been emitted.  The running
source must not be modified, and Git therefore still proves that the operator,
column partitions, complete convention list, and rank formulas preceded the
result.

## Consequence

The first artifact remains a preregistered **STRUCTURAL** modular census, but
must not be advertised as perfectly target-blind.  Any later numerical
coincidence with 119, 120, 121, or their combinations is at most **PATTERN**
until a mechanically independent adversarial construction proves the relevant
subspaces directly.  Restarting the same calculation cannot repair blindness,
because the target is now known; it would only hide the chronology.

This note changes no source, matrix, criterion, or outcome.  It strengthens the
acceptance gate for any subsequent carrier claim.
