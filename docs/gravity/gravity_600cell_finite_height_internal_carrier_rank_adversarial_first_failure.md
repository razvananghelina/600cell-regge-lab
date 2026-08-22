# First adversarial internal-rank execution failure

Date: 2026-08-22.

Implementation commit executed: `5065a7c`.

The first targeted execution passed provenance, registry, geometry,
finite-height background, independent derivative-branch, stationarity,
physical row-split, synthetic-classifier and hostile-corruption controls for
the even schedule.  Its direct global census printed

```text
diagonal rank = 119
full rank     = 239
```

It then completed the even nonlinear secant evaluations but stopped before
forming their verdict, before the odd reconstruction and before reading the
primary scientific artifact.  The exact exception was

```text
TypeError: '<' not supported between instances of 'str' and 'float'
```

The cause was operational: the computed prediction error had been formatted
as a JSON string before applying the already frozen numeric `1e-5` gate.
The correction retains the value as a float through classification and
formats it only while constructing the public record.  No formula, step,
threshold, rank rule, matrix or physical interpretation is changed.
