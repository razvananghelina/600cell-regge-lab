# Preregistered correction: directional diagnostic complex parser

Date: 2026-08-21

First diagnostic crash preserved in commit `dd7e05a`.

Add one parser for the frozen textual form `(real +/- imagj)`:

1. remove the outer parentheses;
2. split from the right into the exact decimal real field, sign and imaginary
   field;
3. remove `j`, attach the sign, and parse both fields independently with
   `mp.mpf`;
4. return `mp.mpc(real,imag)`.

Use it only for the frozen first-Richardson reproduction control.  Test it on
both `+` and `-` synthetic strings before use.  Do not modify or re-evaluate
the scientific ladder definition, steps, precision, envelopes, ratios,
corruption, outcome hierarchy or scope.

