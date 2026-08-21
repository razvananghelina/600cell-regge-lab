# First execution of the adversarial threshold correction

Date: 2026-08-21

Implementation commit `24560c8` changed the componentwise primary comparison
to the preregistered frozen schedule envelopes.  The run correctly reported

```text
24/24 component rows inside their envelopes,
maximum error/envelope ratio = 0.00066666667.
```

It nevertheless retained the formal `DISAGREEMENT` outcome because one
downstream assignment still read

```text
rank_one = nonzero and all(error < 1e-68 ...),
```

rather than reusing the corrected `primary_match` boolean.  Thus the old
impossible threshold survived in exactly one classifier line.  The artifact
is preserved with SHA-256

```text
fe8890cf7bc8b3a393a366a7b5ca49786bcbe2ea017e46ff1f1da371ef0b64c8.
```

The already preregistered correction in `a92ab1f` uniquely fixes this:

```text
rank_one = nonzero and primary_match.
```

No further protocol or numerical change is allowed.

