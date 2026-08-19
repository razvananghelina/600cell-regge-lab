# Disclosed correction after the first dynamic-shift extension run

Date: 2026-08-19

First-result commit: `fcb3a14`.

First artifact SHA-256:

```text
134f8e68335f3acdd40eb909d4dcb4fae4361329c4540a65182fe061affb499a
```

The first registered execution returned

```text
12/13 PASS
BRANCH_GEOMETRY_OPEN.
```

The failed gate printed

```text
c=(-4*a*s2-4*a*s3,
   -4*a*s1-4*a*s3,
   -4*a*s1-4*a*s2),
E^-1*c/2=(a*s1,a*s2,a*s3).
```

An isolated exact inspection found that the two compared objects were a
mutable and an immutable SymPy matrix with algebraically equal but
structurally different entries:

```text
-4*a*s2-4*a*s3
4*a*(-s2-s3).
```

Container equality returned false, while all three simplified entrywise
differences were exactly zero.  The independently tested transformed identity
`E^-1*c/2=a*s` had already passed in the same failed run.

The correction is frozen before rerun:

1. replace container equality by the exact condition that every simplified
   entry of `c-2*a*E*s` equals zero;
2. retain the transformed identity and every other gate unchanged;
3. add no tolerance and change no predicted verdict.

This is a mechanical symbolic-comparison repair, not a relaxation of the
scientific condition.  The failed artifact remains in git history.
