# Independent binary E6 audit

This directory is intentionally separate from the earlier E6 implementation.
It uses row-major 16-bit matrices, with entry `(i,j)` stored in bit `4*i+j`.

`audit_reps.cpp` reconstructs the Delsarte duals of the field representative
and the two published non-field representatives, checks their exact rank
distributions, and searches for explicit left/right maps to their transposes.

`classify_spreads.cpp` independently enumerates every normalized 4-dimensional
binary full-rank matrix code through the identity as a 3-space in
`M_4(F_2)/<I>`, then computes its complete left/right equivalence orbits and the
transpose action. This audits the claim that the three named representatives
exhaust all classes.

`kernel_lift_exhaust.cpp` enumerates every codimension-two kernel in a 12-space
as a unique 2-by-12 RREF annihilator. A lift on a kernel is represented by an
arbitrary extension to the full 12-space, using 48 independent Boolean variables;
this avoids choosing or trusting a kernel basis. For each kernel it constructs
the rank-two lift constraints and decides them with a rollback-free affine-equation DPLL:
each forbidden row-space condition is the clause `(linear form 1 = 1) OR
(linear form 2 = 1)`. The two recursive branches are the disjoint exhaustive
partition `form1=1` and `form1=0, form2=1`. A built-in randomized small-instance
test compares this solver to brute force before every run.

Build with, for example:

```sh
c++ -std=c++20 -O3 -pthread audit_reps.cpp -o audit_reps
c++ -std=c++20 -O3 -pthread classify_spreads.cpp -o classify_spreads
c++ -std=c++20 -O3 -pthread kernel_lift_exhaust.cpp -o kernel_lift_exhaust
```

Reproduce the structural checks with `./audit_reps` and
`./classify_spreads`.  Reproduce a complete class-I or class-II solver record
with, for example,

```sh
./kernel_lift_exhaust --class 1 --mode orbits --start 0 --end 3240 --threads 4
./kernel_lift_exhaust --class 2 --mode orbits --start 0 --end 26643 --threads 4
```

For class III use the four half-open ranges recorded in
`orbit_v3_full_class3.aggregate.json`.

The full kernel count is
`[12 choose 2]_2 = 2,794,155`. Runs may be deterministically sharded by the
half-open global index interval `--start N --end M`.

In `--mode orbits`, the program independently enumerates every left/right
automorphism of the chosen spread (testing all 20,160 invertible right factors
and 15 possible images of the identity), verifies closure, derives the
contragredient action on the dual code, and explicitly partitions all 2,794,155
kernels. The same 48-variable solver checks one representative per orbit; the
sum of the reported orbit masses is an exact coverage check.

Version 3 uses the standard FNV-1a 64-bit offset basis
`14695981039346656037` for its deterministic record digest.  Earlier v2 logs
used a custom (mistyped-offset) FNV recurrence and are retained only as
redundant semantic-run evidence.

## Completed v3 audit

The three spread classes have respectively 900, 108, and 18 left/right
automorphisms.  Orbit reduction leaves 3,240, 26,643, and 156,679 kernels.
Every one of these 186,562 representatives is UNSAT.  The summed orbit masses
are 2,794,155 in each class, so the audit covers all 8,382,465 class/kernel
pairs.  There were no SAT results or checker errors.

The class-I and class-II records are `orbit_v3_full_class1.json` and
`orbit_v3_full_class2.json`.  Class III was run in four deterministic shards;
`orbit_v3_full_class3.aggregate.json` gives their totals and binds each shard
by SHA-256.  The aggregate search counters are 230,848,282 DPLL nodes and
530,284,775 forced affine equations.
