# A counterexample to the Etzion--Silberstein conjecture

This repository contains the paper, exact certificate, and independent
verification programs for the binary Ferrers diagram

```text
column heights: (5, 5, 5, 5, 1, 1)
minimum rank distance: 3
```

The Etzion--Silberstein bound is 12, but the largest possible dimension is
exactly 11.  Thus this instance disproves the conjecture.

The paper also proves a row-cone propagation identity.  Iterating it from
this seed gives explicit binary counterexamples with Singleton-type bound 12
and exact optimum 11 at every minimum rank distance `d >= 3` (and likewise
for the transposed diagrams).

## Fast verification of the lower bound

The explicit 11-dimensional code is stored as eleven generators in
`verification/dim11/e6_dim11_certificate.json`.  Two separate verifiers check
the Ferrers support, linear independence, and the ranks of all 2,048 words.

```bash
python3 verification/dim11/verify_e6_dim11.py

c++ -O3 -std=c++20 \
  verification/dim11/verify_e6_dim11_independent.cpp \
  -o /tmp/verify_e6_dim11
/tmp/verify_e6_dim11 \
  verification/dim11/e6_dim11_certificate.json
```

Both report minimum nonzero rank 3 and rank distribution

```text
rank 0:    1
rank 3:  605
rank 4: 1098
rank 5:  344
```

## Upper-bound proof

A hypothetical 12-dimensional code projects injectively onto a binary
`[4 x 4,12,2]` MRD code `U`.  The two short columns define a codimension-two
kernel `K <= U`; the remaining top row must avoid the row space of every
rank-two matrix in `K`.

The Frobenius dual `U^perp` is a four-dimensional space in which every
nonzero matrix is invertible.  There are exactly three left--right equivalence
classes of these binary spread codes.  Each class has

```text
[12 choose 2]_2 = 2,794,155
```

possible kernels.  A short structural lemma says that any feasible kernel
must contain exactly 105 rank-two matrices.  Exact orbit enumeration then
eliminates both nonfield classes and all but four field-class kernel orbits:

```text
class       kernel orbits   minimum rank-2 count   survivors   raw mass
field               3,240                    105           4         70
II                 26,643                    117           0          0
III               156,679                    113           0          0
```

The four surviving exact lift systems are all unsatisfiable.  They are also
provided as ordinary DIMACS instances with independently checked proof
records in `verification/structural`.

The package contains three deliberately redundant full checks:

1. `verification/primary`: an orbit-reduced 40-variable solver.  It checks
   186,562 orbit representatives whose recorded masses cover every kernel.
2. `verification/independent`: a separately written 48-variable solver that
   independently reconstructs the classification, automorphism groups, and
   kernel orbits.
3. `verification/raw_cuda`: a 48-variable checker that enumerates every raw
   kernel in every class, without orbit reduction.  It reports 8,382,465
   `UNSAT`, zero `SAT`, and zero errors.

The hardest representative in each class was also translated independently
to ordinary DIMACS and checked by both CaDiCaL and CryptoMiniSat; those inputs
and logs are in `verification/cnf_samples`.

All rank computations, parity systems, and searches use exact arithmetic over
`F_2`.  No floating-point computation is involved.

The same result disproves the binary `(n,d)=(5,3)` case of the
puncturing--inclusion MRD conjecture formulated by Couvée and Neri in 2026.

## Counterexamples at every distance

For a diagram `D=(c_1,...,c_n)` at distance `d`, let
`b=nu_min(D,d)>0`.  Its row cone is

```text
R(D) = (c_1+1, ..., c_n+1, 1 repeated b times).
```

The paper proves exactly

```text
nu_min(R(D),d+1) = nu_min(D,d)
kappa_F(R(D),d+1) = kappa_F(D,d),
```

where `kappa_F` is the largest possible dimension over the field `F`.
Starting from `(5,5,5,5,1,1)` over `F_2` therefore preserves the one-unit
gap while raising the distance at each step.

Two portable checks accompany the proof:

```bash
python3 verification/family/verify_family.py
python3 verification/family/bruteforce_tiny.py
```

The first constructs eleven stages and checks every word, support condition,
Singleton cut, and shifted rank distribution.  The second independently
enumerates every binary subspace in a tiny exhaustive test range.

## Paper

The manuscript source is [`main.tex`](main.tex), with a compiled copy in
[`main.pdf`](main.pdf).

## Reproducibility and integrity

`verification/audit_manifest.json` records all class counts, coverage totals,
terminal statuses, deterministic result digests, and SHA-256 hashes.  Run

```bash
python3 verification/verify_release.py
```

to validate the complete release.  To compile the bundled independent proof
checker and replay all four compressed DRAT certificates, also run

```bash
python3 verification/structural/verify_drat_proofs.py
```

## Reference

Jitendra Prajapati, *A counterexample to the Etzion--Silberstein conjecture*,
2026.
