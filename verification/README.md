# Verification package

All mathematical decisions in this directory use exact arithmetic over
`F_2`.  The package has four independent layers.

- `dim11/` contains the explicit lower-bound certificate and two verifiers.
- `family/` contains the exact `E_4` and `E_5` certificates, the recursive
  row-cone verifier, and an independent exhaustive tiny-instance check.
- `structural/` contains the four residual field-class instances left by the
  rank-two filter, together with independent SAT logs and proof certificates.
- `primary/` contains the complete orbit-reduced 40-variable enumeration.
- `independent/` reconstructs the spread classification and kernel orbits in
  a separately written 48-variable implementation.
- `raw_cuda/` enumerates all 8,382,465 raw kernels without orbit reduction.
- `cnf_samples/` contains three additional worst-node cross-checks, one from
  each MRD class.

Run the standard-library audit from the repository root:

```bash
python3 verification/verify_release.py
```

This recomputes the 2,048-word lower-bound rank distribution, streams and
cross-checks every primary orbit/result record, checks the independent shard
coverage, validates all raw-kernel terminal summaries and deterministic
digests, checks the residual solver records, and verifies every SHA-256 hash
in `audit_manifest.json`.

Run the family checks separately:

```bash
python3 verification/family/verify_small_examples.py
python3 verification/family/verify_family.py
python3 verification/family/bruteforce_tiny.py
```

The standard audit validates the compressed proof files and their recorded
terminal checks.  A full replay of all four proofs is available separately:

```bash
python3 verification/structural/verify_drat_proofs.py
```

The C++ and CUDA sources are included so that each computation can be rebuilt.
Compiled binaries are intentionally excluded from the release.
