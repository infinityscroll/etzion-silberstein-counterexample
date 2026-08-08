# Row-cone family verification

The paper proves that the row cone

`R(D) = (c_1 + 1, ..., c_n + 1, 1^b)`, where `b = nu_min(D,d) > 0`,

preserves both the Singleton-type bound and the exact optimum dimension while
raising the minimum distance from `d` to `d + 1`.

Run the two independent finite checks from the repository root:

```bash
python3 verification/family/verify_family.py
python3 verification/family/bruteforce_tiny.py
```

`verify_family.py` starts from the published 11-dimensional certificate,
constructs the first eleven row cones, and enumerates all 2,048 codewords at
each stage. It checks support, every Singleton cut, and the shifted rank
distribution.

`bruteforce_tiny.py` does not use the published certificate. It enumerates
every binary subspace in its stated tiny range via unique RREF bases, checks
the Gaussian-binomial subspace counts, and independently confirms the exact
row-cone identity for every feasible diagram-distance pair in that range.

These programs audit the construction and small cases. The general identity
is proved symbolically in the paper.
