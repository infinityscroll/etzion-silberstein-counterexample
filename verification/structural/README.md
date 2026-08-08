# Four residual kernel instances

For every two-plane `N <= F_2^4`, a feasible lift forces
`dim(K intersect V_N)=2`.  Hence a feasible kernel must contain exactly
`35 * (2^2 - 1) = 105` rank-two matrices.  The exact orbit files in
`../primary` show that both nonfield MRD classes have larger minima and that
only four field-class orbits meet 105.  Their total raw orbit mass is 70.

Each `residual_*.cnf` is the exact lift system for one surviving orbit.  Both
CaDiCaL 3.0.1 and CryptoMiniSat 5.14.7 returned `UNSATISFIABLE` with exit code
20.  CaDiCaL's proof was independently checked by `drat-trim`; the reduced
proof was checked a second time.  The terminal logs contain `s VERIFIED`.

The reduced proofs are gzip-compressed.  Replay all four with the bundled
`drat-trim` source:

```bash
python3 verify_drat_proofs.py
```

Or replay one with an existing `drat-trim` executable:

```bash
gzip -dc residual_1_2.core.drat.gz | \
  drat-trim residual_1_2.cnf
```

The bundled checker source is from
[drat-trim](https://github.com/marijnheule/drat-trim) and retains its MIT
license in `tools/drat-trim/LICENSE`.

`generate_fixed_kernel_cnf.py` and its dependency `generate_cnf.py` regenerate
the instances.  Space-padded DIMACS headers are intentional: they remain
seekable during generation and are parsed consistently by independent proof
checkers.
