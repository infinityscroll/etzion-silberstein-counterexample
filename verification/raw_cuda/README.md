# Raw CUDA audit

This checker enumerates all `2,794,155` kernels in each of the three MRD
classes, without orbit reduction.  It requires a CUDA-capable compiler and
C++20 because the host code uses `std::popcount` and `std::countl_zero`.

```bash
nvcc -O3 -std=c++20 cuda_all_kernels.cu -o cuda_all_kernels
nvcc -O3 -std=c++20 cuda_dpll_selftest.cu -o cuda_dpll_selftest

./cuda_dpll_selftest selftest.json
./cuda_all_kernels 0 0 2794155 class_field.json
./cuda_all_kernels 2 0 2794155 class_II.json
./cuda_all_kernels 3 0 2794155 class_III.json
```

The class numbers `0`, `2`, and `3` denote the field, II, and III classes,
respectively.  Status codes in the result files are `0 = UNSAT`, `1 = SAT`,
`2 = constraint overflow`, and `3 = solver overflow`.  The released complete
runs contain only `UNSAT` statuses.
