# Smaller members of the E_n family

The two JSON files give explicit optimal binary codes for
`E_4=(3,3,1,1)` and `E_5=(4,4,4,1,1)` at minimum rank distance three.
Run the exact standard-library verifier from the repository root:

```bash
python3 verification/small_examples/verify_small_examples.py
```

It checks the Ferrers support, Singleton cuts, linear independence, and rank
of every codeword.
