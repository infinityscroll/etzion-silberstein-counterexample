#!/usr/bin/env python3
"""Family questions around the E_6 counterexample:
1. ES bound of the transpose diagram F_6 = E_6^T = (6,4,4,4,4) at d=3
   (if 12, the E_6 result kills F_6 by transposition, free).
2. Do E_4=(3,3,1,1) (bound 2) and E_5=(4,4,4,1,1) (bound 6) attain their
   bounds over F_2? Randomized greedy search for explicit codes.
"""
import random

def es_bound(heights, d):
    # CN convention: nu_j = sum_{i=j+1..n} max(0, c_i - d + j + 1),
    # heights nonincreasing
    c = sorted(heights, reverse=True)
    return min(sum(max(0, ci - d + j + 1) for ci in c[j:]) for j in range(d))

print("nu_min(E_6=(5,5,5,5,1,1), 3) =", es_bound([5,5,5,5,1,1], 3))
print("nu_min(F_6=(6,4,4,4,4), 3)   =", es_bound([6,4,4,4,4], 3))
print("nu_min(E_5=(4,4,4,1,1), 3)   =", es_bound([4,4,4,1,1], 3))
print("nu_min(E_4=(3,3,1,1), 3)     =", es_bound([3,3,1,1], 3))

def cells_of(heights):
    return [(r, c) for c, h in enumerate(heights) for r in range(h)]

def rank_f2(colmasks, nrows, ncols):
    rows = []
    for r in range(nrows):
        row = 0
        for c in range(ncols):
            if colmasks[c] >> r & 1:
                row |= 1 << c
        rows.append(row)
    rank = 0
    for bit in range(ncols - 1, -1, -1):
        piv = next((i for i, row in enumerate(rows) if row >> bit & 1), None)
        if piv is None:
            continue
        pr = rows.pop(piv)
        rows = [x ^ pr if x >> bit & 1 else x for x in rows]
        rank += 1
    return rank

def search_code(heights, target_dim, d, tries=200000, seed=1):
    rng = random.Random(seed)
    cells = cells_of(heights)
    nrows, ncols = max(heights), len(heights)
    def rand_matrix():
        cm = [0] * ncols
        for (r, c) in cells:
            if rng.random() < 0.5:
                cm[c] |= 1 << r
        return tuple(cm)
    def ok(basis):
        n = len(basis)
        for coeffs in range(1, 1 << n):
            cm = [0] * ncols
            for i in range(n):
                if coeffs >> i & 1:
                    cm = [a ^ b for a, b in zip(cm, basis[i])]
            if rank_f2(cm, nrows, ncols) < d:
                return False
        return True
    basis = []
    for t in range(tries):
        cand = rand_matrix()
        if ok(basis + [cand]):
            basis.append(cand)
            if len(basis) == target_dim:
                return basis
        # occasional restart pruning
        if t % 20000 == 19999 and basis:
            basis.pop(rng.randrange(len(basis)))
    return None

for name, heights, k in (("E_4", [3,3,1,1], 2), ("E_5", [4,4,4,1,1], 6)):
    code = search_code(heights, k, 3, seed=42)
    if code:
        print(f"{name}: FOUND [{name},{k},3]_2 code — bound attained. Basis (column masks):")
        for g in code:
            print("   ", list(g))
    else:
        print(f"{name}: no dim-{k} code found by randomized search (NOT conclusive)")
