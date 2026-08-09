#!/usr/bin/env python3
"""EXACT settle of E_5 = (4,4,4,1,1), d=3, q=2: does a [E_5, 6, 3]_2 code exist?

Reduction (mirror of the E_6 paper's Lemma 1, dimensions checked):
[E_5,6,3]_2 exists  <=>  exists binary [3x3,6,2] MRD code U, codim-2 K <= U,
linear R: K -> F_2^3 with R(M) not in row(M) for every rank-2 M in K.

Over F_2 the dual of U is a [3x3,3,3] spread code <-> semifield of order 8,
and the only semifield of order 8 is F_8. So U is unique up to equivalence:
U = (F_8 multiplication matrices)^perp under the Frobenius form. We enumerate
ALL 651 codim-2 subspaces K and decide each affine parity system exactly.
"""
from itertools import combinations

# --- F_8 = F_2[x]/(x^3+x+1); multiplication matrices in basis 1,x,x^2 ---
def mul_f8(a, b):
    r = 0
    for i in range(3):
        if b >> i & 1:
            r ^= a << i
    for i in (5, 4, 3):
        if r >> i & 1:
            r ^= (0b1011 << (i - 3))
    return r & 7

def mult_matrix(a):
    # columns: images of basis vectors 1, x, x^2; encode matrix as 9-bit int,
    # bit 3*i+j = entry (row i, col j)
    m = 0
    for j, e in enumerate((1, 2, 4)):
        col = mul_f8(a, e)
        for i in range(3):
            if col >> i & 1:
                m |= 1 << (3 * i + j)
    return m

spread = [mult_matrix(a) for a in range(1, 8)]  # nonzero mult maps
# spread spans a 3-dim space; get basis
def span_basis(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return basis

S_basis = span_basis(spread)
assert len(S_basis) == 3

def frob(a, b):  # Frobenius form over F_2 = parity of AND
    return bin(a & b).count("1") & 1

# U = S^perp in 9-bit matrix space
U_basis = []
for v in range(1, 512):
    if all(frob(v, s) == 0 for s in S_basis):
        w = v
        for b in U_basis:
            w = min(w, w ^ b)
        if w:
            U_basis.append(w)
            U_basis.sort(reverse=True)
assert len(U_basis) == 6, len(U_basis)

def rank3x3(m):
    rows = [(m >> (3 * i)) & 7 for i in range(3)]
    rank = 0
    for bit in (2, 1, 0):
        piv = next((i for i, r in enumerate(rows) if r >> bit & 1), None)
        if piv is None:
            continue
        pr = rows.pop(piv)
        rows = [r ^ pr if r >> bit & 1 else r for r in rows]
        rank += 1
    return rank

# sanity: U is [3x3,6,2] MRD — no rank-1 nonzero elements
def elems(basis):
    n = len(basis)
    for c in range(1, 1 << n):
        v = 0
        for i in range(n):
            if c >> i & 1:
                v ^= basis[i]
        yield c, v

assert all(rank3x3(v) >= 2 for _, v in elems(U_basis))

def right_kernel_vec(m):
    # unique nonzero n in F_2^3 with M n = 0 for rank-2 m
    rows = [(m >> (3 * i)) & 7 for i in range(3)]
    for n in range(1, 8):
        if all(bin(r & n).count("1") % 2 == 0 for r in rows):
            return n
    return None

# enumerate all 651 codim-2 subspaces K of U via 2-dim subspaces of the dual
# coordinate space F_2^6 (coordinates w.r.t. U_basis)
def coords_iter():
    seen = set()
    for a in range(1, 64):
        for b in range(a + 1, 64):
            key = frozenset({a, b, a ^ b})
            if key in seen:
                continue
            seen.add(key)
            yield (a, b)

feasible_kernels = 0
total = 0
witness = None
for (a, b) in coords_iter():
    total += 1
    # K = {c in F_2^6 : parity(c&a)=0, parity(c&b)=0}
    kbasis = []
    for c in range(1, 64):
        if bin(c & a).count("1") % 2 == 0 and bin(c & b).count("1") % 2 == 0:
            w = c
            for bb in kbasis:
                w = min(w, w ^ bb)
            if w:
                kbasis.append(w)
                kbasis.sort(reverse=True)
    assert len(kbasis) == 4
    # rank-2 elements of K and their kernels
    constraints = []  # (coeff vector over 12 unknowns, rhs=1)
    okK = True
    for c in range(1, 16):
        cm = 0
        cvec = 0
        for i in range(4):
            if c >> i & 1:
                cvec ^= kbasis[i]
        for i in range(6):
            if cvec >> i & 1:
                cm ^= U_basis[i]
        r = rank3x3(cm)
        if r == 2:
            n = right_kernel_vec(cm)
            # R(M) = sum over basis elts in combination c of R_i; unknowns:
            # R_i in F_2^3 for i=0..3 -> 12 bits. Condition: R(M).n == 1.
            # coefficient of bit (3*i + j) is [c_i] * n_j
            coeff = 0
            for i in range(4):
                if c >> i & 1:
                    for j in range(3):
                        if n >> j & 1:
                            coeff |= 1 << (3 * i + j)
            constraints.append(coeff)
    # solve affine system: coeff . x = 1 for all constraints (over F_2)
    # Gaussian elimination on (coeff | 1) rows
    rows = [(cf << 1) | 1 for cf in constraints]  # bit0 = rhs
    consistent = True
    pivots = []
    for row in rows:
        for p in pivots:
            if row & p & ~1:
                top_row = row >> 1
                top_p = p >> 1
                # reduce by matching leading bit
        # simpler full reduction
    # redo reduction cleanly
    basis_rows = []
    for row in rows:
        cur = row
        changed = True
        while changed:
            changed = False
            for p in basis_rows:
                if (cur >> 1) and (p >> 1) and ((cur >> 1).bit_length() == (p >> 1).bit_length()):
                    cur ^= p
                    changed = True
        if cur >> 1:
            basis_rows.append(cur)
            basis_rows.sort(key=lambda r: (r >> 1).bit_length(), reverse=True)
        elif cur & 1:
            consistent = False
            break
    if consistent:
        feasible_kernels += 1
        if witness is None:
            witness = (a, b, len(constraints))

print(f"kernels enumerated: {total} (expect 651)")
print(f"kernels with consistent lift system: {feasible_kernels}")
if feasible_kernels:
    print(f"E_5 ATTAINS its bound over F_2: [E_5,6,3]_2 exists "
          f"(first witness kernel annihilator {witness[:2]}, "
          f"{witness[2]} rank-2 constraints)")
else:
    print("E_5 FAILS over F_2 — smaller counterexample than E_6!")
