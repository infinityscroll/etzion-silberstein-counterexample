# Mathematical audit for the binary E6 instance

## Reduction

Use the top-left-aligned diagram with column heights
`E6=(5,5,5,5,1,1)`.  A supported matrix has the block form

```text
[ r(A) | t(A) ]   (one row, widths 4 and 2)
[  A   |  0   ]   (four rows)
```

For distance three, the three Etzion--Silberstein deletion counts are all 12:
delete two top rows; delete the top row and first column; or delete the first
two columns.  Hence the conjectured maximum dimension is 12.

Suppose a 12-dimensional binary code `C` on this support has minimum rank
three.  Projection onto the bottom `4 x 4` block is injective, because an
element of its kernel has rank at most one.  Its image `U` therefore has
dimension 12.  Deleting the top row lowers rank by at most one, so every
nonzero element of `U` has rank at least two.  Thus `U` is a binary
`[4 x 4,12,2]` MRD code.

The tail map `t:U -> F_2^2` has a kernel `K`.  The corresponding tail-zero
subcode is a `5 x 4` rank-distance-three code, so the rank-metric Singleton
bound gives `dim K <= 10`.  Since `dim ker(t) >= 12-2=10`, the tail map is
surjective and `dim K=10`.

For `A in K`, minimum rank three says exactly

```text
rank([r(A); A]) >= 3.
```

Only rank-two `A` impose a condition, namely `r(A)` must not belong to the row
space of `A`.  Conversely, given a `[4 x 4,12,2]` code `U`, a codimension-two
subspace `K`, and a linear `r:K -> F_2^4` satisfying this condition, choose a
surjection `t:U -> F_2^2` with kernel `K` and extend `r` linearly to `U`.
If `t(A)` is nonzero, the top row is automatically independent of the bottom
four rows, so the resulting supported 12-space has minimum rank three.  This
proves both directions of the computational reduction.

## Exact lift clauses

For a rank-two matrix `A`, let `a,b` be a basis of the two-dimensional
orthogonal complement of its row space.  Then

```text
r(A) not in rowspace(A)
    iff (a dot r(A) = 1) OR (b dot r(A) = 1).
```

The checker represents `r` by an arbitrary extension `U -> F_2^4`, using 48
Boolean variables.  Such an extension always exists for every map on `K`.
Each displayed condition is therefore an OR of two linear equations.  The
solver maintains an exact affine row-echelon system and uses the disjoint,
exhaustive split

```text
first equation = 1
or
first equation = 0 and second equation = 1.
```

Before any research instance, 300 deterministic pseudorandom systems with at
most ten variables are compared against complete truth-table enumeration.

Every codimension-two kernel is represented once by the unique RREF basis of
its two-dimensional annihilator.  The exact number is

```text
[12 choose 2]_2 = (4096-1)(4096-2)/((4-1)(4-2)) = 2,794,155.
```

## Left, right, dual, and transpose invariance

For invertible `P,Q`, put `B=P A Q`, `K'=P K Q`, and
`r'(B)=r(A)Q`.  Then

```text
[r'(B); B] = diag(1,P) [r(A); A] Q,
```

so the lift condition and all ranks are preserved.  Hence feasibility is
invariant under oriented left/right equivalence, without using transpose.

For the Delsarte product `<A,B>=Tr(A B^T)`, if `S'=P S Q` then

```text
(S')^perp = P^(-T) S^perp Q^(-T).
```

Thus left/right classes of four-dimensional full-rank codes correspond to
left/right classes of their 12-dimensional dual MRD codes.

Transposition of `[r(A);A]` changes a top-row lift into a left-column lift and
transposes the Ferrers diagram: `E6^T=F6`.  This alone would not justify
collapsing oriented classes.  The independent enumeration instead finds three
left/right orbits already, and each is fixed by transpose.  Explicit
left/right maps from each named spread transpose back to itself were also
found.  Therefore no orientation is lost by checking the three named duals.

## Classification audit

The orthogonal of every binary `[4 x 4,12,2]` MRD code is a
`[4 x 4,4,4]` MRD code.  Here is a direct argument.  If a nonzero matrix `A`
in the four-dimensional orthogonal had rank at most three, choose nonzero `x`
with `x^T A=0`.  The four-space `{x y^T : y in F_2^4}` consists of rank-at-most
one matrices and lies, together with the 12-space `U`, in the 15-dimensional
hyperplane `A^perp`.  The two spaces must intersect nontrivially, giving a
rank-one word in `U`, a contradiction.  Thus the orthogonal consists of zero
together with 15 invertible matrices.  (The converse is also the standard
Delsarte duality statement.)

After normalizing a full-rank four-space to contain the identity, quotient by
`<I>`.  A nonzero quotient vector is admissible exactly when both members of
its coset are invertible.  Exhaustive enumeration gives 2,912 admissible
vectors and 19,936 normalized spreads.  Exhaustive left/right action partitions
them into exactly three orbits, of masses 16,800, 2,800, and 336.  The named
class-3, class-2, and field spreads lie in those three distinct orbits,
respectively, and transpose fixes every orbit.

This is an independent finite classification, rather than reliance on the
published MAGMA statement alone.
