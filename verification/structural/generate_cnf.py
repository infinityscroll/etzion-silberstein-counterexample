#!/usr/bin/env python3
"""Generate compact SAT instances for a binary MFD code on E_6.

The support has column heights (5,5,5,5,1,1).  A systematic basis has ten
generators with zero tail and two generators with tails 01 and 10.  For a
zero-tail codeword its 5x4 core must have rank at least three.  For a nonzero
tail, adjoining the fixed column e_1 raises the rank by one, so it is enough
and necessary that the bottom 4x4 projection of the core have rank at least
two.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path


ROWS = 5
CORE_COLUMNS = 4
DIMENSION = 12
KERNEL_DIMENSION = 10


class CNF:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.file = path.open("w")
        self.variables = 0
        self.clauses = 0
        # Keep the fixed-width header seekable, but pad with spaces rather
        # than zeroes.  Some independent DRAT checkers interpret very long
        # zero-padded DIMACS counts incorrectly even though SAT solvers accept
        # them.
        self.header = "p cnf {variables:20d} {clauses:20d}\n"
        self.file.write(self.header.format(variables=0, clauses=0))

    def var(self) -> int:
        self.variables += 1
        return self.variables

    def clause(self, *literals: int) -> None:
        assert literals and all(literal != 0 for literal in literals)
        self.file.write(" ".join(map(str, literals)) + " 0\n")
        self.clauses += 1

    def finish(self) -> None:
        self.file.flush()
        self.file.seek(0)
        self.file.write(self.header.format(variables=self.variables, clauses=self.clauses))
        self.file.close()


def encode_equiv(cnf: CNF, output: int, source: int) -> None:
    cnf.clause(-output, source)
    cnf.clause(output, -source)


def encode_xor(cnf: CNF, output: int, left: int, right: int) -> None:
    cnf.clause(-left, -right, -output)
    cnf.clause(-left, right, output)
    cnf.clause(left, -right, output)
    cnf.clause(left, right, -output)


def encode_and(cnf: CNF, output: int, inputs: tuple[int, ...]) -> None:
    for source in inputs:
        cnf.clause(-output, source)
    cnf.clause(output, *(-source for source in inputs))


def encode_parity(cnf: CNF, inputs: list[int]) -> int:
    assert inputs
    if len(inputs) == 1:
        output = cnf.var()
        encode_equiv(cnf, output, inputs[0])
        return output
    current = inputs[0]
    for source in inputs[1:]:
        output = cnf.var()
        encode_xor(cnf, output, current, source)
        current = output
    return current


def encode_determinant(cnf: CNF, entries: list[list[int]]) -> int:
    size = len(entries)
    assert size in (2, 3) and all(len(row) == size for row in entries)
    products: list[int] = []
    for permutation in itertools.permutations(range(size)):
        product = cnf.var()
        encode_and(cnf, product, tuple(entries[row][permutation[row]] for row in range(size)))
        products.append(product)
    # Signs disappear in characteristic two, so determinant is the XOR of
    # all permutation monomials.
    return encode_parity(cnf, products)


def canonical_core(orbit: str) -> int:
    columns = {
        "rank3_contains_e1": (1, 2, 4, 0),
        "rank3_avoids_e1": (2, 4, 8, 0),
        "rank4_contains_e1": (1, 2, 4, 8),
        "rank4_avoids_e1": (2, 4, 8, 16),
    }[orbit]
    return sum(column << (ROWS * index) for index, column in enumerate(columns))


def generate(output: Path, orbit: str) -> dict[str, object]:
    cnf = CNF(output)
    generators = [[cnf.var() for _ in range(ROWS * CORE_COLUMNS)] for _ in range(DIMENSION)]

    fixed: int | None = None
    if orbit == "systematic":
        # Projection of the 10-dimensional zero-tail kernel to any two core
        # columns is injective: a word vanishing there would be supported on
        # only two columns and have rank at most two.  It is therefore an
        # isomorphism.  Choose the kernel basis to make the first two columns
        # the 10x10 identity.  Adding kernel words to either tail generator
        # then makes its first two columns zero.  This loses no solutions.
        for generator in range(KERNEL_DIMENSION):
            for bit in range(2 * ROWS):
                variable = generators[generator][bit]
                cnf.clause(variable if bit == generator else -variable)
        for generator in range(KERNEL_DIMENSION, DIMENSION):
            for bit in range(2 * ROWS):
                cnf.clause(-generators[generator][bit])
    else:
        fixed = canonical_core(orbit)
        for bit, variable in enumerate(generators[0]):
            cnf.clause(variable if ((fixed >> bit) & 1) else -variable)

    word_variables: dict[int, list[int]] = {}
    determinant_counts = {"3x3": 0, "2x2": 0}
    for coefficients in range(1, 1 << DIMENSION):
        selected = [index for index in range(DIMENSION) if (coefficients >> index) & 1]
        word = [
            encode_parity(cnf, [generators[index][bit] for index in selected])
            for bit in range(ROWS * CORE_COLUMNS)
        ]
        word_variables[coefficients] = word

        tail = ((coefficients >> KERNEL_DIMENSION) & 0b11)
        determinants: list[int] = []
        if tail == 0:
            for row_set in itertools.combinations(range(ROWS), 3):
                for column_set in itertools.combinations(range(CORE_COLUMNS), 3):
                    entries = [
                        [word[column * ROWS + row] for column in column_set]
                        for row in row_set
                    ]
                    determinants.append(encode_determinant(cnf, entries))
                    determinant_counts["3x3"] += 1
        else:
            # Quotient by the span of e_1: retain rows 1,...,4.
            for row_set in itertools.combinations(range(1, ROWS), 2):
                for column_set in itertools.combinations(range(CORE_COLUMNS), 2):
                    entries = [
                        [word[column * ROWS + row] for column in column_set]
                        for row in row_set
                    ]
                    determinants.append(encode_determinant(cnf, entries))
                    determinant_counts["2x2"] += 1
        cnf.clause(*determinants)

    variable_count = cnf.variables
    clause_count = cnf.clauses
    cnf.finish()
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    metadata = {
        "status": "GENERATED",
        "field": "F2",
        "support_column_heights": [5, 5, 5, 5, 1, 1],
        "dimension": DIMENSION,
        "minimum_rank": 3,
        "orbit": orbit,
        "canonical_first_generator_core": fixed,
        "variables": variable_count,
        "clauses": clause_count,
        "determinants": determinant_counts,
        "generator_variables": generators,
        "word_variables": word_variables,
        "cnf_sha256": digest,
    }
    output.with_suffix(".map.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbit", required=True, choices=(
        "systematic",
        "rank3_contains_e1", "rank3_avoids_e1",
        "rank4_contains_e1", "rank4_avoids_e1",
    ))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata = generate(args.output, args.orbit)
    print(json.dumps({key: metadata[key] for key in (
        "orbit", "variables", "clauses", "determinants", "cnf_sha256"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
