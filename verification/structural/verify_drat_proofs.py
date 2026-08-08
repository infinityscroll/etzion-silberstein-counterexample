#!/usr/bin/env python3
"""Replay the four compressed DRAT certificates with the bundled checker."""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path


STEMS = (
    "residual_1_2",
    "residual_2_5",
    "residual_256_512",
    "residual_256_3072",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--checker",
        type=Path,
        help="use an existing drat-trim executable instead of compiling the bundled source",
    )
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="e6-drat-") as temporary:
        work = Path(temporary)
        checker = args.checker
        if checker is None:
            checker = work / "drat-trim"
            subprocess.run(
                [
                    "cc",
                    "-O2",
                    "-DNDEBUG",
                    str(here / "tools" / "drat-trim" / "drat-trim.c"),
                    "-o",
                    str(checker),
                ],
                check=True,
            )
        checker = checker.resolve()

        for stem in STEMS:
            compressed = here / f"{stem}.core.drat.gz"
            proof = work / f"{stem}.core.drat"
            with gzip.open(compressed, "rb") as source, proof.open("wb") as target:
                shutil.copyfileobj(source, target)
            result = subprocess.run(
                [str(checker), str(here / f"{stem}.cnf"), str(proof)],
                check=False,
                capture_output=True,
                text=True,
            )
            transcript = result.stdout + result.stderr
            if result.returncode != 0 or "s VERIFIED" not in transcript:
                raise SystemExit(
                    f"{stem}: proof replay failed (exit {result.returncode})\n{transcript}"
                )
            print(f"{stem}: VERIFIED")


if __name__ == "__main__":
    main()
