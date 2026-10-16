#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        fieldnames = None
        rows = []
        for path in args.input:
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if not reader.fieldnames:
                    raise ValueError(f"{path}: missing header")
                if fieldnames is None:
                    fieldnames = reader.fieldnames
                elif reader.fieldnames != fieldnames:
                    raise ValueError(f"{path}: header differs from the first input")
                rows.extend(reader)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, ValueError) as exc:
        print(f"TSV concatenation failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

