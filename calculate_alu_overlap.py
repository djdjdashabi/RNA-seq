#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path


BIN_SIZE = 100_000
CALLERS = ("reditools", "sprint", "jacusa2")
REQUIRED = {"sample_id", "caller", "contig", "start0", "end1", "a2i_candidate"}


class AluError(ValueError):
    pass


def normalize_contig(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return "MT" if value.upper() == "M" else value


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_alu(path: Path):
    bins = defaultdict(list)
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise AluError(f"Alu BED line {line_number}: expected at least 4 columns")
            try:
                start0, end1 = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise AluError(f"Alu BED line {line_number}: invalid coordinates") from exc
            # This BED is already filtered to RepeatMasker repClass=SINE and
            # repFamily=Alu.  Valid Alu monomers include FAM, FLAM_A, FLAM_C,
            # and FRAM, so the repeat name itself need not start with "Alu".
            if start0 < 0 or end1 <= start0 or not fields[3].strip():
                raise AluError(f"Alu BED line {line_number}: invalid Alu interval")
            contig = normalize_contig(fields[0])
            for bin_number in range(start0 // BIN_SIZE, (end1 - 1) // BIN_SIZE + 1):
                bins[(contig, bin_number)].append((start0, end1, fields[3]))
    return bins


def read_sites(paths: list[Path]):
    sites = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = sorted(REQUIRED - set(reader.fieldnames or []))
            if missing:
                raise AluError(f"{path}: canonical input missing columns: {', '.join(missing)}")
            for row in reader:
                marker = row["a2i_candidate"].strip().lower()
                if marker not in {"true", "false"}:
                    raise AluError(
                        f"{path}:{reader.line_num}: invalid a2i_candidate {marker!r}"
                    )
                if marker != "true":
                    continue
                try:
                    start0, end1 = int(row["start0"]), int(row["end1"])
                except ValueError as exc:
                    raise AluError(f"{path}:{reader.line_num}: invalid coordinate") from exc
                if start0 < 0 or end1 != start0 + 1:
                    raise AluError(f"{path}:{reader.line_num}: invalid single-base interval")
                key = (row["sample_id"], row["caller"], normalize_contig(row["contig"]), start0, end1)
                sites[key] = None
    return sites


def calculate(sites, bins, expected_pairs=None):
    expected_pairs = expected_pairs or set()
    memberships = []
    counts = defaultdict(lambda: [0, 0])
    for sample, caller, contig, start0, end1 in sorted(sites):
        matches = {
            name
            for alu_start, alu_end, name in bins.get((contig, start0 // BIN_SIZE), [])
            if alu_start < end1 and start0 < alu_end
        }
        counts[(sample, caller)][0] += 1
        if matches:
            counts[(sample, caller)][1] += 1
        memberships.append(
            {
                "sample_id": sample,
                "caller": caller,
                "contig": contig,
                "start0": str(start0),
                "end1": str(end1),
                "in_alu": "true" if matches else "false",
                "alu_names": ",".join(sorted(matches)),
            }
        )
    summary = []
    for sample, caller in sorted(set(counts) | expected_pairs):
        total, alu = counts[(sample, caller)]
        summary.append(
            {
                "sample_id": sample,
                "caller": caller,
                "scope": "a2i",
                "total_sites": str(total),
                "alu_sites": str(alu),
                "alu_fraction": format(alu / total if total else 0.0, ".12g"),
            }
        )
    return summary, memberships


def write_tsv(path: Path, fields: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--alu-bed", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--expected", action="append", default=[], metavar="SAMPLE:CALLER")
    return parser.parse_args(argv)


def parse_expected(values):
    pairs = set()
    for value in values:
        try:
            sample, caller = value.split(":", 1)
        except ValueError as exc:
            raise AluError(f"invalid --expected value: {value!r}") from exc
        if not sample or caller not in CALLERS:
            raise AluError(f"invalid --expected value: {value!r}")
        pairs.add((sample, caller))
    return pairs


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        summary, membership = calculate(
            read_sites(args.input), read_alu(args.alu_bed), parse_expected(args.expected)
        )
        write_tsv(args.summary, ["sample_id", "caller", "scope", "total_sites", "alu_sites", "alu_fraction"], summary)
        write_tsv(args.membership, ["sample_id", "caller", "contig", "start0", "end1", "in_alu", "alu_names"], membership)
    except (OSError, AluError) as exc:
        print(f"Alu overlap failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
