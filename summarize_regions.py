#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


REGIONS = (
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "noncoding_exon",
    "intron",
    "intergenic",
)
CALLERS = ("reditools", "sprint", "jacusa2")
REQUIRED = {"sample_id", "caller", "contig", "start0", "end1", "a2i_candidate", "region"}


class RegionSummaryError(ValueError):
    pass


def read_sites(paths: list[Path]):
    sites = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = sorted(REQUIRED - set(reader.fieldnames or []))
            if missing:
                raise RegionSummaryError(
                    f"{path}: annotated input missing columns: {', '.join(missing)}"
                )
            for row in reader:
                try:
                    start0, end1 = int(row["start0"]), int(row["end1"])
                except ValueError as exc:
                    raise RegionSummaryError(
                        f"{path}:{reader.line_num}: invalid coordinate"
                    ) from exc
                region = row["region"]
                if start0 < 0 or end1 <= start0 or region not in REGIONS:
                    raise RegionSummaryError(
                        f"{path}:{reader.line_num}: invalid interval or region {region!r}"
                    )
                marker = row["a2i_candidate"].lower()
                if marker not in {"true", "false"}:
                    raise RegionSummaryError(
                        f"{path}:{reader.line_num}: invalid a2i_candidate {marker!r}"
                    )
                key = (row["sample_id"], row["caller"], row["contig"], start0, end1)
                existing = sites.get(key)
                if existing and existing["region"] != region:
                    raise RegionSummaryError(
                        f"{path}:{reader.line_num}: duplicate site has conflicting regions"
                    )
                if existing:
                    existing["a2i"] = existing["a2i"] or marker == "true"
                else:
                    sites[key] = {"region": region, "a2i": marker == "true"}
    return sites


def summarize(sites, expected_pairs=None):
    expected_pairs = expected_pairs or set()
    grouped = defaultdict(list)
    for (sample, caller, _contig, _start0, _end1), value in sites.items():
        grouped[(sample, caller)].append(value)
    output = []
    for sample, caller in sorted(set(grouped) | expected_pairs):
        values = grouped[(sample, caller)]
        for scope in ("all", "a2i"):
            selected = values if scope == "all" else [value for value in values if value["a2i"]]
            counts = Counter(value["region"] for value in selected)
            total = len(selected)
            for region in REGIONS:
                count = counts[region]
                output.append(
                    {
                        "sample_id": sample,
                        "caller": caller,
                        "scope": scope,
                        "region": region,
                        "site_count": str(count),
                        "fraction": format(count / total if total else 0.0, ".12g"),
                    }
                )
    return output


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected", action="append", default=[], metavar="SAMPLE:CALLER")
    return parser.parse_args(argv)


def parse_expected(values):
    pairs = set()
    for value in values:
        try:
            sample, caller = value.split(":", 1)
        except ValueError as exc:
            raise RegionSummaryError(f"invalid --expected value: {value!r}") from exc
        if not sample or caller not in CALLERS:
            raise RegionSummaryError(f"invalid --expected value: {value!r}")
        pairs.add((sample, caller))
    return pairs


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        rows = summarize(read_sites(args.input), parse_expected(args.expected))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["sample_id", "caller", "scope", "region", "site_count", "fraction"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, RegionSummaryError) as exc:
        print(f"region summary failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
