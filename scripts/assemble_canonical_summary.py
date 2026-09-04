#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


FIELDS = ["sample_id", "caller", "scope", "metric", "group", "value", "significant", "source"]
CALLERS = ("reditools", "sprint", "jacusa2")


class AssemblyError(ValueError):
    pass


def read_rows(path: Path, required: set[str], label: str):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise AssemblyError(f"{label} missing columns: {', '.join(missing)}")
        yield from reader


def emit(sample, caller, scope, metric, group, value, source, significant=""):
    return {
        "sample_id": sample,
        "caller": caller,
        "scope": scope,
        "metric": metric,
        "group": group,
        "value": value,
        "significant": significant,
        "source": source,
    }


def assemble(args):
    output = []
    expected_pairs = set()
    for row in read_rows(
        args.candidate_summary,
        {"sample_id", "caller", "scope", "metric", "group", "value"},
        "candidate summary",
    ):
        expected_pairs.add((row["sample_id"], row["caller"]))
        output.append(emit(row["sample_id"], row["caller"], row["scope"], row["metric"], row["group"], row["value"], "candidate_summary"))

    for row in read_rows(
        args.overlap,
        {"sample_id", "caller_a", "caller_b", "scope", "intersection", "union", "jaccard"},
        "overlap summary",
    ):
        pair = f"{row['caller_a']}|{row['caller_b']}"
        output.extend(
            [
                emit(row["sample_id"], pair, row["scope"], "pairwise_intersection", "", row["intersection"], "overlap"),
                emit(row["sample_id"], pair, row["scope"], "pairwise_union", "", row["union"], "overlap"),
                emit(row["sample_id"], pair, row["scope"], "pairwise_jaccard", "", row["jaccard"], "overlap"),
            ]
        )

    combinations = Counter()
    upset_groups = (
        ("reditools",),
        ("sprint",),
        ("jacusa2",),
        ("reditools", "sprint"),
        ("reditools", "jacusa2"),
        ("sprint", "jacusa2"),
        CALLERS,
    )
    for sample in {sample for sample, _caller in expected_pairs}:
        for scope in ("all", "a2i"):
            for callers in upset_groups:
                combinations[(sample, scope, "&".join(callers))] = 0
    for row in read_rows(args.support, {"sample_id", "scope", "callers"}, "site support"):
        raw_callers = list(filter(None, row["callers"].split(",")))
        unknown = sorted(set(raw_callers) - set(CALLERS))
        if unknown:
            raise AssemblyError(
                "site support contains unknown callers: " + ", ".join(unknown)
            )
        callers = sorted(raw_callers, key=lambda caller: CALLERS.index(caller))
        combination = "&".join(callers)
        combinations[(row["sample_id"], row["scope"], combination)] += 1
    for (sample, scope, combination), count in sorted(combinations.items()):
        output.append(emit(sample, "", scope, "upset_site_count", combination, str(count), "site_support"))

    for row in read_rows(
        args.alu,
        {"sample_id", "caller", "scope", "total_sites", "alu_sites", "alu_fraction"},
        "Alu summary",
    ):
        for metric in ("total_sites", "alu_sites", "alu_fraction"):
            output.append(emit(row["sample_id"], row["caller"], row["scope"], metric, "Alu", row[metric], "alu_summary"))

    for row in read_rows(
        args.regions,
        {"sample_id", "caller", "scope", "region", "site_count", "fraction"},
        "region summary",
    ):
        output.extend(
            [
                emit(
                    row["sample_id"], row["caller"], row["scope"],
                    "region_site_count", row["region"], row["site_count"], "region_summary",
                ),
                emit(
                    row["sample_id"], row["caller"], row["scope"],
                    "region_fraction", row["region"], row["fraction"], "region_summary",
                ),
            ]
        )

    for go_path in args.go:
        for row in read_rows(
            go_path,
            {
                "sample_id", "caller", "strand_mode", "go_id", "go_term", "category",
                "overlap_size", "selected_size", "fdr_bh", "significant",
            },
            "GO enrichment",
        ):
            try:
                overlap_size = int(row["overlap_size"])
                selected_size = int(row["selected_size"])
            except ValueError as exc:
                raise AssemblyError("GO enrichment contains a non-integer gene count") from exc
            group = f"{row['go_id']}|{row['category']}|{row['go_term']}"
            scope = f"a2i:{row['strand_mode']}"
            common = (
                row["sample_id"], row["caller"], scope, group,
                "go_enrichment", row["significant"],
            )
            output.extend(
                [
                    emit(common[0], common[1], common[2], "go_fdr_bh", common[3], row["fdr_bh"], common[4], common[5]),
                    emit(common[0], common[1], common[2], "go_overlap_size", common[3], str(overlap_size), common[4], common[5]),
                    emit(
                        common[0], common[1], common[2], "go_gene_ratio", common[3],
                        format(overlap_size / selected_size if selected_size else 0.0, ".12g"),
                        common[4], common[5],
                    ),
                ]
            )
    return output


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-summary", type=Path, required=True)
    parser.add_argument("--overlap", type=Path, required=True)
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--alu", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--go", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        rows = assemble(args)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, AssemblyError) as exc:
        print(f"canonical summary assembly failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
