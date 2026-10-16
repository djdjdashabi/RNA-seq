#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import itertools
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


CALLER_ORDER = ("reditools", "sprint", "jacusa2")
EDITING_LOG10_MIN = -5.0
EDITING_LOG10_MAX = 0.0
EDITING_LOG10_BINS = 50
REQUIRED_COLUMNS = {
    "sample_id",
    "caller",
    "contig",
    "start0",
    "end1",
    "ref",
    "alt",
    "substitution",
    "a2i_candidate",
    "editing_level",
}


class SummaryError(ValueError):
    pass


def compact_number(value: float) -> str:
    return format(value, ".12g")


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows_by_key = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames or []))
            if missing:
                raise SummaryError(
                    f"{path}: standardized input missing columns: {', '.join(missing)}"
                )
            for row in reader:
                try:
                    start0 = int(row["start0"])
                    end1 = int(row["end1"])
                    editing_level = float(row["editing_level"])
                except ValueError as exc:
                    raise SummaryError(
                        f"{path}:{reader.line_num}: invalid standardized numeric value"
                    ) from exc
                marker = row["a2i_candidate"].strip().lower()
                if (
                    not row["sample_id"].strip()
                    or row["caller"] not in CALLER_ORDER
                    or not row["contig"].strip()
                    or start0 < 0
                    or end1 != start0 + 1
                    or marker not in {"true", "false"}
                    or not math.isfinite(editing_level)
                    or not 0 <= editing_level <= 1
                    or row["substitution"] != f"{row['ref']}>{row['alt']}"
                ):
                    raise SummaryError(
                        f"{path}:{reader.line_num}: invalid standardized single-base record"
                    )
                key = (
                    row["sample_id"],
                    row["caller"],
                    row["contig"],
                    start0,
                    end1,
                    row["ref"],
                    row["alt"],
                )
                normalized = dict(row)
                normalized["start0"] = str(start0)
                normalized["end1"] = str(end1)
                normalized["a2i_candidate"] = marker
                normalized["editing_level"] = compact_number(editing_level)
                existing = rows_by_key.get(key)
                if existing is not None and any(
                    existing[field] != normalized[field]
                    for field in ("substitution", "a2i_candidate", "editing_level")
                ):
                    raise SummaryError(
                        f"{path}:{reader.line_num}: duplicate allele has conflicting values"
                    )
                rows_by_key.setdefault(key, normalized)
    return list(rows_by_key.values())


def caller_sort_key(caller: str) -> tuple[int, str]:
    try:
        return CALLER_ORDER.index(caller), caller
    except ValueError:
        return len(CALLER_ORDER), caller


def scope_rows(rows: list[dict[str, str]], scope: str) -> list[dict[str, str]]:
    if scope == "all":
        return rows
    return [row for row in rows if row["a2i_candidate"].lower() == "true"]


def site_key(row: dict[str, str]) -> tuple[str, int, int]:
    return row["contig"], int(row["start0"]), int(row["end1"])


def editing_level_histogram(levels: list[float]):
    width = (EDITING_LOG10_MAX - EDITING_LOG10_MIN) / EDITING_LOG10_BINS
    counts = [0] * EDITING_LOG10_BINS
    for level in levels:
        if not 0 <= level <= 1:
            raise SummaryError(f"editing level outside [0, 1]: {level}")
        log_level = math.log10(max(level, 10 ** EDITING_LOG10_MIN))
        index = min(
            EDITING_LOG10_BINS - 1,
            max(0, int((log_level - EDITING_LOG10_MIN) / width)),
        )
        counts[index] += 1
    denominator = len(levels)
    return [
        (
            EDITING_LOG10_MIN + (index + 0.5) * width,
            count / denominator if denominator else 0.0,
        )
        for index, count in enumerate(counts)
    ]


def build_summary(
    rows: list[dict[str, str]], expected_pairs: set[tuple[str, str]] | None = None
) -> list[dict[str, str]]:
    expected_pairs = expected_pairs or set()
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample_id"], row["caller"])].append(row)
    output = []
    for sample, caller in sorted(
        set(grouped) | expected_pairs,
        key=lambda key: (key[0], caller_sort_key(key[1])),
    ):
        native_rows = grouped[(sample, caller)]
        for scope in ("all", "a2i"):
            selected = scope_rows(native_rows, scope)
            sites = {site_key(row) for row in selected}
            levels = [float(row["editing_level"]) for row in selected]
            scalars = [
                ("site_count", len(sites)),
                ("allele_count", len(selected)),
                ("editing_level_mean", statistics.fmean(levels) if levels else 0.0),
                ("editing_level_median", statistics.median(levels) if levels else 0.0),
            ]
            for metric, value in scalars:
                output.append(
                    {
                        "sample_id": sample,
                        "caller": caller,
                        "scope": scope,
                        "metric": metric,
                        "group": "",
                        "value": str(value) if isinstance(value, int) else compact_number(value),
                    }
                )
            if scope == "a2i":
                for bin_center, fraction in editing_level_histogram(levels):
                    output.append(
                        {
                            "sample_id": sample,
                            "caller": caller,
                            "scope": scope,
                            "metric": "editing_level_log10_fraction",
                            "group": compact_number(bin_center),
                            "value": compact_number(fraction),
                        }
                    )
            if scope == "all" and selected:
                counts = Counter(row["substitution"] for row in selected)
                total = sum(counts.values())
                for substitution in sorted(counts):
                    output.append(
                        {
                            "sample_id": sample,
                            "caller": caller,
                            "scope": scope,
                            "metric": "substitution_fraction",
                            "group": substitution,
                            "value": compact_number(counts[substitution] / total),
                        }
                    )
    return output


def build_overlap(
    rows: list[dict[str, str]], expected_pairs: set[tuple[str, str]] | None = None
) -> list[dict[str, str]]:
    expected_pairs = expected_pairs or set()
    samples = sorted({row["sample_id"] for row in rows} | {pair[0] for pair in expected_pairs})
    output = []
    for sample in samples:
        sample_rows = [row for row in rows if row["sample_id"] == sample]
        callers = sorted(
            {row["caller"] for row in sample_rows}
            | {caller for expected_sample, caller in expected_pairs if expected_sample == sample},
            key=caller_sort_key,
        )
        for scope in ("all", "a2i"):
            sites = {
                caller: {
                    site_key(row)
                    for row in scope_rows(sample_rows, scope)
                    if row["caller"] == caller
                }
                for caller in callers
            }
            for caller_a, caller_b in itertools.combinations(callers, 2):
                left, right = sites[caller_a], sites[caller_b]
                intersection = len(left & right)
                union = len(left | right)
                output.append(
                    {
                        "sample_id": sample,
                        "caller_a": caller_a,
                        "caller_b": caller_b,
                        "scope": scope,
                        "n_a": str(len(left)),
                        "n_b": str(len(right)),
                        "intersection": str(intersection),
                        "union": str(union),
                        "jaccard": compact_number(intersection / union if union else 0.0),
                    }
                )
    return output


def build_support(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    samples = sorted({row["sample_id"] for row in rows})
    output = []
    for sample in samples:
        sample_rows = [row for row in rows if row["sample_id"] == sample]
        for scope in ("all", "a2i"):
            support = defaultdict(set)
            for row in scope_rows(sample_rows, scope):
                support[site_key(row)].add(row["caller"])
            for (contig, start0, end1), callers in sorted(
                support.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])
            ):
                ordered_callers = sorted(callers, key=caller_sort_key)
                output.append(
                    {
                        "sample_id": sample,
                        "scope": scope,
                        "contig": contig,
                        "start0": str(start0),
                        "end1": str(end1),
                        "callers": ",".join(ordered_callers),
                        "support_count": str(len(ordered_callers)),
                    }
                )
    return output


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--overlap", required=True, type=Path)
    parser.add_argument("--support", required=True, type=Path)
    parser.add_argument(
        "--expected",
        action="append",
        default=[],
        metavar="SAMPLE:CALLER",
        help="Expected sample/caller pair; emits explicit zero summaries when empty",
    )
    return parser.parse_args(argv)


def parse_expected(values: list[str]) -> set[tuple[str, str]]:
    pairs = set()
    for value in values:
        try:
            sample, caller = value.split(":", 1)
        except ValueError as exc:
            raise SummaryError(f"invalid --expected value: {value!r}") from exc
        if not sample or caller not in CALLER_ORDER:
            raise SummaryError(f"invalid --expected value: {value!r}")
        pairs.add((sample, caller))
    return pairs


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        rows = read_rows(args.input)
        expected_pairs = parse_expected(args.expected)
        write_tsv(
            args.summary,
            ["sample_id", "caller", "scope", "metric", "group", "value"],
            build_summary(rows, expected_pairs),
        )
        write_tsv(
            args.overlap,
            ["sample_id", "caller_a", "caller_b", "scope", "n_a", "n_b", "intersection", "union", "jaccard"],
            build_overlap(rows, expected_pairs),
        )
        write_tsv(
            args.support,
            ["sample_id", "scope", "contig", "start0", "end1", "callers", "support_count"],
            build_support(rows),
        )
    except (OSError, SummaryError) as exc:
        print(f"summary failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
