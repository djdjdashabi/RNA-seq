#!/usr/bin/env python3

from __future__ import annotations

import argparse
import ast
import csv
import re
import sys
from pathlib import Path


COLUMNS = [
    "sample_id",
    "caller",
    "contig",
    "start0",
    "end1",
    "pos1",
    "native_strand",
    "ref",
    "alt",
    "substitution",
    "canonical_ref",
    "canonical_alt",
    "a2i_candidate",
    "strand_consistent",
    "coverage",
    "alt_count",
    "editing_level",
    "native_source",
    "native_line",
]

BASES = ("A", "C", "G", "T")


class StandardizationError(ValueError):
    """Raised when native output does not satisfy its documented schema."""


def normalize_contig(value: str) -> str:
    contig = value.strip()
    if contig.lower().startswith("chr"):
        contig = contig[3:]
    if contig.upper() == "M":
        return "MT"
    if not contig:
        raise StandardizationError("empty contig")
    return contig


def parse_base_counts(value: str, line_number: int) -> dict[str, int]:
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError) as exc:
        raise StandardizationError(
            f"line {line_number}: invalid BaseCount[A,C,G,T]: {value!r}"
        ) from exc
    if not isinstance(parsed, (list, tuple)) or len(parsed) != 4:
        raise StandardizationError(
            f"line {line_number}: BaseCount[A,C,G,T] must contain four counts"
        )
    try:
        counts = [int(item) for item in parsed]
    except (TypeError, ValueError) as exc:
        raise StandardizationError(
            f"line {line_number}: base counts must be integers"
        ) from exc
    if any(count < 0 for count in counts):
        raise StandardizationError(f"line {line_number}: negative base count")
    return dict(zip(BASES, counts))


def parse_substitutions(value: str, ref: str, line_number: int) -> list[str]:
    if value.strip() in {"", "-", "."}:
        return []
    matches = re.findall(r"([ACGT])(?:>|2)?([ACGT])", value.upper())
    substitutions = []
    for native_ref, alt in matches:
        if native_ref != ref:
            raise StandardizationError(
                f"line {line_number}: AllSubs reference {native_ref} disagrees with {ref}"
            )
        if alt != ref and alt not in substitutions:
            substitutions.append(alt)
    if not substitutions:
        raise StandardizationError(
            f"line {line_number}: cannot parse AllSubs value {value!r}"
        )
    return substitutions


def canonical_orientation(ref: str, alt: str) -> tuple[str, str, bool]:
    if (ref, alt) == ("T", "C"):
        return "A", "G", True
    return ref, alt, (ref, alt) == ("A", "G")


def standard_row(
    *,
    sample: str,
    caller: str,
    contig: str,
    start0: int,
    end1: int,
    strand: str,
    ref: str,
    alt: str,
    coverage: int,
    alt_count: int,
    source: str,
    native_line: int,
) -> dict[str, object]:
    canonical_ref, canonical_alt, is_a2i = canonical_orientation(ref, alt)
    level = alt_count / coverage if coverage else 0.0
    return {
        "sample_id": sample,
        "caller": caller,
        "contig": contig,
        "start0": start0,
        "end1": end1,
        "pos1": end1,
        "native_strand": strand,
        "ref": ref,
        "alt": alt,
        "substitution": f"{ref}>{alt}",
        "canonical_ref": canonical_ref,
        "canonical_alt": canonical_alt,
        "a2i_candidate": "true" if is_a2i else "false",
        "strand_consistent": "",
        "coverage": coverage,
        "alt_count": alt_count,
        "editing_level": format(level, ".12g"),
        "native_source": source,
        "native_line": native_line,
    }


def parse_reditools(path: Path, sample: str) -> list[dict[str, object]]:
    required = {
        "Region",
        "Position",
        "Reference",
        "Strand",
        "Coverage",
        "BaseCount[A,C,G,T]",
        "AllSubs",
    }
    output = []
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        present = set(reader.fieldnames or [])
        missing = sorted(required - present)
        if missing:
            raise StandardizationError(
                "REDItools output missing required columns: " + ", ".join(missing)
            )
        for record in reader:
            line_number = reader.line_num
            try:
                contig = normalize_contig(record["Region"])
                pos1 = int(record["Position"])
                ref = record["Reference"].strip().upper()
                strand = record["Strand"].strip() or "."
                coverage = int(record["Coverage"])
            except (TypeError, ValueError) as exc:
                raise StandardizationError(
                    f"line {line_number}: invalid REDItools coordinate or coverage"
                ) from exc
            if pos1 < 1 or ref not in BASES or coverage < 0:
                raise StandardizationError(f"line {line_number}: invalid REDItools values")
            counts = parse_base_counts(record["BaseCount[A,C,G,T]"], line_number)
            if sum(counts.values()) != coverage:
                raise StandardizationError(
                    f"line {line_number}: base counts do not sum to coverage"
                )
            for alt in parse_substitutions(record["AllSubs"], ref, line_number):
                key = (sample, "reditools", contig, pos1 - 1, pos1, ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                output.append(
                    standard_row(
                        sample=sample,
                        caller="reditools",
                        contig=contig,
                        start0=pos1 - 1,
                        end1=pos1,
                        strand=strand,
                        ref=ref,
                        alt=alt,
                        coverage=coverage,
                        alt_count=counts[alt],
                        source=path.name,
                        native_line=line_number,
                    )
                )
    return output


def parse_sprint(path: Path, sample: str) -> list[dict[str, object]]:
    expected_header = [
        "#Chrom", "Start(0base)", "End(1base)", "Type",
        "Supporting_reads", "Strand", "AD:DP",
    ]
    output = []
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        records = [(line_number, fields) for line_number, fields in enumerate(reader, 1) if any(fields)]
        if not records:
            raise StandardizationError("SPRINT output is empty")
        has_header = records[0][1] == expected_header
        start_index = 1 if has_header else 0
        if records[0][1][0].startswith("#") and not has_header:
            raise StandardizationError("SPRINT output has an unknown header schema")
        for line_number, fields in records[start_index:]:
            if len(fields) != 7:
                raise StandardizationError(
                    f"line {line_number}: SPRINT regular output must contain 7 columns"
                )
            record = dict(zip(expected_header, fields))
            try:
                contig = normalize_contig(record["#Chrom"])
                start0 = int(record["Start(0base)"])
                end1 = int(record["End(1base)"])
                strand = record["Strand"].strip() or "."
            except (TypeError, ValueError) as exc:
                raise StandardizationError(
                    f"line {line_number}: invalid SPRINT coordinate"
                ) from exc
            if start0 < 0 or end1 != start0 + 1:
                raise StandardizationError(
                    f"line {line_number}: SPRINT site must be a single-base BED interval"
                )
            mutation = record["Type"].strip().upper().replace(">", "").replace("2", "")
            if not re.fullmatch(r"[ACGT]{2}", mutation) or mutation[0] == mutation[1]:
                raise StandardizationError(
                    f"line {line_number}: invalid SPRINT substitution {record['Type']!r}"
                )
            depth_match = re.fullmatch(r"(\d+):(\d+)", record["AD:DP"].strip())
            if not depth_match:
                raise StandardizationError(
                    f"line {line_number}: invalid SPRINT AD:DP value {record['AD:DP']!r}"
                )
            alt_count, coverage = (int(value) for value in depth_match.groups())
            if coverage < 1 or alt_count > coverage:
                raise StandardizationError(f"line {line_number}: invalid SPRINT depth")
            ref, alt = mutation
            key = (sample, "sprint", contig, start0, end1, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            output.append(
                standard_row(
                    sample=sample,
                    caller="sprint",
                    contig=contig,
                    start0=start0,
                    end1=end1,
                    strand=strand,
                    ref=ref,
                    alt=alt,
                    coverage=coverage,
                    alt_count=alt_count,
                    source=path.name,
                    native_line=line_number,
                )
            )
    return output


def parse_jacusa2(path: Path, sample: str) -> list[dict[str, object]]:
    output = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise StandardizationError(
                    f"line {line_number}: JACUSA2 call-1 output must contain 10 columns"
                )
            contig_raw, start_raw, end_raw, _name, _score, strand, counts_raw, _info, _filter, ref = fields
            try:
                contig = normalize_contig(contig_raw)
                start0 = int(start_raw)
                end1 = int(end_raw)
                counts_list = [int(value) for value in counts_raw.split(",")]
            except (TypeError, ValueError) as exc:
                raise StandardizationError(
                    f"line {line_number}: invalid JACUSA2 coordinate or base counts"
                ) from exc
            if start0 < 0 or end1 != start0 + 1:
                raise StandardizationError(
                    f"line {line_number}: JACUSA2 site must be a single-base BED interval"
                )
            ref = ref.strip().upper()
            if ref not in BASES or len(counts_list) != 4 or any(value < 0 for value in counts_list):
                raise StandardizationError(f"line {line_number}: invalid JACUSA2 base vector")
            counts = dict(zip(BASES, counts_list))
            coverage = sum(counts_list)
            if coverage < 1:
                raise StandardizationError(f"line {line_number}: zero JACUSA2 coverage")
            for alt in BASES:
                if alt == ref or counts[alt] == 0:
                    continue
                key = (sample, "jacusa2", contig, start0, end1, ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                output.append(
                    standard_row(
                        sample=sample,
                        caller="jacusa2",
                        contig=contig,
                        start0=start0,
                        end1=end1,
                        strand=strand.strip() or ".",
                        ref=ref,
                        alt=alt,
                        coverage=coverage,
                        alt_count=counts[alt],
                        source=path.name,
                        native_line=line_number,
                    )
                )
    return output


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", required=True, choices=("reditools", "sprint", "jacusa2"))
    parser.add_argument("--sample", required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> None:
    if args.caller == "reditools":
        rows = parse_reditools(args.input, args.sample)
    elif args.caller == "sprint":
        rows = parse_sprint(args.input, args.sample)
    elif args.caller == "jacusa2":
        rows = parse_jacusa2(args.input, args.sample)
    else:
        raise StandardizationError(f"caller parser not implemented yet: {args.caller}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except (OSError, StandardizationError) as exc:
        print(f"standardization failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

