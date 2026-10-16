#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path


BIN_SIZE = 100_000
REGION_ORDER = {
    "CDS": 0,
    "five_prime_UTR": 1,
    "three_prime_UTR": 2,
    "noncoding_exon": 3,
    "intron": 4,
    "intergenic": 5,
}
REQUIRED = {
    "sample_id",
    "caller",
    "contig",
    "start0",
    "end1",
    "ref",
    "alt",
    "a2i_candidate",
}


class AnnotationError(ValueError):
    pass


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def normalize_contig(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return "MT" if value.upper() == "M" else value


def parse_attributes(value: str, line_number: int) -> dict[str, str]:
    attributes = {}
    for part in value.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        match = re.fullmatch(r'(\S+)\s+"([^"]*)"', part)
        if not match:
            raise AnnotationError(f"GTF line {line_number}: malformed attribute {part!r}")
        attributes[match.group(1)] = match.group(2)
    return attributes


def parse_gtf(path: Path):
    genes: dict[str, dict[str, object]] = {}
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise AnnotationError(f"GTF line {line_number}: expected 9 columns, found {len(fields)}")
            contig, _, feature, start, end, _, strand, _, raw_attributes = fields
            try:
                start0, end1 = int(start) - 1, int(end)
            except ValueError as exc:
                raise AnnotationError(f"GTF line {line_number}: invalid coordinates") from exc
            if start0 < 0 or end1 <= start0 or strand not in {"+", "-", "."}:
                raise AnnotationError(f"GTF line {line_number}: invalid interval or strand")
            attributes = parse_attributes(raw_attributes, line_number)
            gene_id = attributes.get("gene_id", "")
            if not gene_id:
                raise AnnotationError(f"GTF line {line_number}: missing gene_id")
            gene = genes.setdefault(
                gene_id,
                {
                    "gene_id": gene_id,
                    "gene_name": attributes.get("gene_name", ""),
                    "gene_biotype": attributes.get("gene_biotype", attributes.get("gene_type", "")),
                    "contig": normalize_contig(contig),
                    "strand": strand,
                    "bounds": None,
                    "features": defaultdict(list),
                },
            )
            if gene["contig"] != normalize_contig(contig) or gene["strand"] != strand:
                raise AnnotationError(f"GTF line {line_number}: inconsistent gene coordinates for {gene_id}")
            if attributes.get("gene_name") and not gene["gene_name"]:
                gene["gene_name"] = attributes["gene_name"]
            if feature == "gene":
                gene["bounds"] = (start0, end1)
            elif feature in {"CDS", "five_prime_utr", "three_prime_utr", "exon"}:
                gene["features"][feature].append((start0, end1))

    bins = defaultdict(set)
    for gene_id, gene in genes.items():
        bounds = gene["bounds"]
        if bounds is None:
            all_intervals = [interval for intervals in gene["features"].values() for interval in intervals]
            if not all_intervals:
                raise AnnotationError(f"GTF gene {gene_id}: no gene bounds or features")
            bounds = (min(item[0] for item in all_intervals), max(item[1] for item in all_intervals))
            gene["bounds"] = bounds
        for bin_number in range(bounds[0] // BIN_SIZE, (bounds[1] - 1) // BIN_SIZE + 1):
            bins[(gene["contig"], bin_number)].add(gene_id)
    return genes, bins


def overlaps(interval: tuple[int, int], start0: int, end1: int) -> bool:
    return interval[0] < end1 and start0 < interval[1]


def gene_region(gene: dict[str, object], start0: int, end1: int) -> str:
    features = gene["features"]
    for feature, label in (
        ("CDS", "CDS"),
        ("five_prime_utr", "five_prime_UTR"),
        ("three_prime_utr", "three_prime_UTR"),
        ("exon", "noncoding_exon"),
    ):
        if any(overlaps(interval, start0, end1) for interval in features.get(feature, [])):
            return label
    return "intron"


def strand_consistency(ref: str, alt: str, gene_strand: str) -> str:
    if (ref, alt) == ("A", "G"):
        return "true" if gene_strand == "+" else "false"
    if (ref, alt) == ("T", "C"):
        return "true" if gene_strand == "-" else "false"
    return "false"


def read_candidates(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        missing = sorted(REQUIRED - set(columns))
        if missing:
            raise AnnotationError(f"{path}: canonical input missing columns: {', '.join(missing)}")
        rows = list(reader)
    return columns, rows


def annotate(columns, rows, genes, bins):
    annotated = []
    links = []
    for input_line, row in enumerate(rows, 2):
        try:
            start0, end1 = int(row["start0"]), int(row["end1"])
        except ValueError as exc:
            raise AnnotationError(f"candidate line {input_line}: invalid coordinates") from exc
        if start0 < 0 or end1 <= start0:
            raise AnnotationError(f"candidate line {input_line}: invalid interval")
        contig = normalize_contig(row["contig"])
        candidates = []
        for gene_id in bins.get((contig, start0 // BIN_SIZE), set()):
            gene = genes[gene_id]
            if overlaps(gene["bounds"], start0, end1):
                candidates.append(gene)
        candidates.sort(key=lambda gene: gene["gene_id"])
        regions = []
        consistencies = []
        for gene in candidates:
            region = gene_region(gene, start0, end1)
            consistent = strand_consistency(row["ref"], row["alt"], gene["strand"])
            regions.append(region)
            consistencies.append(consistent)
            links.append(
                {
                    "sample_id": row["sample_id"],
                    "caller": row["caller"],
                    "contig": contig,
                    "start0": str(start0),
                    "end1": str(end1),
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "a2i_candidate": row["a2i_candidate"],
                    "gene_id": gene["gene_id"],
                    "gene_name": gene["gene_name"],
                    "gene_biotype": gene["gene_biotype"],
                    "gene_strand": gene["strand"],
                    "region": region,
                    "strand_consistent": consistent,
                }
            )
        new_row = dict(row)
        new_row["contig"] = contig
        new_row["gene_ids"] = ",".join(gene["gene_id"] for gene in candidates)
        new_row["gene_names"] = ",".join(gene["gene_name"] for gene in candidates)
        new_row["region"] = min(regions, key=REGION_ORDER.get) if regions else "intergenic"
        new_row["strand_consistent"] = (
            "true" if "true" in consistencies else "false" if consistencies else ""
        )
        annotated.append(new_row)
    return annotated, links


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--sites-output", type=Path, required=True)
    parser.add_argument("--links-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        genes, bins = parse_gtf(args.gtf)
        columns, rows = read_candidates(args.input)
        annotated, links = annotate(columns, rows, genes, bins)
        write_tsv(args.sites_output, columns + ["gene_ids", "gene_names", "region"], annotated)
        write_tsv(
            args.links_output,
            [
                "sample_id", "caller", "contig", "start0", "end1", "ref", "alt",
                "a2i_candidate", "gene_id", "gene_name", "gene_biotype", "gene_strand",
                "region", "strand_consistent",
            ],
            links,
        )
    except (OSError, AnnotationError) as exc:
        print(f"annotation failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

