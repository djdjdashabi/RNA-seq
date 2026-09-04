#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import sys
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_DOWNLOADS = Path.home() / "essay" / "resources" / "downloads"

SCHEMAS = {
    "gene_info": {"#tax_id", "GeneID", "Symbol", "dbXrefs"},
    "gene2go": {"#tax_id", "GeneID", "GO_ID", "Qualifier", "GO_term", "Category"},
}


class BuildError(ValueError):
    """Raised when an NCBI input does not satisfy its documented schema."""


def resolve_default_source(filename: str) -> Path:
    candidates = [PROJECT_ROOT / "annotation" / filename, DEFAULT_DOWNLOADS / filename]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    locations = ", ".join(str(path) for path in candidates)
    raise BuildError(f"cannot find {filename}; checked: {locations}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_source(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def write_human_subset(dataset: str, source: Path, output: Path):
    if not source.is_file():
        raise BuildError(f"{dataset} source does not exist: {source}")
    if source.resolve() == output.resolve():
        raise BuildError(f"{dataset} source and output must differ")
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{output.name}.", dir=output.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    row_count = 0
    try:
        with open_source(source) as input_handle:
            reader = csv.reader(input_handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise BuildError(f"{dataset} is empty") from exc
            if len(header) != len(set(header)):
                raise BuildError(f"{dataset} contains duplicate column names")
            missing = sorted(SCHEMAS[dataset] - set(header))
            if missing:
                raise BuildError(f"{dataset} missing columns: {', '.join(missing)}")
            tax_index = header.index("#tax_id")
            if output.suffix == ".gz":
                output_handle = gzip.open(temporary, "wt", encoding="utf-8", newline="")
            else:
                output_handle = temporary.open("w", encoding="utf-8", newline="")
            with output_handle:
                writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                for line_number, row in enumerate(reader, start=2):
                    if len(row) != len(header):
                        raise BuildError(
                            f"{dataset} line {line_number} has {len(row)} columns; expected {len(header)}"
                        )
                    if row[tax_index] == "9606":
                        writer.writerow(row)
                        row_count += 1
        if row_count == 0:
            raise BuildError(f"{dataset} contains no human tax_id 9606 rows")
        return temporary, row_count, sha256_file(source)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def write_manifest(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["dataset", "source", "source_sha256", "output", "human_rows"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-info-source", type=Path)
    parser.add_argument("--gene2go-source", type=Path)
    parser.add_argument(
        "--gene-info-output",
        type=Path,
        default=PROJECT_ROOT / "annotation" / "gene_info_human.tsv.gz",
    )
    parser.add_argument(
        "--gene2go-output",
        type=Path,
        default=PROJECT_ROOT / "annotation" / "gene2go_human.tsv.gz",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=PROJECT_ROOT / "annotation" / "gene2go_build_manifest.tsv",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    temporary_outputs = []
    try:
        gene_info_source = args.gene_info_source or resolve_default_source("Homo_sapiens.gene_info.gz")
        gene2go_source = args.gene2go_source or resolve_default_source("gene2go.gz")
        info_tmp, info_count, info_hash = write_human_subset(
            "gene_info", gene_info_source, args.gene_info_output
        )
        temporary_outputs.append(info_tmp)
        go_tmp, go_count, go_hash = write_human_subset("gene2go", gene2go_source, args.gene2go_output)
        temporary_outputs.append(go_tmp)
        os.replace(info_tmp, args.gene_info_output)
        temporary_outputs.remove(info_tmp)
        os.replace(go_tmp, args.gene2go_output)
        temporary_outputs.remove(go_tmp)
        write_manifest(
            args.manifest,
            [
                {
                    "dataset": "gene_info",
                    "source": str(gene_info_source.resolve()),
                    "source_sha256": info_hash,
                    "output": str(args.gene_info_output.resolve()),
                    "human_rows": info_count,
                },
                {
                    "dataset": "gene2go",
                    "source": str(gene2go_source.resolve()),
                    "source_sha256": go_hash,
                    "output": str(args.gene2go_output.resolve()),
                    "human_rows": go_count,
                },
            ],
        )
        print(f"gene_info: {info_count} human rows -> {args.gene_info_output}")
        print(f"gene2go: {go_count} human rows -> {args.gene2go_output}")
    except (OSError, BuildError, csv.Error) as exc:
        for temporary in temporary_outputs:
            temporary.unlink(missing_ok=True)
        print(f"gene2go build failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
