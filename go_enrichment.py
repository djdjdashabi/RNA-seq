#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path

try:
    from scipy.stats import hypergeom
except ImportError:  # Tiny fixture fallback; production environment pins SciPy.
    hypergeom = None


LINK_REQUIRED = {"sample_id", "caller", "gene_id", "gene_name", "a2i_candidate", "strand_consistent"}
GENE_INFO_REQUIRED = {"#tax_id", "GeneID", "Symbol", "dbXrefs"}
GENE2GO_REQUIRED = {"#tax_id", "GeneID", "GO_ID", "Qualifier", "GO_term", "Category"}


class GoError(ValueError):
    pass


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_dicts(path: Path, required: set[str], label: str):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise GoError(f"{label} missing columns: {', '.join(missing)}")
        yield from reader


def read_gene_info(path: Path):
    ensembl_ids = defaultdict(set)
    symbols = defaultdict(set)
    for row in read_dicts(path, GENE_INFO_REQUIRED, "gene_info"):
        if row["#tax_id"] != "9606":
            continue
        entrez = row["GeneID"]
        symbols[row["Symbol"].upper()].add(entrez)
        for dbxref in row["dbXrefs"].split("|"):
            if dbxref.startswith("Ensembl:"):
                identifier = dbxref.split(":", 1)[1].split(".", 1)[0]
                ensembl_ids[identifier].add(entrez)
    ensembl = {
        identifier: next(iter(ids))
        for identifier, ids in ensembl_ids.items()
        if len(ids) == 1
    }
    ambiguous_ensembl = {
        identifier: ids
        for identifier, ids in ensembl_ids.items()
        if len(ids) > 1
    }
    unique_symbols = {symbol: next(iter(ids)) for symbol, ids in symbols.items() if len(ids) == 1}
    return ensembl, unique_symbols, ambiguous_ensembl


def read_gene2go(path: Path):
    terms = defaultdict(set)
    metadata = {}
    background = set()
    for row in read_dicts(path, GENE2GO_REQUIRED, "gene2go"):
        if row["#tax_id"] != "9606" or "NOT" in row["Qualifier"].split("|"):
            continue
        entrez, go_id = row["GeneID"], row["GO_ID"]
        if not entrez or not go_id:
            raise GoError("gene2go contains an empty GeneID or GO_ID")
        background.add(entrez)
        terms[go_id].add(entrez)
        current = (row["GO_term"], row["Category"])
        if go_id in metadata and metadata[go_id] != current:
            raise GoError(f"gene2go has inconsistent metadata for {go_id}")
        metadata[go_id] = current
    if not background:
        raise GoError("gene2go has no eligible human annotations")
    return terms, metadata, background


def read_links(path: Path, strand_mode: str):
    groups = defaultdict(dict)
    for row in read_dicts(path, LINK_REQUIRED, "site-gene links"):
        if row["a2i_candidate"].lower() != "true":
            continue
        if strand_mode == "consistent" and row["strand_consistent"].lower() != "true":
            continue
        key = (row["sample_id"], row["caller"])
        groups[key][(row["gene_id"].split(".", 1)[0], row["gene_name"])] = None
    return groups


def hypergeom_sf(x: int, population: int, successes: int, draws: int) -> float:
    if hypergeom is not None:
        return float(hypergeom.sf(x - 1, population, successes, draws))
    denominator = math.comb(population, draws)
    upper = min(successes, draws)
    return sum(
        math.comb(successes, hits) * math.comb(population - successes, draws - hits)
        for hits in range(x, upper + 1)
        if 0 <= draws - hits <= population - successes
    ) / denominator


def bh_adjust(pvalues: list[float]) -> list[float]:
    count = len(pvalues)
    adjusted = [1.0] * count
    running = 1.0
    ordered = sorted(enumerate(pvalues), key=lambda item: item[1])
    for rank_index in range(count - 1, -1, -1):
        original_index, pvalue = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, pvalue * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def analyse(
    groups,
    ensembl,
    symbols,
    ambiguous_ensembl,
    terms,
    metadata,
    background,
    strand_mode,
    fdr_threshold,
):
    mapping_rows = []
    enrichment_rows = []
    for (sample, caller), identifiers in sorted(groups.items()):
        mapped = set()
        for gene_id, gene_name in sorted(identifiers):
            if gene_id in ensembl:
                entrez, method = ensembl[gene_id], "ensembl_dbxref"
            elif gene_id in ambiguous_ensembl:
                symbol_match = symbols.get(gene_name.upper())
                if symbol_match in ambiguous_ensembl[gene_id]:
                    entrez, method = symbol_match, "unique_symbol_after_ambiguous_ensembl"
                else:
                    entrez, method = "", "ambiguous_ensembl"
            elif gene_name.upper() in symbols:
                entrez, method = symbols[gene_name.upper()], "unique_symbol"
            else:
                entrez, method = "", "unmapped"
            if entrez:
                mapped.add(entrez)
            mapping_rows.append(
                {
                    "sample_id": sample, "caller": caller, "strand_mode": strand_mode,
                    "gene_id": gene_id, "gene_name": gene_name, "entrez_id": entrez,
                    "mapping_method": method, "metric": "", "value": "",
                }
            )
        selected = mapped & background
        metrics = {
            "input_gene_count": len(identifiers),
            "mapped_gene_count": len(mapped),
            "annotated_gene_count": len(selected),
            "mapping_rate": len(mapped) / len(identifiers) if identifiers else 0.0,
            "annotated_mapping_rate": len(selected) / len(identifiers) if identifiers else 0.0,
        }
        for metric, value in metrics.items():
            mapping_rows.append(
                {
                    "sample_id": sample, "caller": caller, "strand_mode": strand_mode,
                    "gene_id": "", "gene_name": "", "entrez_id": "", "mapping_method": "",
                    "metric": metric,
                    "value": str(value) if isinstance(value, int) else format(value, ".12g"),
                }
            )
        tested = []
        if selected:
            for go_id, term_genes in sorted(terms.items()):
                overlap = selected & term_genes
                pvalue = hypergeom_sf(len(overlap), len(background), len(term_genes), len(selected))
                tested.append((go_id, term_genes, overlap, pvalue))
        adjusted = bh_adjust([item[3] for item in tested])
        for (go_id, term_genes, overlap, pvalue), fdr in zip(tested, adjusted):
            term, category = metadata[go_id]
            enrichment_rows.append(
                {
                    "sample_id": sample, "caller": caller, "strand_mode": strand_mode,
                    "go_id": go_id, "go_term": term, "category": category,
                    "background_size": str(len(background)), "term_size": str(len(term_genes)),
                    "selected_size": str(len(selected)), "overlap_size": str(len(overlap)),
                    "overlap_entrez_ids": ",".join(sorted(overlap, key=int)),
                    "p_value": format(pvalue, ".12g"), "fdr_bh": format(fdr, ".12g"),
                    "significant": "true" if fdr <= fdr_threshold else "false",
                }
            )
    return mapping_rows, enrichment_rows


def write_tsv(path: Path, fields: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--gene-info", type=Path, required=True)
    parser.add_argument("--gene2go", type=Path, required=True)
    parser.add_argument("--enrichment-output", type=Path, required=True)
    parser.add_argument("--mapping-output", type=Path, required=True)
    parser.add_argument("--strand-mode", choices=("all", "consistent"), default="all")
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not 0 < args.fdr <= 1:
        print("GO enrichment failed: --fdr must be in (0, 1]", file=sys.stderr)
        return 2
    try:
        groups = read_links(args.links, args.strand_mode)
        ensembl, symbols, ambiguous_ensembl = read_gene_info(args.gene_info)
        terms, metadata, background = read_gene2go(args.gene2go)
        mapping, enrichment = analyse(
            groups, ensembl, symbols, ambiguous_ensembl, terms, metadata,
            background, args.strand_mode, args.fdr
        )
        write_tsv(
            args.mapping_output,
            ["sample_id", "caller", "strand_mode", "gene_id", "gene_name", "entrez_id", "mapping_method", "metric", "value"],
            mapping,
        )
        write_tsv(
            args.enrichment_output,
            ["sample_id", "caller", "strand_mode", "go_id", "go_term", "category", "background_size", "term_size", "selected_size", "overlap_size", "overlap_entrez_ids", "p_value", "fdr_bh", "significant"],
            enrichment,
        )
    except (OSError, GoError, ValueError) as exc:
        print(f"GO enrichment failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
