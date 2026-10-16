#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mpl_colors
from matplotlib import ticker


REQUIRED = {"sample_id", "caller", "scope", "metric", "group", "value", "significant", "source"}
CALLERS = ("reditools", "sprint", "jacusa2")
CALLER_LABELS = {"reditools": "REDItools", "sprint": "SPRINT", "jacusa2": "JACUSA2"}
CALLER_COLORS = {"reditools": "#4477AA", "sprint": "#CC6677", "jacusa2": "#228833"}
SAMPLE_STYLES = ("-", "--", ":", "-.")
SAMPLE_MARKERS = ("o", "s", "^", "D")

SUBSTITUTIONS = ("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
SUBSTITUTION_COLORS = {
    "A>C": "#A6CEE3", "A>G": "#1F78B4", "A>T": "#B2DF8A",
    "C>A": "#66B885", "C>G": "#F4B6B2", "C>T": "#D95F5F",
    "G>A": "#F2C879", "G>C": "#E69F45", "G>T": "#CAB2D6",
    "T>A": "#8C6BB1", "T>C": "#17A9B8", "T>G": "#B8B8B8",
}
REGIONS = ("CDS", "five_prime_UTR", "three_prime_UTR", "noncoding_exon", "intron", "intergenic")
REGION_LABELS = {
    "CDS": "CDS", "five_prime_UTR": "5′UTR", "three_prime_UTR": "3′UTR",
    "noncoding_exon": "Non-coding exon", "intron": "Intron", "intergenic": "Intergenic",
}
REGION_COLORS = {
    "CDS": "#4C78A8", "five_prime_UTR": "#72B7B2", "three_prime_UTR": "#54A24B",
    "noncoding_exon": "#B8DE9C", "intron": "#F2CF5B", "intergenic": "#D9D9D9",
}
CATEGORY_ORDER = ("Process", "Function", "Component")
CATEGORY_LABELS = {"Process": "BP", "Function": "MF", "Component": "CC"}


class PlotError(ValueError):
    pass


def configure_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 6.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "grid.color": "#E6E6E6",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_summary(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED - set(reader.fieldnames or []))
        if missing:
            raise PlotError("canonical summary missing columns: " + ", ".join(missing))
        rows = list(reader)
    for line_number, row in enumerate(rows, 2):
        try:
            float(row["value"])
        except ValueError as exc:
            raise PlotError(f"canonical summary line {line_number}: nonnumeric value") from exc
    return rows


def short_sample(sample: str) -> str:
    return sample.replace("SRR370435", "") if sample.startswith("SRR370435") else sample


def ordered_samples(rows):
    return sorted({row["sample_id"] for row in rows})


def panel_label(axis, label: str):
    axis.text(-0.12, 1.06, label, transform=axis.transAxes, fontsize=11, fontweight="bold", va="top")


def soften_axis(axis, *, xgrid=False, ygrid=True):
    axis.set_axisbelow(True)
    axis.grid(axis="x" if xgrid and not ygrid else "y" if ygrid and not xgrid else "both")


def save_figure(fig, output_dir: Path, number: int):
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"figure{number}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"figure{number}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure_method_features(rows, output_dir):
    distributions = [row for row in rows if row["metric"] == "editing_level_log10_fraction" and row["scope"] == "a2i"]
    spectra = [row for row in rows if row["metric"] == "substitution_fraction" and row["scope"] == "all"]
    regions = [row for row in rows if row["metric"] == "region_fraction" and row["scope"] == "a2i"]
    alu = [row for row in rows if row["metric"] == "alu_fraction" and row["scope"] == "a2i"]
    if not distributions or not spectra or not regions or not alu:
        raise PlotError("figure 1 requires editing distributions, substitution spectrum, region fractions, and Alu fractions")

    samples = ordered_samples(distributions)
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    for caller in CALLERS:
        for sample_index, sample in enumerate(samples):
            selected = sorted(
                (row for row in distributions if row["caller"] == caller and row["sample_id"] == sample),
                key=lambda row: float(row["group"]),
            )
            if not selected:
                continue
            ax_a.plot(
                [float(row["group"]) for row in selected],
                [float(row["value"]) for row in selected],
                color=CALLER_COLORS[caller],
                linestyle=SAMPLE_STYLES[sample_index % len(SAMPLE_STYLES)],
                linewidth=1.7,
                label=f"{CALLER_LABELS[caller]} · {short_sample(sample)}",
            )
    ax_a.set_xlabel("log$_{10}$(editing level)")
    ax_a.set_ylabel("Fraction of potential A-to-I alleles")
    ax_a.set_title("Editing-level distribution", loc="left")
    ax_a.set_xlim(-5, 0)
    ax_a.legend(frameon=False, ncol=2, handlelength=2.6)
    soften_axis(ax_a)
    panel_label(ax_a, "A")

    bar_groups = [(sample, caller) for sample in samples for caller in CALLERS]
    x = np.arange(len(bar_groups))
    bottom = np.zeros(len(bar_groups))
    for substitution in SUBSTITUTIONS:
        values = np.array(
            [
                next(
                    (
                        float(row["value"])
                        for row in spectra
                        if row["sample_id"] == sample and row["caller"] == caller and row["group"] == substitution
                    ),
                    0.0,
                )
                for sample, caller in bar_groups
            ]
        )
        ax_b.bar(x, values, bottom=bottom, width=0.76, color=SUBSTITUTION_COLORS[substitution], label=substitution, linewidth=0)
        bottom += values
    ax_b.set_xticks(x, [f"{short_sample(sample)}\n{CALLER_LABELS[caller]}" for sample, caller in bar_groups], rotation=35, ha="right")
    ax_b.set_ylim(0, 1)
    ax_b.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    ax_b.set_ylabel("Allele fraction")
    ax_b.set_title("Substitution spectrum", loc="left")
    ax_b.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.28), columnspacing=0.8, handlelength=1.2)
    soften_axis(ax_b)
    panel_label(ax_b, "B")

    bottom = np.zeros(len(bar_groups))
    for region in REGIONS:
        values = np.array(
            [
                next(
                    (
                        float(row["value"])
                        for row in regions
                        if row["sample_id"] == sample and row["caller"] == caller and row["group"] == region
                    ),
                    0.0,
                )
                for sample, caller in bar_groups
            ]
        )
        ax_c.bar(x, values, bottom=bottom, width=0.76, color=REGION_COLORS[region], label=REGION_LABELS[region], linewidth=0)
        bottom += values
    ax_c.set_xticks(x, [f"{short_sample(sample)}\n{CALLER_LABELS[caller]}" for sample, caller in bar_groups], rotation=35, ha="right")
    ax_c.set_ylim(0, 1)
    ax_c.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    ax_c.set_ylabel("Potential A-to-I site fraction")
    ax_c.set_title("Genomic-region composition", loc="left")
    ax_c.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.28), columnspacing=0.8, handlelength=1.2)
    soften_axis(ax_c)
    panel_label(ax_c, "C")

    max_alu = max(float(row["value"]) for row in alu)
    alu_label_offset = max(max_alu * 0.025, 0.002)
    for caller_index, caller in enumerate(CALLERS):
        caller_values = []
        for sample_index, sample in enumerate(samples):
            value = next(
                float(row["value"])
                for row in alu
                if row["sample_id"] == sample and row["caller"] == caller
            )
            caller_values.append(value)
            offset = (sample_index - (len(samples) - 1) / 2) * 0.14
            ax_d.scatter(
                caller_index + offset,
                value,
                s=42,
                marker=SAMPLE_MARKERS[sample_index % len(SAMPLE_MARKERS)],
                color=CALLER_COLORS[caller],
                edgecolor="white",
                linewidth=0.8,
                zorder=3,
                label=short_sample(sample) if caller_index == 0 else None,
            )
            ax_d.text(caller_index + offset, value + alu_label_offset, f"{value:.1%}", ha="center", va="bottom", fontsize=6)
        if len(caller_values) > 1:
            offsets = [(index - (len(samples) - 1) / 2) * 0.14 for index in range(len(samples))]
            ax_d.plot(np.array(offsets) + caller_index, caller_values, color="#A0A0A0", linewidth=0.8, zorder=1)
        ax_d.scatter(caller_index, float(np.mean(caller_values)), s=58, marker="D", facecolor="white", edgecolor=CALLER_COLORS[caller], linewidth=1.2, zorder=4)
    ax_d.set_xticks(range(len(CALLERS)), [CALLER_LABELS[caller] for caller in CALLERS])
    ax_d.set_ylim(0, max(max_alu * 1.18, 0.05))
    ax_d.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    ax_d.set_ylabel("Potential A-to-I sites in Alu")
    ax_d.set_title("Alu overlap by biological replicate", loc="left")
    ax_d.legend(title="Replicate", frameon=False, loc="upper left")
    soften_axis(ax_d)
    panel_label(ax_d, "D")

    save_figure(fig, output_dir, 1)
    return distributions + spectra + regions + alu


def figure_upset(rows, output_dir):
    intersections = [row for row in rows if row["metric"] == "upset_site_count" and row["scope"] == "a2i"]
    set_sizes = [row for row in rows if row["metric"] == "site_count" and row["scope"] == "a2i"]
    if not intersections or not set_sizes:
        raise PlotError("figure 2 requires a2i upset_site_count and site_count rows")
    samples = ordered_samples(intersections)
    combinations = (
        frozenset(("reditools",)), frozenset(("sprint",)), frozenset(("jacusa2",)),
        frozenset(("reditools", "sprint")), frozenset(("reditools", "jacusa2")),
        frozenset(("sprint", "jacusa2")), frozenset(CALLERS),
    )
    fig = plt.figure(figsize=(9.0, 3.7 * len(samples)), constrained_layout=True)
    outer = fig.add_gridspec(len(samples), 2, width_ratios=(1.15, 3.6), hspace=0.24, wspace=0.12)

    for sample_index, sample in enumerate(samples):
        size_axis = fig.add_subplot(outer[sample_index, 0])
        right = outer[sample_index, 1].subgridspec(2, 1, height_ratios=(3.1, 1.15), hspace=0.04)
        bar_axis = fig.add_subplot(right[0])
        matrix_axis = fig.add_subplot(right[1], sharex=bar_axis)

        sizes = [
            next(float(row["value"]) for row in set_sizes if row["sample_id"] == sample and row["caller"] == caller)
            for caller in CALLERS
        ]
        y = np.arange(len(CALLERS))
        size_axis.barh(y, sizes, color=[CALLER_COLORS[caller] for caller in CALLERS], height=0.62)
        size_axis.set_yticks(y, [CALLER_LABELS[caller] for caller in CALLERS])
        size_axis.invert_xaxis()
        size_axis.invert_yaxis()
        if all(value > 0 for value in sizes):
            size_axis.set_xscale("log")
            size_axis.set_xlabel("Set size (log scale)")
        else:
            size_axis.set_xscale("symlog", linthresh=1)
            size_axis.set_xlabel("Set size (symlog scale)")
        size_axis.set_title("Caller sets", loc="left")
        size_axis.grid(axis="x")
        largest_size = max(sizes, default=0)
        for index, value in enumerate(sizes):
            if value == 0:
                size_axis.text(
                    0.98,
                    index,
                    "0",
                    transform=size_axis.get_yaxis_transform(),
                    va="center",
                    ha="right",
                    color="#2A2A2A",
                    fontsize=6.5,
                    fontweight="bold",
                )
                continue
            is_small = value < largest_size * 0.10
            size_axis.text(
                value * (1.18 if is_small else 0.88),
                index,
                f"{int(value):,}",
                va="center",
                ha="right" if is_small else "left",
                color="#2A2A2A" if is_small else "white",
                fontsize=6.5,
                fontweight="bold",
            )

        count_by_combo = {
            frozenset(filter(None, row["group"].split("&"))): int(float(row["value"]))
            for row in intersections
            if row["sample_id"] == sample
        }
        counts = [count_by_combo.get(combination, 0) for combination in combinations]
        x = np.arange(len(combinations))
        colors = [
            CALLER_COLORS[next(iter(combination))] if len(combination) == 1 else "#6F5D7B"
            for combination in combinations
        ]
        bar_axis.bar(x, counts, color=colors, width=0.72)
        positive = [value for value in counts if value > 0]
        if positive and max(positive) / min(positive) > 100:
            bar_axis.set_yscale("log")
            bar_axis.set_ylabel("Exclusive intersection\nsize (log scale)")
        else:
            bar_axis.set_ylabel("Exclusive intersection size")
        bar_axis.set_title(f"{sample}: potential A-to-I caller intersections", loc="left", fontweight="bold")
        bar_axis.tick_params(axis="x", labelbottom=False)
        bar_axis.grid(axis="y")
        for xpos, value in zip(x, counts):
            if value > 0:
                bar_axis.annotate(f"{value:,}", (xpos, value), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=6)

        for xpos, combination in zip(x, combinations):
            selected_indices = []
            for caller_index, caller in enumerate(CALLERS):
                active = caller in combination
                matrix_axis.scatter(
                    xpos,
                    caller_index,
                    s=32 if active else 20,
                    color=CALLER_COLORS[caller] if active else "#D9D9D9",
                    edgecolor="white" if active else "none",
                    linewidth=0.5,
                    zorder=3,
                )
                if active:
                    selected_indices.append(caller_index)
            if len(selected_indices) > 1:
                matrix_axis.plot([xpos, xpos], [min(selected_indices), max(selected_indices)], color="#555555", linewidth=1.2, zorder=2)
        matrix_axis.set_yticks(range(len(CALLERS)), [CALLER_LABELS[caller] for caller in CALLERS])
        matrix_axis.set_ylim(len(CALLERS) - 0.5, -0.5)
        matrix_axis.set_xticks(x, [str(index + 1) for index in x])
        matrix_axis.set_xlabel("Exclusive caller combination")
        matrix_axis.grid(False)
        matrix_axis.spines["left"].set_visible(False)
        matrix_axis.spines["bottom"].set_visible(False)
        panel_label(bar_axis, chr(ord("A") + sample_index))

    save_figure(fig, output_dir, 2)
    return intersections + set_sizes


def figure_jaccard(rows, output_dir):
    jaccard_rows = [row for row in rows if row["metric"] == "pairwise_jaccard" and row["scope"] == "a2i"]
    intersection_rows = [row for row in rows if row["metric"] == "pairwise_intersection" and row["scope"] == "a2i"]
    if not jaccard_rows or not intersection_rows:
        raise PlotError("figure 3 requires a2i pairwise_jaccard and pairwise_intersection rows")
    samples = ordered_samples(jaccard_rows)
    fig, axes = plt.subplots(1, len(samples), figsize=(4.2 * len(samples) + 0.7, 4.0), constrained_layout=True, squeeze=False)
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("white")
    image = None
    for sample_index, (axis, sample) in enumerate(zip(axes[0], samples)):
        matrix = np.full((len(CALLERS), len(CALLERS)), np.nan)
        intersections = np.zeros((len(CALLERS), len(CALLERS)), dtype=int)
        for row in jaccard_rows:
            if row["sample_id"] != sample:
                continue
            left, right = row["caller"].split("|", 1)
            i, j = CALLERS.index(left), CALLERS.index(right)
            matrix[max(i, j), min(i, j)] = float(row["value"])
        for row in intersection_rows:
            if row["sample_id"] != sample:
                continue
            left, right = row["caller"].split("|", 1)
            i, j = CALLERS.index(left), CALLERS.index(right)
            intersections[max(i, j), min(i, j)] = int(float(row["value"]))
        image = axis.pcolormesh(
            np.arange(len(CALLERS) + 1) - 0.5,
            np.arange(len(CALLERS) + 1) - 0.5,
            np.ma.masked_invalid(matrix),
            vmin=0,
            vmax=1,
            cmap=cmap,
            shading="flat",
            rasterized=False,
        )
        for index in range(len(CALLERS)):
            axis.add_patch(plt.Rectangle((index - 0.5, index - 0.5), 1, 1, facecolor="#F2F2F2", edgecolor="white"))
            axis.text(index, index, "1.000\nself", ha="center", va="center", color="#666666", fontsize=7)
        for i in range(len(CALLERS)):
            for j in range(i):
                value = matrix[i, j]
                color = "white" if value >= 0.48 else "#222222"
                axis.text(j, i, f"{value:.3f}\nn={intersections[i, j]:,}", ha="center", va="center", color=color, fontsize=7)
        axis.set_xticks(range(len(CALLERS)), [CALLER_LABELS[caller] for caller in CALLERS], rotation=35, ha="right")
        axis.set_yticks(range(len(CALLERS)), [CALLER_LABELS[caller] for caller in CALLERS])
        axis.set_xlim(-0.5, len(CALLERS) - 0.5)
        axis.set_ylim(len(CALLERS) - 0.5, -0.5)
        axis.set_aspect("equal")
        axis.set_title(sample, fontweight="bold")
        axis.spines["top"].set_visible(True)
        axis.spines["right"].set_visible(True)
        panel_label(axis, chr(ord("A") + sample_index))
    fig.colorbar(image, ax=axes.ravel().tolist(), label="Jaccard index", shrink=0.8, pad=0.03)
    fig.suptitle("Concordance of potential A-to-I candidate sets", fontsize=10, fontweight="bold")
    save_figure(fig, output_dir, 3)
    return jaccard_rows + intersection_rows


def parse_go_group(group: str):
    parts = group.split("|", 2)
    if len(parts) != 3:
        raise PlotError(f"invalid GO group encoding: {group!r}")
    return parts[0], parts[1], parts[2]


def figure_go(rows, output_dir):
    all_fdr_rows = [
        row for row in rows
        if row["metric"] == "go_fdr_bh" and row["scope"] == "a2i:all"
    ]
    all_overlap_rows = [
        row for row in rows
        if row["metric"] == "go_overlap_size" and row["scope"] == "a2i:all"
    ]
    for row in all_fdr_rows + all_overlap_rows:
        if row["significant"].lower() not in {"true", "false"}:
            raise PlotError("GO rows contain an invalid significant flag")
    fdr_rows = [row for row in all_fdr_rows if row["significant"].lower() == "true"]
    fig, axis = plt.subplots(figsize=(9.0, 6.8), constrained_layout=True)
    if not fdr_rows:
        axis.axis("off")
        axis.text(0.5, 0.5, "No GO terms significant at BH-FDR 0.05", ha="center", va="center", fontsize=11)
        selected_fdr = []
        selected_overlap = []
    else:
        best_by_term = {}
        for row in fdr_rows:
            go_id, category, term = parse_go_group(row["group"])
            key = (go_id, category, term)
            best_by_term[key] = min(best_by_term.get(key, 1.0), float(row["value"]))
        selected_terms = []
        for category in CATEGORY_ORDER:
            category_terms = sorted(
                (item for item in best_by_term if item[1] == category),
                key=lambda item: (best_by_term[item], item[0]),
            )[:4]
            selected_terms.extend(category_terms)
        if not selected_terms:
            raise PlotError("significant GO rows use unknown ontology categories")

        selected_groups = {"|".join(term) for term in selected_terms}
        selected_fdr = [row for row in all_fdr_rows if row["group"] in selected_groups]
        selected_overlap = [row for row in all_overlap_rows if row["group"] in selected_groups]
        overlap_lookup = {}
        for row in selected_overlap:
            key = (row["sample_id"], row["caller"], row["group"])
            if key in overlap_lookup:
                raise PlotError(f"duplicate GO overlap row: {key}")
            overlap_lookup[key] = (
                float(row["value"]),
                row["significant"].lower(),
            )
        samples = ordered_samples(selected_fdr)
        columns = [(sample, caller) for sample in samples for caller in CALLERS]
        x_lookup = {key: index for index, key in enumerate(columns)}
        y_lookup = {"|".join(term): index for index, term in enumerate(selected_terms)}

        points = []
        fdr_lookup = {}
        for row in selected_fdr:
            key = (row["sample_id"], row["caller"], row["group"])
            if key in fdr_lookup:
                raise PlotError(f"duplicate GO FDR row: {key}")
            fdr_lookup[key] = row["significant"].lower()
        if set(fdr_lookup) != set(overlap_lookup):
            missing = sorted(set(fdr_lookup) - set(overlap_lookup))
            orphaned = sorted(set(overlap_lookup) - set(fdr_lookup))
            raise PlotError(
                f"GO FDR/overlap key mismatch; missing={missing[:3]}, orphaned={orphaned[:3]}"
            )
        for row in selected_fdr:
            key = (row["sample_id"], row["caller"], row["group"])
            overlap, overlap_significant = overlap_lookup[key]
            if overlap_significant != fdr_lookup[key]:
                raise PlotError(f"GO FDR/overlap significance mismatch: {key}")
            points.append(
                (
                    x_lookup[(row["sample_id"], row["caller"])],
                    y_lookup[row["group"]],
                    -math.log10(max(float(row["value"]), np.finfo(float).tiny)),
                    overlap,
                )
            )
        max_overlap = max((point[3] for point in points), default=1.0)
        color_values = [point[2] for point in points]
        norm = mpl_colors.Normalize(vmin=min(color_values), vmax=max(color_values) if max(color_values) > min(color_values) else min(color_values) + 1)
        sizes = [25 + 190 * math.sqrt(point[3] / max_overlap) for point in points]
        scatter = axis.scatter(
            [point[0] for point in points],
            [point[1] for point in points],
            c=color_values,
            s=sizes,
            cmap="cividis",
            norm=norm,
            edgecolor="white",
            linewidth=0.7,
        )
        labels = [
            f"[{CATEGORY_LABELS[category]}] {textwrap.fill(term, 34)}\n{go_id}"
            for go_id, category, term in selected_terms
        ]
        axis.set_xticks(
            range(len(columns)),
            [f"{short_sample(sample)}\n{CALLER_LABELS[caller]}" for sample, caller in columns],
            rotation=30,
            ha="right",
        )
        axis.set_yticks(range(len(selected_terms)), labels)
        axis.set_ylim(len(selected_terms) - 0.5, -0.5)
        axis.set_xlabel("Biological replicate and caller")
        axis.set_title("Exploratory GO enrichment of potential A-to-I candidates", loc="left", fontweight="bold")
        axis.grid(color="#ECECEC", linewidth=0.6)
        colorbar = fig.colorbar(scatter, ax=axis, pad=0.02, shrink=0.78)
        colorbar.set_label("−log$_{10}$(BH-FDR)")
        legend_values = sorted({max(1, int(round(max_overlap * fraction))) for fraction in (0.25, 0.5, 1.0)})
        handles = [
            axis.scatter([], [], s=25 + 190 * math.sqrt(value / max_overlap), color="#777777", edgecolor="white", label=f"{value:,}")
            for value in legend_values
        ]
        axis.legend(handles=handles, title="Overlap genes", frameon=False, loc="upper left", bbox_to_anchor=(1.13, 1.0), labelspacing=1.2)
        panel_label(axis, "A")
    save_figure(fig, output_dir, 4)
    return selected_fdr + selected_overlap


def write_figure_data(path: Path, summary_path: Path, groups):
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    fields = ["figure", "sample_id", "caller", "scope", "metric", "group", "value", "significant", "source", "canonical_sha256"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for number, rows in enumerate(groups, 1):
            seen = set()
            for row in rows:
                key = tuple(row[field] for field in ("sample_id", "caller", "scope", "metric", "group", "value", "significant", "source"))
                if key in seen:
                    continue
                seen.add(key)
                writer.writerow({"figure": str(number), **row, "canonical_sha256": digest})


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    configure_style()
    try:
        rows = read_summary(args.summary)
        groups = [
            figure_method_features(rows, args.output_dir),
            figure_upset(rows, args.output_dir),
            figure_jaccard(rows, args.output_dir),
            figure_go(rows, args.output_dir),
        ]
        write_figure_data(args.data_output, args.summary, groups)
    except (OSError, PlotError, ValueError, StopIteration) as exc:
        print(f"plotting failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
