#!/usr/bin/env python3
"""Visualize S-LDSC enrichment results for the sldsc_gsea panel."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULT_DIR = Path(
    os.environ.get("TRACECB_ENRICHMENT_DIR", REPO_ROOT / "results/enrichment")
)
METADATA_PATH = REPO_ROOT / "src/figures/metadata.json"

INCREMENTAL_ORDER = ("Original", "traceC increment", "traceCB increment")
ANNOTATION_SHORT_LABELS = {
    "Original": "Original",
    "traceC increment": "traceC\ninc.",
    "traceCB increment": "traceCB\ninc.",
}
CELLTYPE_LABELS = {
    "Monocytes": "Mono",
    "CD4+T_cells": "CD4+ T",
    "CD8+T_cells": "CD8+ T",
    "B_cells": "B",
    "NK_cells": "NK",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    return parser.parse_args()


def load_visual_metadata() -> dict[str, object]:
    with METADATA_PATH.open() as handle:
        return json.load(handle)


def significance_from_p(value: float) -> str:
    if pd.isna(value):
        return ""
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return ""


def significance_from_z(value: float) -> str:
    if pd.isna(value):
        return ""
    abs_value = abs(float(value))
    if abs_value >= 3.29:
        return "***"
    if abs_value >= 2.58:
        return "**"
    if abs_value >= 1.96:
        return "*"
    return ""


def value_annotations(values: pd.DataFrame, stars: pd.DataFrame | None = None) -> pd.DataFrame:
    annotations = values.copy().astype(object)
    for row_index, row in enumerate(values.index):
        for col_index, column in enumerate(values.columns):
            value = values.iat[row_index, col_index]
            if pd.isna(value):
                annotations.iat[row_index, col_index] = ""
                continue
            suffix = "" if stars is None else stars.iat[row_index, col_index]
            annotations.iat[row_index, col_index] = f"{value:.2f}{suffix}"
    return annotations


def save_heatmap(
    values: pd.DataFrame,
    annotations: pd.DataFrame,
    out_path: Path,
    *,
    cmap: str,
    cbar_label: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    sns.set_theme(style="white", font="DejaVu Sans")
    width = max(13.5, 0.54 * values.shape[1] + 5.0)
    height = max(6.4, 0.45 * values.shape[0] + 2.5)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        values,
        annot=annotations,
        fmt="",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.45,
        linecolor="white",
        cbar_kws={"label": cbar_label, "shrink": 0.70},
        annot_kws={"fontsize": 7.5},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0, labelsize=8)
    ax.tick_params(axis="y", labelrotation=0, labelsize=9)
    fig.subplots_adjust(left=0.22, bottom=0.16, right=0.96, top=0.86)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def study_label_from_metadata(qtdid: str, metadata: dict[str, object]) -> str:
    qtdids = list(metadata["QTDids"])
    index = qtdids.index(qtdid)
    celltype = list(metadata["Celltypes"])[index]
    name = str(list(metadata["Names"])[index]).replace("_2014", "").replace("_2017", "").replace("_2021", "")
    return f"{CELLTYPE_LABELS.get(celltype, celltype)} | {name}"


def ordered_incremental_matrix(
    frame: pd.DataFrame,
    value_col: str,
    metadata: dict[str, object],
    include_average: bool = True,
) -> pd.DataFrame:
    frame = frame[frame["TraitLabel"] != "Drug allergy"].copy()
    traits = (
        frame[["TraitOrder", "TraitLabel", "TraitGroup"]]
        .drop_duplicates()
        .sort_values("TraitOrder")
        .reset_index(drop=True)
    )
    study_ids = [study for study in metadata["QTDids"] if study in set(frame["Study"])]
    study_labels = {study: study_label_from_metadata(study, metadata) for study in study_ids}
    columns: list[tuple[str, str]] = []
    for _, trait in traits.iterrows():
        for annotation in INCREMENTAL_ORDER:
            columns.append((trait["TraitLabel"], annotation))

    matrix = pd.DataFrame(
        index=[study_labels[study] for study in study_ids],
        columns=pd.MultiIndex.from_tuples(columns),
    )
    for _, row in frame.iterrows():
        key = (row["TraitLabel"], row["AnnotationLabel"])
        label = study_labels.get(row["Study"])
        if label in matrix.index and key in matrix.columns:
            matrix.loc[label, key] = row[value_col]
    matrix = matrix.astype(float)
    if include_average:
        matrix.loc["Average"] = matrix.mean(axis=0, skipna=True)
    matrix.columns = pd.MultiIndex.from_tuples(matrix.columns)
    return matrix


def flatten_annotation_columns(matrix: pd.DataFrame) -> pd.DataFrame:
    out = matrix.copy()
    out.columns = [ANNOTATION_SHORT_LABELS[annotation] for _, annotation in out.columns]
    return out


def add_top_trait_labels(ax: plt.Axes, matrix: pd.DataFrame) -> None:
    traits = [trait for trait, _ in matrix.columns]
    seen_traits = list(dict.fromkeys(traits))
    top = ax.secondary_xaxis("top")
    centers = []
    labels = []
    for trait in seen_traits:
        indices = [idx for idx, value in enumerate(traits) if value == trait]
        centers.append((min(indices) + max(indices) + 1) / 2)
        labels.append(trait)
    top.set_xticks(centers)
    top.set_xticklabels(labels, fontsize=10)
    top.tick_params(length=0, pad=8)
    top.spines["top"].set_visible(False)
    for boundary in range(3, len(traits), 3):
        ax.axvline(boundary, color="white", linewidth=1.8)


def plot_incremental_enrichment(frame: pd.DataFrame, out_dir: Path) -> None:
    metadata = load_visual_metadata()
    incremental = frame[frame["Model"] == "incremental"].copy()
    enrichment = ordered_incremental_matrix(incremental, "Enrichment", metadata)
    p_values = ordered_incremental_matrix(incremental, "Enrichment_p", metadata, include_average=False)
    p_values.loc["Average"] = np.nan
    enrichment_stars = p_values.reindex_like(enrichment).map(significance_from_p)

    plot_values = flatten_annotation_columns(enrichment)
    plot_annotations = value_annotations(plot_values, flatten_annotation_columns(enrichment_stars))

    sns.set_theme(style="white", font="DejaVu Sans")
    width = max(13.5, 0.54 * plot_values.shape[1] + 5.0)
    height = max(6.4, 0.45 * plot_values.shape[0] + 2.5)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        plot_values,
        annot=plot_annotations,
        fmt="",
        cmap="viridis",
        vmin=0.9,
        vmax=1.3,
        linewidths=0.45,
        linecolor="white",
        cbar_kws={"label": "Enrichment", "shrink": 0.70, "extend": "both"},
        annot_kws={"fontsize": 7.5},
        ax=ax,
    )
    add_top_trait_labels(ax, enrichment)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0, labelsize=7.5)
    ax.tick_params(axis="y", labelrotation=0, labelsize=9)
    fig.subplots_adjust(left=0.22, bottom=0.14, right=0.96, top=0.84)
    fig.savefig(out_dir / "incremental_enrichment_heatmap.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "incremental_enrichment_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(out_dir: Path) -> None:
    body = """# S-LDSC GSEA Visualization

Files in this directory summarize eGene-interval SNP S-LDSC enrichment screens.
The eGene sets and cis intervals are EAS/BBJ-derived. Their binary annotations
and custom LD scores were regenerated on 1000 Genomes Phase 3 EUR SNPs. Each
custom annotation was run separately as `EUR baseline-LD v2.2 + one custom
annotation`, not as a three-column custom joint model. EUR HapMap3 non-MHC
weights and EUR allele frequencies were used.

- `incremental_enrichment_heatmap.*`: study-specific annotations from
  original eGenes, traceC-increment eGenes, and traceCB-increment eGenes. SNPs
  are annotated if they overlap the selected eGene cis intervals. Drug allergy
  is omitted from this preview figure. The bottom row is the study-wise average
  enrichment for each trait/annotation column. The color scale is fixed at
  0.9--1.3 for cross-panel comparison; the colorbar extensions mark clipped
  values, while cell labels retain the estimates.
"""
    (out_dir / "README.md").write_text(body)


def main() -> None:
    args = parse_args()
    summary_dir = args.result_dir / "results" / "summary"
    out_dir = args.result_dir / "visualization"
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(summary_dir / "master_sldsc_gsea_results.csv")
    plot_incremental_enrichment(frame, out_dir)
    for stale in (
        "incremental_tau_z_heatmap.pdf",
        "incremental_tau_z_heatmap.png",
        "overall_enrichment_95CI.pdf",
        "overall_enrichment_95CI.png",
    ):
        path = out_dir / stale
        if path.exists():
            path.unlink()
    write_readme(out_dir)
    print(f"[done] wrote visualizations under {out_dir}", flush=True)


if __name__ == "__main__":
    main()
