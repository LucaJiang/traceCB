"""Shared GSEApy ORA utilities for traceCB enrichment figures."""

from __future__ import annotations

import json
import math
import re
import sys
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import gseapy as gp
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from gseapy.parser import read_gmt
from matplotlib import colors as mcolors
from matplotlib import transforms as mtransforms
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402

SRC_DIR = Path(__file__).resolve().parents[1]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from visual.utils import geneid2name, load_all_summary, meta_data, p2z  # noqa: E402

MIN_P = 1e-300
DEFAULT_RESULT_ROOT = Path(
    "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea"
)
DEFAULT_GMT_DIR = Path("/home/wjiang49/group/wjiang49/data/gsea_gmt")

CELLTYPE_ORDER = tuple(dict.fromkeys(meta_data.get("Celltypes", []))) + ("Other",)
CELLTYPE_RANK = {celltype: i for i, celltype in enumerate(CELLTYPE_ORDER)}
QTD_RANK = {qtdid: i for i, qtdid in enumerate(meta_data.get("QTDids", []))}
CELLTYPE_LABELS = {
    "Monocytes": "Mono",
    "CD4+T_cells": "CD4+ T",
    "CD8+T_cells": "CD8+ T",
    "B_cells": "B",
    "NK_cells": "NK",
    "Other": "Other",
}
CELLTYPE_COLORS = {
    **meta_data.get("celltype_colors", {}),
    "Other": "#d0d0d0",
}

GROUP_LABELS = {
    "original": "Original",
    "traceC_increment": "traceC inc.",
    "traceCB_increment": "traceCB inc.",
    "traceC_full": "traceC",
    "traceCB_full": "traceCB",
    "Heritability Significant": "Heritability significant",
    "Heritability Significant EAS": "EAS heritability significant",
    "Heritability Significant EUR": "EUR heritability significant",
    "Correlation Significant": "Correlation Significant",
}

LIBRARY_LABELS = {
    "h.all.v2026.1.Hs.symbols": "Hallmark",
    "c5.go.bp.v2026.1.Hs.symbols": "GO:BP",
    "c2.cp.reactome.v2026.1.Hs.symbols": "Reactome",
    "c7.all.v2026.1.Hs.symbols": "C7",
    "c8.all.v2026.1.Hs.symbols": "C8",
    "hpa_blood_immune_2025": "HPA blood/immune",
    "GTEx_Tissues_V8_2023": "GTEx",
    "GWAS_Catalog_2025": "GWAS Catalog",
    "KEGG_2021_Human": "KEGG",
    "c8_immune_blood_subset": "C8 immune/blood",
    "go_bp_immune_subset": "GO BP immune",
    "reactome_immune_subset": "Reactome immune",
    "gtex_whole_blood_subset": "GTEx Whole Blood",
    "gwas_immune_subset": "GWAS immune/inflammatory",
}


@dataclass(frozen=True)
class LibrarySpec:
    path: Path
    library: str
    label: str
    analysis_tier: str | None = None


@dataclass(frozen=True)
class OraGeneSet:
    qtdid: str
    study: str
    celltype: str
    query_group: str
    query_group_label: str
    gene_symbols: tuple[str, ...]
    background_symbols: tuple[str, ...]
    gene_set_mode: str | None = None


def sanitize_filename(value: object) -> str:
    text = re.sub(r"[^\w.-]+", "_", str(value).strip())
    return text.strip("_") or "result"


def unique_sorted(values: Iterable[object]) -> tuple[str, ...]:
    return tuple(sorted({str(v) for v in values if isinstance(v, str) and v}))


def clip_p(value: object) -> float:
    if pd.isna(value):
        return np.nan
    return float(np.clip(float(value), MIN_P, 1.0))


def build_gene_symbol_map() -> dict[str, str]:
    converter = geneid2name()
    return dict(zip(converter.gtex_df["GENE_ID"], converter.gtex_df["GENE_NAME"]))


def gene_id_to_symbol(gene_id: object, gene_symbol_map: dict[str, str]) -> str | None:
    if not isinstance(gene_id, str) or not gene_id:
        return None
    return gene_symbol_map.get(gene_id.split(".")[0])


def gmt_library_id(path: Path) -> str:
    return path.name.removesuffix(".gmt")


def library_spec(
    path: Path, analysis_tier: str | None = None, label: str | None = None
) -> LibrarySpec:
    library = gmt_library_id(path)
    return LibrarySpec(
        path=path,
        library=library,
        label=label or LIBRARY_LABELS.get(library, library),
        analysis_tier=analysis_tier,
    )


def load_gmt_term_count(path: Path) -> int:
    return len(read_gmt(str(path)))


def write_gmt(path: Path, term_to_genes: dict[str, list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for term in sorted(term_to_genes):
            genes = sorted(set(term_to_genes[term]))
            if genes:
                handle.write(f"{term}\t\t" + "\t".join(genes) + "\n")


def subset_gmt(in_path: Path, out_path: Path, pattern: re.Pattern[str]) -> Path:
    gmt = read_gmt(str(in_path))
    subset = {
        term: genes
        for term, genes in gmt.items()
        if pattern.search(term.replace("_", " "))
    }
    write_gmt(out_path, subset)
    return out_path


def prepare_method_gene_table(
    target_qtdids: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    _, all_df = load_all_summary()
    gene_symbol_map = build_gene_symbol_map()

    df = all_df[all_df["QTDid"].isin(target_qtdids)].copy()
    df["GENE"] = df["GENE"].astype(str).str.split(".").str[0]
    df["gene_symbol"] = df["GENE"].map(
        lambda gene: gene_id_to_symbol(gene, gene_symbol_map)
    )
    df = df[df["gene_symbol"].notna()].copy()

    df["original"] = pd.to_numeric(df["TAR_SeSNP"], errors="coerce").fillna(0) > 0
    df["traceC_full"] = pd.to_numeric(df["TAR_CeSNP"], errors="coerce").fillna(0) > 0
    df["traceCB_full"] = pd.to_numeric(df["TAR_TeSNP"], errors="coerce").fillna(0) > 0
    df["traceC_increment"] = df["traceC_full"] & ~df["original"]
    df["traceCB_increment"] = df["traceCB_full"] & ~df["original"]
    df["study"] = df["QTDid"].map(meta_data["id2name"])
    df["celltype"] = df["QTDid"].map(meta_data["id2celltype"])

    rows = []
    for qtdid, study_df in df.groupby("QTDid", sort=False):
        rows.append(
            {
                "QTDid": qtdid,
                "Study": meta_data["id2name"].get(qtdid, qtdid),
                "CellType": meta_data["id2celltype"].get(qtdid, "Other"),
                "Tested_Genes": study_df["gene_symbol"].nunique(),
                "Original": int(
                    study_df.loc[study_df["original"], "gene_symbol"].nunique()
                ),
                "traceC_inc": int(
                    study_df.loc[study_df["traceC_increment"], "gene_symbol"].nunique()
                ),
                "traceCB_inc": int(
                    study_df.loc[study_df["traceCB_increment"], "gene_symbol"].nunique()
                ),
                "traceC_full": int(
                    study_df.loc[study_df["traceC_full"], "gene_symbol"].nunique()
                ),
                "traceCB_full": int(
                    study_df.loc[study_df["traceCB_full"], "gene_symbol"].nunique()
                ),
            }
        )
    return df, pd.DataFrame(rows)


def build_method_ora_gene_sets(study_df: pd.DataFrame, mode: str) -> list[OraGeneSet]:
    if mode == "incremental":
        groups = ("original", "traceC_increment", "traceCB_increment")
    elif mode == "full":
        groups = ("original", "traceC_full", "traceCB_full")
    else:
        raise ValueError(f"Unsupported gene-set mode: {mode}")

    gene_sets: list[OraGeneSet] = []
    for qtdid, group_df in study_df.groupby("QTDid", sort=False):
        background = unique_sorted(group_df["gene_symbol"])
        for group in groups:
            gene_sets.append(
                OraGeneSet(
                    qtdid=qtdid,
                    study=meta_data["id2name"].get(qtdid, qtdid),
                    celltype=meta_data["id2celltype"].get(qtdid, "Other"),
                    query_group=group,
                    query_group_label=GROUP_LABELS.get(group, group),
                    gene_symbols=unique_sorted(
                        group_df.loc[group_df[group], "gene_symbol"]
                    ),
                    background_symbols=background,
                    gene_set_mode=mode,
                )
            )
    return gene_sets


ANCESTRY_H2_COLUMNS = {
    "EAS": ("H1SQ", "H1SQSE"),
    "EUR": ("H2SQ", "H2SQSE"),
}


def map_gene_ids_to_symbols(
    gene_ids: Iterable[object], gene_symbol_map: dict[str, str]
) -> tuple[str, ...]:
    symbols = []
    for gene_id in gene_ids:
        symbol = gene_id_to_symbol(gene_id, gene_symbol_map)
        if symbol:
            symbols.append(symbol)
    return unique_sorted(symbols)


def filter_heritability_significant(
    summary_df: pd.DataFrame, p_threshold: float, ancestries: tuple[str, ...]
) -> pd.DataFrame:
    z_threshold = p2z(p_threshold)
    mask = pd.Series(True, index=summary_df.index)
    for ancestry in ancestries:
        h2_col, se_col = ANCESTRY_H2_COLUMNS[ancestry]
        h2 = pd.to_numeric(summary_df[h2_col], errors="coerce")
        se = pd.to_numeric(summary_df[se_col], errors="coerce")
        mask &= (h2 / se) > z_threshold
        mask &= h2 > 1e-12
    return summary_df.loc[mask].copy()


def build_significant_ora_gene_sets(
    target_qtdids: list[str],
    h2_p_threshold: float,
) -> tuple[list[OraGeneSet], pd.DataFrame]:
    summary_sign_df_all, summary_df_all = load_all_summary()
    gene_symbol_map = build_gene_symbol_map()
    gene_sets: list[OraGeneSet] = []
    summary_rows = []
    definitions = [
        ("Heritability Significant", ("EAS", "EUR")),
        ("Heritability Significant EAS", ("EAS",)),
        ("Heritability Significant EUR", ("EUR",)),
    ]

    for qtdid in target_qtdids:
        study_df = summary_df_all[summary_df_all["QTDid"] == qtdid].copy()
        corr_df = summary_sign_df_all[summary_sign_df_all["QTDid"] == qtdid].copy()
        background = map_gene_ids_to_symbols(study_df["GENE"], gene_symbol_map)
        study = meta_data["id2name"].get(qtdid, qtdid)
        celltype = meta_data["id2celltype"].get(qtdid, "Other")
        for query_group, ancestries in definitions:
            h2_df = filter_heritability_significant(
                study_df, h2_p_threshold, ancestries
            )
            genes = map_gene_ids_to_symbols(h2_df["GENE"], gene_symbol_map)
            gene_sets.append(
                OraGeneSet(
                    qtdid=qtdid,
                    study=study,
                    celltype=celltype,
                    query_group=query_group,
                    query_group_label=GROUP_LABELS[query_group],
                    gene_symbols=genes,
                    background_symbols=background,
                )
            )
            summary_rows.append(
                {
                    "QTDid": qtdid,
                    "Study": study,
                    "CellType": celltype,
                    "Query_Group": query_group,
                    "Query_Genes": len(genes),
                    "Background_Genes": len(background),
                }
            )
        corr_genes = map_gene_ids_to_symbols(corr_df["GENE"], gene_symbol_map)
        gene_sets.append(
            OraGeneSet(
                qtdid=qtdid,
                study=study,
                celltype=celltype,
                query_group="Correlation Significant",
                query_group_label=GROUP_LABELS["Correlation Significant"],
                gene_symbols=corr_genes,
                background_symbols=background,
            )
        )
        summary_rows.append(
            {
                "QTDid": qtdid,
                "Study": study,
                "CellType": celltype,
                "Query_Group": "Correlation Significant",
                "Query_Genes": len(corr_genes),
                "Background_Genes": len(background),
            }
        )
    return gene_sets, pd.DataFrame(summary_rows)


def filter_gmt_for_background(
    gmt_path: Path,
    background: set[str],
    min_size: int,
    max_size: int,
) -> dict[str, list[str]]:
    gmt = read_gmt(str(gmt_path))
    filtered: dict[str, list[str]] = {}
    for term, genes in gmt.items():
        members = sorted(set(genes) & background)
        if min_size <= len(members) <= max_size:
            filtered[term] = members
    return filtered


def parse_overlap_size(value: object) -> int:
    try:
        return int(str(value).split("/")[0])
    except (TypeError, ValueError):
        return 0


def _run_one_ora(
    gene_set: OraGeneSet,
    lib: LibrarySpec,
    min_set_size: int,
    max_set_size: int,
) -> pd.DataFrame | None:
    query = list(gene_set.gene_symbols)
    background = list(gene_set.background_symbols)
    if not query or not background:
        return None
    filtered_gmt = filter_gmt_for_background(
        lib.path,
        set(background),
        min_set_size,
        max_set_size,
    )
    if not filtered_gmt:
        return None
    try:
        result = gp.enrich(
            gene_list=query,
            gene_sets=filtered_gmt,
            background=background,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
            verbose=False,
        )
    except (LookupError, ValueError):
        return None
    if result.res2d is None or result.res2d.empty:
        return None

    res_df = result.res2d.copy()
    res_df.rename(
        columns={
            "Term": "term",
            "Overlap": "overlap",
            "P-value": "p_value",
            "Adjusted P-value": "adjusted_p_value",
            "Odds Ratio": "odds_ratio",
            "Combined Score": "combined_score",
            "Genes": "overlap_genes",
        },
        inplace=True,
    )
    res_df["analysis_tier"] = lib.analysis_tier
    res_df["library"] = lib.library
    res_df["library_label"] = lib.label
    res_df["QTDid"] = gene_set.qtdid
    res_df["Study"] = gene_set.study
    res_df["CellType"] = gene_set.celltype
    res_df["query_group"] = gene_set.query_group
    res_df["query_group_label"] = gene_set.query_group_label
    res_df["gene_set_mode"] = gene_set.gene_set_mode
    res_df["query_group_size"] = len(query)
    res_df["background_size"] = len(background)
    res_df["overlap_size"] = res_df["overlap"].map(parse_overlap_size)
    keep_cols = [
        "analysis_tier",
        "gene_set_mode",
        "library",
        "library_label",
        "QTDid",
        "Study",
        "CellType",
        "query_group",
        "query_group_label",
        "query_group_size",
        "background_size",
        "term",
        "overlap",
        "overlap_size",
        "p_value",
        "adjusted_p_value",
        "odds_ratio",
        "combined_score",
        "overlap_genes",
    ]
    return res_df[[c for c in keep_cols if c in res_df.columns]]


def run_gseapy_ora(
    gene_sets: list[OraGeneSet],
    libraries: list[LibrarySpec],
    min_set_size: int,
    max_set_size: int,
    workers: int = 8,
) -> pd.DataFrame:
    tasks = [(gene_set, lib) for gene_set in gene_sets for lib in libraries]
    if not tasks:
        return pd.DataFrame()

    rows = []
    workers = max(1, min(int(workers), len(tasks)))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(_run_one_ora, gene_set, lib, min_set_size, max_set_size)
            for gene_set, lib in tasks
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            if row is not None and not row.empty:
                rows.append(row)
            if completed % 50 == 0 or completed == len(tasks):
                print(f"Completed {completed}/{len(tasks)} ORA tasks", flush=True)
    if not rows:
        return pd.DataFrame()
    sort_cols = [
        c
        for c in [
            "analysis_tier",
            "gene_set_mode",
            "query_group",
            "library",
            "adjusted_p_value",
            "p_value",
        ]
        if c in rows[0].columns
    ]
    return pd.concat(rows, ignore_index=True).sort_values(sort_cols)


def study_order_and_celltypes(df: pd.DataFrame) -> tuple[list[str], dict[str, str]]:
    records = {}
    for pos, row in enumerate(df.to_dict("records")):
        study = row.get("Study")
        qtdid = row.get("QTDid")
        if not isinstance(study, str) or not study:
            continue
        if isinstance(qtdid, str) and qtdid in meta_data.get("id2celltype", {}):
            celltype = meta_data["id2celltype"][qtdid]
            rank = QTD_RANK.get(qtdid, pos)
        else:
            celltype = row.get("CellType", "Other")
            rank = pos
        records.setdefault(
            study,
            {
                "celltype": celltype if isinstance(celltype, str) else "Other",
                "rank": rank,
                "first_seen": pos,
            },
        )

    def sort_key(item: tuple[str, dict[str, object]]) -> tuple[int, int, int, str]:
        study, values = item
        celltype = str(values["celltype"])
        return (
            CELLTYPE_RANK.get(celltype, CELLTYPE_RANK["Other"]),
            int(values["rank"]),
            int(values["first_seen"]),
            study,
        )

    ordered = sorted(records.items(), key=sort_key)
    studies = [study for study, _ in ordered]
    celltypes = {study: str(values["celltype"]) for study, values in ordered}
    return studies, celltypes


def celltype_ranges(
    study_order: list[str], study_celltypes: dict[str, str]
) -> list[tuple[str, int, int]]:
    ranges = []
    start = 0
    while start < len(study_order):
        celltype = study_celltypes.get(study_order[start], "Other")
        end = start + 1
        while (
            end < len(study_order)
            and study_celltypes.get(study_order[end], "Other") == celltype
        ):
            end += 1
        ranges.append((celltype, start, end))
        start = end
    return ranges


def add_celltype_annotations_top(
    ax, study_order: list[str], study_celltypes: dict[str, str]
) -> None:
    if len(study_order) <= 1:
        return
    band_transform = ax.get_xaxis_transform() + mtransforms.ScaledTranslation(
        0,
        4 / 72,
        ax.figure.dpi_scale_trans,
    )
    text_transform = ax.get_xaxis_transform() + mtransforms.ScaledTranslation(
        0,
        9 / 72,
        ax.figure.dpi_scale_trans,
    )
    for celltype, start, end in celltype_ranges(study_order, study_celltypes):
        color = CELLTYPE_COLORS.get(celltype, CELLTYPE_COLORS["Other"])
        label = CELLTYPE_LABELS.get(celltype, celltype)
        left = start - 0.5
        right = end - 0.5
        center = (start + end - 1) / 2
        ax.axvspan(left, right, facecolor=color, alpha=0.20, edgecolor="none", zorder=0)
        ax.plot(
            [left + 0.08, right - 0.08],
            [1.0, 1.0],
            transform=band_transform,
            color=color,
            linewidth=4,
            solid_capstyle="butt",
            clip_on=False,
            zorder=3,
        )
        ax.text(
            center,
            1.0,
            label,
            transform=text_transform,
            ha="center",
            va="bottom",
            fontsize=10.5,
            fontweight="bold",
            color="black",
            clip_on=False,
            zorder=4,
        )


def shorten_term(term: object, max_len: int = 50, max_lines: int = 2) -> str:
    text = str(term).replace("_", " ")
    if len(text) <= max_len:
        return text
    shortened = textwrap.shorten(text, width=max_len * max_lines, placeholder="...")
    return "\n".join(
        textwrap.wrap(
            shortened,
            width=max_len,
            break_long_words=True,
            break_on_hyphens=False,
        )
    )


def label_line_total(labels: Iterable[object]) -> int:
    return sum(str(label).count("\n") + 1 for label in labels)


def label_text_width(labels: Iterable[object]) -> int:
    return max(
        (len(line) for label in labels for line in str(label).splitlines()), default=0
    )


def ora_plot_layout(n_terms: int) -> dict[str, object]:
    if n_terms >= 120:
        return {
            "row_height": 0.17,
            "line_height": 0.045,
            "base_height": 3.1,
            "y_fontsize": 5.4,
            "x_fontsize": 8.5,
            "marker_sizes": (8, 80),
            "legend_fontsize": 8.5,
            "legend_title_fontsize": 9.5,
        }
    if n_terms >= 80:
        return {
            "row_height": 0.20,
            "line_height": 0.055,
            "base_height": 3.2,
            "y_fontsize": 6.2,
            "x_fontsize": 9,
            "marker_sizes": (10, 95),
            "legend_fontsize": 9,
            "legend_title_fontsize": 10,
        }
    return {
        "row_height": 0.43,
        "line_height": 0.13,
        "base_height": 3.4,
        "y_fontsize": 10,
        "x_fontsize": 10,
        "marker_sizes": (28, 210),
        "legend_fontsize": 10,
        "legend_title_fontsize": 11,
    }


def scale_marker_areas(
    values: np.ndarray,
    marker_sizes: tuple[float, float],
    value_range: tuple[float, float] | None = None,
) -> np.ndarray:
    """Map overlap counts linearly onto scatter-marker areas (points squared)."""
    values = np.asarray(values, dtype=float)
    min_area, max_area = (float(value) for value in marker_sizes)
    if value_range is None:
        value_min = float(np.nanmin(values))
        value_max = float(np.nanmax(values))
    else:
        value_min, value_max = (float(value) for value in value_range)
    if math.isclose(value_min, value_max):
        return np.full_like(values, (min_area + max_area) / 2, dtype=float)
    return min_area + (values - value_min) * (max_area - min_area) / (
        value_max - value_min
    )


def representative_overlap_values(values: np.ndarray) -> list[int]:
    """Return up to three readable reference counts for a marker-size legend."""
    values = np.asarray(values, dtype=float)
    value_min = float(np.nanmin(values))
    value_max = float(np.nanmax(values))
    if math.isclose(value_min, value_max):
        return [int(round(value_min))]

    locator = MaxNLocator(
        nbins=6,
        integer=True,
        steps=[1, 2, 2.5, 5, 10],
    )
    ticks = sorted(
        {
            int(round(tick))
            for tick in locator.tick_values(value_min, value_max)
            if value_min <= tick <= value_max
        }
    )
    if not ticks:
        ticks = [int(round(value_min)), int(round(value_max))]
    if len(ticks) > 3:
        positions = np.linspace(0, len(ticks) - 1, num=3).round().astype(int)
        ticks = [ticks[position] for position in positions]
    return list(dict.fromkeys(ticks))


def pick_top_terms(
    df: pd.DataFrame,
    group_cols: list[str],
    top_n: int,
    alpha: float,
) -> pd.DataFrame:
    if df.empty:
        return df
    rows = []
    for _, group_df in df.groupby(group_cols, dropna=False, sort=False):
        sig_df = group_df[group_df["adjusted_p_value"] <= alpha]
        source_df = sig_df if not sig_df.empty else group_df
        rows.append(
            source_df.sort_values(["adjusted_p_value", "p_value"]).head(top_n).copy()
        )
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def select_terms_for_ora_plot(
    reference_df: pd.DataFrame,
    alpha: float,
    top_terms_per_library: int,
) -> pd.DataFrame:
    """Select significant study/pathway rows for a split.

    If a split has no significant rows, fall back to top terms per library so
    the corresponding group still has a diagnostic plot. The significant case
    intentionally does not plot non-significant rows from other studies,
    because that makes all study columns look artificially similar.
    """
    sig_df = reference_df[reference_df["adjusted_p_value"] <= alpha].copy()
    if not sig_df.empty:
        return sig_df
    return pick_top_terms(reference_df, ["library"], top_terms_per_library, alpha)


def plot_ora_dotplots(
    ora_df: pd.DataFrame,
    out_dir: Path,
    alpha: float,
    top_terms_per_library: int,
    title_map: dict[tuple[object, ...], str],
    split_cols: list[str],
    filename_prefix: str,
) -> list[Path]:
    if ora_df.empty:
        return []
    out_dir.mkdir(parents=True, exist_ok=True)
    df = ora_df.copy()
    df["adjusted_p_value"] = pd.to_numeric(df["adjusted_p_value"], errors="coerce")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["overlap_size"] = pd.to_numeric(df["overlap_size"], errors="coerce")
    df = df[
        df["adjusted_p_value"].notna() & df["overlap_size"].notna()
    ].copy()
    if df.empty:
        return []

    saved: list[Path] = []
    split_values = df[split_cols].drop_duplicates().to_dict("records")
    for values in split_values:
        reference_mask = pd.Series(True, index=df.index)
        key = []
        for col in split_cols:
            value = values[col]
            key.append(value)
            reference_mask &= df[col].eq(value)
        reference_df = df.loc[reference_mask].copy()
        if reference_df.empty:
            continue
        sub_df = select_terms_for_ora_plot(reference_df, alpha, top_terms_per_library)
        study_order, study_celltypes = study_order_and_celltypes(reference_df)
        sub_df = sub_df[sub_df["Study"].isin(study_order)].copy()
        if sub_df.empty:
            continue
        sub_df["display_term"] = sub_df["term"].map(shorten_term)
        sub_df["neg_log10_fdr"] = -np.log10(
            sub_df["adjusted_p_value"].clip(lower=MIN_P)
        )

        term_order = (
            sub_df.groupby("display_term")["adjusted_p_value"]
            .min()
            .sort_values(ascending=False)
            .index.tolist()
        )
        sub_df["display_term"] = pd.Categorical(
            sub_df["display_term"], categories=term_order, ordered=True
        )
        sub_df["Study"] = pd.Categorical(
            sub_df["Study"], categories=study_order, ordered=True
        )
        sub_df = sub_df.sort_values(["display_term", "Study"])

        total_y_lines = label_line_total(term_order)
        max_study_width = label_text_width(study_order)
        layout = ora_plot_layout(len(term_order))
        height = max(
            6.4,
            float(layout["row_height"]) * len(term_order)
            + float(layout["line_height"]) * total_y_lines
            + float(layout["base_height"]),
        )
        width = max(15.5, 1.03 * len(study_order) + 0.06 * max_study_width + 6.3)

        fig, ax = plt.subplots(figsize=(width, height))
        fdr_scores = sub_df["neg_log10_fdr"].to_numpy(dtype=float)
        overlap_counts = sub_df["overlap_size"].to_numpy(dtype=float)
        marker_areas = scale_marker_areas(
            overlap_counts,
            layout["marker_sizes"],
        )
        fdr_min = float(np.nanmin(fdr_scores))
        fdr_max = float(np.nanmax(fdr_scores))
        if math.isclose(fdr_min, fdr_max):
            fdr_min -= 0.5
            fdr_max += 0.5
        color_norm = mcolors.Normalize(vmin=fdr_min, vmax=fdr_max)
        study_positions = {study: i for i, study in enumerate(study_order)}
        term_positions = {term: i for i, term in enumerate(term_order)}
        scatter = ax.scatter(
            sub_df["Study"].map(study_positions).to_numpy(dtype=float),
            sub_df["display_term"].map(term_positions).to_numpy(dtype=float),
            c=fdr_scores,
            s=marker_areas,
            cmap="viridis",
            norm=color_norm,
            edgecolor="black",
            linewidth=0.25,
            zorder=2,
        )
        add_celltype_annotations_top(ax, study_order, study_celltypes)
        ax.set_xlabel("")
        ax.set_ylabel("")
        title = title_map.get(tuple(key), "ORA pathway enrichment")
        ax.set_title(title, fontsize=15, fontweight="bold", pad=36)
        ax.set_xticks(range(len(study_order)))
        ax.set_xticklabels(study_order)
        ax.tick_params(axis="x", rotation=45, labelsize=layout["x_fontsize"])
        ax.tick_params(axis="y", labelsize=layout["y_fontsize"])
        for label in ax.get_xticklabels():
            label.set_ha("right")
        ax.set_yticks(range(len(term_order)))
        ax.set_yticklabels(term_order)
        ax.set_ylim(len(term_order) - 0.5, -0.5)
        ax.set_xlim(-0.5, len(study_order) - 0.5)

        colorbar_ax = fig.add_axes([0.82, 0.58, 0.014, 0.18])
        colorbar = fig.colorbar(scatter, cax=colorbar_ax)
        colorbar.set_label(
            r"$-\log_{10}(P_{\mathrm{adj}})$",
            fontsize=layout["legend_title_fontsize"],
            labelpad=8,
        )
        colorbar.ax.set_title(
            "FDR-adjusted\n$P$ value",
            fontsize=layout["legend_title_fontsize"],
            pad=7,
        )
        colorbar.ax.tick_params(labelsize=layout["legend_fontsize"])

        overlap_values = representative_overlap_values(overlap_counts)
        overlap_areas = scale_marker_areas(
            np.asarray(overlap_values, dtype=float),
            layout["marker_sizes"],
            value_range=(
                float(np.nanmin(overlap_counts)),
                float(np.nanmax(overlap_counts)),
            ),
        )
        size_handles = [
            Line2D(
                [],
                [],
                linestyle="none",
                marker="o",
                markersize=math.sqrt(area),
                markerfacecolor="#6f6f6f",
                markeredgecolor="black",
                markeredgewidth=0.25,
            )
            for area in overlap_areas
        ]
        ax.legend(
            handles=size_handles,
            labels=[str(value) for value in overlap_values],
            title="Overlapping genes, $n$",
            loc="upper left",
            bbox_to_anchor=(1.04, 0.52),
            frameon=False,
            fontsize=layout["legend_fontsize"],
            title_fontsize=layout["legend_title_fontsize"],
            borderaxespad=0,
            handletextpad=1.0,
            labelspacing=0.8,
        )
        fig.subplots_adjust(top=0.84, right=0.78, bottom=0.18, left=0.34)

        suffix = "_".join(sanitize_filename(v) for v in key)
        out_path = out_dir / f"{filename_prefix}_{suffix}.pdf"
        fig.savefig(out_path, bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
        saved.append(out_path)
    return saved


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def clean_output_dir(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.pdf", "*.png", "*.tsv", "*.json"):
        for path in out_dir.glob(pattern):
            path.unlink()
    readme = out_dir / "README.md"
    if readme.exists():
        readme.unlink()
