"""Run two-tier GSEApy ORA for original, traceC, and traceCB eGene groups."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

import ora_utils as ou

DEFAULT_OUT_DIR = ou.DEFAULT_RESULT_ROOT / "ora" / "two_tier_enrichment"

IMMUNE_TERM_PATTERN = re.compile(
    r"(immune|immunolog|inflam|hemat|hemato|blood|leuko|lymph|mono|macro|"
    r"neutro|eosin|baso|eryth|platelet|megakaryo|myeloid|b[_ -]?cell|"
    r"t[_ -]?cell|nk[_ -]?cell|dendritic|interferon|interleukin|cytokine|"
    r"antigen|antibody|phagocyt|granul|mast[_ -]?cell|plasma[_ -]?cell|"
    r"reg[_ -]?t|treg|th1|th2|th17|bcr|tcr|innate|adaptive|autoimmune)",
    re.IGNORECASE,
)
WHOLE_BLOOD_PATTERN = re.compile(r"whole blood", re.IGNORECASE)
GWAS_IMMUNE_PATTERN = re.compile(
    r"(immune|inflamm|autoimmune|hemat|blood|anemia|lymph|leuk|myelo|"
    r"monocyte|neutroph|eosinoph|basoph|platelet|erythro|asthma|allergy|"
    r"allergic|eczema|dermatitis|psoriasis|crohn|ulcerative colitis|ibd|"
    r"lupus|multiple sclerosis|rheumatoid|ankylosing|celiac|thyroiditis|"
    r"sjogren|vasculitis|scleroderma|gout)",
    re.IGNORECASE,
)

TIER_LABELS = {
    "blood_immune_focused": "Blood/Immune Focused ORA",
    "broad_pathway_sensitivity": "Broad Pathway Sensitivity ORA",
}
MODE_LABELS = {
    "incremental": "Incremental eGenes",
    "full": "Full eGenes",
}
QUERY_GROUP_ORDER = {
    "original": 0,
    "traceC_increment": 1,
    "traceCB_increment": 2,
    "traceC_full": 1,
    "traceCB_full": 2,
}


def prepare_tier_libraries(out_dir: Path) -> tuple[list[ou.LibrarySpec], pd.DataFrame]:
    subset_dir = out_dir / "gmt_subsets"
    subset_dir.mkdir(parents=True, exist_ok=True)

    specs = [
        ou.library_spec(
            ou.DEFAULT_GMT_DIR / "c7.all.v2026.1.Hs.symbols.gmt",
            analysis_tier="blood_immune_focused",
            label="MSigDB C7 immunologic signatures",
        ),
        ou.library_spec(
            ou.subset_gmt(
                ou.DEFAULT_GMT_DIR / "c8.all.v2026.1.Hs.symbols.gmt",
                subset_dir / "c8_immune_blood_subset.gmt",
                IMMUNE_TERM_PATTERN,
            ),
            analysis_tier="blood_immune_focused",
            label="MSigDB C8 immune/blood subset",
        ),
        ou.library_spec(
            ou.subset_gmt(
                ou.DEFAULT_GMT_DIR / "c5.go.bp.v2026.1.Hs.symbols.gmt",
                subset_dir / "go_bp_immune_subset.gmt",
                IMMUNE_TERM_PATTERN,
            ),
            analysis_tier="blood_immune_focused",
            label="GO BP immune subset",
        ),
        ou.library_spec(
            ou.subset_gmt(
                ou.DEFAULT_GMT_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
                subset_dir / "reactome_immune_subset.gmt",
                IMMUNE_TERM_PATTERN,
            ),
            analysis_tier="blood_immune_focused",
            label="Reactome immune subset",
        ),
        ou.library_spec(
            ou.subset_gmt(
                ou.DEFAULT_GMT_DIR / "GTEx_Tissues_V8_2023.gmt",
                subset_dir / "gtex_whole_blood_subset.gmt",
                WHOLE_BLOOD_PATTERN,
            ),
            analysis_tier="blood_immune_focused",
            label="GTEx Whole Blood subset",
        ),
        ou.library_spec(
            ou.DEFAULT_GMT_DIR / "hpa_blood_immune_2025.gmt",
            analysis_tier="blood_immune_focused",
            label="HPA blood/immune",
        ),
        ou.library_spec(
            ou.subset_gmt(
                ou.DEFAULT_GMT_DIR / "GWAS_Catalog_2025.gmt",
                subset_dir / "gwas_immune_subset.gmt",
                GWAS_IMMUNE_PATTERN,
            ),
            analysis_tier="blood_immune_focused",
            label="GWAS immune/inflammatory subset",
        ),
        ou.library_spec(
            ou.DEFAULT_GMT_DIR / "h.all.v2026.1.Hs.symbols.gmt",
            analysis_tier="broad_pathway_sensitivity",
            label="Hallmark",
        ),
        ou.library_spec(
            ou.DEFAULT_GMT_DIR / "c5.go.bp.v2026.1.Hs.symbols.gmt",
            analysis_tier="broad_pathway_sensitivity",
            label="GO BP all",
        ),
        ou.library_spec(
            ou.DEFAULT_GMT_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
            analysis_tier="broad_pathway_sensitivity",
            label="Reactome all",
        ),
    ]

    rows = []
    for spec in specs:
        rows.append(
            {
                "analysis_tier": spec.analysis_tier,
                "library": spec.library,
                "library_label": spec.label,
                "gmt_path": str(spec.path),
                "pathway_count": ou.load_gmt_term_count(spec.path),
            }
        )
    return specs, pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--target-qtdids", nargs="+", default=ou.meta_data["QTDids"])
    parser.add_argument("--min-set-size", type=int, default=5)
    parser.add_argument("--max-set-size", type=int, default=1000)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--top-terms-per-library", type=int, default=3)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--no-clean", dest="clean", action="store_false")
    parser.set_defaults(clean=True)
    return parser.parse_args()


def build_tier_summary(ora_df: pd.DataFrame, library_summary: pd.DataFrame, alpha: float) -> pd.DataFrame:
    total_terms = (
        library_summary.groupby("analysis_tier", as_index=False)["pathway_count"]
        .sum()
        .rename(columns={"pathway_count": "input_pathways"})
    )
    if ora_df.empty:
        total_terms["gene_set_mode"] = ""
        total_terms["query_group"] = ""
        total_terms["query_group_label"] = ""
        total_terms["ora_rows"] = 0
        total_terms["tested_pathways_after_background_filter"] = 0
        total_terms["significant_pathways_any_study"] = 0
        total_terms["significant_ora_rows"] = 0
        return total_terms
    rows = []
    group_cols = ["analysis_tier", "gene_set_mode", "query_group", "query_group_label"]
    for (tier, mode, query_group, query_label), sub_df in ora_df.groupby(group_cols, dropna=False):
        sig = sub_df[sub_df["adjusted_p_value"] <= alpha]
        rows.append(
            {
                "analysis_tier": tier,
                "gene_set_mode": mode,
                "query_group": query_group,
                "query_group_label": query_label,
                "ora_rows": len(sub_df),
                "tested_pathways_after_background_filter": sub_df["term"].nunique(),
                "significant_pathways_any_study": sig["term"].nunique(),
                "significant_ora_rows": len(sig),
            }
        )
    summary = pd.DataFrame(rows)
    return summary.merge(total_terms, on="analysis_tier", how="left")


def write_readme(
    out_dir: Path,
    args: argparse.Namespace,
    count_df: pd.DataFrame,
    library_summary: pd.DataFrame,
    tier_summary: pd.DataFrame,
    ora_df: pd.DataFrame,
    saved_figures: list[Path],
) -> None:
    tier_count = library_summary.groupby("analysis_tier", as_index=False)["pathway_count"].sum()
    tier_lines = [
        f"| {TIER_LABELS.get(row.analysis_tier, row.analysis_tier)} | {int(row.pathway_count):,} |"
        for row in tier_count.itertuples(index=False)
    ]
    significant_lines = []
    sort_df = tier_summary.copy()
    sort_df["query_group_rank"] = sort_df["query_group"].map(QUERY_GROUP_ORDER).fillna(99)
    for row in sort_df.sort_values(["analysis_tier", "gene_set_mode", "query_group_rank"]).itertuples(index=False):
        significant_lines.append(
            "| "
            f"{TIER_LABELS.get(row.analysis_tier, row.analysis_tier)} | "
            f"{row.gene_set_mode} | "
            f"{row.query_group_label} | "
            f"{int(row.input_pathways):,} | "
            f"{int(row.tested_pathways_after_background_filter):,} | "
            f"{int(row.significant_pathways_any_study):,} | "
            f"{int(row.significant_ora_rows):,} |"
        )
    figure_lines = [f"- `{path.name}`" for path in saved_figures]

    body = f"""# ORA: Two-Tier Enrichment

This directory contains local-GMT ORA results generated with `gseapy.enrich()`.
Each ORA test uses the matching study's tested eGenes as the background.

## eGene Group Definitions

| Mode | Group | Definition |
|---|---|---|
| incremental | Original | eGenes discovered by the original single-ancestry eQTL method |
| incremental | traceC inc. | traceC eGenes minus original eGenes |
| incremental | traceCB inc. | traceCB eGenes minus original eGenes |
| full | Original | eGenes discovered by the original single-ancestry eQTL method |
| full | traceC | All traceC eGenes |
| full | traceCB | All traceCB eGenes |

## Tier Definitions

`blood_immune_focused` is treated as the immune-focused tier. It contains full
MSigDB C7 immunologic signatures plus blood/immune term subsets from C8, GO BP,
Reactome, GTEx Whole Blood, HPA blood/immune, and GWAS Catalog. The subsets are
selected by explicit case-insensitive term-name patterns for immune, blood,
hematopoietic, leukocyte, cytokine, antigen, antibody, and canonical immune
disease keywords. C7 and HPA blood/immune are included directly because those
libraries are already immune/blood focused.

`broad_pathway_sensitivity` is the broad sensitivity tier. It contains full
Hallmark, full GO BP, and full Reactome libraries, without immune-term filtering.

## Pathway Counts by Tier

| Tier | Input pathways before study-background filtering |
|---|---:|
{chr(10).join(tier_lines)}

## Significant Pathway Counts

Significant means at least one study has ORA FDR <= {args.alpha} for the
specified tier, mode, and eGene group.

| Tier | Mode | eGene group | Input pathways | Tested pathways after background filter | Significant pathways in at least one study | Significant ORA rows |
|---|---|---|---:|---:|---:|---:|
{chr(10).join(significant_lines)}

## Study Gene Counts

`group_size_summary.tsv` reports per-study query sizes. The same study-specific
tested-gene universe is used as background for all groups within that study.

## Result Summary

- ORA result rows: {len(ora_df):,}
- ORA rows with FDR <= {args.alpha}: {int((ora_df["adjusted_p_value"] <= args.alpha).sum()) if not ora_df.empty else 0:,}

## Figures

{chr(10).join(figure_lines)}

## Files

| File | Description |
|---|---|
| `group_size_summary.tsv` | Per-study tested genes and original/traceC/traceCB group sizes |
| `tier_library_summary.tsv` | GMT libraries used in each tier and their pathway counts |
| `tier_pathway_summary.tsv` | Tier-, mode-, and eGene-group-level tested and significant pathway counts |
| `ora_results.tsv` | Full GSEApy ORA table |
| `ora_top_terms.tsv` | Top ORA terms per tier, mode, group, library, and study |
| `ora_two_tier_*.pdf` | Publication-style ORA dotplots; each plot includes all significant study-pathway ORA rows for that tier/mode/group, with a top-term fallback only when the group has no significant rows |
"""
    (out_dir / "README.md").write_text(body)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if args.clean:
        ou.clean_output_dir(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ou.write_json(out_dir / "run_config.json", vars(args))

    study_df, count_df = ou.prepare_method_gene_table(args.target_qtdids)
    count_df.to_csv(out_dir / "group_size_summary.tsv", sep="\t", index=False)

    metric_cols = [
        "QTDid",
        "study",
        "celltype",
        "GENE",
        "gene_symbol",
        "original",
        "traceC_increment",
        "traceCB_increment",
        "traceC_full",
        "traceCB_full",
    ]
    study_df[[c for c in metric_cols if c in study_df.columns]].to_csv(
        out_dir / "gene_level_method_groups.tsv",
        sep="\t",
        index=False,
    )

    libraries, library_summary = prepare_tier_libraries(out_dir)
    library_summary.to_csv(out_dir / "tier_library_summary.tsv", sep="\t", index=False)

    gene_sets = (
        ou.build_method_ora_gene_sets(study_df, mode="incremental")
        + ou.build_method_ora_gene_sets(study_df, mode="full")
    )
    ora_df = ou.run_gseapy_ora(
        gene_sets=gene_sets,
        libraries=libraries,
        min_set_size=args.min_set_size,
        max_set_size=args.max_set_size,
        workers=args.workers,
    )
    ora_df.to_csv(out_dir / "ora_results.tsv", sep="\t", index=False)
    if not ora_df.empty:
        (
            ora_df.sort_values(
                [
                    "analysis_tier",
                    "gene_set_mode",
                    "query_group",
                    "library",
                    "adjusted_p_value",
                    "p_value",
                ]
            )
            .groupby(
                ["analysis_tier", "gene_set_mode", "query_group", "library", "QTDid"],
                as_index=False,
            )
            .head(20)
            .to_csv(out_dir / "ora_top_terms.tsv", sep="\t", index=False)
        )
    tier_summary = build_tier_summary(ora_df, library_summary, args.alpha)
    tier_summary.to_csv(out_dir / "tier_pathway_summary.tsv", sep="\t", index=False)

    title_map = {}
    if not ora_df.empty:
        for row in ora_df[["analysis_tier", "gene_set_mode", "query_group", "query_group_label"]].drop_duplicates().itertuples(index=False):
            title_map[(row.analysis_tier, row.gene_set_mode, row.query_group)] = (
                f"{TIER_LABELS.get(row.analysis_tier, row.analysis_tier)}: "
                f"{row.query_group_label} {MODE_LABELS.get(row.gene_set_mode, row.gene_set_mode)}"
            )
    saved_figures = ou.plot_ora_dotplots(
        ora_df=ora_df,
        out_dir=out_dir,
        alpha=args.alpha,
        top_terms_per_library=args.top_terms_per_library,
        title_map=title_map,
        split_cols=["analysis_tier", "gene_set_mode", "query_group"],
        filename_prefix="ora_two_tier",
    )
    write_readme(out_dir, args, count_df, library_summary, tier_summary, ora_df, saved_figures)
    print(f"[done] two-tier ORA written to {out_dir}")


if __name__ == "__main__":
    main()
