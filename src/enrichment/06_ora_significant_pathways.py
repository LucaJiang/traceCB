"""Run GSEApy ORA for heritability- and correlation-significant genes."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import ora_utils as ou

DEFAULT_OUT_DIR = ou.DEFAULT_RESULT_ROOT / "ora" / "significant_pathways"
DEFAULT_GMT_PATHS = [
    ou.DEFAULT_GMT_DIR / "h.all.v2026.1.Hs.symbols.gmt",
    ou.DEFAULT_GMT_DIR / "c5.go.bp.v2026.1.Hs.symbols.gmt",
    ou.DEFAULT_GMT_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
    ou.DEFAULT_GMT_DIR / "c7.all.v2026.1.Hs.symbols.gmt",
    ou.DEFAULT_GMT_DIR / "c8.all.v2026.1.Hs.symbols.gmt",
    ou.DEFAULT_GMT_DIR / "GTEx_Tissues_V8_2023.gmt",
    ou.DEFAULT_GMT_DIR / "GWAS_Catalog_2025.gmt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--target-qtdids", nargs="+", default=ou.meta_data["QTDids"])
    parser.add_argument("--gmt", nargs="+", default=[str(path) for path in DEFAULT_GMT_PATHS])
    parser.add_argument("--h2-p-threshold", type=float, default=0.05)
    parser.add_argument("--min-set-size", type=int, default=5)
    parser.add_argument("--max-set-size", type=int, default=1000)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--top-terms-per-library", type=int, default=4)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--no-clean", dest="clean", action="store_false")
    parser.set_defaults(clean=True)
    return parser.parse_args()


def write_readme(
    out_dir: Path,
    args: argparse.Namespace,
    gene_set_summary: pd.DataFrame,
    ora_df: pd.DataFrame,
    saved_figures: list[Path],
) -> None:
    fdr_rows = int((ora_df["adjusted_p_value"] <= args.alpha).sum()) if not ora_df.empty else 0
    significant_terms = (
        ora_df.loc[ora_df["adjusted_p_value"] <= args.alpha, "term"].nunique()
        if not ora_df.empty
        else 0
    )
    query_summary = (
        gene_set_summary.groupby("Query_Group", as_index=False)["Query_Genes"]
        .agg(["min", "median", "max"])
        .reset_index()
    )
    query_lines = []
    for row in query_summary.itertuples(index=False):
        query_lines.append(
            f"| {row.Query_Group} | {int(row.min)} | {int(row.median)} | {int(row.max)} |"
        )
    figure_lines = [f"- `{path.name}`" for path in saved_figures]
    body = f"""# ORA: Significant Pathways

This directory contains local-GMT ORA results generated with `gseapy.enrich()`.
The tested-gene background is study-specific: all genes tested in the same
eQTL study are used as the ORA background for that study.

## Query Gene Definitions

| Query group | Definition |
|---|---|
| Heritability Significant | Genes with significant positive EAS and EUR local heritability estimates at P < {args.h2_p_threshold} |
| Heritability Significant EAS | Genes with significant positive EAS local heritability estimates at P < {args.h2_p_threshold} |
| Heritability Significant EUR | Genes with significant positive EUR local heritability estimates at P < {args.h2_p_threshold} |
| Correlation Significant | Genes in the traceCB correlation-significant summary table |

## Gene Counts Across Studies

| Query group | Min genes | Median genes | Max genes |
|---|---:|---:|---:|
{chr(10).join(query_lines)}

## Result Summary

- ORA result rows: {len(ora_df):,}
- ORA rows with FDR <= {args.alpha}: {fdr_rows:,}
- Unique pathways with at least one study/query FDR <= {args.alpha}: {significant_terms:,}

## Figures

{chr(10).join(figure_lines)}

## Files

| File | Description |
|---|---|
| `significant_gene_sets.tsv` | Per-study query gene counts and backgrounds |
| `pathway_enrichment_results.tsv` | Full GSEApy ORA table |
| `significant_pathways.tsv` | ORA rows with adjusted P-value/FDR <= {args.alpha} |
| `top_pathways_by_gene_set.tsv` | Top terms per query group, library, and study |
| `ora_significant_*.pdf` | Publication-style ORA dotplots with cell-type labels above the study axis |
"""
    (out_dir / "README.md").write_text(body)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if args.clean:
        ou.clean_output_dir(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ou.write_json(out_dir / "run_config.json", vars(args))

    gene_sets, gene_set_summary = ou.build_significant_ora_gene_sets(
        target_qtdids=args.target_qtdids,
        h2_p_threshold=args.h2_p_threshold,
    )
    gene_set_summary.to_csv(out_dir / "significant_gene_sets.tsv", sep="\t", index=False)

    libraries = [ou.library_spec(Path(path)) for path in args.gmt]
    ora_df = ou.run_gseapy_ora(
        gene_sets=gene_sets,
        libraries=libraries,
        min_set_size=args.min_set_size,
        max_set_size=args.max_set_size,
        workers=args.workers,
    )
    ora_df.to_csv(out_dir / "pathway_enrichment_results.tsv", sep="\t", index=False)
    if not ora_df.empty:
        ora_df[ora_df["adjusted_p_value"] <= args.alpha].to_csv(
            out_dir / "significant_pathways.tsv",
            sep="\t",
            index=False,
        )
        (
            ora_df.sort_values(["query_group", "library", "adjusted_p_value", "p_value"])
            .groupby(["query_group", "library", "QTDid"], as_index=False)
            .head(20)
            .to_csv(out_dir / "top_pathways_by_gene_set.tsv", sep="\t", index=False)
        )

    title_map = {
        (group,): f"ORA of {ou.GROUP_LABELS.get(group, group)} Genes"
        for group in sorted(ora_df["query_group"].dropna().unique())
    } if not ora_df.empty else {}
    saved_figures = ou.plot_ora_dotplots(
        ora_df=ora_df,
        out_dir=out_dir,
        alpha=args.alpha,
        top_terms_per_library=args.top_terms_per_library,
        title_map=title_map,
        split_cols=["query_group"],
        filename_prefix="ora_significant",
    )
    write_readme(out_dir, args, gene_set_summary, ora_df, saved_figures)
    print(f"[done] significant-pathway ORA written to {out_dir}")


if __name__ == "__main__":
    main()
