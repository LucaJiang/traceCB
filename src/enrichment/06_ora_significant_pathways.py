"""Run GSEApy ORA for heritability- and correlation-significant genes."""

from __future__ import annotations

import argparse
from pathlib import Path
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
    parser.add_argument(
        "--gmt", nargs="+", default=[str(path) for path in DEFAULT_GMT_PATHS]
    )
    parser.add_argument("--h2-p-threshold", type=float, default=0.05)
    parser.add_argument("--min-set-size", type=int, default=5)
    parser.add_argument("--max-set-size", type=int, default=1000)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--top-terms-per-library", type=int, default=4)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--no-clean", dest="clean", action="store_false")
    parser.set_defaults(clean=True)
    return parser.parse_args()


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
    gene_set_summary.to_csv(
        out_dir / "significant_gene_sets.tsv", sep="\t", index=False
    )

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
            ora_df.sort_values(
                ["query_group", "library", "adjusted_p_value", "p_value"]
            )
            .groupby(["query_group", "library", "QTDid"], as_index=False)
            .head(20)
            .to_csv(out_dir / "top_pathways_by_gene_set.tsv", sep="\t", index=False)
        )

    title_map = (
        {
            (group,): f"ORA of {ou.GROUP_LABELS.get(group, group)} Genes"
            for group in sorted(ora_df["query_group"].dropna().unique())
        }
        if not ora_df.empty
        else {}
    )
    _ = ou.plot_ora_dotplots(
        ora_df=ora_df,
        out_dir=out_dir,
        alpha=args.alpha,
        top_terms_per_library=args.top_terms_per_library,
        title_map=title_map,
        split_cols=["query_group"],
        filename_prefix="ora_significant",
    )


if __name__ == "__main__":
    main()
