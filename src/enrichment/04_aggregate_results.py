#!/usr/bin/env python3
"""Aggregate S-LDSC results for the sldsc_gsea panel."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULT_DIR = Path(
    os.environ.get("TRACECB_ENRICHMENT_DIR", REPO_ROOT / "results/enrichment")
)

STUDY_ORDER = {
    "QTD000021": 0,
    "QTD000069": 1,
    "QTD000081": 2,
    "QTD000031": 3,
    "QTD000067": 4,
    "QTD000371": 5,
    "QTD000066": 6,
    "QTD000372": 7,
    "QTD000073": 8,
    "QTD000115": 9,
    "overall": 99,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    return parser.parse_args()


def read_result(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    required = {
        "Category",
        "Prop._SNPs",
        "Prop._h2",
        "Prop._h2_std_error",
        "Enrichment",
        "Enrichment_std_error",
        "Enrichment_p",
        "Coefficient",
        "Coefficient_std_error",
        "Coefficient_z-score",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    return frame


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    out = np.full_like(p, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    ranked = p[order] * len(valid) / np.arange(1, len(valid) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return pd.Series(out, index=values.index)


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir
    trait_manifest = pd.read_csv(result_dir / "metadata" / "trait_manifest.tsv", sep="\t")
    annot_manifest = pd.read_csv(result_dir / "metadata" / "annotation_manifest.tsv", sep="\t")
    raw_dir = result_dir / "results" / "raw"
    out_dir = result_dir / "results" / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    numeric_columns = [
        "Prop._SNPs",
        "Prop._h2",
        "Prop._h2_std_error",
        "Enrichment",
        "Enrichment_std_error",
        "Enrichment_p",
        "Coefficient",
        "Coefficient_std_error",
        "Coefficient_z-score",
    ]

    for _, trait in trait_manifest.sort_values("TraitOrder").iterrows():
        for annot_id, annot_df in annot_manifest.groupby("AnnotID", sort=False):
            annot_df = annot_df.sort_values("AnnotationOrder").reset_index(drop=True)
            path = raw_dir / f"{trait['Trait']}__{annot_id}.results"
            try:
                frame = numeric(read_result(path), numeric_columns)
                custom = frame.tail(len(annot_df)).reset_index(drop=True)
                if len(custom) != len(annot_df):
                    raise ValueError(f"{path}: expected {len(annot_df)} custom rows, found {len(custom)}")
                for idx, custom_row in custom.iterrows():
                    annot = annot_df.iloc[idx]
                    rows.append(
                        {
                            "Trait": trait["Trait"],
                            "TraitLabel": trait["TraitLabel"],
                            "TraitGroup": trait["TraitGroup"],
                            "TraitOrder": trait["TraitOrder"],
                            "AnnotID": annot["AnnotID"],
                            "Model": annot["Model"],
                            "Study": annot["Study"],
                            "StudyLabel": annot["StudyLabel"],
                            "StudyOrder": STUDY_ORDER.get(str(annot["Study"]), 98),
                            "AnnotationOrder": annot["AnnotationOrder"],
                            "Annotation": annot["Annotation"],
                            "AnnotationLabel": annot["AnnotationLabel"],
                            "InputGenes": annot.get("InputGenes", np.nan),
                            "InputIntervals": annot.get("InputIntervals", np.nan),
                            "InputSNPs": annot["InputSNPs"],
                            "LDSC_Category": custom_row["Category"],
                            "Prop._SNPs": custom_row["Prop._SNPs"],
                            "Prop._h2": custom_row["Prop._h2"],
                            "Prop._h2_std_error": custom_row["Prop._h2_std_error"],
                            "Enrichment": custom_row["Enrichment"],
                            "Enrichment_std_error": custom_row["Enrichment_std_error"],
                            "Enrichment_p": custom_row["Enrichment_p"],
                            "Enrichment_CI95_low": custom_row["Enrichment"]
                            - 1.96 * custom_row["Enrichment_std_error"],
                            "Enrichment_CI95_high": custom_row["Enrichment"]
                            + 1.96 * custom_row["Enrichment_std_error"],
                            "Tau_Coefficient": custom_row["Coefficient"],
                            "Tau_Coefficient_SE": custom_row["Coefficient_std_error"],
                            "Tau_Z_score": custom_row["Coefficient_z-score"],
                            "source_file": str(path),
                        }
                    )
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{path}: {exc}")

    if errors:
        error_path = out_dir / "aggregate_errors.txt"
        error_path.write_text("\n".join(errors) + "\n")
        raise RuntimeError(f"{len(errors)} result files failed to parse; see {error_path}")

    result = pd.DataFrame(rows)
    result = result.sort_values(
        ["Model", "TraitOrder", "StudyOrder", "AnnotationOrder"], kind="stable"
    ).reset_index(drop=True)
    result["Enrichment_FDR"] = bh_adjust(result["Enrichment_p"])
    result.to_csv(out_dir / "master_sldsc_gsea_results.csv", index=False)
    result[result["Model"] == "incremental"].to_csv(
        out_dir / "incremental_sldsc_results.csv", index=False
    )
    result[result["Model"] == "overall"].to_csv(
        out_dir / "overall_sldsc_results.csv", index=False
    )

    summary = (
        result.groupby(["Model", "TraitGroup", "AnnotationLabel"], observed=True)
        .agg(
            n=("Enrichment", "size"),
            mean_enrichment=("Enrichment", "mean"),
            median_enrichment=("Enrichment", "median"),
            mean_tau_z=("Tau_Z_score", "mean"),
            n_nominal_enrichment_p_lt_0_05=("Enrichment_p", lambda x: int((x < 0.05).sum())),
            n_enrichment_fdr_lt_0_05=("Enrichment_FDR", lambda x: int((x < 0.05).sum())),
            mean_prop_snps=("Prop._SNPs", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "sldsc_gsea_group_summary.csv", index=False)
    print(f"[done] aggregated {len(result)} custom annotation rows", flush=True)


if __name__ == "__main__":
    main()
