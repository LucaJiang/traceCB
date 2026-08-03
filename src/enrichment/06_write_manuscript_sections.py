#!/usr/bin/env python3
"""Write data-driven manuscript text for the EUR-reference S-LDSC rerun."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RESULT_DIR = (
    Path(__file__).resolve().parents[2]
    / "output"
    / "sldsc_gsea_eur_release_matched"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=DEFAULT_RESULT_DIR,
    )
    return parser.parse_args()


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
    result_path = result_dir / "results/summary/master_sldsc_gsea_results.csv"
    validation_path = result_dir / "metadata/validation_summary.json"
    frame = pd.read_csv(result_path)
    if not validation_path.is_file():
        raise FileNotFoundError(validation_path)
    validation = json.loads(validation_path.read_text())
    if validation.get("status") != "PASS":
        raise ValueError(f"Reference-stack validation did not pass: {validation_path}")

    expected_fdr = bh_adjust(frame["Enrichment_p"])
    if "Enrichment_FDR" not in frame:
        raise ValueError(
            f"{result_path} has no Enrichment_FDR column; rerun 04_aggregate_results.py"
        )
    observed_fdr = pd.to_numeric(frame["Enrichment_FDR"], errors="coerce")
    if not np.allclose(
        observed_fdr.to_numpy(),
        expected_fdr.to_numpy(),
        equal_nan=True,
        rtol=1e-12,
        atol=1e-15,
    ):
        raise ValueError(
            f"Enrichment_FDR values in {result_path} do not match global BH adjustment"
        )

    n_tests = len(frame)
    n_nominal = int((frame["Enrichment_p"] < 0.05).sum())
    n_fdr = int((frame["Enrichment_FDR"] < 0.05).sum())
    median_by_annotation = (
        frame.groupby("AnnotationLabel", sort=False)["Enrichment"].median().round(2)
    )
    medians = "; ".join(f"{key}: {value:.2f}" for key, value in median_by_annotation.items())
    strongest = frame.sort_values("Enrichment_p").iloc[0]

    body = f"""# Manuscript revision: EUR ancestry-matched univariate S-LDSC

## Methods

We performed univariate stratified LD score regression (S-LDSC) for European-ancestry Pan-UK Biobank GWAS summary statistics. The biological annotations remained the EAS/BBJ-derived eGene sets used in the original analysis: for each selected eGene, the cis interval was defined as the minimum and maximum tested SNP position in the corresponding study, without changing the eGene membership or interval boundaries. We projected these fixed GRCh37/hg19 intervals onto 1000 Genomes Project Phase 3 European (EUR) SNPs to generate binary custom annotations. We recomputed each custom annotation's LD scores from the 1000 Genomes Phase 3 EUR PLINK reference using a 1-cM window. Each annotation was tested separately in a model containing the 97-category 1000 Genomes Phase 3 EUR baseline-LD v2.2 model and one custom annotation, with overlapping-annotation enrichment enabled. We used EUR HapMap3 non-MHC regression weights, EUR allele-frequency files, and the compatible HapMap3 non-MHC regression SNP list. Because this LDSC implementation requires identical regression-SNP rows across reference prefixes, we row-restricted the official EUR baseline-LD v2.2 LD-score files to that exact common HapMap3 non-MHC list; no baseline LD-score values, annotation definitions, or model categories were changed. We verified GRCh37 coordinate concordance and rsID compatibility among the EAS/BBJ interval-source variants, EUR PLINK panel, baseline-LD files, weights, frequency files, custom annotations, and GWAS summary statistics before regression.

The EUR S-LDSC reference files were obtained from version 4 of the S-LDSC reference-file archive (Zenodo DOI: 10.5281/zenodo.10515792): `1000G_Phase3_baselineLD_v2.2_ldscores.tgz`, `1000G_Phase3_weights_hm3_no_MHC.tgz`, EUR 1000 Genomes Phase 3 PLINK/frequency files, and `hm3_no_MHC.list.txt`. Published MD5 checksums were verified before extraction, and SHA-256 checksums were recorded for the exact local inputs.

## Results

The ancestry-matched rerun evaluated {n_tests} custom annotation–trait models. At nominal P < 0.05, {n_nominal} models showed evidence of enrichment; {n_fdr} remained significant after Benjamini–Hochberg correction across all {n_tests} tests. Median enrichment estimates by annotation family were {medians}. The smallest enrichment P value was observed for {strongest['AnnotationLabel']} in {strongest['StudyLabel']} for {strongest['TraitLabel']} (enrichment = {strongest['Enrichment']:.2f}, 95% CI {strongest['Enrichment_CI95_low']:.2f}–{strongest['Enrichment_CI95_high']:.2f}, P = {strongest['Enrichment_p']:.3g}, FDR = {strongest['Enrichment_FDR']:.3g}). These estimates supersede the earlier ancestry-mismatched S-LDSC results; the underlying EAS/BBJ-derived eGene sets and cis-interval definitions were unchanged.

## Supplementary figure caption

**Supplementary Figure — Enrichment of EAS/BBJ-derived eGene annotations using EUR ancestry-matched S-LDSC references.** Heatmap of univariate S-LDSC enrichment estimates for fixed EAS/BBJ-derived eGene cis-interval annotations projected onto 1000 Genomes Phase 3 EUR SNPs. Custom annotation LD scores were recomputed with the EUR PLINK reference, and each annotation was fitted separately with the 97-category EUR baseline-LD v2.2 model using EUR HapMap3 non-MHC regression weights, EUR allele frequencies, and an identical compatible HapMap3 non-MHC regression-SNP list across reference prefixes. Asterisks denote nominal enrichment P values (*P < 0.05, **P < 0.01, ***P < 0.001). The bottom row gives the mean across studies. Drug allergy is omitted from the displayed preview heatmap but is retained in the complete enrichment tables.

## Reproducibility and data availability draft

European-ancestry Pan-UK Biobank GWAS summary statistics and phenotype documentation are available from the Pan-UK Biobank resource (https://pan.ukbb.broadinstitute.org/). The EUR S-LDSC reference resources are available from Zenodo (version 4; DOI: 10.5281/zenodo.10515792). Analysis code, the trait configuration, and SHA-256 checksums that lock the EAS/BBJ-derived annotation definitions are maintained in the traceCB repository (https://github.com/LucaJiang/traceCB). Generated annotations, validation reports, LD scores, enrichment tables, and figures are written under the configured result directory (`{result_dir}`). **Before publication, replace this sentence with the DOI or accession of the archival deposit containing the derived outputs required by the journal's reproducibility policy.**
"""
    out_dir = result_dir / "manuscript"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "sldsc_eur_reference_revision.md"
    out_path.write_text(body)
    print(f"[done] wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
