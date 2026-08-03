#!/usr/bin/env python3
"""Prepare eGene-interval SNP annotations for S-LDSC.

This pipeline does not directly annotate the discovered eSNP rsID lists. It
first defines method-specific eGene sets from each study's GMM summary table,
then marks every EUR 1000G reference SNP as 1 if its physical position falls
inside the tested cis interval of any selected eGene.

Each custom annotation is written as a separate one-column thin `.annot.gz`
prefix. Downstream S-LDSC runs therefore estimate `baselineLD + one custom
annotation` at a time.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = Path(__file__).with_name("config")
DEFAULT_RESULT_DIR = REPO_ROOT / "output" / "sldsc_gsea_eur_release_matched"
DEFAULT_TRAIT_CONFIG = CONFIG_DIR / "traits.tsv"

STUDIES = (
    "QTD000021",
    "QTD000069",
    "QTD000081",
    "QTD000031",
    "QTD000067",
    "QTD000371",
    "QTD000066",
    "QTD000372",
    "QTD000073",
    "QTD000115",
)

STUDY_LABELS = {
    "QTD000021": "Mono | BLUEPRINT(191)",
    "QTD000069": "Mono | CEDAR(286)",
    "QTD000081": "Mono | Fairfax(420)",
    "QTD000031": "CD4+ T | BLUEPRINT(167)",
    "QTD000067": "CD4+ T | CEDAR(290)",
    "QTD000371": "CD4+ T | Kasela(280)",
    "QTD000066": "CD8+ T | CEDAR(277)",
    "QTD000372": "CD8+ T | Kasela(269)",
    "QTD000073": "B | CEDAR(262)",
    "QTD000115": "NK | Gilchrist(247)",
}

@dataclass(frozen=True)
class AnnotationSpec:
    annot_id: str
    model: str
    study: str
    column: str
    label: str
    genes: frozenset[str]
    intervals_by_chr: dict[int, tuple[tuple[int, int], ...]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--gwas-root", type=Path)
    parser.add_argument(
        "--trait-config", type=Path, default=DEFAULT_TRAIT_CONFIG
    )
    parser.add_argument(
        "--bim-prefix",
        type=Path,
        help="EUR PLINK prefix; defaults to RESULT_DIR/reference/.../1000G.EUR.QC.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_summary(study_dir: Path, study: str) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    for chrom in range(1, 23):
        path = study_dir / study / "GMM" / f"chr{chrom}" / "summary.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(
            path,
            usecols=["GENE", "TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"],
        )
        frame["CHR"] = chrom
        pieces.append(frame)
    out = pd.concat(pieces, ignore_index=True)
    out["GENE"] = out["GENE"].astype(str).str.split(".").str[0]
    for column in ("TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"):
        out[column] = pd.to_numeric(out[column], errors="coerce")
    out["TAR_CeSNP"] = out["TAR_CeSNP"].fillna(out["TAR_SeSNP"])
    out["TAR_TeSNP"] = out["TAR_TeSNP"].fillna(out["TAR_CeSNP"])
    return out


def egenes_from_summary(summary: pd.DataFrame) -> dict[str, set[str]]:
    original = set(summary.loc[summary["TAR_SeSNP"] > 0, "GENE"])
    tracec = set(summary.loc[summary["TAR_CeSNP"] > 0, "GENE"])
    tracecb = set(summary.loc[summary["TAR_TeSNP"] > 0, "GENE"])
    return {
        "original": original,
        "traceC_increment": tracec - original,
        "traceCB_increment": tracecb - original,
        "traceC_overall": tracec,
        "traceCB_overall": tracecb,
    }


def read_gene_intervals_for_study(
    study_dir: Path,
    study: str,
    target_genes: set[str],
) -> dict[str, dict[int, tuple[int, int]]]:
    intervals: dict[str, dict[int, tuple[int, int]]] = {}
    if not target_genes:
        return intervals

    for chrom in range(1, 23):
        path = study_dir / study / "INFO" / f"chr{chrom}.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path, usecols=["GENE", "POS"])
        frame["GENE"] = frame["GENE"].astype(str).str.split(".").str[0]
        frame = frame[frame["GENE"].isin(target_genes)].copy()
        if frame.empty:
            continue
        frame["POS"] = pd.to_numeric(frame["POS"], errors="coerce")
        frame = frame.dropna(subset=["POS"])
        grouped = frame.groupby("GENE", sort=False)["POS"].agg(["min", "max"])
        for gene, row in grouped.iterrows():
            intervals.setdefault(gene, {})[chrom] = (int(row["min"]), int(row["max"]))
    return intervals


def intervals_for_genes(
    genes: set[str],
    interval_source: dict[str, dict[int, tuple[int, int]]],
) -> dict[int, tuple[tuple[int, int], ...]]:
    by_chr: dict[int, list[tuple[int, int]]] = {chrom: [] for chrom in range(1, 23)}
    for gene in genes:
        for chrom, interval in interval_source.get(gene, {}).items():
            by_chr[chrom].append(interval)
    return {
        chrom: tuple(sorted(intervals))
        for chrom, intervals in by_chr.items()
    }


def merge_interval_sources(
    *sources: dict[str, dict[int, tuple[int, int]]]
) -> dict[str, dict[int, tuple[int, int]]]:
    merged: dict[str, dict[int, tuple[int, int]]] = {}
    for source in sources:
        for gene, by_chr in source.items():
            for chrom, interval in by_chr.items():
                if gene not in merged:
                    merged[gene] = {}
                if chrom not in merged[gene]:
                    merged[gene][chrom] = interval
                else:
                    old_start, old_end = merged[gene][chrom]
                    new_start, new_end = interval
                    merged[gene][chrom] = (min(old_start, new_start), max(old_end, new_end))
    return merged


def build_annotation_specs(study_dir: Path) -> list[AnnotationSpec]:
    study_gene_sets: dict[tuple[str, str], set[str]] = {}
    study_intervals: dict[str, dict[str, dict[int, tuple[int, int]]]] = {}

    for study in STUDIES:
        print(f"[intervals] {study}", flush=True)
        summary = read_summary(study_dir, study)
        gene_sets = egenes_from_summary(summary)
        study_gene_sets[(study, "original")] = gene_sets["original"]
        study_gene_sets[(study, "traceC_increment")] = gene_sets["traceC_increment"]
        study_gene_sets[(study, "traceCB_increment")] = gene_sets["traceCB_increment"]
        study_gene_sets[(study, "traceC_overall")] = gene_sets["traceC_overall"]
        study_gene_sets[(study, "traceCB_overall")] = gene_sets["traceCB_overall"]
        all_needed = (
            gene_sets["original"]
            | gene_sets["traceC_increment"]
            | gene_sets["traceCB_increment"]
            | gene_sets["traceC_overall"]
            | gene_sets["traceCB_overall"]
        )
        study_intervals[study] = read_gene_intervals_for_study(study_dir, study, all_needed)

    specs: list[AnnotationSpec] = []
    for study in STUDIES:
        per_study_source = study_intervals[study]
        for order_key, label in (
            ("original", "Original"),
            ("traceC_increment", "traceC increment"),
            ("traceCB_increment", "traceCB increment"),
        ):
            genes = study_gene_sets[(study, order_key)]
            specs.append(
                AnnotationSpec(
                    annot_id=f"incremental_{study}_{order_key}",
                    model="incremental",
                    study=study,
                    column=order_key,
                    label=label,
                    genes=frozenset(genes),
                    intervals_by_chr=intervals_for_genes(genes, per_study_source),
                )
            )

    overall_definitions = (
        ("overall_original", "overall", "overall", "original_overall", "Original overall", "original"),
        ("overall_traceC", "overall", "overall", "traceC_overall", "traceC overall", "traceC_overall"),
        ("overall_traceCB", "overall", "overall", "traceCB_overall", "traceCB overall", "traceCB_overall"),
    )
    for annot_id, model, study_name, column, label, source_key in overall_definitions:
        genes = set().union(*(study_gene_sets[(study, source_key)] for study in STUDIES))
        source = merge_interval_sources(*(study_intervals[study] for study in STUDIES))
        specs.append(
            AnnotationSpec(
                annot_id=annot_id,
                model=model,
                study=study_name,
                column=column,
                label=label,
                genes=frozenset(genes),
                intervals_by_chr=intervals_for_genes(genes, source),
            )
        )
    return specs


def read_bim(bim_prefix: Path, chrom: int) -> pd.DataFrame:
    path = Path(f"{bim_prefix}{chrom}.bim")
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        usecols=[1, 3],
        names=["SNP", "POS"],
        dtype={"SNP": str, "POS": int},
    )


def annotate_positions(positions: np.ndarray, intervals: tuple[tuple[int, int], ...]) -> np.ndarray:
    if not intervals:
        return np.zeros(len(positions), dtype=np.uint8)
    starts = np.fromiter((start for start, _ in intervals), dtype=np.int64)
    ends = np.fromiter((end for _, end in intervals), dtype=np.int64)
    starts.sort()
    ends.sort()
    open_count = np.searchsorted(starts, positions, side="right") - np.searchsorted(
        ends, positions, side="left"
    )
    return (open_count > 0).astype(np.uint8)


def write_trait_manifest(
    result_dir: Path, trait_config: Path, gwas_root: Path | None
) -> None:
    path = result_dir / "metadata" / "trait_manifest.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    traits = pd.read_csv(trait_config, sep="\t", dtype=str)
    required = {
        "Trait",
        "TraitLabel",
        "TraitGroup",
        "GWASAncestry",
        "GenomeBuild",
        "GWASSource",
        "SumstatsRelativePath",
    }
    missing = required.difference(traits.columns)
    if missing:
        raise ValueError(f"{trait_config} missing columns: {sorted(missing)}")
    if traits["Trait"].duplicated().any():
        raise ValueError(f"Duplicate trait IDs in {trait_config}")

    rows: list[dict[str, object]] = []
    for index, trait in traits.iterrows():
        sumstats_path = Path(trait["SumstatsRelativePath"])
        if not sumstats_path.is_absolute():
            if gwas_root is None:
                raise ValueError(
                    "--gwas-root is required when trait paths are relative"
                )
            sumstats_path = gwas_root / sumstats_path
        sumstats_path = sumstats_path.resolve()
        if not sumstats_path.exists():
            raise FileNotFoundError(sumstats_path)
        row = {
            "TraitOrder": index,
            **{column: trait[column] for column in required - {"SumstatsRelativePath"}},
            "SumstatsPath": str(sumstats_path),
        }
        rows.append(row)
    with path.open("w", newline="") as handle:
        fieldnames = [
            "TraitOrder",
            "Trait",
            "TraitLabel",
            "TraitGroup",
            "GWASAncestry",
            "GenomeBuild",
            "GWASSource",
            "SumstatsPath",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_annotation_files(
    result_dir: Path,
    specs: list[AnnotationSpec],
    bim_prefix: Path,
    overwrite: bool,
) -> None:
    annot_root = result_dir / "annotations" / "ldscores"
    manifest_rows: list[dict[str, object]] = []
    bim_by_chr = {chrom: read_bim(bim_prefix, chrom) for chrom in range(1, 23)}

    for spec in specs:
        print(f"[annotation] {spec.annot_id}", flush=True)
        out_dir = annot_root / spec.annot_id
        out_dir.mkdir(parents=True, exist_ok=True)
        annotated_snps: list[str] = []
        annotated_count = 0
        for chrom, bim in bim_by_chr.items():
            out_path = out_dir / f"{spec.annot_id}.{chrom}.annot.gz"
            values = annotate_positions(
                bim["POS"].to_numpy(dtype=np.int64),
                spec.intervals_by_chr.get(chrom, ()),
            )
            annotated_count += int(values.sum())
            if values.any():
                annotated_snps.extend(bim.loc[values.astype(bool), "SNP"].tolist())
            if out_path.exists() and not overwrite:
                continue
            with gzip.open(out_path, "wt", compresslevel=6) as handle:
                handle.write(f"{spec.column}\n")
                handle.write("\n".join(map(str, values.tolist())))
                handle.write("\n")

        set_path = result_dir / "snp_sets" / f"{spec.annot_id}.txt"
        set_path.parent.mkdir(parents=True, exist_ok=True)
        if overwrite or not set_path.exists():
            set_path.write_text("\n".join(sorted(set(annotated_snps))) + "\n")

        interval_count = sum(len(v) for v in spec.intervals_by_chr.values())
        definition_payload = {
            "AnnotID": spec.annot_id,
            "Genes": sorted(spec.genes),
            "IntervalsByChromosome": {
                str(chrom): [list(interval) for interval in intervals]
                for chrom, intervals in sorted(spec.intervals_by_chr.items())
            },
        }
        definition_sha256 = hashlib.sha256(
            json.dumps(definition_payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        manifest_rows.append(
            {
                "AnnotID": spec.annot_id,
                "Model": spec.model,
                "Study": spec.study,
                "StudyLabel": STUDY_LABELS.get(spec.study, "Overall union"),
                "AnnotationOrder": {"original": 0, "traceC_increment": 1, "traceCB_increment": 2,
                                    "original_overall": 0, "traceC_overall": 1, "traceCB_overall": 2}[
                    spec.column
                ],
                "Annotation": spec.column,
                "AnnotationLabel": spec.label,
                "InputGenes": len(spec.genes),
                "InputIntervals": interval_count,
                "DefinitionSHA256": definition_sha256,
                "InputSNPs": annotated_count,
                "AnnotationDerivation": "EAS/BBJ-derived eGene cis intervals",
                "SNPReferencePopulation": "EUR",
                "GenomeBuild": "GRCh37",
                "AnnotPrefix": str(annot_root / spec.annot_id / f"{spec.annot_id}."),
            }
        )

    manifest_path = result_dir / "metadata" / "annotation_manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    gene_set_path = result_dir / "metadata" / "annotation_egene_sets.tsv.gz"
    with gzip.open(gene_set_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["AnnotID", "Gene"], delimiter="\t")
        writer.writeheader()
        for spec in specs:
            for gene in sorted(spec.genes):
                writer.writerow({"AnnotID": spec.annot_id, "Gene": gene})

    interval_path = result_dir / "metadata" / "annotation_cis_intervals.tsv.gz"
    with gzip.open(interval_path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["AnnotID", "CHR", "Start", "End"],
            delimiter="\t",
        )
        writer.writeheader()
        for spec in specs:
            for chrom, intervals in sorted(spec.intervals_by_chr.items()):
                for start, end in intervals:
                    writer.writerow(
                        {
                            "AnnotID": spec.annot_id,
                            "CHR": chrom,
                            "Start": start,
                            "End": end,
                        }
                    )


def write_readme(result_dir: Path) -> None:
    body = """# S-LDSC GSEA-Oriented eGene-Interval SNP Annotation Panel

This directory contains S-LDSC inputs and outputs for eGene-interval SNP
annotations. The annotation unit used by LDSC is a SNP, but SNPs are selected
by eGene intervals rather than by discovered eSNP rsID membership.

Definition:

1. eGenes are called from each study's `GMM/chr*/summary.csv`:
   - `original`: `TAR_SeSNP > 0`.
   - `traceC_increment`: `TAR_CeSNP > 0` and not `original`.
   - `traceCB_increment`: `TAR_TeSNP > 0` and not `original`.
   Missing traceC/traceCB eSNP counts are filled from the previous method,
   matching the plotting/GSEA utilities used elsewhere in the project.
2. For every selected eGene, its cis interval is the min and max `POS` among
   tested SNPs in the corresponding study's `INFO/chr*.csv`.
3. The eGene sets and cis-interval definitions are EAS/BBJ-derived and are not
   changed by the ancestry matching of the S-LDSC reference stack.
4. A 1000 Genomes Phase 3 EUR reference SNP is annotated as 1 if its GRCh37 BIM position overlaps at
   least one selected eGene interval; all other reference SNPs are 0.

Each custom annotation is run separately as `EUR baseline-LD v2.2 + one
custom annotation` with `--overlap-annot`. This avoids estimating original,
traceC_increment, and traceCB_increment in the same custom joint model.

All LD scores, regression weights, allele frequencies, and regression SNPs
used downstream are EUR ancestry-matched. The baseline model is the 97-category
1000 Genomes Phase 3 EUR baseline-LD model v2.2.

Annotation families:

- Study-specific incremental annotations:
  `original`, `traceC_increment`, and `traceCB_increment`, one prefix per
  study and annotation.
- Overall annotations:
  `original_overall`, `traceC_overall`, and `traceCB_overall`, each using the
  union of the corresponding eGenes across all 10 studies.
"""
    (result_dir / "README.md").write_text(body)


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir.resolve()
    study_dir = args.study_dir.resolve()
    result_dir.mkdir(parents=True, exist_ok=True)
    bim_prefix = args.bim_prefix or (
        result_dir
        / "reference/1000G_EUR_Phase3_plink/1000G.EUR.QC."
    )
    bim_prefix = bim_prefix.resolve()
    specs = build_annotation_specs(study_dir)
    write_trait_manifest(result_dir, args.trait_config, args.gwas_root)
    write_annotation_files(result_dir, specs, bim_prefix, args.overwrite)
    write_readme(result_dir)
    print(
        f"[done] wrote {len(specs)} single-column annotations under {result_dir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
