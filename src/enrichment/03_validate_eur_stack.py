#!/usr/bin/env python3
"""Step 03: validate ancestry, build, identifiers, and EUR cache isolation."""

from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path

import pandas as pd


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

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULT_DIR = Path(
    os.environ.get("TRACECB_ENRICHMENT_DIR", REPO_ROOT / "results/enrichment")
)
DEFAULT_EXPECTED_DEFINITIONS = (
    Path(__file__).with_name("config") / "annotation_definition_checksums.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=DEFAULT_RESULT_DIR,
    )
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument(
        "--bfile-prefix",
        type=Path,
        help="EUR PLINK prefix; defaults under RESULT_DIR/reference.",
    )
    parser.add_argument(
        "--frq-prefix",
        type=Path,
        help="EUR allele-frequency prefix; defaults under RESULT_DIR/reference.",
    )
    parser.add_argument("--check-custom-ldscores", action="store_true")
    parser.add_argument(
        "--expected-definitions",
        type=Path,
        default=DEFAULT_EXPECTED_DEFINITIONS,
        help="Tracked annotation-definition checksums used to detect drift.",
    )
    return parser.parse_args()


def require(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def prefix_path(prefix: Path, chrom: int, suffix: str) -> Path:
    return Path(f"{prefix}{chrom}{suffix}")


def selected_genes(study_dir: Path, study: str) -> set[str]:
    pieces = []
    for chrom in range(1, 23):
        path = study_dir / study / "GMM" / f"chr{chrom}" / "summary.csv"
        pieces.append(
            pd.read_csv(
                require(path),
                usecols=["GENE", "TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"],
            )
        )
    frame = pd.concat(pieces, ignore_index=True)
    frame["GENE"] = frame["GENE"].astype(str).str.split(".").str[0]
    for column in ("TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["TAR_CeSNP"] = frame["TAR_CeSNP"].fillna(frame["TAR_SeSNP"])
    frame["TAR_TeSNP"] = frame["TAR_TeSNP"].fillna(frame["TAR_CeSNP"])
    keep = (frame[["TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"]] > 0).any(axis=1)
    return set(frame.loc[keep, "GENE"])


def read_snp_set(path: Path) -> set[str]:
    with require(path).open() as handle:
        values = [line.strip().split()[0] for line in handle if line.strip()]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate SNP identifiers in {path}")
    if values and values[0].upper() == "SNP":
        raise ValueError(f"Unexpected header in one-column SNP list {path}")
    return set(values)


def compare_annotation_definitions(current: pd.DataFrame, expected_path: Path) -> int:
    expected = pd.read_csv(require(expected_path), sep="\t", dtype={"DefinitionSHA256": str})
    columns = [
        "AnnotID",
        "InputGenes",
        "InputIntervals",
        "DefinitionSHA256",
    ]
    missing_current = set(columns).difference(current.columns)
    missing_expected = set(columns).difference(expected.columns)
    if missing_current or missing_expected:
        raise ValueError(
            "Definition checksum columns are missing: "
            f"current={sorted(missing_current)}, expected={sorted(missing_expected)}"
        )
    left = current[columns].sort_values("AnnotID").reset_index(drop=True)
    right = expected[columns].sort_values("AnnotID").reset_index(drop=True)
    if not left.equals(right):
        comparison = left.merge(
            right,
            on="AnnotID",
            how="outer",
            suffixes=("_observed", "_expected"),
        )
        raise ValueError(
            "EAS/BBJ eGene membership or cis-interval definitions changed:\n"
            + comparison.to_string(index=False)
        )
    return len(left)


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir.resolve()
    bfile_prefix = args.bfile_prefix or (
        result_dir / "reference/1000G_EUR_Phase3_plink/1000G.EUR.QC."
    )
    frq_prefix = args.frq_prefix or (
        result_dir / "reference/1000G_Phase3_frq/1000G.EUR.QC."
    )
    metadata_dir = result_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    if "sldsc_gsea_eur" not in result_dir.name.lower():
        raise ValueError(f"EUR rerun directory is not explicitly EUR-labelled: {result_dir}")

    annot_manifest_path = require(metadata_dir / "annotation_manifest.tsv")
    trait_manifest_path = require(metadata_dir / "trait_manifest.tsv")
    annots = pd.read_csv(annot_manifest_path, sep="\t")
    traits = pd.read_csv(trait_manifest_path, sep="\t")
    preserved_definitions = compare_annotation_definitions(
        annots, args.expected_definitions
    )
    if set(annots["SNPReferencePopulation"]) != {"EUR"}:
        raise ValueError("Custom annotations are not labelled as EUR-reference annotations")
    if set(annots["AnnotationDerivation"]) != {
        "EAS/BBJ-derived eGene cis intervals"
    }:
        raise ValueError("Annotation derivation provenance is missing or incorrect")
    if set(traits["GWASAncestry"]) != {"EUR"} or set(traits["GenomeBuild"]) != {
        "GRCh37"
    }:
        raise ValueError("GWAS ancestry/build provenance is not EUR/GRCh37")

    baseline_prefix = (
        result_dir
        / "reference/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."
    )
    weights_prefix = (
        result_dir
        / "reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
    )
    print_snps_path = result_dir / "reference/hm3_no_MHC.list.txt"
    regression_snps = read_snp_set(print_snps_path)
    if not regression_snps:
        raise ValueError(f"Empty regression SNP list: {print_snps_path}")

    report: list[dict[str, object]] = []
    reference_hm3: set[str] = set()
    weights_all: set[str] = set()
    source_position_matches = 0
    source_position_mismatches = 0
    selected_by_study = {
        study: selected_genes(args.study_dir, study) for study in STUDIES
    }
    first_annot_prefix = Path(str(annots.iloc[0]["AnnotPrefix"]))

    baseline_categories: int | None = None
    fam_count: int | None = None
    for chrom in range(1, 23):
        bfile = prefix_path(bfile_prefix, chrom, ".bim")
        require(prefix_path(bfile_prefix, chrom, ".bed"))
        fam = require(prefix_path(bfile_prefix, chrom, ".fam"))
        require(prefix_path(frq_prefix, chrom, ".frq"))
        current_fam_count = sum(1 for line in fam.open() if line.strip())
        if fam_count is None:
            fam_count = current_fam_count
        elif current_fam_count != fam_count:
            raise ValueError("EUR PLINK sample count differs across chromosomes")

        bim = pd.read_csv(
            require(bfile),
            sep=r"\s+",
            header=None,
            usecols=[0, 1, 3, 4, 5],
            names=["CHR", "SNP", "BP", "A1", "A2"],
            dtype={"SNP": str},
        )
        if (bim["CHR"] != chrom).any() or bim["SNP"].duplicated().any():
            raise ValueError(f"Invalid chromosome or duplicate SNP IDs in {bfile}")
        bim_positions = dict(zip(bim["SNP"], bim["BP"], strict=True))
        in_hm3 = bim["SNP"].isin(regression_snps)
        reference_hm3.update(bim.loc[in_hm3, "SNP"])
        if chrom == 6:
            mhc = bim.loc[in_hm3 & bim["BP"].between(25_000_000, 34_000_000), "SNP"]
            if not mhc.empty:
                raise ValueError("The HapMap3 non-MHC list contains MHC variants")

        baseline_ld = prefix_path(baseline_prefix, chrom, ".l2.ldscore.gz")
        baseline_annot = prefix_path(baseline_prefix, chrom, ".annot.gz")
        require(prefix_path(baseline_prefix, chrom, ".l2.M"))
        require(prefix_path(baseline_prefix, chrom, ".l2.M_5_50"))
        require(baseline_annot)
        with gzip.open(baseline_annot, "rt") as handle:
            next(handle)
            baseline_annot_rows = sum(1 for _ in handle)
        if baseline_annot_rows != len(bim):
            raise ValueError(
                f"Release mismatch: baseline annotation has {baseline_annot_rows} rows "
                f"but EUR BIM has {len(bim)} rows on chr{chrom}"
            )
        header = pd.read_csv(require(baseline_ld), sep=r"\s+", nrows=0)
        n_categories = len(header.columns) - 3
        if baseline_categories is None:
            baseline_categories = n_categories
        if n_categories != 97 or baseline_categories != n_categories:
            raise ValueError(
                f"Expected EUR baseline-LD v2.2 with 97 annotations; "
                f"found {n_categories} on chromosome {chrom}"
            )
        baseline = pd.read_csv(
            baseline_ld,
            sep=r"\s+",
            usecols=["CHR", "SNP", "BP"],
            dtype={"SNP": str},
        )
        expected_bp = baseline["SNP"].map(bim_positions)
        baseline_shared = expected_bp.notna()
        if (expected_bp[baseline_shared] != baseline.loc[baseline_shared, "BP"]).any():
            raise ValueError(f"EUR baseline-LD/BIM identifier-position mismatch on chr{chrom}")

        weight_path = prefix_path(weights_prefix, chrom, ".l2.ldscore.gz")
        weights = pd.read_csv(
            require(weight_path),
            sep=r"\s+",
            usecols=["CHR", "SNP", "BP"],
            dtype={"SNP": str},
        )
        expected_bp = weights["SNP"].map(bim_positions)
        weights_shared = expected_bp.notna()
        if (expected_bp[weights_shared] != weights.loc[weights_shared, "BP"]).any():
            raise ValueError(f"EUR weights/BIM identifier-position mismatch on chr{chrom}")
        weights_all.update(weights["SNP"])

        frq = pd.read_csv(
            prefix_path(frq_prefix, chrom, ".frq"),
            sep=r"\s+",
            usecols=["CHR", "SNP", "A1", "A2", "MAF"],
            dtype={"SNP": str},
        )
        if len(frq) != len(bim) or set(frq["SNP"]) != set(bim["SNP"]):
            raise ValueError(f"EUR frequency/BIM SNP mismatch on chr{chrom}")

        # Every thin annotation is generated in the same EUR BIM order. Check
        # all headers and one full row count per chromosome; LDSC checks every
        # annotation's row count when custom LD scores are computed.
        for prefix in annots["AnnotPrefix"].drop_duplicates():
            annot_path = Path(f"{prefix}{chrom}.annot.gz")
            with gzip.open(require(annot_path), "rt") as handle:
                if not handle.readline().strip():
                    raise ValueError(f"Empty annotation header: {annot_path}")
        first_annot = Path(f"{first_annot_prefix}{chrom}.annot.gz")
        with gzip.open(first_annot, "rt") as handle:
            next(handle)
            annot_rows = sum(1 for _ in handle)
        if annot_rows != len(bim):
            raise ValueError(f"Annotation/BIM row mismatch on chr{chrom}")

        for study in STUDIES:
            info_path = args.study_dir / study / "INFO" / f"chr{chrom}.csv"
            info = pd.read_csv(
                require(info_path),
                usecols=["GENE", "RSID", "POS"],
                dtype={"GENE": str, "RSID": str},
            )
            info["GENE"] = info["GENE"].str.split(".").str[0]
            info = info[info["GENE"].isin(selected_by_study[study])]
            expected = info["RSID"].map(bim_positions)
            shared = expected.notna()
            source_position_matches += int((expected[shared] == info.loc[shared, "POS"]).sum())
            source_position_mismatches += int((expected[shared] != info.loc[shared, "POS"]).sum())

        if args.check_custom_ldscores:
            for prefix in annots["AnnotPrefix"].drop_duplicates():
                path = Path(f"{prefix}{chrom}.l2.ldscore.gz")
                custom_header = pd.read_csv(require(path), sep=r"\s+", nrows=0)
                if list(custom_header.columns[:3]) != ["CHR", "SNP", "BP"]:
                    raise ValueError(f"Malformed custom LD score: {path}")

        report.append(
            {
                "chromosome": chrom,
                "eur_reference_snps": len(bim),
                "eur_hm3_non_mhc_snps": int(in_hm3.sum()),
                "baseline_ld_snps": len(baseline),
                "baseline_ld_snps_in_eur_bim": int(baseline_shared.sum()),
                "baseline_annotation_rows": baseline_annot_rows,
                "weight_snps": len(weights),
                "weight_snps_in_eur_bim": int(weights_shared.sum()),
                "eur_reference_individuals": current_fam_count,
            }
        )

    if source_position_mismatches:
        raise ValueError(
            f"EAS/BBJ interval-source SNPs disagree with EUR GRCh37 positions: "
            f"{source_position_mismatches} mismatches"
        )
    if source_position_matches < 1_000_000:
        raise ValueError("Too few shared interval-source/EUR SNPs to establish build compatibility")
    if reference_hm3 != regression_snps:
        missing_reference = len(regression_snps - reference_hm3)
        raise ValueError(
            "HapMap3 list and release-matched EUR PLINK panel are incompatible: "
            f"{missing_reference} absent from the EUR panel"
        )
    weight_overlap = len(regression_snps & weights_all)
    if weight_overlap < 1_000_000:
        raise ValueError("Insufficient compatibility with EUR regression weights")

    gwas_rows = []
    for trait in traits.itertuples(index=False):
        sumstats_path = require(Path(trait.SumstatsPath))
        sumstats = pd.read_csv(sumstats_path, sep=r"\s+", usecols=["SNP"], dtype=str)
        if sumstats["SNP"].duplicated().any():
            raise ValueError(f"Duplicate SNP IDs in {sumstats_path}")
        overlap = int(sumstats["SNP"].isin(regression_snps).sum())
        if overlap < 100_000:
            raise ValueError(f"Insufficient GWAS/EUR regression SNP overlap for {trait.Trait}")
        gwas_rows.append(
            {
                "trait": trait.Trait,
                "gwas_snps": len(sumstats),
                "eur_hm3_non_mhc_overlap": overlap,
                "overlap_fraction": overlap / len(sumstats),
            }
        )

    report_frame = pd.DataFrame(report)
    report_frame.to_csv(metadata_dir / "reference_compatibility_by_chr.tsv", sep="\t", index=False)
    pd.DataFrame(gwas_rows).to_csv(
        metadata_dir / "gwas_reference_compatibility.tsv", sep="\t", index=False
    )
    summary = {
        "status": "PASS",
        "annotation_derivation": "EAS/BBJ-derived eGene cis intervals (unchanged)",
        "preserved_annotation_definitions": preserved_definitions,
        "expected_definition_checksums": str(args.expected_definitions.resolve()),
        "sldsc_reference_population": "EUR",
        "genome_build": "GRCh37/hg19",
        "baseline_model": "1000 Genomes Phase 3 EUR baseline-LD v2.2",
        "baseline_annotation_count": baseline_categories,
        "eur_reference_individuals": fam_count,
        "interval_source_position_matches": source_position_matches,
        "interval_source_position_mismatches": source_position_mismatches,
        "eur_hm3_non_mhc_regression_snps": len(regression_snps),
        "eur_weight_overlap_snps": weight_overlap,
        "checked_custom_ldscores": args.check_custom_ldscores,
        "old_eas_artifacts_reused": False,
        "result_directory": str(result_dir),
    }
    (metadata_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
