#!/usr/bin/env python3
"""Validate newly computed custom EUR LD-score files before S-LDSC regression."""

from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULT_DIR = Path(
    os.environ.get("TRACECB_ENRICHMENT_DIR", REPO_ROOT / "results/enrichment")
)


def validate_custom_file(path: Path, chrom: int, expected_snps: list[str]) -> int:
    """Validate one custom LD-score file against the exact regression-SNP order."""
    if not path.is_file():
        raise FileNotFoundError(path)

    with gzip.open(path, "rt") as handle:
        header = handle.readline().split()
        if header[:3] != ["CHR", "SNP", "BP"] or len(header) != 4:
            raise ValueError(f"Malformed single-annotation LD-score header: {path}")

        for row_index, expected_snp in enumerate(expected_snps, start=1):
            line = handle.readline()
            if not line:
                raise ValueError(
                    f"Truncated custom LD score at row {row_index}: {path}"
                )
            fields = line.split()
            if len(fields) != 4:
                raise ValueError(
                    f"Malformed custom LD-score row {row_index}: {path}"
                )
            if fields[0] != str(chrom) or fields[1] != expected_snp:
                raise ValueError(
                    f"Custom LD-score SNP/order mismatch at row {row_index}: {path}"
                )

        if handle.readline():
            raise ValueError(f"Extra custom LD-score rows: {path}")

    for suffix in (".l2.M", ".l2.M_5_50"):
        companion = Path(str(path).removesuffix(".l2.ldscore.gz") + suffix)
        if not companion.is_file() or companion.stat().st_size == 0:
            raise FileNotFoundError(companion)
    return len(expected_snps)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=DEFAULT_RESULT_DIR,
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(
        args.result_dir / "metadata/annotation_manifest.tsv", sep="\t"
    )
    prefixes = [Path(value) for value in manifest["AnnotPrefix"].drop_duplicates()]
    regression_snps = set(
        pd.read_csv(
            args.result_dir / "reference/hm3_no_MHC.list.txt", header=None, dtype=str
        )[0]
    )
    total_files = 0
    total_checked_snps = 0
    all_custom_snps_checked = 0
    for chrom in range(1, 23):
        baseline = pd.read_csv(
            args.result_dir
            / (
                "reference/1000G_Phase3_baselineLD_v2.2_exact_hm3/"
                f"baselineLD.{chrom}.l2.ldscore.gz"
            ),
            sep=r"\s+",
            usecols=["SNP"],
            dtype=str,
        )
        expected = baseline.loc[baseline["SNP"].isin(regression_snps), "SNP"].tolist()
        total_checked_snps += len(expected)
        for prefix in prefixes:
            path = Path(f"{prefix}{chrom}.l2.ldscore.gz")
            all_custom_snps_checked += validate_custom_file(path, chrom, expected)
            total_files += 1

    if total_checked_snps != len(regression_snps):
        raise ValueError(
            f"Filtered baseline contains {total_checked_snps} rows, but the exact "
            f"regression list contains {len(regression_snps)} SNPs"
        )

    summary = {
        "status": "PASS",
        "reference_population": "EUR",
        "new_custom_ldscore_files": total_files,
        "expected_custom_ldscore_files": len(prefixes) * 22,
        "exact_stack_regression_snps": len(regression_snps),
        "baseline_snps_checked": total_checked_snps,
        "all_custom_snps_checked": all_custom_snps_checked,
        "cached_eas_ldscores_reused": False,
    }
    out = args.result_dir / "metadata/custom_ldscore_validation.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
