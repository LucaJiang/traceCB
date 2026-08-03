#!/usr/bin/env python3
"""Filter EUR baseline-LD v2.2 rows to the exact custom-LD regression SNPs."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-prefix", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--regression-snps", type=Path, required=True)
    parser.add_argument("--custom-prefix", type=Path, required=True)
    return parser.parse_args()


def snp_column(path: Path) -> list[str]:
    values = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split()
        snp_index = header.index("SNP")
        for line in handle:
            values.append(line.split()[snp_index])
    return values


def main() -> None:
    args = parse_args()
    keep = {
        line.strip().split()[0]
        for line in args.regression_snps.open()
        if line.strip()
    }
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    all_written: set[str] = set()
    for chrom in range(1, 23):
        source_ld = Path(f"{args.source_prefix}{chrom}.l2.ldscore.gz")
        output_ld = Path(f"{args.output_prefix}{chrom}.l2.ldscore.gz")
        written = 0
        with gzip.open(source_ld, "rt") as source, gzip.open(
            output_ld, "wt", compresslevel=6
        ) as output:
            header = source.readline()
            output.write(header)
            columns = header.rstrip("\n").split()
            snp_index = columns.index("SNP")
            if len(columns) - 3 != 97:
                raise ValueError(f"Not a 97-annotation baseline-LD file: {source_ld}")
            for line in source:
                snp = line.split(None, 2)[snp_index]
                if snp in keep:
                    output.write(line)
                    all_written.add(snp)
                    written += 1

        for suffix in (".annot.gz", ".l2.M", ".l2.M_5_50"):
            shutil.copy2(
                Path(f"{args.source_prefix}{chrom}{suffix}"),
                Path(f"{args.output_prefix}{chrom}{suffix}"),
            )
        filtered_snps = snp_column(output_ld)
        custom_snps = snp_column(Path(f"{args.custom_prefix}{chrom}.l2.ldscore.gz"))
        if filtered_snps != custom_snps:
            raise ValueError(
                f"Filtered baseline and custom LD-score SNP rows differ on chr{chrom}"
            )
        rows.append({"chromosome": chrom, "retained_rows": written})

    if all_written != keep:
        raise ValueError(
            f"Filtered baseline retained {len(all_written)} of {len(keep)} requested SNPs"
        )
    summary = {
        "status": "PASS",
        "source_model": "1000 Genomes Phase 3 EUR baseline-LD v2.2",
        "annotation_count": 97,
        "operation": "row filtering only; LD scores were not recomputed or altered",
        "retained_regression_snps": len(all_written),
        "custom_row_order_identical": True,
        "chromosomes": rows,
    }
    report = args.output_prefix.parent / "FILTERING_PROVENANCE.json"
    report.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
