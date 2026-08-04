#!/usr/bin/env python3
"""Step 01: build the exact HapMap3 non-MHC intersection for the EUR stack."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-list", type=Path, required=True)
    parser.add_argument("--bfile-prefix", type=Path, required=True)
    parser.add_argument("--baseline-prefix", type=Path, required=True)
    parser.add_argument("--weights-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = set(pd.read_csv(args.source_list, header=None, dtype=str)[0])
    retained: list[str] = []
    rows = []
    for chrom in range(1, 23):
        bim = pd.read_csv(
            f"{args.bfile_prefix}{chrom}.bim",
            sep=r"\s+",
            header=None,
            usecols=[1, 3],
            names=["SNP", "BP"],
            dtype={"SNP": str},
        )
        baseline = pd.read_csv(
            f"{args.baseline_prefix}{chrom}.l2.ldscore.gz",
            sep=r"\s+",
            usecols=["SNP", "BP"],
            dtype={"SNP": str},
        )
        weights = pd.read_csv(
            f"{args.weights_prefix}{chrom}.l2.ldscore.gz",
            sep=r"\s+",
            usecols=["SNP", "BP"],
            dtype={"SNP": str},
        )
        bim_bp = dict(zip(bim["SNP"], bim["BP"], strict=True))
        for name, frame in (("baseline", baseline), ("weights", weights)):
            expected = frame["SNP"].map(bim_bp)
            shared = expected.notna()
            if (expected[shared] != frame.loc[shared, "BP"]).any():
                raise ValueError(f"{name}/EUR BIM position mismatch on chromosome {chrom}")

        common = source & set(bim["SNP"]) & set(baseline["SNP"]) & set(weights["SNP"])
        if chrom == 6:
            mhc = set(bim.loc[bim["BP"].between(25_000_000, 34_000_000), "SNP"])
            common -= mhc
        ordered = [snp for snp in baseline["SNP"] if snp in common]
        if len(ordered) != len(common):
            raise ValueError(f"Duplicate or missing intersected SNPs on chromosome {chrom}")
        retained.extend(ordered)
        rows.append(
            {
                "chromosome": chrom,
                "source_hm3_non_mhc": int(bim["SNP"].isin(source).sum()),
                "baseline_ld": len(baseline),
                "weights": len(weights),
                "exact_stack_intersection": len(ordered),
            }
        )

    if len(retained) < 1_000_000 or len(retained) != len(set(retained)):
        raise ValueError(f"Unexpected compatible regression SNP count: {len(retained)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(retained) + "\n")
    pd.DataFrame(rows).to_csv(
        args.output.with_name("hm3_compatibility_by_chr.tsv"), sep="\t", index=False
    )
    print(f"[done] wrote {len(retained)} exact-stack HapMap3 non-MHC SNPs", flush=True)


if __name__ == "__main__":
    main()
