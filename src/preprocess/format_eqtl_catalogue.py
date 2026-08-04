"""Format eQTL Catalogue studies as chromosome-level traceCB inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


STUDY_SAMPLE_SIZES = {
    "QTD000021": 191,
    "QTD000031": 167,
    "QTD000066": 277,
    "QTD000067": 290,
    "QTD000069": 286,
    "QTD000073": 262,
    "QTD000081": 420,
    "QTD000115": 247,
    "QTD000371": 280,
    "QTD000372": 269,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--studies",
        nargs="+",
        choices=tuple(STUDY_SAMPLE_SIZES),
        default=list(STUDY_SAMPLE_SIZES),
    )
    return parser.parse_args()


def format_study(path: Path, output_dir: Path, sample_size: int) -> None:
    frame = pd.read_csv(path, sep="\t", compression="gzip", dtype={"chromosome": str})
    frame = frame.loc[(frame["chromosome"] != "X") & (frame["type"] == "SNP")].copy()
    frame["chromosome"] = frame["chromosome"].astype(int)
    frame["n"] = sample_size

    for chromosome in range(1, 23):
        chromosome_frame = (
            frame.loc[frame["chromosome"] == chromosome]
            .sort_values("molecular_trait_id")
            .drop_duplicates(subset=["rsid", "gene_id"])
            .copy()
        )
        chromosome_frame["z"] = chromosome_frame["beta"] / chromosome_frame["se"]
        result = chromosome_frame[
            [
                "chromosome",
                "rsid",
                "gene_id",
                "position",
                "alt",
                "ref",
                "beta",
                "se",
                "pvalue",
                "z",
                "n",
            ]
        ].copy()
        result.columns = [
            "CHR",
            "RSID",
            "GENE",
            "POS",
            "A1",
            "A2",
            "BETA",
            "SE",
            "PVAL",
            "Z",
            "N",
        ]
        result.to_csv(output_dir / f"chr{chromosome}.csv", index=False)


def main() -> None:
    args = parse_args()
    for study in args.studies:
        source = args.input_dir / f"{study}.all.tsv.gz"
        if not source.exists():
            raise FileNotFoundError(source)
        study_output = args.output_dir / study
        study_output.mkdir(parents=True, exist_ok=True)
        format_study(source, study_output, STUDY_SAMPLE_SIZES[study])
        print(f"Prepared {study}")


if __name__ == "__main__":
    main()
