"""Prepare chr22 SNP-to-gene annotation for whole-chromosome simulation.

The UKB chr22 test data are in PLINK coordinates that are typically GRCh37/hg19,
so the default source is GENCODE v19 (GRCh37.p13). The script writes:

1. chr22_gencode_v19_genes.tsv: one row per chr22 gene interval.
2. <bim stem>.chr22.snp_gene.tsv: one row per BIM SNP with nearest/overlap gene.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import urllib.request
from pathlib import Path


DEFAULT_GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_19/gencode.v19.annotation.gtf.gz"
)
DEFAULT_SIM_DIR = Path("data/simulation")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create chr22 SNP-to-gene annotation from a PLINK BIM file."
    )
    parser.add_argument(
        "--bim",
        default="/Users/lucajiang/learn/CityU/UKBheight/height_ukb_50k_chr22.bim",
        help="Input PLINK .bim file.",
    )
    parser.add_argument(
        "--gtf",
        default=str(DEFAULT_SIM_DIR / "gencode.v19.annotation.gtf.gz"),
        help="GENCODE GTF path. Downloaded when missing and --download is set.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download the default GENCODE GTF when --gtf does not exist.",
    )
    parser.add_argument(
        "--gencode_url",
        default=DEFAULT_GENCODE_URL,
        help="GENCODE GTF download URL.",
    )
    parser.add_argument(
        "--chromosome",
        default="22",
        help="Chromosome to annotate. Accepts values like 22 or chr22.",
    )
    parser.add_argument(
        "--gene_type",
        default=None,
        help="Optional gene_type filter, for example protein_coding.",
    )
    parser.add_argument(
        "--out_dir",
        default=str(DEFAULT_SIM_DIR),
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default=None,
        help="Output prefix. Defaults to the BIM stem.",
    )
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom


def maybe_download_gtf(gtf_path: Path, url: str, download: bool) -> None:
    if gtf_path.exists():
        return
    if not download:
        raise FileNotFoundError(
            f"{gtf_path} does not exist. Re-run with --download or provide --gtf."
        )
    gtf_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, gtf_path)


GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def parse_gtf_attributes(attr_text: str) -> dict[str, str]:
    return {key: value for key, value in GTF_ATTR_RE.findall(attr_text)}


def strip_version(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def load_gene_intervals(
    gtf_path: Path, chromosome: str, gene_type: str | None
) -> list[dict[str, str | int]]:
    target_chrom = normalize_chrom(chromosome)
    genes: list[dict[str, str | int]] = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = normalize_chrom(fields[0])
            feature = fields[2]
            if chrom != target_chrom or feature != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            row_gene_type = attrs.get("gene_type", attrs.get("gene_biotype", "NA"))
            if gene_type is not None and row_gene_type != gene_type:
                continue
            genes.append(
                {
                    "chr": target_chrom,
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gene_id": strip_version(attrs.get("gene_id", "NA")),
                    "gene_id_version": attrs.get("gene_id", "NA"),
                    "gene_name": attrs.get("gene_name", "NA"),
                    "gene_type": row_gene_type,
                }
            )
    genes.sort(key=lambda row: (int(row["start"]), int(row["end"])))
    return genes


def read_bim(bim_path: Path, chromosome: str) -> list[dict[str, str | int | float]]:
    target_chrom = normalize_chrom(chromosome)
    snps: list[dict[str, str | int | float]] = []
    with bim_path.open() as handle:
        for idx, line in enumerate(handle):
            fields = line.split()
            if len(fields) < 6:
                continue
            chrom = normalize_chrom(fields[0])
            if chrom != target_chrom:
                continue
            snps.append(
                {
                    "snp_index": idx,
                    "chr": chrom,
                    "snp": fields[1],
                    "cm": float(fields[2]),
                    "bp": int(fields[3]),
                    "a1": fields[4],
                    "a2": fields[5],
                }
            )
    return snps


def distance_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    if pos < start:
        return start - pos
    return pos - end


def annotate_snps(
    snps: list[dict[str, str | int | float]],
    genes: list[dict[str, str | int]],
) -> list[dict[str, str | int | float]]:
    rows: list[dict[str, str | int | float]] = []

    for snp in snps:
        pos = int(snp["bp"])
        overlap_ids = [
            idx
            for idx, gene in enumerate(genes)
            if int(gene["start"]) <= pos <= int(gene["end"])
        ]

        if overlap_ids:
            primary_idx = min(
                overlap_ids,
                key=lambda j: (
                    abs(pos - (int(genes[j]["start"]) + int(genes[j]["end"])) // 2),
                    int(genes[j]["start"]),
                ),
            )
        else:
            primary_idx = min(
                range(len(genes)),
                key=lambda j: distance_to_interval(
                    pos, int(genes[j]["start"]), int(genes[j]["end"])
                ),
            )

        primary = genes[primary_idx]
        distance = distance_to_interval(
            pos, int(primary["start"]), int(primary["end"])
        )
        overlap_gene_ids = ";".join(
            str(genes[j]["gene_id"]) for j in sorted(overlap_ids)
        )
        rows.append(
            {
                **snp,
                "gene_id": primary["gene_id"],
                "gene_name": primary["gene_name"],
                "gene_type": primary["gene_type"],
                "gene_start": primary["start"],
                "gene_end": primary["end"],
                "gene_strand": primary["strand"],
                "distance_to_gene": distance,
                "is_gene_body": int(distance == 0),
                "overlap_gene_ids": overlap_gene_ids,
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    bim_path = Path(args.bim)
    gtf_path = Path(args.gtf)
    out_dir = Path(args.out_dir)
    out_prefix = args.out_prefix or bim_path.stem

    maybe_download_gtf(gtf_path, args.gencode_url, args.download)
    genes = load_gene_intervals(gtf_path, args.chromosome, args.gene_type)
    if not genes:
        raise ValueError(f"No genes found for chromosome {args.chromosome}.")
    snps = read_bim(bim_path, args.chromosome)
    if not snps:
        raise ValueError(f"No chr{args.chromosome} SNPs found in {bim_path}.")

    gene_path = out_dir / "chr22_gencode_v19_genes.tsv"
    write_tsv(
        gene_path,
        genes,
        [
            "chr",
            "start",
            "end",
            "strand",
            "gene_id",
            "gene_id_version",
            "gene_name",
            "gene_type",
        ],
    )

    annotated = annotate_snps(snps, genes)
    snp_gene_path = out_dir / f"{out_prefix}.chr22.snp_gene.tsv"
    write_tsv(
        snp_gene_path,
        annotated,
        [
            "snp_index",
            "chr",
            "snp",
            "cm",
            "bp",
            "a1",
            "a2",
            "gene_id",
            "gene_name",
            "gene_type",
            "gene_start",
            "gene_end",
            "gene_strand",
            "distance_to_gene",
            "is_gene_body",
            "overlap_gene_ids",
        ],
    )

    n_gene_body = sum(int(row["is_gene_body"]) for row in annotated)
    print(f"Wrote {len(genes)} genes to {gene_path}")
    print(f"Wrote {len(annotated)} SNP annotations to {snp_gene_path}")
    print(f"{n_gene_body}/{len(annotated)} SNPs overlap at least one gene body")


if __name__ == "__main__":
    main()
