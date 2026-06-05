"""Prepare aligned chr22 eQTL simulation inputs.

This prepares the shared SNP universe across pop1/pop2 simulation genotypes and
1000G EAS/EUR LD panels, then creates gene-window annotations where each gene is
extended by 500 kb on both sides.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_SIM_DIR = Path("data/simulation")
DEFAULT_GTF = DEFAULT_SIM_DIR / "gencode.v19.annotation.gtf.gz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare chr22 eQTL simulation inputs.")
    parser.add_argument(
        "--pop1_prefix",
        default="data/simulation/wg_ukb/height_merge_qc2_chr22",
        help="Population 1 PLINK prefix.",
    )
    parser.add_argument(
        "--pop2_prefix",
        default="data/simulation/ukb/height_ukb_50k_chr22",
        help="Population 2 PLINK prefix.",
    )
    parser.add_argument(
        "--ref1_prefix",
        default="data/simulation/1000G/1000G.EAS.QC.maf.22",
        help="1000G population 1 PLINK prefix for LD.",
    )
    parser.add_argument(
        "--ref2_prefix",
        default="data/simulation/1000G/1000G.EUR.QC.maf.22",
        help="1000G population 2 PLINK prefix for LD.",
    )
    parser.add_argument("--gtf", default=str(DEFAULT_GTF))
    parser.add_argument("--chromosome", default="22")
    parser.add_argument("--window", type=int, default=500_000)
    parser.add_argument("--min_snps_per_gene", type=int, default=5)
    parser.add_argument("--gene_type", default=None)
    parser.add_argument("--out_dir", default=str(DEFAULT_SIM_DIR / "chr22_eqtl"))
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom


def read_bim(prefix: str | Path, chromosome: str) -> list[dict]:
    target_chrom = normalize_chrom(chromosome)
    rows = []
    with Path(f"{prefix}.bim").open() as handle:
        for idx, line in enumerate(handle):
            fields = line.split()
            if len(fields) < 6:
                continue
            chrom = normalize_chrom(fields[0])
            if chrom != target_chrom:
                continue
            rows.append(
                {
                    "idx": idx,
                    "chr": chrom,
                    "snp": fields[1],
                    "cm": fields[2],
                    "bp": int(fields[3]),
                    "a1": fields[4],
                    "a2": fields[5],
                }
            )
    return rows


def unique_by_snp(rows: list[dict]) -> dict[str, dict]:
    counts = Counter(row["snp"] for row in rows)
    return {row["snp"]: row for row in rows if counts[row["snp"]] == 1}


GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def parse_gtf_attributes(attr_text: str) -> dict[str, str]:
    return {key: value for key, value in GTF_ATTR_RE.findall(attr_text)}


def strip_version(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def load_genes(
    gtf_path: str | Path, chromosome: str, window: int, gene_type: str | None
) -> list[dict]:
    target_chrom = normalize_chrom(chromosome)
    genes = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if normalize_chrom(fields[0]) != target_chrom or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            row_gene_type = attrs.get("gene_type", attrs.get("gene_biotype", "NA"))
            if gene_type is not None and row_gene_type != gene_type:
                continue
            start = int(fields[3])
            end = int(fields[4])
            genes.append(
                {
                    "chr": target_chrom,
                    "gene_start": start,
                    "gene_end": end,
                    "window_start": max(1, start - window),
                    "window_end": end + window,
                    "strand": fields[6],
                    "gene_id": strip_version(attrs.get("gene_id", "NA")),
                    "gene_id_version": attrs.get("gene_id", "NA"),
                    "gene_name": attrs.get("gene_name", "NA"),
                    "gene_type": row_gene_type,
                }
            )
    genes.sort(key=lambda row: (row["window_start"], row["window_end"], row["gene_id"]))
    return genes


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_list(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for value in values:
            handle.write(f"{value}\n")


def allele_action(target: dict, query: dict) -> str:
    target_pair = (target["a1"], target["a2"])
    query_pair = (query["a1"], query["a2"])
    if query_pair == target_pair:
        return "same"
    if query_pair == (target_pair[1], target_pair[0]):
        return "flip"
    if set(query_pair) == set(target_pair):
        return "flip"
    return "mismatch"


def write_annot_gz(path: Path, snps: list[dict], genes: list[dict], memberships: dict[str, set[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        header = ["CHR", "BP", "SNP", "CM", "base"] + [gene["gene_id"] for gene in genes]
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for snp in snps:
            snp_genes = memberships.get(snp["snp"], set())
            writer.writerow(
                [
                    snp["chr"],
                    snp["bp"],
                    snp["snp"],
                    snp["cm"],
                    1,
                    *[1 if gene["gene_id"] in snp_genes else 0 for gene in genes],
                ]
            )


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pop1 = unique_by_snp(read_bim(args.pop1_prefix, args.chromosome))
    pop2 = unique_by_snp(read_bim(args.pop2_prefix, args.chromosome))
    ref1 = unique_by_snp(read_bim(args.ref1_prefix, args.chromosome))
    ref2 = unique_by_snp(read_bim(args.ref2_prefix, args.chromosome))
    common_candidate_ids = sorted(
        set(pop1) & set(pop2) & set(ref1) & set(ref2), key=lambda s: pop1[s]["bp"]
    )
    common_ids = []
    pop2_flip = []
    ref1_flip = []
    ref2_flip = []
    mismatches = []
    for snp_id in common_candidate_ids:
        pop2_action = allele_action(pop1[snp_id], pop2[snp_id])
        ref1_action = allele_action(pop1[snp_id], ref1[snp_id])
        ref2_action = allele_action(pop1[snp_id], ref2[snp_id])
        if "mismatch" in (pop2_action, ref1_action, ref2_action):
            mismatches.append(snp_id)
            continue
        common_ids.append(snp_id)
        if pop2_action == "flip":
            pop2_flip.append(snp_id)
        if ref1_action == "flip":
            ref1_flip.append(snp_id)
        if ref2_action == "flip":
            ref2_flip.append(snp_id)
    if not common_ids:
        raise ValueError("No common non-duplicated SNP IDs across the four BIM files.")

    common_rows = []
    for snp_id in common_ids:
        row = pop1[snp_id]
        common_rows.append(
            {
                "snp": snp_id,
                "chr": row["chr"],
                "bp": row["bp"],
                "cm": row["cm"],
                "a1": row["a1"],
                "a2": row["a2"],
                "pop1_idx": pop1[snp_id]["idx"],
                "pop2_idx": pop2[snp_id]["idx"],
                "ref1_idx": ref1[snp_id]["idx"],
                "ref2_idx": ref2[snp_id]["idx"],
            }
        )
    genes = load_genes(args.gtf, args.chromosome, args.window, args.gene_type)
    gene_rows = []
    long_rows = []
    memberships: dict[str, set[str]] = defaultdict(set)
    for gene in genes:
        gene_snps = [
            snp
            for snp in common_rows
            if gene["window_start"] <= int(snp["bp"]) <= gene["window_end"]
        ]
        if len(gene_snps) < args.min_snps_per_gene:
            continue
        gene_rows.append({**gene, "nsnp": len(gene_snps)})
        for snp in gene_snps:
            memberships[snp["snp"]].add(gene["gene_id"])
            long_rows.append(
                {
                    **snp,
                    "gene_id": gene["gene_id"],
                    "gene_name": gene["gene_name"],
                    "gene_type": gene["gene_type"],
                    "gene_start": gene["gene_start"],
                    "gene_end": gene["gene_end"],
                    "window_start": gene["window_start"],
                    "window_end": gene["window_end"],
                    "gene_strand": gene["strand"],
                }
            )

    retained_snps = [row for row in common_rows if row["snp"] in memberships]
    retained_snp_ids = [row["snp"] for row in retained_snps]
    retained_gene_ids = {row["gene_id"] for row in gene_rows}
    gene_rows = [row for row in gene_rows if row["gene_id"] in retained_gene_ids]

    write_list(out_dir / "chr22.common_snps.txt", common_ids)
    write_list(out_dir / "chr22.pop2_flip_to_pop1.txt", pop2_flip)
    write_list(out_dir / "chr22.ref1_flip_to_pop1.txt", ref1_flip)
    write_list(out_dir / "chr22.ref2_flip_to_pop1.txt", ref2_flip)
    write_list(out_dir / "chr22.allele_mismatch_excluded.txt", mismatches)
    write_list(out_dir / "chr22.print_snps.txt", retained_snp_ids)
    write_tsv(
        out_dir / "chr22.common_snps.tsv",
        common_rows,
        ["snp", "chr", "bp", "cm", "a1", "a2", "pop1_idx", "pop2_idx", "ref1_idx", "ref2_idx"],
    )
    write_tsv(
        out_dir / "chr22.genes_500kb.tsv",
        gene_rows,
        [
            "chr",
            "gene_start",
            "gene_end",
            "window_start",
            "window_end",
            "strand",
            "gene_id",
            "gene_id_version",
            "gene_name",
            "gene_type",
            "nsnp",
        ],
    )
    write_tsv(
        out_dir / "chr22.snp_gene_500kb.tsv",
        long_rows,
        [
            "snp",
            "chr",
            "bp",
            "cm",
            "a1",
            "a2",
            "pop1_idx",
            "pop2_idx",
            "ref1_idx",
            "ref2_idx",
            "gene_id",
            "gene_name",
            "gene_type",
            "gene_start",
            "gene_end",
            "window_start",
            "window_end",
            "gene_strand",
        ],
    )
    write_annot_gz(out_dir / "22.annot.gz", retained_snps, gene_rows, memberships)
    print(f"Common non-duplicated SNP candidates: {len(common_candidate_ids)}")
    print(f"Common allele-harmonizable SNPs: {len(common_ids)}")
    print(f"Excluded allele mismatches: {len(mismatches)}")
    print(f"pop2/ref1/ref2 flips to pop1 allele order: {len(pop2_flip)}/{len(ref1_flip)}/{len(ref2_flip)}")
    print(f"Retained SNPs in at least one 500kb gene window: {len(retained_snps)}")
    print(f"Retained genes with >= {args.min_snps_per_gene} SNPs: {len(gene_rows)}")
    print(f"SNP-gene memberships: {len(long_rows)}")
    print(f"Wrote outputs under {out_dir}")


if __name__ == "__main__":
    main()
