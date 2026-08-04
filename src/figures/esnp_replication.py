#!/usr/bin/env python3
"""Measure eSNP replication for Original, traceC, and traceCB.

The analysis unit is a significant SNP-gene pair (p < 1e-5 by default), which
matches the way TAR_*eSNP is counted in the GMM summary files.  Replication is
evaluated against:

* hum0197: the matching broad cell type (Mono, CD4T, CD8T, B, or NK)
* hum0343: the same bulk whole-blood eQTL data for every study
* CIMA: matching broad cell type, using significant lead cis-eQTL records

Alleles are harmonized before beta signs are compared.  Same, swapped,
strand-complemented, and strand-complemented/swapped allele pairs are handled;
unharmonizable pairs are retained in the summary as allele-mismatch QC counts.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import polars as pl
import pyarrow.parquet as pq


DEFAULT_THRESHOLD = 1e-5

METHODS = ("Original", "traceC", "traceCB")
REPLICATES = ("hum0197", "hum0343", "CIMA")
METHOD_COLUMNS = {
    "Original": ("TAR_SBETA", "TAR_SPVAL"),
    "traceC": ("TAR_CBETA", "TAR_CPVAL"),
    "traceCB": ("TAR_TBETA", "TAR_TPVAL"),
}
SUMMARY_COLUMNS = {
    "Original": "TAR_SeSNP",
    "traceC": "TAR_CeSNP",
    "traceCB": "TAR_TeSNP",
}
CELLTYPE_TO_HUM0197 = {
    "Monocytes": "Mono",
    "CD4+T_cells": "CD4T",
    "CD8+T_cells": "CD8T",
    "B_cells": "B",
    "NK_cells": "NK",
}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

DETAIL_COLUMNS = (
    "QTDid",
    "study_name",
    "cell_type",
    "replicate",
    "replicate_cell_type",
    "replicate_data_level",
    "replicate_variant_scope",
    "method",
    "GENE",
    "RSID",
    "CHR",
    "POS",
    "my_A1",
    "my_A2",
    "my_beta",
    "my_pval",
    "replicate_variant_id",
    "replicate_A1",
    "replicate_A2",
    "replicate_beta",
    "replicate_pval",
    "allele_alignment",
    "aligned_replicate_beta",
    "same_sign",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--hum0197-dir", type=Path, required=True)
    parser.add_argument("--hum0343-esnp", type=Path, required=True)
    parser.add_argument("--cima-lead-eqtl", type=Path, required=True)
    parser.add_argument("--gene-annotation", type=Path, required=True)
    parser.add_argument("--gtex-lookup", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="Number of studies processed concurrently (default: 2).",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Rebuild filtered replicate parquet caches.",
    )
    parser.add_argument(
        "--study-ids",
        nargs="+",
        help="Optional subset of QTDid values (default: all in metadata.json).",
    )
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        type=int,
        default=list(range(1, 23)),
        help="Optional chromosome subset, mainly for validation runs.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {message}", flush=True)


def load_metadata() -> dict:
    metadata_path = Path(__file__).with_name("metadata.json")
    with metadata_path.open() as handle:
        return json.load(handle)


def source_signature(source: Path, threshold: float, kind: str) -> dict:
    stat = source.stat()
    return {
        "kind": kind,
        "source": str(source.resolve()),
        "source_size": stat.st_size,
        "source_mtime_ns": stat.st_mtime_ns,
        "threshold": threshold,
        "cache_version": 1,
    }


def cache_is_current(parquet_path: Path, signature: dict) -> bool:
    metadata_path = parquet_path.with_suffix(".json")
    if not parquet_path.exists() or not metadata_path.exists():
        return False
    try:
        with metadata_path.open() as handle:
            observed = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return False
    return observed == signature


def write_json_atomic(data: dict, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def deduplicate_cache(
    raw_path: Path, final_path: Path, subset: list[str], signature: dict
) -> int:
    deduplicated = final_path.with_suffix(".deduplicated.tmp.parquet")
    (
        pl.scan_parquet(raw_path)
        .unique(subset=subset, keep="first", maintain_order=False)
        .sink_parquet(
            deduplicated,
            compression="zstd",
            row_group_size=100_000,
            statistics=True,
        )
    )
    os.replace(deduplicated, final_path)
    raw_path.unlink()
    count = (
        pl.scan_parquet(final_path)
        .select(pl.len().alias("n"))
        .collect()
        .item()
    )
    write_json_atomic(signature, final_path.with_suffix(".json"))
    return int(count)


def prepare_gtex_lookup_cache(
    source: Path,
    cache_path: Path,
    refresh: bool,
) -> int:
    signature = source_signature(source, -1.0, "gtex_b38_rsid_lookup")
    if not refresh and cache_is_current(cache_path, signature):
        count = pl.scan_parquet(cache_path).select(pl.len()).collect().item()
        log(f"Using GTEx rsID cache {cache_path.name}: {count:,} variants")
        return int(count)

    log(f"Building GRCh38-to-rsID cache from {source.name}")
    temporary = cache_path.with_suffix(".tmp.parquet")
    chromosome = (
        pl.col("chr")
        .str.strip_prefix("chr")
        .cast(pl.Int16, strict=False)
        .alias("replicate_CHR")
    )
    query = (
        pl.scan_csv(source, separator="\t")
        .select(
            pl.col("variant_id")
            .str.strip_suffix("_b38")
            .alias("replicate_variant_id"),
            pl.col("rs_id_dbSNP151_GRCh38p7").alias("replicate_RSID"),
            chromosome,
        )
        .filter(
            pl.col("replicate_CHR").is_between(1, 22)
            & pl.col("replicate_RSID").is_not_null()
            & (pl.col("replicate_RSID") != ".")
        )
    )
    query.sink_parquet(
        temporary,
        compression="zstd",
        row_group_size=100_000,
        statistics=True,
    )
    os.replace(temporary, cache_path)
    count = int(pl.scan_parquet(cache_path).select(pl.len()).collect().item())
    write_json_atomic(signature, cache_path.with_suffix(".json"))
    log(f"Built {cache_path.name}: {count:,} GRCh38 variants with rsID")
    return count


def prepare_hum0197_raw_cache(
    source: Path,
    raw_cache_path: Path,
    threshold: float,
    refresh: bool,
) -> int:
    signature = source_signature(source, threshold, "hum0197_raw")
    if not refresh and cache_is_current(raw_cache_path, signature):
        count = pl.scan_parquet(raw_cache_path).select(pl.len()).collect().item()
        log(f"Using hum0197 cache {raw_cache_path.name}: {count:,} eSNP pairs")
        return int(count)

    log(f"Building hum0197 raw cache from {source.name}")
    temporary_raw = raw_cache_path.with_suffix(".raw.tmp.parquet")
    variant = pl.col("variant_id").str.extract_groups(
        r"^chr(?P<CHR>[0-9]+)_(?P<POS>[0-9]+)_"
        r"(?P<REF>[^_]+)_(?P<ALT>[^_]+)$"
    )
    query = (
        pl.scan_csv(source, separator="\t")
        .filter(pl.col("pval_nominal") < threshold)
        .select(
            pl.col("phenotype_id")
            .str.split(".")
            .list.first()
            .alias("GENE"),
            pl.col("variant_id").alias("replicate_variant_id"),
            pl.col("slope").cast(pl.Float64).alias("replicate_beta"),
            pl.col("pval_nominal").cast(pl.Float64).alias("replicate_pval"),
            variant.alias("variant"),
        )
        .unnest("variant")
        .with_columns(
            pl.col("CHR").cast(pl.Int16),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
        )
        .filter(pl.col("CHR").is_between(1, 22))
        .rename(
            {
                "CHR": "replicate_CHR",
                "POS": "replicate_POS",
                "ALT": "replicate_A1",
                "REF": "replicate_A2",
            }
        )
    )
    query.sink_parquet(
        temporary_raw,
        compression="zstd",
        row_group_size=100_000,
        statistics=True,
    )
    count = deduplicate_cache(
        temporary_raw,
        raw_cache_path,
        ["GENE", "replicate_variant_id"],
        signature,
    )
    log(f"Built {raw_cache_path.name}: {count:,} unique eSNP pairs")
    return count


def prepare_hum0197_mapped_cache(
    raw_cache_path: Path,
    gtex_cache_path: Path,
    mapped_cache_path: Path,
    signature: dict,
    refresh: bool,
) -> int:
    if not refresh and cache_is_current(mapped_cache_path, signature):
        count = pl.scan_parquet(mapped_cache_path).select(pl.len()).collect().item()
        log(
            f"Using hum0197 rsID cache {mapped_cache_path.name}: "
            f"{count:,} mapped eSNP pairs"
        )
        return int(count)

    log(f"Mapping {raw_cache_path.name} to dbSNP151 rsIDs")
    work_dir = Path(
        tempfile.mkdtemp(prefix=".hum0197_map_", dir=mapped_cache_path.parent)
    )
    part_paths: list[Path] = []
    try:
        for chromosome in range(1, 23):
            part_path = work_dir / f"chr{chromosome}.parquet"
            oasis = pl.scan_parquet(raw_cache_path).filter(
                pl.col("replicate_CHR") == chromosome
            )
            lookup = (
                pl.scan_parquet(gtex_cache_path)
                .filter(pl.col("replicate_CHR") == chromosome)
                .drop("replicate_CHR")
            )
            (
                oasis.join(
                    lookup,
                    on="replicate_variant_id",
                    how="inner",
                    validate="m:1",
                )
                .sink_parquet(
                    part_path,
                    compression="zstd",
                    row_group_size=100_000,
                    statistics=True,
                )
            )
            part_paths.append(part_path)

        temporary = mapped_cache_path.with_suffix(".tmp.parquet")
        (
            pl.scan_parquet([str(path) for path in part_paths])
            .sort("replicate_CHR", "replicate_POS", "GENE")
            .sink_parquet(
                temporary,
                compression="zstd",
                row_group_size=100_000,
                statistics=True,
            )
        )
        os.replace(temporary, mapped_cache_path)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    count = int(pl.scan_parquet(mapped_cache_path).select(pl.len()).collect().item())
    write_json_atomic(signature, mapped_cache_path.with_suffix(".json"))
    log(f"Built {mapped_cache_path.name}: {count:,} rsID-mapped eSNP pairs")
    return count


def prepare_hum0343_cache(
    source: Path,
    cache_path: Path,
    threshold: float,
    refresh: bool,
) -> int:
    signature = source_signature(source, threshold, "hum0343")
    if not refresh and cache_is_current(cache_path, signature):
        count = pl.scan_parquet(cache_path).select(pl.len()).collect().item()
        log(f"Using hum0343 cache {cache_path.name}: {count:,} eSNP pairs")
        return int(count)

    log(f"Building hum0343 cache from {source.name}")
    raw_path = cache_path.with_suffix(".raw.tmp.parquet")
    query = (
        pl.scan_csv(source)
        .filter(pl.col("pval") < threshold)
        .select(
            pl.col("gene").str.split(".").list.first().alias("GENE"),
            pl.col("rsid").alias("replicate_RSID"),
            pl.col("variant_id").alias("replicate_variant_id"),
            pl.col("chr").cast(pl.Int16).alias("replicate_CHR"),
            pl.col("pos").cast(pl.Int64).alias("replicate_POS"),
            pl.col("alt").str.to_uppercase().alias("replicate_A1"),
            pl.col("ref").str.to_uppercase().alias("replicate_A2"),
            pl.col("beta").cast(pl.Float64).alias("replicate_beta"),
            pl.col("pval").cast(pl.Float64).alias("replicate_pval"),
        )
        .filter(
            pl.col("replicate_CHR").is_between(1, 22)
            & pl.col("replicate_RSID").is_not_null()
        )
    )
    query.sink_parquet(
        raw_path,
        compression="zstd",
        row_group_size=100_000,
        statistics=True,
    )
    count = deduplicate_cache(
        raw_path,
        cache_path,
        ["GENE", "replicate_RSID"],
        signature,
    )
    log(f"Built {cache_path.name}: {count:,} unique eSNP pairs")
    return count


def starts_with_any(column: str, prefixes: tuple[str, ...]) -> pl.Expr:
    return pl.any_horizontal(
        [pl.col(column).str.starts_with(prefix) for prefix in prefixes]
    )


def cima_broad_celltype() -> pl.Expr:
    celltype = pl.col("celltype")
    return (
        pl.when(
            celltype.is_in(["CD4", "CD4T"])
            | starts_with_any("celltype", ("CD4_",))
        )
        .then(pl.lit("CD4+T_cells"))
        .when(
            celltype.is_in(["CD8", "CD8T"])
            | starts_with_any("celltype", ("CD8_",))
        )
        .then(pl.lit("CD8+T_cells"))
        .when(
            (celltype == "B")
            | starts_with_any(
                "celltype",
                (
                    "Bn_",
                    "Transitional_B_",
                    "Switched_Bm_",
                    "Unswitched_Bm_",
                    "pre-Switched_Bm_",
                    "Atypical_Bm_",
                ),
            )
        )
        .then(pl.lit("B_cells"))
        .when(
            celltype.is_in(["Mono", "Monocyte"])
            | starts_with_any("celltype", ("cMono_", "ncMono_", "intMono_"))
        )
        .then(pl.lit("Monocytes"))
        .when(
            (celltype == "NK")
            | starts_with_any(
                "celltype",
                (
                    "Mature_NK_",
                    "Terminal_NK_",
                    "Transitional_NK_",
                    "NK_bright_",
                    "Inflamed_NK_",
                    "Proliferative_NK_",
                ),
            )
        )
        .then(pl.lit("NK_cells"))
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias("replicate_cell_type")
    )


def cache_counts_by_celltype(cache_path: Path) -> dict[str, int]:
    counts = (
        pl.scan_parquet(cache_path)
        .group_by("replicate_cell_type")
        .agg(pl.len().alias("n"))
        .collect()
    )
    observed = {
        row["replicate_cell_type"]: int(row["n"])
        for row in counts.iter_rows(named=True)
    }
    return {celltype: observed.get(celltype, 0) for celltype in CELLTYPE_TO_HUM0197}


def prepare_cima_raw_cache(
    source: Path,
    cache_path: Path,
    threshold: float,
    refresh: bool,
) -> dict[str, int]:
    signature = source_signature(source, threshold, "cima_lead_cis_eqtl_raw")
    if not refresh and cache_is_current(cache_path, signature):
        counts = cache_counts_by_celltype(cache_path)
        log(
            f"Using CIMA raw cache {cache_path.name}: "
            f"{sum(counts.values()):,} lead eSNP pairs"
        )
        return counts

    log(f"Building CIMA lead cis-eQTL cache from {source.name}")
    temporary = cache_path.with_suffix(".tmp.parquet")
    query = (
        pl.scan_csv(source)
        .filter(
            (pl.col("analysis") == "cis-eQTL")
            & (pl.col("pval_nominal") < threshold)
        )
        .with_columns(cima_broad_celltype())
        .filter(pl.col("replicate_cell_type").is_not_null())
        .select(
            "phenotype_id",
            pl.col("celltype").alias("replicate_subtype"),
            "replicate_cell_type",
            (
                pl.col("variant_id")
                + "_"
                + pl.col("A2(REF)").str.to_uppercase()
                + "_"
                + pl.col("A1(ALT/effect allele)").str.to_uppercase()
            ).alias("replicate_variant_id"),
            pl.col("variant_id")
            .str.extract(r"^chr([0-9]+)_", 1)
            .cast(pl.Int16, strict=False)
            .alias("replicate_CHR"),
            pl.col("variant_id")
            .str.extract(r"^chr[0-9]+_([0-9]+)$", 1)
            .cast(pl.Int64, strict=False)
            .alias("replicate_POS"),
            pl.col("A1(ALT/effect allele)")
            .str.to_uppercase()
            .alias("replicate_A1"),
            pl.col("A2(REF)").str.to_uppercase().alias("replicate_A2"),
            pl.col("slope").cast(pl.Float64).alias("replicate_beta"),
            pl.col("pval_nominal").cast(pl.Float64).alias("replicate_pval"),
        )
        .filter(
            pl.col("replicate_CHR").is_between(1, 22)
            & pl.col("replicate_POS").is_not_null()
            & pl.col("replicate_A1").is_not_null()
            & pl.col("replicate_A2").is_not_null()
        )
        .sort("replicate_pval")
        .unique(
            subset=[
                "replicate_cell_type",
                "phenotype_id",
                "replicate_variant_id",
            ],
            keep="first",
            maintain_order=True,
        )
    )
    query.sink_parquet(
        temporary,
        compression="zstd",
        row_group_size=100_000,
        statistics=True,
    )
    os.replace(temporary, cache_path)
    write_json_atomic(signature, cache_path.with_suffix(".json"))
    counts = cache_counts_by_celltype(cache_path)
    log(
        f"Built {cache_path.name}: {sum(counts.values()):,} unique lead eSNP pairs"
    )
    return counts


def load_gene_name_map(
    gene_annotation: Path,
    hum0343_esnp: Path,
) -> pl.DataFrame:
    records: list[tuple[str, str]] = []
    with gene_annotation.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gene_id = re.search(r'gene_id "([^"]+)"', fields[8])
            gene_name = re.search(r'gene_name "([^"]+)"', fields[8])
            if gene_id and gene_name:
                records.append(
                    (gene_name.group(1), gene_id.group(1).split(".")[0])
                )

    gtf = (
        pl.DataFrame(records, schema=["phenotype_id", "GENE"], orient="row")
        .unique()
        .group_by("phenotype_id")
        .agg(
            pl.col("GENE").n_unique().alias("n_gene_ids"),
            pl.col("GENE").first(),
        )
        .filter(pl.col("n_gene_ids") == 1)
        .drop("n_gene_ids")
        .with_columns(pl.lit(0).alias("priority"))
    )
    hum0343 = (
        pl.scan_csv(hum0343_esnp)
        .select(
            pl.col("gene_name").alias("phenotype_id"),
            pl.col("gene").str.split(".").list.first().alias("GENE"),
        )
        .drop_nulls()
        .unique()
        .group_by("phenotype_id")
        .agg(
            pl.col("GENE").n_unique().alias("n_gene_ids"),
            pl.col("GENE").first(),
        )
        .filter(pl.col("n_gene_ids") == 1)
        .drop("n_gene_ids")
        .with_columns(pl.lit(1).alias("priority"))
        .collect(engine="streaming")
    )
    return (
        pl.concat([gtf, hum0343])
        .sort("priority")
        .unique(subset="phenotype_id", keep="first", maintain_order=True)
        .drop("priority")
    )


def prepare_cima_mapped_cache(
    raw_cache_path: Path,
    gtex_cache_path: Path,
    mapped_cache_path: Path,
    gene_annotation: Path,
    hum0343_esnp: Path,
    signature: dict,
    refresh: bool,
) -> dict[str, int]:
    if not refresh and cache_is_current(mapped_cache_path, signature):
        counts = cache_counts_by_celltype(mapped_cache_path)
        log(
            f"Using CIMA mapped cache {mapped_cache_path.name}: "
            f"{sum(counts.values()):,} pairs"
        )
        return counts

    log("Mapping CIMA gene symbols and GRCh38 variants to Ensembl IDs and rsIDs")
    gene_map = load_gene_name_map(gene_annotation, hum0343_esnp)
    mapped_genes = pl.scan_parquet(raw_cache_path).join(
        gene_map.lazy(), on="phenotype_id", how="inner", validate="m:1"
    )
    work_dir = Path(
        tempfile.mkdtemp(prefix=".cima_map_", dir=mapped_cache_path.parent)
    )
    part_paths: list[Path] = []
    try:
        for chromosome in range(1, 23):
            part_path = work_dir / f"chr{chromosome}.parquet"
            cima = mapped_genes.filter(pl.col("replicate_CHR") == chromosome)
            lookup = (
                pl.scan_parquet(gtex_cache_path)
                .filter(pl.col("replicate_CHR") == chromosome)
                .drop("replicate_CHR")
            )
            mapped_part = (
                cima.join(
                    lookup,
                    on="replicate_variant_id",
                    how="inner",
                    validate="m:1",
                )
                .sort("replicate_pval")
                .unique(
                    subset=["replicate_cell_type", "GENE", "replicate_RSID"],
                    keep="first",
                    maintain_order=True,
                )
                .collect(engine="streaming")
            )
            if not mapped_part.is_empty():
                mapped_part.write_parquet(
                    part_path,
                    compression="zstd",
                    statistics=True,
                )
                part_paths.append(part_path)

        temporary = mapped_cache_path.with_suffix(".tmp.parquet")
        if not part_paths:
            raise RuntimeError("No CIMA records could be mapped to Ensembl ID and rsID")
        (
            pl.scan_parquet([str(path) for path in part_paths])
            .sort("replicate_cell_type", "replicate_CHR", "replicate_POS", "GENE")
            .sink_parquet(
                temporary,
                compression="zstd",
                row_group_size=100_000,
                statistics=True,
            )
        )
        os.replace(temporary, mapped_cache_path)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    write_json_atomic(signature, mapped_cache_path.with_suffix(".json"))
    counts = cache_counts_by_celltype(mapped_cache_path)
    log(
        f"Built {mapped_cache_path.name}: {sum(counts.values()):,} "
        "Ensembl/rsID-mapped lead eSNP pairs"
    )
    return counts


def prepare_replicate_caches(
    hum0197_dir: Path,
    hum0343_esnp: Path,
    cima_lead_eqtl: Path,
    gene_annotation: Path,
    gtex_lookup: Path,
    cache_dir: Path,
    threshold: float,
    refresh: bool,
) -> tuple[dict[str, Path], dict[str, int], dict[str, int]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    counts: dict[str, int] = {}
    mapped_counts: dict[str, int] = {}
    if not gtex_lookup.exists():
        raise FileNotFoundError(gtex_lookup)
    gtex_cache_path = cache_dir / "gtex_b38_rsid.parquet"
    prepare_gtex_lookup_cache(gtex_lookup, gtex_cache_path, refresh)
    with gtex_cache_path.with_suffix(".json").open() as handle:
        gtex_signature = json.load(handle)

    for cell_type, alias in CELLTYPE_TO_HUM0197.items():
        source = hum0197_dir / f"{alias}_PC15_MAF0.05.cis_nominal.txt.gz"
        if not source.exists():
            raise FileNotFoundError(source)
        raw_cache_path = cache_dir / f"hum0197_{alias}_all_significant.parquet"
        mapped_cache_path = cache_dir / f"hum0197_{alias}_rsid.parquet"
        raw_count = prepare_hum0197_raw_cache(
            source, raw_cache_path, threshold, refresh
        )
        mapped_signature = {
            "kind": "hum0197_rsid_mapped",
            "raw_source_signature": source_signature(
                source, threshold, "hum0197_raw"
            ),
            "gtex_source_signature": gtex_signature,
            "cache_version": 1,
        }
        mapped_count = prepare_hum0197_mapped_cache(
            raw_cache_path,
            gtex_cache_path,
            mapped_cache_path,
            mapped_signature,
            refresh,
        )
        paths[f"hum0197:{cell_type}"] = mapped_cache_path
        counts[f"hum0197:{cell_type}"] = raw_count
        mapped_counts[f"hum0197:{cell_type}"] = mapped_count

    if not hum0343_esnp.exists():
        raise FileNotFoundError(hum0343_esnp)
    hum0343_cache = cache_dir / "hum0343.parquet"
    paths["hum0343"] = hum0343_cache
    counts["hum0343"] = prepare_hum0343_cache(
        hum0343_esnp, hum0343_cache, threshold, refresh
    )
    mapped_counts["hum0343"] = counts["hum0343"]

    if not cima_lead_eqtl.exists():
        raise FileNotFoundError(cima_lead_eqtl)
    if not gene_annotation.exists():
        raise FileNotFoundError(gene_annotation)
    cima_raw_cache = cache_dir / "CIMA_lead_cis_eqtl_all_significant.parquet"
    cima_mapped_cache = cache_dir / "CIMA_lead_cis_eqtl_rsid.parquet"
    cima_raw_counts = prepare_cima_raw_cache(
        cima_lead_eqtl, cima_raw_cache, threshold, refresh
    )
    cima_signature = {
        "kind": "cima_lead_cis_eqtl_ensembl_rsid_mapped",
        "raw_source_signature": source_signature(
            cima_lead_eqtl, threshold, "cima_lead_cis_eqtl_raw"
        ),
        "gene_annotation_signature": source_signature(
            gene_annotation, -1.0, "gene_annotation"
        ),
        "hum0343_gene_map_signature": source_signature(
            hum0343_esnp, -1.0, "hum0343_gene_map"
        ),
        "gtex_source_signature": gtex_signature,
        "celltype_mapping_version": 1,
        "cache_version": 1,
    }
    cima_mapped_counts = prepare_cima_mapped_cache(
        cima_raw_cache,
        gtex_cache_path,
        cima_mapped_cache,
        gene_annotation,
        hum0343_esnp,
        cima_signature,
        refresh,
    )
    for cell_type in CELLTYPE_TO_HUM0197:
        key = f"CIMA:{cell_type}"
        paths[key] = cima_mapped_cache
        counts[key] = cima_raw_counts[cell_type]
        mapped_counts[key] = cima_mapped_counts[cell_type]
    return paths, counts, mapped_counts


def expected_counts(summary_path: Path) -> dict[str, int]:
    summary = pl.read_csv(summary_path)
    summary = summary.with_columns(
        pl.col("TAR_CeSNP").fill_null(pl.col("TAR_SeSNP")),
    ).with_columns(
        pl.col("TAR_TeSNP").fill_null(pl.col("TAR_CeSNP")),
    )
    result: dict[str, int] = {}
    for method, column in SUMMARY_COLUMNS.items():
        value = summary.select(pl.col(column).fill_null(0).sum()).item()
        result[method] = int(round(value))
    return result


def method_structs(beta_source: str | None = None, pval_source: str | None = None) -> pl.Expr:
    expressions = []
    for method, (beta_column, pval_column) in METHOD_COLUMNS.items():
        beta = pl.col(beta_source or beta_column).cast(pl.Float64).alias("my_beta")
        pval = pl.col(pval_source or pval_column).cast(pl.Float64).alias("my_pval")
        expressions.append(
            pl.struct(pl.lit(method).alias("method"), beta, pval)
        )
    return pl.concat_list(*expressions).alias("association")


def load_my_esnps(
    study_dir: Path,
    cell_type: str,
    chromosome: int,
    threshold: float,
) -> tuple[pl.DataFrame, dict[str, int]]:
    gmm_dir = study_dir / "GMM" / f"chr{chromosome}"
    summary_path = gmm_dir / "summary.csv"
    gmm_files = sorted(gmm_dir.glob("*.parquet"))
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)

    summary = pl.read_csv(summary_path)
    gmm_genes = [path.stem for path in gmm_files]
    fallback_genes = (
        summary.filter(
            pl.col("TAR_SeSNP").is_not_null()
            & ~pl.col("GENE").is_in(gmm_genes)
        )
        .get_column("GENE")
        .to_list()
    )

    frames: list[pl.LazyFrame] = []
    if gmm_files:
        any_significant = None
        for _, pval_column in METHOD_COLUMNS.values():
            condition = pl.col(pval_column) < threshold
            any_significant = (
                condition if any_significant is None else any_significant | condition
            )
        gmm = (
            pl.scan_parquet(
                str(gmm_dir / "*.parquet"), include_file_paths="source_file"
            )
            .filter(any_significant)
            .with_columns(
                pl.col("source_file")
                .str.extract(r"/([^/]+)\.parquet$", 1)
                .alias("GENE"),
                method_structs(),
            )
            .explode("association")
            .unnest("association")
            .filter(pl.col("my_pval") < threshold)
            .select("GENE", "RSID", "method", "my_beta", "my_pval")
        )
        frames.append(gmm)

    if fallback_genes:
        tar_path = study_dir / f"TAR_{cell_type}" / f"chr{chromosome}.csv"
        fallback = (
            pl.scan_csv(tar_path)
            .filter(
                pl.col("GENE").is_in(fallback_genes)
                & (pl.col("PVAL") < threshold)
            )
            .with_columns(method_structs("BETA", "PVAL"))
            .explode("association")
            .unnest("association")
            .select("GENE", "RSID", "method", "my_beta", "my_pval")
        )
        frames.append(fallback)

    if not frames:
        raise RuntimeError(f"No analyzable associations found in {gmm_dir}")

    info_path = study_dir / "INFO" / f"chr{chromosome}.csv"
    info = (
        pl.scan_csv(info_path)
        .select("GENE", "RSID", "POS", "A1", "A2")
        .with_columns(
            pl.lit(chromosome).cast(pl.Int16).alias("CHR"),
            pl.col("POS").cast(pl.Int64),
            pl.col("A1").str.to_uppercase().alias("my_A1"),
            pl.col("A2").str.to_uppercase().alias("my_A2"),
        )
        .drop("A1", "A2")
    )
    associations = (
        pl.concat(frames)
        .join(info, on=["GENE", "RSID"], how="inner")
        .unique(subset=["GENE", "RSID", "method"], maintain_order=False)
        .collect(engine="streaming")
    )

    observed = {
        row["method"]: int(row["len"])
        for row in associations.group_by("method").len().iter_rows(named=True)
    }
    observed = {method: observed.get(method, 0) for method in METHODS}
    expected = expected_counts(summary_path)
    if observed != expected:
        raise RuntimeError(
            f"eSNP count mismatch for {study_dir.name} chr{chromosome}: "
            f"observed={observed}, expected={expected}"
        )
    return associations, observed


def complement(column: str) -> pl.Expr:
    return pl.col(column).replace_strict(COMPLEMENT, default=None)


def harmonize_candidates(candidates: pl.DataFrame) -> pl.DataFrame:
    same = (
        (pl.col("my_A1") == pl.col("replicate_A1"))
        & (pl.col("my_A2") == pl.col("replicate_A2"))
    )
    swapped = (
        (pl.col("my_A1") == pl.col("replicate_A2"))
        & (pl.col("my_A2") == pl.col("replicate_A1"))
    )
    strand_same = (
        (pl.col("my_A1") == complement("replicate_A1"))
        & (pl.col("my_A2") == complement("replicate_A2"))
    )
    strand_swapped = (
        (pl.col("my_A1") == complement("replicate_A2"))
        & (pl.col("my_A2") == complement("replicate_A1"))
    )
    alignment = (
        pl.when(same)
        .then(pl.lit("same"))
        .when(swapped)
        .then(pl.lit("swapped"))
        .when(strand_same)
        .then(pl.lit("strand_same"))
        .when(strand_swapped)
        .then(pl.lit("strand_swapped"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_alignment")
    )
    result = candidates.with_columns(alignment).filter(
        pl.col("allele_alignment") != "mismatch"
    )
    result = result.with_columns(
        pl.when(pl.col("allele_alignment").is_in(["swapped", "strand_swapped"]))
        .then(-pl.col("replicate_beta"))
        .otherwise(pl.col("replicate_beta"))
        .alias("aligned_replicate_beta")
    )
    return result.with_columns(
        (
            (pl.col("my_beta") * pl.col("aligned_replicate_beta")) > 0
        ).alias("same_sign")
    )


def count_by_method(frame: pl.DataFrame) -> dict[str, int]:
    if frame.is_empty():
        return {method: 0 for method in METHODS}
    observed = {
        row["method"]: int(row["len"])
        for row in frame.group_by("method").len().iter_rows(named=True)
    }
    return {method: observed.get(method, 0) for method in METHODS}


def compare_with_replicate(
    my_esnps: pl.DataFrame,
    replicate_name: str,
    replicate_cache: Path,
    chromosome: int,
    replicate_cell_type: str,
) -> tuple[dict[str, dict[str, int]], pl.DataFrame]:
    replicate_query = pl.scan_parquet(replicate_cache).filter(
        pl.col("replicate_CHR") == chromosome
    )
    if replicate_name == "CIMA":
        replicate_query = replicate_query.filter(
            pl.col("replicate_cell_type") == replicate_cell_type
        )
    replicate = replicate_query.collect(engine="streaming")
    if replicate_name in REPLICATES:
        candidates = my_esnps.join(
            replicate,
            left_on=["GENE", "RSID"],
            right_on=["GENE", "replicate_RSID"],
            how="inner",
        )
    else:
        raise ValueError(replicate_name)

    candidates = candidates.unique(
        subset=["method", "GENE", "RSID", "replicate_variant_id"],
        maintain_order=False,
    )
    candidate_counts = count_by_method(candidates)
    harmonized = harmonize_candidates(candidates)
    harmonized_counts = count_by_method(harmonized)
    evaluable = harmonized.filter(
        pl.col("my_beta").is_finite()
        & pl.col("replicate_beta").is_finite()
        & (pl.col("my_beta") != 0)
        & (pl.col("replicate_beta") != 0)
    )
    evaluable_counts = count_by_method(evaluable)
    same_sign_counts = count_by_method(evaluable.filter(pl.col("same_sign")))
    counts = {
        method: {
            "candidate": candidate_counts[method],
            "harmonized": harmonized_counts[method],
            "evaluable": evaluable_counts[method],
            "same_sign": same_sign_counts[method],
        }
        for method in METHODS
    }
    return counts, harmonized


def add_detail_metadata(
    details: pl.DataFrame,
    qtdid: str,
    study_name: str,
    cell_type: str,
    replicate: str,
    replicate_cell_type: str,
    replicate_data_level: str,
    replicate_variant_scope: str,
) -> pl.DataFrame:
    if details.is_empty():
        return details
    return (
        details.with_columns(
            pl.lit(qtdid).alias("QTDid"),
            pl.lit(study_name).alias("study_name"),
            pl.lit(cell_type).alias("cell_type"),
            pl.lit(replicate).alias("replicate"),
            pl.lit(replicate_cell_type).alias("replicate_cell_type"),
            pl.lit(replicate_data_level).alias("replicate_data_level"),
            pl.lit(replicate_variant_scope).alias("replicate_variant_scope"),
        )
        .select(*DETAIL_COLUMNS)
        .sort("replicate", "method", "CHR", "POS", "GENE", "RSID")
    )


def write_detail_chunk(
    writer: pq.ParquetWriter | None,
    details: pl.DataFrame,
    path: Path,
) -> pq.ParquetWriter | None:
    if details.is_empty():
        return writer
    table = details.to_arrow()
    if writer is None:
        writer = pq.ParquetWriter(path, table.schema, compression="zstd")
    elif table.schema != writer.schema:
        table = table.cast(writer.schema)
    writer.write_table(table)
    return writer


def safe_ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def replicate_context(replicate: str, cell_type: str) -> tuple[str, str, str, str]:
    if replicate == "hum0197":
        return (
            f"hum0197:{cell_type}",
            cell_type,
            "celltype",
            "all_nominal_cis_eQTL_pairs",
        )
    if replicate == "hum0343":
        return (
            "hum0343",
            "Whole_blood",
            "tissue",
            "all_nominal_cis_eQTL_pairs",
        )
    if replicate == "CIMA":
        return (
            f"CIMA:{cell_type}",
            cell_type,
            "celltype",
            "lead_cis_eQTL_only",
        )
    raise ValueError(replicate)


def process_study(
    qtdid: str,
    study_root: Path,
    study_name: str,
    cell_type: str,
    chromosomes: Iterable[int],
    threshold: float,
    replicate_paths: dict[str, Path],
    replicate_counts: dict[str, int],
    replicate_mapped_counts: dict[str, int],
    work_detail_path: Path,
) -> list[dict]:
    study_dir = study_root / qtdid
    my_counts = {method: 0 for method in METHODS}
    comparison_counts = {
        replicate: {
            method: {"candidate": 0, "harmonized": 0, "evaluable": 0, "same_sign": 0}
            for method in METHODS
        }
        for replicate in REPLICATES
    }
    detail_writer: pq.ParquetWriter | None = None
    try:
        for chromosome in chromosomes:
            my_esnps, chromosome_counts = load_my_esnps(
                study_dir, cell_type, chromosome, threshold
            )
            for method in METHODS:
                my_counts[method] += chromosome_counts[method]

            for replicate in REPLICATES:
                (
                    cache_key,
                    replicate_cell_type,
                    replicate_data_level,
                    replicate_variant_scope,
                ) = replicate_context(replicate, cell_type)
                counts, details = compare_with_replicate(
                    my_esnps,
                    replicate,
                    replicate_paths[cache_key],
                    chromosome,
                    replicate_cell_type,
                )
                for method in METHODS:
                    for metric in comparison_counts[replicate][method]:
                        comparison_counts[replicate][method][metric] += counts[method][metric]
                details = add_detail_metadata(
                    details,
                    qtdid,
                    study_name,
                    cell_type,
                    replicate,
                    replicate_cell_type,
                    replicate_data_level,
                    replicate_variant_scope,
                )
                detail_writer = write_detail_chunk(
                    detail_writer, details, work_detail_path
                )
            log(
                f"{qtdid} chr{chromosome}: "
                + ", ".join(
                    f"{method}={chromosome_counts[method]:,}" for method in METHODS
                )
            )
    finally:
        if detail_writer is not None:
            detail_writer.close()

    rows: list[dict] = []
    for replicate in REPLICATES:
        (
            cache_key,
            replicate_cell_type,
            replicate_data_level,
            replicate_variant_scope,
        ) = replicate_context(replicate, cell_type)
        replicate_count = replicate_counts[cache_key]
        replicate_mapped_count = replicate_mapped_counts[cache_key]
        for method in METHODS:
            counts = comparison_counts[replicate][method]
            rows.append(
                {
                    "QTDid": qtdid,
                    "study_name": study_name,
                    "cell_type": cell_type,
                    "replicate": replicate,
                    "replicate_cell_type": replicate_cell_type,
                    "replicate_data_level": replicate_data_level,
                    "replicate_variant_scope": replicate_variant_scope,
                    "method": method,
                    "pvalue_threshold": threshold,
                    "replicate_esnp_count": replicate_count,
                    "replicate_rsid_mapped_esnp_count": replicate_mapped_count,
                    "my_esnp_count": my_counts[method],
                    "candidate_key_overlap_count": counts["candidate"],
                    "allele_mismatch_count": counts["candidate"] - counts["harmonized"],
                    "esnp_intersection_count": counts["harmonized"],
                    "sign_evaluable_intersection_count": counts["evaluable"],
                    "same_sign_count": counts["same_sign"],
                    "same_sign_proportion": safe_ratio(
                        counts["same_sign"], counts["evaluable"]
                    ),
                    "same_sign_percentage": (
                        100 * counts["same_sign"] / counts["evaluable"]
                        if counts["evaluable"]
                        else None
                    ),
                    "esnp_overlap_proportion_of_my_data": safe_ratio(
                        counts["harmonized"], my_counts[method]
                    ),
                    "esnp_overlap_proportion_of_replicate": safe_ratio(
                        counts["harmonized"], replicate_count
                    ),
                    "esnp_overlap_proportion_of_rsid_mapped_replicate": safe_ratio(
                        counts["harmonized"], replicate_mapped_count
                    ),
                }
            )
    log(f"Completed {qtdid}")
    return rows


def combine_detail_files(paths: list[Path], output_path: Path) -> int:
    existing = [path for path in paths if path.exists()]
    if not existing:
        raise RuntimeError("No harmonized intersection detail rows were produced")
    temporary = output_path.with_suffix(".tmp.parquet")
    (
        pl.scan_parquet([str(path) for path in existing])
        .sort("QTDid", "replicate", "method", "CHR", "POS", "GENE", "RSID")
        .sink_parquet(
            temporary,
            compression="zstd",
            row_group_size=100_000,
            statistics=True,
        )
    )
    os.replace(temporary, output_path)
    return int(pl.scan_parquet(output_path).select(pl.len()).collect().item())


def write_readme(output_dir: Path) -> None:
    content = """eSNP replication outputs
========================

Statistical unit
----------------
An eSNP is counted as one significant SNP-gene association, not one unique
rsID. The default significance threshold is p < 1e-5, matching traceCB's
TAR_*eSNP definition. Genes for which GMM was not run inherit the Original
association for traceC and traceCB, matching visual/utils.py::load_all_summary.

Replicates
----------
hum0197 uses the matching broad cell-level result: Mono, CD4T, CD8T, B, or NK.
hum0343 uses its whole-blood result for all studies. The hum0197 slope and the
hum0343 beta are effects of ALT; traceCB A1 is also the effect allele.
CIMA is celltype-level. Its L4 subtypes are grouped into the same five broad
cell types, and duplicate lead SNP-gene pairs within a broad group retain the
record with the smallest nominal p-value. CIMA slope uses its documented ALT
effect allele.

Scope warning
-------------
The downloaded CIMA resource contains lead cis-eQTL records only, whereas
hum0197 and hum0343 contain full nominal cis-eQTL results. CIMA counts and
overlaps therefore represent significant lead SNP-gene pairs and are not
directly comparable to the full eSNP counts from hum0197 or hum0343. See the
replicate_variant_scope column in the summary and detail outputs.

Allele and sign definitions
---------------------------
The intersection contains SNP-gene pairs whose alleles can be harmonized.
Same, swapped, strand-complemented, and strand-complemented/swapped pairs are
supported. A replicate beta is negated for swapped orientations before signs
are compared. same_sign_proportion uses sign_evaluable_intersection_count as
its denominator; this is normally identical to esnp_intersection_count.

Files
-----
esnp_replication_summary.csv: one study x method x replicate per row. For
hum0197 and CIMA it reports both all significant pairs and the identifier-
mapped subset used for intersection.
esnp_replication_intersections.parquet: auditable harmonized pair-level detail.
esnp_replication_run_metadata.json: input paths, parameters, and output counts.
_esnp_replication_cache/: filtered, deduplicated replicate data for fast reruns.
"""
    path = output_dir / "esnp_replication_README.txt"
    temporary = path.with_suffix(".txt.tmp")
    temporary.write_text(content)
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    if args.threshold <= 0 or args.threshold >= 1:
        raise ValueError("--threshold must be between 0 and 1")
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    chromosomes = sorted(set(args.chromosomes))
    if not chromosomes or chromosomes[0] < 1 or chromosomes[-1] > 22:
        raise ValueError("--chromosomes must contain values from 1 through 22")

    metadata = load_metadata()
    study_ids = args.study_ids or metadata["QTDids"]
    unknown = sorted(set(study_ids) - set(metadata["QTDids"]))
    if unknown:
        raise ValueError(f"Unknown QTDid values: {unknown}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = args.output_dir / "_esnp_replication_cache"
    replicate_paths, replicate_counts, replicate_mapped_counts = prepare_replicate_caches(
        args.hum0197_dir,
        args.hum0343_esnp,
        args.cima_lead_eqtl,
        args.gene_annotation,
        args.gtex_lookup,
        cache_dir,
        args.threshold,
        args.refresh_cache,
    )

    summary_rows: list[dict] = []
    temporary_work = Path(
        tempfile.mkdtemp(prefix=".esnp_replication_work_", dir=args.output_dir)
    )
    detail_paths = {qtdid: temporary_work / f"{qtdid}.parquet" for qtdid in study_ids}
    try:
        worker_count = min(args.workers, len(study_ids))
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(
                    process_study,
                    qtdid,
                    args.study_dir,
                    metadata["id2name"][qtdid],
                    metadata["id2celltype"][qtdid],
                    chromosomes,
                    args.threshold,
                    replicate_paths,
                    replicate_counts,
                    replicate_mapped_counts,
                    detail_paths[qtdid],
                ): qtdid
                for qtdid in study_ids
            }
            for future in as_completed(futures):
                summary_rows.extend(future.result())

        summary = pl.DataFrame(summary_rows).sort("QTDid", "replicate", "method")
        summary_path = args.output_dir / "esnp_replication_summary.csv"
        temporary_summary = summary_path.with_suffix(".tmp.csv")
        summary.write_csv(temporary_summary, float_scientific=False)
        os.replace(temporary_summary, summary_path)

        detail_path = args.output_dir / "esnp_replication_intersections.parquet"
        detail_count = combine_detail_files(list(detail_paths.values()), detail_path)
        write_readme(args.output_dir)

        run_metadata = {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "study_dir": str(args.study_dir.resolve()),
            "hum0197_dir": str(args.hum0197_dir.resolve()),
            "hum0343_esnp": str(args.hum0343_esnp.resolve()),
            "cima_lead_eqtl": str(args.cima_lead_eqtl.resolve()),
            "gene_annotation": str(args.gene_annotation.resolve()),
            "gtex_lookup": str(args.gtex_lookup.resolve()),
            "output_dir": str(args.output_dir.resolve()),
            "pvalue_threshold": args.threshold,
            "study_ids": study_ids,
            "chromosomes": chromosomes,
            "workers": worker_count,
            "hum0197_celltype_mapping": CELLTYPE_TO_HUM0197,
            "cima_celltype_mapping": (
                "L4 subtypes mapped to Monocytes, CD4+T_cells, CD8+T_cells, "
                "B_cells, or NK_cells; mixed/adjacent populations excluded"
            ),
            "replicate_esnp_counts": replicate_counts,
            "replicate_rsid_mapped_esnp_counts": replicate_mapped_counts,
            "summary_row_count": summary.height,
            "intersection_detail_row_count": detail_count,
        }
        write_json_atomic(
            run_metadata, args.output_dir / "esnp_replication_run_metadata.json"
        )
        log(f"Saved {summary_path}")
        log(f"Saved {detail_path} ({detail_count:,} rows)")
    finally:
        shutil.rmtree(temporary_work, ignore_errors=True)


if __name__ == "__main__":
    main()
