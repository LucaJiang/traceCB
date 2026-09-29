"""Shared figure paths exported by ``source scripts/config.sh``.

Repository-relative fallbacks allow imports without access to the full datasets.
Input loaders validate the files they actually need, with configuration guidance.
"""

import os
from pathlib import Path


def configured_path(variable, default):
    return Path(os.environ.get(variable) or default).expanduser().resolve()


REPO_ROOT = configured_path("TRACECB_REPO_ROOT", Path(__file__).resolve().parents[2])
DATA_ROOT = configured_path("TRACECB_DATA_ROOT", REPO_ROOT / "data")
OUTPUT_ROOT = configured_path("TRACECB_OUTPUT_ROOT", REPO_ROOT / "results")
STUDY_ROOT = configured_path("TRACECB_STUDY_ROOT", DATA_ROOT / "traceCB")
POPULATION = os.environ.get("TARGET_POPULATION", "EAS")
TISSUE_SOURCE = os.environ.get("TISSUE_SOURCE", "eQTLGen")
STUDY_DIR = configured_path(
    "TRACECB_STUDY_DIR",
    STUDY_ROOT / f"{POPULATION}_{TISSUE_SOURCE}",
)
FIGURE_DIR = configured_path(
    "TRACECB_FIGURE_DIR", OUTPUT_ROOT / "figures" / f"{POPULATION}_{TISSUE_SOURCE}"
)
METADATA_FILE = configured_path(
    "TRACECB_FIGURE_METADATA",
    REPO_ROOT / "src/figures/metadata.json",
)
GENE_ANNOTATION = configured_path(
    "TRACECB_GTEX_GENE_ANNOTATION",
    DATA_ROOT / "GTEx/gencode.v26.GRCh38.genes.gtf",
)
GTEX_LOOKUP = configured_path(
    "TRACECB_GTEX_LOOKUP",
    DATA_ROOT / "GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz",
)
ONEK1K_FILE = configured_path(
    "TRACECB_ONEK1K_FILE",
    STUDY_ROOT / "onek1k_supp/onek1k_esnp.csv",
)
OASIS_DIR = configured_path(
    "TRACECB_OASIS_DIR",
    DATA_ROOT / "hum0197/eQTL_summary_statistics",
)
CIMA_DIR = configured_path("TRACECB_CIMA_DIR", DATA_ROOT / "CIMA")
CIMA_LEAD_EQTL = configured_path(
    "TRACECB_CIMA_LEAD_EQTL",
    CIMA_DIR / "xQTL/CIMA_Lead_cis-xQTL.csv",
)
CIMA_CELL_TYPES = configured_path(
    "TRACECB_CIMA_CELL_TYPES",
    CIMA_DIR / "Cell_Atlas/CIMA_Cell_Type_Level_and_Marker.xlsx",
)
REPLICATION_EGENES = configured_path(
    "TRACECB_REPLICATION_EGENES",
    DATA_ROOT / "hum0343/hum0343_eGene.csv",
)
REPLICATION_ESNPS = configured_path(
    "TRACECB_REPLICATION_ESNPS",
    DATA_ROOT / "hum0343/hum0343_eSNP.csv",
)
CELL_PROPORTIONS = configured_path(
    "TRACECB_CELL_PROPORTIONS",
    STUDY_ROOT / "cell_type_proportion/ind_celltype_proportion.csv",
)
AFR_STUDY_DIR = configured_path(
    "TRACECB_AFR_STUDY_DIR",
    STUDY_ROOT / f"AFR_{TISSUE_SOURCE}",
)
AFR_FIGURE_DIR = configured_path("TRACECB_AFR_FIGURE_DIR", FIGURE_DIR / "afr")
ESNP_REPLICATION_DIR = configured_path(
    "TRACECB_ESNP_REPLICATION_DIR",
    FIGURE_DIR / "esnp_replication",
)
TIMING_FILE = configured_path(
    "TRACECB_TIMING_FILE",
    REPO_ROOT / "tmp/timing/summary_timing.csv",
)
COLOC_INPUT_DIR = configured_path("TRACECB_COLOC_INPUT_DIR", STUDY_ROOT / "coloc")
COLOC_DIR = configured_path(
    "TRACECB_COLOC_DIR",
    OUTPUT_ROOT / f"{POPULATION}_{TISSUE_SOURCE}/coloc",
)
COLOC_GENES = configured_path(
    "TRACECB_COLOC_GENES",
    COLOC_INPUT_DIR / "bcx/bcx_mon.closest.protein_coding.bed",
)
COLOC_REPLICATION = configured_path(
    "TRACECB_COLOC_REPLICATION",
    COLOC_DIR / "replication.csv",
)
COLOC_FIGURE_DIR = configured_path("TRACECB_COLOC_FIGURE_DIR", FIGURE_DIR / "coloc")


def require_files(paths, variable):
    """Reject incomplete datasets before treating absent records as negatives."""
    missing = [str(path) for path in paths if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(
            f"Missing input file(s) for {variable}:\n  " + "\n  ".join(missing)
            + "\nRun 'source scripts/config.sh' from the repository root, then rerun "
            "this Python script in the same shell. If files are still missing, "
            f"check {variable} in scripts/config.sh."
        )


def require_file(path, variable):
    require_files([path], variable)
    return path
