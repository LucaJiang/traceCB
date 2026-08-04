import importlib.util
from pathlib import Path
import tomllib

import numpy as np
import pandas as pd

import traceCB
from traceCB.gmm import GMM
from traceCB.ldsc import Run_Cross_LDSC


REPO_ROOT = Path(__file__).resolve().parents[1]
TOY_DIR = REPO_ROOT / "data" / "toy_example"
PREPARE_LOCI_PATH = REPO_ROOT / "src" / "coloc" / "prepare_loci.py"


def load_prepare_loci_module():
    spec = importlib.util.spec_from_file_location("prepare_loci", PREPARE_LOCI_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_package_version_matches_project_metadata():
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        project = tomllib.load(handle)
    assert traceCB.__version__ == project["project"]["version"]


def test_toy_data_runs_cross_population_model():
    target_gene = "ENSG00000025708"
    eas = pd.read_csv(TOY_DIR / "eas_summary_statistics.csv")
    eur = pd.read_csv(TOY_DIR / "eur_summary_statistics.csv")
    eas_ld = pd.read_csv(TOY_DIR / "eas_ld.csv")
    eur_ld = pd.read_csv(TOY_DIR / "eur_ld.csv")
    cross_ld = pd.read_csv(TOY_DIR / "cross_ld.csv")

    aligned = (
        eas.merge(
            eur.drop(columns=["GENE"]),
            on="RSID",
            suffixes=("_tar", "_aux"),
        )
        .merge(
            eas_ld.rename(columns={target_gene: "LD_tar", "SNP": "RSID"})[
                ["RSID", "LD_tar"]
            ],
            on="RSID",
        )
        .merge(
            eur_ld.rename(columns={target_gene: "LD_aux", "SNP": "RSID"})[
                ["RSID", "LD_aux"]
            ],
            on="RSID",
        )
        .merge(
            cross_ld.rename(columns={target_gene: "LD_x", "SNP": "RSID"})[
                ["RSID", "LD_x"]
            ],
            on="RSID",
        )
    )

    omega, omega_se = Run_Cross_LDSC(
        (aligned.BETA_tar / aligned.SE_tar).to_numpy(),
        aligned.N_tar.to_numpy(),
        aligned.LD_tar.to_numpy(),
        (aligned.BETA_aux / aligned.SE_aux).to_numpy(),
        aligned.N_aux.to_numpy(),
        aligned.LD_aux.to_numpy(),
        aligned.LD_x.to_numpy(),
        np.array([1.0, 1.0, 0.0]),
    )

    assert len(aligned) == 2072
    assert np.all(np.isfinite(omega))
    assert np.all(np.isfinite(omega_se))
    assert np.allclose(omega, omega.T)

    row = aligned.iloc[0]
    estimates = GMM(
        omega,
        np.eye(2),
        row.BETA_tar,
        row.SE_tar,
        row.LD_tar,
        row.BETA_aux,
        row.SE_aux,
        row.LD_aux,
        row.LD_x,
    )
    assert np.all(np.isfinite(estimates))
    assert estimates[1] > 0
    assert estimates[3] > 0


def test_lead_variant_bed_uses_valid_half_open_coordinates(tmp_path):
    prepare_loci = load_prepare_loci_module()
    ldlink_path = tmp_path / "trait_ldlink.tsv"
    pd.DataFrame(
        {
            "GWAS_snp": ["rs1"],
            "GWAS_snp_pos": ["chr1:100"],
            "locus_name": ["chr1:90_110"],
            "chrom": ["chr1"],
            "start": [90],
            "end": [110],
        }
    ).to_csv(ldlink_path, sep="\t", index=False)

    prepare_loci.create_index_snp_bed(ldlink_path, "trait", tmp_path)

    bed = pd.read_csv(
        tmp_path / "trait_GWAS_index_snps.bed", sep="\t", header=None
    )
    assert bed.iloc[0].tolist() == ["chr1", 99, 100, "rs1", 0]
