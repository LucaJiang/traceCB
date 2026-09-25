import importlib.util
from pathlib import Path

import polars as pl


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "src/figures/esnp_replication.py"
SPEC = importlib.util.spec_from_file_location("esnp_replication", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_harmonize_candidates_handles_allele_orientation_and_strand():
    candidates = pl.DataFrame(
        {
            "case": ["same", "swapped", "strand_same", "strand_swapped", "mismatch"],
            "my_A1": ["A", "G", "T", "C", "A"],
            "my_A2": ["G", "A", "C", "T", "C"],
            "replicate_A1": ["A", "A", "A", "A", "A"],
            "replicate_A2": ["G", "G", "G", "G", "T"],
            "my_beta": [0.2, 0.3, 0.2, 0.3, 0.4],
            "replicate_beta": [0.1, -0.2, 0.1, -0.4, 0.3],
        }
    )

    harmonized = MODULE.harmonize_candidates(candidates).sort("case")
    observed = {
        row["case"]: row
        for row in harmonized.select(
            "case", "allele_alignment", "aligned_replicate_beta", "same_sign"
        ).iter_rows(named=True)
    }

    assert set(observed) == {"same", "swapped", "strand_same", "strand_swapped"}
    assert observed["same"]["allele_alignment"] == "same"
    assert observed["swapped"]["allele_alignment"] == "swapped"
    assert observed["strand_same"]["allele_alignment"] == "strand_same"
    assert observed["strand_swapped"]["allele_alignment"] == "strand_swapped"
    assert observed["swapped"]["aligned_replicate_beta"] == 0.2
    assert observed["strand_swapped"]["aligned_replicate_beta"] == 0.4
    assert all(row["same_sign"] for row in observed.values())
