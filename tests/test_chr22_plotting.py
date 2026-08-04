import importlib.util
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "simulation"
    / "chr22"
    / "plot_results.py"
)
spec = importlib.util.spec_from_file_location("chr22_plot_results", MODULE_PATH)
chr22_plot = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(chr22_plot)


def test_default_error_bars_use_gene_level_repeats():
    args = chr22_plot.parse_args([])

    assert args.error_unit == "gene"


def test_gene_error_bars_use_genes_not_replicates_as_units():
    df = pd.DataFrame(
        {
            "h1sq": [0.1, 0.1, 0.1, 0.1],
            "h2sq": [0.1, 0.1, 0.1, 0.1],
            "gc": [0.7, 0.7, 0.7, 0.7],
            "n1": [100, 100, 100, 100],
            "n2": [400, 400, 400, 400],
            "nt": [5000, 5000, 5000, 5000],
            "propt": [0.2, 0.2, 0.2, 0.2],
            "pcausal": [0.005, 0.005, 0.005, 0.005],
            "method": ["traceCB", "traceCB", "traceCB", "traceCB"],
            "gene_id": ["gene_a", "gene_a", "gene_b", "gene_b"],
            "rep": [0, 1, 0, 1],
            "nsnp_gene": [10, 10, 20, 20],
            "power_shared": [0.2, 0.4, 0.6, 0.8],
        }
    )

    summary = chr22_plot.summarize_metric(df, "power_shared", "gene_mean", "gene")

    assert summary.shape[0] == 2
    assert summary["gene_id"].tolist() == ["gene_a", "gene_b"]
    assert summary["value"].tolist() == pytest.approx([0.3, 0.7])


def test_chr22_visualizer_uses_seaborn_ci_instead_of_t_critical_table():
    assert not hasattr(chr22_plot, "T_CRIT_975")
    assert not hasattr(chr22_plot, "t_critical_975")


def test_alpha_null_default_ymax_matches_chr22_plot_target():
    assert chr22_plot.METRIC_DEFAULT_YMAX["alpha_null"] == 0.4


def test_errorbar_option_accepts_sd():
    args = chr22_plot.parse_args(["--errorbar", "sd"])

    assert args.errorbar == "sd"


def test_barplot_errorbar_uses_requested_mode(monkeypatch):
    seen_errorbars = []

    def fake_barplot(**kwargs):
        seen_errorbars.append(kwargs["errorbar"])
        ax = kwargs["ax"]
        ax.bar([0], [0.5], label="traceCB")

    monkeypatch.setattr(chr22_plot.sns, "barplot", fake_barplot)
    plot_df = pd.DataFrame(
        {
            "run_prefix": ["baseline", "baseline"],
            "h2sq": [0.1, 0.1],
            "n2": [100, 100],
            "nt": [5000, 5000],
            "propt": [0.2, 0.2],
            "method": ["traceCB", "traceCB"],
            "value": [0.2, 0.8],
            "weight": [1.0, 1.0],
        }
    )

    fig, ax = chr22_plot.plt.subplots()
    chr22_plot.plot_metric_axis(
        ax,
        plot_df,
        "power_shared",
        ["baseline"],
        0.0,
        1.0,
        errorbar="sd",
    )
    chr22_plot.plt.close(fig)

    assert seen_errorbars == ["sd"]
