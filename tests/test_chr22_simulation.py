import importlib.util
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "simulation"
    / "chr22"
    / "simulate.py"
)
spec = importlib.util.spec_from_file_location("chr22_simulate", MODULE_PATH)
chr22_simulate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(chr22_simulate)


def test_grid_cli_accepts_multi_value_parameters_with_grid_defaults():
    args = chr22_simulate.parse_args(
        [
            "--run_prefix",
            "mix",
            "--h2sq",
            "0.1",
            "0.2",
            "--n2",
            "100",
            "400",
        ]
    )

    assert args.h2sq == [0.1, 0.2]
    assert args.n2 == [100, 400]
    assert args.out_dir == "bench/result/chr22_eqtl_mixture"
    assert args.max_snps_per_gene == 0
    assert args.estimate_omega is True
    assert args.architecture_probs == [0.25, 0.25, 0.25, 0.25]
    assert chr22_simulate.should_use_grid_runname(args)


def test_single_value_cli_stays_single_run_with_original_defaults():
    args = chr22_simulate.parse_args(["--runname", "one", "--h1sq", "0.2"])

    assert not chr22_simulate.should_use_grid_runname(args)
    single_args = chr22_simulate.build_setting_args(
        args,
        h1sq=0.2,
        h2sq=0.1,
        gc=0.7,
        n1=100,
        n2=400,
        nt=1000,
        propt=0.2,
        pcausal=0.005,
    )
    assert single_args.runname == "one"
    assert single_args.h1sq == 0.2
    assert single_args.out_dir == "bench/result"
    assert single_args.max_snps_per_gene == 1000
    assert single_args.estimate_omega is False


def test_grid_run_names_are_stable_across_parameter_settings():
    args = chr22_simulate.parse_args(
        [
            "--run_prefix",
            "mix",
            "--h1sq",
            "0.1",
            "--h2sq",
            "0.1",
            "0.2",
            "--gc",
            "0.7",
            "--n1",
            "100",
            "--n2",
            "100",
            "--nt",
            "5000",
            "--propt",
            "0.01",
            "--pcausal",
            "0.005",
            "--estimate_omega",
        ]
    )

    runs = list(chr22_simulate.iter_setting_args(args))

    assert [run.runname for run in runs] == [
        (
            "mix_h1sq_0.1_h2sq_0.1_gc_0.7_n1_100_n2_100_"
            "nt_5000_propt_0.01_pcausal_0.005_omega_False"
        ),
        (
            "mix_h1sq_0.1_h2sq_0.2_gc_0.7_n1_100_n2_100_"
            "nt_5000_propt_0.01_pcausal_0.005_omega_False"
        ),
    ]
    assert [run.h2sq for run in runs] == [0.1, 0.2]
    assert all(run.dry_run is False for run in runs)


def test_gene_architectures_are_randomly_assigned_from_probabilities():
    groups = [(f"gene_{i}", [{"gene_start": i, "gene_end": i + 1}]) for i in range(8)]

    architecture_map = chr22_simulate.assign_gene_architectures(
        groups, [0.25, 0.25, 0.25, 0.25], seed=7
    )

    counts = {
        architecture: list(architecture_map.values()).count(architecture)
        for architecture in chr22_simulate.ARCHITECTURES
    }

    assert sum(counts.values()) == 8
    assert counts != {architecture: 2 for architecture in chr22_simulate.ARCHITECTURES}


def test_gene_snp_count_summary_has_one_row_per_gene():
    groups = [
        (
            "gene_a",
            [
                {"gene_name": "A", "gene_type": "protein_coding"},
                {"gene_name": "A", "gene_type": "protein_coding"},
            ],
        ),
        ("gene_b", [{"gene_name": "B", "gene_type": "lncRNA"}]),
    ]

    rows = chr22_simulate.summarize_gene_snp_counts(groups)

    assert rows == [
        {
            "gene_id": "gene_a",
            "gene_name": "A",
            "gene_type": "protein_coding",
            "nsnp_gene": 2,
        },
        {
            "gene_id": "gene_b",
            "gene_name": "B",
            "gene_type": "lncRNA",
            "nsnp_gene": 1,
        },
    ]


def test_pop1_specific_genes_do_not_generate_tissue_effects():
    assert not chr22_simulate.has_tissue_genetic_effect("pop1_specific")
    assert not chr22_simulate.has_tissue_genetic_effect("null")
    assert chr22_simulate.has_tissue_genetic_effect("shared")
    assert chr22_simulate.has_tissue_genetic_effect("pop2_specific")
