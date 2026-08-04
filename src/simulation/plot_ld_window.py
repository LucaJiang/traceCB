#!/usr/bin/env python
"""Visualize the LD window used by the small-window simulations."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EAS = ROOT / "data/simulation/EAS_n5000_chr22_loci29.npy"
DEFAULT_EUR = ROOT / "data/simulation/EUR_n20000_chr22_loci29.npy"
DEFAULT_EAS_TXT = ROOT / "data/simulation/EAS_n5000_chr22_loci29.txt"
DEFAULT_OUT_DIR = ROOT / "results/simulation/ld_window"
DEFAULT_BIMS = [
    ROOT / "data/simulation/1000G/1000G.EAS.QC.maf.22.bim",
    ROOT / "data/simulation/1000G/1000G.EUR.QC.maf.22.bim",
    ROOT / "data/simulation/ukb/height_ukb_50k_chr22.bim",
    ROOT / "data/simulation/wg_ukb/height_merge_qc2_chr22.bim",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pop1-geno", type=Path, default=DEFAULT_EAS)
    parser.add_argument("--pop2-geno", type=Path, default=DEFAULT_EUR)
    parser.add_argument("--snp-txt", type=Path, default=DEFAULT_EAS_TXT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--start", type=int, default=1000)
    parser.add_argument("--end", type=int, default=3000)
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def read_snp_ids(path: Path) -> list[str]:
    with path.open() as handle:
        return handle.readline().strip().split()


def read_positions(bim_paths: list[Path]) -> dict[str, int]:
    positions: dict[str, int] = {}
    for bim_path in bim_paths:
        if not bim_path.exists():
            continue
        with bim_path.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 4:
                    continue
                snp_id = fields[1]
                if snp_id in positions:
                    continue
                try:
                    positions[snp_id] = int(fields[3])
                except ValueError:
                    continue
    return positions


def standardize(genotype: np.ndarray) -> np.ndarray:
    x = np.asarray(genotype, dtype=np.float32).copy()
    x -= x.mean(axis=0, dtype=np.float64).astype(np.float32)
    sd = x.std(axis=0, dtype=np.float64).astype(np.float32)
    sd[sd == 0] = 1.0
    x /= sd
    return x


def corr_matrix(genotype_path: Path, start: int, end: int) -> np.ndarray:
    genotype = np.load(genotype_path, mmap_mode="r")
    if end > genotype.shape[1]:
        raise ValueError(f"Window end {end} exceeds {genotype_path} SNP count {genotype.shape[1]}")
    window = np.asarray(genotype[:, start:end], dtype=np.float32)
    x = standardize(window)
    corr = (x.T @ x) / np.float32(x.shape[0])
    corr = np.clip(corr, -1.0, 1.0).astype(np.float32, copy=False)
    np.fill_diagonal(corr, 1.0)
    return corr


def offdiag_values(matrix: np.ndarray) -> np.ndarray:
    mask = ~np.eye(matrix.shape[0], dtype=bool)
    return matrix[mask]


def add_heatmap(ax, matrix: np.ndarray, title: str, cmap: str, vmin: float, vmax: float, start: int, end: int):
    im = ax.imshow(
        matrix,
        origin="upper",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[start, end, end, start],
        rasterized=True,
    )
    ax.set_title(title)
    ticks = [start, start + 500, start + 1000, start + 1500, end - 1]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel("Global SNP index")
    ax.set_ylabel("Global SNP index")
    for spine in ax.spines.values():
        spine.set_edgecolor("#c1121f")
        spine.set_linewidth(1.2)
    return im


def plot_ld(
    snp_ids: list[str],
    positions: dict[str, int],
    start: int,
    end: int,
    corr1: np.ndarray,
    corr2: np.ndarray,
    out_path: Path,
    dpi: int,
) -> None:
    r2_pop1 = corr1 * corr1
    r2_pop2 = corr2 * corr2
    cross = corr1 * corr2

    all_pos = np.array([positions.get(snp_id, np.nan) for snp_id in snp_ids], dtype=float) / 1_000_000
    window_idx = np.arange(start, end)
    window_pos = all_pos[start:end]

    fig = plt.figure(figsize=(18, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(nrows=2, ncols=3, height_ratios=[0.36, 1.0])
    context_ax = fig.add_subplot(grid[0, :])
    axes = [fig.add_subplot(grid[1, i]) for i in range(3)]

    context_ax.scatter(np.arange(len(snp_ids)), all_pos, s=3, color="#9aa0a6", alpha=0.5, linewidths=0)
    context_ax.scatter(window_idx, window_pos, s=4, color="#c1121f", alpha=0.9, linewidths=0)
    context_ax.axvspan(start, end - 1, color="#c1121f", alpha=0.08)
    context_ax.axvline(start, color="#c1121f", linestyle="--", linewidth=1.0)
    context_ax.axvline(end - 1, color="#c1121f", linestyle="--", linewidth=1.0)
    context_ax.set_xlim(0, len(snp_ids) - 1)
    context_ax.set_ylabel("Position (Mb)")
    context_ax.set_xlabel("Global SNP index in regenerated simulation panel")
    context_ax.set_title(
        f"Simulation LD window: columns {start}:{end} "
        f"({end - start} SNPs; 1-based ordinal {start + 1}-{end})"
    )
    context_ax.text(
        start,
        np.nanmax(all_pos),
        f"{snp_ids[start]}\n{int(positions[snp_ids[start]]):,}",
        ha="left",
        va="top",
        fontsize=8,
        color="#6f1d1b",
    )
    context_ax.text(
        end - 1,
        np.nanmax(all_pos),
        f"{snp_ids[end - 1]}\n{int(positions[snp_ids[end - 1]]):,}",
        ha="right",
        va="top",
        fontsize=8,
        color="#6f1d1b",
    )

    im1 = add_heatmap(axes[0], r2_pop1, "EAS LD ($r^2$)", "viridis", 0.0, 1.0, start, end)
    im2 = add_heatmap(axes[1], r2_pop2, "EUR LD ($r^2$)", "viridis", 0.0, 1.0, start, end)
    im3 = add_heatmap(
        axes[2],
        cross,
        "Cross-pop shared LD ($r_{EAS} \\times r_{EUR}$)",
        "coolwarm",
        -1.0,
        1.0,
        start,
        end,
    )
    fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def write_summary(
    path: Path,
    snp_ids: list[str],
    positions: dict[str, int],
    start: int,
    end: int,
    corr1: np.ndarray,
    corr2: np.ndarray,
) -> None:
    r2_pop1 = corr1 * corr1
    r2_pop2 = corr2 * corr2
    cross = corr1 * corr2
    rows = [
        ("total_snps", len(snp_ids)),
        ("window_start_0_based", start),
        ("window_end_exclusive_0_based", end),
        ("window_snp_count", end - start),
        ("window_ordinal_1_based", f"{start + 1}-{end}"),
        ("first_window_snp", snp_ids[start]),
        ("first_window_pos", positions[snp_ids[start]]),
        ("last_window_snp", snp_ids[end - 1]),
        ("last_window_pos", positions[snp_ids[end - 1]]),
        ("window_pos_min", min(positions[snp_id] for snp_id in snp_ids[start:end])),
        ("window_pos_max", max(positions[snp_id] for snp_id in snp_ids[start:end])),
    ]
    for label, matrix in [
        ("eas_r2", r2_pop1),
        ("eur_r2", r2_pop2),
        ("cross_r_product", cross),
    ]:
        values = offdiag_values(matrix)
        rows.extend(
            [
                (f"{label}_offdiag_mean", float(np.mean(values))),
                (f"{label}_offdiag_p95", float(np.percentile(values, 95))),
                (f"{label}_offdiag_p99", float(np.percentile(values, 99))),
                (f"{label}_offdiag_max", float(np.max(values))),
                (f"{label}_offdiag_min", float(np.min(values))),
            ]
        )
    with path.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    snp_ids = read_snp_ids(args.snp_txt)
    if args.start < 0 or args.end <= args.start or args.end > len(snp_ids):
        raise ValueError(f"Invalid window {args.start}:{args.end} for {len(snp_ids)} SNPs")
    positions = read_positions(DEFAULT_BIMS)
    missing_pos = [snp_id for snp_id in snp_ids[args.start : args.end] if snp_id not in positions]
    if missing_pos:
        raise ValueError(f"Missing positions for {len(missing_pos)} SNPs, first: {missing_pos[:5]}")

    print(f"Loading and correlating EAS window {args.start}:{args.end}")
    corr1 = corr_matrix(args.pop1_geno, args.start, args.end)
    print(f"Loading and correlating EUR window {args.start}:{args.end}")
    corr2 = corr_matrix(args.pop2_geno, args.start, args.end)

    figure_path = args.out_dir / f"simulation_ld_window_{args.start}_{args.end}.png"
    summary_path = args.out_dir / f"simulation_ld_window_{args.start}_{args.end}_summary.tsv"
    plot_ld(snp_ids, positions, args.start, args.end, corr1, corr2, figure_path, args.dpi)
    write_summary(summary_path, snp_ids, positions, args.start, args.end, corr1, corr2)
    print(f"figure={figure_path}")
    print(f"summary={summary_path}")


if __name__ == "__main__":
    main()
