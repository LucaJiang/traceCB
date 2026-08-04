"""Summarize posterior-probability thresholds across colocalization results."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


COUNT_COLUMNS = (
    "p_original",
    "p_traceC",
    "p_traceCB",
    "pp_h3_original",
    "pp_h3_traceC",
    "pp_h3_traceCB",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing *_coloc.csv result files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Summary CSV path. Defaults to <input-dir>/coloc_summary.csv.",
    )
    parser.add_argument("--threshold", type=float, default=0.7)
    return parser.parse_args()


def summarize(input_dir: Path, threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path in sorted(input_dir.glob("*_coloc.csv")):
        frame = pd.read_csv(path)
        missing = set(COUNT_COLUMNS) - set(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        row: dict[str, object] = {
            "result": path.name.removesuffix("_coloc.csv"),
            "total_loci": len(frame),
        }
        row.update(
            {
                f"{column}_gt_{threshold:g}": int((frame[column] > threshold).sum())
                for column in COUNT_COLUMNS
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    output = args.output or args.input_dir / "coloc_summary.csv"
    result = summarize(args.input_dir, args.threshold)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(f"Wrote {len(result)} result summaries to {output}")


if __name__ == "__main__":
    main()
