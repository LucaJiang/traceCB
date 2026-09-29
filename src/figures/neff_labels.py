"""Shared, non-overlapping labels for the three traceCB Neff comparisons.

Save this file as src/figures/neff_labels.py. Only label candidates are filtered;
the scatter-point data are never modified. Call label_neff_scatter AFTER the last
layout/axis-limit change. Uses NumPy, pandas and Matplotlib; no adjustText needed.
"""
from typing import Callable, Optional, Sequence

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.transforms import Bbox

PAIRS = (
    ("TAR_SNEFF", "TAR_CNEFF"),
    ("TAR_SNEFF", "TAR_TNEFF"),
    ("TAR_CNEFF", "TAR_TNEFF"),
)
NEFF_COLS = ("TAR_SNEFF", "TAR_CNEFF", "TAR_TNEFF")
RATIO_COLS = ("traceC_Original", "traceCB_Original", "traceCB_traceC")


def select_shared_genes(
    df: pd.DataFrame,
    gene_name: Callable[[str], Optional[str]],
    max_genes: int = 5,
    ranking_mode: str = "all_three",
) -> pd.DataFrame:
    """Rank one study's shared labels by selection_score, largest first.

    all_three (default): min(traceC/Original, traceCB/Original, traceCB/traceC).
    first_two: min(traceC/Original, traceCB/Original); traceCB/traceC does not
    affect ranking or ratio eligibility, including when it is <=1 or infinite.
    Tie-break: traceCB/Original, then GENE. Each GENE must identify one point.
    The result has at most five distinct gene symbols; unknown symbols are skipped.
    Coordinates must be finite and positive; only ratios used by ranking_mode
    must be finite. These checks affect annotation candidates, never scatter data.
    min_ratio always records the minimum of all three ratios, independently of
    ranking_mode and selection_score. All three ratios are retained in the result.
    """
    if not isinstance(max_genes, (int, np.integer)) or not 0 <= max_genes <= 5:
        raise ValueError("max_genes must be an integer between 0 and 5")
    if ranking_mode not in ("all_three", "first_two"):
        raise ValueError("ranking_mode must be 'all_three' or 'first_two'")
    ranking_ratios = RATIO_COLS if ranking_mode == "all_three" else RATIO_COLS[:2]
    required = ["GENE", *NEFF_COLS]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    work = df.loc[:, required].copy()
    work = work.loc[work["GENE"].notna()].copy()
    work["GENE"] = work["GENE"].astype(str)
    for col in NEFF_COLS:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    values = work[list(NEFF_COLS)].to_numpy(dtype=float, na_value=np.nan)
    valid = np.isfinite(values).all(axis=1) & (values > 0).all(axis=1)
    work = work.loc[valid].drop_duplicates(required).copy()
    if work["GENE"].duplicated().any():
        raise ValueError("Conflicting rows for the same GENE; pass ONE study at a time")
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for (x_col, y_col), ratio_col in zip(PAIRS, RATIO_COLS):
            work[ratio_col] = work[y_col] / work[x_col]
    work = work.loc[np.isfinite(work[list(ranking_ratios)]).all(axis=1)].copy()
    work["min_ratio"] = work[list(RATIO_COLS)].min(axis=1)
    work["ranking_mode"] = ranking_mode
    work["selection_score"] = work[list(ranking_ratios)].min(axis=1)
    work = work.sort_values(
        ["selection_score", "traceCB_Original", "GENE"],
        ascending=[False, False, True], kind="mergesort",
    )
    # Lookup in ranking order; do not choose five and THEN discard unknown names.
    rows, seen_symbols = [], set()
    if max_genes:
        for _, row in work.iterrows():
            symbol = gene_name(row["GENE"])
            if not isinstance(symbol, str) or not symbol.strip():
                continue
            symbol = symbol.strip()
            if symbol in seen_symbols:
                continue
            row = row.copy()
            row["GENE_NAME"] = symbol
            rows.append(row)
            seen_symbols.add(symbol)
            if len(rows) == max_genes:
                break
    return pd.DataFrame(rows, columns=[*work.columns, "GENE_NAME"]).reset_index(drop=True)


def _contains(outer: Bbox, inner: Bbox) -> bool:
    return (inner.x0 >= outer.x0 and inner.x1 <= outer.x1
            and inner.y0 >= outer.y0 and inner.y1 <= outer.y1)


def _positions(ax, selected, x_col, y_col, points, fontsize, gap_pt):
    """Deterministic, greedy placement using measured glyph bounding boxes."""
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    px_per_pt = fig.dpi / 72.0
    pad = gap_pt * px_per_pt
    area = ax.get_window_extent(renderer).padded(-2 * px_per_pt)
    targets = ax.transData.transform(selected[[x_col, y_col]].to_numpy(float))
    obstacles = []
    legend = ax.get_legend()
    if legend is not None and legend.get_visible():
        obstacles.append(legend.get_window_extent(renderer).padded(pad))
    for text in ax.texts:
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if np.isfinite(box.extents).all():
                obstacles.append(box.padded(pad))
    # Protect the actual selected points as well as the other labels.
    radius = 2 * px_per_pt
    obstacles += [Bbox.from_extents(x-radius, y-radius, x+radius, y+radius)
                  for x, y in targets]
    point_px = ax.transData.transform(points) if len(points) else np.empty((0, 2))
    placed, output = [], []
    directions = ((1, 1), (-1, 1), (1, -1), (-1, -1),
                  (1, 0), (-1, 0), (0, 1), (0, -1))
    for (_, row), anchor in zip(selected.iterrows(), targets):
        probe = ax.text(0.5, 0.5, row["GENE_NAME"], transform=ax.transAxes,
                        fontsize=fontsize, fontstyle="italic", ha="center", va="center")
        measured = probe.get_window_extent(renderer)
        # Keep a margin for differences between screen/PDF text metrics.
        width, height = measured.width + 2*pad, measured.height + 2*pad
        probe.remove()
        if width >= area.width or height >= area.height:
            raise RuntimeError(f"Label {row['GENE_NAME']} does not fit; lower label_fontsize")
        candidates = []
        for gap in (4, 8, 12, 18, 26, 38, 55, 80):
            for dx, dy in directions:
                candidates.append((anchor[0] + dx*(width/2 + gap*px_per_pt),
                                   anchor[1] + dy*(height/2 + gap*px_per_pt)))
        # A grid fallback finds free space when top-ranked genes cluster together.
        for cy in np.linspace(area.y0+height/2, area.y1-height/2, 25):
            for cx in np.linspace(area.x0+width/2, area.x1-width/2, 25):
                candidates.append((cx, cy))
        best, best_score = None, np.inf
        for cx, cy in candidates:
            box = Bbox.from_bounds(cx-width/2, cy-height/2, width, height)
            if not _contains(area, box) or any(box.overlaps(b) for b in obstacles + placed):
                continue
            # Prefer short leaders and, softly, less densely populated locations.
            covered = ((point_px[:, 0] >= box.x0) & (point_px[:, 0] <= box.x1)
                       & (point_px[:, 1] >= box.y0) & (point_px[:, 1] <= box.y1)).sum()
            distance_pt = np.linalg.norm(np.array([cx, cy])-anchor) / px_per_pt
            score = distance_pt + 0.2 * covered
            if score < best_score:
                best, best_score = box, score
        if best is None:
            raise RuntimeError(
                f"Cannot place {row['GENE_NAME']} without overlap; "
                "try max_genes=4 or a smaller label_fontsize. No genes were silently dropped."
            )
        placed.append(best)
        output.append((row, anchor, best))
    return output


def label_neff_scatter(
    axes: Sequence[Axes],
    df: pd.DataFrame,
    gene_name: Callable[[str], Optional[str]],
    max_genes: int = 5,
    label_fontsize: float = 9.0,
    gap_pt: float = 1.5,
    ranking_mode: str = "all_three",
) -> pd.DataFrame:
    """Add the SAME <=5 gene labels to a finished three-panel scatter figure.

    Uses true renderer text bounds, not adjustText version-specific parameters.
    Raises instead of saving overlapping labels if the requested labels do not fit.
    Text positions are relative to each axes; never move data or change axis limits.
    Do not call tight_layout/change limits/resize axes after this function.
    Select once with ranking_mode (all_three by default, or first_two), then
    label the same genes in every panel, even on/below the third panel's diagonal.
    The returned table records all ratios, ranking_mode, and selection_score;
    min_ratio retains its meaning as the minimum of all three ratios.
    """
    axes = list(axes)
    if len(axes) != 3 or len({id(ax.figure) for ax in axes}) != 1:
        raise ValueError("Pass three axes from the same figure, in Original/C/CB comparison order")
    if not np.isfinite(label_fontsize) or label_fontsize <= 0 or gap_pt < 0:
        raise ValueError("label_fontsize must be positive and gap_pt nonnegative")
    if any(ax.get_xscale() != "linear" or ax.get_yscale() != "linear" for ax in axes):
        raise ValueError("This helper is intended for the existing LINEAR Neff scatter plots")
    fig = axes[0].figure
    fig.canvas.draw()  # Finalize axes geometry and renderer before measuring text.
    if not hasattr(fig.canvas, "get_renderer"):
        raise RuntimeError("Use an Agg-compatible canvas (e.g. MPLBACKEND=Agg on the server)")
    # Common visibility filter applies ONLY to candidate labels, never to scatter points.
    numeric = df.copy()
    for col in NEFF_COLS:
        numeric[col] = pd.to_numeric(numeric[col], errors="coerce")
    visible = np.ones(len(numeric), dtype=bool)
    for ax, (x_col, y_col) in zip(axes, PAIRS):
        xmin, xmax = sorted(ax.get_xlim())
        ymin, ymax = sorted(ax.get_ylim())
        visible &= (numeric[x_col].between(xmin, xmax)
                    & numeric[y_col].between(ymin, ymax)).fillna(False).to_numpy(bool)
    selected = select_shared_genes(
        numeric.loc[visible], gene_name, max_genes, ranking_mode=ranking_mode
    )
    if selected.empty:
        return selected
    plans = []
    # Compute all three plans before adding any final labels (no per-panel dropping).
    for ax, (x_col, y_col) in zip(axes, PAIRS):
        points = numeric[[x_col, y_col]].to_numpy(dtype=float, na_value=np.nan)
        points = points[np.isfinite(points).all(axis=1)]
        plans.append(_positions(ax, selected, x_col, y_col, points, label_fontsize, gap_pt))
    all_artists, label_groups = [], []
    try:
        for ax, (x_col, y_col), plan in zip(axes, PAIRS, plans):
            texts = []
            for row, anchor, box in plan:
                center = np.array([(box.x0+box.x1)/2, (box.y0+box.y1)/2])
                tx, ty = ax.transAxes.inverted().transform(center)
                # Start the leader at the padded box boundary, NOT through the letters.
                end = np.clip(anchor, [box.x0, box.y0], [box.x1, box.y1])
                ex, ey = ax.transAxes.inverted().transform(end)
                arrow = ax.annotate(
                    "", xy=(row[x_col], row[y_col]), xycoords="data",
                    xytext=(ex, ey), textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="-", color="0.4", lw=0.5,
                                    shrinkA=0, shrinkB=1.5), zorder=4,
                )
                text = ax.text(tx, ty, row["GENE_NAME"], transform=ax.transAxes,
                               fontsize=label_fontsize, fontstyle="italic",
                               ha="center", va="center", color="black", zorder=5)
                text.set_gid("neff-gene-label")
                for artist in (arrow, text):
                    artist.set_in_layout(False)
                    all_artists.append(artist)
                texts.append(text)
            label_groups.append(texts)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for ax, texts in zip(axes, label_groups):
            boxes = [text.get_window_extent(renderer) for text in texts]
            area = ax.get_window_extent(renderer)
            legend = ax.get_legend()
            legend_box = (legend.get_window_extent(renderer)
                          if legend is not None and legend.get_visible() else None)
            for i, box in enumerate(boxes):
                if (not _contains(area, box)
                        or any(box.overlaps(other) for other in boxes[:i])
                        or (legend_box is not None and box.overlaps(legend_box))):
                    raise RuntimeError("Final label overlap/boundary check failed; decrease label_fontsize")
    except Exception:
        for artist in reversed(all_artists):
            artist.remove()
        raise
    return selected
