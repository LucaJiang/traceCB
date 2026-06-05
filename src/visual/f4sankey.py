import argparse
import glob
import json
import os
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from upsetplot import UpSet, from_memberships


REPO_ROOT = Path(__file__).resolve().parents[2]
LOCAL_DATA_DIR = REPO_ROOT / "data" / "coloc"
LOCAL_SAVE_DIR = REPO_ROOT / "data" / "img" / "eas_eqtlgen"
SERVER_DATA_DIR = Path("/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc")
SERVER_SAVE_DIR = Path(
    "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/coloc"
)
METADATA_PATH = Path(__file__).with_name("metadata.json")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Figure 4 colocalization Sankey and study-level gene UpSet."
    )
    parser.add_argument("--prefix", default="bcx", help="Coloc file prefix.")
    parser.add_argument(
        "--prefixes",
        nargs="+",
        default=["bcx", "bbj"],
        help="Prefixes to compare in the combined UpSet plot.",
    )
    parser.add_argument(
        "--data-dir",
        default=None,
        help="Directory containing *_coloc.csv files. Defaults to local data/coloc if present.",
    )
    parser.add_argument(
        "--save-dir",
        default=None,
        help="Output directory. Defaults to local data/img/eas_eqtlgen if present.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.7,
        help="Posterior probability threshold for H3/H4 calls.",
    )
    parser.add_argument(
        "--upset-method",
        default="traceCB",
        choices=["original", "traceC", "traceCB", "any"],
        help="Method used to define colocalized genes in the UpSet plot.",
    )
    parser.add_argument(
        "--min-subset-size",
        type=int,
        default=1,
        help="Minimum UpSet intersection size to display.",
    )
    parser.add_argument(
        "--max-intersections",
        type=int,
        default=20,
        help="Maximum number of study-intersection patterns per method subplot.",
    )
    return parser.parse_args()


def resolve_default_path(user_path, local_path, server_path):
    if user_path is not None:
        return Path(user_path)
    if local_path.exists():
        return local_path
    return server_path


def load_metadata():
    with open(METADATA_PATH, "r") as f:
        return json.load(f)


def parse_qtdid(file_path):
    match = re.search(r"(QTD\d+)", os.path.basename(file_path))
    if match is None:
        raise ValueError(f"Could not parse QTDid from file name: {file_path}")
    return match.group(1)


def parse_trait(file_path):
    filename = os.path.basename(file_path)
    return filename.split("_eQTLGen_")[0].removesuffix("_coloc")


def load_coloc_data(data_dir, prefix):
    coloc_files = sorted(glob.glob(str(data_dir / f"{prefix}*_coloc.csv")))
    print(f"Found {len(coloc_files)} coloc files in {data_dir}")
    if not coloc_files:
        raise FileNotFoundError(
            f"No coloc files matched {data_dir}/{prefix}*_coloc.csv"
        )

    df_list = []
    for file_path in coloc_files:
        try:
            tmp_df = pd.read_csv(file_path)
        except pd.errors.EmptyDataError:
            continue
        if tmp_df.empty:
            continue
        tmp_df["QTDid"] = parse_qtdid(file_path)
        tmp_df["trait"] = parse_trait(file_path)
        tmp_df["source_file"] = os.path.basename(file_path)
        df_list.append(tmp_df)

    if not df_list:
        raise ValueError(f"All matched coloc files were empty: {data_dir}")
    df = pd.concat(df_list, ignore_index=True)
    print(f"Total rows before filtering: {len(df)}")
    return df


def has_any_h3_h4_signal(df, threshold):
    return (
        (df["pp_h3_original"] > threshold)
        | (df["p_original"] > threshold)
        | (df["pp_h3_traceC"] > threshold)
        | (df["p_traceC"] > threshold)
        | (df["pp_h3_traceCB"] > threshold)
        | (df["p_traceCB"] > threshold)
    )


def get_category(df, h3_col, h4_col, threshold):
    """Return category for each row: 0=COLOC(H4), 1=Independent(H3), 2=Other."""
    is_h3 = df[h3_col] > threshold
    is_h4 = df[h4_col] > threshold
    category = pd.Series(2, index=df.index)
    category[is_h3] = 1
    category[is_h4] = 0
    return category


def color_to_rgba(color, alpha=0.5):
    if color.startswith("rgb("):
        rgb_values = color[4:-1].split(",")
        r, g, b = [int(v.strip()) for v in rgb_values]
    else:
        hex_color = color.lstrip("#")
        r, g, b = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r}, {g}, {b}, {alpha})"


def plot_sankey(df, save_dir, prefix, threshold):
    df = df.loc[has_any_h3_h4_signal(df, threshold)].copy()
    print(f"Total rows after filtering (at least one h3/h4 > {threshold}): {len(df)}")

    h3_original = (df["pp_h3_original"] > threshold).sum()
    h3_traceC = (df["pp_h3_traceC"] > threshold).sum()
    h3_traceCB = (df["pp_h3_traceCB"] > threshold).sum()
    h4_original = (df["p_original"] > threshold).sum()
    h4_traceC = (df["p_traceC"] > threshold).sum()
    h4_traceCB = (df["p_traceCB"] > threshold).sum()
    other_original = (
        (df["pp_h3_original"] <= threshold) & (df["p_original"] <= threshold)
    ).sum()
    other_traceC = (
        (df["pp_h3_traceC"] <= threshold) & (df["p_traceC"] <= threshold)
    ).sum()
    other_traceCB = (
        (df["pp_h3_traceCB"] <= threshold) & (df["p_traceCB"] <= threshold)
    ).sum()

    print("=== Statistics ===")
    print(f"Original: H3={h3_original}, H4={h4_original}, Other={other_original}")
    print(f"traceC:   H3={h3_traceC}, H4={h4_traceC}, Other={other_traceC}")
    print(f"traceCB:  H3={h3_traceCB}, H4={h4_traceCB}, Other={other_traceCB}")

    cat_original = get_category(df, "pp_h3_original", "p_original", threshold)
    cat_traceC = get_category(df, "pp_h3_traceC", "p_traceC", threshold)
    cat_traceCB = get_category(df, "pp_h3_traceCB", "p_traceCB", threshold)

    flow_orig_to_traceC = pd.crosstab(cat_original, cat_traceC)
    flow_traceC_to_traceCB = pd.crosstab(cat_traceC, cat_traceCB)
    print("\n=== Flow: Original -> traceC ===")
    print(flow_orig_to_traceC)
    print("\n=== Flow: traceC -> traceCB ===")
    print(flow_traceC_to_traceCB)

    orig_counts = cat_original.value_counts().sort_index()
    traceC_counts = cat_traceC.value_counts().sort_index()
    traceCB_counts = cat_traceCB.value_counts().sort_index()
    nodes = [
        f"COLOC ({orig_counts.get(0, 0)})",
        f"Independent ({orig_counts.get(1, 0)})",
        f"Undetermined ({orig_counts.get(2, 0)})",
        f"COLOC ({traceC_counts.get(0, 0)})",
        f"Independent ({traceC_counts.get(1, 0)})",
        f"Undetermined ({traceC_counts.get(2, 0)})",
        f"COLOC ({traceCB_counts.get(0, 0)})",
        f"Independent ({traceCB_counts.get(1, 0)})",
        f"Undetermined ({traceCB_counts.get(2, 0)})",
    ]

    coloc_color = "rgb(124, 200, 124)"
    independent_color = "rgb(150, 181, 248)"
    undetermined_color = "rgb(235, 160, 152)"
    category_colors = [coloc_color, independent_color, undetermined_color]

    sources = []
    targets = []
    values = []
    link_colors = []
    for src_cat in range(3):
        for tgt_cat in range(3):
            if (
                src_cat in flow_orig_to_traceC.index
                and tgt_cat in flow_orig_to_traceC.columns
            ):
                val = flow_orig_to_traceC.loc[src_cat, tgt_cat]
                if val > 0:
                    sources.append(src_cat)
                    targets.append(tgt_cat + 3)
                    values.append(val)
                    link_colors.append(color_to_rgba(category_colors[src_cat]))

    for src_cat in range(3):
        for tgt_cat in range(3):
            if (
                src_cat in flow_traceC_to_traceCB.index
                and tgt_cat in flow_traceC_to_traceCB.columns
            ):
                val = flow_traceC_to_traceCB.loc[src_cat, tgt_cat]
                if val > 0:
                    sources.append(src_cat + 3)
                    targets.append(tgt_cat + 6)
                    values.append(val)
                    link_colors.append(color_to_rgba(category_colors[src_cat]))

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=20,
                    thickness=25,
                    line=dict(color="white", width=1),
                    label=nodes,
                    color=[
                        coloc_color,
                        independent_color,
                        undetermined_color,
                        coloc_color,
                        independent_color,
                        undetermined_color,
                        coloc_color,
                        independent_color,
                        undetermined_color,
                    ],
                    x=[0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.99, 0.99, 0.99],
                    y=[0.15, 0.5, 0.85, 0.15, 0.5, 0.85, 0.15, 0.5, 0.85],
                ),
                link=dict(source=sources, target=targets, value=values, color=link_colors),
            )
        ]
    )

    yloc = 1.04
    fig.update_layout(
        title={
            "text": "Colocalization Pattern Flow",
            "x": 0.5,
            "xanchor": "center",
            "font": {"size": 16, "family": "Arial", "color": "#333"},
        },
        font=dict(size=13, family="Arial", color="#333"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        width=550,
        height=380,
        margin=dict(l=5, r=5, t=70, b=10),
        annotations=[
            dict(
                text="Original",
                x=-0.01,
                y=yloc,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14, color="#555", family="Arial"),
            ),
            dict(
                text="traceC",
                x=0.5,
                y=yloc,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14, color="#555", family="Arial"),
            ),
            dict(
                text="traceCB",
                x=1.01,
                y=yloc,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14, color="#555", family="Arial"),
            ),
        ],
    )

    html_path = save_dir / f"{prefix}_sankey_by_category.html"
    pdf_path = save_dir / f"{prefix}_sankey_by_category.pdf"
    fig.write_html(html_path)
    print(f"Sankey diagram saved to {html_path}")
    try:
        fig.write_image(pdf_path, width=550, height=380, scale=2)
        print(f"PDF saved to {pdf_path}")
    except Exception as e:
        print(f"Could not save Sankey PDF: {e}")


def get_coloc_gene_sets_by_study(df, meta_data, method, threshold):
    if method == "any":
        coloc_mask = (
            (df["p_original"] > threshold)
            | (df["p_traceC"] > threshold)
            | (df["p_traceCB"] > threshold)
        )
    else:
        coloc_mask = df[f"p_{method}"] > threshold

    coloc_df = df.loc[coloc_mask, ["QTDid", "gene"]].dropna().copy()
    coloc_df["gene"] = coloc_df["gene"].astype(str)

    gene_sets = {qtdid: set() for qtdid in meta_data["QTDids"]}
    for qtdid, group_df in coloc_df.groupby("QTDid"):
        gene_sets.setdefault(qtdid, set()).update(group_df["gene"])
    return gene_sets


def write_upset_tables(gene_sets, meta_data, save_dir, prefix, method):
    all_genes = sorted(set().union(*gene_sets.values())) if gene_sets else []
    membership_rows = []
    pattern_counter = Counter()
    pattern_genes = {}
    for gene in all_genes:
        member_qtdids = tuple(
            qtdid for qtdid in meta_data["QTDids"] if gene in gene_sets[qtdid]
        )
        pattern_counter[member_qtdids] += 1
        pattern_genes.setdefault(member_qtdids, []).append(gene)
        row = {"gene": gene, "n_studies": len(member_qtdids)}
        for qtdid in meta_data["QTDids"]:
            row[meta_data["id2name"][qtdid]] = gene in gene_sets[qtdid]
        membership_rows.append(row)

    membership_df = pd.DataFrame(membership_rows)
    membership_path = save_dir / f"{prefix}_{method}_coloc_gene_by_study.csv"
    membership_df.to_csv(membership_path, index=False)

    intersection_rows = []
    for member_qtdids, count in pattern_counter.most_common():
        study_names = [meta_data["id2name"][qtdid] for qtdid in member_qtdids]
        intersection_rows.append(
            {
                "studies": ";".join(study_names),
                "qtdids": ";".join(member_qtdids),
                "n_studies": len(member_qtdids),
                "n_genes": count,
                "genes": ";".join(pattern_genes[member_qtdids]),
            }
        )
    intersection_df = pd.DataFrame(intersection_rows)
    intersection_path = save_dir / f"{prefix}_{method}_coloc_gene_intersections.csv"
    intersection_df.to_csv(intersection_path, index=False)
    print(f"Gene-by-study table saved to {membership_path}")
    print(f"Gene intersection table saved to {intersection_path}")


def plot_study_gene_upset(
    df, meta_data, save_dir, prefix, method, threshold, min_subset_size
):
    gene_sets = get_coloc_gene_sets_by_study(df, meta_data, method, threshold)
    study_labels = {
        qtdid: meta_data["id2name"].get(qtdid, qtdid) for qtdid in meta_data["QTDids"]
    }
    print(
        f"\n=== Unique H4-colocalized genes for UpSet "
        f"({method}, p > {threshold}) ==="
    )
    for qtdid in meta_data["QTDids"]:
        print(f"{study_labels[qtdid]}: {len(gene_sets[qtdid])}")

    all_genes = sorted(set().union(*gene_sets.values()))
    if not all_genes:
        print("No colocalized genes found for UpSet; skip plotting.")
        return

    memberships = [
        [study_labels[qtdid] for qtdid in meta_data["QTDids"] if gene in gene_sets[qtdid]]
        for gene in all_genes
    ]
    upset_data = from_memberships(memberships)

    plt.figure(figsize=(12, 6))
    upset = UpSet(
        upset_data,
        subset_size="count",
        sort_by="cardinality",
        sort_categories_by="-input",
        min_subset_size=min_subset_size,
        show_counts=True,
    )
    upset.plot()
    plt.suptitle(
        f"Unique colocalized genes shared across studies ({method}, PP.H4 > {threshold})",
        fontsize=12,
    )
    output_path = save_dir / f"{prefix}_{method}_coloc_gene_upset.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"UpSet plot saved to {output_path}")

    write_upset_tables(gene_sets, meta_data, save_dir, prefix, method)


def count_intersections(gene_sets, meta_data):
    all_genes = sorted(set().union(*gene_sets.values())) if gene_sets else []
    pattern_counter = Counter()
    for gene in all_genes:
        pattern = tuple(
            qtdid for qtdid in meta_data["QTDids"] if gene in gene_sets[qtdid]
        )
        if pattern:
            pattern_counter[pattern] += 1
    return pattern_counter


def select_patterns(prefix_counters, min_subset_size, max_intersections):
    all_patterns = set()
    for counter in prefix_counters.values():
        all_patterns.update(counter.keys())

    eligible_patterns = [
        pattern
        for pattern in all_patterns
        if max(prefix_counters[prefix].get(pattern, 0) for prefix in prefix_counters)
        >= min_subset_size
    ]
    eligible_patterns.sort(
        key=lambda pattern: (
            -sum(prefix_counters[prefix].get(pattern, 0) for prefix in prefix_counters),
            -max(prefix_counters[prefix].get(pattern, 0) for prefix in prefix_counters),
            -len(pattern),
            pattern,
        )
    )
    return eligible_patterns[:max_intersections]


def write_combined_upset_tables(
    method_gene_sets, meta_data, save_dir, prefixes, methods, threshold
):
    set_rows = []
    intersection_rows = []
    for method in methods:
        for prefix in prefixes:
            gene_sets = method_gene_sets[method][prefix]
            for qtdid in meta_data["QTDids"]:
                set_rows.append(
                    {
                        "method": method,
                        "prefix": prefix,
                        "qtdid": qtdid,
                        "study": meta_data["id2name"][qtdid],
                        "n_genes": len(gene_sets[qtdid]),
                    }
                )

            counter = count_intersections(gene_sets, meta_data)
            for pattern, count in counter.most_common():
                intersection_rows.append(
                    {
                        "method": method,
                        "prefix": prefix,
                        "qtdids": ";".join(pattern),
                        "studies": ";".join(meta_data["id2name"][q] for q in pattern),
                        "n_studies": len(pattern),
                        "n_genes": count,
                    }
                )

    set_path = save_dir / f"bcx_bbj_coloc_gene_set_sizes_pph4_gt_{threshold}.csv"
    intersection_path = (
        save_dir / f"bcx_bbj_coloc_gene_intersections_pph4_gt_{threshold}.csv"
    )
    pd.DataFrame(set_rows).to_csv(set_path, index=False)
    pd.DataFrame(intersection_rows).to_csv(intersection_path, index=False)
    print(f"Combined UpSet set sizes saved to {set_path}")
    print(f"Combined UpSet intersections saved to {intersection_path}")


def draw_upset_panel(
    subfig,
    method,
    prefix_gene_sets,
    meta_data,
    prefixes,
    colors,
    min_subset_size,
    max_intersections,
):
    prefix_counters = {
        prefix: count_intersections(prefix_gene_sets[prefix], meta_data)
        for prefix in prefixes
    }
    patterns = select_patterns(prefix_counters, min_subset_size, max_intersections)
    study_labels = [meta_data["id2name"][qtdid] for qtdid in meta_data["QTDids"]]
    grid = subfig.add_gridspec(
        2,
        2,
        width_ratios=[1.1, 5.9],
        height_ratios=[2.3, 1.45],
        wspace=0.0,
        hspace=0.03,
    )
    set_ax = subfig.add_subplot(grid[1, 0])
    bar_ax = subfig.add_subplot(grid[0, 1])
    matrix_ax = subfig.add_subplot(grid[1, 1], sharex=bar_ax)
    subfig.add_subplot(grid[0, 0]).axis("off")

    y_positions = list(range(len(meta_data["QTDids"])))
    bar_height = 0.34
    for i, prefix in enumerate(prefixes):
        offset = (i - (len(prefixes) - 1) / 2) * bar_height
        sizes = [len(prefix_gene_sets[prefix][qtdid]) for qtdid in meta_data["QTDids"]]
        set_ax.barh(
            [y + offset for y in y_positions],
            sizes,
            height=bar_height,
            color=colors[prefix],
            label=prefix,
            alpha=0.9,
        )

    set_ax.set_yticks(y_positions)
    set_ax.set_yticklabels(study_labels, fontsize=8)
    set_ax.invert_yaxis()
    max_set_size = max(
        len(prefix_gene_sets[prefix][qtdid])
        for prefix in prefixes
        for qtdid in meta_data["QTDids"]
    )
    set_ax.set_xlim(max_set_size * 1.08, 0)
    set_ax.set_xlabel("Set size", fontsize=8)
    set_ax.tick_params(axis="x", labelsize=7)
    for spine in ["top", "right"]:
        set_ax.spines[spine].set_visible(False)

    x_positions = list(range(len(patterns)))
    group_width = 0.72
    single_width = group_width / max(len(prefixes), 1)
    for i, prefix in enumerate(prefixes):
        offset = (i - (len(prefixes) - 1) / 2) * single_width
        counts = [prefix_counters[prefix].get(pattern, 0) for pattern in patterns]
        bars = bar_ax.bar(
            [x + offset for x in x_positions],
            counts,
            width=single_width * 0.9,
            color=colors[prefix],
            label=prefix,
            alpha=0.92,
        )
        bar_ax.bar_label(
            bars, labels=[str(c) if c > 0 else "" for c in counts], fontsize=6
        )

    bar_ax.set_title(method, fontsize=12, pad=8)
    bar_ax.set_ylabel("Intersection size", fontsize=9)
    bar_ax.legend(frameon=False, fontsize=8, ncol=len(prefixes), loc="upper right")
    bar_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    bar_ax.tick_params(axis="y", labelsize=8)
    for spine in ["top", "right"]:
        bar_ax.spines[spine].set_visible(False)

    matrix_ax.set_ylim(-0.5, len(meta_data["QTDids"]) - 0.5)
    matrix_ax.invert_yaxis()
    matrix_ax.set_yticks(y_positions)
    matrix_ax.set_yticklabels([])
    matrix_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    matrix_ax.tick_params(axis="y", length=0)
    for x, pattern in zip(x_positions, patterns):
        matrix_ax.scatter(
            [x] * len(y_positions),
            y_positions,
            s=13,
            color="#d5d5d5",
            zorder=1,
        )
        member_positions = [
            idx for idx, qtdid in enumerate(meta_data["QTDids"]) if qtdid in pattern
        ]
        if member_positions:
            matrix_ax.plot(
                [x] * len(member_positions),
                member_positions,
                color="#303030",
                linewidth=1.1,
                zorder=2,
            )
            matrix_ax.scatter(
                [x] * len(member_positions),
                member_positions,
                s=18,
                color="#303030",
                zorder=3,
            )
    matrix_ax.set_xlim(-0.7, max(len(patterns) - 0.3, 0.3))
    matrix_ax.set_xlabel("Study intersection pattern", fontsize=8)
    for spine in ["top", "right"]:
        matrix_ax.spines[spine].set_visible(False)


def plot_combined_prefix_method_upset(
    prefix_dfs,
    meta_data,
    save_dir,
    prefixes,
    threshold,
    min_subset_size,
    max_intersections,
):
    methods = ["original", "traceC", "traceCB"]
    colors = {"bcx": "#4C78A8", "bbj": "#F58518"}
    for prefix in prefixes:
        colors.setdefault(prefix, f"C{len(colors)}")

    method_gene_sets = {
        method: {
            prefix: get_coloc_gene_sets_by_study(
                prefix_dfs[prefix], meta_data, method, threshold
            )
            for prefix in prefixes
        }
        for method in methods
    }

    fig = plt.figure(figsize=(13, 14), constrained_layout=True)
    subfigs = fig.subfigures(len(methods), 1, hspace=0.06)
    if len(methods) == 1:
        subfigs = [subfigs]
    for subfig, method in zip(subfigs, methods):
        draw_upset_panel(
            subfig,
            method,
            method_gene_sets[method],
            meta_data,
            prefixes,
            colors,
            min_subset_size,
            max_intersections,
        )

    fig.suptitle(
        f"Unique colocalized genes shared across studies (PP.H4 > {threshold})",
        fontsize=14,
    )
    output_path = save_dir / "bcx_bbj_methods_coloc_gene_upset.pdf"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Combined bcx/bbj method UpSet plot saved to {output_path}")

    write_combined_upset_tables(
        method_gene_sets, meta_data, save_dir, prefixes, methods, threshold
    )


def main():
    args = parse_args()
    data_dir = resolve_default_path(args.data_dir, LOCAL_DATA_DIR, SERVER_DATA_DIR)
    save_dir = resolve_default_path(args.save_dir, LOCAL_SAVE_DIR, SERVER_SAVE_DIR)
    save_dir.mkdir(parents=True, exist_ok=True)

    meta_data = load_metadata()
    df = load_coloc_data(data_dir, args.prefix)
    plot_sankey(df, save_dir, args.prefix, args.threshold)
    plot_study_gene_upset(
        df,
        meta_data,
        save_dir,
        args.prefix,
        args.upset_method,
        args.threshold,
        args.min_subset_size,
    )

    prefix_dfs = {}
    for prefix in args.prefixes:
        prefix_dfs[prefix] = load_coloc_data(data_dir, prefix)
    plot_combined_prefix_method_upset(
        prefix_dfs,
        meta_data,
        save_dir,
        args.prefixes,
        args.threshold,
        args.min_subset_size,
        args.max_intersections,
    )


if __name__ == "__main__":
    main()
