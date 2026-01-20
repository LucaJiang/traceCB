# Plot correlation density for all studies combined
from visual.utils import *
import warnings

warnings.filterwarnings("ignore")

# target_qtdids to include
target_qtdids = [
    "QTD000021",
    "QTD000031",
    "QTD000066",
    "QTD000067",
    "QTD000069",
    "QTD000073",
    "QTD000081",
    "QTD000115",
    "QTD000371",
    "QTD000372",
]


def plot_cor_density_all_studies(summary_df_all, summary_sign_df_all):
    """
    Plot correlation density for all studies combined.
    Applies the same filtering logic as individual studies.
    """
    # Apply filtering: only include genes with significant H1 and H2 heritability
    z_threshold = p2z(0.1)
    target_summary_df = summary_df_all.copy()

    # Filter for significant heritability
    target_summary_df.loc[:, "H1_SIGN"] = (
        target_summary_df.H1SQ / target_summary_df.H1SQSE
    ) > z_threshold
    target_summary_df.loc[:, "H2_SIGN"] = (
        target_summary_df.H2SQ / target_summary_df.H2SQSE
    ) > z_threshold
    both_sign_index = (
        target_summary_df.H1_SIGN
        & target_summary_df.H2_SIGN
        & (target_summary_df.H1SQ > 1e-12)
        & (target_summary_df.H2SQ > 1e-12)
    )
    target_summary_df = target_summary_df.loc[both_sign_index, :]

    # Use all correlation significant data (no additional filtering needed)
    target_summary_sign_df = summary_sign_df_all.copy()

    # 添加分组标签
    target_summary_df.loc[:, "Type"] = "Heritability Significant"
    target_summary_sign_df.loc[:, "Type"] = "Correlation Significant"

    # 合并数据
    combined_df = pd.concat(
        [target_summary_df, target_summary_sign_df], ignore_index=True
    )

    print(f"Heritability Significant genes: {len(target_summary_df)}")
    print(f"Correlation Significant genes: {len(target_summary_sign_df)}")
    print(f"Combined data points: {len(combined_df)}")

    # 使用 displot 并添加 rug
    g = sns.displot(
        data=combined_df,
        x="COR_X_ORI",
        hue="Type",
        kind="hist",
        bins=200,
        stat="count",
        kde=True,
        alpha=0.5,
        palette={
            "Heritability Significant": "#120eff",
            "Correlation Significant": "#db1717f8",
        },
        hue_order=["Heritability Significant", "Correlation Significant"],
        height=4,
        aspect=2.0,
    )

    # 设置图例
    sns.move_legend(
        g,
        "upper left",
        bbox_to_anchor=(0.16, 0.85),
        title="Significance Type",
        frameon=True,
        fontsize=8,
        title_fontsize=8,
    )

    g.set_axis_labels("Genetic Correlation", "Count")
    g.fig.suptitle("All Studies Combined")
    g.despine(top=False, right=False)

    plt.tight_layout()
    plt.savefig(
        os.path.join(save_path, "single_study", "f3cor_density_all_studies.pdf"),
        bbox_inches="tight",
    )
    plt.close()
    print(f"Saved: {save_path}/single_study/f3cor_density_all_studies.pdf")


if __name__ == "__main__":
    # Create output directory if it doesn't exist
    if not os.path.exists(f"{save_path}/single_study/"):
        os.makedirs(f"{save_path}/single_study/")

    # Load all summary data
    print("Loading summary data...")
    summary_sign_df_all, summary_df_all = load_all_summary()

    # Filter to only include target studies
    summary_df_all = summary_df_all[summary_df_all.QTDid.isin(target_qtdids)]
    summary_sign_df_all = summary_sign_df_all[
        summary_sign_df_all.QTDid.isin(target_qtdids)
    ]

    print(f"Total studies included: {summary_df_all.QTDid.nunique()}")
    print(f"Total genes (heritability): {len(summary_df_all)}")
    print(f"Total genes (correlation significant): {len(summary_sign_df_all)}")

    # Plot correlation density for all studies combined
    plot_cor_density_all_studies(summary_df_all, summary_sign_df_all)
