# plot f5 egene upset for afr+gtex
from utils import *

# import upsetplot
from upsetplot import UpSet, plot

target_qtdids = ["QTD000067", "QTD000371", "QTD000069", "QTD000081"]
study_path_main = "/gpfs1/scratch/wjiang49/xpmm/AFR_GTEx"
save_path = "/gpfs1/scratch/wjiang49/xpmm/AFR_GTEx/results"


def format_eGene_df(df):
    gene_indicators = pd.DataFrame(
        {
            meta_data["method_name"][0]: df["TAR_SeSNP"] > 0,
            meta_data["method_name"][1]: df["TAR_CeSNP"] > 0,
            meta_data["method_name"][2]: df["TAR_TeSNP"] > 0,
        }
    )
    result = gene_indicators.value_counts().sort_index()
    # remove rows with all False
    result = result.drop((False, False, False))
    return result.to_frame(name="value")


def plot_upset(df, title, save_file):
    plt.figure(figsize=(6, 6))
    upset = UpSet(
        df["value"],
        sort_categories_by="-input",
        show_counts=True,
    )
    for i in range(3):
        upset.style_categories(
            meta_data["method_name"][i],
            bar_facecolor=meta_data["Colors"][meta_data["method_name"][i]],
        )
    upset.plot()
    plt.suptitle(title)
    plt.savefig(save_file, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    summary_sign_df_all, _ = load_all_summary(study_path_main)
    for target_qtdid in target_qtdids:
        df = summary_sign_df_all[summary_sign_df_all.QTDid == target_qtdid].copy()
        eGene_df = format_eGene_df(df)
        title = "eGene for " + meta_data["id2name"][target_qtdid]
        save_file = f"{save_path}/f6egene_upset_{target_qtdid}.png"
        plot_upset(eGene_df, title, save_file)
