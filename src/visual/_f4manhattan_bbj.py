from utils import *

# Input: gene_id, chr_num, gene_name, pos_range
plot_gene_info = {
    ("ENSG00000118369", 11, "USP35", "77,899,858-77,925,757"), # gene, chr, gene_name, pos_range
}
# Global variables
gwas_path = "/home/wjiang49/group/wjiang49/data/xpmm/coloc/bcx/bcx_mon_GWAS.csv"
MIN_PVAL = 1e-20
EXTEND_RANGE = 100_000

def load_qtd_gene_data(qtdid, gene_id, chr_num):
    """
    Load gene data from specific QTD study
    """
    celltype = meta_data["id2celltype"][qtdid]
    file_path = f"{study_path_main}/{qtdid}/AUX_{celltype}/chr{chr_num}.csv"
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    # Load data
    df = pd.read_csv(file_path)
    
    # Filter for specific gene
    gene_data = df[df['GENE'] == gene_id].copy()
    
    if gene_data.empty:
        print(f"Gene {gene_id} not found in {qtdid}")
        return None
    
    # Add study info
    gene_data.loc[:, 'QTDid'] = qtdid
    gene_data.loc[:, 'Study_Name'] = meta_data["id2name"][qtdid]
    gene_data.loc[:, 'Celltype'] = celltype
    
    return gene_data

def create_combined_manhattan_plot(gene_id, chr_num, gene_name, pos_range):
    """
    Create combined Manhattan plot with GWAS on top and eQTL studies below
    """
    # Parse position range and extend by EXTEND_RANGE on both sides
    start_str, end_str = pos_range.split('-')
    gene_start = int(start_str.replace(',', ''))
    gene_end = int(end_str.replace(',', ''))
    
    # Extend range
    extended_start = gene_start - EXTEND_RANGE
    extended_end = gene_end + EXTEND_RANGE
    
    print(f"Gene range: {gene_start:,} - {gene_end:,}")
    print(f"Extended range: {extended_start:,} - {extended_end:,}")
    
    # Load GWAS data
    gwas_df = pd.read_csv(gwas_path)
    gwas_df = gwas_df[gwas_df["CHR"] == chr_num]
    gwas_df.loc[:, "PVAL_GWAS"] = gwas_df.Z.apply(z2p)
    gwas_df["PVAL_GWAS"] = gwas_df["PVAL_GWAS"].clip(lower=MIN_PVAL)
    
    # Filter GWAS data to extended range
    gwas_df = gwas_df[(gwas_df["POS"] >= extended_start) & (gwas_df["POS"] <= extended_end)]
    print(f"GWAS SNPs in extended range: {len(gwas_df)}")
    
    # Collect eQTL data from all QTD studies
    all_qtd_data = []
    for qtdid in meta_data["QTDids"]:
        gene_data = load_qtd_gene_data(qtdid, gene_id, chr_num)
        if gene_data is not None:
            # Merge with GWAS to get positions
            merged_data = pd.merge(gene_data, gwas_df[['SNP', 'POS']], 
                                 left_on='RSID', right_on='SNP', how='inner')
            if not merged_data.empty:
                # Filter to extended range
                merged_data = merged_data[(merged_data["POS"] >= extended_start) & 
                                        (merged_data["POS"] <= extended_end)]
                if not merged_data.empty:
                    merged_data['PVAL'] = merged_data['PVAL'].clip(lower=MIN_PVAL)
                    all_qtd_data.append(merged_data)
                    print(f"Loaded {qtdid} ({meta_data['id2name'][qtdid]}): {len(merged_data)} SNPs")
    
    if not all_qtd_data:
        print(f"No eQTL data found for gene {gene_id} in any QTD study")
        return None
    
    # Sort data by position
    gwas_df = gwas_df.sort_values(by="POS")
    for i in range(len(all_qtd_data)):
        all_qtd_data[i] = all_qtd_data[i].sort_values(by="POS")
    
    # Determine number of subplots (1 for GWAS + number of QTD studies with data)
    n_studies = len(all_qtd_data)
    n_plots = n_studies + 1  # GWAS + eQTL studies
    
    # Create figure with subplots
    fig, axes = plt.subplots(n_plots, 1, figsize=(12, 2*n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]
    elif n_plots == 2:
        axes = list(axes)
    
    # Plot 1: GWAS Manhattan plot
    if not gwas_df.empty:
        gwas_log10p = -np.log10(gwas_df["PVAL_GWAS"])
        gwas_pos = gwas_df["POS"] / 1e6  # Convert to Mb
        
        # Create significance level masks
        mask_8 = gwas_log10p > 8  # p < 1e-8
        mask_5 = (gwas_log10p > 5) & (gwas_log10p <= 8)  # 1e-8 <= p < 1e-5
        mask_ns = gwas_log10p <= 5  # p >= 1e-5
        
        # Plot points with different colors for significance levels
        axes[0].scatter(gwas_pos[mask_ns], gwas_log10p[mask_ns], c="gray", s=3, alpha=0.6)
        axes[0].scatter(gwas_pos[mask_5], gwas_log10p[mask_5], c="blue", s=3)
        axes[0].scatter(gwas_pos[mask_8], gwas_log10p[mask_8], c="red", s=3)
        
        # Add significance threshold lines
        max_gwas_log10p = max(gwas_log10p) if len(gwas_log10p) > 0 else 0
        if max_gwas_log10p > 8:
            axes[0].axhline(y=8, color="red", linestyle="--", alpha=0.5, linewidth=1)
        if max_gwas_log10p > 5:
            axes[0].axhline(y=5, color="blue", linestyle="--", alpha=0.5, linewidth=1)
    
    axes[0].set_ylabel("-log₁₀(P)", fontsize=10)
    axes[0].set_title(f"GWAS - {gene_name}", fontsize=10, fontweight="bold")
    axes[0].grid(True, alpha=0.3)
    axes[0].spines["right"].set_visible(False)
    axes[0].spines["top"].set_visible(False)
    
    # Plot eQTL Manhattan plots for each study
    for i, qtd_data in enumerate(all_qtd_data):
        ax_idx = i + 1  # Skip GWAS plot
        
        if not qtd_data.empty:
            eqtl_log10p = -np.log10(qtd_data["PVAL"])
            eqtl_pos = qtd_data["POS"] / 1e6  # Convert to Mb
            
            # Create significance level masks
            mask_8 = eqtl_log10p > 8  # p < 1e-8
            mask_5 = (eqtl_log10p > 5) & (eqtl_log10p <= 8)  # 1e-8 <= p < 1e-5
            mask_ns = eqtl_log10p <= 5  # p >= 1e-5
            
            # Plot points with different colors for significance levels
            axes[ax_idx].scatter(eqtl_pos[mask_ns], eqtl_log10p[mask_ns], c="gray", s=3, alpha=0.6)
            axes[ax_idx].scatter(eqtl_pos[mask_5], eqtl_log10p[mask_5], c="blue", s=3)
            axes[ax_idx].scatter(eqtl_pos[mask_8], eqtl_log10p[mask_8], c="red", s=3)
            
            # Add significance threshold lines
            max_eqtl_log10p = max(eqtl_log10p) if len(eqtl_log10p) > 0 else 0
            if max_eqtl_log10p > 8:
                axes[ax_idx].axhline(y=8, color="red", linestyle="--", alpha=0.5, linewidth=1)
            if max_eqtl_log10p > 5:
                axes[ax_idx].axhline(y=5, color="blue", linestyle="--", alpha=0.5, linewidth=1)
        
        # Get study info
        study_name = qtd_data['Study_Name'].iloc[0] if not qtd_data.empty else "No data"
        celltype = qtd_data['Celltype'].iloc[0] if not qtd_data.empty else ""
        
        axes[ax_idx].set_ylabel("-log₁₀(P)", fontsize=10)
        axes[ax_idx].set_title(f"{study_name} ({celltype})", fontsize=10, fontweight="bold")
        axes[ax_idx].grid(True, alpha=0.3)
        axes[ax_idx].spines["right"].set_visible(False)
        axes[ax_idx].spines["top"].set_visible(False)
    
    # Set x-axis label only for the bottom plot
    axes[-1].set_xlabel(f"Chr {chr_num} Position (Mb)", fontsize=12)
    
    # Add main title
    fig.suptitle(f"Manhattan Plots for Gene {gene_name} ({gene_id}) - Chromosome {chr_num}", 
                fontsize=14, fontweight="bold")
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)
    
    # Save plot
    output_path = f"{save_path}/manhattan_{gene_name}_{gene_id}_chr{chr_num}.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Manhattan plot saved to: {output_path}")
    
    plt.show()
    
    return all_qtd_data

# Main execution
if __name__ == "__main__":
    # Create output directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)
    
    # Process each gene in plot_gene_info
    for gene_id, chr_num, gene_name, pos_range in plot_gene_info:
        print(f"\nProcessing gene: {gene_name} ({gene_id}) on chromosome {chr_num}")
        
        # Create combined Manhattan plot (GWAS + all eQTL studies)
        qtd_data_list = create_combined_manhattan_plot(gene_id, chr_num, gene_name, pos_range)
        
        if qtd_data_list is not None:
            # Print summary statistics
            total_snps = sum(len(data) for data in qtd_data_list)
            print(f"\nSummary for {gene_name}:")
            print(f"Studies with data: {len(qtd_data_list)}")
            print(f"Total SNPs across all eQTL studies: {total_snps}")
            
            for i, data in enumerate(qtd_data_list):
                if not data.empty:
                    study_name = data['Study_Name'].iloc[0]
                    max_log10p = -np.log10(data['PVAL']).max()
                    n_sig = (data['PVAL'] < 5e-8).sum()
                    print(f"  {study_name}: {len(data)} SNPs, max -log10(P): {max_log10p:.2f}, significant: {n_sig}")