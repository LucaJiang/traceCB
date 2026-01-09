import pandas as pd
import glob

# 定义路径参数
input_dir = f"/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc"
output_file = f"{input_dir}/coloc_statistics_summary.csv"


# 获取所有结果文件（匹配*_coloc.csv）
file_list = glob.glob(f"{input_dir}/*_coloc.csv")
# chr, gene, n_snp_coloc, nsnp_eqtl, nsnp_gwas, p_original, p_traceC, p_traceCB, ...
# 存储统计结果
statistics = []

for file in file_list:
    df = pd.read_csv(file)
    # # unique gene
    # df = df.drop_duplicates(subset=["gene"])
    # ngene = df.shape[0]

    # 统计各p值列大于0.7的数量
    # updated column names based on run_coloc.r changes
    p_original_count = (df["p_original"] > 0.7).sum()
    p_traceC_count = (df["p_traceC"] > 0.7).sum()
    p_traceCB_count = (df["p_traceCB"] > 0.7).sum()
    pp_h3_original_count = (df["pp_h3_original"] > 0.7).sum()
    pp_h3_traceC_count = (df["pp_h3_traceC"] > 0.7).sum()
    pp_h3_traceCB_count = (df["pp_h3_traceCB"] > 0.7).sum()

    # 记录结果（包含文件名）
    statistics.append(
        {
            "filename": file.split("/")[-1].split("_coloc.")[0],
            "total_loci": df.shape[0],
            "p_original_gt0.7": p_original_count,
            "p_traceC_gt0.7": p_traceC_count,
            "p_traceCB_gt0.7": p_traceCB_count,
            "pp_h3_original_gt0.7": pp_h3_original_count,
            "pp_h3_traceC_gt0.7": pp_h3_traceC_count,
            "pp_h3_traceCB_gt0.7": pp_h3_traceCB_count,
        }
    )

# 转换为DataFrame并保存结果
result_df = pd.DataFrame(statistics).sort_values(by="filename")
result_df.to_csv(output_file, index=False)
print(f"统计完成！结果已保存至: {output_file}")
# print all results
print(result_df)
