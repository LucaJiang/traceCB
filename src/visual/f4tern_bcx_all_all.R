# plt f4 ternary with labeled genes

library(ggtern)
library(gridExtra)
library(cowplot)
library(dplyr)

target_trait <- "bcx"
target_trait_name <- "BCX Monocyte Counts"
target_celltype <- ""
label_gene <- ""
label_gene_name <- ""

result_parh <-
  "/Users/lucajiang/learn/CityU/xpmm/coloc/data"
meta_data <-
  jsonlite::fromJSON("/Users/lucajiang/learn/CityU/xpmm/data/metadata.json")
save_path <-
  "/Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/other"

all_files <-
  list.files(result_parh, pattern = "coloc.csv", full.names = TRUE)

all_results <- data.frame()
for (file in all_files) {
  current_data <- read.csv(file)
  if (nrow(current_data) == 0) {
    next
  }
  #   "/Users/lucajiang/learn/CityU/xpmm/coloc/data/bcx_mon_QTD000115_NK_cells_coloc.csv"
  filename <- basename(file)
  if (target_celltype != "") {
    match_result <- regexec(
      paste0(
        "^",
        target_trait,
        "_(\\w*)_(QTD\\d+)_",
        target_celltype,
        "_coloc\\.csv$"
      ),
      filename
    )
  } else {
    match_result <- regexec(
      paste0("^", target_trait, "_(\\w*)_(QTD\\d+)\\S*_coloc\\.csv$"),
      filename
    )
  }
  
  if (match_result[[1]][1] != -1) {
    # split the filename by _
    filename_parts <- strsplit(filename, "_")[[1]]
    current_data$trait <- filename_parts[2]
    current_data$qtdid <- filename_parts[3]
    if (length(all_results) == 0) {
      all_results <- current_data
    } else {
      all_results <- rbind(all_results, current_data)
    }
  }
}
head(all_results)
#  [1] "chr"                    "gene"                   "start"                  "end"
#  [5] "nsnp_eqtl"              "nsnp_gwas"              "n_snp_coloc"            "p_single"
#  [9] "p_gmm_cross"            "p_gmm_cross_tissue"     "pp_h0_single"           "pp_h1_single"
# [13] "pp_h2_single"           "pp_h3_single"           "pp_h0_gmm_cross"        "pp_h1_gmm_cross"
# [17] "pp_h2_gmm_cross"        "pp_h3_gmm_cross"        "pp_h0_gmm_cross_tissue" "pp_h1_gmm_cross_tissue"
# [21] "pp_h2_gmm_cross_tissue" "pp_h3_gmm_cross_tissue" "qtdid"

all_results$h4_improve <-
  all_results$p_gmm_cross - all_results$p_single
all_results$h4_tissue_improve <-
  all_results$p_gmm_cross_tissue - all_results$p_gmm_cross

all_results$h3_improve <-
  all_results$pp_h3_gmm_cross - all_results$pp_h3_single
all_results$h3_improve_tissue <-
  all_results$pp_h3_gmm_cross_tissue - all_results$pp_h3_gmm_cross

all_results$qtdname <- meta_data$id2name[all_results$qtdid]
all_results$qtdname <- as.factor(unlist(all_results$qtdname))
all_results$celltype <- meta_data$id2celltype[all_results$qtdid]
all_results$celltype <- as.factor(unlist(all_results$celltype))

# find cell type specific h4 improvement
all_celltypes <- unique(all_results$celltype)

# 找到细胞类型特异性改善的位点，比较最大和第二大改善
find_celltype_specific_improvement <- function(data, min_difference = 0.05) {
  
  # 为每个基因计算细胞类型特异性改善
  genes <- unique(data$gene)
  
  result_df <- data.frame(
    gene = character(),
    trait = character(),
    a = character(),
    a.improve = numeric(),
    b = character(),
    b.improve = numeric(),
    difference = numeric(),
    qtdname_a = character(),
    qtdname_b = character(),
    stringsAsFactors = FALSE
  )
  
  for (current_gene in genes) {
    # 获取当前基因在所有细胞类型中的数据
    gene_data <- data[data$gene == current_gene, ]
    
    if (length(unique(gene_data$celltype)) > 1) {  # 确保在多个细胞类型中都有数据
      
      # 计算每个细胞类型的平均h4_improve
      celltype_improvements <- aggregate(h4_improve ~ celltype + qtdname, 
                                       data = gene_data, 
                                       FUN = mean)
      
      # 按改善程度排序
      celltype_improvements <- celltype_improvements[order(celltype_improvements$h4_improve, decreasing = TRUE), ]
      
      if (nrow(celltype_improvements) >= 2) {
        # 获取最大和第二大改善
        max_celltype <- celltype_improvements$celltype[1]
        max_improvement <- celltype_improvements$h4_improve[1]
        max_qtdname <- celltype_improvements$qtdname[1]
        max_trait <- celltype_improvements$trait[1]
        
        second_celltype <- celltype_improvements$celltype[2]
        second_improvement <- celltype_improvements$h4_improve[2]
        second_qtdname <- celltype_improvements$qtdname[2]
        
        improvement_difference <- max_improvement - abs(second_improvement)
        
        # 检查是否满足最小差异条件
        if (improvement_difference > min_difference) {
          result_df <- rbind(result_df, data.frame(
            gene = current_gene,
            a = as.character(max_celltype),
            a.improve = max_improvement,
            b = as.character(second_celltype),
            b.improve = second_improvement,
            difference = improvement_difference,
            qtdname_a = as.character(max_qtdname),
            qtdname_b = as.character(second_qtdname),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # 按差异排序
  result_df <- result_df[order(result_df$difference, decreasing = TRUE), ]
  
  return(result_df)
}

# 查找细胞类型特异性改善的位点
celltype_specific_results <- find_celltype_specific_improvement(all_results, min_difference = 0.05)

# 显示结果
cat("发现", nrow(celltype_specific_results), "个细胞类型特异性改善的基因:\n\n")

if (nrow(celltype_specific_results) > 0) {
  print(celltype_specific_results)
  
  # 显示前几个结果的详细信息
  cat("\n前20个差异最大的结果:\n")
  for (i in 1:min(20, nrow(celltype_specific_results))) {
    row <- celltype_specific_results[i, ]
    cat("基因:", row$gene, "\n")
    cat("  最大改善细胞类型:", row$a, "(", row$qtdname_a, ")", "- 改善:", round(row$a.improve, 4), "\n")
    cat("  第二大改善细胞类型:", row$b, "(", row$qtdname_b, ")", "- 改善:", round(row$b.improve, 4), "\n")
    cat("  差异:", round(row$difference, 4), "\n\n")
  }
}
all_results[all_results$gene == "ENSG00000100298", ]
all_results[all_results$gene == "ENSG00000162104", ]
all_results[all_results$gene == "ENSG00000168575", ]
all_results[all_results$gene == "ENSG00000178878", ]
all_results[all_results$gene == "ENSG00000090621", ]
all_results[all_results$gene == "ENSG00000115919", ]

# 前20个差异最大的结果:
# 基因: ENSG00000100298 APOBEC3H, mon, 7
#   最大改善细胞类型: CD4+T_cells ( CEDAR(290) ) - 改善: 0.9909 
#   第二大改善细胞类型: CD8+T_cells ( CEDAR(277) ) - 改善: 0.1689 
#   差异: 0.822 

# 基因: ENSG00000162104 ADCY9, mon, 4
#   最大改善细胞类型: Monocytes ( Fairfax_2014(420) ) - 改善: 0.7826 
#   第二大改善细胞类型: NK_cells ( Gilchrist_2021(247) ) - 改善: 0 
#   差异: 0.7826 

# 基因: ENSG00000168575 SLC20A2, mcv, 3
#   最大改善细胞类型: NK_cells ( Gilchrist_2021(247) ) - 改善: 0.7564 
#   第二大改善细胞类型: CD8+T_cells ( Kasela_2017(269) ) - 改善: 0.109 
#   差异: 0.6474 

# 基因: ENSG00000178878 APOLD1, mon&wbc, 10
#   最大改善细胞类型: Monocytes ( BLUEPRINT(191) ) - 改善: 0.7241 
#   第二大改善细胞类型: CD8+T_cells ( CEDAR(277) ) - 改善: 0.084 
#   差异: 0.6401

# 基因: ENSG00000090621 PABPC4, mch&mcv, 6 
#   最大改善细胞类型: CD8+T_cells ( CEDAR(277) ) - 改善: 0.6822 
#   第二大改善细胞类型: Monocytes ( CEDAR(286) ) - 改善: 0.1273  #improve in tissue
#   差异: 0.5549 

# 基因: ENSG00000115919 KYNU, mon, 4 unmatched
#   最大改善细胞类型: Monocytes ( Fairfax_2014(420) ) - 改善: 0.6984 
#   第二大改善细胞类型: Monocytes ( CEDAR(286) ) - 改善: 0.161 
#   差异: 0.5375 

apply_tern_settings <- function(p, title, x_col, y_col) {
  
  result <- p + # can not avoid warning
    stat_density_tern(
      aes(alpha = after_stat(level)),
      bdl = 0.15,
      bdl.val = 0.15,
      h = c(2, 2),
      expand = c(2, 2),
      color = "#456aff",
    ) +
    geom_point(size = 1) +
    theme_light() +
    theme_nomask() +
    Tlab("", "Coloc") +
    Llab("", "Independent") +
    Rlab("", "Undetermined") +
    theme_showarrows() +
    labs(title = title,
         x = "",
         y = "",
         z = "") +
    theme(plot.title = element_text(size = 11),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"))
    
  return(result)
}

all_results <- all_results[all_results$gene == "ENSG00000162104", ]

p1 <- ggtern(data = all_results,
             aes(
               x = pp_h3_single,
               y = p_single,
               z = 1 - pp_h3_single - p_single,
               color = celltype
             )) %>%
  apply_tern_settings("Original", "pp_h3_single", "p_single") +
  theme(legend.position = "none")
p2 <- ggtern(
  data = all_results,
  aes(
    x = pp_h3_gmm_cross,
    y = p_gmm_cross,
    z = 1 - pp_h3_gmm_cross - p_gmm_cross,
    color = celltype
  )
) %>%
  apply_tern_settings("traceC", "pp_h3_gmm_cross", "p_gmm_cross") +
  theme(legend.position = "none")
p3 <- ggtern(
  data = all_results,
  aes(
    x = pp_h3_gmm_cross_tissue,
    y = p_gmm_cross_tissue,
    z = 1 - pp_h3_gmm_cross_tissue - p_gmm_cross_tissue,
    color = celltype
  )
) %>%
  apply_tern_settings("traceCB",
                      "pp_h3_gmm_cross_tissue",
                      "p_gmm_cross_tissue") +
  theme(legend.position = "none")

print(p1)
print(p2)
print(p3)

# plots_row <- grid.arrange(p1, p2, p3, ncol = 3)
plots_row <- plot_grid(
  p1,
  p2,
  p3,
  ncol = 3,
  align = "hv",
  rel_widths = c(1, 1, 1)
)

# 创建研究图例
standalone_legend <-
  ggplot(all_results, aes(x = 1, y = 1, color = celltype)) +
  geom_point() +
  scale_color_discrete(name = "Study") +
  guides(color = guide_legend(
    title = "Cell Type",
    nrow = 1,
    byrow = TRUE,
    title.position = "left"
  )) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "lines"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

study_legend <- ggpubr::get_legend(standalone_legend)


# 主标题设置
title <- ggdraw() +
  draw_label(
    paste0(
      "Colocalization Results for ",
      target_trait_name,
      " in All Study"
    ),
    size = 14,
    x = 0.5,
    hjust = 0.5,
  )

final_plot <- plot_grid(
  title,
  # 顶部标题
  study_legend,
  # 合并的图例
  plots_row,
  # 三个图表
  ncol = 1,
  # 垂直排列
  greedy = TRUE,
  rel_heights = c(0.04, 0.06, 0.40)  # 调整高度比例
)
print(final_plot)
ggsave(
  filename = paste0(
    save_path,
    "/",
    target_celltype,
    "_",
    target_trait,
    "_coloc_ternary_all.png"
  ),
  plot = final_plot,
  width = 10,
  height = 4,
  dpi = 300,
  units = "in"
)

