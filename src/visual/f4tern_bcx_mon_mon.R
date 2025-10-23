# plt f4 ternary with labeled genes

library(ggtern)
library(gridExtra)
library(cowplot)
library(dplyr)
library(colorspace)

target_trait <- "bcx_mon"
target_trait_name <- "BCX Monocyte Counts"
target_celltype <- "Monocytes"
label_gene <- "ENSG00000172543"
label_gene_name <- "CTSW"

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
  match_result <- regexec(
    paste0(
      "^",
      target_trait,
      "_(QTD\\d+)_",
      target_celltype,
      "_coloc\\.csv$"
    ),
    filename
  )
  
  if (match_result[[1]][1] != -1) {
    # split the filename by _
    filename_parts <- strsplit(filename, "_")[[1]]
    qtdid <- filename_parts[3]
    current_data$qtdid <- qtdid
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

# find h3 and h4 improve less than MIN_IMPROVE_THRESHOLD
MIN_IMPROVE_THRESHOLD <- 0.01
all_results <- all_results %>%
  filter(
    h3_improve > MIN_IMPROVE_THRESHOLD |
      h4_improve > MIN_IMPROVE_THRESHOLD
  )

all_results$qtdname <- meta_data$id2name[all_results$qtdid]
all_results$qtdname <- as.factor(unlist(all_results$qtdname))

# 创建基因ID到基因名称的映射
gene_name_mapping <- setNames(label_gene_name, label_gene)

# 筛选出要标记的基因
label_data <- all_results[all_results$gene == label_gene,]
label_data$gene_name <- gene_name_mapping[label_data$gene]

# 添加标签名称映射
label_name_shorten <- list(
  "Monocytes" = "Monocyte",
  "CD4+T_cells" = "CD4T",
  "CD8+T_cells" = "CD8T",
  "B_cells" = "B",
  "NK_cells" = "NK"
)

apply_tern_settings <- function(p, title, x_col, y_col) {
  # 为标记数据创建临时列，包括z坐标
  temp_label_data <- label_data
  temp_label_data$plot_x <- temp_label_data[[x_col]]
  temp_label_data$plot_y <- temp_label_data[[y_col]]
  temp_label_data$plot_z <-
    1 - temp_label_data$plot_x - temp_label_data$plot_y
  
  result <- p + # can not avoid warning
    stat_density_tern(
      aes(alpha = after_stat(level)),
      bdl = 0.15,
      bdl.val = 0.15,
      h = c(2, 2),
      expand = c(2, 2),
      color = "#456aff",
    ) +
    geom_point(size = 0.8) +
    # 添加标记基因的空心圆圈 - 用固定颜色
    geom_point(
      data = temp_label_data[temp_label_data$gene_name == label_gene_name, ],
      aes(x = plot_x, y = plot_y, z = plot_z),
      color = "red",
      size = 2,
      shape = 21,
      stroke = 1,
      fill = NA,
      inherit.aes = FALSE
    ) +
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


p1 <- ggtern(data = all_results,
             aes(
               x = pp_h3_single,
               y = p_single,
               z = 1 - pp_h3_single - p_single,
               color = qtdname
             )) %>%
  apply_tern_settings("Original", "pp_h3_single", "p_single") +
  theme(legend.position = "none")
p2 <- ggtern(
  data = all_results,
  aes(
    x = pp_h3_gmm_cross,
    y = p_gmm_cross,
    z = 1 - pp_h3_gmm_cross - p_gmm_cross,
    color = qtdname
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
    color = qtdname
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
  ggplot(all_results, aes(x = 1, y = 1, color = qtdname)) +
  geom_point() +
  guides(color = guide_legend(
    title = "Study",
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

# 创建基因标记图例
gene_legend_plot <- ggplot(data.frame(gene = label_gene_name, x = 1, y = 1), 
                          aes(x = x, y = y, color = gene)) +
  geom_point(size = 3, shape = 21, stroke = 1.5, fill = NA) +
  scale_color_manual(values = setNames("red", label_gene_name), 
                    name = "Labeled Gene") +
  guides(color = guide_legend(
    title = "Gene",
    nrow = 1,
    byrow = TRUE,
    title.position = "left",
    override.aes = list(color = "red", size = 3, shape = 21, stroke = 1.5, fill = NA)
  )) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10, face = "italic"),
    legend.key.size = unit(0.6, "lines"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
# 提取图例
study_legend <- ggpubr::get_legend(standalone_legend)
gene_legend <- ggpubr::get_legend(gene_legend_plot)

# 合并图例
combined_legend <- plot_grid(
  study_legend, 
  gene_legend, 
  ncol = 2, 
  rel_widths = c(0.7, 0.3),  # 调整宽度比例
  align = "h",               # 水平对齐
  axis = "tb"              # 顶部和底部对齐
)

# 主标题设置
title <- ggdraw() +
  draw_label(
    paste0(
      "Colocalization Results for ",
      target_trait_name,
      " in ",
      target_celltype,
      " Study"
    ),
    size = 14,
    x = 0.5,
    hjust = 0.5,
  )

final_plot <- plot_grid(
  title,
  # 顶部标题
  combined_legend,
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
    "_coloc_ternary.png"
  ),
  plot = final_plot,
  width = 10,
  height = 4,
  dpi = 300,
  units = "in"
)

