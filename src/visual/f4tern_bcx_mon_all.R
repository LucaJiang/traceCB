# plt f4 ternary with labeled genes

library(ggtern)
library(gridExtra)
library(cowplot)
library(dplyr)
library(colorspace)

target_trait <- "bcx_mon"
target_trait_name <- "BCX Monocyte Counts"
target_celltype <- ""
label_gene <- ""
label_gene_name <- ""

result_parh <-
  "/Users/lucajiang/learn/CityU/traceCB/data/coloc"
meta_data <-
  jsonlite::fromJSON("/Users/lucajiang/learn/CityU/traceCB/src/visual/metadata.json")
save_path <-
  "/Users/lucajiang/learn/CityU/traceCB/data/img/eas_eqtlgen"

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
        "_eQTLGen_(QTD\\d+)_",
        target_celltype,
        "_coloc\\.csv$"
      ),
      filename
    )
  } else {
    match_result <- regexec(
      paste0("^", target_trait, "_eQTLGen_(QTD\\d+)\\S*_coloc\\.csv$"),
      filename
    )
  }

  if (match_result[[1]][1] != -1) {
    # split the filename by _
    filename_parts <- strsplit(filename, "_")[[1]]
    qtdid <- filename_parts[4]
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
  all_results$p_traceC - all_results$p_original
all_results$h4_tissue_improve <-
  all_results$p_traceCB - all_results$p_traceC

all_results$h3_improve <-
  all_results$pp_h3_traceC - all_results$pp_h3_original
all_results$h3_improve_tissue <-
  all_results$pp_h3_traceCB - all_results$pp_h3_traceC

all_results$qtdname <- meta_data$id2name[all_results$qtdid]
all_results$qtdname <- as.factor(unlist(all_results$qtdname))
all_results$celltype <- meta_data$id2celltype[all_results$qtdid]
all_results$celltype <- as.factor(unlist(all_results$celltype))

# find h3 and h4 improve less than MIN_IMPROVE_THRESHOLD
MIN_IMPROVE_THRESHOLD <- 0.01
all_results <- all_results %>%
  filter(
    h3_improve > MIN_IMPROVE_THRESHOLD |
      h4_improve > MIN_IMPROVE_THRESHOLD
  )


all_results$qtdname <- meta_data$id2name[all_results$qtdid]
all_results$qtdname <- as.factor(unlist(all_results$qtdname))
all_results$celltype <- meta_data$id2celltype[all_results$qtdid]
all_results$celltype <- as.factor(unlist(all_results$celltype))

# 添加标签名称映射
label_name_shorten <- list(
  "Monocytes" = "Monocytes",
  "CD4+T_cells" = "CD4^'+'~T~cells",
  "CD8+T_cells" = "CD8^'+'~T~cells",
  "B_cells" = "B~cells",
  "NK_cells" = "NK~cells"
)

# 创建缩短的标签列
all_results$celltype_short <- factor(
  sapply(all_results$celltype, function(x) label_name_shorten[[as.character(x)]]),
  levels = unname(unlist(label_name_shorten))
)

celltype_colors <- meta_data$celltype_colors
for (i in seq_along(celltype_colors)) {
  celltype_colors[[i]] <- darken(celltype_colors[[i]], amount = 0.3)
}
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
    geom_point(size = 0.3, alpha = 1) +
    scale_color_manual(values = celltype_colors) +
    theme_light() +
    theme_nomask() +
    Tlab("", "Coloc") +
    Llab("", "Independent") +
    Rlab("", "Undetermined") +
    theme_showarrows() +
    labs(
      title = title,
      x = "",
      y = "",
      z = ""
    ) +
    theme(
      plot.title = element_text(size = 11),
      plot.margin = margin(-10, -10, -10, -10, unit = "pt"),
      tern.panel.background = element_rect(fill = NA, colour = NA)
    )

  return(result)
}


p1 <- ggtern(
  data = all_results,
  aes(
    x = pp_h3_original,
    y = p_original,
    z = 1 - pp_h3_original - p_original,
    color = celltype
  )
) %>%
  apply_tern_settings("Original", "pp_h3_original", "p_original") +
  theme(legend.position = "none")
p2 <- ggtern(
  data = all_results,
  aes(
    x = pp_h3_traceC,
    y = p_traceC,
    z = 1 - pp_h3_traceC - p_traceC,
    color = celltype
  )
) %>%
  apply_tern_settings("traceC", "pp_h3_traceC", "p_traceC") +
  theme(legend.position = "none")
p3 <- ggtern(
  data = all_results,
  aes(
    x = pp_h3_traceCB,
    y = p_traceCB,
    z = 1 - pp_h3_traceCB - p_traceCB,
    color = celltype
  )
) %>%
  apply_tern_settings(
    "traceCB",
    "pp_h3_traceCB",
    "p_traceCB"
  ) +
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
  rel_widths = c(1, 1, 1),
  axis = "tb",
  hjust = 0,
  vjust = 0
)

# 创建研究图例
standalone_legend <-
  ggplot(all_results, aes(x = 1, y = 1, color = celltype)) +
  geom_point() +
  scale_color_manual(
    values = celltype_colors,
    labels = parse(text = unname(unlist(label_name_shorten)))
  ) +
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
  rel_heights = c(0.04, 0.06, 0.40) # 调整高度比例
)
print(final_plot)
ggsave(
  filename = paste0(
    save_path,
    "/bcx_moncount_coloc_ternary.pdf"
  ),
  plot = final_plot,
  width = 10,
  height = 4,
  dpi = 300,
  units = "in",
  device = "pdf"
)

