# plot Locuszoom use Grch 37/hg19

# ----- H3K27ac -----
library(rtracklayer)
library(ggplot2)
library(gggenes)
library(colorspace)
library(gridExtra)
library(GenomicFeatures)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(ggtext) # 添加ggtext包用于支持富文本
data(SLE_gwas_sub) # data from locuszoomr package

save.path <-
  "/Users/lucajiang/learn/CityU/traceCB/data/img/eas_eqtlgen"
# 导入bigWig数据
base.path <- "/Users/lucajiang/learn/CityU/xpmm/data/Locuszoom/"
meta_data <- jsonlite::fromJSON("/Users/lucajiang/learn/CityU/traceCB/src/visual/metadata.json")
celltype_colors <- meta_data$celltype_colors
color_mapping <- c(
  "B" = darken(celltype_colors$B_cells, amount = 0.4),
  "Mon" = darken(celltype_colors$Monocytes, amount = 0.4),
  "NK" = darken(celltype_colors$NK_cells, amount = 0.4),
  "CD4T" = darken(celltype_colors$`CD4+T_cells`, amount = 0.4),
  "CD8T" = darken(celltype_colors$`CD8+T_cells`, amount = 0.4)
)


plot_locuszoom <- function(gene_names, gene_infos, chrs, start_positions, end_positions, track_maxrow=6) {
  window_size <- 100000 # 100kb window size
  regions <-
    paste0(
      "chr",
      chrs,
      ":",
      start_positions - window_size,
      "-",
      end_positions + window_size
    )

  # load meta_data.json

  # color_mapping <- c(
  #   "B" = celltype_colors$B_cells,
  #   "Mon" = celltype_colors$Monocytes,
  #   "NK" = celltype_colors$NK_cells,
  #   "CD4T" = celltype_colors$`CD4+T_cells`,
  #   "CD8T" = celltype_colors$`CD8+T_cells`
  # )

  # 创建细胞类型标签映射（使用expression）
  celltype_labels <- c(
    "B" = "B~cells",
    "Mon" = "Monocytes",
    "NK" = "NK~cells",
    "CD4T" = "CD4^'+'~T~cells",
    "CD8T" = "CD8^'+'~T~cells"
  )

  for (i in 1:length(regions)) {
    region <- GRanges(regions[i])
    gene_name <- gene_names[i]
    gene_info <- gene_infos[i]
    ### --------------------------
    ### 第一部分：信号峰图
    ### --------------------------
    # From ENCODE https://www.encodeproject.org/search/?type=File&searchTerm=H3K27ac+CD4+positive+T+cell&file_type=bigWig&biosample_ontology.cell_slims=CD4%2B+T+cell
    b <-
      import(paste0(base.path, "B_ENCFF701BIL.bigWig"), which = region)
    monocytes <-
      import(paste0(base.path, "MON_ENCFF840HBF.bigWig"), which = region)
    cd4 <-
      import(paste0(base.path, "CD4_ENCFF357NOB.bigWig"), which = region)
    cd8 <-
      import(paste0(base.path, "CD8_ENCFF455UVQ.bigWig"), which = region)
    nk <-
      import(paste0(base.path, "NK_ENCFF473CXT.bigWig"), which = region)

    # 创建绘图数据框
    plot_data <- rbind(
      data.frame(
        position = start(b),
        signal = score(b),
        cell = "B"
      ),
      data.frame(
        position = start(monocytes),
        signal = score(monocytes),
        cell = "Mon"
      ),
      data.frame(
        position = start(nk),
        signal = score(nk),
        cell = "NK"
      ),
      data.frame(
        position = start(cd4),
        signal = score(cd4),
        cell = "CD4T"
      ),
      data.frame(
        position = start(cd8),
        signal = score(cd8),
        cell = "CD8T"
      )
    )
    plot_data$cell <-
      factor(plot_data$cell, levels = names(color_mapping))
    signal_plot <-
      ggplot(plot_data, aes(x = position / 1e6, y = signal)) +
      geom_line(aes(color = cell), linewidth = 0.6, alpha = 1) +
      scale_color_manual(values = color_mapping) +
      facet_grid(cell ~ .,
        scales = "fixed", switch = "y",
        labeller = as_labeller(celltype_labels, default = label_parsed)
      ) +
      xlab("") +
      ylab("H3K27ac Signal Intensity") +
      # theme_minimal() +
      ggtitle(gene_name) +
      theme(
        strip.text.y.left = element_text(
          angle = 0,
          hjust = 0.5,
          margin = margin(r = 0, l = 0) # 减小标签右侧间隙
        ),
        strip.background.y = element_blank(),
        # 减少标签的左右边距),
        legend.position = "none",
        # 移除图例，因为facet已经显示了标签
        panel.spacing.y = unit(0.2, "lines"),
        # 减少子图间距
        axis.text.y = element_blank(),
        # 去掉y轴刻度文字
        axis.ticks.y = element_blank(),
        # 去掉y轴刻度线
        # axis.text.x = element_blank(),
        # 去掉x轴刻度文字
        # axis.ticks.x = element_blank(),
        # 去掉x轴刻度线
        # strip.background = element_blank(),
        strip.switch.pad.grid = unit(0, "pt"),
        # 减少switch模式下的间距
        panel.grid.minor = element_blank(),
        # 去掉次要网格线
        panel.grid.major = element_blank(),
        # 去掉主要网格线
        panel.background = element_blank(),
        # 去掉面板背景
        plot.margin = margin(
          t = 2,
          r = 5,
          b = 0,
          l = 5,
          unit = "pt"
        )
      )
    print(signal_plot)
    ggsave(
      paste0(save.path, "/", gene_name, "_H3K27ac_signal.pdf"),
      plot = signal_plot,
      width = 10,
      height = 3.5, # 减少高度
      dpi = 300,
      device = "pdf"
    )
    # ---------------------
    ### 第二部分：基因结构图
    ### ---------------------
    # 使用locuszoomr包绘制基因结构图
    loc <-
      locus(
        SLE_gwas_sub,
        gene = gene_name,
        flank = window_size,
        LD = "r2",
        ens_db = "EnsDb.Hsapiens.v75"
      )
    track <-
      gg_genetracks(
        loc,
        maxrows = track_maxrow,
        filter_gene_biotype = "protein_coding",
        cex.text = 1.2
      ) +
      theme(
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.margin = margin(t = 0, r = 5, b = 2, l = 5, unit = "pt") # 减少顶部边距
      )

    # # 组合信号峰图和基因结构图
    # combined_plot <- grid.arrange(
    #   signal_plot,
    #   track,
    #   ncol = 1,
    #   heights = c(4, 0.8)  # 减少第二个图的高度从1到0.8
    # )

    # 保存图像
    ggsave(
      paste0(save.path, "/", gene_name, "_genetrack.pdf"),
      plot = track,
      width = 12,
      height = 2, # 适当减少总高度
      dpi = 300,
      device = "pdf"
    )
  }
}

# genes to plot
gene_names <- c("CTSW")
gene_infos <- c(
  "ENSG00000172543, 11: 65,647,280-65,651,212"
)
chrs <- c(11) # chromosomes for the genes
start_positions <- c(65647280) # start positions for the genes
end_positions <- c(65651212) # end positions for the genes
plot_locuszoom(gene_names, gene_infos, chrs, start_positions, end_positions, 5)

# genes to plot
gene_names <- c("RABGAP1")
gene_infos <- c(
  "ENSG00000011454, 9: 125,703,112-125,867,145"
)
chrs <- c(9) # chromosomes for the genes
start_positions <- c(125703112) # start positions for the genes
end_positions <- c(125867145) # end positions for the genes
plot_locuszoom(gene_names, gene_infos, chrs, start_positions, end_positions)

# genes to plot
gene_names <- c("ADD3")
gene_infos <- c(
  "ENSG00000148700, 10: 111,756,126-111,895,323"
)
chrs <- c(10) # chromosomes for the genes
start_positions <- c(111756126) # start positions for the genes
end_positions <- c(111895323) # end positions for the genes
plot_locuszoom(gene_names, gene_infos, chrs, start_positions, end_positions)

# genes to plot
gene_names <- c("KYNU")
gene_infos <- c(
  "ENSG00000115919, 2: 143,635,067-143,799,890"
)
chrs <- c(2) # chromosomes for the genes
start_positions <- c(143635067) # start positions for the genes
end_positions <- c(143799890) # end positions for the genes
plot_locuszoom(gene_names, gene_infos, chrs, start_positions, end_positions)

# genes to plot
gene_names <- c("USP35")
gene_infos <- c(
  "ENSG00000118369, 11:77,899,858-77,925,757"
)
chrs <- c(11) # chromosomes for the genes
start_positions <- c(77899858) # start positions for the genes
end_positions <- c(77925757) # end positions for the genes
plot_locuszoom(gene_names, gene_infos, chrs, start_positions, end_positions)

# genes to plot
gene_names <- c("WDR48")
gene_infos <- c(
  "ENSG00000114742, 3: 39,093,489-39,138,155"
)
chrs <- c(3) # chromosomes for the genes
start_positions <- c(39093489) # start positions for the genes
end_positions <- c(39138155) # end positions for the genes
plot_locuszoom(gene_names, gene_infos, chrs, start_positions, end_positions)
