#!/usr/bin/env Rscript

# Load required libraries
library(locuszoomr)
library(readr)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(readr)
library(data.table)
library(EnsDb.Hsapiens.v75)
library(dplyr)
library(rtracklayer)
library(colorspace)


# Configuration
study_path_main <- "/Users/lucajiang/learn/CityU/traceCB/data/EAS_eQTLGen" # Update this path
save_path <- "/Users/lucajiang/learn/CityU/traceCB/data/img/eas_eqtlgen" # Update this path
MIN_PVAL <- 1e-40
MAX_RANGE <- 200000
gwas_path <- "/Users/lucajiang/learn/CityU/XeQTL/data/coloc/bcx_mon/bcx_mon_GWAS.csv"
meta_data <- jsonlite::fromJSON("/Users/lucajiang/learn/CityU/traceCB/src/visual/metadata.json")
epi_track_path <- "/Users/lucajiang/learn/CityU/xpmm/data/Locuszoom/"

token <- "72edb9cc22c9" # LD token for LDlink API

# Plot gene info
plot_gene_info <- list(
  # CTSW (chromosome 11)
  # c("ENSG00000172543", "QTD000021", 11, "CTSW"),
  c("ENSG00000172543", "QTD000069", 11, "CTSW"),
  c("ENSG00000172543", "QTD000081", 11, "CTSW"),
  c("ENSG00000172543", "QTD000031", 11, "CTSW"),
  c("ENSG00000172543", "QTD000067", 11, "CTSW"),
  c("ENSG00000172543", "QTD000371", 11, "CTSW"),
  c("ENSG00000172543", "QTD000066", 11, "CTSW"),
  c("ENSG00000172543", "QTD000372", 11, "CTSW"),
  c("ENSG00000172543", "QTD000073", 11, "CTSW"),
  c("ENSG00000172543", "QTD000115", 11, "CTSW"),
  # KCNN4 (chromosome 19)
  c("ENSG00000104783", "QTD000021", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000069", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000081", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000031", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000067", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000371", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000066", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000372", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000073", 19, "KCNN4"),
  c("ENSG00000104783", "QTD000115", 19, "KCNN4"),
  # MDFIC (chromosome 7)
  c("ENSG00000135272", "QTD000021", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000069", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000081", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000031", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000067", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000371", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000066", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000372", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000073", 7, "MDFIC"),
  c("ENSG00000135272", "QTD000115", 7, "MDFIC"),
  # CTDP1 (chromosome 18)
  c("ENSG00000060069", "QTD000081", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000031", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000067", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000371", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000066", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000372", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000073", 18, "CTDP1"),
  c("ENSG00000060069", "QTD000115", 18, "CTDP1")
)


gene_max_range <- list(
  "CTSW" = 2e5,
  "KCNN4" = 2e5,
  "MDFIC" = 5e5,
  "CTDP1" = 10e5
)

# Color mapping for cell types
celltype_colors <- meta_data$celltype_colors
color_mapping <- c(
  "B" = darken(celltype_colors$B_cells, amount = 0.4),
  "Mon" = darken(celltype_colors$Monocytes, amount = 0.4),
  "NK" = darken(celltype_colors$NK_cells, amount = 0.4),
  "CD4T" = darken(celltype_colors$`CD4+T_cells`, amount = 0.4),
  "CD8T" = darken(celltype_colors$`CD8+T_cells`, amount = 0.4)
)

# Function to convert Z-score to p-value
z2p <- function(z) {
  2 * pnorm(-abs(z))
}

# Main plotting function
plot_manhattan_locuszoom <- function(gene_id, qtdid, chromosome, gene_name, max_range = 2e5) {
  cat(sprintf("Plotting manhattan for %s (%s) in %s, chr%d\n", gene_name, gene_id, qtdid, chromosome))

  # Load data
  eqtl_path <- sprintf("%s/%s/%s.csv", study_path_main, qtdid, gene_id)
  if (!file.exists(eqtl_path)) {
    cat(sprintf("eQTL file not found: %s\n", eqtl_path))
    return(NULL)
  }
  eqtl_df <- as_tibble(fread(eqtl_path, na.strings = c("", "NA", "nan", "NaN")))
  gwas_df <- as_tibble(fread(gwas_path, na.strings = c("", "NA", "nan", "NaN")))
  gwas_df$PVAL <- z2p(gwas_df$Z)

  # Merge data
  merged_df <- inner_join(eqtl_df, gwas_df, by = c("RSID" = "SNP"))

  # Clip position range if too large
  min_pos <- min(merged_df$POS)
  max_pos <- max(merged_df$POS)
  if (max_pos - min_pos > max_range) {
    mid_pos <- (min_pos + max_pos) / 2
    new_min_pos <- mid_pos - max_range / 2
    new_max_pos <- mid_pos + max_range / 2
    merged_df <- merged_df %>%
      dplyr::filter(POS >= new_min_pos & POS <= new_max_pos)
    cat(sprintf(
      "Clipped position range for %s in %s: %g - %g from %g - %g\n",
      gene_name, qtdid, new_min_pos, new_max_pos, min_pos, max_pos
    ))
    min_pos <- new_min_pos
    max_pos <- new_max_pos
  }
  merged_df <- merged_df %>% arrange(POS)

  # Define target methods and labels
  target_methods <- c("PVAL", "TAR_SPVAL", "TAR_CPVAL", "TAR_TPVAL")
  method_labels <- c(
    "GWAS",
    meta_data$method_name[1],
    meta_data$method_name[2],
    meta_data$method_name[3]
  )
  cell_type <- meta_data$id2celltype[[qtdid]]
  # Save plot to PDF with layered structure
  save_name <- sprintf(
    "%s/%s_%s_%s_%s_manhattan_LD.pdf",
    save_path, gene_name, cell_type, qtdid, gene_id
  )
  if (file.exists(save_name)) {
    file.remove(save_name)
    # return(c(save_name, min_pos, max_pos))
  }
  pdf(save_name, width = 8, height = 7)

  # Set up layered plot with 4 layers (4 manhattan plots + 1 gene track)
  oldpar <- set_layers(4)

  # Create plots for each method
  for (i in seq_along(target_methods)) {
    method <- target_methods[i]
    label <- method_labels[i]
    # Prepare data for locuszoomr
    plot_data <- merged_df %>%
      dplyr::select(rsid = RSID, chrom = CHR, pos = POS, p = all_of(method)) %>%
      dplyr::filter(!is.na(p)) %>%
      dplyr::mutate(p = pmax(p, MIN_PVAL)) %>% # Apply MIN_PVAL as floor
      as.data.table()
    print(head(plot_data))

    # Identify index SNP
    index_snp <- plot_data$rsid[which.min(plot_data$p)]

    # Create locus object
    loc <- locus(
      data = plot_data,
      ens_db = "EnsDb.Hsapiens.v75",
      gene = gene_name,
      index_snp = index_snp,
      seqname = chromosome
    )

    # Link LD data using your token
    loc <- link_LD(loc, token = token, pop = "EAS") # Using East Asian population

    # Create scatter plot for each method
    if (i == 1) {
      legend_pos <- "topleft"
    } else {
      legend_pos <- NULL
    }

    scatter_plot(loc,
      labels = "index",
      xticks = FALSE,
      legend_pos = legend_pos,
    )
    mtext(label, side = 3, line = 0.5, adj = 0.98, cex = 1.2, font = 1)
  }

  # Add gene track at the bottom
  # Use the last locus object for gene track
  genetracks(loc,
    filter_gene_biotype = "protein_coding",
    maxrows = 5, highlight = gene_name, italics = TRUE
  )

  # mtext(
  #   sprintf(
  #     "Manhattan Plot of %s in GWAS and %s %s",
  #     gene_name, meta_data$id2name[[qtdid]], cell_type
  #   ),
  #   side = 1, line = 5, cex = 1.2, font = 2
  # )
  # Revert par() settings
  par(oldpar)
  dev.off()

  cat(sprintf("Manhattan plot saved to: %s\n", save_name))
  return(c(save_name, min_pos, max_pos))
}

# Function for H3K27ac signal plot
plot_signal_plot <- function(chr, start_pos, end_pos, gene_name) {
  # start_pos <- start_pos - 8000 # add more space to align with Manhattan
  # end_pos <- end_pos + 8000
  save_name <- sprintf("%s/%s_H3K27ac_signal.pdf", save_path, gene_name)
  if (file.exists(save_name)) {
    file.remove(save_name)
    # return(NULL)
  }
  # Create region for signal data
  region <- GRanges(paste0("chr", chr, ":", start_pos, "-", end_pos))

  # Import bigWig files
  b <- import(paste0(epi_track_path, "B_ENCFF701BIL.bigWig"), which = region)
  monocytes <- import(paste0(epi_track_path, "MON_ENCFF840HBF.bigWig"), which = region)
  cd4 <- import(paste0(epi_track_path, "CD4_ENCFF357NOB.bigWig"), which = region)
  cd8 <- import(paste0(epi_track_path, "CD8_ENCFF455UVQ.bigWig"), which = region)
  nk <- import(paste0(epi_track_path, "NK_ENCFF473CXT.bigWig"), which = region)

  # Create plotting data frame
  plot_data <- rbind(
    data.frame(position = start(b), signal = score(b), cell = "B"),
    data.frame(position = start(monocytes), signal = score(monocytes), cell = "Mon"),
    data.frame(position = start(nk), signal = score(nk), cell = "NK"),
    data.frame(position = start(cd4), signal = score(cd4), cell = "CD4T"),
    data.frame(position = start(cd8), signal = score(cd8), cell = "CD8T")
  )

  plot_data$cell <- factor(plot_data$cell, levels = names(color_mapping))

  # Create signal plot
  signal_plot <- ggplot(plot_data, aes(x = position / 1e6, y = signal)) +
    geom_line(aes(color = cell), linewidth = 0.6, alpha = 1) +
    scale_color_manual(values = color_mapping) +
    facet_grid(cell ~ ., scales = "fixed", switch = "y") +
    scale_x_continuous(
      limits = c(start_pos / 1e6, end_pos / 1e6),
      expand = c(0, 0)
    ) +
    xlab(paste0("Chromosome ", chr, " position (Mb)")) +
    ylab("H3K27ac Signal") +
    # ggtitle(paste0(gene_name, " - H3K27ac Signal")) +
    theme_minimal() +
    theme(
      strip.background.y = element_blank(),
      legend.position = "none",
      panel.spacing.y = unit(0.2, "lines"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.switch.pad.grid = unit(0, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      axis.text.x = element_text(),
      axis.ticks.x = element_line()
    )
  ggsave(
    save_name,
    plot = signal_plot,
    width = 10,
    height = 3.5, # 减少高度
    # dpi = 300
  )
  return(save_name)
}

# Alternative function using pre-computed LD data (if available)
plot_manhattan_with_precomputed_LD <- function(gene_id, qtdid, chromosome, gene_name, ld_file = NULL) {
  cat(sprintf("Plotting manhattan with pre-computed LD for %s (%s) in %s, chr%d\n", gene_name, gene_id, qtdid, chromosome))

  # Load eQTL data
  eqtl_path <- sprintf("%s/%s/GMM/chr%d/%s.csv", study_path_main, qtdid, chromosome, gene_id)
  if (!file.exists(eqtl_path)) {
    cat(sprintf("eQTL file not found: %s\n", eqtl_path))
    return(NULL)
  }
  eqtl_df <- read_csv(eqtl_path, show_col_types = FALSE)

  # Load GWAS data
  gwas_df <- read_csv(gwas_path, show_col_types = FALSE)
  gwas_df$PVAL <- z2p(gwas_df$Z)

  # Merge data
  merged_df <- inner_join(eqtl_df, gwas_df, by = c("RSID" = "SNP"))

  # If LD file is provided, load it
  if (!is.null(ld_file) && file.exists(ld_file)) {
    ld_df <- read_csv(ld_file, show_col_types = FALSE)
    merged_df <- left_join(merged_df, ld_df, by = "RSID")
  } else {
    # error handling if LD file is not found
    stop("LD file not found or not provided.")
  }

  # Position clipping and sorting
  min_pos <- min(merged_df$POS)
  max_pos <- max(merged_df$POS)
  if (max_pos - min_pos > MAX_RANGE) {
    mid_pos <- (min_pos + max_pos) / 2
    new_min_pos <- mid_pos - MAX_RANGE / 2
    new_max_pos <- mid_pos + MAX_RANGE / 2
    merged_df <- merged_df %>%
      filter(POS >= new_min_pos & POS <= new_max_pos)
  }

  merged_df <- merged_df %>% arrange(POS)

  # Define target methods and labels
  target_methods <- c("PVAL", "TAR_SPVAL", "TAR_CPVAL", "TAR_TPVAL")
  method_labels <- c("GWAS", meta_data$method_name[1], meta_data$method_name[2], meta_data$method_name[3])

  # Clip p-values
  for (col in target_methods) {
    merged_df[[col]] <- pmax(merged_df[[col]], MIN_PVAL)
  }

  # Create plots list
  plots <- list()

  for (i in seq_along(target_methods)) {
    method <- target_methods[i]
    label <- method_labels[i]

    # Find index SNP
    index_snp <- find_index_snp(merged_df, method)

    # Prepare data
    plot_data <- merged_df %>%
      select(rsid = RSID, chrom = CHR, pos = POS, pval = all_of(method), r2) %>%
      mutate(
        chrom = as.character(chrom),
        logp = -log10(pval),
        ld_category = case_when(
          rsid == index_snp ~ "Index SNP",
          r2 >= 0.8 ~ "r² ≥ 0.8",
          r2 >= 0.6 ~ "0.6 ≤ r² < 0.8",
          r2 >= 0.4 ~ "0.4 ≤ r² < 0.6",
          r2 >= 0.2 ~ "0.2 ≤ r² < 0.4",
          TRUE ~ "r² < 0.2"
        )
      ) %>%
      as.data.table()

    # Use locuszoomr with LD
    loc <- locus(
      data = plot_data,
      chrom = chromosome,
      start = min(plot_data$pos),
      end = max(plot_data$pos),
      index_snp = index_snp,
      LD = "r2"
    )

    p <- locus_plot(loc,
      labels = index_snp,
      LD = "r2"
    ) +
      labs(
        title = label,
        y = "-log10(P)",
        x = if (i == length(target_methods)) sprintf("Chr %d Position (Mb)", chromosome) else ""
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        axis.title.x = if (i != length(target_methods)) element_blank() else element_text(),
        axis.text.x = if (i != length(target_methods)) element_blank() else element_text(),
        legend.position = if (i == 1) "right" else "none",
        panel.grid.major.y = element_line(color = "gray90", size = 0.3)
      )

    plots[[i]] <- p
  }

  # Combine and save
  combined_plot <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots)))
  final_plot <- combined_plot +
    plot_annotation(
      # title = sprintf("Manhattan Plot of %s in GWAS and %s (with LD)",
      #                 gene_name, meta_data$id2name[[qtdid]]),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  save_name <- sprintf(
    "%s/%s_%s_%s_manhattan_precomputed_LD.png",
    save_path, qtdid, gene_id, gene_name
  )

  ggsave(
    save_name,
    final_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  cat(sprintf("Manhattan plot with pre-computed LD saved to: %s\n", save_name))
  return(save_name)
}

# Main execution
for (gene_info in plot_gene_info) {
  gene_info <- unlist(gene_info)
  gene_id <- gene_info[1]
  qtdid <- gene_info[2]
  chromosome <- as.numeric(gene_info[3])
  gene_name <- gene_info[4]

  # Get max_range from the dictionary
  max_range <- gene_max_range[[gene_name]]
  if (is.null(max_range)) {
    max_range <- 2e5 # Default value if gene not found
  }

  ret <- plot_manhattan_locuszoom(gene_id, qtdid, chromosome, gene_name, max_range = max_range)

  if (is.null(ret)) {
    cat(sprintf("Failed to plot for %s (%s) in %s\n", gene_name, gene_id, qtdid))
  } else {
    cat(sprintf("Plot saved: %s\n", ret[1]))
    # signal_save_name <- plot_signal_plot(chromosome, as.numeric(ret[2]), as.numeric(ret[3]), gene_name)
    if (!is.null(signal_save_name)) {
      cat(sprintf("Signal plot saved: %s\n", signal_save_name))
    }
    cat(sprintf("Finished processing %s (%s) in %s\n", gene_name, gene_id, qtdid))
  }
}
