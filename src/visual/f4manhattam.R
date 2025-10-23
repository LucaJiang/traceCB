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


study_path_main <- "/Users/lucajiang/learn/CityU/xpmm/data/EAS_GTEx"
save_path <- "/Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/locuszoom"
epi_track_path <- "/Users/lucajiang/learn/CityU/xpmm/data/Locuszoom/"

# Configuration
MIN_PVAL <- 1e-20
gwas_path <- "/Users/lucajiang/learn/CityU/XeQTL/data/coloc/bcx_mon/bcx_mon_GWAS.csv"
meta_data <- jsonlite::fromJSON("/Users/lucajiang/learn/CityU/xpmm/data/metadata.json")
token <- "72edb9cc22c9" # LD token for LDlink API

# Plot gene info
plot_gene_info <- list(
  # # CTSW (chromosome 11)
  # c("ENSG00000172543", "QTD000021", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000069", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000081", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000031", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000067", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000371", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000066", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000372", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000073", 11, "CTSW"),
  # c("ENSG00000172543", "QTD000115", 11, "CTSW"),
  # 
  # # RABGAP1 (chromosome 9)
  # c("ENSG00000011454", "QTD000021", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000069", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000081", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000031", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000067", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000371", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000066", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000372", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000073", 9, "RABGAP1"),
  # c("ENSG00000011454", "QTD000115", 9, "RABGAP1"),
  # 
  # # ADD3 (chromosome 10)
  # c("ENSG00000148700", "QTD000021", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000069", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000081", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000031", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000067", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000371", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000066", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000372", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000073", 10, "ADD3"),
  # c("ENSG00000148700", "QTD000115", 10, "ADD3"),
  # 
  # # USP35 (chromosome 11)
  # c("ENSG00000118369", "QTD000021", 11, "USP35"),
  # c("ENSG00000118369", "QTD000069", 11, "USP35"),
  # c("ENSG00000118369", "QTD000081", 11, "USP35"),
  # c("ENSG00000118369", "QTD000031", 11, "USP35"),
  # c("ENSG00000118369", "QTD000067", 11, "USP35"),
  # c("ENSG00000118369", "QTD000371", 11, "USP35"),
  # c("ENSG00000118369", "QTD000066", 11, "USP35"),
  # c("ENSG00000118369", "QTD000372", 11, "USP35"),
  # c("ENSG00000118369", "QTD000073", 11, "USP35"),
  # c("ENSG00000118369", "QTD000115", 11, "USP35"),
  
  # # WDR48
  # c("ENSG00000114742", "QTD000021", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000031", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000066", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000067", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000069", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000073", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000081", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000115", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000371", 3, "WDR48"),
  # c("ENSG00000114742", "QTD000372", 3, "WDR48"),

  # KYNU (chromosome 2)
  # c("ENSG00000115919", "QTD000021", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000031", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000066", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000067", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000069", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000073", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000081", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000115", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000371", 2, "KYNU"),
  # c("ENSG00000115919", "QTD000372", 2, "KYNU"),

  # # CD14 (chromosome 5)
  # c("ENSG00000170458", "QTD000021", 5, "CD14"),
  # c("ENSG00000170458", "QTD000031", 5, "CD14"),
  # c("ENSG00000170458", "QTD000066", 5, "CD14"),
  # c("ENSG00000170458", "QTD000067", 5, "CD14"),
  # c("ENSG00000170458", "QTD000069", 5, "CD14"),
  # c("ENSG00000170458", "QTD000073", 5, "CD14"),
  # c("ENSG00000170458", "QTD000081", 5, "CD14"),
  # c("ENSG00000170458", "QTD000115", 5, "CD14"),
  # c("ENSG00000170458", "QTD000371", 5, "CD14"),
  # c("ENSG00000170458", "QTD000372", 5, "CD14"),

  # ZYX (chromosome 7)
  c("ENSG00000159840", "QTD000081", 7, "ZYX")

  # # ADCY9 (chromosome 16)
  # c("ENSG00000162104", "QTD000021", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000031", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000066", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000067", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000069", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000073", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000081", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000115", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000371", 16, "ADCY9"),
  # c("ENSG00000162104", "QTD000372", 16, "ADCY9")
)

# # Single gene for testing
# gene_info <- c("ENSG00000172543", "QTD000021", 11, "CTSW")

gene_max_range <- list(
  "CTSW" = 2e5,
  "RABGAP1" = 6e5,
  "ADD3" = 6e5,  
  "USP35" = 3e5,
  "WDR48" = 4e5,
  "KYNU" = 6e5,
  "CD14" = 2e5,
  "ZYX" = 5e5,
  "ADCY9" = 6e5
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
  eqtl_path <- sprintf("%s/%s/GMM/chr%d/%s.csv", study_path_main, qtdid, chromosome, gene_id)
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
    cat(sprintf("Clipped position range for %s in %s: %g - %g from %g - %g\n", 
                gene_name, qtdid, new_min_pos, new_max_pos, min_pos, max_pos))
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
  save_name <- sprintf("%s/%s_%s_%s_%s_manhattan_LD.pdf", 
                      save_path, gene_name, cell_type, qtdid, gene_id)
  if (file.exists(save_name)) {
    file.remove(save_name)
    # return(c(save_name, min_pos, max_pos))
  }
  pdf(save_name, width = 8, height = 7)
  
  # Set up layered plot with 5 layers (4 manhattan plots + 1 gene track)
  oldpar <- set_layers(5)
  
  # Create plots for each method
  for (i in seq_along(target_methods)) {
    method <- target_methods[i]
    label <- method_labels[i]
    # Prepare data for locuszoomr
    plot_data <- merged_df %>%
      dplyr::select(rsid = RSID, chrom = CHR, pos = POS, p = all_of(method)) %>%
      dplyr::filter(!is.na(p)) %>%
      dplyr::mutate(p = pmax(p, MIN_PVAL)) %>%  # Apply MIN_PVAL as floor
      as.data.table()
    print(head(plot_data))
    # Create locus object
    loc <- locus(
      data = plot_data,
      ens_db = "EnsDb.Hsapiens.v75",
      gene = gene_name
    )
    
    # Link LD data using your token
    loc <- link_LD(loc, token = token, pop = "EAS")  # Using East Asian population
    
    # Create scatter plot for each method
    if (i == 1) {
      legend_pos = "topleft"
    } else {
      legend_pos = NULL
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
  genetracks(loc, filter_gene_biotype = "protein_coding", 
            maxrows = 2, highlight = gene_name, italics = TRUE)
  
  mtext(sprintf("Manhattan Plot of %s in GWAS and %s %s", 
                  gene_name, meta_data$id2name[[qtdid]], cell_type),
          side = 1, line = 5, cex = 1.2, font = 2)
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
    ggtitle(paste0(gene_name, " - H3K27ac Signal")) +
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
    height = 3.5,  # 减少高度
    # dpi = 300
  )
  return(save_name)
}

# Main execution
for (gene_info in plot_gene_info) {
  gene_id <- gene_info[1]
  qtdid <- gene_info[2] 
  chromosome <- as.numeric(gene_info[3])
  gene_name <- gene_info[4]
  
  # Get max_range from the dictionary
  max_range <- gene_max_range[[gene_name]]
  if (is.null(max_range)) {
    max_range <- 2e5  # Default value if gene not found
  }
  
  ret <- plot_manhattan_locuszoom(gene_id, qtdid, chromosome, gene_name, max_range = max_range)
  
  if (is.null(ret)) {
    cat(sprintf("Failed to plot for %s (%s) in %s\n", gene_name, gene_id, qtdid))
  } else {
    cat(sprintf("Plot saved: %s\n", ret[1]))
    signal_save_name <- plot_signal_plot(chromosome, as.numeric(ret[2]), as.numeric(ret[3]), gene_name)
    if (!is.null(signal_save_name)) {
      cat(sprintf("Signal plot saved: %s\n", signal_save_name))
    }
    cat(sprintf("Finished processing %s (%s) in %s\n", gene_name, gene_id, qtdid))
  }
}

