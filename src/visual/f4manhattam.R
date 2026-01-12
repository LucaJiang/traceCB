#!/usr/bin/env Rscript

# Load required libraries
library(locuszoomr)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(readr)
library(data.table)

# Configuration (adapt from your Python utils)
study_path_main <- "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen"  # Update this path
save_path <- "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/manhattan"        # Update this path
# Configuration
MIN_PVAL <- 1e-40
MAX_RANGE <- 200000
gwas_path <- "/home/wjiang49/group/wjiang49/data/traceCB/coloc/bcx/bcx_mon_GWAS.csv"
meta_data <- jsonlite::fromJSON("/home/wjiang49/traceCB/src/visual/metadata.json")

# LD token for LDlink API
token <- "72edb9cc22c9"

# Plot gene info (same as your Python version)
plot_gene_info <- list(
  c("ENSG00000172543", "QTD000021", 11, "CTSW"),
  c("ENSG00000172543", "QTD000081", 11, "CTSW"),
  c("ENSG00000172543", "QTD000069", 11, "CTSW")
)

# Function to convert Z-score to p-value
z2p <- function(z) {
  2 * pnorm(-abs(z))
}

# Function to find index SNP (most significant SNP)
find_index_snp <- function(df, pval_col) {
  min_pval_idx <- which.min(df[[pval_col]])
  return(df$RSID[min_pval_idx])
}

# Main plotting function
plot_manhattan_locuszoom <- function(gene_id, qtdid, chromosome, gene_name) {
  cat(sprintf("Plotting manhattan for %s (%s) in %s, chr%d\n", gene_name, gene_id, qtdid, chromosome))
  
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
  
  # Clip position range if too large
  min_pos <- min(merged_df$POS)
  max_pos <- max(merged_df$POS)
  if (max_pos - min_pos > MAX_RANGE) {
    mid_pos <- (min_pos + max_pos) / 2
    new_min_pos <- mid_pos - MAX_RANGE / 2
    new_max_pos <- mid_pos + MAX_RANGE / 2
    merged_df <- merged_df %>% 
      filter(POS >= new_min_pos & POS <= new_max_pos)
    cat(sprintf("Clipped position range for %s in %s: %g - %g from %g - %g\n", 
                gene_name, qtdid, new_min_pos, new_max_pos, min_pos, max_pos))
  }
  
  # Sort by position
  merged_df <- merged_df %>% arrange(POS)
  
  # Define target methods and labels
  target_methods <- c("PVAL", "TAR_SPVAL", "TAR_CPVAL", "TAR_TPVAL")
  method_labels <- c(
    "GWAS",
    meta_data$method_name[1],
    meta_data$method_name[2], 
    meta_data$method_name[3]
  )
  
  # Clip p-values to avoid log10(0)
  for (col in target_methods) {
    merged_df[[col]] <- pmax(merged_df[[col]], MIN_PVAL)
  }
  
  # Create plots list
  plots <- list()
  
  for (i in seq_along(target_methods)) {
    method <- target_methods[i]
    label <- method_labels[i]
    
    # Find index SNP for this method
    index_snp <- find_index_snp(merged_df, method)
    cat(sprintf("Index SNP for %s: %s\n", label, index_snp))
    
    # Prepare data for locuszoomr
    plot_data <- merged_df %>%
      select(rsid = RSID, chrom = CHR, pos = POS, pval = all_of(method)) %>%
      mutate(
        chrom = as.character(chrom),
        logp = -log10(pval)
      ) %>%
      as.data.table()
    
    # Initialize p as NULL
    p <- NULL
    
    # Create locus object
    tryCatch({
      # Try to get LD data from LDlink
      loc <- locuszoomr::locus(
        data = plot_data,
        gene = gene_name,
        flank = 50000,
        LD = index_snp,
        ens_db = "EnsDb.Hsapiens.v75"
      )
      
      # Link LD data using your token
      loc <- locuszoomr::link_LD(loc, token = token, pop = "EAS")
      
      # Create the plot with LD coloring
      p <- locuszoomr::locus_plot(loc, 
                     labels = index_snp,
                     LD_scheme = c("blue", "skyblue", "green", "orange", "red")) +
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
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")
        )
      
    }, error = function(e) {
      cat(sprintf("Failed to get LD data for %s, using fallback method: %s\n", label, e$message))
      
      # Fallback: create plot without LD data
      p <<- ggplot(plot_data, aes(x = pos / 1e6, y = logp)) +
        geom_point(size = 0.8, alpha = 0.7, color = "steelblue") +
        # Highlight index SNP
        geom_point(data = plot_data[plot_data$rsid == index_snp, ], 
                  aes(x = pos / 1e6, y = logp), 
                  color = "red", size = 2, shape = 18) +
        geom_hline(yintercept = 5, linetype = "dashed", color = "blue", alpha = 0.5, linewidth = 0.5) +
        geom_hline(yintercept = 8, linetype = "dashed", color = "red", alpha = 0.5, linewidth = 0.5) +
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
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")
        ) +
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", x),
          breaks = scales::pretty_breaks(n = 6)
        )
    })
    
    # Check if p was created successfully
    if (is.null(p)) {
      cat(sprintf("Warning: Plot creation failed for %s, skipping...\n", label))
      next
    }
    
    plots[[i]] <- p
  }
  
  # Check if we have any plots
  if (length(plots) == 0) {
    cat("Error: No plots were created successfully\n")
    return(NULL)
  }
  
  # Combine plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots)))
  
  # Add overall title
  final_plot <- combined_plot + 
    plot_annotation(
      title = sprintf("Manhattan Plot of %s in GWAS and %s", 
                     gene_name, meta_data$id2name[[qtdid]]),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Save plot
  save_name <- sprintf("%s/%s_%s_%s_manhattan_LD.pdf", 
                      save_path, qtdid, gene_id, gene_name)
  
  ggsave(
    save_name,
    final_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat(sprintf("Manhattan plot saved to: %s\n", save_name))
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
    # Create dummy r2 values (you would replace this with real LD data)
    merged_df$r2 <- runif(nrow(merged_df), 0, 1)
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
                   LD = "r2") +
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
      title = sprintf("Manhattan Plot of %s in GWAS and %s (with LD)", 
                     gene_name, meta_data$id2name[[qtdid]]),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  save_name <- sprintf("%s/%s_%s_%s_manhattan_precomputed_LD.pdf", 
                      save_path, qtdid, gene_id, gene_name)
  
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

# Create save directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Main execution
for (gene_info in plot_gene_info) {
  gene_id <- gene_info[1]
  qtdid <- gene_info[2] 
  chromosome <- as.numeric(gene_info[3])
  gene_name <- gene_info[4]
  
  # Try plotting with error handling
  tryCatch({
    save_name <- plot_manhattan_locuszoom(gene_id, qtdid, chromosome, gene_name)
    if (!is.null(save_name)) {
      cat(sprintf("Successfully created plot: %s\n", save_name))
    }
  }, error = function(e) {
    cat(sprintf("Failed to create plot for %s: %s\n", gene_name, e$message))
    cat(sprintf("Traceback: %s\n", paste(sys.calls(), collapse = "\n")))
  })
}

cat("Manhattan plot generation completed\n")
