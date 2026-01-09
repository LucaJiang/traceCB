library(coloc)
library(data.table)
library(dplyr)
library(stringr)

# cmd: Rscript src/run_coloc.r /gpfs1/scratch/wjiang49/XeQTL/coloc/bcx/bcx_mon_loci.csv /gpfs1/scratch/wjiang49/XeQTL/coloc/bcx/bcx_mon_GWAS.csv /gpfs1/scratch/wjiang49/XeQTL/BLUEPRINT/gmm/Monocytes /gpfs1/scratch/wjiang49/XeQTL/coloc/bcx/bcx_mon_coloc.csv > /gpfs1/scratch/wjiang49/XeQTL/coloc/bcx/bcx_mon.log

args <- commandArgs(trailingOnly = TRUE)
loci.path <- args[1]
gwas.path <- args[2]
eqtl.path <- args[3]
save.path <- args[4]

posterior.prob <- 0.70

# read loci
loci <- fread(loci.path) %>%
    mutate(
        chr = str_replace(chrom, "chr", ""),
        snp = GWAS_snp,
        gene = ensembl
    ) %>%
    dplyr::select(chr, gene, start, end)
n.loci <- nrow(loci)

# read gwas data
gwas <- fread(gwas.path)
#  SNP Z CHR POS MAF BETA SE

# initialize results df
# Using new nomenclature: original (S), traceC (C), traceCB (T)
results <- tibble(
    chr = character(),
    gene = character(),
    start = numeric(),
    end = numeric(),
    nsnp_eqtl = numeric(),
    nsnp_gwas = numeric(),
    n_snp_coloc = numeric(),
    p_original = numeric(),
    p_traceC = numeric(),
    p_traceCB = numeric(),
    pp_h0_original = numeric(),
    pp_h1_original = numeric(),
    pp_h2_original = numeric(),
    pp_h3_original = numeric(),
    pp_h0_traceC = numeric(),
    pp_h1_traceC = numeric(),
    pp_h2_traceC = numeric(),
    pp_h3_traceC = numeric(),
    pp_h0_traceCB = numeric(),
    pp_h1_traceCB = numeric(),
    pp_h2_traceCB = numeric(),
    pp_h3_traceCB = numeric()
)

n.coloc.genes <- 0
for (i in 1:n.loci) {
    loci.data <- loci[i, ]

    # read eqtl data of target gene
    if (!file.exists(paste0(eqtl.path, "/chr", loci.data$chr, "/", loci.data$gene, ".csv"))) {
        next
    }
    eqtl.gene <-
        fread(paste0(eqtl.path, "/chr", loci.data$chr, "/", loci.data$gene, ".csv"))
    n.snp.eqtl <- nrow(eqtl.gene)

    # get gwas data of target locus
    gwas.locus <-
        gwas %>% filter(CHR == loci.data$chr & !is.na(BETA) & !is.na(SE))
    n.snp.gwas <- nrow(gwas.locus)

    # get common snps
    snp.common <- intersect(eqtl.gene$RSID, gwas.locus$SNP)
    n.snp.coloc <- length(snp.common)

    # coloc
    if (n.snp.coloc > 10) {
        n.coloc.genes <- n.coloc.genes + 1
        results <- results %>%
            add_row(
                chr = loci.data$chr,
                gene = loci.data$gene,
                start = loci.data$start,
                end = loci.data$end,
                nsnp_eqtl = n.snp.eqtl,
                nsnp_gwas = n.snp.gwas,
                n_snp_coloc = n.snp.coloc
            )
        gwas.coloc <-
            gwas.locus %>%
            filter(SNP %in% snp.common) %>%
            distinct(SNP, .keep_all = TRUE) %>%
            arrange(SNP)

        eqtl.coloc <-
            eqtl.gene %>%
            filter(RSID %in% snp.common) %>%
            distinct(RSID, .keep_all = TRUE) %>%
            arrange(RSID)

        # Mapping: original -> S, traceC -> C, traceCB -> T
        targets <- c("original", "traceC", "traceCB")
        suffixes <- c("S", "C", "T")

        for (j in seq_along(targets)) {
            target <- targets[j]
            suffix <- suffixes[j]
            
            target.coloc.res <- coloc.abf(
                dataset1 = list(
                    beta = gwas.coloc$BETA,
                    varbeta = gwas.coloc$SE^2,
                    snp = gwas.coloc$SNP,
                    type = "cc"
                ),
                dataset2 = list(
                    beta = eqtl.coloc[[paste0("TAR_", suffix, "BETA")]],
                    varbeta = eqtl.coloc[[paste0("TAR_", suffix, "SE")]]^2,
                    snp = eqtl.coloc$RSID,
                    sdY = 1,
                    type = "quant"
                )
            )
            
            results[n.coloc.genes, paste0("p_", target)] <-
                target.coloc.res$summary[["PP.H4.abf"]]
            results[n.coloc.genes, paste0("pp_h0_", target)] <-
                target.coloc.res$summary[["PP.H0.abf"]]
            results[n.coloc.genes, paste0("pp_h1_", target)] <-
                target.coloc.res$summary[["PP.H1.abf"]]
            results[n.coloc.genes, paste0("pp_h2_", target)] <-
                target.coloc.res$summary[["PP.H2.abf"]]
            results[n.coloc.genes, paste0("pp_h3_", target)] <-
                target.coloc.res$summary[["PP.H3.abf"]]
        }
    }
}
# delete rows with NA
results <- results %>% filter(!is.na(p_original))

# find genes with high posterior probability in 3 methods
for (method in c("original", "traceC", "traceCB")) {
    n.method.coloc <-
        sum(results[[paste0("p_", method)]] > posterior.prob)
    print(paste0(method, ": ", n.method.coloc))
}
# find number of row where original != traceC
n.cross.sign <- sum(results$p_original != results$p_traceC)
print(paste0("candidate genes: ", n.loci, ", available genes: ", n.coloc.genes, ", cross-significant genes: ", n.cross.sign))
print(
    results %>%
        select(chr, gene, n_snp_coloc, p_original, p_traceC, p_traceCB),
    n = Inf
)

results.loci <- results %>%
    mutate(chr = as.integer(chr)) %>%
    arrange(chr, gene)

write.csv(results.loci, save.path, row.names = FALSE, quote = FALSE)
