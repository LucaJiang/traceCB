library(coloc)
library(data.table)
library(dplyr)
library(stringr)

# cmd: Rscript src/run_coloc_replicate.r /gpfs1/scratch/wjiang49/XeQTL/coloc/bcx/bcx_mon_loci.csv /gpfs1/scratch/wjiang49/XeQTL/coloc/bcx/bcx_mon_GWAS.csv  /gpfs1/scratch/wjiang49/XeQTL/coloc/replicate/replicate_eqtl.csv /gpfs1/scratch/wjiang49/XeQTL/coloc/replicate/bcx_mon_replicate_tissue_coloc.csv > /gpfs1/scratch/wjiang49/XeQTL/coloc/replicate/bcx_mon_replicate.log

# set path
# setwd("/gpfs1/scratch/wjiang49/XeQTL")
# loci.path <- "coloc/eur_ms/eur_ms_loci.csv"
# gwas.path <- "coloc/eur_ms/eur_ms_GWAS.csv" # SNP,CHR,POS,Z
# replicate.path <- "coloc/replicate/replicate_eqtl.csv"
# save.path <- "coloc/replicate/eur_ms_replicate_tissue_coloc.csv"

args <- commandArgs(trailingOnly = TRUE)
loci.path <- args[1]
gwas.path <- args[2]
replicate.path <- args[3]
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

# read replicate data
replicate <- fread(replicate.path) %>%
    # capitalize column names
    setnames(names(.), toupper(names(.))) %>%
    # rename columns
    rename(SNP = RSID) %>%
    # select columns
    select(CHR, POS, GENE, SNP, BETA, SE) %>%
    # remove rows with NA
    filter(!is.na(BETA) & !is.na(SE))
# chr,pos,variant_id,ref,alt,gene,gene_name,beta,se,pval,pip_susie,rsid
# 1,100000012,1_100000012_G_T,G,T,ENSG00000162688,AGL,0.0505631,0.0160197,0.00164814,0.0011359478923747,rs10875231

# initialize results df
results <- tibble(
    chr = character(),
    gene = character(),
    start = numeric(),
    end = numeric(),
    nsnp_gwas = numeric(),
    nsnp_replicate = numeric(),
    n_snp_coloc = numeric(),
    p_h0 = numeric(),
    p_h1 = numeric(),
    p_h2 = numeric(),
    p_h3 = numeric(),
    p_h4 = numeric(),
)

n.coloc.genes <- 0
for (i in 1:n.loci) {
    loci.data <- loci[i, ]

    # get replicate data of target locus
    replicate.locus <-
        replicate %>% filter(
            GENE == loci.data$gene &
                CHR == loci.data$chr &
                !is.na(BETA) & !is.na(SE)
        )
    n.snp.replicate <- nrow(replicate.locus)

    if (n.snp.replicate < 10) {
        next
    }

    # get gwas data of target locus
    gwas.locus <-
        gwas %>% filter(CHR == loci.data$chr & !is.na(BETA) & !is.na(SE))
    n.snp.gwas <- nrow(gwas.locus)

    # get common snps
    snp.common <- intersect(gwas.locus$SNP, replicate.locus$SNP)
    n.snp.coloc <- length(snp.common)
    # random select 1/10 snps

    # coloc
    if (n.snp.coloc > 10) {
        n.coloc.genes <- n.coloc.genes + 1
        results <- results %>%
            add_row(
                chr = loci.data$chr,
                gene = loci.data$gene,
                start = loci.data$start,
                end = loci.data$end,
                nsnp_gwas = n.snp.gwas,
                nsnp_replicate = n.snp.replicate,
                n_snp_coloc = n.snp.coloc
            )
        replicate.coloc <-
            replicate.locus %>%
            filter(SNP %in% snp.common) %>%
            distinct(SNP, .keep_all = TRUE) %>% # 筛选唯一的SNP，并保留所有其他列
            arrange(SNP)

        gwas.coloc <-
            gwas.locus %>%
            filter(SNP %in% snp.common) %>%
            distinct(SNP, .keep_all = TRUE) %>% # 筛选唯一的SNP，并保留所有其他列
            arrange(SNP)

        coloc.results <- coloc.abf(
            dataset1 = list(
                beta = gwas.coloc$BETA,
                varbeta = gwas.coloc$SE^2,
                snp = gwas.coloc$SNP,
                type = "cc"
            ),
            dataset2 = list(
                beta = replicate.coloc$BETA,
                varbeta = replicate.coloc$SE^2,
                snp = replicate.coloc$SNP,
                sdY = 1,
                type = "quant"
            )
        )
        results[n.coloc.genes, "p_h0"] <- coloc.results$summary[["PP.H0.abf"]]
        results[n.coloc.genes, "p_h1"] <- coloc.results$summary[["PP.H1.abf"]]
        results[n.coloc.genes, "p_h2"] <- coloc.results$summary[["PP.H2.abf"]]
        results[n.coloc.genes, "p_h3"] <- coloc.results$summary[["PP.H3.abf"]]
        results[n.coloc.genes, "p_h4"] <- coloc.results$summary[["PP.H4.abf"]]
    }
}
# find genes with high posterior probability in 3 methods
n.coloc <- sum(results$p_h4 > posterior.prob)
print(
    paste0(
        "candidate loci: ",
        n.loci,
        ", available loci: ",
        n.coloc.genes,
        ", coloc loci: ",
        n.coloc
    )
)

results.loci <- results %>%
    mutate(chr = as.integer(chr)) %>%
    arrange(chr, gene)

write.csv(results.loci,
    save.path,
    row.names = FALSE,
    quote = FALSE
)

head(results.loci)
