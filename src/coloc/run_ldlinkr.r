library(LDlinkR)
library(dplyr)
library(readr)

# Load the data from command line arguments
args <- commandArgs(trailingOnly = TRUE)
leading_snp_file <- args[1]
output_path_filename <- args[2]

token <- "" # replace with your LDlink token !!!!!!!!!!!!!!
print(token)

# Read the leading SNP data
leading_snp <- read.csv(leading_snp_file, header = TRUE, sep = ",")
get_proxy_results <- function(snp, token) {
    proxy_result <- LDlinkR::LDproxy(snp,
        pop = "JPT", r2d = "r2", genome_build = "grch37",
        token = token
    )
    proxy_result[proxy_result$R2 >= 0.1, ]
}

proxy_results <- lapply(leading_snp$SNP, get_proxy_results, token = token)
# save proxy_results
# save(proxy_results, file = "proxy_results.RData")
loci <- lapply(proxy_results, function(x) {
    if (is.data.frame(x) && nrow(x) > 1) {
        chr <- x$Coord %>%
            gsub(":.+", "", .) %>%
            unique()
        pos <- x$Coord %>%
            gsub("chr[0-9]{1,2}:", "", .) %>%
            as.numeric()
        locus_name <- paste0(chr, ":", min(pos), "_", max(pos))
        top_snp_df <- dplyr::filter(x, Distance == 0)
        top_snp <- top_snp_df %>% pull(RS_Number)
        top_snp_pos <- top_snp_df %>% pull(Coord)
        locus <- tibble(
            GWAS_snp = top_snp,
            GWAS_snp_pos = top_snp_pos,
            locus_name = locus_name,
            chrom = chr,
            start = min(pos),
            end = max(pos),
            pval = leading_snp$PVALUE[leading_snp$SNP == top_snp],
        )
        return(locus)
    } else {
        return(NULL)
    }
}) %>%
    bind_rows() %>%
    write_tsv(., output_path_filename)
