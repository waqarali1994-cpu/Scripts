
getwd()
setwd('~/Documents/PhD work/Yield_project/Bugeater_snp_effect/')

# Path to your VCF file (can be gzipped)
vcf_file <- "2021_snps_v2.vcf"  # adjust as needed

# Prefix for PLINK output
plink_prefix <- "plink_with_indels_21"



# Construct PLINK command (make sure 'plink' is in your PATH)
plink_cmd <- paste(
  "plink",
  "--vcf", vcf_file,
  "--make-bed",
  "--double-id",          # This line fixes the problem!
  "--out", plink_prefix
)


# Run the command
system(plink_cmd)









setwd('~/Documents/PhD work/kmer_gwas/kmers_p value issue/')


library(rMVP)
library(data.table)
library(bigmemory)
library(parallel)



setwd('~/Documents/PhD work/Rice_project/')
# Full-featured function (Recommended)
MVP.Data(fileVCF="IR64_merge_GWAS_filtered.vcf",
         #filePhe="Phenotype.txt",
         fileKin=TRUE,
         filePC=TRUE,
         out="IR64_filtered_mvp.vcf"
)

# Only convert genotypes
#MVP.Data.VCF2MVP("myVCF.vcf", out='mvp') # the genotype data should be fully imputed before using this function


# Detect available cores
num_cpus <- parallel::detectCores()

# Define traits to analyze
traits <- c("Control")

# Output directory
results_dir <- "results/MLM_pvalue_rice"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Read in data
genotype      <- attach.big.matrix("IR64_filtered_mvp.geno.desc")
phenotype <- fread("Pheno_Control_filtered.csv")
map           <- read.table("IR64_filtered_mvp.geno.map", header = TRUE)
Kinship       <- attach.big.matrix("IR64_filtered_mvp.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("IR64_filtered_mvp.pc.desc"))

# Bonferroni threshold
n_snps <- nrow(map)
bonf_p <- 0.05 / n_snps

# GWAS loop
for (trait in traits) {
  cat("Running MLM for trait:", trait, "\n")
  
  phe_trait <- phenotype[, c("genotype", trait), with = FALSE]
  colnames(phe_trait)[1] <- "genotype"  # rMVP expects "Taxa" as the sample ID column
  
  imMVP <- MVP(
    phe         = phe_trait,
    geno        = genotype,
    map         = map,
    K           = Kinship,
    CV.MLM      = Covariates_PC,   # Use external PCs
    # nPC.MLM    = 3,              # DO NOT use this with CV.MLM
    vc.method   = "EMMA",
    method      = "MLM",
    maxLoop     = 10,
    file.output = FALSE,
    p.threshold = bonf_p,
    ncpus       = num_cpus
  )
  
  # Save results
  output_file <- file.path(results_dir, paste0("GWAS_IR64", trait, ".csv.gz"))
  fwrite(cbind(imMVP$map, imMVP$mlm.results), file = output_file)
  
  cat("✔️ GWAS for", trait, "done →", output_file, "\n")
}



















library(rMVP)
library(data.table)
library(bigmemory)

setwd('~/Documents/PhD work/Yield project/PCA_MLM//')
# Attach necessary data from mvp.plink files
genotype <- attach.big.matrix("mvp.plink.geno.desc")
phenotype <- read.csv("Phenotype_matched.csv", header = TRUE, sep = ",")
map <- read.table("mvp.plink.geno.map", head = TRUE)
Kinship <- attach.big.matrix("mvp.plink.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.plink.pc.desc"))

# Define traits to analyze
traits <- c("PC1", "PC2", "PC3")  # Update with actual trait names if different

# Define results directory
results_dir <- "results/MLM"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Run GWAS using MLM Model for each trait
for (trait in traits) {
  imMVP <- MVP(
    phe = phenotype[, c("Genotype", trait)],
    geno = genotype,
    map = map,
    nPC.MLM = 3,
    K = Kinship,
    vc.method = "EMMA",
    ncpus = 16,
    maxLoop = 10,
    method = "MLM",
    file.output = FALSE,
    p.threshold = 0.05 / 14570744
  )
  
  # Save GWAS results
  output_file <- file.path(results_dir, paste0("GWAS_", trait, ".csv.gz"))
  fwrite(cbind(imMVP$map, imMVP$mlm.results), file = output_file)
  
  print(paste("GWAS for", trait, "completed. Results saved in", output_file))
}


mlm <- fread('GWAS_NumberRows.csv.gz')
mlm$NumberRows.MLM







library(dplyr)
library(data.table)
library(ggplot2)
setwd('~/Documents/PhD work/Rice_project/results/MLM_pvalue_rice//')

# -------------------------------
# Load GWAS results
# -------------------------------
mlm <- fread("GWAS_SN.csv.gz")


#chr3_pval <- min(mlm[CHROM == "chr3", AN.MLM], na.rm = TRUE)
#chr8_pval <- min(mlm[CHROM == "chr8", AN.MLM], na.rm = TRUE)

colnames(mlm)


df <- mlm %>%
  mutate(
    CHROM = as.numeric(gsub("chr", "", CHROM)),  # ensure chromosome numeric
    logP  = -log10(SN.MLM)
  ) %>%
  arrange(CHROM, POS)

# -------------------------------
# Filter SNPs above threshold
# -------------------------------
threshold <- 2   # i.e. p < 1e-6
df_sig <- df %>% filter(logP >= threshold)


# Compute cumulative base pair positions
data_cum <- df %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0))

df_sig <- inner_join(df_sig, data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add)

axis_set <- df %>%
  inner_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add) %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum))

chr_limits <- df %>%
  inner_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add) %>%
  group_by(CHROM) %>%
  summarise(start = min(bp_cum), end = max(bp_cum))

# -------------------------------
# Manhattan Plot
# -------------------------------
p <- ggplot(df_sig, aes(x = bp_cum, y = logP)) +
  geom_rect(data = chr_limits,
            aes(xmin = start, xmax = end, ymin = 0, ymax = 0),
            fill = rep(c("darkgrey", "lightgrey"), length.out = nrow(chr_limits)),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_point(color = "red", alpha = 0.8, size = 2) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
  scale_x_continuous(breaks = axis_set$center,
                     labels = paste0("Chr", sprintf("%02d", axis_set$CHROM))) +
  # Option 1: Remove expand (recommended)
  scale_y_continuous(
    breaks = seq(0, ceiling(max(df_sig$logP)), by = 1),
    minor_breaks = seq(0, ceiling(max(df_sig$logP)), by = 0.5),
    expand = expansion(mult = c(0, 0.05))  # Only 5% expansion at top
  ) +
  labs(x = "SN", y = expression(-log[10](italic(p)))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title  = element_text(size = 22)
  )

print(p)





library(dplyr)
library(data.table)
library(ggplot2)
setwd('~/Documents/PhD work/kmer_gwas/kmers_p value issue/results/MLM_pvalue_flowering//')

# -------------------------------
# Load GWAS results
# -------------------------------
mlm <- fread("GWAS_AN.csv.gz")
colnames(mlm)

df <- mlm %>%
  mutate(
    CHROM = as.numeric(gsub("chr", "", CHROM)),  # ensure chromosome numeric
    logP  = -log10(AN.MLM)
  ) %>%
  arrange(CHROM, POS)

# -------------------------------
# Calculate Bonferroni threshold
# -------------------------------
n_snps <- nrow(df)  # Total number of SNPs tested
bonferroni_pval <- 0.05 / n_snps
bonferroni_threshold <- -log10(bonferroni_pval)

cat("Total SNPs tested:", n_snps, "\n")
cat("Bonferroni p-value:", bonferroni_pval, "\n")
cat("Bonferroni -log10(p) threshold:", bonferroni_threshold, "\n")

# -------------------------------
# Filter SNPs above threshold
# -------------------------------
threshold <- 2   # suggestive threshold (p < 0.01)
df_sig <- df %>% filter(logP >= threshold)

# Compute cumulative base pair positions
data_cum <- df %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0))

df_sig <- inner_join(df_sig, data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add)

axis_set <- df %>%
  inner_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add) %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum))

chr_limits <- df %>%
  inner_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add) %>%
  group_by(CHROM) %>%
  summarise(start = min(bp_cum), end = max(bp_cum))

# -------------------------------
# Manhattan Plot
# -------------------------------
p <- ggplot(df_sig, aes(x = bp_cum, y = logP)) +
  geom_rect(data = chr_limits,
            aes(xmin = start, xmax = end, ymin = 0, ymax = 0),
            fill = rep(c("darkgrey", "lightgrey"), length.out = nrow(chr_limits)),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_point(color = "red", alpha = 0.8, size = 2) +
  #geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", linewidth = 1) +  # Suggestive threshold
  geom_hline(yintercept = bonferroni_threshold, linetype = "solid", color = "red", linewidth = 1) +  # Bonferroni threshold
  scale_x_continuous(breaks = axis_set$center,
                     labels = paste0("Chr", sprintf("%02d", axis_set$CHROM))) +
  scale_y_continuous(
    breaks = seq(0, ceiling(max(df_sig$logP)), by = 1),
    minor_breaks = seq(0, ceiling(max(df_sig$logP)), by = 0.5),
    expand = expansion(mult = c(0, 0.05))  # Only 5% expansion at top
  ) +
  labs(x = "AN", y = expression(-log[10](italic(p)))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title  = element_text(size = 22)
  )

print(p)

# Optional: Save the plot
ggsave("Manhattan_plot_SN_bonferroni.png", p, width = 14, height = 6, dpi = 300)







