# Load required libraries
library(hierfstat)
library(adegenet)
library(dplyr)
library(poppr)
library(ggplot2)

# Working directory to script location
setwd("msc-codes-fmb_genetic_diversity/amova-Ho-He/scripts/")

# Output directory
output_dir <- "../results/haplopes"
dir.create(output_dir, showWarnings = FALSE)

# Load input data - Used previously for structure analysis
population_data <- read.delim("../data/pop_file.txt", header = FALSE)
genotype_data <- read.delim("../data/174300_r_snps.txt", header = FALSE)

# ----------------------
# Filter to 166 haplotypes only * Check UGA set and update 
# ----------------------

haplotypes_166 <- c(
  "E1", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E2", "E20", "E21", "E22", "E25", "E27", "E28", "E29", "E3",
  "E30", "E31", "E33", "E37", "E38", "E39", "E4", "E40", "E41", "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49",
  "E50", "E51", "E52", "E53", "E54", "E57", "E58", "E61", "E63",
  "K1", "K13", "K14", "K15", "K16", "K17", "K18", "K20", "K21", "K23", "K24", "K26", "K27", "K28", "K29", "K31", "K32",
  "K33", "K34", "K35", "K37", "K38", "K39", "K4", "K40", "K41", "K42", "K44", "K6", "K7", "K8", "K9",
  "T11", "T12", "T15", "T16", "T17", "T18", "T19", "T2", "T20", "T22", "T24", "T25", "T26", "T29", "T3", "T30", "T31",
  "T33", "T34", "T35", "T37", "T38", "T4", "T40", "T42", "T43", "T45", "T46", "T47", "T48", "T49", "T5", "T50", "T51",
  "T53", "T54", "T55", "T56", "T57", "T58", "T6", "T7", "T8", "T9",
  "U1", "U13", "U14", "U15", "U2", "U21", "U22", "U24", "U25", "U26", "U27", "U3", "U31", "U32", "U34", "U35", "U37",
  "U38", "U4", "U40", "U42", "U44", "U45", "U46", "U5", "U51", "U56", "U57", "U8",
  "UK10", "UK11", "UK12", "UK13", "UK14", "UK15", "UK2", "UK21", "UK22", "UK24", "UK25", "UK3", "UK4", "UK5", "UK6",
  "UK7", "UK9"
)

# Filter population data
population_data <- population_data %>% filter(V2 %in% haplotypes_166)

# Filter genotype data to match individuals retained in population_data
genotype_data <- genotype_data %>% filter(V1 %in% population_data$V1)

# ----------------------

# Metadata
individual_ids <- as.character(population_data$V1)
population_labels <- as.character(population_data$V2)

# Prepare SNP genotype matrix (excluding label column)
snp_matrix <- select(genotype_data, -V1)

# Convert to genind object (haploid)
genind_obj <- df2genind(snp_matrix,
                        ploidy = 1,
                        ind.names = individual_ids,
                        pop = population_labels,
                        sep = "\t",
                        strata = data.frame(population_labels))

# Save summary statistics of genind object
genind_summary <- summary(genind_obj)
capture.output(genind_summary,
               file = file.path(output_dir, "genind_summary.txt"))

# Convert genind to hierfstat format
hierfstat_obj <- genind2hierfstat(genind_obj)

# Compute basic diversity statistics (Ho, He, Fst)
basic_stats <- basic.stats(hierfstat_obj, diploid = FALSE)
write.csv(basic_stats$Ho, file = file.path(output_dir, "Ho_values.csv"))
write.csv(basic_stats$Hs, file = file.path(output_dir, "He_values.csv"))
write.csv(basic_stats$overall, file = file.path(output_dir, "overall_stats.csv"))

# Plot Ho and He distributions
pdf(file = file.path(output_dir, "He_Ho_distribution.pdf"))
par(mfrow = c(1, 2))
hist(basic_stats$Ho, main = "Observed Heterozygosity (Ho)", xlab = "Ho", col = "skyblue", breaks = 30)
hist(basic_stats$Hs, main = "Expected Heterozygosity (He)", xlab = "He", col = "salmon", breaks = 30)
dev.off()

# Compute global Fst (Weir and Cockerham)
fst_wc <- wc(hierfstat_obj, diploid = FALSE)
capture.output(fst_wc,
               file = file.path(output_dir, "fst_wc.txt"))

# Pairwise Fst estimates using different methods
pairwise_fst_wc84 <- genet.dist(hierfstat_obj, method = "WC84", diploid = FALSE)
write.csv(as.matrix(pairwise_fst_wc84),
          file = file.path(output_dir, "pairwise_fst_wc84.csv"))

pairwise_fst_nei87 <- genet.dist(hierfstat_obj, method = "Nei87", diploid = FALSE)
write.csv(as.matrix(pairwise_fst_nei87),
          file = file.path(output_dir, "pairwise_fst_nei87.csv"))

# AMOVA analysis ----
# Recreate genind object for clone-correction
genind_raw <- df2genind(snp_matrix,
                        ploidy = 1,
                        ind.names = individual_ids,
                        pop = population_labels,
                        sep = "\t")

# Convert to genclone object and apply clone correction
genclone_obj <- as.genclone(genind_raw, strata = population_data)
genclone_obj_cc <- clonecorrect(genclone_obj, strata = ~V2)

# Perform AMOVA
amova_result <- poppr.amova(genclone_obj_cc, ~V2, cutoff = 0.95)
capture.output(amova_result,
               file = file.path(output_dir, "amova_result.txt"))

# Permutation test for AMOVA significance
amova_test <- randtest(amova_result, nrepet = 1000)
capture.output(amova_test,
               file = file.path(output_dir, "amova_randtest.txt"))

# Plot permutation test result
pdf(file = file.path(output_dir, "amova_randtest_plot.pdf"))
plot(amova_test, main = "AMOVA Significance Test (Permutation)")
dev.off()

message("Genetic diversity analysis complete for 166 haplotypes. All outputs saved in: ", output_dir)