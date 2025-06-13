# Load required libraries
library(hierfstat)
library(adegenet)
library(dplyr)
library(poppr)
library(ggplot2)

# Working directory to script location
setwd("msc-codes-fmb_genetic_diversity/amova-Ho-He/scripts/")

# Output directory
output_dir <- "../results"
dir.create(output_dir, showWarnings = FALSE)

# Load input data - Used previously for structure analysis
population_data <- read.delim("../data/pop_file.txt", header = FALSE)
genotype_data <- read.delim("../data/174300_r_snps.txt", header = FALSE)

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
# # Recreate genind object for clone-correction
# genind_raw <- df2genind(snp_matrix,
#                         ploidy = 1,
#                         ind.names = individual_ids,
#                         pop = population_labels,
#                         sep = "\t")

# # Convert to genclone object and apply clone correction
# genclone_obj <- as.genclone(genind_raw, strata = population_data)
# genclone_obj_cc <- clonecorrect(genclone_obj, strata = ~V2)

# Add population labels as a data.frame with column name "Pop"
strata_info <- data.frame(Pop = population_labels)
rownames(strata_info) <- individual_ids  # important to match rownames

# Create genind object and assign strata
genind_raw <- df2genind(snp_matrix,
                        ploidy = 1,
                        ind.names = individual_ids,
                        pop = population_labels,
                        sep = "\t")

# Assign strata
strata(genind_raw) <- strata_info

# Convert to genclone and retain strata
genclone_obj <- as.genclone(genind_raw)

# Clone-correct using the correct strata formula
genclone_obj_cc <- clonecorrect(genclone_obj, strata = ~Pop)

# Perform AMOVA
amova_result <- poppr.amova(genclone_obj_cc, ~Pop, cutoff = 0.95)
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

message("Genetic diversity analysis complete. All outputs saved in: ", output_dir)
