# Genetic diversity analysis

#"hierfstat"----
library("hierfstat")
library("adegenet")
library("dplyr")

setwd('msc-codes-fmb_genetic_diversity/amova-Ho-He/scripts/')

# My data
pop_file <- read.delim("../data/pop_file.txt", header = F) # population lables

df <- read.delim("../data/174300_r_snps.txt", header = F) # read the data file

loci <- select(df, -c(V1)) # exclude the first column (lables) so that we are able to convert the data to a matri
isolates <- as.character(pop_file$V1) # individuals labels
all_coutries <- as.character(pop_file$V2) # country labels

df_2 <- df2genind(loci, ploidy = 1, ind.names = isolates, pop = all_coutries, sep = "\t", strata = data.frame(all_coutries))

summary(df_2) # save output

df_3 <- genind2hierfstat(df_2)
stats_hier <- basic.stats(df_3, diploid = F) # Ho and He analyses
stats_hier
# basic.stats(df_3, diploid =F) # Ho and He analyses
wc(df_3, diploid = F) # Fst

# Pairwise Fst
wc84 <- genet.dist(df_3, method = "WC84", diploid = F)
wc84
nei87 <- genet.dist(df_3, method = "Nei87", diploid = F)
nei87

# 166 haplotypes
# haplotypes_166 <- c(
#   "E1", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E2", "E20", "E21", "E22", "E25", "E27", "E28", "E29", "E3",
#   "E30", "E31", "E33", "E37", "E38", "E39", "E4", "E40", "E41", "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49",
#   "E50", "E51", "E52", "E53", "E54", "E57", "E58", "E61", "E63",
#   "K1", "K13", "K14", "K15", "K16", "K17", "K18", "K20", "K21", "K23", "K24", "K26", "K27", "K28", "K29", "K31", "K32",
#   "K33", "K34", "K35", "K37", "K38", "K39", "K4", "K40", "K41", "K42", "K44", "K6", "K7", "K8", "K9",
#   "T11", "T12", "T15", "T16", "T17", "T18", "T19", "T2", "T20", "T22", "T24", "T25", "T26", "T29", "T3", "T30", "T31",
#   "T33", "T34", "T35", "T37", "T38", "T4", "T40", "T42", "T43", "T45", "T46", "T47", "T48", "T49", "T5", "T50", "T51",
#   "T53", "T54", "T55", "T56", "T57", "T58", "T6", "T7", "T8", "T9",
#   "U1", "U13", "U14", "U15", "U2", "U21", "U22", "U24", "U25", "U26", "U27", "U3", "U31", "U32", "U34", "U35", "U37",
#   "U38", "U4", "U40", "U42", "U44", "U45", "U46", "U5", "U51", "U56", "U57", "U8",
#   "UK10", "UK11", "UK12", "UK13", "UK14", "UK15", "UK2", "UK21", "UK22", "UK24", "UK25", "UK3", "UK4", "UK5", "UK6",
#   "UK7", "UK9"
# )

# df_2 <- df2genind(loci, ploidy = 1, ind.names = isolates, pop = all_coutries, sep = "\t") # need these for the optimized code to add stratification by population


#strata(df_2) <- data.frame(pop(df_2)) # These are initial approches optimized code is below

#df_2 <- df2genind(loci, ploidy = 1, ind.names = isolates, pop = all_coutries, sep = "\t", strata = data.frame(pop(df_2)))
# df_2 <- df2genind(loci, ploidy = 1, ind.names = isolates, pop = all_coutries, sep = "\t")
# addStrata(df)

#div <- summary(df_2) #Get Genetic diversity (observed and expected heterozygosity)
#names(div)
#plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", main="Observed heterozygosity per locus")
#plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", main="Expected heterozygosity as a function of observed heterozygosity per locus")

#"poppr - AMOVA test"----
library(poppr)

df_2 <- df2genind(loci, ploidy = 1, ind.names = isolates, pop = all_coutries, sep = "\t")
# strata(df_2) <- data.frame(pop(df_2))##

df_2_genclone = as.genclone(df_2)
pop_file <- read.delim("../data/pop_file.txt", header = F)

df_2_genclone_strat <- as.genclone(df_2_genclone, strata=pop_file)

## AMOVA TESTS (not including outgroup individuals)
amova_results <- poppr.amova(df_2_genclone_strat, pop, cutoff=0.95)#, clonecorrect = TRUE)
amova_results
# test for significance 
amova_test <- randtest(amova_results, nrepet = 1000)
amova_test
plot(amova_test)

