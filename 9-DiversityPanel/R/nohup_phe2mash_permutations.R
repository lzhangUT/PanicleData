library(tidyverse)
library(bigsnpr)
library(snpdiver)

metadata <- readRDS("~/Github/PanicleData/9-DiversityPanel/data/metadata.rds")
subpops <- readRDS("~/Github/pvdiv-phenology-gxe/data/subpop_color_coding.rds")

ncores <- 10

effectspid <- big_attach("~/Github/PanicleData/9-DiversityPanel/analysis/gwas/PID_PERM/gwas_effects_all_panicle_phe_perm_pid.rds")
metapid <- read_csv("~/Github/PanicleData/9-DiversityPanel/analysis/gwas/PID_PERM/gwas_effects_all_panicle_phe_perm_pid_associated_metadata.csv")
snp05 <- snp_attach("~/Github/PanicleData/9-DiversityPanel/data/Pvirgatum_V5_GWAS_381g_PanicleData_subset_maf_0.05.rds")
all_phe <- readRDS(file = file.path("~/Github/PanicleData/9-DiversityPanel/",
                                    "data", "all_panicle_phe_perm_pid.rds"))

# Set up SNP subset given the PLANT_ID to include, the maf, and reduce to a 
# small set of SNPs to test
plants <- snp05$fam$sample.ID
panicle_381 <- which(plants %in% all_phe$PLANT_ID)
maf_thresh10 <-snp_MAF(snp05$genotypes, ind.row = panicle_381)
snp_subset_maf10 <- which(maf_thresh10 > 0.1)
snp_subset_maf10_10k <- c(1:floor(length(snp_subset_maf10)/1000))*1000
snp_subset <- snp_subset_maf10[snp_subset_maf10_10k]

# subset SNPs & effects given the above snp_subset of snp05
snp05_sub <- subset(snp05, ind.col = snp_subset)
snp_05s <- snp_attach(snp05_sub)
effects_sub <- big_copy(effectspid, ind.row = snp_subset)

# pick the phenotype subset
phe_rep <- c(10, 18, 22, 28, 30, 39, 44, 47, 51)  # PKLE A D CLMB KBSM A D C , D A E , E C B ,
# just use one phenotype rep from each condition

mtest <- dive_effects2mash(effects = effects_sub, snp = snp_05s, 
                           metadata = metapid, phe = phe_rep,
                           suffix = "test_REPs_maf10_10ksubset", 
                           outputdir = "~/Github/PanicleData/9-DiversityPanel/analysis/gwas/PID_PERM")

saveRDS(mtest, file = file.path("~/Github/PanicleData/9-DiversityPanel/analysis/gwas/", 
                                "PID_PERM/Mash_output_all_panicle_phe_perm_pid_REPs_mash_maf10_10ksubset.rds"))
