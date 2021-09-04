devtools::load_all() # if working in snpdiver folder
library(tidyverse)
library(bigsnpr)

workingdir <- file.path("~", "Github", "PanicleData", "9-DiversityPanel")
datadir <- file.path(workingdir, "data")
metadata <- readRDS(file.path(datadir, "metadata.rds"))
phe_v <- c("PAN_LEN", "PRIM_BN", "SEC_BN")

bigsnp_inputs <- read_delim(file.path(workingdir, "analysis", 
                                      "bigsnp_inputs.txt"),  delim = " ", 
                            col_names = "SNPfiles") 
subpop_v2 <- list(Gulf_and_Midwest = c("Gulf", "Midwest"), Midwest = "Midwest", 
                  Gulf = "Gulf", Atlantic = "Atlantic", 
                  Atlantic_and_Midwest = c("Atlantic", "Midwest"),
                  Atlantic_and_Gulf = c("Atlantic", "Gulf"), 
                  All = c("Atlantic", "Gulf", "Midwest", "4X"))

outputdir <- file.path(workingdir, "analysis", "mash")
if (!dir.exists(outputdir)) {
  dir.create(outputdir)
}

kinshipdir <- file.path(workingdir, "analysis", "heritability")
prefix1 <- "pvdiv_panicles_2019_BLUPs_kinship_geno_subset_"
prefix2 <- "pvdiv_panicles_2019_BLUPs_kinship_"

plantiddir <- file.path(workingdir, "analysis", "PLANT_ID_BLUPs")
prefix3 <- "pvdiv_panicles_2019_BLUPs_PLANT_ID_geno_subset_"
prefix4 <- "pvdiv_panicles_2019_BLUPs_PLANT_ID_"

phe_gwas1 <- read_csv(file = file.path(kinshipdir, 
                                       paste0(prefix1, "BLUP_phenotypes",
                                              ".csv"))) %>%
  rename(sample.ID = .data$PLANT_ID)
phe_gwas2 <- read_csv(file = file.path(kinshipdir, 
                                       paste0(prefix2, "BLUP_phenotypes",
                                              ".csv"))) %>%
  rename(sample.ID = .data$PLANT_ID)
phe_gwas3 <- read_csv(file = file.path(plantiddir, 
                                       paste0(prefix3, "BLUP_phenotypes",
                                              ".csv"))) %>%
  rename(sample.ID = .data$PLANT_ID)
phe_gwas4 <- read_csv(file = file.path(plantiddir, 
                                       paste0(prefix4, "BLUP_phenotypes",
                                              ".csv"))) %>%
  rename(sample.ID = .data$PLANT_ID)

phe_gwas_list <- list(PLANT_ID_All = phe_gwas4,
                      PLANT_ID_geno_subset = phe_gwas3,
                      kinship_geno_subset = phe_gwas1,
                      kinship_All = phe_gwas2)

k = 7
suffix_svd <- paste0("big_three_garden_SNPs", names(subpop_v2)[k])
pavir_svd <- readRDS(file = file.path(outputdir, paste0(suffix_svd,
                                                        "_svd.rds")))

for (i in seq_along(phe_gwas_list)) {
  pavir_pan <- phe_gwas_list[[i]]
  suffix = paste0(names(phe_gwas_list)[i])
  pavir_snp <- snp_attach(bigsnp_inputs$SNPfiles[k])
  
  m_pan <- dive_phe2mash(df = pavir_pan, snp = pavir_snp, svd = pavir_svd,
                         save.plots = FALSE,
                         type = "linear", suffix = suffix, 
                         outputdir = outputdir, num.strong = 5000, 
                         num.random = 100000)
  saveRDS(m_pan, file = file.path(outputdir, paste0("Full_mash_model_",
                                                    "5000", "_SNPs_U_ed_and_",
                                                    "100000", "_SNPs_mash_fit_",
                                                    suffix, ".rds")))
  mash_plot_covar(m_pan, saveoutput = TRUE, suffix = suffix)
  
  manhattan <- mash_plot_manhattan_by_condition(m_pan, snp = pavir_snp, 
                                                saveoutput = TRUE)
  mash_clumps <- snp_clumping(pavir_snp$genotypes, 
                              infos.chr = markers$CHRN, thr.r2 = 0.2, 
                              infos.pos = markers$POS, 
                              S = manhattan$ggman_df$log10BayesFactor)
  mash_df_clumped <- manhattan$ggman_df[mash_clumps,]
  write_csv(mash_df_clumped, file = file.path(outputdir,
                                              paste0("Clumped_mash_output_df_",
                                                     suffix, ".csv")))
}
