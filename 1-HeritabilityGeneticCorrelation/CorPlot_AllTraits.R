setwd('c:/Users/Li Zhang/Desktop/PanicleData/1-HeritabilityGeneticCorrelation/')
pdf('PairwiseGeneticCorrelationBetweenSitesForPL_PBN_SBN.pdf',onefile=T)
M = read.csv('GeneticCorrelation_between_site_PL.csv', row.names = 1)
M= as.matrix(M)
corrplot::corrplot.mixed(M,lower.col = 'blue')

M = read.csv('GeneticCorrelation_between_site_PBN.csv', row.names = 1)
M= as.matrix(M)
corrplot::corrplot.mixed(M,lower.col = 'blue')

M = read.csv('GeneticCorrelation_between_site_SBN.csv', row.names = 1)
M= as.matrix(M)
corrplot::corrplot.mixed(M,lower.col = 'blue',use = "pairwise.complete.obs")
dev.off()
