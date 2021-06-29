
rm(list=ls())
options(warn=-1)
library(readxl)
library(tidyverse)

###1. get the panicle data
setwd('C:/Users/Li Zhang/Desktop/PanicleData/0-RawDataAndGrandParentPlotting/')
getwd()

PanicleFiles=intersect(list.files(pattern='2016'),list.files(pattern='Panicle'))
Panicle = lapply(PanicleFiles, function(x){
  tmp = read_excel(x)
  tmp = tmp %>% dplyr::select(c(PLOT_GL, PAN, SAMP_ID, LEN, PRM, SEC))
})
Panicle = do.call(rbind, Panicle)

PlantList = read_excel('NSF_4WCR_Master Plant List_Final.xlsx')

##remove outliers--- after first try of boxplot
phenotype = PlantList %>% left_join(Panicle) %>% dplyr::select(SITE, LINE,CROSS, LEN, PRM, SEC) %>% filter(SITE!='BFLG') %>%
        mutate (CROSS=str_replace(CROSS,'VWAD', 'F2')) %>% mutate (CROSS=str_replace(CROSS,'ADVW', 'F2')) %>%
  filter(CROSS!='F1')%>% mutate_at(vars(matches('SITE|LINE|CROSS')), list(as.character)) %>%
  mutate_at(vars(matches('LEN|PRM|SEC')),list(as.numeric)) %>%
  mutate(LEN=replace(LEN, LEN>200, NA)) %>% mutate(PRM=replace(PRM, PRM>75, NA)) %>% mutate(SEC=replace(SEC, SEC>50, NA)) %>%
  group_by(SITE, LINE, CROSS) %>% summarise_all(list(mean),na.rm=T) %>% ungroup()

write.csv(phenotype,'PaniclePhenotypeData.csv', row.names = F)
file.copy('PaniclePhenotypeData.csv', '../1-GenstatFormatting/', overwrite = T)

###test the survivorship at each site using histogram
tiff("Histogram of Survivorship.tif",width = 900, height = 700,res=90 )
par(mfrow=c(3,4))
h = hist(phenotype$LEN, ylim = c(0, 1200), main='All Sites', xlab = "Panicle Length")
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))
sites = c("KING", "PKLE", "TMPL", "OVTN", "STIL", "CLMB", "MNHT", "LINC", "KBSM", "BRKG")
for (s in sites){
  tmp = phenotype %>% filter(SITE!=s)
  h = hist(tmp$LEN, ylim = c(0, 1200), main=paste0("Removing Site_", s), xlab = "Panicle Length")
  text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))
}

dev.off()



###genotype number smaller than 405 is common across the 10 sites
phenotype_all = phenotype %>% mutate_at(vars('LINE'), list(as.numeric)) %>% filter(LINE<405)%>% dplyr::rename(PL=LEN, PBN=PRM,SBN=SEC)


###plotting the phenotype data across the 10 sites,violin plot for F2 and points for F0
F0 = phenotype %>% filter(CROSS=='F0') %>% dplyr::select(-CROSS) %>% 
  mutate(LINE=str_replace(LINE, 'DAC6','DAC')) %>% mutate(LINE=str_replace(LINE, 'WBC3','WBC')) %>%
  gather (TRAIT, VALUE, -c(LINE, SITE)) %>% mutate(TRAIT=str_replace(TRAIT,'LEN','PL'))%>% 
  mutate(TRAIT=str_replace(TRAIT,'PRM','PBN'))%>% mutate(TRAIT=str_replace(TRAIT,'SEC','SBN'))
F2 = phenotype_all %>% filter(CROSS=='F2') %>% dplyr::select(-CROSS) %>% 
  gather (TRAIT, VALUE, -c(LINE, SITE)) %>% mutate(TRAIT=str_replace(TRAIT,'LEN','PL'))%>% 
  mutate(TRAIT=str_replace(TRAIT,'PRM','PBN'))%>% mutate(TRAIT=str_replace(TRAIT,'SEC','SBN'))
###Site information
sites = c("KING", "PKLE", "TMPL", "OVTN", "STIL", "CLMB", "MNHT", "LINC", "KBSM", "BRKG")
F0$SITE = factor(F0$SITE,levels= sites)
F2$SITE = factor(F2$SITE,levels= sites)
traits = c('PL','PBN','SBN')
F0$TRAIT = factor(F0$TRAIT, levels = traits)
F2$TRAIT = factor(F2$TRAIT, levels = traits)

pdf('ViolinPlot.pdf',height = 8)
ggplot(F2,aes(SITE, VALUE)) + geom_violin(color="black",fill='lightgrey') + xlab("")+theme_bw()+
  facet_wrap(TRAIT~., scales = 'free',strip.position = "left",nrow = 3) + ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 10, angle = 90))+
   geom_jitter(data=F0, aes(SITE, VALUE, shape=LINE), color='black', width = 0.05,size=2) + 
  theme(axis.title=element_text(size=14,face="bold"))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title=element_blank())
dev.off()

pdf('PhenotypicCorrelation.pdf')
M = cor(phenotype_all[,4:6], use = 'complete.obs')
 corrplot::corrplot.mixed(M,lower.col = 'blue')
dev.off()

PhenotypicCorrelation = NULL
for (s in unique(phenotype_all$SITE)){
  tmp = phenotype_all %>% filter(SITE==s)
  M = cor(tmp[,4:6],use = 'complete.obs')
  M1 = cbind(M,s)
  PhenotypicCorrelation = rbind(PhenotypicCorrelation, M1)
}
write.csv(PhenotypicCorrelation,'PhenotypicCorrelationAtEachSiteBetweenTraits.csv',row.names = F)



