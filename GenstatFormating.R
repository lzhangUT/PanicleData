
rm(list=ls())
#options(warn=-1)
library(tidyverse)
library(pegas)

### get the panicle data, phenotype and genotype
setwd('C:/Users/Li Zhang/Desktop/PanicleData/1-GenstatFormatting/')
getwd()

phenotype = read.csv('PaniclePhenotypeData.csv')

genotype = read.loci('4way.qua',skip = 4)
genotype = genotype %>% dplyr::select(id) %>% rename(LINE=id) %>% 
  mutate_at(vars(LINE),list(as.character)) %>% mutate_at(vars(LINE),list(as.numeric))

###seperate the data for modeling building and model evaluation.
phenotype_ModelBuilding = phenotype %>% mutate_at(vars('LINE'), list(as.character))%>% 
  mutate_at(vars('LINE'), list(as.numeric)) %>% filter(LINE<405)

counts_overlap = phenotype_ModelBuilding %>% count(LINE)%>% count(n)
ggplot(counts_overlap, aes(n, nn))+geom_col()+
  scale_x_reverse(breaks=seq(10,1,-1),expand=c(0,0))+
  scale_y_continuous(limits = c(0,215),expand = c(0,0))+
  theme_bw()+xlab("Number of Field Sites")+
  ylab("Number of Overlapping Genotypes")+
  geom_text(aes(label=nn), vjust=-0.5, size=3)+
  theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title = element_text(size = 10), axis.text = element_text(size =8),
  axis.line.x = element_line(size = 0.35, colour = 'grey50'),
  axis.line.y = element_line(size = 0.35, colour = 'grey50'),
  axis.ticks = element_line(size = 0.25, colour = 'grey50')
)
ggsave('OverlappingGenotypesAcrossSites.tiff')

### Formatting for Genstat  
out_ModelBuilding = vector('list', 0)
for (s in unique(phenotype_ModelBuilding$SITE)){
  tmp = phenotype_ModelBuilding %>% dplyr::select(-CROSS) %>% filter(SITE==s)
  out_ModelBuilding[[s]] = tmp %>% right_join(genotype)%>% mutate(SITE = s)
}
out_ModelBuilding =bind_rows(out_ModelBuilding)

LineRemove_ModelBuilding = out_ModelBuilding %>% group_by(LINE) %>% summarise(SUM = sum(PRM,na.rm=T)) %>% filter(SUM==0) 
Genstat_ModelBuilding = filter(out_ModelBuilding,!LINE %in% LineRemove_ModelBuilding$LINE)

Genstat1 = Genstat_ModelBuilding %>% dplyr::mutate(SITE2=SITE)
Genstat1[is.na(Genstat1)] = "*"
write.csv(Genstat1,'PanicleData_ModelBuilding.csv',row.names = F)

###read the loci file and remove the genotype information for these lines which had no phenotyping data
geno=read.loci('4way_reduced_marker.loc',header = F)
colnames(geno) = c('Marker','Phase','NoIdea',as.character(genotype$LINE))
geno_model_building = geno[,!colnames(geno) %in% LineRemove_ModelBuilding$LINE]
write.loci(geno_model_building,'geno_model_building.loc',loci.sep = '\t',quote=F,row.name=F,col.names=F)

Genotype_ModelBuilding= data.frame(line=as.numeric(colnames(geno_model_building)[-c(1:3)]))
write.csv(Genotype_ModelBuilding,'Genotype_ModelBuilding.csv',row.names = F)

### put the loci and genotype_sub file together for Genstat format
y = readLines('geno_model_building.loc')
z = readLines('Genotype_ModelBuilding.csv')

cat('name=cross',file='4way_genotype_ModelBuilding.loc',sep = '\n',append = F)
cat('popt=CP',file='4way_genotype_ModelBuilding.loc',sep = '\n',append = T)
cat(paste0('nloc=', nrow(geno)),file='4way_genotype_ModelBuilding.loc',sep = '\n',append = T)
cat(paste0('nind = ',length(Genotype_ModelBuilding$line)),file='4way_genotype_ModelBuilding.loc',sep = '\n',append = T)
cat(y,file='4way_genotype_ModelBuilding.loc',sep = '\n',append = T)
cat('individual names:',file='4way_genotype_ModelBuilding.loc',sep = '\n',append = T)
cat(z[-1],file='4way_genotype_ModelBuilding.loc',sep = '\n',append = T)

file.remove(c('geno_model_building.loc','Genotype_ModelBuilding.csv'))
file.copy(c('4way_genotype_ModelBuilding.loc','4way_male_reduced_marker.map','PanicleData_ModelBuilding.csv'),"c:/Users/Li Zhang/Desktop/PanicleData/2-GenstatRunning/",overwrite = T)

#########subset data for model evaluation
phenotype_ModelEvaluation = phenotype %>% mutate_at(vars('LINE'), list(as.character))%>%
  mutate_at(vars('LINE'), list(as.numeric)) %>% 
  filter(LINE >= 405) %>% filter(SITE=='CLMB'|SITE=='KBSM'|SITE=='PKLE')

out_ModelEvaluation = vector('list', 0)
for (s in unique(phenotype_ModelEvaluation$SITE)){
  tmp = phenotype_ModelEvaluation %>% dplyr::select(-CROSS) %>% filter(SITE==s)
  out_ModelEvaluation[[s]] = tmp %>% right_join(genotype)%>% mutate(SITE = s)
}
out_ModelEvaluation =bind_rows(out_ModelEvaluation)

LineRemove_ModelEvaluation = out_ModelEvaluation %>% group_by(LINE) %>% summarise(SUM = sum(PRM,na.rm=T)) %>% filter(SUM==0) 
Genstat_ModelEvaluation = filter(out_ModelEvaluation,!LINE %in% LineRemove_ModelEvaluation$LINE)
write.csv(Genstat_ModelEvaluation ,'PanicleData_ModelEvluation.csv',row.names = F)

#####################genotype file for the evaluation dataset
geno_model_evaluation = geno[,!colnames(geno) %in% LineRemove_ModelEvaluation$LINE]
write.loci(geno_model_evaluation,'4way_genotype_ModelEvaluation.loc',loci.sep = '\t',quote=F,row.name=F,col.names=F)

file.copy(c('4way_genotype_ModelEvaluation.loc','PanicleData_ModelEvluation.csv'),"c:/Users/Li Zhang/Desktop/PanicleData/4-GenstatModelEvaluation/",overwrite = T)
###End
#############Extra of FL50 and HT
FL50_HT_RAW = read.csv('FL50_HT.csv')

FL50_HT_RAW = FL50_HT_RAW %>% mutate_at(vars('LINE'), list(as.character))%>% mutate_at(vars('LINE'), list(as.numeric)) %>% filter(LINE<405)
out = vector('list', 0)
for (s in unique(FL50_HT_RAW$SITE)){
  tmp = FL50_HT_RAW  %>% mutate_at(vars(BIOMASS),as.character) %>%mutate_at(vars(BIOMASS),as.numeric) %>%
    dplyr::select(-YEAR) %>%  filter(SITE==s) %>% group_by(SITE, LINE)%>%
    summarise(FL50 = mean(FL50,na.rm = T), TC=mean(TC, na.rm = T), BIO=mean(BIOMASS, na.rm=T)) %>% ungroup()
  out[[s]] = tmp %>%  right_join(genotype)%>% mutate(SITE = s)
}
out =bind_rows(out)

LineRemove_out = out%>% group_by(LINE) %>% summarise(SUM = sum(TC,na.rm=T)) %>% filter(SUM==0) 
Genstat_out = filter(out,!LINE %in% LineRemove_out$LINE)

Genstat_FL50 = Genstat_out %>% dplyr::mutate(SITE2=SITE)
Genstat_FL50[is.na(Genstat_FL50)] = "*"
write.csv(Genstat_FL50,'FL50_BIO_TC_2016_Genstat.csv',row.names = F)

###USING geno_model_building.loc and FL50_HT_2016_Genstat.csv for Genstat running.