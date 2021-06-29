
rm(list=ls())
options(warn=-1)
library(tidyverse)
library(pegas)
library(sommer)

setwd('C:/Users/Li Zhang/Desktop/PanicleData/1-HeritabilityGeneticCorrelation/')
getwd()

###Formatting map file
genotype = read.loci('4way.qua',skip = 4)
genotype = genotype %>% dplyr::select(id) %>% rename(LINE=id) %>% 
  mutate_at(vars(LINE),list(as.character)) %>% mutate_at(vars(LINE),list(as.numeric))

genomap = read.delim('4way_male_reduced_marker.map',header = F, sep = '') #read raw map 
CPmap = genomap %>% dplyr::filter(!str_detect(V1, 'group')) %>% separate(V1,c('Chr','Locus'), sep = '_') %>% 
  rowwise %>% dplyr::mutate(Chrom=substr(Chr,5,6)) %>% dplyr::select(Locus,V2,Chrom) %>% dplyr::rename(Position=V2)
  
###read the loci file and remove the genotype information for these lines which had no phenotyping data
geno=read.loci('4way_reduced_marker.loc',header = F)
geno2 = t(geno)
geno2 = data.frame(geno2[-c(1:3),])
colnames(geno2) = CPmap$Locus
CPgeno1 = cbind(LINE=genotype$LINE,geno2)

CPgeno2 = c()
for (i in 1:(ncol(CPgeno1)-1)){
  df = CPgeno1[,c(1,i+1)]
  df1 = cbind(df,df[,2],df[,2],df[,2])
  df1[,2:5] = lapply(df1[,2:5], as.character)
  colnames(df1)[2:5] = c(rep(colnames(df1)[2],4))
  
  for (j in 1:nrow(df1)){
    if(df1[j,2]=="ac" |df1[j,2]=="ad"){df1[j,2] = 1}else{df1[j,2] = 0}
    if(df1[j,3]=="bc" |df1[j,3]=="bd"){df1[j,3] = 1}else{df1[j,3] = 0}
    if(df1[j,4]=="ac" |df1[j,4]=="bc"){df1[j,4] = 1}else{df1[j,4] = 0}
    if(df1[j,5]=="ad" |df1[j,5]=="bd"){df1[j,5] = 1}else{df1[j,5] = 0}
  }
  X = df1[,-1]
  X = data.matrix(X, rownames.force = NA)
  qr = qr(X)
  y = X[, qr$pivot[1:qr$rank]]
  CPgeno2 = cbind(CPgeno2,y)
}

rownames(CPgeno2) = genotype$LINE
CPgeno = data.matrix(CPgeno2)

A = A.mat(CPgeno) # additive relationship matrix 
#D = D.mat(CPgeno) # dominance relationship matrix
#E = E.mat(CPgeno) # epistatic relationship matrix
save.image('DataFilesNeededForSommer.Rdata')
###################################################
load("DataFilesNeededForSommer.Rdata")
phen = read.csv('PaniclePhenotypeData.csv')
pheno2 = phen %>% mutate_at(vars(LINE), list(as.character)) %>% mutate_at(vars(LINE), list(as.numeric)) %>% 
   dplyr::select(-CROSS) %>% filter(LINE<405) 

out = vector('list', 0)
for (s in unique(pheno2$SITE)){
  tmp = pheno2 %>% filter(SITE==s)
  out[[s]] = tmp %>% right_join(genotype)%>% mutate(SITE = s)
}
out =bind_rows(out)
CPpheno = out %>% mutate_at(vars(LINE), list(as.factor))
  
## calculate the heritability and variance component for each trait at each location
pheid = colnames(CPpheno)[-c(1:2)]
heritability = c()
varComp = c()

for (s in unique(CPpheno$SITE)){
  df1 = CPpheno %>% dplyr::filter(SITE==s) 
  for (phe in pheid){
    if (sum(df1[,phe],na.rm = T)==0) {h2 = data.frame(SITE=s,TRAIT=phe,h2=NA,h2_SE=NA)
    Vave = cbind(data.frame(SITE=rep(s,2),TRAIT=rep(phe,2),VARIANCE=c('Va','Ve'), 
                            VarComp=rep(NA,2),VarCompSE=rep(NA,2),Zratio=rep(NA,2)))}else{
                              form =as.formula(paste(phe," ~ 1"))
                              A_phe = mmer(form, random=~vs(LINE, Gu=A), rcov=~vs(us(SITE),units),data=df1, verbose = F)
                              VaVe = cbind(data.frame(SITE=rep(s,2), TRAIT=rep(phe,2),VARIANCE=c('Va','Ve')),data.frame(summary(A_phe)$varcomp))
                              ###if interested in the random effect ('blups') as breeding value for each line:
                              Eff=data.frame(Eff=unlist(A_phe$U$`u:LINE`)+A_phe$sigma$`u:LINE`)
                              ##the SE for random effect:
                              SE=sqrt(diag(A_phe$PevU$`u:LINE`[[1]]))  
                              #################
                              h2 = data.frame(c(s, phe,round(pin(A_phe, h2 ~ V1 / ( V1 + V2 )),2)))
                              colnames(h2) = c('SITE','TRAIT','h2','h2_SE')}
    heritability = rbind(heritability,h2)
    varComp = rbind(varComp,VaVe)}
}
###plot heritability
sites = c("KING", "PKLE", "TMPL", "OVTN", "STIL", "CLMB", "MNHT", "LINC", "KBSM", "BRKG")
heritability$SITE = factor(heritability$SITE, levels = sites)
ggplot(heritability, aes(SITE,h2)) + geom_col() + facet_grid(TRAIT~.)+ ggtitle('h2 at each site')

###plot variance components
vc = varComp
vc$SITE = factor(vc$SITE,levels=sites)
vc_va = subset(vc, VARIANCE=='Va')
vc_ve = subset(vc, VARIANCE=='Ve')
vc_ve = vc_ve %>% dplyr::mutate(VarComp=-VarComp)

ggplot(vc_va,aes(SITE,VarComp,fill=VARIANCE))+geom_col() + facet_grid(TRAIT~.,scales = 'free')+
  geom_hline(yintercept=0, color = "black",lwd=1.2)+ylab('')+xlab('')+
  geom_col(data=vc_ve,aes(SITE,VarComp,fill=VARIANCE))+facet_grid(TRAIT~.,scales = 'free')+theme_bw()
ggsave('Va_Ve.png')
######

###genetic correlation among trait for each site
Rg <- c()
for (i in unique(CPpheno$SITE)){
  df <- subset(CPpheno,SITE==i)
  print(i)
  rg1 = mmer(cbind(LEN, PRM, SEC)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df) 
  x1 = as.data.frame(cov2cor(rg1$sigma$`u:LINE`)) %>% rownames_to_column('trait') %>% dplyr::mutate(SITE=i)
  Rg = rbind(Rg,x1)
}

write.csv(Rg,'GeneticCorrelation_each_site.csv',row.names=F)
###genetic correlation across sites
rg1 = mmer(cbind(LEN, PRM, SEC)~1, random=~vs(LINE, Gu=A) + SITE, rcov=~units,data=CPpheno) 
rg = cov2cor(rg1$sigma$`u:LINE`)

save.image('heritability_genetic_correlation.RData')


###pairwised genetic correlation between site for one trait at a time
CPpheno2 = CPpheno %>% gather('traits','values',-LINE,-SITE)%>% spread(SITE, values)
Rg2 <- c()
for (i in unique(CPpheno2$traits)){
  df <- subset(CPpheno2,traits==i)
  print(i)
  rg1 = mmer(cbind(BRKG, KBSM)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg2 = mmer(cbind(BRKG, LINC)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg3 = mmer(cbind(BRKG, MNHT)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg4 = mmer(cbind(BRKG, CLMB)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg5 = mmer(cbind(BRKG, STIL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg6 = mmer(cbind(BRKG, OVTN)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg7 = mmer(cbind(BRKG, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg8 = mmer(cbind(BRKG, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg9 = mmer(cbind(BRKG, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg10 = mmer(cbind(KBSM, LINC)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg11 = mmer(cbind(KBSM, MNHT)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg12 = mmer(cbind(KBSM, CLMB)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg13 = mmer(cbind(KBSM, STIL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg14 = mmer(cbind(KBSM, OVTN)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg15 = mmer(cbind(KBSM, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg16 = mmer(cbind(KBSM, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg17 = mmer(cbind(KBSM, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg18 = mmer(cbind(LINC, MNHT)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg19 = mmer(cbind(LINC, CLMB)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg20 = mmer(cbind(LINC, STIL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg21 = mmer(cbind(LINC, OVTN)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg22 = mmer(cbind(LINC, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg23 = mmer(cbind(LINC, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg24 = mmer(cbind(LINC, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg25 = mmer(cbind(MNHT, CLMB)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg26 = mmer(cbind(MNHT, STIL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg27 = mmer(cbind(MNHT, OVTN)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg28 = mmer(cbind(MNHT, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg29 = mmer(cbind(MNHT, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg30 = mmer(cbind(MNHT, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg31 = mmer(cbind(CLMB, STIL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg32 = mmer(cbind(CLMB, OVTN)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg33 = mmer(cbind(CLMB, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg34 = mmer(cbind(CLMB, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg35 = mmer(cbind(CLMB, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg36 = mmer(cbind(STIL, OVTN)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg37 = mmer(cbind(STIL, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg38 = mmer(cbind(STIL, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg39 = mmer(cbind(STIL, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg40 = mmer(cbind(OVTN, TMPL)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg41 = mmer(cbind(OVTN, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg42 = mmer(cbind(OVTN, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg43 = mmer(cbind(TMPL, PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg44 = mmer(cbind(TMPL, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  rg45 = mmer(cbind(PKLE, KING)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df)
  
  x1 = as.data.frame(cov2cor(rg1$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x1)[2:3]=c('SITE1','SITE2')
  x2 = as.data.frame(cov2cor(rg2$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x2)[2:3]=c('SITE1','SITE2')
  x3 = as.data.frame(cov2cor(rg3$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x3)[2:3]=c('SITE1','SITE2')
  x4 = as.data.frame(cov2cor(rg4$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x4)[2:3]=c('SITE1','SITE2')
  x5 = as.data.frame(cov2cor(rg5$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x5)[2:3]=c('SITE1','SITE2')
  x6 = as.data.frame(cov2cor(rg6$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x6)[2:3]=c('SITE1','SITE2')
  x7 = as.data.frame(cov2cor(rg7$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x7)[2:3]=c('SITE1','SITE2')
  x8 = as.data.frame(cov2cor(rg8$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x8)[2:3]=c('SITE1','SITE2')
  x9 = as.data.frame(cov2cor(rg9$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x9)[2:3]=c('SITE1','SITE2')
  x10 = as.data.frame(cov2cor(rg10$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x10)[2:3]=c('SITE1','SITE2')
  
  x11 = as.data.frame(cov2cor(rg11$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x11)[2:3]=c('SITE1','SITE2')
  x12 = as.data.frame(cov2cor(rg12$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x12)[2:3]=c('SITE1','SITE2')
  x13 = as.data.frame(cov2cor(rg13$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x13)[2:3]=c('SITE1','SITE2')
  x14 = as.data.frame(cov2cor(rg14$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x14)[2:3]=c('SITE1','SITE2')
  x15 = as.data.frame(cov2cor(rg15$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x15)[2:3]=c('SITE1','SITE2')
  x16 = as.data.frame(cov2cor(rg16$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x16)[2:3]=c('SITE1','SITE2')
  x17 = as.data.frame(cov2cor(rg17$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x17)[2:3]=c('SITE1','SITE2')
  x18 = as.data.frame(cov2cor(rg18$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x18)[2:3]=c('SITE1','SITE2')
  x19 = as.data.frame(cov2cor(rg19$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x19)[2:3]=c('SITE1','SITE2')
  x20 = as.data.frame(cov2cor(rg20$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x20)[2:3]=c('SITE1','SITE2')
  
  x21 = as.data.frame(cov2cor(rg21$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x21)[2:3]=c('SITE1','SITE2')
  x22 = as.data.frame(cov2cor(rg22$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x22)[2:3]=c('SITE1','SITE2')
  x23 = as.data.frame(cov2cor(rg23$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x23)[2:3]=c('SITE1','SITE2')
  x24 = as.data.frame(cov2cor(rg24$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x24)[2:3]=c('SITE1','SITE2')
  x25 = as.data.frame(cov2cor(rg25$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x25)[2:3]=c('SITE1','SITE2')
  x26 = as.data.frame(cov2cor(rg26$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x26)[2:3]=c('SITE1','SITE2')
  x27 = as.data.frame(cov2cor(rg27$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x27)[2:3]=c('SITE1','SITE2')
  x28 = as.data.frame(cov2cor(rg28$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x28)[2:3]=c('SITE1','SITE2')
  x29 = as.data.frame(cov2cor(rg29$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x29)[2:3]=c('SITE1','SITE2')
  x30 = as.data.frame(cov2cor(rg30$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x30)[2:3]=c('SITE1','SITE2')
  
  x31 = as.data.frame(cov2cor(rg31$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x31)[2:3]=c('SITE1','SITE2')
  x32 = as.data.frame(cov2cor(rg32$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x32)[2:3]=c('SITE1','SITE2')
  x33 = as.data.frame(cov2cor(rg33$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x33)[2:3]=c('SITE1','SITE2')
  x34 = as.data.frame(cov2cor(rg34$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x34)[2:3]=c('SITE1','SITE2')
  x35 = as.data.frame(cov2cor(rg35$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x35)[2:3]=c('SITE1','SITE2')
  x36 = as.data.frame(cov2cor(rg36$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x36)[2:3]=c('SITE1','SITE2')
  x37 = as.data.frame(cov2cor(rg37$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x37)[2:3]=c('SITE1','SITE2')
  x38 = as.data.frame(cov2cor(rg38$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x38)[2:3]=c('SITE1','SITE2')
  x39 = as.data.frame(cov2cor(rg39$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x39)[2:3]=c('SITE1','SITE2')
  x40 = as.data.frame(cov2cor(rg40$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x40)[2:3]=c('SITE1','SITE2')
  
  x41 = as.data.frame(cov2cor(rg41$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x41)[2:3]=c('SITE1','SITE2')
  x42 = as.data.frame(cov2cor(rg42$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x42)[2:3]=c('SITE1','SITE2')
  x43 = as.data.frame(cov2cor(rg43$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x43)[2:3]=c('SITE1','SITE2')
  x44 = as.data.frame(cov2cor(rg44$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x44)[2:3]=c('SITE1','SITE2')
  x45 = as.data.frame(cov2cor(rg45$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(traits=i)
  colnames(x45)[2:3]=c('SITE1','SITE2')
  
  x= rbind(x1, x2,x3,x4,x5,x6,x7,x8,x9,x10,
           x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,
           x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,
           x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,
           x41,x42,x43,x44,x45)
  Rg2 = rbind(Rg2,x)
}

write.csv(Rg2,'GeneticCorrelation_between_site.csv',row.names=F)

########genetic correlation between panicle traits and fl50, tc, biomass
FL50 = read.csv('FL50_BIO_TC_2016.csv')
CPpheno3 = CPpheno %>% left_join(FL50 %>% mutate_at(vars(LINE), as.character) %>% mutate_at(vars(LINE),as.factor))

Rg3 <- c()
for (i in unique(CPpheno3$SITE)){
  df <- subset(CPpheno3,SITE==i)
  print(i)
  rg1 = mmer(cbind(LEN, PRM, SEC, FL50, TC, BIO)~1, random=~vs(LINE, Gu=A), rcov=~units,data=df) 
  x1 = as.data.frame(cov2cor(rg1$sigma$`u:LINE`)) %>% rownames_to_column('trait') %>% dplyr::mutate(SITE=i)
  Rg3 = rbind(Rg3,x1)
}

write.csv(Rg3,'GeneticCorrelationBetweenPanicleTraitsFL50_each_site.csv',row.names=F)

save.image(file = 'GeneticCorrelation_All.RData')




###
heritability = heritability %>% mutate(TRAIT=str_replace(TRAIT,'LEN','PL'))%>% 
  mutate(TRAIT=str_replace(TRAIT,'PRM','PBN'))%>% mutate(TRAIT=str_replace(TRAIT,'SEC','SBN'))
heritability$TRAIT = factor(heritability$TRAIT, levels=c('PL','PBN','SBN'))

load('heritability_genetic_correlation.RData')
pdf('Heritability.pdf')
ggplot(heritability, aes(SITE,h2)) + geom_col(alpha=0.6) + facet_grid(TRAIT~.)+theme_bw()+ ylab('Heritability')+xlab('')+
  geom_errorbar(aes(ymin=h2-h2_SE, ymax=h2+h2_SE), width=0.3)+
  theme(text = element_text(face = "bold", size = 12))+
  theme(axis.text = element_text(face = "bold", size = 12),axis.title = element_text(face = "bold", size = 14))+
  theme(strip.text = element_text(face="bold", size=12))
dev.off() 
 
pdf('GeneticCorrelation.pdf')
corrplot::corrplot.mixed(rg, lower.col = 'blue')
dev.off()


# write.csv(rg,'GeneticCorrelationAmongTraitsAcrossSites.csv')


