rm(list=ls())
library(BMTME)
library(tidyverse)
library(pegas)
library(hydroGOF)
setwd('c://Users/Li Zhang/Desktop/PanicleData/8-GenomicSelection/')

kin = readRDS(file ="Kinship_van_Raden_750g_fourway_cross_withInbred.rds")
genotype = read.loci('4way.qua',skip = 4)
genotype = genotype %>% dplyr::select(id) %>% rename(LINE=id) %>% 
  mutate_at(vars(LINE),list(as.character)) %>% mutate_at(vars(LINE),list(as.numeric))

phenotype = read.csv('PaniclePhenotypeData.csv')

phenotype = phenotype %>%dplyr::select(-CROSS) %>%  
  mutate_at(vars("LINE"), as.character) %>%
  mutate_at(vars("LINE"), as.numeric) 

phenos = NULL
for ( s in unique(phenotype$SITE)){
  tmp = genotype %>% left_join( phenotype %>% filter(SITE==s) )%>% 
    mutate(SITE=s) %>% arrange(LINE)
  phenos = rbind(phenos,tmp)
}
phenos = phenos %>% rename(ID=LINE, PL=LEN, PBN= PRM, SBN=SEC)

LG = cholesky(kin)
ZG = model.matrix(~0+as.factor(phenos$ID))
Z.G = ZG %*% LG
Z.E = model.matrix(~0 + as.factor(phenos$SITE))
ZEG = model.matrix(~0 + as.factor(phenos$ID):as.factor(phenos$SITE))
G2 = kronecker(diag(length(unique(phenos$SITE))), data.matrix(kin))
LG2 = cholesky(G2)
Z.EG = ZEG %*% LG2

Y = as.matrix(phenos[,-c(1,2)])
fm = BMTME(Y=Y, X=Z.E, Z1=Z.G,Z2=Z.EG, nIter= 1500, burnIn=1000, thin=2, bs=50)

load(("GS_panicle_4way.RData"))
load("GS_panicle_4way_WithoutSTIL.RData")

Observed = fm$Y
Predicted = fm$yHat
re1 = lm(Observed[,1]~Predicted[,1])
summary(re1)

re2 = lm(Observed[,2]~Predicted[,2])
summary(re2)

re3 = lm(Observed[,3]~Predicted[,3])
summary(re3)

Observed = Observed %>% as.data.frame()%>% gather('Traits','Observed')
Predicted = Predicted  %>% as.data.frame()%>% gather('Traits','Predicted')
Pre_Obs = Observed %>% bind_cols(Predicted)%>% select(-Traits1)

p1= ggplot((Pre_Obs %>% filter(Traits=='PL')), aes(Observed, Predicted))+geom_point()+
  theme_bw()+ggtitle('BMTME fitting for PL' )+
  geom_abline()+xlim(10, 80)+ylim(10,80)+
  annotate(geom="text",x=20,y= 60, label= bquote(R^2~'='~0.47))

p2= ggplot((Pre_Obs %>% filter(Traits=='PBN')), aes(Observed, Predicted))+geom_point()+
  theme_bw()+ggtitle('BMTME fitting for PBN' )+
  geom_abline()+xlim(10,50)+ylim(10,50)+
  annotate(geom="text",x=15,y= 40, label=bquote(R^2~'='~0.53))


p3= ggplot((Pre_Obs %>% filter(Traits=='SBN')), aes(Observed, Predicted))+geom_point()+
  theme_bw()+ggtitle('BMTME fitting for SBN' )+
  geom_abline()+xlim(0,25)+ylim(0,25)+
  annotate(geom="text",x=5,y= 20, label=bquote(R^2~'='~0.53))

##PLOTTING for the pearson correlation
x=summary(pm_BMORS)
write.csv(x,'PredictionAccuracy_4WAY_panicle.csv',row.names = F)
x$Environment= factor(x$Environment, 
                      levels=c('KING','PKLE','TMPL','OVTN',
                               'STIL','CLMB','MNHT','LINC','KBSM','BRKG'))
p4= ggplot((x%>%filter(Trait=='PL')), aes(Environment, Pearson))+geom_point()+
  ggtitle('BMORS Evaluation for PL')+
  geom_errorbar(aes(ymin=Pearson-SE_Pearson, ymax=Pearson+SE_Pearson), width=0.3)+
  theme_bw()+ylab("Pearson's Correlation")+
  theme(axis.text.x = element_text(size=7))+xlab("")

p5= ggplot((x%>%filter(Trait=='PBN')), aes(Environment, Pearson))+geom_point()+
  ggtitle('BMORS Evaluation for PBN')+
  geom_errorbar(aes(ymin=Pearson-SE_Pearson, ymax=Pearson+SE_Pearson), width=0.3)+
  theme_bw()+ylab("Pearson's Correlation")+
  theme(axis.text.x = element_text(size=7))+xlab("")

p6= ggplot((x%>%filter(Trait=='SBN')), aes(Environment, Pearson))+geom_point()+
  ggtitle('BMORS Evaluation for SBN')+
  geom_errorbar(aes(ymin=Pearson-SE_Pearson, ymax=Pearson+SE_Pearson), width=0.3)+
  theme_bw()+ylab("Pearson's Correlation")+
  theme(axis.text.x = element_text(size=7))+xlab("")


g=gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
ggsave(file='BMTME Model Fitting & Evaluation for Panicle_4way.tiff', g, width = 10.2, height = 5)

###########Genomic selection and evaluation on the three sites
load("GS_panicle_4way.RData")

Observed = read.csv('PanicleData_ModelEvluation.csv')
Observed = Observed %>% rename(ID=LINE,PL= LEN, PBN=PRM,SBN=SEC)
Predicted_all = fm$yHat
Predicted_all = phenos %>% select(ID, SITE)%>% 
  bind_cols(as.data.frame(Predicted_all))

Predicted= Predicted_all %>% filter(ID>=405) %>% 
  filter(SITE=='CLMB'|SITE=='KBSM'|SITE=='PKLE')

Observed = Observed %>% as.data.frame()%>% gather('Trait','Obs',-SITE,-ID)
Predicted = Predicted  %>% as.data.frame()%>% gather('Trait','Pred', -SITE,-ID)
Obs_Pre = Observed %>% bind_cols(Predicted)%>% select(-Trait1, -ID1, -SITE1)

Stats1 = Obs_Pre %>% group_by(Trait, SITE) %>% 
  summarize(Pbias =pbias(Pred, Obs),
            r=cor(Pred, Obs, use='complete.obs'))


ggplot(Pre_Obs,aes(Observed, Predicted, color=SITE))+geom_point()+facet_grid(SITE~Traits, scales = 'free')+
  theme_bw()+
  geom_text(data = eq_lm, 
            aes(x = 20, y = 20, label = eq_lm, family = "serif"), 
            size=6, parse = TRUE) +
  geom_text(data=eq_r,aes(x=20,y=30, label=paste0('r=', round(eq_r$r, 2)), family='serif'),
            size=6)+
  geom_abline()+
  theme(text = element_text(size=8))

theme_Li = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5),
  axis.title = element_text(size = 10), axis.text = element_text(size =8 ),
  axis.line.x = element_line(size = 0.35, colour = 'grey50'),
  axis.line.y = element_line(size = 0.35, colour = 'grey50'),
  axis.ticks = element_line(size = 0.25, colour = 'grey50'),
  strip.background = element_blank(),
  strip.text = element_text(hjust = 0.5, size = 14 ,vjust = 0)
)

Obs_Pre$SITE = factor(Obs_Pre$SITE,levels=c('KBSM','CLMB','PKLE'))

tmp = Obs_Pre %>% filter(Trait=='PL')
rs1 = Stats1 %>% filter(Trait=='PL')
myplots1 =
  ggplot(tmp,aes(Obs, Pred)) + geom_point()+
  facet_grid(SITE~., scales = 'free')+xlab('')+
  ylab('Predicted') + geom_abline()+theme_bw()+
  theme_Li+ theme(legend.position = 'none')+ggtitle("PL")+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank())+
  geom_text(data=rs1, aes(x=25, 60, label=paste0("r=",round(rs1$r, 2)),family = "serif"))+
  geom_text(data=rs1, aes(x=25, 55, label=paste0("%bias=",round(rs1$Pbias, 2),"%"),family = "serif"))+
  ylim(min(min(tmp$Obs, na.rm = T), min(tmp$Pred, na.rm = T)), 
       max(max(tmp$Obs, na.rm = T), max(tmp$Pred, na.rm = T)))+
  xlim(min(min(tmp$Obs, na.rm = T), min(tmp$Pred, na.rm = T)), 
       max(max(tmp$Obs, na.rm = T), max(tmp$Pred, na.rm = T)))


tmp = Obs_Pre %>% filter(Trait=='PBN')
rs2 = Stats1 %>% filter(Trait=='PBN')
myplots2 =
  ggplot(tmp,aes(Obs, Pred)) + geom_point()+
  facet_grid(SITE~., scales = 'free')+xlab('Observed')+
  ylab('') + geom_abline()+theme_bw()+
  theme_Li+ theme(legend.position = 'none')+ggtitle("PBN")+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank())+
  geom_text(data=rs2, aes(x=15, 45, label=paste0("r=",round(rs2$r, 2)),family = "serif"))+
  geom_text(data=rs2, aes(x=15, 40, label=paste0("%bias=",round(rs2$Pbias, 2),"%"),family = "serif"))+
  ylim(min(min(tmp$Obs, na.rm = T), min(tmp$Pred, na.rm = T)), 
       max(max(tmp$Obs, na.rm = T), max(tmp$Pred, na.rm = T)))+
  xlim(min(min(tmp$Obs, na.rm = T), min(tmp$Pred, na.rm = T)), 
       max(max(tmp$Obs, na.rm = T), max(tmp$Pred, na.rm = T)))

tmp = Obs_Pre %>% filter(Trait=='SBN')
rs3 = Stats1 %>% filter(Trait=='SBN')
myplots3 =
  ggplot(tmp,aes(Obs, Pred)) + geom_point()+
  facet_grid(SITE~., scales = 'free')+xlab('')+
  ylab('') + geom_abline()+theme_bw()+
  theme_Li+ theme(legend.position = 'none')+ggtitle("SBN")+
  geom_text(data=rs3, aes(x=7, 25, label=paste0("r=",round(rs3$r, 2)),family = "serif"))+
  geom_text(data=rs3, aes(x=7, 22, label=paste0("%bias=",round(rs3$Pbias, 2),"%"),family = "serif"))+
  ylim(min(min(tmp$Obs, na.rm = T), min(tmp$Pred, na.rm = T)), 
       max(max(tmp$Obs, na.rm = T), max(tmp$Pred, na.rm = T)))+
  xlim(min(min(tmp$Obs, na.rm = T), min(tmp$Pred, na.rm = T)), 
       max(max(tmp$Obs, na.rm = T), max(tmp$Pred, na.rm = T)))

F5 = cowplot::plot_grid(myplots1, myplots2, myplots3 , nrow=1)
ggsave(file='GenomicPrediction.tiff',F5, height =5, width = 8)




p1= ggplot((Pre_Obs %>% filter(Traits=='PL')), aes(Observed, Predicted, color=SITE))+geom_point()+
  theme_bw()+ggtitle('Genomic Prediction for PL' )+
  geom_abline()+xlim(10, 80)+ylim(10,80)+
  theme(text = element_text(size=8))+
  annotate(geom="text",x=20,y= 60, label= bquote(R^2~'='~0.36))+
  annotate(geom="text",x=20,y= 70, label= "r=0.61")

p2= ggplot((Pre_Obs %>% filter(Traits=='PBN')), aes(Observed, Predicted, color=SITE))+geom_point()+
  theme_bw()+ggtitle('Genomic Prediction for PBN' )+
  geom_abline()+xlim(10,50)+ylim(10,50)+
  theme(text = element_text(size=8))+
  annotate(geom="text",x=18,y= 40, label=bquote(R^2~'='~0.45))+
  annotate(geom="text",x=18,y= 45, label='r=0.67')


p3= ggplot((Pre_Obs %>% filter(Traits=='SBN')), aes(Observed, Predicted, color=SITE))+geom_point()+
  theme_bw()+ggtitle('Genomic Prediction for for SBN' )+
  geom_abline()+xlim(0,25)+ylim(0,25)+
  theme(text = element_text(size=8))+
  annotate(geom="text",x=4,y= 18, label=bquote(R^2~'='~0.54))+
  annotate(geom="text",x=4,y= 21, label='r=0.74')

g=gridExtra::grid.arrange(p1,p2,p3,ncol=3)
ggsave(file='BMTME Model Evaluation for Panicle_4way.tiff', g, width =10, height = 3)




