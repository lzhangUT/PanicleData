rm(list=ls())
library(tidyverse)
library(pegas)
library(xlsx)
library(hydroGOF)

setwd('c://Users/Li Zhang/Desktop/PanicleData/4-GenstatModelEvaluation/')

geno_validation = read.loci('4way_genotype_ModelEvaluation.loc',header = F)
Observed = read.csv('PanicleData_ModelEvluation.csv')
ValidationLine = data.frame(LINE=unique(Observed$LINE))

traits =c('PL','PBN','SBN')
constants = c(46.77, 30.25,12.76)

Stats = NULL
Obs_Pre_All=NULL
i=0

for (t in traits){
  i=i+1
  Eff = read.xlsx(paste0('EFF_',t,'.xlsx'), sheetIndex = 1)
  geno_QTLs = geno_validation %>% slice(Eff$POS) %>% dplyr::select(-V2, -V3) %>% mutate(POS = Eff$POS) 
  colnames(geno_QTLs) = c('MarkerName', unique(ValidationLine$LINE),'POS')
    
  ef = Eff %>% dplyr::select(POS, grep("EFF",colnames(Eff)))%>% gather('Envi','eff',-POS)%>% 
    separate(Envi, c("A", "Cross","SITE"))%>%  mutate(eff=replace_na(eff, 0)) %>% filter(SITE %in% unique(Observed$SITE))
  
  SiteMean = read.csv(paste0('SiteMean_',t,'.csv'))
  Eff_Mean = ef %>% left_join(SiteMean) %>% arrange(SITE, POS)%>% mutate(Cross=str_replace(Cross,'add2','2'))

  ###rebuild the Genstat model
  Pred = NULL
  for (s in unique(Eff_Mean$SITE)){
    tmp1 = Eff_Mean %>% filter(SITE==s)
    for (l in 1:nrow(ValidationLine)){
      tmp2 = geno_QTLs[,c(1,(l+1),ncol(geno_QTLs))]
      print(c(t, s,l))
      colnames(tmp2)[2] ='LINE'
      tmp = tmp1 %>% left_join(tmp2) %>% mutate_at('LINE', funs(as.character)) %>%
        mutate(Code = case_when(
          stringr::str_detect(Cross,'add') & str_detect(LINE,'ac|ad' ) ~ -1,
          stringr::str_detect(Cross,'add') & str_detect(LINE,'bc|bd' ) ~ 1,
          stringr::str_detect(Cross,'2') & str_detect(LINE,'ac|bc' ) ~ -1,
          stringr::str_detect(Cross,'2') & str_detect(LINE,'ad|bd' ) ~ 1,
          stringr::str_detect(Cross,'dom') & str_detect(LINE,'ac|bd' ) ~ -1,
          stringr::str_detect(Cross,'dom') & str_detect(LINE,'bc|ad' ) ~ 1
        )) %>% mutate(EffQTL = eff*Code)
      PredValue = constants[i] + mean(tmp$mu) + sum(tmp$EffQTL,na.rm = T)
      PredValue2 = cbind(data.frame(Trait=t), data.frame(SITE=s),data.frame(LINE=ValidationLine$LINE[l]),data.frame(Pred=PredValue))
      Pred = rbind(Pred,PredValue2)
    }
  } 
  Obs_Pre = Observed %>% dplyr::select(SITE,LINE,t)  %>% left_join(Pred) 
  colnames(Obs_Pre)[3] ='Obs'
  Obs_Pre_All = rbind(Obs_Pre_All, Obs_Pre)
}
write.csv(Obs_Pre_All,'Observed_Predicted.csv', row.names = F)

# model_lm = Obs_Pre_All %>% split(list(.$Trait, .$SITE)) %>% map(~lm(Pred~Obs,data=.))
# R2_lm = model_lm %>% map(summary) %>% map_dbl('r.squared')
# eq_lm = paste0(expression(R^2),'==', round(R2_lm,2))
# eq_lm = data.frame(eq_lm,TRAITS=names(R2_lm)) 
# eq_lm = eq_lm %>% separate(TRAITS,c('Trait','SITE'))
# eq_lm$eq_lm = as.character(eq_lm$eq_lm)

Stats1 = Obs_Pre_All %>% group_by(Trait, SITE) %>% 
  summarize(RMSE= rmse(Pred, Obs), Pbias =pbias(Pred, Obs),
            r=cor(Pred, Obs, use='complete.obs'))

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


Obs_Pre_All$SITE = factor(Obs_Pre_All$SITE,levels=c('KBSM','CLMB','PKLE'))

  tmp = Obs_Pre_All %>% filter(Trait=='PL')
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

  
  tmp = Obs_Pre_All %>% filter(Trait=='PBN')
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
  
  tmp = Obs_Pre_All %>% filter(Trait=='SBN')
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
ggsave(file='Model_Evaluation_Figure5.tiff',F5, height =5, width = 8)

