rm(list=ls())
setwd('c://Users/Li Zhang/Desktop/PanicleData/6-EnvironmentalData/Weather/')

library(tidyverse)
library(MASS)
library(leaps)
library(bestglm)

eff = read.csv('QTLEffectWithInteraction.csv')
eff = eff %>% filter(CROSS != 'dom')
env = read.csv('TmeanRainDLSRAD_avg.csv')

env =env %>%  mutate(TR1 =TMEAN*SRAD, TR2 = TMEAN*RAIN, TR3=SRAD*RAIN)

Res = NULL
for (t in unique(eff$TRAIT)){
  for (i in unique(eff$CROSS)){
    tmp = eff %>% filter(TRAIT==t & CROSS==i)
    for (j in unique(tmp$INDEX)){
      print(c(i,j,t))
      tmp2 = tmp %>% filter(INDEX ==j) %>% left_join(env) %>% dplyr::select(TMEAN, RAIN, SRAD,DL,EFF) 
       # mutate(TR1 =TMEAN*SRAD, TR2 = TMEAN*RAIN, TR3=SRAD*RAIN) %>% dplyr::select(TMEAN, TR1,RAIN,SRAD,  TR2,TR3,EFF)
      rg = bestglm(tmp2, IC='BIC')
     x = summary(rg$BestModel)
     y = data.frame(x$coefficients)
     y = y %>% rownames_to_column('terms')%>% mutate(TRAIT=t, CROSS=i,INDEX=j)
     Res = rbind(Res,y)
    }
  }
}

QTL = read.csv('QTLlistForAllTraits_flankmarker.csv')
Res = Res %>% left_join(QTL) 
write.csv(Res,'VariableSelection.csv',row.names = F)
