### Formatting Genstat Results##
rm(list=ls())
library(tidyverse)
library(reshape2)
library(xlsx)
setwd("c://Users/Li Zhang/Desktop/PanicleData/3-GenstatOutput/")

###formatting EFF files
QTLEffectFiles=intersect(list.files(pattern='EFF_WithoutSTIL'),list.files(pattern='xlsx'))
qtlEff = lapply(QTLEffectFiles, function(x){
  tmp = read.xlsx(x, sheetIndex = 1)
  tmp2 = tmp %>% gather('ENVI','VALUE',-POS) %>% separate(ENVI, c('TYPE','CROSS','SITE')) %>% dplyr::rename(INDEX=POS)%>%
    spread(TYPE,VALUE) %>% mutate(TRAIT = strsplit(strsplit(x,"_",fixed=T)[[1]][3],".", fixed= T)[[1]][1] )
})
qtlEff = do.call(rbind, qtlEff)

###formating LOD files
genomap = read.delim('4way_male_reduced_marker.map',sep = '', header = F)
genomap = genomap %>% filter(!str_detect(V1,'group')) %>% dplyr::rename(MARKER=V1, POS = V2) %>% 
  mutate_at(vars(POS, MARKER),list(as.character)) %>% mutate_at(vars(POS),list(as.numeric))

qtlLod = lapply(intersect(list.files(pattern = "LOD_WithoutSTIL"), list.files(pattern='.xlsx')), function(x){
  tmp = read.xlsx(x, sheetIndex = 1)
  tmp2 = tmp %>% mutate(MARKER = genomap$MARKER, POS=genomap$POS, TRAIT = strsplit(strsplit(x,"_",fixed=T)[[1]][3],".", fixed= T)[[1]][1] ) %>%
    dplyr::rename(LOD= qtl_minlogp_) %>% dplyr::select(MARKER, POS, LOD, TRAIT) %>% rownames_to_column('INDEX') %>% mutate_at(vars(INDEX),list(as.numeric))
})
qtlLod = do.call(rbind, qtlLod)

QTLlist = qtlEff %>% left_join(qtlLod) %>% dplyr::select(INDEX,TRAIT, MARKER,POS, LOD) %>% unique() %>% 
  mutate_at(vars(POS, MARKER),list(as.character))  %>% rowwise() %>% mutate(CHR = substr(MARKER, 5,6))%>%
  mutate_at(vars(POS),list(as.numeric)) %>% mutate(QTL = paste0(CHR,'@', round(POS,2))) 

### add the flank marker for these QTLs
QTL_flank = NULL
for (phe in unique(QTLlist$TRAIT)){
  tmp1 = QTLlist %>% filter(TRAIT==phe) 
  tmp2 = qtlLod %>% filter(TRAIT==phe) 

  for (i in 1:nrow(tmp1)){
    tmp3 = tmp2 %>% filter(MARKER== as.character(tmp1[i,'MARKER']))
    threshold = tmp3$LOD - 1.5
    m = tmp1$INDEX[i] -1 ##one marker lower than the marker with the highest LOD
    n = tmp1$INDEX[i] +1 ##one marker higher than the marker with the highest LOD
    
    while(tmp2$LOD[m] > threshold) {m=m-1}
    tmp1$lowposition[i] = tmp2$POS[m]
    tmp1$down_marker[i] = tmp2$MARKER[m]

    while(tmp2$LOD[n] > threshold){n=n+1}
    tmp1$highposition[i] = tmp2$POS[n]
    tmp1$up_marker[i] = tmp2$MARKER[n]
    
    if (tmp1$lowposition[i]> tmp1$POS[i]){
      tmp1$lowposition[i] = tmp2$POS[m+1]
      tmp1$down_marker[i] = tmp2$MARKER[m+1]
    }

      if (tmp1$highposition[i]==0){
      tmp1$highposition[i] = tmp2$POS[n-1]
      tmp1$up_marker[i] = tmp2$MARKER[n-1]
      }
  }
  QTL_flank = rbind(QTL_flank,tmp1)
}

####check if there is g x e, using qtlEFF file
for (i in 1:nrow(QTLlist)){
  tmp = qtlEff %>% filter(INDEX==QTLlist$INDEX[i] & TRAIT==QTLlist$TRAIT[i] & CROSS=='add')
  if (tmp$EFF[1] != tmp$EFF[2]){QTL_flank$QxE[i] = 'Yes'}
  if (tmp$EFF[1] == tmp$EFF[2]){QTL_flank$QxE[i] = 'No'}
}

QTL_flank = QTL_flank %>% mutate(ChrNo = str_extract(CHR, '[1-9]'), ChrAnnot = str_extract(CHR, '[A-Z]')) %>%
  mutate_at(vars('ChrNo'), list(as.numeric)) %>%  mutate(Dist = case_when(
    str_detect(ChrAnnot,'K') ~ (ChrNo*2-1),
    str_detect(ChrAnnot,'N') ~ ChrNo*2
  ))

write.csv(QTL_flank,'QTLlistForAllTraits_WithoutSTIL_flankmarker.csv',row.names = F)
file.copy('QTLlistForAllTraits_WithoutSTIL_flankmarker.csv', 
          'c:/Users/Li Zhang/Desktop/PanicleData/5-QTLSearch/',overwrite = T)

QxE = qtlEff %>% filter(INDEX %in% QTL_flank$INDEX[str_detect(QTL_flank$QxE,'Yes')])
write.csv(QxE,'QTLEffectWithInteraction_WithoutSTIL.csv', row.names = F)
file.copy('QTLEffectWithInteraction.csv', 'c:/Users/Li Zhang/Desktop/PanicleData/6-EnvironmentalData/Weather/',overwrite = T)

###plot summary of QTLs
load('sexAveraged.cross.rda')
QTL_flank$TRAIT = factor(QTL_flank$TRAIT, levels = c('PL','PBN','SBN','FL50','TC', 'BIO'))
#cols  =  c("goldenrod1","blue",'purple','brown','green','deeppink2',)
cols=palette(rainbow(6))
library(qtlTools)
# png('QTL_summary.png', width = 900)
pdf('QTL_SUMMARY.pdf', width = 9)
with(QTL_flank, segmentsOnMap(cross, phe = TRAIT, chr = CHR, 
                                                         l = lowposition, h = highposition,
                                                         peaklod = LOD, peakcM = POS,  
                                                         showPeaks = TRUE,
                                                         chrBuffer = c(0.15,0.1), 
                                                         tick.width=0.05, lwd="byLod",
                                                         col = cols,
                                                         leg.inset=.57, legendCex=1, 
                                                         legendPosition="topleft"))
#mtext('Summary of QTLs',side = 3)
for (phe in unique(QTL_flank$TRAIT)){
    tmp = subset(QTL_flank, TRAIT==phe)
   for (i in 1:nrow(tmp)){
     if (tmp$QxE[i]=='Yes')
    text(x=(tmp$Dist[i]+0.4), y=(tmp$lowposition[i]-2),'*',col = 'red',cex = 1.5)
   }
}

dev.off()

### plotting QTL effects and SE
qtlEff2 = qtlEff %>% filter(CROSS != 'dom')%>% mutate(CROSS = str_replace(CROSS,'add2','C x D')) %>% 
  mutate(CROSS = str_replace(CROSS,'add','A x B'))  %>% left_join(QTLlist) 

sites = c("KING", "PKLE", "TMPL", "OVTN",  "CLMB", "MNHT", "LINC", "KBSM", "BRKG")
qtlEff2$SITE = factor(qtlEff2$SITE, levels = sites)

#######plotting the pleiotropic QTL between traits
pdf('OverlappingQTLsForDifferentTraits_Pleitropy_freeXaxis.pdf',onefile = T)
qtlEff3 = qtlEff2 %>% filter(QTL=='2K@77.89'|QTL=='2K@70.12'|QTL=='2K@74.02')

ggplot(qtlEff3, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + ggtitle('2K@70.12/2K@74.02/2K@77.89')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff4 = qtlEff2 %>% filter(QTL=='3N@62.06'|QTL=='3N@63.97')

ggplot(qtlEff4, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + ggtitle('3N@63.97/3N@62.06')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff5 = qtlEff2 %>% filter(QTL=='4K@26.26'|QTL=='4K@3.89')

ggplot(qtlEff5, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + ggtitle('4K@3.89/4K@26.26')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff6 = qtlEff2 %>% filter(QTL=='9N@38.02'|QTL=='9N@22.11'|QTL=='9N@42.2'|QTL=='9N@47.97')

ggplot(qtlEff6, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + 
  ggtitle('9N@42.2/9N@22.11/9N@38.02/9N@38.02/9N@47.97')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff7 = qtlEff2 %>% filter(QTL=='2N@66.12'|QTL=='2N@70.02')

ggplot(qtlEff7, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + 
  ggtitle('2N@70.02/2N@66.12')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff8 = qtlEff2 %>% filter(QTL=='3K@38'|QTL=='3K@41.85'|QTL=='3K@57.6')

ggplot(qtlEff8, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + 
  ggtitle('3K@41.85/3K@38/3K@57.6')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff9 = qtlEff2 %>% filter(QTL=='5N@84.04'|QTL=='5N@79.93')

ggplot(qtlEff9, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + 
  ggtitle('5N@84.04/5N@84.04/5N@79.93')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff10 = qtlEff2 %>% filter(QTL=='7N@54.06')

ggplot(qtlEff10, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + 
  ggtitle('7N@54.06')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

qtlEff11 = qtlEff2 %>% filter(QTL=='9N@26.03'|QTL=='9N@22.11')
ggplot(qtlEff11, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~TRAIT, scales = 'free_x') + 
  ggtitle('9N@22.11/9N@26.03')+
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 15))+
  theme(strip.text = element_text(face="bold", size=15))+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()

dev.off()

###plot it

for (phe in unique(qtlEff2$TRAIT)){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  pdf(paste0('QTL_Effect_Plot_WithoutSTIL_',phe,'.pdf'), width = 14)
  print(
  ggplot(tmp1, aes(SITE, EFF)) + geom_point(size=2) + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme(axis.text = element_text(face = "bold", size = 15))+
    theme(strip.text = element_text(face="bold", size=15))+
   geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.3) + ylab('QTL Effect') +xlab("") +coord_flip()
  )
  dev.off()
}





