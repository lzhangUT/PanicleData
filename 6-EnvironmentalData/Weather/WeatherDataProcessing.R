rm(list=ls())
library(tidyverse)
library(readxl)

setwd('c://Users/Li Zhang/Desktop/PanicleData/6-EnvironmentalData/Weather/')


##function to calculate daylength
day.length.hrs <- function(day, latitude) {
  daylength.coeff <- asin(0.39795*cos(0.2163108 + 2*atan(0.9671396*tan(0.00860*(day - 186)))))
  daylength.hrs <- 24 - (24/pi)*acos((sin(0.8333*pi/180) + sin(latitude*pi/180)*sin(daylength.coeff))/(cos(latitude*pi/180)*cos(daylength.coeff)))
  return(daylength.hrs)
}

latfile = read_excel('NSF_4WCR_Site Locations.xlsx')
latfile = latfile[-8,]

wthfiles = paste0(latfile$SiteCode,'.csv')

wthfile = lapply(wthfiles, function(x){
  print(x)
  tmp = read_csv(x)
  s = strsplit(x,'\\.')[[1]][1]
  lat = latfile %>% filter(SiteCode==s)
  tmp1 = tmp %>% rename(SRAD=`SRAD_MJ/m2`) %>% 
    mutate_at(vars('TMAX_C', 'TMIN_C','RAIN_mm','SRAD'), list(as.numeric))  %>% 
    mutate(DOY = 1:nrow(.), TMEAN_C = ((TMAX_C+TMIN_C)/2)) %>% mutate(DL = day.length.hrs(DOY, lat$Latitude)) %>% 
    dplyr::select (YEAR, MONTH, DAY, DOY, TMAX_C, TMIN_C, TMEAN_C, RAIN_mm, SRAD, DL)%>% mutate(SITE=s)
  return(tmp1)
 # write.csv(tmp1,paste0(s,'_withDL.csv'), row.names = F)
})

wth = do.call(rbind, wthfile)
wth_avg = wth %>% group_by(SITE) %>% summarise(TMEAN= mean(TMEAN_C,na.rm = T), RAIN=sum(RAIN_mm), 
                                               DL= mean(DL, na.rm = T), SRAD=mean(SRAD, na.rm = T) )
write.csv(wth_avg, 'TmeanRainDLSRAD_avg.csv',row.names = F)

sites = c("KING", "PKLE", "TMPL", "OVTN", "STIL", "CLMB", "MNHT", "LINC", "KBSM", "BRKG")
wth$SITE = factor(wth$SITE, levels = sites)
p1= ggplot(wth) + geom_boxplot(aes(SITE, TMEAN_C))+ylab(expression(paste('Mean Temperature (',~degree,'C)',sep='')))+xlab("")+
  theme_bw()+ theme(text = element_text(face = "bold", size = 13),axis.title=element_text(size=13,face="bold"))+
  theme(axis.text = element_text(face = "bold", size = 13))

rain = wth %>% group_by(SITE) %>% summarize(RAIN=sum(RAIN_mm))
p2 =ggplot(rain) + geom_col(aes(SITE, RAIN), alpha=0.7)+ylab('Total Rainfall (mm)')+xlab("")+
  theme_bw()+ theme(text = element_text(face = "bold", size = 13),axis.title=element_text(size=13,face="bold"))+
  theme(axis.text = element_text(face = "bold", size = 13))

pdf('Weather.pdf')
print(gridExtra::grid.arrange(p1,p2, nrow=2))

dev.off()
