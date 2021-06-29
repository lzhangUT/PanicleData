
rm(list=ls())
library(tidyverse)
library(readxl)
library(ggmap)
library(mapdata)
library(rgeos)
library(sp)
library(leaflet)
library(maptools)
library(mapview)
setwd('c://Users/Li Zhang/Desktop/PanicleData/0-SiteMap/')

unzip('SwitchgrassOccurance.zip', overwrite = T)
Occurances = read.csv('0018466-190621201848488.csv', sep = '\t')

tmp = Occurances %>% filter(countryCode=='US') %>% dplyr::select(decimalLatitude, decimalLongitude) %>% filter(complete.cases(.)) %>%
  filter(decimalLatitude<50 & decimalLongitude< (-60))

coordinates(tmp) = ~decimalLongitude + decimalLatitude

buf1 = gBuffer(tmp, width=.5, byid = T)
buf2 = gUnaryUnion(buf1)

mapStates = map("state", fill = TRUE, plot = FALSE)

sites = read_excel('NSF_4WCR_Site Locations, Plant and Harv Dates.xlsx')
sites = sites[-8,]
sites$SiteCode = factor(sites$SiteCode, levels=rev(c("KING", "PKLE", "TMPL", "OVTN", "STIL", "CLMB", "MNHT", "LINC", "KBSM", "BRKG")))


m = leaflet(data=mapStates) %>% addTiles()%>% setView(-90, 38,4.2) %>% addPolygons(fillColor = 'grey', color='black', weight = 0.7,opacity = 0.4)%>%
  addPolygons(data=buf2, color='green', weight = 1)%>% 
  addMarkers(data = sites,~Longitude, ~Latitude, label=~SiteCode,
                   labelOptions = labelOptions(noHide = T, textOnly = T,
                                               style = list('color'='blue','font-size'='18px')))

mapshot(m,file = "SiteMap2.pdf",remove_controls = c("zoomControl","homeButton", "layersControl"))



