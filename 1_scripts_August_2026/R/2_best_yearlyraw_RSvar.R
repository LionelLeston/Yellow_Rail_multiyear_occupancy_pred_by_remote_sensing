library(unmarked)#for running dynamic occu. models based on maximum likelihood
#run models in unmarked package in R. Pick a best model and rerun in JAGS
library(dplyr)#for manipulating data frames
library(tidyr)#for manipulating data frames
library(sf)#for manipulating shapefile spatial data
library(terra)#for manipulating raster spatial data
library(MuMIn)#for likelihood-based model selection
library(AICcmodavg)#for Mackenzie-Bailey goodness-of-fit
library(ggplot2)#for plotting
library(grid)#for multipanel plots
library(gridExtra)#for multipanel plots
library(rjags)
library(runjags)
library(coda)
library(ggridges)
library(viridisLite)
theme_set(theme_minimal())

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

#get observations
finalYERAsitesShapefile<-st_read("0_data/0_raw/point counts/allYERAsites_June2025.shp")
str(finalYERAsitesShapefile)
#need clusterIDs corrected so that a couple too-large
#clusters are split up among wetlands
finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  mutate(clusterIDfinal=clusterID)%>%
  mutate(clusterIDfinal=case_when(locationC %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) ~ "2_8b",
                                  locationC %in% c(64,65,58) ~ "7_8b",
                                  locationC %in% c(60,61,62,63) ~ "7_8c",
                                  locationC %in% c(53,54,55,56) ~ "7_8d",
                                  !locationC %in% c(1:16,53:56,58,60:65) ~ clusterID))
levels(as.factor(finalYERAsitesShapefile$clusterIDfinal))


obs<-read.csv("0_data/0_raw/point counts/YERA_2025_COLEXT_DATA_DETECTIONHISTORYcombined.csv", header=TRUE)
#manually went through to combine detection histories stored on multiple rows
#for same sites, condensed detection histories to single row per site.

#This analysis will use 2016 - 2021 and sites
#from the whole YERA study area. 2016 is start year
#because that is earliest year we have Sentinel data for.
#2021 is the end year because southern sites were
#not systematically surveyed beyond 2021
obsFilt<-obs%>%
  select("yera2016_1","yera2016_2","yera2016_3","yera2016_4", 
         "yera2016_5","yera2017_1","yera2017_2","yera2017_3","yera2017_4", 
         "yera2017_5","yera2018_1","yera2018_2","yera2018_3","yera2018_4", 
         "yera2018_5","yera2019_1","yera2019_2","yera2019_3","yera2019_4", 
         "yera2019_5","yera2021_1","yera2021_2","yera2021_3","yera2021_4", 
         "yera2021_5")%>%
  mutate(numVisits=rowSums(!is.na(.)))%>%
  mutate(locationC=obs$locationC)

#,"yera2022_1","yera2022_2","yera2022_3","yera2022_4", 
#"yera2022_5","yera2023_1","yera2023_2","yera2023_3","yera2023_4", 
#"yera2023_5"

m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits
str(m2)
table(m2$clusterID)#125 sites, 8 clusters
#Note: not all sites might be covered by Hedley layer
#2_10  2_8  3_1  3_6  3_8  4_5  7_5  7_8 
#   5   15   12    7    7   36   29   14

#get covariates affecting initial occu.
initFilt<-m2%>%
  #select(locationC, clusterID)%>%
  arrange(locationC)
#123 sites with at least 3 years of data from 2016 - 2021

#Hedley et al. layer
Hedley2019<-rast("0_data/1_processed/6_HedleyEtAl2019/YERA_COLL_SUM_FINAL11.tif")
plot(Hedley2019$YERA_COLL_SUM_FINAL11)
raster_bbox<-st_bbox(Hedley2019)
bbox_polygon <- st_as_sfc(raster_bbox)
bbox_sf<-st_as_sf(bbox_polygon)

#Reproject points to Hedley layer's CRS
initFilt.rp<-initFilt%>%
  st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")%>%
  st_transform(crs=st_crs(Hedley2019))
#Extract values from Hedley layer to points
habSuitability<-terra::extract(Hedley2019, initFilt.rp)
hist(habSuitability$YERA_COLL_SUM_FINAL11)
table(habSuitability$YERA_COLL_SUM_FINAL11)
#115 sites with 3+ years of data
#overlapping Hedley layer

#Value:   2  3  4  5  6  7  8 
#nsite:   1  5  3 18 28 30 30

initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11
initFilt.rp$HabSuit4Class<-ifelse(initFilt.rp$HabSuitOrig<6, "5orless", initFilt.rp$HabSuitOrig)
#turns it to charactor
initFilt.rp$HabSuit4Class<-as.factor(initFilt.rp$HabSuit4Class)

obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(yera2016_1,yera2016_2,yera2016_3,yera2016_4,yera2016_5,
         yera2017_1,yera2017_2,yera2017_3,yera2017_4,yera2017_5,
         yera2018_1,yera2018_2,yera2018_3,yera2018_4,yera2018_5,
         yera2019_1,yera2019_2,yera2019_3,yera2019_4,yera2019_5,
         yera2021_1,yera2021_2,yera2021_3,yera2021_4,yera2021_5)
obsFiltFinal$geometry<-NULL
obsFiltFinal<-as.matrix(obsFiltFinal)

#Only use sites that have values from Hedley et al. layer
siteCovFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID, clusterIDfinal, HabSuit4Class)%>%
  mutate(HabSuit4ClassNumeric=as.numeric(HabSuit4Class))
siteCovFinal$geometry<-NULL

HabSuit4ClassNumeric<-siteCovFinal$HabSuit4ClassNumeric
yearlyHabSuit4Class<-as.matrix(cbind(HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric))

#get yearly raw RS variables at locations in siteCovFinal
#Sentinel-1
#DPOL = cross polarization ration
sentinel1_DPOL_2016<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2016DPOL_summer.csv", header=TRUE)
sentinel1_DPOL_2016Filt<-sentinel1_DPOL_2016%>%
  mutate(year=2016,
         DPOL_150summer=mean_mean_150)%>%
  select(locationC,year,DPOL_150summer)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)

sentinel1_DPOL_2017<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2017DPOL_summer.csv", header=TRUE)
sentinel1_DPOL_2017Filt<-sentinel1_DPOL_2017%>%
  mutate(year=2017,
         DPOL_150summer=mean_mean_150)%>%
  select(locationC,year,DPOL_150summer)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)

sentinel1_DPOL_2018<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2018DPOL_summer.csv", header=TRUE)
sentinel1_DPOL_2018Filt<-sentinel1_DPOL_2018%>%
  mutate(year=2018,
         DPOL_150summer=mean_mean_150)%>%
  select(locationC,year,DPOL_150summer)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)

sentinel1_DPOL_2019<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2019DPOL_summer.csv", header=TRUE)
sentinel1_DPOL_2019Filt<-sentinel1_DPOL_2019%>%
  mutate(year=2019,
         DPOL_150summer=mean_mean_150)%>%
  select(locationC,year,DPOL_150summer)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)

#N/A: data available but no YERA point counts in this year
#sentinel1_DPOL_2020<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2020DPOL_summer.csv", header=TRUE)
#sentinel1_DPOL_2020Filt<-sentinel1_DPOL_2020%>%
#  mutate(year=2020,
#         DPOL_150summer=mean_mean_150)%>%
#  select(locationC,year,DPOL_150summer)%>%
#  filter(locationC %in% levels(as.factor(m4$locationC)))%>%
#  arrange(locationC)

sentinel1_DPOL_2021<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2021DPOL_summer.csv", header=TRUE)
sentinel1_DPOL_2021Filt<-sentinel1_DPOL_2021%>%
  mutate(year=2021,
         DPOL_150summer=mean_mean_150)%>%
  select(locationC,year,DPOL_150summer)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)

sentinel1_all<-bind_rows(sentinel1_DPOL_2016Filt,
                         sentinel1_DPOL_2017Filt,
                         sentinel1_DPOL_2018Filt,
                         sentinel1_DPOL_2019Filt,
                         #sentinel1_DPOL_2020Filt,
                         sentinel1_DPOL_2021Filt)

#get variables in site major, year minor order
yearlyDPOLsummer<-sentinel1_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,DPOL_150summer)%>%
  pivot_wider(names_from=year,values_from=DPOL_150summer)#208 rows: 7 years at each of 208 sites

siteCovFinal$DPOLsummer<-yearlyDPOLsummer$`2016`

#VH = summer VH minus fall VH
sentinel1_deltaVH_2016<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2016VHdiff.csv", header=TRUE)
sentinel1_deltaVH_2016Filt<-sentinel1_deltaVH_2016%>%
  mutate(year=2016,
         deltaVH_150=mean_mean_150)%>%
  select(locationC,year,deltaVH_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)%>%
  mutate(deltaVH_150=ifelse(is.na(deltaVH_150), mean(deltaVH_150, na.rm=T), deltaVH_150))

sentinel1_deltaVH_2017<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2017VHdiff.csv", header=TRUE)
sentinel1_deltaVH_2017Filt<-sentinel1_deltaVH_2017%>%
  mutate(year=2017,
         deltaVH_150=mean_mean_150)%>%
  select(locationC,year,deltaVH_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)%>%
  mutate(deltaVH_150=ifelse(is.na(deltaVH_150), mean(deltaVH_150, na.rm=T), deltaVH_150))

sentinel1_deltaVH_2018<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2018VHdiff.csv", header=TRUE)
sentinel1_deltaVH_2018Filt<-sentinel1_deltaVH_2018%>%
  mutate(year=2018,
         deltaVH_150=mean_mean_150)%>%
  select(locationC,year,deltaVH_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)%>%
  mutate(deltaVH_150=ifelse(is.na(deltaVH_150), mean(deltaVH_150, na.rm=T), deltaVH_150))

sentinel1_deltaVH_2019<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2019VHdiff.csv", header=TRUE)
sentinel1_deltaVH_2019Filt<-sentinel1_deltaVH_2019%>%
  mutate(year=2019,
         deltaVH_150=mean_mean_150)%>%
  select(locationC,year,deltaVH_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)%>%
  mutate(deltaVH_150=ifelse(is.na(deltaVH_150), mean(deltaVH_150, na.rm=T), deltaVH_150))

sentinel1_deltaVH_2021<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2021VHdiff.csv", header=TRUE)
sentinel1_deltaVH_2021Filt<-sentinel1_deltaVH_2021%>%
  mutate(year=2021,
         deltaVH_150=mean_mean_150)%>%
  select(locationC,year,deltaVH_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)%>%
  mutate(deltaVH_150=ifelse(is.na(deltaVH_150), mean(deltaVH_150, na.rm=T), deltaVH_150))

sentinel1_deltaVH<-bind_rows(sentinel1_deltaVH_2016Filt,
                             sentinel1_deltaVH_2017Filt,
                             sentinel1_deltaVH_2018Filt,
                             sentinel1_deltaVH_2019Filt,
                             #sentinel1_deltaVH_2020Filt,
                             sentinel1_deltaVH_2021Filt)

#get variables in site major, year minor order
yearlyDeltaVH<-sentinel1_deltaVH%>%
  arrange(locationC,year)%>%
  select(locationC,year,deltaVH_150)%>%
  pivot_wider(names_from=year,values_from=deltaVH_150)#71 rows: 7 years at each of 65 sites

siteCovFinal$DeltaVH<-yearlyDeltaVH$`2016`

#Sentinel-2
sentinel2_2016<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2016.csv", header=TRUE)
sentinel2_2016Filt<-sentinel2_2016%>%
  filter(month==6)%>%
  select(locationC,year,ARI_mean_150,BSI_mean_150,   
         EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
         NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
         NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
         SI_mean_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))

sentinel2_2017<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2017.csv", header=TRUE)
sentinel2_2017Filt<-sentinel2_2017%>%
  filter(month==6)%>%
  select(locationC,year,ARI_mean_150,BSI_mean_150,   
         EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
         NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
         NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
         SI_mean_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))

sentinel2_2018<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2018.csv", header=TRUE)
sentinel2_2018Filt<-sentinel2_2018%>%
  filter(month==6)%>%
  select(locationC,year,ARI_mean_150,BSI_mean_150,   
         EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
         NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
         NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
         SI_mean_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))

sentinel2_2019<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2019.csv", header=TRUE)
sentinel2_2019Filt<-sentinel2_2019%>%
  filter(month==6)%>%
  select(locationC,year,ARI_mean_150,BSI_mean_150,   
         EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
         NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
         NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
         SI_mean_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))

#sentinel2_2020<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2020.csv", header=TRUE)
#sentinel2_2020Filt<-sentinel2_2020%>%
#  filter(month==6)%>%
#  select(locationC,year,ARI_mean_150,BSI_mean_150,   
#         EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
#         NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
#         NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
#         SI_mean_150)%>%
#  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))

sentinel2_2021<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2021.csv", header=TRUE)
sentinel2_2021Filt<-sentinel2_2021%>%
  filter(month==6)%>%
  select(locationC,year,ARI_mean_150,BSI_mean_150,   
         EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
         NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
         NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
         SI_mean_150)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))

sentinel2_all<-bind_rows(sentinel2_2016Filt,
                         sentinel2_2017Filt,
                         sentinel2_2018Filt,
                         sentinel2_2019Filt,
                         #sentinel2_2020Filt,
                         sentinel2_2021Filt)

#get variables in site major, year minor order
yearlyARI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,ARI_mean_150)%>%
  pivot_wider(names_from=year,values_from=ARI_mean_150)#7 years at each of 71 sites

siteCovFinal$ARI<-yearlyARI$`2016`

yearlyREIP<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,REIP_mean_150)%>%
  pivot_wider(names_from=year,values_from=REIP_mean_150)#7 years at each of 71 sites

siteCovFinal$REIP<-yearlyREIP$`2016`

yearlyNDMI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDMI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDMI_mean_150)#7 years at each of 71 sites

siteCovFinal$NDMI<-yearlyNDMI$`2016`

yearlyNDVI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDVI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDVI_mean_150)#7 years at each of 71 sites

siteCovFinal$NDVI<-yearlyNDVI$`2016`

yearlyNDWI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDWI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDWI_mean_150)#7 years at each of 71 sites

siteCovFinal$NDWI<-yearlyNDWI$`2016`

yearlyBSI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,BSI_mean_150)%>%
  pivot_wider(names_from=year,values_from=BSI_mean_150)#7 years at each of 71 sites

siteCovFinal$BSI<-yearlyBSI$`2016`

yearlyNBR<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NBR_mean_150)%>%
  pivot_wider(names_from=year,values_from=NBR_mean_150)#7 years at each of 71 sites

siteCovFinal$NBR<-yearlyNBR$`2016`

umf<-unmarkedMultFrame(y=obsFiltFinal,
                       siteCovs = siteCovFinal,
                       obsCovs=list(time=matrix(4,115,25)),#8 years*5 rounds
                       yearlySiteCovs=list(Y_habSuit=yearlyHabSuit4Class[,2:6],
                                           Y_ARI=yearlyARI[,2:6],
                                           Y_BSI=yearlyBSI[,2:6],
                                           Y_DeltaVH=yearlyDeltaVH[,2:6],
                                           Y_DPOLsummer=yearlyDPOLsummer[,2:6],
                                           Y_NBR=yearlyNBR[,2:6],
                                           Y_NDMI=yearlyNDMI[,2:6],
                                           Y_NDVI=yearlyNDVI[,2:6],
                                           Y_NDWI=yearlyNDWI[,2:6],
                                           Y_REIP=yearlyREIP[,2:6]),
                       numPrimary=5) 

#Scale predictors and create nonlinear terms
umf@siteCovs$DPOLsummer.s<-scale(umf@siteCovs$DPOLsummer, center = TRUE, scale = TRUE)
umf@siteCovs$DPOLsummer.s2<-umf@siteCovs$DPOLsummer.s^2
umf@siteCovs$DPOLsummer.sqrt<-sqrt(umf@siteCovs$DPOLsummer)

umf@siteCovs$DeltaVH.s<-scale(umf@siteCovs$DeltaVH, center = TRUE, scale = TRUE)
umf@siteCovs$DeltaVH.s2<-umf@siteCovs$DeltaVH.s^2
umf@siteCovs$DeltaVH.sqrt<-sqrt(umf@siteCovs$DeltaVH)

umf@siteCovs$ARI.s<-scale(umf@siteCovs$ARI, center = TRUE, scale = TRUE)
umf@siteCovs$ARI.s2<-umf@siteCovs$ARI.s^2
umf@siteCovs$ARI.sqrt<-sqrt(umf@siteCovs$ARI)

umf@siteCovs$BSI.s<-scale(umf@siteCovs$BSI, center = TRUE, scale = TRUE)
umf@siteCovs$BSI.s2<-umf@siteCovs$BSI.s^2
umf@siteCovs$BSI.sqrt<-sqrt(umf@siteCovs$BSI)

umf@siteCovs$NBR.s<-scale(umf@siteCovs$NBR, center = TRUE, scale = TRUE)
umf@siteCovs$NBR.s2<-umf@siteCovs$NBR.s^2
umf@siteCovs$NBR.sqrt<-sqrt(umf@siteCovs$NBR)

umf@siteCovs$REIP.s<-scale(umf@siteCovs$REIP, center = TRUE, scale = TRUE)
umf@siteCovs$REIP.s2<-umf@siteCovs$REIP.s^2
umf@siteCovs$REIP.sqrt<-sqrt(umf@siteCovs$REIP)

umf@siteCovs$NDMI.s<-scale(umf@siteCovs$NDMI, center = TRUE, scale = TRUE)
umf@siteCovs$NDMI.s2<-umf@siteCovs$NDMI.s^2
umf@siteCovs$NDMI.sqrt<-sqrt(umf@siteCovs$NDMI)

umf@siteCovs$NDVI.s<-scale(umf@siteCovs$NDVI, center = TRUE, scale = TRUE)
umf@siteCovs$NDVI.s2<-umf@siteCovs$NDVI.s^2
umf@siteCovs$NDVI.sqrt<-sqrt(umf@siteCovs$NDVI)

umf@siteCovs$NDWI.s<-scale(umf@siteCovs$NDWI, center = TRUE, scale = TRUE)
umf@siteCovs$NDWI.s2<-umf@siteCovs$NDWI.s^2
umf@siteCovs$NDWI.sqrt<-sqrt(umf@siteCovs$NDWI)

umf@yearlySiteCovs$DPOLsummer.s<-scale(umf@yearlySiteCovs$Y_DPOLsummer, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$DPOLsummer.s2<-umf@yearlySiteCovs$DPOLsummer.s^2
umf@yearlySiteCovs$DPOLsummer.sqrt<-sqrt(umf@yearlySiteCovs$Y_DPOLsummer)

umf@yearlySiteCovs$DeltaVH.s<-scale(umf@yearlySiteCovs$Y_DeltaVH, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$DeltaVH.s2<-umf@yearlySiteCovs$DeltaVH.s^2
umf@yearlySiteCovs$DeltaVH.sqrt<-sqrt(umf@yearlySiteCovs$Y_DeltaVH)

umf@yearlySiteCovs$ARI.s<-scale(umf@yearlySiteCovs$Y_ARI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$ARI.s2<-umf@yearlySiteCovs$ARI.s^2
umf@yearlySiteCovs$ARI.sqrt<-sqrt(umf@yearlySiteCovs$Y_ARI)

umf@yearlySiteCovs$BSI.s<-scale(umf@yearlySiteCovs$Y_BSI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$BSI.s2<-umf@yearlySiteCovs$BSI.s^2
umf@yearlySiteCovs$BSI.sqrt<-sqrt(umf@yearlySiteCovs$Y_BSI)

umf@yearlySiteCovs$NBR.s<-scale(umf@yearlySiteCovs$Y_NBR, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$NBR.s2<-umf@yearlySiteCovs$NBR.s^2
umf@yearlySiteCovs$NBR.sqrt<-sqrt(umf@yearlySiteCovs$Y_NBR)

umf@yearlySiteCovs$REIP.s<-scale(umf@yearlySiteCovs$Y_REIP, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$REIP.s2<-umf@yearlySiteCovs$REIP.s^2
umf@yearlySiteCovs$REIP.sqrt<-sqrt(umf@yearlySiteCovs$Y_REIP)

umf@yearlySiteCovs$NDMI.s<-scale(umf@yearlySiteCovs$Y_NDMI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$NDMI.s2<-umf@yearlySiteCovs$NDMI.s^2
umf@yearlySiteCovs$NDMI.sqrt<-sqrt(umf@yearlySiteCovs$Y_NDMI)

umf@yearlySiteCovs$NDVI.s<-scale(umf@yearlySiteCovs$Y_NDVI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$NDVI.s2<-umf@yearlySiteCovs$NDVI.s^2
umf@yearlySiteCovs$NDVI.sqrt<-sqrt(umf@yearlySiteCovs$Y_NDVI)

umf@yearlySiteCovs$NDWI.s<-scale(umf@yearlySiteCovs$Y_NDWI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$NDWI.s2<-umf@yearlySiteCovs$NDWI.s^2
umf@yearlySiteCovs$NDWI.sqrt<-sqrt(umf@yearlySiteCovs$Y_NDWI)


umf@siteCovs$HabSuit4ClassNumeric<-as.factor(umf@siteCovs$HabSuit4ClassNumeric)
umf@yearlySiteCovs$Y_habSuit<-as.factor(umf@yearlySiteCovs$Y_habSuit)

model1<-colext(psiformula= ~HabSuit4ClassNumeric, 
               gammaformula =  ~ Y_habSuit, 
               epsilonformula = ~ Y_habSuit,
               pformula = ~ 1, umf)
summary(model1)#AIC: 1682.768 (Null model)

NDVIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+NDVI.s, 
                      gammaformula =  ~ Y_habSuit+Y_NDVI, 
                      epsilonformula = ~ Y_habSuit+Y_NDVI,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1684.926
NDVIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+NDVI.s+NDVI.s2, 
                      gammaformula =  ~ Y_habSuit+NDVI.s+NDVI.s2, 
                      epsilonformula = ~ Y_habSuit+NDVI.s+NDVI.s2,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1689.587
NDVIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+NDVI.sqrt, 
                      gammaformula =  ~ Y_habSuit+NDVI.sqrt, 
                      epsilonformula = ~ Y_habSuit+NDVI.sqrt,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1684.911
#NDVI doesn't improve model

NDMIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.s, 
                   gammaformula =  ~ Y_habSuit+Y_NDMI, 
                   epsilonformula = ~ Y_habSuit+Y_NDMI,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1671.625
NDMIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.s+NDMI.s2, 
                   gammaformula =  ~ Y_habSuit+NDMI.s+NDMI.s2, 
                   epsilonformula = ~ Y_habSuit+NDMI.s+NDMI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1670.168
NDMIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.sqrt, 
                   gammaformula =  ~ Y_habSuit+NDMI.sqrt, 
                   epsilonformula = ~ Y_habSuit+NDMI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)
#Error in optim(unlist(tmb_pars), fn = tmb_obj$fn, gr = tmb_obj$gr, method = method,  : 
#initial value in 'vmmin' is not finite

#NDMI potentially improves model

NDWIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+NDWI.s, 
                   gammaformula =  ~ Y_habSuit+Y_NDWI, 
                   epsilonformula = ~ Y_habSuit+Y_NDWI,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1684.279
NDWIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+NDWI.s+NDWI.s2, 
                   gammaformula =  ~ Y_habSuit+NDWI.s+NDWI.s2, 
                   epsilonformula = ~ Y_habSuit+NDWI.s+NDWI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1683.668
NDWIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+NDWI.sqrt, 
                   gammaformula =  ~ Y_habSuit+NDWI.sqrt, 
                   epsilonformula = ~ Y_habSuit+NDWI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1683.533
#NDWI doesn't appear to improve model

REIPModel1<-colext(psiformula= ~HabSuit4ClassNumeric+REIP.s, 
                   gammaformula =  ~ Y_habSuit+Y_REIP, 
                   epsilonformula = ~ Y_habSuit+Y_REIP,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1674.281: NaN for error
REIPModel2<-colext(psiformula= ~HabSuit4ClassNumeric+REIP.s+REIP.s2, 
                   gammaformula =  ~ Y_habSuit+REIP.s+REIP.s2, 
                   epsilonformula = ~ Y_habSuit+REIP.s+REIP.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1679.066
REIPModel3<-colext(psiformula= ~HabSuit4ClassNumeric+REIP.sqrt, 
                   gammaformula =  ~ Y_habSuit+REIP.sqrt, 
                   epsilonformula = ~ Y_habSuit+REIP.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1674.84 Large or missing SE values
#REIP might improve model predictions but with caveats

ARIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+ARI.s, 
                   gammaformula =  ~ Y_habSuit+Y_ARI, 
                   epsilonformula = ~ Y_habSuit+Y_ARI,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1687.527 Large or missing SE values
ARIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+ARI.s+ARI.s2, 
                   gammaformula =  ~ Y_habSuit+ARI.s+ARI.s2, 
                   epsilonformula = ~ Y_habSuit+ARI.s+ARI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1685.894
ARIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+ARI.sqrt, 
                   gammaformula =  ~ Y_habSuit+ARI.sqrt, 
                   epsilonformula = ~ Y_habSuit+ARI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)
#Error in optim(unlist(tmb_pars), fn = tmb_obj$fn, gr = tmb_obj$gr, method = method,  : 
#initial value in 'vmmin' is not finite

#ARI doesn't appear to improve model

NBRModel1<-colext(psiformula= ~HabSuit4ClassNumeric+NBR.s, 
                  gammaformula =  ~ Y_habSuit+Y_NBR, 
                  epsilonformula = ~ Y_habSuit+Y_NBR,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1676.073
NBRModel2<-colext(psiformula= ~HabSuit4ClassNumeric+NBR.s+NBR.s2, 
                  gammaformula =  ~ Y_habSuit+NBR.s+NBR.s2, 
                  epsilonformula = ~ Y_habSuit+NBR.s+NBR.s2,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1679.067
NBRModel3<-colext(psiformula= ~HabSuit4ClassNumeric+NBR.sqrt, 
                  gammaformula =  ~ Y_habSuit+NBR.sqrt, 
                  epsilonformula = ~ Y_habSuit+NBR.sqrt,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1676.181 Large or missing SE values
#NBR appears to improve model

BSIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+BSI.s, 
                  gammaformula =  ~ Y_habSuit+Y_BSI, 
                  epsilonformula = ~ Y_habSuit+Y_BSI,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1675.123
BSIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+BSI.s+BSI.s2, 
                  gammaformula =  ~ Y_habSuit+BSI.s+BSI.s2, 
                  epsilonformula = ~ Y_habSuit+BSI.s+BSI.s2,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1677.044
BSIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+BSI.sqrt, 
                  gammaformula =  ~ Y_habSuit+BSI.sqrt, 
                  epsilonformula = ~ Y_habSuit+BSI.sqrt,
                  pformula = ~ 1, umf, na.action = na.fail)
#Error in seq.default(T, M * T, by = T) : wrong sign in 'by' argument

#BSI might improve model. Note: BSI, NDMI, and NBR cannot be together in same model

DPOLsummerPsiModel1<-colext(psiformula= ~HabSuit4ClassNumeric+DPOLsummer.s, 
                            gammaformula =  ~ Y_habSuit+Y_DPOLsummer, 
                            epsilonformula = ~ Y_habSuit+Y_DPOLsummer,
                            pformula = ~ 1, umf)#AIC: 1678.499
DPOLsummerPsiModel2<-colext(psiformula= ~HabSuit4ClassNumeric+DPOLsummer.s+DPOLsummer.s2, 
                            gammaformula =  ~ Y_habSuit+DPOLsummer.s+DPOLsummer.s2, 
                            epsilonformula = ~ Y_habSuit+DPOLsummer.s+DPOLsummer.s2,
                            pformula = ~ 1, umf)#AIC: 1682.81
DPOLsummerPsiModel3<-colext(psiformula= ~HabSuit4ClassNumeric+DPOLsummer.sqrt, 
                            gammaformula =  ~ Y_habSuit+DPOLsummer.sqrt, 
                            epsilonformula = ~ Y_habSuit+DPOLsummer.sqrt,
                            pformula = ~ 1, umf)#AIC: 1678.612 Large or missing SE values
#DPOLsummer might improve model

DeltaVHPsiModel1<-colext(psiformula= ~HabSuit4ClassNumeric+DeltaVH.s, 
                         gammaformula =  ~ Y_habSuit+Y_DeltaVH, 
                         epsilonformula = ~ Y_habSuit+Y_DeltaVH,
                         pformula = ~ 1, umf)#AIC: 1660.899
DeltaVHPsiModel2<-colext(psiformula= ~HabSuit4ClassNumeric+DeltaVH.s+DeltaVH.s2, 
                         gammaformula =  ~ Y_habSuit+DeltaVH.s+DeltaVH.s2, 
                         epsilonformula = ~ Y_habSuit+DeltaVH.s+DeltaVH.s2,
                         pformula = ~ 1, umf)#AIC: 1649.548
DeltaVHPsiModel3<-colext(psiformula= ~HabSuit4ClassNumeric+DeltaVH.sqrt, 
                         gammaformula =  ~ Y_habSuit+DeltaVH.sqrt, 
                         epsilonformula = ~ Y_habSuit+DeltaVH.sqrt,
                         pformula = ~ 1, umf)#
#Error in optim(unlist(tmb_pars), fn = tmb_obj$fn, gr = tmb_obj$gr, method = method,  : 
#initial value in 'vmmin' is not finite 

cor(siteCovFinal$NDMI, siteCovFinal$NDWI)
#NDVI and NDWI not correlated (0.02)
cor(siteCovFinal$NDMI, siteCovFinal$NBR)
#NDMI and NBR strongly positively correlated (0.85)
cor(siteCovFinal$NDMI, siteCovFinal$BSI)
#NDMI and BSI strongly negatively correlated (-0.91)
cor(siteCovFinal$NBR, siteCovFinal$BSI)
#NBR and BSI strongly negatively correlated (-0.95)
#pick one (NBR, BSI, NDMI)

#parallel processing and using dredge didn't work, so dropping it

#RS model without Hedley layer
#RSonlyModel<-colext(psiformula= ~NDMI.s+NDMI.s2+DPOLsummer.s+DeltaVH.s+DeltaVH.s2, 
#                         gammaformula =  ~ 1, 
#                         epsilonformula = ~ NDMI.s+DeltaVH.s+DeltaVH.s2,
#                         pformula = ~ 1, umf)#AIC: 1637.077
#DPOLsummer.s could be dropped

RSonlyModel<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+DeltaVH.s2, 
                    gammaformula =  ~ 1, 
                    epsilonformula = ~ NDMI.s+DeltaVH.s+DeltaVH.s2,
                    pformula = ~ 1, umf)#AIC: 1635.325

#"Full Model"
#FullModel<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.s+NDMI.s2+DPOLsummer.s+DeltaVH.s+DeltaVH.s2, 
#                         gammaformula =  ~ Y_habSuit, #n.s.
#                         epsilonformula = ~ Y_habSuit+NDMI.s+DeltaVH.s+DeltaVH.s2,
#                         pformula = ~ 1, umf)#AIC: 1643.738

#FullModel<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.s+NDMI.s2+DPOLsummer.s+DeltaVH.s+DeltaVH.s2, 
#                  gammaformula =  ~ 1, 
#                  epsilonformula = ~ Y_habSuit+NDMI.s+DeltaVH.s+DeltaVH.s2,
#                  pformula = ~ 1, umf)#AIC: 1638.379
#DPOLsummer.s could be dropped

FullModel<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.s+NDMI.s2+DeltaVH.s+DeltaVH.s2, 
                  gammaformula =  ~ 1, 
                  epsilonformula = ~ Y_habSuit+NDMI.s+DeltaVH.s+DeltaVH.s2,
                  pformula = ~ 1, umf)#AIC: 1624.751
#Initial (logit-scale):
# Estimate    SE      z  P(>|z|)
# (Intercept)              2.962 1.048  2.825 0.004724
# HabSuit4ClassNumeric2   -2.270 1.104 -2.057 0.039657
# HabSuit4ClassNumeric3   -0.982 1.073 -0.915 0.360139
# HabSuit4ClassNumeric4   -0.344 1.081 -0.319 0.749996
# NDMI.s                  -1.092 0.427 -2.556 0.010601
# NDMI.s2                 -1.227 0.388 -3.165 0.001549
# DeltaVH.s               -2.021 0.553 -3.653 0.000259
# DeltaVH.s2              -0.840 0.357 -2.356 0.018496

# Colonization (logit-scale):
#   Estimate    SE     z  P(>|z|)
# -0.947 0.174 -5.44 5.43e-08

# Extinction (logit-scale):
#   Estimate    SE      z P(>|z|)
# (Intercept)   -1.027 0.437 -2.349 0.01883
# Y_habSuit2     1.066 0.653  1.633 0.10246
# Y_habSuit3     1.511 0.586  2.579 0.00991
# Y_habSuit4    -0.388 0.548 -0.708 0.47882
# NDMI.s         0.544 0.272  2.003 0.04516
# DeltaVH.s      0.749 0.257  2.907 0.00364
# DeltaVH.s2     0.647 0.213  3.040 0.00237

#Detection (logit-scale):
#  Estimate     SE    z P(>|z|)
#0.478 0.0722 6.62 3.6e-11

mb.gof.test(FullModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Number of seasons:  5 

#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
# 87.9864  83.5710  36.2483  35.4254  58.1469 
#Total chi-square = 301.378 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#110  153  168  186  296 
#Estimate of c-hat = 1.76  

FullModelWOHedley<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+DeltaVH.s2, 
                  gammaformula =  ~ 1, 
                  epsilonformula = ~ NDMI.s+DeltaVH.s+DeltaVH.s2,
                  pformula = ~ 1, umf)
#Initial (logit-scale):
#Estimate    SE     z  P(>|z|)
#(Intercept)    1.798 0.570  3.16 0.001602
#NDMI.s        -0.786 0.326 -2.41 0.015852
#NDMI.s2       -1.009 0.348 -2.90 0.003722
#DeltaVH.s     -1.938 0.512 -3.79 0.000152
#DeltaVH.s2    -0.701 0.324 -2.16 0.030615
#Colonization (logit-scale):
#  Estimate    SE     z  P(>|z|)
#-0.953 0.174 -5.46 4.67e-08
#Extinction (logit-scale):
#  Estimate    SE     z P(>|z|)
#(Intercept)   -0.513 0.246 -2.09 0.03697
#NDMI.s         0.487 0.243  2.00 0.04515
#DeltaVH.s      0.799 0.244  3.28 0.00105
#DeltaVH.s2     0.550 0.200  2.75 0.00604
#Detection (logit-scale):
#  Estimate     SE   z  P(>|z|)
#0.477 0.0723 6.6 4.06e-11
#AIC: 1635.325 
#Number of sites: 115

mb.gof.test(FullModelWOHedley, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#87.1814  79.7869  37.2350  40.6519  59.1408 
#Total chi-square = 303.996 
#Number of bootstrap samples = 1000
#P-value = 0.001
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#117  154  170  189  307 
#Estimate of c-hat = 1.75 

#Get estimates of uncertainty in yearly occu. probability via 
#nonparametric bootstrapping
#use bootstrap to get estimates of psi (occu.) from derived parameter estimates

bootModel100x<-nonparboot(FullModel, B=100)
str(bootModel100x@projected)
#num [1:2, 1:8, 1:146] 0.46 0.54 0.694 0.306 0.688 ...
#this is a 3d array whose number is equal to 146(#sites)*8(#years)*2(1st vector for prob. absence,
#2nd vector for prob. occu.)
#http://www.statmethods.net/advstats/bootstrapping.html for a description of the bootstrap object structure
#boot( ) calls the statistic function R times. Each time, it generates a set of 
#random indices, with replacement, from the integers 1:nrow(data). These indices 
#are used within the statistic function to select a sample. The statistics are 
#calculated on the sample and the results are accumulated in the bootobject. The 
#bootobject structure includes

#basic bootstrap object element description 
#t0           The observed values of k statistics applied to the orginal data.  
#t            An R x k matrix where each row is a bootstrap replicate of the k statistics. 

#You can access these as bootobject$t0 and bootobject$t.

#Once you generate the bootstrap samples, print(bootobject) and plot(bootobject) can be used to examine the 
#results. If the results look reasonable, you can use boot.ci( ) function to obtain confidence intervals for 
#the statistic(s). 

bootModel100x@projected[2,,1]
#this is how you call up predicted probability of occu. based on the bootstrap samples for site 1
#in the colext unmarked Frame in each of the 7 years
#Yellow Rail predicted occu. probability has declined at this site over time.
bootModel100x@projected[1,,1]
#this is how you call up predicted probability that site is unoccupied based on the original data for site 1
#in the colext unmarked Frame
#These results are equivalent to 1 minus bootModel100x@projected[2,,1]

bootModel100x@projected[2,,9]#this is how you call up predicted probability of occu. based on the bootstrap samples for site 9
#in the colext unmarked Frame in each of the 7 years
bootModel100x@projected[2,,]#this is how you call up predicted probability of occu. based on the bootstrap samples 
#in the colext unmarked Frame, for all sites in each of the 7 years. 
#a 208 x 7 matrix
smoothed(bootModel100x)[2,] 
#average occu. probability in each year across all 208 sites
bootModel100x@projected.mean.bsse#the standard errors for projected occu.. Note that SE is same for probability of
#presence and probability of absence
bootModel100x@projected.mean.bsse[2,]#Also note that there is 1 SE calculated for each year but not for each site
#unlike probability of occu.

#Get regional projected changes in occu.
umf@siteCovs%>%
  select(locationC,clusterID)
#row numbers indicate how sites in a specific
#region are numbered in bootModel100x@projected

#Cluster 1 (formerly 3_6)
cluster3_6List<-list()#Audet
for (j in c(110:115)){
  siteRownum<-j
  m.cluster3_6List<-bootModel100x@projected[2,,j]
  m.cluster3_6List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster3_6List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster3_6List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster3_6List,
                                  l95=m.cluster3_6List_l95,
                                  u95=m.cluster3_6List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster3_6<-do.call(rbind, cluster3_6List)

Occucluster3_6 <- ggplot(data=cluster3_6, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 1")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(-0.25,1.25)+xlim(2016,2021)

#Cluster 2
cluster4_5List<-list()#West Mine
for (j in c(45:59,61,75,78,80:81,86:89,93,97:99,101:105,107:109)){
  siteRownum<-j
  m.cluster4_5List<-bootModel100x@projected[2,,j]
  m.cluster4_5List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster4_5List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster4_5List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster4_5List,
                                  l95=m.cluster4_5List_l95,
                                  u95=m.cluster4_5List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster4_5<-do.call(rbind, cluster4_5List)

Occucluster4_5 <- ggplot(data=cluster4_5, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 2")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster4_5 (West Mine)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(0,1)+xlim(2016,2021)

umf@siteCovs%>%
  select(locationC,clusterID)
#row numbers indicate how sites in a specific
#region are numbered in bootModel100x@projected

#Cluster 3
cluster7_5List<-list()#East Mine
for (j in c(60,62:74,76:77,79,82:85,90:92,94:96,100,106)){
  siteRownum<-j
  m.cluster7_5List<-bootModel100x@projected[2,,j]
  m.cluster7_5List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster7_5List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster7_5List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster7_5List,
                                  l95=m.cluster7_5List_l95,
                                  u95=m.cluster7_5List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster7_5<-do.call(rbind, cluster7_5List)

Occucluster7_5 <- ggplot(data=cluster7_5, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 3")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster7_5 (East Mine)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(0,1)+xlim(2016,2021)

#Cluster 4 (No sites with 3+ years of visits between 2016 - 2021)
# cluster2_10List<-list()#Gipsy Lake
# for (j in c(47:49)){
#   siteRownum<-j
#   m.cluster2_10List<-bootModel100x@projected[2,,j]
#   m.cluster2_10List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
#   m.cluster2_10List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
#   cluster2_10List[[j]]<-data.frame(site=siteRownum,
#                                    occ=m.cluster2_10List,
#                                    l95=m.cluster2_10List_l95,
#                                    u95=m.cluster2_10List_u95,
#                                    year=c(2013,2014,2015,2016,2017,2018,2019,2021))
# }
# cluster2_10<-do.call(rbind, cluster2_10List)
# 
# Occucluster2_10 <- ggplot(data=cluster2_10, aes(x=year, y=occ, group=site)) +
#   geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
#   geom_line(aes(x=year, y=occ, group=site), colour="blue")+
#   my.theme+ggtitle("Cluster 4")+xlab("Year")+ylab("Predicted occu.")+
#   #my.theme+ggtitle("Cluster2_10 (Gipsy Lake)")+xlab("Year")+ylab("Predicted occu.")+
#   ylim(-0.25,1.25)+xlim(2013,2021)


#Cluster 5
cluster3_1List<-list()#881 between Quigley and Chard
for (j in c(33:44)){
  siteRownum<-j
  m.cluster3_1List<-bootModel100x@projected[2,,j]
  m.cluster3_1List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster3_1List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster3_1List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster3_1List,
                                  l95=m.cluster3_1List_l95,
                                  u95=m.cluster3_1List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster3_1<-do.call(rbind, cluster3_1List)

Occucluster3_1 <- ggplot(data=cluster3_1, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+ 
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 5")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster3_1 (Highway 881)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(-0.25,1.25)+xlim(2016,2021)


#Cluster 6

cluster7_8List<-list()#Lac La Biche
for (j in c(23:32)){
  siteRownum<-j
  m.cluster7_8List<-bootModel100x@projected[2,,j]
  m.cluster7_8List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster7_8List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster7_8List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster7_8List,
                                  l95=m.cluster7_8List_l95,
                                  u95=m.cluster7_8List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster7_8<-do.call(rbind, cluster7_8List)

Occucluster7_8 <- ggplot(data=cluster7_8, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 6")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster7_8 (Lac La Biche)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(-0.25,1.25)+xlim(2016,2021)


#Cluster 7
cluster3_8List<-list()#North of La Corey
for (j in c(16:22)){
  siteRownum<-j
  m.cluster3_8List<-bootModel100x@projected[2,,j]
  m.cluster3_8List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster3_8List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster3_8List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster3_8List,
                                  l95=m.cluster3_8List_l95,
                                  u95=m.cluster3_8List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster3_8<-do.call(rbind, cluster3_8List)

Occucluster3_8 <- ggplot(data=cluster3_8, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 7")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster3_8 (North of La Corey)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(-0.25,1.25)+xlim(2016,2021)

umf@siteCovs%>%
  select(locationC,clusterID)
#Cluster 8
cluster2_8List<-list()#Imperial Oil/Cold Lake
for (j in c(1:15)){
  siteRownum<-j
  m.cluster2_8List<-bootModel100x@projected[2,,j]
  m.cluster2_8List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster2_8List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster2_8List[[j]]<-data.frame(site=siteRownum,
                                  occ=m.cluster2_8List,
                                  l95=m.cluster2_8List_l95,
                                  u95=m.cluster2_8List_u95,
                                  year=c(2016,2017,2018,2019,2021))
}
cluster2_8<-do.call(rbind, cluster2_8List)

Occucluster2_8 <- ggplot(data=cluster2_8, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), alpha = 0.2, fill = "blue")+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 8")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster2_8 (Imperial Oil/Cold Lake)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(-0.25,1.25)+xlim(2016,2021)

#Need cluster 9? (4_8) still


multipanelPlot<-grid.arrange(Occucluster3_6,#cluster 1 (Audet Lake)
                             Occucluster4_5,#cluster 2 (West Mine)
                             Occucluster7_5,#cluster 3 (East Mine)
                             #Occucluster2_10,#cluster 4 (Gipsy Lake)
                             Occucluster3_1,#cluster 5 (Highway 881)
                             Occucluster7_8,#cluster 6 (Lac La Biche)
                             Occucluster3_8,#cluster 7 (N of La Corey)
                             Occucluster2_8,#cluster 8 (Imperial Oil/Cold Lake)
                             nrow=3,ncol=3)
ggsave(filename="2_outputs/model output/unmarked/best yearly raw/occu.ByRegion2016_2021_HedleyXYearlyRaw.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)
#ggsave(filename="2_outputs/1_unmarked_outputs/November/Hedley 2016 2021/occu.ByRegion2016_2021_HedleyXYearlyRaw.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)

