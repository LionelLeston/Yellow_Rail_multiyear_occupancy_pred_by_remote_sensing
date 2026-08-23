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
#125 sites with at least 3 years of data from 2016 - 2021

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
#116 sites with 3+ years of data
#overlapping Hedley layer

#Value:   2  3  4  5  6  7  8 
#nsite:   1  5  3 18 29 30 30

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
                                     HabSuit4ClassNumeric))

#get static raw RS variables at locations in siteCovFinal

#sentinel1
#DPOL = cross polarization ration
sentinel1_DPOL_2016<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/Sentinel1_Indices_YERA_2016DPOL_summer.csv", header=TRUE)
sentinel1_DPOL_2016Filt<-sentinel1_DPOL_2016%>%
  mutate(year=2016,
         DPOL_150summer=mean_mean_150)%>%
  select(locationC,year,DPOL_150summer)%>%
  filter(locationC %in% levels(as.factor(siteCovFinal$locationC)))%>%
  arrange(locationC)#71

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
#These are year-specific values. Get an averaged value across years for each site
staticDPOLsummer<-sentinel1_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,DPOL_150summer)%>%
  group_by(locationC)%>%
  summarise(sDPOLsummer=mean(DPOL_150summer))

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
#These are year-specific values. Get an averaged value across years for each site
staticDeltaVH<-sentinel1_deltaVH%>%
  arrange(locationC,year)%>%
  select(locationC,year,deltaVH_150)%>%
  group_by(locationC)%>%
  summarise(sDeltaVH=mean(deltaVH_150))



#sentinel2
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

# sentinel2_2020<-read.csv("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/Sentinel2_Indices_YERASites_2020.csv", header=TRUE)
# sentinel2_2020Filt<-sentinel2_2020%>%
#   filter(month==6)%>%
#   select(locationC,year,ARI_mean_150,BSI_mean_150,   
#          EVI_mean_150,LAI_mean_150,MSAVI2_mean_150,
#          NBR_mean_150,NDMI_mean_150,NDRE_mean_150,NDVI_mean_150,  
#          NDWI_mean_150,REIP_mean_150,SAVI_mean_150,  
#          SI_mean_150)%>%
#   filter(locationC %in% levels(as.factor(m4$locationC)))

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
#These are year-specific values. Get an averaged value across years for each site
staticARI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,ARI_mean_150)%>%
  group_by(locationC)%>%
  summarise(sARI=mean(ARI_mean_150))

yearlyREIP<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,REIP_mean_150)%>%
  pivot_wider(names_from=year,values_from=REIP_mean_150)#7 years at each of 71 sites
#These are year-specific values. Get an averaged value across years for each site
staticREIP<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,REIP_mean_150)%>%
  group_by(locationC)%>%
  summarise(sREIP=mean(REIP_mean_150))

yearlyNDMI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDMI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDMI_mean_150)#7 years at each of 71 sites
#These are year-specific values. Get an averaged value across years for each site
staticNDMI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDMI_mean_150)%>%
  group_by(locationC)%>%
  summarise(sNDMI=mean(NDMI_mean_150))

yearlyNDVI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDVI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDVI_mean_150)#7 years at each of 71 sites
#These are year-specific values. Get an averaged value across years for each site
staticNDVI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDVI_mean_150)%>%
  group_by(locationC)%>%
  summarise(sNDVI=mean(NDVI_mean_150))

yearlyNDWI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDWI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDWI_mean_150)#7 years at each of 71 sites
#These are year-specific values. Get an averaged value across years for each site
staticNDWI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDWI_mean_150)%>%
  group_by(locationC)%>%
  summarise(sNDWI=mean(NDWI_mean_150))

yearlyBSI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,BSI_mean_150)%>%
  pivot_wider(names_from=year,values_from=BSI_mean_150)#7 years at each of 71 sites
#These are year-specific values. Get an averaged value across years for each site
staticBSI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,BSI_mean_150)%>%
  group_by(locationC)%>%
  summarise(sBSI=mean(BSI_mean_150))

yearlyNBR<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NBR_mean_150)%>%
  pivot_wider(names_from=year,values_from=NBR_mean_150)#7 years at each of 71 sites
#These are year-specific values. Get an averaged value across years for each site
staticNBR<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NBR_mean_150)%>%
  group_by(locationC)%>%
  summarise(sNBR=mean(NBR_mean_150))

siteCovFinal$sDPOLsummer<-staticDPOLsummer[,2]
siteCovFinal$sDeltaVH<-staticDeltaVH[,2]
siteCovFinal$sARI<-staticARI[,2]
siteCovFinal$sBSI<-staticBSI[,2]
siteCovFinal$sNBR<-staticNBR[,2]
siteCovFinal$sREIP<-staticREIP[,2]
siteCovFinal$sNDMI<-staticNDMI[,2]
siteCovFinal$sNDVI<-staticNDVI[,2]
siteCovFinal$sNDWI<-staticNDWI[,2]

sDPOLsummer<-staticDPOLsummer[,2]
sDeltaVH<-staticDeltaVH[,2]
sARI<-staticARI[,2]
sBSI<-staticBSI[,2]
sNBR<-staticNBR[,2]
sREIP<-staticREIP[,2]
sNDMI<-staticNDMI[,2]
sNDVI<-staticNDVI[,2]
sNDWI<-staticNDWI[,2]

umf<-unmarkedMultFrame(y=obsFiltFinal,
                       siteCovs = siteCovFinal,
                       obsCovs=list(time=matrix(4,115,25)),#5 years*5 rounds
                       yearlySiteCovs=list(Y_habSuit=yearlyHabSuit4Class,
                                           sDPOLsummer=as.matrix(cbind(sDPOLsummer,sDPOLsummer,sDPOLsummer,
                                                                       sDPOLsummer,sDPOLsummer)),
                                           sDeltaVH=as.matrix(cbind(sDeltaVH,sDeltaVH,sDeltaVH,
                                                                    sDeltaVH,sDeltaVH)),
                                           sNDMI=as.matrix(cbind(sNDMI,sNDMI,sNDMI,
                                                                 sNDMI,sNDMI)),
                                           sNDVI=as.matrix(cbind(sNDVI,sNDVI,sNDVI,
                                                                 sNDVI,sNDVI)),
                                           sNDWI=as.matrix(cbind(sNDWI,sNDWI,sNDWI,
                                                                 sNDWI,sNDWI)),
                                           sARI=as.matrix(cbind(sARI,sARI,sARI,
                                                                sARI,sARI)),
                                           sBSI=as.matrix(cbind(sBSI,sBSI,sBSI,
                                                                sBSI,sBSI)),
                                           sNBR=as.matrix(cbind(sNBR,sNBR,sNBR,
                                                                sNBR,sNBR)),
                                           sREIP=as.matrix(cbind(sREIP,sREIP,sREIP,
                                                                 sREIP,sREIP))),
                       numPrimary=5) 

#Scale predictors and create nonlinear terms
umf@siteCovs$sNDMI.s<-scale(umf@siteCovs$sNDMI, center = TRUE, scale = TRUE)
umf@siteCovs$sNDMI.s2<-umf@siteCovs$sNDMI.s^2
umf@siteCovs$sNDMI.sqrt<-sqrt(umf@siteCovs$sNDMI)

umf@siteCovs$sNDVI.s<-scale(umf@siteCovs$sNDVI, center = TRUE, scale = TRUE)
umf@siteCovs$sNDVI.s2<-umf@siteCovs$sNDVI.s^2
umf@siteCovs$sNDVI.sqrt<-sqrt(umf@siteCovs$sNDVI)

umf@siteCovs$sNDWI.s<-scale(umf@siteCovs$sNDWI, center = TRUE, scale = TRUE)
umf@siteCovs$sNDWI.s2<-umf@siteCovs$sNDWI.s^2
umf@siteCovs$sNDWI.sqrt<-sqrt(umf@siteCovs$sNDWI)

umf@siteCovs$sARI.s<-scale(umf@siteCovs$sARI, center = TRUE, scale = TRUE)
umf@siteCovs$sARI.s2<-umf@siteCovs$sARI.s^2
umf@siteCovs$sARI.sqrt<-sqrt(umf@siteCovs$sARI)

umf@siteCovs$sREIP.s<-scale(umf@siteCovs$sREIP, center = TRUE, scale = TRUE)
umf@siteCovs$sREIP.s2<-umf@siteCovs$sREIP.s^2
umf@siteCovs$sREIP.sqrt<-sqrt(umf@siteCovs$sREIP)

umf@siteCovs$sNBR.s<-scale(umf@siteCovs$sNBR, center = TRUE, scale = TRUE)
umf@siteCovs$sNBR.s2<-umf@siteCovs$sNBR.s^2
umf@siteCovs$sNBR.sqrt<-sqrt(umf@siteCovs$sNBR)

umf@siteCovs$sBSI.s<-scale(umf@siteCovs$sBSI, center = TRUE, scale = TRUE)
umf@siteCovs$sBSI.s2<-umf@siteCovs$sBSI.s^2
umf@siteCovs$sBSI.sqrt<-sqrt(umf@siteCovs$sBSI)

umf@siteCovs$sDPOLsummer.s<-scale(umf@siteCovs$sDPOLsummer, center = TRUE, scale = TRUE)
umf@siteCovs$sDPOLsummer.s2<-umf@siteCovs$sDPOLsummer.s^2
umf@siteCovs$sDPOLsummer.sqrt<-sqrt(umf@siteCovs$sDPOLsummer)

umf@siteCovs$sDeltaVH.s<-scale(umf@siteCovs$sDeltaVH, center = TRUE, scale = TRUE)
umf@siteCovs$sDeltaVH.s2<-umf@siteCovs$sDeltaVH.s^2
umf@siteCovs$sDeltaVH.sqrt<-sqrt(umf@siteCovs$sDeltaVH)


umf@yearlySiteCovs$sNDMI.s<-scale(umf@yearlySiteCovs$sNDMI, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sNDMI.s2<-umf@yearlySiteCovs$sNDMI.s^2
umf@yearlySiteCovs$sNDMI.sqrt<-sqrt(umf@yearlySiteCovs$sNDMI)

umf@yearlySiteCovs$sNDVI.s<-scale(umf@yearlySiteCovs$sNDVI, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sNDVI.s2<-umf@yearlySiteCovs$sNDVI.s^2
umf@yearlySiteCovs$sNDVI.sqrt<-sqrt(umf@yearlySiteCovs$sNDVI)

umf@yearlySiteCovs$sNDWI.s<-scale(umf@yearlySiteCovs$sNDWI, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sNDWI.s2<-umf@yearlySiteCovs$sNDWI.s^2
umf@yearlySiteCovs$sNDWI.sqrt<-sqrt(umf@yearlySiteCovs$sNDWI)

umf@yearlySiteCovs$sARI.s<-scale(umf@yearlySiteCovs$sARI, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sARI.s2<-umf@yearlySiteCovs$sARI.s^2
umf@yearlySiteCovs$sARI.sqrt<-sqrt(umf@yearlySiteCovs$sARI)

umf@yearlySiteCovs$sREIP.s<-scale(umf@yearlySiteCovs$sREIP, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sREIP.s2<-umf@yearlySiteCovs$sREIP.s^2
umf@yearlySiteCovs$sREIP.sqrt<-sqrt(umf@yearlySiteCovs$sREIP)

umf@yearlySiteCovs$sNBR.s<-scale(umf@yearlySiteCovs$sNBR, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sNBR.s2<-umf@yearlySiteCovs$sNBR.s^2
umf@yearlySiteCovs$sNBR.sqrt<-sqrt(umf@yearlySiteCovs$sNBR)

umf@yearlySiteCovs$sBSI.s<-scale(umf@yearlySiteCovs$sBSI, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sBSI.s2<-umf@yearlySiteCovs$sBSI.s^2
umf@yearlySiteCovs$sBSI.sqrt<-sqrt(umf@yearlySiteCovs$sBSI)

umf@yearlySiteCovs$sDPOLsummer.s<-scale(umf@yearlySiteCovs$sDPOLsummer, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$sDPOLsummer.s2<-umf@yearlySiteCovs$sDPOLsummer.s^2
umf@yearlySiteCovs$sDPOLsummer.sqrt<-sqrt(umf@yearlySiteCovs$sDPOLsummer)

umf@yearlySiteCovs$sDeltaVH.s<-scale(umf@yearlySiteCovs$sDeltaVH, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$sDeltaVH.s2<-umf@yearlySiteCovs$sDeltaVH.s^2
umf@yearlySiteCovs$sDeltaVH.sqrt<-sqrt(umf@yearlySiteCovs$sDeltaVH)


umf@siteCovs$HabSuit4ClassNumeric<-as.factor(umf@siteCovs$HabSuit4ClassNumeric)
umf@yearlySiteCovs$Y_habSuit<-as.factor(umf@yearlySiteCovs$Y_habSuit)

model1<-colext(psiformula= ~HabSuit4ClassNumeric, 
               gammaformula =  ~ Y_habSuit, 
               epsilonformula = ~ Y_habSuit,
               pformula = ~ 1, umf)
summary(model1)#AIC: 1682.768 (Null model)

NDMIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sNDMI.s, 
                      gammaformula =  ~ Y_habSuit+sNDMI.s, 
                      epsilonformula = ~ Y_habSuit+sNDMI.s,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1673.527
NDMIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sNDMI.s+sNDMI.s2, 
                      gammaformula =  ~ Y_habSuit+sNDMI.s+sNDMI.s2, 
                      epsilonformula = ~ Y_habSuit+sNDMI.s+sNDMI.s2,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1665.456
NDMIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sNDMI.sqrt, 
                      gammaformula =  ~ Y_habSuit+sNDMI.sqrt, 
                      epsilonformula = ~ Y_habSuit+sNDMI.sqrt,
                      pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#NDMI improves model, potentially more than yearly raw NDMI

NDVIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sNDVI.s, 
                   gammaformula =  ~ Y_habSuit+sNDVI.s, 
                   epsilonformula = ~ Y_habSuit+sNDVI.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1682.453
NDVIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sNDVI.s+sNDVI.s2, 
                   gammaformula =  ~ Y_habSuit+sNDVI.s+sNDVI.s2, 
                   epsilonformula = ~ Y_habSuit+sNDVI.s+sNDVI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1689.336
NDVIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sNDVI.sqrt, 
                   gammaformula =  ~ Y_habSuit+sNDVI.sqrt, 
                   epsilonformula = ~ Y_habSuit+sNDVI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#NDVI doesn't improve model

NDWIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sNDWI.s, 
                   gammaformula =  ~ Y_habSuit+sNDWI.s, 
                   epsilonformula = ~ Y_habSuit+sNDWI.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1679.271
NDWIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sNDWI.s+sNDWI.s2, 
                   gammaformula =  ~ Y_habSuit+sNDWI.s+sNDWI.s2, 
                   epsilonformula = ~ Y_habSuit+sNDWI.s+sNDWI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1684.913
NDWIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sNDWI.sqrt, 
                   gammaformula =  ~ Y_habSuit+sNDWI.sqrt, 
                   epsilonformula = ~ Y_habSuit+sNDWI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#NDWI improves model, more than yearly NDWI

ARIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sARI.s, 
                   gammaformula =  ~ Y_habSuit+sARI.s, 
                   epsilonformula = ~ Y_habSuit+sARI.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1682.018
ARIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sARI.s+sARI.s2, 
                   gammaformula =  ~ Y_habSuit+sARI.s+sARI.s2, 
                   epsilonformula = ~ Y_habSuit+sARI.s+sARI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1674.134
ARIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sARI.sqrt, 
                   gammaformula =  ~ Y_habSuit+sARI.sqrt, 
                   epsilonformula = ~ Y_habSuit+sARI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#ARI improves model, better than yearly raw

REIPModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sREIP.s, 
                  gammaformula =  ~ Y_habSuit+sREIP.s, 
                  epsilonformula = ~ Y_habSuit+sREIP.s,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1688.188
REIPModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sREIP.s+sREIP.s2, 
                  gammaformula =  ~ Y_habSuit+sREIP.s+sREIP.s2, 
                  epsilonformula = ~ Y_habSuit+sREIP.s+sREIP.s2,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1691.766
REIPModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sREIP.sqrt, 
                  gammaformula =  ~ Y_habSuit+sREIP.sqrt, 
                  epsilonformula = ~ Y_habSuit+sREIP.sqrt,
                  pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#REIP doesn't improve model

BSIModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sBSI.s, 
                   gammaformula =  ~ Y_habSuit+sBSI.s, 
                   epsilonformula = ~ Y_habSuit+sBSI.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1676.912
BSIModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sBSI.s+sBSI.s2, 
                   gammaformula =  ~ Y_habSuit+sBSI.s+sBSI.s2, 
                   epsilonformula = ~ Y_habSuit+sBSI.s+sBSI.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1679.043
BSIModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sBSI.sqrt, 
                   gammaformula =  ~ Y_habSuit+sBSI.sqrt, 
                   epsilonformula = ~ Y_habSuit+sBSI.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#BSI improves model, not as well as yearly raw

NBRModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sNBR.s, 
                   gammaformula =  ~ Y_habSuit+sNBR.s, 
                   epsilonformula = ~ Y_habSuit+sNBR.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1676.624
NBRModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sNBR.s+sNBR.s2, 
                   gammaformula =  ~ Y_habSuit+sNBR.s+sNBR.s2, 
                   epsilonformula = ~ Y_habSuit+sNBR.s+sNBR.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1677.287
NBRModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sNBR.sqrt, 
                   gammaformula =  ~ Y_habSuit+sNBR.sqrt, 
                   epsilonformula = ~ Y_habSuit+sNBR.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#NBR improves model, not as well as yearly NBR

DPOLsummerModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sDPOLsummer.s, 
                   gammaformula =  ~ Y_habSuit+sDPOLsummer.s, 
                   epsilonformula = ~ Y_habSuit+sDPOLsummer.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1662.664
DPOLsummerModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sDPOLsummer.s+sDPOLsummer.s2, 
                   gammaformula =  ~ Y_habSuit+sDPOLsummer.s+sDPOLsummer.s2, 
                   epsilonformula = ~ Y_habSuit+sDPOLsummer.s+sDPOLsummer.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1666.131
DPOLsummerModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sDPOLsummer.sqrt, 
                   gammaformula =  ~ Y_habSuit+sDPOLsummer.sqrt, 
                   epsilonformula = ~ Y_habSuit+sDPOLsummer.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#DPOLsummer improves model, better than yearly raw DPOLsummer

DeltaVHModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sDeltaVH.s, 
                         gammaformula =  ~ Y_habSuit+sDeltaVH.s, 
                         epsilonformula = ~ Y_habSuit+sDeltaVH.s,
                         pformula = ~ 1, umf, na.action = na.fail)#AIC: 1677.264
DeltaVHModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sDeltaVH.s+sDeltaVH.s2, 
                         gammaformula =  ~ Y_habSuit+sDeltaVH.s+sDeltaVH.s2, 
                         epsilonformula = ~ Y_habSuit+sDeltaVH.s+sDeltaVH.s2,
                         pformula = ~ 1, umf, na.action = na.fail)#AIC: 1681.018
DeltaVHModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sDeltaVH.sqrt, 
                         gammaformula =  ~ Y_habSuit+sDeltaVH.sqrt, 
                         epsilonformula = ~ Y_habSuit+sDeltaVH.sqrt,
                         pformula = ~ 1, umf, na.action = na.fail)#Error (sqrt of neg values I think)
#DeltaVH improves model, not as well as yearly raw DeltaVH

cor(umf@siteCovs$sNDMI, umf@siteCovs$sNBR)
#0.92
cor(umf@siteCovs$sNDMI, umf@siteCovs$sBSI)
#-0.97
cor(umf@siteCovs$sNBR, umf@siteCovs$sBSI)
#-0.98



#parallel processing and using dredge didn't work, so dropping it

#RS model without Hedley layer
#RSonlyModel<-colext(psiformula= ~grass.s+bare.s+shrub.s, 
#                         gammaformula =  ~ 1, 
#                         epsilonformula = ~ Y_DW_shrub.s,
#                         pformula = ~ 1, umf)#AIC: 1681.127
#grass.s and shrub.s could be dropped from psi model

#RSonlyModel<-colext(psiformula= ~bare.sqrt+grass.sqrt+shrub.sqrt+wet565.sqrt+elev.sqrt, #drop bare, shrub
#                    gammaformula =  ~ Y_distUpland.sqrt+Y_mrvbf.sqrt+Y_swi.sqrt, #drop MRVBF, distUpland?
#                    epsilonformula = ~ Y_tpi.s+Y_elev.sqrt+Y_wet565.sqrt, #drop wet565?
#                    pformula = ~ 1, umf)#AIC: 1648.279  Large or missing SE values

#RSonlyModel<-colext(psiformula= ~sDeltaVH.s+sNDMI.s+sNDMI.s2+sNDWI.s+sBSI.s, #drop sNDMI (keep sNBR for now), sNDWI
#                    gammaformula =  ~ sDPOLsummer.s+sNDWI.s+sARI.s, #drop sNDWI, sARI?
#                    epsilonformula = ~ sDeltaVH.s+sDPOLsummer.s+sNDMI.s+sNDMI.s2, #drop sNDMI?
#                    pformula = ~ 1, umf)#AIC: 1647.372

RSonlyModel<-colext(psiformula= ~sDeltaVH.s+sNBR.s, #drop sNDMI, sBSI, sNDWI
                    gammaformula =  ~ sDPOLsummer.s, #drop sNDWI, sARI?
                    epsilonformula = ~ sDeltaVH.s+sDPOLsummer.s, #drop sNDMI?
                    pformula = ~ 1, umf)#1633.023

FullModel<-colext(psiformula= ~HabSuit4ClassNumeric+sDeltaVH.s+sNBR.s, 
                  gammaformula =  ~ Y_habSuit+sDPOLsummer.s, 
                  epsilonformula = ~ Y_habSuit+sDeltaVH.s+sDPOLsummer.s,
                  pformula = ~ 1, umf)#AIC: 1635.211 #Large or missing SE values
#Basically the RS variables by themselves do a good enough job without the Hedley layer

mb.gof.test(RSonlyModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#86.5644 105.8298  34.9330  31.2366  63.9188 
#Total chi-square = 322.4826 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#115  157  171  185  282 
#Estimate of c-hat = 1.87   

mb.gof.test(FullModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
#86.9039 105.9319  35.0104  30.4049  62.4272 
#Total chi-square = 320.6783 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#117  161  176  194  310 
#Estimate of c-hat = 1.79 

#Get estimates of uncertainty in yearly occu. probability via 
#nonparametric bootstrapping
#use bootstrap to get estimates of psi (occu.) from derived parameter estimates

bootModel100x<-nonparboot(RSonlyModel, B=100)
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
ggsave(filename="2_outputs/model outputs/unmarked/best static raw/occu.ByRegion2016_2021_HedleyXStaticRaw.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)
#ggsave(filename="2_outputs/1_unmarked_outputs/November/Hedley 2016 2021/occu.ByRegion2016_2021_HedleyXStaticRaw.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)

