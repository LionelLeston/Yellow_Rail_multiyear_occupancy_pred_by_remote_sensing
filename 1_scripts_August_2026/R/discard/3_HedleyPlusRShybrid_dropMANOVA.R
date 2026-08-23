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


#Hedley plus RS hybrid models
#get observations
finalYERAsitesShapefile<-st_read("0_data/0_raw/point counts/allYERAsites_June2025.shp")
finalYERAsitesShapefile$lat<-st_coordinates(finalYERAsitesShapefile)[,2]
finalYERAsitesShapefile$long<-st_coordinates(finalYERAsitesShapefile)[,1]
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

siteCovFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, lat, long, clusterID, clusterIDfinal, HabSuit4Class)%>%
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

Latitude<-siteCovFinal$lat
yearlyLatitude<-as.matrix(cbind(Latitude,Latitude,Latitude,Latitude,Latitude))

Longitude<-siteCovFinal$long
yearlyLongitude<-as.matrix(cbind(Longitude,Longitude,Longitude,Longitude,Longitude))

#get top yearly raw and static RS variables at locations in siteCovFinal
#Sentinel-1
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
#since static DPOLsummer did slightly better than yearly DPOLsummer AIC-wise
staticDPOLsummer<-sentinel1_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,DPOL_150summer)%>%
  group_by(locationC)%>%
  summarise(sDPOLsummer=mean(DPOL_150summer))

siteCovFinal$sDPOLsummer<-staticDPOLsummer[,2]
#siteCovFinal$DPOLsummer<-yearlyDPOLsummer$`2016`
sDPOLsummer<-staticDPOLsummer[,2]

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
  pivot_wider(names_from=year,values_from=deltaVH_150)

siteCovFinal$DeltaVH<-yearlyDeltaVH$`2016`
#We won't bother with static DeltaVH since it didn't perform as well AIC-wise as yearly 

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
#static NDMI and NBR were top static raw predictors
#static NDMI may do better AIC-wise than yearly NDMI
#yearly NBR may do better AIC-wise than static NBR
yearlyARI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,ARI_mean_150)%>%
  pivot_wider(names_from=year,values_from=ARI_mean_150)

siteCovFinal$ARI<-yearlyARI$`2016`

yearlyREIP<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,REIP_mean_150)%>%
  pivot_wider(names_from=year,values_from=REIP_mean_150)

siteCovFinal$REIP<-yearlyREIP$`2016`

yearlyNDMI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDMI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDMI_mean_150)
#These are year-specific values. Get an averaged value across years for each site
staticNDMI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDMI_mean_150)%>%
  group_by(locationC)%>%
  summarise(sNDMI=mean(NDMI_mean_150))

siteCovFinal$NDMI<-yearlyNDMI$`2016`
siteCovFinal$sNDMI<-staticNDMI[,2]

yearlyNDVI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDVI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDVI_mean_150)#7 years at each of 71 sites

siteCovFinal$NDVI<-yearlyNDVI$`2016`

yearlyNDWI<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NDWI_mean_150)%>%
  pivot_wider(names_from=year,values_from=NDWI_mean_150)

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
#These are year-specific values. Get an averaged value across years for each site
staticNBR<-sentinel2_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,NBR_mean_150)%>%
  group_by(locationC)%>%
  summarise(sNBR=mean(NBR_mean_150))
siteCovFinal$sNBR<-staticNBR[,2]

siteCovFinal$NBR<-yearlyNBR$`2016`

#get yearly and static derived RS variables at locations in siteCovFinal
#Dynamic World Layers
dw2016<-read.csv("0_data/1_processed/5_DynamicWorld/ExtractedData2016.csv", header=TRUE)
str(dw2016)
dw2016Filt<-dw2016%>%
  select(locationC,bare,built,crops,flooded_vegetation,trees,water,grass,shrub_and_scrub,snow_and_ice)%>%
  mutate(year=2016,
         woody=trees+shrub_and_scrub)%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  arrange(locationC)

dw2017<-read.csv("0_data/1_processed/5_DynamicWorld/ExtractedData2017.csv", header=TRUE)
str(dw2017)
dw2017Filt<-dw2017%>%
  select(locationC,bare,built,crops,flooded_vegetation,trees,water,grass,shrub_and_scrub,snow_and_ice)%>%
  mutate(year=2017,
         woody=trees+shrub_and_scrub)%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  arrange(locationC)

dw2018<-read.csv("0_data/1_processed/5_DynamicWorld/ExtractedData2018.csv", header=TRUE)
str(dw2018)
dw2018Filt<-dw2018%>%
  select(locationC,bare,built,crops,flooded_vegetation,trees,water,grass,shrub_and_scrub,snow_and_ice)%>%
  mutate(year=2018,
         woody=trees+shrub_and_scrub)%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  arrange(locationC)

dw2019<-read.csv("0_data/1_processed/5_DynamicWorld/ExtractedData2019.csv", header=TRUE)
str(dw2019)
dw2019Filt<-dw2019%>%
  select(locationC,bare,built,crops,flooded_vegetation,trees,water,grass,shrub_and_scrub,snow_and_ice)%>%
  mutate(year=2019,
         woody=trees+shrub_and_scrub)%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  arrange(locationC)

dw2021<-read.csv("0_data/1_processed/5_DynamicWorld/ExtractedData2021.csv", header=TRUE)
str(dw2021)
dw2021Filt<-dw2021%>%
  select(locationC,bare,built,crops,flooded_vegetation,trees,water,grass,shrub_and_scrub,snow_and_ice)%>%
  mutate(year=2021,
         woody=trees+shrub_and_scrub)%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  arrange(locationC)


dw_all<-bind_rows(dw2016Filt,
                  dw2017Filt,
                  dw2018Filt,
                  dw2019Filt,
                  dw2021Filt)

#get variables in site major, year minor order
#bare ground (static) and wet565, elevation, tpi, swi made it into best static derived
#model and were more parsimonious without than with Hedley layer
#bare ground (dynamic), shrub (dynamic), and grass (dynamic) made it into best yearly
#derived model. Bare (static) performed better AIC-wise than bare (dynamic)
yearlyDW_flooded<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,flooded_vegetation)%>%
  pivot_wider(names_from=year,values_from=flooded_vegetation)

yearlyDW_water<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,water)%>%
  pivot_wider(names_from=year,values_from=water)

yearlyDW_woody<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,woody)%>%
  pivot_wider(names_from=year,values_from=woody)

yearlyDW_bare<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,bare)%>%
  pivot_wider(names_from=year,values_from=bare)
#These are year-specific values. Get an averaged value across years for each site
staticDW_bare<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,bare)%>%
  group_by(locationC)%>%
  summarise(sbare=mean(bare))

yearlyDW_grass<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,grass)%>%
  pivot_wider(names_from=year,values_from=grass)
#These are year-specific values. Get an averaged value across years for each site
staticDW_grass<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,grass)%>%
  group_by(locationC)%>%
  summarise(sgrass=mean(grass))

yearlyDW_shrub<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,shrub_and_scrub)%>%
  pivot_wider(names_from=year,values_from=shrub_and_scrub)


#get more siteCovs for initial occupancy models
siteCovFinal$grass<-yearlyDW_grass$`2016`
siteCovFinal$shrub<-yearlyDW_shrub$`2016`
siteCovFinal$bare<-yearlyDW_bare$`2016`

sbare<-staticDW_bare$sbare
sgrass<-staticDW_grass$sgrass
siteCovFinal$sbare<-staticDW_bare[,2]
siteCovFinal$sgrass<-staticDW_grass[,2]
sDPOLsummer<-staticDPOLsummer[,2]
sNDMI<-staticNDMI[,2]
sNBR<-staticNBR[,2]
sbareSY<-as.matrix(cbind(sbare,
                         sbare,
                         sbare,
                         sbare,
                         sbare))
sgrassSY<-as.matrix(cbind(sgrass,
                          sgrass,
                          sgrass,
                          sgrass,
                          sgrass))
#topography
topo<-read.csv("0_data/1_processed/2_topography/topography150m_June2025.csv", header=TRUE)
topoFilt<-topo%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  mutate(elev=ALOS_dem_raster..MEAN.,
         swi=Topographic.Wetness.Index..MEAN.,
         tpi=Topographic.Position.Index..MEAN.,
         mrvbf=MRVBF..MEAN.)%>%
  select(locationC, elev, swi, tpi, mrvbf)

meanElev<-mean(topoFilt$elev, na.rm=TRUE)
meanSWI<-mean(topoFilt$swi, na.rm=TRUE)
meanTPI<-mean(topoFilt$tpi, na.rm=TRUE)
meanMRVBF<-mean(topoFilt$mrvbf, na.rm=TRUE)

topoFilt<-topoFilt%>%
  mutate(elev=case_when(is.na(elev) ~ meanElev,
                        !is.na(elev) ~ elev))%>%
  mutate(swi=case_when(is.na(swi) ~ meanSWI,
                       !is.na(swi) ~ swi))%>%
  mutate(tpi=case_when(is.na(tpi) ~ meanTPI,
                       !is.na(tpi) ~ tpi))%>%
  mutate(mrvbf=case_when(is.na(mrvbf) ~ meanMRVBF,
                         !is.na(mrvbf) ~ mrvbf))

siteCovFinal$elev<-topoFilt$elev
siteCovFinal$swi<-topoFilt$swi
siteCovFinal$tpi<-topoFilt$tpi
siteCovFinal$mrvbf<-topoFilt$mrvbf

Elev<-siteCovFinal$elev
Y_elev<-as.matrix(cbind(Elev,Elev,Elev,Elev,Elev))
SWI<-siteCovFinal$swi
Y_swi<-as.matrix(cbind(SWI,SWI,SWI,SWI,SWI))
TPI<-siteCovFinal$tpi
Y_tpi<-as.matrix(cbind(TPI,TPI,TPI,TPI,TPI))
MRVBF<-siteCovFinal$mrvbf
Y_mrvbf<-as.matrix(cbind(MRVBF,MRVBF,MRVBF,MRVBF,MRVBF))


wet<-read.csv("0_data/1_processed/3_wetlands/wetlandProb_uplandDist.csv", header=TRUE)
wet<-wet%>%
  filter(locationC %in% siteCovFinal$locationC)%>%
  select(locationC, wetlandProb_150m, wetlandProb_565m, distUpland)%>%
  arrange(locationC)
cor(wet$wetlandProb_150m, wet$distUpland)#0.65

siteCovFinal$wetlandProb_150m<-wet$wetlandProb_150m
siteCovFinal$wetlandProb_565m<-wet$wetlandProb_565m
siteCovFinal$distUpland<-wet$distUpland

Wet565<-wet$wetlandProb_565m
Y_wet565<-as.matrix(cbind(Wet565,Wet565,Wet565,Wet565,Wet565))



umf<-unmarkedMultFrame(y=obsFiltFinal,
                       siteCovs = siteCovFinal,
                       obsCovs=list(time=matrix(4,115,25)),#8 years*5 rounds
                       yearlySiteCovs=list(Y_habSuit=yearlyHabSuit4Class[,2:6],
                                           Y_latitude=yearlyLatitude,
                                           Y_longitude=yearlyLongitude,
                                           Y_elev=Y_elev,
                                           Y_swi=Y_swi,
                                           Y_tpi=Y_tpi,
                                           Y_mrvbf=Y_mrvbf,
                                           #Y_wet150=Y_wet150,
                                           Y_wet565=Y_wet565,
                                           Y_DW_bare=yearlyDW_bare[,2:6],
                                           Y_DW_grass=yearlyDW_grass[,2:6],
                                           Y_DW_flooded=yearlyDW_flooded[,2:6],
                                           Y_DW_water=yearlyDW_water[,2:6],
                                           Y_DW_woody=yearlyDW_woody[,2:6],
                                           S_DW_grass=sgrassSY,
                                           Y_DW_shrub=yearlyDW_shrub[,2:6],
                                           Y_ARI=yearlyARI[,2:6],
                                           Y_BSI=yearlyBSI[,2:6],
                                           Y_DeltaVH=yearlyDeltaVH[,2:6],
                                           sDPOLsummer=as.matrix(cbind(sDPOLsummer,sDPOLsummer,sDPOLsummer,
                                                                       sDPOLsummer,sDPOLsummer)),
                                           sNDMI=as.matrix(cbind(sNDMI,sNDMI,sNDMI,
                                                                 sNDMI,sNDMI)),
                                           Y_DPOLsummer=yearlyDPOLsummer[,2:6],
                                           Y_NBR=yearlyNBR[,2:6],
                                           sNBR=as.matrix(cbind(sNBR,sNBR,sNBR,
                                                                 sNBR,sNBR)),
                                           Y_NDMI=yearlyNDMI[,2:6],
                                           Y_NDVI=yearlyNDVI[,2:6],
                                           Y_NDWI=yearlyNDWI[,2:6],
                                           Y_REIP=yearlyREIP[,2:6]),
                       numPrimary=5) 


#Scale predictors and create nonlinear terms
umf@siteCovs$elev.s<-scale(umf@siteCovs$elev, center = TRUE, scale = TRUE)
umf@siteCovs$elev.s2<-umf@siteCovs$elev.s^2
umf@siteCovs$elev.sqrt<-sqrt(umf@siteCovs$elev)

umf@siteCovs$lat.s<-scale(umf@siteCovs$lat, center = TRUE, scale = TRUE)
umf@siteCovs$lat.s2<-umf@siteCovs$lat.s^2
umf@siteCovs$lat.sqrt<-sqrt(umf@siteCovs$lat)
umf@siteCovs$North<-ifelse(umf@siteCovs$lat<57,0,1)

umf@siteCovs$swi.s<-scale(umf@siteCovs$swi, center = TRUE, scale = TRUE)
umf@siteCovs$swi.s2<-umf@siteCovs$swi.s^2
umf@siteCovs$swi.sqrt<-sqrt(umf@siteCovs$swi)

umf@siteCovs$tpi.s<-scale(umf@siteCovs$tpi, center = TRUE, scale = TRUE)
umf@siteCovs$tpi.s2<-umf@siteCovs$tpi.s^2
umf@siteCovs$tpi.sqrt<-sqrt(umf@siteCovs$tpi)

umf@siteCovs$wet565.s<-scale(umf@siteCovs$wetlandProb_565m, center = TRUE, scale = TRUE)
umf@siteCovs$wet565.s2<-umf@siteCovs$wet565.s^2
umf@siteCovs$wet565.sqrt<-sqrt(umf@siteCovs$wetlandProb_565m)

#umf@siteCovs$DPOLsummer.s<-scale(umf@siteCovs$DPOLsummer, center = TRUE, scale = TRUE)[,1]
#umf@siteCovs$DPOLsummer.s2<-umf@siteCovs$DPOLsummer.s^2
#umf@siteCovs$DPOLsummer.sqrt<-sqrt(umf@siteCovs$DPOLsummer)

umf@siteCovs$sDPOLsummer.s<-scale(umf@siteCovs$sDPOLsummer, center = TRUE, scale = TRUE)
umf@siteCovs$sDPOLsummer.s2<-umf@siteCovs$sDPOLsummer.s^2
umf@siteCovs$sDPOLsummer.sqrt<-sqrt(umf@siteCovs$sDPOLsummer)

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

umf@siteCovs$sNBR.s<-scale(umf@siteCovs$sNBR, center = TRUE, scale = TRUE)
umf@siteCovs$sNBR.s2<-umf@siteCovs$sNBR.s^2
umf@siteCovs$sNBR.sqrt<-sqrt(umf@siteCovs$sNBR)

umf@siteCovs$REIP.s<-scale(umf@siteCovs$REIP, center = TRUE, scale = TRUE)
umf@siteCovs$REIP.s2<-umf@siteCovs$REIP.s^2
umf@siteCovs$REIP.sqrt<-sqrt(umf@siteCovs$REIP)

umf@siteCovs$NDMI.s<-scale(umf@siteCovs$NDMI, center = TRUE, scale = TRUE)
umf@siteCovs$NDMI.s2<-umf@siteCovs$NDMI.s^2
umf@siteCovs$NDMI.sqrt<-sqrt(umf@siteCovs$NDMI)

umf@siteCovs$sNDMI.s<-scale(umf@siteCovs$sNDMI, center = TRUE, scale = TRUE)
umf@siteCovs$sNDMI.s2<-umf@siteCovs$sNDMI.s^2
umf@siteCovs$sNDMI.sqrt<-sqrt(umf@siteCovs$sNDMI)

umf@siteCovs$NDVI.s<-scale(umf@siteCovs$NDVI, center = TRUE, scale = TRUE)
umf@siteCovs$NDVI.s2<-umf@siteCovs$NDVI.s^2
umf@siteCovs$NDVI.sqrt<-sqrt(umf@siteCovs$NDVI)

umf@siteCovs$NDWI.s<-scale(umf@siteCovs$NDWI, center = TRUE, scale = TRUE)
umf@siteCovs$NDWI.s2<-umf@siteCovs$NDWI.s^2
umf@siteCovs$NDWI.sqrt<-sqrt(umf@siteCovs$NDWI)

umf@siteCovs$bare.s<-scale(umf@siteCovs$bare, center = TRUE, scale = TRUE)
umf@siteCovs$bare.s2<-umf@siteCovs$bare.s^2
umf@siteCovs$bare.sqrt<-sqrt(umf@siteCovs$bare)

umf@siteCovs$grass.s<-scale(umf@siteCovs$grass, center = TRUE, scale = TRUE)
umf@siteCovs$grass.s2<-umf@siteCovs$grass.s^2
umf@siteCovs$grass.sqrt<-sqrt(umf@siteCovs$grass)

umf@siteCovs$sgrass.s<-scale(umf@siteCovs$sgrass, center = TRUE, scale = TRUE)
umf@siteCovs$sgrass.s2<-umf@siteCovs$sgrass.s^2
umf@siteCovs$sgrass.sqrt<-sqrt(umf@siteCovs$sgrass)

umf@siteCovs$shrub.s<-scale(umf@siteCovs$shrub, center = TRUE, scale = TRUE)
umf@siteCovs$shrub.s2<-umf@siteCovs$shrub.s^2
umf@siteCovs$shrub.sqrt<-sqrt(umf@siteCovs$shrub)

#umf@yearlySiteCovs$Y_distUpland.s<-scale(umf@yearlySiteCovs$Y_distUpland, center = TRUE, scale = TRUE)
#umf@yearlySiteCovs$Y_distUpland.s2<-umf@yearlySiteCovs$Y_distUpland.s^2
#umf@yearlySiteCovs$Y_distUpland.sqrt<-sqrt(umf@yearlySiteCovs$Y_distUpland)

#umf@yearlySiteCovs$Y_wet150.s<-scale(umf@yearlySiteCovs$Y_wet150, center = TRUE, scale = TRUE)
#umf@yearlySiteCovs$Y_wet150.s2<-umf@yearlySiteCovs$Y_wet150.s^2
#umf@yearlySiteCovs$Y_wet150.sqrt<-sqrt(umf@yearlySiteCovs$Y_wet150)

umf@yearlySiteCovs$Y_wet565.s<-scale(umf@yearlySiteCovs$Y_wet565, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_wet565.s2<-umf@yearlySiteCovs$Y_wet565.s^2
umf@yearlySiteCovs$Y_wet565.sqrt<-sqrt(umf@yearlySiteCovs$Y_wet565)

umf@yearlySiteCovs$Y_elev.s<-scale(umf@yearlySiteCovs$Y_elev, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_elev.s2<-umf@yearlySiteCovs$Y_elev.s^2
umf@yearlySiteCovs$Y_elev.sqrt<-sqrt(umf@yearlySiteCovs$Y_elev)

umf@yearlySiteCovs$Y_latitude.s<-scale(umf@yearlySiteCovs$Y_latitude, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_latitude.s2<-umf@yearlySiteCovs$Y_latitude.s^2
umf@yearlySiteCovs$Y_latitude.sqrt<-sqrt(umf@yearlySiteCovs$Y_latitude)
nlevels(as.factor(umf@yearlySiteCovs$Y_latitude))
range(umf@yearlySiteCovs$Y_latitude)#54.43446 57.69631
hist(umf@yearlySiteCovs$Y_latitude)
umf@yearlySiteCovs$North<-ifelse(umf@yearlySiteCovs$Y_latitude<57,0,1)

umf@yearlySiteCovs$Y_swi.s<-scale(umf@yearlySiteCovs$Y_swi, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_swi.s2<-umf@yearlySiteCovs$Y_swi.s^2
umf@yearlySiteCovs$Y_swi.sqrt<-sqrt(umf@yearlySiteCovs$Y_swi)

umf@yearlySiteCovs$Y_tpi.s<-scale(umf@yearlySiteCovs$Y_tpi, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_tpi.s2<-umf@yearlySiteCovs$Y_tpi.s^2
umf@yearlySiteCovs$Y_tpi.sqrt<-sqrt(umf@yearlySiteCovs$Y_tpi)

#umf@yearlySiteCovs$Y_mrvbf.s<-scale(umf@yearlySiteCovs$Y_mrvbf, center = TRUE, scale = TRUE)
#umf@yearlySiteCovs$Y_mrvbf.s2<-umf@yearlySiteCovs$Y_mrvbf.s^2
#umf@yearlySiteCovs$Y_mrvbf.sqrt<-sqrt(umf@yearlySiteCovs$Y_mrvbf)

umf@yearlySiteCovs$sDPOLsummer.s<-scale(umf@yearlySiteCovs$sDPOLsummer, center = TRUE, scale = TRUE)[,1]
umf@yearlySiteCovs$sDPOLsummer.s2<-umf@yearlySiteCovs$sDPOLsummer.s^2
umf@yearlySiteCovs$sDPOLsummer.sqrt<-sqrt(umf@yearlySiteCovs$sDPOLsummer)

umf@yearlySiteCovs$DeltaVH.s<-scale(umf@yearlySiteCovs$Y_DeltaVH, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$DeltaVH.s2<-umf@yearlySiteCovs$DeltaVH.s^2
umf@yearlySiteCovs$DeltaVH.sqrt<-sqrt(umf@yearlySiteCovs$Y_DeltaVH)

#umf@yearlySiteCovs$ARI.s<-scale(umf@yearlySiteCovs$Y_ARI, center = TRUE, scale = TRUE)[,1]
#umf@yearlySiteCovs$ARI.s2<-umf@yearlySiteCovs$ARI.s^2
#umf@yearlySiteCovs$ARI.sqrt<-sqrt(umf@yearlySiteCovs$Y_ARI)

umf@yearlySiteCovs$BSI.s<-scale(umf@yearlySiteCovs$Y_BSI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$BSI.s2<-umf@yearlySiteCovs$BSI.s^2
umf@yearlySiteCovs$BSI.sqrt<-sqrt(umf@yearlySiteCovs$Y_BSI)

umf@yearlySiteCovs$NBR.s<-scale(umf@yearlySiteCovs$Y_NBR, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$NBR.s2<-umf@yearlySiteCovs$NBR.s^2
umf@yearlySiteCovs$NBR.sqrt<-sqrt(umf@yearlySiteCovs$Y_NBR)

umf@yearlySiteCovs$sNBR.s<-scale(umf@yearlySiteCovs$sNBR, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$sNBR.s2<-umf@yearlySiteCovs$sNBR.s^2
umf@yearlySiteCovs$sNBR.sqrt<-sqrt(umf@yearlySiteCovs$sNBR)

umf@yearlySiteCovs$REIP.s<-scale(umf@yearlySiteCovs$Y_REIP, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$REIP.s2<-umf@yearlySiteCovs$REIP.s^2
umf@yearlySiteCovs$REIP.sqrt<-sqrt(umf@yearlySiteCovs$Y_REIP)

umf@yearlySiteCovs$NDMI.s<-scale(umf@yearlySiteCovs$Y_NDMI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$NDMI.s2<-umf@yearlySiteCovs$NDMI.s^2
umf@yearlySiteCovs$NDMI.sqrt<-sqrt(umf@yearlySiteCovs$Y_NDMI)

umf@yearlySiteCovs$sNDMI.s<-scale(umf@yearlySiteCovs$sNDMI, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$sNDMI.s2<-umf@yearlySiteCovs$sNDMI.s^2
umf@yearlySiteCovs$sNDMI.sqrt<-sqrt(umf@yearlySiteCovs$sNDMI)

umf@yearlySiteCovs$bare.s<-scale(umf@yearlySiteCovs$Y_DW_bare, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$bare.s2<-umf@yearlySiteCovs$bare.s^2
umf@yearlySiteCovs$bare.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_bare)

umf@yearlySiteCovs$grass.s<-scale(umf@yearlySiteCovs$Y_DW_grass, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$grass.s2<-umf@yearlySiteCovs$grass.s^2
umf@yearlySiteCovs$grass.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_grass)

umf@yearlySiteCovs$sgrass.s<-scale(umf@yearlySiteCovs$S_DW_grass, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$sgrass.s2<-umf@yearlySiteCovs$sgrass.s^2
umf@yearlySiteCovs$sgrass.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_grass)

umf@yearlySiteCovs$shrub.s<-scale(umf@yearlySiteCovs$Y_DW_shrub, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$shrub.s2<-umf@yearlySiteCovs$shrub.s^2
umf@yearlySiteCovs$shrub.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_shrub)


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

mb.gof.test(model1, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
#87.0703 104.1541  35.4833  31.4974  61.1523 
#Total chi-square = 319.3575 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#112  155  169  187  315 
#Estimate of c-hat = 1.85 

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
RSonlyModel<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+DeltaVH.s2+sNBR.s+bare.s+sgrass.s+sgrass.s2+lat.sqrt, #either keep NDMI or NBR
                    gammaformula =  ~ sDPOLsummer.s+Y_swi.sqrt, 
                    epsilonformula = ~ NDMI.s+DeltaVH.s+DeltaVH.s2+sDPOLsummer.s+Y_tpi.s+Y_latitude.sqrt+shrub.s,
                    pformula = ~ 1, umf)#AIC: 1576.832
#Drop eps ~ NDMI, TPI, shrub, gamma ~ DPOL, psi ~ NBR, DeltaVHs2, bare

#square-root function of latitude
RSonlyModel<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+sgrass.s+sgrass.s2+lat.sqrt, #either keep NDMI or NBR
                    gammaformula =  ~ Y_swi.sqrt, 
                    epsilonformula = ~ DeltaVH.s+DeltaVH.s2+sDPOLsummer.s+Y_latitude.sqrt,
                    pformula = ~ 1, umf)#AIC: 1575.717

mb.gof.test(RSonlyModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#87.7963  87.0103  35.1685  32.5430  62.2708 
#Total chi-square = 304.7889 
#Number of bootstrap samples = 1000
#P-value = 0.009
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#121  154  169  187  451 
#Estimate of c-hat = 1.74 

#linear function of latitude
RSonlyModelB<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+sgrass.s+sgrass.s2+lat.s, #either keep NDMI or NBR
                    gammaformula =  ~ Y_swi.sqrt, 
                    epsilonformula = ~ DeltaVH.s+DeltaVH.s2+sDPOLsummer.s+Y_latitude.s,
                    pformula = ~ 1, umf)
#AIC: 1575.505 -effect of latitude on psi +effect on extinction

mb.gof.test(RSonlyModelB, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#87.8395  87.6681  35.1417  32.8354  62.1547 
#Total chi-square = 305.6394 
#Number of bootstrap samples = 1000
#P-value = 0.005
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#114  153  167  182  457 
#Estimate of c-hat = 1.78 

#binary function of latitude
RSonlyModelC<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+sgrass.s+sgrass.s2+North, #either keep NDMI or NBR
                     gammaformula =  ~ Y_swi.sqrt, 
                     epsilonformula = ~ DeltaVH.s+DeltaVH.s2+sDPOLsummer.s+North,
                     pformula = ~ 1, umf)#AIC: 1578.716

mb.gof.test(RSonlyModelC, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#87.9167  85.6383  35.1850  33.1167  61.8632 
#Total chi-square = 303.72 
#Number of bootstrap samples = 1000
#P-value = 0.007
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#114  154  168  186  384 
#Estimate of c-hat = 1.74 

NullModel<-colext(psiformula= ~1, 
                    gammaformula =  ~ 1, 
                    epsilonformula = ~ 1,
                    pformula = ~ 1, umf)#AIC: 1690.285

mb.gof.test(NullModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#86.3005 101.9738  36.0336  32.9264  62.6683 
#Total chi-square = 319.9025 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#121  158  173  188  262 
#Estimate of c-hat = 1.84 

#"Full Model"
FullModel<-colext(psiformula= ~HabSuit4ClassNumeric+NDMI.s+NDMI.s2+DeltaVH.s+sgrass.s+sgrass.s2+lat.sqrt, 
                  gammaformula =  ~ Y_swi.sqrt, 
                  epsilonformula = ~ Y_habSuit+DeltaVH.s+DeltaVH.s2+sDPOLsummer.s+Y_latitude.sqrt,
                  pformula = ~ 1, umf)#AIC: 1570.577
#Drop psi ~ HabSuit4ClassNumeric?

FullModel<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+sgrass.s+sgrass.s2+lat.sqrt, 
                  gammaformula =  ~ Y_swi.sqrt, 
                  epsilonformula = ~ Y_habSuit+DeltaVH.s+DeltaVH.s2+sDPOLsummer.s+Y_latitude.sqrt,
                  pformula = ~ 1, umf)#AIC: 1571.606

#Call:
#colext(psiformula = ~NDMI.s + NDMI.s2 + DeltaVH.s + sgrass.s + 
#         sgrass.s2 + lat.sqrt, gammaformula = ~Y_swi.sqrt, epsilonformula = ~Y_habSuit + 
#         DeltaVH.s + DeltaVH.s2 + sDPOLsummer.s + Y_latitude.sqrt, 
#       pformula = ~1, data = umf)

#Initial (logit-scale):
#  Estimate     SE     z P(>|z|)
#(Intercept)    98.92 51.172  1.93  0.0532
#NDMI.s         -1.09  0.567 -1.92  0.0549
#NDMI.s2        -1.45  0.581 -2.49  0.0126
#DeltaVH.s      -1.12  0.491 -2.28  0.0226
#sgrass.s        2.79  1.204  2.31  0.0206
#sgrass.s2       3.29  1.344  2.45  0.0144
#lat.sqrt      -13.09  6.779 -1.93  0.0535

#Colonization (logit-scale):
#  Estimate   SE     z  P(>|z|)
#(Intercept)    -17.1 5.19 -3.30 0.000963
#Y_swi.sqrt       5.2 1.66  3.14 0.001688

#Extinction (logit-scale):
#  Estimate     SE      z  P(>|z|)
#(Intercept)      -66.139 26.101 -2.534 0.011277
#Y_habSuit2         1.330  0.726  1.832 0.066917
#Y_habSuit3         1.722  0.664  2.592 0.009545
#Y_habSuit4         0.226  0.677  0.333 0.738786
#DeltaVH.s          0.652  0.300  2.173 0.029802
#DeltaVH.s2         1.014  0.279  3.631 0.000282
#sDPOLsummer.s      1.391  0.398  3.500 0.000465
#Y_latitude.sqrt    8.602  3.451  2.492 0.012693

#Detection (logit-scale):
#  Estimate     SE    z  P(>|z|)
#0.472 0.0727 6.49 8.46e-11

#AIC: 1571.606 
#Number of sites: 115

mb.gof.test(FullModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
#88.0588  88.9277  34.9721  31.8075  62.6698 
#Total chi-square = 306.4358 
#Number of bootstrap samples = 1000
#P-value = 0.003
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#113  153  168  186  391 
#Estimate of c-hat = 1.76  

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
ggsave(filename="2_outputs/1_unmarked_outputs/November/Hedley 2016 2021/occu.ByRegion2016_2021_HedleyXHybridRS.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)

FullModelWithoutLat<-colext(psiformula= ~NDMI.s+NDMI.s2+DeltaVH.s+sgrass.s+sgrass.s2, 
                  gammaformula =  ~ Y_swi.sqrt, 
                  epsilonformula = ~ Y_habSuit+DeltaVH.s+DeltaVH.s2+sDPOLsummer.s,
                  pformula = ~ 1, umf)#AIC: 1579.463
#Initial (logit-scale):
#Estimate    SE      z  P(>|z|)
#(Intercept)    0.211 0.691  0.305 0.760550
#NDMI.s        -0.690 0.440 -1.568 0.116853
#NDMI.s2       -1.049 0.464 -2.262 0.023683
#DeltaVH.s     -1.562 0.431 -3.621 0.000293
#sgrass.s       2.494 1.016  2.455 0.014071
#sgrass.s2      3.299 1.299  2.540 0.011074

#Colonization (logit-scale):
#  Estimate   SE     z P(>|z|)
#(Intercept)   -16.69 5.15 -3.24 0.00120
#Y_swi.sqrt      5.06 1.65  3.07 0.00211

#Extinction (logit-scale):
#  Estimate    SE      z  P(>|z|)
#(Intercept)     -1.428 0.493 -2.895 0.003786
#Y_habSuit2       1.502 0.696  2.158 0.030924
#Y_habSuit3       2.014 0.666  3.025 0.002485
#Y_habSuit4       0.340 0.668  0.509 0.610758
#DeltaVH.s        0.871 0.285  3.054 0.002258
#DeltaVH.s2       0.988 0.255  3.870 0.000109
#sDPOLsummer.s    1.524 0.411  3.710 0.000207

#Detection (logit-scale):
#  Estimate     SE    z  P(>|z|)
#0.478 0.0723 6.62 3.69e-11
mb.gof.test(FullModelWithoutLat, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
#88.0585  91.0656  35.1793  33.1377  62.2055 
#Total chi-square = 309.6466 
#Number of bootstrap samples = 1000
#P-value = 0.009
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#111  153  167  183  413 
#Estimate of c-hat = 1.8

#Get result in JAGS

#cluster ID for random effect
clusterID<-as.integer(umf@siteCovs$clusterID)
#reconfigure site and site-year covariates as dummy variable matrices
class2<-umf@siteCovs$HabSuit4ClassNumeric
class2<-ifelse(class2==2,1,0)
class3<-umf@siteCovs$HabSuit4ClassNumeric
class3<-ifelse(class3==3,1,0)
class4<-umf@siteCovs$HabSuit4ClassNumeric
class4<-ifelse(class4==4,1,0)

#class2SY<-umf@yearlySiteCovs$habSuit
#class2SY<-ifelse(class2SY==2,1,0)
#class3SY<-umf@yearlySiteCovs$habSuit
#class3SY<-ifelse(class3SY==2,1,0)
#class4SY<-umf@yearlySiteCovs$habSuit
#class4SY<-ifelse(class4SY==2,1,0)

class2SY<-as.matrix(cbind(class2,class2,class2,class2,class2))
class3SY<-as.matrix(cbind(class3,class3,class3,class3,class3))
class4SY<-as.matrix(cbind(class4,class4,class4,class4,class4))

#make y into a 3-dimensional array
obsFiltFinal.jags<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC,yera2016_1,yera2016_2,yera2016_3,yera2016_4,yera2016_5,
         yera2017_1,yera2017_2,yera2017_3,yera2017_4,yera2017_5,
         yera2018_1,yera2018_2,yera2018_3,yera2018_4,yera2018_5,
         yera2019_1,yera2019_2,yera2019_3,yera2019_4,yera2019_5,
         yera2021_1,yera2021_2,yera2021_3,yera2021_4,yera2021_5)
obsFiltFinal.jags$geometry<-NULL

detectionsLong<-obsFiltFinal.jags%>%
  pivot_longer(cols=yera2016_1:yera2021_5, names_to="year_round", values_to="detected")%>%
  mutate(year=substr(year_round,5,8))%>%
  mutate(round=substr(year_round,10,10))%>%
  arrange(locationC,year,round)
nlevels(as.factor(detectionsLong$locationC))
#115 sites

#get detections and put in 3-d array
R <- nlevels(as.factor(detectionsLong$locationC))#146 sites
T <- nlevels(as.factor(detectionsLong$round))#up to 5 visits per site per year
S <- nlevels(as.factor(detectionsLong$year))#8 years (no visits to any sites in 2020 due to COVID)

nsite<-R
nrep<-T
nyear<-S

sites<-levels(as.factor(detectionsLong$locationC))
rounds<-levels(as.factor(detectionsLong$round))
years<-levels(as.factor(detectionsLong$year))

y <-array(NA,
          dim=c(nsite,nrep,nyear),
          dimnames=list(sites, rounds, years)
)
for (i in sites){
  for (j in rounds){
    for (k in years){
      try(placeHolder<-detectionsLong%>%
            filter(locationC==i)%>%
            filter(round==j)%>%
            filter(year==k)%>%
            select(detected))
      try(y[i,j,k]<-placeHolder$detected) 
    }
  }
}
str(y)
#NO DETECTION SUBMODEL COVARIATES IN THIS ANALYSIS
#THEREFORE NO DETECTION ARRAYS NEEDED

#INCLUDING LATITUDE
#####################################################################################################################
#                                                                                                                   #
#                                                                                                                   #
#                       Create Model File and Initialize Parameters - Modelling Individual Footprint Types          #
#                                                                                                                   #
#                                                                                                                   #
#                                                                                                                   #
#                                                                                                                   #
#####################################################################################################################

ndmi<-umf@siteCovs$sNDMI.s[,1]
ndmi2<-umf@siteCovs$sNDMI.s2[,1]
deltaVH<-umf@siteCovs$DeltaVH.s[,1]
sgrass<-umf@siteCovs$sgrass.s[,1]
sgrass2<-umf@siteCovs$sgrass.s2[,1]
latsqrt<-umf@siteCovs$lat.sqrt

Yswi<-scale(Y_swi, center=TRUE, scale=TRUE)
YdeltaVH<-scale(yearlyDeltaVH[,2:6], center=TRUE, scale=TRUE)
YdeltaVH2<-YdeltaVH^2

sDPOL<-as.matrix(cbind(sDPOLsummer,sDPOLsummer,sDPOLsummer,
                            sDPOLsummer,sDPOLsummer))
YsDPOLsummer<-scale(sDPOL, center=TRUE, scale=TRUE)

Ylatsqrt<-scale(sqrt(yearlyLatitude), center=TRUE, scale=TRUE)

#Run multi-season occu. model in JAGS
# Specify model in JAGS language
model_code <-
  "model
{
  # Specify priors
  alpha.psi ~ dunif(-3, 3) 
  beta1.psi ~ dunif(-3, 3)
  beta2.psi ~ dunif(-3, 3)
  beta3.psi ~ dunif(-3, 3) 
  beta4.psi ~ dunif(-3, 3)
  beta5.psi ~ dunif(-3, 3)
  beta6.psi ~ dunif(-3, 3) 

  alpha.gamma ~ dunif(-3, 3)
  beta1.gamma ~ dunif(-3, 3) 

  alpha.phi ~ dunif(-3, 3)  
  beta1.phi ~ dunif(-3, 3)
  beta2.phi ~ dunif(-3, 3)
  beta3.phi ~ dunif(-3, 3)
  beta4.phi ~ dunif(-3, 3)
  beta5.phi ~ dunif(-3, 3)
  beta6.phi ~ dunif(-3, 3)
  beta7.phi ~ dunif(-3, 3)

  alpha.p ~ dunif(-3, 3) #no detection covariates for now
  
  #random effects on psi, gamma and phi 
  #115 sites, 7 cluster IDs
  for (m in 1:ncluster) {
    u_psi[m] ~ dnorm(0, tau_psi)
    u_phi[m] ~ dnorm(0, tau_phi) #Wait, why dnorm and not dunif?
    u_gamma[m] ~ dnorm(0, tau_gamma)
  }
  tau_psi <- pow(sigma_psi, -2)
  sigma_psi ~ dunif(0, 5)

  tau_phi <- pow(sigma_phi, -2)
  sigma_phi ~ dunif(0, 5)

  tau_gamma <- pow(sigma_gamma, -2)
  sigma_gamma ~ dunif(0, 5)
  

  #assign values generated for u_phi and u_gamma to a vector of size 93 based
  #on clusterID?


  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){
    z0[i] ~dbern(psi1[i]) #True occu. at site i in season 1
    logit(psi1[i])<-alpha.psi +beta1.psi*ndmi[i] +beta2.psi*ndmi2[i] +beta3.psi*deltaVH[i] +beta4.psi*sgrass[i] +beta5.psi*sgrass2[i] +beta6.psi*latsqrt[i] + u_psi[cluster[i]]

    z[i,1]<-z0[i]
    for (k in 2:nyear){
      #survival probability from season k-1 to season k
      phi[i,k-1] <-psurv[i,k-1]
      logit(psurv[i,k-1])<-alpha.phi + u_phi[cluster[i]] +beta1.phi*class2SY[i,k-1] +beta2.phi*class3SY[i,k-1] +beta3.phi*class4SY[i,k-1] +beta4.phi*YdeltaVH[i,k-1] +beta5.phi*YdeltaVH2[i,k-1] +beta6.phi*YsDPOLsummer[i,k-1] +beta7.phi*Ylatsqrt[i,k-1] 

      #colonization probability from season k-1 to season k
      gamma[i,k-1] <-pcol[i,k-1]
      logit(pcol[i,k-1])<-alpha.gamma  +u_gamma[cluster[i]] +beta1.gamma*Yswi[i,k-1] 

      muZ[i,k]<-z[i,k-1]*phi[i,k-1] +(1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~dbern(muZ[i,k])
    } #k
  } #i
  # Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (k in 1:nyear){
        p[i,j,k]<-det[i,j,k]
        logit(det[i,j,k])<-alpha.p    #(no detection covariates for now)    
        muy[i,j,k]<-z[i,k]*p[i,j,k]
        y[i,j,k]~dbern(muy[i,j,k])
        log_likelihood_contrib[i,j,k]<-y[i,j,k]*log(muy[i,j,k]) + (1-y[i,j,k])*log(1-muy[i,j,k])
      } #k
    } #j
  } #i
  # Derived parameters: Sample and population occu., growth rate and turnover
  n.occ[1]<-sum(z[1:nsite,1])
  for (i in 1:nsite){
    psi[i,1] <-psi1[i]
    for (k in 2:nyear){
      psi[i,k] <-psi[i,k-1]*phi[i,k-1] +(1-psi[i,k-1])*gamma[i,k-1]
    }
  }
  for (k in 2:nyear){
      n.occ[k] <-sum(z[1:nsite,k])
      growthr[k-1]<-mean(psi[1:nsite,k-1])/mean(psi[1:nsite,k-1])
      #no longer growthr[k]<-psi[k]/psi[k-1] or psi[i,k]/psi[i,k-1]
      turnover[k-1] <-(1 - mean(psi[1:nsite,k-1])) *mean(gamma[1:nsite,k-1])/mean(psi[1:nsite,k])
      #turnover[k-1] <-(1 - psi[i,k-1]) *gamma[k-1]/psi[i,k]
  }
  
}"
  
  
# Bundle data
nsite<-115
nyear<-5
nrep<-5
data_jags <-list(y=y,
                   #thr=thr,
                   nsite=nsite,
                   nrep=nrep,
                   nyear=nyear,
                   #psi covariates
                   ndmi=ndmi,
                   ndmi2=ndmi2,
                   deltaVH=deltaVH,
                   sgrass=sgrass,
                   sgrass2=sgrass2,
                   latsqrt=latsqrt,
                   #site-year covariates (gamma, phi)
                   class2SY=class2SY,
                   class3SY=class3SY,
                   class4SY=class4SY,
                   Yswi=Yswi,
                   YdeltaVH=YdeltaVH,
                   YdeltaVH2=YdeltaVH2,
                   YsDPOLsummer=YsDPOLsummer,
                   Ylatsqrt=Ylatsqrt,
                   #cluster information
                   ncluster=7,
                   cluster=clusterID)
#no detection covariates (p)
  
# Initial values
zst <-apply(y,c(1,3),max,na.rm=TRUE)
#Observed occurrence as inits for z
zst[zst==-Inf]<-0
  
#Run model via parallel processing:
initsP = parallel.seeds("base::BaseRNG", 3) # a list of lists
initsP[[1]]$z=zst
initsP[[1]]$alpha.psi=runif(1, -3, 3)
initsP[[1]]$beta1.psi=runif(1, -3, 3)
initsP[[1]]$beta2.psi=runif(1, -3, 3)
initsP[[1]]$beta3.psi=runif(1, -3, 3)
initsP[[1]]$beta4.psi=runif(1, -3, 3)
initsP[[1]]$beta5.psi=runif(1, -3, 3)
initsP[[1]]$beta6.psi=runif(1, -3, 3)
initsP[[1]]$u_psi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.phi=runif(1, -3, 3)
initsP[[1]]$beta1.phi =runif(1, -3, 3)
initsP[[1]]$beta2.phi =runif(1, -3, 3)
initsP[[1]]$beta3.phi =runif(1, -3, 3)
initsP[[1]]$beta4.phi =runif(1, -3, 3)
initsP[[1]]$beta5.phi =runif(1, -3, 3)
initsP[[1]]$beta6.phi =runif(1, -3, 3)
initsP[[1]]$beta7.phi =runif(1, -3, 3)
initsP[[1]]$u_phi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.gamma=runif(1, -3, 3)
initsP[[1]]$beta1.gamma =runif(1, -3, 3)
initsP[[1]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.p=runif(1,-3,3)
  
initsP[[2]]$z=zst
initsP[[2]]$alpha.psi=runif(1, -3, 3)
initsP[[2]]$beta1.psi=runif(1, -3, 3)
initsP[[2]]$beta2.psi=runif(1, -3, 3)
initsP[[2]]$beta3.psi=runif(1, -3, 3)
initsP[[2]]$beta4.psi=runif(1, -3, 3)
initsP[[2]]$beta5.psi=runif(1, -3, 3)
initsP[[2]]$beta6.psi=runif(1, -3, 3)
initsP[[2]]$u_psi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.phi=runif(1, -3, 3)
initsP[[2]]$beta1.phi =runif(1, -3, 3)
initsP[[2]]$beta2.phi =runif(1, -3, 3)
initsP[[2]]$beta3.phi =runif(1, -3, 3)
initsP[[2]]$beta4.phi =runif(1, -3, 3)
initsP[[2]]$beta5.phi =runif(1, -3, 3)
initsP[[2]]$beta6.phi =runif(1, -3, 3)
initsP[[2]]$beta7.phi =runif(1, -3, 3)
initsP[[2]]$u_phi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.gamma=runif(1, -3, 3)
initsP[[2]]$beta1.gamma =runif(1, -3, 3)
initsP[[2]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.p=runif(1,-3,3)
  
initsP[[3]]$z=zst
initsP[[3]]$alpha.psi=runif(1, -3, 3)
initsP[[3]]$beta1.psi=runif(1, -3, 3)
initsP[[3]]$beta2.psi=runif(1, -3, 3)
initsP[[3]]$beta3.psi=runif(1, -3, 3)
initsP[[3]]$beta4.psi=runif(1, -3, 3)
initsP[[3]]$beta5.psi=runif(1, -3, 3)
initsP[[3]]$beta6.psi=runif(1, -3, 3)
initsP[[3]]$u_psi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.phi=runif(1, -3, 3)
initsP[[3]]$beta1.phi =runif(1, -3, 3)
initsP[[3]]$beta2.phi =runif(1, -3, 3)
initsP[[3]]$beta3.phi =runif(1, -3, 3)
initsP[[3]]$beta4.phi =runif(1, -3, 3)
initsP[[3]]$beta5.phi =runif(1, -3, 3)
initsP[[3]]$beta6.phi =runif(1, -3, 3)
initsP[[3]]$beta7.phi =runif(1, -3, 3)
initsP[[3]]$u_phi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.gamma=runif(1, -3, 3)
initsP[[3]]$beta1.gamma =runif(1, -3, 3)
initsP[[3]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.p=runif(1,-3,3)
  
# Parameters monitored
params <-c("psi",
             "alpha.psi",
             "beta1.psi",
             "beta2.psi",
             "beta3.psi",
             "beta4.psi",
             "beta5.psi",
             "beta6.psi",
             "u_psi",
             "phi",
             "alpha.phi",
             "beta1.phi",
             "beta2.phi",
             "beta3.phi",
             "beta4.phi",
             "beta5.phi",
             "beta6.phi",
             "beta7.phi",
             "u_phi",
             "gamma",
             "alpha.gamma",
             "beta1.gamma",
             "u_gamma",
             "alpha.p",
             "p","z","muZ",
             "log_likelihood_contrib",
             "n.occ",
             "growthr",
             "turnover")#,"sens","fpr","auc"
  
  
# Running the model via parallel processing using runjags
ni <-60000#instead of 40000
nt <-10
nb <-50000#instead of 30000
nc <-3
#so there should be 1000 posterior draws
  
niter<-ni
nchains<-nc
thinning<-nt
  
results <- runjags::run.jags(
    model = model_code,
    monitor = params,#keep_track_of,
    data = data_jags,#datalist
    n.chains = nchains,
    inits=initsP,
    adapt = niter/10,
    burnin = nb, 
    sample = niter/thinning, 
    thin = thinning,
    method = "parallel"
  )
cat("Model run completed\n")
print(Sys.time())
  
# Extracting results
posterior_draws <- as.mcmc.list(results)
cat("Posterior draws completed\n")
smod <- summary(posterior_draws)
  
save(detectionsLong, 
       data_jags, 
       results, 
       posterior_draws, 
       smod, 
       nc, 
       ni, 
       nb, 
       nt, 
       params, 
       file = "2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_HedleyPlusHybridWithLat.Rdata")
  
load("2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_HedleyPlusHybridWithLat.Rdata")
# Extract posterior samples
#Process the MCMC output:
int.mcmc.mat <- as.matrix(posterior_draws)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)
names(int.mcmc.dat)
int.mcmc.dat[,"phi[64,1]"]#18000 estimates for this parameter
  
  
# Diagnostics:                                        Point est. Upper C.I.
gelman.diag(posterior_draws[1:3][,"alpha.psi"])#[1,]          1       1       
gelman.diag(posterior_draws[1:3][,"beta1.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta2.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta3.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta4.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta5.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta6.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.phi"])#[1,]          1          1       
gelman.diag(posterior_draws[1:3][,"beta1.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta2.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta3.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta4.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta5.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta6.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta7.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.gamma"])#[1,]          1          1       
gelman.diag(posterior_draws[1:3][,"beta1.gamma"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.p"])#[1,]          1          1           
  
range(
    effectiveSize(posterior_draws[1:3][,"alpha.psi"]),        
    effectiveSize(posterior_draws[1:3][,"beta1.psi"]),           
    effectiveSize(posterior_draws[1:3][,"beta2.psi"]),           
    effectiveSize(posterior_draws[1:3][,"beta3.psi"]),           
    effectiveSize(posterior_draws[1:3][,"beta4.psi"]),           
    effectiveSize(posterior_draws[1:3][,"beta5.psi"]),           
    effectiveSize(posterior_draws[1:3][,"beta6.psi"]),           
    effectiveSize(posterior_draws[1:3][,"alpha.phi"]),        
    effectiveSize(posterior_draws[1:3][,"beta1.phi"]),           
    effectiveSize(posterior_draws[1:3][,"beta2.phi"]),           
    effectiveSize(posterior_draws[1:3][,"beta3.phi"]),           
    effectiveSize(posterior_draws[1:3][,"beta4.phi"]),           
    effectiveSize(posterior_draws[1:3][,"beta5.phi"]),           
    effectiveSize(posterior_draws[1:3][,"beta6.phi"]),           
    effectiveSize(posterior_draws[1:3][,"beta7.phi"]),           
    effectiveSize(posterior_draws[1:3][,"alpha.gamma"]),        
    effectiveSize(posterior_draws[1:3][,"beta1.gamma"])           
  )
#effective sample size varies from 1120.647 16351.098
  
heidel.diag(posterior_draws[1:3][,"alpha.psi"])#passed first, failed second        
heidel.diag(posterior_draws[1:3][,"beta1.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta2.psi"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"beta3.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta4.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta5.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta6.psi"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"alpha.phi"])#passed both        
heidel.diag(posterior_draws[1:3][,"beta1.phi"])#passed both          
heidel.diag(posterior_draws[1:3][,"beta2.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta3.phi"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"beta4.phi"])#passed both          
heidel.diag(posterior_draws[1:3][,"beta5.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta6.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta7.phi"])#passed both          
heidel.diag(posterior_draws[1:3][,"alpha.gamma"])#passed both        
heidel.diag(posterior_draws[1:3][,"beta1.gamma"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.p"])#passed both           
#All passed the stationary test (p>0.1), failed the halfwidth test
  
gd <- gelman.diag(posterior_draws[1:3][,"alpha.psi"])
es <- effectiveSize(posterior_draws[1:3][,"alpha.psi"])#   
hd <- heidel.diag(posterior_draws[1:3][,"alpha.psi"])#
  
#Trace Plots and Density Plots
plot(posterior_draws[1:3][,"beta1.psi"],density=FALSE)#- effect on psi
densplot(posterior_draws[[1]][,"beta1.psi"], main="beta1.psi (NDMI) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)
  
plot(posterior_draws[1:3][,"beta2.psi"],density=FALSE)#nonsig effect on psi
densplot(posterior_draws[[1]][,"beta2.psi"], main="beta2.psi (NDMI2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)
  
plot(posterior_draws[1:3][,"beta3.psi"],density=FALSE)#- effect on psi
densplot(posterior_draws[[1]][,"beta3.psi"], main="beta3.psi (DeltaVH) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)
  
plot(posterior_draws[1:3][,"beta4.psi"],density=FALSE)#+ effect on psi
densplot(posterior_draws[[1]][,"beta4.psi"], main="beta4.psi (sgrass) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta5.psi"],density=FALSE)#+ effect on psi
densplot(posterior_draws[[1]][,"beta5.psi"], main="beta5.psi (sgrass2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta6.psi"],density=FALSE)#nonsig effect on psi
densplot(posterior_draws[[1]][,"beta6.psi"], main="beta6.psi (latsqrt) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta1.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta1.phi"], main="beta1.phi (habitat category 2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)
  
plot(posterior_draws[1:3][,"beta2.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta2.phi"], main="beta2.phi (habitat category 3) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)
  
plot(posterior_draws[1:3][,"beta3.phi"],density=FALSE)#nonsig effect on phi
densplot(posterior_draws[[1]][,"beta3.phi"], main="beta3.phi (habitat category 4) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta4.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta4.phi"], main="beta4.phi (DeltaVH) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta5.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta5.phi"], main="beta5.phi (DeltaVH2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta6.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta6.phi"], main="beta6.phi (sDPOLsummer) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta7.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta7.phi"], main="beta7.phi (latsqrt) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta1.gamma"],density=FALSE)#+ effect on gamma
densplot(posterior_draws[[1]][,"beta1.gamma"], main="beta1.gamma (SWI) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

#Extract Model Coefficients
alpha.psi.PostMedian<-smod$quantiles["alpha.psi","50%"]
alpha.psi.Post.2.5<-smod$quantiles["alpha.psi","2.5%"]
alpha.psi.Post.97.5<-smod$quantiles["alpha.psi","97.5%"]

beta1.psi.PostMedian<-smod$quantiles["beta1.psi","50%"]
beta1.psi.Post.2.5<-smod$quantiles["beta1.psi","2.5%"]
beta1.psi.Post.97.5<-smod$quantiles["beta1.psi","97.5%"]

beta2.psi.PostMedian<-smod$quantiles["beta2.psi","50%"]
beta2.psi.Post.2.5<-smod$quantiles["beta2.psi","2.5%"]
beta2.psi.Post.97.5<-smod$quantiles["beta2.psi","97.5%"]

beta3.psi.PostMedian<-smod$quantiles["beta3.psi","50%"]
beta3.psi.Post.2.5<-smod$quantiles["beta3.psi","2.5%"]
beta3.psi.Post.97.5<-smod$quantiles["beta3.psi","97.5%"]

beta4.psi.PostMedian<-smod$quantiles["beta4.psi","50%"]
beta4.psi.Post.2.5<-smod$quantiles["beta4.psi","2.5%"]
beta4.psi.Post.97.5<-smod$quantiles["beta4.psi","97.5%"]

beta5.psi.PostMedian<-smod$quantiles["beta5.psi","50%"]
beta5.psi.Post.2.5<-smod$quantiles["beta5.psi","2.5%"]
beta5.psi.Post.97.5<-smod$quantiles["beta5.psi","97.5%"]

beta6.psi.PostMedian<-smod$quantiles["beta6.psi","50%"]
beta6.psi.Post.2.5<-smod$quantiles["beta6.psi","2.5%"]
beta6.psi.Post.97.5<-smod$quantiles["beta6.psi","97.5%"]

alpha.phi.PostMedian<-smod$quantiles["alpha.phi","50%"]
alpha.phi.Post.2.5<-smod$quantiles["alpha.phi","2.5%"]
alpha.phi.Post.97.5<-smod$quantiles["alpha.phi","97.5%"]

beta1.phi.PostMedian<-smod$quantiles["beta1.phi","50%"]
beta1.phi.Post.2.5<-smod$quantiles["beta1.phi","2.5%"]
beta1.phi.Post.97.5<-smod$quantiles["beta1.phi","97.5%"]

beta2.phi.PostMedian<-smod$quantiles["beta2.phi","50%"]
beta2.phi.Post.2.5<-smod$quantiles["beta2.phi","2.5%"]
beta2.phi.Post.97.5<-smod$quantiles["beta2.phi","97.5%"]

beta3.phi.PostMedian<-smod$quantiles["beta3.phi","50%"]
beta3.phi.Post.2.5<-smod$quantiles["beta3.phi","2.5%"]
beta3.phi.Post.97.5<-smod$quantiles["beta3.phi","97.5%"]

beta4.phi.PostMedian<-smod$quantiles["beta4.phi","50%"]
beta4.phi.Post.2.5<-smod$quantiles["beta4.phi","2.5%"]
beta4.phi.Post.97.5<-smod$quantiles["beta4.phi","97.5%"]

beta5.phi.PostMedian<-smod$quantiles["beta5.phi","50%"]
beta5.phi.Post.2.5<-smod$quantiles["beta5.phi","2.5%"]
beta5.phi.Post.97.5<-smod$quantiles["beta5.phi","97.5%"]

beta6.phi.PostMedian<-smod$quantiles["beta6.phi","50%"]
beta6.phi.Post.2.5<-smod$quantiles["beta6.phi","2.5%"]
beta6.phi.Post.97.5<-smod$quantiles["beta6.phi","97.5%"]

beta7.phi.PostMedian<-smod$quantiles["beta7.phi","50%"]
beta7.phi.Post.2.5<-smod$quantiles["beta7.phi","2.5%"]
beta7.phi.Post.97.5<-smod$quantiles["beta7.phi","97.5%"]

alpha.gamma.PostMedian<-smod$quantiles["alpha.gamma","50%"]
alpha.gamma.Post.2.5<-smod$quantiles["alpha.gamma","2.5%"]
alpha.gamma.Post.97.5<-smod$quantiles["alpha.gamma","97.5%"]

beta1.gamma.PostMedian<-smod$quantiles["beta1.gamma","50%"]
beta1.gamma.Post.2.5<-smod$quantiles["beta1.gamma","2.5%"]
beta1.gamma.Post.97.5<-smod$quantiles["beta1.gamma","97.5%"]

#create ggridges plots
# Combine the results in a single data frame: gamma
post_sim_gamma <- data.frame(samples = c(posterior_draws[[1]][,"beta1.gamma"],posterior_draws[[2]][,"beta1.gamma"],posterior_draws[[3]][,"beta1.gamma"]),
                               posteriors = rep(c("SAGA Wetness Index"), each = 18000),
                               submodel = rep(c("gamma"), each = 18000))
  
ggplot(post_sim_gamma, aes(x = samples, fill = posteriors)) + 
    geom_density(alpha = 0.5) + theme_ridges()
  
# Closed polygons
p0<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors)) + 
    geom_density_ridges2(fill = "#00AFBB")+xlim(-1.25,2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/gammaHedleyClassHybridInclLatitude.png",
         plot = p0,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
ggr_gamma<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors))+
    geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
    #scale_fill_gradientn(colours = terrain.colors(10))
    scale_fill_gradientn(colours =  c("red", "white", "blue"))+xlim(-5,5)+
    xlab("Posterior distrib.")+xlim(-1.25,2)+
    ylab("Colonization predictors")+ facet_grid(. ~ submodel)+
    theme(legend.position = "none")+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/gammaHedleyHybridInclLatitude_Coefsridgeplot.png",
         plot = ggr_gamma,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
# Combine the results in a single data frame: phi
post_sim_phi <- data.frame(samples = c(posterior_draws[[1]][,"beta1.phi"],posterior_draws[[2]][,"beta1.phi"],posterior_draws[[3]][,"beta1.phi"],
                                       posterior_draws[[1]][,"beta2.phi"],posterior_draws[[2]][,"beta2.phi"],posterior_draws[[3]][,"beta2.phi"],
                                       posterior_draws[[1]][,"beta3.phi"],posterior_draws[[2]][,"beta3.phi"],posterior_draws[[3]][,"beta3.phi"],
                                       posterior_draws[[1]][,"beta4.phi"],posterior_draws[[2]][,"beta4.phi"],posterior_draws[[3]][,"beta4.phi"],
                                       posterior_draws[[1]][,"beta5.phi"],posterior_draws[[2]][,"beta5.phi"],posterior_draws[[3]][,"beta5.phi"],
                                       posterior_draws[[1]][,"beta6.phi"],posterior_draws[[2]][,"beta6.phi"],posterior_draws[[3]][,"beta6.phi"],
                                       posterior_draws[[1]][,"beta7.phi"],posterior_draws[[2]][,"beta7.phi"],posterior_draws[[3]][,"beta7.phi"]),
                             posteriors = rep(c("Class 2 (Hedley 6)","Class 3 (Hedley 7)","Class 4 (Hedley 8)","Delta VH","Delta VH^2","static DPOL","Latitude^0.5"), each = 18000),
                             submodel = rep(c("phi","phi","phi","phi","phi","phi","phi"), each = 18000))
  
ggplot(post_sim_phi, aes(x = samples, fill = posteriors)) + 
    geom_density(alpha = 0.5) + theme_ridges()
  
  # Closed polygons
p0<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors)) + 
    geom_density_ridges2(fill = "#00AFBB")+xlim(-3,3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/phiHedleyClassHybridInclLatitude.png",
         plot = p0,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
ggr_phi<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors))+
    geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
    #scale_fill_gradientn(colours = terrain.colors(10))
    scale_fill_gradientn(colours =  c("red", "white", "blue"))+
    xlab("Posterior distrib.")+xlim(-3,3)+
    ylab("Persistence predictors")+ facet_grid(. ~ submodel)+
    theme(legend.position = "none")+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/phiHedleyHybridInclLatitude_Coefsridgeplot.png",
         plot = ggr_phi,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
  
# Combine the results in a single data frame: psi
post_sim_psi <- data.frame(samples = c(posterior_draws[[1]][,"beta1.psi"],posterior_draws[[2]][,"beta1.psi"],posterior_draws[[3]][,"beta1.psi"],
                                       posterior_draws[[1]][,"beta2.psi"],posterior_draws[[2]][,"beta2.psi"],posterior_draws[[3]][,"beta2.psi"],
                                       posterior_draws[[1]][,"beta3.psi"],posterior_draws[[2]][,"beta3.psi"],posterior_draws[[3]][,"beta3.psi"],
                                       posterior_draws[[1]][,"beta4.psi"],posterior_draws[[2]][,"beta4.psi"],posterior_draws[[3]][,"beta4.psi"],
                                       posterior_draws[[1]][,"beta5.psi"],posterior_draws[[2]][,"beta5.psi"],posterior_draws[[3]][,"beta5.psi"],
                                       posterior_draws[[1]][,"beta6.psi"],posterior_draws[[2]][,"beta6.psi"],posterior_draws[[3]][,"beta6.psi"]),
                             posteriors = rep(c("NDMI","NDMI^2","DeltaVH","static Grass150m","static Grass150m^2","Latitude^0.5"), each = 18000),
                             submodel = rep(c("psi","psi","psi","psi","psi","psi"), each = 18000))
  
ggplot(post_sim_psi, aes(x = samples, fill = posteriors)) + 
    geom_density(alpha = 0.5) + theme_ridges()
  
  # Closed polygons
p0<-ggplot(post_sim_psi, aes(x = samples, y = posteriors, group = posteriors)) + 
    geom_density_ridges2(fill = "#00AFBB")+xlim(-3,3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psiHedleyClassHybridInclLatitude.png",
         plot = p0,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
ggr_psi<-ggplot(post_sim_psi, aes(x = samples, y = posteriors, group = posteriors))+
    geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
    #scale_fill_gradientn(colours = terrain.colors(10))
    scale_fill_gradientn(colours =  c("red", "white", "blue"))+
    xlab("Posterior distrib.")+xlim(-3,3)+
    ylab("Initial occu. predictors")+ facet_grid(. ~ submodel)+
    theme(legend.position = "none")+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psiHedleyHybridInclLatitude_Coefsridgeplot.png",
         plot = ggr_psi,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
ggr_gamma2<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors))+
    geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
    #scale_fill_gradientn(colours = terrain.colors(10))
    scale_fill_gradientn(colours =  c("red", "white", "blue"))+xlim(-5,5)+
    xlab("Posterior distrib.")+xlim(-1.25,2)+
    ylab("Colonization predictors")+ facet_grid(. ~ submodel)+
    theme(legend.position = "none")+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )+
    theme(axis.text.y = element_blank())
ggr_phi2<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors))+
    geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
    #scale_fill_gradientn(colours = terrain.colors(10))
    scale_fill_gradientn(colours =  c("red", "white", "blue"))+
    xlab("Posterior distrib.")+xlim(-3,3)+
    ylab("Persistence predictors")+ facet_grid(. ~ submodel)+
    theme(legend.position = "none")+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )+
    theme(axis.text.y = element_blank())
  
g3<-grid.arrange(ggr_psi,
                   ggr_gamma,
                   ggr_phi,
                   nrow=1,ncol=3)
  
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/allHedleyHybridInclLatitude_Coefsridgeplot.png",
         plot = g3,
         width = 16,
         height = 8,
         units = c("in"),
         dpi = 300)
  
#Get WAIC
which(substr(colnames(posterior_draws[[1]]),1,23)=="log_likelihood_contrib[") #4375:6649 #There are 2275 (=65 sites*7 years*up to 5 visits per year)
colnames(posterior_draws[[1]])[5445]#"log_likelihood_contrib[1,1,1]"
posterior_draws[[1]][5445]
traceplot(posterior_draws[1:3][,5445], main="log_likelihood_contrib[1,1,1] Trace Plot", xlab="Iteration", ylab="Value", lwd=1)
densplot(posterior_draws[1:3][,"log_likelihood_contrib[1,1,1]"], main="log_likelihood_contrib[1,1,1] Density Plot", xlab="Value", ylab="Density", lwd=1)
smod$quantiles["log_likelihood_contrib[1,1,1]","97.5%"]
  
library(loo)
#waic(array): An I by C by N array, where I is the number of MCMC iterations per chain,
#C is the number of chains, and N is the number of data points.
#So in my case, a 60000 (or it is 6000=niter/nthin) * 3 * 8484
myArray<-array(NA,
                 dim=c(6000,3,2874))
  
for (i in 1:6000){
  for (j in 1:3){
    for (k in 1:2874){
      try(myArray[i,j,k]<-posterior_draws[[j]][i,(k+5444)])
    }
  }
}
waicMyArray<-waic(myArray)
#Computed from 18000 by 2874 log-likelihood matrix.
#Estimate   SE
#elpd_waic   -892.5 17.1
#p_waic        89.7  2.6
#waic        1784.9 34.1
looMyArray<-loo(myArray)
#Computed from 18000 by 2874 log-likelihood matrix.
#Estimate   SE
#elpd_loo   -892.6 17.0
#p_loo        89.8  2.5
#looic      1785.1 34.1
#------
#MCSE of elpd_loo is NA.
#MCSE and ESS estimates assume independent draws (r_eff=1).
#Pareto k diagnostic values:
#Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     1765  61.4%   14554   
#(0.7, 1]   (bad)         0   0.0%   <NA>    
#(1, Inf)   (very bad) 1109  38.6%   <NA>    
#See help('pareto-k-diagnostic') for details.
save(myArray, waicMyArray, looMyArray, file="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/WAIC_HedleyHybridInclLatitude2016_2021.Rdata")
  
#Get Estimated Number of Occupied Sites (out of 115)
which(substr(colnames(posterior_draws[[1]]),1,6)=="n.occ[")
# [1] 9014 9015 9016 9017 9018 9019 9020
colnames(posterior_draws[[1]])[9014]#n.occ[1]
  
n.occ1.PostMedian<-smod$quantiles["n.occ[1]","50%"]
n.occ1.Post.2.5<-smod$quantiles["n.occ[1]","2.5%"]
n.occ1.Post.97.5<-smod$quantiles["n.occ[1]","97.5%"]
  
n.occ2.PostMedian<-smod$quantiles["n.occ[2]","50%"]
n.occ2.Post.2.5<-smod$quantiles["n.occ[2]","2.5%"]
n.occ2.Post.97.5<-smod$quantiles["n.occ[2]","97.5%"]
  
n.occ3.PostMedian<-smod$quantiles["n.occ[3]","50%"]
n.occ3.Post.2.5<-smod$quantiles["n.occ[3]","2.5%"]
n.occ3.Post.97.5<-smod$quantiles["n.occ[3]","97.5%"]
  
n.occ4.PostMedian<-smod$quantiles["n.occ[4]","50%"]
n.occ4.Post.2.5<-smod$quantiles["n.occ[4]","2.5%"]
n.occ4.Post.97.5<-smod$quantiles["n.occ[4]","97.5%"]
  
n.occ5.PostMedian<-smod$quantiles["n.occ[5]","50%"]
n.occ5.Post.2.5<-smod$quantiles["n.occ[5]","2.5%"]
n.occ5.Post.97.5<-smod$quantiles["n.occ[5]","97.5%"]
  
#need to know #sites visited in each year:
nrow(data_jags$y[,,1])
year1<-data.frame(data_jags$y[,,1])
nrow(year1[!is.na(rowSums(year1)),])#82 
visited_year1<-year1[!is.na(rowSums(year1)),]#82
nrow(visited_year1[(rowSums(visited_year1)>0),])#47 
  
year2<-data.frame(data_jags$y[,,2])
nrow(year2[!is.na(rowSums(year2)),])#101
visited_year2<-year2[!is.na(rowSums(year2)),]#101
nrow(visited_year2[(rowSums(visited_year2)>0),])#58
  
year3<-data.frame(data_jags$y[,,3])
nrow(year3[!is.na(rowSums(year3)),])#62
visited_year3<-year3[!is.na(rowSums(year3)),]#62
nrow(visited_year3[(rowSums(visited_year3)>0),])#12
  
year4<-data.frame(data_jags$y[,,4])
nrow(year4[!is.na(rowSums(year4)),])#84
visited_year4<-year4[!is.na(rowSums(year4)),]#84
nrow(visited_year4[(rowSums(visited_year4)>0),])#30
  
year5<-data.frame(data_jags$y[,,5])
nrow(year5[!is.na(rowSums(year5)),])#79
visited_year5<-year5[!is.na(rowSums(year5)),]#79
nrow(visited_year5[(rowSums(visited_year5)>0),])#25
  
n.occDF<-data.frame(year=c(2016,2017,2018,2019,2021),
                      nsitesVisited=c(82,101,62,84,79),
                      nsitesDetected=c(47,58,12,30,25),
                      n.occMedian=c(n.occ1.PostMedian,
                                    n.occ2.PostMedian,
                                    n.occ3.PostMedian,
                                    n.occ4.PostMedian,
                                    n.occ5.PostMedian),
                      n.occ.2.5=c(n.occ1.Post.2.5,
                                  n.occ2.Post.2.5,
                                  n.occ3.Post.2.5,
                                  n.occ4.Post.2.5,
                                  n.occ5.Post.2.5),
                      n.occ.97.5=c(n.occ1.Post.97.5,
                                   n.occ2.Post.97.5,
                                   n.occ3.Post.97.5,
                                   n.occ4.Post.97.5,
                                   n.occ5.Post.97.5))
  
cor(n.occDF$nsitesVisited, n.occDF$n.occMedian)
#[1] 0.7058976
n.occDF$propDetected<-n.occDF$nsitesDetected/n.occDF$nsitesVisited
  
p1<-ggplot(n.occDF, aes(x=year, y=n.occMedian))+
    geom_point(aes(x=year, y=nsitesVisited), col="blue")+
    geom_line(aes(x=year, y=nsitesVisited), col="blue")+
    geom_point(aes(x=year, y=n.occMedian), col="red")+
    geom_line(aes(x=year, y=n.occMedian), col="red")+
    geom_line(aes(x=year, y=n.occ.2.5), linetype="dashed")+
    geom_line(aes(x=year, y=n.occ.97.5), linetype="dashed")+
    xlab("Year")+ylab("Number of Occupied Sites")+
    xlim(2016,2021)+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )+
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))
  
  
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridInclLatitude_2016_2021.png",
         plot = p1,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
p2<-ggplot(n.occDF, aes(x=year, y=propDetected))+
    geom_point(aes(x=year, y=propDetected), col="blue")+
    geom_line(aes(x=year, y=propDetected), col="blue")+
    xlab("Year")+ylab("Proportion of Visited Sites With Detections")+
    xlim(2016,2021)+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )+
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))
  
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridInclLatitude_2016_2021_SitesWithDetections.png",
         plot = p2,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
p3<-grid.arrange(p1,p2,nrow=1,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridInclLatitude_NumberOccupiedSitesPlusRaw_2016_2023.png",
         plot = p3,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
  
#Get Turnover
which(substr(colnames(posterior_draws[[1]]),1,9)=="turnover[")
#[1]8329 8330 8331 8332
colnames(posterior_draws[[1]])[8329]#turnover[1]
  
turnover1.PostMedian<-smod$quantiles["turnover[1]","50%"]
turnover1.Post.2.5<-smod$quantiles["turnover[1]","2.5%"]
turnover1.Post.97.5<-smod$quantiles["turnover[1]","97.5%"]
  
turnover2.PostMedian<-smod$quantiles["turnover[2]","50%"]
turnover2.Post.2.5<-smod$quantiles["turnover[2]","2.5%"]
turnover2.Post.97.5<-smod$quantiles["turnover[2]","97.5%"]
  
turnover3.PostMedian<-smod$quantiles["turnover[3]","50%"]
turnover3.Post.2.5<-smod$quantiles["turnover[3]","2.5%"]
turnover3.Post.97.5<-smod$quantiles["turnover[3]","97.5%"]
  
turnover4.PostMedian<-smod$quantiles["turnover[4]","50%"]
turnover4.Post.2.5<-smod$quantiles["turnover[4]","2.5%"]
turnover4.Post.97.5<-smod$quantiles["turnover[4]","97.5%"]
  

turnoverDF<-data.frame(year=c(2017,2018,2019,2021),
                         turnoverMedian=c(turnover1.PostMedian,
                                          turnover2.PostMedian,
                                          turnover3.PostMedian,
                                          turnover4.PostMedian),
                         turnover.2.5=c(turnover1.Post.2.5,
                                        turnover2.Post.2.5,
                                        turnover3.Post.2.5,
                                        turnover4.Post.2.5),
                         turnover.97.5=c(turnover1.Post.97.5,
                                         turnover2.Post.97.5,
                                         turnover3.Post.97.5,
                                         turnover4.Post.97.5))
  
p4<-ggplot(turnoverDF, aes(x=year, y=turnoverMedian))+
    geom_point(aes(x=year, y=turnoverMedian), col="red")+
    geom_line(aes(x=year, y=turnoverMedian), col="red")+
    geom_line(aes(x=year, y=turnover.2.5), linetype="dashed")+
    geom_line(aes(x=year, y=turnover.97.5), linetype="dashed")+
    xlab("Year")+ylab("Turnover from Previous Year")+
    theme(
      axis.text = element_text(size = 14), # Adjust all axis tick label sizes
      axis.title = element_text(size = 14) # Adjust all axis title sizes
    )
  
  
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridInclLatitude_2016_2021_PredictedTurnover.png",
         plot = p4,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  

#Get individual predictor plots



  
#Get AUC
  
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
colnames(posterior_draws[[1]])[4410]#"z[1,1]"
colnames(posterior_draws[[1]])[4525]#"z[1,2]"
colnames(posterior_draws[[1]])[4984]#"z[115,5]"
  
which(substr(colnames(posterior_draws[[1]]),1,4)=="muZ[")
colnames(posterior_draws[[1]])[4985]#"muZ[1,2]"
colnames(posterior_draws[[1]])[5444]#"muZ[115,5]"
str(posterior_draws[[1]][,"muZ[2,2]"])#no estimates for season 1
posterior_draws[[1]][,"muZ[1,2]"]#no estimates for season 1
  
niter<-ni
nthin<-nt
nchain<-nc
nsite<-115
nyear<-5
  
library(ggmcmc)
allmuZs<-posterior_draws[1:3][,4985:5444]#
ggd <- ggs(allmuZs)
allZs<-posterior_draws[1:3][,4525:4984]
ggd2 <- ggs(allZs)#starting from year 2
  
nstore <- niter / nthin
psi.pred <- array(dim=c(nsite, nstore*nchain))
Z.est <- array(dim=c(nsite, nstore*nchain))

#2017
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 2, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7861593
median(AUClist)#[1] 0.7886154
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6966912 0.7886154 0.8623279

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~97% above 0.7 

#2018
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 3, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.8366093
median(AUClist)#[1] 0.8386876
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7696102 0.8386876 0.8917023

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~100% above 0.7 

#2019
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 4, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.8147827
median(AUClist)#[1] 0.816414
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7420091 0.8164140 0.8780527

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~100% above 0.7 

#2021
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", nyear, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", nyear, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  
  
  
require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i
  
AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.801509
median(AUClist)#[1] 0.8031746
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7250337 0.8031746 0.8665194
  
# Specify the value to compare against
value_threshold <- 0.7
  
# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~99% above 0.7 
  
#Model is adequate
  
#Predicted yearly occu. for individual clusters
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
#575 entries = 115 sites * 5 years
library(ggmcmc)
allZs<-posterior_draws[1:3][,4410:4984]#
ggd <- ggs(allZs)
  
nsite<-115
niter <-60000#instead of 40000
nthin <-10
nchain <-3

nstore <- niter / nthin

#recreate obsFiltFinal dataset but with location and clusterID
#instead of YERA detections
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
#from the whole YERA study area. 2016 is first year
#for which we have Sentinel data.
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


m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits
str(m2)
table(m2$clusterID)#136 sites
# 2_1 2_10  2_8  3_1  3_6  3_8  4_5  4_8  7_5  7_8  9_1
#   3    5   15   13    7   11   45    2   33   21    5

table(m2$clusterIDfinal)
#2_1 2_10 2_8b  3_1  3_6  3_8  4_5  4_8  7_5 
#  3    5   15   13    7   11   45    2   33 
#7_8 7_8b 7_8c 7_8d  9_1 
# 13    1    3    4    5  

#get covariates affecting initial occu.
initFilt<-m2%>%
  #select(locationC, clusterID)%>%
  arrange(locationC)

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
#146 sites with 3+ years of data
#overlapping Hedley layer

#Value:  2  3  4  5  6  7  8 
#nsite:  3  7  3 22 35 36 40

initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11

obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID)
obsFiltFinal$geometry<-NULL

#Probabilities of occu. at year 1 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 1, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2016<-rowMeans(z.pred)
OccEst2016<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2016$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2016bycluster<-OccEst2016%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2016bycluster_quantiles <- OccEst2016bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2016)

#Probabilities of occu. at year 2 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2017<-rowMeans(z.pred)
OccEst2017<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2017$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2017bycluster<-OccEst2017%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2017bycluster_quantiles <- OccEst2017bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2017)


#Probabilities of occu. at year 3 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2018<-rowMeans(z.pred)
OccEst2018<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2018$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2018bycluster<-OccEst2018%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2018bycluster_quantiles <- OccEst2018bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2018)


#Probabilities of occu. at year 4 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2019<-rowMeans(z.pred)
OccEst2019<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2019$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2019bycluster<-OccEst2019%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2019bycluster_quantiles <- OccEst2019bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2019)

#Probabilities of occu. at year 5 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 5, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2021<-rowMeans(z.pred)
OccEst2021<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2021$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2021bycluster<-OccEst2021%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2021bycluster_quantiles <- OccEst2021bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2021)

#meanOccEstAllYears<-cbind(meanOccEst2016,
#                          meanOccEst2017,
#                          meanOccEst2018,
#                          meanOccEst2019,
#                          meanOccEst2021)

finalOccDF<-rbind(OccEst2016bycluster_quantiles,
                  OccEst2017bycluster_quantiles,
                  OccEst2018bycluster_quantiles,
                  OccEst2019bycluster_quantiles,
                  OccEst2021bycluster_quantiles)


#finalOccDF<-cbind(obsFiltFinal, meanOccEstAllYears)
#finalOccDF$clusterID
#write.csv(finalOccDF,file="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021_PredictedOcc.csv")


#Cluster 1 (formerly 3_6)
finalOccDF_cluster1<-finalOccDF%>%
  filter(clusterID=="3_6")

# finalOccDF_cluster1long<-finalOccDF_cluster1%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster1)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_6 <- ggplot(data=finalOccDF_cluster1long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 1: ",nrow(finalOccDF_cluster1)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_6 <- ggplot(data=finalOccDF_cluster1, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 1: ",nrow(obsFiltFinal%>%filter(clusterID=="3_6"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 2
finalOccDF_cluster2<-finalOccDF%>%
  filter(clusterID=="4_5")

# finalOccDF_cluster2long<-finalOccDF_cluster2%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster2)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster4_5 <- ggplot(data=finalOccDF_cluster2long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 2: ",nrow(finalOccDF_cluster2)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster4_5 <- ggplot(data=finalOccDF_cluster2, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 2: ",nrow(obsFiltFinal%>%filter(clusterID=="4_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 3
finalOccDF_cluster3<-finalOccDF%>%
  filter(clusterID=="7_5")

# finalOccDF_cluster3long<-finalOccDF_cluster3%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster3)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

#Occucluster7_5 <- ggplot(data=finalOccDF_cluster3long, aes(x=Year, y=medYear)) +
#  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#  geom_line(aes(x=Year, y=medYear), colour="blue")+
#  my.theme+ggtitle(paste0("Cluster 3: ",nrow(finalOccDF_cluster3)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#  ylim(-0.5,1.25)+xlim(2016,2021)+
#  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster7_5 <- ggplot(data=finalOccDF_cluster3, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 3: ",nrow(obsFiltFinal%>%filter(clusterID=="7_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 4 (just 3 sites at GG Provincial Park)
#No sites with 3 or more years of obs from 2016-2021 in cluster 4
# finalOccDF_cluster4<-finalOccDF%>%
#     filter(clusterID=="2_1")
#   
# finalOccDF_cluster4long<-finalOccDF_cluster4%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster4)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster2_1 <- ggplot(data=finalOccDF_cluster4long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 4: ",nrow(finalOccDF_cluster4)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 5
finalOccDF_cluster5<-finalOccDF%>%
  filter(clusterID=="3_1")

# finalOccDF_cluster5long<-finalOccDF_cluster5%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster5)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_1 <- ggplot(data=finalOccDF_cluster5long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 5: ",nrow(finalOccDF_cluster5)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_1 <- ggplot(data=finalOccDF_cluster5, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 5: ",nrow(obsFiltFinal%>%filter(clusterID=="3_1"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 6
finalOccDF_cluster6<-finalOccDF%>%
  filter(clusterID=="7_8")

# finalOccDF_cluster6long<-finalOccDF_cluster6%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster6)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster7_8 <- ggplot(data=finalOccDF_cluster6long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 6: ",nrow(finalOccDF_cluster6)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


Occucluster7_8 <- ggplot(data=finalOccDF_cluster6, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 6: ",nrow(obsFiltFinal%>%filter(clusterID=="7_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 7
finalOccDF_cluster7<-finalOccDF%>%
  filter(clusterID=="3_8")

# finalOccDF_cluster7long<-finalOccDF_cluster7%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster7)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_8 <- ggplot(data=finalOccDF_cluster7long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 7: ",nrow(finalOccDF_cluster7)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_8 <- ggplot(data=finalOccDF_cluster7, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 7: ",nrow(obsFiltFinal%>%filter(clusterID=="3_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 8
finalOccDF_cluster8<-finalOccDF%>%
  filter(clusterID=="2_8")

# finalOccDF_cluster8long<-finalOccDF_cluster8%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster8)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster2_8 <- ggplot(data=finalOccDF_cluster8long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 8: ",nrow(finalOccDF_cluster8)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster2_8 <- ggplot(data=finalOccDF_cluster8, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 8: ",nrow(obsFiltFinal%>%filter(clusterID=="2_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 9 No sites with at least 3 years of data from 2016-2021 in this cluster
# finalOccDF_cluster9<-finalOccDF%>%
#     filter(clusterID=="4_8")
#   
# finalOccDF_cluster9long<-finalOccDF_cluster9%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster9)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster4_8 <- ggplot(data=finalOccDF_cluster9long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 9: ",nrow(finalOccDF_cluster9)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

p9<-grid.arrange(Occucluster3_6,
                 Occucluster4_5,
                 Occucluster7_5,
                 #Occucluster2_1,
                 Occucluster3_1,
                 Occucluster7_8,
                 Occucluster3_8,
                 Occucluster2_8,
                 #Occucluster4_8,
                 nrow=3,ncol=3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridInclLatitude_MedianPlus95BCIoccupancy_2016_2023Dec8.png",
         plot = p9,
         width = 12,
         height = 8,
         units = c("in"),
         dpi = 300)
  
  
#EXCLUDING LATITUDE
#####################################################################################################################
#                                                                                                                   #
#                                                                                                                   #
#                       Create Model File and Initialize Parameters - Modelling Individual Footprint Types          #
#                                                                                                                   #
#                                                                                                                   #
#                                                                                                                   #
#                                                                                                                   #
#####################################################################################################################

#Run multi-season occu. model in JAGS
#Same code as above, just with one less predictor for psi and phi
# Specify model in JAGS language
model_code <-
  "model
{
  # Specify priors
  alpha.psi ~ dunif(-3, 3) 
  beta1.psi ~ dunif(-3, 3)
  beta2.psi ~ dunif(-3, 3)
  beta3.psi ~ dunif(-3, 3) 
  beta4.psi ~ dunif(-3, 3)
  beta5.psi ~ dunif(-3, 3)

  alpha.gamma ~ dunif(-3, 3)
  beta1.gamma ~ dunif(-3, 3) 

  alpha.phi ~ dunif(-3, 3)  
  beta1.phi ~ dunif(-3, 3)
  beta2.phi ~ dunif(-3, 3)
  beta3.phi ~ dunif(-3, 3)
  beta4.phi ~ dunif(-3, 3)
  beta5.phi ~ dunif(-3, 3)
  beta6.phi ~ dunif(-3, 3)

  alpha.p ~ dunif(-3, 3) #no detection covariates for now
  
  #random effects on psi, gamma and phi 
  #115 sites, 7 cluster IDs
  for (m in 1:ncluster) {
    u_psi[m] ~ dnorm(0, tau_psi)
    u_phi[m] ~ dnorm(0, tau_phi) #Wait, why dnorm and not dunif?
    u_gamma[m] ~ dnorm(0, tau_gamma)
  }
  tau_psi <- pow(sigma_psi, -2)
  sigma_psi ~ dunif(0, 5)

  tau_phi <- pow(sigma_phi, -2)
  sigma_phi ~ dunif(0, 5)

  tau_gamma <- pow(sigma_gamma, -2)
  sigma_gamma ~ dunif(0, 5)
  

  #assign values generated for u_phi and u_gamma to a vector of size 93 based
  #on clusterID?


  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){
    z0[i] ~dbern(psi1[i]) #True occu. at site i in season 1
    logit(psi1[i])<-alpha.psi +beta1.psi*ndmi[i] +beta2.psi*ndmi2[i] +beta3.psi*deltaVH[i] +beta4.psi*sgrass[i] +beta5.psi*sgrass2[i] + u_psi[cluster[i]]

    z[i,1]<-z0[i]
    for (k in 2:nyear){
      #survival probability from season k-1 to season k
      phi[i,k-1] <-psurv[i,k-1]
      logit(psurv[i,k-1])<-alpha.phi + u_phi[cluster[i]] +beta1.phi*class2SY[i,k-1] +beta2.phi*class3SY[i,k-1] +beta3.phi*class4SY[i,k-1] +beta4.phi*YdeltaVH[i,k-1] +beta5.phi*YdeltaVH2[i,k-1] +beta6.phi*YsDPOLsummer[i,k-1] 

      #colonization probability from season k-1 to season k
      gamma[i,k-1] <-pcol[i,k-1]
      logit(pcol[i,k-1])<-alpha.gamma  +u_gamma[cluster[i]] +beta1.gamma*Yswi[i,k-1] 

      muZ[i,k]<-z[i,k-1]*phi[i,k-1] +(1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~dbern(muZ[i,k])
    } #k
  } #i
  # Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (k in 1:nyear){
        p[i,j,k]<-det[i,j,k]
        logit(det[i,j,k])<-alpha.p    #(no detection covariates for now)    
        muy[i,j,k]<-z[i,k]*p[i,j,k]
        y[i,j,k]~dbern(muy[i,j,k])
        log_likelihood_contrib[i,j,k]<-y[i,j,k]*log(muy[i,j,k]) + (1-y[i,j,k])*log(1-muy[i,j,k])
      } #k
    } #j
  } #i
  # Derived parameters: Sample and population occu., growth rate and turnover
  n.occ[1]<-sum(z[1:nsite,1])
  for (i in 1:nsite){
    psi[i,1] <-psi1[i]
    for (k in 2:nyear){
      psi[i,k] <-psi[i,k-1]*phi[i,k-1] +(1-psi[i,k-1])*gamma[i,k-1]
    }
  }
  for (k in 2:nyear){
      n.occ[k] <-sum(z[1:nsite,k])
      growthr[k-1]<-mean(psi[1:nsite,k-1])/mean(psi[1:nsite,k-1])
      #no longer growthr[k]<-psi[k]/psi[k-1] or psi[i,k]/psi[i,k-1]
      turnover[k-1] <-(1 - mean(psi[1:nsite,k-1])) *mean(gamma[1:nsite,k-1])/mean(psi[1:nsite,k])
      #turnover[k-1] <-(1 - psi[i,k-1]) *gamma[k-1]/psi[i,k]
  }
  
}"
  
  
# Bundle data
nsite<-115
nyear<-5
nrep<-5
data_jags <-list(y=y,
                   #thr=thr,
                   nsite=nsite,
                   nrep=nrep,
                   nyear=nyear,
                   #psi covariates
                   ndmi=ndmi,
                   ndmi2=ndmi2,
                   deltaVH=deltaVH,
                   sgrass=sgrass,
                   sgrass2=sgrass2,
                   #dropped latsqrt
                   #site-year covariates (gamma, phi)
                   class2SY=class2SY,
                   class3SY=class3SY,
                   class4SY=class4SY,
                   Yswi=Yswi,
                   YdeltaVH=YdeltaVH,
                   YdeltaVH2=YdeltaVH2,
                   YsDPOLsummer=YsDPOLsummer,
                   #dropped Ylatsqrt
                   #cluster information
                   ncluster=7,
                   cluster=clusterID)
#no detection covariates (p)
  
# Initial values
zst <-apply(y,c(1,3),max,na.rm=TRUE)
#Observed occurrence as inits for z
zst[zst==-Inf]<-0
  
#Run model via parallel processing:
initsP = parallel.seeds("base::BaseRNG", 3) # a list of lists
initsP[[1]]$z=zst
initsP[[1]]$alpha.psi=runif(1, -3, 3)
initsP[[1]]$beta1.psi=runif(1, -3, 3)
initsP[[1]]$beta2.psi=runif(1, -3, 3)
initsP[[1]]$beta3.psi=runif(1, -3, 3)
initsP[[1]]$beta4.psi=runif(1, -3, 3)
initsP[[1]]$beta5.psi=runif(1, -3, 3)
initsP[[1]]$u_psi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.phi=runif(1, -3, 3)
initsP[[1]]$beta1.phi =runif(1, -3, 3)
initsP[[1]]$beta2.phi =runif(1, -3, 3)
initsP[[1]]$beta3.phi =runif(1, -3, 3)
initsP[[1]]$beta4.phi =runif(1, -3, 3)
initsP[[1]]$beta5.phi =runif(1, -3, 3)
initsP[[1]]$beta6.phi =runif(1, -3, 3)
initsP[[1]]$u_phi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.gamma=runif(1, -3, 3)
initsP[[1]]$beta1.gamma =runif(1, -3, 3)
initsP[[1]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.p=runif(1,-3,3)
  
initsP[[2]]$z=zst
initsP[[2]]$alpha.psi=runif(1, -3, 3)
initsP[[2]]$beta1.psi=runif(1, -3, 3)
initsP[[2]]$beta2.psi=runif(1, -3, 3)
initsP[[2]]$beta3.psi=runif(1, -3, 3)
initsP[[2]]$beta4.psi=runif(1, -3, 3)
initsP[[2]]$beta5.psi=runif(1, -3, 3)
initsP[[2]]$u_psi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.phi=runif(1, -3, 3)
initsP[[2]]$beta1.phi =runif(1, -3, 3)
initsP[[2]]$beta2.phi =runif(1, -3, 3)
initsP[[2]]$beta3.phi =runif(1, -3, 3)
initsP[[2]]$beta4.phi =runif(1, -3, 3)
initsP[[2]]$beta5.phi =runif(1, -3, 3)
initsP[[2]]$beta6.phi =runif(1, -3, 3)
initsP[[2]]$u_phi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.gamma=runif(1, -3, 3)
initsP[[2]]$beta1.gamma =runif(1, -3, 3)
initsP[[2]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.p=runif(1,-3,3)
  
initsP[[3]]$z=zst
initsP[[3]]$alpha.psi=runif(1, -3, 3)
initsP[[3]]$beta1.psi=runif(1, -3, 3)
initsP[[3]]$beta2.psi=runif(1, -3, 3)
initsP[[3]]$beta3.psi=runif(1, -3, 3)
initsP[[3]]$beta4.psi=runif(1, -3, 3)
initsP[[3]]$beta5.psi=runif(1, -3, 3)
initsP[[3]]$u_psi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.phi=runif(1, -3, 3)
initsP[[3]]$beta1.phi =runif(1, -3, 3)
initsP[[3]]$beta2.phi =runif(1, -3, 3)
initsP[[3]]$beta3.phi =runif(1, -3, 3)
initsP[[3]]$beta4.phi =runif(1, -3, 3)
initsP[[3]]$beta5.phi =runif(1, -3, 3)
initsP[[3]]$beta6.phi =runif(1, -3, 3)
initsP[[3]]$u_phi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.gamma=runif(1, -3, 3)
initsP[[3]]$beta1.gamma =runif(1, -3, 3)
initsP[[3]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.p=runif(1,-3,3)
  
# Parameters monitored
params <-c("psi",
             "alpha.psi",
             "beta1.psi",
             "beta2.psi",
             "beta3.psi",
             "beta4.psi",
             "beta5.psi",
             "u_psi",
             "phi",
             "alpha.phi",
             "beta1.phi",
             "beta2.phi",
             "beta3.phi",
             "beta4.phi",
             "beta5.phi",
             "beta6.phi",
             "u_phi",
             "gamma",
             "alpha.gamma",
             "beta1.gamma",
             "u_gamma",
             "alpha.p",
             "p","z","muZ",
             "log_likelihood_contrib",
             "n.occ",
             "growthr",
             "turnover")#,"sens","fpr","auc"
  
  
# Running the model via parallel processing using runjags
ni <-60000#instead of 40000
nt <-10
nb <-50000#instead of 30000
nc <-3
#so there should be 1000 posterior draws
  
niter<-ni
nchains<-nc
thinning<-nt
  
results <- runjags::run.jags(
    model = model_code,
    monitor = params,#keep_track_of,
    data = data_jags,#datalist
    n.chains = nchains,
    inits=initsP,
    adapt = niter/10,
    burnin = nb, 
    sample = niter/thinning, 
    thin = thinning,
    method = "parallel"
  )
cat("Model run completed\n")
print(Sys.time())
  
# Extracting results
posterior_draws <- as.mcmc.list(results)
cat("Posterior draws completed\n")
smod <- summary(posterior_draws)
  
save(detectionsLong, 
       data_jags, 
       results, 
       posterior_draws, 
       smod, 
       nc, 
       ni, 
       nb, 
       nt, 
       params, 
       file = "2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_HedleyPlusHybridExclLat.Rdata")
  
load("2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_HedleyPlusHybridExclLat.Rdata")
#Process the MCMC output:
int.mcmc.mat <- as.matrix(posterior_draws)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)
names(int.mcmc.dat)
int.mcmc.dat[,"phi[64,1]"]#18000 estimates for this parameter


# Diagnostics:                                        Point est. Upper C.I.
gelman.diag(posterior_draws[1:3][,"alpha.psi"])#[1,]          1       1       
gelman.diag(posterior_draws[1:3][,"beta1.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta2.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta3.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta4.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta5.psi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.phi"])#[1,]          1          1       
gelman.diag(posterior_draws[1:3][,"beta1.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta2.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta3.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta4.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta5.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta6.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.gamma"])#[1,]          1          1       
gelman.diag(posterior_draws[1:3][,"beta1.gamma"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.p"])#[1,]          1          1           

range(
  effectiveSize(posterior_draws[1:3][,"alpha.psi"]),        
  effectiveSize(posterior_draws[1:3][,"beta1.psi"]),           
  effectiveSize(posterior_draws[1:3][,"beta2.psi"]),           
  effectiveSize(posterior_draws[1:3][,"beta3.psi"]),           
  effectiveSize(posterior_draws[1:3][,"beta4.psi"]),           
  effectiveSize(posterior_draws[1:3][,"beta5.psi"]),           
  effectiveSize(posterior_draws[1:3][,"alpha.phi"]),        
  effectiveSize(posterior_draws[1:3][,"beta1.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta2.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta3.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta4.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta5.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta6.phi"]),           
  effectiveSize(posterior_draws[1:3][,"alpha.gamma"]),        
  effectiveSize(posterior_draws[1:3][,"beta1.gamma"])           
)
#effective sample size varies from 3921.372 16423.842

heidel.diag(posterior_draws[1:3][,"alpha.psi"])#passed both        
heidel.diag(posterior_draws[1:3][,"beta1.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta2.psi"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"beta3.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta4.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta5.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.phi"])#passed first 2 chains, failed third        
heidel.diag(posterior_draws[1:3][,"beta1.phi"])#passed both          
heidel.diag(posterior_draws[1:3][,"beta2.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta3.phi"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"beta4.phi"])#passed both          
heidel.diag(posterior_draws[1:3][,"beta5.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta6.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.gamma"])#passed both        
heidel.diag(posterior_draws[1:3][,"beta1.gamma"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.p"])#passed both           
#All passed the stationary test (p>0.1), failed the halfwidth test

gd <- gelman.diag(posterior_draws[1:3][,"alpha.psi"])
es <- effectiveSize(posterior_draws[1:3][,"alpha.psi"])#   
hd <- heidel.diag(posterior_draws[1:3][,"alpha.psi"])#

#Trace Plots and Density Plots
plot(posterior_draws[1:3][,"beta1.psi"],density=FALSE)#- effect on psi
densplot(posterior_draws[[1]][,"beta1.psi"], main="beta1.psi (NDMI) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta2.psi"],density=FALSE)#nonsig effect on psi
densplot(posterior_draws[[1]][,"beta2.psi"], main="beta2.psi (NDMI2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta3.psi"],density=FALSE)#- effect on psi
densplot(posterior_draws[[1]][,"beta3.psi"], main="beta3.psi (DeltaVH) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta4.psi"],density=FALSE)#+ effect on psi
densplot(posterior_draws[[1]][,"beta4.psi"], main="beta4.psi (sgrass) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta5.psi"],density=FALSE)#+ effect on psi
densplot(posterior_draws[[1]][,"beta5.psi"], main="beta5.psi (sgrass2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta1.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta1.phi"], main="beta1.phi (habitat category 2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta2.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta2.phi"], main="beta2.phi (habitat category 3) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta3.phi"],density=FALSE)#nonsig effect on phi
densplot(posterior_draws[[1]][,"beta3.phi"], main="beta3.phi (habitat category 4) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta4.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta4.phi"], main="beta4.phi (DeltaVH) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta5.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta5.phi"], main="beta5.phi (DeltaVH2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta6.phi"],density=FALSE)#- effect on phi
densplot(posterior_draws[[1]][,"beta6.phi"], main="beta6.phi (sDPOLsummer) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta1.gamma"],density=FALSE)#+ effect on gamma
densplot(posterior_draws[[1]][,"beta1.gamma"], main="beta1.gamma (SWI) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)


#create ggridges plots
# Combine the results in a single data frame: gamma
post_sim_gamma <- data.frame(samples = c(posterior_draws[[1]][,"beta1.gamma"],posterior_draws[[2]][,"beta1.gamma"],posterior_draws[[3]][,"beta1.gamma"]),
                             posteriors = rep(c("SAGA Wetness Index"), each = 18000),
                             submodel = rep(c("gamma"), each = 18000))

ggplot(post_sim_gamma, aes(x = samples, fill = posteriors)) + 
  geom_density(alpha = 0.5) + theme_ridges()

# Closed polygons
p0<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors)) + 
  geom_density_ridges2(fill = "#00AFBB")+xlim(-1.25,2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/gammaHedleyClassHybridExclLatitude.png",
       plot = p0,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

ggr_gamma<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors))+
  geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
  #scale_fill_gradientn(colours = terrain.colors(10))
  scale_fill_gradientn(colours =  c("red", "white", "blue"))+xlim(-5,5)+
  xlab("Posterior distrib.")+xlim(-1.25,2)+
  ylab("Colonization predictors")+ facet_grid(. ~ submodel)+
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/gammaHedleyHybridExclLatitude_Coefsridgeplot.png",
       plot = ggr_gamma,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

# Combine the results in a single data frame: phi
post_sim_phi <- data.frame(samples = c(posterior_draws[[1]][,"beta1.phi"],posterior_draws[[2]][,"beta1.phi"],posterior_draws[[3]][,"beta1.phi"],
                                       posterior_draws[[1]][,"beta2.phi"],posterior_draws[[2]][,"beta2.phi"],posterior_draws[[3]][,"beta2.phi"],
                                       posterior_draws[[1]][,"beta3.phi"],posterior_draws[[2]][,"beta3.phi"],posterior_draws[[3]][,"beta3.phi"],
                                       posterior_draws[[1]][,"beta4.phi"],posterior_draws[[2]][,"beta4.phi"],posterior_draws[[3]][,"beta4.phi"],
                                       posterior_draws[[1]][,"beta5.phi"],posterior_draws[[2]][,"beta5.phi"],posterior_draws[[3]][,"beta5.phi"],
                                       posterior_draws[[1]][,"beta6.phi"],posterior_draws[[2]][,"beta6.phi"],posterior_draws[[3]][,"beta6.phi"]),
                           posteriors = rep(c("Class 2 (Hedley 6)","Class 3 (Hedley 7)","Class 4 (Hedley 8)","Delta VH","Delta VH^2","static DPOL"), each = 18000),
                           submodel = rep(c("phi","phi","phi","phi","phi","phi"), each = 18000))

ggplot(post_sim_phi, aes(x = samples, fill = posteriors)) + 
  geom_density(alpha = 0.5) + theme_ridges()

# Closed polygons
p0<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors)) + 
  geom_density_ridges2(fill = "#00AFBB")+xlim(-3,3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/phiHedleyClassHybridExclLatitude.png",
       plot = p0,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)
ggr_phi<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors))+
  geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
  #scale_fill_gradientn(colours = terrain.colors(10))
  scale_fill_gradientn(colours =  c("red", "white", "blue"))+
  xlab("Posterior distrib.")+xlim(-3,3)+
  ylab("Persistence predictors")+ facet_grid(. ~ submodel)+
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/phiHedleyHybridExclLatitude_Coefsridgeplot.png",
       plot = ggr_phi,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


# Combine the results in a single data frame: psi
post_sim_psi <- data.frame(samples = c(posterior_draws[[1]][,"beta1.psi"],posterior_draws[[2]][,"beta1.psi"],posterior_draws[[3]][,"beta1.psi"],
                                       posterior_draws[[1]][,"beta2.psi"],posterior_draws[[2]][,"beta2.psi"],posterior_draws[[3]][,"beta2.psi"],
                                       posterior_draws[[1]][,"beta3.psi"],posterior_draws[[2]][,"beta3.psi"],posterior_draws[[3]][,"beta3.psi"],
                                       posterior_draws[[1]][,"beta4.psi"],posterior_draws[[2]][,"beta4.psi"],posterior_draws[[3]][,"beta4.psi"],
                                       posterior_draws[[1]][,"beta5.psi"],posterior_draws[[2]][,"beta5.psi"],posterior_draws[[3]][,"beta5.psi"]),
                           posteriors = rep(c("NDMI","NDMI^2","DeltaVH","static Grass150m","static Grass150m^2"), each = 18000),
                           submodel = rep(c("psi","psi","psi","psi","psi"), each = 18000))

ggplot(post_sim_psi, aes(x = samples, fill = posteriors)) + 
  geom_density(alpha = 0.5) + theme_ridges()

# Closed polygons
p0<-ggplot(post_sim_psi, aes(x = samples, y = posteriors, group = posteriors)) + 
  geom_density_ridges2(fill = "#00AFBB")+xlim(-3,3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psiHedleyClassHybridExclLatitude.png",
       plot = p0,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)
ggr_psi<-ggplot(post_sim_psi, aes(x = samples, y = posteriors, group = posteriors))+
  geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
  #scale_fill_gradientn(colours = terrain.colors(10))
  scale_fill_gradientn(colours =  c("red", "white", "blue"))+
  xlab("Posterior distrib.")+xlim(-3,3)+
  ylab("Initial occu. predictors")+ facet_grid(. ~ submodel)+
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psiHedleyHybridExclLatitude_Coefsridgeplot.png",
       plot = ggr_psi,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

ggr_gamma2<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors))+
  geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
  #scale_fill_gradientn(colours = terrain.colors(10))
  scale_fill_gradientn(colours =  c("red", "white", "blue"))+xlim(-5,5)+
  xlab("Posterior distrib.")+xlim(-1.25,2)+
  ylab("Colonization predictors")+ facet_grid(. ~ submodel)+
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  theme(axis.text.y = element_blank())
ggr_phi2<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors))+
  geom_density_ridges_gradient(aes(fill=after_stat(x)), scale = 3, rel_min_height = 0.01) +
  #scale_fill_gradientn(colours = terrain.colors(10))
  scale_fill_gradientn(colours =  c("red", "white", "blue"))+
  xlab("Posterior distrib.")+xlim(-3,3)+
  ylab("Persistence predictors")+ facet_grid(. ~ submodel)+
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  theme(axis.text.y = element_blank())

g3<-grid.arrange(ggr_psi,
                 ggr_gamma,
                 ggr_phi,
                 nrow=1,ncol=3)

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/allHedleyHybridExclLatitude_Coefsridgeplot.png",
       plot = g3,
       width = 16,
       height = 8,
       units = c("in"),
       dpi = 300)







#Get WAIC
which(substr(colnames(posterior_draws[[1]]),1,23)=="log_likelihood_contrib[") #4375:6649 #There are 2275 (=65 sites*7 years*up to 5 visits per year)
colnames(posterior_draws[[1]])[5443]#"log_likelihood_contrib[1,1,1]"
posterior_draws[[1]][5443]
traceplot(posterior_draws[1:3][,5443], main="log_likelihood_contrib[1,1,1] Trace Plot", xlab="Iteration", ylab="Value", lwd=1)
densplot(posterior_draws[1:3][,"log_likelihood_contrib[1,1,1]"], main="log_likelihood_contrib[1,1,1] Density Plot", xlab="Value", ylab="Density", lwd=1)
smod$quantiles["log_likelihood_contrib[1,1,1]","97.5%"]

library(loo)
#waic(array): An I by C by N array, where I is the number of MCMC iterations per chain,
#C is the number of chains, and N is the number of data points.
#So in my case, a 60000 (or it is 6000=niter/nthin) * 3 * 8484
myArray<-array(NA,
               dim=c(6000,3,2874))

for (i in 1:6000){
  for (j in 1:3){
    for (k in 1:2874){
      try(myArray[i,j,k]<-posterior_draws[[j]][i,(k+5442)])
    }
  }
}
waicMyArray<-waic(myArray)
#Computed from 18000 by 2874 log-likelihood matrix.
#Estimate   SE
#elpd_waic   -885.9 17.0
#p_waic        91.0  2.6
#waic        1771.7 34.1
looMyArray<-loo(myArray)
#Computed from 18000 by 2874 log-likelihood matrix.
#Estimate   SE
#elpd_loo   -886.0 17.0
#p_loo        91.2  2.6
#looic      1772.1 34.0
#MCSE of elpd_loo is NA.
#MCSE and ESS estimates assume independent draws (r_eff=1).
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     1758  61.2%   14538   
#(0.7, 1]   (bad)         2   0.1%   <NA>    
#(1, Inf)   (very bad) 1114  38.8%   <NA>    
#See help('pareto-k-diagnostic') for details.
save(myArray, waicMyArray, looMyArray, file="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/WAIC_HedleyHybridExclLatitude2016_2021.Rdata")

#Get Estimated Number of Occupied Sites (out of 115)
which(substr(colnames(posterior_draws[[1]]),1,6)=="n.occ[")
# [1] 8318 8319 8320 8321 8322
colnames(posterior_draws[[1]])[8318]#n.occ[1]

n.occ1.PostMedian<-smod$quantiles["n.occ[1]","50%"]
n.occ1.Post.2.5<-smod$quantiles["n.occ[1]","2.5%"]
n.occ1.Post.97.5<-smod$quantiles["n.occ[1]","97.5%"]

n.occ2.PostMedian<-smod$quantiles["n.occ[2]","50%"]
n.occ2.Post.2.5<-smod$quantiles["n.occ[2]","2.5%"]
n.occ2.Post.97.5<-smod$quantiles["n.occ[2]","97.5%"]

n.occ3.PostMedian<-smod$quantiles["n.occ[3]","50%"]
n.occ3.Post.2.5<-smod$quantiles["n.occ[3]","2.5%"]
n.occ3.Post.97.5<-smod$quantiles["n.occ[3]","97.5%"]

n.occ4.PostMedian<-smod$quantiles["n.occ[4]","50%"]
n.occ4.Post.2.5<-smod$quantiles["n.occ[4]","2.5%"]
n.occ4.Post.97.5<-smod$quantiles["n.occ[4]","97.5%"]

n.occ5.PostMedian<-smod$quantiles["n.occ[5]","50%"]
n.occ5.Post.2.5<-smod$quantiles["n.occ[5]","2.5%"]
n.occ5.Post.97.5<-smod$quantiles["n.occ[5]","97.5%"]

#need to know #sites visited in each year:
nrow(data_jags$y[,,1])
year1<-data.frame(data_jags$y[,,1])
nrow(year1[!is.na(rowSums(year1)),])#82 
visited_year1<-year1[!is.na(rowSums(year1)),]#82
nrow(visited_year1[(rowSums(visited_year1)>0),])#47 

year2<-data.frame(data_jags$y[,,2])
nrow(year2[!is.na(rowSums(year2)),])#101
visited_year2<-year2[!is.na(rowSums(year2)),]#101
nrow(visited_year2[(rowSums(visited_year2)>0),])#58

year3<-data.frame(data_jags$y[,,3])
nrow(year3[!is.na(rowSums(year3)),])#62
visited_year3<-year3[!is.na(rowSums(year3)),]#62
nrow(visited_year3[(rowSums(visited_year3)>0),])#12

year4<-data.frame(data_jags$y[,,4])
nrow(year4[!is.na(rowSums(year4)),])#84
visited_year4<-year4[!is.na(rowSums(year4)),]#84
nrow(visited_year4[(rowSums(visited_year4)>0),])#30

year5<-data.frame(data_jags$y[,,5])
nrow(year5[!is.na(rowSums(year5)),])#79
visited_year5<-year5[!is.na(rowSums(year5)),]#79
nrow(visited_year5[(rowSums(visited_year5)>0),])#25

n.occDF<-data.frame(year=c(2016,2017,2018,2019,2021),
                    nsitesVisited=c(82,101,62,84,79),
                    nsitesDetected=c(47,58,12,30,25),
                    n.occMedian=c(n.occ1.PostMedian,
                                  n.occ2.PostMedian,
                                  n.occ3.PostMedian,
                                  n.occ4.PostMedian,
                                  n.occ5.PostMedian),
                    n.occ.2.5=c(n.occ1.Post.2.5,
                                n.occ2.Post.2.5,
                                n.occ3.Post.2.5,
                                n.occ4.Post.2.5,
                                n.occ5.Post.2.5),
                    n.occ.97.5=c(n.occ1.Post.97.5,
                                 n.occ2.Post.97.5,
                                 n.occ3.Post.97.5,
                                 n.occ4.Post.97.5,
                                 n.occ5.Post.97.5))

cor(n.occDF$nsitesVisited, n.occDF$n.occMedian)
#[1] 0.7058976
n.occDF$propDetected<-n.occDF$nsitesDetected/n.occDF$nsitesVisited

p1<-ggplot(n.occDF, aes(x=year, y=n.occMedian))+
  geom_point(aes(x=year, y=nsitesVisited), col="blue")+
  geom_line(aes(x=year, y=nsitesVisited), col="blue")+
  geom_point(aes(x=year, y=n.occMedian), col="red")+
  geom_line(aes(x=year, y=n.occMedian), col="red")+
  geom_line(aes(x=year, y=n.occ.2.5), linetype="dashed")+
  geom_line(aes(x=year, y=n.occ.97.5), linetype="dashed")+
  xlab("Year")+ylab("Number of Occupied Sites")+
  xlim(2016,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021.png",
       plot = p1,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p2<-ggplot(n.occDF, aes(x=year, y=propDetected))+
  geom_point(aes(x=year, y=propDetected), col="blue")+
  geom_line(aes(x=year, y=propDetected), col="blue")+
  xlab("Year")+ylab("Proportion of Visited Sites With Detections")+
  xlim(2016,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021_SitesWithDetections.png",
       plot = p2,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p3<-grid.arrange(p1,p2,nrow=1,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_NumberOccupiedSitesPlusRaw_2016_2023.png",
       plot = p3,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Get Turnover
which(substr(colnames(posterior_draws[[1]]),1,9)=="turnover[")
#[1]8329 8330 8331 8332
colnames(posterior_draws[[1]])[8329]#turnover[1]

turnover1.PostMedian<-smod$quantiles["turnover[1]","50%"]
turnover1.Post.2.5<-smod$quantiles["turnover[1]","2.5%"]
turnover1.Post.97.5<-smod$quantiles["turnover[1]","97.5%"]

turnover2.PostMedian<-smod$quantiles["turnover[2]","50%"]
turnover2.Post.2.5<-smod$quantiles["turnover[2]","2.5%"]
turnover2.Post.97.5<-smod$quantiles["turnover[2]","97.5%"]

turnover3.PostMedian<-smod$quantiles["turnover[3]","50%"]
turnover3.Post.2.5<-smod$quantiles["turnover[3]","2.5%"]
turnover3.Post.97.5<-smod$quantiles["turnover[3]","97.5%"]

turnover4.PostMedian<-smod$quantiles["turnover[4]","50%"]
turnover4.Post.2.5<-smod$quantiles["turnover[4]","2.5%"]
turnover4.Post.97.5<-smod$quantiles["turnover[4]","97.5%"]


turnoverDF<-data.frame(year=c(2017,2018,2019,2021),
                       turnoverMedian=c(turnover1.PostMedian,
                                        turnover2.PostMedian,
                                        turnover3.PostMedian,
                                        turnover4.PostMedian),
                       turnover.2.5=c(turnover1.Post.2.5,
                                      turnover2.Post.2.5,
                                      turnover3.Post.2.5,
                                      turnover4.Post.2.5),
                       turnover.97.5=c(turnover1.Post.97.5,
                                       turnover2.Post.97.5,
                                       turnover3.Post.97.5,
                                       turnover4.Post.97.5))

p4<-ggplot(turnoverDF, aes(x=year, y=turnoverMedian))+
  geom_point(aes(x=year, y=turnoverMedian), col="red")+
  geom_line(aes(x=year, y=turnoverMedian), col="red")+
  geom_line(aes(x=year, y=turnover.2.5), linetype="dashed")+
  geom_line(aes(x=year, y=turnover.97.5), linetype="dashed")+
  xlab("Year")+ylab("Turnover from Previous Year")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021_PredictedTurnover.png",
       plot = p4,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Get AUC

which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
colnames(posterior_draws[[1]])[4408]#"z[1,1]"
colnames(posterior_draws[[1]])[4523]#"z[1,2]"
colnames(posterior_draws[[1]])[4982]#"z[115,5]"

which(substr(colnames(posterior_draws[[1]]),1,4)=="muZ[")
colnames(posterior_draws[[1]])[4983]#"muZ[1,2]"
colnames(posterior_draws[[1]])[5442]#"muZ[115,5]"
str(posterior_draws[[1]][,"muZ[2,2]"])#no estimates for season 1
posterior_draws[[1]][,"muZ[1,2]"]#no estimates for season 1

niter<-ni
nthin<-nt
nchain<-nc
nsite<-115
nyear<-5

library(ggmcmc)
allmuZs<-posterior_draws[1:3][,4983:5442]#
ggd <- ggs(allmuZs)
allZs<-posterior_draws[1:3][,4523:4982]
ggd2 <- ggs(allZs)#starting from year 2

nstore <- niter / nthin
psi.pred <- array(dim=c(nsite, nstore*nchain))
Z.est <- array(dim=c(nsite, nstore*nchain))

#2017
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 2, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7823828
median(AUClist)#[1] 0.7844774
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6973039 0.7844774 0.8556923
#slightly lower AUC scores on average than when latitude
#is included as an effect on psi and phi.

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~97% above 0.7 

#2018
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 3, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.8255165
median(AUClist)#[1] 0.8275978
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7526667 0.8275978 0.8870733
#slightly lower AUC scores on average than when latitude
#is included as an effect on psi and phi.

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~100% above 0.7 

#2019
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 4, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.8120102
median(AUClist)#[1] 0.8137129
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7370000 0.8137129 0.8776487
#slightly lower AUC scores on average than when latitude
#is included as an effect on psi and phi.

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~100% above 0.7 


#2021
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", nyear, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", nyear, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7872093
median(AUClist)#[1] 0.789225
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7066579 0.7892250 0.8555556
#slightly lower AUC scores on average than when latitude
#is included as an effect on psi and phi.

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~99% above 0.7 

#Model is adequate

#Predicted yearly occu. for individual clusters
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
#575 entries = 115 sites * 5 years
library(ggmcmc)
allZs<-posterior_draws[1:3][,4408:4982]#
ggd <- ggs(allZs)

nsite<-115
niter <-60000#instead of 40000
nthin <-10
nchain <-3

nstore <- niter / nthin

#recreate obsFiltFinal dataset but with location and clusterID
#instead of YERA detections
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
#from the whole YERA study area. 2016 is first year
#for which we have Sentinel data.
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


m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits
str(m2)
table(m2$clusterID)#136 sites
# 2_1 2_10  2_8  3_1  3_6  3_8  4_5  4_8  7_5  7_8  9_1
#   3    5   15   13    7   11   45    2   33   21    5

table(m2$clusterIDfinal)
#2_1 2_10 2_8b  3_1  3_6  3_8  4_5  4_8  7_5 
#  3    5   15   13    7   11   45    2   33 
#7_8 7_8b 7_8c 7_8d  9_1 
# 13    1    3    4    5  

#get covariates affecting initial occu.
initFilt<-m2%>%
  #select(locationC, clusterID)%>%
  arrange(locationC)

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
#146 sites with 3+ years of data
#overlapping Hedley layer

#Value:  2  3  4  5  6  7  8 
#nsite:  3  7  3 22 35 36 40

initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11

obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID)
obsFiltFinal$geometry<-NULL

#Probabilities of occu. at year 1 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 1, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2016<-rowMeans(z.pred)
OccEst2016<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2016$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2016bycluster<-OccEst2016%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2016bycluster_quantiles <- OccEst2016bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2016)

#Probabilities of occu. at year 2 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2017<-rowMeans(z.pred)
OccEst2017<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2017$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2017bycluster<-OccEst2017%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2017bycluster_quantiles <- OccEst2017bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2017)


#Probabilities of occu. at year 3 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2018<-rowMeans(z.pred)
OccEst2018<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2018$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2018bycluster<-OccEst2018%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2018bycluster_quantiles <- OccEst2018bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2018)


#Probabilities of occu. at year 4 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2019<-rowMeans(z.pred)
OccEst2019<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2019$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2019bycluster<-OccEst2019%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2019bycluster_quantiles <- OccEst2019bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2019)

#Probabilities of occu. at year 5 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 5, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2021<-rowMeans(z.pred)
OccEst2021<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2021$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2021bycluster<-OccEst2021%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2021bycluster_quantiles <- OccEst2021bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2021)

#meanOccEstAllYears<-cbind(meanOccEst2016,
#                          meanOccEst2017,
#                          meanOccEst2018,
#                          meanOccEst2019,
#                          meanOccEst2021)

finalOccDF<-rbind(OccEst2016bycluster_quantiles,
                  OccEst2017bycluster_quantiles,
                  OccEst2018bycluster_quantiles,
                  OccEst2019bycluster_quantiles,
                  OccEst2021bycluster_quantiles)


#finalOccDF<-cbind(obsFiltFinal, meanOccEstAllYears)
#finalOccDF$clusterID
#write.csv(finalOccDF,file="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021_PredictedOcc.csv")


#Cluster 1 (formerly 3_6)
finalOccDF_cluster1<-finalOccDF%>%
  filter(clusterID=="3_6")

# finalOccDF_cluster1long<-finalOccDF_cluster1%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster1)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_6 <- ggplot(data=finalOccDF_cluster1long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 1: ",nrow(finalOccDF_cluster1)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_6 <- ggplot(data=finalOccDF_cluster1, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 1: ",nrow(obsFiltFinal%>%filter(clusterID=="3_6"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 2
finalOccDF_cluster2<-finalOccDF%>%
  filter(clusterID=="4_5")

# finalOccDF_cluster2long<-finalOccDF_cluster2%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster2)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster4_5 <- ggplot(data=finalOccDF_cluster2long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 2: ",nrow(finalOccDF_cluster2)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster4_5 <- ggplot(data=finalOccDF_cluster2, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 2: ",nrow(obsFiltFinal%>%filter(clusterID=="4_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 3
finalOccDF_cluster3<-finalOccDF%>%
  filter(clusterID=="7_5")

# finalOccDF_cluster3long<-finalOccDF_cluster3%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster3)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

#Occucluster7_5 <- ggplot(data=finalOccDF_cluster3long, aes(x=Year, y=medYear)) +
#  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#  geom_line(aes(x=Year, y=medYear), colour="blue")+
#  my.theme+ggtitle(paste0("Cluster 3: ",nrow(finalOccDF_cluster3)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#  ylim(-0.5,1.25)+xlim(2016,2021)+
#  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster7_5 <- ggplot(data=finalOccDF_cluster3, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 3: ",nrow(obsFiltFinal%>%filter(clusterID=="7_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 4 (just 3 sites at GG Provincial Park)
#No sites with 3 or more years of obs from 2016-2021 in cluster 4
# finalOccDF_cluster4<-finalOccDF%>%
#     filter(clusterID=="2_1")
#   
# finalOccDF_cluster4long<-finalOccDF_cluster4%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster4)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster2_1 <- ggplot(data=finalOccDF_cluster4long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 4: ",nrow(finalOccDF_cluster4)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 5
finalOccDF_cluster5<-finalOccDF%>%
  filter(clusterID=="3_1")

# finalOccDF_cluster5long<-finalOccDF_cluster5%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster5)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_1 <- ggplot(data=finalOccDF_cluster5long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 5: ",nrow(finalOccDF_cluster5)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_1 <- ggplot(data=finalOccDF_cluster5, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 5: ",nrow(obsFiltFinal%>%filter(clusterID=="3_1"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 6
finalOccDF_cluster6<-finalOccDF%>%
  filter(clusterID=="7_8")

# finalOccDF_cluster6long<-finalOccDF_cluster6%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster6)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster7_8 <- ggplot(data=finalOccDF_cluster6long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 6: ",nrow(finalOccDF_cluster6)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


Occucluster7_8 <- ggplot(data=finalOccDF_cluster6, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 6: ",nrow(obsFiltFinal%>%filter(clusterID=="7_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 7
finalOccDF_cluster7<-finalOccDF%>%
  filter(clusterID=="3_8")

# finalOccDF_cluster7long<-finalOccDF_cluster7%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster7)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_8 <- ggplot(data=finalOccDF_cluster7long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 7: ",nrow(finalOccDF_cluster7)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_8 <- ggplot(data=finalOccDF_cluster7, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 7: ",nrow(obsFiltFinal%>%filter(clusterID=="3_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 8
finalOccDF_cluster8<-finalOccDF%>%
  filter(clusterID=="2_8")

# finalOccDF_cluster8long<-finalOccDF_cluster8%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster8)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster2_8 <- ggplot(data=finalOccDF_cluster8long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 8: ",nrow(finalOccDF_cluster8)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster2_8 <- ggplot(data=finalOccDF_cluster8, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 8: ",nrow(obsFiltFinal%>%filter(clusterID=="2_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 9 No sites with at least 3 years of data from 2016-2021 in this cluster
# finalOccDF_cluster9<-finalOccDF%>%
#     filter(clusterID=="4_8")
#   
# finalOccDF_cluster9long<-finalOccDF_cluster9%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster9)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster4_8 <- ggplot(data=finalOccDF_cluster9long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 9: ",nrow(finalOccDF_cluster9)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

p9<-grid.arrange(Occucluster3_6,
                 Occucluster4_5,
                 Occucluster7_5,
                 #Occucluster2_1,
                 Occucluster3_1,
                 Occucluster7_8,
                 Occucluster3_8,
                 Occucluster2_8,
                 #Occucluster4_8,
                 nrow=3,ncol=3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_MedianPlus95BCIoccupancy_2016_2023Dec8.png",
       plot = p9,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Null Model and Predictions
# Specify model in JAGS language
model_code <-
  "model
{
  # Specify priors
  alpha.psi ~ dunif(-3, 3) 

  alpha.gamma ~ dunif(-3, 3)

  alpha.phi ~ dunif(-3, 3)  

  alpha.p ~ dunif(-3, 3) #no detection covariates for now
  
  #random effects on psi, gamma and phi 
  #115 sites, 7 cluster IDs
  for (m in 1:ncluster) {
    u_psi[m] ~ dnorm(0, tau_psi)
    u_phi[m] ~ dnorm(0, tau_phi) #Wait, why dnorm and not dunif?
    u_gamma[m] ~ dnorm(0, tau_gamma)
  }
  tau_psi <- pow(sigma_psi, -2)
  sigma_psi ~ dunif(0, 5)

  tau_phi <- pow(sigma_phi, -2)
  sigma_phi ~ dunif(0, 5)

  tau_gamma <- pow(sigma_gamma, -2)
  sigma_gamma ~ dunif(0, 5)
  

  #assign values generated for u_phi and u_gamma to a vector of size 93 based
  #on clusterID?


  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){
    z0[i] ~dbern(psi1[i]) #True occu. at site i in season 1
    logit(psi1[i])<-alpha.psi + u_psi[cluster[i]]

    z[i,1]<-z0[i]
    for (k in 2:nyear){
      #survival probability from season k-1 to season k
      phi[i,k-1] <-psurv[i,k-1]
      logit(psurv[i,k-1])<-alpha.phi + u_phi[cluster[i]] 

      #colonization probability from season k-1 to season k
      gamma[i,k-1] <-pcol[i,k-1]
      logit(pcol[i,k-1])<-alpha.gamma  +u_gamma[cluster[i]] 

      muZ[i,k]<-z[i,k-1]*phi[i,k-1] +(1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~dbern(muZ[i,k])
    } #k
  } #i
  # Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (k in 1:nyear){
        p[i,j,k]<-det[i,j,k]
        logit(det[i,j,k])<-alpha.p    #(no detection covariates for now)    
        muy[i,j,k]<-z[i,k]*p[i,j,k]
        y[i,j,k]~dbern(muy[i,j,k])
        log_likelihood_contrib[i,j,k]<-y[i,j,k]*log(muy[i,j,k]) + (1-y[i,j,k])*log(1-muy[i,j,k])
      } #k
    } #j
  } #i
  # Derived parameters: Sample and population occu., growth rate and turnover
  n.occ[1]<-sum(z[1:nsite,1])
  for (i in 1:nsite){
    psi[i,1] <-psi1[i]
    for (k in 2:nyear){
      psi[i,k] <-psi[i,k-1]*phi[i,k-1] +(1-psi[i,k-1])*gamma[i,k-1]
    }
  }
  for (k in 2:nyear){
      n.occ[k] <-sum(z[1:nsite,k])
      growthr[k-1]<-mean(psi[1:nsite,k-1])/mean(psi[1:nsite,k-1])
      #no longer growthr[k]<-psi[k]/psi[k-1] or psi[i,k]/psi[i,k-1]
      turnover[k-1] <-(1 - mean(psi[1:nsite,k-1])) *mean(gamma[1:nsite,k-1])/mean(psi[1:nsite,k])
      #turnover[k-1] <-(1 - psi[i,k-1]) *gamma[k-1]/psi[i,k]
  }
  
}"
  
  
# Bundle data
nsite<-115
nyear<-5
nrep<-5
data_jags <-list(y=y,
                   #thr=thr,
                   nsite=nsite,
                   nrep=nrep,
                   nyear=nyear,
                   #psi covariates
                   #site-year covariates (gamma, phi)
                   #cluster information
                   ncluster=7,
                   cluster=clusterID)
#no detection covariates (p)
  
# Initial values
zst <-apply(y,c(1,3),max,na.rm=TRUE)
#Observed occurrence as inits for z
zst[zst==-Inf]<-0
  
#Run model via parallel processing:
initsP = parallel.seeds("base::BaseRNG", 3) # a list of lists
initsP[[1]]$z=zst
initsP[[1]]$alpha.psi=runif(1, -3, 3)
initsP[[1]]$u_psi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.phi=runif(1, -3, 3)
initsP[[1]]$u_phi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.gamma=runif(1, -3, 3)
initsP[[1]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.p=runif(1,-3,3)
  
initsP[[2]]$z=zst
initsP[[2]]$alpha.psi=runif(1, -3, 3)
initsP[[2]]$u_psi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.phi=runif(1, -3, 3)
initsP[[2]]$u_phi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.gamma=runif(1, -3, 3)
initsP[[2]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.p=runif(1,-3,3)
  
initsP[[3]]$z=zst
initsP[[3]]$alpha.psi=runif(1, -3, 3)
initsP[[3]]$u_psi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.phi=runif(1, -3, 3)
initsP[[3]]$u_phi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.gamma=runif(1, -3, 3)
initsP[[3]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.p=runif(1,-3,3)
  
# Parameters monitored
params <-c("psi",
             "alpha.psi",
             "u_psi",
             "phi",
             "alpha.phi",
             "u_phi",
             "gamma",
             "alpha.gamma",
             "u_gamma",
             "alpha.p",
             "p","z","muZ",
             "log_likelihood_contrib",
             "n.occ",
             "growthr",
             "turnover")#,"sens","fpr","auc"
  
  
# Running the model via parallel processing using runjags
ni <-60000#instead of 40000
nt <-10
nb <-50000#instead of 30000
nc <-3
#so there should be 1000 posterior draws
  
niter<-ni
nchains<-nc
thinning<-nt
  
results <- runjags::run.jags(
    model = model_code,
    monitor = params,#keep_track_of,
    data = data_jags,#datalist
    n.chains = nchains,
    inits=initsP,
    adapt = niter/10,
    burnin = nb, 
    sample = niter/thinning, 
    thin = thinning,
    method = "parallel"
)
cat("Model run completed\n")
print(Sys.time())
  
# Extracting results
posterior_draws <- as.mcmc.list(results)
cat("Posterior draws completed\n")
smod <- summary(posterior_draws)
  
save(detectionsLong, 
       data_jags, 
       results, 
       posterior_draws, 
       smod, 
       nc, 
       ni, 
       nb, 
       nt, 
       params, 
       file = "2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_NullModel.Rdata")

load("2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_NullModel.Rdata")  

#Get WAIC
which(substr(colnames(posterior_draws[[1]]),1,23)=="log_likelihood_contrib[") #4375:6649 #There are 2275 (=65 sites*7 years*up to 5 visits per year)
colnames(posterior_draws[[1]])[5431]#"log_likelihood_contrib[1,1,1]"
posterior_draws[[1]][5431]
traceplot(posterior_draws[1:3][,5431], main="log_likelihood_contrib[1,1,1] Trace Plot", xlab="Iteration", ylab="Value", lwd=1)
densplot(posterior_draws[1:3][,"log_likelihood_contrib[1,1,1]"], main="log_likelihood_contrib[1,1,1] Density Plot", xlab="Value", ylab="Density", lwd=1)
smod$quantiles["log_likelihood_contrib[1,1,1]","97.5%"]

library(loo)
#waic(array): An I by C by N array, where I is the number of MCMC iterations per chain,
#C is the number of chains, and N is the number of data points.
#So in my case, a 60000 (or it is 6000=niter/nthin) * 3 * 8484
myArray<-array(NA,
               dim=c(6000,3,2875))

for (i in 1:6000){
  for (j in 1:3){
    for (k in 1:2875){
      try(myArray[i,j,k]<-posterior_draws[[j]][i,(k+5430)])
    }
  }
}

waicMyArray<-waic(myArray)
#Computed from 18000 by 2875 log-likelihood matrix.
#Estimate   SE
#elpd_waic   -880.5 16.6
#p_waic       105.9  2.9
#waic        1760.9 33.2
looMyArray<-loo(myArray)
#Computed from 18000 by 2875 log-likelihood matrix.
#Estimate   SE
#elpd_loo   -880.4 16.6
#p_loo       105.8  2.8
#looic      1760.9 33.2
#------
#MCSE of elpd_loo is NA.
#MCSE and ESS estimates assume independent draws (r_eff=1).
#Pareto k diagnostic values:
#Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     1740  60.5%   15098   
#(0.7, 1]   (bad)         0   0.0%   <NA>    
#(1, Inf)   (very bad) 1135  39.5%   <NA>    
#See help('pareto-k-diagnostic') for details.

#Predicted yearly occu. for individual clusters
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
#575 entries = 115 sites * 5 years
library(ggmcmc)
allZs<-posterior_draws[1:3][,4396:4970]#
ggd <- ggs(allZs)

nsite<-115
niter <-60000#instead of 40000
nthin <-10
nchain <-3

nstore <- niter / nthin

#recreate obsFiltFinal dataset but with location and clusterID
#instead of YERA detections
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
#from the whole YERA study area. 2016 is first year
#for which we have Sentinel data.
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


m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits
str(m2)
table(m2$clusterID)#136 sites
# 2_1 2_10  2_8  3_1  3_6  3_8  4_5  4_8  7_5  7_8  9_1
#   3    5   15   13    7   11   45    2   33   21    5

table(m2$clusterIDfinal)
#2_1 2_10 2_8b  3_1  3_6  3_8  4_5  4_8  7_5 
#  3    5   15   13    7   11   45    2   33 
#7_8 7_8b 7_8c 7_8d  9_1 
# 13    1    3    4    5  

#get covariates affecting initial occu.
initFilt<-m2%>%
  #select(locationC, clusterID)%>%
  arrange(locationC)

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

#Value:  2  3  4  5  6  7  8 
#nsite:  3  7  3 22 35 36 40

initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11

obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID)
obsFiltFinal$geometry<-NULL

#Probabilities of occu. at year 1 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 1, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2016<-rowMeans(z.pred)
OccEst2016<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2016$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2016bycluster<-OccEst2016%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2016bycluster_quantiles <- OccEst2016bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2016)

#Probabilities of occu. at year 2 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2017<-rowMeans(z.pred)
OccEst2017<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2017$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2017bycluster<-OccEst2017%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2017bycluster_quantiles <- OccEst2017bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2017)


#Probabilities of occu. at year 3 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2018<-rowMeans(z.pred)
OccEst2018<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2018$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2018bycluster<-OccEst2018%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2018bycluster_quantiles <- OccEst2018bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2018)


#Probabilities of occu. at year 4 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2019<-rowMeans(z.pred)
OccEst2019<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2019$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2019bycluster<-OccEst2019%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2019bycluster_quantiles <- OccEst2019bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2019)

#Probabilities of occu. at year 5 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 5, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2021<-rowMeans(z.pred)
OccEst2021<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2021$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2021bycluster<-OccEst2021%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2021bycluster_quantiles <- OccEst2021bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2021)

#meanOccEstAllYears<-cbind(meanOccEst2016,
#                          meanOccEst2017,
#                          meanOccEst2018,
#                          meanOccEst2019,
#                          meanOccEst2021)

finalOccDF<-rbind(OccEst2016bycluster_quantiles,
                  OccEst2017bycluster_quantiles,
                  OccEst2018bycluster_quantiles,
                  OccEst2019bycluster_quantiles,
                  OccEst2021bycluster_quantiles)


#finalOccDF<-cbind(obsFiltFinal, meanOccEstAllYears)
#finalOccDF$clusterID
#write.csv(finalOccDF,file="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021_PredictedOcc.csv")


#Cluster 1 (formerly 3_6)
finalOccDF_cluster1<-finalOccDF%>%
  filter(clusterID=="3_6")

# finalOccDF_cluster1long<-finalOccDF_cluster1%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster1)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_6 <- ggplot(data=finalOccDF_cluster1long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 1: ",nrow(finalOccDF_cluster1)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_6 <- ggplot(data=finalOccDF_cluster1, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 1: ",nrow(obsFiltFinal%>%filter(clusterID=="3_6"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 2
finalOccDF_cluster2<-finalOccDF%>%
  filter(clusterID=="4_5")

# finalOccDF_cluster2long<-finalOccDF_cluster2%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster2)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster4_5 <- ggplot(data=finalOccDF_cluster2long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 2: ",nrow(finalOccDF_cluster2)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster4_5 <- ggplot(data=finalOccDF_cluster2, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 2: ",nrow(obsFiltFinal%>%filter(clusterID=="4_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 3
finalOccDF_cluster3<-finalOccDF%>%
  filter(clusterID=="7_5")

# finalOccDF_cluster3long<-finalOccDF_cluster3%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster3)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

#Occucluster7_5 <- ggplot(data=finalOccDF_cluster3long, aes(x=Year, y=medYear)) +
#  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#  geom_line(aes(x=Year, y=medYear), colour="blue")+
#  my.theme+ggtitle(paste0("Cluster 3: ",nrow(finalOccDF_cluster3)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#  ylim(-0.5,1.25)+xlim(2016,2021)+
#  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster7_5 <- ggplot(data=finalOccDF_cluster3, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 3: ",nrow(obsFiltFinal%>%filter(clusterID=="7_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 4 (just 3 sites at GG Provincial Park)
#No sites with 3 or more years of obs from 2016-2021 in cluster 4
# finalOccDF_cluster4<-finalOccDF%>%
#     filter(clusterID=="2_1")
#   
# finalOccDF_cluster4long<-finalOccDF_cluster4%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster4)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster2_1 <- ggplot(data=finalOccDF_cluster4long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 4: ",nrow(finalOccDF_cluster4)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 5
finalOccDF_cluster5<-finalOccDF%>%
  filter(clusterID=="3_1")

# finalOccDF_cluster5long<-finalOccDF_cluster5%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster5)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_1 <- ggplot(data=finalOccDF_cluster5long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 5: ",nrow(finalOccDF_cluster5)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_1 <- ggplot(data=finalOccDF_cluster5, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 5: ",nrow(obsFiltFinal%>%filter(clusterID=="3_1"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 6
finalOccDF_cluster6<-finalOccDF%>%
  filter(clusterID=="7_8")

# finalOccDF_cluster6long<-finalOccDF_cluster6%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster6)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster7_8 <- ggplot(data=finalOccDF_cluster6long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 6: ",nrow(finalOccDF_cluster6)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


Occucluster7_8 <- ggplot(data=finalOccDF_cluster6, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 6: ",nrow(obsFiltFinal%>%filter(clusterID=="7_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 7
finalOccDF_cluster7<-finalOccDF%>%
  filter(clusterID=="3_8")

# finalOccDF_cluster7long<-finalOccDF_cluster7%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster7)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_8 <- ggplot(data=finalOccDF_cluster7long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 7: ",nrow(finalOccDF_cluster7)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_8 <- ggplot(data=finalOccDF_cluster7, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 7: ",nrow(obsFiltFinal%>%filter(clusterID=="3_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 8
finalOccDF_cluster8<-finalOccDF%>%
  filter(clusterID=="2_8")

# finalOccDF_cluster8long<-finalOccDF_cluster8%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster8)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster2_8 <- ggplot(data=finalOccDF_cluster8long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 8: ",nrow(finalOccDF_cluster8)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster2_8 <- ggplot(data=finalOccDF_cluster8, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 8: ",nrow(obsFiltFinal%>%filter(clusterID=="2_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 9 No sites with at least 3 years of data from 2016-2021 in this cluster
# finalOccDF_cluster9<-finalOccDF%>%
#     filter(clusterID=="4_8")
#   
# finalOccDF_cluster9long<-finalOccDF_cluster9%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster9)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster4_8 <- ggplot(data=finalOccDF_cluster9long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 9: ",nrow(finalOccDF_cluster9)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

p9<-grid.arrange(Occucluster3_6,
                 Occucluster4_5,
                 Occucluster7_5,
                 #Occucluster2_1,
                 Occucluster3_1,
                 Occucluster7_8,
                 Occucluster3_8,
                 Occucluster2_8,
                 #Occucluster4_8,
                 nrow=3,ncol=3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/NullModel_MedianPlus95BCIoccupancy_2016_2023Dec8.png",
       plot = p9,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Get Estimated Number of Occupied Sites (out of 115)
which(substr(colnames(posterior_draws[[1]]),1,6)=="n.occ[")
# [1] 8306 8307 8308 8309 8310
colnames(posterior_draws[[1]])[8306]#n.occ[1]

n.occ1.PostMedian<-smod$quantiles["n.occ[1]","50%"]
n.occ1.Post.2.5<-smod$quantiles["n.occ[1]","2.5%"]
n.occ1.Post.97.5<-smod$quantiles["n.occ[1]","97.5%"]

n.occ2.PostMedian<-smod$quantiles["n.occ[2]","50%"]
n.occ2.Post.2.5<-smod$quantiles["n.occ[2]","2.5%"]
n.occ2.Post.97.5<-smod$quantiles["n.occ[2]","97.5%"]

n.occ3.PostMedian<-smod$quantiles["n.occ[3]","50%"]
n.occ3.Post.2.5<-smod$quantiles["n.occ[3]","2.5%"]
n.occ3.Post.97.5<-smod$quantiles["n.occ[3]","97.5%"]

n.occ4.PostMedian<-smod$quantiles["n.occ[4]","50%"]
n.occ4.Post.2.5<-smod$quantiles["n.occ[4]","2.5%"]
n.occ4.Post.97.5<-smod$quantiles["n.occ[4]","97.5%"]

n.occ5.PostMedian<-smod$quantiles["n.occ[5]","50%"]
n.occ5.Post.2.5<-smod$quantiles["n.occ[5]","2.5%"]
n.occ5.Post.97.5<-smod$quantiles["n.occ[5]","97.5%"]

#need to know #sites visited in each year:
nrow(data_jags$y[,,1])
year1<-data.frame(data_jags$y[,,1])
nrow(year1[!is.na(rowSums(year1)),])#82 
visited_year1<-year1[!is.na(rowSums(year1)),]#82
nrow(visited_year1[(rowSums(visited_year1)>0),])#47 

year2<-data.frame(data_jags$y[,,2])
nrow(year2[!is.na(rowSums(year2)),])#101
visited_year2<-year2[!is.na(rowSums(year2)),]#101
nrow(visited_year2[(rowSums(visited_year2)>0),])#58

year3<-data.frame(data_jags$y[,,3])
nrow(year3[!is.na(rowSums(year3)),])#62
visited_year3<-year3[!is.na(rowSums(year3)),]#62
nrow(visited_year3[(rowSums(visited_year3)>0),])#12

year4<-data.frame(data_jags$y[,,4])
nrow(year4[!is.na(rowSums(year4)),])#84
visited_year4<-year4[!is.na(rowSums(year4)),]#84
nrow(visited_year4[(rowSums(visited_year4)>0),])#30

year5<-data.frame(data_jags$y[,,5])
nrow(year5[!is.na(rowSums(year5)),])#79
visited_year5<-year5[!is.na(rowSums(year5)),]#79
nrow(visited_year5[(rowSums(visited_year5)>0),])#25

n.occDF<-data.frame(year=c(2016,2017,2018,2019,2021),
                    nsitesVisited=c(82,101,62,84,79),
                    nsitesDetected=c(47,58,12,30,25),
                    n.occMedian=c(n.occ1.PostMedian,
                                  n.occ2.PostMedian,
                                  n.occ3.PostMedian,
                                  n.occ4.PostMedian,
                                  n.occ5.PostMedian),
                    n.occ.2.5=c(n.occ1.Post.2.5,
                                n.occ2.Post.2.5,
                                n.occ3.Post.2.5,
                                n.occ4.Post.2.5,
                                n.occ5.Post.2.5),
                    n.occ.97.5=c(n.occ1.Post.97.5,
                                 n.occ2.Post.97.5,
                                 n.occ3.Post.97.5,
                                 n.occ4.Post.97.5,
                                 n.occ5.Post.97.5))

cor(n.occDF$nsitesVisited, n.occDF$n.occMedian)
#[1] 0.714214
n.occDF$propDetected<-n.occDF$nsitesDetected/n.occDF$nsitesVisited

p1<-ggplot(n.occDF, aes(x=year, y=n.occMedian))+
  geom_point(aes(x=year, y=nsitesVisited), col="blue")+
  geom_line(aes(x=year, y=nsitesVisited), col="blue")+
  geom_point(aes(x=year, y=n.occMedian), col="red")+
  geom_line(aes(x=year, y=n.occMedian), col="red")+
  geom_line(aes(x=year, y=n.occ.2.5), linetype="dashed")+
  geom_line(aes(x=year, y=n.occ.97.5), linetype="dashed")+
  xlab("Year")+ylab("Number of Occupied Sites")+
  xlim(2016,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/NullModel_2016_2021.png",
       plot = p1,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p2<-ggplot(n.occDF, aes(x=year, y=propDetected))+
  geom_point(aes(x=year, y=propDetected), col="blue")+
  geom_line(aes(x=year, y=propDetected), col="blue")+
  xlab("Year")+ylab("Proportion of Visited Sites With Detections")+
  xlim(2016,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/NullModel_2016_2021_SitesWithDetections.png",
       plot = p2,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p3<-grid.arrange(p1,p2,nrow=1,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/NullModel_NumberOccupiedSitesPlusRaw_2016_2023.png",
       plot = p3,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Get Turnover
which(substr(colnames(posterior_draws[[1]]),1,9)=="turnover[")
#[1]8315 8316 8317 8318
colnames(posterior_draws[[1]])[8315]#turnover[1]

turnover1.PostMedian<-smod$quantiles["turnover[1]","50%"]
turnover1.Post.2.5<-smod$quantiles["turnover[1]","2.5%"]
turnover1.Post.97.5<-smod$quantiles["turnover[1]","97.5%"]

turnover2.PostMedian<-smod$quantiles["turnover[2]","50%"]
turnover2.Post.2.5<-smod$quantiles["turnover[2]","2.5%"]
turnover2.Post.97.5<-smod$quantiles["turnover[2]","97.5%"]

turnover3.PostMedian<-smod$quantiles["turnover[3]","50%"]
turnover3.Post.2.5<-smod$quantiles["turnover[3]","2.5%"]
turnover3.Post.97.5<-smod$quantiles["turnover[3]","97.5%"]

turnover4.PostMedian<-smod$quantiles["turnover[4]","50%"]
turnover4.Post.2.5<-smod$quantiles["turnover[4]","2.5%"]
turnover4.Post.97.5<-smod$quantiles["turnover[4]","97.5%"]


turnoverDF<-data.frame(year=c(2017,2018,2019,2021),
                       turnoverMedian=c(turnover1.PostMedian,
                                        turnover2.PostMedian,
                                        turnover3.PostMedian,
                                        turnover4.PostMedian),
                       turnover.2.5=c(turnover1.Post.2.5,
                                      turnover2.Post.2.5,
                                      turnover3.Post.2.5,
                                      turnover4.Post.2.5),
                       turnover.97.5=c(turnover1.Post.97.5,
                                       turnover2.Post.97.5,
                                       turnover3.Post.97.5,
                                       turnover4.Post.97.5))

p4<-ggplot(turnoverDF, aes(x=year, y=turnoverMedian))+
  geom_point(aes(x=year, y=turnoverMedian), col="red")+
  geom_line(aes(x=year, y=turnoverMedian), col="red")+
  geom_line(aes(x=year, y=turnover.2.5), linetype="dashed")+
  geom_line(aes(x=year, y=turnover.97.5), linetype="dashed")+
  xlab("Year")+ylab("Turnover from Previous Year")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/NullModel_2016_2021_PredictedTurnover.png",
       plot = p4,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Get AUC

which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
colnames(posterior_draws[[1]])[4396]#"z[1,1]"
colnames(posterior_draws[[1]])[4511]#"z[1,2]"
colnames(posterior_draws[[1]])[4970]#"z[115,5]"

which(substr(colnames(posterior_draws[[1]]),1,4)=="muZ[")
colnames(posterior_draws[[1]])[4971]#"muZ[1,2]"
colnames(posterior_draws[[1]])[5430]#"muZ[115,5]"
str(posterior_draws[[1]][,"muZ[2,2]"])#no estimates for season 1
posterior_draws[[1]][,"muZ[1,2]"]#no estimates for season 1

niter<-ni
nthin<-nt
nchain<-nc
nsite<-115
nyear<-5

library(ggmcmc)
allmuZs<-posterior_draws[1:3][,4971:5430]#
ggd <- ggs(allmuZs)
allZs<-posterior_draws[1:3][,4511:4970]
ggd2 <- ggs(allZs)#starting from year 2

nstore <- niter / nthin
psi.pred <- array(dim=c(nsite, nstore*nchain))
Z.est <- array(dim=c(nsite, nstore*nchain))

#2017
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 2, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.5963986
median(AUClist)#[1] 0.5967692
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.4921857 0.5967692 0.6999422
#much lower AUC scores on average 

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~2.5% above 0.7 

#2018
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 3, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7227006
median(AUClist)#[1] 0.7248787
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6200994 0.7248787 0.8128576 
#much lower AUC scores on average 

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~69% above 0.7 

#2019
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 4, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.767458
median(AUClist)#[1] 0.771483
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6570364 0.7714830 0.8543352 
#much lower AUC scores on average 

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~95% above 0.7 

 
#2021
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", nyear, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", nyear, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7088674 instead of 0.7872093
median(AUClist)#[1] 0.7137871 instead of 0.789225
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.5932810 0.7137871 0.8029030
#much lower AUC scores on average than when latitude
#is included as an effect on psi and phi.

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~60% above 0.7 

#Model is adequate

#RS hybrid model for whole study area, minus Hedley but include latitude

#cluster ID for random effect
clusterID<-as.integer(umf@siteCovs$clusterID)
#reconfigure site and site-year covariates as dummy variable matrices
class2<-umf@siteCovs$HabSuit4ClassNumeric
class2<-ifelse(class2==2,1,0)
class3<-umf@siteCovs$HabSuit4ClassNumeric
class3<-ifelse(class3==3,1,0)
class4<-umf@siteCovs$HabSuit4ClassNumeric
class4<-ifelse(class4==4,1,0)

#class2SY<-umf@yearlySiteCovs$habSuit
#class2SY<-ifelse(class2SY==2,1,0)
#class3SY<-umf@yearlySiteCovs$habSuit
#class3SY<-ifelse(class3SY==2,1,0)
#class4SY<-umf@yearlySiteCovs$habSuit
#class4SY<-ifelse(class4SY==2,1,0)

class2SY<-as.matrix(cbind(class2,class2,class2,class2,class2))
class3SY<-as.matrix(cbind(class3,class3,class3,class3,class3))
class4SY<-as.matrix(cbind(class4,class4,class4,class4,class4))

#make y into a 3-dimensional array
obsFiltFinal.jags<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC,yera2016_1,yera2016_2,yera2016_3,yera2016_4,yera2016_5,
         yera2017_1,yera2017_2,yera2017_3,yera2017_4,yera2017_5,
         yera2018_1,yera2018_2,yera2018_3,yera2018_4,yera2018_5,
         yera2019_1,yera2019_2,yera2019_3,yera2019_4,yera2019_5,
         yera2021_1,yera2021_2,yera2021_3,yera2021_4,yera2021_5)
obsFiltFinal.jags$geometry<-NULL

detectionsLong<-obsFiltFinal.jags%>%
  pivot_longer(cols=yera2016_1:yera2021_5, names_to="year_round", values_to="detected")%>%
  mutate(year=substr(year_round,5,8))%>%
  mutate(round=substr(year_round,10,10))%>%
  arrange(locationC,year,round)
nlevels(as.factor(detectionsLong$locationC))
#115 sites

#get detections and put in 3-d array
R <- nlevels(as.factor(detectionsLong$locationC))#146 sites
T <- nlevels(as.factor(detectionsLong$round))#up to 5 visits per site per year
S <- nlevels(as.factor(detectionsLong$year))#8 years (no visits to any sites in 2020 due to COVID)

nsite<-R
nrep<-T
nyear<-S

sites<-levels(as.factor(detectionsLong$locationC))
rounds<-levels(as.factor(detectionsLong$round))
years<-levels(as.factor(detectionsLong$year))

y <-array(NA,
          dim=c(nsite,nrep,nyear),
          dimnames=list(sites, rounds, years)
)
for (i in sites){
  for (j in rounds){
    for (k in years){
      try(placeHolder<-detectionsLong%>%
            filter(locationC==i)%>%
            filter(round==j)%>%
            filter(year==k)%>%
            select(detected))
      try(y[i,j,k]<-placeHolder$detected) 
    }
  }
}
str(y)

#NO DETECTION SUBMODEL COVARIATES IN THIS ANALYSIS
#THEREFORE NO DETECTION ARRAYS NEEDED

#INCLUDING LATITUDE BUT NOT HEDLEY LAYER
#CONSIDER INCLUDING CLUSTER 4 (GIPSY CORDON PROVINCIAL PARK)
#####################################################################################################################
#                                                                                                                   #
#                                                                                                                   #
#                       Create Model File and Initialize Parameters - Modelling Individual Footprint Types          #
#                                                                                                                   #
#                                                                                                                   #
#                                                                                                                   #
#                                                                                                                   #
#####################################################################################################################

ndmi<-umf@siteCovs$sNDMI.s[,1]
ndmi2<-umf@siteCovs$sNDMI.s2[,1]
deltaVH<-umf@siteCovs$DeltaVH.s[,1]
sgrass<-umf@siteCovs$sgrass.s[,1]
sgrass2<-umf@siteCovs$sgrass.s2[,1]
latsqrt<-umf@siteCovs$lat.sqrt

Yswi<-scale(Y_swi, center=TRUE, scale=TRUE)
YdeltaVH<-scale(yearlyDeltaVH[,2:6], center=TRUE, scale=TRUE)
YdeltaVH2<-YdeltaVH^2

sDPOL<-as.matrix(cbind(sDPOLsummer,sDPOLsummer,sDPOLsummer,
                       sDPOLsummer,sDPOLsummer))
YsDPOLsummer<-scale(sDPOL, center=TRUE, scale=TRUE)

Ylatsqrt<-scale(sqrt(yearlyLatitude), center=TRUE, scale=TRUE)

#Run multi-season occu. model in JAGS
# Specify model in JAGS language
model_code <-
  "model
{
  # Specify priors
  alpha.psi ~ dunif(-3, 3) 
  beta1.psi ~ dunif(-3, 3)
  beta2.psi ~ dunif(-3, 3)
  beta3.psi ~ dunif(-3, 3) 
  beta4.psi ~ dunif(-3, 3)
  beta5.psi ~ dunif(-3, 3)
  beta6.psi ~ dunif(-3, 3) 

  alpha.gamma ~ dunif(-3, 3)
  beta1.gamma ~ dunif(-3, 3) 

  alpha.phi ~ dunif(-3, 3)  
  beta1.phi ~ dunif(-3, 3)
  beta2.phi ~ dunif(-3, 3)
  beta3.phi ~ dunif(-3, 3)
  beta4.phi ~ dunif(-3, 3)

  alpha.p ~ dunif(-3, 3) #no detection covariates for now
  
  #random effects on psi, gamma and phi 
  #115 sites, 7 cluster IDs
  for (m in 1:ncluster) {
    u_psi[m] ~ dnorm(0, tau_psi)
    u_phi[m] ~ dnorm(0, tau_phi) #Wait, why dnorm and not dunif?
    u_gamma[m] ~ dnorm(0, tau_gamma)
  }
  tau_psi <- pow(sigma_psi, -2)
  sigma_psi ~ dunif(0, 5)

  tau_phi <- pow(sigma_phi, -2)
  sigma_phi ~ dunif(0, 5)

  tau_gamma <- pow(sigma_gamma, -2)
  sigma_gamma ~ dunif(0, 5)
  

  #assign values generated for u_phi and u_gamma to a vector of size 93 based
  #on clusterID?


  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){
    z0[i] ~dbern(psi1[i]) #True occu. at site i in season 1
    logit(psi1[i])<-alpha.psi +beta1.psi*ndmi[i] +beta2.psi*ndmi2[i] +beta3.psi*deltaVH[i] +beta4.psi*sgrass[i] +beta5.psi*sgrass2[i] +beta6.psi*latsqrt[i] + u_psi[cluster[i]]

    z[i,1]<-z0[i]
    for (k in 2:nyear){
      #survival probability from season k-1 to season k
      phi[i,k-1] <-psurv[i,k-1]
      logit(psurv[i,k-1])<-alpha.phi + u_phi[cluster[i]] +beta1.phi*YdeltaVH[i,k-1] +beta2.phi*YdeltaVH2[i,k-1] +beta3.phi*YsDPOLsummer[i,k-1] +beta4.phi*Ylatsqrt[i,k-1] 

      #colonization probability from season k-1 to season k
      gamma[i,k-1] <-pcol[i,k-1]
      logit(pcol[i,k-1])<-alpha.gamma  +u_gamma[cluster[i]] +beta1.gamma*Yswi[i,k-1] 

      muZ[i,k]<-z[i,k-1]*phi[i,k-1] +(1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~dbern(muZ[i,k])
    } #k
  } #i
  # Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (k in 1:nyear){
        p[i,j,k]<-det[i,j,k]
        logit(det[i,j,k])<-alpha.p    #(no detection covariates for now)    
        muy[i,j,k]<-z[i,k]*p[i,j,k]
        y[i,j,k]~dbern(muy[i,j,k])
        log_likelihood_contrib[i,j,k]<-y[i,j,k]*log(muy[i,j,k]) + (1-y[i,j,k])*log(1-muy[i,j,k])
      } #k
    } #j
  } #i
  # Derived parameters: Sample and population occu., growth rate and turnover
  n.occ[1]<-sum(z[1:nsite,1])
  for (i in 1:nsite){
    psi[i,1] <-psi1[i]
    for (k in 2:nyear){
      psi[i,k] <-psi[i,k-1]*phi[i,k-1] +(1-psi[i,k-1])*gamma[i,k-1]
    }
  }
  for (k in 2:nyear){
      n.occ[k] <-sum(z[1:nsite,k])
      growthr[k-1]<-mean(psi[1:nsite,k-1])/mean(psi[1:nsite,k-1])
      #no longer growthr[k]<-psi[k]/psi[k-1] or psi[i,k]/psi[i,k-1]
      turnover[k-1] <-(1 - mean(psi[1:nsite,k-1])) *mean(gamma[1:nsite,k-1])/mean(psi[1:nsite,k])
      #turnover[k-1] <-(1 - psi[i,k-1]) *gamma[k-1]/psi[i,k]
  }
  
}"
  
  
# Bundle data
nsite<-115
nyear<-5
nrep<-5
data_jags <-list(y=y,
                   #thr=thr,
                   nsite=nsite,
                   nrep=nrep,
                   nyear=nyear,
                   #psi covariates
                   ndmi=ndmi,
                   ndmi2=ndmi2,
                   deltaVH=deltaVH,
                   sgrass=sgrass,
                   sgrass2=sgrass2,
                   latsqrt=latsqrt,
                   #site-year covariates (gamma, phi)
                   Yswi=Yswi,
                   YdeltaVH=YdeltaVH,
                   YdeltaVH2=YdeltaVH2,
                   YsDPOLsummer=YsDPOLsummer,
                   Ylatsqrt=Ylatsqrt,
                   #cluster information
                   ncluster=7,
                   cluster=clusterID)
#no detection covariates (p)
  
# Initial values
zst <-apply(y,c(1,3),max,na.rm=TRUE)
#Observed occurrence as inits for z
zst[zst==-Inf]<-0
  
#Run model via parallel processing:
initsP = parallel.seeds("base::BaseRNG", 3) # a list of lists
initsP[[1]]$z=zst
initsP[[1]]$alpha.psi=runif(1, -3, 3)
initsP[[1]]$beta1.psi=runif(1, -3, 3)
initsP[[1]]$beta2.psi=runif(1, -3, 3)
initsP[[1]]$beta3.psi=runif(1, -3, 3)
initsP[[1]]$beta4.psi=runif(1, -3, 3)
initsP[[1]]$beta5.psi=runif(1, -3, 3)
initsP[[1]]$beta6.psi=runif(1, -3, 3)
initsP[[1]]$u_psi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.phi=runif(1, -3, 3)
initsP[[1]]$beta1.phi =runif(1, -3, 3)
initsP[[1]]$beta2.phi =runif(1, -3, 3)
initsP[[1]]$beta3.phi =runif(1, -3, 3)
initsP[[1]]$beta4.phi =runif(1, -3, 3)
initsP[[1]]$u_phi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.gamma=runif(1, -3, 3)
initsP[[1]]$beta1.gamma =runif(1, -3, 3)
initsP[[1]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.p=runif(1,-3,3)
  
initsP[[2]]$z=zst
initsP[[2]]$alpha.psi=runif(1, -3, 3)
initsP[[2]]$beta1.psi=runif(1, -3, 3)
initsP[[2]]$beta2.psi=runif(1, -3, 3)
initsP[[2]]$beta3.psi=runif(1, -3, 3)
initsP[[2]]$beta4.psi=runif(1, -3, 3)
initsP[[2]]$beta5.psi=runif(1, -3, 3)
initsP[[2]]$beta6.psi=runif(1, -3, 3)
initsP[[2]]$u_psi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.phi=runif(1, -3, 3)
initsP[[2]]$beta1.phi =runif(1, -3, 3)
initsP[[2]]$beta2.phi =runif(1, -3, 3)
initsP[[2]]$beta3.phi =runif(1, -3, 3)
initsP[[2]]$beta4.phi =runif(1, -3, 3)
initsP[[2]]$u_phi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.gamma=runif(1, -3, 3)
initsP[[2]]$beta1.gamma =runif(1, -3, 3)
initsP[[2]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.p=runif(1,-3,3)
  
initsP[[3]]$z=zst
initsP[[3]]$alpha.psi=runif(1, -3, 3)
initsP[[3]]$beta1.psi=runif(1, -3, 3)
initsP[[3]]$beta2.psi=runif(1, -3, 3)
initsP[[3]]$beta3.psi=runif(1, -3, 3)
initsP[[3]]$beta4.psi=runif(1, -3, 3)
initsP[[3]]$beta5.psi=runif(1, -3, 3)
initsP[[3]]$beta6.psi=runif(1, -3, 3)
initsP[[3]]$u_psi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.phi=runif(1, -3, 3)
initsP[[3]]$beta1.phi =runif(1, -3, 3)
initsP[[3]]$beta2.phi =runif(1, -3, 3)
initsP[[3]]$beta3.phi =runif(1, -3, 3)
initsP[[3]]$beta4.phi =runif(1, -3, 3)
initsP[[3]]$u_phi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.gamma=runif(1, -3, 3)
initsP[[3]]$beta1.gamma =runif(1, -3, 3)
initsP[[3]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.p=runif(1,-3,3)
  
# Parameters monitored
params <-c("psi",
             "alpha.psi",
             "beta1.psi",
             "beta2.psi",
             "beta3.psi",
             "beta4.psi",
             "beta5.psi",
             "beta6.psi",
             "u_psi",
             "phi",
             "alpha.phi",
             "beta1.phi",
             "beta2.phi",
             "beta3.phi",
             "beta4.phi",
             "u_phi",
             "gamma",
             "alpha.gamma",
             "beta1.gamma",
             "u_gamma",
             "alpha.p",
             "p","z","muZ",
             "log_likelihood_contrib",
             "n.occ",
             "growthr",
             "turnover")#,"sens","fpr","auc"
  
  
# Running the model via parallel processing using runjags
ni <-60000#instead of 40000
nt <-10
nb <-50000#instead of 30000
nc <-3
#so there should be 1000 posterior draws
  
niter<-ni
nchains<-nc
thinning<-nt
  
results <- runjags::run.jags(
    model = model_code,
    monitor = params,#keep_track_of,
    data = data_jags,#datalist
    n.chains = nchains,
    inits=initsP,
    adapt = niter/10,
    burnin = nb, 
    sample = niter/thinning, 
    thin = thinning,
    method = "parallel"
  )
cat("Model run completed\n")
print(Sys.time())
  
# Extracting results
posterior_draws <- as.mcmc.list(results)
cat("Posterior draws completed\n")
smod <- summary(posterior_draws)
  
save(detectionsLong, 
       data_jags, 
       results, 
       posterior_draws, 
       smod, 
       nc, 
       ni, 
       nb, 
       nt, 
       params, 
       file = "2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_HybridWithLatMinusHedley.Rdata")

load("2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/JAGS_2016_2021_HybridWithLatMinusHedley.Rdata")  
int.mcmc.mat <- as.matrix(posterior_draws)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)

#Get WAIC
which(substr(colnames(posterior_draws[[1]]),1,23)=="log_likelihood_contrib[") #4375:6649 #There are 2275 (=65 sites*7 years*up to 5 visits per year)
colnames(posterior_draws[[1]])[5442]#"log_likelihood_contrib[1,1,1]"
posterior_draws[[1]][5442]
traceplot(posterior_draws[1:3][,5442], main="log_likelihood_contrib[1,1,1] Trace Plot", xlab="Iteration", ylab="Value", lwd=1)
densplot(posterior_draws[1:3][,"log_likelihood_contrib[1,1,1]"], main="log_likelihood_contrib[1,1,1] Density Plot", xlab="Value", ylab="Density", lwd=1)
smod$quantiles["log_likelihood_contrib[1,1,1]","97.5%"]

library(loo)
#waic(array): An I by C by N array, where I is the number of MCMC iterations per chain,
#C is the number of chains, and N is the number of data points.
#So in my case, a 60000 (or it is 6000=niter/nthin) * 3 * 8484
myArray<-array(NA,
               dim=c(6000,3,2875))

for (i in 1:6000){
  for (j in 1:3){
    for (k in 1:2875){
      try(myArray[i,j,k]<-posterior_draws[[j]][i,(k+5441)])
    }
  }
}
waicMyArray<-waic(myArray)
#Computed from 18000 by 2875 log-likelihood matrix.
#Estimate   SE
#elpd_waic   -890.9 17.0
#p_waic        91.7  2.6
#waic        1781.9 34.1
looMyArray<-loo(myArray)
#Computed from 18000 by 2875 log-likelihood matrix.
#elpd_loo   -891.0 17.0
#p_loo        91.8  2.6
#looic      1782.0 34.0
#------
#MCSE of elpd_loo is NA.
#MCSE and ESS estimates assume independent draws (r_eff=1).
#Pareto k diagnostic values:
#Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     1770  61.6%   14530   
#(0.7, 1]   (bad)         0   0.0%   <NA>    
#(1, Inf)   (very bad) 1105  38.4%   <NA>    
#See help('pareto-k-diagnostic') for details.

alpha.psi.PostMedian<-smod$quantiles["alpha.psi","50%"]
alpha.psi.Post.2.5<-smod$quantiles["alpha.psi","2.5%"]
alpha.psi.Post.97.5<-smod$quantiles["alpha.psi","97.5%"]

beta1.psi.PostMedian<-smod$quantiles["beta1.psi","50%"]
beta1.psi.Post.2.5<-smod$quantiles["beta1.psi","2.5%"]
beta1.psi.Post.97.5<-smod$quantiles["beta1.psi","97.5%"]

beta2.psi.PostMedian<-smod$quantiles["beta2.psi","50%"]
beta2.psi.Post.2.5<-smod$quantiles["beta2.psi","2.5%"]
beta2.psi.Post.97.5<-smod$quantiles["beta2.psi","97.5%"]

beta3.psi.PostMedian<-smod$quantiles["beta3.psi","50%"]
beta3.psi.Post.2.5<-smod$quantiles["beta3.psi","2.5%"]
beta3.psi.Post.97.5<-smod$quantiles["beta3.psi","97.5%"]

beta4.psi.PostMedian<-smod$quantiles["beta4.psi","50%"]
beta4.psi.Post.2.5<-smod$quantiles["beta4.psi","2.5%"]
beta4.psi.Post.97.5<-smod$quantiles["beta4.psi","97.5%"]

beta5.psi.PostMedian<-smod$quantiles["beta5.psi","50%"]
beta5.psi.Post.2.5<-smod$quantiles["beta5.psi","2.5%"]
beta5.psi.Post.97.5<-smod$quantiles["beta5.psi","97.5%"]

beta6.psi.PostMedian<-smod$quantiles["beta6.psi","50%"]
beta6.psi.Post.2.5<-smod$quantiles["beta6.psi","2.5%"]
beta6.psi.Post.97.5<-smod$quantiles["beta6.psi","97.5%"]

alpha.phi.PostMedian<-smod$quantiles["alpha.phi","50%"]
alpha.phi.Post.2.5<-smod$quantiles["alpha.phi","2.5%"]
alpha.phi.Post.97.5<-smod$quantiles["alpha.phi","97.5%"]

beta1.phi.PostMedian<-smod$quantiles["beta1.phi","50%"]
beta1.phi.Post.2.5<-smod$quantiles["beta1.phi","2.5%"]
beta1.phi.Post.97.5<-smod$quantiles["beta1.phi","97.5%"]

beta2.phi.PostMedian<-smod$quantiles["beta2.phi","50%"]
beta2.phi.Post.2.5<-smod$quantiles["beta2.phi","2.5%"]
beta2.phi.Post.97.5<-smod$quantiles["beta2.phi","97.5%"]

beta3.phi.PostMedian<-smod$quantiles["beta3.phi","50%"]
beta3.phi.Post.2.5<-smod$quantiles["beta3.phi","2.5%"]
beta3.phi.Post.97.5<-smod$quantiles["beta3.phi","97.5%"]

beta4.phi.PostMedian<-smod$quantiles["beta4.phi","50%"]
beta4.phi.Post.2.5<-smod$quantiles["beta4.phi","2.5%"]
beta4.phi.Post.97.5<-smod$quantiles["beta4.phi","97.5%"]

#get predictor plots
#psi covariates
#ndmi=ndmi,
#ndmi2=ndmi2,
#deltaVH=deltaVH,
#sgrass=sgrass,
#sgrass2=sgrass2,
#latsqrt=latsqrt,
#site-year covariates (gamma, phi)
#Yswi=Yswi,
#YdeltaVH=YdeltaVH,
#YdeltaVH2=YdeltaVH2,
#YsDPOLsummer=YsDPOLsummer,
#Ylatsqrt=Ylatsqrt,

#psi~alpha.psi +beta1.psi*ndmi[i] +beta2.psi*ndmi2[i] +beta3.psi*deltaVH[i] +beta4.psi*sgrass[i] +beta5.psi*sgrass2[i] +beta6.psi*latsqrt[i] + u_psi[cluster[i]]
#phi~beta1.phi*YdeltaVH[i,k-1] +beta2.phi*YdeltaVH2[i,k-1] +beta3.phi*YsDPOLsummer[i,k-1] +beta4.phi*Ylatsqrt[i,k-1] 
#gamma~beta1.gamma*Yswi[i,k-1] 
library(ggplot2)
my.theme2 <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

#phi vs delta VH  
newdat1<-data.frame(DeltaVHSY.s=YdeltaVH[,5],
                    DeltaVHSY.s2=YdeltaVH2[,5],
                    sDPOLsummer.SY.s=mean(YsDPOLsummer[,5]),
                    latsqrt=mean(Ylatsqrt))

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat1)){#instead of length(newdat1)
  int.sim[, i] <- int.mcmc.dat$alpha.phi +
    int.mcmc.dat$beta1.phi * newdat1[i,c("DeltaVHSY.s")]+
    int.mcmc.dat$beta2.phi * newdat1[i,c("DeltaVHSY.s2")]+
    int.mcmc.dat$beta3.phi * newdat1[i,c("sDPOLsummer.SY.s")]+
    int.mcmc.dat$beta4.phi * newdat1[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat1,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(DeltaVHSY.s)

#Plotting Predictions
p_phi_deltaVH<-ggplot(plot.dat, aes(x = DeltaVHSY.s, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = DeltaVHSY.s, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = DeltaVHSY.s, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("Delta VH (2021)")+ylab("Prob. Persistence")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/Phi_VS_DeltaVH_modelwithlat_notHedley.png",
       plot = p_phi_deltaVH,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#phi vs. sDPOL summer
newdat2<-data.frame(DeltaVHSY.s=mean(YdeltaVH[,5]),
                    DeltaVHSY.s2=mean(YdeltaVH2[,5]),
                    sDPOLsummer.SY.s=YsDPOLsummer[,5],
                    latsqrt=mean(Ylatsqrt))

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat2)){#instead of length(newdat2)
  int.sim[, i] <- int.mcmc.dat$alpha.phi +
    int.mcmc.dat$beta1.phi * newdat2[i,c("DeltaVHSY.s")]+
    int.mcmc.dat$beta2.phi * newdat2[i,c("DeltaVHSY.s2")]+
    int.mcmc.dat$beta3.phi * newdat2[i,c("sDPOLsummer.SY.s")]+
    int.mcmc.dat$beta4.phi * newdat2[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat2,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(sDPOLsummer.SY.s)

#Plotting Predictions
p_phi_dpol<-ggplot(plot.dat, aes(x = sDPOLsummer.SY.s, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = sDPOLsummer.SY.s, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = sDPOLsummer.SY.s, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("Mean Summer DPOL (2016-2021)")+ylab("Prob. Persistence")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/Phi_VS_meansummerDPOL_modelwithlat_notHedley.png",
       plot = p_phi_dpol,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#phi vs. lat
newdat3<-data.frame(DeltaVHSY.s=mean(YdeltaVH[,5]),
                    DeltaVHSY.s2=mean(YdeltaVH2[,5]),
                    sDPOLsummer.SY.s=mean(YsDPOLsummer[,5]),
                    latsqrt=Ylatsqrt[,5],
                    lat=umf@siteCovs$lat)

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat3)){#instead of length(newdat3)
  int.sim[, i] <- int.mcmc.dat$alpha.phi +
    int.mcmc.dat$beta1.phi * newdat3[i,c("DeltaVHSY.s")]+
    int.mcmc.dat$beta2.phi * newdat3[i,c("DeltaVHSY.s2")]+
    int.mcmc.dat$beta3.phi * newdat3[i,c("sDPOLsummer.SY.s")]+
    int.mcmc.dat$beta4.phi * newdat3[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat3,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(lat)

#Plotting Predictions
p_phi_lat<-ggplot(plot.dat, aes(x = lat, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = lat, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = lat, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("Latitude")+ylab("Prob. Persistence")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/Phi_VS_sqrtlat_modelwithlat_notHedley.png",
       plot = p_phi_lat,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#gamma: SWI
newdat4<-data.frame(SWI.s=Yswi[,5])

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat4)){#instead of length(newdat2)
  int.sim[, i] <- int.mcmc.dat$alpha.gamma +
    int.mcmc.dat$beta1.gamma * newdat4[i,c("SWI.s")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat4,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(SWI.s)

#Plotting Predictions
p_gamma_SWI<-ggplot(plot.dat, aes(x = SWI.s, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = SWI.s, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = SWI.s, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("SAGA Wetness Index")+ylab("Prob. Colonization")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/Gamma_VS_SWI_modelwithlat_notHedley.png",
       plot = p_gamma_SWI,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#psi vs deltaVH  
newdat5<-data.frame(ndmi=ndmi,
                    ndmi2=ndmi2,
                    deltaVH=mean(deltaVH),
                    sgrass=mean(sgrass),
                    sgrass2=mean(sgrass2),
                    latsqrt=mean(latsqrt),
                    lat=mean(umf@siteCovs$lat))

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat5)){#instead of length(newdat5)
  int.sim[, i] <- int.mcmc.dat$alpha.psi +
    int.mcmc.dat$beta1.psi * newdat5[i,c("ndmi")]+
    int.mcmc.dat$beta2.psi * newdat5[i,c("ndmi2")]+
    int.mcmc.dat$beta3.psi * newdat5[i,c("deltaVH")]+
    int.mcmc.dat$beta4.psi * newdat5[i,c("sgrass")]+
    int.mcmc.dat$beta5.psi * newdat5[i,c("sgrass2")]+
    int.mcmc.dat$beta6.psi * newdat5[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat5,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(ndmi)

#Plotting Predictions
p_psi_ndmi<-ggplot(plot.dat, aes(x = ndmi, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = ndmi, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = ndmi, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("NDMI (2021)")+ylab("Prob. Init. Occupancy")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psi_VS_NDMI_modelwithlat_notHedley.png",
       plot = p_psi_ndmi,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#psi vs Delta VH 
newdat6<-data.frame(ndmi=mean(ndmi),
                    ndmi2=mean(ndmi2),
                    deltaVH=deltaVH,
                    sgrass=mean(sgrass),
                    sgrass2=mean(sgrass2),
                    latsqrt=mean(latsqrt),
                    lat=mean(umf@siteCovs$lat))

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat6)){#instead of length(newdat6)
  int.sim[, i] <- int.mcmc.dat$alpha.psi +
    int.mcmc.dat$beta1.psi * newdat6[i,c("ndmi")]+
    int.mcmc.dat$beta2.psi * newdat6[i,c("ndmi2")]+
    int.mcmc.dat$beta3.psi * newdat6[i,c("deltaVH")]+
    int.mcmc.dat$beta4.psi * newdat6[i,c("sgrass")]+
    int.mcmc.dat$beta5.psi * newdat6[i,c("sgrass2")]+
    int.mcmc.dat$beta6.psi * newdat6[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat6,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(deltaVH)

#Plotting Predictions
p_psi_deltaVH<-ggplot(plot.dat, aes(x = deltaVH, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = deltaVH, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = deltaVH, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("Delta VH (2021)")+ylab("Prob. Init. Occupancy")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psi_VS_deltaVH_modelwithlat_notHedley.png",
       plot = p_psi_deltaVH,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#psi vs grass cover  
newdat7<-data.frame(ndmi=mean(ndmi),
                    ndmi2=mean(ndmi2),
                    deltaVH=mean(deltaVH),
                    sgrass=sgrass,
                    sgrass2=sgrass2,
                    latsqrt=mean(latsqrt),
                    lat=mean(umf@siteCovs$lat))

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat7)){#instead of length(newdat7)
  int.sim[, i] <- int.mcmc.dat$alpha.psi +
    int.mcmc.dat$beta1.psi * newdat7[i,c("ndmi")]+
    int.mcmc.dat$beta2.psi * newdat7[i,c("ndmi2")]+
    int.mcmc.dat$beta3.psi * newdat7[i,c("deltaVH")]+
    int.mcmc.dat$beta4.psi * newdat7[i,c("sgrass")]+
    int.mcmc.dat$beta5.psi * newdat7[i,c("sgrass2")]+
    int.mcmc.dat$beta6.psi * newdat7[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat7,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(sgrass)

#Plotting Predictions
p_psi_grass<-ggplot(plot.dat, aes(x = sgrass, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = sgrass, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = sgrass, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("Grass cover (2016-2021)")+ylab("Prob. Init. Occupancy")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psi_VS_grasscover_modelwithlat_notHedley.png",
       plot = p_psi_grass,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#psi vs latitude 
newdat8<-data.frame(ndmi=mean(ndmi),
                    ndmi2=mean(ndmi2),
                    deltaVH=mean(deltaVH),
                    sgrass=mean(sgrass),
                    sgrass2=mean(sgrass2),
                    latsqrt=latsqrt,
                    lat=umf@siteCovs$lat)

int.sim <- matrix(rep(NA,
                      nrow(int.mcmc.dat)),
                  nrow = nrow(int.mcmc.dat), ncol = nrow(int.mcmc.dat))
for(i in 1:nrow(newdat8)){#instead of length(newdat8)
  int.sim[, i] <- int.mcmc.dat$alpha.psi +
    int.mcmc.dat$beta1.psi * newdat8[i,c("ndmi")]+
    int.mcmc.dat$beta2.psi * newdat8[i,c("ndmi2")]+
    int.mcmc.dat$beta3.psi * newdat8[i,c("deltaVH")]+
    int.mcmc.dat$beta4.psi * newdat8[i,c("sgrass")]+
    int.mcmc.dat$beta5.psi * newdat8[i,c("sgrass2")]+
    int.mcmc.dat$beta6.psi * newdat8[i,c("latsqrt")]
}
int.sim <- plogis(int.sim)
#Note: the variance now comes from the posterior, not the vcov matrix.
bayes.c.eff.median <- apply(int.sim, 2, median)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025), na.rm=TRUE))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975), na.rm=TRUE))

#Combine estimates into one dataframe for plotting:
plot.dat <- data.frame(newdat8,
                       bayes.c.eff.median[1:115],
                       bayes.c.eff.lower[1:115],
                       bayes.c.eff.upper[1:115])
plot.dat<-plot.dat%>%
  arrange(lat)

#Plotting Predictions
p_psi_lat<-ggplot(plot.dat, aes(x = lat, y = bayes.c.eff.median.1.115.)) +
  geom_point() +
  geom_line(data = plot.dat, aes(x = lat, y = bayes.c.eff.median.1.115.), color = "blue") +
  geom_ribbon(data = plot.dat, aes(x = lat, y = bayes.c.eff.median.1.115., ymin = bayes.c.eff.lower.1.115., ymax = bayes.c.eff.upper.1.115.), alpha = 0.2, fill = "blue") +
  xlab("Latitude")+ylab("Prob. Init. Occupancy")+my.theme2

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psi_VS_latitude_modelwithlat_notHedley.png",
       plot = p_psi_lat,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

# Define a layout matrix with NA for blank space
# Layout: p_psi_deltaVH, p_psi_ndmi, p_psi_grass, p_psi_lat,
#         p_gamma_SWI, blank, blank, blank,
#         p_phi_deltaVH, p_phi_dpol, p_phi_lat, blank
layout <- matrix(c(1, 2, 3, 4, 5, NA, NA, NA, 6, 7, 8, NA), nrow = 3, byrow = TRUE)

# Arrange plots using the layout matrix
g8<-grid.arrange(p_psi_deltaVH, p_psi_ndmi, p_psi_grass, p_psi_lat, 
             p_gamma_SWI, p_phi_deltaVH, p_phi_dpol, p_phi_lat,
             layout_matrix = layout)

#g8<-grid.arrange(p_psi_deltaVH,
#                 p_psi_ndmi,
#                 p_psi_grass,
#                 p_psi_lat,
#                 p_gamma_SWI,
#                 p_phi_deltaVH,
#                 p_phi_dpol,
#                 p_phi_lat,
#                 nrow=3,ncol=4)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/psi_gamma_phi_predplots_alt.png",
       plot = g8,
       width = 14,
       height = 10,
       units = c("in"),
       dpi = 300)

#Predicted yearly occu. for individual clusters
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
#575 entries = 115 sites * 5 years
library(ggmcmc)
allZs<-posterior_draws[1:3][,4407:4981]#
ggd <- ggs(allZs)

nsite<-115
niter <-60000#instead of 40000
nthin <-10
nchain <-3

nstore <- niter / nthin

#recreate obsFiltFinal dataset but with location and clusterID
#instead of YERA detections
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
#from the whole YERA study area. 2016 is first year
#for which we have Sentinel data.
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


m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits
str(m2)
table(m2$clusterID)#136 sites
# 2_1 2_10  2_8  3_1  3_6  3_8  4_5  4_8  7_5  7_8  9_1
#   3    5   15   13    7   11   45    2   33   21    5

table(m2$clusterIDfinal)
#2_1 2_10 2_8b  3_1  3_6  3_8  4_5  4_8  7_5 
#  3    5   15   13    7   11   45    2   33 
#7_8 7_8b 7_8c 7_8d  9_1 
# 13    1    3    4    5  

#get covariates affecting initial occu.
initFilt<-m2%>%
  #select(locationC, clusterID)%>%
  arrange(locationC)

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

#Value:  2  3  4  5  6  7  8 
#nsite:  3  7  3 22 35 36 40

initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11

obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID)%>%
  mutate(NScluster=case_when(clusterID %in% c("3_6","4_5","7_5") ~ "N",
                             clusterID %in% c("3_1","7_8","3_8","2_8") ~ "S"))
obsFiltFinal$geometry<-NULL

#Probabilities of occu. at year 1 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 1, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2016<-rowMeans(z.pred)
OccEst2016<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2016$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2016bycluster<-OccEst2016%>%
  select(-NScluster)%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2016bycluster_quantiles <- OccEst2016bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2016)

OccEst2016byNS<-OccEst2016%>%
  select(-clusterID)%>%
  group_by(NScluster)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2016byNS_quantiles <- OccEst2016byNS %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(NScluster, q2.5,q50,q97.5)%>%
  mutate(year=2016)
OccEst2016byNS_quantiles$geometry<-NULL

#Probabilities of occu. at year 2 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2017<-rowMeans(z.pred)
OccEst2017<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2017$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2017bycluster<-OccEst2017%>%
  select(-NScluster)%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2017bycluster_quantiles <- OccEst2017bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2017)


OccEst2017byNS<-OccEst2017%>%
  select(-clusterID)%>%
  group_by(NScluster)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2017byNS_quantiles <- OccEst2017byNS %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(NScluster, q2.5,q50,q97.5)%>%
  mutate(year=2017)
OccEst2017byNS_quantiles$geometry<-NULL

#Probabilities of occu. at year 3 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2018<-rowMeans(z.pred)
OccEst2018<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2018$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2018bycluster<-OccEst2018%>%
  select(-NScluster)%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2018bycluster_quantiles <- OccEst2018bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2018)


OccEst2018byNS<-OccEst2018%>%
  select(-clusterID)%>%
  group_by(NScluster)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2018byNS_quantiles <- OccEst2018byNS %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(NScluster, q2.5,q50,q97.5)%>%
  mutate(year=2018)
OccEst2018byNS_quantiles$geometry<-NULL


#Probabilities of occu. at year 4 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2019<-rowMeans(z.pred)
OccEst2019<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2019$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2019bycluster<-OccEst2019%>%
  select(-NScluster)%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2019bycluster_quantiles <- OccEst2019bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2019)

OccEst2019byNS<-OccEst2019%>%
  select(-clusterID)%>%
  group_by(NScluster)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2019byNS_quantiles <- OccEst2019byNS %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(NScluster, q2.5,q50,q97.5)%>%
  mutate(year=2019)
OccEst2019byNS_quantiles$geometry<-NULL


#Probabilities of occu. at year 5 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 5, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2021<-rowMeans(z.pred)
OccEst2021<-cbind(obsFiltFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2021$locationC<-NULL 
#drop before summarizing everything by clusterID
OccEst2021bycluster<-OccEst2021%>%
  select(-NScluster)%>%
  group_by(clusterID)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2021bycluster_quantiles <- OccEst2021bycluster %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`X1_mean`:`X18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(clusterID, q2.5,q50,q97.5)%>%
  mutate(year=2021)


OccEst2021byNS<-OccEst2021%>%
  select(-clusterID)%>%
  group_by(NScluster)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2021byNS_quantiles <- OccEst2021byNS %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(NScluster, q2.5,q50,q97.5)%>%
  mutate(year=2021)
OccEst2021byNS_quantiles$geometry<-NULL


#meanOccEstAllYears<-cbind(meanOccEst2016,
#                          meanOccEst2017,
#                          meanOccEst2018,
#                          meanOccEst2019,
#                          meanOccEst2021)

finalOccDF<-rbind(OccEst2016bycluster_quantiles,
                  OccEst2017bycluster_quantiles,
                  OccEst2018bycluster_quantiles,
                  OccEst2019bycluster_quantiles,
                  OccEst2021bycluster_quantiles)

finalOccDF_byNS<-rbind(OccEst2016byNS_quantiles,
                  OccEst2017byNS_quantiles,
                  OccEst2018byNS_quantiles,
                  OccEst2019byNS_quantiles,
                  OccEst2021byNS_quantiles)


#finalOccDF<-cbind(obsFiltFinal, meanOccEstAllYears)
#finalOccDF$clusterID
#write.csv(finalOccDF,file="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/HedleyHybridExclLatitude_2016_2021_PredictedOcc.csv")


#Cluster 1 (formerly 3_6)
finalOccDF_cluster1<-finalOccDF%>%
  filter(clusterID=="3_6")

# finalOccDF_cluster1long<-finalOccDF_cluster1%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster1)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_6 <- ggplot(data=finalOccDF_cluster1long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 1: ",nrow(finalOccDF_cluster1)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_6 <- ggplot(data=finalOccDF_cluster1, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 1: ",nrow(obsFiltFinal%>%filter(clusterID=="3_6"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 2
finalOccDF_cluster2<-finalOccDF%>%
  filter(clusterID=="4_5")

# finalOccDF_cluster2long<-finalOccDF_cluster2%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster2)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster4_5 <- ggplot(data=finalOccDF_cluster2long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 2: ",nrow(finalOccDF_cluster2)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster4_5 <- ggplot(data=finalOccDF_cluster2, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 2: ",nrow(obsFiltFinal%>%filter(clusterID=="4_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 3
finalOccDF_cluster3<-finalOccDF%>%
  filter(clusterID=="7_5")

# finalOccDF_cluster3long<-finalOccDF_cluster3%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster3)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

#Occucluster7_5 <- ggplot(data=finalOccDF_cluster3long, aes(x=Year, y=medYear)) +
#  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#  geom_line(aes(x=Year, y=medYear), colour="blue")+
#  my.theme+ggtitle(paste0("Cluster 3: ",nrow(finalOccDF_cluster3)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#  ylim(-0.5,1.25)+xlim(2016,2021)+
#  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster7_5 <- ggplot(data=finalOccDF_cluster3, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 3: ",nrow(obsFiltFinal%>%filter(clusterID=="7_5"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 4 (just 3 sites at GG Provincial Park)
#No sites with 3 or more years of obs from 2016-2021 in cluster 4
# finalOccDF_cluster4<-finalOccDF%>%
#     filter(clusterID=="2_1")
#   
# finalOccDF_cluster4long<-finalOccDF_cluster4%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster4)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster2_1 <- ggplot(data=finalOccDF_cluster4long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 4: ",nrow(finalOccDF_cluster4)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 5
finalOccDF_cluster5<-finalOccDF%>%
  filter(clusterID=="3_1")

# finalOccDF_cluster5long<-finalOccDF_cluster5%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster5)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_1 <- ggplot(data=finalOccDF_cluster5long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 5: ",nrow(finalOccDF_cluster5)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_1 <- ggplot(data=finalOccDF_cluster5, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 5: ",nrow(obsFiltFinal%>%filter(clusterID=="3_1"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 6
finalOccDF_cluster6<-finalOccDF%>%
  filter(clusterID=="7_8")

# finalOccDF_cluster6long<-finalOccDF_cluster6%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster6)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster7_8 <- ggplot(data=finalOccDF_cluster6long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 6: ",nrow(finalOccDF_cluster6)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


Occucluster7_8 <- ggplot(data=finalOccDF_cluster6, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 6: ",nrow(obsFiltFinal%>%filter(clusterID=="7_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 7
finalOccDF_cluster7<-finalOccDF%>%
  filter(clusterID=="3_8")

# finalOccDF_cluster7long<-finalOccDF_cluster7%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster7)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster3_8 <- ggplot(data=finalOccDF_cluster7long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 7: ",nrow(finalOccDF_cluster7)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster3_8 <- ggplot(data=finalOccDF_cluster7, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 7: ",nrow(obsFiltFinal%>%filter(clusterID=="3_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Cluster 8
finalOccDF_cluster8<-finalOccDF%>%
  filter(clusterID=="2_8")

# finalOccDF_cluster8long<-finalOccDF_cluster8%>%
#   pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                names_to = "year",
#                values_to = "occ")%>%
#   mutate(Year=substr(year,11,14))%>%
#   group_by(Year)%>%
#   summarize(medYear=median(occ),
#             sdYear=sd(occ),
#             seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster8)),
#             bci2.5=medYear-1.96*seYear,
#             bci97.5=medYear+1.96*seYear)%>%
#   mutate(Year=as.integer(Year))

# Occucluster2_8 <- ggplot(data=finalOccDF_cluster8long, aes(x=Year, y=medYear)) +
#   geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#   geom_line(aes(x=Year, y=medYear), colour="blue")+
#   my.theme+ggtitle(paste0("Cluster 8: ",nrow(finalOccDF_cluster8)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#   ylim(-0.5,1.25)+xlim(2016,2021)+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster2_8 <- ggplot(data=finalOccDF_cluster8, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 8: ",nrow(obsFiltFinal%>%filter(clusterID=="2_8"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))


#Cluster 9 No sites with at least 3 years of data from 2016-2021 in this cluster
# finalOccDF_cluster9<-finalOccDF%>%
#     filter(clusterID=="4_8")
#   
# finalOccDF_cluster9long<-finalOccDF_cluster9%>%
#     pivot_longer(cols=meanOccEst2016:meanOccEst2021,
#                  names_to = "year",
#                  values_to = "occ")%>%
#     mutate(Year=substr(year,11,14))%>%
#     group_by(Year)%>%
#     summarize(medYear=median(occ),
#               sdYear=sd(occ),
#               seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster9)),
#               bci2.5=medYear-1.96*seYear,
#               bci97.5=medYear+1.96*seYear)%>%
#     mutate(Year=as.integer(Year))
#   
# Occucluster4_8 <- ggplot(data=finalOccDF_cluster9long, aes(x=Year, y=medYear)) +
#     geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
#     geom_line(aes(x=Year, y=medYear), colour="blue")+
#     my.theme+ggtitle(paste0("Cluster 9: ",nrow(finalOccDF_cluster9)," sites"))+xlab("Year")+ylab("Pred. occu.")+
#     ylim(-0.5,1.25)+xlim(2016,2021)+
#     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

p9<-grid.arrange(Occucluster3_6,
                 Occucluster4_5,
                 Occucluster7_5,
                 #Occucluster2_1,
                 Occucluster3_1,
                 Occucluster7_8,
                 Occucluster3_8,
                 Occucluster2_8,
                 #Occucluster4_8,
                 nrow=3,ncol=3)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_MedianPlus95BCIoccupancy_2016_2023Dec8.png",
       plot = p9,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Compare Occupancy at Northern (1-3) versus Southern Sites (5-8)
#Northern Clusters 3_6, 4_5, 7_5
finalOccDF_north<-finalOccDF_byNS%>%
  filter(NScluster%in%c("N"))

Occucluster_North <- ggplot(data=finalOccDF_north, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme2+ggtitle(paste0("North clusters: ",nrow(obsFiltFinal%>%filter(NScluster=="N"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(0,1)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Southern Clusters 3_1, 3_8, 2_8, 7_8
finalOccDF_south<-finalOccDF_byNS%>%
  filter(NScluster%in%c("S"))

Occucluster_South <- ggplot(data=finalOccDF_south, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme2+ggtitle(paste0("South clusters: ",nrow(obsFiltFinal%>%filter(NScluster=="S"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(0,1)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

Occucluster_NvsS<-grid.arrange(Occucluster_North,
                 Occucluster_South,
                 nrow=1,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_MedianPlus95BCIoccupancy_2016_2023_NvsS.png",
       plot = Occucluster_NvsS,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Get Estimated Number of Occupied Sites (out of 115)
which(substr(colnames(posterior_draws[[1]]),1,6)=="n.occ[")
# [1] 8317 8318 8319 8320 8321
colnames(posterior_draws[[1]])[8317]#n.occ[1]

n.occ1.PostMedian<-smod$quantiles["n.occ[1]","50%"]
n.occ1.Post.2.5<-smod$quantiles["n.occ[1]","2.5%"]
n.occ1.Post.97.5<-smod$quantiles["n.occ[1]","97.5%"]

n.occ2.PostMedian<-smod$quantiles["n.occ[2]","50%"]
n.occ2.Post.2.5<-smod$quantiles["n.occ[2]","2.5%"]
n.occ2.Post.97.5<-smod$quantiles["n.occ[2]","97.5%"]

n.occ3.PostMedian<-smod$quantiles["n.occ[3]","50%"]
n.occ3.Post.2.5<-smod$quantiles["n.occ[3]","2.5%"]
n.occ3.Post.97.5<-smod$quantiles["n.occ[3]","97.5%"]

n.occ4.PostMedian<-smod$quantiles["n.occ[4]","50%"]
n.occ4.Post.2.5<-smod$quantiles["n.occ[4]","2.5%"]
n.occ4.Post.97.5<-smod$quantiles["n.occ[4]","97.5%"]

n.occ5.PostMedian<-smod$quantiles["n.occ[5]","50%"]
n.occ5.Post.2.5<-smod$quantiles["n.occ[5]","2.5%"]
n.occ5.Post.97.5<-smod$quantiles["n.occ[5]","97.5%"]

#Get Mean Yearly Colonization Probability Across Sites (out of 115)
which(substr(colnames(posterior_draws[[1]]),1,6)=="gamma[")
# [1] 1062..1521 (460 estimates = 115 sites * 4 years)
colnames(posterior_draws[[1]])[1062]#gamma[1,1]
colnames(posterior_draws[[1]])[1063]#gamma[2,1]
colnames(posterior_draws[[1]])[1521]#gamma[115,4]
#ordered by time period then by site within time period

gamma1.1.PostMedian<-smod$quantiles["gamma[1,1]","50%"]
gamma1.1.Post.2.5<-smod$quantiles["gamma[1,1]","2.5%"]
gamma1.1.Post.97.5<-smod$quantiles["gamma[1,1]","97.5%"]

gammaList_year1<-list()
for (i in 1:115) {
  gamma_y1_PostMedian<-smod$quantiles[paste0("gamma[",i,",1]"),"50%"]
  gamma_y1_Post.2.5<-smod$quantiles[paste0("gamma[",i,",1]"),"2.5%"]
  gamma_y1_Post.97.5<-smod$quantiles[paste0("gamma[",i,",1]"),"97.5%"]
  gammaList_year1[[i]]<-data.frame(gamma_y1_PostMedian=gamma_y1_PostMedian,
                                   gamma_y1_Post.2.5=gamma_y1_Post.2.5,
                                   gamma_y1_Post.97.5=gamma_y1_Post.97.5)
}
gamma_year1<-do.call(rbind, gammaList_year1)
gamma_year1median<-gamma_year1%>%summarize(gamma_y1_PostMedian=median(gamma_y1_PostMedian),
                        gamma_y1_Post.2.5=median(gamma_y1_Post.2.5),
                        gamma_y1_Post.97.5=median(gamma_y1_Post.97.5))

gammaList_year2<-list()
for (i in 1:115) {
  gamma_y2_PostMedian<-smod$quantiles[paste0("gamma[",i,",2]"),"50%"]
  gamma_y2_Post.2.5<-smod$quantiles[paste0("gamma[",i,",2]"),"2.5%"]
  gamma_y2_Post.97.5<-smod$quantiles[paste0("gamma[",i,",2]"),"97.5%"]
  gammaList_year2[[i]]<-data.frame(gamma_y2_PostMedian=gamma_y2_PostMedian,
                                   gamma_y2_Post.2.5=gamma_y2_Post.2.5,
                                   gamma_y2_Post.97.5=gamma_y2_Post.97.5)
}
gamma_year2<-do.call(rbind, gammaList_year2)
gamma_year2median<-gamma_year2%>%summarize(gamma_y2_PostMedian=median(gamma_y2_PostMedian),
                                           gamma_y2_Post.2.5=median(gamma_y2_Post.2.5),
                                           gamma_y2_Post.97.5=median(gamma_y2_Post.97.5))

gammaList_year3<-list()
for (i in 1:115) {
  gamma_y3_PostMedian<-smod$quantiles[paste0("gamma[",i,",3]"),"50%"]
  gamma_y3_Post.2.5<-smod$quantiles[paste0("gamma[",i,",3]"),"2.5%"]
  gamma_y3_Post.97.5<-smod$quantiles[paste0("gamma[",i,",3]"),"97.5%"]
  gammaList_year3[[i]]<-data.frame(gamma_y3_PostMedian=gamma_y3_PostMedian,
                                   gamma_y3_Post.2.5=gamma_y3_Post.2.5,
                                   gamma_y3_Post.97.5=gamma_y3_Post.97.5)
}
gamma_year3<-do.call(rbind, gammaList_year3)
gamma_year3median<-gamma_year3%>%summarize(gamma_y3_PostMedian=median(gamma_y3_PostMedian),
                                           gamma_y3_Post.2.5=median(gamma_y3_Post.2.5),
                                           gamma_y3_Post.97.5=median(gamma_y3_Post.97.5))


gammaList_year4<-list()
for (i in 1:115) {
  gamma_y4_PostMedian<-smod$quantiles[paste0("gamma[",i,",4]"),"50%"]
  gamma_y4_Post.2.5<-smod$quantiles[paste0("gamma[",i,",4]"),"2.5%"]
  gamma_y4_Post.97.5<-smod$quantiles[paste0("gamma[",i,",4]"),"97.5%"]
  gammaList_year4[[i]]<-data.frame(gamma_y4_PostMedian=gamma_y4_PostMedian,
                                   gamma_y4_Post.2.5=gamma_y4_Post.2.5,
                                   gamma_y4_Post.97.5=gamma_y4_Post.97.5)
}
gamma_year4<-do.call(rbind, gammaList_year4)
gamma_year4median<-gamma_year4%>%summarize(gamma_y4_PostMedian=median(gamma_y4_PostMedian),
                                           gamma_y4_Post.2.5=median(gamma_y4_Post.2.5),
                                           gamma_y4_Post.97.5=median(gamma_y4_Post.97.5))

#gamma_y1_PostMedian gamma_y1_Post.2.5 gamma_y1_Post.97.5
#          0.2706265         0.1625311          0.4260332
#gamma_y2_PostMedian gamma_y2_Post.2.5 gamma_y2_Post.97.5
#          0.2706265         0.1625311          0.4260332
#gamma_y3_PostMedian gamma_y3_Post.2.5 gamma_y3_Post.97.5
#          0.2706265         0.1625311          0.4260332
#gamma_y4_PostMedian gamma_y4_Post.2.5 gamma_y4_Post.97.5
#          0.2706265         0.1625311          0.4260332

#Get Mean Yearly Persistence Probability Across Sites (out of 115)
which(substr(colnames(posterior_draws[[1]]),1,4)=="phi[")
# [1] 590..1049 (460 estimates = 115 sites * 4 years)
colnames(posterior_draws[[1]])[590]#phi[1,1]
colnames(posterior_draws[[1]])[591]#phi[2,1]
colnames(posterior_draws[[1]])[1049]#phi[115,4]
#ordered by time period then by site within time period

phi1.1.PostMedian<-smod$quantiles["phi[1,1]","50%"]
phi1.1.Post.2.5<-smod$quantiles["phi[1,1]","2.5%"]
phi1.1.Post.97.5<-smod$quantiles["phi[1,1]","97.5%"]

phiList_year1<-list()
for (i in 1:115) {
  phi_y1_PostMedian<-smod$quantiles[paste0("phi[",i,",1]"),"50%"]
  phi_y1_Post.2.5<-smod$quantiles[paste0("phi[",i,",1]"),"2.5%"]
  phi_y1_Post.97.5<-smod$quantiles[paste0("phi[",i,",1]"),"97.5%"]
  phiList_year1[[i]]<-data.frame(phi_y1_PostMedian=phi_y1_PostMedian,
                                   phi_y1_Post.2.5=phi_y1_Post.2.5,
                                   phi_y1_Post.97.5=phi_y1_Post.97.5)
}
phi_year1<-do.call(rbind, phiList_year1)
phi_year1median<-phi_year1%>%summarize(phi_y1_PostMedian=median(phi_y1_PostMedian),
                                           phi_y1_Post.2.5=median(phi_y1_Post.2.5),
                                           phi_y1_Post.97.5=median(phi_y1_Post.97.5))

phiList_year2<-list()
for (i in 1:115) {
  phi_y2_PostMedian<-smod$quantiles[paste0("phi[",i,",2]"),"50%"]
  phi_y2_Post.2.5<-smod$quantiles[paste0("phi[",i,",2]"),"2.5%"]
  phi_y2_Post.97.5<-smod$quantiles[paste0("phi[",i,",2]"),"97.5%"]
  phiList_year2[[i]]<-data.frame(phi_y2_PostMedian=phi_y2_PostMedian,
                                   phi_y2_Post.2.5=phi_y2_Post.2.5,
                                   phi_y2_Post.97.5=phi_y2_Post.97.5)
}
phi_year2<-do.call(rbind, phiList_year2)
phi_year2median<-phi_year2%>%summarize(phi_y2_PostMedian=median(phi_y2_PostMedian),
                                           phi_y2_Post.2.5=median(phi_y2_Post.2.5),
                                           phi_y2_Post.97.5=median(phi_y2_Post.97.5))

phiList_year3<-list()
for (i in 1:115) {
  phi_y3_PostMedian<-smod$quantiles[paste0("phi[",i,",3]"),"50%"]
  phi_y3_Post.2.5<-smod$quantiles[paste0("phi[",i,",3]"),"2.5%"]
  phi_y3_Post.97.5<-smod$quantiles[paste0("phi[",i,",3]"),"97.5%"]
  phiList_year3[[i]]<-data.frame(phi_y3_PostMedian=phi_y3_PostMedian,
                                   phi_y3_Post.2.5=phi_y3_Post.2.5,
                                   phi_y3_Post.97.5=phi_y3_Post.97.5)
}
phi_year3<-do.call(rbind, phiList_year3)
phi_year3median<-phi_year3%>%summarize(phi_y3_PostMedian=median(phi_y3_PostMedian),
                                           phi_y3_Post.2.5=median(phi_y3_Post.2.5),
                                           phi_y3_Post.97.5=median(phi_y3_Post.97.5))


phiList_year4<-list()
for (i in 1:115) {
  phi_y4_PostMedian<-smod$quantiles[paste0("phi[",i,",4]"),"50%"]
  phi_y4_Post.2.5<-smod$quantiles[paste0("phi[",i,",4]"),"2.5%"]
  phi_y4_Post.97.5<-smod$quantiles[paste0("phi[",i,",4]"),"97.5%"]
  phiList_year4[[i]]<-data.frame(phi_y4_PostMedian=phi_y4_PostMedian,
                                   phi_y4_Post.2.5=phi_y4_Post.2.5,
                                   phi_y4_Post.97.5=phi_y4_Post.97.5)
}
phi_year4<-do.call(rbind, phiList_year4)
phi_year4median<-phi_year4%>%summarize(phi_y4_PostMedian=median(phi_y4_PostMedian),
                                           phi_y4_Post.2.5=median(phi_y4_Post.2.5),
                                           phi_y4_Post.97.5=median(phi_y4_Post.97.5))

#phi_y1_PostMedian phi_y1_Post.2.5 phi_y1_Post.97.5
#        0.4170265       0.2106755        0.6780975
#phi_y2_PostMedian phi_y2_Post.2.5 phi_y2_Post.97.5
#        0.3743985       0.1722623        0.6471765
#phi_y3_PostMedian phi_y3_Post.2.5 phi_y3_Post.97.5
#        0.3659255        0.177599        0.6407346
#phi_y4_PostMedian phi_y4_Post.2.5 phi_y4_Post.97.5
#        0.467077       0.2457541        0.7195799

n.occDF<-data.frame(year=c(2016,2017,2018,2019,2021),
                    nsitesVisited=c(82,101,62,84,79),
                    nsitesDetected=c(47,58,12,30,25),
                    phiMedian=c(NA,0.42,0.37,0.37,0.48),
                    phi.2.5=c(NA,0.21,0.17,0.18,0.25),
                    phi.97.5=c(NA,0.68,0.65,0.64,0.72),
                    gammaMedian=c(NA,0.27,0.27,0.27,0.27),
                    gamma.2.5=c(NA,0.16,0.16,0.16,0.16),
                    gamma.97.5=c(NA,0.43,0.43,0.43,0.43),
                    n.occMedian=c(n.occ1.PostMedian,
                                  n.occ2.PostMedian,
                                  n.occ3.PostMedian,
                                  n.occ4.PostMedian,
                                  n.occ5.PostMedian),
                    n.occ.2.5=c(n.occ1.Post.2.5,
                                n.occ2.Post.2.5,
                                n.occ3.Post.2.5,
                                n.occ4.Post.2.5,
                                n.occ5.Post.2.5),
                    n.occ.97.5=c(n.occ1.Post.97.5,
                                 n.occ2.Post.97.5,
                                 n.occ3.Post.97.5,
                                 n.occ4.Post.97.5,
                                 n.occ5.Post.97.5))

gammaphiDF<-data.frame(year=c(2017,2018,2019,2021),
                    phiMedian=c(0.42,0.37,0.37,0.48),
                    phi.2.5=c(0.21,0.17,0.18,0.25),
                    phi.97.5=c(0.68,0.65,0.64,0.72),
                    gammaMedian=c(0.27,0.27,0.27,0.27),
                    gamma.2.5=c(0.16,0.16,0.16,0.16),
                    gamma.97.5=c(0.43,0.43,0.43,0.43))



cor(n.occDF$nsitesVisited, n.occDF$n.occMedian)
#[1] 0.7058976
n.occDF$propDetected<-n.occDF$nsitesDetected/n.occDF$nsitesVisited

p1<-ggplot(n.occDF, aes(x=year, y=n.occMedian))+
  #geom_point(aes(x=year, y=nsitesVisited), col="blue")+
  #geom_line(aes(x=year, y=nsitesVisited), col="blue")+
  #geom_point(aes(x=year, y=n.occMedian), col="red")+
  #geom_line(aes(x=year, y=n.occMedian), col="red")+
  #geom_line(aes(x=year, y=n.occ.2.5), linetype="dashed")+
  #geom_line(aes(x=year, y=n.occ.97.5), linetype="dashed")+
  geom_ribbon(aes(x=year, ymin=n.occ.2.5, ymax=n.occ.97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=n.occMedian), colour="blue")+
  xlab("Year")+ylab("Number of Occupied Sites")+
  xlim(2016,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))+
  my.theme2


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_2016_2021.png",
       plot = p1,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p2<-ggplot(n.occDF, aes(x=year, y=propDetected))+
  geom_point(aes(x=year, y=propDetected), col="blue")+
  geom_line(aes(x=year, y=propDetected), col="blue")+
  xlab("Year")+ylab("Proportion of Visited Sites With Detections")+
  xlim(2016,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))

ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_2016_2021_SitesWithDetections.png",
       plot = p2,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p3<-grid.arrange(p1,p2,nrow=1,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_NumberOccupiedSitesPlusRaw_2016_2023.png",
       plot = p3,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p4<-ggplot(gammaphiDF, aes(x=year, y=phiMedian))+
  geom_ribbon(aes(x=year, ymin=phi.2.5, ymax=phi.97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=phiMedian), colour="blue")+
  xlab("Year")+ylab("Prob. Persistence")+
  xlim(2017,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2017,2021)))))+
  my.theme2


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_2016_2021_yearlyPersistence.png",
       plot = p4,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


p5<-ggplot(gammaphiDF, aes(x=year, y=gammaMedian))+
  geom_ribbon(aes(x=year, ymin=gamma.2.5, ymax=gamma.97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=gammaMedian), colour="blue")+
  xlab("Year")+ylab("Prob. Colonization")+
  xlim(2017,2021)+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2017,2021)))))+
  my.theme2


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_2016_2021_yearlyColonization.png",
       plot = p5,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Get Turnover
which(substr(colnames(posterior_draws[[1]]),1,9)=="turnover[")
#[1]8326 8327 8328 8329
colnames(posterior_draws[[1]])[8326]#turnover[1]

turnover1.PostMedian<-smod$quantiles["turnover[1]","50%"]
turnover1.Post.2.5<-smod$quantiles["turnover[1]","2.5%"]
turnover1.Post.97.5<-smod$quantiles["turnover[1]","97.5%"]

turnover2.PostMedian<-smod$quantiles["turnover[2]","50%"]
turnover2.Post.2.5<-smod$quantiles["turnover[2]","2.5%"]
turnover2.Post.97.5<-smod$quantiles["turnover[2]","97.5%"]

turnover3.PostMedian<-smod$quantiles["turnover[3]","50%"]
turnover3.Post.2.5<-smod$quantiles["turnover[3]","2.5%"]
turnover3.Post.97.5<-smod$quantiles["turnover[3]","97.5%"]

turnover4.PostMedian<-smod$quantiles["turnover[4]","50%"]
turnover4.Post.2.5<-smod$quantiles["turnover[4]","2.5%"]
turnover4.Post.97.5<-smod$quantiles["turnover[4]","97.5%"]


turnoverDF<-data.frame(year=c(2017,2018,2019,2021),
                       turnoverMedian=c(turnover1.PostMedian,
                                        turnover2.PostMedian,
                                        turnover3.PostMedian,
                                        turnover4.PostMedian),
                       turnover.2.5=c(turnover1.Post.2.5,
                                      turnover2.Post.2.5,
                                      turnover3.Post.2.5,
                                      turnover4.Post.2.5),
                       turnover.97.5=c(turnover1.Post.97.5,
                                       turnover2.Post.97.5,
                                       turnover3.Post.97.5,
                                       turnover4.Post.97.5))

p6<-ggplot(turnoverDF, aes(x=year, y=turnoverMedian))+
  #geom_point(aes(x=year, y=turnoverMedian), col="red")+
  #geom_line(aes(x=year, y=turnoverMedian), col="red")+
  #geom_line(aes(x=year, y=turnover.2.5), linetype="dashed")+
  #geom_line(aes(x=year, y=turnover.97.5), linetype="dashed")+
  geom_ribbon(aes(x=year, ymin=turnover.2.5, ymax=turnover.97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=turnoverMedian), colour="blue")+
  xlab("Year")+ylab("Turnover from Previous Year")+
  theme(
    axis.text = element_text(size = 14), # Adjust all axis tick label sizes
    axis.title = element_text(size = 14) # Adjust all axis title sizes
  )+my.theme2


ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_2016_2021_PredictedTurnover.png",
       plot = p6,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p7<-grid.arrange(p1,p4,p5,p6,nrow=2,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_nov2025/Hedley 2016 2021/RShybridMinusHedley_2016_2021_PredictedOccPhiGammaTurnover.png",
       plot = p7,
       width = 12,
       height = 12,
       units = c("in"),
       dpi = 300)

#Get AUC

which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
colnames(posterior_draws[[1]])[4407]#"z[1,1]"
colnames(posterior_draws[[1]])[4522]#"z[1,2]"
colnames(posterior_draws[[1]])[4981]#"z[115,5]"

which(substr(colnames(posterior_draws[[1]]),1,4)=="muZ[")
colnames(posterior_draws[[1]])[4982]#"muZ[1,2]"
colnames(posterior_draws[[1]])[5441]#"muZ[115,5]"
str(posterior_draws[[1]][,"muZ[2,2]"])#no estimates for season 1
posterior_draws[[1]][,"muZ[1,2]"]#no estimates for season 1

niter<-ni
nthin<-nt
nchain<-nc
nsite<-115
nyear<-5

library(ggmcmc)
allmuZs<-posterior_draws[1:3][,4982:5441]#
ggd <- ggs(allmuZs)
allZs<-posterior_draws[1:3][,4522:4981]
ggd2 <- ggs(allZs)#starting from year 2

nstore <- niter / nthin
psi.pred <- array(dim=c(nsite, nstore*nchain))
Z.est <- array(dim=c(nsite, nstore*nchain))

#2017
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 2, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7510286
median(AUClist)#[1] 0.7527779
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6681985 0.7527779 0.8252308

value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~96% above 0.7 

#2018
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 3, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.8302872
median(AUClist)#[1] 0.8323643
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7596200 0.8323643 0.8896669

value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~100% above 0.7 

#2019
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", 4, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.8132914
median(AUClist)#[1] 0.8149225
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7377858 0.8149225 0.8798920

value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~100% above 0.7 

#2021
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("muZ[", j, ",", nyear, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
  indx <- which(ggd2$Parameter == paste("z[", j, ",", nyear, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}  


require(ROCR)
#One iteration: looks good
i=100
psi.vals <- psi.pred[, i]#1:115
Z.vals <- Z.est[, i]#1:115
pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
#A prediction instance with 115 data points
perf <- performance(pred, "auc") 
#A performance instance 'Area under the ROC curve'
AUC <- perf@y.values[[1]]
perf <- performance(pred, "tpr","fpr") 
plot(perf)
fpr <- perf@x.values[[1]]#  values plotted on x axis
#length of perf@x.values[[1]] can vary with i
tpr <- perf@y.values[[1]]#  values plotted on y axis  
#length of perf@y.values[[1]] can vary with i

AUC1 <-list()
for (i in 1:(nstore*nchain)){
  psi.vals <- psi.pred[, i]#
  Z.vals <- Z.est[, i]#
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #Prediction instances
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7874757 instead of 0.7872093
median(AUClist)#[1] 0.7897165 instead of 0.789225
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.7045003 0.7897165 0.8568167

value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~98% above 0.7 
