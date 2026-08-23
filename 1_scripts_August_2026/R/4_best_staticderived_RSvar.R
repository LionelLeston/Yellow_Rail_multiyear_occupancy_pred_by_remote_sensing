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
finalYERAsitesShapefile$lat<-st_coordinates(finalYERAsitesShapefile)[,2]
finalYERAsitesShapefile$long<-st_coordinates(finalYERAsitesShapefile)[,1]
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
#"yera2022_5"
#,"yera2023_1","yera2023_2","yera2023_3","yera2023_4", 
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
table(m2$clusterID)#123 sites, 8 clusters
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
#113 sites with 3+ years of data
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
  select(locationC, lat, long, clusterID, clusterIDfinal, HabSuit4Class)%>%
  mutate(HabSuit4ClassNumeric=as.numeric(HabSuit4Class))
siteCovFinal$geometry<-NULL

HabSuit4ClassNumeric<-siteCovFinal$HabSuit4ClassNumeric
yearlyHabSuit4Class<-as.matrix(cbind(HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric))
Latitude<-siteCovFinal$lat
yearlyLatitude<-as.matrix(cbind(Latitude,Latitude,Latitude,Latitude,Latitude))

Longitude<-siteCovFinal$long
yearlyLongitude<-as.matrix(cbind(Longitude,Longitude,Longitude,Longitude,Longitude))

#get static derived RS variables at locations in siteCovFinal

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

Wet150<-wet$wetlandProb_150m
Y_wet150<-as.matrix(cbind(Wet150,Wet150,Wet150,Wet150,Wet150))

Wet565<-wet$wetlandProb_565m
Y_wet565<-as.matrix(cbind(Wet565,Wet565,Wet565,Wet565,Wet565))

DistUpland<-wet$distUpland
Y_distUpland<-as.matrix(cbind(DistUpland,DistUpland,DistUpland,DistUpland,DistUpland))

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
yearlyDW_flooded<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,flooded_vegetation)%>%
  pivot_wider(names_from=year,values_from=flooded_vegetation)
#These are year-specific values. Get an averaged value across years for each site
staticDW_flooded<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,flooded_vegetation)%>%
  group_by(locationC)%>%
  summarise(sflooded_vegetation=mean(flooded_vegetation))

yearlyDW_water<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,water)%>%
  pivot_wider(names_from=year,values_from=water)
#These are year-specific values. Get an averaged value across years for each site
staticDW_water<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,water)%>%
  group_by(locationC)%>%
  summarise(swater=mean(water))

yearlyDW_woody<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,woody)%>%
  pivot_wider(names_from=year,values_from=woody)
#These are year-specific values. Get an averaged value across years for each site
staticDW_woody<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,woody)%>%
  group_by(locationC)%>%
  summarise(swoody=mean(woody))

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

yearlyDW_shrub<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,shrub_and_scrub)%>%
  pivot_wider(names_from=year,values_from=shrub_and_scrub)
#These are year-specific values. Get an averaged value across years for each site
staticDW_shrub<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,shrub_and_scrub)%>%
  group_by(locationC)%>%
  summarise(sshrub_and_scrub=mean(shrub_and_scrub))


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

swoody<-staticDW_woody$swoody
swater<-staticDW_water$swater
sflooded<-staticDW_flooded$sflooded_vegetation
sbare<-staticDW_bare$sbare
sgrass<-staticDW_grass$sgrass
sshrub<-staticDW_shrub$sshrub_and_scrub

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
sshrubSY<-as.matrix(cbind(sshrub,
                          sshrub,
                          sshrub,
                          sshrub,
                          sshrub))
swoodySY<-as.matrix(cbind(swoody,
                          swoody,
                          swoody,
                          swoody,
                          swoody))
swaterSY<-as.matrix(cbind(swater,
                          swater,
                          swater,
                          swater,
                          swater))
sfloodedSY<-as.matrix(cbind(sflooded,
                            sflooded,
                            sflooded,
                            sflooded,
                            sflooded))

#get more siteCovs for initial occupancy models
siteCovFinal$swoody<-staticDW_woody$swoody
siteCovFinal$swater<-staticDW_water$swater
siteCovFinal$sflooded<-staticDW_flooded$sflooded_vegetation
siteCovFinal$sbare<-staticDW_bare$sbare
siteCovFinal$sgrass<-staticDW_grass$sgrass
siteCovFinal$sshrub<-staticDW_shrub$sshrub_and_scrub


umf<-unmarkedMultFrame(y=obsFiltFinal,
                       siteCovs = siteCovFinal,
                       obsCovs=list(time=matrix(4,115,25)),#8 years*5 rounds
                       yearlySiteCovs=list(Y_habSuit=yearlyHabSuit4Class,
                                           #same as "static" HabSuit4Class would be
                                           Y_latitude=yearlyLatitude,
                                           Y_longitude=yearlyLongitude,
                                           Y_elev=Y_elev,
                                           Y_swi=Y_swi,
                                           Y_tpi=Y_tpi,
                                           Y_mrvbf=Y_mrvbf,
                                           Y_wet150=Y_wet150,
                                           Y_wet565=Y_wet565,
                                           Y_distUpland=Y_distUpland,
                                           S_DW_bare=sbareSY,
                                           S_DW_flooded=sfloodedSY,
                                           S_DW_grass=sgrassSY,
                                           S_DW_shrub=sshrubSY,
                                           S_DW_water=swaterSY,
                                           S_DW_woody=swoodySY),
                       numPrimary=5) 

#Scale predictors and create nonlinear terms
umf@siteCovs$elev.s<-scale(umf@siteCovs$elev, center = TRUE, scale = TRUE)
umf@siteCovs$elev.s2<-umf@siteCovs$elev.s^2
umf@siteCovs$elev.sqrt<-sqrt(umf@siteCovs$elev)

umf@siteCovs$lat.s<-scale(umf@siteCovs$lat, center = TRUE, scale = TRUE)
umf@siteCovs$lat.s2<-umf@siteCovs$lat.s^2
umf@siteCovs$lat.sqrt<-sqrt(umf@siteCovs$lat)

umf@siteCovs$swi.s<-scale(umf@siteCovs$swi, center = TRUE, scale = TRUE)
umf@siteCovs$swi.s2<-umf@siteCovs$swi.s^2
umf@siteCovs$swi.sqrt<-sqrt(umf@siteCovs$swi)

umf@siteCovs$tpi.s<-scale(umf@siteCovs$tpi, center = TRUE, scale = TRUE)
umf@siteCovs$tpi.s2<-umf@siteCovs$tpi.s^2
umf@siteCovs$tpi.sqrt<-sqrt(umf@siteCovs$tpi)

umf@siteCovs$mrvbf.s<-scale(umf@siteCovs$mrvbf, center = TRUE, scale = TRUE)
umf@siteCovs$mrvbf.s2<-umf@siteCovs$mrvbf.s^2
umf@siteCovs$mrvbf.sqrt<-sqrt(umf@siteCovs$mrvbf)

umf@siteCovs$wet150.s<-scale(umf@siteCovs$wetlandProb_150m, center = TRUE, scale = TRUE)
umf@siteCovs$wet150.s2<-umf@siteCovs$wet150.s^2
umf@siteCovs$wet150.sqrt<-sqrt(umf@siteCovs$wetlandProb_150m)

umf@siteCovs$wet565.s<-scale(umf@siteCovs$wetlandProb_565m, center = TRUE, scale = TRUE)
umf@siteCovs$wet565.s2<-umf@siteCovs$wet565.s^2
umf@siteCovs$wet565.sqrt<-sqrt(umf@siteCovs$wetlandProb_565m)

umf@siteCovs$distUpland.s<-scale(umf@siteCovs$distUpland, center = TRUE, scale = TRUE)
umf@siteCovs$distUpland.s2<-umf@siteCovs$distUpland.s^2
umf@siteCovs$distUpland.sqrt<-sqrt(umf@siteCovs$distUpland)

umf@siteCovs$sbare.s<-scale(umf@siteCovs$sbare, center = TRUE, scale = TRUE)
umf@siteCovs$sbare.s2<-umf@siteCovs$sbare.s^2
umf@siteCovs$sbare.sqrt<-sqrt(umf@siteCovs$sbare)

umf@siteCovs$sflooded.s<-scale(umf@siteCovs$sflooded, center = TRUE, scale = TRUE)
umf@siteCovs$sflooded.s2<-umf@siteCovs$sflooded.s^2
umf@siteCovs$sflooded.sqrt<-sqrt(umf@siteCovs$sflooded)

umf@siteCovs$sgrass.s<-scale(umf@siteCovs$sgrass, center = TRUE, scale = TRUE)
umf@siteCovs$sgrass.s2<-umf@siteCovs$sgrass.s^2
umf@siteCovs$sgrass.sqrt<-sqrt(umf@siteCovs$sgrass)

umf@siteCovs$sshrub.s<-scale(umf@siteCovs$sshrub, center = TRUE, scale = TRUE)
umf@siteCovs$sshrub.s2<-umf@siteCovs$sshrub.s^2
umf@siteCovs$sshrub.sqrt<-sqrt(umf@siteCovs$sshrub)

umf@siteCovs$swater.s<-scale(umf@siteCovs$swater, center = TRUE, scale = TRUE)
umf@siteCovs$swater.s2<-umf@siteCovs$swater.s^2
umf@siteCovs$swater.sqrt<-sqrt(umf@siteCovs$swater)

umf@siteCovs$swoody.s<-scale(umf@siteCovs$swoody, center = TRUE, scale = TRUE)
umf@siteCovs$swoody.s2<-umf@siteCovs$swoody.s^2
umf@siteCovs$swoody.sqrt<-sqrt(umf@siteCovs$swoody)

umf@yearlySiteCovs$Y_distUpland.s<-scale(umf@yearlySiteCovs$Y_distUpland, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_distUpland.s2<-umf@yearlySiteCovs$Y_distUpland.s^2
umf@yearlySiteCovs$Y_distUpland.sqrt<-sqrt(umf@yearlySiteCovs$Y_distUpland)

umf@yearlySiteCovs$Y_wet150.s<-scale(umf@yearlySiteCovs$Y_wet150, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_wet150.s2<-umf@yearlySiteCovs$Y_wet150.s^2
umf@yearlySiteCovs$Y_wet150.sqrt<-sqrt(umf@yearlySiteCovs$Y_wet150)

umf@yearlySiteCovs$Y_wet565.s<-scale(umf@yearlySiteCovs$Y_wet565, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_wet565.s2<-umf@yearlySiteCovs$Y_wet565.s^2
umf@yearlySiteCovs$Y_wet565.sqrt<-sqrt(umf@yearlySiteCovs$Y_wet565)

umf@yearlySiteCovs$Y_elev.s<-scale(umf@yearlySiteCovs$Y_elev, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_elev.s2<-umf@yearlySiteCovs$Y_elev.s^2
umf@yearlySiteCovs$Y_elev.sqrt<-sqrt(umf@yearlySiteCovs$Y_elev)

umf@yearlySiteCovs$Y_latitude.s<-scale(umf@yearlySiteCovs$Y_latitude, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_latitude.s2<-umf@yearlySiteCovs$Y_latitude.s^2
umf@yearlySiteCovs$Y_latitude.sqrt<-sqrt(umf@yearlySiteCovs$Y_latitude)


umf@yearlySiteCovs$Y_swi.s<-scale(umf@yearlySiteCovs$Y_swi, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_swi.s2<-umf@yearlySiteCovs$Y_swi.s^2
umf@yearlySiteCovs$Y_swi.sqrt<-sqrt(umf@yearlySiteCovs$Y_swi)

umf@yearlySiteCovs$Y_tpi.s<-scale(umf@yearlySiteCovs$Y_tpi, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_tpi.s2<-umf@yearlySiteCovs$Y_tpi.s^2
umf@yearlySiteCovs$Y_tpi.sqrt<-sqrt(umf@yearlySiteCovs$Y_tpi)

umf@yearlySiteCovs$Y_mrvbf.s<-scale(umf@yearlySiteCovs$Y_mrvbf, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_mrvbf.s2<-umf@yearlySiteCovs$Y_mrvbf.s^2
umf@yearlySiteCovs$Y_mrvbf.sqrt<-sqrt(umf@yearlySiteCovs$Y_mrvbf)


umf@yearlySiteCovs$S_DW_bare.s<-scale(umf@yearlySiteCovs$S_DW_bare, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$S_DW_bare.s2<-umf@yearlySiteCovs$S_DW_bare.s^2
umf@yearlySiteCovs$S_DW_bare.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_bare)

umf@yearlySiteCovs$S_DW_flooded.s<-scale(umf@yearlySiteCovs$S_DW_flooded, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$S_DW_flooded.s2<-umf@yearlySiteCovs$S_DW_flooded.s^2
umf@yearlySiteCovs$S_DW_flooded.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_flooded)

umf@yearlySiteCovs$S_DW_grass.s<-scale(umf@yearlySiteCovs$S_DW_grass, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$S_DW_grass.s2<-umf@yearlySiteCovs$S_DW_grass.s^2
umf@yearlySiteCovs$S_DW_grass.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_grass)

umf@yearlySiteCovs$S_DW_shrub.s<-scale(umf@yearlySiteCovs$S_DW_shrub, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$S_DW_shrub.s2<-umf@yearlySiteCovs$S_DW_shrub.s^2
umf@yearlySiteCovs$S_DW_shrub.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_shrub)

umf@yearlySiteCovs$S_DW_water.s<-scale(umf@yearlySiteCovs$S_DW_water, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$S_DW_water.s2<-umf@yearlySiteCovs$S_DW_water.s^2
umf@yearlySiteCovs$S_DW_water.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_water)

umf@yearlySiteCovs$S_DW_woody.s<-scale(umf@yearlySiteCovs$S_DW_woody, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$S_DW_woody.s2<-umf@yearlySiteCovs$S_DW_woody.s^2
umf@yearlySiteCovs$S_DW_woody.sqrt<-sqrt(umf@yearlySiteCovs$S_DW_woody)


umf@siteCovs$HabSuit4ClassNumeric<-as.factor(umf@siteCovs$HabSuit4ClassNumeric)
umf@yearlySiteCovs$Y_habSuit<-as.factor(umf@yearlySiteCovs$Y_habSuit)

model1<-colext(psiformula= ~HabSuit4ClassNumeric, 
               gammaformula =  ~ Y_habSuit, 
               epsilonformula = ~ Y_habSuit,
               pformula = ~ 1, umf)
summary(model1)#AIC: 1694.376 (Null model)

WoodyModel1<-colext(psiformula= ~HabSuit4ClassNumeric+swoody.s, 
                      gammaformula =  ~ Y_habSuit+S_DW_woody.s, 
                      epsilonformula = ~ Y_habSuit+S_DW_woody.s,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1696.852
WoodyModel2<-colext(psiformula= ~HabSuit4ClassNumeric+swoody.s+swoody.s2, 
                      gammaformula =  ~ Y_habSuit+S_DW_woody.s+S_DW_woody.s2, 
                      epsilonformula = ~ Y_habSuit+S_DW_woody.s+S_DW_woody.s2,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1697.391
WoodyModel3<-colext(psiformula= ~HabSuit4ClassNumeric+swoody.sqrt, 
                      gammaformula =  ~ Y_habSuit+S_DW_woody.sqrt, 
                      epsilonformula = ~ Y_habSuit+S_DW_woody.sqrt,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1696.384
#Woody cover doesn't improve model

ShrubModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sshrub.s, 
                    gammaformula =  ~ Y_habSuit+S_DW_shrub.s, 
                    epsilonformula = ~ Y_habSuit+S_DW_shrub.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1686.358
ShrubModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sshrub.s+sshrub.s2, 
                    gammaformula =  ~ Y_habSuit+S_DW_shrub.s+S_DW_shrub.s2, 
                    epsilonformula = ~ Y_habSuit+S_DW_shrub.s+S_DW_shrub.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1690.488
ShrubModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sshrub.sqrt, 
                    gammaformula =  ~ Y_habSuit+S_DW_shrub.sqrt, 
                    epsilonformula = ~ Y_habSuit+S_DW_shrub.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1687.445
#Shrub cover does improve model; AIC values not as good as for yearly shrub cover

BareModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sbare.s, 
                    gammaformula =  ~ Y_habSuit+S_DW_bare.s, 
                    epsilonformula = ~ Y_habSuit+S_DW_bare.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1693.476
BareModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sbare.s+sbare.s2, 
                    gammaformula =  ~ Y_habSuit+S_DW_bare.s+S_DW_bare.s2, 
                    epsilonformula = ~ Y_habSuit+S_DW_bare.s+S_DW_bare.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1698.009
BareModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sbare.sqrt, 
                    gammaformula =  ~ Y_habSuit+S_DW_bare.sqrt, 
                    epsilonformula = ~ Y_habSuit+S_DW_bare.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1693.876 Large or missing SE values
#Bare Ground only appears to slightly improve model; AICs not as good as when using yearly values

GrassModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sgrass.s, 
                    gammaformula =  ~ Y_habSuit+S_DW_grass.s, 
                    epsilonformula = ~ Y_habSuit+S_DW_grass.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1679.785
GrassModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sgrass.s+sgrass.s2, 
                    gammaformula =  ~ Y_habSuit+S_DW_grass.s+S_DW_grass.s2, 
                    epsilonformula = ~ Y_habSuit+S_DW_grass.s+S_DW_grass.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1668.914
GrassModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sgrass.sqrt, 
                    gammaformula =  ~ Y_habSuit+S_DW_grass.sqrt, 
                    epsilonformula = ~ Y_habSuit+S_DW_grass.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1681.76 large SE values
#grass cover does improve model; lower AICs than when using yearly values

FloodedModel1<-colext(psiformula= ~HabSuit4ClassNumeric+sflooded.s, 
                    gammaformula =  ~ Y_habSuit+S_DW_flooded.s, 
                    epsilonformula = ~ Y_habSuit+S_DW_flooded.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1696.418
FloodedModel2<-colext(psiformula= ~HabSuit4ClassNumeric+sflooded.s+sflooded.s2, 
                    gammaformula =  ~ Y_habSuit+S_DW_flooded.s+S_DW_flooded.s2, 
                    epsilonformula = ~ Y_habSuit+S_DW_flooded.s+S_DW_flooded.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1700.481
FloodedModel3<-colext(psiformula= ~HabSuit4ClassNumeric+sflooded.sqrt, 
                    gammaformula =  ~ Y_habSuit+S_DW_flooded.sqrt, 
                    epsilonformula = ~ Y_habSuit+S_DW_flooded.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1696.847
#flooded land cover doesn't improve model fit; AIC values higher than when using yearly values for flooded land

WaterModel1<-colext(psiformula= ~HabSuit4ClassNumeric+swater.s, 
                      gammaformula =  ~ Y_habSuit+S_DW_water.s, 
                      epsilonformula = ~ Y_habSuit+S_DW_water.s,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1698.469
WaterModel2<-colext(psiformula= ~HabSuit4ClassNumeric+swater.s+swater.s2, 
                      gammaformula =  ~ Y_habSuit+S_DW_water.s+S_DW_water.s2, 
                      epsilonformula = ~ Y_habSuit+S_DW_water.s+S_DW_water.s2,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1699.978
WaterModel3<-colext(psiformula= ~HabSuit4ClassNumeric+swater.sqrt, 
                      gammaformula =  ~ Y_habSuit+S_DW_water.sqrt, 
                      epsilonformula = ~ Y_habSuit+S_DW_water.sqrt,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1698.885 Large or missing SE values
#water cover doesn't improve model fit; AIC values higher than when using yearly values for open water

ElevModel1<-colext(psiformula= ~HabSuit4ClassNumeric+elev.s, 
                    gammaformula =  ~ Y_habSuit+Y_elev.s, 
                    epsilonformula = ~ Y_habSuit+Y_elev.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1682.707
ElevModel2<-colext(psiformula= ~HabSuit4ClassNumeric+elev.s+elev.s2, 
                    gammaformula =  ~ Y_habSuit+Y_elev.s+Y_elev.s2, 
                    epsilonformula = ~ Y_habSuit+Y_elev.s+Y_elev.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1686.607
ElevModel3<-colext(psiformula= ~HabSuit4ClassNumeric+elev.sqrt, 
                    gammaformula =  ~ Y_habSuit+Y_elev.sqrt, 
                    epsilonformula = ~ Y_habSuit+Y_elev.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1682.368 
#Elevation improves model fit

cor(umf@siteCovs$elev, umf@siteCovs$lat)#-0.9646267
LatModel1<-colext(psiformula= ~HabSuit4ClassNumeric+lat.s, 
                  gammaformula =  ~ Y_habSuit+Y_latitude.s, 
                  epsilonformula = ~ Y_habSuit+Y_latitude.s,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1666.605
#When latitude is included, - effect on psi, + effect on epsilon
#Lower AIC than elevation
LatModel2<-colext(psiformula= ~HabSuit4ClassNumeric+lat.s+lat.s2, 
                  gammaformula =  ~ Y_habSuit+Y_latitude.s+Y_latitude.s2, 
                  epsilonformula = ~ Y_habSuit+Y_latitude.s+Y_latitude.s2,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1668.426
#When latitude is included, marginal effects
#Lower AIC than elevation
LatModel3<-colext(psiformula= ~HabSuit4ClassNumeric+lat.sqrt, 
                  gammaformula =  ~ Y_habSuit+Y_latitude.sqrt, 
                  epsilonformula = ~ Y_habSuit+Y_latitude.sqrt,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1665.642 
#When latitude is included, - effect on psi, + effect on epsilon
#Lower AIC than elevation

LatElevModel1<-colext(psiformula= ~HabSuit4ClassNumeric+elev.s+lat.s, 
                   gammaformula =  ~ Y_habSuit+Y_elev.s+Y_latitude.s, 
                   epsilonformula = ~ Y_habSuit+Y_elev.s+Y_latitude.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1661.793
#When latitude is included, effect of elevation on psi and epsilon goes away and is negative for gamma
LatElevModel2<-colext(psiformula= ~HabSuit4ClassNumeric+elev.s+elev.s2+lat.s+lat.s2, 
                   gammaformula =  ~ Y_habSuit+Y_elev.s+Y_elev.s2+Y_latitude.s+Y_latitude.s2, 
                   epsilonformula = ~ Y_habSuit+Y_elev.s+Y_elev.s2+Y_latitude.s+Y_latitude.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1669.094
#When latitude is included, effect of elevation is zero
LatElevModel3<-colext(psiformula= ~HabSuit4ClassNumeric+elev.sqrt+lat.sqrt, 
                   gammaformula =  ~ Y_habSuit+Y_elev.sqrt+Y_latitude.sqrt, 
                   epsilonformula = ~ Y_habSuit+Y_elev.sqrt+Y_latitude.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1662.388 
#When latitude is included, effect of elevation on psi and epsilon goes away and is negative for gamma
#Include latitude, not elevation in final model

mrvbfModel1<-colext(psiformula= ~HabSuit4ClassNumeric+mrvbf.s, 
                   gammaformula =  ~ Y_habSuit+Y_mrvbf.s, 
                   epsilonformula = ~ Y_habSuit+Y_mrvbf.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1692.054
mrvbfModel2<-colext(psiformula= ~HabSuit4ClassNumeric+mrvbf.s+mrvbf.s2, 
                   gammaformula =  ~ Y_habSuit+Y_mrvbf.s+Y_mrvbf.s2, 
                   epsilonformula = ~ Y_habSuit+Y_mrvbf.s+Y_mrvbf.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1694.075
mrvbfModel3<-colext(psiformula= ~HabSuit4ClassNumeric+mrvbf.sqrt, 
                   gammaformula =  ~ Y_habSuit+Y_mrvbf.sqrt, 
                   epsilonformula = ~ Y_habSuit+Y_mrvbf.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1691.831 
#Some improvement

tpiModel1<-colext(psiformula= ~HabSuit4ClassNumeric+tpi.s, 
                   gammaformula =  ~ Y_habSuit+Y_tpi.s, 
                   epsilonformula = ~ Y_habSuit+Y_tpi.s,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1689.146
tpiModel2<-colext(psiformula= ~HabSuit4ClassNumeric+tpi.s+tpi.s2, 
                   gammaformula =  ~ Y_habSuit+Y_tpi.s+Y_tpi.s2, 
                   epsilonformula = ~ Y_habSuit+Y_tpi.s+Y_tpi.s2,
                   pformula = ~ 1, umf, na.action = na.fail)#AIC: 1692.36
tpiModel3<-colext(psiformula= ~HabSuit4ClassNumeric+tpi.sqrt, 
                   gammaformula =  ~ Y_habSuit+Y_tpi.sqrt, 
                   epsilonformula = ~ Y_habSuit+Y_tpi.sqrt,
                   pformula = ~ 1, umf, na.action = na.fail)#missing values 
#TPI improves model fit

swiModel1<-colext(psiformula= ~HabSuit4ClassNumeric+swi.s, 
                  gammaformula =  ~ Y_habSuit+Y_swi.s, 
                  epsilonformula = ~ Y_habSuit+Y_swi.s,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1687.455
swiModel2<-colext(psiformula= ~HabSuit4ClassNumeric+swi.s+swi.s2, 
                  gammaformula =  ~ Y_habSuit+Y_swi.s+Y_swi.s2, 
                  epsilonformula = ~ Y_habSuit+Y_swi.s+Y_swi.s2,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1691.162
swiModel3<-colext(psiformula= ~HabSuit4ClassNumeric+swi.sqrt, 
                  gammaformula =  ~ Y_habSuit+Y_swi.sqrt, 
                  epsilonformula = ~ Y_habSuit+Y_swi.sqrt,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1687.26 missing values 

wet150Model1<-colext(psiformula= ~HabSuit4ClassNumeric+wet150.s, 
                  gammaformula =  ~ Y_habSuit+Y_wet150.s, 
                  epsilonformula = ~ Y_habSuit+Y_wet150.s,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1695.321
wet150Model2<-colext(psiformula= ~HabSuit4ClassNumeric+wet150.s+wet150.s2, 
                  gammaformula =  ~ Y_habSuit+Y_wet150.s+Y_wet150.s2, 
                  epsilonformula = ~ Y_habSuit+Y_wet150.s+Y_wet150.s2,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1691.584
wet150Model3<-colext(psiformula= ~HabSuit4ClassNumeric+wet150.sqrt, 
                  gammaformula =  ~ Y_habSuit+Y_wet150.sqrt, 
                  epsilonformula = ~ Y_habSuit+Y_wet150.sqrt,
                  pformula = ~ 1, umf, na.action = na.fail)#AIC: 1694.348

wet565Model1<-colext(psiformula= ~HabSuit4ClassNumeric+wet565.s, 
                     gammaformula =  ~ Y_habSuit+Y_wet565.s, 
                     epsilonformula = ~ Y_habSuit+Y_wet565.s,
                     pformula = ~ 1, umf, na.action = na.fail)#AIC: 1689.629
wet565Model2<-colext(psiformula= ~HabSuit4ClassNumeric+wet565.s+wet565.s2, 
                     gammaformula =  ~ Y_habSuit+Y_wet565.s+Y_wet565.s2, 
                     epsilonformula = ~ Y_habSuit+Y_wet565.s+Y_wet565.s2,
                     pformula = ~ 1, umf, na.action = na.fail)#AIC: 1689.639
wet565Model3<-colext(psiformula= ~HabSuit4ClassNumeric+wet565.sqrt, 
                     gammaformula =  ~ Y_habSuit+Y_wet565.sqrt, 
                     epsilonformula = ~ Y_habSuit+Y_wet565.sqrt,
                     pformula = ~ 1, umf, na.action = na.fail)#AIC: 1689.142
#wet565 does a better job than wet150 but effect is opposite of expected

distUplandModel1<-colext(psiformula= ~HabSuit4ClassNumeric+distUpland.s, 
                     gammaformula =  ~ Y_habSuit+Y_distUpland.s, 
                     epsilonformula = ~ Y_habSuit+Y_distUpland.s,
                     pformula = ~ 1, umf, na.action = na.fail)#AIC: 1693.132
distUplandModel2<-colext(psiformula= ~HabSuit4ClassNumeric+distUpland.s+distUpland.s2, 
                     gammaformula =  ~ Y_habSuit+Y_distUpland.s+Y_distUpland.s2, 
                     epsilonformula = ~ Y_habSuit+Y_distUpland.s+Y_distUpland.s2,
                     pformula = ~ 1, umf, na.action = na.fail)#AIC: 1695.593
distUplandModel3<-colext(psiformula= ~HabSuit4ClassNumeric+distUpland.sqrt, 
                     gammaformula =  ~ Y_habSuit+Y_distUpland.sqrt, 
                     epsilonformula = ~ Y_habSuit+Y_distUpland.sqrt,
                     pformula = ~ 1, umf, na.action = na.fail)#AIC: 1692.3


cor(umf@siteCovs$grass.s, umf@siteCovs$elev.s)
#weak correlation (0.004156551)
cor(umf@siteCovs$tpi.s, umf@siteCovs$elev.s)
#weak correlation (-0.3939487)


#parallel processing and using dredge didn't work, so dropping it

#RS model without Hedley layer
RSonlyModel<-colext(psiformula= ~sgrass.s+sgrass.s2+sshrub.s+wet565.sqrt+lat.sqrt, #drop shrub, wet565?
                    gammaformula =  ~ Y_distUpland.sqrt+Y_mrvbf.sqrt+Y_swi.sqrt, #drop MRVBF, distUpland?
                    epsilonformula = ~ Y_tpi.s+Y_latitude.sqrt+Y_wet565.sqrt, #drop wet565?
                    pformula = ~ 1, umf)#AIC: 1623.48  Large or missing SE values
#shrub and wet565 could be dropped

RSonlyModel<-colext(psiformula= ~sgrass.s+sgrass.s2+lat.sqrt, #drop wet565, shrub
                    gammaformula =  ~ Y_swi.sqrt, #drop MRVBF, distUpland?
                    epsilonformula = ~ Y_tpi.s+Y_latitude.sqrt, #drop wet565?
                    pformula = ~ 1, umf)#AIC: 1617.925

FullModel<-colext(psiformula= ~sgrass.s+sgrass.s2+lat.sqrt, 
                  gammaformula =  ~ Y_swi.sqrt, 
                  epsilonformula = ~ Y_habSuit+Y_tpi.s+Y_latitude.sqrt,
                  pformula = ~ 1, umf)#AIC: 1617.634
#Basically the RS variables by themselves do a good enough job without the Hedley layer

mb.gof.test(RSonlyModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
# 89.3010  98.9516  34.5826  31.0401  64.2118 
#Total chi-square = 318.0871 
#Number of bootstrap samples = 1000
#P-value = 0.002
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#113  155  168  184  332 
#Estimate of c-hat = 1.84 

mb.gof.test(FullModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
#89.2598  99.3165  34.4851  30.6841  63.8252 
#Total chi-square = 317.5708 
#Number of bootstrap samples = 1000
#P-value = 0.001
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#119  154  168  185  345 
#Estimate of c-hat = 1.84 

FullModelWOHedley<-colext(psiformula= ~sgrass.s+sgrass.s2+lat.sqrt, 
                  gammaformula =  ~ Y_swi.sqrt, 
                  epsilonformula = ~ Y_tpi.s+Y_latitude.sqrt,
                  pformula = ~ 1, umf)#AIC: 1617.925
#Basically the RS variables by themselves do a good enough job without the Hedley layer

mb.gof.test(FullModelWOHedley, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)

#Goodness-of-fit for dynamic occupancy model
#Number of seasons:  5 
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#89.3010  98.9516  34.5826  31.0401  64.2118 
#Total chi-square = 318.0871 
#Number of bootstrap samples = 1000
#P-value = 0.003
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#118  156  171  189  356 
#Estimate of c-hat = 1.82 

#Get estimates of uncertainty in yearly occu. probability via 
#nonparametric bootstrapping
#use bootstrap to get estimates of psi (occu.) from derived parameter estimates

bootModel100x<-nonparboot(RSonlyModel, B=100)
str(bootModel100x@projected)
#num [1:2, 1:5, 1:115] 0.369 0.631 0.434 0.566 0.448 
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
ggsave(filename="2_outputs/model outputs/unmarked/best static derived/occu.ByRegion2016_2021_HedleyXStaticDerived.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)
#ggsave(filename="2_outputs/1_unmarked_outputs/November/Hedley 2016 2021/occu.ByRegion2016_2021_HedleyXStaticDerived.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)

