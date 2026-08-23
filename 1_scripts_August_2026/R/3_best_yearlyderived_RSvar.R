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
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric,
                                     HabSuit4ClassNumeric))

#get yearly derived RS variables at locations in siteCovFinal
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

yearlyDW_grass<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,grass)%>%
  pivot_wider(names_from=year,values_from=grass)

yearlyDW_shrub<-dw_all%>%
  arrange(locationC,year)%>%
  select(locationC,year,shrub_and_scrub)%>%
  pivot_wider(names_from=year,values_from=shrub_and_scrub)


#get more siteCovs for initial occupancy models
siteCovFinal$woody<-yearlyDW_woody$`2016`
siteCovFinal$water<-yearlyDW_water$`2016`
siteCovFinal$flooded<-yearlyDW_flooded$`2016`
siteCovFinal$bare<-yearlyDW_bare$`2016`
siteCovFinal$grass<-yearlyDW_grass$`2016`
siteCovFinal$shrub<-yearlyDW_shrub$`2016`


umf<-unmarkedMultFrame(y=obsFiltFinal,
                       siteCovs = siteCovFinal,
                       obsCovs=list(time=matrix(4,115,25)),#8 years*5 rounds
                       yearlySiteCovs=list(Y_habSuit=yearlyHabSuit4Class[,2:6],
                                           Y_DW_bare=yearlyDW_bare[,2:6],
                                           Y_DW_flooded=yearlyDW_flooded[,2:6],
                                           Y_DW_grass=yearlyDW_grass[,2:6],
                                           Y_DW_shrub=yearlyDW_shrub[,2:6],
                                           Y_DW_water=yearlyDW_water[,2:6],
                                           Y_DW_woody=yearlyDW_woody[,2:6]),
                       numPrimary=5) 

#Scale predictors and create nonlinear terms
umf@siteCovs$bare.s<-scale(umf@siteCovs$bare, center = TRUE, scale = TRUE)
umf@siteCovs$bare.s2<-umf@siteCovs$bare.s^2
umf@siteCovs$bare.sqrt<-sqrt(umf@siteCovs$bare)

umf@siteCovs$flooded.s<-scale(umf@siteCovs$flooded, center = TRUE, scale = TRUE)
umf@siteCovs$flooded.s2<-umf@siteCovs$flooded.s^2
umf@siteCovs$flooded.sqrt<-sqrt(umf@siteCovs$flooded)

umf@siteCovs$grass.s<-scale(umf@siteCovs$grass, center = TRUE, scale = TRUE)
umf@siteCovs$grass.s2<-umf@siteCovs$grass.s^2
umf@siteCovs$grass.sqrt<-sqrt(umf@siteCovs$grass)

umf@siteCovs$shrub.s<-scale(umf@siteCovs$shrub, center = TRUE, scale = TRUE)
umf@siteCovs$shrub.s2<-umf@siteCovs$shrub.s^2
umf@siteCovs$shrub.sqrt<-sqrt(umf@siteCovs$shrub)

umf@siteCovs$water.s<-scale(umf@siteCovs$water, center = TRUE, scale = TRUE)
umf@siteCovs$water.s2<-umf@siteCovs$water.s^2
umf@siteCovs$water.sqrt<-sqrt(umf@siteCovs$water)

umf@siteCovs$woody.s<-scale(umf@siteCovs$woody, center = TRUE, scale = TRUE)
umf@siteCovs$woody.s2<-umf@siteCovs$woody.s^2
umf@siteCovs$woody.sqrt<-sqrt(umf@siteCovs$woody)

umf@yearlySiteCovs$Y_DW_bare.s<-scale(umf@yearlySiteCovs$Y_DW_bare, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_DW_bare.s2<-umf@yearlySiteCovs$Y_DW_bare.s^2
umf@yearlySiteCovs$Y_DW_bare.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_bare)

umf@yearlySiteCovs$Y_DW_flooded.s<-scale(umf@yearlySiteCovs$Y_DW_flooded, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_DW_flooded.s2<-umf@yearlySiteCovs$Y_DW_flooded.s^2
umf@yearlySiteCovs$Y_DW_flooded.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_flooded)

umf@yearlySiteCovs$Y_DW_grass.s<-scale(umf@yearlySiteCovs$Y_DW_grass, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_DW_grass.s2<-umf@yearlySiteCovs$Y_DW_grass.s^2
umf@yearlySiteCovs$Y_DW_grass.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_grass)

umf@yearlySiteCovs$Y_DW_shrub.s<-scale(umf@yearlySiteCovs$Y_DW_shrub, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_DW_shrub.s2<-umf@yearlySiteCovs$Y_DW_shrub.s^2
umf@yearlySiteCovs$Y_DW_shrub.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_shrub)

umf@yearlySiteCovs$Y_DW_water.s<-scale(umf@yearlySiteCovs$Y_DW_water, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_DW_water.s2<-umf@yearlySiteCovs$Y_DW_water.s^2
umf@yearlySiteCovs$Y_DW_water.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_water)

umf@yearlySiteCovs$Y_DW_woody.s<-scale(umf@yearlySiteCovs$Y_DW_woody, center = TRUE, scale = TRUE)
umf@yearlySiteCovs$Y_DW_woody.s2<-umf@yearlySiteCovs$Y_DW_woody.s^2
umf@yearlySiteCovs$Y_DW_woody.sqrt<-sqrt(umf@yearlySiteCovs$Y_DW_woody)


umf@siteCovs$HabSuit4ClassNumeric<-as.factor(umf@siteCovs$HabSuit4ClassNumeric)
umf@yearlySiteCovs$Y_habSuit<-as.factor(umf@yearlySiteCovs$Y_habSuit)

model1<-colext(psiformula= ~HabSuit4ClassNumeric, 
               gammaformula =  ~ Y_habSuit, 
               epsilonformula = ~ Y_habSuit,
               pformula = ~ 1, umf)
summary(model1)#AIC: 1682.768 (Null model)

WoodyModel1<-colext(psiformula= ~HabSuit4ClassNumeric+woody.s, 
                      gammaformula =  ~ Y_habSuit+Y_DW_woody.s, 
                      epsilonformula = ~ Y_habSuit+Y_DW_woody.s,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1684.869
WoodyModel2<-colext(psiformula= ~HabSuit4ClassNumeric+woody.s+woody.s2, 
                      gammaformula =  ~ Y_habSuit+Y_DW_woody.s+Y_DW_woody.s2, 
                      epsilonformula = ~ Y_habSuit+Y_DW_woody.s+Y_DW_woody.s2,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1688.854
WoodyModel3<-colext(psiformula= ~HabSuit4ClassNumeric+woody.sqrt, 
                      gammaformula =  ~ Y_habSuit+Y_DW_woody.sqrt, 
                      epsilonformula = ~ Y_habSuit+Y_DW_woody.sqrt,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1684.888
#Woody cover doesn't improve model

ShrubModel1<-colext(psiformula= ~HabSuit4ClassNumeric+shrub.s, 
                    gammaformula =  ~ Y_habSuit+Y_DW_shrub.s, 
                    epsilonformula = ~ Y_habSuit+Y_DW_shrub.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1673.399
ShrubModel2<-colext(psiformula= ~HabSuit4ClassNumeric+shrub.s+shrub.s2, 
                    gammaformula =  ~ Y_habSuit+Y_DW_shrub.s+Y_DW_shrub.s2, 
                    epsilonformula = ~ Y_habSuit+Y_DW_shrub.s+Y_DW_shrub.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1673.913
ShrubModel3<-colext(psiformula= ~HabSuit4ClassNumeric+shrub.sqrt, 
                    gammaformula =  ~ Y_habSuit+Y_DW_shrub.sqrt, 
                    epsilonformula = ~ Y_habSuit+Y_DW_shrub.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1674.411
#Shrub cover does improve model
#Number values don't seem to be correlated with habitat suitability
cor(umf@siteCovs$shrub.s, siteCovFinal$HabSuit4ClassNumeric)#0.06105979
siteCovFinal%>%
  group_by(HabSuit4ClassNumeric)%>%
  summarize(shrub=mean(shrub))
#  HabSuit4ClassNumeric shrub
#<dbl> <dbl>
#1                    1 0.147
#2                    2 0.154
#3                    3 0.122
#4                    4 0.172

BareModel1<-colext(psiformula= ~HabSuit4ClassNumeric+bare.s, 
                    gammaformula =  ~ Y_habSuit+Y_DW_bare.s, 
                    epsilonformula = ~ Y_habSuit+Y_DW_bare.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1674.124
BareModel2<-colext(psiformula= ~HabSuit4ClassNumeric+bare.s+bare.s2, 
                    gammaformula =  ~ Y_habSuit+Y_DW_bare.s+Y_DW_bare.s2, 
                    epsilonformula = ~ Y_habSuit+Y_DW_bare.s+Y_DW_bare.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1676.403
BareModel3<-colext(psiformula= ~HabSuit4ClassNumeric+bare.sqrt, 
                    gammaformula =  ~ Y_habSuit+Y_DW_bare.sqrt, 
                    epsilonformula = ~ Y_habSuit+Y_DW_bare.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1673.579 Large or missing SE values
#Bare Ground appears to improve model

GrassModel1<-colext(psiformula= ~HabSuit4ClassNumeric+grass.s, 
                    gammaformula =  ~ Y_habSuit+Y_DW_grass.s, 
                    epsilonformula = ~ Y_habSuit+Y_DW_grass.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1672.806
GrassModel2<-colext(psiformula= ~HabSuit4ClassNumeric+grass.s+grass.s2, 
                    gammaformula =  ~ Y_habSuit+Y_DW_grass.s+Y_DW_grass.s2, 
                    epsilonformula = ~ Y_habSuit+Y_DW_grass.s+Y_DW_grass.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1673.028
GrassModel3<-colext(psiformula= ~HabSuit4ClassNumeric+grass.sqrt, 
                    gammaformula =  ~ Y_habSuit+Y_DW_grass.sqrt, 
                    epsilonformula = ~ Y_habSuit+Y_DW_grass.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1673.847
#grass cover does improve model

FloodedModel1<-colext(psiformula= ~HabSuit4ClassNumeric+flooded.s, 
                    gammaformula =  ~ Y_habSuit+Y_DW_flooded.s, 
                    epsilonformula = ~ Y_habSuit+Y_DW_flooded.s,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1678.833
FloodedModel2<-colext(psiformula= ~HabSuit4ClassNumeric+flooded.s+flooded.s2, 
                    gammaformula =  ~ Y_habSuit+Y_DW_flooded.s+Y_DW_flooded.s2, 
                    epsilonformula = ~ Y_habSuit+Y_DW_flooded.s+Y_DW_flooded.s2,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1677.713
FloodedModel3<-colext(psiformula= ~HabSuit4ClassNumeric+flooded.sqrt, 
                    gammaformula =  ~ Y_habSuit+Y_DW_flooded.sqrt, 
                    epsilonformula = ~ Y_habSuit+Y_DW_flooded.sqrt,
                    pformula = ~ 1, umf, na.action = na.fail)#AIC: 1679.399
#flooded land cover improves model fit: colonization +

WaterModel1<-colext(psiformula= ~HabSuit4ClassNumeric+water.s, 
                      gammaformula =  ~ Y_habSuit+Y_DW_water.s, 
                      epsilonformula = ~ Y_habSuit+Y_DW_water.s,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1685.356
WaterModel2<-colext(psiformula= ~HabSuit4ClassNumeric+water.s+water.s2, 
                      gammaformula =  ~ Y_habSuit+Y_DW_water.s+Y_DW_water.s2, 
                      epsilonformula = ~ Y_habSuit+Y_DW_water.s+Y_DW_water.s2,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1686.252
WaterModel3<-colext(psiformula= ~HabSuit4ClassNumeric+water.sqrt, 
                      gammaformula =  ~ Y_habSuit+Y_DW_water.sqrt, 
                      epsilonformula = ~ Y_habSuit+Y_DW_water.sqrt,
                      pformula = ~ 1, umf, na.action = na.fail)#AIC: 1685.593 Large or missing SE values


cor(siteCovFinal$grass, siteCovFinal$woody)
#weak negative correlation (-0.3260957)
cor(siteCovFinal$grass, siteCovFinal$shrub)
#somewhat strongly positively correlated (0.6028941)
cor(siteCovFinal$grass, siteCovFinal$bare)
#weakly positively correlated (0.3665885)
cor(siteCovFinal$shrub, siteCovFinal$bare)
#weakly positively correlated (0.4810908)

#parallel processing and using dredge didn't work, so dropping it

#RS model without Hedley layer
RSonlyModel<-colext(psiformula= ~grass.s+bare.s+shrub.s, 
                         gammaformula =  ~ Y_DW_flooded.s, 
                         epsilonformula = ~ Y_DW_shrub.s,
                         pformula = ~ 1, umf)#AIC: 1681.127
#grass.s and shrub.s could be dropped from psi model

RSonlyModel<-colext(psiformula= ~bare.s, 
                    gammaformula =  ~ Y_DW_flooded.s, 
                    epsilonformula = ~ Y_DW_shrub.s,
                    pformula = ~ 1, umf)#AIC: 1667.803 flooded might be dropped

RSonlyModel<-colext(psiformula= ~bare.s, 
                    gammaformula =  ~ 1, 
                    epsilonformula = ~ Y_DW_shrub.s,
                    pformula = ~ 1, umf)#AIC: 1670.515

FullModel<-colext(psiformula= ~HabSuit4ClassNumeric+bare.s, 
                  gammaformula =  ~ Y_DW_flooded.s, 
                  epsilonformula = ~ Y_habSuit+shrub.s,
                  pformula = ~ 1, umf)#AIC: 1663.081

mb.gof.test(FullModel, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Number of seasons:  5 

#Chi-square statistic:
#  Season 1 Season 2 Season 3 Season 4 Season 5 
# 87.1882  99.3747  34.9360  31.8261  62.1324 
#Total chi-square = 315.4575 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#122  157  173  189  276 
#Estimate of c-hat = 1.81 

FullModelWOHedley<-colext(psiformula= ~bare.s, 
                  gammaformula =  ~ Y_DW_flooded.s, 
                  epsilonformula = ~ shrub.s,
                  pformula = ~ 1, umf)
#Initial (logit-scale):
#Estimate    SE    z  P(>|z|)
#(Intercept)    0.492 0.273 1.80 0.071238
#bare.s         1.302 0.378 3.44 0.000573
#Colonization (logit-scale):
#  Estimate    SE     z  P(>|z|)
#(Intercept)      -0.908 0.179 -5.09 3.61e-07
#Y_DW_flooded.s    0.260 0.146  1.78 7.55e-02
#Extinction (logit-scale):
#  Estimate    SE      z P(>|z|)
#(Intercept)   0.0296 0.181  0.163  0.8705
#shrub.s      -0.3563 0.170 -2.094  0.0362
#Detection (logit-scale):
#  Estimate     SE   z  P(>|z|)
#0.477 0.0723 6.6 4.14e-11
#AIC: 1668.324 
#Number of sites: 115

mb.gof.test(FullModelWOHedley, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5 
#85.6847  97.0455  35.3407  32.8538  63.2435 
#Total chi-square = 314.1681 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#109  158  172  190  289 
#Estimate of c-hat = 1.8 

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
ggsave(filename="2_outputs/model outputs/unmarked/best yearly derived/occu.ByRegion2016_2021_HedleyXYearlyDerived.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)
#ggsave(filename="2_outputs/1_unmarked_outputs/November/Hedley 2016 2021/occu.ByRegion2016_2021_HedleyXYearlyDerived.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)

