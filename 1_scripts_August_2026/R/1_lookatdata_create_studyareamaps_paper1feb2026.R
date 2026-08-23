library(tidyverse)
library(sf)
library(ggplot2)
library(grid)
library(gridExtra)
library(terra)
library(ggspatial)
library(tidyterra)
library(RColorBrewer)
library(gtable)
library(egg)                            
library(ggpubr)
library(ggthemes)
library(svglite)
library(spatialEco)
library(leaflet)

finalYERAsitesShapefile<-st_read("0_data/0_raw/point counts/allYERAsites_June2025.shp")
finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  mutate(clusterIDfinal=clusterID)%>%
  mutate(clusterIDfinal=case_when(locationC %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) ~ "2_8b",
                                  locationC %in% c(64,65,58) ~ "7_8b",
                                  locationC %in% c(60,61,62,63) ~ "7_8c",
                                  locationC %in% c(53,54,55,56) ~ "7_8d",
                                  !locationC %in% c(1:16,53:56,58,60:65) ~ clusterID))
levels(as.factor(finalYERAsitesShapefile$clusterIDfinal))

finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")%>%
  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#st_transform(crs=st_crs(majorRoads))

studyArea<-st_bbox(finalYERAsitesShapefile)
#   xmin       ymin     xmax       ymax 
#568323.8 6040557.5  820183.4 6400476.0


obs<-read.csv("0_data/0_raw/point counts/YERA_2025_COLEXT_DATA_DETECTIONHISTORYcombined.csv", header=TRUE)
#manually went through to combine detection histories stored on multiple rows
#for same sites, condensed detection histories to single row per site.
obsFilt<-obs%>%
  select("yera2016_1","yera2016_2","yera2016_3","yera2016_4", 
         "yera2016_5","yera2017_1","yera2017_2","yera2017_3","yera2017_4", 
         "yera2017_5","yera2018_1","yera2018_2","yera2018_3","yera2018_4", 
         "yera2018_5","yera2019_1","yera2019_2","yera2019_3","yera2019_4", 
         "yera2019_5","yera2021_1","yera2021_2","yera2021_3","yera2021_4", 
         "yera2021_5")%>% #January 25, 2026: drop 2013-2015, 2022-2023
  mutate(numVisits=rowSums(!is.na(.)))%>%
  mutate(locationC=obs$locationC)

m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  filter(!clusterID=="2_10")%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits


#Provinces
provs <-  st_read("0_data/0_raw/spatial/BCR_Terrestrial_master.shp")%>%
  filter(COUNTRY=="CANADA")%>%  #limit to Canada
  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")%>%  
  group_by(PROVINCE_S)%>%
  summarize()#dissolve provinces within Canadian provinces

alberta<-provs%>%
  filter(PROVINCE_S=="ALBERTA")  #limit to Alberta 
#already dissolved polygons within provinces, so no need to dissolve them here
plot(alberta["PROVINCE_S"])

#Alberta label
Alberta.Pt<- data.frame(ID = "ALBERTA" , lat = 57, lon = -113.4937)
Alberta.sf <- st_as_sf(x = Alberta.Pt,                         
                       coords = c("lon", "lat"),
                       crs = "+proj=longlat +ellps=WGS84 
+datum=WGS84 +no_defs +towgs84=0,0,0")%>%
  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

#Restrict to those points within the boundaries of the predictive Yellow Rail
#raster layer described in "Hedley, R. W., McLeod, L. J., Yip, D. A., Farr, D., 
#Knaga, P., Drake, K. L., & Bayne, E. (2020). Modeling the occurrence of the 
#Yellow Rail (Coturnicops noveboracensis) in the context of ongoing resource 
#development in the oil sands region of Alberta. Avian Conservation & Ecology, 
#15(1), 137."
Hedley2019<-rast("0_data/1_processed/6_HedleyEtAl2019/YERA_COLL_SUM_FINAL11.tif")
plot(Hedley2019$YERA_COLL_SUM_FINAL11)
HedleyLayer<-Hedley2019$YERA_COLL_SUM_FINAL11
HedleyLayer[!is.na(HedleyLayer)]<-1
#add cities
City <- data.frame(ID = c("Edmonton","Calgary","Fort McMurray"), 
                   lat = c(53.5461,51.0447,56.726389), 
                   lon = c(-113.4937,-114.0719,-111.380278))
city.sf <- st_as_sf(x = City,                         
                    coords = c("lon", "lat"),
                    crs = "+proj=longlat +ellps=WGS84
+datum=WGS84 +no_defs +towgs84=0,0,0")
city.sf<-city.sf%>%
  st_transform(crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
str(city.sf)

#Transform points to same projection as footprint
#finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
#  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")


#smaller map (study area)
st_bbox(m2)
#    xmin      ymin      xmax      ymax 
#668887.0 6040557.5  820183.4 6399753.5 

m2<-m2%>%mutate(clusterID1_8=case_when(clusterID=="3_6"~"1",
                                       clusterID=="4_5"~"2",
                                       clusterID=="7_5"~"3",
                                       #clusterID=="2_10"~"4",
                                       clusterID=="3_1"~"4",
                                       clusterID=="7_8"~"5",
                                       clusterID=="3_8"~"6",
                                       clusterID=="2_8"~"7"))
#normally cluster 4 would be 2_10 but that cluster was not analyzed
#since it is outside of Hedley layer
#just centroids of each cluster
m2_centroids<-m2%>%
  group_by(clusterID1_8)%>%
  mutate(Number=n())%>%
  select(clusterID1_8,Number)%>%
  group_by(clusterID1_8,Number) %>%
  summarise(st_union(geometry)) %>%
  st_centroid()

#plot larger map (Canada with Alberta highlighted)
map0<-ggplot() + 
  geom_sf(data = provs, size = 0.5, colour="black", fill="white")+
  geom_sf(data = alberta, size = 0.5, colour="grey", fill="lightyellow")+
  #geom_spatraster(data = HedleyLayer, show.legend = FALSE)+
  scale_fill_viridis_c(na.value = "transparent", option = "viridis")+
  #guides(fill=NA)+
  geom_sf(data = m2_centroids, 
          #shape=22, #square symbol instead of solid circle(default)
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  geom_sf(data = city.sf, 
          shape=22, #square symbol instead of solid circle(default)
          size = 5, 
          color = "black", 
          fill = "black") +
  geom_sf_text(data=city.sf, 
               aes(label = ID), 
               size=5, 
               nudge_y = -50000, nudge_x = 5000)+
  labs(y=" ", x=" ")+   
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white"))
#GitHub repos
ggsave("2_outputs/study area maps/CanadaMap.png", map0, width=6, height=8, units="in", dpi=300)
#Formerly
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/CanadaMap.png", map0, width=6, height=8, units="in", dpi=300)

#plot larger map (Alberta)
map1<-ggplot() + 
  geom_sf(data = alberta, size = 0.5, colour="grey", fill="lightyellow")+
  #geom_spatraster(data = HedleyLayer, show.legend = FALSE)+
  scale_fill_viridis_c(na.value = "transparent", option = "viridis")+
  #guides(fill=NA)+
  geom_sf(data = m2_centroids, 
          #shape=22, #square symbol instead of solid circle(default)
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  geom_sf(data = city.sf, 
          shape=22, #square symbol instead of solid circle(default)
          size = 5, 
          color = "black", 
          fill = "black") +
  geom_sf_text(data=city.sf, 
               aes(label = ID), 
               size=5, 
               nudge_y = -50000, nudge_x = 5000)+
  labs(y=" ", x=" ")+   
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white"))
#GitHub repos
ggsave("2_outputs/study area maps/AlbertaMap.png", map1, width=6, height=8, units="in", dpi=300)
#Formerly
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/AlbertaMap.png", map1, width=6, height=8, units="in", dpi=300)

map2<-ggplot() + 
  geom_sf(data = alberta, size = 0.5, colour="grey", fill=NA)+
  geom_spatraster(data = HedleyLayer, show.legend = FALSE)+
  scale_fill_viridis_c(na.value = "transparent", option = "viridis")+
  geom_sf(data = m2_centroids, 
          #shape=22, #square symbol instead of solid circle(default)
          aes(col = clusterID1_8, size=Number)) +
  scale_colour_manual(name="Cluster",
                      labels=c("1","2","3",
                               "4","5","6",
                               "7"),
                      #labels=c("2_10","2_8","3_1",
                      #         "3_6","3_8","4_5",
                      #         "7_5","7_8"),
                      values=c("#E69F00",
                               "#56B4E9",
                               "brown",
                               "#F0E442",
                               "#0072B2",
                               "#CC79A7",
                               "#DC3220"),
                      na.value="transparent")+
  geom_sf(data = city.sf, 
          shape=22, #square symbol instead of solid circle(default)
          size = 5, 
          color = "black", 
          fill = "black") +
  geom_sf_text(data=city.sf, aes(label = ID), size=6, nudge_y = 20000, nudge_x = 8000)+
  labs(y=" ", x=" ")+   
  xlim(600000,830000)+ylim(6000000,6405000)+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  annotation_scale(location="bl")+ #add scale bar
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white"))
#GitHub repos
ggsave("2_outputs/study area maps/paper1_WholeStudyAreaPlusHedley.png", map2, width=6, height=8, units="in", dpi=300)
#Former
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_WholeStudyAreaPlusHedley.png", map2, width=6, height=8, units="in", dpi=300)


  
m2.latlong<-m2%>%
  st_transform(crs="+proj=longlat +datum=WGS84")
pal <- colorFactor("viridis", levels = m2$clusterID)

map2b<-leaflet() %>%
  addTiles() %>%
  addProviderTiles("CartoDB.VoyagerLabelsUnder")%>%
  #addProviderTiles(providers$Esri.WorldImagery) %>%
  setView(lng = -114, lat = 56, zoom = 6) %>%
  #setView(lng = -111.3166667, lat = 56.4166667, zoom = 6) %>%
  addMarkers(lng = -3.7, lat = 40.4)%>%
  addCircleMarkers(data = m2.latlong, radius=~2, color=~ifelse(m2$clusterID=="2_10", "#000000",
                                                               ifelse(m2$clusterID=="2_8", "#E69F00",
                                                                      ifelse(m2$clusterID=="3_1", "#56B4E9", 
                                                                             ifelse(m2$clusterID=="3_6", "#009E73", 
                                                                                    ifelse(m2$clusterID=="3_8", "#F0E442", 
                                                                                           ifelse(m2$clusterID=="4_5", "#0072B2", 
                                                                                                  ifelse(m2$clusterID=="7_5", "#CC79A7","#DC3220"))))))),)

#Suppose I took cluster 2 or 3, made a 1km buffer around the points,
#made a minimum convex polygon around that buffered area, then used it to
#clip different raster layers that were used to predict YERA occupancy
#in this study. I then overlay the clipped rasters in a study area figure.
m2_justcluster7_5<-m2%>%filter(clusterID=="7_5")
m2_justcluster7_5_10km<-m2_justcluster7_5%>%
  st_buffer(10000)%>%
  st_union(by_feature=FALSE)%>%st_as_sf()
plot(m2_justcluster7_5_10km)

#E.g., Sentinel-1 predictor DPOL
#Import rasters 
dpol2016<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/DPOL_summer_2016.tif")
dpol2016_just7_5<-dpol2016%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
plot(dpol2016_just7_5$VH_Gamma)
range(values(dpol2016_just7_5$VH_Gamma), na.rm=TRUE)
hist(values(dpol2016_just7_5$VH_Gamma))
dpol2016_just7_5[dpol2016_just7_5>3]=3
dpol2016_just7_5[dpol2016_just7_5<(-3)]=(-3)
plot(dpol2016_just7_5$VH_Gamma)
range(values(dpol2016_just7_5$VH_Gamma), na.rm=TRUE)
hist(values(dpol2016_just7_5$VH_Gamma))

dpol2017<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/DPOL_summer_2017.tif")
dpol2017_just7_5<-dpol2017%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
plot(dpol2017_just7_5$VH_Gamma)
hist(values(dpol2017_just7_5$VH_Gamma))
dpol2017_just7_5[dpol2017_just7_5>3]=3
dpol2017_just7_5[dpol2017_just7_5<(-3)]=(-3)
plot(dpol2017_just7_5$VH_Gamma)


dpol2018<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/DPOL_summer_2018.tif")
dpol2018_just7_5<-dpol2018%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
plot(dpol2018_just7_5$VH_Gamma)
hist(values(dpol2018_just7_5$VH_Gamma))
dpol2018_just7_5[dpol2018_just7_5>3]=3
dpol2018_just7_5[dpol2018_just7_5<(-3)]=(-3)
plot(dpol2018_just7_5$VH_Gamma)

dpol2019<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/DPOL_summer_2019.tif")
dpol2019_just7_5<-dpol2019%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
plot(dpol2019_just7_5$VH_Gamma)
hist(values(dpol2019_just7_5$VH_Gamma))
dpol2019_just7_5[dpol2019_just7_5>3]=3
dpol2019_just7_5[dpol2019_just7_5<(-3)]=(-3)
plot(dpol2019_just7_5$VH_Gamma)

dpol2021<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 1 data/DPOL_summer_2021.tif")
dpol2021_just7_5<-dpol2021%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
plot(dpol2021_just7_5$VH_Gamma)
hist(values(dpol2021_just7_5$VH_Gamma))
dpol2021_just7_5[dpol2021_just7_5>3]=3
dpol2021_just7_5[dpol2021_just7_5<(-3)]=(-3)
plot(dpol2021_just7_5$VH_Gamma)

# 1. Combine your layers into a single multi-layer SpatRaster object
# (Can be a list of file paths, or individual SpatRaster objects combined with c())
r_stack <- c(dpol2016_just7_5, 
                  dpol2017_just7_5, 
                  dpol2018_just7_5, 
                  dpol2019_just7_5, 
                  dpol2021_just7_5) 

# 2. Calculate the cell-wise mean
dpol2016_2021_just7_5 <- mean(r_stack, na.rm = TRUE)
plot(dpol2016_2021_just7_5)

map_dpol2016_7_5<-ggplot() + 
  geom_spatraster(data=dpol2016_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="DPOL 2016") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_DPOL2016_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_DPOL2016_cluster7_5.png", 
       map_dpol2016_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, bg="transparent")

map_dpol2017_7_5<-ggplot() + 
  geom_spatraster(data=dpol2017_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="DPOL 2017") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_DPOL2017_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_DPOL2017_cluster7_5.png", 
       map_dpol2017_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300)

map_dpol2018_7_5<-ggplot() + 
  geom_spatraster(data=dpol2018_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="DPOL 2018") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_DPOL2018_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_DPOL2018_cluster7_5.png", 
       map_dpol2018_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, bg="transparent")

map_dpol2019_7_5<-ggplot() + 
  geom_spatraster(data = dpol2019_just7_5) +
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="DPOL 2019") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_DPOL2019_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_DPOL2019_cluster7_5.png", 
       map_dpol2019_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg = "transparent" # Crucial step!
       )

map_dpol2021_7_5<-ggplot() + 
  geom_spatraster(data=dpol2021_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="DPOL 2021") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_DPOL2021_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_DPOL2021_cluster7_5.png", 
       map_dpol2021_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg="transparent")



map_dpol2016_2021_7_5<-ggplot() + 
  geom_spatraster(data=dpol2016_2021_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Mean DPOL 2016-2021") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )
ggsave("2_outputs/study area maps/paper1_meanDPOL2016_2021_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_meanDPOL2016_2021_cluster7_5.png", 
       map_dpol2016_2021_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg="transparent")

#Now do Sentinel-2 predictor: NDMI
#Import rasters 
sentinel2_2016<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/YERASites-May2025-S2_monthly_2016_6.tif")
ndmi2016<-sentinel2_2016$NDMI
ndmi2016_just7_5<-ndmi2016%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(ndmi2016_just7_5$NDMI)

sentinel2_2017<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/YERASites-May2025-S2_monthly_2017_6.tif")
ndmi2017<-sentinel2_2017$NDMI
ndmi2017_just7_5<-ndmi2017%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(ndmi2017_just7_5$NDMI)


sentinel2_2018<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/YERASites-May2025-S2_monthly_2018_6.tif")
ndmi2018<-sentinel2_2018$NDMI
ndmi2018_just7_5<-ndmi2018%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(ndmi2018_just7_5$NDMI)

sentinel2_2019<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/YERASites-May2025-S2_monthly_2019_6.tif")
ndmi2019<-sentinel2_2019$NDMI
ndmi2019_just7_5<-ndmi2019%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(ndmi2019_just7_5$NDMI)

sentinel2_2021<-rast("0_data/1_processed/1_Sentinel_Landsat/Sentinel 2 data/YERASites-May2025-S2_monthly_2021_6.tif")
ndmi2021<-sentinel2_2021$NDMI
ndmi2021_just7_5<-ndmi2021%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(ndmi2021_just7_5$NDMI)

# 1. Combine your layers into a single multi-layer SpatRaster object
# (Can be a list of file paths, or individual SpatRaster objects combined with c())
r_stack <- c(ndmi2016_just7_5, 
             ndmi2017_just7_5, 
             ndmi2018_just7_5, 
             ndmi2019_just7_5, 
             ndmi2021_just7_5) 

# 2. Calculate the cell-wise mean
ndmi2016_2021_just7_5 <- mean(r_stack, na.rm = TRUE)
#plot(ndmi2016_2021_just7_5)

map_ndmi2016_7_5<-ggplot() + 
  geom_spatraster(data=ndmi2016_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="NDMI 2016") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_ndmi2016_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_ndmi2016_cluster7_5.png", 
       map_ndmi2016_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, bg="transparent")

map_ndmi2017_7_5<-ggplot() + 
  geom_spatraster(data=ndmi2017_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="NDMI 2017") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_ndmi2017_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_ndmi2017_cluster7_5.png", 
       map_ndmi2017_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300)

map_ndmi2018_7_5<-ggplot() + 
  geom_spatraster(data=ndmi2018_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="NDMI 2018") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_ndmi2018_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_ndmi2018_cluster7_5.png", 
       map_ndmi2018_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, bg="transparent")

map_ndmi2019_7_5<-ggplot() + 
  geom_spatraster(data = ndmi2019_just7_5) +
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="NDMI 2019") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_ndmi2019_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_ndmi2019_cluster7_5.png", 
       map_ndmi2019_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg = "transparent" # Crucial step!
)

map_ndmi2021_7_5<-ggplot() + 
  geom_spatraster(data=ndmi2021_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="NDMI 2021") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_ndmi2021_cluster7_5.png", 
#ggsave("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/paper1_ndmi2021_cluster7_5.png", 
       map_ndmi2021_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg="transparent")



map_ndmi2016_2021_7_5<-ggplot() + 
  geom_spatraster(data=ndmi2016_2021_just7_5)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Mean\nNDMI\n2016-2021") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_meanndmi2016_2021_cluster7_5.png", 
       map_ndmi2016_2021_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg="transparent")

#Now do Dynamic World predictor: grass
#Import rasters 
DynamicWorld_2016<-rast("0_data/1_processed/5_DynamicWorld/DynamicWorld2016_AOI.tif")
grass2016_just7_5<-DynamicWorld_2016%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(grass2016_just7_5$grass)

DynamicWorld_2017<-rast("0_data/1_processed/5_DynamicWorld/DynamicWorld2017_AOI.tif")
grass2017_just7_5<-DynamicWorld_2017%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(grass2017_just7_5$grass)

DynamicWorld_2018<-rast("0_data/1_processed/5_DynamicWorld/DynamicWorld2018_AOI.tif")
grass2018_just7_5<-DynamicWorld_2018%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(grass2018_just7_5$grass)

DynamicWorld_2019<-rast("0_data/1_processed/5_DynamicWorld/DynamicWorld2019_AOI.tif")
grass2019_just7_5<-DynamicWorld_2019%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(grass2019_just7_5$grass)

DynamicWorld_2021<-rast("0_data/1_processed/5_DynamicWorld/DynamicWorld2021_AOI.tif")
grass2021_just7_5<-DynamicWorld_2021%>%
  project(m2_justcluster7_5_10km)%>%
  crop(m2_justcluster7_5_10km, mask=TRUE)
#plot(grass2021_just7_5$grass)

# 1. Combine your layers into a single multi-layer SpatRaster object
# (Can be a list of file paths, or individual SpatRaster objects combined with c())
r_stack <- c(grass2016_just7_5, 
             grass2017_just7_5, 
             grass2018_just7_5, 
             grass2019_just7_5, 
             grass2021_just7_5) 

# 2. Calculate the cell-wise mean
grass2016_2021_just7_5 <- mean(r_stack, na.rm = TRUE)
#plot(grass2016_2021_just7_5)

map_grass2016_7_5<-ggplot() + 
  geom_spatraster(data=grass2016_just7_5$grass)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Grass 2016") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_grass2016_cluster7_5.png", 
       map_grass2016_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, bg="transparent")

map_grass2017_7_5<-ggplot() + 
  geom_spatraster(data=grass2017_just7_5$grass)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Grass 2017") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_grass2017_cluster7_5.png", 
       map_grass2017_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300)

map_grass2018_7_5<-ggplot() + 
  geom_spatraster(data=grass2018_just7_5$grass)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Grass 2018") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_grass2018_cluster7_5.png", 
       map_grass2018_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, bg="transparent")

map_grass2019_7_5<-ggplot() + 
  geom_spatraster(data = grass2019_just7_5$grass) +
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Grass 2019") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_grass2019_cluster7_5.png", 
       map_grass2019_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg = "transparent" # Crucial step!
)

map_grass2021_7_5<-ggplot() + 
  geom_spatraster(data=grass2021_just7_5$grass)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Grass 2021") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_grass2021_cluster7_5.png", 
       map_grass2021_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg="transparent")



map_grass2016_2021_7_5<-ggplot() + 
  geom_spatraster(data=grass2016_2021_just7_5$mean)+
  # Force NA values to be completely transparent
  scale_fill_viridis_c(na.value = "transparent", option = "viridis") +
  geom_sf(data = m2_justcluster7_5, 
          size = 2, 
          color = "darkred", 
          fill = "darkred") +
  labs(y = " ", x = " ", fill="Mean\ngrass\n2016-2021") +   
  # 1. Base theme goes FIRST
  theme_bw() +
  # 2. Overrides to strip out ALL backgrounds, borders, and grids go SECOND
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # Strip everything out for absolute transparency
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave("2_outputs/study area maps/paper1_meangrass2016_2021_cluster7_5.png", 
       map_grass2016_2021_7_5, 
       width=6, 
       height=8, 
       units="in", 
       dpi=300, 
       bg="transparent")
