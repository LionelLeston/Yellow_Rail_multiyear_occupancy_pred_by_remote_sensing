library(sf)
library(tidyverse)
library(terra)
finalYERAsitesShapefile<-st_read("0_data/0_raw/point counts/allYERAsites_June2025.shp")
nrow(finalYERAsitesShapefile)#283
finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  st_transform("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#reprojected to same CRS as human footprint layers
box<-st_bbox(finalYERAsitesShapefile)
box[1]#xmin
box[2]#ymin
box[3]#xmax
box[4]#ymax

#https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r
#create a global object with same extent as point count layer
#(actually used conventional seismic line since its extent is a bit bigger)
hf2022_bb <- matrix(c((box[1]-10000), (box[4]+10000),
                      (box[3]+10000), (box[4]+10000),
                      (box[3]+10000), (box[2]-10000),
                      (box[1]-10000), (box[2]-10000),
                      (box[1]-10000), (box[4]+10000)), byrow = TRUE, ncol = 2) %>%
  list() %>% 
  st_polygon() %>% 
  st_sfc(., crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

st_write(hf2022_bb, "0_data/0_raw/point counts/OSM_YERA_colext_studyarea.shp", append=FALSE)
#this will be our study area

studyArea<-st_read("0_data/0_raw/point counts/OSM_YERA_colext_studyarea.shp")
st_bbox(studyArea)
#   xmin       ymin     xmax       ymax 
#568323.8 6040557.5  820183.4 6400476.0

wetlandProb<-rast("0_data/0_raw/wetlands/WetlandProbabilityOS.tif")
#using probability threshold of 0.35 to classify upland versus wetland
wetlandProb<-wetlandProb%>%
  crop(studyArea)
plot(wetlandProb["WetlandProbabilityOS"])

upland<-wetlandProb$WetlandProbabilityOS
upland[upland>0.35]<-NA
upland[!upland>0.35]<-1
distanceUpland<-terra::distance(upland)

#create buffers around points
b150<-st_buffer(finalYERAsitesShapefile,150)#150 m circle buffer
b565<-st_buffer(finalYERAsitesShapefile,565)#565 m circle buffer

#extract mean harvest year value in each buffer
pointList<-list()
for (i in 1:nrow(finalYERAsitesShapefile)){
  pointData<-finalYERAsitesShapefile[i,]
  #Extracted and summarized within 150 m and 565 m of point
  try(wetlandProb_150m <- terra::extract(wetlandProb$WetlandProbabilityOS, b150[i,], fun="mean", na.rm=TRUE))
  try(pointData$wetlandProb_150m<-wetlandProb_150m$WetlandProbabilityOS)
  try(wetlandProb_565m <- terra::extract(wetlandProb$WetlandProbabilityOS, b565[i,], fun="mean", na.rm=TRUE))
  try(pointData$wetlandProb_565m<-wetlandProb_565m$WetlandProbabilityOS)

  try(distUpland <- terra::extract(distanceUpland, finalYERAsitesShapefile[i,]))
  try(pointData$distUpland<-distUpland$WetlandProbabilityOS)
  
  try(pointList[[i]]<-pointData)
  print(paste0("Line ",i," done."))
}
pointWetlandProb<-do.call(rbind,pointList)
summaries<-data.frame(pointWetlandProb)
write.csv(summaries, file="0_data/1_processed/3_wetlands/wetlandProb_uplandDist.csv")


hf2022_bb<-st_read("0_data/0_raw/point counts/OSM_YERA_colext_studyarea.shp")
#then make a grid whose cells (~100x100 km) can be used to intersect shapefile and make tiles
hf2022_grid_10x10 <- st_make_grid(studyArea, n = c(10, 10), 
                                  crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs", 
                                  what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))


#wetland polygons extracted from CWCS wetland layer
wetlands<-st_read("0_data/0_raw/wetlands/AESRD_CWCSmergedwetlandinven2018layer.shp")

openWater<-wetlands%>%
  filter(CWCS_Class %in% c("Open Water"))%>%
  st_transform(crs=st_crs(studyArea))%>%
  st_intersection(studyArea)
bog<-wetlands%>%
  filter(CWCS_Class %in% c("Bog"))%>%
  st_transform(crs=st_crs(studyArea))%>%
  st_intersection(studyArea)
fen<-wetlands%>%
  filter(CWCS_Class %in% c("Fen"))%>%
  st_transform(crs=st_crs(studyArea))%>%
  st_intersection(studyArea)
marsh<-wetlands%>%
  filter(CWCS_Class %in% c("Marsh"))%>%
  st_transform(crs=st_crs(studyArea))%>%
  st_intersection(studyArea)

rm(wetlands)
gc()

#wrap in a function
areanoyeardist<-function(pointfile, polyfile){
  pkey.sf.summary<-list()
  #pkey.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    pkey.sf.i<-pointfile[i,]#pkey.sf.demo[i,]
    SS.i<-pkey.sf.i$locationC
    #YEAR.i<-pkey.sf.i$year
    nSitesUsed.i<-pkey.sf.i$nSitesUsed
    clusterID1.i<-pkey.sf.i$clusterID1
    clusterID2.i<-pkey.sf.i$clusterID2
    clusterID.i<-pkey.sf.i$clusterID
    
    #1. create 150-m and 565-m buffers
    b150<-st_buffer(pkey.sf.i,150)
    p.150<-st_intersection(polyfile, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection

    #area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Recalc_Area))
    
    b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
    p.565<-st_intersection(polyfile, b565)
    p.565$Recalc_Area<-st_area(p.565)#recalculate areas of clipped polygons
    area.565 <- ifelse((nrow(p.565)==0), 0, (p.565 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    #area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Recalc_Area))

    #to filter amount of footprint type within 150 and 565 m
    pkey.sf.summary[[i]]<-data.frame(locationC=SS.i,
                                     #year=YEAR.i,
                                     nSitesUsed=nSitesUsed.i,
                                     clusterID1=clusterID1.i,
                                     clusterID2=clusterID2.i,
                                     clusterID=clusterID.i,
                                     AREA150=area.150,
                                     PROP150=area.150/(3.14*150*150),
                                     AREA565=area.565,
                                     PROP565=area.565/(3.14*565*565))
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(rbind,pkey.sf.summary)
  return(summaries)
}

#now go through each cell in the grid and use it to intersect part of the
#the point counts. For each group of point counts, generate a buffer 
#equal to maxdist and use this multiple-point-counts-buffer to intersect
#wetland layers (fen, bog, marsh, and open water);
#then use this final intersection and point-counts-group
#in the areanoyeardist function
#pointcountgroup<-list()

#fen
fen<-fen%>%
  st_transform(crs=st_crs(hf2022_cell))
#Uses areanoyeardist function
maxdist<-1000
finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  st_transform(crs=st_crs(hf2022_cell))
for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(finalYERAsitesShapefile, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(fen, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areanoyeardist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areanoyeardist(pointfile,polyfile)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/3_wetlands/fen/area_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of features calculated within cell ",j))
  }
}


#bog
bog<-bog%>%
  st_transform(crs=st_crs(hf2022_cell))
#Uses areanoyeardist function
maxdist<-1000
finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  st_transform(crs=st_crs(hf2022_cell))
for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(finalYERAsitesShapefile, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(bog, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areanoyeardist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areanoyeardist(pointfile,polyfile)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/3_wetlands/bog/area_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of features calculated within cell ",j))
  }
}

#marsh
marsh<-marsh%>%
  st_transform(crs=st_crs(hf2022_cell))
#Uses areanoyeardist function
maxdist<-1000
for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(finalYERAsitesShapefile, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(marsh, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areanoyeardist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areanoyeardist(pointfile,polyfile)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/3_wetlands/marsh/area_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of features calculated within cell ",j))
  }
}

#open water
openWater<-openWater%>%
  st_transform(crs=st_crs(hf2022_cell))
#Uses areanoyeardist function
maxdist<-1000
finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
  st_transform(crs=st_crs(hf2022_cell))
for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(finalYERAsitesShapefile, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(openWater, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areanoyeardist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areanoyeardist(pointfile,polyfile)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/3_wetlands/open water/area_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of features calculated within cell ",j))
  }
}


#Now for each wetland type, combine the smaller files into one large file
#for each wetland type's area within 150 and 565 m of
#each point count survey. 

#fen
names1<-list.files(path = "0_data/1_processed/3_wetlands/fen/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/3_wetlands/fen/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/3_wetlands/fen-area-allsites.csv")

#bog
names1<-list.files(path = "0_data/1_processed/3_wetlands/bog/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/3_wetlands/bog/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/3_wetlands/bog-area-allsites.csv")

#marsh
names1<-list.files(path = "0_data/1_processed/3_wetlands/marsh/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/3_wetlands/marsh/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/3_wetlands/marsh-area-allsites.csv")

#open water
names1<-list.files(path = "0_data/1_processed/3_wetlands/open water/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/3_wetlands/open water/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/3_wetlands/open-water-area-allsites.csv")

#merge files (code assumes same order of location and year among all files)
dense.linear<-read.csv("0_data/1_processed/0_spatial data extracted/human footprint/dense-linear-areadist.age-allsurveys.csv", header=TRUE)
dense.linear$P150.denseLinear<-dense.linear$PROP150
dense.linear$P565.denseLinear<-dense.linear$PROP565
dense.linear$DIST.denseLinear<-dense.linear$NEAR.DIST
str(dense.linear)
m1<-dense.linear%>%
  select(locationC,year,P150.denseLinear,P565.denseLinear,DIST.denseLinear)


road<-read.csv("0_data/1_processed/0_spatial data extracted/human footprint/road-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.Road<-road$PROP150
m1$P565.Road<-road$PROP565
m1$DIST.Road<-road$NEAR.DIST


well.high<-read.csv("0_data/1_processed/0_spatial data extracted/human footprint/well-high-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.wellHigh<-well.high$PROP150
m1$P565.wellHigh<-well.high$PROP565
m1$DIST.wellHigh<-well.high$NEAR.DIST


well.low<-read.csv("0_data/1_processed/0_spatial data extracted/human footprint/well-low-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.wellLow<-well.low$PROP150
m1$P565.wellLow<-well.low$PROP565
m1$DIST.wellLow<-well.low$NEAR.DIST


plants.mines<-read.csv("0_data/1_processed/0_spatial data extracted/human footprint/plant-mine-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.plantMines<-plants.mines$PROP150
m1$P565.plantMines<-plants.mines$PROP565
m1$DIST.plantMines<-plants.mines$NEAR.DIST

pipelines.powerlines<-read.csv("0_data/1_processed/0_spatial data extracted/human footprint/pipeline-powerline-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.wide.linear<-pipelines.powerlines$PROP150
m1$P565.wide.linear<-pipelines.powerlines$PROP565
m1$DIST.wide.linear<-pipelines.powerlines$NEAR.DIST

write.csv(m1, file="0_data/1_processed/0_spatial data extracted/0_all-human-footprint.csv")

