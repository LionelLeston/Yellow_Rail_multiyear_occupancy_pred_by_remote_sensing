library(sf)
library(tidyverse)
#use a file with site and year data
#for dynamic occupancy models, need all site-year combinations, even if there are no 
#visits to a given site in some of those years
finalYERAsitesShapefile<-st_read("0_data/0_raw/point counts/allYERAsites_June2025.shp")
nrow(finalYERAsitesShapefile)#283
sites<-bind_rows(finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile,
                 finalYERAsitesShapefile)
year<-c(rep(2013,283),
        rep(2014,283),
        rep(2015,283),
        rep(2016,283),
        rep(2017,283),
        rep(2018,283),
        rep(2019,283),
        rep(2021,283),
        rep(2022,283),
        rep(2023,283))#2380
siteYears<-cbind(sites,year)
SSYR<-siteYears%>%
  st_transform("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#reprojected to same CRS as human footprint layers
box<-st_bbox(SSYR)
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

#then make a grid whose cells (~100x100 km) can be used to intersect shapefile and make tiles
hf2022_grid_10x10 <- st_make_grid(hf2022_bb, n = c(10, 10), 
                                  crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs", 
                                  what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))


#human footprint shapefiles extracted from ABMI geodatabase

#wrap in a function
areadist<-function(pointfile, polyfile, maxdist){
  pkey.sf.summary<-list()
  #pkey.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    pkey.sf.i<-pointfile[i,]#pkey.sf.demo[i,]
    SS.i<-pkey.sf.i$locationC
    YEAR.i<-pkey.sf.i$year
    nSitesUsed.i<-pkey.sf.i$nSitesUsed
    clusterID1.i<-pkey.sf.i$clusterID1
    clusterID2.i<-pkey.sf.i$clusterID2
    clusterID.i<-pkey.sf.i$clusterID
    
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(YEAR<YEAR.i)
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(pkey.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection

    #area.150<-ifelse((nrow(p.150)==0), 0, sum(p.150$Recalc_Area))
    
    b565<-st_buffer(pkey.sf.i,565)#1 square kilometre
    p.565<-st_intersection(polyfile.y, b565)
    p.565$Recalc_Area<-st_area(p.565)#recalculate areas of clipped polygons
    area.565 <- ifelse((nrow(p.565)==0), 0, (p.565 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    #area.565<-ifelse((nrow(p.565)==0), 0, sum(p.565$Recalc_Area))

    #to filter amount of footprint type within 150 and 565 m
    #2. select a maximum buffer distance around the point (X), beyond which footprint features are exceedingly unlikely to influence bird abundance within the point, then filter to only the footprint polygons that fall within that buffer distance.
    b.maxdist<-st_buffer(pkey.sf.i,maxdist)
    p.maxdist<-st_intersection(polyfile.y, b.maxdist)
    #3. estimate minimum distance among the time-and-buffer-filtered polygons
    #4. if nrow(time-and-buffer-filtered polygons) == 0, i.e. there is no footprint of a particular type within X m of a point count in the year of the survey, nearest distance to that footprint is capped at X.
    nearest<-ifelse((nrow(p.maxdist)==0), maxdist, 
                    min(st_distance(pkey.sf.i,p.maxdist)))
    pkey.sf.summary[[i]]<-data.frame(locationC=SS.i,
                                     year=YEAR.i,
                                     nSitesUsed=nSitesUsed.i,
                                     clusterID1=clusterID1.i,
                                     clusterID2=clusterID2.i,
                                     clusterID=clusterID.i,
                                     AREA150=area.150,
                                     PROP150=area.150/(3.14*150*150),
                                     AREA565=area.565,
                                     PROP565=area.565/(3.14*565*565),
                                     NEAR.DIST=nearest)
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(rbind,pkey.sf.summary)
  return(summaries)
}

#now go through each cell in the grid and use it to intersect part of the
#the point counts. For each group of point counts, generate a buffer 
#equal to maxdist and use this multiple-point-counts-buffer to intersect
#BADR footprint layers (dense linear, road, plant-mine, well-high, and 
#well-low); then use this final intersection and point-counts-group
#in the areadist function
#pointcountgroup<-list()

#Dense Linear
#conventional seismic lines
c.seismic<-st_read("0_data/0_raw/human footprint 2019/conventional.seismic2019corr.shp")
c.seismic<-st_transform(c.seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
#Low-impact seismic lines
li.seismic<-st_read("0_data/0_raw/human footprint 2019/low.impact.seismic2019corr.shp")
li.seismic<-st_transform(li.seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

dense.linear<-rbind(c.seismic, li.seismic)

maxdist<-1000
oldestyear<-min(dense.linear$YEAR[dense.linear$YEAR>0])
for (j in 1:nlevels(as.factor(hf2018_grid_10x10$ID))){
  hf2018_cell<-hf2018_grid_10x10[hf2018_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2018_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(dense.linear, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist.age function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist.age(pointfile,polyfile,maxdist,oldestyear)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/human footprint/dense linear/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

st_layers("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb")
#Industrial
industrial<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o08_Industrials_OSM_HFIe_2022")
industrial<-st_transform(industrial, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-8000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in 97:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(industrial, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/industrial areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}



#Mines
mines<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o07_Mines_OSM_HFIe_2022")
mines<-st_transform(mines, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-30000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in 98:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(mines, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/mines areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}



#Abandoned Wells
abandonedWells<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o16_WellsAbnd_OSM_HFIe_2022")
abandonedWells<-st_transform(abandonedWells, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-7000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in 97:nlevels(as.factor(hf2022_grid_10x10$ID))){
#for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(abandonedWells, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/abandoned wells areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Pipelines
pipelines<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o19_Pipelines_Centerlines_HFIe_2022")
pipelines<-st_transform(pipelines, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-50000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(pipelines, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/pipelines areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Transmission lines
powerlines<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o13_TransmissionLines_Centerlines_HFIe_2022")
powerlines<-st_transform(powerlines, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")


maxdist<-300000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in 9:19){
#for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(powerlines, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/powerlines areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Roads
roads<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o03_Roads_OSM_HFIe_2022")
roads<-roads%>%
  # filter(FEATURE_TY %in% c("ROAD-GRAVEL-1L",      
  #                          "ROAD-GRAVEL-2L",
  #                          "ROAD-PAVED-1L",       
  #                          "ROAD-PAVED-2L",       
  #                          "ROAD-PAVED-3L",       
  #                          "ROAD-PAVED-4L",       
  #                          "ROAD-PAVED-5L",       
  #                          "ROAD-PAVED-DIV"))%>%
  st_transform("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-8000
#oldestyear<-min(all.roads$YEAR[all.roads$YEAR>0], na.rm=TRUE)
for (j in c(86:87)){
#for (j in 26:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(roads, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/roads areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Active Wells
activeWells<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o09_WellsActive_OSM_HFIe_2022")
activeWells<-st_transform(activeWells, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-15000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in c(97)){
#for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(activeWells, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/active wells areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}

#Conventional Seismic
seismic<-st_read("0_data/0_raw/HFIeOSA_2022_v1_1.gdb/HFIeOSA_2022_v1_1.gdb", layer="o20_Seismic_Lines_and_Trails_OSM_HFIe_2022")
seismic<-st_transform(seismic, "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

maxdist<-3000#adjust this upward or downward depending on if/how many site-footprint distances exceed this limit
#Uses areadist function
for (j in 97:97){
#for (j in 1:nlevels(as.factor(hf2022_grid_10x10$ID))){
  hf2022_cell<-hf2022_grid_10x10[hf2022_grid_10x10$ID==j,]
  #Then isolate the points that fall within a given cell.
  pointfile<-st_intersection(SSYR, hf2022_cell)
  if (nrow(pointfile)>0){
    print(paste0("Points in cell ",j," isolated"))
    #buffer those cells by maxdist
    pointfilebuffer<-st_buffer(pointfile, maxdist)
    print(paste0("Points in cell ",j," buffered by maximum distance"))
    polyfile<-st_intersection(seismic, pointfilebuffer)
    print(paste0("Features within ",maxdist," m of points in cell ",j," isolated"))
    #Then run the areadist function, saving a separate file for each chunk.
    print(paste0("Area of /distance to features within ",maxdist," m of points in cell ",j," calculated"))
    summaries<-areadist(pointfile,polyfile,maxdist)
    summaries<-data.frame(summaries)
    write.csv(summaries, file=paste0("0_data/1_processed/0_human_footprint/seismic areadist 2025/areadist_",j,".csv"))
    #Combine the separate files.
  }
  else if(nrow(pointfile)==0){
    print(paste0("No points present and no areas of/distances to features calculated within cell ",j))
  }
}


#Now for each footprint, combine the smaller files into one large file
#for each footprint type's area, weighted mean age, and nearest distance
#to each point count survey. There should be 83439 observations in each 
#footprint file.

#seismic lines
names1<-list.files(path = "0_data/1_processed/0_human_footprint/seismic areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/seismic areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/seismic-areadist-allsurveys.csv")


#roads
names1<-list.files(path = "0_data/1_processed/0_human_footprint/roads areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/roads areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/roads-areadist-allsurveys.csv")


#active wells
names1<-list.files(path = "0_data/1_processed/0_human_footprint/active wells areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/active wells areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/active-wells-areadist-allsurveys.csv")

#abandoned wells
names1<-list.files(path = "0_data/1_processed/0_human_footprint/abandoned wells areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/abandoned wells areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/abandoned-wells-areadist-allsurveys.csv")

#pipelines
names1<-list.files(path = "0_data/1_processed/0_human_footprint/pipelines areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/pipelines areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/pipelines-areadist-allsurveys.csv")

#powerlines
names1<-list.files(path = "0_data/1_processed/0_human_footprint/powerlines areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/powerlines areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/powerlines-areadist-allsurveys.csv")

#mines
names1<-list.files(path = "0_data/1_processed/0_human_footprint/mines areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/mines areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/mines-areadist-allsurveys.csv")

#industrial
names1<-list.files(path = "0_data/1_processed/0_human_footprint/industrial areadist 2025/")
modelresults<-list()#empty list each time we draw a new number of samples
for (i in names1){
  #read in model coefficients from best model for each species
  try(tdf<-read.csv(paste0("0_data/1_processed/0_human_footprint/industrial areadist 2025/",i),header=TRUE))
  #add to model list
  modelresults[[i]]<-tdf#append temporary data frame at each loop iteration to the list
  print(paste0("Data frame ",i," added to list"))
  
}
drawnresults = do.call(bind_rows, modelresults)#bind data frames together
write.csv(drawnresults, file="0_data/1_processed/0_human_footprint/industrial-areadist-allsurveys.csv")



#merge files (code assumes same order of location and year among all files)
dense.linear<-read.csv("0_data/1_processed/0_human_footprint/human footprint/dense-linear-areadist.age-allsurveys.csv", header=TRUE)
dense.linear$P150.denseLinear<-dense.linear$PROP150
dense.linear$P565.denseLinear<-dense.linear$PROP565
dense.linear$DIST.denseLinear<-dense.linear$NEAR.DIST
str(dense.linear)
m1<-dense.linear%>%
  select(locationC,year,P150.denseLinear,P565.denseLinear,DIST.denseLinear)


road<-read.csv("0_data/1_processed/0_human_footprint/human footprint/road-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.Road<-road$PROP150
m1$P565.Road<-road$PROP565
m1$DIST.Road<-road$NEAR.DIST


well.high<-read.csv("0_data/1_processed/0_human_footprint/human footprint/well-high-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.wellHigh<-well.high$PROP150
m1$P565.wellHigh<-well.high$PROP565
m1$DIST.wellHigh<-well.high$NEAR.DIST


well.low<-read.csv("0_data/1_processed/0_human_footprint/human footprint/well-low-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.wellLow<-well.low$PROP150
m1$P565.wellLow<-well.low$PROP565
m1$DIST.wellLow<-well.low$NEAR.DIST


plants.mines<-read.csv("0_data/1_processed/0_human_footprint/human footprint/plant-mine-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.plantMines<-plants.mines$PROP150
m1$P565.plantMines<-plants.mines$PROP565
m1$DIST.plantMines<-plants.mines$NEAR.DIST

pipelines.powerlines<-read.csv("0_data/1_processed/0_human_footprint/human footprint/pipeline-powerline-areadist.age-allsurveys.csv", header=TRUE)
m1$P150.wide.linear<-pipelines.powerlines$PROP150
m1$P565.wide.linear<-pipelines.powerlines$PROP565
m1$DIST.wide.linear<-pipelines.powerlines$NEAR.DIST

write.csv(m1, file="0_data/1_processed/0_human_footprint/0_all-human-footprint.csv")

