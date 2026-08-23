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
library(terra)
theme_set(theme_minimal())

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(margin=margin(16,0,0,0)),
        axis.title.y=element_text(margin=margin(0,16,0,0)),
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
levels(as.factor(finalYERAsitesShapefile$clusterID))


obs<-read.csv("0_data/0_raw/point counts/YERA_2025_COLEXT_DATA_DETECTIONHISTORYcombined.csv", header=TRUE)
#manually went through to combine detection histories stored on multiple rows
#for same sites, condensed detection histories to single row per site.

obsFilt<-obs%>%
  select("yera2016_1","yera2016_2","yera2016_3","yera2016_4","yera2016_5",
         "yera2017_1","yera2017_2","yera2017_3","yera2017_4", 
         "yera2017_5","yera2018_1","yera2018_2","yera2018_3","yera2018_4", 
         "yera2018_5","yera2019_1","yera2019_2","yera2019_3","yera2019_4", 
         "yera2019_5","yera2021_1","yera2021_2","yera2021_3","yera2021_4", 
         "yera2021_5")%>%
  mutate(numVisits=rowSums(!is.na(.)))%>%
  mutate(locationC=obs$locationC)

m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
str(m1)#274 sites

#both mining and non-mining regions, but limit to sites
#with at least 3 years of visits
m2<-m1%>%
  filter(numVisits>10)%>%
  arrange(locationC)%>%
  unique()#at least 3 years of visits
str(m2)
table(m2$clusterID)#123 sites
#2_10  2_8  3_1  3_6  3_8  4_5  7_5  7_8 
#   4   15   12    6    7   36   29   14 

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

library(tidyterra)
quickmap<-ggplot() + 
  geom_spatraster(data = Hedley2019)+
  geom_sf(data = initFilt.rp, 
          #shape=22, #square symbol instead of solid circle(default)
          size = 1, 
          aes(col = clusterID))
#Gipsy Gordon points are not covered by Hedley layer

#I think I ended up deciding not to do this analysis
habSuitability<-terra::extract(Hedley2019, initFilt.rp)
hist(habSuitability$YERA_COLL_SUM_FINAL11)
table(habSuitability$YERA_COLL_SUM_FINAL11)
#152 sites with 3+ years of data
#overlapping Hedley layer

#Value:  2  3  4  5  6  7  8 
#nsite:  1  5  3 18 28 30 30 

initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11
initFilt.rp$HabSuit4Class<-ifelse(initFilt.rp$HabSuitOrig<6, "A", ifelse(initFilt.rp$HabSuitOrig<7, "B", ifelse(initFilt.rp$HabSuitOrig<8, "C", "D")))
#turns it to charactor
initFilt.rp$HabSuit4Class<-as.factor(initFilt.rp$HabSuit4Class)
table(initFilt.rp$HabSuit4Class)
# A  B  C  D 
#27 28 30 30 

#make habitat quality sampling frame
# visited2013<-initFilt.rp[,c("yera2013_1",
#                             "yera2013_2",
#                             "yera2013_3",
#                             "yera2013_4",
#                             "yera2013_5")]
# visited2013$geometry<-NULL
# initFilt.rp$visited2013<-rowMeans(visited2013, na.rm=TRUE)
# initFilt.rp$visited2013<-ifelse(initFilt.rp$visited2013=="NaN",NA,initFilt.rp$visited2013)
# initFilt.rp$hqsf2013<-ifelse(!is.na(initFilt.rp$visited2013),initFilt.rp$HabSuit4Class,NA)
# 
# visited2014<-initFilt.rp[,c("yera2014_1",
#                             "yera2014_2",
#                             "yera2014_3",
#                             "yera2014_4",
#                             "yera2014_5")]
# visited2014$geometry<-NULL
# initFilt.rp$visited2014<-rowMeans(visited2014, na.rm=TRUE)
# initFilt.rp$visited2014<-ifelse(initFilt.rp$visited2014=="NaN",NA,initFilt.rp$visited2014)
# initFilt.rp$hqsf2014<-ifelse(!is.na(initFilt.rp$visited2014),initFilt.rp$HabSuit4Class,NA)
# 
# visited2015<-initFilt.rp[,c("yera2015_1",
#                             "yera2015_2",
#                             "yera2015_3",
#                             "yera2015_4",
#                             "yera2015_5")]
# visited2015$geometry<-NULL
# initFilt.rp$visited2015<-rowMeans(visited2015, na.rm=TRUE)
# initFilt.rp$visited2015<-ifelse(initFilt.rp$visited2015=="NaN",NA,initFilt.rp$visited2015)
# initFilt.rp$hqsf2015<-ifelse(!is.na(initFilt.rp$visited2015),initFilt.rp$HabSuit4Class,NA)

visited2016<-initFilt.rp[,c("yera2016_1",
                            "yera2016_2",
                            "yera2016_3",
                            "yera2016_4",
                            "yera2016_5")]
visited2016$geometry<-NULL
initFilt.rp$visited2016<-rowMeans(visited2016, na.rm=TRUE)
initFilt.rp$visited2016<-ifelse(initFilt.rp$visited2016=="NaN",NA,initFilt.rp$visited2016)
initFilt.rp$hqsf2016<-ifelse(!is.na(initFilt.rp$visited2016),initFilt.rp$HabSuit4Class,NA)

visited2017<-initFilt.rp[,c("yera2017_1",
                            "yera2017_2",
                            "yera2017_3",
                            "yera2017_4",
                            "yera2017_5")]
visited2017$geometry<-NULL
initFilt.rp$visited2017<-rowMeans(visited2017, na.rm=TRUE)
initFilt.rp$visited2017<-ifelse(initFilt.rp$visited2017=="NaN",NA,initFilt.rp$visited2017)
initFilt.rp$hqsf2017<-ifelse(!is.na(initFilt.rp$visited2017),initFilt.rp$HabSuit4Class,NA)

visited2018<-initFilt.rp[,c("yera2018_1",
                            "yera2018_2",
                            "yera2018_3",
                            "yera2018_4",
                            "yera2018_5")]
visited2018$geometry<-NULL
initFilt.rp$visited2018<-rowMeans(visited2018, na.rm=TRUE)
initFilt.rp$visited2018<-ifelse(initFilt.rp$visited2018=="NaN",NA,initFilt.rp$visited2018)
initFilt.rp$hqsf2018<-ifelse(!is.na(initFilt.rp$visited2018),initFilt.rp$HabSuit4Class,NA)

visited2019<-initFilt.rp[,c("yera2019_1",
                            "yera2019_2",
                            "yera2019_3",
                            "yera2019_4",
                            "yera2019_5")]
visited2019$geometry<-NULL
initFilt.rp$visited2019<-rowMeans(visited2019, na.rm=TRUE)
initFilt.rp$visited2019<-ifelse(initFilt.rp$visited2019=="NaN",NA,initFilt.rp$visited2019)
initFilt.rp$hqsf2019<-ifelse(!is.na(initFilt.rp$visited2019),initFilt.rp$HabSuit4Class,NA)

visited2021<-initFilt.rp[,c("yera2021_1",
                            "yera2021_2",
                            "yera2021_3",
                            "yera2021_4",
                            "yera2021_5")]
visited2021$geometry<-NULL
initFilt.rp$visited2021<-rowMeans(visited2021, na.rm=TRUE)
initFilt.rp$visited2021<-ifelse(initFilt.rp$visited2021=="NaN",NA,initFilt.rp$visited2021)
initFilt.rp$hqsf2021<-ifelse(!is.na(initFilt.rp$visited2021),initFilt.rp$HabSuit4Class,NA)

# visited2022<-initFilt.rp[,c("yera2022_1",
#                             "yera2022_2",
#                             "yera2022_3",
#                             "yera2022_4",
#                             "yera2022_5")]
# visited2022$geometry<-NULL
# initFilt.rp$visited2022<-rowMeans(visited2022, na.rm=TRUE)
# initFilt.rp$visited2022<-ifelse(initFilt.rp$visited2022=="NaN",NA,initFilt.rp$visited2022)
# initFilt.rp$hqsf2022<-ifelse(!is.na(initFilt.rp$visited2022),initFilt.rp$HabSuit4Class,NA)
# 
# visited2023<-initFilt.rp[,c("yera2023_1",
#                             "yera2023_2",
#                             "yera2023_3",
#                             "yera2023_4",
#                             "yera2023_5")]
# visited2023$geometry<-NULL
# initFilt.rp$visited2023<-rowMeans(visited2023, na.rm=TRUE)
# initFilt.rp$visited2023<-ifelse(initFilt.rp$visited2023=="NaN",NA,initFilt.rp$visited2023)
# initFilt.rp$hqsf2023<-ifelse(!is.na(initFilt.rp$visited2023),initFilt.rp$HabSuit4Class,NA)

#plot(initFilt.rp["hqsf2013"])
#plot(initFilt.rp["hqsf2014"])
#plot(initFilt.rp["hqsf2015"])
plot(initFilt.rp["hqsf2016"])
plot(initFilt.rp["hqsf2017"])
plot(initFilt.rp["hqsf2018"])
plot(initFilt.rp["hqsf2019"])
plot(initFilt.rp["hqsf2021"])
#plot(initFilt.rp["hqsf2022"])
#plot(initFilt.rp["hqsf2023"])

initFilt.rp$geometry<-NULL
#Number of sites with a Hedley rank
ranked<-initFilt.rp%>%
  select(clusterID,locationC,
         hqsf2016,hqsf2017,hqsf2018,hqsf2019,hqsf2021)%>%
  pivot_longer(cols=hqsf2016:hqsf2021, names_to="year", values_to="rating")%>%
  filter(!is.na(rating))%>%
  mutate(year=substr(year,5,8))
nlevels(as.factor(ranked$locationC))#115 sites

#Plotting
hqsf<-initFilt.rp%>%
  select(clusterID,locationC,
         hqsf2016,hqsf2017,hqsf2018,hqsf2019,hqsf2021)%>%
  pivot_longer(cols=hqsf2016:hqsf2021, names_to="year", values_to="rating")%>%
  filter(!is.na(rating))%>%
  mutate(year=substr(year,5,8))%>%
  select(year, rating)%>%
  group_by(year)%>%
  count(rating)%>%
  mutate(ratingLetter=case_when(rating==1 ~ "A",
                                rating==2 ~ "B",
                                rating==3 ~ "C", 
                                rating==4 ~ "D",
                                is.na(rating) ~ NA))
print(hqsf)

hqsf_wide<-hqsf%>%
  pivot_wider(names_from = ratingLetter, values_from = n)
hqsf$Rating<-as.factor(hqsf$ratingLetter)
hqsf$year<-as.integer(hqsf$year)
ratingXtime<-ggplot(hqsf, aes(x = year, y = n, group = Rating, colour=Rating)) +
  geom_line(size = 1) + # Adjust line thickness as needed
  labs(
    title = "",
    x = "Year",
    y = "Number of Sites"
  ) + xlim(2016,2021) +my.theme
  #theme_minimal()

ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/sitetypesvisitedXyear.png",
       plot = ratingXtime,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Do not drop northern sites outside the Hedley et al. layer for now
obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(yera2016_1,yera2016_2,yera2016_3,yera2016_4,yera2016_5,
         yera2017_1,yera2017_2,yera2017_3,yera2017_4,yera2017_5,
         yera2018_1,yera2018_2,yera2018_3,yera2018_4,yera2018_5,
         yera2019_1,yera2019_2,yera2019_3,yera2019_4,yera2019_5,
         yera2021_1,yera2021_2,yera2021_3,yera2021_4,yera2021_5)#,
         #yera2022_1,yera2022_2,yera2022_3,yera2022_4,yera2022_5,
         #yera2023_1,yera2023_2,yera2023_3,yera2023_4,yera2023_5)
obsFiltFinal$geometry<-NULL
obsFiltFinal<-as.matrix(obsFiltFinal)#115 sites

siteCovFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID, clusterIDfinal, HabSuit4Class)%>%
  mutate(HabSuit4ClassNumeric=as.numeric(HabSuit4Class))
siteCovFinal$geometry<-NULL

HabSuit4ClassNumeric<-siteCovFinal$HabSuit4ClassNumeric
yearlySiteCovFinal<-as.matrix(cbind(HabSuit4ClassNumeric,
                                    HabSuit4ClassNumeric,
                                    HabSuit4ClassNumeric,
                                    HabSuit4ClassNumeric,
                                    HabSuit4ClassNumeric))

         
umf<-unmarkedMultFrame(y=obsFiltFinal,
                       siteCovs = siteCovFinal,
                       obsCovs=list(time=matrix(4,115,25)),#5 years*5 rounds
                       yearlySiteCovs=list(habSuit=yearlySiteCovFinal),
                       numPrimary=5) 

umf@siteCovs$HabSuit4ClassNumeric<-as.factor(umf@siteCovs$HabSuit4ClassNumeric)
umf@yearlySiteCovs$habSuit<-as.factor(umf@yearlySiteCovs$habSuit)

model1<-colext(psiformula= ~HabSuit4ClassNumeric, 
               gammaformula =  ~ habSuit, 
               epsilonformula = ~ habSuit,
               pformula = ~ 1, umf)
summary(model1)
#Including a y-intercept
#Initial (logit-scale):
#                      Estimate    SE      z P(>|z|)
#(Intercept)             0.3555 0.424  0.838   0.402
#HabSuit4ClassNumeric2  -1.0269 0.631 -1.628   0.104
#HabSuit4ClassNumeric3  -0.0836 0.646 -0.129   0.897
#HabSuit4ClassNumeric4   1.2562 0.786  1.599   0.110

#Colonization (logit-scale):
#  Estimate    SE      z  P(>|z|)
#(Intercept)  -0.9749 0.366 -2.661 0.00779
#habSuit2      0.1890 0.491  0.385 0.70047
#habSuit3     -0.0836 0.513 -0.163 0.87049
#habSuit4      0.3667 0.523  0.701 0.48322

#Extinction (logit-scale):
#  Estimate    SE      z  P(>|z|)
#(Intercept)   -0.381 0.343 -1.110  0.2671
#habSuit2       1.035 0.576  1.796  0.0725
#habSuit3       1.235 0.529  2.333  0.0196
#habSuit4      -0.236 0.470 -0.502  0.6160

#Detection (logit-scale):
#Estimate     SE    z  P(>|z|)
#   0.478 0.0722 6.62 3.5e-11

#AIC: 1682.768 
#Number of sites: 115
#In contrast to when more years are used
#Initial occupancy & colonization do not vary with Hedley rank
#B & C have higher extinction probability than A, probably because
#Yellow Rails are never in A to begin with

#Setting y-intercept to zero
model1b<-colext(psiformula= ~HabSuit4ClassNumeric-1, 
               gammaformula =  ~ habSuit-1, 
               epsilonformula = ~ habSuit-1,
               pformula = ~ 1, umf)
summary(model1b)
#Initial (logit-scale):
#                      Estimate    SE      z P(>|z|)
#HabSuit4ClassNumeric1    0.355 0.424  0.837  0.4024
#HabSuit4ClassNumeric2   -0.671 0.467 -1.437  0.1508
#HabSuit4ClassNumeric3    0.272 0.487  0.559  0.5759
#HabSuit4ClassNumeric4    1.611 0.661  2.436  0.0148

#Colonization (logit-scale):
#         Estimate    SE     z  P(>|z|)
#habSuit1   -0.975 0.366 -2.66 0.00777
#habSuit2   -0.786 0.327 -2.40 0.01633
#habSuit3   -1.059 0.359 -2.95 0.00317
#habSuit4   -0.608 0.373 -1.63 0.10339

#Extinction (logit-scale):
#  Estimate    SE      z  P(>|z|)
#habSuit1   -0.381 0.343 -1.11  0.2669
#habSuit2    0.653 0.463  1.41  0.1579
#habSuit3    0.854 0.403  2.12  0.0341
#habSuit4   -0.617 0.321 -1.92  0.0550

#Detection (logit-scale):
#Estimate     SE    z  P(>|z|)
#   0.478 0.0722 6.62 3.52e-11

#AIC: 1682.768 
#Number of sites: 115
#Sites ranked as an 8 by Hedley initially more likely to be occupied.
#YERA slightly less likely to go extinct at colonized sites ranked as an 8 by Hedley 
#YERA more likely to go extinct at colonized sites ranked as an 7 by Hedley 
#YERA less likely to colonize sites ranked as an 7 or less

mb.gof.test(model1, nsim = 1000, plot.hist = TRUE,
            report = NULL, parallel = TRUE,
            ncores=2, cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, plot.seasons = FALSE)
#Goodness-of-fit for dynamic occu. model

#Number of seasons:  5 

#Chi-square statistic:
#Season 1 Season 2 Season 3 Season 4 Season 5  
# 87.0703 104.1541  35.4833  31.4974  61.1523 

#Total chi-square = 319.3575 
#Number of bootstrap samples = 1000
#P-value = 0
#Quantiles of bootstrapped statistics:
#  0%  25%  50%  75% 100% 
#110  156  171  188  274 
#Estimate of c-hat = 1.84  

#Get estimates of uncertainty in yearly occu. probability via 
#nonparametric bootstrapping
#use bootstrap to get estimates of psi (occu.) from derived parameter estimates
bootModel100x<-nonparboot(model1, B=100)
str(bootModel100x@projected)
#num [1:2, 1:5, 1:115] 0.432 0.568 0.719 0.281 0.731 ...
#this is a 3d array whose number is equal to 115(#sites)*5(#years)*2(1st vector for prob. absence,
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
#in the colext unmarked Frame in each of the 5 years
#Yellow Rail predicted occu. probability has declined at this site over time.
bootModel100x@projected[1,,1]
#this is how you call up predicted probability that site is unoccupied based on the original data for site 1
#in the colext unmarked Frame
#These results are equivalent to 1 minus bootModel100x@projected[2,,1]

bootModel100x@projected[2,,9]#this is how you call up predicted probability of occu. based on the bootstrap samples for site 9
#in the colext unmarked Frame in each of the 5 years
bootModel100x@projected[2,,]#this is how you call up predicted probability of occu. based on the bootstrap samples 
#in the colext unmarked Frame, for all sites in each of the 5 years. 
#a 208 x 7 matrix
smoothed(bootModel100x)[2,] 
#average occu. probability in each year across all 115 sites
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
umf@siteCovs%>%
  select(locationC,clusterID)
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


#Cluster 3
umf@siteCovs%>%
  select(locationC,clusterID)
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

#Cluster 4: no sites within Hedley layer
umf@siteCovs%>%
  select(locationC,clusterID)
cluster2_10List<-list()#Gipsy Lake
for (j in c(47:49)){
  siteRownum<-j
  m.cluster2_10List<-bootModel100x@projected[2,,j]
  m.cluster2_10List_l95<-bootModel100x@projected[2,,j]-1.96*bootModel100x@projected.mean.bsse[2,]
  m.cluster2_10List_u95<-bootModel100x@projected[2,,j]+1.96*bootModel100x@projected.mean.bsse[2,]
  cluster2_10List[[j]]<-data.frame(site=siteRownum,
                                   occ=m.cluster2_10List,
                                   l95=m.cluster2_10List_l95,
                                   u95=m.cluster2_10List_u95,
                                   year=c(2013,2014,2015,2016,2017,2018,2019,2021))
}
cluster2_10<-do.call(rbind, cluster2_10List)

Occucluster2_10 <- ggplot(data=cluster2_10, aes(x=year, y=occ, group=site)) +
  geom_ribbon(aes(x=year, ymin=l95, ymax=u95, group=site), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=occ, group=site), colour="blue")+
  my.theme+ggtitle("Cluster 4")+xlab("Year")+ylab("Predicted occu.")+
  #my.theme+ggtitle("Cluster2_10 (Gipsy Lake)")+xlab("Year")+ylab("Predicted occu.")+
  ylim(-0.25,1.25)+xlim(2013,2021)


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
umf@siteCovs%>%
  select(locationC,clusterID)
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
ggsave(filename="2_outputs/1_unmarked_outputs/February 2026/occu.ByRegion2016_2021_JustHedleyRaster.png", plot=multipanelPlot, width=16, height=12, units="in", dpi=300)


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
R <- nlevels(as.factor(detectionsLong$locationC))#115 sites
T <- nlevels(as.factor(detectionsLong$round))#up to 5 visits per site per year
S <- nlevels(as.factor(detectionsLong$year))#5 years (no visits to any sites in 2020 due to COVID)

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
# Specify model in JAGS language
model_code <-
  "model
{
  # Specify priors
  alpha.psi ~ dunif(-3, 3) 
  beta1.psi ~ dunif(-3, 3)
  beta2.psi ~ dunif(-3, 3)
  beta3.psi ~ dunif(-3, 3) 

  alpha.gamma ~ dunif(-3, 3)
  beta1.gamma ~ dunif(-3, 3) 
  beta2.gamma ~ dunif(-3, 3)
  beta3.gamma ~ dunif(-3, 3)

  alpha.phi ~ dunif(-3, 3)  
  beta1.phi ~ dunif(-3, 3)
  beta2.phi ~ dunif(-3, 3)
  beta3.phi ~ dunif(-3, 3)

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
    logit(psi1[i])<-alpha.psi +beta1.psi*class2[i] +beta2.psi*class3[i] +beta3.psi*class4[i] + u_psi[cluster[i]]

    z[i,1]<-z0[i]
    for (k in 2:nyear){
      #survival probability from season k-1 to season k
      phi[i,k-1] <-psurv[i,k-1]
      logit(psurv[i,k-1])<-alpha.phi + u_phi[cluster[i]] +beta1.phi*class2SY[i,k-1] +beta2.phi*class3SY[i,k-1] +beta3.phi*class4SY[i,k-1] 

      #colonization probability from season k-1 to season k
      gamma[i,k-1] <-pcol[i,k-1]
      logit(pcol[i,k-1])<-alpha.gamma  +u_gamma[cluster[i]] +beta1.gamma*class2SY[i,k-1] +beta2.gamma*class3SY[i,k-1] +beta3.gamma*class4SY[i,k-1] 

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
                   class2=class2,
                   class3=class3,
                   class4=class4,
                   #site-year covariates (gamma, phi)
                   class2SY=class2SY,
                   class3SY=class3SY,
                   class4SY=class4SY,
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
initsP[[1]]$u_psi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.phi=runif(1, -3, 3)
initsP[[1]]$beta1.phi =runif(1, -3, 3)
initsP[[1]]$beta2.phi =runif(1, -3, 3)
initsP[[1]]$beta3.phi =runif(1, -3, 3)
initsP[[1]]$u_phi=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.gamma=runif(1, -3, 3)
initsP[[1]]$beta1.gamma =runif(1, -3, 3)
initsP[[1]]$beta2.gamma =runif(1, -3, 3)
initsP[[1]]$beta3.gamma =runif(1, -3, 3)
initsP[[1]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[1]]$alpha.p=runif(1,-3,3)
  
initsP[[2]]$z=zst
initsP[[2]]$alpha.psi=runif(1, -3, 3)
initsP[[2]]$beta1.psi=runif(1, -3, 3)
initsP[[2]]$beta2.psi=runif(1, -3, 3)
initsP[[2]]$beta3.psi=runif(1, -3, 3)
initsP[[2]]$u_psi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.phi=runif(1, -3, 3)
initsP[[2]]$beta1.phi =runif(1, -3, 3)
initsP[[2]]$beta2.phi =runif(1, -3, 3)
initsP[[2]]$beta3.phi =runif(1, -3, 3)
initsP[[2]]$u_phi=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.gamma=runif(1, -3, 3)
initsP[[2]]$beta1.gamma =runif(1, -3, 3)
initsP[[2]]$beta2.gamma =runif(1, -3, 3)
initsP[[2]]$beta3.gamma =runif(1, -3, 3)
initsP[[2]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[2]]$alpha.p=runif(1,-3,3)
  
initsP[[3]]$z=zst
initsP[[3]]$alpha.psi=runif(1, -3, 3)
initsP[[3]]$beta1.psi=runif(1, -3, 3)
initsP[[3]]$beta2.psi=runif(1, -3, 3)
initsP[[3]]$beta3.psi=runif(1, -3, 3)
initsP[[3]]$u_psi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.phi=runif(1, -3, 3)
initsP[[3]]$beta1.phi =runif(1, -3, 3)
initsP[[3]]$beta2.phi =runif(1, -3, 3)
initsP[[3]]$beta3.phi =runif(1, -3, 3)
initsP[[3]]$u_phi=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.gamma=runif(1, -3, 3)
initsP[[3]]$beta1.gamma =runif(1, -3, 3)
initsP[[3]]$beta2.gamma =runif(1, -3, 3)
initsP[[3]]$beta3.gamma =runif(1, -3, 3)
initsP[[3]]$u_gamma=rnorm(7, 0, 5)
  
initsP[[3]]$alpha.p=runif(1,-3,3)
  
# Parameters monitored
params <-c("psi",
             "alpha.psi",
             "beta1.psi",
             "beta2.psi",
             "beta3.psi",
             "u_psi",
             "phi",
             "alpha.phi",
             "beta1.phi",
             "beta2.phi",
             "beta3.phi",
             "u_phi",
             "gamma",
             "alpha.gamma",
             "beta1.gamma",
             "beta2.gamma",
             "beta3.gamma",
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
       file = "2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/JAGS_2016_2021_HedleyRaster.Rdata")
  
load("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/JAGS_2016_2021_HedleyRaster.Rdata")
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
gelman.diag(posterior_draws[1:3][,"alpha.phi"])#[1,]          1          1       
gelman.diag(posterior_draws[1:3][,"beta1.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta2.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta3.phi"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.gamma"])#[1,]          1          1       
gelman.diag(posterior_draws[1:3][,"beta1.gamma"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta2.gamma"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"beta3.gamma"])#[1,]          1          1           
gelman.diag(posterior_draws[1:3][,"alpha.p"])#[1,]          1          1           

range(
  effectiveSize(posterior_draws[1:3][,"alpha.psi"]),        
  effectiveSize(posterior_draws[1:3][,"beta1.psi"]),           
  effectiveSize(posterior_draws[1:3][,"beta2.psi"]),           
  effectiveSize(posterior_draws[1:3][,"beta3.psi"]),           
  effectiveSize(posterior_draws[1:3][,"alpha.phi"]),        
  effectiveSize(posterior_draws[1:3][,"beta1.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta2.phi"]),           
  effectiveSize(posterior_draws[1:3][,"beta3.phi"]),           
  effectiveSize(posterior_draws[1:3][,"alpha.gamma"]),        
  effectiveSize(posterior_draws[1:3][,"beta1.gamma"]),           
  effectiveSize(posterior_draws[1:3][,"beta2.gamma"]),           
  effectiveSize(posterior_draws[1:3][,"beta3.gamma"])           
)
#effective sample size varies from 3931.714 13650.778

heidel.diag(posterior_draws[1:3][,"alpha.psi"])#passed both        
heidel.diag(posterior_draws[1:3][,"beta1.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta2.psi"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"beta3.psi"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.phi"])#passed first, failed second        
heidel.diag(posterior_draws[1:3][,"beta1.phi"])#passed both          
heidel.diag(posterior_draws[1:3][,"beta2.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta3.phi"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.gamma"])#passed both        
heidel.diag(posterior_draws[1:3][,"beta1.gamma"])#passed both           
heidel.diag(posterior_draws[1:3][,"beta2.gamma"])#passed first, failed second           
heidel.diag(posterior_draws[1:3][,"beta3.gamma"])#passed both           
heidel.diag(posterior_draws[1:3][,"alpha.p"])#passed both           
#All passed the stationary test (p>0.1), failed the halfwidth test

gd <- gelman.diag(posterior_draws[1:3][,"alpha.psi"])
es <- effectiveSize(posterior_draws[1:3][,"alpha.psi"])#   
hd <- heidel.diag(posterior_draws[1:3][,"alpha.psi"])#

#Trace Plots and Density Plots
plot(posterior_draws[1:3][,"beta1.psi"],density=FALSE)#nonsig (mostly -) effect on psi
densplot(posterior_draws[[1]][,"beta1.psi"], main="beta1.psi (habitat category 2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta2.psi"],density=FALSE)#nonsig effect on psi
densplot(posterior_draws[[1]][,"beta2.psi"], main="beta2.psi (habitat category 3) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta3.psi"],density=FALSE)#mostly + effect on psi
densplot(posterior_draws[[1]][,"beta3.psi"], main="beta3.psi (habitat category 4) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta1.phi"],density=FALSE)#nonsig (mostly -) effect on phi
densplot(posterior_draws[[1]][,"beta1.phi"], main="beta1.phi (habitat category 2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta2.phi"],density=FALSE)#negative effect on phi
densplot(posterior_draws[[1]][,"beta2.phi"], main="beta2.phi (habitat category 3) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta3.phi"],density=FALSE)#nonsig (more +) effect on phi
densplot(posterior_draws[[1]][,"beta3.phi"], main="beta3.phi (habitat category 4) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta1.gamma"],density=FALSE)#nonsig effect on gamma
densplot(posterior_draws[[1]][,"beta1.gamma"], main="beta1.gamma (habitat category 2) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta2.gamma"],density=FALSE)#nonsig (mostly +) effect on gamma
densplot(posterior_draws[[1]][,"beta2.gamma"], main="beta2.gamma (habitat category 3) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)

plot(posterior_draws[1:3][,"beta3.gamma"],density=FALSE)#nonsig (mostly +) effect on gamma
densplot(posterior_draws[[1]][,"beta3.gamma"], main="beta3.gamma (habitat category 4) Density Plot", xlab="Value", ylab="Density", col="red", lwd=2)


#create ggridges plots
# Combine the results in a single data frame: gamma
post_sim_gamma <- data.frame(samples = c(posterior_draws[[1]][,"beta1.gamma"],posterior_draws[[2]][,"beta1.gamma"],posterior_draws[[3]][,"beta1.gamma"],
                                         posterior_draws[[1]][,"beta2.gamma"],posterior_draws[[2]][,"beta2.gamma"],posterior_draws[[3]][,"beta2.gamma"],
                                         posterior_draws[[1]][,"beta3.gamma"],posterior_draws[[2]][,"beta3.gamma"],posterior_draws[[3]][,"beta3.gamma"]),
                             posteriors = rep(c("Class 2 (Hedley 6)","Class 3 (Hedley 7)","Class 4 (Hedley 8)"), each = 18000),
                             submodel = rep(c("gamma","gamma","gamma"), each = 18000))

ggplot(post_sim_gamma, aes(x = samples, fill = posteriors)) + 
  geom_density(alpha = 0.5) + theme_ridges()

# Closed polygons
p0<-ggplot(post_sim_gamma, aes(x = samples, y = posteriors, group = posteriors)) + 
  geom_density_ridges2(fill = "#00AFBB")+xlim(-1.25,2)
ggsave(filename="2_outputs/2_jags_outputs_feb2026/gammaHedleyClass.png",
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
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/gammaHedleyClass_Coefsridgeplot.png",
       plot = ggr_gamma,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

# Combine the results in a single data frame: phi
post_sim_phi <- data.frame(samples = c(posterior_draws[[1]][,"beta1.phi"],posterior_draws[[2]][,"beta1.phi"],posterior_draws[[3]][,"beta1.phi"],
                                       posterior_draws[[1]][,"beta2.phi"],posterior_draws[[2]][,"beta2.phi"],posterior_draws[[3]][,"beta2.phi"],
                                       posterior_draws[[1]][,"beta3.phi"],posterior_draws[[2]][,"beta3.phi"],posterior_draws[[3]][,"beta3.phi"]),
                           posteriors = rep(c("Class 2 (Hedley 6)","Class 3 (Hedley 7)","Class 4 (Hedley 8)"), each = 18000),
                           submodel = rep(c("phi","phi","phi"), each = 18000))

ggplot(post_sim_phi, aes(x = samples, fill = posteriors)) + 
  geom_density(alpha = 0.5) + theme_ridges()

# Closed polygons
p0<-ggplot(post_sim_phi, aes(x = samples, y = posteriors, group = posteriors)) + 
  geom_density_ridges2(fill = "#00AFBB")+xlim(-3,3)
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/phiHedleyClass.png",
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
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/phiHedleyClass_Coefsridgeplot.png",
       plot = ggr_phi,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


# Combine the results in a single data frame: psi
post_sim_psi <- data.frame(samples = c(posterior_draws[[1]][,"beta1.psi"],posterior_draws[[2]][,"beta1.psi"],posterior_draws[[3]][,"beta1.psi"],
                                       posterior_draws[[1]][,"beta2.psi"],posterior_draws[[2]][,"beta2.psi"],posterior_draws[[3]][,"beta2.psi"],
                                       posterior_draws[[1]][,"beta3.psi"],posterior_draws[[2]][,"beta3.psi"],posterior_draws[[3]][,"beta3.psi"]),
                           posteriors = rep(c("Class 2 (Hedley 6)","Class 3 (Hedley 7)","Class 4 (Hedley 8)"), each = 18000),
                           submodel = rep(c("psi","psi","psi"), each = 18000))

ggplot(post_sim_psi, aes(x = samples, fill = posteriors)) + 
  geom_density(alpha = 0.5) + theme_ridges()

# Closed polygons
p0<-ggplot(post_sim_psi, aes(x = samples, y = posteriors, group = posteriors)) + 
  geom_density_ridges2(fill = "#00AFBB")+xlim(-3,3)
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/psiHedleyClass.png",
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
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/psiHedleyClass_Coefsridgeplot.png",
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
                  ggr_gamma2,
                  ggr_phi2,
                  nrow=1,ncol=3)

ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/allHedleyClass_Coefsridgeplot.png",
       plot = g3,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Make box plots for different Hedley class types A(1-5), B(6), C(7), D(8)
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

#psi: A
psiA.median<-exp(alpha.psi.PostMedian)/(1+exp(alpha.psi.PostMedian))
psiA.2.5<-exp(alpha.psi.Post.2.5)/(1+exp(alpha.psi.Post.2.5))
psiA.97.5<-exp(alpha.psi.Post.97.5)/(1+exp(alpha.psi.Post.97.5))
#psi: B
psiB.median<-exp(alpha.psi.PostMedian+beta1.psi.PostMedian)/(1+exp(alpha.psi.PostMedian+beta1.psi.PostMedian))
psiB.2.5<-exp(alpha.psi.Post.2.5+beta1.psi.Post.2.5)/(1+exp(alpha.psi.Post.2.5+beta1.psi.Post.2.5))
psiB.97.5<-exp(alpha.psi.Post.97.5+beta1.psi.Post.97.5)/(1+exp(alpha.psi.Post.97.5+beta1.psi.Post.97.5))
#psi: C
psiC.median<-exp(alpha.psi.PostMedian+beta2.psi.PostMedian)/(1+exp(alpha.psi.PostMedian+beta2.psi.PostMedian))
psiC.2.5<-exp(alpha.psi.Post.2.5+beta2.psi.Post.2.5)/(1+exp(alpha.psi.Post.2.5+beta2.psi.Post.2.5))
psiC.97.5<-exp(alpha.psi.Post.97.5+beta2.psi.Post.97.5)/(1+exp(alpha.psi.Post.97.5+beta2.psi.Post.97.5))
#psi: D
psiD.median<-exp(alpha.psi.PostMedian+beta3.psi.PostMedian)/(1+exp(alpha.psi.PostMedian+beta3.psi.PostMedian))
psiD.2.5<-exp(alpha.psi.Post.2.5+beta3.psi.Post.2.5)/(1+exp(alpha.psi.Post.2.5+beta3.psi.Post.2.5))
psiD.97.5<-exp(alpha.psi.Post.97.5+beta3.psi.Post.97.5)/(1+exp(alpha.psi.Post.97.5+beta3.psi.Post.97.5))

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

#phi: A
phiA.median<-exp(alpha.phi.PostMedian)/(1+exp(alpha.phi.PostMedian))
phiA.2.5<-exp(alpha.phi.Post.2.5)/(1+exp(alpha.phi.Post.2.5))
phiA.97.5<-exp(alpha.phi.Post.97.5)/(1+exp(alpha.phi.Post.97.5))
#phi: B
phiB.median<-exp(alpha.phi.PostMedian+beta1.phi.PostMedian)/(1+exp(alpha.phi.PostMedian+beta1.phi.PostMedian))
phiB.2.5<-exp(alpha.phi.Post.2.5+beta1.phi.Post.2.5)/(1+exp(alpha.phi.Post.2.5+beta1.phi.Post.2.5))
phiB.97.5<-exp(alpha.phi.Post.97.5+beta1.phi.Post.97.5)/(1+exp(alpha.phi.Post.97.5+beta1.phi.Post.97.5))
#phi: C
phiC.median<-exp(alpha.phi.PostMedian+beta2.phi.PostMedian)/(1+exp(alpha.phi.PostMedian+beta2.phi.PostMedian))
phiC.2.5<-exp(alpha.phi.Post.2.5+beta2.phi.Post.2.5)/(1+exp(alpha.phi.Post.2.5+beta2.phi.Post.2.5))
phiC.97.5<-exp(alpha.phi.Post.97.5+beta2.phi.Post.97.5)/(1+exp(alpha.phi.Post.97.5+beta2.phi.Post.97.5))
#phi: D
phiD.median<-exp(alpha.phi.PostMedian+beta3.phi.PostMedian)/(1+exp(alpha.phi.PostMedian+beta3.phi.PostMedian))
phiD.2.5<-exp(alpha.phi.Post.2.5+beta3.phi.Post.2.5)/(1+exp(alpha.phi.Post.2.5+beta3.phi.Post.2.5))
phiD.97.5<-exp(alpha.phi.Post.97.5+beta3.phi.Post.97.5)/(1+exp(alpha.phi.Post.97.5+beta3.phi.Post.97.5))


alpha.gamma.PostMedian<-smod$quantiles["alpha.gamma","50%"]
alpha.gamma.Post.2.5<-smod$quantiles["alpha.gamma","2.5%"]
alpha.gamma.Post.97.5<-smod$quantiles["alpha.gamma","97.5%"]

beta1.gamma.PostMedian<-smod$quantiles["beta1.gamma","50%"]
beta1.gamma.Post.2.5<-smod$quantiles["beta1.gamma","2.5%"]
beta1.gamma.Post.97.5<-smod$quantiles["beta1.gamma","97.5%"]

beta2.gamma.PostMedian<-smod$quantiles["beta2.gamma","50%"]
beta2.gamma.Post.2.5<-smod$quantiles["beta2.gamma","2.5%"]
beta2.gamma.Post.97.5<-smod$quantiles["beta2.gamma","97.5%"]

beta3.gamma.PostMedian<-smod$quantiles["beta3.gamma","50%"]
beta3.gamma.Post.2.5<-smod$quantiles["beta3.gamma","2.5%"]
beta3.gamma.Post.97.5<-smod$quantiles["beta3.gamma","97.5%"]

#gamma: A
gammaA.median<-exp(alpha.gamma.PostMedian)/(1+exp(alpha.gamma.PostMedian))
gammaA.2.5<-exp(alpha.gamma.Post.2.5)/(1+exp(alpha.gamma.Post.2.5))
gammaA.97.5<-exp(alpha.gamma.Post.97.5)/(1+exp(alpha.gamma.Post.97.5))
#gamma: B
gammaB.median<-exp(alpha.gamma.PostMedian+beta1.gamma.PostMedian)/(1+exp(alpha.gamma.PostMedian+beta1.gamma.PostMedian))
gammaB.2.5<-exp(alpha.gamma.Post.2.5+beta1.gamma.Post.2.5)/(1+exp(alpha.gamma.Post.2.5+beta1.gamma.Post.2.5))
gammaB.97.5<-exp(alpha.gamma.Post.97.5+beta1.gamma.Post.97.5)/(1+exp(alpha.gamma.Post.97.5+beta1.gamma.Post.97.5))
#gamma: C
gammaC.median<-exp(alpha.gamma.PostMedian+beta2.gamma.PostMedian)/(1+exp(alpha.gamma.PostMedian+beta2.gamma.PostMedian))
gammaC.2.5<-exp(alpha.gamma.Post.2.5+beta2.gamma.Post.2.5)/(1+exp(alpha.gamma.Post.2.5+beta2.gamma.Post.2.5))
gammaC.97.5<-exp(alpha.gamma.Post.97.5+beta2.gamma.Post.97.5)/(1+exp(alpha.gamma.Post.97.5+beta2.gamma.Post.97.5))
#gamma: D
gammaD.median<-exp(alpha.gamma.PostMedian+beta3.gamma.PostMedian)/(1+exp(alpha.gamma.PostMedian+beta3.gamma.PostMedian))
gammaD.2.5<-exp(alpha.gamma.Post.2.5+beta3.gamma.Post.2.5)/(1+exp(alpha.gamma.Post.2.5+beta3.gamma.Post.2.5))
gammaD.97.5<-exp(alpha.gamma.Post.97.5+beta3.gamma.Post.97.5)/(1+exp(alpha.gamma.Post.97.5+beta3.gamma.Post.97.5))

data_summary <- data.frame(
  habcat = c("A (Hedley 0-5)", "B (Hedley 6)", "C (Hedley 7)", "D (Hedley 8)"),
  mean_psi = c(psiA.median, psiB.median, psiC.median, psiD.median),
  ci_lower_psi = c(psiA.2.5, psiB.2.5, psiC.2.5, psiD.2.5),
  ci_upper_psi = c(psiA.97.5, psiB.97.5, psiC.97.5, psiD.97.5),
  mean_phi = c(phiA.median, phiB.median, phiC.median, phiD.median),
  ci_lower_phi = c(phiA.2.5, phiB.2.5, phiC.2.5, phiD.2.5),
  ci_upper_phi = c(phiA.97.5, phiB.97.5, phiC.97.5, phiD.97.5),
  mean_gamma = c(gammaA.median, gammaB.median, gammaC.median, gammaD.median),
  ci_lower_gamma = c(gammaA.2.5, gammaB.2.5, gammaC.2.5, gammaD.2.5),
  ci_upper_gamma = c(gammaA.97.5, gammaB.97.5, gammaC.97.5, gammaD.97.5)
)

crossbar_occ<-ggplot(data_summary, aes(x = habcat, y = mean_psi)) +
  geom_crossbar(aes(ymin = ci_lower_psi, ymax = ci_upper_psi), width = 0.5, fill = "lightblue") +
  labs(x = "Habitat Suitability",
       y = "Occupancy") +
  theme_minimal()

crossbar_col<-ggplot(data_summary, aes(x = habcat, y = mean_gamma)) +
  geom_crossbar(aes(ymin = ci_lower_gamma, ymax = ci_upper_gamma), width = 0.5, fill = "lightblue") +
  labs(x = "Habitat Suitability",
       y = "Colonization") +
  theme_minimal()


crossbar_pers<-ggplot(data_summary, aes(x = habcat, y = mean_phi)) +
  geom_crossbar(aes(ymin = ci_lower_phi, ymax = ci_upper_phi), width = 0.5, fill = "lightblue") +
  labs(x = "Habitat Suitability",
       y = "Persistence") +
  theme_minimal()

g3alternate<-grid.arrange(crossbar_occ,
                          crossbar_col,
                          crossbar_pers,
                          nrow=1,ncol=3)

ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/allHedleyClass_Probs_CrossbarPlot.png",
       plot = g3alternate,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

finalbarandhistplots<-grid.arrange(g3,
                          g3alternate,
                          nrow=2,ncol=1)

ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/Hist_Plus_CrossbarPlot.png",
       plot = finalbarandhistplots,
       width = 12,
       height = 16,
       units = c("in"),
       dpi = 300)


#Get WAIC
which(substr(colnames(posterior_draws[[1]]),1,23)=="log_likelihood_contrib[") #5440:8314 #There are 2275 (=65 sites*7 years*up to 5 visits per year)
colnames(posterior_draws[[1]])[5440]#"log_likelihood_contrib[1,1,1]"
posterior_draws[[1]][5440]
traceplot(posterior_draws[1:3][,5440], main="log_likelihood_contrib[1,1,1] Trace Plot", xlab="Iteration", ylab="Value", lwd=1)
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
      try(myArray[i,j,k]<-posterior_draws[[j]][i,(k+5439)])
    }
  }
}
waicMyArray<-waic(myArray)
#Computed from 18000 by 5840 log-likelihood matrix.
#          Estimate   SE
#elpd_waic   -883.7 16.7
#p_waic       103.7  2.8
#waic        1767.3 33.4
looMyArray<-loo(myArray)
#Warning message:
#Some Pareto k diagnostic values are too high. 
#See help('pareto-k-diagnostic') for details.
#Computed from 18000 by 5840 log-likelihood matrix.
#         Estimate   SE
#elpd_loo   -883.7 16.7
#p_loo       103.7  2.8
#looic      1767.4 33.4
#------
#MCSE of elpd_loo is NA.
#MCSE and ESS estimates assume independent draws (r_eff=1).
#Pareto k diagnostic values:
#                        Count   Pct. Min. ESS
#(-Inf, 0.7]   (good)     1755  61.1%   15053   
#(0.7, 1]   (bad)         0   0.0%   <NA>    
#(1, Inf)   (very bad) 1119  38.9%   <NA>    
#See help('pareto-k-diagnostic') for details.
save(myArray, waicMyArray, looMyArray, file="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/WAIC_HedleyModelNorth2013_2021.Rdata")

#Get Estimated Number of Occupied Sites (out of 136)
which(substr(colnames(posterior_draws[[1]]),1,6)=="n.occ[")
# [1] 8315 8316 8317 8318 8319
colnames(posterior_draws[[1]])[8315]#n.occ[1]

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
#[1] 0.7011897
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


ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModel_2016_2021_predOcc.png",
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

ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModel_2016_2021_SitesWithDetections.png",
       plot = p2,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

p3<-grid.arrange(p1,p2,nrow=1,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/NumberOccupiedSitesPlusRaw_2016_2021.png",
       plot = p3,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Get Turnover
which(substr(colnames(posterior_draws[[1]]),1,9)=="turnover[")
#[1]8324 8325 8326 8327
colnames(posterior_draws[[1]])[8324]#turnover[1]

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

ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModel_2016_2021_PredictedTurnover.png",
       plot = p4,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


#Get AUC

which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")
colnames(posterior_draws[[1]])[4405]#"z[1,1]"
colnames(posterior_draws[[1]])[4520]#"z[1,2]"
colnames(posterior_draws[[1]])[4979]#"z[115,5]"

which(substr(colnames(posterior_draws[[1]]),1,4)=="muZ[")
colnames(posterior_draws[[1]])[4980]#"muZ[1,2]"
colnames(posterior_draws[[1]])[5439]#"muZ[115,5]"
str(posterior_draws[[1]][,"muZ[2,2]"])#no estimates for season 1
posterior_draws[[1]][,"muZ[1,2]"]#no estimates for season 1

niter<-ni
nthin<-nt
nchain<-nc
nsite<-115
nyear<-5

library(ggmcmc)
allmuZs<-posterior_draws[1:3][,4980:5439]#
ggd <- ggs(allmuZs)
allZs<-posterior_draws[1:3][,4520:4979]
ggd2 <- ggs(allZs)#starting from year 2

nstore <- niter / nthin
psi.pred <- array(dim=c(nsite, nstore*nchain))
Z.est <- array(dim=c(nsite, nstore*nchain))


#AUC in 2017
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
  psi.vals <- psi.pred[, i]
  Z.vals <- Z.est[, i]
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #A prediction instance with 115 data points
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.5966939
median(AUClist)#[1] 0.5966451
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.5059504 0.5966451 0.6881571 

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~3.60% above 0.7 

#AUC in 2018
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
  psi.vals <- psi.pred[, i]
  Z.vals <- Z.est[, i]
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #A prediction instance with 115 data points
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7254322
median(AUClist)#[1] 0.7271018
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6324821 0.7271018 0.8091076

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~72.34% above 0.7 

#AUC in 2019
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
  psi.vals <- psi.pred[, i]
  Z.vals <- Z.est[, i]
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #A prediction instance with 115 data points
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7716757
median(AUClist)#[1] 0.7743333
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6728505 0.7743333 0.8546674

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~96.68% above 0.7 

#AUC in 2021
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
  psi.vals <- psi.pred[, i]
  Z.vals <- Z.est[, i]
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  #A prediction instance with 115 data points
  perf <- performance(pred, "auc") 
  #A performance instance 'Area under the ROC curve'
  AUC1[i] <- perf@y.values[[1]]
}
AUClist<-do.call(rbind, AUC1)
mean(AUClist)#[1] 0.7347415
median(AUClist)#[1] 0.7363242
quantile(AUClist, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.6389981 0.7363242 0.8229684 

# Specify the value to compare against
value_threshold <- 0.7

# Create the histogram
hist_result <- hist(AUClist, plot = FALSE) # plot=FALSE prevents plotting the histogram
relative_frequencies <- hist_result$counts / sum(hist_result$counts)
indices_above_threshold <- which(hist_result$breaks[-1] > value_threshold)
proportion_above_threshold <- sum(relative_frequencies[indices_above_threshold])
print(proportion_above_threshold)#~76.87% above 0.7 

#Model is adequate

#Predicted yearly occu. for individual clusters
which(substr(colnames(posterior_draws[[1]]),1,4)=="muZ[")
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[")

#575 entries = 115 sites * 5 years
library(ggmcmc)
allZs<-posterior_draws[1:3][,4405:4979]#
ggd <- ggs(allZs)

nsite<-115
nstore <- niter / nthin

#Probabilities of occu. at year 1 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 1, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
meanOccEst2016<-rowMeans(z.pred)

#Probabilities of occu. at year 2 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
meanOccEst2017<-rowMeans(z.pred)


#Probabilities of occu. at year 3 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
meanOccEst2018<-rowMeans(z.pred)

#Probabilities of occu. at year 4 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
meanOccEst2019<-rowMeans(z.pred)

#Probabilities of occu. at year 5 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 5, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
meanOccEst2021<-rowMeans(z.pred)


meanOccEstAllYears<-cbind(meanOccEst2016,
      meanOccEst2017,
      meanOccEst2018,
      meanOccEst2019,
      meanOccEst2021)

#recreate obsFiltFinal dataset but with location and clusterID
#instead of YERA detections
# finalYERAsitesShapefile<-st_read("0_data/0_raw/point counts/allYERAsites_June2025.shp")
# str(finalYERAsitesShapefile)
# #need clusterIDs corrected so that a couple too-large
# #clusters are split up among wetlands
# finalYERAsitesShapefile<-finalYERAsitesShapefile%>%
#   mutate(clusterIDfinal=clusterID)%>%
#   mutate(clusterIDfinal=case_when(locationC %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) ~ "2_8b",
#                                   locationC %in% c(64,65,58) ~ "7_8b",
#                                   locationC %in% c(60,61,62,63) ~ "7_8c",
#                                   locationC %in% c(53,54,55,56) ~ "7_8d",
#                                   !locationC %in% c(1:16,53:56,58,60:65) ~ clusterID))
# levels(as.factor(finalYERAsitesShapefile$clusterIDfinal))
# 
# 
# obs<-read.csv("0_data/0_raw/point counts/YERA_2025_COLEXT_DATA_DETECTIONHISTORYcombined.csv", header=TRUE)
# #manually went through to combine detection histories stored on multiple rows
# #for same sites, condensed detection histories to single row per site.
# 
# #This analysis will use 2013 - 2021 and sites
# #from the whole YERA study area. 2013-2015 are
# #included because we are not using Sentinel data.
# #2021 is the end year because southern sites were
# #not systematically surveyed beyond 2021
# obsFilt<-obs%>%
#   filter(!is.na(HabSuit4ClassNumeric))%>%
#   select("yera2016_1","yera2016_2","yera2016_3","yera2016_4", 
#          "yera2016_5","yera2017_1","yera2017_2","yera2017_3","yera2017_4", 
#          "yera2017_5","yera2018_1","yera2018_2","yera2018_3","yera2018_4", 
#          "yera2018_5","yera2019_1","yera2019_2","yera2019_3","yera2019_4", 
#          "yera2019_5","yera2021_1","yera2021_2","yera2021_3","yera2021_4", 
#          "yera2021_5")%>%
#   mutate(numVisits=rowSums(!is.na(.)))%>%
#   mutate(locationC=obs$locationC)
# 
# #,"yera2023_1","yera2023_2","yera2023_3","yera2023_4", 
# #"yera2023_5"
# 
# m1<-merge(finalYERAsitesShapefile, obsFilt, by=c("locationC"))
# str(m1)
# 
# #both mining and non-mining regions, but limit to sites
# #with at least 3 years of visits
# m2<-m1%>%
#   filter(numVisits>10)%>%
#   arrange(locationC)%>%
#   unique()#at least 3 years of visits
# str(m2)
# table(m2$clusterID)#136 sites
# # 2_10  2_8  3_1  3_6  3_8  4_5  7_5  7_8 
# #4   15   12    6    7   36   29   14 
# 
# table(m2$clusterIDfinal)
# #2_10 2_8b  3_1  3_6  3_8  4_5  7_5  7_8 7_8b 7_8c 7_8d 
# #4   15   12    6    7   36   29    7    1    2    4   
# 
# #get covariates affecting initial occu.
# initFilt<-m2%>%
#   #select(locationC, clusterID)%>%
#   arrange(locationC)
# #160 sites with at least 3 years of data from 2013 - 2021
# 
# #Hedley et al. layer
# Hedley2019<-rast("0_data/1_processed/6_HedleyEtAl2019/YERA_COLL_SUM_FINAL11.tif")
# plot(Hedley2019$YERA_COLL_SUM_FINAL11)
# raster_bbox<-st_bbox(Hedley2019)
# bbox_polygon <- st_as_sfc(raster_bbox)
# bbox_sf<-st_as_sf(bbox_polygon)
# 
# #Reproject points to Hedley layer's CRS
# initFilt.rp<-initFilt%>%
#   st_set_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")%>%
#   st_transform(crs=st_crs(Hedley2019))
# #Extract values from Hedley layer to points
# habSuitability<-terra::extract(Hedley2019, initFilt.rp)
# hist(habSuitability$YERA_COLL_SUM_FINAL11)
# table(habSuitability$YERA_COLL_SUM_FINAL11)
# #115 sites with 3+ years of data
# #overlapping Hedley layer
# 
# #Value:  2  3  4  5  6  7  8 
# #nsite:  1  5  3 18 28 30 30
# 
# initFilt.rp$HabSuitOrig<-habSuitability$YERA_COLL_SUM_FINAL11

obsFiltFinal<-initFilt.rp%>%
  filter(!is.na(HabSuitOrig))%>%
  select(locationC, clusterID)
obsFiltFinal$geometry<-NULL

finalOccDF<-cbind(obsFiltFinal, meanOccEstAllYears)
finalOccDF$clusterID
write.csv(finalOccDF,file="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModel_2016_2021_PredictedOcc.csv")

finalOccDF<-read.csv("2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModel_2016_2021_PredictedOcc.csv", header=TRUE)
levels(as.factor(finalOccDF$clusterID))
#[1] "2_8" "3_1" "3_6" "3_8" "4_5" "7_5" "7_8"

#Cluster 1 (formerly 3_6)
finalOccDF_cluster1<-finalOccDF%>%
  filter(clusterID=="3_6")

finalOccDF_cluster1long<-finalOccDF_cluster1%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster1)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster3_6 <- ggplot(data=finalOccDF_cluster1long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 1: ",nrow(finalOccDF_cluster1)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  #my.theme+ggtitle("Cluster3_6 (Audet Lake)")+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))

#Cluster 2
finalOccDF_cluster2<-finalOccDF%>%
  filter(clusterID=="4_5")

finalOccDF_cluster2long<-finalOccDF_cluster2%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster2)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster4_5 <- ggplot(data=finalOccDF_cluster2long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 2: ",nrow(finalOccDF_cluster2)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


#Cluster 3
finalOccDF_cluster3<-finalOccDF%>%
  filter(clusterID=="7_5")

finalOccDF_cluster3long<-finalOccDF_cluster3%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster3)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster7_5 <- ggplot(data=finalOccDF_cluster3long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 3: ",nrow(finalOccDF_cluster3)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


#Cluster 4 (just 3 sites at GG Provincial Park, not in Hedley layer)
finalOccDF_cluster4<-finalOccDF%>%
  filter(clusterID=="2_1")

finalOccDF_cluster4long<-finalOccDF_cluster4%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster4)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster2_1 <- ggplot(data=finalOccDF_cluster4long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 4: ",nrow(finalOccDF_cluster4)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


#Cluster 4 (would be 5 if I did not restrict to Hedley sites)
finalOccDF_cluster5<-finalOccDF%>%
  filter(clusterID=="3_1")

finalOccDF_cluster5long<-finalOccDF_cluster5%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster5)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster3_1 <- ggplot(data=finalOccDF_cluster5long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 4: ",nrow(finalOccDF_cluster5)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


#Cluster 5 (would be 6 if I did not restrict to Hedley sites)
finalOccDF_cluster6<-finalOccDF%>%
  filter(clusterID=="7_8")

finalOccDF_cluster6long<-finalOccDF_cluster6%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster6)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster7_8 <- ggplot(data=finalOccDF_cluster6long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 5: ",nrow(finalOccDF_cluster6)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


#Cluster 6 (would be 7 if I did not restrict to Hedley sites)
finalOccDF_cluster7<-finalOccDF%>%
  filter(clusterID=="3_8")

finalOccDF_cluster7long<-finalOccDF_cluster7%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster7)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster3_8 <- ggplot(data=finalOccDF_cluster7long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 6: ",nrow(finalOccDF_cluster7)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))


#Cluster 7 (would be 8 if I did not restrict to Hedley sites)
finalOccDF_cluster8<-finalOccDF%>%
  filter(clusterID=="2_8")

finalOccDF_cluster8long<-finalOccDF_cluster8%>%
  pivot_longer(cols=meanOccEst2016:meanOccEst2021,
               names_to = "year",
               values_to = "occ")%>%
  mutate(Year=substr(year,11,14))%>%
  group_by(Year)%>%
  summarize(medYear=median(occ),
            sdYear=sd(occ),
            seYear=sd(occ)/sqrt(nrow(finalOccDF_cluster8)),
            bci2.5=medYear-1.96*seYear,
            bci97.5=medYear+1.96*seYear)%>%
  mutate(Year=as.integer(Year))

Occucluster2_8 <- ggplot(data=finalOccDF_cluster8long, aes(x=Year, y=medYear)) +
  geom_ribbon(aes(x=Year, ymin=bci2.5, ymax=bci97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=Year, y=medYear), colour="blue")+
  my.theme+ggtitle(paste0("Cluster 7: ",nrow(finalOccDF_cluster8)," sites"))+xlab("Year")+ylab("Pred. occu.")+
  xlim(2016,2021)+#ylim(-0.5,1.25)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2013,2021)))))

p7<-grid.arrange(Occucluster3_6,
                 Occucluster4_5,
                 Occucluster7_5,
                 Occucluster3_1,
                 Occucluster7_8,
                 Occucluster3_8,
                 Occucluster2_8,
                 nrow=3,ncol=3)
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModelNorth_MedianPlus95BCIoccupancy_2013_2023.png",
       plot = p7,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)

#Still to do
#Predicted yearly occu. for individual clusters
which(substr(colnames(posterior_draws[[1]]),1,2)=="z[") 
#575 entries = 115 sites * 5 years
library(ggmcmc)
allZs<-posterior_draws[1:3][,6389:7216]#
ggd <- ggs(allZs)

nsite<-115
nstore <- niter / nthin
nchain<-3

# obsFiltFinal<-initFilt.rp%>%
#   filter(!is.na(HabSuitOrig))%>%
#   select(locationC, clusterID, HabSuit4Class)
# obsFiltFinal$geometry<-NULL
#use siteCovFinal

#Probabilities of occu. at year 1 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 1, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2016<-rowMeans(z.pred)
# OccEst2016<-cbind(siteCovFinal, z.pred)
# #I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
# #in columns 3-18002
# OccEst2016$locationC<-NULL 
# OccEst2016$HabSuit4Class<-NULL
# #drop before summarizing everything by clusterID
# OccEst2016bycluster<-OccEst2016%>%
#   group_by(clusterID)%>%
#   summarize(
#     across(
#       .cols = everything(), # Selects all columns not in group_by
#       .fns = mean,         # Applies the mean function
#       .names = "{.col}_mean" # Renames the new columns
#     )
#   )
# #get median, 2.5%, and 97.5% quantile
# OccEst2016bycluster_quantiles <- OccEst2016bycluster %>%
#   rowwise() %>%
#   mutate(
#     q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
#     q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
#     q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
#   ) %>%
#   ungroup()%>%
#   select(clusterID, q2.5,q50,q97.5)%>%
#   mutate(year=2016)


OccEst2016Xhabitat<-cbind(siteCovFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2016Xhabitat$locationC<-NULL 
OccEst2016Xhabitat$clusterID<-NULL 
OccEst2016Xhabitat$HabSuit4ClassNumeric<-NULL 
#drop before summarizing everything by habitat class
OccEst2016byhabitat<-OccEst2016Xhabitat%>%
  group_by(HabSuit4Class)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2016byhabitat_quantiles <- OccEst2016byhabitat %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(HabSuit4Class, q2.5,q50,q97.5)%>%
  mutate(year=2016)


#Probabilities of occu. at year 2 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 2, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
# #meanOccEst2017<-rowMeans(z.pred)
# OccEst2017<-cbind(obsFiltFinal, z.pred)
# #I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
# #in columns 3-18002
# OccEst2017$locationC<-NULL 
# OccEst2017$HabSuit4Class<-NULL
# #drop before summarizing everything by clusterID
# OccEst2017bycluster<-OccEst2017%>%
#   group_by(clusterID)%>%
#   summarize(
#     across(
#       .cols = everything(), # Selects all columns not in group_by
#       .fns = mean,         # Applies the mean function
#       .names = "{.col}_mean" # Renames the new columns
#     )
#   )
# #get median, 2.5%, and 97.5% quantile
# OccEst2017bycluster_quantiles <- OccEst2017bycluster %>%
#   rowwise() %>%
#   mutate(
#     q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
#     q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
#     q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
#   ) %>%
#   ungroup()%>%
#   select(clusterID, q2.5,q50,q97.5)%>%
#   mutate(year=2017)


OccEst2017Xhabitat<-cbind(siteCovFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2017Xhabitat$locationC<-NULL 
OccEst2017Xhabitat$clusterID<-NULL 
#drop before summarizing everything by habitat class
OccEst2017byhabitat<-OccEst2017Xhabitat%>%
  group_by(HabSuit4Class)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2017byhabitat_quantiles <- OccEst2017byhabitat %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(HabSuit4Class, q2.5,q50,q97.5)%>%
  mutate(year=2017)


#Probabilities of occu. at year 3 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 3, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2018<-rowMeans(z.pred)
# OccEst2018<-cbind(obsFiltFinal, z.pred)
# #I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
# #in columns 3-18002
# OccEst2018$locationC<-NULL 
# OccEst2018$HabSuit4Class<-NULL
# #drop before summarizing everything by clusterID
# OccEst2018bycluster<-OccEst2018%>%
#   group_by(clusterID)%>%
#   summarize(
#     across(
#       .cols = everything(), # Selects all columns not in group_by
#       .fns = mean,         # Applies the mean function
#       .names = "{.col}_mean" # Renames the new columns
#     )
#   )
# #get median, 2.5%, and 97.5% quantile
# OccEst2018bycluster_quantiles <- OccEst2018bycluster %>%
#   rowwise() %>%
#   mutate(
#     q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
#     q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
#     q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
#   ) %>%
#   ungroup()%>%
#   select(clusterID, q2.5,q50,q97.5)%>%
#   mutate(year=2018)

OccEst2018Xhabitat<-cbind(siteCovFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2018Xhabitat$locationC<-NULL 
OccEst2018Xhabitat$clusterID<-NULL 
#drop before summarizing everything by habitat class
OccEst2018byhabitat<-OccEst2018Xhabitat%>%
  group_by(HabSuit4Class)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2018byhabitat_quantiles <- OccEst2018byhabitat %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(HabSuit4Class, q2.5,q50,q97.5)%>%
  mutate(year=2018)


#Probabilities of occu. at year 4 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 4, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2019<-rowMeans(z.pred)
# OccEst2019<-cbind(obsFiltFinal, z.pred)
# #I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
# #in columns 3-18002
# OccEst2019$locationC<-NULL 
# OccEst2019$HabSuit4Class<-NULL
# #drop before summarizing everything by clusterID
# OccEst2019bycluster<-OccEst2019%>%
#   group_by(clusterID)%>%
#   summarize(
#     across(
#       .cols = everything(), # Selects all columns not in group_by
#       .fns = mean,         # Applies the mean function
#       .names = "{.col}_mean" # Renames the new columns
#     )
#   )
# #get median, 2.5%, and 97.5% quantile
# OccEst2019bycluster_quantiles <- OccEst2019bycluster %>%
#   rowwise() %>%
#   mutate(
#     q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
#     q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
#     q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
#   ) %>%
#   ungroup()%>%
#   select(clusterID, q2.5,q50,q97.5)%>%
#   mutate(year=2019)

OccEst2019Xhabitat<-cbind(siteCovFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2019Xhabitat$locationC<-NULL 
OccEst2019Xhabitat$clusterID<-NULL 
#drop before summarizing everything by habitat class
OccEst2019byhabitat<-OccEst2019Xhabitat%>%
  group_by(HabSuit4Class)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2019byhabitat_quantiles <- OccEst2019byhabitat %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(HabSuit4Class, q2.5,q50,q97.5)%>%
  mutate(year=2019)



#Probabilities of occu. at year 5 (second number index for stored phi and gamma values=6)
z.pred <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("z[", j, ",", 5, "]", sep=""))
  z.pred[j, ] <- ggd$value[indx]#Error in psi.pred[j, ] <- ggd$value[indx] : replacement has length zero
}  
str(z.pred)
#meanOccEst2021<-rowMeans(z.pred)
# OccEst2021<-cbind(obsFiltFinal, z.pred)
# #I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
# #in columns 3-18002
# OccEst2021$locationC<-NULL 
# OccEst2021$HabSuit4Class<-NULL
# #drop before summarizing everything by clusterID
# OccEst2021bycluster<-OccEst2021%>%
#   group_by(clusterID)%>%
#   summarize(
#     across(
#       .cols = everything(), # Selects all columns not in group_by
#       .fns = mean,         # Applies the mean function
#       .names = "{.col}_mean" # Renames the new columns
#     )
#   )
# #get median, 2.5%, and 97.5% quantile
# OccEst2021bycluster_quantiles <- OccEst2021bycluster %>%
#   rowwise() %>%
#   mutate(
#     q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
#     q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
#     q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
#   ) %>%
#   ungroup()%>%
#   select(clusterID, q2.5,q50,q97.5)%>%
#   mutate(year=2021)

OccEst2021Xhabitat<-cbind(siteCovFinal, z.pred)
#I have 18000 estimates (1's and 0's) that I want to summarize by clusterID
#in columns 3-18002
OccEst2021Xhabitat$locationC<-NULL 
OccEst2021Xhabitat$clusterID<-NULL 
#drop before summarizing everything by habitat class
OccEst2021byhabitat<-OccEst2021Xhabitat%>%
  group_by(HabSuit4Class)%>%
  summarize(
    across(
      .cols = everything(), # Selects all columns not in group_by
      .fns = mean,         # Applies the mean function
      .names = "{.col}_mean" # Renames the new columns
    )
  )
#get median, 2.5%, and 97.5% quantile
OccEst2021byhabitat_quantiles <- OccEst2021byhabitat %>%
  rowwise() %>%
  mutate(
    q2.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.025),
    q50 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.50),
    q97.5 = quantile(c_across(`1_mean`:`18000_mean`), probs = 0.975)
  ) %>%
  ungroup()%>%
  select(HabSuit4Class, q2.5,q50,q97.5)%>%
  mutate(year=2021)



finalOccDFxhabitat<-rbind(OccEst2016byhabitat_quantiles,
                          OccEst2017byhabitat_quantiles,
                          OccEst2018byhabitat_quantiles,
                          OccEst2019byhabitat_quantiles,
                          OccEst2021byhabitat_quantiles)
#Habitat A
finalOccDF_habitatA<-finalOccDFxhabitat%>%
  filter(HabSuit4Class=="A")

OccuhabclassA <- ggplot(data=finalOccDF_habitatA, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Group A: ",nrow(siteCovFinal%>%dplyr::filter(HabSuit4Class=="A"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Habitat B
finalOccDF_habitatB<-finalOccDFxhabitat%>%
  filter(HabSuit4Class=="B")

OccuhabclassB <- ggplot(data=finalOccDF_habitatB, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Group B: ",nrow(siteCovFinal%>%dplyr::filter(HabSuit4Class=="B"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Habitat C
finalOccDF_habitatC<-finalOccDFxhabitat%>%
  filter(HabSuit4Class=="C")

OccuhabclassC <- ggplot(data=finalOccDF_habitatC, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Group C: ",nrow(siteCovFinal%>%dplyr::filter(HabSuit4Class=="C"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

#Habitat D
finalOccDF_habitatD<-finalOccDFxhabitat%>%
  filter(HabSuit4Class=="D")

OccuhabclassD <- ggplot(data=finalOccDF_habitatD, aes(x=year, y=q50)) +
  geom_ribbon(aes(x=year, ymin=q2.5, ymax=q97.5), fill="blue", alpha=0.2)+
  geom_line(aes(x=year, y=q50), colour="blue")+
  my.theme+ggtitle(paste0("Group D: ",nrow(siteCovFinal%>%dplyr::filter(HabSuit4Class=="D"))," sites"))+xlab("Year")+ylab("Pred. occu.")+
  ylim(-0.5,1.25)+xlim(2016,2021)+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(2016,2021)))))

p4<-grid.arrange(OccuhabclassA,
                 OccuhabclassB,
                 OccuhabclassC,
                 OccuhabclassD,
                 nrow=2,ncol=2)
ggsave(filename="2_outputs/2_jags_outputs_feb2026/HedleyModel_2016_2021/HedleyModel_MedianPlus95BCIoccupancy_2016_2021xhabitat.png",
       plot = p4,
       width = 12,
       height = 8,
       units = c("in"),
       dpi = 300)


