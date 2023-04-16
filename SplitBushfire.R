#Ronja Seitz
#07.04.2023

###Burn Severity Assesment of a bushfire in 2017 near Split Croatia and Analysis of Regrowth the following years

library(terra)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)

setwd("/Uni/Programming/") #set working directory

##Load data
#Load Pre Fire Image
preFire <- "Material/S2A_MSIL2A_20170707T095031_N0205_R079_T33TXJ_20170707T095257.SAFE" #path to one sentinel2 image
l <- list.files(path = preFire, recursive = TRUE, full.names=TRUE, pattern ="20m.jp2$") #search for all the bands
length(l) #how many files do we have

preFire <- rast(l[2:10]) #rasterinzing the multisepctral bands
names(preFire) <- c("B02", "B03","B04","B05", "B06", "B07", "B11", "B12","B08a") #renaming

#Load Post Fire Image
postFire <- "Material/S2A_MSIL2A_20170806T095031_N0205_R079_T33TXJ_20170806T095744.SAFE" #path to one sentinel2 image
l <- list.files(path = postFire, recursive = TRUE, full.names=TRUE, pattern ="20m.jp2$") #search for all the bands
length(l) #how many files do we have

postFire <- rast(l[2:10]) #rasterizing the multisepctral bands
names(postFire) <- c("B02", "B03","B04","B05", "B06", "B07", "B11", "B12","B08a") #renaming

##load aoi and process it 
aoi <- st_read("BushfireSplit_AOI.gpkg")# Create a geopackage file
aoi <- vect(aoi)
aoi <- terra::project(aoi, "EPSG:32633") #project aoi to crs of sentinel2 scene

#run second time if it did not work the first time
preFire <- crop(preFire,aoi) #crop sentinel pre fire scene to aoi
plot(preFire$B02, main = "Pre Fire Band 2") #plot one band

postFire <- crop(postFire,aoi) #crop sentinel post fire scene to aoi
plot(postFire$B02, main = "Post Fire Band 2") #plot one band

pal <- colorRampPalette(rev(brewer.pal(11, 'PiYG')))(200) #define color palelette

plotRGB(preFire, r="B04", g="B03", b="B02", stretch = "lin") #plot true color image pre fire
plotRGB(preFire, r="B12", g="B08a", b="B04", stretch = "lin") #plot false color image of pre fire

plotRGB(postFire, r="B04", g="B03", b="B02", stretch = "lin") #plot true color image post fire
plotRGB(postFire, r="B12", g="B08a", b="B04", stretch = "lin") #plot false color image post fire

##create water mask, because indices are sensitive to water
water <- rast("ASTWBDV001_N43E016/ASTWBDV001_N43E016_att.tif") #load water dataset
unique(water) #check values of water dataset

water <- terra::project(water, "EPSG:32633") #project water to crs of sentinel2 scene
#water <- classify(x = water, cbind(1, NA))
water <- as.polygons(water) #change from raster to vector
water <- water[1] #index raster to get everything except ocean
water <- crop(water,aoi) #crop water scene to aoi

par(mfrow=c(2,1)) #two images in one column
preFire <- mask(x = preFire, mask = water, inverse = F) #mask pre fire scene with water mask
plot(preFire$B02, main = "Pre Fire Band 2 With Extent of Water Mask") #check if masking worked with one band
plot(water, add = T) #add water mask to plot

postFire <- mask(x = postFire, mask = water)
plot(postFire$B02, main = "Post Fire Band 2 With Extent of Water Mask")
plot(water, add = T)

##calculate NBR+, an advanced burned area index
NBRPlus_fun <- function(img,B12,B08a,B03,B02){
  NBR_Plus <- (B12-B08a-B03-B02)/(B12+B08a+B03+B02)
}

postF_NBRPlus <- NBRPlus_fun(img = postFire,B12 = postFire$B12, B08a = postFire$B08a,
                       B03 = postFire$B03, B02 = postFire$B02) 


##burn severity nbr+
preF_NBRPlus <- NBRPlus_fun(img = preFire, B12 = preFire$B12, B08a = preFire$B08a, 
                               B03 = preFire$B03, B02 = preFire$B02)
DNBR_Plus <- preF_NBRPlus - postF_NBRPlus
plot(DNBR_Plus, main = "Burn Severity NBR+", col = pal)

##calculate NBR, an established and often used fire index
NBR_fun <- function(img,B08a,B12){
  NBR <- (B08a-B12)/(B08a+B12)
}

preFireNBR <- NBR_fun(img = preFire, B08a = preFire$B08a, B12 = preFire$B12)
postFireNBR <- NBR_fun(img = postFire, B08a = postFire$B08a, B12 = postFire$B12)

#Burn severity
DNBR <- preFireNBR - postFireNBR 
DNBR_breaks <- c(-0.500,-0.251,-0.101,0.99,0.269,0.439,0.659,1.3)

plot(DNBR, breaks = DNBR_breaks, main = "Burn Severity NBR", col = pal)
DNBR_df <- as.data.frame(x = DNBR, xy = T) #as dataframe
#plot with ggplot
ggplot(DNBR_df)+geom_raster( aes(x = x, y = y,fill = B08a ))+
  scale_fill_distiller(palette = "Spectral")

#get only burned area
dnbr_mask <-  DNBR
dnbr_mask[dnbr_mask$B08a<0.100] <- -0.5
plot(dnbr_mask)

##calculate NDVI (NIR & RED)
NDVI_fun <- function(img, B08a, B4){
  NDVI <- (B08a-B4)/(B08a+B4)
}

par(mfrow=c(1,1)) #define how many figures are shown on one page

NDVI_2017pre <- NDVI_fun(img = preFire, B08a = preFire$B08a, B4 = preFire$B04) #calculate NDVI before the bushfire
names(NDVI_2017pre) <- c("NDVI")
plot(NDVI_2017pre, main = "Pre Fire NDVI") #plot
NDVI_2017post <- NDVI_fun(img = postFire, B08a = postFire$B08a, B4 = postFire$B04) #calcilate NDVI after the bushfire
names(NDVI_2017post) <- c("NDVI")
plot(NDVI_2017post, main = "Post Fire NDVI") #plot

DNDVI_2017 <- NDVI_2017post - NDVI_2017pre #calculate change
plot(DNDVI_2017, main = "Difference Post & Pre Fire NDVI") #plot

##get burn area for later comparison
#reclassify burn severity of NBR to get burned and non-burned area
m <- c(-1,0.2,NA, #determine burned and non-burned areas
       0.2,1.4,1) #burned areas have value 1 and non-burned areas NA
rcmat <- matrix(m, ncol = 3, byrow = T)
DNBR_rc <- classify(DNBR, rcmat)  #reclassify DNBR

plotRGB(postFire, r="B12", g="B08a", b="B04", stretch = "lin") #plot post fire false color image
plot(DNBR_rc,col = "blue", add = T) #plot reclassified DNBR as overlay

#get non-burned area
m_non <- c(-1,0.2,1, #determine burned and non-burned areas
       0.2,1.4,NA) #burned areas have value NA and non-burned areas 1
rcmat_non <- matrix(m_non, ncol = 3, byrow = T)
DNBR_rc_non <- classify(DNBR,rcmat_non)  #reclassify DNBR

plotRGB(postFire, r="B12", g="B08a", b="B04", stretch = "lin") #plot post fire false color image
plot(DNBR_rc_non,col = "blue", add = T) #plot reclassified DNBR as overlay

##Get areas in km2
#burned
burned_area <- as.data.frame(DNBR_rc, xy = T) %>%
  summarise(pixels = n()) %>%
  mutate(`area m^2` = pixels * 20 * 20, #multiply by pixel area
                  `area km^2` = round(`area m^2`/1e6,2))
#non-burned
nonburned_area <- as.data.frame(DNBR_rc_non, xy = T) %>%
  summarise(pixels = n()) %>%
  mutate(`area m^2` = pixels * 20 * 20, #multiply by pixel area
         `area km^2` = round(`area m^2`/1e6,2))

print(burned_area)
print(nonburned_area)

##get NDVI values within burned area and outside
#within
NDVI_2017a <- data.frame()
NDVI_2017a <- rbind(terra::zonal(NDVI_2017pre,DNBR_rc,fun = mean))
NDVI_2017in_sum <- terra::zonal(NDVI_2017post,DNBR_rc,fun = mean)
NDVI_2017a <- rbind(NDVI_2017a, NDVI_2017in_sum )
NDVI_2017a <- NDVI_2017a[2]

#outside
NDVI_2017b <- data.frame()
NDVI_2017b <- rbind(terra::zonal(NDVI_2017pre,DNBR_rc_non,fun = mean))
NDVI_2017out_sum <- terra::zonal(NDVI_2017post,DNBR_rc_non,fun = mean)
NDVI_2017b <- rbind(NDVI_2017b, NDVI_2017out_sum)
NDVI_2017b <- NDVI_2017b[2]

#merge both df
NDVI_2017 <- data.frame()
NDVI_2017 <- cbind(NDVI_2017a,NDVI_2017b)
names(NDVI_2017) <- c("IN","OUT")
row.names(NDVI_2017) <- c("preFire","postFire")
NDVI_2017


##calculate ratio of regrowth over the following years using the NDVI
scenes <- list.files(pattern = ".*SAFE")
scenes <- scenes[4:length(scenes)]
rasters <- list()

for (i in 1:length(scenes)) {
  files <- list.files(scenes[i], recursive = TRUE, full.names=TRUE, pattern ="B8A_20m.jp2$|4_20m.jp2$")%>%
    rast()
  
  rasters[[i]] <- files
}

rst <- rast(rasters)
rst <- crop(rst,aoi)    
rst <- mask(x = rst, mask = water)

#calculate NDVI of multiple yearsn
NDVI_2018 <- NDVI_fun(img = rst, B08a = rst[[2]], B4 = rst[[1]])
NDVI_2019 <- NDVI_fun(img = rst, B08a = rst[[4]], B4 = rst[[3]])
NDVI_2020 <- NDVI_fun(img = rst, B08a = rst[[6]], B4 = rst[[5]])
NDVI_2021 <- NDVI_fun(img = rst, B08a = rst[[8]], B4 = rst[[7]])
NDVI_2022 <- NDVI_fun(img = rst, B08a = rst[[10]], B4 = rst[[9]])
NDVI_2023 <- NDVI_fun(img = rst, B08a = rst[[12]], B4 = rst[[11]])

#ndvi stats inside and outside of burned area
#2018
#within
NDVI_2018in_df <- data.frame()
NDVI_2018in_df <- rbind(terra::zonal(NDVI_2018,DNBR_rc,fun = mean))
NDVI_2018in_df <- NDVI_2018in_df[2]
names(NDVI_2018in_df) <- c("NDVI_IN")

#outside
NDVI_2018out_df <- data.frame()
NDVI_2018out_df <- rbind(terra::zonal(NDVI_2018,DNBR_rc_non,fun = mean))
NDVI_2018out_df <- NDVI_2018out_df[2]
names(NDVI_2018out_df) <- c("NDVI_OUT")

#2019
#within
NDVI_2019in_df <- data.frame()
NDVI_2019in_df <- rbind(terra::zonal(NDVI_2019,DNBR_rc,fun = mean))
NDVI_2019in_df <- NDVI_2019in_df[2]
names(NDVI_2019in_df) <- c("NDVI_IN")

#outside
NDVI_2019out_df <- data.frame()
NDVI_2019out_df <- rbind(terra::zonal(NDVI_2019,DNBR_rc_non,fun = mean))
NDVI_2019out_df <- NDVI_2019out_df[2]
names(NDVI_2019out_df) <- c("NDVI_OUT")

#2020
#within
NDVI_2020in_df <- data.frame()
NDVI_2020in_df <- rbind(terra::zonal(NDVI_2020,DNBR_rc,fun = mean))
NDVI_2020in_df <- NDVI_2020in_df[2]
names(NDVI_2020in_df) <- c("NDVI_IN")

#outside
NDVI_2020out_df <- data.frame()
NDVI_2020out_df <- rbind(terra::zonal(NDVI_2020,DNBR_rc_non,fun = mean))
NDVI_2020out_df <- NDVI_2020out_df[2]
names(NDVI_2020out_df) <- c("NDVI_OUT")

#2021
#within
NDVI_2021in_df <- data.frame()
NDVI_2021in_df <- rbind(terra::zonal(NDVI_2021,DNBR_rc,fun = mean))
NDVI_2021in_df <- NDVI_2021in_df[2]
names(NDVI_2021in_df) <- c("NDVI_IN")

#outside
NDVI_2021out_df <- data.frame()
NDVI_2021out_df <- rbind(terra::zonal(NDVI_2021,DNBR_rc_non,fun = mean))
NDVI_2021out_df <- NDVI_2021out_df[2]
names(NDVI_2021out_df) <- c("NDVI_OUT")

#2022
#within
NDVI_2022in_df <- data.frame()
NDVI_2022in_df <- rbind(terra::zonal(NDVI_2022,DNBR_rc,fun = mean))
NDVI_2022in_df <- NDVI_2022in_df[2]
names(NDVI_2022in_df) <- c("NDVI_IN")

#outside
NDVI_2022out_df <- data.frame()
NDVI_2022out_df <- rbind(terra::zonal(NDVI_2022,DNBR_rc_non,fun = mean))
NDVI_2022out_df <- NDVI_2022out_df[2]
names(NDVI_2022out_df) <- c("NDVI_OUT")

#2023
#within
NDVI_2023in_df <- data.frame()
NDVI_2023in_df <- rbind(terra::zonal(NDVI_2023,DNBR_rc,fun = mean))
NDVI_2023in_df <- NDVI_2023in_df[2]
names(NDVI_2023in_df) <- c("NDVI_IN")

#outside
NDVI_2023out_df <- data.frame()
NDVI_2023out_df <- rbind(terra::zonal(NDVI_2023,DNBR_rc_non,fun = mean))
NDVI_2023out_df <- NDVI_2023out_df[2]
names(NDVI_2023out_df) <- c("NDVI_OUT")

#merge df from 2018 to 2023
NDVI_multiY <- data.frame()
NDVI_multiY <- rbind(NDVI_2018out_df,NDVI_2019out_df,NDVI_2020out_df,NDVI_2021out_df,NDVI_2022out_df,NDVI_2023out_df)
NDVI_multiY$NDVI_IN <- rbind(NDVI_2018in_df,NDVI_2019in_df,NDVI_2020in_df,NDVI_2021in_df,NDVI_2022in_df,NDVI_2023in_df)
NDVI_multiY$DIF <- NDVI_multiY$NDVI_OUT - NDVI_multiY$NDVI_IN
print(NDVI_multiY)

#get rid of weird naming of columns and plotting
NDVI_multiY_mat <- as.matrix(NDVI_multiY)
NDVI_multiY_mat <- matrix(NDVI_multiY_mat, ncol = ncol(NDVI_multiY), dimnames = NULL)

par(mar = c(5.25, 4.25, 4.25, 4.25))
matplot(NDVI_multiY_mat,
  type = "l",
  lty = 1,
  main = "Change of NDVI over the years 2018 to 2023",
  ylab = "NDVI")
legend("bottomleft", legend = c("Out", "IN","DIF"), lty = 1, col = 1:2)
