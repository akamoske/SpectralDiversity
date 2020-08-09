# load the needed library
library(dplyr)
library(vegan)
library(raster)
library(rgdal)
library(rgeos)
library(geometry)
library(matrixStats) 

# read in reflectance data
all.ref <- read.csv(paste0("C:/Users/Aaron Kamoske/Dropbox/Publications/Kamoske_Dissertation_Ch3/",
                           "data/reflectance/plot_reflectance_20200412.csv"))

#---------------------------------------------------------------------------------------#
# PCA based metrics
#---------------------------------------------------------------------------------------#

##### RUN A PCA #####

# set the data up
dim(all.ref)
all.ref <- na.omit(all.ref)
dim(all.ref)
pca.data <- all.ref[,-1]

# run a pca
site.pca <- prcomp(pca.data,
                   center = TRUE,
                   scale = FALSE)

# see how well this did
summary(site.pca)

# the first two components account for 97.51% of the variation in the dataset so we will use those
# so we need to save the loading values for the first 2 components
pc1.loadings <- site.pca$rotation[,1]
pc2.loadings <- site.pca$rotation[,2]

# apply the loadings to each of the datasets
pca.data2 <- scale(pca.data, center = TRUE, scale = FALSE)
pc.df <- all.ref

pc1 <- pca.data2 * pc1.loadings
pc.df$pc1 <- rowSums(pc1)

pc2 <- pca.data2 * pc2.loadings
pc.df$pc2 <- rowSums(pc2)

pc.df <- pc.df[,c(1,320,321)]
colnames(pc.df) <- c("plotID", "PC1", "PC2")

##### CONVEX HULL VOLUME #####

# calculate species indices for each plot
ch <- as.data.frame(pc.df %>%
                      group_by(plotID) %>%
                      summarise(convexHull = convhulln(na.omit(data.frame(PC1, PC2)), options = "FA")$vol))

# make a final dataset to save all subsequent variables too
final.data <- ch

##### SUM OF SQUARES #####

# rename the dataframe
ss.df <- pc.df

# create new columns with the centered data
ss.df$pc1sc <- raster::scale(ss.df$PC1, center = TRUE, scale = FALSE)^2
ss.df$pc2sc <- raster::scale(ss.df$PC2, center = TRUE, scale = FALSE)^2

# sum the rows across these two columns
ss.df$ss <- ss.df$pc1sc + ss.df$pc2sc

# set 0 to NA
ss.df$ss[ss.df$ss <= 0] <- NA

# write spectral diversity function
spec.div.sum <- function(x) {
  if (all(is.na(x)) == TRUE) {
    return(NA)
  } else {
    sp.di <- sum(na.omit(x)) / (length(na.omit(x)) - 1)
    return(sp.di)
  }
}

# calculate species indices for each plot
sos <- as.data.frame(ss.df %>%
                       group_by(plotID) %>%
                       summarise(sumSquares_sum = spec.div.sum(ss))) 

# add the data to the final dataset
final.data <- merge(final.data, sos, by = "plotID")

##### SUM OF VARIANCE #####

# calculate species indices for each plot
sv <- as.data.frame(pc.df %>%
                      group_by(plotID) %>%
                      summarise(variance_sum = sum(var(PC1) + var(PC2))))

# add the data to the final dataset
final.data <- merge(final.data, sv, by = "plotID")

##### Principal Components #####

# calculate species indices for each plot
prin.c1 <- as.data.frame(pc.df %>%
                           group_by(plotID) %>%
                           summarise(pc1_mean = mean(PC1, na.rm = TRUE),
                                     pc1_sd = sd(PC1, na.rm = TRUE),
                                     pc1_max = max(PC1, na.rm = TRUE),
                                     pc1_min = min(PC1, na.rm = TRUE),
                                     pc1_range = pc1_max - pc1_min))

# add the data to the final dataset
final.data <- cbind(final.data, prin.c1[,2:ncol(prin.c1)])

# calculate species indices for each plot
prin.c2 <- as.data.frame(pc.df %>%
                           group_by(plotID) %>%
                           summarise(pc2_mean = mean(PC2, na.rm = TRUE),
                                     pc2_sd = sd(PC2, na.rm = TRUE),
                                     pc2_max = max(PC2, na.rm = TRUE),
                                     pc2_min = min(PC2, na.rm = TRUE),
                                     pc2_range = pc2_max - pc2_min))

# add the data to the final dataset
final.data <- cbind(final.data, prin.c2[,2:ncol(prin.c2)])
                                  
#---------------------------------------------------------------------------------------#
# reflectance based metrics
#---------------------------------------------------------------------------------------#

##### COEFFICIENT OF VARIATION #####

# resave the reflectance data
co.va <- all.ref
colnames(co.va)[1] <- "plotID"

# set up data for saving
final.data$coefVariantion <- 0

# calcalute coefficient of variation for each 
for (i in 1:length(unique(co.va$plotID))) {
  
  # pull out the plot id
  plot.id <- unique(co.va$plotID)[i]
  
  # subset the data
  df.sub <- co.va[co.va$plotID == plot.id,]
  
  # find column means
  mean.col <- colMeans(df.sub[,-1])
  
  # find column sd
  sd.col <- co.va[,-1] %>%
    summarise_each(funs(sd(., na.rm=TRUE)))
  
  # calculate coefficient of variation for each wavelength
  coef.var <- rowMeans(sd.col / mean.col)
  
  # save data to final dataset
  final.data[i,ncol(final.data)] <- coef.var
}

##### NDVI #####

# calculate a narrow band NDVI for each pixel
co.va$NDVI <- (co.va$nm832.300781 - co.va$nm676.982178) / (co.va$nm832.300781 + co.va$nm676.982178)

# only take the columns we need
ndvi <- co.va[,c(1,321)]

# calculate species indices for each plot
ndvi.sd <- as.data.frame(ndvi %>%
                           group_by(plotID) %>%
                           summarise(NDVI_mean = mean(NDVI, na.rm = TRUE),
                                     NDVI_sd = sd(NDVI, na.rm = TRUE),
                                     NDVI_max = max(NDVI, na.rm = TRUE),
                                     NDVI_min = min(NDVI, na.rm = TRUE),
                                     NDVI_range = NDVI_max - NDVI_min))
# add the data to the final dataset
final.data <- cbind(final.data, ndvi.sd[,2:ncol(ndvi.sd)])

##### PRI #####

# calculate a narrow band PRI for each pixel
co.va$PRI <- (co.va$nm531.684082 - co.va$nm571.766296) / (co.va$nm531.684082 + co.va$nm571.766296)

# only take the columns we need
pri <- co.va[,c(1,322)]

# calculate species indices for each plot
pri.sd <- as.data.frame(pri %>%
                          group_by(plotID) %>%
                          summarise(PRI_mean = mean(PRI, na.rm = TRUE),
                                    PRI_sd = sd(PRI, na.rm = TRUE),
                                    PRI_max = max(PRI, na.rm = TRUE),
                                    PRI_min = min(PRI, na.rm = TRUE),
                                    PRI_range = PRI_max - PRI_min))

# add the data to the final dataset
final.data <- cbind(final.data, pri.sd[,2:ncol(pri.sd)])

##### RVSI #####

# calculate a narrow band RVSI for each pixel
# ((714+752) / 2) -733
co.va$RVSI <- ((co.va$nm712.054199 + co.va$nm752.136414) / 2) - co.va$nm732.095276

# only take the columns we need
rvsi <- co.va[,c(1,323)]

# calculate species indices for each plot
rvsi.sd <- as.data.frame(rvsi %>%
                           group_by(plotID) %>%
                           summarise(RVSI_mean = mean(RVSI, na.rm = TRUE),
                                     RVSI_sd = sd(RVSI, na.rm = TRUE),
                                     RVSI_max = max(RVSI, na.rm = TRUE),
                                     RVSI_min = min(RVSI, na.rm = TRUE),
                                     RVSI_range = RVSI_max - RVSI_min))

# add the data to the final dataset
final.data <- cbind(final.data, rvsi.sd[,2:ncol(rvsi.sd)])

##### RED EDGE NDVI #####
# (750-705) / (750+705) - Ollinger 2011

co.va$reNDVI <- (co.va$nm752.136414 - co.va$nm707.043884) / (co.va$nm752.136414 + co.va$nm707.043884)

# only take the columns we need
rendvi <- co.va[,c(1,324)]

# calculate species indices for each plot
rendvi.sd <- as.data.frame(rendvi %>%
                           group_by(plotID) %>%
                           summarise(reNDVI_mean = mean(reNDVI, na.rm = TRUE),
                                     reNDVI_sd = sd(reNDVI, na.rm = TRUE),
                                     reNDVI_max = max(reNDVI, na.rm = TRUE),
                                     reNDVI_min = min(reNDVI, na.rm = TRUE),
                                     reNDVI_range = reNDVI_max - reNDVI_min))

# add the data to the final dataset
final.data <- cbind(final.data, rendvi.sd[,2:ncol(rendvi.sd)])

##### SWIR 1 - 1500nm-1800nm#####

co.va$SWIR1 <- rowMeans(co.va[,c(181:240)])

# only take the columns we need
swir1 <- co.va[,c(1,325)]

# calculate species indices for each plot
swir1.sd <- as.data.frame(swir1 %>%
                             group_by(plotID) %>%
                             summarise(SWIR1_mean = mean(SWIR1, na.rm = TRUE),
                                       SWIR1_sd = sd(SWIR1, na.rm = TRUE),
                                       SWIR1_max = max(SWIR1, na.rm = TRUE),
                                       SWIR1_min = min(SWIR1, na.rm = TRUE),
                                       SWIR1_range = SWIR1_max - SWIR1_min))

# add the data to the final dataset
final.data <- cbind(final.data, swir1.sd[,2:ncol(swir1.sd)])

##### SWIR 2 - 2000nm-2400nm#####

co.va$SWIR2 <- rowMeans(co.va[,c(241:319)])

# only take the columns we need
swir2 <- co.va[,c(1,326)]

# calculate species indices for each plot
swir2.sd <- as.data.frame(swir2 %>%
                            group_by(plotID) %>%
                            summarise(SWIR2_mean = mean(SWIR2, na.rm = TRUE),
                                      SWIR2_sd = sd(SWIR2, na.rm = TRUE),
                                      SWIR2_max = max(SWIR2, na.rm = TRUE),
                                      SWIR2_min = min(SWIR2, na.rm = TRUE),
                                      SWIR2_range = SWIR2_max - SWIR2_min))

# add the data to the final dataset
final.data <- cbind(final.data, swir2.sd[,2:ncol(swir2.sd)])

##### NIR - 800nm-1350nm #####

co.va$NIR <- rowMeans(co.va[,c(62:171)])

# only take the columns we need
nir <- co.va[,c(1,327)]

# calculate species indices for each plot
nir.sd <- as.data.frame(nir %>%
                            group_by(plotID) %>%
                            summarise(NIR_mean = mean(NIR, na.rm = TRUE),
                                      NIR_sd = sd(NIR, na.rm = TRUE),
                                      NIR_max = max(NIR, na.rm = TRUE),
                                      NIR_min = min(NIR, na.rm = TRUE),
                                      NIR_range = NIR_max - NIR_min))

# add the data to the final dataset
final.data <- cbind(final.data, nir.sd[,2:ncol(nir.sd)])

#---------------------------------------------------------------------------------------#
# export final data
#---------------------------------------------------------------------------------------#

# write this data to disc
write.csv(final.data,
          paste0("C:/Users/Aaron Kamoske/Dropbox/Publications/Kamoske_Dissertation_Ch3/",
                 "/data/spectral_diversity/spectral_diversity_20200423.csv"),
          row.names = FALSE)




