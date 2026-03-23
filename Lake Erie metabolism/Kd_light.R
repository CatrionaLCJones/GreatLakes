########################## Kd (light attenuation) ##############################

############################################################################
##################### MODEL COMPONDENT Kd light attenuation #######################
############################################################################


# Kd (light attenuation coefficient)
# Kd = 1/z*ln (E 0 /E z ) where E 0 is the surface incident irradiance at 
# 0 m and E z is the irradiance at z m depth.

# Components
# Surface PAR (E0) - shallowest par measurement on ctd profile
# Sub-surface PAR (Ez) - deepest par measurement on ctd profile
# Depth (z) - difference between shallowest and deepest measurements

setwd("C:/Users/Katie/Documents/Research/Purdue_post_doc/Lake Erie Metabolism model/Datasets/csv files")
ctd <- read.csv('NOAA_ctd_201718.csv')

################################################################################
################################################################################

#This dataframe contains 2 years, deal with them separately

ctd_2018 <- ctd[str_detect(ctd$datetime, "2018"), ] 
ctd_2017 <- ctd[str_detect(ctd$datetime, "2017"), ] 

################################################################################
############################### 2018 ########################################
################################################################################

#Assign loop list
site_list <- c("CHRP1","CHRP2","CHRP3","CHRP4","CHRP5","CHRP7","CHRP8","CHRP9")

for (g in site_list) {
  
  subset <- ctd_2018 %>% filter(site==g)
  
  site <- subset[1,1]
  datetime <- subset[1,4]
  lat <- subset[1,2]
  long <- subset[1,3]
  
  # need to filter on loop list value
  max_depth <- max(subset$depth)
  min_depth <- min(subset$depth)
  k0 <- subset %>% slice_head(n=1) %>% pull(par)
  kz <-subset %>% slice_tail(n=1) %>% pull(par)
  z = max_depth - min_depth
  
  Kd = 1/z*log(k0/kz)
  
  Light_atten <- data.frame(site, lat, long, datetime, Kd)
  
  assign(paste0("Light_atten_", g), Light_atten)    
  
}

Light_atten_2018 <- rbind(Light_atten_CHRP1, Light_atten_CHRP2, Light_atten_CHRP3,
                          Light_atten_CHRP4, Light_atten_CHRP5,  Light_atten_CHRP7,
                          Light_atten_CHRP8, Light_atten_CHRP9)

write.csv(Light_atten_2018, file="C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/Kd_2018.csv")

################################################################################
############################### 2017 ########################################
################################################################################

#Assign loop list
site_list <- c("CHRP1","CHRP3","CHRP4","CHRP5","CHRP7","CHRP8","CHRP9")

for (g in site_list) {
  
  subset <- ctd_2017 %>% filter(site==g)
  
  site <- subset[1,1]
  datetime <- subset[1,4]
  lat <- subset[1,2]
  long <- subset[1,3]
  
  # need to filter on loop list value
  max_depth <- max(subset$depth)
  min_depth <- min(subset$depth)
  k0 <- subset %>% slice_head(n=1) %>% pull(par)
  kz <-subset %>% slice_tail(n=1) %>% pull(par)
  z = max_depth - min_depth
  
  Kd = 1/z*log(k0/kz)
  
  Light_atten <- data.frame(site, lat, long, datetime, Kd)
  
  assign(paste0("Light_atten_", g), Light_atten)    
  
}

Light_atten_2017 <- rbind(Light_atten_CHRP1, Light_atten_CHRP3,
                          Light_atten_CHRP4, Light_atten_CHRP5,  Light_atten_CHRP7,
                          Light_atten_CHRP8, Light_atten_CHRP9)

write.csv(Light_atten_2017, file="C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/Kd_2017.csv")
