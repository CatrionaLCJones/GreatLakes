################# MIXING DEPTH ESTIMATION AND GAP FILLING ######################

############################################################################
##################### MODEL COMPONDENT z-mix ##############################
############################################################################
#First step is to determine the mixing depth. Use the temperature only data for this as it
#is collected at a wider range of depths than the DO (for the most part) so this gives more
#data for the function to work with. This loops for all sites/years
###############################################################################

library(LakeMetabolizer)
library(rLakeAnalyzer)
library(lubridate)
library(tidyverse)
library(data.table)
library(stringr)
library(anytime)
library(zoo)

Temperature <- list.files(pattern = '*temp.txt')
Loop_list_z <- substr(Temperature, 1, 10)

for (i in Loop_list_z) { 
  
  #Set the file path for the text file. 
  #file path for main data = "C:/Users/Katie/Documents/Research/Purdue_post_doc/Lake Erie Metabolism model/Datasets/"
  FilePath <- (paste0("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/", i, "_temp.txt"))
  
  #Load the text file in, specifying timezone as EST - the default is GMT
  Temp <- load.ts(FilePath, tz = "EST")
  
  #Remove all columns except the datetime and each temperature column
  Temp <- select(Temp, -c(latitude, longitude))
  
  #Use temperature data to get mixing depths
  Mixing_depth <- ts.thermo.depth(Temp,  na.rm=TRUE)[,2]
  #Pull the datatime out of the Temp dataframe and combine with mixing depth to get a
  #mixing depth for each datetime
  datetime <- Temp$datetime
  Mixing_depth <- data.frame(datetime, Mixing_depth)
  
  #Capture the loop output
  assign(paste0("Mixing_depth_", i), Mixing_depth)
  
}
  
  sites <- c("CHRP1",
            "CHRP2",
            "CHRP3",
            "CHRP4",
         "CHRP5", #exclude from 2018
            "CHRP7", 
            "CHRP8",
            "CHRP9") 
  
  years <- c("2017")
  
  for (y in years) {
  
  for (s in sites) {
    
    Mixing <- get(paste0("Mixing_depth_",s,"_",y))

    #interpolate y value based on x list of empty datetimes
    new_y = approx(Mixing$datetime, Mixing$Mixing_depth , xout=Mixing$datetime)
    
    new_y <- new_y %>% as.data.frame(col.names=c("datetime", "Mixing_depth"))
    
    assign((paste0("Mixing_depth_",s,"_",y)), new_y) 
  }

}
