library(LakeMetabolizer)
library(rLakeAnalyzer)
library(lubridate)
library(tidyverse)
library(data.table)
library(stringr)
library(anytime)
library(zoo)

#Read in new bubble coefficient function
source("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Scripts/Wave height/k.dieke.mel.R")

############################################################################
##################### MODEL COMPONDENT k gas flux coefficient##################
############################################################################

#There is no data for CHRP5 in 2018. Run script once with CHRP5 included in sites
#and 2018 excluded from year, then run again with CHRP5 excluded and 2018 included

years <- c("2017","2018","2019")


#For CHRP2 in 2017, no data but seems to be closest to CHRP8 in 2018/19 so use that

Light_atten_2017 <- read.csv("Kd_2017.csv")
Light_atten_2018 <- read.csv("Kd_2018.csv")
Light_atten_2019 <- read.csv("Kd_2018.csv")

for (y in years) {
  
  FilePath <- (paste("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/Weather_", y,".txt", sep=""))
  #Load the text file in, specifying timezone as EST - the default is GMT
  Climate<- load.ts(FilePath, tz = "EST")
  
  sites <- c("CHRP1","CHRP2","CHRP3","CHRP4","CHRP5","CHRP7","CHRP8","CHRP9")
  
  for (s in sites) {
    
    try({
      
      FilePath <- (paste("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/",s,"_",y,"_temp.txt", sep=""))
      #Load the text file in, specifying timezone as EST - the default is GMT
      Temp <- load.ts(FilePath, tz = "EST")
      
      #Create a single dataframe with all relevant parameters to ensure matching number of rows
      #for model
      
      #These are all single values which can be applied across all rows of the dataframe
      Kd <- as.numeric(get(paste0("Light_atten_",y)) %>% filter(site==s) %>% select(Kd))
      Weather <- Climate %>% filter(site==s)
      Weather <- Weather %>% mutate(Kd = Kd)
      Weather <- Weather %>% mutate(lake.area=25740)
      Weather <- Weather %>% mutate(wnd.z = 10)
      datetime <- get(paste0("Light_atten_",y)) %>% filter(site==s) %>% select(datetime)
      
      #For mixing depth, there are multiple values based on datetime and more rows in this dataframe
      #than in the 'weather' dataframe (3 hourly measurements vs 10 minute measurements).
      #So we need to join the following mixing depth columns: datetime, mixing depth and surface temp 
      #with all columns from the Weather dataframe. A left join will keep all the rows for Weather and only
      #bring in mixing depth columns where there's a match on our joining colum (datetime)
      #First subset the columns from mixing depth that we want to keep
      
      #Get the surfce water temperature from the Temp profile
      wtr <- Temp[,c(1,4)] #This relies on their being a 0-1m reading
      
      #Then join the dateframes by datetime
      kgas <- left_join(Weather, get(paste0("Mixing_depth_",s,"_",y)), by = c("datetime" = "datetime"))
      kgas <- left_join(kgas, wtr, by = c("datetime" = "datetime"))
      
      #Remove rows where either surface temp or mixing depth are NA
      kgas_NA <- kgas[complete.cases(kgas[ , 3]),]
      kgas_NA <- kgas_NA[complete.cases(kgas_NA[ , 17]),]
      
      #Move datetime so that it's the first column
      kgas_NA <- kgas_NA %>% relocate(datetime)
      kgas_NA <- kgas_NA %>% dplyr::rename(Ts = 18)
      kgas_NA <- kgas_NA %>% dplyr::rename(atm.press = 7)
      
      kgas_NA$datetime <- format(as.POSIXct(kgas_NA$datetime,
                                            format = "%Y-%m-%d %H:%M:%S"),
                                 format = "%Y-%m-%d %H:%M")
      
      as.POSIXct(kgas_NA$datetime, format = "%Y-%m-%d %H:%M")
    
      #Produces a timeseries vector of gas exchange coefficient
      k600_dm = k.dieke.mel(kgas_NA, wnd.z=10, 
                            Kd=kgas_NA$Kd, 
                            atm.press=kgas_NA$atm.press, 
                            lat=kgas_NA$latitude,
                            datetime=kgas_NA$datetime,
                            z.aml=kgas_NA$mixing_depth, 
                            airt=kgas_NA$airt,
                            wnd=kgas_NA$wnd, 
                            rh=kgas_NA$rh, 
                            sw=kgas_NA$sw, 
                            lwnet=kgas_NA$lwnet,
                            wvht=kgas_NA$wvht,
                            fetch=kgas_NA$fetch,
                            Ts=kgas_NA$Ts)
      
      #Then convert k600 values to gas flux
      kgas.dm.wvht = k600.2.kGAS.base(k600_dm, kgas_NA$Ts, 'O2')
      
      kgas_dm_wvht <- kgas_NA %>% mutate(kgas.dm = kgas.dm.wvht)
      
      kgas_dm_wvht[sapply(kgas_dm_wvht, is.infinite)] <- NA
      
      assign(paste0("kgas_dm_wvht_",s,"_",y), kgas_dm_wvht)    
      
      silent=TRUE})
  }
} 
