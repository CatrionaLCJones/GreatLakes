library(LakeMetabolizer)
library(rLakeAnalyzer)
library(lubridate)
library(tidyverse)
library(plyr)
library(zoo)

##################### DOES ONE YEAR AT A TIME - CHANGE MANUALLY #############################

FilePath <- ("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/MeteorologicalData_2021.txt")
#Load the text file in, specifying timezone as EST - the default is GMT
Weather <- load.ts(FilePath, tz = "EST")

FilePath <- ("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/Wind_waves_2021.txt")
#Load the text file in, specifying timezone as EST - the default is GMT
Waves <- load.ts(FilePath, tz = "EST")

sites<- c("CHRP1",
           "CHRP2",
          "CHRP3",
          "CHRP4",
        "CHRP5", #exclude from 2018
    #   "CHRP7", #exclude from 2021
          "CHRP8",
          "CHRP9") 

for (s in sites) {

FilePath <- paste0("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/",s,"_2021_temp.txt")
Temp <- load.ts(FilePath, tz = "EST")

#sub_temp <- Temp %>% filter(site==s)
sub_weather <- Weather %>% filter(site==s)
sub_wave <- Waves %>% filter(site==s)


datetime <- format(seq(as.POSIXct("2021-05-22 00:00:00", tz="EST"), 
                       length.out=nrow(unique(Temp)), by='10 min'), '%Y-%m-%d %H:%M:%S')

datetime <- as.POSIXct(datetime, tz="EST")

Gap_filled_weather <- as.data.frame(datetime)

Gap_filled <- merge(Gap_filled_weather, sub_weather, all=TRUE)
Gap_filled <- merge(sub_wave, Gap_filled, all=TRUE)

Gap_filled <- subset(Gap_filled, datetime %in% Gap_filled_weather$datetime)

Gap_filled <- Gap_filled %>% replace_na(list(site = Gap_filled[2,2], 
                                    latitude = Gap_filled[2,3],
                                    longitude = Gap_filled[2,4])) 

###########################################
colnumbers <- c(3,4,5,8,11,12,13,15,18,19,20)
###########################################

 for (i in colnumbers) {  
  
#interpolate y value based on x list of empty datetimes
new_y = approx(Gap_filled$datetime, Gap_filled[,i] , xout=datetime)

new_y <- new_y %>% as.data.frame(col.names=c("datetime", paste0(i)))

assign((paste0("output_",i)), new_y) 

}
##############################################
df_list = list(output_3, output_4, output_5, output_8, output_11, output_12, output_13, output_15, output_18,output_19,output_20) 
###############################################

Weather_out <- df_list %>% reduce(inner_join, by='datetime')
names(Weather_out) <- c("datetime","wnd.dir","wvht","fetch","rh","wnd","pres", "sw","swnet","lwnet","cloud","airt")

Weather_out <- Weather_out %>% mutate(site = s)

assign((paste0(s,"_Weather")), Weather_out)
}

Weather_2021 <- rbind(CHRP1_Weather,
                      CHRP2_Weather, 
                      CHRP3_Weather,
                      CHRP4_Weather,
                   CHRP5_Weather, #exlcude from 2018
               #    CHRP7_Weather, #exclude from 2021
                      CHRP8_Weather,
                      CHRP9_Weather) 

Weather_2021 <- drop_na(Weather_2021)

write.table(Weather_2021, "C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/Weather_2021.txt", 
            sep=",",row.names=FALSE)



##############################################################################################
#Read in temperature file
FilePath <- ("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/CHRP9_2017_temp.txt")
Temp <- load.ts(FilePath, tz = "EST")
head(Temp)
#Fill gaps in temp data
#Temp$wtr_1.0 <- na.approx(Temp$wtr_1.0, rule=2, na.rm=FALSE)
#Temp$wtr_2.1 <- na.approx(Temp$wtr_2.1, rule=2, na.rm=FALSE)
#Temp$wtr_3.0 <- na.approx(Temp$wtr_3.0, rule=2, na.rm=FALSE)
#Temp$wtr_4.1 <- na.approx(Temp$wtr_4.1, rule=2, na.rm=FALSE)
Temp$wtr_5.1 <- na.approx(Temp$wtr_5.1, rule=2, na.rm=FALSE)
#Temp$wtr_6.0 <- na.approx(Temp$wtr_6.0, rule=2, na.rm=FALSE)
Temp$wtr_7.1 <- na.approx(Temp$wtr_7.1, rule=2, na.rm=FALSE)
Temp$wtr_8.0 <- na.approx(Temp$wtr_8.0, rule=2, na.rm=FALSE)
#Temp$wtr_9.1 <- na.approx(Temp$wtr_9.1, rule=2, na.rm=FALSE)
#Temp$wtr_9.5 <- na.approx(Temp$wtr_9.5, rule=2, na.rm=FALSE)
Temp$wtr_10.0 <- na.approx(Temp$wtr_10.0, rule=2, na.rm=FALSE)
#Temp$wtr_10.5 <- na.approx(Temp$wtr_10.5, rule=2, na.rm=FALSE)
Temp$wtr_11.5 <- na.approx(Temp$wtr_11.5, rule=2, na.rm=FALSE)
#Temp$wtr_12.0 <- na.approx(Temp$wtr_12.0, rule=2, na.rm=FALSE)
#Temp$wtr_13.5 <- na.approx(Temp$wtr_13.5, rule=2, na.rm=FALSE)
#Temp$wtr_14.9 <- na.approx(Temp$wtr_14.9, rule=2, na.rm=FALSE)
Temp$wtr_15.0 <- na.approx(Temp$wtr_15.0, rule=2, na.rm=FALSE)
Temp$wtr_16.0 <- na.approx(Temp$wtr_16.0, rule=2, na.rm=FALSE)
Temp$wtr_17.0 <- na.approx(Temp$wtr_17.0, rule=2, na.rm=FALSE)
Temp$wtr_18.0 <- na.approx(Temp$wtr_18.0, rule=2, na.rm=FALSE)
#Temp$wtr_19.6 <- na.approx(Temp$wtr_19.6, rule=2, na.rm=FALSE)
Temp$wtr_20.4 <- na.approx(Temp$wtr_20.4, rule=2, na.rm=FALSE)
##Temp$wtr_21.3 <- na.approx(Temp$wtr_21.3, rule=2, na.rm=FALSE)
#Temp$wtr_22.0 <- na.approx(Temp$wtr_22.0, rule=2, na.rm=FALSE)
#Temp$wtr_23.0 <- na.approx(Temp$wtr_23.0, rule=2, na.rm=FALSE)
#Temp$wtr_23.8 <- na.approx(Temp$wtr_23.8, rule=2, na.rm=FALSE)

write.table(Temp, "C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/CHRP9_2017_temp.txt", 
            sep=",",row.names=FALSE)