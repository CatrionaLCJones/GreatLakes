library(LakeMetabolizer)
library(rLakeAnalyzer)
library(lubridate)
library(tidyverse)
library(data.table)
library(stringr)
library(anytime)
library(zoo)

###############################################################################
# In order to obtain model components, run the following scripts in order:

# 1. Mixing depth 
# 2. Kd (light)
# 3. Gap filling
# 4. Kgas
# 5. Metabolism model (this script)

############################################################################
##################### MAIN METABOLISM MODEL #################################
############################################################################
#Create loop lists
#There is no data for CHRP5 in 2018. Run script once with CHRP5 included in sites
#and 2018 excluded from year, then run again with CHRP5 excluded and 2018 included

site <- c("CHRP1", "CHRP2","CHRP3","CHRP4","CHRP5", "CHRP7","CHRP8","CHRP9")
years <- c("2017","2019")

for (s in site)   {
  for (y in years) {

    FilePath <- (paste0("C:/Users/Katie/Documents/Research/Purdue_post-doc/Lake Erie Metabolism model/Datasets/",s,"_",y,"_DO_tidy.txt"))
    
    #Load the text file in, specifying timezone as EST - the default is GMT
    DO <- load.ts(FilePath, tz = "EST")
    DO <- DO %>% dplyr::select(-c("date","time"))
    par = sw.to.par(get(paste0("kgas_dm_wvht_",s,"_",y)), sw.col='swnet', coeff=2.114)
    
    PAR <- par %>% group_by(datetime) %>% filter(site==s) %>% dplyr::select(datetime, par)
    
    PAR$datetime <- as.POSIXlt(PAR$datetime, tz="EST")
    Kgas <- get(paste0("kgas_dm_wvht_",s,"_",y))
    Kgas$datetime <- as.POSIXlt(Kgas$datetime, tz="EST")
    
    Metabolism <- inner_join(DO, Kgas,
                             by=c('datetime', 'site')) %>%
      inner_join(., PAR, by='datetime') 
    
    colnames(Metabolism)[6] ="do.obs"
    
    Metabolism$do.sat = o2.at.sat.base(Metabolism[,5], altitude=174)
    
    Metabolism <- Metabolism %>% dplyr::select(-site)
    
    Metabolism$datetime <- as.POSIXct(Metabolism$datetime, tz="EST")
    
    Metabolism <- Metabolism %>% mutate(date = (format(as.POSIXct(datetime,
                                                                  format = '%Y-%m-%d %H:%M:%S'),
                                               format = '%Y-%m-%d')))
    


    #Create a datarame to capture the loop output
    output <- data.frame(matrix(ncol=5))
    
    #Give the model a list of dates to loop on
    date_list <- unique(Metabolism$date)

    # Loop through each day.
   for (i in date_list) {
         
      assign(paste0("OLS_",s,"_",y), data.frame())
      
      #Subset the dataframe by date
      subset <- Metabolism %>% filter(date == i)
      
      subset <- subset %>% dplyr::select(-date)
      
      #Gap fill remaining columns
      cols <- 4:ncol(subset)
      
      for (c in cols) {
        tryCatch({
          subset[,c]<- na.approx(subset[,c], rule=2, na.rm=FALSE)
        }, error=function(e){}) }
      
      subset$do.sat<- na.approx(subset$do.sat, rule=2, na.rm=FALSE)
      subset$do.obs<- na.approx(subset[,5], rule=2, na.rm=FALSE)
      
      
      #within the date loop, we want r to ignore dataframes that have
      #less than 144 observations as this means it is an incomplete
      #24 hour period and will stop the script
      
      # Check if the dataframe has 144 rows
      if (nrow(subset) == 144) {
        # Process the dataframe
    
      metabolism_ols <- metab.ols(datetime=subset$datetime, subset$do.obs, subset$do.sat,
                                  subset$kgas.dm,subset$mixing_depth, subset$par, 
                                  subset[,4], error.type="OE" )
      
#Collects the following output: date, GPP, R, NEP, and model residuals
      output[,1] = i
      output[,2] = metabolism_ols$metab[1]
      output[,3] = metabolism_ols$metab[2] 
      output[,4] = metabolism_ols$metab[3]
      output[,5] = metabolism_ols$mod[2]
      
      write.table(output, file=paste0("Wvht_OLS_",s,"_",y,"_singledepth.txt"), col.names=FALSE,     
                  row.names=FALSE,append=TRUE)
      
        } else {
          # Ignore the dataframe if it doesn't have 144 rows
          print(paste("Ignoring", subset, "as it does not have 144 rows"))
        }
  }
}
}
#######################################################################################
################# Tidying the MLE files ################################################
###########################################################################################
OLS <- list.files(pattern = '*.txt')
Loop_list_z <- substr(OLS, 10, 19)

for (z in Loop_list_z) {
  
  data <- read.table(paste0("Wvht_OLS_",z,"_singledepth.txt"))
  
  colnames(data)[1] ="date"
  colnames(data)[2] ="GPP"
  colnames(data)[3] ="R"
  colnames(data)[4] ="NEP"
  colnames(data)[5] = "Residuals"
  
  data <- data %>%
    pivot_longer(
      cols = c("GPP","R","NEP", "Residuals"),
      names_to = "Measure"
    )  
  
  write.table(data, file=paste0("Wvht_Metabolism_",z,"_singledepth.txt"), col.names=TRUE,     
              row.names=FALSE)
  
}


