library(tidyverse)

#Wave height model
CHRP1_2017_WVHT <- read.table("Wvht_OLS_CHRP1_2017_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP1_2018_WVHT <- read.table("Wvht_OLS_CHRP1_2018_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP1_2019_WVHT <- read.table("Wvht_OLS_CHRP1_2019_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP2_2017_WVHT <- read.table("Wvht_OLS_CHRP2_2017_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP2_2018_WVHT <- read.table("Wvht_OLS_CHRP2_2018_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP2_2019_WVHT <- read.table("Wvht_OLS_CHRP2_2019_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP3_2017_WVHT <- read.table("Wvht_OLS_CHRP3_2017_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP3_2018_WVHT <- read.table("Wvht_OLS_CHRP3_2018_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP3_2019_WVHT <- read.table("Wvht_OLS_CHRP3_2019_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP4_2017_WVHT <- read.table("Wvht_OLS_CHRP4_2017_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP4_2018_WVHT <- read.table("Wvht_OLS_CHRP4_2018_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP4_2019_WVHT <- read.table("Wvht_OLS_CHRP4_2019_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP5_2017_WVHT <- read.table("Wvht_OLS_CHRP5_2017_singledepth.txt") %>% mutate(site = "CHRP5")
CHRP5_2019_WVHT <- read.table("Wvht_OLS_CHRP5_2019_singledepth.txt") %>% mutate(site = "CHRP5")
CHRP7_2017_WVHT <- read.table("Wvht_OLS_CHRP7_2017_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP7_2018_WVHT <- read.table("Wvht_OLS_CHRP7_2018_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP7_2019_WVHT <- read.table("Wvht_OLS_CHRP7_2019_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP8_2017_WVHT <- read.table("Wvht_OLS_CHRP8_2017_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP8_2018_WVHT <- read.table("Wvht_OLS_CHRP8_2018_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP8_2019_WVHT <- read.table("Wvht_OLS_CHRP8_2019_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP9_2017_WVHT <- read.table("Wvht_OLS_CHRP9_2017_singledepth.txt") %>% mutate(site = "CHRP9")
CHRP9_2018_WVHT <- read.table("Wvht_OLS_CHRP9_2018_singledepth.txt") %>% mutate(site = "CHRP9")
CHRP9_2019_WVHT <- read.table("Wvht_OLS_CHRP9_2019_singledepth.txt") %>% mutate(site = "CHRP9")

Waveheight <- rbind(CHRP1_2017_WVHT, CHRP1_2018_WVHT, CHRP1_2019_WVHT,
                    CHRP2_2017_WVHT, CHRP2_2018_WVHT, CHRP2_2019_WVHT,
                    CHRP3_2017_WVHT, CHRP3_2018_WVHT, CHRP3_2019_WVHT,
                    CHRP4_2017_WVHT, CHRP4_2018_WVHT, CHRP4_2019_WVHT,
                    CHRP5_2017_WVHT, CHRP5_2019_WVHT,
                    CHRP7_2017_WVHT, CHRP7_2018_WVHT, CHRP7_2019_WVHT,
                    CHRP8_2017_WVHT, CHRP8_2018_WVHT, CHRP8_2019_WVHT,
                    CHRP9_2017_WVHT, CHRP9_2018_WVHT, CHRP9_2019_WVHT)

Waveheight <- Waveheight %>% mutate(model = "Wvht")

#Wave height model
CHRP1_2017_WVHT_ND <- read.table("Wvht_OLS_CHRP1_2017_neardeep.txt") %>% mutate(site = "CHRP1")
CHRP1_2018_WVHT_ND <- read.table("Wvht_OLS_CHRP1_2018_neardeep.txt") %>% mutate(site = "CHRP1")
CHRP1_2019_WVHT_ND <- read.table("Wvht_OLS_CHRP1_2019_neardeep.txt") %>% mutate(site = "CHRP1")
CHRP3_2017_WVHT_ND <- read.table("Wvht_OLS_CHRP3_2017_neardeep.txt") %>% mutate(site = "CHRP3")
CHRP3_2018_WVHT_ND <- read.table("Wvht_OLS_CHRP3_2018_neardeep.txt") %>% mutate(site = "CHRP3")
CHRP3_2019_WVHT_ND <- read.table("Wvht_OLS_CHRP3_2019_neardeep.txt") %>% mutate(site = "CHRP3")
CHRP5_2017_WVHT_ND <- read.table("Wvht_OLS_CHRP5_2017_neardeep.txt") %>% mutate(site = "CHRP5")
CHRP5_2019_WVHT_ND <- read.table("Wvht_OLS_CHRP5_2019_neardeep.txt") %>% mutate(site = "CHRP5")
CHRP7_2017_WVHT_ND <- read.table("Wvht_OLS_CHRP7_2017_neardeep.txt") %>% mutate(site = "CHRP7")
CHRP7_2018_WVHT_ND <- read.table("Wvht_OLS_CHRP7_2018_neardeep.txt") %>% mutate(site = "CHRP7")
CHRP7_2019_WVHT_ND <- read.table("Wvht_OLS_CHRP7_2019_neardeep.txt") %>% mutate(site = "CHRP7")
CHRP8_2017_WVHT_ND <- read.table("Wvht_OLS_CHRP8_2017_neardeep.txt") %>% mutate(site = "CHRP8")
CHRP8_2018_WVHT_ND <- read.table("Wvht_OLS_CHRP8_2018_neardeep.txt") %>% mutate(site = "CHRP8")
CHRP8_2019_WVHT_ND <- read.table("Wvht_OLS_CHRP8_2019_neardeep.txt") %>% mutate(site = "CHRP8")


Waveheight_Deep <- rbind(CHRP1_2017_WVHT_ND, CHRP1_2018_WVHT_ND, CHRP1_2019_WVHT_ND,
                      CHRP3_2017_WVHT_ND, CHRP3_2018_WVHT_ND, CHRP3_2019_WVHT_ND,
                           CHRP5_2017_WVHT_ND, CHRP5_2019_WVHT_ND,
                    CHRP7_2017_WVHT_ND, CHRP7_2018_WVHT_ND, CHRP7_2019_WVHT_ND,
                    CHRP8_2017_WVHT_ND, CHRP8_2018_WVHT_ND, CHRP8_2019_WVHT_ND)

Waveheight_Deep <- Waveheight_Deep %>% mutate(model = "Deep") 

#Fetch model
CHRP1_2017_Fetch <- read.table("Fetch_OLS_CHRP1_2017_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP1_2018_Fetch <- read.table("Fetch_OLS_CHRP1_2018_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP1_2019_Fetch <- read.table("Fetch_OLS_CHRP1_2019_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP2_2017_Fetch <- read.table("Fetch_OLS_CHRP2_2017_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP2_2018_Fetch <- read.table("Fetch_OLS_CHRP2_2018_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP2_2019_Fetch <- read.table("Fetch_OLS_CHRP2_2019_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP3_2017_Fetch <- read.table("Fetch_OLS_CHRP3_2017_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP3_2018_Fetch <- read.table("Fetch_OLS_CHRP3_2018_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP3_2019_Fetch <- read.table("Fetch_OLS_CHRP3_2019_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP4_2017_Fetch <- read.table("Fetch_OLS_CHRP4_2017_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP4_2018_Fetch <- read.table("Fetch_OLS_CHRP4_2018_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP4_2019_Fetch <- read.table("Fetch_OLS_CHRP4_2019_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP5_2017_Fetch <- read.table("Fetch_OLS_CHRP5_2017_singledepth.txt") %>% mutate(site = "CHRP5")
CHRP5_2019_Fetch <- read.table("Fetch_OLS_CHRP5_2019_singledepth.txt") %>% mutate(site = "CHRP5")
CHRP7_2017_Fetch <- read.table("Fetch_OLS_CHRP7_2017_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP7_2018_Fetch <- read.table("Fetch_OLS_CHRP7_2018_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP7_2019_Fetch <- read.table("Fetch_OLS_CHRP7_2019_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP8_2017_Fetch <- read.table("Fetch_OLS_CHRP8_2017_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP8_2018_Fetch <- read.table("Fetch_OLS_CHRP8_2018_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP8_2019_Fetch <- read.table("Fetch_OLS_CHRP8_2019_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP9_2017_Fetch <- read.table("Fetch_OLS_CHRP9_2017_singledepth.txt") %>% mutate(site = "CHRP9")
CHRP9_2018_Fetch <- read.table("Fetch_OLS_CHRP9_2018_singledepth.txt") %>% mutate(site = "CHRP9")
CHRP9_2019_Fetch <- read.table("Fetch_OLS_CHRP9_2019_singledepth.txt") %>% mutate(site = "CHRP9")

Fetch <- rbind(CHRP1_2017_Fetch, CHRP1_2018_Fetch, CHRP1_2019_Fetch,
               CHRP2_2017_Fetch, CHRP2_2018_Fetch, CHRP2_2019_Fetch,
               CHRP3_2017_Fetch, CHRP3_2018_Fetch, CHRP3_2019_Fetch,
               CHRP4_2017_Fetch, CHRP4_2018_Fetch, CHRP4_2019_Fetch,
               CHRP5_2017_Fetch, CHRP5_2019_Fetch,
               CHRP7_2017_Fetch, CHRP7_2018_Fetch, CHRP7_2019_Fetch,
               CHRP8_2017_Fetch, CHRP8_2018_Fetch, CHRP8_2019_Fetch,
               CHRP9_2017_Fetch, CHRP9_2018_Fetch, CHRP9_2019_Fetch)

Fetch <- Fetch %>% mutate(model = "Fetch")


#Read model
CHRP1_2017_Read <- read.table("Read_OLS_CHRP1_2017_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP1_2018_Read <- read.table("Read_OLS_CHRP1_2018_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP1_2019_Read <- read.table("Read_OLS_CHRP1_2019_singledepth.txt") %>% mutate(site = "CHRP1")
CHRP2_2017_Read <- read.table("Read_OLS_CHRP2_2017_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP2_2018_Read <- read.table("Read_OLS_CHRP2_2018_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP2_2019_Read <- read.table("Read_OLS_CHRP2_2019_singledepth.txt") %>% mutate(site = "CHRP2")
CHRP3_2017_Read <- read.table("Read_OLS_CHRP3_2017_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP3_2018_Read <- read.table("Read_OLS_CHRP3_2018_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP3_2019_Read <- read.table("Read_OLS_CHRP3_2019_singledepth.txt") %>% mutate(site = "CHRP3")
CHRP4_2017_Read <- read.table("Read_OLS_CHRP4_2017_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP4_2018_Read <- read.table("Read_OLS_CHRP4_2018_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP4_2019_Read <- read.table("Read_OLS_CHRP4_2019_singledepth.txt") %>% mutate(site = "CHRP4")
CHRP5_2017_Read <- read.table("Read_OLS_CHRP5_2017_singledepth.txt") %>% mutate(site = "CHRP5")
CHRP5_2019_Read <- read.table("Read_OLS_CHRP5_2019_singledepth.txt") %>% mutate(site = "CHRP5")
CHRP7_2017_Read <- read.table("Read_OLS_CHRP7_2017_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP7_2018_Read <- read.table("Read_OLS_CHRP7_2018_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP7_2019_Read <- read.table("Read_OLS_CHRP7_2019_singledepth.txt") %>% mutate(site = "CHRP7")
CHRP8_2017_Read <- read.table("Read_OLS_CHRP8_2017_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP8_2018_Read <- read.table("Read_OLS_CHRP8_2018_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP8_2019_Read <- read.table("Read_OLS_CHRP8_2019_singledepth.txt") %>% mutate(site = "CHRP8")
CHRP9_2017_Read <- read.table("Read_OLS_CHRP9_2017_singledepth.txt") %>% mutate(site = "CHRP9")
CHRP9_2018_Read <- read.table("Read_OLS_CHRP9_2018_singledepth.txt") %>% mutate(site = "CHRP9")
CHRP9_2019_Read <- read.table("Read_OLS_CHRP9_2019_singledepth.txt") %>% mutate(site = "CHRP9")

Read <- rbind(CHRP1_2017_Read, CHRP1_2018_Read, CHRP1_2019_Read,
              CHRP2_2017_Read, CHRP2_2018_Read, CHRP2_2019_Read,
              CHRP3_2017_Read, CHRP3_2018_Read, CHRP3_2019_Read,
              CHRP4_2017_Read, CHRP4_2018_Read, CHRP4_2019_Read,
              CHRP5_2017_Read, CHRP5_2019_Read,
              CHRP7_2017_Read, CHRP7_2018_Read, CHRP7_2019_Read,
              CHRP8_2017_Read, CHRP8_2018_Read, CHRP8_2019_Read,
              CHRP9_2017_Read, CHRP9_2018_Read, CHRP9_2019_Read)

Read <- Read %>% mutate(model = "Read")

##################################################################
#Combine into a single OLS dataframe
Metabolism <- rbind(Waveheight, Fetch, Read)

Deep <- rbind(Waveheight, Waveheight_Deep)

#Change column names
colnames(Metabolism)[1] ="date"
colnames(Metabolism)[2] ="GPP"
colnames(Metabolism)[3] ="R"
colnames(Metabolism)[4] ="NEP"
colnames(Metabolism)[5] = "Residuals"

#Create a OLS .csv file for future analysis
write.csv(Metabolism, file="Metabolism.csv", row.names=FALSE)

colnames(Deep)[1] ="date"
colnames(Deep)[2] ="GPP"
colnames(Deep)[3] ="R"
colnames(Deep)[4] ="NEP"
colnames(Deep)[5] = "Residuals"

#Create a OLS .csv file for future analysis
write.csv(Depth, file="Deep.csv", row.names=FALSE)
