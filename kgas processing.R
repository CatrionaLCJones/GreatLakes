library(lubridate)
library(tidyverse)
library(xts)
library(viridis)
library(ggfortify)
library(car)

kgas_dm_fetch <- rbind(kgas_dm_fetch_CHRP1_2017, kgas_dm_fetch_CHRP1_2018, 
                       kgas_dm_fetch_CHRP1_2019, kgas_dm_fetch_CHRP2_2017,
                       kgas_dm_fetch_CHRP2_2018, kgas_dm_fetch_CHRP2_2019,
                       kgas_dm_fetch_CHRP3_2017, kgas_dm_fetch_CHRP3_2018, 
                       kgas_dm_fetch_CHRP3_2019, kgas_dm_fetch_CHRP4_2017, 
                       kgas_dm_fetch_CHRP4_2018, kgas_dm_fetch_CHRP4_2019,
                       kgas_dm_fetch_CHRP5_2017, kgas_dm_fetch_CHRP5_2019,
                       kgas_dm_fetch_CHRP7_2017, kgas_dm_fetch_CHRP7_2018, 
                       kgas_dm_fetch_CHRP7_2019, kgas_dm_fetch_CHRP8_2017, 
                       kgas_dm_fetch_CHRP8_2018, kgas_dm_fetch_CHRP8_2019,
                       kgas_dm_fetch_CHRP9_2017, kgas_dm_fetch_CHRP9_2018, 
                       kgas_dm_fetch_CHRP9_2019)

write.csv(kgas_dm_fetch, file="kgas_dm_fetch.csv", row.names=FALSE)

####################################################################################

kgas_dm_wvht <- rbind(kgas_dm_wvht_CHRP1_2017, kgas_dm_wvht_CHRP1_2018, 
                       kgas_dm_wvht_CHRP1_2019, kgas_dm_wvht_CHRP2_2017,
                       kgas_dm_wvht_CHRP2_2018, kgas_dm_wvht_CHRP2_2019,
                       kgas_dm_wvht_CHRP3_2017, kgas_dm_wvht_CHRP3_2018, 
                       kgas_dm_wvht_CHRP3_2019, kgas_dm_wvht_CHRP4_2017, 
                       kgas_dm_wvht_CHRP4_2018, kgas_dm_wvht_CHRP4_2019,
                       kgas_dm_wvht_CHRP5_2017, kgas_dm_wvht_CHRP5_2019,
                       kgas_dm_wvht_CHRP7_2017, kgas_dm_wvht_CHRP7_2018, 
                       kgas_dm_wvht_CHRP7_2019, kgas_dm_wvht_CHRP8_2017, 
                       kgas_dm_wvht_CHRP8_2018, kgas_dm_wvht_CHRP8_2019,
                       kgas_dm_wvht_CHRP9_2017, kgas_dm_wvht_CHRP9_2018, 
                       kgas_dm_wvht_CHRP9_2019)

write.csv(kgas_dm_wvht, file="kgas_dm_wvht.csv", row.names=FALSE)

#################################################################################

kgas_read_soloviev <- rbind(kgas_read_soloviev_CHRP1_2017, kgas_read_soloviev_CHRP1_2018, 
                       kgas_read_soloviev_CHRP1_2019, kgas_read_soloviev_CHRP2_2017,
                       kgas_read_soloviev_CHRP2_2018, kgas_read_soloviev_CHRP2_2019,
                       kgas_read_soloviev_CHRP3_2017, kgas_read_soloviev_CHRP3_2018, 
                       kgas_read_soloviev_CHRP3_2019, kgas_read_soloviev_CHRP4_2017, 
                       kgas_read_soloviev_CHRP4_2018, kgas_read_soloviev_CHRP4_2019,
                       kgas_read_soloviev_CHRP5_2017, kgas_read_soloviev_CHRP5_2019,
                       kgas_read_soloviev_CHRP7_2017, kgas_read_soloviev_CHRP7_2018, 
                       kgas_read_soloviev_CHRP7_2019, kgas_read_soloviev_CHRP8_2017, 
                       kgas_read_soloviev_CHRP8_2018, kgas_read_soloviev_CHRP8_2019,
                       kgas_read_soloviev_CHRP9_2017, kgas_read_soloviev_CHRP9_2018, 
                       kgas_read_soloviev_CHRP9_2019)

write.csv(kgas_read_soloviev, file="kgas_read_soloviev.csv", row.names=FALSE)

##################################################################################

Sol_read <- read.csv("kgas_read_soloviev.csv")
Wvht <- read.csv("kgas_dm_wvht.csv")
Fetch <- read.csv("kgas_dm_fetch.csv")

Wvht <- Wvht %>% mutate(model = "Wvht") %>% rename(kgas=kgas.dm)
Fetch <- Fetch %>% mutate(model = "Fetch") %>% rename(kgas=kgas.dm)
Sol_read <- Sol_read %>% mutate(model = "Read") %>% rename(kgas=kgas_read)


Kgas <- rbind(Sol_read, Wvht, Fetch)

Kgas <- Kgas %>% mutate(Julien = as.POSIXlt(datetime)$yday) 

Kgas_means <- Kgas_clean %>%
  group_by(model, year(datetime)) %>%
  summarize(avg_kgas = mean(kgas),
            se_kgas = sd(kgas, na.rm=TRUE)/sum(!is.na(kgas)),
            Ci_kgas = 1.96*se_kgas,
            var_kgas = var(kgas, na.rm=TRUE))


Kgas_means <- Kgas_means[is.finite(Kgas_means$avg_kgas),]

Kgas_2017 <- Kgas_means %>% filter(year == 2017) 
Kgas_2018 <- Kgas_means %>% filter(year == 2018)
Kgas_2019 <- Kgas_means %>% filter(year == 2019)

ggplot(data=Kgas_2017, aes(x = Julien, y = avg_kgas, color = Model)) +
  geom_line() +
  geom_line(aes(x=Julien, y = avg_wvht)) +
 # geom_point() +
  geom_errorbar(aes(ymin=avg_kgas-Ci_kgas, ymax=avg_kgas+Ci_kgas), 
                width = 0.1)+
  scale_fill_viridis() +
  facet_wrap(~site) +
  labs(y = "Gas flux velocity (m d-1)",
       x = "Julian day") + theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c("#7AD151ff","#FDE725FF","#440154FF")) 

Kgas_clean <- Kgas[is.finite(Kgas$kgas),]

mod <- aov(kgas~model, Kgas_clean)
summary(mod)
Anova(mod)
TukeyHSD(mod)
