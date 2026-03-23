library(ggpubr)
library(egg)
library(lubridate)
library(viridis)
library(tidyverse)
library(ggfortify)
library(lme4)
library(lmerTest)
library(report)
library(emmeans)

###########################################################################
Res <- read.csv("Deep.csv")

Res$Date <- as.POSIXct(Res$Date, tz="EST")
Res <- Res %>% mutate(my = (format(as.Date(Date), "%Y-%m")))
Res <- Res %>% mutate(md = (format(as.Date(Date), "%m-%d")))  

Res$Year <- format(Res$Date, "%Y")

Res <- Res %>% filter(site!="CHRP2", site!="CHRP4", site!="CHRP9")

CHRP3 <- Res %>% filter(site=="CHRP3")

CHRP8 <- Res %>% filter(site=="CHRP8")


#########################################################################

Res_wide <- read.csv("Deep_untidy.csv")

Res_wide$Date <- as.POSIXct(Res_wide$Date, tz="EST")
Res_wide <- Res_wide %>% mutate(my = (format(as.Date(Date), "%Y-%m")))

Res_wide$Year <- format(Res_wide$Date, "%Y")

Res_wide <- Res_wide %>% filter(site!="CHRP2", site!="CHRP4", site!="CHRP9")

Res_long <- Res_wide %>%
  pivot_longer(
    cols = c(Shallow_GPP, Deep_GPP, Shallow_R, Deep_R, Shallow_NEP, Deep_NEP),
    names_to = c("model", ".value"),
    names_pattern = "(Shallow|Deep)_(.*)"
  ) %>%
  mutate(model = as.factor(model),
         site = as.factor(site))
##############################################################################
#Remove nonsense results

Res_long$GPP[Res_long$GPP < 0] <- 0       
Res_long$R[Res_long$R > 0] <- 0

# Mixed effects models to compare two sets of results:
lmer_GPP <- lmer(GPP ~ model + (1|site) + (1|Date), data = Res_long)
summary(lmer_GPP)

lmer_R <- lmer(R ~ model + (1|site) + (1|Date), data = Res_long)
summary(lmer_R)

lmer_NEP <- lmer(NEP ~ model + (1|site) + (1|Date), data = Res_long)
summary(lmer_NEP)

#Plots showing difference between two methods

Sum_Metab_site <- Res_long %>%
  group_by(model, site, my) %>%
  summarize(avg_GPP = mean(GPP),
            se_GPP = sd(GPP, na.rm=TRUE)/sum(!is.na(GPP)),
            Ci_GPP = 1.96*se_GPP,
            avg_R = mean(R),
            se_R = sd(R, na.rm=TRUE)/sum(!is.na(R)),
            Ci_R = 1.96*se_R,
            avg_NEP = mean(NEP),
            se_NEP = sd(NEP, na.rm=TRUE)/sum(!is.na(NEP)),
            Ci_NEP = 1.96*se_NEP)

Sum_Metab <- Res_long %>%
  group_by(model, my) %>%
  summarize(avg_GPP = mean(GPP),
            se_GPP = sd(GPP, na.rm=TRUE)/sum(!is.na(GPP)),
            Ci_GPP = 1.96*se_GPP,
            avg_R = mean(R),
            se_R = sd(R, na.rm=TRUE)/sum(!is.na(R)),
            Ci_R = 1.96*se_R,
            avg_NEP = mean(NEP),
            se_NEP = sd(NEP, na.rm=TRUE)/sum(!is.na(NEP)),
            Ci_NEP = 1.96*se_NEP)

DP_Resp <- ggplot(aes(x = model, y = avg_R, color=site), data=Sum_Metab_site) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1)) +
  geom_hline(yintercept=0) +
  xlab("Model") + ylab("Mean Resp. (mg/O2/L/day)") +
  theme_bw() +
  scale_color_viridis(discrete = TRUE)+
  scale_x_discrete(labels = c("Deep" = "Submerged", "Shallow" = "Full")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  annotate("text", x = 0.5, y = -0.125, size = 6, label="B")


DP_GPP <- ggplot(aes(x = model, y = avg_GPP, color=site), data=Sum_Metab_site) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1)) +
  geom_hline(yintercept=0) +
  xlab("Model") + ylab("Mean GPP (mg/O2/L/day)") +
  theme_bw() +
  scale_color_viridis(discrete = TRUE) +
  scale_x_discrete(labels = c("Deep" = "Submerged", "Shallow" = "Full")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  annotate("text", x = 0.5, y = 1, size = 6, label="A")

DP_NEP <- ggplot(aes(x = model, y = avg_NEP, color=site), data=Sum_Metab_site) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1)) +
  geom_hline(yintercept=0) +
  xlab("Model") + ylab("Mean NEP. (mg/O2/L/day)") +
  theme_bw() +
  scale_color_viridis(discrete = TRUE)+
  scale_x_discrete(labels = c("Deep" = "Submerged", "Shallow" = "Full")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10)) +
annotate("text", x = 0.5, y = 0.5, size = 6, label="C")

ggpubr::ggarrange(DP_GPP, DP_Resp, DP_NEP, ncol=1,
                  common.legend = TRUE, 
                  legend = "right")

#Export 1300 x 650