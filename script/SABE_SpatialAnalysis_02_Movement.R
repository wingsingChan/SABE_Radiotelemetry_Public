# Data Analysis - Movement ----
## Load library ----
library(tidyverse)
library(car)
library(emmeans)
library(PerformanceAnalytics)
library(usdm)
library(lme4)
library(glmmTMB)
library(lmerTest)
library(MuMIn)
library(DHARMa)
library(ggplot2)
library(ggpubr)
library(ggeffects)
library(scales)
library(broom.mixed)

## Load Data ----
SABE.MCPresults <- read.csv("data/SABE_movement_mcp.csv", header = TRUE)
SABE.MCPresults$sex <- as.factor(SABE.MCPresults$sex)
SABE.seasonMCPresults <- read.csv("data/SABE_movement_mcp_season.csv", header = TRUE)
SABE.seasonMCPresults$sex <- as.factor(SABE.seasonMCPresults$sex)
SABE.seasonMCPresults$season <- factor(SABE.seasonMCPresults$season, 
                                       levels = c("2020 Wet Season", "2020 Mating Season", 
                                                  "2020 Dry Season", "2021 Nesting Season", 
                                                  "2021 Wet Season", "2021 Mating Season"))
SABE.seasonMCPresults$season.4 <- sub("\\d+\\s", "", SABE.seasonMCPresults$season)
SABE.seasonMCPresults$season.4 <- factor(SABE.seasonMCPresults$season.4, 
                                         levels = c("Nesting Season", "Wet Season", 
                                                    "Mating Season", "Dry Season"))

SABE.trajDf <- read.csv("data/SABE_movement_dist.csv", header = TRUE)
SABE.trajDf$sex <- as.factor(SABE.trajDf$sex)
SABE.trajDf$season <- factor(SABE.trajDf$season, 
                             levels = c("2020 Wet Season", "2020 Mating Season", 
                                        "2020 Dry Season", "2021 Nesting Season", 
                                        "2021 Wet Season", "2021 Mating Season", 
                                        "2021 Dry Season"))
SABE.trajDf$season.4 <- sub("\\d+\\s", "", SABE.trajDf$season)
SABE.trajDf$season.4 <- factor(SABE.trajDf$season.4, levels = c("Nesting Season", "Wet Season", 
                                                                "Mating Season", "Dry Season"))

bio <- SABE.trajDf %>% 
  dplyr::select(c(turtle.id, sex, cl, pl, wt)) %>% 
  distinct()

## Summary Statistics ----
### Summary Table for MCP ----
## Population Average
SABE.MCPresults %>%
  subset(MCP_Lv == 100 & Pts > 5) %>% 
  dplyr::summarise(meanMCP = mean(Area_m2),
                   sdMCP = sd(Area_m2))

## Sexes
SABE.MCPsummary <- SABE.MCPresults %>%
  subset(MCP_Lv == 100 & Pts > 5) %>% 
  dplyr::group_by(sex) %>%
  dplyr::summarise(meanMCP = mean(Area_m2),
                   sdMCP = sd(Area_m2),
                   varMCP = var(Area_m2), 
                   n = n()) %>% 
  mutate(se = sdMCP / sqrt(n), 
         lwr = meanMCP - qt(1 - (0.05/2), n - 1) * se, 
         upr = meanMCP + qt(1 - (0.05/2), n - 1) * se)

## Seasons
SABE.seasonMCPsummary <- SABE.seasonMCPresults %>% 
  subset(MCP_Lv == 100 & Pts > 5) %>% 
  dplyr::group_by(season.4) %>%
  dplyr::summarise(meanMCP = mean(Area_m2),
                   sdMCP = sd(Area_m2),
                   varMCP = var(Area_m2),
                   n = n()) %>% 
  mutate(se = sdMCP / sqrt(n), 
         lwr = meanMCP - qt(1 - (0.05/2), n - 1) * se, 
         upr = meanMCP + qt(1 - (0.05/2), n - 1) * se)

### Summary Table for Displacement Distance ----
## Population Average
SABE.trajDf %>% 
  summarise(AvgDist = mean(dailyDist),
            SdDist = sd(dailyDist))

## Seasons
SABE.seasonDailyDistResults <- SABE.trajDf %>% 
  group_by(season.4) %>%  
  summarise(AvgDist = mean(dailyDist),
            SdDist = sd(dailyDist),
            q1 = quantile(dailyDist, .25),
            q3 = quantile(dailyDist, .75),
            median = median(dailyDist), 
            MaxDist = max(dailyDist), 
            n = n()) %>% 
  mutate(se = SdDist / sqrt(n), 
         lwr = AvgDist - qt(1 - (0.05/2), n - 1) * se, 
         upr = AvgDist + qt(1 - (0.05/2), n - 1) * se)

## Sexes & Seasons
SABE.seasonDailyDistResults <- SABE.trajDf %>% 
  group_by(sex, season.4) %>%  
  summarise(AvgDist = mean(dailyDist),
            SdDist = sd(dailyDist),
            q1 = quantile(dailyDist, .25),
            q3 = quantile(dailyDist, .75),
            median = median(dailyDist), 
            MaxDist = max(dailyDist), 
            n = n()) %>% 
  mutate(se = SdDist / sqrt(n), 
         lwr = AvgDist - qt(1 - (0.05/2), n - 1) * se, 
         upr = AvgDist + qt(1 - (0.05/2), n - 1) * se)

## Statistical Tests ---- 
### Linear Mixed Model: MCP ----
mcpMd.global <- lmer(log(Area_m2 + 1) ~ # Biometric vars
                                         sex + poly(scale(cl), 2) + 
                                        # Seasonal response
                                         season.4 + 
                                        # Reproductive drive
                                         sex * season.4 + 
                                        # Thermoregulative drive
                                         poly(scale(cl), 2) * season.4+ 
                                        # Sampling efforts 
                                         scale(Pts) +
                                        # Individual var.
                                         (1|turtle.id),
                     data = SABE.seasonMCPresults, na.action = "na.fail")
step(mcpMd.global)

mcpMd1 <- lmer(log(Area_m2 + 1) ~ # Biometric vars
                                   sex + 
                                  # Seasonal response
                                   season.4 + 
                                  # Sampling efforts 
                                   scale(Pts) +
                                  # Individual var.
                                   (1|turtle.id),
               data = SABE.seasonMCPresults, na.action = "na.fail")
summary(mcpMd1)
Anova(mcpMd1, type = 2)
emmeans(mcpMd1, list(pairwise ~ sex), adjust = "tukey")
emmeans(mcpMd1, list(pairwise ~ season.4), adjust = "tukey")

qqPlot(resid(mcpMd1))
plot(mcpMd1)

mcpMd1.sim <- simulateResiduals(mcpMd1)
plot(mcpMd1.sim)

### Plotting -- MCP vs Sex & Season 
### Figure 2
ggplot(SABE.seasonMCPresults, aes(x = sex, y = Area_m2, fill = sex, col = sex)) + 
  facet_grid(~season.4) + 
  geom_boxplot(outlier.shape = NA, alpha = .2) + 
  geom_point(alpha = .6) +
  scale_x_discrete(labels = c("Female", "Juvenile", "Male")) + 
  ylab(expression(paste("Home range size  ", (m^2)))) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12))



### Linear Mixed Model: Displacement Distance ----
SABE.seasonDistResults <- SABE.trajDf %>%
  group_by(turtle.id, season) %>%
  summarise(AvgDist = mean(dailyDist),
            sd = sd(dailyDist),
            var = var(dailyDist),
            n = n(),
            TotDist = sum(dailyDist),
            Median = median(dailyDist)) %>%
  filter(n > 4) %>%
  mutate(se = sd / sqrt(n),
         lwr = AvgDist - qt(1 - (0.05/2), n - 1) * se,
         upr = AvgDist + qt(1 - (0.05/2), n - 1) * se)

SABE.seasonDistResults$year.2 <- sub("\\s.*", "", SABE.seasonDistResults$season)
SABE.seasonDistResults$year.2 <- as.factor(SABE.seasonDistResults$year.2)
SABE.seasonDistResults$season.4 <- sub("\\d+\\s", "", SABE.seasonDistResults$season)
SABE.seasonDistResults$season.4 <- factor(SABE.seasonDistResults$season.4, c("Nesting Season", "Wet Season", "Mating Season", "Dry Season"))
SABE.seasonDistResults <- merge(SABE.seasonDistResults, bio,
                                by.x = c("turtle.id"), by.y = c("turtle.id"))

distMd.global <- lmer(AvgDist ~ # Biometric vars
                                 sex + poly(scale(cl),2) + 
                                # Seasonal response
                                 season.4 + 
                                # Reproductive drive
                                 sex * season.4 + 
                                # Thermoregulative drive
                                  poly(scale(cl),2) * season.4 + 
                                # Individual var.
                                  (1|turtle.id),
                      data = SABE.seasonDistResults, na.action = "na.fail")
step(distMd.global)

distMd1 <- lmer(AvgDist ~ # Biometric vars
                          sex + 
                          # Seasonal response
                          season.4 + 
                          # Reproductive drive
                          sex * season.4 + 
                          # Individual var.
                          (1|turtle.id),
                data = SABE.seasonDistResults, na.action = "na.fail")
summary(distMd1)
Anova(distMd1, type = 3)

distMdEmm.season <- emmeans(distMd1, ~ season.4, adjust = "tukey")
summary(distMdEmm.season, type = "response")
pairs(distMdEmm.season)

distMdEmm.sex <- emmeans(distMd1, ~ sex | season.4, adjust = "tukey")
summary(distMdEmm.sex, type = "response")
pairs(distMdEmm.sex)

distMd1.sim <- simulateResiduals(distMd1)
plot(distMd1.sim)

### Plotting -- Displacement Distance & Sex | Season ----
### Figure 1
ggplot(SABE.seasonDailyDistResults, aes(x = season.4, y = AvgDist, col = sex)) +
  geom_point(aes(pch =  sex), size = 2, position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = .1, position = position_dodge(.6)) + 
  ylab("Measured daily displacement distance (m)") + 
  theme_classic() +  
  theme(legend.position = c(.075,.9)) + 
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank()) + 
  scale_color_discrete(labels = c("Female", "Juvenile", "Male")) +
  scale_fill_discrete(labels = c("Female", "Juvenile", "Male")) +
  scale_shape(labels = c("Female", "Juvenile", "Male")) + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 12))

  

