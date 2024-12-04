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
                                         sex + scale(cl) + scale(I(cl^2)) + 
                                        # Seasonal response
                                         season.4 + 
                                        # Reproductive drive
                                         sex * season.4 + 
                                        # Thermoregulative drive
                                         scale(cl) * season.4 + 
                                         scale(I(cl^2)) * season.4+ 
                                        # Individual var.
                                         (1|turtle.id),
                     data = SABE.seasonMCPresults, na.action = "na.fail")
step(mcpMd.global)

mcpMd1 <- lmer(log(Area_m2 + 1) ~ # Biometric vars
                                   sex + scale(cl) + scale(I(cl^2)) + 
                                  # Seasonal response
                                   season.4 + 
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

### Plotting -- MCP vs Sex ----
#### Figure S5. 
ggplot(SABE.seasonMCPresults, aes(x = sex, y = Area_m2, fill = sex)) + 
  geom_boxplot(outlier.shape = NA, width = .5) + 
  geom_jitter(width = .1, alpha = .3) +
  scale_x_discrete(labels = c("Female", "Juvenile", "Male")) + 
  ylab(expression(paste("Home range size  ", (m^2)))) + 
  geom_signif(comparisons = list(c("M", "J"), c("M", "F")), 
              map_signif_level = TRUE,
              y_position = c(34000, 36000),
              annotation = c("*", "**"), 
              tip_length = 0.01, 
              col = "black") + 
  scale_y_continuous(breaks = seq(0, 30000, by = 10000)) +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + 
  theme(legend.position = "none")

### Plotting -- MCP vs Season 
### Figure 3
ggplot(SABE.seasonMCPresults, aes(x = season.4, y = Area_m2, fill = season.4)) + 
  geom_boxplot(outlier.shape = NA, width = .5) + 
  geom_jitter(width = .1, alpha = .3) +
  scale_x_discrete(labels = c("Nesting", "Wet", "Mating", "Dry")) + 
  ylab(expression(paste("Home range size  ", (m^2)))) + 
  geom_signif(comparisons = list(c("Nesting Season", "Wet Season")), 
              map_signif_level = TRUE,
              y_position = c(34000),
              annotation = c("**"), 
              tip_length = 0.01, 
              col = "black") + 
  scale_y_continuous(breaks = seq(0, 30000, by = 10000)) +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + 
  theme(legend.position = "none")

### Plotting -- Predictive Model ----
### Figure 2
maleCA <- bio %>% filter(sex == "M") %>% summarise(seq(floor(min(cl)), floor(max(cl)), .1))
maleCA <- as.vector(t(maleCA))
predMale <- data.frame(sex = "M", season.4 = rep(c("Dry Season", "Nesting Season", "Wet Season", "Mating Season"), length(maleCA)), 
                       cl = maleCA)
femaleCA <- bio %>% filter(sex == "F") %>% summarise(seq(floor(min(cl)), floor(max(cl)), .1))
femaleCA <- as.vector(t(femaleCA))
predFemale <- data.frame(sex = "F", season.4 = rep(c("Dry Season", "Nesting Season", "Wet Season", "Mating Season"), length(femaleCA)), 
                         cl = femaleCA)
juvenileCA <- bio %>% filter(sex == "J") %>% summarise(seq(floor(min(cl)), floor(max(cl)), .1))
juvenileCA <- as.vector(t(juvenileCA))
predJuvenile <- data.frame(sex = "J", season.4 = rep(c("Dry Season", "Nesting Season", "Wet Season", "Mating Season"), length(juvenileCA)), 
                           cl = juvenileCA)
predDf <- rbind(predFemale, predJuvenile, predMale)
predDf$sex <- factor(predDf$sex, c("F", "J", "M"))
predDf$season.4 <- factor(predDf$season.4, c("Dry Season", "Nesting Season", "Wet Season", "Mating Season"))
str(predDf)

easyPredCI <- function(model, newdata, alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model, newdata=newdata, re.form=NA)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2], newdata)
  beta <- fixef(model) ## fixed-effects coefficients 
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## construct 95% Normal CIs on the link scale and
  ## transform back to the response (probability) scale:
  ## crit <- -qnorm(alpha/2)
  (cbind(fit = pred0,
         lwr=pred0-1*pred.se,
         upr=pred0+1*pred.se))
}
predMCP <- cbind(predDf, easyPredCI(mcpMd1, predDf))

SABE.seasonMCPresults$response <- log(SABE.seasonMCPresults$Area_m2 + 1)
ggplot(predMCP) + 
  facet_wrap(~season.4, nrow = 1) + 
  geom_line(aes(cl, fit, fill = factor(sex), col = factor(sex)), size = .4) +
  geom_ribbon(aes(cl, ymin = lwr, ymax = upr, fill = factor(sex)), alpha = .2) +
  geom_point(data = SABE.seasonMCPresults, 
             aes(cl, y = response, color = factor(sex), shape = factor(sex)), 
             size = 2, alpha = .8) + 
  theme_bw() +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.grid = element_blank(), 
        legend.title = element_blank()) +
  xlab("Carapace length (mm)") + 
  ylab(expression(paste("Predicted values of log(Home Range + 1) ", (m^2)))) + 
  coord_cartesian(ylim=c(0, 12))



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
                                 sex + scale(cl) + scale(I(cl^2)) + 
                                # Seasonal response
                                 season.4 + 
                                # Reproductive drive
                                 sex * season.4 + 
                                # Thermoregulative drive
                                 (scale(cl) + scale(I(cl^2))) * season.4 + 
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
  theme(legend.position = c(.05,.65)) + 
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank())

