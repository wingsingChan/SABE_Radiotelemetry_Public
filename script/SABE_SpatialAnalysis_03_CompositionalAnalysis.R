# Data Analysis - Compositional Analysis ---- 
## Load Libraries ----
library(adehabitatHS)
library(tidyverse)
library(reshape2)
library(MVN)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)

checkNormality <- function(data) {
  plot(stats::density(data))
  hist(data)
  print(paste('Difference between mean and trimmed mean = ',
              mean(data) - mean(data, trim=0.05)))
  print(paste('Standard error = ',sqrt(var(data)/length(data))))
  library(moments)
  print(paste('Skewness = ',skewness(data)))
  print(paste('Kurtosis = ',kurtosis(data)))
  print(shapiro.test(data))
  print(ks.test(data,"pnorm", mean=mean(data), sd=sd(data)))
  qqnorm(data); qqline(data)  
}

## Load Data ----
studyarea <- read.csv("data/SABE_habitatProp_studyarea.csv", header = TRUE)
mcp <- read.csv("data/SABE_habitatProp_mcp.csv", header = TRUE)
locs <- read.csv("data/SABE_habitatProp_loc.csv", header = TRUE)

rownames(studyarea) <- studyarea$turtle.id
studyarea$turtle.id <- NULL
rownames(mcp) <- mcp$turtle.id
mcp$turtle.id <- NULL
rownames(locs) <- locs$turtle.id
locs$turtle.id <- NULL

studyareaSeasonDf <- read.csv("data/SABE_habitatProp_studyarea_season.csv", header = TRUE)
mcpSeasonDf <- read.csv("data/SABE_habitatProp_mcp_season.csv", header = TRUE)
locsSeasonDf <- read.csv("data/SABE_habitatProp_loc_season.csv", header = TRUE)

seasonLv <- c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season", 
              "2021 Wet Season")
studyareaSeasonDf$season <- factor(studyareaSeasonDf$season, levels = seasonLv)
studyareaSeason = list()
for(i in 1:length(levels(studyareaSeasonDf$season))){
  studyareaSeason[[i]] <- studyareaSeasonDf %>% filter(season == seasonLv[i]) %>% dplyr::select(-season)
  rownames(studyareaSeason[[i]]) <- studyareaSeason[[i]]$turtle.id
  studyareaSeason[[i]]$turtle.id <- NULL
  i = i + 1
}
mcpSeasonDf$season <- factor(mcpSeasonDf$season, levels = seasonLv)
mcpSeason = list()
for(i in 1:length(levels(mcpSeasonDf$season))){
  mcpSeason[[i]] <- mcpSeasonDf %>% filter(season == seasonLv[i]) %>% dplyr::select(-season)
  rownames(mcpSeason[[i]]) <- mcpSeason[[i]]$turtle.id
  mcpSeason[[i]]$turtle.id <- NULL
  i = i + 1
}
locsSeasonDf$season <- factor(locsSeasonDf$season, levels = seasonLv)
locsSeason = list()
for(i in 1:length(levels(locsSeasonDf$season))){
  locsSeason[[i]] <- locsSeasonDf %>% filter(season == seasonLv[i]) %>% dplyr::select(-season)
  rownames(locsSeason[[i]]) <- locsSeason[[i]]$turtle.id
  locsSeason[[i]]$turtle.id <- NULL
  i = i + 1
}

## Second-order Habitat Selection ----
mcp2 <- mcp
mcp2[mcp2 == "0"] <- 0.01
a <- apply(mcp2, 2, function(x) {log(x/mcp2$Riffle)})
PoolRuna <- a[,1] - a[,3]

std2 <- studyarea
std2[std2 == "0"] <- 0.01
a0 <- apply(std2, 2, function(x) {log(x/std2$Riffle)})
PoolRuna0 <- a0[,1] - a0[,3]

da <- a - a0
da <- as.data.frame(da)
da <- da %>% dplyr::select(-Riffle)
PoolRunda <- PoolRuna - PoolRuna0
da <- cbind(da, "PoolRun" = PoolRunda)

results_a <- mvn(data = da[,-3], mvnTest = "energy")
results_a$multivariateNormality

checkNormality(da[,1])
t.test(da$Pool, mu = 0, alternative = "two.sided")
checkNormality(da[,2])
t.test(da$Run, mu = 0, alternative = "two.sided")
checkNormality(da[,3])
t.test(da$PoolRun, mu = 0, alternative = "two.sided")

habitatSec <- compana(mcp, studyarea, test = "parametric", alpha = 0.05)
habitatSec
habitatSec[["rm"]] # Preferences
habitatSec[["rmv"]] # t values
2*pt(-abs(habitatSec[["rmv"]]), habitatSec[["rmnb"]]-1) # p value

## Second-order Habitat Selection by Seasons ----
## Multivariate normality were violated for 1. Dry Season// 2. Mating Season// 4. Wet Season 2020// 5. Wet Season
habitatSecSeason = list()
daSeason = list()
for(i in 1:length(mcpSeason)){
  print(seasonLv[i])
  
  mcp2 <- mcpSeason[[i]]
  mcp2[mcp2 == "0"] <- 0.01
  a <- apply(mcp2, 2, function(x) {log(x/mcp2$Riffle)})
  PoolRuna <- a[,1] - a[,3]
  
  std2 <- studyareaSeason[[i]]
  std2[std2 == "0"] <- 0.01
  a0 <- apply(std2, 2, function(x) {log(x/std2$Riffle)})
  PoolRuna0 <- a0[,1] - a0[,3]
  
  daSeason[[i]] <- a - a0
  daSeason[[i]] <- as.data.frame(daSeason[[i]])
  daSeason[[i]] <- daSeason[[i]] %>% dplyr::select(-Riffle)
  
  PoolRunda <- PoolRuna - PoolRuna0
  daSeason[[i]] <- cbind(daSeason[[i]], "PoolRun" = PoolRunda)
  
  results_a <- mvn(data = daSeason[[i]][,-3], mvnTest = "energy")
  print(results_a$multivariateNormality)
  
    habitatSecSeason[[i]] <- compana(mcpSeason[[i]], studyareaSeason[[i]], test = "randomisation", nrep = 10000, alpha = 0.05)
    print(habitatSecSeason[[i]])
    print(habitatSecSeason[[i]][["rank"]])

    print(wilcox.test(daSeason[[i]]$Pool, mu = 0, alternative = "two.sided"))
    print(wilcox.test(daSeason[[i]]$Run, mu = 0, alternative = "two.sided"))
    print(wilcox.test(daSeason[[i]]$PoolRun, mu = 0, alternative = "two.sided"))

  i = i + 1
}

## Third-order Habitat Selection ----
locs3 <- locs
locs3[locs3 == "0"] <- 0.01
b <- apply(locs3, 2, function(x) {log(x/locs3$Riffle)})
PoolRunb <- b[,1] - b[,3]

mcp3 <- mcp
mcp3[mcp3 == "0"] <- 0.01
b0 <- apply(mcp3, 2, function(x) {log(x/mcp3$Riffle)})
PoolRunb0 <- b0[,1] - b0[,3]

db <- b - b0
db <- as.data.frame(db)
db <- db %>% dplyr::select(-Riffle)
PoolRundb <- PoolRunb - PoolRunb0
db <- cbind(db, "PoolRun" = PoolRundb)

results_b <- mvn(data = db[,-3], mvnTest = "energy")
results_b$multivariateNormality

checkNormality(log(db[,1]+10))
wilcox.test(db$Pool, mu = 0, alternative = "two.sided")
checkNormality(log(db[,2]+10))
wilcox.test(db$Run, mu = 0, alternative = "two.sided")
checkNormality(log(db[,3]+10))
wilcox.test(db$PoolRun, mu = 0, alternative = "two.sided")

habitatThrd <- compana(locs, mcp, test = "randomisation", nrep = 10000, alpha = 0.05) 
habitatThrd
habitatThrd[["rm"]]
habitatThrd[["rmv"]] # mean difference between the log-ratio of used and available

## Third-order Habitat Selection by Seasons ----
habitatThrdSeason = list()
dbSeason = list()
for(i in 1:length(locsSeason)){
  print(seasonLv[i])
  
  locs3 <- locsSeason[[i]]
  locs3[locs3 == "0"] <- 0.01
  b <- apply(locs3, 2, function(x) {log(x/locs3$Riffle)})
  PoolRunb <- b[,1] - b[,3]
  
  mcp3 <- mcpSeason[[i]]
  mcp3[mcp3 == "0"] <- 0.01
  b0 <- apply(mcp3, 2, function(x) {log(x/mcp3$Riffle)})
  PoolRunb0 <- b0[,1] - b0[,3]
  
  dbSeason[[i]] <- b - b0
  dbSeason[[i]] <- as.data.frame(dbSeason[[i]])
  dbSeason[[i]] <- dbSeason[[i]] %>% dplyr::select(-Riffle)
  PoolRundb <- PoolRunb - PoolRunb0
  dbSeason[[i]] <- cbind(dbSeason[[i]], "PoolRun" = PoolRundb)
  
  results_b <- mvn(data = dbSeason[[i]][,-3], mvnTest = "energy")
  print(results_b$multivariateNormality)
  
    habitatThrdSeason[[i]] <- compana(locsSeason[[i]], mcpSeason[[i]], test = "randomisation", nrep = 10000, alpha = 0.05)
    print(habitatThrdSeason[[i]])
    print(habitatThrdSeason[[i]][["rank"]])
    
    print(wilcox.test(dbSeason[[i]]$Pool, mu = 0, alternative = "two.sided"))
    print(wilcox.test(dbSeason[[i]]$Run, mu = 0, alternative = "two.sided"))
    print(wilcox.test(dbSeason[[i]]$PoolRun, mu = 0, alternative = "two.sided"))

  i = i + 1
}

## Graphical Output ----
daData <- da
daData <- daData %>% 
  dplyr::rename("Pool:Riffle" = Pool,
                "Run:Riffle" = Run,
                "Pool:Run" = PoolRun)
daData <- data.frame(melt(daData), lv = "Second Order")

dbData <- db
dbData <- dbData %>% 
  dplyr::rename("Pool:Riffle" = Pool,
                "Run:Riffle" = Run,
                "Pool:Run" = PoolRun)
dbData <- data.frame(melt(dbData), lv = "Third Order")

dData <- rbind(daData, dbData)

daSeasonData <- do.call(rbind, lapply(1:length(daSeason), function(i) data.frame(season = seasonLv[i], daSeason[[i]])))
daSeasonData <- daSeasonData %>%
  dplyr::rename("Pool:Riffle" = Pool,
                "Run:Riffle" = Run,
                "Pool:Run" = PoolRun)
daSeasonData <- data.frame(melt(daSeasonData, id.vars = "season"), lv = "Second Order")

dbSeasonData <- do.call(rbind, lapply(1:length(dbSeason), function(i) data.frame(season = seasonLv[i], dbSeason[[i]])))
dbSeasonData <- dbSeasonData %>%
  dplyr::rename("Pool:Riffle" = Pool,
                "Run:Riffle" = Run,
                "Pool:Run" = PoolRun)
dbSeasonData <- data.frame(melt(dbSeasonData, id.vars = "season"), lv = "Third Order")

dSeasonData <- rbind(daSeasonData, dbSeasonData)
dSeasonData$season <- factor(dSeasonData$season, levels = seasonLv)

daSummary <- daData %>% dplyr::select(-lv)
daSummary <- daSummary %>%
  group_by(variable) %>%
  dplyr::summarise(mean = mean(value), 
                   sd = sd(value), 
                   n = n()) %>% 
  dplyr::mutate(se = sd / sqrt(n), 
                lwr = mean - qt(1 - (0.05/2), n - 1) * se, 
                upr = mean + qt(1 - (0.05/2), n - 1) * se)
daSummary <- as.data.frame(daSummary)

dbSummary <- dbData %>% dplyr::select(-lv)
dbSummary <- dbSummary %>%
  group_by(variable) %>%
  dplyr::summarise(mean = mean(value), 
                   sd = sd(value), 
                   n = n()) %>% 
  dplyr::mutate(se = sd / sqrt(n), 
                lwr = mean - qt(1 - (0.05/2), n - 1) * se, 
                upr = mean + qt(1 - (0.05/2), n - 1) * se)
dbSummary <- as.data.frame(dbSummary)

daSummary <- data.frame(daSummary, lv = "Second Order")
dbSummary <- data.frame(dbSummary, lv = "Third Order")
dSummary <- rbind(daSummary, dbSummary)

daSeasonSummary <- daSeasonData %>% dplyr::select(-lv)
daSeasonSummary <- daSeasonSummary %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value), 
                   sd = sd(value), 
                   n = n()) %>% 
  dplyr::mutate(se = sd / sqrt(n), 
                lwr = mean - qt(1 - (0.05/2), n - 1) * se, 
                upr = mean + qt(1 - (0.05/2), n - 1) * se)
daSeasonSummary <- as.data.frame(daSeasonSummary)

dbSeasonSummary <- dbSeasonData %>% dplyr::select(-lv)
dbSeasonSummary <- dbSeasonSummary %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value), 
                   sd = sd(value), 
                   n = n()) %>% 
  dplyr::mutate(se = sd / sqrt(n), 
                lwr = mean - qt(1 - (0.05/2), n - 1) * se, 
                upr = mean + qt(1 - (0.05/2), n - 1) * se)
dbSeasonSummary <- as.data.frame(dbSeasonSummary)

daSeasonSummary <- data.frame(daSeasonSummary, lv = "Second Order")
dbSeasonSummary <- data.frame(dbSeasonSummary, lv = "Third Order")
dSeasonSummary <- rbind(daSeasonSummary, dbSeasonSummary)
dSeasonSummary$season <- factor(dSeasonSummary$season, levels = seasonLv)

relativeHS2 <- ggplot(filter(dSummary, lv == "Second Order"), aes(x = variable, y = mean, col = lv)) +
  geom_point(size = 4, position = position_dodge(.2)) +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se), width = 0.1, position = position_dodge(.2)) + 
  geom_jitter(data = filter(dData, lv == "Second Order"), aes(x = variable, y = value, fill = lv), size = 2, alpha = 0.2, position = position_jitter(.05)) + 
  geom_hline(yintercept = 0) + 
  ylab(expression(paste("Relative use of habitat ", italic(" i"), " (", italic("d"["i"]), ")"))) +
  scale_color_manual(values = c("dark red")) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank()) + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 12))
relativeHS2

relativeSeasonHS2 <- ggplot(filter(dSeasonSummary, lv == "Second Order"), aes(x = variable, y = mean, col = lv)) +
  facet_wrap(~ season, nrow = 3) + 
  geom_point(size = 4, position = position_dodge(.2)) +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se), width = 0.05, position = position_dodge(.2)) + 
  geom_jitter(data = filter(dSeasonData, lv == "Second Order"), aes(x = variable, y = value, fill = lv), size = 2, alpha = 0.2, position = position_jitter(.05)) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c("dark red")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "gainsboro")) + 
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(), 
        panel.grid = element_blank(), 
        legend.position = "none") + 
  theme(strip.text = element_text(size = 18), 
        axis.text = element_text(size = 14))
relativeSeasonHS2

relativeHS3 <- ggplot(filter(dSummary, lv == "Third Order"), aes(x = variable, y = mean, col = lv)) +
  geom_point(size = 4, position = position_dodge(.2)) +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se), width = 0.1, position = position_dodge(.2)) + 
  geom_jitter(data = filter(dData, lv == "Third Order"), aes(x = variable, y = value, fill = lv), size = 2, alpha = 0.2, position = position_jitter(.05)) + 
  geom_hline(yintercept = 0) + 
  ylab(expression(paste("Relative use of habitat ", italic(" i"), " (", italic("d"["i"]), ")"))) +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)) +
  scale_color_manual(values = c("dark blue")) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank()) + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 12))
relativeHS3

relativeSeasonHS3 <- ggplot(filter(dSeasonSummary, lv == "Third Order"), aes(x = variable, y = mean, col = lv)) +
  facet_wrap(~ season, nrow = 3) + 
  geom_point(size = 4, position = position_dodge(.2)) +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se), width = 0.05, position = position_dodge(.2)) + 
  geom_jitter(data = filter(dSeasonData, lv == "Third Order"), aes(x = variable, y = value, fill = lv), size = 2, alpha = 0.2, position = position_jitter(.05)) + 
  geom_hline(yintercept = 0) + 
  scale_y_continuous(breaks = seq(-10, 10, by = 2)) +
  scale_color_manual(values = c("dark blue")) + 
  theme_bw() +
  theme(panel.background = element_rect(fill = "gainsboro")) + 
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(), 
        panel.grid = element_blank(), 
        legend.position = "none") + 
  theme(strip.text = element_text(size = 18), 
        axis.text = element_text(size = 14))
relativeSeasonHS3

## Figure 4
relativeHS2 <- relativeHS2 + 
  theme(axis.title.y = element_text(angle = 90)) + 
  labs(title = "(a) Second-order selection") + 
  theme(plot.title = element_text(size = 12))
relativeHS3 <- relativeHS3 + 
  theme(axis.title.y = element_text(angle = 90)) + 
  labs(title = "(b) Third-order selection") + 
  theme(plot.title = element_text(size = 12))
grid.arrange(relativeHS2, relativeHS3)

## Figure S6
relativeSeasonHS2 <- relativeSeasonHS2 + 
  ylab(expression(paste("Relative use of habitat ", italic(" i"), " (", italic("d"["i"]), ")"))) +
  theme(axis.title.y = element_text(size = 20, angle = 90)) + 
  labs(title = "(a) Second-order selection") + 
  theme(plot.title = element_text(size = 20))
relativeSeasonHS3 <- relativeSeasonHS3 + 
  ylab(expression(paste("Relative use of habitat ", italic(" i"), " (", italic("d"["i"]), ")"))) +
  theme(axis.title.y = element_text(size = 20, angle = 90)) + 
  labs(title = "(b) Third-order selection") + 
  theme(plot.title = element_text(size = 20))
grid.arrange(relativeSeasonHS2, relativeSeasonHS3, nrow=1)
