# Resource Selection Function ----
## Load Libraries ----
library(tidyverse)
library(PerformanceAnalytics)
library(usdm)
library(mvnormtest)
library(MANOVA.RM)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(car)
library(lmerTest)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(broom.mixed)
library(dotwhisker)
library(ggeffects)
library(scales)

## Required dataset ----
bio <- read.csv("data/SABE_biometry.csv", header = TRUE)
str(bio)

turtlesRSF <- read.csv("data/SABE_habitatselection.csv", header = TRUE, na = "NA")
str(turtlesRSF)

turtlesRSF$turtle.id <- as.factor(as.character(turtlesRSF$turtle.id))
turtlesRSF$frequency <- as.factor(as.character(turtlesRSF$frequency))
turtlesRSF$date <- as.Date(turtlesRSF$date, format = "%Y-%m-%d")
turtlesRSF$time <- as.POSIXct(turtlesRSF$time, format = "%Y-%m-%d %H:%M:%S")
turtlesRSF$session <- as.factor(turtlesRSF$session)
turtlesRSF$month <- factor(turtlesRSF$month, c("January", "February", "March", "April", "May", "June", 
                                               "July", "August", "September", "October", "November", "December"))
turtlesRSF$year <- as.factor(turtlesRSF$year)
turtlesRSF$season <- factor(turtlesRSF$season, c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season", 
                                                 "2021 Wet Season", "2021 Mating Season", "2021 Dry Season"))
turtlesRSF$behaviour <- as.factor(turtlesRSF$behaviour)
turtlesRSF$habitat <- as.factor(turtlesRSF$habitat)

## Check pairwise correlation between explanatory variables
corrCheck <- turtlesRSF %>% 
  dplyr::select(waterDepth, streamWidth, gravel, pebble, cobble, boulder, litter, canopy)
par(mfrow = c(3,3))
apply(corrCheck, 2, hist)
apply(log(corrCheck + 1), 2, shapiro.test)

chart.Correlation(corrCheck,
                  method="spearman",
                  histogram=TRUE,
                  pch=16) ## correlation coefficient

usdm::vif(corrCheck) ## vif
vifstep(corrCheck, th=4)

## MANOVA ----
## Preliminary analysis to compare if habitat composition differ between session, season and sex
select <- turtlesRSF %>% filter(status == 1)
select <- merge(select, bio, by = "turtle.id")

predictor <- c("waterDepth", "streamWidth", "gravel", "pebble", "cobble", "litter")

select %>% dplyr::select(predictor) %>% 
  summarise(mean = apply(., 2, mean), 
            sd = apply(., 2, sd))

selectD <- select %>% filter(session == "Day")
selectN <- select %>% filter(session == "Night")
select %>% group_by(turtle.id, session) %>% 
            summarise_at(vars(predictor), 
                         list(mean))

### Compare if habitat composition differ between day and night
mshapiro.test(t(select[, names(select) %in% predictor]))

selectDielDiff <- MANOVA.wide(cbind(waterDepth, streamWidth, gravel, pebble, cobble, litter) ~ session, 
                              data = select, subject = "turtle.id", resampling = "paramBS", 
                              iter = 10000, seed = 101)
summary(selectDielDiff) ## significant difference of habitat use between daytime and nighttime

### Compare if habitat composition differ seasonally and among reproductive classes
#### During daytime
mshapiro.test(t(selectD[, names(selectD) %in% predictor]))

#### Seasonal difference
selectDSsnDiff <- MANOVA.wide(cbind(waterDepth, streamWidth, gravel, pebble, cobble, litter, canopy) ~ factor(season), 
                              data = filter(selectD, season != "2021 Dry Season"), subject = "turtle.id", resampling = "paramBS", 
                              iter = 10000, seed = 101)
summary(selectDSsnDiff) ## significant divergence in habitat use between seasons in daytime
simCI(selectDSsnDiff, contrast = "pairwise", type = "Tukey")

#### Sexual difference 
selectDSexDiff <- MANOVA.wide(cbind(waterDepth, streamWidth, gravel, pebble, cobble, litter, canopy) ~ sex, 
                              data = selectD, subject = "turtle.id", resampling = "paramBS", 
                              iter = 10000, seed = 101)
summary(selectDSexDiff) 
simCI(selectDSexDiff, contrast = "pairwise", type = "Tukey")

#### During nighttime
mshapiro.test(t(selectN[, names(selectN) %in% predictor]))

#### Seasonal difference
selectNSsnDiff <- MANOVA.wide(cbind(waterDepth, streamWidth, gravel, pebble, cobble, litter) ~ factor(season), 
                              data = selectN, subject = "turtle.id", resampling = "paramBS", 
                              iter = 10000, seed = 101)
summary(selectNSsnDiff) ## significant divergence in habitat use between season in nighttime
simCI(selectNSsnDiff, contrast = "pairwise", type = "Tukey")

#### Sexual difference
selectNSexDiff <- MANOVA.wide(cbind(waterDepth, streamWidth, gravel, pebble, cobble, litter) ~ sex, 
                              data = selectN, subject = "turtle.id", resampling = "paramBS", 
                              iter = 10000, seed = 101)
summary(selectNSexDiff) 
simCI(selectNSexDiff, contrast = "pairwise", type = "Tukey")

## NMDS ----
selectH <- as.matrix(select[,c("waterDepth", "streamWidth", "gravel", "pebble", "cobble", "boulder", "litter")]) # response variables
selectH <- scale(selectH)
selectDH <- as.matrix(selectD[,c("waterDepth", "streamWidth", "gravel", "pebble", "cobble", "boulder", "litter", "canopy")]) 
selectDH <- scale(selectDH)
selectNH <- as.matrix(selectN[,c("waterDepth", "streamWidth", "gravel", "pebble", "cobble", "boulder", "litter")]) 
selectNH <- scale(selectNH)

select.dist <- vegdist(selectH, method = "euclidean")
select.dist

select.div <- adonis2(select.dist ~ as.factor(select$session), data = select, permutation = 9999) 
select.div

## Nmds plot (diurnal)
selectD.NMDS <- metaMDS(selectDH, distance = "euclidean", k = 3, trymax = 100, autotransform = TRUE)
selectD.NMDS
stressplot(selectD.NMDS)

NMDS1 <- selectD.NMDS$points[,1]
NMDS2 <- selectD.NMDS$points[,2]
selectD.NMDSplot <- cbind(selectD, NMDS1, NMDS2)
selectD.NMDSplot <- selectD.NMDSplot[,!duplicated(colnames(selectD.NMDSplot))]

fit <- envfit(selectD.NMDS, selectDH)
arrows <- data.frame(fit$vectors$arrows, R = fit$vectors$r, P = fit$vectors$pvals)
arrows$FG <- rownames(arrows)
arrows.p <- arrows[arrows$P < 0.05,]
ggplot(data = selectD.NMDSplot, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = selectD.NMDSplot, size = 3, alpha = 0.2, col = "steelblue") + 
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = NMDS1*2, yend = NMDS2*2), 
               size = 1, alpha = 0.5, colour = "black", 
               arrow = arrow(length = unit(.2, "cm")*arrows$R)) +
  geom_text_repel(data = arrows, aes(x = NMDS1*2, y = NMDS2*2, label = FG), 
                  colour = "grey30", fontface = "bold") + 
  theme_classic() + 
  annotate("text", x = 4, y = 4, label = paste('Stress =', round(selectD.NMDS$stress,3)))

## Nmds plot (nocturnal) 
selectN.NMDS <- metaMDS(selectNH, distance = "euclidean", k = 3, trymax = 100, autotransform = TRUE)
selectN.NMDS
stressplot(selectN.NMDS)

NMDS1 <- selectN.NMDS$points[,1]
NMDS2 <- selectN.NMDS$points[,2]
selectN.NMDSplot <- cbind(selectN, NMDS1, NMDS2)
selectN.NMDSplot <- selectN.NMDSplot[,!duplicated(colnames(selectN.NMDSplot))]

fit <- envfit(selectN.NMDS, selectNH)
arrows <- data.frame(fit$vectors$arrows, R = fit$vectors$r, P = fit$vectors$pvals)
arrows$FG <- rownames(arrows)
arrows.p <- arrows[arrows$P < 0.05,]
ggplot(data = selectN.NMDSplot, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = selectN.NMDSplot, size = 3, alpha = 0.2, col = "steelblue") + 
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = NMDS1*2, yend = NMDS2*2), 
               size = 1, alpha = 0.5, colour = "black", 
               arrow = arrow(length = unit(.2, "cm")*arrows$R)) +
  geom_text_repel(data = arrows, aes(x = NMDS1*2, y = NMDS2*2, label = FG), 
                  colour = "grey30", fontface = "bold") + 
  theme_classic() + 
  annotate("text", x = 4, y = 4, label = paste('Stress =', round(selectN.NMDS$stress,3)))



## Resource Selection Functions ----
turtlesRSFD <- turtlesRSF %>% filter(session == "Day")
turtlesRSFN <- turtlesRSF %>% filter(session == "Night")

### Diurnal ----
turtlesRSFD.Mod <- glmmTMB(status ~(scale(waterDepth) + scale(I(waterDepth^2)) + 
                                      scale(streamWidth) + scale(I(streamWidth^2)) + 
                                      scale(gravel) + scale(I(gravel^2)) + 
                                      scale(pebble) + scale(I(pebble^2)) + 
                                      scale(cobble) + scale(I(cobble^2)) + 
                                      scale(litter) + scale(I(litter^2)) + 
                                      scale(canopy) + scale(I(canopy^2))) : (season) + 
                             (0+scale(waterDepth)|turtle.id) + (0+scale(I(waterDepth^2))|turtle.id) + 
                             (0+scale(streamWidth)|turtle.id) + (0+scale(I(streamWidth^2))|turtle.id) + 
                             (0+scale(gravel)|turtle.id) + (0+scale(I(gravel^2))|turtle.id) + 
                             (0+scale(pebble)|turtle.id) + (0+scale(I(pebble^2))|turtle.id) + 
                             (0+scale(cobble)|turtle.id) + (0+scale(I(cobble^2))|turtle.id) +
                             (0+scale(litter)|turtle.id) + (0+scale(I(litter^2))|turtle.id) + 
                             (0+scale(canopy)|turtle.id) + (0+scale(I(canopy^2))|turtle.id), 
                           family = binomial(), data = filter(turtlesRSFD, season != "2021 Dry Season"))
summary(turtlesRSFD.Mod)
Anova(turtlesRSFD.Mod, type = 3)

View(coef(turtlesRSFD.Mod)[["cond"]][["turtle.id"]]) ## species-specific effect on slope

### Check Model Assumptions
turtlesRSFD.Mod.sim <- simulateResiduals(turtlesRSFD.Mod)
plot(turtlesRSFD.Mod.sim)

### Nocturnal ----
turtlesRSFN.Mod <- glmmTMB(status ~(scale(waterDepth) + scale(I(waterDepth^2)) + 
                                        scale(streamWidth) + scale(I(streamWidth^2)) + 
                                        scale(gravel) + scale(I(gravel^2)) + 
                                        scale(pebble) + scale(I(pebble^2)) + 
                                        scale(cobble) + scale(I(cobble^2)) + 
                                        scale(litter) + scale(I(litter^2))) : (season) + 
                                     (0+scale(waterDepth)|turtle.id) + (0+scale(I(waterDepth^2))|turtle.id) + 
                                     (0+scale(streamWidth)|turtle.id) + (0+scale(I(streamWidth^2))|turtle.id) + 
                                     (0+scale(gravel)|turtle.id) + (0+scale(I(gravel^2))|turtle.id) + 
                                     (0+scale(pebble)|turtle.id) + (0+scale(I(pebble^2))|turtle.id) + 
                                     (0+scale(cobble)|turtle.id) + (0+scale(I(cobble^2))|turtle.id) +
                                     (0+scale(litter)|turtle.id) + (0+scale(I(litter^2))|turtle.id),
                                   family = binomial(), data = turtlesRSFN)
summary(turtlesRSFN.Mod)
Anova(turtlesRSFN.Mod, type = 3)

View(coef(turtlesRSFN.Mod)[["cond"]][["turtle.id"]])

### Check Model Assumptions
turtlesRSFN.Mod.sim <- simulateResiduals(turtlesRSFN.Mod)
plot(turtlesRSFN.Mod.sim)

### Plotting ----
#### Plots of Coefficient Estimates ----
#### Figure 5
turtlesRSFD.Mod.coef <- tidy(turtlesRSFD.Mod, conf.int = TRUE)
turtlesRSFD.Mod.coef <- turtlesRSFD.Mod.coef %>% filter(effect == "fixed")
turtlesRSFD.Mod.coef$session <- "Day"

turtlesRSFN.Mod.coef <- tidy(turtlesRSFN.Mod, conf.int = TRUE)
turtlesRSFN.Mod.coef <- turtlesRSFN.Mod.coef %>% filter(effect == "fixed")
turtlesRSFN.Mod.coef$session <- "Night"

turtlesRSF.Mod.coef <- rbind(turtlesRSFD.Mod.coef, turtlesRSFN.Mod.coef)
turtlesRSF.Mod.coef$group <- sub(".*season", "", turtlesRSF.Mod.coef$term)
turtlesRSF.Mod.coef$group <- sub(":.*", "", turtlesRSF.Mod.coef$group)
turtlesRSF.Mod.coef$term <- sub(".*scale", "", turtlesRSF.Mod.coef$term)
turtlesRSF.Mod.coef$term <- sub(":.*", "", turtlesRSF.Mod.coef$term)
turtlesRSF.Mod.coef$term <- sub(paste0(".*", "[(]", "I", "[(]"), "", turtlesRSF.Mod.coef$term)
turtlesRSF.Mod.coef$term <- sub(paste0("[)]", "[)]"), "", turtlesRSF.Mod.coef$term)
turtlesRSF.Mod.coef$term <- gsub("\\(|)", "", turtlesRSF.Mod.coef$term)

turtlesRSF.Mod.coef$fun <- ifelse(grepl('^2', turtlesRSF.Mod.coef$term, fixed=TRUE), "quadratic", "linear")
turtlesRSF.Mod.coef$variable <- gsub("\\^2", "", turtlesRSF.Mod.coef$term)
turtlesRSF.Mod.coef$variable <- gsub("([a-z])([A-Z])","\\1 \\2", turtlesRSF.Mod.coef$variable)
turtlesRSF.Mod.coef$variable <- str_to_title(turtlesRSF.Mod.coef$variable)
turtlesRSF.Mod.coef$variable[turtlesRSF.Mod.coef$variable == "Litter"] <- "Leaf Litter"
turtlesRSF.Mod.coef$variable[turtlesRSF.Mod.coef$variable == "Canopy"] <- "Canopy Cover"

turtlesRSF.Mod.coef <- turtlesRSF.Mod.coef %>% rename(model = group)
turtlesRSF.Mod.coef[!turtlesRSF.Mod.coef$model %in% c("2020 Mating Season", "2020 Dry Season", "2021 Nesting Season", 
                                                      "2021 Wet Season", "2021 Mating Season", "(Intercept)"), "model"] <- "2020 Wet Season"

seasonLv <- c("2020 Wet Season", "2021 Wet Season", "2020 Mating Season", "2021 Mating Season", 
              "2020 Dry Season", "2021 Dry Season", "2021 Nesting Season", "2022 Nesting Season")
emptyRSF <- data.frame(expand.grid(model = factor(seasonLv, levels = seasonLv), 
                                   variable = levels(factor(turtlesRSF.Mod.coef$variable)), 
                                   session = levels(factor(turtlesRSF.Mod.coef$session))))
turtlesRSF.Mod.coef <- merge(emptyRSF, turtlesRSF.Mod.coef, 
                             by = c("model", "variable", "session"), 
                             all.x = TRUE)
turtlesRSF.Mod.coef$effect <- "fixed"
turtlesRSF.Mod.coef$component <- "cond"
turtlesRSF.Mod.coef <- turtlesRSF.Mod.coef %>% filter(variable != "Intercept")
turtlesRSF.Mod.coef$variable <- factor(turtlesRSF.Mod.coef$variable, 
                                       levels = c("Water Depth", "Stream Width",  
                                                  "Gravel", "Pebble", "Cobble",  
                                                  "Leaf Litter", "Canopy Cover")) 
turtlesRSF.Mod.coef$year <- as.numeric(substr(turtlesRSF.Mod.coef$model, 1, 4))
turtlesRSF.Mod.coef$year <- as.factor(turtlesRSF.Mod.coef$year)
turtlesRSF.Mod.coef$season <- sub(".*^\\d+\\s", "", turtlesRSF.Mod.coef$model)
turtlesRSF.Mod.coef$season <- as.factor(turtlesRSF.Mod.coef$season)

turtlesRSF.Mod.coef %>% 
  arrange(fun) %>% 
  filter(model %in% c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season")) %>%
  mutate(sign = ifelse(p.value<0.05, ifelse(estimate>0, "pos", "neg"), "ns")) %>% 
  drop_na(estimate) %>% 
  transform(estimate = abs(estimate), 
            conf.low = abs(estimate) - 1.96*std.error, 
            conf.high = abs(estimate) + 1.96*std.error) %>% 
  ggplot(aes(x = estimate, y = factor(variable, 
                                      levels = rev(levels(factor(variable)))), 
             shape = factor(fun, 
                            levels = rev(levels(factor(fun)))), 
             col = factor(sign, 
                          levels = rev(c("pos", "neg", "ns"))))) + 
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), 
                size = 1, width = 0, 
                position = position_dodge(.8)) + 
  geom_point(size = 4, 
             position = position_dodge(.8)) + 
  theme_bw() +
  xlab("| Coefficient Estimate |") + 
  ylab("") + 
  # scale_y_discrete(labels = rev(y.labels)) +
  scale_shape_manual(values = c("linear" = 16, 
                                "quadratic" = 17), 
                     guide = guide_legend(reverse = TRUE), 
                     labels = rev(c("Linear response", "Quadratic response"))) + 
  scale_colour_manual(values = c("pos" = "#EE6677", 
                                 "neg" = "#4477AA", 
                                 "ns" = "grey"), 
                      guide = guide_legend(reverse = TRUE), 
                      labels = rev(c("Postively significant", "Negatively significant", "Non-significance"))) +
  
  theme(legend.position = "right", 
        legend.title = element_blank()) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  facet_grid(session ~ factor(season, 
                              levels = c("Wet Season", "Mating Season", "Dry Season", "Nesting Season")), 
             scales = "free_y") + 
  theme(text = element_text(size = 14), 
        axis.text = element_text(size = 14), 
        strip.text = element_text(size = 14), 
        legend.text = element_text(size = 14))



#### Plots of Marginal Effects ----
#### Right Panel: Predicted Effect 
turtlesRSFD.waterDepth <- ggpredict(turtlesRSFD.Mod, term = c("waterDepth [0:100, by = 1]", "season"))
turtlesRSFD.waterDepth$variable <- "waterDepth"
turtlesRSFD.streamWidth <- ggpredict(turtlesRSFD.Mod, term = c("streamWidth [all]", "season"))
turtlesRSFD.streamWidth$variable <- "streamWidth"
turtlesRSFD.gravel <- ggpredict(turtlesRSFD.Mod, term = c("gravel [0:100, by = 1]", "season"))
turtlesRSFD.gravel$variable <- "gravel"
turtlesRSFD.pebble <- ggpredict(turtlesRSFD.Mod, term = c("pebble [0:100, by = 1]", "season"))
turtlesRSFD.pebble$variable <- "pebble"
turtlesRSFD.cobble <- ggpredict(turtlesRSFD.Mod, term = c("cobble [0:100, by = 1]", "season"))
turtlesRSFD.cobble$variable <- "cobble"
turtlesRSFD.litter <- ggpredict(turtlesRSFD.Mod, term = c("litter [0:100, by = 1]", "season"))
turtlesRSFD.litter$variable <- "litter"
turtlesRSFD.canopy <- ggpredict(turtlesRSFD.Mod, term = c("canopy [0:100, by = 1]", "season"))
turtlesRSFD.canopy$variable <- "canopy"
turtlesRSFD.ME <- rbind(turtlesRSFD.waterDepth, turtlesRSFD.streamWidth, 
                        turtlesRSFD.gravel, turtlesRSFD.pebble, turtlesRSFD.cobble, 
                        turtlesRSFD.litter, turtlesRSFD.canopy)
turtlesRSFD.ME$session <- "Day"

turtlesRSFN.waterDepth <- ggpredict(turtlesRSFN.Mod, term = c("waterDepth [0:100, by = 1]", "season"))
turtlesRSFN.waterDepth$variable <- "waterDepth"
turtlesRSFN.streamWidth <- ggpredict(turtlesRSFN.Mod, term = c("streamWidth [all]", "season"))
turtlesRSFN.streamWidth$variable <- "streamWidth"
turtlesRSFN.gravel <- ggpredict(turtlesRSFN.Mod, term = c("gravel [0:100, by = 1]", "season"))
turtlesRSFN.gravel$variable <- "gravel"
turtlesRSFN.pebble <- ggpredict(turtlesRSFN.Mod, term = c("pebble [0:100, by = 1]", "season"))
turtlesRSFN.pebble$variable <- "pebble"
turtlesRSFN.cobble <- ggpredict(turtlesRSFN.Mod, term = c("cobble [0:100, by = 1]", "season"))
turtlesRSFN.cobble$variable <- "cobble"
turtlesRSFN.litter <- ggpredict(turtlesRSFN.Mod, term = c("litter [0:100, by = 1]", "season"))
turtlesRSFN.litter$variable <- "litter"
turtlesRSFN.ME <- rbind(turtlesRSFN.waterDepth, turtlesRSFN.streamWidth, 
                        turtlesRSFN.gravel, turtlesRSFN.pebble, turtlesRSFN.cobble, 
                        turtlesRSFN.litter)
turtlesRSFN.ME$session <- "Night"

turtlesRSF.ME <- rbind(turtlesRSFD.ME, turtlesRSFN.ME)
turtlesRSF.ME$year <- substr(turtlesRSF.ME$group, 1, 4)
turtlesRSF.ME$year <- as.factor(turtlesRSF.ME$year)
turtlesRSF.ME$season <- sub(".*^\\d+\\s", "", turtlesRSF.ME$group)
turtlesRSF.ME$season <- factor(turtlesRSF.ME$season, 
                               levels = c("Wet Season", "Mating Season", 
                                          "Dry Season", "Nesting Season"))
turtlesRSF.ME$variable <- as.factor(turtlesRSF.ME$variable)
turtlesRSF.ME$session <- as.factor(turtlesRSF.ME$session)

turtlesRSF.ME.new <- data.frame(matrix(ncol = length(names(turtlesRSF.ME)), nrow = 0))
colnames(turtlesRSF.ME.new) <- names(turtlesRSF.ME)
for(i in 1:length(levels(turtlesRSF.ME$variable))){
  variableName = levels(turtlesRSF.ME$variable)[i]
  for(j in 1:length(levels(turtlesRSF.ME$group))){
    groupName = levels(turtlesRSF.ME$group)[j]
    min = turtlesRSF %>% filter(season == groupName) %>% dplyr::select(variableName) %>% min(., na.rm = TRUE)
    max = turtlesRSF %>% filter(season == groupName) %>% dplyr::select(variableName) %>% max(., na.rm = TRUE)
    selRange <- turtlesRSF.ME %>% filter(group == groupName) %>% filter(variable == variableName) %>% filter(x >= min & x <= max)
    turtlesRSF.ME.new <- rbind(turtlesRSF.ME.new, selRange)
    j = j + 1
  }
  i = i + 1
}
turtlesRSF.ME <- turtlesRSF.ME.new
remove(turtlesRSF.ME.new)

turtlesRSF.ME$variable <- fct_rev(turtlesRSF.ME$variable)
turtlesRSF.ME$variable <- factor(turtlesRSF.ME$variable, levels = c("waterDepth", "streamWidth", 
                                                                    "gravel", "pebble", "cobble", 
                                                                    "litter", "canopy"))

turtlesRSF.MEPlot = list()
for(i in 1:I(length(levels(turtlesRSF.ME$variable)))){
  turtlesRSF.MEPlot[[i]] <- turtlesRSF.ME %>% 
    filter(group %in% c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season")) %>%
    filter(variable != "canopy" | session != "Night") %>%
    filter(variable %in% c(levels(turtlesRSF.ME$variable)[[i]])) %>%
    plyr::mutate(ci_range = conf.high - conf.low) %>%
    filter(ci_range < 0.9) %>%
    ggplot(., aes(x = x, y = predicted, group = group)) +
    geom_line(aes(col = group), size = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) + 
    theme_classic() + 
    ylab("Predicted probabilities") + 
    xlab(str_to_title(gsub("([A-Z])", " \\1", levels(turtlesRSF.ME$variable)[[i]]))) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
    scale_colour_manual(values = c("2020 Wet Season" = colorspace::qualitative_hcl(4, "Dark3")[3],
                                   "2020 Mating Season" = colorspace::qualitative_hcl(4, "Dark3")[4],
                                   "2020 Dry Season" = colorspace::qualitative_hcl(4, "Dark3")[2],
                                   "2021 Nesting Season" = colorspace::qualitative_hcl(4, "Dark3")[1])) +
    scale_fill_manual(values = c("2020 Wet Season" = colorspace::qualitative_hcl(4, "Dark3")[3],
                                 "2020 Mating Season" = colorspace::qualitative_hcl(4, "Dark3")[4],
                                 "2020 Dry Season" = colorspace::qualitative_hcl(4, "Dark3")[2],
                                 "2021 Nesting Season" = colorspace::qualitative_hcl(4, "Dark3")[1])) +
    theme(strip.background.y  = element_blank(),
          strip.text.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom") +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) + 
    theme(text = element_text(size = 16),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16)) +
    facet_grid(factor(season, levels = c("Wet Season", "Mating Season", 
                                         "Dry Season", "Nesting Season"))~session);
  
  if(levels(turtlesRSF.ME$variable)[[i]] != "canopy"){
    turtlesRSF.MEPlot[[i]] <- turtlesRSF.MEPlot[[i]] +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                data.frame(session = "Night"), fill = "darkblue", alpha = .05, inherit.aes = FALSE)
  };
  
  i = i + 1;
}

#### Left Panel: Coefficient Values
turtlesRSF.Mod.coef$term <- factor(turtlesRSF.Mod.coef$term, 
                                   levels = c("waterDepth", "waterDepth^2", "streamWidth", "streamWidth^2", 
                                              "gravel", "gravel^2", "pebble", "pebble^2", "cobble", "cobble^2", 
                                              "litter", "litter^2", "canopy", "canopy^2")) 

turtlesRSF.coefPlot = list()
a = 1
for(i in 1:I(length(levels(turtlesRSF.Mod.coef$term))/2)){
  ## Find x-axis range
  max = turtlesRSF.Mod.coef %>% 
    filter(model %in% c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season")) %>%
    filter(term %in% c(levels(turtlesRSF.Mod.coef$term)[[a]], levels(turtlesRSF.Mod.coef$term)[[a+1]])) %>%
    summarise(max = max(abs(c(estimate, conf.low, conf.high)), na.rm = TRUE));
  
  ## Plotting
  turtlesRSF.coefPlot[[i]] <- turtlesRSF.Mod.coef %>% 
    filter(model %in% c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season")) %>%
    filter(term != "canopy" | session != "Night") %>% 
    filter(term != "canopy^2" | session != "Night") %>%
    filter(term %in% c(levels(turtlesRSF.Mod.coef$term)[[a]], levels(turtlesRSF.Mod.coef$term)[[a+1]])) %>%
    arrange(factor(season, levels = c("Wet Season", "Mating Season", 
                                      "Dry Season", "Nesting Season"))) %>%
    dwplot(whisker_args = list(aes(col = season), size = .6), 
           dot_args = list(aes(col = season, pch = fun), size = 4, fill = "white"), 
           dodge_size = 1, 
           show_intercept = TRUE) + 
    theme_classic() +
    xlab("Coefficient Estimate") + 
    ylab("") + 
    scale_x_continuous(limits = c(-ceiling((as.numeric(max))), ceiling(as.numeric(max)))) + 
    ylab(str_to_title(gsub("([A-Z])", " \\1", levels(turtlesRSF.Mod.coef$term)[[a]]))) + 
    scale_colour_manual(values = c("Wet Season" = colorspace::qualitative_hcl(4, "Dark3")[3], 
                                   "Mating Season" = colorspace::qualitative_hcl(4, "Dark3")[4], 
                                   "Dry Season" = colorspace::qualitative_hcl(4, "Dark3")[2], 
                                   "Nesting Season" = colorspace::qualitative_hcl(4, "Dark3")[1])) +
    scale_shape_discrete(labels = c("Linear Function", "Quadratic Function")) +
    theme(legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          strip.text.y = element_blank()) +
    theme(text = element_text(size = 16),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16)) +
    guides(colour = "none", 
           shape = guide_legend(nrow = 2, byrow = FALSE)) + 
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    facet_grid(factor(season, levels = c("Wet Season", "Mating Season", 
                                         "Dry Season", "Nesting Season"))~session);
  i = i + 1;
  a = a + 2
  
}

## Water Depth
ggarrange(turtlesRSF.coefPlot[[1]], 
          turtlesRSF.MEPlot[[1]], 
          ncol = 2, widths = c(2,3), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
## Stream Width
ggarrange(turtlesRSF.coefPlot[[2]], 
          turtlesRSF.MEPlot[[2]], 
          ncol = 2, widths = c(2,3), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
## Gravel
ggarrange(turtlesRSF.coefPlot[[3]], 
          turtlesRSF.MEPlot[[3]], 
          ncol = 2, widths = c(2,3), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
## Pebble
ggarrange(turtlesRSF.coefPlot[[4]], 
          turtlesRSF.MEPlot[[4]], 
          ncol = 2, widths = c(2,3), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
## Cobble
ggarrange(turtlesRSF.coefPlot[[5]], 
          turtlesRSF.MEPlot[[5]], 
          ncol = 2, widths = c(2,3), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
## Leaf Litter
ggarrange(turtlesRSF.coefPlot[[6]], 
          turtlesRSF.MEPlot[[6]], 
          ncol = 2, widths = c(4,7), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
## Canopy Cover
ggarrange(turtlesRSF.coefPlot[[7]], 
          turtlesRSF.MEPlot[[7]], 
          ncol = 2, widths = c(4,7), 
          labels = c("A", "B"), 
          font.label = list(size = 16))
