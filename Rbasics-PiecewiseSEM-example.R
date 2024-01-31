library(tidyverse) # Main graphing package w/ ggplot
library(cowplot) # Additional graphics feature
library(lme4) # mixed models: Use for generalized linear mixed models (glmer)
library(nlme) # mixed models: Use for general linear mixed effects models (lme)
library(emmeans) # use for extracting ls means, pairwise comparisons, and linear contrasts
library(piecewiseSEM) # Use for path analysis and SEM
library(MuMIn) # use calculates marginal & conditional R2 (Hagawatha) Code = r.squaredGLMM(AnovaModelName)
library(ape) # use for partitioning variance components with varcomp commant 
library(effects) # extracts fitted estimates; code = effect("expvar1:expvar2", model) 
library(nlmeU) # power analysis for mixed models
library(ggplot2)
library(predictmeans)
library(ggeffects)
library(multcompView)
library(multcomp)
library(dplyr)
library(partR2)
library(ggsignif)
library(car)


#oo <- options(repos = "https://cran.r-project.org/")
#install.packages("Matrix")
#install.packages("lme4")
#options(oo)

###########################################################
###########################################################
# PIECEWISE SEM WORKED EXAMPLE
###########################################################
###########################################################
# PIECEWISE SEM: Fit models for each response and then piece together inferences 
# using tests of direct separation by fitting missing relationships to test whether
# path coefficients are significantly different from zero. Chi-square test statistic
# compares observed vs. estimated variance-covariance matrices. Significance indicates
# that missing information (i.e., path) would improve model fit.

# STEP 1: Local estimation of models (hierarchcial structure, generalized linear, etc.)
# STEP 2: Structured hypothesis - Identify set of missing relationships using "basis set"
# STEP 3: Test whether effect is not significantly different from zero (p > 0.05)
#         when controlling for covariates already specified in model; i.e, the test
#         considers the partial effect of one variable on response if either or both
#         are already connect to other variables in the model
# STEP 4: Identify best fit model & report path diagram w/ standardized path coefficients,
#         significance-levels, marginal & conditional r2, correlation coefficients
# STEP 5: Use SEM to extract partial residuals (effects) and model separately when appropriate

# REFER TO WEB BOOK ** https://jslefche.github.io/sem_book/index.html **
###########################################################
###########################################################

###########################################################
# DATA WRANGLING
###########################################################
#psa <-  read.csv("//Users/jmw309/Documents/R/PSA-CE2-Weeds-Yield-Rye-Compiled.csv", na.strings = "na")
psa<-read.csv("~/Desktop/PSA-CE2-Weeds-Yield-Rye-Compiled.csv", na.strings='na')
psa$Rye.kgha <- as.integer(psa$Rye.kgha)
psa$CC.win <- as.integer(psa$CC.win) # create continuous variable for "CC termination window; 0 to 21 days"
psa$RR.tot <- psa$Total_density/psa$UTC.tot.density.rep # create response ratio standardized to no cover control at replicate level
psa$LRR.tot <- log((psa$Total_density/psa$UTC.tot.density.rep)+0.01) # log response ratio
psa <- subset(psa, Herb == "POST")
psa <- subset(psa, Crop == "Corn"|Crop == "Beans")
psa <- subset(psa, CC == "1-3 DAP"|CC == "14-21 DPP"|CC == "3-7 DPP")

###########################################################
# STEP 1 LOCAL ESTIMATION: Model 1 cover crop response to planting date
###########################################################
# 1. Inspect univariate statistics
library(tidyverse)
# cheack data by using summary (ie length is right, no missing data
###CC.win is the window of cover crop termiantion timings )
summary <- psa %>% # name new dataset (summary) and call raw data set to summarize (data)
  group_by (CC, CC.win, CropYear, Block) %>% # identify explanatory variables you want to summarize over
  summarise_at(vars(Rye.kgha), funs(mean, length)) # identify predictor variables (can multiple) and type of summary stat
print.data.frame(summary) # output summary table

# 2a. LME Regression (CHECK HOMOGENEITY OF VARIANCE)
######check assumptions of normaility and homogeneity to see if you need to transform data
m1 <- lmer(Rye.kgha      ~ CC.win + (1|CropYear), data = psa, na.action = na.omit) 
plot(m1) 
m2 <- lmer(log(Rye.kgha) ~ CC.win + (1|CropYear), data = psa, na.action = na.omit)
plot(m2)  # inspect residuals to inspect homogeneity of variance (choose best fit)
anova (m1, m2)

# 2b. LME Regression (RANDOM INTERCEPT MODEL VS. RANDOM INTERCEPT/SLOPE MODEL)
m1 <- lmer(log(Rye.kgha) ~ CC.win + (1|CropYear), data = psa, na.action = na.omit) # Random intercept model
m2 <- lmer(log(Rye.kgha) ~ CC.win + (CC.win|CropYear), data = psa, na.action = na.omit) # Random intercept and slope model
anova(m1, m2) # check AICs

# 3. Extract test-statistics and parameter estimates of best fit model
m0 <- lmer(log(Rye.kgha) ~ (CC.win|CropYear), data = psa, na.action = na.omit) 
m1 <- lmer(log(Rye.kgha) ~ CC.win + (CC.win|CropYear), data = psa, na.action = na.omit)
plot(m1)      # inspect residuals to inspect homogeneity of variance (choose best fit)
summary(m1)   # inspect parameter estimates
anova(m0, m1) # model significance: extract Wald Chi Sq & p-value
library(piecewiseSEM) # package to extract marginal & conditional r2 values
#####gives you 2 r2 coefficients. marginal tells you what proportion variance of rye bio is explained by fixed effect.
#######conditional tells you proportion of variance in response that is explained by fixed AND random effect. what's left over the is residual error. 

rsquared(m1, method=NULL) # use to extract conditional r2 value 

# 4a. Extract and plot fitted model at the RANDOM EFFECTS LEVEL
library(merTools)
psa$fit.ranef = predict(m1) 
psa$fit.ranef
fig.rand <- ggplot(psa, aes(x = CC.win, y = log(Rye.kgha), shape = CropYear)) +
  theme_light() +
  geom_point(size=2.0, alpha=0.9) + # point features
  geom_smooth(method="lm", aes(x = CC.win, y = fit.ranef, lty = CropYear), size = 0.5, alpha = 0.2) +
  facet_grid (. ~ CropYear) +
  ylab ("log" ~ (kg ~ ha^{-1})) +
  xlab("Termination window (0-21 d)") +
  ggtitle ("Cereal rye biomass") +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12, face = "bold"))
fig.rand
# extra code  
#####
  #theme(axis.title.y = element_blank()) +
  #theme(panel.grid.minor.x = element_blank()) +
  #theme(panel.grid.major.x = element_blank()) +
  #scale_fill_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
  #scale_color_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
  #scale_x_continuous (limits = c(0, 1)) + 
  #scale_y_continuous (limits = c(-99, 20)) + # setting axis range; scale_y_continuous or scale_x_continuous
#####

# 4b. Extract and plot fitted model at the POPULATION LEVEL
##### model averaged over regression coeffifients
psa$fit.popef = predict(m1, re.form = NA) 
fig.popest <- ggplot(psa, aes(x = CC.win, y = log(Rye.kgha), shape = CropYear)) +
  theme_light() +
  geom_point(size=2.0, alpha=0.9) + # point features
  geom_smooth(method="lm", aes(x = CC.win, y = fit.popef), size = 0.5, alpha = 0.2) +
  ylab ("log" ~ (kg ~ ha^{-1})) +
  xlab("Termination window (0-21 d)") +
  ggtitle ("Cereal rye biomass") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12, face = "bold"))
fig.popest  
  #####
  #theme(axis.title.y = element_blank()) +
  #theme(panel.grid.minor.x = element_blank()) +
  #theme(panel.grid.major.x = element_blank()) +
  #scale_fill_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
  #scale_color_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
  #scale_x_continuous (limits = c(0, 1)) + 
  #scale_y_continuous (limits = c(-99, 20)) + # setting axis range; scale_y_continuous or scale_x_continuous
  #####
  ggsave(file="/Users/johnwallace/Documents/Documents/R/PSA.f1.jpg", scale=1.5, width=4.0, height=2.0, units="in", dpi=900)

# 5a. Extract and plot INTERCEPTS for random effects (i.e., BLUPs or conditional means)
library(broom.mixed) # package for extracting random effects of interest w/ CIs
r1 <- broom.mixed::tidy(m1, effects = "ran_vals", conf.int = TRUE, response = TRUE)
as.data.frame(r1) ### conditional means (interclass correlation coefficient (ICC)). conditioned on pop level estimate (relative to 0)

intercept <- subset(r1, term == "(Intercept)")
fig.int <-ggplot(intercept, aes (x = reorder(level, -estimate), y = estimate, ymin = conf.low, ymax = conf.high, fill = level, color = level, shape = level)) + # extracted from summary data
  geom_hline(yintercept = 0, color = "gray85", size=0.75) +
  geom_errorbar(width = 0.0, size=0.5, position=position_dodge(0.9)) + 
  geom_point(size=3.0, shape=21, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle ("Intercepts") + # optional title in quotes
  theme_light () + # several themes available, other: theme_minimal, theme_grey, theme_classic
  xlab (NULL) + # use this to remove axis title
  ylab ("Conditional mean biomass (95% CI)") + #Label syntax is complicated; but is a typical layout
  theme(legend.position = "NULL") + #legend placement, can be: right, left, top
  theme(legend.title=element_blank()) + #removes legend title
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  #theme(panel.grid.major.x = element_blank()) +
  #theme(panel.grid.major.y = element_blank()) +
  #theme(panel.grid.minor.x = element_blank()) +
  theme(plot.title = element_text(size = 12, face = "bold"))
fig.int
#### 0 is population level. shows deviation from overall average
#ggsave(file="/Users/johnwallace/Documents/Documents/R/PSA.f1.jpg", scale=1.5, width=4.0, height=2.0, units="in", dpi=900)

# 5b. Extract and plot SLOPES for random effects (i.e., BLUPs or conditional means)
r1 <- broom.mixed::tidy(m1, effects = "ran_vals", conf.int = TRUE, response = TRUE)
as.data.frame(r1) 
slope <- subset(r1, term == "CC.win")
fig.int <-ggplot(slope, aes (x = reorder(level, -estimate), y = estimate, ymin = conf.low, ymax = conf.high, fill = level, color = level, shape = level)) + # extracted from summary data
  geom_hline(yintercept = 0, color = "gray85", size=0.75) +
  geom_errorbar(width = 0.0, size=0.5, position=position_dodge(0.9)) + 
  geom_point(size=3.0, shape=21, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle ("Slope") + # optional title in quotes
  theme_light () + # several themes available, other: theme_minimal, theme_grey, theme_classic
  xlab (NULL) + # use this to remove axis title
  ylab ("Conditional mean biomass (95% CI)") + #Label syntax is complicated; but is a typical layout
  #scale_y_continuous (limits = c(-60, 60)) + # setting axis range; scale_y_continuous or scale_x_continuous
  theme(legend.position = "NULL") + #legend placement, can be: right, left, top
  theme(legend.title=element_blank()) + #removes legend title
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  #theme(panel.grid.major.x = element_blank()) +
  #theme(panel.grid.major.y = element_blank()) +
  #theme(panel.grid.minor.x = element_blank()) +
  theme(plot.title = element_text(size = 12, face = "bold"))

fig.int
#ggsave(file="/Users/johnwallace/Documents/Documents/R/PSA.f1.jpg", scale=1.5, width=4.0, height=2.0, units="in", dpi=900)

# 6. ICC (Intraclass Correlation Coefficients); proportion of total variance explained at each grouping level *for nested structures*
icc <-  VarCorr(m1) %>% 
  as_data_frame() %>%
  mutate(icc=vcov/sum(vcov))
icc.df <- as.data.frame(icc)
icc.df
var <- c("Variable", "Variable")
grp <- c("GroupA", "Groupb")
icc <- c(00, 00) # input these ICC scores icc.df data frame using percentages
df <- data.frame(var, grp, icc) # create a data frame
df$grp = factor(df$grp, levels=c("GroupA", "GroupB")) # Order explanatory variables in order of ICC score (descending)
cbPalette <- c("#999999", "#E69F00")
ICC.plot <-ggplot(df, aes (x = var, y = icc)) + # extracted from summary data
  #geom_bar(aes(fill=grp, position="fill", stat="identity")) + # use this when you have more than one stacked bar
  geom_col(aes(fill=grp)) + # use this if plotting a single stacked bar
  coord_flip() + 
  ggtitle ("Intraclass correlation coefficient") + # optional title in quotes
  theme_light() + # several themes available, other: theme_minimal, theme_grey, theme_classic
  xlab (NULL) + # use this to remove axis title
  ylab ("variance explained (%)") + #Label syntax is complicated; but is a typical layout
  theme(legend.position = "bottom") + #legend placement, can be: right, left, top
  theme(legend.title=element_blank()) + #removes legend title
  #scale_y_continuous (limits = c(0, 100)) + 
  # scale_x_break(c(0.75, 0.95)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.text = element_text(size=10)) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values=c("GroupA" = "#CCCCCC", "GroupB" = "#E69F00")) 
ICC.plot
ggsave(file="/Users/johnwallace/Documents/Documents/R/vt-f1b.jpg", scale=1.33, width=3.5, height=1.5, units="in", dpi=900)

###########################################################
# STEP 1 LOCAL ESTIMATION: Model 2 cover crop response to planting date
###########################################################
#### ie doing the same thing as before but for the 2nd component model for PSEM
# 1. Inspect univariate statistics
library(tidyverse)
summary <- psa %>% # name new dataset (summary) and call raw data set to summarize (data)
  group_by (Rye.kgha, CC.win, CropYear, Block) %>% # identify explanatory variables you want to summarize over
  summarise_at(vars(RR.tot, LRR.tot), funs(mean, length)) # identify predictor variables (can multiple) and type of summary stat
print.data.frame(summary) # output summary table

# 2a. LME Regression (CHECK FOR HOMOGENEITY OF VARIANCE)
m1 <- lmer(RR.tot      ~ Rye.kgha + (1|CropYear), data = psa, na.action = na.omit) #response ratio
plot(m1)
m2 <- lmer(LRR.tot     ~ Rye.kgha + (1|CropYear), data = psa, na.action = na.omit) #log response ratio
plot(m2)      # inspect residuals to inspect homogeneity of variance (choose best fit)

# 2b. LME Regression (RANDOM INTERCEPT MODEL VS. RANDOM SLOPE/INTERCEPT MODEL; check AIC)
m1 <- lmer(LRR.tot      ~ Rye.kgha + (1|CropYear), data = psa, na.action = na.omit) 
m2 <- lmer(LRR.tot      ~ Rye.kgha + (Rye.kgha|CropYear), data = psa, na.action = na.omit) 
anova(m1, m2) # check AICs ### random slop did not improve fit here

#### NOTE: when calculation response ratios without spiked weeds, should standardize weeds to overall weedy check (not at replicate level)

# 3. Extract test-statistics and parameter estimates of best fit model
m0 <- lmer(LRR.tot ~            (1|CropYear), data = psa, na.action = na.omit) 
m1 <- lmer(LRR.tot ~ Rye.kgha + (1|CropYear), data = psa, na.action = na.omit)
plot(m1)      # inspect residuals to inspect homogeneity of variance (choose best fit)
summary(m1)   # inspect parameter estimates
anova(m0, m1) # model significance: extract Wald Chi Sq & p-value
library(piecewiseSEM) # package to extract marginal & conditional r2 values
rsquared(m1, method=NULL) # use to extract conditional r2 value 

# 4a. Extract and plot fitted model at the random effects level
library(merTools)
psa$fit.ranef = predict(m1) 
fig.rand <- ggplot(psa, aes(x = Rye.kgha, y = LRR.tot, shape = CropYear)) +
  theme_light() +
  geom_point(size=2.0, alpha=0.9) + # point features
  geom_smooth(method="lm", aes(x = Rye.kgha, y = fit.ranef, lty = CropYear), size = 0.5, alpha = 0.2) +
  facet_grid (. ~ CropYear) +
  xlab ("Cereal rye biomass" ~ (kg ~ ha^{-1})) +
  ylab("LRR") +
  ggtitle ("Weed density relative to control") +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=10)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12, face = "bold"))
fig.rand ### rand intercept model: diff interecepts, same slope
#####
#theme(axis.title.y = element_blank()) +
#theme(panel.grid.minor.x = element_blank()) +
#theme(panel.grid.major.x = element_blank()) +
#scale_fill_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
#scale_color_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
#scale_x_continuous (limits = c(0, 1)) + 
#scale_y_continuous (limits = c(-99, 20)) + # setting axis range; scale_y_continuous or scale_x_continuous
#####

# 4B. Extract and plot fitted model at the population-estimate level
psa$fit.popef = predict(m1, re.form = NA) 
fig.popest <- ggplot(psa, aes(x = Rye.kgha, y = LRR.tot, shape = CropYear)) +
  theme_light() +
  geom_point(size=2.0, alpha=0.9) + # point features
  geom_smooth(method="lm", aes(x = Rye.kgha, y = fit.popef), size = 0.5, alpha = 0.2) +
  xlab ("Cereal rye biomass" ~ (kg ~ ha^{-1})) +
  ylab("LRR") +
  ggtitle ("Weed density relative to control") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(plot.title = element_text(size = 12, face = "bold"))
fig.popest  ### gaining biomass gradient by combining site years. but not explaining lots of variation in the model
#####
#theme(axis.title.y = element_blank()) +
#theme(panel.grid.minor.x = element_blank()) +
#theme(panel.grid.major.x = element_blank()) +
#scale_fill_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
#scale_color_manual(values=c("14-21 DPP" = "peru", "3-7 DPP"="gold", "1-3 DAP"="sea green")) +
#scale_x_continuous (limits = c(0, 1)) + 
#scale_y_continuous (limits = c(-99, 20)) + # setting axis range; scale_y_continuous or scale_x_continuous
#####
ggsave(file="/Users/johnwallace/Documents/Documents/R/PSA.f1.jpg", scale=1.5, width=4.0, height=2.0, units="in", dpi=900)

# 5a. Extract and plot intercepts for random effects (i.e., BLUPs or conditional means)
library(broom.mixed) # package for extracting random effects of interest w/ CIs
r1 <- broom.mixed::tidy(m1, effects = "ran_vals", conf.int = TRUE, response = TRUE)
as.data.frame(r1) 
intercept <- subset(r1, term == "(Intercept)")
fig.int <-ggplot(intercept, aes (x = reorder(level, -estimate), y = estimate, ymin = conf.low, ymax = conf.high, fill = level, color = level, shape = level)) + # extracted from summary data
  geom_hline(yintercept = 0, color = "gray85", size=0.75) +
  geom_errorbar(width = 0.0, size=0.5, position=position_dodge(0.9)) + 
  geom_point(size=3.0, shape=21, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle ("Intercepts") + # optional title in quotes
  theme_light () + # several themes available, other: theme_minimal, theme_grey, theme_classic
  xlab (NULL) + # use this to remove axis title
  ylab ("Conditional mean (95% CI)") + #Label syntax is complicated; but is a typical layout
  theme(legend.position = "NULL") + #legend placement, can be: right, left, top
  theme(legend.title=element_blank()) + #removes legend title
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  #theme(panel.grid.major.x = element_blank()) +
  #theme(panel.grid.major.y = element_blank()) +
  #theme(panel.grid.minor.x = element_blank()) +
  theme(plot.title = element_text(size = 12, face = "bold"))
fig.int
ggsave(file="/Users/johnwallace/Documents/Documents/R/PSA.f1.jpg", scale=1.5, width=4.0, height=2.0, units="in", dpi=900)


# 6. ICC (Intraclass Correlation Coefficients); proportion of total variance explained at each grouping level *for nested structures*
icc <-  VarCorr(m1) %>% 
  as_data_frame() %>%
  mutate(icc=vcov/sum(vcov))
icc.df <- as.data.frame(icc)
icc.df
var <- c("Variable", "Variable")
grp <- c("GroupA", "Groupb")
icc <- c(00, 00) # input these ICC scores icc.df data frame using percentages
df <- data.frame(var, grp, icc) # create a data frame
df$grp = factor(df$grp, levels=c("GroupA", "GroupB")) # Order explanatory variables in order of ICC score (descending)
cbPalette <- c("#999999", "#E69F00")
ICC.plot <-ggplot(df, aes (x = var, y = icc)) + # extracted from summary data
  #geom_bar(aes(fill=grp, position="fill", stat="identity")) + # use this when you have more than one stacked bar
  geom_col(aes(fill=grp)) + # use this if plotting a single stacked bar
  coord_flip() + 
  ggtitle ("Intraclass correlation coefficient") + # optional title in quotes
  theme_light() + # several themes available, other: theme_minimal, theme_grey, theme_classic
  xlab (NULL) + # use this to remove axis title
  ylab ("variance explained (%)") + #Label syntax is complicated; but is a typical layout
  theme(legend.position = "bottom") + #legend placement, can be: right, left, top
  theme(legend.title=element_blank()) + #removes legend title
  #scale_y_continuous (limits = c(0, 100)) + 
  # scale_x_break(c(0.75, 0.95)) + 
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.text = element_text(size=10)) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values=c("GroupA" = "#CCCCCC", "GroupB" = "#E69F00")) 
ICC.plot
ggsave(file="/Users/johnwallace/Documents/Documents/R/vt-f1b.jpg", scale=1.33, width=3.5, height=1.5, units="in", dpi=900)


###########################################################
# STEPS 2-5: Hypothesized diagrams & tests of independence in piecewiseSEM package
###########################################################

#####everything so far was coming up with best fit model for these submodels below:
### need to avoid having oversaturated, need a missing path
basis.psa = psem(
  lmer(Rye.kgha ~ CC.win   + (CC.win|CropYear), data = psa, na.action = na.omit),
  lmer(LRR.tot  ~ Rye.kgha + (1|CropYear),      data = psa, na.action = na.omit)
)
basisSet(basis.psa) # Identifies indepedence claims in basis set
dSep (basis.psa)    # Test of directed separation for individual claim; (should i include another term?)
                    # p < 0.05 indicates it should be included to improve model fit
fisherC(basis.psa)  # Model-wide test-statistic and P-value when multiple dSep tests (did my data fit the SEM?)
                    # p > 0.05 indicates that data support the hypothesized model structure
summary(basis.psa)  # Produces all test statistics, reg & correlation coefficients
                    # No chi-sq estimate or warnings may result from convergence issue in submodel
                    # Either tweak model parameters to converge or rely on Fisher's C
### use Std.Estimate (standardized regression coefficient), not "estimate" to make inference about strngth of the relationship
##### put those on the final SEM fig.
part.resid <- partialResid(LRR.tot ~ CC.win, basis.psa) # Extract and inspect PARTIAL RESIDUALS 
                                                    # Removes variation explained by x1 so partial of effect of x2 can be modeled
with(part.resid, plot(xresid, yresid))              # plot partial residuals
as.data.frame(part.resid)                           # extract as data.frame
### use to test if something else you're not modeling is going on (or if it's not linear)

###########################################################













