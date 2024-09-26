####################################################
#                                                  #
# R file for statistical evaluation                #
#                                                  #
####################################################


# packages and data preparation


# packages
library(xlsx); library(gplots); library(nlme); library(piecewiseSEM)
library(multcomp); library(lsmeans); library(car); library(MASS)
library(arm); library(pbkrtest)


# data import:
DATA <- read.table("Data_Pathogenicity.txt", header = TRUE, sep = "\t")

# factor declarations:
DATA$Pythium_species <- as.factor(DATA$Pythium_species)
str(DATA)


# Germination
# new variables
DATA$GerminationAbs <- as.integer(DATA$Germination*3)
DATA$AnzPfl <- 3
DATA$Non <- as.integer(DATA$AnzPfl - DATA$GerminationAbs)

# model and anova
mod1 <- bayesglm(cbind(DATA$GerminationAbs,DATA$Non) ~ Pythium_species,
  data=DATA, family=binomial(link="logit"))
mod2 <- update(mod1, . ~ 0 + Pythium_species)

# model estimates:
# summary(mod2)
logit <- coef(mod2)
exp(logit) / (1 + exp(logit)) 

# anova:
Anova(mod1, type="II")


# Which levels of Pythium_species differ?
summary(glht(mod1, mcp(Pythium_species="Dunnett")))



# Shoot_weight_fresh_per_plant 
# boxplot:
windows(7,7); par(mfrow=c(1,1), mar=c(8, 4, 4, 2) + 0.1)
boxplot2(Shoot_weight_fresh_per_plant ~ Pythium_species, data=DATA, main="Shoot_weight_fresh_per_plant", xlab="", ylab="", las=3)

# model and anova
mod1 <- gls(Shoot_weight_fresh_per_plant ~ Pythium_species, data=DATA)
mod2 <- update(mod1, . ~ 0 + Pythium_species)

# residual analysis:
windows(); plot(mod1, main="Shoot_weight_fresh_per_plant")

# model estimates:
summary(mod2)

# pseudo R^2:
rsquared(mod1)

# anova:
anova(mod1)

# Which levels of Pythium_species differ?
summary(glht(mod1, mcp(Pythium_species="Dunnett"), df=mod2$dims$N - mod2$dims$p))



# Root_weight_dry_per_plant
# boxplot:
windows(7,7); par(mfrow=c(1,1), mar=c(8, 4, 4, 2) + 0.1)
boxplot2(Root_weight_dry_per_plant ~ Pythium_species, data=DATA, main="Root_weight_dry_per_plant", xlab="", ylab="", las=3)

# model and anova
mod1 <- gls(Root_weight_dry_per_plant ~ Pythium_species, data=DATA)
mod2 <- update(mod1, . ~ 0 + Pythium_species)

# residual analysis:
windows(); plot(mod1, main="Root_weight_dry_per_plant")

# model estimates:
summary(mod2)

# pseudo R^2:
rsquared(mod1)

# anova:
anova(mod1)

# Which levels of Pythium_species differ?
summary(glht(mod1, mcp(Pythium_species="Dunnett"), df=mod2$dims$N - mod2$dims$p))