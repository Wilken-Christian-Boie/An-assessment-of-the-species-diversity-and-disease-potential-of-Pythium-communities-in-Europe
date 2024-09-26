####################################################
#                                                  #
# R file for statistical evaluation                #
#                                                  #
####################################################


# packages and data preparation

# packages:
library(xlsx); library(gplots); library(nlme); library(piecewiseSEM)
library(multcomp); library(lsmeans); library(sandwich)
library(broom); library(dplyr)

export_stats <- function(stats, filename) {
  write.table(broom::tidy(stats) %>% dplyr::mutate_if(is.numeric, round, 5), file = paste0("../Statistics/", filename, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# data_import
DATA <- read.table("../Data/Pythium_matrix_relativ.txt", header = TRUE, sep = "\t")
str(DATA)

# factor declarations:
for (c in c("Location","Year","Country","Compartment")) {
  DATA[[c]] <- as.factor(DATA[[c]])
}
DATA$PF <- as.factor(paste(DATA$Year,DATA$Country))
str(DATA)

# subset:
DATA <- droplevels(subset(DATA, Compartment=="Soil"))


# Globisporangium attrantheridium
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.attrantheridium ~ Year+Country, data=DATA, main="Globoisporangium attrantheridium", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.attrantheridium ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Globisporangium attrantheridium")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_attrantheridium_country")


# Globisporangium heterothallicum
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.heterothallicum ~ Year+Country, data=DATA, main="Globisporangium heterothallicum", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.heterothallicum ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Globisporangium.heterothallicum")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_heterothallicum_country")


# Globisporangium sylvaticum
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.sylvaticum ~ Year+Country, data=DATA, main="Globisporangium sylvaticum", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.sylvaticum ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="s__sp")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_sylvaticum_country")


# Pythium aff hydnosporum
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Pythium.aff._hydnosporum ~ Year+Country, data=DATA, main="Pythium aff.hydnosporum", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Pythium.aff._hydnosporum ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Pythium aff hydnosporum")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Pythium_aff_hydnosporum_country")


# Pythium monospermum
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Pythium.monospermum ~ Year+Country, data=DATA, main="Pythium monospermum", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Pythium.monospermum ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Pythium monospermum")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Pythium_monospermum_country")


# Pythium arrhenomanes
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Pythium.arrhenomanes ~ Year+Country, data=DATA, main="Pythium arrhenomanes", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Pythium.arrhenomanes ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Pythium arrhenomanes")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Pythium_arrhenomanes_country")


# Globisporangium rostratifingens
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.rostratifingens ~ Year+Country, data=DATA, main="Globisporangium rostratifingens", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.rostratifingens ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Globisporangium rostratifingens")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_rostratifingens_country")


# Globisporangium ultimum var ultimum
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.ultimum_var._ultimum ~ Year+Country, data=DATA, main="Globisporangium ultimum var ultimum", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.ultimum_var._ultimum ~ PF,
  data=DATA, 
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Globisporangium ultimum var ultimum")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_ultimum_var_ultimum_country")


# Globisporangium intermedium
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.intermedium ~ Year+Country, data=DATA, main="Globisporangium intermedium", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.intermedium ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Globisporangium intermedium")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_intermedium_country")


# Globisporangium apiculatum
# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Globisporangium.apiculatum ~ Year+Country, data=DATA, main="Globisporangium apiculatum", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))

# model and anova
mod1 <- gls(Globisporangium.apiculatum ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Globisporangium apiculatum")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which countries differ?
comp2 <- lsmeans(mod3, specs="Country", contr="eff")$contrasts
summary(as.glht(comp2, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp2, by=NULL, alternative="two.sided")), "Globisporangium_apiculatum_country")
