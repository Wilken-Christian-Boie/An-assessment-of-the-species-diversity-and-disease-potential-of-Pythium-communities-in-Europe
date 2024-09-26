####################################################
#                                                  #
# R file for statistical evaluation                #
#                                                  #
####################################################


# packages and data preparation


# packages:
library(xlsx); library(gplots); library(nlme); library(piecewiseSEM)
library(multcomp); library(lsmeans); library(car); library(MASS); library(SimComp)


# data import:
DATA <- read.table("../Data/qRT_PCR_data.txt", header = TRUE, sep = "\t")
str(DATA)

# factor declarations:
for (c in c("Target","Species","Plant")) {
  DATA[[c]] <- as.factor(DATA[[c]])
}
str(DATA)

DATA$PF <- as.factor(paste(DATA$Target,DATA$Species))
DATA$WDH <- as.factor(paste(DATA$Plant,DATA$Species))
str(DATA)

# boxplot:
#windows(14,7); par(mfrow=c(1,1), mar=c(6, 4, 4, 2) + 0.1)
boxplot2(Expression ~ Species+Target, data=DATA, main="Expression", xlab="", ylab="log", las=3,
  col=rainbow(length(levels(DATA$Species))))

#windows(14,7); par(mfrow=c(1,1), mar=c(6, 4, 4, 2) + 0.1)
boxplot2(Expression ~ Species+Target, data=DATA, main="Expression", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Species))), log="y")

# model and anova
mod1 <- gls(log(Expression) ~ Species*Target,
  data=DATA, 
  correlation=corSymm(form=~as.integer(Target)|WDH),
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="optim"))
mod2 <- update(mod1, . ~ 0 + Species:Target)
AIC(mod1,mod2)


# residual analysis:
#windows(); plot(mod1, main="Expression")

# model estimates:
summary(mod2)
exp(summary(mod2)$tTable[,1])

# pseudo R^2:
rsquared(mod1)

# anova:
anova(mod1, type="sequential")

# Which levels of Species differ?
comp1 <- lsmeans(mod1, specs="Species", by=c("Target"), contr="pairwise")$contrasts
summary(as.glht(comp1, by=NULL, alternative="two.sided"))
