####################################################
#                                                  #
# R file for statistical evaluation                #
#                                                  #
####################################################


# packages and data preparation

# packages:
library(xlsx); library(gplots); library(nlme); library(piecewiseSEM)
library(multcomp); library(lsmeans); library(sandwich); library(dplyr); library(broom)

export_stats <- function(stats, filename) {
  write.table(broom::tidy(stats) %>% dplyr::mutate_if(is.numeric, round, 5), file = paste0("../Statistics/", filename, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# data import:
pythium_species <- read.table("../Data/ASV_table_pythium_species.txt", header = TRUE, sep = "\t")
pythium_species <- pythium_species %>% mutate_if(is.numeric, ~1 * (. != 0))
pyth_metadata <- read.table("../../../output/Metadata_pythium.txt", header = TRUE, sep = "\t")

pythium_matrix <- merge(pyth_metadata, pythium_species, by.x = "Id", by.y = "Id")

df <- tidyr::gather(data = pythium_matrix,
                    key = Taxa,
                    value = Abundance,
                    Globisporangium.sylvaticum:Globisporangium.intermedium)

df <- df %>% filter(Abundance != "0")
df <- df %>% group_by(Year, Location, Country, Soilweight) %>% summarise(n_taxa = n_distinct(Taxa))
DATA <- df %>% group_by(Year, Location, Country, Soilweight) %>% summarise(mean_taxa_count = mean(n_taxa))

# factor declarations:
for (c in c("Location", "Year", "Country", "Soilweight")) {
  DATA[[c]] <- as.factor(DATA[[c]])
}
DATA$PF <- as.factor(paste(DATA$Year, DATA$Country, DATA$Soilweight))
str(DATA)

# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(mean_taxa_count ~ Year+Country, data=DATA, main="mean_taxa_count", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))


# number_of_species_pro_soilweight
# model and anova
mod1 <- gls(mean_taxa_count ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year + Soilweight)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="mean_taxa_count")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, spec="Soilweight"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year","Soilweight")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which soilweights differ?
comp3 <- lsmeans(mod3, specs="Soilweight", contr="eff")$contrasts
summary(as.glht(comp3, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp3, by=NULL, alternative="two.sided")), "Soilweight_mean_taxa_count")

# abundance_per_soilweight
# import data
pythium_species <- read.table("../Data/ASV_table_pythium_species.txt", header = TRUE, sep = "\t")
pythium_matrix <- merge(pyth_metadata, pythium_species, by.x = "Id", by.y = "Id")

df <- tidyr::gather(data = pythium_matrix,
                    key = Taxa,
                    value = Abundance,
                    Globisporangium.sylvaticum:Globisporangium.intermedium)

df <- df %>% filter(Abundance != "0")
DATA <- df %>% group_by(Year, Location, Country, Soilweight) %>% summarise(Abundance = sum(Abundance))


# factor declarations:
for (c in c("Location","Year","Country","Soilweight")) {
  DATA[[c]] <- as.factor(DATA[[c]])
}
DATA$PF <- as.factor(paste(DATA$Year, DATA$Country, DATA$Soilweight))
str(DATA)

# boxplot:
#windows(10,7); par(mfrow=c(1,1), mar=c(10, 4, 4, 2) + 0.1)
boxplot2(Abundance ~ Year+Country, data=DATA, main="Abundance", xlab="", ylab="", las=3,
  col=rainbow(length(levels(DATA$Year))))


# model and anova
mod1 <- gls(Abundance ~ PF,
  data=DATA,
  na.action=na.exclude,
  control=list(maxIter=500,msMaxIter=500,niterEM=500,msMaxEval=500,opt="nlminb"))
mod2 <- update(mod1, . ~ 0 + PF)
AIC(mod1,mod2)

mod3 <- update(mod1, . ~ Country + Year + Soilweight)
anova(mod1,mod3)

# residual analysis:
#windows(); plot(mod3, main="Abundance")

# model estimates:
summary(mod3)
as.data.frame(lsmeans(mod3, specs=~1))
as.data.frame(lsmeans(mod3, specs="Country"))
as.data.frame(lsmeans(mod3, specs=~"Year"))
as.data.frame(lsmeans(mod3, spec="Soilweight"))
as.data.frame(lsmeans(mod3, specs=c("Country","Year","Soilweight")))

# pseudo R^2:
rsquared(mod3)

# anova:
anova(mod3)

# Which soilweights differ?
comp3 <- lsmeans(mod3, specs="Soilweight", contr="eff")$contrasts
summary(as.glht(comp3, by=NULL, alternative="two.sided"))
export_stats(summary(as.glht(comp3, by=NULL, alternative="two.sided")), "Soilweight_Abundance")
