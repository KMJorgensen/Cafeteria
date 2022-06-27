# R-script to reproduce analyses for:
#"A critical test of ectomycorrhizal exploration types as predictors of mycelial foraging"
# Jörgensen, K., Clemmensen, KE., Wallander, H., Lindahl, BD. 

# Contact about R script: karolina.jorgensen@slu.se                          

# This has been tested under:
# R version 4.0.1 (2022-06-27)


#ABOUT THE SCRIPT ####
# This script needs the .csv file generated in the script called
# "Cafeteria_Jorgensen_eta_logratios" (also available on GitHub)

# Rationale: Statistical tests ####
# This script describes the statistical approach used in the paper:
# 1. Evaluate whether ectomycorrhizal exploration type is a good predictor for 
#   growth patterns. Specifically, test the effect of exploration type on log-ratio 
#   between roots/bags and soilbags/sandbags. 
# 2. Since exploration type was not a good predictor, evaluate the 
#   growth patterns between roots/bags and soilbags/sandbags on genus level. 
# 3. Test whether soil substrate types (F, RH, P, M soil type) affects community composition.
# 4. Test which genera colonise  substrates differently.
# 5. Test which substrates specific genera prefer. 
# 6. Same approach as (3.) but compare differences between pure sand bags and bags amended with apatite.


# Required packages
# Statistical analyses
library(lme4) # ver 1.1.26
library(lmerTest) # ver 3.1.3
library(stats) # ver 4.0.1
library(emmeans) # ver 1.5.3
library(vegan) # ver 2.5-7

# Plotting and data handling
library(dplyr) # ver 1.0.3
library(tidyr) # ver 1.1.2
library(ggplot2) # ver 3.3.3
library(ggthemes) # ver 4.2.4
library(ggpubr) # ver 0.4.0
library(stringr) # ver 1.4.0

# Load data sets ####
setwd("~/docs/Dokument/Projekt/VR-projekt/Cafeteriaexperiment/Submission_prep/New Phytologist - revision/For GitHub")

data <- read.csv2("genera_logratios.csv", header = TRUE) #from Cafeteria_Jorgensen_etal_logrations.R

# Substrate data
abund.data <- read.csv2("Rel_abundances.csv", header = T)
abund.data$Site.no <- as.factor(abund.data$Site.no) # Site ID
abund.data$Cafeteria.ID <- as.factor(abund.data$Cafeteria.ID) # Cafeteria ID
abund.data$Caf <- as.factor(abund.data$Caf) # Cafeteria replicate, within site
abund.data$Substrate <- as.factor(abund.data$Substrate) # Either root or one of six bag substrates
abund.data$Bagtype <- as.factor(abund.data$Bagtype) # Either root, sandbag (S, A) or soilbag (F, M, RH, P)
abund.data$Bag_root <- as.factor(abund.data$Bag_root) # Root or ingrowth bag

# Add information about exploration types 
data <- data %>%
  mutate(Expl.type = case_when(
    Genus == "Amanita" ~ "MDS",
    Genus == "Amphinema" ~ "MDF",
    Genus == "Cenococcum" ~ "SD",
    Genus == "Cortinarius" ~ "MDF",
    Genus == "Hyaloscypha" ~ "C",
    Genus == "Hygrophorus" ~ "SD",
    Genus == "Lactarius" ~ "C",
    Genus == "Piloderma" ~ "MDF",
    Genus == "Pseudotomentella" ~ "MDS",
    Genus == "Russula" ~ "C",
    Genus == "Tomentella_Thelephora" ~ "MDS",
    Genus == "Tylospora" ~ "SD"
  ))

data$Site <- as.character(data$Site)
data$Cafeteria.ID <- as.character(data$Cafeteria.ID)


# Make separate DFs for analyses
# Root vs. bags (Comparison = Root_bag)
# Soilbags vs. sandbags (Comparison = bagtype)

root_bag <- data %>% filter(Comparison == "Root_bag")
bagtype <- data %>% filter(Comparison == "Bagtype")

# MODELS ####
# Test 1. Exploration type ####
# Evaluate whether ectomycorrhizal exploration type is a good predictor for 
#   growth patterns. Specifically, test the effect of exploration type on log-ratio 
#   between roots/bags and soilbags/sandbags

bagroot.expltype.mod <- lmerTest::lmer(Log_ratio ~ Expl.type + (1|Genus) + (1|Site/Cafeteria.ID), data = root_bag)
plot(bagroot.expltype.mod)
hist(resid(bagroot.expltype.mod))
anova(bagroot.expltype.mod) #p=0.72

bagtype.expltype.mod <- lmerTest::lmer(Log_ratio ~ Expl.type + (1|Genus) + (1|Site/Cafeteria.ID), data = bagtype)
plot(bagtype.expltype.mod)
hist(resid(bagtype.expltype.mod))
anova(bagtype.expltype.mod) #p=0.32

# Test 2. Single genus models ####
# Since exploration type was not a good predictor of colonisation patterns, 
# The preference off each genus will be tested separately
# The Benjamini-Hochberg method will be used to correct for multiple testing

# Make vectors to store p.values for each analysis - this is used to make the p-value corrections
p.value.rootbag <- numeric() # Comparison: root vs. bag preference
p.value.bagtype <- numeric() # Comparison: soil- vs. sandbag preference

# _Amanita ####
amanita.rootbag <- root_bag %>% filter(Genus == "Amanita")
amanita.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = amanita.rootbag)
summary(amanita.mod1)
hist(resid(amanita.mod1))
p.value.rootbag[1] <- coef(summary(amanita.mod1))[1,5]

amanita.bagtype <- bagtype %>% filter(Genus == "Amanita")
amanita.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = amanita.bagtype)
summary(amanita.mod2)
hist(resid(amanita.mod2))
p.value.bagtype[1] <- coef(summary(amanita.mod2))[1,5]

# _Amphinema ####
amphinema.rootbag <- root_bag %>% filter(Genus == "Amphinema")
amphinema.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = amphinema.rootbag)
summary(amphinema.mod1)
hist(resid(amphinema.mod1))
p.value.rootbag[2] <- coef(summary(amphinema.mod1))[1,5]

amphinema.bagtype <- bagtype %>% filter(Genus == "Amphinema")
amphinema.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = amphinema.bagtype)
summary(amphinema.mod2)
hist(resid(amphinema.mod2))
p.value.bagtype[2] <- coef(summary(amphinema.mod2))[1,5]

# _Cenococcum  ####
cenococcum.rootbag <- root_bag %>% filter(Genus == "Cenococcum")
cenococcum.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = cenococcum.rootbag)
summary(cenococcum.mod1)
hist(resid(cenococcum.mod1))
p.value.rootbag[3] <- coef(summary(cenococcum.mod1))[1,5]

cenococcum.bagtype <- bagtype %>% filter(Genus == "Cenococcum")
cenococcum.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = cenococcum.bagtype)
summary(cenococcum.mod2)
hist(resid(cenococcum.mod2))
p.value.bagtype[3] <- coef(summary(cenococcum.mod2))[1,5]

# _Cortinarius ####
cortinarius.rootbag <- root_bag %>% filter(Genus == "Cortinarius")
cortinarius.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = cortinarius.rootbag)
summary(cortinarius.mod1)
hist(resid(cortinarius.mod1))
p.value.rootbag[4] <- coef(summary(cortinarius.mod1))[1,5]

cortinarius.bagtype <- bagtype %>% filter(Genus == "Cortinarius")
cortinarius.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = cortinarius.bagtype)
summary(cortinarius.mod2)
hist(resid(cortinarius.mod2))
p.value.bagtype[4] <- coef(summary(cortinarius.mod2))[1,5]


# _Hyaloscypha ####
hyaloscypha.rootbag <- root_bag %>% filter(Genus == "Hyaloscypha")
hyaloscypha.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = hyaloscypha.rootbag)
summary(hyaloscypha.mod1)
hist(resid(hyaloscypha.mod1))
p.value.rootbag[5] <- coef(summary(hyaloscypha.mod1))[1,5]

hyaloscypha.bagtype <- bagtype %>% filter(Genus == "Hyaloscypha")
hyaloscypha.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = hyaloscypha.bagtype)
summary(hyaloscypha.mod2)
hist(resid(hyaloscypha.mod2))
p.value.bagtype[5] <- coef(summary(hyaloscypha.mod2))[1,5]

# _Hygrophorus ####
hygrophorus.rootbag <- root_bag %>% filter(Genus == "Hygrophorus")
hygrophorus.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = hygrophorus.rootbag)
summary(hygrophorus.mod1)
hist(resid(hygrophorus.mod1))
p.value.rootbag[6] <- coef(summary(hygrophorus.mod1))[1,5]

hygrophorus.bagtype <- bagtype %>% filter(Genus == "Hygrophorus")
hygrophorus.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = hygrophorus.bagtype)
summary(hygrophorus.mod2)
hist(resid(hygrophorus.mod2))
p.value.bagtype[6] <- coef(summary(hygrophorus.mod2))[1,5]

# _Lactarius ####
lactarius.rootbag <- root_bag %>% filter(Genus == "Lactarius")
lactarius.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = lactarius.rootbag)
summary(lactarius.mod1)
hist(resid(lactarius.mod1))
p.value.rootbag[7] <- coef(summary(lactarius.mod1))[1,5]

lactarius.bagtype <- bagtype %>% filter(Genus == "Lactarius")
lactarius.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = lactarius.bagtype)
summary(lactarius.mod2)
hist(resid(lactarius.mod2))
p.value.bagtype[7] <- coef(summary(lactarius.mod2))[1,5]

# _Piloderma ####
piloderma.rootbag <- root_bag %>% filter(Genus == "Piloderma")
piloderma.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = piloderma.rootbag)
summary(piloderma.mod1)
hist(resid(piloderma.mod1))
p.value.rootbag[8] <- coef(summary(piloderma.mod1))[1,5]

piloderma.bagtype <- bagtype %>% filter(Genus == "Piloderma")
piloderma.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = piloderma.bagtype)
summary(piloderma.mod2)
hist(resid(piloderma.mod2))
p.value.bagtype[8] <- coef(summary(piloderma.mod2))[1,5]

# _Pseudotomentella ####
pseudotomentella.rootbag <- root_bag %>% filter(Genus == "Pseudotomentella")
pseudotomentella.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = pseudotomentella.rootbag)
summary(pseudotomentella.mod1)
hist(resid(pseudotomentella.mod1))
p.value.rootbag[9] <- coef(summary(pseudotomentella.mod1))[1,5]

pseudotomentella.bagtype <- bagtype %>% filter(Genus == "Pseudotomentella")
pseudotomentella.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = pseudotomentella.bagtype)
summary(pseudotomentella.mod2)
hist(resid(pseudotomentella.mod2))
p.value.bagtype[9] <- coef(summary(pseudotomentella.mod2))[1,5]

# _Russula ####
russula.rootbag <- root_bag %>% filter(Genus == "Russula")
russula.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = russula.rootbag)
summary(russula.mod1)
hist(resid(russula.mod1))
p.value.rootbag[10] <- coef(summary(russula.mod1))[1,5]

russula.bagtype <- bagtype %>% filter(Genus == "Russula")
russula.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = russula.bagtype)
summary(russula.mod2)
hist(resid(russula.mod2))
p.value.bagtype[10] <- coef(summary(russula.mod2))[1,5]

# _Tomentella ####
tomentella.rootbag <- root_bag %>% filter(Genus == "Tomentella_Thelephora")
tomentella.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = tomentella.rootbag)
summary(tomentella.mod1)
hist(resid(tomentella.mod1))
p.value.rootbag[11] <- coef(summary(tomentella.mod1))[1,5]

tomentella.bagtype <- bagtype %>% filter(Genus == "Tomentella_Thelephora")
tomentella.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = tomentella.bagtype)
summary(tomentella.mod2)
hist(resid(tomentella.mod2))
p.value.bagtype[11] <- coef(summary(tomentella.mod2))[1,5]

# _Tylospora ####
tylospora.rootbag <- root_bag %>% filter(Genus == "Tylospora")
tylospora.mod1 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = tylospora.rootbag)
summary(tylospora.mod1)
hist(resid(tylospora.mod1))
p.value.rootbag[12] <- coef(summary(tylospora.mod1))[1,5]

tylospora.bagtype <- bagtype %>% filter(Genus == "Tylospora")
tylospora.mod2 <- lmerTest::lmer(Log_ratio ~ 1 + (1|Site), data = tylospora.bagtype)
summary(tylospora.mod2)
hist(resid(tylospora.mod2))
p.value.bagtype[12] <- coef(summary(tylospora.mod2))[1,5]


# __Export .txt with model p-values ####
# Use "p.adjust" to correct for multiple testing

genera <- c("Amanita", "Amphinema", "Cenococcum", "Cortinarius",
            "Hyaloscypha", "Hygrophorus", "Lactarius","Piloderma",
            "Pseudotomentella", "Russula", "Tomentella", "Tylospora")


bags.roots <- data.frame(Genus = genera, p.value = p.value.rootbag,
                         p.value.BH = p.adjust(p.value.rootbag, method="BH"))
bags.roots$sign <- ifelse(bags.roots$p.value.BH<=0.05, "*", "")
write.table(bags.roots, file = "bagroot_stats.txt", sep = ",", quote = FALSE, row.names = F)

bagtypes <- data.frame(Genus = genera, p.value = p.value.bagtype,
                         p.value.BH = p.adjust(p.value.bagtype, method="BH"))
bagtypes$sign <- ifelse(bagtypes$p.value.BH<=0.05, "*", "")
write.table(bagtypes, file = "bagtype_stats.txt", sep = ",", quote = FALSE, row.names = F)

# Test 3. Soil substrates ####
# PERMANOVA was used to test whether the different soil substrates in bags 
# affected ectomycorrhizal community composition

# make new df with soil bags only
soil.substrates <- abund.data %>% filter(Substrate %in% c("F", "RH", "P", "M"))
  str(soil.substrates)

soil.substrates$Site.no <- as.factor(soil.substrates$Site.no)
soil.substrates$Cafeteria.ID <- as.factor(soil.substrates$Cafeteria.ID)
soil.substrates$Caf <- as.factor(soil.substrates$Caf)
soil.substrates$Substrate <- as.factor(soil.substrates$Substrate)

# Calculate rel.abundance of exploration types
soil.substrates$C_type <- soil.substrates$Hyaloscypha+soil.substrates$Lactarius+soil.substrates$Russula
soil.substrates$SD_type <- soil.substrates$Cenococcum+soil.substrates$Hygrophorus+soil.substrates$Tylospora
soil.substrates$MDS_type <- soil.substrates$Amanita+soil.substrates$Pseudotomentella+soil.substrates$Tomentella_Thelephora
soil.substrates$MDF_type <- soil.substrates$Amphinema+soil.substrates$Cortinarius+soil.substrates$Piloderma

# Sqrt-transform relative abundances,
# this creates Hellinger transformed data
soil.substrates[,8:23] <- sqrt(soil.substrates[,8:23])

# Subset only sqrt-transformed abundances
spp <- soil.substrates[,8:19]
spp_expl <- soil.substrates[,20:23]
env <- soil.substrates[,1:5]

summary(env) 

# Run PERMANOVA, with Cafeteria.ID as strata to constrain permutations within each cafeteria

perm <- with(env, how(nperm = 1000, blocks = Cafeteria.ID))
perm

permanova.soil_expl <-adonis2(formula = spp_expl ~ Substrate+Site.no+Cafeteria.ID, data = env,
                         method = "bray", permutations = perm)
permanova.soil_expl # p=0.001

permanova.soil <-adonis2(formula = spp ~ Substrate+Site.no+Cafeteria.ID, data = env,
            method = "bray", permutations = perm)

permanova.soil # p=0.001

# Test 4. Single genus models ####
# Since test 3 showed that soil substrate type mattered
# for community composition, individual lme-models can be 
# run to see which genera prefer which substrate.
# p-values will be corrected for multiple testing (Benjamini Hochberg)

# Individual models for genus preference - one per genus (Hellinger transformed)

# Make vector to store p-values
p.values.soil <- numeric()

# Run individual models per species.
# _Amanita ####
aman.soil <- read.csv2("amanita_soil.csv", header = TRUE)
aman.soil$Site.no <- as.factor(aman.soil$Site.no)
aman.soil$Cafeteria.ID <- as.factor(aman.soil$Cafeteria.ID)
aman.soil$Caf <- as.factor(aman.soil$Caf)

aman.soil$abund.sqrt <- sqrt(aman.soil$Amanita)
aman.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Caf), data = aman.soil)
hist(resid(aman.soil.mod))
anova(aman.soil.mod)
p.values.soil[1] <- anova(aman.soil.mod)$`Pr(>F)`

# _Amphinema ####
amph.soil <- read.csv2("amphinema_soil.csv", header = TRUE)
amph.soil$Site.no <- as.factor(amph.soil$Site.no)
amph.soil$Cafeteria.ID <- as.factor(amph.soil$Cafeteria.ID)

amph.soil$abund.sqrt <- sqrt(amph.soil$Amphinema)
amph.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = amph.soil)
hist(resid(amph.soil.mod))
anova(amph.soil.mod) 
p.values.soil[2] <- anova(amph.soil.mod)$`Pr(>F)`

# _Cenococcum ####
ceno.soil <- read.csv2("cenococcum_soil.csv", header = TRUE)
ceno.soil$Site.no <- as.factor(ceno.soil$Site.no)
ceno.soil$Cafeteria.ID <- as.factor(ceno.soil$Cafeteria.ID)

ceno.soil$abund.sqrt <- sqrt(ceno.soil$Cenococcum)
ceno.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = ceno.soil)
hist(resid(ceno.soil.mod))
anova(ceno.soil.mod) 
p.values.soil[3] <- anova(ceno.soil.mod)$`Pr(>F)`

# _Cortinarius ####
cort.soil <- read.csv2("cortinarius_soil.csv", header = TRUE)
cort.soil$Site.no <- as.factor(cort.soil$Site.no)
cort.soil$Cafeteria.ID <- as.factor(cort.soil$Cafeteria.ID)

cort.soil$abund.sqrt <- sqrt(cort.soil$Cortinarius)
cort.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = cort.soil)
hist(resid(cort.soil.mod))
anova(cort.soil.mod) 
p.values.soil[4] <- anova(cort.soil.mod)$`Pr(>F)`

# _Hyaloscypha ####
hyal.soil <- read.csv2("hyaloscypha_soil.csv", header = TRUE)
hyal.soil$Site.no <- as.factor(hyal.soil$Site.no)
hyal.soil$Cafeteria.ID <- as.factor(hyal.soil$Cafeteria.ID)

hyal.soil$abund.sqrt <- sqrt(hyal.soil$Hyaloscypha)
hyal.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = hyal.soil)
hist(resid(hyal.soil.mod))
anova(hyal.soil.mod) 
p.values.soil[5] <- anova(hyal.soil.mod)$`Pr(>F)`

# _Hygrophorus ####
hyg.soil <- read.csv2("hygrophorus_soil.csv", header = TRUE)
hyg.soil$Site.no <- as.factor(hyg.soil$Site.no)
hyg.soil$Cafeteria.ID <- as.factor(hyg.soil$Cafeteria.ID)

hyg.soil$abund.sqrt <- sqrt(hyg.soil$Hygrophorus)
hyg.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = hyg.soil)
hist(resid(hyg.soil.mod))
anova(hyg.soil.mod) 
p.values.soil[6] <- anova(hyg.soil.mod)$`Pr(>F)`

# _Lactarius ####
lact.soil <- read.csv2("lactarius_soil.csv", header = TRUE)
lact.soil$Site.no <- as.factor(lact.soil$Site.no)
lact.soil$Cafeteria.ID <- as.factor(lact.soil$Cafeteria.ID)

lact.soil$abund.sqrt <- sqrt(lact.soil$Lactarius)
lact.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = lact.soil)
hist(resid(lact.soil.mod))
anova(lact.soil.mod) 
p.values.soil[7] <- anova(lact.soil.mod)$`Pr(>F)`

# _Piloderma ####
pilo.soil <- read.csv2("piloderma_soil.csv", header = TRUE)
pilo.soil$Site.no <- as.factor(pilo.soil$Site.no)
pilo.soil$Cafeteria.ID <- as.factor(pilo.soil$Cafeteria.ID)

pilo.soil$abund.sqrt <- sqrt(pilo.soil$Piloderma)
pilo.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = pilo.soil)
hist(resid(pilo.soil.mod))
anova(pilo.soil.mod) 
p.values.soil[8] <- anova(pilo.soil.mod)$`Pr(>F)`

# _Pseudotomentella ####
pseu.soil <- read.csv2("pseudotomentella_soil.csv", header = TRUE)
pseu.soil$Site.no <- as.factor(pseu.soil$Site.no)
pseu.soil$Cafeteria.ID <- as.factor(pseu.soil$Cafeteria.ID)

pseu.soil$abund.sqrt <- sqrt(pseu.soil$Pseudotomentella)
pseu.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = pseu.soil)
hist(resid(pseu.soil.mod))
anova(pseu.soil.mod) 
p.values.soil[9] <- anova(pseu.soil.mod)$`Pr(>F)`

# _Russula ####
russ.soil <- read.csv2("russula_soil.csv", header = TRUE)
russ.soil$Site.no <- as.factor(russ.soil$Site.no)
russ.soil$Cafeteria.ID <- as.factor(russ.soil$Cafeteria.ID)

russ.soil$abund.sqrt <- sqrt(russ.soil$Russula)
russ.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = russ.soil)
hist(resid(russ.soil.mod))
anova(russ.soil.mod) 
p.values.soil[10] <- anova(russ.soil.mod)$`Pr(>F)`

# _Tomentella ####
tomen.soil <- read.csv2("tomentella_soil.csv", header = TRUE)
tomen.soil$Site.no <- as.factor(tomen.soil$Site.no)
tomen.soil$Cafeteria.ID <- as.factor(tomen.soil$Cafeteria.ID)

tomen.soil$abund.sqrt <- sqrt(tomen.soil$Tomentella_Thelephora)
tomen.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = tomen.soil)
hist(resid(tomen.soil.mod))
anova(tomen.soil.mod)
summary(tomen.soil.mod)
p.values.soil[11] <- anova(tomen.soil.mod)$`Pr(>F)`

# _Tylospora ####
tylo.soil <- read.csv2("tylospora_soil.csv", header = TRUE)
tylo.soil$Site.no <- as.factor(tylo.soil$Site.no)
tylo.soil$Cafeteria.ID <- as.factor(tylo.soil$Cafeteria.ID)

tylo.soil$abund.sqrt <- sqrt(tylo.soil$Tylospora)
tylo.soil.mod <- lmerTest::lmer(abund.sqrt ~ Substrate + (1|Site.no/Cafeteria.ID), data = tylo.soil)
hist(resid(tylo.soil.mod))
anova(tylo.soil.mod) 
summary(tylo.soil.mod)
p.values.soil[12] <- anova(tylo.soil.mod)$`Pr(>F)`

genera <- c("Amanita", "Amphinema", "Cenococcum", "Cortinarius",
            "Hyaloscypha", "Hygrophorus", "Lactarius","Piloderma",
            "Pseudotomentella", "Russula", "Tomentella", "Tylospora")

substrates.summ <- data.frame(Genus = genera, p.value = p.values.soil,
                              p.value.BH = p.adjust(p.values.soil, method="BH"))
substrates.summ$sign <- ifelse(substrates.summ$p.value.BH<=0.05, "*", "")
write.table(substrates.summ, file="substrates_stats.txt", sep=",", quote=FALSE, rownames=FALSE)

# Test 5. Substrate preferences for genera #### 
# lme-models showed that Amphinema, Cenococcum and Tomentella 
# have substrate specific colonisation patterns

# Post-hoc test to see differences in preference
amph.posthoc <- emmeans(amph.soil.mod, list(pairwise ~ Substrate), adjust = "tukey")
ceno.posthoc <- emmeans(ceno.soil.mod, list(pairwise ~ Substrate), adjust = "tukey")
tomen.posthoc <- emmeans(tomen.soil.mod, list(pairwise ~ Substrate), adjust = "tukey")


# Test 6. Sand substrates ####
# PERMANOVA was used to test whether pure sand or apatite amended sand 
# affected ectomycorrhizal community composition

# make new df with sand bags only
sand.substrates <- abund.data %>% filter(Substrate %in% c("1_S", "2_A"))
str(sand.substrates)

sand.substrates$Site.no <- as.factor(sand.substrates$Site.no)
sand.substrates$Cafeteria.ID <- as.factor(sand.substrates$Cafeteria.ID)
sand.substrates$Caf <- as.factor(sand.substrates$Caf)
sand.substrates$Substrate <- as.factor(sand.substrates$Substrate)

# Sqrt-transform relative abundances,
# this creates Hellinger transformed data
sand.substrates[,8:19] <- sqrt(sand.substrates[,8:19])

spp.sand <- sand.substrates[,8:19]
env.sand <- sand.substrates[,1:5]


# Run PERMANOVA, with Cafeteria.ID as strata to constrain permutations within each cafeteria
perm2 <- with(env.sand, how(nperm = 1000,blocks = Cafeteria.ID))

permanova.sand <- adonis2(formula = spp.sand ~ Substrate+Site.no+Cafeteria.ID, data = env.sand,
                                           method = "bray", permutations = perm2)
permanova.sand # p=0.46

# No significant substrate effect -> no post-hoc testing

# Figure 1 ####
# Plot model estimates and SE as error bars for each genus 

estimate.rootbag <- numeric()
estimate.rootbag[1] <- coef(summary(amanita.mod1))[1,1]
estimate.rootbag[2] <- coef(summary(amphinema.mod1))[1,1]
estimate.rootbag[3] <- coef(summary(cenococcum.mod1))[1,1]
estimate.rootbag[4] <- coef(summary(cortinarius.mod1))[1,1]
estimate.rootbag[5] <- coef(summary(hyaloscypha.mod1))[1,1]
estimate.rootbag[6] <- coef(summary(hygrophorus.mod1))[1,1]
estimate.rootbag[7] <- coef(summary(lactarius.mod1))[1,1]
estimate.rootbag[8] <- coef(summary(piloderma.mod1))[1,1]
estimate.rootbag[9] <- coef(summary(pseudotomentella.mod1))[1,1]
estimate.rootbag[10] <- coef(summary(russula.mod1))[1,1]
estimate.rootbag[11] <- coef(summary(tomentella.mod1))[1,1]
estimate.rootbag[12] <- coef(summary(tylospora.mod1))[1,1]

SE.rootbag <- numeric()
SE.rootbag[1] <- coef(summary(amanita.mod1))[1,2]
SE.rootbag[2] <- coef(summary(amphinema.mod1))[1,2]
SE.rootbag[3] <- coef(summary(cenococcum.mod1))[1,2]
SE.rootbag[4] <- coef(summary(cortinarius.mod1))[1,2]
SE.rootbag[5] <- coef(summary(hyaloscypha.mod1))[1,2]
SE.rootbag[6] <- coef(summary(hygrophorus.mod1))[1,2]
SE.rootbag[7] <- coef(summary(lactarius.mod1))[1,2]
SE.rootbag[8] <- coef(summary(piloderma.mod1))[1,2]
SE.rootbag[9] <- coef(summary(pseudotomentella.mod1))[1,2]
SE.rootbag[10] <- coef(summary(russula.mod1))[1,2]
SE.rootbag[11] <- coef(summary(tomentella.mod1))[1,2]
SE.rootbag[12] <- coef(summary(tylospora.mod1))[1,2]

estimate.bagtype <- numeric()
estimate.bagtype[1] <- coef(summary(amanita.mod2))[1,1]
estimate.bagtype[2] <- coef(summary(amphinema.mod2))[1,1]
estimate.bagtype[3] <- coef(summary(cenococcum.mod2))[1,1]
estimate.bagtype[4] <- coef(summary(cortinarius.mod2))[1,1]
estimate.bagtype[5] <- coef(summary(hyaloscypha.mod2))[1,1]
estimate.bagtype[6] <- coef(summary(hygrophorus.mod2))[1,1]
estimate.bagtype[7] <- coef(summary(lactarius.mod2))[1,1]
estimate.bagtype[8] <- coef(summary(piloderma.mod2))[1,1]
estimate.bagtype[9] <- coef(summary(pseudotomentella.mod2))[1,1]
estimate.bagtype[10] <- coef(summary(russula.mod2))[1,1]
estimate.bagtype[11] <- coef(summary(tomentella.mod2))[1,1]
estimate.bagtype[12] <- coef(summary(tylospora.mod2))[1,1]

SE.bagtype <- numeric()
SE.bagtype[1] <- coef(summary(amanita.mod2))[1,2]
SE.bagtype[2] <- coef(summary(amphinema.mod2))[1,2]
SE.bagtype[3] <- coef(summary(cenococcum.mod2))[1,2]
SE.bagtype[4] <- coef(summary(cortinarius.mod2))[1,2]
SE.bagtype[5] <- coef(summary(hyaloscypha.mod2))[1,2]
SE.bagtype[6] <- coef(summary(hygrophorus.mod2))[1,2]
SE.bagtype[7] <- coef(summary(lactarius.mod2))[1,2]
SE.bagtype[8] <- coef(summary(piloderma.mod2))[1,2]
SE.bagtype[9] <- coef(summary(pseudotomentella.mod2))[1,2]
SE.bagtype[10] <- coef(summary(russula.mod2))[1,2]
SE.bagtype[11] <- coef(summary(tomentella.mod2))[1,2]
SE.bagtype[12] <- coef(summary(tylospora.mod2))[1,2]


expl.types <- c("MDS", "MDF", "SD", "MDF", "C", "SD", "C", "MDF",
                "MDS", "C", "MDS", "SD")
#Mean abundance of genera, used for point size in Fig. 1
adj.means <- read.csv2("Adjusted.means.csv", header=TRUE)

fig1.dat <- tibble(Genus = genera, expl.types, estimate.rootbag, SE.rootbag, estimate.bagtype, SE.bagtype)
fig1.dat <- inner_join(fig1.dat, adj.means, by = "Genus")


ggplot(fig1.dat, aes(x = estimate.rootbag, xmax = estimate.rootbag+SE.rootbag, xmin = estimate.rootbag-SE.rootbag,
                     y = estimate.bagtype, ymax = estimate.bagtype+SE.bagtype, ymin = estimate.bagtype-SE.bagtype,
                     color = expl.types))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey50", size=0.4)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size =0.4)+
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_text(size=13),
    legend.text = element_text(size =11)
  )+
  labs(x = "Log-ratio (ingrowth bags / roots)",
       y = "Log-ratio (soil bags / sand bags)")+
  geom_point(aes(size=Adj.mean))+
  scale_size_continuous(range=c(1,7))+
  geom_errorbar(width = 0.1)+
  geom_errorbarh(width = 0.1)+
  geom_text(label = fig1.dat$Genus, nudge_x = 1, nudge_y = -0.1, check_overlap = F, color = "black")+
  scale_color_manual(name = "Exploration type",
                     values=c("#034363","#18878B", "#A35257", "#E07800"),
                     breaks=c("C", "SD", "MDS", "MDF"),
                     labels= c("C", "SD", "MDS", "MDF"))


# The figure was edited in Adobe Illustrator to add symbols where lme-models 
# were significant + improve the general appearance.

# Figure 2 ####

#In Figure 2, the mean relative abundance +/- SE of each genus is plotted.
#  F = organic topsoil, low N
#  M = mull soil
#  RH = organic topsoil, high N
#  P = organic topsoil, P fertilised

M.dat <- abund.data %>% filter(Substrate == "M")
F.dat <- abund.data %>% filter(Substrate== "F")
RH.dat <- abund.data %>% filter(Substrate == "RH")
P.dat <- abund.data %>% filter(Substrate == "P")

# Gather data into long dataset to be able to calculate mean abundance
M_long <- gather(M.dat, "Genus", "Rel_abund", 8:19)
F_long <- gather(F.dat, "Genus", "Rel_abund", 8:19)
RH_long <- gather(RH.dat, "Genus", "Rel_abund", 8:19)
P_long <- gather(P.dat, "Genus", "Rel_abund", 8:19)

# Calculate mean rel abundance on roots, per genus

M.mean <- M_long %>%
  group_by(Genus) %>%
  summarise( 
    n.M=n(),
    mean.M=mean(Rel_abund, na.rm = TRUE),
    sd.M=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.M=sd.M/sqrt(n.M))

F.mean <- F_long %>%
  group_by(Genus) %>%
  summarise( 
    n.F=n(),
    mean.F=mean(Rel_abund, na.rm = TRUE),
    sd.F=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.F=sd.F/sqrt(n.F))

RH.mean <- RH_long %>%
  group_by(Genus) %>%
  summarise( 
    n.RH=n(),
    mean.RH=mean(Rel_abund, na.rm = TRUE),
    sd.RH=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.RH=sd.RH/sqrt(n.RH))

P.mean <- P_long %>%
  group_by(Genus) %>%
  summarise( 
    n.P=n(),
    mean.P=mean(Rel_abund, na.rm = TRUE),
    sd.P=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.P=sd.P/sqrt(n.P))


substrate.summ <- left_join(M.mean, F.mean, by="Genus")
substrate.summ <- left_join(substrate.summ, RH.mean, by="Genus")
substrate.summ <- left_join(substrate.summ, P.mean, by="Genus")


substrate.mean.vals <- tibble(matrix(nrow=12))
substrate.mean.vals$Genus <- substrate.summ$Genus
substrate.mean.vals$F <- substrate.summ$mean.F
substrate.mean.vals$RH <- substrate.summ$mean.RH
substrate.mean.vals$P <- substrate.summ$mean.P
substrate.mean.vals$M <- substrate.summ$mean.M
substrate.mean.vals <- substrate.mean.vals[,-1]

# Calculate rel.abund of "other" (not the 12 most frequent) genera to include in fig. 2

fig2.dat <-substrate.mean.vals

F.other <- 1-sum(fig2.dat[1:12,2])
RH.other <- 1-sum(fig2.dat[1:12,3])
P.other <- 1-sum(fig2.dat[1:12,4])
M.other <- 1-sum(fig2.dat[1:12,5])

fig2.dat[13,1] <- "Other"

fig2.dat[13,2] <- F.other
fig2.dat[13,3] <- RH.other
fig2.dat[13,4] <- P.other
fig2.dat[13,5] <- M.other

fig2.dat.long <- gather(fig2.dat, "Substrate", "Mean", 2:5)

fig2.dat.long <- fig2.dat.long %>%
  mutate(order = case_when(
    Genus == "Amanita" ~ "h",
    Genus == "Amphinema" ~ "k",
    Genus == "Cenococcum" ~ "e",
    Genus == "Cortinarius" ~ "l",
    Genus == "Hyaloscypha" ~ "b",
    Genus == "Hygrophorus" ~ "f",
    Genus == "Lactarius" ~ "c",
    Genus == "Piloderma" ~ "m",
    Genus == "Pseudotomentella" ~ "i",
    Genus == "Russula" ~ "d",
    Genus == "Tomentella_Thelephora" ~ "j",
    Genus == "Tylospora" ~ "g",
    Genus == "Other" ~ "a"
  ))

fig2.dat.long$genus.order <- str_c(fig2.dat.long$order,"_", fig2.dat.long$Genus)

fig2.dat.long <- fig2.dat.long %>%
  mutate(subs.order = case_when(
    Substrate == "F" ~ "1_F",
    Substrate == "RH" ~ "2_RH",
    Substrate == "P" ~ "3_P",
    Substrate == "M" ~ "4_M",
  ))

fig2.dat.long$Mean <- fig2.dat.long$Mean*100

ggplot(fig2.dat.long, aes(x= subs.order, y = Mean, fill = genus.order))+
  theme_classic2()+
  geom_bar(stat="identity")+
  theme(legend.position = "right")+
  scale_fill_manual(name="Genus",
                    values=c("grey80", "#034363","#056B9E","#0793DA","#18878B","#21BBC0","#3FD9DE",
                             "#A35257","#BB777C","#CFA0A3","#E07800","#FF961F","#FFB35C"),
                    breaks=c("a_Other","b_Hyaloscypha", "c_Lactarius", "d_Russula",
                             "e_Cenococcum", "f_Hygrophorus", "g_Tylospora",
                             "h_Amanita", "i_Pseudotomentella", "j_Tomentella_Thelephora",
                             "k_Amphinema", "l_Cortinarius", "m_Piloderma"),
                    labels=c("Other", "Hyaloscypha", "Lactarius", "Russula",
                             "Cenococcum", "Hygrophorus", "Tylospora",
                             "Amanita", "Pseudotomentella", "Tomentella",
                             "Amphinema", "Cortinarius", "Piloderma"))



