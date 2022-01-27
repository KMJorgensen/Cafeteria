# R-script to reproduce analyses for:
#"Mycelial foraging strategies and phenotypic plasticity of ectomycorrhizal fungi in boreal Picea abies forests"
# Jörgensen, K., Clemmensen, KE., Wallander, H., Lindahl, BD. 


# Contact about R script: karolina.jorgensen@slu.se                          

# This has been tested under:
# R version 4.0.1 (2022-01-27)


# Required packages
# Statistical analyses
library(lme4) # ver 1.1.26
library(lmerTest) # ver 3.1.3
library(stats) # ver 4.0.1
library(emmeans) # ver 1.5.3

# Plotting and data handling
library(dplyr) # ver 1.0.3
library(tidyr) # ver 1.1.2
library(ggplot2) # ver 3.3.3
library(sjPlot) # ver 2.8.7 
library(ggpubr) # ver 0.4.0

setwd("~/docs/Dokument/Projekt/VR-projekt/Cafeteriaexperiment/Submission_prep")

# Load dataset with untransformed relative abundance of ECM community of 
# 15 most frequent genera in cafeterias.
# Abbreviations for substrate types:
#  F = organic topsoil, low N
#  M = mull soil
#  RH = organic topsoil, high N
#  P = organic topsoil, P fertilised
#  S = sand
#  A = sand + apatite

data <- read.csv2("Rel_abundances.csv", header = T)
data$Site.no <- as.character(data$Site.no) # Site ID
data$Cafeteria.ID <- as.character(data$Cafeteria.ID) # Cafeteria ID
data$Caf <- as.character(data$Caf) # Cafeteria replicate, within site
data$Substrate <- as.factor(data$Substrate) # Either root or one of six bag substrates
data$Bagtype <- as.factor(data$Bagtype) # Either root, sandbag (S, A) or soilbag (F, M, RH, P)
data$Bag_root <- as.factor(data$Bag_root) # Root or ingrowth bag

data[is.na(data)]<-0 # Set abundance to 0 in bags where ECM genera are not present
colnames(data) # fungal genera in columns 8:21

# Square-root transform relative abundances and inorganic N.
# This creates Hellinger-transformed community data
# This dataset is used for the lme-models for ECM genera

data.sq <- data
data.sq[,8:21] <- sqrt(data.sq[,8:21])

# Analyses: Mixed models ####
# Strategy is to:
#  Use sqrt-transformed (Hellinger transformed) relative abundances of each genus. 
#  Use nested structure to account for random effects (1|Site / Cafeteria)
#  Use Benjamini & Hochberg correction to account for multiple testing (15 genera)

# Make "Root" the reference treatment in bag vs root analysis
contrasts(data.sq$Bag_root) <- contr.treatment(levels(data.sq$Bag_root),
                                            base = which(levels(data.sq$Bag_root)== "Root"))

# Goal is to investigate whether the relative abundance of a genus differ between: 
# 1. Bags and roots? - models called NN.bag
# 2. Soil bags and sand bags? - models called NN.bagtype
# 3. Apatite bags and sand bags? - models called NN.sand
# 4. Different soil bags? - models called NN.soil

# For differences between soil substrates, a posthoc Tukey test was done using the "emmeans" package

# For each genus, only cafeterias where the genus is present is filtered out 
# this filtering creates an individual dataset for each genus e.g. Amphinema = "amph.data"

# Make vectors to store p.values for each analysis - this is used to make the p-value corrections
p.value.bag <- numeric()
p.value.bagtype <- numeric()
p.value.sand <- numeric()
p.value.soil <- numeric()

# Models ####

# ___Amphinema ####

# Filter out cafeterias where genus is present
amph.data <- data.sq %>% filter(Cafeteria.ID %in% c("9","14","16","17","20","22","23","24","35","37",
                                                 "38","40","47","50","1","2","4","8","10","13","18",
                                                 "19","21","32","39","43","45","48","25","31","33",
                                                 "36","3","12","15","28","41","6","7","29","30","42"
                                                 ,"5","11","26","27","44"))

amphinema.bag <- lmerTest::lmer(Amphinema ~ Bag_root + (1|Site.no/Cafeteria.ID), data = amph.data)
summary(amphinema.bag)
hist(resid(amphinema.bag))
plot(amphinema.bag)
p.value.bag[1] <- coef(summary(amphinema.bag))[2,5]

amph.bagtype <- amph.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
amphinema.bagtype <- lmerTest::lmer(Amphinema ~ Bagtype + (1|Site.no/Cafeteria.ID), data = amph.bagtype)
summary(amphinema.bagtype)
hist(resid(amphinema.bagtype))
p.value.bagtype[1] <- coef(summary(amphinema.bagtype))[2,5]

amph.sand <- amph.data %>% filter(Substrate %in% c("2_A", "1_S"))
amphinema.sand <- lmerTest::lmer(Amphinema ~ Substrate + (1|Site.no/Cafeteria.ID), data = amph.sand)
summary(amphinema.sand)
hist(resid(amphinema.sand))
p.value.sand[1] <- coef(summary(amphinema.sand))[2,5]

amph.soil <- amph.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
amphinema.soil <- lmerTest::lmer(Amphinema ~ Substrate + (1|Site.no/Cafeteria.ID), data = amph.soil)
anova(amphinema.soil)
summary(amphinema.soil)
hist(resid(amphinema.soil))
p.value.soil[1] <- anova(amphinema.soil)$`Pr(>F)`
amph <- emmeans(amphinema.soil, list(pairwise ~ Substrate), adjust = "tukey")

 

#___Tomentella (incl. Thelephora) ####

toth.data <- data.sq %>% filter(Cafeteria.ID %in% c("1","2","4","9","13","16","20","22","24","28","35",
                                                 "37","47","50","11","12","14","17","18","21","23",
                                                 "26","29","36","38","39","43","6","7","8","10","15",
                                                 "25","33","40","41","42","44","45","3","5","27","32",
                                                 "48","19","31","30"))

tomentella_thelephora.bag <- lmerTest::lmer(Tomentella_Thelephora ~ Bag_root + (1|Site.no/Cafeteria.ID), data = toth.data)
summary(tomentella_thelephora.bag)
hist(resid(tomentella_thelephora.bag))
p.value.bag[2] <- coef(summary(tomentella_thelephora.bag))[2,5]

toth.bagtype <- toth.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
tomentella_thelephora.bagtype <- lmerTest::lmer(Tomentella_Thelephora ~ Bagtype + (1|Site.no/Cafeteria.ID), data = toth.bagtype)
summary(tomentella_thelephora.bagtype)
hist(resid(tomentella_thelephora.bagtype))
p.value.bagtype[2] <- coef(summary(tomentella_thelephora.bagtype))[2,5]

toth.sand <- toth.data %>% filter(Substrate %in% c("2_A", "1_S"))
tomentella_thelephora.sand <- lmerTest::lmer(Tomentella_Thelephora ~ Substrate + (1|Site.no/Cafeteria.ID), data = toth.sand)
summary(tomentella_thelephora.sand)
hist(resid(tomentella_thelephora.sand))
p.value.sand[2] <- coef(summary(tomentella_thelephora.sand))[2,5]

toth.soil <- toth.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
tomentella_thelephora.soil <- lmerTest::lmer(Tomentella_Thelephora ~ Substrate + (1|Site.no/Cafeteria.ID), data = toth.soil)
anova(tomentella_thelephora.soil)
summary(tomentella_thelephora.soil)
hist(resid(tomentella_thelephora.soil))
p.value.soil[2] <- anova(tomentella_thelephora.soil)$`Pr(>F)`
tom_thel <- emmeans(tomentella_thelephora.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Cenococcum ####
ceno.data <- data.sq %>% filter(Cafeteria.ID %in% c("17","22","29","47","8","9","13","14","20","32","2",
                                                 "18","23","27","28","41","44","1","3","4","7","16",
                                                 "21","24","25","31","33","35","42","45","48","50","6",
                                                 "10","11","12","19","26","36","37","39","43","5","15",
                                                 "30","38","40"))

cenococcum.bag <- lmerTest::lmer(Cenococcum ~ Bag_root + (1|Site.no/Cafeteria.ID), data = ceno.data)
summary(cenococcum.bag)
hist(resid(cenococcum.bag))
p.value.bag[3] <- coef(summary(cenococcum.bag))[2,5]

ceno.bagtype <- ceno.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
cenococcum.bagtype <- lmerTest::lmer(Cenococcum ~ Bagtype + (1|Site.no/Cafeteria.ID), data = ceno.bagtype)
summary(cenococcum.bagtype)
hist(resid(cenococcum.bagtype))
p.value.bagtype[3] <- coef(summary(cenococcum.bagtype))[2,5]

ceno.sand <- ceno.data %>% filter(Substrate %in% c("2_A", "1_S"))
cenococcum.sand <- lmerTest::lmer(Cenococcum ~ Substrate + (1|Site.no/Cafeteria.ID), data = ceno.sand)
summary(cenococcum.sand)
hist(resid(cenococcum.sand))
p.value.sand[3] <- coef(summary(cenococcum.sand))[2,5]

ceno.soil <- ceno.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
cenococcum.soil <- lmerTest::lmer(Cenococcum ~ Substrate + (1|Site.no/Cafeteria.ID), data = ceno.soil)
anova(cenococcum.soil)
summary(cenococcum.soil)
hist(resid(cenococcum.soil))
p.value.soil[3] <- anova(cenococcum.soil)$`Pr(>F)`
ceno <- emmeans(cenococcum.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Russula ####
russ.data <- data.sq %>% filter(Cafeteria.ID %in% c("24","50","14","16","17","22","28","11","13","15","20","25",
                                                 "37","38","41","43","44","1","2","5","6","10","18","21","26",
                                                 "30","32","36","3","7","8","12","19","27","35","39","40","42",
                                                 "47","4","9","23","33","45","48","29","31"))

russula.bag <- lmerTest::lmer(Russula ~ Bag_root + (1|Site.no/Cafeteria.ID), data = russ.data)
summary(russula.bag)
plot(russula.bag)
hist(resid(russula.bag))
p.value.bag[4] <- coef(summary(russula.bag))[2,5]

russ.bagtype <- russ.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
russula.bagtype <- lmerTest::lmer(Russula ~ Bagtype + (1|Site.no/Cafeteria.ID), data = russ.bagtype)
summary(russula.bagtype)
hist(resid(russula.bagtype))
p.value.bagtype[4] <- coef(summary(russula.bagtype))[2,5]

russ.sand <- russ.data %>% filter(Substrate %in% c("2_A", "1_S"))
russula.sand <- lmerTest::lmer(Russula ~ Substrate + (1|Site.no/Cafeteria.ID), data = russ.sand)
summary(russula.sand)
hist(resid(russula.sand))
p.value.sand[4] <- coef(summary(russula.sand))[2,5]

russ.soil <- russ.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
russula.soil <- lmerTest::lmer(Russula ~ Substrate + (1|Site.no/Cafeteria.ID), data = russ.soil)
anova(russula.soil)
summary(russula.soil)
hist(resid(russula.soil))
p.value.soil[4] <- anova(russula.soil)$`Pr(>F)`
russ <- emmeans(russula.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Piloderma ####
pilo.data <- data.sq %>% filter(Cafeteria.ID %in% c("24","50","14","16","17","22","28","11","13","15","20","25",
                                                 "37","38","41","43","44","1","2","5","6","10","18","21","26",
                                                 "30","32","36","3","7","8","12","19","27","35","39","40","42",
                                                 "47","4","9","23","33","45","48","29","31"))

piloderma.bag <- lmerTest::lmer(Piloderma ~ Bag_root + (1|Site.no/Cafeteria.ID), data = pilo.data)
summary(piloderma.bag)
hist(resid(piloderma.bag))
p.value.bag[5] <- coef(summary(piloderma.bag))[2,5]

pilo.bagtype <- pilo.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
piloderma.bagtype <- lmerTest::lmer(Piloderma ~ Bagtype + (1|Site.no/Cafeteria.ID), data = pilo.bagtype)
summary(piloderma.bagtype)
hist(resid(piloderma.bagtype))
p.value.bagtype[5] <- coef(summary(piloderma.bagtype))[2,5]

pilo.sand <- pilo.data %>% filter(Substrate %in% c("2_A", "1_S"))
piloderma.sand <- lmerTest::lmer(Piloderma ~ Substrate + (1|Site.no/Cafeteria.ID), data = pilo.sand)
summary(piloderma.sand)
hist(resid(piloderma.sand))
p.value.sand[5] <- coef(summary(piloderma.sand))[2,5]

pilo.soil <- pilo.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
piloderma.soil <- lmerTest::lmer(Piloderma ~ Substrate + (1|Site.no/Cafeteria.ID), data = pilo.soil)
anova(piloderma.soil)
summary(piloderma.soil)
p.value.soil[5] <- anova(piloderma.soil)$`Pr(>F)`
pilo <- emmeans(piloderma.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Tylospora ####
tylo.data <- data.sq %>% filter(Cafeteria.ID %in% c("18","39","40","50","13","16","20","26","27","30","31",
                                                 "38","1","4","29","35","41","43","44","2","15","17","21",
                                                 "28","36","37","48","10","12","19","23","45","5","6","7",
                                                 "11","14","24","32","42","3","9","22","25","33","47"))

tylospora.bag <- lmerTest::lmer(Tylospora ~ Bag_root + (1|Site.no/Cafeteria.ID), data = tylo.data)
summary(tylospora.bag)
hist(resid(tylospora.bag))
p.value.bag[6] <- coef(summary(tylospora.bag))[2,5]

tylo.bagtype <- tylo.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
tylospora.bagtype <- lmerTest::lmer(Tylospora ~ Bagtype + (1|Site.no/Cafeteria.ID), data = tylo.bagtype)
summary(tylospora.bagtype)
hist(resid(tylospora.bagtype))
p.value.bagtype[6] <- coef(summary(tylospora.bagtype))[2,5]

tylo.sand <- tylo.data %>% filter(Substrate %in% c("2_A", "1_S"))
tylospora.sand <- lmerTest::lmer(Tylospora ~ Substrate + (1|Site.no/Cafeteria.ID), data = tylo.sand)
summary(tylospora.sand)
hist(resid(tylospora.sand))
p.value.sand[6] <- coef(summary(tylospora.sand))[2,5]

tylo.soil <- tylo.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
tylospora.soil <- lmerTest::lmer(Tylospora ~ Substrate + (1|Site.no/Cafeteria.ID), data = tylo.soil)
anova(tylospora.soil)
summary(tylospora.soil)
hist(resid(tylospora.soil))
p.value.soil[6] <- anova(tylospora.soil)$`Pr(>F)`
tylo <- emmeans(tylospora.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Cortinarius ####
cort.data <- data.sq %>% filter(Cafeteria.ID %in% c("35","1","3","4","5","6","20","28","43","2","7","10",
                                                 "25","29","32","36","37","39","44","12","22","24","27",
                                                 "30","33","38","42","8","9","13","15","18","19","23",
                                                 "31","47","48","50","17","21","26","40","41"))

cortinarius.bag <- lmerTest::lmer(Cortinarius ~ Bag_root + (1|Site.no/Cafeteria.ID), data = cort.data)
summary(cortinarius.bag)
hist(resid(cortinarius.bag))
p.value.bag[7] <- coef(summary(cortinarius.bag))[2,5]

cort.bagtype <- cort.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
cortinarius.bagtype <- lmerTest::lmer(Cortinarius ~ Bagtype + (1|Site.no/Cafeteria.ID), data = cort.bagtype)
summary(cortinarius.bagtype)
hist(resid(cortinarius.bagtype))
p.value.bagtype[7] <- coef(summary(cortinarius.bagtype))[2,5]

cort.sand <- cort.data %>% filter(Substrate %in% c("2_A", "1_S"))
cortinarius.sand <- lmerTest::lmer(Cortinarius ~ Substrate + (1|Site.no/Cafeteria.ID), data = cort.sand)
summary(cortinarius.sand)
hist(resid(cortinarius.sand))
p.value.sand[7] <- coef(summary(cortinarius.sand))[2,5]

cort.soil <- cort.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
cortinarius.soil <- lmerTest::lmer(Cortinarius ~ Substrate + (1|Site.no/Cafeteria.ID), data = cort.soil)
anova(cortinarius.soil)
summary(cortinarius.soil)
hist(resid(cortinarius.soil))
p.value.soil[7] <- anova(cortinarius.soil)$`Pr(>F)`
cort <- emmeans(cortinarius.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Hyaloscypha ####
hyal.data <- data.sq %>% filter(Cafeteria.ID %in% c("38","50","2","4","24","25","39","1","3","5","11",
                                                 "21","23","28","32","35","36","40","41","42","43",
                                                 "44","48","6","7","12","13","14","18","20","27","29",
                                                 "30","45","47","10","15","17","22","26","31","37"))

hyaloscypha.bag <- lmerTest::lmer(Hyaloscypha ~ Bag_root + (1|Site.no/Cafeteria.ID), data = hyal.data)
summary(hyaloscypha.bag)
hist(resid(hyaloscypha.bag))
p.value.bag[8] <- coef(summary(hyaloscypha.bag))[2,5]

hyal.bagtype <- hyal.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
hyaloscypha.bagtype <- lmerTest::lmer(Hyaloscypha ~ Bagtype + (1|Site.no/Cafeteria.ID), data = hyal.bagtype)
summary(hyaloscypha.bagtype)
hist(resid(hyaloscypha.bagtype))
p.value.bagtype[8] <- coef(summary(hyaloscypha.bagtype))[2,5]

hyal.sand <- hyal.data %>% filter(Substrate %in% c("2_A", "1_S"))
hyaloscypha.sand <- lmerTest::lmer(Hyaloscypha ~ Substrate + (1|Site.no/Cafeteria.ID), data = hyal.sand)
summary(hyaloscypha.sand)
hist(resid(hyaloscypha.sand))
p.value.sand[8] <- coef(summary(hyaloscypha.sand))[2,5]

hyal.soil <- hyal.data %>% filter(Substrate %in% c("RH", "F", "P", "M")) 
hyaloscypha.soil <- lmerTest::lmer(Hyaloscypha ~ Substrate + (1|Site.no/Cafeteria.ID), data = hyal.soil)
anova(hyaloscypha.soil)
summary(hyaloscypha.soil)
hist(resid(hyaloscypha.soil))
p.value.soil[8] <- anova(hyaloscypha.soil)$`Pr(>F)`
hyal <- emmeans(hyaloscypha.soil, list(pairwise ~ Substrate), adjust = "tukey")


# ___Pseudotomentella ####
pseu.data <- data.sq %>% filter(Cafeteria.ID %in% c("18","20","39","50","8","16","36","9","14","44","19","21",
                                                 "33","38","2","4","7","10","11","17","24","41","43","5","6",
                                                 "13","15","25","30","35","37","40","42","45","47","48")) 

pseudotomentella.bag <- lmerTest::lmer(Pseudotomentella ~ Bag_root + (1|Site.no/Cafeteria.ID), data = pseu.data)
summary(pseudotomentella.bag)
hist(resid(pseudotomentella.bag))
p.value.bag[9] <- coef(summary(pseudotomentella.bag))[2,5]

pseu.bagtype <- pseu.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
pseudotomentella.bagtype <- lmerTest::lmer(Pseudotomentella ~ Bagtype + (1|Site.no/Cafeteria.ID), data = pseu.bagtype)
summary(pseudotomentella.bagtype)
hist(resid(pseudotomentella.bagtype))
p.value.bagtype[9] <- coef(summary(pseudotomentella.bagtype))[2,5]

pseu.sand <- pseu.data %>% filter(Substrate %in% c("2_A", "1_S"))
pseudotomentella.sand <- lmerTest::lmer(Pseudotomentella ~ Substrate + (1|Site.no/Cafeteria.ID), data = pseu.sand)
summary(pseudotomentella.sand)
hist(resid(pseudotomentella.sand))
p.value.sand[9] <- coef(summary(pseudotomentella.sand))[2,5]

pseu.soil <- pseu.data %>% filter(Substrate %in% c("RH", "F", "P", "M"))
pseudotomentella.soil <- lmerTest::lmer(Pseudotomentella ~ Substrate + (1|Site.no/Cafeteria.ID), data = pseu.soil)
anova(pseudotomentella.soil)
summary(pseudotomentella.soil)
hist(resid(pseudotomentella.soil))
p.value.soil[9] <- anova(pseudotomentella.soil)$`Pr(>F)`
pseudo <- emmeans(pseudotomentella.soil, list(pairwise ~ Substrate), adjust = "tukey")


# ___Lactarius ####
lact.data <- data.sq %>% filter(Cafeteria.ID %in% c("40","31","39","5","6","23","41","43","44","1","2","7",
                                                 "9","10","14","24","29","32","33","38","50","3","4","8",
                                                 "16","18","20","21","22","37","48","12","13","15","17",
                                                 "25","26","28","36","42","45","47")) 

lactarius.bag <- lmerTest::lmer(Lactarius ~ Bag_root + (1|Site.no/Cafeteria.ID), data = lact.data)
summary(lactarius.bag)
hist(resid(lactarius.bag))
p.value.bag[10] <- coef(summary(lactarius.bag))[2,5]

lact.bagtype <- lact.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
lactarius.bagtype <- lmerTest::lmer(Lactarius ~ Bagtype + (1|Site.no/Cafeteria.ID), data = lact.bagtype)
summary(lactarius.bagtype)
hist(resid(lactarius.bagtype))
p.value.bagtype[10] <- coef(summary(lactarius.bagtype))[2,5]

lact.soil <- lact.data %>% filter(Substrate %in% c("RH", "F", "P", "M"))
lactarius.soil <- lmerTest::lmer(Lactarius ~ Substrate + (1|Site.no/Cafeteria.ID), data = lact.soil)
anova(lactarius.soil)
summary(lactarius.soil)
hist(resid(lactarius.soil))
p.value.soil[10] <- anova(lactarius.soil)$`Pr(>F)`

lact <- emmeans(lactarius.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Amanita ####
aman.data <- data.sq %>% filter(Cafeteria.ID %in% c("27","5","6","2","7","12","26","38","4","11","23",
                                                 "24","29","30","31","40","41","13","15","21","22",
                                                 "28","32","33","36","43","1","3","8","14","20","25",
                                                 "35","37","39","42","44","47")) 

amanita.bag <- lmerTest::lmer(Amanita ~ Bag_root + (1|Site.no/Cafeteria.ID), data = data)
summary(amanita.bag)
hist(resid(amanita.bag))
p.value.bag[11] <- coef(summary(amanita.bag))[2,5]

aman.bagtype <- aman.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
amanita.bagtype <- lmerTest::lmer(Amanita ~ Bagtype + (1|Site.no/Cafeteria.ID), data = aman.bagtype)
summary(amanita.bagtype)
hist(resid(amanita.bagtype))
p.value.bagtype[11] <- coef(summary(amanita.bagtype))[2,5]

aman.sand <- aman.data %>% filter(Substrate %in% c("2_A", "1_S"))
amanita.sand <- lmerTest::lmer(Amanita ~ Substrate + (1|Site.no/Cafeteria.ID), data = aman.sand)
summary(amanita.sand)
hist(resid(amanita.sand))
p.value.sand[10] <- coef(summary(amanita.sand))[2,5]

aman.soil <- aman.data %>% filter(Substrate %in% c("RH", "F", "P", "M"))
amanita.soil <- lmerTest::lmer(Amanita ~ Substrate + (1|Site.no/Cafeteria.ID), data = aman.soil)
anova(amanita.soil)
summary(amanita.soil)
hist(resid(amanita.soil))
p.value.soil[11] <- anova(amanita.soil)$`Pr(>F)`
amanita <- emmeans(amanita.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Tomentellopsis ####
topsis.data <- data.sq %>% filter(Cafeteria.ID %in% c("27","29","10","20","35","9","11","24","30","32",
                                                   "33","1","3","26","8","14","22","36","50","4","17",
                                                   "21","23","31","39","40","41","42","47")) 

tomentellopsis.bag <- lmerTest::lmer(Tomentellopsis ~ Bag_root + (1|Site.no/Cafeteria.ID), data = topsis.data)
summary(tomentellopsis.bag)
hist(resid(tomentellopsis.bag))
p.value.bag[12] <- coef(summary(tomentellopsis.bag))[2,5]

topsis.bagtype <- topsis.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
tomentellopsis.bagtype <- lmerTest::lmer(Tomentellopsis ~ Bagtype + (1|Site.no/Cafeteria.ID), data = topsis.bagtype)
summary(tomentellopsis.bagtype)
hist(resid(tomentellopsis.bagtype))
p.value.bagtype[12] <- coef(summary(tomentellopsis.bagtype))[2,5]

topsis.sand <- topsis.data %>% filter(Substrate %in% c("2_A", "1_S"))
tomentellopsis.sand <- lmerTest::lmer(Tomentellopsis ~ Substrate + (1|Site.no/Cafeteria.ID), data = topsis.sand)
summary(tomentellopsis.sand)
hist(resid(tomentellopsis.sand))
p.value.sand[11] <- coef(summary(tomentellopsis.sand))[2,5]

topsis.soil <- topsis.data %>% filter(Substrate %in% c("RH", "F", "P", "M"))
tomentellopsis.soil <- lmerTest::lmer(Tomentellopsis ~ Substrate + (1|Site.no/Cafeteria.ID), data = topsis.soil)
anova(tomentellopsis.soil)
summary(tomentellopsis.soil)
hist(resid(tomentellopsis.soil))
p.value.soil[12] <- anova(tomentellopsis.soil)$`Pr(>F)`
topsis <- emmeans(tomentellopsis.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Trichophaea ####
trich.data <- data.sq %>% filter(Cafeteria.ID %in% c("20","47","18","39","40","8","50","9","36","48","41",
                                                  "1","12","17","19","22","24","27","32","38","42","43")) 

trichophaea.bag <- lmerTest::lmer(Trichophaea ~ Bag_root + (1|Site.no/Cafeteria.ID), data = trich.data)
summary(trichophaea.bag)
hist(resid(trichophaea.bag))
p.value.bag[13] <- coef(summary(trichophaea.bag))[2,5]

trich.bagtype <- trich.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
trichophaea.bagtype <- lmerTest::lmer(Trichophaea ~ Bagtype + (1|Site.no/Cafeteria.ID), data = trich.bagtype)
summary(trichophaea.bagtype)
hist(resid(trichophaea.bagtype))
p.value.bagtype[13] <- coef(summary(trichophaea.bagtype))[2,5]

trich.sand <- trich.data %>% filter(Substrate %in% c("2_A", "1_S"))
trichophaea.sand <- lmerTest::lmer(Trichophaea ~ Substrate + (1|Site.no/Cafeteria.ID), data = trich.sand)
summary(trichophaea.sand)
hist(resid(trichophaea.sand))
p.value.sand[12] <- coef(summary(trichophaea.sand))[2,5]

trich.soil <- trich.data %>% filter(Substrate %in% c("RH", "F", "P", "M"))
trichophaea.soil <- lmerTest::lmer(Trichophaea ~ Substrate + (1|Site.no/Cafeteria.ID), data = trich.soil)
anova(trichophaea.soil)
summary(trichophaea.soil)
hist(resid(trichophaea.soil))
p.value.soil[13] <- anova(trichophaea.soil)$`Pr(>F)`
trich <- emmeans(trichophaea.soil, list(pairwise ~ Substrate), adjust = "tukey")

# ___Paxillus ####
pax.data <- data.sq %>% filter(Cafeteria.ID %in% c("38","5","6","28","7","25","21","36","40","4","39",
                                                "2","14","19","20","23","30","41","44","45","50"))

paxillus.bag <- lmerTest::lmer(Paxillus ~ Bag_root + (1|Site.no/Cafeteria.ID), data = pax.data)
summary(paxillus.bag)
hist(resid(paxillus.bag))
p.value.bag[14] <- coef(summary(paxillus.bag))[2,5]

pax.bagtype <- pax.data %>% filter(Bagtype %in% c("Sandbag", "Soilbag"))
paxillus.bagtype <- lmerTest::lmer(Paxillus ~ Bagtype + (1|Site.no/Cafeteria.ID), data = pax.bagtype)
summary(paxillus.bagtype)
hist(resid(paxillus.bagtype))
p.value.bagtype[14] <- coef(summary(paxillus.bagtype))[2,5]

pax.sand <- pax.data %>% filter(Substrate %in% c("2_A", "1_S"))
paxillus.sand <- lmerTest::lmer(Paxillus ~ Substrate + (1|Site.no/Cafeteria.ID), data = pax.sand)
summary(paxillus.sand)
hist(resid(paxillus.sand))
p.value.sand[13] <- coef(summary(paxillus.sand))[2,5]

pax.soil <- pax.data %>% filter(Substrate %in% c("RH", "F", "P", "M"))
paxillus.soil <- lmerTest::lmer(Paxillus ~ Substrate + (1|Site.no/Cafeteria.ID), data = pax.soil)
anova(paxillus.soil)
summary(paxillus.soil)
hist(resid(paxillus.soil))
p.value.soil[14] <- anova(paxillus.soil)$`Pr(>F)`
pax <- emmeans(paxillus.soil, list(pairwise ~ Substrate), adjust = "tukey")
#________________________________________________________________________________________________________________####
# P-value corrections ####
# Make a df with genera, un-adjusted p-values and BH-adjusted P-values
# Indicate with "*" if model is statistically significant (P<0.05) also
# after P-value correction

genera <- c("Amphinema", "Tomentella_Thelephora", "Cenococcum", "Russula", "Piloderma", "Tylospora",
            "Cortinarius", "Hyaloscypha", "Pseudotomentella", "Lactarius", "Amanita",
            "Tomentellopsis", "Trichophaea", "Paxillus")
genera2 <- c("Amphinema", "Tomentella_Thelephora", "Cenococcum", "Russula", "Piloderma", "Tylospora",
             "Cortinarius", "Hyaloscypha", "Pseudotomentella", "Amanita",
             "Tomentellopsis", "Trichophaea", "Paxillus") # no Lactarius in sandbags, removed from string

# _Bags vs. roots ####

bags.roots <- data.frame(Genus = genera, p.value = p.value.bag,
                         p.value.BH = p.adjust(p.value.bag, method="BH"))
bags.roots$sign <- ifelse(bags.roots$p.value.BH<=0.05, "*", "")

# _Soil bags vs. sand bags ####
soil.sand <- data.frame(Genus = genera, p.value = p.value.bagtype,
                         p.value.BH = p.adjust(p.value.bagtype, method="BH"))
soil.sand$sign <- ifelse(soil.sand$p.value.BH<=0.05, "*", "")

# _Sand vs. Apatite ####
sand.apatite <- data.frame(Genus = genera2, p.value = p.value.sand,
                           p.value.BH = p.adjust(p.value.sand, method="BH"))
sand.apatite$sign <- ifelse(sand.apatite$p.value.BH<=0.05, "*", "")

# _Different soil substrates ####
soil.substrates <- data.frame(Genus=genera, p.value=p.value.soil,
                              p.value.BH=p.adjust(p.value.soil, method="BH"))
soil.substrates$sign <- ifelse(soil.substrates$p.value.BH<=0.05, "*", "")

#_______________________________________________________________________________####

#Produce Figure 1 ####
# In Figure 1, the "relative preference" of ECM genera towards colonising bags relative to 
# the abundance on roots (x-axis) and to colonise soil-bags relative to sand bags (y-axis)

# "Relative preference" is calculated as the abundance ratio of bags/roots or soil bags/sand bags. 
# Use un-transformed relative abundances 

# Size of points reflect adjusted relative abundance of each genus 
# (mean abund on roots + mean abundance in bags)/2

#_Calculate mean abundances####
#colnames(data) # Fungi in 8:21 

# Bags relative to roots ####
# Filter separate datasets for bags and roots
bag.dat <- data %>% filter(Bag_root == "Bag")
root.dat <- data %>% filter(Bag_root == "Root")
sand.dat <- data %>% filter(Bagtype == "Sandbag")
soil.dat <- data %>% filter(Bagtype == "Soilbag")

# Gather data into long dataset to be able to calculate mean abundance
bags_long <- gather(bag.dat, "Genus", "Rel_abund", 8:21)
root_long <- gather(root.dat, "Genus", "Rel_abund", 8:21)
sand_long <- gather(sand.dat, "Genus", "Rel_abund", 8:21)
soil_long <- gather(soil.dat, "Genus", "Rel_abund", 8:21)

# Calculate mean rel abundance on roots, per genus

root.mean <- root_long %>%
  group_by(Genus) %>%
  summarise( 
    n.roots=n(),
    mean.roots=mean(Rel_abund, na.rm = TRUE),
    sd.roots=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.roots=sd.roots/sqrt(n.roots))

bag.mean <- bags_long %>%
  group_by(Genus) %>%
  summarise( 
    n.bags=n(),
    mean.bags=mean(Rel_abund, na.rm = TRUE),
    sd.bags=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.bags=sd.bags/sqrt(n.bags))

sand.mean <- sand_long %>%
  group_by(Genus) %>%
  summarise( 
    n.sandbags=n(),
    mean.sandbags=mean(Rel_abund, na.rm = TRUE),
    sd.sandbags=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.sandbags=sd.sandbags/sqrt(n.sandbags))

soil.mean <- soil_long %>%
  group_by(Genus) %>%
  summarise( 
    n.soilbags=n(),
    mean.soilbags=mean(Rel_abund, na.rm = TRUE),
    sd.soilbags=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.soilbags=sd.soilbags/sqrt(n.soilbags))

means <- left_join(root.mean, bag.mean, by="Genus")
means <- left_join(means, sand.mean, by="Genus")
means <- left_join(means, soil.mean, by="Genus")

means$adjusted.abund <- (means$mean.roots+means$mean.bags)/2
means$bagpref <- means$mean.bags/means$mean.roots
means$bagtypepref <- means$mean.soilbags/means$mean.sandbags

# Add info about exploration types
#Amphinema=MDF
#Tomentella=MDS
#Cenococcum=SD
#Russula=C
#Piloderma=MNF
#Tylospora=SD
#Cortinarius=MDF
#Hyalsocypha=C
#Pseudotomentella=MDS
#Lactarius=C
#Amanita=MDS
#Tomentellopsis=MDS
#Trichophaea=SD
#Paxillus=LD

expl_types <- c("MDF","MDS","SD","C","MDF","SD","MDF",
                    "C", "MDS","C","MDS","MDS","SD","LD")
exploration <- data.frame("Genus"=genera, "exploration_type"=expl_types)

rel_preference <- left_join(means, exploration, by = "Genus")

#_Drawing the figure####

Figure_1 <- ggplot(rel_preference, aes(x = bagpref, y = bagtypepref, color = exploration_type))+
  theme_classic()+
  geom_point(aes(size=adjusted.abund))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  scale_size_continuous(range = c(1, 14))+
  geom_hline(yintercept = 1, linetype="dotted", color = "black", size=0.6)+
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", size =0.6)+
  geom_text(label = rel_preference$Genus, nudge_x = 0.1, nudge_y = -0.1, check_overlap = F, color = "black")+
  scale_color_manual(name = "Exploration type",
                     values=c("#823038", "#D65721", "#F59C28", "#33658A", "#2F4858"),
                     breaks=c("C", "SD", "MDS", "MDF", "LD"),
                     labels= c("C", "SD", "MDS", "MDF", "LD"))

# The figure was edited in Adobe Illustrator to add symbols where lme-models 
# were significant + improve the general appearance. 

#_______________________________________________________________________________________####
# Produce Figure 2 ####

#In Figure 2, the relative preference for different soil substrates is plotted for all genera
# Relative preference is caluclated as: rel abundance in specific substrate/mean abundance in all substrates
#--> e.g. Mull soil preference = M/mean(M+P+RH+F)
#  F = organic topsoil, low N
#  M = mull soil
#  RH = organic topsoil, high N
#  P = organic topsoil, P fertilised

M.dat <- data %>% filter(Substrate == "M")
F.dat <- data %>% filter(Substrate== "F")
RH.dat <- data %>% filter(Substrate == "RH")
P.dat <- data %>% filter(Substrate == "P")

# Gather data into long dataset to be able to calculate mean abundance
M_long <- gather(M.dat, "Genus", "Rel_abund", 8:21)
F_long <- gather(F.dat, "Genus", "Rel_abund", 8:21)
RH_long <- gather(RH.dat, "Genus", "Rel_abund", 8:21)
P_long <- gather(P.dat, "Genus", "Rel_abund", 8:21)

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

substrate.means <- left_join(M.mean, F.mean, by="Genus")
substrate.means <- left_join(substrate.means, RH.mean, by="Genus")
substrate.means <- left_join(substrate.means, P.mean, by="Genus")

substrate.means$a_M.rel <- substrate.means$mean.M/((substrate.means$mean.M+substrate.means$mean.F+
                                                       substrate.means$mean.RH+substrate.means$mean.P)/4)
substrate.means$d_F.rel <- substrate.means$mean.F/((substrate.means$mean.M+substrate.means$mean.F+
                                                       substrate.means$mean.RH+substrate.means$mean.P)/4)
substrate.means$c_RH.rel <- substrate.means$mean.RH/((substrate.means$mean.M+substrate.means$mean.F+
                                                       substrate.means$mean.RH+substrate.means$mean.P)/4)
substrate.means$b_P.rel <- substrate.means$mean.P/((substrate.means$mean.M+substrate.means$mean.F+
                                                       substrate.means$mean.RH+substrate.means$mean.P)/4)



# In the figure, genera are ordered (from top to bottom) by their relative preference for
# bags over roots (X-axis in Figure 1) -->

# Add info about exploration types
#Amphinema=k
#Tomentella=g
#Cenococcum=d
#Russula=e
#Piloderma=c
#Tylospora=j
#Cortinarius=b
#Hyalsocypha=a
#Pseudotomentella=h
#Lactarius=f
#Amanita=l
#Tomentellopsis=m
#Trichophaea=i
#Paxillus=n

x.order <- c("k","g","d","e","c","j","b","a","h","f","l","m","i","n")
plot.order <- data.frame("Genus"=genera, "Order"=x.order)
substrate.means <- left_join(substrate.means, plot.order, by="Genus")

substrate.means$genus.order <- paste(substrate.means$Order, "_", substrate.means$Genus)

subs <- subset(substrate.means, select=c(genus.order, a_M.rel, d_F.rel, c_RH.rel, b_P.rel))
subs_long <- subs %>% gather(Substrate,Rel_abund, 2:5)

figure2 <- ggplot(subs_long, aes(x = genus.order, y = Rel_abund, color = Substrate))+
  theme_classic()+
  geom_point(position = position_dodge(0.9), size = 3)+
  coord_flip()+
  geom_hline(yintercept = 1, linetype="dotted", 
             color = "black", size=0.4)+
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "bottom",
    axis.text = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank())+
  scale_color_manual(name = "Substrate",
                     values=c("#002642", "#840032", "#E59500", "#A6A6A6"),
                     breaks=c("a_M.rel", "b_P.rel","c_RH.rel", "d_F.rel"),
                     labels = c("Mull soil", "Organic topsoil - P fertilised","Organic topsoil - high N deposition", "Organic topsoil - low N deposition"))


# The figure was edited in Adobe Illustrator to add symbols where lme-models 
# were significant + improve the general appearance. 