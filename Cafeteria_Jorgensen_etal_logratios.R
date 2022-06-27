# R-script to reproduce analyses for:
#"A critical test of ectomycorrhizal exploration types as predictors of mycelial foraging"
# Jörgensen, K., Clemmensen, KE., Wallander, H., Lindahl, BD. 

# Contact about R script: karolina.jorgensen@slu.se                          

# This has been tested under:
# R version 4.0.1 (2022-06-27)


# Plotting and data handling
library(dplyr) # ver 1.0.3
library(tidyr) # ver 1.1.2


# ABUT THE SCRIPT ####
# This script is used to calculate log ratios of abundance of ectomycorrhizal genera
# Ratios are calculated as:
# Roots vs. ingrowth bags = log(mean(all bags)+µ/roots)
# Soil bags vs. sand bags = log((mean(soil)+µ)/(mean(sand)+µ))

setwd("~/docs/Dokument/Projekt/VR-projekt/Cafeteriaexperiment/Submission_prep/New Phytologist - revision/For GitHub")

# Load dataset with relative abundance (of ECM community) of 
# 12 genera that are present on roots in at least 10 cafeterias.
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

colnames(data) # fungal genera in columns 8:19

# Note about µ ####
# Average sequencing depth is 1881 counts per sample, there are 6 substrate types in each cafeteria.
# Create a constant (µ) that will be added to abundances in bags to avoid 0 division
µ <- 1/(6*1881)

# Calculate log-ratios ####

# _Amanita ####
# Make new dataset with only one genus
aman.data <- data[,c(1:7, 8)]
# Filter out cafeterias where genus is present on roots 
aman.root.pres <- data %>% filter(Bag_root == "Root" & Amanita >0)
# Make a vector of cafeteria numbers where genus is present
aman.root.pres.vector <- aman.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
aman.data <- aman.data %>% filter(Cafeteria.ID %in% c(aman.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
aman.roots <- aman.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
aman.bag.dat <- aman.data %>% filter(Bag_root == "Bag")
aman.bag.long <- gather(aman.bag.dat, "Genus", "Rel_abund", 8)

aman.bag.mean <- aman.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

aman.bag.root <- left_join(aman.roots, aman.bag.mean, by = "Cafeteria.ID")
aman.bag.root$log.ratio <- log((aman.bag.root$mean.bags + µ)/aman.bag.root$Amanita)

# ___Sandbag-soilbag comparison ####
aman.soil <- aman.data %>% filter(Bagtype == "Soilbag")
aman.sand <- aman.data %>% filter(Bagtype == "Sandbag")

aman.soil.long <- gather(aman.soil, "Genus", "Rel_abund", 8)
aman.sand.long <- gather(aman.sand, "Genus", "Rel_abund", 8)

aman.soil.mean <- aman.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

aman.sand.mean <- aman.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

aman.bagtypes <- left_join(aman.soil.mean, aman.sand.mean, by = "Cafeteria.ID")
aman.bagtypes$log.ratio <- log((aman.bagtypes$mean.soilbags+µ)/(aman.bagtypes$mean.sandbags+µ))


#___Gather into one DF ####
amanita.logratios <- data_frame(matrix(NA, nrow=nrow(aman.bag.root)))
amanita.logratios$Site <- aman.bag.root$Site.no
amanita.logratios$Cafeteria.ID <- aman.bag.root$Cafeteria.ID
amanita.logratios$Genus <- aman.bag.root$Genus
amanita.logratios$Root_bag <- aman.bag.root$log.ratio
amanita.logratios$Bagtype <- aman.bagtypes$log.ratio
amanita.logratios <- amanita.logratios[,-1]

amanita.logratios.long <- gather(amanita.logratios, "Comparison", "Log_ratio", 4:5)

# _Amphinema ####
# Make new dataset with only one genus
amph.data <- data[,c(1:7, 9)]
# Filter out cafeterias where genus is present on roots 
amph.root.pres <- data %>% filter(Bag_root == "Root" & Amphinema >0)
# Make a vector of cafeteria numbers where genus is present
amph.root.pres.vector <- amph.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
amph.data <- amph.data %>% filter(Cafeteria.ID %in% c(amph.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
amph.roots <- amph.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
amph.bag.dat <- amph.data %>% filter(Bag_root == "Bag")
amph.bag.long <- gather(amph.bag.dat, "Genus", "Rel_abund", 8)

amph.bag.mean <- amph.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

amph.bag.root <- left_join(amph.roots, amph.bag.mean, by = "Cafeteria.ID")
amph.bag.root$log.ratio <- log((amph.bag.root$mean.bags + µ)/amph.bag.root$Amphinema)

# ___Sandbag-soilbag comparison ####
amph.soil <- amph.data %>% filter(Bagtype == "Soilbag")
amph.sand <- amph.data %>% filter(Bagtype == "Sandbag")

amph.soil.long <- gather(amph.soil, "Genus", "Rel_abund", 8)
amph.sand.long <- gather(amph.sand, "Genus", "Rel_abund", 8)

amph.soil.mean <- amph.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

amph.sand.mean <- amph.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

amph.bagtypes <- left_join(amph.soil.mean, amph.sand.mean, by = "Cafeteria.ID")
amph.bagtypes$log.ratio <- log((amph.bagtypes$mean.soilbags+µ)/(amph.bagtypes$mean.sandbags+µ))


#___Gather into one DF ####
amphinema.logratios <- data_frame(matrix(NA, nrow=nrow(amph.bag.root)))
amphinema.logratios$Site <- amph.bag.root$Site.no
amphinema.logratios$Cafeteria.ID <- amph.bag.root$Cafeteria.ID
amphinema.logratios$Genus <- amph.bag.root$Genus
amphinema.logratios$Root_bag <- amph.bag.root$log.ratio
amphinema.logratios$Bagtype <- amph.bagtypes$log.ratio
amphinema.logratios <- amphinema.logratios[,-1]

amphinema.logratios.long <- gather(amphinema.logratios, "Comparison", "Log_ratio", 4:5)

# _Cenococcum ####
# Make new dataset with only one genus
ceno.data <- data[,c(1:7, 10)]
# Filter out cafeterias where genus is present on roots 
ceno.root.pres <- data %>% filter(Bag_root == "Root" & Cenococcum >0)
# Make a vector of cafeteria numbers where genus is present
ceno.root.pres.vector <- ceno.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
ceno.data <- ceno.data %>% filter(Cafeteria.ID %in% c(ceno.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
ceno.roots <- ceno.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
ceno.bag.dat <- ceno.data %>% filter(Bag_root == "Bag")
ceno.bag.long <- gather(ceno.bag.dat, "Genus", "Rel_abund", 8)

ceno.bag.mean <- ceno.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

ceno.bag.root <- left_join(ceno.roots, ceno.bag.mean, by = "Cafeteria.ID")
ceno.bag.root$log.ratio <- log((ceno.bag.root$mean.bags + µ)/ceno.bag.root$Cenococcum)

# ___Sandbag-soilbag comparison ####
ceno.soil <- ceno.data %>% filter(Bagtype == "Soilbag")
ceno.sand <- ceno.data %>% filter(Bagtype == "Sandbag")

ceno.soil.long <- gather(ceno.soil, "Genus", "Rel_abund", 8)
ceno.sand.long <- gather(ceno.sand, "Genus", "Rel_abund", 8)

ceno.soil.mean <- ceno.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

ceno.sand.mean <- ceno.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

ceno.bagtypes <- left_join(ceno.soil.mean, ceno.sand.mean, by = "Cafeteria.ID")
ceno.bagtypes$log.ratio <- log((ceno.bagtypes$mean.soilbags+µ)/(ceno.bagtypes$mean.sandbags+µ))


#___Gather into one DF ####
cenococcum.logratios <- data_frame(matrix(NA, nrow=nrow(ceno.bag.root)))
cenococcum.logratios$Site <- ceno.bag.root$Site.no
cenococcum.logratios$Cafeteria.ID <- ceno.bag.root$Cafeteria.ID
cenococcum.logratios$Genus <- ceno.bag.root$Genus
cenococcum.logratios$Root_bag <- ceno.bag.root$log.ratio
cenococcum.logratios$Bagtype <- ceno.bagtypes$log.ratio
cenococcum.logratios <- cenococcum.logratios[,-1]

cenococcum.logratios.long <- gather(cenococcum.logratios, "Comparison", "Log_ratio", 4:5)

# _Cortinarius ####
# Make new dataset with only one genus
cort.data <- data[,c(1:7, 11)]
# Filter out cafeterias where genus is present on roots 
cort.root.pres <- data %>% filter(Bag_root == "Root" & Cortinarius >0)
# Make a vector of cafeteria numbers where genus is present
cort.root.pres.vector <- cort.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
cort.data <- cort.data %>% filter(Cafeteria.ID %in% c(cort.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
cort.roots <- cort.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
cort.bag.dat <- cort.data %>% filter(Bag_root == "Bag")
cort.bag.long <- gather(cort.bag.dat, "Genus", "Rel_abund", 8)

cort.bag.mean <- cort.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

cort.bag.root <- left_join(cort.roots, cort.bag.mean, by = "Cafeteria.ID")
cort.bag.root$log.ratio <- log((cort.bag.root$mean.bags + µ)/cort.bag.root$Cortinarius)

# ___Sandbag-soilbag comparison ####
cort.soil <- cort.data %>% filter(Bagtype == "Soilbag")
cort.sand <- cort.data %>% filter(Bagtype == "Sandbag")

cort.soil.long <- gather(cort.soil, "Genus", "Rel_abund", 8)
cort.sand.long <- gather(cort.sand, "Genus", "Rel_abund", 8)

cort.soil.mean <- cort.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

cort.sand.mean <- cort.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

cort.bagtypes <- left_join(cort.soil.mean, cort.sand.mean, by = "Cafeteria.ID")
cort.bagtypes$log.ratio <- log((cort.bagtypes$mean.soilbags+µ)/(cort.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
cortinarius.logratios <- data_frame(matrix(NA, nrow=nrow(cort.bag.root)))
cortinarius.logratios$Site <- cort.bag.root$Site.no
cortinarius.logratios$Cafeteria.ID <- cort.bag.root$Cafeteria.ID
cortinarius.logratios$Genus <- cort.bag.root$Genus
cortinarius.logratios$Root_bag <- cort.bag.root$log.ratio
cortinarius.logratios$Bagtype <- cort.bagtypes$log.ratio
cortinarius.logratios <- cortinarius.logratios[,-1]

cortinarius.logratios.long <- gather(cortinarius.logratios, "Comparison", "Log_ratio", 4:5)

# _Hyaloscypha ####
# Make new dataset with only one genus
hyal.data <- data[,c(1:7, 12)]
# Filter out cafeterias where genus is present on roots 
hyal.root.pres <- data %>% filter(Bag_root == "Root" & Hyaloscypha >0)
# Make a vector of cafeteria numbers where genus is present
hyal.root.pres.vector <- hyal.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
hyal.data <- hyal.data %>% filter(Cafeteria.ID %in% c(hyal.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
hyal.roots <- hyal.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
hyal.bag.dat <- hyal.data %>% filter(Bag_root == "Bag")
hyal.bag.long <- gather(hyal.bag.dat, "Genus", "Rel_abund", 8)

hyal.bag.mean <- hyal.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

hyal.bag.root <- left_join(hyal.roots, hyal.bag.mean, by = "Cafeteria.ID")
hyal.bag.root$log.ratio <- log((hyal.bag.root$mean.bags + µ)/hyal.bag.root$Hyaloscypha)

# ___Sandbag-soilbag comparison ####
hyal.soil <- hyal.data %>% filter(Bagtype == "Soilbag")
hyal.sand <- hyal.data %>% filter(Bagtype == "Sandbag")

hyal.soil.long <- gather(hyal.soil, "Genus", "Rel_abund", 8)
hyal.sand.long <- gather(hyal.sand, "Genus", "Rel_abund", 8)

hyal.soil.mean <- hyal.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

hyal.sand.mean <- hyal.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

hyal.bagtypes <- left_join(hyal.soil.mean, hyal.sand.mean, by = "Cafeteria.ID")
hyal.bagtypes$log.ratio <- log((hyal.bagtypes$mean.soilbags+µ)/(hyal.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
hyaloscypha.logratios <- data_frame(matrix(NA, nrow=nrow(hyal.bag.root)))
hyaloscypha.logratios$Site <- hyal.bag.root$Site.no
hyaloscypha.logratios$Cafeteria.ID <- hyal.bag.root$Cafeteria.ID
hyaloscypha.logratios$Genus <- hyal.bag.root$Genus
hyaloscypha.logratios$Root_bag <- hyal.bag.root$log.ratio
hyaloscypha.logratios$Bagtype <- hyal.bagtypes$log.ratio
hyaloscypha.logratios <- hyaloscypha.logratios[,-1]

hyaloscypha.logratios.long <- gather(hyaloscypha.logratios, "Comparison", "Log_ratio", 4:5)

# _Hygrophorus ####
# Make new dataset with only one genus
hyg.data <- data[,c(1:7, 13)]
# Filter out cafeterias where genus is present on roots 
hyg.root.pres <- data %>% filter(Bag_root == "Root" & Hygrophorus >0)
# Make a vector of cafeteria numbers where genus is present
hyg.root.pres.vector <- hyg.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
hyg.data <- hyg.data %>% filter(Cafeteria.ID %in% c(hyg.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
hyg.roots <- hyg.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
hyg.bag.dat <- hyg.data %>% filter(Bag_root == "Bag")
hyg.bag.long <- gather(hyg.bag.dat, "Genus", "Rel_abund", 8)

hyg.bag.mean <- hyg.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

hyg.bag.root <- left_join(hyg.roots, hyg.bag.mean, by = "Cafeteria.ID")
hyg.bag.root$log.ratio <- log((hyg.bag.root$mean.bags + µ)/hyg.bag.root$Hygrophorus)

# ___Sandbag-soilbag comparison ####
hyg.soil <- hyg.data %>% filter(Bagtype == "Soilbag")
hyg.sand <- hyg.data %>% filter(Bagtype == "Sandbag")

hyg.soil.long <- gather(hyg.soil, "Genus", "Rel_abund", 8)
hyg.sand.long <- gather(hyg.sand, "Genus", "Rel_abund", 8)

hyg.soil.mean <- hyg.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

hyg.sand.mean <- hyg.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

hyg.bagtypes <- left_join(hyg.soil.mean, hyg.sand.mean, by = "Cafeteria.ID")
hyg.bagtypes$log.ratio <- log((hyg.bagtypes$mean.soilbags+µ)/(hyg.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
hygrophorus.logratios <- data_frame(matrix(NA, nrow=nrow(hyg.bag.root)))
hygrophorus.logratios$Site <- hyg.bag.root$Site.no
hygrophorus.logratios$Cafeteria.ID <- hyg.bag.root$Cafeteria.ID
hygrophorus.logratios$Genus <- hyg.bag.root$Genus
hygrophorus.logratios$Root_bag <- hyg.bag.root$log.ratio
hygrophorus.logratios$Bagtype <- hyg.bagtypes$log.ratio
hygrophorus.logratios <- hygrophorus.logratios[,-1]

hygrophorus.logratios.long <- gather(hygrophorus.logratios, "Comparison", "Log_ratio", 4:5)

# _Lactarius ####
# Make new dataset with only one genus
lact.data <- data[,c(1:7, 14)]
# Filter out cafeterias where genus is present on roots 
lact.root.pres <- data %>% filter(Bag_root == "Root" & Lactarius >0)
# Make a vector of cafeteria numbers where genus is present
lact.root.pres.vector <- lact.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
lact.data <- lact.data %>% filter(Cafeteria.ID %in% c(lact.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
lact.roots <- lact.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
lact.bag.dat <- lact.data %>% filter(Bag_root == "Bag")
lact.bag.long <- gather(lact.bag.dat, "Genus", "Rel_abund", 8)

lact.bag.mean <- lact.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

lact.bag.root <- left_join(lact.roots, lact.bag.mean, by = "Cafeteria.ID")
lact.bag.root$log.ratio <- log((lact.bag.root$mean.bags + µ)/lact.bag.root$Lactarius)

# ___Sandbag-soilbag comparison ####
lact.soil <- lact.data %>% filter(Bagtype == "Soilbag")
lact.sand <- lact.data %>% filter(Bagtype == "Sandbag")

lact.soil.long <- gather(lact.soil, "Genus", "Rel_abund", 8)
lact.sand.long <- gather(lact.sand, "Genus", "Rel_abund", 8)

lact.soil.mean <- lact.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

lact.sand.mean <- lact.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

lact.bagtypes <- left_join(lact.soil.mean, lact.sand.mean, by = "Cafeteria.ID")
lact.bagtypes$log.ratio <- log((lact.bagtypes$mean.soilbags+µ)/(lact.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
lactarius.logratios <- data_frame(matrix(NA, nrow=nrow(lact.bag.root)))
lactarius.logratios$Site <- lact.bag.root$Site.no
lactarius.logratios$Cafeteria.ID <- lact.bag.root$Cafeteria.ID
lactarius.logratios$Genus <- lact.bag.root$Genus
lactarius.logratios$Root_bag <- lact.bag.root$log.ratio
lactarius.logratios$Bagtype <- lact.bagtypes$log.ratio
lactarius.logratios <- lactarius.logratios[,-1]

lactarius.logratios.long <- gather(lactarius.logratios, "Comparison", "Log_ratio", 4:5)

# _Piloderma ####
# Make new dataset with only one genus
pilo.data <- data[,c(1:7, 15)]
# Filter out cafeterias where genus is present on roots 
pilo.root.pres <- data %>% filter(Bag_root == "Root" & Piloderma >0)
# Make a vector of cafeteria numbers where genus is present
pilo.root.pres.vector <- pilo.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
pilo.data <- pilo.data %>% filter(Cafeteria.ID %in% c(pilo.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
pilo.roots <- pilo.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
pilo.bag.dat <- pilo.data %>% filter(Bag_root == "Bag")
pilo.bag.long <- gather(pilo.bag.dat, "Genus", "Rel_abund", 8)

pilo.bag.mean <- pilo.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

pilo.bag.root <- left_join(pilo.roots, pilo.bag.mean, by = "Cafeteria.ID")
pilo.bag.root$log.ratio <- log((pilo.bag.root$mean.bags + µ)/pilo.bag.root$Piloderma)

# ___Sandbag-soilbag comparison ####
pilo.soil <- pilo.data %>% filter(Bagtype == "Soilbag")
pilo.sand <- pilo.data %>% filter(Bagtype == "Sandbag")

pilo.soil.long <- gather(pilo.soil, "Genus", "Rel_abund", 8)
pilo.sand.long <- gather(pilo.sand, "Genus", "Rel_abund", 8)

pilo.soil.mean <- pilo.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

pilo.sand.mean <- pilo.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

pilo.bagtypes <- left_join(pilo.soil.mean, pilo.sand.mean, by = "Cafeteria.ID")
pilo.bagtypes$log.ratio <- log((pilo.bagtypes$mean.soilbags+µ)/(pilo.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
piloderma.logratios <- data_frame(matrix(NA, nrow=nrow(pilo.bag.root)))
piloderma.logratios$Site <- pilo.bag.root$Site.no
piloderma.logratios$Cafeteria.ID <- pilo.bag.root$Cafeteria.ID
piloderma.logratios$Genus <- pilo.bag.root$Genus
piloderma.logratios$Root_bag <- pilo.bag.root$log.ratio
piloderma.logratios$Bagtype <- pilo.bagtypes$log.ratio
piloderma.logratios <- piloderma.logratios[,-1]

piloderma.logratios.long <- gather(piloderma.logratios, "Comparison", "Log_ratio", 4:5)

# _Pseudotomentella ####
# Make new dataset with only one genus
pseu.data <- data[,c(1:7, 16)]
# Filter out cafeterias where genus is present on roots 
pseu.root.pres <- data %>% filter(Bag_root == "Root" & Pseudotomentella >0)
# Make a vector of cafeteria numbers where genus is present
pseu.root.pres.vector <- pseu.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
pseu.data <- pseu.data %>% filter(Cafeteria.ID %in% c(pseu.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
pseu.roots <- pseu.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
pseu.bag.dat <- pseu.data %>% filter(Bag_root == "Bag")
pseu.bag.long <- gather(pseu.bag.dat, "Genus", "Rel_abund", 8)

pseu.bag.mean <- pseu.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

pseu.bag.root <- left_join(pseu.roots, pseu.bag.mean, by = "Cafeteria.ID")
pseu.bag.root$log.ratio <- log((pseu.bag.root$mean.bags + µ)/pseu.bag.root$Pseudotomentella)

# ___Sandbag-soilbag comparison ####
pseu.soil <- pseu.data %>% filter(Bagtype == "Soilbag")
pseu.sand <- pseu.data %>% filter(Bagtype == "Sandbag")

pseu.soil.long <- gather(pseu.soil, "Genus", "Rel_abund", 8)
pseu.sand.long <- gather(pseu.sand, "Genus", "Rel_abund", 8)

pseu.soil.mean <- pseu.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

pseu.sand.mean <- pseu.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

pseu.bagtypes <- left_join(pseu.soil.mean, pseu.sand.mean, by = "Cafeteria.ID")
pseu.bagtypes$log.ratio <- log((pseu.bagtypes$mean.soilbags+µ)/(pseu.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
pseudotomentella.logratios <- data_frame(matrix(NA, nrow=nrow(pseu.bag.root)))
pseudotomentella.logratios$Site <- pseu.bag.root$Site.no
pseudotomentella.logratios$Cafeteria.ID <- pseu.bag.root$Cafeteria.ID
pseudotomentella.logratios$Genus <- pseu.bag.root$Genus
pseudotomentella.logratios$Root_bag <- pseu.bag.root$log.ratio
pseudotomentella.logratios$Bagtype <- pseu.bagtypes$log.ratio
pseudotomentella.logratios <- pseudotomentella.logratios[,-1]

pseudotomentella.logratios.long <- gather(pseudotomentella.logratios, "Comparison", "Log_ratio", 4:5)

# _Russula ####
# Make new dataset with only one genus
russ.data <- data[,c(1:7, 17)]
# Filter out cafeterias where genus is present on roots 
russ.root.pres <- data %>% filter(Bag_root == "Root" & Russula >0)
# Make a vector of cafeteria numbers where genus is present
russ.root.pres.vector <- russ.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
russ.data <- russ.data %>% filter(Cafeteria.ID %in% c(russ.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
russ.roots <- russ.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
russ.bag.dat <- russ.data %>% filter(Bag_root == "Bag")
russ.bag.long <- gather(russ.bag.dat, "Genus", "Rel_abund", 8)

russ.bag.mean <- russ.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

russ.bag.root <- left_join(russ.roots, russ.bag.mean, by = "Cafeteria.ID")
russ.bag.root$log.ratio <- log((russ.bag.root$mean.bags + µ)/russ.bag.root$Russula)

# ___Sandbag-soilbag comparison ####
russ.soil <- russ.data %>% filter(Bagtype == "Soilbag")
russ.sand <- russ.data %>% filter(Bagtype == "Sandbag")

russ.soil.long <- gather(russ.soil, "Genus", "Rel_abund", 8)
russ.sand.long <- gather(russ.sand, "Genus", "Rel_abund", 8)

russ.soil.mean <- russ.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

russ.sand.mean <- russ.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

russ.bagtypes <- left_join(russ.soil.mean, russ.sand.mean, by = "Cafeteria.ID")
russ.bagtypes$log.ratio <- log((russ.bagtypes$mean.soilbags+µ)/(russ.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
russula.logratios <- data_frame(matrix(NA, nrow=nrow(russ.bag.root)))
russula.logratios$Site <- russ.bag.root$Site.no
russula.logratios$Cafeteria.ID <- russ.bag.root$Cafeteria.ID
russula.logratios$Genus <- russ.bag.root$Genus
russula.logratios$Root_bag <- russ.bag.root$log.ratio
russula.logratios$Bagtype <- russ.bagtypes$log.ratio
russula.logratios <- russula.logratios[,-1]

russula.logratios.long <- gather(russula.logratios, "Comparison", "Log_ratio", 4:5)

# _Tomentella (incl. Thelephora) ####
# Make new dataset with only one genus
tomen.data <- data[,c(1:7, 18)]
# Filter out cafeterias where genus is present on roots 
tomen.root.pres <- data %>% filter(Bag_root == "Root" & Tomentella_Thelephora >0)
# Make a vector of cafeteria numbers where genus is present
tomen.root.pres.vector <- tomen.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
tomen.data <- tomen.data %>% filter(Cafeteria.ID %in% c(tomen.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
tomen.roots <- tomen.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
tomen.bag.dat <- tomen.data %>% filter(Bag_root == "Bag")
tomen.bag.long <- gather(tomen.bag.dat, "Genus", "Rel_abund", 8)

tomen.bag.mean <- tomen.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

tomen.bag.root <- left_join(tomen.roots, tomen.bag.mean, by = "Cafeteria.ID")
tomen.bag.root$log.ratio <- log((tomen.bag.root$mean.bags + µ)/tomen.bag.root$Tomentella_Thelephora)

# ___Sandbag-soilbag comparison ####
tomen.soil <- tomen.data %>% filter(Bagtype == "Soilbag")
tomen.sand <- tomen.data %>% filter(Bagtype == "Sandbag")

tomen.soil.long <- gather(tomen.soil, "Genus", "Rel_abund", 8)
tomen.sand.long <- gather(tomen.sand, "Genus", "Rel_abund", 8)

tomen.soil.mean <- tomen.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

tomen.sand.mean <- tomen.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

tomen.bagtypes <- left_join(tomen.soil.mean, tomen.sand.mean, by = "Cafeteria.ID")
tomen.bagtypes$log.ratio <- log((tomen.bagtypes$mean.soilbags+µ)/(tomen.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
tomentella.logratios <- data_frame(matrix(NA, nrow=nrow(tomen.bag.root)))
tomentella.logratios$Site <- tomen.bag.root$Site.no
tomentella.logratios$Cafeteria.ID <- tomen.bag.root$Cafeteria.ID
tomentella.logratios$Genus <- tomen.bag.root$Genus
tomentella.logratios$Root_bag <- tomen.bag.root$log.ratio
tomentella.logratios$Bagtype <- tomen.bagtypes$log.ratio
tomentella.logratios <- tomentella.logratios[,-1]

tomentella.logratios.long <- gather(tomentella.logratios, "Comparison", "Log_ratio", 4:5)

# _Tylospora ####
# Make new dataset with only one genus
tylo.data <- data[,c(1:7, 19)]
# Filter out cafeterias where genus is present on roots 
tylo.root.pres <- data %>% filter(Bag_root == "Root" & Tylospora >0)
# Make a vector of cafeteria numbers where genus is present
tylo.root.pres.vector <- tylo.root.pres$Cafeteria.ID
# Use vector to filter out cafeterias (all substrates) where the genus is present on roots
tylo.data <- tylo.data %>% filter(Cafeteria.ID %in% c(tylo.root.pres.vector))

# ___Root-bag comparison ####
#Subset roots only
tylo.roots <- tylo.data %>% filter(Bag_root == "Root")

# Calculate mean abundance in bags 
tylo.bag.dat <- tylo.data %>% filter(Bag_root == "Bag")
tylo.bag.long <- gather(tylo.bag.dat, "Genus", "Rel_abund", 8)

tylo.bag.mean <- tylo.bag.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.bags=mean(Rel_abund, na.rm = FALSE))

tylo.bag.root <- left_join(tylo.roots, tylo.bag.mean, by = "Cafeteria.ID")
tylo.bag.root$log.ratio <- log((tylo.bag.root$mean.bags + µ)/tylo.bag.root$Tylospora)

# ___Sandbag-soilbag comparison ####
tylo.soil <- tylo.data %>% filter(Bagtype == "Soilbag")
tylo.sand <- tylo.data %>% filter(Bagtype == "Sandbag")

tylo.soil.long <- gather(tylo.soil, "Genus", "Rel_abund", 8)
tylo.sand.long <- gather(tylo.sand, "Genus", "Rel_abund", 8)

tylo.soil.mean <- tylo.soil.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.soilbags=mean(Rel_abund, na.rm = FALSE))

tylo.sand.mean <- tylo.sand.long %>%
  group_by(Genus, Cafeteria.ID) %>%
  summarise( 
    mean.sandbags=mean(Rel_abund, na.rm = FALSE))

tylo.bagtypes <- left_join(tylo.soil.mean, tylo.sand.mean, by = "Cafeteria.ID")
tylo.bagtypes$log.ratio <- log((tylo.bagtypes$mean.soilbags+µ)/(tylo.bagtypes$mean.sandbags+µ))

#___Gather into one DF ####
tylospora.logratios <- data_frame(matrix(NA, nrow=nrow(tylo.bag.root)))
tylospora.logratios$Site <- tylo.bag.root$Site.no
tylospora.logratios$Cafeteria.ID <- tylo.bag.root$Cafeteria.ID
tylospora.logratios$Genus <- tylo.bag.root$Genus
tylospora.logratios$Root_bag <- tylo.bag.root$log.ratio
tylospora.logratios$Bagtype <- tylo.bagtypes$log.ratio
tylospora.logratios <- tylospora.logratios[,-1]

tylospora.logratios.long <- gather(tylospora.logratios, "Comparison", "Log_ratio", 4:5)

# DF log-ratios all genera ####
# Make one big DF with all genera 
# Make a .csv to use in a separate stats script

all.genera.logratios <- rbind(amanita.logratios.long, 
                              amphinema.logratios.long,
                              cenococcum.logratios.long,
                              cortinarius.logratios.long,
                              hyaloscypha.logratios.long,
                              hygrophorus.logratios.long,
                              lactarius.logratios.long,
                              piloderma.logratios.long,
                              pseudotomentella.logratios.long,
                              russula.logratios.long,
                              tomentella.logratios.long,
                              tylospora.logratios.long)

write.csv2(all.genera.logratios, "genera_logratios.csv", row.names = FALSE)


# Export .csv for soilbags and sandbags ####
# These files are needed to run the models for each genus in
# Test 4 in the "Statistics script.R" 

write.csv2(aman.soil, "amanita_soil.csv", row.names=FALSE)
write.csv2(amph.soil, "amphinema_soil.csv", row.names=FALSE)
write.csv2(ceno.soil, "cenococcum_soil.csv", row.names=FALSE)
write.csv2(cort.soil, "cortinarius_soil.csv", row.names=FALSE)
write.csv2(hyal.soil, "hyaloscypha_soil.csv", row.names=FALSE)
write.csv2(hyg.soil, "hygrophorus_soil.csv", row.names=FALSE)
write.csv2(lact.soil, "lactarius_soil.csv", row.names=FALSE)
write.csv2(pilo.soil, "piloderma_soil.csv", row.names=FALSE)
write.csv2(pseu.soil, "pseudotomentella_soil.csv", row.names=FALSE)
write.csv2(russ.soil, "russula_soil.csv", row.names=FALSE)
write.csv2(tomen.soil, "tomentella_soil.csv", row.names=FALSE)
write.csv2(tylo.soil, "tylospora_soil.csv", row.names=FALSE)


# Calculate adjusted abundances####
# Adjusted abundance is used for point sizes in Figure 1
# adjusted.abund=(mean(root)+mean(bags))/2

aman.adj <- inner_join(aman.roots[,c(3,8)], aman.bag.mean[,2:3], by = "Cafeteria.ID")
amph.adj <- inner_join(amph.roots[,c(3,8)], amph.bag.mean[,2:3], by = "Cafeteria.ID")
ceno.adj <- inner_join(ceno.roots[,c(3,8)], ceno.bag.mean[,2:3], by = "Cafeteria.ID")
cort.adj <- inner_join(cort.roots[,c(3,8)], cort.bag.mean[,2:3], by = "Cafeteria.ID")
hyal.adj <- inner_join(hyal.roots[,c(3,8)], hyal.bag.mean[,2:3], by = "Cafeteria.ID")
hyg.adj <- inner_join(hyg.roots[,c(3,8)], hyg.bag.mean[,2:3], by = "Cafeteria.ID")
lact.adj <- inner_join(lact.roots[,c(3,8)], lact.bag.mean[,2:3], by = "Cafeteria.ID")
pilo.adj <- inner_join(pilo.roots[,c(3,8)], pilo.bag.mean[,2:3], by = "Cafeteria.ID")
pseu.adj <- inner_join(pseu.roots[,c(3,8)], pseu.bag.mean[,2:3], by = "Cafeteria.ID")
russ.adj <- inner_join(russ.roots[,c(3,8)], russ.bag.mean[,2:3], by = "Cafeteria.ID")
tomen.adj <- inner_join(tomen.roots[,c(3,8)], tomen.bag.mean[,2:3], by = "Cafeteria.ID")
tylo.adj <- inner_join(tylo.roots[,c(3,8)], tylo.bag.mean[,2:3], by = "Cafeteria.ID")

aman.adj$adj.mean <- (aman.adj$Amanita+aman.adj$mean.bags)/2
amph.adj$adj.mean <- (amph.adj$Amphinema+amph.adj$mean.bags)/2
ceno.adj$adj.mean <- (ceno.adj$Cenococcum+ceno.adj$mean.bags)/2
cort.adj$adj.mean <- (cort.adj$Cortinarius+cort.adj$mean.bags)/2
hyal.adj$adj.mean <- (hyal.adj$Hyaloscypha+hyal.adj$mean.bags)/2
hyg.adj$adj.mean <- (hyg.adj$Hygrophorus+hyg.adj$mean.bags)/2
lact.adj$adj.mean <- (lact.adj$Lactarius+lact.adj$mean.bags)/2
pilo.adj$adj.mean <- (pilo.adj$Piloderma+pilo.adj$mean.bags)/2
pseu.adj$adj.mean <- (pseu.adj$Pseudotomentella+pseu.adj$mean.bags)/2
russ.adj$adj.mean <- (russ.adj$Russula+russ.adj$mean.bags)/2
tomen.adj$adj.mean <- (tomen.adj$Tomentella+tomen.adj$mean.bags)/2
tylo.adj$adj.mean <- (tylo.adj$Tylospora+tylo.adj$mean.bags)/2


adj.means <- numeric()
adj.means[1] <- mean(aman.adj$adj.mean)
adj.means[2] <- mean(amph.adj$adj.mean)
adj.means[3] <- mean(ceno.adj$adj.mean)
adj.means[4] <- mean(cort.adj$adj.mean)
adj.means[5] <- mean(hyal.adj$adj.mean)
adj.means[6] <- mean(hyg.adj$adj.mean)
adj.means[7] <- mean(lact.adj$adj.mean)
adj.means[8] <- mean(pilo.adj$adj.mean)
adj.means[9] <- mean(pseu.adj$adj.mean)
adj.means[10] <- mean(russ.adj$adj.mean)
adj.means[11] <- mean(tomen.adj$adj.mean)
adj.means[12] <- mean(tylo.adj$adj.mean)


genera <- c("Amanita", "Amphinema", "Cenococcum", "Cortinarius",
            "Hyaloscypha", "Hygrophorus", "Lactarius","Piloderma",
            "Pseudotomentella", "Russula", "Tomentella", "Tylospora")

adj.means.df <- tibble(matrix(NA, nrow=12))
adj.means.df$Genus <- genera
adj.means.df$Adj.mean <- adj.means
adj.means.df <- adj.means.df[,-1]

write.csv2(adj.means.df, "Adjusted.means.csv", row.names=TRUE)

