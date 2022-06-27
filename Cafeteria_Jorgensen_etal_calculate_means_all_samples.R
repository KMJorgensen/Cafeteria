# Calculate mean values in different substrates based on all samples


# Plotting and data handling
library(dplyr) # ver 1.0.3
library(tidyr) # ver 1.1.2
library(tibble)

data <- read.csv2("Rel_abundances.csv", header = T)
data$Site.no <- as.character(data$Site.no) # Site ID
data$Cafeteria.ID <- as.character(data$Cafeteria.ID) # Cafeteria ID
data$Caf <- as.character(data$Caf) # Cafeteria replicate, within site
data$Substrate <- as.factor(data$Substrate) # Either root or one of six bag substrates
data$Bagtype <- as.factor(data$Bagtype) # Either root, sandbag (S, A) or soilbag (F, M, RH, P)
data$Bag_root <- as.factor(data$Bag_root) # Root or ingrowth bag


colnames(data) # fungal genera in columns 8:19


#_Calculate mean abundances####
#colnames(data) # Fungi in 8:21 

# Bags relative to roots ####
# Filter separate datasets for bags and roots
bag.dat <- data %>% filter(Bag_root == "Bag")
root.dat <- data %>% filter(Bag_root == "Root")
sand.dat <- data %>% filter(Bagtype == "Sandbag")
soil.dat <- data %>% filter(Bagtype == "Soilbag")
M.dat <- data %>% filter(Substrate == "M")
F.dat <- data %>% filter(Substrate== "F")
RH.dat <- data %>% filter(Substrate == "RH")
P.dat <- data %>% filter(Substrate == "P")
S.dat <- data %>% filter(Substrate == "1_S")
A.dat <- data %>% filter(Substrate == "2_A")

# Gather data into long dataset to be able to calculate mean abundance
bags_long <- gather(bag.dat, "Genus", "Rel_abund", 8:19)
root_long <- gather(root.dat, "Genus", "Rel_abund", 8:19)
sand_long <- gather(sand.dat, "Genus", "Rel_abund", 8:19)
soil_long <- gather(soil.dat, "Genus", "Rel_abund", 8:19)
M_long <- gather(M.dat, "Genus", "Rel_abund", 8:19)
F_long <- gather(F.dat, "Genus", "Rel_abund", 8:19)
RH_long <- gather(RH.dat, "Genus", "Rel_abund", 8:19)
P_long <- gather(P.dat, "Genus", "Rel_abund", 8:19)
S_long <- gather(S.dat, "Genus", "Rel_abund", 8:19)
A_long <- gather(A.dat, "Genus", "Rel_abund", 8:19)

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

S.mean <- S_long %>%
  group_by(Genus) %>%
  summarise( 
    n.S=n(),
    mean.S=mean(Rel_abund, na.rm = TRUE),
    sd.S=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.S=sd.S/sqrt(n.S))

A.mean <- A_long %>%
  group_by(Genus) %>%
  summarise( 
    n.A=n(),
    mean.A=mean(Rel_abund, na.rm = TRUE),
    sd.A=sd(Rel_abund, na.rm = TRUE)
  ) %>%
  mutate(se.A=sd.A/sqrt(n.A))

root.mean <- root.mean[,c(-2,-4)]
bag.mean <- bag.mean[,c(-2,-4)]
soil.mean <- soil.mean[,c(-2,-4)]
sand.mean <- sand.mean[,c(-2,-4)]
F.mean <- F.mean[,c(-2,-4)]
RH.mean <- RH.mean[,c(-2,-4)]
P.mean <- P.mean[,c(-2,-4)]
M.mean <- M.mean[,c(-2,-4)]
S.mean <- S.mean[,c(-2,-4)]
A.mean <- A.mean[,c(-2,-4)]

means <- left_join(root.mean, bag.mean, by = "Genus")
means <- means %>% add_column(Adjusted.mean=NA, .after=1)
means$Adjusted.mean <- (means$mean.roots+means$mean.bags)/2
means <- left_join(means, soil.mean, by="Genus")
means <- left_join(means, sand.mean, by="Genus")
means <- left_join(means, F.mean, by="Genus")
means <- left_join(means, RH.mean, by="Genus")
means <- left_join(means, P.mean, by="Genus")
means <- left_join(means, M.mean, M.mean, by="Genus")
means <- left_join(means, S.mean, by="Genus")
means <- left_join(means, A.mean, by="Genus")

write.table(means, "means_SE.txt", sep = ",", quote = FALSE, row.names = F)

