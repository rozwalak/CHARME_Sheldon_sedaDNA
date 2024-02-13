library(readxl)
library(tidyverse)
library(rioja)
library(vegan)
library(car)
library(ggpubr)

#Loading data (available in input folder in github repo or in Supplementary Table associated to Article)
sheldon_db <- read_excel("input/Sheldon_data.xlsx")

#Below we computed hierarchical clustering (CONISS) on different combination of sea-ice proxies
all_proxy <- sheldon_db %>%
  select("Age", "Diatoms_rel_abundance", "Dinoflagellata_rel_abundance", "Polarella_rel_to_Dinoflagellata", "Islandinium_rel_to_Dinoflagellata","IPSO25_ug_to_g_sed","HBI_III_ug_to_g_sed", "BiogenicSiO2", "TOC") %>%
  mutate(Age = as.integer(Age)) %>%
  column_to_rownames(var = "Age")

all_proxy.dist <- vegdist(all_proxy, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
all_proxy.chclust <- chclust(all_proxy.dist , method="coniss")
plot(all_proxy.chclust, main = "Constrained Hierarchical Clustering Dendrogram - All Proxy")
#################################
biomarkers_proxy <- sheldon_db %>%
  select("Age","IPSO25_ug_to_g_sed","HBI_III_ug_to_g_sed") %>%
  mutate(Age = as.integer(Age)) %>%
  column_to_rownames(var = "Age")

biomarkers_proxy.dist <- vegdist(biomarkers_proxy, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
biomarkers_proxy.chclust <- chclust(biomarkers_proxy.dist, method="coniss")
plot(biomarkers_proxy.chclust, main = "Constrained Hierarchical Clustering Dendrogram - Biomarkers")
#################################
DNA_proxy <- sheldon_db %>%
  select("Age", "Diatoms_rel_abundance", "Dinoflagellata_rel_abundance", "Polarella_rel_to_Dinoflagellata", "Islandinium_rel_to_Dinoflagellata") %>%
  mutate(Age = as.integer(Age)) %>%
  column_to_rownames(var = "Age")

DNA_proxy.dist <- vegdist(DNA_proxy, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
DNA_proxy.chclust <- chclust(DNA_proxy.dist , method="coniss")
plot(DNA_proxy.chclust, main = "Constrained Hierarchical Clustering Dendrogram - DNA proxy")

#To objectively select number of clusters on hierarchical trees we calculated stick-bone models
#Significant number of clusters is indicated by black line above red line. 

all_proxy_bstick <- bstick(all_proxy.chclust)
all_proxy_bstick

biomarkers_bstick <- bstick(biomarkers_proxy.chclust)
biomarkers_bstick

DNA_proxy_bstick <- bstick(DNA_proxy.chclust)
DNA_proxy_bstick

#To test homogeneity of sea-ice proxies variance across samples we used Levene's Test:
leveneTest(Polarella_rel_to_Dinoflagellata ~ Coniss_all_proxy, sheldon_db)
leveneTest(Islandinium_rel_to_Dinoflagellata ~ Coniss_all_proxy, sheldon_db)
leveneTest(Dinoflagellata_rel_abundance ~ Coniss_all_proxy, sheldon_db)
leveneTest(Diatoms_rel_abundance ~ Coniss_all_proxy, sheldon_db)
leveneTest(IPSO25_ug_to_g_sed ~ Coniss_all_proxy, sheldon_db)
leveneTest(HBI_III_ug_to_g_sed ~ Coniss_all_proxy, sheldon_db)
leveneTest(BiogenicSiO2 ~ Coniss_all_proxy, sheldon_db)
leveneTest(TOC ~ Coniss_all_proxy, sheldon_db)

# To test significance of difference between sea-ice proxies and previously selected clusters, here so-called: "warmer" and "colder" 
kruskal.test(Polarella_rel_to_Dinoflagellata ~ Coniss_all_proxy, sheldon_db)
kruskal.test(Islandinium_rel_to_Dinoflagellata  ~ Coniss_all_proxy, sheldon_db)
kruskal.test(Dinoflagellata_rel_abundance ~ Coniss_all_proxy, sheldon_db)
kruskal.test(Diatoms_rel_abundance ~ Coniss_all_proxy, sheldon_db)
kruskal.test(IPSO25_ug_to_g_sed~ Coniss_all_proxy, sheldon_db)
kruskal.test(HBI_III_ug_to_g_sed ~ Coniss_all_proxy, sheldon_db)
kruskal.test(BiogenicSiO2 ~ Coniss_all_proxy, sheldon_db)
kruskal.test(TOC ~ Coniss_all_proxy, sheldon_db)

#To visualise differences between sea-ice proxies in "warmer" and "colder" clusters
Polarella <- sheldon_db %>%
  select("Coniss_all_proxy", "Polarella_rel_to_Dinoflagellata_rarefied") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

Islandinium <- sheldon_db %>%
  select("Coniss_all_proxy", "Islandinium_rel_to_Dinoflagellata_rarefied") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

Islandinium <- sheldon_db %>%
  select("Coniss_all_proxy", "Islandinium_rel_to_Dinoflagellata_rarefied") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

Dinoflagellates <- sheldon_db %>%
  select("Coniss_all_proxy", "Dinoflagellata_rel_abundance_rarefied") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

Diatoms <- sheldon_db %>%
  select("Coniss_all_proxy", "Diatoms_rel_abundance_rarefied") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

IPSO25 <- sheldon_db %>%
  select("Coniss_all_proxy", "IPSO25_ug_to_g_sed") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

HBI_III <- sheldon_db %>%
  select("Coniss_all_proxy", "HBI_III_ug_to_g_sed") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

BiogenicSiO2 <- sheldon_db %>%
  select("Coniss_all_proxy", "BiogenicSiO2") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

TOC <- sheldon_db %>%
  select("Coniss_all_proxy", "TOC") %>%
  pivot_longer(cols = -Coniss_all_proxy, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=Coniss_all_proxy)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("warmer" = "#ffcc99ff", "colder" = "#80b3ffff")) +
  theme_minimal()

combined_plots <- ggarrange(Polarella, Islandinium, Dinoflagellates,Diatoms, BiogenicSiO2,TOC, IPSO25,HBI_III, nrow=1, ncol=8, common.legend = TRUE, align = "hv")
combined_plots
