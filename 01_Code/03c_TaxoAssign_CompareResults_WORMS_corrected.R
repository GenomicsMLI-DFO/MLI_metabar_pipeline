# Info --------------------------------------------------------------------

# Compare all taxonomic assignment that were performed
# and create a final taxonomic assignment result
#
# Audrey Bourret
# 2024
#

# Library -----------------------------------------------------------------

library(stringr)
library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
`%nin%` = Negate(`%in%`)

# Load assignment results -------------------------------------------------

# This script will check across all folders available
res.dirs <- list.dirs("02_Results/03_TaxoAssign", full.names = F)[-1] 
res.dirs

# Compile the results
res <- tibble()

for(i in res.dirs){

res.file <- list.files(file.path("02_Results/03_TaxoAssign",i), pattern = "RES.all", full.names = T) %>% 
              str_subset("wSamples", negate = T)  

res.int <- read_csv(res.file) %>% mutate(Folder = i)

if("QueryAccVer" %in% names(res.int)){
  res.int <- res.int %>% dplyr::rename(ESV = QueryAccVer )
}

res <- bind_rows(res, res.int)

}

write_csv(res , "02_Results/03_TaxoAssign/Assignements.ALL.csv")

# Graphics ----------------------------------------------------------------

graph.ESV <- res %>% dplyr::filter(Levels %in% c("species")) %>% 
  dplyr::group_by(Taxon, phylum, Levels, Loci,  Method, Threshold, Folder, RefSeq ) %>% 
  dplyr::summarize(N_ESV = n()) %>%
  mutate(CAT = paste(Method, Threshold)) %>% 
  ggplot(aes(y = Taxon, x = CAT))+
  geom_bin2d(aes(fill = N_ESV)) +
  scale_fill_viridis_c() +
  facet_grid(phylum~RefSeq + Method, scale = "free", space = "free")+
theme_bw() +
  labs(title = "Visual comparison detected species",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.ESV

graph.loci <- res %>% dplyr::filter(Levels %in% c("species")) %>% 
  dplyr::group_by(Taxon, phylum, Levels, Method, Threshold, Folder, RefSeq ) %>% 
  summarise(N_loci = length(unique(Loci))) %>% 
  mutate(CAT = paste(Method, Threshold)) %>% 
  ggplot(aes(y = Taxon, x = CAT))+
  geom_bin2d(aes(fill = N_loci)) +
  scale_fill_viridis_c() +
  facet_grid(phylum~RefSeq + Method, scale = "free", space = "free")+
  theme_bw() +
  labs(title = "Visual comparison detected species",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.loci



ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_ALL_nESV.png"), plot =  graph.ESV, height = 12, width = 10)

ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_ALL_nLoci.png"), plot =  graph.loci, height = 12, width = 10)


# Final dataset -----------------------------------------------------------

# The first method that will be used (usually a local DB)
method.1 <- get.value("assign.1")
method.1

# The second method that will be used to "fill" the hole
method.2 <- get.value("assign.2")
method.2

# Subset the dataset with the preferred assignment method
res.1 <- res %>% dplyr::filter(Method ==  stringr::str_split(method.1, pattern = ";")[[1]][2] ,
                               Threshold == stringr::str_split(method.1, pattern = ";")[[1]][3],
                               Script  == stringr::str_split(method.1, pattern = ";")[[1]][1]) %>% 
                 mutate(Assignment.method = method.1)

res.1 %>% group_by(Loci) %>% summarise(N = n())

# Subset the dataset with the second assignment method
res.2 <- res %>% dplyr::filter(Method ==  stringr::str_split(method.2, pattern = ";")[[1]][2] ,
                               Threshold == stringr::str_split(method.2, pattern = ";")[[1]][3],
                               Script  == stringr::str_split(method.2, pattern = ";")[[1]][1]) %>% 
  mutate(Assignment.method = method.2)

res.2 %>% group_by(Loci) %>% summarise(N = n())

res.2 %>% dplyr::filter(ESV %nin% res.1$ESV) %>%  group_by(Loci) %>% summarise(N = n())


res.final <- bind_rows(res.1, 
                       res.2 %>% dplyr::filter(ESV %nin% res.1$ESV) )
res.final
# Check that no duplicated values persisted
table(duplicated(res.final$ESV))

# Save the results

write_csv(res.final,  "02_Results/03_TaxoAssign/Assignements.not.corrected.csv")



graph.ESV.final <- res.final %>% dplyr::filter(Levels %in% c("species", "genus")) %>% 
  dplyr::group_by(Taxon, phylum, Levels, Loci, Assignment.method ) %>% 
  summarise(N_ESV = n()) %>% 
  dplyr::mutate(Assignment.method = factor(Assignment.method, levels = c(method.1, method.2))) %>% 
   ggplot(aes(y = Taxon, x = Assignment.method))+
  geom_bin2d(aes(fill = N_ESV)) +
  scale_fill_viridis_c() +
  facet_grid(phylum~Loci, scale = "free", space = "free")+
  theme_bw() +
  labs(title = "Detected species combining 2 methods",
       x = "Assignement method")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph.ESV.final

ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_Final_nESV.png"), plot =  graph.ESV.final, height = 12, width = 10)


# Check and correct taxonomy for WoRMS ------------------------------
#import function validateWORMS
source(file.path(here::here(), "01_Code", "Functions", "fct_WoRMS.R"))

#load Final assignation file
assignations <- readr::read_csv(file.path("02_Results/03_TaxoAssign/Assignements.not.corrected.csv"))

#remove duplicate Taxon and NA
assignations.unique <- assignations %>% distinct(Taxon, .keep_all = T) %>% 
  select(species,Taxon, Levels, genus, family, order, class, phylum, kingdom, ESV)
assignations.unique <- assignations.unique %>% drop_na(Taxon)

#check if a taxon or a vector of taxon complies with WORMS ID
list.taxon <- assignations.unique$Taxon

valid.taxo <- data.frame()

sep <- seq(from = 1, to = length(list.taxon), by = 10) #separate list to avoid error with taxa without entry in WORMS
sep

for(x in sep){
  cat("\nWe now process row number", x)
  y <- x + 9
  
  if(y > length(list.taxon)){
    y <- length(list.taxon)
  }
  
  valid.taxo.int <- validateWORMS(list.taxon[c(x:y)])
  
  valid.taxo <- bind_rows(valid.taxo, valid.taxo.int)
  
}

write.csv(valid.taxo, file = file.path(here::here(), "02_Results/03_TaxoAssign/Assignements_WoRMS_valid.csv"))

#select and check only the taxa with the "unaccepted" status
unacc.taxo <- valid.taxo %>% filter(status %in% c("unaccepted", "accepted | unaccepted"))

valid.taxo.f <- valid.taxo %>% separate(valid_name, c("valid_name", "B"),sep = "\\|", remove = TRUE) %>% #sep is a regular expression so the | needs to be escaped with \\
                               separate(kingdom,  c("kingdom", "C"),sep = "\\|",remove=TRUE) %>%
                               separate(phylum,  c("phylum", "D"),sep = "\\|", remove =TRUE) %>%
                               separate(class,  c("class", "E"),sep = "\\|", remove=TRUE) %>%
                               separate(order,  c("order", "F"),sep = "\\|", remove=TRUE) %>%
                               separate(family, c("family", "G"),sep = "\\|", remove =TRUE) %>%
                               select(-c(B,C,D,E,F,G)) %>%
                               mutate(species = case_when(rank == "Species" ~ valid_name, rank != "Species" ~ NA) ) %>%
  dplyr::rename(Taxon = scientificname)


#add WORMS into the original Assignement file
assignations.valid <- assignations %>% left_join(valid.taxo.f %>% select(Taxon, AphiaID, authority, status, valid_name, kingdom, phylum, class,order,family,genus, species), by = "Taxon")
assignations.valid <- assignations.valid  %>% separate(valid_name, c("A",NA), remove = F) %>% 
  mutate(genus = case_when(Levels == "species" ~ A, 
                           Levels == "genus" ~ A, 
                           Levels == "family"~ NA,
                           Levels == "order"~ NA,
                           Levels == "class"~ NA,
                           Levels == "phylum"~ NA,
                           Levels == "kingdom"~ NA) ) %>% 
  select(-c(genus.y,Taxon,A)) %>%
  dplyr::rename(kingdom_NCBI = kingdom.x, 
         phylum_NCBI = phylum.x,
         class_NCBI = class.x,
         order_NCBI = order.x,
         family_NCBI = family.x,
         genus_NCBI = genus.x,
         species_NCBI =species.x,
         kingdom = kingdom.y,
         phylum =  phylum.y,
         class = class.y,
         order = order.y,
         family = family.y,
         Taxon = valid_name,
         species = species.y
         #genus = genus.y,
  )

write.csv(assignations.valid, file = file.path(here::here(), "02_Results/03_TaxoAssign/Assignements.Final.csv"))


# https://www.marinespecies.org/aphia.php?p=manual to check the significance of WoRMS status. 







