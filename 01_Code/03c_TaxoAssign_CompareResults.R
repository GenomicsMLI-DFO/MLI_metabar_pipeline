# Info --------------------------------------------------------------------

# Compare all taxonomic assignment that were performed
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


# Load assignment results -------------------------------------------------

res.dirs <- list.dirs("02_Results/03_TaxoAssign", full.names = F)[-1] 


res.dirs

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
  summarise(N_ESV = n()) %>% 
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



ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_ALL_nESV.png"), plot =  graph.ESV, height = 8, width = 10)

ggsave(filename = file.path("02_Results/03_TaxoAssign", "Comparison_ALL_nLoci.png"), plot =  graph.loci, height = 8, width = 10)

