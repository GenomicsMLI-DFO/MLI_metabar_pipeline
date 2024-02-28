# Info --------------------------------------------------------------------

# Taxonomic assignment with RDP Classifier
#  
# Audrey Bourret 
# 2023-12-18#

# Library -----------------------------------------------------------------

#library(readr)
library(magrittr)
library(dplyr)
library(readr)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
source(file.path(here::here(), "01_Code", "Functions", "rdp.R"))

# Set output directory ----------------------------------------------------

res.path <- file.path(here::here(), "02_Results/03_TaxoAssign/02_RDP")
if(!file.exists(res.path)) dir.create(res.path)

# Dataset -----------------------------------------------------------------

rdp.path <- get.value("rdp.path")

LOCUS <- stringr::str_split(get.value("Loci"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analyse (", paste(LOCUS, collapse = ", "),")\nTheses parameters can be changed with the file Option.txt", sep = " ")

numCores <- as.numeric(as.character(get.value("NumCores")))

cat(numCores, "core(s) will be used by the script",
    "\nThis parameter can be changed with the file Option.txt", sep = " ")

if(!file.exists(rdp.path)){
  cat("Please set a valid rdp path in the file Option.txt", sep = " ") 
} else("RDP path is valid :)")

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )

# Run RDP -------------------------------------------------------------

# Add python env to this specific project
Sys.setenv(PATH = paste(c("/home/genobiwan/miniconda3/envs/OpenDNA/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))


system2("rdp_classifier", "--help")

PARAM.RDP <- readr::read_tsv(file.path(here::here(), "01_Code/Parameters/rdp_param.tsv"))
PARAM.RDP



# Perform rdp

for(l in LOCUS){
  
  cat("\nWorking on" , l, "\n")
  
  quick.rdp(fasta.file = file.path(here::here(),"00_Data/03c_ESV", paste0("ESV.",l,".fasta")), 
            out.path = file.path(res.path, paste0("rdp.",l, ".out")),
            hier.path= file.path(res.path, paste0("hierarchy.",l, ".out")),
            training.path = file.path(rdp.path ,PARAM.RDP %>% dplyr::filter(Locus == l) %>% dplyr::pull(TrainingSet)))
  
}


# Load results

for(l in LOCUS){
  
  assign(x = paste0("RES.",l,".rdp"), 
         value = readr::read_tsv(file.path(res.path, paste0("rdp.",l, ".out")), col_names = F)
         
  )
}

# Compile at different threshold ---------------------------------------------------------

RES.all.rdp <- tibble()

for(l in LOCUS){
  
  for(t in c(30,50,80)){
    
    rdp.int <- get(paste0("RES.",l,".rdp")) %>% taxon.thr.RDP(threshold = t) %>% 
      dplyr::mutate(Loci = l,
                    Script = "RDP",
                    Method = "RDP",
                    Threshold = t,
                    RefSeq = PARAM.RDP %>% dplyr::filter(Locus == l) %>% dplyr::pull(TrainingSet) %>% stringr::str_remove("/rRNAClassifier.properties"))
    
    readr::write_csv(rdp.int, file = file.path(res.path, paste0("RDP.", t, ".", l, ".csv")))
    rdp.int$confidence_order = as.numeric(as.character(rdp.int$confidence_order))
    rdp.int$confidence_genus = as.numeric(as.character(rdp.int$confidence_genus))
    rdp.int$confidence_species = as.numeric(as.character(rdp.int$confidence_species))
    RES.all.rdp <- bind_rows(RES.all.rdp, rdp.int)  
    
  }
}

readr::write_csv(RES.all.rdp, file = file.path(res.path, paste0("RES.all.rdp.csv")))


# LOAD
ESV.taxo.ALL <- readr::read_csv(file = file.path(here::here(), "02_Results/03_TaxoAssign", paste0("ESV.taxo.ALL.csv")))


# Combine all datasets, but dataset by dataset to be sure to keep unassigned
FINAL_RES <- dplyr::tibble()

for(m in c("RDP") ){
  
  for(t in c(30,50, 80)){
    
    RefSeq.int <- RES.all.rdp %>% filter(Method == m, Threshold == t) %>% pull(RefSeq) %>% unique()
    
    RES.all.rdp.int <- RES.all.rdp %>% filter(Method == m, Threshold == t) %>% select(-c(Loci,Method, Threshold)) %>% distinct(.keep_all = T)
    
    FINAL_RES.int <- ESV.taxo.ALL %>% left_join(RES.all.rdp.int,
                                                by =  c("ESV" = "ESV")) %>% 
      mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon),
             Script = "RDP",
             Method = m, 
             Threshold = t,
             RefSeq = RefSeq.int)
    
    FINAL_RES <- bind_rows(FINAL_RES, FINAL_RES.int)
    
  }
  
}

readr::write_csv(FINAL_RES, file = file.path(res.path, paste0("RES.all.rdp.wSamples.csv")))

save(list = c("RES.all.rdp"),
   file = file.path(res.path, "ESVtab_assign.Rdata"  ))


# Basic figures -----------------------------------------------------------

library(ggplot2)

# Differencre between LCA and Tophit at different threshold
graph1 <- FINAL_RES %>% mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon)) %>% 
  group_by(Loci, phylum, Taxon, Method, Threshold) %>% summarise(Nreads = sum(Nreads)) %>% 
  filter(Nreads > 0) %>% 
  mutate(Method.thresh =paste0(Method, Threshold)) %>% 
  #filter(Method.thresh != "NANA") %>% 
  ggplot(aes(x = Method.thresh, y = Taxon, fill = Nreads)) + 
  geom_bin2d() +
  scale_fill_distiller(palette = "Spectral", trans = "log10", na.value = "white") +
  facet_grid(phylum ~Loci, scale = "free", space = "free") +
  theme_bw() +
  labs(title = "Visual comparison of BLAST assignments")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  )

graph1

ggsave(filename = file.path(res.path, "Comparison_RDP.png"), plot =  graph1, height = 8, width = 10)



# Write a final log

cat("\nEND of 03_TaxoAssign_RDP_classifier.R script\n",
    date(),
    "\n-------------------------\n", 
    
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("Biostrings", packageVersion("Biostrings"), sep = ": "),     
    
    "\n~ External programs ~",
    paste("RDP", system2("which", "rdp_classifier", stdout=T, stderr=T), collapse = (": ")),     
    paste("with the reference db:", paste(PARAM.RDP %>% dplyr::pull(TrainingSet) %>% str_remove("/rRNAClassifier.properties") %>% unique(), collapse = ", ")),
    
    # Add it to the log file
    file = file.path(res.path, "TaxoAssign_RDP.log"), 
    append = F, sep = "\n")

