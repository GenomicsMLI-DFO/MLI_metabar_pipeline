# Info --------------------------------------------------------------------

# Simplest taxonomic assignment
# Blast + both top hit vs LCA
# 
# Audrey Bourret
# 2022-02-20
#

# Library -----------------------------------------------------------------

library(readr)
library(magrittr)
library(dplyr)
library(Biostrings)


source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
source(file.path(here::here(), "01_Code", "Functions", "blast.R"))

# Dataset -----------------------------------------------------------------

NCBI.path <- get.value("NCBI.path")

LOCUS <- stringr::str_split(get.value("Loci"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analyse (", paste(LOCUS, collapse = ", "),")\nTheses parameters can be changed with the file Option.txt", sep = " ")

numCores <- as.numeric(as.character(get.value("NumCores")))
numCores <- 20

cat(numCores, "core(s) will be used by the script",
    "\nThis parameter can be changed with the file Option.txt", sep = " ")

if(!file.exists(NCBI.path)){
cat("Please set a valid NCBI path in the file Option.txt", sep = " ") 
}

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

# Load taxonomy data

ncbi.tax <- readr::read_tsv(file.path(NCBI.path,"rankedlineage.dmp"), 
                     col_names = c("id", "name", "species", "genus", "family", "order", "class","phylum", "kingdom", "superkingdom"), 
                     col_types=("i-c-c-c-c-c-c-c-c-c-"))

ncbi.tax %>% head()

# BLAST query -------------------------------------------------------------

# Check that we can run blastX

system2("blastn", "-help")

#PARAM.BLAST <- readr::read_tsv(file.path(here::here(), "01_Code/Parameters/blast_param.tsv"))
#PARAM.BLAST
res.path <- file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast")

if(!file.exists(res.path)) dir.create(res.path)
    
for(l in LOCUS){

  cat("\nWorking on " , l, "\n")
  
    assign(x = paste0("RES.",l,".ncbi"), 
         value = quick.blastn(fasta.file = file.path(here::here(),"00_Data/03c_ESV", paste0("ESV.",l,".fasta")), 
                              out.file = file.path(res.path, paste0("Blast.",l, ".raw.out")),
                              perc_perc_identity = 95, 
                              qcov_hsp_perc = 95, 
                              max_target_seqs = 500, 
                              evalue = "1e-50",
                              NCBI.path = NCBI.path,
                              n.cores = numCores,
                              ncbi.tax = ncbi.tax)
         
  )
  
  
}


# Compute TOPHIT + LCA ---------------------------------------------------------

library(dplyr)

RES.all <- tibble()

for(l in LOCUS){

for(t in c(95,97,99)){
  # TOP HIT
  TOP.int <- get(paste0("RES.",l,".ncbi")) %>% BLAST_TOPHIT(threshold = t) %>% 
                sum.BLAST() %>% 
                dplyr::mutate(Loci = l,
                              Method = "TOP",
                              Threshold = t)
  
  write_csv(TOP.int, file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast", paste0("TopHit.", t, ".", l, ".csv")))

  # LCA
  LCA.int <- get(paste0("RES.",l,".ncbi")) %>% BLAST_TOPHIT(threshold = t) %>% 
                sum.BLAST() %>% 
                dplyr::mutate(Loci = l,
                              Method = "LCA",
                              Threshold = t)

  readr::write_csv(LCA.int, file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast", paste0("LCA.", t, ".", l,  ".csv")))

RES.all <- bind_rows(RES.all, TOP.int, LCA.int)  
  
  }  
  
}

# Load all results




#Loading seq table

for(l in LOCUS){
  
load(file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("ESVtab.noCHIM.", l, ".Rdata")  ))

assign(x = paste0("DNA.", l), 
       value =Biostrings::readDNAStringSet(file.path(here::here(), "00_Data", "03c_ESV",paste0("ESV.", l, ".fasta")))
)
  
}


tidy.ESV <- function(ESVtab, DNA.seq) {
 DNA.tidy <- tibble(ESV = names(DNA.seq), SEQ =  DNA.seq %>% as.character())

ESV.tidy <- ESVtab %>% as_tibble() %>% 
                             dplyr::mutate(ID_labo = row.names(ESVtab)) %>%
                             tidyr::pivot_longer(cols = !ID_labo, names_to = "SEQ", values_to = "Nreads") %>% 
                             dplyr::left_join(DNA.tidy)
 
 return(ESV.tidy) 
}


ESV_RES <- tibble()

for(l in LOCUS){

  ESV_RES <- bind_rows(ESV_RES, 
                       tidy.ESV( get(paste0("ESVtab.",l)),   get(paste0("DNA.",l))) %>% mutate(Loci = l))
  
  
}





FINAL_RES <- ESV_RES %>% left_join(RES.all %>% select(-Loci), by =  c("ESV" = "QueryAccVer"))

FINAL_RES

#save(list = c("FINAL_RES", "ESV_RES", "RES.all"),
#     file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast/", "ESVtab_assign.Rdata"  ))


# Basic figures -----------------------------------------------------------




library(ggplot2)

graph1 <- FINAL_RES %>% mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon)) %>% 
group_by(ID_labo, Loci, phylum, class, genus, Taxon, Method, Threshold) %>% summarise(Nreads = sum(Nreads)) %>% 
  filter(Nreads > 0) %>% 
  left_join(data.info) %>%
  #filter(Type_echantillon %in% c("PPC", "PNC")) %>% 
  ggplot(aes(x = ID_labo, y = Taxon, fill = Nreads)) + 
  geom_bin2d() +
  scale_fill_distiller(palette = "Spectral", trans = "log10", na.value = "white") +
  facet_grid(phylum ~Loci + Method + Threshold, scale = "free", space = "free") +
  theme_bw() +
  #labs(title = "PCR control results")+ 
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

#ggsave(filename = file.path(here::here(), "02_Results/03_Assignments/01_BlastX", "PC_MiFish.png"), plot =  graph3, height = 8, width = 10)

