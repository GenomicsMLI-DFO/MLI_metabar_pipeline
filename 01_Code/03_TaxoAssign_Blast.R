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
#numCores <- 20

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

PARAM.BLAST <- readr::read_tsv(file.path(here::here(), "01_Code/Parameters/blast_param.tsv"))
PARAM.BLAST$evalue <- as.character(PARAM.BLAST$evalue )
PARAM.BLAST

res.path <- file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast")

if(!file.exists(res.path)) dir.create(res.path)

for(l in LOCUS){
  
  cat("\nWorking on " , l, "\n")
  
  assign(x = paste0("RES.",l,".ncbi"), 
         value = quick.blastn(fasta.file = file.path(here::here(),"00_Data/03c_ESV", paste0("ESV.",l,".fasta")), 
                              out.file = file.path(res.path, paste0("Blast.",l, ".raw.out")),
                              perc_identity = get.blast.value(l, "perc_identity", PARAM.BLAST), 
                              qcov_hsp_perc = get.blast.value(l, "qcov_hsp_perc", PARAM.BLAST), 
                              max_target_seqs = get.blast.value(l, "max_target_seqs", PARAM.BLAST),
                              evalue = get.blast.value(l, "evalue", PARAM.BLAST),
                              NCBI.path = NCBI.path,
                              n.cores = numCores,
                              ncbi.tax = ncbi.tax)
         
  )
  
  
}


# Compute TOPHIT + LCA ---------------------------------------------------------

library(dplyr)

RES.all.ncbi <- tibble()

for(l in LOCUS){

for(t in c(95,97,99)){
  # TOP HIT
  TOP.int <- get(paste0("RES.",l,".ncbi")) %>% BLAST_TOPHIT(threshold = t) %>% 
                sum.BLAST() %>% 
                dplyr::mutate(Loci = l,
                              Method = "TOP",
                              Threshold = t)
  
  readr::write_csv(TOP.int, file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast", paste0("TopHit.", t, ".", l, ".csv")))

  # LCA
  LCA.int <- get(paste0("RES.",l,".ncbi")) %>% BLAST_TOPHIT(threshold = t) %>% 
                sum.BLAST() %>% 
                dplyr::mutate(Loci = l,
                              Method = "LCA",
                              Threshold = t)

  readr::write_csv(LCA.int, file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast", paste0("LCA.", t, ".", l,  ".csv")))

RES.all.ncbi <- bind_rows(RES.all.ncbi, TOP.int, LCA.int)  
  
  }  
  
}


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


ESV.taxo.ALL <- tibble()

for(l in LOCUS){

  ESV.taxo.ALL <- bind_rows(ESV.taxo.ALL, 
                       tidy.ESV( get(paste0("ESVtab.",l)),   get(paste0("DNA.",l))) %>% mutate(Loci = l) 
                                     )
  
  
}


# Combine all datasets, but dataset by dataset to be sure to keep unassigned

FINAL_RES <- dplyr::tibble()

for(m in c("LCA", "TOP") ){
  
  for(t in c(95,97,99)){

    RES.all.ncbi.int <- RES.all.ncbi %>% filter(Method == m, Threshold == t) %>% select(-c(Loci,Method, Threshold)) %>% distinct(.keep_all = T)

    FINAL_RES.int <- ESV.taxo.ALL %>% left_join(RES.all.ncbi.int,
                                              by =  c("ESV" = "QueryAccVer")) %>% 
                                      mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon),
                                             Method = m, 
                                             Threshold = t)

   FINAL_RES <- bind_rows(FINAL_RES, FINAL_RES.int)

  }
  
}


#FINAL_RES %>% group_by(ESV) %>% summarise(N = n()) %>% arrange(desc(N))

save(list = c("FINAL_RES", "ESV.taxo.ALL", "RES.all.ncbi"),
     file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast/", "ESVtab_assign.Rdata"  ))


# Basic figures -----------------------------------------------------------

library(ggplot2)

# Differencre between LCA and Tophit at different threshold
graph1 <- FINAL_RES %>% mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon)) %>% 
  group_by(Loci, phylum, class, genus, Taxon, Method, Threshold) %>% summarise(Nreads = sum(Nreads)) %>% 
  filter(Nreads > 0) %>% 
   mutate(Method.thresh =paste0(Method, Threshold) ) %>% 
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

ggsave(filename = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast", "Comparison_Blast.png"), plot =  graph1, height = 8, width = 10)

