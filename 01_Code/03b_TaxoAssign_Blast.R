# Info --------------------------------------------------------------------

# Simplest taxonomic assignment
# Designed to work with NCBI nt but could work with a local DB
# Blast + both top hit vs LCA
# 
# Audrey Bourret
# 2022-2023
#

# Library -----------------------------------------------------------------

library(readr)
library(magrittr)
library(dplyr)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
source(file.path(here::here(), "01_Code", "Functions", "blast.R"))


res.path <- file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast_nt")

if(!file.exists(res.path)) dir.create(res.path)

# Dataset -----------------------------------------------------------------

NCBI.path <- get.value("NCBI.nt.path")
NCBI.path

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



# BLAST query -------------------------------------------------------------

# Check that we can run blastX

system2("blastn", "-help")
system2('blastn', '-version')

PARAM.BLAST <- readr::read_tsv(file.path(here::here(), "01_Code/Parameters/blast_param.tsv"))
PARAM.BLAST$evalue <- as.character(PARAM.BLAST$evalue )
PARAM.BLAST

# Load taxonomy data

if(PARAM.BLAST$db %>% unique() == "nt"){
   ncbi.tax <- readr::read_tsv(file.path(NCBI.path,"rankedlineage.dmp"), 
   col_names = c("id", "name", "species", "genus", "family", "order", "class","phylum", "kingdom", "superkingdom"), 
   col_types=("i-c-c-c-c-c-c-c-c-c-"))
  
} else{
  ncbi.tax <- readr::read_tsv(file.path(NCBI.path,"rankedlineage.dmp"))
  }



ncbi.tax %>% head()

# Perform blast

for(l in LOCUS){
  
  cat("\nWorking on " , l, "\n")
  
  quick.blastn(fasta.file = file.path(here::here(),"00_Data/03c_ESV", paste0("ESV.",l,".fasta")), 
               db = get.blast.value(l, "db", PARAM.BLAST),
               out.file = file.path(res.path, paste0("Blast.",l, ".raw.out")),
               perc_identity = get.blast.value(l, "perc_identity", PARAM.BLAST), 
               qcov_hsp_perc = get.blast.value(l, "qcov_hsp_perc", PARAM.BLAST), 
               max_target_seqs = get.blast.value(l, "max_target_seqs", PARAM.BLAST),
               evalue = get.blast.value(l, "evalue", PARAM.BLAST),
               NCBI.path = NCBI.path,
               n.cores = numCores)
  
  
  
  
}

# Load results and add taxonomical information

for(l in LOCUS){
  
  assign(x = paste0("RES.",l,".ncbi"), 
         value = load.blast(out.file = file.path(res.path, paste0("Blast.",l, ".raw.out")),
                            ncbi.tax = ncbi.tax,
                            db = get.blast.value(l, "db", PARAM.BLAST))
         
  )
  
  
}


# Compute TOPHIT + LCA ---------------------------------------------------------

RES.all.ncbi <- tibble()

for(l in LOCUS){
  
  for(t in c(95,97,99)){
    # TOP HIT
    TOP.int <- get(paste0("RES.",l,".ncbi")) %>% BLAST_TOPHIT(threshold = t) %>% 
      sum.BLAST() %>% 
      dplyr::mutate(Loci = l,
                    Script = "Blast",
                    Method = "TOP",
                    Threshold = t,
                    RefSeq = get.blast.value(l, "db", PARAM.BLAST))
    
    readr::write_csv(TOP.int, file = file.path( res.path, paste0("TopHit.", t, ".", l, ".csv")))
    
    # LCA
    LCA.int <- get(paste0("RES.",l,".ncbi")) %>% BLAST_LCA(threshold = t) %>% 
      sum.BLAST() %>% 
      dplyr::mutate(Loci = l,
                    Script = "Blast",
                    Method = "LCA",
                    Threshold = t,
                    RefSeq = get.blast.value(l, "db", PARAM.BLAST))
    
    readr::write_csv(LCA.int, file = file.path( res.path, paste0("LCA.", t, ".", l,  ".csv")))
    
    RES.all.ncbi <- bind_rows(RES.all.ncbi, TOP.int, LCA.int)  
    
  }  
  
}

readr::write_csv(RES.all.ncbi, file = file.path(res.path, paste0("RES.all.ncbi.csv")))



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

readr::write_csv(ESV.taxo.ALL, file = file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast_nt", paste0("ESV.taxo.ALL.csv")))

# Combine all datasets, but dataset by dataset to be sure to keep unassigned

FINAL_RES <- dplyr::tibble()

for(m in c("LCA", "TOP") ){
  
  for(t in c(95,97,99)){
    
    RefSeq.int <- RES.all.ncbi %>% filter(Method == m, Threshold == t) %>% pull(RefSeq) %>% unique()
    
    RES.all.ncbi.int <- RES.all.ncbi %>% filter(Method == m, Threshold == t) %>% select(-c(Loci,Method, Threshold, RefSeq)) %>% distinct(.keep_all = T)
    
    FINAL_RES.int <- ESV.taxo.ALL %>% left_join(RES.all.ncbi.int,
                                                by =  c("ESV" = "QueryAccVer")) %>% 
      mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon),
             Script = "Blast",
             Method = m, 
             Threshold = t,
             RefSeq = RefSeq.int)
    
    FINAL_RES <- bind_rows(FINAL_RES, FINAL_RES.int)
    
  }
  
}

readr::write_csv(FINAL_RES, file = file.path(res.path, paste0("RES.all.ncbi.wSamples.csv")))


#FINAL_RES %>% group_by(ESV) %>% summarise(N = n()) %>% arrange(desc(N))

save(list = c("RES.all.ncbi"),
     file = file.path(res.path, "ESVtab_assign.Rdata"  ))


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

ggsave(filename = file.path(res.path, "Comparison_Blast.png"), plot =  graph1, height = 8, width = 10)

# Write a final log

cat("\nEND of 03_TaxoAssign_Blast.R script\n",
    paste("Pipeline version:", get.value("MLI.version")),
    date(),
    "\n-------------------------\n", 
    
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("Biostrings", packageVersion("Biostrings"), sep = ": "),     
   
    "\n~ External programs ~",
    paste(system2("blastn", "-version", stdout=T, stderr=T), collapse = ("; ")),     
    paste("with the reference db:", get.blast.value(l, "db", PARAM.BLAST)),
     # Add it to the log file
    file = file.path( res.path, "TaxoAssign_Blast.log"), 
    append = F, sep = "\n")


