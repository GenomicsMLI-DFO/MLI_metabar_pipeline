# Info --------------------------------------------------------------------

# Prepare a tidy ESV table to be used with both BLAST and RDP
#  
# Audrey Bourret 
# 2023-12-18

# Library -----------------------------------------------------------------

library(readr)
library(magrittr)
library(dplyr)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))

# Dataset -----------------------------------------------------------------

LOCUS <- stringr::str_split(get.value("Loci"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analyse (", paste(LOCUS, collapse = ", "),")\nTheses parameters can be changed with the file Option.txt", sep = " ")


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

readr::write_csv(ESV.taxo.ALL, file = file.path(here::here(), "02_Results/03_TaxoAssign/ESV.taxo.ALL.csv"))



