
# Info --------------------------------------------------------------------

# eDNA pipeline using DADA2
# Template pipeline
# 
# Audrey Bourret
# 2021-09-11
#


# Library -----------------------------------------------------------------

# Data manipulation

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(here)
library(parallel)
library(ggplot2)

library(dada2)
library(Biostrings)

#library(tidyverse) # includes ggplot2 dplyr and stringr

#library(gtools)    # for mixedsort
#library(readxl)

#library(devtools)
#devtools::install_github("kassambara/ggpubr")
#library(ggpubr)    # on github - for nice graphs

# Fastq and fasta manipulation


# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# BiocManager::install("dada2")

#library(Biostrings)
#library(dada2); packageVersion("dada2") # Faire mettre cette info dans le log

# Internal functions
source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
source(file.path(here::here(), "01_Code", "Functions", "fastqc.R"))
source(file.path(here::here(), "01_Code", "Functions", "cutadapt.R"))
source(file.path(here::here(), "01_Code", "Functions", "dada2.R"))

# Add python env to this specific project
#Sys.setenv(PATH = paste(c("/home/genobiwan/Documents/PythonVenv/GenoBaseEnv/bin",
#                          Sys.getenv("PATH")),
#                        collapse = .Platform$path.sep))

system2("cutadapt", "--help")
system2("multiqc", "--help")

#auto.folder()

# Data --------------------------------------------------------------------

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

LOCUS <- stringr::str_split(get.value("Loci"), pattern = ";")[[1]]
SENS  <- stringr::str_split(get.value("Sens"), pattern = ";")[[1]]
RUN   <- stringr::str_split(get.value("Run"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analyse (", paste(LOCUS, collapse = ", "), ") from the run(s)", paste(RUN, collapse = ", "),
    "\nTheses parameters can be changed with the file Option.txt", sep = " ")

numCores <- as.numeric(as.character(get.value("NumCores")))

cat(numCores, "core(s) will be used by the script",
    "\nThis parameter can be changed with the file Option.txt", sep = " ")

# Setting up the log file

# Functions ---------------------------------------------------------------

# Check that fastqc and mutliqc are found on this computer

system2("fastqc", "--help")
system2("multiqc", "--help")

# 1. RAW quality assesment (fastqc) ---------------------------------------

# FastQC
fastqc(folder.in = file.path(here::here(), "00_Data", "01b_RawData_rename"),
       folder.out = file.path(here::here(), "02_Results", "01_FastQC", "01_Raw"))

multiqc(folder.out = file.path(here::here(), "02_Results", "01_FastQC", "01_Raw"),
        loci = LOCUS, 
        sens = SENS)

# 2. RAW to FILT (cutadapt + dada2) ------------------------------------------------------------

# 2.1 CUTADAPT

# We need to remove the adaptors, and discard reads untrimmed
# this function can work with more than one loci

system2("cutadapt", "--help")

cutadapt(folder.in = file.path(here::here(), "00_Data", "01b_RawData_rename"), 
         folder.out = file.path(here::here(), "00_Data", "02a_Cutadapt"), 
         loci = LOCUS, 
         sens = SENS, 
         numCores = numCores,
         novaseq = FALSE) 

# Running another fastqc following cutadapt
fastqc(folder.in = file.path(here::here(), "00_Data", "02a_Cutadapt"),
       folder.out = file.path(here::here(), "02_Results", "01_FastQC", "02_Cutadapt"),
       numCores = numCores)

multiqc(folder.out = file.path(here::here(), "02_Results", "01_FastQC", "02_Cutadapt"),
        loci = LOCUS, 
        sens = SENS)


# 2.2 DADA2

# Check if these parameters seem good, if not change them
PARAM.DADA2 <- readr::read_tsv(file.path(here::here(), "01_Code/Functions/dada2_param.tsv"))
PARAM.DADA2

dada2.filter (folder.in = file.path(here::here(), "00_Data", "02a_Cutadapt"), 
         folder.out = file.path(here::here(), "00_Data", "02b_Filtered_dada2"), 
         loci = LOCUS, 
         sens = SENS, 
         param.dada2 = PARAM.DADA2,
         numCores = numCores
)
         
COI.res <- readr::read_csv(file.path("00_Data/02b_Filtered_dada2/log/COI_dada2_summary.csv")) 

COI.res <- COI.res %>% mutate(Keep = reads.in,
                   Remove = reads.in - reads.out) %>% 
  select(-c(reads.in, reads.out)) %>% 
  tidyr::pivot_longer(cols = c(Keep, Remove), names_to = "Reads", values_to = "N")  
  
COI.res %>%  mutate(Reads = factor(Reads, levels = c("Remove", "Keep")),
                    ID_labo = ID %>% str_remove("_COI")) %>%
  left_join(data.info) %>% 
  ggplot(aes(x = ID_labo, y = N, fill = Reads)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Type_echantillon, scale = "free")

# QUALITY ASSEMENT

# Running another fastqc following cutadapt
fastqc(folder.in = file.path(here::here(), "00_Data", "02b_Filtered_dada2"),
       folder.out = file.path(here::here(), "02_Results", "01_FastQC", "03_Dada2"),
       numCores = numCores)

multiqc(folder.out = file.path(here::here(), "02_Results", "01_FastQC", "03_Dada2"),
        loci = LOCUS, 
        sens = SENS)


# 3. FILT to ASV (denoise with DADA2) --------------------------------------------------

# 3.1 Compute error rate

# Be careful, MUST include samples from the same RUN (cause errors can be run specific)
# I was used to compute it by loci, but I don't think it's necessary
# except if the error rate have something to do with the position and all markers don't have the same length
# Some part are weird simply because it was in a loop previously

  # create files lists
  filesF.temp <- list.files(file.path(here::here(), "00_Data", "02b_Filtered_dada2"), full.name =T, pattern = ".fastq") %>%
                 str_subset(SENS[1])
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "files were found\n")  
  
  #err.F.temp <- learnErrors(filesF.temp)
  assign(x = paste0("err.F"), value = learnErrors(filesF.temp,
                                                     nbases = 1e8, # 1e8 is the default - increase to sample more samples
                                                     randomize = T,
                                                     MAX_CONSIST = 10, # default - can be more if not reach
                                                     multithread = ifelse(numCores > 1, T, F)))
  
  #if(nrow(PARAM.temp == 2)){
    filesR.temp <-  list.files(file.path(here::here(), "00_Data", "02b_Filtered_dada2"), full.name =T, pattern = ".fastq") %>%
                    str_subset(SENS[2])
    
    # Add a message to be sure that the right number of files were used
    cat(length(filesR.temp), "files were found\n")  
    
    assign(x = paste0("err.R"), value = learnErrors(filesR.temp,
                                                       nbases = 1e8, # 1e8 is the default
                                                       randomize = T,
                                                       MAX_CONSIST = 10, # default - can be more if not reach
                                                       multithread = ifelse(numCores > 1, T, F)))
    
  #} else {err.R.temp <- vector()}
  
  #  cat("\nPlotting error rate results for" , l, "\n")  
  
  # Print a PDF of error rate
    
  pdf(file.path(here::here(), "00_Data", "03a_ErrorRate_dada2", "ErrorsRate.pdf")) 
    print(plotErrors(get(paste0("err.F")), nominalQ=TRUE))
    print(plotErrors(get(paste0("err.R")), nominalQ=TRUE))
  dev.off()

  #cat("\nSaving error rate results for" , l, "\n") 
    
save(list = paste0("err.F"),
     file = file.path(here::here(), "00_Data", "03a_ErrorRate_dada2", paste0("err.F.Rdata")))
  
#if(nrow(PARAM.temp == 2)){ 
  save(list = paste0("err.R"),
       file = file.path(here::here(), "00_Data", "03a_ErrorRate_dada2",paste0("err.R.Rdata")))
#  }  


# 3.2 Dereplication and sample inference

# Since we are working with BIG DATA, we don't have enough memory
# to run it by locus most of the time. So this code run it by sample directly.

# To reload error rate if necessary  
for(x in list.files(file.path(here::here(), "00_Data", "03a_ErrorRate_dada2"), full.names = T, pattern = ".Rdata")){
  load(x)
}

for(l in LOCUS){
  
  cat("\nWorking on " , l, "\n")

  # create files lists
  filesF.temp <- list.files(file.path(here::here(), "00_Data", "02b_Filtered_dada2"), full.name =T, pattern = ".fastq") %>%
    str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
    str_subset(SENS[1] %>% as.character())
  
  #filesR.temp <- filesF.temp %>% str_replace(SENS[1], SENS[2])
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "files were found\n")  

  mergers <- vector("list", length(filesF.temp))
  names(mergers) <- filesF.temp

  # Set a progress bar
  pb <- txtProgressBar(min = 0, max = length(filesF.temp), style = 3)
  
  for(i in seq_along(filesF.temp)){
    
    samF <- filesF.temp[i]
    samR <- samF %>% str_replace(SENS[1] %>% as.character(), SENS[2] %>% as.character())
    
    cat("\nProcessing:", samF ,"\n" )

    cat("Dereplication\n" )
    
    # Dereplication
    derep.F <-  derepFastq(samF)
    derep.R <-  derepFastq(samR)    
    
    cat("Sample inference\n" )
    
    #Sample inference
    dada.F <-  dada(derep.F, 
                     err =  get(paste0("err.F")), 
                     multithread = ifelse(numCores > 1, T, F),
                     pool=TRUE)
    
    dada.R <-  dada(derep.R, 
                     err =  get(paste0("err.R")), 
                     multithread = ifelse(numCores > 1, T, F),
                     pool=TRUE)

    cat("Merge paires\n" )
    
    merger <- mergePairs(dadaF = dada.F, 
                              derepF = derep.F, 
                              dadaR = dada.R, 
                              derepR = derep.R, 
                              minOverlap = 30, 
                              maxMismatch = 0,
                              returnRejects = FALSE,
                              verbose=TRUE)
    
    
      mergers[[samF]] <- merger
  
  setTxtProgressBar(pb, i)
      
  }
 
  rm(list = c("derep.F", "derep.R", "dada.F", "dada.R", "merger"))
  
  names(mergers) <-   names(mergers) %>% str_remove(file.path(here::here(), "00_Data", "02b_Filtered_dada2")) %>% 
                                         str_remove(paste0("_", l, "_R1")) %>% 
                                         str_remove("_cutadapt") %>% 
                                         str_remove(".fastq.gz") %>% 
                                         str_remove("/")
  
  cat("\nMaking SeqTab for" , l, "\n")
  
  #seqtab <- makeSequenceTable(mergers)
  
  assign(x = paste0("seqtab.", l, ".int"), value = makeSequenceTable(mergers))
  
  cat("\nSaving SeqTab for" , l, "\n") 
  
  save(list = paste0("seqtab.", l, ".int"),
       file = file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("seqtab.wCHIM.", l, ".Rdata")  ))
  
  close(pb)
  
 }



# 3.6 Remove chimera 

for(l in LOCUS){
  
  cat("\nRemoving chimera for" , l, "\n")
  
  assign(x = paste0("ESVtab.", l), value = removeBimeraDenovo(get(paste0("seqtab.",l,".int")), method = "consensus", 
                                                              multithread = ifelse(numCores > 1, T, F), verbose = TRUE)
  )
  
  cat("\nSaving SeqTab without chimera for" , l, "\n") 
  
  save(list = paste0("ESVtab.", l),
       file = file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("ESVtab.noCHIM.", l, ".Rdata")  ))

  
}


# Save ESV table and fasta file --------------------------------------


for(l in LOCUS){

  write.dada2.res(ESVtab = get(paste0("ESVtab.",l)), 
                  loci = l, 
                  folder = file.path(here::here(), "00_Data", "03c_ESV"))

}

# Write a final log
  
cat("\nEND of 02_Process_RAW.R script\n",
    date(),
    "\n-------------------------\n", 
   
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("dada2", packageVersion("dada2"), sep = ": "),     
    #paste("fastqcr", packageVersion("fastqcr"), sep = ": "),     
    paste("Biostrings", packageVersion("Biostrings"), sep = ": "),   

    "\n~ External programs ~",
    paste("fastqc", system2("fastqc", "-v", stdout=T, stderr=T), sep = ": "),     
    paste("multiqc", system2("multiqc", "--version", stdout=T, stderr=T), sep = ": "),     
    paste("cutadapt", system2("cutadapt", "--version", stdout=T, stderr=T), sep = ": "),  
    
    # Add it to the log file
    file = file.path(here::here(), "00_Data", "03c_ESV", "Process_RAW.log"), 
    append = F, sep = "\n")
