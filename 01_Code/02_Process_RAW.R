
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
library(here)
library(parallel)
library(ggplot2)

library(dada2)
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
Sys.setenv(PATH = paste(c("/home/genleia/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

system2("cutadapt", "--help")
system2("multiqc", "--help")

#auto.folder()


# Data --------------------------------------------------------------------

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

LOCUS <- stringr::str_split(get.value("Loci"), pattern = ";")[[1]]
SENS  <- stringr::str_split(get.value("Sens"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analysed:", LOCUS, 
    "\nThis parameter can be changed with the file Option.txt", sep = " ")

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

# Check if these paramaters seem good, if not change them
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
  
COI.res %>%   ggplot(aes(x = ID, y = N, fill = Reads)) +
  geom_bar()
# QUALITY ASSEMENT

fastqc(get.value("filt_dada2.path"), get.value("result.FQdada2.path"))

multiqc(get.value("result.FQdada2.path"))

# 3. FILT to ASV (denoise with DADA2) --------------------------------------------------

# 3.1 Calcul du taux d'erreur

# Be careful, MUST include samples from the same RUN (cause error can be run specific)


  # create files lists
  filesF.temp <- list.files(get.value("filt_dada2.path"), full.name =T, pattern = ".fastq") %>%
    #str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
    str_subset(paste0(SENS[1],"_cut"))
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "files were found\n")  
  
  #err.F.temp <- learnErrors(filesF.temp)
  assign(x = paste0("err.F"), value = learnErrors(filesF.temp,
                                                     nbases = 1e8, # 1e8 is the default - increase to sample more samples
                                                     randomize = T,
                                                     MAX_CONSIST = 10, # default - can be more if not reach
                                                     multithread = ifelse(numCores > 1, T, F)))
  
  #if(nrow(PARAM.temp == 2)){
    filesR.temp <- list.files(get.value("filt_dada2.path"), full.name =T, pattern = ".fastq") %>% 
      #str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
      str_subset(paste0(SENS[2] %>% as.character(),"_cut"))
  
    # Add a message to be sure that the right number of files were used
    cat(length(filesR.temp), "files were found\n")  
    
    assign(x = paste0("err.R"), value = learnErrors(filesR.temp,
                                                       nbases = 1e8, # 1e8 is the default
                                                       randomize = T,
                                                       MAX_CONSIST = 10, # default - can be more if not reach
                                                       multithread = ifelse(numCores > 1, T, F)))
    
  #} else {err.R.temp <- vector()}
  
  #  cat("\nPlotting error rate results for" , l, "\n")  
    
  pdf(file.path(get.value("error.dada2.path"), paste0("ErrorsRate.dada2_june2021.pdf"))) 
    print(plotErrors(get(paste0("err.F")), nominalQ=TRUE))
    print(plotErrors(get(paste0("err.R")), nominalQ=TRUE))
  dev.off()

  #cat("\nSaving error rate results for" , l, "\n") 
    
save(list = paste0("err.F"),
     file = file.path(get.value("error.dada2.path"),paste0("err.F.data")  ))
  
#if(nrow(PARAM.temp == 2)){ 
  save(list = paste0("err.R"),
       file = file.path(get.value("error.dada2.path"),paste0("err.R.data")  ))
#  }  



cat("Dada2 error rate assessment was performed and graph were saved:",
    "01_Results/ErrorsRate.dada2.pdf",
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# 3.2 Dereplication and sample inference

# Since we are working with BIG DATA, we don't have enought memory
# to run by locus. This is a first test by sample

for(x in list.files(get.value("error.dada2.path"), full.names = T, pattern = ".data")){
  load(x)
}

for(l in LOCUS){
  
  cat("\nWorking on " , l, "\n")

  # create files lists
  filesF.temp <- list.files(get.value("filt_dada2.path"), full.name =T, pattern = ".fastq") %>%
    str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
    str_subset(paste0(SENS[1] %>% as.character(),"_cut"))
  
  #filesR.temp <- filesF.temp %>% str_replace(SENS[1], SENS[2])
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "files were found\n")  

  mergers <- vector("list", length(filesF.temp))
  names(mergers) <- filesF.temp

  # Set a progress bar
  pb <- txtProgressBar(min = 0, max = length(filesF.temp), style = 3)
  
  for(i in seq_along(filesF.temp)){
    
    samF <- filesF.temp[i]
    samR <- samF %>% str_replace(paste0(SENS[1] %>% as.character(),"_cut"), paste0(SENS[2] %>% as.character(),"_cut"))
    
    cat("\nProcessing:", samF %>%  str_remove(get.value("filt_dada2.path")) %>% 
          str_remove(paste0("_", l, "_R1_cut.fastq.gz")) %>% 
          str_remove("/"), "\n" )

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
  
  names(mergers) <-   names(mergers) %>% str_remove(get.value("filt_dada2.path")) %>% 
                                         str_remove(paste0("_", l, "_R1_cut.fastq.gz")) %>% 
                                         str_remove("/")
  
  cat("\nMaking SeqTab for" , l, "\n")
  
  #seqtab <- makeSequenceTable(mergers)
  
  assign(x = paste0("seqtab.", l, ".int"), value = makeSequenceTable(mergers))
  
  cat("\nSaving SeqTab for" , l, "\n") 
  
  save(list = paste0("seqtab.", l, ".int"),
       file = file.path(get.value("seqtab.dada2.path"),paste0("seqtab.wCHIM.", l, ".data")  ))
  
  close(pb)
  
 }

cat("Dada2 dereplication, sample inference and seqence table were performed.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# 3.6 Remove chimera ------------------------------------------------------

# Why false for multithread?

for(l in LOCUS){
  
  cat("\nRemoving chimera for" , l, "\n")
  
  assign(x = paste0("ESVtab.", l), value = removeBimeraDenovo(get(paste0("seqtab.",l,".int")), method = "consensus", 
                                                              multithread = FALSE, verbose = TRUE)
  )
  
  cat("\nSaving SeqTab without chimera for" , l, "\n") 
  
  save(list = paste0("ESVtab.", l),
       file = file.path(get.value("chimera.dada2.path"),paste0("ESVtab.noCHIM.", l, ".data")  ))

  
}

cat("Chimera were removed.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 



# Write results

write.dada2.res <- function(tab, name, folder){

  file1 <- file.path(folder, paste0("all.", name, "_ESV.fasta"))
  file2 <- file1 %>% str_replace(".fasta", "_ESVtable.txt")
  
  DNA <- DNAStringSet(getSequences(tab))
  names(DNA) <- paste0("ESV_", 1:length(DNA))
  
  writeXStringSet(DNA, file1)
  write.table(tab, file2)
  
}

for(l in LOCUS){

  write.dada2.res(get(paste0("ESVtab.",l)), l, get.value("ASV.dada2.path"))

}

save(file = get.value("ASVtable.data"), 
     list = ls(pattern = "ESVtab."))

cat("Data were saved:",
    get.value("ASVtable.data"),
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# END of the script -------------------------------------------------------

# save.image(file.path(log.path,"Process_RAW.Rdata"))

end.time <- round(Sys.time() - start.time,2) 

cat("\nEND of the raw data processing!",
    #paste0("Rdata saved: ", file.path(log.path,"Process_RAW.Rdata")),
    paste0("\nTime to run this script: ", end.time, " ",units(end.time)),
    "\n-------------------------\n", 
   
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("dada2", packageVersion("dada2"), sep = ": "),     
    #paste("fastqcr", packageVersion("fastqcr"), sep = ": "),     
    paste("Biostrings", packageVersion("Biostrings"), sep = ": "),   

    "\n~ External programs ~",
    
    
    # Add it to the log file
    file=get.value("Raw.log"), 
    append = T, sep = "\n")


if(get_os() %in% c("os","linux")){ # to run only on os and not windows
  
  cat(paste("fastqc", system2("fastqc", "-v", stdout=T, stderr=T), sep = ": "),     
      paste("cutadapt", system2("cutadapt", "--version", stdout=T, stderr=T), sep = ": "),     
      #paste("vsearch", system2("vsearch", "-v", stdout=T, stderr=T)[1] , sep = ": "),   
      #paste("usearch", system2("usearch", "--version", stdout=T, stderr=T) , sep = ": "),      
      "\n-------------------------\n", 
            # Add it to the log file
      file=get.value("Raw.log"), 
      append = T, sep = "\n")


}

# END OF THE SCRIPT

seqtab.Pvit.int %>% str()
