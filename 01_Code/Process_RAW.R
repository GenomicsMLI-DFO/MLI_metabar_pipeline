
# Info --------------------------------------------------------------------

# eDNA pipeline using DADA2
# Test rapide pour phoca 2020 - avec Taq multi
# 
# Audrey Bourret
# 2021-06-23
#

start.time <- Sys.time()

#

# Library -----------------------------------------------------------------

# Data manipulation

library(tidyverse) # includes ggplot2 dplyr and stringr

library(gtools)    # for mixedsort
library(readxl)

#library(devtools)
#devtools::install_github("kassambara/ggpubr")
#library(ggpubr)    # on github - for nice graphs

# Fastq and fasta manipulation


# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# BiocManager::install("dada2")

library(Biostrings)
library(dada2); packageVersion("dada2") # Faire mettre cette info dans le log

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}

# Add pythom env to this specific project
Sys.setenv(PATH = paste(c("/home/genleia/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

system2("cutadapt", "--help")
system2("multiqc", "--help")

auto.folder()


# Data --------------------------------------------------------------------

data.info <- read_csv(file.path(get.value("info.path"), "SeqInfo.csv") )
data.info

# FR.Fish <- data.info %>% filter(Locus == "Fish", ID_projet == "Frankeinfish") %>% 
#                          pull(ID_GQ)
# 
# FR.COI <- data.info %>% filter(Locus == "COI", ID_projet == "Frankeinfish") %>% 
#   pull(ID_GQ)
# 
# PP3.COI <- data.info %>% filter(Locus == "COI", ID_projet == "PP3") %>% 
#   pull(ID_GQ)
# 
# paste(PP3.COI, collapse = "|")
# https://benjjneb.github.io/dada2/bigdata_paired.html

LOCUS <- str_split(get.value("locus"), pattern = ";")[[1]]
SENS  <- str_split(get.value("sens"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analysed:", LOCUS, "\n", sep = " ")

numCores <- if(get_os() %in% c("os","linux")){
  detectCores() # Utilise le max de coeurs si sur linux
} else 1

cat("There are", numCores, "cores available on this computer", sep = " ")

# Setting up the log file

if(file.exists(get.value("Raw.log"))){

  switch(menu(title = "Do you want to erase the previous RAW log files?", graphics = T, 
            choice = c("yes", "no")
            )+1,
       # Answer 0
       cat("Nothing done\n"),
       # Answer 1 (yes)
      {cat("\nA new log file Process_RAW.log.txt was created (the previous RAW log was erased)\n\n")
       cat("\n-------------------------\n", 
           "Process raw eDNA data\n",
           date(),
           "\n-------------------------\n", 
           file=get.value("Raw.log"), 
           append = FALSE, sep = "\n")
           }, 
       # Answer 2 (no) 
      {cat ("\nInformation will be append to the log file Process_RAW.log.txt\n\n")
       cat("\n\n-------------------------\n",
            "NEW ANALYSIS PERFORMED", 
            date(),
            "\n-------------------------\n", 
            file=get.value("Raw.log"), 
            append = T, sep = "\n")
      }
        )

} else {
  cat ("\nThe log file Process_RAW.log.txt was created\n\n")  
  cat("\n-------------------------\n", 
      "Process raw eDNA data\n",
      date(),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = FALSE, sep = "\n")
}

# Functions ---------------------------------------------------------------

# Check that fastqc is found on this computer

#system2("fastqc", "--help")

# Wrapper around FastQC
fastqc <- function(folder.in, folder.out) {
  
  if(get_os() %in% c("os","linux")){ # to run only on os and not windows
    
    Nfiles <- length(list.files(folder.in))
    
    cat("Performing a FastQC analysis on", Nfiles,"files \n")
    
    file.remove(list.files(folder.out, full.name = T, pattern ="fastqc"))
    
    temp <-  mclapply(list.files(folder.in, full.name = T, pattern = "fastq"),
                      FUN = function(i){cmd <- paste("--outdir", folder.out, i, "-q")
                      system2("fastqc", cmd)
                      } ,
                      mc.cores = numCores
    )  
    
    # Save info on the log file
    cat("FastQC analysis was performed:",
        get.value("result.FQcutadapt.path"),
        "\n-------------------------\n", 
        file=get.value("Raw.log"), 
        append = T, sep = "\n") 
    
  } else {cat("Cannot perform FastQC on windows yet -- sorry!!")}
} # End of my function


# Creating a FastQC report with MultiQC
multiqc <- function(folder.out){
  for(l in LOCUS){
    print(l)
    for(s in SENS){
      print(s)
      cmd <- paste(list.files(folder.out, full.names = T) %>%
                     str_subset(paste0("_",l,"_")) %>% # update 2020-06-12 for FishAB vs Fish ...
                     str_subset(paste0("_",s)) %>%
                     str_subset(".zip"),
                   "--outdir", file.path(folder.out, "MultiQC_report"),
                   "--filename", paste0("multiqc_report_",l, "_", s, ".html"),
                   "-f" # to rewrite on previous data
      )

      system2("multiqc", cmd)

    }
  }
  
  
}



# 1. RAW quality assesment (fastqc) ---------------------------------------

# FastQC
get.value("raw_rename.path")
file.exists(get.value("raw_rename.path"))

get.value("result.FQraw.path")
file.exists(get.value("result.FQraw.path"))

fastqc(get.value("raw_rename.path"), get.value("result.FQraw.path"))

list.files(get.value("raw_rename.path"))[1:10]

# MultiQC for the others? Or by project rather then by loci?

multiqc(get.value("result.FQraw.path"))

# 2. RAW to FILT (cutadapt + dada2) ------------------------------------------------------------

# 2.1 CUTADAPT

list.files(get.value("filt_cutadapt.path"), full.name = T, pattern =".fastq") %>% 
  str_subset(pattern = "_12S_")


system2("cutadapt", "--help")

# CHECK THE PATTERN FOR LOCUS - NOW "_LOCUS_"
# NOW WITH NOVASEQ poly-G problem resolve with --nextseq-trim


if(get_os() %in% c("os","linux")){ 

  if(length(SENS) == 2){
    
    # Remove old files
    file.remove(list.files(get.value("filt_cutadapt.path"), full.name = T, pattern =".fastq"))
    file.remove(list.files(get.value("filt_cutadapt.log"), full.name = T, pattern ="_log.txt"))
    
    for(l in LOCUS){
    
      print(l)
      
      F.primer <- str_split(get.value(paste(l, "primers", sep=".")), pattern = ";")[[1]][1]
      R.primer <- str_split(get.value(paste(l, "primers", sep=".")), pattern = ";")[[1]][2] 
      
      F.primer <- F.primer %>% str_replace_all("I", "N")
      R.primer <- R.primer %>% str_replace_all("I", "N")
      
      mclapply(list.files(get.value("raw_rename.path"), pattern = paste0("_",l,"_"), full.names=T) %>% str_subset(paste0(SENS[1],".fastq")),
               FUN = function(file1){file2 <- file1 %>% str_replace(paste0(SENS[1],".fastq"), paste0(SENS[2],".fastq")) 
                                     new.file1 <- file1 %>% str_replace(get.value("raw_rename.path"), get.value("filt_cutadapt.path")) %>%
                                                            str_replace (".fastq", "_cut.fastq") 
                                     new.file2 <- new.file1 %>% str_replace(paste0(SENS[1],"_cut.fastq"), paste0(SENS[2],"_cut.fastq")) 
                 
                                     cmd <- paste("-g", paste0("^", F.primer),
                                                  "-G", paste0("^", R.primer),
                                                  "-o", new.file1, 
                                                  "--paired-output", new.file2, 
                                                  file1,
                                                  file2,
                                                  #"-f", "fastq",      
                                                  "--discard-untrimmed", 
                                                  #"--nextseq-trim=20", # To trim when quality go down just before poly-Q tail
                                                  #"--report=minimal",
                                                  #"-j", numCores, #Ncore
                                                  sep = " ") # forward adapter
                                     
                                     A <- system2("cutadapt", cmd, stdout=T, stderr=T) # to R console
                                     #system2("cutadapt", cmd)
                                     # TO UPDaTE
                                     
                                     # save a file log
                                     cat(file = new.file1 %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_cutadapt.log")) %>% 
                                           str_replace(".fastq.gz","_log.txt") %>% 
                                           str_remove(paste0("_",SENS[1])),
                                         A, # what to put in my file
                                         append= F, sep = "\n")
                                     
                                 },
               mc.cores = numCores
      )  
      
      
    }
    
    
  } else {print("No code for unpaired reads")}
  

} else {cat("Cannot perform cutadapt on windows yet -- sorry!!")}



# Creating a FastQC report with MultiQC
# THIS IS NOT WORKING WELL, we should add something here
#for(l in LOCUS){
#  print(l)
#  cmd <- paste(list.files(get.value("filt_cutadapt.log"), full.names = T) %>% 
#                 str_subset(l),
#               "--outdir", get.value("result.Cutadapt"),
#               #"--module", "cutadapt",
#               "--filename", paste0("multiqc_report_",l, ".html"),
#               "-f", # to rewrite on previous data
#               "-v -s -d"
#  )
#  
#  system2("multiqc", cmd)
#  
#}


# Running another fastqc following cutadapt
fastqc(get.value("filt_cutadapt.path"), get.value("result.FQcutadapt.path"))

multiqc(get.value("result.FQcutadapt.path"))


# 2.2 DADA2

PARAM.DADA2 <- expand.grid(Locus = LOCUS, Sens = SENS) %>% as.data.frame()
PARAM.DADA2$truncQ    <- 10
PARAM.DADA2$truncLen  <- 0
PARAM.DADA2$trimLeft  <- 0
PARAM.DADA2$trimRight <- 0
PARAM.DADA2$maxLen    <- Inf
PARAM.DADA2$minLen    <- 50
PARAM.DADA2$minQ      <- 0
PARAM.DADA2$maxEE     <- 1 # was 2 before

#PARAM.DADA2$truncLen[which(PARAM.DADA2$Locus == "12S")] <- 100
PARAM.DADA2$truncLen[which(PARAM.DADA2$Locus == "12S")] <- 150
#PARAM.DADA2$truncLen[which(PARAM.DADA2$Locus == "FishAB")] <- 160

PARAM.DADA2$minLen[which(PARAM.DADA2$Locus == "COI")]   <- 150
PARAM.DADA2$minLen[which(PARAM.DADA2$Locus == "Pvit")]   <- 150

PARAM.DADA2

# Remove old files
file.remove(list.files(get.value("filt_dada2.path"), full.name = T, pattern =".fastq"))

for(l in LOCUS){

cat("\nFiltering",l, "\n")  
  
PARAM.temp <- PARAM.DADA2 %>% filter(Locus == l)  

# create files lists
filesF.temp <- list.files(get.value("filt_cutadapt.path"), full.name =T, pattern = ".fastq") %>%
                str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
                str_subset(paste0(PARAM.temp$Sens[1] %>% as.character(),"_cut"))


# Add a message to be sure that the right number of files were used
cat(length(filesF.temp), "F files were found\n")  

if(nrow(PARAM.temp == 2)){
  filesR.temp <- list.files(get.value("filt_cutadapt.path"), full.name =T, pattern = ".fastq") %>% 
    str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
    str_subset(paste0(PARAM.temp$Sens[2] %>% as.character(),"_cut"))
  

  
} else {filesR.temp <- vector()}

# Add a message to be sure that the right number of files were used
cat(length(filesR.temp), "R files were found\n") 

filesF.filt.temp <- filesF.temp  %>% str_replace(get.value("filt_cutadapt.path"),get.value("filt_dada2.path"))
filesR.filt.temp <- filesR.temp  %>% str_replace(get.value("filt_cutadapt.path"),get.value("filt_dada2.path"))

filter.summary.temp <- filterAndTrim(fwd = filesF.temp ,
                                    filt = filesF.filt.temp,
                                    rev = filesR.temp,
                                    filt.rev = filesR.filt.temp,
                                    truncQ= PARAM.temp$truncQ, # minimum Q score, 10 = 90% base call accuracy
                                    truncLen = PARAM.temp$truncLen, # Taille min/max des reads
                                    trimLeft= PARAM.temp$trimLeft, # Deja enlevé avec cutadapt, sinon c(18,18) 
                                    trimRight= PARAM.temp$trimRight,
                                    maxLen = PARAM.temp$maxLen,
                                    minLen = PARAM.temp$minLen, # after trimming and truncation
                                    minQ = PARAM.temp$minQ,
                                    maxEE=PARAM.temp$maxEE,
                                    compress = TRUE, # Voir si c'est OK de compresser spour JAMP
                                    multithread=ifelse(numCores > 1, T, F), # TRUE on linux
                                    verbose = TRUE) 

cat(l, ":\n",
    #Data
    "Date: ", date(), "\n",
    # Dada2
    "Dada2 v", packageVersion("dada2") %>% as.character() , "\n",
    #N samples
    nrow(filter.summary.temp),
    " samples were process\n",

    #Prop samples
    round(sum(filter.summary.temp[,2])/sum(filter.summary.temp[,1]),3)*100,
    "% of reads remains after trimming (",
    sum(filter.summary.temp[,2]),
    " reads)\n",
    #Parameters:
    "Parameters: truncQ = ", paste(PARAM.temp$truncQ[1], PARAM.temp$truncQ[2], sep="-"), 
                ", truncLen = ", paste(PARAM.temp$truncLen[1], PARAM.temp$truncLen[2], sep="-"), # Taille min/max des reads
                ", trimLeft = ", paste(PARAM.temp$trimLeft[1], PARAM.temp$trimLeft[2], sep="-"), # Deja enlevé avec cutadapt, sinon c(18,18) 
                ", maxLen = ", paste(PARAM.temp$maxLen[1], PARAM.temp$maxLen[2], sep="-"),
                ", minLen = ", paste(PARAM.temp$minLen[1], PARAM.temp$minLen[2], sep="-"), # after trimming and truncation
                ", maxEE = ", paste(PARAM.temp$maxEE[1], PARAM.temp$maxEE[2], sep="-"),"\n",
    "\n",
    sep = "",
    file = file.path(get.value("filt_dada2.log"), paste0(l,"_dada2_filtering.log")),
    append = TRUE
    )

cat(readLines(file.path(get.value("filt_dada2.log"), paste0(l,"_dada2_filtering.log"))), sep = "\n")

}

cat("Quality filtering was performed with dada2:",
    get.value("filt_dada2.path"),
    get.value("dada2.filt.data"),
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 

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
