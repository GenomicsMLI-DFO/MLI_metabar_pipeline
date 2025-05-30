
# Info --------------------------------------------------------------------

# eDNA pipeline using DADA2
# Template pipeline
# Samples demultiplexed but loci multiplexed

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

# Internal functions
source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
source(file.path(here::here(), "01_Code", "Functions", "fastqc.R"))
source(file.path(here::here(), "01_Code", "Functions", "cutadapt.R"))
source(file.path(here::here(), "01_Code", "Functions", "dada2.R"))

# Add python env to this specific project
PythonENV.path <- get.value("PythonENV.path")
PythonENV.path

file.exists(PythonENV.path)

Sys.setenv(PATH = paste(c(PythonENV.path,
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))


# Check that fastqc, mutliqc and fastp are found on this computer

system2("cutadapt", "--help")
system2("multiqc", "--help")
system2("fastp", "--help")

# Data --------------------------------------------------------------------

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

# What is in the data.info
data.info %>% group_by(Loci) %>% summarise(N = n())
data.info %>% group_by(Run) %>% summarise(N = n())

# What is in the options
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

# 1. RAW quality assesment (fastp) ---------------------------------------

files.to.use <- list.files("00_Data/01b_RawData_rename",
                           pattern = "R1.fastq.gz",
                           full.names = T)
files.to.use %>% length()

# Progress bar
stepi <- 0
pb <- txtProgressBar(min = 0, style = 3, max = length(files.to.use), initial = stepi ) 


for(x in files.to.use){
  
  cat("\nProcessing:", x, "\n")
  
  stepi <- stepi + 1
  setTxtProgressBar(pb,stepi)
  
  #settings the files
  file1 <- x
  file2 <- file1 %>% stringr::str_replace("_R1", "_R2")
  
  file1.out <- file1 %>%  stringr::str_replace("01b_RawData_rename", "01c_RawData_Fastp") 
  
  file2.out <- file1.out %>% stringr::str_replace("_R1", "_R2")
  out.json <- file1 %>% str_replace("00_Data/01b_RawData_rename/", "02_Results/01_FastQC/01_Raw/" ) %>% 
    str_replace("_R1.fastq.gz", ".json") 
  
  out.html <- out.json %>% stringr::str_replace(".json", ".html")
  out.log  <- out.json %>% stringr::str_replace(".json", ".log")
  
  # The command
  
  if(file.exists(x)){
    
    cmd1 <- paste("-w", numCores, # Number of cores
                  "-i", file1,
                  "-I", file2,
                  "-o", file1.out,
                  "-O", file2.out,
                  "-j", out.json,
                  "-h", out.html,
                  sep = " ")
   cmd1
    
    A1 <- system2("fastp", cmd1, stdout=T, stderr=T) # to R console
    A1
    
    # save a file log
    cat(file = out.log,
        "Running Fastp", cmd1, A1,
        append= T, sep = "\n\n")
  
    gc()   

  }
}

close(pb) # Close the connection

# Internal function
multiqc.fastp("02_Results/01_fastp/")

# 2. RAW to FILT (cutadapt + dada2) ------------------------------------------------------------

# 2.1 CUTADAPT

# We need to remove the adaptors, and discard reads untrimmed
# this function can work with more than one loci
# The option novaseq TRUE/FALSE allowed to used to option -nextseq-trim=20

system2("cutadapt", "--help")

cutadapt.multiplex(folder.in = file.path(here::here(), "00_Data", "01c_RawData_Fastp"), 
                   folder.out = file.path(here::here(), "00_Data", "02a_Cutadapt"), 
                   loci = LOCUS, 
                   sens = SENS, 
                   numCores = numCores,
                   novaseq = FALSE) 

# Extract cutadapt res

cutadapt.res <- read_csv(file.path(here::here(), "00_Data", "02a_Cutadapt", "log", "Cutadapt_Stats.csv"))

# Copy in another folder
write_csv(cutadapt.res, file = file.path(here::here(), "02_Results", "02_Filtrations","Cutadapt_Stats_Nreads.csv"))

graph.cutadapt.1 <- cutadapt.res %>% pivot_longer(c(Raw, Adapt), names_to = "Step", values_to = "Nreads") %>%
  mutate(Nreads1 = Nreads + 1 ,
         Step = factor(Step, levels = c("Raw", "Adapt"))) %>% 
  left_join(data.info %>% select(ID_sample, Loci, ID_project, Sample_type), 
            by = c("ID_sample", "Loci")) %>% 
  ggplot(aes(x = Step, y = Nreads1, col = Sample_type, group = ID_sample)) +
  #geom_jitter(height = 0) +
  geom_point() +
  geom_line() + 
  #geom_boxplot() +
  scale_y_continuous(trans="log10") +
  labs(y = "N reads + 1 (log)", x = "Pipeline step")+ 
  facet_grid(~Loci) +
  theme_bw()+
  theme(legend.position = "bottom")
graph.cutadapt.1

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "Cutadapt_filt.png"), plot =  graph.cutadapt.1, height = 5, width = 6)


graph.cutadapt.2 <- cutadapt.res %>%  
  mutate(Perc = Adapt / Raw) %>% 
  left_join(data.info %>% select(ID_sample, Loci, ID_project, Sample_type), 
            by = c("ID_sample", "Loci")) %>% 
  ggplot(aes(x = ID_sample, y = Perc, fill = Loci)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.5, lty = "dashed")+
  geom_hline(yintercept = 1, col = "red")+
  #geom_jitter(height = 0) +
  labs(y = "% of read recovery", x = "Library")+ 
  facet_grid(~Sample_type, scale = "free", space = "free") +
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
graph.cutadapt.2

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "Cutadapt_prop_multiplex.png"), plot =  graph.cutadapt.2, height = 4, width = 10)




# Running another fastqc following cutadapt
fastqc(folder.in = file.path(here::here(), "00_Data", "02a_Cutadapt"),
       folder.out = file.path(here::here(), "02_Results", "01_FastQC", "02_Cutadapt"),
       numCores = numCores)

multiqc(folder.out = file.path(here::here(), "02_Results", "01_FastQC", "02_Cutadapt"),
        loci = LOCUS, 
        sens = SENS)


# 2.2 DADA2

# Check if these parameters seem good, if not change them
PARAM.DADA2 <- readr::read_tsv(file.path(here::here(), "01_Code/Parameters/dada2_param.tsv"))
PARAM.DADA2

dada2.filter (folder.in = file.path(here::here(), "00_Data", "02a_Cutadapt"), 
         folder.out = file.path(here::here(), "00_Data", "02b_Filtered_dada2"), 
         loci = LOCUS, 
         sens = SENS, 
         param.dada2 = PARAM.DADA2,
         numCores = numCores
)
         
# Extract dada2 summary
dada2.summary <- data.frame()

for(l in LOCUS){
  
  dada2.int <- readr::read_csv(file.path("00_Data", "02b_Filtered_dada2", "log", paste0(l, "_dada2_summary.csv"))) 
  dada2.int <- dada2.int %>% mutate(Loci = l)
  dada2.summary <- bind_rows(dada2.summary, dada2.int)
  
}

dada2.summary


dada2.summary <- dada2.summary %>% mutate(Dada2 = reads.out,
                                          Adapt = reads.in,
                                          Remove = reads.in - reads.out) %>% 
  select(-c(reads.in, reads.out)) %>% 
  tidyr::pivot_longer(cols = c(Dada2, Adapt), names_to = "Reads", values_to = "N")  

write_csv(dada2.summary, file = file.path(here::here(), "02_Results", "02_Filtrations", "Dada2Filt_Stats_Nreads.csv"))

graph.dada2 <- dada2.summary %>%  mutate(N1 = N + 1,
                                         #Reads = factor(Reads, levels = c("Remove", "Keep")),
                                         ID_sample = ID %>% str_remove(paste(paste0("_", LOCUS), collapse = "|"))) %>%
  left_join(data.info) %>% 
  ggplot(aes(x = Reads, y = N1, col = Sample_type, group = ID_sample)) +
  geom_point() +
  geom_line() + 
  #geom_boxplot() +
  scale_y_continuous(trans="log10") +
  labs(y = "N reads + 1 (log)", x = "Pipeline step")+ 
  facet_grid(~Loci) +
  theme_bw()+
  theme(legend.position = "bottom")

graph.dada2

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "DADA2_filt.png"), plot =  graph.dada2, height = 5, width = 6)

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
# Unusure if it's alway necessary to compute error rate by amplicon, but given that
# we run amplicon of various length, I think it's necessary

for(l in LOCUS){
  
  cat("\nCalculting error rate for" , l, "- F\n")
  
  # create files lists
  filesF.temp <- list.files(file.path(here::here(), "00_Data", "02b_Filtered_dada2"), full.name =T, pattern = ".fastq") %>%
    str_subset(paste0("_", SENS[1], "_")) %>% 
    str_subset(paste0("_",l,"_")) 
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "files were found\n")  
  
  #err.F.temp <- learnErrors(filesF.temp)
  assign(x = paste0("err.F.", l), value = learnErrors(filesF.temp,
                                                      nbases = 1e8, # 1e8 is the default - increase to sample more samples
                                                      randomize = T,
                                                      MAX_CONSIST = 10, # default - can be more if not reach
                                                      multithread = ifelse(numCores > 1, numCores, F)))
  
  #if(nrow(PARAM.temp == 2)){
  filesR.temp <-  list.files(file.path(here::here(), "00_Data", "02b_Filtered_dada2"), full.name =T, pattern = ".fastq") %>%
    str_subset(paste0("_", SENS[2], "_")) %>% 
    str_subset(paste0("_",l,"_")) 
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesR.temp), "files were found\n")  
  
  assign(x = paste0("err.R.", l), value = learnErrors(filesR.temp,
                                                      nbases = 1e8, # 1e8 is the default
                                                      randomize = T,
                                                      MAX_CONSIST = 10, # default - can be more if not reach
                                                      multithread = ifelse(numCores > 1, numCores, F)))
  
  #} else {err.R.temp <- vector()}
  
  #  cat("\nPlotting error rate results for" , l, "\n")  
  
  # Print a PDF of error rate
  
  pdf(file.path(here::here(), "00_Data", "03a_ErrorRate_dada2", paste0("ErrorsRate.", l, ".pdf"))) 
  suppressWarnings(print(plotErrors(get(paste0("err.F.",l)), nominalQ=TRUE)))
  suppressWarnings(print(plotErrors(get(paste0("err.R.",l)), nominalQ=TRUE)))
  dev.off()
  
  cat("\nSaving error rate results for" , l, "\n") 
  
  save(list = paste0("err.F.",l),
       file = file.path(here::here(), "00_Data", "03a_ErrorRate_dada2", paste0("err.F.",l,".Rdata")))
  
  #if(nrow(PARAM.temp == 2)){ 
  save(list = paste0("err.R.",l),
       file = file.path(here::here(), "00_Data", "03a_ErrorRate_dada2",paste0("err.R.",l,".Rdata")))
}  


# 3.2 Dereplication and sample inference

# Since we are working with BIG DATA, we don't have enough memory
# to run it by locus most of the time. So this code run it by sample directly.

# To reload error rate if necessary  
if(length(str_subset(ls(), "err.F.")) != length(LOCUS)){
  
  for(x in list.files(file.path(here::here(), "00_Data", "03a_ErrorRate_dada2"), full.names = T, pattern = ".Rdata")){
    load(x)
  }
  
}

for(l in LOCUS){
  
  cat("\nWorking on " , l, "\n")
  
  # create files lists
  filesF.temp <- list.files(file.path(here::here(), "00_Data", "02b_Filtered_dada2"), full.name =T, pattern = ".fastq") %>%
    str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
    str_subset(paste0("_", SENS[1], "_"))
  
  #filesR.temp <- filesF.temp %>% str_replace(SENS[1], SENS[2])
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "files were found\n")  
  
  mergers <- vector("list", length(filesF.temp))
  names(mergers) <- filesF.temp

  # Load overlap threshold
   minOverlap <- PARAM.DADA2 %>% dplyr::filter(Locus == l, Sens == "R1") %>% pull(minOverlap)  
  
  # Set a progress bar
  pb <- txtProgressBar(min = 0, max = length(filesF.temp), style = 3)
  
  for(i in seq_along(filesF.temp)){
    
    samF <- filesF.temp[i]
    samR <- samF %>% str_replace( paste0("_", SENS[1], "_"),  paste0("_", SENS[2], "_"))
    
    cat("\nProcessing:", samF ,"\n" )
    
    cat("Dereplication\n" )
    
    # Dereplication
    derep.F <-  derepFastq(samF)
    derep.R <-  derepFastq(samR)    
    
    cat("Sample inference\n" )
    
    #Sample inference
    dada.F <-  dada(derep.F, 
                    err =  get(paste0("err.F.",l)), 
                    multithread = ifelse(numCores > 1, numCores, F),
                    pool=TRUE)
    
    dada.R <-  dada(derep.R, 
                    err =  get(paste0("err.R.",l)), 
                    multithread = ifelse(numCores > 1, numCores, F),
                    pool=TRUE)
    
    cat("Merge paires\n" )
    
    merger <- mergePairs(dadaF = dada.F, 
                         derepF = derep.F, 
                         dadaR = dada.R, 
                         derepR = derep.R, 
                         minOverlap = minOverlap, 
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
       file = file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("seqtab.raw.", l, ".Rdata")  ))
  
  close(pb)
  
}
# 3.6 Length filter

# To reload seqtab if necessary  
if(length(str_subset(ls(), "seqtab.")) != length(LOCUS)){
  
  for(x in list.files(file.path(here::here(), "00_Data", "03b_SeqTab_dada2"), full.names = T, pattern = ".Rdata")){
    load(x)
  }
  
}


length.stat <- tibble()

for(l in LOCUS){
  
  cat("\nLength filtratrion for" , l, "\n")
  
  MINLEN <- PARAM.DADA2 %>% dplyr::filter(Locus == l, Sens == "R1") %>% pull(minESVLen)
  MAXLEN <- PARAM.DADA2 %>% dplyr::filter(Locus == l, Sens == "R1") %>% pull(maxESVLen)
  
  seqlens <- nchar(getSequences(get(paste0("seqtab.",l,".int"))))
  
  length.stat  <-  bind_rows(length.stat, tibble(Locus = l,ESVlength = seqlens))
  
  assign(x = paste0("seqtab.length.", l), value = get(paste0("seqtab.",l,".int"))[, seqlens >= MINLEN & seqlens <= MAXLEN]  )

  cat("\nKept",  ncol(get(paste0("seqtab.length.", l))), "ESV out of",   ncol(get(paste0("seqtab.",l,".int"))), "\n") 

  cat("\nSaving SeqTab after length filtration for" , l, "\n") 
  
  save(list = paste0("seqtab.length.", l),
       file = file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("seqtab.length.", l, ".Rdata")  ))
  
}

gg.prefilt <- length.stat  %>% ggplot(aes(x =ESVlength )) +
  geom_histogram() +
  geom_vline(data = PARAM.DADA2 %>% dplyr::filter(Locus %in% LOCUS, Sens == "R1"), aes(xintercept =minESVLen ), col = "darkred", lty = "dashed") +
  geom_vline(data = PARAM.DADA2 %>% dplyr::filter(Locus %in% LOCUS, Sens == "R1"), aes(xintercept =maxESVLen ), col = "darkred", lty = "dashed") +
  facet_wrap(~Locus, scale = "free", nrow = 1) +
  ggtitle("ESV filtration based on sequence length") +
  theme_bw(base_size = 8)

gg.prefilt

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "ESV_length_filt.png"), 
       plot =  gg.prefilt, height = 5, width = 6)

# 3.6 Remove chimera 

# To reload seqtab if necessary  
if(length(str_subset(ls(), "seqtab.")) != length(LOCUS)*2){
  
  for(x in list.files(file.path(here::here(), "00_Data", "03b_SeqTab_dada2"), full.names = T, pattern = ".Rdata")){
    load(x)
  }
  
}


for(l in LOCUS){
  
  cat("\nRemoving chimera for" , l, "\n")
  
  assign(x = paste0("ESVtab.", l), value = removeBimeraDenovo(get(paste0("seqtab.length.",l)), method = "consensus", 
                                                              multithread = ifelse(numCores > 1, numCores, F), verbose = TRUE)
  )
  
  cat("\nSaving SeqTab without chimera for" , l, "\n") 
  
  save(list = paste0("ESVtab.", l),
       file = file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("ESVtab.noCHIM.", l, ".Rdata")  ))
  
  
}

# Compute stats before and after chimera removal
# From the N reads perspective

Nread.summary <- data.frame(ID_sample = character(),
                            Loci = character(),
                            Merge = numeric(),
                            Final = numeric(),
                            stringsAsFactors = F)  

for(l in LOCUS){
  
  RES <- data.frame(ID_sample = rowSums(get(paste0("seqtab.",l,".int"))) %>% names(),
                    Loci = l,
                    Merge = rowSums(get(paste0("seqtab.",l,".int"))),
                    ESVlength =  rowSums(get(paste0("seqtab.length.",l))),
                    Final =  rowSums(get(paste0("ESVtab.",l))) ,
                    stringsAsFactors = F)  
  
  
  Nread.summary   <- bind_rows(Nread.summary, RES)
  
}  

Nread.summary

write_csv(Nread.summary, file = file.path(here::here(), "02_Results", "02_Filtrations", "ESVtab_Stats_Nreads.csv"))


graph.Nread <- Nread.summary %>% tidyr::pivot_longer(cols = c(Merge, ESVlength, Final), names_to = "Reads", values_to = "N")  %>%
  mutate(Reads = factor(Reads, levels = c("Merge", "ESVlength", "Final")),
         N1 = N + 1) %>%   
  left_join(data.info) %>% 
  ggplot(aes(x = Reads, y = N1, col = Sample_type, group = ID_sample)) +
  geom_point() +
  geom_line() + 
  #geom_boxplot() +
  scale_y_continuous(trans="log10") +
  labs(y = "N reads + 1 (log)", x = "Pipeline step")+ 
  facet_grid(~Loci) +
  theme_bw()+
  theme(legend.position = "bottom")

graph.Nread

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "Nreads_filt.png"), plot =  graph.Nread, height = 5, width = 6)

# From the N ESV perspective    

ESV.summary <- data.frame(ID_sample = character(),
                          Loci = character(),
                          Merge = numeric(),
                          ESVlength = numeric(),
                          Final = numeric(),
                          stringsAsFactors = F)  

for(l in LOCUS){
  
  RES <- data.frame(ID_sample = rowSums(get(paste0("seqtab.",l,".int"))) %>% names(),
                    Loci = l,
                    Merge = apply(get(paste0("seqtab.",l,".int")), MARGIN = 1, FUN = function(x){length(x[x>0])}),
                    ESVlength = apply(get(paste0("seqtab.length.",l)), MARGIN = 1, FUN = function(x){length(x[x>0])}),      
                    Final =  apply(get(paste0("ESVtab.",l)), MARGIN = 1, FUN = function(x){length(x[x>0])}) ,
                    stringsAsFactors = F)  
  
  
  ESV.summary   <- bind_rows(ESV.summary, RES)
  
}  

ESV.summary

write_csv(ESV.summary, file = file.path(here::here(), "02_Results", "02_Filtrations", "ESVtab_Stats_NESV.csv"))


graph.ESV <- ESV.summary %>% tidyr::pivot_longer(cols = c(Merge,ESVlength, Final), names_to = "Reads", values_to = "N")  %>%
  mutate(Reads = factor(Reads, levels = c("Merge", "ESVlength", "Final"))) %>%   
  left_join(data.info) %>% 
  ggplot(aes(x = Reads, y = N, col = Sample_type, group = ID_sample)) +
  geom_point() +
  geom_line() + 
  #geom_boxplot() +
  scale_y_continuous(trans="log10") +
  labs(y = "N ESV", x = "Pipeline step")+ 
  facet_grid(~Loci) +
  theme_bw()+
  theme(legend.position = "bottom")

graph.ESV

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "ESV_filt.png"), plot =  graph.ESV, height = 5, width = 6)


# Save ESV table and fasta file --------------------------------------


for(l in LOCUS){
  
  write.dada2.res(ESVtab = get(paste0("ESVtab.",l)), 
                  loci = l, 
                  folder = file.path(here::here(), "00_Data", "03c_ESV"))
  
}

# Compute ESV lengh -------------------------------------------------------

seq.df <- data.frame(Loci = character(),
                     width = numeric(),
                     Nreads = numeric())

for(l in LOCUS){
  

DNA.int <- readDNAStringSet(file.path(here::here(), "00_Data", "03c_ESV", paste0("ESV.",l,".fasta" )))

seq.df.int <-  data.frame(Loci = l,
                          width = DNA.int@ranges@width,
                          Nreads = get(paste0("ESVtab.",l))  %>% colSums() )


seq.df <- bind_rows(seq.df, seq.df.int)


}

gg.seq1 <- seq.df %>% ggplot(aes(x = width)) +
  geom_bar() +
  facet_wrap(~Loci, scale = "free", ncol= length(LOCUS)) + 
  labs(x = "ESV length (pb)", y = "N ESV") +
  theme_bw()
gg.seq1

gg.seq2 <- seq.df %>% ggplot(aes(x = width, y = Nreads)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Loci, scale = "free", ncol = length(LOCUS)) + 
  labs(x = "ESV length (pb)", y = "N reads") +
  theme_bw()
gg.seq2


gg.seq <- ggpubr::ggarrange(gg.seq1, gg.seq2,
                  nrow = 2,
                  labels = LETTERS,
                  align = "hv")

gg.seq

ggsave(filename = file.path(here::here(), "02_Results/02_Filtrations/", "ESV_length.png"), plot =  gg.seq, height = 5, width = 8)


# Write a final log

cat("\nEND of 02_Process_RAW.R script\n",
    paste("Pipeline version:", get.value("MLI.version")),
    date(),
    "\n-------------------------\n", 
    
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("dada2", packageVersion("dada2"), sep = ": "),     
    paste("Biostrings", packageVersion("Biostrings"), sep = ": "),   
    
    "\n~ External programs ~",
    paste("fastp", system2("fastqc", "-v", stdout=T, stderr=T), sep = ": "),       
    paste("fastqc", system2("fastqc", "-v", stdout=T, stderr=T), sep = ": "),     
    paste("multiqc", system2("multiqc", "--version", stdout=T, stderr=T), sep = ": "),     
    paste("cutadapt", system2("cutadapt", "--version", stdout=T, stderr=T), sep = ": "),  
    
    # Add it to the log file
    file = file.path(here::here(), "00_Data", "03c_ESV", "Process_RAW.log"), 
    append = T, sep = "\n")
