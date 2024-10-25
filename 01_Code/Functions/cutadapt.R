
cutadapt <- function(folder.in, folder.out, loci, sens, numCores, novaseq = FALSE) {
  
    if(length(sens) == 2){
      
    cat("\nWorking on both directions, F and R\n")  
      
      for(l in loci){
        
        cat("\nProcessing", l, "\n")  
        
        F.primer <- str_split(get.value(paste(l, "primers", sep=".")), pattern = ";")[[1]][1]
        R.primer <- str_split(get.value(paste(l, "primers", sep=".")), pattern = ";")[[1]][2] 
        
        F.primer <- F.primer %>% str_replace_all("I", "N")
        R.primer <- R.primer %>% str_replace_all("I", "N")
        
        parallel::mclapply(list.files(folder.in, pattern = paste0("_",l,"_"), full.names=T) %>% str_subset(paste0(SENS[1],".fastq")),
                 FUN = function(file1){file2 <- file1 %>% str_replace(paste0(SENS[1],".fastq"), paste0(SENS[2],".fastq")) 
                 new.file1 <- file1 %>% str_replace(folder.in, folder.out) %>%
                   str_replace (".fastq", "_cutadapt.fastq") 
                 new.file2 <- new.file1 %>% str_replace(paste0(SENS[1],"_cutadapt.fastq"), paste0(SENS[2],"_cutadapt.fastq")) 
                 
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
                 # Add nexttrim
                 if(novaseq == T){
                 cmd <- paste(cmd, 
                              "--nextseq-trim=20", # To trim when quality go down just before poly-Q tail

                              sep = " ")  
                   
                 }
                 
                 A <- system2("cutadapt", cmd, stdout=T, stderr=T) # to R console
                 #system2("cutadapt", cmd)
                 # TO UPDaTE
                 
                 # save a file log
                 cat(file = new.file1 %>% str_replace(folder.out, file.path(folder.out, "log")) %>% 
                       str_replace(".fastq.gz","_log.txt") %>% 
                       str_remove(paste0("_",SENS[1])),
                     A, # what to put in my file
                     append= F, sep = "\n")
                 
                 },
                 mc.cores = numCores
        )  
        
        
      }
      
      
    } else {print("No code for unpaired reads")}
    
  
  cat("\nTrimming is over, files were saved in", folder.out, "\n")
  
  # Get stats
  
  extract_cutadapt_stats(folder.in = file.path(folder.out, "log"), 
                         loci = loci)
    
}


extract_cutadapt_line <- function(TAB, TEXT){
  res <- TAB %>% str_subset(TEXT) %>% 
    str_remove(TEXT) %>% 
    str_remove_all(" ") %>%
    str_remove_all(",") #%>% 
  #as.numeric()
  return(res)
}


extract_cutadapt_stats <- function(folder.in, loci){
  # List log files
  cutadapt.log <- list.files(folder.in, full.names = F) %>% str_subset("log.txt")
  
   
  # New dataframe 
  RES <- data.frame(ID_sample = character(),
                    Loci = character(),
                    Raw = numeric(),
                    Adapt = numeric(),
                    stringsAsFactors = F)
  # Loop over the files
  
  for(x in seq_along(cutadapt.log)){
    
    cat(cutadapt.log[x], sep="\n")
 
    log.file <- readLines(file.path(folder.in, cutadapt.log[x]))
    
       
    sample_full <-  cutadapt.log[x] %>% str_remove("_cutadapt_log.txt")
    sample <- sample_full %>% stringr::str_remove(paste(paste0("_", loci), collapse = "|"))
    locus  <- sample_full %>% stringr::str_remove(paste0(sample, "_"))  
    N.start <- extract_cutadapt_line(log.file, "Total read pairs processed:" ) %>% as.numeric()
    N.cut <- sapply(str_split(extract_cutadapt_line(log.file, "Pairs written \\(passing filters\\):"), "\\("), `[`,1 ) %>% as.numeric()
    
    RES[x,] <- c(sample, locus, N.start, N.cut)
    
  
  
}

  Reads.sum <- RES %>% mutate(Raw = as.numeric(as.character(Raw)),
                              Adapt = as.numeric(as.character(Adapt)))
  
  write_csv(Reads.sum, file = file.path(folder.in, "Cutadapt_Stats.csv"))
  
  cat("\nStats can be found in", file.path(folder.in, "Cutadapt_Stats.csv"), "\n")
                              
}




# For multiplex data

cutadapt.multiplex <- function(folder.in, folder.out, loci, sens, numCores, novaseq = FALSE) {
  
  if(length(sens) == 2){
    
    cat("\nWorking on both directions, F and R\n")  
    
    for(l in loci){
      
      cat("\nProcessing", l, "\n")  
      
      F.primer <- str_split(get.value(paste(l, "primers", sep=".")), pattern = ";")[[1]][1]
      R.primer <- str_split(get.value(paste(l, "primers", sep=".")), pattern = ";")[[1]][2] 
      
      F.primer <- F.primer %>% str_replace_all("I", "N")
      R.primer <- R.primer %>% str_replace_all("I", "N")
      
      parallel::mclapply(list.files(folder.in, pattern = paste0("_multi_"), full.names=T) %>% str_subset(paste0(SENS[1],".fastq")),
                         FUN = function(file1){file2 <- file1 %>% str_replace(paste0(SENS[1],".fastq"), paste0(SENS[2],".fastq")) 
                         new.file1 <- file1 %>% str_replace(folder.in, folder.out) %>%
                           str_replace (".fastq", "_cutadapt.fastq") %>% 
                           str_replace ("_multi_", paste0("_", l, "_")) 
                         new.file2 <- new.file1 %>% str_replace(paste0(SENS[1],"_cutadapt.fastq"), paste0(SENS[2],"_cutadapt.fastq")) 
                         
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
                         # Add nexttrim
                         if(novaseq == T){
                           cmd <- paste(cmd, 
                                        "--nextseq-trim=20", # To trim when quality go down just before poly-Q tail
                                        
                                        sep = " ")  
                           
                         }
                         
                         A <- system2("cutadapt", cmd, stdout=T, stderr=T) # to R console
                         #system2("cutadapt", cmd)
                         # TO UPDaTE
                         
                         # save a file log
                         cat(file = new.file1 %>% str_replace(folder.out, file.path(folder.out, "log")) %>% 
                               str_replace(".fastq.gz","_log.txt") %>% 
                               str_remove(paste0("_",SENS[1])),
                             A, # what to put in my file
                             append= F, sep = "\n")
                         
                         },
                         mc.cores = numCores
      )  
      
      
    }
    
    
  } else {print("No code for unpaired reads")}
  
  
  cat("\nTrimming is over, files were saved in", folder.out, "\n")
  
  # Get stats
  
  extract_cutadapt_stats(folder.in = file.path(folder.out, "log"), 
                         loci = loci)
  
}


