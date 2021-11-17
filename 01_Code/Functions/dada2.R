
dada2.filter <- function(folder.in, folder.out, loci, sens, param.dada2, numCores) {
  
for(l in loci){
  
  cat("\nFiltering",l, "\n")  
  
  PARAM.temp <- param.dada2 %>% dplyr::filter(Locus == l)  
  
  # create files lists
  filesF.temp <- list.files(folder.in, full.name =T, pattern = ".fastq") %>%
    str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
    str_subset(PARAM.temp$Sens[1] %>% as.character())
  
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesF.temp), "F files were found\n")  
  
  if(nrow(PARAM.temp == 2)){
    filesR.temp <- list.files(folder.in, full.name =T, pattern = ".fastq") %>% 
      str_subset(paste0("_",l,"_")) %>% # Updated 2020-06-12 for FishAB vs Fish 
      str_subset(PARAM.temp$Sens[2])
    
    
    
  } else {filesR.temp <- vector()}
  
  # Add a message to be sure that the right number of files were used
  cat(length(filesR.temp), "R files were found\n") 
  
  filesF.filt.temp <- filesF.temp  %>% str_replace(folder.in, folder.out)
  filesR.filt.temp <- filesR.temp  %>% str_replace(folder.in, folder.out)
  
  filter.summary.temp <- dada2::filterAndTrim(fwd = filesF.temp ,
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
  
  filter.summary <- data.frame(ID = row.names(filter.summary.temp) %>% str_remove("_R1_cutadapt.fastq.gz"),
                               reads.in = filter.summary.temp[,1],
                               reads.out = filter.summary.temp[,2]
                               )
  
  readr::write_csv(filter.summary, file = file.path(folder.out, "log", paste0(l,"_dada2_summary.csv")))
  
  cat(l, ":\n",
      #Data
      "Date: ", date(), "\n",
      # Dada2
      "Dada2 v", packageVersion("dada2") %>% as.character() , "\n",
      #N samples
      nrow(filter.summary.temp),
      " samples were processed\n",
      
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
      file = file.path(folder.out, "log", paste0(l,"_dada2_filtering.log")),
      append = TRUE
  )
  
  cat(readLines(file.path(folder.out, "log", paste0(l,"_dada2_filtering.log"))), sep = "\n")
  
}

cat("Quality filtering was performed with dada2",
    "Files were saved in", folder.out, sep = "\n") 

}


write.dada2.res <- function(ESVtab, loci, folder){
  
  file1 <- file.path(folder, paste0("ESV.", loci, ".fasta"))
  file2 <- file1 %>% str_replace(".fasta", "_table.txt")
  
  DNA <- DNAStringSet(getSequences(ESVtab))
  names(DNA) <- paste0("ESV_", loci, "_", 1:length(DNA))
  
  writeXStringSet(DNA, file1)
  write.table(ESVtab, file2)
  
}
