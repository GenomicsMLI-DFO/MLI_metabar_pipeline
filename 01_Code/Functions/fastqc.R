#system2("fastqc", "--help")

# Wrapper around FastQC
fastqc <- function(folder.in, folder.out, numCores = 1) {
  
    Nfiles <- length(list.files(folder.in, pattern = ".fastq.gz"))
    
    cat("Performing a FastQC analysis on", Nfiles,"files \n")
    
    #file.remove(list.files(folder.out, full.name = T, pattern ="fastqc"))
    
    temp <-  parallel::mclapply(list.files(folder.in, full.name = T, pattern = "fastq"),
                      FUN = function(i){cmd <- paste("--outdir", folder.out, i, "-q")
                      system2("fastqc", cmd)
                      } ,
                      mc.cores = numCores
    )  
    
   # # Save info on the log file
   #  cat("FastQC analysis was performed:",
   #     get.value("result.FQcutadapt.path"),
  #      "\n-------------------------\n", 
   #     file=get.value("Raw.log"), 
    #    append = T, sep = "\n") 
    
    cat("\nFastQC is over, reports were saved in", folder.out, "\n")  
    
    
} # End of my function


# Creating a FastQC report with MultiQC
multiqc <- function(folder.out, loci, sens){
  for(l in loci){
    print(l)
    for(s in sens){
      print(s)
      cmd <- paste(paste(list.files(folder.out, full.names = T) %>%
                     str_subset(paste0("_",l,"_")) %>% # update 2020-06-12 for FishAB vs Fish ...
                     str_subset(paste0("_",s)) %>%
                     str_subset(".zip"), collapse = " "),
                   "--outdir", file.path(folder.out, "MultiQC_report"),
                   "--filename", paste0("multiqc_report_",l, "_", s, ".html"),
                   "-f" # to rewrite on previous data
      )
      
      system2("multiqc", cmd)
      
    }
  }
  
  cat(paste("\nMultiQC is over, reports were saved in", paste0(folder.out,"/MultiQC_report")),
      
      "\nYou should check:",
      "1. Read length", "2. Quality drops", "3. Adaptor content", sep = "\n"
      )  
  
}


# Creating a FastQC report with MultiQC
# but not by loci
multiqc.multiplex <- function(folder.out, sens){
  #for(l in loci){
  #  print(l)
    for(s in sens){
      print(s)
      cmd <- paste(paste(list.files(folder.out, full.names = T) %>%
                     #str_subset(paste0("_",l,"_")) %>% # update 2020-06-12 for FishAB vs Fish ...
                     str_subset(paste0("_",s)) %>%
                     str_subset(".zip"), collapse = " "),
                   "--outdir", file.path(folder.out, "MultiQC_report"),
                   "--filename", paste0("multiqc_report_", s, ".html"),
                   "-f" # to rewrite on previous data
      )
      
      system2("multiqc", cmd)
      
   # }
  }
  
  cat(paste("\nMultiQC is over, reports were saved in", paste0(folder.out,"/MultiQC_report")),
      
      "\nYou should check:",
      "1. Read length", "2. Quality drops", "3. Adaptor content", sep = "\n"
  )  
  
}


