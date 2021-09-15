#system2("fastqc", "--help")

# Wrapper around FastQC
fastqc <- function(folder.in, folder.out) {
  
    Nfiles <- length(list.files(folder.in))
    
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

