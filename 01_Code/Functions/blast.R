library(stringr)
library(readr)

# Function to call blastn

quick.blastn <- function(fasta.file, out.file, 
                         db = "nt",
                         perc_identity = 95, 
                         qcov_hsp_perc = 95, 
                         max_target_seqs = 500, 
                         evalue = "1e-50",
                         NCBI.path = NCBI.path,
                         n.cores = numCores
                         ){
  
  cat("\nPerforming Blast over NCBI", db, "\n")
  cmd1 <- paste("-db", db, # echo $BLASTDB
                
                "-query",  fasta.file,
                "-evalue", evalue,
                "-qcov_hsp_perc", qcov_hsp_perc, 
                "-out", out.file, 
                "-perc_identity", perc_identity,
                "-num_threads", n.cores,
                "-max_target_seqs", max_target_seqs, 
                sep = " ")# forward adapter
  

  if(db %in% c("nt")){
  cmd1 <- paste(cmd1, 
        "-outfmt", "\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"", 
         sep = " ")  
  } else {
  cmd1 <-  paste(cmd1, 
          "-outfmt", "7", 
          sep = " ")           
         }
  
  
  A <-system2("blastn", cmd1, stdout=T, stderr=T,
              env = paste0("BLASTDB=", NCBI.path))
  
}


load.blast <- function(out.file, 
                       ncbi.tax = ncbi.tax,
                       db = "nt"){
  cat("\nLoading Blast results\n")
  
  RES.ncbi <- readr::read_tsv(out.file, comment = "#", col_names = F)
  
  if(db %in% c("nt")){
  
  names(RES.ncbi) <- c("QueryAccVer", "SubjectAccVer", "TaxoId","SciName", "SKindom", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
  RES.ncbi <- RES.ncbi %>% dplyr::left_join(ncbi.tax, by = c("TaxoId" = "id")) 
  } else{
    names(RES.ncbi) <- c("QueryAccVer", "SubjectAccVer",  "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
    RES.ncbi <- RES.ncbi %>% dplyr::left_join(ncbi.tax, by = c("SubjectAccVer" = "id")) 
  }

  return(RES.ncbi)
}

# Performed LCA at a given threshold, and return everything
BLAST_LCA <- function(RES, threshold = 0.97){
  DF <- data.frame()
  
  RES.OK <- RES %>% filter(#AlignmentLength >= .95 * width,
    Identity >= threshold) %>% 
    mutate(Taxon = NA,
           Levels = NA)
  
  # if from ncbi
  if(stringr::str_count(names(RES.OK), "SciName") %>% sum() == 1){
    RES.OK <- RES.OK %>% filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
    RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                            sapply(str_split(RES.OK$SciName, " "),`[`,2))
    RES.OK$genus <- sapply(str_split(RES.OK$species, " "),`[`,1)
    RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))
  }
  
  ASV <- RES.OK %>% pull(QueryAccVer) %>% unique()
  for(x in seq_along(ASV)){
    RES.INT <- RES.OK %>% filter(QueryAccVer == ASV[x])
    
    # loops around ranks
    for(y in c("species", "genus", "family", "order", "class", "phylum", "kingdom")) {
      
      N.LCA <- RES.INT %>% filter(!is.na(y)) %>%  pull(y) %>% unique() %>% str_subset("NA", negate = T) %>% length()
      
      if(N.LCA==1){
        RES.INT$Taxon <-  RES.INT %>% filter(!is.na(y)) %>% pull(y) %>% unique()  %>% str_subset("NA", negate = T)
        RES.INT$Levels <-  y
        break;
      }else{
        RES.INT[,y] <- NA
      }
    } # END of the loop around columns
    #RES.INT <- RES.INT %>% select(QueryAccVer, kingdom, phylum, class, order, family, genus, species, Taxon, Levels) %>% 
    #                     distinct(.keep_all = T)
    
    if(nrow(RES.INT[which(!is.na(RES.INT[,y])),])>0){
      DF <- bind_rows(DF, RES.INT[which(!is.na(RES.INT[,y])),])      
    }
    
    
  }
  return(DF)  
  
}   

BLAST_TOPHIT <- function(RES, threshold = 0.95){
  DF <- data.frame()
  
  RES.OK <- RES %>% filter(#AlignmentLength >= .95 * width,
    Identity >= threshold) %>% 
    mutate(Taxon = NA,
           Levels = NA)
  
  # if from ncbi
  if(str_count(names(RES.OK), "SciName") %>% sum() == 1){
    RES.OK <- RES.OK %>% filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
    RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                            sapply(str_split(RES.OK$SciName, " "),`[`,2))
    RES.OK$genus <- sapply(str_split(RES.OK$species, " "),`[`,1)
    RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))
  }
  
  ASV <- RES.OK %>% pull(QueryAccVer) %>% unique()
  for(x in seq_along(ASV)){
    RES.INT <- RES.OK %>% filter(QueryAccVer == ASV[x])
    
    evalue.min <- min(RES.INT$evalue)
    RES.INT <- RES.INT %>% filter(evalue == evalue.min)
    
    #identity.max <- max(RES.INT$Identity)
    #RES.INT <- RES.INT %>% filter(Identity == identity.max)
    
    # loops around ranks
    for(y in c("species", "genus", "family", "order", "class", "phylum", "kingdom")) {
      
      N.LCA <- RES.INT %>% filter(!is.na(y)) %>%  pull(y) %>% unique() %>% str_subset("NA", negate = T) %>% length()
      
      if(N.LCA==1){
        RES.INT$Taxon <-  RES.INT %>% filter(!is.na(y)) %>% pull(y) %>% unique()  %>% str_subset("NA", negate = T)
        RES.INT$Levels <-  y
        break;
      }else{
        RES.INT[,y] <- NA
      }
    } # END of the loop around columns
    #RES.INT <- RES.INT %>% select(QueryAccVer, kingdom, phylum, class, order, family, genus, species, Taxon, Levels) %>% 
    #                     distinct(.keep_all = T)
    
    if(nrow(RES.INT[which(!is.na(RES.INT[,y])),])>0){
      DF <- bind_rows(DF, RES.INT[which(!is.na(RES.INT[,y])),])      
    }
    
  }
  return(DF)  
  
}   

sum.BLAST <- function(DF){
  RES <- DF %>% select(QueryAccVer, Taxon, Levels, species, genus, family, order, class, phylum, kingdom) %>% unique() 
  
  return(RES)
}

# Extract blast parameters
get.blast.value <- function(Loci, value, df = PARAM.BLAST){
  df %>% dplyr::filter(Locus == Loci) %>%  dplyr::pull(value)
  
}
