quick.rdp <- function(fasta.file, 
                      out.path,
                      hier.path,
                      training.path){
  
  cmd <- paste("-Xmx40g",
               "classify" ,
               "-t", training.path,
               "-o", out.path,
               '-hier_outfile', hier.path, #tab-delimited output file containing the assignment count for each taxon in the hierarchical format.
               fasta.file
  )
  
  A <- system2("rdp_classifier", cmd, stdout = TRUE, stderr = TRUE)
  A
  
}


# A function

taxon.thr.RDP <- function(rdp.out, threshold = 80){
  
  names(rdp.out)  <- c("ESV", "nothing", "cellular", "Level1", "confidence_cellular", 
                       # "superkingdom", "Level2", "confidence_superkingdom",
                       "kingdom", "Level3", "confidence_kingdom", 
                       "phylum", "Level4", "confidence_phylum", 
                       "class", "Level5", "confidence_class", 	
                       "order", "Level6", "confidence_order",
                       "family", "Level7", "confidence_family",
                       "genus", "Level8", "confidence_genus",
                       "species", "Level9", "confidence_species")
  
  
  rdp.out$Taxon <- NA
  rdp.out$Levels <- NA
  
  # Go one ESV at the time
  for(x in 1:nrow(rdp.out)){
    
    # Reset at each row
    Taxon.int <- NULL
    Levels.int <- NULL
    
    for(y in c("confidence_species", "confidence_genus", "confidence_family", "confidence_order", "confidence_class", "confidence_phylum", "confidence_kingdom", #"confidence_superkingdom",
               "confidence_cellular")) {
      
      level.thr <- rdp.out[x,] %>% pull(y)
      
      if(level.thr >= threshold/100){
        Taxon.int <- rdp.out[x,] %>% pull(stringr::str_remove(y, "confidence_")) %>% stringr::str_replace_all("_", " ")
        Levels.int <-  stringr::str_remove(y, "confidence_")
        
        
        break;
      }
      
    }
    
    if(is.null(Taxon.int)){
      
      Taxon.int <- "Not assigned"
      Levels.int <- NA
      
    }
    
    rdp.out[x, "Taxon"] <- Taxon.int
    rdp.out[x, "Levels"] <- Levels.int
    
    
    
  }
  
  rdp.final <- rdp.out %>% dplyr::select(ESV, Taxon, Levels, 
                                         cellular, confidence_cellular, 
                                         # superkingdom, confidence_superkingdom,
                                         kingdom, confidence_kingdom, 
                                         phylum, confidence_phylum, 
                                         class, confidence_class, 	
                                         order, confidence_order,
                                         family, confidence_family,
                                         genus, confidence_genus,
                                         species, confidence_species)
  
  return(rdp.final)
  
}

