
# Info --------------------------------------------------------------------

# Simple code to rename Raw files into something a little bit 
# easier to work with

# Audrey Bourret
# 2021-09-08

# Library -----------------------------------------------------------------

library(parallel)
library(here)
library(magrittr)
library(stringr)

# Parameters --------------------------------------------------------------------

numCores <- 1


# Create new files names --------------------------------------------------

# Get the current names of zipped files

old.names <- list.files(file.path(here::here(), "00_Data", "01a_RawData"),
                        full.name = T, 
                        pattern = ".fastq") %>% 
                        stringr::str_subset(".md5", negate = T) %>% 
                        stringr::str_subset("i1|i2", negate = T)

head(old.names)
length(old.names)

# Pattern to remove (new pattern can be added):

# MiSeq GenomeQuebec pattern
pat.MI <- paste0("MI.M", paste(rep("[:digit:]",5), collapse=""), "_", paste(rep("[:digit:]",4), collapse=""), ".", paste(rep("[:digit:]",3), collapse=""), ".FLD", paste(rep("[:digit:]",4), collapse=""), ".")
# NovaSeq GenomeQuebec pattern + index
pat.NS <- c(# Pattern in 2020
            paste0("NS.", paste(rep("[:digit:]",4), collapse=""), ".", paste(rep("[:digit:]",3), collapse=""), ".FLD", paste(rep("[:digit:]",4), collapse=""), ".", paste(rep("[:digit:]",4), collapse=""), "---PE1-CS1-IDT_i5_[:digit:]."),
            # Pattern in 2021
            paste0("NS.", paste(rep("[:digit:]",4), collapse=""), ".FLD", paste(rep("[:digit:]",4), collapse=""), "---PE1-CS1-IDT_i5_[:digit:].")
            )

pat.rm <- c(pat.MI, pat.NS)
pat.rm 

# Doit être TRUE partout, sinon changer le patron

name.test <- old.names %>% stringr::str_detect(paste(pat.rm, collapse = "|"))
table(name.test)


rename.raw <- function(file.name) {
         file.name %>% stringr::str_replace("01a_RawData",  "01b_RawData_rename") %>% 
                       stringr::str_remove_all(paste(pat.rm, collapse = "|")) 
  }

new.names <- rename.raw(old.names)
new.names 

# Change files names ------------------------------------------------------

# Parallel version, up to 6 times faster
#system.time(

parallel::mclapply(old.names,
        FUN = function(i) file.copy(i, 
                                    rename.raw(i)),
          mc.cores = numCores
          )  
#)
