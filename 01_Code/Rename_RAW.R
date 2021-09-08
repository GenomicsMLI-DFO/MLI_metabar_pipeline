
# Info --------------------------------------------------------------------

# Simple code to rename Raw files into something a little bit 
# easier to work with

# Special cases for run M03992 and M05812

# Audrey Bourret
# 2019-04-03

# Library -----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(parallel)

#BiocManager::install("Biostrings")

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}

# Data --------------------------------------------------------------------

numCores <- if(get_os() %in% c("os","linux")){
               detectCores() # Utilise le max de coeurs si sur linux
            } else 1

cat("There are", numCores, "cores available on this computer", sep = " ")
#info.path <- get.value("info.path")
#Sample.xl <- get.value("Sample.xl")

#raw_unz.path <- get.value("raw_unz.path")


# Create new files names --------------------------------------------------

# Get the current names of zipped files

# Pattern to remove :
list.files(get.value("raw.path"), full.name = T, pattern = ".fastq")[1:10]

pat.rm <- "MI.M[:digit:][:digit:][:digit:][:digit:][:digit:]_[:digit:][:digit:][:digit:][:digit:].[:digit:][:digit:][:digit:].FLD[:digit:][:digit:][:digit:][:digit:]."
pat.rm.extra <- "--PE1-CS1-IDT_i5_10."
# Doit être TRUE partout, sinon changer le patron
name.test <- old.names %>% str_detect(pat.rm)
name.test <- old.names %>% str_detect(pat.rm.extra)

table(name.test)

list.files(get.value("raw.path"), full.name = T, pattern = ".fastq")[!name.test]

old.names <- list.files(get.value("raw.path"), full.name = T, pattern = ".fastq") %>% str_subset(".md5", negate = T) %>% 
                                                                                      str_subset("i1|i2", negate = T)       
old.names

new.names <- old.names %>% 
                str_replace(get.value("raw.path"),  get.value("raw_rename.path")) %>% 
                str_remove(pat.rm) %>% 
                str_remove(pat.rm.extra) %>%       
                str_replace("MiFish", "12S") %>% 
                str_replace("dloop", "Pvit")

old.names[1:10]
new.names

# Change files names ------------------------------------------------------

# Remove old files

#file.remove(list.files(get.value("raw_rename.path"), full.name = T, pattern =".fastq"))

# Copier unz -> unz_rename (can take some times)

# system.time(
# file.copy(old.names,
#           new.names
#           )
# )

# Parallele version, 6 times faster with 96 cores
#system.time(
mclapply(old.names,
        FUN = function(i) file.copy(i, 
                                    i %>% 
                                      str_replace(get.value("raw.path"),  get.value("raw_rename.path")) %>% 
                                      str_remove(pat.rm)%>% 
                                      str_remove(pat.rm.extra) %>% 
                                      str_replace("MiFish", "12S") %>% 
                                      str_replace("dloop", "Pvit")),
          mc.cores = numCores
          )  
#)
