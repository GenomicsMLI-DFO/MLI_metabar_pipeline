
# Info --------------------------------------------------------------------

# Simple code to rename Raw files into something a little bit 
# easier to work with
# multiplex files need _multi_ pattern instead of locus name

# Audrey Bourret
# 2022-04-14

# Library -----------------------------------------------------------------

library(parallel)
library(here)
library(magrittr)
library(stringr)

# Parameters --------------------------------------------------------------------

numCores <- 1


# Create new files names --------------------------------------------------

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

# Get the current names of zipped files

old.names <- list.files("./00_Data/01a_RawData",
                        full.name = T, 
                        pattern = ".fastq") %>% 
                        stringr::str_subset(".md5", negate = T) %>% 
                        stringr::str_subset("i1|i2", negate = T)

head(old.names)
length(old.names)

R1.old <- paste0("./00_Data/01a_RawData/", data.info$File, "_R1.fastq.gz")
R2.old <- paste0("./00_Data/01a_RawData/", data.info$File, "_R2.fastq.gz")

# There should be as many times duplicated values as multiplesex loci 
R1.old %>% length()
R1.old %>% unique() %>% length()

# Check that all is detected and TRUE
old.names %>% str_detect(paste(paste(unique(R1.old), collapse = "|"), paste(unique(R2.old), collapse = "|"), sep = "|")) %>% table()

R1.new <- paste0("./00_Data/01b_RawData_rename/", data.info$ID_labo, "_multi_R1.fastq.gz")  %>% unique()
R2.new <- paste0("./00_Data/01b_RawData_rename/",data.info$ID_labo, "_multi_R2.fastq.gz")  %>% unique()

R1.new %>% head()

length(R1.old) == length(R1.new)
 
# Change files names ------------------------------------------------------

for(i in seq_along(R1.old)){
  file.copy(from = R1.old, to = R1.new)
  file.copy(from = R2.old, to = R2.new)  
}

