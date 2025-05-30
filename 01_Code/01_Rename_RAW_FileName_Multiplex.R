
# Info --------------------------------------------------------------------

# Simple code to rename Raw files into something a little bit 
# easier to work with
# multiplex files need _multi_ pattern instead of locus name

# Library -----------------------------------------------------------------

library(parallel)
library(here)
library(magrittr)
library(stringr)

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

# There should be as many times duplicated values as multiplexed loci 
R1.old %>% length()
R1.old %>% unique() %>% length()

R2.old %>% length()
R2.old %>% unique() %>% length()

# Keep only unique R1 and R2 old value

R1.old <- R1.old %>% unique()
R2.old <- R2.old %>% unique()

# Check that all is detected and TRUE
old.names %>% str_detect(paste(paste(unique(R1.old), collapse = "|"), paste(unique(R2.old), collapse = "|"), sep = "|")) %>% table()

# check which ones are FALSE if any
old.names %>% str_subset(paste(paste(unique(R1.old), collapse = "|"), paste(unique(R2.old), collapse = "|"), sep = "|"),negate=TRUE) %>% table()

R1.new <- paste0("./00_Data/01b_RawData_rename/", data.info$ID_sample, "_multi_R1.fastq.gz")  %>% unique()
R2.new <- paste0("./00_Data/01b_RawData_rename/",data.info$ID_sample, "_multi_R2.fastq.gz")  %>% unique()

R1.new %>% head()

# The old and new files should be the same length

length(R1.old) == length(R1.new)
length(R2.old) == length(R2.new)
 
# Change files names ------------------------------------------------------

for(i in seq_along(R1.old)){
  file.copy(from = R1.old[i], to = R1.new[i])
  file.copy(from = R2.old[i], to = R2.new[i])  
}

## It could be faster to simply rename the files instead of copying them :) 
#for(i in seq_along(R1.old)){
#  file.rename(from = R1.old[i], to = R1.new[i])
#  file.rename(from = R2.old[i], to = R2.new[i])  
#}
