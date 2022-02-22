
# Info --------------------------------------------------------------------

#
# Simple code to add TAGS 
#
# Audrey Bourret
# 2021-06-29
#

# Library -----------------------------------------------------------------

library(here)
#library(readxl)
library(tidyverse)

# MetaData --------------------------------------------------------------------

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info


# Index files (TAGS)
index_i7 <- read.delim(here::here("00_Data/00_FileInfos/GQ_fluidigm_index.txt"), sep = " ")
index_i7 %>% head()

index_i5 <- read.delim(here::here("00_Data/00_FileInfos/GQ_fluidigm_dual_index.txt"), sep = " ") %>%
  mutate(Dual_Mix = as.character(Dual_Mix))
index_i5 %>% head()


# Joining -----------------------------------------------------------------

data.info.modif <- data.info %>%  mutate(
    Key_i7 = sapply(str_split(File, "[.]"), `[`, 4),
    Key_i5 = sapply(str_split(File, "_"), `[`, 3)
  ) %>%
  left_join(index_i7 %>% select(Key_i7 = Barcode_Name, tag_fwd = Barcode_Sequence)) %>%
  left_join(index_i5 %>% select(Key_i5 = Dual_Mix, tag_rev = Direction)) %>% 
  mutate(tag_rev = ifelse(is.na(tag_rev), "AAAAAAAA", tag_rev))

data.info.modif %>% select(File, Key_i5, Key_i7, tag_fwd, tag_rev)


# Saving ------------------------------------------------------------------

readr::write_csv(data.info.modif, file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv"))

