
# Info --------------------------------------------------------------------

# Corrections for negative control using metabar package
# https://metabarfactory.github.io/metabaR/
# Template pipeline
# 
# Audrey Bourret
# 2022-02-21

# Library -----------------------------------------------------------------

library("metabaR")
library(magrittr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))
source(file.path(here::here(), "01_Code", "Functions", "metabar.R"))

# Dataset -----------------------------------------------------------------

LOCUS <- stringr::str_split(get.value("Loci"), pattern = ";")[[1]]

cat(length(LOCUS), "loci will be analyse (", paste(LOCUS, collapse = ", "),")\nTheses parameters can be changed with the file Option.txt", sep = " ")

SUBGROUP <- stringr::str_split(get.value("group.metabaR"), pattern = ";")[[1]]

cat(length(SUBGROUP), "subgroup will be considered(", paste(SUBGROUP, collapse = ", "),")\nTheses parameters can be changed with the file Option.txt", sep = " ")

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

# Create a list of which PCR to be considered in each SUBGROUP

SUBGROUP.ls <- list()

for(x in SUBGROUP){
  
  subgroup <- c(data.info %>% dplyr::filter(ID_subprojet == x) %>% pull(ID_subprojet) %>% unique(),
                data.info %>% dplyr::filter(ID_subprojet == x) %>% pull(ID_projet) %>% unique(),                                                         
                "ALL", "All", "all") %>% unique()
  
  id <- data.info %>% dplyr::filter(ID_subprojet %in% subgroup)  %>% pull(ID_labo) %>% unique()
  
  SUBGROUP.ls[[x]] <- id
  
  
}

# Data conversion for metabaR ---------------------------------------------

# Copied from: https://metabarfactory.github.io/metabaR/articles/metabaRF-vignette.html?fbclid=IwAR04iNKTcIkKtsoOY66NOsWiAU2WSkopyqAy6ifUO5ibaOaVLfPfZp-3OxQ
#
# The basic data format used in metabaR is a metabarlist, a list of four tables:
#   
# - reads a table of class matrix consisting of PCRs as rows, and molecular operational taxonomic units (MOTUs) as columns. The number of reads for each MOTU in each PCR is given in each cell, with 0 corresponding to no reads.
# 
# - motus a table of class data.frame where MOTUs are listed as rows, and their attributes as columns. A mandatory field in this table is the field “sequence”, i.e. the DNA sequence representative of the MOTU. Examples of other attributes that can be included in this table are the MOTU taxonomic information, the taxonomic assignment scores, MOTU sequence GC content, MOTU total abundance in the dataset, etc.
# 
# - pcrs a table of class data.frame consisting of PCRs as rows, and PCR attributes as columns. This table is particularly important in metabaR, as it contains all the information related to the extraction, PCR and sequencing protocols and design that are necessary to assess and improve the quality of metabarcoding data (Taberlet et al. 2018; Zinger et al. 2019). This table can also include information relating to the PCR design, such as the well coordinates and plate of each PCR, the tag combinations specific of the PCR, the primers used, etc. Mandatory fields are:
#   sample_id: a vector indicating the biological sample origin of each PCR (e.g. the sample name)
# type : the type of PCR, either a sample or an experimental control amplification. Only two values allowed: "sample" or "control".
# control_type : the type of control. Only five values are possible in this field:
#   NA if type="sample", i.e. for any PCR obtained from a biological sample.
# "extraction" for DNA extraction negative controls, i.e. PCR amplification of an extraction where the DNA template was replaced by extraction buffer.
# "pcr" for PCR negative controls, i.e. pcr amplification where the DNA template was replaced by PCR buffer or sterile water.
# "sequencing" for sequencing negative controls, i.e. unused tag/library index combinations.
# "positive" for DNA extraction or PCR positive controls, i.e. pcr amplifications of known biological samples or DNA template (e.g. a mock community).
# 
# - samples a table of class data.frame consisting of biological samples as rows, and associated information as columns. Such information includes e.g. geographic coordinates, abiotic parameters, experimental treatment, etc. This table does not include information on the DNA metabarcoding experimental controls, which can only be found in pcrs


# Prepare metabar format

# Reads 
for(l in LOCUS){
  
  load(file.path(here::here(), "00_Data", "03b_SeqTab_dada2",paste0("ESVtab.noCHIM.", l, ".Rdata")  ))
  
  
  assign(x = paste0("DNA.", l), 
         value =Biostrings::readDNAStringSet(file.path(here::here(), "00_Data", "03c_ESV",paste0("ESV.", l, ".fasta")))
  )
  
}

tidy.ESV <- function(ESVtab, DNA.seq) {
  DNA.tidy <- tibble(ESV = names(DNA.seq), SEQ =  DNA.seq %>% as.character())
  
  ESV.tidy <- ESVtab %>% as_tibble() %>% 
    dplyr::mutate(ID_labo = row.names(ESVtab)) %>%
    tidyr::pivot_longer(cols = !ID_labo, names_to = "SEQ", values_to = "Nreads") %>% 
    dplyr::left_join(DNA.tidy)
  
  return(ESV.tidy) 
}


for(l in LOCUS){
  
  reads.int <- tidy.ESV( get(paste0("ESVtab.",l)),   get(paste0("DNA.",l))) %>% 
    dplyr::select(ID_labo, ESV, Nreads) %>% tidyr::pivot_wider(names_from = ESV, values_from = Nreads)
  
  reads <- as.matrix(reads.int %>% select(-ID_labo))
  dimnames(reads)[[1]] <- reads.int %>% pull(ID_labo)
  
  assign(x = paste0("reads.", l), 
         value = reads)
  
}

# Motus

load(file.path(here::here(), "02_Results/03_TaxoAssign/01_Blast/", "ESVtab_assign.Rdata"))

for(l in LOCUS){
  
  motus.int <- data.frame(sequence =  as.vector( get(paste0("DNA.",l))),
                          QueryAccVer = names(get(paste0("DNA.",l)))) 
  
  motus <-   motus.int %>% left_join(RES.all.ncbi %>% filter(Loci == l, Method == "TOP", Threshold == 95) %>% distinct(QueryAccVer, .keep_all = T) )
  row.names(motus) <- names(get(paste0("DNA.",l)))
  
  assign(x = paste0("motus.", l), 
         value = motus)
  
}

# PCRs

data.info

for(l in LOCUS){
  
  pcr.int <- data.info  %>%  dplyr::filter(Loci == l) %>% 
    dplyr::rename(sample_id = ID_labo,
                  plate_no = ID_plaque,
                  project = ID_subprojet) %>% 
    dplyr::mutate (type = ifelse(Type_echantillon %in% c("Echantillon", "ECH"), "sample", "control"),
                   control_type = ifelse(type == "sample", NA,
                                         ifelse(Type_echantillon %in% c("Neg_PCR", "PNC"),"pcr",
                                                ifelse(Type_echantillon %in% c("NTC"), "sequencing",
                                                       ifelse(Type_echantillon %in% c("PPC"), "positive", 
                                                              "extraction")))),
                   plate_row= str_sub(ID_puit, 1, 1),
                   plate_col= str_sub(ID_puit, 2,3) ,
                   plate_col = as.numeric(as.character(plate_col)),
                   #tag_fwd = Sequence_i7,
                   #tag_rev = Sequence_i5,
                   primer_fwd = sapply(stringr::str_split(get.value(paste0(l,".primers")), pattern = ";"), `[`, 1),
                   primer_rev = sapply(stringr::str_split(get.value(paste0(l,".primers")), pattern = ";"), `[`, 2),
    ) %>% as.data.frame()
  
  pcr.int
  row.names(pcr.int) <-  pcr.int$sample_id
  
  assign(x = paste0("pcr.", l), 
         value = pcr.int)
  
}


# samples

for(l in LOCUS){
  
  samples.int <- data.frame(sample_id =  data.info %>% dplyr::filter(Loci == l) %>% pull(ID_labo), info = NA) 
  row.names(samples.int )<- samples.int$sample_id
  
  assign(x = paste0("samples.", l), 
         value = samples.int)
  
}


## Create metabar objects

for(l in LOCUS){
  
  metabarlist.int <- metabarlist_generator(reads = get(paste0("reads.", l)), 
                                           motus = get(paste0("motus.", l)), 
                                           pcrs = get(paste0("pcr.", l)), 
                                           samples = get(paste0("samples.", l)))
  
  
  assign(x = paste0("metabarlist.ori.", l), 
         value = metabarlist.int )
  
}

# Load metabar threshold

metabar.param <- readr::read_tsv(file = file.path(here::here(), "01_Code/Parameters/metabar_param.tsv"))

metabar.param

# Tag jump ----------------------------------------------------------------

# Threshold tag can be defined in the file:  "01_Code/Parameters/metabar_param.tsv"
# If you change the parameters, just rerun this part

metabar.param <- readr::read_tsv(file = file.path(here::here(), "01_Code/Parameters/metabar_param.tsv"))
metabar.param %>% filter(Locus %in% LOCUS) %>%  select(Locus, tag.threshold)

# Define a vector of thresholds to test
thresholds.tag.test <- c(0, 0.0001, 0.001, 0.003, 0.005, 0.01, 0.03, 0.05) 

# LOOP over all LOCI
for(l in LOCUS){
  
  cat("\nLooking at different tag jump threshold for", l, "\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  
  # Compute basic stats and saved the resuls
  metabarlist.int$pcrs$nb_reads <- rowSums(metabarlist.int$reads)
  metabarlist.int$pcrs$nb_motus <- rowSums(metabarlist.int$reads>0)
  metabarlist.int$motus$count   <- colSums(metabarlist.int$reads)
  
  assign(x = paste0("metabarlist.ori.", l), 
         value = metabarlist.int )
  
  # Load threshold
  thresholds.tag <-  metabar.param %>% filter(Locus == l) %>% pull(tag.threshold)
  
  cat("Current threshold :", thresholds.tag , "\n")  
  
  #remove empty MOTUs
  #metabarlist.int.clean <- subset_metabarlist(metabarlist.int, "reads", 
  #                                             indices = (rowSums(metabarlist.int$reads)>0))
  
  # Run the tests and stores the results in a list
  tests.tagjump <- lapply(thresholds.tag.test, function(x) tagjumpslayer(metabarlist.int, x))
  # method = "cut" vs. method = "substract"
  # tests2020_substract <-lapply(thresholds, function(x) tagjumpslayer(Fish2020_clean,x, method = "substract"))
  
  # NOTE: Method = "substract" ne changerait rien ici parce que l'on utilise des présences / absences pour les différents MOTU. Si au dessus du seuil, ça ne changerait rien d'enlever des reads.
  
  
  names(tests.tagjump) <- paste("t_", thresholds.tag.test, sep="")
  
  # Format the data for ggplot with amount of reads at each threshold
  tests.tagjump.long <- reshape2::melt(as.matrix(do.call("rbind", lapply(tests.tagjump, function(x) rowSums(x$reads)))))
  colnames(tests.tagjump.long) <- c("threshold", "sample", "abundance")
  
  # Add richness in MOTUs at each threshold
  tests.tagjump.long$richness <-
    reshape2::melt(as.matrix(do.call("rbind", lapply(tests.tagjump, function(x) {
      rowSums(x$reads > 0)
    }))))$value
  
  # Add control type information on pcrs and make data curation threshold numeric
  tests.tagjump.long$controls <- metabarlist.int$pcrs$control[match(tests.tagjump.long$sample, rownames(metabarlist.int$pcrs))]
  tests.tagjump.long$threshold <- as.numeric(gsub("t_", "", tests.tagjump.long$threshold))
  
  # New table formatting for ggplot
  tests.tagjump.long.2 <- reshape2::melt(tests.tagjump.long, id.vars=colnames(tests.tagjump.long)[-grep("abundance|richness", colnames(tests.tagjump.long))])
  
  tag.gg <- tests.tagjump.long.2 %>% mutate(controls = ifelse(is.na(controls), "sample", controls)) %>% 
    ggplot(aes(x=as.factor(threshold), y=value + 1)) + 
    geom_jitter(aes(color=controls), height = 0, alpha=0.5) + 
    geom_boxplot(color="grey40", fill = "white", alpha = 0) + 
    geom_vline(xintercept =  factor(thresholds.tag), col="orange", lty=2) + 
    scale_color_manual(values = c("brown", "red", "cyan4","pink","green", "black","yellow","purple"), na.value = "darkgrey") +
    facet_grid(variable ~ controls, scale="free") + 
    theme_bw() + 
    scale_y_log10() +
    labs(x="MOTU pcr : total abundance filtering threshold", y="MOTUs + 1",
         title = paste("Validating tag jumping threshold for", l)) + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          axis.text.x = element_text(angle=40, h=1), 
          legend.position = "none")
  
  print(tag.gg)
  
  n.sample <- tag.gg$data$controls %>% unique() %>% length()
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("00_tagjump.threshold_",l, ".png")), 
         plot = tag.gg ,
         width = 2*n.sample,
         height = 4,
         units = c("in"))
  
  # Put in an object to be able to export it to the automatic report
  assign(x = paste0("tag.gg.", l), 
         value = tag.gg)
  
  # Check
  cat("Threshold MUST be validated with the graph ", paste0("02_Results/04_ESVtable_correction/00_tagjump.threshold_",l, ".png") , "\n")  
  
  metabarlist.int.clean <- tagjumpslayer(metabarlist.int, thresholds.tag )
  
  metabarlist.int.clean$pcrs$nb_reads.tagjump <- rowSums(metabarlist.int.clean$reads)
  metabarlist.int.clean$pcrs$nb_motus.tagjump <- rowSums(metabarlist.int.clean$reads>0)
  metabarlist.int.clean$motus$count.tagjump   <- colSums(metabarlist.int.clean$reads)
  
  assign(x = paste0("metabarlist.tagclean.", l), 
         value = metabarlist.int.clean)
  #
  # Save a summary to compute stat later on
  
  summary.int <- metabarlist.int.clean$pcrs %>% select(sample_id, nb_reads, nb_motus, nb_reads.tagjump, nb_motus.tagjump) %>% 
    mutate(Loci = l) 
  
  if(!file.exists(file.path(here::here(), "00_Data", "04_ESVcorrected"))){dir.create(file.path(here::here(), "00_Data", "04_ESVcorrected"))}  
  readr::write_csv(summary.int, file.path(here::here(), "00_Data", "04_ESVcorrected", paste0("ESVtab_Stats_postTagjump_", l, ".csv")))
  cat("Summary stats were saved here:", paste0("00_Data/04_ESVcorrected/ESVtab_Stats_postTagjump_",l, ".csv") , "\n")  
  
  # identify occurrence of the most abundant OTU
  idx <- which.max(metabarlist.int$motus$count)
  p1 <- ggpcrplate.modif(metabarlist.int,
                         legend_title = "# reads",
                         FUN = function(m) {
                           m$reads[, idx]
                         }
  )
  
  # same on clean data
  p2 <- ggpcrplate.modif(metabarlist.int.clean,
                         legend_title = "# reads",
                         FUN = function(m) {
                           m$reads[, idx]
                         }
  )
  
  plate.tag.gg <- ggpubr::ggarrange(p1 + scale_size(limits = c(1, max(metabarlist.int$reads[, idx]))) +
                                      ggtitle("Most abundant MOTU ori"),
                                    p2 + scale_size(limits = c(1, max(metabarlist.int$reads[, idx]))) +
                                      ggtitle("Most abundant MOTU after tagjumpslayer"),
                                    nrow = 1, ncol = 2 , common.legend = T, legend = "right")
  
  
  n.plate <- p1$data$plate_no %>% unique() %>% length()
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("00_tagjump.plate_",l, ".png")), 
         plot = plate.tag.gg ,
         width = 8,
         height = 2.5 * n.plate,
         units = c("in"), bg = "white")
  
  
} 

# Save important graph to an R object (to export it to the automatic report)
save(list = ls(pattern = "tag.gg."),
     file = file.path(here::here(), "02_Results/04_ESVtable_correction", "00_tag.gg.Rdata") )

# Compute N read By plate -------------------------------------------------------------

# Print stats and add a few count columns
for(l in LOCUS){
  
  cat("\nStatistic for", l, "original\n")  
  metabarlist.int       <- get(paste0("metabarlist.ori.",l))
  print(summary_metabarlist(metabarlist.int))
  
  cat("\nStatistic for", l, "after tag jump correction \n")  
  metabarlist.int.clean <- get(paste0("metabarlist.tagclean.",l))
  print(summary_metabarlist(metabarlist.int.clean))
  
  # Save image
  
  plate.ori.gg <- ggpcrplate.modif(metabarlist.int, legend_title = "N reads")
  
  plate.clean.gg <- ggpcrplate.modif(metabarlist.int.clean, legend_title = "N reads")
  
  plate.gg <- ggpubr::ggarrange(plate.ori.gg + scale_size(limits = c(1, max(metabarlist.int$reads[, ]))) +
                                  ggtitle("N reads observed ori"),
                                plate.clean.gg + scale_size(limits = c(1, max(metabarlist.int$reads[, ]))) +
                                  ggtitle("N reads observed after tagjumpslayer"),
                                nrow = 1, ncol = 2 , common.legend = T, legend = "right"
  )
  
  n.plate <- plate.ori.gg$data$plate_no %>% unique() %>% length()
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("01_plate_Nreads_",l, ".png")), 
         plot = plate.gg ,
         width = 8,
         height = 2.5 * n.plate,
         units = c("in"), bg = "white")
  
  #ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("00_plate_Nreads_",l, ".png")), 
  #       plot = plate.gg + ggtitle(paste("N reads by well for", l)),
  #       width = 5,
  #       height = 3 * n.plate,
  #       units = c("in"))
  
  cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/01_plate_Nreads_",l,".png"), "\n")
  
}

# Flag low read depth outlier ---------------------------------------------

for(l in LOCUS){
  
  cat("\nSequencing depth for", l, "\n")  
  
  depth.threshold <- metabar.param %>% filter(Locus == l) %>% pull(depth.threshold)
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  
  metabarlist.int.clean  <- get(paste0("metabarlist.tagclean.",l))
  
  depth.ori.gg <-  metabarlist.int$pcrs %>%  mutate(control_type = ifelse(is.na(control_type), "sample", control_type)) %>% 
    ggplot(aes(x = nb_reads+1, fill = control_type)) +
    geom_histogram(bins=40, color="grey") + 
    geom_vline(xintercept = depth.threshold, lty=2, color="orange") + # threshold
    scale_x_log10() + 
    scale_fill_manual(breaks = c("sample", "positive", "extraction", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
    labs(x="N reads + 1 (log)", 
         y="N samples") +
    labs(title = paste("Sequencing depth for", l),
         subtitle = "original" ) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_text(size = 8),
          title = element_text(size = 10),
          legend.title = element_blank())
  
  depth.gg.ori.project <-  metabarlist.int$pcrs %>%  mutate(control_type = ifelse(is.na(control_type), "sample",
                                                                                  ifelse(control_type == "extraction", "field/lab",      
                                                                                         control_type))) %>% 
    ggplot(aes(x = nb_reads+1, fill = control_type)) +
    geom_histogram(bins=40, color="grey") + 
    geom_vline(xintercept = depth.threshold, lty=2, color="orange") + # threshold
    scale_x_log10() + 
    scale_fill_manual(breaks = c("sample", "positive", "field/lab", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
    labs(x="N reads + 1 (log)", 
         y="N samples") +
    labs(title = paste("Sequencing depth for", l),
         subtitle = "original" ) +
    theme_bw() + 
    facet_grid(project ~ ., scale = "free_y") +
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_text(size = 8),
          title = element_text(size = 10),
          strip.text.y = element_text(angle = 0) ,
          legend.title = element_blank())
  
  #ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("01_depth_",l,".png")), 
  #       plot = depth.gg,
  #       width = 5,
  #       height = 3,
  #       units = c("in"))
  
  #ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("01_depth_byProject_",l,".png")), 
  #       plot = depth.gg.project,
  #       width = 6,
  #       height = 5,
  #      units = c("in"))
  
  depth.clean.gg <-  metabarlist.int.clean$pcrs %>%  mutate(control_type = ifelse(is.na(control_type), "sample", control_type)) %>% 
    ggplot(aes(x = nb_reads+1, fill = control_type)) +
    geom_histogram(bins=40, color="grey") + 
    geom_vline(xintercept = depth.threshold, lty=2, color="orange") + # threshold
    scale_x_log10() + 
    scale_fill_manual(breaks = c("sample", "positive", "extraction", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
    labs(x="N reads + 1 (log)", 
         y="N samples") +
    labs(title = paste("Sequencing depth for", l),
         subtitle = "after tagjumpslayer" ) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_text(size = 8),
          title = element_text(size = 10),
          legend.title = element_blank())
  
  depth.gg.clean.project <-  metabarlist.int.clean$pcrs %>%  mutate(control_type = ifelse(is.na(control_type), "sample",
                                                                                          ifelse(control_type == "extraction", "field/lab",      
                                                                                                 control_type))) %>% 
    ggplot(aes(x = nb_reads+1, fill = control_type)) +
    geom_histogram(bins=40, color="grey") + 
    geom_vline(xintercept = depth.threshold, lty=2, color="orange") + # threshold
    scale_x_log10() + 
    scale_fill_manual(breaks = c("sample", "positive", "field/lab", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
    labs(x="N reads + 1 (log)", 
         y="N samples") +
    labs(title = paste("Sequencing depth for", l),
         subtitle = "after tagjumpslayer" ) +
    theme_bw() + 
    facet_grid(project ~ ., scale = "free_y") +
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_text(size = 8),
          title = element_text(size = 10),
          strip.text.y = element_text(angle = 0) ,
          legend.title = element_blank())
  
  depth.gg <- ggpubr::ggarrange(depth.ori.gg,
                                depth.clean.gg,
                                nrow = 1, ncol = 2 , common.legend = T, legend = "bottom"
  )
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("02_depth_",l,".png")), 
         plot = depth.gg,
         width = 8,
         height = 3,
         units = c("in"),
         bg = "white")
  
  depth.gg.project <- ggpubr::ggarrange(depth.gg.ori.project,
                                        depth.gg.clean.project,
                                        nrow = 1, ncol = 2 , common.legend = T, legend = "bottom"
  )
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("02_depth_byProject_",l,".png")), 
         plot = depth.gg.project,
         width = 10,
         height = 8,
         units = c("in"),
         bg = "white")
  
  #ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("01_depth_byProject_",l,".png")), 
  #       plot = depth.gg.project,
  #       width = 6,
  #       height = 5,
  #      units = c("in"))
  
  cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/02_depth_",l,".png"), "\n")
  
  # Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
  metabarlist.int$pcrs$seqdepth_ok <- ifelse(metabarlist.int$pcrs$nb_reads < depth.threshold, F, T)
  metabarlist.int.clean$pcrs$seqdepth_ok <- ifelse(metabarlist.int.clean$pcrs$nb_reads.tagjump < depth.threshold, F, T)
  
  assign(x = paste0("metabarlist.ori.", l), 
         value = metabarlist.int )
  
  assign(x = paste0("metabarlist.tagclean.", l), 
         value = metabarlist.int.clean )
  
}


# Flag contaminants -------------------------------------------------------

# Will be performed both by project and overall, to have an idea of the difference 
# between both

for(l in LOCUS){
  
  cat("\nLooking at contaminants for", l, "\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  metabarlist.int.clean <- get(paste0("metabarlist.tagclean.",l))
  
  
  for(x in SUBGROUP){
    
    cat("Running on", x, "subgroup\n")  
    
    id.int <-  SUBGROUP.ls[[x]]
    
    metabarlist.int.sub <- subset_metabarlist(metabarlist.int, table="pcrs",
                                              indices = metabarlist.int$pcrs$sample_id %in% id.int)
    
    metabarlist.int.clean.sub <- subset_metabarlist(metabarlist.int.clean, table="pcrs",
                                                    indices = metabarlist.int.clean$pcrs$sample_id %in% id.int)
    
    # Recompute basic stats by subgroup
    metabarlist.int.sub$pcrs$nb_reads.subgroup <- rowSums(metabarlist.int.sub$reads)
    metabarlist.int.sub$pcrs$nb_motus.subgroup <- rowSums(metabarlist.int.sub$reads>0)
    metabarlist.int.sub$motus$count.subgroup   <- colSums(metabarlist.int.sub$reads)
    
    metabarlist.int.clean.sub$pcrs$nb_reads.tagjump.subgroup <- rowSums(metabarlist.int.clean.sub$reads)
    metabarlist.int.clean.sub$pcrs$nb_motus.tagjump.subgroup <- rowSums(metabarlist.int.clean.sub$reads>0)
    metabarlist.int.clean.sub$motus$count.tagjump.subgroup   <- colSums(metabarlist.int.clean.sub$reads)
    
    # Contaslayer
    
    metabarlist.int.sub  <- contaslayer(metabarlist.int.sub , method="max",
                                        control_types = c("pcr", "extraction"),
                                        output_col = "not_a_max_conta")
    
    metabarlist.int.clean.sub  <- contaslayer(metabarlist.int.clean.sub , method="max",
                                              control_types = c("pcr", "extraction"),
                                              output_col = "not_a_max_conta")
    
    common_contam.ori.max <- metabarlist.int.sub$motus %>% filter(not_a_max_conta == F) %>% select(Taxon, Levels, count = count.subgroup) %>% mutate(method = "sub.ori.max") %>% arrange(desc(count))
    common_contam.ori.max$ESV <- row.names(common_contam.ori.max)
    
    common_contam.clean.max <- metabarlist.int.clean.sub$motus %>% filter(not_a_max_conta == F) %>% select(Taxon, Levels, count = count.tagjump.subgroup) %>% mutate(method = "sub.tagclean.max") %>% arrange(desc(count))
    common_contam.clean.max$ESV <- row.names(common_contam.clean.max)
    
    common_contam <- bind_rows(common_contam.ori.max, common_contam.clean.max) %>% 
      pivot_wider(names_from = method, values_from = count)
    
    readr::write_csv(common_contam, 
                     file = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("03_conta_",l,"_",x, ".csv")))
    
    cat("The list of identified contaminant is here:", paste0("02_Results/04_ESVtable_correction/03_conta_",l,"_",x, ".csv"), "\n")
    
    assign(x = paste0("metabarlist.ori.", l, ".", x), 
           value = metabarlist.int.sub )
    
    assign(x = paste0("metabarlist.tagclean.", l, ".", x), 
           value = metabarlist.int.clean.sub)
    
    
    
  }  
  
  cat("Running on overall dataset\n")  
  
  # Run contaslayer on the overall dataset
  
  metabarlist.int  <- contaslayer(metabarlist.int, method = "max",
                                  control_types = c("pcr", "extraction"),
                                  output_col = "not_a_max_conta")
  
  metabarlist.int.clean  <- contaslayer(metabarlist.int.clean, method = "max",
                                        control_types = c("pcr", "extraction"),
                                        output_col = "not_a_max_conta")
  
  common_contam.ori.max <- metabarlist.int$motus %>% filter(not_a_max_conta == F) %>% select(Taxon, Levels, count = count) %>% mutate(method = "ori.max") %>% arrange(desc(count))
  
  common_contam.ori.max$ESV <- row.names(common_contam.ori.max)
  
  common_contam.clean.max <- metabarlist.int.clean$motus %>% filter(not_a_max_conta == F) %>% select(Taxon, Levels, count =count.tagjump) %>% mutate(method = "tagclean.max") %>% arrange(desc(count))
  
  common_contam.clean.max$ESV <- row.names(common_contam.clean.max)
  
  common_contam <- bind_rows( common_contam.ori.max, common_contam.clean.max) %>% 
    pivot_wider(names_from = method, values_from = count)
  
  readr::write_csv(common_contam, 
                   file = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("03_conta_",l, "_Overall.csv")))
  
  cat("The list of identified contaminant is here:", paste0("02_Results/04_ESVtable_correction/03_conta_",l,"_Overvall.csv"), "\n")
  
  assign(x = paste0("metabarlist.ori.", l), 
         value = metabarlist.int )
  
  assign(x = paste0("metabarlist.tagclean.", l), 
         value = metabarlist.int.clean)
  
# NE pas faire les figures si il n'y a pas assez de contaminants
if( (nrow(metabarlist.int$motus %>% filter(not_a_max_conta == F) ) & nrow(metabarlist.int.clean$motus %>% filter(not_a_max_conta == F) ))> 2) {
  
  conta.ori.gg <- ggpcrplate.cont(metabarlist.int, N = Inf)

  conta.clean.gg <- ggpcrplate.cont(metabarlist.int.clean, N = Inf)

  n.plate <- conta.ori.gg $data$plate_no %>% unique() %>% length()

  conta.gg <- ggpubr::ggarrange(conta.ori.gg,
                                conta.clean.gg,
                                nrow = 1, ncol = 2 , common.legend = T, legend = "right"
  )

}
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("03_plate_conta_",l,".png")), 
         plot = conta.gg,
         width = 8,
         height = 1.5 * n.plate,
         units = c("in"),
         bg = "white")
  
  
}

# Flag high % contaminant samples -------------------------------------------------

for(l in LOCUS){
  
  cat("\nLooking at samples with high contaminants for", l, "\n")  
  
  threshold.prop.cont <- metabar.param %>% filter(Locus == l) %>% pull(prop.cont.threshold)
  
  cat("Running on overall dataset\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  metabarlist.int.clean <- get(paste0("metabarlist.tagclean.",l))
  
  if(#is.null(nrow(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"])) |
    is.null(nrow(metabarlist.int$reads[,metabarlist.int$motus$not_a_max_conta=="FALSE"]))) {     
    
    Rel.conta.prop.ori <- data.frame(#conta.all.prop = metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"] / rowSums(metabarlist.int$reads),
      conta.max.prop = metabarlist.int$reads[,metabarlist.int$motus$not_a_max_conta=="FALSE"] / rowSums(metabarlist.int$reads))
    
  } else {
    Rel.conta.prop.ori <- data.frame(#conta.all.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"]) / rowSums(metabarlist.int$reads),
      conta.max.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_max_conta=="FALSE"]) / rowSums(metabarlist.int$reads))
  }
  
  if(#is.null(nrow(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"])) |
    is.null(nrow(metabarlist.int.clean$reads[,metabarlist.int.clean$motus$not_a_max_conta=="FALSE"])))  {     
    
    Rel.conta.prop.clean <- data.frame(#conta.all.prop = metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"] / rowSums(metabarlist.int$reads),
      conta.max.prop = metabarlist.int.clean$reads[,metabarlist.int.clean$motus$not_a_max_conta=="FALSE"] / rowSums(metabarlist.int.clean$reads))
    
  } else {
    Rel.conta.prop.clean <- data.frame(#conta.all.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"]) / rowSums(metabarlist.int$reads),
      conta.max.prop = rowSums(metabarlist.int.clean$reads[,metabarlist.int.clean$motus$not_a_max_conta=="FALSE"]) / rowSums(metabarlist.int.clean$reads))
  }
  
  
  # Add information on control types
  Rel.conta.prop.ori$control_type <- metabarlist.int$pcrs$control_type[match(rownames(Rel.conta.prop.ori), rownames(metabarlist.int$pcrs))]
  # Add information on control types
  Rel.conta.prop.clean$control_type <- metabarlist.int.clean$pcrs$control_type[match(rownames(Rel.conta.prop.clean), rownames(metabarlist.int.clean$pcrs))]
  
  conta.gg.ori <-  Rel.conta.prop.ori %>% pivot_longer(cols = c(conta.max.prop)) %>% 
    mutate(control_type = ifelse(is.na(control_type), "sample", control_type)) %>% 
    ggplot(aes(x=control_type, y=value, color=control_type)) + 
    geom_boxplot() + geom_jitter(alpha=0.5) +
    geom_hline(yintercept = threshold.prop.cont, lty = "dashed", color="orange") +
    scale_color_manual(breaks = c("sample", "positive", "extraction", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
    labs(x=NULL, y="Prop. reads flagged as contaminant",
         title = paste("Prop. of flagged MOTUs for", l),
         subtitle = "original") + 
    facet_grid(.~name) +
    theme_bw() +
    theme(axis.title = element_text(size = 8),
          title = element_text(size = 10))
  
  conta.gg.clean <-  Rel.conta.prop.clean %>% pivot_longer(cols = c(conta.max.prop)) %>% 
    mutate(control_type = ifelse(is.na(control_type), "sample", control_type)) %>% 
    ggplot(aes(x=control_type, y=value, color=control_type)) + 
    geom_boxplot() + geom_jitter(alpha=0.5) +
    geom_hline(yintercept = threshold.prop.cont, lty = "dashed", color="orange") +
    scale_color_manual(breaks = c("sample", "positive", "extraction", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
    labs(x=NULL, y="Prop. reads flagged as contaminant",
         title = paste("Prop. of flagged MOTUs for", l),
         subtitle = "after tagjumpslayer" ) + 
    facet_grid(.~name) +
    theme_bw() +
    theme(axis.title = element_text(size = 8),
          title = element_text(size = 10))
  
  conta.gg <- ggpubr::ggarrange(conta.gg.ori,
                                conta.gg.clean,
                                nrow = 1, ncol = 2 , common.legend = T, legend = "bottom"
  )
  
  ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("04_conta.prop_",l,"_ALL.png")), 
         plot = conta.gg ,
         width = 7,
         height = 3,
         bg = "white",
         units = c("in"))
  
  cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/04_conta.prop_",l,"_ALL.png"), "\n")
  
  metabarlist.int$pcrs$low_max_conta_level_10 <- ifelse(Rel.conta.prop.ori$conta.max.prop[match(rownames(metabarlist.int$pcrs), rownames(Rel.conta.prop.ori))]>0.1,  F, T)
  metabarlist.int.clean$pcrs$low_max_conta_level_10 <- ifelse(Rel.conta.prop.clean$conta.max.prop[match(rownames(metabarlist.int.clean$pcrs), rownames(Rel.conta.prop.clean))]>0.1,  F, T)
  
  metabarlist.int$pcrs$low_max_conta_level_10[metabarlist.int$pcrs$nb_reads == 0] <- FALSE
  metabarlist.int.clean$pcrs$low_max_conta_level_10[metabarlist.int.clean$pcrs$nb_reads.tagjump == 0] <- FALSE
  
  assign(x = paste0("metabarlist.ori.", l), 
         value = metabarlist.int )
  
  assign(x = paste0("metabarlist.tagclean.", l), 
         value = metabarlist.int.clean )
  
  
  for(x in SUBGROUP){
    
    cat("Running on", x, "subgroup\n")  
    
    metabarlist.int.sub <- get(paste0("metabarlist.ori.",l, ".", x))
    metabarlist.int.clean.sub <- get(paste0("metabarlist.tagclean.",l, ".", x))
    
    if(#is.null(nrow(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"])) |
      is.null(nrow(metabarlist.int.sub$reads[,metabarlist.int.sub$motus$not_a_max_conta=="FALSE"]))) {     
      
      Rel.conta.prop.ori.sub <- data.frame(#conta.all.prop = metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"] / rowSums(metabarlist.int$reads),
        conta.max.prop = metabarlist.int.sub$reads[,metabarlist.int.sub$motus$not_a_max_conta=="FALSE"] / rowSums(metabarlist.int.sub$reads))
      
    } else {
      Rel.conta.prop.ori.sub <- data.frame(#conta.all.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"]) / rowSums(metabarlist.int$reads),
        conta.max.prop = rowSums(metabarlist.int.sub$reads[,metabarlist.int.sub$motus$not_a_max_conta=="FALSE"]) / rowSums(metabarlist.int.sub$reads))
    }
    
    if(#is.null(nrow(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"])) |
      is.null(nrow(metabarlist.int.clean.sub$reads[,metabarlist.int.clean.sub$motus$not_a_max_conta=="FALSE"])))  {     
      
      Rel.conta.prop.clean.sub <- data.frame(#conta.all.prop = metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"] / rowSums(metabarlist.int$reads),
        conta.max.prop = metabarlist.int.clean.sub$reads[,metabarlist.int.clean.sub$motus$not_a_max_conta=="FALSE"] / rowSums(metabarlist.int.clean$reads))
      
    } else {
      Rel.conta.prop.clean.sub <- data.frame(#conta.all.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"]) / rowSums(metabarlist.int$reads),
        conta.max.prop = rowSums(metabarlist.int.clean.sub$reads[,metabarlist.int.clean.sub$motus$not_a_max_conta=="FALSE"]) / rowSums(metabarlist.int.clean.sub$reads))
    }
    
    # Add information on control types
    Rel.conta.prop.ori.sub$control_type <- metabarlist.int.sub$pcrs$control_type[match(rownames(Rel.conta.prop.ori.sub), rownames(metabarlist.int.sub$pcrs))]
    # Add information on control types
    Rel.conta.prop.clean.sub$control_type <- metabarlist.int.clean.sub$pcrs$control_type[match(rownames(Rel.conta.prop.clean.sub), rownames(metabarlist.int.clean.sub$pcrs))]
    
    metabarlist.int.sub$pcrs$low_max_conta_level_10 <- ifelse(Rel.conta.prop.ori.sub$conta.max.prop[match(rownames(metabarlist.int.sub$pcrs), rownames(Rel.conta.prop.ori.sub))]>0.1,  F, T)
    metabarlist.int.clean.sub$pcrs$low_max_conta_level_10 <- ifelse(Rel.conta.prop.clean.sub$conta.max.prop[match(rownames(metabarlist.int.clean.sub$pcrs), rownames(Rel.conta.prop.clean.sub))]>0.1,  F, T)
    
    
    metabarlist.int.sub$pcrs$low_max_conta_level_10[metabarlist.int.sub$pcrs$nb_reads.subgroup == 0] <- FALSE
    metabarlist.int.clean.sub$pcrs$low_max_conta_level_10[metabarlist.int.clean.sub$pcrs$nb_reads.tagjump.subgroup == 0] <- FALSE
    
    assign(x = paste0("metabarlist.ori.", l, ".", x), 
           value = metabarlist.int.sub )
    
    assign(x = paste0("metabarlist.tagclean.", l, ".", x), 
           value = metabarlist.int.clean.sub)
    
    conta.gg.ori.sub <-  Rel.conta.prop.ori.sub %>% pivot_longer(cols = c(conta.max.prop)) %>% 
      mutate(control_type = ifelse(is.na(control_type), "sample", control_type)) %>% 
      ggplot(aes(x=control_type, y=value, color=control_type)) + 
      geom_boxplot() + geom_jitter(alpha=0.5) +
      geom_hline(yintercept = threshold.prop.cont, lty = "dashed", color="orange") +
      scale_color_manual(breaks = c("sample", "positive", "extraction", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
      labs(x=NULL, y="Prop. reads flagged as contaminant",
           title = paste("Prop. of flagged MOTUs for", l),
           subtitle = "original") + 
      facet_grid(.~name) +
      theme_bw() +
      theme(axis.title = element_text(size = 8),
            title = element_text(size = 10))
    
    
    conta.gg.clean.sub <-  Rel.conta.prop.clean.sub %>% pivot_longer(cols = c(conta.max.prop)) %>% 
      mutate(control_type = ifelse(is.na(control_type), "sample", control_type)) %>% 
      ggplot(aes(x=control_type, y=value, color=control_type)) + 
      geom_boxplot() + geom_jitter(alpha=0.5) +
      geom_hline(yintercept = threshold.prop.cont, lty = "dashed", color="orange") +
      scale_color_manual(breaks = c("sample", "positive", "extraction", "pcr", "sequencing"), values = c("darkgray", "cyan4", "brown", "red", "pink"), na.value = "darkgrey") +
      labs(x=NULL, y="Prop. reads flagged as contaminant",
           title = paste("Prop. of flagged MOTUs for", l),
           subtitle = "after tagjumpslayer" ) + 
      facet_grid(.~name) +
      theme_bw() +
      theme(axis.title = element_text(size = 8),
            title = element_text(size = 10))
    
    conta.gg.sub <- ggpubr::ggarrange(conta.gg.ori.sub,
                                      conta.gg.clean.sub,
                                      nrow = 1, ncol = 2 , common.legend = T, legend = "bottom"
    )
    
    ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("04_conta.prop_",l,"_", x,".png")), 
           plot = conta.gg.sub ,
           width = 7,
           height = 3,
           bg = "white",
           units = c("in"))
    
    cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/04_conta.prop_",l,"_", x,".png"), "\n")  
    
  }
}


# Class MOTUs and Samples (pcr) in different categories --------------------------------------------------------

# Perform this by subgroups

for(l in LOCUS){
  
  cat("\nSummerizing contamination problems for", l, "\n\n")  
  
  for(x in SUBGROUP){
    
    cat("Running on", x, "subgroup\n")  
    
    metabarlist.int <- get(paste0("metabarlist.ori.",l, ".", x))
    metabarlist.int.clean <- get(paste0("metabarlist.tagclean.",l, ".", x))
    
    
    metabarlist.int$motus <- metabarlist.int$motus %>% mutate(artefact_type = ifelse(not_a_max_conta == FALSE, "Contamination",
                                                                                     ifelse(not_a_max_conta == TRUE, "Not artefactual", "Undefined??")))
    
    metabarlist.int.clean$motus <- metabarlist.int.clean$motus %>% mutate(artefact_type = ifelse(not_a_max_conta == FALSE, "Contamination",
                                                                                                 ifelse(not_a_max_conta == TRUE, "Not artefactual", "Undefined??")))
    
    
    summary.artefact.motus_N.ori <- metabarlist.int$motus %>%  
      group_by(artefact_type) %>% summarise(N = n()) %>% 
      mutate(SUM  = sum(N),
             prop = N / sum(N),
             dataset = "original",
             level = "motus")
    
    summary.artefact.motus_reads.ori <- metabarlist.int$motus %>%  
      group_by(artefact_type) %>% summarise(N = sum(count.subgroup)) %>% 
      mutate(SUM  = sum(N),
             prop = N / sum(N),
             dataset = "original",
             level = "reads")
    
    summary.artefact.motus_N.clean <- metabarlist.int.clean$motus %>%  
      group_by(artefact_type) %>% summarise(N = n()) %>% 
      mutate(SUM  = sum(N),
             prop = N / sum(N),
             dataset = "tagjump",
             level = "motus")
    
    summary.artefact.motus_reads.clean <- metabarlist.int.clean$motus %>%  
      group_by(artefact_type) %>% summarise(N = sum(count.tagjump.subgroup)) %>% 
      mutate(SUM  = sum(N),
             prop = N / sum(N),
             dataset = "tagjump",
             level = "reads")
    
    
    readr::write_csv(bind_rows(summary.artefact.motus_N.ori, summary.artefact.motus_reads.ori, summary.artefact.motus_N.clean, summary.artefact.motus_reads.clean), 
                     file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("05_Artefact_MOTUs_",l, "_", x,".csv")))
    
    
    
    graph.artefact.motus_N.ori <- summary.artefact.motus_N.ori %>% 
      ggplot(aes(x=1, y = prop, fill=artefact_type)) +
      geom_bar(stat = "identity") +#  xlim(0, 1) +
      labs(fill="Artifact type") + 
      coord_polar(theta="y") + theme_void() + 
      scale_fill_manual(limits = c("Not artefactual", "Contamination"), values = c("deepskyblue1", "darkorange1" )) + 
      geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 5) +
      ggtitle(paste("Prop. MOTUs for", l))
    
    graph.artefact.motus_reads.ori <- summary.artefact.motus_reads.ori %>% 
      ggplot(aes(x=1, y = prop, fill=artefact_type)) +
      geom_bar(stat = "identity") +#  xlim(0, 1) +
      labs(fill="Artifact type") + 
      coord_polar(theta="y") + theme_void() + 
      scale_fill_manual(limits = c("Not artefactual", "Contamination"), values = c("deepskyblue1", "darkorange1" )) + 
      geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 5) +
      ggtitle(paste("Prop. reads for", l))
    
    graph.artefact.motus_N.clean <- summary.artefact.motus_N.clean %>% 
      ggplot(aes(x=1, y = prop, fill=artefact_type)) +
      geom_bar(stat = "identity") +#  xlim(0, 1) +
      labs(fill="Artifact type") + 
      coord_polar(theta="y") + theme_void() + 
      scale_fill_manual(limits = c("Not artefactual", "Contamination"), values = c("deepskyblue1", "darkorange1" )) + 
      geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 5) +
      ggtitle(paste("Prop. MOTUs for", l))
    
    graph.artefact.motus_reads.clean <- summary.artefact.motus_reads.clean %>% 
      ggplot(aes(x=1, y = prop, fill=artefact_type)) +
      geom_bar(stat = "identity") +#  xlim(0, 1) +
      labs(fill="Category") + 
      coord_polar(theta="y") + theme_void() + 
      scale_fill_manual(limits = c("Not artefactual", "Contamination"), values = c("deepskyblue1", "darkorange1" )) + 
      geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 5) +
      ggtitle(paste("Prop. reads for", l))
    
    
    graph.artefact.motus.ori <- ggpubr::ggarrange(graph.artefact.motus_N.ori,
                                                  graph.artefact.motus_reads.ori,
                                                  nrow = 1, ncol = 2 , common.legend = T, legend = "bottom"
    )
    
    graph.artefact.motus.clean <- ggpubr::ggarrange(graph.artefact.motus_N.clean,
                                                    graph.artefact.motus_reads.clean,
                                                    nrow = 1, ncol = 2 , common.legend = T, legend = "bottom"
    )
    
    ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction",  paste0("05_Artefact_MOTUs_original_",l, "_", x,".png")), 
           plot = graph.artefact.motus.ori,
           width = 7,
           height = 3,
           bg = "white",
           units = c("in"))
    
    cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/05_Artefact_MOTUs_original_",l, "_", x,".png"), "\n") 
    
    ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction",  paste0("05_Artefact_MOTUs_tagjump_",l, "_", x,".png")), 
           plot = graph.artefact.motus.clean ,
           width = 7,
           height = 3,
           bg = "white",
           units = c("in"))
    
    cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/05_Artefact_MOTUs_tagjump_",l, "_", x,".png"), "\n")  
    
    metabarlist.int$pcrs <-  metabarlist.int$pcrs %>% mutate(artefact_type = ifelse(seqdepth_ok == FALSE & low_max_conta_level_10 == TRUE, "Low sequencing depth",
                                                                                    ifelse(low_max_conta_level_10 == FALSE& seqdepth_ok == TRUE, "Contamination > 10%",
                                                                                           ifelse(low_max_conta_level_10 == FALSE & seqdepth_ok == FALSE, "Contamination > 10% and low sequencing depth",        
                                                                                                  ifelse(low_max_conta_level_10 == TRUE & low_max_conta_level_10 == TRUE, "Not artefactual", "Undefined??")))))
    
    metabarlist.int.clean$pcrs <-  metabarlist.int.clean$pcrs %>% mutate(artefact_type = ifelse(seqdepth_ok == FALSE & low_max_conta_level_10 == TRUE, "Low sequencing depth",
                                                                                                ifelse(low_max_conta_level_10 == FALSE& seqdepth_ok == TRUE, "Contamination > 10%",
                                                                                                       ifelse(low_max_conta_level_10 == FALSE & seqdepth_ok == FALSE, "Contamination > 10% and low sequencing depth",        
                                                                                                              ifelse(low_max_conta_level_10 == TRUE & low_max_conta_level_10 == TRUE, "Not artefactual", "Undefined??")))))
    
    
    metabarlist.int.clean$pcrs %>% pull(artefact_type) %>% table()
    
    summary.artefact.pcr.ori <- metabarlist.int.clean$pcrs %>% dplyr::filter(type == "sample") %>% 
      group_by(artefact_type) %>% summarise(N = n()) %>% 
      mutate(SUM = sum(N),
             prop = N / sum(N)) %>% 
      mutate(dataset = "original")
    
    graph.artefact.pcr.ori <-  summary.artefact.pcr.ori  %>% 
      ggplot(aes(x=1, y = prop, fill=artefact_type)) +
      geom_bar(stat = "identity") +  #xlim(0, 1) +
      labs(fill="Category") + 
      coord_polar(theta="y") + theme_void() + 
      scale_fill_manual(limits = c("Not artefactual", "Low sequencing depth", "Contamination > 10%", "Contamination > 10% and low sequencing depth"), values = c("deepskyblue1", "darkgoldenrod1", "darkorange1", "darkred" )) + 
      geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 5) +
      ggtitle(paste("Prop. sample (pcr) for", l))
    
    summary.artefact.pcr.clean <- metabarlist.int.clean$pcrs %>% dplyr::filter(type == "sample") %>% 
      group_by(artefact_type) %>% summarise(N = n()) %>% 
      mutate(SUM = sum(N),
             prop = N / sum(N)) %>% 
      mutate(dataset = "tagjump")
    
    graph.artefact.pcr.clean <-  summary.artefact.pcr.clean  %>% 
      ggplot(aes(x=1, y = prop, fill=artefact_type)) +
      geom_bar(stat = "identity") +  #xlim(0, 1) +
      labs(fill="Category") + 
      coord_polar(theta="y") + theme_void() + 
      scale_fill_manual(limits = c("Not artefactual", "Low sequencing depth", "Contamination > 10%", "Contamination > 10% and low sequencing depth"), values = c("deepskyblue1", "darkgoldenrod1", "darkorange1", "darkred" )) + 
      geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 5) +
      ggtitle(paste("Prop. sample (pcr) for", l))
    graph.artefact.pcr.clean
    
    # Save summary for the automatic report
    readr::write_csv(bind_rows(summary.artefact.pcr.ori, summary.artefact.pcr.clean), 
                     file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("05_Artefact_PCRs_",l, "_", x,".csv")))
    
    ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("05_Artefact_PCRs_original_",l, "_", x,".png")), 
           plot = graph.artefact.pcr.ori,
           width = 5,
           height = 4,
           units = c("in"), bg = "white")
    
    cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/05_Artefact_PCRs_original_",l, "_", x,".png"), "\n")  
    
    
    ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("05_Artefact_PCRs_tagjump_",l, "_", x,".png")), 
           plot = graph.artefact.pcr.clean,
           width = 5,
           height = 4,
           units = c("in"), bg = "white")
    
    cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/05_Artefact_PCRs_tagjump_",l, "_", x,".png"), "\n")  
    
    # Export to keep the categories
    
    assign(x = paste0("metabarlist.ori.", l, ".", x), 
           value = metabarlist.int)
    
    assign(x = paste0("metabarlist.tagclean.", l, ".", x), 
           value = metabarlist.int.clean)
    
    # Export!!
    
    metabarlist.int$reads %>% write.csv( file = file.path(here::here(), "00_Data/04_ESVcorrected", paste0("ESVtab.raw_",l,"_",x, ".csv")))
    metabarlist.int$pcrs %>% readr::write_csv(file.path(here::here(),"00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.raw_",l, "_",x, ".csv")))
    metabarlist.int$motus %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("MOTUs.Metabarinfo.raw_",l, "_",x, ".csv")))
    metabarlist.int$pcrs %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.raw_",l, "_",x, ".csv")))
    
    metabarlist.int.clean$reads %>% write.csv( file = file.path(here::here(), "00_Data/04_ESVcorrected", paste0("ESVtab.postTagjump_",l,"_",x, ".csv")))
    metabarlist.int.clean$pcrs %>% readr::write_csv(file.path(here::here(),"00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.postTagjump_",l, "_",x, ".csv")))
    metabarlist.int.clean$motus %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("MOTUs.Metabarinfo.postTagjump_",l, "_",x, ".csv")))
    metabarlist.int.clean$pcrs %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.postTagjump_",l, "_",x, ".csv")))
    
    cat("Data (.csv) are exported here:", paste0("00_Data/04_ESVcorrected/") , "\n")  
    
    
    
    
  }
}


# Apply the corrections ----------------------------------------------------

#

#metabarlist.int.clean <- subset_metabarlist(metabarlist.int, "reads", 

for(l in LOCUS){
  
  cat("\nFinal step for", l, "\n\n")  
  
  # Load parameters
  
  thresholds.tag <-  metabar.param %>% filter(Locus == l) %>% pull(tag.threshold)
  motus.correct <- metabar.param %>% filter(Locus == l) %>% pull(motus.correct)
  pcr.correct <- metabar.param %>% filter(Locus == l) %>% pull(pcr.correct)
  
  tag.correct <-  metabar.param %>% filter(Locus == l) %>% pull(tag.correct)
  
  
  for(x in SUBGROUP){
    
    cat("Running on", x, "subgroup\n")  
    
    # id.int <-  SUBGROUP.ls[[x]]
    
    metabarlist.int.sub <- get(paste0("metabarlist.ori.",l, ".", x))
    
    #metabarlist.int.sub$pcrs$type
    
    metabarlist.int.sub <- subset_metabarlist(  metabarlist.int.sub, table="pcrs", 
                                                indices = (metabarlist.int.sub$pcrs$type == "sample" ) ) 
    
    
    
    # Choose the right dataset
    if(tag.correct == T){
      
      cat("Working on the tag jump clean dataset\n")  
      metabarlist.correct.int <- get(paste0("metabarlist.tagclean.",l, ".", x))
    } else {
      cat("Working on the original dataset (no tag jump correction)\n")       
      metabarlist.correct.int <- get(paste0("metabarlist.ori.",l, ".", x))
    } 
    
    if(motus.correct == T){ 
      
      # Subset on MOTUs and SAMPLE 
      cat("Keeping only not artefactual MOTUs\n")  
      metabarlist.correct.int <- subset_metabarlist(metabarlist.correct.int, table="motus", 
                                                    indices = (metabarlist.correct.int$motus$artefact_type == "Not artefactual") )
      
    } else{
      cat("No filtration on MOTUs \n") 
    }  
    
    if(pcr.correct == T){   
      
      metabarlist.correct.int <- subset_metabarlist(metabarlist.correct.int, table="pcrs", 
                                                    indices = (metabarlist.correct.int$pcrs$artefact_type == "Not artefactual"))
      
    }  else {
      cat("No filtration on samples (pcrs) \n") 
    }  
    
    
    if(nrow(metabarlist.correct.int$pcrs) >0) { # Stop here if no sample were keep
      
      # Keep only sample
      metabarlist.correct.int <- subset_metabarlist(metabarlist.correct.int, table="pcrs", 
                                                    indices = (metabarlist.correct.int$pcrs$type == "sample" ) )  
      
      
      
      if(nrow(metabarlist.correct.int$pcrs) > 1) { # Stop here if no sample were keep
        
        # Add stats on filtration
        
        metabarlist.correct.int$motus$count_postmetbaR = colSums(metabarlist.correct.int$reads)
        metabarlist.correct.int$pcrs$nb_reads_postmetabaR = rowSums(metabarlist.correct.int$reads)
        metabarlist.correct.int$pcrs$nb_motus_postmetabaR = rowSums(ifelse(metabarlist.correct.int$reads>0, T, F))
        
        
        summary.int <- metabarlist.correct.int$pcrs %>% select(sample_id, nb_motus, nb_reads, nb_reads_postmetabaR, nb_motus_postmetabaR) %>% 
          mutate(Loci = l,
                 ID_subproject = x)
        
        
        readr::write_csv(summary.int, file.path(here::here(), "00_Data", "04_ESVcorrected", paste0("ESVtab_Stats_postMetabaR_", l, "_",x, ".csv")))
        cat("Summary stats were saved here:", paste0("00_Data/04_ESVcorrected/ESVtab_Stats_postMetabaR_", l, "_",x, ".csv") , "\n")  
        
        # 
        check.correction <- reshape2::melt(metabarlist.correct.int$pcrs[,c("nb_reads", "nb_reads_postmetabaR", 
                                                                           "nb_motus", "nb_motus_postmetabaR")])
        check.correction$type <- ifelse(grepl("motus", check.correction$variable), "richness", "abundance")
        
        post.correction.gg <- ggplot(data = check.correction, aes(x = variable, y = value)) +
          geom_boxplot( color = "darkgrey") +
          geom_jitter(alpha=0.1, color = "darkgrey") +
          theme_bw() +
          facet_wrap(~type, scales = "free", ncol = 5) +
          theme(axis.text.x = element_text(angle=45, h=1)) +
          ggtitle(paste("Comparison before/after correction for", l))
        
        
        ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("06_summary.postcorrection_",l, "_",x, ".png")), 
               plot = post.correction.gg,
               width = 5,
               height = 3,
               units = c("in"), bg = "white")
        
        
        cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/06_summary.postcorrection_",l, "_", x,".png"), "\n")  
        
        # Export corrected the results 
        metabarlist.correct.int$reads %>% write.csv( file = file.path(here::here(), "00_Data/04_ESVcorrected", paste0("ESVtab.corrected_",l,"_",x, ".csv")))
        metabarlist.correct.int$pcrs %>% readr::write_csv(file.path(here::here(),"00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.corrected_",l, "_",x, ".csv")))
        metabarlist.correct.int$motus %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("MOTUs.Metabarinfo.corrected_",l, "_",x, ".csv")))
        metabarlist.correct.int$pcrs %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.corrected_",l, "_",x, ".csv")))
        cat("Data (.csv) are exported here:", paste0("00_Data/04_ESVcorrected/") , "\n")  
        
        # Final visualisation
        
        read.correct.tidy <- metabarlist.correct.int$reads %>% as.data.frame() %>% 
          mutate(ID_labo = row.names(metabarlist.correct.int$reads)) %>% 
          pivot_longer(-ID_labo, names_to = "QueryAccVer", values_to = "Nreads") %>% 
          left_join(metabarlist.correct.int$motus %>% select(QueryAccVer, Taxon, genus, phylum)) %>% 
          mutate(Taxon = ifelse(is.na(Taxon), "Unassigned", Taxon)) %>% 
          group_by(ID_labo, Taxon, phylum) %>% summarise(Nreads = sum(Nreads)) %>% 
          left_join(data.info %>% dplyr::filter(Loci == l)) %>% filter(Type_echantillon %in% c("Echantillon", "ECH"))
        
        read.ori.tidy <- metabarlist.int.sub$reads %>% as.data.frame() %>% 
          mutate(ID_labo = row.names(metabarlist.int.sub$reads)) %>% 
          pivot_longer(-ID_labo, names_to = "QueryAccVer", values_to = "Nreads") %>% 
          left_join(metabarlist.int.sub$motus %>% select(QueryAccVer, Taxon, genus, phylum)) %>% 
          mutate(Taxon = ifelse(is.na(Taxon), "Unassigned", Taxon)) %>% 
          group_by(ID_labo, Taxon, phylum) %>% summarise(Nreads = sum(Nreads)) %>% 
          left_join(data.info %>% dplyr::filter(Loci == l)) %>% dplyr::filter(Type_echantillon %in% c("Echantillon", "ECH"),
                                                                              Nreads > 0)
        
        if(nrow( read.correct.tidy) > 0){
          final.correction.gg  <- read.correct.tidy %>%  ggplot(aes(fill = Nreads, x = ID_labo, y = Taxon)) +
            labs(x= "", y = "") + 
            geom_bin2d(color = "darkgray")+
            scale_fill_distiller(trans = "log10",
                                 palette = "Spectral",
                                 na.value = "gray95"#,
                                 #breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000), labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000", "10,000,000")
            ) +
            theme_minimal()+
            facet_grid(phylum ~ ID_projet, scale = "free", space = "free") + 
            ggtitle(paste("Overall visualisation after correction for", l)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  legend.position = "right")
        } else{final.correct.gg <-  ggplot()}
        
        if(nrow( read.ori.tidy) > 0 ){
          final.ori.gg <- read.ori.tidy %>%  ggplot(aes(fill = Nreads, x = ID_labo, y = Taxon)) +
            labs(x= "", y = "") + 
            geom_bin2d(color = "darkgray")+
            scale_fill_distiller(trans = "log10",
                                 palette = "Spectral",
                                 na.value = "gray95"#,
                                 #breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000), labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000", "10,000,000")
            ) +
            theme_minimal() +
            facet_grid(phylum ~ ID_projet, scale = "free", space = "free") +   
            
            ggtitle(paste("Overall visualisation before correction for", l)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  legend.position = "right")
        } else{final.ori.gg <- ggplot()}
        
        
        ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("06_overall.postcorrection_",l,"_",x,  ".png")), 
               plot = final.correction.gg,
               width = 8,
               height = 8,
               units = c("in"), bg = "white")
        
        cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/06_overall.postcorrection_",l, "_", x,".png"), "\n")  
        
        ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("06_overall.original_",l, "_",x, ".png")), 
               plot = final.ori.gg,
               width = 8,
               height = 8,
               units = c("in"), bg = "white")
        
        cat("Figure saved:", paste0("02_Results/04_ESVtable_correction/06_overall.original_",l, "_", x,".png"), "\n")  
        
        assign(x = paste0("metabarlist.correct.", l, ".", x), 
               value = metabarlist.correct.int )
        
      }
      
    }
    
  } # End of project loop
} # End of overall loop  

# Save overall file --------------------------------------------------------------


for(l in LOCUS){
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  metabarlist.int.clean <- get(paste0("metabarlist.tagclean.",l))
  
  
  metabarlist.int$reads %>% write.csv( file = file.path(here::here(), "00_Data/04_ESVcorrected", paste0("ESVtab.raw_",l,"_ALL.csv")))
  metabarlist.int$pcrs %>% readr::write_csv(file.path(here::here(),"00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.raw_",l, "_ALL.csv")))
  metabarlist.int$motus %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("MOTUs.Metabarinfo.raw_",l, "_ALL.csv")))
  metabarlist.int$pcrs %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.raw_",l, "_ALL.csv")))
  
  metabarlist.int.clean$reads %>% write.csv( file = file.path(here::here(), "00_Data/04_ESVcorrected", paste0("ESVtab.postTagjump_",l,"_ALL.csv")))
  metabarlist.int.clean$pcrs %>% readr::write_csv(file.path(here::here(),"00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.postTagjump_",l, "_ALL.csv")))
  metabarlist.int.clean$motus %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("MOTUs.Metabarinfo.postTagjump_",l, "_ALL.csv")))
  metabarlist.int.clean$pcrs %>% readr::write_csv(file.path(here::here(), "00_Data/04_ESVcorrected", paste0("Samples.Metabarinfo.postTagjump_",l, "_ALL.csv")))
  
  cat("Data (.csv) are exported here:", paste0("00_Data/04_ESVcorrected/") , "\n")  
  
}



# Write a final log

cat("\nEND of 04_ESVtable_correction.R script\n",
    date(),
    "\n-------------------------\n", 
    
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("metabaR", packageVersion("metabaR"), sep = ": "),     
    # Add it to the log file
    file = file.path(here::here(), "00_Data", "04_ESVcorrected", "ESVtab_correction.log"), 
    append = F, sep = "\n")
