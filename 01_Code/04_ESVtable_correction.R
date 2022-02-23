
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

data.info <- readr::read_csv(file.path(here::here(), "00_Data", "00_FileInfos", "SeqInfo.csv") )
data.info

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
  
  motus<-   motus.int %>% left_join(RES.all.ncbi %>% filter(Loci == l, Method == "TOP", Threshold == 95))
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
                                         project = ID_projet) %>% 
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


# RUN METABAR -------------------------------------------------------------


# Print stats and add a few count columns
for(l in LOCUS){

cat("\nStatistic for", l, "\n\n")  
  
metabarlist.int <- get(paste0("metabarlist.ori.",l))

print(summary_metabarlist(metabarlist.int))

metabarlist.int$pcrs$nb_reads <- rowSums(metabarlist.int$reads)
metabarlist.int$pcrs$nb_motus <- rowSums(metabarlist.int$reads>0)
metabarlist.int$motus$count   <- colSums(metabarlist.int$reads)

# Save image

plate.gg <- ggpcrplate.modif(metabarlist.int, legend_title = "N reads")

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("00_plate_ori_",l, ".png")), 
       plot = plate.gg ,
       width = 6,
       height = 6,
       units = c("in"))



assign(x = paste0("metabarlist.ori.", l), 
       value = metabarlist.int )

}

# Flag low read depth outlier ---------------------------------------------

depth.threshold <- 1000

depth.gg <- list()

for(l in LOCUS){
  
  cat("\nSequencing depth for", l, "\n\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  

gg.int <- ggplot(metabarlist.int$pcrs, aes(x = nb_reads+1, fill = control_type)) +
  geom_histogram(bins=40, color="grey") + 
  geom_vline(xintercept = depth.threshold, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  scale_fill_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x="N reads + 1 (log)", 
       y="N samples") +
  labs(title = paste("Sequencing depth for", l)) +
  theme_bw() + 
  theme(panel.grid = element_blank())

depth.gg[[2]] <- gg.int

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
metabarlist.int$pcrs$seqdepth_ok <- ifelse(metabarlist.int$pcrs$nb_reads < depth.threshold, F, T)


assign(x = paste0("metabarlist.ori.", l), 
       value = metabarlist.int )

}


graph.depth <- ggpubr::ggarrange(plotlist = depth.gg)

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", "01_depth_allsamples.png"), 
       plot = graph.depth,
       width = 6,
       height = 6,
       units = c("in"))


# Flag contaminants -------------------------------------------------------

#conta.gg <- list()

for(l in LOCUS){

  cat("\nLooking at contaminants for", l, "\n\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))

  # Run contaslayer
    
metabarlist.int  <- contaslayer(metabarlist.int , method="all",
                        control_types = c("pcr", "sequencing", "extraction"),
                        output_col = "not_a_all_conta")

metabarlist.int  <- contaslayer(metabarlist.int , method="max",
                                control_types = c("pcr", "sequencing", "extraction"),
                                output_col = "not_a_max_conta")

common_contam.all <- metabarlist.int$motus %>% filter(not_a_all_conta == F) %>% select(Taxon, Levels, count) %>% mutate(method = "all") %>% arrange(desc(count)) 
common_contam.max <- metabarlist.int$motus %>% filter(not_a_max_conta == F) %>% select(Taxon, Levels, count) %>% mutate(method = "max") %>% arrange(desc(count))

common_contam.all$ESV <- row.names(common_contam.all)
common_contam.max$ESV <- row.names(common_contam.max)

readr::write_csv(bind_rows(common_contam.all, common_contam.max) %>% arrange(desc(count)), 
    file = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("02_conta_",l, ".csv")))

conta.gg <- ggpcrplate.cont(metabarlist.int, N = 20)

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("02_plate_conta_",l, ".png")), 
       plot = conta.gg ,
       width = 6,
       height = 6,
       units = c("in"))

assign(x = paste0("metabarlist.ori.", l), 
       value = metabarlist.int )

}

# Flag high % contaminant samples -------------------------------------------------

threshold.prop.cont <- 0.1

conta.gg <- list()

for(l in LOCUS){
  
  cat("\nLooking at samples with high contaminants for", l, "\n\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  
Rel.conta.prop <- data.frame(conta.all.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_all_conta=="FALSE"]) / rowSums(metabarlist.int$reads),
                             conta.max.prop = rowSums(metabarlist.int$reads[,metabarlist.int$motus$not_a_max_conta=="FALSE"]) / rowSums(metabarlist.int$reads))

# Add information on control types
Rel.conta.prop$control_type <- metabarlist.int$pcrs$control_type[match(rownames(Rel.conta.prop), rownames(metabarlist.int$pcrs))]

conta.gg[[l]] <-  Rel.conta.prop %>% pivot_longer(cols = c(conta.max.prop, conta.all.prop)) %>% 
  ggplot(aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + geom_jitter(alpha=0.5) +
  geom_hline(yintercept = threshold.prop.cont, lty = "dashed", color="orange") +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. reads flagged as contaminant by sample",
       title = paste("Abundance of flagged MOTUs for", l)) + 
  facet_grid(.~name) +
  theme_bw() 

metabarlist.int$pcrs$low_all_conta_level_10 <- ifelse(Rel.conta.prop$conta.all.prop[match(rownames(metabarlist.int$pcrs), rownames(Rel.conta.prop))]>0.1,  F, T)
metabarlist.int$pcrs$low_max_conta_level_10 <- ifelse(Rel.conta.prop$conta.max.prop[match(rownames(metabarlist.int$pcrs), rownames(Rel.conta.prop))]>0.1,  F, T)

assign(x = paste0("metabarlist.ori.", l), 
       value = metabarlist.int )

}

graph.prop.conta <- ggpubr::ggarrange(plotlist = conta.gg)

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", "03_conta.prop_allsamples.png"), 
       plot = graph.prop.conta ,
       width = 6,
       height = 6,
       units = c("in"))



# Tag jump ----------------------------------------------------------------

# Define a vector of thresholds to test
thresholds.tag.test <- c(0, 0.0001, 0.001, 0.01, 0.03, 0.05) 
thresholds.tag <- 0.01

for(l in LOCUS){
  
  cat("\nLooking at different tag jump threshold for", l, "\n\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))


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

tag.gg <- ggplot(tests.tagjump.long.2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept =  factor(thresholds.tag), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), height = 0, alpha=0.5) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink","green", "black","yellow","purple"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=2) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="MOTUs",
       title = paste("Tested threshold for", l)) + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("04_tagjump.threshold_",l, ".png")), 
       plot = tag.gg ,
       width = 6,
       height = 6,
       units = c("in"))

# Check

metabarlist.int.clean2 <- tagjumpslayer(metabarlist.int, 0.01)

# identify occurrence of the most abundant OTU
idx <- which.max(metabarlist.int$motus$count)
p1 <- ggpcrplate.modif(metabarlist.int,
                 legend_title = "# reads",
                 FUN = function(m) {
                   m$reads[, idx]
                 }
)



# same on clean data
p2 <- ggpcrplate.modif(metabarlist.int.clean2,
                 legend_title = "# reads",
                 FUN = function(m) {
                   m$reads[, idx]
                 }
)

p2 + scale_size(limits = c(1, max(metabarlist.int$reads[, idx]))) +
  ggtitle("Distribution of the most abundant MOTU after curation")

plate.tag.gg <- ggpubr::ggarrange(p1 + scale_size(limits = c(1, max(metabarlist.int$reads[, idx]))) +
                                    ggtitle("Most abundant MOTU ori"),
                                  p2 + scale_size(limits = c(1, max(metabarlist.int$reads[, idx]))) +
                                    ggtitle("Most abundant MOTU after tagjumpslayer"),
                                  nrow = 1, ncol = 2 , common.legend = T
                                  )


ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("04_tagjump.plate_",l, ".png")), 
       plot = plate.tag.gg ,
       width = 8,
       height = 8,
       units = c("in"), bg = "white")


}

# Class MOTUs and Samples (pcr) in different categories --------------------------------------------------------

for(l in LOCUS){
  
  cat("\nSummerizing contamination problems for", l, "\n\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))

metabarlist.int$motus <- metabarlist.int$motus %>% mutate(artefact_type = ifelse(not_a_all_conta == FALSE & not_a_max_conta == FALSE, "Control contamination (max and all)",
                                                                          ifelse(not_a_all_conta == FALSE & not_a_max_conta == TRUE, "Control contamination (all)",
                                                                          ifelse(not_a_all_conta == TRUE & not_a_max_conta == FALSE, "Control contamination (max)",
                                                                          ifelse(not_a_all_conta == TRUE & not_a_max_conta == TRUE, "Not artefactual", "Undefined??")))))

graph.artefact.motus <- metabarlist.int$motus %>%  ggplot(aes(x=1, fill=artefact_type)) +
  geom_bar() + # xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  #theme(legend.direction = "vertical") + 
  ggtitle(paste("MOTUs artefacts overview for", l))


graph.artefact.motus

metabarlist.int$pcrs <-  metabarlist.int$pcrs %>% mutate(artefact_type = ifelse(seqdepth_ok == FALSE, "Too low sequencing depth",
                                                           ifelse(low_all_conta_level_10 == FALSE & low_max_conta_level_10 == FALSE, "Contamination > 10% (max and all)",
                                                                                ifelse(low_all_conta_level_10 == FALSE & low_max_conta_level_10 == TRUE, "Contamination > 10%  (all)",
                                                                                       ifelse(low_all_conta_level_10 == TRUE & low_max_conta_level_10 == FALSE, "Contamination > 10%  (max)",
                                                                                              ifelse(low_all_conta_level_10 == TRUE & low_max_conta_level_10 == TRUE, "Not artefactual", "Undefined??"))))))


graph.artefact.pcr <-  metabarlist.int$pcrs %>% dplyr::filter(type == "sample") %>% 
  ggplot(aes(x=1, fill=artefact_type)) +
  geom_bar() +#  xlim(0, 1) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set2") + 
  #theme(legend.direction = "vertical") + 
  ggtitle(paste("Sample (pcr) artefacts overview for", l))

graph.artefact.pcr

graph.artefact.pcr + facet_wrap(~ project)

graph.artefact <- ggpubr::ggarrange(graph.artefact.motus + theme(legend.position = "bottom"), 
                  graph.artefact.pcr  + theme(legend.position = "bottom"), align = c("hv"))


graph.artefact.pcr.project <- graph.artefact.pcr + facet_wrap(~ project)

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("05_summary.conta_",l, ".png")), 
       plot = graph.artefact  ,
       width = 8,
       height = 5,
       units = c("in"), bg = "white")

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("05_summary.conta_byProject_",l, ".png")), 
       plot = graph.artefact.pcr.project  ,
       width = 8,
       height = 8,
       units = c("in"), bg = "white")

assign(x = paste0("metabarlist.ori.", l), 
       value = metabarlist.int )

}



# Apply the corrections ----------------------------------------------------

# Step 1 : Take the right tagjumpslayer dataset


#metabarlist.int.clean <- subset_metabarlist(metabarlist.int, "reads", 
#                                            indices = (rowSums(metabarlist.int$reads)>0))
thresholds.tag

for(l in LOCUS){
  
  cat("\nFinal step for", l, "\n\n")  
  
  metabarlist.int <- get(paste0("metabarlist.ori.",l))
  


# Run the tests and stores the results in a list
metabarlist.correct.int <- tagjumpslayer(metabarlist.int, thresholds.tag)

summary_metabarlist(metabarlist.int)

summary_metabarlist(metabarlist.correct.int)


# Subset on MOTUs and SAMPLE 

metabarlist.correct.int <- subset_metabarlist(metabarlist.correct.int, table="motus", 
                                              indices = (metabarlist.correct.int$motus$artefact_type == "Not artefactual") )

metabarlist.correct.int <- subset_metabarlist(metabarlist.correct.int, table="pcrs", 
                                              indices = (metabarlist.correct.int$pcrs$artefact_type == "Not artefactual" &
                                                           metabarlist.correct.int$pcrs$type == "sample" ) )

# Add stats on filtration

metabarlist.correct.int$motus$count = colSums(metabarlist.correct.int$reads)
metabarlist.correct.int$pcrs$nb_reads_postmetabaR = rowSums(metabarlist.correct.int$reads)
metabarlist.correct.int$pcrs$nb_motus_postmetabaR = rowSums(ifelse(metabarlist.correct.int$reads>0, T, F))



summary_metabarlist(metabarlist.int)
summary_metabarlist(metabarlist.correct.int)


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


ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("06_summary.postcorrection_",l, ".png")), 
       plot = post.correction.gg,
       width = 6,
       height = 6,
       units = c("in"), bg = "white")

# Export the results 
metabarlist.correct.int$reads %>% write.csv( file = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("ESVtab.corrected_",l, ".csv")))

metabarlist.int$motus %>% readr::write_csv(file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("MOTUs.Metabarinfo_",l, ".csv")))

metabarlist.int$pcrs %>% readr::write_csv(file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("Samples.Metabarinfo.all_",l, ".csv")))
metabarlist.correct.int$pcrs %>% readr::write_csv(file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("Samples.Metabarinfo.corrected_",l, ".csv")))

# Final visualisation

read.correct.tidy <- metabarlist.correct.int$reads %>% as.data.frame() %>% 
                 mutate(ID_labo = row.names(metabarlist.correct.int$reads)) %>% 
                 pivot_longer(-ID_labo, names_to = "QueryAccVer", values_to = "Nreads") %>% 
                 left_join(metabarlist.correct.int$motus %>% select(QueryAccVer, Taxon, phylum)) %>% 
                 mutate(Taxon = ifelse(is.na(Taxon), "Unassigned", Taxon)) %>% 
                 group_by(ID_labo, Taxon, phylum) %>% summarise(Nreads = sum(Nreads)) %>% 
                 left_join(data.info)

read.ori.tidy <- metabarlist.int$reads %>% as.data.frame() %>% 
                 mutate(ID_labo = row.names(metabarlist.int$reads)) %>% 
                 pivot_longer(-ID_labo, names_to = "QueryAccVer", values_to = "Nreads") %>% 
                 left_join(metabarlist.int$motus %>% select(QueryAccVer, Taxon, phylum)) %>% 
                 mutate(Taxon = ifelse(is.na(Taxon), "Unassigned", Taxon)) %>% 
                 group_by(ID_labo, Taxon, phylum) %>% summarise(Nreads = sum(Nreads)) %>% 
                 left_join(data.info) %>% dplyr::filter(Type_echantillon %in% c("Echantillon", "ECH"))


final.correction.gg  <- read.correct.tidy %>%  ggplot(aes(fill = Nreads, x = ID_labo, y = Taxon)) +
  labs(x= "", y = "") + 
  geom_bin2d(color = "darkgray")+
  scale_fill_distiller(trans = "log10",
                       palette = "Spectral",
                       na.value = "gray95",
                       breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000), labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000", "10,000,000")) +
  theme_minimal()+
  facet_grid(phylum ~ ID_projet, scale = "free", space = "free") + 
  ggtitle(paste("Overall visualisation after correction for", l)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "right")

final.ori.gg <- read.ori.tidy %>%  ggplot(aes(fill = Nreads, x = ID_labo, y = Taxon)) +
  labs(x= "", y = "") + 
  geom_bin2d(color = "darkgray")+
  scale_fill_distiller(trans = "log10",
                       palette = "Spectral",
                       na.value = "gray95",
                       breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000), labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000", "10,000,000")) +
  theme_minimal()+
  facet_grid(phylum ~ ID_projet, scale = "free", space = "free") + 
  ggtitle(paste("Overall visualisation before correction for", l)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "right")

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("06_overall.postcorrection_",l, ".png")), 
       plot = final.correction.gg,
       width = 8,
       height = 8,
       units = c("in"), bg = "white")

ggsave(filename = file.path(here::here(), "02_Results/04_ESVtable_correction", paste0("06_overall.original_",l, ".png")), 
       plot = final.ori.gg,
       width = 8,
       height = 8,
       units = c("in"), bg = "white")



assign(x = paste0("metabarlist.correct.", l), 
       value = metabarlist.correct.int )

}



# Save Rdata --------------------------------------------------------------

rm(list = c("metabarlist.int", "metabarlist.correct.int", "metabarlist.int.clean2"))

save(list = ls(pattern = "metabarlist"), 
     file = file.path(here::here(), "02_Results/04_ESVtable_correction","metabarlist.Rdata"))



