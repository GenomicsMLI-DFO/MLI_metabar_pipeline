# MLI_metabar_pipeline
Template repository to perform metabarcoding analysis. 


This pipeline is intended to run in R (and Rstudio), but need external programs such as FastQC (REF) and cutadapt (REF). 


## How to use MLI_metabar_pipeline

### Pre-requisite

- R and Rstudio
- FastQC
- MultiQC
- cutadapt

Check that this work :
system2("fastqc", "--help")

MultiQC and cutadapt from python environment - should be added to the R path. Can be done as:

# Add python env to this specific project
Sys.setenv(PATH = paste(c("/home/genleia/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))



### Before starting an analysis

- Put raw read files in 00_Data/01a_RawData (see examples)
- Put metadata in 00_Data/00_FileInfos (see examples)

SeqInfo.csv

- Install the depending R package : `...`

This can be done all at once with this command line in R :

```{r}
install.packages(c("readr", "tidyr", "magrittr", "dplyr", "stringr", "here", "parallel", "ggplot2"))

```
Somes are from biostrings

```{r}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
BiocManager::install("Biostrings")
BiocManager::install("dada2")

```


### Rename raw files

Use **01_Rename_RAW.R** within *01_Code* folder to rename zipped fastq files. It will remove these patterns :

MiSeq : "MI.M*00000*_*0000*.*000*.FLD*0000*."
NovaSeq : "NS.*0000*.*000*.FLD*0000*.*0000*---PE1-CS1-IDT_i5_*0*."

Other patterns or transformation can be implemented

### From **raw reads** to **ESV table**

Use the file **02_Process_RAW.R** within *01_Code* folder to transform raw reads into an ESV table. 

These are the steps:
1. fastQC and multiQC on raw reads
2. Cutdapt to check for and remove adaptors (followed by a second fastQC/multiQC)
3. dada2 filtering (followed by a third fastQC/multiQC)
4. dada2 error rate assessment
5. dada2 dereplication, sample inference and merging
6. dada2 chimera removal

### Deal with negative samples

## Example



## References
