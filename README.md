# MLI metabarcoding pipeline
Template repository to perform metabarcoding analysis. 

This pipeline is intended to run in R (and Rstudio), but need external programs (see pre-requisite section). 

## How to use MLI_metabar_pipeline

### Pre-requisite

- R and Rstudio
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

To be sure that the external command are found by R, try to run these commands first:

```{r}
system2("fastqc", "--help")
system2("multiqc", "--help")
system2("cutadapt", "--help")
```

MultiQC and cutadapt can be installed in a python environment that should be added to the R path. Can be done as:

```{r}
Sys.setenv(PATH = paste(c("/path/to/PythonVenv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))
```
### Before starting an analysis

- Put raw read files in 00_Data/01a_RawData (see examples)
- Put metadata in 00_Data/00_FileInfos (see examples)

SeqInfo.csv

- Install the depending R package : `readr`, `tidyr`, `magrittr`,`dplyr`,`stringr`,`here`,`parallel`, `ggplot2`

This can be done all at once with this command line in R :

```{r}
install.packages(c("readr", "tidyr", "magrittr", "dplyr", "stringr", "here", "parallel", "ggplot2"))
```
Somes are from biostrings : `dada2`, `Biostrings`

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
BiocManager::install("Biostrings")
BiocManager::install("dada2")
```


### Rename raw files

To be process, fastq files must be rename as **SAMPLENAME_MARKER_R1or2.fastq.gz**. This can be done with one of the two following scripts:

#### Rename option 1 : Only remove sequencer/run id pattern

Use **01_Rename_RAW_SeqPattern.R** within *01_Code* folder to rename zipped fastq files. It will remove these patterns :

MiSeq

"MI.M*00000*_*0000*.*000*.FLD*0000*."

NovaSeq 

"NS.*0000*.*000*.FLD*0000*.*0000*---PE1-CS1-IDT_i5_*0*."

"NS.*0000*.FLD*0000*---PE1-CS1-IDT_i5_*0*."

Other patterns or transformation can be implemented

#### Rename option 2 : Start from file ID name 

Use **01_Rename_RAW_FileName.R** within *01_Code* folder to rename zipped fastq files. It will used the metadata file (*SeqInfo.csv*) to create a new name. File name must exclude the **_R1orR2.fastq.gz** part.

### From **raw reads** to **ESV table**

Use the file **02_Process_RAW.R** within *01_Code* folder to transform raw reads into an ESV table. 

These are the steps:
1. fastQC and multiQC on raw reads
2. Cutdapt to check for and remove adaptors (followed by a second fastQC/multiQC)
3. dada2 filtering (followed by a third fastQC/multiQC)
4. dada2 error rate assessment
5. dada2 dereplication, sample inference and merging
6. dada2 chimera removal

Then you can use **03_Raw_report.Rmd** within *01_Code* folder to get a first view on what have been done.

### Deal with negative samples

## Example

A test dataset is available with the template pipeline.

These files are withe the **00_Data/01a_Raw_Data** folder.

4 samples at the COI marker, both F (R1) and R (R2).

## References
