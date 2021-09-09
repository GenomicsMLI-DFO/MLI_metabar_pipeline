# MLI_metabar_pipeline
Template repository to perform metabarcoding analysis. 


This pipeline is intended to run in R (and Rstudio), but need external programs such as FastQC (REF) and cutadapt (REF). 


## How to use MLI_metabar_pipeline

### Pre-requisite

- R and Rstudio
- FastQC
- cutadapt

### Before starting an analysis

- Put raw read files in 00_Data/01a_RawData (see examples)
- Put metadata in 00_Data/00_FileInfos (see examples)

- Install the depending R package : `...`

This can be done all at once with this command line in R :

```{r}
install.packages(c("..."))
```

### Rename raw files

Use **Rename_RAW.R** within *01_Code* folder to rename zipped fastq files. It will remove these patterns :

MiSeq : "MI.M*00000*_*0000*.*000*.FLD*0000*."
NovaSeq : "NS.*0000*.*000*.FLD*0000*.*0000*---PE1-CS1-IDT_i5_*0*."

Other patterns or transformation can be implemented

### From **raw reads** to **ESV table**

Use the file **Process_RAW.R** within *01_Code* folder to transform raw reads into an ESV table. 

### Deal with negative samples

## Example



## References
