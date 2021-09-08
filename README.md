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
- Use **Rename_RAW.R** to rename zipped fastq files
- Install the depending R package : `...`

This can be done all at once with this command line in R :

```{r}
install.packages(c("..."))
```


### From **raw reads** to **ESV table**

Use the file **Process_RAW.R** to transform raw reads into an ESV table. 

### Deal with negative samples

## Example



## References
