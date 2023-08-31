

# Info --------------------------------------------------------------------

# Script to check that all needed external programs can be reach by R
# If you reach an error message, it's because the program is not install 
# system wide
# 
# Audrey Bourret
# 2023-08-31


# Library -----------------------------------------------------------------

# Internal functions
source(file.path(here::here(), "01_Code", "Functions", "get.value.R"))


# Step 0: Python virtual environment --------------------------------------

# If some program (e.g., multiqc and cutadapt) are installed in a python 
# virtual environment, you could add the virtual environment in the R path env
# You can add your specific python virtual environment in the Option.txt file

PythonENV.path <- get.value("PythonENV.path")
PythonENV.path

file.exists(PythonENV.path)

# Add python env to this specific project
Sys.setenv(PATH = paste(c(PythonENV.path,
                          Sys.getenv("PATH")),
                          collapse = .Platform$path.sep))


# Step 1: Process Raw -----------------------------------------------------

# FASTQC
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

system2("fastqc", "--help")

# MultiQC
# https://multiqc.info/

system2("multiqc", "--help")

# Cutadapt
# https://cutadapt.readthedocs.io/en/stable/

system2("cutadapt", "--help")


# Step 2: Blast -----------------------------------------------------------

# You need the refere to the folder where the taxonomic data are
NCBI.path <- get.value("NCBI.path")
NCBI.path
file.exists(NCBI.path)

# Blastn
# https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata

system2("blastn",  "-help", 
        env = paste0("BLASTDB=", NCBI.path))

