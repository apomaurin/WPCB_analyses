###
#Title : "Removal of primet"
#Author : Apolline Maurin
# Date : 17/02/2023
###


# 1. R SET UP -------------------------------------------------------------


# 1.1. Install the required packages -----------

install.packages("dada2")
install.packages("plyr")
install.packages("R.utils")
install.packages("ShortRead")

# 1.2. Load the required packages --------------

library(dada2) ; packageVersion("dada2")
library(plyr)
library(R.utils)
library(ShortRead)


# 1.3. Define the working directory ------------

### First define where you're gonna work : 
setwd("C:/Users/apoll/OneDrive/Documents/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/")
### Then, check if you are where you want to be : 
getwd()


# 1.4. Create folders for your output files ----

dir.create("Output_Files_Seq_2_Methodo_Bacteria")

dir.create("C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Output_Files_Seq_2_Methodo_Bacteria/Figures")

dir.create("C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Output_Files_Seq_2_Methodo_Bacteria/Files")

dir.create("C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Output_Files_Seq_2_Methodo_Bacteria/Files/Output_Files")

dir.create("Output_Files_Seq_2_Temp_Bacteria")

dir.create("C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Output_Files_Seq_2_Temp_Bacteria/Figures")

dir.create("C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Output_Files_Seq_2_Temp_Bacteria/Files")

dir.create("C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Output_Files_Seq_2_Temp_Bacteria/Files/Output_Files")



# 2. CUTADAPT -------------------------------------------------------------


# 2.1. Unzip the files -------------------------

# Define the directory :
path <- "C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Data/1_Raw/Bacteria/Temp_Bacteria"
list.files(path)

# Get all the zipped files
zipF <- list.files(path = "C:/Users/apoll/OneDrive/Documents/Th?se/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Data/1_Raw/Bacteria/Temp_Bacteria", pattern = "*.gz", full.names = TRUE)

# Unzip all your files
ldply(.data = zipF, .fun = gunzip, overwrite = TRUE)

# Check if it's alright
list.files(path)


# 2.2. Filtering and Trimming -------------------


# Sorting out the files ----

###there are 2 types of files: Forward (R1) and Reverse (R2) 
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))


# Define the primers WITHOUT Illumina adaptator ----

FWD <- "GTGCCAGCMGCCGCGGTAA" ## Change the sequence to your forward primer sequence
REV <- "GGGACTACHVGGGTWTCTAAT"  ## Change the sequence to your reverse primer sequence


# Definition of the forward and reverse primer complement ----

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather
                               #than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# Sorting and trimming ---- 

### Define the subdirectory for the filtered files : "filtN"
###Then, put them inside
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)


# Primerhits ----
### Count the number of Foward and reverse reads in which the primers are found

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


# Concatenation of the sequences ----
#### In a files : one per Foward and reverse sample

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# 2.3. Cutadapt ---------------------------------

# Define the path to get to cutadapt ----

# Rely on this path to found yours, it should be in a similar place after the downloading

cutadapt <- "C:/Users/apoll/AppData/Local/Programs/Python/Python310/Scripts/cutadapt.exe" # Change it with the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


# Creation of the cutadapt file ----

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


# Cut the primers off with cutadapt ----

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 


# Run Cutadapt ----

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             "--discard-untrimmed", "--minimum-length",length(200),
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


# Concatenation in a file ---- 

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))



