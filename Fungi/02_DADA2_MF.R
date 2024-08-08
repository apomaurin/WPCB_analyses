###
#Title : "DADA2_Methode_Fungi"
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
library(stats)


# 1.3. Define the working directory ------------

### First define where you're gonna work : 
setwd("C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi")
### Then, check if you are where you want to be : 
getwd()

# 2- DADA2 ----------------------------------------------------------------


# 2.1. Define the directory -------------------------

###Put the directory containing the Fastq files from cutadapt
path <- "C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Data/Raw/cutadapt"
list.files(path)


# 2.2. Manipulation of the Fastq files --------------
### the objective is to matched the list of the foward Fastq files with the reverse one

cutFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

# Extract sample names : 
cutFs.name=sub("S503.","x.",cutFs) #Remplace the barcoded primer ID by "x." to 
#manage the sample names
cutFs.name=sub("S505.","x.",cutFs.name)
cutFs.name=sub("S506.","x.",cutFs.name)
cutFs.name=sub("S507.","x.",cutFs.name)
cutFs.name=sub("S510.","x.",cutFs.name)
cutFs.name=sub("S511.","x.",cutFs.name)
cutFs.name=sub("S513.","x.",cutFs.name)


get.sample.name <- function(fname) strsplit(basename(fname), "x.")[[1]][2]
sample.names <- unname(sapply(cutFs.name, get.sample.name), 
                       sapply(cutRs.name, get.sample.name))
head(sample.names)


# 2.3. Quality profile ------------------------------

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:4])

### What to do with it: ###

#### Legend:  
# --> In gray: the heat map of the frequency of each quality score at each 
#base point  
# --> In green: the mean quality score at each position 
# --> In orange: the quartiles of the quality distribution 
# --> Red line: the scale of reads that extend to at least that position  

#### Interpretation: ###
# --> Flat line means same length 
# --> If the red line fall, it means there is a difference in length  

#### The advice is to truncate the last 10 nucletoides to avoid less 
#well-controlled errors that can arise there  

### Keep in mind: ###
# --> That the quality is often worse for the reverse reads, it is common in 
     #Illumina sequencing  
# --> DADA2 incorporates quality information into its error model which makes 
     #the algorithm robust to lower quality sequence 


# 2.4. Filter & Trim --------------------------------

# Assign the filenames for the filtered fastq.gz files ----

path <- "C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Data/Clean"

filtFs <- file.path(path, "Methode_Fungi_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "Methode_Fungi_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Define the filtering parameters and filter ---- 

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)



### Good to know: ### 

#### It is possible to change the filtering parameters: 
# --> If too few reads are passing the filter: relax maxEE (ex: maxEE=c(2,5)) 
# --> To remove low quality tail: reduce truncLen 

#### It is not recommending to truncate reads when using ITS because of the 
    #large length variation


# 2.5. Learn the error rates ------------------------ 

# Unzipp the files ----

# Define the directory :
path <- "C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Data/Clean/Methode_Fungi_filtered"

# Get all the zipped files
zipF <- list.files(path = "C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Data/Clean/Methode_Fungi_filtered",
                   pattern = "*.gz", full.names = TRUE)

# Unzip all your files
ldply(.data = zipF, .fun = gunzip, overwrite = TRUE)

# Check if it's alright
list.files(path)

# DADA2 algorithm : parametric error model (err)  

### LearnErrors : learns error model from the data by alternating estimation of
#the error rates and inference of sample composition until they converge on 
#a jointly consistent solution

### Tool to visualize the frequency of error rate as a function of quality score 

filtFs <- file.path(path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, paste0(sample.names, "_R_filt.fastq"))

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)


plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


### What to do with it: ###

### These plots show the error rates for each possible transition (A->C, A->G, etc.)  

### Legend: 
# --> Points = observed error rates for each consensus quality score 
# --> Black lines = estimate error after convergence of the matrice-learning algorithm 
# --> Red lines = shows the error rates expected under the nominal def of the Q-score  

### Basically, as long as points and black line are a fit + error rates drop, it's good ! 

# 2.5. Dereplication --------------------------------

### Dereplication combines all identical sequencing reads into into a unique 
#sequences with a corresponding abundance equal to the number of reads with 
#that unique sequence. 

### Dereplication substantially reduces computation time by eliminating 
#redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# 2.6 Sample Inference ------------------------------
### We are now ready to apply the core sample inference algorithm to the filtered
#and trimmed sequence data.

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

# Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

#DADA2 algorithm inferred x true sequence variants from x unique sequences in the first sample
#help("dada-class") #if needed


# 2.7. Merge paired reads ---------------------------

# Mergeing of the foward and reverse reads together to obtain a full denoised sequences 
#--> Merging is performed by aligning the denoised F and R reads and construct 
#a merged "contig" sequences 
# Merdge only if at least 12 bases overlap and are identical to each other in 
#the overlap region That can be changed via function argument if wanted

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# The mergers object is a list of data.frame from each sample 
# ach data.frame contains : $sequence, $abundance, indices of the $foward and 
#$reverse sequence variants that were merged 



# 2.8. Construct sequence table ---------------------

# Construction of an amplicon sequence variant table (ASV) table ----

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths ----

table(nchar(getSequences(seqtab)))

# This sequence table = matrix 
#--> Rows = samples 
#--> Columns = sequence variants  

# If sequences are much longer or shorter than expected
#--> It may be the r! of non-specific priming 
#--> Possibility to remove non-target-length sequences from your sequence table : 

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
#dim(seqtab2)
#table(nchar(getSequences(seqtab2))



# 2.9. Remove chimeras -----------------------------

# The core dada method corrects substitution and indel errors but chimeras remain.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
write.table(seqtab.nochim, file="C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_seqtab.nochim.tsv", 
            sep="\t", quote=F, col.names=NA)


# 2.10. Track reads through pipeline ----------------

# As a final check of our progress, we'll look at the number of reads that made 
#it through each step in the pipeline: 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace 
#sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


write.table(track, file="C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_Summary.txt")


# 2.11. Assign taxonomy -----------------------------

# Preparation ----

#Download on zenodo.org the file with the format: 
#- silva_nr99_v138.1_train_set.fa.gz 
#- silva_species_assignment_v138.1.fa.gz 
#--> https://benjjneb.github.io/dada2/training.html  


#Place them in the directory with the fastq files Filtered and extract it  

# Assign Taxonomy and species ----

path <- "C:/Users/Léo-Paul/Desktop/Apo/Seq_2/Resultats_Seq_2/Methode_Fungi/Data/Clean/Methode_Fungi_filtered"
list.files(path)

#- Taxon: 
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Léo-Paul/Desktop/Apo/Taxonomic_Assignation/sh_general_release_dynamic_16.10.2022.fasta",
                       multithread=FALSE, tryRC = TRUE)

getwd()
write.table(taxa, file="./Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_Taxa_solo.tsv",
            sep="\t", quote=F, col.names=NA)

# Inspection of the taxonomic/species assignments ----

taxa.print <- taxa 
rownames(taxa.print) <- NULL # Removing sequence rownames for display only
head(taxa.print)


# 2.12. Standard output files ----------------------

# This code has been modify from Thibault's one (dada2_Thibault.R) # 


# Giving our seq headers more manageable names (ASV_1, ASV_2...) ----

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV_Seq2B", i, sep="_") #Change the name !
}


# Output files ----  

# Fasta:

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta,"./Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_ASVs.fa") #Change the name !


# Count table:

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "./Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_ASVs_counts.tsv", sep="\t", quote=F, 
            col.names=NA)#Change the name !


# Tax table:

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, file="./Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_ASV_taxonomy.tsv", sep="\t", quote=F, 
            col.names=NA)#Change the name !


# Export results in a matrix of samples (rows) and AVS counts (columns), in descending order:

ASV.sample <- seqtab.nochim
colnames(ASV.sample) <- sub(">", "", asv_headers)

ASV.sample.sort <- ASV.sample
for (i in ncol(ASV.sample.sort):1) {
  ASV.sample.sort <- ASV.sample.sort[order(-ASV.sample.sort[,i]),]
}

# write.table(ASV.sample.sort, "Seq2B_ASV_samples.txt", sep="\t", quote=F):  
# Change the name !

write.table(ASV.sample.sort, file="./Output_Files_Seq_2_Methode_Fungi/Files/Output_Files/Seq2_MF_ASV_samples.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !


