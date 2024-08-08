###
#Title : "Data_Treatment_Methode_Bacteria"
#Author : Apolline Maurin
# Date : 06/03/2023
###


# 1. R SET UP -------------------------------------------------------------


# 1.1. Install the required packages -------------------

###Install packages if you don't have it already. otherwise, don't load it. 

install.packages("phyloseq")
install.packages("Biostrings")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("vegan")

# 1.2. Load the required packages ----------------------

library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(base)

library(BiocGenerics)

# 1.3. Define the working directory --------------------

### First define where you're gonna work : 
setwd("C:/Users/apoll/OneDrive/Documents/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Bacteria")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 1.4. Create folders for your output files ------------

dir.create("./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files")


# 2. Load the data --------------------------------------------------------

count_tab <- read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Output_Files/Seq2_MB_ASVs_counts.tsv", header=T, row.names = 1, check.names=F, sep="\t")
tax_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Output_Files/Seq2_MB_ASV_taxonomy.tsv", header=T, row.names = 1, check.names=F, sep="\t"))
ASV_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Output_Files/Seq2_MB_ASVs.fa", header=F, check.names=F, sep="\t"))
seqtab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Output_Files/Seq2_MB_seqtab.nochim.tsv", header=T, row.names = 1, check.names=F, sep="\t"))



# 3. Prepare the data -----------------------------------------------------

# 3.1. Preparation of the ASV table -------------------

### The table ASV is organised as a FASTA file: 
# ">ASV_X
#  DNA sequence" 

### As we don't want to lose this structure for the analyses to come, but we 
#want to be able to filter at 0.005% and remove any chlorophyll or 
#Eukaryota, we are going to add a column with the corresponding ASV and 
#filter according to it : the title table (or name). 

### It will contain the names of the ASV twice as two line count for the same 
#ASV (the first, the line and the second, the sequence)

# 1/ Create 2 data frame: 

### -> One with the ASV's names
samples.out <- data.frame(A = rownames(count_tab))

### -> One with as many ASV (here 2788) and fill with the number 2
#Why do you ask ? You'll understand in a few step!
df_rep <- rep(2, 2788)

# 2/ Bind the two dataframe together: 
df <- data.frame(cbind("A" = samples.out , "B" = df_rep))

# 3/ create the title table:

### This is the step where you need your table full of 2. We are going to 
#replicate the line of the column A the number of time that is indicate on 
#the column B (i.e. 2 times for each).

title <- data.frame(name = rep(df$A,times=df$B))

# 4/ Create the ASV table containing this new column:
ASV_tab2 <- cbind(title, ASV_tab)

### -> Check, just to be sure! 

# 3.2. Preparation of the seqtab table -------------------

seqtab2 <- cbind(samples.out, t(seqtab))

# 3.3. Filter at 0.005 % -------------------------------

### Here, we want to remove ASVs that are too rare. To do so, we'll only keep 
#the one that are present in 0,005% of the sample 
#(=0,00005% * total Sum of ASVs)

# Sum the ASVs ----

Sum_ASV <- cbind(count_tab, Sum = rowSums(count_tab))
Sum_total <- sum(Sum_ASV$Sum)

# Filter the count table ----

Sum_ASV_Filter <- filter(Sum_ASV, Sum_ASV$Sum >= (0.00005*Sum_total))

# Create a new file without the sum column
Count_tab_filter <- Sum_ASV_Filter %>% select(-(Sum))

# Filter the tax table and the ASV table ----

### The tax and ASV tables doesn't have the ASV count, therefor, we'll use the count 
#table filtered to keep the taxonomy needed.

tax_tab_filter <- tax_tab[rownames(Count_tab_filter),]

### For the ASV table, we can not use the rownames. We will then create a 
#"filter" data frame which will contain the names of the ASV we want to keep.
# Then, we'll keep only those fir the ASV. 
filter <- row.names(Count_tab_filter)
ASV_tab_filter <- ASV_tab2[ASV_tab2$name %in% filter, ]

# Filter the tax table and the seqtab table ----
seq_tab_filter <- seqtab2[seqtab2$name %in% filter, ]


### How the function works : 
# --> tax_tab : the name of your data frame
# --> [] : inside the data frame
# --> [rownames(Count_tab_filter),] : only keep the rows with the following 
#ASVs (the one of interest, after filtering), keep all the columns

### If you wanted to keep some of the column, it's possible to precise it after "," :
# --> Example : [rownames(Count_tab_filter), colnames(patatipatata)]


# 3.4. Remove the unwanted ASVs ------------------------

### Within the ASVs, some are Chloroplast or Eukaryota, which are not of 
#interest. therefor we have to remove them before doing any work on it.

# Remove the chloroplast and eukaryota ----
tax_tab_clean <- subset(tax_tab_filter, !(Order %in% c('Chloroplast' , 'Eukaryota')))
write.table(tax_tab_clean, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Tax_Tab.tsv", sep="\t", quote=F, col.names=NA)

# Check if they have been properly removed ----

Chloroplast <- which(tax_tab_clean == "Chloroplast", arr.ind = TRUE)
Chloroplast

Eucaryote <- which(tax_tab_clean == "Eukaryota", arr.ind = TRUE)
Eucaryote

# Remove the according ASVs from the count table, the ASV table and the seqtab ----

Count_tab_clean <- Count_tab_filter[rownames(tax_tab_clean),]
write.table(Count_tab_clean, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Count_Tab.tsv", sep="\t", quote=F, col.names=NA)

filter2 <- row.names(tax_tab_clean)
ASV_tab_clean <- ASV_tab_filter[ASV_tab_filter$name %in% filter2, ]
### Remove the column used to filter (the on named "name") and save:
ASV_tab_final <- subset(ASV_tab_clean, select = -name)
write.table(ASV_tab_final, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_ASV_Tab.tsv", sep="\t", quote=F, row.names = F, col.names=F)
write.table(ASV_tab_final, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_ASV_Tab.fa", sep="\t", quote=F, row.names = F, col.names=F)

seq_tab_clean <- seq_tab_filter[seqtab2$name %in% filter2, ]
### Remove the column used to filter (the on named "name") and save:
seq_tab_final <- subset(seq_tab_clean, select = -name)
write.table(seq_tab_final, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_seqtab.tsv", sep="\t", quote=F, row.names = T, col.names=NA)
write.table(seq_tab_clean, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_seqtab_name.tsv", sep="\t", quote=F, row.names = T, col.names=NA)

# 4. Create a data frame with all your information and count --------------

# 4.1. Create a table of information on your samples ---

### The info table will contain your variables according to the name of your sample.
# --> For each, I indicate, the number of tract and which tree. The 2 
#variables need to be of the same size as the subject (= your sample)
# --> Here, my variables are the number of tracts and the tree.

samples.out <- rownames(t(Count_tab_clean))
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
tract <- c("2_tracts", "3_tracts", "4_tracts","2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts", "2_tracts", "2_tracts","2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts","2_tracts", "2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts" ,"Control","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery")
tree <- c("Tree_1","Tree_1", "Tree_1", "Tree_2","Tree_2","Tree_2", "Tree_3","Tree_3","Tree_3","Tree_4","Tree_4","Tree_4","Tree_5","Tree_6","Tree_7","Tree_7","Tree_7","Tree_8","Tree_8","Tree_8","Tree_9","Tree_11","Tree_11","Tree_11","Tree_12", "Tree_12", "Tree_12", "Temoin", "Tree_1", "Tree_2","Tree_3","Tree_4","Tree_5","Tree_6","Tree_7","Tree_8","Tree_9","Tree_11","Tree_12")
status <- c("Tract", "Tract", "Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract", "Tract", "Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract", "Tract","Tract","Tract","Tract","Tract","Tract" ,"Control","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery")
sample_info <- data.frame(Subject=subject, Tract=tract, Tree=tree, Status=status)
head(sample_info)
write.table(sample_info, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Sample_info.tsv", sep="\t", quote=F, col.names=NA)
rownames(sample_info) <- samples.out


# 4.2. Create the data frame ---------------------------

### First, turn the row name into a column (which will be the common column 
#for the next step):
count_tabr <- t(Count_tab_clean)
count_tabr_subj_tot <- data.frame(Subject = rownames(count_tabr), count_tabr, row.names=NULL)

### then, create the data.frame:
Info_count_tot <- full_join(sample_info, count_tabr_subj_tot, by="Subject")

###If you want to remove the common column, use this:
#df = data.frame(Info_count[,-1], row.names = sample_info$Subject)

### Finally, reverse again to get the ASV on the column:
final_table_tot <- t(Info_count_tot)

### Save the file:
write.table(final_table_tot, file = "./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Final_Table.tsv", sep="\t", quote=F, col.names = FALSE)

