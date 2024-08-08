###
#Title : "Stats_Alpha_Div_Bacteria"
#Author : Apolline Maurin
# Date : 16/03/2023
###

# 1. R SET UP -------------------------------------------------------------


# 1.1. Install the required packages -------------------

###Install packages if you don't have it already. otherwise, don't load it. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")


BiocManager::install("phyloseq")
install.packages("Biostrings")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("vegan")
install.packages("dendextend")
install.packages('devtools')
BiocManager::install("DECIPHER")
install.packages("phangorn")
install.packages("ggpubr")
BiocManager::install("ANCOMBC", force = TRUE)
install.packages("DT")


# 1.2. Load the required packages ----------------------

library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(dendextend) ;packageVersion("dendextend")
library(devtools)
library(DECIPHER)
library(dada2) ; packageVersion("dada2")
library(phangorn)
library(ggpubr)

# 1.3. Define the working directory --------------------

### First define where you're gonna work : 
setwd("D:/Thèse/7_Labo/Sequencage/Seq_3_Pretest_Objectif3")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 2. Load the data --------------------------------------------------------

count_tab <- as.data.frame(read.table("./Output_Files_Seq_3_PTB/Files/Stats_Files/Seq3_PTB_Count_Tab.tsv", header=T, row.names = 1, check.names=F, sep="\t"))
tax_tab <- as.matrix(read.table("./Output_Files_Seq_3_PTB/Files/Stats_Files/Seq3_PTB_Tax_Tab.tsv", header=T, row.names = 1, check.names=F, sep="\t"))
ASV_tab <- as.data.frame(read.table("./Output_Files_Seq_3_PTB/Files/Stats_Files/Seq3_PTB_ASV_Tab.fa", header=T, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("./Output_Files_Seq_3_PTB/Files/Stats_Files/Seq3_PTB_Sample_info.tsv", header=T, row.names = 1, check.names=F, sep="\t"))
seq_tab <- as.data.frame(read.table("./Output_Files_Seq_3_PTB/Files/Stats_Files/Seq3_PTB_seqtab.tsv", header=T,row.names=1, check.names=F, sep="\t"))

# 3. Transforming data (Hellinger) --------------------------

count_tab_ATB_SW <- select(count_tab, -contains("_N"))
seq_tab_ATB_SW <- select(seq_tab, -contains("_N"))
sample_info_ATB_SW <- subset(sample_info, !(Treatment %in% c('N')))

count_tab_ATB_N <- select(count_tab, -contains("_SW"))
seq_tab_ATB_N <- select(seq_tab, -contains("_SW"))
sample_info_ATB_N <- subset(sample_info, !(Treatment %in% c('SW')))

count_tab_SW_N <- select(count_tab, -contains("_ATB"))
seq_tab_SW_N <- select(seq_tab, -contains("_ATB"))
sample_info_SW_N <- subset(sample_info, !(Treatment %in% c('ATB')))

### All
data_stand <- decostand(count_tab, method = "hellinger")
seqtab <- decostand(seq_tab, method = "hellinger")

data_stand_ATB_N <- decostand(count_tab_ATB_N, method = "hellinger")
seqtab_ATB_N <- decostand(seq_tab, method = "hellinger")

data_stand_ATB_SW <- decostand(count_tab_ATB_SW, method = "hellinger")
seqtab_ATB_SW <- decostand(seq_tab, method = "hellinger")

data_stand_SW_N <- decostand(count_tab_SW_N, method = "hellinger")
seqtab_SW_N <- decostand(seq_tab, method = "hellinger")

# 4. Creating a phylogenetic tree (UNIFRAC) -----

### To calculate a distance which the phylogenetic distance is taking into 
#account  we use the UNIFRAC distance. To do so, we need to do a phylogenetic
#tree with the fasta file.   

# 3.3.1. Alignement ----

### First, we read the fasta file as DNA strings:
#Reminder: the file need to start with ">". Make sure that when you 
#registered your file, you do it so the name are ">ASV"

fastaFile <- readDNAStringSet("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_files/Seq2_MB_ASV_tab.fa")
fastaFile

### Then, we align the sequences:
seqs <- getSequences(fastaFile)
names(fastaFile) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(fastaFile), anchor=NA)

# 3.3.2. Tree construction ----

### Transformation from the DNA format to the phyDat format used to construct 
#the phylogenetic tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

### Computing the distances under different substitution models. 
#this will quantify th dissimilarity or evolutionary divergence between 
#pair of sequences. 
dm <- dist.ml(phang.align)

### Construction of the tree, using the neighbor-joining algorithm. 
#It's an iterative algorithm that builds a tree by joining pairs of 
#sequences or taxa based on their pairwise distances.
treeNJ <- NJ(dm) # Note, tip order != sequence order

### Then, we computes the likelihood of a phylogenetic tree given a sequence 
#alignment and a model
fit = pml(treeNJ, data=phang.align) #negative edges length changed to 0!

### Fitting of the phylogenetic tree using a Generalized time reversible with 
#Gama rate variation (GTR+G+I) maximum likelihood tree using the 
#neighbor-joining tree as a starting point 
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
#plot(fitGTR)


# 5. Creating a Phyloseq Object -----------------

### Phyloseq object all
samples.out <- rownames(t(data_stand))
rownames(sample_info) <- samples.out

ps_trans_all <- phyloseq(otu_table(data_stand, taxa_are_rows=TRUE), 
                         sample_data(sample_info), 
                         tax_table(tax_tab))

samples.out1 <- rownames(t(data_stand_ATB_N))
rownames(sample_info_ATB_N) <- samples.out1

ps_trans_ATB_N <- phyloseq(otu_table(data_stand_ATB_N, taxa_are_rows=TRUE), 
                           sample_data(sample_info_ATB_N), 
                           tax_table(tax_tab))

samples.out2 <- rownames(t(data_stand_ATB_SW))
rownames(sample_info_ATB_SW) <- samples.out2

ps_trans_ATB_SW <- phyloseq(otu_table(data_stand_ATB_SW, taxa_are_rows=TRUE), 
                            sample_data(sample_info_ATB_SW), 
                            tax_table(tax_tab))

samples.out3 <- rownames(t(data_stand_SW_N))
rownames(sample_info_SW_N) <- samples.out3

ps_trans_SW_N <- phyloseq(otu_table(data_stand_SW_N, taxa_are_rows=TRUE), 
                          sample_data(sample_info_SW_N), 
                          tax_table(tax_tab))


### Phyloseq object all (UNIFRAC)
samples.out3 <- rownames(t(seq_tab))
sample_info_UNIFRAC <- sample_info
rownames(sample_info_UNIFRAC) <- samples.out3

samples.out4 <- rownames(seq_tab)
taxa_tab_UNIFRAC <- taxa_tab
rownames(taxa_tab_UNIFRAC) <- samples.out4

ps_UNIFRAC_all <- phyloseq(otu_table(seqtab_all, taxa_are_rows=TRUE), 
                           sample_data(sample_info_UNIFRAC), 
                           tax_table(taxa_tab_UNIFRAC),
                           phy_tree(fitGTR$tree))

# 6. Distance calculation -----------------------

# 6.1. Bray-Curtis ----------------------------------------------------

bc_rel_asv_all <- data.frame(phyloseq::otu_table(ps_trans_all))
bc_rel_asv_all <- t(bc_rel_asv_all)
bc_dist_all <- as.matrix(vegan::vegdist(bc_rel_asv_all, method = "bray")) 
as.matrix(bc_dist_all)[1:5, 1:5]  

bc_rel_asv_ATB_N <- data.frame(phyloseq::otu_table(ps_trans_ATB_N))
bc_rel_asv_ATB_N <- t(bc_rel_asv_ATB_N)
bc_dist_ATB_N <- as.matrix(vegan::vegdist(bc_rel_asv_ATB_N, method = "bray")) 
as.matrix(bc_dist_ATB_N)[1:5, 1:5]  

bc_rel_asv_ATB_SW <- data.frame(phyloseq::otu_table(ps_trans_ATB_SW))
bc_rel_asv_ATB_SW <- t(bc_rel_asv_ATB_SW)
bc_dist_ATB_SW <- as.matrix(vegan::vegdist(bc_rel_asv_ATB_SW, method = "bray")) 
as.matrix(bc_dist_ATB_SW)[1:5, 1:5]  

bc_rel_asv_SW_N <- data.frame(phyloseq::otu_table(ps_trans_SW_N))
bc_rel_asv_SW_N <- t(bc_rel_asv_SW_N)
bc_dist_SW_N <- as.matrix(vegan::vegdist(bc_rel_asv_SW_N, method = "bray")) 
as.matrix(bc_dist_SW_N)[1:5, 1:5]  

# 6.2. Hellinger ------------------------------------------------------

eucli_rel_asv_all <- data.frame(phyloseq::otu_table(ps_trans_all))
eucli_rel_asv_all <- t(eucli_rel_asv_all)
eucli_dist_all <- as.matrix(vegan::vegdist(eucli_rel_asv_all, method = "euclidean")) 
as.matrix(eucli_dist_all)[1:5, 1:5]

eucli_rel_asv_ATB_N <- data.frame(phyloseq::otu_table(ps_trans_ATB_N))
eucli_rel_asv_ATB_N <- t(eucli_rel_asv_ATB_N)
eucli_dist_ATB_N <- as.matrix(vegan::vegdist(bc_rel_asv_ATB_N, method = "bray")) 
as.matrix(bc_dist_ATB_N)[1:5, 1:5]  

eucli_rel_asv_ATB_SW <- data.frame(phyloseq::otu_table(ps_trans_ATB_SW))
eucli_rel_asv_ATB_SW <- t(eucli_rel_asv_ATB_SW)
eucli_dist_ATB_SW <- as.matrix(vegan::vegdist(bc_rel_asv_ATB_SW, method = "bray")) 
as.matrix(bc_dist_ATB_SW)[1:5, 1:5]  

eucli_rel_asv_SW_N <- data.frame(phyloseq::otu_table(ps_trans_SW_N))
eucli_rel_asv_SW_N <- t(bc_rel_asv_SW_N)
eucli_dist_SW_N <- as.matrix(vegan::vegdist(bc_rel_asv_SW_N, method = "bray")) 
as.matrix(bc_dist_SW_N)[1:5, 1:5]  

# 6.3. UNIFRAC --------------------------------------------------------

# Weighted ---- 
### = takes into account the relative abundance of species/taxa shared between samples
### Note : it will randomly assign the root of the phylogenetic tree
wunifrac_dis_all <- as.matrix(UniFrac(ps_UNIFRAC_all, weighted = TRUE, normalized = TRUE))
wunifrac_dis_tract <- as.matrix(UniFrac(ps_UNIFRAC_tract, weighted = TRUE, normalized = TRUE))
wunifrac_dis_4t <- as.matrix(UniFrac(ps_UNIFRAC_4t, weighted = TRUE, normalized = TRUE))
wunifrac_dis_3t <- as.matrix(UniFrac(ps_UNIFRAC_3t, weighted = TRUE, normalized = TRUE))
wunifrac_dis_2t <- as.matrix(UniFrac(ps_UNIFRAC_2t, weighted = TRUE, normalized = TRUE))

# Unweighted ----
### = only considers presence/absence
### Note : it will randomly assign the root of the phylogenetic tree
unifrac_dis_all <- as.matrix(UniFrac(ps_UNIFRAC_all, weighted = FALSE, normalized = TRUE))
unifrac_dis_tract <- as.matrix(UniFrac(ps_UNIFRAC_tract,  weighted = FALSE, normalized = TRUE))


# 7. PCoA ---------------------------------------

# 7.1. Bray-Curtis ----------------------------------------------------

### All
theme_set(theme_bw())
ordu = ordinate(ps_trans_ATB_SW, "PCoA", "bray")
plot_ordination(ps_trans_ATB_SW, ordu)
bc_all = plot_ordination(ps_trans_ATB_SW, ordu, color="Temperature", shape="Treatment", axes=1:2)
bc_all_geom = bc_all + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  stat_ellipse(aes(group = Treatment), linetype = 2)
bc_all_geom


# 7.2. Hellinger ------------------------------------------------------

### All
theme_set(theme_bw())
ordu3 = ordinate(ps_trans_all, "PCoA", "euclidean")
plot_ordination(ps_trans_all, ordu3)
eucli_all = plot_ordination(ps_trans_all, ordu3, color="Temperature", shape="Treatment", axes=1:2)
eucli_all_geom = eucli_all + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Treatment), linetype = 2)
eucli_all_geom


# 7.3. UNIFRAC --------------------------------------------------------

# Weighted ---- 
### All
theme_set(theme_bw())
ordu5 = ordinate(ps_UNIFRAC_all, "PCoA", "wunifrac")
plot_ordination(ps_UNIFRAC_all, ordu5)
wunifrac_all = plot_ordination(ps_UNIFRAC_all, ordu5, color="Tract", shape="Tract", axes=1:2)
wunifrac_all_geom = wunifrac_all + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
wunifrac_all_geom


### Two graph in one
ggarrange(wunifrac_all_geom, wunifrac_tract_geom + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Unweighted ---- 
### All
theme_set(theme_bw())
ordu7 = ordinate(ps_UNIFRAC_all, "PCoA", "unifrac")
plot_ordination(ps_UNIFRAC_all, ordu7)
unifrac_all = plot_ordination(ps_UNIFRAC_all, ordu7, color="Tract", shape="Bloc", axes=1:2)
unifrac_all_geom = unifrac_all + geom_point(size=5) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
unifrac_all_geom


### All distances - All
ggarrange(bc_all_geom, eucli_all_geom, wunifrac_all_geom, unifrac_all_geom + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "right")


### The selected distances 
#All
final_plot <- ggarrange(bc_all_geom, wunifrac_all_geom, 
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1,
                        common.legend = TRUE, legend = "right")
final_plot

ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/PCoA_bacteria.jpeg", plot = final_plot, width = 12, height = 6, units = "in", dpi = 300)


# 8. Permanova ----------------------------------

# 8.1. Bray-Curtis ----------------------------------------------------
perma_bc_all <-adonis2(bc_rel_asv_all~Treatment+Temperature, data=sample_info, permutations=9999)
perma_bc_all
perma_bc_all2 <- mutate(perma_bc_all, Distance = "BC_all")

perma_bc_ATB_N <-adonis2(bc_dist_ATB_N~Treatment, data=sample_info_ATB_N, permutations=9999)
perma_bc_ATB_N
perma_bc_ATB_N2 <- mutate(perma_bc_ATB_N, Distance = "BC_ATB_N")

perma_bc_ATB_SW <-adonis2(bc_dist_ATB_SW~Treatment+Temperature, data=sample_info_ATB_SW, permutations=9999)
perma_bc_ATB_SW
perma_bc_ATB_SW2 <- mutate(perma_bc_ATB_SW, Distance = "BC_ATB_SW")

perma_bc_SW_N <-adonis2(bc_dist_SW_N~Treatment+Temperature, data=sample_info_SW_N, permutations=9999)
perma_bc_SW_N
perma_bc_SW_N2 <- mutate(perma_bc_SW_N, Distance = "BC_SW_N")

# 8.2. Hellinger ------------------------------------------------------

perma_eucli_all <-adonis2(eucli_dist_ATB_N~Treatment+Temperature, data=sample_info_ATB_N, permutations=9999)
perma_eucli_all
perma_eucli_all2 <- mutate(perma_eucli_all, Distance = "Eucli_all")

# 8.3. Weighted UNIFRAC -----------------------------------------------

perma_wuni_all <- adonis2(wunifrac_dis_all~Bloc+Tract+Tree+Status, data=sample_info_UNIFRAC, permutations=9999)
perma_wuni_all
perma_wuni_all2 <- mutate(perma_wuni_all, Distance = "Wunifrac_all")

# 8.4. Unweighted UNIFRAC  --------------------------------------------

perma_uni_all <-adonis2(unifrac_dis_all~Bloc+Tract+Tree, data=sample_info_UNIFRAC, permutations=9999)
perma_uni_all
perma_uni_all2 <- mutate(perma_uni_all, Distance = "Unifrac_all")

# 8.5. Creating a table with all the permanova  -----------------------

### Create the table
perma_all <- bind_rows(perma_bc_all2, perma_bc_tract2, perma_eucli_all2, 
                       perma_eucli_tract2, perma_wuni_all2, perma_wuni_tract2, 
                       perma_uni_all2, perma_uni_tract2)

write.table(perma_all, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_permanova_tot.tsv", sep="\t", quote=F, col.names=NA)


### Only keep the p.value and the distance name and remove the 
#Total and residual line
perma_pr <- select(perma_all, "Distance", p.value = "Pr(>F)")
perma_pr2 <- subset(perma_pr, !grepl(c("Total|Residual"), rownames(perma_pr)))

### Create a new column to store the Tract or Tree information
perma_pr2$Group <- rownames(perma_pr2)

### Remove the "..." and the number from the row names
perma_pr2$Group <- gsub("\\.\\.\\.\\d+", "", perma_pr2$Group)

### Pivot the table to reshape it
reshaped_perma <- pivot_wider(perma_pr2, names_from = Distance, values_from = p.value)

### Set the row names to "Tract" and "Tree"
perma_name <- reshaped_perma$Group
# or: perma_name <- c("Tract", "Tree")

### Delete the Group column 
reshaped_perma$Group <- NULL

### Set the row names to "Tract" and "Tree"
rownames(reshaped_perma) <- perma_name

write.table(reshaped_perma, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_permanova.tsv", sep="\t", quote=F, col.names=NA)

## GLobalement, différence ATB et SW, pas de différence ATB et rien, mais aussi effet SW et N, donc je ne sais aps trop quoi conclure... 