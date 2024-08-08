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

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

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
setwd("C:/Users/apoll/OneDrive/Documents/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Fungi")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 2. Load the data --------------------------------------------------------

count_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Count_Tab_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.matrix(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Tax_Tab.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Sample_info_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))
ASV_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_ASV_Tab.fa", header=F, check.names=F, sep="\t"))
seq_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_seqtab_No_Conta.tsv", header=T,row.names=1, check.names=F, sep="\t"))

# 3. Transforming data --------------------------

# 3.1. Create table with only the tract -------------------------------

### Create table with only the tracts
count_tab_tract <- select(count_tab, -contains("_G"))
write.table(count_tab_tract, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Count_Tab_Tract_Only.tsv", sep="\t", quote=F, col.names=NA)

seq_tab_tract <- select(seq_tab, -contains("_G"))
write.table(seq_tab_tract, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Seqtab_Tract_Only.tsv", sep="\t", quote=F, col.names=NA)

sample_info_tract <- subset(sample_info, !(Status %in% c('Gallery')))
write.table(sample_info_tract, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Sample_Info_Tract_Only.tsv", sep="\t", quote=F, col.names=NA)

### Create table with only 2t, 3t or 4t and Galleries

#### 4t
count_tab_4t <- select(count_tab, -contains("2t_champignons"))
count_tab_4t <- select(count_tab_4t, -contains("3t_champignons"))

seq_tab_4t <- select(seq_tab, -contains("2t_champignons"))
seq_tab_4t <- select(seq_tab_4t, -contains("3t_champignons"))

sample_info_4t <- subset(sample_info, !(Tract %in% c('2_tracts','3_tracts')))

#### 3t
count_tab_3t <- select(count_tab, -contains("4t_champignons"))
count_tab_3t <- select(count_tab_3t, -contains("2t_champignons"))

seq_tab_3t <- select(seq_tab, -contains("2t_champignons"))
seq_tab_3t <- select(seq_tab_3t, -contains("4t_champignons"))

sample_info_3t <- subset(sample_info, !(Tract %in% c('2_tracts','4_tracts')))

#### 2t
count_tab_2t <- select(count_tab, -contains("3t_champignons"))
count_tab_2t <- select(count_tab_2t, -contains("4t_champignons"))

seq_tab_2t <- select(seq_tab, -contains("4t_champignons"))
seq_tab_2t <- select(seq_tab_2t, -contains("3t_champignons"))

sample_info_2t <- subset(sample_info, !(Tract %in% c('4_tracts','3_tracts')))

# 3.2. Data transformation --------------------------------------------

# 3.2.1. Hellinger transformation ----

## All
data_stand_all <- decostand(count_tab, method = "hellinger")
seqtab_all <- decostand(seq_tab, method = "hellinger")

### Only tracts
data_stand_tract <- decostand(count_tab_tract, method = "hellinger")
seqtab_tract <- decostand(seq_tab_tract, method = "hellinger")

### 4t
data_stand_4t <- decostand(count_tab_4t, method = "hellinger")
seqtab_4t <- decostand(seq_tab_4t, method = "hellinger")

### 3t
data_stand_3t <- decostand(count_tab_3t, method = "hellinger")
seqtab_3t <- decostand(seq_tab_3t, method = "hellinger")

### 2t
data_stand_2t <- decostand(count_tab_2t, method = "hellinger")
seqtab_2t <- decostand(seq_tab_2t, method = "hellinger")


# 4. Creating a phylogenetic tree (UNIFRAC) -----

### To calculate a distance which the phylogenetic distance is taking into 
#account we use the UNIFRAC distance. To do so, we need to do a phylogenetic
#tree with the fasta file.   

# 3.3.1. Alignement ----

### First, we read the fasta file as DNA strings:
#Reminder: the file need to start with ">". Make sure that when you 
#registered your file, you do it so the name are ">ASV"

fastaFile <- readDNAStringSet("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/Seq2_MF_ASV_tab.fa")
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
samples.out <- rownames(t(data_stand_all))
rownames(sample_info) <- samples.out

ps_trans_all <- phyloseq(otu_table(data_stand_all, taxa_are_rows=TRUE), 
                         sample_data(sample_info), 
                         tax_table(taxa_tab))

### Phyloseq object tract
samples.out2 <- rownames(t(data_stand_tract))
rownames(sample_info_tract) <- samples.out2

ps_trans_tract <- phyloseq(otu_table(data_stand_tract, taxa_are_rows=TRUE), 
                           sample_data(sample_info_tract), 
                           tax_table(taxa_tab))

### Phyloseq object all (UNIFRAC)
samples.out3 <- rownames(t(seq_tab))
samples.out3
sample_info_UNIFRAC <- sample_info
rownames(sample_info_UNIFRAC) <- samples.out3

samples.out4 <- rownames(seq_tab)
taxa_tab_UNIFRAC <- taxa_tab
rownames(taxa_tab_UNIFRAC) <- samples.out4

ps_UNIFRAC_all <- phyloseq(otu_table(seqtab_all, taxa_are_rows=TRUE), 
                           sample_data(sample_info_UNIFRAC), 
                           tax_table(taxa_tab_UNIFRAC),
                           phy_tree(fitGTR$tree))

### Phyloseq object tract (UNIFRAC)
samples.out5 <- rownames(t(seq_tab_tract))
sample_info_tract_UNIFRAC <- sample_info_tract
rownames(sample_info_tract_UNIFRAC) <- samples.out5

ps_UNIFRAC_tract <- phyloseq(otu_table(seqtab_tract, taxa_are_rows=TRUE), 
                             sample_data(sample_info_tract_UNIFRAC), 
                             tax_table(taxa_tab_UNIFRAC),
                             phy_tree(fitGTR$tree))

samples.out6 <- rownames(t(seq_tab_4t))
sample_info_4t_UNIFRAC <- sample_info_4t
rownames(sample_info_4t_UNIFRAC) <- samples.out6

ps_UNIFRAC_4t <- phyloseq(otu_table(seqtab_4t, taxa_are_rows=TRUE), 
                          sample_data(sample_info_4t_UNIFRAC), 
                          tax_table(taxa_tab_UNIFRAC),
                          phy_tree(fitGTR$tree))

samples.out7 <- rownames(t(seq_tab_3t))
sample_info_3t_UNIFRAC <- sample_info_3t
rownames(sample_info_3t_UNIFRAC) <- samples.out7

ps_UNIFRAC_3t <- phyloseq(otu_table(seqtab_3t, taxa_are_rows=TRUE), 
                          sample_data(sample_info_3t_UNIFRAC), 
                          tax_table(taxa_tab_UNIFRAC),
                          phy_tree(fitGTR$tree))

samples.out8 <- rownames(t(seq_tab_2t))
sample_info_2t_UNIFRAC <- sample_info_2t
rownames(sample_info_2t_UNIFRAC) <- samples.out8

ps_UNIFRAC_2t <- phyloseq(otu_table(seqtab_2t, taxa_are_rows=TRUE), 
                          sample_data(sample_info_2t_UNIFRAC), 
                          tax_table(taxa_tab_UNIFRAC),
                          phy_tree(fitGTR$tree))

# 6. Distance calculation -----------------------

# 6.1. Bray-Curtis ----------------------------------------------------

bc_rel_asv_all <- data.frame(phyloseq::otu_table(ps_trans_all))
bc_rel_asv_all <- t(bc_rel_asv_all)
bc_dist_all <- as.matrix(vegan::vegdist(bc_rel_asv_all, method = "bray")) 
as.matrix(bc_dist_all)[1:5, 1:5]  

bc_rel_asv_tract <- data.frame(phyloseq::otu_table(ps_trans_tract))
bc_rel_asv_tract <- t(bc_rel_asv_tract)
bc_dist_tract <- as.matrix(vegan::vegdist(bc_rel_asv_tract, method = "bray")) 
as.matrix(bc_dist_tract)[1:5, 1:5]  

bc_rel_asv_4t <- t(count_tab_4t)
bc_dist_4t <- as.matrix(vegan::vegdist(bc_rel_asv_4t, method = "bray")) 
as.matrix(bc_dist_4t)[1:5, 1:5]  

bc_rel_asv_3t <- t(count_tab_3t)
bc_dist_3t <- as.matrix(vegan::vegdist(bc_rel_asv_3t, method = "bray")) 
as.matrix(bc_dist_3t)[1:5, 1:5]  

bc_rel_asv_2t <- t(count_tab_2t)
bc_dist_2t <- as.matrix(vegan::vegdist(bc_rel_asv_2t, method = "bray")) 
as.matrix(bc_dist_2t)[1:5, 1:5] 

# 6.2. Hellinger ------------------------------------------------------

eucli_rel_asv_all <- data.frame(phyloseq::otu_table(ps_trans_all))
eucli_rel_asv_all <- t(eucli_rel_asv_all)
eucli_dist_all <- as.matrix(vegan::vegdist(eucli_rel_asv_all, method = "euclidean")) 
as.matrix(eucli_dist_all)[1:5, 1:5]

eucli_rel_asv_tract <- data.frame(phyloseq::otu_table(ps_trans_tract))
eucli_rel_asv_tract <- t(eucli_rel_asv_tract)
eucli_dist_tract <- as.matrix(vegan::vegdist(eucli_rel_asv_tract, method = "euclidean")) 
as.matrix(eucli_dist_tract)[1:5, 1:5]

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
ordu = ordinate(ps_trans_all, "PCoA", "bray")
plot_ordination(ps_trans_all, ordu)
bc_all = plot_ordination(ps_trans_all, ordu, color="Tract", shape="Tract", axes=1:2)
bc_all_geom = bc_all + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)

### Tract
theme_set(theme_bw())
ordu2 = ordinate(ps_trans_tract, "PCoA", "bray")
plot_ordination(ps_trans_tract, ordu2)
bc_tract = plot_ordination(ps_trans_tract, ordu2, color="Tract", shape="Tract", axes=1:2)
bc_tract_geom = bc_tract + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)

### Two graph in one
ggarrange(bc_all_geom, bc_tract_geom + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# 7.2. Hellinger ------------------------------------------------------

### All
theme_set(theme_bw())
ordu3 = ordinate(ps_trans_all, "PCoA", "euclidean")
plot_ordination(ps_trans_all, ordu3)
eucli_all = plot_ordination(ps_trans_all, ordu3, color="Tract", shape="Tract", axes=1:2)
eucli_all_geom = eucli_all + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
eucli_all_geom

### Tract
theme_set(theme_bw())
ordu4 = ordinate(ps_trans_tract, "PCoA", "euclidean")
plot_ordination(ps_trans_tract, ordu4)
eucli_tract = plot_ordination(ps_trans_tract, ordu4, color="Tract", shape="Tract", axes=1:2)
eucli_tract_geom = eucli_tract + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
eucli_tract_geom

### Two graph in one
ggarrange(eucli_all_geom, eucli_tract_geom + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

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


### Tract
theme_set(theme_bw())
ordu6 = ordinate(ps_UNIFRAC_tract, "PCoA", "wunifrac")
plot_ordination(ps_UNIFRAC_tract, ordu6)
wunifrac_tract = plot_ordination(ps_UNIFRAC_tract, ordu6, color="Tract", shape="Tract", axes=1:2)
wunifrac_tract_geom = wunifrac_tract + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
wunifrac_tract_geom

### Two graph in one
ggarrange(wunifrac_all_geom, wunifrac_tract_geom + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Unweighted ---- 
### All
theme_set(theme_bw())
ordu7 = ordinate(ps_UNIFRAC_all, "PCoA", "unifrac")
plot_ordination(ps_UNIFRAC_all, ordu7)
unifrac_all = plot_ordination(ps_UNIFRAC_all, ordu7, color="Tract", shape="Tract", axes=1:2)
unifrac_all_geom = unifrac_all + geom_point(size=5) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
unifrac_all_geom

### Tract
theme_set(theme_bw())
ordu8 = ordinate(ps_UNIFRAC_tract, "PCoA", "unifrac")
plot_ordination(ps_UNIFRAC_tract, ordu8)
unifrac_tract = plot_ordination(ps_UNIFRAC_tract, ordu8, color="Tract", shape="Tract", axes=1:2)
unifrac_tract_geom = unifrac_tract + geom_point(size=3) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  stat_ellipse(aes(group = Tract), linetype = 2)
unifrac_tract_geom

### Two graph in one
ggarrange(unifrac_all_geom, unifrac_tract_geom + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

### All distances - All
ggarrange(bc_all_geom, eucli_all_geom, wunifrac_all_geom, unifrac_all_geom + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "right")

### All distances - Tract
ggarrange(bc_tract_geom, eucli_tract_geom, wunifrac_tract_geom, unifrac_tract_geom + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "right")

### Distance of interest
#All
final_plot <- ggarrange(bc_all_geom, wunifrac_all_geom + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = TRUE, legend = "right")
final_plot
ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/PCoA_Fungi.tiff", plot = final_plot, width = 10, height = 6, units = "in", dpi = 300)


#Tract
ggarrange(bc_tract_geom, wunifrac_tract_geom + rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1,
                         common.legend = TRUE, legend = "right")
final_table


# 8. Permanova ----------------------------------

# 8.1. Bray-Curtis ----------------------------------------------------

perma_bc_all <-adonis2(bc_dist_all~Bloc+Tract+Tree, data=sample_info, permutations=9999)
perma_bc_all
perma_bc_all2 <- mutate(perma_bc_all, Distance = "BC_all")

perma_bc_tract <-adonis2(bc_dist_tract~Bloc+Tract+Tree, data=sample_info_tract, permutations=9999)
perma_bc_tract
perma_bc_tract2 <- mutate(perma_bc_tract, Distance = "BC_tract")

perma_bc_4t <-adonis2(bc_dist_4t~Bloc+Tract+Tree, data=sample_info_4t, permutations=9999)
perma_bc_4t
perma_bc_4t2 <- mutate(perma_bc_4t, Distance = "BC_tract")

perma_bc_3t <-adonis2(bc_dist_3t~Bloc+Tract+Tree, data=sample_info_3t, permutations=9999)
perma_bc_3t
perma_bc_3t2 <- mutate(perma_bc_3t, Distance = "BC_tract")

perma_bc_2t <-adonis2(bc_dist_2t~Bloc+Tract+Tree, data=sample_info_2t, permutations=9999)
perma_bc_2t
perma_bc_2t2 <- mutate(perma_bc_2t, Distance = "BC_tract")

# 8.2. Hellinger ------------------------------------------------------

perma_eucli_all <-adonis2(eucli_dist_all~Tract+Tree, data=sample_info, permutations=9999)
perma_eucli_all
perma_eucli_all2 <- mutate(perma_eucli_all, Distance = "Eucli_all")

perma_eucli_tract <-adonis2(eucli_dist_tract~Tract+Tree, data=sample_info_tract, permutations=9999)
perma_eucli_tract
perma_eucli_tract2 <- mutate(perma_eucli_tract, Distance = "Eucli_tract")

# 8.3. Weighted UNIFRAC -----------------------------------------------

perma_wuni_all <-adonis2(wunifrac_dis_all~Bloc+Tract+Tree, data=sample_info_UNIFRAC, permutations=9999)
perma_wuni_all
perma_wuni_all2 <- mutate(perma_wuni_all, Distance = "Wunifrac_all")

perma_wuni_tract <-adonis2(wunifrac_dis_tract~Bloc+Tract+Tree, data=sample_info_tract_UNIFRAC, permutations=9999)
perma_wuni_tract
perma_wuni_tract2 <- mutate(perma_wuni_tract, Distance = "Wunifrac_tract")

perma_wuni_4t <-adonis2(wunifrac_dis_4t~Bloc+Tract+Tree, data=sample_info_4t_UNIFRAC, permutations=9999)
perma_wuni_4t
perma_wuni_4t2 <- mutate(perma_wuni_4t, Distance = "Wunifrac_tract")

perma_wuni_3t <-adonis2(wunifrac_dis_3t~Bloc+Tract+Tree, data=sample_info_3t_UNIFRAC, permutations=9999)
perma_wuni_3t
perma_wuni_3t2 <- mutate(perma_wuni_3t, Distance = "Wunifrac_tract")

perma_wuni_2t <-adonis2(wunifrac_dis_2t~Bloc+Tract+Tree, data=sample_info_2t_UNIFRAC, permutations=9999)
perma_wuni_2t
perma_wuni_2t2 <- mutate(perma_wuni_2t, Distance = "Wunifrac_tract")

# 8.4. Unweighted UNIFRAC  --------------------------------------------

perma_uni_all <-adonis2(unifrac_dis_all~Tract+Tree, data=sample_info_UNIFRAC, permutations=9999)
perma_uni_all
perma_uni_all2 <- mutate(perma_uni_all, Distance = "Unifrac_all")

perma_uni_tract <-adonis2(unifrac_dis_tract~Tract+Tree, data=sample_info_tract_UNIFRAC, permutations=9999)
perma_uni_tract
perma_uni_tract2 <- mutate(perma_uni_tract, Distance = "Unifrac_tract")

# 8.5. Creating a table with all the permanova  -----------------------

### Create the table
perma_all <- bind_rows(perma_bc_all2, perma_bc_tract2, perma_eucli_all2, 
                       perma_eucli_tract2, perma_wuni_all2, perma_wuni_tract2, 
                       perma_uni_all2, perma_uni_tract2)

write.table(perma_all, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_permanova_tot.tsv", sep="\t", quote=F, col.names=NA)


### Only keep the p.value and the distance name and remove the 
#Total and residual line
perma_pr <- select(perma_all, "Distance", p.value = "Pr(>F)")
perma_pr2 <- subset(perma_pr, !grepl(c("Total|Residual"), rownames(perma_pr)))

### Remove the "..." and the number from the row names
rownames(perma_pr2) <- gsub("\\.\\.\\.\\d+", "", rownames(perma_pr2))

### Create a new column to store the Tract or Tree information
perma_pr2$Group <- ifelse(grepl("Tract", rownames(perma_pr2)), "Tract", "Tree")

### Pivot the table to reshape it
reshaped_perma <- pivot_wider(perma_pr2, names_from = Distance, values_from = p.value)

### Set the row names to "Tract" and "Tree"
perma_name <- reshaped_perma$Group
# or: perma_name <- c("Tract", "Tree")

### Delete the Group column 
reshaped_perma$Group <- NULL

### Set the row names to "Tract" and "Tree"
rownames(reshaped_perma) <- perma_name

write.table(reshaped_perma, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_permanova.tsv", sep="\t", quote=F, col.names=NA)
