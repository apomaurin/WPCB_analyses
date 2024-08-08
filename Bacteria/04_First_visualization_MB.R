###
#Title : "First_Visualization_Methode_Bacteria"
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
install.packages("dendextend")

# 1.2. Load the required packages ----------------------

library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(dendextend) ;packageVersion("dendextend")

# 1.3. Define the working directory --------------------

### First define where you're gonna work : 
setwd("C:/Users/apoll/OneDrive/Documents/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Bacteria")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 2. Load the data --------------------------------------------------------

count_tab <- as.matrix(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Count_Tab.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.matrix(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Tax_Tab.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))

#creation de sample info: 

samples.out <- rownames(t(count_tab))
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
tract <- c("2_tracts", "3_tracts", "4_tracts","2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts", "2_tracts", "2_tracts","2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts","2_tracts", "2_tracts","3_tracts","4_tracts","2_tracts","3_tracts","4_tracts" ,"Control","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery")
tree <- c("Tree_1","Tree_1", "Tree_1", "Tree_2","Tree_2","Tree_2", "Tree_3","Tree_3","Tree_3","Tree_4","Tree_4","Tree_4","Tree_5","Tree_6","Tree_7","Tree_7","Tree_7","Tree_8","Tree_8","Tree_8","Tree_9","Tree_11","Tree_11","Tree_11","Tree_12", "Tree_12", "Tree_12", "Temoin", "Tree_1", "Tree_2","Tree_3","Tree_4","Tree_5","Tree_6","Tree_7","Tree_8","Tree_9","Tree_11","Tree_12")
status <- c("Tract", "Tract", "Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract", "Tract", "Tract","Tract","Tract","Tract","Tract","Tract","Tract","Tract", "Tract","Tract","Tract","Tract","Tract","Tract" ,"Control","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery","Gallery")
bloc <- c("B2", "B2", "B2","B2","B2","B2","B2","B2","B2","B2","B2","B2","B2","B2","B6","B6","B6","B6","B6","B6","B6","B6","B6","B6","B6","B6","B6","Control","B2","B2","B2","B2","B2","B2","B6","B6","B6","B6","B6")
sample_info <- data.frame(Subject=subject, Tract=tract, Tree=tree, Status=status, Bloc=bloc)
head(sample_info)
write.table(sample_info, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Sample_info.tsv", sep="\t", quote=F, col.names=NA)
rownames(sample_info) <- samples.out

# 3. Create a phyloseq object ---------------------------------------------

samples.out <- rownames(t(count_tab))
rownames(sample_info) <- samples.out

phylo_data <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))


# 3. First visualization of the data  -------------------------------------

# 3.1. Rarefaction curve ---------------------------------------

S <- specnumber(t(count_tab)) # observed number of species
raremax <- min(rowSums(t(count_tab)))
Srare <- rarefy(t(count_tab), raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(t(count_tab), step = 20, sample = raremax, col = "orange", cex = 0.8)

# 3.2. Species richness --------------

count.pre <- rowSums(count_tab)
count.pre
barplot(count.pre, main = "Species richness", ylim=c(0,80000),
        xlab = "ASV", ylab = "Number of species", col = "grey ", las = 1, cex.axis= 0.7, cex.names = 0.5)



# 4. Visualization of relative abundance ----------------------------------

### Here we'll agglomerate the reads to the phylum-level using Phyloseq plot 
    #and the relative abundance by Status


# 4.1. Count and transformation --------------------------------

# Get count phyla ----
table(phyloseq::tax_table(phylo_data)[, "Phylum"])

# Transformation of the data ----

### Conversion to relative abundance
### When doing the transformation on my Phyloseq object, it seems a problem occur. 
    #Therefor, I'll do the transformation on my count_tab and then put in in my 
    #Phyloseq object. 

# Transformation: 
### For the transformation, we'll use the Hellinger method. It's the one that 
    #is widely used in microbial ecology as it's particularly suited to species 
    #abundance date. In fact, this transformation gives low weights to variable 
    #with low counts and many zero. (2001, Legendre)

data_stand <- decostand(count_tab, method = "hellinger")
#View(data_stand)

# Phyloseq Object:
samples.out <- rownames(t(data_stand))
rownames(sample_info) <- samples.out

phylo_data_trans <- phyloseq(otu_table(data_stand, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))

# Checking everything is fine:
phyloseq::otu_table(phylo_data_trans)[1:5, 1:5]


# 4.2. Plot ----------------------------------------------------

# BarPlot ----

### If you want to look at a particularly type of sample. 
    #In the example I wanted to plot only the sample belonging to the tree 4
#tree4 <- subset_samples(phylo_data_trans, Tree == "Tree_4")


phyloseq::plot_bar(phylo_data_trans, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Tract , scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

### It seems that there is more stability in the community present in the 
    #gallery than in the tract. 
### Moreover, it seems clear that the control have been contaminated by a 
    #gallery. Somme of the sample might have been too. 

phyloseq::plot_bar(phylo_data_trans, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

### From the tree perspective, I think that some of the insect sample have 
    #been contaminated. 
### Would be interested to check how many tree were used in the composite for
    #tree 5.


# BoxPlot ----

### Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_data, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

### Melt and plot
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Tract, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

### this allows us to better see which community is shared between the Gallery 
    #and the insect. 

### For example it seems that the Proteobacteria is more prominent in the 
    #insect.Similarly, almost all the other community are in majority in the 
    #gallery or equal.


# Conclusion ----

### To conclude, some of the sample are clearly contaminated. I evaluated that 
#about 4 to 5 sample are contaminated. 

### I will remove them from the sample, but first, I'd like to do a hierarchical
#clustering to validate see if my suspicions are correct. 



# 5. Hierarchical clustering ----------------------------------------------

# 5.1. Distance calculation ------------------------------------

### To examine how samples cluster, we use some measure of taxonomic 
    #(dis)similarity. For example, Bray-Curtis dissimilarity is a popular one 
    #in microbial ecology.

### I'm not sure of the distance I'll use, so I'll repeat the hierarchical 
    #clustering with different calculation distance. I'll try and decide. 

### Once again, the calculation is not done is the phyloseq object. 

# Bray-Curtis ----
bc_rel_asv <- data.frame(phyloseq::otu_table(phylo_data_trans))
bc_rel_asv <- t(bc_rel_asv)
bc_dist <- vegan::vegdist(bc_rel_asv, method = "bray") 
as.matrix(bc_dist)[1:5, 1:5]  

# Hellinger ----
eucli_rel_asv <- data.frame(phyloseq::otu_table(phylo_data_trans))
eucli_rel_asv <- t(eucli_rel_asv)
eucli_dist <- vegan::vegdist(eucli_rel_asv, method = "euclidean") 
as.matrix(eucli_dist)[1:5, 1:5]

# Robust Aitchison ----
ra_rel_asv <- data.frame(phyloseq::otu_table(phylo_data_trans))
ra_rel_asv <- t(ra_rel_asv)
ra_dist <- vegan::vegdist(ra_rel_asv, method = "robust.aitchison") 
as.matrix(ra_dist)[1:5, 1:5]


# 5.2. Plot ----------------------------------------------------

# Bray-Curtis ----

### Create the hierarchical clustering:
    #--> as.dendrogram convert the hierearchical clustering object into a dendrogram 
   #object that can be manipulated (put color, etc.)
bc_clust <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
### Put some color and adjust the margin:
meta <- data.frame(phyloseq::sample_data(phylo_data_trans))
colorCode <- c(Control = "tomato", Tract = "olivedrab", Gallery = "sandybrown")
labels_colors(bc_clust) <- colorCode[meta$Status][order.dendrogram(bc_clust)]
par(mar = c(10, 5, 4, 2) + 0.1) #adjust the margin so the name fit in the image
### Plot:
plot(bc_clust)

heatmap(as.matrix(bc_dist), Rowv = bc_clust, symm = TRUE, margin = c(3,3))

# Hellinger ----

### Create the hierarchical clustering:
eucli_clust <- as.dendrogram(hclust(eucli_dist, method = "ward.D2"))
### Put some color and adjust the margin:
labels_colors(eucli_clust) <- colorCode[meta$Status][order.dendrogram(eucli_clust)]
par(mar = c(10, 5, 4, 2) + 0.1) #adjust the margin so the name fit in the image
### Plot:
plot(eucli_clust)

# Robust Aitchison ----

### Create the hierarchical clustering:
ra_clust <- as.dendrogram(hclust(ra_dist, method = "ward.D2"))
### Put some color and adjust the margin:
labels_colors(ra_clust) <- colorCode[meta$Status][order.dendrogram(ra_clust)]
par(mar = c(10, 5, 4, 2) + 0.1) #adjust the margin so the name fit in the image
### Plot:
plot(ra_clust)

### In a sure way, the one contaminated are: 
    # --> A2 - 3t 
    # --> A7 - 2t 
### In a less sure way:
    # --> A3 - 4t (only BC)
    # --> A6 - 2t (only BC)
    # --> A4 - 2t (only RA)

### I'm not sure weather I should remove them all or just the 2 sure one, or
   #only the one from BC. I'll have a tendency for BC or all.  



# 6. Remove the contaminated samples --------------------------------------

# 6.1. Found the line ------------------------------------------

### First, make sure your table are a data.frame. If not, it won't work anyway.  
count_tab_filter <- as.data.frame(count_tab)

### Then, remove the contaminated samples and save them in a new file: 

# Column name:
count_tab_filter <- count_tab_filter[,!names(count_tab_filter) %in% c("A2_3t_bacteries", "A7_2t_bacteries", "A3_4t_bacteries", "A6_2t_bacteries", "Temoin_1_bacteries")]
write.table(count_tab_filter, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Count_Tab_No_Conta.tsv", sep="\t", quote=F, col.names=NA)

# Row name:
sample_info_filter <- sample_info[!(sample_info$Subject) %in% c("A2_3t_bacteries", "A7_2t_bacteries", "A3_4t_bacteries", "A6_2t_bacteries", "Temoin_1_bacteries"),]
write.table(sample_info_filter, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Seq2_MB_Sample_Info_No_Conta.tsv", sep="\t", quote=F, col.names=NA)



# 6.2. Save the  ---------------------------------------

samples.out <- rownames(t(count_tab))
rownames(sample_info) <- samples.out

phylo_data <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))

##Today : 
# Remove the sample contaminated : okay
# pursue the pipeline : 
  # alpha diversity
  # beta diversity 
  # ASV rarefaction curve

 










