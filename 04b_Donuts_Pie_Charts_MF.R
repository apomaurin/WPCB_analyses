###
#Title : "Donuts_Pie_Charts_Fungi"
#Author : Apolline Maurin
# Date : 15/05/2024
###

# 1.1. Install the required packages -------------------
  
  ###Install packages if you don't have it already. otherwise, don't load it. 
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dplyr")
install.packages("ggplot2")
install.packages("forcats")
install.packages("RColorBrewer")
install.packages("tidyverse")
install.packages("vegan")
install.packages("ggtext")
install.packages("eulerr")
install.packages("microbiome")
devtools::install_github('microsud/microbiomeutilities')
BiocManager::install("phyloseq")

# 1.2. Load the required packages ----------------------

library(dplyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(vegan)
library(ggtext)
library(eulerr)
library(microbiome)
library(microbiomeutilities)
library(phyloseq)

# 1.3. Define the working directory --------------------

### First define where you're gonna work : 
setwd("D:/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Fungi")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 2. Load the data --------------------------------------------------------

count_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Count_Tab_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.matrix(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Tax_Tab.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Sample_info_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))

# 3. Donuts and pie charts ------------------------------------------------

## Create phyloseq object ----
samples.out <- rownames(t(count_tab))
rownames(sample_info) <- samples.out

ps <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
               sample_data(sample_info), 
               tax_table(taxa_tab))

## rarefy without replacement ----
ps_rarefied <- rarefy_even_depth(ps, rngseed=1, replace=F)
count_rarefied <- as.data.frame(otu_table(ps_rarefied))

## Transformation of the table ----

### Reverse the table and put the rowname in a column name Subject or ASV
count_tabr <- t(count_rarefied)
count_tabr_subj_tot <- data.frame(Subject = rownames(count_tabr), count_tabr, row.names=NULL)

taxa_tab_asv <- data.frame(ASV = rownames(taxa_tab), taxa_tab, row.names=NULL)

### Load the data, give a name to the values in the count tab and remove the NA

# Process the otu_counts
otu_counts <- count_tabr_subj_tot %>%
  select(Subject, starts_with("ASV")) %>%
  pivot_longer(-Subject, names_to = "ASV", values_to = "count")

# Process the taxonomy
taxonomy <- taxa_tab_asv %>%
  drop_na(Family)

### Join the count, sample info and taxonomy + Compute percentage + define
#the variables considered for each chart (here = Status)
otu_rel_abund <- inner_join(sample_info, otu_counts, by = "Subject") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Subject) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to = "level",
               values_to = "taxon") %>%
  mutate(Status = factor(Status, levels = c("Gallery",
                                            "Tract")))

### Select a taxonomic level + define the percentage
taxon_rel_abund <- otu_rel_abund %>%
  filter(level == "Family") %>%
  group_by(Status, Subject, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(Status, taxon) %>%
  summarize(mean_rel_abund = 100 * mean(rel_abund), .groups = "drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "*\\1*"))


### Conserve percentage superior at 3%, class under "Other" the rest
taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 3,
            mean = mean(mean_rel_abund),
            .groups = "drop")

## Join and plot ----
donuts_plots <-
  inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Relative abundance < 3%", taxon)) %>%
  group_by(Status, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(x=Status, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 0.75) +
  scale_fill_manual(name=NULL,
                    breaks=c("*Helotiales_fam_Incertae_sedis*", "*Chionosphaeraceae*", 
                             "*Debaryomycetaceae*", "*Byssocorticiaceae*",
                             "*Hyaloscyphaceae*", "*Nectriaceae*",
                             "Relative abundance < 3%"),
                    values = c(brewer.pal(6, "RdBu"), "gray")) +
  scale_x_discrete(breaks=c("Tract",
                            "Gallery"),
                   labels=c("Tract",
                            "Gallery")) +
  coord_polar(theta = "y") +
  scale_y_continuous(breaks = NULL) +
  labs(x=NULL, 
       y=NULL) +
  theme_classic() +
  theme(legend.key.size = unit(10, "pt"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_markdown())
donuts_plots

### Save
ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/1_Microbiome_Composition/Donuts_Pie_Charts_Family_Fungi.jpeg", plot = donuts_plots, width = 12, height = 6, units = "in", dpi = 300)
ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/1_Microbiome_Composition/Donuts_Pie_Charts_Family_Fungi.pdf", plot = donuts_plots, width = 12, height = 6, units = "in", dpi = 300)

#Phylum
#scale_fill_manual(name=NULL,
#                  breaks=c("*Basidiomycota*", "*Ascomycota*",
#                           "Other"),
#                  values = c(brewer.pal(10, "RdBu"), brewer.pal(3, "OrRd"), "gray")) +
#Family
#scale_fill_manual(name=NULL,
#                  breaks=c("*Helotiales_fam_Incertae_sedis*", "*Chionosphaeraceae*", 
#                           "*Debaryomycetaceae*", "*Byssocorticiaceae*",
#                           "*Hyaloscyphaceae*", "*Nectriaceae*",
#                           "Relative abundance < 3%"),
#                  values = c(brewer.pal(6, "RdBu"), "gray")) +

## Statistics ----

# Conserve percentage superior at 3%, class under "Other" the rest
taxon_pool2 <- taxon_rel_abund %>%
  group_by(Status, taxon) %>%
  summarize(pool = max(mean_rel_abund) < 3,
            mean = mean(mean_rel_abund),
            .groups = "drop")

# Calculate stats for each status
stats_by_status <- taxon_pool2 %>%
  group_by(Status) %>%
  summarise(
    total_asvs = n(),
    num_asvs_under_3 = sum(pool),
    percentage_under_3 = (sum(pool) / n()) * 100,
    num_asvs_above_3 = n() - sum(pool),
    percentage_above_3 = ((n() - sum(pool)) / n()) * 100
  )

# Print the results
stats_by_status %>%
  rowwise() %>%
  mutate(
    print_output = paste0(
      "Status: ", Status, "\n",
      "  Total ASVs: ", total_asvs, "\n",
      "  Number of ASVs with mean relative abundance < 3: ", num_asvs_under_3, "\n",
      "  Percentage of ASVs with mean relative abundance < 3: ", round(percentage_under_3, 2), "%\n",
      "  Number of ASVs with mean relative abundance >= 3: ", num_asvs_above_3, "\n",
      "  Percentage of ASVs with mean relative abundance >= 3: ", round(percentage_above_3, 2), "%\n"
    )
  ) %>%
  pull(print_output) %>%
  cat(sep = "\n\n")

# Additional calculations for shared ASVs

shared_asvs <- taxon_rel_abund %>%
  group_by(taxon) %>%
  filter(n_distinct(Status) > 1) %>%
  summarise(
    mean_rel_abund = mean(mean_rel_abund),
    pool = max(mean_rel_abund) < 3,
    .groups = "drop"
  )

shared_stats <- shared_asvs %>%
  summarise(
    total_asvs = n(),
    num_asvs_under_3 = sum(pool),
    percentage_under_3 = (sum(pool) / n()) * 100,
    num_asvs_above_3 = n() - sum(pool),
    percentage_above_3 = ((n() - sum(pool)) / n()) * 100
  ) %>%
  mutate(Status = "Shared")

# Print shared stats
shared_stats %>%
  rowwise() %>%
  mutate(
    print_output = paste0(
      "Status: ", Status, "\n",
      "  Total ASVs: ", total_asvs, "\n",
      "  Number of ASVs with mean relative abundance < 3: ", num_asvs_under_3, "\n",
      "  Percentage of ASVs with mean relative abundance < 3: ", round(percentage_under_3, 2), "%\n",
      "  Number of ASVs with mean relative abundance >= 3: ", num_asvs_above_3, "\n",
      "  Percentage of ASVs with mean relative abundance >= 3: ", round(percentage_above_3, 2), "%\n"
    )
  ) %>%
  pull(print_output) %>%
  cat(sep = "\n\n")



mean_tract <- subset(taxon_pool2, !(Status %in% c('Gallery')))
test1 = sum(mean_tract$mean)

test = 4.440939e+01 + 4.369997e+01 + 5.997797e+00

test/test1*100

# 3. Venn diagram ---------------------------------------------------------

## Hellinger transformation ----
data_stand <- decostand(count_tab, method = "hellinger")

## Create phyloseq object ----
samples.out <- rownames(t(data_stand))
rownames(sample_info) <- samples.out

ps_trans <- phyloseq(otu_table(data_stand, taxa_are_rows=TRUE), 
                     sample_data(sample_info), 
                     tax_table(taxa_tab))

## Aggregate data to family level ----
ps_family <- tax_glom(ps_rarefied, taxrank = "Family")

# Verify the new taxonomic level
taxa_names(ps_family)[1:5]

## Prepare the names to easily check wich ASV and microorganisms are found ----
### first remove "ASV_"
taxa_names(ps_family) <- gsub( "ASV_","", taxa_names(ps_family))
taxa_names(ps_family)[1:5]
### format names
ps_family_name <- format_to_besthit(ps_family)
### check names
taxa_names(ps_family_name)[1:5]

## Identify the levels of each group (here found in Status) ----
status <- unique(as.character(meta(ps_family_name)$Status))
print(status)

## Initialize an empty list to store core taxa information ----
list_core <- list()

## Loop through each status to identify core taxa at the family level ----
for (n in status) {
  ps.sub <- subset_samples(ps_family_name, Status == n)
  
  core_m <- core_members(ps.sub, 
                         detection = 0.001, 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

print(list_core)

## Identify common core microbiome across all groups ----
common_core <- Reduce(intersect, list_core)
print(common_core)

## Combine results into a single data frame ----
core_table <- data.frame(Taxon = unique(unlist(list_core)))

for (n in status) {
  core_table[[n]] <- core_table$Taxon %in% list_core[[n]]
}

core_table$CommonCore <- core_table$Taxon %in% common_core

# Save the core microbiome table to a CSV file
write.table(core_table, "./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Microbiome_Composition/Seq2_MF_core_microbiome_Phylum.tsv", sep="\t", quote=F, col.names=NA)


## Specify colors and plot Venn diagram ----
mycols <- brewer.pal(6, "RdBu")
venn_diagram <- plot(venn(list_core),
                     fills = mycols)
venn_diagram

# Save
ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/1_Microbiome_Composition/Venn_Diagram_Family_Fungi.jpeg", plot = venn_diagram, width = 12, height = 6, units = "in", dpi = 300)
ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/1_Microbiome_Composition/Venn_Diagram_Family_Fungi.pdf", plot = venn_diagram, width = 12, height = 6, units = "in", dpi = 300)

