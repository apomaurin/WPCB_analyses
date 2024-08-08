###
#Title: "Ancombc_Fungi"
#Author: Apolline Maurin
# Date: 14/12/2023
###

# This code uses ANCOMBC in order to determine taxa whose absolute abundances
#are significantly different with changes in the covariate of interest.

# Sources: 
#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#https://www.bioconductor.org/packages/release/bioc/manuals/ANCOMBC/man/ANCOMBC.pdf

# 1. R SET UP -------------------------------------------------------------


# 1.1. Install the required packages -------------------

###Install packages if you don't have it already. otherwise, don't load it. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")


BiocManager::install("phyloseq")
BiocManager::install("ANCOMBC")
BiocManager::install("mia")
BiocManager::install("microbiome")
install.packages("dplyr")
install.packages("DT")

# 1.2. Load the required packages ----------------------

library(phyloseq) ;packageVersion("phyloseq")
library(ANCOMBC);packageVersion("ANCOMBC")
library(mia);packageVersion("mia")
library(microbiome);packageVersion("microbiome")
library(dplyr);packageVersion("dplyr")
library(DT);packageVersion("DT")


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


# 3. ANCOM-BC --------------------------------

# 3.1. Create a phyloseq object with the count_table ------------------
## ANCOM-BC do a normalization
samples.out <- rownames(t(count_tab))
rownames(sample_info) <- samples.out

phylo_data <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))

# 3.2. Agglomerate under the desired taxon ----------------------------

#List of comparisons
comparisons <- c("Gallery", "2_tracts", "2_tracts", "4_tracts")
output_dir <- "./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/ANCOMBC/"

# Loop through comparisons
for (comparison in comparisons) {
  # Set seed for reproducibility
  set.seed(123)
  
  ## 1- Option 1:
  #--> Taxa are agglomerated when the genus is the same 
  #--> limit = Unknown taxon are grouped together
  
  ### Aggregate to the desired taxonomy level
  ps_aggrete = aggregate_taxa(phylo_data, "Family")
  
  
  # Convert the tract variable to a factor and reorder for the current comparison
  sample_data(ps_aggrete)$Tract <- as.factor(sample_data(ps_aggrete)$Tract)
  sample_data(ps_aggrete)$Tract <- relevel(sample_data(ps_aggrete)$Tract, comparison)
  
  # Run ANCOM-BC
  out_aggregate <- ancombc(data = NULL, tax_level = "Family", phyloseq = ps_aggrete,
                           formula = "Tract",
                           p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,
                           group = "Tract", struc_zero = TRUE, neg_lb = TRUE,
                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                           alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
  
  # Save result to a variable with a corresponding name
  assign(paste("res_aggre_", comparison, sep = ""), out_aggregate$res)
  
  # Save result to a file with a corresponding name
  filename <- paste(output_dir, "Ancom_BC_Aggregate_", comparison, "_Genus.tsv", sep = "")
  write.table(out_aggregate$res, file = filename, sep = "\t", quote = FALSE, col.names = NA)
}


### Here is the same thing, but with the option 2. if you want to use this 
#method, replace it in the previous loop.

## 2- Option 2:
#--> Taxa are agglomerate with the genus or the family 
#--> Unknown taxa are not agglomerate together
#--> Unknown taxa are agglomerate according to the family

### Converts phyloseq objects into TreeSummarizedExperiment objects
#tse = mia::makeTreeSummarizedExperimentFromPhyloseq(phylo_data)
### Aggregate to the desired taxonomy level
#tse_agglo = agglomerateByRank(tse, "Genus")
### Converts TreeSummarizedExperiment objects into phyloseq objects
#ps_agglomerate = makePhyloseqFromTreeSummarizedExperiment(tse_agglo)

# Order for comparison 

#sample_data(ps_agglomerate)$Tract <- as.factor(sample_data(ps_agglomerate)$Tract)
#sample_data(ps_agglomerate)$Tract <- relevel(sample_data(ps_agglomerate)$Tract, "4_tracts") 

# ANCOMBC
#set.seed(123) 
#out_agglomerate = ancombc(data = NULL, tax_level = "Genus", phyloseq = ps_agglomerate,
  #                        formula = "Bloc",
  #                        p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0,
  #                        group = "Bloc", struc_zero = TRUE, neg_lb = TRUE,
  #                        tol = 1e-5, max_iter = 100, conserve = TRUE,
  #                        alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)

#res_agglo = out_agglomerate$res
#write.table(res_agglo, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/ANCOMBC/Ancom_BC_Agglomerate_Bloc_MF.tsv", sep="\t", quote=F, col.names=NA)

# - Formula: how the abundance of a taxon depend of a variable
##--> Example : How the abundance depend on the number of tract
# - prv_cut: kind of a threshold, Taxa with prevalence (the proportion of
##samples in which the taxon is present) less than prv_cut will be 
##excluded in the analysis
# - Lib_cut: cut the library according to a threshold
##--> If you don't want to cut anything, put 0
# - Group: the group variable in metadata = the one you compare with the other
##--> If you do define it beforehand, it will take the first character 
##in the group cited (here Status --> 2_tracts)
# - Struc_zero: Whether to detect structural zeros based on group (put TRUE)
# - Neg_lb: whether to classify a taxon as a structural zero using its asymptotic
#lower bound (put TRUE)
# - Conserve: whether to use a conservative variance estimator for the test 
##statistic (put TRUE)
# - Global: whether to perform the global test (put TRUE)
# - Verbose: if you want the level of output or messages displayed during the 
## execution of ancombc

#For the parameters that are not described, the default parameters is used.


# 3.5. Look at the results --------------------------------------------
## Change the name to look at the desired option

# A) Log fold change (describe the ratio between two values)
## --> lfc > 0 = more of this taxon in the tract
## --> lfc < 0 = more of this taxon in the gallery
## --> lfc = 0 = same amount in both

res_aggre_4_tracts$lfc
tab_lfc4 = res_aggre_4_tracts$lfc
col_name = c("Taxon", "Intercept", "2 tracts vs 4 tracts", "3 tracts vs 4 tracts", "Galleries vs 4 tracts") #Change the name with yours
colnames(tab_lfc4) = col_name
tab_lfc4 %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

lfc_all <- full_join(tab_lfc2, tab_lfc4, by = "Taxon")
lfc_ultime <- full_join(lfc_all, tab_lfcG, by = "Taxon")
write.table(lfc_ultime, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/ANCOMBC/lfc_tab_Fungi.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !


# B) Standard error (standard deviation of the distribution)
## = measures the amount of discrepancy that can be expected in a sample 
#estimate compared to the true value in the population.
## --> Closer to 0, the better

tab_se4 = res_aggre_4_tracts$se
colnames(tab_se4) = col_name
tab_se4 %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
## --> SE seems pretty hight for my results

se_all <- full_join(tab_se2, tab_se4, by = "Taxon")
se_ultime <- full_join(se_all, tab_seG, by = "Taxon")
write.table(se_ultime, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/ANCOMBC/se_tab_Fungi.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !

# C) Test Statistics (lfc/SE)
## --> Essentially, it is a count of the number of sub-hypotheses that have 
#passed for a given species
## --> So if you have 1000 species, and W=60, for ASV k, then H0k is rejected 60 
#times.This basically means, that the ratio ASV k and 60 other species were 
#detected to be significantly different across the x and y groups
# https://forum.qiime2.org/t/specify-w-cutoff-for-anacom/1844/10 

tab_w4 = res_aggre_4_tracts$W
colnames(tab_w4) = col_name
tab_w4 %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

w_all <- full_join(tab_w2, tab_w4, by = "Taxon")
w_ultime <- full_join(w_all, tab_wG, by = "Taxon")
write.table(w_ultime, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/ANCOMBC/w_tab_Fungi.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !
# D) P-value
# --> p-value < 0.05 generally accepted

tab_p4 = res_aggre_4_tracts$p_val
colnames(tab_p4) = col_name
tab_p4 %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

p_all <- full_join(tab_p2, tab_p4, by = "Taxon")
p_ultime <- full_join(p_all, tab_pG, by = "Taxon")
write.table(p_ultime, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/ANCOMBC/pvalue_tab_Fungi.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !

# E) Adjusted p-value
## As the matrix has more ASVs than samples, it's possible to have a significant
#p-value even if it's not. 
## The adjusted p-value is a kind of penalties to be more rigorous. 

tab_q4 = res_aggre_4_tracts$q
colnames(tab_q4) = col_name
tab_q4 %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

q_all <- full_join(tab_q2, tab_q4, by = "Taxon")
q_ultime <- full_join(q_all, tab_qG, by = "Taxon")
write.table(q_ultime, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/ANCOMBC/adjusted_pvalue_tab_Fungi.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !

# F) Differentially abundant taxa
## True = difference 
## False = no difference

tab_diff4 = res_aggre_4_tracts$diff_abn
colnames(tab_diff4) = col_name
tab_diff4 %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

diff_all <- full_join(tab_diff2, tab_diff4, by = "Taxon")
diff_ultime <- full_join(diff_all, tab_diffG, by = "Taxon")
write.table(diff_ultime, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_files/ANCOMBC/diff_abundance_tab_Fungi.tsv", sep="\t", quote=F, 
            col.names=NA) #Change the name !

# D) What to with that ? 
## Globally, the data we'll mostly use are the log fold change and the adjusted 
#p-value. 
## It's must not be a surprise if we don't have that much ASVs that are significant

# 3.6. Log transformation of abundance --------------------------------
## This step is optional and allows us have a look at the "absolute" abundance
#by transforming the abundance data

# Save the abundance table
samp_frac = out_aggregate$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out_aggregate$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")


# 3.7 Plot ------------------------------------------------------------

## First, create a dataframe with only the lfc that differ:

# Directory where files will be saved
output_dir <- "./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/ANCOMBC/"

# List of comparisons
comparisons <- c("Gallery", "2_tracts", "4_tracts")

# Loop through comparisons
for (comparison in comparisons) {
  # Load ANCOM-BC result
  res_aggre <- get(paste("res_aggre_", comparison, sep = ""))
  
  # Calculate log fold change
  df_lfc <- data.frame(res_aggre$lfc[, -1] * res_aggre$diff_abn[, -1], check.names = FALSE) %>%
    mutate(taxon_id = res_aggre$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
  
  # Save log fold change to file
  filename_lfc <- paste(output_dir, "LogFoldChange_", comparison, "_Genus.tsv", sep = "")
  write.table(df_lfc, file = filename_lfc, sep = "\t", quote = FALSE, col.names = NA)
  
  # Save log fold change to variable with corresponding name
  assign(paste("df_lfc_", comparison, sep = ""), df_lfc)
  
  # Calculate standard error
  df_se <- data.frame(res_aggre$se[, -1] * res_aggre$diff_abn[, -1], check.names = FALSE) %>% 
    mutate(taxon_id = res_aggre$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
  colnames(df_se)[-1] <- paste0(colnames(df_se)[-1], "SE")
  
  # Save standard error to file
  filename_se <- paste(output_dir, "StandardError_", comparison, "_Genus.tsv", sep = "")
  write.table(df_se, file = filename_se, sep = "\t", quote = FALSE, col.names = NA)
  
  # Save standard error to variable with corresponding name
  assign(paste("df_se_", comparison, sep = ""), df_se)
  
  # Calculate q-values
  df_q <- data.frame(res_aggre$q_val[, -1] * res_aggre$diff_abn[, -1], check.names = FALSE) %>% 
    mutate(taxon_id = res_aggre$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
  colnames(df_q)[-1] <- paste0(colnames(df_q)[-1], "PV")
  
  # Save q-values to file
  filename_q <- paste(output_dir, "QValues_", comparison, "_Genus.tsv", sep = "")
  write.table(df_q, file = filename_q, sep = "\t", quote = FALSE, col.names = NA)
  
  # Save q-values to variable with corresponding name
  assign(paste("df_q_", comparison, sep = ""), df_q)
}

## Then, merge the lfc and se dataframe, eliminate 0,  Sorts the data frame in 
#descending order, add a new column with + or - lfc.

# Define the comparisons
comparisons <- c("2_tracts", "3_tracts", "4_tracts", "Gallery")

# Create an empty list to store the results
result_list <- list()

# Loop over the comparisons
for (comp in comparisons) {
  # Left join the log fold change and standard error data frames
  result <- left_join(get(paste0("df_lfc_", comp)), 
                      get(paste0("df_se_", comp)), 
                      by = "taxon_id")
  
  # Transmute and filter the data
  result <- result %>%
    transmute(taxon_id, 
              Tract = get(paste0("Tract", comp)), 
              SE = get(paste0("Tract", comp, "SE"))) %>%
    filter(Tract != 0) %>%
    arrange(desc(Tract)) %>%
    mutate(direct = ifelse(Tract > 0, "Positive LFC", "Negative LFC"))
  
  # Set factor levels
  result$taxon_id <- factor(result$taxon_id, levels = result$taxon_id)
  result$direct <- factor(result$direct, levels = c("Positive LFC", "Negative LFC"))
  
  # Save the result to the list
  result_list[[paste0("df_fig_tract", comp, "vsG")]] <- result
}

# Print the list structure
str(result_list)



### Tract 2 ----
df_fig_tract2vs3 = df_lfc_3_tracts %>% 
  dplyr::left_join(df_se_3_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract2_tracts, Tract2_tractsSE) %>%
  dplyr::filter(Tract2_tracts != 0) %>% 
  dplyr::arrange(desc(Tract2_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract2_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract2vs3$taxon_id = factor(df_fig_tract2vs3$taxon_id, levels = df_fig_tract2vs3$taxon_id)
df_fig_tract2vs3$direct = factor(df_fig_tract2vs3$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))
df_fig_tract2vs4 = df_lfc_4_tracts %>% 
  dplyr::left_join(df_se_4_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract2_tracts, Tract2_tractsSE) %>%
  dplyr::filter(Tract2_tracts != 0) %>% 
  dplyr::arrange(desc(Tract2_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract2_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract2vs4$taxon_id = factor(df_fig_tract2vs4$taxon_id, levels = df_fig_tract2vs4$taxon_id)
df_fig_tract2vs4$direct = factor(df_fig_tract2vs4$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

df_fig_tract2vsG = df_lfc_Gallery %>% 
  dplyr::left_join(df_se_Gallery, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract2_tracts, Tract2_tractsSE) %>%
  dplyr::filter(Tract2_tracts != 0) %>% 
  dplyr::arrange(desc(Tract2_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract2_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract2vsG$taxon_id = factor(df_fig_tract2vsG$taxon_id, levels = df_fig_tract2vsG$taxon_id)
df_fig_tract2vsG$direct = factor(df_fig_tract2vsG$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

### Tract 3 ----
df_fig_tract3vs2 = df_lfc_2_tracts %>% 
  dplyr::left_join(df_se_2_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract3_tracts, Tract3_tractsSE) %>%
  dplyr::filter(Tract3_tracts != 0) %>% 
  dplyr::arrange(desc(Tract3_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract3_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract3vs2$taxon_id = factor(df_fig_tract3vs2$taxon_id, levels = df_fig_tract3vs2$taxon_id)
df_fig_tract3vs2$direct = factor(df_fig_tract3vs2$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

df_fig_tract3vs4 = df_lfc_4_tracts %>% 
  dplyr::left_join(df_se_4_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract3_tracts, Tract3_tractsSE) %>%
  dplyr::filter(Tract3_tracts != 0) %>% 
  dplyr::arrange(desc(Tract3_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract3_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract3vs4$taxon_id = factor(df_fig_tract3vs4$taxon_id, levels = df_fig_tract3vs4$taxon_id)
df_fig_tract3vs4$direct = factor(df_fig_tract3vs4$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

df_fig_tract3vsG = df_lfc_Gallery %>% 
  dplyr::left_join(df_se_Gallery, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract3_tracts, Tract3_tractsSE) %>%
  dplyr::filter(Tract3_tracts != 0) %>% 
  dplyr::arrange(desc(Tract3_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract3_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract3vsG$taxon_id = factor(df_fig_tract3vsG$taxon_id, levels = df_fig_tract3vsG$taxon_id)
df_fig_tract3vsG$direct = factor(df_fig_tract3vsG$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

### Tract 4 ----
df_fig_tract4vs2 = df_lfc_2_tracts %>% 
  dplyr::left_join(df_se_2_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract4_tracts, Tract4_tractsSE) %>%
  dplyr::filter(Tract4_tracts != 0) %>% 
  dplyr::arrange(desc(Tract4_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract4_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract4vs2$taxon_id = factor(df_fig_tract4vs2$taxon_id, levels = df_fig_tract4vs2$taxon_id)
df_fig_tract4vs2$direct = factor(df_fig_tract4vs2$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))
df_fig_tract4vs3 = df_lfc_3_tracts %>% 
  dplyr::left_join(df_se_3_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract4_tracts, Tract4_tractsSE) %>%
  dplyr::filter(Tract4_tracts != 0) %>% 
  dplyr::arrange(desc(Tract4_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract4_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract4vs3$taxon_id = factor(df_fig_tract4vs3$taxon_id, levels = df_fig_tract4vs3$taxon_id)
df_fig_tract4vs3$direct = factor(df_fig_tract4vs3$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

df_fig_tract4vsG = df_lfc_Gallery %>% 
  dplyr::left_join(df_se_Gallery, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Tract4_tracts, Tract4_tractsSE) %>%
  dplyr::filter(Tract4_tracts != 0) %>% 
  dplyr::arrange(desc(Tract4_tracts)) %>% 
  dplyr::mutate(direct = ifelse(Tract4_tracts > 0, "Positive LFC", "Negative LFC"))
df_fig_tract4vsG$taxon_id = factor(df_fig_tract4vsG$taxon_id, levels = df_fig_tract4vsG$taxon_id)
df_fig_tract4vsG$direct = factor(df_fig_tract4vsG$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

### Galleries ----
df_fig_tractGvs2 = df_lfc_2_tracts %>% 
  dplyr::left_join(df_se_2_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, TractGallery, TractGallerySE) %>%
  dplyr::filter(TractGallery != 0) %>% 
  dplyr::arrange(desc(TractGallery)) %>% 
  dplyr::mutate(direct = ifelse(TractGallery > 0, "Positive LFC", "Negative LFC"))
df_fig_tractGvs2$taxon_id = factor(df_fig_tractGvs2$taxon_id, levels = df_fig_tractGvs2$taxon_id)
df_fig_tractGvs2$direct = factor(df_fig_tractGvs2$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))
df_fig_tractGvs3 = df_lfc_3_tracts %>% 
  dplyr::left_join(df_se_3_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, TractGallery, TractGallerySE) %>%
  dplyr::filter(TractGallery != 0) %>% 
  dplyr::arrange(desc(TractGallery)) %>% 
  dplyr::mutate(direct = ifelse(TractGallery > 0, "Positive LFC", "Negative LFC"))
df_fig_tractGvs3$taxon_id = factor(df_fig_tractGvs3$taxon_id, levels = df_fig_tractGvs3$taxon_id)
df_fig_tractGvs3$direct = factor(df_fig_tractGvs3$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))
df_fig_tractGvs4 = df_lfc_4_tracts %>% 
  dplyr::left_join(df_se_4_tracts, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, TractGallery, TractGallerySE) %>%
  dplyr::filter(TractGallery != 0) %>% 
  dplyr::arrange(desc(TractGallery)) %>% 
  dplyr::mutate(direct = ifelse(TractGallery > 0, "Positive LFC", "Negative LFC"))
df_fig_tractGvs4$taxon_id = factor(df_fig_tractGvs4$taxon_id, levels = df_fig_tractGvs4$taxon_id)
df_fig_tractGvs4$direct = factor(df_fig_tractGvs4$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))


prout <- full_join(df_fig_tract3vs2, df_fig_tract3vs4, by = "taxon_id")
prout2 <- full_join(prout,df_fig_tract4vs2, by = "taxon_id")

proutG <- full_join(df_fig_tractGvs2,df_fig_tract3vsG,df_fig_tractGvs4, by = "taxon_id")
proutG2 <- full_join(proutG, df_fig_tractGvs4, by = "taxon_id")


proutultime <- full_join(prout2, proutG2, by = "taxon_id")

## Finally, plot !
### Bar blot - Waterfall plot ----
p_tract = ggplot(data = df_fig_tract3, 
                 aes(x = taxon_id, y = Tract3_tracts, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Tract3_tracts - Tract3_tractsSE, ymax = Tract3_tracts + Tract3_tractsSE),
                width = 0.2, position = position_dodge(0.05), color = "black") +
  geom_text(aes(label = Tract3_tracts), vjust = -1) + 
  coord_flip() +
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes between insects with four tracts and two tracts") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

p_tract

### Bubble plot ----

# Let's assume you have a data frame 'ancombc_results' with the following columns:
# 'Taxon' for the taxon name, 'Significance' for the significance result from ANCOM-BC,
# 'Abundance' for the mean relative abundance, and 'Location' for the locality.

# Replace 'ancombc_results' with your actual ANCOM-BC results data frame.
# You also need to have 'mean_abundance' and 'sd_abundance' for each taxon, which represents
# the mean and standard deviation of its relative abundance across samples.

# Example of what 'ancombc_results' could look like:
# ancombc_results <- data.frame(
#   Taxon = c('Taxon1', 'Taxon2', ...),
#   Significance = c(TRUE, FALSE, ...), # TRUE if significant, FALSE if not
#   Mean_Abundance = c(10, 5, ...), # Replace with your actual data
#   SD_Abundance = c(2, 1, ...), # Replace with your actual data
#   Location = c('Portugal', 'Madagascar', ...)
# )

# Creating the plot

bubbleplot_all <- ggplot(proutultime, aes(x = taxon_id,)) +
  geom_point(aes(y = "2 tracts vs 3 tracts", size = abs(Tract3_tracts.x), color = direct.x.x), alpha = 0.7, position = position_dodge(width = 0.5)) + # Add transparency for overlap
  geom_point(aes(y = "4 tracts vs 3 tracts", size = abs(Tract3_tracts.y), color = direct.y.x), alpha = 0.7, position = position_dodge(width = 0.01)) +
  geom_point(aes(y = "4 tracts vs 2 tracts", size = abs(Tract4_tracts), color = direct.x.x.x), alpha = 0.7, position = position_dodge(width = 0.5)) + # Add transparency for overlap
  geom_point(aes(y = "Gallery vs 2 tracts", size = abs(TractGallery.x), color = direct.x.y), alpha = 0.7, position = position_dodge(width = 0.01)) +
  geom_point(aes(y = "Gallery vs 3 tracts", size = abs(Tract3_tracts), color = direct.y.y), alpha = 0.7, position = position_dodge(width = 0.01)) +
  geom_point(aes(y = "Gallery vs 4 tracts", size = abs(TractGallery.y), color = direct.y.y.y), alpha = 0.7, position = position_dodge(width = 0.01)) +
  coord_flip() + 
  scale_size_continuous(range = c(0.5, 6), breaks = c(1,3,5)) + # Adjust the size scale to your data
  scale_color_manual(values = c('tan2', 'steelblue'), na.translate =F) + # Color significant taxa differently
  theme_minimal() +
  labs(title = "ANCOM-BC Results", x = "Family", y = NULL, size = "Log fold change", color = "Positive or negative log fold change") +
  theme(legend.position = "right")
bubbleplot_all


ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/3_ANCOMBC/ANCOMBC_BubblePlot_Fungi_all_family.jpeg", plot = bubbleplot_all, width = 12, height = 6, units = "in", dpi = 300)
ggsave("C:/Users/apoll/OneDrive/Documents/Thèse/14_Articles/1_Methode/Figures/3_ANCOMBC/ANCOMBC_BubblePlot_Fungi_all_family.pdf", plot = bubbleplot_all, width = 12, height = 6, units = "in", dpi = 300)



### Heat map ----

## The clr-transformed values of a taxon might give a contradictory answer, it 
#is therefor recommended to use the W statistics as the input values for the 
#heatmap (source: https://github.com/FrederickHuangLin/ANCOMBC/issues/10)


df_fig_bmi = df_lfc %>% 
  filter(Tract3_tracts != 0 | Tract2_tracts != 0 | TractGallery != 0) %>%
  transmute(taxon_id, 
            `4 tracts vs. 2 tracts` = round(Tract2_tracts, 3),
            `4 tracts vs. 3 tracts` = round(Tract3_tracts, 3), 
            `4 tracts vs. Gallery` = round(TractGallery, 3)) %>%
  pivot_longer(cols = `4 tracts vs. 2 tracts`:`4 tracts vs. 3 tracts` : `4 tracts vs. Gallery`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = 0


p_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "steelblue", high = "tan2", mid = "whitesmoke", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_bmi
























