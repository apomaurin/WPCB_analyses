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

BiocManager::install("phyloseq")
install.packages("Biostrings")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("vegan")
install.packages("dendextend")
install.packages("iNEXT")
install.packages('devtools')
install.packages("ggbreak")
install.packages("aplot")
install.packages("stats")
install.packages("FSA")

# 1.2. Load the required packages ----------------------

library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(dendextend) ;packageVersion("dendextend")
library(devtools);packageVersion("devtools")
library(iNEXT);packageVersion("iNext")
library(ggbreak);packageVersion("ggbreak")
library(aplot);packageVersion("aplot")
library(ggbreak);packageVersion("ggbreak")
library(aplot);packageVersion("aplot")
library(stats)
library(FSA)
library(ggsignif)

# 1.3. Define the working directory --------------------

### First define where you're gonna work : 
setwd("C:/Users/apoll/OneDrive/Documents/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Bacteria")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 2. Load the data --------------------------------------------------------

count_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Count_Tab_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.matrix(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Tax_Tab.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Sample_info_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))

# 3. Create a phyloseq object ---------------------------------------------

samples.out <- rownames(t(count_tab))
rownames(sample_info) <- samples.out

ps <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))

# 4. Alpha-diversity ------------------------------------------------------

#"Alpha-diversity is the diversity in a single ecosystem or sample. 
  #The simplest measure is richness, the number of species (or OTUs) observed in 
  #the sample. Other metrics consider the abundances (frequencies) of the OTUs, 
  #for example to give lower weight to lower-abundance OTUs." (Robert Edgar)

# 4.1. Read and observed richness ----------------------

ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(ps),
                         "observed" = phyloseq::estimate_richness(ps, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

faire p.value stqt regression linearire 


### Surprisingly, the more read I have, the less diversity I found... 
    #I have no idea how to interpret that. 
    #I guess it's correlated, no in the expected way, but still.


# 4.2. Standardization ---------------------------------------

# iNext ----

### Before doing any alpha-diversity analyses, I have to normalize/standardize 
    #my data. I cannot use the raw one because of the difference in total reads. 

### I could just rarefied my sample to have the same number of reads everywhere
    #however, this means cutoff some of the reads and I don't think it's a good 
    #idea. One of the reason is that I might cutoff some important count without 
    #even knowing it.  

### Another way is to use the package iNext with the order of Hill. 
    #--> It will define a sample reference from the sample according to the 
         #sample size or the coverage that will help to extrapolation or 
         #rarefied our sample. 
    #--> If the sample is inf to the sample reference, it will implement an 
         #extrapolated estimator of diversity. 

#ech en haut un ech, par colone et les ASV la 1er colonne
alpha_div <- iNEXT(count_tab, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95,nboot=50)
alpha_div

### Legend: 
# Knots = generally 40, it is the size of the sample, between 1 and the 
#endpoint. Each knot represents a particular sample size for which 
#diversity estimates will be calculated (2016, Hsieh)
# Endpoint = sample size which is the the endpoint of the R/E calculation. 
#If NULL = double the reference sample size

# Keep the data frame ----

### Saving the originals 
## Data info
ad_DataInfo <- alpha_div$DataInfo
write.table(ad_DataInfo, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_ad_DataInfo.tsv", sep="\t", quote=F, col.names=NA)

## Estimation
ad_mid <- alpha_div$iNextEst

# Size_based
ad_est_size <- ad_mid$size_based
write.table(ad_est_size, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_ad_est_size.tsv", sep="\t", quote=F, col.names=NA)

# Coverage_based
ad_est_cover <- ad_mid$coverage_based
write.table(ad_est_cover, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_ad_est_cover.tsv", sep="\t", quote=F, col.names=NA)

## Alpha diversity
ad_solo <- alpha_div$AsyEst
write.table(ad_solo, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_AsyEst.tsv", sep="\t", quote=F, col.names=NA)


## Organizing the alpha_diversity data frame

ad_solo_filter <- ad_solo[,-c(4,5,6,7)] #only keep the name and value
ad_only <- pivot_wider(ad_solo_filter, names_from = "Assemblage", 
                       values_from = "Observed") #order by diversity measurement
ad_reverse <- as.data.frame(t(ad_only)[-1,]) #turn the table
colnames(ad_reverse) <- c("Richness","Shannon", "Simpson") #put the name you like
alpha_diversity <- data.frame(Subject = rownames(ad_reverse), ad_reverse) #add a column with the samples name

write.table(alpha_diversity, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_alpha_diversity.tsv", sep="\t", quote=F, col.names=NA)


# 4.3. Alpha diversity visualization -------------------------

# Plotting the alpha diversity ----
## Loading data
### First define where you're gonna work : 
setwd("D:/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Bacteria")
### Then, check if you are where you want to be : 
getwd()
list.files()

alpha_div <- read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_alpha_diversity.tsv", header=T, row.names=1, check.names=F, sep="\t")
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Sample_info_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))

## Data frame with all information ----

alpha_div <- merge(alpha_div,sample_info, by = "Subject", all = TRUE)

# Plotting ----

### Creation of a function to create boxplot with the information data (here 
    #alpha_div) and metric (here Richness, Shannon or Simpson)

create_boxplot <- function(data, metric) {
  data %>%
    gather(key = "metric", value = "value", !!metric) %>%
    mutate(metric = factor(metric, levels = unique(.$metric))) %>%
    ggplot(aes(x = Tract, y = value)) +
    geom_boxplot() +
    geom_jitter(aes(color = Tract), height = 0, width = .2, show.legend = FALSE) +
    labs(x = "", y = "") +
    facet_wrap(~ metric, scales = "free") +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) +
    theme_bw() +
    geom_signif(comparisons = list(c("2_tracts", "3_tracts"), 
                                   c("2_tracts", "4_tracts"),
                                   c("3_tracts", "4_tracts"),
                                   c("3_tracts", "Gallery"),
                                   c("4_tracts", "Gallery")),
                map_signif_level=TRUE,
                size = 0.4,
                y_position = 200)
}

### Legend for futur me: 
#Gather(): used to gather columns together into key-value pairs. 
#mutate(): create a new column to a factor and set the order of the level
#geom_jitter(): create jittered points for each level of Tract

## All
### Plot ! For Shannon and Simpson, a scale break has been added for a 
    #better appreciation of the results. 
plot_richness <- create_boxplot(alpha_div, "Richness")
plot_richness

plot_shannon <- create_boxplot(alpha_div, "Shannon") + scale_y_break(c(5,15), scales = "free")
plot_shannon

plot_simpson <- create_boxplot(alpha_div, "Simpson") + scale_y_break(c(2.5,15), scales = "free")
plot_simpson

### Finally, Plot them all together ! 
plot_alpha_div <- plot_list(gglist=list(plot_richness, plot_shannon, plot_simpson),labels = NULL)
plot_alpha_div
ggsave("D:/Thèse/14_Articles/1_Methode/Figures/Alpha_Diversity_bacteria2.pdf", plot = plot_alpha_div, width = 12, height = 6, units = "in", dpi = 300)


### It seems, once again, that there is a difference between the tracts and the 
   #gallery.

## Only the tracts
###Removing the gallery: 
alpha_div_tract <- subset(alpha_div, !(Status %in% c('Gallery')))

### Then plot ! For Shannon and Simpson, a scale break has been added for a 
#better appreciation of the results. 
plot_richness_tract <- create_boxplot(alpha_div_tract, "Richness")
plot_richness_tract

plot_shannon_tract <- create_boxplot(alpha_div_tract, "Shannon")
plot_shannon_tract

plot_simpson_tract <- create_boxplot(alpha_div_tract, "Simpson")
plot_simpson_tract

### Finally, Plot them all together ! 
plot_list(gglist=list(plot_richness_tract, plot_shannon_tract, plot_simpson_tract),labels = NULL)

#Is there a real difference ? I'll do an ANOVA next to check. In the Simpson 
 #diversity, it seems that there is a little influence..


# 4.4. Statistical test --------------------------------------


### 4.4.1. Summarize the median ----
ad_median <- alpha_div_tract %>%
              group_by(Tract) %>%
              dplyr::summarise(median_Richness = median(Richness),
                               median_shannon = median(Shannon),
                               median_Simpson = median(Simpson))
ad_median

write.table(ad_median, file="./Output_Files_Seq_2_Methode_bacteria/Files/Stats_Files/Alpha_Diversity/Seq2_MB_alpha_div_median.tsv", sep="\t", quote=F, col.names=NA)

# Visually, I'm not sure, I'll see with the test, but first, are the data normal?


### 4.4.2. Linear model ----

# Creation of a linear model:
lm_rich <- lm(tract4G$Richness ~ tract4G$Bloc)
lm_sha <- lm(tract4G$Shannon ~ tract4G$Bloc)
lm_sim <- lm(tract4G$Simpson ~ tract4G$Bloc)

# Diagnostics graphs:
plot(lm_rich)
plot(lm_sha)
plot(lm_sim)

# Normality test for the linear model: 
### How it works: 
    #--> The hypothesis H0 is that the data are normal 
    #--> If the p-value is less than or equal to 0.05, we reject H0 
         #= The data are not normal 
    #--> Failing the test allow us to state with 95% confidence that the data 
         #doesn't fit the normal distribution.

shapiro.test(residuals(lm_rich))
#p-value = 0.3081 -> the data seems to follow a normal distribution. 
shapiro.test(residuals(lm_sha))
#p-value = 0,01435 -> the data doesn't fit the normal distribution.
shapiro.test(residuals(lm_sim))
#p-value = 0,0006244 -> the data doesn't fit the normal distribution.



### 4.4.3. Normality ----

# Before doing any test, we need look at the normality of each population. Here, 
 #I'm going to look at the normality of 2 tracts, 3 tracts, 4 tracts, galleries and all

# Creation of the subtable needed:
tract4 <- subset(alpha_div, !(Tract %in% c('2_tracts',"3_tracts","Gallery")))
tract3 <- subset(alpha_div, !(Tract %in% c('2_tracts',"4_tracts","Gallery")))
tract2 <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts","Gallery")))
galleries <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts","2_tracts")))

# Normality evaluation:

# Create an empty data frame to store the results
results_df <- data.frame(
  Sample = character(),
  Variable = character(),
  Normality = logical(),
  p_value = numeric(),
  After_Transformation_Normality = logical(),
  After_Transformation_p_value = numeric(),
  stringsAsFactors = FALSE
)

# Assuming you have loaded your data frames and they are named tract2, tract4, and tract3

samples <- list(tract2, tract3, tract4, alpha_div_tract, galleries)
sample_names <- c("tract2", "tract3", "tract4","alpha_div_tract", "galleries")

normality_results <- list()

for (i in seq_along(samples)) {
  sample <- samples[[i]]
  sample_name <- sample_names[i]
  
  normality_results[[sample_name]] <- list(before_transformation = list(), after_transformation = list())
  
  # Normality check loop - First test
  for (var_name in c("Richness", "Shannon", "Simpson")) {
    shapiro_test_result_before <- shapiro.test(sample[[var_name]])
    normality_results[[sample_name]]$before_transformation[[var_name]] <- list(normal = shapiro_test_result_before$p.value > 0.05, p_value = shapiro_test_result_before$p.value)
    
    # Apply log transformation if necessary
    if (!normality_results[[sample_name]]$before_transformation[[var_name]]$normal) {
      sample[[var_name]] <- log(sample[[var_name]])
    }
  }
  
  # Normality check after transformation - Second test
  for (var_name in c("Richness", "Shannon", "Simpson")) {
    shapiro_test_result_after <- shapiro.test(sample[[var_name]])
    normality_results[[sample_name]]$after_transformation[[var_name]] <- list(normal = shapiro_test_result_after$p.value > 0.05, p_value = shapiro_test_result_after$p.value)
    
    # Add the results to the data frame
    results_df <- rbind(results_df, data.frame(
      Sample = sample_name,
      Variable = var_name,
      Before_Transformation_Normal = normality_results[[sample_name]]$before_transformation[[var_name]]$normal,
      Before_Transformation_p_value = normality_results[[sample_name]]$before_transformation[[var_name]]$p_value,
      After_Transformation_Normal = normality_results[[sample_name]]$after_transformation[[var_name]]$normal,
      After_Transformation_p_value = normality_results[[sample_name]]$after_transformation[[var_name]]$p_value
    ))
  }
}

# Print the data frame
print(results_df)

## Conclusion: 
   # If the data are normal, use the parametric test ANOVA or t test
   # If not after the log transformation, use a non parametric test Mann-Whitney 
    #or Kruskal-Wallisfor the ANOVA, I'll use a parametric test for the specific 
               #richness and non parametric test for Shannon and Simpson index.

### 4.4.4. Statistical test ----

# Keep in mind: 
  ## if data are normal: 
     #- 2 groups: t test
     #- 3 groups or more: ANOVA test
  ## if data are not normal:
     #- 2 groups: Mann-Whitney = wilcoxon
     #- 3 groups or more: Kruskal-Wallis 
    #Mann-Whitney = wilcoxon - 2 groups

#### Bloc vs tracts ----

# Create an empty data frame to store the results
statistical_results_bloc <- data.frame(
  Sample = character(),
  Variable = character(),
  Test = character(),
  Comparison = character(),
  Test_Statistic = numeric(),  # Add a column for test statistic
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each row of results_df
for (i in 1:nrow(results_df)) {
  # Extract relevant information from results_df
  sample_name <- results_df$Sample[i]
  variable <- results_df$Variable[i]
  before_transformation_normal <- results_df$Before_Transformation_Normal[i]
  after_transformation_normal <- results_df$After_Transformation_Normal[i]
  p_value <- results_df$Before_Transformation_p_value[i]
  
  # Extract the sample data
  sample_data <- switch(sample_name, 
                        "tract2" = tract2,
                        "tract3" = tract3,
                        "tract4" = tract4,
                        "galleries" = galleries)
  
  # Check if the data are normal and do not need transformation
  if (before_transformation_normal) {
    # Use ANOVA or t test
    if (length(unique(sample_data$Bloc)) == 2) {
      # Use t test
      test_result <- t.test(sample_data[[variable]] ~ sample_data$Bloc)
      comparison <- paste(names(table(sample_data$Bloc))[1], "vs", names(table(sample_data$Bloc))[2])
      test_type <- "t-test"
    } else {
      # Use ANOVA
      test_result <- aov(sample_data[[variable]] ~ sample_data$Bloc)
      comparison <- "ANOVA"
      test_type <- "ANOVA"
    }
  } else {
    # Check if the transformed data are normal
    if (after_transformation_normal) {
      # Use transformed data and ANOVA or t test
      transformed_data <- log(sample_data[[variable]])
      if (length(unique(sample_data$Bloc)) == 2) {
        # Use t test
        test_result <- t.test(transformed_data ~ sample_data$Bloc)
        comparison <- paste(names(table(sample_data$Bloc))[1], "vs", names(table(sample_data$Bloc))[2])
        test_type <- "t-test"
      } else {
        # Use ANOVA
        test_result <- aov(transformed_data ~ sample_data$Bloc)
        comparison <- "ANOVA"
        test_type <- "ANOVA"
      }
    } else {
      # Use data before transformation and Mann-Whitney or Kruskal-Wallis test
      if (length(unique(sample_data$Bloc)) == 2) {
        # Use Mann-Whitney test
        test_result <- wilcox.test(sample_data[[variable]] ~ sample_data$Bloc)
        comparison <- paste(names(table(sample_data$Bloc))[1], "vs", names(table(sample_data$Bloc))[2])
        test_type <- "Mann-Whitney"
      } else {
        # Use Kruskal-Wallis test
        test_result <- kruskal.test(sample_data[[variable]] ~ sample_data$Bloc)
        comparison <- "Kruskal-Wallis"
        test_type <- "Kruskal-Wallis"
      }
    }
  }
  
  # Extract test statistic or t-value
  if (test_type %in% c("t-test", "Mann-Whitney")) {
    test_statistic <- test_result$statistic
  } else {
    test_statistic <- NA  # For ANOVA or Kruskal-Wallis
  }
  
  # Store the test results in the statistical_results data frame
  statistical_results_bloc <- rbind(statistical_results_bloc, 
                               data.frame(Sample = sample_name, 
                                          Variable = variable, 
                                          Test = test_type, 
                                          Comparison = comparison, 
                                          Test_Statistic = test_statistic,  # Include test statistic
                                          p_value = test_result$p.value,
                                          stringsAsFactors = FALSE))
}

# Print the statistical results
print(statistical_results_bloc)

#### Tracts vs tracts ----

# Create an empty data frame to store the results
statistical_results_tracts <- data.frame(
  Sample1 = character(),
  Sample2 = character(),
  Variable = character(),
  Test = character(),
  Comparison = character(),
  Test_Statistic = numeric(),  # Add a column for test statistic
  p_value = numeric(),
  stringsAsFactors = FALSE
)

sample_levels <- c("tract2", "tract3", "tract4", "galleries")

# Loop through each pair of levels
for (i in 1:(length(sample_levels)-1)) {
  for (j in (i+1):length(sample_levels)) {
    # Extract relevant information
    sample1_name <- sample_levels[i]
    sample2_name <- sample_levels[j]
    
    # Extract the sample data
    sample1_data <- switch(sample1_name, 
                           "tract2" = tract2,
                           "tract3" = tract3,
                           "tract4" = tract4,
                           "galleries" = galleries)
    sample2_data <- switch(sample2_name, 
                           "tract2" = tract2,
                           "tract3" = tract3,
                           "tract4" = tract4,
                           "galleries" = galleries)
    
    # Loop through each variable
    for (variable in c("Richness", "Shannon", "Simpson")) {
      # Check normality for sample1 before transformation
      normality_sample1_before <- results_df[results_df$Sample == sample1_name & 
                                               results_df$Variable == variable, "Before_Transformation_Normal"]
      # Check normality for sample1 after transformation
      normality_sample1_after <- results_df[results_df$Sample == sample1_name & 
                                              results_df$Variable == variable, "After_Transformation_Normal"]
      
      # Check normality for sample2 before transformation
      normality_sample2_before <- results_df[results_df$Sample == sample2_name & 
                                               results_df$Variable == variable, "Before_Transformation_Normal"]
      # Check normality for sample2 after transformation
      normality_sample2_after <- results_df[results_df$Sample == sample2_name & 
                                              results_df$Variable == variable, "After_Transformation_Normal"]
      
      # Perform appropriate test based on normality
      if (!is.na(normality_sample1_before) && !is.na(normality_sample1_after) && 
          !is.na(normality_sample2_before) && !is.na(normality_sample2_after)) {
        if ((normality_sample1_before || normality_sample1_after) && 
            (normality_sample2_before || normality_sample2_after)) {
          # Use t-test
          test_result <- t.test(sample1_data[[variable]], sample2_data[[variable]])
          test_type <- "t-test"
        } else {
          # Use Mann-Whitney U test
          test_result <- wilcox.test(sample1_data[[variable]], sample2_data[[variable]])
          test_type <- "Mann-Whitney"
        }
      } else {
        # If normality information is missing, default to Mann-Whitney U test
        test_result <- wilcox.test(sample1_data[[variable]], sample2_data[[variable]])
        test_type <- "Mann-Whitney"
      }
      
      # Extract test statistic or U value
      test_statistic <- ifelse(test_type == "t-test", test_result$statistic, test_result$statistic)
      
      # Store the test results in the statistical_results data frame
      statistical_results_tracts <- rbind(statistical_results_tracts, 
                                          data.frame(Sample1 = sample1_name, 
                                                     Sample2 = sample2_name,
                                                     Variable = variable, 
                                                     Test = test_type, 
                                                     Comparison = paste(sample1_name, "vs", sample2_name), 
                                                     Test_Statistic = test_statistic,  # Include test statistic
                                                     p_value = test_result$p.value,
                                                     stringsAsFactors = FALSE))
    }
  }
}

# Print the statistical results
print(statistical_results_tracts)


### All the p-value are above 0.05, we cannot reject the H0, meaning that, it 
#is 95% sure that the there is no difference between of the number of tract.
print(normality_sample1)
print(normality_sample2)

#### All tracts ---- 
#I didn't want to take anymore time to automatized this part so here we are:

## Richness: Normal - no transformation
all_tract_richness <- aov(Richness ~ Tract + Bloc + Tree, data = alpha_div_tract)
summary(all_tract_richness)

## Shannon: Normal - with transformation
sample_info_tract <- subset(sample_info, !(Status %in% c('Gallery')))
log_shannon <- log(alpha_div_tract$Shannon)
log_alpha <- cbind(log_shannon, sample_info_tract)

all_tract_shannon <- aov(log_shannon ~ Tract + Bloc + Tree, data = alpha_div_tract)
summary(all_tract_shannon)

## Simpson: Not normal
all_tract_simpson_Tract <- kruskal.test(Simpson ~ Tract, data = alpha_div_tract)
all_tract_simpson_Tract
all_tract_simpson_Tree <- kruskal.test(Simpson ~ Tree, data = alpha_div_tract)
all_tract_simpson_Tree
all_tract_simpson_Bloc <- kruskal.test(Simpson ~ Bloc, data = alpha_div_tract)
all_tract_simpson_Bloc

#### All sample ----
## Richness, Shannon and Simpson are not normal

# Richness
all_richness_Tract <- kruskal.test(Richness ~ Tract, data = alpha_div)
all_richness_Tract
dunnTest(Richness ~ Tract, data=alpha_div, method="bonferroni")
# Differences are with the galleries only.
all_richness_Tree <- kruskal.test(Richness ~ Tree, data = alpha_div)
all_richness_Tree
all_richness_Bloc <- kruskal.test(Richness ~ Bloc, data = alpha_div)
all_richness_Bloc

# Shannon
all_shannon_Tract <- kruskal.test(Shannon ~ Tract, data = alpha_div)
all_shannon_Tract
dunnTest(Shannon ~ Tract, data=alpha_div, method="bonferroni")
## Differences are with the galleries only.
all_shannon_Tree <- kruskal.test(Shannon ~ Tree, data = alpha_div)
all_shannon_Tree
all_shannon_Bloc <- kruskal.test(Shannon ~ Bloc, data = alpha_div)
all_shannon_Bloc

# Simpson
all_simpson_Tract <- kruskal.test(Simpson ~ Tract, data = alpha_div)
all_simpson_Tract
dunnTest(Simpson ~ Tract, data=alpha_div, method="bonferroni")
## Differences are with the galleries only.
all_simpson_Tree <- kruskal.test(Simpson ~ Tree, data = alpha_div)
all_simpson_Tree
all_simpson_Bloc <- kruskal.test(Simpson ~ Bloc, data = alpha_div)
all_simpson_Bloc

### 4.4.5 Power ----

### I'm gonna calculate the power which is 1-Beta = the probability that 
    #we accept H0 and that it's not correct.  
### To do so, I use GPower which needs 
    #--> Effect size: which need the standard deviation within each group 
         #(I guess the mean SD) : 
alpha_div_tract %>%
  group_by(Tract) %>%
  dplyr::summarise(sd_Richness = sd(Richness),
                   sd_shannon = sd(Shannon),
                   sd_Simpson = sd(Simpson))

#     Tract     sd_Richness sd_shannon sd_Simpson
#   2_tracts        45.4      0.650      0.317
#   3_tracts        47.5      1.27       0.806
#   4_tracts        61.3      0.884      0.467
#mean of 51.4

    #The mean of each group: 
alpha_div %>%
  group_by(Tract) %>%
  dplyr::summarise(mean_Richness = mean(Richness),
                   mean_shannon = mean(Shannon),
                   mean_Simpson = mean(Simpson))
# Tract    mean_Richness mean_shannon mean_Simpson
# 2_tracts          88           1.84         1.28
# 3_tracts          63.6         2.25         1.58
# 4_tracts          81.9         2.12         1.52

    #The size of each group:
tract4 <- subset(alpha_div, !(Tract %in% c('2_tracts',"3_tracts","Gallery")))
tract3 <- subset(alpha_div, !(Tract %in% c('2_tracts',"4_tracts","Gallery")))
tract2 <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts","Gallery")))
Galleries <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts","2_tracts")))
    
#--> The alpha error : 0.05
    #--> The total sample size : 23
    #--> The number of group : 3

#For the richness: the power is of 0.1150223 (11%)
#For the Shannon: 0.0798820 (7%)
#For the Simpson: 0.1566392 (15%)

### It is acceptable to accept the H0 if the power is under 20%. 
    #It also mean that it seems that O have enough replicate.
