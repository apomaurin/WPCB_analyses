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

# 1.2. Load the required packages ----------------------

library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(dendextend) ;packageVersion("dendextend")
library(devtools)
library(iNEXT)
library(ggbreak);packageVersion("ggbreak")
library(aplot);packageVersion("aplot")
library(ggsignif)

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
write.table(ad_DataInfo, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_ad_DataInfo.tsv", sep="\t", quote=F, col.names=NA)
                                 
## Estimation
ad_mid <- alpha_div$iNextEst

# Size_based
ad_est_size <- ad_mid$size_based
write.table(ad_est_size, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_ad_est_size.tsv", sep="\t", quote=F, col.names=NA)

# Coverage_based
ad_est_cover <- ad_mid$coverage_based
write.table(ad_est_cover, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_ad_est_cover.tsv", sep="\t", quote=F, col.names=NA)

## Alpha diversity
ad_solo <- alpha_div$AsyEst
write.table(ad_solo, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_AsyEst.tsv", sep="\t", quote=F, col.names=NA)


## Organizing the alpha_diversity data frame

ad_solo_filter <- ad_solo[,-c(4,5,6,7)] #only keep the name and value
ad_only <- pivot_wider(ad_solo_filter, names_from = "Assemblage", 
                       values_from = "Observed") #order by diversity measurement
ad_reverse <- as.data.frame(t(ad_only)[-1,]) #turn the table
colnames(ad_reverse) <- c("Richness","Shannon", "Simpson") #put the name you like
alpha_diversity <- data.frame(Subject = rownames(ad_reverse), ad_reverse) #add a column with the samples name

write.table(alpha_diversity, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_alpha_diversity.tsv", sep="\t", quote=F, col.names=NA)


# 4.3. Alpha diversity visualization -------------------------

# Plotting the alpha diversity ----
### Loading data

setwd("D:/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Fungi")
### Then, check if you are where you want to be : 
getwd()
list.files()

alpha_div <- read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_alpha_diversity.tsv", header=T, row.names=1, check.names=F, sep="\t")
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Seq2_MF_Sample_info_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))

alpha_div_all <- read.table("./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MFall_alpha_diversity.tsv", header=T, row.names=1, check.names=F, sep="\t")

### Data frame with all information ----

alpha_div <- merge(alpha_div,sample_info, by = "Subject", all = TRUE)
alpha_div_all <- merge(alpha_div_all,sample_info, by = "Subject", all = TRUE)

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
                y_position = 70)
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

plot_shannon <- create_boxplot(alpha_div, "Shannon")
plot_shannon

plot_simpson <- create_boxplot(alpha_div, "Simpson")
plot_simpson

### Finally, Plot them all together ! 
plot_alpha_div <- plot_list(gglist=list(plot_richness, plot_shannon, plot_simpson),labels = NULL)
plot_alpha_div
ggsave("D:/Thèse/14_Articles/1_Methode/Figures/Alpha_Diversity_fungi2.pdf", plot = plot_alpha_div, width = 12, height = 6, units = "in", dpi = 300)

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

#Old code, just in case
#alpha_div %>%
#gather(key = metric, value = value, c("Richness", "Shannon", "Simpson"))  %>%
#mutate(metric = factor(metric, levels = c("Richness", "Shannon", "Simpson"))) %>%
#  ggplot(aes(x = Tract, y = value)) +
#  geom_boxplot() +
#  geom_jitter(aes(color = Tract), height = 0, width = .2) +
#  labs(x = "", y = "") +
#  facet_wrap(~ metric, scales = "free") +
#  theme(legend.position="none")


# 4.4. Statistical test --------------------------------------

### 4.4.1. Summarize the median ----
ad_median <- alpha_div_all %>%
              group_by(Tract) %>%
              dplyr::summarise(median_Richness = median(Richness),
                               median_shannon = median(Shannon),
                               median_Simpson = median(Simpson))
ad_median

write.table(ad_median, file="./Output_Files_Seq_2_Methode_Fungi/Files/Stats_Files/Alpha_Diversity/Seq2_MF_alpha_div_median.tsv", sep="\t", quote=F, col.names=NA)

# Visually, I'm not sure, I'll see with the test, but first, are the data normal?
alpha_div2 <- subset(alpha_div, !(Tract %in% c('3_tracts','4_tracts')))
alpha_div3 <- subset(alpha_div, !(Tract %in% c('4_tracts', '2_tracts')))
alpha_div4 <- subset(alpha_div, !(Tract %in% c('3_tracts','2_tracts')))


### 4.4.2. Linear model ----

# Creation of a linear model:
lm_rich <- lm(alpha_div$Richness ~ alpha_div$Tract)
lm_sha <- lm(alpha_div$Shannon ~ alpha_div$Tract)
lm_sim <- lm(alpha_div$Simpson ~ alpha_div$Tract)

# Diagnostics graphs:
plot(lm_rich)
plot(lm_sha)
plot(lm_sim)

# Normality test: 
### How it works: 
    #--> The hypothesis H0 is that the data are normal 
    #--> If the p-value is less than or equal to 0.05, we reject H0 
         #= The data are not normal 
    #--> Failing the test allow us to state with 95% confidence that the data 
         #doesn't fit the normal distribution.

shapiro.test(residuals(lm_rich))
#p-value = 0.6003 -> the data seems to follow a normal distribution. 
shapiro.test(residuals(lm_sha))
#p-value = 0.5015 -> the data seems to follow a normal distribution. 
shapiro.test(residuals(lm_sim))
#p-value = 0.5332 -> the data seems to follow a normal distribution. 

## Conclusion: for the ANOVA, I'll use a parametric test for the all the 
   #Hills numbers

### 4.4.3. Normality ----

# Before doing any test, we need look at the normality of each population. Here, 
#I'm going to look at the normality of 2 tracts, 3 tracts, 4 tracts, galleries and all

# Creation of the subtable needed:
tract4 <- subset(alpha_div, !(Tract %in% c('2_tracts',"3_tracts")))
tract3 <- subset(alpha_div, !(Tract %in% c('2_tracts',"4_tracts")))
tract2 <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts")))
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

samples <- list(tract2, tract4, alpha_div_tract, galleries, alpha_div)
sample_names <- c("tract2", "tract4","alpha_div_tract", "galleries", "alpha_div")

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

# Loop through each pair of levels
for (i in 1:(length(sample_levels)-1)) {
  for (j in (i+1):length(sample_levels)) {
    # Extract relevant information
    sample1_name <- sample_levels[i]
    sample2_name <- sample_levels[j]
    
    # Extract the sample data
    sample1_data <- switch(sample1_name, 
                           "tract2" = tract2,
                           "tract4" = tract4,
                           "galleries" = galleries)
    sample2_data <- switch(sample2_name, 
                           "tract2" = tract2,
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

### All tracts ---- 
#I didn't want to take anymore time to automatized this part so here we are:

## Richness: Normal - no transformation
all_tract_richness <- aov(Richness ~ Tract + Bloc + Tree, data = alpha_div_tract)
summary(all_tract_richness)

## Shannon: Normal - no transformation
all_tract_shannon <- aov(Shannon ~ Tract + Bloc + Tree, data = alpha_div_tract)
summary(all_tract_shannon)

## Shannon: Normal - no transformation
all_tract_simpson <- aov(Simpson ~ Tract + Bloc + Tree, data = alpha_div_tract)
summary(all_tract_simpson)


#### All sample ----

#Richness - no transformation:
anov_rich <- aov(Richness ~ Tract+Bloc+Tree, data = alpha_div)
summary(anov_rich)
TukeyHSD(anov_rich)
# It seems to be a significant differences between the gallery and 2 tracts (p = 0.39).

#Shannon - no transformation:
anov_sha <- aov(Shannon ~ Tract+Bloc+Tree, data = alpha_div)
summary(anov_sha)
TukeyHSD(anov_sha)

#Simpson - transformation:
log_simpson <- log(alpha_div$Simpson)
log_simpson <- cbind(log_simpson, sample_info)
anov_sim <- aov(log_simpson ~ Tract+Bloc+Tree, data = alpha_div)
summary(anov_sim) 
TukeyHSD(anov_sim)

### All the p-value are above 0.05, we cannot reject the H0, meaning that, it 
    #is 95% sure that the there is no difference between of the number of tract.

### Power ----

### I'm gonna calculate the power which is 1-Beta = the probability that 
    #with accept H0 and that it's not correct.  
### To do so, I use GPower :
    ## Test Family - F tests 
    ## Statistical test : ANOVA : Fixed effects, omnibus, one-way
    ## For the test, I'll need :
       #--> Effect size: which need the standard deviation within each group 
         #(I guess the mean SD) : 

alpha_div %>%
  group_by(Tract) %>%
  dplyr::summarise(sd_Richness = sd(Richness),
                   sd_shannon = sd(Shannon),
                   sd_Simpson = sd(Simpson))


count_tab_tract <- subset(alpha_div_all, !(Tract %in% c('Gallery')))
count_tab_tract <- subset(count_tab_tract, !(Tract %in% c('2_tracts')))
count_tab_tract <- subset(count_tab_tract, !(Subject %in% c('Temoin_1_champignons')))
count_tab_tract <- subset(count_tab_tract, !(Subject %in% c('Temoin_2_champignons')))


alpha_div %>%
  group_by(Tract) %>%
  dplyr::summarise(sd_Richness = sd(Richness),
                   sd_shannon = sd(Shannon),
                   sd_Simpson = sd(Simpson))

count_tab_tract <- select(alpha_div_tract, -contains("2_tracts"))
count_tab_tract <- subset(alpha_div_tract, !(Tract %in% c('2_tracts')))

       #The mean of each group: 
alpha_div %>%
  group_by(Tract) %>%
  dplyr::summarise(mean_Richness = mean(Richness),
                   mean_shannon = mean(Shannon),
                   mean_Simpson = mean(Simpson))
 
       #The size of each group:
tract4 <- subset(alpha_div, !(Tract %in% c('2_tracts',"3_tracts","Gallery", "Control")))
tract3 <- subset(alpha_div, !(Tract %in% c('2_tracts',"4_tracts","Gallery", "Control")))
tract2 <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts","Gallery", "Control")))
Gallery <- subset(alpha_div, !(Tract %in% c('4_tracts',"3_tracts","2_tracts","Control")))

       #--> The alpha error : 0.05
       #--> The total sample size : 25 (or 23)
       #--> The number of group : 5 (or 4)

#Between tract : 
#For the richness: the power is of 0.3040084 (30%)
#For the Shannon: 0.0908675 (9%) (without control)
#For the Simpson: 0.0624725 (6%) (without control)

### It is acceptable to accept the H0 if the power is under 20%. 
    #It also mean that it seems that we have enough replicate.

## The power is way over 20%, not reliable ? To much variaiton in the sd ? 

#Tract-Galleries
##Richness : 
   ### Criterion : alpha = 0.000741986 = 0.07%
   ### Post-hoc :beta = 0.0004336 (1-B = 0.9999563) = 0.04% 
## Shanon : 
   ### Criterion : alpha = 0.9187290 = 92% 
   ### Post-hoc :beta = 0.8918648 (1-B = 0.108135) = 89% 
## Simpson4 : 
   ### Criterion : alpha = 0.9152055 = 92%
   ### Post-hoc :beta = 0,8860504 (1-B = 0.1139496) = 88% 

## Beta, not good, except for the richness... 


##Pie chart ---- 
##data

count_tab_tract <- select(count_tab, -contains("_G"))
sample_info_tract <- subset(sample_info, !(Status %in% c('Gallery')))

##transformation of the data
data_stand <- decostand(count_tab, method = "hellinger")

count_tabr <- t(data_stand)
count_tabr_subj_tot <- data.frame(Subject = rownames(count_tabr), count_tabr, row.names=NULL)

taxa_tab_asv <- data.frame(ASV = rownames(taxa_tab), taxa_tab, row.names=NULL)

### then, create the data.frame:

metadata <- sample_info %>%
  select(Subject, Status)
  
  otu_counts <- count_tabr_subj_tot %>%
  select(Subject, starts_with("ASV")) %>%
  pivot_longer(-Subject, names_to="ASV", values_to = "count")

taxonomy <- taxa_tab_asv %>%
  drop_na(Family)

##Compute percentage

otu_rel_abund <- inner_join(sample_info, otu_counts, by="Subject") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(Subject) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon") %>%
  mutate(Status = factor(Status,
                         levels=c("Tract",
                                  "Gallery")))
otu_rel_abund

taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="Family", ) %>%
  group_by(Status, Subject, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Status, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))


taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 3,
            mean = mean(mean_rel_abund),
            .groups="drop")

inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Status, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(x=Status, y=mean_rel_abund, fill=taxon)) +
  geom_col()+
  scale_fill_manual(name=NULL,
                    breaks=c("*Nectriaceae*", "*Chionosphaeraceae*",
                             "*Helotiales_fam_Incertae_sedis*", "*Chrysozymaceae*", 
                             "*Mytilinidiaceae*", "*Cucurbitariaceae*",
                             "*Malasseziaceae*", "*Abrothallaceae*",
                             "*Colacogloeaceae*", "*Ceratobasidiaceae*",
                             "*Tubeufiaceae*", "*Hamatocanthoscyphaceae*",
                             "*Serendipitaceae*","Other"),
                    values = c(brewer.pal(10, "RdBu"),brewer.pal(3, "OrRd") , "gray")) +
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

