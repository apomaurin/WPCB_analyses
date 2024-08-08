###
#Title : "Accumuation curve"
#Author : Apolline Maurin
# Date : 01/09/2023
###

# Intro -------------------------------------------------------------------

# The objective is to create a species accumulation curve or the number of species 
 #for a certain number of sampled sites or individuals.

# We'll use the function specaccum from the vegan package for the species 
  #accumulation adn ggplot to plot it. 

# Species accumulation curves (SAC) are used to compare diversity properties of 
  #community data sets using different accumulator functions. The classic method 
  #is "random" which finds the mean SAC and its standard deviation from random 
  #permutations of the data, or subsampling without replacement 
  #(Gotelli & Colwell 2001).

# Sources: 
  # -> https://rdrr.io/rforge/vegan/man/specaccum.html
  # -> https://search.r-project.org/CRAN/refmans/vegan/html/specpool.html

# 1. R SET UP -------------------------------------------------------------

# 1.1. Install the required packages -------------------

###Install packages if you don't have it already. otherwise, don't load it. 

install.packages("vegan")
install.packages("ggplot2")


# 1.2. Load the required packages ----------------------

library(vegan)
library(ggplot2)


# 1.3. Define the working directory --------------------
  
### First define where you're gonna work : 
setwd("D:/Thèse/7_Labo/Sequencage/Seq_2/Resultats_Seq_2/Methode_Bacteria")
### Then, check if you are where you want to be : 
getwd()
list.files()


# 2. Load the data --------------------------------------------------------

count_tab <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Count_Tab_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.matrix(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Tax_Tab.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("./Output_Files_Seq_2_Methode_Bacteria/Files/Stats_Files/Seq2_MB_Sample_info_No_Conta.tsv", header=T, row.names=1, check.names=F, sep="\t"))


# 3. Prepare the data -----------------------------------------------------

# 3.1. ASV in columns and names in rows -------------------------------

count_tab2 <- t(count_tab)
count_tab2 <- data.frame(Subject = rownames(count_tab2), count_tab2)  

# 3.2. filtered data frame --------------------------------------------
### To filter, I need to add the column subject in the dataframe:
ct_subject <- merge(count_tab2,sample_info, by = "Subject", all = TRUE)

samples.out <- rownames(count_tab2)
rownames(ct_subject) <- samples.out

# 3.2.1. Data frame with only tracts ----
ct_tracts <- subset(ct_subject, !(Status %in% c('Gallery')),
                    (-c(Subject, Tract, Tree, Status, Bloc)))

ct_2t <- subset(ct_subject, !(Tract %in% c('3_tracts', '4_tracts', 'Gallery')),
                (-c(Subject, Tract, Tree, Status, Bloc)))
ct_3t <- subset(ct_subject, !(Tract %in% c('2_tracts', '4_tracts','Gallery')),
                (-c(Subject, Tract, Tree, Status, Bloc)))
ct_4t <- subset(ct_subject, !(Tract %in% c('2_tracts', '3_tracts', 'Gallery')),
                (-c(Subject, Tract, Tree, Status, Bloc)))

# 3.2.2. Data frame with only gallery ----
ct_gal <- subset(ct_subject, !(Tract %in% c('2_tracts', '3_tracts', '4_tracts')),
                 (-c(Subject, Tract, Tree, Status, Bloc)))

# 4. Accumulation curve ---------------------------------------------------

## function details:
### The function "specaccum" finds species accumulation curves or the number of 
#species for a certain number of sampled sites or individuals.

### The Method "random" will add sites in a random order.

# 4.1. Tracts ---------------------------------------------------------

layout(matrix(c(1, 2, 3, 4), nrow = 2,
              ncol = 2, byrow = TRUE))

# 2t ----
sp1_2t <- specaccum(ct_2t)
sp2_2t <- specaccum(ct_2t, "random")
sp2_2t
summary(sp2_2t)


plot_2t <- plot(sp1_2t, ci.type="poly", col="tan1", lwd=2, ci.lty=0, ci.col="moccasin", xlab = "Number of tree", ylab = "Species", cex.axis = 1, cex.lab = 1.5)
plot_2t <- boxplot(sp2_2t, col="olivedrab1", add=TRUE, pch="+") 


par(mgp=c(4,1,0))
par(mar=c(6,6,4,1) + 0.1)

save.image(plot_2t, file="./Output_Files_Seq_2_Methode_bacteria/Figures/Accumulation_curve")

mtext("A", side = 3, at = -25, line = 1, cex = 1, cex.axis = 4)


# 3t ----
sp1_3t <- specaccum(ct_3t)
sp2_3t <- specaccum(ct_3t, "random")
sp2_3t
summary(sp2_3t)

plot_3t <- plot(sp1_3t, ci.type="poly", col="tan1", lwd=2, ci.lty=0, ci.col="moccasin", xlab = "Number of tree", ylab = "Species", cex.axis = 1, cex.lab = 1.5)
plot_3t <- boxplot(sp2_3t, col="olivedrab1", add=TRUE, pch="+") 


ggplot(data = sp1_2t_df) +
  geom_point(mapping = aes(x = sites, y = richness))


sp1_2t_df <- cbind(sites = sp1_2t$sites, richness = sp1_2t$richness)

mtext("B", side = 3, at = -25, line = 1, cex = 1)

# 4t ----
sp1_4t <- specaccum(ct_4t)
sp2_4t <- specaccum(ct_4t, "random")
sp2_4t
summary(sp2_4t)

plot_4t <- plot(sp1_4t, ci.type="poly", col="tan1", lwd=2, ci.lty=0, ci.col="moccasin", xlab = "Number of tree", ylab = "Species", cex.axis = 1, cex.lab = 1.5)
plot_4t <- boxplot(sp2_4t, col="olivedrab1", add=TRUE, pch="+") 

mtext("C", side = 3, at = -25, line = 1, cex = 1)

# Load necessary library
library(vegan) # For specaccum function

# Set up the plotting area
par(mgp=c(4,1,0))
par(mar=c(6,6,4,1) + 0.1)

# Plot for 2_tracts
sp1_2t <- specaccum(ct_2t)
sp2_2t <- specaccum(ct_2t, "random")

plot(sp1_2t, ci.type="poly", col="tan1", lwd=2, ci.lty=0, ci.col="moccasin", 
     xlab = "Number of trees", ylab = "Species", cex.axis = 1, cex.lab = 1.5,
     main = "Accumulation Curves")

# Add the boxplots for 2_tracts
boxplot(sp2_2t, col="olivedrab1", add=TRUE, pch="+", at=sp1_2t$sites)

# Add additional accumulation curves for 3_tracts and 4_tracts
sp1_3t <- specaccum(ct_3t)
sp2_3t <- specaccum(ct_3t, "random")
sp1_4t <- specaccum(ct_4t)
sp2_4t <- specaccum(ct_4t, "random")

# Add the 3_tracts curve
lines(sp1_3t$sites, sp1_3t$species, type="l", col="blue", lwd=2)
boxplot(sp2_3t, col="cyan", add=TRUE, pch="+", at=sp1_3t$sites)

# Add the 4_tracts curve
lines(sp1_4t$sites, sp1_4t$species, type="l", col="green", lwd=2)
boxplot(sp2_4t, col="yellow", add=TRUE, pch="+", at=sp1_4t$sites)

# Add legend to the plot
legend("bottomright", legend=c("2_tracts", "3_tracts", "4_tracts"),
       col=c("tan1", "blue", "green"), lwd=2, bty="n")

# Add subplot labels
mtext("A", side = 3, line = 1, cex = 1)

# Optionally, save the plot
# pdf("./Output_Files_Seq_2_Methode_bacteria/Figures/Accumulation_curves_combined.pdf")
# dev.off()


# 4.2. Gallery --------------------------------------------------------
# function details:
### -> Method random will add sites in a random order. 

sp1_gal <- specaccum(ct_gal)
sp2_gal <- specaccum(ct_gal, "random")
sp2_gal
summary(sp2_gal)

plot_gal <- plot(sp1_gal, ci.type="poly", col="tan1", lwd=2, ci.lty=0, ci.col="moccasin", xlab = "Number of tree", ylab = "Species", cex.axis = 1, cex.lab = 1.5)
plot_gal <- boxplot(sp2_gal, col="olivedrab1", add=TRUE, pch="+") 

par(mgp=c(4,1,0))
par(mar=c(6,6,4,1) + 0.1)


###Another way to have all the plot on one graph
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)


# Specaccum objects
sp1_2t <- specaccum(ct_2t)
sp2_2t <- specaccum(ct_2t, "random")
sp1_3t <- specaccum(ct_3t)
sp2_3t <- specaccum(ct_3t, "random")
sp1_4t <- specaccum(ct_4t)
sp2_4t <- specaccum(ct_4t, "random")
sp1_gal <- specaccum(ct_gal)
sp2_gal <- specaccum(ct_gal, "random")

# Function to convert specaccum objects to data frames
specaccum_to_df <- function(specaccum_obj, type) {
  data.frame(
    sites = specaccum_obj$sites,
    richness = specaccum_obj$richness,
    type = type,
    ci_lower = specaccum_obj$sd,
    ci_upper = specaccum_obj$sd 
  )
}


# Convert specaccum objects to data frames

df_2t <- specaccum_to_df(sp1_2t, "2 Tracts")
df_random_2t <- specaccum_to_df(sp2_2t, "2 Tracts Random")

df_3t <- specaccum_to_df(sp1_3t, "3 Tracts")
df_random_3t <- specaccum_to_df(sp2_3t, "3 Tracts Random")

df_4t <- specaccum_to_df(sp1_4t, "4 Tracts")
df_random_4t <- specaccum_to_df(sp2_4t, "4 Tracts Random")

df_gal <-specaccum_to_df(sp1_gal, "Gallery")
df_random_gal <- specaccum_to_df(sp1_gal, "Gallery Random")  

# Combine data frames
combined_df <- bind_rows(
  df_2t,
  df_3t,
  df_4t,
  df_gal
)

accumu_plot <-
  ggplot(combined_df, aes(x = sites, y = richness, color = type, fill = type, shape=type)) +
  geom_line() +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = richness - ci_lower, ymax = richness + ci_upper), alpha = 0.2) +
  labs(
    x = "Number of tree",
    y = "Number of species",
    title = "Species Accumulation Curves"
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

ggsave("D:/Thèse/14_Articles/1_Methode/Figures/Accumulation_curve_bacteria.jpeg", plot = accumu_plot, width = 12, height = 6, units = "in", dpi = 300)
ggsave("D:/Thèse/14_Articles/1_Methode/Figures/Accumulation_curve_bacteria.pdf", plot = accumu_plot, width = 12, height = 6, units = "in", dpi = 300)





# 5. Extrapolation --------------------------------------------------------

# The functions estimate the extrapolated species richness in a species pool, 
  #or the number of unobserved species. Function specpool is based on incidences 
  #in sample sites, and gives a single estimate for a collection of sample sites 
  #(matrix).

# 5.1. Tracts ---------------------------------------------------------

# 2t ----
specpool(ct_2t)
pool_2t <- poolaccum(ct_2t)
pool_2t
plot(pool_2t)

# 4tree : 53% according to the extrapolation 
# 5tree : 59%
# 6tree : 65 %
# 7 tree : 69%
# 9 tree : 76%

# 3t ----
specpool(ct_3t)
pool_3t <- poolaccum(ct_3t)
pool_3t
plot(pool_3t)

# 4tree : 48% according to the extrapolation 
# 5tree : 56%
# 6tree : 62%
# 7 tree : 67%

# 4t ----
specpool(ct_4t)
pool_4t <- poolaccum(ct_4t)
pool_4t
plot(pool_4t)

# 4tree : 54% according to the extrapolation 
# 5tree : 60%
# 6tree : 66%
# 7 tree : 72%

# 5.2. Gallery --------------------------------------------------------
specpool(ct_gal)
pool_gal <- poolaccum(ct_gal)
pool_gal
plot(pool_gal)

# 4tree : 87% according to the extrapolation 
# 5tree : 90%
# 6tree : 92%
# 7 tree : 93%
# 9 tree : 95%
# 11 tree : 96%

# I would say that it's not possible to extrapolate what I have obtain from the 
  #gallery to the insect. It seems that the type of bacteria found in the 
  #galleries are more stable throughout the plantation than the insect's tracts.

# However, between 2, 3 of 4 tracts, I don't have that much differences. 
  #And around 7 trees, it seems that we have around 70% of the diversity 
  #and at 9 trees almost 80%


# 5. How many more ? ------------------------------------------------------

### To do so, I'll use the Excel table give by Chao and al.(2009) in the 
    #supplemental data of the article "Sufficient sampling for asymptotic 
    #minimum species richness estimators".

### To use it, I need to know how many time an ASV is count one or two time
    #= number of singleton and doubleton.The following code will help up know that.

# Count the number of samples in which each ASV is present
  #Change the name for rowSums with your count table of interest
presence_count <- rowSums(t(ct_gal) > 0)

# Count the number of ASVs that are present in either one or two samples
count_one_samples <- sum(presence_count == 1)
count_two_samples <- sum(presence_count == 2)

# Print the count
print(count_one_samples)
print(count_two_samples)

### If you want to chek the lines where the singleton and the doubleton are found

# Identify the rows where ASVs have counts of one or two
which_one_samples <- which(presence_count == 1)
which_two_samples <- which(presence_count == 2)

# Print the row numbers and corresponding ASV names or identifiers
cat("ASVs with counts of one or two are found in rows:", which_one_samples, "\n")
cat("Corresponding ASV names or identifiers:", rownames(t(ct_2t))[which_one_samples], "\n")

cat("ASVs with counts of one or two are found in rows:", which_two_samples, "\n")
cat("Corresponding ASV names or identifiers:", rownames(t(ct_2t))[which_two_samples], "\n")

### For 2t: 
    #-> 114 singleton (f1)
    #-> 69 doubleton (f2)
## According to the Excel spreadsheet, in order to successfully have all the 
  #communities, I need an addition of 57 samples (total of 66 samples)
## To have 95% of th communities, I would need 12 more samples (total of 21 samples)
## To have 98% of th communities, I would need 17 more samples (total of 26 samples)
## It seems I have already 76% of the communities ?

### It is important to note that When There is only a few sample (like I have),
    #The estimation is often exaggerated. But It do give us an idea. 

### For the galleries
    #-> 25 singleton (f1)
    #-> 31 doubleton (f2)
## According to the Excel spreadsheet, in order to successfully have all the 
#communities, I need an addition of 23 samples (total of 34 samples)
## To have 98% of th communities, I would need 1 more samples (total of 12 samples)
## It seems that I already have 97% of the communities
    

#Insect : around 22 trees (mean of MB and MF)
#Tree : around 17 trees (Mean of MB and MF)




