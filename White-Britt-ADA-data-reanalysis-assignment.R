#Data reanalysis of Mahler et al. 2016 

#Ecomorphological Affinities
##Goal: Quantify the ecomorphological similarity of new species to other anoles
###Using pPCA Revell (2009) on a covariance matrix of log-transformed and averaged data for body size and 
###shape residuals with the Bayesian MCC chronogram they calculated independent "size" and "shape" axes.  
###load log-transformed and averaged data for body size and shape data 
df <- "Traits.csv"

###load necessary packages for data tidying and transformation
library(tidyverse)

###rename dataframe and read csv into a tibble
traits <- read_csv(df, col_names = TRUE)

###see the variables in the tibble
glimpse(traits)

###The data file they provided in Dryad contains 29 more species than the data table they produced so...
m_sp_col <- read_csv("Mahler_species_col.csv", col_names = TRUE)
traits_2 <- read_csv("Traits_2.csv", col_names = TRUE)
traits_reduced <- m_sp_col %>% left_join(traits_2)

###load NEXUS tree file - they specify this NEXUS tree in their README file 
tree <- read.nexus("5Loc_25milburn_100kthin.MCC.median.nex")

###The number of taxa in the data file and tree do not match
tip_label <-tibble(tree$tip.label)

dropped_tips <- tip_label %>% 
anti_join(traits_reduced, c("tree$tip.label" = "species"))

###the above gives a dataset of 27 tips to drop from an OG data set of 97, leaving 70. Two sp are now missing in the morph set.
missing_dropped_tips <- traits_reduced %>% 
  anti_join(tip_label, c("species" = "tree$tip.label"))

###centralis and guazuma appear to not be in the tree files but are reported in the loadings table for the pPCA in the paper.
###With this in mind I suppose I will continue with the analyses without the morph data for those species.

#tree<-(drop.tip(tree, dropped_tips))
#plot(tree)

###now we will remove missing data from the trait morph
trait70sp <- traits_reduced %>% 
semi_join(tip_label, c("species" = "tree$tip.label"))

###Nevermind, found this in supplements: 
  #"Thus, for our comparative analyses, measurements of A. guazuma correspond to A. garridoi in the
  #phylogeny, and measurements of A. centralis correspond to A. terueli in the phylogeny"
traits_3 <- read_csv("Traits_3.csv", col_names = TRUE)
traits_reduced <- m_sp_col %>% left_join(traits_3)
trait72sp <- traits_3 %>% 
  semi_join(tip_label, c("species" = "tree$tip.label"))

###load necessary package for phylogenetic PCA
library(phytools)

###build the model of phylogenetic relationships to the covariance matrix of traits to determine if 
### covariance occurs at a rate greater than expected by Brownian Motion
pPCA <- phyl.pca(tree, trait70sp, method="BM", mode="cov")

##visualize in screeplot the PCAs
screeplot(pPCA)
scatter(pPCA, useLag=TRUE)
plot(pPCA, useLag=TRUE)

###Step 1. Calculate the Euclidean distance to every other species in a four-dimensional morphospace of traits.

###Step 2. Test whether the new species forms a distinct phenotypic cluster with a specific clade (Chamaeleinides) by conducting a linear discriminant analysis (LDA) with 7 categories (1 represents the new species, 1 Chamaeleonides, then the remaining ecomorphs).  