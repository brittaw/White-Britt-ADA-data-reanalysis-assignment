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
#m_sp_col <- read_csv("Mahler_species_col.csv", col_names = TRUE)
#traits_2 <- read_csv("Traits_2.csv", col_names = TRUE)
#traits_reduced <- m_sp_col %>% left_join(traits_2)

###load NEXUS tree file - they specify this NEXUS tree in their README file 
tree <- read.nexus("5Loc_25milburn_100kthin.MCC.median.nex")

###now we will remove missing data from the trait morph
#trait70sp <- traits_reduced %>% 
#semi_join(tip_label, c("species" = "tree$tip.label"))

###Nevermind, found this in supplements: 
  #"Thus, for our comparative analyses, measurements of A. guazuma correspond to A. garridoi in the
  #phylogeny, and measurements of A. centralis correspond to A. terueli in the phylogeny"
traits_3 <- read_csv("Traits_3.csv", col_names = TRUE)
tip_label <-tibble(tree$tip.label)
trait72sp <- traits_3 %>% 
  semi_join(tip_label, c("species" = "tree$tip.label"))

#get a list of taxa missing
missing_tips <- tip_label %>% 
  anti_join(trait72sp, c("tree$tip.label" = "species"))

#create list
miss <- c("Leiocephalus", "Polychrus", "acutus", "aeneus", "carpenteri", "ferreus", "fuscoauratus", 
           "gingivinus", "griseus", "heterodermus", "leachii", "limifrons", "luciae", "marmoratus", 
           "tandai", "ortonii", "pogus", "punctatus", "rejectus", "roquet", "schwartzi", "trachyderma", 
           "transversalis", "trinitatis", "wattsii")

##Now the tree has many species the data file does not
tree<-(drop.tip(tree, miss))
plot(tree)

###load necessary package for phylogenetic PCA
library(phytools)

#rename column sepcies to match tree
trait72sp <- trait72sp %>% rename("tree$tip.label" = species)

#maybe this to fix error below 
X<-as.matrix(trait72sp, row.names=1)

###build the model of phylogenetic relationships to the covariance matrix of traits to determine if 
### covariance occurs at a rate greater than expected by Brownian Motion

pPCA <- phyl.pca(tree, X, method="BM", mode="cov")

##visualize in screeplot the PCAs
screeplot(pPCA)
scatter(pPCA, useLag=TRUE)
plot(pPCA, useLag=TRUE)

##Calculate the Euclidean distance to every other species in by the first four axes of a phylogenetic PCA on traits. 
RePCA <- read_csv("Scores.retained.csv", col_names = TRUE)
as_matrix(RePCA)
library(philentropy)

# compute the Euclidean Distance with default parameters
distance(RePCA, method = "euclidean")

##Conducted a 7-category linear discriminant analysis (LDA) on traditional ecomorph class for the four cryptic giants in our sample (A. landestoyi, A. chamaeleonides, A. guamuhaya, A. porcus) plus all species assigned to one of the six traditional ecomorph classes
