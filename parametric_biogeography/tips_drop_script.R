##It is important and recoomended to use a dated tree for biogeographic analysis because it doesn;t make much sense to trace the
#spatial or distributional history of a species (or popuation) without considering the time-scale. And of course, geographic 
#distribtuions are so often shaped by geoloigcal features, which in turn change with time. 

#This script is a small tool in the process of paramteric biogeographic analysis. Usually, phylogenetic trees are calulated with outgroups
#for the purposes of rooting. While using the tree for biogeographic analysis, the outgroup taxa may be irrelevant and you might want to drop them 
#fom your tree. This R script is written for this purpose: to drop tips that are not required from a phylogenetic tree. 

library(ape)
install.packages("diversitree")
library(diversitree)
library(phytools)
library(phylobase)
library(phangorn)
########################################################################
#If you still want to compile from sources, you could first install 
#the fftw system library from the command line

#brew install fftw

#and then run 
install.packages("fftwtools")
########################################################################
### tree manupulations
########################################################################
# package "ape"
tree <- read.nexus("best tree")
tree
tree$tip.label
par(mfcol=c(1,1))
plot(tree)
plot(tree, cex=0.5, show.node.label=TRUE)
# zoom(tree, 80:110)

#dropping taxa that haves indices 48 and 49 in the tree$tip.label vector. 

?drop.tip #need to get rid of 2 taxa
tip <- c(48, 49) # excluding taxa 48, 49
tr_droptip <- drop.tip(tree, tip, trim.internal = TRUE)
plot(tr_droptip, cex=0.7)

# write droptip tree as .nex and .tre
write.tree(tr_droptip, file="reduced_species_2_c90s-5_2_RAxML_tips_dropped.tre")
write.nexus(tr_droptip, file="reduced_species_2_c90s-5_2_RAxML_tips_dropped.nex")


######################################################################
#Checking if the tree is bifurcating and ultrametric or not
######################################################################
tree <- read.newick("reduced_species_2_c90s-5_2_RAxML_tips_dropped.tre")
## ultrametric or not
is.ultrametric(tree)
## binary or not
is.binary(tree)
#####################################################################
