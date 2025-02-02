#I used the package ape in R to reconstruct ancestral discrete states. My analysis involved a discrete binary character (with two states 0 and 1)
#There are three models for discrete character evolution: equal rates or ER model (rate of transition from 0 to 1 is same as that from 1 to 0), all rates
#different or ARD model (rate of transition from 0 to 1 not equal to the rate of transition from 1 to 0). Now, for a binary character, these two models 
#are enough (k(k-1) parameters where k is the number of discrete states). However, if you have a character that has three or more states, the symmetric model or SYM model 
#can be implemented as well. If my memory serves right, for a binary character, this model is the same as ER model. 
#All models are implemented in the maximum likelihood framework. Check out the documentation of the ape package for more functions. 


install.packages("ape")
library(ape)
#ancestral state reconstruction
tree = read.nexus("reduced_species_2_c90s-5_2_RAxML_tips_dropped.nex")
tree
tree$tip.label
plot(tree)
characters = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1,1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
dataf = data.frame(tree$tip.label, characters)
length(characters)
names(characters) = tree$tip.label
character
characters
ER = ace(characters, tree, type = "discrete", model = "ER")
SYM = ace(characters, tree, type = "discrete", model = "SYM")
ARD = ace(characters, tree, type = "discrete", model = "ARD")
ER$loglik
SYM$loglik
ARD$loglik
ER$rates
SYM$rates
ARD$rates
1-pchisq(2*abs(ER$loglik - ARD$loglik), 1)
#plot tree
install.packages("corHMM")
library(corHMM)
plotRECON(tree, ER$lik.anc, piecolors = NULL, cex = 0.5, pie.cex=0.75, file = NULL)
ER$lik.anc
plotRECON(tree, ARD$lik.anc, piecolors = NULL, cex = 0.5, pie.cex=0.75, file = NULL)
ARD$lik.anc
ARD$rates

#References
#Paradis, E., Blomberg, S., Bolker, B., Brown, J., Claude, J., Cuong, H. S., ... & Didier, G. (2019). Package ‘ape’. 
#Analyses of phylogenetics and evolution, version, 2(4), 47.
