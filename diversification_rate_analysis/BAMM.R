
library(ape)
library(BAMMtools)
tree = read.tree("reduced_species_tree.tre")
is.ultrametric(tree)
is.binary(tree)
plot(tree)
min(tree$edge.length)
tree$edge.length
tree$Nnode
plot.phylo(tree, type = "phylogram", use.edge.length = TRUE, show.node.label = TRUE)
install.packages("BAMMtools")
setBAMMpriors(tree)
setBAMMpriors(read.tree("reduced_species_2_tree.tre"))
generateControlFile(file = "trial.text", type = "diversification")
#analysing BAMM output

tree <- read.tree("reduced_species_tree.tre")
plot(tree)
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)

#accessing whether MCMC has converged

mcmcout = read.csv("BAMM_complete_splitter.txt", header = T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart = floor(0.1*nrow(mcmcout))
postburn = mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
postburn$N_shifts

#how many rate shifts in our data? posterior probabilites of models

post_probs = table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
shift_probs = summary(edata)

#computing Bayes Factors for model comparisons

postfile = "BAMM_complete_splitter.txt"
bfmat = computeBayesFactors(postfile, expectedNumberOfShifts = 1, burnin = 0.1)
bfmat
plotPrior(postfile, expectedNumberOfShifts = 1)
summary(edata)

########mean phylorate plot to plot diversification rates#######

plot.bammdata(edata, lwd=2, legend = T)

#credible sets of shift configurations

css = credibleShiftSet(edata, expectedNumberOfShifts = 1, threshold = 5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

#####plotting rgw maximum a posteriori probability configuration#####

best = getBestShiftConfiguration(edata, expectedNumberOfShifts = 1)
edata
plot.bammdata(best, lwd = 2, legend = T)

##getting speciation rates across phylogeny

allrates = getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

##getting extinction rates across phylogeny

allrates = getCladeRates(edata)
mean(allrates$mu)
quantile(allrates$mu, c(0.05, 0.95))

#to get node numbers

plot.phylo(tree, show.node.label = T)
nodelabels()

#clade specific rates-clade 1 -- speciation

myr <- getCladeRates(edata, node = 54, nodetype = "include")
mean(myr$lambda)
quantile(myr$lambda, c(0.05, 0.95))

#clade specific extinction rates - clade1

myr <- getCladeRates(edata, node = 54)
mean(myr$mu)
quantile(myr$mu, c(0.05, 0.95))

#clade2 -- speciation rate

pachy <- getCladeRates(edata, node = 62, nodetype = "include")
mean(pachy$lambda)

quantile(pachy$lambda, c(0.05, 0.95))

#clade2 -- extinction rate

pachy <- getCladeRates(edata, node = 62, nodetype = "include")
mean(pachy$mu)

quantile(pachy$mu, c(0.05, 0.95))

#clade3 -- speciation rate
ban <- getCladeRates(edata, node = 65, nodetype = "include")
mean(ban$lambda)
quantile(ban$lambda, c(0.05, 0.95))

#clade3-- extinction rate

ban <- getCladeRates(edata, node = 65, nodetype = "include")
mean(ban$mu)
quantile(ban$mu, c(0.05, 0.95))

#rate through time analysis - overall

plotRateThroughTime(edata, ratetype = "speciation", intervals = seq(from = 0.05, 0.95, by = 0.01), ylim = c(0,1), yticks = 10)
text(x=30, y= 0.8, label="Overall speciation rate", font=4, cex=2.0, pos=4)


#rate through time for clade1

sp1 = "species6985"
sp2 = "species7239"
tipnode1 <- which(tree$tip.label == sp1)
tipnode2 <- which(tree$tip.label == sp2)
mrca <- getMRCA(tree, tip = c(tipnode1, tipnode2))
plotRateThroughTime(edata,node = 54, ratetype = "speciation", intervals = seq(from = 0.05, 0.95, by = 0.01), ylim = c(0,1), yticks = 10)
plotRateThroughTime(edata,node = 54, nodetype = "exclude", ratetype = "speciation", intervals = seq(from = 0.05, 0.95, by = 0.01), ylim = c(0,1), yticks = 10)
plot(edata,node = 54, ratetype = "speciation", intervals = seq(from = 0.05, 0.95, by = 0.01))
text(x=50, y= 0.8, label="myrmecophytic clade", font=2, cex=2.0, pos=4)


#rate through time for clad3

sp3 = "species6483"
sp4 = "species303"
tipnode3 <- which(tree$tip.label == sp3)
tipnode4 <- which(tree$tip.label == sp4)
mrca2  <- getMRCA(tree, tip = c(tipnode3, tipnode4))
plotRateThroughTime(edata,node = 65, ratetype = "speciation", intervals = seq(from = 0.05, 0.95, by = 0.01), ylim = c(0,1), yticks = 10)


#refer to
#Rabosky, D.L., 2014. Automatic detection of key innovations, 
#rate shifts, and diversitydependence on phylogenetic trees. PLoS One 9, e89543. h
#Rabosky, D.L., Grundler, M., Anderson, C., Title, P., Shi, J.J., Brown, J.W., Huang, H.,
#Larson, J.G., 2014. BAMM tools: an R package for the analysis of evolutionary
#dynamics on phylogenetic trees. Methods Ecol. Evol. 5, 701–707.
