library(diversitree)
tree = read.tree("reduced_species_tree.tre")
is.ultrametric(tree)
plot(tree)
tree$tip.label
data = read.table("character_data.txt", row.names=1)


#convert data table into a vector that is a component of the tree file

data.v = data[,1]
names(data.v) = row.names(data)
tree$tip.state = data.v

###########################################################################################

#checking for character-specific diversification rates with 8 different models# 

###########################################################################################

f = c(0.5, 0.931) #this vector contains details on sampling. c(0.5, 0.931) means 50% of all species in state 0 are sampled
                  #while 93.1% of all species in state 1 are sampled. 

## MODEL 1 --> FULL MODEL (lambda1 != lambda0, mu1 != mu0, q01 != q10) ##

model1 = make.bisse(tree=tree,states=tree$tip.state, sampling.f = f)   #make BiSSE

par1 = starting.point.bisse(tree)                   #determine starting point for parameter values

model1.fit <-find.mle(model1, par1)                       #do the actual ML search

model1.fit$lnLik

coefficients(model1.fit)

## MODEL 2 (lambda1 != lambda0, mu1 != mu0, q01 = q10) ##

model2 = constrain(model1, q01 ~ q10)        #constrain model1 to get model2

par2 = p[argnames(model2)]                   #drops the redundant q01 term from the starting parameter vector

model2.fit <- find.mle(model2, par2)         #do the ML search

model2.fit$lnLik                            

coefficients(model2.fit)

## MODEL 3 (lambda1 != lambda0, mu1 = mu0, q01 = q10) ##

model3 = constrain(model1, q01 ~ q10, mu0 ~ mu1)    #constrain model1 (q and mu) to get model3

par3 = p[argnames(model3)]                          #drops the redundant mu2 and q01 from the starting parameter vector

model3.fit = find.mle(model3, par3)                 #do the ML search

model3.fit$lnLik

coefficients(model3.fit)

## MODEL 4 (lambda1 != lambda0, mu1 = mu0, q01 != q10) ##

model4 = constrain(model1, mu0 ~ mu1)          #constrain model1 (mu) to get model4

par4 = p[argnames(model4)]                     #drops the redundant mu term from the starting parameter vector

model4.fit = find.mle(model4, par4)           #do the ML search

model4.fit$lnLik

coefficients(model4.fit)

## MODEL 5 (lambda1 = lambda0, mu1 != mu0, q01 = q10) ##

model5 = constrain(model1, lambda1 ~ lambda0, q01 ~ q10)          #constrain model1 (lambda, q) to get model5

par5 = p[argnames(model5)]                     #drops the redundant mu term from the starting parameter vector

model5.fit = find.mle(model5, par5)           #do the ML search

model5.fit$lnLik

coefficients(model5.fit)

## MODEL 6 (lambda1 = lambda0, mu1 = mu0, q01 != q10) ##

model6 = constrain(model1, lambda1 ~ lambda0, mu1 ~ mu0)          #constrain model1 (lambda, mu) to get model6

par6 = p[argnames(model6)]                     #drops the redundant lambda, mu term from the starting parameter vector

model6.fit = find.mle(model6, par6)           #do the ML search

model6.fit$lnLik

coefficients(model6.fit)

##  MODEL 7 (lambda1 = lambda0, mu1 != mu0, q01 != q10) ##

model7 = constrain(model1, lambda1 ~ lambda0)          #constrain model1 (lambda) to get model7

par7 = p[argnames(model7)]                     #drops the redundant lambda term from the starting parameter vector

model7.fit = find.mle(model7, par7)           #do the ML search

model7.fit$lnLik

coefficients(model7.fit)

##  MODEL 8 (lambda1 = lambda0, mu1 = mu0, q01 = q10) ##

model8 = constrain(model1, lambda1 ~ lambda0, mu1 ~ mu0, q01 ~ q10)          #constrain model1 (lambda, mu and q) to get model8

par8 = p[argnames(model8)]                     #drops the redundant lambda term from the starting parameter vector

model8.fit = find.mle(model8, par8)           #do the ML search

model8.fit$lnLik

coefficients(model8.fit)

#comparing models with anova 

anova(model1.fit, model2.fit, model3.fit, model4.fit, model5.fit, model6.fit, model7.fit, model8.fit)

## analysis with Bayesian MCMC

prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))    # prior as 1/2r where r is the character independent net diversification rate (lambda 0 - mu0)   

## step size tuning

set.seed(1)
tmp = mcmc(model1,model1.fit$par, nsteps = 100, prior = prior, lower = 0, w= rep(1,6), print.every = 0)
w <- diff(sapply(tmp[2:7], range))
w

## run the chain for 10000 generations

samples <- mcmc(model1, model1.fit$par, nsteps=10000, w=w, lower=0, prior=prior,
                print.every=0)

##plotting posterior distributions of speciation rates for the two states

col <- c("#004165", "#eaab00")
profiles.plot(samples[c("lambda0", "lambda1")], col.line=col, las=1,
              xlab="Speciation rate", legend="topright")


#refer to 
#Analysing diversification with diversitree Rich FitzJohn, with GeoSSE by Emma Goldberg 25 March 2012, version 0.9-2
#FitzJohn, R. G., Goldberg, E. E., & Magnuson-Ford, K. (2012). diversitree: Comparative phylogenetic analyses of diversification. R package version 0.9–1.
