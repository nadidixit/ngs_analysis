install.packages("hisse")
library(hisse)
library(diversitree)

#load tree and character dataset

tree = read.tree("reduced_species_tree.tre")
plot(tree)
is.ultrametric(tree)
data = read.table("character_data.txt", row.names=1)

#convert data table into a vector that is a component of the tree file

data.v = data[,1]
names(data.v) = row.names(data)
tree$tip.state = data.v                  
data = data.frame(names(tree$tip.state), tree$tip.state)

f = c(0.072, 0.931) #sampling proportions as explained in the BiSSE.R script


#let's start by setting up a full BiSSE (bisse) model with HiSSE

t.bisse = c(1,2)         #turnover rates different acorss both observed traits
e.bisse = c(1,2)        #extinction rates different across both observed traits

matrix.bisse = TransMatMakerHiSSE(hidden.traits = 0)
print(matrix.bisse)

bisse = hisse(phy = tree, data = data, f=f, turnover = t.bisse, eps = e.bisse, hidden.states = FALSE, trans.rate = matrix.bisse)
bisse_ = hisse(phy = tree, data = data, f= f, turnover = t.bisse, eps = c(1,1), hidden.states = FALSE, trans.rate = matrix.bisse)

#setting up a null BiSSE model (bisse2) where turnover and extinction rates are equal (but transition rates are different)

t.bisse2 = c(1,1)             #turnover rates same acrosss both traits
e.bisse2 = c(1,1)             #extinction rates same across both traits

matrix.bisse2 = TransMatMakerHiSSE(hidden.traits = 0)
print(matrix.bisse2)

bisse2 = hisse(phy = tree, data = data, f = f, turnover = t.bisse2, eps = e.bisse2, hidden.states = FALSE, trans.rate = matrix.bisse2)
bisse2_ = hisse(phy = tree, data = data, turnover = t.bisse2, eps = c(1,1), hidden.states = FALSE, trans.rate = matrix.bisse2)

#setting up full HiSSE model

t.hisse = c(1,2,3,4)
e.hisse = c(1,2,3,4)

matrix.hisse = TransMatMakerHiSSE(hidden.traits = 1)
print(matrix.hisse)

matrix.hisse[8] = 6
matrix.hisse[9] = 7       #modifying trasition matrix to form a full HiSSE model. complicated transitions involving changes in both characters are removed
matrix.hisse[14] = 8

hisse = hisse(phy = tree, data = data, f= f, turnover = t.hisse, eps = e.hisse, hidden.states = TRUE, trans.rate = matrix.hisse)
hisse_ = hisse(phy = tree, data = data, f =f,  turnover = t.hisse, eps = c(1,1,1,1), hidden.states = TRUE, trans.rate = matrix.hisse)

#setting up HiSSe where only one observed states are associated with both hidden states (only 0A, 1A, 1B) -- hisse2

t.hisse2 = c(1,2,0,3)
e.hisse2 = c(1,2,0,3)

matrix.hisse2 = TransMatMakerHiSSE(hidden.traits = 1)
print(matrix.hisse2)

matrix.hisse2[3] = 0
matrix.hisse2[8] = 4
matrix.hisse2[9] = 0        #modifying transition matrix to include only those transitions that match the model description 
matrix.hisse2[12] = 0
matrix.hisse2[14] = 3
matrix.hisse2[15] = 0

hisse2 = hisse(phy = tree, data = data, f=f, turnover = t.hisse2, eps = e.hisse2, hidden.states = TRUE, trans.rate = matrix.hisse2)
hisse2_ =  hisse(phy = tree, data = data, f=f, turnover = t.hisse2, eps = c(1,1,0,1), hidden.states = TRUE, trans.rate = matrix.hisse2)


#cid-2 as given in documentation

t.cid2_d = c(1,1,2,2)
e.cid2_d = c(1,1,2,2)

matrix.cid2_d = TransMatMakerHiSSE(hidden.traits = 1, make.null = TRUE)
print(matrix.cid2_d)

cid2_d = hisse(phy = tree, data = data, f=f, turnover = t.cid2_d, eps = e.cid2_d, hidden.states = TRUE, trans.rate = matrix.cid2_d)
cid2_d_ = hisse(phy = tree, data = data, f=f, turnover = t.cid2_d, eps = c(1,1,1,1), hidden.states = TRUE, trans.rate = matrix.cid2_d)

#setting up cid-4 as in documentation

t.cid4 = c(1,1,2,2,3,3,4,4)
e.cid4 = rep(1,8)

matrix.cid4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
print(matrix.cid4)

cid4 = hisse(phy = tree, data = data, f=f, turnover = t.cid4, eps = e.cid4, hidden.states = TRUE, trans.rate = matrix.cid4)
