#I used the analysis.pca module within the ipyrad package for principal component analysis. It is very easy and very informative. This script is mostly unmodified from that on the ipyrad tutorial page, but 
#I have made some changes, especially in the SNP filtration step. 
#!/usr/bin/env python
# coding: utf-8

# In[1]:
import ipyrad.analysis as ipa
import pandas as pd
import toyplot
import toyplot.pdf
import toyplot.png


# In[2]:

#plug in the path to your ipyrad assembly .snps.hdf5 output file.
data = "./beccariana_outfiles/beccariana.snps.hdf5"


# In[3]:
#As in my structure analysis, I have chosen not to use imap dictionaries for filtering SNPs for the same reason as I mentioned before: assuming populations before the analysis can potentially bias your results.
#I have initiated an imap dictionary later (see below) for color coding. This is helful of course, because you can see clearly when specimens from two different presumed populations cluster close togther. 
#PCA is very sensitive to missing sites. Either going with a very high mincov value to greatly reduce missing sites or using an imputation technique is highly recommended. 
#One way to impute is called "sample", which is done as below. I don;t like this method, because missing sites are imputed with randomly selected genotypes based on allele frequencies in previously
#defined populations in the imap dictionary. Again, biased. 

pca = ipa.pca(
       data = data,
       mincov=0.75,
       impute_method="Sample",
)

#In[4]:
#OR imputation can be done with the k-means method where genotypes are randomly sampled based on frequency of alleles in kmeans cluster generated populations. 
#here you don;t depend on user-defined popualtions but as with all kmeans clustering, you assume a specific number of populations beforehand. If you go with this, it would be a 
#good idea to try out differnt k values. This imputation technique can be done as follows. Here, k=5 populations are assumed. 

pca = ipa.pca(
       data = data,
       mincov=0.75,
       impute_method=5,
)

#In[5]
#If you don't want to impute missing sites at all, which is usually not recommended, you make sure mincov is high enough to greatly reduce missing sites.
#Trying out all imputation techniques to see what fits the best is a good way to go

pca = ipa.pca(
       data = data,
       mincov=0.95,
       impute_method=None,
)


# In[6]:
#initiating imap dictionary for coloring
imap = {
      "pop1": ["sample1", "sample2"],
    "pop2" : ["sample3", "sample4", "sample5", "sample6", "sample7", "sample9"],
    "pop3" : ["sample8", "sample11", "sample12", "sample22", "sample15", "sample16"],
    "pop4" : ["sample10", "sample13", "sample14", "sample17", "sample18"],
    "pop5" : ["sample19", "sample20", "sample21", "sample23", "sample24"],
}

pca.imap = imap


# In[7]:
#running PCA with multiple replicates is a good idea to avoid the effect of linkage among SNPs. Each replicate is generated by randomly sampling unliked SNPs and then running PCA on it. 

pca.run(nreplicates = 100, seed = 12345)


# In[8]:
# store the PC axes as a dataframe
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)

# write the PC axes to a CSV file
df.to_csv("pca_analysis.csv")

# show the first ten samples and the first 10 PC axes
df.iloc[:10, :10].round(2)


# In[9]:
#It is a good idea to plot differnt pairs of PC axes to get an idea of your data, but of course, the gratest variation is seen along 0 and 1 axes. 

canvas, axes = pca.draw(1,0)


# In[91]:


toyplot.pdf.render(canvas, "pca(0,1).pdf")
toyplot.pdf.render(canvas, "pca(0,1).png")
