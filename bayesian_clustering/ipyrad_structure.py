#For population clustering analysis, I used the structure module within the analysis subpackage of ipyrad. This script is mostly unmodified from that on the ipyrad site, but I made a few changes to the way 
#SNPs are filtered using the imap dictionary. Although ipyrad.analysis.structure is a very useful tool, I strongly recommend reading the STRUCTURE documentation (cited at the end as comment)
#available on the Pritchard Lab site to familiarize yourself with the different parameter and model options. 
#Best K vaule is determined with the Evanno method here, but since this method has its own shortcomings, I recommend checking the plots for multiple K values and going with the 
#one that seems biologically the most informative, regardless of the deltaK values. I would also suggest going with a fist batch of assumed K values (for example 2, 3, 4, 5, 6) 
#and then checking how the plots change with differnt K values. You can stop to include higher K values when the plots remain unchanged and start showing 
#a "horizontal patterning" (you will know it when you see it)

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ipyrad.analysis as ipa
import toyplot
import toyplot.pdf
import toyplot.png


# In[2]:
#plug in the path to your ipyrad assembly .snps.hdf5 output file here
data = "./ipyrad_output_file.snps.hdf5"


# In[3]:
#intiate a structure object. In the ipyrad tutorial, two parameters imap and minmap are also used to filter SNPs. I had an issue with this method and I left it out. 
#imap is essentially a dcitionary where you group your samples into "presumed populations" (or even species). minmap is a percentage value (say x) which is used 
#to include SNPs that are present in at least x proportion of the samples in each group (as defined in imap). Since this filtration is based
#on assumed populations and the point of STRUCTURE is to discern distinct populations, this step can potentially bias your results by including
#SNPs common to predefined populations. I found this problematic, and hence only included mincov which is used to filter out SNPs that are present 
#in less than a certain proportion of all samples (not dependent on any predefined grouping). I did initiate an imap dictionary later on, but this was used only to 
#determine the order of samples in the plot
struct1 = ipa.structure(
name = "analysis1.str",
data = data,
mincov = 0.75)


# In[4]:
#In the STRUCTURE documentation, Pritchard suggets that a run length of 100,000 is usually sufficent. As of now, it is not possible to check for convergence
#with the ipyrad.analysis.structure tool
struct1.mainparams.burnin = 50000
struct1.mainparams.numreps = 100000


# In[5]:
#You can print the paramter values with the following two commands. Most values are deafult vlaues recommended in the STRUCTURE docuemtation (cited at the end as comment) with the
#"admixture" ancestry model chosen. You can, of course, change them if you like. The default values and the admixture model were appropriate for my 
#data and they worked well for me. 

print(struct1.mainparams)
print(struct1.extraparams)


# In[6]:
#Running the analysis
struct1.run(nreps=10, kpop=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], auto=True)


# In[7]:
#I initiated the imap dictionary here, but only for plotting purposes, not for filtering SNPs
imap = {
    "pop1": ["sample1", "sample2"],
    "pop2" : ["sample3", "sample4", "sample5", "sample6", "sample7", "sample9"],
    "pop3" : ["sample8", "sample11", "sample12", "sample22", "sample15", "sample16"],
    "pop4" : ["sample10", "sample13", "sample14", "sample17", "sample18"],
    "pop5" : ["sample19", "sample20", "sample21", "sample23", "sample24"],
  
}


# In[8]:

etable = struct1.get_evanno_table([2,3,4,5, 6, 7, 8, 9, 10, 11, 12])
etable.astype('float').dtypes
etable.dtypes
etable


# In[9]:
#plotting the deltaK values. Reading Evanno et al.(2005) before you interpret this plot is a very good idea

# get canvas object and set size
canvas = toyplot.Canvas(width=400, height=300)

# plot the mean log probability of the models in red
axes = canvas.cartesian(ylabel="estLnProbMean")
axes.plot(etable.estLnProbMean * -1, color="darkred", marker="o")
axes.y.spine.style = {"stroke": "darkred"}

# plot delta K with its own scale bar of left side and in blue
axes = axes.share("x", ylabel="deltaK", ymax=etable.deltaK.max() + etable.deltaK.max() * .25)
axes.plot(etable.deltaK, color="steelblue", marker="o");
axes.y.spine.style = {"stroke": "steelblue"}

# set x labels
axes.x.ticks.locator = toyplot.locator.Explicit(range(len(etable.index)), etable.index)
axes.x.label.text = "K (N ancestral populations)"


# In[10]:
#saving the plot as a .png file. You can also save it as a pdf by using toyplot.pdf.render(), but you need to import the toyplot.pdf module first.
toyplot.png.render(canvas, "deltaK.png")


# In[11]:
#selecting a k value for your structure plot.
k = 2
table = struct1.get_clumpp_table(k)

# In[12]:
# sort list by columns
table.sort_values(by=list(range(k)), inplace=True)

# or, sort by a list of names (here taken from imap)
import itertools
onames = list(itertools.chain(*imap.values()))
table = table.loc[onames]


# In[13]:
# build barplot
canvas = toyplot.Canvas(width=1000, height=400)
axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
axes.bars(table)

# add labels to x-axis
ticklabels = [i for i in table.index.tolist()]
axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
axes.x.ticks.labels.angle = -60
axes.x.ticks.show = True
axes.x.ticks.labels.offset = 10
axes.x.ticks.labels.style = {"font-size": "12px"}


# In[14]:
#I would recommend to plot barplots for many k values and checking which one seems more relevant and informative, regardless of the deltak vlaues
toyplot.png.render(canvas, "k2_plot.png")



#References
#Eaton, D. A. R., and Overcast, I. (2020). Ipyrad: interactive assembly and analysis of RADseq datasets. Bioinformatics 36, 2592–2594.
#Evanno, G., Regnaut, S., and Goudet, J. (2005). Detecting the number of clusters of individuals using the software structure: a simulation study.
#Pritchard, J. K., Stephens, M., and Donnelly, P. (2000). Inference of population structure using multilocus genotype data. Genetics 155, 945–959.
#Pritchard, J. K., Wen, W., and Falush, D. (2010). Documentation for STRUCTURE Software: Version 2. University of Chicago, Chicago, IL.
