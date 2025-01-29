#I used the raxml module available in the analysis subpackage of ipyrad to calculate trees in the maximum likelihood framework. ipyrad.analysis.raxml() automates the process of generating
#RAxML command line strings and running them through python code. This script is largely unmodified from that available on the ipyrad site. 
#!/usr/bin/env python
# coding: utf-8
# In[1]:

import ipyrad.analysis as ipa
import toytree
import toyplot.pdf

# In[2]:
#you plug in the path to your phylip formatted output file from the ipyrad assembly here. 
data1 = "./the_phylip_formatted_output_file.phy"

# In[3]:
#I used 100 bootstrap replicates, which is decent enough. This command by default chooses the -f a option elucidated in the RAxMl documentation that conducts a rapid Bootstrap analysis and search for the best-scoring ML tree in
#one single program run while choosing the generalized time reversible (GTR) model as the substitution model along with the GAMMA model for rate heterogeneity (GTRGAMMA).  
rax1 = ipa.raxml(T = 12, N = 100, data = data1, n = "run_name")

# In[4]:
#you can check the parameters chosen (by default or set by you) by running this command
print(rax1.command)

# In[5]:

rax1.run(block = True)

# In[6]:

tre = toytree.tree(rax1.trees.bipartitions)

# In[7]:
#this roots the tree with the two chosen samples
rtre = tre.root(names = ["Sample1", "Sample2"])

# In[8]:
#more options to plot the tree, including font and size can be explored in the toytree documentation (https://eaton-lab.org/toytree/quick_guide/)
canvas, axes = rtre.draw(tip_labels_align = True, node_labels = "support", node_sizes = 10, node_markers = "r2x1.25", node_labels_style = {"font-size": "9px"}, width = 1100, height = 1000)

# In[9]:
#saving your plotted tree as pdf
toyplot.pdf.render(canvas, "plotted_tree.pdf")

# In[10]:
#I just wanted to take a look at the bootstrap values for this tree as an array
bootstrap1 = rtre.get_node_values("support")
bootstrap1

#Refer to Eaton, D.A., Overcast, I., 2020. ipyrad: interactive assembly and analysis of RADseqdatasets. Bioinformatics 36, 2592–2594.
#Eaton, D.A., 2020. Toytree: a minimalist tree visualization and manipulation library for Python. Methods Ecol. Evol. 11, 187–191
#Stamatakis, A., 2014. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30, 1312–1313.


