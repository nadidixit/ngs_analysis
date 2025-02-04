# ngs_analysis
Contains Python and R code (written/borrowed/modified) and further information on bioinformatic methodology that was employed in the analysis of NGS data, specifically GBS, to infer the spatial-temporal evolutionary history of myrmecophytism in the plant genus Macaranga.

The idea behind this repository is to create a concise collection and documentation of tools and methodologies commonly applied in evolutionary studies that make use of NGS genomic datasets, specifically those generated using restriction-enzyme digestion. I have added suggestions and comments to help the reader through the files. 

Deeper description: All files in this repository are connected to the doctoral thesis titled "Spatial-temporal evolution of the Southeast Asian Macaranga (Euphorbiaceae) ant-plant lineages: a next-generation sequencing-based analysis", that was authored by me and supervised by Dr. Daniela Guicking. The goal of the doctoral project was to gain deeper insights into the evolution of myrmecophytism, an obligate mutualistic association between plants and ants, by using the plant genus Macaranga of the family Euphorbiaceae as a model system. The project encompassed methodologies such as NGS data assembly, SNP genotyping, phylogenetic tree calculations, molecular dating, parametric biogeography, population structure analysis, PCA, ancestral state reconstruction, and diversification rate analysis methods. These methodologies heavily relied on the ipyrad python package (Eaton & Overcast, 2020) and the modules and functions within. Additional methodologies relied on specific tools such as BEAST2 (Bouckaert et al., 2019), RASP (Yu et al., 2020), and Mesquite (Maddison & Maddison, 2021), and several packages on R. 

The application of aforementioned tools is peer-reviewed. Analyses, methodologies, and interpretations can be found in Dixit et al. (2023) and Dixit & Guicking (2024). 


References:
Eaton, D.A., Overcast, I., 2020. ipyrad: interactive assembly and analysis of RADseq datasets. Bioinformatics 36, 2592–2594. https://doi.org/10.1093/bioinformatics/btz966.
Bouckaert, R., Vaughan, T.G., Barido-Sottani, J., Duchêne, S., Fourment, M.,Gavryushkina, A., Heled, J., Jones, G., Kühnert, D., De maio, N., et al., 2019. BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. PLoS Comput. Biol. 15, e1006650
Yu, Y., Blair, C., He, X., 2020. RASP 4: ancestral state reconstruction tool for multiple genes and characters. Mol. Biol. Evol. 37, 604–606.
Maddison, W.P., Maddison, D.R., 2021. Mesquite: a modular system for evolutionary analysis. Version 3, 70.
Dixit, N.M., Zirpel, M., Slik, J.W., Jamsari, J., Weising, K., Guicking, D., 2023. Biogeography of the Sunda Shelf revisited: insights from Macaranga section Pruinosae (Euphorbiaceae). Front. Ecol. Evol. 10, 1049243.
Dixit, N. M., & Guicking, D., 2024. Exploring the evolutionary dynamics of myrmecophytism: Perspectives from the Southeast Asian Macaranga ant-plant symbiosis. Molecular Phylogenetics and Evolution 194, 108028.
