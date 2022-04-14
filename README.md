# PASI description

PASI (Pathway Analysis for Sample-level Information) is a pathway analysis tool that provides pathway values for each analyzed sample and pathway separately. Details about the algorithm are available in the original open access publication "PASI: A novel pathway method to identify delicate group effects" [1] (please cite it if you utilize PASI in your research). In case of bugs, missing documentation, or problems with installation, please contact us: maria.jaakkola@utu.fi

PASI can be installed by opening R and typing devtools::install_github("elolab/PASI") (requires package devtools to be installed).

### PASI Input and output

The only mandatory PASI inputs from the user are **data** and **grouplabels**.
| Input | Description |
| ----------- | ----------- |
| data | A data frame of gene expression data in log2 scale (rows: Entrez genes, cols: samples) |
| grouplabels | An integer vector indicating sample groups. If argument **score** is set to "deregulation", control samples should be indicated with label 0. |
| pathwayadress | NULL (default) or a path to custom pathways provided by the user (see format details below). If this is set to NULL, pathway structures are accessed automatically from KEGG. |
| useKEGG | TRUE (default) or FALSE indicating, if the KEGG pathways should be automatically accessed from API. Reguires internet connection. |
| score | Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect (more details below). |
| nodemin | An integer defining the minimum for how many measured nodes a pathway should include. Pathways with less measured nodes are excluded from the analysis. |

The argument **grouplabels** should be an integer vector with as many elements as columns (i.e. samples) in the gene expression data. Number 0 indicates a control sample. In case the expression data contains no sample groups, **grouplabels** can be set to dummy value of only zeros rep(0,ncol(data).

File format for the user defined pathways in **pathwayadress** is a .txt file (one file per pathway) and the directory should not include other .txt files than pathway files. Each pathway file should include the pathway's name as the first row and the following rows including nodes and relations. Node-lines should include only the Entrez gene id of the node. Relation-lines include three elemnts separated with a space: Entrez gene id of the start node (first), Entrez gene id of the end node (second), and either + or - indicating activation or inhibition, respectively (third). A toy example with five nodes (Entrez gene ids 5269, 8828, 10938, 1203, 8824) and three relations (5269 and 8828 activating 10938, and 8828 inhibiting 1203) is given below.

Pathway name here  
5269  
8828  
10938  
1203  
8824  
5269 10938 +  
8828 10938 +  
8828 1203 -

There are two options for argument **score**. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If **score** is set to "deregulation" (only this option was available in versions 1.0 and 1.1 of PASI and it is the approach in the original publication), the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.

Function PASI returns a data frame where columns are samples, rows are pathways and elements are sample specific pathway scores. The interpretation of a score is dependent on argument **score** (see above). The same results are written into a file named as "PASI YYYY-MM-DD.txt" in the current working directory (type getwd() to see what is your working directory at the moment).

### FilterNoise Input and output

Besides the main function PASI doing the pathway analysis, the R package offers a function FilterNoise for exluding genes with low expression in all sample groups.

| Input | Description |
| ----------- | ----------- |
| data | A data frame of gene expression data (rows: genes, cols: samples) |
| grouplabels | An integer vector indicating sample groups.|
| datatype | Either "microarray" or "rnaseq" (default) |
| noisedefault | The default cutoff (numeric value or character ”automatic” (default)) for real signal. |

Function FilterNoise returns a data frame similar to data, but without genes (rows) with median expression below a cutoff in all sample groups. The cutoff for real signal and noise is searched based on the provided data using either the value provided in argument **noisedefault** or the data type specific default (3.3 for RNAseq and 6 for microarray) as a starting point. If no cutoff in which normal distribution for noise and real signal overlap is found, the default value is used. 

##### References

[1] Maria K Jaakkola, Aidan J McGlinchey, Riku Klén, and Laura L Elo:
PASI: A novel pathway method to identify delicate group effects.
PLoS one 2018; 13(7):e0199991.