# PASI description

PASI (Pathway Analysis for Sample-level Information) is a pathway analysis tool that provides pathway values for each analyzed sample and pathway separately. Details about the algorithm are available in the original open access publication "PASI: A novel pathway method to identify delicate group effects" [1] (please cite it if you utilize PASI in your research). In case of bugs, missing documentation, or problems with installation, please contact us: maria.jaakkola@utu.fi

PASI can be installed by opening R and typing devtools::install_github("elolab/PASI") (requires package devtools to be installed).

### Input and output

The only mandatory input from the user are **data** and **grouplabels**.
| Input | Description |
| ----------- | ----------- |
| data | A data frame of gene expression data in log2 scale (rows: Entrez genes, cols: samples) |
| grouplabels | An integer vector indicating sample groups. If argument **score** is set to "deregulation", control samples should be indicated with label 0. |
| pathwayadress | NULL (default) or a path to .xml pathways downloaded from KEGG. If this is set to NULL, pathway structures are accessed automatically from KEGG. |
| datatype | Either "microarray" (default) or "rnaseq" |
| noisedefault | The default cutoff (numeric value or character ”automatic” (default)) for real signal |
| score | Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect (more details below). |

The argument **grouplabels** should be an integer vector with as many elements as columns (i.e. samples) in the gene expression data. Number 0 indicates a control sample. In case the expression data contains no sample groups, **grouplabels** can be set to dummy value of only zeros rep(0,ncol(data).

There are two options for argument **score**. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If **score** is set to "deregulation" (only this option was available in versions 1.0 and 1.1 of PASI and it is the approach in the original publication), the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.

The main function (and the only one meant for external use) PASI returns a data frame where columns are samples, rows are pathways and elements are sample specific pathway scores. The interpretation of a score is dependent on argument **score** (see above). The same results are written into a file named as "PASI YYYY-MM-DD.txt" in the current working directory (type getwd() to see what is your working directory at the moment).


##### References

[1] Maria K Jaakkola, Aidan J McGlinchey, Riku Klén, and Laura L Elo:
PASI: A novel pathway method to identify delicate group effects.
PLoS one 2018; 13(7):e0199991.