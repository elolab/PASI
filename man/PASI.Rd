\name{PASI}
\alias{PASI}

\title{
Pathway analysis for sample-level information
}
\description{
Main function of PASI returns pathway values in data frame and writes them into a .txt file.
}
\usage{
PASI(data, grouplabels, pathwayadress=NULL, useKEGG=TRUE, score="activity", nodemin=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data frame including gene expression data (rows: Entrez genes, cols: samples) as log2 values. }
  \item{grouplabels}{ An integer vector indicating sample groups, control samples should be indicated with label 0. }
  \item{pathwayadress}{ NULL (default) or a path to user-defined pathway files. If this is set to NULL, only pathway structures accessed automatically from KEGG are used.}
  \item{useKEGG}{ TRUE (default) or FALSE indicating, if the KEGG pathways should be automatically accessed from API. Reguires internet connection.}
  \item{score}{ Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect}
  \item{nodemin}{ Pathways with less than ''nodemin'' (integer) measured nodes will be excluded from the analysis}
}
\details{ 
The argument ''grouplabels'' should correspond to the columns in file ''data''.

File format for the user defined pathways in ''pathwayadress'' is a .txt file (one file per pathway) and the directory should not include other .txt files than pathway files. Each pathway file should include the pathway's name as the first row and the following rows including nodes and relations. Node-lines should include only the Entrez gene id of the node. Relation-lines include three elemnts separated with a space: Entrez gene id of the start node (first), Entrez gene id of the end node (second), and either + or - indicating activation or inhibition, respectively (third). A toy example with five nodes (Entrez gene ids 5269, 8828, 10938, 1203, 8824) and three relations (5269 and 8828 activating 10938, and 8828 inhibiting 1203) is given below.

Pathway name here  
5269  
8828  
10938  
1203  
8824  
5269 10938 +  
8828 10938 +  
8828 1203 -
}

\value{
A data frame where rows are pathways, columns are samples and elements are calculated pathway values.
}

\references{
MK Jaakkola et al. PASI: A novel pathway method to identify delicate group effects. PLOS ONE 2018, 13(7):e0199991.
}

\author{
Maria K. Jaakkola
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
# Run PASI for data frame mydata, which includes samples from two groups to obtain pathway activity scores.
# pasiexample = PASI(mydata, grouplabels=c(rep(0,19),rep(1,15)))
}

\keyword{ design }
