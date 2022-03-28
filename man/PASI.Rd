\name{PASI}
\alias{PASI}

\title{
Pathway analysis for sample-level information
}
\description{
Main function (and the only one for external use) of PASI returns pathway values in data frame and writes them into a .txt file.
}
\usage{
PASI(data, grouplabels, pathwayadress=NULL, datatype="microarray", noisedefault="automatic", score="activity", nodemin=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data frame including gene expression data (rows: Entrez genes, cols: samples) as log2 values }
  \item{grouplabels}{ An integer vector indicating sample groups. If argument 'score' is set to "deregulation", control samples should be indicated with label 0. }
  \item{pathwayadress}{ A path to .xml pathways downloaded from KEGG, if NULL (default), pathways are automatically accessed via KEGG api}
  \item{datatype}{ Either "microarray" or "rnaseq" (default) }
  \item{noisedefault}{ The default cutoff (numeric value or character "automatic") for real signal }
  \item{score}{ Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect}
  \item{nodemin}{ Pathways with less than **nodemin** measured nodes will be excluded from the analysis}
}
\details{ 
The argument ''grouplabels'' should correspond to the columns in file ''data''.
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
# Run PASI for data frame mydata (microarray data), which includes samples from two groups to obtain pathway activity scores.
# pasiexample = PASI(mydata, grouplabels=c(rep(0,19),rep(1,15)))
}

\keyword{ design }
