\name{FilterNoise}
\alias{FilterNoise}

\title{
FilterNoise
}
\description{
Removes genes with median expression below cutoff in all sample groups from given data. The cutoff for real signal and noise is searched based on the provided data using either the value provided in argument ''noisedefault'' or the data type specific default (3.3 for RNAseq and 6 for microarray) as a starting point. If no cutoff in which normal distribution for noise and real signal overlap is found, the default value is used.
}
\usage{
FilterNoise(data, grouplabels, datatype="rnaseq", noisedefault="automatic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data frame including gene expression data (rows: Entrez genes, cols: samples) as log2 values }
  \item{grouplabels}{ An integer vector indicating sample groups. If argument 'score' is set to "deregulation", control samples should be indicated with label 0. }
  \item{datatype}{ Either "microarray" or "rnaseq" (default) }
  \item{noisedefault}{ The default cutoff (numeric value or character "automatic") for real signal in log2 scale }
}
\details{ 
The argument ''grouplabels'' should correspond to the columns in file ''data''.
}

\value{
A data frame similar to ''data'', but without genes (rows) with median expression below a cutoff in all sample groups.
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
%%
}

\keyword{ design }
