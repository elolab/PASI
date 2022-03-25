# Input: 'data' is a numeric matrix (cols: samples, rows: pathways) with wanted effects neutralized
#
#
#
PickImportant = function(data, mainfeature, seed=252525){
  
  # Pick samples with the main feature defined (cases)
  index = which(!is.na(mainfeature))
  feature = mainfeature[index]
  usedata = data[,index]
  
  # Initialize coefficient matrix
  if(is.numeric(feature)){
    cols = 1
    coefnames = "feature"
  } else{
    cols = length(unique(feature))
    coefnames = paste("feature", unique(feature), sep="")
  }
  coefficient = matrix(0, nrow=nrow(data), ncol=cols)
  rownames(coefficient) = rownames(data)
  colnames(coefficient) = coefnames
  
  # Fit model for those samples
  for(i in 1:nrow(data)){
    modelres = rlm(as.numeric(usedata[i,])~feature)
    coefficient[i,coefnames] = modelres$coefficients[coefnames]
  }
  
  # Initialize coefficient matrix for randomized time labels
  coefs = array(0, dim=c(1000, nrow(usedata), ncol(coefficient)))
  featurenames = paste("temp",colnames(coefficient),sep="")
  dimnames(coefs) = list(paste("R",1:1000,sep=""),rownames(usedata),featurenames)
  
  # Calculate coefficients with randomly sampled features (takes forever)
  for(i in 1:nrow(coefs)){
    set.seed(seed+i)
    tempfeature = sample(feature, size=length(feature), replace=FALSE)
    for(j in 1:nrow(usedata)){
      coefs[i,j,featurenames] = rlm(as.numeric(usedata[j,])~tempfeature)$coefficients[featurenames]
    }
  }
  
  # Initialize a matrix for p-value calculation
  realgreater = matrix(0, nrow=nrow(usedata), ncol=ncol(coefficient))
  rownames(realgreater) = rownames(usedata)
  colnames(realgreater) = gsub("feature","",colnames(coefficient))
  
  # Calculate p-values and convert them into fdr
  for(n in 1:dim(coefs)[1]) realgreater = realgreater + (coefficient >= as.matrix(coefs[n,,]))
  realsmaller = nrow(coefs) - realgreater
  pvals = 2*pmin(realgreater/nrow(coefs), realsmaller/nrow(coefs))
  fdr = apply(pvals, 2, p.adjust, method="fdr")
  
  # Combine p-values and fdr into one significance matrix
  significance = cbind(pvals, fdr)
  rownames(significance) = rownames(pvals)
  colnames(significance) = c(paste("pval",colnames(pvals),sep="_"), paste("fdr",colnames(pvals),sep="_"))

  return(significance)
}