# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector (starting from 0) indicating which samples belong to the same group,
#        "pathwayadress" is path to .xml pathways downloaded from KEGG, 
#        "datatype" is a character like "rnaseq_tpm" describing the type of genomic input data in file "data",
#        "noisedefault" is the default cutoff (numeric value) for real signal, and
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# OUTPUT: Checks that all arguments are valid (gives an error if not) and sets some default values. 

CheckInput = function(data, grouplabels, pathwayadress, datatype, noisedefault, score){
  
  datatypes = c("microarray", "rnaseq")
  message = "Argument datatype should be either 'microarray' or 'rnaseq'."
  if(!(datatype %in% datatypes)) stop(message)
  
  if(noisedefault!="automatic"){
    message = "Argument noisedefault should be either a numeric value, NA, or 'automatic' (default)"
    if(!is.numeric(noisedefault) & (!is.na(noisedefault))) stop(message)
  }
  
  if(noisedefault == "automatic"){
    noisedefault = c(6, 3.3)[datatypes %in% datatype]
  }
  
  if(!is.null(pathwayadress)){
    if(!dir.exists(pathwayadress)) stop("Argument pathwayadress is not null or a proper adress. Windows users: remember to change \ into /")
  }
  
  if(!(score %in% c("activity", "deregulation"))) stop("Argument score should be either 'activity' or 'deregulation'.")
  
  if(!is.data.frame(data)){
    stop("Argument 'data' should be a data frame")
  }
  
  # If it seems that the data is not in log scale, do the log transformation
  if(max(data, na.rm=T) > 25) data = log2(data+1)
  
  # Check that all columns in data have a lable
  if(length(grouplabels) != ncol(data)) stop("Number of samples in 'data' and 'grouplabels' are different.")
  
  # Check that in case of deregulation score, some samples are defined as controls
  if(score=="deregulation"){
    controls = sum(grouplabels == 0)
    message = "Deregulation score requires control samples (indicated by 0 in 'grouplabels')."
    if(controls < 1) stop(message)
    if(controls < 10){
      warning("Deregulation scores are more reliable with at least 10-15 control samples (0 in grouplabels)")
    }
  } 
  
  return(list(data, noisedefault))
}