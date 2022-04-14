# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector indicating samples groups (control: 0),
#        "datatype" is should be either "microarray" or "rnaseq",
#        "noisedefault" is the default cutoff (numeric value) for real signal
# OUTPUT: Returns a data frame (cols: samples, rows: pathways) of pathway scores. 

FilterNoise = function(data, grouplabels, datatype="rnaseq", noisedefault="automatic"){

  # Check that argument datatype is ok
  datatypes = c("microarray", "rnaseq")
  message = "Argument datatype should be either 'microarray' or 'rnaseq'."
  if(!(datatype %in% datatypes)) stop(message)
  
  # Check that argument noisedefault is ok
  message = "Argument noisedefault should be either a numeric value or 'automatic' (default)"
  if(!is.numeric(noisedefault) & (noisedefault != "automatic")) stop(message)
  
  # Detect noise
  if(is.numeric(noisedefault)){
    noiselevel = noisedefault
  } else{
    noisedefault = c(6, 3.3)[which(datatypes %in% datatype)]
    noiselevel = DetectNoise(data, noisedefault)
  } 
  
  # Remove genes that are just noise in all groups
  realsignal = apply(data, 1, function(x){
    groupsignal = unlist(lapply(split(x, labels), median, na.rm=T))
    keep = T
    if(all(groupsignal < noise)) keep = F
    return(keep)
  })
  
  filtereddata = data[realsignal,] 
  return(filtereddata)
}