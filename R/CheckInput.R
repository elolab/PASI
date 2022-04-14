# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector (starting from 0) indicating which samples belong to the same group,
#        "pathwayadress" is path to .xml pathways downloaded from KEGG, 
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# OUTPUT: Checks that all arguments are valid (gives an error if not) and sets some default values. 

CheckInput = function(data, grouplabels, pathwayadress, useKEGG, score, nodemin){ 

  # If no pathway adress is given set useKEGG to TRUE as some pathways need to be accessed
  if(is.null(pathwayadress)){
    if(!useKEGG) print("Argument useKEGG is set to TRUE bacause no custom pathways were provided and some pathways are needed for analysis")
    useKEGG = TRUE
  }

  # Check that all columns in data have a label
  if(length(grouplabels) != ncol(data)){
    stop("Number of samples in 'data' and 'grouplabels' are different.")
  } else{
    prefix = "Argument 'grouplabels' should indicate control samples with 0."
    suffix = "If the data contains only one sample group, use 0 for all samples in 'grouplabels'."
    message = paste(prefix, suffix, sep=" ")
    if(!any(grouplabels == 0)) stop(message)
  }
  
  if(!is.null(pathwayadress)){
    if(!dir.exists(pathwayadress)) stop("Argument pathwayadress is not null or a proper adress. Windows users: remember to change \ into /")
  }
  
  if(!(score %in% c("activity", "deregulation"))) stop("Argument score should be either 'activity' or 'deregulation'.")
  
  if(!is.data.frame(data)){
    stop("Argument 'data' should be a data frame. Tip: function as.data.frame() converts a matrix into a data frame.")
  }
  
  # If it seems that the data is not in log scale, do the log transformation
  if(max(data, na.rm=T) > 50) data = log2(data+1)
  
  # Check that in case of deregulation score, some samples are defined as controls
  if(score=="deregulation"){
    controls = sum(grouplabels == 0)
    message = "Deregulation score requires control samples (indicated by 0 in 'grouplabels')."
    if(controls < 1) stop(message)
    if(controls < 10){
      warning("Deregulation scores are more reliable with at least 10-15 control samples (0 in grouplabels)")
    }
  } 
  
  return(data)
}