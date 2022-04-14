# Input: "nodevalues" is a list (pathways) of matrices with cols as samples and rows as pathway nodes,
#        "relationvalues" is a list similar to "nodevalues", but matrix rows are relations,
#        "importancefactors" is a list (pathways) of scalars for each node in a pathway.
# Output: A matrix (cols: samples, rows: pathways) of pathway activity scores

CalculatePathwayValues_Activity = function(nodevalues, relationvalues, importancefactors){
  
  numberofsamples = ncol(nodevalues[[1]])
  
  # Calculate average node value scaled with importances (list of pathways)
  nodemeans = mapply(function(x, importance){
      scalednodes = importance*x
      averagenodevals = colMeans(scalednodes, na.rm=TRUE) 
      return(averagenodevals)},
    x = nodevalues, importance = importancefactors, SIMPLIFY=FALSE)
  
  # Calculate average relation value (list of pathways)
  relationmeans = lapply(relationvalues, function(y){
    if(!all(is.na(y))){ averagerelationvals = colMeans(y, na.rm=TRUE)
    } else averagerelationvals = rep(NA, numberofsamples)
    return(averagerelationvals)
  })
  
  # Calculate pathway specific weights for node and relation values based on their ratio
  nodeweights = mapply(function(nodeval, relval){
    nnod = nrow(nodeval)                             # SHOULD I USE ONLY THE NUMBER OF MEASURED NODES/RELATIONS?
    if(!all(is.na(relval))){ nrel = nrow(relval)
    } else nrel = 0
    return(nnod/(nrel+nnod))
  }, nodeval = nodevalues,  relval = relationvalues)
  relationweights = 1 - nodeweights
  
  # Combine relation and node scores into pahway scores
  nodeeffect = as.matrix(data.frame(nodemeans)) # rows: samples, cols: pathways
  relationeffect = as.matrix(data.frame(relationmeans))
  relationeffect[is.na(relationeffect)] = 0
  pathwayscores = nodeweights*nodeeffect + relationweights*relationeffect
  pathwayscores = t(pathwayscores)
  
  return(pathwayscores)
}



# Input: "nodevalues" is a list (pathways) of matrices with cols as samples and rows as pathway nodes, 
#        "importancefactors" is a list (pathways) of scalars for each node in a pathway.
# Output: A matrix (cols: samples, rows: pathways) of pathway deregulation scores   

CalculatePathwayValues_Deregulation = function(nodevalues, importancefactors){
  
  # Initialize results
  results = mat.or.vec(length(nodevalues), ncol(nodevalues[[1]]))
  colnames(results) = colnames(nodevalues[[1]])
  rownames(results) = names(nodevalues)
  
  for(i in 1:length(nodevalues)){ # i = pathway
    
    pathway = nodevalues[[i]]
    importance = importancefactors[[i]]
    
    # Calculate scaled ranks
    scores = t(apply(pathway, 1, function(x){
      ranks = rank(x, na.last="keep")
      values = seq(1, ncol(pathway), length.out=length(ranks))
      scaledranks = values[ranks]
      return(scaledranks)
    }))
    
    # Multiply each node with its importance
    scaledscores = importance*scores
    
    # Calculate and save pathway scores for pathway i
    pathwayscores = colMeans(scaledscores, na.rm=TRUE)
    results[i,] = pathwayscores
  }
  
  return(results)
}